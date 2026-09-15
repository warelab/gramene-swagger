'use strict';

// Variant resolution, I/O half (genotyping spec §2.3, §2.4, §2.7 rules 7-14, §3.2, §3.6-§3.9, §4.15): orchestration
// over the assembly resolver, one FASTA window per operation, the Ensembl client (client.js) and the pure
// normalizer (normalize.js).
//
//   hasVariation(system_name, deps)           -> Promise<boolean>
//   variationInfo(system_name, resolved, cfg) -> {available, source, release}          pure; GET /primers/genomes
//   listVariants(query, deps)                 -> PrimerVariantList                      GET /primers/variants
//   lookupVariant(query, deps)                -> PrimerVariantLookup                    GET /primers/variants/{variant_id}
//   resolveDesignVariant(request, deps)       -> the design's variant, assembly and sequence read (§2.7 rules 7-14)
//   neighboursFor(request, deps)              -> the design's known neighbours and their summary (§4.15)
//
// deps (all optional): {cfg, client (a createVariationClient instance), resolve(system_name) | assemblies |
//   (catalog, mongo), sequence {fetch, regionLength}, signal, log}. Without deps.client the process-wide client is
//   used, built from require('../config').get() (never from deps.cfg).
//
// Every Ensembl call goes through that client (limiter, single flight, breaker, caches) and nothing here holds a
// design semaphore slot, so the genotyping design calls resolveDesignVariant and neighboursFor before
// semaphore.acquire (§4.2). Overlap records are always fetched in the client's fixed chunk_bp chunks and filtered
// here, so a listing, a later design over the same locus and the UI's neighbour track share cache entries (§3.3).
// Each operation reads the genome once, max_shift + 250 bases beyond the span it needs (§3.6), and reads a wider
// window only when normalization reaches the edge (SequenceWindowError).

const { PrimerHttpError, isPrimerHttpError } = require('../errors');
const normalize = require('./normalize');
const { createVariationClient, VariationError, ID_RE } = require('./client');

const KINDS = Object.freeze(['snv', 'mnv', 'insertion', 'deletion', 'complex']);
const SPECIES_RE = /^[a-z0-9_]+$/;
const MAX_REGION_NAME = 255;
const SOURCE_NAME = 'ensembl';
// §3.6: every read extends max_shift + WINDOW_MARGIN bases beyond the span the operation needs.
const WINDOW_MARGIN = 250;
// Wider re-reads after a SequenceWindowError before the error is let through.
const MAX_WIDENINGS = 3;
const FEATURE_DISABLED_RETRY_AFTER_S = 300;
// Fallbacks for a partial test config (config.js DEFAULTS carry the real values).
const DEFAULT_TEMPLATE_FLANK = 400;
const DEFAULT_DENSE_WINDOW = 30;
const DEFAULT_DENSE_COUNT = 2;

let sharedClient = null;

// ---- config and small helpers -----------------------------------------------------------------------------

function configOf(deps) {
  return deps.cfg || require('../config').get();
}

// primers.variation over the normalizer's defaults, so a partial test config still normalizes.
function variationCfg(cfg) {
  return Object.assign({}, normalize.DEFAULTS, cfg && cfg.variation);
}

function variationEnabled(cfg) {
  return !!(cfg && cfg.variation && cfg.variation.enabled !== false);
}

// The Ensembl species segment configured for system_name, or null. Enabled or not.
function speciesConfigured(cfg, systemName) {
  const species = cfg && cfg.variation && cfg.variation.species;
  if (!species || typeof systemName !== 'string' || !Object.prototype.hasOwnProperty.call(species, systemName)) return null;
  const s = species[systemName];
  return typeof s === 'string' && SPECIES_RE.test(s) ? s : null;
}

function releaseOf(cfg) {
  const r = cfg && cfg.variation && cfg.variation.release;
  return r === undefined || r === null ? null : String(r);
}

function isPosInt(v) {
  return Number.isSafeInteger(v) && v >= 1;
}

function intValue(v) {
  if (Number.isSafeInteger(v)) return v;
  if (typeof v === 'string' && /^\d{1,15}$/.test(v)) return Number(v);
  return null;
}

function spanOverlaps(a, b, lo, hi) {
  return Math.min(a, b) <= hi && Math.max(a, b) >= lo;
}

function warning(code, message, details) {
  return { code: code, message: message, details: details };
}

function plural(n, one, many) {
  return n === 1 ? one : many;
}

function clientOf(deps) {
  if (deps.client) return deps.client;
  if (!sharedClient) sharedClient = createVariationClient({ cfg: require('../config').get(), log: console });
  return sharedClient;
}

async function resolveAssembly(systemName, deps) {
  if (typeof deps.resolve === 'function') return deps.resolve(systemName);
  const opts = {};
  ['cfg', 'mongo', 'catalog', 'log'].forEach(function (k) {
    if (deps[k] !== undefined) opts[k] = deps[k];
  });
  return (deps.assemblies || require('../assemblies')).resolve(systemName, opts);
}

// ---- errors ------------------------------------------------------------------------------------------------

function featureDisabled() {
  return new PrimerHttpError(503, 'FEATURE_DISABLED', 'known-variant lookups are disabled on this server',
    { retry_after_s: FEATURE_DISABLED_RETRY_AFTER_S });
}

function noVariationData(systemName) {
  return new PrimerHttpError(422, 'NO_VARIATION_DATA',
    systemName + ' has no known-variant data; enter the variant as region, position, ref and alt', { system_name: systemName });
}

function noSequence(systemName) {
  return new PrimerHttpError(422, 'NO_SEQUENCE', 'no genome sequence is available for ' + systemName, { system_name: systemName });
}

function invalidRequest(field, message) {
  return new PrimerHttpError(400, 'INVALID_REQUEST', message, { field: field });
}

function outOfBounds(region, details, regionLength) {
  return new PrimerHttpError(400, 'REGION_OUT_OF_BOUNDS', 'the requested coordinates are beyond the end of region ' + region +
    ' (' + regionLength + ' bp)', Object.assign({ region: region }, details, { length: regionLength }));
}

function isUnavailable(err) {
  return !!err && err.code === 'VARIATION_SOURCE_UNAVAILABLE';
}

async function requireRegionLength(sequence, fasta, region) {
  const len = await sequence.regionLength(fasta, region);
  if (!isPosInt(len)) {
    throw new PrimerHttpError(404, 'UNKNOWN_REGION', 'region ' + JSON.stringify(String(region)) + ' is not in this assembly',
      { region: String(region) });
  }
  return len;
}

// ---- warnings ----------------------------------------------------------------------------------------------

function duplicateIdsWarning(entry) {
  return warning('DUPLICATE_VARIANT_IDS', entry.ids.length + ' Ensembl ids describe the same event ' + entry.key + '; they were merged',
    { key: entry.key, ids: entry.ids.slice() });
}

function pushSkipped(warnings, skipped) {
  if (!skipped || skipped.count === 0) return;
  warnings.push(warning('VARIATION_RECORDS_SKIPPED', skipped.count + ' malformed Ensembl variation ' + plural(skipped.count, 'record was', 'records were') +
    ' skipped', { count: skipped.count, reasons: skipped.reasons }));
}

function pushRefMismatches(warnings, entries) {
  const count = entries.filter(function (e) { return e.ref_verified === false; }).length;
  if (count === 0) return;
  warnings.push(warning('REF_MISMATCHES', count + ' ' + plural(count, 'variant has a reference allele that differs', 'variants have reference alleles that differ') +
    ' from this genome', { count: count }));
}

// {count, reasons} over the client's skipped records inside [lo, hi] (or of unknown position) plus the
// normalizer's own skips. Records of two overlapping chunks are counted once.
function skippedSummary(clientSkipped, lo, hi, normSkipped) {
  const out = { count: 0, reasons: {} };
  const add = function (reason, n) {
    out.count += n;
    out.reasons[reason] = (out.reasons[reason] || 0) + n;
  };
  const seen = new Set();
  (clientSkipped || []).forEach(function (s) {
    if (s.start !== null && s.end !== null) {
      if (!spanOverlaps(s.start, s.end, lo, hi)) return;
      const key = s.start + '|' + s.end + '|' + s.reason;
      if (seen.has(key)) return;
      seen.add(key);
    }
    add(s.reason, 1);
  });
  if (normSkipped) Object.keys(normSkipped.reasons).forEach(function (r) { add(r, normSkipped.reasons[r]); });
  return out;
}

// ---- Ensembl records and the sequence read -------------------------------------------------------------------

function chunkStarts(lo, hi, chunkBp) {
  const out = [];
  for (let s = Math.floor((lo - 1) / chunkBp) * chunkBp + 1; s <= hi; s += chunkBp) out.push(s);
  return out;
}

// Records of the chunks covering [lo, hi] (clamped to the region), fetched concurrently under the client's
// limiter. -> {records, skipped, start, end}, where [start, end] is the span the chunks cover.
async function fetchChunks(ctx, region, regionLength, lo, hi) {
  lo = Math.max(1, lo);
  hi = Math.min(regionLength, hi);
  const starts = chunkStarts(lo, hi, ctx.v.chunk_bp);
  const chunks = await Promise.all(starts.map(function (s) {
    return ctx.client.overlapChunk(ctx.species, region, s, { signal: ctx.signal, regionLength: regionLength });
  }));
  const records = [];
  const skipped = [];
  chunks.forEach(function (c) {
    Array.prototype.push.apply(records, c.records);
    Array.prototype.push.apply(skipped, c.skipped || []);
  });
  return { records: records, skipped: skipped, start: starts[0], end: Math.min(regionLength, starts[starts.length - 1] + ctx.v.chunk_bp - 1) };
}

async function readWindow(ctx, region, regionLength, start, end) {
  const maxLength = ctx.cfg.design && ctx.cfg.design.max_fetch_length;
  const seq = await ctx.sequence.fetch(ctx.fasta, region, start, end, 1, maxLength ? { maxLength: maxLength } : undefined);
  return normalize.sequenceWindow(seq, start, regionLength);
}

// fn(genome) over one sequence read of [win.start, win.end] (clamped to the region), or over win.genome when the
// caller already holds a read of this region. A SequenceWindowError re-reads a window widened to the failed read
// plus max_shift + WINDOW_MARGIN, at most MAX_WIDENINGS times. -> {genome, value}
async function withGenome(ctx, region, regionLength, win, fn) {
  let genome = win.genome || null;
  let start = Math.max(1, win.start);
  let end = Math.min(regionLength, win.end);
  for (let attempt = 0; ; attempt++) {
    if (!genome) genome = await readWindow(ctx, region, regionLength, start, end);
    try {
      return { genome: genome, value: await fn(genome) };
    } catch (err) {
      if (!(err instanceof normalize.SequenceWindowError) || attempt >= MAX_WIDENINGS) throw err;
      const pad = ctx.v.max_shift + WINDOW_MARGIN;
      start = Math.max(1, Math.min(genome.start, err.needed.start - pad));
      end = Math.min(regionLength, Math.max(genome.end, err.needed.end + pad));
      genome = null;
    }
  }
}

function normalizeOptions(v, region, extra) {
  return Object.assign({ max_allele_length: v.max_allele_length, max_shift: v.max_shift, ems_source_pattern: v.ems_source_pattern, region: region },
    extra || {});
}

function makeContext(cfg, deps, resolved, species, fasta) {
  return {
    cfg: cfg,
    v: variationCfg(cfg),
    resolved: resolved,
    species: species,
    fasta: fasta,
    sequence: deps.sequence || require('../sequence'),
    client: species ? clientOf(deps) : null,
    signal: deps.signal,
    log: deps.log || console
  };
}

// The assembly, species segment and FASTA of a genome with variation data. Errors: 404 UNKNOWN_GENOME,
// 422 AMBIGUOUS_ASSEMBLY, 422 NO_VARIATION_DATA, 422 NO_SEQUENCE, 503 MONGO_UNAVAILABLE.
async function variationContext(systemName, deps, cfg) {
  const resolved = await resolveAssembly(systemName, deps);
  const species = speciesConfigured(cfg, systemName);
  if (species === null) throw noVariationData(systemName);
  const fasta = resolved && resolved.fasta && resolved.fasta.dna;
  if (!fasta) throw noSequence(systemName);
  return makeContext(cfg, deps, resolved, species, fasta);
}

// §2.3: the entries whose minimal span overlaps [start, end]. Records come from the chunks covering
// [start - 1, end + 1]; the records of those chunks up to max_shift past the end are normalized too, so a
// right-shifted description of an in-window event still merges into its entry. -> {entries, skipped, genome}
async function windowEntries(ctx, region, regionLength, start, end, genome) {
  const fetched = await fetchChunks(ctx, region, regionLength, start - 1, end + 1);
  const records = fetched.records.filter(function (r) { return spanOverlaps(r.start, r.end, start - 1, end + 1 + ctx.v.max_shift); });
  const pad = ctx.v.max_shift + WINDOW_MARGIN;
  const out = await withGenome(ctx, region, regionLength, { start: start - pad, end: end + pad, genome: genome }, function (g) {
    return normalize.recordsToEntries(records, g, normalizeOptions(ctx.v, region));
  });
  return {
    entries: out.value.entries.filter(function (e) { return normalize.inWindow(e.minimal, start, end); }),
    skipped: skippedSummary(fetched.skipped, start - 1, end + 1, out.value.skipped),
    genome: out.genome
  };
}

// Mappings onto sequences of this assembly: none -> 422 VARIANT_NOT_ON_ASSEMBLY, several -> 422 AMBIGUOUS_VARIANT_MAPPING.
async function assemblyMapping(ctx, id, mappings) {
  const usable = [];
  const seen = new Set();
  for (const m of mappings) {
    const len = await ctx.sequence.regionLength(ctx.fasta, m.seq_region_name);
    // an insertion after the last base has start = length + 1
    if (!isPosInt(len) || Math.min(m.start, m.end) < 1 || m.end > len || m.start > len + 1) continue;
    const key = [m.seq_region_name, m.start, m.end, m.allele_string].join('\t');
    if (seen.has(key)) continue;
    seen.add(key);
    usable.push({ region: m.seq_region_name, start: m.start, end: m.end, allele_string: m.allele_string, region_length: len });
  }
  if (usable.length === 0) {
    throw new PrimerHttpError(422, 'VARIANT_NOT_ON_ASSEMBLY', 'variant ' + id + ' has no mapping onto a sequence of this assembly', { id: id });
  }
  if (usable.length > 1) {
    throw new PrimerHttpError(422, 'AMBIGUOUS_VARIANT_MAPPING', 'variant ' + id + ' maps to ' + usable.length + ' places on this assembly', {
      id: id,
      mappings: usable.map(function (m) { return { region: m.region, start: m.start, end: m.end, allele_string: m.allele_string }; })
    });
  }
  return usable[0];
}

// §2.4 and §2.7 rule 9: one id to its entries, one per alt, with every alias of the same event merged. The lookup
// gives the mapping; the overlap chunk at the mapping gives the Ensembl record itself (with its short source code,
// §3.7) and every other id describing the same event. `flank` widens the one sequence read for a design, whose
// template and neighbour normalization then use the same read.
// -> {entries (those holding id), skipped, region, region_length, allele_string, genome}
async function resolveId(ctx, systemName, id, flank) {
  let lookup;
  try {
    lookup = await ctx.client.variation(ctx.species, id, { signal: ctx.signal });
  } catch (err) {
    if (err && err.code === 'UNKNOWN_VARIANT') {
      throw new PrimerHttpError(404, 'UNKNOWN_VARIANT', 'unknown variant ' + id + ' for ' + systemName, { id: id, system_name: systemName });
    }
    throw err;
  }
  const mapping = await assemblyMapping(ctx, id, lookup.mappings);
  const region = mapping.region;
  const lo = Math.min(mapping.start, mapping.end);
  const hi = Math.max(mapping.start, mapping.end);
  const fetched = await fetchChunks(ctx, region, mapping.region_length, lo - 1, hi + 1);
  let own = fetched.records.filter(function (r) { return r.id === id && spanOverlaps(r.start, r.end, lo - 1, hi + 1); });
  if (own.length === 0) {
    // Not in the overlap (a lookup-only id): the mapping itself, without a source code.
    own = [{ id: id, seq_region_name: region, start: mapping.start, end: mapping.end, alleles: mapping.allele_string.split('/'),
      source: null, consequence_type: lookup.consequence }];
  }
  const opts = normalizeOptions(ctx.v, region);
  const pad = flank + ctx.v.max_shift + WINDOW_MARGIN;
  const out = await withGenome(ctx, region, mapping.region_length, { start: lo - pad, end: hi + pad }, async function (genome) {
    // The requested records alone give each event's shift, which bounds where other descriptions of it can lie.
    let spanLo = lo - 1;
    let spanHi = hi + 1;
    normalize.recordsToEntries(own, genome, opts).entries.forEach(function (e) {
      spanLo = Math.min(spanLo, e.vcf.position - 1);
      spanHi = Math.max(spanHi, e.vcf.position + e.vcf.ref.length + e.shift);
    });
    const pool = spanLo < fetched.start || spanHi > fetched.end
      ? await fetchChunks(ctx, region, mapping.region_length, spanLo, spanHi) : fetched;
    const candidates = own.concat(pool.records.filter(function (r) { return spanOverlaps(r.start, r.end, spanLo, spanHi); }));
    const res = normalize.recordsToEntries(candidates, genome, Object.assign({}, opts, { requested_id: id, synonyms: lookup.synonyms }));
    return {
      entries: res.entries.filter(function (e) { return e.ids.indexOf(id) >= 0; }),
      skipped: skippedSummary(pool.skipped, spanLo, spanHi, res.skipped)
    };
  });
  return {
    entries: out.value.entries,
    skipped: out.value.skipped,
    region: region,
    region_length: mapping.region_length,
    allele_string: mapping.allele_string,
    genome: out.genome
  };
}

// ---- GET /primers/genomes --------------------------------------------------------------------------------------

// §3.2 and §2.2, pure: variation data is available when primers.variation is enabled, lists the genome's species
// segment and the assembly resolved with a FASTA. `variation` is never null: {available: false, source: null,
// release: null} without data.
function variationInfo(systemName, resolved, cfg) {
  cfg = cfg || require('../config').get();
  const available = variationEnabled(cfg) && speciesConfigured(cfg, systemName) !== null &&
    !!(resolved && resolved.fasta && resolved.fasta.dna);
  return { available: available, source: available ? SOURCE_NAME : null, release: available ? releaseOf(cfg) : null };
}

// hasVariation(system_name, deps) -> Promise<boolean>. Never touches the filesystem for a genome without a
// configured species; a genome that does not resolve (4xx) has no variation data.
async function hasVariation(systemName, deps) {
  deps = deps || {};
  const cfg = configOf(deps);
  if (!variationEnabled(cfg) || speciesConfigured(cfg, systemName) === null) return false;
  let resolved;
  try {
    resolved = await resolveAssembly(systemName, deps);
  } catch (err) {
    if (isPrimerHttpError(err) && err.status < 500) return false;
    throw err;
  }
  return variationInfo(systemName, resolved, cfg).available;
}

// ---- GET /primers/variants --------------------------------------------------------------------------------------

// Handler-side query rules (swagger checks the types first). Values may also arrive as query strings.
function listQuery(query, v) {
  query = query || {};
  const region = query.region;
  if (typeof region !== 'string' || region === '' || region.length > MAX_REGION_NAME) {
    throw invalidRequest('region', 'region must be a non-empty string of at most ' + MAX_REGION_NAME + ' characters');
  }
  const start = intValue(query.start);
  const end = intValue(query.end);
  if (start === null || end === null || start < 1 || end < start) {
    throw new PrimerHttpError(400, 'REGION_OUT_OF_BOUNDS', 'start and end must be integers with 1 <= start <= end',
      { region: region, start: start, end: end });
  }
  const length = end - start + 1;
  if (length > v.max_window) {
    throw new PrimerHttpError(400, 'VARIANT_WINDOW_TOO_LONG', 'the window is ' + length + ' bp; the limit is ' + v.max_window + ' bp',
      { length: length, max: v.max_window });
  }
  let limit = v.list_limit_default;
  if (query.limit !== undefined && query.limit !== null) {
    limit = intValue(query.limit);
    if (limit === null || limit < 1 || limit > v.list_limit_max) throw invalidRequest('limit', 'limit must be an integer from 1 to ' + v.list_limit_max);
  }
  let types = KINDS.slice();
  if (query.types !== undefined && query.types !== null) {
    const list = (Array.isArray(query.types) ? query.types : String(query.types).split(','))
      .map(function (t) { return String(t).trim(); }).filter(function (t) { return t !== ''; });
    if (list.some(function (t) { return KINDS.indexOf(t) < 0; })) {
      throw invalidRequest('types', 'types must be a comma-separated list of ' + KINDS.join(', '));
    }
    if (list.length) types = list;
  }
  let includeEms = true;
  if (query.include_ems !== undefined && query.include_ems !== null) {
    if (query.include_ems === false || query.include_ems === 'false') includeEms = false;
    else if (query.include_ems !== true && query.include_ems !== 'true') throw invalidRequest('include_ems', 'include_ems must be true or false');
  }
  return { system_name: query.system_name, region: region, start: start, end: end, limit: limit, types: types, include_ems: includeEms };
}

// listVariants({system_name, region, start, end, types, include_ems, limit}, deps) -> §2.3 body
// Errors: 400 REGION_OUT_OF_BOUNDS, VARIANT_WINDOW_TOO_LONG, INVALID_REQUEST; 404 UNKNOWN_GENOME, UNKNOWN_REGION;
// 422 NO_VARIATION_DATA, NO_SEQUENCE, AMBIGUOUS_ASSEMBLY; 503 VARIATION_SOURCE_UNAVAILABLE, FEATURE_DISABLED, MONGO_UNAVAILABLE.
async function listVariants(query, deps) {
  deps = deps || {};
  const cfg = configOf(deps);
  if (!variationEnabled(cfg)) throw featureDisabled();
  const v = variationCfg(cfg);
  const q = listQuery(query, v);
  const ctx = await variationContext(q.system_name, deps, cfg);
  const regionLength = await requireRegionLength(ctx.sequence, ctx.fasta, q.region);
  if (q.end > regionLength) throw outOfBounds(q.region, { start: q.start, end: q.end }, regionLength);

  const win = await windowEntries(ctx, q.region, regionLength, q.start, q.end, null);
  const matched = win.entries.filter(function (e) { return q.types.indexOf(e.kind) >= 0 && (q.include_ems || !e.ems); });
  const variants = matched.slice(0, q.limit);
  const warnings = [];
  pushSkipped(warnings, win.skipped);
  pushRefMismatches(warnings, matched);
  if (matched.length > variants.length) {
    warnings.push(warning('VARIANTS_TRUNCATED', 'showing ' + variants.length + ' of ' + matched.length + ' variants; narrow the window',
      { returned: variants.length, total: matched.length, limit: q.limit }));
  }
  return {
    system_name: q.system_name,
    region: q.region,
    start: q.start,
    end: q.end,
    source: { name: SOURCE_NAME, release: releaseOf(cfg) },
    total: matched.length,
    returned: variants.length,
    truncated: matched.length > variants.length,
    variants: variants,
    warnings: warnings
  };
}

// ---- GET /primers/variants/{variant_id} ----------------------------------------------------------------------------

// lookupVariant({system_name, variant_id}, deps) -> §2.4 body. Every entry of the id is returned, the
// non-designable ones with their issues (as in a listing).
// Errors: 400 INVALID_REQUEST; 404 UNKNOWN_VARIANT, UNKNOWN_GENOME; 422 NO_VARIATION_DATA, NO_SEQUENCE,
// AMBIGUOUS_VARIANT_MAPPING, VARIANT_NOT_ON_ASSEMBLY; 503 VARIATION_SOURCE_UNAVAILABLE, FEATURE_DISABLED.
async function lookupVariant(query, deps) {
  deps = deps || {};
  query = query || {};
  const cfg = configOf(deps);
  if (!variationEnabled(cfg)) throw featureDisabled();
  const id = query.variant_id;
  if (typeof id !== 'string' || !ID_RE.test(id)) throw invalidRequest('variant_id', 'variant_id must match ' + ID_RE.source);
  const ctx = await variationContext(query.system_name, deps, cfg);
  const site = await resolveId(ctx, query.system_name, id, 0);
  const warnings = [];
  site.entries.forEach(function (e) { if (e.ids.length > 1) warnings.push(duplicateIdsWarning(e)); });
  pushSkipped(warnings, site.skipped);
  pushRefMismatches(warnings, site.entries);
  return {
    requested_id: id,
    system_name: query.system_name,
    source: { name: SOURCE_NAME, release: releaseOf(cfg) },
    variants: site.entries,
    warnings: warnings
  };
}

// ---- genotyping design ------------------------------------------------------------------------------------------

function invalidVariant(reason, message) {
  return new PrimerHttpError(400, 'INVALID_VARIANT', message, { reason: reason });
}

function unsupportedAllele(allele) {
  return new PrimerHttpError(400, 'UNSUPPORTED_ALLELE', 'allele ' + allele + ' cannot be designed', { allele: allele });
}

// The entry, or the 400 its first blocking issue stands for (STAR_ALLELE does not block).
function assertDesignable(entry) {
  const issue = entry.issues.find(function (i) { return i.code !== 'STAR_ALLELE'; });
  if (!issue) return entry;
  const d = issue.details || {};
  switch (issue.code) {
    case 'REF_MISMATCH':
      throw new PrimerHttpError(400, 'REF_MISMATCH', 'the reference allele ' + d.given + ' does not match the genome ' +
        (String(d.genome).length === 1 ? 'base ' : 'bases ') + d.genome + ' at ' + d.region + ':' + d.position,
      { region: d.region, position: d.position, given: d.given, genome: d.genome });
    case 'UNSUPPORTED_ALLELE':
      throw unsupportedAllele(d.allele);
    case 'ALLELE_TOO_LONG':
      throw invalidVariant('allele_too_long', 'alleles are limited to ' + d.max + ' nt');
    default: // REPEAT_TOO_LONG
      throw new PrimerHttpError(400, 'VARIANT_TOO_REPETITIVE', issue.message,
        { region: entry.region, position: entry.vcf.position, shift: entry.shift, max: d.max });
  }
}

// §2.7 rules 10-12 and 14 over the entries of one id. alt matches an entry's minimal or VCF allele.
function chooseEntry(entries, id, alt, alleleString) {
  if (entries.length === 0) throw unsupportedAllele(alleleString);
  let alleles = entries[0].alleles;
  entries.forEach(function (e) { if (e.alleles.length > alleles.length) alleles = e.alleles; });
  if (alt !== undefined && alt !== null) {
    const a = String(alt).toUpperCase();
    if (a === '*' || a.indexOf('N') >= 0) throw unsupportedAllele(a);
    const match = entries.find(function (e) { return e.minimal.alt === a || e.vcf.alt === a; });
    if (!match) {
      throw new PrimerHttpError(400, 'ALT_NOT_AT_SITE', a + ' is not an alternative allele of ' + id, { id: id, alt: a, alleles: alleles.slice() });
    }
    return assertDesignable(match);
  }
  const designable = entries.filter(function (e) { return e.designable; });
  if (designable.length > 1) {
    const rank = function (e) {
      const i = alleles.indexOf(e.minimal.alt);
      return i < 0 ? alleles.length : i;
    };
    const alts = designable.slice().sort(function (x, y) { return rank(x) - rank(y); }).map(function (e) { return e.minimal.alt; });
    throw new PrimerHttpError(400, 'ALT_REQUIRED', id + ' has more than one alternative allele; choose one with variant.alt', { id: id, alts: alts });
  }
  return assertDesignable(designable[0] || entries[0]);
}

// The entry with requested_id placed after key, as in the §2.9 variant block.
function withRequestedId(entry, requestedId) {
  const out = { key: entry.key, requested_id: requestedId };
  Object.keys(entry).forEach(function (k) { if (k !== 'key') out[k] = entry[k]; });
  return out;
}

// resolveDesignVariant(request, deps): §2.7 rules 7-14, all before any design semaphore slot is taken (§4.2).
// Rules 2-3 are re-checked here (genotyping/request.js owns rules 1-6).
//   request: {system_name, variant: {id [, alt]} | {region, position, ref, alt}, flank}
//     flank: the §4.3 template flank F (default genotyping.template_flank); the one sequence read spans the variant
//     +- (flank + max_shift + 250), so the template and the neighbour normalization can reuse it.
//   -> {variant: the designable, REF-verified entry with requested_id (the id, or null for manual input) after key,
//         and without submission_sequence; for manual input ids/records stay empty until neighboursFor fills them,
//       resolved (assemblies.resolve), fasta, region, region_length,
//       species: the Ensembl species segment, or null when the genome has no (enabled) variation data,
//       has_variation, genome: the normalize.sequenceWindow accessor of the read,
//       warnings: [DUPLICATE_VARIANT_IDS] for an id whose event has aliases}
// Errors: 400 INVALID_VARIANT, UNSUPPORTED_ALLELE, ALT_REQUIRED, ALT_NOT_AT_SITE, REF_MISMATCH, VARIANT_TOO_REPETITIVE,
//   REGION_OUT_OF_BOUNDS; 404 UNKNOWN_GENOME, UNKNOWN_REGION, UNKNOWN_VARIANT; 422 NO_SEQUENCE, AMBIGUOUS_ASSEMBLY,
//   NO_VARIATION_DATA, AMBIGUOUS_VARIANT_MAPPING, VARIANT_NOT_ON_ASSEMBLY; 503 VARIATION_SOURCE_UNAVAILABLE,
//   FEATURE_DISABLED (an id while primers.variation.enabled is false), MONGO_UNAVAILABLE.
async function resolveDesignVariant(request, deps) {
  deps = deps || {};
  request = request || {};
  const cfg = configOf(deps);
  const v = variationCfg(cfg);
  const input = request.variant !== null && typeof request.variant === 'object' && !Array.isArray(request.variant) ? request.variant : {};
  const byId = input.id !== undefined && input.id !== null;
  let manual = null;
  if (byId) {
    const mixed = ['region', 'position', 'ref'].some(function (k) { return input[k] !== undefined && input[k] !== null; });
    if (typeof input.id !== 'string' || !ID_RE.test(input.id) || mixed) {
      throw invalidVariant('id_or_manual', 'give either variant.id (with an optional alt) or region, position, ref and alt');
    }
  } else {
    manual = normalize.parseManual(input, v);
  }
  const flank = Number.isSafeInteger(request.flank) && request.flank >= 0 ? request.flank
    : (cfg.genotyping && Number.isSafeInteger(cfg.genotyping.template_flank) ? cfg.genotyping.template_flank : DEFAULT_TEMPLATE_FLANK);

  const systemName = request.system_name;
  const resolved = await resolveAssembly(systemName, deps);
  const fasta = resolved && resolved.fasta && resolved.fasta.dna;
  if (!fasta) throw noSequence(systemName);
  const configured = speciesConfigured(cfg, systemName);
  const species = configured !== null && variationEnabled(cfg) ? configured : null;
  const ctx = makeContext(cfg, deps, resolved, species, fasta);

  if (byId) {
    if (configured === null) throw noVariationData(systemName);
    if (species === null) throw featureDisabled();
    const site = await resolveId(ctx, systemName, input.id, flank);
    const entry = chooseEntry(site.entries, input.id, input.alt, site.allele_string);
    return {
      variant: withRequestedId(entry, input.id),
      resolved: resolved,
      fasta: fasta,
      region: site.region,
      region_length: site.region_length,
      species: species,
      has_variation: true,
      genome: site.genome,
      warnings: entry.ids.length > 1 ? [duplicateIdsWarning(entry)] : []
    };
  }

  const region = manual.seq_region_name;
  const regionLength = await requireRegionLength(ctx.sequence, fasta, region);
  const beyond = { position: manual.start };
  if (manual.start > regionLength || manual.end > regionLength) throw outOfBounds(region, beyond, regionLength);
  const lo = Math.min(manual.start, manual.end);
  const hi = Math.max(manual.start, manual.end);
  const pad = flank + v.max_shift + WINDOW_MARGIN;
  const out = await withGenome(ctx, region, regionLength, { start: lo - pad, end: hi + pad }, function (genome) {
    return normalize.recordsToEntries([manual], genome, normalizeOptions(v, region));
  });
  // No entry: the event needs an anchor base before position 1.
  if (out.value.entries.length === 0) throw outOfBounds(region, beyond, regionLength);
  const entry = assertDesignable(out.value.entries[0]);
  return {
    variant: withRequestedId(entry, null),
    resolved: resolved,
    fasta: fasta,
    region: region,
    region_length: regionLength,
    species: species,
    has_variation: species !== null,
    genome: out.genome,
    warnings: []
  };
}

// A manual variant takes its Ensembl identity from the window entry of the same key (§3.9).
function mergeTarget(variant, target) {
  const out = Object.assign({}, variant);
  ['ids', 'records', 'ems', 'consequence', 'alleles', 'multiallelic'].forEach(function (k) { out[k] = target[k]; });
  const star = target.issues.filter(function (i) { return i.code === 'STAR_ALLELE'; });
  if (star.length && !variant.issues.some(function (i) { return i.code === 'STAR_ALLELE'; })) out.issues = variant.issues.concat(star);
  return out;
}

// neighboursFor(request, deps): §4.15 known variants around a design, from the same cached chunks as every other
// call (never one fetch per variant), normalized over the design's own sequence read.
//   request: {system_name, variant (resolveDesignVariant().variant), window: {start, end} (the template window),
//             resolved, region_length, genome (from resolveDesignVariant; looked up again when absent),
//             strict (default: variant.requested_id is a string, i.e. a design by id)}
//   -> {data: 'ensembl' | 'none' | 'unavailable',
//       variant: the request variant; for manual input (requested_id null) with a window entry of the same key,
//                its ids, records, ems, consequence, alleles, multiallelic and STAR_ALLELE issue come from that entry,
//       window,
//       entries: the neighbour entries (minimal span overlapping the window; the target key and entries sharing a
//                target id excluded; EMS included), sorted by position then key,
//       target: the window entry with the target key, or null,
//       summary: PrimerGenotypingNeighbourSummary {data, window, variants, non_ems, ems, dense_non_ems},
//       skipped: {count, reasons} of malformed records (not a design warning),
//       warnings: [NO_VARIATION_DATA {system_name}] | [NEIGHBOURS_UNAVAILABLE {reason}] |
//                 DUPLICATE_VARIANT_IDS (manual input merged into an aliased event) and DENSE_NEIGHBOURS {count, window},
//       genome: the sequence accessor used (possibly widened), or the request's genome}
// A genome without (enabled) variation data gives data 'none'. Ensembl unavailable gives data 'unavailable' with
// zero counts, unless strict, which rethrows the 503 (§3.8).
async function neighboursFor(request, deps) {
  deps = deps || {};
  request = request || {};
  const cfg = configOf(deps);
  const g = cfg.genotyping || {};
  const variant = request.variant;
  const systemName = request.system_name;
  const w = request.window || {};
  if (!variant || typeof variant.key !== 'string' || !isPosInt(w.start) || !isPosInt(w.end) || w.end < w.start) {
    throw new TypeError('neighboursFor: variant and window {start, end} are required');
  }
  const window = { start: w.start, end: w.end };
  const strict = request.strict !== undefined ? request.strict === true : typeof variant.requested_id === 'string';
  const denseWindow = Number.isSafeInteger(g.dense_window) ? g.dense_window : DEFAULT_DENSE_WINDOW;
  const denseCount = Number.isSafeInteger(g.dense_count) ? g.dense_count : DEFAULT_DENSE_COUNT;
  const summary = function (data, counts) {
    return Object.assign({ data: data, window: { start: window.start, end: window.end } }, counts || { variants: 0, non_ems: 0, ems: 0, dense_non_ems: 0 });
  };
  const nothing = function (data, warn) {
    return { data: data, variant: variant, window: window, entries: [], target: null, summary: summary(data), skipped: { count: 0, reasons: {} },
      warnings: [warn], genome: request.genome || null };
  };

  const enabled = variationEnabled(cfg);
  const species = enabled ? speciesConfigured(cfg, systemName) : null;
  if (species === null) {
    // Switched off and "genome without data" share the warning code (§2.15 has no separate one); only the message
    // says which, so a disabled server does not claim the genome lacks data.
    const message = enabled
      ? systemName + ' has no known-variant data; neighbouring variants were not screened'
      : 'known-variant lookups are disabled on this server; neighbouring variants were not screened';
    return nothing('none', warning('NO_VARIATION_DATA', message, { system_name: systemName }));
  }
  const resolved = request.resolved || await resolveAssembly(systemName, deps);
  const fasta = resolved && resolved.fasta && resolved.fasta.dna;
  if (!fasta) throw noSequence(systemName);
  const ctx = makeContext(cfg, deps, resolved, species, fasta);
  const region = variant.region;
  const regionLength = isPosInt(request.region_length) ? request.region_length : await requireRegionLength(ctx.sequence, fasta, region);

  let win;
  try {
    win = await windowEntries(ctx, region, regionLength, window.start, window.end, request.genome || null);
  } catch (err) {
    if (strict || !isUnavailable(err)) throw err;
    return nothing('unavailable', warning('NEIGHBOURS_UNAVAILABLE',
      'known variants near the primers could not be fetched from Ensembl and were not screened', { reason: err.details.reason }));
  }

  const target = win.entries.find(function (e) { return e.key === variant.key; }) || null;
  const merged = target && typeof variant.requested_id !== 'string' ? mergeTarget(variant, target) : variant;
  const targetIds = new Set(merged.ids || []);
  const entries = win.entries.filter(function (e) {
    return e.key !== variant.key && !e.ids.some(function (id) { return targetIds.has(id); });
  });
  const nonEms = entries.filter(function (e) { return !e.ems; });
  const a = variant.vcf.position - denseWindow;
  const b = variant.vcf.position + variant.vcf.ref.length - 1 + denseWindow;
  const dense = nonEms.filter(function (e) { return normalize.inWindow(e.minimal, a, b); }).length;

  const warnings = [];
  if (merged !== variant && merged.ids.length > 1) warnings.push(duplicateIdsWarning(merged));
  if (dense > denseCount) {
    warnings.push(warning('DENSE_NEIGHBOURS', dense + ' known non-EMS variants lie within ' + denseWindow + ' bp of the variant',
      { count: dense, window: denseWindow }));
  }
  return {
    data: 'ensembl',
    variant: merged,
    window: window,
    entries: entries,
    target: target,
    summary: summary('ensembl', { variants: entries.length, non_ems: nonEms.length, ems: entries.length - nonEms.length, dense_non_ems: dense }),
    skipped: win.skipped,
    warnings: warnings,
    genome: win.genome
  };
}

function _resetClient() {
  sharedClient = null;
}

module.exports = {
  hasVariation,
  variationInfo,
  listVariants,
  lookupVariant,
  resolveDesignVariant,
  neighboursFor,
  VariationError,
  KINDS,
  ID_PATTERN: ID_RE.source,
  WINDOW_MARGIN,
  _resetClient,
  _internal: { chunkStarts, fetchChunks, withGenome, skippedSummary, listQuery, chooseEntry, assertDesignable, mergeTarget }
};
