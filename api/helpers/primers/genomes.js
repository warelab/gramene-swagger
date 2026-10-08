'use strict';

// Genome catalog for the /primers endpoints (spec §A.2.2, §A.4.1 item 1).
//
// Catalog = mongo `maps` docs of type 'genome' joined with the pruned NCBI `taxonomy`.
// Species of a genome = the first entry of its taxon's `ancestors` (ordered self -> root) whose
// taxonomy doc has rank 'species'. Never taxon_id/1000: rice genomes hang under 39946, 1736656 and
// 1736659, which all descend from 4530 (Oryza sativa); 4558006 -> 4558, not 'Sorghum pan'.
//
// The catalog is cached per mongo handle for 10 minutes. A failed reload keeps serving the previous
// catalog (maps only change at release time) and is retried at most every 30 s; with no catalog at
// all the caller gets 503 MONGO_UNAVAILABLE. Mongo lookups are equality / $in on ids only.

const { PrimerHttpError, isPrimerHttpError, redactPaths } = require('./errors');

const SYSTEM_NAME_RE = /^[a-z0-9_]+$/;
const SYSTEM_NAME_MAX_LENGTH = 128;
const CATALOG_TTL_MS = 10 * 60 * 1000;
const CATALOG_RETRY_MS = 30 * 1000;
const MONGO_TIMEOUT_MS = 10 * 1000;
const MONGO_RETRY_AFTER_S = 30;

const MAP_FIELDS = { _id: 1, system_name: 1, display_name: 1, taxon_id: 1, type: 1, regions: 1 };
const TAXON_FIELDS = { _id: 1, name: 1, rank: 1, ancestors: 1 };

function isValidSystemName(name) {
  return typeof name === 'string' && name.length > 0 && name.length <= SYSTEM_NAME_MAX_LENGTH &&
    SYSTEM_NAME_RE.test(name);
}

function unknownGenomeError(name) {
  // Only echo names that passed the pattern: they cannot carry markup, paths or control characters.
  if (isValidSystemName(name)) {
    return new PrimerHttpError(404, 'UNKNOWN_GENOME', 'unknown genome: ' + name, { system_name: name });
  }
  return new PrimerHttpError(404, 'UNKNOWN_GENOME', 'unknown genome', {});
}

function logSafe(log, level, msg) {
  try {
    const fn = (log && (log[level] || log.error || log.log)) || function () {};
    fn.call(log, msg);
  } catch (e) { /* logging must never throw */ }
}

function compareText(a, b) {
  return a < b ? -1 : a > b ? 1 : 0;
}

function compareGenomes(a, b) {
  return a.display_name.localeCompare(b.display_name, 'en', { numeric: true, sensitivity: 'base' }) ||
    compareText(a.system_name, b.system_name);
}

// ---- pure core -------------------------------------------------------------------------------

function toTaxonIndex(taxonomy) {
  if (taxonomy instanceof Map) return taxonomy;
  const index = new Map();
  (taxonomy || []).forEach(function (doc) {
    if (doc && doc._id !== undefined && doc._id !== null && Number.isFinite(Number(doc._id))) {
      index.set(Number(doc._id), doc);
    }
  });
  return index;
}

// speciesOf(taxonId, taxonomy) -> {taxon_id, name} | null
// taxonomy: array of taxonomy docs {_id, rank, name, ancestors} or a Map(id -> doc).
function speciesOf(taxonId, taxonomy) {
  const index = toTaxonIndex(taxonomy);
  const id = Number(taxonId);
  const self = index.get(id);
  if (!self) return null;
  const chain = Array.isArray(self.ancestors) ? self.ancestors.slice() : [];
  if (Number(chain[0]) !== id) chain.unshift(id);
  for (let i = 0; i < chain.length; i++) {
    const doc = index.get(Number(chain[i]));
    if (doc && doc.rank === 'species') {
      return Object.freeze({ taxon_id: Number(doc._id), name: typeof doc.name === 'string' ? doc.name : String(doc._id) });
    }
  }
  return null;
}

function normalizeRegions(regions) {
  const names = regions && Array.isArray(regions.names) ? regions.names : [];
  const lengths = regions && Array.isArray(regions.lengths) ? regions.lengths : [];
  const n = Math.min(names.length, lengths.length);
  return Object.freeze({
    names: Object.freeze(names.slice(0, n).map(String)),
    lengths: Object.freeze(lengths.slice(0, n).map(Number))
  });
}

function speciesKey(genome) {
  return genome.species ? 'species:' + genome.species.taxon_id : 'genome:' + genome.system_name;
}

// buildCatalog(mapDocs, taxonomyDocs[, {loaded_at}]) -> frozen catalog
//   {loaded_at, genomes[] (sorted by display_name), bySystemName: Map, bySpecies: Map(key -> genome[])}
// genome: {system_name, display_name, taxon_id, map_id, regions {names[], lengths[]}, species {taxon_id, name}|null}
// Docs that are not type 'genome' or whose system_name fails the pattern are skipped; for a duplicated
// system_name the doc with the smallest _id wins.
function buildCatalog(mapDocs, taxonomyDocs, opts) {
  opts = opts || {};
  const index = toTaxonIndex(taxonomyDocs);
  const docs = (mapDocs || []).filter(function (doc) {
    return doc && (doc.type === undefined || doc.type === 'genome') && doc._id !== undefined && doc._id !== null &&
      isValidSystemName(doc.system_name);
  }).sort(function (a, b) { return compareText(String(a._id), String(b._id)); });

  const bySystemName = new Map();
  docs.forEach(function (doc) {
    if (bySystemName.has(doc.system_name)) return;
    const taxonId = Number(doc.taxon_id);
    const hasTaxon = doc.taxon_id !== undefined && doc.taxon_id !== null && Number.isFinite(taxonId);
    bySystemName.set(doc.system_name, Object.freeze({
      system_name: doc.system_name,
      display_name: typeof doc.display_name === 'string' && doc.display_name ? doc.display_name : doc.system_name,
      taxon_id: hasTaxon ? taxonId : null,
      map_id: String(doc._id),
      regions: normalizeRegions(doc.regions),
      species: hasTaxon ? speciesOf(taxonId, index) : null
    }));
  });

  const list = Array.from(bySystemName.values()).sort(compareGenomes);
  const bySpecies = new Map();
  list.forEach(function (genome) {
    const key = speciesKey(genome);
    if (!bySpecies.has(key)) bySpecies.set(key, []);
    bySpecies.get(key).push(genome);
  });
  bySpecies.forEach(function (members) { Object.freeze(members); });

  return Object.freeze({
    loaded_at: opts.loaded_at || new Date().toISOString(),
    genomes: Object.freeze(list),
    bySystemName: bySystemName,
    bySpecies: bySpecies
  });
}

// ---- mongo loading and cache ------------------------------------------------------------------

function uniqueNumbers(values) {
  const out = new Set();
  values.forEach(function (v) {
    if (v === undefined || v === null || v === '') return;
    const n = Number(v);
    if (Number.isFinite(n)) out.add(n);
  });
  return Array.from(out);
}

function withTimeout(promise, ms, what) {
  let timer;
  const timeout = new Promise(function (resolve, reject) {
    timer = setTimeout(function () { reject(new Error(what + ' timed out after ' + ms + ' ms')); }, ms);
    if (timer.unref) timer.unref();
  });
  return Promise.race([promise, timeout]).finally(function () { clearTimeout(timer); });
}

function defaultMongo() {
  return require('gramene-mongodb-config');
}

// loadCatalog(deps) -> catalog (uncached). deps: {mongo, mongo_timeout_ms, log}
async function loadCatalog(deps) {
  deps = deps || {};
  const mongo = deps.mongo || defaultMongo();
  const log = deps.log || console;
  const timeoutMs = deps.mongo_timeout_ms || MONGO_TIMEOUT_MS;
  try {
    return await withTimeout((async function () {
      // gramene-mongodb-config resolves mongoCollection() to undefined when the connection failed.
      const maps = await mongo.maps.mongoCollection();
      const taxonomy = await mongo.taxonomy.mongoCollection();
      if (!maps || !taxonomy) throw new Error('mongo collections maps/taxonomy are unavailable');
      const mapDocs = await maps.find({ type: 'genome' }, { fields: Object.assign({}, MAP_FIELDS) }).toArray();
      if (!Array.isArray(mapDocs) || mapDocs.length === 0) throw new Error('maps has no genome documents');
      const leafIds = uniqueNumbers(mapDocs.map(function (d) { return d && d.taxon_id; }));
      const leaves = leafIds.length
        ? await taxonomy.find({ _id: { $in: leafIds } }, { fields: Object.assign({}, TAXON_FIELDS) }).toArray()
        : [];
      const have = new Set(leaves.map(function (d) { return Number(d._id); }));
      const ancestorIds = uniqueNumbers([].concat.apply([], leaves.map(function (d) {
        return Array.isArray(d.ancestors) ? d.ancestors : [];
      }))).filter(function (id) { return !have.has(id); });
      const ancestors = ancestorIds.length
        ? await taxonomy.find({ _id: { $in: ancestorIds } }, { fields: Object.assign({}, TAXON_FIELDS) }).toArray()
        : [];
      return buildCatalog(mapDocs, leaves.concat(ancestors));
    })(), timeoutMs, 'genome catalog query');
  } catch (err) {
    if (isPrimerHttpError(err)) throw err;
    logSafe(log, 'error', 'primers: genome catalog unavailable: ' + (err && err.message ? err.message : String(err)));
    throw new PrimerHttpError(503, 'MONGO_UNAVAILABLE', 'the genome catalog is temporarily unavailable',
      { retry_after_s: MONGO_RETRY_AFTER_S });
  }
}

let caches = new WeakMap();

function cacheFor(mongo) {
  let cache = caches.get(mongo);
  if (!cache) {
    cache = { catalog: null, loadedAt: 0, lastAttempt: -Infinity, inflight: null };
    caches.set(mongo, cache);
  }
  return cache;
}

function startLoad(cache, mongo, deps) {
  const now = deps.now || Date.now;
  cache.lastAttempt = now();
  const p = loadCatalog(Object.assign({}, deps, { mongo: mongo })).then(function (catalog) {
    cache.catalog = catalog;
    cache.loadedAt = now();
    return catalog;
  });
  cache.inflight = p;
  const clear = function () { if (cache.inflight === p) cache.inflight = null; };
  p.then(clear, clear);
  return p;
}

// getCatalog(deps) -> Promise<catalog>
// deps: {catalog (use as-is, no mongo), mongo, now, catalog_ttl_ms, catalog_retry_ms, mongo_timeout_ms, log}
function getCatalog(deps) {
  deps = deps || {};
  if (deps.catalog) return Promise.resolve(deps.catalog);
  const mongo = deps.mongo || defaultMongo();
  const now = (deps.now || Date.now)();
  const ttl = deps.catalog_ttl_ms !== undefined ? deps.catalog_ttl_ms : CATALOG_TTL_MS;
  const retry = deps.catalog_retry_ms !== undefined ? deps.catalog_retry_ms : CATALOG_RETRY_MS;
  const cache = cacheFor(mongo);
  if (cache.catalog && now - cache.loadedAt < ttl) return Promise.resolve(cache.catalog);
  if (cache.catalog) {
    // Stale: answer from the previous catalog and refresh in the background.
    if (!cache.inflight && now - cache.lastAttempt >= retry) {
      const log = deps.log || console;
      startLoad(cache, mongo, deps).catch(function (err) {
        logSafe(log, 'warn', 'primers: genome catalog reload failed (' + (err && err.code) + '); serving the previous catalog');
      });
    }
    return Promise.resolve(cache.catalog);
  }
  return cache.inflight || startLoad(cache, mongo, deps);
}

function clearCache() {
  caches = new WeakMap();
}

// ---- lookups ------------------------------------------------------------------------------------

// byName(system_name, deps) -> genome | null
async function byName(systemName, deps) {
  if (!isValidSystemName(systemName)) return null;
  const catalog = await getCatalog(deps);
  return catalog.bySystemName.get(systemName) || null;
}

// requireGenome(system_name, deps) -> genome, else 404 UNKNOWN_GENOME
async function requireGenome(systemName, deps) {
  const genome = await byName(systemName, deps);
  if (!genome) throw unknownGenomeError(systemName);
  return genome;
}

// sameSpecies(system_name, deps) -> genome[] of the same species, INCLUDING the genome itself,
// sorted by display_name. Unknown -> 404 UNKNOWN_GENOME. A genome without a species-rank
// ancestor is its own group.
async function sameSpecies(systemName, deps) {
  if (!isValidSystemName(systemName)) throw unknownGenomeError(systemName);
  const catalog = await getCatalog(deps);
  const genome = catalog.bySystemName.get(systemName);
  if (!genome) throw unknownGenomeError(systemName);
  return catalog.bySpecies.get(speciesKey(genome)) || Object.freeze([genome]);
}

// ---- GET /primers/genomes -----------------------------------------------------------------------

// variationFor(system_name, resolved) -> {available, source, release} (variation/index.js variationInfo)
function genomeEntry(genome, outcome, isQuery, log, variationFor) {
  const entry = {
    system_name: genome.system_name,
    display_name: genome.display_name,
    taxon_id: genome.taxon_id,
    map_id: genome.map_id,
    is_query: isQuery,
    has_sequence: false,
    has_blastdb: false,
    has_cdna_blastdb: false,
    has_variation: false,
    repeat_masking: 'absent',
    total_bases: null,
    warnings: []
  };
  const resolved = outcome && outcome.resolved;
  if (resolved) {
    entry.has_sequence = !!(resolved.fasta && resolved.fasta.dna);
    entry.has_blastdb = !!(resolved.blastdb && resolved.blastdb.dna);
    entry.has_cdna_blastdb = !!(resolved.blastdb && resolved.blastdb.cdna);
    entry.has_variation = variationFor(genome.system_name, resolved).available;
    entry.repeat_masking = resolved.repeat_masking || 'absent';
    entry.total_bases = typeof resolved.total_bases === 'number' ? resolved.total_bases : null;
    entry.warnings = (resolved.warnings || []).map(function (w) {
      return { code: w.code, message: redactPaths(w.message) };
    });
    return entry;
  }
  const err = outcome && outcome.error;
  if (isPrimerHttpError(err) && err.status < 500) {
    entry.warnings = [{ code: err.code, message: redactPaths(err.message) }];
  } else {
    logSafe(log, 'error', 'primers: could not resolve ' + genome.system_name + ': ' + (err && err.stack ? err.stack : String(err)));
    entry.warnings = [{
      code: isPrimerHttpError(err) ? err.code : 'INTERNAL',
      message: 'the assembly files for this genome could not be inspected'
    }];
  }
  return entry;
}

// genomesResponse(system_name, deps) -> A.2.2 body, plus the genotyping spec §2.2 variation fields
//   {system_name, species {taxon_id, name}|null, variation {available, source 'ensembl'|null, release|null},
//    counts {total, with_blastdb, with_cdna_blastdb},
//    genomes [{system_name, display_name, taxon_id, map_id, is_query, has_sequence, has_blastdb,
//              has_cdna_blastdb, has_variation, repeat_masking, total_bases, warnings[{code, message}]}]}
// Query genome first, then the rest of its species by display_name. Never contains filesystem paths.
// `variation` describes the query genome and is never null; has_variation needs primers.variation to be enabled,
// to list the genome's species and the assembly to resolve with a FASTA (§3.2). No Ensembl call is made.
// deps: everything getCatalog/assemblies.resolveMany accept (cfg included), plus {assemblies} to stub the resolver.
async function genomesResponse(systemName, deps) {
  deps = deps || {};
  if (!isValidSystemName(systemName)) throw unknownGenomeError(systemName);
  const catalog = await getCatalog(deps);
  const query = catalog.bySystemName.get(systemName);
  if (!query) throw unknownGenomeError(systemName);
  const members = catalog.bySpecies.get(speciesKey(query)) || [query];
  const ordered = [query].concat(members.filter(function (g) { return g !== query; }));
  const assemblies = deps.assemblies || require('./assemblies');
  const outcomes = await assemblies.resolveMany(ordered.map(function (g) { return g.system_name; }),
    Object.assign({}, deps, { catalog: catalog }));
  const log = deps.log || console;
  const cfg = deps.cfg || require('./config').get();
  const variation = require('./variation'); // lazily: variation/index.js reaches assemblies.js, which requires this module
  const variationFor = function (name, resolved) { return variation.variationInfo(name, resolved, cfg); };
  const list = ordered.map(function (genome, i) { return genomeEntry(genome, outcomes[i], genome === query, log, variationFor); });
  return {
    system_name: systemName,
    species: query.species ? { taxon_id: query.species.taxon_id, name: query.species.name } : null,
    variation: variationFor(systemName, outcomes[0] && outcomes[0].resolved),
    counts: {
      total: list.length,
      with_blastdb: list.filter(function (g) { return g.has_blastdb; }).length,
      with_cdna_blastdb: list.filter(function (g) { return g.has_cdna_blastdb; }).length
    },
    genomes: list
  };
}

module.exports = {
  SYSTEM_NAME_RE,
  SYSTEM_NAME_MAX_LENGTH,
  CATALOG_TTL_MS,
  isValidSystemName,
  unknownGenomeError,
  speciesOf,
  speciesKey,
  compareGenomes,
  buildCatalog,
  loadCatalog,
  getCatalog,
  clearCache,
  byName,
  requireGenome,
  sameSpecies,
  genomesResponse
};
