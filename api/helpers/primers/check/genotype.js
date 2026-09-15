'use strict';

// The genotyping block of POST /primers/check (genotyping spec §2.11, §5.1-§5.3, §5.6). One module for the API and the worker.
//
// Submit side, called by check/normalize.js in the §5.1 order:
//   validateShape(genotyping)                         step 1: the PrimerCheckGenotyping definition   400 INVALID_REQUEST {field}
//   validateLinks(genotyping, pairs)                  step 4: §5.2 part A, pure -> links              400 GENOTYPING_SET_INVALID
//   validateSets(genotyping, pairs, reference, deps)  step 6: one FASTA read of windowFor(), then prepare()
//   prepare(genotyping, pairs, genome, opts)          pure: §5.2 part B and the per-job values of §5.6 over that window
//   requestBlock(prepared)                            step 8: request.genotyping, with the client's keys only (§5.3)
// Step 6 errors: 400 REF_MISMATCH {region, position, given, genome}, VARIANT_TOO_REPETITIVE {region, position, shift, max},
// GENOTYPING_SET_INVALID {set_id, reason, pair_ids, primer?, sequence?, mm_pos?}, REGION_OUT_OF_BOUNDS
// {region, start, end, length}; 404 UNKNOWN_REGION; 422 NO_SEQUENCE.
//
// Worker side (check/run.js): milestone M8 adds emptyResults, callGenome, semiGlobal, mergeCopies, gapCompressedIdentity,
// megablastCopies, predict, writeResults and unavailableEntry here. They start from prepare()'s object, which normalize stores
// as resolved.genotyping (the job doc only, never returned to clients); validateSets() with
// deps.sequence = {fetch: fetchWindow, regionLength} re-derives the same object from the stored request and pairs.
//
// Orientation and the allele side come from the shared primer, never from the client, and are never stored in the request:
// a shared right primer makes the left primers allele-specific (forward), a shared left primer the right primers (reverse).
// Coordinates are 1-based inclusive on the reference plus strand. An ALT-haplotype coordinate numbers the ALT haplotype like
// the reference up to the variant (variation/normalize.js altToRef).
//
// prepare() -> {
//   algorithm_version: 'g1',
//   variant: {key, region, position, ref, alt (left-aligned, uppercase), region_length, shift, tract {start, end} | null,
//             zone {start, end}, discriminating {forward, reverse}, flank (K of §5.6),
//             core {start, end, ref, alt}, haplotypes {start, end, ref, alt}}   ref/alt are null when the span runs past the read
//   window: {region, start, end}                                                the reference span read at submit
//   sets: [{id, ref_pair, alt_pair, orientation, primers {as_ref, as_alt, common}, expected {region, start, end},
//           footprints {as_ref, common: {start, end, strand}}, deliberate_mismatch {as_ref, as_alt: 2 | 3 | null},
//           deliberate_mismatch_positions (as_ref's then as_alt's, nulls left out), in_shift_tract {as_ref, as_alt},
//           segment {start, end, sequence | null}}]                             segment: R_s of §5.6, expected ± 50 bp
// }

const { PrimerHttpError } = require('../errors');
const variation = require('../variation/normalize');
const realign = require('./realign');

// The caller and prediction version. Genotyping jobs hash with ALGORITHM_VERSION + '+' + this ('2+g1', §5.3), and results
// report it as results.genotyping.algorithm_version, so changing the caller invalidates genotyping jobs alone.
const GENOTYPING_VERSION = 'g1';
// Only these modes honour expected product locations (specificity.js), which the sets rely on (§5.1 step 3).
const MODES = Object.freeze(['gene', 'region']);
const MAX_SETS = 5;
const KEYS = Object.freeze(['variant', 'sets']);
const VARIANT_KEYS = Object.freeze(['region', 'position', 'ref', 'alt']);
const SET_KEYS = Object.freeze(['id', 'ref_pair', 'alt_pair']);
const REGION_MAX = 255;
const SET_ID_RE = /^[A-Za-z0-9_.:-]{1,64}$/;
const ALLELE_RE = /^[ACGTacgt]{1,50}$/;
// The submit-time read (§3.6, §5.2 part B): [position - 1,700, position + len(ref) + 1,700], clamped to the region.
const WINDOW_PAD = 1700;
// three_prime_window of the part B alignments; only mm_pos, gaps and the ends are read.
const THREE_PRIME_WINDOW = 5;
// An allele-specific primer may carry one undeclared extra mismatch this far from its 3' end: a deliberate mismatch (§5.2).
const DELIBERATE_DISTANCES = Object.freeze([2, 3]);
// GENOTYPING_SET_INVALID reasons (§2.14), in the order the checks run.
const REASONS = Object.freeze(['unknown_pair', 'same_pair', 'pair_reused', 'expected_required', 'expected_differs',
  'no_shared_common', 'not_at_variant', 'alleles_swapped', 'too_many_edits', 'common_in_zone', 'expected_mismatch']);
// Config fallbacks (primers.variation.max_shift, primers.check.genotype_flank_min / genotype_amplicon_pad, §3.1).
const DEFAULTS = Object.freeze({ max_shift: 1000, genotype_flank_min: 15, genotype_amplicon_pad: 50 });

const ROLES = Object.freeze({
  as_ref: Object.freeze({ pair: 'ref_pair', own: 'ref', other: 'alt', label: 'REF-specific primer' }),
  as_alt: Object.freeze({ pair: 'alt_pair', own: 'alt', other: 'ref', label: 'ALT-specific primer' })
});
const ALLELE_NAME = Object.freeze({ ref: 'REF', alt: 'ALT' });

function isPlainObject(v) {
  return v !== null && typeof v === 'object' && !Array.isArray(v) &&
    (Object.getPrototypeOf(v) === Object.prototype || Object.getPrototypeOf(v) === null);
}

function invalid(message, details) {
  return new PrimerHttpError(400, 'INVALID_REQUEST', message, details || {});
}

function checkKeys(obj, allowed, where) {
  Object.keys(obj).forEach(function (k) {
    if (allowed.indexOf(k) < 0) throw invalid('unknown field ' + where + k, { field: where + k });
  });
}

function setInvalid(setId, reason, pairIds, message, extra) {
  return new PrimerHttpError(400, 'GENOTYPING_SET_INVALID', 'set ' + setId + ': ' + message,
    Object.assign({ set_id: setId, reason: reason, pair_ids: pairIds }, extra || {}));
}

function count(v, fallback) {
  return Number.isSafeInteger(v) && v >= 0 ? v : fallback;
}

// ---- step 1: shape (mirrors PrimerCheckGenotyping; swagger normally rejects first) --------------------------------------

// Beyond the definition: ref and alt must differ, and set ids must be unique (results are keyed by set id).
function validateShape(genotyping) {
  const g = genotyping;
  if (!isPlainObject(g)) throw invalid('genotyping must be an object with variant and sets', { field: 'genotyping' });
  checkKeys(g, KEYS, 'genotyping.');
  const v = g.variant;
  if (!isPlainObject(v)) {
    throw invalid('genotyping.variant must be an object with region, position, ref and alt', { field: 'genotyping.variant' });
  }
  checkKeys(v, VARIANT_KEYS, 'genotyping.variant.');
  if (typeof v.region !== 'string' || v.region === '' || v.region.length > REGION_MAX) {
    throw invalid('genotyping.variant.region is not valid', { field: 'genotyping.variant.region' });
  }
  if (!Number.isSafeInteger(v.position) || v.position < 1) {
    throw invalid('genotyping.variant.position must be a positive integer', { field: 'genotyping.variant.position' });
  }
  ['ref', 'alt'].forEach(function (k) {
    if (typeof v[k] !== 'string' || !ALLELE_RE.test(v[k])) {
      throw invalid('genotyping.variant.' + k + ' must be 1-50 A, C, G or T bases', { field: 'genotyping.variant.' + k });
    }
  });
  if (v.ref.toUpperCase() === v.alt.toUpperCase()) {
    throw invalid('genotyping.variant.ref and alt must be different alleles', { field: 'genotyping.variant.alt' });
  }
  if (!Array.isArray(g.sets) || g.sets.length < 1 || g.sets.length > MAX_SETS) {
    throw invalid('genotyping.sets must be an array of 1-' + MAX_SETS + ' sets', { field: 'genotyping.sets', max: MAX_SETS });
  }
  const ids = new Set();
  g.sets.forEach(function (s, i) {
    const where = 'genotyping.sets[' + i + ']';
    if (!isPlainObject(s)) throw invalid(where + ' must be an object', { field: where });
    checkKeys(s, SET_KEYS, where + '.');
    SET_KEYS.forEach(function (k) {
      if (typeof s[k] !== 'string' || !SET_ID_RE.test(s[k])) throw invalid(where + '.' + k + ' is not valid', { field: where + '.' + k });
    });
    if (ids.has(s.id)) throw invalid('genotyping set id ' + s.id + ' is used more than once', { field: where + '.id' });
    ids.add(s.id);
  });
}

// ---- step 4: links between the sets and the pairs (§5.2 part A, pure) ----------------------------------------------------

function pairMap(pairs) {
  if (pairs instanceof Map) return pairs;
  const m = new Map();
  (Array.isArray(pairs) ? pairs : []).forEach(function (p) {
    if (p && typeof p.id === 'string' && !m.has(p.id)) m.set(p.id, p);
  });
  return m;
}

// pairs: the normalized pairs (an array, or a Map by id). For each set, in order: both pairs exist and differ, belong to no
// other set, carry identical expected locations, and share exactly one primer on the same side.
// -> [{id, ref_pair, alt_pair, orientation, primers {as_ref, as_alt, common} (uppercase), expected {region, start, end}}]
function validateLinks(genotyping, pairs) {
  const byId = pairMap(pairs);
  const usedBy = new Map();
  return genotyping.sets.map(function (s) {
    const ids = [s.ref_pair, s.alt_pair];
    const missing = ids.filter(function (id, i) { return !byId.has(id) && ids.indexOf(id) === i; });
    if (missing.length) throw setInvalid(s.id, 'unknown_pair', missing, 'no pair has the id ' + missing.join(' or '));
    if (s.ref_pair === s.alt_pair) {
      throw setInvalid(s.id, 'same_pair', [s.ref_pair], 'ref_pair and alt_pair must name two different pairs, not ' + s.ref_pair + ' twice');
    }
    const reused = ids.filter(function (id) { return usedBy.has(id); });
    if (reused.length) {
      throw setInvalid(s.id, 'pair_reused', reused, 'pair ' + reused[0] + ' already belongs to set ' + usedBy.get(reused[0]));
    }
    ids.forEach(function (id) { usedBy.set(id, s.id); });

    const ref = byId.get(s.ref_pair);
    const alt = byId.get(s.alt_pair);
    const lacking = [ref, alt].filter(function (p) { return !isPlainObject(p.expected); }).map(function (p) { return p.id; });
    if (lacking.length) {
      throw setInvalid(s.id, 'expected_required', lacking,
        lacking.join(' and ') + (lacking.length === 1 ? ' needs' : ' need') + ' an expected product location');
    }
    const a = ref.expected;
    const b = alt.expected;
    if (a.region !== b.region || a.start !== b.start || a.end !== b.end) {
      throw setInvalid(s.id, 'expected_differs', ids, s.ref_pair + ' and ' + s.alt_pair + ' must have the same expected product location');
    }
    const left = [String(ref.left).toUpperCase(), String(alt.left).toUpperCase()];
    const right = [String(ref.right).toUpperCase(), String(alt.right).toUpperCase()];
    if ((left[0] === left[1]) === (right[0] === right[1])) {
      throw setInvalid(s.id, 'no_shared_common', ids,
        s.ref_pair + ' and ' + s.alt_pair + ' must share exactly one primer, the common primer, on the same side');
    }
    const forward = right[0] === right[1];
    return {
      id: s.id,
      ref_pair: s.ref_pair,
      alt_pair: s.alt_pair,
      orientation: forward ? 'forward' : 'reverse',
      primers: forward ? { as_ref: left[0], as_alt: left[1], common: right[0] } : { as_ref: right[0], as_alt: right[1], common: left[0] },
      expected: { region: a.region, start: a.start, end: a.end }
    };
  });
}

// ---- step 6: the variant and the positions of each set (§5.2 part B) --------------------------------------------------------

function windowFor(variant, regionLength) {
  const end = variant.position + variant.ref.length + WINDOW_PAD;
  return {
    start: Math.max(1, variant.position - WINDOW_PAD),
    end: Number.isSafeInteger(regionLength) ? Math.min(regionLength, end) : end
  };
}

function tooRepetitive(region, position, shift, max) {
  return new PrimerHttpError(400, 'VARIANT_TOO_REPETITIVE',
    'the variant at ' + region + ':' + position + ' can slide too far to have a stable discriminating base',
    { region: region, position: position, shift: shift, max: max });
}

function outOfBounds(v, regionLength) {
  const end = v.position + v.ref.length - 1;
  return new PrimerHttpError(400, 'REGION_OUT_OF_BOUNDS',
    'the genotyping variant ' + v.region + ':' + v.position + '-' + end + ' reaches beyond region ' + v.region +
    (regionLength === null ? '' : ' (' + regionLength + ' bp)'),
    { region: v.region, start: v.position, end: end, length: regionLength });
}

// A read error while left-aligning or scanning the variant. OUTSIDE_REGION: the event reaches the region's end.
// SEQUENCE_WINDOW (it slides past the ±1,700 bp read, so far beyond max_shift) and NO_DISCRIMINATING_BASE: too repetitive.
function scanError(err, v, regionLength, shift, max) {
  const code = err && err.code;
  if (code === 'OUTSIDE_REGION') return outOfBounds(v, regionLength);
  if (code === 'SEQUENCE_WINDOW' || code === 'NO_DISCRIMINATING_BASE') return tooRepetitive(v.region, v.position, shift, max);
  return err;
}

// The reference span [start, end] read on both haplotypes; the strings are null when the span runs past the window read.
function haplotypeSpan(v, genome, start, end) {
  try {
    const h = variation.haplotypes(v, genome, start, end);
    return { start: start, end: end, ref: h.ref, alt: h.alt };
  } catch (e) {
    if (e && e.code === 'SEQUENCE_WINDOW') return { start: start, end: end, ref: null, alt: null };
    throw e;
  }
}

// realign's 3'-anchored alignment of a primer on a haplotype string that starts at hapStart, with its expected 3' end at
// anchor: over the whole haplotype, or only over anchor ± span.
function alignOn(primer, face, hap, hapStart, anchor, span) {
  let seq = hap;
  let start = hapStart;
  if (span !== undefined) {
    start = Math.max(hapStart, anchor - span);
    seq = hap.slice(start - hapStart, Math.min(hapStart + hap.length - 1, anchor + span) - hapStart + 1);
  }
  return realign.realignSiteOnWindow({ key: 'genotype', primer: primer, face: face, p3: anchor }, seq, start, THREE_PRIME_WINDOW);
}

function edits(a) {
  return 'edits at ' + a.mm_pos.join(', ') + ' nt from its 3\' end' + (a.gaps > 0 ? ', ' + a.gaps + ' of them indel bases' : '');
}

// §5.2 part B for one set, in order: each allele-specific primer ends at the variant on its own haplotype
// (not_at_variant); neither fits the other haplotype better (alleles_swapped); each carries at most one mismatch, at
// distance 2 or 3, and no indel (too_many_edits); each fits its own haplotype strictly better (alleles_swapped); the expected
// product is on the variant's region and the common primer, placed by expected, lies beyond the zone ± 1 (common_in_zone);
// the REF-specific primer's 5' end is expected.start (forward) or expected.end (reverse) (expected_mismatch).
function checkSet(link, c) {
  const forward = link.orientation === 'forward';
  const face = forward ? 'F' : 'R';
  const disc = c.d[link.orientation];
  // The allele-specific 3' base on each haplotype: forward primers keep their indices, reverse ones shift by delta (§4.9).
  const anchor = { ref: disc.position, alt: forward ? disc.position : disc.position + c.delta };
  const at = c.region + ':' + disc.position;
  const p = link.primers;
  const both = [link.ref_pair, link.alt_pair];
  const roles = Object.keys(ROLES);
  function primerError(role, reason, message, extra) {
    return setInvalid(link.id, reason, [link[ROLES[role].pair]], message,
      Object.assign({ primer: ROLES[role].pair, sequence: p[role] }, extra || {}));
  }

  const own = {};
  const other = {};
  roles.forEach(function (role) {
    const r = ROLES[role];
    own[role] = alignOn(p[role], face, c.hap[r.own], c.win.start, anchor[r.own]);
    if (!own[role] || own[role].p3 !== anchor[r.own]) {
      throw primerError(role, 'not_at_variant', 'the ' + r.label + ' ' + p[role] + ' does not end at the variant (' + at + ')');
    }
  });
  // The other haplotype is read at its own anchor only, so a copy of the primer elsewhere in the window cannot pass for it.
  roles.forEach(function (role) {
    const r = ROLES[role];
    other[role] = alignOn(p[role], face, c.hap[r.other], c.win.start, anchor[r.other], p[role].length + Math.abs(c.delta) + realign.PAD);
    if (other[role].mm < own[role].mm) {
      throw primerError(role, 'alleles_swapped', 'the ' + r.label + ' ' + p[role] + ' matches the ' + ALLELE_NAME[r.other] +
        ' allele better than the ' + ALLELE_NAME[r.own] + ' allele at ' + at + '; are ref_pair and alt_pair swapped?');
    }
  });
  roles.forEach(function (role) {
    const a = own[role];
    const deliberate = a.mm_pos.length === 1 && DELIBERATE_DISTANCES.indexOf(a.mm_pos[0]) >= 0;
    if (a.gaps > 0 || (a.mm_pos.length > 0 && !deliberate)) {
      throw primerError(role, 'too_many_edits', 'the ' + ROLES[role].label + ' ' + p[role] + ' may differ from the ' +
        ALLELE_NAME[ROLES[role].own] + ' allele at ' + at + ' only by one mismatch 2 or 3 nt from its 3\' end, but has ' + edits(a),
      { mm_pos: a.mm_pos.slice() });
    }
  });
  roles.forEach(function (role) {
    const r = ROLES[role];
    if (own[role].mm >= other[role].mm) {
      throw primerError(role, 'alleles_swapped', 'the ' + r.label + ' ' + p[role] + ' does not match the ' + ALLELE_NAME[r.own] +
        ' allele better than the ' + ALLELE_NAME[r.other] + ' allele at ' + at);
    }
  });

  const e = link.expected;
  if (e.region !== c.region) {
    throw setInvalid(link.id, 'expected_mismatch', both, 'the expected product is on region ' + e.region +
      ', not on the variant\'s region ' + c.region, { primer: 'ref_pair', sequence: p.as_ref });
  }
  const n = p.common.length;
  const common = forward ? { start: e.end - n + 1, end: e.end, strand: -1 } : { start: e.start, end: e.start + n - 1, strand: 1 };
  if (forward ? common.start <= c.zone.end + 1 : common.end >= c.zone.start - 1) {
    throw setInvalid(link.id, 'common_in_zone', both, 'the common primer ' + p.common + ' at ' + c.region + ':' + common.start + '-' +
      common.end + ' must lie ' + (forward ? 'after' : 'before') + ' the variant zone ' + c.region + ':' + c.zone.start + '-' + c.zone.end +
      ' and the base beside it', { primer: 'common', sequence: p.common });
  }
  const five = own.as_ref.p5;
  const want = forward ? e.start : e.end;
  if (five !== want) {
    throw setInvalid(link.id, 'expected_mismatch', both, 'the REF-specific primer ' + p.as_ref + ' has its 5\' end at ' + c.region + ':' +
      five + ', but expected.' + (forward ? 'start' : 'end') + ' is ' + want, { primer: 'ref_pair', sequence: p.as_ref });
  }

  const inTract = function (pos) { return c.tract !== null && pos !== null && pos >= c.tract.start && pos <= c.tract.end; };
  const mismatch = {
    as_ref: own.as_ref.mm_pos.length ? own.as_ref.mm_pos[0] : null,
    as_alt: own.as_alt.mm_pos.length ? own.as_alt.mm_pos[0] : null
  };
  const segment = { start: Math.max(1, e.start - c.pad), end: c.clampEnd(e.end + c.pad), sequence: null };
  try {
    segment.sequence = c.genome.bases(segment.start, segment.end);
  } catch (err) {
    if (!err || err.code !== 'SEQUENCE_WINDOW') throw err;
  }
  return {
    id: link.id,
    ref_pair: link.ref_pair,
    alt_pair: link.alt_pair,
    orientation: link.orientation,
    primers: { as_ref: p.as_ref, as_alt: p.as_alt, common: p.common },
    expected: { region: e.region, start: e.start, end: e.end },
    footprints: {
      as_ref: forward ? { start: five, end: anchor.ref, strand: 1 } : { start: anchor.ref, end: five, strand: -1 },
      common: common
    },
    deliberate_mismatch: mismatch,
    deliberate_mismatch_positions: [mismatch.as_ref, mismatch.as_alt].filter(function (x) { return x !== null; }),
    // §4.12: the 3' base of each allele-specific primer, on its own haplotype, inside the shift tract.
    in_shift_tract: {
      as_ref: inTract(disc.position),
      as_alt: c.tract !== null && (disc.alt_maps_to === null || inTract(disc.alt_maps_to))
    },
    segment: segment
  };
}

// Pure. genome: a variation/normalize.js sequenceWindow over the reference region of genotyping.variant (it must cover
// windowFor(), and carries start, end and regionLength). opts: {cfg (primers config), links (validateLinks' result)}.
function prepare(genotyping, pairs, genome, opts) {
  const o = opts || {};
  const cfg = o.cfg || {};
  const maxShift = count((cfg.variation || {}).max_shift, DEFAULTS.max_shift);
  const flankMin = count((cfg.check || {}).genotype_flank_min, DEFAULTS.genotype_flank_min);
  const pad = count((cfg.check || {}).genotype_amplicon_pad, DEFAULTS.genotype_amplicon_pad);
  const links = o.links || validateLinks(genotyping, pairs);
  const region = genotyping.variant.region;
  const regionLength = Number.isSafeInteger(genome.regionLength) ? genome.regionLength : null;
  const clampEnd = function (x) { return regionLength === null ? x : Math.min(regionLength, x); };
  const given = {
    region: region,
    position: genotyping.variant.position,
    ref: genotyping.variant.ref.toUpperCase(),
    alt: genotyping.variant.alt.toUpperCase()
  };

  // REF verification on the client's own representation, then the design's pure normalization (§3.5).
  let genomeRef;
  try {
    genomeRef = genome.bases(given.position, given.position + given.ref.length - 1);
  } catch (e) {
    throw e && e.code === 'OUTSIDE_REGION' ? outOfBounds(given, regionLength) : e;
  }
  if (genomeRef !== given.ref) {
    throw new PrimerHttpError(400, 'REF_MISMATCH', 'the reference allele ' + given.ref + ' does not match the genome ' +
      (genomeRef.length === 1 ? 'base ' : 'bases ') + genomeRef + ' at ' + region + ':' + given.position,
    { region: region, position: given.position, given: given.ref, genome: genomeRef });
  }
  let v = null;
  let shift = null;
  let d;
  try {
    v = variation.leftNormalize(given, genome);
    shift = variation.shiftOf(v, genome, { max_shift: maxShift });
    if (shift > maxShift) throw tooRepetitive(region, v.position, shift, maxShift);
    d = variation.discriminating(v, genome);
  } catch (e) {
    throw scanError(e, v ? { region: region, position: v.position, ref: v.ref } : given, regionLength, shift, maxShift);
  }
  const zone = variation.zoneOf(v, d);
  const win = { start: genome.start, end: genome.end };
  const ctx = {
    region: region,
    genome: genome,
    d: d,
    zone: zone,
    tract: variation.tractOf(v, shift),
    delta: v.alt.length - v.ref.length,
    hap: variation.haplotypes(v, genome, win.start, win.end),
    win: win,
    pad: pad,
    clampEnd: clampEnd
  };
  const sets = links.map(function (link) { return checkSet(link, ctx); });

  // §5.6: the core is the left-aligned VCF span with its anchor base, extended right by shift, ± 1 bp; K flanks it.
  const core = haplotypeSpan(v, genome, Math.max(1, v.position - 1), clampEnd(v.position + v.ref.length + shift));
  const flank = Math.max(flankMin, (zone.end - zone.start + 1) + Math.max(v.ref.length, v.alt.length));
  return {
    algorithm_version: GENOTYPING_VERSION,
    variant: {
      key: region + ':' + v.position + ':' + v.ref + ':' + v.alt,
      region: region,
      position: v.position,
      ref: v.ref,
      alt: v.alt,
      region_length: regionLength,
      shift: shift,
      tract: ctx.tract,
      zone: zone,
      discriminating: d,
      flank: flank,
      core: core,
      haplotypes: haplotypeSpan(v, genome, Math.max(1, core.start - flank), clampEnd(core.end + flank))
    },
    window: { region: region, start: win.start, end: win.end },
    sets: sets
  };
}

// reference: the resolved reference assembly (fasta.dna). deps: {cfg, sequence {regionLength, fetch} (default ../sequence)}.
// Reads windowFor(variant) once and returns prepare()'s object.
async function validateSets(genotyping, pairs, reference, deps) {
  deps = deps || {};
  const cfg = deps.cfg || require('../config').get();
  const sequence = deps.sequence || require('../sequence');
  const links = validateLinks(genotyping, pairs);
  const v = genotyping.variant;
  const systemName = reference && reference.system_name;
  const fasta = reference && reference.fasta && reference.fasta.dna;
  if (!fasta) {
    throw new PrimerHttpError(422, 'NO_SEQUENCE', 'no genome sequence is available for ' + systemName, { system_name: systemName });
  }
  const regionLength = await sequence.regionLength(fasta, v.region);
  if (!Number.isSafeInteger(regionLength) || regionLength < 1) {
    throw new PrimerHttpError(404, 'UNKNOWN_REGION', 'region ' + JSON.stringify(String(v.region)) + ' is not in this assembly',
      { region: String(v.region) });
  }
  if (v.position + v.ref.length - 1 > regionLength) throw outOfBounds(v, regionLength);
  const win = windowFor(v, regionLength);
  const seq = await sequence.fetch(fasta, v.region, win.start, win.end, 1, { maxLength: win.end - win.start + 1 });
  return prepare(genotyping, pairs, variation.sequenceWindow(seq, win.start, regionLength), { cfg: cfg, links: links });
}

// ---- step 8 ---------------------------------------------------------------------------------------------------------------------

// request.genotyping: the client's keys only, with the variant uppercased and left-aligned (§5.3).
function requestBlock(prepared) {
  const v = prepared.variant;
  return {
    variant: { region: v.region, position: v.position, ref: v.ref, alt: v.alt },
    sets: prepared.sets.map(function (s) { return { id: s.id, ref_pair: s.ref_pair, alt_pair: s.alt_pair }; })
  };
}

module.exports = {
  GENOTYPING_VERSION,
  MODES,
  MAX_SETS,
  WINDOW_PAD,
  REASONS,
  validateShape,
  validateLinks,
  windowFor,
  validateSets,
  prepare,
  requestBlock
};
