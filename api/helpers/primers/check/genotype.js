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
// Worker side (check/run.js, §5.6-§5.9), below the submit side: fillSequences, megablastQuery, emptyResults, callGenome,
// semiGlobal, mergeCopies, gapCompressedIdentity, megablastCopies, predict, writeResults and unavailableEntry. They start from
// prepare()'s object, which normalize stores as resolved.genotyping (the job doc only, never returned to clients); the worker
// re-derives it with validateSets() over deps.sequence = {fetch: fetchWindow, regionLength} and fills in, with fillSequences(),
// any sequence that ran past the submit-time window. BLAST is never loaded here: run.js runs the megablast fallback.
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
const classify = require('./classify');

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

// ==== worker side: the per-assembly allele caller (§5.6), the amplification prediction (§5.7), results (§2.12, §5.8) =========
//
// run.js prepares once per job (validateSets, fillSequences), calls callGenome for the reference (the control) and for every
// pan-genome genome with a genome target, and rebuilds results.genotyping with writeResults before each flush.
//   products (per set, in prepared.sets order): {ref: {amplicons, unlikely}, alt: {amplicons, unlikely}}, the REF and ALT pairs'
//     amplicons.finalize lists of that task; left_* describes the pair's left primer, so the allele-specific site is left_* for
//     a forward set and right_* for a reverse one
//   entry = {genome: PrimerGenotypeGenome, sets: [PrimerGenotypeSetGenome per set]}
// Alignments share one shape: {ops (M X I D per column; I = a query base against a gap, D = a genome base against a gap),
// qToT (genome index of each query base, -1 on a gap), qCol (column of each query base), qFirst, qLast (the aligned query
// span), target (the genome bases, on the reference strand)}.

// Config fallbacks (primers.check.genotype_*, §3.1).
const CALLER_DEFAULTS = Object.freeze({
  genotype_amplicon_pad: 50,
  genotype_ortholog_min_identity: 95,
  genotype_ortholog_size_tolerance: 0.2,
  genotype_max_copies: 10,
  genotype_max_anchors: 50,
  genotype_megablast_min_identity: 95,
  genotype_megablast_min_query_cover: 0.8,
  genotype_megablast_min_bitscore_frac: 0.9,
  genotype_offlocus_max_mismatches: 2
});
// WEAK_OFF_TARGETS details list at most this many example products.
const WEAK_OFF_TARGET_EXAMPLES = 5;
// The megablast query is ref[zone.start - 200 - K, zone.end + 200 + K] (§5.6 step 8).
const MEGABLAST_FLANK = 200;
// An anchor's genome window is the product ± WINDOW_PADS x genotype_amplicon_pad, not ± 1 pad as §5.6 step 2 words it: the
// genome ends are free, so the extra bases add no edit, while an indel inside the pad would otherwise cut R_s short (the plus
// copy of pi180348 has a 1 bp insertion at 1:11050, inside the left pad of set S2 of §2.13).
const WINDOW_PADS = 2;
// paralog_copies counts the copies that fail the ortholog test with a gap-compressed identity of at least this. Below it an
// anchor is an unrelated off-target product, not a copy of the locus: a real check yields 11-16 of those per assembly at
// rs871475760, all at 66-77 %, while pi180348's 1:72.9 Mb paralog of the locus aligns at 90 %.
const PARALOG_MIN_IDENTITY = 80;
// Anchor windows on one region at most this far apart are read with one fetch of at most FETCH_MAX_SPAN bases.
const FETCH_MERGE_GAP = 1000;
const FETCH_MAX_SPAN = 20000;

const ALLELES = Object.freeze(['ref', 'alt', 'other', 'ambiguous', 'missing', 'unavailable']);
const COPY_CALLS = Object.freeze(['ref', 'alt', 'other', 'missing']);
// Best first: the order in which a primer's status wins over its products and copies (§5.7).
const STATUSES = Object.freeze(['match', 'weak', 'terminal_mismatch', 'uncertain', 'unknown', 'blocked', 'no_product']);
const PREDICTIONS = Object.freeze(['ref', 'alt', 'both', 'none', 'no_call', 'unknown']);
const GENOME_REASONS = Object.freeze(['variant_not_covered', 'no_orthologous_copy', 'fallback_failed', 'fallback_budget',
  'db_unavailable', 'blast_error', 'call_failed']);
// Also the order in which a set call lists its reasons.
const SET_REASONS = Object.freeze(['common_primer_3p_mismatch', 'common_primer_weak', 'ref_signal_off_locus', 'alt_signal_off_locus',
  'shift_tract_uncertain', 'approx_alignment', 'no_orthologous_copy', 'third_allele']);
const CONTROL_STATUSES = Object.freeze(['pass', 'warn', 'fail']);
const CONTROL_REASONS = Object.freeze(['allele_not_ref', 'prediction_not_ref', 'weak', 'off_locus_products']);

const DIAG = 1;
const UP = 2;
const LEFT = 3;

function cfgNumber(ccfg, key) {
  const v = ccfg ? ccfg[key] : undefined;
  return typeof v === 'number' && Number.isFinite(v) && v >= 0 ? v : CALLER_DEFAULTS[key];
}

function cfgInteger(ccfg, key) {
  const v = ccfg ? ccfg[key] : undefined;
  return Number.isInteger(v) && v >= 0 ? v : CALLER_DEFAULTS[key];
}

// 100 x part / whole, rounded half away from zero to 2 decimals in integer arithmetic (§2.1); null without a whole.
function percent(part, whole) {
  if (!(whole > 0)) return null;
  return Math.floor((20000 * part + whole) / (2 * whole)) / 100;
}

function isPaired(product) {
  return product.orientation === 'LR' || product.orientation === 'RL';
}

function listOf(side) {
  return side ? (side.amplicons || []).concat(side.unlikely || []) : [];
}

// ---- sequences for the worker ------------------------------------------------------------------------------------------------

// Fills in, with fetch(region, start, end) -> plus-strand bases, every sequence of prepare()'s object that ran past the
// submit-time window (null): set segments, the core and the haplotypes. Mutates and returns prepared.
async function fillSequences(prepared, fetch) {
  const v = prepared.variant;
  for (const s of prepared.sets) {
    if (s.segment.sequence === null) s.segment.sequence = String(await fetch(v.region, s.segment.start, s.segment.end)).toUpperCase();
  }
  const spans = [v.core, v.haplotypes].filter(function (sp) { return sp.ref === null || sp.alt === null; });
  if (spans.length) {
    const start = Math.min.apply(null, spans.map(function (sp) { return sp.start; }));
    const end = Math.max.apply(null, spans.map(function (sp) { return sp.end; }));
    const genome = variation.sequenceWindow(String(await fetch(v.region, start, end)), start, v.region_length);
    spans.forEach(function (sp) {
      const h = variation.haplotypes({ position: v.position, ref: v.ref, alt: v.alt }, genome, sp.start, sp.end);
      sp.ref = h.ref;
      sp.alt = h.alt;
    });
  }
  return prepared;
}

// The megablast fallback query, clamped to the region: {region, start, end}.
function megablastQuery(prepared) {
  const v = prepared.variant;
  const end = v.zone.end + MEGABLAST_FLANK + v.flank;
  return {
    region: v.region,
    start: Math.max(1, v.zone.start - MEGABLAST_FLANK - v.flank),
    end: Number.isSafeInteger(v.region_length) ? Math.min(v.region_length, end) : end
  };
}

// ---- alignment and the core (§5.6 steps 2-4) ---------------------------------------------------------------------------------

// Semi-global unit-cost alignment of query against target: the query is fully consumed and the target's ends are free. The
// end is the leftmost target index of least cost; the traceback prefers the diagonal, then a query base against a gap, then a
// target base against a gap, except that it stays in a gap run while that is still optimal, so that one indel is one run for
// gap_compressed_identity (a 35 bp deletion would otherwise be split into many runs). -> an alignment (see above) plus cost.
function semiGlobal(query, target) {
  const q = String(query).toUpperCase();
  const t = String(target).toUpperCase();
  const n = q.length;
  const m = t.length;
  if (n === 0) throw new TypeError('semiGlobal: the query is empty');
  const cols = m + 1;
  const D = new Int32Array((n + 1) * cols);
  const dir = new Uint8Array((n + 1) * cols);
  for (let i = 1; i <= n; i++) {
    const row = i * cols;
    const prev = row - cols;
    const qc = q.charCodeAt(i - 1);
    D[row] = i;
    dir[row] = UP;
    for (let j = 1; j <= m; j++) {
      const diag = D[prev + j - 1] + (qc === t.charCodeAt(j - 1) ? 0 : 1);
      const up = D[prev + j] + 1;
      const left = D[row + j - 1] + 1;
      let best = diag;
      let d = DIAG;
      if (up < best) { best = up; d = UP; }
      if (left < best) { best = left; d = LEFT; }
      D[row + j] = best;
      dir[row + j] = d;
    }
  }
  const last = n * cols;
  let end = 0;
  for (let j = 1; j <= m; j++) if (D[last + j] < D[last + end]) end = j;
  const ops = [];
  const qToT = new Int32Array(n).fill(-1);
  let i = n;
  let j = end;
  let previous = 0;
  while (i > 0) {
    let d = j === 0 ? UP : dir[i * cols + j];
    if (previous === UP && d !== UP && D[i * cols + j] === D[(i - 1) * cols + j] + 1) d = UP;
    else if (previous === LEFT && d !== LEFT && D[i * cols + j] === D[i * cols + j - 1] + 1) d = LEFT;
    previous = d;
    if (d === DIAG) {
      ops.push(q.charCodeAt(i - 1) === t.charCodeAt(j - 1) ? 'M' : 'X');
      qToT[i - 1] = j - 1;
      i--;
      j--;
    } else if (d === UP) {
      ops.push('I');
      i--;
    } else {
      ops.push('D');
      j--;
    }
  }
  ops.reverse();
  const qCol = new Int32Array(n).fill(-1);
  for (let k = 0, qi = 0; k < ops.length; k++) if (ops[k] !== 'D') qCol[qi++] = k;
  return { cost: D[last + end], ops: ops.join(''), qToT: qToT, qCol: qCol, qFirst: 0, qLast: n - 1, target: t };
}

// Columns, edits, gap columns and gap runs (a run of I or of D columns) of an ops string.
function opsStats(ops) {
  let edits = 0;
  let gapColumns = 0;
  let gapRuns = 0;
  let prev = '';
  for (let k = 0; k < ops.length; k++) {
    const o = ops[k];
    if (o !== 'M') edits++;
    if ((o === 'I' || o === 'D') && o !== prev) gapRuns++;
    if (o === 'I' || o === 'D') gapColumns++;
    prev = o;
  }
  return { columns: ops.length, edits: edits, matches: ops.length - edits, gapColumns: gapColumns, gapRuns: gapRuns };
}

// §2.12: percent identity with each gap run counted as one column and one edit, 2 decimals. alignment: an ops string or {ops}.
function gapCompressedIdentity(alignment) {
  const st = opsStats(typeof alignment === 'string' ? alignment : alignment.ops);
  return percent(st.matches, st.columns - st.gapColumns + st.gapRuns);
}

// §5.6 step 3: the query indices qs..qe of the core mapped through an alignment. An end on a gap moves outward (the start
// left, the end right) to the nearest query base aligned to a genome base; an end outside the aligned query span, or no such
// base, gives null. -> {lo, hi (genome indices), c0, c1 (columns)}
function readCore(al, qs, qe) {
  if (qs < al.qFirst || qe > al.qLast || qs > qe) return null;
  let a = qs;
  while (a >= al.qFirst && al.qToT[a] < 0) a--;
  let b = qe;
  while (b <= al.qLast && al.qToT[b] < 0) b++;
  if (a < al.qFirst || b > al.qLast || al.qToT[b] < al.qToT[a]) return null;
  return { lo: al.qToT[a], hi: al.qToT[b], c0: al.qCol[a], c1: al.qCol[b] };
}

// One alignment of the reference span starting at segStart (an anchor's R_s, or the megablast query): its figures, the
// exact-core call (§5.6 step 4), the edits outside the core and the genome coordinate aligned to vcf.position.
// coord(genome index) -> genome coordinate.
function readAlignment(al, segStart, variant, coord) {
  const st = opsStats(al.ops);
  const core = readCore(al, variant.core.start - segStart, variant.core.end - segStart);
  const qv = variant.position - segStart;
  const vi = qv >= al.qFirst && qv <= al.qLast ? al.qToT[qv] : -1;
  const out = {
    columns: st.columns,
    matches: st.matches,
    gc_columns: st.columns - st.gapColumns + st.gapRuns,
    identity: percent(st.matches, st.columns),
    gap_compressed_identity: percent(st.matches, st.columns - st.gapColumns + st.gapRuns),
    aligned_length: st.columns,
    observed: null,
    flank_edits: null,
    call: 'missing',
    variant_position: vi >= 0 ? coord(vi) : null,
    locus: null
  };
  if (core) {
    let inCore = 0;
    for (let k = core.c0; k <= core.c1; k++) if (al.ops[k] !== 'M') inCore++;
    out.observed = al.target.slice(core.lo, core.hi + 1);
    out.call = out.observed === variant.core.ref ? 'ref' : out.observed === variant.core.alt ? 'alt' : 'other';
    out.flank_edits = st.edits - inCore;
  }
  out.locus = out.variant_position !== null ? out.variant_position : core ? coord(core.lo) : null;
  return out;
}

// ---- anchors, copies and orthologs (§5.6 steps 1, 5, 6) --------------------------------------------------------------------

// Every LR/RL product of each set's two pairs, likely ones first, de-duplicated by (region, start, end) within the set, at most
// maxAnchors per set. Approximate products are left out: a site that was not re-aligned (over the mismatch cap, or without
// FASTA) has extrapolated ends, so its product can span a secondary site of one primer beside the locus (S1 1:11068-11173 in the
// reference, next to the real 1:11068-11132). -> [{set, product}]
function anchorsOf(prepared, products, maxAnchors) {
  const out = [];
  prepared.sets.forEach(function (set, si) {
    const p = products[si] || {};
    const seen = new Set();
    [p.ref && p.ref.amplicons, p.alt && p.alt.amplicons, p.ref && p.ref.unlikely, p.alt && p.alt.unlikely].forEach(function (list) {
      (list || []).forEach(function (a) {
        const key = a.region + '\t' + a.start + '\t' + a.end;
        if (seen.size >= maxAnchors || !isPaired(a) || a.approx === true || seen.has(key)) return;
        seen.add(key);
        out.push({ set: si, product: a });
      });
    });
  });
  return out;
}

// Reads [start - pad, end + pad] of every anchor (pad: WINDOW_PADS x genotype_amplicon_pad), clamped to its region, sharing
// fetches between nearby windows. Sets
// anchor.lo and anchor.window (uppercase plus strand), or anchor.window = null when it cannot be read. -> failed reads
async function readWindows(anchors, pad, deps) {
  const lengths = new Map();
  for (const a of anchors) {
    const region = a.product.region;
    if (!lengths.has(region)) {
      let len;
      try {
        len = typeof deps.regionLength === 'function' ? await deps.regionLength(region) : undefined;
      } catch (e) {
        len = undefined;
      }
      lengths.set(region, len);
    }
    const len = lengths.get(region);
    a.lo = Math.max(1, a.product.start - pad);
    a.hi = Number.isSafeInteger(len) ? Math.min(len, a.product.end + pad) : a.product.end + pad;
    a.window = null;
  }
  const sorted = anchors.filter(function (a) { return a.hi >= a.lo; }).sort(function (x, y) {
    return x.product.region < y.product.region ? -1 : x.product.region > y.product.region ? 1 : x.lo - y.lo;
  });
  const chunks = [];
  sorted.forEach(function (a) {
    const c = chunks[chunks.length - 1];
    if (c && c.region === a.product.region && a.lo <= c.hi + FETCH_MERGE_GAP && Math.max(c.hi, a.hi) - c.lo + 1 <= FETCH_MAX_SPAN) {
      c.hi = Math.max(c.hi, a.hi);
      c.anchors.push(a);
    } else {
      chunks.push({ region: a.product.region, lo: a.lo, hi: a.hi, anchors: [a] });
    }
  });
  async function read(region, lo, hi) {
    try {
      const seq = await deps.fetch(region, lo, hi);
      return typeof seq === 'string' && seq.length === hi - lo + 1 ? seq.toUpperCase() : null;
    } catch (e) {
      if (e && (e.code === 'ABORTED' || e.name === 'AbortError')) throw e;
      return null;
    }
  }
  let failed = 0;
  for (const c of chunks) {
    const seq = c.anchors.length > 1 ? await read(c.region, c.lo, c.hi) : null;
    for (const a of c.anchors) {
      a.window = seq !== null ? seq.slice(a.lo - c.lo, a.hi - c.lo + 1) : await read(c.region, a.lo, a.hi);
      if (a.window === null) failed++;
    }
  }
  return failed;
}

// §5.6 step 5. records (in product order): readAlignment() results plus {region, start, end, strand, anchor (a product, not an
// HSP), ortholog, size, ref_size}. Records on one region and strand whose loci (the aligned vcf.position, else the first core
// base) differ by at most opts.zoneLength are one copy: start/end the envelope, anchors the number of product records, and the
// figures of the member with the most columns (the first on ties). -> [{copy (§2.12 PrimerGenotypeCopy), best, members, region,
// strand, locus}]
function sameLocus(group, record, tolerance) {
  return group.locus !== null && record.locus !== null && group.region === record.region && group.strand === record.strand &&
    Math.abs(group.locus - record.locus) <= tolerance;
}

function mergeCopies(records, opts) {
  const o = opts || {};
  const tolerance = Number.isFinite(o.zoneLength) ? o.zoneLength : 1;
  const groups = [];
  records.forEach(function (r) {
    let g = null;
    if (r.locus !== null) {
      g = groups.find(function (x) { return sameLocus(x, r, tolerance); }) || null;
    }
    if (!g) {
      g = { region: r.region, strand: r.strand, locus: r.locus, members: [] };
      groups.push(g);
    }
    g.members.push(r);
  });
  return groups.map(function (g) {
    let best = g.members[0];
    g.members.forEach(function (m) { if (m.columns > best.columns) best = m; });
    const anchors = g.members.filter(function (m) { return m.anchor; });
    const flags = anchors.map(function (m) { return m.ortholog; });
    const variantPosition = best.variant_position !== null ? best.variant_position
      : (g.members.find(function (m) { return m.variant_position !== null; }) || { variant_position: null }).variant_position;
    return {
      copy: {
        region: g.region,
        start: Math.min.apply(null, g.members.map(function (m) { return m.start; })),
        end: Math.max.apply(null, g.members.map(function (m) { return m.end; })),
        strand: g.strand,
        variant_position: variantPosition,
        identity: best.identity,
        gap_compressed_identity: best.gap_compressed_identity,
        aligned_length: best.aligned_length,
        observed: best.observed,
        flank_edits: best.flank_edits,
        call: best.call,
        anchors: anchors.length,
        ortholog: flags.some(function (f) { return f === true; }) ? true
          : flags.length && flags.every(function (f) { return f === false; }) ? false : null,
        source: anchors.length ? 'amplicon' : 'megablast'
      },
      best: best,
      members: g.members,
      region: g.region,
      strand: g.strand,
      locus: g.locus
    };
  });
}

function sizeOk(record, tolerance) {
  return !record.anchor || !(record.ref_size > 0) || Math.abs(record.size - record.ref_size) <= tolerance * record.ref_size;
}

// §5.6 step 6: annotated as an ortholog, a megablast copy (its HSP passed the fallback filters), or gap-compressed identity of
// the longest member >= genotype_ortholog_min_identity with every anchor's size within the tolerance of its set's product.
function isOrthologous(group, ccfg) {
  if (group.copy.ortholog === true || group.copy.source === 'megablast') return true;
  const minIdentity = cfgNumber(ccfg, 'genotype_ortholog_min_identity');
  const tolerance = cfgNumber(ccfg, 'genotype_ortholog_size_tolerance');
  if (!(100 * group.best.matches >= minIdentity * group.best.gc_columns)) return false;
  return group.members.every(function (m) { return sizeOk(m, tolerance); });
}

// The alignment of one megablast row over a query of qlen bases, or null when its columns disagree with its coordinates.
function hspAlignment(row, qlen) {
  const q = row.qseq.toUpperCase();
  const s = row.sseq.toUpperCase();
  const qToT = new Int32Array(qlen).fill(-1);
  const qCol = new Int32Array(qlen).fill(-1);
  let ops = '';
  let target = '';
  let qi = row.qstart - 1;
  for (let k = 0; k < q.length; k++) {
    if (q[k] === '-') {
      if (s[k] === '-') return null;
      ops += 'D';
      target += s[k];
      continue;
    }
    if (qi >= qlen) return null;
    qCol[qi] = k;
    if (s[k] === '-') {
      ops += 'I';
    } else {
      ops += q[k] === s[k] ? 'M' : 'X';
      qToT[qi] = target.length;
      target += s[k];
    }
    qi++;
  }
  if (qi !== row.qend || target.length !== Math.abs(row.send - row.sstart) + 1) return null;
  return { ops: ops, qToT: qToT, qCol: qCol, qFirst: row.qstart - 1, qLast: row.qend - 1, target: target };
}

// §5.6 step 8 over parsed megablast rows (check/blast.js parseMegablastLine) of the query ref[query.start, query.end]: HSPs with
// bitscore >= frac x the best, pident >= the minimum and query cover >= the minimum whose query columns cover the core are read
// like anchors (sseq is already on the query strand) and merged. -> {copies: mergeCopies groups, hsps, kept, covering}
function megablastCopies(rows, query, prepared, ccfg) {
  const v = prepared.variant;
  const qlen = query.end - query.start + 1;
  const valid = (rows || []).filter(function (r) {
    return r && typeof r.qseq === 'string' && typeof r.sseq === 'string' && Number.isFinite(r.bitscore) && Number.isFinite(r.pident) &&
      r.qstart >= 1 && r.qend <= qlen;
  });
  const best = valid.reduce(function (b, r) { return Math.max(b, r.bitscore); }, 0);
  const kept = valid.filter(function (r) {
    return r.bitscore >= cfgNumber(ccfg, 'genotype_megablast_min_bitscore_frac') * best &&
      r.pident >= cfgNumber(ccfg, 'genotype_megablast_min_identity') &&
      r.qend - r.qstart + 1 >= cfgNumber(ccfg, 'genotype_megablast_min_query_cover') * qlen;
  });
  const qs = v.core.start - query.start + 1;
  const qe = v.core.end - query.start + 1;
  const covering = kept.filter(function (r) { return r.qstart <= qs && r.qend >= qe; });
  const records = [];
  covering.forEach(function (r) {
    const al = hspAlignment(r, qlen);
    if (!al) return;
    const plus = r.strand === 1;
    const read = readAlignment(al, query.start, v, function (k) { return plus ? r.sstart + k : r.sstart - k; });
    records.push(Object.assign(read, {
      region: r.sseqid, start: Math.min(r.sstart, r.send), end: Math.max(r.sstart, r.send), strand: plus ? 1 : -1,
      anchor: false, ortholog: null, size: null, ref_size: null
    }));
  });
  return { copies: mergeCopies(records, { zoneLength: v.zone.end - v.zone.start + 1 }), hsps: valid.length, kept: kept.length, covering: covering.length };
}

// ---- amplification prediction (§5.7) ---------------------------------------------------------------------------------------

function primerCall(status, product, mmPos, residual) {
  return {
    status: status,
    likelihood: product ? product.likelihood : null,
    mm_pos: Array.isArray(mmPos) ? mmPos.slice() : null,
    residual_mm_pos: Array.isArray(residual) ? residual.slice() : null
  };
}

const NO_PRODUCT = Object.freeze({ status: 'no_product', likelihood: null, mm_pos: null, residual_mm_pos: null });

// mm_pos without the primer's own declared deliberate-mismatch position (once); null stays null.
function residualOf(mmPos, declared) {
  if (!Array.isArray(mmPos)) return null;
  const out = mmPos.slice();
  const i = Number.isInteger(declared) ? out.indexOf(declared) : -1;
  if (i >= 0) out.splice(i, 1);
  return out;
}

function siteOf(product, side) {
  return { mm: product[side + '_mm'], mm_pos: product[side + '_mm_pos'], mm_3p: product[side + '_3p_mm'] };
}

// The common primer's site alone: never a product likelihood, never a declared mismatch.
function commonStatus(site, params) {
  if (classify.isIgnored(site, params) || classify.isOverAmplifyingCap(site, params) || classify.isBlocked(site, params)) return 'blocked';
  if (!Array.isArray(site.mm_pos)) return 'unknown';
  if (site.mm_pos.indexOf(1) >= 0) return 'terminal_mismatch';
  if (site.mm_pos.indexOf(2) >= 0 || site.mm_pos.indexOf(3) >= 0) return 'weak';
  return 'match';
}

// An allele-specific primer in its own pair's product.
function allelePrimerCall(product, side, declared, inShiftTract) {
  const mm = product[side + '_mm_pos'];
  const residual = residualOf(mm, declared);
  let status;
  if (product.likelihood === 'unlikely') status = 'blocked';
  else if (!Array.isArray(mm)) status = 'unknown';
  else if (residual.indexOf(1) >= 0) status = inShiftTract ? 'uncertain' : 'terminal_mismatch';
  else if (residual.indexOf(2) >= 0 || residual.indexOf(3) >= 0) status = 'weak';
  else status = 'match';
  return primerCall(status, product, mm, residual);
}

// The best call: status order, then likelihood, then the first in product order; no_product without candidates.
function bestCall(calls) {
  let best = null;
  calls.forEach(function (c) {
    if (best === null) {
      best = c;
      return;
    }
    const d = STATUSES.indexOf(c.status) - STATUSES.indexOf(best.status);
    if (d < 0 || (d === 0 && classify.LIKELIHOOD_RANK[c.likelihood] < classify.LIKELIHOOD_RANK[best.likelihood])) best = c;
  });
  return best || Object.assign({}, NO_PRODUCT);
}

function unknownCall() {
  return { status: 'unknown', likelihood: null, mm_pos: null, residual_mm_pos: null };
}

function unknownPrediction() {
  return { ref_primer: unknownCall(), alt_primer: unknownCall(), common_primer: unknownCall(), predicted: 'unknown', strength: null,
    agrees: null, reasons: [], off_locus_products: 0 };
}

const amplifies = function (s) { return s === 'match' || s === 'weak'; };

// The mismatches of a product's allele-specific site, without the primer's own declared deliberate mismatch.
function allelePrimerMismatches(product, side, declared) {
  const mm = product[side + '_mm'];
  const pos = product[side + '_mm_pos'];
  return Array.isArray(pos) && Number.isInteger(declared) && pos.indexOf(declared) >= 0 ? mm - 1 : mm;
}

// §5.7 for one set on one genome. products: {ref: {amplicons, unlikely}, alt: {...}}; genome: {allele, copies (the orthologous
// PrimerGenotypeCopy objects)}; params: the check's classify params. opts: {max_offlocus_mismatches (default 2), weak (an array)}:
// an off-locus product changes the prediction only when each of its primers has at most max_offlocus_mismatches mismatches
// (M8b); one the rule would otherwise count is pushed to opts.weak as {which: 'ref' | 'alt', product} instead.
// -> {ref_primer, alt_primer, common_primer, predicted, strength, agrees, reasons, off_locus_products}
function predict(set, products, genome, params, opts) {
  const o = opts || {};
  const maxMismatches = Number.isInteger(o.max_offlocus_mismatches) && o.max_offlocus_mismatches >= 0
    ? o.max_offlocus_mismatches : CALLER_DEFAULTS.genotype_offlocus_max_mismatches;
  if (genome.allele === 'unavailable') return unknownPrediction();
  const p = products || {};
  const asSide = set.orientation === 'forward' ? 'left' : 'right';
  const commonSide = asSide === 'left' ? 'right' : 'left';
  const declared = set.deliberate_mismatch || {};
  const inTract = set.in_shift_tract || {};
  const copies = genome.copies || [];
  const onLocus = function (a) {
    return copies.some(function (c) { return String(c.region) === String(a.region) && a.start <= c.end && a.end >= c.start; });
  };
  const split = function (side) {
    const out = { on: [], off: [] };
    listOf(side).forEach(function (a) { if (isPaired(a)) (onLocus(a) ? out.on : out.off).push(a); });
    return out;
  };
  const ref = split(p.ref);
  const alt = split(p.alt);

  const refCall = bestCall(ref.on.map(function (a) { return allelePrimerCall(a, asSide, declared.as_ref, inTract.as_ref === true); }));
  const altCall = bestCall(alt.on.map(function (a) { return allelePrimerCall(a, asSide, declared.as_alt, inTract.as_alt === true); }));
  const commonCall = bestCall(ref.on.concat(alt.on).map(function (a) {
    const site = siteOf(a, commonSide);
    return primerCall(commonStatus(site, params), a, site.mm_pos, site.mm_pos);
  }));
  const statuses = [refCall.status, altCall.status, commonCall.status];

  const reasons = [];
  let predicted;
  let strength = null;
  if (!ref.on.length && !alt.on.length) {
    predicted = 'none';
    if (genome.allele === 'missing') reasons.push('no_orthologous_copy');
  } else if (commonCall.status === 'terminal_mismatch' || commonCall.status === 'blocked') {
    predicted = 'no_call';
    reasons.push('common_primer_3p_mismatch');
  } else if (statuses.indexOf('unknown') >= 0 || statuses.indexOf('uncertain') >= 0) {
    predicted = 'unknown';
    if (statuses.indexOf('unknown') >= 0) reasons.push('approx_alignment');
    if (statuses.indexOf('uncertain') >= 0) reasons.push('shift_tract_uncertain');
  } else {
    const r = amplifies(refCall.status);
    const a = amplifies(altCall.status);
    predicted = r && a ? 'both' : r ? 'ref' : a ? 'alt' : 'none';
    if (predicted !== 'none') {
      const weak = (r && refCall.status === 'weak') || (a && altCall.status === 'weak') || commonCall.status === 'weak';
      strength = weak ? 'weak' : 'normal';
      if (commonCall.status === 'weak') reasons.push('common_primer_weak');
      if (predicted === 'both' && genome.allele === 'other') reasons.push('third_allele');
    }
  }

  // Off-locus products add only the dye that is missing; unknown never changes. A product with more than maxMismatches
  // mismatches in either primer (up to max_amplifying_mismatches, so the check still lists it) is only reported.
  let offLocus = 0;
  const signal = { ref: false, alt: false };
  [['ref', ref.off, declared.as_ref], ['alt', alt.off, declared.as_alt]].forEach(function (row) {
    row[1].forEach(function (a) {
      if (!classify.countsAsAmplicon(a.likelihood)) return;
      const residual = residualOf(a[asSide + '_mm_pos'], row[2]);
      if (residual === null || residual.indexOf(1) >= 0 || !amplifies(commonStatus(siteOf(a, commonSide), params))) return;
      if (allelePrimerMismatches(a, asSide, row[2]) > maxMismatches || a[commonSide + '_mm'] > maxMismatches) {
        if (Array.isArray(o.weak)) o.weak.push({ which: row[0], product: a });
        return;
      }
      offLocus++;
      signal[row[0]] = true;
    });
  });
  if (signal.ref) reasons.push('ref_signal_off_locus');
  if (signal.alt) reasons.push('alt_signal_off_locus');
  if (offLocus && predicted !== 'unknown') {
    const has = predicted === 'none' || predicted === 'no_call'
      ? { ref: signal.ref, alt: signal.alt }
      : { ref: predicted === 'ref' || predicted === 'both' || signal.ref, alt: predicted === 'alt' || predicted === 'both' || signal.alt };
    predicted = has.ref && has.alt ? 'both' : has.ref ? 'ref' : has.alt ? 'alt' : predicted;
  }

  let agrees = null;
  if (['ambiguous', 'unavailable'].indexOf(genome.allele) < 0 && predicted !== 'unknown' && predicted !== 'no_call' &&
      statuses.indexOf('uncertain') < 0) {
    agrees = (genome.allele === 'ref' && predicted === 'ref') || (genome.allele === 'alt' && predicted === 'alt') ||
      ((genome.allele === 'other' || genome.allele === 'missing') && predicted === 'none');
  }
  reasons.sort(function (x, y) { return SET_REASONS.indexOf(x) - SET_REASONS.indexOf(y); });
  return { ref_primer: refCall, alt_primer: altCall, common_primer: commonCall, predicted: predicted, strength: strength,
    agrees: agrees, reasons: reasons, off_locus_products: offLocus };
}

// ---- one genome (§5.6 steps 1-9) --------------------------------------------------------------------------------------------

// input: {prepared, system_name, display_name, is_reference, products (per set), cfg (primers.check), params (classify params)}
// deps: {fetch(region, start, end) -> plus-strand bases, regionLength(region) -> length | undefined,
//        megablast() -> {status: 'ok', rows, query {start, end}} | {status: 'budget'} | {status: 'failed'}   (optional)}
// -> {genome, sets, megablast: null | 'ok' | 'budget' | 'failed', anchors, failed_reads}
async function callGenome(input, deps) {
  const prepared = input.prepared;
  const v = prepared.variant;
  const ccfg = input.cfg || {};
  const products = input.products || [];
  const d = deps || {};
  const zoneLength = v.zone.end - v.zone.start + 1;

  const anchors = anchorsOf(prepared, products, cfgNumber(ccfg, 'genotype_max_anchors'));
  const failedReads = await readWindows(anchors, WINDOW_PADS * cfgNumber(ccfg, 'genotype_amplicon_pad'), d);
  const records = [];
  anchors.forEach(function (a) {
    const set = prepared.sets[a.set];
    if (a.window === null || !set.segment.sequence) return;
    const minus = a.product.orientation === 'RL';
    const hi = a.lo + a.window.length - 1;
    const al = semiGlobal(set.segment.sequence, minus ? realign.revcomp(a.window) : a.window);
    const read = readAlignment(al, set.segment.start, v, minus ? function (k) { return hi - k; } : function (k) { return a.lo + k; });
    records.push(Object.assign(read, {
      region: a.product.region, start: a.product.start, end: a.product.end, strand: minus ? -1 : 1, anchor: true,
      ortholog: a.product.ortholog === undefined ? null : a.product.ortholog, size: a.product.size,
      ref_size: set.expected.end - set.expected.start + 1
    }));
  });
  // An anchor outside the size tolerance that lies at the locus of a right-sized copy is another product of one of its primers
  // there and adds nothing to that copy; elsewhere it makes a copy of its own, which the size test keeps out of the orthologs.
  const tolerance = cfgNumber(ccfg, 'genotype_ortholog_size_tolerance');
  const sized = mergeCopies(records.filter(function (r) { return sizeOk(r, tolerance); }), { zoneLength: zoneLength });
  const strays = records.filter(function (r) {
    return !sizeOk(r, tolerance) && !sized.some(function (g) { return sameLocus(g, r, zoneLength); });
  });
  const groups = sized.concat(mergeCopies(strays, { zoneLength: zoneLength }));
  const orthologs = groups.filter(function (g) { return isOrthologous(g, ccfg); });
  const called = function (gs) { return gs.filter(function (g) { return g.copy.call !== 'missing'; }); };

  let used = orthologs;
  let source = 'amplicon';
  let reason = null;
  let megablast = null;
  if (!called(orthologs).length) {
    let uncovered = orthologs.length > 0;
    if (typeof d.megablast === 'function') {
      const mb = await d.megablast();
      megablast = mb && mb.status ? mb.status : 'failed';
      if (megablast === 'ok') {
        const res = megablastCopies(mb.rows, mb.query, prepared, ccfg);
        if (called(res.copies).length) {
          used = res.copies;
          source = 'megablast';
        } else {
          uncovered = uncovered || res.kept > 0;
        }
      } else {
        reason = megablast === 'budget' ? 'fallback_budget' : 'fallback_failed';
      }
    }
    if (source !== 'megablast' && reason === null) reason = uncovered ? 'variant_not_covered' : 'no_orthologous_copy';
  }

  const calls = called(used).map(function (g) { return g.copy; });
  let allele = 'missing';
  let observed = null;
  if (calls.length) {
    allele = calls.every(function (c) { return c.call === calls[0].call; }) ? calls[0].call : 'ambiguous';
    if (allele !== 'ambiguous' && calls.every(function (c) { return c.observed === calls[0].observed; })) observed = calls[0].observed;
    reason = null;
  }
  const copies = used.map(function (g) { return g.copy; });
  const genome = {
    system_name: input.system_name,
    display_name: input.display_name != null ? input.display_name : input.system_name,
    is_reference: input.is_reference === true,
    allele: allele,
    observed: observed,
    source: allele === 'missing' ? null : source,
    copies: copies.slice(0, cfgNumber(ccfg, 'genotype_max_copies')),
    orthologous_copies: used.length,
    paralog_copies: groups.filter(function (grp) {
      return orthologs.indexOf(grp) < 0 && 100 * grp.best.matches >= PARALOG_MIN_IDENTITY * grp.best.gc_columns;
    }).length,
    reason: allele === 'missing' ? reason : null
  };
  const maxMismatches = cfgInteger(ccfg, 'genotype_offlocus_max_mismatches');
  const weakOffTargets = [];
  const sets = prepared.sets.map(function (set, si) {
    const weak = [];
    const call = predict(set, products[si], { allele: allele, copies: copies }, input.params, { max_offlocus_mismatches: maxMismatches, weak: weak });
    const asSide = set.orientation === 'forward' ? 'left' : 'right';
    weak.forEach(function (w) {
      const declared = (set.deliberate_mismatch || {})[w.which === 'ref' ? 'as_ref' : 'as_alt'];
      const mm = { left: w.product.left_mm, right: w.product.right_mm };
      mm[asSide] = allelePrimerMismatches(w.product, asSide, declared);
      weakOffTargets.push({ set: si, set_id: set.id, pair_id: w.which === 'ref' ? set.ref_pair : set.alt_pair, which: w.which, product: w.product, left_mm: mm.left, right_mm: mm.right });
    });
    return Object.assign({ system_name: input.system_name }, call);
  });
  return { genome: genome, sets: sets, megablast: megablast, anchors: anchors.length, failed_reads: failedReads, weak_off_targets: weakOffTargets };
}

// WEAK_OFF_TARGETS (M8b): the off-locus products left out of the predictions because a primer has more than
// genotype_offlocus_max_mismatches mismatches, over a whole job. details: {count, max_mismatches, examples [{system_name, set_id,
// pair_id, allele, region, start, end, size, orientation, left_mm, right_mm}]}, mismatches as the rule counts them.
function weakOffTargetDetails(ccfg) {
  return { count: 0, max_mismatches: cfgInteger(ccfg, 'genotype_offlocus_max_mismatches'), examples: [] };
}

// Adds one genome's callGenome weak_off_targets to details (mutated and returned). The examples kept are the first
// WEAK_OFF_TARGET_EXAMPLES in (genome, set, REF before ALT, region, start) order, where order.genomes lists the system names (the
// reference first) and order.sets the set ids, so they do not depend on the order in which genomes finish.
function addWeakOffTargets(details, systemName, items, order) {
  details.count += items.length;
  const rank = function (e) {
    return [order.genomes.indexOf(e.system_name), order.sets.indexOf(e.set_id), e.allele === 'ref' ? 0 : 1];
  };
  const examples = details.examples.concat(items.map(function (w) {
    const p = w.product;
    return { system_name: systemName, set_id: w.set_id, pair_id: w.pair_id, allele: w.which, region: p.region, start: p.start, end: p.end,
      size: p.size, orientation: p.orientation, left_mm: w.left_mm, right_mm: w.right_mm };
  }));
  examples.sort(function (a, b) {
    const ra = rank(a);
    const rb = rank(b);
    for (let k = 0; k < ra.length; k++) if (ra[k] !== rb[k]) return ra[k] - rb[k];
    if (a.region !== b.region) return String(a.region) < String(b.region) ? -1 : 1;
    return a.start - b.start;
  });
  details.examples = examples.slice(0, WEAK_OFF_TARGET_EXAMPLES);
  return details;
}

// A genome the check could not search (db_unavailable, blast_error) or whose call threw (call_failed): allele unavailable and
// every prediction unknown. info: {system_name, display_name, is_reference}
function unavailableEntry(prepared, info, reason) {
  const o = info || {};
  return {
    genome: {
      system_name: o.system_name,
      display_name: o.display_name != null ? o.display_name : o.system_name,
      is_reference: o.is_reference === true,
      allele: 'unavailable',
      observed: null,
      source: null,
      copies: [],
      orthologous_copies: 0,
      paralog_copies: 0,
      reason: GENOME_REASONS.indexOf(reason) >= 0 ? reason : 'call_failed'
    },
    sets: prepared.sets.map(function () { return Object.assign({ system_name: o.system_name }, unknownPrediction()); })
  };
}

// ---- results.genotyping (§2.12, §5.8) ---------------------------------------------------------------------------------------

function genomeSummary(genomes) {
  const s = { genomes_total: 0, ref: 0, alt: 0, other: 0, ambiguous: 0, missing: 0, unavailable: 0 };
  genomes.forEach(function (g) {
    s.genomes_total++;
    s[g.allele]++;
  });
  return s;
}

function setSummary(calls) {
  const s = { genomes_total: 0, predicted_ref: 0, predicted_alt: 0, both: 0, none: 0, no_call: 0, unknown: 0, weak: 0, agree: 0, disagree: 0, not_comparable: 0 };
  calls.forEach(function (c) {
    s.genomes_total++;
    s[c.predicted === 'ref' || c.predicted === 'alt' ? 'predicted_' + c.predicted : c.predicted]++;
    if (c.strength === 'weak') s.weak++;
    s[c.agrees === true ? 'agree' : c.agrees === false ? 'disagree' : 'not_comparable']++;
  });
  return s;
}

// A pair's reference verdict is consistent with the allele when a blocked (or absent) allele-specific primer goes with
// on_target_missing and any other status with an amplifying verdict (§5.5).
function consistentWithAllele(verdict, call) {
  if (!call) return false;
  return call.status === 'blocked' || call.status === 'no_product' ? verdict === 'on_target_missing'
    : verdict === 'specific' || verdict === 'off_targets';
}

// §5.8: pass when the reference is ref and the set predicts ref on it; warn when that prediction is weak or off-locus products
// appear; fail otherwise. reasons: the control's own codes, then the set call's reasons.
function referenceControl(genome, call) {
  const own = [];
  if (genome.allele !== 'ref') own.push('allele_not_ref');
  if (call.predicted !== 'ref') own.push('prediction_not_ref');
  let status = own.length ? 'fail' : 'pass';
  if (status === 'pass') {
    if (call.strength === 'weak') own.push('weak');
    if (call.off_locus_products > 0) own.push('off_locus_products');
    if (own.length) status = 'warn';
  }
  return { status: status, allele: genome.allele, reasons: own.concat(call.reasons) };
}

function publicVariant(v) {
  return {
    key: v.key,
    region: v.region,
    position: v.position,
    ref: v.ref,
    alt: v.alt,
    shift: v.shift,
    zone: { start: v.zone.start, end: v.zone.end },
    flank: v.flank,
    core: { ref: v.core.ref, alt: v.core.alt },
    haplotypes: { ref: v.haplotypes.ref, alt: v.haplotypes.alt }
  };
}

// The block run.js creates before the reference stage; writeResults fills it in.
function emptyResults(prepared) {
  return {
    algorithm_version: prepared.algorithm_version || GENOTYPING_VERSION,
    variant: publicVariant(prepared.variant),
    summary: genomeSummary([]),
    genomes: [],
    sets: prepared.sets.map(function (s) {
      return {
        id: s.id,
        ref_pair: s.ref_pair,
        alt_pair: s.alt_pair,
        orientation: s.orientation,
        deliberate_mismatch_positions: s.deliberate_mismatch_positions.slice(),
        specificity: null,
        control: null,
        reference: null,
        summary: setSummary([]),
        genomes: []
      };
    })
  };
}

// Rebuilds block (emptyResults) in place for the next flush. input: {reference: entry | null, genomes: [entry] (finished
// pan-genome genomes in request order), specificity: [{ref_verdict, alt_verdict, off_target_count}] per set | null}.
// Summaries count pan-genome genomes only. -> block
function writeResults(block, input) {
  const o = input || {};
  const reference = o.reference || null;
  const entries = o.genomes || [];
  block.summary = genomeSummary(entries.map(function (e) { return e.genome; }));
  block.genomes = (reference ? [reference.genome] : []).concat(entries.map(function (e) { return e.genome; }));
  block.sets.forEach(function (s, si) {
    const call = reference ? reference.sets[si] : null;
    const spec = o.specificity ? o.specificity[si] : null;
    s.specificity = call && spec ? {
      ref_pair: { verdict: spec.ref_verdict, consistent_with_allele: consistentWithAllele(spec.ref_verdict, call.ref_primer) },
      alt_pair: { verdict: spec.alt_verdict, consistent_with_allele: consistentWithAllele(spec.alt_verdict, call.alt_primer) },
      off_target_count: spec.off_target_count
    } : null;
    s.control = call ? referenceControl(reference.genome, call) : null;
    s.reference = call;
    s.genomes = entries.map(function (e) { return e.sets[si]; });
    s.summary = setSummary(s.genomes);
  });
  return block;
}

module.exports = {
  GENOTYPING_VERSION,
  MODES,
  MAX_SETS,
  WINDOW_PAD,
  REASONS,
  CALLER_DEFAULTS,
  ALLELES,
  COPY_CALLS,
  STATUSES,
  PREDICTIONS,
  GENOME_REASONS,
  SET_REASONS,
  CONTROL_STATUSES,
  CONTROL_REASONS,
  validateShape,
  validateLinks,
  windowFor,
  validateSets,
  prepare,
  requestBlock,
  fillSequences,
  megablastQuery,
  emptyResults,
  callGenome,
  semiGlobal,
  mergeCopies,
  gapCompressedIdentity,
  megablastCopies,
  predict,
  writeResults,
  unavailableEntry,
  weakOffTargetDetails,
  addWeakOffTargets,
  WEAK_OFF_TARGET_EXAMPLES,
  _internal: { readCore, readAlignment, opsStats, percent, anchorsOf, isOrthologous, hspAlignment, residualOf, commonStatus, referenceControl }
};
