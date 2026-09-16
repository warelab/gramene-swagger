'use strict';

// Genotyping sets (spec §4.8-§4.18), pure: the §4.8 pair filters, ALT primer derivation, discrimination with the check's
// own re-alignment and classifier, known variants under each primer, issues, quality, the §4.16 score and ranking, set
// keys, check fragments and the proposed check request. Primer3 and ntthal results arrive as arguments (scoring.js).
//
// Coordinates: template t is 1-based on the plus strand (genotyping/template.js), genomic g = template.start + t - 1.
// alt_seq indexes like seq up to the variant and is shifted by alt_offset after it. Every reported number is rounded
// once from exact decimals (decimal.js, §2.1); every threshold test compares exact values.

const crypto = require('crypto');
const guard = require('./guard');
const decimal = require('./decimal');
const mismatch = require('./mismatch');
const order = require('./order');
const normalize = require('../variation/normalize');
const realign = require('../check/realign');
const classify = require('../check/classify');
const { revcomp } = require('../sequence');

const ROLES = Object.freeze(['as_ref', 'as_alt', 'common']);
const HAPLOTYPE = Object.freeze({ as_ref: 'ref', as_alt: 'alt', common: 'both' });
// §2.15: AS_TM_IMBALANCE is high above 2.0 °C (warn above genotyping.as_tm_diff_warn).
const AS_TM_DIFF_HIGH = '2.0';
// §4.16 S6: priced by S3 and S4, so never counted as issues.
const UNCOUNTED = new Set(['AS_TM_IMBALANCE', 'COMMON_TM_OUT_OF_RANGE']);
// classify.classifyAmplicon partner for a single-primer question: a perfect site.
const PERFECT_PARTNER = Object.freeze({ mm: 0, mm_3p: 0, terminal_mm: false, gaps: 0, mm_pos: [] });
// The six tailed pairings of §4.14, in thermo.tailed order.
const PAIRINGS = Object.freeze([
  ['ref_alt_any_th', 'as_ref', 'as_alt', 'any'], ['ref_alt_end_th', 'as_ref', 'as_alt', 'end'],
  ['ref_common_any_th', 'as_ref', 'common', 'any'], ['ref_common_end_th', 'as_ref', 'common', 'end'],
  ['alt_common_any_th', 'as_alt', 'common', 'any'], ['alt_common_end_th', 'as_alt', 'common', 'end']
]);
const STRUCTURE_METRICS = Object.freeze(['hairpin_th', 'self_any_th', 'self_end_th']);
const METRIC_LABEL = Object.freeze({ hairpin_th: 'hairpin', self_any_th: 'self-dimer', self_end_th: 'self-dimer (3′ end)' });
const DEG = ' °C';

// ---- neighbours under a primer (§4.15) ----------------------------------------------------------------------------------

// `count` reference coordinates under a primer, from its 3' base at p3: a plus-strand primer reads leftwards, a
// minus-strand primer rightwards.
function coordsFrom3p(p3, count, strand) {
  const out = [];
  for (let k = 0; k < count; k++) out.push(strand === 1 ? p3 - k : p3 + k);
  return out;
}

function neighbourHit(n, distance) {
  return {
    key: n.key,
    ids: (n.ids || []).slice(),
    label: n.label,
    start: n.minimal.start,
    end: n.minimal.end,
    alleles: (n.alleles || [n.minimal.ref, n.minimal.alt]).join('/'),
    ems: n.ems === true,
    distance_from_3p: distance
  };
}

function byDistance(x, y) {
  return x.distance_from_3p - y.distance_from_3p || (x.key < y.key ? -1 : x.key > y.key ? 1 : 0);
}

// Known variants under reference coordinates listed from a primer's 3' end, as PrimerNeighbourHit objects:
// distance_from_3p is 1 + the index of the first listed coordinate inside the neighbour's minimal span (an insertion
// covers both flanking bases). Sorted by distance, then key. neighbours: canonical entries; targetKey is skipped;
// coords may hold null (an inserted base, which still takes a position along the primer).
function neighbourHits(targetKey, neighbours, coords) {
  const real = coords.filter(function (c) { return c !== null && c !== undefined; });
  if (real.length === 0) return [];
  const lo = Math.min.apply(null, real);
  const hi = Math.max.apply(null, real);
  const hits = [];
  (neighbours || []).forEach(function (n) {
    if (!n || n.key === targetKey || !n.minimal) return;
    const a = Math.min(n.minimal.start, n.minimal.end);
    const b = Math.max(n.minimal.start, n.minimal.end);
    if (b < lo || a > hi) return;
    const i = coords.findIndex(function (c) { return c !== null && c !== undefined && c >= a && c <= b; });
    if (i >= 0) hits.push(neighbourHit(n, i + 1));
  });
  return hits.sort(byDistance);
}

function nonEms(hits) {
  return hits.filter(function (h) { return !h.ems; });
}

// ---- §4.8 steps 1-4 -------------------------------------------------------------------------------------------------------

// §4.8 steps 1-3 for one design-run pair (boulder.extractPairs, template coordinates).
//   s: {orientation, as_3p (template position of the forced 3' end), zone (template), delta (alt_offset),
//       template_start, target_key, neighbours (entries), window (neighbour_3p_window), block (reject on a 3' neighbour)}
// -> {reject: 'force' | 'overlap' | 'common_neighbour_3p'} | {candidate}
//   candidate: {orientation, pair_index, pair, as, common (the pair's oligos), common_neighbours_3p (non-EMS neighbours
//               inside the common primer's 3' window; present when they did not block, e.g. for an EMS target)}
function screenPair(pair, s) {
  const forward = s.orientation === 'forward';
  const as = forward ? pair.left : pair.right;
  const common = forward ? pair.right : pair.left;
  if ((forward ? as.end : as.start) !== s.as_3p) return { reject: 'force' };
  const set = { orientation: s.orientation, as: { start: as.start, end: as.end }, common: { start: common.start, end: common.end } };
  if (!guard.commonPrimerOk(set, s.zone, { delta: s.delta })) return { reject: 'overlap' };
  // The common primer is the minus-strand right primer of a forward set (3' base at its start) and the plus-strand left
  // primer of a reverse set (3' base at its end).
  const p3 = s.template_start - 1 + (forward ? common.start : common.end);
  const hits = nonEms(neighbourHits(s.target_key, s.neighbours, coordsFrom3p(p3, s.window, forward ? -1 : 1)));
  if (hits.length > 0 && s.block) return { reject: 'common_neighbour_3p' };
  return {
    candidate: { orientation: s.orientation, pair_index: pair.rank, pair: pair, as: as, common: common, common_neighbours_3p: hits }
  };
}

// candidatesFromPairs(pairs, s, seen) -> one entry per pair, in Primer3 rank order: {reject} or {candidate}, after the
// force guard, the zone guard on both haplotypes, the common-primer 3' neighbour and duplicates of an (allele-specific,
// common) sequence pair already seen in this orientation (`seen`, shared across ladder levels). The scoring cap and the
// budget reservation (steps 5-6) belong to the orientation runner.
function candidatesFromPairs(pairs, s, seen) {
  const set = seen || new Set();
  return pairs.map(function (pair) {
    const r = screenPair(pair, s);
    if (r.reject) return r;
    const dup = r.candidate.as.seq + '|' + r.candidate.common.seq;
    if (set.has(dup)) return { reject: 'duplicate' };
    set.add(dup);
    return r;
  });
}

// ---- §4.9 ALT derivation and products ------------------------------------------------------------------------------------

// Ascending reference blocks of per-base coordinates (null = inserted base). -> {blocks, inserted_bases}
function blocksOf(coords) {
  const blocks = [];
  let inserted = 0;
  coords.forEach(function (c) {
    if (c === null || c === undefined) {
      inserted++;
      return;
    }
    const last = blocks[blocks.length - 1];
    if (last && last.end + 1 === c) last.end = c;
    else blocks.push({ start: c, end: c });
  });
  return { blocks: blocks, inserted_bases: inserted };
}

function genomicOf(region, coords, strand) {
  const b = blocksOf(coords);
  return {
    genomic: {
      region: region,
      start: b.blocks.length ? b.blocks[0].start : null,
      end: b.blocks.length ? b.blocks[b.blocks.length - 1].end : null,
      strand: strand,
      blocks: b.blocks
    },
    inserted_bases: b.inserted_bases
  };
}

function altCoords(variant, vt, start, end) {
  const out = [];
  for (let t = start; t <= end; t++) out.push(normalize.altToRef(variant.vcf, vt.start + t - 1));
  return out;
}

// deriveAlt(candidate, vt, variant) -> the ALT allele-specific primer of a candidate (§4.9): the same 5' anchor and
// length as the REF primer, its 3' base at the first haplotype difference. Forward: the same template span on alt_seq;
// reverse: the span shifted by alt_offset.
//   -> {start, end (alt_seq span), matched_seq, coords (reference coordinate per base, template order),
//       genomic {region, start, end, strand, blocks}, inserted_bases}
function deriveAlt(candidate, vt, variant) {
  const forward = candidate.orientation === 'forward';
  const delta = vt.features.alt_offset;
  const start = candidate.as.start + (forward ? 0 : delta);
  const end = candidate.as.end + (forward ? 0 : delta);
  if (start < 1 || end > vt.alt_seq.length) throw new RangeError('the ALT primer footprint lies outside the ALT template');
  const plus = vt.alt_seq.slice(start - 1, end);
  const coords = altCoords(variant, vt, start, end);
  const g = genomicOf(variant.region, coords, forward ? 1 : -1);
  return { start: start, end: end, matched_seq: forward ? plus : revcomp(plus), coords: coords, genomic: g.genomic, inserted_bases: g.inserted_bases };
}

// productSpans(candidate, vt, variant, altSize) -> {ref, alt} PrimerGenotypingProduct. altSize: the ALT check_primers
// run's product size (defaults to the ALT span length, which equals ref size + alt_offset).
function productSpans(candidate, vt, variant, altSize) {
  const forward = candidate.orientation === 'forward';
  const delta = vt.features.alt_offset;
  const as = candidate.as;
  const common = candidate.common;
  const ref = forward ? { start: as.start, end: common.end } : { start: common.start, end: as.end };
  const alt = forward ? { start: as.start, end: common.end + delta } : { start: common.start, end: as.end + delta };
  const gs = vt.start + ref.start - 1;
  const ge = vt.start + ref.end - 1;
  const altGenomic = genomicOf(variant.region, altCoords(variant, vt, alt.start, alt.end), 1);
  return {
    ref: {
      size: candidate.pair.product_size,
      template: { start: ref.start, end: ref.end, sequence: 'ref' },
      genomic: { region: variant.region, start: gs, end: ge, strand: 1, blocks: [{ start: gs, end: ge }] },
      inserted_bases: 0
    },
    alt: {
      size: altSize === undefined || altSize === null ? alt.end - alt.start + 1 : altSize,
      template: { start: alt.start, end: alt.end, sequence: 'alt' },
      genomic: altGenomic.genomic,
      inserted_bases: altGenomic.inserted_bases
    }
  };
}

// ---- §4.12 discrimination ---------------------------------------------------------------------------------------------------

// The check's own 3'-anchored alignment of `primer` on one haplotype window [fp.start - 3, fp.end + 3], classified against
// a perfect partner with the check's default params. -> PrimerLikelihoodAt
function alignAt(primer, orientation, hap, fp) {
  const forward = orientation === 'forward';
  const lo = Math.max(1, fp.start - realign.PAD);
  const hi = Math.min(hap.length, fp.end + realign.PAD);
  const aln = realign.realignSiteOnWindow({ key: 'discrimination', primer: primer, face: forward ? 'F' : 'R', p3: forward ? fp.end : fp.start },
    hap.slice(lo - 1, hi), lo, classify.DEFAULT_PARAMS.three_prime_window);
  return { mm_pos: aln.mm_pos.slice(), likelihood: classify.classifyAmplicon(aln, PERFECT_PARTNER) };
}

// in_shift_tract of both allele-specific primers (§3.5, §4.12), from the anchor-excluded shift tract: as_ref when its 3'
// position lies inside; as_alt when alt_maps_to lies inside or is an inserted base (null) while the variant has a tract.
function shiftTractFlags(variant, orientation) {
  const tract = normalize.tractOf(variant.vcf, variant.shift);
  const d = variant.discriminating[orientation];
  const inside = function (p) { return tract !== null && p !== null && p !== undefined && p >= tract.start && p <= tract.end; };
  return { as_ref: inside(d.position), as_alt: tract !== null && (d.alt_maps_to === null || inside(d.alt_maps_to)) };
}

// discrimination({orientation, variant, template, as_ref, as_alt}) -> {as_ref, as_alt} PrimerGenotypingDiscrimination
//   as_ref / as_alt: {target_seq (as ordered), matched_seq, start, end (own haplotype span)}
// The other haplotype's span: forward the same indices; reverse shifted by +alt_offset (REF primer on alt_seq) or
// -alt_offset (ALT primer on seq).
function discrimination(input) {
  const o = input.orientation;
  const vt = input.template;
  const v = input.variant;
  const delta = vt.features.alt_offset;
  const forward = o === 'forward';
  const d = v.discriminating[o];
  const flags = shiftTractFlags(v, o);
  const ref = input.as_ref;
  const alt = input.as_alt;
  const refOther = forward ? { start: ref.start, end: ref.end } : { start: ref.start + delta, end: ref.end + delta };
  const altOther = forward ? { start: alt.start, end: alt.end } : { start: alt.start - delta, end: alt.end - delta };
  return {
    as_ref: {
      own_allele: alignAt(ref.target_seq, o, vt.seq, ref),
      other_allele: alignAt(ref.target_seq, o, vt.alt_seq, refOther),
      terminal_mismatch_class: mismatch.terminalClass(ref.matched_seq, o, d.alt_base),
      in_shift_tract: flags.as_ref
    },
    as_alt: {
      own_allele: alignAt(alt.target_seq, o, vt.alt_seq, alt),
      other_allele: alignAt(alt.target_seq, o, vt.seq, altOther),
      terminal_mismatch_class: mismatch.terminalClass(alt.matched_seq, o, d.ref_base),
      in_shift_tract: flags.as_alt
    }
  };
}

// ---- §4.15 neighbours of the three primers ---------------------------------------------------------------------------------

// primerNeighbours({orientation, variant, template, as, alt (deriveAlt), common, neighbours}) -> {as_ref, as_alt, common}
// hit lists, each measured along its own primer from its own 3' end: the REF allele-specific primer over its REF
// footprint, the ALT one over its mapped ALT bases (inserted bases count), the common primer over its footprint.
function primerNeighbours(input) {
  const forward = input.orientation === 'forward';
  const vt = input.template;
  const key = input.variant.key;
  const g = function (t) { return vt.start + t - 1; };
  const as = input.as;
  const common = input.common;
  const asCoords = forward ? coordsFrom3p(g(as.end), as.end - as.start + 1, 1) : coordsFrom3p(g(as.start), as.end - as.start + 1, -1);
  const altCoords3 = forward ? input.alt.coords.slice().reverse() : input.alt.coords.slice();
  const commonCoords = forward ? coordsFrom3p(g(common.start), common.end - common.start + 1, -1)
    : coordsFrom3p(g(common.end), common.end - common.start + 1, 1);
  return {
    as_ref: neighbourHits(key, input.neighbours, asCoords),
    as_alt: neighbourHits(key, input.neighbours, altCoords3),
    common: neighbourHits(key, input.neighbours, commonCoords)
  };
}

// Distinct non-EMS neighbour keys under the three primers.
function neighbourSites(hits) {
  const keys = new Set();
  ROLES.forEach(function (r) { nonEms(hits[r]).forEach(function (h) { keys.add(h.key); }); });
  return keys.size;
}

// ---- §4.13 Tm balance ---------------------------------------------------------------------------------------------------------

// tmBalance({as_ref_tm, as_alt_tm, as_ref_nominal, as_alt_nominal, common_tm}, g) -> exact micro units
//   {as_tm_diff = |as_ref.tm - as_alt.tm|, common_minus_as = common.tm - max(nominal Tm), common_out = distance outside
//    [common_tm_low, common_tm_high]}. tm is the reported Tm (the duplex Tm for a mismatch primer); nominal is Primer3's.
function tmBalance(t, g) {
  const diff = decimal.abs(decimal.micro(t.as_ref_tm) - decimal.micro(t.as_alt_tm));
  const cma = decimal.micro(t.common_tm) - decimal.max(decimal.micro(t.as_ref_nominal), decimal.micro(t.as_alt_nominal));
  const low = decimal.micro(g.common_tm_low);
  const high = decimal.micro(g.common_tm_high);
  return { as_tm_diff: diff, common_minus_as: cma, common_out: cma < low ? low - cma : cma > high ? cma - high : 0n };
}

// ---- §2.15 issues ----------------------------------------------------------------------------------------------------------------

function fixed2(v) {
  return decimal.round(v, 2).toFixed(2);
}

function thresholdText(x) {
  const n = Number(x);
  return Number.isInteger(n) ? String(n) : String(x);
}

function signedText(x) {
  const n = decimal.round(decimal.micro(x), 1);
  return (n < 0 ? '−' : '+') + Math.abs(n).toFixed(1);
}

// " GC content too low; Temperature too low;" -> "GC content too low, temperature too low"
function problemPhrases(problems) {
  return String(problems).split(';').map(function (p) { return p.trim(); }).filter(Boolean).map(function (p) {
    return /^[A-Z][a-z]/.test(p) ? p.charAt(0).toLowerCase() + p.slice(1) : p;
  }).join(', ');
}

function whereOf(role) {
  return role === 'common' ? 'the common primer' : 'the allele-specific primers';
}

function neighbourList(hits) {
  return hits.map(function (h) { return (h.ids[0] || h.key) + ' (' + h.distance_from_3p + ' nt)'; }).join(', ');
}

// NEIGHBOUR_IN_PRIMER (outside the 3' window, warn) and NEIGHBOUR_AT_3P (inside it, high) for one primer site.
function neighbourIssues(role, hits, window, add) {
  const outside = hits.filter(function (h) { return h.distance_from_3p > window; });
  const inside = hits.filter(function (h) { return h.distance_from_3p <= window; });
  const details = function (list) {
    return { role: role, ids: list.map(function (h) { return h.ids[0] || h.key; }), distances: list.map(function (h) { return h.distance_from_3p; }) };
  };
  if (outside.length) {
    const one = outside.length === 1;
    const d = details(outside);
    add('NEIGHBOUR_IN_PRIMER', 'warn',
      one ? d.ids[0] + ' in ' + whereOf(role) + ', ' + d.distances[0] + ' nt from its 3′ end'
        : outside.length + ' known variants in ' + whereOf(role) + ': ' + neighbourList(outside) + ' from the 3′ end',
      one ? 'known variant ' + d.ids[0] + ' lies in ' + whereOf(role) + ', ' + d.distances[0] + ' nt from its 3′ end'
        : 'known variants lie in ' + whereOf(role) + ': ' + neighbourList(outside) + ' from the 3′ end',
      d);
  }
  if (inside.length) {
    const d = details(inside);
    add('NEIGHBOUR_AT_3P', 'high',
      neighbourList(inside) + ' inside the ' + window + '-nt 3′ window of ' + whereOf(role),
      'known variant' + (inside.length === 1 ? ' ' : 's ') + neighbourList(inside) + ' inside the last ' + window + ' nt of ' + whereOf(role) +
        ' would weaken or silence priming on genomes that carry ' + (inside.length === 1 ? 'it' : 'them'),
      d);
  }
}

// issuesFor(x) -> the issues of one set, in the §2.15 order, each {code, severity, message (issue text), warning
// (warning text), details, warning_details?, value_exact?}.
//   x: {g (primers.genotyping), variant, orientation,
//       problems {as_ref, as_alt}: echoed PRIMER_*_PROBLEMS of the ordered derived oligos (as_ref only with a mismatch),
//       balance (tmBalance), dyes {as_ref, as_alt},
//       tailed {as_ref, as_alt: {hairpin_th, self_any_th, self_end_th}} | null, tailed_pairs {<PAIRINGS key>} | null
//         (exact ntthal decimals),
//       mismatch_not_applicable: the position | null,
//       mismatch_structures {as_ref, as_alt: {hairpin_th, self_any_th, self_end_th}} | null (the mismatch runs),
//       neighbours {as_ref, as_alt, common} (primerNeighbours), discrimination {as_ref, as_alt}}
function issuesFor(x) {
  const g = x.g;
  const out = [];
  const warnTh = decimal.micro(g.structure_warn_th);
  const highTh = decimal.micro(g.structure_high_th);
  const add = function (code, severity, message, warning, details, extra) {
    out.push(Object.assign({ code: code, severity: severity, message: message, warning: warning, details: details }, extra || {}));
  };

  ['as_alt', 'as_ref'].forEach(function (role) {
    const p = x.problems && x.problems[role];
    if (typeof p !== 'string' || p.trim() === '') return;
    const phrases = problemPhrases(p);
    add('ALT_PRIMER_SUBOPTIMAL', 'warn', role + ': Primer3 reports ' + phrases,
      'the derived ' + role + ' primer would not have been chosen de novo (Primer3: ' + phrases + '); it clears the hard floors',
      { oligo: role, problems: p });
  });

  const b = x.balance;
  if (b.as_tm_diff > decimal.micro(g.as_tm_diff_warn)) {
    const diff = decimal.round(b.as_tm_diff, 2);
    add('AS_TM_IMBALANCE', b.as_tm_diff > decimal.micro(AS_TM_DIFF_HIGH) ? 'high' : 'warn',
      'as_ref and as_alt Tm differ by ' + fixed2(b.as_tm_diff) + DEG,
      'as_ref and as_alt Tm differ by ' + fixed2(b.as_tm_diff) + DEG + ' (limit ' + Number(g.as_tm_diff_warn).toFixed(1) + ')',
      { diff: diff });
  }
  if (b.common_out > 0n) {
    const dir = b.common_minus_as < 0n ? 'below' : 'above';
    const abs = fixed2(decimal.abs(b.common_minus_as));
    add('COMMON_TM_OUT_OF_RANGE', 'warn', 'common primer Tm ' + abs + DEG + ' ' + dir + ' the warmer allele-specific primer',
      'the common primer Tm is ' + abs + DEG + ' ' + dir + ' the warmer allele-specific primer (wanted ' + signedText(g.common_tm_low) +
        ' to ' + signedText(g.common_tm_high) + DEG + ')',
      { value: decimal.round(b.common_minus_as, 2) });
  }

  if (x.tailed) {
    ['as_ref', 'as_alt'].forEach(function (role) {
      const t = x.tailed[role];
      const dye = x.dyes[role];
      STRUCTURE_METRICS.forEach(function (metric) {
        const v = decimal.micro(t[metric]);
        if (v < warnTh) return;
        const high = metric === 'hairpin_th' && v >= highTh;
        const severity = high ? 'high' : 'warn';
        const value = decimal.round(v, 2);
        const text = metric === 'hairpin_th'
          ? dye + '-tailed ' + role + ': hairpin Tm ' + fixed2(v) + DEG + ' is ' + (high
            ? 'at or above ' + thresholdText(g.structure_high_th) + DEG + ', the final KASP annealing temperature'
            : 'above ' + thresholdText(g.structure_warn_th) + DEG)
          : dye + '-tailed ' + role + ': ' + METRIC_LABEL[metric] + ' Tm ' + fixed2(v) + DEG + ' is at or above ' + thresholdText(g.structure_warn_th) + DEG;
        add('TAILED_STRUCTURE', severity, dye + '-tailed ' + role + ' ' + METRIC_LABEL[metric] + ' ' + fixed2(v) + DEG, text,
          { oligo: role, metric: metric, value: value },
          { warning_details: { oligo: role, metric: metric, value: value, severity: severity }, value_exact: t[metric] });
      });
    });
    PAIRINGS.forEach(function (p) {
      const v = decimal.micro(x.tailed_pairs[p[0]]);
      if (v < warnTh) return;
      const value = decimal.round(v, 2);
      const what = (p[1] === 'common' ? 'common' : x.dyes[p[1]] + '-tailed ' + p[1]) + ' × ' +
        (p[2] === 'common' ? 'common' : x.dyes[p[2]] + '-tailed ' + p[2]) + ' dimer' + (p[3] === 'end' ? ' (3′ end)' : '');
      add('TAILED_STRUCTURE', 'warn', what + ' ' + fixed2(v) + DEG,
        what + ': Tm ' + fixed2(v) + DEG + ' is at or above ' + thresholdText(g.structure_warn_th) + DEG,
        { pair: [p[1], p[2]], metric: p[0], value: value },
        { warning_details: { pair: [p[1], p[2]], metric: p[0], value: value, severity: 'warn' }, value_exact: x.tailed_pairs[p[0]] });
    });
  }

  if (x.mismatch_not_applicable) {
    add('MISMATCH_NOT_APPLICABLE', 'warn',
      'no deliberate mismatch: the base at −' + x.mismatch_not_applicable + ' differs between the allele-specific primers',
      'the base at −' + x.mismatch_not_applicable + ' differs between the REF and ALT allele-specific primers, so no deliberate mismatch was applied',
      { position: x.mismatch_not_applicable });
  }
  if (x.mismatch_structures) {
    ['as_ref', 'as_alt'].forEach(function (role) {
      STRUCTURE_METRICS.forEach(function (metric) {
        const raw = x.mismatch_structures[role][metric];
        if (raw === null || raw === undefined) return;
        const v = decimal.micro(raw);
        if (v < warnTh) return;
        const severity = v >= highTh ? 'high' : 'warn';
        const value = decimal.round(v, 2);
        add('MISMATCH_STRUCTURE', severity, role + ' ' + metric + ' ' + fixed2(v) + DEG,
          'deliberate-mismatch ' + role + ': ' + METRIC_LABEL[metric] + ' Tm ' + fixed2(v) + DEG + ' is at or above ' +
            thresholdText(severity === 'high' ? g.structure_high_th : g.structure_warn_th) + DEG,
          { oligo: role, metric: metric, value: value },
          { warning_details: { oligo: role, metric: metric, value: value, severity: severity } });
      });
    });
  }

  // One issue per code for the allele-specific site (the union of both primers' non-EMS neighbours, each with its smaller
  // distance) and one for the common primer (§4.15).
  const window = g.neighbour_3p_window;
  const union = new Map();
  nonEms(x.neighbours.as_ref.concat(x.neighbours.as_alt)).forEach(function (h) {
    const cur = union.get(h.key);
    if (!cur || cur.distance_from_3p > h.distance_from_3p) union.set(h.key, h);
  });
  neighbourIssues('as_ref', Array.from(union.values()).sort(byDistance), window, add);
  neighbourIssues('common', nonEms(x.neighbours.common), window, add);

  ['as_ref', 'as_alt'].forEach(function (role) {
    const other = x.discrimination[role].other_allele;
    if (other.likelihood !== 'likely') return;
    add('WEAK_DISCRIMINATION', 'high', role + ' is predicted to amplify the other allele (mm_pos ' + JSON.stringify(other.mm_pos) + ')',
      'the ' + role + ' primer is predicted to amplify the other allele as well, so it cannot discriminate', { primer: role, mm_pos: other.mm_pos.slice() });
  });
  ['as_ref', 'as_alt'].forEach(function (role) {
    if (!x.discrimination[role].in_shift_tract) return;
    add('SHIFT_TRACT_DISCRIMINATION', 'warn', role + ' 3′ base lies in the ' + x.variant.shift + '-base shift tract',
      'the ' + role + ' primer\'s 3′ base lies inside the shift tract of the indel; on the other allele it may still prime through a bulge',
      { primer: role, shift: x.variant.shift });
  });
  if (x.discrimination.as_ref.terminal_mismatch_class === 'weak' && x.discrimination.as_alt.terminal_mismatch_class === 'weak') {
    add('WEAK_TERMINAL_CLASS', 'info', 'both allele-specific 3′ mismatches are in the weak class (A/G or C/T site)', null, {});
  }
  return out;
}

// ---- §4.16 score, quality and ranking ---------------------------------------------------------------------------------------

// scoreSet({penalty, relaxation_level, balance, issues, kasp, product_size, g}) -> {exact (micro units), score, quality,
// counted {high, warn}, terms}. S1 the candidate's own PRIMER_PAIR_i_PENALTY; S2 5 x level; S3 2 x max(0, as_tm_diff -
// as_tm_diff_warn); S4 2 x the distance of common_minus_as outside its range; S5 (kasp only) max(0, size - 100) / 25;
// S6 3 x counted high + 1 x counted warn issues (never the two Tm issues, never info). S7: exact sum, rounded once.
function scoreSet(t) {
  const g = t.g;
  const counted = t.issues.filter(function (i) { return i.severity !== 'info' && !UNCOUNTED.has(i.code); });
  const high = counted.filter(function (i) { return i.severity === 'high'; }).length;
  const warn = counted.filter(function (i) { return i.severity === 'warn'; }).length;
  const terms = {
    penalty: decimal.micro(t.penalty),
    relaxation: decimal.fromInt(5 * t.relaxation_level),
    as_tm: 2n * decimal.max(0n, t.balance.as_tm_diff - decimal.micro(g.as_tm_diff_warn)),
    common_tm: 2n * t.balance.common_out,
    product: t.kasp ? BigInt(Math.max(0, t.product_size - 100)) * (decimal.UNIT / 25n) : 0n,
    issues: decimal.fromInt(3 * high + warn)
  };
  const exact = Object.keys(terms).reduce(function (sum, k) { return sum + terms[k]; }, 0n);
  const quality = t.issues.some(function (i) { return i.severity === 'high'; }) ? 'poor'
    : t.issues.some(function (i) { return i.severity === 'warn'; }) || t.relaxation_level >= 1 ? 'usable' : 'good';
  return { exact: exact, score: decimal.round(exact, 2), quality: quality, counted: { high: high, warn: warn }, terms: terms };
}

// The exact score, then the key.
function compareSets(a, b) {
  return a.exact < b.exact ? -1 : a.exact > b.exact ? 1 : a.key < b.key ? -1 : a.key > b.key ? 1 : 0;
}

// rankSets(sets, numSets) -> the returned sets in rank order: the best set of each orientation first (ordered between
// themselves by compareSets), then every remaining set by compareSets; cut to numSets. sets: [{orientation, exact, key}].
function rankSets(sets, numSets) {
  const bests = ['forward', 'reverse'].map(function (o) {
    return sets.filter(function (s) { return s.orientation === o; }).sort(compareSets)[0];
  }).filter(Boolean).sort(compareSets);
  const rest = sets.filter(function (s) { return bests.indexOf(s) < 0; }).sort(compareSets);
  return bests.concat(rest).slice(0, numSets);
}

// First 12 hex of sha256(orientation|as_ref.target_seq|as_alt.target_seq|common.target_seq) (§2.1, D5).
function setKey(orientation, refTarget, altTarget, commonTarget) {
  return crypto.createHash('sha256').update(orientation + '|' + refTarget + '|' + altTarget + '|' + commonTarget).digest('hex').slice(0, 12);
}

// ---- oligos, fragments and the check request (§4.17, §4.18) ------------------------------------------------------------------

function rounded(v) {
  return v === null || v === undefined ? null : decimal.roundValue(v, 2);
}

// oligo(fields) -> PrimerGenotypingOligo in the §2.9 key order. Numbers are exact tool values, rounded here.
//   fields: {role, allele, target_seq, matched_seq, dye, tail_seq, tm, tm_method, matched_tm, hairpin_th, self_any_th,
//            self_end_th, end_stability, primer3_problems, template, genomic, inserted_bases, deliberate_mismatch,
//            discrimination, tailed {hairpin_th, self_any_th, self_end_th} | null, neighbours}
function oligo(f) {
  const orderSeq = (f.tail_seq || '') + f.target_seq;
  return {
    role: f.role,
    allele: f.allele,
    three_prime_base: f.target_seq[f.target_seq.length - 1],
    haplotype: HAPLOTYPE[f.role],
    target_seq: f.target_seq,
    matched_seq: f.matched_seq,
    dye: f.dye || null,
    tail_seq: f.tail_seq || null,
    order_seq: orderSeq,
    len: f.target_seq.length,
    order_len: orderSeq.length,
    tm: rounded(f.tm),
    tm_method: f.tm_method,
    matched_tm: rounded(f.matched_tm),
    gc: decimal.gcPercent(f.target_seq),
    hairpin_th: rounded(f.hairpin_th),
    self_any_th: rounded(f.self_any_th),
    self_end_th: rounded(f.self_end_th),
    end_stability: rounded(f.end_stability),
    primer3_problems: typeof f.primer3_problems === 'string' && f.primer3_problems !== '' ? f.primer3_problems : null,
    template: f.template,
    genomic: f.genomic,
    inserted_bases: f.inserted_bases,
    deliberate_mismatch: f.deliberate_mismatch || null,
    discrimination: f.discrimination || null,
    tailed: f.tailed ? { hairpin_th: rounded(f.tailed.hairpin_th), self_any_th: rounded(f.tailed.self_any_th), self_end_th: rounded(f.tailed.self_end_th) } : null,
    neighbours: f.neighbours
  };
}

// checkFragment(set) -> PrimerGenotypingCheckFragment: the REF and ALT pairs of the set with the REF product as expected.
function checkFragment(set) {
  const forward = set.orientation === 'forward';
  const g = set.products.ref.genomic;
  const common = set.primers.common.target_seq;
  const pair = function (suffix, as) {
    return {
      id: set.id + '_' + suffix,
      left: forward ? as.target_seq : common,
      right: forward ? common : as.target_seq,
      expected: { region: g.region, start: g.start, end: g.end }
    };
  };
  return {
    set: { id: set.id, ref_pair: set.id + '_REF', alt_pair: set.id + '_ALT' },
    pairs: [pair('REF', set.primers.as_ref), pair('ALT', set.primers.as_alt)]
  };
}

// proposedCheckRequest(sets, {system_name, variant, max_sets, max_pairs, max_unique_primers}) -> PrimerGenotypingCheckBlock,
// or null without sets. Sets are added in rank order while the request stays inside max_sets (check_max_sets 5),
// max_pairs (check.max_pairs 10) and max_unique_primers (check_max_unique_primers 13); from the first set that does not
// fit, the rest are omitted_set_ids.
function proposedCheckRequest(sets, opts) {
  if (!sets || sets.length === 0) return null;
  const included = [];
  const omitted = [];
  const primers = new Set();
  let pairs = [];
  let open = true;
  sets.forEach(function (s) {
    if (open) {
      const next = new Set(primers);
      s.check.pairs.forEach(function (p) {
        next.add(p.left.toUpperCase());
        next.add(p.right.toUpperCase());
      });
      if (included.length < opts.max_sets && pairs.length + s.check.pairs.length <= opts.max_pairs && next.size <= opts.max_unique_primers) {
        included.push(s);
        pairs = pairs.concat(s.check.pairs.map(function (p) { return { id: p.id, left: p.left, right: p.right, expected: Object.assign({}, p.expected) }; }));
        next.forEach(function (p) { primers.add(p); });
        return;
      }
      open = false;
    }
    omitted.push(s.id);
  });
  const v = opts.variant;
  return {
    request: {
      system_name: opts.system_name,
      mode: 'region',
      checks: ['specificity', 'pangenome'],
      pairs: pairs,
      genotyping: {
        variant: { region: v.region, position: v.vcf.position, ref: v.vcf.ref, alt: v.vcf.alt },
        sets: included.map(function (s) { return Object.assign({}, s.check.set); })
      }
    },
    set_ids: included.map(function (s) { return s.id; }),
    unique_primers: primers.size,
    omitted_set_ids: omitted
  };
}

// finalizeSet(scored, rank, {label, variant}) -> PrimerGenotypingSet in the §2.9 key order: id S{rank+1}, warnings (the
// warn and high issues), issues, the check fragment and the order rows.
function finalizeSet(scored, rank, opts) {
  const set = {
    id: 'S' + (rank + 1),
    key: scored.key,
    rank: rank,
    orientation: scored.orientation,
    relaxation_level: scored.relaxation_level,
    quality: scored.quality,
    score: scored.score,
    primers: scored.primers,
    products: scored.products,
    thermo: scored.thermo,
    tm_balance: scored.tm_balance,
    neighbour_sites: scored.neighbour_sites,
    primer3_penalty: scored.primer3_penalty,
    warnings: scored.issues.filter(function (i) { return i.severity !== 'info'; }).map(function (i) {
      return { code: i.code, message: i.warning, details: i.warning_details || i.details };
    }),
    issues: scored.issues.map(function (i) { return { code: i.code, severity: i.severity, message: i.message, details: i.details }; }),
    check: null,
    order: null
  };
  set.check = checkFragment(set);
  set.order = order.orderRows({ id: set.id, key: set.key, orientation: set.orientation, primers: set.primers, products: set.products, issues: scored.issues },
    { label: opts.label, variant_key: opts.variant.key });
  return set;
}

module.exports = {
  coordsFrom3p,
  neighbourHits,
  nonEms,
  screenPair,
  candidatesFromPairs,
  blocksOf,
  deriveAlt,
  productSpans,
  alignAt,
  shiftTractFlags,
  discrimination,
  primerNeighbours,
  neighbourSites,
  tmBalance,
  problemPhrases,
  issuesFor,
  scoreSet,
  compareSets,
  rankSets,
  setKey,
  oligo,
  checkFragment,
  proposedCheckRequest,
  finalizeSet,
  PAIRINGS,
  AS_TM_DIFF_HIGH,
  UNCOUNTED
};
