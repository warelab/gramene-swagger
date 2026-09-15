'use strict';

// Scoring one genotyping candidate (spec §4.9-§4.16): the Primer3 check_primers runs of §4.10, the deliberate mismatch
// of §4.11 with its ntthal duplex Tm, the §4.6 hard floors, the tailed structures of §4.14, then discrimination,
// neighbours, issues and the score through the pure helpers of sets.js.
//
// scorePair(haplotype, left, right, params, deps) runs one check_primers record: PRIMER_TASK=check_primers,
// SEQUENCE_PRIMER / SEQUENCE_PRIMER_REVCOMP, PRIMER_PICK_ANYWAY=1 (only so Primer3 never silently returns nothing for an
// already chosen primer; the floors are ours), and the level's params with the pair-level constraints widened (sizes
// 15-36, max_tm_diff 30, product range [36, upper + 100]) while primer-level constraints keep their values, so
// PRIMER_*_PROBLEMS stays meaningful. A PRIMER_ERROR (e.g. "Specified right primer not in sequence" when the common
// primer is absent from the ALT haplotype, §4.4) or no returned pair never fails the request: the candidate is dropped
// as alt_scoring_failed. Every run is counted in the design budget when it starts.
//
// scoreCandidate(candidate, ctx) is genotyping/design.js's ctx.scoreCandidate. Order of work, which the design budget
// counts: the ALT run; with a deliberate mismatch the REF and ALT mismatch runs and their two duplex calls; the floors;
// the 15 tailed ntthal calls. -> {set: candidate (with candidate.scored)} | {dropped: 'alt_scoring_failed' | 'below_floor'}

const boulder = require('../boulder');
const design = require('../design');
const { PrimerHttpError } = require('../errors');
const { revcomp } = require('../sequence');
const decimal = require('./decimal');
const mismatch = require('./mismatch');
const sets = require('./sets');

// KASP tails (§4.14). They never reach Primer3 or the check.
const TAILS = Object.freeze({ FAM: 'GAAGGTGACCAAGTTCATGCT', HEX: 'GAAGGTCGGAGTCAACGGATT' });
const DYES = Object.freeze({
  none: Object.freeze({ as_ref: null, as_alt: null }),
  ref_fam_alt_hex: Object.freeze({ as_ref: 'FAM', as_alt: 'HEX' }),
  ref_hex_alt_fam: Object.freeze({ as_ref: 'HEX', as_alt: 'FAM' })
});
// §4.10 widened pair-level constraints of a scoring run.
const WIDENED = Object.freeze({ min_size: 15, max_size: 36, max_tm_diff: 30, product_min: 36, product_pad: 100 });

function deadlineError(ms) {
  return new PrimerHttpError(504, 'DEADLINE_EXCEEDED', 'the design did not finish within ' + ms + ' ms', { deadline_ms: ms });
}

function rounded(v) {
  return v === null || v === undefined ? null : decimal.roundValue(v, 2);
}

// The level params with the §4.10 pair-level constraints widened, one pair returned.
function scoringParams(params) {
  const upper = Math.max.apply(null, params.product_size_ranges.map(function (r) { return r[1]; }));
  return Object.assign({}, params, {
    min_size: WIDENED.min_size,
    max_size: WIDENED.max_size,
    max_tm_diff: WIDENED.max_tm_diff,
    product_size_ranges: [[WIDENED.product_min, upper + WIDENED.product_pad]],
    num_return: 1
  });
}

// The check_primers record: design.buildRecord (unmasked haplotype, no target or regions), then PRIMER_TASK replaced,
// PRIMER_PICK_ANYWAY and the two primers appended.
function scoringRecord(id, seq, left, right, params) {
  const tags = design.buildRecord({ id: id, seq: seq }, {
    params: scoringParams(params), target: null, included: null, excluded: [], junctions: null, avoid_repeats: false
  });
  tags.PRIMER_TASK = 'check_primers';
  tags.PRIMER_PICK_ANYWAY = 1;
  tags.SEQUENCE_PRIMER = left;
  tags.SEQUENCE_PRIMER_REVCOMP = right;
  return tags;
}

// scorePair(haplotype {id, seq}, left, right, params, deps {cfg, primer3, budget, deadline, log})
//   -> {pair (boulder.extractPairs pair 0, with left.problems and right.problems), tags} | {failed: 'alt_scoring_failed', error}
async function scorePair(haplotype, left, right, params, deps) {
  const tags = scoringRecord(haplotype.id, haplotype.seq, left, right, params);
  const cfg = deps.cfg;
  const dl = deps.deadline;
  const remaining = Math.floor(dl.remaining());
  if (remaining <= 0) throw deadlineError(cfg.design.deadline_ms);
  deps.budget.countPrimer3();
  const result = await dl.within(deps.primer3.run(tags, {
    timeoutMs: Math.min(cfg.design.primer3_timeout_ms, remaining), signal: dl.signal, log: deps.log, tmpDir: cfg.tmp_dir
  }));
  if (result.error) return { failed: 'alt_scoring_failed', error: result.error };
  const pairs = boulder.extractPairs(result.tags);
  if (pairs.length < 1) return { failed: 'alt_scoring_failed', error: null };
  const pair = pairs[0];
  pair.left.problems = result.tags.PRIMER_LEFT_0_PROBLEMS || null;
  pair.right.problems = result.tags.PRIMER_RIGHT_0_PROBLEMS || null;
  return { pair: pair, tags: result.tags };
}

// §4.6 hard floors, exact: every allele-specific primer's matched (Primer3) Tm >= as_min_tm; every derived primer's
// target GC >= as_min_gc. oligos: {as_ref, as_alt: {matched_tm, target_seq, derived}}. -> true when a floor fails.
function belowFloors(oligos, floors) {
  return ['as_ref', 'as_alt'].some(function (role) {
    const o = oligos[role];
    if (o.matched_tm === null || o.matched_tm === undefined || decimal.micro(o.matched_tm) < decimal.micro(floors.as_min_tm)) return true;
    return o.derived === true && !decimal.gcAtLeast(o.target_seq, floors.as_min_gc);
  });
}

function structuresOf(side) {
  return { hairpin_th: side.hairpin_th, self_any_th: side.self_any_th, self_end_th: side.self_end_th };
}

function pairThermo(pair) {
  return { compl_any_th: rounded(pair.compl_any_th), compl_end_th: rounded(pair.compl_end_th) };
}

async function scoreCandidate(candidate, ctx) {
  const req = ctx.req;
  const assay = req.assay;
  const g = ctx.cfg.genotyping;
  const vt = ctx.template;
  const v = ctx.variant;
  const o = candidate.orientation;
  const forward = o === 'forward';
  const d = v.discriminating[o];
  const as = candidate.as;
  const common = candidate.common;
  const params = candidate.params;
  const id = vt.id + '_' + o + '_L' + candidate.level + '_p' + candidate.pair_index;
  const primersWith = function (asSeq) { return forward ? [asSeq, common.seq] : [common.seq, asSeq]; };
  const asSide = function (pair) { return forward ? pair.left : pair.right; };
  const alt = sets.deriveAlt(candidate, vt, v);

  // The ALT primer with the common primer on the ALT haplotype.
  const altPrimers = primersWith(alt.matched_seq);
  const altRun = await scorePair({ id: id + '_alt', seq: vt.alt_seq }, altPrimers[0], altPrimers[1], params, ctx);
  if (altRun.failed) return { dropped: 'alt_scoring_failed' };
  const altOut = asSide(altRun.pair);

  // §4.11: one mismatch run per haplotype, on a copy whose template base opposite the mismatch is substituted, and the
  // duplex Tm against the unmodified own-allele footprint.
  let mm = null;
  let notApplicable = null;
  if (assay.deliberate_mismatch === 'auto') {
    const k = assay.mismatch_position;
    if (!mismatch.applicable(as.seq, alt.matched_seq, k)) {
      notApplicable = k;
    } else {
      const refBlock = mismatch.chooseMismatch(as.seq, o, d.alt_base, k);
      const altBlock = mismatch.chooseMismatch(alt.matched_seq, o, d.ref_base, k);
      const refSeq = mismatch.applyMismatch(as.seq, refBlock);
      const altSeq = mismatch.applyMismatch(alt.matched_seq, altBlock);
      const refPrimers = primersWith(refSeq);
      const refRun = await scorePair({ id: id + '_mmref', seq: mismatch.substituteTemplate(vt.seq, as, o, refBlock) },
        refPrimers[0], refPrimers[1], params, ctx);
      if (refRun.failed) return { dropped: 'alt_scoring_failed' };
      const altMmPrimers = primersWith(altSeq);
      const altMmRun = await scorePair({ id: id + '_mmalt', seq: mismatch.substituteTemplate(vt.alt_seq, alt, o, altBlock) },
        altMmPrimers[0], altMmPrimers[1], params, ctx);
      if (altMmRun.failed) return { dropped: 'alt_scoring_failed' };
      const duplex = await ctx.deadline.within(Promise.all([
        ctx.thermo.duplex(refSeq, revcomp(as.seq)),
        ctx.thermo.duplex(altSeq, revcomp(alt.matched_seq))
      ]));
      mm = {
        ref: { block: refBlock, seq: refSeq, run: refRun, out: asSide(refRun.pair), duplex: duplex[0] },
        alt: { block: altBlock, seq: altSeq, run: altMmRun, out: asSide(altMmRun.pair), duplex: duplex[1] }
      };
    }
  }

  const refTarget = mm ? mm.ref.seq : as.seq;
  const altTarget = mm ? mm.alt.seq : alt.matched_seq;
  if (belowFloors({
    as_ref: { matched_tm: as.tm, target_seq: refTarget, derived: mm !== null },
    as_alt: { matched_tm: altOut.tm, target_seq: altTarget, derived: true }
  }, req.floors)) {
    return { dropped: 'below_floor' };
  }

  const dyes = DYES[assay.tails] || DYES.none;
  const tailOf = function (role) { return dyes[role] ? TAILS[dyes[role]] : null; };
  let tailed = null;
  let tailedPairs = null;
  if (dyes.as_ref || dyes.as_alt) {
    const R = (tailOf('as_ref') || '') + refTarget;
    const A = (tailOf('as_alt') || '') + altTarget;
    const C = common.seq;
    const t = ctx.thermo;
    const r = await ctx.deadline.within(Promise.all([
      t.hairpin(R), t.selfAny(R), t.selfEnd(R), t.hairpin(A), t.selfAny(A), t.selfEnd(A), t.cross(R, A), t.cross(R, C), t.cross(A, C)
    ]));
    tailed = { as_ref: { hairpin_th: r[0], self_any_th: r[1], self_end_th: r[2] }, as_alt: { hairpin_th: r[3], self_any_th: r[4], self_end_th: r[5] } };
    tailedPairs = {
      ref_alt_any_th: r[6].any, ref_alt_end_th: r[6].end,
      ref_common_any_th: r[7].any, ref_common_end_th: r[7].end,
      alt_common_any_th: r[8].any, alt_common_end_th: r[8].end
    };
  }

  const disc = sets.discrimination({
    orientation: o, variant: v, template: vt,
    as_ref: { target_seq: refTarget, matched_seq: as.seq, start: as.start, end: as.end },
    as_alt: { target_seq: altTarget, matched_seq: alt.matched_seq, start: alt.start, end: alt.end }
  });
  const hits = sets.primerNeighbours({ orientation: o, variant: v, template: vt, as: as, alt: alt, common: common, neighbours: ctx.neighbours });
  const balance = sets.tmBalance({
    as_ref_tm: mm ? mm.ref.duplex : as.tm, as_alt_tm: mm ? mm.alt.duplex : altOut.tm,
    as_ref_nominal: as.tm, as_alt_nominal: altOut.tm, common_tm: common.tm
  }, g);
  const issues = sets.issuesFor({
    g: g, variant: v, orientation: o,
    problems: { as_ref: mm ? mm.ref.out.problems : null, as_alt: mm ? mm.alt.out.problems : altOut.problems },
    balance: balance, dyes: dyes, tailed: tailed, tailed_pairs: tailedPairs,
    mismatch_not_applicable: notApplicable,
    mismatch_structures: mm ? { as_ref: structuresOf(mm.ref.out), as_alt: structuresOf(mm.alt.out) } : null,
    neighbours: hits, discrimination: disc
  });
  const score = sets.scoreSet({
    penalty: candidate.pair.penalty, relaxation_level: candidate.level, balance: balance, issues: issues,
    kasp: assay.type === 'kasp', product_size: candidate.pair.product_size, g: g
  });

  const gpos = function (t) { return vt.start + t - 1; };
  const span = function (s, e) { return { region: v.region, start: gpos(s), end: gpos(e) }; };
  const refSide = mm ? mm.ref.out : as;
  const altSideOut = mm ? mm.alt.out : altOut;
  const refGenomic = Object.assign(span(as.start, as.end), { strand: forward ? 1 : -1, blocks: [{ start: gpos(as.start), end: gpos(as.end) }] });
  const commonGenomic = Object.assign(span(common.start, common.end), { strand: forward ? -1 : 1, blocks: [{ start: gpos(common.start), end: gpos(common.end) }] });
  const primers = {
    as_ref: sets.oligo({
      role: 'as_ref', allele: v.vcf.ref, target_seq: refTarget, matched_seq: as.seq, dye: dyes.as_ref, tail_seq: tailOf('as_ref'),
      tm: mm ? mm.ref.duplex : as.tm, tm_method: mm ? 'ntthal_duplex' : 'primer3', matched_tm: mm ? as.tm : null,
      hairpin_th: refSide.hairpin_th, self_any_th: refSide.self_any_th, self_end_th: refSide.self_end_th, end_stability: refSide.end_stability,
      primer3_problems: mm ? mm.ref.out.problems : null,
      template: { start: as.start, end: as.end, sequence: 'ref' }, genomic: refGenomic, inserted_bases: 0,
      deliberate_mismatch: mm ? mm.ref.block : null, discrimination: disc.as_ref, tailed: tailed && tailed.as_ref, neighbours: hits.as_ref
    }),
    as_alt: sets.oligo({
      role: 'as_alt', allele: v.vcf.alt, target_seq: altTarget, matched_seq: alt.matched_seq, dye: dyes.as_alt, tail_seq: tailOf('as_alt'),
      tm: mm ? mm.alt.duplex : altOut.tm, tm_method: mm ? 'ntthal_duplex' : 'primer3', matched_tm: mm ? altOut.tm : null,
      hairpin_th: altSideOut.hairpin_th, self_any_th: altSideOut.self_any_th, self_end_th: altSideOut.self_end_th, end_stability: altSideOut.end_stability,
      primer3_problems: altSideOut.problems,
      template: { start: alt.start, end: alt.end, sequence: 'alt' }, genomic: alt.genomic, inserted_bases: alt.inserted_bases,
      deliberate_mismatch: mm ? mm.alt.block : null, discrimination: disc.as_alt, tailed: tailed && tailed.as_alt, neighbours: hits.as_alt
    }),
    common: sets.oligo({
      role: 'common', allele: null, target_seq: common.seq, matched_seq: common.seq, dye: null, tail_seq: null,
      tm: common.tm, tm_method: 'primer3', matched_tm: null,
      hairpin_th: common.hairpin_th, self_any_th: common.self_any_th, self_end_th: common.self_end_th, end_stability: common.end_stability,
      primer3_problems: null,
      template: { start: common.start, end: common.end, sequence: 'ref' }, genomic: commonGenomic, inserted_bases: 0,
      deliberate_mismatch: null, discrimination: null, tailed: null, neighbours: hits.common
    })
  };
  const thermo = {
    // The ordered pairs: the mismatch runs when a deliberate mismatch is applied (§4.10), else the design and ALT runs.
    ref_common: pairThermo(mm ? mm.ref.run.pair : candidate.pair),
    alt_common: pairThermo(mm ? mm.alt.run.pair : altRun.pair),
    tailed: tailedPairs ? Object.keys(tailedPairs).reduce(function (acc, k) { acc[k] = rounded(tailedPairs[k]); return acc; }, {}) : null
  };

  candidate.scored = {
    key: sets.setKey(o, refTarget, altTarget, common.seq),
    orientation: o,
    relaxation_level: candidate.level,
    pair_index: candidate.pair_index,
    exact: score.exact,
    score: score.score,
    quality: score.quality,
    counted: score.counted,
    terms: score.terms,
    primers: primers,
    products: sets.productSpans(candidate, vt, v, altRun.pair.product_size),
    thermo: thermo,
    tm_balance: { as_tm_diff: decimal.round(balance.as_tm_diff, 2), common_minus_as: decimal.round(balance.common_minus_as, 2) },
    neighbour_sites: sets.neighbourSites(hits),
    primer3_penalty: decimal.roundValue(candidate.pair.penalty, 4),
    issues: issues
  };
  return { set: candidate };
}

module.exports = {
  scorePair,
  scoreCandidate,
  scoringParams,
  scoringRecord,
  belowFloors,
  TAILS,
  DYES,
  WIDENED
};
