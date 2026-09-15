'use strict';

// POST /primers/genotyping/design (spec §4.1-§4.8).
//
// runOrientation(orientation, ctx) is step 9 of §4.1: the neighbour blocker check, then one Primer3 design run per
// relaxation level, with the allele-specific 3' end forced and the SEQUENCE_TARGET guard beside the zone (§4.4, §4.7).
// Each returned pair goes, in Primer3 rank order, through the §4.8 filters: force guard, zone guard on both
// haplotypes, common-primer 3' neighbour, duplicates, the per-orientation scoring cap and the budget reservation. A
// pair that passes is floored on its designed allele-specific Tm and handed to the injected scorer. Every attempt
// keeps Primer3's explain data, and pairs_returned always equals rejected + not_scored + sets. The ladder stops at the
// first level that keeps a set, or when the budget stops it.
//
// designGenotyping(body, deps) orchestrates steps 1 and 6-9: request, deadline, the injected variant resolution (steps
// 2-5, before the semaphore), semaphore, template, repeat mask and both orientations. Sets (steps 10-16: ALT
// derivation, scoring, ranking, order rows) are built on the candidates runOrientation returns and are not produced here.

const boulder = require('../boulder');
const design = require('../design');
const { PrimerHttpError } = require('../errors');
const guard = require('./guard');
const request = require('./request');
const templates = require('./template');

const ORIENTATIONS = Object.freeze(['forward', 'reverse']);
const GENOTYPING_DESIGN_VERSION = '1';
const REJECTED_KEYS = Object.freeze(['force', 'overlap', 'common_neighbour_3p', 'alt_scoring_failed', 'below_floor', 'duplicate']);
const REQUIRED_CONFIG = Object.freeze(['template_flank', 'num_return_per_run', 'max_scored_per_orientation', 'max_primer3_runs',
  'max_thermo_calls', 'guard_gap', 'mask_exempt_pad', 'neighbour_3p_window']);

// Worst-case cost reserved for one candidate before it is scored (§4.8 step 6): one check_primers run, or three with a
// deliberate mismatch; 15 ntthal calls with tails, plus 2 mismatch duplexes.
const RESERVE = Object.freeze({ runs: 1, mismatch_runs: 3, tailed_calls: 15, mismatch_calls: 2 });

function logSafe(log, level, msg) {
  try {
    const l = log || console;
    (l[level] || l.error || l.log).call(l, msg);
  } catch (e) { /* logging must never throw */ }
}

function warning(code, message, details) {
  return { code: code, message: message, details: details || {} };
}

function deadlineError(ms) {
  return new PrimerHttpError(504, 'DEADLINE_EXCEEDED', 'the design did not finish within ' + ms + ' ms', { deadline_ms: ms });
}

// cfg.genotyping with every key this module reads (§3.1).
function genotypingConfig(cfg) {
  const g = cfg && cfg.genotyping;
  REQUIRED_CONFIG.forEach(function (k) {
    if (!g || !Number.isSafeInteger(g[k]) || g[k] < 0) throw new TypeError('primers config: genotyping.' + k + ' must be a non-negative integer');
  });
  return g;
}

// ---- exact decimals ------------------------------------------------------------------------------------

function decimalParts(x) {
  const s = typeof x === 'number' ? (Number.isFinite(x) ? String(x) : '') : String(x).trim();
  const m = /^(-?)(\d+)(?:\.(\d+))?$/.exec(s);
  return m ? { sign: m[1], int: m[2], frac: m[3] || '' } : null;
}

// Exact comparison of two finite decimals in plain notation (a value as Primer3 prints it, a config threshold), as the
// §2.1 rule requires of threshold tests. -> -1, 0 or 1
function compareDecimal(a, b) {
  const pa = decimalParts(a);
  const pb = decimalParts(b);
  if (!pa || !pb) throw new TypeError('not a plain decimal: ' + String(pa ? b : a));
  const places = Math.max(pa.frac.length, pb.frac.length);
  const va = BigInt(pa.sign + pa.int + pa.frac.padEnd(places, '0'));
  const vb = BigInt(pb.sign + pb.int + pb.frac.padEnd(places, '0'));
  return va < vb ? -1 : va > vb ? 1 : 0;
}

// ---- budget --------------------------------------------------------------------------------------------

// createBudget(limits, opts) -> the Primer3 run and ntthal call counters of one design request (§4.8 step 6, §4.16).
//   limits: {max_primer3_runs, max_thermo_calls} (primers.genotyping)
//   opts.thermoCalls(): the distinct ntthal calls made so far (thermo.js memoizes); without it, the countThermo() total
// Every Primer3 run is counted when it starts, design runs included; only candidates are checked against the caps.
// exhausted lists, in order, the orientations whose ladder the budget stopped.
function createBudget(limits, opts) {
  const l = limits || {};
  ['max_primer3_runs', 'max_thermo_calls'].forEach(function (k) {
    if (!Number.isSafeInteger(l[k]) || l[k] < 0) throw new TypeError(k + ' must be a non-negative integer');
  });
  const thermoCalls = opts && typeof opts.thermoCalls === 'function' ? opts.thermoCalls : null;
  let thermo = 0;
  return {
    max_primer3_runs: l.max_primer3_runs,
    max_thermo_calls: l.max_thermo_calls,
    primer3_runs: 0,
    get thermo_calls() { return thermoCalls ? thermoCalls() : thermo; },
    exhausted: [],
    countPrimer3: function (n) { this.primer3_runs += n === undefined ? 1 : n; },
    countThermo: function (n) { thermo += n === undefined ? 1 : n; },
    fits: function (runs, calls) {
      return this.primer3_runs + runs <= this.max_primer3_runs && this.thermo_calls + calls <= this.max_thermo_calls;
    }
  };
}

// The worst-case cost one candidate reserves for the effective assay. -> {primer3_runs, thermo_calls}
function reservation(assay) {
  const mismatch = assay.deliberate_mismatch === 'auto';
  return {
    primer3_runs: mismatch ? RESERVE.mismatch_runs : RESERVE.runs,
    thermo_calls: (assay.tails !== 'none' ? RESERVE.tailed_calls : 0) + (mismatch ? RESERVE.mismatch_calls : 0)
  };
}

function emptyRejected() {
  const r = {};
  REJECTED_KEYS.forEach(function (k) { r[k] = 0; });
  return r;
}

// ---- neighbours ----------------------------------------------------------------------------------------

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

// Known variants under reference coordinates listed from a primer's 3' end (§4.15), as PrimerNeighbourHit objects:
// distance_from_3p is 1 + the index of the first listed coordinate inside the neighbour's minimal span (an insertion
// covers both flanking bases). Sorted by distance, then key. neighbours: canonical entries; targetKey is skipped;
// coords may hold null (an inserted base).
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
  return hits.sort(function (x, y) {
    return x.distance_from_3p - y.distance_from_3p || (x.key < y.key ? -1 : x.key > y.key ? 1 : 0);
  });
}

function nonEms(hits) {
  return hits.filter(function (h) { return !h.ems; });
}

// ---- one pair ------------------------------------------------------------------------------------------

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

// ---- Primer3 -------------------------------------------------------------------------------------------

// The design-run record of §4.7: design.buildRecord plus SEQUENCE_FORCE_LEFT_END or SEQUENCE_FORCE_RIGHT_END, the
// SEQUENCE_TARGET guard and PRIMER_NUM_RETURN, appended in that order. No PRIMER_PICK_ANYWAY, no included or excluded
// region, no client tag.
//   vt: the variant template; o: {seq (the orientation's Primer3 template), params (effective level params),
//       avoid_repeats, repeat_mask_mode, guard_gap, num_return}
function designRecord(vt, orientation, level, o) {
  const target = guard.targetTag(vt.features.zone, orientation, vt, o.guard_gap);
  if (target === null) throw new RangeError('no template base is left for the ' + orientation + ' SEQUENCE_TARGET guard');
  const tags = design.buildRecord({ id: vt.id + '_' + orientation + '_L' + level, seq: o.seq }, {
    params: Object.assign({}, o.params, { num_return: o.num_return }),
    target: null,
    included: null,
    excluded: [],
    junctions: null,
    avoid_repeats: o.avoid_repeats,
    repeat_mask_mode: o.repeat_mask_mode
  });
  if (orientation === 'forward') tags.SEQUENCE_FORCE_LEFT_END = vt.features.discriminating.forward;
  else tags.SEQUENCE_FORCE_RIGHT_END = vt.features.discriminating.reverse;
  tags.SEQUENCE_TARGET = target;
  tags.PRIMER_NUM_RETURN = o.num_return;
  return tags;
}

// One design run under the deadline, counted in the budget. A PRIMER_ERROR of a design run is 400 PRIMER3_INPUT_ERROR.
async function runDesignRecord(ctx, tags) {
  const cfg = ctx.cfg;
  const dl = ctx.deadline;
  const remaining = Math.floor(dl.remaining());
  if (remaining <= 0) throw deadlineError(cfg.design.deadline_ms);
  ctx.budget.countPrimer3();
  const result = await dl.within(ctx.primer3.run(tags, {
    timeoutMs: Math.min(cfg.design.primer3_timeout_ms, remaining), signal: dl.signal, log: ctx.log, tmpDir: cfg.tmp_dir
  }));
  if (result.error) {
    throw new PrimerHttpError(400, 'PRIMER3_INPUT_ERROR', 'Primer3 rejected the input: ' + result.error, { primer3_error: result.error });
  }
  return result;
}

// The default scorer keeps every candidate as it is; the scoring modules replace it (ctx.scoreCandidate).
async function keepCandidate(candidate) {
  return { set: candidate };
}

// The hard Tm floor of the designed allele-specific (REF) primer (§4.6), exact. The derived primers are floored by
// the scorer.
function belowFloor(oligo, floors) {
  return oligo.tm === null || oligo.tm === undefined || compareDecimal(oligo.tm, floors.as_min_tm) < 0;
}

// ---- one orientation -----------------------------------------------------------------------------------

// runOrientation(orientation, ctx) -> Promise<{orientation, sets, budget_exhausted, warnings}>
//   ctx: {cfg, req (request.normalize), variant (canonical entry), template (buildVariantTemplate), maskedSeq (the
//         repeat-masked template, or null), neighbours (canonical entries of the template window), budget
//         (createBudget), deadline (design.createDeadline), primer3 {run}, scoreCandidate(candidate, ctx) ->
//         {set} | {dropped: 'alt_scoring_failed' | 'below_floor'} (default keepCandidate), log}
//   orientation: PrimerGenotypingOrientation {status, reason, discriminating_position, relaxation_level, sets_found,
//                blockers, attempts[{level, changes, explain {left, right, pair}, pairs_returned, rejected, not_scored,
//                sets}]}
//   sets: what the scorer kept, in Primer3 order (each candidate also carries level and params)
//   budget_exhausted: the budget stopped this orientation's ladder (also pushed to budget.exhausted)
//   warnings: PRIMER3_WARNING and VARIANT_IN_REPEAT of this orientation
async function runOrientation(orientation, ctx) {
  if (ORIENTATIONS.indexOf(orientation) < 0) throw new TypeError('orientation must be forward or reverse');
  const req = ctx.req;
  const g = genotypingConfig(ctx.cfg);
  const vt = ctx.template;
  const v = ctx.variant;
  const d = v.discriminating[orientation];
  const pub = {
    status: 'ok', reason: null, discriminating_position: d.position, relaxation_level: null, sets_found: 0, blockers: [], attempts: []
  };
  const out = { orientation: pub, sets: [], budget_exhausted: false, warnings: [] };
  const skip = req.orientations.indexOf(orientation) < 0 ? 'not_requested' : vt.skip[orientation];
  if (skip) {
    pub.status = 'skipped';
    pub.reason = skip;
    return out;
  }

  // §4.15: a non-EMS neighbour in the allele-specific primer's last bases blocks a natural target's orientation.
  const block = req.assay.neighbour_policy === 'avoid_3p' && v.ems !== true;
  const at3p = nonEms(neighbourHits(v.key, ctx.neighbours, coordsFrom3p(d.position, g.neighbour_3p_window, orientation === 'forward' ? 1 : -1)));
  if (block && at3p.length > 0) {
    pub.status = 'blocked';
    pub.reason = 'neighbour_at_3p';
    pub.blockers = at3p;
    return out;
  }

  const p3 = templates.orientationTemplate(vt, ctx.maskedSeq, orientation);
  if (p3.masked_bases > 0) {
    out.warnings.push(warning('VARIANT_IN_REPEAT', orientation + ': ' + p3.masked_bases +
      ' bases of the allele-specific primer window are masked as repeat and were exempted from the mask',
    { orientation: orientation, masked_bases: p3.masked_bases }));
  }
  const screen = {
    orientation: orientation, as_3p: vt.features.discriminating[orientation], zone: vt.features.zone, delta: vt.features.alt_offset,
    template_start: vt.start, target_key: v.key, neighbours: ctx.neighbours, window: g.neighbour_3p_window, block: block
  };
  const reserve = reservation(req.assay);
  const score = ctx.scoreCandidate || keepCandidate;
  const seen = new Set();
  let scored = 0;

  for (let level = 0; level <= req.assay.max_relaxation; level++) {
    const params = req.levels[level];
    const tags = designRecord(vt, orientation, level, {
      seq: p3.seq, params: params, avoid_repeats: req.avoid_repeats, repeat_mask_mode: req.repeat_mask_mode,
      guard_gap: g.guard_gap, num_return: g.num_return_per_run
    });
    const result = await runDesignRecord(ctx, tags);
    if (result.warning) out.warnings.push(warning('PRIMER3_WARNING', result.warning, { orientation: orientation, level: level }));
    const pairs = boulder.extractPairs(result.tags);
    const attempt = {
      level: level,
      changes: level === 0 ? {} : req.ladder[level - 1].changes,
      explain: boulder.extractExplain(result.tags),
      pairs_returned: pairs.length,
      rejected: emptyRejected(),
      not_scored: 0,
      sets: 0
    };
    pub.attempts.push(attempt);

    let stopped = false;
    for (let i = 0; i < pairs.length; i++) {
      const screened = screenPair(pairs[i], screen);
      if (screened.reject) {
        attempt.rejected[screened.reject]++;
        continue;
      }
      const candidate = screened.candidate;
      const dup = candidate.as.seq + '|' + candidate.common.seq;
      if (seen.has(dup)) {
        attempt.rejected.duplicate++;
        continue;
      }
      seen.add(dup);
      if (scored >= g.max_scored_per_orientation) {
        attempt.not_scored++;
        continue;
      }
      if (!ctx.budget.fits(reserve.primer3_runs, reserve.thermo_calls)) {
        // This pair and every later pair of the level stay unscored, and the ladder stops (§4.8 step 6).
        attempt.not_scored += pairs.length - i;
        stopped = true;
        break;
      }
      scored++;
      candidate.level = level;
      candidate.params = params;
      const s = belowFloor(candidate.as, req.floors) ? { dropped: 'below_floor' } : await score(candidate, ctx);
      if (s.dropped) {
        if (REJECTED_KEYS.indexOf(s.dropped) < 0) throw new TypeError('scoreCandidate dropped a candidate for an unknown reason: ' + s.dropped);
        attempt.rejected[s.dropped]++;
        continue;
      }
      attempt.sets++;
      out.sets.push(s.set);
    }
    if (stopped) {
      out.budget_exhausted = true;
      ctx.budget.exhausted.push(orientation);
    }
    if (out.sets.length > 0) {
      pub.relaxation_level = level;
      break;
    }
    if (stopped) break;
  }

  pub.sets_found = out.sets.length;
  if (out.sets.length === 0) {
    pub.status = 'no_sets';
    if (out.budget_exhausted) pub.reason = 'budget_exhausted';
    else if (pub.attempts.some(function (a) { return a.rejected.below_floor > 0; })) pub.reason = 'below_floor';
  }
  return out;
}

// Response warnings about the orientations (§2.15): ORIENTATION_BLOCKED, ORIENTATION_SKIPPED (not for an orientation
// the client did not request), ORIENTATION_NO_SETS, RELAXED_CONSTRAINTS, and NO_SETS when neither orientation kept a set.
function orientationWarnings(results) {
  const out = [];
  ORIENTATIONS.forEach(function (o) {
    const r = results[o].orientation;
    if (r.status === 'blocked') {
      const ids = r.blockers.map(function (b) { return b.ids[0] || b.key; });
      const distances = r.blockers.map(function (b) { return b.distance_from_3p; });
      out.push(warning('ORIENTATION_BLOCKED', o + ': ' + r.blockers.map(function (b, i) {
        return 'known variant ' + ids[i] + ' (' + b.alleles + ') lies ' + b.distance_from_3p + ' nt from the allele-specific primer\'s 3′ end';
      }).join('; '), { orientation: o, ids: ids, distances: distances }));
    } else if (r.status === 'skipped' && r.reason !== 'not_requested') {
      out.push(warning('ORIENTATION_SKIPPED', o + ': ' + (r.reason === 'too_close_to_end'
        ? 'too close to the end of the sequence' : 'N bases in the allele-specific primer window'), { orientation: o, reason: r.reason }));
    } else if (r.status === 'no_sets') {
      out.push(warning('ORIENTATION_NO_SETS', o + ': no primer set within ' + r.attempts.length + ' relaxation level' +
        (r.attempts.length === 1 ? '' : 's'), { orientation: o, levels_tried: r.attempts.length }));
    } else if (r.status === 'ok' && r.relaxation_level >= 1) {
      out.push(warning('RELAXED_CONSTRAINTS', o + ': sets need relaxation level ' + r.relaxation_level,
        { orientation: o, level: r.relaxation_level, changes: r.attempts[r.relaxation_level].changes }));
    }
  });
  if (ORIENTATIONS.every(function (o) { return results[o].orientation.sets_found === 0; })) {
    out.push(warning('NO_SETS', 'no allele-specific primer set was found in either orientation', {}));
  }
  return out;
}

function versionOrNull(primer3, log) {
  return Promise.resolve().then(function () { return primer3.version(); }).catch(function (err) {
    logSafe(log, 'warn', 'primers: primer3 version unavailable: ' + (err && err.code));
    return null;
  });
}

// ---- orchestration -------------------------------------------------------------------------------------

// designGenotyping(body, deps) -> Promise<response>
//   deps: {cfg, log, signal (client abort), now, semaphore {acquire}, primer3 {run, version}, repeatMask {repeatMask},
//          sequence, scoreCandidate, and
//          resolveVariant(req, {signal, deadline}) -> {assembly (assemblies.resolve), variant (canonical entry, REF
//            verified), neighbours {data, entries}, warnings [], variation_source}: the variant and neighbour
//            resolution of §4.1 steps 2-5, which runs before the semaphore}
// response: {variant, template, assay, orientations (null for template_only), sets [], check null, settings, engine,
//            warnings, candidates {forward, reverse} (the kept candidates; not part of the public response)}
// template_only skips the semaphore unless avoid_repeats needs the masker (§4.2).
async function designGenotyping(body, deps) {
  deps = deps || {};
  const cfg = deps.cfg || require('../config').get();
  if (cfg.enabled === false) throw new PrimerHttpError(503, 'FEATURE_DISABLED', 'primer design is disabled on this server', {});
  const log = deps.log || console;
  const req = request.normalize(body, cfg);
  genotypingConfig(cfg);
  const dl = design.createDeadline(cfg.design.deadline_ms, deps.signal, deps.now);
  let release = null;
  try {
    if (typeof deps.resolveVariant !== 'function') {
      throw new PrimerHttpError(500, 'INTERNAL', 'genotyping variant resolution is not configured', {});
    }
    const resolved = await dl.within(deps.resolveVariant(req, { signal: dl.signal, deadline: dl }));
    if (!req.template_only || req.avoid_repeats) {
      const semaphore = deps.semaphore || require('../semaphore').designSemaphore();
      release = await semaphore.acquire({ signal: dl.signal });
    }
    const tdeps = Object.assign({}, deps, { cfg: cfg, req: req, signal: dl.signal, log: log });
    const vt = await dl.within(templates.buildVariantTemplate(resolved.variant, resolved.assembly, tdeps));
    const warnings = (resolved.warnings || []).concat(vt.warnings);

    let mask = { masked: false, mask_source: null, mask: [], masked_fraction: 0 };
    let maskedSeq = null;
    if (req.avoid_repeats) {
      const masker = deps.repeatMask || require('../repeat_mask');
      const m = await dl.within(masker.repeatMask(vt.region_template, { mode: req.repeat_mask_mode, signal: dl.signal, deadline: dl.deadlineAt }, tdeps));
      mask = { masked: m.masked, mask_source: m.mask_source, mask: m.mask, masked_fraction: m.masked_fraction };
      maskedSeq = m.seq;
      m.warnings.forEach(function (w) { warnings.push(w); });
    }

    const primer3 = deps.primer3 || require('../primer3');
    const response = {
      variant: resolved.variant,
      template: templates.publicVariantTemplate(vt, mask),
      assay: Object.assign({}, req.assay, { ems_target: resolved.variant.ems === true }),
      orientations: null,
      sets: [],
      check: null,
      settings: { preset: req.preset, params: req.params, pinned: req.pinned, ladder: req.ladder, floors: req.floors },
      engine: { primer3: null, thermo: null, genotyping_design: GENOTYPING_DESIGN_VERSION, variation_source: resolved.variation_source || null },
      warnings: warnings
    };
    if (req.template_only) return response;

    response.engine.primer3 = await dl.within(versionOrNull(primer3, log));
    const ctx = {
      cfg: cfg,
      req: req,
      variant: resolved.variant,
      template: vt,
      maskedSeq: maskedSeq,
      neighbours: resolved.neighbours && Array.isArray(resolved.neighbours.entries) ? resolved.neighbours.entries : [],
      budget: deps.budget || createBudget(cfg.genotyping),
      deadline: dl,
      primer3: primer3,
      scoreCandidate: deps.scoreCandidate,
      log: log
    };
    const results = {};
    for (const o of ORIENTATIONS) {
      results[o] = await runOrientation(o, ctx);
      results[o].warnings.forEach(function (w) { warnings.push(w); });
    }
    response.orientations = { forward: results.forward.orientation, reverse: results.reverse.orientation };
    orientationWarnings(results).forEach(function (w) { warnings.push(w); });
    Object.defineProperty(response, 'candidates', {
      value: { forward: results.forward.sets, reverse: results.reverse.sets }, enumerable: false
    });
    Object.defineProperty(response, 'budget', { value: ctx.budget, enumerable: false });
    return response;
  } finally {
    dl.dispose();
    if (release) release();
  }
}

module.exports = {
  designGenotyping,
  runOrientation,
  screenPair,
  designRecord,
  neighbourHits,
  coordsFrom3p,
  createBudget,
  reservation,
  keepCandidate,
  compareDecimal,
  orientationWarnings,
  emptyRejected,
  ORIENTATIONS,
  REJECTED_KEYS,
  GENOTYPING_DESIGN_VERSION
};
