'use strict';

require('../../../api/helpers/primers/node_compat');

const test = require('node:test');
const should = require('should');
const fs = require('fs');

const gd = require('../../../api/helpers/primers/genotyping/design');
const gtemplate = require('../../../api/helpers/primers/genotyping/template');
const scoring = require('../../../api/helpers/primers/genotyping/scoring');
const request = require('../../../api/helpers/primers/genotyping/request');
const design = require('../../../api/helpers/primers/design');
const template = require('../../../api/helpers/primers/template');
const boulder = require('../../../api/helpers/primers/boulder');
const normalize = require('../../../api/helpers/primers/variation/normalize');
const repeatMask = require('../../../api/helpers/primers/repeat_mask');
const { revcomp } = require('../../../api/helpers/primers/sequence');
const stubs = require('../fixtures/design/stubs');
const fx = require('../fixtures/catalog/fsfixture');
const cases = require('../fixtures/primer3/genotyping/cases');

// Spec §4.1-§4.18 and the genotyping_design rows of §7.3: records, explain strings, pairs_returned, not_scored, the
// scoring cap and budget reservation, the guard, floors and mask exemption (M5); the complete §2.9 response, the given
// values of the §2.10 excerpts, the budget warning, template_only, Ensembl outages and the deadline (M6). No
// primer3_core or ntthal is spawned outside the PRIMERS_REALDATA=1 tests.
// Recorded fixtures: test/primers/fixtures/primer3/genotyping/ and test/primers/fixtures/thermo/genotyping.json
// (record_genotyping.js, primer3_core and ntthal 2.6.1, 2026-09-15). Apart from SEQUENCE_ID, every recorded Primer3 input
// and output is byte-identical to the design fix pass recordings of the spec (genotyping-design/scratch-fix2/score2/
// default/p3), which the §2.9-§2.10 values were computed from. expected_<case>.json hold the §2.9 response and the
// §2.10 excerpts exactly as the spec prints them ("…" marks an elided value).

const REALDATA = process.env.PRIMERS_REALDATA === '1';
const PRIMER3_BIN = '/home/olson/bin/primer3_core';
const NTTHAL_BIN = '/home/olson/primer3-2.6.1/bin/ntthal';

function expected(name) {
  return JSON.parse(fs.readFileSync(cases.DIR + '/expected_' + name + '.json', 'utf8'));
}

function rejected(over) {
  return Object.assign({ force: 0, overlap: 0, common_neighbour_3p: 0, alt_scoring_failed: 0, below_floor: 0, duplicate: 0 }, over || {});
}

function counts(a) {
  return { level: a.level, pairs_returned: a.pairs_returned, rejected: a.rejected, not_scored: a.not_scored, sets: a.sets };
}

function raws(a) {
  return { left: a.explain.left.raw, right: a.explain.right.raw, pair: a.explain.pair.raw };
}

async function rejectsWith(promise, status, code) {
  let err = null;
  try {
    await promise;
  } catch (e) {
    err = e;
  }
  should.exist(err, 'expected ' + code);
  err.should.have.property('code', code);
  err.should.have.property('status', status);
  return err;
}

async function replay(name, extra, body) {
  const calls = [];
  const thermoCalls = [];
  const x = Object.assign({ primer3: cases.recordedPrimer3(calls) }, extra || {});
  if (x.thermo === undefined) x.thermo = cases.recordedThermo(x.cfg || cases.cfg(), thermoCalls);
  const deps = cases.deps(name, x);
  const res = await gd.designGenotyping(body || cases.body(name), deps);
  return { res: res, calls: calls, thermoCalls: thermoCalls, deps: deps };
}

// A scorer that keeps every candidate unscored and counts `runs` check_primers runs for it, as scoring would.
function countingScorer(runs) {
  return async function (candidate, ctx) {
    ctx.budget.countPrimer3(runs);
    return { set: candidate };
  };
}

// Every value the spec gives in an excerpt equals the response ("…" keys and values are skipped). -> differences
function subsetDiff(got, want, p, out) {
  out = out || [];
  p = p || '$';
  if (want === '…') return out;
  if (Array.isArray(want)) {
    if (!Array.isArray(got)) out.push(p + ': not an array');
    else {
      if (got.length !== want.length) out.push(p + ': length ' + got.length + ', want ' + want.length);
      want.forEach(function (w, i) { subsetDiff(got[i], w, p + '[' + i + ']', out); });
    }
    return out;
  }
  if (want !== null && typeof want === 'object') {
    if (got === null || typeof got !== 'object') out.push(p + ': got ' + JSON.stringify(got));
    else {
      Object.keys(want).forEach(function (k) {
        if (k === '…') return;
        if (!(k in got)) out.push(p + '.' + k + ': missing');
        else subsetDiff(got[k], want[k], p + '.' + k, out);
      });
    }
    return out;
  }
  if (got !== want) out.push(p + ': got ' + JSON.stringify(got) + ', want ' + JSON.stringify(want));
  return out;
}

// The candidate of design-run pair `index`; e: {level, as [seq, start, end], common [seq, start, end] (genomic), size,
// penalty, as_tm, common_tm}.
function expectCandidate(res, orientation, index, e) {
  const c = res.candidates[orientation].find(function (x) { return x.pair_index === index; });
  should.exist(c, orientation + ' pair ' + index);
  const t = function (g) { return g - res.template.start + 1; };
  c.orientation.should.equal(orientation);
  c.level.should.equal(e.level);
  c.as.seq.should.equal(e.as[0]);
  [c.as.start, c.as.end].should.eql([t(e.as[1]), t(e.as[2])]);
  c.common.seq.should.equal(e.common[0]);
  [c.common.start, c.common.end].should.eql([t(e.common[1]), t(e.common[2])]);
  c.pair.product_size.should.equal(e.size);
  c.pair.penalty.should.equal(e.penalty);
  // The spec quotes Tm to 2 decimals, rounded half away from zero (60.425 -> 60.43).
  c.as.tm.should.be.approximately(e.as_tm, 0.0051);
  c.common.tm.should.be.approximately(e.common_tm, 0.0051);
  return c;
}

// Primer3-like output on a template: pairs [{left: [pos, len, tm], right: [pos (5' end), len, tm]}].
function synthetic(seq, pairs, warning) {
  const tags = {
    PRIMER_LEFT_EXPLAIN: 'considered 13, ok 2', PRIMER_RIGHT_EXPLAIN: 'considered 100, ok 50', PRIMER_PAIR_EXPLAIN: 'considered 20, ok ' + pairs.length,
    PRIMER_FIRST_BASE_INDEX: '1', PRIMER_PAIR_NUM_RETURNED: String(pairs.length)
  };
  pairs.forEach(function (p, i) {
    tags['PRIMER_PAIR_' + i + '_PENALTY'] = '1.500000';
    tags['PRIMER_LEFT_' + i] = p.left[0] + ',' + p.left[1];
    tags['PRIMER_RIGHT_' + i] = p.right[0] + ',' + p.right[1];
    tags['PRIMER_LEFT_' + i + '_SEQUENCE'] = seq.slice(p.left[0] - 1, p.left[0] - 1 + p.left[1]);
    tags['PRIMER_RIGHT_' + i + '_SEQUENCE'] = revcomp(seq.slice(p.right[0] - p.right[1], p.right[0]));
    tags['PRIMER_LEFT_' + i + '_TM'] = p.left[2] || '58.000';
    tags['PRIMER_RIGHT_' + i + '_TM'] = p.right[2] || '58.000';
    tags['PRIMER_PAIR_' + i + '_PRODUCT_SIZE'] = String(p.right[0] - p.left[0] + 1);
  });
  if (warning) tags.PRIMER_WARNING = warning;
  return { tags: tags, error: null, warning: warning || null, exitCode: 0 };
}

function primer3Stub(result, calls) {
  return {
    run: async function (tags) {
      if (calls) calls.push(tags);
      return typeof result === 'function' ? result(tags) : result;
    },
    version: async function () { return '2.6.1'; }
  };
}

// A thermo that must never be used.
const NO_THERMO = Object.freeze({
  calls: 0,
  run: function () { throw new Error('ntthal must not run'); },
  hairpin: function () { throw new Error('ntthal must not run'); },
  selfAny: function () { throw new Error('ntthal must not run'); },
  selfEnd: function () { throw new Error('ntthal must not run'); },
  cross: function () { throw new Error('ntthal must not run'); },
  duplex: function () { throw new Error('ntthal must not run'); }
});

// ctx for runOrientation on one case; over: {body (merged into the case body), ctx (merged into ctx)}. Call dispose().
async function contextFor(name, over) {
  over = over || {};
  const deps = cases.deps(name);
  const req = request.normalize(Object.assign(cases.body(name), over.body || {}), deps.cfg);
  const resolved = await cases.resolve(req, deps);
  const vt = await gtemplate.buildVariantTemplate(resolved.variant, resolved.assembly, Object.assign({}, deps, { req: req }));
  const dl = design.createDeadline(45000);
  const ctx = Object.assign({
    cfg: deps.cfg, req: req, variant: resolved.variant, template: vt, maskedSeq: null, neighbours: resolved.neighbours.entries,
    budget: gd.createBudget(deps.cfg.genotyping), deadline: dl, thermo: deps.thermo, log: stubs.silentLog
  }, over.ctx || {});
  ctx.dispose = function () { dl.dispose(); };
  return ctx;
}

const BUDGETS = Object.freeze({
  rs871475760_kasp: [18, 123], tmp_1_11193_C_T_kasp: [20, 150], rs5413864115_kasp: [20, 141], tmp_1_11502_C_CGT_kasp: [9, 66], rs871475760_as_pcr: [29, 4]
});

// ---- the recordings --------------------------------------------------------------------------------------------

test('recorded fixtures: every example regenerates its Primer3 records and ntthal calls byte for byte; §4.8 budgets 18/123 … 29/4', async function () {
  const index = cases.loadIndex();
  index.primer3_version.should.equal('2.6.1');
  Object.keys(index.runs).should.have.length(96);
  for (const name of Object.keys(cases.CASES)) {
    const r = await replay(name);
    r.calls.map(function (c) { return c.name; }).should.eql(index.cases[name]);
    r.calls.forEach(function (c) {
      c.input.should.equal(cases.readRecording(c.name + '.input.txt'));
      index.runs[cases.sha256(c.input)].name.should.equal(c.name);
    });
    [r.res.budget.primer3_runs, r.res.budget.thermo_calls].should.eql(BUDGETS[name], name);
    r.thermoCalls.should.have.length(BUDGETS[name][1]);
    r.res.warnings.map(function (w) { return w.code; }).should.not.containEql('DESIGN_BUDGET_EXHAUSTED');
  }
  // Any other input is refused.
  const stub = cases.recordedPrimer3();
  const tags = boulder.parse(cases.readRecording('rs871475760_kasp_forward_L0.input.txt'));
  tags.PRIMER_NUM_RETURN = '10';
  let refused = null;
  try { await stub.run(tags); } catch (e) { refused = e; }
  should.exist(refused);
  refused.code.should.equal('UNRECORDED_PRIMER3_INPUT');
  refused = null;
  try { await cases.recordedThermo().hairpin('ACGTACGTACGTACGTACGT'); } catch (e) { refused = e; }
  should.exist(refused);
  refused.code.should.equal('UNRECORDED_THERMO_INPUT');
});

test('§4.7 design records: buildRecord, then the forced 3′ end, the SEQUENCE_TARGET guard and PRIMER_NUM_RETURN 20; never PICK_ANYWAY', function () {
  const want = {
    rs871475760_kasp_forward_L0: ['SEQUENCE_FORCE_LEFT_END', '401', '402,10', 'sorghum_bicolor_1_10709-11509_forward_L0', '61-120', 801],
    rs871475760_kasp_reverse_L0: ['SEQUENCE_FORCE_RIGHT_END', '401', '391,10', 'sorghum_bicolor_1_10709-11509_reverse_L0', '61-120', 801],
    tmp_1_11193_C_T_kasp_forward_L0: ['SEQUENCE_FORCE_LEFT_END', '401', '402,10', 'sorghum_bicolor_1_10793-11593_forward_L0', '61-120', 801],
    tmp_1_11193_C_T_kasp_forward_L1: ['SEQUENCE_FORCE_LEFT_END', '401', '402,10', 'sorghum_bicolor_1_10793-11593_forward_L1', '65-150', 801],
    tmp_1_11193_C_T_kasp_reverse_L0: ['SEQUENCE_FORCE_RIGHT_END', '401', '391,10', 'sorghum_bicolor_1_10793-11593_reverse_L0', '61-120', 801],
    tmp_1_11193_C_T_kasp_reverse_L1: ['SEQUENCE_FORCE_RIGHT_END', '401', '391,10', 'sorghum_bicolor_1_10793-11593_reverse_L1', '65-150', 801],
    rs5413864115_kasp_forward_L0: ['SEQUENCE_FORCE_LEFT_END', '403', '405,10', 'sorghum_bicolor_1_10883-11685_forward_L0', '61-120', 803],
    rs5413864115_kasp_forward_L1: ['SEQUENCE_FORCE_LEFT_END', '403', '405,10', 'sorghum_bicolor_1_10883-11685_forward_L1', '65-150', 803],
    rs5413864115_kasp_reverse_L0: ['SEQUENCE_FORCE_RIGHT_END', '401', '390,10', 'sorghum_bicolor_1_10883-11685_reverse_L0', '61-120', 803],
    rs5413864115_kasp_reverse_L1: ['SEQUENCE_FORCE_RIGHT_END', '401', '390,10', 'sorghum_bicolor_1_10883-11685_reverse_L1', '65-150', 803],
    tmp_1_11502_C_CGT_kasp_reverse_L0: ['SEQUENCE_FORCE_RIGHT_END', '401', '391,10', 'sorghum_bicolor_1_11102-11903_reverse_L0', '61-120', 802],
    rs871475760_as_pcr_forward_L0: ['SEQUENCE_FORCE_LEFT_END', '401', '402,10', 'sorghum_bicolor_1_10709-11509_forward_L0', '150-300', 801],
    rs871475760_as_pcr_reverse_L0: ['SEQUENCE_FORCE_RIGHT_END', '401', '391,10', 'sorghum_bicolor_1_10709-11509_reverse_L0', '150-300', 801]
  };
  const index = cases.loadIndex();
  const runs = Object.keys(index.runs).map(function (h) { return index.runs[h]; });
  const designRuns = runs.filter(function (r) { return /_L\d$/.test(r.name); });
  designRuns.map(function (r) { return r.name; }).sort().should.eql(Object.keys(want).sort());
  designRuns.forEach(function (run) {
    const w = want[run.name];
    const text = cases.readRecording(run.input);
    const tags = boulder.parse(text);
    const keys = text.split('\n').filter(Boolean).map(function (l) { return l.split('=')[0]; });
    keys.slice(0, 2).should.eql(['SEQUENCE_ID', 'SEQUENCE_TEMPLATE']);
    keys.slice(-3).should.eql([w[0], 'SEQUENCE_TARGET', ''], run.name);
    tags[w[0]].should.equal(w[1]);
    tags.SEQUENCE_TARGET.should.equal(w[2]);
    tags.SEQUENCE_ID.should.equal(w[3]);
    tags.PRIMER_PRODUCT_SIZE_RANGE.should.equal(w[4]);
    tags.SEQUENCE_TEMPLATE.should.have.length(w[5]);
    tags.should.containEql({ PRIMER_TASK: 'generic', PRIMER_FIRST_BASE_INDEX: '1', PRIMER_EXPLAIN_FLAG: '1', PRIMER_NUM_RETURN: '20' });
    ['PRIMER_PICK_ANYWAY', 'SEQUENCE_INCLUDED_REGION', 'SEQUENCE_EXCLUDED_REGION', 'PRIMER_MAX_NS_ACCEPTED', 'PRIMER_LOWERCASE_MASKING']
      .forEach(function (k) { tags.should.not.have.property(k); });
    if (/_L1$/.test(run.name)) tags.should.containEql({ PRIMER_MAX_SIZE: '32', PRIMER_MIN_TM: '55', PRIMER_MAX_TM: '65', PRIMER_MIN_GC: '20', PRIMER_MAX_GC: '80' });
    cases.readRecording(run.output).should.match(/\n=\n$/);
  });
});

test('§4.10 scoring records: check_primers, PICK_ANYWAY, the two primers, widened pair constraints; one base substituted for a mismatch run', function () {
  const index = cases.loadIndex();
  const scoringRuns = Object.keys(index.runs).map(function (h) { return index.runs[h]; }).filter(function (r) { return /_p\d+_(alt|mmref|mmalt)$/.test(r.name); });
  scoringRuns.should.have.length(96 - 13);
  const kinds = { alt: 0, mmref: 0, mmalt: 0 };
  scoringRuns.forEach(function (run) {
    const text = cases.readRecording(run.input);
    const tags = boulder.parse(text);
    const keys = text.split('\n').filter(Boolean).map(function (l) { return l.split('=')[0]; });
    const kind = /_(alt|mmref|mmalt)$/.exec(run.name)[1];
    kinds[kind]++;
    keys.slice(-4).should.eql(['PRIMER_PICK_ANYWAY', 'SEQUENCE_PRIMER', 'SEQUENCE_PRIMER_REVCOMP', ''], run.name);
    tags.should.containEql({ PRIMER_TASK: 'check_primers', PRIMER_PICK_ANYWAY: '1', PRIMER_NUM_RETURN: '1', PRIMER_MIN_SIZE: '15', PRIMER_MAX_SIZE: '36',
      PRIMER_PAIR_MAX_DIFF_TM: '30' });
    tags.PRIMER_PRODUCT_SIZE_RANGE.should.equal(/as_pcr/.test(run.name) ? '36-400' : /_L0_/.test(run.name) ? '36-220' : '36-250');
    ['SEQUENCE_TARGET', 'SEQUENCE_FORCE_LEFT_END', 'SEQUENCE_FORCE_RIGHT_END', 'PRIMER_MAX_NS_ACCEPTED'].forEach(function (k) { tags.should.not.have.property(k); });
    tags.SEQUENCE_ID.should.match(new RegExp('_p\\d+_' + kind + '$'));
    const out = boulder.parse(cases.readRecording(run.output));
    should(out.PRIMER_ERROR).be.undefined();
    out.PRIMER_PAIR_NUM_RETURNED.should.equal('1');
  });
  kinds.should.eql({ alt: 16 + 16 + 16 + 8 + 9, mmref: 9, mmalt: 9 });
  // the mismatch run's template differs from the REF haplotype at exactly the base opposite the -2 mismatch
  const ref = boulder.parse(cases.readRecording('rs871475760_as_pcr_reverse_L0.input.txt')).SEQUENCE_TEMPLATE;
  const mm = boulder.parse(cases.readRecording('rs871475760_as_pcr_reverse_L0_p0_mmref.input.txt'));
  const diffs = [];
  for (let i = 0; i < ref.length; i++) if (ref[i] !== mm.SEQUENCE_TEMPLATE[i]) diffs.push([i + 1, ref[i], mm.SEQUENCE_TEMPLATE[i]]);
  diffs.should.eql([[402, 'T', 'C']]);
  mm.SEQUENCE_PRIMER_REVCOMP.should.equal('ATCTTTGACTAGCGAGAAATTCGG');
  scoring.scoringParams({ max_size: 30, product_size_ranges: [[61, 120], [200, 250]], opt_tm: 60 }).should.eql(
    { max_size: 36, product_size_ranges: [[36, 350]], opt_tm: 60, min_size: 15, max_tm_diff: 30, num_return: 1 });
});

// ---- §2.9 ------------------------------------------------------------------------------------------------------

test('§2.9 rs871475760 kasp: the complete response deep-equals the spec (template sequences elided as the spec prints them)', async function () {
  const { res, deps } = await replay('rs871475760_kasp');
  const got = JSON.parse(JSON.stringify(res));
  const window = cases.windowSeq();
  got.template.seq.should.equal(window.slice(10709 - 10500, 11509 - 10500 + 1));
  got.template.alt_seq.should.equal(got.template.seq.slice(0, 400) + 'A' + got.template.seq.slice(401));
  got.template.seq = got.template.seq.slice(0, 13) + '…' + got.template.seq.slice(-8);
  got.template.alt_seq = got.template.alt_seq.slice(0, 13) + '…' + got.template.alt_seq.slice(-8);
  got.should.eql(expected('rs871475760_kasp'));
  // key order as printed
  Object.keys(got).should.eql(Object.keys(expected('rs871475760_kasp')));
  Object.keys(got.sets[0]).should.eql(['id', 'key', 'rank', 'orientation', 'relaxation_level', 'quality', 'score', 'primers', 'products', 'thermo',
    'tm_balance', 'neighbour_sites', 'primer3_penalty', 'warnings', 'issues', 'check', 'order']);
  Object.keys(got.sets[0].primers.as_alt).should.eql(Object.keys(expected('rs871475760_kasp').sets[0].primers.as_alt));
  Object.keys(res).should.not.containEql('candidates');
  deps.semaphore.order.should.eql(['acquire', 'release']);
});

test('§2.9 rs871475760 kasp: the pairs of S1 (reverse, Primer3 pair 0) and S2 (forward, pair 2), and forward pair 0', async function () {
  const { res } = await replay('rs871475760_kasp');
  res.candidates.forward.map(function (c) { return c.pair_index; }).should.eql([0, 1, 2, 3, 4, 5, 6, 7]);
  res.candidates.reverse.map(function (c) { return c.pair_index; }).should.eql([0, 1, 2, 3, 4, 5, 6, 8]);
  const s1 = expectCandidate(res, 'reverse', 0, { level: 0, as: ['ATCTTTGACTAGCGAGAAATTCAG', 11109, 11132],
    common: ['AGCTTCTCTAAGTGGTTATCCGA', 11068, 11090], size: 65, penalty: 7.179713, as_tm: 57.10, common_tm: 58.72 });
  s1.as.end_stability.should.equal(3.02);
  s1.common.hairpin_th.should.be.approximately(34.99, 0.005);
  s1.common.end_stability.should.equal(4.55);
  s1.params.should.equal(res.settings.params);
  s1.scored.key.should.equal('f9df650ad116');
  const s2 = expectCandidate(res, 'forward', 2, { level: 0, as: ['GGTTATCCGAATATAGTCATACTCTATTC', 11081, 11109],
    common: ['TCTTTGTCTACTGAGAAATCCAGA', 11149, 11172], size: 92, penalty: 14.634369, as_tm: 57.28, common_tm: 57.08 });
  s2.as.hairpin_th.should.be.approximately(39.2, 0.005);
  s2.as.self_any_th.should.be.approximately(12.01, 0.005);
  s2.common.hairpin_th.should.be.approximately(35.39, 0.005);
  s2.scored.key.should.equal('1accc54c262d');
  // §4.16 S1: pair 0 of the forward run is the 30-mer with penalty 14.376148, not S2's pair.
  const p0 = res.candidates.forward[0];
  p0.as.seq.should.equal('TGGTTATCCGAATATAGTCATACTCTATTC');
  p0.pair.penalty.should.equal(14.376148);
  res.candidates.reverse.forEach(function (c) { c.common_neighbours_3p.should.eql([]); });
  res.budget.primer3_runs.should.equal(18);
  res.budget.exhausted.should.eql([]);
});

// ---- §2.10 -----------------------------------------------------------------------------------------------------

test('§2.10(a) tmp_1_11193_C_T: an EMS target needing level 1 in both orientations, with the four explain strings and the spec values', async function () {
  const { res } = await replay('tmp_1_11193_C_T_kasp');
  subsetDiff(JSON.parse(JSON.stringify(res)), expected('tmp_1_11193_C_T_kasp')).should.eql([]);
  res.template.should.containEql({ region: '1', start: 10793, end: 11593, length: 801 });
  res.assay.ems_target.should.be.true();
  res.variant.records.should.eql([{ id: 'tmp_1_11193_C_T', source: 'EMS_PMID38100514_Jiao', ems: true }]);
  const L1 = { max_size: 32, min_tm: 55, max_tm: 65, min_gc: 20, max_gc: 80, product_size_ranges: [[65, 150]] };
  ['forward', 'reverse'].forEach(function (o) {
    res.orientations[o].should.containEql({ status: 'ok', reason: null, discriminating_position: 11193, relaxation_level: 1, sets_found: 8, blockers: [] });
    res.orientations[o].attempts.map(counts).should.eql([
      { level: 0, pairs_returned: 0, rejected: rejected(), not_scored: 0, sets: 0 },
      { level: 1, pairs_returned: 20, rejected: rejected(), not_scored: 12, sets: 8 }
    ]);
    res.orientations[o].attempts.map(function (a) { return a.changes; }).should.eql([{}, L1]);
  });
  res.orientations.forward.attempts.map(raws).should.eql([
    { left: 'considered 13, GC content failed 8, low tm 3, ok 2', right: 'considered 4771, GC content failed 1241, low tm 1100, high tm 927, ok 1503',
      pair: 'considered 3006, unacceptable product size 3006, ok 0' },
    { left: 'considered 15, low tm 9, ok 6', right: 'considered 5490, GC content failed 293, low tm 1517, high tm 807, ok 2873',
      pair: 'considered 6436, unacceptable product size 6413, ok 23' }
  ]);
  res.orientations.reverse.attempts.map(raws).should.eql([
    { left: 'considered 4771, GC content failed 950, low tm 1792, high tm 507, ok 1522', right: 'considered 13, GC content failed 13, ok 0',
      pair: 'considered 0, ok 0' },
    { left: 'considered 5490, GC content failed 104, low tm 1994, high tm 435, ok 2957', right: 'considered 15, GC content failed 1, low tm 12, ok 2',
      pair: 'considered 1364, unacceptable product size 1341, tm diff too large 2, ok 21' }
  ]);
  res.warnings.should.eql([
    { code: 'EMS_TARGET', message: 'the target is an EMS mutation, private to BTx623-background mutant lines; natural neighbours are reported as warnings rather than blocking an orientation',
      details: { sources: ['EMS_PMID38100514_Jiao'] } },
    { code: 'RELAXED_CONSTRAINTS', message: 'forward: sets need relaxation level 1', details: { orientation: 'forward', level: 1, changes: L1 } },
    { code: 'RELAXED_CONSTRAINTS', message: 'reverse: sets need relaxation level 1', details: { orientation: 'reverse', level: 1, changes: L1 } }
  ]);
  res.neighbours.should.eql({ data: 'ensembl', window: { start: 10793, end: 11593 }, variants: 54, non_ems: 45, ems: 9, dense_non_ems: 1 });
  res.variant.submission_sequence.should.equal('GCATATTCTGGATTTCTCWGTAGACAAAGATAGATAACARAAATAGCTCT[C/T]TAGAGTATACACAATATATTAAAAAGTTGTTAGAGAGTGAAAATATATAG');
  expectCandidate(res, 'forward', 0, { level: 1, as: ['ACAAAGATAGATAACAAAAATAGCTCTC', 11166, 11193],
    common: ['TCACACTTTGAATCATCATTTGGA', 11268, 11291], size: 126, penalty: 14.470497, as_tm: 56.43, common_tm: 57.10 });
  expectCandidate(res, 'reverse', 4, { level: 1, as: ['ACAACTTTTTAATATATTGTGTATACTCTAG', 11193, 11223],
    common: ['TGAATTTCTCGCTAGTCAAAGA', 11110, 11131], size: 114, penalty: 18.278417, as_tm: 55.21, common_tm: 55.51 });
});

test('§2.10(b) rs5413864115 entered in Ensembl style: the shiftable deletion, ids filled from Ensembl, and the spec values', async function () {
  const { res } = await replay('rs5413864115_kasp');
  subsetDiff(JSON.parse(JSON.stringify(res)), expected('rs5413864115_kasp')).should.eql([]);
  res.variant.should.containEql({ key: '1:11282:CA:C', requested_id: null, ids: ['rs5413864115'], label: '1:11283-11283 A/-', kind: 'deletion', shift: 2 });
  res.variant.vcf.should.eql({ position: 11282, ref: 'CA', alt: 'C' });
  res.variant.zone.should.eql({ start: 11282, end: 11286 });
  res.template.should.containEql({ start: 10883, end: 11685, length: 803, alt_length: 802 });
  res.template.features.should.eql({ variant: { start: 400, end: 401 }, zone: { start: 400, end: 404 },
    discriminating: { forward: 403, reverse: 401 }, alt_offset: -1, exempt: [[367, 38], [400, 38]] });
  res.template.alt_seq.should.equal(res.template.seq.slice(0, 399) + 'C' + res.template.seq.slice(401));
  res.orientations.forward.should.containEql({ status: 'ok', discriminating_position: 11285, relaxation_level: 1, sets_found: 8 });
  res.orientations.reverse.should.containEql({ status: 'ok', discriminating_position: 11283, relaxation_level: 1, sets_found: 8 });
  res.orientations.forward.attempts.map(counts).should.eql([
    { level: 0, pairs_returned: 0, rejected: rejected(), not_scored: 0, sets: 0 },
    { level: 1, pairs_returned: 20, rejected: rejected({ common_neighbour_3p: 11 }), not_scored: 1, sets: 8 }
  ]);
  res.orientations.reverse.attempts.map(counts).should.eql([
    { level: 0, pairs_returned: 0, rejected: rejected(), not_scored: 0, sets: 0 },
    { level: 1, pairs_returned: 20, rejected: rejected(), not_scored: 12, sets: 8 }
  ]);
  res.warnings.map(function (w) { return w.code; }).should.eql(['SHIFTABLE_INDEL', 'DENSE_NEIGHBOURS', 'RELAXED_CONSTRAINTS', 'RELAXED_CONSTRAINTS']);
  res.warnings[0].should.eql({ code: 'SHIFTABLE_INDEL',
    message: 'the deletion can slide 2 bases inside AAA; the forward and reverse primers end at different bases and the wrong-allele primer may still prime through a 1-nt bulge',
    details: { shift: 2, forward_position: 11285, reverse_position: 11283, zone: { start: 11282, end: 11286 } } });
  res.warnings[1].should.eql({ code: 'DENSE_NEIGHBOURS', message: '4 known non-EMS variants lie within 30 bp of the variant', details: { count: 4, window: 30 } });
  res.variant.submission_sequence.should.equal('AAAATATATAGAAAACAATTTTATACAGATGATTTTCCAAATGATGATTC[A/]AAGTGTGAAATTTGRAAAGWCTCTTRGASATGMTYTAAGTGGAAGGAACA');
  expectCandidate(res, 'forward', 1, { level: 1, as: ['ACAGATGATTTTCCAAATGATGATTCAAA', 11257, 11285],
    common: ['CCCCATGTTTTTGTTCCTTCCA', 11323, 11344], size: 88, penalty: 8.47932, as_tm: 59.22, common_tm: 59.30 });
  expectCandidate(res, 'reverse', 0, { level: 1, as: ['AGAGTCTTTTCAAATTTCACACTTT', 11283, 11307],
    common: ['ACAAAAATAGCTCTCTAGAGTATACACA', 11179, 11206], size: 129, penalty: 14.793151, as_tm: 56.09, common_tm: 58.12 });
});

test('§2.10(c) tmp_1_11502_C_CGT: forward blocked by rs5413863234 at distance 4, reverse at level 0; ids merged; the spec values', async function () {
  const { res, calls } = await replay('tmp_1_11502_C_CGT_kasp');
  subsetDiff(JSON.parse(JSON.stringify(res)), expected('tmp_1_11502_C_CGT_kasp')).should.eql([]);
  res.variant.ids.should.eql(['tmp_1_11502_C_CGT', 'rs5413863549']);
  res.template.should.containEql({ start: 11102, end: 11903, length: 802, alt_length: 804 });
  res.orientations.forward.should.eql({
    status: 'blocked', reason: 'neighbour_at_3p', discriminating_position: 11503, relaxation_level: null, sets_found: 0,
    blockers: [{ key: '1:11500:G:A', ids: ['rs5413863234'], label: '1:11500 G/A', start: 11500, end: 11500, alleles: 'G/A', ems: false, distance_from_3p: 4 }],
    attempts: []
  });
  res.orientations.reverse.should.containEql({ status: 'ok', reason: null, discriminating_position: 11502, relaxation_level: 0, sets_found: 8, blockers: [] });
  res.orientations.reverse.attempts.map(counts).should.eql([
    { level: 0, pairs_returned: 20, rejected: rejected({ common_neighbour_3p: 6 }), not_scored: 6, sets: 8 }
  ]);
  // one design run (no forward run at all), then the ALT scoring run of each of the 8 scored reverse candidates
  calls.map(function (c) { return c.name; }).should.eql(['tmp_1_11502_C_CGT_kasp_reverse_L0'].concat(res.candidates.reverse.map(function (c) {
    return 'tmp_1_11502_C_CGT_kasp_reverse_L0_p' + c.pair_index + '_alt';
  })));
  calls.should.have.length(9);
  res.warnings.should.eql([
    { code: 'DUPLICATE_VARIANT_IDS', message: '2 Ensembl ids describe the same event 1:11502:C:CGT; they were merged',
      details: { key: '1:11502:C:CGT', ids: ['tmp_1_11502_C_CGT', 'rs5413863549'] } },
    { code: 'ORIENTATION_BLOCKED', message: 'forward: known variant rs5413863234 (G/A) lies 4 nt from the allele-specific primer\'s 3′ end',
      details: { orientation: 'forward', ids: ['rs5413863234'], distances: [4] } },
    { code: 'DENSE_NEIGHBOURS', message: '5 known non-EMS variants lie within 30 bp of the variant', details: { count: 5, window: 30 } }
  ]);
  res.sets[0].warnings[0].should.eql({ code: 'AS_TM_IMBALANCE', message: 'as_ref and as_alt Tm differ by 1.16 °C (limit 1.0)', details: { diff: 1.16 } });
  res.variant.submission_sequence.should.equal('ATGTTAGGATCTTTGCAACCCWGTGTTGCGTGCAATCTCGGTATCTCRCC[/GT]ATATGAYGTTAGGWTTTCTTWTCCTGCAACTGCCAAGAGAAYAAATATAT');
  expectCandidate(res, 'reverse', 3, { level: 0, as: ['GCAGGAAAAGAAATCCTAACATCATATG', 11502, 11529],
    common: ['AGGATCTTTGCAACCCTGTGTT', 11458, 11479], size: 72, penalty: 7.237494, as_tm: 59.19, common_tm: 60.43 });
});

test('§2.10(d) rs871475760 as_pcr: preset, attempts, the mismatch blocks of S1; the EMS neighbour at the common 3′ end does not reject', async function () {
  const { res } = await replay('rs871475760_as_pcr');
  subsetDiff(JSON.parse(JSON.stringify(res)), expected('rs871475760_as_pcr')).should.eql([]);
  res.settings.params.should.eql({ opt_size: 24, min_size: 18, max_size: 30, opt_tm: 60, min_tm: 57, max_tm: 63, min_gc: 30, max_gc: 70,
    max_tm_diff: 3, max_poly_x: 4, product_size_ranges: [[150, 300]] });
  res.assay.should.eql({ type: 'as_pcr', orientation: 'both', tails: 'none', deliberate_mismatch: 'auto', mismatch_position: 2,
    num_sets: 1, max_relaxation: 2, neighbour_policy: 'avoid_3p', ems_target: false, kasp_mix: null });
  res.orientations.forward.attempts.map(counts).should.eql([
    { level: 0, pairs_returned: 20, rejected: rejected({ common_neighbour_3p: 19 }), not_scored: 0, sets: 1 }
  ]);
  res.orientations.reverse.attempts.map(counts).should.eql([
    { level: 0, pairs_returned: 20, rejected: rejected({ common_neighbour_3p: 5 }), not_scored: 7, sets: 8 }
  ]);
  const s1 = expectCandidate(res, 'reverse', 0, { level: 0, as: ['ATCTTTGACTAGCGAGAAATTCAG', 11109, 11132],
    common: ['TGCATCAACAAATGTGCTATGTGT', 10880, 10903], size: 253, penalty: 2.922723, as_tm: 57.10, common_tm: 60.02 });
  s1.common_neighbours_3p.should.eql([]);
  res.sets[0].primers.as_ref.should.containEql({ tm_method: 'ntthal_duplex', tailed: null, dye: null, order_seq: 'ATCTTTGACTAGCGAGAAATTCGG' });
  should(res.sets[0].thermo.tailed).be.null();
  // §4.8 table: 2 design runs + 9 candidates x 3 scoring runs; 2 duplex calls per set, memoized across candidates.
  [res.budget.primer3_runs, res.budget.thermo_calls].should.eql([29, 4]);
  const counted = await replay('rs871475760_as_pcr', { scoreCandidate: countingScorer(3), thermo: NO_THERMO });
  counted.res.budget.primer3_runs.should.equal(29);
  counted.res.sets.should.eql([]);
  should(counted.res.check).be.null();
});

// ---- the scoring cap and the budget (§4.8) ----------------------------------------------------------------------------

test('§4.8 budget: max_primer3_runs 12 stops reverse after 2 sets with DESIGN_BUDGET_EXHAUSTED in a 200; the returned keys stay', async function () {
  const capped = await replay('rs871475760_kasp', { cfg: cases.cfg({ max_primer3_runs: 12 }) });
  const res = capped.res;
  res.orientations.forward.attempts.map(function (a) { return [a.pairs_returned, a.not_scored, a.sets]; }).should.eql([[20, 12, 8]]);
  res.orientations.reverse.attempts.map(function (a) { return [a.pairs_returned, a.not_scored, a.sets, a.rejected.common_neighbour_3p]; })
    .should.eql([[20, 18, 2, 0]]);
  res.orientations.reverse.should.containEql({ status: 'ok', relaxation_level: 0, sets_found: 2 });
  res.candidates.reverse.map(function (c) { return c.pair_index; }).should.eql([0, 1]);
  res.budget.should.containEql({ primer3_runs: 12, max_primer3_runs: 12, max_thermo_calls: 272 });
  res.budget.exhausted.should.eql(['reverse']);
  res.warnings.map(function (w) { return [w.code, w.details]; }).should.eql([['DESIGN_BUDGET_EXHAUSTED',
    { orientation: 'reverse', primer3_runs: 12, thermo_calls: 87, max_primer3_runs: 12, max_thermo_calls: 272, not_scored: 18, sets_returned: 2 }]]);
  res.sets.map(function (s) { return s.key; }).should.eql(['f9df650ad116', '1accc54c262d']);

  // the ntthal reservation: max_thermo_calls 20 scores one forward candidate (15 calls) and nothing else
  const thin = (await replay('rs871475760_kasp', { cfg: cases.cfg({ max_thermo_calls: 20 }) })).res;
  thin.orientations.forward.attempts.map(counts).should.eql([{ level: 0, pairs_returned: 20, rejected: rejected(), not_scored: 19, sets: 1 }]);
  thin.orientations.reverse.should.containEql({ status: 'no_sets', reason: 'budget_exhausted', sets_found: 0 });
  thin.orientations.reverse.attempts.map(counts).should.eql([{ level: 0, pairs_returned: 20, rejected: rejected(), not_scored: 20, sets: 0 }]);
  thin.sets.map(function (s) { return [s.id, s.orientation, s.score]; }).should.eql([['S1', 'forward', 22.3]]);
  thin.warnings.map(function (w) { return [w.code, w.details]; }).should.eql([
    ['ORIENTATION_NO_SETS', { orientation: 'reverse', levels_tried: 1 }],
    ['DESIGN_BUDGET_EXHAUSTED', { orientation: 'forward', primer3_runs: 3, thermo_calls: 15, max_primer3_runs: 54, max_thermo_calls: 20, not_scored: 19, sets_returned: 1 }],
    ['DESIGN_BUDGET_EXHAUSTED', { orientation: 'reverse', primer3_runs: 3, thermo_calls: 15, max_primer3_runs: 54, max_thermo_calls: 20, not_scored: 20, sets_returned: 1 }]
  ]);
  thin.check.set_ids.should.eql(['S1']);
});

test('reservation, createBudget and exact decimal comparison', function () {
  gd.reservation({ tails: 'ref_fam_alt_hex', deliberate_mismatch: 'none' }).should.eql({ primer3_runs: 1, thermo_calls: 15 });
  gd.reservation({ tails: 'none', deliberate_mismatch: 'auto' }).should.eql({ primer3_runs: 3, thermo_calls: 2 });
  gd.reservation({ tails: 'ref_hex_alt_fam', deliberate_mismatch: 'auto' }).should.eql({ primer3_runs: 3, thermo_calls: 17 });
  gd.reservation({ tails: 'none', deliberate_mismatch: 'none' }).should.eql({ primer3_runs: 1, thermo_calls: 0 });

  const b = gd.createBudget({ max_primer3_runs: 54, max_thermo_calls: 272 });
  b.countPrimer3(53);
  b.fits(1, 272).should.be.true();
  b.fits(2, 0).should.be.false();
  b.countThermo(258);
  b.fits(1, 14).should.be.true();
  b.fits(1, 15).should.be.false();
  let calls = 0;
  const shared = gd.createBudget({ max_primer3_runs: 54, max_thermo_calls: 20 }, { thermoCalls: function () { return calls; } });
  calls = 6;
  shared.thermo_calls.should.equal(6);
  shared.fits(0, 15).should.be.false();
  (function () { gd.createBudget({ max_primer3_runs: 54 }); }).should.throw(TypeError);

  gd.compareDecimal('51.999', 52).should.equal(-1);
  gd.compareDecimal(52, '52.000').should.equal(0);
  gd.compareDecimal(57.099, 52).should.equal(1);
  gd.compareDecimal('-0.125', '-0.12').should.equal(-1);
  (function () { gd.compareDecimal('1e3', 2); }).should.throw(TypeError);
});

// ---- the pair filters, floors and ladder on synthetic Primer3 output ------------------------------------------------------------

test('§4.8 filters in rank order: force, overlap, floor, duplicate, a scoring drop; pairs_returned = rejected + not_scored + sets', async function () {
  const ctx = await contextFor('rs871475760_kasp');
  try {
    const seq = ctx.template.seq;
    const scored = [];
    ctx.primer3 = primer3Stub(synthetic(seq, [
      { left: [372, 29], right: [464, 24] }, // 3' end at 400, not the forced 401
      { left: [373, 29], right: [423, 24] }, // common primer 400-423 = 1:11108-11131 covers the variant
      { left: [373, 29, '51.999'], right: [464, 24] }, // allele-specific primer under the 52 C floor
      { left: [373, 29, '58.100'], right: [464, 24] }, // the same sequences again
      { left: [372, 30, '52.000'], right: [464, 24] }, // the scorer drops it
      { left: [374, 28], right: [464, 24] }
    ]));
    ctx.scoreCandidate = async function (c) {
      scored.push(c.pair_index);
      return scored.length === 1 ? { dropped: 'alt_scoring_failed' } : { set: c };
    };
    const out = await gd.runOrientation('forward', ctx);
    out.orientation.should.containEql({ status: 'ok', reason: null, relaxation_level: 0, sets_found: 1 });
    out.orientation.attempts.map(counts).should.eql([{ level: 0, pairs_returned: 6,
      rejected: rejected({ force: 1, overlap: 1, below_floor: 1, duplicate: 1, alt_scoring_failed: 1 }), not_scored: 0, sets: 1 }]);
    scored.should.eql([4, 5]);
    out.sets.map(function (c) { return [c.pair_index, c.level, c.as.len]; }).should.eql([[5, 0, 28]]);
    out.budget_exhausted.should.be.false();

    ctx.scoreCandidate = async function () { return { dropped: 'too_cold' }; };
    ctx.budget = gd.createBudget(ctx.cfg.genotyping);
    let err = null;
    try {
      await gd.runOrientation('forward', ctx);
    } catch (e) {
      err = e;
    }
    should(err).be.instanceOf(TypeError);
    err.message.should.match(/too_cold/);
  } finally {
    ctx.dispose();
  }
});

test('§4.7/§4.10: a check_primers PRIMER_ERROR drops the candidate as alt_scoring_failed, never a 500; the run is counted', async function () {
  const ctx = await contextFor('rs871475760_kasp');
  try {
    const design1 = synthetic(ctx.template.seq, [{ left: [373, 29], right: [464, 24] }]);
    const seen = [];
    ctx.primer3 = primer3Stub(function (tags) {
      seen.push(tags.PRIMER_TASK);
      return tags.PRIMER_TASK === 'check_primers'
        ? { tags: { PRIMER_ERROR: 'Specified right primer not in sequence' }, error: 'Specified right primer not in sequence', warning: null }
        : design1;
    });
    ctx.thermo = NO_THERMO;
    ctx.scoreCandidate = scoring.scoreCandidate;
    const out = await gd.runOrientation('forward', ctx);
    out.orientation.attempts.map(counts).should.eql([
      { level: 0, pairs_returned: 1, rejected: rejected({ alt_scoring_failed: 1 }), not_scored: 0, sets: 0 },
      { level: 1, pairs_returned: 1, rejected: rejected({ duplicate: 1 }), not_scored: 0, sets: 0 },
      { level: 2, pairs_returned: 1, rejected: rejected({ duplicate: 1 }), not_scored: 0, sets: 0 }
    ]);
    out.orientation.should.containEql({ status: 'no_sets', reason: null });
    seen.should.eql(['generic', 'check_primers', 'generic', 'generic']);
    ctx.budget.primer3_runs.should.equal(4);
  } finally {
    ctx.dispose();
  }
});

test('§4.6 floors and ladder: every candidate under the floor runs L0-L2 and ends no_sets below_floor', async function () {
  const ctx = await contextFor('rs871475760_kasp');
  try {
    const calls = [];
    ctx.primer3 = primer3Stub(synthetic(ctx.template.seq, [
      { left: [373, 29, '51.900'], right: [464, 24] },
      { left: [372, 30, '45.400'], right: [464, 24] }
    ]), calls);
    const out = await gd.runOrientation('forward', ctx);
    out.orientation.should.containEql({ status: 'no_sets', reason: 'below_floor', relaxation_level: null, sets_found: 0 });
    out.orientation.attempts.map(counts).should.eql([
      { level: 0, pairs_returned: 2, rejected: rejected({ below_floor: 2 }), not_scored: 0, sets: 0 },
      { level: 1, pairs_returned: 2, rejected: rejected({ duplicate: 2 }), not_scored: 0, sets: 0 },
      { level: 2, pairs_returned: 2, rejected: rejected({ duplicate: 2 }), not_scored: 0, sets: 0 }
    ]);
    out.orientation.attempts.map(function (a) { return a.changes; }).should.eql([{}, ctx.req.ladder[0].changes, { min_tm: 52, max_tm_diff: 6 }]);
    calls.map(function (t) { return [t.PRIMER_MIN_TM, t.PRIMER_MAX_SIZE, t.PRIMER_PAIR_MAX_DIFF_TM, t.PRIMER_PRODUCT_SIZE_RANGE]; })
      .should.eql([[57, 30, 3, '61-120'], [55, 32, 3, '65-150'], [52, 32, 6, '65-150']]);
    calls.forEach(function (t) { t.should.not.have.property('PRIMER_PICK_ANYWAY'); });
    gd.orientationWarnings({ forward: out, reverse: { orientation: { status: 'skipped', reason: 'not_requested', sets_found: 0 } } })
      .map(function (w) { return [w.code, w.details]; })
      .should.eql([['ORIENTATION_NO_SETS', { orientation: 'forward', levels_tried: 3 }], ['NO_SETS', {}]]);
  } finally {
    ctx.dispose();
  }
});

test('§4.8 step 6: a candidate that does not fit the budget stops the ladder; no sets gives reason budget_exhausted', async function () {
  const ctx = await contextFor('rs871475760_kasp', { ctx: { budget: gd.createBudget({ max_primer3_runs: 1, max_thermo_calls: 272 }) } });
  try {
    ctx.primer3 = primer3Stub(synthetic(ctx.template.seq, [
      { left: [372, 29], right: [464, 24] },
      { left: [373, 29], right: [464, 24] },
      { left: [374, 28], right: [464, 24] }
    ]));
    const out = await gd.runOrientation('forward', ctx);
    out.orientation.should.containEql({ status: 'no_sets', reason: 'budget_exhausted', sets_found: 0 });
    out.orientation.attempts.map(counts).should.eql([{ level: 0, pairs_returned: 3, rejected: rejected({ force: 1 }), not_scored: 2, sets: 0 }]);
    out.budget_exhausted.should.be.true();
    ctx.budget.exhausted.should.eql(['forward']);
    ctx.budget.primer3_runs.should.equal(1);
  } finally {
    ctx.dispose();
  }
});

test('§4.15 blockers: avoid_3p blocks a natural target without a Primer3 run; ignore and an EMS target do not', async function () {
  const noRun = primer3Stub(function () { throw new Error('Primer3 must not run for a blocked orientation'); });
  let ctx = await contextFor('tmp_1_11502_C_CGT_kasp', { ctx: { primer3: noRun } });
  try {
    const out = await gd.runOrientation('forward', ctx);
    out.orientation.status.should.equal('blocked');
    out.orientation.blockers.map(function (b) { return [b.ids[0], b.distance_from_3p]; }).should.eql([['rs5413863234', 4]]);
  } finally {
    ctx.dispose();
  }
  const empty = synthetic('', []);
  for (const over of [{ body: { assay: { neighbour_policy: 'ignore', num_sets: 1 } } }, { ctx: {} }]) {
    const calls = [];
    ctx = await contextFor('tmp_1_11502_C_CGT_kasp', over);
    try {
      ctx.primer3 = primer3Stub(empty, calls);
      if (!over.body) ctx.variant = Object.assign({}, ctx.variant, { ems: true });
      const out = await gd.runOrientation('forward', ctx);
      out.orientation.should.containEql({ status: 'no_sets', reason: null, blockers: [] });
      calls.should.have.length(3);
      calls[0].SEQUENCE_FORCE_LEFT_END.should.equal(402);
      calls[0].SEQUENCE_TARGET.should.equal('403,10');
    } finally {
      ctx.dispose();
    }
  }
});

test('§4.8 step 3: neighbour_policy ignore keeps the pairs a common 3′ neighbour would reject (rs871475760 reverse)', async function () {
  const b = cases.body('rs871475760_kasp');
  b.assay.neighbour_policy = 'ignore';
  const { res } = await replay('rs871475760_kasp', { scoreCandidate: countingScorer(1), thermo: NO_THERMO }, b);
  res.orientations.reverse.attempts.map(counts).should.eql([{ level: 0, pairs_returned: 20, rejected: rejected(), not_scored: 12, sets: 8 }]);
  res.candidates.reverse.map(function (c) { return c.pair_index; }).should.eql([0, 1, 2, 3, 4, 5, 6, 7]);
  const kept = res.candidates.reverse.find(function (c) { return c.pair_index === 7; });
  kept.common_neighbours_3p.length.should.be.above(0);
  kept.common_neighbours_3p.forEach(function (n) {
    n.ems.should.be.false();
    n.distance_from_3p.should.be.within(1, 5);
  });
  res.orientations.forward.attempts.map(counts).should.eql([{ level: 0, pairs_returned: 20, rejected: rejected(), not_scored: 12, sets: 8 }]);
});

// ---- orchestration ---------------------------------------------------------------------------------------------------

function maskStub(calls, runs, fraction) {
  return {
    repeatMask: async function (tpl, opts) {
      calls.push({ template: tpl, opts: opts });
      return { mask: runs, mask_source: 'blast_depth', masked: true, masked_fraction: fraction, seq: repeatMask.applyMask(tpl.seq, runs, opts.mode),
        mode: opts.mode, tags: repeatMask.maskTags(opts.mode), warnings: [{ code: 'BLAST_DEPTH_MASK', message: repeatMask.MESSAGES.BLAST_DEPTH_MASK }] };
    }
  };
}

test('§4.5 mask exemption: each orientation\'s Primer3 template has its allele-specific window unmasked; VARIANT_IN_REPEAT', async function () {
  const maskCalls = [];
  const p3 = [];
  const b = Object.assign(cases.body('rs871475760_kasp'), { avoid_repeats: true });
  const deps = cases.deps('rs871475760_kasp', { repeatMask: maskStub(maskCalls, [[300, 201]], 0.2509), primer3: primer3Stub(synthetic('', []), p3), thermo: NO_THERMO });
  const res = await gd.designGenotyping(b, deps);
  maskCalls.should.have.length(1);
  maskCalls[0].template.should.containEql({ mode: 'region', region: '1', start: 10709, end: 11509, strand: 1 });
  maskCalls[0].opts.mode.should.equal('n_mask');
  res.template.should.containEql({ masked: true, mask_source: 'blast_depth', mask: [[300, 201]], masked_fraction: 0.2509 });
  p3.should.have.length(6);
  const seq = res.template.seq;
  const allN = function (s) { return /^N+$/.test(s); };
  const fwd = p3[0].SEQUENCE_TEMPLATE;
  p3[0].SEQUENCE_ID.should.equal('sorghum_bicolor_1_10709-11509_forward_L0');
  p3[0].PRIMER_MAX_NS_ACCEPTED.should.equal(0);
  fwd.slice(0, 299).should.equal(seq.slice(0, 299));
  allN(fwd.slice(299, 364)).should.be.true();
  fwd.slice(364, 401).should.equal(seq.slice(364, 401));
  allN(fwd.slice(401, 500)).should.be.true();
  fwd.slice(500).should.equal(seq.slice(500));
  const rev = p3[3].SEQUENCE_TEMPLATE;
  p3[3].SEQUENCE_ID.should.equal('sorghum_bicolor_1_10709-11509_reverse_L0');
  allN(rev.slice(299, 400)).should.be.true();
  rev.slice(400, 437).should.equal(seq.slice(400, 437));
  allN(rev.slice(437, 500)).should.be.true();
  res.warnings.map(function (w) { return w.code; })
    .should.eql(['BLAST_DEPTH_MASK', 'VARIANT_IN_REPEAT', 'VARIANT_IN_REPEAT', 'ORIENTATION_NO_SETS', 'ORIENTATION_NO_SETS', 'NO_SETS']);
  res.warnings[1].details.should.eql({ orientation: 'forward', masked_bases: 37 });
  res.warnings[2].details.should.eql({ orientation: 'reverse', masked_bases: 37 });
  deps.semaphore.order.should.eql(['acquire', 'release']);
  res.sets.should.eql([]);
  should(res.check).be.null();

  // three_prime: lowercase outside the exemption, uppercase inside
  const p3b = [];
  const b2 = Object.assign(cases.body('rs871475760_kasp'), { avoid_repeats: true, repeat_mask_mode: 'three_prime', assay: { orientation: 'forward' } });
  await gd.designGenotyping(b2, cases.deps('rs871475760_kasp', { repeatMask: maskStub([], [[380, 10]], 0.0125), primer3: primer3Stub(synthetic('', []), p3b), thermo: NO_THERMO }));
  p3b.should.have.length(3);
  p3b[0].PRIMER_LOWERCASE_MASKING.should.equal(1);
  p3b[0].SEQUENCE_TEMPLATE.should.equal(seq);

  const vt = { seq: 'ACGTACGTAC', length: 10, features: { exempt: [[3, 4], [5, 3]] } };
  gtemplate.orientationTemplate(vt, 'acgtNNNNac', 'forward').should.eql({ seq: 'acGTACNNac', exempt: [3, 4], masked_bases: 4 });
  gtemplate.orientationTemplate(vt, null, 'reverse').should.eql({ seq: 'ACGTACGTAC', exempt: [5, 3], masked_bases: 0 });
  (function () { gtemplate.orientationTemplate(vt, 'ACGT', 'forward'); }).should.throw(RangeError);
});

test('template_only: no Primer3, no ntthal and no semaphore; neighbours and the submission string are still resolved (§4.1, §4.2)', async function () {
  const noRun = primer3Stub(function () { throw new Error('Primer3 must not run for template_only'); });
  let deps = cases.deps('rs871475760_kasp', { primer3: noRun, thermo: NO_THERMO });
  let res = await gd.designGenotyping(Object.assign(cases.body('rs871475760_kasp'), { template_only: true }), deps);
  should(res.orientations).be.null();
  res.sets.should.eql([]);
  should(res.check).be.null();
  res.engine.should.eql({ primer3: null, thermo: null, genotyping_design: '1', variation_source: 'ensembl 115' });
  res.template.length.should.equal(801);
  res.settings.preset.should.equal('kasp');
  res.neighbours.should.eql({ data: 'ensembl', window: { start: 10709, end: 11509 }, variants: 57, non_ems: 45, ems: 12, dense_non_ems: 0 });
  res.variant.submission_sequence.should.equal(expected('rs871475760_kasp').variant.submission_sequence);
  res.assay.kasp_mix.stock_uM.should.equal(100);
  res.warnings.should.eql([]);
  deps.semaphore.order.should.eql([]);

  const maskCalls = [];
  deps = cases.deps('rs871475760_kasp', { primer3: noRun, thermo: NO_THERMO, repeatMask: maskStub(maskCalls, [[1, 200]], 0.2497) });
  res = await gd.designGenotyping(Object.assign(cases.body('rs871475760_kasp'), { template_only: true, avoid_repeats: true }), deps);
  maskCalls.should.have.length(1);
  deps.semaphore.order.should.eql(['acquire', 'release']);
  res.template.should.containEql({ masked: true, mask: [[1, 200]], masked_fraction: 0.2497 });
  res.warnings.map(function (w) { return w.code; }).should.eql(['BLAST_DEPTH_MASK']);
});

test('Ensembl unavailable: a manual design is 200 with NEIGHBOURS_UNAVAILABLE and no ids; a design by id is 503 and never takes a slot', async function () {
  const manual = await replay('rs5413864115_kasp', { ensembl: 'down', scoreCandidate: countingScorer(1), thermo: NO_THERMO });
  const res = manual.res;
  res.variant.should.containEql({ key: '1:11282:CA:C', requested_id: null, ids: [], records: [] });
  res.neighbours.should.eql({ data: 'unavailable', window: { start: 10883, end: 11685 }, variants: 0, non_ems: 0, ems: 0, dense_non_ems: 0 });
  res.warnings.map(function (w) { return w.code; }).should.eql(['SHIFTABLE_INDEL', 'NEIGHBOURS_UNAVAILABLE', 'RELAXED_CONSTRAINTS', 'RELAXED_CONSTRAINTS']);
  res.warnings[1].details.should.eql({ reason: 'transport' });
  res.engine.variation_source.should.equal('ensembl 115');
  // without neighbours every flank base is the reference base: the §4.17 string with each IUPAC code resolved
  const IUPAC = { R: 'AG', Y: 'CT', S: 'CG', W: 'AT', K: 'GT', M: 'AC' };
  const withCodes = 'AAAATATATAGAAAACAATTTTATACAGATGATTTTCCAAATGATGATTC[A/]AAGTGTGAAATTTGRAAAGWCTCTTRGASATGMTYTAAGTGGAAGGAACA';
  const plain = res.variant.submission_sequence;
  plain.should.not.match(/[RYSWKM]/);
  plain.should.have.length(withCodes.length);
  for (let i = 0; i < withCodes.length; i++) {
    if (IUPAC[withCodes[i]]) IUPAC[withCodes[i]].should.containEql(plain[i]);
    else plain[i].should.equal(withCodes[i]);
  }
  manual.deps.semaphore.order.should.eql(['acquire', 'release']);

  const deps = cases.deps('rs871475760_kasp', { ensembl: 'down', primer3: primer3Stub(function () { throw new Error('Primer3 must not run'); }), thermo: NO_THERMO });
  const err = await rejectsWith(gd.designGenotyping(cases.body('rs871475760_kasp'), deps), 503, 'VARIATION_SOURCE_UNAVAILABLE');
  err.details.should.match({ reason: 'transport' });
  deps.semaphore.order.should.eql([]);

  // a genome without variation data (no primers.variation.species entry; same sequence here): an id is 422 before any
  // Primer3 run, a manual variant designs with neighbours data "none" and warning NO_VARIATION_DATA
  const rio = async function (name) { return stubs.resolvedStub({ system_name: name }); };
  const noData = cases.deps('rs871475760_kasp', { resolve: rio, primer3: primer3Stub(function () { throw new Error('Primer3 must not run'); }), thermo: NO_THERMO });
  const e422 = await rejectsWith(gd.designGenotyping(Object.assign(cases.body('rs871475760_kasp'), { system_name: 'sorghum_rio' }), noData), 422, 'NO_VARIATION_DATA');
  e422.details.should.eql({ system_name: 'sorghum_rio' });
  noData.semaphore.order.should.eql([]);
  const noDataManual = cases.deps('rs871475760_kasp', { resolve: rio, thermo: NO_THERMO });
  const m = await gd.designGenotyping({ system_name: 'sorghum_rio', variant: { region: '1', position: 11109, ref: 'C', alt: 'A' }, template_only: true }, noDataManual);
  m.neighbours.should.eql({ data: 'none', window: { start: 10709, end: 11509 }, variants: 0, non_ems: 0, ems: 0, dense_non_ems: 0 });
  m.warnings.map(function (w) { return [w.code, w.details]; }).should.eql([['NO_VARIATION_DATA', { system_name: 'sorghum_rio' }]]);
  m.variant.should.containEql({ key: '1:11109:C:A', ids: [], requested_id: null });
  should(m.engine.variation_source).be.null();
});

test('orchestration errors: deadline 504, PRIMER_ERROR 400, PRIMER3_WARNING, disabled 503, unknown keys 400; the slot is always released', async function () {
  let deps = cases.deps('rs871475760_kasp', { cfg: cases.cfg(null, { design: { deadline_ms: 100 } }), thermo: NO_THERMO,
    primer3: { run: function () { return new Promise(function () {}); }, version: async function () { return '2.6.1'; } } });
  await rejectsWith(gd.designGenotyping(cases.body('rs871475760_kasp'), deps), 504, 'DEADLINE_EXCEEDED');
  deps.semaphore.order.should.eql(['acquire', 'release']);

  deps = cases.deps('rs871475760_kasp', { thermo: NO_THERMO, primer3: primer3Stub({ tags: { PRIMER_ERROR: 'SEQUENCE_TARGET beyond end of sequence' },
    error: 'SEQUENCE_TARGET beyond end of sequence', warning: null }) });
  const err = await rejectsWith(gd.designGenotyping(cases.body('rs871475760_kasp'), deps), 400, 'PRIMER3_INPUT_ERROR');
  err.details.should.eql({ primer3_error: 'SEQUENCE_TARGET beyond end of sequence' });
  deps.semaphore.order.should.eql(['acquire', 'release']);

  const warned = [];
  deps = cases.deps('rs871475760_kasp', { thermo: NO_THERMO, primer3: primer3Stub(synthetic('', [], 'unrecognized tag')) });
  const res = await gd.designGenotyping(Object.assign(cases.body('rs871475760_kasp'), { assay: { orientation: 'forward', max_relaxation: 0 } }), deps);
  res.warnings.forEach(function (w) { warned.push([w.code, w.details]); });
  warned.should.eql([
    ['PRIMER3_WARNING', { orientation: 'forward', level: 0 }],
    ['ORIENTATION_NO_SETS', { orientation: 'forward', levels_tried: 1 }],
    ['NO_SETS', {}]
  ]);
  res.orientations.reverse.should.containEql({ status: 'skipped', reason: 'not_requested', attempts: [] });

  await rejectsWith(gd.designGenotyping(cases.body('rs871475760_kasp'), cases.deps('rs871475760_kasp', { cfg: cases.cfg(null, { enabled: false }) })),
    503, 'FEATURE_DISABLED');
  await rejectsWith(gd.designGenotyping(Object.assign(cases.body('rs871475760_kasp'), { mode: 'region' }), cases.deps('rs871475760_kasp')),
    400, 'INVALID_REQUEST');

  // an injected resolver replaces the variation path
  const own = cases.deps('rs871475760_kasp', { thermo: NO_THERMO, primer3: primer3Stub(synthetic('', [])) });
  const seenReq = [];
  own.resolveVariant = async function (req, opts) {
    seenReq.push([req.system_name, opts.signal instanceof AbortSignal]);
    return cases.resolve(req, cases.deps('rs871475760_kasp'));
  };
  (await gd.designGenotyping(Object.assign(cases.body('rs871475760_kasp'), { template_only: true }), own)).variant.key.should.equal('1:11109:C:A');
  seenReq.should.eql([['sorghum_bicolor', true]]);
});

// ---- the template (§4.3) ------------------------------------------------------------------------------------------------

const C_DNA = '/fake/synthetic/dna/Syn.dna.toplevel.fa.gz';

// A synthetic region 'c' with an SNV at `position`; mutate(seq) may edit the random sequence first.
function syntheticWorld(length, position, mutate) {
  let seq = fx.randomSeq(length, 11).toUpperCase();
  if (mutate) seq = mutate(seq);
  const ref = seq[position - 1];
  const alt = ref === 'A' ? 'C' : 'A';
  const rec = normalize.parseManual({ region: 'c', position: position, ref: ref, alt: alt });
  const variant = normalize.recordsToEntries([rec], normalize.sequenceWindow(seq, 1, length), { region: 'c' }).entries[0];
  const files = {};
  files[C_DNA] = { c: { length: length, windows: [{ start: 1, seq: seq }] } };
  const cfg = cases.cfg();
  const req = request.normalize({ system_name: 'synthetic', variant: { region: 'c', position: position, ref: ref, alt: alt } }, cfg);
  return {
    seq: seq, variant: variant, resolved: stubs.resolvedStub({ system_name: 'synthetic', fasta: { dna: C_DNA, dna_sm: null } }),
    deps: { cfg: cfg, req: req, sequence: stubs.sequenceStub(files), log: stubs.silentLog }
  };
}

test('buildVariantTemplate: flank, ALT haplotype, features; one-sided room skips an orientation, no room at all is 400', async function () {
  template.buildRegionTemplate.should.be.a.Function();
  template.requireRegionLength.should.be.a.Function();
  const w = syntheticWorld(1000, 500);
  const vt = await gtemplate.buildVariantTemplate(w.variant, w.resolved, w.deps);
  vt.should.containEql({ system_name: 'synthetic', region: 'c', start: 100, end: 900, strand: 1, length: 801, alt_length: 801, region_length: 1000,
    product_min: 61, skip: { forward: null, reverse: null } });
  vt.seq.should.equal(w.seq.slice(99, 900));
  vt.alt_seq.should.equal(vt.seq.slice(0, 400) + w.variant.vcf.alt + vt.seq.slice(401));
  vt.region_template.should.containEql({ mode: 'region', id: 'synthetic_c_100-900', start: 100, end: 900 });
  gtemplate.templateWindow(w.variant, 1000, gtemplate.templateFlank(w.deps.req.levels, 400)).should.eql({ start: 100, end: 900 });

  // flank = max(400, largest product-range upper bound + 40)
  const wide = syntheticWorld(3000, 1500);
  wide.deps.req = request.normalize({ system_name: 'synthetic', variant: { region: 'c', position: 1500, ref: wide.variant.vcf.ref, alt: wide.variant.vcf.alt },
    params: { product_size_ranges: [[100, 700]] } }, wide.deps.cfg);
  const wt = await gtemplate.buildVariantTemplate(wide.variant, wide.resolved, wide.deps);
  [wt.start, wt.end].should.eql([760, 2240]);
  gtemplate.templateFlank(wide.deps.req.levels, 400).should.equal(740);

  const nearEnd = syntheticWorld(540, 500);
  const ne = await gtemplate.buildVariantTemplate(nearEnd.variant, nearEnd.resolved, nearEnd.deps);
  [ne.start, ne.end].should.eql([100, 540]);
  ne.skip.should.eql({ forward: 'too_close_to_end', reverse: null });
  const ctx = { cfg: nearEnd.deps.cfg, req: nearEnd.deps.req, variant: nearEnd.variant, template: ne, neighbours: [],
    budget: gd.createBudget(nearEnd.deps.cfg.genotyping), deadline: design.createDeadline(1000), primer3: primer3Stub(synthetic('', [])) };
  try {
    (await gd.runOrientation('forward', ctx)).orientation.should.containEql({ status: 'skipped', reason: 'too_close_to_end', attempts: [] });
  } finally {
    ctx.deadline.dispose();
  }

  const tiny = syntheticWorld(100, 50);
  const err = await rejectsWith(gtemplate.buildVariantTemplate(tiny.variant, tiny.resolved, tiny.deps), 400, 'VARIANT_TOO_CLOSE_TO_END');
  err.details.should.eql({ region: 'c', position: 50, region_length: 100, needed: 61 });
  tiny.deps.sequence.calls.should.eql([]);
});

test('buildVariantTemplate: an N in the allele-specific window skips that orientation; REF, repetitive and no-sequence errors', async function () {
  const withN = syntheticWorld(1000, 500, function (s) { return s.slice(0, 479) + 'N' + s.slice(480); });
  (await gtemplate.buildVariantTemplate(withN.variant, withN.resolved, withN.deps)).skip.should.eql({ forward: 'n_in_primer_window', reverse: null });
  const far = syntheticWorld(1000, 500, function (s) { return s.slice(0, 460) + 'N' + s.slice(461); });
  (await gtemplate.buildVariantTemplate(far.variant, far.resolved, far.deps)).skip.should.eql({ forward: null, reverse: null });

  const w = syntheticWorld(1000, 500);
  const wrong = JSON.parse(JSON.stringify(w.variant));
  wrong.vcf.ref = wrong.vcf.ref === 'G' ? 'T' : 'G';
  const e1 = await rejectsWith(gtemplate.buildVariantTemplate(wrong, w.resolved, w.deps), 400, 'REF_MISMATCH');
  e1.details.should.eql({ region: 'c', position: 500, given: wrong.vcf.ref, genome: w.variant.vcf.ref });
  const repetitive = Object.assign({}, w.variant, { zone: null, discriminating: null, shift: 1001 });
  const e2 = await rejectsWith(gtemplate.buildVariantTemplate(repetitive, w.resolved, w.deps), 400, 'VARIANT_TOO_REPETITIVE');
  e2.details.should.eql({ region: 'c', position: 500, shift: 1001, max: 1000 });
  await rejectsWith(gtemplate.buildVariantTemplate(w.variant, stubs.resolvedStub({ fasta: { dna: null } }), w.deps), 422, 'NO_SEQUENCE');
  const other = Object.assign({}, w.variant, { region: 'd' });
  await rejectsWith(gtemplate.buildVariantTemplate(other, w.resolved, w.deps), 404, 'UNKNOWN_REGION');
});

// ---- real binaries -----------------------------------------------------------------------------------------------------

const realPrimer3 = { skip: REALDATA && fs.existsSync(PRIMER3_BIN) ? false : 'set PRIMERS_REALDATA=1 (needs ' + PRIMER3_BIN + ')' };
const realBoth = { skip: REALDATA && fs.existsSync(PRIMER3_BIN) && fs.existsSync(NTTHAL_BIN) ? false : 'set PRIMERS_REALDATA=1 (needs primer3_core and ntthal)' };

test('real primer3_core reproduces every recorded output (PRIMERS_REALDATA=1)', realPrimer3, async function () {
  const primer3 = require('../../../api/helpers/primers/primer3');
  const index = cases.loadIndex();
  for (const hash of Object.keys(index.runs)) {
    const run = index.runs[hash];
    const r = await primer3.run(cases.readRecording(run.input), { bin: PRIMER3_BIN, timeoutMs: 30000, maxStdoutBytes: 2000000, log: stubs.silentLog });
    should(r.error).be.null();
    r.tags.should.eql(boulder.parse(cases.readRecording(run.output)), run.name);
  }
});

test('real binaries: the five designs reproduce their set keys and scores; check_primers on REF reproduces the design run (±0.01)', realBoth, async function () {
  const primer3 = require('../../../api/helpers/primers/primer3');
  const thermo = require('../../../api/helpers/primers/thermo');
  const want = { rs871475760_kasp: [['f9df650ad116', 9.18], ['1accc54c262d', 19.63]], tmp_1_11193_C_T_kasp: [['f5a5c6f2ebcf', 22.51], ['c791a4956fd3', 24.84]],
    rs5413864115_kasp: [['a5277232d8ab', 15.48], ['cb3ef66afd37', 23.95]], tmp_1_11502_C_CGT_kasp: [['53942cb55348', 9.56]], rs871475760_as_pcr: [['7f9af6b1c938', 7.92]] };
  for (const name of Object.keys(want)) {
    const config = cases.cfg(null, { primer3_core: PRIMER3_BIN, ntthal: NTTHAL_BIN });
    const started = Date.now();
    const res = await gd.designGenotyping(cases.body(name), cases.deps(name, { cfg: config, primer3: primer3, thermo: thermo.createThermo({ cfg: config }) }));
    (Date.now() - started).should.be.below(3000);
    res.sets.map(function (s) { return [s.key, s.score]; }).should.eql(want[name], name);
    if (name === 'rs871475760_as_pcr') continue;
    for (const c of res.candidates.forward.concat(res.candidates.reverse)) {
      const forward = c.orientation === 'forward';
      const dl = design.createDeadline(30000);
      try {
        const run = await scoring.scorePair({ id: 'equivalence', seq: res.template.seq }, forward ? c.as.seq : c.common.seq, forward ? c.common.seq : c.as.seq,
          c.params, { cfg: config, primer3: primer3, budget: gd.createBudget({ max_primer3_runs: 99, max_thermo_calls: 0 }), deadline: dl, log: stubs.silentLog });
        ['left', 'right'].forEach(function (side) {
          const d = c.pair[side];
          const s = run.pair[side];
          ['tm', 'hairpin_th', 'self_any_th', 'self_end_th', 'end_stability'].forEach(function (k) { s[k].should.be.approximately(d[k], 0.01, name + ' ' + side + ' ' + k); });
        });
        run.pair.product_size.should.equal(c.pair.product_size);
      } finally {
        dl.dispose();
      }
    }
  }
});
