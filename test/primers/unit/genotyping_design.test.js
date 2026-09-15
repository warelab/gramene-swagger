'use strict';

require('../../../api/helpers/primers/node_compat');

const test = require('node:test');
const should = require('should');
const fs = require('fs');

const gd = require('../../../api/helpers/primers/genotyping/design');
const gtemplate = require('../../../api/helpers/primers/genotyping/template');
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

// Spec §4.1-§4.8 and the M5 rows of genotyping_design in §7.3: records, explain strings, pairs_returned, not_scored,
// the scoring cap and budget reservation, the guard, floors and mask exemption. No primer3_core is spawned outside the
// PRIMERS_REALDATA=1 test.
// Recorded fixtures: test/primers/fixtures/primer3/genotyping/ (record_genotyping.js, primer3_core 2.6.1, 2026-09-15).
// Apart from SEQUENCE_ID, every recorded input and output is byte-identical to the design fix pass recordings of the
// spec (genotyping-design/scratch-fix2/score2/default/p3), which the §2.9-§2.10 values were computed from; the attempt
// numbers pinned below for §2.10(b) and (d), which the spec elides, are that pass's.

const REALDATA = process.env.PRIMERS_REALDATA === '1';
const PRIMER3_BIN = '/home/olson/bin/primer3_core';

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
  const deps = cases.deps(name, Object.assign({ primer3: cases.recordedPrimer3(calls) }, extra || {}));
  const res = await gd.designGenotyping(body || cases.body(name), deps);
  return { res: res, calls: calls, deps: deps };
}

// A scorer that keeps every candidate and counts `runs` check_primers runs for it, as scoring would.
function countingScorer(runs) {
  return async function (candidate, ctx) {
    ctx.budget.countPrimer3(runs);
    return { set: candidate };
  };
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

// ctx for runOrientation on one case; over: {body (merged into the case body), ctx (merged into ctx)}. Call dispose().
async function contextFor(name, over) {
  over = over || {};
  const deps = cases.deps(name);
  const req = request.normalize(Object.assign(cases.body(name), over.body || {}), deps.cfg);
  const resolved = await deps.resolveVariant(req);
  const vt = await gtemplate.buildVariantTemplate(resolved.variant, resolved.assembly, Object.assign({}, deps, { req: req }));
  const dl = design.createDeadline(45000);
  const ctx = Object.assign({
    cfg: deps.cfg, req: req, variant: resolved.variant, template: vt, maskedSeq: null, neighbours: resolved.neighbours.entries,
    budget: gd.createBudget(deps.cfg.genotyping), deadline: dl, log: stubs.silentLog
  }, over.ctx || {});
  ctx.dispose = function () { dl.dispose(); };
  return ctx;
}

// ---- the recordings --------------------------------------------------------------------------------------------

test('recorded fixtures: every example regenerates its Primer3 design records byte for byte', async function () {
  const index = cases.loadIndex();
  index.primer3_version.should.equal('2.6.1');
  Object.keys(index.runs).should.have.length(13);
  for (const name of Object.keys(cases.CASES)) {
    const r = await replay(name);
    r.calls.map(function (c) { return c.name; }).should.eql(index.cases[name]);
    r.calls.forEach(function (c) {
      c.input.should.equal(cases.readRecording(c.name + '.input.txt'));
      index.runs[cases.sha256(c.input)].name.should.equal(c.name);
    });
  }
  // Any other input is refused.
  const stub = cases.recordedPrimer3();
  const tags = boulder.parse(cases.readRecording('rs871475760_kasp_forward_L0.input.txt'));
  tags.PRIMER_NUM_RETURN = '10';
  let refused = null;
  try { await stub.run(tags); } catch (e) { refused = e; }
  should.exist(refused);
  refused.code.should.equal('UNRECORDED_PRIMER3_INPUT');
});

test('§4.7 records: buildRecord, then the forced 3′ end, the SEQUENCE_TARGET guard and PRIMER_NUM_RETURN 20; never PICK_ANYWAY', function () {
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
  Object.keys(index.runs).map(function (h) { return index.runs[h].name; }).sort().should.eql(Object.keys(want).sort());
  Object.keys(index.runs).forEach(function (hash) {
    const run = index.runs[hash];
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

// ---- §2.9 ------------------------------------------------------------------------------------------------------

test('§2.9 rs871475760 kasp: template, settings, assay and both orientations with their Primer3 explain data', async function () {
  const { res } = await replay('rs871475760_kasp');
  res.template.should.containEql({ system_name: 'sorghum_bicolor', region: '1', start: 10709, end: 11509, strand: 1, length: 801,
    alt_length: 801, masked: false, mask_source: null, mask: [], masked_fraction: 0 });
  res.template.seq.slice(0, 13).should.equal('GCGAGTTCTCAAG');
  res.template.seq.slice(-8).should.equal('CATATGAT');
  res.template.seq[400].should.equal('C');
  res.template.alt_seq.should.equal(res.template.seq.slice(0, 400) + 'A' + res.template.seq.slice(401));
  res.template.features.should.eql({ variant: { start: 401, end: 401 }, zone: { start: 401, end: 401 },
    discriminating: { forward: 401, reverse: 401 }, alt_offset: 0, exempt: [[365, 37], [401, 37]] });
  res.assay.should.eql({ type: 'kasp', orientation: 'both', tails: 'ref_fam_alt_hex', deliberate_mismatch: 'none', mismatch_position: 2,
    num_sets: 2, max_relaxation: 2, neighbour_policy: 'avoid_3p', ems_target: false });
  res.settings.should.eql({
    preset: 'kasp',
    params: { opt_size: 22, min_size: 18, max_size: 30, opt_tm: 60, min_tm: 57, max_tm: 63, min_gc: 30, max_gc: 70,
      max_tm_diff: 3, max_poly_x: 5, product_size_ranges: [[61, 120]] },
    pinned: [],
    ladder: [
      { level: 1, changes: { max_size: 32, min_tm: 55, max_tm: 65, min_gc: 20, max_gc: 80, product_size_ranges: [[65, 150]] } },
      { level: 2, changes: { min_tm: 52, max_tm_diff: 6 } }
    ],
    floors: { as_min_tm: 52, as_min_gc: 15 }
  });
  res.orientations.should.eql({
    forward: {
      status: 'ok', reason: null, discriminating_position: 11109, relaxation_level: 0, sets_found: 8, blockers: [],
      attempts: [
        { level: 0, changes: {},
          explain: {
            left: { raw: 'considered 13, GC content failed 5, low tm 6, ok 2', considered: 13, 'GC content failed': 5, 'low tm': 6, ok: 2 },
            right: { raw: 'considered 4771, GC content failed 1718, low tm 1123, high tm 777, ok 1153', considered: 4771, 'GC content failed': 1718, 'low tm': 1123, 'high tm': 777, ok: 1153 },
            pair: { raw: 'considered 1463, unacceptable product size 1443, ok 20', considered: 1463, 'unacceptable product size': 1443, ok: 20 } },
          pairs_returned: 20,
          rejected: { force: 0, overlap: 0, common_neighbour_3p: 0, alt_scoring_failed: 0, below_floor: 0, duplicate: 0 },
          not_scored: 12,
          sets: 8 }
      ]
    },
    reverse: {
      status: 'ok', reason: null, discriminating_position: 11109, relaxation_level: 0, sets_found: 8, blockers: [],
      attempts: [
        { level: 0, changes: {},
          explain: {
            left: { raw: 'considered 4771, GC content failed 598, low tm 1550, high tm 1031, ok 1592', considered: 4771, 'GC content failed': 598, 'low tm': 1550, 'high tm': 1031, ok: 1592 },
            right: { raw: 'considered 13, low tm 6, ok 7', considered: 13, 'low tm': 6, ok: 7 },
            pair: { raw: 'considered 905, unacceptable product size 883, tm diff too large 1, ok 21', considered: 905, 'unacceptable product size': 883, 'tm diff too large': 1, ok: 21 } },
          pairs_returned: 20,
          rejected: { force: 0, overlap: 0, common_neighbour_3p: 6, alt_scoring_failed: 0, below_floor: 0, duplicate: 0 },
          not_scored: 6,
          sets: 8 }
      ]
    }
  });
  res.warnings.should.eql([]);
  res.engine.should.eql({ primer3: '2.6.1', thermo: null, genotyping_design: '1', variation_source: 'ensembl 115' });
  res.sets.should.eql([]);
  should(res.check).be.null();
  res.variant.key.should.equal('1:11109:C:A');
  Object.keys(res).should.not.containEql('candidates');
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
  const s2 = expectCandidate(res, 'forward', 2, { level: 0, as: ['GGTTATCCGAATATAGTCATACTCTATTC', 11081, 11109],
    common: ['TCTTTGTCTACTGAGAAATCCAGA', 11149, 11172], size: 92, penalty: 14.634369, as_tm: 57.28, common_tm: 57.08 });
  s2.as.hairpin_th.should.be.approximately(39.2, 0.005);
  s2.as.self_any_th.should.be.approximately(12.01, 0.005);
  s2.common.hairpin_th.should.be.approximately(35.39, 0.005);
  // §4.16 S1: pair 0 of the forward run is the 30-mer with penalty 14.376148, not S2's pair.
  const p0 = res.candidates.forward[0];
  p0.as.seq.should.equal('TGGTTATCCGAATATAGTCATACTCTATTC');
  p0.pair.penalty.should.equal(14.376148);
  res.candidates.reverse.forEach(function (c) { c.common_neighbours_3p.should.eql([]); });
  res.budget.primer3_runs.should.equal(2);
  res.budget.exhausted.should.eql([]);
});

// ---- §2.10 -----------------------------------------------------------------------------------------------------

test('§2.10(a) tmp_1_11193_C_T: an EMS target needing level 1 in both orientations, with the four explain strings', async function () {
  const { res } = await replay('tmp_1_11193_C_T_kasp');
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
    { code: 'RELAXED_CONSTRAINTS', message: 'forward: sets need relaxation level 1', details: { orientation: 'forward', level: 1, changes: L1 } },
    { code: 'RELAXED_CONSTRAINTS', message: 'reverse: sets need relaxation level 1', details: { orientation: 'reverse', level: 1, changes: L1 } }
  ]);
  expectCandidate(res, 'forward', 0, { level: 1, as: ['ACAAAGATAGATAACAAAAATAGCTCTC', 11166, 11193],
    common: ['TCACACTTTGAATCATCATTTGGA', 11268, 11291], size: 126, penalty: 14.470497, as_tm: 56.43, common_tm: 57.10 });
  expectCandidate(res, 'reverse', 4, { level: 1, as: ['ACAACTTTTTAATATATTGTGTATACTCTAG', 11193, 11223],
    common: ['TGAATTTCTCGCTAGTCAAAGA', 11110, 11131], size: 114, penalty: 18.278417, as_tm: 55.21, common_tm: 55.51 });
});

test('§2.10(b) rs5413864115 entered in Ensembl style: the shiftable deletion template and level-1 pairs', async function () {
  const { res } = await replay('rs5413864115_kasp');
  res.variant.should.containEql({ key: '1:11282:CA:C', ids: ['rs5413864115'], label: '1:11283-11283 A/-', kind: 'deletion', shift: 2 });
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
  res.warnings.map(function (w) { return [w.code, w.details.orientation, w.details.level]; })
    .should.eql([['RELAXED_CONSTRAINTS', 'forward', 1], ['RELAXED_CONSTRAINTS', 'reverse', 1]]);
  expectCandidate(res, 'forward', 1, { level: 1, as: ['ACAGATGATTTTCCAAATGATGATTCAAA', 11257, 11285],
    common: ['CCCCATGTTTTTGTTCCTTCCA', 11323, 11344], size: 88, penalty: 8.47932, as_tm: 59.22, common_tm: 59.30 });
  expectCandidate(res, 'reverse', 0, { level: 1, as: ['AGAGTCTTTTCAAATTTCACACTTT', 11283, 11307],
    common: ['ACAAAAATAGCTCTCTAGAGTATACACA', 11179, 11206], size: 129, penalty: 14.793151, as_tm: 56.09, common_tm: 58.12 });
});

test('§2.10(c) tmp_1_11502_C_CGT: forward blocked by rs5413863234 at distance 4, reverse at level 0', async function () {
  const { res, calls } = await replay('tmp_1_11502_C_CGT_kasp');
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
  calls.map(function (c) { return c.name; }).should.eql(['tmp_1_11502_C_CGT_kasp_reverse_L0']);
  res.warnings.should.eql([{ code: 'ORIENTATION_BLOCKED',
    message: 'forward: known variant rs5413863234 (G/A) lies 4 nt from the allele-specific primer\'s 3′ end',
    details: { orientation: 'forward', ids: ['rs5413863234'], distances: [4] } }]);
  expectCandidate(res, 'reverse', 3, { level: 0, as: ['GCAGGAAAAGAAATCCTAACATCATATG', 11502, 11529],
    common: ['AGGATCTTTGCAACCCTGTGTT', 11458, 11479], size: 72, penalty: 7.237494, as_tm: 59.19, common_tm: 60.43 });
  res.budget.primer3_runs.should.equal(1);
});

test('§2.10(d) rs871475760 as_pcr: preset, attempts and the pair of S1; the EMS neighbour at the common 3′ end does not reject', async function () {
  const { res } = await replay('rs871475760_as_pcr', { scoreCandidate: countingScorer(3) });
  res.settings.params.should.eql({ opt_size: 24, min_size: 18, max_size: 30, opt_tm: 60, min_tm: 57, max_tm: 63, min_gc: 30, max_gc: 70,
    max_tm_diff: 3, max_poly_x: 4, product_size_ranges: [[150, 300]] });
  res.assay.should.eql({ type: 'as_pcr', orientation: 'both', tails: 'none', deliberate_mismatch: 'auto', mismatch_position: 2,
    num_sets: 1, max_relaxation: 2, neighbour_policy: 'avoid_3p', ems_target: false });
  res.orientations.forward.attempts.map(counts).should.eql([
    { level: 0, pairs_returned: 20, rejected: rejected({ common_neighbour_3p: 19 }), not_scored: 0, sets: 1 }
  ]);
  res.orientations.reverse.attempts.map(counts).should.eql([
    { level: 0, pairs_returned: 20, rejected: rejected({ common_neighbour_3p: 5 }), not_scored: 7, sets: 8 }
  ]);
  const s1 = expectCandidate(res, 'reverse', 0, { level: 0, as: ['ATCTTTGACTAGCGAGAAATTCAG', 11109, 11132],
    common: ['TGCATCAACAAATGTGCTATGTGT', 10880, 10903], size: 253, penalty: 2.922723, as_tm: 57.10, common_tm: 60.02 });
  s1.common_neighbours_3p.should.eql([]);
  // §4.8 table: 2 design runs + 9 candidates x 3 scoring runs.
  res.budget.primer3_runs.should.equal(29);
});

// ---- the scoring cap and the budget (§4.8) ----------------------------------------------------------------------------

test('§4.8 budget: default caps give 18 Primer3 runs; max_primer3_runs 12 stops reverse after 2 sets without an error', async function () {
  const full = await replay('rs871475760_kasp', { scoreCandidate: countingScorer(1) });
  full.res.budget.primer3_runs.should.equal(18);
  full.res.budget.exhausted.should.eql([]);

  const capped = await replay('rs871475760_kasp', { cfg: cases.cfg({ max_primer3_runs: 12 }), scoreCandidate: countingScorer(1) });
  const res = capped.res;
  res.orientations.forward.attempts.map(function (a) { return [a.pairs_returned, a.not_scored, a.sets]; }).should.eql([[20, 12, 8]]);
  res.orientations.reverse.attempts.map(function (a) { return [a.pairs_returned, a.not_scored, a.sets, a.rejected.common_neighbour_3p]; })
    .should.eql([[20, 18, 2, 0]]);
  res.orientations.reverse.should.containEql({ status: 'ok', relaxation_level: 0, sets_found: 2 });
  res.candidates.reverse.map(function (c) { return c.pair_index; }).should.eql([0, 1]);
  res.budget.should.containEql({ primer3_runs: 12, max_primer3_runs: 12, max_thermo_calls: 272 });
  res.budget.exhausted.should.eql(['reverse']);
  res.warnings.should.eql([]);
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
  const { res } = await replay('rs871475760_kasp', {}, b);
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
  const deps = cases.deps('rs871475760_kasp', { repeatMask: maskStub(maskCalls, [[300, 201]], 0.2509), primer3: primer3Stub(synthetic('', []), p3) });
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

  // three_prime: lowercase outside the exemption, uppercase inside
  const p3b = [];
  const b2 = Object.assign(cases.body('rs871475760_kasp'), { avoid_repeats: true, repeat_mask_mode: 'three_prime', assay: { orientation: 'forward' } });
  await gd.designGenotyping(b2, cases.deps('rs871475760_kasp', { repeatMask: maskStub([], [[380, 10]], 0.0125), primer3: primer3Stub(synthetic('', []), p3b) }));
  p3b.should.have.length(3);
  p3b[0].PRIMER_LOWERCASE_MASKING.should.equal(1);
  p3b[0].SEQUENCE_TEMPLATE.should.equal(seq);

  const vt = { seq: 'ACGTACGTAC', length: 10, features: { exempt: [[3, 4], [5, 3]] } };
  gtemplate.orientationTemplate(vt, 'acgtNNNNac', 'forward').should.eql({ seq: 'acGTACNNac', exempt: [3, 4], masked_bases: 4 });
  gtemplate.orientationTemplate(vt, null, 'reverse').should.eql({ seq: 'ACGTACGTAC', exempt: [5, 3], masked_bases: 0 });
  (function () { gtemplate.orientationTemplate(vt, 'ACGT', 'forward'); }).should.throw(RangeError);
});

test('template_only: no Primer3 run and no semaphore, unless avoid_repeats needs the masker (§4.1, §4.2)', async function () {
  const noRun = primer3Stub(function () { throw new Error('Primer3 must not run for template_only'); });
  let deps = cases.deps('rs871475760_kasp', { primer3: noRun });
  let res = await gd.designGenotyping(Object.assign(cases.body('rs871475760_kasp'), { template_only: true }), deps);
  should(res.orientations).be.null();
  res.sets.should.eql([]);
  should(res.check).be.null();
  should(res.engine.primer3).be.null();
  res.template.length.should.equal(801);
  res.settings.preset.should.equal('kasp');
  deps.semaphore.order.should.eql([]);

  const maskCalls = [];
  deps = cases.deps('rs871475760_kasp', { primer3: noRun, repeatMask: maskStub(maskCalls, [[1, 200]], 0.2497) });
  res = await gd.designGenotyping(Object.assign(cases.body('rs871475760_kasp'), { template_only: true, avoid_repeats: true }), deps);
  maskCalls.should.have.length(1);
  deps.semaphore.order.should.eql(['acquire', 'release']);
  res.template.should.containEql({ masked: true, mask: [[1, 200]], masked_fraction: 0.2497 });
});

test('orchestration errors: deadline 504, PRIMER_ERROR 400, disabled 503, no resolver 500; the slot is always released', async function () {
  let deps = cases.deps('rs871475760_kasp', { cfg: cases.cfg(null, { design: { deadline_ms: 100 } }),
    primer3: { run: function () { return new Promise(function () {}); }, version: async function () { return '2.6.1'; } } });
  await rejectsWith(gd.designGenotyping(cases.body('rs871475760_kasp'), deps), 504, 'DEADLINE_EXCEEDED');
  deps.semaphore.order.should.eql(['acquire', 'release']);

  deps = cases.deps('rs871475760_kasp', { primer3: primer3Stub({ tags: { PRIMER_ERROR: 'SEQUENCE_TARGET beyond end of sequence' },
    error: 'SEQUENCE_TARGET beyond end of sequence', warning: null }) });
  const err = await rejectsWith(gd.designGenotyping(cases.body('rs871475760_kasp'), deps), 400, 'PRIMER3_INPUT_ERROR');
  err.details.should.eql({ primer3_error: 'SEQUENCE_TARGET beyond end of sequence' });
  deps.semaphore.order.should.eql(['acquire', 'release']);

  const warned = [];
  deps = cases.deps('rs871475760_kasp', { primer3: primer3Stub(synthetic('', [], 'unrecognized tag')) });
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
  deps = cases.deps('rs871475760_kasp', { resolveVariant: undefined });
  await rejectsWith(gd.designGenotyping(cases.body('rs871475760_kasp'), deps), 500, 'INTERNAL');
  deps.semaphore.order.should.eql([]);
  await rejectsWith(gd.designGenotyping(Object.assign(cases.body('rs871475760_kasp'), { mode: 'region' }), cases.deps('rs871475760_kasp')),
    400, 'INVALID_REQUEST');
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

  // flank = max(400, largest product-range upper bound + 40)
  const wide = syntheticWorld(3000, 1500);
  wide.deps.req = request.normalize({ system_name: 'synthetic', variant: { region: 'c', position: 1500, ref: wide.variant.vcf.ref, alt: wide.variant.vcf.alt },
    params: { product_size_ranges: [[100, 700]] } }, wide.deps.cfg);
  const wt = await gtemplate.buildVariantTemplate(wide.variant, wide.resolved, wide.deps);
  [wt.start, wt.end].should.eql([760, 2240]);

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

// ---- real binary -----------------------------------------------------------------------------------------------------

test('real primer3_core reproduces every recorded output (PRIMERS_REALDATA=1)',
  { skip: REALDATA && fs.existsSync(PRIMER3_BIN) ? false : 'set PRIMERS_REALDATA=1 (needs ' + PRIMER3_BIN + ')' }, async function () {
    const primer3 = require('../../../api/helpers/primers/primer3');
    const index = cases.loadIndex();
    for (const hash of Object.keys(index.runs)) {
      const run = index.runs[hash];
      const r = await primer3.run(cases.readRecording(run.input), { bin: PRIMER3_BIN, timeoutMs: 30000, maxStdoutBytes: 2000000, log: stubs.silentLog });
      should(r.error).be.null();
      r.tags.should.eql(boulder.parse(cases.readRecording(run.output)), run.name);
    }
  });
