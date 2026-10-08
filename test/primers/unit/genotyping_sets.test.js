'use strict';

require('../../../api/helpers/primers/node_compat');

const { test, before } = require('node:test');
const should = require('should');

const gd = require('../../../api/helpers/primers/genotyping/design');
const sets = require('../../../api/helpers/primers/genotyping/sets');
const scoring = require('../../../api/helpers/primers/genotyping/scoring');
const order = require('../../../api/helpers/primers/genotyping/order');
const decimal = require('../../../api/helpers/primers/genotyping/decimal');
const cases = require('../fixtures/primer3/genotyping/cases');
const EXPECTED_29 = require('../fixtures/primer3/genotyping/expected_rs871475760_kasp.json');

// Spec §4.8-§4.18 and the genotyping_sets rows of §7.3: ALT derivation, set keys, discrimination, floors, issues, quality
// and score (rules S1-S7), exact arithmetic and rounding, the shift tract, ranking, the proposed check request and the
// order rows. The example sets come from replaying the recorded Primer3 and ntthal runs of the five §2.9-§2.10 designs
// (test/primers/fixtures/primer3/genotyping); nothing is spawned.

const R = {};

before(async function () {
  for (const name of Object.keys(cases.CASES)) {
    R[name] = await gd.designGenotyping(cases.body(name), cases.deps(name, { primer3: cases.recordedPrimer3() }));
  }
});

function candidate(res, orientation, pairIndex) {
  const c = res.candidates[orientation].find(function (x) { return x.pair_index === pairIndex; });
  should.exist(c, orientation + ' pair ' + pairIndex);
  return c;
}

// ---- §4.9 --------------------------------------------------------------------------------------------------------

test('deriveAlt and productSpans: the five §4.9 rows (sequences, blocks, inserted bases, product sizes)', function () {
  const rows = [
    ['rs871475760_kasp', 'reverse', 0, 'ATCTTTGACTAGCGAGAAATTCAT', [[11109, 11132]], 0, 65, 65, 0],
    ['rs871475760_kasp', 'forward', 2, 'GGTTATCCGAATATAGTCATACTCTATTA', [[11081, 11109]], 0, 92, 92, 0],
    ['rs5413864115_kasp', 'forward', 1, 'ACAGATGATTTTCCAAATGATGATTCAAG', [[11257, 11282], [11284, 11286]], 0, 88, 87, 0],
    ['rs5413864115_kasp', 'reverse', 0, 'AGAGTCTTTTCAAATTTCACACTTG', [[11282, 11282], [11284, 11307]], 0, 129, 128, 0],
    ['tmp_1_11502_C_CGT_kasp', 'reverse', 3, 'GCAGGAAAAGAAATCCTAACATCATATA', [[11503, 11529]], 1, 72, 74, 2]
  ];
  rows.forEach(function (row) {
    const res = R[row[0]];
    const c = candidate(res, row[1], row[2]);
    const alt = sets.deriveAlt(c, res.template, res.variant);
    alt.matched_seq.should.equal(row[3], row[0]);
    alt.genomic.blocks.map(function (b) { return [b.start, b.end]; }).should.eql(row[4]);
    alt.genomic.strand.should.equal(row[1] === 'forward' ? 1 : -1);
    alt.inserted_bases.should.equal(row[5]);
    alt.matched_seq.length.should.equal(c.as.seq.length);
    const p = sets.productSpans(c, res.template, res.variant);
    [p.ref.size, p.alt.size, p.alt.inserted_bases].should.eql([row[6], row[7], row[8]]);
    // the ALT check_primers run confirms the size: ref + alt_offset
    c.scored.products.alt.size.should.equal(row[7]);
    c.scored.products.ref.size.should.equal(row[6]);
  });
  const del = candidate(R.rs5413864115_kasp, 'forward', 1).scored.products;
  del.alt.genomic.blocks.should.eql([{ start: 11257, end: 11282 }, { start: 11284, end: 11344 }]);
  del.ref.genomic.should.eql({ region: '1', start: 11257, end: 11344, strand: 1, blocks: [{ start: 11257, end: 11344 }] });
});

test('setKey: the eight example keys', function () {
  [
    ['reverse', 'ATCTTTGACTAGCGAGAAATTCAG', 'ATCTTTGACTAGCGAGAAATTCAT', 'AGCTTCTCTAAGTGGTTATCCGA', 'f9df650ad116'],
    ['forward', 'GGTTATCCGAATATAGTCATACTCTATTC', 'GGTTATCCGAATATAGTCATACTCTATTA', 'TCTTTGTCTACTGAGAAATCCAGA', '1accc54c262d'],
    ['forward', 'ACAAAGATAGATAACAAAAATAGCTCTC', 'ACAAAGATAGATAACAAAAATAGCTCTT', 'TCACACTTTGAATCATCATTTGGA', 'f5a5c6f2ebcf'],
    ['reverse', 'ACAACTTTTTAATATATTGTGTATACTCTAG', 'ACAACTTTTTAATATATTGTGTATACTCTAA', 'TGAATTTCTCGCTAGTCAAAGA', 'c791a4956fd3'],
    ['forward', 'ACAGATGATTTTCCAAATGATGATTCAAA', 'ACAGATGATTTTCCAAATGATGATTCAAG', 'CCCCATGTTTTTGTTCCTTCCA', 'a5277232d8ab'],
    ['reverse', 'AGAGTCTTTTCAAATTTCACACTTT', 'AGAGTCTTTTCAAATTTCACACTTG', 'ACAAAAATAGCTCTCTAGAGTATACACA', 'cb3ef66afd37'],
    ['reverse', 'GCAGGAAAAGAAATCCTAACATCATATG', 'GCAGGAAAAGAAATCCTAACATCATATA', 'AGGATCTTTGCAACCCTGTGTT', '53942cb55348'],
    ['reverse', 'ATCTTTGACTAGCGAGAAATTCGG', 'ATCTTTGACTAGCGAGAAATTCGT', 'TGCATCAACAAATGTGCTATGTGT', '7f9af6b1c938']
  ].forEach(function (r) { sets.setKey(r[0], r[1], r[2], r[3]).should.equal(r[4]); });
  [R.rs871475760_kasp, R.tmp_1_11193_C_T_kasp, R.rs5413864115_kasp, R.tmp_1_11502_C_CGT_kasp, R.rs871475760_as_pcr].map(function (res) {
    return res.sets.map(function (s) { return s.key; });
  }).should.eql([['f9df650ad116', '1accc54c262d'], ['f5a5c6f2ebcf', 'c791a4956fd3'], ['a5277232d8ab', 'cb3ef66afd37'], ['53942cb55348'], ['7f9af6b1c938']]);
});

// ---- §4.12 --------------------------------------------------------------------------------------------------------

test('discrimination with realign + classify: [] likely and [1] likely_weak; with a -2 mismatch [2] likely and [2,1] unlikely', function () {
  const res = R.rs871475760_kasp;
  const span = { start: 401, end: 424 };
  const plain = sets.discrimination({
    orientation: 'reverse', variant: res.variant, template: res.template,
    as_ref: Object.assign({ target_seq: 'ATCTTTGACTAGCGAGAAATTCAG', matched_seq: 'ATCTTTGACTAGCGAGAAATTCAG' }, span),
    as_alt: Object.assign({ target_seq: 'ATCTTTGACTAGCGAGAAATTCAT', matched_seq: 'ATCTTTGACTAGCGAGAAATTCAT' }, span)
  });
  ['as_ref', 'as_alt'].forEach(function (role) {
    plain[role].should.eql({ own_allele: { mm_pos: [], likelihood: 'likely' }, other_allele: { mm_pos: [1], likelihood: 'likely_weak' },
      terminal_mismatch_class: 'max', in_shift_tract: false });
  });
  const mm = sets.discrimination({
    orientation: 'reverse', variant: res.variant, template: res.template,
    as_ref: Object.assign({ target_seq: 'ATCTTTGACTAGCGAGAAATTCGG', matched_seq: 'ATCTTTGACTAGCGAGAAATTCAG' }, span),
    as_alt: Object.assign({ target_seq: 'ATCTTTGACTAGCGAGAAATTCGT', matched_seq: 'ATCTTTGACTAGCGAGAAATTCAT' }, span)
  });
  ['as_ref', 'as_alt'].forEach(function (role) {
    mm[role].own_allele.should.eql({ mm_pos: [2], likelihood: 'likely' });
    mm[role].other_allele.should.eql({ mm_pos: [2, 1], likelihood: 'unlikely' });
  });
  // the AS-PCR design reports exactly these
  R.rs871475760_as_pcr.sets[0].primers.as_ref.discrimination.should.eql(mm.as_ref);
  // every KASP allele-specific primer of the examples, the homopolymer deletion included: other allele [1] likely_weak
  ['rs871475760_kasp', 'tmp_1_11193_C_T_kasp', 'rs5413864115_kasp', 'tmp_1_11502_C_CGT_kasp'].forEach(function (name) {
    R[name].sets.forEach(function (s) {
      ['as_ref', 'as_alt'].forEach(function (role) {
        s.primers[role].discrimination.other_allele.should.eql({ mm_pos: [1], likelihood: 'likely_weak' }, name + ' ' + s.id + ' ' + role);
      });
    });
  });
});

test('shift tract (§3.5, §4.12): rs5413864115 S1 and S2 as_ref inside, as_alt outside, one SHIFT_TRACT_DISCRIMINATION each; S2 23.95', function () {
  const res = R.rs5413864115_kasp;
  res.sets.map(function (s) {
    return [s.id, s.primers.as_ref.discrimination.in_shift_tract, s.primers.as_alt.discrimination.in_shift_tract,
      s.issues.filter(function (i) { return i.code === 'SHIFT_TRACT_DISCRIMINATION'; }).map(function (i) { return i.details; })];
  }).should.eql([['S1', true, false, [{ primer: 'as_ref', shift: 2 }]], ['S2', true, false, [{ primer: 'as_ref', shift: 2 }]]]);
  sets.shiftTractFlags(res.variant, 'forward').should.eql({ as_ref: true, as_alt: false });
  sets.shiftTractFlags(res.variant, 'reverse').should.eql({ as_ref: true, as_alt: false });
  sets.shiftTractFlags(R.tmp_1_11502_C_CGT_kasp.variant, 'reverse').should.eql({ as_ref: false, as_alt: false });
  res.sets[1].score.should.equal(23.95);
  decimal.format(candidate(res, 'reverse', 0).scored.exact).should.equal('23.953151');
});

// ---- §4.6 floors ----------------------------------------------------------------------------------------------------

test('floors: an ALT primer at 51.9 °C is dropped; at 54.89 °C with 19.35 % GC it is kept, with ALT_PRIMER_SUBOPTIMAL', function () {
  const floors = { as_min_tm: 52, as_min_gc: 15 };
  const s2alt = 'ACAACTTTTTAATATATTGTGTATACTCTAA';
  decimal.gcPercent(s2alt).should.equal(19.35);
  const ref = { matched_tm: 55.21, target_seq: 'ACAACTTTTTAATATATTGTGTATACTCTAG', derived: false };
  scoring.belowFloors({ as_ref: ref, as_alt: { matched_tm: 51.9, target_seq: s2alt, derived: true } }, floors).should.be.true();
  scoring.belowFloors({ as_ref: ref, as_alt: { matched_tm: '51.999999', target_seq: s2alt, derived: true } }, floors).should.be.true();
  scoring.belowFloors({ as_ref: ref, as_alt: { matched_tm: 52, target_seq: s2alt, derived: true } }, floors).should.be.false();
  scoring.belowFloors({ as_ref: ref, as_alt: { matched_tm: 54.89, target_seq: s2alt, derived: true } }, floors).should.be.false();
  // GC floor for derived primers only: 4/31 = 12.9 %
  const lowGc = 'ACAACTTTTTAATATATTATATATAATTTAA';
  scoring.belowFloors({ as_ref: ref, as_alt: { matched_tm: 54.89, target_seq: lowGc, derived: true } }, floors).should.be.true();
  scoring.belowFloors({ as_ref: { matched_tm: 55, target_seq: lowGc, derived: false }, as_alt: { matched_tm: 54.89, target_seq: s2alt, derived: true } }, floors)
    .should.be.false();

  const s2 = R.tmp_1_11193_C_T_kasp.sets[1];
  s2.primers.as_alt.should.containEql({ target_seq: s2alt, tm: 54.89, gc: 19.35, primer3_problems: ' GC content too low; Temperature too low;' });
  s2.warnings.should.eql([{ code: 'ALT_PRIMER_SUBOPTIMAL',
    message: 'the derived as_alt primer would not have been chosen de novo (Primer3: GC content too low, temperature too low); it clears the hard floors',
    details: { oligo: 'as_alt', problems: ' GC content too low; Temperature too low;' } }]);
  sets.problemPhrases(' Temperature too low;').should.equal('temperature too low');
});

// ---- §4.16 issues, quality and score ----------------------------------------------------------------------------------

test('issues, quality and score of the eight example sets (§4.16 table)', function () {
  const rows = [
    ['rs871475760_kasp', 'S1', 'f9df650ad116', 'reverse', 0, 'usable', 9.18, '9.179713', 7.1797, 0, 2],
    ['rs871475760_kasp', 'S2', '1accc54c262d', 'forward', 0, 'poor', 19.63, '19.634369', 14.6344, 1, 2],
    ['tmp_1_11193_C_T_kasp', 'S1', 'f5a5c6f2ebcf', 'forward', 1, 'usable', 22.51, '22.510497', 14.4705, 0, 2],
    ['tmp_1_11193_C_T_kasp', 'S2', 'c791a4956fd3', 'reverse', 1, 'usable', 24.84, '24.838417', 18.2784, 0, 1],
    ['rs5413864115_kasp', 'S1', 'a5277232d8ab', 'forward', 1, 'usable', 15.48, '15.479320', 8.4793, 0, 2],
    ['rs5413864115_kasp', 'S2', 'cb3ef66afd37', 'reverse', 1, 'usable', 23.95, '23.953151', 14.7932, 0, 3],
    ['tmp_1_11502_C_CGT_kasp', 'S1', '53942cb55348', 'reverse', 0, 'usable', 9.56, '9.561494', 7.2375, 0, 2],
    ['rs871475760_as_pcr', 'S1', '7f9af6b1c938', 'reverse', 0, 'usable', 7.92, '7.922723', 2.9227, 0, 5]
  ];
  rows.forEach(function (r) {
    const res = R[r[0]];
    const s = res.sets.find(function (x) { return x.id === r[1]; });
    [s.key, s.orientation, s.relaxation_level, s.quality, s.score, s.primer3_penalty].should.eql([r[2], r[3], r[4], r[5], r[6], r[8]], r[0] + ' ' + r[1]);
    const c = res.candidates[s.orientation].find(function (x) { return x.scored && x.scored.key === s.key; });
    decimal.format(c.scored.exact).should.equal(r[7]);
    c.scored.counted.should.eql({ high: r[9], warn: r[10] });
  });
  const codes = function (name, id) {
    return R[name].sets.find(function (x) { return x.id === id; }).issues.map(function (i) { return i.code + '/' + i.severity + (i.details.role ? ':' + i.details.role : ''); });
  };
  codes('rs871475760_kasp', 'S1').should.eql(['ALT_PRIMER_SUBOPTIMAL/warn', 'TAILED_STRUCTURE/warn']);
  codes('rs871475760_kasp', 'S2').should.eql(['ALT_PRIMER_SUBOPTIMAL/warn', 'TAILED_STRUCTURE/high', 'NEIGHBOUR_IN_PRIMER/warn:common']);
  // AS_TM_IMBALANCE is an issue but is priced by S3, not counted by S6
  codes('tmp_1_11502_C_CGT_kasp', 'S1').should.eql(['AS_TM_IMBALANCE/warn', 'NEIGHBOUR_IN_PRIMER/warn:as_ref', 'NEIGHBOUR_IN_PRIMER/warn:common', 'WEAK_TERMINAL_CLASS/info']);
  const s11502 = R.tmp_1_11502_C_CGT_kasp.sets[0];
  s11502.tm_balance.should.eql({ as_tm_diff: 1.16, common_minus_as: 1.24 });
  const c11502 = candidate(R.tmp_1_11502_C_CGT_kasp, 'reverse', 3).scored;
  decimal.format(c11502.terms.as_tm).should.equal('0.324000');
  // WEAK_TERMINAL_CLASS (A/G and C/T sites) is info: in issues only, never in warnings, never scored
  [['tmp_1_11193_C_T_kasp', 'S1'], ['tmp_1_11193_C_T_kasp', 'S2'], ['rs5413864115_kasp', 'S1'], ['tmp_1_11502_C_CGT_kasp', 'S1']].forEach(function (x) {
    const s = R[x[0]].sets.find(function (y) { return y.id === x[1]; });
    s.issues.filter(function (i) { return i.code === 'WEAK_TERMINAL_CLASS'; }).should.eql([{ code: 'WEAK_TERMINAL_CLASS', severity: 'info',
      message: 'both allele-specific 3′ mismatches are in the weak class (A/G or C/T site)', details: {} }]);
    s.warnings.map(function (w) { return w.code; }).should.not.containEql('WEAK_TERMINAL_CLASS');
  });
  // quality: high -> poor; warn or relaxation -> usable; otherwise good
  sets.scoreSet({ penalty: 1, relaxation_level: 0, balance: { as_tm_diff: 0n, common_out: 0n }, issues: [], kasp: true, product_size: 90, g: cases.GENOTYPING }).quality.should.equal('good');
  sets.scoreSet({ penalty: 1, relaxation_level: 1, balance: { as_tm_diff: 0n, common_out: 0n }, issues: [], kasp: true, product_size: 90, g: cases.GENOTYPING }).quality.should.equal('usable');
  sets.scoreSet({ penalty: 1, relaxation_level: 0, balance: { as_tm_diff: 0n, common_out: 0n }, kasp: true, product_size: 90, g: cases.GENOTYPING,
    issues: [{ code: 'COMMON_TM_OUT_OF_RANGE', severity: 'warn' }, { code: 'AS_TM_IMBALANCE', severity: 'high' }] }).should.containEql({ quality: 'poor', score: 1 });
});

test('score rules S1 and S6: the own pair penalty; one neighbour issue per primer site; ALT_PRIMER_SUBOPTIMAL only from ordered oligos', function () {
  // S1: rs871475760 S2 is Primer3's pair 2 (14.634369); pair 0 of that run has 14.376148
  const s2 = candidate(R.rs871475760_kasp, 'forward', 2);
  decimal.format(s2.scored.terms.penalty).should.equal('14.634369');
  candidate(R.rs871475760_kasp, 'forward', 0).pair.penalty.should.equal(14.376148);
  // S6 grouping: tmp_1_11193_C_T S1 has one allele-specific and one common NEIGHBOUR_IN_PRIMER (22.51, not 23.51)
  const t1 = R.tmp_1_11193_C_T_kasp.sets[0];
  t1.warnings.map(function (w) { return [w.code, w.details]; }).should.eql([
    ['NEIGHBOUR_IN_PRIMER', { role: 'as_ref', ids: ['rs873026643'], distances: [12] }],
    ['NEIGHBOUR_IN_PRIMER', { role: 'common', ids: ['rs5413864115'], distances: [16] }]
  ]);
  t1.primers.as_alt.neighbours.map(function (n) { return [n.ids[0], n.distance_from_3p]; }).should.eql([['rs873026643', 12]]);
  t1.score.should.equal(22.51);
  t1.neighbour_sites.should.equal(2);
  // rs5413864115 S1: a single common issue listing three neighbours
  R.rs5413864115_kasp.sets[0].warnings[0].details.should.eql({ role: 'common', ids: ['tmp_1_11334_A_G', 'tmp_1_11335_A_G', 'tmp_1_11340_TG_T'], distances: [12, 13, 19] });
  // rs871475760 AS-PCR S1: two ALT_PRIMER_SUBOPTIMAL, from the two mismatch runs; the matched ALT primer's own run also has
  // problems, which are not counted (7.92, not 8.92)
  const as = R.rs871475760_as_pcr.sets[0];
  as.issues.filter(function (i) { return i.code === 'ALT_PRIMER_SUBOPTIMAL'; }).map(function (i) { return i.details; }).should.eql([
    { oligo: 'as_alt', problems: ' Hairpin stability too high;' }, { oligo: 'as_ref', problems: ' Hairpin stability too high;' }
  ]);
  const altRun = cases.readRecording('rs871475760_as_pcr_reverse_L0_p0_alt.output.txt');
  altRun.should.match(/PRIMER_RIGHT_0_PROBLEMS= \S/);
  as.score.should.equal(7.92);
  as.primers.as_alt.primer3_problems.should.equal(' Hairpin stability too high;');
});

test('the eight scored forward rs871475760 candidates rank 19.63 (Primer3 pair 2) first and pair 0 (22.30) last', function () {
  const fwd = R.rs871475760_kasp.candidates.forward.map(function (c) { return c.scored; }).sort(sets.compareSets);
  fwd.map(function (s) { return s.score; }).should.eql([19.63, 19.79, 19.93, 21.53, 21.67, 21.76, 21.9, 22.3]);
  [fwd[0].pair_index, fwd[7].pair_index].should.eql([2, 0]);
  fwd[7].primers.as_ref.target_seq.should.equal('TGGTTATCCGAATATAGTCATACTCTATTC');
  fwd[7].issues.filter(function (i) { return i.code === 'TAILED_STRUCTURE' && i.severity === 'high'; }).should.have.length(2);
});

test('rankSets: the best set of each orientation first, then the rest by exact score and key; cut to num_sets', function () {
  const mk = function (o, exact, key) { return { orientation: o, exact: BigInt(exact), key: key }; };
  const f1 = mk('forward', 19634369, 'b');
  const f2 = mk('forward', 9000000, 'z');
  const f3 = mk('forward', 9000000, 'a');
  const r1 = mk('reverse', 9179713, 'c');
  const r2 = mk('reverse', 9200000, 'd');
  sets.rankSets([f1, f2, f3, r1, r2], 10).should.eql([f3, r1, f2, r2, f1]);
  sets.rankSets([f1, r1, r2], 2).should.eql([r1, f1]);
  sets.rankSets([r2, r1], 6).should.eql([r1, r2]);
  // §2.9: S1 reverse (9.18) before S2 forward (19.63), although forward has seven sets scoring below the eighth reverse set
  R.rs871475760_kasp.sets.map(function (s) { return [s.id, s.rank, s.orientation]; }).should.eql([['S1', 0, 'reverse'], ['S2', 1, 'forward']]);
});

// ---- §2.1 exact arithmetic -------------------------------------------------------------------------------------------

test('exact arithmetic and the one rounding rule (§2.1)', function () {
  const t = candidate(R.tmp_1_11193_C_T_kasp, 'forward', 0).scored;
  t.tm_balance.common_minus_as.should.equal(0.68);
  decimal.round(decimal.micro('57.102') - decimal.micro('56.427'), 2).should.equal(0.68);
  (Math.round((57.102 - 56.427) * 100) / 100).should.equal(0.67); // the float shortcut the rule forbids
  decimal.round(decimal.micro('0.125'), 2).should.equal(0.13);
  decimal.round(decimal.micro('-0.125'), 2).should.equal(-0.13);
  decimal.round(decimal.micro('2.675'), 2).should.equal(2.68);
  Number((2.675).toFixed(2)).should.equal(2.67); // binary floating point: 2.675 is stored as 2.67499999...
  decimal.round(decimal.micro('9.179713'), 2).should.equal(9.18);
  decimal.round(decimal.micro('23.953151'), 2).should.equal(23.95);
  decimal.ratio(600, 31, 2).should.equal(19.35);
  decimal.gcPercent('ACAACTTTTTAATATATTGTGTATACTCTAA').should.equal(19.35);
  decimal.roundValue(7.179713, 4).should.equal(7.1797);
  decimal.round(decimal.micro('-0.004'), 2).should.equal(0);
  decimal.round(decimal.micro('51.502395'), 1).should.equal(51.5);
  (decimal.micro('4.5500000') === 4550000n).should.be.true();
  (decimal.micro(-1.0) === -1000000n).should.be.true();
  (function () { decimal.micro('1e3'); }).should.throw(TypeError);
  (function () { decimal.micro('0.0000001'); }).should.throw(TypeError);
  (function () { decimal.micro(NaN); }).should.throw(TypeError);
  // 6/31 = 19.3548387...: exactly, not by its rounded 19.35
  decimal.gcAtLeast('ACAACTTTTTAATATATTGTGTATACTCTAA', '19.354839').should.be.false();
  decimal.gcAtLeast('ACAACTTTTTAATATATTGTGTATACTCTAA', '19.354838').should.be.true();
  decimal.gcAtLeast('ACAACTTTTTAATATATTGTGTATACTCTAA', 15).should.be.true();
});

// ---- §4.17-§4.18 ------------------------------------------------------------------------------------------------------

test('proposedCheckRequest deep-equals the §2.9 check block and stops at 5 sets, 10 pairs and 13 primers', function () {
  R.rs871475760_kasp.check.should.eql(EXPECTED_29.check);
  const variant = R.rs871475760_kasp.variant;
  const caps = { system_name: 'sorghum_bicolor', variant: variant, max_sets: 5, max_pairs: 10, max_unique_primers: 13 };
  // a distinct 20-mer per counter value (base 4 over 16 nt)
  const primer = function (x) {
    let s = '';
    for (let i = 0; i < 16; i++, x = Math.floor(x / 4)) s += 'ACGT'[x % 4];
    return s + 'GGCC';
  };
  let n = 1;
  const fake = function (id, shareCommon) {
    const common = shareCommon ? primer(0) : primer(n++);
    const s = { id: id, orientation: 'reverse', products: { ref: { genomic: { region: '1', start: 1, end: 100 } } },
      primers: { as_ref: { target_seq: primer(n++) }, as_alt: { target_seq: primer(n++) }, common: { target_seq: common } } };
    s.check = sets.checkFragment(s);
    return s;
  };
  // 6 sets sharing one common primer: 1 + 2 x 5 = 11 primers, but only 5 sets and 10 pairs fit
  const shared = ['S1', 'S2', 'S3', 'S4', 'S5', 'S6'].map(function (id) { return fake(id, true); });
  const a = sets.proposedCheckRequest(shared, caps);
  a.set_ids.should.eql(['S1', 'S2', 'S3', 'S4', 'S5']);
  a.omitted_set_ids.should.eql(['S6']);
  a.request.pairs.should.have.length(10);
  a.unique_primers.should.equal(11);
  // distinct primers: 4 sets = 12 primers; the fifth would make 15 > 13, and every later set is omitted with it
  const distinct = ['S1', 'S2', 'S3', 'S4', 'S5', 'S6'].map(function (id) { return fake(id, false); });
  const b = sets.proposedCheckRequest(distinct, caps);
  b.set_ids.should.eql(['S1', 'S2', 'S3', 'S4']);
  b.omitted_set_ids.should.eql(['S5', 'S6']);
  b.unique_primers.should.equal(12);
  b.request.genotyping.sets.map(function (s) { return s.id; }).should.eql(b.set_ids);
  should(sets.proposedCheckRequest([], caps)).be.null();
  sets.proposedCheckRequest(distinct, Object.assign({}, caps, { max_pairs: 4 })).set_ids.should.eql(['S1', 'S2']);
});

test('orderRows, oligo names and labels (§2.1, §4.17)', function () {
  const res = R.rs871475760_kasp;
  res.sets[0].order.map(function (r) { return r.name; }).should.eql(['rs871475760_S1_REF_FAM', 'rs871475760_S1_ALT_HEX', 'rs871475760_S1_COM']);
  res.sets[1].order.map(function (r) { return r.notes; }).should.eql(['', 'tail hairpin 66.2 °C', '']);
  R.rs871475760_as_pcr.sets[0].order.map(function (r) { return [r.name, r.dye, r.notes]; }).should.eql([
    ['rs871475760_S1_REF', null, 'deliberate mismatch −2 A→G (Little 1995 Table 9.8.1)'],
    ['rs871475760_S1_ALT', null, 'deliberate mismatch −2 A→G (Little 1995 Table 9.8.1)'],
    ['rs871475760_S1_COM', null, '']
  ]);
  order.sanitizeLabel('1:11109:C:A').should.equal('1_11109_C_A');
  order.labelFor(null, { ids: [], key: '1:11109:C:A' }).should.equal('1_11109_C_A');
  order.labelFor(null, { ids: ['tmp_1_13549_TTA_T,*'], key: 'k' }).should.equal('tmp_1_13549_TTA_T__');
  order.labelFor('my.assay-1', { ids: ['rs1'], key: 'k' }).should.equal('my.assay-1');
  order.labelFor(null, { ids: ['x'.repeat(60)], key: 'k' }).should.have.length(40);
  order.oligoName('L', 'S3', 'as_alt', 'HEX').should.equal('L_S3_ALT_HEX');
  order.oligoName('L', 'S3', 'common', 'HEX').should.equal('L_S3_COM');
  should(order.kaspMix('as_pcr')).be.null();
  order.kaspMix('kasp').should.eql(EXPECTED_29.assay.kasp_mix);
  R.rs5413864115_kasp.sets[0].order[0].name.should.equal('rs5413864115_S1_REF_FAM');
});

test('candidatesFromPairs: rank order, the four §4.8 rejections, duplicates shared across levels', function () {
  const res = R.rs871475760_kasp;
  const c = candidate(res, 'reverse', 0);
  const screen = { orientation: 'reverse', as_3p: 401, zone: res.template.features.zone, delta: 0, template_start: res.template.start,
    target_key: res.variant.key, neighbours: [], window: 5, block: true };
  const seen = new Set();
  const moved = JSON.parse(JSON.stringify(c.pair));
  moved.right.start = 400;
  const out = sets.candidatesFromPairs([c.pair, c.pair, moved], screen, seen);
  out.map(function (x) { return x.reject || 'candidate'; }).should.eql(['candidate', 'duplicate', 'force']);
  sets.candidatesFromPairs([c.pair], screen, seen).map(function (x) { return x.reject; }).should.eql(['duplicate']);
  sets.candidatesFromPairs([c.pair], screen, new Set())[0].candidate.pair_index.should.equal(0);
});
