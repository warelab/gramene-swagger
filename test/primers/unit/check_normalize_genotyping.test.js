'use strict';

// check/normalize.js with a genotyping block (genotyping spec §2.11, §5.1-§5.3, §7.3): the §2.11 body and the example sets of
// the other variants are accepted, request.genotyping keeps only the client's keys, resolved.genotyping holds check/genotype.js
// prepare(), every GENOTYPING_SET_INVALID reason and variant error fires, and genotyping jobs hash with algorithm '2+g1'.
// Offline: the sorghum_v11 catalog fixture, stub assemblies and a sequence stub over the recorded window sorghum_bicolor
// 1:9000-15500 (test/primers/fixtures/check_core/genotype/stubs.js).

const test = require('node:test');
const should = require('should');

const stubs = require('../fixtures/check_core/genotype/stubs');
const { normalize, validateShape, BODY_KEYS } = require('../../../api/helpers/primers/check/normalize');
const genotype = require('../../../api/helpers/primers/check/genotype');
const variation = require('../../../api/helpers/primers/variation/normalize');
const check = require('../../../api/helpers/primers/check');
const jobs = require('../../../api/helpers/primers/jobs');
const { createMemoryStore } = require('../../../api/helpers/primers/jobs/memory_store');
const { PrimerHttpError } = require('../../../api/helpers/primers/errors');

const clone = stubs.clone;
const B = stubs.BODY_2_11;
const CHR1_LENGTH = 80884392;

function mutate(body, fn) {
  const b = clone(body);
  fn(b);
  return b;
}

// The error normalize() rejects with; the test fails if it resolves.
async function refusal(body, deps) {
  try {
    await normalize(clone(body), deps || stubs.checkDeps());
  } catch (e) {
    return e;
  }
  throw new Error('normalize accepted ' + JSON.stringify(body).slice(0, 300));
}

// deps whose catalog lookup and every later stage throw, so a refusal proves that the rule ran before them.
function pureDeps() {
  const unreachable = function () { throw Object.assign(new Error('normalize went past the pure genotyping rules'), { code: 'UNREACHABLE' }); };
  return {
    cfg: stubs.makeCfg(),
    genomes: { getCatalog: unreachable },
    assemblies: { resolve: unreachable },
    sequence: { regionLength: unreachable, fetch: unreachable },
    mongo: {},
    log: stubs.quiet
  };
}

// ---- acceptance ------------------------------------------------------------------------------------------------------------

test('§2.11 body: accepted; request.genotyping is the client block with no added key, from one FASTA read, estimate 95 CPU-s', async function () {
  const deps = stubs.checkDeps();
  const norm = await normalize(clone(B), deps);
  Object.keys(norm.request).should.eql(['system_name', 'mode', 'checks', 'genomes', 'params', 'pairs', 'genotyping']);
  norm.request.genotyping.should.eql(B.genotyping);
  norm.request.pairs.should.eql(B.pairs);
  JSON.stringify(norm.request).should.not.match(/orientation|algorithm|deliberate|shift/);
  norm.request.genomes.should.eql(B.genomes);
  norm.kind.should.equal('pangenome');
  norm.estimate.should.eql({ cpu_s: 95, total: 4 });
  norm.warnings.should.eql([]);
  deps.sequence.calls.should.eql([['1', 9409, 12810]]); // [position - 1,700, position + len(ref) + 1,700]
  // the stored request still passes the handler's own PrimerCheckRequest shape rules
  (function () { validateShape(clone(norm.request), stubs.makeCfg().check); }).should.not.throw();
});

test('§2.11 body: resolved.genotyping holds the derived orientations, the §5.6 core and haplotypes, and each set\'s R_s segment', async function () {
  const g = (await normalize(clone(B), stubs.checkDeps())).resolved.genotyping;
  JSON.parse(JSON.stringify(g)).should.eql(g); // plain JSON for the job doc
  g.algorithm_version.should.equal('g1');
  g.window.should.eql({ region: '1', start: 9409, end: 12810 });
  g.variant.should.match({ key: '1:11109:C:A', region: '1', position: 11109, ref: 'C', alt: 'A', region_length: CHR1_LENGTH,
    shift: 0, tract: null, zone: { start: 11109, end: 11109 }, flank: 15 });
  g.variant.discriminating.should.eql({
    forward: { position: 11109, ref_base: 'C', alt_base: 'A', alt_maps_to: 11109 },
    reverse: { position: 11109, ref_base: 'C', alt_base: 'A', alt_maps_to: 11109 }
  });
  g.variant.core.should.eql({ start: 11108, end: 11110, ref: 'TCT', alt: 'TAT' });
  g.variant.haplotypes.should.eql({ start: 11093, end: 11125, ref: 'ATAGTCATACTCTATTCTGAATTTCTCGCTAGT', alt: 'ATAGTCATACTCTATTATGAATTTCTCGCTAGT' });
  g.sets.should.have.length(2);
  g.sets[0].should.eql({
    id: 'S1', ref_pair: 'S1_REF', alt_pair: 'S1_ALT', orientation: 'reverse',
    primers: { as_ref: 'ATCTTTGACTAGCGAGAAATTCAG', as_alt: 'ATCTTTGACTAGCGAGAAATTCAT', common: 'AGCTTCTCTAAGTGGTTATCCGA' },
    expected: { region: '1', start: 11068, end: 11132 },
    footprints: { as_ref: { start: 11109, end: 11132, strand: -1 }, common: { start: 11068, end: 11090, strand: 1 } },
    deliberate_mismatch: { as_ref: null, as_alt: null },
    deliberate_mismatch_positions: [],
    in_shift_tract: { as_ref: false, as_alt: false },
    segment: { start: 11018, end: 11182, sequence: stubs.bases(11018, 11182) }
  });
  g.sets[1].should.match({
    id: 'S2', orientation: 'forward', deliberate_mismatch_positions: [], in_shift_tract: { as_ref: false, as_alt: false },
    primers: { as_ref: 'GGTTATCCGAATATAGTCATACTCTATTC', as_alt: 'GGTTATCCGAATATAGTCATACTCTATTA', common: 'TCTTTGTCTACTGAGAAATCCAGA' },
    footprints: { as_ref: { start: 11081, end: 11109, strand: 1 }, common: { start: 11149, end: 11172, strand: -1 } },
    segment: { start: 11031, end: 11222 }
  });
});

test('without genotyping the same body normalizes exactly as before: no request or resolved key, same estimate and dbs', async function () {
  const withG = await normalize(clone(B), stubs.checkDeps());
  const plainBody = mutate(B, function (b) { delete b.genotyping; });
  const deps = stubs.checkDeps();
  const plain = await normalize(plainBody, deps);
  deps.sequence.calls.should.eql([]);
  should(plain.request).not.have.property('genotyping');
  Object.keys(plain.resolved).should.eql(['assemblies', 'gene', 'species']);
  const expected = clone(withG);
  delete expected.request.genotyping;
  delete expected.resolved.genotyping;
  plain.should.eql(expected); // §5.4: 94.16 and 94.96 CPU-s both round up to 95
});

test('a lowercase variant that is not left-aligned is stored uppercase and left-aligned, and shares the left-aligned body\'s job id', async function () {
  const aligned = await normalize(stubs.body(stubs.VARIANTS.rs5413864115, stubs.SETS.rs5413864115), stubs.checkDeps());
  const shifted = await normalize(stubs.body({ region: '1', position: 11284, ref: 'aa', alt: 'a' }, stubs.SETS.rs5413864115), stubs.checkDeps());
  shifted.request.genotyping.variant.should.eql({ region: '1', position: 11282, ref: 'CA', alt: 'C' });
  shifted.request.should.eql(aligned.request);
  jobs.jobId(shifted.request, shifted.dbs, check.algorithmVersionFor(shifted.request))
    .should.equal(jobs.jobId(aligned.request, aligned.dbs, check.algorithmVersionFor(aligned.request)));
  shifted.resolved.genotyping.window.should.eql({ region: '1', start: 9584, end: 12986 }); // read around the given position

  const g = aligned.resolved.genotyping;
  g.variant.should.match({ key: '1:11282:CA:C', shift: 2, tract: { start: 11283, end: 11285 }, zone: { start: 11282, end: 11286 }, flank: 15 });
  g.variant.core.should.eql({ start: 11281, end: 11286, ref: 'TCAAAG', alt: 'TCAAG' });
  g.variant.haplotypes.should.match({ start: 11266, end: 11301 });
  g.variant.haplotypes.ref.should.have.length(36);
  g.variant.haplotypes.alt.should.have.length(35);
  // §4.12: for a deletion only as_ref can end inside the tract (S1 11285, S2 11283; as_alt 11286 and the anchor 11282)
  g.sets.map(function (s) { return [s.id, s.orientation, s.in_shift_tract, s.footprints.as_ref]; }).should.eql([
    ['S1', 'forward', { as_ref: true, as_alt: false }, { start: 11257, end: 11285, strand: 1 }],
    ['S2', 'reverse', { as_ref: true, as_alt: false }, { start: 11283, end: 11307, strand: -1 }]
  ]);
});

test('the insertion tmp_1_11502_C_CGT and the TGG repeat rs5413863413: §5.6 cores, flanks and haplotypes, §4.12 shift-tract flags', async function () {
  const ins = (await normalize(stubs.body(stubs.VARIANTS.tmp_1_11502_C_CGT, stubs.SETS.tmp_1_11502_C_CGT), stubs.checkDeps())).resolved.genotyping;
  ins.variant.should.match({ key: '1:11502:C:CGT', shift: 0, tract: null, zone: { start: 11502, end: 11503 }, flank: 15 });
  ins.variant.core.should.eql({ start: 11501, end: 11503, ref: 'CCA', alt: 'CCGTA' });
  ins.variant.haplotypes.should.match({ start: 11486, end: 11518 });
  [ins.variant.haplotypes.ref.length, ins.variant.haplotypes.alt.length].should.eql([33, 35]);
  ins.sets.map(function (s) { return [s.orientation, s.in_shift_tract, s.footprints]; }).should.eql([
    ['reverse', { as_ref: false, as_alt: false }, { as_ref: { start: 11502, end: 11529, strand: -1 }, common: { start: 11458, end: 11479, strand: 1 } }]
  ]);

  const rep = (await normalize(stubs.body(stubs.VARIANTS.rs5413863413, stubs.SETS.rs5413863413), stubs.checkDeps())).resolved.genotyping;
  rep.variant.should.match({ key: '1:13735:TTGG:T', shift: 13, tract: { start: 13736, end: 13751 }, zone: { start: 13735, end: 13752 }, flank: 22 });
  rep.variant.discriminating.should.eql({
    forward: { position: 13749, ref_base: 'G', alt_base: 'A', alt_maps_to: 13752 },
    reverse: { position: 13738, ref_base: 'G', alt_base: 'T', alt_maps_to: 13735 }
  });
  rep.variant.core.should.eql({ start: 13734, end: 13752, ref: 'ATTGGTGGTGGTGGTGGTA', alt: 'ATTGGTGGTGGTGGTA' });
  rep.variant.haplotypes.should.match({ start: 13712, end: 13774 });
  [rep.variant.haplotypes.ref.length, rep.variant.haplotypes.alt.length].should.eql([63, 60]);
  // §5.7: both sets have as_ref inside the tract and as_alt outside it
  rep.sets.map(function (s) { return [s.id, s.orientation, s.in_shift_tract]; }).should.eql([
    ['F', 'forward', { as_ref: true, as_alt: false }],
    ['R', 'reverse', { as_ref: true, as_alt: false }]
  ]);
});

test('the §2.10(d) AS-PCR primers, each with an undeclared -2 mismatch, are accepted: deliberate_mismatch_positions [2, 2]', async function () {
  const norm = await normalize(stubs.body(stubs.VARIANTS.rs871475760, stubs.SETS.rs871475760_as_pcr), stubs.checkDeps());
  norm.request.genotyping.sets.should.eql([{ id: 'S1', ref_pair: 'S1_REF', alt_pair: 'S1_ALT' }]);
  norm.resolved.genotyping.sets[0].should.match({
    orientation: 'reverse',
    deliberate_mismatch: { as_ref: 2, as_alt: 2 },
    deliberate_mismatch_positions: [2, 2],
    footprints: { as_ref: { start: 11109, end: 11132, strand: -1 }, common: { start: 10880, end: 10903, strand: 1 } }
  });
});

test('gene mode is accepted like region mode', async function () {
  const doc = {
    _id: 'SORBI_3001G000200', system_name: 'sorghum_bicolor', taxon_id: 4558006,
    location: { region: '1', start: 11180, end: 14899, strand: -1, map: 'GCA_000003195.3' },
    gene_structure: { canonical_transcript: 'SORBI_3001G000200.1', transcripts: [{ id: 'SORBI_3001G000200.1' }] }
  };
  const mongo = {
    genes: {
      mongoCollection: async function () {
        return { find: function (q) { return { toArray: async function () { return q._id === doc._id ? [clone(doc)] : []; } }; } };
      }
    }
  };
  const body = stubs.body(stubs.VARIANTS.tmp_1_11502_C_CGT, stubs.SETS.tmp_1_11502_C_CGT, { mode: 'gene', gene_id: 'SORBI_3001G000200' });
  const norm = await normalize(body, stubs.checkDeps({ mongo: mongo }));
  norm.request.should.match({ mode: 'gene', gene_id: 'SORBI_3001G000200', transcript_id: 'SORBI_3001G000200.1' });
  norm.request.genotyping.should.eql(body.genotyping);
});

// ---- rejections --------------------------------------------------------------------------------------------------------------

test('genotyping in transcript or sequence mode -> 400 INVALID_REQUEST {field: "genotyping", reason: "mode"}, before any catalog lookup', async function () {
  genotype.MODES.should.eql(['gene', 'region']);
  for (const extra of [{ mode: 'transcript', gene_id: 'SORBI_3001G000200' }, { mode: 'sequence' }]) {
    const e = await refusal(Object.assign(clone(B), extra), pureDeps());
    e.should.be.instanceOf(PrimerHttpError);
    e.should.match({ status: 400, code: 'INVALID_REQUEST' });
    e.details.should.eql({ field: 'genotyping', reason: 'mode' });
  }
});

test('step 1 mirrors PrimerCheckGenotyping: 400 INVALID_REQUEST {field}, before any catalog lookup', async function () {
  BODY_KEYS.should.containEql('genotyping');
  (function () { validateShape(clone(B), stubs.makeCfg().check); }).should.not.throw();
  const V = B.genotyping.variant;
  const S = B.genotyping.sets;
  const variant = function (patch) { return { variant: Object.assign({}, V, patch), sets: S }; };
  const set0 = function (patch) { return { variant: V, sets: [Object.assign({}, S[0], patch)] }; };
  const rows = [
    [null, 'genotyping'],
    [[], 'genotyping'],
    [{ variant: V, sets: S, orientation: 'reverse' }, 'genotyping.orientation'],
    [{ sets: S }, 'genotyping.variant'],
    [{ variant: V }, 'genotyping.sets'],
    [variant({ strand: 1 }), 'genotyping.variant.strand'],
    [variant({ region: '' }), 'genotyping.variant.region'],
    [variant({ region: 'r'.repeat(256) }), 'genotyping.variant.region'],
    [variant({ position: 0 }), 'genotyping.variant.position'],
    [variant({ position: 11109.5 }), 'genotyping.variant.position'],
    [variant({ position: '11109' }), 'genotyping.variant.position'],
    [variant({ ref: 'N' }), 'genotyping.variant.ref'],
    [variant({ alt: '-' }), 'genotyping.variant.alt'],
    [variant({ alt: 'A'.repeat(51) }), 'genotyping.variant.alt'],
    [variant({ alt: 'c' }), 'genotyping.variant.alt'], // the same allele as ref
    [{ variant: V, sets: [] }, 'genotyping.sets'],
    [{ variant: V, sets: [0, 1, 2, 3, 4, 5].map(function (i) { return { id: 'X' + i, ref_pair: 'S1_REF', alt_pair: 'S1_ALT' }; }) }, 'genotyping.sets'],
    [{ variant: V, sets: ['S1'] }, 'genotyping.sets[0]'],
    [set0({ orientation: 'reverse' }), 'genotyping.sets[0].orientation'],
    [{ variant: V, sets: [{ id: 'S1', ref_pair: 'S1_REF' }] }, 'genotyping.sets[0].alt_pair'],
    [set0({ id: 'S 1' }), 'genotyping.sets[0].id'],
    [set0({ ref_pair: 'p'.repeat(65) }), 'genotyping.sets[0].ref_pair'],
    [{ variant: V, sets: [S[0], Object.assign({}, S[1], { id: 'S1' })] }, 'genotyping.sets[1].id'] // set ids are unique
  ];
  for (const row of rows) {
    const e = await refusal(Object.assign(clone(B), { genotyping: row[0] }), pureDeps());
    e.should.match({ status: 400, code: 'INVALID_REQUEST' }, JSON.stringify(row[0]));
    e.details.field.should.equal(row[1], JSON.stringify(row[0]));
  }
});

// [case, body, details]; part A runs before any catalog lookup, part B after the reference resolves (one FASTA read).
const PART_A = [
  ['unknown_pair', mutate(B, function (b) { b.genotyping.sets[0].ref_pair = 'NOPE'; }),
    { set_id: 'S1', reason: 'unknown_pair', pair_ids: ['NOPE'] }],
  ['same_pair', mutate(B, function (b) { b.genotyping.sets[0].alt_pair = 'S1_REF'; }),
    { set_id: 'S1', reason: 'same_pair', pair_ids: ['S1_REF'] }],
  ['pair_reused', mutate(B, function (b) { b.genotyping.sets[1].ref_pair = 'S1_REF'; }),
    { set_id: 'S2', reason: 'pair_reused', pair_ids: ['S1_REF'] }],
  ['expected_required (missing on one pair)', mutate(B, function (b) { delete b.pairs[1].expected; }),
    { set_id: 'S1', reason: 'expected_required', pair_ids: ['S1_ALT'] }],
  ['expected_differs', mutate(B, function (b) { b.pairs[1].expected.end = 11133; }),
    { set_id: 'S1', reason: 'expected_differs', pair_ids: ['S1_REF', 'S1_ALT'] }],
  ['no_shared_common (no primer shared)', mutate(B, function (b) { b.pairs[1].left = 'AGCTTCTCTAAGTGGTTATCCGT'; }),
    { set_id: 'S1', reason: 'no_shared_common', pair_ids: ['S1_REF', 'S1_ALT'] }],
  ['no_shared_common (both primers shared)', mutate(B, function (b) { b.pairs[1].right = b.pairs[0].right; }),
    { set_id: 'S1', reason: 'no_shared_common', pair_ids: ['S1_REF', 'S1_ALT'] }]
];

const COMMON_OVER_VARIANT = [['C', 'forward', 'GGTTATCCGAATATAGTCATACTCTATTC', 'GGTTATCCGAATATAGTCATACTCTATTA', 'TCTTTGACTAGCGAGAAATTCAGA', 11081, 11131]];
const SWAPPED_DELETION = mutate(stubs.body(stubs.VARIANTS.rs5413864115, stubs.SETS.rs5413864115), function (b) {
  b.genotyping.sets[1] = { id: 'S2', ref_pair: 'S2_ALT', alt_pair: 'S2_REF' };
});
const PART_B = [
  ['not_at_variant (§2.11: a REF-specific primer 1 nt off)', mutate(B, function (b) { b.pairs[0].right = 'TCTTTGACTAGCGAGAAATTCAGA'; }),
    { set_id: 'S1', reason: 'not_at_variant', pair_ids: ['S1_REF'], primer: 'ref_pair', sequence: 'TCTTTGACTAGCGAGAAATTCAGA' }],
  ['not_at_variant (an ALT-specific primer 1 nt the other way)', mutate(B, function (b) { b.pairs[1].right = 'CATCTTTGACTAGCGAGAAATTCA'; }),
    { set_id: 'S1', reason: 'not_at_variant', pair_ids: ['S1_ALT'], primer: 'alt_pair', sequence: 'CATCTTTGACTAGCGAGAAATTCA' }],
  ['alleles_swapped (SNV)', mutate(B, function (b) { b.genotyping.sets[0] = { id: 'S1', ref_pair: 'S1_ALT', alt_pair: 'S1_REF' }; }),
    { set_id: 'S1', reason: 'alleles_swapped', pair_ids: ['S1_ALT'], primer: 'ref_pair', sequence: 'ATCTTTGACTAGCGAGAAATTCAT' }],
  ['alleles_swapped (reverse deletion set)', SWAPPED_DELETION,
    { set_id: 'S2', reason: 'alleles_swapped', pair_ids: ['S2_ALT'], primer: 'ref_pair', sequence: 'AGAGTCTTTTCAAATTTCACACTTG' }],
  ['too_many_edits (2 edits: -3 and -2)', mutate(B, function (b) { b.pairs[0].right = 'ATCTTTGACTAGCGAGAAATTGGG'; }),
    { set_id: 'S1', reason: 'too_many_edits', pair_ids: ['S1_REF'], primer: 'ref_pair', sequence: 'ATCTTTGACTAGCGAGAAATTGGG', mm_pos: [3, 2] }],
  ['too_many_edits (one mismatch at -5)', mutate(B, function (b) { b.pairs[0].right = 'ATCTTTGACTAGCGAGAAAATCAG'; }),
    { set_id: 'S1', reason: 'too_many_edits', pair_ids: ['S1_REF'], primer: 'ref_pair', sequence: 'ATCTTTGACTAGCGAGAAAATCAG', mm_pos: [5] }],
  ['common_in_zone (a common primer over the variant, 1:11108-11131)', stubs.body(stubs.VARIANTS.rs871475760, COMMON_OVER_VARIANT),
    { set_id: 'C', reason: 'common_in_zone', pair_ids: ['C_REF', 'C_ALT'], primer: 'common', sequence: 'TCTTTGACTAGCGAGAAATTCAGA' }],
  ['expected_mismatch (expected.end 1 bp past the REF-specific primer)', mutate(B, function (b) { b.pairs[0].expected.end = 11133; b.pairs[1].expected.end = 11133; }),
    { set_id: 'S1', reason: 'expected_mismatch', pair_ids: ['S1_REF', 'S1_ALT'], primer: 'ref_pair', sequence: 'ATCTTTGACTAGCGAGAAATTCAG' }],
  ['expected_mismatch (expected on another region)', mutate(B, function (b) { b.pairs[0].expected.region = '2'; b.pairs[1].expected.region = '2'; }),
    { set_id: 'S1', reason: 'expected_mismatch', pair_ids: ['S1_REF', 'S1_ALT'], primer: 'ref_pair', sequence: 'ATCTTTGACTAGCGAGAAATTCAG' }]
];

test('the rejection cases cover every GENOTYPING_SET_INVALID reason of §2.14', function () {
  const reasons = new Set(PART_A.concat(PART_B).map(function (row) { return row[2].reason; }));
  Array.from(reasons).sort().should.eql(genotype.REASONS.slice().sort());
  genotype.REASONS.should.have.length(11);
});

test('§5.2 part A: unknown_pair, same_pair, pair_reused, expected_required, expected_differs, no_shared_common, before any catalog lookup', async function () {
  for (const row of PART_A) {
    const e = await refusal(row[1], pureDeps());
    e.should.match({ status: 400, code: 'GENOTYPING_SET_INVALID' }, row[0] + ': ' + e.message);
    e.details.should.eql(row[2], row[0]);
    e.message.should.match(/^set S\d: /, row[0]);
  }
  (await refusal(PART_A[5][1], pureDeps())).message
    .should.equal('set S1: S1_REF and S1_ALT must share exactly one primer, the common primer, on the same side'); // §2.11
  // expected is required on both pairs of a set, whichever lacks it
  (await refusal(mutate(B, function (b) { delete b.pairs[2].expected; delete b.pairs[3].expected; }), pureDeps())).details
    .should.eql({ set_id: 'S2', reason: 'expected_required', pair_ids: ['S2_REF', 'S2_ALT'] });
});

test('§5.2 part B: not_at_variant (§2.11), alleles_swapped, too_many_edits, common_in_zone, expected_mismatch, each from one FASTA read', async function () {
  for (const row of PART_B) {
    const deps = stubs.checkDeps();
    const e = await refusal(row[1], deps);
    e.should.match({ status: 400, code: 'GENOTYPING_SET_INVALID' }, row[0] + ': ' + e.message);
    e.details.should.eql(row[2], row[0]);
    deps.sequence.calls.should.have.length(1, row[0]);
  }
  (await refusal(PART_B[0][1])).message
    .should.equal('set S1: the REF-specific primer TCTTTGACTAGCGAGAAATTCAGA does not end at the variant (1:11109)'); // §2.11
});

test('variant errors: REF_MISMATCH, VARIANT_TOO_REPETITIVE, UNKNOWN_REGION, REGION_OUT_OF_BOUNDS, NO_SEQUENCE', async function () {
  const refMismatch = await refusal(mutate(B, function (b) { b.genotyping.variant.ref = 'A'; b.genotyping.variant.alt = 'C'; }));
  refMismatch.should.match({ status: 400, code: 'REF_MISMATCH', message: 'the reference allele A does not match the genome base C at 1:11109' });
  refMismatch.details.should.eql({ region: '1', position: 11109, given: 'A', genome: 'C' });

  const repetitive = await refusal(stubs.body(stubs.VARIANTS.rs5413864115, stubs.SETS.rs5413864115),
    stubs.checkDeps({ cfg: stubs.makeCfg({ variation: { max_shift: 1 } }) }));
  repetitive.should.match({ status: 400, code: 'VARIANT_TOO_REPETITIVE' });
  repetitive.details.should.eql({ region: '1', position: 11282, shift: 2, max: 1 });

  const unknownRegion = await refusal(mutate(B, function (b) { b.genotyping.variant.region = 'Z'; }));
  unknownRegion.should.match({ status: 404, code: 'UNKNOWN_REGION' });
  unknownRegion.details.should.eql({ region: 'Z' });

  const beyond = await refusal(mutate(B, function (b) { b.genotyping.variant.position = CHR1_LENGTH; b.genotyping.variant.ref = 'CA'; }));
  beyond.should.match({ status: 400, code: 'REGION_OUT_OF_BOUNDS' });
  beyond.details.should.eql({ region: '1', start: CHR1_LENGTH, end: CHR1_LENGTH + 1, length: CHR1_LENGTH });

  const noFasta = await refusal(B, stubs.checkDeps({ assemblies: stubs.stubAssemblies({ sorghum_bicolor: { fasta: { dna: null, dna_sm: null } } }) }));
  noFasta.should.match({ status: 422, code: 'NO_SEQUENCE' });
  noFasta.details.should.eql({ system_name: 'sorghum_bicolor' });

  const unreadable = new PrimerHttpError(422, 'NO_SEQUENCE', 'sequence files are not available for this genome', {});
  (await refusal(B, stubs.checkDeps({ sequence: stubs.sequenceStub({ fail: unreadable }) }))).should.equal(unreadable);
});

test('prepare maps reads past the region to REGION_OUT_OF_BOUNDS and a slide past the window to VARIANT_TOO_REPETITIVE', function () {
  const links = genotype.validateLinks(B.genotyping, B.pairs);
  const throwsFrom = function (variant, genome) {
    try {
      genotype.prepare({ variant: variant, sets: B.genotyping.sets }, B.pairs, genome, { cfg: stubs.makeCfg(), links: links });
    } catch (e) {
      return e;
    }
    throw new Error('prepare accepted ' + JSON.stringify(variant));
  };
  // left-aligning AA>A at position 2 of a poly-A region start needs base 0
  const start = throwsFrom({ region: '1', position: 2, ref: 'AA', alt: 'A' }, variation.sequenceWindow('AAAACGTACGTACGT', 1, 15));
  start.should.match({ status: 400, code: 'REGION_OUT_OF_BOUNDS' });
  start.details.should.eql({ region: '1', start: 2, end: 3, length: 15 });
  // a window that starts at the given position of a right-shifted deletion cannot left-align it
  const win = throwsFrom({ region: '1', position: 11284, ref: 'AA', alt: 'A' }, variation.sequenceWindow(stubs.bases(11284, 13000), 11284, CHR1_LENGTH));
  win.should.match({ status: 400, code: 'VARIANT_TOO_REPETITIVE' });
  win.details.should.eql({ region: '1', position: 11284, shift: null, max: 1000 });
  genotype.windowFor({ position: 11109, ref: 'C' }, CHR1_LENGTH).should.eql({ start: 9409, end: 12810 });
  genotype.windowFor({ position: 100, ref: 'CA' }, 1000).should.eql({ start: 1, end: 1000 });
});

// ---- job ids ---------------------------------------------------------------------------------------------------------------

test('job ids: a genotyping request hashes with algorithm 2+g1 through jobs.submit; without genotyping it keeps algorithm 2', async function () {
  check.ALGORITHM_VERSION.should.equal('2');
  genotype.GENOTYPING_VERSION.should.equal('g1');
  check.algorithmVersionFor({ system_name: 'sorghum_bicolor', pairs: [] }).should.equal('2');
  check.algorithmVersionFor(undefined).should.equal('2');

  const norm = await normalize(clone(B), stubs.checkDeps());
  check.algorithmVersionFor(norm.request).should.equal('2+g1');
  const store = createMemoryStore({ cfg: stubs.makeCfg(), siteKey: 'check_genotyping_test', shared: { slots: new Map(), panSlots: new Map() } });
  const submitted = await jobs.submit(clone(B), stubs.checkDeps({ store: store, check: check }));
  submitted.job_id.should.equal(jobs.jobId(norm.request, norm.dbs, '2+g1'));
  submitted.job_id.should.not.equal(jobs.jobId(norm.request, norm.dbs, '2'));
  submitted.should.match({ status: 'queued', kind: 'pangenome', estimate: { cpu_s: 95 }, created: true });
  const doc = await store.getJob(submitted.job_id);
  doc.request.should.eql(norm.request);
  doc.resolved.genotyping.sets.map(function (s) { return s.orientation; }).should.eql(['reverse', 'forward']);
  // the status body returns request.genotyping as sent, never resolved
  const status = await jobs.status(submitted.job_id, stubs.checkDeps({ store: store }));
  status.request.genotyping.should.eql(B.genotyping);
  should(status).not.have.property('resolved');

  const plainBody = mutate(B, function (b) { delete b.genotyping; });
  const plain = await normalize(clone(plainBody), stubs.checkDeps());
  const plainStore = createMemoryStore({ cfg: stubs.makeCfg(), siteKey: 'check_genotyping_plain', shared: { slots: new Map(), panSlots: new Map() } });
  (await jobs.submit(plainBody, stubs.checkDeps({ store: plainStore, check: check }))).job_id
    .should.equal(jobs.jobId(plain.request, plain.dbs, check.ALGORITHM_VERSION));
});
