'use strict';

require('../../../api/helpers/primers/node_compat');

const test = require('node:test');
const should = require('should');
const fs = require('fs');

const design = require('../../../api/helpers/primers/design');
const coords = require('../../../api/helpers/primers/coords');
const template = require('../../../api/helpers/primers/template');
const boulder = require('../../../api/helpers/primers/boulder');
const { revcomp } = require('../../../api/helpers/primers/sequence');
const { PrimerHttpError } = require('../../../api/helpers/primers/errors');
const stubs = require('../fixtures/design/stubs');

const CFG = stubs.cfg();
const G200 = 'SORBI_3001G000200';

function norm(body, cfg) {
  return design.normalize(body, cfg || CFG);
}

function throwsCode(fn, code, status) {
  let err = null;
  try {
    fn();
  } catch (e) {
    err = e;
  }
  should.exist(err, 'expected ' + code);
  err.should.have.property('code', code);
  if (status) err.should.have.property('status', status);
  return err;
}

async function rejectsCode(promise, code, status) {
  let err = null;
  try {
    await promise;
  } catch (e) {
    err = e;
  }
  should.exist(err, 'expected ' + code);
  err.should.have.property('code', code);
  if (status) err.should.have.property('status', status);
  return err;
}

const REGION = { region: '1', start: 11080, end: 15099, strand: -1 };

// ---- normalize: request shape --------------------------------------------------------------------

test('normalize: the body must be an object with a known mode and no unknown fields', function () {
  throwsCode(function () { norm(null); }, 'INVALID_REQUEST', 400);
  throwsCode(function () { norm([]); }, 'INVALID_REQUEST', 400);
  throwsCode(function () { norm({}); }, 'INVALID_REQUEST', 400).details.field.should.equal('mode');
  throwsCode(function () { norm({ mode: 'protein' }); }, 'INVALID_REQUEST', 400);
  throwsCode(function () { norm({ mode: 'gene', gene_id: 'X', bogus: 1 }); }, 'INVALID_REQUEST', 400).details.field.should.equal('bogus');
});

test('normalize: mode-required fields', function () {
  throwsCode(function () { norm({ mode: 'gene' }); }, 'INVALID_REQUEST').details.field.should.equal('gene_id');
  throwsCode(function () { norm({ mode: 'transcript' }); }, 'INVALID_REQUEST').details.field.should.equal('gene_id');
  throwsCode(function () { norm({ mode: 'region', region: REGION }); }, 'INVALID_REQUEST').details.field.should.equal('system_name');
  throwsCode(function () { norm({ mode: 'region', system_name: 'sorghum_bicolor' }); }, 'INVALID_REQUEST').details.field.should.equal('region');
  throwsCode(function () { norm({ mode: 'sequence' }); }, 'INVALID_REQUEST').details.field.should.equal('sequence');
  norm({ mode: 'gene', gene_id: G200 }).gene_id.should.equal(G200);
  norm({ mode: 'transcript', gene_id: G200, transcript_id: 'SORBI_3001G000200.1' }).transcript_id.should.equal('SORBI_3001G000200.1');
  norm({ mode: 'region', system_name: 'sorghum_bicolor', region: REGION }).region.should.eql(REGION);
  norm({ mode: 'sequence', sequence: '>x\nACGT' }).sequence.should.equal('>x\nACGT');
});

test('normalize: field types - gene_id must be a string (no operator objects reach mongo), patterns, ranges, enums, intervals', function () {
  throwsCode(function () { norm({ mode: 'gene', gene_id: { $gt: '' } }); }, 'INVALID_REQUEST').details.field.should.equal('gene_id');
  throwsCode(function () { norm({ mode: 'gene', gene_id: 'x'.repeat(256) }); }, 'INVALID_REQUEST');
  throwsCode(function () { norm({ mode: 'transcript', gene_id: 'G', transcript_id: 7 }); }, 'INVALID_REQUEST');
  throwsCode(function () { norm({ mode: 'region', system_name: '../etc', region: REGION }); }, 'INVALID_REQUEST').details.field.should.equal('system_name');
  throwsCode(function () { norm({ mode: 'region', system_name: 'Sorghum', region: REGION }); }, 'INVALID_REQUEST');
  throwsCode(function () { norm({ mode: 'gene', gene_id: 'G', flank_up: 10001 }); }, 'INVALID_REQUEST').details.field.should.equal('flank_up');
  throwsCode(function () { norm({ mode: 'gene', gene_id: 'G', flank_down: -1 }); }, 'INVALID_REQUEST');
  throwsCode(function () { norm({ mode: 'gene', gene_id: 'G', flank_down: 1.5 }); }, 'INVALID_REQUEST');
  throwsCode(function () { norm({ mode: 'gene', gene_id: 'G', avoid_repeats: 'yes' }); }, 'INVALID_REQUEST');
  throwsCode(function () { norm({ mode: 'gene', gene_id: 'G', template_only: 1 }); }, 'INVALID_REQUEST');
  throwsCode(function () { norm({ mode: 'gene', gene_id: 'G', repeat_mask_mode: 'hard' }); }, 'INVALID_REQUEST');
  throwsCode(function () { norm({ mode: 'gene', gene_id: 'G', target: [0, 10] }); }, 'INVALID_REQUEST').details.field.should.equal('target');
  throwsCode(function () { norm({ mode: 'gene', gene_id: 'G', included: [5] }); }, 'INVALID_REQUEST');
  const many = [];
  for (let i = 0; i < 51; i++) many.push([i + 1, 1]);
  throwsCode(function () { norm({ mode: 'gene', gene_id: 'G', excluded: many }); }, 'INVALID_REQUEST');
  throwsCode(function () { norm({ mode: 'gene', gene_id: 'G', excluded: [[1, 1], [2, 0]] }); }, 'INVALID_REQUEST').details.field.should.equal('excluded[1]');
  throwsCode(function () { norm({ mode: 'region', system_name: 's', region: { region: '1', start: 1, end: 100, strand: 0 } }); }, 'INVALID_REQUEST');
  throwsCode(function () { norm({ mode: 'region', system_name: 's', region: { chr: '1', start: 1, end: 100 } }); }, 'INVALID_REQUEST');
  throwsCode(function () { norm({ mode: 'region', system_name: 's', region: { region: '1', start: 0, end: 100 } }); }, 'INVALID_REQUEST');
  throwsCode(function () { norm({ mode: 'sequence', sequence: 'A'.repeat(60001) }); }, 'INVALID_REQUEST');
  norm({ mode: 'region', system_name: 's', region: { region: '1', start: 1, end: 100 } }).region.strand.should.equal(1);
  norm({ mode: 'gene', gene_id: 'G', excluded: [[1, 1], [5, 2]] }).excluded.should.eql([[1, 1], [5, 2]]);
});

test('normalize: region start > end is REGION_OUT_OF_BOUNDS and an oversized region TEMPLATE_TOO_LONG, before any I/O', function () {
  throwsCode(function () { norm({ mode: 'region', system_name: 's', region: { region: '1', start: 100, end: 99 } }); }, 'REGION_OUT_OF_BOUNDS', 400);
  const err = throwsCode(function () { norm({ mode: 'region', system_name: 's', region: { region: '1', start: 1, end: 50001 } }); }, 'TEMPLATE_TOO_LONG', 400);
  err.details.should.eql({ length: 50001, max: 50000 });
  norm({ mode: 'region', system_name: 's', region: { region: '1', start: 1, end: 50000 } }).region.end.should.equal(50000);
});

// ---- presets and params ---------------------------------------------------------------------------

test('presets: qpcr for transcripts (junction_spanning defaults to true), pcr otherwise; client params override preset values', function () {
  const t = norm({ mode: 'transcript', gene_id: 'G' });
  t.preset.should.equal('qpcr');
  t.params.should.eql(design.PRESETS.qpcr);
  t.params.should.not.equal(design.PRESETS.qpcr);
  [t.junction_spanning, t.avoid_repeats, t.repeat_mask_mode, t.template_only, t.flank_up].should.eql([true, false, null, false, 0]);
  norm({ mode: 'transcript', gene_id: 'G', junction_spanning: false }).junction_spanning.should.be.false();
  norm({ mode: 'transcript', gene_id: 'G', flank_up: 500 }).flank_up.should.equal(0);

  const g = norm({ mode: 'gene', gene_id: 'G', flank_up: 200, params: { max_size: 30, product_size_ranges: [[300, 800]] } });
  g.preset.should.equal('pcr');
  g.junction_spanning.should.be.false();
  g.flank_up.should.equal(200);
  g.params.should.eql(Object.assign({}, design.PRESETS.pcr, { max_size: 30, product_size_ranges: [[300, 800]] }));
  design.PRESETS.pcr.max_size.should.equal(25);
  Object.isFrozen(design.PRESETS.pcr.product_size_ranges).should.be.true();
  norm({ mode: 'region', system_name: 's', region: REGION }).preset.should.equal('pcr');
  norm({ mode: 'sequence', sequence: 'ACGT' }).preset.should.equal('pcr');

  design.PRESETS.pcr.should.eql({ opt_size: 20, min_size: 18, max_size: 25, opt_tm: 60, min_tm: 57, max_tm: 63, min_gc: 30, max_gc: 70,
    max_tm_diff: 3, max_poly_x: 4, product_size_ranges: [[100, 1000]], num_return: 5 });
  design.PRESETS.qpcr.should.eql({ opt_size: 20, min_size: 18, max_size: 24, opt_tm: 60, min_tm: 58, max_tm: 62, min_gc: 35, max_gc: 65,
    max_tm_diff: 2, max_poly_x: 4, product_size_ranges: [[70, 150]], num_return: 5,
    min_3_prime_overlap_of_junction: 4, min_5_prime_overlap_of_junction: 7 });
});

test('avoid_repeats: n_mask (default) forces max_ns 0 in the effective params; three_prime keeps it', function () {
  const n = norm({ mode: 'region', system_name: 's', region: REGION, avoid_repeats: true, params: { max_ns: 3 } });
  n.repeat_mask_mode.should.equal('n_mask');
  n.params.max_ns.should.equal(0);
  const p = norm({ mode: 'region', system_name: 's', region: REGION, avoid_repeats: true, repeat_mask_mode: 'three_prime', params: { max_ns: 3 } });
  p.repeat_mask_mode.should.equal('three_prime');
  p.params.max_ns.should.equal(3);
  should(norm({ mode: 'region', system_name: 's', region: REGION, repeat_mask_mode: 'three_prime' }).repeat_mask_mode).be.null();
});

test('params: closed set with type and range checks -> INVALID_PARAMS', function () {
  const bad = [
    { primer_opt_size: 20 }, { min_size: 14 }, { max_size: 37 }, { max_size: 20.5 }, { num_return: 21 }, { num_return: 0 },
    { opt_tm: '60' }, { opt_tm: 29.9 }, { max_gc: 101 }, { max_tm_diff: 31 }, { max_poly_x: 11 }, { gc_clamp: 6 },
    { max_ns: 6 }, { salt_monovalent: 1001 }, { dna_conc: -1 }, { min_3_prime_overlap_of_junction: 21 },
    { min_5_prime_overlap_of_junction: 0 }, { product_size_ranges: [] }, { product_size_ranges: [[30]] },
    { product_size_ranges: [[10, 100]] }, { product_size_ranges: [[100, 50001]] }, { product_size_ranges: 'x' }
  ];
  bad.forEach(function (params) {
    throwsCode(function () { norm({ mode: 'gene', gene_id: 'G', params: params }); }, 'INVALID_PARAMS', 400);
  });
  const eleven = [];
  for (let i = 0; i < 11; i++) eleven.push([100 + i, 200 + i]);
  throwsCode(function () { norm({ mode: 'gene', gene_id: 'G', params: { product_size_ranges: eleven } }); }, 'INVALID_PARAMS');
  throwsCode(function () { norm({ mode: 'gene', gene_id: 'G', params: [] }); }, 'INVALID_PARAMS');
  throwsCode(function () { norm({ mode: 'gene', gene_id: 'G', params: { __proto__: { x: 1 }, bogus: 1 } }); }, 'INVALID_PARAMS');
  const ok = norm({ mode: 'gene', gene_id: 'G', params: { opt_gc: 50, gc_clamp: 1, max_end_stability: 9, salt_divalent: 1.5, dntp_conc: 0.6 } });
  ok.params.should.have.properties({ opt_gc: 50, gc_clamp: 1, max_end_stability: 9, salt_divalent: 1.5, dntp_conc: 0.6 });
});

test('cross-field: min <= opt <= max for size, Tm and GC', function () {
  throwsCode(function () { norm({ mode: 'gene', gene_id: 'G', params: { min_size: 26 } }); }, 'INVALID_PARAMS').details.param.should.equal('min_size');
  throwsCode(function () { norm({ mode: 'gene', gene_id: 'G', params: { opt_size: 26 } }); }, 'INVALID_PARAMS').details.param.should.equal('opt_size');
  throwsCode(function () { norm({ mode: 'gene', gene_id: 'G', params: { opt_tm: 56 } }); }, 'INVALID_PARAMS').details.param.should.equal('opt_tm');
  throwsCode(function () { norm({ mode: 'gene', gene_id: 'G', params: { min_tm: 64 } }); }, 'INVALID_PARAMS').details.param.should.equal('min_tm');
  throwsCode(function () { norm({ mode: 'gene', gene_id: 'G', params: { min_gc: 80 } }); }, 'INVALID_PARAMS').details.param.should.equal('min_gc');
  throwsCode(function () { norm({ mode: 'gene', gene_id: 'G', params: { opt_gc: 20 } }); }, 'INVALID_PARAMS').details.param.should.equal('opt_gc');
  norm({ mode: 'gene', gene_id: 'G', params: { min_size: 15, opt_size: 15, max_size: 15 } }).params.max_size.should.equal(15);
});

test('cross-field: product ranges need min < max, and max_size may not exceed the smallest product size', function () {
  throwsCode(function () { norm({ mode: 'gene', gene_id: 'G', params: { product_size_ranges: [[300, 300]] } }); }, 'INVALID_PARAMS');
  throwsCode(function () { norm({ mode: 'gene', gene_id: 'G', params: { product_size_ranges: [[100, 200], [500, 400]] } }); }, 'INVALID_PARAMS')
    .details.index.should.equal(1);
  throwsCode(function () { norm({ mode: 'gene', gene_id: 'G', params: { product_size_ranges: [[20, 100]] } }); }, 'INVALID_PARAMS')
    .details.param.should.equal('max_size');
  norm({ mode: 'gene', gene_id: 'G', params: { max_size: 36, product_size_ranges: [[36, 100], [200, 300]] } }).params.max_size.should.equal(36);
});

test('junction overlaps must be <= floor(max_size / 2) when primers must span a junction (Primer3 aborts globally otherwise)', function () {
  const err = throwsCode(function () {
    norm({ mode: 'transcript', gene_id: G200, params: { max_size: 24, min_5_prime_overlap_of_junction: 13 } });
  }, 'INVALID_PARAMS', 400);
  err.details.should.eql({ param: 'min_5_prime_overlap_of_junction', max: 12 });
  norm({ mode: 'transcript', gene_id: G200, params: { max_size: 24, min_5_prime_overlap_of_junction: 12 } }).params.min_5_prime_overlap_of_junction.should.equal(12);
  norm({ mode: 'transcript', gene_id: G200, params: { max_size: 25, min_5_prime_overlap_of_junction: 12 } });
  throwsCode(function () { norm({ mode: 'transcript', gene_id: G200, params: { max_size: 25, min_5_prime_overlap_of_junction: 13 } }); }, 'INVALID_PARAMS');
  throwsCode(function () { norm({ mode: 'transcript', gene_id: G200, params: { max_size: 24, min_3_prime_overlap_of_junction: 13 } }); }, 'INVALID_PARAMS')
    .details.param.should.equal('min_3_prime_overlap_of_junction');
  // no junction list is sent: not Primer3's global error
  norm({ mode: 'transcript', gene_id: G200, junction_spanning: false, params: { max_size: 24, min_5_prime_overlap_of_junction: 13 } });
  norm({ mode: 'gene', gene_id: G200, params: { min_5_prime_overlap_of_junction: 20 } });
  // the qpcr defaults (7/4) fit every allowed max_size (>= 15)
  norm({ mode: 'transcript', gene_id: G200, params: { min_size: 15, opt_size: 15, max_size: 15 } });
});

// ---- template-dependent checks ---------------------------------------------------------------------

test('validateAgainstTemplate: target, included and excluded inside the template; a product range that can fit', function () {
  const t = { length: 1000 };
  const e1 = throwsCode(function () { design.validateAgainstTemplate(norm({ mode: 'sequence', sequence: 'A', target: [900, 102] }), t); },
    'INTERVAL_OUT_OF_BOUNDS', 400);
  e1.details.should.eql({ field: 'target', interval: [900, 102], template_length: 1000 });
  design.validateAgainstTemplate(norm({ mode: 'sequence', sequence: 'A', target: [901, 100] }), t);
  throwsCode(function () { design.validateAgainstTemplate(norm({ mode: 'sequence', sequence: 'A', included: [1001, 1] }), t); }, 'INTERVAL_OUT_OF_BOUNDS');
  throwsCode(function () { design.validateAgainstTemplate(norm({ mode: 'sequence', sequence: 'A', excluded: [[1, 10], [995, 10]] }), t); },
    'INTERVAL_OUT_OF_BOUNDS').details.field.should.equal('excluded[1]');
  throwsCode(function () {
    design.validateAgainstTemplate(norm({ mode: 'sequence', sequence: 'A', params: { product_size_ranges: [[1001, 2000]] } }), t);
  }, 'INVALID_PARAMS');
  design.validateAgainstTemplate(norm({ mode: 'sequence', sequence: 'A', params: { product_size_ranges: [[1001, 2000], [100, 1000]] } }), t);
});

test('junctionList: transcripts with junction_spanning only; SINGLE_EXON_TRANSCRIPT; JUNCTIONS_TRUNCATED keeps 200', function () {
  const w = [];
  should(design.junctionList(norm({ mode: 'gene', gene_id: 'G' }), { junctions: [10] }, w)).be.null();
  should(design.junctionList(norm({ mode: 'transcript', gene_id: 'G', junction_spanning: false }), { junctions: [10] }, w)).be.null();
  w.should.eql([]);
  design.junctionList(norm({ mode: 'transcript', gene_id: 'G' }), { junctions: [10, 20] }, w).should.eql([10, 20]);
  should(design.junctionList(norm({ mode: 'transcript', gene_id: 'G' }), { junctions: [] }, w)).be.null();
  w.map(function (x) { return x.code; }).should.eql(['SINGLE_EXON_TRANSCRIPT']);
  const many = [];
  for (let i = 1; i <= 250; i++) many.push(i * 10);
  const w2 = [];
  design.junctionList(norm({ mode: 'transcript', gene_id: 'G' }), { junctions: many }, w2).should.eql(many.slice(0, 200));
  w2.map(function (x) { return x.code; }).should.eql(['JUNCTIONS_TRUNCATED']);
  const w3 = [];
  const inside = design.junctionList(norm({ mode: 'transcript', gene_id: 'G', included: [2001, 300] }), { junctions: many }, w3);
  inside.should.eql(many.filter(function (j) { return j >= 2001 && j < 2300; }));
  w3.map(function (x) { return x.code; }).should.eql(['JUNCTIONS_TRUNCATED']);
});

test('FASTA cleaning for sequence mode (template.cleanSequence)', function () {
  template.cleanSequence('>seq1\nACGT ACGT\n10 ACGTACGT\r\nacgtacgt\n').seq.should.equal('ACGTACGTACGTACGTACGTACGT');
  const r = template.cleanSequence('>a\nACGTACGTRYACGTACGTACGT\n');
  r.seq.should.equal('ACGTACGTNNACGTACGTACGT');
  r.iupac.should.equal(2);
  throwsCode(function () { template.cleanSequence('>a\nACGTACGTACGTACGTACGTXACGT'); }, 'INVALID_SEQUENCE', 400);
});

// ---- Boulder record ------------------------------------------------------------------------------

test('buildRecord reproduces the captured Primer3 input records byte for byte (transcript qPCR and 50 kb PCR)', function () {
  const tj = stubs.primer3Tags('transcript_junction.input.txt');
  const junctions = tj.SEQUENCE_OVERLAP_JUNCTION_LIST.split(' ').map(Number);
  boulder.serialize(design.buildRecord({ id: 'SORBI_3001G000200.1', seq: tj.SEQUENCE_TEMPLATE }, { params: design.PRESETS.qpcr, junctions: junctions }))
    .should.equal(stubs.primer3Text('transcript_junction.input.txt'));
  const g = stubs.primer3Tags('genomic_50kb.input.txt');
  boulder.serialize(design.buildRecord({ id: 'sorghum_bicolor_4_7400001-7450000', seq: g.SEQUENCE_TEMPLATE },
    { params: Object.assign({}, design.PRESETS.pcr, { num_return: 20 }) }))
    .should.equal(stubs.primer3Text('genomic_50kb.input.txt'));
});

test('buildRecord: intervals, optional params, masking modes and junction tags', function () {
  const params = Object.assign({}, design.PRESETS.pcr, { opt_gc: 50, gc_clamp: 1, max_end_stability: 9, max_ns: 2, salt_monovalent: 50,
    salt_divalent: 1.5, dntp_conc: 0.6, dna_conc: 50 });
  const tags = design.buildRecord({ id: 'a b\nc', seq: 'ACGT' },
    { params: params, target: [499, 50], included: [100, 3000], excluded: [[1000, 40], [2000, 10]], avoid_repeats: false, repeat_mask_mode: null });
  Object.keys(tags).slice(0, 5).should.eql(['SEQUENCE_ID', 'SEQUENCE_TEMPLATE', 'SEQUENCE_TARGET', 'SEQUENCE_INCLUDED_REGION', 'SEQUENCE_EXCLUDED_REGION']);
  tags.should.have.properties({
    SEQUENCE_ID: 'a_b_c', SEQUENCE_TARGET: '499,50', SEQUENCE_INCLUDED_REGION: '100,3000', SEQUENCE_EXCLUDED_REGION: '1000,40 2000,10',
    PRIMER_OPT_GC_PERCENT: 50, PRIMER_GC_CLAMP: 1, PRIMER_MAX_END_STABILITY: 9, PRIMER_MAX_NS_ACCEPTED: 2, PRIMER_SALT_MONOVALENT: 50,
    PRIMER_SALT_DIVALENT: 1.5, PRIMER_DNTP_CONC: 0.6, PRIMER_DNA_CONC: 50, PRIMER_PRODUCT_MIN_TM: 0, PRIMER_PRODUCT_MAX_TM: 150,
    PRIMER_LIBERAL_BASE: 1, PRIMER_FIRST_BASE_INDEX: 1, PRIMER_EXPLAIN_FLAG: 1, PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT: 1
  });
  tags.should.not.have.property('SEQUENCE_OVERLAP_JUNCTION_LIST');
  tags.should.not.have.property('PRIMER_LOWERCASE_MASKING');
  boulder.serialize(tags).should.match(/\nPRIMER_DNTP_CONC=0\.6\n/);

  const n = design.buildRecord({ id: 'x', seq: 'ACGTNNNN' }, { params: params, avoid_repeats: true, repeat_mask_mode: 'n_mask' });
  n.PRIMER_MAX_NS_ACCEPTED.should.equal(0);
  n.should.not.have.property('PRIMER_LOWERCASE_MASKING');
  const withoutNs = design.buildRecord({ id: 'x', seq: 'ACGT' }, { params: design.PRESETS.pcr, avoid_repeats: true, repeat_mask_mode: 'n_mask' });
  withoutNs.PRIMER_MAX_NS_ACCEPTED.should.equal(0);
  const l = design.buildRecord({ id: 'x', seq: 'ACGTacgt' }, { params: params, avoid_repeats: true, repeat_mask_mode: 'three_prime' });
  l.PRIMER_LOWERCASE_MASKING.should.equal(1);
  l.PRIMER_MAX_NS_ACCEPTED.should.equal(2);

  const j = design.buildRecord({ id: 'x', seq: 'ACGT' },
    { params: Object.assign({}, design.PRESETS.qpcr, { min_3_prime_overlap_of_junction: 5, min_5_prime_overlap_of_junction: 9 }), junctions: [10, 20] });
  j.should.have.properties({
    SEQUENCE_OVERLAP_JUNCTION_LIST: '10 20', PRIMER_MIN_3_PRIME_OVERLAP_OF_JUNCTION: 5, PRIMER_MIN_5_PRIME_OVERLAP_OF_JUNCTION: 9,
    PRIMER_INTERNAL_MIN_3_PRIME_OVERLAP_OF_JUNCTION: 5, PRIMER_INTERNAL_MIN_5_PRIME_OVERLAP_OF_JUNCTION: 9
  });
  (function () { boulder.serialize(design.buildRecord({ id: 'x', seq: 'ACGT\nPRIMER_TASK=pick_sequencing_primers' }, { params: design.PRESETS.pcr })); })
    .should.throw({ code: 'BOULDER_INVALID_TAG' });
});

// ---- design() orchestration (no primer3 binary, no mongo) -----------------------------------------

function fixtureResult(name) {
  const tags = stubs.primer3Tags(name);
  return { tags: tags, error: tags.PRIMER_ERROR || null, warning: tags.PRIMER_WARNING || null, exitCode: 0 };
}

function primer3Stub(result, calls) {
  return {
    run: async function (tags, opts) {
      if (calls) calls.push({ tags: tags, opts: opts });
      return typeof result === 'function' ? result(tags, opts) : result;
    },
    version: async function () { return '2.6.1'; }
  };
}

function semaphoreStub(order) {
  return {
    acquire: async function () {
      order.push('acquire');
      return function () { order.push('release'); };
    }
  };
}

// Primer3-like output for pairs [{left: [start, len], right: [end, len]}] on a template sequence.
function syntheticResult(tseq, pairs) {
  const tags = { PRIMER_FIRST_BASE_INDEX: '1', PRIMER_PAIR_NUM_RETURNED: String(pairs.length),
    PRIMER_LEFT_EXPLAIN: 'considered 10, ok 2', PRIMER_RIGHT_EXPLAIN: 'considered 10, ok 2', PRIMER_PAIR_EXPLAIN: 'considered 1, ok 1' };
  pairs.forEach(function (p, i) {
    tags['PRIMER_PAIR_' + i + '_PENALTY'] = '0.5';
    tags['PRIMER_LEFT_' + i] = p.left[0] + ',' + p.left[1];
    tags['PRIMER_RIGHT_' + i] = p.right[0] + ',' + p.right[1];
    tags['PRIMER_LEFT_' + i + '_SEQUENCE'] = tseq.slice(p.left[0] - 1, p.left[0] - 1 + p.left[1]);
    tags['PRIMER_RIGHT_' + i + '_SEQUENCE'] = revcomp(tseq.slice(p.right[0] - p.right[1], p.right[0]));
    tags['PRIMER_LEFT_' + i + '_TM'] = '60.1';
    tags['PRIMER_RIGHT_' + i + '_TM'] = '59.9';
    tags['PRIMER_PAIR_' + i + '_PRODUCT_SIZE'] = String(p.right[0] - p.left[0] + 1);
    tags['PRIMER_PAIR_' + i + '_PRODUCT_TM'] = '80.5';
  });
  return { tags: tags, error: null, warning: null, exitCode: 0 };
}

test('design (transcript): record equals the captured input; pairs carry junctions, genome-verified blocks and genomic_size', async function () {
  const calls = [];
  const order = [];
  const deps = stubs.sorghumDeps({ primer3: primer3Stub(fixtureResult('transcript_junction.output.txt'), calls), semaphore: semaphoreStub(order) });
  const r = await design.design({ mode: 'transcript', gene_id: G200 }, deps);
  calls.length.should.equal(1);
  boulder.serialize(calls[0].tags).should.equal(stubs.primer3Text('transcript_junction.input.txt'));
  calls[0].opts.timeoutMs.should.equal(30000);
  should.exist(calls[0].opts.signal);
  order.should.eql(['acquire', 'release']);

  Object.keys(r).should.eql(['template', 'pairs', 'explain', 'settings', 'engine', 'warnings']);
  Object.keys(r.template).should.eql(['mode', 'system_name', 'gene_id', 'transcript_id', 'region', 'start', 'end', 'strand', 'length', 'seq',
    'masked', 'mask_source', 'mask', 'masked_fraction', 'features']);
  r.template.should.have.properties({ mode: 'transcript', system_name: 'sorghum_bicolor', gene_id: G200, transcript_id: 'SORBI_3001G000200.1',
    region: '1', start: 11180, end: 14899, strand: -1, length: 1982, masked: false, mask_source: null, masked_fraction: 0 });
  r.template.mask.should.eql([]);
  r.template.seq.should.equal(stubs.cdna200());
  r.template.features.junctions.slice(0, 3).should.eql([397, 493, 597]);

  r.pairs.length.should.equal(5);
  r.pairs.forEach(function (p, i) {
    p.rank.should.equal(i);
    (p.left.junction !== null || p.right.junction !== null).should.be.true();
    p.product_size.should.be.within(70, 150);
    p.product_tm.should.be.a.Number();
    let lo = Infinity;
    let hi = -Infinity;
    ['left', 'right'].forEach(function (side) {
      const o = p[side];
      o.seq.should.equal(o.seq.toUpperCase());
      o.seq.length.should.equal(o.len);
      r.template.seq.slice(o.start - 1, o.end).should.equal(side === 'left' ? o.seq : revcomp(o.seq));
      o.genomic.region.should.equal('1');
      o.genomic.strand.should.equal(side === 'left' ? -1 : 1);
      stubs.blocksSeq(G200, o.genomic.blocks, o.genomic.strand).should.equal(o.seq);
      o.genomic.start.should.equal(o.genomic.blocks[0].start);
      o.genomic.end.should.equal(o.genomic.blocks[o.genomic.blocks.length - 1].end);
      if (o.junction) {
        o.junction.overlap_5p.should.be.aboveOrEqual(7);
        o.junction.overlap_3p.should.be.aboveOrEqual(4);
        o.genomic.blocks.length.should.equal(2);
      }
      o.genomic.blocks.forEach(function (b) { lo = Math.min(lo, b.start); hi = Math.max(hi, b.end); });
    });
    p.product.genomic.should.eql({ region: '1', start: lo, end: hi, strand: -1 });
    p.product.genomic_size.should.equal(hi - lo + 1);
  });
  r.settings.should.eql({ preset: 'qpcr', junction_spanning: true, avoid_repeats: false, repeat_mask_mode: null, params: design.PRESETS.qpcr });
  r.engine.should.eql({ primer3: '2.6.1' });
  r.warnings.should.eql([]);
  r.explain.left.raw.should.be.a.String();
  r.explain.pair.should.have.property('ok');
});

test('design (gene, - strand): genomic mapping of both primers and the product; no genomic_size; product_tm only when trusted', async function () {
  const calls = [];
  const deps = stubs.sorghumDeps({
    primer3: primer3Stub(function (tags) { return syntheticResult(tags.SEQUENCE_TEMPLATE, [{ left: [500, 20], right: [719, 20] }]); }, calls),
    semaphore: semaphoreStub([])
  });
  const r = await design.design({ mode: 'gene', gene_id: G200, flank_up: 200, flank_down: 100 }, deps);
  calls[0].tags.SEQUENCE_ID.should.equal(G200);
  calls[0].tags.should.not.have.property('SEQUENCE_OVERLAP_JUNCTION_LIST');
  calls[0].tags.SEQUENCE_TEMPLATE.should.equal(revcomp(stubs.windowSeq(G200)));
  const p = r.pairs[0];
  p.left.genomic.should.eql({ region: '1', start: 14581, end: 14600, strand: -1, blocks: [{ start: 14581, end: 14600 }] });
  p.right.genomic.should.eql({ region: '1', start: 14381, end: 14400, strand: 1, blocks: [{ start: 14381, end: 14400 }] });
  stubs.blocksSeq(G200, p.left.genomic.blocks, -1).should.equal(p.left.seq);
  stubs.blocksSeq(G200, p.right.genomic.blocks, 1).should.equal(p.right.seq);
  p.product.should.eql({ start: 500, end: 719, genomic: { region: '1', start: 14381, end: 14600, strand: -1 } });
  p.product.should.not.have.property('genomic_size');
  should(p.left.junction).be.null();
  should(p.product_tm).be.null();
  r.template.features.cds.start.should.equal(499);
  r.settings.should.have.properties({ preset: 'pcr', junction_spanning: false });
});

test('design: PRIMER_ERROR -> 400 PRIMER3_INPUT_ERROR with details.primer3_error; the semaphore is released', async function () {
  const order = [];
  const deps = stubs.sorghumDeps({ primer3: primer3Stub(fixtureResult('error_junction_overlap.output.txt')), semaphore: semaphoreStub(order) });
  const err = await rejectsCode(design.design({ mode: 'transcript', gene_id: G200 }, deps), 'PRIMER3_INPUT_ERROR', 400);
  err.details.primer3_error.should.equal('PRIMER_MIN_5_PRIME_OVERLAP_OF_JUNCTION > PRIMER_MAX_SIZE / 2');
  order.should.eql(['acquire', 'release']);
});

test('design: PRIMER_WARNING -> PRIMER3_WARNING; zero pairs -> 200 with NO_PAIRS and explain', async function () {
  const warned = fixtureResult('transcript_junction.output.txt');
  warned.warning = 'Unrecognized base in input sequence';
  const r1 = await design.design({ mode: 'transcript', gene_id: G200 }, stubs.sorghumDeps({ primer3: primer3Stub(warned), semaphore: semaphoreStub([]) }));
  r1.warnings.should.eql([{ code: 'PRIMER3_WARNING', message: 'Unrecognized base in input sequence' }]);
  const r2 = await design.design({ mode: 'transcript', gene_id: G200 },
    stubs.sorghumDeps({ primer3: primer3Stub(fixtureResult('no_pairs.output.txt')), semaphore: semaphoreStub([]) }));
  r2.pairs.should.eql([]);
  r2.warnings.map(function (w) { return w.code; }).should.eql(['NO_PAIRS']);
  r2.explain.left.raw.should.be.a.String();
});

test('design: template_only returns the template and skips Primer3', async function () {
  const calls = [];
  const r = await design.design({ mode: 'gene', gene_id: G200, template_only: true },
    stubs.sorghumDeps({ primer3: primer3Stub(fixtureResult('no_pairs.output.txt'), calls), semaphore: semaphoreStub([]) }));
  calls.length.should.equal(0);
  r.pairs.should.eql([]);
  should(r.explain).be.null();
  r.engine.should.eql({ primer3: '2.6.1' });
  r.template.length.should.equal(3720);
  r.warnings.should.eql([]);
});

test('design: the semaphore is acquired before the template is built; BUSY rejects before any lookup', async function () {
  const order = [];
  const deps = stubs.sorghumDeps({
    primer3: primer3Stub(fixtureResult('no_pairs.output.txt')),
    semaphore: semaphoreStub(order),
    template: { buildTemplate: async function (req, d) { order.push('build'); return template.buildTemplate(req, d); } }
  });
  await design.design({ mode: 'gene', gene_id: G200 }, deps);
  order.should.eql(['acquire', 'build', 'release']);

  let looked = 0;
  const busy = stubs.sorghumDeps({
    findGene: async function () { looked++; return null; },
    semaphore: { acquire: async function () { throw new PrimerHttpError(503, 'BUSY', 'busy', { retry_after_s: 5 }); } }
  });
  await rejectsCode(design.design({ mode: 'gene', gene_id: G200 }, busy), 'BUSY', 503);
  looked.should.equal(0);
});

test('design: validation errors and FEATURE_DISABLED are raised before the semaphore', async function () {
  const order = [];
  await rejectsCode(design.design({ mode: 'transcript', gene_id: G200, params: { max_size: 24, min_5_prime_overlap_of_junction: 13 } },
    stubs.sorghumDeps({ semaphore: semaphoreStub(order) })), 'INVALID_PARAMS', 400);
  await rejectsCode(design.design({ mode: 'gene', gene_id: G200 }, stubs.sorghumDeps({ cfg: stubs.cfg({ enabled: false }), semaphore: semaphoreStub(order) })),
    'FEATURE_DISABLED', 503);
  order.should.eql([]);
});

test('design: one deadline covers the template build; 504 DEADLINE_EXCEEDED releases the semaphore', async function () {
  const order = [];
  const deps = stubs.sorghumDeps({
    cfg: stubs.cfg({ design: { deadline_ms: 60 } }),
    semaphore: semaphoreStub(order),
    template: { buildTemplate: function () { return new Promise(function () {}); } }
  });
  const started = Date.now();
  const err = await rejectsCode(design.design({ mode: 'gene', gene_id: G200 }, deps), 'DEADLINE_EXCEEDED', 504);
  (Date.now() - started).should.be.below(2000);
  err.details.deadline_ms.should.equal(60);
  order.should.eql(['acquire', 'release']);
});

test('design: Primer3 gets min(primer3_timeout_ms, time left) and the shared abort signal', async function () {
  const calls = [];
  await design.design({ mode: 'gene', gene_id: G200 },
    stubs.sorghumDeps({ cfg: stubs.cfg({ design: { deadline_ms: 5000 } }), primer3: primer3Stub(fixtureResult('no_pairs.output.txt'), calls), semaphore: semaphoreStub([]) }));
  calls[0].opts.timeoutMs.should.be.within(4000, 5000);
  const order = [];
  const hang = {
    run: function (tags, opts) {
      return new Promise(function (resolve, reject) {
        opts.signal.addEventListener('abort', function () { reject(opts.signal.reason); });
      });
    },
    version: async function () { return '2.6.1'; }
  };
  await rejectsCode(design.design({ mode: 'gene', gene_id: G200 },
    stubs.sorghumDeps({ cfg: stubs.cfg({ design: { deadline_ms: 150 } }), primer3: hang, semaphore: semaphoreStub(order) })), 'DEADLINE_EXCEEDED', 504);
  order.should.eql(['acquire', 'release']);
});

test('design (sequence): IUPAC codes are converted (IUPAC_CONVERTED), never sent to Primer3; genomic fields are null', async function () {
  const calls = [];
  const raw = '>amp\n' + stubs.cdna200().slice(0, 100) + 'RY' + stubs.cdna200().slice(100, 500) + 'K';
  const r = await design.design({ mode: 'sequence', sequence: raw },
    stubs.sorghumDeps({ primer3: primer3Stub(function (tags) { return syntheticResult(tags.SEQUENCE_TEMPLATE, [{ left: [10, 20], right: [300, 20] }]); }, calls),
      semaphore: semaphoreStub([]) }));
  calls[0].tags.SEQUENCE_TEMPLATE.should.match(/^[ACGTN]+$/);
  calls[0].tags.SEQUENCE_TEMPLATE.slice(100, 102).should.equal('NN');
  calls[0].tags.SEQUENCE_ID.should.equal('amp');
  r.warnings.map(function (w) { return w.code; }).should.eql(['IUPAC_CONVERTED']);
  [r.template.system_name, r.template.region, r.template.start, r.template.strand].should.eql([null, null, null, null]);
  should(r.pairs[0].left.genomic).be.null();
  should(r.pairs[0].right.genomic).be.null();
  should(r.pairs[0].product.genomic).be.null();
  r.pairs[0].product.should.not.have.property('genomic_size');
});

test('design: avoid_repeats passes the masked sequence to Primer3 (N + PRIMER_MAX_NS_ACCEPTED=0) and returns the mask', async function () {
  const calls = [];
  const seen = [];
  const masker = {
    repeatMask: async function (t, opts, d) {
      seen.push(opts);
      return { mask: [[1, 10]], mask_source: 'softmask', masked: true, masked_fraction: 0.0027, seq: 'N'.repeat(10) + t.seq.slice(10).toUpperCase(),
        mode: opts.mode, tags: {}, warnings: [{ code: 'MOSTLY_REPEAT', message: 'x' }] };
    }
  };
  const r = await design.design({ mode: 'gene', gene_id: G200, avoid_repeats: true, params: { max_ns: 2 } },
    stubs.sorghumDeps({ repeatMask: masker, primer3: primer3Stub(fixtureResult('no_pairs.output.txt'), calls), semaphore: semaphoreStub([]) }));
  seen[0].mode.should.equal('n_mask');
  should.exist(seen[0].signal);
  seen[0].deadline.should.be.a.Number();
  // no flanks: the template is the 3720 bp gene span; Primer3 gets the masked copy, the response the unmasked one
  r.template.length.should.equal(3720);
  calls[0].tags.SEQUENCE_TEMPLATE.should.equal('NNNNNNNNNN' + r.template.seq.slice(10));
  calls[0].tags.PRIMER_MAX_NS_ACCEPTED.should.equal(0);
  r.template.should.have.properties({ masked: true, mask_source: 'softmask', masked_fraction: 0.0027 });
  r.template.mask.should.eql([[1, 10]]);
  // the captured window is 1:11080-15099; on the - strand the gene (11180-14899) starts 200 bases in
  r.template.seq.should.equal(revcomp(stubs.windowSeq(G200)).slice(200, 3920));
  r.settings.should.have.properties({ avoid_repeats: true, repeat_mask_mode: 'n_mask' });
  r.settings.params.max_ns.should.equal(0);
  r.warnings.map(function (w) { return w.code; }).should.eql(['MOSTLY_REPEAT', 'NO_PAIRS']);
});

test('design: single-exon transcript drops the junction constraint with SINGLE_EXON_TRANSCRIPT; intervals beyond the cDNA are rejected', async function () {
  const doc = stubs.geneDoc(G200);
  doc.gene_structure.transcripts = [{ id: 'S.1', length: 397, exons: ['EER93047-1'] }];
  doc.gene_structure.canonical_transcript = 'S.1';
  const calls = [];
  const r = await design.design({ mode: 'transcript', gene_id: G200 },
    stubs.sorghumDeps({ findGene: async function () { return doc; }, primer3: primer3Stub(fixtureResult('no_pairs.output.txt'), calls), semaphore: semaphoreStub([]) }));
  calls[0].tags.should.not.have.property('SEQUENCE_OVERLAP_JUNCTION_LIST');
  r.warnings.map(function (w) { return w.code; }).should.eql(['SINGLE_EXON_TRANSCRIPT', 'NO_PAIRS']);
  r.settings.junction_spanning.should.be.true();

  const order = [];
  const calls2 = [];
  await rejectsCode(design.design({ mode: 'transcript', gene_id: G200, target: [1900, 100] },
    stubs.sorghumDeps({ primer3: primer3Stub(fixtureResult('no_pairs.output.txt'), calls2), semaphore: semaphoreStub(order) })), 'INTERVAL_OUT_OF_BOUNDS', 400);
  calls2.length.should.equal(0);
  order.should.eql(['acquire', 'release']);
});

test('design: lookup errors from the template surface unchanged (UNKNOWN_GENE 404, SYSTEM_NAME_MISMATCH 400)', async function () {
  await rejectsCode(design.design({ mode: 'gene', gene_id: 'NOPE' }, stubs.sorghumDeps({ semaphore: semaphoreStub([]) })), 'UNKNOWN_GENE', 404);
  await rejectsCode(design.design({ mode: 'gene', gene_id: G200, system_name: 'sorghum_rio' }, stubs.sorghumDeps({ semaphore: semaphoreStub([]) })),
    'SYSTEM_NAME_MISMATCH', 400);
});

// ---- template_only, multi-record FASTA, n_mask product Tm -------------------------------------------

test('template_only: product ranges are not checked against the template length; interval and param checks still apply', async function () {
  const t = { length: 80 };
  design.validateAgainstTemplate(norm({ mode: 'sequence', sequence: 'A', template_only: true }), t);
  throwsCode(function () { design.validateAgainstTemplate(norm({ mode: 'sequence', sequence: 'A' }), t); }, 'INVALID_PARAMS', 400);
  throwsCode(function () { design.validateAgainstTemplate(norm({ mode: 'sequence', sequence: 'A', template_only: true, target: [70, 20] }), t); },
    'INTERVAL_OUT_OF_BOUNDS', 400);
  throwsCode(function () { norm({ mode: 'sequence', sequence: 'A', template_only: true, params: { product_size_ranges: [[20, 100]] } }); },
    'INVALID_PARAMS', 400);

  const seq80 = stubs.cdna200().slice(0, 80);
  const calls = [];
  const deps = function () {
    return stubs.sorghumDeps({ primer3: primer3Stub(fixtureResult('no_pairs.output.txt'), calls), semaphore: semaphoreStub([]) });
  };
  const r = await design.design({ mode: 'sequence', sequence: seq80, template_only: true }, deps());
  calls.length.should.equal(0);
  r.template.length.should.equal(80);
  r.template.seq.should.equal(seq80);
  r.pairs.should.eql([]);
  should(r.explain).be.null();
  r.warnings.should.eql([]);
  r.settings.params.product_size_ranges.should.eql([[100, 1000]]);
  await rejectsCode(design.design({ mode: 'sequence', sequence: seq80 }, deps()), 'INVALID_PARAMS', 400);
  await rejectsCode(design.design({ mode: 'sequence', sequence: seq80, template_only: true, excluded: [[75, 10]] }, deps()), 'INTERVAL_OUT_OF_BOUNDS', 400);
  calls.length.should.equal(0);
});

test('design (sequence): a multi-record FASTA still designs, with warning MULTIPLE_RECORDS', async function () {
  const cdna = stubs.cdna200();
  const calls = [];
  const r = await design.design({ mode: 'sequence', sequence: '>recA\n' + cdna.slice(0, 300) + '\n>recB\n' + cdna.slice(600, 900) + '\n' },
    stubs.sorghumDeps({ primer3: primer3Stub(function (tags) { return syntheticResult(tags.SEQUENCE_TEMPLATE, [{ left: [250, 20], right: [400, 20] }]); }, calls),
      semaphore: semaphoreStub([]) }));
  calls.length.should.equal(1);
  calls[0].tags.SEQUENCE_TEMPLATE.should.equal(cdna.slice(0, 300) + cdna.slice(600, 900));
  r.template.length.should.equal(600);
  r.warnings.should.eql([{ code: 'MULTIPLE_RECORDS', message: '2 FASTA records were joined into one template; primers may span the joins' }]);
  r.pairs.length.should.equal(1);
});

test('longSeqTm reproduces primer3_core\'s PRIMER_PAIR_i_PRODUCT_TM in the captured records (default salts); divalent <= dNTP adds nothing', function () {
  let n = 0;
  ['genomic_50kb', 'transcript_junction'].forEach(function (name) {
    const input = stubs.primer3Tags(name + '.input.txt');
    boulder.extractPairs(stubs.primer3Tags(name + '.output.txt')).forEach(function (p) {
      const tm = design.longSeqTm(input.SEQUENCE_TEMPLATE.slice(p.product.start - 1, p.product.end).toUpperCase(), {});
      Math.abs(tm - p.product_tm).should.be.below(0.0001);
      n++;
    });
  });
  n.should.equal(25);
  const low = design.longSeqTm('GGCCAATT', { salt_monovalent: 50, salt_divalent: 0.5, dntp_conc: 0.6 });
  low.should.equal(design.longSeqTm('GGCCAATT', { salt_monovalent: 50, salt_divalent: 0, dntp_conc: 0 }));
  design.longSeqTm('GGCCAATT', {}).should.be.approximately(low + 16.6 * Math.log10((50 + 120 * Math.sqrt(0.9)) / 50), 1e-9);
  design.longSeqTm('ggccNNNN', {}).should.equal(design.longSeqTm('AAAAAAAA', {})); // only uppercase G/C count, like Primer3
  should(design.longSeqTm('ACGT', { salt_monovalent: 0, salt_divalent: 0 })).be.null();
});

test('design (n_mask): product_tm of a product covering masked bases is recomputed on the unmasked template; others untouched', async function () {
  const cdna = stubs.cdna200();
  const seq = cdna.slice(0, 300) + 'gcggccgcgc'.repeat(10) + cdna.slice(400, 700); // user mask 301-400, all G/C
  const params = { salt_monovalent: 40, salt_divalent: 2, dntp_conc: 0.8 };
  const PAIRS = [{ left: [100, 20], right: [600, 20] }, { left: [420, 20], right: [650, 20] }];
  // what primer3_core does: long_seq_tm on the template it was given (N counts as A/T), with the record's salts
  function simulatedPrimer3(tags) {
    const res = syntheticResult(tags.SEQUENCE_TEMPLATE, PAIRS);
    const salts = { salt_monovalent: tags.PRIMER_SALT_MONOVALENT, salt_divalent: tags.PRIMER_SALT_DIVALENT, dntp_conc: tags.PRIMER_DNTP_CONC };
    PAIRS.forEach(function (p, i) {
      res.tags['PRIMER_PAIR_' + i + '_PRODUCT_TM'] = design.longSeqTm(tags.SEQUENCE_TEMPLATE.slice(p.left[0] - 1, p.right[0]).toUpperCase(), salts).toFixed(4);
      res.tags['PRIMER_PAIR_' + i + '_PRODUCT_TM_OLIGO_TM_DIFF'] = '20.0000';
    });
    return res;
  }
  const calls = [];
  const r = await design.design({ mode: 'sequence', sequence: seq, avoid_repeats: true, params: params },
    stubs.sorghumDeps({ primer3: primer3Stub(simulatedPrimer3, calls), semaphore: semaphoreStub([]) }));
  r.template.mask_source.should.equal('user_lowercase');
  r.template.mask.should.eql([[301, 100]]);
  calls[0].tags.SEQUENCE_TEMPLATE.slice(300, 400).should.equal('N'.repeat(100));
  calls[0].tags.PRIMER_SALT_DIVALENT.should.equal(2);
  const upper = seq.toUpperCase();
  const unmasked = design.longSeqTm(upper.slice(99, 600), params);
  const asReported = Number(calls[0].tags.SEQUENCE_TEMPLATE && design.longSeqTm(calls[0].tags.SEQUENCE_TEMPLATE.slice(99, 600), params).toFixed(4));
  (unmasked - asReported).should.be.approximately(41 * 100 / 501, 0.001); // every masked G/C was counted as A/T
  r.pairs[0].product_tm.should.be.approximately(unmasked, 0.00005);
  r.pairs[0].product_tm.should.equal(Math.round(unmasked * 10000) / 10000);
  r.pairs[1].product_tm.should.equal(Number(design.longSeqTm(upper.slice(419, 650), params).toFixed(4))); // no masked base: Primer3's value

  const fixed = function (tags) {
    const res = syntheticResult(tags.SEQUENCE_TEMPLATE, PAIRS.slice(0, 1));
    res.tags.PRIMER_PAIR_0_PRODUCT_TM = '77.7777';
    res.tags.PRIMER_PAIR_0_PRODUCT_TM_OLIGO_TM_DIFF = '17.7777';
    return res;
  };
  const three = await design.design({ mode: 'sequence', sequence: seq, avoid_repeats: true, repeat_mask_mode: 'three_prime' },
    stubs.sorghumDeps({ primer3: primer3Stub(fixed), semaphore: semaphoreStub([]) }));
  three.pairs[0].product_tm.should.equal(77.7777); // Primer3 upcases lowercase itself: nothing to correct
  const noMask = await design.design({ mode: 'sequence', sequence: seq },
    stubs.sorghumDeps({ primer3: primer3Stub(fixed), semaphore: semaphoreStub([]) }));
  noMask.pairs[0].product_tm.should.equal(77.7777);
  const untrusted = await design.design({ mode: 'sequence', sequence: seq, avoid_repeats: true },
    stubs.sorghumDeps({ primer3: primer3Stub(function (tags) { return syntheticResult(tags.SEQUENCE_TEMPLATE, PAIRS.slice(0, 1)); }), semaphore: semaphoreStub([]) }));
  should(untrusted.pairs[0].product_tm).be.null();
});

const REAL_PRIMER3 = fs.existsSync(CFG.primer3_core);

test('design (n_mask, real primer3_core): product_tm equals primer3_core\'s own PRODUCT_TM for the same pair on the unmasked template (<= 0.01 C)',
  { skip: !REAL_PRIMER3 && 'primer3_core not installed' }, async function () {
    const primer3 = require('../../../api/helpers/primers/primer3');
    const win = stubs.windowSeq(G200); // sorghum_bicolor 1:11080-15099 (+), 4020 bp
    const seq = win.slice(0, 1499) + win.slice(1499, 1799).toLowerCase() + win.slice(1799); // user mask 1500-1799
    const params = { product_size_ranges: [[400, 900]], salt_monovalent: 40, salt_divalent: 2.5, dntp_conc: 0.8 };
    const r = await design.design({ mode: 'sequence', sequence: seq, avoid_repeats: true, target: [1600, 50], params: params },
      { cfg: CFG, log: stubs.silentLog, semaphore: semaphoreStub([]) });
    r.template.mask.should.eql([[1500, 300]]);
    r.pairs.length.should.be.above(0);
    const upper = seq.toUpperCase();
    const nMasked = upper.slice(0, 1499) + 'N'.repeat(300) + upper.slice(1799);
    for (const p of r.pairs) {
      coords.overlapsRuns(p.product.start, p.product.end, r.template.mask).should.be.true();
      const forced = await primer3.run(stubs.forcedPairRecord(upper, p, params), { timeoutMs: 30000 });
      should(forced.error).be.null();
      forced.tags.PRIMER_LEFT_0.should.equal(p.left.start + ',' + p.left.len);
      forced.tags.PRIMER_RIGHT_0.should.equal(p.right.end + ',' + p.right.len);
      const own = Number(forced.tags.PRIMER_PAIR_0_PRODUCT_TM);
      Math.abs(p.product_tm - own).should.be.belowOrEqual(0.01);
      (own - design.longSeqTm(nMasked.slice(p.product.start - 1, p.product.end), params)).should.be.above(1); // the N-masked value was wrong
    }
  });
