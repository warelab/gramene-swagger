'use strict';

require('../../../api/helpers/primers/node_compat');

const test = require('node:test');
const should = require('should');

const request = require('../../../api/helpers/primers/genotyping/request');
const presets = require('../../../api/helpers/primers/genotyping/presets');
const design = require('../../../api/helpers/primers/design');
const cases = require('../fixtures/primer3/genotyping/cases');

// Spec §2.7 handler rules 1-6, §4.6 presets and ladder, and the genotyping_request rows of §7.3. Pure, offline.

const CFG = cases.cfg();

function body(extra) {
  return Object.assign({ system_name: 'sorghum_bicolor', variant: { id: 'rs871475760', alt: 'A' } }, extra || {});
}

function norm(b, cfg) {
  return request.normalize(b, cfg || CFG);
}

function throwsCode(fn, code, details) {
  let err = null;
  try {
    fn();
  } catch (e) {
    err = e;
  }
  should.exist(err, 'expected ' + code);
  err.should.have.property('code', code);
  err.should.have.property('status', 400);
  if (details) err.details.should.containEql(details);
  return err;
}

// ---- presets --------------------------------------------------------------------------------------------

test('presets: frozen tables of §4.6 and the §3.1 floors', function () {
  Object.isFrozen(presets.PRESETS.kasp).should.be.true();
  Object.isFrozen(presets.LADDER.as_pcr[1]).should.be.true();
  Object.isFrozen(presets.FLOORS).should.be.true();
  presets.FLOORS.should.eql({ as_min_tm: 52, as_min_gc: 15 });
  presets.levelParams('kasp', 0).should.eql({ opt_size: 22, min_size: 18, max_size: 30, opt_tm: 60, min_tm: 57, max_tm: 63,
    min_gc: 30, max_gc: 70, max_tm_diff: 3, max_poly_x: 5, product_size_ranges: [[61, 120]] });
  presets.levelParams('kasp', 1).should.eql({ opt_size: 22, min_size: 18, max_size: 32, opt_tm: 60, min_tm: 55, max_tm: 65,
    min_gc: 20, max_gc: 80, max_tm_diff: 3, max_poly_x: 5, product_size_ranges: [[65, 150]] });
  presets.levelParams('kasp', 2).should.containEql({ max_size: 32, min_tm: 52, max_tm: 65, max_tm_diff: 6, product_size_ranges: [[65, 150]] });
  presets.levelParams('as_pcr', 0).should.eql({ opt_size: 24, min_size: 18, max_size: 30, opt_tm: 60, min_tm: 57, max_tm: 63,
    min_gc: 30, max_gc: 70, max_tm_diff: 3, max_poly_x: 4, product_size_ranges: [[150, 300]] });
  presets.levelParams('as_pcr', 2).should.containEql({ max_size: 32, min_tm: 52, max_tm: 65, min_gc: 20, max_gc: 80, max_tm_diff: 6,
    product_size_ranges: [[150, 300]] });
  // A fresh object every time: mutating it never reaches the frozen tables.
  const p = presets.levelParams('kasp', 0);
  p.product_size_ranges[0][0] = 1;
  presets.levelParams('kasp', 0).product_size_ranges.should.eql([[61, 120]]);
  (function () { presets.levelParams('pcr', 0); }).should.throw(TypeError);
  (function () { presets.levelParams('kasp', 3); }).should.throw(RangeError);
});

test('presets: paramChanges lists what a level changes, in key order', function () {
  presets.paramChanges(presets.levelParams('kasp', 0), presets.levelParams('kasp', 1))
    .should.eql({ max_size: 32, min_tm: 55, max_tm: 65, min_gc: 20, max_gc: 80, product_size_ranges: [[65, 150]] });
  presets.paramChanges(presets.levelParams('kasp', 1), presets.levelParams('kasp', 2)).should.eql({ min_tm: 52, max_tm_diff: 6 });
  presets.paramChanges(presets.levelParams('kasp', 0), presets.levelParams('kasp', 0)).should.eql({});
});

// ---- defaults and the §2.9 settings --------------------------------------------------------------------------

test('normalize: kasp defaults and the §2.9 settings echo', function () {
  const req = norm(body({ assay: { type: 'kasp', num_sets: 2 } }));
  req.variant_input.should.equal('id');
  req.variant.should.eql({ id: 'rs871475760', alt: 'A' });
  req.assay.should.eql({ type: 'kasp', orientation: 'both', tails: 'ref_fam_alt_hex', deliberate_mismatch: 'none',
    mismatch_position: 2, num_sets: 2, max_relaxation: 2, neighbour_policy: 'avoid_3p' });
  Object.keys(req.assay).should.eql(['type', 'orientation', 'tails', 'deliberate_mismatch', 'mismatch_position', 'num_sets',
    'max_relaxation', 'neighbour_policy']);
  req.orientations.should.eql(['forward', 'reverse']);
  req.preset.should.equal('kasp');
  req.params.should.eql({ opt_size: 22, min_size: 18, max_size: 30, opt_tm: 60, min_tm: 57, max_tm: 63, min_gc: 30, max_gc: 70,
    max_tm_diff: 3, max_poly_x: 5, product_size_ranges: [[61, 120]] });
  req.pinned.should.eql([]);
  req.ladder.should.eql([
    { level: 1, changes: { max_size: 32, min_tm: 55, max_tm: 65, min_gc: 20, max_gc: 80, product_size_ranges: [[65, 150]] } },
    { level: 2, changes: { min_tm: 52, max_tm_diff: 6 } }
  ]);
  req.floors.should.eql({ as_min_tm: 52, as_min_gc: 15 });
  req.levels.should.have.length(3);
  req.levels[0].should.equal(req.params);
  req.should.containEql({ avoid_repeats: false, repeat_mask_mode: null, template_only: false, label: null });
});

test('normalize: defaults per assay type (as_pcr -> tails none, mismatch auto); rule 6 allows any combination', function () {
  norm(body()).assay.should.containEql({ type: 'kasp', tails: 'ref_fam_alt_hex', deliberate_mismatch: 'none', num_sets: 6 });
  const as = norm(body({ assay: { type: 'as_pcr' } }));
  as.assay.should.eql({ type: 'as_pcr', orientation: 'both', tails: 'none', deliberate_mismatch: 'auto', mismatch_position: 2,
    num_sets: 6, max_relaxation: 2, neighbour_policy: 'avoid_3p' });
  as.params.product_size_ranges.should.eql([[150, 300]]);
  as.ladder[0].changes.should.eql({ max_size: 32, min_tm: 55, max_tm: 65, min_gc: 20, max_gc: 80 });
  norm(body({ assay: { type: 'as_pcr', tails: 'ref_hex_alt_fam' } })).assay.should.containEql({ tails: 'ref_hex_alt_fam', deliberate_mismatch: 'auto' });
  norm(body({ assay: { type: 'kasp', deliberate_mismatch: 'auto', mismatch_position: 3 } })).assay
    .should.containEql({ tails: 'ref_fam_alt_hex', deliberate_mismatch: 'auto', mismatch_position: 3 });
  request.ASSAY_DEFAULTS.kasp.should.containEql({ tails: 'ref_fam_alt_hex', deliberate_mismatch: 'none' });
  request.ASSAY_DEFAULTS.as_pcr.should.containEql({ tails: 'none', deliberate_mismatch: 'auto' });
});

test('normalize: num_sets default and maximum from genotyping config; orientation and max_relaxation', function () {
  norm(body(), cases.cfg({ num_sets_default: 4 })).assay.num_sets.should.equal(4);
  throwsCode(function () { norm(body({ assay: { num_sets: 11 } })); }, 'INVALID_REQUEST', { field: 'assay.num_sets' });
  throwsCode(function () { norm(body({ assay: { num_sets: 0 } })); }, 'INVALID_REQUEST', { field: 'assay.num_sets' });
  norm(body({ assay: { orientation: 'reverse' } })).orientations.should.eql(['reverse']);
  const r0 = norm(body({ assay: { max_relaxation: 0 } }));
  r0.levels.should.have.length(1);
  r0.ladder.should.eql([]);
  norm(body({ assay: { max_relaxation: 1 } })).ladder.map(function (l) { return l.level; }).should.eql([1]);
  throwsCode(function () { norm(body({ assay: { max_relaxation: 3 } })); }, 'INVALID_REQUEST', { field: 'assay.max_relaxation' });
  throwsCode(function () { norm(body({ assay: { type: 'tetra_arms' } })); }, 'INVALID_REQUEST', { field: 'assay.type' });
  throwsCode(function () { norm(body({ assay: { mismatch_position: 1 } })); }, 'INVALID_REQUEST', { field: 'assay.mismatch_position' });
  throwsCode(function () { norm(body({ assay: { neighbour_policy: 'block' } })); }, 'INVALID_REQUEST', { field: 'assay.neighbour_policy' });
  norm(body(), cases.cfg({ as_min_tm: 53, as_min_gc: 18 })).floors.should.eql({ as_min_tm: 53, as_min_gc: 18 });
});

// ---- rule 1: unknown keys --------------------------------------------------------------------------------------

test('rule 1: unknown top-level or nested keys -> INVALID_REQUEST {field}', function () {
  throwsCode(function () { norm(body({ mode: 'region' })); }, 'INVALID_REQUEST', { field: 'mode' });
  throwsCode(function () { norm(body({ variant: { id: 'rs871475760', strand: 1 } })); }, 'INVALID_REQUEST', { field: 'variant.strand' });
  throwsCode(function () { norm(body({ assay: { tail: 'none' } })); }, 'INVALID_REQUEST', { field: 'assay.tail' });
  // Design params outside the genotyping subset, even valid ones, are unknown here.
  throwsCode(function () { norm(body({ params: { num_return: 5 } })); }, 'INVALID_REQUEST', { field: 'params.num_return' });
  throwsCode(function () { norm(body({ params: { max_ns: 0 } })); }, 'INVALID_REQUEST', { field: 'params.max_ns' });
  throwsCode(function () { norm(null); }, 'INVALID_REQUEST');
  throwsCode(function () { norm([]); }, 'INVALID_REQUEST');
  request.PARAM_KEYS.should.have.length(18);
  request.PARAM_KEYS.forEach(function (k) { design.PARAM_SPECS.should.have.property(k); });
});

test('system_name, label, repeat_mask_mode, booleans', function () {
  throwsCode(function () { norm(body({ system_name: undefined })); }, 'INVALID_REQUEST', { field: 'system_name' });
  throwsCode(function () { norm(body({ system_name: 'Sorghum' })); }, 'INVALID_REQUEST', { field: 'system_name' });
  throwsCode(function () { norm(body({ label: 'S1 REF' })); }, 'INVALID_REQUEST', { field: 'label' });
  throwsCode(function () { norm(body({ label: 'x'.repeat(41) })); }, 'INVALID_REQUEST', { field: 'label' });
  norm(body({ label: 'my_marker-1.a' })).label.should.equal('my_marker-1.a');
  // rule 4: repeat_mask_mode without avoid_repeats is ignored, as in /primers/design
  should(norm(body({ repeat_mask_mode: 'three_prime' })).repeat_mask_mode).be.null();
  norm(body({ avoid_repeats: true })).repeat_mask_mode.should.equal('n_mask');
  norm(body({ avoid_repeats: true, repeat_mask_mode: 'three_prime' })).repeat_mask_mode.should.equal('three_prime');
  throwsCode(function () { norm(body({ avoid_repeats: true, repeat_mask_mode: 'soft' })); }, 'INVALID_REQUEST', { field: 'repeat_mask_mode' });
  throwsCode(function () { norm(body({ template_only: 'yes' })); }, 'INVALID_REQUEST', { field: 'template_only' });
  norm(body({ template_only: true })).template_only.should.be.true();
});

// ---- rules 2 and 3: the variant ------------------------------------------------------------------------------------

test('rule 2: id together with manual fields, or neither -> INVALID_VARIANT id_or_manual', function () {
  throwsCode(function () { norm(body({ variant: { id: 'rs871475760', position: 11109 } })); }, 'INVALID_VARIANT', { reason: 'id_or_manual' });
  throwsCode(function () { norm(body({ variant: { id: 'rs871475760', region: '1' } })); }, 'INVALID_VARIANT', { reason: 'id_or_manual' });
  throwsCode(function () { norm(body({ variant: { id: 'rs871475760', ref: 'C' } })); }, 'INVALID_VARIANT', { reason: 'id_or_manual' });
  throwsCode(function () { norm(body({ variant: { region: '1', position: 11109, ref: 'C' } })); }, 'INVALID_VARIANT', { reason: 'id_or_manual' });
  throwsCode(function () { norm(body({ variant: { alt: 'A' } })); }, 'INVALID_VARIANT', { reason: 'id_or_manual' });
  throwsCode(function () { norm(body({ variant: {} })); }, 'INVALID_VARIANT', { reason: 'id_or_manual' });
  throwsCode(function () { norm(body({ variant: undefined })); }, 'INVALID_REQUEST', { field: 'variant' });
  // Real ids with ',' and '*' pass the §2.4 pattern; a leading '.' does not.
  norm(body({ variant: { id: 'tmp_1_13549_TTA_T,*' } })).variant.should.eql({ id: 'tmp_1_13549_TTA_T,*', alt: null });
  throwsCode(function () { norm(body({ variant: { id: '../x' } })); }, 'INVALID_REQUEST', { field: 'variant.id' });
  norm(body({ variant: { id: 'rs5413864115', alt: '-' } })).variant.should.eql({ id: 'rs5413864115', alt: '-' });
  throwsCode(function () { norm(body({ variant: { id: 'rs871475760', alt: 'X' } })); }, 'INVALID_REQUEST', { field: 'variant.alt' });
});

test('rule 3: manual alleles -> INVALID_VARIANT alleles / allele_too_long; * and N wait for rule 12', function () {
  const manual = function (ref, alt) { return body({ variant: { region: '1', position: 11109, ref: ref, alt: alt } }); };
  throwsCode(function () { norm(manual('-', '-')); }, 'INVALID_VARIANT', { reason: 'alleles' });
  throwsCode(function () { norm(manual('c', 'C')); }, 'INVALID_VARIANT', { reason: 'alleles' });
  throwsCode(function () { norm(manual('-', '*')); }, 'INVALID_VARIANT', { reason: 'alleles' });
  throwsCode(function () { norm(manual('CN', '-')); }, 'INVALID_VARIANT', { reason: 'alleles' });
  throwsCode(function () { norm(manual('A', 'X')); }, 'INVALID_VARIANT', { reason: 'alleles' });
  throwsCode(function () { norm(manual('A', 'A'.repeat(51))); }, 'INVALID_VARIANT', { reason: 'allele_too_long' });
  throwsCode(function () { norm(body({ variant: { region: '1', position: 0, ref: 'C', alt: 'A' } })); }, 'INVALID_REQUEST', { field: 'variant.position' });
  throwsCode(function () { norm(body({ variant: { region: '', position: 1, ref: 'C', alt: 'A' } })); }, 'INVALID_REQUEST', { field: 'variant.region' });
  // Both input styles, uppercased; '*' and N alleles are left for UNSUPPORTED_ALLELE after the catalog lookup.
  const ens = norm(manual('a', '-'));
  ens.variant_input.should.equal('manual');
  ens.variant.should.eql({ region: '1', position: 11109, ref: 'A', alt: '-' });
  norm(body({ variant: { region: '1', position: 11282, ref: 'ca', alt: 'c' } })).variant.should.eql({ region: '1', position: 11282, ref: 'CA', alt: 'C' });
  norm(manual('C', 'N')).variant.alt.should.equal('N');
  norm(manual('C', '*')).variant.alt.should.equal('*');
});

// ---- rule 5 and the ladder ---------------------------------------------------------------------------------------

test('pinned params: params.min_tm 58 is never relaxed and settings.pinned lists it', function () {
  const req = norm(body({ params: { min_tm: 58 } }));
  req.pinned.should.eql(['min_tm']);
  req.levels.map(function (p) { return p.min_tm; }).should.eql([58, 58, 58]);
  req.levels.map(function (p) { return p.max_tm; }).should.eql([63, 65, 65]);
  req.ladder.should.eql([
    { level: 1, changes: { max_size: 32, max_tm: 65, min_gc: 20, max_gc: 80, product_size_ranges: [[65, 150]] } },
    { level: 2, changes: { max_tm_diff: 6 } }
  ]);
  const pinned = norm(body({ params: { product_size_ranges: [[80, 200]], salt_divalent: 2.5, max_tm_diff: 3 } }));
  pinned.pinned.should.eql(['max_tm_diff', 'salt_divalent', 'product_size_ranges']);
  pinned.levels.map(function (p) { return p.product_size_ranges; }).should.eql([[[80, 200]], [[80, 200]], [[80, 200]]]);
  pinned.levels.map(function (p) { return p.max_tm_diff; }).should.eql([3, 3, 3]);
  pinned.params.salt_divalent.should.equal(2.5);
});

test('product minimum: max_size 30 -> [[61,120]]; user [[40,90]] with max_size 25 -> [[51,90]]', function () {
  norm(body()).params.product_size_ranges.should.eql([[61, 120]]);
  const req = norm(body({ params: { product_size_ranges: [[40, 90]], max_size: 25 } }));
  req.levels.map(function (p) { return p.product_size_ranges; }).should.eql([[[51, 90]], [[51, 90]], [[51, 90]]]);
  // A pinned range still follows a relaxed max_size.
  norm(body({ params: { product_size_ranges: [[40, 90]] } })).levels.map(function (p) { return p.product_size_ranges; })
    .should.eql([[[61, 90]], [[65, 90]], [[65, 90]]]);
});

test('rule 5: param bounds and cross-field rules -> INVALID_PARAMS, at every level the design may run', function () {
  throwsCode(function () { norm(body({ params: { min_tm: 64 } })); }, 'INVALID_PARAMS', { param: 'min_tm' });
  throwsCode(function () { norm(body({ params: { min_tm: 61 } })); }, 'INVALID_PARAMS', { param: 'opt_tm' });
  throwsCode(function () { norm(body({ params: { opt_size: 40 } })); }, 'INVALID_PARAMS', { param: 'opt_size' });
  throwsCode(function () { norm(body({ params: { gc_clamp: 19, min_size: 18 } })); }, 'INVALID_PARAMS', { param: 'gc_clamp' });
  throwsCode(function () { norm(body({ params: { product_size_ranges: [[100, 90]] } })); }, 'INVALID_PARAMS', { param: 'product_size_ranges' });
  throwsCode(function () { norm(body({ params: { product_size_ranges: [[60, 70], [70, 80], [80, 90], [90, 100], [100, 110]] } })); },
    'INVALID_PARAMS', { param: 'product_size_ranges' });
  throwsCode(function () { norm(body({ params: { product_size_ranges: [[100, 1001]] } })); }, 'INVALID_PARAMS', { param: 'product_size_ranges', index: 0 });
  throwsCode(function () { norm(body({ params: 'fast' })); }, 'INVALID_PARAMS');
  // [[20, 62]] is [[61, 62]] at level 0 but [[65, 62]] once level 1 raises max_size to 32.
  throwsCode(function () { norm(body({ params: { product_size_ranges: [[20, 62]] } })); }, 'INVALID_PARAMS', { param: 'product_size_ranges', level: 1 });
  norm(body({ params: { product_size_ranges: [[20, 62]] }, assay: { max_relaxation: 0 } })).params.product_size_ranges.should.eql([[61, 62]]);
});

test('design.parseParams is exported unchanged', function () {
  design.parseParams.should.be.a.Function();
  design.parseParams({ min_tm: 58, product_size_ranges: [[80, 200]] }).should.eql({ min_tm: 58, product_size_ranges: [[80, 200]] });
  should(design.parseParams(undefined)).eql({});
});
