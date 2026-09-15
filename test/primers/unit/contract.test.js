'use strict';

// Swagger contract for /primers (spec §A.10, §A.12): api/swagger/swagger.yaml loads and validates with sway,
// the request definitions accept the documented request bodies and every fixture under
// test/primers/fixtures/contract/requests/ (written by gramene-primers `npm run fixtures`), reject legacy
// and malformed bodies, and stay in sync with the handlers' own closed sets (design.js, check/normalize.js).
//
// Fixture files (*.json, searched recursively). Either a wrapper
//   {"method": "POST", "path": "/primers/design", "body": {...}, "query": {...}, "expect": "valid" | "invalid"}
// ("path" with or without the /sorghum_v11 basePath; method defaults to POST, expect to valid), or a bare
// request body whose endpoint comes from the file name (design*.json -> POST /primers/design,
// check*.json -> POST /primers/check) or, failing that, from its shape (a "pairs" array -> check, else design).
// variants-*.json (GET /primers/variants[/{variant_id}]) must be wrappers {method, path, query}; their method defaults to GET.
// check-genotyping-*.json are check bodies with a genotyping block (genotyping spec §7.7): routed by the check* name like any check
// fixture, and run through check/normalize.js up to the catalog stub.
// genotyping-design-*.json bare bodies go to POST /primers/genotyping/design.

require('../../../api/helpers/primers/node_compat');

const fs = require('fs');
const path = require('path');
const { test, before } = require('node:test');
const should = require('should');
const sway = require('sway');

const ROOT = path.resolve(__dirname, '..', '..', '..');
const SWAGGER_YAML = path.join(ROOT, 'api', 'swagger', 'swagger.yaml');
// PRIMERS_CONTRACT_FIXTURES points the fixture test at another directory (e.g. a gramene-primers checkout).
const FIXTURE_DIR = process.env.PRIMERS_CONTRACT_FIXTURES
  ? path.resolve(process.env.PRIMERS_CONTRACT_FIXTURES)
  : path.join(ROOT, 'test', 'primers', 'fixtures', 'contract', 'requests');

const design = require(path.join(ROOT, 'api/helpers/primers/design'));
const checkNormalize = require(path.join(ROOT, 'api/helpers/primers/check/normalize'));
const checkGenotype = require(path.join(ROOT, 'api/helpers/primers/check/genotype'));
const config = require(path.join(ROOT, 'api/helpers/primers/config'));
const grequest = require(path.join(ROOT, 'api/helpers/primers/genotyping/request'));
const gpresets = require(path.join(ROOT, 'api/helpers/primers/genotyping/presets'));

const P2_L = 'GGACAGCTCCACAACATATCAG';
const P2_R = 'GGACATTTGAAGCCCATGGCC';

let api = null;
let basePath = null;

before(async function () {
  api = await sway.create({ definition: SWAGGER_YAML });
  basePath = api.definition.basePath;
});

// sway request validation for (method, url) with a JSON body / query; returns the flat list of error codes.
function validate(method, url, opts) {
  opts = opts || {};
  const full = url.indexOf(basePath) === 0 ? url : basePath + url;
  const op = api.getOperation({ url: full, method: method.toLowerCase() });
  should.exist(op, 'no swagger operation for ' + method + ' ' + full);
  const q = full.indexOf('?') >= 0 ? Object.fromEntries(new URLSearchParams(full.slice(full.indexOf('?') + 1))) : {};
  const headers = opts.headers || { 'content-type': 'application/json' };
  const result = op.validateRequest({ url: full, method: method.toLowerCase(), headers: headers, body: opts.body, query: opts.query || q });
  return { op: op, result: result, codes: codesOf(result.errors) };
}

function codesOf(errors) {
  const out = [];
  (errors || []).forEach(function (e) {
    if (e.code) out.push(e.code);
    if (Array.isArray(e.errors)) codesOf(e.errors).forEach(function (c) { out.push(c); });
  });
  return out;
}

function valid(method, url, body, opts) {
  const v = validate(method, url, Object.assign({ body: body }, opts || {}));
  v.result.errors.should.eql([], method + ' ' + url + ' should be valid: ' + JSON.stringify(v.result.errors).slice(0, 600));
  return v;
}

function invalid(method, url, body, expectCode, opts) {
  const v = validate(method, url, Object.assign({ body: body }, opts || {}));
  v.result.errors.length.should.be.above(0, method + ' ' + url + ' should be rejected: ' + String(JSON.stringify(body)).slice(0, 300));
  if (expectCode) v.codes.should.containEql(expectCode);
  return v;
}

function definitions() {
  return api.definition.definitions;
}

// ---- swagger document ------------------------------------------------------------------------------

test('swagger.yaml validates with sway: no errors or warnings besides the pre-existing duplicate operationId', function () {
  const r = api.validate();
  const errors = r.errors.filter(function (e) { return e.code !== 'DUPLICATE_OPERATIONID'; });
  errors.should.eql([]);
  const primersWarnings = r.warnings.filter(function (w) { return /primers|Primer/.test(JSON.stringify(w.path || [])); });
  primersWarnings.should.eql([]);
});

test('the /primers operations exist before /{collection}, tagged "Primer design", controller primers, JSON only', function () {
  const paths = Object.keys(api.definition.paths);
  const coll = paths.indexOf('/{collection}');
  const expected = [['/primers/design', 'post', 'designPrimers'], ['/primers/genomes', 'get', 'primerGenomes'],
    ['/primers/check', 'post', 'submitPrimerCheck'], ['/primers/check/{job_id}', 'get', 'getPrimerCheck'],
    ['/primers/variants', 'get', 'listPrimerVariants'], ['/primers/variants/{variant_id}', 'get', 'getPrimerVariant'],
    ['/primers/genotyping/design', 'post', 'designGenotypingPrimers']];
  expected.forEach(function (row) {
    const idx = paths.indexOf(row[0]);
    idx.should.be.aboveOrEqual(0, row[0]);
    idx.should.be.below(coll, row[0] + ' must precede /{collection}');
    const p = api.definition.paths[row[0]];
    p['x-swagger-router-controller'].should.equal('primers');
    Object.keys(p).filter(function (k) { return !/^x-/.test(k); }).should.eql([row[1]]);
    const op = p[row[1]];
    op.operationId.should.equal(row[2]);
    op.tags.should.eql(['Primer design']);
    op.consumes.should.eql(['application/json']);
    op.produces.should.eql(['application/json']);
  });
  // the operations the controller implements
  const ctrl = require(path.join(ROOT, 'api/controllers/primers'));
  expected.forEach(function (row) { ctrl[row[2]].should.be.a.Function(); });
  // sway resolves concrete urls to the primers paths, not to /{collection}
  api.getPath({ url: basePath + '/primers/design' }).path.should.equal('/primers/design');
  api.getPath({ url: basePath + '/primers/check/0123456789abcdef0123456789abcdef' }).path.should.equal('/primers/check/{job_id}');
  api.getPath({ url: basePath + '/primers/variants' }).path.should.equal('/primers/variants');
  api.getPath({ url: basePath + '/primers/variants/tmp_1_13549_TTA_T%2C%2A' }).path.should.equal('/primers/variants/{variant_id}');
  api.getPath({ url: basePath + '/primers/genotyping/design' }).path.should.equal('/primers/genotyping/design');
});

test('request definitions are strict (additionalProperties false); response definitions are documented', function () {
  const d = definitions();
  ['PrimerDesignRequest', 'PrimerDesignParams', 'PrimerRegion', 'PrimerCheckRequest',
    'PrimerGenotypingRequest', 'PrimerVariantInput', 'PrimerGenotypingAssay', 'PrimerGenotypingParams'].forEach(function (name) {
    should.exist(d[name], name);
    d[name].additionalProperties.should.equal(false, name);
  });
  d.PrimerGenotypingRequest.required.should.eql(['system_name', 'variant']);
  d.PrimerInterval.minItems.should.equal(2);
  d.PrimerInterval.maxItems.should.equal(2);
  d.PrimerCheckRequest.properties.params.additionalProperties.should.equal(false);
  d.PrimerCheckRequest.properties.pairs.items.additionalProperties.should.equal(false);
  d.PrimerCheckRequest.properties.pairs.items.properties.expected.additionalProperties.should.equal(false);
  ['PrimerError', 'PrimerDesignResponse', 'PrimerGenomesResponse', 'PrimerCheckJob', 'PrimerCheckResults',
    'PrimerVariantList', 'PrimerVariantLookup', 'PrimerVariant'].forEach(function (name) {
    should.exist(d[name], name);
  });
});

// Genotyping spec §2.2, §2.6 and §7.7: every new definition is referenced (no UNUSED_DEFINITION), and the additive
// response fields are documented.
test('variants definitions (§2.6) are all present and referenced; PrimerWarning.details and the genomes variation fields (§2.2)', function () {
  const r = api.validate();
  r.warnings.filter(function (w) { return w.code === 'UNUSED_DEFINITION' && /^Primer/.test(String((w.path || [])[1])); }).should.eql([]);
  const d = definitions();
  ['PrimerVariantRecord', 'PrimerVariantIssue', 'PrimerAlleleSite', 'PrimerVariantZone', 'PrimerVariantVcf', 'PrimerVariantMinimal',
    'PrimerVariantMultiallelic', 'PrimerVariant', 'PrimerVariantSource', 'PrimerVariantList', 'PrimerVariantLookup'].forEach(function (name) {
    should.exist(d[name], name);
    should.exist(d[name].properties, name + ' has properties');
  });
  d.PrimerWarning.properties.details.type.should.equal('object');
  d.PrimerGenomesResponse.properties.variation.properties.should.have.keys('available', 'source', 'release');
  d.PrimerGenomesResponse.properties.variation.properties.source.enum.should.eql(['ensembl']);
  d.PrimerGenomesResponse.properties.genomes.items.properties.has_variation.type.should.equal('boolean');
  d.PrimerVariantIssue.properties.code.enum.should.eql(['REF_MISMATCH', 'STAR_ALLELE', 'ALLELE_TOO_LONG', 'UNSUPPORTED_ALLELE', 'REPEAT_TOO_LONG']);
  // the places M3's normalizer can put a null
  [d.PrimerAlleleSite.properties.alt_maps_to, d.PrimerVariantZone, d.PrimerVariantMultiallelic, d.PrimerVariant.properties.discriminating,
    d.PrimerVariant.properties.consequence, d.PrimerVariantRecord.properties.source].forEach(function (s) { s['x-nullable'].should.equal(true); });
});

// Genotyping spec §2.7-§2.8 and §7.7: the design response definitions are all present, have properties and are referenced,
// and every place the design can put a null is x-nullable.
test('genotyping design definitions (§2.7-§2.8) are present and referenced; the response nulls are x-nullable', function () {
  const r = api.validate();
  r.warnings.filter(function (w) { return w.code === 'UNUSED_DEFINITION' && /^Primer/.test(String((w.path || [])[1])); }).should.eql([]);
  const d = definitions();
  ['PrimerNeighbourHit', 'PrimerLikelihoodAt', 'PrimerDeliberateMismatch', 'PrimerGenotypingDiscrimination', 'PrimerTailedStructures',
    'PrimerOligoTemplateSpan', 'PrimerGenotypingOligo', 'PrimerGenotypingOrderRow', 'PrimerCheckPairFragment', 'PrimerGenotypingSetRef',
    'PrimerGenotypingCheckFragment', 'PrimerGenotypingProduct', 'PrimerGenotypingPairThermo', 'PrimerGenotypingTailedThermo',
    'PrimerGenotypingSetIssue', 'PrimerGenotypingSet', 'PrimerGenotypingExplainSide', 'PrimerGenotypingRejected', 'PrimerGenotypingAttempt',
    'PrimerGenotypingOrientation', 'PrimerGenotypingTemplateFeatures', 'PrimerGenotypingTemplate', 'PrimerKaspMix',
    'PrimerGenotypingAssayEffective', 'PrimerGenotypingNeighbourSummary', 'PrimerGenotypingLadderStep', 'PrimerGenotypingSettings',
    'PrimerGenotypingEngine', 'PrimerGenotypingCheckBlock', 'PrimerGenotypingResponse'].forEach(function (name) {
    should.exist(d[name], name);
    should.exist(d[name].properties, name + ' has properties');
  });
  const post = api.definition.paths['/primers/genotyping/design'].post;
  post.responses['200'].schema.$ref.should.equal('#/definitions/PrimerGenotypingResponse');
  post.parameters[0].schema.$ref.should.equal('#/definitions/PrimerGenotypingRequest');
  ['400', '404', '422', '500', '503', '504'].forEach(function (s) { post.responses[s].schema.$ref.should.equal('#/definitions/PrimerError'); });
  [d.PrimerNeighbourHit.properties.distance_from_3p, d.PrimerDeliberateMismatch, d.PrimerGenotypingDiscrimination, d.PrimerTailedStructures,
    d.PrimerGenotypingOligo.properties.allele, d.PrimerGenotypingOligo.properties.dye, d.PrimerGenotypingOligo.properties.tail_seq,
    d.PrimerGenotypingOligo.properties.matched_tm, d.PrimerGenotypingOligo.properties.primer3_problems, d.PrimerGenotypingOrderRow.properties.allele,
    d.PrimerGenotypingOrderRow.properties.dye, d.PrimerGenotypingOrderRow.properties.tail_seq, d.PrimerGenotypingTailedThermo,
    d.PrimerGenotypingOrientation.properties.reason, d.PrimerGenotypingOrientation.properties.relaxation_level, d.PrimerGenotypingTemplate.properties.mask_source,
    d.PrimerKaspMix, d.PrimerGenotypingEngine.properties.primer3, d.PrimerGenotypingEngine.properties.thermo, d.PrimerGenotypingEngine.properties.variation_source,
    d.PrimerGenotypingCheckBlock, d.PrimerGenotypingResponse.properties.orientations].forEach(function (s) { s['x-nullable'].should.equal(true); });
});

// ---- documented bodies -----------------------------------------------------------------------------

const DESIGN_EXAMPLES = [
  // §A.2.1 examples
  { mode: 'gene', gene_id: 'SORBI_3001G000200', flank_up: 200, flank_down: 100, params: { product_size_ranges: [[300, 800]] } },
  { mode: 'transcript', gene_id: 'SORBI_3001G000200', junction_spanning: true },
  { mode: 'region', system_name: 'sorghum_bicolor', region: { region: '1', start: 11080, end: 15099, strand: -1 }, target: [499, 50], excluded: [[1000, 40]] },
  { mode: 'sequence', sequence: '>amp\nATGGCCRYT...', system_name: 'sorghum_bicolor', avoid_repeats: true },
  // §10.3 rows
  { mode: 'transcript', gene_id: 'SORBI_3001G000200' },
  { mode: 'gene', gene_id: 'SORBI_3001G000700' },
  { mode: 'gene', gene_id: 'SORBI_3004G087700', avoid_repeats: true, template_only: true },
  { mode: 'region', system_name: 'sorghum_tx436pac', region: { region: '4', start: 7547610, end: 7564601 }, avoid_repeats: true },
  {
    mode: 'region', system_name: 'sorghum_bicolor', region: { region: '1', start: 11080, end: 15099, strand: -1 },
    target: [499, 50], included: [100, 3000], excluded: [[1000, 40]],
    params: { min_size: 19, max_size: 22, opt_size: 20, min_tm: 58, max_tm: 61, opt_tm: 59.5, min_gc: 40, max_gc: 60, max_tm_diff: 2, product_size_ranges: [[200, 300]] }
  },
  // every field and every parameter at once
  {
    mode: 'transcript', gene_id: 'SORBI_3001G000200', transcript_id: 'SORBI_3001G000200.1', system_name: 'sorghum_bicolor',
    region: { region: '1', start: 1, end: 2, strand: 1 }, sequence: 'ACGT', flank_up: 0, flank_down: 10000,
    target: [1, 1], included: [1, 100], excluded: [], avoid_repeats: false, repeat_mask_mode: 'three_prime',
    junction_spanning: false, template_only: false,
    params: {
      opt_size: 20, min_size: 18, max_size: 24, opt_tm: 60, min_tm: 58, max_tm: 62, opt_gc: 50, min_gc: 35, max_gc: 65,
      max_tm_diff: 2, max_poly_x: 4, gc_clamp: 1, max_end_stability: 9, max_ns: 0, salt_monovalent: 50, salt_divalent: 1.5,
      dntp_conc: 0.6, dna_conc: 50, num_return: 20, min_3_prime_overlap_of_junction: 4, min_5_prime_overlap_of_junction: 7,
      product_size_ranges: [[70, 150], [150, 300]]
    }
  }
];

const CHECK_EXAMPLES = [
  // §10.4 P2 and the §A.2.3 shape with every field
  { system_name: 'sorghum_bicolor', mode: 'gene', gene_id: 'SORBI_3004G087700', pairs: [{ id: 'P2', left: P2_L, right: P2_R, expected: { region: '4', start: 7423537, end: 7423746 } }] },
  { system_name: 'sorghum_bicolor', pairs: [{ id: 'P2', left: P2_L.toLowerCase(), right: P2_R }] },
  { system_name: 'sorghum_bicolor', mode: 'transcript', gene_id: 'SORBI_3004G087700', pairs: [{ id: 'J', left: 'CCAACAAAGTCATGGATGCACT', right: 'GTGAACATCATGCTGCCCGATG' }] },
  {
    system_name: 'sorghum_bicolor', mode: 'gene', gene_id: 'SORBI_3004G087700', transcript_id: 'SORBI_3004G087700.1',
    checks: ['specificity', 'pangenome'], genomes: ['sorghum_353', 'sorghum_grassl', 'sorghum_leoti'],
    params: { max_product_size: 4000, ignore_mismatches: 6, min_total_mismatches: 2, min_3p_mismatches: 2, three_prime_window: 5, include_unlikely: false, repeat_site_threshold: 5 },
    pairs: [{ id: 'P3', left: 'GATATCAGTGGAATCATAAGACCG', right: 'CATCGATATCAGGATCTGGCTT', expected: { region: '4', start: 7422482, end: 7423061 } }]
  },
  { system_name: 'sorghum_bicolor', mode: 'sequence', checks: ['pangenome'], pairs: Array.from({ length: 10 }, function (_, i) { return { id: 'p.' + i + ':x-y_z', left: 'ACGTACGTACGTACG', right: 'ACGTACGTACGTACGTACGTACGTACGTACGTACGT' }; }) }
];

test('documented design request bodies validate against PrimerDesignRequest', function () {
  DESIGN_EXAMPLES.forEach(function (body) { valid('POST', '/primers/design', body); });
});

test('documented check request bodies validate against PrimerCheckRequest', function () {
  CHECK_EXAMPLES.forEach(function (body) { valid('POST', '/primers/check', body); });
});

// Genotyping spec §2.9-§2.10: the documented design requests, and a body with every field.
const GENOTYPING_EXAMPLES = [
  { system_name: 'sorghum_bicolor', variant: { id: 'rs871475760', alt: 'A' }, assay: { type: 'kasp', num_sets: 2 } },
  { system_name: 'sorghum_bicolor', variant: { id: 'tmp_1_11193_C_T' }, assay: { num_sets: 2 } },
  { system_name: 'sorghum_bicolor', variant: { region: '1', position: 11283, ref: 'A', alt: '-' }, assay: { num_sets: 2 } },
  { system_name: 'sorghum_bicolor', variant: { id: 'tmp_1_11502_C_CGT' }, assay: { num_sets: 1 } },
  { system_name: 'sorghum_bicolor', variant: { id: 'rs871475760' }, assay: { type: 'as_pcr', num_sets: 1 } },
  { system_name: 'sorghum_bicolor', variant: { region: '1', position: 11109, ref: 'A', alt: 'C' }, template_only: true },
  { system_name: 'sorghum_bicolor', variant: { id: 'rs5413863494' } },
  { system_name: 'sorghum_bicolor', variant: { id: 'tmp_1_13549_TTA_T,*', alt: 't' } },
  { system_name: 'sorghum_bicolor', variant: { region: '1', position: 11282, ref: 'ca', alt: 'c' } },
  { system_name: 'sorghum_bicolor', variant: { region: '1', position: 11503, ref: '-', alt: 'GT' } },
  {
    system_name: 'sorghum_bicolor', variant: { region: '1', position: 11109, ref: 'C', alt: 'A' },
    assay: { type: 'as_pcr', orientation: 'reverse', tails: 'ref_hex_alt_fam', deliberate_mismatch: 'auto', mismatch_position: 3, num_sets: 10,
      max_relaxation: 0, neighbour_policy: 'ignore' },
    avoid_repeats: true, repeat_mask_mode: 'three_prime', template_only: false, label: 'my.assay-1',
    params: {
      opt_size: 22, min_size: 18, max_size: 30, opt_tm: 60, min_tm: 57, max_tm: 63, opt_gc: 50, min_gc: 30, max_gc: 70, max_tm_diff: 3,
      max_poly_x: 5, gc_clamp: 1, max_end_stability: 9, salt_monovalent: 50, salt_divalent: 1.5, dntp_conc: 0.6, dna_conc: 50,
      product_size_ranges: [[61, 120], [150, 300]]
    }
  }
];

test('documented genotyping design bodies validate against PrimerGenotypingRequest and pass genotyping/request.normalize', function () {
  const cfg = cfgForTests();
  GENOTYPING_EXAMPLES.forEach(function (body) {
    valid('POST', '/primers/genotyping/design', body);
    (function () { grequest.normalize(JSON.parse(JSON.stringify(body)), cfg); }).should.not.throw(JSON.stringify(body));
  });
});

test('malformed genotyping design bodies are rejected by the validator', function () {
  const base = { system_name: 'sorghum_bicolor', variant: { id: 'rs871475760' } };
  const withVariant = function (v) { return Object.assign({}, base, { variant: v }); };
  const bad = [
    [Object.assign({ mode: 'region' }, base), 'OBJECT_ADDITIONAL_PROPERTIES'],
    [{ system_name: 'sorghum_bicolor' }, 'OBJECT_MISSING_REQUIRED_PROPERTY'],
    [{ variant: { id: 'rs871475760' } }, 'OBJECT_MISSING_REQUIRED_PROPERTY'],
    [Object.assign({}, base, { system_name: '../etc' }), 'PATTERN'],
    [withVariant({ id: 'rs871475760', name: 'x' }), 'OBJECT_ADDITIONAL_PROPERTIES'],
    [withVariant({ id: 'a/b' }), 'PATTERN'],
    [withVariant({ id: '.hidden' }), 'PATTERN'],
    [withVariant({ id: 'a'.repeat(256) }), 'PATTERN'],
    [withVariant({ region: '1', position: 0, ref: 'C', alt: 'A' }), 'MINIMUM'],
    [withVariant({ region: '1', position: 1.5, ref: 'C', alt: 'A' }), 'INVALID_TYPE'],
    [withVariant({ region: '1', position: 11109, ref: 'N', alt: 'A' }), 'PATTERN'],
    [withVariant({ region: '1', position: 11109, ref: '*', alt: 'A' }), 'PATTERN'],
    [withVariant({ region: '1', position: 11109, ref: 'C', alt: 'A'.repeat(51) }), 'PATTERN'],
    [withVariant({ region: 'r'.repeat(256), position: 1, ref: 'C', alt: 'A' }), 'MAX_LENGTH'],
    [Object.assign({}, base, { assay: { type: 'tetra_arms' } }), 'ENUM_MISMATCH'],
    [Object.assign({}, base, { assay: { tails: 'fam' } }), 'ENUM_MISMATCH'],
    [Object.assign({}, base, { assay: { mismatch_position: 4 } }), 'ENUM_MISMATCH'],
    [Object.assign({}, base, { assay: { num_sets: 11 } }), 'MAXIMUM'],
    [Object.assign({}, base, { assay: { num_sets: 0 } }), 'MINIMUM'],
    [Object.assign({}, base, { assay: { max_relaxation: 3 } }), 'MAXIMUM'],
    [Object.assign({}, base, { assay: { neighbour_policy: 'block' } }), 'ENUM_MISMATCH'],
    [Object.assign({}, base, { assay: { dye: 'FAM' } }), 'OBJECT_ADDITIONAL_PROPERTIES'],
    [Object.assign({}, base, { avoid_repeats: 'yes' }), 'INVALID_TYPE'],
    [Object.assign({}, base, { repeat_mask_mode: 'hard' }), 'ENUM_MISMATCH'],
    [Object.assign({}, base, { label: 'my label' }), 'PATTERN'],
    [Object.assign({}, base, { label: 'x'.repeat(41) }), 'PATTERN'],
    [Object.assign({}, base, { params: { num_return: 5 } }), 'OBJECT_ADDITIONAL_PROPERTIES'],
    [Object.assign({}, base, { params: { max_ns: 1 } }), 'OBJECT_ADDITIONAL_PROPERTIES'],
    [Object.assign({}, base, { params: { opt_size: 14 } }), 'MINIMUM'],
    [Object.assign({}, base, { params: { product_size_ranges: [[61, 120], [61, 120], [61, 120], [61, 120], [61, 120]] } }), 'ARRAY_LENGTH_LONG'],
    [Object.assign({}, base, { params: { product_size_ranges: [[61, 1001]] } }), 'MAXIMUM'],
    [Object.assign({}, base, { params: { product_size_ranges: [[19, 120]] } }), 'MINIMUM'],
    [Object.assign({}, base, { variant: null }), 'INVALID_TYPE']
  ];
  bad.forEach(function (row) { invalid('POST', '/primers/genotyping/design', row[0], row[1]); });
  invalid('POST', '/primers/genotyping/design', GENOTYPING_EXAMPLES[0], 'INVALID_CONTENT_TYPE', { headers: { 'content-type': 'text/plain' } });
});

// Genotyping spec §2.11 and §7.7 (M7): PrimerCheckRequest.genotyping, its three strict definitions, the documented bodies (the
// §2.11 example and the design's check.request of §2.9) and malformed genotyping blocks.
const CHECK_GENOTYPING_PAIRS = [
  { id: 'S1_REF', left: 'AGCTTCTCTAAGTGGTTATCCGA', right: 'ATCTTTGACTAGCGAGAAATTCAG', expected: { region: '1', start: 11068, end: 11132 } },
  { id: 'S1_ALT', left: 'AGCTTCTCTAAGTGGTTATCCGA', right: 'ATCTTTGACTAGCGAGAAATTCAT', expected: { region: '1', start: 11068, end: 11132 } },
  { id: 'S2_REF', left: 'GGTTATCCGAATATAGTCATACTCTATTC', right: 'TCTTTGTCTACTGAGAAATCCAGA', expected: { region: '1', start: 11081, end: 11172 } },
  { id: 'S2_ALT', left: 'GGTTATCCGAATATAGTCATACTCTATTA', right: 'TCTTTGTCTACTGAGAAATCCAGA', expected: { region: '1', start: 11081, end: 11172 } }
];
const CHECK_GENOTYPING_BLOCK = {
  variant: { region: '1', position: 11109, ref: 'C', alt: 'A' },
  sets: [{ id: 'S1', ref_pair: 'S1_REF', alt_pair: 'S1_ALT' }, { id: 'S2', ref_pair: 'S2_REF', alt_pair: 'S2_ALT' }]
};
const CHECK_GENOTYPING_EXAMPLES = [
  { system_name: 'sorghum_bicolor', mode: 'region', checks: ['specificity', 'pangenome'], genomes: ['sorghum_bicolorv5', 'sorghum_pi180348', 'sorghum_pi329250'],
    pairs: CHECK_GENOTYPING_PAIRS, genotyping: CHECK_GENOTYPING_BLOCK },
  { system_name: 'sorghum_bicolor', mode: 'region', checks: ['specificity', 'pangenome'], pairs: CHECK_GENOTYPING_PAIRS, genotyping: CHECK_GENOTYPING_BLOCK }
];

test('PrimerCheckRequest.genotyping references three strict definitions (§2.11)', function () {
  const d = definitions();
  d.PrimerCheckRequest.properties.genotyping.should.eql({ $ref: '#/definitions/PrimerCheckGenotyping' });
  ['PrimerCheckGenotyping', 'PrimerCheckGenotypingVariant', 'PrimerCheckGenotypingSet'].forEach(function (name) {
    should.exist(d[name], name);
    d[name].additionalProperties.should.equal(false, name);
  });
  d.PrimerCheckGenotyping.required.should.eql(['variant', 'sets']);
  d.PrimerCheckGenotyping.properties.sets.should.match({ type: 'array', minItems: 1, maxItems: 5 });
  d.PrimerCheckGenotypingVariant.required.should.eql(['region', 'position', 'ref', 'alt']);
  d.PrimerCheckGenotypingSet.required.should.eql(['id', 'ref_pair', 'alt_pair']);
});

test('documented genotyping check bodies validate against PrimerCheckRequest', function () {
  CHECK_GENOTYPING_EXAMPLES.forEach(function (body) { valid('POST', '/primers/check', body); });
  // lowercase alleles, a variant that is not left-aligned and gene mode are the handler's business, not the contract's
  valid('POST', '/primers/check', Object.assign({}, CHECK_GENOTYPING_EXAMPLES[1], { mode: 'gene', gene_id: 'SORBI_3001G000200',
    genotyping: { variant: { region: '1', position: 11284, ref: 'aa', alt: 'a' }, sets: [CHECK_GENOTYPING_BLOCK.sets[0]] } }));
});

test('malformed genotyping blocks are rejected', function () {
  const base = CHECK_GENOTYPING_EXAMPLES[1];
  const g = function (block) { return Object.assign({}, base, { genotyping: block }); };
  const variant = function (patch) { return g({ variant: Object.assign({}, CHECK_GENOTYPING_BLOCK.variant, patch), sets: CHECK_GENOTYPING_BLOCK.sets }); };
  const set0 = function (patch) { return g({ variant: CHECK_GENOTYPING_BLOCK.variant, sets: [Object.assign({}, CHECK_GENOTYPING_BLOCK.sets[0], patch)] }); };
  const bad = [
    [g(null), 'INVALID_TYPE'],
    [g(Object.assign({ orientation: 'reverse' }, CHECK_GENOTYPING_BLOCK)), 'OBJECT_ADDITIONAL_PROPERTIES'],
    [g({ variant: CHECK_GENOTYPING_BLOCK.variant }), 'OBJECT_MISSING_REQUIRED_PROPERTY'],
    [g({ sets: CHECK_GENOTYPING_BLOCK.sets }), 'OBJECT_MISSING_REQUIRED_PROPERTY'],
    [g({ variant: CHECK_GENOTYPING_BLOCK.variant, sets: [] }), 'ARRAY_LENGTH_SHORT'],
    [g({ variant: CHECK_GENOTYPING_BLOCK.variant, sets: Array(6).fill(CHECK_GENOTYPING_BLOCK.sets[0]) }), 'ARRAY_LENGTH_LONG'],
    [set0({ orientation: 'reverse' }), 'OBJECT_ADDITIONAL_PROPERTIES'],
    [g({ variant: CHECK_GENOTYPING_BLOCK.variant, sets: [{ id: 'S1', ref_pair: 'S1_REF' }] }), 'OBJECT_MISSING_REQUIRED_PROPERTY'],
    [set0({ id: 'S 1' }), 'PATTERN'],
    [set0({ ref_pair: 'p'.repeat(65) }), 'PATTERN'],
    [set0({ alt_pair: '' }), 'PATTERN'],
    [variant({ strand: 1 }), 'OBJECT_ADDITIONAL_PROPERTIES'],
    [variant({ position: 0 }), 'MINIMUM'],
    [variant({ position: '11109' }), 'INVALID_TYPE'],
    [variant({ region: 'r'.repeat(256) }), 'MAX_LENGTH'],
    [variant({ ref: 'N' }), 'PATTERN'],
    [variant({ alt: '-' }), 'PATTERN'],
    [variant({ alt: 'A'.repeat(51) }), 'PATTERN']
  ];
  bad.forEach(function (row) { invalid('POST', '/primers/check', row[0], row[1]); });
});

test('legacy and malformed check bodies are rejected', function () {
  const pair = { id: 'P2', left: P2_L, right: P2_R };
  const base = { system_name: 'sorghum_bicolor', pairs: [pair] };
  const bad = [
    [Object.assign({ max_mismatches: 3 }, base), 'OBJECT_ADDITIONAL_PROPERTIES'], // legacy top-level max_mismatches
    [Object.assign({}, base, { checks: ['transcriptome'] }), 'ENUM_MISMATCH'], // legacy check name
    [Object.assign({}, base, { checks: ['specificity', 'transcriptome'] }), 'ENUM_MISMATCH'],
    [Object.assign({}, base, { params: { max_mismatches: 3 } }), 'OBJECT_ADDITIONAL_PROPERTIES'],
    [Object.assign({}, base, { targets: ['genome'] }), 'OBJECT_ADDITIONAL_PROPERTIES'],
    [{ pairs: [pair] }, 'OBJECT_MISSING_REQUIRED_PROPERTY'],
    [{ system_name: 'sorghum_bicolor' }, 'OBJECT_MISSING_REQUIRED_PROPERTY'],
    [Object.assign({}, base, { system_name: '../etc' }), 'PATTERN'],
    [Object.assign({}, base, { mode: 'protein' }), 'ENUM_MISMATCH'],
    [Object.assign({}, base, { checks: [] }), 'ARRAY_LENGTH_SHORT'],
    [Object.assign({}, base, { genomes: ['sorghum_353', 'Sorghum_Rio'] }), 'PATTERN'],
    [Object.assign({}, base, { genomes: Array(151).fill('sorghum_353') }), 'ARRAY_LENGTH_LONG'],
    [Object.assign({}, base, { pairs: [] }), 'ARRAY_LENGTH_SHORT'],
    [Object.assign({}, base, { pairs: Array(11).fill(pair) }), 'ARRAY_LENGTH_LONG'],
    [Object.assign({}, base, { pairs: [{ id: 'P2', left: P2_L }] }), 'OBJECT_MISSING_REQUIRED_PROPERTY'],
    [Object.assign({}, base, { pairs: [Object.assign({ tm: 60 }, pair)] }), 'OBJECT_ADDITIONAL_PROPERTIES'],
    [Object.assign({}, base, { pairs: [Object.assign({}, pair, { left: 'ACGTNACGTACGTACGTACG' })] }), 'PATTERN'],
    [Object.assign({}, base, { pairs: [Object.assign({}, pair, { left: 'ACGTACGTACGTAC' })] }), 'PATTERN'], // 14 nt
    [Object.assign({}, base, { pairs: [Object.assign({}, pair, { right: 'A'.repeat(37) })] }), 'PATTERN'],
    [Object.assign({}, base, { pairs: [Object.assign({}, pair, { id: 'P 2' })] }), 'PATTERN'],
    [Object.assign({}, base, { pairs: [Object.assign({}, pair, { expected: { region: '4', start: 1 } })] }), 'OBJECT_MISSING_REQUIRED_PROPERTY'],
    [Object.assign({}, base, { pairs: [Object.assign({}, pair, { expected: { region: '4', start: 1, end: 2, strand: 1 } })] }), 'OBJECT_ADDITIONAL_PROPERTIES'],
    [Object.assign({}, base, { pairs: [Object.assign({}, pair, { expected: { region: '4', start: 0, end: 2 } })] }), 'MINIMUM'],
    [Object.assign({}, base, { params: { max_product_size: 10001 } }), 'MAXIMUM'],
    [Object.assign({}, base, { params: { ignore_mismatches: 7 } }), 'MAXIMUM'],
    [Object.assign({}, base, { params: { include_unlikely: 'yes' } }), 'INVALID_TYPE'],
    [Object.assign({}, base, { params: { max_amplifying_mismatches: 6 } }), 'MAXIMUM'],
    [Object.assign({}, base, { params: { max_amplifying_mismatches: -1 } }), 'MINIMUM'],
    [Object.assign({}, base, { params: { max_amplifying_mismatches: 2.5 } }), 'INVALID_TYPE'],
    [Object.assign({}, base, { params: { max_amplifying_mismatches: '3' } }), 'INVALID_TYPE']
  ];
  bad.forEach(function (row) { invalid('POST', '/primers/check', row[0], row[1]); });
});

// Stringency cap (plan decision; algorithm version 2): integer 0-5 in the contract; "below ignore_mismatches" is a
// cross-field rule the swagger schema cannot express, so the validator accepts it and the handler refuses it.
test('params.max_amplifying_mismatches: 0-5 passes swagger; a value not below ignore_mismatches is 400 INVALID_PARAMS in the handler', async function () {
  const body = function (params) {
    return { system_name: 'sorghum_bicolor', params: params, pairs: [{ id: 'P2', left: P2_L, right: P2_R }] };
  };
  [0, 1, 2, 3, 4, 5].forEach(function (cap) { valid('POST', '/primers/check', body({ max_amplifying_mismatches: cap, ignore_mismatches: 6 })); });
  valid('POST', '/primers/check', body({ max_amplifying_mismatches: 3, ignore_mismatches: 3 }));
  const def = definitions().PrimerCheckRequest.properties.params.properties.max_amplifying_mismatches;
  def.should.match({ type: 'integer', minimum: 0, maximum: 5 });

  // The handler checks params before any catalog or filesystem work; the stubs make a later stage fail loudly.
  const unreachable = function () { throw Object.assign(new Error('normalize went past the params rules'), { code: 'UNREACHABLE' }); };
  const deps = { cfg: cfgForTests(), genomes: { getCatalog: unreachable }, assemblies: { resolve: unreachable }, mongo: {}, log: { info() {}, warn() {}, error() {}, log() {} } };
  const refusal = async function (params) {
    try {
      await checkNormalize.normalize(body(params), deps);
    } catch (e) {
      return e;
    }
    return null;
  };
  for (const [params, cap, ignore] of [[{ max_amplifying_mismatches: 3, ignore_mismatches: 3 }, 3, 3], [{ max_amplifying_mismatches: 5, ignore_mismatches: 4 }, 5, 4]]) {
    const e = await refusal(params);
    should.exist(e, JSON.stringify(params) + ' must be refused');
    e.should.match({ status: 400, code: 'INVALID_PARAMS' });
    e.details.should.eql({ field: 'params.max_amplifying_mismatches', max_amplifying_mismatches: cap, ignore_mismatches: ignore });
  }
  // a cap below ignore_mismatches, or an omitted cap with ignore_mismatches 3 (default lowered to 2), passes the rule
  for (const params of [{ max_amplifying_mismatches: 2, ignore_mismatches: 3 }, { ignore_mismatches: 3 }, {}]) {
    (await refusal(params)).should.match({ code: 'UNREACHABLE' });
  }
});

test('malformed design bodies are rejected by the validator', function () {
  const bad = [
    [{ mode: 'gene', gene_id: 'X', bogus: 1 }, 'OBJECT_ADDITIONAL_PROPERTIES'],
    [{ gene_id: 'X' }, 'OBJECT_MISSING_REQUIRED_PROPERTY'],
    [{ mode: 'protein' }, 'ENUM_MISMATCH'],
    [{ mode: 'region', system_name: '../etc', region: { region: '1', start: 1, end: 100 } }, 'PATTERN'],
    [{ mode: 'region', system_name: 'sorghum_bicolor', region: { region: '1', start: 1, end: 100, strand: 2 } }, 'ENUM_MISMATCH'],
    [{ mode: 'region', system_name: 'sorghum_bicolor', region: { region: '1', start: 0, end: 100 } }, 'MINIMUM'],
    [{ mode: 'region', system_name: 'sorghum_bicolor', region: { region: '1', start: 1 } }, 'OBJECT_MISSING_REQUIRED_PROPERTY'],
    [{ mode: 'region', system_name: 'sorghum_bicolor', region: { chr: '1', start: 1, end: 2 } }, 'OBJECT_ADDITIONAL_PROPERTIES'],
    [{ mode: 'gene', gene_id: 'X', flank_up: 10001 }, 'MAXIMUM'],
    [{ mode: 'gene', gene_id: 'X', flank_down: -1 }, 'MINIMUM'],
    [{ mode: 'gene', gene_id: 'X', target: [1] }, 'ARRAY_LENGTH_SHORT'],
    [{ mode: 'gene', gene_id: 'X', included: [1, 2, 3] }, 'ARRAY_LENGTH_LONG'],
    [{ mode: 'gene', gene_id: 'X', target: [0, 10] }, 'MINIMUM'],
    [{ mode: 'gene', gene_id: 'X', excluded: Array(51).fill([1, 1]) }, 'ARRAY_LENGTH_LONG'],
    [{ mode: 'sequence', sequence: 'A'.repeat(60001) }, 'MAX_LENGTH'],
    [{ mode: 'gene', gene_id: 'X'.repeat(256) }, 'MAX_LENGTH'],
    [{ mode: 'gene', gene_id: 'X', repeat_mask_mode: 'hard' }, 'ENUM_MISMATCH'],
    [{ mode: 'gene', gene_id: 'X', avoid_repeats: 'true' }, 'INVALID_TYPE'],
    [{ mode: 'gene', gene_id: 'X', gene_id_list: [] }, 'OBJECT_ADDITIONAL_PROPERTIES'],
    [{ mode: 'gene', gene_id: 'X', params: { primer_opt_size: 20 } }, 'OBJECT_ADDITIONAL_PROPERTIES'],
    [{ mode: 'gene', gene_id: 'X', params: { num_return: 21 } }, 'MAXIMUM'],
    [{ mode: 'gene', gene_id: 'X', params: { opt_size: 14 } }, 'MINIMUM'],
    [{ mode: 'gene', gene_id: 'X', params: { max_size: 20.5 } }, 'INVALID_TYPE'],
    [{ mode: 'gene', gene_id: 'X', params: { product_size_ranges: [] } }, 'ARRAY_LENGTH_SHORT'],
    [{ mode: 'gene', gene_id: 'X', params: { product_size_ranges: [[10, 100]] } }, 'MINIMUM'],
    [{ mode: 'gene', gene_id: 'X', params: { product_size_ranges: [[100, 200, 300]] } }, 'ARRAY_LENGTH_LONG'],
    [{ mode: 'gene', gene_id: 'X', params: { product_size_ranges: Array(11).fill([100, 200]) } }, 'ARRAY_LENGTH_LONG'],
    [{ mode: 'gene', gene_id: 'X', sequence: null }, 'INVALID_TYPE']
  ];
  bad.forEach(function (row) { invalid('POST', '/primers/design', row[0], row[1]); });
});

test('POST without Content-Type: application/json is INVALID_CONTENT_TYPE; path and query parameters are validated', function () {
  invalid('POST', '/primers/design', { mode: 'transcript', gene_id: 'SORBI_3001G000200' }, 'INVALID_CONTENT_TYPE', { headers: {} });
  invalid('POST', '/primers/design', { mode: 'transcript', gene_id: 'SORBI_3001G000200' }, 'INVALID_CONTENT_TYPE', { headers: { 'content-type': 'application/x-www-form-urlencoded' } });
  invalid('POST', '/primers/check', CHECK_EXAMPLES[0], 'INVALID_CONTENT_TYPE', { headers: { 'content-type': 'text/plain' } });
  valid('POST', '/primers/design', DESIGN_EXAMPLES[1], { headers: { 'content-type': 'application/json; charset=utf-8' } });

  valid('GET', '/primers/genomes?system_name=sorghum_bicolor', undefined, { headers: {} });
  invalid('GET', '/primers/genomes', undefined, 'REQUIRED', { headers: {}, query: {} });
  invalid('GET', '/primers/genomes?system_name=../etc', undefined, 'PATTERN', { headers: {} });
  invalid('GET', '/primers/genomes?system_name=' + 'a'.repeat(129), undefined, 'MAX_LENGTH', { headers: {} });

  valid('GET', '/primers/check/0123456789abcdef0123456789abcdef', undefined, { headers: {} });
  invalid('GET', '/primers/check/0123456789ABCDEF0123456789ABCDEF', undefined, 'PATTERN', { headers: {} });
  invalid('GET', '/primers/check/0123', undefined, 'PATTERN', { headers: {} });
});

// Genotyping spec §2.3-§2.5 and §7.6: the variant id pattern accepts every real id shape (',' '*' ':' '.', 255
// characters) once URL-decoded, and rejects '/', a leading '.' or '-', and 256 characters.
test('GET /primers/variants and /primers/variants/{variant_id}: query and path parameters as documented', function () {
  const opts = { headers: {} };
  const list = '/primers/variants?system_name=sorghum_bicolor&region=1&start=11180&end=11290';
  valid('GET', list, undefined, opts);
  valid('GET', list + '&types=snv,deletion&include_ems=false&limit=5000', undefined, opts);
  valid('GET', '/primers/variants', undefined, { headers: {}, query: { system_name: 'sorghum_bicolor', region: '1', start: '11180', end: '11290', types: 'insertion' } });
  invalid('GET', '/primers/variants?system_name=sorghum_bicolor&region=1&start=11180', undefined, 'REQUIRED', opts);
  invalid('GET', '/primers/variants?region=1&start=1&end=2', undefined, 'REQUIRED', opts);
  invalid('GET', list.replace('start=11180', 'start=0'), undefined, 'MINIMUM', opts);
  invalid('GET', list.replace('start=11180', 'start=abc'), undefined, 'INVALID_TYPE', opts);
  invalid('GET', list + '&limit=5001', undefined, 'MAXIMUM', opts);
  invalid('GET', list + '&limit=0', undefined, 'MINIMUM', opts);
  invalid('GET', list + '&types=snp', undefined, 'ENUM_MISMATCH', opts);
  invalid('GET', list + '&include_ems=maybe', undefined, 'INVALID_TYPE', opts);
  invalid('GET', list.replace('sorghum_bicolor', '../etc'), undefined, 'PATTERN', opts);
  invalid('GET', list.replace('region=1', 'region=' + 'r'.repeat(256)), undefined, 'MAX_LENGTH', opts);

  const lookup = function (id) { return '/primers/variants/' + id + '?system_name=sorghum_bicolor'; };
  ['rs871475760', 'tmp_1_11502_C_CGT', 'tmp_1_13549_TTA_T%2C%2A', 'tmp_1_13549_TTA_T,*', 'a'.repeat(255), 'X.1:2-3'].forEach(function (id) {
    valid('GET', lookup(id), undefined, opts);
  });
  const op = api.getOperation({ url: basePath + lookup('tmp_1_13549_TTA_T%2C%2A'), method: 'get' });
  op.getParameter('variant_id').getValue({ url: basePath + lookup('tmp_1_13549_TTA_T%2C%2A') }).value.should.equal('tmp_1_13549_TTA_T,*');
  ['a%2Fb', '.hidden', '-x', 'a'.repeat(256), 'a%20b'].forEach(function (id) { invalid('GET', lookup(id), undefined, 'PATTERN', opts); });
  invalid('GET', '/primers/variants/rs871475760', undefined, 'REQUIRED', { headers: {}, query: {} });
});

test('the variants parameters match variation/index.js and config: kinds, limit, id pattern', function () {
  const variation = require(path.join(ROOT, 'api/helpers/primers/variation'));
  const param = function (p, name) { return api.definition.paths[p].get.parameters.find(function (x) { return x.name === name; }); };
  param('/primers/variants', 'types').items.enum.should.eql(variation.KINDS.slice());
  definitions().PrimerVariant.properties.kind.enum.should.eql(variation.KINDS.slice());
  param('/primers/variants', 'limit').maximum.should.equal(config.DEFAULTS.variation.list_limit_max);
  param('/primers/variants/{variant_id}', 'variant_id').pattern.should.equal(variation.ID_PATTERN);
  param('/primers/variants', 'system_name').should.match({ pattern: '^[a-z0-9_]+$', maxLength: 128, required: true });
});

// ---- the swagger definitions and the handlers accept the same closed sets ------------------------------

function cfgForTests() {
  return config._build({ env: {}, fileConfig: {} }).config;
}

test('PrimerDesignParams matches design.js PARAM_SPECS: keys, integer vs number, bounds, product_size_ranges', function () {
  const props = definitions().PrimerDesignParams.properties;
  Object.keys(props).sort().should.eql(Object.keys(design.PARAM_SPECS).sort());
  Object.keys(design.PARAM_SPECS).forEach(function (key) {
    const spec = design.PARAM_SPECS[key];
    const p = props[key];
    if (spec.ranges) {
      p.type.should.equal('array');
      p.minItems.should.equal(spec.minItems);
      p.maxItems.should.equal(spec.maxItems);
      p.items.minItems.should.equal(2);
      p.items.maxItems.should.equal(2);
      p.items.items.type.should.equal('integer');
      p.items.items.minimum.should.equal(spec.min);
      p.items.items.maximum.should.equal(spec.max);
    } else {
      p.type.should.equal(spec.int ? 'integer' : 'number', key);
      p.minimum.should.equal(spec.min, key);
      p.maximum.should.equal(spec.max, key);
    }
  });
  definitions().PrimerDesignRequest.properties.mode.enum.should.eql(design.MODES.slice());
  definitions().PrimerDesignRequest.properties.repeat_mask_mode.enum.should.eql(design.MASK_MODES.slice());
});

test('every PrimerDesignRequest / PrimerRegion property is a field design.normalize knows (no "unknown field")', function () {
  const cfg = cfgForTests();
  const unknownField = function (fn) {
    try {
      fn();
      return null;
    } catch (e) {
      return e && e.code === 'INVALID_REQUEST' && /^unknown field/.test(e.message) ? e.details.field : null;
    }
  };
  Object.keys(definitions().PrimerDesignRequest.properties).forEach(function (k) {
    const body = { mode: 'sequence', sequence: 'ACGTACGTACGTACGTACGTACGT' };
    if (!(k in body)) body[k] = null; // present but absent-valued: only the key check applies
    should(unknownField(function () { design.normalize(body, cfg); })).equal(null, 'design.normalize does not know ' + k);
  });
  unknownField(function () { design.normalize({ mode: 'sequence', sequence: 'ACGT', bogus: 1 }, cfg); }).should.equal('bogus');
  Object.keys(definitions().PrimerRegion.properties).forEach(function (k) {
    const region = { region: '1', start: 1, end: 100 };
    if (!(k in region)) region[k] = null;
    should(unknownField(function () { design.normalize({ mode: 'region', system_name: 'sorghum_bicolor', region: region }, cfg); }))
      .equal(null, 'design.normalize does not know region.' + k);
  });
});

test('PrimerCheckRequest matches check/normalize.js: params rules, modes, checks, and a full body passes both layers', function () {
  const req = definitions().PrimerCheckRequest;
  const params = req.properties.params.properties;
  Object.keys(params).sort().should.eql(Object.keys(checkNormalize.PARAM_RULES).sort());
  Object.keys(checkNormalize.PARAM_RULES).forEach(function (k) {
    const rule = checkNormalize.PARAM_RULES[k];
    params[k].type.should.equal(rule.type, k);
    if (rule.type === 'integer') {
      params[k].minimum.should.equal(rule.min, k);
      params[k].maximum.should.equal(rule.max, k);
    }
  });
  req.properties.mode.enum.should.eql(checkNormalize.MODES.slice());
  req.properties.checks.items.enum.should.eql(checkNormalize.CHECKS.slice());
  const ccfg = cfgForTests().check;
  CHECK_EXAMPLES.forEach(function (body) {
    (function () { checkNormalize.validateShape(JSON.parse(JSON.stringify(body)), ccfg); }).should.not.throw();
  });
  // a key the swagger contract lacks is also refused by the handler
  (function () { checkNormalize.validateShape({ system_name: 'sorghum_bicolor', max_mismatches: 3, pairs: [{ id: 'P', left: P2_L, right: P2_R }] }, ccfg); })
    .should.throw({ code: 'INVALID_REQUEST' });
  // genotyping spec §7.7 (M7): every PrimerCheckRequest property is a handler body key, genotyping included; the documented
  // genotyping bodies pass the handler's shape rules, which refuse a key their definition lacks; the definitions match check/genotype.js
  Object.keys(req.properties).sort().should.eql(checkNormalize.BODY_KEYS.slice().sort());
  checkNormalize.BODY_KEYS.should.containEql('genotyping');
  CHECK_GENOTYPING_EXAMPLES.forEach(function (body) {
    (function () { checkNormalize.validateShape(JSON.parse(JSON.stringify(body)), ccfg); }).should.not.throw();
  });
  (function () { checkNormalize.validateShape(Object.assign({}, CHECK_GENOTYPING_EXAMPLES[1], { genotyping: Object.assign({ orientation: 'reverse' }, CHECK_GENOTYPING_BLOCK) }), ccfg); })
    .should.throw({ code: 'INVALID_REQUEST', details: { field: 'genotyping.orientation' } });
  const gd = definitions();
  gd.PrimerCheckGenotyping.properties.sets.maxItems.should.equal(checkGenotype.MAX_SETS);
  Object.keys(gd.PrimerCheckGenotyping.properties).sort().should.eql(['sets', 'variant']);
  Object.keys(gd.PrimerCheckGenotypingVariant.properties).sort().should.eql(['alt', 'position', 'ref', 'region']);
  Object.keys(gd.PrimerCheckGenotypingSet.properties).sort().should.eql(['alt_pair', 'id', 'ref_pair']);
});

// Genotyping spec §7.7 (sync tests, design half): PrimerGenotypingParams is the closed subset of design.PARAM_SPECS with
// equal bounds (product_size_ranges has its own, genotyping/request.js PRODUCT_RANGES); every request, variant, assay and
// params property is a field genotyping/request.normalize knows; the enums equal the module's tables.
test('PrimerGenotypingRequest matches genotyping/request.js: params subset and bounds, known fields, assay enums, id and label patterns', function () {
  const d = definitions();
  const props = d.PrimerGenotypingParams.properties;
  Object.keys(props).should.eql(grequest.PARAM_KEYS.slice());
  Object.keys(props).forEach(function (key) {
    const spec = design.PARAM_SPECS[key];
    should.exist(spec, 'design.PARAM_SPECS has ' + key);
    const p = props[key];
    if (spec.ranges) {
      p.type.should.equal('array');
      p.minItems.should.equal(1);
      p.maxItems.should.equal(grequest.PRODUCT_RANGES.maxItems);
      p.items.minItems.should.equal(2);
      p.items.maxItems.should.equal(2);
      p.items.items.type.should.equal('integer');
      p.items.items.minimum.should.equal(grequest.PRODUCT_RANGES.min);
      p.items.items.maximum.should.equal(grequest.PRODUCT_RANGES.max);
    } else {
      p.type.should.equal(spec.int ? 'integer' : 'number', key);
      p.minimum.should.equal(spec.min, key);
      p.maximum.should.equal(spec.max, key);
    }
  });
  Object.keys(d.PrimerGenotypingRequest.properties).should.eql(grequest.TOP_LEVEL_FIELDS.slice());
  Object.keys(d.PrimerVariantInput.properties).should.eql(grequest.VARIANT_FIELDS.slice());
  Object.keys(d.PrimerGenotypingAssay.properties).should.eql(grequest.ASSAY_FIELDS.slice());
  Object.keys(grequest.ENUMS).forEach(function (k) {
    d.PrimerGenotypingAssay.properties[k].enum.should.eql(grequest.ENUMS[k].slice(), k);
  });
  d.PrimerGenotypingAssay.properties.num_sets.maximum.should.equal(config.DEFAULTS.genotyping.max_sets);
  d.PrimerGenotypingAssay.properties.max_relaxation.maximum.should.equal(gpresets.MAX_LEVEL);
  d.PrimerGenotypingRequest.properties.repeat_mask_mode.enum.should.eql(design.MASK_MODES.slice());
  d.PrimerVariantInput.properties.id.pattern.should.equal(require(path.join(ROOT, 'api/helpers/primers/variation')).ID_PATTERN);
  d.PrimerVariantInput.properties.ref.pattern.should.equal('^([ACGTacgt]{1,' + config.DEFAULTS.variation.max_allele_length + '}|-)$');
  d.PrimerGenotypingRequest.properties.label.pattern.should.equal('^[A-Za-z0-9_.-]{1,40}$');
  d.PrimerGenotypingAssayEffective.properties.tails.enum.should.eql(grequest.ENUMS.tails.slice());

  // present but absent-valued keys: only the unknown-field rule applies, and it knows every documented property
  const cfg = cfgForTests();
  const unknownField = function (body) {
    try {
      grequest.normalize(body, cfg);
      return null;
    } catch (e) {
      return e && e.code === 'INVALID_REQUEST' && /^unknown field/.test(e.message) ? e.details.field : null;
    }
  };
  const base = function () { return { system_name: 'sorghum_bicolor', variant: { id: 'rs871475760' } }; };
  Object.keys(d.PrimerGenotypingRequest.properties).forEach(function (k) {
    const body = base();
    if (!(k in body)) body[k] = null;
    should(unknownField(body)).equal(null, 'request.normalize does not know ' + k);
  });
  [['variant', d.PrimerVariantInput], ['assay', d.PrimerGenotypingAssay], ['params', d.PrimerGenotypingParams]].forEach(function (row) {
    Object.keys(row[1].properties).forEach(function (k) {
      const body = base();
      body[row[0]] = Object.assign({}, row[0] === 'variant' ? body.variant : {});
      if (!(k in body[row[0]])) body[row[0]][k] = null;
      should(unknownField(body)).equal(null, 'request.normalize does not know ' + row[0] + '.' + k);
    });
  });
  unknownField(Object.assign(base(), { bogus: 1 })).should.equal('bogus');
  unknownField(Object.assign(base(), { assay: { dye: 'FAM' } })).should.equal('assay.dye');
});

// ---- request fixtures from gramene-primers ----------------------------------------------------------------

function listJson(dir) {
  const out = [];
  fs.readdirSync(dir, { withFileTypes: true }).forEach(function (ent) {
    const p = path.join(dir, ent.name);
    if (ent.isDirectory()) listJson(p).forEach(function (x) { out.push(x); });
    else if (ent.isFile() && /\.json$/i.test(ent.name)) out.push(p);
  });
  return out.sort();
}

function fixtureCase(file, json) {
  const name = path.basename(file);
  const variants = /^variants-/i.test(name);
  if (json && typeof json === 'object' && !Array.isArray(json) && typeof json.path === 'string') {
    return {
      method: String(json.method || (variants ? 'GET' : 'POST')).toUpperCase(),
      url: json.path.indexOf(basePath) === 0 ? json.path.slice(basePath.length) : json.path,
      body: json.body,
      query: json.query,
      expect: json.expect === 'invalid' ? 'invalid' : 'valid'
    };
  }
  if (variants) throw new Error(file + ': variants-*.json fixtures must be {method, path, query} wrappers (GET requests have no body)');
  let url = null;
  if (/^genotyping-design/i.test(name)) url = '/primers/genotyping/design';
  else if (/^design/i.test(name)) url = '/primers/design';
  else if (/^check/i.test(name)) url = '/primers/check';
  else if (json && Array.isArray(json.pairs)) url = '/primers/check';
  else if (json && typeof json.mode === 'string') url = '/primers/design';
  if (!url) throw new Error('cannot tell which endpoint ' + file + ' is for; use design*/check* names or a {path, body} wrapper');
  return { method: 'POST', url: url, body: json, query: undefined, expect: 'valid' };
}

const HAVE_FIXTURES = fs.existsSync(FIXTURE_DIR);

test('every request fixture in test/primers/fixtures/contract/requests validates (gramene-primers npm run fixtures)',
  { skip: HAVE_FIXTURES ? false : 'no fixtures yet: ' + path.relative(ROOT, FIXTURE_DIR) + ' does not exist' },
  function (t) {
    const files = listJson(FIXTURE_DIR);
    files.length.should.be.above(0, 'the fixture directory exists but holds no .json files');
    const counts = { valid: 0, invalid: 0 };
    files.forEach(function (file) {
      const json = JSON.parse(fs.readFileSync(file, 'utf8'));
      const c = fixtureCase(file, json);
      const headers = c.method === 'GET' ? {} : { 'content-type': 'application/json' };
      const v = validate(c.method, c.url, { body: c.body, query: c.query, headers: headers });
      const rel = path.relative(FIXTURE_DIR, file);
      if (c.expect === 'valid') {
        v.result.errors.should.eql([], rel + ': ' + JSON.stringify(v.result.errors).slice(0, 600));
      } else {
        v.result.errors.length.should.be.above(0, rel + ' was expected to be rejected');
      }
      counts[c.expect]++;
    });
    t.diagnostic('validated ' + files.length + ' fixture(s): ' + JSON.stringify(counts));
  });

// The regenerated gramene-primers fixtures must also pass the handlers' own layers, and together send every check
// param the server knows (so a server param the package does not send, or sends but swagger lacks, shows up here).
test('request fixtures pass the handler shape rules too, and the check fixtures send every check param (incl. max_amplifying_mismatches)',
  { skip: HAVE_FIXTURES ? false : 'no fixtures yet: ' + path.relative(ROOT, FIXTURE_DIR) + ' does not exist' },
  async function (t) {
    const cfg = cfgForTests();
    const unreachable = function () { throw Object.assign(new Error('past the request rules'), { code: 'UNREACHABLE' }); };
    const deps = { cfg: cfg, genomes: { getCatalog: unreachable }, assemblies: { resolve: unreachable }, mongo: {}, log: { info() {}, warn() {}, error() {}, log() {} } };
    const paramKeys = new Set();
    let checks = 0;
    let designs = 0;
    let genotyping = 0;
    for (const file of listJson(FIXTURE_DIR)) {
      const json = JSON.parse(fs.readFileSync(file, 'utf8'));
      const c = fixtureCase(file, json);
      const rel = path.relative(FIXTURE_DIR, file);
      if (c.expect !== 'valid' || c.method !== 'POST') continue;
      if (c.url === '/primers/check') {
        checks++;
        (function () { checkNormalize.validateShape(JSON.parse(JSON.stringify(c.body)), cfg.check); }).should.not.throw(rel);
        Object.keys(c.body.params || {}).forEach(function (k) { paramKeys.add(k); });
        // every request rule before the catalog lookup (pairs, expected, params incl. the cap < ignore_mismatches rule)
        let err = null;
        try { await checkNormalize.normalize(JSON.parse(JSON.stringify(c.body)), deps); } catch (e) { err = e; }
        should.exist(err, rel);
        err.code.should.equal('UNREACHABLE', rel + ': ' + err.code + ' ' + err.message);
      } else if (c.url === '/primers/design') {
        designs++;
        (function () { design.normalize(JSON.parse(JSON.stringify(c.body)), cfg); }).should.not.throw(rel);
      } else if (c.url === '/primers/genotyping/design') {
        genotyping++;
        // §2.7 rules 1-6 (pure; rules 7-15 need the catalog, Ensembl or the FASTA)
        (function () { grequest.normalize(JSON.parse(JSON.stringify(c.body)), cfg); }).should.not.throw(rel);
      }
    }
    checks.should.be.above(0);
    designs.should.be.above(0);
    Array.from(paramKeys).sort().should.eql(Object.keys(checkNormalize.PARAM_RULES).sort());
    t.diagnostic('handler layers: ' + checks + ' check, ' + designs + ' design and ' + genotyping + ' genotyping design fixture(s); check params sent: ' +
      Array.from(paramKeys).sort().join(','));
  });

// Genotyping spec §7.7 (M7): the check-genotyping-*.json request fixtures exist and pass check/normalize.js's pure genotyping
// rules (shape, gene/region mode, set links) up to the catalog stub. The positional rules need the reference FASTA and are unit-tested
// in check_normalize_genotyping.test.js.
test('check-genotyping-*.json request fixtures pass the genotyping rules of check/normalize.js up to the catalog stub',
  { skip: HAVE_FIXTURES ? false : 'no fixtures yet: ' + path.relative(ROOT, FIXTURE_DIR) + ' does not exist' },
  async function (t) {
    const cfg = cfgForTests();
    const unreachable = function () { throw Object.assign(new Error('past the request rules'), { code: 'UNREACHABLE' }); };
    const deps = { cfg: cfg, genomes: { getCatalog: unreachable }, assemblies: { resolve: unreachable }, mongo: {}, log: { info() {}, warn() {}, error() {}, log() {} } };
    const files = listJson(FIXTURE_DIR).filter(function (f) { return /^check-genotyping-/.test(path.basename(f)); });
    if (!process.env.PRIMERS_CONTRACT_FIXTURES) files.length.should.be.above(0);
    for (const file of files) {
      const body = JSON.parse(fs.readFileSync(file, 'utf8'));
      const rel = path.relative(FIXTURE_DIR, file);
      should.exist(body.genotyping, rel);
      let err = null;
      try { await checkNormalize.normalize(JSON.parse(JSON.stringify(body)), deps); } catch (e) { err = e; }
      should.exist(err, rel);
      err.code.should.equal('UNREACHABLE', rel + ': ' + err.code + ' ' + err.message);
    }
    t.diagnostic('genotyping check fixtures: ' + files.map(function (f) { return path.basename(f); }).join(', '));
  });

// Genotyping spec §2.12, §2.13 and §7.7 (M8): PrimerCheckResults.genotyping and its definitions, in step with check/genotype.js,
// and every field the allele caller writes documented. Responses are not validated by sway, which also ignores x-nullable, so the
// walk below checks declared properties, types, enums and nulls itself.
const GENOTYPING_RESULT_DEFINITIONS = ['PrimerCheckGenotypingResults', 'PrimerGenotypeVariantInfo', 'PrimerGenotypeCopy', 'PrimerGenotypeGenome',
  'PrimerGenotypePrimerCall', 'PrimerGenotypeSetGenome', 'PrimerGenotypeSetSummary', 'PrimerGenotypePairSpecificity', 'PrimerGenotypeReferenceControl',
  'PrimerGenotypeSetResult', 'PrimerGenotypeSummary'];

// Pushes one message per undeclared property, wrong type, value outside an enum, or null without x-nullable.
function undocumentedFields(value, schema, where, errors) {
  let s = schema;
  while (s && s.$ref) s = definitions()[s.$ref.replace('#/definitions/', '')];
  if (value === null) {
    if (s['x-nullable'] !== true) errors.push(where + ': null without x-nullable');
    return errors;
  }
  if (s.enum && s.enum.indexOf(value) < 0) errors.push(where + ': ' + JSON.stringify(value) + ' is not in the enum');
  const types = {
    integer: Number.isInteger(value), number: typeof value === 'number', string: typeof value === 'string', boolean: typeof value === 'boolean',
    array: Array.isArray(value), object: value !== null && typeof value === 'object' && !Array.isArray(value)
  };
  if (s.type && !types[s.type]) errors.push(where + ': not ' + s.type);
  if (s.type === 'array' && Array.isArray(value) && s.items) value.forEach(function (x, i) { undocumentedFields(x, s.items, where + '[' + i + ']', errors); });
  if (s.type === 'object' && types.object && s.properties) {
    Object.keys(value).forEach(function (k) {
      if (!s.properties[k]) errors.push(where + '.' + k + ': not declared');
      else undocumentedFields(value[k], s.properties[k], where + '.' + k, errors);
    });
  }
  return errors;
}

test('PrimerCheckResults.genotyping references the §2.12 definitions; their enums are check/genotype.js\'s', function () {
  const d = definitions();
  const genotype = require(path.join(ROOT, 'api/helpers/primers/check/genotype'));
  d.PrimerCheckResults.properties.genotyping.should.eql({ $ref: '#/definitions/PrimerCheckGenotypingResults' });
  d.PrimerCheckResults.properties.warnings.description.should.match(/WEAK_OFF_TARGETS \{count, max_mismatches, examples/);
  GENOTYPING_RESULT_DEFINITIONS.forEach(function (name) {
    should.exist(d[name], name);
    should.exist(d[name].properties, name + ' has properties');
  });
  api.validate().warnings.filter(function (w) { return w.code === 'UNUSED_DEFINITION' && /^Primer/.test(String((w.path || [])[1])); }).should.eql([]);
  d.PrimerGenotypeGenome.properties.allele.enum.should.eql(genotype.ALLELES);
  d.PrimerGenotypeReferenceControl.properties.allele.enum.should.eql(genotype.ALLELES);
  d.PrimerGenotypeGenome.properties.reason.enum.should.eql(genotype.GENOME_REASONS);
  d.PrimerGenotypeCopy.properties.call.enum.should.eql(genotype.COPY_CALLS);
  d.PrimerGenotypePrimerCall.properties.status.enum.slice().sort().should.eql(genotype.STATUSES.slice().sort());
  d.PrimerGenotypeSetGenome.properties.predicted.enum.should.eql(genotype.PREDICTIONS);
  d.PrimerGenotypeSetGenome.properties.reasons.items.enum.should.eql(genotype.SET_REASONS);
  d.PrimerGenotypeReferenceControl.properties.status.enum.should.eql(genotype.CONTROL_STATUSES);
});

test('results.genotyping as documented: the §2.13 example and blocks written by check/genotype.js use only declared fields, types, enums and nulls', async function () {
  const genotype = require(path.join(ROOT, 'api/helpers/primers/check/genotype'));
  const blast = require(path.join(ROOT, 'api/helpers/primers/check/blast'));
  const classify = require(path.join(ROOT, 'api/helpers/primers/check/classify'));
  const fixtures = path.join(ROOT, 'test/primers/fixtures/check_core/genotype');
  const stubs = require(path.join(fixtures, 'stubs'));
  const P = require(path.join(fixtures, 'products'));
  const schema = { $ref: '#/definitions/PrimerCheckGenotypingResults' };
  undocumentedFields(require(path.join(fixtures, 'results_2_13.json')), schema, 'example', []).should.eql([]);

  const cfg = stubs.makeCfg();
  const prepared = P.prepared(stubs.BODY_2_11, cfg);
  // emptyResults' null specificity, control and reference never reach a client: run.js fills them in before the first flush
  const block = genotype.emptyResults(prepared);
  const reference = P.genome('1', stubs.bases(9000, 15500), 9000);
  const input = function (name, products) {
    return { prepared: prepared, system_name: name, display_name: name, is_reference: name === 'sorghum_bicolor', products: products, cfg: cfg.check, params: classify.DEFAULT_PARAMS };
  };
  const rows = fs.readFileSync(path.join(fixtures, 'megablast_rs871475760.tsv'), 'utf8').split('\n')
    .filter(function (l) { return /^sorghum_pi180348\t/.test(l); }).map(function (l) { return blast.parseMegablastLine(l.slice(l.indexOf('\t') + 1)); });
  genotype.writeResults(block, {
    reference: await genotype.callGenome(input('sorghum_bicolor', prepared.sets.map(function (set) { return P.productsForSet(reference, set, 1); })), P.fetchFrom([reference])),
    genomes: [
      await genotype.callGenome(input('sorghum_pi180348', []), { megablast: async function () { return { status: 'ok', rows: rows, query: genotype.megablastQuery(prepared) }; } }),
      await genotype.callGenome(input('sorghum_is36143', []), { megablast: async function () { return { status: 'budget' }; } }),
      genotype.unavailableEntry(prepared, { system_name: 'sorghum_rio', display_name: 'Rio' }, 'db_unavailable')
    ],
    specificity: prepared.sets.map(function () { return { ref_verdict: 'specific', alt_verdict: 'on_target_missing', off_target_count: 0 }; })
  });
  block.genomes.map(function (x) { return x.allele; }).should.eql(['ref', 'alt', 'missing', 'unavailable']);
  undocumentedFields(block, schema, 'written', []).should.eql([]);
  undocumentedFields({ genotyping: block }, { $ref: '#/definitions/PrimerCheckResults' }, 'results', []).should.eql([]);
});
