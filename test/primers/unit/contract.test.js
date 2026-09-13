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
const config = require(path.join(ROOT, 'api/helpers/primers/config'));

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

test('the four /primers operations exist before /{collection}, tagged "Primer design", controller primers, JSON only', function () {
  const paths = Object.keys(api.definition.paths);
  const coll = paths.indexOf('/{collection}');
  const expected = [['/primers/design', 'post', 'designPrimers'], ['/primers/genomes', 'get', 'primerGenomes'],
    ['/primers/check', 'post', 'submitPrimerCheck'], ['/primers/check/{job_id}', 'get', 'getPrimerCheck']];
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
});

test('request definitions are strict (additionalProperties false); response definitions are documented', function () {
  const d = definitions();
  ['PrimerDesignRequest', 'PrimerDesignParams', 'PrimerRegion', 'PrimerCheckRequest'].forEach(function (name) {
    should.exist(d[name], name);
    d[name].additionalProperties.should.equal(false, name);
  });
  d.PrimerInterval.minItems.should.equal(2);
  d.PrimerInterval.maxItems.should.equal(2);
  d.PrimerCheckRequest.properties.params.additionalProperties.should.equal(false);
  d.PrimerCheckRequest.properties.pairs.items.additionalProperties.should.equal(false);
  d.PrimerCheckRequest.properties.pairs.items.properties.expected.additionalProperties.should.equal(false);
  ['PrimerError', 'PrimerDesignResponse', 'PrimerGenomesResponse', 'PrimerCheckJob', 'PrimerCheckResults'].forEach(function (name) {
    should.exist(d[name], name);
  });
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
  if (json && typeof json === 'object' && !Array.isArray(json) && typeof json.path === 'string') {
    return {
      method: String(json.method || 'POST').toUpperCase(),
      url: json.path.indexOf(basePath) === 0 ? json.path.slice(basePath.length) : json.path,
      body: json.body,
      query: json.query,
      expect: json.expect === 'invalid' ? 'invalid' : 'valid'
    };
  }
  let url = null;
  if (/^design/i.test(name)) url = '/primers/design';
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
      }
    }
    checks.should.be.above(0);
    designs.should.be.above(0);
    Array.from(paramKeys).sort().should.eql(Object.keys(checkNormalize.PARAM_RULES).sort());
    t.diagnostic('handler layers: ' + checks + ' check and ' + designs + ' design fixture(s); check params sent: ' + Array.from(paramKeys).sort().join(','));
  });
