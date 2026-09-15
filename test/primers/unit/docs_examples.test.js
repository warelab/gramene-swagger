'use strict';

// The examples of the "Genotyping primers (KASP / allele-specific PCR)" section of docs/primer_design_api.md (genotyping spec
// §7.7, M9):
//   - every request example validates against api/swagger/swagger.yaml with sway, as contract.test.js validates request
//     bodies, and passes the handlers' own pure request rules (genotyping/request.normalize; check/normalize.js up to the
//     catalog stub, including the genotyping link rules of check/genotype.js);
//   - every response example validates against the swagger response definition of its operation and status. The check is
//     structural: types, enums, patterns and bounds of the members shown;
//   - an example that cites a capture is that recorded request or response, whole or as an explicit excerpt of it, so the
//     documented examples stay real.
//
// Examples are tagged by an HTML comment on the line just before their fence (invisible when the markdown is rendered):
//   <!-- example: request GET /primers/variants capture=FILE#/request/path -->           then ```http with "GET <path>"
//   <!-- example: request POST /primers/genotyping/design capture=FILE#/request/body -->  then ```json with the body
//   <!-- example: response POST /primers/genotyping/design 200 capture=FILE#/response --> then ```json with the body
//   <!-- example: response POST /primers/check 202 illustrative -->                      no capture: structure only
// The path is the swagger path template. FILE is a capture under test/primers/fixtures/docs/, {source, request {method,
// path, body?}, status, response}, and the part after '#' is a JSON pointer into it. Every ```json and ```http block of the
// section must be a tagged example, and every capture file must be cited.
//
// Elision: "…" (U+2026) is the explicit elision marker, as in the design spec. A "…": "…" member marks an object excerpt
// (recorded members were left out), a "…" array item an array excerpt (the documented items are an in-order subsequence of
// the recorded ones), and a string "head…tail" a shortened string (the recorded string starts with head and ends with tail).
// An object without the marker must show every recorded member, and an array without it every recorded item. Markers are
// removed before the structural check, and an object marked as an excerpt may lack members its definition requires.

require('../../../api/helpers/primers/node_compat');

const fs = require('fs');
const path = require('path');
const { test, before } = require('node:test');
const should = require('should');
const sway = require('sway');
const swayHelpers = require('sway/lib/helpers');

const ROOT = path.resolve(__dirname, '..', '..', '..');
const SWAGGER_YAML = path.join(ROOT, 'api', 'swagger', 'swagger.yaml');
const DOC = path.join(ROOT, 'docs', 'primer_design_api.md');
const CAPTURES = path.join(ROOT, 'test', 'primers', 'fixtures', 'docs');
const SECTION = '## Genotyping primers (KASP / allele-specific PCR)';
const RESULTS_HEADING = '### Results: `results.genotyping`';
const ELLIPSIS = '…';
const TAG_RE = /^<!-- example: (request|response) (GET|POST) (\/\S*)(?: (\d{3}))?((?: (?:capture=\S+|illustrative))*) -->$/;

const config = require(path.join(ROOT, 'api/helpers/primers/config'));
const grequest = require(path.join(ROOT, 'api/helpers/primers/genotyping/request'));
const checkNormalize = require(path.join(ROOT, 'api/helpers/primers/check/normalize'));
const checkGenotype = require(path.join(ROOT, 'api/helpers/primers/check/genotype'));

// Each endpoint the section documents has at least one request and one response example.
const ENDPOINTS = [
  ['GET', '/primers/genomes'],
  ['GET', '/primers/variants'],
  ['GET', '/primers/variants/{variant_id}'],
  ['POST', '/primers/genotyping/design'],
  ['POST', '/primers/check']
];

let api = null;
let basePath = null;
let doc = null;
let section = null;

before(async function () {
  api = await sway.create({ definition: SWAGGER_YAML });
  basePath = api.definition.basePath;
  doc = fs.readFileSync(DOC, 'utf8');
  section = parseSection(doc);
});

function clone(v) {
  return JSON.parse(JSON.stringify(v));
}

function has(o, k) {
  return Object.prototype.hasOwnProperty.call(o, k);
}

function cfgForTests() {
  return config._build({ env: {}, fileConfig: {} }).config;
}

// ---- the section and its examples ---------------------------------------------------------------------------------

// -> {lines, offset (0-based line of the heading), examples [{kind, method, path, status, capture, illustrative, info, body,
//     line (1-based line of the fence)}], untagged [{info, line}]}
function parseSection(md) {
  const all = md.split('\n');
  const start = all.indexOf(SECTION);
  if (start < 0) throw new Error('docs/primer_design_api.md has no "' + SECTION + '" section');
  let end = all.length;
  let inFence = false;
  for (let i = start + 1; i < all.length; i++) {
    if (/^```/.test(all[i])) inFence = !inFence;
    if (!inFence && /^## /.test(all[i])) {
      end = i;
      break;
    }
  }
  const lines = all.slice(start, end);
  const out = { lines: lines, offset: start, examples: [], untagged: [] };
  for (let i = 0; i < lines.length; i++) {
    const open = /^```(\w*)\s*$/.exec(lines[i]);
    if (!open) continue;
    let j = i + 1;
    while (j < lines.length && !/^```\s*$/.test(lines[j])) j++;
    if (j >= lines.length) throw new Error('unclosed code fence at line ' + (start + i + 1));
    const line = start + i + 1;
    const previous = i > 0 ? lines[i - 1].trim() : '';
    const tag = TAG_RE.exec(previous);
    if (tag) {
      const attrs = tag[5].trim().split(/\s+/).filter(Boolean);
      const capture = attrs.find(function (a) { return a.indexOf('capture=') === 0; });
      out.examples.push({
        kind: tag[1],
        method: tag[2],
        path: tag[3],
        status: tag[4] ? Number(tag[4]) : null,
        capture: capture ? capture.slice('capture='.length) : null,
        illustrative: attrs.indexOf('illustrative') >= 0,
        info: open[1],
        body: lines.slice(i + 1, j).join('\n'),
        line: line
      });
    } else {
      if (/^<!--\s*example/.test(previous)) throw new Error('malformed example tag before line ' + line + ': ' + previous);
      out.untagged.push({ info: open[1], line: line });
    }
    i = j;
  }
  return out;
}

function where(e) {
  return 'docs/primer_design_api.md:' + e.line + ' (' + e.kind + ' ' + e.method + ' ' + e.path + (e.status ? ' ' + e.status : '') + ')';
}

function captureFile(ref) {
  const file = ref.split('#')[0];
  const full = path.join(CAPTURES, file);
  if (path.dirname(full) !== CAPTURES || !fs.existsSync(full)) throw new Error('no capture ' + file + ' in test/primers/fixtures/docs');
  return JSON.parse(fs.readFileSync(full, 'utf8'));
}

function pointerGet(obj, pointer) {
  if (pointer === '' || pointer === undefined) return obj;
  if (pointer[0] !== '/') throw new Error('not a JSON pointer: ' + pointer);
  return pointer.split('/').slice(1).reduce(function (o, part) {
    const key = part.replace(/~1/g, '/').replace(/~0/g, '~');
    if (o === null || typeof o !== 'object' || !has(o, key)) throw new Error('nothing at ' + pointer);
    return o[key];
  }, obj);
}

function captureValue(ref) {
  const hash = ref.indexOf('#');
  return pointerGet(captureFile(ref), hash < 0 ? '' : ref.slice(hash + 1));
}

// The documented value is the recorded one, or an excerpt of it marked with "…" (see the header). Appends to and returns
// `out`, the list of differences.
function mismatches(docValue, recorded, at, out) {
  out = out || [];
  const fail = function (message) {
    out.push(message);
    return out;
  };
  if (typeof docValue === 'string' && docValue !== ELLIPSIS && docValue.indexOf(ELLIPSIS) >= 0) {
    const i = docValue.indexOf(ELLIPSIS);
    const head = docValue.slice(0, i);
    const tail = docValue.slice(i + ELLIPSIS.length);
    if (typeof recorded !== 'string' || recorded.length <= head.length + tail.length || recorded.indexOf(head) !== 0 ||
        recorded.slice(recorded.length - tail.length) !== tail) {
      fail(at + ': "' + docValue + '" does not shorten the recorded ' + JSON.stringify(recorded).slice(0, 80));
    }
    return out;
  }
  if (Array.isArray(docValue)) {
    if (!Array.isArray(recorded)) return fail(at + ': documented an array, recorded ' + JSON.stringify(recorded).slice(0, 80));
    const excerpt = docValue.indexOf(ELLIPSIS) >= 0;
    const items = docValue.filter(function (x) { return x !== ELLIPSIS; });
    if (!excerpt) {
      if (items.length !== recorded.length) {
        return fail(at + ': ' + items.length + ' items documented, ' + recorded.length + ' recorded (mark an excerpt with "…")');
      }
      items.forEach(function (item, i) { mismatches(item, recorded[i], at + '[' + i + ']', out); });
      return out;
    }
    // an excerpt: each documented item matches a later recorded item than the one before it
    let j = 0;
    items.forEach(function (item, i) {
      while (j < recorded.length && mismatches(item, recorded[j], at, []).length > 0) j++;
      if (j >= recorded.length) fail(at + '[' + i + ']: no later recorded item is ' + JSON.stringify(item).slice(0, 120));
      j++;
    });
    return out;
  }
  if (docValue !== null && typeof docValue === 'object') {
    if (recorded === null || typeof recorded !== 'object' || Array.isArray(recorded)) {
      return fail(at + ': documented an object, recorded ' + JSON.stringify(recorded).slice(0, 80));
    }
    Object.keys(docValue).filter(function (k) { return k !== ELLIPSIS; }).forEach(function (k) {
      if (!has(recorded, k)) fail(at + '.' + k + ': not in the recording');
      else mismatches(docValue[k], recorded[k], at + '.' + k, out);
    });
    if (!has(docValue, ELLIPSIS)) {
      Object.keys(recorded).forEach(function (k) {
        if (!has(docValue, k)) fail(at + '.' + k + ': recorded but not documented (mark an excerpt with "…": "…")');
      });
    }
    return out;
  }
  if (docValue !== recorded) fail(at + ': documented ' + JSON.stringify(docValue) + ', recorded ' + JSON.stringify(recorded));
  return out;
}

function stripElisions(v) {
  if (Array.isArray(v)) return v.filter(function (x) { return x !== ELLIPSIS; }).map(stripElisions);
  if (v !== null && typeof v === 'object') {
    const out = {};
    Object.keys(v).forEach(function (k) { if (k !== ELLIPSIS) out[k] = stripElisions(v[k]); });
    return out;
  }
  return v;
}

// The paths ('/', '/check/request', '/sets/0', ...) of the objects marked "…": "…", as they are after stripElisions and in
// the form schemaErrors reports.
function excerptPaths(v, at, out) {
  out = out || new Set();
  at = at || '';
  if (Array.isArray(v)) {
    v.filter(function (x) { return x !== ELLIPSIS; }).forEach(function (x, i) { excerptPaths(x, at + '/' + i, out); });
  } else if (v !== null && typeof v === 'object') {
    if (has(v, ELLIPSIS)) out.add(at === '' ? '/' : at);
    Object.keys(v).forEach(function (k) { if (k !== ELLIPSIS) excerptPaths(v[k], at + '/' + k, out); });
  }
  return out;
}

// ---- swagger ------------------------------------------------------------------------------------------------------

function operation(method, template) {
  const op = api.getOperations().find(function (o) { return o.pathObject.path === template && o.method === method.toLowerCase(); });
  should.exist(op, 'swagger.yaml has no operation ' + method + ' ' + template);
  return op;
}

function responseSchema(method, template, status) {
  const res = operation(method, template).getResponse(String(status));
  should.exist(res, method + ' ' + template + ' declares no ' + status + ' response');
  String(res.statusCode).should.equal(String(status), method + ' ' + template + ' declares no ' + status + ' response (only default)');
  const schema = res.definitionFullyResolved.schema;
  should.exist(schema, method + ' ' + template + ' ' + status + ' has no response schema');
  return schema;
}

// sway 1.0.0 validates with z-schema, and neither knows the swagger 2.0 vendor extension x-nullable, so a plain validation of
// a real response rejects every null the definitions allow (multiallelic, alt_maps_to, dye, deliberate_mismatch, ...). This
// copy of the resolved schema accepts null exactly where a schema node says x-nullable: true, by widening its type to
// [type, 'null'] and adding null to its enum. Every other node, and every non-null value, is validated unchanged.
function acceptNull(schema, memo) {
  if (schema === null || typeof schema !== 'object' || Array.isArray(schema)) return schema;
  if (memo.has(schema)) return memo.get(schema);
  const out = Object.assign({}, schema);
  memo.set(schema, out);
  if (schema.properties && typeof schema.properties === 'object') {
    out.properties = {};
    Object.keys(schema.properties).forEach(function (k) { out.properties[k] = acceptNull(schema.properties[k], memo); });
  }
  if (Array.isArray(schema.items)) out.items = schema.items.map(function (s) { return acceptNull(s, memo); });
  else if (schema.items && typeof schema.items === 'object') out.items = acceptNull(schema.items, memo);
  if (schema.additionalProperties && typeof schema.additionalProperties === 'object') out.additionalProperties = acceptNull(schema.additionalProperties, memo);
  ['allOf', 'anyOf', 'oneOf'].forEach(function (k) {
    if (Array.isArray(schema[k])) out[k] = schema[k].map(function (s) { return acceptNull(s, memo); });
  });
  if (schema['x-nullable'] === true) {
    if (typeof schema.type === 'string') out.type = [schema.type, 'null'];
    if (Array.isArray(schema.enum) && schema.enum.indexOf(null) < 0) out.enum = schema.enum.concat([null]);
  }
  return out;
}

function schemaErrors(schema, value) {
  return swayHelpers.validateAgainstSchema(swayHelpers.getJSONSchemaValidator(), schema, value).errors.map(function (e) {
    return { code: e.code, path: '/' + (e.path || []).join('/'), message: e.message };
  });
}

function responseErrors(method, template, status, value) {
  return schemaErrors(acceptNull(responseSchema(method, template, status), new WeakMap()), value);
}

// The structural errors of a documented response body (with its "…" markers). An object marked as an excerpt may leave out
// members its definition requires (the check request inside a design response requires pairs); only that error, and only at
// such an object, is ignored.
function documentedResponseErrors(method, template, status, documented) {
  const excerpts = excerptPaths(documented);
  return responseErrors(method, template, status, stripElisions(documented)).filter(function (err) {
    return !(err.code === 'OBJECT_MISSING_REQUIRED_PROPERTY' && excerpts.has(err.path));
  });
}

async function normalizeToCatalogStub(body, cfg) {
  const unreachable = function () { throw Object.assign(new Error('past the request rules'), { code: 'UNREACHABLE' }); };
  const deps = { cfg: cfg, genomes: { getCatalog: unreachable }, assemblies: { resolve: unreachable }, mongo: {}, log: { info() {}, warn() {}, error() {}, log() {} } };
  try {
    await checkNormalize.normalize(clone(body), deps);
  } catch (e) {
    return e;
  }
  return { code: 'NO_ERROR', message: 'normalize returned before the catalog lookup' };
}

// ---- tests --------------------------------------------------------------------------------------------------------

test('the genotyping section tags every json and http block as an example, for every endpoint it documents', function (t) {
  const loose = section.untagged.filter(function (u) { return u.info === 'json' || u.info === 'http'; });
  loose.should.eql([], 'untagged json/http blocks: ' + loose.map(function (u) { return 'line ' + u.line; }).join(', '));
  section.examples.length.should.be.above(0);
  ENDPOINTS.forEach(function (row) {
    ['request', 'response'].forEach(function (kind) {
      section.examples.some(function (e) { return e.kind === kind && e.method === row[0] && e.path === row[1]; })
        .should.equal(true, 'no ' + kind + ' example for ' + row[0] + ' ' + row[1]);
    });
  });
  section.examples.forEach(function (e) {
    if (e.kind === 'response') should.exist(e.status, where(e) + ': a response example needs a status');
    else should(e.status).equal(null, where(e) + ': a request example has no status');
    (e.capture !== null || e.illustrative).should.equal(true, where(e) + ': cite capture=FILE#POINTER or mark the example illustrative');
    e.info.should.equal(e.kind === 'request' && e.method === 'GET' ? 'http' : 'json', where(e) + ': fence type');
  });
  section.lines.indexOf(RESULTS_HEADING).should.be.aboveOrEqual(0, 'the section has a "' + RESULTS_HEADING + '" subsection');
  t.diagnostic(section.examples.length + ' examples: ' + section.examples.filter(function (e) { return e.kind === 'request'; }).length + ' requests');
});

test('every request example validates against swagger.yaml (sway) and passes the handlers\' pure request rules', async function (t) {
  const cfg = cfgForTests();
  const requests = section.examples.filter(function (e) { return e.kind === 'request'; });
  for (const e of requests) {
    let url;
    let body;
    if (e.method === 'GET') {
      const m = /^GET (\/\S+)$/.exec(e.body.trim());
      should.exist(m, where(e) + ': an http example is the single line "GET <path>"');
      url = m[1];
    } else {
      url = basePath + e.path;
      body = JSON.parse(e.body);
    }
    url.indexOf(basePath + '/').should.equal(0, where(e) + ': the request path starts with the base path ' + basePath);
    const op = api.getOperation({ url: url, method: e.method.toLowerCase() });
    should.exist(op, where(e) + ': no swagger operation for ' + url);
    op.pathObject.path.should.equal(e.path, where(e) + ': the request is for ' + op.pathObject.path);
    const query = url.indexOf('?') >= 0 ? Object.fromEntries(new URLSearchParams(url.slice(url.indexOf('?') + 1))) : {};
    const headers = e.method === 'GET' ? {} : { 'content-type': 'application/json' };
    const result = op.validateRequest({ url: url, method: e.method.toLowerCase(), headers: headers, body: body, query: query });
    result.errors.should.eql([], where(e) + ': ' + JSON.stringify(result.errors).slice(0, 600));

    if (e.path === '/primers/genotyping/design') {
      (function () { grequest.normalize(clone(body), cfg); }).should.not.throw(where(e));
    } else if (e.path === '/primers/check') {
      (function () { checkNormalize.validateShape(clone(body), cfg.check); }).should.not.throw(where(e));
      if (body.genotyping) (function () { checkGenotype.validateLinks(body.genotyping, body.pairs); }).should.not.throw(where(e));
      const err = await normalizeToCatalogStub(body, cfg);
      err.code.should.equal('UNREACHABLE', where(e) + ': ' + err.code + ' ' + err.message);
    }

    if (e.capture) {
      const recorded = captureValue(e.capture);
      if (e.method === 'GET') recorded.should.equal(url, where(e) + ': the documented request line is not the recorded one');
      else mismatches(body, recorded, where(e)).should.eql([]);
    }
  }
  t.diagnostic(requests.length + ' request examples validated');
});

test('every response example matches the swagger response definition of its operation and status (x-nullable-aware)', function (t) {
  const responses = section.examples.filter(function (e) { return e.kind === 'response'; });
  responses.forEach(function (e) {
    const errors = documentedResponseErrors(e.method, e.path, e.status, JSON.parse(e.body));
    errors.should.eql([], where(e) + ': ' + JSON.stringify(errors).slice(0, 800));
  });
  t.diagnostic(responses.length + ' response examples validated');
});

test('the structural check is not vacuous: nulls pass only where x-nullable, and wrong values fail', function () {
  const full = section.examples.find(function (e) {
    return e.kind === 'response' && e.path === '/primers/genotyping/design' && e.status === 200 && excerptPaths(JSON.parse(e.body)).size === 0;
  });
  should.exist(full, 'a complete (not excerpted) 200 genotyping design example');
  const body = stripElisions(JSON.parse(full.body));
  const raw = responseSchema('POST', '/primers/genotyping/design', 200);
  // without the x-nullable rewrite, sway's validator rejects the documented nulls
  const plain = schemaErrors(raw, body);
  plain.length.should.be.above(0);
  plain.every(function (err) { return err.code === 'INVALID_TYPE' || err.code === 'ENUM_MISMATCH'; }).should.equal(true, JSON.stringify(plain).slice(0, 400));
  responseErrors('POST', '/primers/genotyping/design', 200, body).should.eql([]);
  const broken = [
    function (b) { b.sets[0].score = null; }, // not x-nullable
    function (b) { b.variant.shift = null; },
    function (b) { b.sets[0].quality = 'excellent'; },
    function (b) { b.sets[0].key = 'F9DF650AD116'; },
    function (b) { b.sets[0].primers.as_ref.dye = 'ROX'; },
    function (b) { b.orientations.forward.reason = 'unknown_reason'; },
    function (b) { b.variant.discriminating.forward.alt_maps_to = '11109'; },
    function (b) { delete b.check.request.pairs; } // required, and this object is not marked as an excerpt
  ];
  broken.forEach(function (breakIt, i) {
    const b = clone(body);
    breakIt(b);
    documentedResponseErrors('POST', '/primers/genotyping/design', 200, b).length.should.be.above(0, 'mutation ' + i + ' must be rejected');
  });
  // the same missing member is accepted only in an object marked as an excerpt
  const excerpt = clone(body);
  delete excerpt.check.request.pairs;
  excerpt.check.request[ELLIPSIS] = ELLIPSIS;
  documentedResponseErrors('POST', '/primers/genotyping/design', 200, excerpt).should.eql([]);
});

test('every example that cites a capture is the recorded request or response, or an explicit excerpt of it', function (t) {
  const cited = section.examples.filter(function (e) { return e.capture !== null; });
  cited.forEach(function (e) {
    const file = captureFile(e.capture);
    ['source', 'request', 'status', 'response'].forEach(function (k) { should.exist(file[k], e.capture + ' has ' + k); });
    const pointer = e.capture.indexOf('#') >= 0 ? e.capture.slice(e.capture.indexOf('#') + 1) : '';
    if (pointer === '/response') {
      file.status.should.equal(e.status, where(e) + ': the capture recorded status ' + file.status);
      file.request.method.should.equal(e.method, where(e));
      const recordedPath = file.request.path.split('?')[0];
      const op = api.getOperation({ url: recordedPath, method: e.method.toLowerCase() });
      should.exist(op, where(e) + ': ' + recordedPath);
      op.pathObject.path.should.equal(e.path, where(e) + ': the capture is for another operation');
    }
    if (e.kind === 'response') mismatches(JSON.parse(e.body), captureValue(e.capture), where(e)).should.eql([]);
  });
  const citedFiles = new Set(cited.map(function (e) { return e.capture.split('#')[0]; }));
  fs.readdirSync(CAPTURES).filter(function (f) { return /\.json$/.test(f); }).forEach(function (f) {
    citedFiles.has(f).should.equal(true, 'test/primers/fixtures/docs/' + f + ' is not cited by any example');
  });
  t.diagnostic(cited.length + ' examples compared with ' + citedFiles.size + ' captures');
});

test('the elision rules of the capture comparison', function () {
  mismatches({ a: 1, b: [1, 2, 3], s: 'ACGTACGT' }, { a: 1, b: [1, 2, 3], s: 'ACGTACGT' }, '$').should.eql([]);
  mismatches({ a: 1, [ELLIPSIS]: ELLIPSIS }, { a: 1, b: 2 }, '$').should.eql([]);
  mismatches({ a: 1 }, { a: 1, b: 2 }, '$').length.should.equal(1); // an unmarked excerpt
  mismatches({ o: { a: 1 } }, { o: { a: 1, b: 2 } }, '$').length.should.equal(1); // nested
  mismatches({ b: [1, 3, ELLIPSIS] }, { b: [1, 2, 3, 4] }, '$').should.eql([]);
  mismatches({ b: [3, 1, ELLIPSIS] }, { b: [1, 2, 3, 4] }, '$').length.should.equal(1); // out of order
  mismatches({ b: [1, 2] }, { b: [1, 2, 3] }, '$').length.should.equal(1);
  mismatches({ b: [{ x: 1 }] }, { b: [{ x: 2 }] }, '$').length.should.equal(1);
  mismatches({ b: 'x' }, { b: ['x'] }, '$').length.should.equal(1);
  mismatches({ s: 'ACG' + ELLIPSIS + 'CGT' }, { s: 'ACGTTTACGT' }, '$').should.eql([]);
  mismatches({ s: 'ACG' + ELLIPSIS + 'AAA' }, { s: 'ACGTTTACGT' }, '$').length.should.equal(1);
  mismatches({ n: 57.1 }, { n: 57.10000001 }, '$').length.should.equal(1);
  mismatches({ n: null }, { n: 0 }, '$').length.should.equal(1);
  stripElisions({ a: [1, ELLIPSIS], [ELLIPSIS]: ELLIPSIS, s: 'A' + ELLIPSIS + 'B' }).should.eql({ a: [1], s: 'A' + ELLIPSIS + 'B' });
  Array.from(excerptPaths({ [ELLIPSIS]: ELLIPSIS, a: [ELLIPSIS, { b: { [ELLIPSIS]: ELLIPSIS } }], c: {} })).sort().should.eql(['/', '/a/0/b']);
});

test('the configuration section documents every environment override and every genotyping config key', function () {
  const start = doc.indexOf('\n### Configuration\n');
  start.should.be.above(0, 'docs/primer_design_api.md has an Operations "### Configuration" subsection');
  const end = doc.indexOf('\n### ', start + 1);
  const text = doc.slice(start, end < 0 ? doc.length : end);
  config.ENV_OVERRIDES.map(function (row) { return row[0]; }).forEach(function (name) {
    text.should.containEql('`' + name + '`', 'environment variable ' + name + ' is documented');
  });
  const d = config.DEFAULTS;
  const keys = ['ntthal']
    .concat(Object.keys(d.variation).map(function (k) { return 'variation.' + k; }))
    .concat(Object.keys(d.genotyping).map(function (k) { return 'genotyping.' + k; }))
    .concat(Object.keys(d.check).filter(function (k) { return /^genotype_/.test(k); }).map(function (k) { return 'check.' + k; }));
  keys.forEach(function (k) { text.should.containEql('`' + k + '`', 'config key ' + k + ' is documented'); });
});
