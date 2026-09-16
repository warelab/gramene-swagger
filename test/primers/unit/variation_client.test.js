'use strict';

// variation/client.js (genotyping spec §3.3 and the variation_client rows of §7.3): URL building, record validation,
// status mapping, caches, the breaker, the limiter and single flight, all against an injected fake fetch and clock.
// Never the network. The 25 kb chunking row drives variation/index.js listVariants over this client.

require('../../../api/helpers/primers/node_compat');

const test = require('node:test');
const should = require('should');
const fs = require('fs');
const path = require('path');

const config = require('../../../api/helpers/primers/config');
const errors = require('../../../api/helpers/primers/errors');
const { createVariationClient, VariationError, _internal } = require('../../../api/helpers/primers/variation/client');
const variation = require('../../../api/helpers/primers/variation');

const { PrimerHttpError } = errors;
const FIX = path.join(__dirname, '..', 'fixtures', 'variation');
const BASE = 'https://data.gramene.org/pansite-ensembl-115';
const SPECIES = 'sorghum_bicolor';
const CHR1_LENGTH = 80884392;

function fixture(name) {
  return JSON.parse(fs.readFileSync(path.join(FIX, name), 'utf8'));
}

function cfgWith(variationOverrides) {
  return config._build({ env: {}, fileConfig: {}, overrides: { variation: variationOverrides || {} } }).config;
}

function quietLog() {
  const lines = [];
  const push = function (m) { lines.push(String(m)); };
  return { lines: lines, log: push, info: push, warn: push, error: push };
}

function clock() {
  const c = { t: 1700000000000 };
  c.now = function () { return c.t; };
  c.advance = function (ms) { c.t += ms; };
  return c;
}

function respond(body, status, headers) {
  return new Response(typeof body === 'string' ? body : JSON.stringify(body), { status: status || 200, headers: headers });
}

function refused() {
  const cause = Object.assign(new Error('connect ECONNREFUSED 127.0.0.1:50199'), { code: 'ECONNREFUSED' });
  return Object.assign(new TypeError('fetch failed'), { cause: cause });
}

// handler(url, init, n) -> Response (or throws); every call is recorded.
function make(opts) {
  opts = opts || {};
  const c = opts.clock || clock();
  const log = quietLog();
  const calls = [];
  const handler = opts.handler || function () { return respond([]); };
  const fetch = async function (url, init) {
    calls.push({ url: url, init: init });
    return handler(url, init, calls.length);
  };
  const client = createVariationClient({ cfg: cfgWith(opts.variation), fetch: fetch, now: c.now, log: log });
  return { client: client, calls: calls, clock: c, log: log };
}

async function rejection(promise) {
  try {
    await promise;
  } catch (e) {
    return e;
  }
  throw new Error('expected a rejection');
}

function deferred() {
  const d = {};
  d.promise = new Promise(function (resolve, reject) { d.resolve = resolve; d.reject = reject; });
  return d;
}

function record(over) {
  return Object.assign({
    feature_type: 'variation', seq_region_name: '1', start: 11109, end: 11109, strand: 1, id: 'rs871475760', source: 'EVA',
    alleles: ['C', 'A'], consequence_type: 'downstream_gene_variant', assembly_name: 'Sorghum_bicolor_NCBIv3', clinical_significance: []
  }, over || {});
}

function unavailableWith(e, reason, retryAfterS) {
  e.should.be.instanceOf(VariationError);
  e.should.be.instanceOf(PrimerHttpError);
  e.status.should.equal(503);
  e.code.should.equal('VARIATION_SOURCE_UNAVAILABLE');
  e.details.reason.should.equal(reason);
  if (retryAfterS !== undefined) e.details.retry_after_s.should.equal(retryAfterS);
}

// ---- URLs and records ---------------------------------------------------------------------------------------

test('overlap: URL, transport options, compact frozen records; malformed records are dropped and counted with their position', async function () {
  const recorded = fixture('overlap_1_10709-11509.json');
  const malformed = [
    'a string',
    record({ feature_type: 'gene' }),
    record({ seq_region_name: '2', start: 11200, end: 11200 }),
    record({ start: 'x' }),
    record({ start: 11300, end: 11300, alleles: ['C'] }),
    record({ start: 11301, end: 11301, alleles: ['C', 'Q'] }),
    record({ start: 11302, end: 11302, alleles: ['A', 'C', 'G', 'T', 'AC', 'AG', 'AT', 'CA', 'CG', 'CT', 'GA'] })
  ];
  const w = make({ handler: function () { return respond(recorded.concat(malformed, [record({ id: 'long_source', source: 'S'.repeat(150), consequence_type: 'c'.repeat(101) })])); } });
  const got = await w.client.overlapChunk(SPECIES, '1', 10001, { regionLength: CHR1_LENGTH });

  w.calls.should.have.length(1);
  w.calls[0].url.should.equal(BASE + '/overlap/region/sorghum_bicolor/1:10001-20000?feature=variation;content-type=application/json');
  const init = w.calls[0].init;
  init.method.should.equal('GET');
  init.redirect.should.equal('error');
  init.headers.should.eql({ Accept: 'application/json', 'User-Agent': 'gramene-swagger-primers' });
  init.signal.should.be.instanceOf(AbortSignal);
  Object.keys(init).sort().should.eql(['headers', 'method', 'redirect', 'signal']);

  got.records.should.have.length(recorded.length + 1);
  got.records[0].should.eql({
    id: recorded[0].id, seq_region_name: '1', start: recorded[0].start, end: recorded[0].end, alleles: recorded[0].alleles,
    source: recorded[0].source, consequence_type: recorded[0].consequence_type
  });
  const long = got.records[got.records.length - 1];
  long.source.should.have.length(100);
  long.consequence_type.should.have.length(100);
  Object.isFrozen(got).should.be.true();
  Object.isFrozen(got.records[0].alleles).should.be.true();
  got.skipped.should.eql([
    { start: null, end: null, reason: 'not_an_object' },
    { start: 11109, end: 11109, reason: 'feature_type' },
    { start: 11200, end: 11200, reason: 'region_mismatch' },
    { start: null, end: null, reason: 'invalid_coordinates' },
    { start: 11300, end: 11300, reason: 'invalid_alleles' },
    { start: 11301, end: 11301, reason: 'invalid_alleles' },
    { start: 11302, end: 11302, reason: 'invalid_alleles' }
  ]);

  // the last chunk of a region stops at its end; a chunk start off the grid never fetches
  await w.client.overlapChunk(SPECIES, '1', 80880001, { regionLength: CHR1_LENGTH });
  w.calls[1].url.should.equal(BASE + '/overlap/region/sorghum_bicolor/1:80880001-80884392?feature=variation;content-type=application/json');
  (function () { w.client.overlapChunk(SPECIES, '1', 10000); }).should.throw(TypeError);
  (function () { w.client.overlapChunk('Sorghum bicolor', '1', 1); }).should.throw(TypeError);
  (function () { w.client.overlapChunk(SPECIES, '1', 1, { regionLength: 0 }); }).should.throw(TypeError);
  w.calls.should.have.length(2);
});

test('ids with "," and "*", and a 255-character id, are accepted and percent-encoded under the configured base; bad ids never fetch', async function () {
  const w = make({ handler: function () { return respond(fixture('variation_rs871475760.json')); } });
  const long = 'a'.repeat(255);
  await w.client.variation(SPECIES, 'tmp_1_13549_TTA_T,*');
  await w.client.variation(SPECIES, long);
  await w.client.variation(SPECIES, 'rs871475760');
  w.calls.map(function (c) { return c.url; }).should.eql([
    BASE + '/variation/sorghum_bicolor/tmp_1_13549_TTA_T%2C*?content-type=application/json',
    BASE + '/variation/sorghum_bicolor/' + long + '?content-type=application/json',
    BASE + '/variation/sorghum_bicolor/rs871475760?content-type=application/json'
  ]);
  w.calls.forEach(function (c) { c.url.indexOf(BASE + '/').should.equal(0); });
  ['a'.repeat(256), '../etc', '.hidden', 'a/b', 'a b', '', 'rs1?x=1', 'rs1#x', null].forEach(function (bad) {
    let thrown = null;
    try {
      w.client.variation(SPECIES, bad);
    } catch (e) {
      thrown = e;
    }
    should(thrown).be.instanceOf(TypeError, String(bad));
  });
  w.calls.should.have.length(3);

  // a trailing slash on the base changes nothing; a region name is encoded as one path segment
  const slash = make({ variation: { base_url: BASE + '/' } });
  await slash.client.overlapChunk(SPECIES, 'scaffold 1/x', 1);
  slash.calls[0].url.should.equal(BASE + '/overlap/region/sorghum_bicolor/scaffold%201%2Fx:1-10000?feature=variation;content-type=application/json');
});

test('record ids with ",", "*" or 200 characters are kept; an id outside even the lenient rule keeps its record with id null', async function () {
  const w = make({
    handler: function () {
      return respond([record({ id: 'tmp_1_13549_TTA_T,*' }), record({ id: 'x'.repeat(200) }), record({ id: 'badid' }),
        record({ id: 'y'.repeat(256) }), record({ id: 7 })]);
    }
  });
  const got = await w.client.overlapChunk(SPECIES, '1', 10001);
  got.records.map(function (r) { return r.id; }).should.eql(['tmp_1_13549_TTA_T,*', 'x'.repeat(200), null, null, null]);
  got.skipped.should.eql([]); // so no VARIATION_RECORDS_SKIPPED
});

test('variation lookup: compact mappings, "." synonyms dropped; a body without mappings is invalid_response', async function () {
  const w = make({
    handler: function (url) {
      if (/rs871475760/.test(url)) return respond(fixture('variation_rs871475760.json'));
      if (/rs5413864115/.test(url)) return respond(fixture('variation_rs5413864115.json'));
      return respond({ name: 'x', mappings: [{ seq_region_name: '1', start: 1 }] });
    }
  });
  (await w.client.variation(SPECIES, 'rs871475760')).should.eql({
    name: 'rs871475760',
    mappings: [{ seq_region_name: '1', start: 11109, end: 11109, allele_string: 'C/A' }],
    synonyms: ['tmp_1_11109_C_A'],
    consequence: 'downstream_gene_variant'
  });
  (await w.client.variation(SPECIES, 'rs5413864115')).mappings.should.eql([{ seq_region_name: '1', start: 11283, end: 11283, allele_string: 'A/-' }]);
  unavailableWith(await rejection(w.client.variation(SPECIES, 'rs1')), 'invalid_response', 30);
  should(_internal.compactLookup([])).be.null();
  should(_internal.compactLookup({ mappings: {} })).be.null();
  should(_internal.compactOverlap({})).be.null();
});

// ---- status mapping -------------------------------------------------------------------------------------------

test('400 "not found" on /variation is 404 UNKNOWN_VARIANT, cached for 5 minutes, and leaves the breaker alone', async function () {
  const w = make({ handler: function () { return respond(fixture('variation_not_found.json'), 400); } });
  const e = await rejection(w.client.variation(SPECIES, 'rs0000000001'));
  e.should.be.instanceOf(VariationError);
  e.should.match({ status: 404, code: 'UNKNOWN_VARIANT', details: { id: 'rs0000000001' } });
  w.client._state().breaker.failures.should.equal(0);
  for (let i = 0; i < 4; i++) (await rejection(w.client.variation(SPECIES, 'rs0000000001'))).code.should.equal('UNKNOWN_VARIANT');
  w.calls.should.have.length(1);
  w.client._state().breaker.should.eql({ failures: 0, open_until: 0 });
  w.clock.advance(300000 - 1);
  await rejection(w.client.variation(SPECIES, 'rs0000000001'));
  w.calls.should.have.length(1);
  w.clock.advance(2);
  await rejection(w.client.variation(SPECIES, 'rs0000000001'));
  w.calls.should.have.length(2);

  // the same body on the overlap endpoint is not an answer about a variant
  unavailableWith(await rejection(w.client.overlapChunk(SPECIES, '1', 1)), 'invalid_response', 30);
});

test('400 other text, HTML 404, 3xx, non-JSON 200, unexpected shape, redirects and bodies over the cap: invalid_response, not cached, breaker +1', async function () {
  const html = '<html><head><title>404 Not Found</title></head><body>' + 'y'.repeat(500) + '</body></html>';
  const cases = [
    ['variation', function () { return respond('{"error":"Cannot fetch slice for this id"}', 400); }],
    ['variation', function () { return respond('Bad request: species not in registry', 400); }],
    ['overlap', function () { return respond(html, 404, { 'content-type': 'text/html' }); }],
    ['overlap', function () { return respond('', 302, { location: 'https://elsewhere.example/' }); }],
    ['overlap', function () { return respond('{"not json', 200); }],
    ['overlap', function () { return respond({ error: 'nope' }, 200); }],
    ['variation', function () { return respond([], 200); }],
    ['overlap', function () { throw Object.assign(new TypeError('fetch failed'), { cause: new Error('unexpected redirect') }); }],
    ['overlap', function () { return respond('[' + '"x",'.repeat(700) + '"x"]', 200); }], // 2.8 KB streamed, no content-length
    ['overlap', function () { return respond('[]', 200, { 'content-length': '999999' }); }],
    ['variation', function () { return respond(JSON.stringify({ mappings: [], pad: 'p'.repeat(1200) }), 200); }]
  ];
  let failures = 0;
  for (const [kind, handler] of cases) {
    const w = make({ handler: handler, variation: { breaker_failures: 100, max_chunk_bytes: 2000, max_lookup_bytes: 1000 } });
    const call = function () { return kind === 'overlap' ? w.client.overlapChunk(SPECIES, '1', 1) : w.client.variation(SPECIES, 'rs1'); };
    unavailableWith(await rejection(call()), 'invalid_response', 30);
    w.client._state().breaker.failures.should.equal(1);
    w.client._state().negative.should.equal(0);
    unavailableWith(await rejection(call()), 'invalid_response', 30);
    w.calls.should.have.length(2); // not cached
    w.client._state().breaker.failures.should.equal(2);
    const logged = w.log.lines.join('\n');
    logged.should.match(/data\.gramene\.org/);
    logged.should.not.match(/pansite|overlap\/region|variation\/sorghum/);
    w.log.lines.forEach(function (l) { l.length.should.be.below(400); });
    failures++;
  }
  failures.should.equal(cases.length);
});

test('503, 429, timeout and ECONNREFUSED: 503 with the matching reason and retry_after_s 30, negatively cached for 30 s', async function () {
  const cases = [
    ['http_5xx', function () { return respond({ error: 'busy' }, 503); }],
    ['rate_limited', function () { return respond('slow down', 429); }],
    ['timeout', function (url, init) {
      return new Promise(function (resolve, reject) {
        init.signal.addEventListener('abort', function () { reject(init.signal.reason); }, { once: true });
      });
    }],
    ['transport', function () { throw refused(); }]
  ];
  for (const [reason, handler] of cases) {
    const w = make({ handler: handler, variation: { timeout_ms: 30, breaker_failures: 100 } });
    unavailableWith(await rejection(w.client.overlapChunk(SPECIES, '1', 10001)), reason, 30);
    w.client._state().negative.should.equal(1);
    w.client._state().breaker.failures.should.equal(1);
    w.clock.advance(20000);
    unavailableWith(await rejection(w.client.overlapChunk(SPECIES, '1', 10001)), reason, 10); // the remaining seconds
    w.calls.should.have.length(1);
    w.clock.advance(10001);
    unavailableWith(await rejection(w.client.overlapChunk(SPECIES, '1', 10001)), reason, 30);
    w.calls.should.have.length(2);
  }
});

test('breaker: 3 failures within 60 s open it; the next call fails with breaker_open and no fetch; cached answers still serve; it closes after 60 s', async function () {
  let down = false;
  const w = make({ handler: function () { if (down) throw refused(); return respond([record()]); } });
  await w.client.overlapChunk(SPECIES, '1', 40001);
  down = true;
  for (const start of [1, 10001, 20001]) unavailableWith(await rejection(w.client.overlapChunk(SPECIES, '1', start)), 'transport', 30);
  w.calls.should.have.length(4);
  w.client._state().breaker.open_until.should.equal(w.clock.t + 60000);

  unavailableWith(await rejection(w.client.overlapChunk(SPECIES, '1', 30001)), 'breaker_open', 60);
  w.calls.should.have.length(4);
  (await w.client.overlapChunk(SPECIES, '1', 40001)).records.should.have.length(1); // a positive cache hit
  w.clock.advance(45000);
  unavailableWith(await rejection(w.client.variation(SPECIES, 'rs1')), 'breaker_open', 15);
  w.calls.should.have.length(4);

  w.clock.advance(15001);
  down = false;
  (await w.client.overlapChunk(SPECIES, '1', 30001)).records.should.have.length(1);
  w.calls.should.have.length(5);
  w.client._state().breaker.should.eql({ failures: 0, open_until: 0 });

  // failures further apart than the window never open it
  const spaced = make({ handler: function () { throw refused(); } });
  for (const start of [1, 10001, 20001, 30001]) {
    unavailableWith(await rejection(spaced.client.overlapChunk(SPECIES, '1', start)), 'transport');
    spaced.clock.advance(30000);
  }
  spaced.calls.should.have.length(4);
});

// ---- limiter, single flight and caches -------------------------------------------------------------------------

test('limiter: with max_concurrent 4 a fifth call waits, then fails with queue_full and retry_after_s 5, without fetching', async function () {
  const gates = [];
  const w = make({
    variation: { queue_wait_ms: 40 },
    handler: function () {
      const d = deferred();
      gates.push(d);
      return d.promise.then(function () { return respond([record()]); });
    }
  });
  const running = [1, 10001, 20001, 30001].map(function (s) { return w.client.overlapChunk(SPECIES, '1', s); });
  const fifth = w.client.overlapChunk(SPECIES, '1', 40001);
  await new Promise(function (resolve) { setImmediate(resolve); });
  w.calls.should.have.length(4);
  w.client._state().should.match({ running: 4, waiting: 1 });
  const e = await rejection(fifth);
  unavailableWith(e, 'queue_full', 5);
  w.calls.should.have.length(4);
  w.client._state().breaker.failures.should.equal(0);
  gates.forEach(function (g) { g.resolve(); });
  (await Promise.all(running)).should.have.length(4);
  w.client._state().running.should.equal(0);
});

test('single flight: concurrent callers of one key share one fetch; an aborted caller stops waiting and the shared request completes', async function () {
  const gate = deferred();
  let fetchSignal = null;
  const w = make({ handler: function (url, init) { fetchSignal = init.signal; return gate.promise.then(function () { return respond([record()]); }); } });
  const ac = new AbortController();
  const first = w.client.overlapChunk(SPECIES, '1', 10001, { signal: ac.signal });
  const second = w.client.overlapChunk(SPECIES, '1', 10001);
  const third = w.client.overlapChunk(SPECIES, '1', 10001, { signal: new AbortController().signal });
  await new Promise(function (resolve) { setImmediate(resolve); });
  w.calls.should.have.length(1);
  w.client._state().inflight.should.equal(1);

  const reason = new PrimerHttpError(400, 'CLIENT_CLOSED_REQUEST', 'the client closed the connection', {});
  ac.abort(reason);
  (await rejection(first)).should.equal(reason);
  gate.resolve();
  (await second).records.should.have.length(1);
  (await third).records.should.have.length(1);
  fetchSignal.aborted.should.be.false();
  w.client._state().should.match({ inflight: 0, cache: 1 });
  (await w.client.overlapChunk(SPECIES, '1', 10001)).records.should.have.length(1);
  w.calls.should.have.length(1);

  // a caller that gave up before calling starts no request; a plain abort becomes 504 DEADLINE_EXCEEDED
  const gone = new AbortController();
  gone.abort();
  (await rejection(w.client.variation(SPECIES, 'rs1', { signal: gone.signal }))).should.match({ status: 504, code: 'DEADLINE_EXCEEDED' });
  w.calls.should.have.length(1);
});

test('caches: an LRU of cache_entries validated values, each positive entry kept for an hour', async function () {
  const w = make({ variation: { cache_entries: 2 }, handler: function () { return respond([record()]); } });
  for (const s of [1, 10001, 20001]) await w.client.overlapChunk(SPECIES, '1', s);
  w.calls.should.have.length(3);
  w.client._state().cache.should.equal(2);
  await w.client.overlapChunk(SPECIES, '1', 20001);
  w.calls.should.have.length(3);
  await w.client.overlapChunk(SPECIES, '1', 1); // evicted as least recently used
  w.calls.should.have.length(4);
  w.clock.advance(3600000 - 1);
  await w.client.overlapChunk(SPECIES, '1', 1);
  w.calls.should.have.length(4);
  w.clock.advance(2);
  await w.client.overlapChunk(SPECIES, '1', 1);
  w.calls.should.have.length(5);
});

test('a 25 kb window fetches exactly 3 chunks; a second overlapping window reuses them', async function () {
  const w = make({ handler: function () { return respond([]); } });
  const cfg = cfgWith({});
  const deps = {
    cfg: cfg,
    client: w.client,
    resolve: async function (name) { return { system_name: name, fasta: { dna: '/fake/' + name + '.fa.gz' } }; },
    sequence: {
      regionLength: async function () { return CHR1_LENGTH; },
      fetch: async function (fasta, region, s, e) { return 'A'.repeat(e - s + 1); }
    }
  };
  const body = await variation.listVariants({ system_name: SPECIES, region: '1', start: 2, end: 25001 }, deps);
  body.should.match({ total: 0, returned: 0, truncated: false, variants: [], warnings: [] });
  w.calls.map(function (c) { return /:(\d+-\d+)\?/.exec(c.url)[1]; }).sort().should.eql(['1-10000', '10001-20000', '20001-30000']);
  await variation.listVariants({ system_name: SPECIES, region: '1', start: 15000, end: 29999 }, deps);
  w.calls.should.have.length(3);
});

test('VariationError is a PrimerHttpError: sendError answers 503 with Retry-After and a client-safe body', function () {
  const e = new VariationError(503, 'VARIATION_SOURCE_UNAVAILABLE', 'the Ensembl variation service is temporarily unavailable', { retry_after_s: 30, reason: 'timeout' });
  errors.isPrimerHttpError(e).should.be.true();
  e.reason.should.equal('timeout');
  const res = {
    headers: {}, statusCode: 200, headersSent: false, body: null,
    set: function (k, v) { this.headers[k.toLowerCase()] = v; return this; },
    status: function (c) { this.statusCode = c; return this; },
    json: function (b) { this.body = b; this.headersSent = true; return this; }
  };
  errors.sendError(res, e, { log: quietLog() });
  res.statusCode.should.equal(503);
  res.headers['retry-after'].should.equal('30');
  res.body.should.eql({ message: 'the Ensembl variation service is temporarily unavailable', code: 'VARIATION_SOURCE_UNAVAILABLE', details: { retry_after_s: 30, reason: 'timeout' } });
});
