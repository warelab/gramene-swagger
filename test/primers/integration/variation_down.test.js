'use strict';

// Integration: Ensembl unavailable (genotyping spec §3.3, §3.8) and the design semaphore under a slow Ensembl (§4.2, the
// slot-starvation regression test of §7.6), against the loopback fake Ensembl (helpers/fake_ensembl.js) and an API that this
// file starts and stops itself (helpers/dev_api.js) on the port of PRIMERS_IT_BASE, with PRIMERS_VARIATION_URL pointing at
// the fake. Skipped unless PRIMERS_IT_BASE is set:
//
//   PRIMERS_IT_BASE=http://127.0.0.1:50112/sorghum_v11 node --test --test-concurrency=1 test/primers/integration/variation_down.test.js
//
// The port must be free: stop the dev API first, or set PRIMERS_IT_DOWN_PORT to another free loopback port. The client's
// circuit breaker (3 failures within 60 s open it for 60 s) and its negative cache (upstream failures, 30 s) live in the API
// process, so each group restarts the API and starts clean with the production defaults. No worker is started and no check
// job is submitted. It takes about two minutes, because the slow groups wait on the fake by design.

const { describe, it, before, after } = require('node:test');
const should = require('should');
const fs = require('fs');
const path = require('path');

const { createFakeEnsembl } = require('./helpers/fake_ensembl');
const { startApi, portOpen } = require('./helpers/dev_api');

const IT_BASE = (process.env.PRIMERS_IT_BASE || '').replace(/\/+$/, '');
const SKIP = IT_BASE ? false : 'set PRIMERS_IT_BASE, e.g. http://127.0.0.1:50112/sorghum_v11 (this file starts its own API on that port)';
const TARGET = IT_BASE ? new URL(IT_BASE) : null;
const PORT = Number(process.env.PRIMERS_IT_DOWN_PORT) || (TARGET ? Number(TARGET.port) : 0);
const BASE_PATH = TARGET ? TARGET.pathname : '/sorghum_v11';
const CAPTURES = path.join(__dirname, '..', 'fixtures', 'docs');
const GROUP_TIMEOUT_MS = 180000;

const CLIENT_TIMEOUT_MS = 8000; // primers.variation.timeout_ms
const SEMAPHORE_WAIT_MS = 10000; // primers.design.wait_timeout_ms: a queued ordinary design gives up with 503 BUSY after this
const CHUNK = 'o|sorghum_bicolor|1|10001'; // the 10 kb overlap chunk holding every example variant
const LOOKUP = function (id) { return 'v|sorghum_bicolor|' + id; };
const LIST = '/primers/variants?system_name=sorghum_bicolor&region=1&start=11180&end=11290';

// Ordinary designs: they never call Ensembl and must never wait behind it.
const ORDINARY = [
  { mode: 'transcript', gene_id: 'SORBI_3001G000200' },
  { mode: 'gene', gene_id: 'SORBI_3001G000700' },
  { mode: 'region', system_name: 'sorghum_bicolor', region: { region: '1', start: 11080, end: 12079 } }
];
const WAVE = 16;
// Genotyping designs by id, one Ensembl lookup key each; their overlap chunk is shared (single flight).
const BY_ID = [
  { system_name: 'sorghum_bicolor', variant: { id: 'rs871475760', alt: 'A' }, assay: { type: 'kasp', num_sets: 2 } },
  { system_name: 'sorghum_bicolor', variant: { id: 'tmp_1_11502_C_CGT' }, assay: { num_sets: 1 } },
  { system_name: 'sorghum_bicolor', variant: { id: 'rs5413864115' }, assay: { num_sets: 1 } },
  { system_name: 'sorghum_bicolor', variant: { id: 'tmp_1_11193_C_T' }, assay: { num_sets: 1 } }
];
const MANUAL_SNV = { system_name: 'sorghum_bicolor', variant: { region: '1', position: 11109, ref: 'C', alt: 'A' }, assay: { type: 'kasp', num_sets: 2 } };
const MANUAL_INSERTION = { system_name: 'sorghum_bicolor', variant: { region: '1', position: 11502, ref: 'C', alt: 'CGT' }, assay: { num_sets: 1 } };

const sleep = (ms) => new Promise((resolve) => setTimeout(resolve, ms));

function capture(name) {
  return JSON.parse(fs.readFileSync(path.join(CAPTURES, 'capture-' + name + '.json'), 'utf8'));
}

function summarize(results, startedAtMs) {
  const ms = results.map(function (r) { return r.ms; }).sort(function (a, b) { return a - b; });
  const statuses = {};
  results.forEach(function (r) {
    const k = r.status + (r.json && r.json.code ? ' ' + r.json.code + (r.json.details && r.json.details.reason ? '/' + r.json.details.reason : '') : '');
    statuses[k] = (statuses[k] || 0) + 1;
  });
  return { n: results.length, fired_at_ms: startedAtMs, statuses: statuses, ms_min: ms[0], ms_median: ms[Math.floor(ms.length / 2)], ms_max: ms[ms.length - 1] };
}

describe('Ensembl unavailable and slow: fake Ensembl and an API of its own (spec §3.8, §4.2)', { skip: SKIP }, function () {
  let fake = null;
  let api = null;
  const timings = {};

  async function http(method, p, opts) {
    opts = opts || {};
    const headers = {};
    let body;
    if (opts.json !== undefined) {
      body = JSON.stringify(opts.json);
      headers['Content-Type'] = 'application/json';
    }
    const t0 = Date.now();
    const res = await fetch(api.base + p, { method: method, headers: headers, body: body, redirect: 'manual' });
    const text = await res.text();
    let json = null;
    try { json = JSON.parse(text); } catch (e) { /* not JSON */ }
    return { status: res.status, headers: res.headers, json: json, text: text, ms: Date.now() - t0 };
  }

  const get = (p) => http('GET', p);
  const genotyping = (body) => http('POST', '/primers/genotyping/design', { json: body });
  const ordinary = (i) => http('POST', '/primers/design', { json: ORDINARY[i % ORDINARY.length] });

  function expectStatus(r, status) {
    r.status.should.equal(status, JSON.stringify(r.json || r.text).slice(0, 500));
    r.headers.get('cache-control').should.equal('no-store');
    return r;
  }

  // 503 VARIATION_SOURCE_UNAVAILABLE with one of `reasons`, retry_after_s and the matching Retry-After header.
  function expectUnavailable(r, reasons) {
    expectStatus(r, 503);
    r.json.code.should.equal('VARIATION_SOURCE_UNAVAILABLE', JSON.stringify(r.json));
    reasons.should.containEql(r.json.details.reason);
    r.json.details.retry_after_s.should.be.a.Number().and.be.aboveOrEqual(1);
    r.headers.get('retry-after').should.equal(String(r.json.details.retry_after_s));
    return r.json.details;
  }

  async function restart(t, mode, options) {
    if (api) await api.stop();
    await fake.setMode(mode, options);
    fake.reset();
    api = await startApi({ port: PORT, basePath: BASE_PATH, env: { PRIMERS_VARIATION_URL: fake.url } });
    should(api.siteKey).equal('sorghum_v11_geno');
    api.log().should.not.match(/PRIMERS_VARIATION_URL must be/);
    t.diagnostic('API pid ' + api.pid + ' on ' + api.base + ' up in ' + api.startup_ms + ' ms; fake Ensembl ' + fake.url + ' mode ' + mode + ' ' + JSON.stringify(options || {}));
  }

  // The slot-starvation scenario: genotyping designs by id wait on the slow fake while waves of ordinary designs run.
  // waves: [{when: 'lookups' | 'overlap' | ms after the genotyping designs were fired}]
  async function starvation(t, label, waves) {
    const t0 = Date.now();
    const geno = BY_ID.map(genotyping);
    await fake.waitFor(function () { return fake.inflight({ kind: 'variation' }) === BY_ID.length; }, 5000, BY_ID.length + ' lookups in flight');
    const waveResults = [];
    for (const w of waves) {
      if (w.when === 'overlap') await fake.waitFor(function () { return fake.inflight({ kind: 'overlap' }) >= 1; }, 3 * CLIENT_TIMEOUT_MS, 'the overlap chunk in flight');
      else if (typeof w.when === 'number') await sleep(Math.max(0, t0 + w.when - Date.now()));
      const firedAt = Date.now() - t0;
      const inflight = { lookups: fake.inflight({ kind: 'variation' }), overlaps: fake.inflight({ kind: 'overlap' }) };
      const results = await Promise.all(Array.from({ length: WAVE }, function (_, i) { return ordinary(i); }));
      const s = Object.assign(summarize(results, firedAt), { ensembl_inflight_when_fired: inflight });
      waveResults.push({ results: results, summary: s });
    }
    const genoResults = await Promise.all(geno);
    const out = {
      waves: waveResults.map(function (w) { return w.summary; }),
      genotyping: genoResults.map(function (r, i) {
        return { id: BY_ID[i].variant.id, status: r.status, code: r.json && r.json.code, reason: r.json && r.json.details && r.json.details.reason, ms: r.ms };
      }),
      ensembl_requests: { lookups: fake.count({ kind: 'variation' }), overlaps: fake.count({ kind: 'overlap' }) }
    };
    timings[label] = out;
    t.diagnostic(label + ' ' + JSON.stringify(out));

    // No ordinary design may be refused, and none may come near the semaphore wait timeout.
    waveResults.forEach(function (w, i) {
      const refused = w.results.filter(function (r) { return r.status !== 200; });
      refused.map(function (r) { return r.status + ' ' + (r.json && r.json.code); }).should.eql([], 'wave ' + (i + 1) + ': ordinary designs refused');
      w.results.filter(function (r) { return r.json && r.json.code === 'BUSY'; }).length.should.equal(0);
      w.results.forEach(function (r) {
        r.headers.get('cache-control').should.equal('no-store');
        r.json.pairs.length.should.be.above(0);
        r.ms.should.be.below(SEMAPHORE_WAIT_MS);
      });
    });
    return { geno: genoResults, out: out };
  }

  before(async function () {
    fake = await createFakeEnsembl().start();
  });

  after(async function () {
    if (api) await api.stop();
    if (fake) await fake.stop();
  });

  it('normal mode: the API uses the fake (PRIMERS_VARIATION_URL accepted) and the recorded bodies reproduce the live captures', { timeout: GROUP_TIMEOUT_MS }, async function (t) {
    await restart(t, 'normal');
    for (const name of ['variants-list-1_11180-11290', 'variants-lookup-tmp_1_11502_C_CGT', 'variants-lookup-unknown-variant', 'genotyping-design-rs871475760-kasp']) {
      const c = capture(name);
      const r = await http(c.request.method, c.request.path.slice(BASE_PATH.length), c.request.body !== undefined ? { json: c.request.body } : {});
      expectStatus(r, c.status);
      r.json.should.eql(c.response, name);
    }
    // One chunk fetch serves the listing, the lookup and the design (cache); one lookup per id.
    fake.count({ kind: 'overlap' }).should.equal(1);
    fake.count({ key: CHUNK }).should.equal(1);
    fake.count({ kind: 'variation' }).should.equal(3);
    fake.count({ kind: 'other' }).should.equal(0);
  });

  it('slot starvation, Ensembl slow 7 s (under the 8 s client timeout): genotyping designs by id wait outside the design semaphore; 2 waves of 16 ordinary designs get no 503 BUSY', { timeout: GROUP_TIMEOUT_MS }, async function (t) {
    const delay = 7000;
    await restart(t, 'slow', { delay_ms: delay });
    const { geno, out } = await starvation(t, 'slow_7s', [{ when: 'lookups' }, { when: 'overlap' }]);
    // Each design by id waited on its lookup and then on the shared overlap chunk, and then designed normally.
    geno.forEach(function (r, i) {
      expectStatus(r, 200);
      r.ms.should.be.aboveOrEqual(2 * delay - 1000, BY_ID[i].variant.id);
      r.json.sets.length.should.be.above(0);
    });
    geno[0].json.sets.map(function (s) { return s.key; }).should.eql(['f9df650ad116', '1accc54c262d']);
    out.ensembl_requests.should.eql({ lookups: BY_ID.length, overlaps: 1 });
    out.waves.forEach(function (w) { (w.ensembl_inflight_when_fired.lookups + w.ensembl_inflight_when_fired.overlaps).should.be.above(0); });
  });

  it('slot starvation, Ensembl slow 9 s (over the 8 s client timeout): genotyping designs by id fail with 503 VARIATION_SOURCE_UNAVAILABLE; 2 waves of 16 ordinary designs get no 503 BUSY', { timeout: GROUP_TIMEOUT_MS }, async function (t) {
    await restart(t, 'slow', { delay_ms: 9000 });
    const { geno, out } = await starvation(t, 'slow_9s', [{ when: 'lookups' }, { when: 5000 }]);
    geno.forEach(function (r, i) {
      expectUnavailable(r, ['timeout', 'breaker_open']);
      r.ms.should.be.within(CLIENT_TIMEOUT_MS - 500, CLIENT_TIMEOUT_MS + 5000, BY_ID[i].variant.id);
    });
    out.ensembl_requests.should.eql({ lookups: BY_ID.length, overlaps: 0 });
  });

  it('5xx: the variants endpoints answer 503 http_5xx; repeats come from the negative cache without calling Ensembl; the breaker opens only after 3 failures on different keys', { timeout: GROUP_TIMEOUT_MS }, async function (t) {
    await restart(t, '5xx');
    // failure 1: the overlap chunk
    let d = expectUnavailable(await get(LIST), ['http_5xx']);
    d.retry_after_s.should.equal(30);
    fake.count({ key: CHUNK }).should.equal(1);
    // the same key again, and another window in the same chunk: negative cache, no request, no breaker failure
    for (const p of [LIST, LIST, '/primers/variants?system_name=sorghum_bicolor&region=1&start=11300&end=11400']) {
      const r = await get(p);
      expectUnavailable(r, ['http_5xx']).retry_after_s.should.be.within(1, 30);
      r.ms.should.be.below(1000);
    }
    fake.count({ key: CHUNK }).should.equal(1);
    // failure 2: a lookup; its repeat is cached too
    expectUnavailable(await get('/primers/variants/rs871475760?system_name=sorghum_bicolor'), ['http_5xx']);
    expectUnavailable(await get('/primers/variants/rs871475760?system_name=sorghum_bicolor'), ['http_5xx']);
    fake.count({ key: LOOKUP('rs871475760') }).should.equal(1);
    // six 503s so far but two upstream failures: the breaker is still closed, so a third key reaches Ensembl (failure 3)
    expectUnavailable(await get('/primers/variants/tmp_1_11502_C_CGT?system_name=sorghum_bicolor'), ['http_5xx']);
    fake.count({ key: LOOKUP('tmp_1_11502_C_CGT') }).should.equal(1);
    // now open: new keys fail at once and never reach Ensembl
    d = expectUnavailable(await get('/primers/variants/rs5413864115?system_name=sorghum_bicolor'), ['breaker_open']);
    d.retry_after_s.should.be.within(1, 60);
    expectUnavailable(await get('/primers/variants?system_name=sorghum_bicolor&region=1&start=30001&end=30100'), ['breaker_open']);
    fake.count({ key: LOOKUP('rs5413864115') }).should.equal(0);
    fake.count({ key: 'o|sorghum_bicolor|1|30001' }).should.equal(0);
    fake.requests.length.should.equal(3);
  });

  it('5xx with the breaker open: a manual design still returns 200 with NEIGHBOURS_UNAVAILABLE; designs by id return 503; ordinary designs and the genomes flags are unaffected', { timeout: GROUP_TIMEOUT_MS }, async function (t) {
    const m = expectStatus(await genotyping(MANUAL_SNV), 200).json;
    m.variant.should.match({ key: '1:11109:C:A', requested_id: null, ids: [] });
    m.neighbours.should.match({ data: 'unavailable', variants: 0 });
    const w = m.warnings.filter(function (x) { return x.code === 'NEIGHBOURS_UNAVAILABLE'; });
    w.should.have.length(1);
    ['http_5xx', 'breaker_open'].should.containEql(w[0].details.reason);
    m.sets.length.should.be.above(0);
    t.diagnostic('manual 1:11109 C/A without neighbours: set keys ' + JSON.stringify(m.sets.map(function (s) { return s.key; })) + ', warnings ' + JSON.stringify(m.warnings.map(function (x) { return x.code; })));

    expectUnavailable(await genotyping(BY_ID[0]), ['http_5xx']);
    expectUnavailable(await genotyping(BY_ID[3]), ['breaker_open']);
    const o = expectStatus(await ordinary(0), 200).json;
    o.pairs.length.should.be.above(0);
    const g = expectStatus(await get('/primers/genomes?system_name=sorghum_bicolor'), 200).json;
    g.variation.should.eql({ available: true, source: 'ensembl', release: '115' });
    fake.requests.length.should.equal(3);
  });

  it('down (connection refused): the list answers 503 transport within 10 s and its repeat at once; a lookup 503; a manual design 200 with NEIGHBOURS_UNAVAILABLE; a design by id 503', { timeout: GROUP_TIMEOUT_MS }, async function (t) {
    await restart(t, 'down');
    fake.listening.should.be.false();
    (await portOpen(fake.port)).should.be.false();
    let r = await get(LIST);
    expectUnavailable(r, ['transport']);
    r.ms.should.be.below(10000);
    const first = r.ms;
    r = await get(LIST);
    expectUnavailable(r, ['transport']);
    r.ms.should.be.below(1000);
    t.diagnostic('list 503 in ' + first + ' ms, repeat in ' + r.ms + ' ms');
    expectUnavailable(await get('/primers/variants/tmp_1_11502_C_CGT?system_name=sorghum_bicolor'), ['transport']);

    const m = expectStatus(await genotyping(MANUAL_INSERTION), 200).json;
    m.variant.should.match({ key: '1:11502:C:CGT', requested_id: null, ids: [] });
    m.neighbours.data.should.equal('unavailable');
    const w = m.warnings.find(function (x) { return x.code === 'NEIGHBOURS_UNAVAILABLE'; });
    should.exist(w, JSON.stringify(m.warnings));
    w.details.should.eql({ reason: 'transport' });
    // no neighbour data, so nothing blocks the forward orientation (rs5413863234 blocks it when Ensembl answers)
    m.orientations.forward.status.should.not.equal('blocked');

    expectUnavailable(await genotyping(BY_ID[0]), ['transport', 'breaker_open']);
    fake.requests.length.should.equal(0);
  });

  it('malformed and HTML 404 bodies: 503 invalid_response, never 404 UNKNOWN_VARIANT; not negatively cached; they count toward the breaker even on one key', { timeout: GROUP_TIMEOUT_MS }, async function (t) {
    await restart(t, 'malformed');
    expectUnavailable(await get(LIST), ['invalid_response']);
    expectUnavailable(await get(LIST), ['invalid_response']);
    // not cached: the same chunk was requested twice
    fake.count({ key: CHUNK }).should.equal(2);

    await fake.setMode('html404');
    const d = expectUnavailable(await get('/primers/variants/rs0000000001?system_name=sorghum_bicolor'), ['invalid_response']);
    d.retry_after_s.should.equal(30);
    fake.count({ key: LOOKUP('rs0000000001') }).should.equal(1);
    // three failures, two on one key: open
    expectUnavailable(await get('/primers/variants/rs871475760?system_name=sorghum_bicolor'), ['breaker_open']);
    fake.count({ key: LOOKUP('rs871475760') }).should.equal(0);
    expectUnavailable(await genotyping(BY_ID[1]), ['breaker_open']);
    fake.requests.length.should.equal(3);
  });

  it('oversize bodies: a padded answer under the 1 MB overlap cap is accepted; over it (streamed, no Content-Length) and over the 256 KB lookup cap (declared) give 503 invalid_response', { timeout: GROUP_TIMEOUT_MS }, async function (t) {
    await restart(t, 'oversize', { bytes: 300000 });
    // the padding is valid JSON: only the byte caps can refuse it
    const ok = expectStatus(await get(LIST), 200).json;
    ok.variants.map(function (v) { return v.key; }).should.eql(['1:11182:A:G', '1:11193:C:T', '1:11203:C:T', '1:11282:CA:C']);

    await fake.setMode('oversize', { bytes: 1500000 });
    let r = await get('/primers/variants?system_name=sorghum_bicolor&region=1&start=30001&end=30100');
    expectUnavailable(r, ['invalid_response']);
    t.diagnostic('1.5 MB streamed chunk refused in ' + r.ms + ' ms');

    await fake.setMode('oversize', { bytes: 300000, declare_length: true });
    r = await get('/primers/variants/rs871475760?system_name=sorghum_bicolor');
    expectUnavailable(r, ['invalid_response']);
    fake.count({ key: LOOKUP('rs871475760') }).should.equal(1);
  });

  it('teardown: the API stops and its port closes; the timings', { timeout: GROUP_TIMEOUT_MS }, async function (t) {
    await api.stop();
    (await portOpen(PORT)).should.be.false();
    api = null;
    await fake.stop();
    fake = null;
    t.diagnostic('timings ' + JSON.stringify(timings));
  });
});
