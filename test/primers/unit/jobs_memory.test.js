'use strict';

// jobs/index.js + memory_store.js: submit/status semantics, caps, stale pruning, requeue, TTL sweep.

const test = require('node:test');
const should = require('should');
const crypto = require('crypto');

const jobs = require('../../../api/helpers/primers/jobs');
const { createMemoryStore } = require('../../../api/helpers/primers/jobs/memory_store');
const workerLib = require('../../../api/helpers/primers/jobs/worker');
const config = require('../../../api/helpers/primers/config');

const quiet = { info() {}, warn() {}, error() {}, log() {} };

function makeCfg(check, top) {
  return config._build({
    env: {},
    fileConfig: {},
    overrides: Object.assign({ check: Object.assign({ global_prefix: 'primers:test_memory:' }, check || {}) }, top || {})
  }).config;
}

function newShared() {
  return { slots: new Map(), panSlots: new Map() };
}

function clock(start) {
  let t = start || 1757700000000;
  const f = function () { return t; };
  f.advance = function (ms) { t += ms; };
  return f;
}

function hexId(n) {
  return crypto.createHash('md5').update('job-' + n).digest('hex');
}

function fakeNorm(o) {
  o = o || {};
  const kind = o.kind || 'specificity';
  return {
    request: {
      system_name: 'sorghum_bicolor', mode: 'region',
      checks: kind === 'pangenome' ? ['pangenome', 'specificity'] : ['specificity'],
      genomes: kind === 'pangenome' ? ['sorghum_353', 'sorghum_rio'] : [],
      params: { max_product_size: 4000 },
      pairs: [{ id: o.pairId || 'P1', left: 'GGACAGCTCCACAACATATCAG', right: 'GGACATTTGAAGCCCATGGCC' }]
    },
    resolved: { assemblies: { sorghum_bicolor: { system_name: 'sorghum_bicolor', fingerprint: 'f'.repeat(40) } }, gene: null, species: null },
    kind: kind,
    warnings: o.warnings || [],
    estimate: { cpu_s: 8, total: kind === 'pangenome' ? 3 : 1 },
    dbs: { sorghum_bicolor: 'f'.repeat(40) }
  };
}

function doc(n, kind, now) {
  return jobs.newJobDoc(hexId(n), fakeNorm({ kind: kind, pairId: 'P' + n }), now || Date.now());
}

function stubCheck() {
  return {
    ALGORITHM_VERSION: '1',
    normalize: async function (body) {
      return fakeNorm({ pairId: body.pairs[0].id, kind: body.kind, warnings: body.warn ? [{ code: 'EXPECTED_IGNORED', message: 'x' }] : [] });
    },
    run: async function () { return {}; }
  };
}

async function rejects(promise, status, code) {
  let err;
  try {
    await promise;
  } catch (e) {
    err = e;
  }
  should.exist(err, 'expected ' + code);
  if (err.name !== 'PrimerHttpError') throw err;
  err.status.should.equal(status);
  err.code.should.equal(code);
  return err;
}

async function start(store, id, worker, now) {
  const cur = await store.getJob(id);
  const started = workerLib.startedDoc(cur, { worker: worker || 'w1', now: now || Date.now() });
  (await store.setJob(started, { expect: ['queued'], worker: '' })).should.be.true();
  return started;
}

test('idempotent submit: 202-style created true, then false with the same job id', async function () {
  const cfg = makeCfg();
  const store = createMemoryStore({ cfg: cfg, siteKey: 'site_a', shared: newShared() });
  const d = { cfg: cfg, store: store, check: stubCheck(), log: quiet };
  const first = await jobs.submit({ pairs: [{ id: 'A' }] }, d);
  first.created.should.be.true();
  first.job_id.should.match(/^[0-9a-f]{32}$/);
  first.status.should.equal('queued');
  first.kind.should.equal('specificity');
  first.queue_position.should.equal(0);
  first.progress.should.eql({ done: 0, total: 1, stage: 'queued', running: [] });
  first.estimate.should.eql({ cpu_s: 8 });
  first.created_at.should.match(/^\d{4}-\d\d-\d\dT/);
  Object.keys(first).sort().should.eql(['created', 'created_at', 'estimate', 'job_id', 'kind', 'progress', 'queue_position', 'status', 'warnings']);

  const again = await jobs.submit({ pairs: [{ id: 'A' }], warn: true }, d);
  again.created.should.be.false();
  again.job_id.should.equal(first.job_id);
  again.created_at.should.equal(first.created_at);
  again.queue_position.should.equal(0);
  again.warnings.should.eql([{ code: 'EXPECTED_IGNORED', message: 'x' }]); // this submission's warnings
  store._state().queues.spec.should.eql([first.job_id]);

  const other = await jobs.submit({ pairs: [{ id: 'B' }], kind: 'pangenome' }, d);
  other.created.should.be.true();
  other.kind.should.equal('pangenome');
  other.queue_position.should.equal(0); // first in the pan queue
  const third = await jobs.submit({ pairs: [{ id: 'C' }] }, d);
  third.queue_position.should.equal(1);
});

test('queue full: > max_queued (both queues) -> 503 QUEUE_FULL retry_after_s 60; an existing job still answers', async function () {
  const cfg = makeCfg({ max_queued: 2 });
  const store = createMemoryStore({ cfg: cfg, siteKey: 'site_a', shared: newShared() });
  const d = { cfg: cfg, store: store, check: stubCheck(), log: quiet };
  await jobs.submit({ pairs: [{ id: 'A' }] }, d);
  await jobs.submit({ pairs: [{ id: 'B' }], kind: 'pangenome' }, d);
  const err = await rejects(jobs.submit({ pairs: [{ id: 'C' }] }, d), 503, 'QUEUE_FULL');
  err.details.retry_after_s.should.equal(60);
  (await jobs.submit({ pairs: [{ id: 'A' }] }, d)).created.should.be.false();
});

test('re-POST of an errored job re-queues it (created true)', async function () {
  const cfg = makeCfg();
  const store = createMemoryStore({ cfg: cfg, siteKey: 'site_a', shared: newShared() });
  const d = { cfg: cfg, store: store, check: stubCheck(), log: quiet };
  const first = await jobs.submit({ pairs: [{ id: 'A' }] }, d);
  (await store.claim(Date.now())).id.should.equal(first.job_id);
  const started = await start(store, first.job_id);
  (await store.fail(workerLib.errorDoc(started, { code: 'CHECK_FAILED', message: 'boom' }, { now: Date.now() }), { worker: 'w1' })).should.be.true();
  const failed = await jobs.status(first.job_id, d);
  failed.status.should.equal('error');
  failed.error.should.eql({ code: 'CHECK_FAILED', message: 'boom' });
  const again = await jobs.submit({ pairs: [{ id: 'A' }] }, d);
  again.created.should.be.true();
  again.status.should.equal('queued');
  const st = await jobs.status(first.job_id, d);
  st.status.should.equal('queued');
  should(st.error).be.null();
  st.attempts.should.equal(0);
  store._state().finished.should.eql([]);
});

test('claim caps: global 2 host-wide, local 2, pan 1; a spec job is claimable while a pan job runs', async function () {
  const cfg = makeCfg({ global_max_jobs: 2, local_max_jobs: 2, pangenome_max_jobs: 1 });
  const shared = newShared();
  const now = clock();
  const a = createMemoryStore({ cfg: cfg, siteKey: 'site_a', shared: shared, now: now });
  const b = createMemoryStore({ cfg: cfg, siteKey: 'site_b', shared: shared, now: now });
  await a.submit(doc(1, 'pangenome', now()));
  await a.submit(doc(2, 'pangenome', now()));
  (await a.claim(now())).should.eql({ id: hexId(1), queue: 'pan' });
  should(await a.claim(now())).be.null(); // pan cap 1: the second pan job waits
  await a.submit(doc(3, 'specificity', now()));
  (await a.claim(now())).should.eql({ id: hexId(3), queue: 'spec' }); // spec job while the pan job runs
  await b.submit(doc(4, 'specificity', now()));
  should(await b.claim(now())).be.null(); // host-wide cap 2 reached by site_a
  a._state().slots.map(function (e) { return e[0]; }).sort().should.eql(['site_a:' + hexId(1), 'site_a:' + hexId(3)].sort());
  a._state().panSlots.map(function (e) { return e[0]; }).should.eql(['site_a:' + hexId(1)]);

  const started = await start(a, hexId(3), 'w1', now());
  (await a.complete(workerLib.doneDoc(started, { now: now() }), '{}', { worker: 'w1' })).should.be.true();
  (await b.claim(now())).should.eql({ id: hexId(4), queue: 'spec' }); // a slot was freed host-wide
  await a.submit(doc(5, 'specificity', now()));
  should(await a.claim(now())).be.null(); // global full again (site_a pan + site_b spec)

  // local cap on its own
  const cfgLocal = makeCfg({ global_max_jobs: 10, local_max_jobs: 2, pangenome_max_jobs: 1 });
  const c = createMemoryStore({ cfg: cfgLocal, siteKey: 'site_c', shared: newShared(), now: now });
  for (let i = 10; i < 13; i++) await c.submit(doc(i, 'specificity', now()));
  should.exist(await c.claim(now()));
  should.exist(await c.claim(now()));
  should(await c.claim(now())).be.null();
  c._state().queues.spec.should.eql([hexId(12)]);
});

test('stale host-wide slots are pruned after stale_ms; heartbeats keep a slot alive', async function () {
  const cfg = makeCfg({ global_max_jobs: 1, local_max_jobs: 2, stale_ms: 60000 });
  const shared = newShared();
  const now = clock();
  const t0 = now();
  const a = createMemoryStore({ cfg: cfg, siteKey: 'site_a', shared: shared, now: now });
  const b = createMemoryStore({ cfg: cfg, siteKey: 'site_b', shared: shared, now: now });
  await a.submit(doc(1, 'specificity', t0));
  await b.submit(doc(2, 'specificity', t0));
  await b.submit(doc(3, 'specificity', t0));
  (await a.claim(t0)).id.should.equal(hexId(1));
  should(await b.claim(t0 + 59999)).be.null();
  (await b.claim(t0 + 60000)).id.should.equal(hexId(2)); // site_a's slot is stale (a hard-killed worker)
  shared.slots.has('site_a:' + hexId(1)).should.be.false();
  a._state().running.map(function (e) { return e[0]; }).should.eql([hexId(1)]); // the site's own running set is untouched

  await start(b, hexId(2), 'wb', t0 + 60000);
  (await b.heartbeat(hexId(2), { queue: 'spec', worker: 'wb', now: t0 + 110000 })).should.be.true();
  await a.submit(doc(4, 'specificity', t0));
  should(await a.claim(t0 + 150000)).be.null(); // heartbeat at +110 s keeps it alive
  (await a.claim(t0 + 170000)).id.should.equal(hexId(4));
});

test('claim drops expired and non-queued ids from the queues', async function () {
  const cfg = makeCfg({ ttl_active_s: 10 });
  const now = clock();
  const store = createMemoryStore({ cfg: cfg, siteKey: 'site_a', shared: newShared(), now: now });
  await store.submit(doc(1, 'specificity', now()));
  await store.submit(doc(2, 'pangenome', now()));
  now.advance(11000);
  should(await store.claim(now())).be.null();
  store._state().queues.should.eql({ spec: [], pan: [] });
  await rejects(jobs.status(hexId(1), { cfg: cfg, store: store, log: quiet }), 404, 'UNKNOWN_JOB');
});

test('requeueRunning: running jobs back to the front; attempts >= max_attempts -> WORKER_LOST', async function () {
  const cfg = makeCfg({ max_attempts: 2, local_max_jobs: 5, global_max_jobs: 5 });
  const now = clock();
  const store = createMemoryStore({ cfg: cfg, siteKey: 'site_a', shared: newShared(), now: now });
  for (let i = 1; i <= 4; i++) await store.submit(doc(i, 'specificity', now()));
  await store.claim(now());
  await start(store, hexId(1), 'dead:1', now()); // attempt 1
  await store.claim(now());
  const second = Object.assign(await store.getJob(hexId(2)), { attempts: 1 });
  (await store.setJob(workerLib.startedDoc(second, { worker: 'dead:1', now: now() }), { expect: ['queued'] })).should.be.true(); // attempt 2
  await store.claim(now()); // hexId(3): claimed, never started
  await store.setPartial(hexId(1), '{"partial":true}', { worker: 'dead:1' });

  const worker = workerLib.createWorker({ store: store, cfg: cfg, siteKey: 'site_a', check: stubCheck(), log: quiet, now: now });
  const r = await worker.requeueRunning();
  r.requeued.slice().sort().should.eql([hexId(1), hexId(3)].sort());
  r.failed.should.eql([hexId(2)]);
  const st = store._state();
  st.queues.spec.slice(0, 2).sort().should.eql([hexId(1), hexId(3)].sort());
  st.queues.spec[2].should.equal(hexId(4));
  st.running.should.eql([]);
  st.slots.should.eql([]);
  const d = { cfg: cfg, store: store, log: quiet };
  const s1 = await jobs.status(hexId(1), d);
  s1.status.should.equal('queued');
  s1.attempts.should.equal(1);
  s1.partial.should.be.false();
  should(await store.getPartial(hexId(1))).be.null();
  s1.progress.stage.should.equal('queued');
  const s2 = await jobs.status(hexId(2), d);
  s2.status.should.equal('error');
  s2.error.code.should.equal('WORKER_LOST');
  s2.attempts.should.equal(2);
  (await jobs.status(hexId(3), d)).status.should.equal('queued');
  // nothing left to recover
  (await worker.requeueRunning()).should.eql({ requeued: [], failed: [], cleaned: [] });
});

test('queued-key TTL sweep refreshes queued jobs and drops ids of expired jobs', async function () {
  const cfg = makeCfg({ ttl_active_s: 100 });
  const now = clock();
  const store = createMemoryStore({ cfg: cfg, siteKey: 'site_a', shared: newShared(), now: now });
  await store.submit(doc(3, 'specificity', now()));
  now.advance(50000);
  await store.submit(doc(1, 'specificity', now()));
  await store.submit(doc(2, 'pangenome', now()));
  now.advance(60000); // job 3 has expired, jobs 1 and 2 have 40 s left
  store._ttlMs('job:' + hexId(1)).should.equal(40000);
  (await store.sweepQueued()).should.eql({ refreshed: 2, removed: 1 });
  store._ttlMs('job:' + hexId(1)).should.equal(100000);
  store._ttlMs('job:' + hexId(2)).should.equal(100000);
  store._state().queues.should.eql({ spec: [hexId(1)], pan: [hexId(2)] });
  await rejects(jobs.status(hexId(3), { cfg: cfg, store: store, log: quiet }), 404, 'UNKNOWN_JOB');
  // the worker's sweep() delegates to the store
  const worker = workerLib.createWorker({ store: store, cfg: cfg, siteKey: 'site_a', check: stubCheck(), log: quiet });
  (await worker.sweep()).should.eql({ refreshed: 2, removed: 0 });
});

test('progress and partial results are visible through status() while running; final results after completion', async function () {
  const cfg = makeCfg();
  const store = createMemoryStore({ cfg: cfg, siteKey: 'site_a', shared: newShared() });
  const d = { cfg: cfg, store: store, check: stubCheck(), log: quiet };
  const sub = await jobs.submit({ pairs: [{ id: 'A' }], kind: 'pangenome' }, d);
  const id = sub.job_id;
  await store.claim(Date.now());
  const started = await start(store, id, 'w1');

  let st = await jobs.status(id, d);
  st.status.should.equal('running');
  st.partial.should.be.false();
  should(st.results).be.null();
  should(st.queue_position).be.null();
  st.attempts.should.equal(1);
  st.started_at.should.be.a.String();
  st.request.should.eql(fakeNorm({ kind: 'pangenome', pairId: 'A' }).request);
  st.should.not.have.property('resolved');
  st.should.not.have.property('worker');

  const progress = { done: 2, total: 3, stage: 'pangenome', running: ['sorghum_rio'] };
  const withProgress = Object.assign({}, started, { progress: progress });
  (await store.setJob(withProgress, { expect: ['running'], worker: 'w1' })).should.be.true();
  (await store.setPartial(id, JSON.stringify({ specificity: { pairs: [] } }), { worker: 'w1' })).should.be.true();
  (await store.setPartial(id, '{"x":1}', { worker: 'intruder' })).should.be.false();
  st = await jobs.status(id, d);
  st.progress.should.eql(progress);
  st.partial.should.be.true();
  st.results.should.eql({ specificity: { pairs: [] } });

  (await store.complete(workerLib.doneDoc(withProgress, { now: Date.now() }), JSON.stringify({ final: true }), { worker: 'w1' })).should.be.true();
  st = await jobs.status(id, d);
  st.status.should.equal('done');
  st.partial.should.be.false();
  st.results.should.eql({ final: true });
  st.progress.should.eql({ done: 3, total: 3, stage: 'done', running: [] });
  st.finished_at.should.be.a.String();
  should(st.error).be.null();
  should(await store.getPartial(id)).be.null();
  (await store.getJob(id)).should.have.property('resolved', null);
  // late writes after completion are refused and leave no phantom running entry
  (await store.heartbeat(id, { queue: 'pan', worker: 'w1' })).should.be.false();
  (await store.setPartial(id, '{}', { worker: 'w1' })).should.be.false();
  (await store.complete(withProgress, '{}', { worker: 'w1' })).should.be.false();
  store._state().running.should.eql([]);
  store._state().panSlots.should.eql([]);
  const again = await jobs.submit({ pairs: [{ id: 'A' }], kind: 'pangenome' }, d);
  again.created.should.be.false();
  again.status.should.equal('done');
  should(again.queue_position).be.null();
});

test('finished jobs are trimmed to max_finished_jobs (oldest job and result keys deleted)', async function () {
  const cfg = makeCfg({ max_finished_jobs: 2, local_max_jobs: 5, global_max_jobs: 5 });
  const now = clock();
  const store = createMemoryStore({ cfg: cfg, siteKey: 'site_a', shared: newShared(), now: now });
  for (let i = 1; i <= 3; i++) {
    await store.submit(doc(i, 'specificity', now()));
    await store.claim(now());
    const started = await start(store, hexId(i), 'w', now());
    now.advance(1000);
    await store.complete(workerLib.doneDoc(started, { now: now() }), '{"n":' + i + '}', { worker: 'w', now: now() });
  }
  should(await store.getJob(hexId(1))).be.null();
  should(await store.getResults(hexId(1))).be.null();
  (await store.getResults(hexId(3))).should.eql({ n: 3 });
  store._state().finished.map(function (e) { return e[0]; }).should.eql([hexId(2), hexId(3)]);
  store._ttlMs('job:' + hexId(3)).should.equal(cfg.check.ttl_done_s * 1000);
  store._ttlMs('result:' + hexId(3)).should.equal(cfg.check.ttl_done_s * 1000);
});

test('status errors: malformed or unknown id -> 404 UNKNOWN_JOB; disabled -> 503 FEATURE_DISABLED', async function () {
  const cfg = makeCfg();
  const store = createMemoryStore({ cfg: cfg, siteKey: 'site_a', shared: newShared() });
  const e1 = await rejects(jobs.status('XYZ', { cfg: cfg, store: store }), 404, 'UNKNOWN_JOB');
  e1.details.should.eql({});
  const e2 = await rejects(jobs.status(hexId(99), { cfg: cfg, store: store }), 404, 'UNKNOWN_JOB');
  e2.details.should.eql({ job_id: hexId(99) });
  const off = makeCfg({}, { enabled: false });
  await rejects(jobs.status(hexId(1), { cfg: off, store: store }), 503, 'FEATURE_DISABLED');
  await rejects(jobs.submit({ pairs: [{ id: 'A' }] }, { cfg: off, store: store, check: stubCheck() }), 503, 'FEATURE_DISABLED');
});

test('a hung or failing store -> 503 JOB_STORE_UNAVAILABLE within the deadline (3 s by default)', async function () {
  jobs.STORE_DEADLINE_MS.should.equal(3000);
  const cfg = makeCfg();
  const store = createMemoryStore({ cfg: cfg, siteKey: 'site_a', shared: newShared() });
  const hung = Object.assign({}, store, { submit: function () { return new Promise(function () {}); } });
  const t0 = Date.now();
  const err = await rejects(jobs.submit({ pairs: [{ id: 'A' }] }, { cfg: cfg, store: hung, check: stubCheck(), storeDeadlineMs: 150, log: quiet }),
    503, 'JOB_STORE_UNAVAILABLE');
  (Date.now() - t0).should.be.below(1000);
  err.details.retry_after_s.should.equal(5);
  const broken = Object.assign({}, store, { getJob: async function () { throw new Error('connect ECONNREFUSED 127.0.0.1:6380'); } });
  const e2 = await rejects(jobs.status(hexId(1), { cfg: cfg, store: broken, log: quiet }), 503, 'JOB_STORE_UNAVAILABLE');
  e2.message.should.not.match(/6380/);
});

test('getStore() follows check.store and the site key; _resetForTests clears it', async function () {
  try {
    config._setForTests({ site_key: 'unit_memory_site', check: { store: 'memory' } });
    await jobs._resetForTests();
    const s = jobs.getStore();
    s.kind.should.equal('memory');
    s.siteKey.should.equal('unit_memory_site');
    jobs.getStore().should.equal(s);
    const override = createMemoryStore({ cfg: makeCfg(), siteKey: 'override', shared: newShared() });
    jobs._setStoreForTests(override);
    jobs.getStore().should.equal(override);
  } finally {
    await jobs._resetForTests();
    config._setForTests(null);
  }
  jobs.siteKeyLogLine('sorghum_v11:sorghum11').should.equal('primers site_key=sorghum_v11:sorghum11');
});
