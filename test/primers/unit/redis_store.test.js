'use strict';

// jobs/redis_store.js + lua.js against Redis 6380 db 1, confined to the dev-test prefixes
//   site keys   primers:<devtest>:*            (site keys '<devtest>' and '<devtest>:b')
//   host slots  primers:<devtest>_global:*     (check.global_prefix)
// <devtest> is PRIMERS_DEVTEST_PREFIX when it matches ^devtest_[a-z0-9_]+$ (so concurrent runs stay apart),
// else devtest_jobs. Every test starts from, and the file ends with, SCAN + DEL of exactly those patterns.
// Skipped when Redis is unreachable. The API-client failure tests use local ports and need no Redis.

const test = require('node:test');
const should = require('should');
const net = require('net');
const zlib = require('zlib');
const path = require('path');
const crypto = require('crypto');
const { spawn } = require('child_process');
const Redis = require('ioredis');

const config = require('../../../api/helpers/primers/config');
const jobs = require('../../../api/helpers/primers/jobs');
const { createRedisStore, keyLayout, API_CLIENT_OPTIONS } = require('../../../api/helpers/primers/jobs/redis_store');
const workerLib = require('../../../api/helpers/primers/jobs/worker');

const REDIS_URL = 'redis://localhost:6380/1';
const DEVTEST = /^devtest_[a-z0-9_]+$/.test(process.env.PRIMERS_DEVTEST_PREFIX || '') ? process.env.PRIMERS_DEVTEST_PREFIX : 'devtest_jobs';
const SITE = DEVTEST;
const SITE_B = DEVTEST + ':b';
const GLOBAL_PREFIX = 'primers:' + DEVTEST + '_global:';
const PATTERNS = ['primers:' + DEVTEST + ':*', 'primers:' + DEVTEST + '_global:*'];
const REALDATA = process.env.PRIMERS_REALDATA === '1';
const WORKER_MAIN = path.resolve(__dirname, '../../../api/helpers/primers/jobs/worker_main.js');
const CONFIG_JS = path.resolve(__dirname, '../../../api/helpers/primers/config.js');
const quiet = { info() {}, warn() {}, error() {}, log() {} };
const FAST = { claimIntervalMs: 10, lockRefreshMs: 50, lockRetryMs: 20, errorSleepMs: 20, shutdownGraceMs: 300, shutdownWriteMs: 150, drainMs: 500 };

let available = false;
let raw = null;
const stores = [];

function makeCfg(check) {
  return config._build({
    env: {},
    fileConfig: {},
    overrides: { check: Object.assign({ redis_url: REDIS_URL, global_prefix: GLOBAL_PREFIX }, check || {}) }
  }).config;
}

function mk(cfg, site, role) {
  const s = createRedisStore({ cfg: cfg, siteKey: site || SITE, role: role || 'api', log: quiet });
  stores.push(s);
  return s;
}

function hexId(n) {
  return crypto.createHash('md5').update('redis-job-' + n).digest('hex');
}

function fakeNorm(o) {
  o = o || {};
  const kind = o.kind || 'specificity';
  return {
    request: {
      system_name: 'sorghum_bicolor', mode: 'region',
      checks: kind === 'pangenome' ? ['pangenome', 'specificity'] : ['specificity'],
      genomes: kind === 'pangenome' ? ['sorghum_353'] : [],
      params: { max_product_size: 4000 },
      pairs: [{ id: o.pairId || 'P1', left: 'GGACAGCTCCACAACATATCAG', right: 'GGACATTTGAAGCCCATGGCC' }]
    },
    resolved: { assemblies: { sorghum_bicolor: { system_name: 'sorghum_bicolor', fingerprint: 'b'.repeat(40) } }, gene: null, species: null },
    kind: kind,
    warnings: [],
    estimate: { cpu_s: 8, total: kind === 'pangenome' ? 2 : 1 },
    dbs: { sorghum_bicolor: 'b'.repeat(40) }
  };
}

function doc(n, kind, now) {
  return jobs.newJobDoc(hexId(n), fakeNorm({ kind: kind, pairId: 'P' + n }), now || Date.now());
}

function stubCheck(run) {
  return {
    ALGORITHM_VERSION: '1',
    normalize: async function (body) { return fakeNorm({ pairId: body.pairs[0].id, kind: body.kind }); },
    run: run || async function () { return {}; }
  };
}

function sleep(ms) {
  return new Promise(function (r) { setTimeout(r, ms); });
}

async function waitFor(fn, timeoutMs, label) {
  const deadline = Date.now() + (timeoutMs || 3000);
  for (;;) {
    const v = await fn();
    if (v) return v;
    if (Date.now() > deadline) throw new Error('timed out waiting for ' + (label || 'condition'));
    await sleep(20);
  }
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

async function probe() {
  const c = new Redis(REDIS_URL, { lazyConnect: true, enableOfflineQueue: false, maxRetriesPerRequest: 0, connectTimeout: 1000, retryStrategy: function () { return null; } });
  c.on('error', function () {});
  try {
    await c.connect();
    return (await c.ping()) === 'PONG';
  } catch (e) {
    return false;
  } finally {
    c.disconnect();
  }
}

async function cleanup() {
  if (!raw) return 0;
  let deleted = 0;
  for (const pattern of PATTERNS) {
    let cursor = '0';
    do {
      const res = await raw.scan(cursor, 'MATCH', pattern, 'COUNT', 500);
      cursor = res[0];
      if (res[1].length) deleted += await raw.del.apply(raw, res[1]);
    } while (cursor !== '0');
  }
  return deleted;
}

test.before(async function () {
  available = await probe();
  if (!available) return;
  raw = new Redis(REDIS_URL);
  raw.on('error', function () {});
  await cleanup();
});

test.after(async function () {
  for (const s of stores) {
    try { await s.close(); } catch (e) { /* ignore */ }
  }
  if (raw) {
    await cleanup();
    let left = 0;
    for (const pattern of PATTERNS) {
      let cursor = '0';
      do {
        const res = await raw.scan(cursor, 'MATCH', pattern, 'COUNT', 500);
        cursor = res[0];
        left += res[1].length;
      } while (cursor !== '0');
    }
    await raw.quit();
    if (left) throw new Error(left + ' dev-test keys left in redis');
  }
});

function rtest(name, fn) {
  test(name, async function (t) {
    if (!available) {
      t.skip('redis unreachable at ' + REDIS_URL);
      return;
    }
    await cleanup();
    await fn(t);
  });
}

// ---- layout ----------------------------------------------------------------------------------------

test('key layout follows spec A.8.1 with the configured global prefix', function () {
  const K = keyLayout(SITE, GLOBAL_PREFIX);
  const site = 'primers:' + DEVTEST + ':';
  K.job('abc').should.equal(site + 'job:abc');
  K.partial('abc').should.equal(site + 'partial:abc');
  K.result('abc').should.equal(site + 'result:abc');
  K.queue('spec').should.equal(site + 'queue:spec');
  K.queue('pangenome').should.equal(site + 'queue:pan');
  K.running.should.equal(site + 'running');
  K.finished.should.equal(site + 'finished');
  K.worker.should.equal(site + 'worker');
  K.slots.should.equal('primers:' + DEVTEST + '_global:slots');
  K.panSlots.should.equal('primers:' + DEVTEST + '_global:pan_slots');
  API_CLIENT_OPTIONS.should.eql({ lazyConnect: true, enableOfflineQueue: false, maxRetriesPerRequest: 1, connectTimeout: 2000 });
  should.throws(function () { createRedisStore({ cfg: makeCfg(), siteKey: 'bad key' }); }, TypeError);
  makeCfg().check.global_prefix.should.equal(GLOBAL_PREFIX);
});

// ---- submit / queues -------------------------------------------------------------------------------

rtest('submit: QUEUED with ttl_active, EXISTS, FULL; queuePosition via LRANGE 0 49; queues have no TTL', async function () {
  const cfg = makeCfg({ max_queued: 3 });
  const s = mk(cfg);
  const K = s.keys;
  const d1 = doc(1);
  (await s.submit(d1)).outcome.should.equal('QUEUED');
  (await raw.ttl(K.job(d1.job_id))).should.be.within(cfg.check.ttl_active_s - 5, cfg.check.ttl_active_s);
  (await raw.lrange(K.queue('spec'), 0, -1)).should.eql([d1.job_id]);
  const again = await s.submit(doc(1, 'specificity', Date.now() + 5000));
  again.outcome.should.equal('EXISTS');
  again.job.created_at.should.equal(d1.created_at);
  (await raw.llen(K.queue('spec'))).should.equal(1);
  (await s.submit(doc(2, 'pangenome'))).outcome.should.equal('QUEUED');
  (await s.submit(doc(3))).outcome.should.equal('QUEUED');
  (await s.submit(doc(4))).outcome.should.equal('FULL');
  (await raw.exists(K.job(hexId(4)))).should.equal(0);
  (await s.queuePosition(hexId(3), 'spec')).should.equal(1);
  (await s.queuePosition(hexId(2), 'pan')).should.equal(0);
  (await s.queuePosition(hexId(9), 'spec')).should.equal(-1);
  (await s.getJob(d1.job_id)).should.eql(d1);
  (await raw.ttl(K.queue('spec'))).should.equal(-1);
  should(await s.getJob(hexId(9))).be.null();
});

rtest('jobs.submit/status through the redis store: 202-style then 200-style, QUEUE_FULL 503', async function () {
  const cfg = makeCfg({ max_queued: 1 });
  const s = mk(cfg);
  const d = { cfg: cfg, store: s, check: stubCheck(), log: quiet };
  const first = await jobs.submit({ pairs: [{ id: 'A' }] }, d);
  first.created.should.be.true();
  first.queue_position.should.equal(0);
  const second = await jobs.submit({ pairs: [{ id: 'A' }] }, d);
  second.created.should.be.false();
  second.job_id.should.equal(first.job_id);
  const err = await rejects(jobs.submit({ pairs: [{ id: 'B' }] }, d), 503, 'QUEUE_FULL');
  err.details.retry_after_s.should.equal(60);
  const st = await jobs.status(first.job_id, d);
  st.status.should.equal('queued');
  st.queue_position.should.equal(0);
  st.should.not.have.property('resolved');
  await rejects(jobs.status(hexId(123), d), 404, 'UNKNOWN_JOB');
});

rtest('an errored job is replaced and re-queued by submit (and leaves the finished set)', async function () {
  const cfg = makeCfg();
  const s = mk(cfg, SITE, 'worker');
  const K = s.keys;
  const d = doc(1);
  await s.submit(d);
  (await s.claim(Date.now())).id.should.equal(d.job_id);
  const started = workerLib.startedDoc(d, { worker: 'w', now: Date.now() });
  (await s.setJob(started, { expect: ['queued'] })).should.be.true();
  (await s.fail(workerLib.errorDoc(started, { code: 'X_FAILED', message: 'm' }, { now: Date.now() }), { worker: 'w' })).should.be.true();
  (await raw.ttl(K.job(d.job_id))).should.be.within(cfg.check.ttl_error_s - 5, cfg.check.ttl_error_s);
  (await raw.exists(K.result(d.job_id))).should.equal(0);
  should.exist(await raw.zscore(K.finished, d.job_id));
  (await s.submit(doc(1))).outcome.should.equal('QUEUED');
  (await s.getJob(d.job_id)).status.should.equal('queued');
  should(await raw.zscore(K.finished, d.job_id)).be.null();
  (await raw.lrange(K.queue('spec'), 0, -1)).should.eql([d.job_id]);
});

// ---- claim -----------------------------------------------------------------------------------------

rtest('claim caps: global 2 across sites, local 2, pan 1; a spec job is claimable while a pan job runs', async function () {
  const cfg = makeCfg({ global_max_jobs: 2, local_max_jobs: 2, pangenome_max_jobs: 1 });
  const a = mk(cfg, SITE, 'worker');
  const b = mk(cfg, SITE_B, 'worker');
  a.keys.slots.should.equal(b.keys.slots);
  const t = Date.now();
  await a.submit(doc(1, 'pangenome'));
  await a.submit(doc(2, 'pangenome'));
  (await a.claim(t)).should.eql({ id: hexId(1), queue: 'pan' });
  should(await a.claim(t)).be.null();
  await a.submit(doc(3));
  (await a.claim(t)).should.eql({ id: hexId(3), queue: 'spec' });
  await b.submit(doc(4));
  should(await b.claim(t)).be.null();
  (await raw.zrange(a.keys.slots, 0, -1)).sort().should.eql([SITE + ':' + hexId(1), SITE + ':' + hexId(3)].sort());
  (await raw.zrange(a.keys.panSlots, 0, -1)).should.eql([SITE + ':' + hexId(1)]);
  (await raw.zrange(a.keys.running, 0, -1)).sort().should.eql([hexId(1), hexId(3)].sort());
  Number(await raw.zscore(a.keys.running, hexId(3))).should.equal(t);

  const started = workerLib.startedDoc(doc(3, 'specificity', t), { worker: 'w', now: t });
  (await a.setJob(started, { expect: ['queued'] })).should.be.true();
  (await a.complete(workerLib.doneDoc(started, { now: t }), '{}', { worker: 'w', now: t })).should.be.true();
  (await b.claim(t)).should.eql({ id: hexId(4), queue: 'spec' });
  (await raw.zrange(b.keys.slots, 0, -1)).sort().should.eql([SITE + ':' + hexId(1), SITE_B + ':' + hexId(4)].sort());
  await a.submit(doc(5));
  should(await a.claim(t)).be.null(); // host-wide cap

  // local cap with room host-wide
  await cleanup();
  const local = makeCfg({ global_max_jobs: 10, local_max_jobs: 2 });
  const c = mk(local, SITE, 'worker');
  for (let i = 10; i < 13; i++) await c.submit(doc(i));
  should.exist(await c.claim(t));
  should.exist(await c.claim(t));
  should(await c.claim(t)).be.null();
  (await raw.lrange(c.keys.queue('spec'), 0, -1)).should.eql([hexId(12)]);
});

rtest('claim prunes stale host-wide slots after stale_ms; heartbeats keep a slot alive', async function () {
  const cfg = makeCfg({ global_max_jobs: 1, local_max_jobs: 2, stale_ms: 60000 });
  const a = mk(cfg, SITE, 'worker');
  const b = mk(cfg, SITE_B, 'worker');
  const t0 = Date.now();
  await a.submit(doc(1));
  await b.submit(doc(2));
  await a.submit(doc(4));
  (await a.claim(t0)).id.should.equal(hexId(1));
  should(await b.claim(t0 + 59999)).be.null();
  (await b.claim(t0 + 60000)).id.should.equal(hexId(2));
  should(await raw.zscore(a.keys.slots, SITE + ':' + hexId(1))).be.null();
  should.exist(await raw.zscore(a.keys.running, hexId(1))); // the site's running set is not pruned
  const started = workerLib.startedDoc(doc(2), { worker: 'wb', now: t0 + 60000 });
  (await b.setJob(started, { expect: ['queued'] })).should.be.true();
  (await b.heartbeat(hexId(2), { queue: 'spec', worker: 'wb', now: t0 + 110000 })).should.be.true();
  should(await a.claim(t0 + 150000)).be.null();
  (await a.claim(t0 + 170000)).id.should.equal(hexId(4));
});

rtest('claim drops ids whose job expired or is no longer queued', async function () {
  const s = mk(makeCfg(), SITE, 'worker');
  const K = s.keys;
  await raw.rpush(K.queue('spec'), hexId(1));
  await raw.set(K.job(hexId(2)), JSON.stringify({ job_id: hexId(2), status: 'running' }), 'EX', 60);
  await raw.rpush(K.queue('spec'), hexId(2));
  await raw.set(K.job(hexId(3)), 'not json', 'EX', 60);
  await raw.rpush(K.queue('pan'), hexId(3));
  await raw.set(K.job(hexId(5)), JSON.stringify({ job_id: hexId(5), status: null }), 'EX', 60);
  await raw.rpush(K.queue('pan'), hexId(5));
  should(await s.claim(Date.now())).be.null();
  (await raw.llen(K.queue('spec'))).should.equal(0);
  (await raw.llen(K.queue('pan'))).should.equal(0);
  (await raw.zcard(K.running)).should.equal(0);
});

// ---- job lifecycle ---------------------------------------------------------------------------------

rtest('lifecycle: CAS start, heartbeat, gzip partial, complete (TTLs, slots, finished), no late writes', async function () {
  const cfg = makeCfg();
  const s = mk(cfg, SITE, 'worker');
  const K = s.keys;
  const d = doc(1, 'pangenome');
  const id = d.job_id;
  const member = SITE + ':' + id;
  await s.submit(d);
  const t0 = Date.now();
  (await s.claim(t0)).should.eql({ id: id, queue: 'pan' });
  const started = workerLib.startedDoc(d, { worker: 'w', now: t0 });
  (await s.setJob(started, { expect: ['running'] })).should.be.false(); // still queued
  (await s.setJob(started, { expect: ['queued'] })).should.be.true();

  await raw.expire(K.job(id), 100);
  (await s.heartbeat(id, { queue: 'pan', worker: 'w', now: t0 + 5000 })).should.be.true();
  Number(await raw.zscore(K.running, id)).should.equal(t0 + 5000);
  Number(await raw.zscore(K.slots, member)).should.equal(t0 + 5000);
  Number(await raw.zscore(K.panSlots, member)).should.equal(t0 + 5000);
  (await raw.ttl(K.job(id))).should.be.within(cfg.check.ttl_active_s - 5, cfg.check.ttl_active_s);
  (await s.heartbeat(id, { queue: 'pan', worker: 'intruder', now: t0 + 6000 })).should.be.false();

  const partial = { specificity: { pairs: [{ id: 'P1', verdict: 'specific' }] } };
  (await s.setPartial(id, JSON.stringify(partial), { worker: 'w' })).should.be.true();
  const pbuf = await raw.getBuffer(K.partial(id));
  pbuf[0].should.equal(0x1f);
  pbuf[1].should.equal(0x8b);
  JSON.parse(zlib.gunzipSync(pbuf).toString()).should.eql(partial);
  (await raw.ttl(K.partial(id))).should.be.within(cfg.check.ttl_active_s - 5, cfg.check.ttl_active_s);
  (await s.getPartial(id)).should.eql(partial);
  (await s.setPartial(id, '{}', { worker: 'intruder' })).should.be.false();

  const progressed = Object.assign({}, started, { progress: { done: 1, total: 2, stage: 'pangenome', running: ['sorghum_353'] } });
  (await s.setJob(progressed, { expect: ['running'], worker: 'w' })).should.be.true();
  (await s.getJob(id)).progress.running.should.eql(['sorghum_353']); // arrays survive (never re-encoded in Lua)

  const results = { engine: { algorithm_version: '1' }, pangenome: { pairs: [] }, blob: 'x'.repeat(10000) };
  (await s.complete(workerLib.doneDoc(progressed, { now: t0 + 9000 }), JSON.stringify(results), { worker: 'w', now: t0 + 9000 })).should.be.true();
  (await raw.ttl(K.job(id))).should.be.within(cfg.check.ttl_done_s - 5, cfg.check.ttl_done_s);
  (await raw.ttl(K.result(id))).should.be.within(cfg.check.ttl_done_s - 5, cfg.check.ttl_done_s);
  (await raw.exists(K.partial(id))).should.equal(0);
  should(await raw.zscore(K.running, id)).be.null();
  should(await raw.zscore(K.slots, member)).be.null();
  should(await raw.zscore(K.panSlots, member)).be.null();
  Number(await raw.zscore(K.finished, id)).should.equal(t0 + 9000);
  const rbuf = await raw.getBuffer(K.result(id));
  rbuf[0].should.equal(0x1f);
  rbuf.length.should.be.below(1000);
  (await s.getResults(id)).should.eql(results);
  const done = await s.getJob(id);
  done.status.should.equal('done');
  should(done.resolved).be.null();
  done.progress.should.eql({ done: 2, total: 2, stage: 'done', running: [] });

  (await s.heartbeat(id, { queue: 'pan', worker: 'w', now: t0 + 10000 })).should.be.false();
  should(await raw.zscore(K.running, id)).be.null(); // no phantom running entry
  (await s.setPartial(id, '{}', { worker: 'w' })).should.be.false();
  (await raw.exists(K.partial(id))).should.equal(0);
  (await s.complete(progressed, '{}', { worker: 'w' })).should.be.false();
});

rtest('finished is trimmed to max_finished_jobs; a re-queued job with a stale finished entry is kept', async function () {
  const cfg = makeCfg({ max_finished_jobs: 2, local_max_jobs: 5, global_max_jobs: 5 });
  const s = mk(cfg, SITE, 'worker');
  const K = s.keys;
  const t = Date.now();
  const active = doc(9, 'pangenome');
  await s.submit(active);
  await raw.zadd(K.finished, t - 100000, active.job_id);
  for (let i = 1; i <= 3; i++) {
    const d = doc(i);
    await s.submit(d);
    (await s.claim(t)).id.should.equal(d.job_id);
    const st = workerLib.startedDoc(d, { worker: 'w', now: t });
    (await s.setJob(st, { expect: ['queued'] })).should.be.true();
    (await s.complete(workerLib.doneDoc(st, { now: t + i }), '{"n":' + i + '}', { worker: 'w', now: t + i })).should.be.true();
  }
  (await raw.exists(K.job(active.job_id))).should.equal(1);
  (await raw.exists(K.job(hexId(1)))).should.equal(0);
  (await raw.exists(K.result(hexId(1)))).should.equal(0);
  (await raw.zrange(K.finished, 0, -1)).should.eql([hexId(2), hexId(3)]);
  (await s.getResults(hexId(3))).should.eql({ n: 3 });
});

rtest('requeue: front (or back) of the queue, slots released, partial deleted, worker CAS', async function () {
  const cfg = makeCfg();
  const s = mk(cfg, SITE, 'worker');
  const K = s.keys;
  const d1 = doc(1);
  const d2 = doc(2);
  await s.submit(d1);
  await s.submit(d2);
  (await s.claim(Date.now())).id.should.equal(d1.job_id);
  const started = workerLib.startedDoc(d1, { worker: 'w', now: Date.now() });
  (await s.setJob(started, { expect: ['queued'] })).should.be.true();
  await s.setPartial(d1.job_id, '{"p":1}', { worker: 'w' });
  (await s.requeue(workerLib.requeuedDoc(started), { worker: 'intruder' })).should.be.false();
  (await s.requeue(workerLib.requeuedDoc(started), { worker: 'w', front: true })).should.be.true();
  (await raw.lrange(K.queue('spec'), 0, -1)).should.eql([d1.job_id, d2.job_id]);
  (await raw.zcard(K.running)).should.equal(0);
  (await raw.zcard(K.slots)).should.equal(0);
  (await raw.exists(K.partial(d1.job_id))).should.equal(0);
  const q = await s.getJob(d1.job_id);
  q.status.should.equal('queued');
  q.attempts.should.equal(1);
  (await raw.ttl(K.job(d1.job_id))).should.be.within(cfg.check.ttl_active_s - 5, cfg.check.ttl_active_s);

  (await s.claim(Date.now())).id.should.equal(d1.job_id);
  const again = workerLib.startedDoc(q, { worker: 'w', now: Date.now() });
  (await s.setJob(again, { expect: ['queued'] })).should.be.true();
  (await s.requeue(workerLib.requeuedDoc(again), { worker: 'w', front: false })).should.be.true();
  (await raw.lrange(K.queue('spec'), 0, -1)).should.eql([d2.job_id, d1.job_id]);
  (await s.getJob(d1.job_id)).attempts.should.equal(2);
});

rtest('sweepQueued refreshes queued job TTLs and drops ids of missing jobs', async function () {
  const cfg = makeCfg();
  const s = mk(cfg, SITE, 'worker');
  const K = s.keys;
  const d1 = doc(1);
  const d2 = doc(2, 'pangenome');
  await s.submit(d1);
  await s.submit(d2);
  await raw.expire(K.job(d1.job_id), 100);
  await raw.del(K.job(d2.job_id));
  await raw.rpush(K.queue('spec'), d1.job_id);
  (await s.sweepQueued()).should.eql({ refreshed: 2, removed: 1 });
  (await raw.ttl(K.job(d1.job_id))).should.be.within(cfg.check.ttl_active_s - 5, cfg.check.ttl_active_s);
  (await raw.lrange(K.queue('spec'), 0, -1)).should.eql([d1.job_id, d1.job_id]);
  (await raw.llen(K.queue('pan'))).should.equal(0);
});

rtest('worker lock: SET NX PX, refresh (1), re-acquire after expiry (2), lost to another token (0), compare-and-delete', async function () {
  const s = mk(makeCfg(), SITE, 'worker');
  const K = s.keys;
  (await s.acquireWorkerLock('tokA', 30000)).should.be.true();
  (await s.acquireWorkerLock('tokB', 30000)).should.be.false();
  (await raw.pttl(K.worker)).should.be.within(29000, 30000);
  (await s.refreshWorkerLock('tokB', 30000)).should.equal(0);
  await raw.pexpire(K.worker, 1000);
  (await s.refreshWorkerLock('tokA', 30000)).should.equal(1);
  (await raw.pttl(K.worker)).should.be.above(20000);
  await raw.del(K.worker);
  (await s.refreshWorkerLock('tokA', 30000)).should.equal(2);
  (await raw.get(K.worker)).should.equal('tokA');
  (await s.releaseWorkerLock('tokB')).should.be.false();
  (await s.releaseWorkerLock('tokA')).should.be.true();
  (await raw.exists(K.worker)).should.equal(0);
});

rtest('requeueRunning on redis: WORKER_LOST after max_attempts, requeue below it, phantom entries cleaned', async function () {
  const cfg = makeCfg({ max_attempts: 2, local_max_jobs: 5, global_max_jobs: 5 });
  const s = mk(cfg, SITE, 'worker');
  const K = s.keys;
  const d1 = doc(1);
  const d3 = doc(3);
  await s.submit(d1);
  await s.submit(d3);
  (await s.claim(Date.now())).id.should.equal(d1.job_id);
  (await s.setJob(workerLib.startedDoc(Object.assign({}, d1, { attempts: 1 }), { worker: 'dead', now: Date.now() }), { expect: ['queued'] })).should.be.true();
  (await s.claim(Date.now())).id.should.equal(d3.job_id);
  (await s.setJob(workerLib.startedDoc(d3, { worker: 'dead', now: Date.now() }), { expect: ['queued'] })).should.be.true();
  const phantom = hexId(77);
  await raw.zadd(K.running, Date.now(), phantom);
  await raw.zadd(K.slots, Date.now(), SITE + ':' + phantom);
  const w = workerLib.createWorker({ store: s, cfg: cfg, siteKey: SITE, check: stubCheck(), log: quiet });
  const r = await w.requeueRunning();
  r.should.eql({ requeued: [d3.job_id], failed: [d1.job_id], cleaned: [phantom] });
  const lost = await s.getJob(d1.job_id);
  lost.status.should.equal('error');
  lost.error.code.should.equal('WORKER_LOST');
  (await raw.ttl(K.job(d1.job_id))).should.be.within(cfg.check.ttl_error_s - 5, cfg.check.ttl_error_s);
  (await raw.zcard(K.running)).should.equal(0);
  (await raw.zcard(K.slots)).should.equal(0);
  (await raw.lrange(K.queue('spec'), 0, -1)).should.eql([d3.job_id]);
});

rtest('runningEntries returns heartbeat scores; staleOnly reconciliation skips fresh claims', async function () {
  const cfg = makeCfg({ local_max_jobs: 5, global_max_jobs: 5 });
  const s = mk(cfg, SITE, 'worker');
  const old = doc(1);
  const fresh = doc(2);
  await s.submit(old);
  await s.submit(fresh);
  const t0 = Date.now();
  (await s.claim(t0 - 120000)).id.should.equal(old.job_id);
  (await s.claim(t0)).id.should.equal(fresh.job_id);
  (await s.runningEntries()).should.eql([[old.job_id, t0 - 120000], [fresh.job_id, t0]]);
  const w = workerLib.createWorker({ store: s, cfg: cfg, siteKey: SITE, check: stubCheck(), log: quiet, now: function () { return t0; } });
  (await w.requeueRunning({ staleOnly: true })).should.eql({ requeued: [old.job_id], failed: [], cleaned: [] });
  (await raw.zrange(s.keys.running, 0, -1)).should.eql([fresh.job_id]);
  (await raw.lrange(s.keys.queue('spec'), 0, -1)).should.eql([old.job_id]);
});

// ---- worker against redis --------------------------------------------------------------------------

rtest('end to end: jobs.submit -> worker (redis store) -> partial visible -> done; identical re-POST shares the job', async function () {
  const cfg = makeCfg({ progress_min_interval_ms: 0, heartbeat_ms: 50 });
  const api = mk(cfg, SITE, 'api');
  const ws = mk(cfg, SITE, 'worker');
  const check = stubCheck(async function (request, ctx) {
    ctx.progress({ done: 0, total: 1, stage: 'reference', running: ['sorghum_bicolor'] });
    ctx.partial({ specificity: null, note: 'partial' });
    await sleep(300);
    return { specificity: { pairs: [{ id: request.pairs[0].id, verdict: 'specific' }] } };
  });
  const d = { cfg: cfg, store: api, check: check, log: quiet };
  const sub = await jobs.submit({ pairs: [{ id: 'E2E' }] }, d);
  sub.created.should.be.true();
  const w = workerLib.createWorker({ store: ws, cfg: cfg, siteKey: SITE, check: check, log: quiet, timing: FAST });
  w.start();
  try {
    const running = await waitFor(async function () {
      const x = await jobs.status(sub.job_id, d);
      return x.status === 'running' && x.partial && x;
    }, 5000, 'running with partial results');
    running.results.should.eql({ specificity: null, note: 'partial' });
    running.progress.running.should.eql(['sorghum_bicolor']);
    const done = await waitFor(async function () {
      const x = await jobs.status(sub.job_id, d);
      return x.status === 'done' && x;
    }, 5000, 'done');
    done.results.should.eql({ specificity: { pairs: [{ id: 'E2E', verdict: 'specific' }] } });
    done.attempts.should.equal(1);
    const again = await jobs.submit({ pairs: [{ id: 'E2E' }] }, d);
    again.created.should.be.false();
    again.status.should.equal('done');
  } finally {
    await w.stop();
  }
  (await raw.exists(ws.keys.worker)).should.equal(0);
});

rtest('worker shutdown on redis: the running job returns to the front within 1 s and the lock is released', async function () {
  const cfg = makeCfg({ local_max_jobs: 1 });
  const ws = mk(cfg, SITE, 'worker');
  const K = ws.keys;
  const check = stubCheck(function (request, ctx) {
    return new Promise(function (resolve, reject) {
      ctx.signal.addEventListener('abort', function () { reject(ctx.signal.reason); });
    });
  });
  const d1 = doc(1);
  const d2 = doc(2);
  await ws.submit(d1);
  await ws.submit(d2);
  const w = workerLib.createWorker({ store: ws, cfg: cfg, siteKey: SITE, check: check, log: quiet, timing: FAST });
  w.start();
  await waitFor(async function () { return (await ws.getJob(d1.job_id)).status === 'running'; }, 5000, 'running');
  const t0 = Date.now();
  await w.stop();
  (Date.now() - t0).should.be.below(1000);
  const j = await ws.getJob(d1.job_id);
  j.status.should.equal('queued');
  j.attempts.should.equal(1);
  (await raw.lrange(K.queue('spec'), 0, -1)).should.eql([d1.job_id, d2.job_id]);
  (await raw.zcard(K.running)).should.equal(0);
  (await raw.zcard(K.slots)).should.equal(0);
  (await raw.exists(K.worker)).should.equal(0);
});

function childEnv() {
  const env = Object.assign({}, process.env, {
    PRIMERS_SITE_KEY: SITE, PRIMERS_GLOBAL_PREFIX: GLOBAL_PREFIX, PRIMERS_REDIS_URL: REDIS_URL, PRIMERS_JOB_STORE: 'redis', PRIMERS_ENABLED: '1'
  });
  delete env.PRIMERS_GLOBAL_MAX_JOBS;
  delete env.PRIMERS_MAX_QUEUED;
  return env;
}

function track(child) {
  child.out = '';
  child.err = '';
  child.stdout.on('data', function (c) { child.out += c; });
  child.stderr.on('data', function (c) { child.err += c; });
  child.exited = new Promise(function (resolve) {
    child.on('exit', function (code, signal) { resolve({ code: code, signal: signal, at: Date.now() }); });
  });
  return child;
}

// worker_main.main() in a child process with a stub check (local_max_jobs 1).
function startStubWorker(mode) {
  const script = [
    "'use strict';",
    'const config = require(' + JSON.stringify(CONFIG_JS) + ');',
    "const cfg = config._build({ overrides: { check: { local_max_jobs: 1 } } }).config;",
    'const wm = require(' + JSON.stringify(WORKER_MAIN) + ');',
    'const mode = ' + JSON.stringify(mode) + ';',
    'const run = mode === "block"',
    '  ? function (req, ctx) { return new Promise(function (res, rej) { ctx.signal.addEventListener("abort", function () { rej(ctx.signal.reason); }); }); }',
    '  : async function (req, ctx) { return { finished_by: "restart", attempt: ctx.attempt }; };',
    'wm.main({ config: { get: function () { return cfg; }, siteKey: config.siteKey, basePath: config.basePath },',
    '  probeMongo: false, check: { ALGORITHM_VERSION: "1", run: run } })',
    '  .catch(function (e) { console.error(e); process.exit(1); });'
  ].join('\n');
  return track(spawn(process.execPath, ['-e', script], { cwd: '/', env: childEnv(), stdio: ['ignore', 'pipe', 'pipe'] }));
}

function killIfAlive(child) {
  if (child && child.exitCode === null && child.signalCode === null) child.kill('SIGKILL');
}

rtest('worker_main process: SIGTERM requeues the running job to the front, releases the lock, exits 0 within 1 s; a restart finishes it', async function () {
  const cfg = makeCfg();
  const s = mk(cfg, SITE, 'api');
  const K = s.keys;
  const d1 = doc(1);
  const d2 = doc(2);
  await s.submit(d1);
  await s.submit(d2);
  let child = startStubWorker('block');
  try {
    await waitFor(async function () { const j = await s.getJob(d1.job_id); return j && j.status === 'running'; }, 15000, 'job running in the child');
    child.out.split('\n')[0].should.equal('primers site_key=' + SITE);
    (await raw.exists(K.worker)).should.equal(1);
    (await s.getJob(d2.job_id)).status.should.equal('queued');
    const t0 = Date.now();
    child.kill('SIGTERM');
    const ex = await child.exited;
    const ms = ex.at - t0;
    console.log('# worker_main SIGTERM -> exit in ' + ms + ' ms (code ' + ex.code + ')');
    ex.code.should.equal(0);
    ms.should.be.below(1000);
    const j = await s.getJob(d1.job_id);
    j.status.should.equal('queued');
    j.attempts.should.equal(1);
    (await raw.lrange(K.queue('spec'), 0, -1)).should.eql([d1.job_id, d2.job_id]);
    (await raw.zcard(K.running)).should.equal(0);
    (await raw.zcard(K.slots)).should.equal(0);
    (await raw.exists(K.worker)).should.equal(0);

    child = startStubWorker('finish');
    const done = await waitFor(async function () {
      const x = await jobs.status(d1.job_id, { cfg: cfg, store: s, log: quiet });
      return x.status === 'done' && x;
    }, 15000, 'finished after restart');
    done.attempts.should.equal(2);
    done.results.should.eql({ finished_by: 'restart', attempt: 2 });
    await waitFor(async function () { return (await s.getJob(d2.job_id)).status === 'done'; }, 15000, 'second job');
    child.kill('SIGTERM');
    (await child.exited).code.should.equal(0);
  } finally {
    killIfAlive(child);
  }
});

rtest('real data: node worker_main.js probes mongo, logs site_key, takes the lock and exits 0 on SIGTERM', async function (t) {
  if (!REALDATA) {
    t.skip('set PRIMERS_REALDATA=1');
    return;
  }
  const K = keyLayout(SITE, GLOBAL_PREFIX);
  const child = track(spawn(process.execPath, [WORKER_MAIN], { cwd: '/', env: childEnv(), stdio: ['ignore', 'pipe', 'pipe'] }));
  try {
    await waitFor(function () { return child.out.indexOf('primers site_key=' + SITE) >= 0; }, 15000, 'site_key line');
    await waitFor(async function () { return (await raw.exists(K.worker)) === 1; }, 20000, 'worker lock (after the mongo probe)');
    const t0 = Date.now();
    child.kill('SIGTERM');
    const ex = await child.exited;
    console.log('# real worker_main SIGTERM -> exit in ' + (ex.at - t0) + ' ms');
    ex.code.should.equal(0);
    (ex.at - t0).should.be.below(1000);
    (await raw.exists(K.worker)).should.equal(0);
  } finally {
    killIfAlive(child);
  }
});

// ---- API client failure modes (no redis needed) ------------------------------------------------------

function freePort() {
  return new Promise(function (resolve, reject) {
    const srv = net.createServer();
    srv.on('error', reject);
    srv.listen(0, '127.0.0.1', function () {
      const port = srv.address().port;
      srv.close(function () { resolve(port); });
    });
  });
}

test('API client: connection refused -> 503 JOB_STORE_UNAVAILABLE immediately', async function () {
  const port = await freePort();
  const cfg = makeCfg({ redis_url: 'redis://127.0.0.1:' + port + '/1' });
  const s = createRedisStore({ cfg: cfg, siteKey: SITE, role: 'api', log: quiet });
  try {
    const t0 = Date.now();
    await rejects(jobs.submit({ pairs: [{ id: 'A' }] }, { cfg: cfg, store: s, check: stubCheck(), log: quiet }), 503, 'JOB_STORE_UNAVAILABLE');
    (Date.now() - t0).should.be.below(3000);
    await rejects(jobs.status(hexId(1), { cfg: cfg, store: s, log: quiet }), 503, 'JOB_STORE_UNAVAILABLE');
  } finally {
    await s.close();
  }
});

test('API client: a server that accepts but never answers -> 503 within the 3 s deadline', async function () {
  const sockets = [];
  const server = net.createServer(function (sock) {
    sockets.push(sock);
    sock.on('error', function () {});
  });
  await new Promise(function (r) { server.listen(0, '127.0.0.1', r); });
  const cfg = makeCfg({ redis_url: 'redis://127.0.0.1:' + server.address().port + '/1' });
  const s = createRedisStore({ cfg: cfg, siteKey: SITE, role: 'api', log: quiet });
  try {
    const t0 = Date.now();
    await rejects(jobs.submit({ pairs: [{ id: 'A' }] }, { cfg: cfg, store: s, check: stubCheck(), log: quiet }), 503, 'JOB_STORE_UNAVAILABLE');
    const ms = Date.now() - t0;
    console.log('# black-hole redis -> 503 after ' + ms + ' ms');
    ms.should.be.below(3100);
    const t1 = Date.now();
    await rejects(jobs.status(hexId(1), { cfg: cfg, store: s, log: quiet }), 503, 'JOB_STORE_UNAVAILABLE');
    (Date.now() - t1).should.be.below(3100);
  } finally {
    await s.close();
    sockets.forEach(function (x) { x.destroy(); });
    await new Promise(function (r) { server.close(r); });
  }
});
