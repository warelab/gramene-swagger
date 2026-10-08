'use strict';

// Redis job store (spec §A.8.1 layout, §A.8.2 clients). Same interface as memory_store.js.
//
// Keys, P = primers:<site_key>:
//   P job:<id>        string JSON job doc          queued/running ttl_active_s (refreshed); done ttl_done_s; error ttl_error_s
//   P partial:<id>    gzip JSON partial results    ttl_active_s, deleted on completion
//   P result:<id>     gzip JSON final results      ttl_done_s
//   P queue:spec|pan  list of job ids              no TTL (ids only)
//   P running         zset id -> heartbeat ms
//   P finished        zset id -> finish ms, trimmed to max_finished_jobs by deleting the oldest job/result keys
//   P worker          worker lock (SET NX PX 30000, refreshed every 10 s)
//   <global_prefix>slots, <global_prefix>pan_slots   host-wide zsets '<site_key>:<id>' -> heartbeat ms
//
// Clients:
//   role 'api'    lazyConnect, enableOfflineQueue false, maxRetriesPerRequest 1, connectTimeout 2000; connect()
//                 once and reset on failure. Callers (jobs/index.js) add the 3 s deadline -> 503 JOB_STORE_UNAVAILABLE.
//   role 'worker' ioredis defaults (auto reconnect, offline queue).

const zlib = require('zlib');
const util = require('util');
const lua = require('./lua');
const { redactPaths } = require('../errors');

const gzip = util.promisify(zlib.gzip);
const gunzip = util.promisify(zlib.gunzip);

const API_CLIENT_OPTIONS = Object.freeze({ lazyConnect: true, enableOfflineQueue: false, maxRetriesPerRequest: 1, connectTimeout: 2000 });
const WORKER_CLIENT_OPTIONS = Object.freeze({ lazyConnect: true });
const CONNECT_DEADLINE_MS = 2500;
const CLOSE_TIMEOUT_MS = 1000;
const QUEUE_POSITION_WINDOW = 50;
const ERROR_LOG_INTERVAL_MS = 30000;
const SITE_KEY_RE = /^[A-Za-z0-9_.:-]{1,128}$/;

function queueOf(kind) {
  return kind === 'pangenome' || kind === 'pan' ? 'pan' : 'spec';
}

function keyLayout(siteKey, globalPrefix) {
  const P = 'primers:' + siteKey + ':';
  return Object.freeze({
    prefix: P,
    job: function (id) { return P + 'job:' + id; },
    partial: function (id) { return P + 'partial:' + id; },
    result: function (id) { return P + 'result:' + id; },
    queue: function (q) { return P + 'queue:' + queueOf(q); },
    running: P + 'running',
    finished: P + 'finished',
    worker: P + 'worker',
    jobPrefix: P + 'job:',
    resultPrefix: P + 'result:',
    slots: globalPrefix + 'slots',
    panSlots: globalPrefix + 'pan_slots'
  });
}

function withTimer(promise, ms, message) {
  let timer;
  const timeout = new Promise(function (resolve, reject) {
    timer = setTimeout(function () { reject(new Error(message)); }, ms);
    if (timer.unref) timer.unref();
  });
  return Promise.race([promise, timeout]).finally(function () { clearTimeout(timer); });
}

function parseJson(raw) {
  if (raw === null || raw === undefined) return null;
  return JSON.parse(raw);
}

async function gunzipJson(buf) {
  if (buf === null || buf === undefined) return null;
  return JSON.parse((await gunzip(buf)).toString('utf8'));
}

function createRedisStore(opts) {
  opts = opts || {};
  const cfg = opts.cfg || require('../config').get();
  const ccfg = cfg.check;
  const siteKey = opts.siteKey;
  if (typeof siteKey !== 'string' || !SITE_KEY_RE.test(siteKey)) {
    throw new TypeError('redis job store: siteKey must match ' + SITE_KEY_RE);
  }
  const role = opts.role === 'worker' ? 'worker' : 'api';
  const url = opts.url || ccfg.redis_url;
  const log = opts.log || console;
  const K = keyLayout(siteKey, ccfg.global_prefix || 'primers:global:');
  const maxFinished = Math.max(1, Number(ccfg.max_finished_jobs) || 1);

  let client = null;
  let connecting = null;
  let closed = false;
  let lastErrorLog = 0;

  function makeClient() {
    const Redis = require('ioredis');
    const c = new Redis(url, Object.assign({}, role === 'api' ? API_CLIENT_OPTIONS : WORKER_CLIENT_OPTIONS));
    c.on('error', function (err) {
      const t = Date.now();
      if (t - lastErrorLog < ERROR_LOG_INTERVAL_MS) return;
      lastErrorLog = t;
      try {
        // A refused connection to a dual-stack host is an AggregateError with an empty message.
        const detail = err && (err.message || err.code ||
          (Array.isArray(err.errors) && err.errors.map(function (e) { return e && (e.message || e.code); }).join('; ')));
        (log.error || log.log).call(log, 'primers job store (' + role + ') redis error: ' + redactPaths(String(detail || err)));
      } catch (e) { /* ignore */ }
    });
    lua.define(c);
    return c;
  }

  async function ready() {
    if (closed) throw new Error('the job store is closed');
    if (!client) client = makeClient();
    const c = client;
    if (role === 'worker') return c;
    if (c.status === 'ready') return c;
    if (c.status === 'end') {
      client = null;
      return ready();
    }
    if (c.status === 'wait' && !connecting) {
      connecting = withTimer(c.connect(), CONNECT_DEADLINE_MS, 'redis connect timed out').then(
        function () { connecting = null; },
        function (err) {
          connecting = null;
          if (client === c) client = null;
          try { c.disconnect(); } catch (e) { /* ignore */ }
          throw err;
        });
    }
    if (connecting) {
      await connecting;
      if (c.status === 'ready') return c;
    }
    throw new Error('redis connection is ' + c.status);
  }

  function member(id) {
    return siteKey + ':' + id;
  }

  function csv(list, dflt) {
    return (Array.isArray(list) && list.length ? list : dflt).join(',');
  }

  async function finish(doc, resultJson, ttlS, o) {
    o = o || {};
    const c = await ready();
    const id = doc.job_id;
    const has = resultJson !== null && resultJson !== undefined;
    const buf = has ? await gzip(Buffer.from(String(resultJson), 'utf8')) : '';
    const r = await c.primersFinish(K.job(id), K.result(id), K.partial(id), K.running, K.slots, K.panSlots, K.finished,
      id, member(id), JSON.stringify(doc), ttlS, o.now !== undefined ? o.now : Date.now(), buf, has ? '1' : '0',
      o.worker || '', maxFinished, K.jobPrefix, K.resultPrefix, csv(o.expect, ['running']));
    return Number(r[0]) === 1;
  }

  return {
    kind: 'redis',
    role: role,
    siteKey: siteKey,
    keys: K,

    async submit(doc) {
      const c = await ready();
      const id = doc.job_id;
      const q = queueOf(doc.kind);
      const r = await c.primersSubmit(K.job(id), K.queue(q), K.queue(q === 'pan' ? 'spec' : 'pan'), K.finished,
        JSON.stringify(doc), id, ccfg.ttl_active_s, ccfg.max_queued);
      if (r[0] === 'EXISTS') return { outcome: 'EXISTS', job: parseJson(r[1]) };
      if (r[0] === 'FULL') return { outcome: 'FULL', job: null };
      return { outcome: 'QUEUED', job: doc };
    },

    async getJob(id) {
      const c = await ready();
      return parseJson(await c.get(K.job(id)));
    },

    async getResults(id) {
      const c = await ready();
      return gunzipJson(await c.getBuffer(K.result(id)));
    },

    async getPartial(id) {
      const c = await ready();
      return gunzipJson(await c.getBuffer(K.partial(id)));
    },

    async queuePosition(id, queue) {
      const c = await ready();
      return (await c.lrange(K.queue(queue), 0, QUEUE_POSITION_WINDOW - 1)).indexOf(id);
    },

    async acquireWorkerLock(token, ttlMs) {
      const c = await ready();
      return (await c.set(K.worker, token, 'PX', ttlMs, 'NX')) === 'OK';
    },

    async refreshWorkerLock(token, ttlMs) {
      const c = await ready();
      return Number(await c.primersLockRefresh(K.worker, token, ttlMs));
    },

    async releaseWorkerLock(token) {
      const c = await ready();
      return Number(await c.primersLockRelease(K.worker, token)) === 1;
    },

    async claim(nowMs) {
      const c = await ready();
      const r = await c.primersClaim(K.queue('spec'), K.queue('pan'), K.running, K.slots, K.panSlots,
        nowMs === undefined ? Date.now() : nowMs, ccfg.stale_ms, ccfg.global_max_jobs, ccfg.local_max_jobs,
        ccfg.pangenome_max_jobs, siteKey, K.jobPrefix);
      return Array.isArray(r) ? { id: r[0], queue: r[1] } : null;
    },

    async setJob(doc, o) {
      o = o || {};
      const c = await ready();
      return Number(await c.primersSetJob(K.job(doc.job_id), JSON.stringify(doc), ccfg.ttl_active_s,
        csv(o.expect, ['running']), o.worker || '')) === 1;
    },

    async heartbeat(id, o) {
      o = o || {};
      const c = await ready();
      return Number(await c.primersHeartbeat(K.job(id), K.running, K.slots, K.panSlots, id, member(id),
        o.now !== undefined ? o.now : Date.now(), ccfg.ttl_active_s, queueOf(o.queue) === 'pan' ? '1' : '0',
        o.worker || '')) === 1;
    },

    async setPartial(id, json, o) {
      o = o || {};
      const c = await ready();
      const buf = await gzip(Buffer.from(String(json), 'utf8'));
      return Number(await c.primersSetPartial(K.job(id), K.partial(id), buf, ccfg.ttl_active_s, o.worker || '')) === 1;
    },

    async complete(doc, resultJson, o) {
      return finish(doc, resultJson === undefined ? 'null' : resultJson, ccfg.ttl_done_s, o);
    },

    async fail(doc, o) {
      return finish(doc, null, ccfg.ttl_error_s, o);
    },

    async requeue(doc, o) {
      o = o || {};
      const c = await ready();
      const id = doc.job_id;
      return Number(await c.primersRequeue(K.job(id), K.queue(doc.kind), K.running, K.slots, K.panSlots, K.partial(id),
        id, member(id), JSON.stringify(doc), ccfg.ttl_active_s, o.worker || '', o.front === false ? '0' : '1',
        csv(o.expect, ['running']))) === 1;
    },

    async release(id) {
      const c = await ready();
      await c.multi().zrem(K.running, id).zrem(K.slots, member(id)).zrem(K.panSlots, member(id)).exec();
    },

    async runningIds() {
      const c = await ready();
      return c.zrange(K.running, 0, -1);
    },

    // [[id, heartbeat ms]] by score
    async runningEntries() {
      const c = await ready();
      const flat = await c.zrange(K.running, 0, -1, 'WITHSCORES');
      const out = [];
      for (let i = 0; i + 1 < flat.length; i += 2) out.push([flat[i], Number(flat[i + 1])]);
      return out;
    },

    async sweepQueued() {
      const c = await ready();
      const r = await c.primersSweep(K.queue('spec'), K.queue('pan'), K.jobPrefix, ccfg.ttl_active_s);
      return { refreshed: Number(r[0]), removed: Number(r[1]) };
    },

    async close() {
      closed = true;
      const c = client;
      client = null;
      if (!c) return;
      try {
        if (c.status === 'ready') await withTimer(c.quit(), CLOSE_TIMEOUT_MS, 'redis quit timed out');
        else c.disconnect();
      } catch (e) {
        try { c.disconnect(); } catch (x) { /* ignore */ }
      }
    },

    // Test hook: the connected ioredis client.
    _client: ready
  };
}

module.exports = {
  createRedisStore,
  keyLayout,
  queueOf,
  API_CLIENT_OPTIONS,
  WORKER_CLIENT_OPTIONS,
  QUEUE_POSITION_WINDOW
};
