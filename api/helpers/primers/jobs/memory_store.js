'use strict';

// In-memory job store with the same interface and semantics as redis_store.js (tests, single-process dev).
// TTLs are emulated against an injectable clock. Host-wide slot maps are shared between stores created
// with the same `shared` object (default: one per global_prefix in this process).
//
// Interface (all methods return promises):
//   submit(doc) -> {outcome: 'QUEUED'|'EXISTS'|'FULL', job}
//   getJob(id) -> doc | null;  getResults(id) / getPartial(id) -> object | null
//   queuePosition(id, queue) -> index within the first 50 ids of queue 'spec'|'pan', or -1
//   acquireWorkerLock(token, ttlMs) -> bool;  refreshWorkerLock(token, ttlMs) -> 1 | 2 (re-acquired) | 0 (lost)
//   releaseWorkerLock(token) -> bool
//   claim(nowMs) -> {id, queue: 'spec'|'pan'} | null
//   setJob(doc, {expect: [statuses], worker}) -> bool            (compare-and-set, TTL ttl_active_s)
//   heartbeat(id, {queue, worker, now}) -> bool
//   setPartial(id, json, {worker}) -> bool
//   complete(doc, resultJson, {worker, expect, now}) -> bool      (TTL ttl_done_s, trims finished)
//   fail(doc, {worker, expect, now}) -> bool                      (TTL ttl_error_s, trims finished)
//   requeue(doc, {worker, expect, front}) -> bool                 (back on its queue, slots released)
//   release(id);  runningIds() -> [id];  runningEntries() -> [[id, heartbeatMs]]
//   sweepQueued() -> {refreshed, removed};  close()

const QUEUE_POSITION_WINDOW = 50;
const SHARED = new Map();

function queueOf(kind) {
  return kind === 'pangenome' || kind === 'pan' ? 'pan' : 'spec';
}

function safeParse(raw) {
  if (raw === null || raw === undefined) return null;
  try {
    return JSON.parse(raw);
  } catch (e) {
    return null;
  }
}

function sortedByScore(map) {
  return Array.from(map.entries()).sort(function (a, b) {
    return a[1] - b[1] || (a[0] < b[0] ? -1 : a[0] > b[0] ? 1 : 0);
  });
}

function createMemoryStore(opts) {
  opts = opts || {};
  const cfg = opts.cfg || require('../config').get();
  const ccfg = cfg.check;
  const siteKey = String(opts.siteKey || 'memory');
  const now = opts.now || Date.now;
  const globalPrefix = ccfg.global_prefix || 'primers:global:';
  let shared = opts.shared;
  if (!shared) {
    shared = SHARED.get(globalPrefix);
    if (!shared) {
      shared = { slots: new Map(), panSlots: new Map() };
      SHARED.set(globalPrefix, shared);
    }
  }

  const kv = new Map();
  const queues = { spec: [], pan: [] };
  const running = new Map();
  const finished = new Map();
  let lock = null;

  function getRaw(key) {
    const e = kv.get(key);
    if (!e) return null;
    if (e.expiresAt !== null && e.expiresAt <= now()) {
      kv.delete(key);
      return null;
    }
    return e.value;
  }

  function setRaw(key, value, ttlS) {
    kv.set(key, { value: value, expiresAt: ttlS === null ? null : now() + ttlS * 1000 });
  }

  function expire(key, ttlS) {
    if (getRaw(key) === null) return false;
    kv.get(key).expiresAt = now() + ttlS * 1000;
    return true;
  }

  function statusOf(id) {
    const doc = safeParse(getRaw('job:' + id));
    if (!doc || typeof doc !== 'object') return { status: null, worker: null };
    return {
      status: typeof doc.status === 'string' ? doc.status : null,
      worker: typeof doc.worker === 'string' ? doc.worker : null
    };
  }

  function matches(id, expect, worker) {
    const s = statusOf(id);
    if (s.status === null || (expect || ['running']).indexOf(s.status) < 0) return false;
    return !worker || s.worker === worker;
  }

  function member(id) {
    return siteKey + ':' + id;
  }

  function removeFromQueue(q, id) {
    for (let i = q.length - 1; i >= 0; i--) if (q[i] === id) q.splice(i, 1);
  }

  function lockAlive() {
    if (lock && lock.expiresAt <= now()) lock = null;
    return lock;
  }

  function finish(doc, resultJson, ttlS, o) {
    o = o || {};
    const id = doc.job_id;
    if (!matches(id, o.expect || ['running'], o.worker)) return false;
    if (resultJson !== null && resultJson !== undefined) setRaw('result:' + id, String(resultJson), ttlS);
    else kv.delete('result:' + id);
    setRaw('job:' + id, JSON.stringify(doc), ttlS);
    kv.delete('partial:' + id);
    running.delete(id);
    shared.slots.delete(member(id));
    shared.panSlots.delete(member(id));
    finished.set(id, o.now !== undefined ? o.now : now());
    const maxFinished = Math.max(1, Number(ccfg.max_finished_jobs) || 1);
    const excess = finished.size - maxFinished;
    if (excess > 0) {
      sortedByScore(finished).slice(0, excess).forEach(function (entry) {
        const oid = entry[0];
        const st = statusOf(oid).status;
        if (st === null || st === 'done' || st === 'error') {
          kv.delete('job:' + oid);
          kv.delete('result:' + oid);
        }
        finished.delete(oid);
      });
    }
    return true;
  }

  const store = {
    kind: 'memory',
    siteKey: siteKey,

    async submit(doc) {
      const id = doc.job_id;
      const cur = safeParse(getRaw('job:' + id));
      if (cur && typeof cur.status === 'string' && cur.status !== 'error') return { outcome: 'EXISTS', job: cur };
      if (queues.spec.length + queues.pan.length >= ccfg.max_queued) return { outcome: 'FULL', job: null };
      const json = JSON.stringify(doc);
      setRaw('job:' + id, json, ccfg.ttl_active_s);
      const q = queues[queueOf(doc.kind)];
      removeFromQueue(q, id);
      q.push(id);
      finished.delete(id);
      return { outcome: 'QUEUED', job: JSON.parse(json) };
    },

    async getJob(id) {
      return safeParse(getRaw('job:' + id));
    },

    async getResults(id) {
      return safeParse(getRaw('result:' + id));
    },

    async getPartial(id) {
      return safeParse(getRaw('partial:' + id));
    },

    async queuePosition(id, queue) {
      return queues[queueOf(queue)].slice(0, QUEUE_POSITION_WINDOW).indexOf(id);
    },

    async acquireWorkerLock(token, ttlMs) {
      if (lockAlive()) return false;
      lock = { token: token, expiresAt: now() + ttlMs };
      return true;
    },

    async refreshWorkerLock(token, ttlMs) {
      const cur = lockAlive();
      if (cur && cur.token === token) {
        cur.expiresAt = now() + ttlMs;
        return 1;
      }
      if (!cur) {
        lock = { token: token, expiresAt: now() + ttlMs };
        return 2;
      }
      return 0;
    },

    async releaseWorkerLock(token) {
      const cur = lockAlive();
      if (cur && cur.token === token) {
        lock = null;
        return true;
      }
      return false;
    },

    async claim(nowMs) {
      const t = nowMs === undefined ? now() : nowMs;
      const cutoff = t - ccfg.stale_ms;
      [shared.slots, shared.panSlots].forEach(function (m) {
        m.forEach(function (score, k) { if (score <= cutoff) m.delete(k); });
      });
      if (shared.slots.size >= ccfg.global_max_jobs || running.size >= ccfg.local_max_jobs) return null;
      function pop(q) {
        while (q.length) {
          const id = q.shift();
          if (statusOf(id).status === 'queued') return id;
        }
        return null;
      }
      let id = pop(queues.spec);
      let queue = 'spec';
      if (!id && shared.panSlots.size < ccfg.pangenome_max_jobs) {
        id = pop(queues.pan);
        queue = 'pan';
      }
      if (!id) return null;
      running.set(id, t);
      shared.slots.set(member(id), t);
      if (queue === 'pan') shared.panSlots.set(member(id), t);
      return { id: id, queue: queue };
    },

    async setJob(doc, o) {
      o = o || {};
      if (!matches(doc.job_id, o.expect || ['running'], o.worker)) return false;
      setRaw('job:' + doc.job_id, JSON.stringify(doc), ccfg.ttl_active_s);
      return true;
    },

    async heartbeat(id, o) {
      o = o || {};
      if (!matches(id, ['running'], o.worker)) return false;
      const t = o.now !== undefined ? o.now : now();
      running.set(id, t);
      shared.slots.set(member(id), t);
      if (queueOf(o.queue) === 'pan') shared.panSlots.set(member(id), t);
      expire('job:' + id, ccfg.ttl_active_s);
      return true;
    },

    async setPartial(id, json, o) {
      o = o || {};
      if (!matches(id, ['running'], o.worker)) return false;
      setRaw('partial:' + id, String(json), ccfg.ttl_active_s);
      return true;
    },

    async complete(doc, resultJson, o) {
      return finish(doc, resultJson === undefined ? 'null' : resultJson, ccfg.ttl_done_s, o);
    },

    async fail(doc, o) {
      return finish(doc, null, ccfg.ttl_error_s, o);
    },

    async requeue(doc, o) {
      o = o || {};
      const id = doc.job_id;
      if (!matches(id, o.expect || ['running'], o.worker)) return false;
      setRaw('job:' + id, JSON.stringify(doc), ccfg.ttl_active_s);
      const q = queues[queueOf(doc.kind)];
      removeFromQueue(q, id);
      if (o.front === false) q.push(id);
      else q.unshift(id);
      running.delete(id);
      shared.slots.delete(member(id));
      shared.panSlots.delete(member(id));
      kv.delete('partial:' + id);
      return true;
    },

    async release(id) {
      running.delete(id);
      shared.slots.delete(member(id));
      shared.panSlots.delete(member(id));
    },

    async runningIds() {
      return sortedByScore(running).map(function (e) { return e[0]; });
    },

    // [[id, heartbeat ms]] by score
    async runningEntries() {
      return sortedByScore(running);
    },

    async sweepQueued() {
      let refreshed = 0;
      let removed = 0;
      ['spec', 'pan'].forEach(function (name) {
        const keep = [];
        queues[name].forEach(function (id) {
          if (statusOf(id).status === 'queued') {
            expire('job:' + id, ccfg.ttl_active_s);
            refreshed++;
            keep.push(id);
          } else {
            removed++;
          }
        });
        queues[name] = keep;
      });
      return { refreshed: refreshed, removed: removed };
    },

    async close() {},

    // Test introspection.
    _state() {
      return {
        queues: { spec: queues.spec.slice(), pan: queues.pan.slice() },
        running: sortedByScore(running),
        finished: sortedByScore(finished),
        slots: sortedByScore(shared.slots),
        panSlots: sortedByScore(shared.panSlots),
        lock: lockAlive() ? Object.assign({}, lock) : null
      };
    },

    // Remaining TTL in ms of 'job:<id>', 'partial:<id>' or 'result:<id>' (null when missing, -1 without TTL).
    _ttlMs(key) {
      if (getRaw(key) === null) return null;
      const e = kv.get(key);
      return e.expiresAt === null ? -1 : e.expiresAt - now();
    }
  };
  return store;
}

module.exports = {
  createMemoryStore,
  queueOf,
  QUEUE_POSITION_WINDOW
};
