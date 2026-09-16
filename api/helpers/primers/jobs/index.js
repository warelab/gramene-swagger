'use strict';

// Check jobs, API side (spec §A.2.3, §A.2.4, §A.8).
//   submit(body, deps) -> {job_id, status, kind, queue_position, progress, estimate {cpu_s}, created_at, warnings, created}
//                         created true => HTTP 202, false => 200 (an errored job is re-queued and counts as created)
//   status(job_id, deps) -> {job_id, status, kind, partial, queue_position, progress, created_at, started_at, finished_at,
//                            attempts, request, warnings, estimate, results, error}   (never `resolved` or `worker`)
// Errors: normalize's 4xx/503; 503 QUEUE_FULL (retry_after_s 60); 503 JOB_STORE_UNAVAILABLE (store failure or
// no answer within 3 s); 503 FEATURE_DISABLED; 404 UNKNOWN_JOB.
// deps: {cfg, store, check, now, storeDeadlineMs, log, ...normalize deps}

const crypto = require('crypto');
const { PrimerHttpError, isPrimerHttpError, redactPaths } = require('../errors');

const JOB_ID_RE = /^[0-9a-f]{32}$/;
const JOB_ID_VERSION = 1;
const STORE_DEADLINE_MS = 3000;
const QUEUE_FULL_RETRY_AFTER_S = 60;
const STORE_RETRY_AFTER_S = 5;
const FEATURE_DISABLED_RETRY_AFTER_S = 300;

// JSON with object keys sorted at every level; undefined/function members are skipped like JSON.stringify.
function canonicalJSON(value) {
  if (value === null || typeof value !== 'object') {
    if (typeof value === 'number' && !Number.isFinite(value)) return 'null';
    if (value === undefined || typeof value === 'function' || typeof value === 'symbol') return undefined;
    return JSON.stringify(value);
  }
  if (typeof value.toJSON === 'function') return canonicalJSON(value.toJSON());
  if (Array.isArray(value)) {
    return '[' + value.map(function (v) {
      const s = canonicalJSON(v);
      return s === undefined ? 'null' : s;
    }).join(',') + ']';
  }
  const parts = [];
  Object.keys(value).sort().forEach(function (k) {
    const s = canonicalJSON(value[k]);
    if (s !== undefined) parts.push(JSON.stringify(k) + ':' + s);
  });
  return '{' + parts.join(',') + '}';
}

// sha256(canonicalJSON({v: 1, algo, request: normalized, dbs: {system_name: fingerprint}})).hex[0:32]
function jobId(request, dbs, algorithmVersion) {
  const algo = algorithmVersion !== undefined ? String(algorithmVersion) : require('../check').ALGORITHM_VERSION;
  const text = canonicalJSON({ v: JOB_ID_VERSION, algo: algo, request: request, dbs: dbs || {} });
  return crypto.createHash('sha256').update(text).digest('hex').slice(0, 32);
}

function dbsOf(resolved) {
  const out = {};
  const asms = (resolved && resolved.assemblies) || {};
  Object.keys(asms).forEach(function (sys) {
    out[sys] = asms[sys] && asms[sys].fingerprint !== undefined ? asms[sys].fingerprint : null;
  });
  return out;
}

function queueOf(kind) {
  return kind === 'pangenome' || kind === 'pan' ? 'pan' : 'spec';
}

function newJobDoc(id, norm, nowMs) {
  return {
    job_id: id,
    kind: norm.kind,
    status: 'queued',
    request: norm.request,
    resolved: norm.resolved,
    created_at: new Date(nowMs).toISOString(),
    started_at: null,
    finished_at: null,
    progress: { done: 0, total: (norm.estimate && norm.estimate.total) || 0, stage: 'queued', running: [] },
    attempts: 0,
    worker: null,
    warnings: norm.warnings || [],
    estimate: { cpu_s: norm.estimate ? norm.estimate.cpu_s : null },
    error: null
  };
}

function publicProgress(job) {
  const p = (job && job.progress) || {};
  return {
    done: Number.isInteger(p.done) ? p.done : 0,
    total: Number.isInteger(p.total) ? p.total : 0,
    stage: typeof p.stage === 'string' ? p.stage : (job && job.status) || 'queued',
    running: Array.isArray(p.running) ? p.running : []
  };
}

function logSafe(log, level, msg) {
  try {
    (log[level] || log.log || function () {}).call(log, msg);
  } catch (e) { /* ignore */ }
}

function storeUnavailable() {
  return new PrimerHttpError(503, 'JOB_STORE_UNAVAILABLE', 'the check job store is temporarily unavailable',
    { retry_after_s: STORE_RETRY_AFTER_S });
}

function unknownJob(id) {
  return new PrimerHttpError(404, 'UNKNOWN_JOB', 'unknown or expired check job',
    typeof id === 'string' && JOB_ID_RE.test(id) ? { job_id: id } : {});
}

function featureDisabled() {
  return new PrimerHttpError(503, 'FEATURE_DISABLED', 'primer checks are disabled on this server',
    { retry_after_s: FEATURE_DISABLED_RETRY_AFTER_S });
}

// Run one store operation under the deadline (default 3 s). Any failure -> 503 JOB_STORE_UNAVAILABLE.
async function storeCall(fn, deps) {
  deps = deps || {};
  const ms = deps.storeDeadlineMs > 0 ? deps.storeDeadlineMs : STORE_DEADLINE_MS;
  const log = deps.log || console;
  let timer;
  const work = Promise.resolve().then(fn);
  work.catch(function () { /* a late failure after the deadline is already reported */ });
  const deadline = new Promise(function (resolve, reject) {
    timer = setTimeout(function () {
      const e = new Error('no answer within ' + ms + ' ms');
      e.deadline = true;
      reject(e);
    }, ms);
  });
  try {
    return await Promise.race([work, deadline]);
  } catch (err) {
    if (isPrimerHttpError(err)) throw err;
    logSafe(log, 'error', 'primers job store unavailable: ' + redactPaths(String(err && err.message)));
    throw storeUnavailable();
  } finally {
    clearTimeout(timer);
  }
}

let storeSingleton = null;
let storeOverride = null;

// The process-wide store chosen by cfg.check.store ('redis' | 'memory'), created lazily.
function getStore() {
  if (storeOverride) return storeOverride;
  if (!storeSingleton) {
    const config = require('../config');
    const cfg = config.get();
    const siteKey = config.siteKey();
    storeSingleton = cfg.check.store === 'memory'
      ? require('./memory_store').createMemoryStore({ cfg: cfg, siteKey: siteKey })
      : require('./redis_store').createRedisStore({ cfg: cfg, siteKey: siteKey, role: 'api' });
  }
  return storeSingleton;
}

async function submit(body, deps) {
  deps = deps || {};
  const cfg = deps.cfg || require('../config').get();
  if (!cfg.enabled) throw featureDisabled();
  const check = deps.check || require('../check');
  const norm = await check.normalize(body, Object.assign({}, deps, { cfg: cfg }));
  const dbs = norm.dbs || dbsOf(norm.resolved);
  const id = jobId(norm.request, dbs, check.ALGORITHM_VERSION);
  const doc = newJobDoc(id, norm, (deps.now || Date.now)());
  const store = deps.store || getStore();

  const res = await storeCall(function () { return store.submit(doc); }, deps);
  if (res.outcome === 'FULL') {
    throw new PrimerHttpError(503, 'QUEUE_FULL', 'the check queue is full; try again later',
      { retry_after_s: QUEUE_FULL_RETRY_AFTER_S, max_queued: cfg.check.max_queued });
  }
  const job = res.outcome === 'EXISTS' && res.job ? res.job : doc;
  let position = null;
  if (job.status === 'queued') {
    const p = await storeCall(function () { return store.queuePosition(id, queueOf(job.kind)); }, deps);
    position = Number.isInteger(p) && p >= 0 ? p : null;
  }
  return {
    job_id: id,
    status: job.status,
    kind: job.kind,
    queue_position: position,
    progress: publicProgress(job),
    estimate: { cpu_s: norm.estimate ? norm.estimate.cpu_s : null },
    created_at: job.created_at,
    warnings: norm.warnings || [],
    created: res.outcome === 'QUEUED'
  };
}

async function status(id, deps) {
  deps = deps || {};
  const cfg = deps.cfg || require('../config').get();
  if (!cfg.enabled) throw featureDisabled();
  if (typeof id !== 'string' || !JOB_ID_RE.test(id)) throw unknownJob(id);
  const store = deps.store || getStore();
  const job = await storeCall(function () { return store.getJob(id); }, deps);
  if (!job || typeof job !== 'object') throw unknownJob(id);
  let results = null;
  let position = null;
  if (job.status === 'done') {
    results = await storeCall(function () { return store.getResults(id); }, deps);
  } else if (job.status === 'running') {
    results = await storeCall(function () { return store.getPartial(id); }, deps);
  } else if (job.status === 'queued') {
    const p = await storeCall(function () { return store.queuePosition(id, queueOf(job.kind)); }, deps);
    position = Number.isInteger(p) && p >= 0 ? p : null;
  }
  return {
    job_id: id,
    status: job.status,
    kind: job.kind,
    partial: job.status !== 'done' && results !== null && results !== undefined,
    queue_position: position,
    progress: publicProgress(job),
    created_at: job.created_at || null,
    started_at: job.started_at || null,
    finished_at: job.finished_at || null,
    attempts: Number.isInteger(job.attempts) ? job.attempts : 0,
    request: job.request || null,
    warnings: Array.isArray(job.warnings) ? job.warnings : [],
    estimate: job.estimate || null,
    results: results === undefined ? null : results,
    error: job.status === 'error' ? (job.error || { code: 'CHECK_FAILED', message: 'the check failed' }) : null
  };
}

// The startup line both the API and the worker print, so their site keys can be compared in the logs.
function siteKeyLogLine(siteKey) {
  return 'primers site_key=' + siteKey;
}

function _setStoreForTests(store) {
  storeOverride = store || null;
}

async function _resetForTests() {
  storeOverride = null;
  const s = storeSingleton;
  storeSingleton = null;
  if (s && typeof s.close === 'function') await s.close();
}

module.exports = {
  JOB_ID_RE,
  STORE_DEADLINE_MS,
  canonicalJSON,
  jobId,
  dbsOf,
  queueOf,
  newJobDoc,
  publicProgress,
  storeCall,
  getStore,
  submit,
  status,
  siteKeyLogLine,
  _setStoreForTests,
  _resetForTests
};
