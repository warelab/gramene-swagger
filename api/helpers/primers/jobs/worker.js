'use strict';

// Check worker (spec §A.8.3 worker loop / runJob / shutdown, §A.8.5 ctx). Runs inside the separate pm2 process
// started by worker_main.js; there is no supervisor.
//
// Loop:
//   1. acquire the per-site worker lock (SET NX PX 30 s); another live holder -> retry every 5 s
//   2. after acquiring: requeueRunning (attempts >= max_attempts -> fail WORKER_LOST), queued-key TTL sweep now
//      and every queue_sweep_ms; the lock is refreshed every 10 s (lost to another worker -> abandon jobs, re-acquire)
//   3. every 1 s: claim (Lua; host-wide/local/pan caps) while fewer than local_max_jobs jobs run here, then runJob
//   any loop error -> log and sleep 5 s
// runJob: CAS queued->running (attempts+1), AbortController job timeout, 10 s heartbeat, check.run(request, ctx),
//   progress/partial throttled to one write per progress_min_interval_ms (+ final), result size cap, complete/fail,
//   tmpdir removal. Shutdown (stop()): requeue running jobs to the FRONT of their queue, abort them, release the lock.

const fs = require('fs');
const os = require('os');
const path = require('path');
const crypto = require('crypto');
const { redactPaths } = require('../errors');

const DEFAULT_TIMING = Object.freeze({
  lockTtlMs: 30000,
  lockRefreshMs: 10000,
  lockRetryMs: 5000,
  claimIntervalMs: 1000,
  errorSleepMs: 5000,
  shutdownGraceMs: 700,
  shutdownWriteMs: 250,
  drainMs: 5000,
  finalizeAttempts: 3,
  killGraceMs: 3000
});
// NCBI_DONT_USE_NCBIRC: BLAST+ would otherwise read .ncbirc from its cwd, $HOME, /etc and its own directory.
const CHILD_ENV = Object.freeze({ PATH: '/usr/bin:/bin', NCBI_DONT_USE_NCBIRC: '1' });
const MAX_RUNNING_LIST = 50;
const MAX_STAGE_LENGTH = 40;
const MAX_ERROR_MESSAGE = 500;
const ERROR_CODE_RE = /^[A-Z][A-Z0-9_]{1,63}$/;
const FATAL_CODES = new Set(['MONGO_UNAVAILABLE']);

function codedError(code, message, extra) {
  const e = new Error(message);
  e.code = code;
  if (extra) Object.assign(e, extra);
  return e;
}

function queueOf(kind) {
  return kind === 'pangenome' || kind === 'pan' ? 'pan' : 'spec';
}

function iso(ms) {
  return new Date(ms).toISOString();
}

function totalOf(doc) {
  const t = doc && doc.progress && doc.progress.total;
  return Number.isInteger(t) && t >= 0 ? t : 0;
}

function formatDuration(ms) {
  if (ms >= 60000) return Math.round(ms / 60000) + ' min';
  if (ms >= 1000) return Math.round(ms / 1000) + ' s';
  return Math.max(0, Math.round(ms)) + ' ms';
}

function errMessage(err) {
  return redactPaths(String((err && err.message) || err));
}

// ---- job doc transforms (pure) -------------------------------------------------------------------

function startedDoc(doc, o) {
  return Object.assign({}, doc, {
    status: 'running',
    started_at: iso(o.now),
    finished_at: null,
    attempts: (Number.isInteger(doc.attempts) ? doc.attempts : 0) + 1,
    worker: o.worker,
    error: null,
    progress: { done: 0, total: totalOf(doc), stage: 'starting', running: [] }
  });
}

function requeuedDoc(doc) {
  return Object.assign({}, doc, {
    status: 'queued',
    started_at: null,
    worker: null,
    progress: { done: 0, total: totalOf(doc), stage: 'queued', running: [] }
  });
}

// Finished docs drop `resolved` (only the worker needs it; it can be ~100 KB with a full pan-genome).
function doneDoc(doc, o) {
  const total = totalOf(doc);
  return Object.assign({}, doc, {
    status: 'done',
    finished_at: iso(o.now),
    resolved: null,
    error: null,
    progress: { done: total, total: total, stage: 'done', running: [] }
  });
}

function errorDoc(doc, error, o) {
  const p = doc.progress || {};
  const e = { code: error.code, message: error.message };
  if (error.details && typeof error.details === 'object') e.details = error.details;
  return Object.assign({}, doc, {
    status: 'error',
    finished_at: iso(o.now),
    resolved: null,
    error: e,
    progress: { done: Number.isInteger(p.done) ? p.done : 0, total: totalOf(doc), stage: 'error', running: [] }
  });
}

// {code, message} for a failed run. A timeout always reports JOB_TIMEOUT; Node system errors
// (ENOENT & co.) and codes that are not UPPER_SNAKE become CHECK_FAILED.
function jobError(err, o) {
  o = o || {};
  if (o.timedOut) return { code: 'JOB_TIMEOUT', message: 'check exceeded ' + formatDuration(o.timeoutMs) };
  const code = err && typeof err.code === 'string' && ERROR_CODE_RE.test(err.code) &&
    !err.syscall && typeof err.errno !== 'number' ? err.code : 'CHECK_FAILED';
  return { code: code, message: errMessage(err || 'check failed').slice(0, MAX_ERROR_MESSAGE) || 'check failed' };
}

function isFatal(err) {
  return !!err && (err.fatal === true || FATAL_CODES.has(err.code));
}

function sanitizeProgress(p, prev) {
  p = p || {};
  prev = prev || {};
  const int = function (v, d) { return Number.isInteger(v) && v >= 0 ? v : d; };
  return {
    done: int(p.done, int(prev.done, 0)),
    total: int(p.total, int(prev.total, 0)),
    stage: typeof p.stage === 'string' && p.stage ? p.stage.slice(0, MAX_STAGE_LENGTH) : (prev.stage || 'running'),
    running: Array.isArray(p.running)
      ? p.running.filter(function (s) { return typeof s === 'string'; }).slice(0, MAX_RUNNING_LIST)
      : []
  };
}

function prefixedLog(base, prefix) {
  function out(level) {
    return function () {
      const args = Array.prototype.slice.call(arguments);
      try {
        (base[level] || base.log || function () {}).apply(base, [prefix].concat(args));
      } catch (e) { /* logging must never throw */ }
    };
  }
  return { info: out('info'), warn: out('warn'), error: out('error') };
}

function sleepPlain(ms) {
  return new Promise(function (resolve) { setTimeout(resolve, ms); });
}

function bounded(promise, ms) {
  let timer;
  return Promise.race([
    Promise.resolve(promise).catch(function () {}),
    new Promise(function (resolve) { timer = setTimeout(resolve, ms); })
  ]).finally(function () { clearTimeout(timer); });
}

// ---- worker ----------------------------------------------------------------------------------------

// createWorker({store, cfg, siteKey, check, log, mongo, spawnLines, now, token, timing, onFatal})
//   store     memory_store / redis_store (role 'worker')
//   check     {run(request, ctx)} (default require('../check'))
//   spawnLines(cmd, args, opts) implementation behind ctx.spawnLines (default check/blast.js spawnLinesLocal)
//   onFatal(err) called after a job failed with a fatal error (MONGO_UNAVAILABLE) was requeued
function createWorker(opts) {
  opts = opts || {};
  const cfg = opts.cfg || require('../config').get();
  const ccfg = cfg.check;
  const store = opts.store;
  if (!store) throw new TypeError('createWorker: store is required');
  const siteKey = String(opts.siteKey || '');
  if (!siteKey) throw new TypeError('createWorker: siteKey is required');
  const check = opts.check || require('../check');
  const now = opts.now || Date.now;
  const timing = Object.assign({}, DEFAULT_TIMING, opts.timing || {});
  const token = opts.token || (os.hostname() + ':' + process.pid + ':' + crypto.randomBytes(4).toString('hex'));
  const spawnImpl = opts.spawnLines || function (cmd, args, o) { return require('../check/blast').spawnLinesLocal(cmd, args, o); };
  const onFatal = typeof opts.onFatal === 'function' ? opts.onFatal : null;
  const log = prefixedLog(opts.log || console, 'primers worker [' + siteKey + ']');
  const localMax = Math.max(1, Number(ccfg.local_max_jobs) || 1);
  const minWriteMs = Math.max(0, Number(ccfg.progress_min_interval_ms) || 0);
  const allowed = new Set([cfg.blastn, cfg.blastdbcmd].filter(function (s) { return typeof s === 'string' && s !== ''; }));

  const active = new Map();
  const stats = { claimed: 0, done: 0, failed: 0, requeued: 0, skipped: 0, abandoned: 0, lock_acquisitions: 0 };
  let stopping = false;
  let haveLock = false;
  let loopPromise = null;
  let stopPromise = null;
  let lockTimer = null;
  let sweepTimer = null;
  let wakeFn = null;

  function sleep(ms) {
    return new Promise(function (resolve) {
      let timer = null;
      function done() {
        clearTimeout(timer);
        if (wakeFn === done) wakeFn = null;
        resolve();
      }
      timer = setTimeout(done, ms);
      wakeFn = done;
    });
  }

  function wake() {
    if (wakeFn) wakeFn();
  }

  function clearLoopTimers() {
    if (lockTimer) clearInterval(lockTimer);
    if (sweepTimer) clearInterval(sweepTimer);
    lockTimer = null;
    sweepTimer = null;
  }

  // ---- lock ----

  async function refreshLock() {
    if (!haveLock || stopping) return;
    let r;
    try {
      r = await store.refreshWorkerLock(token, timing.lockTtlMs);
    } catch (err) {
      log.warn('worker lock refresh failed: ' + errMessage(err));
      return;
    }
    if (r === 2) log.warn('worker lock had expired; re-acquired');
    if (r !== 0 || !haveLock) return;
    haveLock = false;
    clearLoopTimers();
    log.error('worker lock taken over by another worker; abandoning ' + active.size + ' running job(s)');
    active.forEach(function (h) { abandon(h); });
    wake();
  }

  async function sweep() {
    const r = await store.sweepQueued();
    if (r && r.removed) log.info('queue sweep: refreshed ' + r.refreshed + ' queued job(s), dropped ' + r.removed + ' stale id(s)');
    return r;
  }

  // Requeue (or fail after max_attempts) every job left in the running set by a dead worker.
  // {staleOnly: true} is the periodic reconciliation while this worker holds the lock: only entries this
  // worker is not running AND whose heartbeat is older than stale_ms are touched (a claim whose start or
  // final write failed on a Redis error), so a job that was just claimed is never taken for an orphan.
  async function requeueRunning(o) {
    o = o || {};
    const out = { requeued: [], failed: [], cleaned: [] };
    const cutoff = now() - Math.max(0, Number(ccfg.stale_ms) || 60000);
    const entries = await store.runningEntries();
    const maxAttempts = Math.max(1, Number(ccfg.max_attempts) || 1);
    for (const entry of entries) {
      const id = entry[0];
      if (active.has(id)) continue;
      if (o.staleOnly && entry[1] > cutoff) continue;
      const doc = await store.getJob(id);
      if (active.has(id)) continue;
      if (!doc || (doc.status !== 'running' && doc.status !== 'queued')) {
        await store.release(id);
        out.cleaned.push(id);
        continue;
      }
      const attempts = Number.isInteger(doc.attempts) ? doc.attempts : 0;
      if (doc.status === 'running' && attempts >= maxAttempts) {
        const failed = errorDoc(doc, {
          code: 'WORKER_LOST',
          message: 'the worker stopped during the check ' + attempts + ' time(s); giving up'
        }, { now: now() });
        if (await store.fail(failed, { expect: ['running', 'queued'], worker: '', now: now() })) {
          out.failed.push(id);
          stats.failed++;
        }
      } else if (await store.requeue(requeuedDoc(doc), { expect: ['running', 'queued'], worker: '', front: true })) {
        out.requeued.push(id);
        stats.requeued++;
      }
    }
    if (out.requeued.length || out.failed.length || out.cleaned.length) {
      log.warn('recovered running jobs: requeued ' + out.requeued.length + ', failed ' + out.failed.length +
        ' (WORKER_LOST), cleaned ' + out.cleaned.length);
    }
    return out;
  }

  // ---- loop ----

  async function loop() {
    while (!stopping) {
      try {
        if (!haveLock) {
          if (active.size > 0) {
            await sleep(timing.claimIntervalMs);
            continue;
          }
          if (!(await store.acquireWorkerLock(token, timing.lockTtlMs))) {
            await sleep(timing.lockRetryMs);
            continue;
          }
          if (stopping) break;
          haveLock = true;
          stats.lock_acquisitions++;
          log.info('acquired the worker lock (' + token + ')');
          clearLoopTimers();
          lockTimer = setInterval(function () { refreshLock(); }, timing.lockRefreshMs);
          try {
            await requeueRunning();
            await sweep();
          } catch (err) {
            // Keep the lock but retry recovery: a failure here must not leave jobs stranded.
            haveLock = false;
            clearLoopTimers();
            try { await store.releaseWorkerLock(token); } catch (e) { /* ignore */ }
            throw err;
          }
          sweepTimer = setInterval(function () {
            sweep()
              .then(function () { return haveLock && !stopping ? requeueRunning({ staleOnly: true }) : null; })
              .catch(function (err) { log.warn('queue sweep failed: ' + errMessage(err)); });
          }, Math.max(1000, Number(ccfg.queue_sweep_ms) || 600000));
        }
        while (!stopping && haveLock && active.size < localMax) {
          const claimed = await store.claim(now());
          if (!claimed) break;
          if (stopping || !haveLock) {
            // Put it back where it was: the claim popped it off the front of its queue. REQUEUE also drops the
            // running entry and slots atomically; release only when there is nothing to requeue. If getJob or
            // requeue throws, the id stays in the running set and requeueRunning recovers it.
            const doc = await store.getJob(claimed.id);
            const requeued = !!doc && doc.status === 'queued' &&
              await store.requeue(doc, { expect: ['queued'], worker: '', front: true });
            if (!requeued) await store.release(claimed.id);
            break;
          }
          stats.claimed++;
          launch(claimed);
        }
        await sleep(timing.claimIntervalMs);
      } catch (err) {
        log.error('loop error: ' + errMessage(err));
        await sleep(timing.errorSleepMs);
      }
    }
  }

  function start() {
    if (!loopPromise) {
      loopPromise = loop().catch(function (err) { log.error('worker loop crashed: ' + errMessage(err)); });
    }
    return loopPromise;
  }

  // ---- jobs ----

  function launch(claimed) {
    const h = {
      id: claimed.id,
      queue: claimed.queue,
      ac: new AbortController(),
      doc: null,
      progress: null,
      partialJson: null,
      tmp: null,
      timedOut: false,
      shutdown: false,
      abandoned: false,
      finalized: false,
      writes: Promise.resolve(),
      throttle: { progress: { last: -Infinity, timer: null, dirty: false }, partial: { last: -Infinity, timer: null, dirty: false } },
      timeoutTimer: null,
      heartbeatTimer: null,
      done: null,
      partialTooLargeLogged: false
    };
    active.set(h.id, h);
    h.done = runJob(h)
      .catch(function (err) { log.error('job ' + h.id + ' crashed: ' + errMessage(err)); })
      .finally(function () {
        clearJobTimers(h);
        removeTmp(h);
        if (active.get(h.id) === h) active.delete(h.id);
        wake();
      });
    return h;
  }

  function clearJobTimers(h) {
    if (h.timeoutTimer) clearTimeout(h.timeoutTimer);
    if (h.heartbeatTimer) clearInterval(h.heartbeatTimer);
    h.timeoutTimer = null;
    h.heartbeatTimer = null;
    ['progress', 'partial'].forEach(function (k) {
      if (h.throttle[k].timer) clearTimeout(h.throttle[k].timer);
      h.throttle[k].timer = null;
    });
  }

  function removeTmp(h) {
    if (!h.tmp) return;
    const dir = h.tmp;
    h.tmp = null;
    try {
      fs.rmSync(dir, { recursive: true, force: true });
    } catch (err) {
      log.warn('could not remove the job tmpdir: ' + errMessage(err));
    }
  }

  function abandon(h) {
    if (h.abandoned) return;
    h.abandoned = true;
    stats.abandoned++;
    clearJobTimers(h);
    h.ac.abort(codedError('LOCK_LOST', 'the worker lock was lost'));
  }

  function enqueueWrite(h, label, fn) {
    const p = h.writes.then(function () {
      if (h.abandoned) return undefined;
      return fn();
    }).catch(function (err) {
      log.warn(label + ' write for job ' + h.id + ' failed: ' + errMessage(err));
    });
    h.writes = p;
    return p;
  }

  function flushThrottled(h, kind) {
    const st = h.throttle[kind];
    if (st.timer) clearTimeout(st.timer);
    st.timer = null;
    if (!st.dirty || h.finalized || h.abandoned) return;
    st.dirty = false;
    st.last = Date.now();
    if (kind === 'progress') {
      enqueueWrite(h, 'progress', function () {
        if (h.finalized) return false;
        h.doc = Object.assign({}, h.doc, { progress: h.progress });
        return store.setJob(h.doc, { expect: ['running'], worker: token });
      });
    } else {
      enqueueWrite(h, 'partial', function () {
        if (h.finalized || h.partialJson === null) return false;
        return store.setPartial(h.id, h.partialJson, { worker: token });
      });
    }
  }

  function scheduleWrite(h, kind) {
    const st = h.throttle[kind];
    st.dirty = true;
    if (st.timer) return;
    const wait = st.last + minWriteMs - Date.now();
    if (wait <= 0) flushThrottled(h, kind);
    else st.timer = setTimeout(function () { st.timer = null; flushThrottled(h, kind); }, wait);
  }

  // The job's private working directory (mkdtemp: mode 0700), created on first use and removed with the job.
  function jobTmpdir(h) {
    if (!h.tmp) {
      const base = typeof cfg.tmp_dir === 'string' && path.isAbsolute(cfg.tmp_dir) ? cfg.tmp_dir : os.tmpdir();
      h.tmp = fs.mkdtempSync(path.join(base, 'primers-check-' + h.id.slice(0, 8) + '-'));
    }
    return h.tmp;
  }

  function buildCtx(h) {
    const jobLog = prefixedLog(opts.log || console, 'primers worker [' + siteKey + '] job ' + h.id);
    return {
      jobId: h.id,
      siteKey: siteKey,
      kind: h.doc.kind,
      attempt: h.doc.attempts,
      signal: h.ac.signal,
      config: cfg,
      request: h.doc.request,
      resolved: h.doc.resolved,
      procs: queueOf(h.doc.kind) === 'pan' ? ccfg.pangenome_job_procs : ccfg.spec_job_procs,
      progress: function (p) {
        if (h.finalized || h.abandoned) return;
        h.progress = sanitizeProgress(p, h.progress || h.doc.progress);
        scheduleWrite(h, 'progress');
      },
      partial: function (results) {
        if (h.finalized || h.abandoned) return;
        // Snapshot now: run() hands over its live results document and keeps mutating it.
        let json;
        try {
          json = JSON.stringify(results === undefined ? null : results);
        } catch (err) {
          log.warn('partial results of job ' + h.id + ' are not serializable: ' + errMessage(err));
          return;
        }
        if (Buffer.byteLength(json) > ccfg.max_result_bytes) {
          if (!h.partialTooLargeLogged) log.warn('partial results of job ' + h.id + ' exceed max_result_bytes; not stored');
          h.partialTooLargeLogged = true;
          return;
        }
        h.partialJson = json;
        scheduleWrite(h, 'partial');
      },
      tmpdir: function () { return jobTmpdir(h); },
      // Mongo loss during the run must throw MONGO_UNAVAILABLE {fatal: true} (requeue + exit 75 via onFatal)
      // instead of finishing with ANNOTATION_UNAVAILABLE: the config lib never reconnects.
      fatalOnMongoUnavailable: true,
      log: jobLog,
      // Allow-list [cfg.blastn, cfg.blastdbcmd]; niced; SIGTERM then SIGKILL after 3 s on job abort.
      // cwd is always the job's private tmpdir, never a shared directory such as /tmp (BLAST+ reads ./.ncbirc).
      spawnLines: function (cmd, args, o) {
        o = o || {};
        if (typeof cmd !== 'string' || !allowed.has(cmd)) {
          return Promise.reject(codedError('SPAWN_NOT_ALLOWED', 'command not allowed: ' + path.basename(String(cmd))));
        }
        if (!Array.isArray(args) || args.some(function (a) { return typeof a !== 'string'; })) {
          return Promise.reject(new TypeError('spawnLines: args must be an array of strings'));
        }
        let cwd;
        try {
          cwd = jobTmpdir(h);
        } catch (err) {
          return Promise.reject(codedError('TMPDIR_FAILED', 'could not create the job working directory: ' + errMessage(err)));
        }
        const signal = o.signal ? AbortSignal.any([h.ac.signal, o.signal]) : h.ac.signal;
        return spawnImpl(cmd, args, {
          stdin: o.stdin,
          timeoutMs: o.timeoutMs,
          onLine: o.onLine,
          signal: signal,
          nice: Number.isInteger(ccfg.nice) ? ccfg.nice : 10,
          env: CHILD_ENV,
          cwd: cwd,
          killGraceMs: timing.killGraceMs
        });
      },
      get mongo() {
        return opts.mongo || require('gramene-mongodb-config');
      }
    };
  }

  async function retryStore(label, fn) {
    let lastErr = null;
    for (let i = 0; i < timing.finalizeAttempts; i++) {
      try {
        return await fn();
      } catch (err) {
        lastErr = err;
        log.warn(label + ' failed (attempt ' + (i + 1) + '): ' + errMessage(err));
        if (i + 1 < timing.finalizeAttempts) await sleepPlain(timing.errorSleepMs);
      }
    }
    log.error(label + ' gave up: ' + errMessage(lastErr));
    return false;
  }

  async function requeueHandle(h, reason) {
    const ok = await retryStore('requeue of job ' + h.id, function () {
      return store.requeue(requeuedDoc(h.doc), { expect: ['running'], worker: token, front: true });
    });
    if (ok) {
      stats.requeued++;
      log.info('job ' + h.id + ' requeued (' + reason + ')');
    }
    return ok;
  }

  async function runJob(h) {
    // 1. queued -> running (compare-and-set), attempts + 1.
    const cur = await store.getJob(h.id);
    if (!cur || cur.status !== 'queued') {
      stats.skipped++;
      await store.release(h.id);
      return 'skipped';
    }
    const started = startedDoc(cur, { worker: token, now: now() });
    if (!(await store.setJob(started, { expect: ['queued'], worker: '' }))) {
      stats.skipped++;
      await store.release(h.id);
      return 'skipped';
    }
    h.doc = started;
    h.progress = started.progress;
    if (h.abandoned) return 'abandoned';
    if (h.shutdown || stopping) {
      h.finalized = true;
      await requeueHandle(h, 'shutdown');
      return 'requeued';
    }
    log.info('job ' + h.id + ' started (' + started.kind + ', attempt ' + started.attempts + ')');

    // 2. timeout + heartbeat.
    const timeoutMs = Math.max(1, Number(ccfg.job_timeout_ms) || 1800000);
    h.timeoutTimer = setTimeout(function () {
      h.timedOut = true;
      h.ac.abort(codedError('JOB_TIMEOUT', 'check exceeded ' + formatDuration(timeoutMs)));
    }, timeoutMs);
    h.heartbeatTimer = setInterval(function () {
      enqueueWrite(h, 'heartbeat', function () {
        return store.heartbeat(h.id, { queue: h.queue, worker: token, now: now() });
      });
    }, Math.max(1, Number(ccfg.heartbeat_ms) || 10000));

    // 3. run; an abort (timeout, shutdown, lock loss) ends the wait even if run() ignores the signal.
    let results;
    let error = null;
    const aborted = new Promise(function (resolve) {
      if (h.ac.signal.aborted) resolve();
      else h.ac.signal.addEventListener('abort', resolve, { once: true });
    });
    try {
      const runP = Promise.resolve().then(function () { return check.run(h.doc.request, buildCtx(h)); });
      runP.catch(function () { /* handled below or ignored after an abort */ });
      const outcome = await Promise.race([
        runP.then(function (r) { return { r: r }; }, function (e) { return { e: e || new Error('check failed') }; }),
        aborted.then(function () { return { aborted: true }; })
      ]);
      if (outcome.aborted) error = h.ac.signal.reason || codedError('ABORTED', 'aborted');
      else if (outcome.e) error = outcome.e;
      else if (h.ac.signal.aborted) error = h.ac.signal.reason || codedError('ABORTED', 'aborted');
      else results = outcome.r;
    } catch (e) {
      error = e;
    }
    clearJobTimers(h);

    if (h.abandoned) return 'abandoned';
    if (h.finalized) return 'shutdown';
    await bounded(h.writes, timing.drainMs);
    if (h.abandoned) return 'abandoned';
    if (h.finalized) return 'shutdown';
    h.finalized = true;

    // 4. finalize.
    if (error && (h.shutdown || stopping) && !h.timedOut) {
      await requeueHandle(h, 'shutdown');
      return 'requeued';
    }
    if (error && isFatal(error)) {
      const maxAttempts = Math.max(1, Number(ccfg.max_attempts) || 1);
      let outcome = 'requeued';
      if ((Number.isInteger(h.doc.attempts) ? h.doc.attempts : 0) >= maxAttempts) {
        // Fatal on every attempt: fail it (error TTL, a re-POST re-runs it) so one job cannot block the queue
        // through endless exit-75 restarts.
        const je = jobError(error, { timedOut: false });
        await retryStore('fail of job ' + h.id, function () {
          return store.fail(errorDoc(h.doc, je, { now: now() }), { expect: ['running'], worker: token, now: now() });
        });
        stats.failed++;
        outcome = 'error';
      } else {
        await requeueHandle(h, error.code || 'fatal');
      }
      log.error('job ' + h.id + ': fatal error ' + (error.code || '') + ' (' + outcome + '): ' + errMessage(error));
      if (onFatal) {
        try { onFatal(error); } catch (e) { /* ignore */ }
      }
      return outcome;
    }
    if (error) {
      const je = jobError(error, { timedOut: h.timedOut, timeoutMs: timeoutMs });
      if (je.code === 'CHECK_FAILED' && error && error.stack) log.error('job ' + h.id + ' failed: ' + redactPaths(error.stack));
      await retryStore('fail of job ' + h.id, function () {
        return store.fail(errorDoc(h.doc, je, { now: now() }), { expect: ['running'], worker: token, now: now() });
      });
      stats.failed++;
      log.warn('job ' + h.id + ' error ' + je.code + ': ' + je.message);
      return 'error';
    }
    let json;
    try {
      json = JSON.stringify(results === undefined ? null : results);
    } catch (e) {
      await retryStore('fail of job ' + h.id, function () {
        return store.fail(errorDoc(h.doc, { code: 'CHECK_FAILED', message: 'the check results could not be serialized' },
          { now: now() }), { expect: ['running'], worker: token, now: now() });
      });
      stats.failed++;
      return 'error';
    }
    const bytes = Buffer.byteLength(json);
    const limit = Number(ccfg.max_result_bytes) || 5000000;
    if (bytes > limit) {
      await retryStore('fail of job ' + h.id, function () {
        return store.fail(errorDoc(h.doc, {
          code: 'RESULT_TOO_LARGE',
          message: 'the check results are ' + bytes + ' bytes, over the limit of ' + limit,
          details: { bytes: bytes, limit: limit }
        }, { now: now() }), { expect: ['running'], worker: token, now: now() });
      });
      stats.failed++;
      log.warn('job ' + h.id + ' error RESULT_TOO_LARGE (' + bytes + ' bytes)');
      return 'error';
    }
    await retryStore('complete of job ' + h.id, function () {
      return store.complete(doneDoc(h.doc, { now: now() }), json, { expect: ['running'], worker: token, now: now() });
    });
    stats.done++;
    log.info('job ' + h.id + ' done (' + bytes + ' bytes)');
    return 'done';
  }

  // ---- shutdown ----

  async function shutdownRequeue(h) {
    h.shutdown = true;
    if (h.finalized || h.abandoned || !h.doc) return;
    h.finalized = true;
    clearJobTimers(h);
    await bounded(h.writes, timing.shutdownWriteMs);
    try {
      if (await store.requeue(requeuedDoc(h.doc), { expect: ['running'], worker: token, front: true })) {
        stats.requeued++;
        log.info('job ' + h.id + ' requeued (shutdown)');
      }
    } catch (err) {
      log.error('could not requeue job ' + h.id + ' on shutdown: ' + errMessage(err) +
        '; the next worker will recover it');
    }
  }

  // Stop claiming, requeue running jobs to the front of their queues, abort them and release the lock.
  // Resolves within roughly shutdownGraceMs + shutdownWriteMs even if run() ignores the abort.
  function stop() {
    if (stopPromise) return stopPromise;
    stopping = true;
    wake();
    stopPromise = (async function () {
      clearLoopTimers();
      const handles = Array.from(active.values());
      await bounded(Promise.all(handles.map(shutdownRequeue)), timing.shutdownGraceMs);
      handles.forEach(function (h) { h.ac.abort(codedError('WORKER_SHUTDOWN', 'the worker is shutting down')); });
      await bounded(Promise.all(handles.map(function (h) { return h.done; })), timing.shutdownWriteMs);
      if (loopPromise) await bounded(loopPromise, timing.shutdownWriteMs);
      if (haveLock) {
        haveLock = false;
        await bounded(store.releaseWorkerLock(token), timing.shutdownWriteMs);
      }
    })();
    return stopPromise;
  }

  return {
    token: token,
    siteKey: siteKey,
    stats: stats,
    start: start,
    stop: stop,
    requeueRunning: requeueRunning,
    sweep: sweep,
    refreshLock: refreshLock,
    activeIds: function () { return Array.from(active.keys()); },
    get haveLock() { return haveLock; },
    get stopping() { return stopping; }
  };
}

module.exports = {
  DEFAULT_TIMING,
  CHILD_ENV,
  createWorker,
  startedDoc,
  requeuedDoc,
  doneDoc,
  errorDoc,
  jobError,
  isFatal,
  sanitizeProgress,
  queueOf
};
