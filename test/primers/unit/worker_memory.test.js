'use strict';

// jobs/worker.js + worker_main.js with the memory store and stub check.run():
// done, error, timeout, shutdown requeue, result cap, throttled progress/partial, spawnLines, locks, recovery.

const test = require('node:test');
const should = require('should');
const fs = require('fs');
const os = require('os');
const path = require('path');
const crypto = require('crypto');
const { spawn } = require('child_process');

const jobs = require('../../../api/helpers/primers/jobs');
const { createMemoryStore } = require('../../../api/helpers/primers/jobs/memory_store');
const workerLib = require('../../../api/helpers/primers/jobs/worker');
const workerMain = require('../../../api/helpers/primers/jobs/worker_main');
const config = require('../../../api/helpers/primers/config');

const WORKER_MAIN = path.resolve(__dirname, '../../../api/helpers/primers/jobs/worker_main.js');
const quiet = { info() {}, warn() {}, error() {}, log() {} };
const FAST = {
  claimIntervalMs: 10, lockRefreshMs: 40, lockRetryMs: 20, errorSleepMs: 20,
  shutdownGraceMs: 300, shutdownWriteMs: 100, drainMs: 500, killGraceMs: 500, finalizeAttempts: 2
};

function makeCfg(check, top) {
  return config._build({
    env: {},
    fileConfig: {},
    overrides: Object.assign({
      tmp_dir: os.tmpdir(),
      blastn: process.execPath,
      blastdbcmd: '/nonexistent/blastdbcmd',
      check: Object.assign({ global_prefix: 'primers:test_worker:', heartbeat_ms: 30, job_timeout_ms: 5000, progress_min_interval_ms: 50 }, check || {})
    }, top || {})
  }).config;
}

function hexId(n) {
  return crypto.createHash('md5').update('worker-job-' + n).digest('hex');
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
    resolved: { assemblies: { sorghum_bicolor: { system_name: 'sorghum_bicolor', fingerprint: 'a'.repeat(40) } }, gene: null, species: null },
    kind: kind,
    warnings: [],
    estimate: { cpu_s: 8, total: kind === 'pangenome' ? 2 : 1 },
    dbs: { sorghum_bicolor: 'a'.repeat(40) }
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
    await sleep(10);
  }
}

function blockUntilAbort(ctx) {
  return new Promise(function (resolve, reject) {
    ctx.signal.addEventListener('abort', function () { reject(ctx.signal.reason); }, { once: true });
  });
}

function setup(o) {
  o = o || {};
  const cfg = makeCfg(o.check, o.top);
  const shared = { slots: new Map(), panSlots: new Map() };
  const store = createMemoryStore({ cfg: cfg, siteKey: 'wtest', shared: shared });
  const check = {
    ALGORITHM_VERSION: '1',
    normalize: async function () { throw new Error('not used'); },
    run: function (request, ctx) { return o.run(request, ctx); }
  };
  const workers = [];
  return {
    cfg: cfg,
    store: store,
    check: check,
    worker: function (extra) {
      const w = workerLib.createWorker(Object.assign({ store: store, cfg: cfg, siteKey: 'wtest', check: check, log: quiet, timing: FAST }, extra || {}));
      workers.push(w);
      return w;
    },
    add: async function (n, kind, pairId) {
      const id = hexId(n);
      const r = await store.submit(jobs.newJobDoc(id, fakeNorm({ kind: kind, pairId: pairId || 'P' + n }), Date.now()));
      r.outcome.should.equal('QUEUED');
      return id;
    },
    status: function (id) {
      return jobs.status(id, { cfg: cfg, store: store, log: quiet });
    },
    stopAll: function () {
      return Promise.all(workers.map(function (w) { return w.stop(); }));
    }
  };
}

async function statusIs(env, id, wanted) {
  const s = await env.status(id);
  return s.status === wanted ? s : null;
}

test('done: results stored, ctx fields per A.8.5, tmpdir removed, slots released', async function () {
  const seen = [];
  const env = setup({
    run: async function (request, ctx) {
      const tmp = ctx.tmpdir();
      fs.writeFileSync(path.join(tmp, 'queries.fa'), '>q0\nACGT\n');
      seen.push({
        jobId: ctx.jobId, siteKey: ctx.siteKey, kind: ctx.kind, procs: ctx.procs, attempt: ctx.attempt,
        resolved: ctx.resolved, request: ctx.request, config: ctx.config, signal: ctx.signal, tmp: tmp,
        tmpAgain: ctx.tmpdir(), log: ctx.log, spawnLines: ctx.spawnLines
      });
      ctx.progress({ done: 1, total: 1, stage: 'reference', running: [] });
      return { ok: true, pairs: request.pairs.map(function (p) { return p.id; }) };
    }
  });
  const specId = await env.add(1, 'specificity', 'S1');
  const w = env.worker();
  w.start();
  try {
    const s = await waitFor(function () { return statusIs(env, specId, 'done'); }, 3000, 'spec job done');
    s.results.should.eql({ ok: true, pairs: ['S1'] });
    s.attempts.should.equal(1);
    s.partial.should.be.false();
    should(s.error).be.null();
    s.progress.should.eql({ done: 1, total: 1, stage: 'done', running: [] });
    const c = seen[0];
    c.jobId.should.equal(specId);
    c.siteKey.should.equal('wtest');
    c.kind.should.equal('specificity');
    c.procs.should.equal(4);
    c.attempt.should.equal(1);
    c.resolved.should.eql(fakeNorm().resolved);
    c.request.should.eql(fakeNorm({ pairId: 'S1' }).request);
    c.config.should.equal(env.cfg);
    c.signal.should.be.instanceOf(AbortSignal);
    c.tmpAgain.should.equal(c.tmp);
    c.tmp.indexOf(os.tmpdir()).should.equal(0);
    c.log.should.have.properties(['info', 'warn', 'error']);
    await waitFor(function () { return !fs.existsSync(c.tmp); }, 2000, 'tmpdir removal');

    const panId = await env.add(2, 'pangenome', 'N1');
    await waitFor(function () { return statusIs(env, panId, 'done'); }, 3000, 'pan job done');
    seen[1].procs.should.equal(8);
    seen[1].kind.should.equal('pangenome');
    await waitFor(function () { return w.activeIds().length === 0; }, 1000, 'idle');
    const st = env.store._state();
    st.running.should.eql([]);
    st.slots.should.eql([]);
    st.panSlots.should.eql([]);
    w.stats.done.should.equal(2);
    w.haveLock.should.be.true();
  } finally {
    await env.stopAll();
  }
  should(env.store._state().lock).be.null();
});

test('error: coded errors keep their code, other errors become CHECK_FAILED; paths are redacted', async function () {
  const env = setup({
    check: { local_max_jobs: 3, global_max_jobs: 3 },
    run: async function (request) {
      const tag = request.pairs[0].id;
      if (tag === 'coded') {
        const e = new Error('blastn failed on /scratch/olson/fasta/sorghum_bicolor/Sb.dna.toplevel');
        e.code = 'REFERENCE_BLAST_FAILED';
        throw e;
      }
      if (tag === 'plain') throw new Error('boom at /home/olson/src/secret.js');
      const e = new Error('ENOENT: no such file or directory');
      e.code = 'ENOENT';
      e.errno = -2;
      e.syscall = 'open';
      throw e;
    }
  });
  const coded = await env.add(1, 'specificity', 'coded');
  const plain = await env.add(2, 'specificity', 'plain');
  const errno = await env.add(3, 'pangenome', 'errno');
  const w = env.worker();
  w.start();
  try {
    const s1 = await waitFor(function () { return statusIs(env, coded, 'error'); }, 3000, 'coded error');
    s1.error.should.eql({ code: 'REFERENCE_BLAST_FAILED', message: 'blastn failed on Sb.dna.toplevel' });
    should(s1.results).be.null();
    s1.partial.should.be.false();
    const s2 = await waitFor(function () { return statusIs(env, plain, 'error'); }, 3000, 'plain error');
    s2.error.code.should.equal('CHECK_FAILED');
    s2.error.message.should.equal('boom at secret.js');
    const s3 = await waitFor(function () { return statusIs(env, errno, 'error'); }, 3000, 'errno error');
    s3.error.code.should.equal('CHECK_FAILED');
    await waitFor(function () { return w.activeIds().length === 0; }, 1000, 'idle');
    env.store._state().running.should.eql([]);
    env.store._state().slots.should.eql([]);
    env.store._state().finished.length.should.equal(3);
    w.stats.failed.should.equal(3);
  } finally {
    await env.stopAll();
  }
});

test('timeout: JOB_TIMEOUT whether run() honours the abort signal or ignores it', async function () {
  const aborted = [];
  const env = setup({
    check: { job_timeout_ms: 150 },
    run: function (request, ctx) {
      ctx.signal.addEventListener('abort', function () { aborted.push(request.pairs[0].id); });
      if (request.pairs[0].id === 'ignores') return new Promise(function () {});
      return blockUntilAbort(ctx);
    }
  });
  const honours = await env.add(1, 'specificity', 'honours');
  const ignores = await env.add(2, 'specificity', 'ignores');
  const w = env.worker();
  w.start();
  try {
    const a = await waitFor(function () { return statusIs(env, honours, 'error'); }, 3000, 'timeout 1');
    a.error.should.eql({ code: 'JOB_TIMEOUT', message: 'check exceeded 150 ms' });
    const b = await waitFor(function () { return statusIs(env, ignores, 'error'); }, 3000, 'timeout 2');
    b.error.code.should.equal('JOB_TIMEOUT');
    aborted.sort().should.eql(['honours', 'ignores']);
    await waitFor(function () { return w.activeIds().length === 0; }, 1000, 'slots freed');
    env.store._state().slots.should.eql([]);
    workerLib.jobError(null, { timedOut: true, timeoutMs: 1800000 }).should.eql({ code: 'JOB_TIMEOUT', message: 'check exceeded 30 min' });
  } finally {
    await env.stopAll();
  }
});

test('shutdown: the running job goes back to the FRONT of its queue within 1 s; a new worker finishes it (attempts 2)', async function () {
  let runs = 0;
  const env = setup({
    check: { local_max_jobs: 1 },
    run: function (request, ctx) {
      runs++;
      if (runs === 1) return blockUntilAbort(ctx);
      return Promise.resolve({ run: runs, attempt: ctx.attempt });
    }
  });
  const first = await env.add(1);
  const w1 = env.worker();
  w1.start();
  await waitFor(function () { return statusIs(env, first, 'running'); }, 3000, 'running');
  const second = await env.add(2);
  const t0 = Date.now();
  await w1.stop();
  const elapsed = Date.now() - t0;
  elapsed.should.be.below(1000);
  const s = await env.status(first);
  s.status.should.equal('queued');
  s.attempts.should.equal(1);
  s.queue_position.should.equal(0);
  should(s.started_at).be.null();
  const st = env.store._state();
  st.queues.spec.should.eql([first, second]);
  st.running.should.eql([]);
  st.slots.should.eql([]);
  should(st.lock).be.null();
  w1.stats.requeued.should.equal(1);

  const w2 = env.worker();
  w2.start();
  try {
    const done = await waitFor(function () { return statusIs(env, first, 'done'); }, 3000, 'finished by the second worker');
    done.attempts.should.equal(2);
    done.results.should.eql({ run: 2, attempt: 2 });
    await waitFor(function () { return statusIs(env, second, 'done'); }, 3000, 'second job');
    runs.should.equal(3);
  } finally {
    await env.stopAll();
  }
});

test('RESULT_TOO_LARGE when the serialized results exceed max_result_bytes', async function () {
  const env = setup({ check: { max_result_bytes: 2000 }, run: async function () { return { blob: 'x'.repeat(5000) }; } });
  const id = await env.add(1);
  const w = env.worker();
  w.start();
  try {
    const s = await waitFor(function () { return statusIs(env, id, 'error'); }, 3000, 'too large');
    s.error.code.should.equal('RESULT_TOO_LARGE');
    s.error.details.limit.should.equal(2000);
    s.error.details.bytes.should.be.above(5000);
    should(await env.store.getResults(id)).be.null();
  } finally {
    await env.stopAll();
  }
});

test('progress and partial writes are throttled to one per interval (plus a trailing write) and visible via status()', async function () {
  let release;
  const gate = new Promise(function (r) { release = r; });
  const env = setup({
    check: { progress_min_interval_ms: 200 },
    run: async function (request, ctx) {
      for (let i = 1; i <= 20; i++) {
        ctx.progress({ done: i, total: 20, stage: 'pangenome', running: ['sorghum_' + i] });
        ctx.partial({ upto: i });
      }
      await gate;
      ctx.progress({ done: 20, total: 20, stage: 'annotate', running: [] });
      return { final: true };
    }
  });
  const writes = { progress: [], partial: [] };
  const setJob = env.store.setJob;
  env.store.setJob = function (doc, o) {
    if (o && o.expect && o.expect[0] === 'running') writes.progress.push(Date.now());
    return setJob(doc, o);
  };
  const setPartial = env.store.setPartial;
  env.store.setPartial = function (id, json, o) {
    writes.partial.push(Date.now());
    return setPartial(id, json, o);
  };
  const id = await env.add(1, 'pangenome');
  const w = env.worker();
  w.start();
  try {
    const running = await waitFor(async function () {
      const s = await env.status(id);
      return s.status === 'running' && s.partial && s;
    }, 3000, 'partial visible');
    running.results.should.have.property('upto');
    await sleep(450);
    writes.progress.length.should.equal(2);
    writes.partial.length.should.equal(2);
    (writes.progress[1] - writes.progress[0]).should.be.aboveOrEqual(190);
    (writes.partial[1] - writes.partial[0]).should.be.aboveOrEqual(190);
    const s = await env.status(id);
    s.progress.should.eql({ done: 20, total: 20, stage: 'pangenome', running: ['sorghum_20'] });
    s.results.should.eql({ upto: 20 });
    release();
    const done = await waitFor(function () { return statusIs(env, id, 'done'); }, 3000, 'done');
    done.results.should.eql({ final: true });
    should(await env.store.getPartial(id)).be.null();
  } finally {
    release();
    await env.stopAll();
  }
});

test('ctx.partial snapshots the results at call time (run() keeps mutating its live document)', async function () {
  let release;
  const gate = new Promise(function (r) { release = r; });
  const env = setup({
    run: async function (request, ctx) {
      const live = { pangenome: { pairs: [{ id: 'P1', summary: { genomes_total: 1 } }] } };
      ctx.partial(live);
      live.pangenome.pairs[0].summary.genomes_total = 999;
      live.half_written = true;
      ctx.partial({ cyclic: (function () { const o = {}; o.self = o; return o; })() }); // unserializable: ignored
      await gate;
      return { final: true };
    }
  });
  const id = await env.add(1, 'pangenome');
  const w = env.worker();
  w.start();
  try {
    const s = await waitFor(async function () {
      const x = await env.status(id);
      return x.status === 'running' && x.partial && x;
    }, 3000, 'partial visible');
    await sleep(120); // past progress_min_interval_ms: no later write may replace the snapshot
    (await env.status(id)).results.should.eql({ pangenome: { pairs: [{ id: 'P1', summary: { genomes_total: 1 } }] } });
    s.results.should.eql({ pangenome: { pairs: [{ id: 'P1', summary: { genomes_total: 1 } }] } });
    release();
    await waitFor(function () { return statusIs(env, id, 'done'); }, 3000, 'done');
  } finally {
    release();
    await env.stopAll();
  }
});

test('ctx.spawnLines: allow-list [blastn, blastdbcmd], niced, child SIGTERMed when the job aborts', async function () {
  const expectedNice = String(Math.min(19, os.getPriority() + 10));
  const out = {};
  let childPid = null;
  const env = setup({
    check: { job_timeout_ms: 1500, nice: 10 },
    run: async function (request, ctx) {
      if (request.pairs[0].id === 'basic') {
        try {
          await ctx.spawnLines('/bin/echo', ['hi'], {});
        } catch (e) {
          out.notAllowed = e.code;
        }
        try {
          await ctx.spawnLines(ctx.config.blastn, ['-e', 1], {});
        } catch (e) {
          out.badArgs = e.name;
        }
        const lines = [];
        const r = await ctx.spawnLines(ctx.config.blastn, ['-e', 'console.log(require("os").getPriority()); console.log(process.env.HOME === undefined)'],
          { onLine: function (l) { lines.push(l); } });
        out.lines = lines;
        out.code = r.code;
        return { ok: true };
      }
      const res = await ctx.spawnLines(ctx.config.blastn, ['-e', 'console.log(process.pid); setInterval(function () {}, 1000)'],
        { onLine: function (l) { childPid = Number(l); } });
      out.abort = res;
      return { unreachable: true };
    }
  });
  const basic = await env.add(1, 'specificity', 'basic');
  const w = env.worker();
  w.start();
  try {
    await waitFor(function () { return statusIs(env, basic, 'done'); }, 5000, 'basic spawn job');
    out.notAllowed.should.equal('SPAWN_NOT_ALLOWED');
    out.badArgs.should.equal('TypeError');
    out.lines.should.eql([expectedNice, 'true']); // niced, minimal env
    out.code.should.equal(0);

    const aborting = await env.add(2, 'specificity', 'abort');
    await waitFor(function () { return childPid; }, 5000, 'child pid');
    const s = await waitFor(function () { return statusIs(env, aborting, 'error'); }, 5000, 'job timeout');
    s.error.code.should.equal('JOB_TIMEOUT');
    await waitFor(function () {
      try {
        process.kill(childPid, 0);
        return false;
      } catch (e) {
        return e.code === 'ESRCH';
      }
    }, 3000, 'child killed');
    await waitFor(function () { return out.abort; }, 3000, 'spawnLines resolved');
    out.abort.aborted.should.be.true();
  } finally {
    await env.stopAll();
  }
});

test('worker lock: a second worker waits; it takes over after the first stops', async function () {
  const env = setup({ run: async function (request, ctx) { return { by: ctx.jobId }; } });
  const w1 = env.worker();
  const w2 = env.worker();
  w1.start();
  await waitFor(function () { return w1.haveLock; }, 2000, 'w1 lock');
  w2.start();
  try {
    await sleep(100);
    w2.haveLock.should.be.false();
    const a = await env.add(1);
    await waitFor(function () { return statusIs(env, a, 'done'); }, 3000, 'job on w1');
    w1.stats.done.should.equal(1);
    await w1.stop();
    await waitFor(function () { return w2.haveLock; }, 2000, 'w2 takes over');
    const b = await env.add(2);
    await waitFor(function () { return statusIs(env, b, 'done'); }, 3000, 'job on w2');
    w2.stats.done.should.equal(1);
  } finally {
    await env.stopAll();
  }
});

test('startup recovery: a job left running by a dead worker is requeued and finished (attempts 2)', async function () {
  const env = setup({ run: async function (request, ctx) { return { recovered: true, attempt: ctx.attempt }; } });
  const id = await env.add(1);
  (await env.store.claim(Date.now())).id.should.equal(id);
  const dead = workerLib.startedDoc(await env.store.getJob(id), { worker: 'deadhost:1:ffff', now: Date.now() });
  (await env.store.setJob(dead, { expect: ['queued'] })).should.be.true();
  const w = env.worker();
  w.start();
  try {
    const s = await waitFor(function () { return statusIs(env, id, 'done'); }, 3000, 'recovered');
    s.attempts.should.equal(2);
    s.results.should.eql({ recovered: true, attempt: 2 });
    w.stats.requeued.should.equal(1);
  } finally {
    await env.stopAll();
  }
});

test('local_max_jobs: at most 2 jobs run at once in one worker', async function () {
  let concurrent = 0;
  let peak = 0;
  const gates = [];
  const env = setup({
    check: { local_max_jobs: 2, global_max_jobs: 5 },
    run: async function () {
      concurrent++;
      peak = Math.max(peak, concurrent);
      await new Promise(function (r) { gates.push(r); });
      concurrent--;
      return {};
    }
  });
  const ids = [await env.add(1), await env.add(2), await env.add(3)];
  const w = env.worker();
  w.start();
  try {
    await waitFor(function () { return gates.length === 2; }, 3000, 'two running');
    await sleep(100);
    gates.length.should.equal(2);
    (await env.status(ids[2])).status.should.equal('queued');
    gates.shift()();
    await waitFor(function () { return gates.length === 2; }, 3000, 'third started');
    gates.shift()();
    gates.shift()();
    for (const id of ids) await waitFor(function () { return statusIs(env, id, 'done'); }, 3000, 'done ' + id);
    peak.should.equal(2);
  } finally {
    gates.forEach(function (g) { g(); });
    await env.stopAll();
  }
});

test('lock lost to another worker: running jobs are abandoned without writes, then recovered', async function () {
  const env = setup({
    run: function (request, ctx) {
      return ctx.attempt === 1 ? blockUntilAbort(ctx) : Promise.resolve({ attempt: ctx.attempt });
    }
  });
  const id = await env.add(1);
  const w = env.worker();
  w.start();
  try {
    await waitFor(function () { return statusIs(env, id, 'running'); }, 3000, 'running');
    (await env.store.releaseWorkerLock(w.token)).should.be.true();
    (await env.store.acquireWorkerLock('thief', 60000)).should.be.true();
    await waitFor(function () { return !w.haveLock && w.activeIds().length === 0; }, 2000, 'abandoned');
    w.stats.abandoned.should.equal(1);
    const s = await env.status(id);
    s.status.should.equal('running'); // no requeue/fail from a worker that lost its lock
    (await env.store.getJob(id)).worker.should.equal(w.token);
    (await env.store.releaseWorkerLock('thief')).should.be.true();
    const done = await waitFor(function () { return statusIs(env, id, 'done'); }, 3000, 'recovered after re-acquiring');
    done.attempts.should.equal(2);
    done.results.should.eql({ attempt: 2 });
  } finally {
    await env.stopAll();
  }
});

test('fatal MONGO_UNAVAILABLE: the job is requeued and onFatal is called', async function () {
  let fatal = null;
  let w = null;
  const env = setup({
    run: async function () {
      const e = new Error('mongoCollection() resolved undefined');
      e.code = 'MONGO_UNAVAILABLE';
      throw e;
    }
  });
  const id = await env.add(1);
  w = env.worker({ onFatal: function (err) { fatal = err; w.stop(); } });
  w.start();
  try {
    await waitFor(function () { return fatal; }, 3000, 'onFatal');
    fatal.code.should.equal('MONGO_UNAVAILABLE');
    const s = await env.status(id);
    s.status.should.equal('queued');
    s.attempts.should.equal(1);
    workerLib.isFatal({ fatal: true }).should.be.true();
    workerLib.isFatal({ code: 'CHECK_FAILED' }).should.be.false();
  } finally {
    await env.stopAll();
  }
});

test('periodic reconciliation (staleOnly) requeues only orphaned running entries whose heartbeat is older than stale_ms', async function () {
  const env = setup({ run: async function () { return {}; } });
  const orphan = await env.add(1);
  const t0 = Date.now();
  (await env.store.claim(t0)).id.should.equal(orphan); // claimed, but the start never happened (e.g. a Redis error)
  try {
    const early = env.worker({ now: function () { return t0 + 30000; } });
    (await early.requeueRunning({ staleOnly: true })).should.eql({ requeued: [], failed: [], cleaned: [] });
    env.store._state().running.map(function (e) { return e[0]; }).should.eql([orphan]);
    const late = env.worker({ now: function () { return t0 + 60000; } });
    (await late.requeueRunning({ staleOnly: true })).should.eql({ requeued: [orphan], failed: [], cleaned: [] });
    const s = await env.status(orphan);
    s.status.should.equal('queued');
    s.queue_position.should.equal(0);
    env.store._state().running.should.eql([]);
  } finally {
    await env.stopAll();
  }
});

test('ctx.spawnLines: env has NCBI_DONT_USE_NCBIRC=1 and cwd is the job\'s private 0700 tmpdir under tmp_dir (created lazily), never tmpdir itself', async function (t) {
  const base = fs.mkdtempSync(path.join(os.tmpdir(), 'primers-wtmp-base-'));
  t.after(function () { fs.rmSync(base, { recursive: true, force: true }); });
  workerLib.CHILD_ENV.should.eql({ PATH: '/usr/bin:/bin', NCBI_DONT_USE_NCBIRC: '1' });
  const out = {};
  const env = setup({
    top: { tmp_dir: base },
    run: async function (request, ctx) {
      out.before = fs.readdirSync(base);
      const lines = [];
      const r = await ctx.spawnLines(ctx.config.blastn, ['-e', 'const fs = require("fs"); console.log(process.cwd()); ' +
        'console.log((fs.statSync(".").mode & 0o777).toString(8)); console.log(Object.keys(process.env).sort().join(",")); ' +
        'console.log(process.env.NCBI_DONT_USE_NCBIRC)'], { onLine: function (l) { lines.push(l); } });
      out.code = r.code;
      out.lines = lines;
      out.tmp = fs.realpathSync(ctx.tmpdir());
      return { ok: true };
    }
  });
  const id = await env.add(1);
  const w = env.worker();
  w.start();
  try {
    await waitFor(function () { return statusIs(env, id, 'done'); }, 5000, 'spawn job');
    out.before.should.eql([]); // nothing created before the first spawn
    out.code.should.equal(0);
    out.lines.should.eql([out.tmp, '700', 'NCBI_DONT_USE_NCBIRC,PATH', '1']);
    out.tmp.should.not.equal(fs.realpathSync(os.tmpdir()));
    path.dirname(out.tmp).should.equal(fs.realpathSync(base));
    path.basename(out.tmp).indexOf('primers-check-' + id.slice(0, 8) + '-').should.equal(0);
    await waitFor(function () { return !fs.existsSync(out.tmp); }, 2000, 'tmpdir removal');
  } finally {
    await env.stopAll();
  }
});

test('fatal on the last attempt (max_attempts): the job fails with its code instead of being requeued again; onFatal still fires', async function () {
  const fatals = [];
  const attempts = [];
  const env = setup({
    check: { max_attempts: 2 },
    run: async function (request, ctx) {
      attempts.push(ctx.attempt);
      const e = new Error('gene annotation (mongo) is unavailable');
      e.code = 'MONGO_UNAVAILABLE';
      e.fatal = true;
      throw e;
    }
  });
  const id = await env.add(1);
  const w = env.worker({ onFatal: function (err) { fatals.push(err.code); } });
  w.start();
  try {
    const s = await waitFor(function () { return statusIs(env, id, 'error'); }, 3000, 'failed after 2 fatal attempts');
    s.error.code.should.equal('MONGO_UNAVAILABLE');
    s.attempts.should.equal(2);
    attempts.should.eql([1, 2]);
    await waitFor(function () { return fatals.length === 2; }, 1000, 'onFatal twice');
    fatals.should.eql(['MONGO_UNAVAILABLE', 'MONGO_UNAVAILABLE']);
    w.stats.requeued.should.equal(1);
    w.stats.failed.should.equal(1);
    env.store._state().queues.spec.should.eql([]);
    env.store._state().running.should.eql([]);
  } finally {
    await env.stopAll();
  }
});

test('claim aborted by stop(): the claimed job goes back to the FRONT of its queue by requeue alone (no release before it)', async function () {
  let ran = 0;
  const env = setup({ run: async function () { ran++; return {}; } });
  const first = await env.add(1);
  const second = await env.add(2);
  const w = env.worker();
  const calls = [];
  const realClaim = env.store.claim;
  const realRelease = env.store.release;
  const realRequeue = env.store.requeue;
  let stopP = null;
  env.store.claim = async function (t) {
    const r = await realClaim.call(env.store, t);
    if (r && !stopP) stopP = w.stop(); // stopping flips while the claim is in flight
    return r;
  };
  env.store.release = async function (id) { calls.push('release ' + id); return realRelease.call(env.store, id); };
  env.store.requeue = async function (doc, o) { calls.push('requeue ' + doc.job_id); return realRequeue.call(env.store, doc, o); };
  w.start();
  await waitFor(function () { return stopP !== null; }, 3000, 'claim + stop');
  await stopP;
  ran.should.equal(0);
  calls.should.eql(['requeue ' + first]);
  const st = env.store._state();
  st.queues.spec.should.eql([first, second]);
  st.running.should.eql([]);
  st.slots.should.eql([]);
  const s = await env.status(first);
  s.status.should.equal('queued');
  s.queue_position.should.equal(0);
});

test('claim aborted by stop() and getJob throws: the job stays in the running set (not stranded) and the next worker runs it', async function () {
  const env = setup({ run: async function (request, ctx) { return { attempt: ctx.attempt }; } });
  const id = await env.add(1);
  const w1 = env.worker();
  const realClaim = env.store.claim;
  const realGetJob = env.store.getJob;
  let stopP = null;
  let failGet = false;
  env.store.claim = async function (t) {
    const r = await realClaim.call(env.store, t);
    if (r && !stopP) {
      failGet = true;
      stopP = w1.stop();
    }
    return r;
  };
  env.store.getJob = async function (jobId) {
    if (failGet) {
      failGet = false;
      throw new Error('the job store is closed');
    }
    return realGetJob.call(env.store, jobId);
  };
  w1.start();
  await waitFor(function () { return stopP !== null; }, 3000, 'claim + stop');
  await stopP;
  const st = env.store._state();
  st.running.map(function (e) { return e[0]; }).should.eql([id]);
  st.queues.spec.should.eql([]);
  (await realGetJob.call(env.store, id)).status.should.equal('queued');
  env.store.claim = realClaim;
  env.store.getJob = realGetJob;
  const w2 = env.worker();
  w2.start();
  try {
    const s = await waitFor(function () { return statusIs(env, id, 'done'); }, 3000, 'recovered by the next worker');
    s.results.should.eql({ attempt: 1 });
    w2.stats.requeued.should.equal(1);
  } finally {
    await env.stopAll();
  }
});

// ---- worker_main -----------------------------------------------------------------------------------

function captureLog() {
  const lines = [];
  const push = function () { lines.push(Array.prototype.slice.call(arguments).join(' ')); };
  return { lines: lines, info: push, warn: push, error: push, log: push };
}

function fakeConfig(cfg, siteKey) {
  return { get: function () { return cfg; }, siteKey: function () { return siteKey; }, basePath: function () { return '/sorghum_v11'; } };
}

test('worker_main: logs site_key first and exits 75 when mongo is unavailable or the probe hangs', async function () {
  const cfg = makeCfg();
  const store = createMemoryStore({ cfg: cfg, siteKey: 'wmain', shared: { slots: new Map(), panSlots: new Map() } });
  let code = null;
  let log = captureLog();
  const r = await workerMain.main({
    config: fakeConfig(cfg, 'wmain'), log: log, exit: function (c) { code = c; }, store: store, signals: false,
    mongo: { genes: { mongoCollection: async function () { return undefined; } } }
  });
  should(r).be.undefined();
  code.should.equal(75);
  log.lines[0].should.equal('primers site_key=wmain');
  log.lines.join('\n').should.match(/mongo is unavailable/);

  code = null;
  log = captureLog();
  await workerMain.main({
    config: fakeConfig(cfg, 'wmain'), log: log, exit: function (c) { code = c; }, store: store, signals: false,
    mongo: { genes: { mongoCollection: function () { return new Promise(function () {}); } } }, mongoProbeTimeoutMs: 50
  });
  code.should.equal(75);
  should(store._state().lock).be.null();
});

test('worker_main: shutdown requeues the running job and exits 0 within 1 s; a fatal job error exits 75', async function () {
  const cfg = makeCfg();
  const shared = { slots: new Map(), panSlots: new Map() };
  const store = createMemoryStore({ cfg: cfg, siteKey: 'wmain', shared: shared });
  const id = hexId(500);
  await store.submit(jobs.newJobDoc(id, fakeNorm(), Date.now()));
  let code = null;
  let exitAt = 0;
  const r = await workerMain.main({
    config: fakeConfig(cfg, 'wmain'), log: quiet, exit: function (c) { code = c; exitAt = Date.now(); },
    store: store, signals: false, timing: FAST,
    mongo: { genes: { mongoCollection: async function () { return {}; } } },
    check: { ALGORITHM_VERSION: '1', run: function (request, ctx) { return blockUntilAbort(ctx); } }
  });
  await waitFor(async function () { const j = await store.getJob(id); return j.status === 'running'; }, 3000, 'running');
  const t0 = Date.now();
  await r.shutdown('SIGTERM', 0);
  code.should.equal(0);
  (exitAt - t0).should.be.below(1000);
  const j = await store.getJob(id);
  j.status.should.equal('queued');
  store._state().queues.spec.should.eql([id]);
  should(store._state().lock).be.null();

  // fatal
  const store2 = createMemoryStore({ cfg: cfg, siteKey: 'wmain2', shared: { slots: new Map(), panSlots: new Map() } });
  const id2 = hexId(501);
  await store2.submit(jobs.newJobDoc(id2, fakeNorm(), Date.now()));
  code = null;
  await workerMain.main({
    config: fakeConfig(cfg, 'wmain2'), log: quiet, exit: function (c) { code = c; }, store: store2, signals: false, timing: FAST,
    mongo: { genes: { mongoCollection: async function () { return {}; } } },
    check: {
      ALGORITHM_VERSION: '1',
      run: async function () { const e = new Error('mongo went away'); e.code = 'MONGO_UNAVAILABLE'; throw e; }
    }
  });
  await waitFor(function () { return code !== null; }, 3000, 'exit 75');
  code.should.equal(75);
  (await store2.getJob(id2)).status.should.equal('queued');
});

test('worker_main: run() throwing {code: MONGO_UNAVAILABLE, fatal: true} under ctx.fatalOnMongoUnavailable requeues the job and exits 75', async function () {
  const cfg = makeCfg();
  const store = createMemoryStore({ cfg: cfg, siteKey: 'wfatal', shared: { slots: new Map(), panSlots: new Map() } });
  const id = hexId(502);
  await store.submit(jobs.newJobDoc(id, fakeNorm(), Date.now()));
  let code = null;
  const flags = [];
  const log = captureLog();
  await workerMain.main({
    config: fakeConfig(cfg, 'wfatal'), log: log, exit: function (c) { code = c; }, store: store, signals: false, timing: FAST,
    mongo: { genes: { mongoCollection: async function () { return {}; } } },
    check: {
      ALGORITHM_VERSION: '2',
      run: async function (request, ctx) {
        flags.push(ctx.fatalOnMongoUnavailable);
        const e = new Error('gene annotation (mongo) is unavailable');
        e.code = 'MONGO_UNAVAILABLE';
        e.fatal = true;
        throw e;
      }
    }
  });
  await waitFor(function () { return code !== null; }, 3000, 'exit 75');
  code.should.equal(75);
  flags.should.eql([true]);
  const j = await store.getJob(id);
  j.status.should.equal('queued');
  j.attempts.should.equal(1);
  store._state().queues.spec.should.eql([id]);
  store._state().running.should.eql([]);
  log.lines.join('\n').should.match(/fatal error MONGO_UNAVAILABLE \(requeued\)/);
});

test('worker_main: disabled or a non-redis store idles instead of running jobs', async function () {
  let code = null;
  const off = makeCfg({}, { enabled: false });
  const r1 = await workerMain.main({ config: fakeConfig(off, 'wmain'), log: quiet, exit: function (c) { code = c; }, signals: false });
  r1.idle.should.be.true();
  r1.stop();
  code.should.equal(0);
  code = null;
  const mem = makeCfg({ store: 'memory' });
  const r2 = await workerMain.main({ config: fakeConfig(mem, 'wmain'), log: quiet, exit: function (c) { code = c; }, signals: false });
  r2.idle.should.be.true();
  r2.stop();
  code.should.equal(0);
});

test('node worker_main.js from cwd / derives site_key sorghum_v11:sorghum11 from swagger.yaml and exits 0 on SIGTERM', async function () {
  const env = Object.assign({}, process.env, { PRIMERS_JOB_STORE: 'memory' });
  delete env.PRIMERS_SITE_KEY;
  delete env.PRIMERS_ENABLED;
  const child = spawn(process.execPath, [WORKER_MAIN], { cwd: '/', env: env, stdio: ['ignore', 'pipe', 'pipe'] });
  let stdout = '';
  let stderr = '';
  child.stdout.on('data', function (c) { stdout += c; });
  child.stderr.on('data', function (c) { stderr += c; });
  const exited = new Promise(function (resolve) { child.on('exit', function (c, s) { resolve({ code: c, signal: s, at: Date.now() }); }); });
  try {
    await waitFor(function () { return /idling/.test(stderr); }, 15000, 'worker_main idle line');
    stdout.split('\n')[0].should.equal('primers site_key=sorghum_v11:sorghum11');
    const t0 = Date.now();
    child.kill('SIGTERM');
    const ex = await exited;
    ex.code.should.equal(0);
    (ex.at - t0).should.be.below(1000);
  } finally {
    if (child.exitCode === null && child.signalCode === null) child.kill('SIGKILL');
  }
});
