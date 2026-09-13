'use strict';

const test = require('node:test');
const should = require('should');

const { Semaphore, designSemaphore, _resetDesignSemaphore } = require('../../../api/helpers/primers/semaphore');
const { PrimerHttpError } = require('../../../api/helpers/primers/errors');
const config = require('../../../api/helpers/primers/config');

function deferred() {
  let resolve;
  const promise = new Promise(function (r) { resolve = r; });
  return { promise: promise, resolve: resolve };
}

const tick = function () { return new Promise(function (r) { setImmediate(r); }); };

test('never runs more than max tasks and serves waiters FIFO', async function () {
  const sem = new Semaphore({ max: 2, maxWaiting: 10, waitTimeoutMs: 5000 });
  let active = 0;
  let peak = 0;
  const started = [];
  const gates = [0, 1, 2, 3, 4].map(function () { return deferred(); });
  const runs = gates.map(function (g, i) {
    return sem.run(async function () {
      started.push(i);
      active++;
      peak = Math.max(peak, active);
      await g.promise;
      active--;
      return i;
    });
  });
  await tick();
  sem.running.should.equal(2);
  sem.waiting.should.equal(3);
  started.should.eql([0, 1]);
  gates[1].resolve();
  await tick();
  started.should.eql([0, 1, 2]);
  gates[0].resolve();
  await tick();
  started.should.eql([0, 1, 2, 3]);
  gates[2].resolve();
  gates[3].resolve();
  gates[4].resolve();
  (await Promise.all(runs)).should.eql([0, 1, 2, 3, 4]);
  peak.should.equal(2);
  sem.running.should.equal(0);
  sem.waiting.should.equal(0);
});

test('a full waiter queue rejects immediately with 503 BUSY retry_after_s 5', async function () {
  const sem = new Semaphore({ max: 1, maxWaiting: 1, waitTimeoutMs: 5000 });
  const release = await sem.acquire();
  const waiter = sem.acquire();
  let err;
  try {
    await sem.acquire();
  } catch (e) {
    err = e;
  }
  err.should.be.instanceOf(PrimerHttpError);
  err.status.should.equal(503);
  err.code.should.equal('BUSY');
  err.details.retry_after_s.should.equal(5);
  err.details.reason.should.equal('queue_full');
  release();
  const release2 = await waiter;
  sem.running.should.equal(1);
  release2();
  sem.running.should.equal(0);
});

test('an expired wait rejects with 503 BUSY and leaves the queue clean', async function () {
  const sem = new Semaphore({ max: 1, maxWaiting: 4, waitTimeoutMs: 30 });
  const release = await sem.acquire();
  const t0 = Date.now();
  let err;
  try {
    await sem.acquire();
  } catch (e) {
    err = e;
  }
  (Date.now() - t0).should.be.aboveOrEqual(25);
  err.code.should.equal('BUSY');
  err.details.reason.should.equal('wait_timeout');
  sem.waiting.should.equal(0);
  sem.running.should.equal(1);
  release();
  sem.running.should.equal(0);
  // per-call timeout override
  const r = await sem.acquire({ timeoutMs: 1 });
  r();
});

test('release is idempotent and run() releases when the task throws', async function () {
  const sem = new Semaphore({ max: 1, maxWaiting: 2, waitTimeoutMs: 1000 });
  const release = await sem.acquire();
  release();
  release();
  sem.running.should.equal(0);
  const r1 = await sem.acquire();
  sem.running.should.equal(1);
  let busy;
  try {
    await sem.acquire({ timeoutMs: 5 });
  } catch (e) {
    busy = e;
  }
  busy.code.should.equal('BUSY');
  r1();
  let taskErr;
  try {
    await sem.run(async function () { throw new Error('task failed'); });
  } catch (e) {
    taskErr = e;
  }
  taskErr.message.should.equal('task failed');
  sem.running.should.equal(0);
  (await sem.run(function () { return 42; })).should.equal(42);
  sem.running.should.equal(0);
});

test('aborting a waiting acquire rejects with the abort reason and removes the waiter', async function () {
  const sem = new Semaphore({ max: 1, maxWaiting: 4, waitTimeoutMs: 5000 });
  const release = await sem.acquire();
  const ac = new AbortController();
  const p = sem.acquire({ signal: ac.signal });
  sem.waiting.should.equal(1);
  const reason = new PrimerHttpError(504, 'DEADLINE_EXCEEDED', 'deadline');
  ac.abort(reason);
  let err;
  try {
    await p;
  } catch (e) {
    err = e;
  }
  err.should.equal(reason);
  sem.waiting.should.equal(0);

  const ac2 = new AbortController();
  ac2.abort();
  let err2;
  try {
    await sem.acquire({ signal: ac2.signal });
  } catch (e) {
    err2 = e;
  }
  err2.should.be.instanceOf(Error);
  release();
  sem.running.should.equal(0);
});

test('constructor validates its limits', function () {
  (function () { return new Semaphore({ max: 0 }); }).should.throw(/max/);
  (function () { return new Semaphore({ maxWaiting: -1 }); }).should.throw(/maxWaiting/);
  (function () { return new Semaphore({ waitTimeoutMs: -5 }); }).should.throw(/waitTimeoutMs/);
  const s = new Semaphore();
  s.max.should.equal(4);
  s.maxWaiting.should.equal(16);
  s.waitTimeoutMs.should.equal(10000);
});

test('designSemaphore is sized from primers.design', function () {
  try {
    _resetDesignSemaphore();
    config._setForTests(null);
    const s = designSemaphore();
    s.max.should.equal(4);
    s.maxWaiting.should.equal(16);
    s.waitTimeoutMs.should.equal(10000);
    designSemaphore().should.equal(s);
    _resetDesignSemaphore();
    config._setForTests({ design: { max_concurrent: 3, max_waiting: 2, wait_timeout_ms: 50 } });
    const t = designSemaphore();
    t.max.should.equal(3);
    t.maxWaiting.should.equal(2);
    t.waitTimeoutMs.should.equal(50);
  } finally {
    _resetDesignSemaphore();
    config._setForTests(null);
  }
});
