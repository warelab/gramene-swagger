'use strict';

// Counting semaphore with a bounded FIFO waiter queue and a wait timeout.
// A full queue or an expired wait rejects with 503 BUSY {retry_after_s}.

const { PrimerHttpError } = require('./errors');

function positiveInt(v, name, allowZero) {
  if (!Number.isInteger(v) || v < (allowZero ? 0 : 1)) {
    throw new TypeError('semaphore ' + name + ' must be an integer >= ' + (allowZero ? 0 : 1));
  }
  return v;
}

function abortError(signal) {
  const r = signal && signal.reason;
  if (r instanceof Error) return r;
  return new PrimerHttpError(504, 'DEADLINE_EXCEEDED', 'request deadline exceeded', {});
}

class Semaphore {
  constructor(opts) {
    opts = opts || {};
    this.max = positiveInt(opts.max === undefined ? 4 : opts.max, 'max');
    this.maxWaiting = positiveInt(opts.maxWaiting === undefined ? 16 : opts.maxWaiting, 'maxWaiting', true);
    const wt = opts.waitTimeoutMs === undefined ? 10000 : opts.waitTimeoutMs;
    if (!(wt === Infinity || (Number.isFinite(wt) && wt >= 0))) throw new TypeError('semaphore waitTimeoutMs must be >= 0');
    this.waitTimeoutMs = wt;
    this.retryAfterS = opts.retryAfterS === undefined ? 5 : opts.retryAfterS;
    this._running = 0;
    this._waiters = [];
  }

  get running() { return this._running; }
  get waiting() { return this._waiters.length; }

  _busy(reason) {
    return new PrimerHttpError(503, 'BUSY', 'too many concurrent requests; retry shortly',
      { retry_after_s: this.retryAfterS, reason: reason });
  }

  _releaser() {
    let released = false;
    return () => {
      if (released) return;
      released = true;
      this._release();
    };
  }

  _release() {
    const w = this._waiters.shift();
    if (w) {
      this._cleanup(w);
      w.resolve(this._releaser()); // the slot passes straight to the next waiter
    } else {
      this._running--;
    }
  }

  _cleanup(w) {
    if (w.timer) clearTimeout(w.timer);
    if (w.signal && w.onAbort) w.signal.removeEventListener('abort', w.onAbort);
  }

  _remove(w) {
    const i = this._waiters.indexOf(w);
    if (i >= 0) this._waiters.splice(i, 1);
    this._cleanup(w);
  }

  // Resolves to an idempotent release() function.
  acquire(opts) {
    opts = opts || {};
    const signal = opts.signal;
    if (signal && signal.aborted) return Promise.reject(abortError(signal));
    if (this._running < this.max && this._waiters.length === 0) {
      this._running++;
      return Promise.resolve(this._releaser());
    }
    if (this._waiters.length >= this.maxWaiting) return Promise.reject(this._busy('queue_full'));
    const waitMs = opts.timeoutMs === undefined ? this.waitTimeoutMs : opts.timeoutMs;
    return new Promise((resolve, reject) => {
      const w = { resolve: resolve, reject: reject, timer: null, signal: signal, onAbort: null };
      if (waitMs !== Infinity) {
        w.timer = setTimeout(() => {
          this._remove(w);
          reject(this._busy('wait_timeout'));
        }, waitMs);
      }
      if (signal) {
        w.onAbort = () => {
          this._remove(w);
          reject(abortError(signal));
        };
        signal.addEventListener('abort', w.onAbort, { once: true });
      }
      this._waiters.push(w);
    });
  }

  // Run fn() while holding a permit; the permit is released however fn settles.
  async run(fn, opts) {
    const release = await this.acquire(opts);
    try {
      return await fn();
    } finally {
      release();
    }
  }
}

let designSem = null;

// Process-wide design semaphore sized from primers.design (created on first use).
function designSemaphore() {
  if (!designSem) {
    const d = require('./config').get().design;
    designSem = new Semaphore({ max: d.max_concurrent, maxWaiting: d.max_waiting, waitTimeoutMs: d.wait_timeout_ms, retryAfterS: 5 });
  }
  return designSem;
}

function _resetDesignSemaphore() {
  designSem = null;
}

module.exports = { Semaphore, designSemaphore, _resetDesignSemaphore };
