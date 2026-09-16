'use strict';

// Spawn primer3_core for one Boulder-IO record, with a timeout, abort signal and output caps.
// No shell, argument array, minimal environment, cwd a private mkdtemp directory (mode 0700) removed afterwards,
// never a shared directory such as /tmp itself.

const { spawn } = require('child_process');
const fs = require('fs');
const os = require('os');
const path = require('path');
const boulder = require('./boulder');
const { PrimerHttpError, redactPaths } = require('./errors');

const STDERR_CAP = 64 * 1024;
const LOG_STDERR_CHARS = 500;
const UNAVAILABLE = new Set(['ENOENT', 'EACCES', 'EPERM', 'ENOEXEC', 'ENOTDIR']);
const TRANSIENT = new Set(['EAGAIN', 'EMFILE', 'ENFILE', 'ENOMEM']);
const SPAWN_ENV = Object.freeze({ PATH: '/usr/bin:/bin' });

function cfg() {
  return require('./config').get();
}

function logError(log, msg) {
  try {
    (log.error || log.log).call(log, msg);
  } catch (e) { /* ignore */ }
}

function spawnError(err, log) {
  const code = err && err.code;
  if (UNAVAILABLE.has(code)) {
    logError(log, 'primer3_core unavailable: ' + code);
    return new PrimerHttpError(503, 'PRIMER3_UNAVAILABLE', 'Primer3 is not available on this server', { retry_after_s: 60 });
  }
  if (TRANSIENT.has(code)) {
    logError(log, 'primer3_core could not start: ' + code);
    return new PrimerHttpError(503, 'PRIMER3_UNAVAILABLE', 'Primer3 could not be started; retry shortly', { retry_after_s: 5 });
  }
  logError(log, 'primer3_core spawn failed: ' + redactPaths(String(err && err.message)));
  return new PrimerHttpError(500, 'PRIMER3_FAILED', 'Primer3 could not be started', {});
}

function deadlineError(timeoutMs) {
  return new PrimerHttpError(504, 'DEADLINE_EXCEEDED', 'Primer3 did not finish within the time limit', { timeout_ms: timeoutMs });
}

function abortError(signal) {
  const r = signal && signal.reason;
  if (r instanceof PrimerHttpError) return r;
  return new PrimerHttpError(504, 'DEADLINE_EXCEEDED', 'request deadline exceeded', {});
}

// Private working directory under base (absolute, else os.tmpdir()); mkdtemp creates it with mode 0700.
function makeWorkDir(base) {
  const root = typeof base === 'string' && path.isAbsolute(base) ? base : os.tmpdir();
  return fs.promises.mkdtemp(path.join(root, 'primers-primer3-'));
}

function removeWorkDir(dir) {
  return fs.promises.rm(dir, { recursive: true, force: true }).catch(function () { /* best effort */ });
}

// run(record, {timeoutMs, maxStdoutBytes, signal, bin, log, tmpDir})
//   record: serialized Boulder-IO string, or a tags object/array (serialized here).
//   tmpDir: absolute base of the private working directory (default os.tmpdir()).
// Resolves {tags, error, warning, exitCode, elapsedMs, stdoutBytes, stderrTail}.
//   error = PRIMER_ERROR (primer3 exits non-zero for global errors but still prints the record;
//   the caller maps it to 400 PRIMER3_INPUT_ERROR); warning = PRIMER_WARNING.
// Rejects PrimerHttpError: 503 PRIMER3_UNAVAILABLE (ENOENT/EACCES...), 500 PRIMER3_FAILED
//   (non-zero exit, signal, truncated record, stdout over the cap, no working directory),
//   504 DEADLINE_EXCEEDED (timeout or abort; the child is SIGKILLed).
function run(record, opts) {
  opts = opts || {};
  const log = opts.log || console;
  let c = null;
  const conf = function () { return c || (c = cfg()); };
  const bin = opts.bin || conf().primer3_core;
  const timeoutMs = opts.timeoutMs === undefined ? conf().design.primer3_timeout_ms : opts.timeoutMs;
  const maxStdout = opts.maxStdoutBytes === undefined ? conf().design.max_stdout_bytes : opts.maxStdoutBytes;
  const signal = opts.signal;

  let input;
  try {
    input = typeof record === 'string' ? record : boulder.serialize(record);
  } catch (e) {
    return Promise.reject(e);
  }
  if (signal && signal.aborted) return Promise.reject(abortError(signal));
  if (!(timeoutMs > 0)) return Promise.reject(deadlineError(timeoutMs));
  const spawnOpts = { timeoutMs: timeoutMs, maxStdout: maxStdout, signal: signal, log: log };
  return makeWorkDir(opts.tmpDir).then(function (cwd) {
    return spawnRecord(bin, input, cwd, spawnOpts).finally(function () { return removeWorkDir(cwd); });
  }, function (err) {
    logError(log, 'primer3_core working directory could not be created: ' + (err && err.code));
    throw new PrimerHttpError(500, 'PRIMER3_FAILED', 'Primer3 could not be started', {});
  });
}

function spawnRecord(bin, input, cwd, o) {
  const log = o.log;
  const timeoutMs = o.timeoutMs;
  const maxStdout = o.maxStdout;
  const signal = o.signal;

  return new Promise(function (resolve, reject) {
    // The signal may have aborted while the working directory was created.
    if (signal && signal.aborted) return reject(abortError(signal));

    const started = Date.now();
    let child;
    try {
      child = spawn(bin, ['-strict_tags'], { cwd: cwd, env: SPAWN_ENV, stdio: ['pipe', 'pipe', 'pipe'], windowsHide: true });
    } catch (e) {
      return reject(spawnError(e, log));
    }

    const chunks = [];
    let outBytes = 0;
    let errTail = '';
    let settled = false;
    let failure = null;
    let timer = null;
    let onAbort = null;

    function cleanup() {
      if (timer) clearTimeout(timer);
      if (signal && onAbort) signal.removeEventListener('abort', onAbort);
    }
    function finish(err, value) {
      if (settled) return;
      settled = true;
      cleanup();
      if (err) {
        try { child.stdout.destroy(); child.stderr.destroy(); } catch (e) { /* ignore */ }
        reject(err);
      } else {
        resolve(value);
      }
    }
    function kill(err) {
      if (!failure) failure = err;
      try { child.kill('SIGKILL'); } catch (e) { /* already gone */ }
    }

    child.on('error', function (err) {
      if (failure) return finish(failure);
      finish(spawnError(err, log));
    });
    child.stdout.on('data', function (chunk) {
      if (failure) return;
      outBytes += chunk.length;
      if (outBytes > maxStdout) {
        logError(log, 'primer3_core stdout exceeded ' + maxStdout + ' bytes; killed');
        kill(new PrimerHttpError(500, 'PRIMER3_FAILED', 'Primer3 output exceeded the size limit', { limit_bytes: maxStdout }));
        return;
      }
      chunks.push(chunk);
    });
    child.stderr.on('data', function (chunk) {
      errTail = (errTail + chunk.toString('utf8')).slice(-STDERR_CAP);
    });
    child.stdin.on('error', function () { /* EPIPE when primer3 exits before reading everything */ });

    // After a kill, don't wait for 'close' (a grandchild could hold the pipes open).
    child.on('exit', function () {
      if (failure) finish(failure);
    });
    child.on('close', function (code, sig) {
      if (failure) return finish(failure);
      const { tags, complete } = boulder.parseRecord(Buffer.concat(chunks));
      const stderrTail = redactPaths(errTail.slice(-LOG_STDERR_CHARS));
      const result = {
        tags: tags,
        error: tags.PRIMER_ERROR || null,
        warning: tags.PRIMER_WARNING || null,
        exitCode: code,
        elapsedMs: Date.now() - started,
        stdoutBytes: outBytes,
        stderrTail: stderrTail
      };
      if (complete && result.error) return finish(null, result);
      if (code !== 0 || sig) {
        logError(log, 'primer3_core exited with ' + (sig ? 'signal ' + sig : 'code ' + code) + ': ' + stderrTail);
        return finish(new PrimerHttpError(500, 'PRIMER3_FAILED', 'Primer3 failed', {}));
      }
      if (!complete) {
        logError(log, 'primer3_core output ended without a record terminator: ' + stderrTail);
        return finish(new PrimerHttpError(500, 'PRIMER3_FAILED', 'Primer3 output was incomplete', {}));
      }
      finish(null, result);
    });

    if (timeoutMs !== Infinity) timer = setTimeout(function () { kill(deadlineError(timeoutMs)); }, timeoutMs);
    if (signal) {
      onAbort = function () { kill(abortError(signal)); };
      signal.addEventListener('abort', onAbort, { once: true });
    }
    child.stdin.end(input);
  });
}

const versionCache = new Map();

function unavailableError() {
  return new PrimerHttpError(503, 'PRIMER3_UNAVAILABLE', 'Primer3 is not available on this server', { retry_after_s: 60 });
}

function about(bin, timeoutMs, tmpDir) {
  return makeWorkDir(tmpDir).then(function (cwd) {
    return aboutIn(bin, timeoutMs, cwd).finally(function () { return removeWorkDir(cwd); });
  }, function () { throw unavailableError(); });
}

function aboutIn(bin, timeoutMs, cwd) {
  return new Promise(function (resolve, reject) {
    const unavailable = unavailableError;
    let child;
    try {
      child = spawn(bin, ['-about'], { cwd: cwd, env: SPAWN_ENV, stdio: ['ignore', 'pipe', 'ignore'], windowsHide: true });
    } catch (e) {
      return reject(unavailable());
    }
    let out = '';
    let done = false;
    const finish = function (err, v) {
      if (done) return;
      done = true;
      clearTimeout(timer);
      if (err) {
        try { child.kill('SIGKILL'); } catch (e) { /* ignore */ }
        reject(err);
      } else {
        resolve(v);
      }
    };
    const timer = setTimeout(function () { finish(unavailable()); }, timeoutMs);
    child.on('error', function () { finish(unavailable()); });
    child.stdout.on('data', function (d) { if (out.length < 4096) out += d.toString('utf8'); });
    child.on('close', function (code) {
      if (code !== 0) return finish(unavailable());
      const m = /release\s+([0-9][^\s]*)/.exec(out);
      const v = m ? m[1] : out.trim().split('\n')[0];
      if (!v) return finish(unavailable());
      finish(null, v);
    });
  });
}

// Primer3 version string (e.g. "2.6.1") from `primer3_core -about`; cached per binary.
// A failure is not cached. Rejects 503 PRIMER3_UNAVAILABLE.
function version(opts) {
  opts = opts || {};
  const bin = opts.bin || cfg().primer3_core;
  if (!versionCache.has(bin)) {
    const p = about(bin, opts.timeoutMs || 2000, opts.tmpDir).catch(function (e) {
      versionCache.delete(bin);
      throw e;
    });
    versionCache.set(bin, p);
  }
  return versionCache.get(bin);
}

function _clearVersionCache() {
  versionCache.clear();
}

// makeWorkDir, removeWorkDir and SPAWN_ENV are shared with thermo.js (ntthal runs the same way).
module.exports = { run, version, makeWorkDir, removeWorkDir, SPAWN_ENV, _clearVersionCache };
