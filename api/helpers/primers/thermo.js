'use strict';

// ntthal (Primer3 2.6.1) for the genotyping design (spec §4.10, §4.11, §4.14): tailed hairpins, tailed self and cross
// dimers, and the duplex Tm of a deliberate-mismatch primer against its own-allele footprint. Nominal Tm values never
// come from here (Primer3 itself supplies them), and oligotm is not used.
//
// createThermo({cfg, spawn, signal, deadline, salts, log}) -> {hairpin, selfAny, selfEnd, cross, duplex, run, calls}
//   cfg: the primers config (ntthal, tmp_dir, genotyping.thermo_concurrency, genotyping.thermo_timeout_ms)
//   spawn(bin, argv, {timeoutMs, signal, log, tmpDir}) -> Promise<{stdout}> (default runBinary; tests inject recordings)
//   signal: abort signal (the design deadline); deadline: design.createDeadline, whose remaining() caps each call
//   salts: {salt_monovalent, salt_divalent, dntp_conc, dna_conc} (the effective design params; Primer3 defaults)
// One instance serves one design request. Each call validates its sequences against ^[ACGT]{1,60}$ before anything is
// spawned (a TypeError otherwise), is memoized per argv (salts, mode, s1, s2), and waits for one of
// thermo_concurrency (4) slots. `calls` counts the distinct calls requested so far: the design budget's thermo_calls.
// Results are exact decimals as ntthal prints them ('42.616970'); "No secondary structure" and negative values are '0'.
//
// ntthal runs like primer3_core: no shell, an argument array, PATH=/usr/bin:/bin, a private mkdtemp working directory
// removed afterwards, a timeout and the abort signal (SIGKILL). Its own defaults (-dv 0 -n 0) are not Primer3's, so the
// salts are always passed: -mv 50 -dv 1.5 -n 0.6 -d 50 -t 37 -r.

const { spawn: childSpawn } = require('child_process');
const primer3 = require('./primer3');
const decimal = require('./genotyping/decimal');
const { PrimerHttpError, redactPaths } = require('./errors');

const SEQ_RE = /^[ACGT]{1,60}$/;
const MODES = Object.freeze(['ANY', 'END1', 'HAIRPIN']);
const TEMPERATURE_C = 37;
// PRIMER_SALT_MONOVALENT, PRIMER_SALT_DIVALENT, PRIMER_DNTP_CONC and PRIMER_DNA_CONC defaults of Primer3 2.6.1.
const PRIMER3_SALTS = Object.freeze({ salt_monovalent: 50, salt_divalent: 1.5, dntp_conc: 0.6, dna_conc: 50 });
const DEFAULT_CONCURRENCY = 4;
const DEFAULT_TIMEOUT_MS = 5000;
const STDOUT_CAP = 4096;
const LOG_STDERR_CHARS = 500;
const BINARY = 'ntthal';
const UNAVAILABLE = new Set(['ENOENT', 'EACCES', 'EPERM', 'ENOEXEC', 'ENOTDIR']);
const TRANSIENT = new Set(['EAGAIN', 'EMFILE', 'ENFILE', 'ENOMEM']);
const NUMBER_RE = /^-?\d+(?:\.\d+)?$/;

function logError(log, msg) {
  try {
    const l = log || console;
    (l.error || l.log).call(l, msg);
  } catch (e) { /* logging must never throw */ }
}

function unavailableError(retryAfter) {
  return new PrimerHttpError(503, 'THERMO_UNAVAILABLE', 'ntthal is not available on this server', { binary: BINARY, retry_after_s: retryAfter });
}

function failedError() {
  return new PrimerHttpError(500, 'THERMO_FAILED', 'ntthal failed', { binary: BINARY });
}

function deadlineError(timeoutMs) {
  return new PrimerHttpError(504, 'DEADLINE_EXCEEDED', 'ntthal did not finish within the time limit', { timeout_ms: timeoutMs });
}

function abortError(signal) {
  const r = signal && signal.reason;
  if (r instanceof PrimerHttpError) return r;
  return new PrimerHttpError(504, 'DEADLINE_EXCEEDED', 'request deadline exceeded', {});
}

function spawnError(err, log) {
  const code = err && err.code;
  if (UNAVAILABLE.has(code)) {
    logError(log, 'ntthal unavailable: ' + code);
    return unavailableError(60);
  }
  if (TRANSIENT.has(code)) {
    logError(log, 'ntthal could not start: ' + code);
    return unavailableError(5);
  }
  logError(log, 'ntthal spawn failed: ' + redactPaths(String(err && err.message)));
  return failedError();
}

// The effective salts: each finite number of `params`, else Primer3's default.
function saltsOf(params) {
  const p = params || {};
  const out = {};
  Object.keys(PRIMER3_SALTS).forEach(function (k) {
    out[k] = typeof p[k] === 'number' && Number.isFinite(p[k]) ? p[k] : PRIMER3_SALTS[k];
  });
  return out;
}

// argv(salts, mode, s1, s2) -> the ntthal argument array; its ' '-joined form is the memo key and the recording key.
function argv(salts, mode, s1, s2) {
  const s = saltsOf(salts);
  const out = ['-mv', String(s.salt_monovalent), '-dv', String(s.salt_divalent), '-n', String(s.dntp_conc), '-d', String(s.dna_conc),
    '-t', String(TEMPERATURE_C), '-r', '-a', mode, '-s1', s1];
  if (s2 !== undefined && s2 !== null) out.push('-s2', s2);
  return out;
}

// parseOutput(stdout) -> exact decimal string, '0' for "No secondary structure" or a negative value, or null when the
// output is not one number.
function parseOutput(stdout) {
  const text = String(stdout === undefined || stdout === null ? '' : stdout).trim();
  if (/^No secondary structure/i.test(text)) return '0';
  if (!NUMBER_RE.test(text)) return null;
  return text.charAt(0) === '-' ? '0' : text;
}

function checkSequence(s, name) {
  if (typeof s !== 'string' || !SEQ_RE.test(s)) {
    const e = new TypeError('thermo: ' + name + ' must match ' + SEQ_RE.source + ' (got ' + String(s).slice(0, 70) + ')');
    e.code = 'THERMO_INVALID_SEQUENCE';
    throw e;
  }
}

// runBinary(bin, args, {timeoutMs, signal, log, tmpDir}) -> Promise<{stdout}>: one ntthal process in a private working
// directory. Rejects 503 THERMO_UNAVAILABLE (missing binary), 500 THERMO_FAILED (non-zero exit, signal, output over the
// cap, no working directory), 504 DEADLINE_EXCEEDED (timeout or abort; the child is SIGKILLed).
function runBinary(bin, args, opts) {
  opts = opts || {};
  const log = opts.log || console;
  if (opts.signal && opts.signal.aborted) return Promise.reject(abortError(opts.signal));
  return primer3.makeWorkDir(opts.tmpDir).then(function (cwd) {
    return spawnIn(bin, args, cwd, opts).finally(function () { return primer3.removeWorkDir(cwd); });
  }, function (err) {
    logError(log, 'ntthal working directory could not be created: ' + (err && err.code));
    throw failedError();
  });
}

function spawnIn(bin, args, cwd, opts) {
  const log = opts.log || console;
  const signal = opts.signal;
  const timeoutMs = opts.timeoutMs > 0 ? opts.timeoutMs : DEFAULT_TIMEOUT_MS;
  return new Promise(function (resolve, reject) {
    if (signal && signal.aborted) return reject(abortError(signal));
    let child;
    try {
      child = childSpawn(bin, args, { cwd: cwd, env: primer3.SPAWN_ENV, stdio: ['ignore', 'pipe', 'pipe'], windowsHide: true });
    } catch (e) {
      return reject(spawnError(e, log));
    }
    let out = '';
    let errTail = '';
    let settled = false;
    let failure = null;
    let onAbort = null;
    const timer = setTimeout(function () { kill(deadlineError(timeoutMs)); }, timeoutMs);

    function finish(err, value) {
      if (settled) return;
      settled = true;
      clearTimeout(timer);
      if (signal && onAbort) signal.removeEventListener('abort', onAbort);
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
      out += chunk.toString('utf8');
      if (out.length > STDOUT_CAP) {
        logError(log, 'ntthal stdout exceeded ' + STDOUT_CAP + ' bytes; killed');
        kill(failedError());
      }
    });
    child.stderr.on('data', function (chunk) {
      errTail = (errTail + chunk.toString('utf8')).slice(-LOG_STDERR_CHARS);
    });
    child.on('exit', function () {
      if (failure) finish(failure);
    });
    child.on('close', function (code, sig) {
      if (failure) return finish(failure);
      if (code !== 0 || sig) {
        logError(log, 'ntthal exited with ' + (sig ? 'signal ' + sig : 'code ' + code) + ': ' + redactPaths(errTail));
        return finish(failedError());
      }
      finish(null, { stdout: out });
    });
    if (signal) {
      onAbort = function () { kill(abortError(signal)); };
      signal.addEventListener('abort', onAbort, { once: true });
    }
  });
}

function posInt(v) {
  return Number.isSafeInteger(v) && v > 0 ? v : null;
}

function createThermo(opts) {
  opts = opts || {};
  const cfg = opts.cfg || require('./config').get();
  const g = cfg.genotyping || {};
  const bin = cfg.ntthal;
  const concurrency = posInt(g.thermo_concurrency) || DEFAULT_CONCURRENCY;
  const timeoutMs = posInt(g.thermo_timeout_ms) || DEFAULT_TIMEOUT_MS;
  const exec = typeof opts.spawn === 'function' ? opts.spawn : runBinary;
  const salts = saltsOf(opts.salts);
  const signal = opts.signal;
  const deadline = opts.deadline;
  const log = opts.log || console;
  const memo = new Map();
  const waiting = [];
  let active = 0;

  function acquire() {
    if (active < concurrency) {
      active++;
      return Promise.resolve();
    }
    return new Promise(function (resolve) { waiting.push(resolve); });
  }
  function release() {
    const next = waiting.shift();
    if (next) next();
    else active--;
  }

  function run(mode, s1, s2) {
    if (MODES.indexOf(mode) < 0) throw new TypeError('thermo: mode must be one of ' + MODES.join(', '));
    checkSequence(s1, 's1');
    if (mode === 'HAIRPIN') {
      if (s2 !== undefined && s2 !== null) throw new TypeError('thermo: HAIRPIN takes one sequence');
    } else {
      checkSequence(s2, 's2');
    }
    const args = argv(salts, mode, s1, s2);
    const key = args.join(' ');
    if (memo.has(key)) return memo.get(key);
    const p = acquire().then(function () {
      if (signal && signal.aborted) throw abortError(signal);
      let t = timeoutMs;
      if (deadline && typeof deadline.remaining === 'function') {
        const left = Math.floor(deadline.remaining());
        if (left <= 0) throw deadlineError(timeoutMs);
        t = Math.min(t, left);
      }
      return exec(bin, args, { timeoutMs: t, signal: signal, log: log, tmpDir: cfg.tmp_dir });
    }).then(function (r) {
      const value = parseOutput(r && r.stdout);
      if (value === null) {
        logError(log, 'ntthal printed an unparseable result: ' + JSON.stringify(String(r && r.stdout).slice(0, 80)));
        throw failedError();
      }
      return value;
    }).finally(release);
    p.catch(function () { /* observed by the callers */ });
    memo.set(key, p);
    return p;
  }

  return {
    run: run,
    hairpin: function (s) { return run('HAIRPIN', s); },
    selfAny: function (s) { return run('ANY', s, s); },
    selfEnd: function (s) { return run('END1', s, s); },
    // Tailed cross dimer of two oligos (§4.14): ANY once, END1 in both directions; end = the larger.
    cross: function (a, b) {
      return Promise.all([run('ANY', a, b), run('END1', a, b), run('END1', b, a)]).then(function (r) {
        return { any: r[0], end: decimal.micro(r[1]) >= decimal.micro(r[2]) ? r[1] : r[2] };
      });
    },
    // Duplex Tm of a primer against its template strand, both 5'->3' (§4.10: target = revcomp(matched_seq)).
    duplex: function (primer, target) { return run('ANY', primer, target); },
    get calls() { return memo.size; },
    salts: salts
  };
}

module.exports = {
  createThermo,
  runBinary,
  argv,
  parseOutput,
  saltsOf,
  PRIMER3_SALTS,
  MODES,
  SEQ_RE,
  BINARY
};
