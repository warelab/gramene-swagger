'use strict';

// Repeat avoidance at design time (spec §B.14, decision: real soft-mask where present, else a
// megablast copy-number mask).
//
//   template mode sequence with lowercase letters  -> the user's own mask ('user_lowercase')
//   assembly repeat_masking === 'soft_masked'      -> lowercase runs of dna_sm ('softmask'); transcript
//                                                     templates splice dna_sm exactly like the cDNA
//   repeat_mask.enabled and a dna BLAST DB         -> megablast depth mask ('blast_depth', warning
//                                                     BLAST_DEPTH_MASK): bases covered by >= min_depth
//                                                     HSPs of >= min_hsp_len at >= perc_identity
//   otherwise                                      -> no mask, warning NO_REPEAT_MASK
// A megablast failure or timeout gives warning REPEAT_MASK_FAILED and an unmasked design; the
// design deadline (opts.signal) still aborts with its own error.

const { spawn } = require('child_process');
const crypto = require('crypto');
const fs = require('fs');
const os = require('os');
const path = require('path');
const coords = require('./coords');
const { PrimerHttpError, redactPaths } = require('./errors');

const NICE_BIN = '/usr/bin/nice';
const NICE_LEVEL = 10;
// NCBI_DONT_USE_NCBIRC: BLAST+ would otherwise read .ncbirc from its cwd, $HOME, /etc and its own directory.
const SPAWN_ENV = Object.freeze({ PATH: '/usr/bin:/bin', NCBI_DONT_USE_NCBIRC: '1' });
const OUTFMT = '6 qseqid qstart qend length';
const MAX_TARGET_SEQS = 500;
const MAX_HSPS = 200;
// Genomic templates are cut into cfg.repeat_mask.chunk windows that overlap by this much, so a
// repeat crossing a chunk border is still seen whole by one query.
const CHUNK_OVERLAP = 1000;
const MOSTLY_REPEAT_FRACTION = 0.8;
const STDERR_TAIL_CHARS = 2000;
const MAX_LINE_CHARS = 64 * 1024;
const MASK_MODES = Object.freeze(['n_mask', 'three_prime']);

const MESSAGES = Object.freeze({
  BLAST_DEPTH_MASK: 'repeats were estimated from genome copy number (megablast depth); multi-copy gene families are masked as repeats',
  REPEAT_MASK_FAILED: 'the repeat mask could not be computed; primers were designed without repeat masking',
  NO_REPEAT_MASK: 'no repeat mask is available for this template; primers were designed without repeat masking',
  MOSTLY_REPEAT: 'more than 80% of the template is masked as repeat'
});

function warning(code, message) {
  return { code: code, message: message || MESSAGES[code] };
}

function logSafe(log, level, msg) {
  try {
    const l = log || console;
    (l[level] || l.error || l.log).call(l, msg);
  } catch (e) { /* logging must never throw */ }
}

function coded(code, message) {
  const e = new Error(message);
  e.code = code;
  return e;
}

function abortReason(signal) {
  const r = signal && signal.reason;
  if (r instanceof PrimerHttpError) return r;
  return new PrimerHttpError(504, 'DEADLINE_EXCEEDED', 'request deadline exceeded', {});
}

// ---- small LRU -----------------------------------------------------------------------------------

function createLru(maxEntries) {
  const map = new Map();
  const max = Math.max(1, maxEntries | 0);
  return {
    get: function (key) {
      if (!map.has(key)) return undefined;
      const v = map.get(key);
      map.delete(key);
      map.set(key, v);
      return v;
    },
    set: function (key, value) {
      map.delete(key);
      map.set(key, value);
      while (map.size > max) map.delete(map.keys().next().value);
    },
    has: function (key) { return map.has(key); },
    clear: function () { map.clear(); },
    get size() { return map.size; }
  };
}

let sharedCache = null;

function cacheFor(deps, cfg) {
  if (deps.cache) return deps.cache;
  if (!sharedCache) sharedCache = createLru(cfg.repeat_mask.cache_entries);
  return sharedCache;
}

// ---- masking primitives ---------------------------------------------------------------------------

// Primer3 template for a mask: n_mask -> masked bases become N; three_prime -> masked bases lowercase.
// Unmasked bases are always uppercase.
function applyMask(seq, runs, mode) {
  const upper = String(seq).toUpperCase();
  const merged = coords.mergeRuns(runs || [], { length: upper.length });
  if (merged.length === 0) return upper;
  const parts = [];
  let pos = 0;
  merged.forEach(function (r) {
    const s = r[0] - 1;
    parts.push(upper.slice(pos, s));
    parts.push(mode === 'three_prime' ? upper.slice(s, s + r[1]).toLowerCase() : 'N'.repeat(r[1]));
    pos = s + r[1];
  });
  parts.push(upper.slice(pos));
  return parts.join('');
}

// Boulder tags a masking mode adds to the record.
function maskTags(mode) {
  return mode === 'three_prime' ? { PRIMER_LOWERCASE_MASKING: 1 } : { PRIMER_MAX_NS_ACCEPTED: 0 };
}

// ---- megablast depth ------------------------------------------------------------------------------

function megablastArgs(db, rm) {
  if (typeof db !== 'string' || db === '' || /\s/.test(db) || db.charAt(0) === '-') {
    throw new TypeError('repeat_mask: BLAST DB path must be a non-empty path without whitespace');
  }
  return ['-task', 'megablast', '-db', db.replace(/\.(nal|nin)$/, ''), '-query', '-', '-dust', 'no', '-soft_masking', 'false',
    '-evalue', String(rm.evalue), '-perc_identity', String(rm.perc_identity),
    '-max_target_seqs', String(MAX_TARGET_SEQS), '-max_hsps', String(MAX_HSPS),
    '-num_threads', String(rm.threads), '-outfmt', OUTFMT];
}

// "q3\t10\t250\t241" -> {query: 3, qstart: 10, qend: 250, length: 241} | null
function parseHspLine(line) {
  const f = String(line).split('\t');
  if (f.length < 4) return null;
  const m = /^q(\d+)$/.exec(f[0]);
  if (!m) return null;
  const qstart = Number(f[1]);
  const qend = Number(f[2]);
  const length = Number(f[3]);
  if (!Number.isSafeInteger(qstart) || !Number.isSafeInteger(qend) || !Number.isSafeInteger(length)) return null;
  return { query: Number(m[1]), qstart: qstart, qend: qend, length: length };
}

// Genomic (or pasted) template -> overlapping chunk queries {seq, offset, lo, hi}:
// query position p (lo..hi) is template position offset + p.
function chunkQueries(seq, chunk, overlap) {
  const upper = String(seq).toUpperCase();
  const n = upper.length;
  const size = Math.max(1, chunk | 0);
  const step = Math.max(1, size - Math.max(0, Math.min(overlap | 0, size - 1)));
  const out = [];
  for (let s = 0; s < n; s += step) {
    const e = Math.min(n, s + size);
    out.push({ seq: upper.slice(s, e), offset: s, lo: 1, hi: e - s });
    if (e >= n) break;
  }
  return out;
}

// Accumulates HSP coverage per query and projects depth >= minDepth onto template runs.
function createDepthCounter(queries, opts) {
  const diffs = queries.map(function (q) { return new Int32Array(q.seq.length + 2); });
  return {
    add: function (h) {
      if (!h || h.query < 0 || h.query >= queries.length || h.length < opts.minLen) return;
      const n = queries[h.query].seq.length;
      const s = Math.max(1, Math.min(h.qstart, h.qend));
      const e = Math.min(n, Math.max(h.qstart, h.qend));
      if (e < s) return;
      diffs[h.query][s]++;
      diffs[h.query][e + 1]--;
    },
    runs: function (templateLength) {
      const flags = new Uint8Array(templateLength);
      queries.forEach(function (q, i) {
        const d = diffs[i];
        let depth = 0;
        for (let p = 1; p <= q.seq.length; p++) {
          depth += d[p];
          if (depth >= opts.minDepth && p >= q.lo && p <= q.hi) {
            const t = q.offset + p;
            if (t >= 1 && t <= templateLength) flags[t - 1] = 1;
          }
        }
      });
      return coords.runsFromFlags(flags);
    }
  };
}

// Pure: HSPs [{query, qstart, qend, length}] over queries -> template mask runs.
function depthRuns(hsps, queries, templateLength, opts) {
  const counter = createDepthCounter(queries, opts);
  (hsps || []).forEach(function (h) { counter.add(h); });
  return counter.runs(templateLength);
}

// Spawn `nice -n 10 blastn ...` with FASTA on stdin; onLine(line) for each stdout line.
// Resolves when blastn exits 0. Rejects: code BLAST_TIMEOUT (timer), BLAST_FAILED (spawn error,
// non-zero exit, overlong line, no working directory), or the abort reason when opts.signal aborts.
// The child is SIGKILLed. Its cwd is a private mkdtemp directory (mode 0700) under opts.tmpDir (absolute,
// else os.tmpdir()), removed afterwards: BLAST+ reads ./.ncbirc, so it never runs in a shared directory.
async function runMegablast(opts) {
  const signal = opts.signal;
  if (signal && signal.aborted) throw abortReason(signal);
  if (!(opts.timeoutMs > 0)) throw coded('BLAST_TIMEOUT', 'no time left for the repeat mask');
  const base = typeof opts.tmpDir === 'string' && path.isAbsolute(opts.tmpDir) ? opts.tmpDir : os.tmpdir();
  let cwd;
  try {
    cwd = await fs.promises.mkdtemp(path.join(base, 'primers-megablast-'));
  } catch (e) {
    throw coded('BLAST_FAILED', 'could not create a working directory for blastn (' + (e && e.code) + ')');
  }
  try {
    return await spawnMegablast(opts, cwd);
  } finally {
    await fs.promises.rm(cwd, { recursive: true, force: true }).catch(function () { /* best effort */ });
  }
}

function spawnMegablast(opts, cwd) {
  return new Promise(function (resolve, reject) {
    const signal = opts.signal;
    // The signal may have aborted while the working directory was created.
    if (signal && signal.aborted) return reject(abortReason(signal));
    let child;
    try {
      child = spawn(opts.niceBin || NICE_BIN, ['-n', String(opts.nice === undefined ? NICE_LEVEL : opts.nice), opts.bin].concat(opts.args),
        { cwd: cwd, env: SPAWN_ENV, stdio: ['pipe', 'pipe', 'pipe'], windowsHide: true });
    } catch (e) {
      return reject(coded('BLAST_FAILED', 'blastn could not be started (' + (e && e.code) + ')'));
    }
    let settled = false;
    let buf = '';
    let errTail = '';
    let timer = null;
    let onAbort = null;

    function finish(err) {
      if (settled) return;
      settled = true;
      if (timer) clearTimeout(timer);
      if (signal && onAbort) signal.removeEventListener('abort', onAbort);
      if (err) {
        try { child.kill('SIGKILL'); } catch (e) { /* already gone */ }
        try { child.stdout.destroy(); child.stderr.destroy(); } catch (e) { /* ignore */ }
        reject(err);
      } else {
        resolve();
      }
    }
    function emit(line) {
      try {
        opts.onLine(line);
      } catch (e) {
        finish(e);
      }
    }

    child.on('error', function (err) {
      finish(coded('BLAST_FAILED', 'blastn could not be started (' + (err && err.code) + ')'));
    });
    child.stdout.setEncoding('utf8');
    child.stdout.on('data', function (chunk) {
      if (settled) return;
      buf += chunk;
      let start = 0;
      let nl;
      while (!settled && (nl = buf.indexOf('\n', start)) >= 0) {
        emit(buf.slice(start, nl));
        start = nl + 1;
      }
      buf = buf.slice(start);
      if (buf.length > MAX_LINE_CHARS) finish(coded('BLAST_FAILED', 'blastn output line too long'));
    });
    child.stderr.on('data', function (d) {
      errTail = (errTail + d.toString('utf8')).slice(-STDERR_TAIL_CHARS);
    });
    child.stdin.on('error', function () { /* EPIPE when blastn exits early */ });
    child.on('close', function (code, sig) {
      if (settled) return;
      if (buf) emit(buf);
      buf = '';
      if (code !== 0 || sig) {
        return finish(coded('BLAST_FAILED', 'blastn exited with ' + (sig ? 'signal ' + sig : 'code ' + code) + ': ' +
          redactPaths(errTail.slice(-500)).trim()));
      }
      finish(null);
    });
    timer = setTimeout(function () {
      finish(coded('BLAST_TIMEOUT', 'megablast did not finish within ' + opts.timeoutMs + ' ms'));
    }, opts.timeoutMs);
    if (signal) {
      onAbort = function () { finish(abortReason(signal)); };
      signal.addEventListener('abort', onAbort, { once: true });
    }
    child.stdin.end(opts.input);
  });
}

// Exon +- exon_pad (genomic) queries for a transcript template, projected onto the cDNA.
async function exonQueries(template, fastaPath, cfg, sequence) {
  const loc = template.location;
  const pad = Math.max(0, cfg.repeat_mask.exon_pad | 0);
  const regionLen = await sequence.regionLength(fastaPath, loc.region);
  if (!Number.isSafeInteger(regionLen)) throw coded('BLAST_FAILED', 'region ' + loc.region + ' is not in the assembly');
  const windows = template.segments.map(function (s) {
    return { start: Math.max(1, s.g_start - pad), end: Math.min(regionLen, s.g_end + pad) };
  });
  const seqs = await require('./template').fetchWindows(sequence, fastaPath, loc.region, loc.strand, windows,
    { maxLength: cfg.design.max_fetch_length });
  return template.segments.map(function (s, i) {
    const w = windows[i];
    const padUp = loc.strand === 1 ? s.g_start - w.start : w.end - s.g_end;
    const n = s.t_end - s.t_start + 1;
    return { seq: seqs[i].toUpperCase(), offset: s.t_start - 1 - padUp, lo: padUp + 1, hi: padUp + n };
  });
}

function sha1(s) {
  return crypto.createHash('sha1').update(s).digest('hex');
}

function blastCacheKey(template, asm, rm) {
  const parts = [asm.system_name, asm.fingerprint || '', asm.blastdb.dna, template.mode, rm.min_hsp_len, rm.perc_identity,
    rm.evalue, rm.min_depth, rm.chunk, rm.exon_pad];
  if (template.mode === 'sequence') {
    parts.push(sha1(template.seq.toUpperCase()));
  } else {
    parts.push(template.region, template.start, template.end, template.strand);
    if (template.segments) parts.push(template.segments.map(function (s) { return s.g_start + '-' + s.g_end; }).join(','));
  }
  return parts.join('|');
}

async function blastDepthMask(template, asm, ctx) {
  const cfg = ctx.cfg;
  const rm = cfg.repeat_mask;
  const queries = template.mode === 'transcript'
    ? await exonQueries(template, asm.fasta.dna, cfg, ctx.sequence)
    : chunkQueries(template.seq, rm.chunk, CHUNK_OVERLAP);
  const counter = createDepthCounter(queries, { minLen: rm.min_hsp_len, minDepth: rm.min_depth });
  const input = queries.map(function (q, i) { return '>q' + i + '\n' + q.seq + '\n'; }).join('');
  const remaining = ctx.deadline ? ctx.deadline - ctx.now() : Infinity;
  const timeoutMs = Math.floor(Math.min(rm.timeout_ms, remaining));
  const started = ctx.now();
  await (ctx.runBlast || runMegablast)({
    bin: cfg.blastn,
    args: megablastArgs(asm.blastdb.dna, rm),
    input: input,
    timeoutMs: timeoutMs,
    signal: ctx.signal,
    tmpDir: cfg.tmp_dir,
    nice: NICE_LEVEL,
    onLine: function (line) { counter.add(parseHspLine(line)); }
  });
  const runs = counter.runs(template.length);
  logSafe(ctx.log, 'info', 'primers: megablast depth mask ' + asm.system_name + ' ' + queries.length + ' queries, ' +
    coords.maskedBases(runs) + ' of ' + template.length + ' bp masked in ' + (ctx.now() - started) + ' ms');
  return runs;
}

// ---- soft mask -----------------------------------------------------------------------------------

async function softMaskRuns(template, asm, ctx) {
  const cfg = ctx.cfg;
  const smPath = asm.fasta.dna_sm;
  let sm;
  if (template.mode === 'transcript') {
    const parts = await require('./template').fetchWindows(ctx.sequence, smPath, template.location.region, template.location.strand,
      template.segments.map(function (s) { return { start: s.g_start, end: s.g_end }; }),
      { maxLength: cfg.design.max_fetch_length });
    sm = parts.join('');
  } else {
    sm = await ctx.sequence.fetch(smPath, template.region, template.start, template.end, template.strand,
      { maxLength: cfg.design.max_fetch_length });
  }
  if (sm.length !== template.length || sm.toUpperCase() !== template.seq.toUpperCase()) {
    throw coded('SOFTMASK_MISMATCH', 'the soft-masked FASTA does not match the template sequence');
  }
  return coords.lowercaseRuns(sm);
}

// ---- entry point ---------------------------------------------------------------------------------

// computeMask(template, opts, deps) -> {runs, source, warnings}
async function computeMask(template, opts, deps) {
  const cfg = deps.cfg;
  const log = deps.log;
  const asm = template.resolved;
  const signal = opts.signal;
  const aborted = function () { return !!(signal && signal.aborted); };

  if (template.mode === 'sequence') {
    const lower = coords.lowercaseRuns(template.seq);
    if (lower.length > 0) return { runs: lower, source: 'user_lowercase', warnings: [] };
  } else if (asm && asm.repeat_masking === 'soft_masked' && asm.fasta && asm.fasta.dna_sm) {
    try {
      return { runs: await softMaskRuns(template, asm, deps), source: 'softmask', warnings: [] };
    } catch (err) {
      if (aborted()) throw abortReason(signal);
      logSafe(log, 'error', 'primers: soft mask of ' + asm.system_name + ' failed: ' + redactPaths(String(err && err.message)));
      return { runs: [], source: null, warnings: [warning('REPEAT_MASK_FAILED')] };
    }
  }

  if (cfg.repeat_mask.enabled && asm && asm.blastdb && asm.blastdb.dna && asm.fasta && asm.fasta.dna) {
    const cache = cacheFor(deps, cfg);
    const key = blastCacheKey(template, asm, cfg.repeat_mask);
    const hit = cache.get(key);
    if (hit) return { runs: hit.map(function (r) { return r.slice(); }), source: 'blast_depth', warnings: [warning('BLAST_DEPTH_MASK')] };
    try {
      const runs = await blastDepthMask(template, asm, {
        cfg: cfg, log: log, sequence: deps.sequence, runBlast: deps.runBlast, signal: signal,
        deadline: opts.deadline, now: deps.now
      });
      cache.set(key, runs.map(function (r) { return r.slice(); }));
      return { runs: runs, source: 'blast_depth', warnings: [warning('BLAST_DEPTH_MASK')] };
    } catch (err) {
      if (aborted()) throw abortReason(signal);
      if (err instanceof PrimerHttpError && err.status === 504) throw err;
      logSafe(log, 'error', 'primers: megablast depth mask of ' + asm.system_name + ' failed: ' +
        (err && err.code ? err.code + ' ' : '') + redactPaths(String(err && err.message)));
      return { runs: [], source: null, warnings: [warning('REPEAT_MASK_FAILED')] };
    }
  }
  return { runs: [], source: null, warnings: [warning('NO_REPEAT_MASK')] };
}

// repeatMask(template, opts, deps) -> {mask, mask_source, masked, masked_fraction, seq, mode, tags, warnings}
//   template: from template.buildTemplate (uses mode, seq, length, region/start/end/strand, location,
//             segments, resolved)
//   opts: {mode: 'n_mask' (default) | 'three_prime', signal (design deadline), deadline (epoch ms)}
//   deps: {cfg, log, sequence {fetch, regionLength}, runBlast(opts) (replaces the megablast spawn),
//          cache (LRU with get/set), now}
//   mask: merged [start, length] template runs; seq: the Primer3 template with the mask applied
//   (uppercase otherwise); tags: Boulder tags for the mode.
async function repeatMask(template, opts, deps) {
  opts = opts || {};
  deps = deps || {};
  const mode = opts.mode === undefined || opts.mode === null ? 'n_mask' : opts.mode;
  if (MASK_MODES.indexOf(mode) < 0) {
    throw new PrimerHttpError(400, 'INVALID_REQUEST', 'repeat_mask_mode must be n_mask or three_prime', { field: 'repeat_mask_mode' });
  }
  const ctx = {
    cfg: deps.cfg || require('./config').get(),
    log: deps.log || console,
    sequence: deps.sequence || require('./sequence'),
    runBlast: deps.runBlast,
    cache: deps.cache,
    now: deps.now || Date.now
  };
  if (opts.signal && opts.signal.aborted) throw abortReason(opts.signal);
  const found = await computeMask(template, opts, ctx);
  const runs = coords.mergeRuns(found.runs, { length: template.length });
  const fraction = coords.maskedFraction(runs, template.length);
  const warnings = found.warnings.slice();
  if (fraction > MOSTLY_REPEAT_FRACTION) warnings.push(warning('MOSTLY_REPEAT'));
  return {
    mask: runs,
    mask_source: found.source,
    masked: runs.length > 0,
    masked_fraction: Math.round(fraction * 10000) / 10000,
    seq: applyMask(template.seq, runs, mode),
    mode: mode,
    tags: maskTags(mode),
    warnings: warnings
  };
}

function _clearCache() {
  if (sharedCache) sharedCache.clear();
  sharedCache = null;
}

function _cacheSize() {
  return sharedCache ? sharedCache.size : 0;
}

module.exports = {
  repeatMask,
  computeMask,
  applyMask,
  maskTags,
  megablastArgs,
  parseHspLine,
  chunkQueries,
  depthRuns,
  runMegablast,
  createLru,
  blastCacheKey,
  MASK_MODES,
  MESSAGES,
  CHUNK_OVERLAP,
  MOSTLY_REPEAT_FRACTION,
  NICE_BIN,
  NICE_LEVEL,
  SPAWN_ENV,
  _clearCache,
  _cacheSize
};
