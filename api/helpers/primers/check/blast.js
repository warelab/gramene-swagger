'use strict';

// BLAST argument arrays, query FASTA, streaming outfmt-6 parser and a local
// niced line-streaming spawner (spec §B.3).

const { spawn } = require('child_process');

const GENOME_OUTFMT = '6 qseqid sseqid qlen qstart qend sstart send sstrand mismatch';
const CDNA_OUTFMT = GENOME_OUTFMT + ' qseq sseq';
const MIN_MAX_TARGET_SEQS = 5000;
const DEFAULT_TIMEOUT_MS = 600000;
const STDERR_TAIL_BYTES = 2048;
const KILL_GRACE_MS = 3000;
const MAX_LINE_CHARS = 1 << 20;
const NICE_BIN = '/usr/bin/nice';
const DEFAULT_NICE = 10;
// NCBI_DONT_USE_NCBIRC: BLAST+ must not read a .ncbirc from the cwd, $HOME, /etc or the binary directory
// (a planted /tmp/.ncbirc could otherwise break or reconfigure every run). Callers also pass a private cwd.
const CHILD_ENV = Object.freeze({ PATH: '/usr/bin:/bin', NCBI_DONT_USE_NCBIRC: '1' });

function codedError(code, message, extra) {
  const err = new Error(message);
  err.code = code;
  if (extra) Object.assign(err, extra);
  return err;
}

// Exactly the §B.3 argument array. `maxTargetSeqs` is the sequence count of the
// target DB (fai line count for genomes, `blastdbcmd -info` for cDNA); the
// argument sent is max(5000, maxTargetSeqs). `db` is the DB path without
// extension (a trailing .nal/.nin is stripped so multi-volume DBs resolve).
function buildArgs(opts) {
  const o = opts || {};
  if (o.target !== 'genome' && o.target !== 'cdna') {
    throw new TypeError("target must be 'genome' or 'cdna'");
  }
  if (!Number.isInteger(o.wordSize) || o.wordSize < 4 || o.wordSize > 64) {
    throw new TypeError('wordSize must be an integer in 4..64');
  }
  if (!Number.isInteger(o.threads) || o.threads < 1 || o.threads > 64) {
    throw new TypeError('threads must be an integer in 1..64');
  }
  if (typeof o.db !== 'string' || o.db === '' || /\s/.test(o.db) || o.db[0] === '-') {
    throw new TypeError('db must be a non-empty path without whitespace');
  }
  const db = o.db.replace(/\.(nal|nin)$/, '');
  const mts = Math.max(MIN_MAX_TARGET_SEQS, Number.isInteger(o.maxTargetSeqs) ? o.maxTargetSeqs : 0);
  return ['-task', 'blastn-short', '-reward', '1', '-penalty', '-1', '-word_size', String(o.wordSize), '-ungapped',
    '-evalue', '30000', '-searchsp', '15000000000', '-dust', 'no', '-soft_masking', 'false',
    '-max_target_seqs', String(mts), '-max_hsps', '100000',
    '-num_threads', String(o.threads), '-db', db, '-query', '-',
    '-outfmt', o.target === 'cdna' ? CDNA_OUTFMT : GENOME_OUTFMT];
}

// Short description for results.engine.reference / .pangenome.
function engineDescription(wordSize) {
  return 'blastn-short r1 p-1 ws' + wordSize + ' ungapped e30000 searchsp1.5e10';
}

function targetFromArgs(args) {
  const i = Array.isArray(args) ? args.indexOf('-outfmt') : -1;
  if (i < 0 || typeof args[i + 1] !== 'string') return 'genome';
  return / qseq( |$)/.test(args[i + 1]) ? 'cdna' : 'genome';
}

// Uppercase, dedupe, keep first-seen order.
function uniquePrimers(seqs) {
  const seen = new Set();
  const out = [];
  for (const s of seqs || []) {
    const u = String(s).toUpperCase();
    if (!seen.has(u)) {
      seen.add(u);
      out.push(u);
    }
  }
  return out;
}

// '>q<i>\n<SEQ>\n' per primer; i is the 0-based index into `uniquePrimers`, so
// hit.q maps straight back to the primer. Rejects non-ACGT input (no FASTA
// injection) and duplicates (the index mapping must be 1:1).
function queryFasta(primers) {
  if (!Array.isArray(primers) || primers.length === 0) {
    throw new TypeError('primers must be a non-empty array');
  }
  const seen = new Set();
  let out = '';
  primers.forEach((s, i) => {
    const u = typeof s === 'string' ? s.toUpperCase() : '';
    if (!/^[ACGT]+$/.test(u)) throw new TypeError('primer ' + i + ' is not an ACGT sequence');
    if (seen.has(u)) throw new TypeError('primer ' + i + ' is a duplicate; pass unique primers');
    seen.add(u);
    out += '>q' + i + '\n' + u + '\n';
  });
  return out;
}

function toInt(field, name, min) {
  // outfmt 6 integers are plain decimal; reject anything else.
  if (!/^\d+$/.test(field)) throw codedError('BLAST_PARSE', 'bad ' + name + ': ' + field);
  const v = Number(field);
  if (v < min) throw codedError('BLAST_PARSE', name + ' below ' + min + ': ' + field);
  return v;
}

// One outfmt-6 line → hit, or null for blank/comment lines. Throws
// err.code 'BLAST_PARSE' on malformed lines.
//   hit = { q, qseqid, sseqid, qlen, qstart, qend, sstart, send, strand: 1|-1, mismatch, qseq?, sseq? }
function parseLine(line, target) {
  if (typeof line !== 'string') return null;
  if (line.length && line.charCodeAt(line.length - 1) === 13) line = line.slice(0, -1);
  if (line === '' || line[0] === '#') return null;
  const f = line.split('\t');
  const cdna = target === 'cdna';
  const want = cdna ? 11 : 9;
  if (f.length !== want) throw codedError('BLAST_PARSE', 'expected ' + want + ' fields, got ' + f.length);
  const m = /^q(\d+)$/.exec(f[0]);
  if (!m) throw codedError('BLAST_PARSE', 'bad qseqid: ' + f[0]);
  if (f[1] === '') throw codedError('BLAST_PARSE', 'empty sseqid');
  const hit = {
    q: Number(m[1]),
    qseqid: f[0],
    sseqid: f[1],
    qlen: toInt(f[2], 'qlen', 1),
    qstart: toInt(f[3], 'qstart', 1),
    qend: toInt(f[4], 'qend', 1),
    sstart: toInt(f[5], 'sstart', 1),
    send: toInt(f[6], 'send', 1),
    strand: 0,
    mismatch: toInt(f[8], 'mismatch', 0)
  };
  if (f[7] === 'plus') hit.strand = 1;
  else if (f[7] === 'minus') hit.strand = -1;
  else throw codedError('BLAST_PARSE', 'bad sstrand: ' + f[7]);
  if (hit.qstart > hit.qend || hit.qend > hit.qlen) throw codedError('BLAST_PARSE', 'bad query range');
  if (hit.strand === 1 ? hit.sstart > hit.send : hit.sstart < hit.send) {
    throw codedError('BLAST_PARSE', 'subject range does not match strand');
  }
  if (cdna) {
    hit.qseq = f[9];
    hit.sseq = f[10];
    if (!hit.qseq || hit.qseq.length !== hit.sseq.length) throw codedError('BLAST_PARSE', 'bad qseq/sseq');
  }
  return hit;
}

// Spawn `cmd args` (under `nice -n 10` by default) and stream stdout lines to
// onLine without buffering the output. Resolves in every process outcome:
//   { code, signal, stderrTail, timedOut, aborted, elapsedMs }
// Rejects only when the process cannot be spawned, when onLine throws (the
// child is killed first) or when a single line exceeds 1 MB.
// opts: { stdin?: string|Buffer, timeoutMs?: number (0/undefined = none), onLine?: fn,
//         signal?: AbortSignal, nice?: number|false, niceBin?, env?, cwd?, killGraceMs? }
// Abort and timeout send SIGTERM, then SIGKILL after killGraceMs (3 s).
function spawnLinesLocal(cmd, args, opts) {
  const o = opts || {};
  return new Promise((resolve, reject) => {
    const started = Date.now();
    const signal = o.signal;
    if (signal && signal.aborted) {
      resolve({ code: null, signal: null, stderrTail: '', timedOut: false, aborted: true, elapsedMs: 0 });
      return;
    }
    const niceness = o.nice === undefined ? DEFAULT_NICE : o.nice;
    const useNice = niceness !== false && niceness !== null;
    const exe = useNice ? (o.niceBin || NICE_BIN) : cmd;
    const argv = useNice ? ['-n', String(niceness), cmd].concat(args || []) : (args || []).slice();
    const graceMs = Number.isInteger(o.killGraceMs) ? o.killGraceMs : KILL_GRACE_MS;

    let child;
    try {
      child = spawn(exe, argv, { stdio: ['pipe', 'pipe', 'pipe'], env: o.env || CHILD_ENV, cwd: o.cwd });
    } catch (e) {
      reject(e);
      return;
    }

    let settled = false;
    let exited = false;
    let timedOut = false;
    let aborted = false;
    let failure = null;
    let pending = '';
    let stderrTail = '';
    let timeoutTimer = null;
    let killTimer = null;

    function terminate() {
      if (exited) return;
      try { child.kill('SIGTERM'); } catch (e) { /* already gone */ }
      if (!killTimer) {
        killTimer = setTimeout(() => {
          if (!exited) {
            try { child.kill('SIGKILL'); } catch (e) { /* already gone */ }
          }
        }, graceMs);
      }
    }

    function onAbort() {
      aborted = true;
      terminate();
    }

    function cleanup() {
      if (timeoutTimer) clearTimeout(timeoutTimer);
      if (killTimer) clearTimeout(killTimer);
      if (signal) signal.removeEventListener('abort', onAbort);
    }

    function deliver(line) {
      if (failure || !o.onLine) return;
      try {
        o.onLine(line);
      } catch (e) {
        failure = e;
        terminate();
      }
    }

    child.stdout.setEncoding('utf8');
    child.stdout.on('data', (chunk) => {
      if (failure) return;
      pending += chunk;
      let start = 0;
      let nl;
      while ((nl = pending.indexOf('\n', start)) !== -1) {
        let line = pending.slice(start, nl);
        if (line.length && line.charCodeAt(line.length - 1) === 13) line = line.slice(0, -1);
        start = nl + 1;
        deliver(line);
        if (failure) return;
      }
      pending = start === 0 ? pending : pending.slice(start);
      if (pending.length > MAX_LINE_CHARS) {
        failure = codedError('LINE_TOO_LONG', 'output line exceeds ' + MAX_LINE_CHARS + ' characters');
        pending = '';
        terminate();
      }
    });

    child.stderr.setEncoding('utf8');
    child.stderr.on('data', (chunk) => {
      stderrTail = (stderrTail + chunk).slice(-STDERR_TAIL_BYTES);
    });

    child.stdin.on('error', () => { /* EPIPE when the child exits early; the exit code tells the story */ });
    if (o.stdin != null) child.stdin.end(o.stdin);
    else child.stdin.end();

    if (signal) signal.addEventListener('abort', onAbort, { once: true });
    if (o.timeoutMs > 0) {
      timeoutTimer = setTimeout(() => {
        timedOut = true;
        terminate();
      }, o.timeoutMs);
    }

    child.on('error', (err) => {
      exited = true;
      cleanup();
      if (settled) return;
      settled = true;
      reject(err);
    });

    child.on('close', (code, sig) => {
      exited = true;
      cleanup();
      if (settled) return;
      settled = true;
      if (!failure && pending !== '') {
        const last = pending;
        pending = '';
        deliver(last);
      }
      if (failure) {
        reject(failure);
        return;
      }
      resolve({ code, signal: sig, stderrTail, timedOut, aborted, elapsedMs: Date.now() - started });
    });
  });
}

// Run one blastn over `primers` (unique uppercase list; hit.primer = primers[hit.q])
// and stream parsed hits to onHit. `spawnLines` is injected (worker ctx.spawnLines
// or spawnLinesLocal); `cwd` (a private directory) is passed through to it.
// Resolves { hits, code, signal, stderrTail, elapsedMs };
// throws err.code ABORTED | BLAST_TIMEOUT | BLAST_FAILED | BLAST_OUTPUT_INVALID
// (err.stderrTail attached). `target` defaults to the one implied by -outfmt.
async function runBlast(opts) {
  const o = opts || {};
  if (typeof o.cmd !== 'string' || !o.cmd) throw new TypeError('cmd is required');
  if (!Array.isArray(o.args)) throw new TypeError('args must be an array');
  const spawnLines = o.spawnLines || spawnLinesLocal;
  const target = o.target || targetFromArgs(o.args);
  const stdin = queryFasta(o.primers);
  const seqs = o.primers.map((s) => s.toUpperCase());
  const started = Date.now();
  let hits = 0;
  let badLines = 0;
  let firstBadLine = null;
  let firstBadReason = null;

  const res = await spawnLines(o.cmd, o.args, {
    stdin,
    timeoutMs: o.timeoutMs > 0 ? o.timeoutMs : DEFAULT_TIMEOUT_MS,
    signal: o.signal,
    cwd: o.cwd,
    onLine(line) {
      let hit;
      try {
        hit = parseLine(line, target);
      } catch (e) {
        badLines++;
        if (firstBadLine === null) {
          firstBadLine = line.slice(0, 200);
          firstBadReason = e.message;
        }
        return;
      }
      if (!hit) return;
      if (hit.q >= seqs.length) {
        badLines++;
        if (firstBadLine === null) {
          firstBadLine = line.slice(0, 200);
          firstBadReason = 'unknown query index';
        }
        return;
      }
      hit.primer = seqs[hit.q];
      hits++;
      if (o.onHit) o.onHit(hit);
    }
  });

  const extra = { stderrTail: res.stderrTail, exitCode: res.code, exitSignal: res.signal, hits };
  if (res.aborted || (res.code !== 0 && o.signal && o.signal.aborted)) {
    throw codedError('ABORTED', 'BLAST aborted', extra);
  }
  if (res.timedOut) {
    throw codedError('BLAST_TIMEOUT', 'BLAST timed out', extra);
  }
  if (res.code !== 0) {
    throw codedError('BLAST_FAILED', 'blastn exited with code ' + res.code + (res.signal ? ' (' + res.signal + ')' : ''), extra);
  }
  if (badLines > 0) {
    throw codedError('BLAST_OUTPUT_INVALID', badLines + ' unparseable BLAST output line(s): ' + firstBadReason,
      Object.assign(extra, { badLines, firstBadLine }));
  }
  return { hits, code: res.code, signal: res.signal, stderrTail: res.stderrTail, elapsedMs: Date.now() - started };
}

// `blastdbcmd -info -db <db>` text → { title, num_sequences, total_bases, volumes[] }.
function parseDbInfo(text) {
  const out = { title: null, num_sequences: null, total_bases: null, volumes: [] };
  const lines = String(text || '').split('\n');
  let inVolumes = false;
  for (const raw of lines) {
    const line = raw.replace(/\r$/, '');
    let m;
    if ((m = /^Database:\s*(.*)$/.exec(line))) out.title = m[1].trim();
    else if ((m = /^\s*([\d,]+) sequences; ([\d,]+) total bases/.exec(line))) {
      out.num_sequences = Number(m[1].replace(/,/g, ''));
      out.total_bases = Number(m[2].replace(/,/g, ''));
    } else if (/^Volumes:/.test(line)) inVolumes = true;
    else if (inVolumes && line.trim()) out.volumes.push(line.trim());
  }
  return out;
}

// Runs blastdbcmd -info through the injected spawnLines. Throws BLASTDBCMD_FAILED
// when it exits non-zero or prints no sequence count.
async function dbInfo(opts) {
  const o = opts || {};
  const spawnLines = o.spawnLines || spawnLinesLocal;
  const lines = [];
  const res = await spawnLines(o.cmd, ['-info', '-db', String(o.db).replace(/\.(nal|nin)$/, '')], {
    timeoutMs: o.timeoutMs > 0 ? o.timeoutMs : 30000,
    signal: o.signal,
    cwd: o.cwd,
    onLine(line) { if (lines.length < 10000) lines.push(line); }
  });
  const info = parseDbInfo(lines.join('\n'));
  if (res.code !== 0 || info.num_sequences === null) {
    throw codedError('BLASTDBCMD_FAILED', 'blastdbcmd -info failed' + (res.code !== 0 ? ' with code ' + res.code : ''),
      { stderrTail: res.stderrTail });
  }
  return info;
}

module.exports = {
  GENOME_OUTFMT,
  CDNA_OUTFMT,
  MIN_MAX_TARGET_SEQS,
  DEFAULT_TIMEOUT_MS,
  STDERR_TAIL_BYTES,
  KILL_GRACE_MS,
  CHILD_ENV,
  buildArgs,
  engineDescription,
  targetFromArgs,
  uniquePrimers,
  queryFasta,
  parseLine,
  spawnLinesLocal,
  runBlast,
  parseDbInfo,
  dbInfo
};
