'use strict';

// Check orchestration (spec §B.1, §B.2, §B.12, §B.15).
//
//   tasks = [reference genome (ws5, -num_threads ctx.procs)]
//         + (mode transcript ? [reference cDNA (ws5)] : [])
//         + (checks has pangenome ? genomes sorted by (total_bases, system_name), ≤ ctx.procs single-thread
//            blastn at ws6, target cdna in transcript mode : [])
//   per task: blast → sites (+ joined gapped sites) → amplicons.candidates → realign plan → realign
//             (FASTA DP | cDNA qseq/sseq) → finalize → annotate (overlaps | transcript→gene groups)
//             → verdicts / pan-genome statuses
//
// ctx (A.8.5): { jobId, siteKey, signal, config (frozen primers cfg), resolved (job.resolved), procs,
//   progress({done,total,stage,running}), partial(results), tmpdir(), log {info,warn,error},
//   spawnLines(cmd, args, {stdin, timeoutMs, onLine, signal, cwd}) → {code, signal, stderrTail}, mongo,
//   fatalOnMongoUnavailable? }
// The results object handed to ctx.partial is the live document; serialize it synchronously.
//
// Working directory: every blastn / blastdbcmd child runs with cwd = a private directory: ctx.tmpdir() when it
// is a directory owned by this user and not group/world-writable, else a mkdtemp (0700) directory created
// here (under the directory ctx.tmpdir() returned, cfg.tmp_dir or os.tmpdir()) and removed at the end.
//
// Re-alignment (genome targets): a site whose gap-aware lower bound (amplicons.realignBound = BLAST mismatches
// + 1 per unaligned tail) exceeds max_amplifying_mismatches cannot be part of an amplifying product, so it is
// not re-aligned against the FASTA; it gets an approximate alignment (realign.realignFromHit with mm = that
// bound, approx: true) and its products are classed unlikely. The max_realign_sites_per_genome cap applies to
// the re-aligned sites only.
//
// Failures: a BLAST (or blastdbcmd) failure is retried once after 2 s; a BLAST_TIMEOUT is not retried. A
// failed reference throws err.code 'REFERENCE_BLAST_FAILED'; a failed pan-genome genome is marked 'error'.
// Abort throws signal.reason when it is an Error with a string code, else err.code 'ABORTED'.
// Mongo: when annotation or transcript→gene mapping is unavailable because a supplied mongo handle failed
// and ctx.fatalOnMongoUnavailable is true (the worker), run() throws err.code 'MONGO_UNAVAILABLE' with
// err.fatal = true; otherwise it adds warning ANNOTATION_UNAVAILABLE and completes.
//
// Genotyping (genotyping spec §5.6-§5.9), only for a request with a genotyping block in gene or region mode:
// before the reference stage, check/genotype.js validateSets re-derives the prepared variant and sets over
// this job's FASTA access and results.genotyping is added after pangenome; the reference stage calls the
// allele caller for the reference (the control) before its flush; every pan-genome genome records an entry
// (unavailable when its DB is missing or its BLAST failed; otherwise the caller, after markOrthologs); each
// flush rebuilds results.genotyping. The megablast fallback runs through ctx.spawnLines in the private cwd:
// genotype_megablast_timeout_ms, no retry, at most genotype_max_megablast per job. Requests without
// genotyping produce exactly the results they always did.

const fs = require('fs');
const os = require('os');
const path = require('path');

const blast = require('./blast');
const sites = require('./sites');
const amplicons = require('./amplicons');
const realign = require('./realign');
const classify = require('./classify');
const sensitivity = require('./sensitivity');
const specificity = require('./specificity');
const pangenome = require('./pangenome');
const annotate = require('./annotate');
const { redactPaths } = require('../errors');

const DEFAULT_ALGORITHM_VERSION = '2';
const BLAST_RETRY_DELAY_MS = 2000;
const REALIGN_MAX_FETCH = 2000000;
const PARAM_KEYS = Object.freeze(Object.keys(classify.DEFAULT_PARAMS));
const NOOP_LOG = Object.freeze({ info() {}, warn() {}, error() {} });
const MAX_SUBJECTS_IN_WARNING = 10;
// Genotyping fallbacks (primers.check.genotype_max_megablast / genotype_megablast_timeout_ms) and the megablast rows kept
// (-max_target_seqs 50 x -max_hsps 10 is at most 500).
const GENOTYPE_MAX_MEGABLAST = 30;
const GENOTYPE_MEGABLAST_TIMEOUT_MS = 20000;
const GENOTYPE_MAX_HSP_ROWS = 1000;

function codedError(code, message, extra) {
  const e = new Error(message);
  e.code = code;
  if (extra) Object.assign(e, extra);
  return e;
}

function posInt(v, dflt) {
  return Number.isInteger(v) && v > 0 ? v : dflt;
}

function pickParams(params) {
  const p = classify.withDefaults(params);
  const out = {};
  for (const k of PARAM_KEYS) out[k] = p[k];
  return out;
}

function normalizePairs(pairs) {
  if (!Array.isArray(pairs) || pairs.length === 0) throw new TypeError('request.pairs must be a non-empty array');
  return pairs.map((p, i) => {
    if (!p || typeof p.left !== 'string' || typeof p.right !== 'string') throw new TypeError('pair ' + i + ' needs left and right');
    return { id: p.id != null ? String(p.id) : 'pair' + (i + 1), left: p.left.toUpperCase(), right: p.right.toUpperCase(), expected: p.expected || null };
  });
}

function algorithmVersion(ctx) {
  if (ctx && ctx.algorithmVersion != null) return String(ctx.algorithmVersion);
  try {
    const index = require('./index');
    if (index && index.ALGORITHM_VERSION != null) return String(index.ALGORITHM_VERSION);
  } catch (e) { /* check/index.js not present (unit tests) */ }
  return DEFAULT_ALGORITHM_VERSION;
}

// A directory BLAST+ may run in: absolute, owned by this user, not group/world-writable.
function isPrivateDir(dir) {
  if (typeof dir !== 'string' || !path.isAbsolute(dir)) return false;
  try {
    const st = fs.statSync(dir);
    if (!st.isDirectory() || (st.mode & 0o022) !== 0) return false;
    return typeof process.getuid !== 'function' || st.uid === process.getuid();
  } catch (e) {
    return false;
  }
}

function isDirectory(dir) {
  try {
    return typeof dir === 'string' && path.isAbsolute(dir) && fs.statSync(dir).isDirectory();
  } catch (e) {
    return false;
  }
}

// Client-safe {code, message} for a failed BLAST/blastdbcmd run (stderr tail, paths → basenames).
function publicError(e) {
  const code = e && typeof e.code === 'string' && /^[A-Z][A-Z0-9_]*$/.test(e.code) ? e.code : 'BLAST_FAILED';
  let message = e && e.message ? String(e.message) : 'BLAST failed';
  if (e && e.stderrTail) {
    const tail = String(e.stderrTail).trim().split('\n').filter(Boolean).slice(-3).join(' | ');
    if (tail) message += ': ' + tail;
  }
  return { code, message: redactPaths(message).slice(0, 500) };
}

const versionCache = new Map();

// 'blastn: 2.13.0+' → '2.13.0' (cached per binary; null when it cannot be determined).
async function blastVersion(cmd, spawnLines, signal, cwd) {
  if (versionCache.has(cmd)) return versionCache.get(cmd);
  let text = '';
  try {
    const res = await spawnLines(cmd, ['-version'], {
      timeoutMs: 15000,
      signal,
      cwd,
      onLine(line) { if (text.length < 2000) text += line + '\n'; }
    });
    const m = /blastn:\s*(\d+\.\d+\.\d+)/.exec(text) || /blast\s+(\d+\.\d+\.\d+)/.exec(text);
    if (res && res.code === 0 && m) {
      versionCache.set(cmd, m[1]);
      return m[1];
    }
  } catch (e) { /* unknown version */ }
  return null;
}

function clearVersionCache() {
  versionCache.clear();
}

// Warnings deduplicated by (code, message); subjects (genome names) are appended to the message. details (the first
// add's) is kept only when given, so {code, message} warnings stay as they were.
class WarningSet {
  constructor() {
    this.map = new Map();
  }

  add(code, message, subject, details) {
    const key = code + ' ' + message;
    let w = this.map.get(key);
    if (!w) {
      w = { code, message, subjects: [] };
      if (details !== undefined) w.details = details;
      this.map.set(key, w);
    }
    if (subject != null && w.subjects.indexOf(subject) < 0) w.subjects.push(subject);
  }

  has(code) {
    for (const w of this.map.values()) if (w.code === code) return true;
    return false;
  }

  toArray() {
    return Array.from(this.map.values()).map((w) => {
      let out;
      if (!w.subjects.length) {
        out = { code: w.code, message: w.message };
      } else {
        const shown = w.subjects.slice(0, MAX_SUBJECTS_IN_WARNING).join(', ');
        const more = w.subjects.length > MAX_SUBJECTS_IN_WARNING ? ' and ' + (w.subjects.length - MAX_SUBJECTS_IN_WARNING) + ' more' : '';
        out = { code: w.code, message: w.message + ' (' + shown + more + ')' };
      }
      if (w.details !== undefined) out.details = w.details;
      return out;
    });
  }
}

// Annotation priority: every selection's head first, then round-robin over the lists.
function priorityList(sets) {
  const out = [];
  for (const s of sets) for (const a of s.head || []) out.push(a);
  const queues = [];
  for (const s of sets) for (const l of s.lists || []) if (l && l.length) queues.push({ l, i: 0 });
  let more = true;
  while (more) {
    more = false;
    for (const q of queues) {
      if (q.i < q.l.length) {
        out.push(q.l[q.i++]);
        more = true;
      }
    }
  }
  return out;
}

function defaultFetchWindow() {
  const sequence = require('../sequence');
  return (fastaPath, region, start1, end1, strand) => sequence.fetch(fastaPath, region, start1, end1, strand, { maxLength: REALIGN_MAX_FETCH });
}

function defaultRegionLength() {
  const sequence = require('../sequence');
  return (fastaPath, region) => Promise.resolve().then(() => sequence.regionLength(fastaPath, region)).catch(() => undefined);
}

class CheckRun {
  constructor(request, ctx, opts) {
    const o = opts || {};
    this.request = request;
    this.ctx = ctx;
    this.cfg = ctx.config || {};
    this.ccfg = this.cfg.check || {};
    this.log = ctx.log && typeof ctx.log.warn === 'function' ? ctx.log : NOOP_LOG;
    this.external = ctx.signal || null;
    this.internal = new AbortController();
    this.signal = this.external ? AbortSignal.any([this.external, this.internal.signal]) : this.internal.signal;
    this.spawnLines = typeof ctx.spawnLines === 'function' ? ctx.spawnLines : blast.spawnLinesLocal;
    this.procs = posInt(ctx.procs, 1);
    this.blastn = this.cfg.blastn || 'blastn';
    this.blastdbcmd = this.cfg.blastdbcmd || 'blastdbcmd';
    this.mode = request.mode || 'region';
    this.transcript = this.mode === 'transcript';
    this.checks = Array.isArray(request.checks) && request.checks.length ? request.checks : ['specificity'];
    this.doPan = this.checks.indexOf('pangenome') >= 0;
    this.params = pickParams(request.params);
    this.resolved = ctx.resolved || {};
    this.assemblies = this.resolved.assemblies || {};
    this.gene = annotate.geneContext(this.resolved.gene);
    this.pairs = normalizePairs(request.pairs);
    this.primers = blast.uniquePrimers([].concat(...this.pairs.map((p) => [p.left, p.right])));
    this.genomes = this.doPan && Array.isArray(request.genomes) ? Array.from(new Set(request.genomes.map(String))) : [];
    this.refWs = posInt(this.ccfg.reference_word_size, 5);
    this.panWs = posInt(this.ccfg.pangenome_word_size, 6);
    this.maxListed = posInt(this.ccfg.max_offtargets_listed, specificity.DEFAULT_MAX_LISTED);
    this.maxCandidates = posInt(this.ccfg.max_candidates_per_pair, amplicons.DEFAULT_MAX_CANDIDATES);
    this.maxRealignSites = posInt(this.ccfg.max_realign_sites_per_genome, amplicons.DEFAULT_MAX_REALIGN_SITES);
    this.realignConcurrency = posInt(this.ccfg.realign_concurrency, 32);
    this.blastTimeoutMs = posInt(this.ccfg.blast_timeout_ms, blast.DEFAULT_TIMEOUT_MS);
    this.retryDelayMs = Number.isInteger(o.retryDelayMs) && o.retryDelayMs >= 0 ? o.retryDelayMs : BLAST_RETRY_DELAY_MS;
    this.fetchWindow = o.fetchWindow || defaultFetchWindow();
    this.regionLength = o.regionLength || defaultRegionLength();
    this.algorithmVersion = algorithmVersion(ctx);
    this.warnings = new WarningSet();
    this.annotator = new annotate.Annotator({
      mongo: ctx.mongo,
      log: this.log,
      maxAnnotated: posInt(this.ccfg.max_annotated_amplicons, annotate.DEFAULT_MAX_ANNOTATED),
      timeoutMs: o.mongoTimeoutMs
    });
    this.dbInfoCache = new Map();
    // Genotyping spec §5.6-§5.9: allele calls for a request with a genotyping block (gene and region modes only).
    this.genotyping = request.genotyping && (this.mode === 'gene' || this.mode === 'region') ? request.genotyping : null;
    this.geno = null;
    this.genoReference = null;
    this.genoSpecificity = null;
    this.genoEntries = new Map();
    this.genoMegablast = 0;
    this.genoQuery = null;
    this.genoWeak = null;
    this.done = 0;
    this.total = 1 + (this.transcript ? 1 : 0) + this.genomes.length;
    this.results = null;
    this.cwd = null;
    this.ownTmp = null;
    this.tag = ctx.jobId ? String(ctx.jobId).slice(0, 12) : 'job';
  }

  // ---- abort / retry -------------------------------------------------------------------------

  abortError() {
    const reason = this.external && this.external.aborted ? this.external.reason : null;
    if (reason instanceof Error && typeof reason.code === 'string') return reason;
    return codedError('ABORTED', 'check aborted', { aborted: true });
  }

  isAbort(e) {
    return (this.external && this.external.aborted) || this.internal.signal.aborted || (e && e.code === 'ABORTED');
  }

  checkAborted() {
    if ((this.external && this.external.aborted) || this.internal.signal.aborted) throw this.abortError();
  }

  sleep(ms) {
    return new Promise((resolve, reject) => {
      if (this.signal.aborted) {
        reject(this.abortError());
        return;
      }
      const onAbort = () => {
        clearTimeout(timer);
        reject(this.abortError());
      };
      const timer = setTimeout(() => {
        this.signal.removeEventListener('abort', onAbort);
        resolve();
      }, ms);
      this.signal.addEventListener('abort', onAbort, { once: true });
    });
  }

  // One retry after retryDelayMs, except for BLAST_TIMEOUT: a stalled run is not repeated (spec §B.10
  // 'error': BLAST failed twice or timed out), so it cannot hold a slot or the job for a second timeout.
  async withRetry(label, fn) {
    let last = null;
    for (let attempt = 1; attempt <= 2; attempt++) {
      this.checkAborted();
      try {
        return await fn(attempt);
      } catch (e) {
        if (this.isAbort(e)) throw this.abortError();
        last = e;
        this.log.warn('primers check ' + this.tag + ': ' + label + ' attempt ' + attempt + ' failed: ' + publicError(e).message);
        if (e && e.code === 'BLAST_TIMEOUT') break;
        if (attempt === 1) await this.sleep(this.retryDelayMs);
      }
    }
    if (last && typeof last === 'object') last.blastStage = true;
    throw last;
  }

  // ---- ctx side channels ----------------------------------------------------------------------

  safeCall(fn, arg) {
    if (typeof fn !== 'function') return;
    try {
      const p = fn(arg);
      if (p && typeof p.catch === 'function') p.catch((e) => this.log.warn('primers check ' + this.tag + ': ctx callback failed: ' + (e && e.message)));
    } catch (e) {
      this.log.warn('primers check ' + this.tag + ': ctx callback failed: ' + (e && e.message));
    }
  }

  progress(stage, running) {
    this.safeCall(this.ctx.progress, { done: this.done, total: this.total, stage, running: running || [] });
  }

  flush(stage, running) {
    this.results.warnings = this.warnings.toArray();
    this.progress(stage, running);
    this.safeCall(this.ctx.partial, this.results);
  }

  addAssemblyWarnings(asm) {
    if (!asm || !Array.isArray(asm.warnings)) return;
    for (const w of asm.warnings) {
      if (w && w.code) this.warnings.add(String(w.code), String(w.message || w.code), asm.system_name);
    }
  }

  // Private cwd for BLAST+ children (see the header). Idempotent.
  async privateCwd() {
    if (this.cwd) return this.cwd;
    let dir = null;
    if (typeof this.ctx.tmpdir === 'function') {
      try {
        dir = await this.ctx.tmpdir();
      } catch (e) {
        this.log.warn('primers check ' + this.tag + ': ctx.tmpdir failed: ' + (e && e.message));
        dir = null;
      }
    }
    if (!isPrivateDir(dir)) {
      const cfgTmp = typeof this.cfg.tmp_dir === 'string' && isDirectory(this.cfg.tmp_dir) ? this.cfg.tmp_dir : os.tmpdir();
      const base = isDirectory(dir) ? dir : cfgTmp;
      dir = fs.mkdtempSync(path.join(base, 'primers-check-'));
      this.ownTmp = dir;
    }
    this.cwd = dir;
    return dir;
  }

  removeOwnTmp() {
    if (!this.ownTmp) return;
    try {
      fs.rmSync(this.ownTmp, { recursive: true, force: true });
    } catch (e) {
      this.log.warn('primers check ' + this.tag + ': could not remove ' + path.basename(this.ownTmp) + ': ' + (e && e.message));
    }
    this.ownTmp = null;
  }

  // ---- one BLAST target ---------------------------------------------------------------------

  async cdnaInfo(asm) {
    const db = asm.blastdb.cdna;
    if (this.dbInfoCache.has(db)) return this.dbInfoCache.get(db);
    const info = await blast.dbInfo({ cmd: this.blastdbcmd, db, spawnLines: this.spawnLines, signal: this.signal, cwd: await this.privateCwd() });
    this.dbInfoCache.set(db, info);
    return info;
  }

  // Alignments for the unique candidate sites of one task → { byKey, realigned, bounded, truncated, fetchErrors }.
  async alignSites(asm, target, lists) {
    const W = this.params.three_prime_window;
    if (target === 'cdna') {
      const uniq = amplicons.uniqueSites(lists, { maxSites: this.maxRealignSites });
      const alns = uniq.sites.map((s) => realign.realignCdnaHit(s, { threePrimeWindow: W }));
      return { byKey: realign.indexAlignments(uniq.sites, alns), realigned: uniq.sites.length, bounded: 0, truncated: uniq.truncated, fetchErrors: 0 };
    }
    const plan = amplicons.realignPlan(lists, { maxSites: this.maxRealignSites, maxAmplifyingMismatches: this.params.max_amplifying_mismatches });
    const fastaPath = asm.fasta && asm.fasta.dna ? asm.fasta.dna : null;
    const alns = await realign.realignGenomeSites(plan.realign, {
      fetchWindow: this.fetchWindow,
      regionLength: this.regionLength,
      fastaPath,
      threePrimeWindow: W,
      concurrency: this.realignConcurrency,
      signal: this.signal
    });
    let fetchErrors = 0;
    for (const a of alns) if (a && a.fetch_error) fetchErrors++;
    const byKey = realign.indexAlignments(plan.realign, alns);
    // Bounded sites keep their extrapolated ends, clamped to the sequence like re-aligned ones.
    const lengths = new Map();
    for (const s of plan.bounded) {
      let regionLength;
      if (fastaPath && typeof this.regionLength === 'function') {
        if (!lengths.has(s.subject)) lengths.set(s.subject, await this.regionLength(fastaPath, s.subject));
        regionLength = lengths.get(s.subject);
      }
      byKey.set(s.key, realign.realignFromHit(s, W, { mm: amplicons.realignBound(s), regionLength }));
    }
    return { byKey, realigned: plan.realign.length, bounded: plan.bounded.length, truncated: plan.truncated, fetchErrors };
  }

  // → { perPair: [{pair, maxSize, amplicons, unlikely, counts, truncated, candidates}], stats }
  async analyze(asm, target, opts) {
    const db = target === 'cdna' ? asm.blastdb.cdna : asm.blastdb.dna;
    const cwd = await this.privateCwd();
    const blasted = await this.withRetry(opts.label, async () => {
      let maxTargetSeqs;
      let totalBases;
      if (target === 'cdna') {
        const info = await this.cdnaInfo(asm);
        maxTargetSeqs = info.num_sequences;
        totalBases = info.total_bases;
      } else {
        maxTargetSeqs = asm.num_sequences;
        totalBases = asm.total_bases;
      }
      const cap = Math.max(blast.MIN_MAX_TARGET_SEQS, Number.isInteger(maxTargetSeqs) ? maxTargetSeqs : 0);
      const store = new sites.SiteStore({ primers: this.primers, totalBases, keepAlignments: target === 'cdna' });
      // Distinct subjects per primer: BLAST prints all HSPs of one subject together, so subject
      // changes in the stream count subjects without keeping a set.
      const lastSubject = new Array(this.primers.length).fill(null);
      const subjects = new Array(this.primers.length).fill(0);
      const res = await blast.runBlast({
        cmd: this.blastn,
        args: blast.buildArgs({ target, wordSize: opts.wordSize, threads: opts.threads, db, maxTargetSeqs }),
        primers: this.primers,
        spawnLines: this.spawnLines,
        timeoutMs: this.blastTimeoutMs,
        signal: this.signal,
        cwd,
        target,
        onHit: (h) => {
          if (h.sseqid !== lastSubject[h.q]) {
            lastSubject[h.q] = h.sseqid;
            subjects[h.q]++;
          }
          store.addHit(h);
        }
      });
      const gapped = store.finalize();
      return { store, cap, subjects, hits: res.hits, gapped, blastMs: res.elapsedMs };
    });
    this.checkAborted();

    const store = blasted.store;
    const t0 = Date.now();
    const perPair = this.pairs.map((pair, i) => ({
      pair,
      maxSize: opts.maxSizes[i],
      cand: amplicons.candidates(pair, store, { maxSize: opts.maxSizes[i], maxCandidates: this.maxCandidates, ignoreMismatches: this.params.ignore_mismatches })
    }));
    const candidatesMs = Date.now() - t0;
    const t1 = Date.now();
    const aligned = await this.alignSites(asm, target, perPair.map((x) => x.cand.candidates));
    const realignMs = Date.now() - t1;
    this.checkAborted();
    if (aligned.fetchErrors > 0) {
      this.warnings.add('NO_FASTA_FOR_REALIGN', 'genome sequence was not available for re-alignment of some sites; their mismatch counts are lower bounds (approx)', asm.system_name);
    }
    for (const x of perPair) {
      const fin = amplicons.finalize(x.cand.candidates, aligned.byKey, this.params);
      const qs = [store.indexOf(x.pair.left), store.indexOf(x.pair.right)];
      const primerTruncated = qs.some((q) => q >= 0 && (store.stats(q).truncated || blasted.subjects[q] >= blasted.cap));
      x.amplicons = fin.amplicons;
      x.unlikely = fin.unlikely;
      x.counts = fin.counts;
      x.candidates = x.cand.candidates.length;
      x.truncated = !!(x.cand.truncated || fin.counts.unaligned > 0 || primerTruncated);
      delete x.cand;
    }
    return {
      store,
      perPair,
      stats: {
        hits: blasted.hits,
        sites: store.totalSites(),
        gapped_sites: blasted.gapped,
        realigned: aligned.realigned,
        bounded: aligned.bounded,
        realign_truncated: aligned.truncated,
        fetch_errors: aligned.fetchErrors,
        blast_ms: blasted.blastMs,
        candidates_ms: candidatesMs,
        realign_ms: realignMs
      }
    };
  }

  async referenceAnalyze(asm, target, label) {
    try {
      return await this.analyze(asm, target, { wordSize: this.refWs, threads: this.procs, maxSizes: this.pairs.map(() => this.params.max_product_size), label });
    } catch (e) {
      if (this.isAbort(e)) throw this.abortError();
      if (!e || !e.blastStage) throw e;
      const pe = publicError(e);
      throw codedError('REFERENCE_BLAST_FAILED', label + ' BLAST failed: ' + pe.message, { details: { target, cause: pe.code } });
    }
  }

  // ---- annotation -----------------------------------------------------------------------------

  annotationUnavailable() {
    if (this.ctx.fatalOnMongoUnavailable === true && this.annotator.mongoFailed) {
      throw codedError('MONGO_UNAVAILABLE', 'gene annotation (mongo) is unavailable', { fatal: true });
    }
    this.warnings.add('ANNOTATION_UNAVAILABLE', 'gene annotation (mongo) is unavailable; genes and ortholog flags are null');
  }

  async annotateGenome(asm, sets) {
    const list = priorityList(sets);
    if (!list.length) return;
    const res = await this.annotator.overlaps(asm, list);
    if (!res.available) this.annotationUnavailable();
  }

  async groupCdna(asm, perPair) {
    const subjects = new Set();
    for (const x of perPair) {
      for (const a of x.amplicons) subjects.add(a.region);
      for (const a of x.unlikely) subjects.add(a.region);
    }
    let seed = null;
    if (this.gene && this.gene.id && asm.system_name === this.request.system_name) {
      seed = new Map(this.gene.transcripts.map((tx) => [tx, this.gene.id]));
    }
    const tg = await this.annotator.transcriptGenes(asm, Array.from(subjects), seed);
    if (!tg.available) this.annotationUnavailable();
    if (tg.unmapped.length) {
      this.warnings.add('TRANSCRIPT_GENE_UNMAPPED', 'some transcripts have no gene in mongo and are grouped by their id without the isoform suffix', asm.system_name);
    }
    return perPair.map((x) => ({ groups: annotate.groupByGene(x.amplicons, tg.map), unlikely: annotate.groupByGene(x.unlikely, tg.map) }));
  }

  // ---- stages ---------------------------------------------------------------------------------

  async referenceStage(ref) {
    const r = this.results;
    this.progress('reference', [ref.system_name]);
    const t0 = Date.now();
    const res = await this.referenceAnalyze(ref, 'genome', 'reference genome');
    for (const seq of this.primers) {
      const summary = res.store.primerSummary(seq, this.params);
      summary.sensitivity = sensitivity.primerSensitivity(seq.length, { referenceWordSize: this.refWs, pangenomeWordSize: this.doPan ? this.panWs : undefined });
      r.primers[seq] = summary;
    }
    const includeUnlikely = !!this.params.include_unlikely;
    const sels = res.perPair.map((x) => specificity.selectGenome({
      amplicons: x.amplicons,
      unlikely: x.unlikely,
      truncated: x.truncated,
      mode: this.mode,
      expected: specificity.honoredExpected(this.mode, x.pair),
      geneLocation: this.gene ? this.gene.location : null
    }));
    this.checkAborted();
    await this.annotateGenome(ref, sels.map((s) => specificity.annotationCandidates(s, { maxListed: this.maxListed, includeUnlikely })));
    r.specificity = {
      target: 'genome',
      pairs: sels.map((s, i) => specificity.genomePairBlock(s, { id: this.pairs[i].id, maxListed: this.maxListed, includeUnlikely }))
    };
    r.timings_ms.reference = Date.now() - t0;
    this.stats.reference = res.stats;
    if (this.geno) await this.genotypeReference(ref, res.perPair, sels);
    this.done++;
    this.flush('reference', []);
  }

  async transcriptomeStage(ref) {
    const r = this.results;
    this.progress('transcriptome', [ref.system_name]);
    const t0 = Date.now();
    const res = await this.referenceAnalyze(ref, 'cdna', 'reference cDNA');
    this.checkAborted();
    const grouped = await this.groupCdna(ref, res.perPair);
    const geneId = this.request.gene_id != null ? String(this.request.gene_id) : this.gene ? this.gene.id : null;
    const includeUnlikely = !!this.params.include_unlikely;
    const sels = grouped.map((g, i) => specificity.selectCdna({ groups: g.groups, unlikelyGroups: g.unlikely, truncated: res.perPair[i].truncated, geneId }));
    r.transcriptome = {
      target: 'cdna',
      gene_id: geneId,
      transcript_id: this.request.transcript_id != null ? String(this.request.transcript_id) : null,
      pairs: sels.map((s, i) => specificity.cdnaPairBlock(s, { id: this.pairs[i].id, maxListed: this.maxListed, includeUnlikely }))
    };
    r.timings_ms.transcriptome = Date.now() - t0;
    this.stats.transcriptome = res.stats;
    this.done++;
    this.flush('transcriptome', []);
  }

  resolveFailure(sys, asm) {
    if (asm && asm.error && asm.error.code) return String(asm.error.code);
    const errs = this.resolved.errors;
    if (errs && errs[sys]) return String(errs[sys].code || errs[sys]);
    return asm ? 'NO_BLASTDB' : 'UNRESOLVED';
  }

  async pangenomeGenome(sys) {
    const asm = this.assemblies[sys] || null;
    const target = this.transcript ? 'cdna' : 'genome';
    const orth = annotate.orthologInfo(this.gene, sys);
    const info = { system_name: sys, display_name: asm && asm.display_name ? asm.display_name : sys, orthologAnnotated: orth.annotated };
    const db = asm && asm.blastdb ? (target === 'cdna' ? asm.blastdb.cdna : asm.blastdb.dna) : null;
    if (!db || (asm && asm.error)) {
      const reason = this.resolveFailure(sys, asm);
      if (this.geno) this.genoEntries.set(sys, require('./genotype').unavailableEntry(this.geno, info, 'db_unavailable'));
      return this.pairs.map(() => pangenome.unavailableEntry(Object.assign({ reason }, info)));
    }
    this.addAssemblyWarnings(asm);
    let res;
    try {
      res = await this.analyze(asm, target, { wordSize: this.panWs, threads: 1, maxSizes: this.panMax, label: 'pan-genome ' + sys + ' ' + target });
    } catch (e) {
      if (this.isAbort(e)) throw this.abortError();
      const err = e && e.blastStage ? publicError(e) : { code: 'CHECK_FAILED', message: 'analysis failed' };
      if (!(e && e.blastStage)) this.log.error('primers check ' + this.tag + ': pan-genome ' + sys + ' failed: ' + (e && e.stack ? e.stack : e));
      if (this.geno) this.genoEntries.set(sys, require('./genotype').unavailableEntry(this.geno, info, e && e.blastStage ? 'blast_error' : 'call_failed'));
      return this.pairs.map(() => pangenome.errorEntry(Object.assign({ error: err }, info)));
    }
    this.stats.pangenome = this.stats.pangenome || {};
    this.stats.pangenome[sys] = res.stats;
    this.checkAborted();
    if (target === 'genome') {
      await this.annotateGenome(asm, res.perPair.map((x) => ({ head: [], lists: [x.amplicons.length ? x.amplicons : x.unlikely.slice(0, 1)] })));
      // Ortholog flags first: the allele caller tests copies for annotated orthologs (genotyping spec §5.9).
      for (const x of res.perPair) {
        pangenome.markOrthologs(x.amplicons, orth.ids);
        pangenome.markOrthologs(x.unlikely, orth.ids);
      }
      if (this.geno) this.genoEntries.set(sys, await this.genotypeGenome(asm, res.perPair, info));
      return res.perPair.map((x, pi) => pangenome.genomeEntry(Object.assign({ products: x.amplicons, unlikely: x.unlikely, referenceSize: this.refSizes[pi], truncated: x.truncated }, info)));
    }
    const grouped = await this.groupCdna(asm, res.perPair);
    return res.perPair.map((x, pi) => {
      pangenome.markOrthologs(grouped[pi].groups, orth.ids);
      pangenome.markOrthologs(grouped[pi].unlikely, orth.ids);
      return pangenome.genomeEntry(Object.assign({ products: grouped[pi].groups, unlikely: grouped[pi].unlikely, referenceSize: this.refSizes[pi], truncated: x.truncated }, info));
    });
  }

  async pangenomeStage() {
    const r = this.results;
    const target = this.transcript ? 'cdna' : 'genome';
    if (this.transcript) {
      this.warnings.add('PANGENOME_TRANSCRIPT_MODELS_ONLY', 'transcript-mode pan-genome checks search annotated transcript models (cDNA BLAST databases) only; copies outside annotated transcripts are not seen');
    }
    const block = target === 'cdna' ? r.transcriptome : r.specificity;
    const transcriptId = this.request.transcript_id != null ? String(this.request.transcript_id) : null;
    this.refSizes = this.pairs.map((p, i) => specificity.referenceSize(block && block.pairs[i], { transcriptId }));
    this.panMax = this.refSizes.map((s) => pangenome.pangenomeMaxSize(this.params.max_product_size, s));
    r.pangenome = {
      target,
      pairs: this.pairs.map((p, i) => ({ id: p.id, reference_size: this.refSizes[i], max_size: this.panMax[i], summary: pangenome.summarize([]), genomes: [] }))
    };
    const entries = this.pairs.map(() => new Array(this.genomes.length).fill(null));
    const tb = (sys) => {
      const a = this.assemblies[sys];
      return a && Number.isFinite(a.total_bases) ? a.total_bases : Infinity;
    };
    const order = this.genomes.map((g, i) => i).sort((a, b) => {
      const da = tb(this.genomes[a]);
      const db = tb(this.genomes[b]);
      if (da !== db) return da < db ? -1 : 1;
      return this.genomes[a] < this.genomes[b] ? -1 : this.genomes[a] > this.genomes[b] ? 1 : 0;
    });
    const running = new Set();
    let next = 0;
    let failure = null;
    let failureIsAbort = false;
    const worker = async () => {
      while (next < order.length && !failure) {
        this.checkAborted();
        const idx = order[next++];
        const sys = this.genomes[idx];
        running.add(sys);
        this.progress('pangenome', Array.from(running));
        const t0 = Date.now();
        let out;
        try {
          out = await this.pangenomeGenome(sys);
        } finally {
          running.delete(sys);
        }
        r.timings_ms[sys] = Date.now() - t0;
        out.forEach((e, pi) => { entries[pi][idx] = e; });
        r.pangenome.pairs.forEach((pb, pi) => {
          pb.genomes = entries[pi].filter(Boolean);
          pb.summary = pangenome.summarize(pb.genomes);
        });
        this.done++;
        if (this.geno) this.writeGenotyping();
        this.flush('pangenome', Array.from(running));
      }
    };
    const n = Math.max(1, Math.min(this.procs, order.length));
    const settled = await Promise.allSettled(Array.from({ length: n }, () => worker().catch((e) => {
      if (!failure) {
        failure = e;
        // Decide before aborting the other workers: the internal abort would make every error look like one.
        failureIsAbort = !!this.isAbort(e);
        if (!failureIsAbort) this.internal.abort();
      }
      throw e;
    })));
    if (this.external && this.external.aborted) throw this.abortError();
    if (failure) throw failureIsAbort ? this.abortError() : failure;
    for (const s of settled) if (s.status === 'rejected') throw s.reason;
  }

  // ---- genotyping (genotyping spec §5.6-§5.9) -------------------------------------------------------

  // Before the reference stage: prepare() re-derived from the stored request over this job's FASTA access, the sequences
  // past the submit-time window filled in, and results.genotyping added after pangenome.
  async genotypeStart(ref) {
    const genotype = require('./genotype');
    const fasta = ref.fasta && ref.fasta.dna ? ref.fasta.dna : null;
    const prepared = await genotype.validateSets(this.genotyping, this.pairs, ref, {
      cfg: this.cfg,
      sequence: {
        regionLength: (fastaPath, region) => this.regionLength(fastaPath, region),
        fetch: (fastaPath, region, start, end, strand) => this.fetchWindow(fastaPath, region, start, end, strand)
      }
    });
    await genotype.fillSequences(prepared, (region, start, end) => this.fetchWindow(fasta, region, start, end, 1));
    this.checkAborted();
    this.geno = prepared;
    const results = {};
    for (const k of Object.keys(this.results)) {
      results[k] = this.results[k];
      if (k === 'pangenome') results.genotyping = genotype.emptyResults(prepared);
    }
    this.results = results;
  }

  // Per set, the REF and ALT pairs' products of one task.
  genotypeProducts(perPair) {
    const byId = new Map();
    this.pairs.forEach((p, i) => byId.set(p.id, perPair[i]));
    const side = (x) => (x ? { amplicons: x.amplicons, unlikely: x.unlikely } : null);
    return this.geno.sets.map((s) => ({ ref: side(byId.get(s.ref_pair)), alt: side(byId.get(s.alt_pair)) }));
  }

  // The allele call and set predictions of one genome with a genome target. A caller exception makes the genome unavailable
  // (call_failed, warning GENOTYPE_FAILED); an abort propagates. stats.genotype[system_name] records the stage's cost.
  async genotypeGenome(asm, perPair, info) {
    const genotype = require('./genotype');
    const t0 = Date.now();
    const cpu0 = process.cpuUsage();
    const fasta = asm && asm.fasta && asm.fasta.dna ? asm.fasta.dna : null;
    let out = null;
    let entry;
    try {
      out = await genotype.callGenome({
        prepared: this.geno,
        system_name: info.system_name,
        display_name: info.display_name,
        is_reference: info.is_reference === true,
        products: this.genotypeProducts(perPair),
        cfg: this.ccfg,
        params: this.params
      }, {
        fetch: (region, start, end) => (fasta ? this.fetchWindow(fasta, region, start, end, 1) : Promise.reject(codedError('NO_SEQUENCE', 'no genome FASTA'))),
        regionLength: (region) => (fasta ? this.regionLength(fasta, region) : Promise.resolve(undefined)),
        megablast: () => this.genotypeMegablast(asm, info.system_name)
      });
      entry = { genome: out.genome, sets: out.sets };
      if (out.weak_off_targets.length) this.genotypeWeakOffTargets(info, out.weak_off_targets);
    } catch (e) {
      if (this.isAbort(e)) throw this.abortError();
      this.log.error('primers check ' + this.tag + ': allele caller failed for ' + info.system_name + ': ' + (e && e.stack ? e.stack : e));
      this.warnings.add('GENOTYPE_FAILED', 'the allele caller failed; these genomes are unavailable', info.system_name);
      entry = genotype.unavailableEntry(this.geno, info, 'call_failed');
    }
    this.checkAborted();
    const cpu = process.cpuUsage(cpu0);
    this.stats.genotype = this.stats.genotype || {};
    this.stats.genotype[info.system_name] = {
      ms: Date.now() - t0,
      cpu_ms: Math.round((cpu.user + cpu.system) / 1000),
      anchors: out ? out.anchors : null,
      failed_reads: out ? out.failed_reads : null,
      megablast: out ? out.megablast : null
    };
    return entry;
  }

  // WEAK_OFF_TARGETS (M8b): one warning for the job naming the genomes with off-locus products left out of the predictions by
  // genotype_offlocus_max_mismatches; its details object is shared by every add and grows with each genome.
  genotypeWeakOffTargets(info, items) {
    const genotype = require('./genotype');
    if (!this.genoWeak) this.genoWeak = genotype.weakOffTargetDetails(this.ccfg);
    const order = { genomes: [this.request.system_name].concat(this.genomes), sets: this.geno.sets.map((s) => s.id) };
    genotype.addWeakOffTargets(this.genoWeak, info.system_name, items, order);
    this.warnings.add('WEAK_OFF_TARGETS', 'off-target products with more than ' + this.genoWeak.max_mismatches + ' mismatches in a primer do not change ' +
      'the allele predictions; the specificity and pan-genome results still list them', info.system_name, this.genoWeak);
  }

  // §5.8 reference control: the reference's own call, its per-set verdicts and REFERENCE_CONTROL_FAILED.
  async genotypeReference(ref, perPair, sels) {
    const index = new Map(this.pairs.map((p, i) => [p.id, i]));
    this.genoSpecificity = this.geno.sets.map((s) => {
      const a = sels[index.get(s.ref_pair)];
      const b = sels[index.get(s.alt_pair)];
      const off = new Set(a.off.concat(b.off).map((x) => x.region + ':' + x.start + ':' + x.end + ':' + x.orientation));
      return { ref_verdict: a.verdict, alt_verdict: b.verdict, off_target_count: off.size };
    });
    const sys = ref.system_name || this.request.system_name;
    this.genoReference = await this.genotypeGenome(ref, perPair, { system_name: sys, display_name: ref.display_name || sys, is_reference: true });
    this.writeGenotyping();
    const failed = this.results.genotyping.sets.filter((s) => s.control && s.control.status === 'fail').map((s) => s.id);
    if (failed.length) {
      const allele = this.genoReference.genome.allele;
      this.warnings.add('REFERENCE_CONTROL_FAILED', 'the reference control failed for ' + (failed.length > 1 ? 'sets ' : 'set ') + failed.join(', ') +
        ': the reference ' + sys + ' must be called ref (it is ' + allele + ') and every set must predict ref on it', null, { allele, sets: failed });
    }
  }

  // §5.6 step 8: one megablast of the reference segment around the variant against asm's genome DB, within the per-job
  // budget, without retry. → {status: 'ok', rows, query} | {status: 'budget'} | {status: 'failed'}
  async genotypeMegablast(asm, sys) {
    const max = Number.isInteger(this.ccfg.genotype_max_megablast) && this.ccfg.genotype_max_megablast >= 0 ? this.ccfg.genotype_max_megablast : GENOTYPE_MAX_MEGABLAST;
    if (this.genoMegablast >= max) {
      this.warnings.add('GENOTYPE_FALLBACK_BUDGET', 'the per-job megablast fallback budget (' + max + ') was reached; these genomes are missing', sys);
      return { status: 'budget' };
    }
    this.genoMegablast++;
    const failed = (message) => {
      this.log.warn('primers check ' + this.tag + ': megablast fallback for ' + sys + ' failed: ' + message);
      this.warnings.add('GENOTYPE_FALLBACK_FAILED', 'the megablast fallback of the allele caller failed; these genomes are missing', sys);
      return { status: 'failed' };
    };
    try {
      const query = await this.genotypeQuery();
      const args = blast.buildMegablastArgs({ db: asm && asm.blastdb ? asm.blastdb.dna : null, threads: 1 });
      const rows = [];
      let bad = 0;
      const res = await this.spawnLines(this.blastn, args, {
        stdin: blast.megablastQueryFasta(query.sequence),
        timeoutMs: posInt(this.ccfg.genotype_megablast_timeout_ms, GENOTYPE_MEGABLAST_TIMEOUT_MS),
        signal: this.signal,
        cwd: await this.privateCwd(),
        onLine: (line) => {
          let row = null;
          try {
            row = blast.parseMegablastLine(line);
          } catch (e) {
            bad++;
          }
          if (row && rows.length < GENOTYPE_MAX_HSP_ROWS) rows.push(row);
        }
      });
      if (res.aborted || this.signal.aborted) throw this.abortError();
      if (res.timedOut) return failed('megablast timed out');
      if (res.code !== 0) return failed(publicError(Object.assign(new Error('blastn exited with code ' + res.code), { stderrTail: res.stderrTail })).message);
      if (bad > 0) return failed(bad + ' unparseable megablast line(s)');
      return { status: 'ok', rows, query };
    } catch (e) {
      if (this.isAbort(e)) throw this.abortError();
      return failed(publicError(e).message);
    }
  }

  // The megablast query {region, start, end, sequence}, read once per job from the reference FASTA.
  genotypeQuery() {
    if (!this.genoQuery) {
      const q = require('./genotype').megablastQuery(this.geno);
      const ref = this.assemblies[this.request.system_name];
      this.genoQuery = Promise.resolve()
        .then(() => this.fetchWindow(ref.fasta.dna, q.region, q.start, q.end, 1))
        .then((seq) => Object.assign({ sequence: String(seq).toUpperCase() }, q));
      this.genoQuery.catch(() => { this.genoQuery = null; });
    }
    return this.genoQuery;
  }

  // results.genotyping for the next flush: the reference, then the finished pan-genome genomes in request order.
  writeGenotyping() {
    require('./genotype').writeResults(this.results.genotyping, {
      reference: this.genoReference,
      specificity: this.genoSpecificity,
      genomes: this.genomes.map((sys) => this.genoEntries.get(sys)).filter(Boolean)
    });
  }

  async execute() {
    const started = Date.now();
    this.checkAborted();
    const ref = this.assemblies[this.request.system_name];
    if (!ref || !ref.blastdb || !ref.blastdb.dna) throw codedError('NO_BLASTDB', 'the reference genome has no dna BLAST database');
    if (this.transcript && !ref.blastdb.cdna) throw codedError('NO_BLASTDB', 'the reference genome has no cDNA BLAST database');
    try {
      const cwd = await this.privateCwd();
      const version = await blastVersion(this.blastn, this.spawnLines, this.signal, cwd);
      this.checkAborted();
      this.stats = {};
      this.results = {
        engine: {
          algorithm_version: this.algorithmVersion,
          blast: version,
          reference: blast.engineDescription(this.refWs),
          pangenome: this.doPan ? blast.engineDescription(this.panWs) : null
        },
        params: Object.assign({}, this.params),
        reference: {
          system_name: ref.system_name || this.request.system_name,
          map_id: ref.map_id != null ? ref.map_id : null,
          total_bases: ref.total_bases != null ? ref.total_bases : null
        },
        sensitivity_note: sensitivity.sensitivityNote({
          maxAmplifyingMismatches: this.params.max_amplifying_mismatches,
          referenceWordSize: this.refWs,
          pangenomeWordSize: this.doPan ? this.panWs : undefined,
          primerLengths: this.primers.map((seq) => seq.length)
        }),
        primers: {},
        specificity: null,
        transcriptome: null,
        pangenome: null,
        warnings: [],
        timings_ms: {}
      };
      this.addAssemblyWarnings(ref);
      if (this.genotyping) await this.genotypeStart(ref);
      try {
        await this.referenceStage(ref);
        if (this.transcript) await this.transcriptomeStage(ref);
        if (this.doPan) await this.pangenomeStage();
      } finally {
        if (!this.internal.signal.aborted) this.internal.abort();
      }
    } finally {
      this.removeOwnTmp();
    }
    this.results.warnings = this.warnings.toArray();
    this.results.timings_ms.total = Date.now() - started;
    return this.results;
  }
}

// run(request, ctx, opts?) → §B.12 results.
//   request: normalized PrimerCheckRequest {system_name, mode, gene_id?, transcript_id?, checks, genomes,
//            params (8 keys), pairs [{id, left, right, expected?}]}
//   opts (tests): { retryDelayMs = 2000, fetchWindow, regionLength, mongoTimeoutMs }
async function run(request, ctx, opts) {
  if (!request || typeof request !== 'object') throw new TypeError('request is required');
  if (!ctx || typeof ctx !== 'object') throw new TypeError('ctx is required');
  const job = new CheckRun(request, ctx, opts);
  return job.execute();
}

// Diagnostics for tests and logs: runs like run() and also returns per-stage statistics.
async function runWithStats(request, ctx, opts) {
  if (!request || typeof request !== 'object') throw new TypeError('request is required');
  if (!ctx || typeof ctx !== 'object') throw new TypeError('ctx is required');
  const job = new CheckRun(request, ctx, opts);
  const results = await job.execute();
  return { results, stats: job.stats };
}

module.exports = {
  DEFAULT_ALGORITHM_VERSION,
  BLAST_RETRY_DELAY_MS,
  run,
  runWithStats,
  _internal: { CheckRun, WarningSet, publicError, priorityList, blastVersion, clearVersionCache, pickParams, normalizePairs, isPrivateDir }
};
