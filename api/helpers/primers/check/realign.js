'use strict';

// Re-alignment of primer sites (spec §B.6).
//
// Genome target: fetch FASTA[subject][min(p5,p3)-3 .. max(p5,p3)+3] (plus strand;
// face R is reverse-complemented here, IUPAC-aware) and align the whole primer
// against it:
//   match 0, mismatch 1, gap 1; leading and trailing window bases are free;
//   the primer is fully consumed; its last base is aligned to a window base (no
//   terminal gap at the 3' end); among equal-cost ends the one nearest the
//   expected 3' end (window index m-3 when unclipped) wins, then the lower index.
//   Traceback prefers diagonal, then primer-base gap, then window-base gap, which
//   places gaps toward the 5' end.
// Edits are reported as distances from the primer 3' end (1 = terminal base):
//   substitution or inserted primer base i → n-i+1; a window base skipped after
//   primer base i → n-i+1. mm_pos is ordered 5' → 3' (descending distance).
//
// cDNA target (the cdna .fa.gz is not indexed): exact positions over the aligned
// part from qseq/sseq, plus lbU-style tail edits (ceil(t/2) per tail, placed at
// the unaligned base next to the HSP and every second base outward);
// approx = (t5+t3 > 0).

const PAD = 3;
const DEFAULT_THREE_PRIME_WINDOW = 5;
const DEFAULT_CONCURRENCY = 32;
const MERGE_GAP = 4096;
const MAX_MERGED_SPAN = 65536;

const DIR_DIAG = 1;
const DIR_UP = 2;   // primer base against a gap (insertion in the primer)
const DIR_LEFT = 3; // window base against a gap (deletion in the primer)

const COMPLEMENT = (() => {
  const t = {};
  const pairs = ['AT', 'CG', 'RY', 'KM', 'SS', 'WW', 'BV', 'DH', 'NN'];
  for (const [a, b] of pairs) {
    t[a] = b; t[b] = a;
    t[a.toLowerCase()] = b.toLowerCase(); t[b.toLowerCase()] = a.toLowerCase();
  }
  t.U = 'A'; t.u = 'a';
  return t;
})();

// IUPAC-aware, case-preserving reverse complement. Unknown characters map to N/n.
function revcomp(seq) {
  let out = '';
  for (let i = seq.length - 1; i >= 0; i--) {
    const c = seq[i];
    const r = COMPLEMENT[c];
    out += r !== undefined ? r : (c === c.toLowerCase() && c !== c.toUpperCase() ? 'n' : (c === '-' || c === '.' ? c : 'N'));
  }
  return out;
}

function countWithin(dists, w) {
  let n = 0;
  for (const d of dists) if (d <= w) n++;
  return n;
}

function codedError(code, message) {
  const e = new Error(message);
  e.code = code;
  return e;
}

// Semi-global DP of `primer` (n) against `window` (m).
// expectedEnd: 1-based window index of the expected 3' end (tie-break target).
// Returns { mm, mm_3p, terminal_mm, gaps, subs, mm_pos, startIdx, endIdx, ops }
// where startIdx/endIdx are 1-based window indices of the primer 5' and 3' ends
// (startIdx is extrapolated over leading primer-base gaps) and ops is the
// alignment as a string of M (match) X (mismatch) I (primer base vs gap)
// D (window base vs gap), 5' → 3'.
function alignPrimer(primer, window, expectedEnd, threePrimeWindow) {
  const p = String(primer).toUpperCase();
  const w = String(window).toUpperCase();
  const n = p.length;
  const m = w.length;
  const W = Number.isInteger(threePrimeWindow) ? threePrimeWindow : DEFAULT_THREE_PRIME_WINDOW;
  if (n === 0) throw new TypeError('empty primer');
  if (m === 0) return null;
  const cols = m + 1;
  const D = new Int32Array((n + 1) * cols);
  const dir = new Uint8Array((n + 1) * cols);
  for (let i = 1; i <= n; i++) {
    D[i * cols] = i;
    dir[i * cols] = DIR_UP;
  }
  // Rows 1..n-1 in full; the last row only needs diagonal ends.
  for (let i = 1; i < n; i++) {
    const pc = p.charCodeAt(i - 1);
    const row = i * cols;
    const prev = (i - 1) * cols;
    for (let j = 1; j <= m; j++) {
      const diag = D[prev + j - 1] + (pc === w.charCodeAt(j - 1) ? 0 : 1);
      const up = D[prev + j] + 1;
      const left = D[row + j - 1] + 1;
      let best = diag;
      let d = DIR_DIAG;
      if (up < best) { best = up; d = DIR_UP; }
      if (left < best) { best = left; d = DIR_LEFT; }
      D[row + j] = best;
      dir[row + j] = d;
    }
  }
  const lastChar = p.charCodeAt(n - 1);
  const prevRow = (n - 1) * cols;
  const target = Number.isFinite(expectedEnd) ? expectedEnd : m - PAD;
  let endIdx = -1;
  let endCost = Infinity;
  for (let j = 1; j <= m; j++) {
    const cost = D[prevRow + j - 1] + (lastChar === w.charCodeAt(j - 1) ? 0 : 1);
    if (cost < endCost || (cost === endCost && Math.abs(j - target) < Math.abs(endIdx - target))) {
      endCost = cost;
      endIdx = j;
    }
  }

  // Traceback from the forced diagonal at (n, endIdx).
  const ops = [];
  const dists = [];
  let subs = 0;
  let gaps = 0;
  if (lastChar === w.charCodeAt(endIdx - 1)) ops.push('M');
  else {
    ops.push('X');
    dists.push(1);
    subs++;
  }
  let i = n - 1;
  let j = endIdx - 1;
  while (i > 0) {
    const d = j === 0 ? DIR_UP : dir[i * cols + j];
    if (d === DIR_DIAG) {
      if (p.charCodeAt(i - 1) === w.charCodeAt(j - 1)) ops.push('M');
      else {
        ops.push('X');
        dists.push(n - i + 1);
        subs++;
      }
      i--;
      j--;
    } else if (d === DIR_UP) {
      ops.push('I');
      dists.push(n - i + 1);
      gaps++;
      i--;
    } else {
      ops.push('D');
      dists.push(n - i + 1);
      gaps++;
      j--;
    }
  }
  ops.reverse();
  let leadingIns = 0;
  while (leadingIns < ops.length && ops[leadingIns] === 'I') leadingIns++;
  const startIdx = j + 1 - leadingIns;
  dists.sort((a, b) => b - a);
  const mm = dists.length;
  if (mm !== endCost) throw new Error('realign traceback cost mismatch');
  return {
    mm,
    mm_3p: countWithin(dists, W),
    terminal_mm: dists.length > 0 && dists[dists.length - 1] === 1,
    gaps,
    subs,
    mm_pos: dists,
    startIdx,
    endIdx,
    ops: ops.join('')
  };
}

// Tail edits implied by lbU: ceil(t/2) per unaligned tail, placed at the base
// adjacent to the HSP and every second base outward. Returns distances.
function tailDistances(t5, t3, n) {
  const d = [];
  for (let k = t3; k >= 1; k -= 2) d.push(k);
  for (let i = t5; i >= 1; i -= 2) d.push(n - i + 1);
  return d;
}

// Published coordinates stay inside the sequence: [1, regionLength] (upper bound only when known).
function clampCoord(v, regionLength) {
  let x = v < 1 ? 1 : v;
  if (Number.isFinite(regionLength) && regionLength >= 1 && x > regionLength) x = regionLength;
  return x;
}

// Alignment of one genome site against its window sequence (plus strand as
// fetched). winStart is the genomic coordinate of seq[0]. Pure; used by
// realignGenomeSites and directly by tests.
// regionLength (optional): p5/p3 are clamped to [1, regionLength]. A primer whose 5' bases overhang the
// sequence start (face F) or end (face R) keeps those overhang edits in mm/mm_pos.
function realignSiteOnWindow(site, seq, winStart, threePrimeWindow, regionLength) {
  if (typeof site.primer !== 'string' || site.primer === '') {
    throw new TypeError('site.primer is required for re-alignment');
  }
  const len = seq.length;
  if (len === 0) return null;
  const winEnd = winStart + len - 1;
  const R = site.face === 'R' || site.strand === -1;
  const w = R ? revcomp(seq) : seq;
  const expectedEnd = R ? winEnd - site.p3 + 1 : site.p3 - winStart + 1;
  const a = alignPrimer(site.primer, w, expectedEnd, threePrimeWindow);
  if (!a) return null;
  const g = (idx) => (R ? winEnd - idx + 1 : winStart + idx - 1);
  return {
    key: site.key,
    mm: a.mm,
    mm_3p: a.mm_3p,
    terminal_mm: a.terminal_mm,
    gaps: a.gaps,
    mm_pos: a.mm_pos,
    p5: clampCoord(g(a.startIdx), regionLength),
    p3: clampCoord(g(a.endIdx), regionLength),
    approx: false,
    ops: a.ops
  };
}

// Fallback when no FASTA is available for a site (warning NO_FASTA_FOR_REALIGN):
// BLAST mismatch count plus lbU tail edits; internal positions are unknown, so
// mm_pos is null and mm_3p/terminal_mm only reflect the tails.
// opts: { mm: a known lower bound on the site's edits, reported instead of the lbU tail bound (run.js passes
//         amplicons.realignBound for sites that are not re-aligned because that bound already exceeds
//         max_amplifying_mismatches); mm_3p is capped at it,
//         regionLength: p5/p3 are clamped to [1, regionLength] }
function realignFromHit(site, threePrimeWindow, opts) {
  const o = opts || {};
  const W = Number.isInteger(threePrimeWindow) ? threePrimeWindow : DEFAULT_THREE_PRIME_WINDOW;
  const n = site.len || (site.primer ? site.primer.length : 0);
  const tails = tailDistances(site.t5, site.t3, n);
  const out = {
    key: site.key,
    mm: site.mm + tails.length,
    mm_3p: countWithin(tails, W),
    terminal_mm: tails.indexOf(1) >= 0,
    gaps: 0,
    mm_pos: null,
    p5: clampCoord(site.p5, o.regionLength),
    p3: clampCoord(site.p3, o.regionLength),
    approx: true
  };
  if (Number.isInteger(o.mm) && o.mm >= 0) {
    out.mm = o.mm;
    out.mm_3p = Math.min(out.mm_3p, o.mm);
    if (o.mm === 0) out.terminal_mm = false;
  }
  return out;
}

// cDNA hit/site with qseq/sseq (query orientation, as BLAST prints them).
function realignCdnaHit(hit, opts) {
  const W = opts && Number.isInteger(opts.threePrimeWindow) ? opts.threePrimeWindow : DEFAULT_THREE_PRIME_WINDOW;
  const n = hit.qlen != null ? hit.qlen : hit.len;
  const qstart = hit.qstart != null ? hit.qstart : hit.t5 + 1;
  const qend = hit.qend != null ? hit.qend : n - hit.t3;
  const t5 = qstart - 1;
  const t3 = n - qend;
  if (typeof hit.qseq !== 'string' || typeof hit.sseq !== 'string' || hit.qseq.length !== hit.sseq.length) {
    throw new TypeError('cDNA hit needs qseq and sseq of equal length');
  }
  const q = hit.qseq.toUpperCase();
  const s = hit.sseq.toUpperCase();
  const dists = tailDistances(t5, t3, n);
  let gaps = 0;
  let qi = qstart;
  for (let k = 0; k < q.length; k++) {
    const qc = q[k];
    const sc = s[k];
    if (qc === '-') {
      dists.push(n - (qi - 1) + 1);
      gaps++;
    } else if (sc === '-') {
      dists.push(n - qi + 1);
      gaps++;
      qi++;
    } else {
      if (qc !== sc) dists.push(n - qi + 1);
      qi++;
    }
  }
  dists.sort((a, b) => b - a);
  let p5 = hit.p5;
  let p3 = hit.p3;
  if (p5 === undefined && hit.sstart !== undefined) {
    const plus = hit.strand === 1 || hit.sstrand === 'plus';
    p5 = plus ? hit.sstart - t5 : hit.sstart + t5;
    p3 = plus ? hit.send + t3 : hit.send - t3;
  }
  return {
    key: hit.key,
    mm: dists.length,
    mm_3p: countWithin(dists, W),
    terminal_mm: dists.length > 0 && dists[dists.length - 1] === 1,
    gaps,
    mm_pos: dists,
    p5,
    p3,
    approx: t5 + t3 > 0
  };
}

function siteWindow(site, pad) {
  const lo = Math.min(site.p5, site.p3) - pad;
  const hi = Math.max(site.p5, site.p3) + pad;
  // Highest subject coordinate BLAST actually aligned (always inside the region).
  const alignedHi = site.face === 'R' || site.strand === -1 ? site.p5 - site.t5 : site.p3 - site.t3;
  return { lo: Math.max(1, lo), hi, alignedHi: Math.max(1, alignedHi) };
}

async function runPool(items, concurrency, fn, signal) {
  let next = 0;
  const workers = [];
  const n = Math.max(1, Math.min(concurrency, items.length));
  for (let k = 0; k < n; k++) {
    workers.push((async () => {
      while (next < items.length) {
        if (signal && signal.aborted) throw codedError('ABORTED', 'realign aborted');
        const item = items[next++];
        await fn(item);
      }
    })());
  }
  await Promise.all(workers);
}

// Re-align genome sites against FASTA windows.
//   sites: hitToSite/SiteStore.site objects (need key, primer, subject, face|strand, p5, p3, t5, t3, mm, len)
//   opts: { fetchWindow(fastaPath, region, start1, end1, strand) → Promise<string>,
//           fastaPath, threePrimeWindow = 5, concurrency = 32, signal?,
//           regionLength?(fastaPath, region) → Promise<number|undefined>, pad = 3 }
// Nearby windows on one subject are fetched together (gap ≤ 4096, span ≤ 65536).
// Returns an array parallel to `sites`:
//   { key, mm, mm_3p, terminal_mm, gaps, mm_pos[], p5, p3, approx: false, ops }
// or, when the window cannot be read, realignFromHit(site) with approx: true and
// fetch_error: <code>. Without fetchWindow/fastaPath every result is the fallback.
async function realignGenomeSites(sites, opts) {
  const o = opts || {};
  const W = Number.isInteger(o.threePrimeWindow) ? o.threePrimeWindow : DEFAULT_THREE_PRIME_WINDOW;
  const pad = Number.isInteger(o.pad) ? o.pad : PAD;
  const results = new Array(sites.length);
  if (!o.fetchWindow || !o.fastaPath) {
    for (let k = 0; k < sites.length; k++) {
      results[k] = Object.assign(realignFromHit(sites[k], W), { fetch_error: 'NO_FASTA_FOR_REALIGN' });
    }
    return results;
  }

  const jobs = sites.map((s, k) => Object.assign({ k, site: s, subject: s.subject }, siteWindow(s, pad)));

  if (o.regionLength) {
    const lengths = new Map();
    for (const job of jobs) {
      if (!lengths.has(job.subject)) lengths.set(job.subject, await o.regionLength(o.fastaPath, job.subject));
      const len = lengths.get(job.subject);
      if (Number.isFinite(len)) {
        job.hi = Math.min(job.hi, len);
        job.regionLength = len;
      }
    }
  }

  jobs.sort((a, b) => (a.subject < b.subject ? -1 : a.subject > b.subject ? 1 : a.lo - b.lo));
  const chunks = [];
  for (const job of jobs) {
    const last = chunks[chunks.length - 1];
    if (last && last.subject === job.subject && job.lo - last.hi <= MERGE_GAP &&
        Math.max(last.hi, job.hi) - last.lo + 1 <= MAX_MERGED_SPAN) {
      last.hi = Math.max(last.hi, job.hi);
      last.jobs.push(job);
    } else {
      chunks.push({ subject: job.subject, lo: job.lo, hi: job.hi, jobs: [job] });
    }
  }

  async function fetchSeq(subject, lo, hi) {
    try {
      const seq = await o.fetchWindow(o.fastaPath, subject, lo, hi, 1);
      return typeof seq === 'string' ? { seq } : { error: 'NO_SEQUENCE' };
    } catch (e) {
      return { error: (e && e.code) || 'FETCH_FAILED' };
    }
  }

  function finish(job, seq, seqStart, error) {
    let r = null;
    if (seq) {
      const from = job.lo - seqStart;
      const slice = seq.slice(Math.max(0, from), job.hi - seqStart + 1);
      if (from >= 0 && slice.length > 0) r = realignSiteOnWindow(job.site, slice, job.lo, W, job.regionLength);
    }
    results[job.k] = r || Object.assign(realignFromHit(job.site, W, { regionLength: job.regionLength }), { fetch_error: error || 'NO_SEQUENCE' });
  }

  async function singleJob(job) {
    let got = await fetchSeq(job.subject, job.lo, job.hi);
    if (got.error || got.seq.length === 0) {
      // The padded/extrapolated end may run past the region end; retry up to the
      // last aligned base, which BLAST proved exists.
      const hi2 = Math.max(job.lo, Math.min(job.hi, job.alignedHi));
      if (hi2 < job.hi) {
        const retry = await fetchSeq(job.subject, job.lo, hi2);
        if (!retry.error && retry.seq.length) {
          job.hi = hi2;
          got = retry;
        }
      }
    }
    finish(job, got.seq, job.lo, got.error);
  }

  await runPool(chunks, o.concurrency > 0 ? o.concurrency : DEFAULT_CONCURRENCY, async (chunk) => {
    const got = chunk.jobs.length > 1 ? await fetchSeq(chunk.subject, chunk.lo, chunk.hi) : { error: 'SINGLE' };
    if (got.error) {
      for (const job of chunk.jobs) await singleJob(job);
      return;
    }
    const seqEnd = chunk.lo + got.seq.length - 1;
    for (const job of chunk.jobs) {
      if (job.hi <= seqEnd) finish(job, got.seq, chunk.lo);
      else await singleJob(job);
    }
  }, o.signal);
  if (o.signal && o.signal.aborted) throw codedError('ABORTED', 'realign aborted');
  return results;
}

// Map site.key → alignment for amplicons.finalize.
function indexAlignments(sites, alignments) {
  const m = new Map();
  for (let k = 0; k < sites.length; k++) {
    if (alignments[k]) m.set(sites[k].key, alignments[k]);
  }
  return m;
}

module.exports = {
  PAD,
  DEFAULT_THREE_PRIME_WINDOW,
  revcomp,
  alignPrimer,
  tailDistances,
  realignSiteOnWindow,
  realignFromHit,
  realignCdnaHit,
  realignGenomeSites,
  indexAlignments
};
