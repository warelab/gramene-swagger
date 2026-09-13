'use strict';

// BLAST hits → primer binding sites, filters and columnar per-primer storage
// (spec §B.4).
//
//   t5 = qstart-1 ; t3 = qlen-qend
//   plus : face 'F'; p5 = sstart - t5; p3 = send + t3     (3' end points to higher coordinates)
//   minus: face 'R'; p5 = sstart + t5; p3 = send - t3
//   lbU  = mismatch + ceil(t5/2) + ceil(t3/2)             (substitution-only lower bound)
//   lbI  = mismatch + (t5>0) + (t3>0)                     (gap-aware lower bound)
//   parse filter (param-independent): keep if lbU <= 5 || (lbI <= 5 && t5+t3 <= 8)
//   job filter (before pairing):      lb = (t5+t3 <= 8) ? min(lbU, lbI) : lbU; drop if lb >= ignore_mismatches
//   near_perfect(primer) = hits with t5==0 && t3==0 && mismatch <= 1
//   site cap per primer = max(50000, ceil(60 * total_bases / 1e6)); over it → truncated + repetitive
//
// Gapped sites. An indel near the middle of a primer of 22 nt or more makes BLAST (ungapped) report the two
// sides as separate HSPs whose unaligned tails are >= 11 nt each, so both fail the parse filter. A hit that
// fails it but has lbI <= 5 and >= 11 aligned bases is held; finalize() (after BLAST) joins held pairs of one
// primer, subject and strand whose extrapolated p5 differ by <= 3 and that are colinear (5' half A before
// 3' half B in both query and subject) into one site:
//   dq = B.qstart - A.qend - 1, ds = subject bases between the halves; k = max(0, -dq, -ds) overlap columns
//   are trimmed from B; p5 = p5(A), p3 = p3(B), t5 = t5(A), t3 = t3(B),
//   mm = mm(A) + mm(B without the trimmed columns) + max(1, |dq+k - (ds+k)|)
//   kept when parseKeep(mm, t5, t3) and the footprint length is within 3 of the primer length.
// cDNA stores join qseq/sseq (primer bases between the halves against N, missing bases as '-').
// Re-alignment then scores the joined site exactly (genome) or from the joined strings (cDNA).

const PARSE_MAX_LB = 5;
const MAX_GAP_TAILS = 8;
const MIN_SITE_CAP = 50000;
const SITES_PER_MB = 60;
const MIN_HALF_ALIGNED = 11;
const MERGE_P5_TOLERANCE = 3;

function lowerBoundU(mm, t5, t3) {
  return mm + Math.ceil(t5 / 2) + Math.ceil(t3 / 2);
}

function lowerBoundI(mm, t5, t3) {
  return mm + (t5 > 0 ? 1 : 0) + (t3 > 0 ? 1 : 0);
}

function parseKeep(mm, t5, t3) {
  return lowerBoundU(mm, t5, t3) <= PARSE_MAX_LB ||
    (lowerBoundI(mm, t5, t3) <= PARSE_MAX_LB && t5 + t3 <= MAX_GAP_TAILS);
}

function jobBound(mm, t5, t3) {
  const u = lowerBoundU(mm, t5, t3);
  return t5 + t3 <= MAX_GAP_TAILS ? Math.min(u, lowerBoundI(mm, t5, t3)) : u;
}

// A hit that fails the parse filter but may be one side of a gapped site.
function holdable(mm, t5, t3, aligned) {
  return lowerBoundI(mm, t5, t3) <= PARSE_MAX_LB && aligned >= MIN_HALF_ALIGNED;
}

function siteKey(primer, subject, face, p5) {
  return primer + ':' + subject + ':' + face + ':' + p5;
}

// hit (blast.parseLine, optionally with hit.primer) → site object
//   { key, q, primer, len, subject, strand: 1|-1, face: 'F'|'R', p5, p3, t5, t3, mm, lbU, lbI, qseq?, sseq? }
function hitToSite(hit) {
  const t5 = hit.qstart - 1;
  const t3 = hit.qlen - hit.qend;
  const plus = hit.strand === 1 || hit.sstrand === 'plus';
  const mm = hit.mismatch;
  const face = plus ? 'F' : 'R';
  const p5 = plus ? hit.sstart - t5 : hit.sstart + t5;
  const site = {
    key: null,
    q: hit.q,
    primer: hit.primer != null ? hit.primer : null,
    len: hit.qlen,
    subject: hit.sseqid,
    strand: plus ? 1 : -1,
    face,
    p5,
    p3: plus ? hit.send + t3 : hit.send - t3,
    t5,
    t3,
    mm,
    lbU: lowerBoundU(mm, t5, t3),
    lbI: lowerBoundI(mm, t5, t3)
  };
  site.key = siteKey(site.primer != null ? site.primer : 'q' + site.q, site.subject, face, p5);
  if (hit.qseq != null) {
    site.qseq = hit.qseq;
    site.sseq = hit.sseq;
    site.qstart = hit.qstart;
    site.qend = hit.qend;
  }
  return site;
}

function passesParseFilter(site) {
  return parseKeep(site.mm, site.t5, site.t3);
}

function jobLowerBound(site) {
  return jobBound(site.mm, site.t5, site.t3);
}

function passesJobFilter(site, ignoreMismatches) {
  return jobLowerBound(site) < ignoreMismatches;
}

function isNearPerfect(site) {
  return site.t5 === 0 && site.t3 === 0 && site.mm <= 1;
}

function isNearPerfectHit(hit) {
  return hit.qstart === 1 && hit.qend === hit.qlen && hit.mismatch <= 1;
}

function siteCap(totalBases) {
  const tb = Number(totalBases) > 0 ? Number(totalBases) : 0;
  return Math.max(MIN_SITE_CAP, Math.ceil(SITES_PER_MB * tb / 1e6));
}

const COLUMNS = ['subj', 'p5', 'p3', 'strand', 't5', 't3', 'mm'];
const HOLD_COLUMNS = ['subj', 'qstart', 'qend', 'sstart', 'send', 'strand', 'mm'];

function newColumns(names, capacity) {
  const c = {};
  for (const name of names) c[name] = new Int32Array(capacity);
  return c;
}

function growColumns(names, cols, capacity) {
  const c = {};
  for (const name of names) {
    c[name] = new Int32Array(capacity);
    c[name].set(cols[name]);
  }
  return c;
}

// Per-task site store: one set of Int32Array columns per primer.
//   new SiteStore({ primers: string[] (index = BLAST q), totalBases, siteCap?, keepAlignments? })
// addHit(hit) counts the hit, the near-perfect counter and applies the parse
// filter and the cap (holding possible halves of gapped sites). Subjects are interned (store.subjects[i]).
// finalize() joins held halves; call it once after BLAST (view() and sites() call it when needed).
class SiteStore {
  constructor(opts) {
    const o = opts || {};
    if (!Array.isArray(o.primers) || o.primers.length === 0) throw new TypeError('primers must be a non-empty array');
    this.primers = o.primers.map((s) => String(s).toUpperCase());
    this.cap = Number.isInteger(o.siteCap) && o.siteCap > 0 ? o.siteCap : siteCap(o.totalBases);
    this.keepAlignments = !!o.keepAlignments;
    this.subjects = [];
    this._subjectIndex = new Map();
    this._primerIndex = new Map();
    this.primers.forEach((s, i) => { if (!this._primerIndex.has(s)) this._primerIndex.set(s, i); });
    const initial = Number.isInteger(o.initialCapacity) && o.initialCapacity > 0 ? o.initialCapacity : 256;
    this._cols = this.primers.map(() => newColumns(COLUMNS, initial));
    this._aln = this.primers.map(() => (this.keepAlignments ? [] : null));
    this._stats = this.primers.map((seq) => ({ seq, len: seq.length, hits: 0, sites: 0, near_perfect_sites: 0, truncated: false }));
    this._held = this.primers.map(() => null);
    this._heldSubjects = [];
    this._heldSubjectIndex = new Map();
    this._pendingHeld = false;
    this.gappedSites = 0;
  }

  indexOf(primerOrIndex) {
    if (Number.isInteger(primerOrIndex)) {
      return primerOrIndex >= 0 && primerOrIndex < this.primers.length ? primerOrIndex : -1;
    }
    const i = this._primerIndex.get(String(primerOrIndex).toUpperCase());
    return i === undefined ? -1 : i;
  }

  _q(primerOrIndex) {
    const q = this.indexOf(primerOrIndex);
    if (q < 0) throw new RangeError('unknown primer ' + primerOrIndex);
    return q;
  }

  _intern(subject) {
    let i = this._subjectIndex.get(subject);
    if (i === undefined) {
      i = this.subjects.length;
      this.subjects.push(subject);
      this._subjectIndex.set(subject, i);
    }
    return i;
  }

  // Appends one site to primer q's columns; false (and truncated) past the cap.
  _store(q, subject, p5, p3, strand, t5, t3, mm, aln) {
    const st = this._stats[q];
    if (st.sites >= this.cap) {
      st.truncated = true;
      return false;
    }
    let cols = this._cols[q];
    const n = st.sites;
    if (n === cols.p5.length) {
      cols = this._cols[q] = growColumns(COLUMNS, cols, n * 2);
    }
    cols.subj[n] = this._intern(subject);
    cols.p5[n] = p5;
    cols.p3[n] = p3;
    cols.strand[n] = strand;
    cols.t5[n] = t5;
    cols.t3[n] = t3;
    cols.mm[n] = mm;
    if (this.keepAlignments) this._aln[q].push(aln || null);
    st.sites = n + 1;
    return true;
  }

  // Returns true when the hit was stored as a site.
  addHit(hit) {
    const q = hit.q;
    const st = this._stats[q];
    if (!st) throw new RangeError('hit for unknown query index ' + q);
    st.hits++;
    const t5 = hit.qstart - 1;
    const t3 = hit.qlen - hit.qend;
    const mm = hit.mismatch;
    if (t5 === 0 && t3 === 0 && mm <= 1) st.near_perfect_sites++;
    if (!parseKeep(mm, t5, t3)) {
      if (holdable(mm, t5, t3, hit.qend - hit.qstart + 1)) this._hold(q, hit);
      return false;
    }
    const plus = hit.strand === 1 || hit.sstrand === 'plus';
    return this._store(q, hit.sseqid, plus ? hit.sstart - t5 : hit.sstart + t5, plus ? hit.send + t3 : hit.send - t3,
      plus ? 1 : -1, t5, t3, mm, hit.qseq != null ? [hit.qseq, hit.sseq] : null);
  }

  _hold(q, hit) {
    let h = this._held[q];
    if (!h) h = this._held[q] = { n: 0, cols: newColumns(HOLD_COLUMNS, 64), aln: this.keepAlignments ? [] : null };
    if (h.n >= this.cap) return;
    if (h.n === h.cols.mm.length) h.cols = growColumns(HOLD_COLUMNS, h.cols, h.n * 2);
    let si = this._heldSubjectIndex.get(hit.sseqid);
    if (si === undefined) {
      si = this._heldSubjects.length;
      this._heldSubjects.push(hit.sseqid);
      this._heldSubjectIndex.set(hit.sseqid, si);
    }
    const c = h.cols;
    const n = h.n;
    c.subj[n] = si;
    c.qstart[n] = hit.qstart;
    c.qend[n] = hit.qend;
    c.sstart[n] = hit.sstart;
    c.send[n] = hit.send;
    c.strand[n] = hit.strand === 1 || hit.sstrand === 'plus' ? 1 : -1;
    c.mm[n] = hit.mismatch;
    if (h.aln) h.aln.push(hit.qseq != null ? [hit.qseq, hit.sseq] : null);
    h.n = n + 1;
    this._pendingHeld = true;
  }

  // Joins held halves into gapped sites (see the header). Returns the number of sites added.
  finalize() {
    let added = 0;
    for (let q = 0; q < this._held.length; q++) {
      const h = this._held[q];
      this._held[q] = null;
      if (h && h.n >= 2) added += this._joinHeld(q, h);
    }
    this._heldSubjects = [];
    this._heldSubjectIndex = new Map();
    this._pendingHeld = false;
    this.gappedSites += added;
    return added;
  }

  _joinHeld(q, h) {
    const c = h.cols;
    const count = h.n;
    const p5 = new Int32Array(count);
    for (let i = 0; i < count; i++) {
      const t5 = c.qstart[i] - 1;
      p5[i] = c.strand[i] === 1 ? c.sstart[i] - t5 : c.sstart[i] + t5;
    }
    const order = new Array(count);
    for (let i = 0; i < count; i++) order[i] = i;
    order.sort((a, b) => c.subj[a] - c.subj[b] || c.strand[a] - c.strand[b] || p5[a] - p5[b] || c.qstart[a] - c.qstart[b]);
    const used = new Uint8Array(count);
    let added = 0;
    for (let x = 0; x < count; x++) {
      const i = order[x];
      if (used[i]) continue;
      for (let y = x + 1; y < count; y++) {
        const j = order[y];
        if (c.subj[j] !== c.subj[i] || c.strand[j] !== c.strand[i] || p5[j] - p5[i] > MERGE_P5_TOLERANCE) break;
        if (used[j]) continue;
        const a = c.qstart[i] <= c.qstart[j] ? i : j;
        const joined = this._join(q, h, a, a === i ? j : i);
        if (!joined) continue;
        used[i] = 1;
        used[j] = 1;
        if (this._store(q, this._heldSubjects[c.subj[i]], joined.p5, joined.p3, c.strand[i], joined.t5, joined.t3, joined.mm, joined.aln)) added++;
        break;
      }
    }
    return added;
  }

  // Held hit a (5' half) + b (3' half) → { p5, p3, t5, t3, mm, aln } or null.
  _join(q, h, a, b) {
    const c = h.cols;
    const st = this._stats[q];
    const len = st.len;
    const plus = c.strand[a] === 1;
    if (!(c.qstart[a] < c.qstart[b] && c.qend[a] < c.qend[b])) return null;
    if (plus ? !(c.sstart[a] < c.sstart[b]) : !(c.sstart[a] > c.sstart[b])) return null;
    // cDNA sites are re-aligned from qseq/sseq, so both halves need them
    if (h.aln && !(h.aln[a] && h.aln[b])) return null;
    const dq = c.qstart[b] - c.qend[a] - 1;
    const ds = plus ? c.sstart[b] - c.send[a] - 1 : c.send[a] - c.sstart[b] - 1;
    const k = Math.max(0, -dq, -ds);
    if (k >= c.qend[b] - c.qstart[b] + 1) return null;
    const dq2 = dq + k;
    const ds2 = ds + k;
    const t5 = c.qstart[a] - 1;
    const t3 = len - c.qend[b];
    let mmB = Math.max(0, c.mm[b] - k);
    let aln = null;
    if (h.aln && h.aln[a] && h.aln[b]) {
      const A = h.aln[a];
      const bq = h.aln[b][0].slice(k);
      const bs = h.aln[b][1].slice(k);
      mmB = 0;
      for (let i = 0; i < bq.length; i++) if (bq[i].toUpperCase() !== bs[i].toUpperCase()) mmB++;
      const midQ = st.seq.slice(c.qend[a], c.qend[a] + dq2);
      const common = Math.min(dq2, ds2);
      const qMid = dq2 > ds2 ? midQ : midQ + '-'.repeat(ds2 - dq2);
      const sMid = 'N'.repeat(common) + (dq2 > ds2 ? '-'.repeat(dq2 - ds2) : 'N'.repeat(ds2 - dq2));
      aln = [A[0] + qMid + bq, A[1] + sMid + bs];
    }
    const mm = c.mm[a] + mmB + Math.max(1, Math.abs(dq2 - ds2));
    if (!parseKeep(mm, t5, t3)) return null;
    const p5 = plus ? c.sstart[a] - t5 : c.sstart[a] + t5;
    const p3 = plus ? c.send[b] + t3 : c.send[b] - t3;
    if (Math.abs(Math.abs(p3 - p5) + 1 - len) > MERGE_P5_TOLERANCE) return null;
    return { p5, p3, t5, t3, mm, aln };
  }

  size(primerOrIndex) {
    return this._stats[this._q(primerOrIndex)].sites;
  }

  // { seq, len, hits, sites, near_perfect_sites, truncated }
  stats(primerOrIndex) {
    return Object.assign({}, this._stats[this._q(primerOrIndex)]);
  }

  // results.primers[seq] core fields (§B.12): { len, near_perfect_sites, repetitive, truncated }
  primerSummary(primerOrIndex, params) {
    const st = this._stats[this._q(primerOrIndex)];
    const threshold = params && params.repeat_site_threshold != null ? params.repeat_site_threshold : 5;
    return {
      len: st.len,
      near_perfect_sites: st.near_perfect_sites,
      repetitive: st.truncated || st.near_perfect_sites > threshold,
      truncated: st.truncated
    };
  }

  // Materialize site i of a primer as a hitToSite-shaped object.
  site(primerOrIndex, i) {
    const q = this._q(primerOrIndex);
    const st = this._stats[q];
    if (!(i >= 0 && i < st.sites)) throw new RangeError('site index out of range');
    const c = this._cols[q];
    const t5 = c.t5[i];
    const t3 = c.t3[i];
    const mm = c.mm[i];
    const plus = c.strand[i] === 1;
    const s = {
      key: null,
      q,
      primer: st.seq,
      len: st.len,
      subject: this.subjects[c.subj[i]],
      strand: plus ? 1 : -1,
      face: plus ? 'F' : 'R',
      p5: c.p5[i],
      p3: c.p3[i],
      t5,
      t3,
      mm,
      lbU: lowerBoundU(mm, t5, t3),
      lbI: lowerBoundI(mm, t5, t3)
    };
    s.key = siteKey(s.primer, s.subject, s.face, s.p5);
    if (this.keepAlignments && this._aln[q][i]) {
      s.qseq = this._aln[q][i][0];
      s.sseq = this._aln[q][i][1];
      s.qstart = t5 + 1;
      s.qend = st.len - t3;
    }
    return s;
  }

  // Columnar view for amplicons.candidates: indices of sites passing the job
  // filter (all sites when ignoreMismatches is not given). Joins pending held halves first.
  //   { q, primer, len, n, idx: Int32Array, subj, p5, strand, lb: Int32Array, subjects }
  view(primerOrIndex, opts) {
    if (this._pendingHeld) this.finalize();
    const q = this._q(primerOrIndex);
    const st = this._stats[q];
    const c = this._cols[q];
    const limit = opts && Number.isInteger(opts.ignoreMismatches) ? opts.ignoreMismatches : Infinity;
    const idx = new Int32Array(st.sites);
    const lb = new Int32Array(st.sites);
    let n = 0;
    for (let i = 0; i < st.sites; i++) {
      const b = jobBound(c.mm[i], c.t5[i], c.t3[i]);
      if (b < limit) {
        idx[n] = i;
        lb[n] = b;
        n++;
      }
    }
    return {
      q,
      primer: st.seq,
      len: st.len,
      n,
      idx: idx.subarray(0, n),
      lb: lb.subarray(0, n),
      subj: c.subj,
      p5: c.p5,
      strand: c.strand,
      subjects: this.subjects
    };
  }

  // Array of materialized sites, optionally job-filtered.
  sites(primerOrIndex, opts) {
    const v = this.view(primerOrIndex, opts);
    const out = new Array(v.n);
    for (let k = 0; k < v.n; k++) out[k] = this.site(v.q, v.idx[k]);
    return out;
  }

  totalSites() {
    let n = 0;
    for (const st of this._stats) n += st.sites;
    return n;
  }
}

module.exports = {
  PARSE_MAX_LB,
  MAX_GAP_TAILS,
  MIN_SITE_CAP,
  SITES_PER_MB,
  MIN_HALF_ALIGNED,
  MERGE_P5_TOLERANCE,
  lowerBoundU,
  lowerBoundI,
  siteKey,
  hitToSite,
  passesParseFilter,
  jobLowerBound,
  passesJobFilter,
  isNearPerfect,
  isNearPerfectHit,
  siteCap,
  SiteStore
};
