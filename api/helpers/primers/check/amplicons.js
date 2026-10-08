'use strict';

// Amplicon calling (spec §B.5): pair forward-facing (F) and reverse-facing (R)
// sites of a primer pair on the same subject, then finalize with re-aligned
// coordinates and Primer-BLAST classes.
//
//   role(seq) = seq == pair.left ? 'L' : 'R'          (identical L/R → 'L')
//   F[], R[] sorted by (subject, p5); for f in F, r in R from
//   lowerBound(r.p5 >= f.p5 + 2 × min primer length - 1) while r.p5 - f.p5 + 1 <= maxSize,
//   keeping size >= len(f) + len(r)                   (the two primer footprints must not overlap)
//   → { subject, fwd: f, rev: r, orientation: role(f)+role(r) }   LR | RL | LL | RR
//   an identical-primer pair (left == right) labels every F×R product 'LR', so it can be the on-target
//   more than maxCandidates (5000) → truncated
// finalize drops products whose re-aligned footprints overlap (counted as invalid).

const { siteKey, jobLowerBound, passesJobFilter, lowerBoundI } = require('./sites');
const classify = require('./classify');

const DEFAULT_MAX_CANDIDATES = 5000;
const DEFAULT_MAX_REALIGN_SITES = 20000;

// Pan-genome maxSize rule from §B.5.
function maxSizeFor(taskKind, maxProductSize, referenceSize) {
  if (taskKind === 'pangenome' && referenceSize > 0) {
    return Math.max(maxProductSize, Math.ceil(1.5 * referenceSize) + 500);
  }
  return maxProductSize;
}

function sitesSourceFor(sitesByPrimer, seq) {
  if (!sitesByPrimer) return null;
  if (typeof sitesByPrimer.view === 'function' && typeof sitesByPrimer.indexOf === 'function') {
    return sitesByPrimer.indexOf(seq) >= 0 ? { store: sitesByPrimer } : null;
  }
  if (sitesByPrimer instanceof Map) {
    const list = sitesByPrimer.get(seq);
    return list ? { list } : null;
  }
  const list = sitesByPrimer[seq];
  return Array.isArray(list) ? { list } : null;
}

// Candidates for one pair.
//   pair: { id?, left, right }
//   sitesByPrimer: SiteStore | Map<SEQ, site[]> | { SEQ: site[] }   (uppercase keys)
//   opts: { maxSize (required), maxCandidates = 5000, ignoreMismatches? (job filter) }
// Returns { candidates: [{ subject, orientation, start, end, size, fwd, rev }], truncated, sites_considered }.
// start/end are the pre-realignment 5' ends (fwd.p5, rev.p5); fwd/rev are site
// objects shared across calls on the same store (same key).
function candidates(pair, sitesByPrimer, opts) {
  const o = opts || {};
  if (!Number.isInteger(o.maxSize) || o.maxSize < 1) throw new TypeError('maxSize must be a positive integer');
  const maxCandidates = Number.isInteger(o.maxCandidates) && o.maxCandidates >= 0 ? o.maxCandidates : DEFAULT_MAX_CANDIDATES;
  const L = String(pair.left).toUpperCase();
  const R = String(pair.right).toUpperCase();
  const identical = L === R;
  const roles = identical ? [[L, 'L']] : [[L, 'L'], [R, 'R']];
  const hasFilter = Number.isInteger(o.ignoreMismatches);

  // Flat entry arrays: subject name, p5, primer length, role, lb, face and a site
  // accessor (store index or list object).
  const subjName = [];
  const p5 = [];
  const plen = [];
  const role = [];
  const lb = [];
  const isF = [];
  const getSite = [];

  for (const [seq, r] of roles) {
    const src = sitesSourceFor(sitesByPrimer, seq);
    if (!src) continue;
    if (src.store) {
      const store = src.store;
      const v = store.view(seq, hasFilter ? { ignoreMismatches: o.ignoreMismatches } : undefined);
      for (let k = 0; k < v.n; k++) {
        const i = v.idx[k];
        subjName.push(v.subjects[v.subj[i]]);
        p5.push(v.p5[i]);
        plen.push(v.len);
        role.push(r);
        lb.push(v.lb[k]);
        isF.push(v.strand[i] === 1);
        getSite.push({ store, q: v.q, i });
      }
    } else {
      for (const s of src.list) {
        if (hasFilter && !passesJobFilter(s, o.ignoreMismatches)) continue;
        if (!s.key) s.key = siteKey(seq, s.subject, s.face, s.p5);
        if (!s.primer) s.primer = seq;
        subjName.push(s.subject);
        p5.push(s.p5);
        plen.push(s.len || seq.length);
        role.push(r);
        lb.push(jobLowerBound(s));
        isF.push(s.face === 'F' || (s.face === undefined && s.strand === 1));
        getSite.push({ site: s });
      }
    }
  }

  const total = p5.length;
  const cmp = (a, b) => {
    if (subjName[a] !== subjName[b]) return subjName[a] < subjName[b] ? -1 : 1;
    if (p5[a] !== p5[b]) return p5[a] - p5[b];
    if (role[a] !== role[b]) return role[a] < role[b] ? -1 : 1;
    return lb[a] - lb[b];
  };
  const fIdx = [];
  const rIdx = [];
  for (let e = 0; e < total; e++) (isF[e] ? fIdx : rIdx).push(e);
  fIdx.sort(cmp);
  rIdx.sort(cmp);
  // Drop duplicate sites (same subject, p5 and role; the best lb sorts first).
  const dedupe = (arr) => arr.filter((e, k) => k === 0 ||
    !(subjName[arr[k - 1]] === subjName[e] && p5[arr[k - 1]] === p5[e] && role[arr[k - 1]] === role[e]));
  const F = dedupe(fIdx);
  const Rv = dedupe(rIdx);

  const cache = new Map();
  const materialize = (e) => {
    const g = getSite[e];
    if (g.site) return g.site;
    const id = g.q + ':' + g.i;
    let s = cache.get(id);
    if (!s) {
      s = g.store.site(g.q, g.i);
      cache.set(id, s);
    }
    return s;
  };

  const out = [];
  let truncated = false;
  const minLen = Math.min(L.length, R.length);
  let fk = 0;
  let rk = 0;
  outer:
  while (fk < F.length) {
    const subject = subjName[F[fk]];
    let fEnd = fk;
    while (fEnd < F.length && subjName[F[fEnd]] === subject) fEnd++;
    while (rk < Rv.length && subjName[Rv[rk]] < subject) rk++;
    let rEnd = rk;
    while (rEnd < Rv.length && subjName[Rv[rEnd]] === subject) rEnd++;
    let lo = rk;
    for (let a = fk; a < fEnd; a++) {
      const f = F[a];
      const fp = p5[f];
      // Smallest product with non-overlapping footprints: size = len(f) + len(r) >= 2 × minLen.
      while (lo < rEnd && p5[Rv[lo]] < fp + 2 * minLen - 1) lo++;
      for (let b = lo; b < rEnd; b++) {
        const r = Rv[b];
        const size = p5[r] - fp + 1;
        if (size > o.maxSize) break;
        if (size < plen[f] + plen[r]) continue;
        if (out.length >= maxCandidates) {
          truncated = true;
          break outer;
        }
        out.push({
          subject,
          orientation: identical ? 'LR' : role[f] + role[r],
          start: fp,
          end: p5[r],
          size,
          fwd: materialize(f),
          rev: materialize(r)
        });
      }
    }
    fk = fEnd;
    rk = rEnd;
  }
  return { candidates: out, truncated, sites_considered: F.length + Rv.length };
}

// Unique sites referenced by candidate lists (one genome), capped.
//   lists: candidate[] or candidate[][]
// Returns { sites, truncated }; candidates whose sites are not included are
// reported as `unaligned` by finalize.
function uniqueSites(lists, opts) {
  const maxSites = opts && Number.isInteger(opts.maxSites) ? opts.maxSites : DEFAULT_MAX_REALIGN_SITES;
  const arrays = lists.length && Array.isArray(lists[0]) ? lists : [lists];
  const seen = new Map();
  let truncated = false;
  for (const arr of arrays) {
    for (const c of arr) {
      for (const s of [c.fwd, c.rev]) {
        if (seen.has(s.key)) continue;
        if (seen.size >= maxSites) {
          truncated = true;
          continue;
        }
        seen.set(s.key, s);
      }
    }
  }
  return { sites: Array.from(seen.values()), truncated };
}

// Gap-aware lower bound on a site's re-aligned edits: BLAST mismatches + 1 per unaligned tail (sites.lbI).
// sites.jobLowerBound is not used here: for tails longer than 8 it is the substitution-only lbU, and a site
// whose unaligned tail is explained by one indel (e.g. P2_L q1-13 with t3 9 on a genome lacking a base: lbU 5,
// re-aligned edits 1) would otherwise be skipped and its product misclassified as unlikely.
function realignBound(s) {
  return lowerBoundI(s.mm, s.t5, s.t3);
}

// Unique sites of candidate lists split for genome re-alignment (run.js):
//   realign  sites whose gap-aware bound (realignBound) <= maxAmplifyingMismatches, capped at maxSites
//            (sites past the cap are left out → their candidates are unaligned → truncated)
//   bounded  sites whose bound already exceeds maxAmplifyingMismatches: no product using them can amplify, so
//            they get an approximate alignment instead of a FASTA re-alignment (not capped)
// opts: { maxSites = 20000, maxAmplifyingMismatches (omitted: every site is re-aligned) }
// Returns { realign, bounded, truncated }.
function realignPlan(lists, opts) {
  const o = opts || {};
  const maxSites = Number.isInteger(o.maxSites) ? o.maxSites : DEFAULT_MAX_REALIGN_SITES;
  const cap = Number.isInteger(o.maxAmplifyingMismatches) ? o.maxAmplifyingMismatches : Infinity;
  const arrays = lists.length && Array.isArray(lists[0]) ? lists : [lists];
  const seen = new Set();
  const realign = [];
  const bounded = [];
  let truncated = false;
  for (const arr of arrays) {
    for (const c of arr) {
      for (const s of [c.fwd, c.rev]) {
        if (seen.has(s.key)) continue;
        if (realignBound(s) > cap) {
          seen.add(s.key);
          bounded.push(s);
        } else if (realign.length >= maxSites) {
          truncated = true;
        } else {
          seen.add(s.key);
          realign.push(s);
        }
      }
    }
  }
  return { realign, bounded, truncated };
}

function orientationStrand(orientation) {
  if (orientation === 'LR') return 1;
  if (orientation === 'RL') return -1;
  return null;
}

// Candidates + alignments (Map key → aln, or function(site) → aln) → classified amplicons.
//   params: classify params (ignore_mismatches, min_total_mismatches, min_3p_mismatches, three_prime_window)
// Returns { amplicons: [likely|likely_weak], unlikely: [...], counts }
//   amplicon = { region, start, end, size, strand, orientation, likelihood,
//                left_mm, right_mm, left_3p_mm, right_3p_mm, left_mm_pos, right_mm_pos,
//                terminal_mismatch, approx }
// left_* describe the L-primer site and right_* the R-primer site; for LL/RR
// products left_* is the forward-facing site and right_* the reverse-facing one.
// Duplicates (same region, orientation, start, end) keep the best class, then
// the fewest mismatches. Lists are sorted by class, total mismatches, region, start.
// A candidate whose re-aligned ends cross, or whose primer footprints overlap (forward [fwd.p5, fwd.p3],
// reverse [rev.p3, rev.p5]; p3 falls back to p5 ± site length - 1), is invalid.
function finalize(cands, alignments, params) {
  const p = classify.withDefaults(params);
  const get = typeof alignments === 'function' ? alignments : (s) => alignments.get(s.key);
  const counts = { candidates: cands.length, likely: 0, likely_weak: 0, unlikely: 0, ignored: 0, unaligned: 0, invalid: 0, duplicates: 0 };
  const best = new Map();
  for (const c of cands) {
    const fa = get(c.fwd);
    const ra = get(c.rev);
    if (!fa || !ra || fa.mm == null || ra.mm == null) {
      counts.unaligned++;
      continue;
    }
    const start = fa.p5;
    const end = ra.p5;
    const fwdEnd = Number.isFinite(fa.p3) ? fa.p3 : start + ((c.fwd && c.fwd.len) || 0) - 1;
    const revStart = Number.isFinite(ra.p3) ? ra.p3 : end - ((c.rev && c.rev.len) || 0) + 1;
    if (!(end >= start) || fwdEnd >= revStart) {
      counts.invalid++;
      continue;
    }
    const likelihood = classify.classifyAmplicon(fa, ra, p);
    if (likelihood === 'ignored') {
      counts.ignored++;
      continue;
    }
    const leftIsFwd = c.orientation !== 'RL';
    const la = leftIsFwd ? fa : ra;
    const rb = leftIsFwd ? ra : fa;
    const amp = {
      region: c.subject,
      start,
      end,
      size: end - start + 1,
      strand: orientationStrand(c.orientation),
      orientation: c.orientation,
      likelihood,
      left_mm: la.mm,
      right_mm: rb.mm,
      left_3p_mm: classify.mismatches3p(la, p.three_prime_window),
      right_3p_mm: classify.mismatches3p(rb, p.three_prime_window),
      left_mm_pos: Array.isArray(la.mm_pos) ? la.mm_pos.slice() : null,
      right_mm_pos: Array.isArray(rb.mm_pos) ? rb.mm_pos.slice() : null,
      terminal_mismatch: !!(fa.terminal_mm || ra.terminal_mm),
      approx: !!(fa.approx || ra.approx)
    };
    const key = amp.region + '|' + amp.orientation + '|' + amp.start + '|' + amp.end;
    const prev = best.get(key);
    if (prev) {
      counts.duplicates++;
      if (compareAmplicons(amp, prev) < 0) best.set(key, amp);
    } else {
      best.set(key, amp);
    }
  }
  const amplicons = [];
  const unlikely = [];
  for (const amp of best.values()) {
    counts[amp.likelihood]++;
    (classify.countsAsAmplicon(amp.likelihood) ? amplicons : unlikely).push(amp);
  }
  amplicons.sort(compareAmplicons);
  unlikely.sort(compareAmplicons);
  return { amplicons, unlikely, counts };
}

function compareAmplicons(a, b) {
  const ra = classify.LIKELIHOOD_RANK[a.likelihood];
  const rb = classify.LIKELIHOOD_RANK[b.likelihood];
  if (ra !== rb) return ra - rb;
  const ma = a.left_mm + a.right_mm;
  const mb = b.left_mm + b.right_mm;
  if (ma !== mb) return ma - mb;
  if (a.region !== b.region) return a.region < b.region ? -1 : 1;
  if (a.start !== b.start) return a.start - b.start;
  if (a.end !== b.end) return a.end - b.end;
  return a.orientation < b.orientation ? -1 : a.orientation > b.orientation ? 1 : 0;
}

module.exports = {
  DEFAULT_MAX_CANDIDATES,
  DEFAULT_MAX_REALIGN_SITES,
  maxSizeFor,
  candidates,
  uniqueSites,
  realignBound,
  realignPlan,
  finalize,
  compareAmplicons,
  orientationStrand
};
