'use strict';

// Primer-BLAST style amplicon classification (spec §B.7, plan stringency decision).
//
//   blocked(s)  = s.mm >= min_total_mismatches && s.mm_3p >= min_3p_mismatches
//   ignored     : fwd.mm >= ignore_mismatches || rev.mm >= ignore_mismatches   (dropped)
//   unlikely    : fwd.mm > max_amplifying_mismatches || rev.mm > max_amplifying_mismatches
//                 || blocked(fwd) || blocked(rev)                              (counted; listed only if include_unlikely)
//   likely_weak : fwd.terminal_mm || rev.terminal_mm
//   likely      : otherwise
//   repetitive(primer) = near_perfect_sites > repeat_site_threshold
// A product amplifies (likely | likely_weak) only when each primer has <= max_amplifying_mismatches (default 3)
// edits and neither primer is 3'-blocked; 4-5 mismatch products are 'unlikely'.
//
// An alignment ("aln") is the object produced by realign.js:
//   { mm, mm_3p, terminal_mm, gaps, mm_pos: number[]|null, approx, ... }
// mm_pos holds distances from the primer 3' end (1 = terminal base). When it is
// present, mm_3p is recomputed for the requested three_prime_window so a single
// re-alignment serves any window.

const DEFAULT_PARAMS = Object.freeze({
  max_product_size: 4000,
  ignore_mismatches: 6,
  min_total_mismatches: 2,
  min_3p_mismatches: 2,
  three_prime_window: 5,
  include_unlikely: false,
  repeat_site_threshold: 5,
  max_amplifying_mismatches: 3
});

const LIKELIHOODS = Object.freeze(['likely', 'likely_weak', 'unlikely', 'ignored']);

// Lower rank is better ("best" likelihood in a group = minimum rank).
const LIKELIHOOD_RANK = Object.freeze({ likely: 0, likely_weak: 1, unlikely: 2, ignored: 3 });

function withDefaults(params) {
  return Object.assign({}, DEFAULT_PARAMS, params || {});
}

function mismatches3p(aln, threePrimeWindow) {
  if (Array.isArray(aln.mm_pos) && Number.isInteger(threePrimeWindow)) {
    let n = 0;
    for (const d of aln.mm_pos) if (d <= threePrimeWindow) n++;
    return n;
  }
  return aln.mm_3p || 0;
}

function isTerminalMismatch(aln) {
  return !!aln.terminal_mm;
}

function isBlocked(aln, params) {
  const p = withDefaults(params);
  return aln.mm >= p.min_total_mismatches && mismatches3p(aln, p.three_prime_window) >= p.min_3p_mismatches;
}

function isIgnored(aln, params) {
  const p = withDefaults(params);
  return aln.mm >= p.ignore_mismatches;
}

// More edits than a primer may have and still count as amplifying.
function isOverAmplifyingCap(aln, params) {
  const p = withDefaults(params);
  return aln.mm > p.max_amplifying_mismatches;
}

// fwd/rev are the alignments of the forward-facing and reverse-facing sites.
function classifyAmplicon(fwd, rev, params) {
  const p = withDefaults(params);
  if (isIgnored(fwd, p) || isIgnored(rev, p)) return 'ignored';
  if (isOverAmplifyingCap(fwd, p) || isOverAmplifyingCap(rev, p)) return 'unlikely';
  if (isBlocked(fwd, p) || isBlocked(rev, p)) return 'unlikely';
  if (isTerminalMismatch(fwd) || isTerminalMismatch(rev)) return 'likely_weak';
  return 'likely';
}

// likely and likely_weak both count as amplicons.
function countsAsAmplicon(likelihood) {
  return likelihood === 'likely' || likelihood === 'likely_weak';
}

function bestLikelihood(list) {
  let best = null;
  for (const l of list) {
    if (best === null || LIKELIHOOD_RANK[l] < LIKELIHOOD_RANK[best]) best = l;
  }
  return best;
}

function isRepetitive(nearPerfectSites, params) {
  const p = withDefaults(params);
  return nearPerfectSites > p.repeat_site_threshold;
}

module.exports = {
  DEFAULT_PARAMS,
  LIKELIHOODS,
  LIKELIHOOD_RANK,
  withDefaults,
  mismatches3p,
  isBlocked,
  isIgnored,
  isOverAmplifyingCap,
  classifyAmplicon,
  countsAsAmplicon,
  bestLikelihood,
  isRepetitive
};
