'use strict';

// BLAST sensitivity guarantee for primer sites (spec §B.3).
//
// blastn-short with -ungapped reports a substitution-only site when it contains
// an exact word of length >= W and its ungapped score (reward 1, penalty -1) is
// about 11 or more (score 12 was reported on squam, 10 was not). A site with k
// substitutions in an L-mer splits the L-k matching bases into k+1 runs, so the
// longest run is at least ceil((L-k)/(k+1)), and the full-length ungapped score
// is L-2k. The largest k that satisfies both bounds is guaranteed to be found.

const MIN_UNGAPPED_SCORE = 12;
const DEFAULT_MAX_AMPLIFYING_MISMATCHES = 3;
const DEFAULT_REFERENCE_WORD_SIZE = 5;
const DEFAULT_PRIMER_LENGTHS = Object.freeze([20, 22]);

function assertPositiveInt(name, v) {
  if (!Number.isInteger(v) || v < 1) {
    throw new TypeError(name + ' must be a positive integer, got ' + v);
  }
}

// Largest k >= 0 with ceil((L-k)/(k+1)) >= W and L-2k >= 12. Both bounds shrink
// as k grows, so the first failing k ends the scan. Returns 0 when even a
// perfect site is not guaranteed (L < W or L < 12), which cannot happen for the
// 15-36 nt primers the API accepts.
function guaranteedMaxMismatches(L, W) {
  assertPositiveInt('L', L);
  assertPositiveInt('W', W);
  let best = 0;
  for (let k = 0; k <= L; k++) {
    const longestRun = Math.ceil((L - k) / (k + 1));
    if (longestRun < W || L - 2 * k < MIN_UNGAPPED_SCORE) break;
    best = k;
  }
  return best;
}

// Per-primer block for results.primers[seq].sensitivity (§B.12).
// opts: { referenceWordSize, pangenomeWordSize? }. The pangenome key is only
// present when a pangenome word size is given.
function primerSensitivity(len, opts) {
  const o = opts || {};
  const out = {};
  if (o.referenceWordSize != null) {
    out.reference = { word_size: o.referenceWordSize, guaranteed_max_mismatches: guaranteedMaxMismatches(len, o.referenceWordSize) };
  }
  if (o.pangenomeWordSize != null) {
    out.pangenome = { word_size: o.pangenomeWordSize, guaranteed_max_mismatches: guaranteedMaxMismatches(len, o.pangenomeWordSize) };
  }
  return out;
}

function span(a, b) {
  return a === b ? String(a) : a + '-' + b;
}

function lengthRange(lengths) {
  const ls = (Array.isArray(lengths) ? lengths : []).filter(function (n) { return Number.isInteger(n) && n > 0; });
  if (!ls.length) return DEFAULT_PRIMER_LENGTHS.slice();
  return [Math.min.apply(null, ls), Math.max.apply(null, ls)];
}

// One sentence about what a search at word size ws is guaranteed to find for primers of minL..maxL nt.
function searchSentence(name, ws, minL, maxL, cap) {
  const gMin = guaranteedMaxMismatches(minL, ws);
  const gMax = guaranteedMaxMismatches(maxL, ws);
  let s = 'The ' + name + ' search (word size ' + ws + ') finds every site with up to ' + span(gMin, gMax) +
    ' mismatches for these ' + span(minL, maxL) + ' nt primers';
  if (gMin >= cap) return s + ', which covers the cap.';
  const missed = gMin + 1 === cap ? String(cap) : (gMin + 1) + '-' + cap;
  return s + ', so products with ' + missed + ' mismatches per primer can be missed' +
    (name === 'pan-genome' ? ' in other genomes' : '') + (gMax >= cap ? ' for the shorter primers' : '') + '.';
}

// Request-specific text for results.sensitivity_note.
// opts: { maxAmplifyingMismatches (default 3), referenceWordSize (default 5),
//         pangenomeWordSize (omit when no pan-genome search runs),
//         primerLengths: lengths of the request's primers (default 20-22 nt) }
function sensitivityNote(opts) {
  const o = opts || {};
  const cap = Number.isInteger(o.maxAmplifyingMismatches) && o.maxAmplifyingMismatches >= 0 ?
    o.maxAmplifyingMismatches : DEFAULT_MAX_AMPLIFYING_MISMATCHES;
  const wsRef = Number.isInteger(o.referenceWordSize) && o.referenceWordSize > 0 ? o.referenceWordSize : DEFAULT_REFERENCE_WORD_SIZE;
  const wsPan = Number.isInteger(o.pangenomeWordSize) && o.pangenomeWordSize > 0 ? o.pangenomeWordSize : null;
  const range = lengthRange(o.primerLengths);
  const parts = [
    'Sites are detected only if they contain an exact match of at least the word size and an ungapped score of about 11 or more.',
    'An indel is found when the primer part on at least one side of it scores that much (both sides are joined into one ' +
      'gapped site); indels near the middle of primers shorter than about 22 nt may be missed.',
    'A product amplifies only when each primer has at most ' + cap + ' mismatch' + (cap === 1 ? '' : 'es') +
      ' (max_amplifying_mismatches).',
    searchSentence('reference', wsRef, range[0], range[1], cap)
  ];
  if (wsPan !== null) parts.push(searchSentence('pan-genome', wsPan, range[0], range[1], cap));
  parts.push('See primers[].sensitivity.');
  return parts.join(' ');
}

// The note for a default request (cap 3, reference word size 5, 20-22 nt primers, no pan-genome search).
const SENSITIVITY_NOTE = sensitivityNote();

module.exports = {
  MIN_UNGAPPED_SCORE,
  DEFAULT_MAX_AMPLIFYING_MISMATCHES,
  SENSITIVITY_NOTE,
  guaranteedMaxMismatches,
  primerSensitivity,
  sensitivityNote
};
