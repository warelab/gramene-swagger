'use strict';

// Verdicts by mode × target (spec §B.8). Pure functions over amplicons.finalize output
// (genome target) and annotate.groupByGene output (cDNA target).
//
// | mode                              | target | on-target                                              | verdicts |
// | gene, region with expected        | genome | likely LR/RL product with region == E.region,          | specific, off_targets, on_target_missing, truncated, error |
// |                                   |        | |start-E.start| <= 2 and |end-E.end| <= 2               |          |
// | gene/region without expected,     | genome | exactly one likely LR/RL product with left_mm == right_mm == 0 | specific, off_targets, unverified_target, truncated, error |
// | sequence                          |        | → on_target_inferred: true                             |          |
// | transcript                        | genome | none; likely products inside the gene span → gdna_products | specific, off_targets, truncated, error (never on_target_missing) |
// | transcript                        | cdna   | best likely group of gene_id                           | specific, off_targets, on_target_missing, truncated, error |
//
// Verdict precedence: no on-target (where one is required) → truncated ? 'truncated' : missing verdict;
// otherwise off-targets present → 'off_targets'; otherwise truncated ? 'truncated' : 'specific'.
// off_target_count is exact; listed off_targets / gdna_products / unlikely are capped (100).
// Without an on-target (unverified_target) every likely product is listed as an off-target.

const ON_TARGET_TOLERANCE = 2;
const DEFAULT_MAX_LISTED = 100;
const VERDICTS = Object.freeze(['specific', 'off_targets', 'on_target_missing', 'unverified_target', 'truncated', 'error']);

function isPairedOrientation(orientation) {
  return orientation === 'LR' || orientation === 'RL';
}

// Public genome amplicon (§B.12 specificity entries).
function publicAmplicon(a) {
  return {
    region: a.region,
    start: a.start,
    end: a.end,
    size: a.size,
    strand: a.strand != null ? a.strand : null,
    orientation: a.orientation,
    likelihood: a.likelihood,
    left_mm: a.left_mm,
    right_mm: a.right_mm,
    left_3p_mm: a.left_3p_mm,
    right_3p_mm: a.right_3p_mm,
    left_mm_pos: a.left_mm_pos != null ? a.left_mm_pos : null,
    right_mm_pos: a.right_mm_pos != null ? a.right_mm_pos : null,
    terminal_mismatch: !!a.terminal_mismatch,
    approx: !!a.approx,
    genes: a.genes !== undefined ? a.genes : null
  };
}

// Public cDNA gene group (annotate.groupByGene output is already public; copy it).
function publicGroup(g) {
  return {
    gene_id: g.gene_id,
    orientation: g.orientation,
    isoforms: g.isoforms.map((i) => Object.assign({}, i)),
    size_min: g.size_min,
    size_max: g.size_max,
    likelihood: g.likelihood,
    left_mm: g.left_mm,
    right_mm: g.right_mm,
    left_3p_mm: g.left_3p_mm,
    right_3p_mm: g.right_3p_mm,
    left_mm_pos: g.left_mm_pos != null ? g.left_mm_pos : null,
    right_mm_pos: g.right_mm_pos != null ? g.right_mm_pos : null,
    terminal_mismatch: !!g.terminal_mismatch,
    approx: !!g.approx
  };
}

function matchesExpected(a, expected, tolerance) {
  const tol = Number.isInteger(tolerance) ? tolerance : ON_TARGET_TOLERANCE;
  return !!expected && isPairedOrientation(a.orientation) && String(a.region) === String(expected.region) &&
    Math.abs(a.start - expected.start) <= tol && Math.abs(a.end - expected.end) <= tol;
}

function insideSpan(a, loc) {
  return !!loc && String(a.region) === String(loc.region) && a.start >= loc.start && a.end <= loc.end;
}

// `expected` is honoured only in gene and region modes (§A.2.3).
function honoredExpected(mode, pair) {
  return (mode === 'gene' || mode === 'region') && pair && pair.expected ? pair.expected : null;
}

function verdictFor(onTargetRequired, onTarget, offCount, truncated, missingVerdict) {
  if (onTargetRequired && !onTarget) return truncated ? 'truncated' : missingVerdict;
  if (offCount > 0) return 'off_targets';
  return truncated ? 'truncated' : 'specific';
}

// Genome-target selection (raw amplicon objects, so the caller can annotate before publishing).
//   input: { amplicons (likely + likely_weak, sorted best first), unlikely, truncated, mode,
//            expected (already honoured or null), geneLocation ({region,start,end}) }
//   → { verdict, onTarget, inferred, off: [], gdna: [], unlikely: [], truncated }
function selectGenome(input) {
  const o = input || {};
  const amps = Array.isArray(o.amplicons) ? o.amplicons : [];
  const unlikely = Array.isArray(o.unlikely) ? o.unlikely : [];
  const truncated = !!o.truncated;
  let onTarget = null;
  let inferred = false;
  let off = [];
  let gdna = [];
  let verdict;
  if (o.mode === 'transcript') {
    for (const a of amps) (insideSpan(a, o.geneLocation) ? gdna : off).push(a);
    verdict = verdictFor(false, null, off.length, truncated);
  } else if (o.expected) {
    const matches = amps.filter((a) => matchesExpected(a, o.expected));
    onTarget = matches.length ? matches[0] : null;
    off = amps.filter((a) => matches.indexOf(a) < 0);
    verdict = verdictFor(true, onTarget, off.length, truncated, 'on_target_missing');
  } else {
    const perfect = amps.filter((a) => isPairedOrientation(a.orientation) && a.left_mm === 0 && a.right_mm === 0);
    if (perfect.length === 1) {
      onTarget = perfect[0];
      inferred = true;
      off = amps.filter((a) => a !== onTarget);
    } else {
      off = amps.slice();
    }
    verdict = verdictFor(true, onTarget, off.length, truncated, 'unverified_target');
  }
  return { verdict, onTarget, inferred, off, gdna, unlikely, truncated };
}

// Amplicons worth annotating for one selection, in priority order: on-target, off-targets and
// gDNA products up to the listing cap, then unlikely ones when they are listed.
function annotationCandidates(sel, opts) {
  const max = opts && Number.isInteger(opts.maxListed) ? opts.maxListed : DEFAULT_MAX_LISTED;
  const out = [];
  if (sel.onTarget) out.push(sel.onTarget);
  const lists = [sel.off.slice(0, max), sel.gdna.slice(0, max)];
  if (opts && opts.includeUnlikely) lists.push(sel.unlikely.slice(0, max));
  return { head: out, lists };
}

// Public §B.12 specificity pair block.
//   opts: { id, maxListed = 100, includeUnlikely = false }
function genomePairBlock(sel, opts) {
  const o = opts || {};
  const max = Number.isInteger(o.maxListed) ? o.maxListed : DEFAULT_MAX_LISTED;
  const block = {
    id: o.id,
    verdict: sel.verdict,
    on_target_inferred: !!sel.inferred,
    truncated: !!sel.truncated,
    on_target: sel.onTarget ? publicAmplicon(sel.onTarget) : null,
    off_target_count: sel.off.length,
    unlikely_count: sel.unlikely.length,
    off_targets: sel.off.slice(0, max).map(publicAmplicon),
    gdna_products: sel.gdna.slice(0, max).map(publicAmplicon)
  };
  if (o.includeUnlikely) block.unlikely = sel.unlikely.slice(0, max).map(publicAmplicon);
  return block;
}

// cDNA-target selection over gene groups (annotate.groupByGene, sorted best first).
//   input: { groups, unlikelyGroups, truncated, geneId }
//   → { verdict, onTarget, off, sameGene (other groups of gene_id), unlikely, truncated }
// The on-target is the best LR/RL group of gene_id (else its best group). Further groups of the
// design gene (another orientation) are not off-targets; they are counted in other_on_target_groups.
function selectCdna(input) {
  const o = input || {};
  const groups = Array.isArray(o.groups) ? o.groups : [];
  const unlikely = Array.isArray(o.unlikelyGroups) ? o.unlikelyGroups : [];
  const truncated = !!o.truncated;
  const mine = o.geneId != null ? groups.filter((g) => g.gene_id === String(o.geneId)) : [];
  const onTarget = mine.find((g) => isPairedOrientation(g.orientation)) || mine[0] || null;
  const off = groups.filter((g) => mine.indexOf(g) < 0);
  const verdict = verdictFor(true, onTarget, off.length, truncated, 'on_target_missing');
  return { verdict, onTarget, off, sameGene: mine.filter((g) => g !== onTarget), unlikely, truncated };
}

// Public transcriptome pair block:
//   { id, verdict, truncated, on_target: group|null, other_on_target_groups, off_target_count,
//     unlikely_count, off_targets: [group] (≤ 100), unlikely?: [group] }
function cdnaPairBlock(sel, opts) {
  const o = opts || {};
  const max = Number.isInteger(o.maxListed) ? o.maxListed : DEFAULT_MAX_LISTED;
  const block = {
    id: o.id,
    verdict: sel.verdict,
    truncated: !!sel.truncated,
    on_target: sel.onTarget ? publicGroup(sel.onTarget) : null,
    other_on_target_groups: sel.sameGene.length,
    off_target_count: sel.off.length,
    unlikely_count: sel.unlikely.length,
    off_targets: sel.off.slice(0, max).map(publicGroup)
  };
  if (o.includeUnlikely) block.unlikely = sel.unlikely.slice(0, max).map(publicGroup);
  return block;
}

// Block for a pair whose target could not be checked (verdict 'error').
function errorPairBlock(target, id, error) {
  const block = {
    id,
    verdict: 'error',
    truncated: false,
    on_target: null,
    off_target_count: 0,
    unlikely_count: 0,
    off_targets: [],
    error: error ? { code: error.code || 'CHECK_FAILED', message: error.message || String(error.code || 'error') } : null
  };
  if (target === 'genome') {
    block.on_target_inferred = false;
    block.gdna_products = [];
  } else {
    block.other_on_target_groups = 0;
  }
  return block;
}

// reference_size for the pan-genome (§B.2): size of the on-target (or inferred) reference product.
// Genome blocks use on_target.size; cDNA blocks use the isoform of transcriptId when present,
// else the group's size_min. null without an on-target.
function referenceSize(block, opts) {
  if (!block || !block.on_target) return null;
  const t = block.on_target;
  if (Array.isArray(t.isoforms)) {
    const tx = opts && opts.transcriptId;
    const iso = tx ? t.isoforms.find((i) => i.transcript_id === tx) : null;
    return iso ? iso.size : t.size_min;
  }
  return Number.isInteger(t.size) ? t.size : null;
}

module.exports = {
  ON_TARGET_TOLERANCE,
  DEFAULT_MAX_LISTED,
  VERDICTS,
  isPairedOrientation,
  publicAmplicon,
  publicGroup,
  matchesExpected,
  insideSpan,
  honoredExpected,
  selectGenome,
  annotationCandidates,
  genomePairBlock,
  selectCdna,
  cdnaPairBlock,
  errorPairBlock,
  referenceSize
};
