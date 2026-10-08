'use strict';

// Small builders for verdict / pan-genome unit tests. Shapes follow amplicons.finalize output
// (genome amplicons) and annotate.groupByGene output (cDNA groups).

const W = 5;

function count3p(pos, w) {
  if (!Array.isArray(pos)) return 0;
  return pos.filter((d) => d <= (w || W)).length;
}

// makeAmplicon({region, start, end, orientation, likelihood, left_mm, right_mm, left_mm_pos, right_mm_pos,
//               terminal_mismatch, approx, genes})
// Mismatch positions default to [] for 0 mismatches; a count without positions leaves them null.
function makeAmplicon(o) {
  const a = o || {};
  const region = a.region != null ? String(a.region) : '4';
  const start = a.start != null ? a.start : 7423537;
  const end = a.end != null ? a.end : 7423746;
  const orientation = a.orientation || 'LR';
  const lpos = a.left_mm_pos !== undefined ? a.left_mm_pos : (a.left_mm ? null : []);
  const rpos = a.right_mm_pos !== undefined ? a.right_mm_pos : (a.right_mm ? null : []);
  const leftMm = a.left_mm != null ? a.left_mm : (Array.isArray(lpos) ? lpos.length : 0);
  const rightMm = a.right_mm != null ? a.right_mm : (Array.isArray(rpos) ? rpos.length : 0);
  const amp = {
    region,
    start,
    end,
    size: end - start + 1,
    strand: orientation === 'LR' ? 1 : orientation === 'RL' ? -1 : null,
    orientation,
    likelihood: a.likelihood || 'likely',
    left_mm: leftMm,
    right_mm: rightMm,
    left_3p_mm: a.left_3p_mm != null ? a.left_3p_mm : count3p(lpos),
    right_3p_mm: a.right_3p_mm != null ? a.right_3p_mm : count3p(rpos),
    left_mm_pos: lpos,
    right_mm_pos: rpos,
    terminal_mismatch: !!a.terminal_mismatch,
    approx: !!a.approx
  };
  if (a.genes !== undefined) amp.genes = a.genes;
  return amp;
}

// makeGroup({gene_id, orientation, sizes: [..] | isoforms: [{transcript_id, size}], likelihood, left_mm, ...})
function makeGroup(o) {
  const g = o || {};
  const geneId = g.gene_id || 'GENE';
  const isoforms = g.isoforms || (g.sizes || [278]).map((size, i) => ({ transcript_id: geneId + '.' + (i + 1), start: 100, end: 100 + size - 1, size, likelihood: g.likelihood || 'likely', left_mm: g.left_mm || 0, right_mm: g.right_mm || 0 }));
  const sizes = isoforms.map((i) => i.size);
  const lpos = g.left_mm_pos !== undefined ? g.left_mm_pos : (g.left_mm ? null : []);
  const rpos = g.right_mm_pos !== undefined ? g.right_mm_pos : (g.right_mm ? null : []);
  return {
    gene_id: geneId,
    orientation: g.orientation || 'LR',
    isoforms,
    size_min: Math.min.apply(null, sizes),
    size_max: Math.max.apply(null, sizes),
    likelihood: g.likelihood || 'likely',
    left_mm: g.left_mm || 0,
    right_mm: g.right_mm || 0,
    left_3p_mm: g.left_3p_mm != null ? g.left_3p_mm : count3p(lpos),
    right_3p_mm: g.right_3p_mm != null ? g.right_3p_mm : count3p(rpos),
    left_mm_pos: lpos,
    right_mm_pos: rpos,
    terminal_mismatch: !!g.terminal_mismatch,
    approx: !!g.approx
  };
}

module.exports = { makeAmplicon, makeGroup };
