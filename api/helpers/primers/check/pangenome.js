'use strict';

// Pan-genome statuses, primary product and per-pair summary (spec §B.10, plan Overrides).
//
// | status          | rule (genome target: distinct amplicons; cDNA target: distinct gene groups) |
// | db_unavailable  | no dna DB (cdna DB in transcript mode), unresolved or ambiguous assembly     |
// | error           | BLAST failed twice or timed out                                              |
// | no_amplicon     | 0 likely; nearest = best unlikely product, if any                            |
// | single_perfect  | exactly 1 likely; left_mm == right_mm == 0                                   |
// | single_mismatch | exactly 1 likely; any mismatch; flags mismatch_in_3p_window, terminal_mismatch |
// | multiple        | >= 2 likely; primary + others[<= 10] + other_amplicons                       |
//
// Primary order: ortholog === true first, then left_mm+right_mm ascending, then |size - reference_size|
// ascending (cDNA groups: the isoform size closest to reference_size), then class, region/gene, start.
// size_delta = size - reference_size (null without a reference size).
// Summary: {genomes_total, single_perfect, single_mismatch, multiple, no_amplicon, db_unavailable, error,
//           amplifies, truncated}; amplifies = single_perfect + single_mismatch + multiple; status counts sum to
//           genomes_total. truncated = genomes whose search hit a site, candidate or re-alignment cap
//           (entry.truncated); it is not a status and not part of that sum. The status of a truncated genome
//           is a lower bound: products may be missing, so no_amplicon / single_* may understate what amplifies.
// maxSize for pan-genome tasks is amplicons.maxSizeFor('pangenome', max_product_size, reference_size)
//   = max(max_product_size, ceil(1.5 × reference_size) + 500).

const classify = require('./classify');
const { maxSizeFor } = require('./amplicons');
const { publicAmplicon, publicGroup } = require('./specificity');

const MAX_OTHERS = 10;
const STATUSES = Object.freeze(['single_perfect', 'single_mismatch', 'multiple', 'no_amplicon', 'db_unavailable', 'error']);
const AMPLIFYING = Object.freeze(['single_perfect', 'single_mismatch', 'multiple']);

function isGroup(p) {
  return Array.isArray(p.isoforms);
}

// Size used for ordering and size_delta.
function productSize(p, referenceSize) {
  if (!isGroup(p)) return p.size;
  if (!(referenceSize > 0)) return p.size_min;
  let best = p.isoforms[0].size;
  for (const iso of p.isoforms) {
    if (Math.abs(iso.size - referenceSize) < Math.abs(best - referenceSize)) best = iso.size;
  }
  return best;
}

function sizeDelta(p, referenceSize) {
  return referenceSize > 0 ? productSize(p, referenceSize) - referenceSize : null;
}

function comparePrimary(referenceSize) {
  return (a, b) => {
    const oa = a.ortholog === true ? 0 : 1;
    const ob = b.ortholog === true ? 0 : 1;
    if (oa !== ob) return oa - ob;
    const ma = a.left_mm + a.right_mm;
    const mb = b.left_mm + b.right_mm;
    if (ma !== mb) return ma - mb;
    if (referenceSize > 0) {
      const da = Math.abs(productSize(a, referenceSize) - referenceSize);
      const db = Math.abs(productSize(b, referenceSize) - referenceSize);
      if (da !== db) return da - db;
    }
    const ra = classify.LIKELIHOOD_RANK[a.likelihood];
    const rb = classify.LIKELIHOOD_RANK[b.likelihood];
    if (ra !== rb) return ra - rb;
    const ka = isGroup(a) ? a.gene_id : a.region;
    const kb = isGroup(b) ? b.gene_id : b.region;
    if (ka !== kb) return ka < kb ? -1 : 1;
    const sa = isGroup(a) ? a.size_min : a.start;
    const sb = isGroup(b) ? b.size_min : b.start;
    if (sa !== sb) return sa - sb;
    return a.orientation < b.orientation ? -1 : a.orientation > b.orientation ? 1 : 0;
  };
}

// Sets p.ortholog on products: genome amplicons → any overlapping gene in the ortholog set
// (null when genes are unknown); cDNA groups → gene_id in the set. null without a design gene.
function markOrthologs(products, orthologIds) {
  for (const p of products || []) {
    if (!orthologIds) p.ortholog = null;
    else if (isGroup(p)) p.ortholog = orthologIds.has(p.gene_id);
    else if (Array.isArray(p.genes)) p.ortholog = p.genes.some((g) => orthologIds.has(String(g.id)));
    else p.ortholog = null;
  }
  return products;
}

// Public product with the pan-genome decorations.
function decorate(p, referenceSize) {
  const base = isGroup(p) ? publicGroup(p) : publicAmplicon(p);
  if (isGroup(p)) base.size = productSize(p, referenceSize);
  base.size_delta = sizeDelta(p, referenceSize);
  base.mismatch_in_3p_window = (p.left_3p_mm || 0) > 0 || (p.right_3p_mm || 0) > 0;
  base.terminal_mismatch = !!p.terminal_mismatch;
  base.ortholog = p.ortholog === undefined ? null : p.ortholog;
  return base;
}

function statusFor(products) {
  if (products.length === 0) return 'no_amplicon';
  if (products.length > 1) return 'multiple';
  const p = products[0];
  return p.left_mm === 0 && p.right_mm === 0 ? 'single_perfect' : 'single_mismatch';
}

// One genome × pair entry.
//   input: { system_name, display_name, products (likely amplicons or groups), unlikely,
//            referenceSize, orthologAnnotated (bool|null), truncated }
// Products should already carry `ortholog` (markOrthologs); missing values count as null.
function genomeEntry(input) {
  const o = input || {};
  const products = (o.products || []).slice();
  const cmp = comparePrimary(o.referenceSize);
  products.sort(cmp);
  const status = statusFor(products);
  let nearest = null;
  if (status === 'no_amplicon' && o.unlikely && o.unlikely.length) {
    nearest = decorate(o.unlikely.slice().sort(cmp)[0], o.referenceSize);
  }
  return {
    system_name: o.system_name,
    display_name: o.display_name != null ? o.display_name : o.system_name,
    status,
    ortholog_annotated: o.orthologAnnotated === undefined ? null : o.orthologAnnotated,
    truncated: !!o.truncated,
    primary: products.length ? decorate(products[0], o.referenceSize) : null,
    other_amplicons: Math.max(0, products.length - 1),
    others: products.slice(1, 1 + MAX_OTHERS).map((p) => decorate(p, o.referenceSize)),
    nearest
  };
}

function emptyEntry(input, status) {
  const o = input || {};
  return {
    system_name: o.system_name,
    display_name: o.display_name != null ? o.display_name : o.system_name,
    status,
    ortholog_annotated: o.orthologAnnotated === undefined ? null : o.orthologAnnotated,
    truncated: false,
    primary: null,
    other_amplicons: 0,
    others: [],
    nearest: null
  };
}

// input: { system_name, display_name, orthologAnnotated, reason }
function unavailableEntry(input) {
  const e = emptyEntry(input, 'db_unavailable');
  e.reason = input && input.reason ? input.reason : 'NO_BLASTDB';
  return e;
}

// input: { system_name, display_name, orthologAnnotated, error: {code, message} }
function errorEntry(input) {
  const e = emptyEntry(input, 'error');
  const err = input && input.error;
  e.error = { code: err && err.code ? String(err.code) : 'BLAST_FAILED', message: err && err.message ? String(err.message) : 'BLAST failed' };
  return e;
}

function summarize(entries) {
  const s = { genomes_total: 0, single_perfect: 0, single_mismatch: 0, multiple: 0, no_amplicon: 0, db_unavailable: 0, error: 0, amplifies: 0, truncated: 0 };
  for (const e of entries || []) {
    if (!e) continue;
    if (STATUSES.indexOf(e.status) < 0) throw new TypeError('unknown pan-genome status ' + e.status);
    s.genomes_total++;
    s[e.status]++;
    if (e.truncated === true) s.truncated++;
  }
  s.amplifies = s.single_perfect + s.single_mismatch + s.multiple;
  return s;
}

function pangenomeMaxSize(maxProductSize, referenceSize) {
  return maxSizeFor('pangenome', maxProductSize, referenceSize);
}

module.exports = {
  MAX_OTHERS,
  STATUSES,
  AMPLIFYING,
  productSize,
  sizeDelta,
  comparePrimary,
  markOrthologs,
  decorate,
  statusFor,
  genomeEntry,
  unavailableEntry,
  errorEntry,
  summarize,
  pangenomeMaxSize
};
