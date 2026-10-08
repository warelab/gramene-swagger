'use strict';

// Pan-genome products for the allele-caller tests (genotyping spec §5.6-§5.7, §7.3), built the way check/run.js builds them:
// every primer site is placed like a BLAST hit (the best ungapped placement of the whole primer), re-aligned with realign.js
// over its footprint ± 3 bp, paired into one LR/RL candidate per pair and classified by amplicons.finalize. The genome is an
// in-memory plus-strand string, so the caller's fetch can read it back.
//
//   genome(region, seq, start = 1)                   -> {region, seq, start, end}
//   productsForSet(genome, set, strand, opts)        -> {ref: {amplicons, unlikely}, alt: {amplicons, unlikely}}
//        one copy of the set's locus on `strand` inside genome (opts.lo/opts.hi limit the search to a window of the genome)
//   mergeProducts(a, b)                              -> the per-pair lists of two copies concatenated and re-sorted
//   fetchFrom(genomes)                               -> deps {fetch(region, start, end), regionLength(region)}
//   prepared(body, cfg)                              -> check/genotype.js prepare() over the recorded sorghum_bicolor window

const path = require('path');

const ROOT = path.resolve(__dirname, '..', '..', '..', '..', '..');
const realign = require(path.join(ROOT, 'api/helpers/primers/check/realign'));
const amplicons = require(path.join(ROOT, 'api/helpers/primers/check/amplicons'));
const classify = require(path.join(ROOT, 'api/helpers/primers/check/classify'));
const genotype = require(path.join(ROOT, 'api/helpers/primers/check/genotype'));
const variation = require(path.join(ROOT, 'api/helpers/primers/variation/normalize'));
const stubs = require('./stubs');

const PAD = realign.PAD;

function genome(region, seq, start) {
  const s = seq.toUpperCase();
  const first = start || 1;
  return { region: String(region), seq: s, start: first, end: first + s.length - 1 };
}

// Best ungapped placement (fewest mismatches, then the leftmost) of primer on the genome window [lo, hi], facing F (the primer
// reads the plus strand) or R; realign.js then re-aligns it over the placed footprint ± 3 bp, as run.js does for a BLAST hit.
function place(g, primer, face, lo, hi) {
  const p = primer.toUpperCase();
  const plus = face === 'F' ? p : realign.revcomp(p);
  const n = p.length;
  let bestAt = -1;
  let bestMm = Infinity;
  for (let x = lo; x + n - 1 <= hi; x++) {
    let mm = 0;
    for (let k = 0; k < n && mm < bestMm; k++) if (g.seq.charCodeAt(x - g.start + k) !== plus.charCodeAt(k)) mm++;
    if (mm < bestMm) {
      bestMm = mm;
      bestAt = x;
    }
  }
  if (bestAt < 0) throw new Error('products: ' + primer + ' does not fit in ' + g.region + ':' + lo + '-' + hi);
  const p5 = face === 'F' ? bestAt : bestAt + n - 1;
  const p3 = face === 'F' ? bestAt + n - 1 : bestAt;
  const wlo = Math.max(g.start, Math.min(p5, p3) - PAD);
  const whi = Math.min(g.end, Math.max(p5, p3) + PAD);
  const site = { key: p + ':' + g.region + ':' + face + ':' + p5, primer: p, len: n, subject: g.region, face: face, p5: p5, p3: p3 };
  const aln = realign.realignSiteOnWindow(site, g.seq.slice(wlo - g.start, whi - g.start + 1), wlo, classify.DEFAULT_PARAMS.three_prime_window, g.end);
  aln.key = site.key;
  return { site: site, aln: aln, ungapped_mm: bestMm };
}

// One pair's product on one copy: in the reference orientation the left primer faces F and the right primer R; on a copy on the
// minus strand both faces flip and the product is RL.
function pairProducts(g, left, right, strand, lo, hi, params) {
  const L = place(g, left, strand === 1 ? 'F' : 'R', lo, hi);
  const R = place(g, right, strand === 1 ? 'R' : 'F', lo, hi);
  const fwd = strand === 1 ? L : R;
  const rev = strand === 1 ? R : L;
  const cand = {
    subject: g.region, orientation: strand === 1 ? 'LR' : 'RL', start: fwd.site.p5, end: rev.site.p5,
    size: rev.site.p5 - fwd.site.p5 + 1, fwd: fwd.site, rev: rev.site
  };
  const fin = amplicons.finalize([cand], new Map([[fwd.site.key, fwd.aln], [rev.site.key, rev.aln]]), params || classify.DEFAULT_PARAMS);
  return { amplicons: fin.amplicons, unlikely: fin.unlikely };
}

// set: a prepare() set (orientation and primers). opts: {lo, hi, params}
function productsForSet(g, set, strand, opts) {
  const o = opts || {};
  const lo = o.lo || g.start;
  const hi = o.hi || g.end;
  const forward = set.orientation === 'forward';
  const pr = set.primers;
  const pair = function (as) { return forward ? [as, pr.common] : [pr.common, as]; };
  const ref = pair(pr.as_ref);
  const alt = pair(pr.as_alt);
  return {
    ref: pairProducts(g, ref[0], ref[1], strand, lo, hi, o.params),
    alt: pairProducts(g, alt[0], alt[1], strand, lo, hi, o.params)
  };
}

function mergeProducts(a, b) {
  const both = function (x, y) { return (x || []).concat(y || []).sort(amplicons.compareAmplicons); };
  return {
    ref: { amplicons: both(a.ref.amplicons, b.ref.amplicons), unlikely: both(a.ref.unlikely, b.ref.unlikely) },
    alt: { amplicons: both(a.alt.amplicons, b.alt.amplicons), unlikely: both(a.alt.unlikely, b.alt.unlikely) }
  };
}

// deps.fetch / deps.regionLength over in-memory genomes (each genome's region starts at its start coordinate).
// calls records [region, start, end].
function fetchFrom(genomes) {
  const byRegion = new Map(genomes.map(function (g) { return [g.region, g]; }));
  const calls = [];
  return {
    calls: calls,
    fetch: async function (region, start, end) {
      calls.push([region, start, end]);
      const g = byRegion.get(String(region));
      if (!g || start < g.start || end > g.end) throw Object.assign(new Error('no bases for ' + region + ':' + start + '-' + end), { code: 'REGION_OUT_OF_BOUNDS' });
      return g.seq.slice(start - g.start, end - g.start + 1);
    },
    regionLength: async function (region) {
      const g = byRegion.get(String(region));
      return g && g.start === 1 ? g.end : undefined;
    }
  };
}

// prepare() over the recorded sorghum_bicolor 1:9000-15500 window for a check body of stubs.js.
function prepared(body, cfg) {
  const v = body.genotyping.variant;
  const len = stubs.regionLengthOf('sorghum_bicolor', v.region);
  const win = genotype.windowFor({ position: v.position, ref: v.ref }, len);
  const start = Math.max(win.start, stubs.WINDOW.start);
  const end = Math.min(win.end, stubs.WINDOW.end);
  return genotype.prepare(body.genotyping, body.pairs, variation.sequenceWindow(stubs.bases(start, end), start, len), { cfg: cfg || stubs.makeCfg() });
}

module.exports = { genome, place, pairProducts, productsForSet, mergeProducts, fetchFrom, prepared };
