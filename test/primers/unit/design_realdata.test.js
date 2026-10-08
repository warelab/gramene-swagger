'use strict';

// Real-data design checks (spec §10.3 rows 1, 2, 4, 5, 11, 12), in process, no server.
//   PRIMERS_REALDATA=1 node --test test/primers/unit/design_realdata.test.js
// Needs mongo (sorghum11), /scratch/olson/fasta, /home/olson/bin/primer3_core and blastn, and the
// fastaIdx service on localhost:8888 (PRIMERS_FASTAIDX) for the genomic block cross-check.
// Runs at most one megablast at a time (2 threads, nice 10).

require('../../../api/helpers/primers/node_compat');

const { test, after } = require('node:test');
const should = require('should');

const REALDATA = process.env.PRIMERS_REALDATA === '1';
const SKIP = REALDATA ? false : 'set PRIMERS_REALDATA=1';
const FASTAIDX = process.env.PRIMERS_FASTAIDX || 'http://localhost:8888';

const design = require('../../../api/helpers/primers/design');
const coords = require('../../../api/helpers/primers/coords');
const boulder = require('../../../api/helpers/primers/boulder');
const stubs = require('../fixtures/design/stubs');

const G200 = 'SORBI_3001G000200';
const G700 = 'SORBI_3001G000700';
const G087700 = 'SORBI_3004G087700';

const latencies = [];

async function timed(label, body) {
  const t0 = process.hrtime.bigint();
  const r = await design.design(body);
  latencies.push({ case: label, ms: Math.round(Number(process.hrtime.bigint() - t0) / 1e6), pairs: r.pairs.length });
  return r;
}

async function rejectsCode(promise, code, status) {
  let err = null;
  try {
    await promise;
  } catch (e) {
    err = e;
  }
  should.exist(err, 'expected ' + code);
  err.should.have.property('code', code);
  if (status) err.should.have.property('status', status);
  return err;
}

const regionCache = new Map();

async function fastaIdx(systemName, region, start, end, strand) {
  const key = [systemName, region, start, end, strand].join(':');
  if (!regionCache.has(key)) {
    const url = FASTAIDX + '/sequence/region/' + systemName + '/' + region + ':' + start + '..' + end + ':' + strand;
    const res = await fetch(url);
    if (!res.ok) throw new Error('fastaIdx HTTP ' + res.status + ' for ' + url);
    const body = await res.json();
    regionCache.set(key, String(body.seq).toUpperCase());
  }
  return regionCache.get(key);
}

// Blocks fetched on the primer's genomic strand, concatenated 5'->3', must equal the primer.
async function verifyPrimer(systemName, primer) {
  const g = primer.genomic;
  should.exist(g);
  g.blocks.length.should.be.aboveOrEqual(1);
  const ordered = g.strand === 1 ? g.blocks : g.blocks.slice().reverse();
  let seq = '';
  for (const b of ordered) seq += await fastaIdx(systemName, g.region, b.start, b.end, g.strand);
  seq.should.equal(primer.seq);
}

async function verifyAllPrimers(systemName, r) {
  r.pairs.length.should.be.above(0);
  for (const p of r.pairs) {
    await verifyPrimer(systemName, p.left);
    await verifyPrimer(systemName, p.right);
  }
}

function checkJunctionPairs(r) {
  r.pairs.forEach(function (p) {
    (p.left.junction !== null || p.right.junction !== null).should.be.true();
    ['left', 'right'].forEach(function (side) {
      const o = p[side];
      if (o.junction) {
        coords.spansJunction(side, o.start, o.end, o.junction.position, 7, 4).should.be.true();
        o.genomic.blocks.length.should.equal(2);
      }
    });
    p.product.genomic_size.should.equal(p.product.genomic.end - p.product.genomic.start + 1);
  });
}

after(function () {
  if (!REALDATA) return;
  console.log('design latencies (ms, in-process): ' + JSON.stringify(latencies));
  require('gramene-mongodb-config').closeMongoDatabase();
});

test('real data: transcript SORBI_3001G000200 (qPCR): 1982 nt, junctions, junction-spanning pairs of 70-150 bp, verified blocks', { skip: SKIP }, async function () {
  const r = await timed('transcript ' + G200, { mode: 'transcript', gene_id: G200 });
  r.template.length.should.equal(1982);
  r.template.features.junctions.slice(0, 3).should.eql([397, 493, 597]);
  r.template.transcript_id.should.equal('SORBI_3001G000200.1');
  r.pairs.length.should.be.within(1, 5);
  r.pairs.forEach(function (p) {
    p.product_size.should.be.within(70, 150);
    p.product_tm.should.be.a.Number();
  });
  checkJunctionPairs(r);
  r.engine.primer3.should.equal('2.6.1');
  r.warnings.should.eql([]);
  // same record as the captured fixture -> same Primer3 output
  const fixture = boulder.extractPairs(stubs.primer3Tags('transcript_junction.output.txt'));
  r.pairs.map(function (p) { return [p.left.seq, p.left.start, p.right.seq, p.right.end]; })
    .should.eql(fixture.map(function (p) { return [p.left.seq, p.left.start, p.right.seq, p.right.end]; }));
  await verifyAllPrimers('sorghum_bicolor', r);
});

test('real data: gene SORBI_3001G000200 flanks 200/100: 1:11080-15099(-), 4020 bp, exon1 201-597, CDS 499, verified blocks', { skip: SKIP }, async function () {
  const r = await timed('gene ' + G200 + ' flanks 200/100', { mode: 'gene', gene_id: G200, flank_up: 200, flank_down: 100 });
  const t = r.template;
  [t.region, t.start, t.end, t.strand, t.length].should.eql(['1', 11080, 15099, -1, 4020]);
  [t.features.exons[0].start, t.features.exons[0].end].should.eql([201, 597]);
  t.features.cds.start.should.equal(499);
  r.pairs.forEach(function (p) { p.product.should.not.have.property('genomic_size'); });
  await verifyAllPrimers('sorghum_bicolor', r);
});

test('real data: + strand SORBI_3001G000700 in gene and transcript modes: every primer\'s blocks verify through fastaIdx', { skip: SKIP }, async function () {
  const g = await timed('gene ' + G700, { mode: 'gene', gene_id: G700 });
  g.template.strand.should.equal(1);
  [g.template.start, g.template.end].should.eql([53781, 63305]);
  await verifyAllPrimers('sorghum_bicolor', g);
  const t = await timed('transcript ' + G700, { mode: 'transcript', gene_id: G700 });
  t.template.transcript_id.should.equal('SORBI_3001G000700.3');
  t.template.features.junctions.should.eql([93, 238, 299, 353, 433, 639]);
  checkJunctionPairs(t);
  await verifyAllPrimers('sorghum_bicolor', t);
});

test('real data: region 1:11080-15099(-) with target, included, excluded and params: every pair obeys every constraint', { skip: SKIP }, async function () {
  const params = { min_size: 19, max_size: 22, min_tm: 58, max_tm: 61, min_gc: 40, max_gc: 60, max_tm_diff: 2, product_size_ranges: [[200, 300]] };
  const r = await timed('region 1:11080-15099(-) constrained', {
    mode: 'region', system_name: 'sorghum_bicolor', region: { region: '1', start: 11080, end: 15099, strand: -1 },
    target: [499, 50], included: [100, 3000], excluded: [[1000, 40]], params: params
  });
  r.template.length.should.equal(4020);
  r.settings.params.should.have.properties(params);
  r.settings.preset.should.equal('pcr');
  r.pairs.length.should.be.above(0);
  const eps = 1e-3;
  r.pairs.forEach(function (p) {
    // covers the target: the primers flank [499, 548]
    p.left.end.should.be.below(499);
    p.right.start.should.be.above(548);
    (p.product.start <= 499 && p.product.end >= 548).should.be.true();
    // inside included [100, 3099], outside excluded [1000, 1039]
    p.left.start.should.be.aboveOrEqual(100);
    p.right.end.should.be.belowOrEqual(3099);
    coords.overlapsRuns(p.left.start, p.left.end, [[1000, 40]]).should.be.false();
    coords.overlapsRuns(p.right.start, p.right.end, [[1000, 40]]).should.be.false();
    p.product_size.should.be.within(200, 300);
    ['left', 'right'].forEach(function (side) {
      const o = p[side];
      o.len.should.be.within(19, 22);
      o.tm.should.be.within(58 - eps, 61 + eps);
      o.gc.should.be.within(40 - eps, 60 + eps);
    });
    Math.abs(p.left.tm - p.right.tm).should.be.belowOrEqual(2 + eps);
  });
  await verifyAllPrimers('sorghum_bicolor', r);
});

test('real data: sequence mode with IUPAC codes -> 200 with IUPAC_CONVERTED, not PRIMER3_INPUT_ERROR', { skip: SKIP }, async function () {
  const cdna = stubs.cdna200();
  const seq = cdna.slice(0, 49) + 'R' + cdna.slice(50, 299) + 'Y' + cdna.slice(300, 449) + 'K' + cdna.slice(450, 600);
  const r = await timed('sequence 600 nt with IUPAC', { mode: 'sequence', sequence: '>iupac_test\n' + seq, system_name: 'sorghum_bicolor' });
  r.warnings.map(function (w) { return w.code; }).should.containEql('IUPAC_CONVERTED');
  r.template.seq.charAt(49).should.equal('N');
  r.template.seq.charAt(299).should.equal('N');
  r.template.length.should.equal(600);
  r.pairs.length.should.be.above(0);
  r.pairs.forEach(function (p) {
    should(p.left.genomic).be.null();
    should(p.product.genomic).be.null();
    p.left.seq.should.not.match(/N/);
  });
});

test('real data: avoid_repeats on sorghum_tx436pac 4:7547610-7564601 uses the real soft-mask; no primer overlaps it', { skip: SKIP }, async function () {
  const assemblies = require('../../../api/helpers/primers/assemblies');
  const sequence = require('../../../api/helpers/primers/sequence');
  const asm = await assemblies.resolve('sorghum_tx436pac');
  asm.repeat_masking.should.equal('soft_masked');
  const r = await timed('region sorghum_tx436pac avoid_repeats (softmask)', {
    mode: 'region', system_name: 'sorghum_tx436pac', region: { region: '4', start: 7547610, end: 7564601 }, avoid_repeats: true
  });
  r.template.mask_source.should.equal('softmask');
  r.template.masked.should.be.true();
  r.template.masked_fraction.should.be.above(0);
  const sm = await sequence.fetch(asm.fasta.dna_sm, '4', 7547610, 7564601, 1);
  r.template.mask.should.eql(coords.lowercaseRuns(sm));
  coords.maskedBases(r.template.mask).should.equal(3238);
  r.warnings.map(function (w) { return w.code; }).should.not.containEql('BLAST_DEPTH_MASK');
  r.settings.params.max_ns.should.equal(0);
  r.pairs.length.should.be.above(0);
  for (const p of r.pairs) {
    for (const side of ['left', 'right']) {
      const o = p[side];
      coords.overlapsRuns(o.start, o.end, r.template.mask).should.be.false();
      const g = o.genomic;
      const ordered = g.strand === 1 ? g.blocks : g.blocks.slice().reverse();
      let s = '';
      for (const b of ordered) s += (await sequence.fetch(asm.fasta.dna, g.region, b.start, b.end, g.strand)).toUpperCase();
      s.should.equal(o.seq);
    }
  }
});

test('real data: n_mask product_tm on sorghum_tx436pac 4:7547610-7564601 target [2483,150] equals primer3_core on the unmasked template (<= 0.01 C)', { skip: SKIP }, async function () {
  const primer3 = require('../../../api/helpers/primers/primer3');
  const body = { mode: 'region', system_name: 'sorghum_tx436pac', region: { region: '4', start: 7547610, end: 7564601 }, avoid_repeats: true,
    target: [2483, 150], params: { product_size_ranges: [[300, 1000]] } };
  const r = await timed('region sorghum_tx436pac n_mask target [2483,150]', body);
  r.template.mask_source.should.equal('softmask');
  r.settings.repeat_mask_mode.should.equal('n_mask');
  r.pairs.length.should.be.above(0);
  const rows = [];
  for (const p of r.pairs) {
    const forced = await primer3.run(stubs.forcedPairRecord(r.template.seq, p, {}), { timeoutMs: 30000 });
    should(forced.error).be.null();
    forced.tags.PRIMER_LEFT_0.should.equal(p.left.start + ',' + p.left.len);
    const own = Number(forced.tags.PRIMER_PAIR_0_PRODUCT_TM);
    rows.push({ rank: p.rank, product: p.product.start + '-' + p.product.end, size: p.product_size,
      masked_bases: stubs.maskedBasesIn(r.template.mask, p.product.start, p.product.end), product_tm: p.product_tm, primer3_unmasked: own });
    Math.abs(p.product_tm - own).should.be.belowOrEqual(0.01);
  }
  rows.filter(function (x) { return x.masked_bases > 0; }).length.should.be.above(0);
  // three_prime (Primer3 upcases the lowercase mask itself) reports the same product Tm for any pair both modes return
  const three = await timed('region sorghum_tx436pac three_prime target [2483,150]', Object.assign({}, body, { repeat_mask_mode: 'three_prime' }));
  let shared = 0;
  three.pairs.forEach(function (q) {
    const p = r.pairs.find(function (x) { return x.left.seq === q.left.seq && x.right.seq === q.right.seq && x.product.start === q.product.start; });
    if (!p) return;
    shared++;
    Math.abs(p.product_tm - q.product_tm).should.be.belowOrEqual(0.01);
  });
  console.log('n_mask product_tm vs primer3_core on the unmasked template: ' + JSON.stringify(rows) + '; pairs shared with three_prime: ' + shared);
});

test('real data: SORBI_3004G087700 avoid_repeats + template_only -> blast_depth mask, BLAST_DEPTH_MASK, pairs []', { skip: SKIP }, async function () {
  const body = { mode: 'gene', gene_id: G087700, avoid_repeats: true, template_only: true };
  const r = await timed('gene ' + G087700 + ' avoid_repeats template_only (megablast)', body);
  r.template.mask_source.should.equal('blast_depth');
  r.warnings.map(function (w) { return w.code; }).should.containEql('BLAST_DEPTH_MASK');
  r.pairs.should.eql([]);
  should(r.explain).be.null();
  r.template.masked.should.be.true();
  r.template.masked_fraction.should.be.above(0);
  r.template.mask.forEach(function (run) {
    run[0].should.be.aboveOrEqual(1);
    (run[0] + run[1] - 1).should.be.belowOrEqual(r.template.length);
  });
  const again = await timed('gene ' + G087700 + ' avoid_repeats template_only (cached mask)', body);
  again.template.mask.should.eql(r.template.mask);
});

test('real data: SORBI_3004G087700 transcript + avoid_repeats projects the exon megablast mask onto the cDNA', { skip: SKIP }, async function () {
  const r = await timed('transcript ' + G087700 + ' avoid_repeats (megablast exons +-100)', { mode: 'transcript', gene_id: G087700, avoid_repeats: true });
  r.template.mask_source.should.equal('blast_depth');
  r.template.length.should.equal(4878);
  r.template.mask.forEach(function (run) {
    (run[0] + run[1] - 1).should.be.belowOrEqual(4878);
  });
  r.pairs.forEach(function (p) {
    coords.overlapsRuns(p.left.start, p.left.end, r.template.mask).should.be.false();
    coords.overlapsRuns(p.right.start, p.right.end, r.template.mask).should.be.false();
  });
  if (r.pairs.length > 0) await verifyAllPrimers('sorghum_bicolor', r);
});

test('real data: error cases (UNKNOWN_GENE, INVALID_PARAMS, UNKNOWN_GENOME, REGION_OUT_OF_BOUNDS, UNKNOWN_TRANSCRIPT)', { skip: SKIP }, async function () {
  await rejectsCode(design.design({ mode: 'gene', gene_id: 'NOPE' }), 'UNKNOWN_GENE', 404);
  await rejectsCode(design.design({ mode: 'transcript', gene_id: G200, params: { max_size: 24, min_5_prime_overlap_of_junction: 13 } }), 'INVALID_PARAMS', 400);
  await rejectsCode(design.design({ mode: 'region', system_name: 'no_such_genome', region: { region: '1', start: 1, end: 100 } }), 'UNKNOWN_GENOME', 404);
  await rejectsCode(design.design({ mode: 'region', system_name: 'sorghum_bicolor', region: { region: '1', start: 80884300, end: 80884400 } }),
    'REGION_OUT_OF_BOUNDS', 400);
  await rejectsCode(design.design({ mode: 'transcript', gene_id: G200, transcript_id: 'SORBI_3001G000200.9' }), 'UNKNOWN_TRANSCRIPT', 404);
});
