'use strict';

// check/cost.js against the spec §B.13 table, with the pan-genome CPU factor measured on the full sorghum panel.

const test = require('node:test');
const should = require('should');

const cost = require('../../../api/helpers/primers/check/cost');
const config = require('../../../api/helpers/primers/config');
const genomes = require('../../../api/helpers/primers/genomes');
const mapsFixture = require('../fixtures/catalog/maps_sorghum_v11.json');
const taxonomyFixture = require('../fixtures/catalog/taxonomy_sorghum_v11.json');

function makeCfg(overrides) {
  return config._build({ env: {}, fileConfig: {}, overrides: overrides }).config;
}

const CFG = makeCfg();
// sorghum_bicolor dna .fai total (verified by the catalog role).
const BICOLOR = { system_name: 'sorghum_bicolor', total_bases: 708735318 };

// The 119 other sorghum_v11 sorghum assemblies, sized by the sum of their map region lengths (~ total_bases).
function otherSorghum() {
  const cat = genomes.buildCatalog(mapsFixture, taxonomyFixture);
  const ref = cat.bySystemName.get('sorghum_bicolor');
  return cat.bySpecies.get(genomes.speciesKey(ref))
    .filter(function (g) { return g.system_name !== 'sorghum_bicolor'; })
    .map(function (g) {
      return { system_name: g.system_name, total_bases: g.regions.lengths.reduce(function (a, b) { return a + b; }, 0) };
    });
}

function within(value, target, fraction) {
  value.should.be.within(target * (1 - fraction), target * (1 + fraction));
}

test('coefficients come from cfg.check.cpu_s_per_primer_gb (5.2 / 2.2 / 1.2)', function () {
  cost.coefficient(5, CFG).should.equal(5.2);
  cost.coefficient(6, CFG).should.equal(2.2);
  cost.coefficient(7, CFG).should.equal(1.2);
  // Unknown word sizes: the largest configured ws <= requested, else the most expensive coefficient.
  cost.coefficient(9, CFG).should.equal(1.2);
  cost.coefficient(4, CFG).should.equal(5.2);
  cost.coefficient(5, makeCfg({ check: { cpu_s_per_primer_gb: { ws5: null, ws6: 3 } } })).should.equal(3);
  cost.CDNA_GB_ESTIMATE.should.equal(0.15);
});

test('B.13: specificity only, 10 primers, ws5 on sorghum_bicolor ~ 37 CPU-s BLAST + 6 CPU-s re-alignment', function () {
  const est = cost.estimate({ unique_primers: 10, mode: 'gene', reference: BICOLOR, pangenome: [], cfg: CFG });
  est.cpu_s.should.equal(43); // ceil(36.85 + 10 × 1 task × 0.6)
  est.total.should.equal(1);
  est.breakdown.should.eql({ reference: 36.9, transcriptome: 0, pangenome: 0, realign: 6 });
  est.word_sizes.should.eql({ reference: 5, pangenome: 6 });
  // region and sequence modes cost the same as gene mode
  cost.estimate({ unique_primers: 10, mode: 'sequence', reference: BICOLOR, cfg: CFG }).cpu_s.should.equal(43);
});

test('re-alignment term: primers × genome tasks × 0.6 CPU-s (configurable); cDNA tasks are not charged', function () {
  cost.DEFAULT_REALIGN_CPU_S_PER_PRIMER_TASK.should.equal(0.6);
  cost.realignCoefficient(CFG).should.equal(0.6);
  const pan3 = [{ total_bases: 7e8 }, { total_bases: 7e8 }, { total_bases: 7e8 }];
  cost.estimate({ unique_primers: 4, mode: 'gene', reference: BICOLOR, pangenome: pan3, cfg: CFG }).breakdown.realign.should.equal(9.6); // 4 × 4 × 0.6
  cost.estimate({ unique_primers: 4, mode: 'transcript', reference: BICOLOR, pangenome: pan3, cfg: CFG }).breakdown.realign.should.equal(2.4); // reference genome only
  const zero = makeCfg({ check: { realign_cpu_s_per_primer_task: 0 } });
  const est = cost.estimate({ unique_primers: 10, mode: 'gene', reference: BICOLOR, cfg: zero });
  est.breakdown.realign.should.equal(0);
  est.cpu_s.should.equal(37);
  cost.estimate({ unique_primers: 10, mode: 'gene', reference: BICOLOR, cfg: makeCfg({ check: { realign_cpu_s_per_primer_task: 1.5 } }) }).breakdown.realign.should.equal(15);
});

test('pan-genome CPU factor: default 2.0 multiplies only the pan-genome BLAST term; configurable, values below 1 fall back', function () {
  cost.DEFAULT_PANGENOME_CPU_FACTOR.should.equal(2);
  cost.pangenomeFactor(CFG).should.equal(2);
  cost.pangenomeFactor(makeCfg({ check: { pangenome_cpu_factor: 1.5 } })).should.equal(1.5);
  cost.pangenomeFactor(makeCfg({ check: { pangenome_cpu_factor: 1 } })).should.equal(1);
  cost.pangenomeFactor(makeCfg({ check: { pangenome_cpu_factor: 0.5 } })).should.equal(2);
  cost.pangenomeFactor(makeCfg({ check: { pangenome_cpu_factor: 'x' } })).should.equal(2);
  const pan3 = [{ total_bases: 7e8 }, { total_bases: 7e8 }, { total_bases: 7e8 }];
  const two = cost.estimate({ unique_primers: 4, mode: 'gene', reference: BICOLOR, pangenome: pan3, cfg: CFG });
  const one = cost.estimate({ unique_primers: 4, mode: 'gene', reference: BICOLOR, pangenome: pan3, cfg: makeCfg({ check: { pangenome_cpu_factor: 1 } }) });
  two.breakdown.pangenome.should.equal(37); // 4 × 2.1 Gb × 2.2 × 2 = 36.96
  one.breakdown.pangenome.should.equal(18.5); // 4 × 2.1 Gb × 2.2 = 18.48
  two.breakdown.reference.should.equal(one.breakdown.reference);
  two.breakdown.realign.should.equal(one.breakdown.realign);
  // The specificity-only estimate has no pan-genome term, so the factor does not change it.
  cost.estimate({ unique_primers: 10, mode: 'gene', reference: BICOLOR, cfg: makeCfg({ check: { pangenome_cpu_factor: 5 } }) }).cpu_s.should.equal(43);
});

test('measured full panel: P3 (2 primers) against all 119 other sorghum genomes estimates above the 758.2 CPU-s it used', function () {
  const others = otherSorghum();
  const est = cost.estimate({ unique_primers: 2, mode: 'gene', reference: BICOLOR, pangenome: others, cfg: CFG });
  est.cpu_s.should.be.above(758);
  est.cpu_s.should.be.below(1000);
  const unfactored = cost.estimate({ unique_primers: 2, mode: 'gene', reference: BICOLOR, pangenome: others, cfg: makeCfg({ check: { pangenome_cpu_factor: 1 } }) });
  unfactored.cpu_s.should.be.below(758);
});

test('B.13 with the factor: 10 primers against all 119 other sorghum assemblies ~ 3,700 CPU-s BLAST (+ 720 re-alignment), under the cap', function () {
  const others = otherSorghum();
  others.length.should.equal(119);
  const est = cost.estimate({ unique_primers: 10, mode: 'gene', reference: BICOLOR, pangenome: others, cfg: CFG });
  within(est.breakdown.pangenome, 3600, 0.05); // the table's ~1,800 single-thread CPU-s × 2.0
  est.breakdown.realign.should.equal(720); // 10 × 120 × 0.6
  est.total.should.equal(120);
  est.cpu_s.should.be.within(Math.floor(est.breakdown.reference + est.breakdown.pangenome + est.breakdown.realign),
    Math.ceil(est.breakdown.reference + est.breakdown.pangenome + est.breakdown.realign) + 1);
  est.cpu_s.should.be.below(6000);
  cost.assertWithinLimit(est, CFG).should.equal(est);
});

test('B.13 with the factor: 20 primers against all 119 other sorghum assemblies exceed the 6,000 CPU-s cap', function () {
  const est = cost.estimate({ unique_primers: 20, mode: 'region', reference: BICOLOR, pangenome: otherSorghum(), cfg: CFG });
  within(est.breakdown.pangenome, 7200, 0.05);
  est.breakdown.realign.should.equal(1440);
  est.cpu_s.should.be.above(6000);
  should.throws(function () { cost.assertWithinLimit(est, CFG); }, /CPU-seconds/);
});

test('B.13: transcript-mode pan-genome (cDNA ~0.1 Gb each), 10 primers ~ 260 CPU-s single-thread, doubled by the factor; the API charges 0.15 Gb', function () {
  const others = otherSorghum();
  const table = cost.estimate({ unique_primers: 10, mode: 'transcript', reference: BICOLOR, pangenome: others, cfg: CFG, cdna_gb: 0.1 });
  table.breakdown.pangenome.should.equal(523.6); // 10 × 119 × 0.1 × 2.2 × 2
  table.breakdown.realign.should.equal(6);
  table.total.should.equal(121);
  const api = cost.estimate({ unique_primers: 10, mode: 'transcript', reference: BICOLOR, pangenome: others, cfg: CFG });
  api.breakdown.pangenome.should.equal(785.4); // 10 × 119 × 0.15 × 2.2 × 2
});

test('cost guard: over max_job_cpu_s -> 422 JOB_TOO_LARGE {estimate_cpu_s, limit}; equal to the limit passes', function () {
  const synthetic = [];
  for (let i = 0; i < 150; i++) synthetic.push({ system_name: 'big_' + i, total_bases: 2e9 });
  const big = cost.estimate({ unique_primers: 20, mode: 'gene', reference: BICOLOR, pangenome: synthetic, cfg: CFG });
  big.cpu_s.should.be.above(6000);
  let err;
  try {
    cost.assertWithinLimit(big, CFG);
  } catch (e) {
    err = e;
  }
  should.exist(err);
  err.status.should.equal(422);
  err.code.should.equal('JOB_TOO_LARGE');
  err.details.should.eql({ estimate_cpu_s: big.cpu_s, limit: 6000 });

  const spec = cost.estimate({ unique_primers: 10, mode: 'gene', reference: BICOLOR, cfg: CFG });
  should.doesNotThrow(function () { cost.assertWithinLimit(spec, makeCfg({ check: { max_job_cpu_s: 43 } })); });
  should.throws(function () { cost.assertWithinLimit(spec, makeCfg({ check: { max_job_cpu_s: 42 } })); }, /CPU-seconds/);
  cost.limitOf(CFG).should.equal(6000);
});

test('word sizes come from config; an assembly without total_bases is charged 1 Gb', function () {
  const est = cost.estimate({
    unique_primers: 1,
    mode: 'gene',
    reference: { system_name: 'x', total_bases: 1e9 },
    pangenome: [],
    cfg: makeCfg({ check: { reference_word_size: 7 } })
  });
  est.breakdown.reference.should.equal(1.2);
  cost.genomeGb({ total_bases: null }).should.equal(cost.FALLBACK_GENOME_GB);
  cost.estimate({ unique_primers: 2, mode: 'gene', reference: { total_bases: null }, cfg: CFG }).cpu_s.should.equal(12); // 2 x 1 Gb x 5.2 + 2 x 0.6
  cost.estimate({ unique_primers: 0, mode: 'gene', reference: BICOLOR, cfg: CFG }).cpu_s.should.equal(0);
});

// Compatibility C4 (genotyping spec §7.2, §5.4): every estimate above that is not already pinned in full (the
// specificity-only test pins its own) is pinned here, so the genotyping term cannot move an existing number; without
// genotyping, breakdown keeps exactly its four keys (breakdown.genotyping appears only when genotyping is requested).
test('C4 compat: without genotyping every estimate above is unchanged and breakdown has exactly its four keys', function () {
  const others = otherSorghum();
  const pan3 = [{ total_bases: 7e8 }, { total_bases: 7e8 }, { total_bases: 7e8 }];
  const synthetic = [];
  for (let i = 0; i < 150; i++) synthetic.push({ system_name: 'big_' + i, total_bases: 2e9 });
  const factor1 = makeCfg({ check: { pangenome_cpu_factor: 1 } });
  // [case, estimate() options, cpu_s, total, breakdown [reference, transcriptome, pangenome, realign], word_sizes]
  const rows = [
    ['sequence mode, 10 primers', { unique_primers: 10, mode: 'sequence', reference: BICOLOR, cfg: CFG }, 43, 1, [36.9, 0, 0, 6]],
    ['gene, 4 primers, 3 pan-genome genomes', { unique_primers: 4, mode: 'gene', reference: BICOLOR, pangenome: pan3, cfg: CFG }, 62, 4, [14.7, 0, 37, 9.6]],
    ['the same with genotyping: false', { unique_primers: 4, mode: 'gene', reference: BICOLOR, pangenome: pan3, cfg: CFG, genotyping: false }, 62, 4, [14.7, 0, 37, 9.6]],
    ['transcript, 4 primers, 3 pan-genome genomes', { unique_primers: 4, mode: 'transcript', reference: BICOLOR, pangenome: pan3, cfg: CFG }, 29, 5, [14.7, 3.1, 7.9, 2.4]],
    ['re-alignment coefficient 0', { unique_primers: 10, mode: 'gene', reference: BICOLOR, cfg: makeCfg({ check: { realign_cpu_s_per_primer_task: 0 } }) }, 37, 1, [36.9, 0, 0, 0]],
    ['re-alignment coefficient 1.5', { unique_primers: 10, mode: 'gene', reference: BICOLOR, cfg: makeCfg({ check: { realign_cpu_s_per_primer_task: 1.5 } }) }, 52, 1, [36.9, 0, 0, 15]],
    ['pan-genome factor 1, 3 genomes', { unique_primers: 4, mode: 'gene', reference: BICOLOR, pangenome: pan3, cfg: factor1 }, 43, 4, [14.7, 0, 18.5, 9.6]],
    ['pan-genome factor 5, specificity only', { unique_primers: 10, mode: 'gene', reference: BICOLOR, cfg: makeCfg({ check: { pangenome_cpu_factor: 5 } }) }, 43, 1, [36.9, 0, 0, 6]],
    ['P3 (2 primers), full panel', { unique_primers: 2, mode: 'gene', reference: BICOLOR, pangenome: others, cfg: CFG }, 889, 120, [7.4, 0, 737.5, 144]],
    ['P3 (2 primers), full panel, factor 1', { unique_primers: 2, mode: 'gene', reference: BICOLOR, pangenome: others, cfg: factor1 }, 521, 120, [7.4, 0, 368.7, 144]],
    ['10 primers, full panel', { unique_primers: 10, mode: 'gene', reference: BICOLOR, pangenome: others, cfg: CFG }, 4445, 120, [36.9, 0, 3687.3, 720]],
    ['20 primers, full panel, region mode', { unique_primers: 20, mode: 'region', reference: BICOLOR, pangenome: others, cfg: CFG }, 8889, 120, [73.7, 0, 7374.7, 1440]],
    ['transcript, 10 primers, full panel, cDNA 0.1 Gb', { unique_primers: 10, mode: 'transcript', reference: BICOLOR, pangenome: others, cfg: CFG, cdna_gb: 0.1 }, 572, 121, [36.9, 5.2, 523.6, 6]],
    ['transcript, 10 primers, full panel, API cDNA 0.15 Gb', { unique_primers: 10, mode: 'transcript', reference: BICOLOR, pangenome: others, cfg: CFG }, 837, 121, [36.9, 7.8, 785.4, 6]],
    ['20 primers, 150 synthetic 2 Gb genomes', { unique_primers: 20, mode: 'gene', reference: BICOLOR, pangenome: synthetic, cfg: CFG }, 28286, 151, [73.7, 0, 26400, 1812]],
    ['reference word size 7, 1 Gb', { unique_primers: 1, mode: 'gene', reference: { system_name: 'x', total_bases: 1e9 }, pangenome: [], cfg: makeCfg({ check: { reference_word_size: 7 } }) }, 2, 1, [1.2, 0, 0, 0.6], { reference: 7, pangenome: 6 }],
    ['reference without total_bases', { unique_primers: 2, mode: 'gene', reference: { total_bases: null }, cfg: CFG }, 12, 1, [10.4, 0, 0, 1.2]],
    ['no primers', { unique_primers: 0, mode: 'gene', reference: BICOLOR, cfg: CFG }, 0, 1, [0, 0, 0, 0]]
  ];
  rows.forEach(function (row) {
    const est = cost.estimate(row[1]);
    Object.keys(est.breakdown).sort().should.eql(['pangenome', 'realign', 'reference', 'transcriptome'], row[0]);
    est.should.eql({
      cpu_s: row[2],
      total: row[3],
      breakdown: { reference: row[4][0], transcriptome: row[4][1], pangenome: row[4][2], realign: row[4][3] },
      word_sizes: row[5] || { reference: 5, pangenome: 6 }
    }, row[0]);
  });
});

// Genotyping spec §5.4 (M7): a genotyping check adds genome_tasks × check.genotype_cpu_s_per_genome (0.2) and
// breakdown.genotyping, only when requested. The spec's table uses the real genome sizes, sorghum_bicolor 708,735,318 bases and
// 119 other assemblies summing to 83,774,241,795; only that sum enters the estimate, so the panel is 119 parts of it.
function realPanel() {
  const total = 83774241795;
  const each = Math.floor(total / 119);
  return Array.from({ length: 119 }, function (_, i) {
    return { system_name: 'genome_' + i, total_bases: each + (i === 0 ? total - each * 119 : 0) };
  });
}

test('genotyping term (§5.4): 120 genome tasks × 0.2 CPU-s and breakdown.genotyping; the spec table over the full panel', function () {
  cost.DEFAULT_GENOTYPE_CPU_S_PER_GENOME.should.equal(0.2);
  cost.genotypeCoefficient(CFG).should.equal(0.2);
  const panel = realPanel();
  // [unique primers, CPU-s without genotyping, with genotyping]
  [[3, 1333, 1357], [6, 2666, 2690], [12, 5332, 5356], [13, 5776, 5800], [14, 6221, 6245]].forEach(function (row) {
    const opts = { unique_primers: row[0], mode: 'region', reference: BICOLOR, pangenome: panel, cfg: CFG };
    const without = cost.estimate(opts);
    const withG = cost.estimate(Object.assign({ genotyping: true }, opts));
    without.cpu_s.should.equal(row[1], row[0] + ' primers');
    withG.cpu_s.should.equal(row[2], row[0] + ' primers with genotyping');
    Object.keys(without.breakdown).should.eql(['reference', 'transcriptome', 'pangenome', 'realign']);
    Object.keys(withG.breakdown).should.eql(['reference', 'transcriptome', 'pangenome', 'realign', 'genotyping']);
    withG.breakdown.should.eql(Object.assign({ genotyping: 24 }, without.breakdown));
    withG.total.should.equal(without.total);
    withG.word_sizes.should.eql(without.word_sizes);
  });
  // 13 primers is the largest full-panel genotyping check (check_max_unique_primers, §4.18); 14 primers are refused
  const opts = function (n) { return { unique_primers: n, mode: 'region', reference: BICOLOR, pangenome: panel, cfg: CFG, genotyping: true }; };
  cost.assertWithinLimit(cost.estimate(opts(13)), CFG).cpu_s.should.equal(5800);
  should.throws(function () { cost.assertWithinLimit(cost.estimate(opts(14)), CFG); }, /estimated at 6245 CPU-seconds/);
});

test('genotyping term: the §2.11 check (6 primers, 3 genomes) stays at 95 CPU-s; the coefficient comes from config', function () {
  const cat = genomes.buildCatalog(mapsFixture, taxonomyFixture);
  const three = ['sorghum_bicolorv5', 'sorghum_pi180348', 'sorghum_pi329250'].map(function (n) {
    return { system_name: n, total_bases: cat.bySystemName.get(n).regions.lengths.reduce(function (a, b) { return a + b; }, 0) };
  });
  const opts = { unique_primers: 6, mode: 'region', reference: BICOLOR, pangenome: three, cfg: CFG };
  cost.estimate(opts).cpu_s.should.equal(95); // 94.16
  const withG = cost.estimate(Object.assign({ genotyping: true }, opts));
  withG.cpu_s.should.equal(95); // 94.96
  withG.breakdown.genotyping.should.equal(0.8);
  cost.estimate({ unique_primers: 3, mode: 'gene', reference: BICOLOR, cfg: CFG, genotyping: true }).breakdown.genotyping.should.equal(0.2);
  cost.estimate({ unique_primers: 0, mode: 'region', reference: BICOLOR, pangenome: three, cfg: CFG, genotyping: true }).cpu_s.should.equal(1);

  const one = makeCfg({ check: { genotype_cpu_s_per_genome: 1 } });
  cost.genotypeCoefficient(one).should.equal(1);
  cost.estimate(Object.assign({}, opts, { cfg: one, genotyping: true })).should.match({ cpu_s: 99, breakdown: { genotyping: 4 } });
  cost.genotypeCoefficient(makeCfg({ check: { genotype_cpu_s_per_genome: 0 } })).should.equal(0);
  cost.genotypeCoefficient(makeCfg({ check: { genotype_cpu_s_per_genome: -1 } })).should.equal(0.2);
  cost.genotypeCoefficient(makeCfg({ check: { genotype_cpu_s_per_genome: 'x' } })).should.equal(0.2);
  // a non-genotyping estimate ignores the coefficient entirely
  cost.estimate(Object.assign({}, opts, { cfg: one })).should.eql(cost.estimate(opts));
});
