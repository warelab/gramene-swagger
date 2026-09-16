'use strict';

// Real-data §10.4 checks through check/run.js (opt-in):
//   PRIMERS_REALDATA=1 node --test --test-concurrency=1 --test-reporter=spec test/primers/unit/run_realdata.test.js
// Read-only: mongo (genes, maps, taxonomy), /scratch/olson/fasta bgzip FASTA and BLAST DBs.
// blastn runs at nice 10 with ctx.procs = 4 (reference: 4 threads; pan-genome: <= 4 single-thread
// processes), one test at a time. Stringency: a product amplifies only with <= max_amplifying_mismatches (3)
// edits per primer and no 3'-blocked primer; 4-5 mismatch products are 'unlikely'.

const { describe, it, before, after } = require('node:test');
const fs = require('fs');
const should = require('should');

const REAL = process.env.PRIMERS_REALDATA === '1';
const PROCS = 4;
const CAP = 3;

describe('real data: §10.4 checks on the sorghum_bicolor reference', { skip: REAL ? false : 'set PRIMERS_REALDATA=1' }, () => {
  let H;
  let P;
  let mongo;
  let cfg;
  let runWithStats;
  const timings = {};

  before(() => {
    H = require('../fixtures/verdicts/realdata_ctx');
    P = H.P;
    mongo = H.mongoConfig();
    cfg = H.primersConfig();
    runWithStats = require('../../../api/helpers/primers/check/run').runWithStats;
  });

  after(() => {
    if (mongo) mongo.closeMongoDatabase();
  });

  async function check(name, body, t) {
    const request = H.normalizedRequest(body);
    const t0 = Date.now();
    const resolved = await H.buildResolved({ cfg, mongo, systemName: request.system_name, genomes: request.genomes, geneId: request.gene_id });
    const resolveMs = Date.now() - t0;
    const warnings = [];
    const { ctx, events, cleanup } = H.makeCtx({ cfg, resolved, mongo, procs: PROCS, jobId: name, log: { info() {}, warn: (m) => warnings.push(m), error: (m) => warnings.push(m) } });
    const calls = [];
    const spawn = ctx.spawnLines;
    ctx.spawnLines = (cmd, args, o) => {
      calls.push({ cmd, args: args.slice(), cwd: o && o.cwd });
      return spawn(cmd, args, o);
    };
    try {
      const t1 = Date.now();
      const { results, stats } = await runWithStats(request, ctx);
      timings[name] = { resolve_ms: resolveMs, run_ms: Date.now() - t1, timings_ms: results.timings_ms, stats };
      t.diagnostic(name + ' ' + JSON.stringify(timings[name]));
      if (warnings.length) t.diagnostic(name + ' log ' + JSON.stringify(warnings));
      // every run: BLAST+ children in the private job directory, summary invariants, a partial after specificity
      const dir = ctx.tmpdir();
      should(calls.length).be.above(0);
      should(calls.every((c) => c.cwd === dir)).be.true();
      should(fs.statSync(dir).mode & 0o777).equal(0o700);
      should(events.partials.length).be.aboveOrEqual(1);
      should(events.partials[0].specificity).be.true();
      should(results.engine).match({ algorithm_version: '2', blast: '2.13.0' });
      should(results.params.max_amplifying_mismatches).equal(CAP);
      if (results.pangenome) {
        for (const pair of results.pangenome.pairs) {
          const s = pair.summary;
          should(s.single_perfect + s.single_mismatch + s.multiple + s.no_amplicon + s.db_unavailable + s.error).equal(s.genomes_total);
          should(s.amplifies).equal(s.single_perfect + s.single_mismatch + s.multiple);
          should(s.genomes_total).equal(request.genomes.length);
          should(s.truncated).equal(pair.genomes.filter((g) => g.truncated).length);
        }
      }
      return { request, results, events, calls, stats };
    } finally {
      cleanup();
    }
  }

  const find = (list, region, start, end, orientation) => (list || []).find((a) => a.region === region && a.start === start && a.end === end && (!orientation || a.orientation === orientation));
  const geneIds = (a) => (a && Array.isArray(a.genes) ? a.genes.map((g) => g.id) : []);
  const brief = (a) => (a ? [a.region + ':' + a.start + '-' + a.end, a.orientation, a.size, a.likelihood, a.left_mm, a.right_mm, a.left_mm_pos, a.right_mm_pos, a.approx] : null);
  const amplifying = (a) => a.likelihood === 'likely' || a.likelihood === 'likely_weak';

  // Every listed amplifying product has each primer within the cap; every unlikely one is over it or 3'-blocked.
  function assertStringency(block) {
    const listed = (block.on_target ? [block.on_target] : []).concat(block.off_targets || [], block.gdna_products || []);
    for (const a of listed) {
      should(amplifying(a)).be.true();
      should(a.left_mm).be.belowOrEqual(CAP);
      should(a.right_mm).be.belowOrEqual(CAP);
    }
    for (const a of block.unlikely || []) should(a.likelihood).equal('unlikely');
  }

  it('gene mode: P1, P2, P3 off-targets; P5_L repetitive; m20_7_14 and del18 sites', { timeout: 600000 }, async (t) => {
    const { results } = await check('gene_mode_6_pairs', {
      system_name: 'sorghum_bicolor',
      mode: 'gene',
      gene_id: 'SORBI_3004G087700',
      params: { include_unlikely: true },
      pairs: [
        { id: 'P1', left: P.P1_L, right: P.P1_R, expected: { region: '4', start: 7422190, end: 7422843 } },
        { id: 'P2', left: P.P2_L, right: P.P2_R, expected: { region: '4', start: 7423537, end: 7423746 } },
        { id: 'P3', left: P.P3_L, right: P.P3_R, expected: { region: '4', start: 7422482, end: 7423061 } },
        { id: 'P5_L', left: P.P5_L, right: P.P2_R },
        { id: 'm20_7_14', left: P.m20_7_14, right: P.P2_R },
        { id: 'del18', left: P.del18, right: P.P2_R }
      ]
    }, t);
    const pair = {};
    for (const p of results.specificity.pairs) pair[p.id] = p;
    t.diagnostic('verdicts [id, verdict, off_target_count, unlikely_count] ' + JSON.stringify(results.specificity.pairs.map((p) => [p.id, p.verdict, p.off_target_count, p.unlikely_count])));
    should(results.reference).eql({ system_name: 'sorghum_bicolor', map_id: 'GCA_000003195.3', total_bases: 708735318 });
    for (const p of results.specificity.pairs) assertStringency(p);

    // P2
    t.diagnostic('P2 off_targets ' + JSON.stringify(pair.P2.off_targets.map(brief)));
    should(pair.P2).match({ verdict: 'off_targets', off_target_count: 2 });
    should(pair.P2.on_target).match({ region: '4', start: 7423537, end: 7423746, orientation: 'LR', left_mm: 0, right_mm: 0, likelihood: 'likely' });
    should(geneIds(pair.P2.on_target)).containEql('SORBI_3004G087700');
    const p2para = find(pair.P2.off_targets, '4', 7437317, 7437526, 'LR');
    should(p2para).match({ left_mm: 0, right_mm: 0, likelihood: 'likely' });
    should(geneIds(p2para)).containEql('SORBI_3004G087800');
    const p2chr5 = find(pair.P2.off_targets, '5', 66890615, 66890824, 'RL');
    should(p2chr5).match({ left_mm_pos: [10, 7], right_mm_pos: [15, 10], left_3p_mm: 0, right_3p_mm: 0, likelihood: 'likely' });
    should(geneIds(p2chr5)).containEql('SORBI_3005G183900');

    // P1: off-targets 4:7435938-7436628 (691 bp) and 5:66891530-66892235 (706 bp)
    should(pair.P1.on_target).match({ region: '4', start: 7422190, end: 7422843, size: 654 });
    const p1a = find(pair.P1.off_targets.concat(pair.P1.unlikely), '4', 7435938, 7436628);
    const p1b = find(pair.P1.off_targets.concat(pair.P1.unlikely), '5', 66891530, 66892235);
    t.diagnostic('P1 §10.4 products ' + JSON.stringify([brief(p1a), brief(p1b)]));
    should(p1a).match({ size: 691 });
    should(p1b).match({ size: 706 });
    // 4:7435938 still amplifies (L 3 / R 2 edits); 5:66891530 has a 4-edit L site: now unlikely (approximate bound)
    should(p1a).match({ likelihood: 'likely', left_mm: 3, right_mm: 2, approx: false });
    should(pair.P1.off_targets).containEql(p1a);
    should(p1b).match({ likelihood: 'unlikely', left_mm: 4, right_mm: 0, approx: true });
    should(pair.P1.unlikely).containEql(p1b);
    should(pair.P1).match({ verdict: 'off_targets', off_target_count: 2 });
    // P1_L site at 4:7432183 is ignored: in no listed product (likely or unlikely) ...
    const near = (a) => a.region === '4' && (Math.abs(a.start - 7432183) <= 3 || Math.abs(a.end - 7432183) <= 3);
    should(pair.P1.off_targets.concat(pair.P1.unlikely).some(near)).be.false();
    // ... because the site re-aligns with >= ignore_mismatches edits (gap-aware DP 6; ungapped 7)
    const sequence = require('../../../api/helpers/primers/sequence');
    const realign = require('../../../api/helpers/primers/check/realign');
    const fasta = '/scratch/olson/fasta/sorghum_bicolor/dna/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel.fa.gz';
    const window = await sequence.fetch(fasta, '4', 7432183 - 22 - 3, 7432183 + 3, -1, { maxLength: 1000 });
    const aln = realign.alignPrimer(P.P1_L, window, window.length - 3, 5);
    t.diagnostic('P1_L at 4:7432183 (minus) DP edits ' + aln.mm + ' ops ' + aln.ops);
    should(aln.mm).be.aboveOrEqual(results.params.ignore_mismatches);
    const plus = (await sequence.fetch(fasta, '4', 7432183 - 22, 7432183, 1, { maxLength: 100 })).toUpperCase();
    const rc = realign.revcomp(P.P1_L);
    let ungapped = 0;
    for (let k = 0; k < rc.length; k++) if (plus[k] !== rc[k]) ungapped++;
    should(ungapped).equal(7);

    // P3: off-target 4:7436221-7436845 likely_weak
    should(pair.P3.on_target).match({ region: '4', start: 7422482, end: 7423061, size: 580 });
    const p3off = find(pair.P3.off_targets.concat(pair.P3.unlikely), '4', 7436221, 7436845);
    t.diagnostic('P3 §10.4 product ' + JSON.stringify(brief(p3off)));
    // the §10.4 likely_weak off-target has a 4-edit L site: now unlikely, so P3 is specific
    should(p3off).match({ likelihood: 'unlikely', left_mm: 4, right_mm: 1, right_mm_pos: [1], terminal_mismatch: true, approx: true });
    should(pair.P3).match({ verdict: 'specific', off_target_count: 0 });
    // those approximate 4-edit bounds are the exact re-aligned edits (the skipped FASTA re-alignment would not amplify either)
    const p1l = await sequence.fetch(fasta, '5', 66892235 - P.P1_L.length + 1 - 3, 66892235 + 3, -1, { maxLength: 1000 });
    should(realign.alignPrimer(P.P1_L, p1l, p1l.length - 3, 5)).match({ mm: 4, mm_pos: [20, 17, 14, 11] });
    const p3l = await sequence.fetch(fasta, '4', 7436221 - 3, 7436221 + P.P3_L.length - 1 + 3, 1, { maxLength: 1000 });
    should(realign.alignPrimer(P.P3_L, p3l, p3l.length - 3, 5)).match({ mm: 4, mm_pos: [22, 12, 6, 1] });

    // P5_L repetitive
    should(results.primers[P.P5_L]).match({ repetitive: true });
    should(results.primers[P.P5_L].near_perfect_sites).be.above(results.params.repeat_site_threshold);
    should(results.primers[P.P2_L]).match({ near_perfect_sites: 2, repetitive: false });

    // m20_7_14: 2-mismatch site at 4:7423537 found at word size 5
    should(find(pair.m20_7_14.off_targets.concat(pair.m20_7_14.on_target ? [pair.m20_7_14.on_target] : []), '4', 7423537, 7423746, 'LR')).match({ left_mm: 2, left_mm_pos: [14, 7], right_mm: 0 });
    // del18: 1-mismatch (gap) site kept, likely
    should(find(pair.del18.off_targets.concat(pair.del18.on_target ? [pair.del18.on_target] : []), '4', 7423537, 7423746, 'LR')).match({ left_mm: 1, likelihood: 'likely' });

    for (const seq of Object.keys(results.primers)) {
      should(results.primers[seq].sensitivity.reference).match({ word_size: 5 });
    }
  });

  it('qPCR transcript mode: J_L + P1_R on SORBI_3004G087700', { timeout: 600000 }, async (t) => {
    const { results, calls } = await check('qpcr_transcript', {
      system_name: 'sorghum_bicolor',
      mode: 'transcript',
      gene_id: 'SORBI_3004G087700',
      transcript_id: 'SORBI_3004G087700.3',
      pairs: [{ id: 'J', left: P.J_L, right: P.P1_R }]
    }, t);
    const g = results.specificity.pairs[0];
    t.diagnostic('genome ' + JSON.stringify([g.verdict, g.off_target_count, g.unlikely_count, g.gdna_products.length]));
    should(g.verdict).not.equal('on_target_missing');
    should(g.on_target).equal(null);
    assertStringency(g);
    const c = results.transcriptome.pairs[0];
    t.diagnostic('transcriptome ' + JSON.stringify([c.verdict, c.off_target_count, c.unlikely_count, c.off_targets.map((x) => [x.gene_id, x.likelihood, x.left_mm, x.right_mm])]));
    should(c.on_target.gene_id).equal('SORBI_3004G087700');
    should(c.on_target.isoforms.map((i) => [i.transcript_id, i.size]).sort()).eql([
      ['SORBI_3004G087700.1', 278], ['SORBI_3004G087700.2', 278], ['SORBI_3004G087700.3', 278]
    ]);
    should(c.off_target_count).equal(2);
    should(c.off_targets.map((x) => x.gene_id).sort()).eql(['SORBI_3004G087800', 'SORBI_3005G183900']);
    for (const x of c.off_targets) {
      should(x.left_mm).be.belowOrEqual(CAP);
      should(x.right_mm).be.belowOrEqual(CAP);
    }
    should(c.off_targets.find((x) => x.gene_id === 'SORBI_3005G183900').isoforms.map((i) => i.transcript_id).sort()).eql(['SORBI_3005G183900.1', 'SORBI_3005G183900.2']);
    // cDNA -max_target_seqs = max(5000, blastdbcmd -info count), fetched once
    const info = calls.filter((x) => x.args[0] === '-info');
    should(info).have.length(1);
    const cdna = calls.find((x) => x.args.indexOf('-db') >= 0 && /cdna\.all$/.test(x.args[x.args.indexOf('-db') + 1]) && x.args[0] !== '-info');
    should(cdna.args[cdna.args.indexOf('-max_target_seqs') + 1]).equal('47110');
    should(c.truncated).be.false();
  });

  it('pan-genome: P3 over sorghum_353, sorghum_grassl, sorghum_leoti', { timeout: 600000 }, async (t) => {
    const { results } = await check('pangenome_p3', {
      system_name: 'sorghum_bicolor',
      mode: 'gene',
      gene_id: 'SORBI_3004G087700',
      checks: ['pangenome'],
      genomes: ['sorghum_353', 'sorghum_grassl', 'sorghum_leoti'],
      params: { include_unlikely: true },
      pairs: [{ id: 'P3', left: P.P3_L, right: P.P3_R, expected: { region: '4', start: 7422482, end: 7423061 } }]
    }, t);
    const pair = results.pangenome.pairs[0];
    t.diagnostic('summary ' + JSON.stringify(pair.summary) + ' statuses ' + JSON.stringify(pair.genomes.map((x) => [x.system_name, x.status, x.other_amplicons, brief(x.primary)])));
    should(pair).match({ reference_size: 580, max_size: 4000 });
    const by = {};
    for (const x of pair.genomes) by[x.system_name] = x;
    should(by.sorghum_353.primary).match({ size: 580, size_delta: 0, terminal_mismatch: true, likelihood: 'likely_weak', left_mm_pos: [1], ortholog: true });
    should(by.sorghum_353.primary).match({ region: '4', start: 7499931, end: 7500510 });
    should(by.sorghum_grassl.primary).match({ size: 580, terminal_mismatch: true, likelihood: 'likely_weak' });
    should(by.sorghum_leoti.primary).match({ size: 581, size_delta: 1, right_mm: 1 });
    should(by.sorghum_353.ortholog_annotated).be.true();
    // the 626 bp (4/1-4/2 mm) and ~3.96 kb (5/4 mm) products no longer amplify, so no genome is 'multiple'
    for (const x of pair.genomes) should(x.status).not.equal('multiple');
    should(pair.genomes.map((x) => x.status)).eql(['single_mismatch', 'single_mismatch', 'single_mismatch']);
    should(pair.summary).match({ amplifies: 3, multiple: 0, truncated: 0 });
    const spec = results.specificity.pairs[0];
    assertStringency(spec);
  });

  it('presence/absence: SORBI_3001G046200 pair over sorghum_353, sorghum_is12661, sorghum_ji2731, sorghum_leoti', { timeout: 900000 }, async (t) => {
    const body = {
      system_name: 'sorghum_bicolor',
      mode: 'gene',
      gene_id: 'SORBI_3001G046200',
      checks: ['pangenome'],
      genomes: ['sorghum_353', 'sorghum_is12661', 'sorghum_ji2731', 'sorghum_leoti'],
      pairs: [{ id: 'PA', left: P.PA_L, right: P.PA_R }]
    };
    const { results, stats } = await check('presence_default_params', body, t);
    const spec = results.specificity.pairs[0];
    should(spec).match({ on_target_inferred: true, on_target: { region: '1', size: 3621 } });
    assertStringency(spec);
    const pair = results.pangenome.pairs[0];
    should(pair).match({ reference_size: 3621, max_size: 5932 });
    const by = {};
    for (const x of pair.genomes) by[x.system_name] = x;
    t.diagnostic('statuses ' + JSON.stringify(pair.genomes.map((x) => [x.system_name, x.status, x.other_amplicons, brief(x.primary), brief(x.nearest)])));
    t.diagnostic('realign per genome ' + JSON.stringify(Object.keys(stats.pangenome).map((s) => [s, stats.pangenome[s].realigned, stats.pangenome[s].bounded, stats.pangenome[s].realign_ms])));
    for (const sys of ['sorghum_353', 'sorghum_is12661']) {
      should(by[sys]).match({ status: 'no_amplicon', ortholog_annotated: false, primary: null });
      should(by[sys].nearest).be.an.Object();
      should(by[sys].nearest.likelihood).equal('unlikely');
      should(Math.max(by[sys].nearest.left_mm, by[sys].nearest.right_mm)).be.above(CAP);
    }
    should(by.sorghum_ji2731).match({ status: 'single_perfect', primary: { size_delta: -177, left_mm: 0, right_mm: 0 } });
    should(by.sorghum_leoti).match({ status: 'single_perfect', primary: { left_mm: 0, right_mm: 0, size_delta: 0, ortholog: true } });
    should(pair.summary).eql({ genomes_total: 4, single_perfect: 2, single_mismatch: 0, multiple: 0, no_amplicon: 2, db_unavailable: 0, error: 0, amplifies: 2, truncated: 0 });
    // most candidate sites are already over the cap by their lower bound and are not re-aligned
    for (const sys of Object.keys(stats.pangenome)) should(stats.pangenome[sys].bounded).be.above(stats.pangenome[sys].realigned);
    t.diagnostic('all timings ' + JSON.stringify(Object.keys(timings).map((k) => [k, timings[k].run_ms])));
  });
});
