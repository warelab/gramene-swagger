'use strict';

// Verdicts by mode × target (spec §B.8, §B.9, §B.12, §B.15).
// Part 1: pure specificity.js rules. Part 2: run() end to end on a synthetic species (fake BLAST that
// reports planted sites, real bgzip FASTA re-alignment through sequence.js, fake mongo).

const { describe, it, before, after } = require('node:test');
const fs = require('fs');
const os = require('os');
const path = require('path');
const should = require('should');

const specificity = require('../../../api/helpers/primers/check/specificity');
const { run, runWithStats, _internal: runInternal } = require('../../../api/helpers/primers/check/run');
const { fakeMongo } = require('../fixtures/verdicts/fake_mongo');
const { makeAmplicon, makeGroup } = require('../fixtures/verdicts/builders');
const W = require('../fixtures/verdicts/fake_world');

const E = { region: '4', start: 7423537, end: 7423746 };
const GENE_SPAN = { region: '4', start: 7421357, end: 7428285, strand: 1, map: 'GCA_000003195.3' };

function amp(region, start, end, o) {
  return makeAmplicon(Object.assign({ region, start, end }, o || {}));
}

describe('specificity: gene/region with expected (genome target)', () => {
  const on = amp('4', 7423537, 7423746);
  const para = amp('4', 7437317, 7437526);
  const chr5 = amp('5', 66890615, 66890824, { orientation: 'RL', left_mm_pos: [10, 7], right_mm_pos: [15, 10] });

  it('the on-target is an LR/RL product within ±2 bp of expected; everything else is an off-target', () => {
    const sel = specificity.selectGenome({ amplicons: [on, para, chr5], unlikely: [], mode: 'gene', expected: { region: '4', start: 7423539, end: 7423744 } });
    should(sel.verdict).equal('off_targets');
    should(sel.onTarget).equal(on);
    should(sel.inferred).be.false();
    should(sel.off).eql([para, chr5]);
    const block = specificity.genomePairBlock(sel, { id: 'P2' });
    should(block).match({ id: 'P2', verdict: 'off_targets', on_target_inferred: false, truncated: false, off_target_count: 2, unlikely_count: 0, gdna_products: [] });
    should(block.off_targets[1]).match({ region: '5', start: 66890615, end: 66890824, orientation: 'RL', left_mm: 2, right_mm: 2, left_mm_pos: [10, 7], right_mm_pos: [15, 10], left_3p_mm: 0, right_3p_mm: 0 });
  });

  it('3 bp away, or an LL/RR product at the expected coordinates, is not the on-target', () => {
    should(specificity.selectGenome({ amplicons: [on], mode: 'region', expected: { region: '4', start: 7423534, end: 7423746 } }).verdict).equal('on_target_missing');
    const rr = amp('4', 7423537, 7423746, { orientation: 'RR' });
    const sel = specificity.selectGenome({ amplicons: [rr], mode: 'gene', expected: E });
    should(sel.verdict).equal('on_target_missing');
    should(sel.off).eql([rr]);
    should(specificity.matchesExpected(amp('4', 7423537, 7423746), { region: 4, start: 7423537, end: 7423746 })).be.true();
  });

  it('specific / off_targets / on_target_missing / truncated precedence', () => {
    const v = (amps, truncated) => specificity.selectGenome({ amplicons: amps, mode: 'gene', expected: E, truncated }).verdict;
    should(v([on], false)).equal('specific');
    should(v([on, para], false)).equal('off_targets');
    should(v([para], false)).equal('on_target_missing');
    should(v([], false)).equal('on_target_missing');
    should(v([para], true)).equal('truncated');
    should(v([on], true)).equal('truncated');
    should(v([on, para], true)).equal('off_targets');
  });

  it('unlikely products are counted, listed only with include_unlikely, and never make the on-target', () => {
    const blocked = amp('4', 7423537, 7423746, { likelihood: 'unlikely', right_mm_pos: [2, 1] });
    const sel = specificity.selectGenome({ amplicons: [para], unlikely: [blocked], mode: 'gene', expected: E });
    should(sel.verdict).equal('on_target_missing');
    const plain = specificity.genomePairBlock(sel, { id: 'x' });
    should(plain.unlikely_count).equal(1);
    should(plain).not.have.property('unlikely');
    const listed = specificity.genomePairBlock(sel, { id: 'x', includeUnlikely: true });
    should(listed.unlikely).have.length(1);
    should(listed.unlikely[0]).match({ likelihood: 'unlikely', right_mm_pos: [2, 1] });
  });

  it('off_target_count is exact while listings are capped at 100; genes default to null', () => {
    const offs = [];
    for (let i = 0; i < 150; i++) offs.push(amp('7', 1000 + i * 10, 1300 + i * 10, { left_mm_pos: [12] }));
    const sel = specificity.selectGenome({ amplicons: [on].concat(offs), mode: 'gene', expected: E });
    const block = specificity.genomePairBlock(sel, { id: 'cap' });
    should(block.off_target_count).equal(150);
    should(block.off_targets).have.length(100);
    should(specificity.genomePairBlock(sel, { id: 'cap', maxListed: 5 }).off_targets).have.length(5);
    should(Object.keys(block)).eql(['id', 'verdict', 'on_target_inferred', 'truncated', 'on_target', 'off_target_count', 'unlikely_count', 'off_targets', 'gdna_products']);
    should(Object.keys(block.on_target)).eql(['region', 'start', 'end', 'size', 'strand', 'orientation', 'likelihood', 'left_mm', 'right_mm', 'left_3p_mm', 'right_3p_mm', 'left_mm_pos', 'right_mm_pos', 'terminal_mismatch', 'approx', 'genes']);
    should(block.on_target.genes).equal(null);
  });
});

describe('specificity: without expected, and sequence mode (genome target)', () => {
  it('exactly one perfect LR/RL product is the inferred on-target', () => {
    const perfect = amp('1', 3406188, 3409808);
    const mm = amp('1', 9000, 9500, { left_mm_pos: [9, 4] });
    const sel = specificity.selectGenome({ amplicons: [perfect, mm], mode: 'gene', expected: null });
    should(sel).match({ verdict: 'off_targets', inferred: true });
    should(sel.onTarget).equal(perfect);
    should(specificity.selectGenome({ amplicons: [perfect], mode: 'sequence' })).match({ verdict: 'specific', inferred: true });
    should(specificity.genomePairBlock(specificity.selectGenome({ amplicons: [perfect], mode: 'region' }), { id: 'r' })).match({ verdict: 'specific', on_target_inferred: true });
  });

  it('zero or several perfect products (or only LL/RR) → unverified_target, and every product is listed', () => {
    const a = amp('4', 7423537, 7423746);
    const b = amp('4', 7437317, 7437526);
    const two = specificity.selectGenome({ amplicons: [a, b], mode: 'region' });
    should(two).match({ verdict: 'unverified_target', inferred: false, onTarget: null });
    should(two.off).eql([a, b]);
    should(specificity.selectGenome({ amplicons: [amp('4', 1, 400, { orientation: 'LL' })], mode: 'sequence' }).verdict).equal('unverified_target');
    should(specificity.selectGenome({ amplicons: [], mode: 'sequence' }).verdict).equal('unverified_target');
    should(specificity.selectGenome({ amplicons: [], mode: 'sequence', truncated: true }).verdict).equal('truncated');
    should(specificity.selectGenome({ amplicons: [amp('4', 1, 400, { left_mm_pos: [3] })], mode: 'region' }).verdict).equal('unverified_target');
  });

  it('expected is honoured only in gene and region modes', () => {
    const pair = { expected: E };
    should(specificity.honoredExpected('gene', pair)).equal(E);
    should(specificity.honoredExpected('region', pair)).equal(E);
    should(specificity.honoredExpected('sequence', pair)).equal(null);
    should(specificity.honoredExpected('transcript', pair)).equal(null);
    should(specificity.honoredExpected('gene', {})).equal(null);
  });
});

describe('specificity: transcript mode, genome target', () => {
  it('products inside the gene span are gdna_products; the verdict is never on_target_missing', () => {
    const inside = amp('4', 7422291, 7422843, { left_mm_pos: [8, 2] });
    const edge = amp('4', 7421357, 7428285);
    const outside = amp('5', 66891530, 66892134, { orientation: 'RL' });
    const acrossEnd = amp('4', 7428000, 7428400);
    const sel = specificity.selectGenome({ amplicons: [inside, edge, outside, acrossEnd], mode: 'transcript', geneLocation: GENE_SPAN });
    should(sel.gdna).eql([inside, edge]);
    should(sel.off).eql([outside, acrossEnd]);
    should(sel.verdict).equal('off_targets');
    const block = specificity.genomePairBlock(sel, { id: 'J' });
    should(block).match({ on_target: null, on_target_inferred: false, off_target_count: 2 });
    should(block.gdna_products).have.length(2);
    for (const amps of [[], [inside], [outside]]) {
      for (const truncated of [false, true]) {
        should(specificity.selectGenome({ amplicons: amps, mode: 'transcript', geneLocation: GENE_SPAN, truncated }).verdict).not.equal('on_target_missing');
      }
    }
    should(specificity.selectGenome({ amplicons: [inside], mode: 'transcript', geneLocation: GENE_SPAN }).verdict).equal('specific');
    should(specificity.selectGenome({ amplicons: [], mode: 'transcript', geneLocation: GENE_SPAN }).verdict).equal('specific');
    should(specificity.selectGenome({ amplicons: [inside], mode: 'transcript', geneLocation: GENE_SPAN, truncated: true }).verdict).equal('truncated');
    // no gene location: every product is an off-target
    should(specificity.selectGenome({ amplicons: [inside], mode: 'transcript', geneLocation: null }).verdict).equal('off_targets');
  });
});

describe('specificity: transcript mode, cDNA target (gene groups)', () => {
  const G = 'SORBI_3004G087700';
  const on = makeGroup({ gene_id: G, sizes: [278, 278, 278] });
  const onLL = makeGroup({ gene_id: G, orientation: 'LL', sizes: [60] });
  const p1 = makeGroup({ gene_id: 'SORBI_3005G183900', sizes: [278, 278], left_mm: 1, left_mm_pos: [3] });
  const p2 = makeGroup({ gene_id: 'SORBI_3004G087800', sizes: [278], left_mm: 1, left_mm_pos: [14], right_mm: 2, right_mm_pos: [17, 4] });

  it('the on-target is the best LR/RL group of gene_id; other genes are off-targets', () => {
    const sel = specificity.selectCdna({ groups: [onLL, on, p1, p2], unlikelyGroups: [], geneId: G });
    should(sel.onTarget).equal(on);
    should(sel.off).eql([p1, p2]);
    should(sel.verdict).equal('off_targets');
    const block = specificity.cdnaPairBlock(sel, { id: 'J' });
    should(block).match({ id: 'J', verdict: 'off_targets', truncated: false, other_on_target_groups: 1, off_target_count: 2, unlikely_count: 0 });
    should(block.on_target.isoforms.map((i) => i.size)).eql([278, 278, 278]);
    should(block.off_targets.map((g) => g.gene_id)).eql(['SORBI_3005G183900', 'SORBI_3004G087800']);
    should(Object.keys(block)).eql(['id', 'verdict', 'truncated', 'on_target', 'other_on_target_groups', 'off_target_count', 'unlikely_count', 'off_targets']);
  });

  it('specific, on_target_missing, truncated', () => {
    should(specificity.selectCdna({ groups: [on], geneId: G }).verdict).equal('specific');
    should(specificity.selectCdna({ groups: [p1], geneId: G }).verdict).equal('on_target_missing');
    should(specificity.selectCdna({ groups: [], geneId: G }).verdict).equal('on_target_missing');
    should(specificity.selectCdna({ groups: [p1], geneId: G, truncated: true }).verdict).equal('truncated');
    should(specificity.selectCdna({ groups: [on], geneId: G, truncated: true }).verdict).equal('truncated');
    const blocked = makeGroup({ gene_id: G, likelihood: 'unlikely' });
    const withUnlikely = specificity.cdnaPairBlock(specificity.selectCdna({ groups: [on], unlikelyGroups: [blocked], geneId: G }), { id: 'u', includeUnlikely: true });
    should(withUnlikely.unlikely_count).equal(1);
    should(withUnlikely.unlikely[0].likelihood).equal('unlikely');
  });
});

describe('specificity.referenceSize and errorPairBlock', () => {
  it('uses the on-target size, the design transcript isoform, else size_min', () => {
    should(specificity.referenceSize({ on_target: { size: 580 } })).equal(580);
    should(specificity.referenceSize({ on_target: null })).equal(null);
    should(specificity.referenceSize(null)).equal(null);
    const group = specificity.publicGroup(makeGroup({ gene_id: 'G', isoforms: [{ transcript_id: 'G.1', start: 1, end: 300, size: 300 }, { transcript_id: 'G.2', start: 1, end: 250, size: 250 }] }));
    should(specificity.referenceSize({ on_target: group }, { transcriptId: 'G.1' })).equal(300);
    should(specificity.referenceSize({ on_target: group }, { transcriptId: 'G.9' })).equal(250);
    should(specificity.referenceSize({ on_target: group })).equal(250);
  });

  it('error blocks carry the target-specific keys', () => {
    should(specificity.errorPairBlock('genome', 'P', { code: 'X', message: 'y' })).match({ id: 'P', verdict: 'error', on_target_inferred: false, gdna_products: [], error: { code: 'X', message: 'y' } });
    should(specificity.errorPairBlock('cdna', 'P', null)).match({ verdict: 'error', other_on_target_groups: 0, error: null });
    should(specificity.VERDICTS).eql(['specific', 'off_targets', 'on_target_missing', 'unverified_target', 'truncated', 'error']);
  });
});

describe('run() on a synthetic reference (fake BLAST, real bgzip FASTA re-alignment, fake mongo)', () => {
  let world;
  let P;
  before(() => {
    world = W.standardWorld();
    P = world.primers;
  });
  after(() => world.cleanup());

  const geneBody = () => ({
    mode: 'gene',
    gene_id: 'REFG1',
    params: { include_unlikely: true },
    pairs: [
      { id: 'P', left: P.L, right: P.R, expected: { region: '1', start: 1001, end: 1300 } },
      { id: 'S', left: P.S_L.toLowerCase(), right: P.S_R, expected: { region: '8', start: 501, end: 700 } },
      { id: 'Q', left: P.Q_L, right: P.Q_R },
      { id: 'MISS', left: P.L, right: P.R, expected: { region: '1', start: 5001, end: 5300 } }
    ]
  });

  it('gene mode: verdicts per pair, primers block, annotation, §B.12 keys, partial after specificity', async () => {
    const { ctx, events } = W.worldCtx(world);
    const r = await run(W.request(geneBody()), ctx, { retryDelayMs: 0 });
    should(Object.keys(r)).eql(['engine', 'params', 'reference', 'sensitivity_note', 'primers', 'specificity', 'transcriptome', 'pangenome', 'warnings', 'timings_ms']);
    should(r.engine).eql({ algorithm_version: '2', blast: '2.13.0', reference: 'blastn-short r1 p-1 ws5 ungapped e30000 searchsp1.5e10', pangenome: null });
    should(r.params).eql({ max_product_size: 4000, ignore_mismatches: 6, min_total_mismatches: 2, min_3p_mismatches: 2, three_prime_window: 5, include_unlikely: true, repeat_site_threshold: 5, max_amplifying_mismatches: 3 });
    should(r.reference).eql({ system_name: 'ref', map_id: 'MAP_REF', total_bases: 27600 });
    should(r.transcriptome).equal(null);
    should(r.pangenome).equal(null);
    should(r.warnings).eql([]);
    should(r.timings_ms).have.properties(['reference', 'total']);
    should(Object.keys(r.primers).sort()).eql([P.L, P.R, P.S_L, P.S_R, P.Q_L, P.Q_R].sort());
    should(r.primers[P.L]).eql({ len: 22, near_perfect_sites: 3, repetitive: false, truncated: false, sensitivity: { reference: { word_size: 5, guaranteed_max_mismatches: 3 } } });

    const byId = {};
    for (const p of r.specificity.pairs) byId[p.id] = p;
    should(r.specificity.target).equal('genome');
    should(byId.P).match({ verdict: 'off_targets', on_target_inferred: false, truncated: false, off_target_count: 2, unlikely_count: 1 });
    should(byId.P.on_target).match({ region: '1', start: 1001, end: 1300, size: 300, orientation: 'LR', likelihood: 'likely', left_mm: 0, right_mm: 0, approx: false });
    should(byId.P.on_target.genes).eql([{ id: 'REFG1', name: 'REFG1', biotype: 'protein_coding', strand: 1 }]);
    should(byId.P.off_targets.map((a) => [a.region, a.start, a.end, a.orientation, a.likelihood, a.left_mm_pos, a.right_mm_pos])).eql([
      ['2', 1001, 1350, 'LR', 'likely', [16, 12], []],
      ['3', 1001, 1300, 'RL', 'likely_weak', [1], []]
    ]);
    should(byId.P.off_targets[0].genes.map((g) => g.id)).eql(['OTHER_G']);
    should(byId.P.off_targets[1].terminal_mismatch).be.true();
    should(byId.P.unlikely.map((a) => [a.region, a.likelihood, a.right_mm_pos, a.right_3p_mm])).eql([['4', 'unlikely', [2, 1], 2]]);
    // the 6-mismatch copy on chromosome 5 is ignored: neither listed nor counted
    should(byId.P.off_targets.concat(byId.P.unlikely).some((a) => a.region === '5')).be.false();
    should(byId.S).match({ verdict: 'specific', off_target_count: 0, on_target: { region: '8', start: 501, end: 700, genes: [{ id: 'REFG8', strand: -1 }] } });
    should(byId.Q).match({ verdict: 'unverified_target', on_target: null, on_target_inferred: false, off_target_count: 2 });
    should(byId.MISS).match({ verdict: 'on_target_missing', on_target: null, off_target_count: 3, unlikely_count: 1 });

    should(events.partials).have.length(1);
    should(events.partials[0].specificity.pairs).have.length(4);
    should(events.partials[0].pangenome).equal(null);
    should(events.progress[0]).eql({ done: 0, total: 1, stage: 'reference', running: ['ref'] });
    should(events.progress[events.progress.length - 1]).eql({ done: 1, total: 1, stage: 'reference', running: [] });

    const blastCalls = world.calls.filter((c) => c.cmd === '/fake/bin/blastn' && c.args[0] !== '-version');
    const args = blastCalls[blastCalls.length - 1].args;
    const arg = (k) => args[args.indexOf(k) + 1];
    should([arg('-task'), arg('-word_size'), arg('-num_threads'), arg('-max_target_seqs'), arg('-db'), arg('-outfmt')]).eql([
      'blastn-short', '5', '4', '5000', world.assemblies.ref.blastdb.dna, '6 qseqid sseqid qlen qstart qend sstart send sstrand mismatch'
    ]);
  });

  it('region mode infers the on-target without expected; sequence mode ignores expected', async () => {
    const pair = { id: 'P', left: P.L, right: P.R, expected: { region: '1', start: 5001, end: 5300 } };
    const region = await run(W.request({ mode: 'region', pairs: [{ id: 'P', left: P.L, right: P.R }] }), W.worldCtx(world).ctx, { retryDelayMs: 0 });
    should(region.specificity.pairs[0]).match({ verdict: 'off_targets', on_target_inferred: true, on_target: { region: '1', start: 1001, end: 1300 } });
    const seq = await run(W.request({ mode: 'sequence', pairs: [pair] }), W.worldCtx(world, { gene: null }).ctx, { retryDelayMs: 0 });
    should(seq.specificity.pairs[0]).match({ verdict: 'off_targets', on_target_inferred: true, on_target: { region: '1', start: 1001 } });
  });

  it('transcript mode: genome gdna_products and cDNA gene groups (§B.9) with TRANSCRIPT_GENE_UNMAPPED', async () => {
    const { ctx, events } = W.worldCtx(world);
    const before = world.calls.length;
    const r = await run(W.request({ mode: 'transcript', gene_id: 'REFG1', transcript_id: 'REFG1.2', pairs: [{ id: 'P', left: P.L, right: P.R }, { id: 'S', left: P.S_L, right: P.S_R }] }), ctx, { retryDelayMs: 0 });
    const [gp, gs] = r.specificity.pairs;
    should(gp).match({ verdict: 'off_targets', on_target: null, off_target_count: 2 });
    should(gp.gdna_products.map((a) => [a.region, a.start, a.end])).eql([['1', 1001, 1300]]);
    should(gs).match({ verdict: 'off_targets', off_target_count: 1 });
    should(r.transcriptome).match({ target: 'cdna', gene_id: 'REFG1', transcript_id: 'REFG1.2' });
    const [tp, ts] = r.transcriptome.pairs;
    should(tp).match({ id: 'P', verdict: 'off_targets', off_target_count: 3, other_on_target_groups: 0 });
    should(tp.on_target).match({ gene_id: 'REFG1', orientation: 'LR', size_min: 300, size_max: 300, left_mm: 0, right_mm: 0, approx: false });
    should(tp.on_target.isoforms.map((i) => i.transcript_id)).eql(['REFG1.1', 'REFG1.2']);
    should(tp.off_targets.map((g) => [g.gene_id, g.left_mm, g.right_mm, g.isoforms.map((i) => i.transcript_id)])).eql([
      ['ODDG', 0, 0, ['odd_tx']],
      ['OTHER_G', 1, 0, ['OTHER_G.1']],
      ['UNMAPPED_X', 0, 1, ['UNMAPPED_X.3']]
    ]);
    should(ts).match({ id: 'S', verdict: 'on_target_missing', on_target: null, off_target_count: 1 });
    should(r.warnings.map((w) => w.code)).eql(['TRANSCRIPT_GENE_UNMAPPED']);
    should(r.timings_ms).have.properties(['reference', 'transcriptome', 'total']);
    should(events.partials.map((x) => [!!x.specificity, !!x.transcriptome])).eql([[true, false], [true, true]]);
    should(events.progress.map((p) => p.stage + ':' + p.done + '/' + p.total)).eql(['reference:0/2', 'reference:1/2', 'transcriptome:1/2', 'transcriptome:2/2']);
    const mine = world.calls.slice(before);
    const info = mine.filter((c) => c.args[0] === '-info');
    should(info).have.length(1);
    should(info[0].cmd).equal('/fake/bin/blastdbcmd');
    const cdnaCall = mine.find((c) => c.args.indexOf(world.assemblies.ref.blastdb.cdna) >= 0 && c.args[0] !== '-info');
    should(cdnaCall.args[cdnaCall.args.indexOf('-outfmt') + 1]).equal('6 qseqid sseqid qlen qstart qend sstart send sstrand mismatch qseq sseq');
    should(cdnaCall.args[cdnaCall.args.indexOf('-word_size') + 1]).equal('5');
  });

  it('a failed reference BLAST is retried once; a second failure throws REFERENCE_BLAST_FAILED without paths', async () => {
    const db = world.dbs.get(world.assemblies.ref.blastdb.dna);
    const body = { mode: 'gene', gene_id: 'REFG1', pairs: [{ id: 'S', left: P.S_L, right: P.S_R, expected: { region: '8', start: 501, end: 700 } }] };
    db.failuresLeft = 1;
    let before = world.blastCalls('ref').length;
    const ok = await run(W.request(body), W.worldCtx(world).ctx, { retryDelayMs: 0 });
    should(ok.specificity.pairs[0].verdict).equal('specific');
    should(world.blastCalls('ref').length - before).equal(2);

    db.failuresLeft = 2;
    before = world.blastCalls('ref').length;
    const err = await run(W.request(body), W.worldCtx(world).ctx, { retryDelayMs: 0 }).then(() => null, (e) => e);
    should(err).be.an.Error();
    should(err.code).equal('REFERENCE_BLAST_FAILED');
    should(err.details).eql({ target: 'genome', cause: 'BLAST_FAILED' });
    should(err.message).match(/reference genome BLAST failed: blastn exited with code 2/);
    should(err.message.indexOf(world.root)).equal(-1);
    should(err.message).match(/ref\.dna\.toplevel\.nal/);
    should(world.blastCalls('ref').length - before).equal(2);
    db.failuresLeft = 0;
  });

  it('abort rejects with the signal reason when it has a code, else with ABORTED', async () => {
    const body = W.request({ mode: 'region', pairs: [{ id: 'S', left: P.S_L, right: P.S_R }] });
    const reason = Object.assign(new Error('check exceeded 30 min'), { code: 'JOB_TIMEOUT' });
    const a = W.worldCtx(world, { onProgress: (p, controller) => { if (p.stage === 'reference') controller.abort(reason); } });
    const e1 = await run(body, a.ctx, { retryDelayMs: 0 }).then(() => null, (e) => e);
    should(e1).equal(reason);
    const b = W.worldCtx(world, { onProgress: (p, controller) => { if (p.stage === 'reference') controller.abort(); } });
    const e2 = await run(body, b.ctx, { retryDelayMs: 0 }).then(() => null, (e) => e);
    should(e2.code).equal('ABORTED');
    const c = W.worldCtx(world);
    c.controller.abort();
    should((await run(body, c.ctx).then(() => null, (e) => e)).code).equal('ABORTED');
  });

  it('mongo unavailable: ANNOTATION_UNAVAILABLE, genes null, verdicts unchanged', async () => {
    const { ctx } = W.worldCtx(world, { mongo: world.mongo({ noCollection: true }).mongo });
    const r = await run(W.request(geneBody()), ctx, { retryDelayMs: 0 });
    should(r.warnings.map((w) => w.code)).eql(['ANNOTATION_UNAVAILABLE']);
    should(r.specificity.pairs.map((p) => p.verdict)).eql(['off_targets', 'specific', 'unverified_target', 'on_target_missing']);
    should(r.specificity.pairs[0].on_target.genes).equal(null);
    should(r.specificity.pairs[0].off_targets.every((a) => a.genes === null)).be.true();
  });

  it('NO_BLASTDB without a reference database; TypeError on malformed input', async () => {
    const { ctx } = W.worldCtx(world);
    ctx.resolved = { assemblies: { ref: Object.assign({}, world.assemblies.ref, { blastdb: { dna: null, cdna: null } }) } };
    should((await run(W.request({ pairs: [{ id: 'S', left: P.S_L, right: P.S_R }] }), ctx).then(() => null, (e) => e)).code).equal('NO_BLASTDB');
    ctx.resolved = { assemblies: { ref: Object.assign({}, world.assemblies.ref, { blastdb: { dna: world.assemblies.ref.blastdb.dna, cdna: null } }) } };
    should((await run(W.request({ mode: 'transcript', gene_id: 'REFG1', pairs: [{ id: 'S', left: P.S_L, right: P.S_R }] }), ctx).then(() => null, (e) => e)).code).equal('NO_BLASTDB');
    should(await run(null, ctx).then(() => null, (e) => e)).be.instanceOf(TypeError);
    should(await run(W.request({ pairs: [] }), W.worldCtx(world).ctx).then(() => null, (e) => e)).be.instanceOf(TypeError);
  });

  it('repeat_site_threshold drives primers[].repetitive', async () => {
    const r = await run(W.request({ mode: 'region', params: { repeat_site_threshold: 2 }, pairs: [{ id: 'P', left: P.L, right: P.R }] }), W.worldCtx(world).ctx, { retryDelayMs: 0 });
    should(r.primers[P.L]).match({ near_perfect_sites: 3, repetitive: true });
    should(r.primers[P.R]).match({ near_perfect_sites: 4, repetitive: true });
  });
});

// Mongo whose overlap or transcript queries can be switched to fail mid-job ('Topology was destroyed').
function switchableMongo(docs) {
  const base = fakeMongo(docs || []);
  const state = { failOverlap: false, failTranscripts: false };
  const failing = () => {
    const cursor = { limit: () => cursor, toArray: () => Promise.reject(new Error('Topology was destroyed')) };
    return cursor;
  };
  const collection = {
    find(query, options) {
      const overlap = Object.prototype.hasOwnProperty.call(query, 'location.map');
      return (overlap ? state.failOverlap : state.failTranscripts) ? failing() : base.collection.find(query, options);
    }
  };
  return { state, mongo: { genes: { mongoCollection: async () => collection } } };
}

describe('run(): stringency, identical primers, BLAST timeouts, private cwd and fatal mongo', () => {
  let world;
  const P = W.PRIMERS;
  before(() => {
    const chr = (name, len, seed) => new W.Seq(name, len, seed);
    world = W.buildWorld({
      genomes: {
        ref: {
          taxon_id: 1001,
          map_id: 'MAP_REF',
          chromosomes: [
            chr('1', 3000, 201).lr(1001, 1300, P.L, P.R),
            chr('2', 3000, 202).lr(1001, 1400, P.L, P.R, { l: [20, 17, 14, 11] }), // 4 mismatches in L, none in the 3' window
            chr('3', 3000, 203).lr(1001, 1500, P.L, P.R, { l: [18, 12, 8] }), // 3 mismatches in L
            chr('4', 3000, 204).forward(501, P.S_L).reverse(800, P.S_L) // one primer on both strands
          ],
          cdna: [chr('G1.1', 800, 205).lr(101, 400, P.L, P.R)]
        },
        pan1: { taxon_id: 1002, map_id: 'MAP_P1', chromosomes: [chr('1', 3000, 206).lr(1001, 1300, P.L, P.R)] }
      }
    });
  });
  after(() => world.cleanup());

  const ctxOf = (o) => W.worldCtx(world, Object.assign({ gene: null, mongo: fakeMongo([]).mongo }, o || {}));
  const pairP = { id: 'P', left: P.L, right: P.R, expected: { region: '1', start: 1001, end: 1300 } };

  it('a product amplifies only with <= max_amplifying_mismatches per primer; sites over the cap are not re-aligned', async () => {
    const body = (params) => W.request({ mode: 'region', params: Object.assign({ include_unlikely: true }, params), pairs: [pairP] });
    const def = await runWithStats(body({}), ctxOf().ctx, { retryDelayMs: 0 });
    const p = def.results.specificity.pairs[0];
    should(def.results.params.max_amplifying_mismatches).equal(3);
    should(p).match({ verdict: 'off_targets', off_target_count: 1, unlikely_count: 1, on_target: { region: '1', start: 1001, end: 1300 } });
    should(p.off_targets[0]).match({ region: '3', start: 1001, end: 1500, left_mm: 3, left_mm_pos: [18, 12, 8], right_mm: 0, likelihood: 'likely', approx: false });
    // the 4-mismatch copy: listed as unlikely with an approximate lower bound, its L site never re-aligned
    should(p.unlikely[0]).match({ region: '2', start: 1001, end: 1400, left_mm: 4, left_mm_pos: null, likelihood: 'unlikely', approx: true });
    should(def.stats.reference).match({ bounded: 1 });

    const four = await runWithStats(body({ max_amplifying_mismatches: 4 }), ctxOf().ctx, { retryDelayMs: 0 });
    const q = four.results.specificity.pairs[0];
    should(q).match({ verdict: 'off_targets', off_target_count: 2, unlikely_count: 0 });
    should(q.off_targets.find((a) => a.region === '2')).match({ left_mm: 4, left_mm_pos: [20, 17, 14, 11], likelihood: 'likely', approx: false });
    should(four.stats.reference.bounded).equal(0);

    const two = (await run(body({ max_amplifying_mismatches: 2 }), ctxOf().ctx, { retryDelayMs: 0 })).specificity.pairs[0];
    should(two).match({ verdict: 'specific', off_target_count: 0, unlikely_count: 2 });
    should(two.unlikely.map((a) => [a.region, a.left_mm, a.approx])).eql([['3', 3, true], ['2', 4, true]]);
  });

  it('an identical-primer pair can match its expected product (chk-identical-lr)', async () => {
    const withExpected = await run(W.request({ mode: 'region', pairs: [{ id: 'X', left: P.S_L, right: P.S_L.toLowerCase(), expected: { region: '4', start: 501, end: 800 } }] }), ctxOf().ctx, { retryDelayMs: 0 });
    should(withExpected.specificity.pairs[0]).match({ verdict: 'specific', off_target_count: 0, on_target: { region: '4', start: 501, end: 800, size: 300, orientation: 'LR', strand: 1, left_mm: 0, right_mm: 0 } });
    const inferred = await run(W.request({ mode: 'sequence', pairs: [{ id: 'X', left: P.S_L, right: P.S_L }] }), ctxOf().ctx, { retryDelayMs: 0 });
    should(inferred.specificity.pairs[0]).match({ verdict: 'specific', on_target_inferred: true, on_target: { orientation: 'LR' } });
  });

  it('a timed-out reference BLAST is not retried; a failed one still is (chk-retry-timeout)', async () => {
    const db = world.dbs.get(world.assemblies.ref.blastdb.dna);
    const body = W.request({ mode: 'region', pairs: [pairP] });
    db.timeoutsLeft = 1;
    let before = world.blastCalls('ref').length;
    const err = await run(body, ctxOf().ctx, { retryDelayMs: 0 }).then(() => null, (e) => e);
    should(err).match({ code: 'REFERENCE_BLAST_FAILED', details: { target: 'genome', cause: 'BLAST_TIMEOUT' } });
    should(world.blastCalls('ref').length - before).equal(1);
    db.failuresLeft = 1;
    before = world.blastCalls('ref').length;
    should((await run(body, ctxOf().ctx, { retryDelayMs: 0 })).specificity.pairs[0].verdict).equal('off_targets');
    should(world.blastCalls('ref').length - before).equal(2);
  });

  it('every blastn / blastdbcmd child runs in a private cwd, never the shared tmp directory (sec-ncbirc-cwd-tmp)', async () => {
    runInternal.clearVersionCache();
    const a = ctxOf();
    let before = world.calls.length;
    await run(W.request({ mode: 'transcript', gene_id: 'G1', checks: ['pangenome'], genomes: ['pan1'], pairs: [{ id: 'P', left: P.L, right: P.R }] }), a.ctx, { retryDelayMs: 0 });
    const mine = world.calls.slice(before);
    should(mine.map((c) => (c.args[0] === '-version' ? 'version' : c.args[0] === '-info' ? 'info' : 'blastn'))).containDeep(['version', 'info', 'blastn']);
    const dir = a.ctx.tmpdir();
    should(mine.every((c) => c.cwd === dir)).be.true();
    should(fs.statSync(dir).mode & 0o777).equal(0o700);

    // ctx.tmpdir() returning the world-writable system tmp: run() makes (and removes) its own private directory
    runInternal.clearVersionCache();
    const b = ctxOf({ tmpdir: () => os.tmpdir() });
    const modes = [];
    const spawn = b.ctx.spawnLines;
    b.ctx.spawnLines = (cmd, args, o) => {
      modes.push(fs.statSync(o.cwd).mode & 0o777);
      return spawn(cmd, args, o);
    };
    before = world.calls.length;
    await run(W.request({ mode: 'region', pairs: [pairP] }), b.ctx, { retryDelayMs: 0 });
    const cwds = Array.from(new Set(world.calls.slice(before).map((c) => c.cwd)));
    should(cwds).have.length(1);
    should(cwds[0]).startWith(path.join(os.tmpdir(), 'primers-check-'));
    should(modes.length).be.aboveOrEqual(2);
    should(modes.every((m) => m === 0o700)).be.true();
    should(fs.existsSync(cwds[0])).be.false();
    should(runInternal.isPrivateDir(os.tmpdir())).be.false();
  });

  it('bounded (not re-aligned) sites keep coordinates inside the sequence (chk-coords-below-one)', async () => {
    const lengths = [];
    const job = new runInternal.CheckRun(W.request({ mode: 'region', pairs: [pairP] }), ctxOf().ctx, {
      retryDelayMs: 0,
      fetchWindow: async () => { throw new Error('bounded sites must not be fetched'); },
      regionLength: async (fp, region) => { lengths.push(region); return 1000; }
    });
    const site = (primer, face, p5, p3, mm) => ({ key: primer + ':c:' + face + ':' + p5, primer, len: primer.length, subject: 'c', face, strand: face === 'F' ? 1 : -1, p5, p3, t5: 0, t3: 0, mm });
    const fwd = site(P.L, 'F', -2, 19, 4); // 5' end extrapolated before the contig start
    const rev = site(P.R, 'R', 1003, 983, 5); // 5' end extrapolated past the contig end (length 1000)
    const out = await job.alignSites(world.assemblies.ref, 'genome', [{ subject: 'c', orientation: 'LR', fwd, rev }]);
    should(out).match({ realigned: 0, bounded: 2, truncated: false, fetchErrors: 0 });
    should(out.byKey.get(fwd.key)).match({ p5: 1, p3: 19, mm: 4, approx: true });
    should(out.byKey.get(rev.key)).match({ p5: 1000, p3: 983, mm: 5, approx: true });
    should(lengths).eql(['c']);
  });

  it('with ctx.fatalOnMongoUnavailable a failed mongo is a fatal MONGO_UNAVAILABLE in every stage (ops-mongo-fatal-unreachable)', async () => {
    const body = W.request({ mode: 'region', pairs: [pairP] });
    const fatal = (ctx) => Object.assign(ctx, { fatalOnMongoUnavailable: true });
    // reference stage, collection unavailable
    const noColl = fatal(ctxOf({ mongo: fakeMongo([], { noCollection: true }).mongo }).ctx);
    should(await run(body, noColl, { retryDelayMs: 0 }).then(() => null, (e) => e)).match({ code: 'MONGO_UNAVAILABLE', fatal: true });
    // without the flag, or without a mongo handle, the check completes with the warning
    should((await run(body, ctxOf({ mongo: fakeMongo([], { noCollection: true }).mongo }).ctx, { retryDelayMs: 0 })).warnings.map((w) => w.code)).eql(['ANNOTATION_UNAVAILABLE']);
    should((await run(body, fatal(ctxOf({ mongo: null }).ctx), { retryDelayMs: 0 })).warnings.map((w) => w.code)).eql(['ANNOTATION_UNAVAILABLE']);

    // transcriptome stage: overlaps work, transcript → gene mapping fails
    const tx = switchableMongo([]);
    tx.state.failTranscripts = true;
    const t = ctxOf({ mongo: tx.mongo });
    fatal(t.ctx);
    const te = await run(W.request({ mode: 'transcript', gene_id: 'G1', pairs: [{ id: 'P', left: P.L, right: P.R }] }), t.ctx, { retryDelayMs: 0 }).then(() => null, (e) => e);
    should(te).match({ code: 'MONGO_UNAVAILABLE', fatal: true });
    should(t.events.partials).have.length(1); // the genome stage had completed

    // pan-genome stage: mongo goes away after the reference stage; the pool reports the error, not ABORTED
    const pan = switchableMongo([]);
    const g = ctxOf({ mongo: pan.mongo, onProgress: (pr) => { if (pr.stage === 'pangenome') pan.state.failOverlap = true; } });
    fatal(g.ctx);
    const ge = await run(W.request({ mode: 'region', checks: ['pangenome'], genomes: ['pan1'], pairs: [pairP] }), g.ctx, { retryDelayMs: 0 }).then(() => null, (e) => e);
    should(ge).match({ code: 'MONGO_UNAVAILABLE', fatal: true });
  });
});
