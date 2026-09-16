'use strict';

// Pan-genome statuses, primary ordering, summary invariant and maxSize (spec §B.10, plan Overrides),
// plus run() over a synthetic species (genome and cDNA targets, retry, error, db_unavailable).

const { describe, it, before, after } = require('node:test');
const should = require('should');

const pangenome = require('../../../api/helpers/primers/check/pangenome');
const { run } = require('../../../api/helpers/primers/check/run');
const { makeAmplicon, makeGroup } = require('../fixtures/verdicts/builders');
const W = require('../fixtures/verdicts/fake_world');

function amp(o) {
  return makeAmplicon(o);
}

const INFO = { system_name: 'sorghum_353', display_name: 'Sb verticilliflorum 353', orthologAnnotated: true };

describe('pangenome.summarize', () => {
  it('counts statuses, amplifies = single_perfect + single_mismatch + multiple, and counts sum to genomes_total', () => {
    const statuses = ['single_perfect', 'single_perfect', 'single_mismatch', 'multiple', 'no_amplicon', 'db_unavailable', 'error', 'error'];
    const s = pangenome.summarize(statuses.map((status) => ({ status })));
    should(s).eql({ genomes_total: 8, single_perfect: 2, single_mismatch: 1, multiple: 1, no_amplicon: 1, db_unavailable: 1, error: 2, amplifies: 4, truncated: 0 });
    should(pangenome.STATUSES.reduce((a, k) => a + s[k], 0)).equal(s.genomes_total);
    should(pangenome.summarize([])).eql({ genomes_total: 0, single_perfect: 0, single_mismatch: 0, multiple: 0, no_amplicon: 0, db_unavailable: 0, error: 0, amplifies: 0, truncated: 0 });
    should(pangenome.summarize([null, { status: 'multiple' }]).genomes_total).equal(1);
    should(() => pangenome.summarize([{ status: 'amplifies' }])).throw(/unknown pan-genome status/);
  });

  it('counts truncated genomes (lower-bound statuses) without changing the status sum (chk-pan-truncated-status)', () => {
    const s = pangenome.summarize([{ status: 'no_amplicon', truncated: true }, { status: 'single_perfect', truncated: true }, { status: 'multiple' }, { status: 'error', truncated: false }]);
    should(s).eql({ genomes_total: 4, single_perfect: 1, single_mismatch: 0, multiple: 1, no_amplicon: 1, db_unavailable: 0, error: 1, amplifies: 2, truncated: 2 });
    should(pangenome.STATUSES.reduce((a, k) => a + s[k], 0)).equal(s.genomes_total);
    const entry = pangenome.genomeEntry({ system_name: 's', products: [], unlikely: [], truncated: true });
    should(entry).match({ status: 'no_amplicon', truncated: true });
    should(pangenome.summarize([entry, pangenome.errorEntry({ system_name: 'e' })])).match({ no_amplicon: 1, error: 1, truncated: 1 });
  });

  it('the invariant holds for random mixes', () => {
    let seed = 7;
    const rand = () => (seed = (seed * 1103515245 + 12345) % 2147483648) / 2147483648;
    for (let t = 0; t < 50; t++) {
      const n = Math.floor(rand() * 20);
      const entries = [];
      for (let i = 0; i < n; i++) entries.push({ status: pangenome.STATUSES[Math.floor(rand() * pangenome.STATUSES.length)] });
      const s = pangenome.summarize(entries);
      should(pangenome.STATUSES.reduce((a, k) => a + s[k], 0)).equal(n);
      should(s.amplifies).equal(s.single_perfect + s.single_mismatch + s.multiple);
    }
  });
});

describe('pangenome.genomeEntry statuses (§B.10)', () => {
  it('no likely product → no_amplicon with the best unlikely product as nearest', () => {
    const near = amp({ likelihood: 'unlikely', right_mm_pos: [2, 1], terminal_mismatch: true });
    const far = amp({ likelihood: 'unlikely', left_mm_pos: [14, 9, 2], right_mm_pos: [5, 1], start: 1, end: 900 });
    const e = pangenome.genomeEntry(Object.assign({ products: [], unlikely: [far, near], referenceSize: 210 }, INFO));
    should(e).match({ system_name: 'sorghum_353', display_name: 'Sb verticilliflorum 353', status: 'no_amplicon', ortholog_annotated: true, primary: null, other_amplicons: 0, others: [] });
    should(e.nearest).match({ likelihood: 'unlikely', right_mm: 2, size_delta: 0, mismatch_in_3p_window: true, terminal_mismatch: true, ortholog: null });
    should(pangenome.genomeEntry(Object.assign({ products: [], unlikely: [] }, INFO)).nearest).equal(null);
  });

  it('one perfect product → single_perfect; one with mismatches → single_mismatch with 3′ flags', () => {
    const perfect = pangenome.genomeEntry(Object.assign({ products: [amp({})], referenceSize: 210 }, INFO));
    should(perfect).match({ status: 'single_perfect', other_amplicons: 0, others: [], nearest: null });
    should(perfect.primary).match({ size: 210, size_delta: 0, mismatch_in_3p_window: false, terminal_mismatch: false });

    const weak = amp({ start: 7499931, end: 7500510, left_mm_pos: [1], likelihood: 'likely_weak', terminal_mismatch: true });
    const e = pangenome.genomeEntry(Object.assign({ products: [weak], referenceSize: 580 }, INFO));
    should(e.status).equal('single_mismatch');
    should(e.primary).match({ start: 7499931, end: 7500510, size: 580, size_delta: 0, likelihood: 'likely_weak', left_mm: 1, left_3p_mm: 1, left_mm_pos: [1], mismatch_in_3p_window: true, terminal_mismatch: true });

    const outside = pangenome.genomeEntry(Object.assign({ products: [amp({ start: 1, end: 581, right_mm_pos: [18] })], referenceSize: 580 }, INFO));
    should(outside.status).equal('single_mismatch');
    should(outside.primary).match({ size: 581, size_delta: 1, right_mm: 1, right_3p_mm: 0, mismatch_in_3p_window: false, terminal_mismatch: false });
  });

  it('two or more → multiple; primary by ortholog, then mismatches, then |size − reference|; others capped at 10', () => {
    const products = [];
    for (let i = 0; i < 14; i++) products.push(amp({ region: 'x' + i, start: 1, end: 600 + i }));
    const mm = amp({ region: 'mm', start: 1, end: 580, left_mm_pos: [9] });
    const orth = amp({ region: 'orth', start: 1, end: 900, left_mm_pos: [9, 8] });
    orth.ortholog = true;
    const all = products.concat([mm, orth]);
    all.forEach((p) => { if (p.ortholog === undefined) p.ortholog = false; });
    const e = pangenome.genomeEntry(Object.assign({ products: all, referenceSize: 600 }, INFO));
    should(e.status).equal('multiple');
    should(e.primary).match({ region: 'orth', ortholog: true, size_delta: 300 });
    should(e.others.map((o) => o.region)).eql(['x0', 'x1', 'x2', 'x3', 'x4', 'x5', 'x6', 'x7', 'x8', 'x9']);
    should(e.other_amplicons).equal(15);

    // without an ortholog: fewest mismatches first, then closest size
    const f = pangenome.genomeEntry(Object.assign({ products: [amp({ region: 'far', start: 1, end: 900 }), amp({ region: 'near', start: 1, end: 605 }), mm] , referenceSize: 600 }, INFO));
    should([f.primary.region].concat(f.others.map((o) => o.region))).eql(['near', 'far', 'mm']);
    // without a reference size the size is not a criterion
    const g = pangenome.genomeEntry(Object.assign({ products: [amp({ region: 'b', start: 1, end: 900 }), amp({ region: 'a', start: 1, end: 605 })] }, INFO));
    should(g.primary).match({ region: 'a', size_delta: null });
  });

  it('entry shapes for db_unavailable and error', () => {
    const u = pangenome.unavailableEntry(Object.assign({ reason: 'AMBIGUOUS_ASSEMBLY' }, INFO));
    should(u).eql({ system_name: 'sorghum_353', display_name: 'Sb verticilliflorum 353', status: 'db_unavailable', ortholog_annotated: true, truncated: false, primary: null, other_amplicons: 0, others: [], nearest: null, reason: 'AMBIGUOUS_ASSEMBLY' });
    const e = pangenome.errorEntry({ system_name: 's', error: { code: 'BLAST_TIMEOUT', message: 'BLAST timed out' } });
    should(e).match({ system_name: 's', display_name: 's', status: 'error', ortholog_annotated: null, error: { code: 'BLAST_TIMEOUT', message: 'BLAST timed out' } });
    should(pangenome.errorEntry({ system_name: 's' }).error).eql({ code: 'BLAST_FAILED', message: 'BLAST failed' });
    should(Object.keys(pangenome.genomeEntry(Object.assign({ products: [amp({})] }, INFO)))).eql(['system_name', 'display_name', 'status', 'ortholog_annotated', 'truncated', 'primary', 'other_amplicons', 'others', 'nearest']);
  });
});

describe('pangenome products: orthologs, cDNA groups, maxSize', () => {
  it('markOrthologs: overlapping gene in the set (genome), gene_id in the set (cDNA), null when unknown', () => {
    const ids = new Set(['353.004G093500']);
    const hit = amp({ genes: [{ id: '353.004G093500' }, { id: 'other' }] });
    const miss = amp({ genes: [{ id: 'other' }] });
    const empty = amp({ genes: [] });
    const unknown = amp({ genes: null });
    const group = makeGroup({ gene_id: '353.004G093500' });
    pangenome.markOrthologs([hit, miss, empty, unknown, group], ids);
    should([hit, miss, empty, unknown, group].map((p) => p.ortholog)).eql([true, false, false, null, true]);
    pangenome.markOrthologs([hit, group], null);
    should([hit.ortholog, group.ortholog]).eql([null, null]);
  });

  it('cDNA groups use the isoform size closest to the reference size', () => {
    const g = makeGroup({ gene_id: 'G', isoforms: [{ transcript_id: 'G.1', start: 1, end: 250, size: 250 }, { transcript_id: 'G.2', start: 1, end: 290, size: 290 }] });
    should(pangenome.productSize(g, 278)).equal(290);
    should(pangenome.sizeDelta(g, 278)).equal(12);
    should(pangenome.productSize(g, null)).equal(250);
    should(pangenome.sizeDelta(g, null)).equal(null);
    const e = pangenome.genomeEntry({ system_name: 's', products: [g], referenceSize: 278, orthologAnnotated: false });
    should(e.primary).match({ gene_id: 'G', size: 290, size_delta: 12, size_min: 250, size_max: 290, ortholog: null });
  });

  it('pan-genome maxSize = max(max_product_size, ceil(1.5 × reference_size) + 500)', () => {
    should(pangenome.pangenomeMaxSize(4000, 580)).equal(4000);
    should(pangenome.pangenomeMaxSize(4000, 3621)).equal(5932);
    should(pangenome.pangenomeMaxSize(4000, 3001)).equal(5002);
    should(pangenome.pangenomeMaxSize(1000, 581)).equal(1372);
    should(pangenome.pangenomeMaxSize(4000, null)).equal(4000);
  });
});

describe('run() pan-genome over a synthetic species', () => {
  let world;
  let P;
  before(() => {
    world = W.standardWorld();
    P = world.primers;
  });
  after(() => world.cleanup());

  const GENOMES = ['pan_x', 'pan_h', 'pan_g', 'pan_f', 'pan_e', 'pan_d', 'pan_c', 'pan_b', 'pan_a'];

  it('genome target: every status, primary ordering, retry, error, maxSize raise, partial invariants', async () => {
    const { ctx, events } = W.worldCtx(world, { errors: { pan_x: { code: 'AMBIGUOUS_ASSEMBLY', message: 'two assemblies' } } });
    const callsG = world.blastCalls('pan_g').length;
    const callsF = world.blastCalls('pan_f').length;
    const r = await run(W.request({
      mode: 'gene',
      gene_id: 'REFG1',
      checks: ['pangenome'],
      genomes: GENOMES,
      pairs: [
        { id: 'P', left: P.L, right: P.R, expected: { region: '1', start: 1001, end: 1300 } },
        { id: 'T', left: P.T_L, right: P.T_R }
      ]
    }), ctx, { retryDelayMs: 0 });

    should(r.engine.pangenome).equal('blastn-short r1 p-1 ws6 ungapped e30000 searchsp1.5e10');
    should(r.primers[P.L].sensitivity).eql({ reference: { word_size: 5, guaranteed_max_mismatches: 3 }, pangenome: { word_size: 6, guaranteed_max_mismatches: 2 } });
    should(r.pangenome.target).equal('genome');
    const [pp, pt] = r.pangenome.pairs;
    should(pp).match({ id: 'P', reference_size: 300, max_size: 4000 });
    should(pp.genomes.map((g) => g.system_name)).eql(GENOMES.slice().sort());
    const st = {};
    for (const g of pp.genomes) st[g.system_name] = g;
    should(pp.genomes.map((g) => g.status)).eql(['single_perfect', 'single_mismatch', 'multiple', 'no_amplicon', 'db_unavailable', 'error', 'single_perfect', 'single_perfect', 'db_unavailable']);
    should(pp.summary).eql({ genomes_total: 9, single_perfect: 3, single_mismatch: 1, multiple: 1, no_amplicon: 1, db_unavailable: 2, error: 1, amplifies: 5, truncated: 0 });

    should(st.pan_a).match({ display_name: 'PAN_A', ortholog_annotated: true, other_amplicons: 0 });
    should(st.pan_a.primary).match({ region: '1', start: 1101, end: 1400, size: 300, size_delta: 0, ortholog: true, genes: [{ id: 'PANA_G1' }] });
    should(st.pan_b).match({ ortholog_annotated: false });
    should(st.pan_b.primary).match({ start: 501, end: 810, size_delta: 10, right_mm: 1, right_mm_pos: [3], right_3p_mm: 1, mismatch_in_3p_window: true, terminal_mismatch: false, likelihood: 'likely', ortholog: false });
    should(st.pan_c.primary).match({ region: '2', start: 1001, end: 1280, size_delta: -20, ortholog: true });
    should(st.pan_c.others.map((o) => [o.region, o.start, o.size_delta, o.ortholog])).eql([['1', 201, 0, false]]);
    should(st.pan_c.other_amplicons).equal(1);
    should(st.pan_d).match({ ortholog_annotated: false, primary: null });
    should(st.pan_d.nearest).match({ likelihood: 'unlikely', right_mm_pos: [2, 1] });
    should(st.pan_e).match({ reason: 'NO_BLASTDB', primary: null });
    should(st.pan_x).match({ status: 'db_unavailable', reason: 'AMBIGUOUS_ASSEMBLY', display_name: 'pan_x' });
    should(st.pan_f.error.code).equal('BLAST_FAILED');
    should(st.pan_f.error.message).match(/pan_f\.dna\.toplevel\.nal/);
    should(st.pan_f.error.message.indexOf(world.root)).equal(-1);
    should(world.blastCalls('pan_f').length - callsF).equal(2);
    should(world.blastCalls('pan_g').length - callsG).equal(2);
    should(st.pan_g.primary).match({ ortholog: true });
    should(st.pan_h.primary).match({ approx: true, left_mm_pos: null, right_mm_pos: null });
    should(r.warnings).eql([{ code: 'NO_FASTA_FOR_REALIGN', message: 'genome sequence was not available for re-alignment of some sites; their mismatch counts are lower bounds (approx) (pan_h)' }]);

    // T: reference 3001 bp → maxSize 5002 finds the 4500 bp product that 4000 would miss
    should(pt).match({ id: 'T', reference_size: 3001, max_size: 5002 });
    should(pt.genomes.find((g) => g.system_name === 'pan_a')).match({ status: 'single_perfect', primary: { size: 4500, size_delta: 1499 } });
    should(pt.summary).eql({ genomes_total: 9, single_perfect: 1, single_mismatch: 0, multiple: 0, no_amplicon: 5, db_unavailable: 2, error: 1, amplifies: 1, truncated: 0 });

    // blastn per pan-genome genome: word size 6, one thread
    const panCall = world.calls.find((c) => c.args.indexOf(world.assemblies.pan_a.blastdb.dna) >= 0);
    should([panCall.args[panCall.args.indexOf('-word_size') + 1], panCall.args[panCall.args.indexOf('-num_threads') + 1]]).eql(['6', '1']);

    // partial documents: specificity first, then genomes appear with a consistent summary
    should(events.partials[0].specificity).be.ok();
    should(events.partials[0].pangenome).equal(null);
    should(events.partials).have.length(1 + GENOMES.length);
    let last = -1;
    for (const doc of events.partials.slice(1)) {
      for (const pair of doc.pangenome.pairs) {
        const s = pair.summary;
        should(s.genomes_total).equal(pair.genomes.length);
        should(pangenome.STATUSES.reduce((a, k) => a + s[k], 0)).equal(s.genomes_total);
      }
      should(doc.pangenome.pairs[0].genomes.length).be.aboveOrEqual(last);
      last = doc.pangenome.pairs[0].genomes.length;
    }
    should(last).equal(9);
    const finalProgress = events.progress[events.progress.length - 1];
    should(finalProgress).eql({ done: 10, total: 10, stage: 'pangenome', running: [] });
    should(events.progress.some((p) => p.stage === 'pangenome' && p.running.length > 1)).be.true();
    should(Object.keys(r.timings_ms)).containDeep(['reference', 'pan_a', 'pan_x', 'total']);
  });

  it('transcript mode: cDNA gene groups per genome and PANGENOME_TRANSCRIPT_MODELS_ONLY', async () => {
    const { ctx } = W.worldCtx(world);
    const r = await run(W.request({
      mode: 'transcript', gene_id: 'REFG1', transcript_id: 'REFG1.2', checks: ['pangenome'], genomes: ['pan_a', 'pan_b', 'pan_e'],
      pairs: [{ id: 'P', left: P.L, right: P.R }]
    }), ctx, { retryDelayMs: 0 });
    should(r.warnings.map((w) => w.code)).eql(['TRANSCRIPT_GENE_UNMAPPED', 'PANGENOME_TRANSCRIPT_MODELS_ONLY']);
    should(r.pangenome.target).equal('cdna');
    const pair = r.pangenome.pairs[0];
    should(pair).match({ reference_size: 300, max_size: 4000 });
    should(pair.summary).eql({ genomes_total: 3, single_perfect: 1, single_mismatch: 1, multiple: 0, no_amplicon: 0, db_unavailable: 1, error: 0, amplifies: 2, truncated: 0 });
    const [a, b, e] = pair.genomes;
    should(a).match({ system_name: 'pan_a', status: 'single_perfect', ortholog_annotated: true });
    should(a.primary).match({ gene_id: 'PANA_G1', orientation: 'LR', size: 300, size_delta: 0, ortholog: true });
    should(a.primary.isoforms.map((i) => i.transcript_id)).eql(['PANA_G1.1']);
    should(b.primary).match({ gene_id: 'PANB_G7', left_mm: 1, left_mm_pos: [2], mismatch_in_3p_window: true, ortholog: false });
    should(b.status).equal('single_mismatch');
    should(e).match({ system_name: 'pan_e', status: 'db_unavailable', reason: 'NO_BLASTDB' });
  });

  it('a timed-out pan-genome BLAST is not retried: the genome is error BLAST_TIMEOUT after one call (chk-retry-timeout)', async () => {
    const before = world.blastCalls('pan_t').length;
    const beforeA = world.blastCalls('pan_a').length;
    const r = await run(W.request({ mode: 'region', checks: ['pangenome'], genomes: ['pan_t', 'pan_a'], pairs: [{ id: 'P', left: P.L, right: P.R }] }), W.worldCtx(world, { gene: null }).ctx, { retryDelayMs: 0 });
    should(world.blastCalls('pan_t').length - before).equal(1);
    should(world.blastCalls('pan_a').length - beforeA).equal(1);
    const byName = {};
    for (const g of r.pangenome.pairs[0].genomes) byName[g.system_name] = g;
    should(byName.pan_t).match({ status: 'error', error: { code: 'BLAST_TIMEOUT' } });
    should(byName.pan_a.status).equal('single_perfect');
  });

  it('summary.truncated counts genomes whose search hit a cap (chk-pan-truncated-status)', async () => {
    const r = await run(W.request({ mode: 'region', checks: ['pangenome'], genomes: ['pan_a', 'pan_b', 'pan_c'], pairs: [{ id: 'P', left: P.L, right: P.R }] }),
      W.worldCtx(world, { gene: null, check: { max_realign_sites_per_genome: 1 } }).ctx, { retryDelayMs: 0 });
    const pair = r.pangenome.pairs[0];
    should(pair.genomes.map((g) => [g.system_name, g.truncated])).eql([['pan_a', true], ['pan_b', true], ['pan_c', true]]);
    should(pair.summary).match({ genomes_total: 3, truncated: 3 });
    should(r.specificity.pairs[0].truncated).be.true();
  });

  it('region mode: ortholog fields are null; abort during the pan-genome stage rejects', async () => {
    const { ctx } = W.worldCtx(world, { gene: null });
    const r = await run(W.request({ mode: 'region', checks: ['pangenome'], genomes: ['pan_a', 'pan_c'], pairs: [{ id: 'P', left: P.L, right: P.R }] }), ctx, { retryDelayMs: 0 });
    const [a, c] = r.pangenome.pairs[0].genomes;
    should(a).match({ ortholog_annotated: null, status: 'single_perfect' });
    should(a.primary.ortholog).equal(null);
    should(c.primary.ortholog).equal(null);
    should(r.pangenome.pairs[0].reference_size).equal(300);

    const reason = Object.assign(new Error('shutdown'), { code: 'WORKER_SHUTDOWN' });
    const aborting = W.worldCtx(world, { gene: null, onProgress: (p, controller) => { if (p.stage === 'pangenome') controller.abort(reason); } });
    const err = await run(W.request({ mode: 'region', checks: ['pangenome'], genomes: ['pan_a', 'pan_b', 'pan_c'], pairs: [{ id: 'P', left: P.L, right: P.R }] }), aborting.ctx, { retryDelayMs: 0 }).then(() => null, (x) => x);
    should(err).equal(reason);
  });
});
