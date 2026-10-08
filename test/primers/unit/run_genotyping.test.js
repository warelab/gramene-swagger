'use strict';

// check/run.js and the genotyping extension (genotyping spec §5.5-§5.9), end to end on the verdicts fake world
// (fake BLAST reporting planted sites, real bgzip FASTA re-alignment, fake mongo).
// First half, "compat" (C5, §7.2; written before any genotyping code): a CheckRun whose request has no `genotyping`
// produces exactly today's Object.keys(results), in the final document and in every partial flush, across the stages
// the feature changes (reference; pan-genome genome, cDNA, db_unavailable and error entries; execute).
// The genotyping cases of §7.3 follow in M8.

const { describe, it, before, after } = require('node:test');
const should = require('should');

const { _internal: runInternal } = require('../../../api/helpers/primers/check/run');
const W = require('../fixtures/verdicts/fake_world');

// The §B.12 results keys, in order, as produced on primer-design d4af3a1 (verdicts.test.js asserts the same list).
const RESULTS_KEYS = ['engine', 'params', 'reference', 'sensitivity_note', 'primers', 'specificity', 'transcriptome', 'pangenome', 'warnings', 'timings_ms'];

describe('compat', () => {
  let world;
  let P;
  before(() => {
    world = W.standardWorld();
    P = world.primers;
  });
  after(() => world.cleanup());

  // A CheckRun over W.request(body), which never carries `genotyping` → { results, partials }.
  const checkRun = async (body) => {
    const request = W.request(body);
    should(request).not.have.property('genotyping');
    const { ctx, events } = W.worldCtx(world);
    const results = await new runInternal.CheckRun(request, ctx, { retryDelayMs: 0 }).execute();
    return { results, partials: events.partials };
  };

  const todaysKeys = (results, partials, flushes) => {
    should(Object.keys(results)).eql(RESULTS_KEYS);
    should(results).not.have.property('genotyping');
    should(partials).have.length(flushes);
    partials.forEach((doc) => should(Object.keys(doc)).eql(RESULTS_KEYS));
  };

  const ON_TARGET_P = { region: '1', start: 1001, end: 1300 };

  it('gene mode, specificity only: exactly today\'s Object.keys(results), no genotyping key', async () => {
    const { results, partials } = await checkRun({
      mode: 'gene',
      gene_id: 'REFG1',
      pairs: [
        { id: 'P', left: P.L, right: P.R, expected: ON_TARGET_P },
        { id: 'S', left: P.S_L, right: P.S_R, expected: { region: '8', start: 501, end: 700 } }
      ]
    });
    todaysKeys(results, partials, 1);
    should(results.specificity.pairs.map((p) => p.verdict)).eql(['off_targets', 'specific']);
  });

  it('gene mode with the pan-genome check (genome, db_unavailable and error entries): the same keys after every flush', async () => {
    const genomes = ['pan_a', 'pan_b', 'pan_e', 'pan_f'];
    const { results, partials } = await checkRun({ mode: 'gene', gene_id: 'REFG1', checks: ['pangenome'], genomes, pairs: [{ id: 'P', left: P.L, right: P.R, expected: ON_TARGET_P }] });
    todaysKeys(results, partials, 1 + genomes.length);
    should(results.pangenome.pairs[0].genomes.map((g) => [g.system_name, g.status])).eql([
      ['pan_a', 'single_perfect'], ['pan_b', 'single_mismatch'], ['pan_e', 'db_unavailable'], ['pan_f', 'error']
    ]);
  });

  it('transcript mode with the pan-genome check (cDNA entries): the same keys after every flush', async () => {
    const genomes = ['pan_a', 'pan_b'];
    const { results, partials } = await checkRun({ mode: 'transcript', gene_id: 'REFG1', transcript_id: 'REFG1.2', checks: ['pangenome'], genomes, pairs: [{ id: 'P', left: P.L, right: P.R }] });
    todaysKeys(results, partials, 2 + genomes.length);
    should(results.pangenome.target).equal('cdna');
    should(results.transcriptome).match({ target: 'cdna', gene_id: 'REFG1' });
  });
});

// ---- M8: the genotyping stages (§5.6-§5.9, §7.3) --------------------------------------------------------------------------------
// The genotyping world (fixtures/check_core/genotype/world.js): a reference whose chromosome 1 carries the recorded
// sorghum_bicolor 1:9000-15500 at its own coordinates, and pan-genome genomes built from the rs871475760 locus 1:10600-11600
// (1:11109 at 2510 of each genome's chromosome). Three sets over 1:11109:C:A: KASP S1 (reverse) and S2 (forward) of §2.9 and the
// gel AS-PCR set A1 of §2.10(d), with a -2 mismatch on both allele-specific primers. check/genotype.js runs for real.

describe('genotyping', () => {
  const G = require('../fixtures/check_core/genotype/world');
  const genotype = require('../../../api/helpers/primers/check/genotype');
  const blast = require('../../../api/helpers/primers/check/blast');
  const EXAMPLE = require('../fixtures/check_core/genotype/results_2_13.json');

  const GENOTYPING_KEYS = ['engine', 'params', 'reference', 'sensitivity_note', 'primers', 'specificity', 'transcriptome', 'pangenome', 'genotyping', 'warnings', 'timings_ms'];
  const MAIN = ['g01_ref', 'g02_alt', 'g03_other', 'g04_missing', 'g05_dup_alt', 'g06_paralog85', 'g07_alt_paralog', 'g08_ref_paralog', 'g09_common3p', 'g10_common2'];
  const UNKNOWN = { status: 'unknown', likelihood: null, mm_pos: null, residual_mm_pos: null };
  const ZERO_SET = { genomes_total: 0, predicted_ref: 0, predicted_alt: 0, both: 0, none: 0, no_call: 0, unknown: 0, weak: 0, agree: 0, disagree: 0, not_comparable: 0 };

  let gw;
  before(() => {
    gw = G.genotypingWorld();
  });
  after(() => gw.cleanup());

  // A CheckRun over a genotyping world → { results, partials, job, ctx }. opts: { world, procs, check, body, megablast }
  const genotypingRun = async (genomes, opts) => {
    const o = opts || {};
    const w = o.world || gw;
    const body = w.body(genomes, o.body);
    const request = W.request(body);
    request.genotyping = body.genotyping;
    const { ctx, events } = W.worldCtx(w.world, { gene: null, procs: o.procs || 4, check: o.check });
    if (o.megablast) ctx.spawnLines = G.withMegablast(w.world, o.megablast);
    const job = new runInternal.CheckRun(request, ctx, { retryDelayMs: 0 });
    const results = await job.execute();
    return { results, partials: events.partials, job, ctx };
  };
  const byName = (list) => {
    const m = {};
    for (const x of list) m[x.system_name] = x;
    return m;
  };
  const prediction = (c) => [c.system_name, c.predicted, c.strength, c.agrees, c.reasons, c.off_locus_products];
  const addsUp = (g) => {
    should(g.summary.ref + g.summary.alt + g.summary.other + g.summary.ambiguous + g.summary.missing + g.summary.unavailable).equal(g.summary.genomes_total);
    should(g.genomes.length).equal(1 + g.summary.genomes_total);
    for (const s of g.sets) {
      const m = s.summary;
      should(m.predicted_ref + m.predicted_alt + m.both + m.none + m.no_call + m.unknown).equal(m.genomes_total);
      should(m.agree + m.disagree + m.not_comparable).equal(m.genomes_total);
      should(m.genomes_total).equal(g.summary.genomes_total);
      should(s.genomes).have.length(m.genomes_total);
    }
  };

  describe('the §7.3 world: REF, ALT, third allele, no locus, two ALT copies, paralogs and private common-primer SNPs', () => {
    let run;
    before(async () => {
      run = await genotypingRun(MAIN);
    });

    it('results gain genotyping right after pangenome; engine algorithm 2, genotyping g1; no warnings; the §2.13 variant block', () => {
      should(Object.keys(run.results)).eql(GENOTYPING_KEYS);
      should(run.results.engine.algorithm_version).equal('2');
      const g = run.results.genotyping;
      should(Object.keys(g)).eql(['algorithm_version', 'variant', 'summary', 'genomes', 'sets']);
      should(g.algorithm_version).equal('g1');
      should(g.variant).eql(EXAMPLE.variant);
      should(run.results.warnings).eql([]);
      should(Object.keys(run.job.stats.genotype).sort()).eql(['ref'].concat(MAIN).sort());
    });

    it('calls: ref, alt, other, missing, alt with 2 copies, ref with the paralog excluded, ref, ref, ref, ref; the reference is its own control', () => {
      const g = run.results.genotyping;
      should(g.genomes.map((x) => [x.system_name, x.is_reference, x.allele, x.observed, x.source, x.orthologous_copies, x.paralog_copies, x.reason])).eql([
        ['ref', true, 'ref', 'TCT', 'amplicon', 1, 0, null],
        ['g01_ref', false, 'ref', 'TCT', 'amplicon', 1, 0, null],
        ['g02_alt', false, 'alt', 'TAT', 'amplicon', 1, 0, null],
        ['g03_other', false, 'other', 'TTT', 'amplicon', 1, 0, null],
        ['g04_missing', false, 'missing', null, null, 0, 0, 'no_orthologous_copy'],
        ['g05_dup_alt', false, 'alt', 'TAT', 'amplicon', 2, 0, null],
        ['g06_paralog85', false, 'ref', 'TCT', 'amplicon', 1, 1, null],
        ['g07_alt_paralog', false, 'ref', 'TCT', 'amplicon', 1, 1, null],
        ['g08_ref_paralog', false, 'ref', 'TCT', 'amplicon', 1, 1, null],
        ['g09_common3p', false, 'ref', 'TCT', 'amplicon', 1, 0, null],
        ['g10_common2', false, 'ref', 'TCT', 'amplicon', 1, 0, null]
      ]);
      // one copy from the three sets' products; its figures from A1's 353 bp segment, the longest alignment
      should(g.genomes[0].copies).eql([{ region: '1', start: 10880, end: 11172, strand: 1, variant_position: 11109, identity: 100, gap_compressed_identity: 100,
        aligned_length: 353, observed: 'TCT', flank_edits: 0, call: 'ref', anchors: 3, ortholog: null, source: 'amplicon' }]);
      const by = byName(g.genomes);
      should(by.g05_dup_alt.copies.map((c) => [c.region, c.start, c.end, c.strand, c.variant_position, c.call, c.anchors])).eql([
        ['1', 2281, 2573, 1, 2510, 'alt', 3], ['2', 2429, 2721, -1, 2492, 'alt', 3]
      ]);
      should(by.g03_other.copies[0]).match({ observed: 'TTT', call: 'other', identity: 99.72, flank_edits: 0 });
      should(g.summary).eql({ genomes_total: 10, ref: 6, alt: 2, other: 1, ambiguous: 0, missing: 1, unavailable: 0 });
      should(g.summary.genomes_total).equal(run.results.pangenome.pairs[0].summary.genomes_total);
      addsUp(g);
    });

    it('predictions, the same for S1, S2 and A1: ref, alt, none, none, alt, ref, both, ref, no_call, weak ref', () => {
      const sets = run.results.genotyping.sets;
      should(sets.map((s) => [s.id, s.ref_pair, s.alt_pair, s.orientation, s.deliberate_mismatch_positions])).eql([
        ['S1', 'S1_REF', 'S1_ALT', 'reverse', []], ['S2', 'S2_REF', 'S2_ALT', 'forward', []], ['A1', 'A1_REF', 'A1_ALT', 'reverse', [2, 2]]
      ]);
      for (const s of sets) {
        should(s.genomes.map(prediction)).eql([
          ['g01_ref', 'ref', 'normal', true, [], 0],
          ['g02_alt', 'alt', 'normal', true, [], 0],
          ['g03_other', 'none', null, true, [], 0],
          ['g04_missing', 'none', null, true, ['no_orthologous_copy'], 0],
          ['g05_dup_alt', 'alt', 'normal', true, [], 0],
          ['g06_paralog85', 'ref', 'normal', true, [], 0],
          ['g07_alt_paralog', 'both', 'normal', false, ['alt_signal_off_locus'], 1],
          ['g08_ref_paralog', 'ref', 'normal', true, ['ref_signal_off_locus'], 1],
          ['g09_common3p', 'no_call', null, null, ['common_primer_3p_mismatch'], 0],
          ['g10_common2', 'ref', 'weak', true, ['common_primer_weak'], 0]
        ], s.id);
        should(s.summary).eql({ genomes_total: 10, predicted_ref: 4, predicted_alt: 2, both: 1, none: 2, no_call: 1, unknown: 0, weak: 1, agree: 8, disagree: 1, not_comparable: 1 });
        should(prediction(s.reference)).eql(['ref', 'ref', 'normal', true, [], 0]);
        should(s.control).eql({ status: 'pass', allele: 'ref', reasons: [] });
      }
      should(sets[0].specificity).eql({ ref_pair: { verdict: 'specific', consistent_with_allele: true }, alt_pair: { verdict: 'specific', consistent_with_allele: true }, off_target_count: 0 });
      should(sets[2].specificity).eql({ ref_pair: { verdict: 'specific', consistent_with_allele: true }, alt_pair: { verdict: 'on_target_missing', consistent_with_allele: true }, off_target_count: 0 });
    });

    it('primer calls: KASP reads a terminal mismatch; the common primer is read from its own site; AS-PCR never turns its deliberate mismatch into no_call', () => {
      const [S1, , A1] = run.results.genotyping.sets;
      const s1 = byName(S1.genomes);
      should(S1.reference).match({
        ref_primer: { status: 'match', likelihood: 'likely', mm_pos: [], residual_mm_pos: [] },
        alt_primer: { status: 'terminal_mismatch', likelihood: 'likely_weak', mm_pos: [1], residual_mm_pos: [1] },
        common_primer: { status: 'match', likelihood: 'likely', mm_pos: [], residual_mm_pos: [] }
      });
      should(s1.g02_alt).match({
        ref_primer: { status: 'terminal_mismatch', likelihood: 'likely_weak', mm_pos: [1], residual_mm_pos: [1] },
        alt_primer: { status: 'match', likelihood: 'likely', mm_pos: [], residual_mm_pos: [] },
        common_primer: { status: 'match', likelihood: 'likely', mm_pos: [], residual_mm_pos: [] }
      });
      should(s1.g09_common3p.common_primer).eql({ status: 'terminal_mismatch', likelihood: 'likely_weak', mm_pos: [1], residual_mm_pos: [1] });
      should(s1.g10_common2.common_primer).eql({ status: 'weak', likelihood: 'likely', mm_pos: [2], residual_mm_pos: [2] });
      // A1 on the reference: the ALT pair's product is unlikely ([2,1]) yet its common site, like the REF pair's, is a perfect match
      should(A1.reference).match({
        ref_primer: { status: 'match', likelihood: 'likely', mm_pos: [2], residual_mm_pos: [] },
        alt_primer: { status: 'blocked', likelihood: 'unlikely', mm_pos: [2, 1], residual_mm_pos: [1] },
        common_primer: { status: 'match', likelihood: 'likely', mm_pos: [], residual_mm_pos: [] },
        predicted: 'ref'
      });
      should(byName(A1.genomes).g02_alt).match({
        ref_primer: { status: 'blocked', likelihood: 'unlikely', mm_pos: [2, 1], residual_mm_pos: [1] },
        alt_primer: { status: 'match', likelihood: 'likely', mm_pos: [2], residual_mm_pos: [] },
        common_primer: { status: 'match', likelihood: 'likely' },
        predicted: 'alt',
        agrees: true
      });
      // the pan-genome block alone cannot tell the alleles apart: both S1 pairs amplify on the ALT genome
      const pan = byName(run.results.pangenome.pairs.find((p) => p.id === 'S1_REF').genomes);
      should(pan.g02_alt.status).equal('single_mismatch');
    });

    it('partial flushes: the first carries only the reference; each later one adds one genome, with summaries that add up', () => {
      const parts = run.partials;
      should(parts).have.length(1 + MAIN.length);
      parts.forEach((doc, k) => {
        should(Object.keys(doc)).eql(GENOTYPING_KEYS);
        should(doc.genotyping.summary.genomes_total).equal(k);
        should(doc.genotyping.genomes[0]).match({ system_name: 'ref', is_reference: true, allele: 'ref' });
        should(doc.pangenome ? doc.pangenome.pairs[0].genomes.length : 0).equal(k);
        addsUp(doc.genotyping);
      });
      should(parts[0].genotyping.sets.map((s) => [s.summary, s.genomes, s.control.status])).eql([[ZERO_SET, [], 'pass'], [ZERO_SET, [], 'pass'], [ZERO_SET, [], 'pass']]);
    });
  });

  it('a specificity-only job: the reference entry alone, every set predicts ref on it, zero summaries', async () => {
    const { results, partials } = await genotypingRun([], { body: { checks: [] } });
    should(Object.keys(results)).eql(GENOTYPING_KEYS);
    const g = results.genotyping;
    should(g.genomes.map((x) => [x.system_name, x.allele])).eql([['ref', 'ref']]);
    should(g.summary).eql({ genomes_total: 0, ref: 0, alt: 0, other: 0, ambiguous: 0, missing: 0, unavailable: 0 });
    for (const s of g.sets) {
      should([s.reference.predicted, s.control.status, s.summary, s.genomes]).eql(['ref', 'pass', ZERO_SET, []]);
    }
    should(partials).have.length(1);
    should(partials[0].genotyping).eql(g);
  });

  it('megablast fallback: a genome without products is called from its HSP; a failed megablast is fallback_failed; no DB or failed BLAST is unavailable', async () => {
    const { results, ctx } = await genotypingRun(['g04_missing', 'm1_mb_alt', 'm2_mb_fail', 'm3_nodb', 'm4_blast_error'], { procs: 1, megablast: { fail: (sys) => sys === 'm2_mb_fail' } });
    const g = byName(results.genotyping.genomes);
    should(g.m1_mb_alt).eql({
      system_name: 'm1_mb_alt', display_name: 'M1_MB_ALT', is_reference: false, allele: 'alt', observed: 'TAT', source: 'megablast',
      copies: [{ region: '1', start: 2295, end: 2725, strand: 1, variant_position: 2510, identity: 96.98, gap_compressed_identity: 96.98, aligned_length: 431,
        observed: 'TAT', flank_edits: 12, call: 'alt', anchors: 0, ortholog: null, source: 'megablast' }],
      orthologous_copies: 1, paralog_copies: 0, reason: null
    });
    should(g.m2_mb_fail).match({ allele: 'missing', observed: null, source: null, copies: [], reason: 'fallback_failed' });
    should(g.g04_missing).match({ allele: 'missing', reason: 'no_orthologous_copy' });
    should(g.m3_nodb).eql({ system_name: 'm3_nodb', display_name: 'M3_NODB', is_reference: false, allele: 'unavailable', observed: null, source: null, copies: [],
      orthologous_copies: 0, paralog_copies: 0, reason: 'db_unavailable' });
    should(g.m4_blast_error).match({ allele: 'unavailable', reason: 'blast_error' });
    should(results.genotyping.summary).eql({ genomes_total: 5, ref: 0, alt: 1, other: 0, ambiguous: 0, missing: 2, unavailable: 2 });
    addsUp(results.genotyping);
    for (const s of results.genotyping.sets) {
      const c = byName(s.genomes);
      should(prediction(c.m1_mb_alt)).eql(['m1_mb_alt', 'none', null, false, [], 0]);
      for (const sys of ['m3_nodb', 'm4_blast_error']) {
        should(c[sys]).eql({ system_name: sys, ref_primer: UNKNOWN, alt_primer: UNKNOWN, common_primer: UNKNOWN, predicted: 'unknown', strength: null, agrees: null, reasons: [], off_locus_products: 0 });
      }
    }
    should(results.warnings).eql([{ code: 'GENOTYPE_FALLBACK_FAILED', message: 'the megablast fallback of the allele caller failed; these genomes are missing (m2_mb_fail)' }]);
    // megablast ran only for the genomes without an orthologous copy, with the §5.6 arguments, in the job's private directory
    const calls = ctx.spawnLines.calls;
    should(calls.map((c) => c.system_name)).eql(['g04_missing', 'm1_mb_alt', 'm2_mb_fail']);
    for (const c of calls) {
      should(c.args).eql(blast.buildMegablastArgs({ db: gw.world.assemblies[c.system_name].blastdb.dna, threads: 1 }));
      should(c.cwd).equal(ctx.tmpdir());
    }
  });

  it('the megablast budget: past genotype_max_megablast a genome is missing with fallback_budget and warning GENOTYPE_FALLBACK_BUDGET', async () => {
    const { results, ctx } = await genotypingRun(['g04_missing', 'm1_mb_alt', 'm2_mb_fail'], { procs: 1, check: { genotype_max_megablast: 1 }, megablast: {} });
    const g = byName(results.genotyping.genomes);
    should([g.g04_missing.reason, g.m1_mb_alt.reason, g.m2_mb_fail.reason]).eql(['no_orthologous_copy', 'fallback_budget', 'fallback_budget']);
    should(g.m1_mb_alt).match({ allele: 'missing', source: null, copies: [] });
    should(ctx.spawnLines.calls.map((c) => c.system_name)).eql(['g04_missing']);
    should(results.warnings).eql([{ code: 'GENOTYPE_FALLBACK_BUDGET', message: 'the per-job megablast fallback budget (1) was reached; these genomes are missing (m1_mb_alt, m2_mb_fail)' }]);
  });

  it('an exception inside the allele caller makes that genome unavailable (call_failed) with GENOTYPE_FAILED; the job completes', async () => {
    const original = genotype.callGenome;
    genotype.callGenome = async (input, deps) => {
      if (input.system_name === 'g02_alt') throw new TypeError('caller bug');
      return original(input, deps);
    };
    let run;
    try {
      run = await genotypingRun(['g01_ref', 'g02_alt'], { procs: 1 });
    } finally {
      genotype.callGenome = original;
    }
    const g = byName(run.results.genotyping.genomes);
    should([g.ref.allele, g.g01_ref.allele, g.g02_alt.allele, g.g02_alt.reason]).eql(['ref', 'ref', 'unavailable', 'call_failed']);
    for (const s of run.results.genotyping.sets) should(prediction(byName(s.genomes).g02_alt)).eql(['g02_alt', 'unknown', null, null, [], 0]);
    should(run.results.warnings).eql([{ code: 'GENOTYPE_FAILED', message: 'the allele caller failed; these genomes are unavailable (g02_alt)' }]);
    should(run.results.pangenome.pairs[0].genomes.map((x) => x.status)).eql(['single_perfect', 'single_mismatch']);
  });

  it('reference control: a reference that also carries the ALT locus is ambiguous, every set fails, REFERENCE_CONTROL_FAILED {allele, sets}', async () => {
    const other = G.genotypingWorld({ referenceAltCopy: true });
    try {
      const { results } = await genotypingRun([], { world: other, body: { checks: [] } });
      const ref = results.genotyping.genomes[0];
      should(ref).match({ system_name: 'ref', is_reference: true, allele: 'ambiguous', observed: null, source: 'amplicon', orthologous_copies: 2, paralog_copies: 0 });
      should(ref.copies.map((c) => [c.region, c.call])).eql([['1', 'ref'], ['2', 'alt']]);
      for (const s of results.genotyping.sets) {
        should(s.reference).match({ predicted: 'both', agrees: null });
        should(s.control).eql({ status: 'fail', allele: 'ambiguous', reasons: ['allele_not_ref', 'prediction_not_ref'] });
      }
      should(results.warnings).eql([{
        code: 'REFERENCE_CONTROL_FAILED',
        message: 'the reference control failed for sets S1, S2, A1: the reference ref must be called ref (it is ambiguous) and every set must predict ref on it',
        details: { allele: 'ambiguous', sets: ['S1', 'S2', 'A1'] }
      }]);
    } finally {
      other.cleanup();
    }
  });
});

// ---- M8b: the off-locus mismatch threshold (check.genotype_offlocus_max_mismatches) -------------------------------------------------
// An off-locus product changes a prediction only when each of its primers has at most genotype_offlocus_max_mismatches (default 2)
// mismatches; weaker ones, which the specificity and pan-genome blocks still list, go into results warning WEAK_OFF_TARGETS.
// g11 and g12 carry the REF locus plus an ALT paralog whose three common sites have 3 (g11) or 2 (g12) mismatches, all at least
// 11 nt from the 3' end; g07's ALT paralog is perfect.

describe('genotyping: the off-locus mismatch threshold (M8b)', () => {
  const G = require('../fixtures/check_core/genotype/world');
  const GENOMES = ['g07_alt_paralog', 'g11_alt_paralog_3mm', 'g12_alt_paralog_2mm'];
  let gw;
  before(() => {
    gw = G.genotypingWorld();
  });
  after(() => gw.cleanup());

  // A CheckRun over the three genomes with ctx.config.check overrides → a JSON copy of the results.
  const thresholdRun = async (check) => {
    const body = gw.body(GENOMES);
    const request = W.request(body);
    request.genotyping = body.genotyping;
    const { ctx } = W.worldCtx(gw.world, { gene: null, procs: 4, check });
    return JSON.parse(JSON.stringify(await new runInternal.CheckRun(request, ctx, { retryDelayMs: 0 }).execute()));
  };
  const predictions = (results) => results.genotyping.sets.map((s) => [s.id, s.control.status].concat(s.genomes.map((x) => [x.system_name, x.predicted, x.agrees, x.reasons, x.off_locus_products])));
  const everySet = (rows) => ['S1', 'S2', 'A1'].map((id) => [id, 'pass'].concat(rows));
  const both = (sys) => [sys, 'both', false, ['alt_signal_off_locus'], 1];
  const ref = (sys) => [sys, 'ref', true, [], 0];
  const example = (sys, set, start, end, size, left, right) => ({ system_name: sys, set_id: set, pair_id: set + '_ALT', allele: 'alt', region: '2', start, end, size, orientation: 'LR', left_mm: left, right_mm: right });
  const G11 = [example('g11_alt_paralog_3mm', 'S1', 2469, 2533, 65, 3, 0), example('g11_alt_paralog_3mm', 'S2', 2482, 2573, 92, 0, 3), example('g11_alt_paralog_3mm', 'A1', 2281, 2533, 253, 3, 0)];

  it('default 2: the ALT paralog with 3 common-site mismatches no longer makes both and is reported; with 2 mismatches it still does', async () => {
    const results = await thresholdRun();
    should(predictions(results)).eql(everySet([both('g07_alt_paralog'), ref('g11_alt_paralog_3mm'), both('g12_alt_paralog_2mm')]));
    should(results.genotyping.genomes.map((x) => [x.system_name, x.allele, x.paralog_copies])).eql([['ref', 'ref', 0], ['g07_alt_paralog', 'ref', 1], ['g11_alt_paralog_3mm', 'ref', 1], ['g12_alt_paralog_2mm', 'ref', 1]]);
    should(results.warnings).eql([{
      code: 'WEAK_OFF_TARGETS',
      message: 'off-target products with more than 2 mismatches in a primer do not change the allele predictions; the specificity and pan-genome results still list them (g11_alt_paralog_3mm)',
      details: { count: 3, max_mismatches: 2, examples: G11 }
    }]);
    // the pan-genome block still counts the weak paralog product as amplifying
    const pan = results.pangenome.pairs.find((p) => p.id === 'S1_ALT').genomes.find((x) => x.system_name === 'g11_alt_paralog_3mm');
    should(pan.status).equal('multiple');
  });

  it('the threshold follows check.genotype_offlocus_max_mismatches: 3 counts the 3-mismatch paralog, 1 drops the 2-mismatch one, 0 keeps a perfect one', async () => {
    const three = await thresholdRun({ genotype_offlocus_max_mismatches: 3 });
    should(predictions(three)).eql(everySet([both('g07_alt_paralog'), both('g11_alt_paralog_3mm'), both('g12_alt_paralog_2mm')]));
    should(three.warnings).eql([]);
    const one = await thresholdRun({ genotype_offlocus_max_mismatches: 1 });
    should(predictions(one)).eql(everySet([both('g07_alt_paralog'), ref('g11_alt_paralog_3mm'), ref('g12_alt_paralog_2mm')]));
    should(one.warnings.map((w) => [w.code, w.details.count, w.details.max_mismatches])).eql([['WEAK_OFF_TARGETS', 6, 1]]);
    should(one.warnings[0].message).endWith('(g11_alt_paralog_3mm, g12_alt_paralog_2mm)');
    should(one.warnings[0].details.examples).eql(G11.concat([example('g12_alt_paralog_2mm', 'S1', 2469, 2533, 65, 2, 0), example('g12_alt_paralog_2mm', 'S2', 2482, 2573, 92, 0, 2)]));
    const zero = await thresholdRun({ genotype_offlocus_max_mismatches: 0 });
    should(predictions(zero).map((row) => row[2])).eql([both('g07_alt_paralog'), both('g07_alt_paralog'), both('g07_alt_paralog')]);
  });
});
