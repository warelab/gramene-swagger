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
