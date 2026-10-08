'use strict';

// LONG TEST (opt-in): the allele caller over the full sorghum panel on real data (genotyping spec §7.5 "full panel", milestone
// M11). About 2.5 min wall and 1,000 CPU-s on squam (the check itself measured 143 s and 990 CPU-s, against an estimate of 1,357),
// with 8 blastn processes at nice 10; the rest of the real-data suite takes about 3 min.
//   PRIMERS_REALDATA=1 node --test --test-concurrency=1 --test-reporter=spec test/primers/unit/genotype_panel_realdata.test.js
// Read-only: mongo (maps, taxonomy, genes), /scratch/olson/fasta bgzip FASTA and BLAST DBs.
// rs871475760 (1:11109:C:A) with KASP set S1 of §2.9 (3 unique primers). The body is the design's check.request
// (fixtures/docs/capture-genotyping-design-rs871475760-kasp.json) cut to S1. check/normalize builds the job as POST /primers/check
// does (genomes omitted: every other genome of the species with a DNA BLAST DB, 119) and check/run.runWithStats runs it in-process
// with a worker-like ctx. The pinned values are those observed in M11 (HEAD 228a76b). The research megablast baseline over the same
// 119 assemblies is 73 C / 44 A / 2 missing (is36143, pi525695) with copies 1:79, 2:37, 4:1: the calls are the same genome by
// genome, and the copy counts differ only where the baseline's near-best bitscore filter drops a diverged copy.

const { describe, it, before, after } = require('node:test');
const should = require('should');

const REAL = process.env.PRIMERS_REALDATA === '1';
const PROCS = 8;

// Observed in M11 (identical in a separate run of S1 with S2 and an AS-PCR set, apart from copy envelopes and paralog counts).
// tx430nano's third copy (Scaffold_20:~32.56 Mb, gap-compressed identity 96.9 % over S1's window) is orthologous by §5.6 but under
// the baseline's 0.9 x best bitscore filter (95 % over 804 bp), so copies are 1:79, 2:36, 3:1, 4:1, not the §7.5 1:79, 2:37, 4:1.
// S1_REF's 293 off-target products over 120 genomes have 3 mismatches in a primer (WEAK_OFF_TARGETS).
const PINS = Object.freeze({
  estimate: { cpu_s: 1357, total: 120 },
  summary: { genomes_total: 119, ref: 73, alt: 44, other: 0, ambiguous: 0, missing: 2, unavailable: 0 },
  alt: ['sorghum_353', 'sorghum_grassl', 'sorghum_grif16309', 'sorghum_is12661', 'sorghum_ji2055t2tagi', 'sorghum_pi154844', 'sorghum_pi154987',
    'sorghum_pi180348', 'sorghum_pi276837', 'sorghum_pi510757', 'sorghum_pi513676', 'sorghum_pi532566', 'sorghum_pi533794', 'sorghum_pi533869',
    'sorghum_pi536008', 'sorghum_pi543243ph256', 'sorghum_pi554647ph355', 'sorghum_pi554649ph410', 'sorghum_pi554650ph449', 'sorghum_pi576434',
    'sorghum_pi595221excels235', 'sorghum_pi596332ph387lm', 'sorghum_pi597980', 'sorghum_pi601555ph308', 'sorghum_pi601556ph309',
    'sorghum_pi601557ph310', 'sorghum_pi601717ph333', 'sorghum_pi601719pha86', 'sorghum_pi601721phb86', 'sorghum_pi601743hp150',
    'sorghum_pi601744r159', 'sorghum_pi601756r145', 'sorghum_pi656015', 'sorghum_pi656027', 'sorghum_pi656029', 'sorghum_pi656050',
    'sorghum_pi656057', 'sorghum_pi656068', 'sorghum_pi656111', 'sorghum_pi660563', 'sorghum_rio', 'sorghum_riouncc', 'sorghum_tx2783pac',
    'sorghum_tx430nano'],
  missing: [['sorghum_is36143', 'no_orthologous_copy'], ['sorghum_pi525695', 'no_orthologous_copy']],
  copies: { 0: 2, 1: 79, 2: 36, 3: 1, 4: 1 },
  many_copies: [['sorghum_pi154987', 4], ['sorghum_tx430nano', 3]],
  paralogs: [],
  megablast: ['sorghum_is36143', 'sorghum_pi525695'],
  set_summary: { genomes_total: 119, predicted_ref: 73, predicted_alt: 44, both: 0, none: 2, no_call: 0, unknown: 0, weak: 0, agree: 119, disagree: 0, not_comparable: 0 },
  control: { status: 'pass', allele: 'ref', reasons: [] },
  not_agreeing: [],
  specificity: [['S1_REF', 'off_targets', 3], ['S1_ALT', 'specific', 0]],
  weak_count: [293]
});

describe('real data: rs871475760 S1 allele calls over the full sorghum panel (§7.5, long)', { skip: REAL ? false : 'set PRIMERS_REALDATA=1' }, () => {
  let H;
  let mongo;
  let cfg;
  let runWithStats;
  let checkNormalize;

  before(() => {
    H = require('../fixtures/verdicts/realdata_ctx');
    mongo = H.mongoConfig();
    cfg = H.primersConfig();
    runWithStats = require('../../../api/helpers/primers/check/run').runWithStats;
    checkNormalize = require('../../../api/helpers/primers/check/normalize');
  });

  after(() => {
    if (mongo) mongo.closeMongoDatabase();
  });

  it('KASP S1 over the reference and 119 other assemblies: allele summary, alt and missing genomes, copies, predictions, control', { timeout: 1800000 }, async (t) => {
    const capture = require('../fixtures/docs/capture-genotyping-design-rs871475760-kasp.json');
    const body = JSON.parse(JSON.stringify(capture.response.check.request));
    body.pairs = body.pairs.filter((p) => /^S1_/.test(p.id));
    body.genotyping.sets = body.genotyping.sets.filter((s) => s.id === 'S1');
    const norm = await checkNormalize.normalize(body, { cfg, mongo });
    should([norm.kind, norm.request.genomes.length, norm.estimate]).eql(['pangenome', 119, PINS.estimate]);

    const logs = [];
    const made = H.makeCtx({ cfg, resolved: norm.resolved, mongo, procs: PROCS, jobId: 'genotype_panel',
      log: { info() {}, warn: (m) => logs.push(m), error: (m) => logs.push(m) } });
    try {
      const t0 = Date.now();
      const { results, stats } = await runWithStats(norm.request, made.ctx);
      const g = results.genotyping;
      const genomes = g.genomes.filter((x) => !x.is_reference);
      const s = g.sets[0];
      t.diagnostic('run ' + (Date.now() - t0) + ' ms; summary ' + JSON.stringify(g.summary) + '; S1 ' + JSON.stringify(s.summary) +
        '; warnings ' + JSON.stringify(results.warnings.map((w) => [w.code, w.details ? w.details.count : null])) + (logs.length ? '; log ' + JSON.stringify(logs) : ''));

      should(results.warnings.filter((w) => /^(GENOTYPE_|REFERENCE_CONTROL)/.test(w.code))).eql([]);
      should(g.summary).eql(PINS.summary);
      should(genomes.filter((x) => x.allele === 'alt').map((x) => x.system_name).sort()).eql(PINS.alt);
      should(genomes.filter((x) => x.allele === 'missing').map((x) => [x.system_name, x.reason])).eql(PINS.missing);
      const copies = {};
      for (const x of genomes) copies[x.orthologous_copies] = (copies[x.orthologous_copies] || 0) + 1;
      should(copies).eql(PINS.copies);
      should(genomes.filter((x) => x.orthologous_copies > 2).map((x) => [x.system_name, x.orthologous_copies])).eql(PINS.many_copies);
      should(genomes.filter((x) => x.paralog_copies > 0).map((x) => [x.system_name, x.paralog_copies])).eql(PINS.paralogs);
      // every call comes from the amplicons; the megablast fallback runs only for the two genomes without an orthologous copy
      should(genomes.filter((x) => x.source === 'amplicon').length).equal(PINS.summary.ref + PINS.summary.alt);
      should(Object.keys(stats.genotype).filter((k) => stats.genotype[k].megablast).sort()).eql(PINS.megablast);
      // §7.5 [U]: pi536008 (scaffold2100 ends 56 bp before the variant) is called alt from its amplicons, without the fallback
      should(genomes.find((x) => x.system_name === 'sorghum_pi536008')).match({ allele: 'alt', observed: 'TAT', source: 'amplicon', orthologous_copies: 1 });

      should([s.summary, s.control, s.reference.predicted]).eql([PINS.set_summary, PINS.control, 'ref']);
      should(s.genomes.filter((c) => c.agrees !== true).map((c) => [c.system_name, c.predicted, c.reasons])).eql(PINS.not_agreeing);
      should(results.specificity.pairs.map((p) => [p.id, p.verdict, p.off_target_count])).eql(PINS.specificity);
      should(results.warnings.filter((w) => w.code === 'WEAK_OFF_TARGETS').map((w) => w.details.count)).eql(PINS.weak_count);
    } finally {
      made.cleanup();
    }
  });
});
