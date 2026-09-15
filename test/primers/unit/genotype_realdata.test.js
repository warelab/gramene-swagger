'use strict';

// Real-data allele calls through check/run.js (genotyping spec §7.5, milestone M8; opt-in):
//   PRIMERS_REALDATA=1 node --test --test-concurrency=1 --test-reporter=spec test/primers/unit/genotype_realdata.test.js
// Read-only: mongo (maps, taxonomy), /scratch/olson/fasta bgzip FASTA and BLAST DBs; blastn at nice 10 with ctx.procs 8.
// rs871475760 (1:11109:C:A) with KASP S1 and S2 of §2.9 and an AS-PCR set on S1's common primer (a -2 mismatch on both
// allele-specific primers), over the 11 research assemblies: the §7.5 allele and copy table, the §2.13 copies and primer calls,
// the reference control and the predictions, plus the measured cost of the genotype stage (each genome's call replayed alone)
// and of one megablast fallback per assembly (whose calls must match). Then pi536008, whose amplicon-window identity is 78.2 %:
// its amplicon call, and what the megablast fallback alone gives.
// Unlike §7.5, which was verified on megablast locus alignments only, S1 and A1 predict both on the 5 ALT genomes: their REF
// pairs have 3 amplifying genome-wide off-targets (see the predictions block).

const { describe, it, before, after } = require('node:test');
const should = require('should');

const REAL = process.env.PRIMERS_REALDATA === '1';
const PROCS = 8;

const COMMON = 'AGCTTCTCTAAGTGGTTATCCGA';
const PAIRS = Object.freeze([
  { id: 'S1_REF', left: COMMON, right: 'ATCTTTGACTAGCGAGAAATTCAG', expected: { region: '1', start: 11068, end: 11132 } },
  { id: 'S1_ALT', left: COMMON, right: 'ATCTTTGACTAGCGAGAAATTCAT', expected: { region: '1', start: 11068, end: 11132 } },
  { id: 'S2_REF', left: 'GGTTATCCGAATATAGTCATACTCTATTC', right: 'TCTTTGTCTACTGAGAAATCCAGA', expected: { region: '1', start: 11081, end: 11172 } },
  { id: 'S2_ALT', left: 'GGTTATCCGAATATAGTCATACTCTATTA', right: 'TCTTTGTCTACTGAGAAATCCAGA', expected: { region: '1', start: 11081, end: 11172 } },
  { id: 'A1_REF', left: COMMON, right: 'ATCTTTGACTAGCGAGAAATTCGG', expected: { region: '1', start: 11068, end: 11132 } },
  { id: 'A1_ALT', left: COMMON, right: 'ATCTTTGACTAGCGAGAAATTCGT', expected: { region: '1', start: 11068, end: 11132 } }
]);
const GENOTYPING = Object.freeze({
  variant: { region: '1', position: 11109, ref: 'C', alt: 'A' },
  sets: [{ id: 'S1', ref_pair: 'S1_REF', alt_pair: 'S1_ALT' }, { id: 'S2', ref_pair: 'S2_REF', alt_pair: 'S2_ALT' }, { id: 'A1', ref_pair: 'A1_REF', alt_pair: 'A1_ALT' }]
});
// §7.5 / research §3: allele and orthologous copies
const TABLE = Object.freeze({
  sorghum_bicolorv5: ['ref', 1], sorghum_austrcf317961: ['ref', 1], sorghum_pi565121: ['ref', 1], sorghum_pi655972: ['ref', 1],
  sorghum_pi180348: ['alt', 2], sorghum_pi276837: ['alt', 2], sorghum_pi656027: ['alt', 2], sorghum_pi510757: ['alt', 1],
  sorghum_rio: ['alt', 2], sorghum_pi329250: ['ref', 1], sorghum_s3691: ['ref', 1]
});

describe('real data: rs871475760 allele calls in the research assemblies (§7.5)', { skip: REAL ? false : 'set PRIMERS_REALDATA=1' }, () => {
  let H;
  let mongo;
  let cfg;
  let runInternal;
  let runWithStats;
  let genotype;
  let EXAMPLE;

  before(() => {
    H = require('../fixtures/verdicts/realdata_ctx');
    mongo = H.mongoConfig();
    cfg = H.primersConfig();
    const run = require('../../../api/helpers/primers/check/run');
    runInternal = run._internal;
    runWithStats = run.runWithStats;
    genotype = require('../../../api/helpers/primers/check/genotype');
    EXAMPLE = require('../fixtures/check_core/genotype/results_2_13.json');
  });

  after(() => {
    if (mongo) mongo.closeMongoDatabase();
  });

  const byName = (list) => {
    const m = {};
    for (const x of list) m[x.system_name] = x;
    return m;
  };

  // One genotyping check through run.js. Every genotypeGenome call is captured, with its job, for the replays below.
  async function check(name, genomes, t) {
    const request = H.normalizedRequest({ system_name: 'sorghum_bicolor', mode: 'region', checks: ['pangenome'], genomes, pairs: PAIRS });
    request.genotyping = JSON.parse(JSON.stringify(GENOTYPING));
    const resolved = await H.buildResolved({ cfg, mongo, systemName: 'sorghum_bicolor', genomes: request.genomes });
    const logs = [];
    const made = H.makeCtx({ cfg, resolved, mongo, procs: PROCS, jobId: name, log: { info() {}, warn: (m) => logs.push(m), error: (m) => logs.push(m) } });
    const proto = runInternal.CheckRun.prototype;
    const original = proto.genotypeGenome;
    const captured = [];
    proto.genotypeGenome = function (asm, perPair, info) {
      captured.push({ job: this, asm, perPair, info });
      return original.call(this, asm, perPair, info);
    };
    const t0 = Date.now();
    let out;
    try {
      out = await runWithStats(request, made.ctx);
    } finally {
      proto.genotypeGenome = original;
    }
    t.diagnostic(name + ': run ' + (Date.now() - t0) + ' ms; timings_ms ' + JSON.stringify(out.results.timings_ms) + (logs.length ? '; log ' + JSON.stringify(logs) : ''));
    t.diagnostic(name + ': in-run genotype stage ' + JSON.stringify(out.stats.genotype));
    const job = captured[0].job;
    // The finished job's internal signal is aborted; a replay needs a live one.
    job.internal = new AbortController();
    job.signal = AbortSignal.any([job.external, job.internal.signal]);
    return { results: out.results, stats: out.stats, captured, job, original, made };
  }

  // One megablast fallback against asm, timed by /usr/bin/time (user + sys of the blastn child).
  async function timedMegablast(job, asm, sys) {
    const spawn = job.spawnLines;
    let cpu = null;
    job.spawnLines = async (cmd, args, o) => {
      const res = await spawn('/usr/bin/time', ['-f', 'GENOTYPE_TIME %e %U %S', cmd].concat(args), o);
      const m = /GENOTYPE_TIME (\S+) (\S+) (\S+)/.exec(res.stderrTail || '');
      if (m) cpu = { wall_s: Number(m[1]), cpu_s: Math.round((Number(m[2]) + Number(m[3])) * 100) / 100 };
      return res;
    };
    try {
      const mb = await job.genotypeMegablast(asm, sys);
      return { mb, cpu };
    } finally {
      job.spawnLines = spawn;
    }
  }

  it('KASP S1 and S2 and an AS-PCR set on S1\'s common primer over the 11 research assemblies: the §7.5 table, the §2.13 copies and calls, the predictions and their cost', { timeout: 900000 }, async (t) => {
    const c = await check('rs871475760_eleven', Object.keys(TABLE), t);
    try {
      const g = c.results.genotyping;
      t.diagnostic('calls ' + JSON.stringify(g.genomes.map((x) => [x.system_name, x.allele, x.observed, x.source, x.orthologous_copies, x.paralog_copies,
        x.copies.map((k) => [k.region, k.start, k.end, k.strand, k.variant_position, k.identity, k.gap_compressed_identity, k.aligned_length, k.flank_edits, k.call, k.anchors])])));
      t.diagnostic('predictions ' + JSON.stringify(g.sets.map((s) => [s.id, s.control.status, s.summary, s.genomes.filter((x) => x.agrees !== true).map((x) => [x.system_name, x.predicted, x.reasons])])));
      should(c.results.warnings.filter((w) => /^(GENOTYPE_|REFERENCE_CONTROL)/.test(w.code))).eql([]);
      const by = byName(g.genomes);
      t.diagnostic('paralog_copies ' + JSON.stringify(g.genomes.map((x) => [x.system_name, x.paralog_copies])));
      for (const sys of Object.keys(TABLE)) {
        should([by[sys].allele, by[sys].orthologous_copies, by[sys].source, by[sys].observed]).eql([TABLE[sys][0], TABLE[sys][1], 'amplicon', TABLE[sys][0] === 'ref' ? 'TCT' : 'TAT'], sys);
      }
      should(g.summary).eql({ genomes_total: 11, ref: 6, alt: 5, other: 0, ambiguous: 0, missing: 0, unavailable: 0 });
      // §2.13 (sets S1 and S2): the reference, bicolorv5, pi180348 and pi329250 entries. The AS-PCR set's product adds one anchor per
      // copy, and pi180348 has one paralog of the locus, at 1:72.9 Mb (S2's product there aligns at 90 %), which the megablast loci
      // behind the example did not include.
      const PARALOGS = { sorghum_bicolor: 0, sorghum_bicolorv5: 0, sorghum_pi180348: 1, sorghum_pi329250: 0 };
      for (const e of EXAMPLE.genomes) {
        should(by[e.system_name]).eql(Object.assign({}, e, { paralog_copies: PARALOGS[e.system_name], copies: e.copies.map((k) => Object.assign({}, k, { anchors: k.anchors + 1 })) }), e.system_name);
      }
      const primerCalls = (x) => [x.system_name, x.ref_primer, x.alt_primer, x.common_primer];
      for (const s of EXAMPLE.sets) {
        const got = g.sets.find((x) => x.id === s.id);
        should(primerCalls(got.reference)).eql(primerCalls(s.reference), s.id);
        const mine = byName(got.genomes);
        for (const e of s.genomes) should(primerCalls(mine[e.system_name])).eql(primerCalls(e), s.id + ' ' + e.system_name);
      }
      // §7.5: S1's own anchor on bicolorv5 is 1:31523-31587
      should(c.results.pangenome.pairs.find((p) => p.id === 'S1_REF').genomes.find((x) => x.system_name === 'sorghum_bicolorv5').primary).match({ region: '1', start: 31523, end: 31587 });
      should(by.sorghum_pi180348.copies.map((k) => [k.start, k.end, k.strand])).eql([[15028, 15132, -1], [38342, 38446, 1]]);
      should(byName(g.sets[2].genomes).sorghum_pi180348).match({ ref_primer: { status: 'blocked', mm_pos: [2, 1] }, alt_primer: { status: 'match', mm_pos: [2] }, common_primer: { status: 'match' } });

      // Predictions (§5.7). On the locus every set reads the genome's allele. Off the locus, S1_REF and A1_REF (both on S1's common
      // primer) amplify 3 genome-wide off-targets each (likely: 3 mismatches per primer, none at a 3' base), which the reference
      // already reports as off_targets; they add FAM on the 5 ALT genomes (both, ref_signal_off_locus) and warn in the control.
      // S2 is specific; on pi180348 its REF pair amplifies the 1:72.9 Mb paralog, which carries C.
      const spec = byName(c.results.specificity.pairs.map((p) => Object.assign({ system_name: p.id }, p)));
      should(['S1_REF', 'S1_ALT', 'S2_REF', 'S2_ALT', 'A1_REF', 'A1_ALT'].map((id) => [spec[id].verdict, spec[id].off_target_count])).eql([
        ['off_targets', 3], ['specific', 0], ['specific', 0], ['specific', 0], ['off_targets', 3], ['on_target_missing', 0]
      ]);
      for (const s of g.sets) {
        const rows = byName(s.genomes);
        for (const sys of Object.keys(TABLE)) {
          const allele = TABLE[sys][0];
          const offLocus = allele === 'alt' && (s.id !== 'S2' || sys === 'sorghum_pi180348');
          const want = offLocus ? ['both', false] : [allele, true];
          should([rows[sys].predicted, rows[sys].agrees]).eql(want, s.id + ' ' + sys);
          if (offLocus) should(rows[sys].reasons).eql(['ref_signal_off_locus'], s.id + ' ' + sys);
        }
        should(s.reference.predicted).equal('ref', s.id);
        should(s.control).eql(s.id === 'S2' ? { status: 'pass', allele: 'ref', reasons: [] } : { status: 'warn', allele: 'ref', reasons: ['off_locus_products', 'ref_signal_off_locus'] }, s.id);
      }
      should(g.sets.map((s) => [s.id, s.summary.agree, s.summary.both])).eql([['S1', 6, 5], ['S2', 10, 1], ['A1', 6, 5]]);

      // Cost of the genotype stage: each genome's call replayed alone (one warm-up, then the mean of 3), and one megablast each.
      const cost = [];
      for (const k of c.captured) {
        await c.original.call(c.job, k.asm, k.perPair, k.info);
        const u0 = process.cpuUsage();
        const w0 = process.hrtime.bigint();
        for (let r = 0; r < 3; r++) await c.original.call(c.job, k.asm, k.perPair, k.info);
        const u = process.cpuUsage(u0);
        const { mb, cpu } = await timedMegablast(c.job, k.asm, k.info.system_name);
        should(mb.status).equal('ok', k.info.system_name);
        const res = genotype.megablastCopies(mb.rows, mb.query, c.job.geno, c.job.ccfg);
        const calls = res.copies.map((x) => x.copy.call);
        // the fallback alone reads the same allele from the same number of copies
        should(calls).eql(new Array(by[k.info.system_name].orthologous_copies).fill(by[k.info.system_name].allele), k.info.system_name);
        cost.push({
          system_name: k.info.system_name,
          caller_cpu_s: Math.round((u.user + u.system) / 3e3) / 1e3,
          caller_wall_ms: Math.round(Number(process.hrtime.bigint() - w0) / 3e6),
          megablast: cpu,
          megablast_copies: calls.length
        });
      }
      t.diagnostic('genotype stage cost per genome ' + JSON.stringify(cost));
      const mean = (xs) => Math.round(1000 * xs.reduce((a, b) => a + b, 0) / xs.length) / 1000;
      t.diagnostic('mean caller CPU-s ' + mean(cost.map((x) => x.caller_cpu_s)) + '; mean megablast CPU-s ' + mean(cost.filter((x) => x.megablast).map((x) => x.megablast.cpu_s)));
      for (const x of cost) should(x.caller_cpu_s).be.below(1, x.system_name);
    } finally {
      c.made.cleanup();
    }
  });

  it('pi536008 (amplicon-window identity 78.2 %, scaffold2100 ends 176 bp past the variant): called alt from its amplicons; the megablast fallback alone keeps no HSP', { timeout: 900000 }, async (t) => {
    const c = await check('rs871475760_pi536008', ['sorghum_pi536008'], t);
    try {
      const e = byName(c.results.genotyping.genomes).sorghum_pi536008;
      t.diagnostic('pi536008 ' + JSON.stringify(e));
      t.diagnostic('pi536008 predictions ' + JSON.stringify(c.results.genotyping.sets.map((s) => [s.id, s.genomes[0].predicted, s.genomes[0].agrees, s.genomes[0].reasons])));
      should(e).match({ allele: 'alt', observed: 'TAT', source: 'amplicon', orthologous_copies: 1, reason: null });
      should(e.copies[0]).match({ region: 'scaffold2100', strand: -1, call: 'alt' });
      // the longest anchor runs off the scaffold end: its raw identity is low, its gap-compressed identity is not
      should(e.copies[0].identity).be.below(95);
      should(e.copies[0].gap_compressed_identity).be.aboveOrEqual(95);
      should(c.stats.genotype.sorghum_pi536008.megablast).equal(null);
      const k = c.captured.find((x) => x.info.system_name === 'sorghum_pi536008');
      const { mb, cpu } = await timedMegablast(c.job, k.asm, 'sorghum_pi536008');
      const res = genotype.megablastCopies(mb.rows, mb.query, c.job.geno, c.job.ccfg);
      t.diagnostic('pi536008 megablast fallback alone ' + JSON.stringify({ status: mb.status, cpu, hsps: res.hsps, kept: res.kept, covering: res.covering,
        rows: mb.rows.map((r) => [r.sseqid, r.sstart, r.send, r.pident, r.length, r.bitscore, r.qstart, r.qend]) }));
      should([mb.status, res.kept, res.copies.length]).eql(['ok', 0, 0]);
      should(mb.rows.find((r) => r.sseqid === 'scaffold2100')).match({ qstart: 160, qend: 431 });
    } finally {
      c.made.cleanup();
    }
  });
});
