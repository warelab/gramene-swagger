'use strict';

// annotate.js: overlap queries, transcript → gene mapping, gene groups (spec §B.9, §B.11).
// Offline: mongo is a fake collection that evaluates the small query subset the module uses.

const { describe, it } = require('node:test');
const should = require('should');

const annotate = require('../../../api/helpers/primers/check/annotate');
const { makeAmplicon } = require('../fixtures/verdicts/builders');
const { fakeMongo } = require('../fixtures/verdicts/fake_mongo');

const BICOLOR = { system_name: 'sorghum_bicolor', map_id: 'GCA_000003195.3', taxon_id: 4558006 };

function geneDoc(id, region, start, end, strand, transcripts, extra) {
  return Object.assign({
    _id: id,
    name: id,
    biotype: 'protein_coding',
    taxon_id: 4558006,
    system_name: 'sorghum_bicolor',
    location: { map: 'GCA_000003195.3', region, start, end, strand },
    gene_structure: { transcripts: (transcripts || []).map((t) => ({ id: t })) }
  }, extra || {});
}

const GENES = [
  geneDoc('SORBI_3004G087700', '4', 7421357, 7428285, 1, ['SORBI_3004G087700.3', 'SORBI_3004G087700.1', 'SORBI_3004G087700.2']),
  geneDoc('SORBI_3004G087800', '4', 7435120, 7442061, 1, ['SORBI_3004G087800.1']),
  geneDoc('SORBI_3005G183900', '5', 66885571, 66892778, -1, ['SORBI_3005G183900.1', 'SORBI_3005G183900.2']),
  // transcript id that does not follow the <gene>.<n> convention: only the taxon query finds it
  geneDoc('ODDGENE', '6', 100, 900, 1, ['odd_model_A']),
  // far upstream (> 1 Mb before) but long: excluded by the indexed window, as specified
  geneDoc('HUGE', '4', 6000000, 7500000, 1, [])
];

describe('annotate query builders', () => {
  it('overlapQuery is the §B.11 indexed window query', () => {
    should(annotate.overlapQuery('GCA_000003195.3', '4', 7423537, 7423746)).eql({
      'location.map': 'GCA_000003195.3',
      'location.region': '4',
      'location.start': { $lte: 7423746, $gte: 6423537 },
      'location.end': { $gte: 7423537 }
    });
  });

  it('transcriptQuery is the §B.9 taxon query', () => {
    should(annotate.transcriptQuery(4558006, ['a.1'])).eql({ taxon_id: 4558006, 'gene_structure.transcripts.id': { $in: ['a.1'] } });
  });

  it('stripIsoform removes a trailing numeric isoform suffix only', () => {
    should(annotate.stripIsoform('SORBI_3004G087700.3')).equal('SORBI_3004G087700');
    should(annotate.stripIsoform('SbiGrassl.01g028090.01')).equal('SbiGrassl.01g028090');
    should(annotate.stripIsoform('odd_model_A')).equal('odd_model_A');
  });

  it('stripIsoform maps the pan-genome transcript id formats to their gene ids (chk-strip-isoform-grouping)', () => {
    // transcript → gene _id pairs as found in the cDNA BLAST DB titles and mongo genes (survey of 120 genomes)
    const cases = [
      ['SbiPI561072.01g062000_T001', 'SbiPI561072.01g062000'],
      ['SbiGRIF16309.01g004960_T008', 'SbiGRIF16309.01g004960'],
      ['SbPI513676.01G013100.11.v1.1', 'SbPI513676.01G013100.v1.1'],
      ['SbPI154987.01G439300.1.v2.1', 'SbPI154987.01G439300.v2.1'],
      ['Sobic.001G000100.1.v5.1', 'Sobic.001G000100.v5.1'],
      ['BTx623.01G000100.t1', 'BTx623.01G000100'],
      ['BTx623.01G511500.t1.1.66327aeb', 'BTx623.01G511500'],
      ['Ji2055.07G218300.t2.2.66328bfb', 'Ji2055.07G218300']
    ];
    for (const [tx, gene] of cases) should(annotate.stripIsoform(tx)).equal(gene);
  });
});

describe('annotate.geneContext / orthologInfo', () => {
  it('normalizes resolved.gene and groups ortholog ids by system_name', () => {
    const g = annotate.geneContext({
      id: 'SORBI_3001G046200',
      location: { region: '1', start: 3406158, end: 3409981, strand: 1, map: 'GCA_000003195.3' },
      transcripts: ['SORBI_3001G046200.1', { id: 'SORBI_3001G046200.2' }],
      orthologs: { sorghum_leoti: ['SbiLeoti.01g006350'], sorghum_grassl: [{ id: 'SbiGrassl.01g006190' }] }
    });
    should(g).eql({
      id: 'SORBI_3001G046200',
      location: { region: '1', start: 3406158, end: 3409981, strand: 1, map: 'GCA_000003195.3' },
      transcripts: ['SORBI_3001G046200.1', 'SORBI_3001G046200.2'],
      orthologs: { sorghum_leoti: ['SbiLeoti.01g006350'], sorghum_grassl: ['SbiGrassl.01g006190'] }
    });
    should(annotate.geneContext({ _id: 'X', transcript_ids: ['X.1'], ortholog_ids: { s: ['Y'] } })).eql({ id: 'X', location: null, transcripts: ['X.1'], orthologs: { s: ['Y'] } });
    should(annotate.geneContext(null)).equal(null);
  });

  it('ortholog_annotated is true/false per genome and null without a design gene', () => {
    const g = annotate.geneContext({ id: 'G', orthologs: { sorghum_leoti: ['L1'], sorghum_353: [] } });
    const leoti = annotate.orthologInfo(g, 'sorghum_leoti');
    should(leoti.annotated).be.true();
    should(leoti.ids.has('L1')).be.true();
    should(annotate.orthologInfo(g, 'sorghum_353').annotated).be.false();
    should(annotate.orthologInfo(g, 'sorghum_is12661').annotated).be.false();
    should(annotate.orthologInfo(null, 'sorghum_leoti')).eql({ annotated: null, ids: null });
  });
});

describe('annotate.groupByGene (§B.9)', () => {
  const amps = [
    makeAmplicon({ region: 'SORBI_3004G087700.1', start: 1001, end: 1278, orientation: 'LR' }),
    makeAmplicon({ region: 'SORBI_3004G087700.3', start: 1090, end: 1367, orientation: 'LR' }),
    makeAmplicon({ region: 'SORBI_3004G087700.2', start: 890, end: 1167, orientation: 'LR' }),
    makeAmplicon({ region: 'SORBI_3005G183900.2', start: 500, end: 777, orientation: 'LR', left_mm: 2, left_mm_pos: [10, 7] }),
    makeAmplicon({ region: 'SORBI_3005G183900.1', start: 520, end: 800, orientation: 'LR', left_mm: 3, left_mm_pos: [12, 10, 7], likelihood: 'likely_weak', terminal_mismatch: true }),
    makeAmplicon({ region: 'SORBI_3004G087800.1', start: 40, end: 317, orientation: 'LR', right_mm: 1, right_mm_pos: [15] }),
    makeAmplicon({ region: 'SORBI_3004G087800.1', start: 60, end: 88, orientation: 'RR' })
  ];

  it('collapses isoforms per (gene, orientation) and describes the group by its best product', () => {
    const map = new Map([
      ['SORBI_3004G087700.1', 'SORBI_3004G087700'], ['SORBI_3004G087700.2', 'SORBI_3004G087700'], ['SORBI_3004G087700.3', 'SORBI_3004G087700'],
      ['SORBI_3005G183900.1', 'SORBI_3005G183900'], ['SORBI_3005G183900.2', 'SORBI_3005G183900'], ['SORBI_3004G087800.1', 'SORBI_3004G087800']
    ]);
    const groups = annotate.groupByGene(amps, map);
    should(groups.map((g) => g.gene_id + '/' + g.orientation)).eql([
      'SORBI_3004G087700/LR', 'SORBI_3004G087800/RR', 'SORBI_3004G087800/LR', 'SORBI_3005G183900/LR'
    ]);
    const on = groups[0];
    should(on.isoforms.map((i) => i.transcript_id)).eql(['SORBI_3004G087700.1', 'SORBI_3004G087700.2', 'SORBI_3004G087700.3']);
    should(on.isoforms.map((i) => i.size)).eql([278, 278, 278]);
    should(on).match({ size_min: 278, size_max: 278, likelihood: 'likely', left_mm: 0, right_mm: 0, approx: false });
    const para = groups[3];
    should(para).match({ size_min: 278, size_max: 281, likelihood: 'likely', left_mm: 2, left_mm_pos: [10, 7], terminal_mismatch: false });
    should(para.isoforms).have.length(2);
  });

  it('without a mapping, _T / .<n>.v<x>.<y> isoforms still collapse into one gene (chk-strip-isoform-grouping)', () => {
    const tx = [];
    for (let i = 1; i <= 10; i++) tx.push('SbiPI561073.01g000100_T' + String(i).padStart(3, '0'));
    const iso = tx.map((id, i) => makeAmplicon({ region: id, start: 100 + i, end: 377 + i }));
    should(annotate.stripIsoform('SbiPI561073.01g000100.1')).equal('SbiPI561073.01g000100'); // the old rule, for contrast
    should(tx.map((id) => id.replace(/\.\d+$/, '')).filter((v, i, a) => a.indexOf(v) === i)).have.length(10);
    should(annotate.groupByGene(iso, new Map()).map((g) => [g.gene_id, g.isoforms.length])).eql([['SbiPI561073.01g000100', 10]]);
    const v = ['SbPI513676.01G013100.1.v1.1', 'SbPI513676.01G013100.11.v1.1'].map((id) => makeAmplicon({ region: id }));
    should(annotate.groupByGene(v, {}).map((g) => [g.gene_id, g.isoforms.length])).eql([['SbPI513676.01G013100.v1.1', 2]]);
  });

  it('unmapped transcripts form their own gene without the isoform suffix; lookups may be objects or functions', () => {
    const byObject = annotate.groupByGene(amps.slice(0, 3), {});
    should(byObject.map((g) => g.gene_id)).eql(['SORBI_3004G087700']);
    const byFn = annotate.groupByGene(amps.slice(3, 5), () => 'PARALOG');
    should(byFn.map((g) => g.gene_id)).eql(['PARALOG']);
    should(annotate.groupByGene([], new Map())).eql([]);
  });
});

describe('Annotator.overlaps (§B.11)', () => {
  it('queries each distinct product once with the indexed window and assigns sorted genes to duplicates', async () => {
    const f = fakeMongo(GENES);
    const ann = new annotate.Annotator({ mongo: f.mongo });
    const a = makeAmplicon({ region: '4', start: 7423537, end: 7423746 });
    const dup = makeAmplicon({ region: '4', start: 7423537, end: 7423746, orientation: 'RR' });
    const b = makeAmplicon({ region: '4', start: 7428000, end: 7435500 });
    const c = makeAmplicon({ region: '9', start: 1, end: 100 });
    const res = await ann.overlaps(BICOLOR, [a, dup, b, c]);
    should(res).eql({ available: true, annotated: 3, skipped: 0 });
    should(f.calls).have.length(3);
    should(f.calls[0].options).eql({ fields: annotate.OVERLAP_FIELDS });
    should(f.calls[0].limit).equal(20);
    should(f.calls[0].query).eql(annotate.overlapQuery('GCA_000003195.3', '4', 7423537, 7423746));
    should(a.genes).eql([{ id: 'SORBI_3004G087700', name: 'SORBI_3004G087700', biotype: 'protein_coding', strand: 1 }]);
    should(dup.genes).equal(a.genes);
    should(b.genes.map((g) => g.id)).eql(['SORBI_3004G087700', 'SORBI_3004G087800']);
    should(c.genes).eql([]);
  });

  it('annotates at most maxAnnotated distinct products, in the given order; the rest get genes: null', async () => {
    const f = fakeMongo(GENES);
    const ann = new annotate.Annotator({ mongo: f.mongo, maxAnnotated: 2 });
    const list = [0, 1, 2, 3].map((k) => makeAmplicon({ region: '4', start: 7423537 + k * 10, end: 7423746 + k * 10 }));
    const res = await ann.overlaps(BICOLOR, list);
    should(res).eql({ available: true, annotated: 2, skipped: 2 });
    should(f.calls).have.length(2);
    should(list[0].genes).have.length(1);
    should(list[1].genes).have.length(1);
    should(list[2].genes).equal(null);
    should(list[3].genes).equal(null);
  });

  it('reports unavailable when there is no collection, and when a query fails or times out', async () => {
    const noMongo = new annotate.Annotator({ mongo: null });
    const a = makeAmplicon({ region: '4', start: 7423537, end: 7423746 });
    should(await noMongo.overlaps(BICOLOR, [a])).match({ available: false, annotated: 0 });
    should(a.genes).equal(null);

    const undef = new annotate.Annotator({ mongo: { genes: { mongoCollection: async () => undefined } } });
    should((await undef.overlaps(BICOLOR, [makeAmplicon({})])).available).be.false();
    should(undef.unavailable).be.true();

    const failing = fakeMongo(GENES, { fail: true });
    const ann = new annotate.Annotator({ mongo: failing.mongo });
    const list = [makeAmplicon({ region: '4', start: 1, end: 100 }), makeAmplicon({ region: '4', start: 7423537, end: 7423746 })];
    should(await ann.overlaps(BICOLOR, list)).match({ available: false, annotated: 0 });
    should(list.map((x) => x.genes)).eql([null, null]);
    // once unavailable, the job does not hit mongo again
    const before = failing.calls.length;
    should((await ann.overlaps(BICOLOR, [makeAmplicon({})])).available).be.false();
    should(failing.calls.length).equal(before);

    const hanging = fakeMongo(GENES, { hang: true });
    const slow = new annotate.Annotator({ mongo: hanging.mongo, timeoutMs: 30 });
    should((await slow.overlaps(BICOLOR, [makeAmplicon({ region: '4', start: 7423537, end: 7423746 })])).available).be.false();
    // a timed-out query is retried once; timing out again is a timeout, not a mongo failure
    should(hanging.calls).have.length(2);
    should([slow.unavailable, slow.mongoFailed, slow.timedOut, slow.timeoutRetries]).eql([true, false, true, 1]);
  });
});

describe('Annotator.transcriptGenes (§B.9)', () => {
  it('maps by _id lookup first, falls back to the taxon query, reports unmapped, caches per genome', async () => {
    const f = fakeMongo(GENES);
    const ann = new annotate.Annotator({ mongo: f.mongo });
    const subjects = ['SORBI_3004G087700.1', 'SORBI_3005G183900.2', 'odd_model_A', 'NOPE_1.1', 'SORBI_3004G087700.1'];
    const res = await ann.transcriptGenes(BICOLOR, subjects);
    should(res.available).be.true();
    should(Array.from(res.map.entries())).eql([
      ['SORBI_3004G087700.1', 'SORBI_3004G087700'], ['SORBI_3005G183900.2', 'SORBI_3005G183900'], ['odd_model_A', 'ODDGENE']
    ]);
    should(res.unmapped).eql(['NOPE_1.1']);
    should(f.calls).have.length(2);
    should(f.calls[0].query).eql({ _id: { $in: ['SORBI_3004G087700', 'SORBI_3005G183900', 'odd_model_A', 'NOPE_1'] } });
    should(f.calls[0].options).eql({ fields: annotate.TRANSCRIPT_FIELDS });
    should(f.calls[1].query).eql(annotate.transcriptQuery(4558006, ['odd_model_A', 'NOPE_1.1']));
    const again = await ann.transcriptGenes(BICOLOR, ['SORBI_3004G087700.1', 'NOPE_1.1']);
    should(f.calls).have.length(2);
    should(again.unmapped).eql(['NOPE_1.1']);
  });

  it('the indexed _id lookup finds _T isoforms without the taxon fallback query', async () => {
    const f = fakeMongo([geneDoc('SbiPI561072.01g062000', '1', 1, 5000, 1, ['SbiPI561072.01g062000_T001', 'SbiPI561072.01g062000_T002'], { taxon_id: 4558064 })]);
    const ann = new annotate.Annotator({ mongo: f.mongo });
    const res = await ann.transcriptGenes({ system_name: 'sorghum_pi561072', taxon_id: 4558064 }, ['SbiPI561072.01g062000_T002', 'SbiPI561072.01g062000_T001']);
    should(Array.from(res.map.entries())).eql([['SbiPI561072.01g062000_T002', 'SbiPI561072.01g062000'], ['SbiPI561072.01g062000_T001', 'SbiPI561072.01g062000']]);
    should(res.unmapped).eql([]);
    should(f.calls).have.length(1);
    should(f.calls[0].query).eql({ _id: { $in: ['SbiPI561072.01g062000'] } });
  });

  it('mongoFailed is set only when a supplied mongo handle fails (ops-mongo-fatal-unreachable)', async () => {
    const none = new annotate.Annotator({ mongo: null });
    await none.overlaps(BICOLOR, [makeAmplicon({})]);
    should([none.unavailable, none.mongoFailed]).eql([true, false]);
    should(none.timedOut).be.false();
    const undef = new annotate.Annotator({ mongo: { genes: { mongoCollection: async () => undefined } } });
    await undef.transcriptGenes(BICOLOR, ['X.1']);
    should([undef.unavailable, undef.mongoFailed]).eql([true, true]);
    should(undef.timedOut).be.false();
    const failing = new annotate.Annotator({ mongo: fakeMongo(GENES, { fail: true }).mongo });
    should((await failing.transcriptGenes(BICOLOR, ['SORBI_3004G087700.1'])).available).be.false();
    should(failing.mongoFailed).be.true();
    should([failing.timedOut, failing.timeoutRetries]).eql([false, 0]);
    const ok = new annotate.Annotator({ mongo: fakeMongo(GENES).mongo });
    should((await ok.overlaps({ system_name: 'x', map_id: null }, [makeAmplicon({})])).available).be.true(); // no map id: harmless
    should([ok.unavailable, ok.mongoFailed]).eql([false, false]);
    should([ok.timedOut, ok.timeoutRetries]).eql([false, 0]);
  });

  it('batches lookups by 500', async () => {
    const f = fakeMongo([]);
    const ann = new annotate.Annotator({ mongo: f.mongo });
    const subjects = [];
    for (let i = 0; i < 1201; i++) subjects.push('G' + i + '.1');
    const res = await ann.transcriptGenes(BICOLOR, subjects);
    should(res.unmapped).have.length(1201);
    should(f.calls.map((c) => (c.query._id ? c.query._id.$in.length : c.query['gene_structure.transcripts.id'].$in.length))).eql([500, 500, 201, 500, 500, 201]);
  });

  it('uses the seed without mongo and reports unavailable instead of unmapped', async () => {
    const ann = new annotate.Annotator({ mongo: null });
    const res = await ann.transcriptGenes(BICOLOR, ['SORBI_3004G087700.2', 'SORBI_3005G183900.1'], new Map([['SORBI_3004G087700.2', 'SORBI_3004G087700']]));
    should(res.available).be.false();
    should(Array.from(res.map.entries())).eql([['SORBI_3004G087700.2', 'SORBI_3004G087700']]);
    should(res.unmapped).eql([]);
  });
});

describe('Annotator query timeouts: retried once, then unavailable but not a mongo failure (ops-mongo-timeout-nonfatal)', () => {
  const T = 20; // timeoutMs: the retry gets 40 ms
  const ON = { region: '4', start: 7423537, end: 7423746 };
  const recordingLog = () => {
    const warns = [];
    return { warns, log: { info() {}, warn: (m) => warns.push(m), error() {} } };
  };
  const flags = (ann) => ({ unavailable: ann.unavailable, mongoFailed: ann.mongoFailed, timedOut: ann.timedOut, timeoutRetries: ann.timeoutRetries });
  // resolves with the matching docs after ms: slower than the first budget, within the retry's
  const late = (ms, docs) => new Promise((resolve) => setTimeout(() => resolve(docs), ms));

  it('overlaps: a timeout then a success on the retry (twice the budget) keeps the annotation', async () => {
    const { warns, log } = recordingLog();
    const f = fakeMongo(GENES, { script: (call) => (call.n === 0 ? 'hang' : call.n === 1 ? late(30, [GENES[0]]) : 'ok') });
    const ann = new annotate.Annotator({ mongo: f.mongo, timeoutMs: T, log });
    const a = makeAmplicon(ON);
    should(await ann.overlaps(BICOLOR, [a])).eql({ available: true, annotated: 1, skipped: 0 });
    should(a.genes.map((g) => g.id)).eql(['SORBI_3004G087700']);
    should(flags(ann)).eql({ unavailable: false, mongoFailed: false, timedOut: false, timeoutRetries: 1 });
    // the same query again, with the same fields and limit
    should(f.calls).have.length(2);
    should(f.calls[1].query).eql(f.calls[0].query);
    should(f.calls[1].options).eql({ fields: annotate.OVERLAP_FIELDS });
    should(f.calls[1].limit).equal(20);
    should(warns).eql(['primers check: mongo genes query timed out after 20 ms; retrying it once with a 40 ms budget']);
    // later queries get the normal budget and no retry
    const b = makeAmplicon({ region: '4', start: 7428000, end: 7435500 });
    should((await ann.overlaps(BICOLOR, [b])).available).be.true();
    should(b.genes.map((g) => g.id)).eql(['SORBI_3004G087700', 'SORBI_3004G087800']);
    should([f.calls.length, ann.timeoutRetries]).eql([3, 1]);
  });

  it('overlaps: two timeouts make the annotation unavailable for the job, with timedOut and without mongoFailed', async () => {
    const { warns, log } = recordingLog();
    const f = fakeMongo(GENES, { hang: true });
    const ann = new annotate.Annotator({ mongo: f.mongo, timeoutMs: T, log, concurrency: 1 });
    const list = [makeAmplicon(ON), makeAmplicon({ region: '4', start: 7428000, end: 7435500 })];
    const res = await ann.overlaps(BICOLOR, list);
    should(res).eql({ available: false, annotated: 0, skipped: 2 });
    should(list.map((x) => x.genes)).eql([null, null]);
    should(flags(ann)).eql({ unavailable: true, mongoFailed: false, timedOut: true, timeoutRetries: 1 });
    should(f.calls).have.length(2); // the first product and its retry; the second product is not queried
    should(warns).eql([
      'primers check: mongo genes query timed out after 20 ms; retrying it once with a 40 ms budget',
      'primers check: annotation unavailable: mongo genes query timed out after 40 ms'
    ]);
    // the rest of the job does not hit mongo and stays a timeout
    should((await ann.overlaps(BICOLOR, [makeAmplicon(ON)])).available).be.false();
    should((await ann.transcriptGenes(BICOLOR, ['SORBI_3004G087700.1'])).available).be.false();
    should(f.calls).have.length(2);
    should(flags(ann)).eql({ unavailable: true, mongoFailed: false, timedOut: true, timeoutRetries: 1 });
  });

  it('transcriptGenes: a timeout then a success on the retry maps the transcripts', async () => {
    const f = fakeMongo(GENES, { script: (call) => (call.n === 0 ? 'hang' : 'ok') });
    const ann = new annotate.Annotator({ mongo: f.mongo, timeoutMs: T });
    const res = await ann.transcriptGenes(BICOLOR, ['SORBI_3004G087700.1', 'odd_model_A']);
    should(res.available).be.true();
    should(Array.from(res.map.entries())).eql([['SORBI_3004G087700.1', 'SORBI_3004G087700'], ['odd_model_A', 'ODDGENE']]);
    should(res.unmapped).eql([]);
    should(f.calls.map((c) => c.query)).eql([
      { _id: { $in: ['SORBI_3004G087700', 'odd_model_A'] } },
      { _id: { $in: ['SORBI_3004G087700', 'odd_model_A'] } },
      annotate.transcriptQuery(4558006, ['odd_model_A'])
    ]);
    should(flags(ann)).eql({ unavailable: false, mongoFailed: false, timedOut: false, timeoutRetries: 1 });
  });

  it('transcriptGenes: two timeouts (here in the taxon fallback query) are unavailable with timedOut, without mongoFailed', async () => {
    const f = fakeMongo(GENES, { script: (call) => (call.n === 0 ? 'ok' : 'hang') });
    const ann = new annotate.Annotator({ mongo: f.mongo, timeoutMs: T });
    const res = await ann.transcriptGenes(BICOLOR, ['SORBI_3004G087700.1', 'odd_model_A']);
    should(res.available).be.false();
    should(res.unmapped).eql([]);
    should(f.calls).have.length(3);
    should(f.calls[2].query).eql(f.calls[1].query);
    should(flags(ann)).eql({ unavailable: true, mongoFailed: false, timedOut: true, timeoutRetries: 1 });
  });

  it('connection-level failures still set mongoFailed: a collection() timeout, a missing collection, a non-timeout error (also on the retry)', async () => {
    const hungColl = new annotate.Annotator({ mongo: { genes: { mongoCollection: () => new Promise(() => {}) } }, timeoutMs: T });
    should((await hungColl.overlaps(BICOLOR, [makeAmplicon(ON)])).available).be.false();
    should(flags(hungColl)).eql({ unavailable: true, mongoFailed: true, timedOut: false, timeoutRetries: 0 });

    const noColl = new annotate.Annotator({ mongo: fakeMongo(GENES, { noCollection: true }).mongo, timeoutMs: T });
    should((await noColl.transcriptGenes(BICOLOR, ['SORBI_3004G087700.1'])).available).be.false();
    should(flags(noColl)).eql({ unavailable: true, mongoFailed: true, timedOut: false, timeoutRetries: 0 });

    const err = fakeMongo(GENES, { script: () => new Error('Topology was destroyed') });
    const errAnn = new annotate.Annotator({ mongo: err.mongo, timeoutMs: T });
    should((await errAnn.overlaps(BICOLOR, [makeAmplicon(ON)])).available).be.false();
    should(err.calls).have.length(1); // not retried
    should(flags(errAnn)).eql({ unavailable: true, mongoFailed: true, timedOut: false, timeoutRetries: 0 });

    const retryErr = fakeMongo(GENES, { script: (call) => (call.n === 0 ? 'hang' : new Error('Topology was destroyed')) });
    const retryAnn = new annotate.Annotator({ mongo: retryErr.mongo, timeoutMs: T });
    should((await retryAnn.transcriptGenes(BICOLOR, ['SORBI_3004G087700.1'])).available).be.false();
    should(flags(retryAnn)).eql({ unavailable: true, mongoFailed: true, timedOut: false, timeoutRetries: 1 });
  });

  it('a timeout and a connection-level failure in one job: mongoFailed ends up true in either order', async () => {
    // timeout first: both products' first attempts time out and are retried (same budget, the ON product dispatched
    // first, so its retry times out first); the other product's retry then fails with a connection error
    const { warns, log } = recordingLog();
    let failOther = null;
    const first = fakeMongo(GENES, {
      script: (call) => (call.query['location.start'].$lte === ON.end || call.n < 2 ? 'hang' : new Promise((resolve, reject) => { failOther = reject; }))
    });
    const timeoutLog = { info() {}, error() {}, warn: (m) => { log.warn(m); if (/annotation unavailable: mongo genes query timed out/.test(m)) failOther(new Error('Topology was destroyed')); } };
    const a = new annotate.Annotator({ mongo: first.mongo, timeoutMs: T, log: timeoutLog });
    should((await a.overlaps(BICOLOR, [makeAmplicon(ON), makeAmplicon({ region: '4', start: 7428000, end: 7435500 })])).available).be.false();
    should(first.calls).have.length(4);
    should(flags(a)).eql({ unavailable: true, mongoFailed: true, timedOut: true, timeoutRetries: 2 });
    should(warns.filter((m) => /annotation unavailable/.test(m))).eql([
      'primers check: annotation unavailable: mongo genes query timed out after 40 ms',
      'primers check: annotation unavailable: Topology was destroyed' // the fatal cause is logged too
    ]);

    // failure first: a query still in flight that then times out is not retried
    const second = fakeMongo(GENES, { script: (call) => (call.n === 0 ? new Error('Topology was destroyed') : 'hang') });
    const b = new annotate.Annotator({ mongo: second.mongo, timeoutMs: T });
    should((await b.overlaps(BICOLOR, [makeAmplicon(ON), makeAmplicon({ region: '4', start: 7428000, end: 7435500 })])).available).be.false();
    should(second.calls).have.length(2);
    should(flags(b)).eql({ unavailable: true, mongoFailed: true, timedOut: true, timeoutRetries: 0 });
  });
});

describe('Annotator query slots: at most MAX_IN_FLIGHT genes queries per mongo handle, the retry race keeps the first attempt (ops-mongo-timeout-nonfatal)', () => {
  const T = 20;
  const ON = { region: '4', start: 7423537, end: 7423746 };
  const RETRYING = 'primers check: mongo genes query timed out after 20 ms; retrying it once with a 40 ms budget';
  const WAITING = 'primers check: mongo genes query timed out after 20 ms waiting for one of the 5 genes query slots; waiting once more with a 40 ms budget';
  const GAVE_UP = 'primers check: annotation unavailable: mongo genes query timed out after 40 ms';
  const recordingLog = () => {
    const warns = [];
    return { warns, log: { info() {}, warn: (m) => warns.push(m), error() {} } };
  };
  const flags = (ann) => ({ unavailable: ann.unavailable, mongoFailed: ann.mongoFailed, timedOut: ann.timedOut, timeoutRetries: ann.timeoutRetries });
  const late = (ms, docs) => new Promise((resolve) => setTimeout(() => resolve(docs), ms));
  const turn = () => new Promise((resolve) => setImmediate(resolve));
  // n distinct products overlapping SORBI_3004G087700
  const products = (n) => Array.from({ length: n }, (_, k) => makeAmplicon({ region: '4', start: 7423537 + k * 10, end: 7423746 + k * 10 }));

  it('the first attempt stays in the race: its answer after timeoutMs is used while the resend is still running', async () => {
    const { warns, log } = recordingLog();
    const f = fakeMongo(GENES, { script: (call) => (call.n === 0 ? late(30, [GENES[0]]) : 'hang') });
    const ann = new annotate.Annotator({ mongo: f.mongo, timeoutMs: T, log });
    const a = makeAmplicon(ON);
    should(await ann.overlaps(BICOLOR, [a])).eql({ available: true, annotated: 1, skipped: 0 });
    should(a.genes.map((g) => g.id)).eql(['SORBI_3004G087700']);
    should(f.calls).have.length(2); // the resend was sent and never answered
    should(f.calls[1].query).eql(f.calls[0].query);
    should(flags(ann)).eql({ unavailable: false, mongoFailed: false, timedOut: false, timeoutRetries: 1 });
    should(warns).eql([RETRYING]);
  });

  it('sends at most 5 genes queries at a time; the others follow in order as answers free the slots', async () => {
    should(annotate.MAX_IN_FLIGHT).equal(5);
    const answers = [];
    const f = fakeMongo(GENES, { script: () => new Promise((resolve) => answers.push(resolve)) });
    const ann = new annotate.Annotator({ mongo: f.mongo, timeoutMs: 5000 });
    const done = ann.overlaps(BICOLOR, products(7));
    await turn();
    should(f.calls).have.length(5);
    answers[0]([]);
    await turn();
    await turn();
    should(f.calls).have.length(6);
    should(f.calls[5].query['location.start'].$lte).equal(7423746 + 50); // the sixth product
    for (let i = 1; i < 7; i++) {
      await turn();
      await turn();
      answers[i]([]);
    }
    should(await done).eql({ available: true, annotated: 7, skipped: 0 });
    should([f.calls.length, ann.queries, ann.timeoutRetries]).eql([7, 7, 0]);
  });

  it('timed-out queries keep their slots until mongo answers: queries still waiting are dropped unsent, and the next job queues instead of adding a backlog', async () => {
    const stuck = [];
    const f = fakeMongo(GENES, { script: (call) => (call.n < 5 ? new Promise((resolve) => stuck.push(resolve)) : 'ok') });

    // job 1: 5 queries sent and never answered, 2 waiting; the 5 resends and the 2 waiting queries are never sent
    const j1 = recordingLog();
    const job1 = new annotate.Annotator({ mongo: f.mongo, timeoutMs: T, log: j1.log });
    should(await job1.overlaps(BICOLOR, products(7))).eql({ available: false, annotated: 0, skipped: 7 });
    should(f.calls).have.length(5);
    should([job1.queries, stuck.length]).eql([5, 5]);
    should(flags(job1)).eql({ unavailable: true, mongoFailed: false, timedOut: true, timeoutRetries: 7 });
    should(j1.warns).eql([RETRYING, RETRYING, RETRYING, RETRYING, RETRYING, WAITING, WAITING, GAVE_UP]);

    // job 2 on the same handle: its query waits for a slot and is never sent
    const j2 = recordingLog();
    const job2 = new annotate.Annotator({ mongo: f.mongo, timeoutMs: T, log: j2.log });
    should((await job2.overlaps(BICOLOR, products(1))).available).be.false();
    should(f.calls).have.length(5);
    should(job2.queries).equal(0);
    should(flags(job2)).eql({ unavailable: true, mongoFailed: false, timedOut: true, timeoutRetries: 1 });
    should(j2.warns).eql([WAITING, GAVE_UP]);

    // mongo answers the abandoned queries: their slots are free again and the next job is annotated
    stuck.forEach((resolve) => resolve([]));
    const job3 = new annotate.Annotator({ mongo: f.mongo, timeoutMs: T });
    const a = makeAmplicon(ON);
    should(await job3.overlaps(BICOLOR, [a])).eql({ available: true, annotated: 1, skipped: 0 });
    should(a.genes.map((g) => g.id)).eql(['SORBI_3004G087700']);
    should(f.calls).have.length(6);
    should(flags(job3)).eql({ unavailable: false, mongoFailed: false, timedOut: false, timeoutRetries: 0 });
  });

  it('a resend waiting for a slot is not sent once the first attempt answers; a query that waited past timeoutMs for a slot still gets its answer', async () => {
    const held = [];
    const f = fakeMongo(GENES, {
      script: (call) => (call.n < 4 || call.n === 5 ? new Promise((resolve) => held.push(resolve)) : call.n === 4 ? late(30, [GENES[0]]) : 'ok')
    });
    // another job holds 4 of the 5 slots
    const other = new annotate.Annotator({ mongo: f.mongo, timeoutMs: 5000 });
    const otherDone = other.overlaps(BICOLOR, products(4));
    await turn();
    should(f.calls).have.length(4);

    const r = recordingLog();
    const ann = new annotate.Annotator({ mongo: f.mongo, timeoutMs: T, log: r.log });
    const a = makeAmplicon(ON);
    should(await ann.overlaps(BICOLOR, [a])).eql({ available: true, annotated: 1, skipped: 0 });
    should(a.genes.map((g) => g.id)).eql(['SORBI_3004G087700']);
    should(f.calls).have.length(5); // its resend waited for a slot and was dropped when the first attempt answered
    should(flags(ann)).eql({ unavailable: false, mongoFailed: false, timedOut: false, timeoutRetries: 1 });
    should(r.warns).eql([RETRYING]);

    // a third long job takes the free slot; the next query waits past its 20 ms, gets a slot at 30 ms and its answer
    const third = new annotate.Annotator({ mongo: f.mongo, timeoutMs: 5000 });
    const thirdDone = third.overlaps(BICOLOR, products(1));
    await turn();
    await turn();
    should(f.calls).have.length(6);
    const w = recordingLog();
    const waiting = new annotate.Annotator({ mongo: f.mongo, timeoutMs: T, log: w.log });
    setTimeout(() => held[0]([]), 30);
    const b = makeAmplicon(ON);
    should(await waiting.overlaps(BICOLOR, [b])).eql({ available: true, annotated: 1, skipped: 0 });
    should(b.genes.map((g) => g.id)).eql(['SORBI_3004G087700']);
    should([f.calls.length, waiting.queries]).eql([7, 1]);
    should(flags(waiting)).eql({ unavailable: false, mongoFailed: false, timedOut: false, timeoutRetries: 1 });
    should(w.warns).eql([WAITING]);

    held.forEach((resolve) => resolve([]));
    should((await otherDone).available).be.true();
    should((await thirdDone).available).be.true();
    should(f.calls).have.length(7);
  });

  it('a connection-level failure drops the queries still waiting for a slot: they are never sent', async () => {
    const held = [];
    const f = fakeMongo(GENES, {
      script: (call) => (call.n < 4 ? new Promise((resolve) => held.push(resolve))
        : call.n === 4 ? new Promise((resolve, reject) => setTimeout(() => reject(new Error('Topology was destroyed')), 5)) : 'ok')
    });
    const ann = new annotate.Annotator({ mongo: f.mongo, timeoutMs: 5000 });
    const done = ann.overlaps(BICOLOR, products(7));
    await late(30);
    should(f.calls).have.length(5);
    should(flags(ann)).eql({ unavailable: true, mongoFailed: true, timedOut: false, timeoutRetries: 0 });
    // the 2 queries that were waiting left the slot queue at once, while the 4 slots are still held
    should([ann._slots.active, ann._slots.waiting.length, ann._waiting.size]).eql([4, 0, 0]);
    held.forEach((resolve) => resolve([]));
    // the dropped queries do not keep the call waiting for their 5 s budget
    const t = Date.now();
    should(await done).eql({ available: false, annotated: 0, skipped: 7 });
    should(Date.now() - t).be.below(1000);
    should([f.calls.length, ann.queries]).eql([5, 5]);
  });

  it('queries of a job still waiting for a slot leave the queue as soon as the job gives up, not at the end of their own budget', async () => {
    const held = [];
    const f = fakeMongo(GENES, { script: (call) => (call.n < 4 ? new Promise((resolve) => held.push(resolve)) : 'hang') });
    const other = new annotate.Annotator({ mongo: f.mongo, timeoutMs: 5000 });
    const otherDone = other.overlaps(BICOLOR, products(4));
    await turn();
    const ann = new annotate.Annotator({ mongo: f.mongo, timeoutMs: 40 });
    // A: sent at 0 ms (the 5th slot), never answered: retried at 40 ms (the resend waits), gives up at 120 ms
    const A = ann.overlaps(BICOLOR, [makeAmplicon(ON)]);
    // B: dispatched at 50 ms, waits for a slot; its own budget would end at 170 ms
    const B = late(50).then(() => ann.transcriptGenes(BICOLOR, ['SORBI_3004G087700.1']));
    should((await A).available).be.false();
    should(flags(ann)).eql({ unavailable: true, mongoFailed: false, timedOut: true, timeoutRetries: 2 });
    should([ann._slots.active, ann._slots.waiting.length, ann._waiting.size]).eql([5, 0, 0]);
    should((await B).available).be.false();
    should([f.calls.length, ann.queries]).eql([5, 1]);
    held.forEach((resolve) => resolve([]));
    should((await otherDone).available).be.true();
  });

  it('once another call made the annotator unavailable, transcriptGenes sends no further query', async () => {
    // the _id lookup answers after 20 ms with nothing, so the taxon fallback query would follow; meanwhile an overlap
    // query of the same job fails
    const f = fakeMongo(GENES, { script: (call) => (call.n === 0 ? late(20, []) : call.n === 1 ? new Error('Topology was destroyed') : 'ok') });
    const ann = new annotate.Annotator({ mongo: f.mongo, timeoutMs: 1000 });
    const tg = ann.transcriptGenes(BICOLOR, ['odd_model_A']);
    should((await ann.overlaps(BICOLOR, [makeAmplicon(ON)])).available).be.false();
    should(await tg).eql({ map: new Map(), unmapped: [], available: false });
    should(f.calls.map((c) => c.query)).eql([{ _id: { $in: ['odd_model_A'] } }, annotate.overlapQuery(BICOLOR.map_id, ON.region, ON.start, ON.end)]);
    should(flags(ann)).eql({ unavailable: true, mongoFailed: true, timedOut: false, timeoutRetries: 0 });
  });

  it('transcriptGenes: a query waiting for a slot past timeoutMs is not re-sent, and a second timeout is timedOut only', async () => {
    const held = [];
    const f = fakeMongo(GENES, { script: () => new Promise((resolve) => held.push(resolve)) });
    const other = new annotate.Annotator({ mongo: f.mongo, timeoutMs: 5000 });
    const otherDone = other.overlaps(BICOLOR, products(5));
    await turn();
    should(f.calls).have.length(5);
    const r = recordingLog();
    const ann = new annotate.Annotator({ mongo: f.mongo, timeoutMs: T, log: r.log });
    const res = await ann.transcriptGenes(BICOLOR, ['SORBI_3004G087700.1']);
    should(res).eql({ map: new Map(), unmapped: [], available: false });
    should([f.calls.length, ann.queries]).eql([5, 0]);
    should(flags(ann)).eql({ unavailable: true, mongoFailed: false, timedOut: true, timeoutRetries: 1 });
    should(r.warns).eql([WAITING, GAVE_UP]);
    held.forEach((resolve) => resolve([]));
    should((await otherDone).available).be.true();
    should(f.calls).have.length(5); // the dropped lookup was never sent
  });
});
