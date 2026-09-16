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
    const undef = new annotate.Annotator({ mongo: { genes: { mongoCollection: async () => undefined } } });
    await undef.transcriptGenes(BICOLOR, ['X.1']);
    should([undef.unavailable, undef.mongoFailed]).eql([true, true]);
    const failing = new annotate.Annotator({ mongo: fakeMongo(GENES, { fail: true }).mongo });
    should((await failing.transcriptGenes(BICOLOR, ['SORBI_3004G087700.1'])).available).be.false();
    should(failing.mongoFailed).be.true();
    const ok = new annotate.Annotator({ mongo: fakeMongo(GENES).mongo });
    should((await ok.overlaps({ system_name: 'x', map_id: null }, [makeAmplicon({})])).available).be.true(); // no map id: harmless
    should([ok.unavailable, ok.mongoFailed]).eql([false, false]);
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
