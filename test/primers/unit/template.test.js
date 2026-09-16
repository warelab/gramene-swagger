'use strict';

require('../../../api/helpers/primers/node_compat');

const test = require('node:test');
const should = require('should');
const fs = require('fs');
const path = require('path');

const template = require('../../../api/helpers/primers/template');
const coords = require('../../../api/helpers/primers/coords');
const { revcomp } = require('../../../api/helpers/primers/sequence');
const stubs = require('../fixtures/design/stubs');
const fx = require('../fixtures/catalog/fsfixture');

const G200 = 'SORBI_3001G000200';
const G700 = 'SORBI_3001G000700';

async function rejectsWith(promise, status, code) {
  let err = null;
  try {
    await promise;
  } catch (e) {
    err = e;
  }
  should.exist(err, 'expected a rejection with ' + code);
  err.should.have.property('code', code);
  err.should.have.property('status', status);
  return err;
}

// ---- synthetic genome: region 'c', 1000 bp ------------------------------------------------------

const C_DNA = '/fake/synthetic/dna/Syn.dna.toplevel.fa.gz';
const C_SEQ = fx.randomSeq(1000, 7);

function syntheticDeps(docs, extra) {
  const files = {};
  files[C_DNA] = { c: { length: 1000, windows: [{ start: 1, seq: C_SEQ }] } };
  return Object.assign({
    cfg: stubs.cfg(),
    log: stubs.silentLog,
    findGene: async function (id) { return docs[id] ? JSON.parse(JSON.stringify(docs[id])) : null; },
    resolve: async function (name) {
      return stubs.resolvedStub({ system_name: name, fasta: { dna: C_DNA, dna_sm: null } });
    },
    sequence: stubs.sequenceStub(files)
  }, extra || {});
}

function plusC(start, end) {
  return C_SEQ.slice(start - 1, end);
}

// Non-coding, 2 exons, no cds and no exon_junctions (AT3G05735.1-like), minus strand.
const NC1 = {
  _id: 'NC1',
  system_name: 'synthetic',
  location: { region: 'c', start: 101, end: 700, strand: -1 },
  gene_structure: {
    exons: [{ id: 'e1', start: 1, end: 150 }, { id: 'e2', start: 401, end: 600 }],
    transcripts: [{ id: 'NC1.1', length: 350, exons: ['e1', 'e2'] }],
    canonical_transcript: 'NC1.1'
  }
};

// ---- transcript mode ------------------------------------------------------------------------------

test('transcript mode (- strand): cDNA equals the captured SORBI_3001G000200.1 template, junctions equal exon_junctions', async function () {
  const deps = stubs.sorghumDeps();
  const t = await template.buildTemplate({ mode: 'transcript', gene_id: G200, transcript_id: null }, deps);
  const doc = stubs.geneDoc(G200);
  t.mode.should.equal('transcript');
  t.length.should.equal(1982);
  t.seq.should.equal(stubs.cdna200());
  t.transcript_id.should.equal('SORBI_3001G000200.1');
  t.junctions.should.eql(doc.gene_structure.transcripts[0].exon_junctions);
  t.features.junctions.should.eql(t.junctions);
  t.junctions.slice(0, 3).should.eql([397, 493, 597]);
  t.segments.should.eql(stubs.transcriptSegments200());
  should(t.features.gene).be.null();
  t.features.exons[0].should.eql({ id: 'EER93047-1', start: 1, end: 397, genomic: { start: 14503, end: 14899 } });
  t.features.cds.should.eql({ start: 299, end: 1651 });
  [t.region, t.start, t.end, t.strand].should.eql(['1', 11180, 14899, -1]);
  t.exon_count.should.equal(11);
  // the whole transcript is read with one envelope fetch
  deps.sequence.calls.length.should.equal(1);
  deps.sequence.calls[0].should.eql({ fastaPath: stubs.SB_DNA, region: '1', start: 11180, end: 14899, strand: -1 });
});

test('transcript mode: the spliced mapper reproduces the verified genomic blocks', async function () {
  const t = await template.buildTemplate({ mode: 'transcript', gene_id: G200 }, stubs.sorghumDeps());
  const left = t.mapper.primer('left', 856, 875);
  left.should.eql({ region: '1', start: 13395, end: 13650, strand: -1, blocks: [{ start: 13395, end: 13400 }, { start: 13637, end: 13650 }] });
  stubs.blocksSeq(G200, left.blocks, -1).should.equal('ATTACATCAAATAGGCCTTG');
  t.seq.slice(855, 875).should.equal('ATTACATCAAATAGGCCTTG');
  const right = t.mapper.primer('right', 960, 979);
  right.blocks.should.eql([{ start: 13291, end: 13310 }]);
  right.strand.should.equal(1);
  stubs.blocksSeq(G200, right.blocks, 1).should.equal('AACTTCTTTGTCGATCCATG');
});

test('transcript mode reads exon by exon when the gene envelope exceeds span_fetch_limit', async function () {
  const deps = stubs.sorghumDeps({ span_fetch_limit: 1000 });
  const t = await template.buildTemplate({ mode: 'transcript', gene_id: G200 }, deps);
  t.seq.should.equal(stubs.cdna200());
  deps.sequence.calls.length.should.equal(11);
});

test('transcript mode (+ strand, SORBI_3001G000700 canonical .3): junctions and cDNA from exon slices', async function () {
  const t = await template.buildTemplate({ mode: 'transcript', gene_id: G700 }, stubs.sorghumDeps());
  const doc = stubs.geneDoc(G700);
  const tr = doc.gene_structure.transcripts.find(function (x) { return x.id === 'SORBI_3001G000700.3'; });
  t.transcript_id.should.equal('SORBI_3001G000700.3');
  t.length.should.equal(tr.length);
  t.junctions.should.eql(tr.exon_junctions);
  const byId = new Map(doc.gene_structure.exons.map(function (e) { return [e.id, e]; }));
  const expected = tr.exons.map(function (id) {
    const g = coords.exonGenomicRange(doc.location, byId.get(id));
    return stubs.blocksSeq(G700, [g], 1);
  }).join('');
  t.seq.should.equal(expected);
  // a primer across the first junction maps to two ascending blocks on the + strand
  const p = t.mapper.primer('left', 80, 99);
  p.blocks.length.should.equal(2);
  p.strand.should.equal(1);
  stubs.blocksSeq(G700, p.blocks, 1).should.equal(t.seq.slice(79, 99));
});

test('explicit transcript_id selects a non-canonical transcript; an unknown one is 404 UNKNOWN_TRANSCRIPT', async function () {
  const t = await template.buildTemplate({ mode: 'transcript', gene_id: G700, transcript_id: 'SORBI_3001G000700.1' }, stubs.sorghumDeps());
  t.transcript_id.should.equal('SORBI_3001G000700.1');
  t.junctions.should.eql([325, 405, 611]);
  const err = await rejectsWith(template.buildTemplate({ mode: 'transcript', gene_id: G700, transcript_id: 'NOPE.1' }, stubs.sorghumDeps()),
    404, 'UNKNOWN_TRANSCRIPT');
  err.details.should.eql({ gene_id: G700, transcript_id: 'NOPE.1' });
});

test('non-coding multi-exon transcript without exon_junctions still gets derived junctions (- strand)', async function () {
  const t = await template.buildTemplate({ mode: 'transcript', gene_id: 'NC1' }, syntheticDeps({ NC1: NC1 }));
  t.length.should.equal(350);
  t.junctions.should.eql([150]);
  should(t.features.cds).be.null();
  t.seq.should.equal(revcomp(plusC(551, 700)) + revcomp(plusC(101, 300)));
  t.features.exons.should.eql([
    { id: 'e1', start: 1, end: 150, genomic: { start: 551, end: 700 } },
    { id: 'e2', start: 151, end: 350, genomic: { start: 101, end: 300 } }
  ]);
});

test('single-exon transcript has no junctions', async function () {
  const doc = JSON.parse(JSON.stringify(NC1));
  doc.gene_structure.transcripts = [{ id: 'NC1.2', length: 150, exons: ['e1'] }];
  doc.gene_structure.canonical_transcript = 'NC1.2';
  const t = await template.buildTemplate({ mode: 'transcript', gene_id: 'NC1' }, syntheticDeps({ NC1: doc }));
  t.junctions.should.eql([]);
  t.exon_count.should.equal(1);
});

test('transcript.length that disagrees with the exons is 500 GENE_STRUCTURE_MISMATCH', async function () {
  const doc = JSON.parse(JSON.stringify(NC1));
  doc.gene_structure.transcripts[0].length = 351;
  await rejectsWith(template.buildTemplate({ mode: 'transcript', gene_id: 'NC1' }, syntheticDeps({ NC1: doc })), 500, 'GENE_STRUCTURE_MISMATCH');
});

test('an exon id missing from gene_structure.exons, or outside the gene, is 500 GENE_STRUCTURE_MISMATCH', async function () {
  const missing = JSON.parse(JSON.stringify(NC1));
  missing.gene_structure.transcripts[0].exons = ['e1', 'e9'];
  await rejectsWith(template.buildTemplate({ mode: 'transcript', gene_id: 'NC1' }, syntheticDeps({ NC1: missing })), 500, 'GENE_STRUCTURE_MISMATCH');
  const outside = JSON.parse(JSON.stringify(NC1));
  outside.gene_structure.exons[1].end = 601;
  outside.gene_structure.transcripts[0].length = 351;
  await rejectsWith(template.buildTemplate({ mode: 'transcript', gene_id: 'NC1' }, syntheticDeps({ NC1: outside })), 500, 'GENE_STRUCTURE_MISMATCH');
});

// ---- gene mode ------------------------------------------------------------------------------------

test('gene mode flanks 200/100 (- strand): 1:11080-15099, 4020 bp, exon1 201-597, CDS start 499 at genomic 14601', async function () {
  const t = await template.buildTemplate({ mode: 'gene', gene_id: G200, flank_up: 200, flank_down: 100 }, stubs.sorghumDeps());
  [t.region, t.start, t.end, t.strand, t.length].should.eql(['1', 11080, 15099, -1, 4020]);
  t.seq.should.equal(revcomp(stubs.windowSeq(G200)));
  t.features.gene.should.eql({ start: 201, end: 3920 });
  t.features.exons[0].should.eql({ id: 'EER93047-1', start: 201, end: 597, genomic: { start: 14503, end: 14899 } });
  t.features.exons.length.should.equal(11);
  t.features.cds.start.should.equal(499);
  t.features.junctions.should.eql([]);
  t.junctions.should.eql([]);
  t.mapper.toGenomic(499).should.equal(14601);
  t.transcript_id.should.equal('SORBI_3001G000200.1');
  t.id.should.equal(G200);
});

test('gene mode (+ strand, SORBI_3001G000700) with flanks 100/100: template, exon overlay and CDS start', async function () {
  const t = await template.buildTemplate({ mode: 'gene', gene_id: G700, flank_up: 100, flank_down: 100 }, stubs.sorghumDeps());
  const doc = stubs.geneDoc(G700);
  [t.start, t.end, t.strand, t.length].should.eql([53681, 63405, 1, 9725]);
  t.seq.should.equal(stubs.windowSeq(G700));
  const tr = doc.gene_structure.transcripts[0];
  const first = doc.gene_structure.exons.find(function (e) { return e.id === tr.exons[0]; });
  t.features.exons[0].should.eql({ id: first.id, start: first.start + 100, end: first.end + 100,
    genomic: { start: 53781 + first.start - 1, end: 53781 + first.end - 1 } });
  // CDS starts at cDNA 1 of .3 -> the first exon's first base
  t.features.cds.start.should.equal(first.start + 100);
  t.mapper.toGenomic(t.features.cds.start).should.equal(53781 + first.start - 1);
});

test('flankExtent clamps to the region on both strands', function () {
  template.flankExtent({ start: 50, end: 400, strand: 1 }, 1000, 200, 700).should.eql({ start: 1, end: 1000, effUp: 49 });
  template.flankExtent({ start: 600, end: 950, strand: -1 }, 1000, 200, 700).should.eql({ start: 1, end: 1000, effUp: 50 });
  template.flankExtent({ start: 11180, end: 14899, strand: -1 }, 80884392, 200, 100).should.eql({ start: 11080, end: 15099, effUp: 200 });
});

test('gene mode clamps flanks at the region ends and shifts features by the effective upstream flank', async function () {
  const docs = {
    GP: { _id: 'GP', system_name: 'synthetic', location: { region: 'c', start: 50, end: 400, strand: 1 },
      gene_structure: { exons: [{ id: 'x1', start: 1, end: 351 }], transcripts: [{ id: 'GP.1', length: 351, exons: ['x1'] }] } },
    GM: { _id: 'GM', system_name: 'synthetic', location: { region: 'c', start: 600, end: 950, strand: -1 },
      gene_structure: { exons: [{ id: 'y1', start: 1, end: 351 }], transcripts: [{ id: 'GM.1', length: 351, exons: ['y1'] }] } }
  };
  const p = await template.buildTemplate({ mode: 'gene', gene_id: 'GP', flank_up: 200, flank_down: 700 }, syntheticDeps(docs));
  [p.start, p.end, p.length].should.eql([1, 1000, 1000]);
  p.features.gene.should.eql({ start: 50, end: 400 });
  p.seq.should.equal(C_SEQ);
  const m = await template.buildTemplate({ mode: 'gene', gene_id: 'GM', flank_up: 200, flank_down: 700 }, syntheticDeps(docs));
  [m.start, m.end, m.length].should.eql([1, 1000, 1000]);
  m.features.gene.should.eql({ start: 51, end: 401 });
  m.seq.should.equal(revcomp(C_SEQ));
  m.mapper.toGenomic(51).should.equal(950);
});

test('gene mode: TEMPLATE_TOO_LONG with a hint when the flanked gene exceeds max_template_length', async function () {
  const deps = stubs.sorghumDeps({ cfg: stubs.cfg({ design: { max_template_length: 4000 } }) });
  const err = await rejectsWith(template.buildTemplate({ mode: 'gene', gene_id: G200, flank_up: 200, flank_down: 100 }, deps),
    400, 'TEMPLATE_TOO_LONG');
  err.details.should.have.properties({ length: 4020, max: 4000 });
  err.details.should.have.property('hint');
  deps.sequence.calls.length.should.equal(0);
});

test('gene lookup errors: UNKNOWN_GENE, SYSTEM_NAME_MISMATCH, NO_SEQUENCE, UNKNOWN_REGION', async function () {
  await rejectsWith(template.buildTemplate({ mode: 'gene', gene_id: 'NOPE' }, stubs.sorghumDeps()), 404, 'UNKNOWN_GENE');
  const mm = await rejectsWith(template.buildTemplate({ mode: 'transcript', gene_id: G200, system_name: 'sorghum_rio' }, stubs.sorghumDeps()),
    400, 'SYSTEM_NAME_MISMATCH');
  mm.details.should.have.properties({ system_name: 'sorghum_rio', gene_system_name: 'sorghum_bicolor' });
  const noSeq = stubs.sorghumDeps({ resolve: async function () { return stubs.resolvedStub({ fasta: { dna: null, dna_sm: null } }); } });
  await rejectsWith(template.buildTemplate({ mode: 'gene', gene_id: G200 }, noSeq), 422, 'NO_SEQUENCE');
  const doc = stubs.geneDoc(G200);
  doc.location.region = 'constructor';
  await rejectsWith(template.buildTemplate({ mode: 'gene', gene_id: G200 },
    stubs.sorghumDeps({ findGene: async function () { return doc; } })), 404, 'UNKNOWN_REGION');
});

test('assembly warnings are carried on the template', async function () {
  const deps = stubs.sorghumDeps({
    resolve: async function () {
      return stubs.resolvedStub({ warnings: [{ code: 'ASSEMBLY_MISMATCH', message: 'only 1 of 11 map regions match' }] });
    }
  });
  const t = await template.buildTemplate({ mode: 'gene', gene_id: G200 }, deps);
  t.warnings.should.eql([{ code: 'ASSEMBLY_MISMATCH', message: 'only 1 of 11 map regions match' }]);
});

test('findGene: string-equality lookup through genes.findOne; non-strings are rejected before mongo; no collection is 503', async function () {
  const calls = [];
  const mongo = {
    genes: {
      mongoCollection: async function () {
        return { findOne: async function (q, opts) { calls.push([q, opts]); return null; } };
      }
    }
  };
  await rejectsWith(template.findGene({ $ne: 'x' }, { mongo: mongo }), 400, 'INVALID_REQUEST');
  await rejectsWith(template.findGene('', { mongo: mongo }), 400, 'INVALID_REQUEST');
  calls.length.should.equal(0);
  should(await template.findGene('SORBI_X', { mongo: mongo })).be.null();
  calls.should.eql([[{ _id: 'SORBI_X' }, { fields: { _id: 1, name: 1, system_name: 1, taxon_id: 1, location: 1, gene_structure: 1 } }]]);
  const down = { genes: { mongoCollection: async function () { return undefined; } } };
  const err = await rejectsWith(template.findGene('SORBI_X', { mongo: down, log: stubs.silentLog }), 503, 'MONGO_UNAVAILABLE');
  err.details.retry_after_s.should.be.above(0);
});

// ---- region mode ----------------------------------------------------------------------------------

test('region mode: + and - strand, genomic mapper, id and features', async function () {
  const plus = await template.buildTemplate({ mode: 'region', system_name: 'synthetic', region: { region: 'c', start: 11, end: 110, strand: 1 } },
    syntheticDeps({}));
  plus.seq.should.equal(plusC(11, 110));
  plus.features.should.eql({});
  plus.id.should.equal('synthetic_c_11-110');
  plus.mapper.toGenomic(1).should.equal(11);
  const minus = await template.buildTemplate({ mode: 'region', system_name: 'synthetic', region: { region: 'c', start: 11, end: 110, strand: -1 } },
    syntheticDeps({}));
  minus.seq.should.equal(revcomp(plusC(11, 110)));
  minus.mapper.primer('left', 1, 20).should.eql({ region: 'c', start: 91, end: 110, strand: -1, blocks: [{ start: 91, end: 110 }] });
  [minus.region, minus.start, minus.end, minus.strand, minus.length].should.eql(['c', 11, 110, -1, 100]);
});

test('region mode errors: REGION_OUT_OF_BOUNDS, UNKNOWN_REGION, TEMPLATE_TOO_LONG', async function () {
  await rejectsWith(template.buildTemplate({ mode: 'region', system_name: 'synthetic', region: { region: 'c', start: 950, end: 1001 } },
    syntheticDeps({})), 400, 'REGION_OUT_OF_BOUNDS');
  await rejectsWith(template.buildTemplate({ mode: 'region', system_name: 'synthetic', region: { region: 'c', start: 50, end: 40 } },
    syntheticDeps({})), 400, 'REGION_OUT_OF_BOUNDS');
  await rejectsWith(template.buildTemplate({ mode: 'region', system_name: 'synthetic', region: { region: 'zz', start: 1, end: 40 } },
    syntheticDeps({})), 404, 'UNKNOWN_REGION');
  const deps = syntheticDeps({}, { cfg: stubs.cfg({ design: { max_template_length: 500 } }) });
  await rejectsWith(template.buildTemplate({ mode: 'region', system_name: 'synthetic', region: { region: 'c', start: 1, end: 501 } }, deps),
    400, 'TEMPLATE_TOO_LONG');
  deps.sequence.calls.length.should.equal(0);
});

// ---- sequence mode --------------------------------------------------------------------------------

test('cleanSequence: FASTA headers, whitespace and digits removed; IUPAC codes become N; case kept only on request', function () {
  const raw = '>amp1 some description\r\n 1 acgtRYKM swbd\n61 hvNNacgtac gtacgtac\n>second\nACGT\n';
  const r = template.cleanSequence(raw, {});
  r.seq.should.equal('ACGTNNNNNNNNNNNNACGTACGTACGTACACGT');
  r.iupac.should.equal(10);
  r.name.should.equal('amp1');
  r.records.should.equal(2);
  const kept = template.cleanSequence(raw, { keepCase: true });
  kept.seq.should.equal('acgtNNNNnnnnnnNNacgtacgtacgtacACGT');
  kept.iupac.should.equal(10);
});

test('cleanSequence counts the FASTA records that carry sequence', function () {
  template.cleanSequence('ACGTACGTACGTACGTACGTACGT').records.should.equal(1);
  template.cleanSequence('>only\nACGTACGTACGTACGTACGTACGT\n').records.should.equal(1);
  // headers with no sequence after them (or only digits and whitespace) do not count
  template.cleanSequence('>a\n>b\n\nACGTACGTACGTACGTACGTACGT\n>trailing\n\n').records.should.equal(1);
  const three = template.cleanSequence('>a\nACGTACGTACGT\n>b\n  12 \n>c\n 1 ACGTACGTACGT\n>d\nacgt\n');
  three.records.should.equal(3);
  three.seq.should.equal('ACGTACGTACGTACGTACGTACGTACGT');
  three.name.should.equal('a');
});

test('sequence mode: a multi-record FASTA is joined (not rejected) with warning MULTIPLE_RECORDS', async function () {
  const a = C_SEQ.slice(0, 30).toUpperCase();
  const b = C_SEQ.slice(500, 530).toUpperCase();
  const t = await template.buildTemplate({ mode: 'sequence', sequence: '>recA\n' + a + '\n>recB\n' + b + '\n' }, syntheticDeps({}));
  t.seq.should.equal(a + b);
  t.length.should.equal(60);
  t.id.should.equal('recA');
  t.warnings.should.eql([{ code: 'MULTIPLE_RECORDS', message: '2 FASTA records were joined into one template; primers may span the joins' }]);
  const withIupac = await template.buildTemplate({ mode: 'sequence', sequence: '>x\n' + a + 'R\n>y\n' + b + '\n>z\n' + a + '\n' }, syntheticDeps({}));
  withIupac.warnings.map(function (w) { return w.code; }).should.eql(['MULTIPLE_RECORDS', 'IUPAC_CONVERTED']);
  withIupac.warnings[0].message.should.equal('3 FASTA records were joined into one template; primers may span the joins');
  (await template.buildTemplate({ mode: 'sequence', sequence: '>recA\n' + a + '\n' }, syntheticDeps({}))).warnings.should.eql([]);
});

test('cleanSequence errors: INVALID_SEQUENCE for letters outside IUPAC DNA, too short or empty; TEMPLATE_TOO_LONG', function () {
  let err = null;
  try { template.cleanSequence('ACGTACGTACGTACGTACGTACGU'); } catch (e) { err = e; }
  err.code.should.equal('INVALID_SEQUENCE');
  err.details.should.eql({ position: 24, character: 'U' });
  (function () { template.cleanSequence('ACGT-ACGTACGTACGTACGTACGT'); }).should.throw({ code: 'INVALID_SEQUENCE' });
  (function () { template.cleanSequence('ACGTACGTACGTACGTACG'); }).should.throw({ code: 'INVALID_SEQUENCE' });
  (function () { template.cleanSequence('>only a header\n\n'); }).should.throw({ code: 'INVALID_SEQUENCE' });
  (function () { template.cleanSequence(42); }).should.throw({ code: 'INVALID_SEQUENCE' });
  (function () { template.cleanSequence('A'.repeat(101), { maxLength: 100 }); }).should.throw({ code: 'TEMPLATE_TOO_LONG' });
});

test('sequence mode: IUPAC_CONVERTED warning, no genomic mapper, optional system_name resolved', async function () {
  const t = await template.buildTemplate({ mode: 'sequence', sequence: '>x\nACGTRYACGTACGTACGTACGTACGT' }, syntheticDeps({}));
  t.seq.should.equal('ACGTNNACGTACGTACGTACGTACGT');
  t.warnings.should.eql([{ code: 'IUPAC_CONVERTED', message: '2 IUPAC ambiguity codes were converted to N' }]);
  should(t.mapper).be.null();
  should(t.resolved).be.null();
  [t.region, t.start, t.end, t.strand].should.eql([null, null, null, null]);
  t.id.should.equal('x');
  const lower = await template.buildTemplate({ mode: 'sequence', sequence: 'acgtACGTACGTACGTACGTacgt', avoid_repeats: true, system_name: 'synthetic' },
    syntheticDeps({}));
  lower.seq.should.equal('acgtACGTACGTACGTACGTacgt');
  lower.system_name.should.equal('synthetic');
  lower.resolved.system_name.should.equal('synthetic');
  const upper = await template.buildTemplate({ mode: 'sequence', sequence: 'acgtACGTACGTACGTACGTacgt' }, syntheticDeps({}));
  upper.seq.should.equal('ACGTACGTACGTACGTACGTACGT');
  const unknown = syntheticDeps({}, {
    resolve: async function (name) { throw new (require('../../../api/helpers/primers/errors').PrimerHttpError)(404, 'UNKNOWN_GENOME', 'x', { system_name: name }); }
  });
  await rejectsWith(template.buildTemplate({ mode: 'sequence', sequence: 'ACGTACGTACGTACGTACGTACGT', system_name: 'nope' }, unknown), 404, 'UNKNOWN_GENOME');
});

test('an unknown mode is 400 INVALID_REQUEST', async function () {
  await rejectsWith(template.buildTemplate({ mode: 'protein' }, {}), 400, 'INVALID_REQUEST');
});

// ---- real sequence.js on a bgzip FASTA written in os.tmpdir() --------------------------------------

test('transcript and region templates through the real sequence.js reader (bgzip FASTA fixture)', async function (t) {
  const root = fx.makeRoot('design-tpl');
  t.after(function () { fx.removeRoot(root); });
  const seq = fx.randomSeq(5000, 11);
  const file = fx.bgzipFasta(path.join(root, 'synthetic'), 'Syn', 'dna', [{ name: 'c', seq: seq }, { name: 'd', seq: fx.randomSeq(300, 3) }]);
  fs.existsSync(file + '.gzi').should.be.true();
  const doc = {
    _id: 'RG', system_name: 'synthetic', location: { region: 'c', start: 1001, end: 4000, strand: -1 },
    gene_structure: {
      exons: [{ id: 'r1', start: 1, end: 120 }, { id: 'r2', start: 1001, end: 1180 }, { id: 'r3', start: 2801, end: 3000 }],
      transcripts: [{ id: 'RG.1', length: 500, exons: ['r1', 'r2', 'r3'], exon_junctions: [120, 300] }],
      canonical_transcript: 'RG.1'
    }
  };
  const deps = {
    cfg: stubs.cfg(),
    log: stubs.silentLog,
    findGene: async function () { return JSON.parse(JSON.stringify(doc)); },
    resolve: async function () { return stubs.resolvedStub({ system_name: 'synthetic', fasta: { dna: file, dna_sm: null } }); }
  };
  const tr = await template.buildTemplate({ mode: 'transcript', gene_id: 'RG' }, deps);
  const plus = function (s, e) { return seq.slice(s - 1, e); };
  tr.seq.should.equal(revcomp(plus(3881, 4000)) + revcomp(plus(2821, 3000)) + revcomp(plus(1001, 1200)));
  tr.junctions.should.eql([120, 300]);
  const rg = await template.buildTemplate({ mode: 'region', system_name: 'synthetic', region: { region: 'c', start: 4901, end: 5000, strand: -1 } }, deps);
  rg.seq.should.equal(revcomp(plus(4901, 5000)));
  await rejectsWith(template.buildTemplate({ mode: 'region', system_name: 'synthetic', region: { region: 'c', start: 4901, end: 5001 } }, deps),
    400, 'REGION_OUT_OF_BOUNDS');
});
