'use strict';

const test = require('node:test');
const should = require('should');

const coords = require('../../../api/helpers/primers/coords');

// Gene docs from GET /sorghum_v11/genes (sorghum_bicolor), canonical transcripts.
// Genomic values below were verified against the genome through fastaIdx (:8888).
const G200 = {
  loc: { region: '1', start: 11180, end: 14899, strand: -1 },
  // SORBI_3001G000200.1 exons in transcript order (gene-relative, 1-based)
  exons: [
    { id: 'EER93047-1', start: 1, end: 397 }, { id: 'EER93047-2', start: 484, end: 579 },
    { id: 'EER93047-3', start: 780, end: 883 }, { id: 'EER93047-4', start: 992, end: 1263 },
    { id: 'EER93047-5', start: 1500, end: 1622 }, { id: 'EER93047-6', start: 1718, end: 1806 },
    { id: 'EER93047-7', start: 1900, end: 2005 }, { id: 'EER93047-8', start: 2114, end: 2211 },
    { id: 'EER93047-9', start: 2748, end: 3008 }, { id: 'EER93047-10', start: 3168, end: 3251 },
    { id: 'EER93047-11', start: 3369, end: 3720 }
  ],
  length: 1982,
  exon_junctions: [397, 493, 597, 869, 992, 1081, 1187, 1285, 1546, 1630]
};

const G700 = {
  loc: { region: '1', start: 53781, end: 63305, strand: 1 },
  // SORBI_3001G000700.3 (canonical) exons in transcript order
  exons: [
    { id: 'OQU90575-1', start: 6310, end: 6402 }, { id: 'OQU90575-2', start: 6527, end: 6671 },
    { id: 'OQU90575-3', start: 6946, end: 7006 }, { id: 'OQU90577-2', start: 7475, end: 7528 },
    { id: 'OQU90577-3', start: 8288, end: 8367 }, { id: 'OQU90577-4', start: 8791, end: 8996 },
    { id: 'OQU90577-5', start: 9069, end: 9525 }
  ],
  length: 1096,
  exon_junctions: [93, 238, 299, 353, 433, 639]
};

test('footprints: left [pos, pos+len-1], right [pos-len+1, pos], product rpos-lpos+1', function () {
  coords.leftFootprint(39528, 20).should.eql({ start: 39528, end: 39547 });
  coords.rightFootprint(40043, 20).should.eql({ start: 40024, end: 40043 });
  coords.footprint('left', 856, 20).should.eql({ start: 856, end: 875 });
  coords.footprint('right', 979, 20).should.eql({ start: 960, end: 979 });
  coords.productSize(39528, 40043).should.equal(516);
  coords.intervalToRange([499, 50]).should.eql({ start: 499, end: 548 });
  coords.rangeToInterval(499, 548).should.eql([499, 50]);
  (function () { coords.footprint('middle', 1, 2); }).should.throw(/side/);
  (function () { coords.leftFootprint(1.5, 2); }).should.throw(/integer/);
});

test('primer genomic strand: left = template strand, right = opposite', function () {
  coords.primerStrand('left', -1).should.equal(-1);
  coords.primerStrand('right', -1).should.equal(1);
  coords.primerStrand('left', 1).should.equal(1);
  coords.primerStrand('right', 1).should.equal(-1);
});

test('SORBI_3001G000200 (-): segments and derived junctions', function () {
  const segs = coords.buildSegments(G200.loc, G200.exons);
  segs[0].should.eql({ id: 'EER93047-1', t_start: 1, t_end: 397, g_start: 14503, g_end: 14899 });
  segs[3].should.eql({ id: 'EER93047-4', t_start: 598, t_end: 869, g_start: 13637, g_end: 13908 });
  segs[4].should.eql({ id: 'EER93047-5', t_start: 870, t_end: 992, g_start: 13278, g_end: 13400 });
  segs[10].should.eql({ id: 'EER93047-11', t_start: 1631, t_end: 1982, g_start: 11180, g_end: 11531 });
  segs[segs.length - 1].t_end.should.equal(G200.length);
  coords.junctionsFromSegments(segs).should.eql(G200.exon_junctions);
  coords.geneRelativeToGenomic(G200.loc, 1).should.equal(14899);
  coords.geneRelativeToGenomic(G200.loc, 3720).should.equal(11180);
  coords.exonGenomicRange(G200.loc, { start: 992, end: 1263 }).should.eql({ start: 13637, end: 13908 });
});

test('SORBI_3001G000200 (-): cDNA 856-875 -> 1:13395-13400 + 13637-13650 (-1); 960-979 -> 13291-13310 (+1)', function () {
  const m = coords.splicedMapper({ region: '1', strand: -1, segments: coords.buildSegments(G200.loc, G200.exons) });
  m.length.should.equal(1982);
  m.junctions.should.eql(G200.exon_junctions);
  m.primer('left', 856, 875).should.eql({
    region: '1', start: 13395, end: 13650, strand: -1,
    blocks: [{ start: 13395, end: 13400 }, { start: 13637, end: 13650 }]
  });
  m.primer('right', 960, 979).should.eql({
    region: '1', start: 13291, end: 13310, strand: 1, blocks: [{ start: 13291, end: 13310 }]
  });
  m.product(856, 979).should.eql({ region: '1', start: 13291, end: 13650, strand: -1 });
  m.genomicSize(856, 979).should.equal(360);
  m.toGenomic(856).should.equal(13650);
  m.toGenomic(870).should.equal(13400);
  m.toTemplate(13650).should.equal(856);
  m.toTemplate(13400).should.equal(870);
  should(m.toTemplate(13500)).equal(null); // intron
  (function () { m.blocks(1980, 1983); }).should.throw(RangeError);
});

test('SORBI_3001G000200 gene-mode template 1:11080-15099(-): t=499 -> 14601', function () {
  // flanks up 200 / down 100 on a - strand gene: gStart = start-100, gEnd = end+200, effUp = 200
  const m = coords.genomicMapper({ region: '1', start: 11080, end: 15099, strand: -1 });
  m.length.should.equal(4020);
  m.toGenomic(499).should.equal(14601);
  m.toTemplate(14601).should.equal(499);
  m.toGenomic(201).should.equal(14899); // exon 1 start (gene start in transcription order)
  m.toGenomic(597).should.equal(14503); // exon 1 end
  should(m.toTemplate(11079)).equal(null);
  m.primer('left', 480, 499).should.eql({ region: '1', start: 14601, end: 14620, strand: -1, blocks: [{ start: 14601, end: 14620 }] });
  m.primer('right', 499, 518).should.eql({ region: '1', start: 14582, end: 14601, strand: 1, blocks: [{ start: 14582, end: 14601 }] });
  m.product(480, 518).should.eql({ region: '1', start: 14582, end: 14620, strand: -1 });
  (function () { m.toGenomic(4021); }).should.throw(RangeError);
});

test('SORBI_3001G000700 (+) equivalents', function () {
  const segs = coords.buildSegments(G700.loc, G700.exons);
  segs[0].should.eql({ id: 'OQU90575-1', t_start: 1, t_end: 93, g_start: 60090, g_end: 60182 });
  segs[1].should.eql({ id: 'OQU90575-2', t_start: 94, t_end: 238, g_start: 60307, g_end: 60451 });
  segs[6].should.eql({ id: 'OQU90577-5', t_start: 640, t_end: 1096, g_start: 62849, g_end: 63305 });
  coords.junctionsFromSegments(segs).should.eql(G700.exon_junctions);
  const m = coords.splicedMapper({ region: '1', strand: 1, segments: segs });
  m.length.should.equal(G700.length);
  // verified: CTTTTGCATTTACAGTAGGA = 1:60169-60182(+) . 1:60307-60312(+)
  m.primer('left', 80, 99).should.eql({
    region: '1', start: 60169, end: 60312, strand: 1,
    blocks: [{ start: 60169, end: 60182 }, { start: 60307, end: 60312 }]
  });
  // verified: GTGCTGAGACTTGAACAGAT = 1:60413-60432(-)
  m.primer('right', 200, 219).should.eql({ region: '1', start: 60413, end: 60432, strand: -1, blocks: [{ start: 60413, end: 60432 }] });
  m.product(80, 219).should.eql({ region: '1', start: 60169, end: 60432, strand: 1 });
  m.genomicSize(80, 219).should.equal(264);
  m.toTemplate(60182).should.equal(93);
  m.toTemplate(60307).should.equal(94);
  coords.geneRelativeToGenomic(G700.loc, 1).should.equal(53781);
  // gene mode, flanks 200/100 on a + strand gene: template 53581-63405, t=201 is the gene start
  const g = coords.genomicMapper({ region: '1', start: 53581, end: 63405, strand: 1 });
  g.toGenomic(201).should.equal(53781);
  g.primer('right', 201, 220).should.eql({ region: '1', start: 53781, end: 53800, strand: -1, blocks: [{ start: 53781, end: 53800 }] });
});

test('junction overlap formulas (A.3)', function () {
  const J = G200.exon_junctions;
  // left [a,b] spans j when j-a+1 >= min5 and b-j >= min3
  coords.junctionOverlap('left', 856, 875, J, { min5: 7, min3: 4 }).should.eql({ position: 869, overlap_5p: 14, overlap_3p: 6 });
  should(coords.junctionOverlap('left', 856, 875, J, { min5: 15, min3: 4 })).equal(null);
  should(coords.junctionOverlap('left', 856, 875, J, { min5: 7, min3: 7 })).equal(null);
  should(coords.junctionOverlap('left', 863, 872, J)).equal(null); // 7 / 3
  coords.junctionOverlap('left', 863, 873, J).should.eql({ position: 869, overlap_5p: 7, overlap_3p: 4 });
  // right [a,b] (5' end at b) spans j when j-a+1 >= min3 and b-j >= min5
  should(coords.junctionOverlap('right', 960, 979, J)).equal(null);
  coords.junctionOverlap('right', 980, 999, J).should.eql({ position: 992, overlap_5p: 7, overlap_3p: 13 });
  coords.junctionOverlap('right', 989, 999, J).should.eql({ position: 992, overlap_5p: 7, overlap_3p: 4 });
  should(coords.junctionOverlap('right', 990, 999, J)).equal(null); // 3' overlap 3 < 4
  should(coords.junctionOverlap('right', 980, 998, J)).equal(null); // 5' overlap 6 < 7
  // the same span read as a left primer: 5' overlap is j-a+1
  should(coords.junctionOverlap('left', 989, 999, J)).equal(null);
  coords.spansJunction('left', 856, 875, 869, 14, 6).should.equal(true);
  coords.spansJunction('left', 856, 875, 869, 15, 6).should.equal(false);
  coords.spansJunction('right', 856, 875, 869, 6, 14).should.equal(true);
  coords.spansJunction('right', 856, 875, 869, 7, 14).should.equal(false);
  should(coords.junctionOverlap('left', 856, 875, [])).equal(null);
});

test('mask runs: merge, flags, lowercase, fraction, overlap', function () {
  coords.mergeRuns([[10, 5], [1, 3], [4, 2], [20, 1], [14, 3], [30, 0]]).should.eql([[1, 5], [10, 7], [20, 1]]);
  coords.mergeRuns([[10, 5], [1, 3], [4, 2], [20, 1], [14, 3]], { length: 18 }).should.eql([[1, 5], [10, 7]]);
  coords.mergeRuns([[-3, 5]], { length: 10 }).should.eql([[1, 1]]);
  coords.mergeRuns([]).should.eql([]);
  coords.runsFromFlags([0, 1, 1, 0, 1]).should.eql([[2, 2], [5, 1]]);
  coords.runsFromFlags(Uint8Array.from([1, 1, 1])).should.eql([[1, 3]]);
  coords.runsFromFlags([]).should.eql([]);
  coords.lowercaseRuns('ACgtNNaaT').should.eql([[3, 2], [7, 2]]);
  coords.lowercaseRuns('ACGT').should.eql([]);
  coords.maskedBases([[1, 5], [3, 5]]).should.equal(7);
  coords.maskedFraction([[1, 5], [3, 5]], 10).should.equal(0.7);
  coords.maskedFraction([], 0).should.equal(0);
  coords.overlapsRuns(5, 9, [[10, 2]]).should.equal(false);
  coords.overlapsRuns(5, 10, [[10, 2]]).should.equal(true);
  coords.overlapsRuns(12, 20, [[10, 2]]).should.equal(false);
  coords.overlapsRuns(11, 20, [[10, 2]]).should.equal(true);
});

test('spliced mapper rejects inconsistent segments', function () {
  (function () {
    coords.splicedMapper({ region: '1', strand: 1, segments: [{ t_start: 1, t_end: 10, g_start: 100, g_end: 109 }, { t_start: 12, t_end: 20, g_start: 200, g_end: 208 }] });
  }).should.throw(/contiguous/);
  (function () {
    coords.splicedMapper({ region: '1', strand: 1, segments: [{ t_start: 1, t_end: 10, g_start: 100, g_end: 120 }] });
  }).should.throw(/contiguous/);
  (function () { coords.splicedMapper({ region: '1', strand: 0, segments: [] }); }).should.throw(/strand/);
});
