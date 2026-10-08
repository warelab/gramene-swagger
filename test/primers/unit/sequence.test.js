'use strict';

require('../../../api/helpers/primers/node_compat');

const test = require('node:test');
const should = require('should');
const fs = require('fs');
const os = require('os');
const path = require('path');
const zlib = require('zlib');
const { execFileSync } = require('child_process');

const sequence = require('../../../api/helpers/primers/sequence');
const { PrimerHttpError } = require('../../../api/helpers/primers/errors');

const REALDATA = process.env.PRIMERS_REALDATA === '1';
const SAMTOOLS = ['/usr/local/bin/samtools', '/usr/bin/samtools'].find(function (p) { return fs.existsSync(p); });

// ---- a tiny hermetic bgzip FASTA (+ .fai + .gzi), written in pure node ---------------------

const BGZF_EOF = Buffer.from('1f8b08040000000000ff0600424302001b0003000000000000000000', 'hex');

function bgzfBlock(data) {
  const cdata = zlib.deflateRawSync(data, { level: 6 });
  const header = Buffer.from([0x1f, 0x8b, 8, 4, 0, 0, 0, 0, 0, 0xff, 6, 0, 0x42, 0x43, 2, 0, 0, 0]);
  const total = header.length + cdata.length + 8;
  header.writeUInt16LE(total - 1, 16);
  const trailer = Buffer.alloc(8);
  trailer.writeUInt32LE(zlib.crc32(data) >>> 0, 0);
  trailer.writeUInt32LE(data.length, 4);
  return Buffer.concat([header, cdata, trailer]);
}

// records: [{name, seq}], lineWidth; blockSize = uncompressed bytes per BGZF block
function writeBgzipFasta(file, records, lineWidth, blockSize) {
  let text = '';
  const fai = [];
  records.forEach(function (r) {
    text += '>' + r.name + ' test record\n';
    const offset = Buffer.byteLength(text);
    for (let i = 0; i < r.seq.length; i += lineWidth) text += r.seq.slice(i, i + lineWidth) + '\n';
    fai.push([r.name, r.seq.length, offset, lineWidth, lineWidth + 1].join('\t'));
  });
  const raw = Buffer.from(text, 'latin1');
  const blocks = [];
  const gzi = [];
  let coff = 0;
  for (let u = 0; u < raw.length; u += blockSize) {
    if (u > 0) gzi.push([coff, u]);
    const b = bgzfBlock(raw.subarray(u, Math.min(raw.length, u + blockSize)));
    blocks.push(b);
    coff += b.length;
  }
  fs.writeFileSync(file, Buffer.concat(blocks.concat([BGZF_EOF])));
  fs.writeFileSync(file + '.fai', fai.join('\n') + '\n');
  const idx = Buffer.alloc(8 + 16 * gzi.length);
  idx.writeBigUInt64LE(BigInt(gzi.length), 0);
  gzi.forEach(function (e, i) {
    idx.writeBigUInt64LE(BigInt(e[0]), 8 + 16 * i);
    idx.writeBigUInt64LE(BigInt(e[1]), 16 + 16 * i);
  });
  fs.writeFileSync(file + '.gzi', idx);
}

function randomSeq(n, alphabet, seed) {
  let x = seed;
  let s = '';
  for (let i = 0; i < n; i++) {
    x = (x * 1103515245 + 12345) & 0x7fffffff;
    s += alphabet[x % alphabet.length];
  }
  return s;
}

const TMP = fs.mkdtempSync(path.join(os.tmpdir(), 'primers-sequence-test-'));
const FA = path.join(TMP, 'Test_genome.v1.dna_sm.toplevel.fa.gz');
const CHR_A = randomSeq(2503, 'ACGTACGTacgtNRY', 7);
const CHR_B = 'ACGTRYKMBVDHSWNacgtrykmbvdhswn' + randomSeq(20, 'ACGT', 3);
writeBgzipFasta(FA, [{ name: 'chrA', seq: CHR_A }, { name: 'chrB', seq: CHR_B }], 60, 1000);

async function rejection(p) {
  try {
    await p;
  } catch (e) {
    return e;
  }
  throw new Error('expected a rejection');
}

test.after(function () {
  sequence._clearCache();
  fs.rmSync(TMP, { recursive: true, force: true });
});

test('revcomp is IUPAC-aware and case-preserving', function () {
  sequence.revcomp('ACGTRYKMBVDHSWNacgtrykmbvdhswn').should.equal('nwsdhbvkmryacgtNWSDHBVKMRYACGT');
  sequence.complement('ACGTRYKMBVDHSWNacgtrykmbvdhswn').should.equal('TGCAYRMKVBHDSWNtgcayrmkvbhdswn');
  sequence.revcomp('').should.equal('');
  sequence.revcomp('A-c.').should.equal('.g-T');
  sequence.revcomp('Uu').should.equal('aA');
  sequence.revcomp(sequence.revcomp(CHR_A)).should.equal(CHR_A);
  sequence.revcomp('GGACAGCTCCACAACATATCAG').should.equal('CTGATATGTTGTGGAGCTGTCC');
  (function () { sequence.revcomp('ACGT☃'); }).should.throw(TypeError);
  (function () { sequence.revcomp(null); }).should.throw(TypeError);
});

test('hermetic fixture is a valid bgzip/faidx file (samtools cross-check)', { skip: !SAMTOOLS && 'samtools not installed' }, function () {
  const out = execFileSync(SAMTOOLS, ['faidx', FA, 'chrA:901-1100', 'chrB:1-30'], { encoding: 'utf8' });
  const recs = out.split('>').filter(Boolean).map(function (r) { return r.split('\n').slice(1).join(''); });
  recs[0].should.equal(CHR_A.slice(900, 1100));
  recs[1].should.equal(CHR_B.slice(0, 30));
});

test('fetch returns 1-based inclusive sequence, case preserved, across BGZF blocks', async function () {
  (await sequence.regionLength(FA, 'chrA')).should.equal(2503);
  (await sequence.regionLength(FA, 'chrB')).should.equal(CHR_B.length);
  should(await sequence.regionLength(FA, 'chrC')).be.undefined();
  should(await sequence.regionLength(FA, 'constructor')).be.undefined();
  should(await sequence.regionLength(FA, '__proto__')).be.undefined();
  should(await sequence.regionLength(FA, '')).be.undefined();

  (await sequence.fetch(FA, 'chrA', 1, 60, 1)).should.equal(CHR_A.slice(0, 60));
  (await sequence.fetch(FA, 'chrA', 900, 2100, 1)).should.equal(CHR_A.slice(899, 2100));
  (await sequence.fetch(FA, 'chrA', 2503, 2503, 1)).should.equal(CHR_A.slice(2502));
  (await sequence.fetch(FA, 'chrA', 1, 2503)).should.equal(CHR_A); // strand defaults to +1
  (await sequence.fetch(FA, 'chrA', 900, 2100, -1)).should.equal(sequence.revcomp(CHR_A.slice(899, 2100)));
  (await sequence.fetch(FA, 'chrB', 1, 30, '-1')).should.equal('nwsdhbvkmryacgtNWSDHBVKMRYACGT');
  (await sequence.fetch(FA, 'chrB', 16, 30, '+')).should.equal('acgtrykmbvdhswn');
});

test('fetch errors: UNKNOWN_REGION 404, REGION_OUT_OF_BOUNDS 400, TEMPLATE_TOO_LONG 400, INVALID_REQUEST 400', async function () {
  let e = await rejection(sequence.fetch(FA, 'chrZ', 1, 10, 1));
  e.should.be.instanceOf(PrimerHttpError);
  e.status.should.equal(404);
  e.code.should.equal('UNKNOWN_REGION');
  e.details.should.eql({ region: 'chrZ' });

  e = await rejection(sequence.fetch(FA, 'constructor', 1, 1, 1));
  e.code.should.equal('UNKNOWN_REGION');

  e = await rejection(sequence.fetch(FA, 'chrA', 2500, 2504, 1));
  e.status.should.equal(400);
  e.code.should.equal('REGION_OUT_OF_BOUNDS');
  e.details.should.eql({ region: 'chrA', start: 2500, end: 2504, length: 2503 });

  for (const bad of [[0, 10], [10, 9], [1.5, 10], ['1', 10], [1, NaN]]) {
    e = await rejection(sequence.fetch(FA, 'chrA', bad[0], bad[1], 1));
    e.code.should.equal('REGION_OUT_OF_BOUNDS', JSON.stringify(bad));
  }

  e = await rejection(sequence.fetch(FA, 'chrA', 1, 100, 1, { maxLength: 99 }));
  e.code.should.equal('TEMPLATE_TOO_LONG');
  e.details.should.eql({ length: 100, max: 99 });

  // the 2 Mb internal cap comes from config and is checked before any I/O
  e = await rejection(sequence.fetch(path.join(TMP, 'absent.fa.gz'), '1', 1, 2000001, 1));
  e.code.should.equal('TEMPLATE_TOO_LONG');
  e.details.max.should.equal(2000000);

  e = await rejection(sequence.fetch(FA, 'chrA', 1, 10, 2));
  e.status.should.equal(400);
  e.code.should.equal('INVALID_REQUEST');
});

test('missing FASTA/index files give 422 NO_SEQUENCE without paths, and the handle is evicted', async function () {
  const absent = path.join(TMP, 'nope', 'Missing.dna.toplevel.fa.gz');
  const e = await rejection(sequence.fetch(absent, '1', 1, 10, 1));
  e.status.should.equal(422);
  e.code.should.equal('NO_SEQUENCE');
  JSON.stringify(e).should.not.match(/nope|Missing|primers-sequence-test/);
  sequence._hasHandle(absent).should.equal(false);

  const noGzi = path.join(TMP, 'NoGzi.fa.gz');
  fs.copyFileSync(FA, noGzi);
  fs.copyFileSync(FA + '.fai', noGzi + '.fai');
  const e2 = await rejection(sequence.fetch(noGzi, 'chrA', 1, 10, 1));
  e2.code.should.equal('NO_SEQUENCE');
  sequence._hasHandle(noGzi).should.equal(false);
  // once the index appears, the next request succeeds (no cached failure)
  fs.copyFileSync(FA + '.gzi', noGzi + '.gzi');
  (await sequence.fetch(noGzi, 'chrA', 1, 10, 1)).should.equal(CHR_A.slice(0, 10));
});

test('handle cache is an LRU of 32 FASTA files keyed by path', async function () {
  sequence._clearCache();
  sequence.MAX_HANDLES.should.equal(32);
  const links = [];
  for (let i = 0; i < 34; i++) {
    const p = path.join(TMP, 'lru' + i + '.fa.gz');
    ['', '.fai', '.gzi'].forEach(function (ext) { fs.symlinkSync(FA + ext, p + ext); });
    links.push(p);
  }
  for (let i = 0; i < 32; i++) (await sequence.fetch(links[i], 'chrB', 1, 4, 1)).should.equal('ACGT');
  sequence._cacheSize().should.equal(32);
  (await sequence.fetch(links[0], 'chrB', 1, 4, 1)).should.equal('ACGT'); // touch: most recently used
  (await sequence.fetch(links[32], 'chrB', 1, 4, 1)).should.equal('ACGT');
  sequence._cacheSize().should.equal(32);
  sequence._hasHandle(links[0]).should.equal(true);
  sequence._hasHandle(links[1]).should.equal(false); // least recently used was evicted
  (await sequence.fetch(links[33], 'chrB', 1, 4, -1)).should.equal('ACGT');
  sequence._hasHandle(links[2]).should.equal(false);
  sequence._cacheSize().should.equal(32);
  sequence._clearCache();
});

// ---- real data (opt-in): PRIMERS_REALDATA=1 --------------------------------------------------

const BICOLOR = '/scratch/olson/fasta/sorghum_bicolor/dna/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel.fa.gz';
const TX436_SM = '/scratch/olson/fasta/sorghum_tx436pac/dna/Sorghum_tx436pac.Sorghum_bicolor-Tx436-Reference-CSHL-USDA-1.0.dna_sm.toplevel.fa.gz';

async function fastaIdx(system, region, start, end, strand) {
  const url = 'http://localhost:8888/sequence/region/' + system + '/' + region + ':' + start + '..' + end + ':' + strand;
  const res = await fetch(url, { signal: AbortSignal.timeout(10000) });
  const body = await res.json();
  should.exist(body.seq, url + ' -> ' + JSON.stringify(body).slice(0, 200));
  return body.seq;
}

test('real data: sequence.fetch matches fastaIdx on + and - strands (sorghum_bicolor)', { skip: !REALDATA && 'set PRIMERS_REALDATA=1' }, async function () {
  const cases = [
    ['1', 13291, 13310, 1], ['1', 13637, 13650, -1], ['1', 13395, 13400, -1],
    ['4', 7400001, 7420000, 1], ['4', 7400001, 7420000, -1], ['1', 11080, 15099, -1]
  ];
  for (const c of cases) {
    const ours = await sequence.fetch(BICOLOR, c[0], c[1], c[2], c[3]);
    const theirs = await fastaIdx('sorghum_bicolor', c[0], c[1], c[2], c[3]);
    ours.length.should.equal(c[2] - c[1] + 1);
    ours.should.equal(theirs, c.join(':'));
  }
  (await sequence.fetch(BICOLOR, '1', 13637, 13650, -1) + await sequence.fetch(BICOLOR, '1', 13395, 13400, -1)).should.equal('ATTACATCAAATAGGCCTTG');
  (await sequence.fetch(BICOLOR, '1', 13291, 13310, 1)).should.equal('AACTTCTTTGTCGATCCATG');
});

test('real data: soft-masked dna_sm keeps lowercase on both strands (sorghum_tx436pac, vs samtools)', { skip: (!REALDATA && 'set PRIMERS_REALDATA=1') || (!SAMTOOLS && 'samtools not installed') }, async function () {
  const sam = execFileSync(SAMTOOLS, ['faidx', TX436_SM, '4:7550080-7550100'], { encoding: 'utf8' }).split('\n').slice(1).join('');
  sam.should.equal('AGGTACTGTGACtaaggatga');
  (await sequence.fetch(TX436_SM, '4', 7550080, 7550100, 1)).should.equal(sam);
  (await sequence.fetch(TX436_SM, '4', 7550080, 7550100, -1)).should.equal('tcatccttaGTCACAGTACCT');
  const len = await sequence.regionLength(TX436_SM, '4');
  len.should.equal(68123459);
  const e = await rejection(sequence.fetch(TX436_SM, '4', len - 5, len + 1, 1));
  e.code.should.equal('REGION_OUT_OF_BOUNDS');
  (await sequence.fetch(TX436_SM, '4', len - 9, len, 1)).length.should.equal(10);
  (await rejection(sequence.fetch(TX436_SM, 'chr4', 1, 10, 1))).code.should.equal('UNKNOWN_REGION');
});

test('real data: a 2 Mb fetch at the internal cap succeeds; one base more is refused', { skip: !REALDATA && 'set PRIMERS_REALDATA=1' }, async function () {
  const s = await sequence.fetch(BICOLOR, '1', 1000001, 3000000, -1);
  s.length.should.equal(2000000);
  s.slice(-20).should.equal(sequence.revcomp(await sequence.fetch(BICOLOR, '1', 1000001, 1000020, 1)));
  (await rejection(sequence.fetch(BICOLOR, '1', 1000001, 3000001, 1))).code.should.equal('TEMPLATE_TOO_LONG');
});
