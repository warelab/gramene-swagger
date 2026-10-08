'use strict';

require('../../../api/helpers/primers/node_compat');

const { test } = require('node:test');
const assert = require('node:assert');
const should = require('should');
const fs = require('fs');
const path = require('path');

const genomes = require('../../../api/helpers/primers/genomes');
const assemblies = require('../../../api/helpers/primers/assemblies');
const fx = require('../fixtures/catalog/fsfixture');

const FIXTURES = path.join(__dirname, '..', 'fixtures', 'catalog');
const SORGHUM_MAPS = JSON.parse(fs.readFileSync(path.join(FIXTURES, 'maps_sorghum_v11.json'), 'utf8'));
const V69_MAPS = JSON.parse(fs.readFileSync(path.join(FIXTURES, 'maps_v69_subset.json'), 'utf8'));
const RESULT_KEYS = ['blastdb', 'dir', 'display_name', 'fasta', 'fingerprint', 'map_id', 'num_sequences', 'prefix',
  'repeat_masking', 'system_name', 'taxon_id', 'total_bases', 'warnings'];
const SILENT = { warn: function () {}, error: function () {}, log: function () {} };
const T2020 = 1596067200; // 2020-07-30
const T2021 = 1637366400; // 2021-11-20

function tmpRoot(t) {
  const root = fx.makeRoot('asm');
  t.after(function () { fx.removeRoot(root); });
  return root;
}

function mapDoc(systemName, id, entries, taxonId) {
  return {
    _id: id, system_name: systemName, display_name: systemName, taxon_id: taxonId || 4558999, type: 'genome',
    regions: { names: entries.map(function (e) { return e[0]; }), lengths: entries.map(function (e) { return e[1]; }) }
  };
}

function depsFor(root, maps, extra) {
  return Object.assign({
    cfg: { fasta_root: root, assembly_overrides: {}, repeat_masking_overrides: {} },
    catalog: genomes.buildCatalog(maps, []),
    log: SILENT,
    cache: false
  }, extra || {});
}

function sorghumMap(systemName) {
  const m = SORGHUM_MAPS.find(function (x) { return x.system_name === systemName; });
  assert.ok(m, systemName);
  return m;
}

function regionsOf(map) {
  return map.regions.names.map(function (n, i) { return [n, map.regions.lengths[i]]; });
}

async function rejectsWith(promise, status, code) {
  let caught = null;
  await assert.rejects(promise, function (err) { caught = err; return true; });
  should(caught.status).equal(status);
  should(caught.code).equal(code);
  return caught;
}

function codes(result) {
  return result.warnings.map(function (w) { return w.code; });
}

function genomeSeqs(seed) {
  return [
    { name: '1', seq: fx.randomSeq(60000, seed + 1) },
    { name: '2', seq: fx.randomSeq(30000, seed + 2) },
    { name: '3', seq: fx.randomSeq(20000, seed + 3) },
    { name: 'scaffold_9', seq: fx.randomSeq(700, seed + 4) }
  ];
}

function softMask(seqs, name, from, to) {
  return seqs.map(function (s) {
    if (s.name !== name) return s;
    return { name: s.name, seq: s.seq.slice(0, from) + s.seq.slice(from, to).toLowerCase() + s.seq.slice(to) };
  });
}

function entriesOf(seqs) {
  return seqs.map(function (s) { return [s.name, s.seq.length]; });
}

// ---- dual-prefix dirs ---------------------------------------------------------------------------

['sorghum_pi534133', 'sorghum_pi576434', 'sorghum_pi656057'].forEach(function (sys) {
  test('dual-prefix dir ' + sys + ': picks the JGI prefix ending in "." + map._id', async function (t) {
    const root = tmpRoot(t);
    const map = sorghumMap(sys);
    const dir = path.join(root, sys);
    const acc = sys.slice('sorghum_'.length);
    const jgi = 'Sorghum_' + acc + '.' + map._id; // Sorghum_pi534133.Sb-PI534133-REFERENCE-JGI-1.0
    const old = 'Sorghum_' + acc + '.' + acc; // Sorghum_pi534133.pi534133
    const chroms = regionsOf(map).filter(function (e) { return e[0] !== 'UNANCHORED'; });
    chroms.length.should.equal(10);

    // Feb 2026 JGI prefix, as on disk: indexed dna, plain unindexed dna_sm/dna_rm .fa, dna/cdna/cds DBs.
    fx.indexOnlyFasta(dir, jgi, 'dna', chroms.concat([['scaffold_21', 269045]]));
    fx.touch(path.join(dir, 'dna', jgi + '.dna_sm.toplevel.fa'), '>1\nACGTacgt\n');
    fx.touch(path.join(dir, 'dna', jgi + '.dna_rm.toplevel.fa'), '>1\nACGTNNNN\n');
    ['dna', 'dna_sm', 'dna_rm'].forEach(function (k) { fx.blastDb(dir, jgi + '.' + k + '.toplevel', { nsqBytes: 300 }); });
    fx.blastDb(dir, jgi + '.cdna.all');
    fx.blastDb(dir, jgi + '.cds.all');
    // Jan 2025 prefix: other chromosome lengths, indexed dna/dna_sm/dna_rm copies and its own DBs.
    const oldChroms = chroms.map(function (e) { return [e[0], e[1] + 1000]; });
    ['dna', 'dna_sm', 'dna_rm'].forEach(function (k) {
      fx.indexOnlyFasta(dir, old, k, oldChroms);
      fx.blastDb(dir, old + '.' + k + '.toplevel', { nsqBytes: 400 });
    });
    fx.blastDb(dir, old + '.cdna.all');

    const scan = assemblies._internal.scanCandidates(fs.readdirSync(dir), fs.readdirSync(path.join(dir, 'dna')));
    scan.prefixes.should.eql([jgi, old].sort());

    const r = await assemblies.resolve(sys, depsFor(root, [map]));
    Object.keys(r).sort().should.eql(RESULT_KEYS);
    r.prefix.should.equal(jgi);
    r.dir.should.equal(dir);
    r.fasta.should.eql({ dna: path.join(dir, 'dna', jgi + '.dna.toplevel.fa.gz'), dna_sm: null });
    r.blastdb.should.eql({ dna: path.join(dir, jgi + '.dna.toplevel'), cdna: path.join(dir, jgi + '.cdna.all') });
    r.repeat_masking.should.equal('absent');
    r.warnings.should.eql([]); // 10 of 10 real map regions (the UNANCHORED bin is not a sequence and not counted)
    r.num_sequences.should.equal(11);
    r.total_bases.should.equal(chroms.reduce(function (s, e) { return s + e[1]; }, 0) + 269045);
    r.fingerprint.should.match(/^[0-9a-f]{40}$/);
    r.map_id.should.equal(map._id);
    r.taxon_id.should.equal(map.taxon_id);
    r.display_name.should.equal(map.display_name);
  });
});

test('the "." + map._id suffix wins over a better region score; a non-unique suffix falls back to scoring', async function (t) {
  const root = tmpRoot(t);
  const entries = [['1', 1000], ['2', 2000], ['3', 3000]];
  const dirA = path.join(root, 'sorghum_sfx');
  fx.indexOnlyFasta(dirA, 'Sorghum_sfx.ACC-2.0', 'dna', [['1', 1], ['2', 2]]);
  fx.indexOnlyFasta(dirA, 'Sorghum_sfx.old', 'dna', entries);
  const dirB = path.join(root, 'sorghum_two');
  fx.indexOnlyFasta(dirB, 'A.ACC-2.0', 'dna', [['1', 1]]);
  fx.indexOnlyFasta(dirB, 'B.ACC-2.0', 'dna', entries);
  const deps = depsFor(root, [mapDoc('sorghum_sfx', 'ACC-2.0', entries), mapDoc('sorghum_two', 'ACC-2.0', entries)]);

  const a = await assemblies.resolve('sorghum_sfx', deps);
  a.prefix.should.equal('Sorghum_sfx.ACC-2.0');
  codes(a).should.eql(['ASSEMBLY_MISMATCH']);
  a.warnings[0].message.should.equal('only 0 of 3 map regions match the assembly index by name and length');

  const b = await assemblies.resolve('sorghum_two', deps);
  b.prefix.should.equal('B.ACC-2.0');
  b.warnings.should.eql([]);
});

// ---- region scoring, sole candidate, ASSEMBLY_MISMATCH -------------------------------------------

test('region scoring picks the unique best .fai; sorghum_bicolor-like 10 of 11 regions gives no ASSEMBLY_MISMATCH', async function (t) {
  const root = tmpRoot(t);
  const map = sorghumMap('sorghum_bicolor');
  const dir = path.join(root, 'sorghum_bicolor');
  const chroms = regionsOf(map).filter(function (e) { return e[0] !== 'UNANCHORED'; });
  const scaffolds = [['super_10', 5000], ['super_11', 4000]];
  fx.indexOnlyFasta(dir, 'Sorghum_bicolor.Sorghum_bicolor_NCBIv3', 'dna', chroms.concat(scaffolds));
  fx.indexOnlyFasta(dir, 'Sorghum_bicolor.Other_v1', 'dna', chroms.slice(0, 3).concat(scaffolds));
  fx.blastDb(dir, 'Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel');
  fx.blastDb(dir, 'Sorghum_bicolor.Other_v1.dna.toplevel');

  const r = await assemblies.resolve('sorghum_bicolor', depsFor(root, [map]));
  r.prefix.should.equal('Sorghum_bicolor.Sorghum_bicolor_NCBIv3');
  r.warnings.should.eql([]);
  const fai = assemblies._internal.parseFai(fs.readFileSync(r.fasta.dna + '.fai', 'utf8'));
  assemblies._internal.regionMatches(fai, map.regions).should.equal(10);
  r.num_sequences.should.equal(12);
  r.blastdb.cdna === null && r.blastdb.dna.should.equal(path.join(dir, 'Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel'));
});

test('ASSEMBLY_MISMATCH below 80% of map regions (the sole candidate is still used); exactly 80% is fine', async function (t) {
  const root = tmpRoot(t);
  const entries = [['1', 100], ['2', 200], ['3', 300], ['4', 400], ['5', 500]];
  fx.indexOnlyFasta(path.join(root, 'sorghum_low'), 'Low.v1', 'dna', entries.slice(0, 3));
  fx.indexOnlyFasta(path.join(root, 'sorghum_edge'), 'Edge.v1', 'dna', entries.slice(0, 4).concat([['5', 501]]));
  const deps = depsFor(root, [mapDoc('sorghum_low', 'L', entries), mapDoc('sorghum_edge', 'E', entries)]);

  const low = await assemblies.resolve('sorghum_low', deps);
  low.prefix.should.equal('Low.v1');
  low.warnings.should.eql([{ code: 'ASSEMBLY_MISMATCH', message: 'only 3 of 5 map regions match the assembly index by name and length' }]);
  const edge = await assemblies.resolve('sorghum_edge', deps);
  edge.prefix.should.equal('Edge.v1');
  edge.warnings.should.eql([]);
});

test('synthetic map bins (UNANCHORED, UNPLACED) are not scored; a bin-only map gets no ASSEMBLY_MISMATCH; the denominator is real regions', async function (t) {
  const root = tmpRoot(t);
  // sorghum_tx430nano / selaginella_moellendorffii: the map is one UNANCHORED bin, the .fai only has contigs
  fx.indexOnlyFasta(path.join(root, 'sorghum_binonly'), 'Bin.v1', 'dna', [['ctg_1', 5000], ['ctg_2', 4000]]);
  // drosophila_melanogaster: 7 of 8 real regions match (the .fai calls the mitochondrion dmel_mitochondrion_genome)
  const fly = [['2L', 100], ['2R', 200], ['3L', 300], ['3R', 400], ['4', 50], ['X', 500], ['Y', 60], ['mitochondrion_genome', 19517]];
  fx.indexOnlyFasta(path.join(root, 'drosophila_like'), 'Fly.v6', 'dna',
    fly.slice(0, 7).concat([['dmel_mitochondrion_genome', 19517], ['211000022278279', 1500]]));
  // 6 of 8 real regions is still a mismatch, reported against the real regions only
  fx.indexOnlyFasta(path.join(root, 'drosophila_low'), 'Low.v6', 'dna', fly.slice(0, 6));
  const deps = depsFor(root, [
    mapDoc('sorghum_binonly', 'B1', [['UNANCHORED', 9000]]),
    mapDoc('drosophila_like', 'F6', fly.concat([['UNANCHORED', 6158518]])),
    mapDoc('drosophila_low', 'L6', fly.concat([['UNPLACED', 1000], ['unanchored', 5]]))
  ]);

  const bin = await assemblies.resolve('sorghum_binonly', deps);
  bin.prefix.should.equal('Bin.v1');
  bin.warnings.should.eql([]);
  const like = await assemblies.resolve('drosophila_like', deps);
  like.prefix.should.equal('Fly.v6');
  like.warnings.should.eql([]);
  const low = await assemblies.resolve('drosophila_low', deps);
  low.warnings.should.eql([{ code: 'ASSEMBLY_MISMATCH', message: 'only 6 of 8 map regions match the assembly index by name and length' }]);

  const internal = assemblies._internal;
  ['UNANCHORED', 'unanchored', 'UNPLACED', 'Unplaced'].forEach(function (n) { internal.isSyntheticRegion(n).should.be.true(); });
  ['Un', 'chrUn', 'Mt', 'Pt', 'UNANCHORED_1', 'mitochondrion_genome', '1'].forEach(function (n) { internal.isSyntheticRegion(n).should.be.false(); });
  internal.realRegionCount({ names: ['1', 'UNANCHORED', 'Un'], lengths: [1, 2, 3] }).should.equal(2);
  // a bin never adds to the score, even if a .fai had a sequence of that name and length
  internal.regionMatches(internal.parseFai(fx.faiText([['UNANCHORED', 9000], ['1', 10]])), { names: ['UNANCHORED', '1'], lengths: [9000, 10] })
    .should.equal(1);
});

// ---- ties ---------------------------------------------------------------------------------------

test('vitis-style tie with identical .fai and .nsq: newest BLAST DB wins with an AMBIGUOUS_ASSEMBLY warning', async function (t) {
  const root = tmpRoot(t);
  const entries = [['1', 23037639], ['2', 18779844], ['Un', 43154196]];
  const dir = path.join(root, 'vitis_vinifera');
  // FASTA mtimes point the other way on purpose: only the BLAST DB mtime decides.
  fx.indexOnlyFasta(dir, 'Vitis_vinifera.12X', 'dna', entries, { mtime: T2021 + 86400 });
  fx.indexOnlyFasta(dir, 'Vitis_vinifera.IGGP_12x', 'dna', entries.slice().reverse(), { mtime: T2020 });
  fx.blastDb(dir, 'Vitis_vinifera.12X.dna.toplevel', { nsqBytes: 1217, mtime: T2020 });
  fx.blastDb(dir, 'Vitis_vinifera.IGGP_12x.dna.toplevel', { nsqBytes: 1217, mtime: T2021 });
  fx.blastDb(dir, 'Vitis_vinifera.12X.cdna.all', { mtime: T2020 });

  const r = await assemblies.resolve('vitis_vinifera', depsFor(root, [mapDoc('vitis_vinifera', 'GCA_000003745.2', entries, 29760001)]));
  r.prefix.should.equal('Vitis_vinifera.IGGP_12x');
  r.warnings.should.eql([{ code: 'AMBIGUOUS_ASSEMBLY', message: '2 identical assemblies found; using the one with the newest BLAST database' }]);
  should(r.blastdb.cdna).be.null(); // files of the other prefix never leak into the result
});

test('a tie that is not identical (.nsq sizes, .fai sets, or nothing matches) is 422 AMBIGUOUS_ASSEMBLY', async function (t) {
  const root = tmpRoot(t);
  const entries = [['1', 1000], ['2', 2000]];
  // A: identical .fai, different .nsq sizes
  const a = path.join(root, 'sorghum_tie_nsq');
  fx.indexOnlyFasta(a, 'P.one', 'dna', entries);
  fx.indexOnlyFasta(a, 'P.two', 'dna', entries);
  fx.blastDb(a, 'P.one.dna.toplevel', { nsqBytes: 1217 });
  fx.blastDb(a, 'P.two.dna.toplevel', { nsqBytes: 1300 });
  // B: same score, but one .fai has an extra scaffold
  const b = path.join(root, 'sorghum_tie_fai');
  fx.indexOnlyFasta(b, 'P.one', 'dna', entries);
  fx.indexOnlyFasta(b, 'P.two', 'dna', entries.concat([['scaffold_1', 5]]));
  fx.blastDb(b, 'P.one.dna.toplevel', { nsqBytes: 1217 });
  fx.blastDb(b, 'P.two.dna.toplevel', { nsqBytes: 1217 });
  // C: neither candidate matches any map region and they differ
  const c = path.join(root, 'sorghum_tie_none');
  fx.indexOnlyFasta(c, 'P.one', 'dna', [['chr1', 5]]);
  fx.indexOnlyFasta(c, 'P.two', 'dna', [['chr1', 6]]);
  // D: one of two identical .fai sets has no BLAST DB
  const d = path.join(root, 'sorghum_tie_nodb');
  fx.indexOnlyFasta(d, 'P.one', 'dna', entries);
  fx.indexOnlyFasta(d, 'P.two', 'dna', entries);
  fx.blastDb(d, 'P.two.dna.toplevel');

  const names = ['sorghum_tie_nsq', 'sorghum_tie_fai', 'sorghum_tie_none', 'sorghum_tie_nodb'];
  const deps = depsFor(root, names.map(function (n) { return mapDoc(n, 'M', entries); }));
  for (const name of names) {
    const err = await rejectsWith(assemblies.resolve(name, deps), 422, 'AMBIGUOUS_ASSEMBLY');
    err.details.should.eql({ system_name: name, candidates: ['P.one', 'P.two'] });
    err.message.should.not.containEql(root);
    err.message.should.containEql('assembly_overrides');
  }
});

// ---- BLAST DB layouts ---------------------------------------------------------------------------

test('multi-volume .nal: alias preferred over .nin, volume files are not candidates, fingerprint follows alias + volumes', async function (t) {
  const root = tmpRoot(t);
  const entries = [['1H', 5000], ['2H', 6000]];
  const dir = path.join(root, 'hordeum_x');
  const prefix = 'Hordeum_x.IBSC_v2';
  fx.indexOnlyFasta(dir, prefix, 'dna', entries);
  fx.blastDb(dir, prefix + '.dna.toplevel', { volumes: 2, nsqBytes: 500 });
  fx.blastDb(dir, prefix + '.dna.toplevel', { nsqBytes: 999 }); // stray single-volume files beside the alias
  fx.blastDb(dir, prefix + '.dna_sm.toplevel', { volumes: 2 });
  fx.blastDb(dir, prefix + '.cdna.all', { volumes: 2 });

  const scan = assemblies._internal.scanCandidates(fs.readdirSync(dir), fs.readdirSync(path.join(dir, 'dna')));
  scan.prefixes.should.eql([prefix]);
  const db = assemblies._internal.blastDbFiles(scan, prefix, 'dna.toplevel');
  db.alias.should.equal(prefix + '.dna.toplevel.nal');
  db.volumes.should.eql([prefix + '.dna.toplevel.00.nsq', prefix + '.dna.toplevel.01.nsq']);
  assemblies._internal.scanCandidates([prefix + '.dna.toplevel.00.nin', prefix + '.dna.toplevel.01.nin'], [])
    .prefixes.should.eql([]);

  const deps = depsFor(root, [mapDoc('hordeum_x', 'IBSC_v2', entries)]);
  const r1 = await assemblies.resolve('hordeum_x', deps);
  r1.prefix.should.equal(prefix);
  r1.blastdb.should.eql({ dna: path.join(dir, prefix + '.dna.toplevel'), cdna: path.join(dir, prefix + '.cdna.all') });

  fx.touch(path.join(dir, prefix + '.dna.toplevel.nsq'), Buffer.alloc(1234), T2021);
  (await assemblies.resolve('hordeum_x', deps)).fingerprint.should.equal(r1.fingerprint);
  fx.touch(path.join(dir, prefix + '.dna.toplevel.01.nsq'), Buffer.alloc(501));
  const r3 = await assemblies.resolve('hordeum_x', deps);
  r3.fingerprint.should.not.equal(r1.fingerprint);
  fx.touch(path.join(dir, prefix + '.dna.toplevel.nal'), 'TITLE changed\n', T2020);
  (await assemblies.resolve('hordeum_x', deps)).fingerprint.should.not.equal(r3.fingerprint);
});

test('a BLAST-only assembly resolves with no FASTA and no sizes; a .nin without .nsq is not a usable DB', async function (t) {
  const root = tmpRoot(t);
  const dir = path.join(root, 'sorghum_blastonly');
  fx.blastDb(dir, 'Sorghum_blastonly.v1.dna.toplevel', { volumes: 3 });
  fx.blastDb(dir, 'Sorghum_blastonly.v1.cdna.all', { nsq: false });
  const r = await assemblies.resolve('sorghum_blastonly', depsFor(root, [mapDoc('sorghum_blastonly', 'v1', [['1', 10]])]));
  r.prefix.should.equal('Sorghum_blastonly.v1');
  r.fasta.should.eql({ dna: null, dna_sm: null });
  r.blastdb.should.eql({ dna: path.join(dir, 'Sorghum_blastonly.v1.dna.toplevel'), cdna: null });
  should(r.total_bases).be.null();
  should(r.num_sequences).be.null();
  r.repeat_masking.should.equal('absent');
  r.warnings.should.eql([]);
  r.fingerprint.should.match(/^[0-9a-f]{40}$/);
});

test('a catalog genome without assembly files resolves to nulls; results are frozen and JSON-safe', async function (t) {
  const root = tmpRoot(t);
  fs.mkdirSync(path.join(root, 'sorghum_empty', 'dna'), { recursive: true });
  fx.touch(path.join(root, 'sorghum_empty', 'dna', 'README'), 'x');
  fx.touch(path.join(root, 'sorghum_empty', 'dna', 'X.dna.toplevel.fa.gz'), 'x'); // no .fai/.gzi
  const maps = [mapDoc('sorghum_empty', 'E', [['1', 1]]), mapDoc('sorghum_nodir', 'N', [['1', 1]])];
  const deps = depsFor(root, maps);
  for (const name of ['sorghum_empty', 'sorghum_nodir']) {
    const r = await assemblies.resolve(name, deps);
    Object.keys(r).sort().should.eql(RESULT_KEYS);
    r.should.eql({
      system_name: name, taxon_id: 4558999, display_name: name, map_id: name === 'sorghum_empty' ? 'E' : 'N',
      prefix: null, dir: path.join(root, name), fasta: { dna: null, dna_sm: null }, blastdb: { dna: null, cdna: null },
      repeat_masking: 'absent', total_bases: null, num_sequences: null, fingerprint: null, warnings: []
    });
    [r, r.fasta, r.blastdb, r.warnings].forEach(function (o) { Object.isFrozen(o).should.be.true(); });
    JSON.parse(JSON.stringify(r)).should.eql(r);
  }
});

// ---- soft-mask detection ------------------------------------------------------------------------

test('the BGZF fixture writer round-trips through @gmod/indexedfasta (multi-block)', async function (t) {
  const root = tmpRoot(t);
  const seqs = genomeSeqs(7);
  const file = fx.bgzipFasta(path.join(root, 'g'), 'G.v1', 'dna', seqs);
  fs.readFileSync(file + '.gzi').readBigUInt64LE(0).should.be.above(0n);
  const { BgzipIndexedFasta } = require('@gmod/indexedfasta');
  const f = new BgzipIndexedFasta({ path: file, faiPath: file + '.fai', gziPath: file + '.gzi' });
  (await f.getSequence('1', 59000, 60000)).should.equal(seqs[0].seq.slice(59000, 60000));
  (await f.getSequence('3', 0, 20000)).should.equal(seqs[2].seq);
  (await f.getSequence('scaffold_9', 100, 200)).should.equal(seqs[3].seq.slice(100, 200));
});

test('repeat_masking: lowercase in sampled dna_sm windows -> soft_masked; uppercase copy -> unmasked_copy; no dna_sm -> absent', async function (t) {
  const root = tmpRoot(t);
  const seqs = genomeSeqs(11);
  const P = 'Sorghum_x.V1';
  const names = ['sorghum_rio', 'sorghum_copy', 'sorghum_plain'];
  const rio = path.join(root, 'sorghum_rio');
  fx.bgzipFasta(rio, P, 'dna', seqs);
  fx.bgzipFasta(rio, P, 'dna_sm', softMask(seqs, '1', 29000, 31000)); // at the midpoint of the longest sequence
  fx.bgzipFasta(rio, P, 'dna_rm', seqs);
  const copy = path.join(root, 'sorghum_copy');
  fx.bgzipFasta(copy, P, 'dna', seqs);
  fx.bgzipFasta(copy, P, 'dna_sm', seqs);
  fx.bgzipFasta(path.join(root, 'sorghum_plain'), P, 'dna', seqs);
  const deps = depsFor(root, names.map(function (n) { return mapDoc(n, 'V1', entriesOf(seqs)); }));

  const r = await assemblies.resolve('sorghum_rio', deps);
  r.repeat_masking.should.equal('soft_masked');
  r.fasta.should.eql({ dna: path.join(rio, 'dna', P + '.dna.toplevel.fa.gz'), dna_sm: path.join(rio, 'dna', P + '.dna_sm.toplevel.fa.gz') });
  r.total_bases.should.equal(110700);
  r.num_sequences.should.equal(4);
  r.warnings.should.eql([]);

  const c = await assemblies.resolve('sorghum_copy', deps);
  c.repeat_masking.should.equal('unmasked_copy');
  c.fasta.dna_sm.should.equal(path.join(copy, 'dna', P + '.dna_sm.toplevel.fa.gz'));

  const p = await assemblies.resolve('sorghum_plain', deps);
  p.repeat_masking.should.equal('absent');
  should(p.fasta.dna_sm).be.null();
});

test('soft-mask sampling reads 20 kb windows from the 3 longest sequences and is cached per file identity', async function (t) {
  const root = tmpRoot(t);
  const seqs = genomeSeqs(21);
  const fai = assemblies._internal.parseFai(fx.faiText(entriesOf(seqs)));
  assemblies._internal.maskSampleWindows(fai).should.eql([
    { name: '1', start: 20000, end: 40000 }, { name: '2', start: 5000, end: 25000 }, { name: '3', start: 0, end: 20000 },
    { name: '1', start: 5000, end: 25000 }, { name: '2', start: 0, end: 20000 },
    { name: '1', start: 35000, end: 55000 }, { name: '2', start: 10000, end: 30000 }
  ]);

  // Lowercase only inside the 1/4 window of sequence 1: found on the 4th read.
  const P = 'Sorghum_q.V1';
  const quarter = path.join(root, 'sorghum_quarter');
  fx.bgzipFasta(quarter, P, 'dna', seqs);
  const smFile = fx.bgzipFasta(quarter, P, 'dna_sm', softMask(seqs, '1', 6000, 6100));
  const copy = path.join(root, 'sorghum_copy2');
  fx.bgzipFasta(copy, P, 'dna', seqs);
  fx.bgzipFasta(copy, P, 'dna_sm', seqs);

  const reads = [];
  const { BgzipIndexedFasta } = require('@gmod/indexedfasta');
  const openFasta = function (file) {
    const real = new BgzipIndexedFasta({ path: file, faiPath: file + '.fai', gziPath: file + '.gzi' });
    return { getSequence: function (name, start, end) { reads.push([path.basename(path.dirname(path.dirname(file))), name, start, end]); return real.getSequence(name, start, end); } };
  };
  const maps = ['sorghum_quarter', 'sorghum_copy2'].map(function (n) { return mapDoc(n, 'V1', entriesOf(seqs)); });
  const deps = depsFor(root, maps, { openFasta: openFasta });

  (await assemblies.resolve('sorghum_quarter', deps)).repeat_masking.should.equal('soft_masked');
  reads.length.should.equal(4);
  reads[3].should.eql(['sorghum_quarter', '1', 5000, 25000]);
  (await assemblies.resolve('sorghum_quarter', deps)).repeat_masking.should.equal('soft_masked');
  reads.length.should.equal(4); // cached although the resolve cache is off

  fs.utimesSync(smFile, T2021, T2021); // file identity changed -> sampled again
  (await assemblies.resolve('sorghum_quarter', deps)).repeat_masking.should.equal('soft_masked');
  reads.length.should.equal(8);

  (await assemblies.resolve('sorghum_copy2', deps)).repeat_masking.should.equal('unmasked_copy');
  reads.length.should.equal(15); // every window read before concluding there is no mask
  reads.slice(8).every(function (r) { return r[0] === 'sorghum_copy2' && r[1] !== 'scaffold_9'; }).should.be.true();
});

test('dna_sm is ignored without .gzi, with an index that differs from dna, or when it cannot be read', async function (t) {
  const root = tmpRoot(t);
  const seqs = genomeSeqs(31);
  const P = 'Sorghum_s.V1';
  const nogzi = path.join(root, 'sorghum_nogzi');
  fx.bgzipFasta(nogzi, P, 'dna', seqs);
  fx.bgzipFasta(nogzi, P, 'dna_sm', softMask(seqs, '1', 20000, 40000), { gzi: false });
  const diff = path.join(root, 'sorghum_diff');
  fx.bgzipFasta(diff, P, 'dna', seqs);
  fx.bgzipFasta(diff, P, 'dna_sm', softMask(seqs, '1', 20000, 40000).slice(0, 3));
  const broken = path.join(root, 'sorghum_broken');
  fx.bgzipFasta(broken, P, 'dna', seqs);
  fx.indexOnlyFasta(broken, P, 'dna_sm', entriesOf(seqs)); // index says the sequences are there; the data is not BGZF
  const logged = [];
  const names = ['sorghum_nogzi', 'sorghum_diff', 'sorghum_broken'];
  const deps = depsFor(root, names.map(function (n) { return mapDoc(n, 'V1', entriesOf(seqs)); }),
    { log: { warn: function (m) { logged.push(m); } } });
  for (const name of names) {
    const r = await assemblies.resolve(name, deps);
    should(r.fasta.dna_sm).be.null();
    r.repeat_masking.should.equal('absent');
    r.fasta.dna.should.equal(path.join(root, name, 'dna', P + '.dna.toplevel.fa.gz'));
  }
  logged.length.should.equal(2);
  logged.join('\n').should.not.containEql(root);
});

test('repeat_masking_overrides win; a masked override without a usable dna_sm stays absent', async function (t) {
  const root = tmpRoot(t);
  const seqs = genomeSeqs(41);
  const P = 'Sorghum_o.V1';
  const masked = path.join(root, 'sorghum_masked');
  fx.bgzipFasta(masked, P, 'dna', seqs);
  fx.bgzipFasta(masked, P, 'dna_sm', softMask(seqs, '1', 20000, 40000));
  const plain = path.join(root, 'sorghum_plain2');
  fx.bgzipFasta(plain, P, 'dna', seqs);
  const bogus = path.join(root, 'sorghum_bogus');
  fx.bgzipFasta(bogus, P, 'dna', seqs);
  fx.bgzipFasta(bogus, P, 'dna_sm', softMask(seqs, '1', 20000, 40000));
  let opened = 0;
  const names = ['sorghum_masked', 'sorghum_plain2', 'sorghum_bogus'];
  const deps = depsFor(root, names.map(function (n) { return mapDoc(n, 'V1', entriesOf(seqs)); }), {
    cfg: {
      fasta_root: root, assembly_overrides: {},
      repeat_masking_overrides: { sorghum_masked: 'unmasked_copy', sorghum_plain2: 'soft_masked', sorghum_bogus: 'bogus' }
    },
    openFasta: function (file) {
      opened++;
      const { BgzipIndexedFasta } = require('@gmod/indexedfasta');
      return new BgzipIndexedFasta({ path: file, faiPath: file + '.fai', gziPath: file + '.gzi' });
    }
  });
  (await assemblies.resolve('sorghum_masked', deps)).repeat_masking.should.equal('unmasked_copy');
  opened.should.equal(0);
  (await assemblies.resolve('sorghum_plain2', deps)).repeat_masking.should.equal('absent');
  (await assemblies.resolve('sorghum_bogus', deps)).repeat_masking.should.equal('soft_masked'); // invalid value ignored
  opened.should.equal(1);
});

// ---- overrides, security, config ----------------------------------------------------------------

test('assembly_overrides wins over the map._id suffix; an override naming no candidate is 422 AMBIGUOUS_ASSEMBLY', async function (t) {
  const root = tmpRoot(t);
  const entries = [['1', 1000]];
  const dir = path.join(root, 'sorghum_ovr');
  fx.indexOnlyFasta(dir, 'Sorghum_ovr.M1', 'dna', entries);
  fx.indexOnlyFasta(dir, 'Sorghum_ovr.legacy', 'dna', entries);
  const maps = [mapDoc('sorghum_ovr', 'M1', entries)];
  const cfg = function (value) {
    return { fasta_root: root, assembly_overrides: { sorghum_ovr: value }, repeat_masking_overrides: {} };
  };
  (await assemblies.resolve('sorghum_ovr', depsFor(root, maps, { cfg: cfg('Sorghum_ovr.legacy') }))).prefix
    .should.equal('Sorghum_ovr.legacy');
  (await assemblies.resolve('sorghum_ovr', depsFor(root, maps, { cfg: cfg(null) }))).prefix.should.equal('Sorghum_ovr.M1');
  for (const bad of ['Sorghum_ovr.missing', '../sorghum_ovr/dna/Sorghum_ovr.legacy', 42]) {
    const err = await rejectsWith(assemblies.resolve('sorghum_ovr', depsFor(root, maps, { cfg: cfg(bad) })), 422, 'AMBIGUOUS_ASSEMBLY');
    err.details.candidates.should.eql(['Sorghum_ovr.M1', 'Sorghum_ovr.legacy']);
  }
});

test('system_name must match ^[a-z0-9_]+$ and be in the catalog before any filesystem access', async function (t) {
  const root = tmpRoot(t);
  fx.indexOnlyFasta(path.join(root, 'sorghum_ok'), 'Ok.v1', 'dna', [['1', 10]]);
  fs.mkdirSync(path.join(root, 'etc'));
  const spy = fx.spyFs();
  const catalog = genomes.buildCatalog([mapDoc('sorghum_ok', 'v1', [['1', 10]]), mapDoc('../etc', 'x', [['1', 1]])], []);
  catalog.bySystemName.has('../etc').should.be.false();
  const deps = { cfg: { fasta_root: root }, catalog: catalog, fs: spy, log: SILENT };
  const bad = ['../x', '../etc', '..', '.', 'a/b', 'sorghum_ok/../sorghum_ok', '/etc', 'Sorghum_ok', 'sorghum-ok', '',
    ' sorghum_ok', 'sorghum_ok\n', 'sorghum_ok ', 'x'.repeat(129), null, undefined, 42, {}, ['sorghum_ok'], 'etc', 'not_in_catalog'];
  for (const name of bad) {
    const err = await rejectsWith(assemblies.resolve(name, deps), 404, 'UNKNOWN_GENOME');
    if (typeof name === 'string' && !/^[a-z0-9_]+$/.test(name)) err.details.should.eql({});
  }
  spy.calls.should.eql([]);

  const r = await assemblies.resolve('sorghum_ok', deps);
  r.prefix.should.equal('Ok.v1');
  spy.calls.length.should.be.above(0);
  spy.calls.forEach(function (c) {
    const rel = path.relative(path.join(root, 'sorghum_ok'), c[1]);
    (rel === '' || (!rel.startsWith('..') && !path.isAbsolute(rel))).should.be.true();
  });
});

test('fasta_root must be configured as an absolute path (500 INTERNAL)', async function (t) {
  const catalog = genomes.buildCatalog([mapDoc('sorghum_ok', 'v1', [['1', 10]])], []);
  for (const cfg of [{ fasta_root: 'relative/fasta' }, { fasta_root: '' }, {}]) {
    await rejectsWith(assemblies.resolve('sorghum_ok', { cfg: cfg, catalog: catalog, log: SILENT }), 500, 'INTERNAL');
  }
});

// ---- caching and resolveMany --------------------------------------------------------------------

test('resolve results are cached for 10 minutes (422 answers too); clearCache() forgets them', async function (t) {
  const root = tmpRoot(t);
  fx.indexOnlyFasta(path.join(root, 'sorghum_c1'), 'C1.v1', 'dna', [['1', 10]]);
  fx.indexOnlyFasta(path.join(root, 'sorghum_amb'), 'A.one', 'dna', [['1', 10]]);
  fx.indexOnlyFasta(path.join(root, 'sorghum_amb'), 'A.two', 'dna', [['1', 11]]);
  let now = 1e6;
  const spy = fx.spyFs();
  const maps = [mapDoc('sorghum_c1', 'v1', [['1', 10]]), mapDoc('sorghum_amb', 'x', [['2', 10]])];
  const deps = depsFor(root, maps, { cache: true, fs: spy, now: function () { return now; } });

  const r1 = await assemblies.resolve('sorghum_c1', deps);
  const n1 = spy.calls.length;
  (await assemblies.resolve('sorghum_c1', deps)).should.equal(r1);
  spy.calls.length.should.equal(n1);

  await rejectsWith(assemblies.resolve('sorghum_amb', deps), 422, 'AMBIGUOUS_ASSEMBLY');
  const n2 = spy.calls.length;
  await rejectsWith(assemblies.resolve('sorghum_amb', deps), 422, 'AMBIGUOUS_ASSEMBLY');
  spy.calls.length.should.equal(n2);

  now += assemblies.RESOLVE_TTL_MS - 1;
  (await assemblies.resolve('sorghum_c1', deps)).should.equal(r1);
  now += 1;
  const r3 = await assemblies.resolve('sorghum_c1', deps);
  r3.should.not.equal(r1);
  r3.should.eql(r1);
  spy.calls.length.should.be.above(n2);

  assemblies.clearCache();
  const n3 = spy.calls.length;
  (await assemblies.resolve('sorghum_c1', deps)).should.not.equal(r3);
  spy.calls.length.should.be.above(n3);
});

test('transient filesystem errors are not cached', async function (t) {
  const root = tmpRoot(t);
  fx.indexOnlyFasta(path.join(root, 'sorghum_eio'), 'E.v1', 'dna', [['1', 10]]);
  let fail = true;
  const spy = fx.spyFs({
    readdir: function (p) {
      if (fail) { const e = new Error('EACCES: permission denied, scandir \'' + p + '\''); e.code = 'EACCES'; return Promise.reject(e); }
      return fs.promises.readdir(p);
    }
  });
  const deps = depsFor(root, [mapDoc('sorghum_eio', 'v1', [['1', 10]])], { cache: true, fs: spy });
  await assert.rejects(assemblies.resolve('sorghum_eio', deps), /EACCES/);
  fail = false;
  (await assemblies.resolve('sorghum_eio', deps)).prefix.should.equal('E.v1');
});

test('resolveMany keeps input order, captures per-genome errors and resolves at most 8 genomes at a time', async function (t) {
  const root = tmpRoot(t);
  const names = [];
  for (let i = 0; i < 20; i++) {
    const name = 'sorghum_g' + String(i).padStart(2, '0');
    names.push(name);
    fx.indexOnlyFasta(path.join(root, name), 'G.v' + i, 'dna', [['1', 10]]);
    if (i === 5) fx.indexOnlyFasta(path.join(root, name), 'G.other', 'dna', [['1', 11]]);
  }
  const inflight = new Map();
  let maxActive = 0;
  const spy = fx.spyFs({
    readdir: async function (p) {
      const g = path.relative(root, p).split(path.sep)[0];
      inflight.set(g, (inflight.get(g) || 0) + 1);
      maxActive = Math.max(maxActive, Array.from(inflight.values()).filter(function (n) { return n > 0; }).length);
      await new Promise(function (r) { setTimeout(r, 10); });
      inflight.set(g, inflight.get(g) - 1);
      return fs.promises.readdir(p);
    }
  });
  const maps = names.map(function (n) { return mapDoc(n, 'nomatch', [['2', 10]]); });
  const input = names.concat(['not_in_catalog', '../etc']);
  const out = await assemblies.resolveMany(input, depsFor(root, maps, { fs: spy }));

  out.map(function (o) { return o.system_name; }).should.eql(input);
  out.forEach(function (o, i) {
    if (i === 5) {
      should(o.resolved).be.null();
      o.error.code.should.equal('AMBIGUOUS_ASSEMBLY');
    } else if (i >= 20) {
      should(o.resolved).be.null();
      o.error.code.should.equal('UNKNOWN_GENOME');
    } else {
      should(o.error).be.null();
      o.resolved.prefix.should.equal('G.v' + i);
    }
  });
  maxActive.should.be.within(2, assemblies.RESOLVE_CONCURRENCY);
  (await assemblies.resolveMany([], {})).should.eql([]);
});

test('resolveMany rejects with 503 MONGO_UNAVAILABLE when the catalog cannot be loaded', async function () {
  genomes.clearCache();
  const mongo = {
    maps: { mongoCollection: async function () { return undefined; } },
    taxonomy: { mongoCollection: async function () { return undefined; } }
  };
  await rejectsWith(assemblies.resolveMany(['sorghum_bicolor'], { mongo: mongo, cfg: { fasta_root: '/nonexistent' }, log: SILENT }),
    503, 'MONGO_UNAVAILABLE');
});

// ---- real data (opt-in) -------------------------------------------------------------------------

const REAL = process.env.PRIMERS_REALDATA === '1';
const FASTA_ROOT = process.env.PRIMERS_FASTA_ROOT || '/scratch/olson/fasta';
const REAL_SKIP = REAL ? false : 'set PRIMERS_REALDATA=1 to resolve against ' + FASTA_ROOT + ' and mongo';
const SOFT_MASKED_SORGHUM = ['sorghum_rio', 'sorghum_tx2783pac', 'sorghum_tx430nano', 'sorghum_tx436pac'];
const DUAL_PREFIX = ['sorghum_pi534133', 'sorghum_pi576434', 'sorghum_pi656057'];

function report(t, label, value) {
  const line = label + ': ' + (typeof value === 'string' ? value : JSON.stringify(value));
  t.diagnostic(line);
  console.log('# ' + line);
}

test('real data: every sorghum_v11 genome map resolves against the fasta root', { skip: REAL_SKIP, timeout: 10 * 60 * 1000 }, async function (t) {
  const mongo = require('gramene-mongodb-config');
  t.after(function () { mongo.closeMongoDatabase(); });
  genomes.clearCache();
  assemblies.clearCache();
  const cfg = { fasta_root: FASTA_ROOT, assembly_overrides: {}, repeat_masking_overrides: {} };
  const catalog = await genomes.getCatalog({ mongo: mongo });
  const names = catalog.genomes.map(function (g) { return g.system_name; });
  report(t, 'catalog genomes', names.length);

  const t0 = Date.now();
  const out = await assemblies.resolveMany(names, { cfg: cfg, catalog: catalog });
  const coldMs = Date.now() - t0;

  const tally = { genomes: out.length, errors: {}, has_sequence: 0, has_blastdb: 0, has_cdna_blastdb: 0,
    repeat_masking: { soft_masked: 0, unmasked_copy: 0, absent: 0 }, warnings: {} };
  const sorghumSoft = [];
  out.forEach(function (o) {
    if (o.error) {
      tally.errors[o.system_name] = o.error.code + ' ' + o.error.message;
      return;
    }
    const r = o.resolved;
    if (r.fasta.dna) tally.has_sequence++;
    if (r.blastdb.dna) tally.has_blastdb++;
    if (r.blastdb.cdna) tally.has_cdna_blastdb++;
    tally.repeat_masking[r.repeat_masking]++;
    r.warnings.forEach(function (w) { (tally.warnings[w.code] = tally.warnings[w.code] || []).push(o.system_name); });
    if (r.repeat_masking === 'soft_masked' && o.system_name.indexOf('sorghum_') === 0) sorghumSoft.push(o.system_name);
  });
  report(t, 'tally', tally);
  const sorghumTally = { genomes: 0, has_sequence: 0, has_blastdb: 0, has_cdna_blastdb: 0,
    repeat_masking: { soft_masked: 0, unmasked_copy: 0, absent: 0 } };
  out.forEach(function (o) {
    if (o.error || o.system_name.indexOf('sorghum_') !== 0) return;
    const r = o.resolved;
    sorghumTally.genomes++;
    if (r.fasta.dna) sorghumTally.has_sequence++;
    if (r.blastdb.dna) sorghumTally.has_blastdb++;
    if (r.blastdb.cdna) sorghumTally.has_cdna_blastdb++;
    sorghumTally.repeat_masking[r.repeat_masking]++;
  });
  report(t, 'sorghum-only tally', sorghumTally);
  report(t, 'ASSEMBLY_MISMATCH details', out.filter(function (o) {
    return o.resolved && codes(o.resolved).indexOf('ASSEMBLY_MISMATCH') >= 0;
  }).map(function (o) {
    return o.system_name + ' (' + o.resolved.prefix + '): ' +
      o.resolved.warnings.find(function (w) { return w.code === 'ASSEMBLY_MISMATCH'; }).message;
  }));
  report(t, 'cold resolve ms', coldMs);
  sorghumTally.genomes.should.equal(120);
  report(t, 'soft_masked non-sorghum', out.filter(function (o) {
    return o.resolved && o.resolved.repeat_masking === 'soft_masked' && o.system_name.indexOf('sorghum_') !== 0;
  }).map(function (o) { return o.system_name; }));
  report(t, 'not fully provisioned', out.filter(function (o) {
    return o.resolved && !(o.resolved.fasta.dna && o.resolved.blastdb.dna && o.resolved.blastdb.cdna);
  }).map(function (o) { return o.system_name + ' seq=' + !!o.resolved.fasta.dna + ' db=' + !!o.resolved.blastdb.dna + ' cdna=' + !!o.resolved.blastdb.cdna; }));

  names.length.should.equal(128);
  tally.errors.should.eql({});
  sorghumSoft.sort().should.eql(SOFT_MASKED_SORGHUM);

  const byName = new Map(out.map(function (o) { return [o.system_name, o.resolved]; }));
  DUAL_PREFIX.forEach(function (sys) {
    const r = byName.get(sys);
    report(t, sys, { prefix: r.prefix, repeat_masking: r.repeat_masking, blastdb: !!r.blastdb.dna, cdna: !!r.blastdb.cdna, warnings: codes(r) });
    r.prefix.should.equal('Sorghum_' + sys.slice(8) + '.' + r.map_id);
    r.blastdb.dna.should.equal(path.join(FASTA_ROOT, sys, r.prefix + '.dna.toplevel'));
  });

  const bicolor = byName.get('sorghum_bicolor');
  const fai = assemblies._internal.parseFai(fs.readFileSync(bicolor.fasta.dna + '.fai', 'utf8'));
  const bicolorRegions = catalog.bySystemName.get('sorghum_bicolor').regions;
  const matched = assemblies._internal.regionMatches(fai, bicolorRegions);
  report(t, 'sorghum_bicolor', { prefix: bicolor.prefix, region_matches: matched + '/' + assemblies._internal.realRegionCount(bicolorRegions) +
    ' real (' + bicolorRegions.names.length + ' map regions)', total_bases: bicolor.total_bases, num_sequences: bicolor.num_sequences,
    repeat_masking: bicolor.repeat_masking, warnings: codes(bicolor) });
  matched.should.equal(10);
  assemblies._internal.realRegionCount(bicolorRegions).should.equal(10);
  bicolor.prefix.should.equal('Sorghum_bicolor.Sorghum_bicolor_NCBIv3');
  codes(bicolor).should.not.containEql('ASSEMBLY_MISMATCH');
  bicolor.repeat_masking.should.equal('unmasked_copy');

  // Synthetic UNANCHORED bins are not regions: bin-only maps and drosophila (7 of 8 real regions) give no ASSEMBLY_MISMATCH.
  const binCases = {};
  ['sorghum_tx430nano', 'selaginella_moellendorffii', 'drosophila_melanogaster'].forEach(function (sys) {
    const r = byName.get(sys);
    const regions = catalog.bySystemName.get(sys).regions;
    const f = assemblies._internal.parseFai(fs.readFileSync(r.fasta.dna + '.fai', 'utf8'));
    binCases[sys] = {
      prefix: r.prefix, map_regions: regions.names.length, real_regions: assemblies._internal.realRegionCount(regions),
      matched: assemblies._internal.regionMatches(f, regions),
      unmatched: regions.names.filter(function (n, i) { return f.lengths.get(n) !== regions.lengths[i]; }), warnings: codes(r)
    };
    codes(r).should.not.containEql('ASSEMBLY_MISMATCH');
  });
  report(t, 'synthetic-bin genomes', binCases);
  binCases.sorghum_tx430nano.real_regions.should.equal(0);
  binCases.selaginella_moellendorffii.real_regions.should.equal(0);
  [binCases.drosophila_melanogaster.matched, binCases.drosophila_melanogaster.real_regions].should.eql([7, 8]);
  binCases.drosophila_melanogaster.unmatched.sort().should.eql(['UNANCHORED', 'mitochondrion_genome']);

  const t1 = Date.now();
  await assemblies.resolveMany(names, { cfg: cfg, catalog: catalog });
  const warmMs = Date.now() - t1;
  report(t, 'cached resolve ms', warmMs);
  warmMs.should.be.below(1000);

  const body = await genomes.genomesResponse('sorghum_bicolor', { cfg: cfg, catalog: catalog });
  report(t, 'genomes response counts', body.counts);
  body.genomes.length.should.equal(120);
  body.genomes[0].system_name.should.equal('sorghum_bicolor');
  body.genomes[0].is_query.should.be.true();
  body.genomes.find(function (g) { return g.system_name === 'sorghum_rio'; }).repeat_masking.should.equal('soft_masked');
  JSON.stringify(body).should.not.containEql(FASTA_ROOT);
  JSON.stringify(body).should.not.containEql('.toplevel');
});

test('real data: multi-volume .nal BLAST DBs resolve (read-only, v69 map docs)', { skip: REAL_SKIP, timeout: 5 * 60 * 1000 }, async function (t) {
  const nal = [];
  fs.readdirSync(FASTA_ROOT).forEach(function (d) {
    let names;
    try { names = fs.readdirSync(path.join(FASTA_ROOT, d)); } catch (e) { return; }
    names.filter(function (n) { return /\.nal$/.test(n); }).forEach(function (n) { nal.push(d + '/' + n); });
  });
  report(t, '.nal files in fasta root', nal);

  const systems = ['aegilops_tauschii', 'hordeum_vulgare', 'hordeum_vulgare_goldenpromise'];
  const maps = V69_MAPS.filter(function (m) { return systems.indexOf(m.system_name) >= 0; });
  maps.length.should.equal(3);
  const catalog = genomes.buildCatalog(maps, []);
  const cfg = { fasta_root: FASTA_ROOT, assembly_overrides: {}, repeat_masking_overrides: {} };
  const out = await assemblies.resolveMany(systems, { cfg: cfg, catalog: catalog, cache: false });
  out.forEach(function (o) {
    report(t, o.system_name, o.error ? { error: o.error.code, details: o.error.details } : {
      prefix: o.resolved.prefix, blastdb: o.resolved.blastdb.dna && path.basename(o.resolved.blastdb.dna),
      cdna: o.resolved.blastdb.cdna && path.basename(o.resolved.blastdb.cdna), repeat_masking: o.resolved.repeat_masking,
      num_sequences: o.resolved.num_sequences, warnings: codes(o.resolved)
    });
  });

  // End to end through the multi-volume prefix of hordeum_vulgare (forced: the v69 map is MorexV3).
  const dir = path.join(FASTA_ROOT, 'hordeum_vulgare');
  const scan = assemblies._internal.scanCandidates(fs.readdirSync(dir), fs.readdirSync(path.join(dir, 'dna')));
  report(t, 'hordeum_vulgare candidates', scan.prefixes);
  scan.prefixes.should.containEql('Hordeum_vulgare.IBSC_v2');
  scan.prefixes.every(function (p) { return !/\.\d\d$/.test(p); }).should.be.true();
  const db = assemblies._internal.blastDbFiles(scan, 'Hordeum_vulgare.IBSC_v2', 'dna.toplevel');
  db.alias.should.equal('Hordeum_vulgare.IBSC_v2.dna.toplevel.nal');
  db.volumes.should.eql(['Hordeum_vulgare.IBSC_v2.dna.toplevel.00.nsq', 'Hordeum_vulgare.IBSC_v2.dna.toplevel.01.nsq']);
  const forced = await assemblies.resolve('hordeum_vulgare', {
    cfg: Object.assign({}, cfg, { assembly_overrides: { hordeum_vulgare: 'Hordeum_vulgare.IBSC_v2' } }), catalog: catalog, cache: false
  });
  forced.blastdb.dna.should.equal(path.join(dir, 'Hordeum_vulgare.IBSC_v2.dna.toplevel'));
  forced.fasta.dna.should.equal(path.join(dir, 'dna', 'Hordeum_vulgare.IBSC_v2.dna.toplevel.fa.gz'));
  forced.fingerprint.should.match(/^[0-9a-f]{40}$/);
  report(t, 'hordeum_vulgare forced IBSC_v2', { repeat_masking: forced.repeat_masking, warnings: codes(forced), num_sequences: forced.num_sequences });

  const aegilops = out.find(function (o) { return o.system_name === 'aegilops_tauschii'; });
  should(aegilops.error).be.null();
  aegilops.resolved.blastdb.dna.should.endWith('.dna.toplevel');
  fs.existsSync(aegilops.resolved.blastdb.dna + '.nal').should.be.true();
});
