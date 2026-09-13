'use strict';

require('../../../api/helpers/primers/node_compat');

const { test } = require('node:test');
const assert = require('node:assert');
const should = require('should');
const fs = require('fs');
const path = require('path');

const genomes = require('../../../api/helpers/primers/genomes');
const fx = require('../fixtures/catalog/fsfixture');

const FIXTURES = path.join(__dirname, '..', 'fixtures', 'catalog');
const V11_MAPS = JSON.parse(fs.readFileSync(path.join(FIXTURES, 'maps_sorghum_v11.json'), 'utf8'));
const V11_TAXONOMY = JSON.parse(fs.readFileSync(path.join(FIXTURES, 'taxonomy_sorghum_v11.json'), 'utf8'));
const SILENT = { warn: function () {}, error: function () {}, log: function () {} };
const ENTRY_KEYS = ['display_name', 'has_blastdb', 'has_cdna_blastdb', 'has_sequence', 'is_query', 'map_id',
  'repeat_masking', 'system_name', 'taxon_id', 'total_bases', 'warnings'];

async function rejectsWith(promise, status, code) {
  let caught = null;
  await assert.rejects(promise, function (err) { caught = err; return true; });
  should(caught.status).equal(status);
  should(caught.code).equal(code);
  return caught;
}

function seqs(seed) {
  return [
    { name: '1', seq: fx.randomSeq(40000, seed) },
    { name: '2', seq: fx.randomSeq(25000, seed + 1) },
    { name: '3', seq: fx.randomSeq(8000, seed + 2) }
  ];
}

function sizeOf(list) {
  return list.reduce(function (s, x) { return s + x.seq.length; }, 0);
}

// The real 128 sorghum_v11 maps; five genomes get files under root (with map regions matching them).
function buildFixture(root) {
  const maps = JSON.parse(JSON.stringify(V11_MAPS));
  const map = function (name) { return maps.find(function (m) { return m.system_name === name; }); };
  const useSeqs = function (name, list) {
    map(name).regions = { names: list.map(function (s) { return s.name; }), lengths: list.map(function (s) { return s.seq.length; }) };
  };
  let dir;
  let P;

  // sorghum_bicolor: dna + uppercase dna_sm copy, dna and cdna BLAST DBs
  const bicolor = seqs(1);
  useSeqs('sorghum_bicolor', bicolor);
  dir = path.join(root, 'sorghum_bicolor');
  P = 'Sorghum_bicolor.Sorghum_bicolor_NCBIv3';
  fx.bgzipFasta(dir, P, 'dna', bicolor);
  fx.bgzipFasta(dir, P, 'dna_sm', bicolor);
  fx.blastDb(dir, P + '.dna.toplevel');
  fx.blastDb(dir, P + '.cdna.all');

  // sorghum_rio: a real soft-mask, dna BLAST DB only
  const rio = seqs(2);
  useSeqs('sorghum_rio', rio);
  dir = path.join(root, 'sorghum_rio');
  P = 'Sorghum_rio.JGI-v2.0';
  fx.bgzipFasta(dir, P, 'dna', rio);
  fx.bgzipFasta(dir, P, 'dna_sm', rio.map(function (s, i) {
    return i ? s : { name: s.name, seq: s.seq.slice(0, 15000) + s.seq.slice(15000, 25000).toLowerCase() + s.seq.slice(25000) };
  }));
  fx.blastDb(dir, P + '.dna.toplevel');

  // sorghum_pi534133: dual prefix; the JGI prefix (map._id suffix) is picked
  const pi = seqs(3);
  useSeqs('sorghum_pi534133', pi);
  dir = path.join(root, 'sorghum_pi534133');
  const jgi = 'Sorghum_pi534133.' + map('sorghum_pi534133')._id;
  fx.indexOnlyFasta(dir, jgi, 'dna', pi.map(function (s) { return [s.name, s.seq.length]; }));
  fx.blastDb(dir, jgi + '.dna.toplevel');
  fx.blastDb(dir, jgi + '.cdna.all');
  fx.indexOnlyFasta(dir, 'Sorghum_pi534133.pi534133', 'dna', [['1', 5]]);
  fx.blastDb(dir, 'Sorghum_pi534133.pi534133.dna.toplevel');

  // sorghum_353: two unrelated prefixes and nothing to choose by -> AMBIGUOUS_ASSEMBLY
  dir = path.join(root, 'sorghum_353');
  fx.indexOnlyFasta(dir, 'Sorghum_353.a', 'dna', [['1', 10]]);
  fx.indexOnlyFasta(dir, 'Sorghum_353.b', 'dna', [['1', 11]]);

  // zea_maysb73: another species, fully provisioned; never listed for sorghum
  const zea = seqs(4);
  useSeqs('zea_maysb73', zea);
  dir = path.join(root, 'zea_maysb73');
  P = 'Zea_mays.Zm-B73-REFERENCE-NAM-5.0';
  fx.bgzipFasta(dir, P, 'dna', zea);
  fx.blastDb(dir, P + '.dna.toplevel');
  fx.blastDb(dir, P + '.cdna.all');

  return { catalog: genomes.buildCatalog(maps, V11_TAXONOMY), sizes: { bicolor: sizeOf(bicolor), rio: sizeOf(rio), pi: sizeOf(pi) } };
}

function depsFor(root, catalog, extra) {
  return Object.assign({
    catalog: catalog,
    cfg: { fasta_root: root, assembly_overrides: {}, repeat_masking_overrides: {} },
    log: SILENT,
    cache: false
  }, extra || {});
}

function stringsIn(value, out) {
  out = out || [];
  if (typeof value === 'string') out.push(value);
  else if (Array.isArray(value)) value.forEach(function (v) { stringsIn(v, out); });
  else if (value && typeof value === 'object') Object.keys(value).forEach(function (k) { stringsIn(value[k], out); });
  return out;
}

test('genomes body: query first, then its species by display_name; flags, counts and no filesystem paths', async function (t) {
  const root = fx.makeRoot('genomes');
  t.after(function () { fx.removeRoot(root); });
  const fixture = buildFixture(root);
  const body = await genomes.genomesResponse('sorghum_rio', depsFor(root, fixture.catalog));

  Object.keys(body).sort().should.eql(['counts', 'genomes', 'species', 'system_name']);
  body.system_name.should.equal('sorghum_rio');
  body.species.should.eql({ taxon_id: 4558, name: 'Sorghum bicolor' });
  body.genomes.length.should.equal(120);
  body.genomes[0].system_name.should.equal('sorghum_rio');
  body.genomes[0].is_query.should.be.true();
  body.genomes.slice(1).every(function (g) { return g.is_query === false; }).should.be.true();
  const rest = body.genomes.slice(1).map(function (g) { return fixture.catalog.bySystemName.get(g.system_name); });
  for (let i = 1; i < rest.length; i++) (genomes.compareGenomes(rest[i - 1], rest[i]) <= 0).should.be.true();
  body.genomes.forEach(function (g) { Object.keys(g).sort().should.eql(ENTRY_KEYS); });
  body.genomes.some(function (g) { return g.system_name === 'zea_maysb73'; }).should.be.false();

  const by = {};
  body.genomes.forEach(function (g) { by[g.system_name] = g; });
  by.sorghum_rio.should.eql({
    system_name: 'sorghum_rio', display_name: 'Sb bicolor PI651496 Rio', taxon_id: 4558116, map_id: 'GCA_015952705.1',
    is_query: true, has_sequence: true, has_blastdb: true, has_cdna_blastdb: false, repeat_masking: 'soft_masked',
    total_bases: fixture.sizes.rio, warnings: []
  });
  by.sorghum_bicolor.should.eql({
    system_name: 'sorghum_bicolor', display_name: 'Sb bicolor BTx623 v3', taxon_id: 4558006, map_id: 'GCA_000003195.3',
    is_query: false, has_sequence: true, has_blastdb: true, has_cdna_blastdb: true, repeat_masking: 'unmasked_copy',
    total_bases: fixture.sizes.bicolor, warnings: []
  });
  by.sorghum_pi534133.should.have.properties({
    has_sequence: true, has_blastdb: true, has_cdna_blastdb: true, repeat_masking: 'absent', total_bases: fixture.sizes.pi, warnings: []
  });
  by.sorghum_353.should.have.properties({ has_sequence: false, has_blastdb: false, has_cdna_blastdb: false, repeat_masking: 'absent', total_bases: null });
  by.sorghum_353.warnings.length.should.equal(1);
  by.sorghum_353.warnings[0].code.should.equal('AMBIGUOUS_ASSEMBLY');
  by.sorghum_tx436pac.should.have.properties({ has_sequence: false, has_blastdb: false, has_cdna_blastdb: false, repeat_masking: 'absent', total_bases: null, warnings: [] });

  body.counts.should.eql({ total: 120, with_blastdb: 3, with_cdna_blastdb: 2 });

  const json = JSON.stringify(body);
  json.should.not.containEql(root);
  json.should.not.containEql('.toplevel');
  json.should.not.containEql('Sorghum_rio.JGI');
  stringsIn(body).forEach(function (s) { s.charAt(0).should.not.equal('/'); });
});

test('a different query genome is listed first and the rest keep display_name order', async function (t) {
  const root = fx.makeRoot('genomes');
  t.after(function () { fx.removeRoot(root); });
  const fixture = buildFixture(root);
  const body = await genomes.genomesResponse('sorghum_bicolor', depsFor(root, fixture.catalog));
  body.genomes[0].system_name.should.equal('sorghum_bicolor');
  body.genomes.filter(function (g) { return g.is_query; }).length.should.equal(1);
  const names = body.genomes.map(function (g) { return g.system_name; });
  names.indexOf('sorghum_rio').should.be.above(0);
  names.slice(1).should.eql(fixture.catalog.bySpecies.get('species:4558')
    .map(function (g) { return g.system_name; }).filter(function (n) { return n !== 'sorghum_bicolor'; }));

  const zea = await genomes.genomesResponse('zea_maysb73', depsFor(root, fixture.catalog));
  // The real sorghum_v11 taxonomy doc for 4577 is named 'Zea maysB73'; the name comes from taxonomy as-is.
  const zeaSpecies = V11_TAXONOMY.find(function (d) { return d._id === 4577; });
  zeaSpecies.rank.should.equal('species');
  zea.species.should.eql({ taxon_id: 4577, name: zeaSpecies.name });
  zea.genomes.map(function (g) { return g.system_name; }).should.eql(['zea_maysb73']);
  zea.counts.should.eql({ total: 1, with_blastdb: 1, with_cdna_blastdb: 1 });
});

test('unknown or invalid system_name is 404 UNKNOWN_GENOME; invalid names never reach mongo or the filesystem', async function (t) {
  const root = fx.makeRoot('genomes');
  t.after(function () { fx.removeRoot(root); });
  const fixture = buildFixture(root);
  const err = await rejectsWith(genomes.genomesResponse('sorghum_nope', depsFor(root, fixture.catalog)), 404, 'UNKNOWN_GENOME');
  err.details.should.eql({ system_name: 'sorghum_nope' });

  const exploding = {};
  Object.defineProperty(exploding, 'maps', { get: function () { throw new Error('mongo must not be used'); } });
  const spy = fx.spyFs();
  for (const bad of ['../etc', 'Sorghum_bicolor', '', 'x'.repeat(129), null, 3]) {
    await rejectsWith(genomes.genomesResponse(bad, { mongo: exploding, fs: spy, cfg: { fasta_root: root }, log: SILENT }), 404, 'UNKNOWN_GENOME');
  }
  spy.calls.should.eql([]);
});

test('mongo unavailable is 503 MONGO_UNAVAILABLE', async function () {
  genomes.clearCache();
  const mongo = {
    maps: { mongoCollection: async function () { return undefined; } },
    taxonomy: { mongoCollection: async function () { return undefined; } }
  };
  const err = await rejectsWith(genomes.genomesResponse('sorghum_bicolor', { mongo: mongo, cfg: { fasta_root: '/nonexistent' }, log: SILENT }),
    503, 'MONGO_UNAVAILABLE');
  err.details.retry_after_s.should.be.above(0);
});

test('unexpected resolver failures become a per-genome INTERNAL warning without paths', async function (t) {
  const root = fx.makeRoot('genomes');
  t.after(function () { fx.removeRoot(root); });
  const fixture = buildFixture(root);
  const logged = [];
  const denied = path.join(root, 'sorghum_bicolor'); // exact dir: sorghum_bicolorv5 etc. must not match
  const spy = fx.spyFs({
    readdir: function (p) {
      if (p === denied || p.indexOf(denied + path.sep) === 0) {
        const e = new Error('EACCES: permission denied, scandir \'' + p + '\'');
        e.code = 'EACCES';
        return Promise.reject(e);
      }
      return fs.promises.readdir(p);
    }
  });
  const body = await genomes.genomesResponse('sorghum_rio', depsFor(root, fixture.catalog, { fs: spy, log: { error: function (m) { logged.push(m); } } }));
  const bicolor = body.genomes.find(function (g) { return g.system_name === 'sorghum_bicolor'; });
  bicolor.should.have.properties({ has_sequence: false, has_blastdb: false, has_cdna_blastdb: false, repeat_masking: 'absent', total_bases: null });
  bicolor.warnings.should.eql([{ code: 'INTERNAL', message: 'the assembly files for this genome could not be inspected' }]);
  JSON.stringify(body).should.not.containEql(root);
  logged.length.should.equal(1);
  body.counts.should.eql({ total: 120, with_blastdb: 2, with_cdna_blastdb: 1 });
});

test('resolveMany is called once, query first, with the loaded catalog; species-less genomes list only themselves', async function () {
  const catalog = genomes.buildCatalog(V11_MAPS.concat([
    { _id: 'O1', system_name: 'orphan_genome', display_name: 'Orphan', taxon_id: 123456789, type: 'genome', regions: {} }
  ]), V11_TAXONOMY);
  const calls = [];
  const stub = {
    resolveMany: async function (names, deps) {
      calls.push({ names: names, catalog: deps.catalog });
      return names.map(function (n) {
        return {
          system_name: n, error: null,
          resolved: { fasta: { dna: '/x', dna_sm: null }, blastdb: { dna: '/y', cdna: null }, repeat_masking: 'unmasked_copy', total_bases: 7, warnings: [{ code: 'ASSEMBLY_MISMATCH', message: 'see /scratch/olson/fasta/x' }] }
        };
      });
    }
  };
  const body = await genomes.genomesResponse('sorghum_353', { catalog: catalog, assemblies: stub, log: SILENT });
  calls.length.should.equal(1);
  calls[0].catalog.should.equal(catalog);
  calls[0].names[0].should.equal('sorghum_353');
  calls[0].names.length.should.equal(120);
  body.genomes[0].warnings.should.eql([{ code: 'ASSEMBLY_MISMATCH', message: 'see x' }]);

  const orphan = await genomes.genomesResponse('orphan_genome', { catalog: catalog, assemblies: stub, log: SILENT });
  should(orphan.species).be.null();
  orphan.genomes.map(function (g) { return g.system_name; }).should.eql(['orphan_genome']);
  orphan.counts.should.eql({ total: 1, with_blastdb: 1, with_cdna_blastdb: 0 });
});
