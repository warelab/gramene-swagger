'use strict';

require('../../../api/helpers/primers/node_compat');

const { test } = require('node:test');
const assert = require('node:assert');
const should = require('should');
const fs = require('fs');
const path = require('path');

const genomes = require('../../../api/helpers/primers/genomes');

// Fixtures are real docs: taxonomy/maps from http://localhost:10069/v69 (rice, barley) and
// http://localhost:50011/sorghum_v11 (all 128 genome maps, the whole 223-doc taxonomy).
const FIXTURES = path.join(__dirname, '..', 'fixtures', 'catalog');
function load(name) { return JSON.parse(fs.readFileSync(path.join(FIXTURES, name), 'utf8')); }
const V69_TAXONOMY = load('taxonomy_v69_rice_barley.json');
const V69_MAPS = load('maps_v69_subset.json');
const V11_TAXONOMY = load('taxonomy_sorghum_v11.json');
const V11_MAPS = load('maps_sorghum_v11.json');
const SILENT = { warn: function () {}, error: function () {}, log: function () {} };

async function rejectsWith(promise, status, code) {
  let caught = null;
  await assert.rejects(promise, function (err) { caught = err; return true; });
  should(caught.status).equal(status);
  should(caught.code).equal(code);
  return caught;
}

function clone(x) { return JSON.parse(JSON.stringify(x)); }

// Minimal stand-in for require('gramene-mongodb-config'): {maps, taxonomy}.mongoCollection().
// mode: ok | undefined (connection failed) | throw | hang | empty
function fakeMongo(mapDocs, taxonomyDocs, mode) {
  const state = { mode: mode || 'ok', calls: { maps: [], taxonomy: [] } };
  function collection(name, docs) {
    return {
      find: function (query, options) {
        state.calls[name].push({ query: clone(query), options: clone(options) });
        if (state.mode === 'throw') throw new Error('connection to localhost:27017 lost');
        let out = state.mode === 'empty' && name === 'maps' ? [] : docs;
        if (query.type !== undefined) out = out.filter(function (d) { return d.type === query.type; });
        if (query._id && query._id.$in) {
          const ids = new Set(query._id.$in);
          out = out.filter(function (d) { return ids.has(d._id); });
        }
        return {
          toArray: function () {
            if (state.mode === 'hang') return new Promise(function () {});
            return Promise.resolve(clone(out));
          }
        };
      }
    };
  }
  state.maps = { mongoCollection: async function () { return state.mode === 'undefined' ? undefined : collection('maps', mapDocs); } };
  state.taxonomy = { mongoCollection: async function () { return state.mode === 'undefined' ? undefined : collection('taxonomy', taxonomyDocs); } };
  return state;
}

async function flush() {
  for (let i = 0; i < 5; i++) await new Promise(function (r) { setImmediate(r); });
}

test('rice genomes under 39946, 1736656 and 1736659 all belong to Oryza sativa (4530)', function () {
  [39946001, 1736656001, 1736659001].forEach(function (id) {
    should(genomes.speciesOf(id, V69_TAXONOMY)).eql({ taxon_id: 4530, name: 'Oryza sativa' });
  });
  // taxon_id/1000 would split rice into three "species", none of them 4530.
  [39946001, 1736656001, 1736659001].map(function (id) { return Math.floor(id / 1000); })
    .should.eql([39946, 1736656, 1736659]);
  // 1736656 sits under Japonica Group 39947, 1736659 under Indica Group 39946: both 'no rank'.
  const idx = new Map(V69_TAXONOMY.map(function (d) { return [d._id, d]; }));
  [39946, 39947, 1736656, 1736659].forEach(function (id) { idx.get(id).rank.should.equal('no rank'); });
});

test('barley 112509001 hangs under subspecies 112509 and belongs to Hordeum vulgare (4513)', function () {
  should(genomes.speciesOf(112509001, V69_TAXONOMY)).eql({ taxon_id: 4513, name: 'Hordeum vulgare' });
  V69_TAXONOMY.find(function (d) { return d._id === 112509; }).rank.should.equal('subspecies');
});

test('sorghum 4558006 belongs to Sorghum bicolor (4558), the first species-ranked ancestor', function () {
  should(genomes.speciesOf(4558006, V11_TAXONOMY)).eql({ taxon_id: 4558, name: 'Sorghum bicolor' });
  // 'Sorghum pan' (455800000) is also ranked species, but it comes after 4558 in the ancestors list.
  const pan = V11_TAXONOMY.find(function (d) { return d._id === 455800000; });
  pan.rank.should.equal('species');
  const leaf = V11_TAXONOMY.find(function (d) { return d._id === 4558006; });
  leaf.ancestors.slice(0, 3).should.eql([4558006, 4558, 455800000]);
});

test('speciesOf: species taxa, string ids, Map indexes, unknown taxa and lineages without a species', function () {
  should(genomes.speciesOf(4530, V69_TAXONOMY)).eql({ taxon_id: 4530, name: 'Oryza sativa' });
  should(genomes.speciesOf('39946001', V69_TAXONOMY)).eql({ taxon_id: 4530, name: 'Oryza sativa' });
  const index = new Map(V69_TAXONOMY.map(function (d) { return [d._id, d]; }));
  should(genomes.speciesOf(1736659001, index)).eql({ taxon_id: 4530, name: 'Oryza sativa' });
  should(genomes.speciesOf(987654321, V69_TAXONOMY)).be.null();
  should(genomes.speciesOf(4527, V69_TAXONOMY)).be.null(); // genus Oryza: nothing species-ranked above it
  should(genomes.speciesOf(5, [{ _id: 5, rank: 'genome' }])).be.null(); // no ancestors field
  should(genomes.speciesOf(5, [{ _id: 5, rank: 'species', name: 'X y' }])).eql({ taxon_id: 5, name: 'X y' });
  should(genomes.speciesOf(10, [{ _id: 10, rank: 'genome', ancestors: [10, 11, 12] }, { _id: 12, rank: 'species', name: 'Z' }]))
    .eql({ taxon_id: 12, name: 'Z' }); // missing intermediate docs are skipped
});

test('buildCatalog on the 128 sorghum_v11 maps: 120 Sorghum bicolor genomes, sorted by display_name, frozen', function () {
  const catalog = genomes.buildCatalog(V11_MAPS, V11_TAXONOMY, { loaded_at: '2026-09-12T00:00:00.000Z' });
  catalog.loaded_at.should.equal('2026-09-12T00:00:00.000Z');
  catalog.genomes.length.should.equal(128);
  catalog.bySpecies.get('species:4558').length.should.equal(120);
  catalog.bySpecies.get('species:4530').map(function (g) { return g.system_name; }).should.eql(['oryza_sativa']);
  catalog.bySpecies.size.should.equal(9);

  const bicolor = catalog.bySystemName.get('sorghum_bicolor');
  bicolor.should.eql({
    system_name: 'sorghum_bicolor',
    display_name: 'Sb bicolor BTx623 v3',
    taxon_id: 4558006,
    map_id: 'GCA_000003195.3',
    regions: {
      names: ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', 'UNANCHORED'],
      lengths: [80884392, 77742459, 74386277, 68658214, 71854669, 61277060, 65505356, 62686529, 59416394, 61233695, 25090273]
    },
    species: { taxon_id: 4558, name: 'Sorghum bicolor' }
  });
  for (let i = 1; i < catalog.genomes.length; i++) {
    (genomes.compareGenomes(catalog.genomes[i - 1], catalog.genomes[i]) <= 0).should.be.true();
  }
  [catalog, catalog.genomes, bicolor, bicolor.regions, bicolor.regions.names, bicolor.species,
    catalog.bySpecies.get('species:4558')].forEach(function (o) { Object.isFrozen(o).should.be.true(); });
});

test('buildCatalog on v69 maps groups the three rice genomes into Oryza sativa and barley into Hordeum vulgare', function () {
  const catalog = genomes.buildCatalog(V69_MAPS, V69_TAXONOMY);
  catalog.bySpecies.get('species:4530').map(function (g) { return g.system_name; }).sort()
    .should.eql(['oryza_indica', 'oryza_sativa_azucena', 'oryza_sativa_n22']);
  catalog.bySystemName.get('hordeum_vulgare').species.should.eql({ taxon_id: 4513, name: 'Hordeum vulgare' });
  // aegilops_tauschii and goldenpromise have no taxonomy docs in this fixture: each is its own group.
  should(catalog.bySystemName.get('aegilops_tauschii').species).be.null();
  catalog.bySpecies.get('genome:aegilops_tauschii').length.should.equal(1);
});

test('buildCatalog skips non-genome maps, invalid system names and duplicates (smallest _id wins)', function () {
  const base = { taxon_id: 4558006, regions: { names: ['1'], lengths: [10] } };
  const catalog = genomes.buildCatalog([
    Object.assign({ _id: 'B', system_name: 'sorghum_dup', display_name: 'second' }, base),
    Object.assign({ _id: 'A', system_name: 'sorghum_dup', display_name: 'first', type: 'genome' }, base),
    Object.assign({ _id: 'C', system_name: 'sorghum_gm', type: 'genetic' }, base),
    Object.assign({ _id: 'D', system_name: '../etc' }, base),
    Object.assign({ _id: 'E', system_name: 'Sorghum_Upper' }, base),
    Object.assign({ system_name: 'no_id' }, base),
    { _id: 'F', system_name: 'no_regions' },
    null
  ], V11_TAXONOMY);
  catalog.genomes.map(function (g) { return g.system_name; }).sort().should.eql(['no_regions', 'sorghum_dup']);
  catalog.bySystemName.get('sorghum_dup').display_name.should.equal('first');
  catalog.bySystemName.get('no_regions').should.eql({
    system_name: 'no_regions', display_name: 'no_regions', taxon_id: null, map_id: 'F',
    regions: { names: [], lengths: [] }, species: null
  });
});

test('isValidSystemName enforces ^[a-z0-9_]+$ and at most 128 characters', function () {
  ['sorghum_bicolor', 'a', '0', 'x'.repeat(128)].forEach(function (n) { genomes.isValidSystemName(n).should.be.true(); });
  ['', 'x'.repeat(129), '../etc', 'a/b', 'Sorghum', 'sorghum-x', 'sorghum.x', 'sorghum_x\n', ' x', null, 7, {}, ['a']]
    .forEach(function (n) { genomes.isValidSystemName(n).should.be.false(); });
});

test('getCatalog loads genome maps and only the needed taxonomy docs, then serves the cache for 10 minutes', async function () {
  genomes.clearCache();
  const extra = { _id: 'gm1', system_name: 'sorghum_genetic', type: 'genetic_map', taxon_id: 4558006 };
  const mongo = fakeMongo(V11_MAPS.concat([extra]), V11_TAXONOMY);
  let now = 5e6;
  const deps = { mongo: mongo, now: function () { return now; }, log: SILENT };

  const first = await Promise.all([genomes.getCatalog(deps), genomes.getCatalog(deps)]);
  first[0].should.equal(first[1]); // concurrent callers share one load
  const c1 = first[0];
  c1.genomes.length.should.equal(128);
  mongo.calls.maps.length.should.equal(1);
  mongo.calls.maps[0].query.should.eql({ type: 'genome' });
  mongo.calls.maps[0].options.fields.should.have.properties(['_id', 'system_name', 'display_name', 'taxon_id', 'regions']);
  mongo.calls.taxonomy.length.should.equal(2);
  const leafIds = mongo.calls.taxonomy[0].query._id.$in;
  leafIds.length.should.equal(128);
  const ancestorIds = mongo.calls.taxonomy[1].query._id.$in;
  ancestorIds.should.containEql(4558);
  ancestorIds.should.not.containEql(4558006);

  now += genomes.CATALOG_TTL_MS - 1;
  (await genomes.getCatalog(deps)).should.equal(c1);
  mongo.calls.maps.length.should.equal(1);

  now += 1; // expired: the previous catalog is served while a reload runs in the background
  (await genomes.getCatalog(deps)).should.equal(c1);
  await flush();
  mongo.calls.maps.length.should.equal(2);
  const c2 = await genomes.getCatalog(deps);
  c2.should.not.equal(c1);
  c2.genomes.length.should.equal(128);
});

test('getCatalog: mongo unavailable with no cached catalog is 503 MONGO_UNAVAILABLE; it recovers once mongo is back', async function () {
  genomes.clearCache();
  const mongo = fakeMongo(V11_MAPS, V11_TAXONOMY, 'undefined');
  const deps = { mongo: mongo, log: SILENT, mongo_timeout_ms: 50 };

  const err = await rejectsWith(genomes.getCatalog(deps), 503, 'MONGO_UNAVAILABLE');
  err.details.retry_after_s.should.equal(30);
  err.message.should.not.match(/27017|localhost/);

  mongo.mode = 'throw';
  await rejectsWith(genomes.getCatalog(deps), 503, 'MONGO_UNAVAILABLE');
  mongo.mode = 'empty';
  await rejectsWith(genomes.getCatalog(deps), 503, 'MONGO_UNAVAILABLE');
  mongo.mode = 'hang';
  const keepAlive = setTimeout(function () {}, 5000); // the query timeout timer is unref'd
  try {
    await rejectsWith(genomes.getCatalog(deps), 503, 'MONGO_UNAVAILABLE');
  } finally {
    clearTimeout(keepAlive);
  }
  await rejectsWith(genomes.byName('sorghum_bicolor', deps), 503, 'MONGO_UNAVAILABLE');
  await rejectsWith(genomes.sameSpecies('sorghum_bicolor', deps), 503, 'MONGO_UNAVAILABLE');

  mongo.mode = 'ok';
  (await genomes.getCatalog(deps)).genomes.length.should.equal(128);
});

test('getCatalog: a failed reload keeps serving the previous catalog and retries at most every 30 s', async function () {
  genomes.clearCache();
  const mongo = fakeMongo(V11_MAPS, V11_TAXONOMY);
  const warnings = [];
  let now = 1e7;
  const deps = { mongo: mongo, now: function () { return now; }, log: { warn: function (m) { warnings.push(m); } } };
  const c1 = await genomes.getCatalog(deps);

  mongo.mode = 'throw';
  now += genomes.CATALOG_TTL_MS;
  (await genomes.getCatalog(deps)).should.equal(c1);
  await flush();
  mongo.calls.maps.length.should.equal(2);
  warnings.length.should.equal(1);
  warnings[0].should.match(/MONGO_UNAVAILABLE/);

  now += 29999;
  (await genomes.getCatalog(deps)).should.equal(c1);
  await flush();
  mongo.calls.maps.length.should.equal(2); // throttled

  mongo.mode = 'ok';
  now += 1;
  (await genomes.getCatalog(deps)).should.equal(c1);
  await flush();
  mongo.calls.maps.length.should.equal(3);
  (await genomes.getCatalog(deps)).should.not.equal(c1);
});

test('deps.catalog bypasses mongo entirely', async function () {
  const catalog = genomes.buildCatalog(V11_MAPS, V11_TAXONOMY);
  const exploding = {};
  Object.defineProperty(exploding, 'maps', { get: function () { throw new Error('mongo must not be used'); } });
  (await genomes.getCatalog({ catalog: catalog, mongo: exploding })).should.equal(catalog);
});

test('byName, requireGenome and sameSpecies', async function () {
  const catalog = genomes.buildCatalog(V11_MAPS.concat([
    { _id: 'X1', system_name: 'orphan_genome', display_name: 'Orphan', taxon_id: 123456789, type: 'genome', regions: {} }
  ]), V11_TAXONOMY);
  const deps = { catalog: catalog };

  (await genomes.byName('sorghum_rio', deps)).map_id.should.equal('GCA_015952705.1');
  should(await genomes.byName('sorghum_nope', deps)).be.null();
  should(await genomes.byName('../etc', deps)).be.null();
  (await genomes.requireGenome('sorghum_rio', deps)).system_name.should.equal('sorghum_rio');
  const unknown = await rejectsWith(genomes.requireGenome('sorghum_nope', deps), 404, 'UNKNOWN_GENOME');
  unknown.details.should.eql({ system_name: 'sorghum_nope' });
  const invalid = await rejectsWith(genomes.sameSpecies('../etc', deps), 404, 'UNKNOWN_GENOME');
  invalid.details.should.eql({});
  invalid.message.should.not.containEql('etc');

  const sorghum = await genomes.sameSpecies('sorghum_bicolor', deps);
  sorghum.length.should.equal(120);
  sorghum.map(function (g) { return g.system_name; }).should.containEql('sorghum_bicolor');
  sorghum.every(function (g) { return g.species.taxon_id === 4558; }).should.be.true();
  (await genomes.sameSpecies('zea_maysb73', deps)).map(function (g) { return g.system_name; }).should.eql(['zea_maysb73']);
  (await genomes.sameSpecies('orphan_genome', deps)).map(function (g) { return g.system_name; }).should.eql(['orphan_genome']);
  await rejectsWith(genomes.sameSpecies('sorghum_nope', deps), 404, 'UNKNOWN_GENOME');
});
