'use strict';

// Deterministic job ids (spec §A.2.3): sha256(canonicalJSON({v: 1, algo, request: normalized, dbs})).slice(0, 32).

const test = require('node:test');
const should = require('should');
const crypto = require('crypto');

const jobs = require('../../../api/helpers/primers/jobs');
const check = require('../../../api/helpers/primers/check');
const { normalize } = require('../../../api/helpers/primers/check/normalize');
const { createMemoryStore } = require('../../../api/helpers/primers/jobs/memory_store');
const config = require('../../../api/helpers/primers/config');
const genomes = require('../../../api/helpers/primers/genomes');
const mapsFixture = require('../fixtures/catalog/maps_sorghum_v11.json');
const taxonomyFixture = require('../fixtures/catalog/taxonomy_sorghum_v11.json');

const CATALOG = genomes.buildCatalog(mapsFixture, taxonomyFixture);
const CFG = config._build({ env: {}, fileConfig: {} }).config;
const quiet = { info() {}, warn() {}, error() {}, log() {} };

function fp(name, salt) {
  return crypto.createHash('sha1').update((salt || 'fp:') + name).digest('hex');
}

function stubAssemblies(fingerprints) {
  fingerprints = fingerprints || {};
  function build(name) {
    const g = CATALOG.bySystemName.get(name);
    if (!g) throw genomes.unknownGenomeError(name);
    return {
      system_name: name, taxon_id: g.taxon_id, display_name: g.display_name, map_id: g.map_id,
      prefix: 'P', dir: '/scratch/olson/fasta/' + name,
      fasta: { dna: '/d/' + name + '.fa.gz', dna_sm: null },
      blastdb: { dna: '/d/' + name + '.dna.toplevel', cdna: '/d/' + name + '.cdna.all' },
      repeat_masking: 'absent', total_bases: 7e8, num_sequences: 10,
      fingerprint: fingerprints[name] || fp(name), warnings: []
    };
  }
  return {
    async resolve(name) { return build(name); },
    async resolveMany(names) { return names.map(function (n) { return { system_name: n, resolved: build(n), error: null }; }); }
  };
}

const GENE = {
  _id: 'SORBI_3004G087700',
  system_name: 'sorghum_bicolor',
  location: { region: '4', start: 7421357, end: 7428285, strand: 1, map: 'GCA_000003195.3' },
  gene_structure: { canonical_transcript: 'SORBI_3004G087700.3', transcripts: [{ id: 'SORBI_3004G087700.3' }, { id: 'SORBI_3004G087700.1' }] },
  homology: { homologous_genes: { ortholog_one2one: ['353.004G093500'] } }
};

function stubMongo() {
  return {
    genes: {
      mongoCollection: async function () {
        return {
          find(query) {
            const out = typeof query._id === 'string'
              ? (query._id === GENE._id ? [GENE] : [])
              : (query._id.$in.indexOf('353.004G093500') >= 0 ? [{ _id: '353.004G093500', system_name: 'sorghum_353' }] : []);
            return { toArray: async function () { return JSON.parse(JSON.stringify(out)); } };
          }
        };
      }
    }
  };
}

function deps(extra) {
  return Object.assign({ cfg: CFG, catalog: CATALOG, assemblies: stubAssemblies(), mongo: stubMongo(), log: quiet }, extra || {});
}

async function idOf(body, d) {
  const norm = await normalize(body, d || deps());
  return jobs.jobId(norm.request, norm.dbs, check.ALGORITHM_VERSION);
}

const L = 'GGACAGCTCCACAACATATCAG';
const R = 'GGACATTTGAAGCCCATGGCC';
const L3 = 'GATATCAGTGGAATCATAAGACCG';
const R3 = 'CATCGATATCAGGATCTGGCTT';

test('canonicalJSON sorts object keys at every level, keeps array order and skips undefined', function () {
  jobs.canonicalJSON({ b: 1, a: { d: [3, { z: 1, y: 2 }], c: undefined }, e: null })
    .should.equal('{"a":{"d":[3,{"y":2,"z":1}]},"b":1,"e":null}');
  jobs.canonicalJSON([undefined, NaN, 'x']).should.equal('[null,null,"x"]');
  jobs.canonicalJSON({ b: 1, a: 2 }).should.equal(jobs.canonicalJSON({ a: 2, b: 1 }));
});

test('jobId is 32 hex characters and follows the documented formula', function () {
  const request = { system_name: 'sorghum_bicolor', pairs: [] };
  const dbs = { sorghum_bicolor: 'abc' };
  const id = jobs.jobId(request, dbs, '2');
  id.should.match(jobs.JOB_ID_RE);
  const expected = crypto.createHash('sha256')
    .update(jobs.canonicalJSON({ v: 1, algo: '2', request: request, dbs: dbs })).digest('hex').slice(0, 32);
  id.should.equal(expected);
  jobs.jobId(request, dbs).should.equal(id); // default algo = check.ALGORITHM_VERSION ('2')
  // '2': max_amplifying_mismatches stringency and discarded overlapping footprints; version-1 jobs are never reused.
  check.ALGORITHM_VERSION.should.equal('2');
  jobs.jobId(request, dbs, '1').should.not.equal(id);
});

test('max_amplifying_mismatches is part of the normalized request: default 3 fills in, other values change the id', async function () {
  const body = { system_name: 'sorghum_bicolor', pairs: [{ id: 'P2', left: L, right: R }] };
  const base = await idOf(body);
  (await normalize(body, deps())).request.params.max_amplifying_mismatches.should.equal(3);
  (await idOf(Object.assign({}, body, { params: { max_amplifying_mismatches: 3 } }))).should.equal(base);
  (await idOf(Object.assign({}, body, { params: { max_amplifying_mismatches: 2 } }))).should.not.equal(base);
  // an omitted cap is lowered below an explicit ignore_mismatches 3; that equals asking for cap 2
  (await idOf(Object.assign({}, body, { params: { ignore_mismatches: 3 } })))
    .should.equal(await idOf(Object.assign({}, body, { params: { ignore_mismatches: 3, max_amplifying_mismatches: 2 } })));
});

test('body key order does not change the id', async function () {
  const a = await idOf({
    system_name: 'sorghum_bicolor', mode: 'gene', gene_id: 'SORBI_3004G087700', checks: ['specificity'],
    params: { ignore_mismatches: 5, max_product_size: 3000 },
    pairs: [{ id: 'P2', left: L, right: R, expected: { region: '4', start: 7423537, end: 7423746 } }]
  });
  const b = await idOf({
    pairs: [{ expected: { end: 7423746, start: 7423537, region: '4' }, right: R, left: L, id: 'P2' }],
    params: { max_product_size: 3000, ignore_mismatches: 5 },
    checks: ['specificity'], gene_id: 'SORBI_3004G087700', mode: 'gene', system_name: 'sorghum_bicolor'
  });
  a.should.equal(b);
});

test('default filling does not change the id (omitted vs explicit defaults, case, checks and genomes order)', async function () {
  const minimal = await idOf({ system_name: 'sorghum_bicolor', pairs: [{ id: 'P2', left: L, right: R }] });
  const explicit = await idOf({
    system_name: 'sorghum_bicolor', mode: 'region', checks: ['specificity', 'specificity'],
    params: { max_product_size: 4000, ignore_mismatches: 6, min_total_mismatches: 2, min_3p_mismatches: 2,
      three_prime_window: 5, include_unlikely: false, repeat_site_threshold: 5 },
    pairs: [{ id: 'P2', left: L.toLowerCase(), right: R }]
  });
  minimal.should.equal(explicit);

  const pan1 = await idOf({ system_name: 'sorghum_bicolor', checks: ['specificity', 'pangenome'], genomes: ['sorghum_rio', 'sorghum_353'],
    pairs: [{ id: 'P2', left: L, right: R }] });
  const pan2 = await idOf({ system_name: 'sorghum_bicolor', checks: ['pangenome'], genomes: ['sorghum_353', 'sorghum_rio', 'sorghum_353', 'sorghum_bicolor'],
    pairs: [{ id: 'P2', left: L, right: R }] });
  pan1.should.equal(pan2);

  // transcript mode: omitted transcript_id defaults to the canonical one; a stripped expected leaves no trace.
  const tx1 = await idOf({ system_name: 'sorghum_bicolor', mode: 'transcript', gene_id: 'SORBI_3004G087700',
    pairs: [{ id: 'Q', left: L, right: R }] });
  const tx2 = await idOf({ system_name: 'sorghum_bicolor', mode: 'transcript', gene_id: 'SORBI_3004G087700', transcript_id: 'SORBI_3004G087700.3',
    pairs: [{ id: 'Q', left: L, right: R, expected: { region: '4', start: 1, end: 100 } }] });
  tx1.should.equal(tx2);

  // a max_product_size raised to cover an expected product equals asking for it explicitly
  const raised = await idOf({ system_name: 'sorghum_bicolor', pairs: [{ id: 'P', left: L, right: R, expected: { region: '4', start: 1, end: 5000 } }] });
  const asked = await idOf({ system_name: 'sorghum_bicolor', params: { max_product_size: 6000 },
    pairs: [{ id: 'P', left: L, right: R, expected: { region: '4', start: 1, end: 5000 } }] });
  raised.should.equal(asked);
});

test('fingerprints, algorithm version and request content change the id', async function () {
  const body = { system_name: 'sorghum_bicolor', checks: ['pangenome'], genomes: ['sorghum_353'], pairs: [{ id: 'P2', left: L, right: R }] };
  const base = await idOf(body);
  (await idOf(body, deps({ assemblies: stubAssemblies({ sorghum_bicolor: fp('sorghum_bicolor', 'rebuilt:') }) }))).should.not.equal(base);
  (await idOf(body, deps({ assemblies: stubAssemblies({ sorghum_353: fp('sorghum_353', 'rebuilt:') }) }))).should.not.equal(base);

  const norm = await normalize(body, deps());
  jobs.jobId(norm.request, norm.dbs, '1').should.not.equal(base);
  jobs.jobId(norm.request, norm.dbs, '2').should.equal(base);

  (await idOf(Object.assign({}, body, { genomes: ['sorghum_rio'] }))).should.not.equal(base);
  (await idOf(Object.assign({}, body, { params: { include_unlikely: true } }))).should.not.equal(base);
  (await idOf(Object.assign({}, body, { pairs: [{ id: 'P2b', left: L, right: R }] }))).should.not.equal(base);
  (await idOf(Object.assign({}, body, { mode: 'sequence' }))).should.not.equal(base);

  const two = [{ id: 'P2', left: L, right: R }, { id: 'P3', left: L3, right: R3 }];
  (await idOf({ system_name: 'sorghum_bicolor', pairs: two }))
    .should.not.equal(await idOf({ system_name: 'sorghum_bicolor', pairs: two.slice().reverse() })); // pair order is kept
});

test('jobs.submit uses exactly that id, and identical bodies share one job', async function () {
  const store = createMemoryStore({ cfg: CFG, siteKey: 'jobid_test', shared: { slots: new Map(), panSlots: new Map() } });
  const d = deps({ store: store, check: { ALGORITHM_VERSION: check.ALGORITHM_VERSION, normalize: normalize } });
  const body = { system_name: 'sorghum_bicolor', pairs: [{ id: 'P2', left: L, right: R }] };
  const first = await jobs.submit(body, d);
  first.job_id.should.equal(await idOf(body));
  first.created.should.be.true();
  const second = await jobs.submit({ pairs: [{ right: R, left: L.toLowerCase(), id: 'P2' }], system_name: 'sorghum_bicolor', mode: 'region' }, d);
  second.job_id.should.equal(first.job_id);
  second.created.should.be.false();
});
