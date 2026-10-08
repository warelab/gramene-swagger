'use strict';

// Records the compatibility snapshots of the genotyping spec §7.2 from the current code:
//   C1 fixtures/contract/compat_job_ids.json           jobs.jobId of every existing check-*.json request fixture
//   C2 fixtures/contract/compat_design_normalize.json  design.normalize of every existing design-*.json request fixture
//   C3 fixtures/contract/compat_extract_pairs.json     boulder.extractPairs of every recorded Primer3 output
// First recorded before any genotyping code, on code identical to primer-design d4af3a1 (genotyping 0a743b7 only
// changes jobs/worker_main.js). Re-record only when a change to the pinned behaviour is intended:
//   node test/primers/tools/record_compat.js           reports which snapshots differ from the code; writes nothing
//   node test/primers/tools/record_compat.js --write   rewrites the snapshots that differ
// Not a test: the unit suite globs test/primers/unit only. compat_ids.test.js, compat_design_normalize.test.js and
// boulder.test.js take the fixture lists, config and offline stubs from here, so they normalize exactly as recorded.

const fs = require('fs');
const path = require('path');
const crypto = require('crypto');

const design = require('../../../api/helpers/primers/design');
const boulder = require('../../../api/helpers/primers/boulder');
const config = require('../../../api/helpers/primers/config');
const genomes = require('../../../api/helpers/primers/genomes');
const jobs = require('../../../api/helpers/primers/jobs');
const check = require('../../../api/helpers/primers/check');
const { normalize } = require('../../../api/helpers/primers/check/normalize');
const mapsFixture = require('../fixtures/catalog/maps_sorghum_v11.json');
const taxonomyFixture = require('../fixtures/catalog/taxonomy_sorghum_v11.json');

const ROOT = path.resolve(__dirname, '..', '..', '..');
const FIXTURES = path.join(ROOT, 'test', 'primers', 'fixtures');
const REQUESTS_DIR = path.join(FIXTURES, 'contract', 'requests');
const PRIMER3_DIR = path.join(FIXTURES, 'primer3');
const SNAPSHOTS = Object.freeze({
  jobIds: path.join(FIXTURES, 'contract', 'compat_job_ids.json'),
  designNormalize: path.join(FIXTURES, 'contract', 'compat_design_normalize.json'),
  extractPairs: path.join(FIXTURES, 'contract', 'compat_extract_pairs.json')
});

// Built-in defaults, no env or config file (contract.test.js cfgForTests).
const CFG = config._build({ env: {}, fileConfig: {} }).config;
const CATALOG = genomes.buildCatalog(mapsFixture, taxonomyFixture);
const quiet = { info() {}, warn() {}, error() {}, log() {} };

// ---- offline check stubs (as in normalize.test.js) ------------------------------------------------------

function fp(name) {
  return crypto.createHash('sha1').update('fp:' + name).digest('hex');
}

// Every catalog genome resolves with sequence, dna and cdna BLAST DBs, 0.7 Gb and fingerprint sha1('fp:' + name).
function stubAssemblies() {
  function build(name) {
    const g = CATALOG.bySystemName.get(name);
    if (!g) throw genomes.unknownGenomeError(name);
    const dir = '/scratch/olson/fasta/' + name;
    return Object.freeze({
      system_name: name,
      taxon_id: g.taxon_id,
      display_name: g.display_name,
      map_id: g.map_id,
      prefix: 'Prefix_' + name,
      dir: dir,
      fasta: { dna: dir + '/dna/Prefix.dna.toplevel.fa.gz', dna_sm: null },
      blastdb: { dna: dir + '/Prefix.dna.toplevel', cdna: dir + '/Prefix.cdna.all' },
      repeat_masking: 'absent',
      total_bases: 700000000,
      num_sequences: 12,
      fingerprint: fp(name),
      warnings: []
    });
  }
  return {
    async resolve(name) {
      return build(name);
    },
    async resolveMany(names) {
      return names.map(function (n) {
        try {
          return { system_name: n, resolved: build(n), error: null };
        } catch (e) {
          return { system_name: n, resolved: null, error: e };
        }
      });
    }
  };
}

// The one gene the check fixtures name (SORBI_3004G087700, canonical transcript .3) and its homologues.
const GENE_DOCS = [
  {
    _id: 'SORBI_3004G087700',
    system_name: 'sorghum_bicolor',
    taxon_id: 4558006,
    location: { region: '4', start: 7421357, end: 7428285, strand: 1, map: 'GCA_000003195.3' },
    gene_structure: {
      canonical_transcript: 'SORBI_3004G087700.3',
      transcripts: [{ id: 'SORBI_3004G087700.3' }, { id: 'SORBI_3004G087700.1' }, { id: 'SORBI_3004G087700.2' }]
    },
    homology: {
      homologous_genes: {
        ortholog_one2one: ['353.004G093500', 'IS12661.004G086400'],
        ortholog_one2many: ['AT1G15520'],
        within_species_paralog: ['SORBI_3004G087800']
      }
    }
  },
  { _id: '353.004G093500', system_name: 'sorghum_353' },
  { _id: 'IS12661.004G086400', system_name: 'sorghum_is12661' },
  { _id: 'AT1G15520', system_name: 'arabidopsis_thaliana' },
  { _id: 'SORBI_3004G087800', system_name: 'sorghum_bicolor' }
];

function stubMongo() {
  const coll = {
    find(query) {
      const out = typeof query._id === 'string'
        ? GENE_DOCS.filter(function (d) { return d._id === query._id; })
        : GENE_DOCS.filter(function (d) { return query._id.$in.indexOf(d._id) >= 0; })
          .map(function (d) { return { _id: d._id, system_name: d.system_name }; });
      return { toArray: function () { return Promise.resolve(JSON.parse(JSON.stringify(out))); } };
    }
  };
  return { genes: { mongoCollection: async function () { return coll; } } };
}

// Fresh check/normalize deps; extra entries (e.g. store, check) are merged over them.
function checkDeps(extra) {
  return Object.assign({ cfg: CFG, catalog: CATALOG, assemblies: stubAssemblies(), mongo: stubMongo(), log: quiet }, extra || {});
}

// ---- fixture lists ------------------------------------------------------------------------------------

function listFiles(dir, re) {
  return fs.readdirSync(dir, { withFileTypes: true })
    .filter(function (ent) { return ent.isFile() && re.test(ent.name); })
    .map(function (ent) { return ent.name; })
    .sort();
}

// The check request fixtures that predate the feature; its check-genotyping-*.json fixtures (§6.3) are not pinned.
function checkFixtureNames() {
  return listFiles(REQUESTS_DIR, /^check-(?!genotyping-).*\.json$/);
}

function designFixtureNames() {
  return listFiles(REQUESTS_DIR, /^design-.*\.json$/);
}

// Top level only: the genotyping recordings live in fixtures/primer3/genotyping/ (§6.3).
function primer3OutputNames() {
  return listFiles(PRIMER3_DIR, /\.output\.txt$/);
}

function readRequest(name) {
  return JSON.parse(fs.readFileSync(path.join(REQUESTS_DIR, name), 'utf8'));
}

function readPrimer3Output(name) {
  return fs.readFileSync(path.join(PRIMER3_DIR, name), 'utf8');
}

// sha256 of the key-order-independent body, so a regenerated fixture with the same content keeps its digest.
function requestDigest(body) {
  return crypto.createHash('sha256').update(jobs.canonicalJSON(body)).digest('hex');
}

// ---- snapshots ----------------------------------------------------------------------------------------

async function jobIdSnapshot() {
  const ids = {};
  const digests = {};
  for (const name of checkFixtureNames()) {
    const body = readRequest(name);
    digests[name] = requestDigest(body);
    const norm = await normalize(body, checkDeps());
    ids[name] = jobs.jobId(norm.request, norm.dbs, check.ALGORITHM_VERSION);
  }
  return {
    about: 'Genotyping spec §7.2 C1: jobs.jobId(norm.request, norm.dbs, check.ALGORITHM_VERSION) of each check-*.json request ' +
      'fixture that predates the feature, normalized offline with the stubs of test/primers/tools/record_compat.js. ' +
      'request_sha256 is sha256(canonicalJSON(body)) of the body each id was recorded from. Re-record only deliberately.',
    algorithm_version: check.ALGORITHM_VERSION,
    ids: ids,
    request_sha256: digests
  };
}

function designNormalizeSnapshot() {
  const normalized = {};
  for (const name of designFixtureNames()) normalized[name] = design.normalize(readRequest(name), CFG);
  return {
    about: 'Genotyping spec §7.2 C2: design.normalize(body, built-in config) of each design-*.json request fixture that ' +
      'predates the feature. Recorded by test/primers/tools/record_compat.js; re-record only deliberately.',
    normalized: normalized
  };
}

function extractPairsSnapshot() {
  const pairs = {};
  for (const name of primer3OutputNames()) {
    const rec = boulder.parseRecord(readPrimer3Output(name));
    if (!rec.complete) throw new Error(name + ' is not a complete Boulder-IO record');
    pairs[name] = boulder.extractPairs(rec.tags);
  }
  return {
    about: 'Genotyping spec §7.2 C3: boulder.extractPairs of each recorded Primer3 output at the top level of ' +
      'test/primers/fixtures/primer3. Recorded by test/primers/tools/record_compat.js; re-record only deliberately.',
    pairs: pairs
  };
}

async function main(argv) {
  const write = argv.indexOf('--write') >= 0;
  const snapshots = [
    [SNAPSHOTS.jobIds, await jobIdSnapshot()],
    [SNAPSHOTS.designNormalize, designNormalizeSnapshot()],
    [SNAPSHOTS.extractPairs, extractPairsSnapshot()]
  ];
  let differ = 0;
  for (const [file, value] of snapshots) {
    const text = JSON.stringify(value, null, 2) + '\n';
    const old = fs.existsSync(file) ? fs.readFileSync(file, 'utf8') : null;
    const rel = path.relative(ROOT, file);
    if (old === text) {
      console.log('unchanged  ' + rel);
      continue;
    }
    differ++;
    if (write) {
      fs.writeFileSync(file, text);
      console.log((old === null ? 'created    ' : 'rewritten  ') + rel);
    } else {
      console.log((old === null ? 'missing    ' : 'differs    ') + rel);
    }
  }
  if (differ > 0 && !write) {
    console.log(differ + ' snapshot(s) differ from the current code; rerun with --write only if that change is intended');
    process.exitCode = 1;
  }
}

if (require.main === module) {
  main(process.argv.slice(2)).catch(function (e) {
    console.error(e && e.stack ? e.stack : e);
    process.exit(2);
  });
}

module.exports = {
  ROOT,
  REQUESTS_DIR,
  PRIMER3_DIR,
  SNAPSHOTS,
  CFG,
  checkDeps,
  checkFixtureNames,
  designFixtureNames,
  primer3OutputNames,
  readRequest,
  readPrimer3Output,
  requestDigest,
  jobIdSnapshot,
  designNormalizeSnapshot,
  extractPairsSnapshot
};
