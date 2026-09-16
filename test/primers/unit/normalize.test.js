'use strict';

// check/normalize.js: spec §A.2.3 handler rules and plan Overrides.
// Offline tests use the real sorghum_v11 catalog fixture, a stub resolver and a stub mongo.
// PRIMERS_REALDATA=1 adds read-only tests against mongo (sorghum11) and /scratch/olson/fasta (no BLAST).

const test = require('node:test');
const { describe, after } = test;
const should = require('should');
const crypto = require('crypto');

const { normalize } = require('../../../api/helpers/primers/check/normalize');
const cost = require('../../../api/helpers/primers/check/cost');
const config = require('../../../api/helpers/primers/config');
const genomes = require('../../../api/helpers/primers/genomes');
const jobs = require('../../../api/helpers/primers/jobs');
const { PrimerHttpError } = require('../../../api/helpers/primers/errors');
const mapsFixture = require('../fixtures/catalog/maps_sorghum_v11.json');
const taxonomyFixture = require('../fixtures/catalog/taxonomy_sorghum_v11.json');

const REALDATA = process.env.PRIMERS_REALDATA === '1';
const CATALOG = genomes.buildCatalog(mapsFixture, taxonomyFixture);
const quiet = { info() {}, warn() {}, error() {}, log() {} };

function makeCfg(overrides) {
  return config._build({ env: {}, fileConfig: {}, overrides: overrides }).config;
}
const CFG = makeCfg();

function fp(name) {
  return crypto.createHash('sha1').update('fp:' + name).digest('hex');
}

// Stub resolver: every catalog genome has sequence, dna and cdna BLAST DBs and 0.7 Gb, unless overridden
// (an override object is merged, an Error is thrown).
function stubAssemblies(overrides) {
  overrides = overrides || {};
  const calls = { resolve: [], resolveMany: [] };
  function build(name) {
    const g = CATALOG.bySystemName.get(name);
    if (!g) throw genomes.unknownGenomeError(name);
    const o = overrides[name];
    if (o instanceof Error) throw o;
    const dir = '/scratch/olson/fasta/' + name;
    return Object.freeze(Object.assign({
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
    }, o || {}));
  }
  return {
    calls: calls,
    async resolve(name) {
      calls.resolve.push(name);
      return build(name);
    },
    async resolveMany(names) {
      calls.resolveMany.push(names.slice());
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

const P2_GENE = {
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
      ortholog_many2many: ['Zm00001eb205480', '353.004G093500'],
      within_species_paralog: ['SORBI_3004G087800']
    }
  }
};
const OTHER_GENES = [
  { _id: '353.004G093500', system_name: 'sorghum_353' },
  { _id: 'IS12661.004G086400', system_name: 'sorghum_is12661' },
  { _id: 'AT1G15520', system_name: 'arabidopsis_thaliana' },
  { _id: 'SORBI_3004G087800', system_name: 'sorghum_bicolor' },
  { _id: 'NOTX_GENE', system_name: 'sorghum_bicolor', gene_structure: { transcripts: [] } }
];

function stubMongo(opts) {
  opts = opts || {};
  const docs = [P2_GENE].concat(OTHER_GENES);
  const calls = [];
  const coll = {
    find(query, options) {
      calls.push({ query: query, options: options });
      let out;
      if (typeof query._id === 'string') {
        out = docs.filter(function (d) { return d._id === query._id; });
      } else {
        out = docs.filter(function (d) { return query._id.$in.indexOf(d._id) >= 0; })
          .map(function (d) { return { _id: d._id, system_name: d.system_name }; });
      }
      return {
        toArray: function () {
          return opts.hang ? new Promise(function () {}) : Promise.resolve(JSON.parse(JSON.stringify(out)));
        }
      };
    }
  };
  return { calls: calls, genes: { mongoCollection: async function () { return opts.down ? undefined : coll; } } };
}

function deps(extra) {
  return Object.assign({ cfg: CFG, catalog: CATALOG, assemblies: stubAssemblies(), mongo: stubMongo(), log: quiet }, extra || {});
}

const P2 = { id: 'P2', left: 'GGACAGCTCCACAACATATCAG', right: 'GGACATTTGAAGCCCATGGCC' };

function body(extra) {
  return Object.assign({ system_name: 'sorghum_bicolor', pairs: [Object.assign({}, P2)] }, extra || {});
}

async function rejects(promise, status, code) {
  let err;
  try {
    await promise;
  } catch (e) {
    err = e;
  }
  should.exist(err, 'expected ' + status + ' ' + code);
  if (err.name !== 'PrimerHttpError') throw err;
  err.status.should.equal(status);
  err.code.should.equal(code);
  return err;
}

// ---- defaults and canonical request ----------------------------------------------------------------

test('a minimal region request is filled with defaults and uppercased primers', async function () {
  const d = deps();
  const norm = await normalize(body({ pairs: [{ id: 'p1', left: 'ggacagctccacaacatatcag', right: 'GGACATTTGAAGCCCATGGCC' }] }), d);
  norm.request.should.eql({
    system_name: 'sorghum_bicolor',
    mode: 'region',
    checks: ['specificity'],
    genomes: [],
    params: {
      max_product_size: 4000, ignore_mismatches: 6, min_total_mismatches: 2, min_3p_mismatches: 2,
      three_prime_window: 5, include_unlikely: false, repeat_site_threshold: 5, max_amplifying_mismatches: 3
    },
    pairs: [{ id: 'p1', left: 'GGACAGCTCCACAACATATCAG', right: 'GGACATTTGAAGCCCATGGCC' }]
  });
  norm.kind.should.equal('specificity');
  norm.warnings.should.eql([]);
  norm.estimate.should.eql({ cpu_s: 9, total: 1 }); // ceil(2 x 0.7 Gb x 5.2 + 2 primers x 1 genome task x 0.6 re-alignment)
  norm.dbs.should.eql({ sorghum_bicolor: fp('sorghum_bicolor') });
  Object.keys(norm.resolved.assemblies).should.eql(['sorghum_bicolor']);
  Object.isFrozen(norm.resolved.assemblies.sorghum_bicolor).should.be.false();
  should(norm.resolved.gene).be.null();
  norm.resolved.species.should.eql({ taxon_id: 4558, name: 'Sorghum bicolor' });
  d.mongo.calls.length.should.equal(0);
  d.assemblies.calls.resolveMany.length.should.equal(0);
});

test('explicit params win over defaults; the config defaults block is used', async function () {
  const norm = await normalize(body({ params: { ignore_mismatches: 4, include_unlikely: true } }),
    deps({ cfg: makeCfg({ check: { defaults: { repeat_site_threshold: 9 } } }) }));
  norm.request.params.ignore_mismatches.should.equal(4);
  norm.request.params.include_unlikely.should.be.true();
  norm.request.params.repeat_site_threshold.should.equal(9);
  norm.request.params.max_product_size.should.equal(4000);
});

// ---- validation --------------------------------------------------------------------------------------

test('duplicate pair ids -> 400 DUPLICATE_PAIR_ID', async function () {
  const err = await rejects(normalize(body({ pairs: [P2, Object.assign({}, P2, { right: 'GTGAACATCATGCTGCCCGATG' })] }), deps()),
    400, 'DUPLICATE_PAIR_ID');
  err.details.should.eql({ pair_id: 'P2' });
});

test('more than 20 distinct primers -> 400 TOO_MANY_PRIMERS; shared and case-variant primers count once', async function () {
  const bases = 'ACGT';
  function primer(n) {
    let s = '';
    for (let i = 0; i < 20; i++) s += bases[(n >> (i % 8)) & 3 ^ (i % 4)];
    return 'AAAA' + s.slice(4) + n.toString(4).replace(/./g, function (c) { return bases[Number(c)]; });
  }
  const pairs = [];
  for (let i = 0; i < 10; i++) pairs.push({ id: 'p' + i, left: primer(2 * i + 1), right: primer(2 * i + 2) });
  new Set(pairs.map(function (p) { return p.left; }).concat(pairs.map(function (p) { return p.right; }))).size.should.equal(20);
  // With the default limits (10 pairs, 20 primers) the cap cannot be exceeded; 20 distinct primers pass.
  (await normalize(body({ pairs: pairs }), deps())).request.pairs.length.should.equal(10);

  const lowCap = makeCfg({ check: { max_unique_primers: 19 } });
  const err = await rejects(normalize(body({ pairs: pairs }), deps({ cfg: lowCap })), 400, 'TOO_MANY_PRIMERS');
  err.details.should.eql({ unique_primers: 20, max: 19 });

  // A primer shared between pairs (in any case) counts once: 11 distinct primers.
  const shared = pairs.map(function (p) { return { id: p.id, left: p.left, right: pairs[0].right.toLowerCase() }; });
  const norm = await normalize(body({ pairs: shared }), deps({ cfg: makeCfg({ check: { max_unique_primers: 11 } }) }));
  norm.request.pairs[3].right.should.equal(pairs[0].right);
  norm.estimate.cpu_s.should.equal(Math.ceil(11 * (0.7 * 5.2 + 0.6)));
});

test('shape errors -> 400 INVALID_REQUEST', async function () {
  const cases = [
    null,
    { pairs: [P2] },
    body({ system_name: '../etc' }),
    body({ bogus: 1 }),
    body({ mode: 'transcriptome' }),
    body({ checks: ['transcriptome'] }),
    body({ checks: [] }),
    body({ pairs: [] }),
    body({ pairs: [Object.assign({}, P2, { left: 'ACGTN' + P2.left })] }),
    body({ pairs: [Object.assign({}, P2, { left: 'ACGTACGTAC' })] }),
    body({ pairs: [Object.assign({}, P2, { extra: true })] }),
    body({ pairs: [Object.assign({}, P2, { expected: { region: '4', start: 0, end: 10 } })] }),
    body({ params: { ignore_mismatches: 7 } }),
    body({ params: { max_mismatches: 3 } }),
    body({ params: { include_unlikely: 'yes' } }),
    body({ params: { max_amplifying_mismatches: 6 } }),
    body({ params: { max_amplifying_mismatches: -1 } }),
    body({ params: { max_amplifying_mismatches: 2.5 } }),
    body({ genomes: ['Sorghum_353'] }),
    body({ gene_id: '' })
  ];
  const eleven = [];
  for (let i = 0; i < 11; i++) eleven.push(Object.assign({}, P2, { id: 'p' + i }));
  cases.push(body({ pairs: eleven }));
  for (const c of cases) await rejects(normalize(c, deps()), 400, 'INVALID_REQUEST');
});

test('max_amplifying_mismatches: default 3, below ignore_mismatches or 400 INVALID_PARAMS; an omitted cap is lowered to fit', async function () {
  (await normalize(body({ params: { max_amplifying_mismatches: 5 } }), deps())).request.params.max_amplifying_mismatches.should.equal(5);
  (await normalize(body({ params: { max_amplifying_mismatches: 0, ignore_mismatches: 3 } }), deps())).request.params
    .should.match({ max_amplifying_mismatches: 0, ignore_mismatches: 3 });
  const e = await rejects(normalize(body({ params: { max_amplifying_mismatches: 4, ignore_mismatches: 4 } }), deps()), 400, 'INVALID_PARAMS');
  e.details.should.eql({ field: 'params.max_amplifying_mismatches', max_amplifying_mismatches: 4, ignore_mismatches: 4 });
  await rejects(normalize(body({ params: { max_amplifying_mismatches: 5, ignore_mismatches: 5 } }), deps()), 400, 'INVALID_PARAMS');
  await rejects(normalize(body({ params: { max_amplifying_mismatches: 4, ignore_mismatches: 3 } }), deps()), 400, 'INVALID_PARAMS');
  (await normalize(body({ params: { ignore_mismatches: 3 } }), deps())).request.params.max_amplifying_mismatches.should.equal(2);
  (await normalize(body({ params: { ignore_mismatches: 4 } }), deps())).request.params.max_amplifying_mismatches.should.equal(3);
  (await normalize(body(), deps({ cfg: makeCfg({ check: { defaults: { max_amplifying_mismatches: 2 } } }) })))
    .request.params.max_amplifying_mismatches.should.equal(2);
});

test('gene and transcript modes require gene_id', async function () {
  await rejects(normalize(body({ mode: 'transcript' }), deps()), 400, 'INVALID_REQUEST');
  await rejects(normalize(body({ mode: 'gene' }), deps()), 400, 'INVALID_REQUEST');
});

// ---- expected / max_product_size -------------------------------------------------------------------

test('transcript and sequence modes strip expected with warning EXPECTED_IGNORED; gene mode keeps it', async function () {
  const expected = { region: '4', start: 7423537, end: 7423746 };
  const withExpected = [Object.assign({}, P2, { expected: expected })];
  const tx = await normalize(body({ mode: 'transcript', gene_id: 'SORBI_3004G087700', pairs: withExpected }), deps());
  should(tx.request.pairs[0].expected).be.undefined();
  tx.warnings.map(function (w) { return w.code; }).should.eql(['EXPECTED_IGNORED']);
  const seq = await normalize(body({ mode: 'sequence', pairs: withExpected }), deps());
  should(seq.request.pairs[0].expected).be.undefined();
  seq.warnings[0].code.should.equal('EXPECTED_IGNORED');
  const gene = await normalize(body({ mode: 'gene', gene_id: 'SORBI_3004G087700', pairs: withExpected }), deps());
  gene.request.pairs[0].expected.should.eql(expected);
  gene.warnings.should.eql([]);
});

test('max_product_size auto-raise: min(10000, max(current, ceil(1.2 x largest expected))) with MAX_PRODUCT_SIZE_RAISED', async function () {
  function pairWith(size, id) {
    return Object.assign({}, P2, { id: id || 'P2', expected: { region: '4', start: 1000, end: 1000 + size - 1 } });
  }
  let n = await normalize(body({ pairs: [pairWith(3000)] }), deps());
  n.request.params.max_product_size.should.equal(4000);
  n.warnings.should.eql([]);

  n = await normalize(body({ pairs: [pairWith(5000), Object.assign(pairWith(2000, 'small'), { left: 'GATCGACAATCCGACGATAGAAG' })] }), deps());
  n.request.params.max_product_size.should.equal(6000);
  n.warnings.map(function (w) { return w.code; }).should.eql(['MAX_PRODUCT_SIZE_RAISED']);

  n = await normalize(body({ pairs: [pairWith(9000)] }), deps());
  n.request.params.max_product_size.should.equal(10000);

  n = await normalize(body({ pairs: [pairWith(10000)] }), deps());
  n.request.params.max_product_size.should.equal(10000);

  // An explicit smaller value is raised only as far as needed.
  n = await normalize(body({ params: { max_product_size: 1000 }, pairs: [pairWith(900)] }), deps());
  n.request.params.max_product_size.should.equal(1080);
  n = await normalize(body({ params: { max_product_size: 8000 }, pairs: [pairWith(3000)] }), deps());
  n.request.params.max_product_size.should.equal(8000);
  n.warnings.should.eql([]);

  const err = await rejects(normalize(body({ pairs: [pairWith(10001)] }), deps()), 400, 'PRODUCT_TOO_LONG_TO_CHECK');
  err.details.should.eql({ pair_id: 'P2', size: 10001, max: 10000 });
  await rejects(normalize(body({ pairs: [Object.assign({}, P2, { expected: { region: '4', start: 20, end: 10 } })] }), deps()),
    400, 'INVALID_REQUEST');
  // Stripped expected (transcript mode) never raises or refuses.
  n = await normalize(body({ mode: 'transcript', gene_id: 'SORBI_3004G087700', pairs: [pairWith(20000)] }), deps());
  n.request.params.max_product_size.should.equal(4000);
});

// ---- gene lookup -------------------------------------------------------------------------------------

test('gene mode: gene doc lookup, canonical transcript default, orthologs grouped with one genes.find', async function () {
  const d = deps();
  const norm = await normalize(body({ mode: 'gene', gene_id: 'SORBI_3004G087700' }), d);
  norm.request.gene_id.should.equal('SORBI_3004G087700');
  norm.request.transcript_id.should.equal('SORBI_3004G087700.3');
  norm.resolved.gene.should.eql({
    id: 'SORBI_3004G087700',
    system_name: 'sorghum_bicolor',
    taxon_id: 4558006,
    location: { region: '4', start: 7421357, end: 7428285, strand: 1, map: 'GCA_000003195.3' },
    transcripts: ['SORBI_3004G087700.3', 'SORBI_3004G087700.1', 'SORBI_3004G087700.2'],
    canonical_transcript: 'SORBI_3004G087700.3',
    transcript_id: 'SORBI_3004G087700.3',
    orthologs: {
      sorghum_353: ['353.004G093500'],
      sorghum_is12661: ['IS12661.004G086400'],
      arabidopsis_thaliana: ['AT1G15520']
    }
  });
  d.mongo.calls.length.should.equal(2);
  d.mongo.calls[0].query.should.eql({ _id: 'SORBI_3004G087700' });
  d.mongo.calls[1].query._id.$in.sort().should.eql(['353.004G093500', 'AT1G15520', 'IS12661.004G086400', 'Zm00001eb205480']);
  d.mongo.calls[1].options.should.eql({ fields: { system_name: 1 } });
});

test('transcript mode: explicit transcript kept, unknown transcript 404, gene without transcripts 404', async function () {
  const norm = await normalize(body({ mode: 'transcript', gene_id: 'SORBI_3004G087700', transcript_id: 'SORBI_3004G087700.1' }), deps());
  norm.request.transcript_id.should.equal('SORBI_3004G087700.1');
  norm.estimate.should.eql({ cpu_s: 11, total: 2 }); // ceil(2 x (0.7 + 0.15) x 5.2 = 8.84, + 2 x 0.6 re-alignment)
  await rejects(normalize(body({ mode: 'transcript', gene_id: 'SORBI_3004G087700', transcript_id: 'SORBI_3004G087700.9' }), deps()),
    404, 'UNKNOWN_TRANSCRIPT');
  await rejects(normalize(body({ mode: 'transcript', gene_id: 'NOTX_GENE' }), deps()), 404, 'UNKNOWN_TRANSCRIPT');
});

test('unknown gene -> 404 UNKNOWN_GENE; gene of another genome -> 400 SYSTEM_NAME_MISMATCH', async function () {
  const e1 = await rejects(normalize(body({ mode: 'gene', gene_id: 'NOPE' }), deps()), 404, 'UNKNOWN_GENE');
  e1.details.should.eql({ gene_id: 'NOPE' });
  const e2 = await rejects(normalize(body({ system_name: 'sorghum_353', mode: 'gene', gene_id: 'SORBI_3004G087700' }), deps()),
    400, 'SYSTEM_NAME_MISMATCH');
  e2.details.should.eql({ gene_id: 'SORBI_3004G087700', system_name: 'sorghum_353', gene_system_name: 'sorghum_bicolor' });
});

test('region and sequence modes ignore gene_id / transcript_id and never touch mongo', async function () {
  const d = deps();
  const norm = await normalize(body({ mode: 'region', gene_id: 'SORBI_3004G087700', transcript_id: 'x' }), d);
  norm.request.should.not.have.property('gene_id');
  norm.request.should.not.have.property('transcript_id');
  d.mongo.calls.length.should.equal(0);
});

test('mongo unavailable (collection undefined or a hung query) -> 503 MONGO_UNAVAILABLE', async function () {
  const e = await rejects(normalize(body({ mode: 'gene', gene_id: 'SORBI_3004G087700' }), deps({ mongo: stubMongo({ down: true }) })),
    503, 'MONGO_UNAVAILABLE');
  e.details.retry_after_s.should.equal(30);
  await rejects(normalize(body({ mode: 'gene', gene_id: 'SORBI_3004G087700' }),
    deps({ mongo: stubMongo({ hang: true }), mongo_timeout_ms: 50 })), 503, 'MONGO_UNAVAILABLE');
});

// ---- reference assembly ------------------------------------------------------------------------------

test('reference: unknown genome 404; no dna BLAST DB 422; transcript mode also needs cdna', async function () {
  await rejects(normalize(body({ system_name: 'sorghum_nope' }), deps()), 404, 'UNKNOWN_GENOME');
  const noDna = stubAssemblies({ sorghum_bicolor: { blastdb: { dna: null, cdna: '/x/y.cdna.all' } } });
  const e1 = await rejects(normalize(body(), deps({ assemblies: noDna })), 422, 'NO_BLASTDB');
  e1.details.should.eql({ system_name: 'sorghum_bicolor', db: 'dna' });
  const noCdna = stubAssemblies({ sorghum_bicolor: { blastdb: { dna: '/x/y.dna.toplevel', cdna: null } } });
  const e2 = await rejects(normalize(body({ mode: 'transcript', gene_id: 'SORBI_3004G087700' }), deps({ assemblies: noCdna })),
    422, 'NO_BLASTDB');
  e2.details.should.eql({ system_name: 'sorghum_bicolor', db: 'cdna' });
  (await normalize(body({ mode: 'gene', gene_id: 'SORBI_3004G087700' }), deps({ assemblies: noCdna }))).kind.should.equal('specificity');
  const ambiguous = stubAssemblies({ sorghum_bicolor: new PrimerHttpError(422, 'AMBIGUOUS_ASSEMBLY', 'x', { candidates: ['a', 'b'] }) });
  await rejects(normalize(body(), deps({ assemblies: ambiguous })), 422, 'AMBIGUOUS_ASSEMBLY');
});

test('reference assembly warnings are passed on, prefixed with the system name', async function () {
  const a = stubAssemblies({ sorghum_bicolor: { warnings: [{ code: 'ASSEMBLY_MISMATCH', message: 'only 1 of 2 map regions match' }] } });
  const norm = await normalize(body(), deps({ assemblies: a }));
  norm.warnings.should.eql([{ code: 'ASSEMBLY_MISMATCH', message: 'sorghum_bicolor: only 1 of 2 map regions match' }]);
});

// ---- pan-genome genomes ------------------------------------------------------------------------------

test('pangenome, genomes omitted: every same-species genome with a dna DB, minus the reference, sorted', async function () {
  const a = stubAssemblies({
    sorghum_rio: { blastdb: { dna: null, cdna: '/x/rio.cdna.all' } },
    sorghum_leoti: { blastdb: { dna: '/x/leoti.dna.toplevel', cdna: null } },
    sorghum_353: new PrimerHttpError(422, 'AMBIGUOUS_ASSEMBLY', 'two assemblies', {})
  });
  const norm = await normalize(body({ checks: ['pangenome'] }), deps({ assemblies: a }));
  norm.kind.should.equal('pangenome');
  norm.request.checks.should.eql(['pangenome', 'specificity']);
  a.calls.resolveMany.length.should.equal(1);
  a.calls.resolveMany[0].length.should.equal(119);
  a.calls.resolveMany[0].should.not.containEql('sorghum_bicolor');
  a.calls.resolveMany[0].should.not.containEql('arabidopsis_thaliana');
  norm.request.genomes.length.should.equal(117);
  norm.request.genomes.should.not.containEql('sorghum_rio');
  norm.request.genomes.should.not.containEql('sorghum_353');
  norm.request.genomes.should.containEql('sorghum_leoti');
  norm.request.genomes.slice().sort().should.eql(norm.request.genomes);
  Object.keys(norm.resolved.assemblies).length.should.equal(118);
  Object.keys(norm.dbs).length.should.equal(118);
  norm.estimate.total.should.equal(118);
  // pan-genome BLAST term × the default pangenome_cpu_factor 2
  norm.estimate.cpu_s.should.equal(Math.ceil(2 * (0.7 * 5.2 + 117 * 0.7 * 2.2 * 2 + 118 * 0.6)));

  // transcript mode needs the cdna DB instead: rio is in, leoti is out.
  const tx = await normalize(body({ mode: 'transcript', gene_id: 'SORBI_3004G087700', checks: ['specificity', 'pangenome'] }),
    deps({ assemblies: a }));
  tx.request.genomes.should.containEql('sorghum_rio');
  tx.request.genomes.should.not.containEql('sorghum_leoti');
  tx.request.genomes.length.should.equal(117);
  tx.estimate.total.should.equal(119);
  tx.estimate.cpu_s.should.equal(Math.ceil(2 * (0.7 * 5.2 + 0.15 * 5.2 + 117 * 0.15 * 2.2 * 2 + 0.6)));
});

test('pangenome over all 119 other genomes (stub sizes) and the A.2.4 progress total 121 in transcript mode', async function () {
  const norm = await normalize(body({ mode: 'transcript', gene_id: 'SORBI_3004G087700', checks: ['pangenome'] }), deps());
  norm.request.genomes.length.should.equal(119);
  norm.estimate.total.should.equal(121);
});

test('pangenome, explicit genomes: deduped, reference dropped, sorted; unknown 404; not checkable 400 with details.genomes', async function () {
  const a = stubAssemblies();
  const norm = await normalize(body({ checks: ['pangenome'], genomes: ['sorghum_rio', 'sorghum_353', 'sorghum_bicolor', 'sorghum_rio'] }),
    deps({ assemblies: a }));
  norm.request.genomes.should.eql(['sorghum_353', 'sorghum_rio']);
  a.calls.resolveMany.should.eql([['sorghum_rio', 'sorghum_353']]);
  norm.estimate.total.should.equal(3);

  const e404 = await rejects(normalize(body({ checks: ['pangenome'], genomes: ['sorghum_353', 'sorghum_nope'] }), deps()),
    404, 'UNKNOWN_GENOME');
  e404.details.should.eql({ system_name: 'sorghum_nope' });

  const other = await rejects(normalize(body({ checks: ['pangenome'], genomes: ['arabidopsis_thaliana', 'sorghum_353'] }), deps()),
    400, 'GENOME_NOT_CHECKABLE');
  other.details.genomes.should.eql(['arabidopsis_thaliana']);
  other.details.reasons.should.eql({ arabidopsis_thaliana: 'other_species' });

  const broken = stubAssemblies({
    sorghum_353: { blastdb: { dna: null, cdna: '/x/353.cdna.all' } },
    sorghum_rio: new PrimerHttpError(422, 'AMBIGUOUS_ASSEMBLY', 'two', {}),
    sorghum_leoti: new Error('EACCES: permission denied')
  });
  const bad = await rejects(normalize(body({ checks: ['pangenome'], genomes: ['sorghum_353', 'sorghum_rio', 'sorghum_leoti', 'sorghum_grassl'] }),
    deps({ assemblies: broken })), 400, 'GENOME_NOT_CHECKABLE');
  bad.details.genomes.should.eql(['sorghum_353', 'sorghum_rio', 'sorghum_leoti']);
  bad.details.reasons.should.eql({ sorghum_353: 'no_blastdb', sorghum_rio: 'ambiguous_assembly', sorghum_leoti: 'assembly_unavailable' });

  const noCdna = stubAssemblies({ sorghum_353: { blastdb: { dna: '/x/353.dna.toplevel', cdna: null } } });
  const tx = await rejects(normalize(body({ mode: 'transcript', gene_id: 'SORBI_3004G087700', checks: ['pangenome'], genomes: ['sorghum_353'] }),
    deps({ assemblies: noCdna })), 400, 'GENOME_NOT_CHECKABLE');
  tx.details.reasons.should.eql({ sorghum_353: 'no_cdna_blastdb' });

  const onlyRef = await rejects(normalize(body({ checks: ['pangenome'], genomes: ['sorghum_bicolor'] }), deps()), 400, 'GENOME_NOT_CHECKABLE');
  onlyRef.details.genomes.should.eql([]);
});

test('genomes without the pangenome check are ignored with warning GENOMES_IGNORED', async function () {
  const d = deps();
  const norm = await normalize(body({ genomes: ['sorghum_353'] }), d);
  norm.request.genomes.should.eql([]);
  norm.kind.should.equal('specificity');
  norm.warnings.map(function (w) { return w.code; }).should.eql(['GENOMES_IGNORED']);
  d.assemblies.calls.resolveMany.length.should.equal(0);
});

test('a species with no other checkable genome -> 400 GENOME_NOT_CHECKABLE', async function () {
  const err = await rejects(normalize(body({ system_name: 'arabidopsis_thaliana', checks: ['pangenome'] }), deps()),
    400, 'GENOME_NOT_CHECKABLE');
  err.details.genomes.should.eql([]);
});

test('cost guard: estimate over max_job_cpu_s -> 422 JOB_TOO_LARGE', async function () {
  const err = await rejects(normalize(body({ checks: ['pangenome'] }), deps({ cfg: makeCfg({ check: { max_job_cpu_s: 100 } }) })),
    422, 'JOB_TOO_LARGE');
  err.details.should.eql({ estimate_cpu_s: Math.ceil(2 * (0.7 * 5.2 + 119 * 0.7 * 2.2 * 2 + 120 * 0.6)), limit: 100 });
});

// ---- real data ---------------------------------------------------------------------------------------

describe('real data: normalize against mongo sorghum11 and /scratch/olson/fasta', { skip: REALDATA ? false : 'set PRIMERS_REALDATA=1' }, function () {
  let realCfg;
  function real() {
    if (!realCfg) realCfg = config._build({ env: {} }).config;
    return { cfg: realCfg, log: quiet };
  }
  const P2_EXPECTED = Object.assign({}, P2, { expected: { region: '4', start: 7423537, end: 7423746 } });

  after(function () {
    require('gramene-mongodb-config').closeMongoDatabase();
  });

  test('P2 gene mode, specificity: gene location, canonical transcript, orthologs by genome, estimate', async function () {
    const norm = await normalize({ system_name: 'sorghum_bicolor', mode: 'gene', gene_id: 'SORBI_3004G087700', pairs: [P2_EXPECTED] }, real());
    norm.kind.should.equal('specificity');
    norm.request.transcript_id.should.equal('SORBI_3004G087700.3');
    norm.resolved.gene.location.should.eql({ region: '4', start: 7421357, end: 7428285, strand: 1, map: 'GCA_000003195.3' });
    norm.resolved.gene.transcripts.slice().sort().should.eql(['SORBI_3004G087700.1', 'SORBI_3004G087700.2', 'SORBI_3004G087700.3']);
    norm.resolved.gene.orthologs.sorghum_353.should.containEql('353.004G093500');
    norm.resolved.gene.orthologs.should.not.have.property('sorghum_bicolor'); // paralogs are not orthologs
    const asm = norm.resolved.assemblies.sorghum_bicolor;
    asm.total_bases.should.equal(708735318);
    asm.blastdb.dna.should.match(/\.dna\.toplevel$/);
    asm.blastdb.cdna.should.match(/\.cdna\.all$/);
    norm.dbs.sorghum_bicolor.should.match(/^[0-9a-f]{40}$/);
    norm.estimate.should.eql({ cpu_s: Math.ceil(2 * (0.708735318 * 5.2 + 0.6)), total: 1 });
    norm.warnings.map(function (w) { return w.code; }).should.not.containEql('ASSEMBLY_MISMATCH');
    console.log('# P2 orthologs in', Object.keys(norm.resolved.gene.orthologs).length, 'genomes; estimate', norm.estimate);
  });

  test('qPCR transcript mode + pangenome (genomes omitted): 119 genomes with a cDNA DB, progress total 121', async function () {
    const norm = await normalize({
      system_name: 'sorghum_bicolor', mode: 'transcript', gene_id: 'SORBI_3004G087700', checks: ['pangenome', 'specificity'],
      pairs: [{ id: 'Q1', left: 'CCAACAAAGTCATGGATGCACT', right: 'GTGAACATCATGCTGCCCGATG', expected: { region: '4', start: 1, end: 278 } }]
    }, real());
    norm.warnings.map(function (w) { return w.code; }).should.containEql('EXPECTED_IGNORED');
    norm.request.genomes.length.should.equal(119);
    norm.request.genomes.should.not.containEql('sorghum_bicolor');
    norm.estimate.total.should.equal(121);
    Object.keys(norm.dbs).length.should.equal(120);
    norm.estimate.cpu_s.should.equal(Math.ceil(2 * (0.708735318 * 5.2 + 0.15 * 5.2 + 119 * 0.15 * 2.2 * 2 + 0.6)));
    console.log('# transcript pangenome estimate', norm.estimate);
  });

  test('B.13 cost on real assembly sizes: specificity 10 primers ~37, pan-genome over 119 genomes ~3,700 CPU-s with the factor; 20 primers exceed the cap', async function () {
    const norm = await normalize({ system_name: 'sorghum_bicolor', mode: 'gene', gene_id: 'SORBI_3004G087700', checks: ['pangenome'],
      pairs: [P2_EXPECTED] }, real());
    norm.request.genomes.length.should.equal(119);
    const asms = norm.resolved.assemblies;
    const pan = norm.request.genomes.map(function (s) { return asms[s]; });
    pan.forEach(function (a) { a.total_bases.should.be.above(4e8); });
    const est = cost.estimate({ unique_primers: 10, mode: 'gene', reference: asms.sorghum_bicolor, pangenome: pan, cfg: real().cfg });
    est.breakdown.reference.should.equal(36.9);
    est.breakdown.pangenome.should.be.within(3200, 4000); // ~1,840 single-thread CPU-s × pangenome_cpu_factor 2
    const est20 = cost.estimate({ unique_primers: 20, mode: 'gene', reference: asms.sorghum_bicolor, pangenome: pan, cfg: real().cfg });
    est20.cpu_s.should.be.above(6000);
    console.log('# real pan-genome Gb', (pan.reduce(function (s, a) { return s + a.total_bases; }, 0) / 1e9).toFixed(2),
      'estimate 10 primers', est.cpu_s, '20 primers', est20.cpu_s);
  });

  test('real errors: UNKNOWN_GENE, SYSTEM_NAME_MISMATCH, GENOME_NOT_CHECKABLE (other species)', async function () {
    await rejects(normalize({ system_name: 'sorghum_bicolor', mode: 'gene', gene_id: 'SORBI_NOPE_GENE', pairs: [P2] }, real()),
      404, 'UNKNOWN_GENE');
    await rejects(normalize({ system_name: 'sorghum_353', mode: 'gene', gene_id: 'SORBI_3004G087700', pairs: [P2] }, real()),
      400, 'SYSTEM_NAME_MISMATCH');
    const err = await rejects(normalize({ system_name: 'sorghum_bicolor', checks: ['pangenome'], genomes: ['arabidopsis_thaliana', 'sorghum_353'],
      pairs: [P2] }, real()), 400, 'GENOME_NOT_CHECKABLE');
    err.details.genomes.should.eql(['arabidopsis_thaliana']);
  });

  test('real job id is stable under key order and default filling', async function () {
    const a = await normalize({ system_name: 'sorghum_bicolor', mode: 'gene', gene_id: 'SORBI_3004G087700', pairs: [P2_EXPECTED] }, real());
    const b = await normalize({
      pairs: [{ right: P2.right.toLowerCase(), expected: { end: 7423746, start: 7423537, region: '4' }, left: P2.left, id: 'P2' }],
      params: { max_product_size: 4000, include_unlikely: false },
      checks: ['specificity'],
      transcript_id: 'SORBI_3004G087700.3',
      gene_id: 'SORBI_3004G087700',
      mode: 'gene',
      system_name: 'sorghum_bicolor'
    }, real());
    jobs.jobId(a.request, a.dbs, '1').should.equal(jobs.jobId(b.request, b.dbs, '1'));
  });
});
