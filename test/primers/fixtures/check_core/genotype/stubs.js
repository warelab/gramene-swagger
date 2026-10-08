'use strict';

// Offline stand-ins for the genotyping check tests (genotyping spec §7.3): the sequence.js contract over one real window of
// sorghum_bicolor chromosome 1, a resolver whose genomes carry their catalog sizes, and the example sets of §2.9-§2.11 and
// §5.7 as check pairs.
// Provenance (recorded 2026-09-15, read-only): sorghum_bicolor_1_9000-15500.plus.txt is sequence.fetch of
//   /scratch/olson/fasta/sorghum_bicolor/dna/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel.fa.gz, region 1, 9000-15500,
//   plus strand, uppercase, no newline inside; where they overlap it equals test/primers/fixtures/design/
//   sorghum_bicolor_1_{10500-12100,11080-15099,13400-14000}.plus.txt. Chromosome 1 is 80,884,392 bp (.fai, catalog fixture).

const fs = require('fs');
const path = require('path');
const crypto = require('crypto');

const ROOT = path.resolve(__dirname, '..', '..', '..', '..', '..');
const config = require(path.join(ROOT, 'api/helpers/primers/config'));
const genomes = require(path.join(ROOT, 'api/helpers/primers/genomes'));
const mapsFixture = require(path.join(ROOT, 'test/primers/fixtures/catalog/maps_sorghum_v11.json'));
const taxonomyFixture = require(path.join(ROOT, 'test/primers/fixtures/catalog/taxonomy_sorghum_v11.json'));

const CATALOG = genomes.buildCatalog(mapsFixture, taxonomyFixture);
const quiet = { info() {}, warn() {}, error() {}, log() {} };
const WINDOW = Object.freeze({ system_name: 'sorghum_bicolor', region: '1', start: 9000, end: 15500 });
const WINDOW_SEQ = fs.readFileSync(path.join(__dirname, 'sorghum_bicolor_1_9000-15500.plus.txt'), 'utf8').trim();
if (WINDOW_SEQ.length !== WINDOW.end - WINDOW.start + 1) throw new Error('genotype stubs: the window fixture has the wrong length');
// sorghum_bicolor's .fai total; the other genomes are sized by their catalog map regions.
const BICOLOR_TOTAL_BASES = 708735318;

function makeCfg(overrides) {
  return config._build({ env: {}, fileConfig: {}, overrides: overrides }).config;
}

function fp(name) {
  return crypto.createHash('sha1').update('fp:' + name).digest('hex');
}

function regionLengthOf(systemName, region) {
  const g = CATALOG.bySystemName.get(systemName);
  const i = g ? g.regions.names.indexOf(region) : -1;
  return i >= 0 ? g.regions.lengths[i] : undefined;
}

// Plus-strand bases of sorghum_bicolor 1:start-end from the fixture window.
function bases(start, end) {
  if (start < WINDOW.start || end > WINDOW.end || end < start) throw new Error('genotype stubs: no fixture bases for 1:' + start + '-' + end);
  return WINDOW_SEQ.slice(start - WINDOW.start, end - WINDOW.start + 1);
}

// The sequence.js contract for the reference (sorghum_bicolor, whatever the path): regionLength from the catalog, fetch from
// the fixture. calls records every fetch as [region, start, end]. opts.fail: an Error that fetch rejects with.
function sequenceStub(opts) {
  opts = opts || {};
  const calls = [];
  return {
    calls: calls,
    async regionLength(fastaPath, region) {
      return regionLengthOf('sorghum_bicolor', region);
    },
    async fetch(fastaPath, region, start, end, strand) {
      calls.push([region, start, end]);
      if (opts.fail) throw opts.fail;
      if (region !== WINDOW.region || strand !== 1) throw new Error('genotype stubs: only plus-strand region 1 is recorded');
      return bases(start, end);
    }
  };
}

// Every catalog genome resolves with sequence and BLAST DBs; sizes as above; overrides[name] is merged in.
function stubAssemblies(overrides) {
  overrides = overrides || {};
  function build(name) {
    const g = CATALOG.bySystemName.get(name);
    if (!g) throw genomes.unknownGenomeError(name);
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
      total_bases: name === 'sorghum_bicolor' ? BICOLOR_TOTAL_BASES : g.regions.lengths.reduce(function (a, b) { return a + b; }, 0),
      num_sequences: g.regions.names.length,
      fingerprint: fp(name),
      warnings: []
    }, overrides[name] || {}));
  }
  return {
    async resolve(name) {
      return build(name);
    },
    async resolveMany(names) {
      return names.map(function (n) { return { system_name: n, resolved: build(n), error: null }; });
    }
  };
}

// check/normalize deps; extra entries are merged over them.
function checkDeps(extra) {
  return Object.assign({ cfg: makeCfg(), catalog: CATALOG, assemblies: stubAssemblies(), sequence: sequenceStub(), mongo: {}, log: quiet }, extra || {});
}

function clone(x) {
  return JSON.parse(JSON.stringify(x));
}

function pair(id, left, right, region, start, end) {
  return { id: id, left: left, right: right, expected: { region: region, start: start, end: end } };
}

// §2.11: design.check.request of §2.9 (rs871475760, KASP S1 reverse and S2 forward) narrowed to three genomes.
const BODY_2_11 = Object.freeze({
  system_name: 'sorghum_bicolor',
  mode: 'region',
  checks: ['specificity', 'pangenome'],
  genomes: ['sorghum_bicolorv5', 'sorghum_pi180348', 'sorghum_pi329250'],
  pairs: [
    pair('S1_REF', 'AGCTTCTCTAAGTGGTTATCCGA', 'ATCTTTGACTAGCGAGAAATTCAG', '1', 11068, 11132),
    pair('S1_ALT', 'AGCTTCTCTAAGTGGTTATCCGA', 'ATCTTTGACTAGCGAGAAATTCAT', '1', 11068, 11132),
    pair('S2_REF', 'GGTTATCCGAATATAGTCATACTCTATTC', 'TCTTTGTCTACTGAGAAATCCAGA', '1', 11081, 11172),
    pair('S2_ALT', 'GGTTATCCGAATATAGTCATACTCTATTA', 'TCTTTGTCTACTGAGAAATCCAGA', '1', 11081, 11172)
  ],
  genotyping: {
    variant: { region: '1', position: 11109, ref: 'C', alt: 'A' },
    sets: [
      { id: 'S1', ref_pair: 'S1_REF', alt_pair: 'S1_ALT' },
      { id: 'S2', ref_pair: 'S2_REF', alt_pair: 'S2_ALT' }
    ]
  }
});

// Region-mode bodies over the example sets of the other variants, one set per [id, orientation, as_ref, as_alt, common,
// expected start, expected end].
function body(variant, sets, extra) {
  const pairs = [];
  sets.forEach(function (s) {
    const forward = s[1] === 'forward';
    pairs.push(pair(s[0] + '_REF', forward ? s[2] : s[4], forward ? s[4] : s[2], variant.region, s[5], s[6]));
    pairs.push(pair(s[0] + '_ALT', forward ? s[3] : s[4], forward ? s[4] : s[3], variant.region, s[5], s[6]));
  });
  return Object.assign({
    system_name: 'sorghum_bicolor',
    mode: 'region',
    pairs: pairs,
    genotyping: { variant: variant, sets: sets.map(function (s) { return { id: s[0], ref_pair: s[0] + '_REF', alt_pair: s[0] + '_ALT' }; }) }
  }, extra || {});
}

const SETS = Object.freeze({
  // §2.10(b) rs5413864115, 1:11282:CA:C: S1 forward (88/87 bp), S2 reverse (129/128 bp).
  rs5413864115: [
    ['S1', 'forward', 'ACAGATGATTTTCCAAATGATGATTCAAA', 'ACAGATGATTTTCCAAATGATGATTCAAG', 'CCCCATGTTTTTGTTCCTTCCA', 11257, 11344],
    ['S2', 'reverse', 'AGAGTCTTTTCAAATTTCACACTTT', 'AGAGTCTTTTCAAATTTCACACTTG', 'ACAAAAATAGCTCTCTAGAGTATACACA', 11179, 11307]
  ],
  // §2.10(c) tmp_1_11502_C_CGT, 1:11502:C:CGT: S1 reverse (72/74 bp).
  tmp_1_11502_C_CGT: [
    ['S1', 'reverse', 'GCAGGAAAAGAAATCCTAACATCATATG', 'GCAGGAAAAGAAATCCTAACATCATATA', 'AGGATCTTTGCAACCCTGTGTT', 11458, 11529]
  ],
  // §2.10(d) rs871475760 gel AS-PCR, both allele-specific primers with a -2 mismatch (253 bp).
  rs871475760_as_pcr: [
    ['S1', 'reverse', 'ATCTTTGACTAGCGAGAAATTCGG', 'ATCTTTGACTAGCGAGAAATTCGT', 'TGCATCAACAAATGTGCTATGTGT', 10880, 11132]
  ],
  // §5.7 rs5413863413, 1:13735:TTGG:T: pair 0 of the level-0 design runs with ALT primers derived by §4.9 (forward 103 bp,
  // reverse 82 bp).
  rs5413863413: [
    ['F', 'forward', 'CCATTGGTGGTGGTGGTG', 'CCATTGGTGGTGGTGGTA', 'CCCCGGTTATAACTGAGGATGG', 13732, 13834],
    ['R', 'reverse', 'GTACTTCTACCACCACCACCAC', 'GTACTTCTACCACCACCACCAA', 'TCCAAGACAAGGTTCATGGAGT', 13678, 13759]
  ]
});

const VARIANTS = Object.freeze({
  rs871475760: { region: '1', position: 11109, ref: 'C', alt: 'A' },
  rs5413864115: { region: '1', position: 11282, ref: 'CA', alt: 'C' },
  tmp_1_11502_C_CGT: { region: '1', position: 11502, ref: 'C', alt: 'CGT' },
  rs5413863413: { region: '1', position: 13735, ref: 'TTGG', alt: 'T' }
});

module.exports = {
  ROOT,
  CATALOG,
  WINDOW,
  BICOLOR_TOTAL_BASES,
  BODY_2_11,
  SETS,
  VARIANTS,
  quiet,
  makeCfg,
  fp,
  bases,
  regionLengthOf,
  sequenceStub,
  stubAssemblies,
  checkDeps,
  clone,
  body
};
