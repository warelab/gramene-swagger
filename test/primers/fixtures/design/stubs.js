'use strict';

// Offline stand-ins for the design tests, an in-memory sequence source with the sequence.js contract,
// a resolver stub and a config builder that ignores config/default.yaml.
// Provenance (captured 2026-09-12, read-only):
//   genes.json: GET localhost:50011/sorghum_v11/genes?idList=SORBI_3001G000200,SORBI_3001G000700,SORBI_3004G087700
//               &fl=_id,name,system_name,taxon_id,location,gene_structure (transcripts trimmed to id, length,
//               exons, exon_junctions, cds)
//   sorghum_bicolor_1_<start>-<end>.plus.txt: samtools faidx of
//               /scratch/olson/fasta/sorghum_bicolor/dna/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel.fa.gz
//               (plus strand, no newlines); chromosome 1 length 80884392 from its .fai

const fs = require('fs');
const path = require('path');
const ROOT = path.resolve(__dirname, '..', '..', '..', '..');
const { revcomp } = require(path.join(ROOT, 'api/helpers/primers/sequence'));
const { PrimerHttpError } = require(path.join(ROOT, 'api/helpers/primers/errors'));
const config = require(path.join(ROOT, 'api/helpers/primers/config'));
const boulder = require(path.join(ROOT, 'api/helpers/primers/boulder'));

const DIR = __dirname;
const PRIMER3_DIR = path.join(DIR, '..', 'primer3');
const SB_DNA = '/fake/sorghum_bicolor/dna/Sb.dna.toplevel.fa.gz';
const SB_CHR1_LENGTH = 80884392;

function readWindow(file) {
  return fs.readFileSync(path.join(DIR, file), 'utf8').trim();
}

// Plus-strand genomic windows of sorghum_bicolor chromosome 1 (samtools faidx of dna.toplevel).
const WINDOWS = Object.freeze({
  SORBI_3001G000200: Object.freeze({ region: '1', start: 11080, end: 15099, file: 'sorghum_bicolor_1_11080-15099.plus.txt' }),
  SORBI_3001G000700: Object.freeze({ region: '1', start: 53681, end: 63405, file: 'sorghum_bicolor_1_53681-63405.plus.txt' })
});

function windowSeq(name) {
  return readWindow(WINDOWS[name].file);
}

let genesCache = null;
function genes() {
  if (!genesCache) genesCache = JSON.parse(fs.readFileSync(path.join(DIR, 'genes.json'), 'utf8'));
  return genesCache;
}

function geneDoc(id) {
  return Object.prototype.hasOwnProperty.call(genes(), id) ? JSON.parse(JSON.stringify(genes()[id])) : null;
}

// In-memory sequence source with the sequence.js contract.
// files: {fastaPath: {regionName: {length, windows: [{start, seq (plus strand)}]}}}
function sequenceStub(files) {
  const calls = [];
  function lookup(fastaPath, region) {
    const f = Object.prototype.hasOwnProperty.call(files, fastaPath) ? files[fastaPath] : null;
    return f && Object.prototype.hasOwnProperty.call(f, region) ? f[region] : null;
  }
  return {
    calls: calls,
    regionLength: async function (fastaPath, region) {
      const r = lookup(fastaPath, region);
      return r ? r.length : undefined;
    },
    fetch: async function (fastaPath, region, start, end, strand) {
      calls.push({ fastaPath: fastaPath, region: region, start: start, end: end, strand: strand });
      const r = lookup(fastaPath, region);
      if (!r) throw new PrimerHttpError(404, 'UNKNOWN_REGION', 'stub: unknown region', { region: String(region) });
      if (end > r.length) throw new PrimerHttpError(400, 'REGION_OUT_OF_BOUNDS', 'stub: beyond region end', { region: region });
      const w = r.windows.find(function (x) { return start >= x.start && end < x.start + x.seq.length; });
      if (!w) throw new Error('stub: ' + region + ':' + start + '-' + end + ' is outside the captured windows');
      const seq = w.seq.slice(start - w.start, end - w.start + 1);
      return strand === -1 ? revcomp(seq) : seq;
    }
  };
}

function resolvedStub(overrides) {
  return Object.assign({
    system_name: 'sorghum_bicolor',
    taxon_id: 4558006,
    display_name: 'Sb bicolor BTx623 v3',
    map_id: 'GCA_000003195.3',
    prefix: 'Sb',
    dir: '/fake/sorghum_bicolor',
    fasta: { dna: SB_DNA, dna_sm: null },
    blastdb: { dna: '/fake/sorghum_bicolor/Sb.dna.toplevel', cdna: null },
    repeat_masking: 'unmasked_copy',
    total_bases: 708735318,
    num_sequences: 1,
    fingerprint: 'f'.repeat(40),
    warnings: []
  }, overrides || {});
}

// Config from DEFAULTS + overrides only (no yaml, no env).
function cfg(overrides) {
  return config._build({ env: {}, fileConfig: {}, overrides: overrides || {} }).config;
}

const silentLog = Object.freeze({ error: function () {}, warn: function () {}, info: function () {}, log: function () {} });

function sorghumSequence() {
  const files = {};
  files[SB_DNA] = {
    '1': {
      length: SB_CHR1_LENGTH,
      windows: Object.keys(WINDOWS).map(function (k) { return { start: WINDOWS[k].start, seq: windowSeq(k) }; })
    }
  };
  return sequenceStub(files);
}

// Offline deps for template.buildTemplate / design.design on the captured sorghum genes.
function sorghumDeps(extra) {
  return Object.assign({
    cfg: cfg(),
    log: silentLog,
    findGene: async function (id) { return geneDoc(id); },
    resolve: async function (name) {
      if (name !== 'sorghum_bicolor') throw new PrimerHttpError(404, 'UNKNOWN_GENOME', 'stub: unknown genome', { system_name: name });
      return resolvedStub();
    },
    sequence: sorghumSequence()
  }, extra || {});
}

function primer3Text(name) {
  return fs.readFileSync(path.join(PRIMER3_DIR, name), 'utf8');
}

function primer3Tags(name) {
  return boulder.parse(primer3Text(name));
}

// SORBI_3001G000200.1 cDNA (1982 nt) from the captured Primer3 input record.
function cdna200() {
  return primer3Tags('transcript_junction.input.txt').SEQUENCE_TEMPLATE;
}

function transcriptSegments200() {
  return JSON.parse(primer3Text('fixtures.json')).transcript_junction.segments;
}

// Plus-strand sequence of genomic blocks [{start, end}] from a captured window; oriented to `strand`.
function blocksSeq(windowName, blocks, strand) {
  const w = WINDOWS[windowName];
  const seq = windowSeq(windowName);
  const plus = blocks.map(function (b) { return seq.slice(b.start - w.start, b.end - w.start + 1); }).join('');
  return strand === -1 ? revcomp(plus) : plus;
}

// A Primer3 record that makes primer3_core evaluate exactly this designed pair on `seq` (SEQUENCE_PRIMER +
// SEQUENCE_PRIMER_REVCOMP, PRIMER_PICK_ANYWAY), so its own PRIMER_PAIR_0_PRODUCT_TM can be compared with the
// design's product_tm. params: salt_monovalent, salt_divalent, dntp_conc (Primer3 defaults when absent).
function forcedPairRecord(seq, pair, params) {
  params = params || {};
  const minLen = Math.min(pair.left.len, pair.right.len);
  const tags = {
    SEQUENCE_ID: 'forced_pair',
    SEQUENCE_TEMPLATE: seq,
    SEQUENCE_PRIMER: pair.left.seq,
    SEQUENCE_PRIMER_REVCOMP: pair.right.seq,
    PRIMER_TASK: 'generic',
    PRIMER_PICK_ANYWAY: 1,
    PRIMER_FIRST_BASE_INDEX: 1,
    PRIMER_NUM_RETURN: 1,
    PRIMER_OPT_SIZE: minLen,
    PRIMER_MIN_SIZE: minLen,
    PRIMER_MAX_SIZE: Math.max(pair.left.len, pair.right.len),
    PRIMER_PRODUCT_SIZE_RANGE: (pair.product_size - 1) + '-' + (pair.product_size + 1),
    PRIMER_PRODUCT_MIN_TM: 0,
    PRIMER_PRODUCT_MAX_TM: 150
  };
  if (params.salt_monovalent !== undefined) tags.PRIMER_SALT_MONOVALENT = params.salt_monovalent;
  if (params.salt_divalent !== undefined) tags.PRIMER_SALT_DIVALENT = params.salt_divalent;
  if (params.dntp_conc !== undefined) tags.PRIMER_DNTP_CONC = params.dntp_conc;
  return tags;
}

// Masked bases of runs [[start, length]] inside [start, end].
function maskedBasesIn(runs, start, end) {
  return (runs || []).reduce(function (n, r) {
    const s = Math.max(start, r[0]);
    const e = Math.min(end, r[0] + r[1] - 1);
    return n + Math.max(0, e - s + 1);
  }, 0);
}

module.exports = {
  ROOT,
  forcedPairRecord,
  maskedBasesIn,
  SB_DNA,
  SB_CHR1_LENGTH,
  WINDOWS,
  windowSeq,
  genes,
  geneDoc,
  sequenceStub,
  resolvedStub,
  cfg,
  silentLog,
  sorghumSequence,
  sorghumDeps,
  primer3Text,
  primer3Tags,
  cdna200,
  transcriptSegments200,
  blocksSeq
};
