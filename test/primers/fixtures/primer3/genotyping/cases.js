'use strict';

// The genotyping design examples of spec §2.9-§2.10 as an offline world, and the replay of their recorded Primer3
// runs (§7.4). Shared by record_genotyping.js (records with the real primer3_core) and the unit tests (replay).
//
// Inputs (read-only, committed):
//   ../../design/sorghum_bicolor_1_10500-12100.plus.txt: sorghum_bicolor chromosome 1, 1:10500-12100 (provenance in
//     test/primers/unit/variation_normalize.test.js); it covers every example template window
//   ../../variation/overlap_1_<start>-<end>.json: live Ensembl 115 overlap bodies of the §4.3 template windows
// variantResolver stands in for variation/index.js: the canonical entry comes from normalize.recordsToEntries over the
// window's overlap records (a manual variant is merged with them, §3.9), and every other entry is a neighbour.

const fs = require('fs');
const path = require('path');
const crypto = require('crypto');
const stubs = require('../../design/stubs');

const ROOT = stubs.ROOT;
const boulder = require(path.join(ROOT, 'api/helpers/primers/boulder'));
const normalize = require(path.join(ROOT, 'api/helpers/primers/variation/normalize'));

const DIR = __dirname;
const FIX = path.join(DIR, '..', '..');
const INDEX = path.join(DIR, 'index.json');
const WINDOW = Object.freeze({ region: '1', start: 10500, end: 12100, file: 'sorghum_bicolor_1_10500-12100.plus.txt' });

// primers.genotyping of spec §3.1, and the primers.variation keys the design reads.
const GENOTYPING = Object.freeze({
  template_flank: 400,
  num_return_per_run: 20,
  num_sets_default: 6,
  max_sets: 10,
  max_scored_per_orientation: 8,
  max_primer3_runs: 54,
  max_thermo_calls: 272,
  thermo_concurrency: 4,
  thermo_timeout_ms: 5000,
  guard_gap: 10,
  mask_exempt_pad: 36,
  as_min_tm: 52,
  as_min_gc: 15,
  structure_warn_th: 47,
  structure_high_th: 55,
  as_tm_diff_warn: 1.0,
  common_tm_low: -1.0,
  common_tm_high: 3.0,
  neighbour_3p_window: 5,
  dense_window: 30,
  dense_count: 2,
  check_max_sets: 5,
  check_max_unique_primers: 13
});
const VARIATION = Object.freeze({ max_allele_length: 50, max_shift: 1000, ems_source_pattern: '^EMS_' });

// The §2.9-§2.10 requests. overlap: the template window whose Ensembl records supply the entry and its neighbours.
const CASES = Object.freeze({
  rs871475760_kasp: Object.freeze({
    section: '2.9', key: '1:11109:C:A', overlap: '10709-11509',
    body: { system_name: 'sorghum_bicolor', variant: { id: 'rs871475760', alt: 'A' }, assay: { type: 'kasp', num_sets: 2 } }
  }),
  tmp_1_11193_C_T_kasp: Object.freeze({
    section: '2.10(a)', key: '1:11193:C:T', overlap: '10793-11593',
    body: { system_name: 'sorghum_bicolor', variant: { id: 'tmp_1_11193_C_T' }, assay: { num_sets: 2 } }
  }),
  rs5413864115_kasp: Object.freeze({
    section: '2.10(b)', key: '1:11282:CA:C', overlap: '10883-11685',
    body: { system_name: 'sorghum_bicolor', variant: { region: '1', position: 11283, ref: 'A', alt: '-' }, assay: { num_sets: 2 } }
  }),
  tmp_1_11502_C_CGT_kasp: Object.freeze({
    section: '2.10(c)', key: '1:11502:C:CGT', overlap: '11102-11903',
    body: { system_name: 'sorghum_bicolor', variant: { id: 'tmp_1_11502_C_CGT' }, assay: { num_sets: 1 } }
  }),
  rs871475760_as_pcr: Object.freeze({
    section: '2.10(d)', key: '1:11109:C:A', overlap: '10709-11509',
    body: { system_name: 'sorghum_bicolor', variant: { id: 'rs871475760' }, assay: { type: 'as_pcr', num_sets: 1 } }
  })
});

function clone(v) {
  return JSON.parse(JSON.stringify(v));
}

// Config from config.js DEFAULTS plus §3.1 genotyping and variation (no yaml, no env). extra: more overrides.
function cfg(genotyping, extra) {
  return stubs.cfg(Object.assign({
    genotyping: Object.assign({}, GENOTYPING, genotyping || {}),
    variation: Object.assign({}, VARIATION)
  }, extra || {}));
}

let windowCache = null;
function windowSeq() {
  if (windowCache === null) windowCache = fs.readFileSync(path.join(FIX, 'design', WINDOW.file), 'utf8').trim();
  return windowCache;
}

function genome() {
  return normalize.sequenceWindow(windowSeq(), WINDOW.start, stubs.SB_CHR1_LENGTH);
}

function overlapRecords(name) {
  return JSON.parse(fs.readFileSync(path.join(FIX, 'variation', 'overlap_1_' + name + '.json'), 'utf8'));
}

// deps.resolveVariant for one case -> {assembly, variant, neighbours {data, entries}, warnings, variation_source}
function variantResolver(name) {
  const c = CASES[name];
  return async function (req) {
    let records = overlapRecords(c.overlap);
    if (req.variant_input === 'manual') records = [normalize.parseManual(req.variant, VARIATION)].concat(records);
    const opts = Object.assign({ region: WINDOW.region }, VARIATION);
    if (req.variant_input === 'id') opts.requested_id = req.variant.id;
    const entries = normalize.recordsToEntries(records, genome(), opts).entries;
    const variant = entries.find(function (e) { return e.key === c.key; });
    if (!variant || (req.variant_input === 'id' && variant.ids.indexOf(req.variant.id) < 0)) {
      throw new Error('cases: ' + name + ' does not resolve to ' + c.key);
    }
    return {
      assembly: stubs.resolvedStub(),
      variant: variant,
      neighbours: { data: 'ensembl', entries: entries.filter(function (e) { return e.key !== c.key; }) },
      warnings: [],
      variation_source: 'ensembl 115'
    };
  };
}

function sequence() {
  const files = {};
  files[stubs.SB_DNA] = { '1': { length: stubs.SB_CHR1_LENGTH, windows: [{ start: WINDOW.start, seq: windowSeq() }] } };
  return stubs.sequenceStub(files);
}

// Offline designGenotyping deps for one case; extra overrides (primer3, cfg, semaphore, scoreCandidate ...).
function deps(name, extra) {
  const order = [];
  return Object.assign({
    cfg: cfg(),
    log: stubs.silentLog,
    sequence: sequence(),
    resolveVariant: variantResolver(name),
    semaphore: {
      order: order,
      acquire: async function () {
        order.push('acquire');
        return function () { order.push('release'); };
      }
    }
  }, extra || {});
}

function body(name) {
  return clone(CASES[name].body);
}

// ---- recordings ----------------------------------------------------------------------------------------

function sha256(s) {
  return crypto.createHash('sha256').update(s).digest('hex');
}

// primer3.run's result for one recorded stdout.
function resultOf(stdout) {
  const r = boulder.parseRecord(stdout);
  return {
    tags: r.tags,
    error: r.tags.PRIMER_ERROR || null,
    warning: r.tags.PRIMER_WARNING || null,
    exitCode: 0,
    elapsedMs: 0,
    stdoutBytes: Buffer.byteLength(stdout),
    stderrTail: ''
  };
}

function loadIndex() {
  return JSON.parse(fs.readFileSync(INDEX, 'utf8'));
}

function readRecording(file) {
  return fs.readFileSync(path.join(DIR, file), 'utf8');
}

// A primer3 {run, version} that answers only byte-identical recorded inputs (sha256 of the serialized record) and
// rejects anything else with code UNRECORDED_PRIMER3_INPUT. calls: receives {name, hash, input} per run.
function recordedPrimer3(calls) {
  const index = loadIndex();
  return {
    run: async function (record) {
      const input = typeof record === 'string' ? record : boulder.serialize(record);
      const hash = sha256(input);
      const hit = Object.prototype.hasOwnProperty.call(index.runs, hash) ? index.runs[hash] : null;
      if (!hit) {
        const e = new Error('unrecorded Primer3 input ' + hash + ' (SEQUENCE_ID=' + boulder.parse(input).SEQUENCE_ID + ')');
        e.code = 'UNRECORDED_PRIMER3_INPUT';
        throw e;
      }
      if (calls) calls.push({ name: hit.name, hash: hash, input: input });
      return resultOf(readRecording(hit.output));
    },
    version: async function () { return index.primer3_version; }
  };
}

module.exports = {
  DIR,
  INDEX,
  WINDOW,
  GENOTYPING,
  VARIATION,
  CASES,
  cfg,
  genome,
  windowSeq,
  overlapRecords,
  variantResolver,
  sequence,
  deps,
  body,
  sha256,
  resultOf,
  loadIndex,
  readRecording,
  recordedPrimer3
};
