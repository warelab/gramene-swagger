'use strict';

// The genotyping design examples of spec §2.9-§2.10 as an offline world, and the replay of their recorded Primer3 and
// ntthal runs (§7.4). Shared by record_genotyping.js (records with the real binaries) and the unit tests (replay).
//
// Inputs (read-only, committed):
//   ../../design/sorghum_bicolor_1_10500-12100.plus.txt: sorghum_bicolor chromosome 1, 1:10500-12100 (provenance in
//     test/primers/unit/variation_normalize.test.js); it covers every example template window. Bases outside it read
//     as N: only the wide normalization reads of variation/index.js reach there.
//   ../../variation/overlap_1_<start>-<end>.json: live Ensembl 115 overlap bodies of the §4.3 template windows
//   ../../variation/variation_<id>.json: live Ensembl 115 lookups; an id without one (tmp_1_11193_C_T) is looked up from
//     its recorded overlap record, with no synonyms (as test/primers/unit/variation_index.test.js does)
//   ../../thermo/genotyping.json: the recorded ntthal results, keyed by argv
// The variant and its neighbours are resolved by the production path (genotyping/design.js variationResolver:
// variation/index.js resolveDesignVariant and neighboursFor) through the real Ensembl client over a fake fetch serving
// those bodies.

const fs = require('fs');
const path = require('path');
const crypto = require('crypto');
const stubs = require('../../design/stubs');

const ROOT = stubs.ROOT;
const boulder = require(path.join(ROOT, 'api/helpers/primers/boulder'));
const normalize = require(path.join(ROOT, 'api/helpers/primers/variation/normalize'));
const { createVariationClient } = require(path.join(ROOT, 'api/helpers/primers/variation/client'));
const { PrimerHttpError } = require(path.join(ROOT, 'api/helpers/primers/errors'));
const thermoModule = require(path.join(ROOT, 'api/helpers/primers/thermo'));
const gdesign = require(path.join(ROOT, 'api/helpers/primers/genotyping/design'));

const DIR = __dirname;
const FIX = path.join(DIR, '..', '..');
const INDEX = path.join(DIR, 'index.json');
const THERMO = path.join(FIX, 'thermo', 'genotyping.json');
const WINDOW = Object.freeze({ region: '1', start: 10500, end: 12100, file: 'sorghum_bicolor_1_10500-12100.plus.txt' });
const OVERLAP_WINDOWS = Object.freeze(['10709-11509', '10793-11593', '10883-11685', '11102-11903']);

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

// The §2.9-§2.10 requests. overlap: the template window of the example.
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
  if (windowCache === null) windowCache = fs.readFileSync(path.join(FIX, 'design', WINDOW.file), 'utf8').trim().toUpperCase();
  return windowCache;
}

function genome() {
  return normalize.sequenceWindow(windowSeq(), WINDOW.start, stubs.SB_CHR1_LENGTH);
}

function overlapRecords(name) {
  return JSON.parse(fs.readFileSync(path.join(FIX, 'variation', 'overlap_1_' + name + '.json'), 'utf8'));
}

// Every recorded overlap record once (the recorded windows overlap).
let recordsCache = null;
function allRecords() {
  if (recordsCache) return recordsCache;
  const seen = new Set();
  recordsCache = [];
  OVERLAP_WINDOWS.forEach(function (w) {
    overlapRecords(w).forEach(function (r) {
      const k = JSON.stringify([r.id, r.start, r.end, r.alleles]);
      if (!seen.has(k)) {
        seen.add(k);
        recordsCache.push(r);
      }
    });
  });
  return recordsCache;
}

function respond(body, status) {
  return new Response(typeof body === 'string' ? body : JSON.stringify(body), { status: status || 200 });
}

// The recorded Ensembl service as a fetch. mode: 'normal' (default), or 'down' (connection refused).
function ensemblFetch(mode) {
  return async function (url) {
    if (mode === 'down') {
      throw Object.assign(new TypeError('fetch failed'), { cause: Object.assign(new Error('connect ECONNREFUSED'), { code: 'ECONNREFUSED' }) });
    }
    const u = new URL(url);
    let m = /\/overlap\/region\/sorghum_bicolor\/([^/]+):(\d+)-(\d+)$/.exec(u.pathname);
    if (m) {
      const region = decodeURIComponent(m[1]);
      const s = Number(m[2]);
      const e = Number(m[3]);
      return respond(allRecords().filter(function (r) {
        return String(r.seq_region_name) === region && Math.min(r.start, r.end) <= e && Math.max(r.start, r.end) >= s;
      }));
    }
    m = /\/variation\/sorghum_bicolor\/([^/]+)$/.exec(u.pathname);
    if (m) {
      const id = decodeURIComponent(m[1]);
      const file = path.join(FIX, 'variation', 'variation_' + id + '.json');
      if (fs.existsSync(file)) return respond(fs.readFileSync(file, 'utf8'));
      const r = allRecords().find(function (x) { return x.id === id; });
      if (r) {
        return respond({ name: id, synonyms: [], most_severe_consequence: r.consequence_type,
          mappings: [{ seq_region_name: r.seq_region_name, start: r.start, end: r.end, allele_string: r.alleles.join('/') }] });
      }
      return respond(fs.readFileSync(path.join(FIX, 'variation', 'variation_not_found.json'), 'utf8'), 400);
    }
    return respond('<html><body>404 Not Found</body></html>', 404);
  };
}

// sequence.js contract over the recorded window, N elsewhere on region 1.
function sequence() {
  const seq = windowSeq();
  const calls = [];
  return {
    calls: calls,
    regionLength: async function (fasta, region) { return fasta === stubs.SB_DNA && region === '1' ? stubs.SB_CHR1_LENGTH : undefined; },
    fetch: async function (fasta, region, start, end) {
      calls.push({ region: region, start: start, end: end });
      if (fasta !== stubs.SB_DNA || region !== '1') throw new PrimerHttpError(404, 'UNKNOWN_REGION', 'stub: unknown region', { region: String(region) });
      let out = '';
      for (let p = start; p <= end; p++) out += p >= WINDOW.start && p <= WINDOW.end ? seq[p - WINDOW.start] : 'N';
      return out;
    }
  };
}

async function resolveAssembly(name) {
  if (name !== 'sorghum_bicolor') throw new PrimerHttpError(404, 'UNKNOWN_GENOME', 'stub: unknown genome', { system_name: name });
  return stubs.resolvedStub();
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

let thermoCache = null;
function loadThermo() {
  if (!thermoCache) thermoCache = JSON.parse(fs.readFileSync(THERMO, 'utf8'));
  return thermoCache;
}

// The spawn of thermo.createThermo over the recorded ntthal results: only a recorded argv is answered; anything else is
// rejected with code UNRECORDED_THERMO_INPUT. calls: receives each argv (' '-joined) spawned.
function recordedNtthal(calls) {
  const fixture = loadThermo();
  return async function (bin, args) {
    const key = args.join(' ');
    if (!Object.prototype.hasOwnProperty.call(fixture.calls, key)) {
      const e = new Error('unrecorded ntthal input: ' + key);
      e.code = 'UNRECORDED_THERMO_INPUT';
      throw e;
    }
    if (calls) calls.push(key);
    return { stdout: fixture.calls[key] + '\n' };
  };
}

// A thermo.js instance (real parsing, memo, pool and counting) over the recorded ntthal results.
function recordedThermo(config, calls) {
  return thermoModule.createThermo({ cfg: config || cfg(), log: stubs.silentLog, spawn: recordedNtthal(calls) });
}

// ---- deps ----------------------------------------------------------------------------------------------

// Offline designGenotyping deps for one case: the production resolver over the recorded Ensembl, the recorded thermo
// and a semaphore that records acquire/release. extra overrides (primer3, cfg, semaphore, scoreCandidate, thermo,
// ensembl: 'down' ...). The variation client, and so its caches, belong to this deps object alone.
function deps(name, extra) {
  const order = [];
  const x = extra || {};
  const out = Object.assign({
    cfg: cfg(),
    log: stubs.silentLog,
    sequence: sequence(),
    resolve: resolveAssembly,
    semaphore: {
      order: order,
      acquire: async function () {
        order.push('acquire');
        return function () { order.push('release'); };
      }
    }
  }, x);
  delete out.ensembl;
  if (!out.variationClient) {
    out.variationClient = createVariationClient({ cfg: out.cfg, log: stubs.silentLog, fetch: ensemblFetch(x.ensembl || 'normal') });
  }
  if (out.thermo === undefined) out.thermo = recordedThermo(out.cfg);
  return out;
}

// The resolved variant of a request through the same resolver designGenotyping uses by default.
function resolve(req, d) {
  return gdesign.variationResolver(d)(req, {});
}

function body(name) {
  return clone(CASES[name].body);
}

module.exports = {
  DIR,
  INDEX,
  THERMO,
  WINDOW,
  GENOTYPING,
  VARIATION,
  CASES,
  cfg,
  genome,
  windowSeq,
  overlapRecords,
  allRecords,
  ensemblFetch,
  sequence,
  resolveAssembly,
  deps,
  resolve,
  body,
  sha256,
  resultOf,
  loadIndex,
  readRecording,
  recordedPrimer3,
  loadThermo,
  recordedNtthal,
  recordedThermo
};
