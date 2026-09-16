'use strict';

// Primer design / check configuration: DEFAULTS, merged with the `primers:` block of
// config/default.yaml (config lib), then environment overrides. The result is deep-frozen.
//
// basePath comes from api/swagger/swagger.yaml (js-yaml safeLoad) so the API and the
// separate check worker process derive the same site_key without talking to each other.

require('./node_compat'); // must precede `config` (config@1.x calls util.isArray etc.)

const fs = require('fs');
const path = require('path');

const REPO_ROOT = path.resolve(__dirname, '..', '..', '..');
const SWAGGER_YAML = path.join(REPO_ROOT, 'api', 'swagger', 'swagger.yaml');

if (!process.env.NODE_CONFIG_DIR) process.env.NODE_CONFIG_DIR = path.join(REPO_ROOT, 'config');

const DEFAULTS = deepFreeze({
  enabled: true,
  site_key: null,
  primer3_core: '/home/olson/bin/primer3_core',
  blastn: '/home/olson/bin/blastn',
  blastdbcmd: '/home/olson/bin/blastdbcmd',
  fasta_root: '/scratch/olson/fasta',
  tmp_dir: '/tmp',
  assembly_overrides: {},
  repeat_masking_overrides: {},
  design: {
    deadline_ms: 45000,
    primer3_timeout_ms: 30000,
    max_concurrent: 4,
    max_waiting: 16,
    wait_timeout_ms: 10000,
    max_template_length: 50000,
    max_flank: 10000,
    max_fetch_length: 2000000,
    max_stdout_bytes: 2000000
  },
  repeat_mask: {
    enabled: true,
    min_hsp_len: 50,
    perc_identity: 85,
    evalue: 1.0e-10,
    min_depth: 3,
    threads: 2,
    chunk: 20000,
    exon_pad: 100,
    timeout_ms: 20000,
    cache_entries: 200
  },
  check: {
    store: 'redis',
    redis_url: 'redis://localhost:6380/1',
    global_prefix: 'primers:global:',
    global_max_jobs: 2,
    pangenome_max_jobs: 1,
    local_max_jobs: 2,
    spec_job_procs: 4,
    pangenome_job_procs: 8,
    nice: 10,
    max_queued: 50,
    max_pairs: 10,
    max_unique_primers: 20,
    max_genomes: 150,
    max_job_cpu_s: 6000,
    cpu_s_per_primer_gb: { ws5: 5.2, ws6: 2.2, ws7: 1.2 },
    reference_word_size: 5,
    pangenome_word_size: 6,
    pangenome_cpu_factor: 2.0,
    blast_timeout_ms: 600000,
    job_timeout_ms: 1800000,
    heartbeat_ms: 10000,
    stale_ms: 60000,
    queue_sweep_ms: 600000,
    progress_min_interval_ms: 5000,
    max_attempts: 2,
    ttl_active_s: 21600,
    ttl_done_s: 86400,
    ttl_error_s: 3600,
    max_finished_jobs: 500,
    max_result_bytes: 5000000,
    max_candidates_per_pair: 5000,
    max_realign_sites_per_genome: 20000,
    realign_concurrency: 32,
    max_offtargets_listed: 100,
    max_annotated_amplicons: 200,
    defaults: {
      max_product_size: 4000,
      ignore_mismatches: 6,
      min_total_mismatches: 2,
      min_3p_mismatches: 2,
      three_prime_window: 5,
      include_unlikely: false,
      repeat_site_threshold: 5,
      max_amplifying_mismatches: 3
    }
  }
});

// [ENV_NAME, dotted config path, type]. Read here, not via custom-environment-variables.yaml,
// so a malformed value is ignored with a warning instead of crashing startup.
const ENV_OVERRIDES = [
  ['PRIMERS_ENABLED', 'enabled', 'bool'],
  ['PRIMERS_SITE_KEY', 'site_key', 'string'],
  ['PRIMER3_CORE', 'primer3_core', 'string'],
  ['BLASTN', 'blastn', 'string'],
  ['BLASTDBCMD', 'blastdbcmd', 'string'],
  ['PRIMERS_FASTA_ROOT', 'fasta_root', 'string'],
  ['PRIMERS_JOB_STORE', 'check.store', 'store'],
  ['PRIMERS_REDIS_URL', 'check.redis_url', 'string'],
  ['PRIMERS_GLOBAL_PREFIX', 'check.global_prefix', 'string'],
  ['PRIMERS_GLOBAL_MAX_JOBS', 'check.global_max_jobs', 'int'],
  ['PRIMERS_MAX_QUEUED', 'check.max_queued', 'int']
];

const UNSAFE_KEYS = new Set(['__proto__', 'constructor', 'prototype']);

function isPlainObject(v) {
  return v !== null && typeof v === 'object' && !Array.isArray(v) &&
    (Object.getPrototypeOf(v) === Object.prototype || Object.getPrototypeOf(v) === null);
}

function clone(v) {
  if (Array.isArray(v)) return v.map(clone);
  if (isPlainObject(v)) {
    const out = {};
    Object.keys(v).forEach(function (k) { if (!UNSAFE_KEYS.has(k)) out[k] = clone(v[k]); });
    return out;
  }
  return v;
}

// Objects merge recursively; arrays, scalars and null replace.
function deepMerge(target, src) {
  if (!isPlainObject(src)) return target;
  Object.keys(src).forEach(function (k) {
    if (UNSAFE_KEYS.has(k)) return;
    if (isPlainObject(src[k]) && isPlainObject(target[k])) target[k] = deepMerge(target[k], src[k]);
    else target[k] = clone(src[k]);
  });
  return target;
}

function deepFreeze(o) {
  if (o !== null && typeof o === 'object' && !Object.isFrozen(o)) {
    Object.freeze(o);
    Object.keys(o).forEach(function (k) { deepFreeze(o[k]); });
  }
  return o;
}

function setPath(obj, dotted, value) {
  const parts = dotted.split('.');
  let cur = obj;
  for (let i = 0; i < parts.length - 1; i++) {
    if (!isPlainObject(cur[parts[i]])) cur[parts[i]] = {};
    cur = cur[parts[i]];
  }
  cur[parts[parts.length - 1]] = value;
}

function applyEnv(cfg, env, warnings) {
  ENV_OVERRIDES.forEach(function (row) {
    const name = row[0];
    const raw = env[name];
    if (raw === undefined || raw === '') return;
    const s = String(raw).trim();
    switch (row[2]) {
      case 'bool':
        setPath(cfg, row[1], /^(1|true|yes)$/i.test(s));
        break;
      case 'int':
        if (/^\d+$/.test(s) && Number.isSafeInteger(Number(s))) setPath(cfg, row[1], Number(s));
        else warnings.push(name + ' is not a non-negative integer; ignored');
        break;
      case 'store':
        if (s === 'redis' || s === 'memory') setPath(cfg, row[1], s);
        else warnings.push(name + ' must be redis or memory; ignored');
        break;
      default:
        setPath(cfg, row[1], s);
    }
  });
}

function loadFileConfig() {
  const config = require('config');
  return config.has('primers') ? config.util.toObject(config.get('primers')) : {};
}

// Pure builder (exported for tests): {env, fileConfig, overrides} -> {config (frozen), warnings}.
function build(opts) {
  opts = opts || {};
  const env = opts.env || process.env;
  const warnings = [];
  const fileCfg = opts.fileConfig !== undefined ? opts.fileConfig : loadFileConfig();
  if (fileCfg !== null && fileCfg !== undefined && !isPlainObject(fileCfg)) warnings.push('primers config block is not a mapping; ignored');
  const cfg = deepMerge(clone(DEFAULTS), isPlainObject(fileCfg) ? fileCfg : {});
  applyEnv(cfg, env, warnings);
  if (opts.overrides) deepMerge(cfg, opts.overrides);
  return { config: deepFreeze(cfg), warnings: warnings };
}

let cached = null;
let basePathOverride = null;
let basePathFromYaml = null;

function get() {
  if (!cached) {
    const r = build();
    r.warnings.forEach(function (w) { console.warn('primers config: ' + w); });
    cached = r.config;
  }
  return cached;
}

function validBasePath(bp) {
  return typeof bp === 'string' && /^\/[A-Za-z0-9_.\-/]*$/.test(bp);
}

function readBasePath(file) {
  const yaml = require('js-yaml');
  const doc = yaml.safeLoad(fs.readFileSync(file || SWAGGER_YAML, 'utf8'));
  const bp = doc && doc.basePath;
  if (!validBasePath(bp)) throw new Error('primers config: swagger.yaml has no usable basePath');
  return bp;
}

function basePath() {
  if (basePathOverride) return basePathOverride;
  if (!basePathFromYaml) basePathFromYaml = readBasePath();
  return basePathFromYaml;
}

function setBasePath(bp) {
  if (!validBasePath(bp)) throw new TypeError('basePath must be a string starting with "/"');
  basePathOverride = bp;
}

// PRIMERS_SITE_KEY || primers.site_key || "<basePath without leading slash>:<mongo db>".
function siteKey() {
  const cfg = get();
  if (cfg.site_key !== null && cfg.site_key !== undefined && String(cfg.site_key) !== '') return String(cfg.site_key);
  const mongo = require('gramene-mongodb-config');
  if (!mongo || typeof mongo.getMongoConfig !== 'function') {
    throw new Error('primers config: gramene-mongodb-config has no getMongoConfig(); set PRIMERS_SITE_KEY');
  }
  const db = mongo.getMongoConfig().db;
  if (!db) throw new Error('primers config: getMongoConfig().db is empty; set PRIMERS_SITE_KEY');
  return basePath().replace(/^\/+/, '') + ':' + db;
}

// Tests: _setForTests(overrides) rebuilds (defaults + yaml + env + overrides);
// _setForTests(null) drops the cache and any basePath override.
function _setForTests(overrides) {
  if (overrides === null || overrides === undefined) {
    cached = null;
    basePathOverride = null;
    basePathFromYaml = null;
    return undefined;
  }
  cached = build({ overrides: overrides }).config;
  return cached;
}

module.exports = {
  get,
  siteKey,
  basePath,
  setBasePath,
  _setForTests,
  _build: build,
  readBasePath,
  DEFAULTS,
  ENV_OVERRIDES: deepFreeze(ENV_OVERRIDES.map(function (r) { return r.slice(); })),
  REPO_ROOT,
  SWAGGER_YAML
};
