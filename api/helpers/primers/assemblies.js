'use strict';

// The single assembly resolver (spec §A.4.1). Maps a catalog genome to the files under
// <fasta_root>/<system_name>/: which assembly prefix to use, its bgzip FASTA and soft-masked copy,
// its BLAST DBs, the repeat-masking status, sizes and a fingerprint of the files.
// Used by design, GET /primers/genomes and check submission; the worker reads job.resolved and
// never re-resolves.
//
// Security: system_name must match ^[a-z0-9_]+$ AND be in the maps catalog before any path.join;
// the joined dir is asserted to sit directly under fasta_root; file names come only from readdir.

const fs = require('fs');
const path = require('path');
const crypto = require('crypto');
const { PrimerHttpError, isPrimerHttpError } = require('./errors');
const genomes = require('./genomes');

const RESOLVE_TTL_MS = 10 * 60 * 1000;
const RESOLVE_CONCURRENCY = 8;
const RESOLVE_CACHE_MAX = 2000;
const REGION_MATCH_MIN_FRACTION = 0.8;
const MASK_CACHE_MAX = 1000;
// Soft-mask sampling: windows of 20 kb from the 3 longest sequences, at their midpoints first
// (pericentromeric, repeat-rich), then at 1/4 and 3/4. Stops at the first lowercase base.
const MASK_SAMPLE = Object.freeze({ sequences: 3, window: 20000, fractions: Object.freeze([0.5, 0.25, 0.75]) });
const REPEAT_MASKING = Object.freeze(['soft_masked', 'unmasked_copy', 'absent']);

const DNA_FASTA_RE = /^(.+)\.dna\.toplevel\.fa\.gz$/;
const DNA_BLASTDB_RE = /^(.+)\.dna\.toplevel\.(?:nal|nin)$/;
const VOLUME_NIN_RE = /\.\d\d\.nin$/;

function hasOwn(obj, key) {
  return !!obj && Object.prototype.hasOwnProperty.call(obj, key);
}

function escapeRegExp(s) {
  return s.replace(/[.*+?^${}()|[\]\\]/g, '\\$&');
}

function logSafe(log, level, msg) {
  try {
    const fn = (log && (log[level] || log.error || log.log)) || function () {};
    fn.call(log, msg);
  } catch (e) { /* logging must never throw */ }
}

function getConfig(deps) {
  if (deps && deps.cfg) return deps.cfg;
  return require('./config').get();
}

function isMissing(err) {
  return !!err && (err.code === 'ENOENT' || err.code === 'ENOTDIR');
}

async function readdirNames(fsp, dir) {
  try {
    return (await fsp.readdir(dir)).filter(function (n) { return typeof n === 'string'; });
  } catch (err) {
    if (isMissing(err)) return [];
    throw err;
  }
}

async function statFile(fsp, file) {
  if (!file) return null;
  try {
    const st = await fsp.stat(file);
    return st.isFile() ? st : null;
  } catch (err) {
    if (isMissing(err)) return null;
    throw err;
  }
}

function statToken(st) {
  return st ? st.size + ':' + Math.trunc(st.mtimeMs) : '-';
}

// ---- pure helpers -------------------------------------------------------------------------------

// parseFai(text) -> {lengths: Map(name -> length), total_bases, num_sequences}
function parseFai(text) {
  const lengths = new Map();
  let total = 0;
  String(text).split('\n').forEach(function (raw) {
    const line = raw.replace(/\r$/, '');
    if (!line) return;
    const cols = line.split('\t');
    const len = Number(cols[1]);
    if (cols.length < 2 || !cols[0] || !Number.isSafeInteger(len) || len < 0 || lengths.has(cols[0])) return;
    lengths.set(cols[0], len);
    total += len;
  });
  return { lengths: lengths, total_bases: total, num_sequences: lengths.size };
}

// Identity of a .fai as a (name, length) set, independent of line order and byte offsets.
function faiIdentity(fai) {
  if (!fai) return null;
  if (!fai.identity) {
    const entries = Array.from(fai.lengths.entries()).map(function (e) { return e[0] + '\t' + e[1]; }).sort();
    fai.identity = crypto.createHash('sha1').update(entries.join('\n')).digest('hex');
  }
  return fai.identity;
}

// Map regions that are bins, not sequences: UNANCHORED (in 124 of the 128 sorghum_v11 maps) collects every
// unanchored scaffold under one name, so it is never in a FASTA index. Bins are left out of the region-match
// score and of the ASSEMBLY_MISMATCH denominator.
const SYNTHETIC_REGION_RE = /^(?:UNANCHORED|UNPLACED|UNASSIGNED|UNLOCALIZED)$/i;

function isSyntheticRegion(name) {
  return SYNTHETIC_REGION_RE.test(String(name));
}

// Number of map regions that are real sequences (not synthetic bins).
function realRegionCount(regions) {
  if (!regions) return 0;
  let n = 0;
  for (let i = 0; i < regions.names.length; i++) {
    if (!isSyntheticRegion(regions.names[i])) n++;
  }
  return n;
}

// Number of real (non-bin) map regions found in the .fai with the same name and length.
function regionMatches(fai, regions) {
  if (!fai || !regions) return 0;
  let n = 0;
  for (let i = 0; i < regions.names.length; i++) {
    if (isSyntheticRegion(regions.names[i])) continue;
    if (fai.lengths.get(String(regions.names[i])) === Number(regions.lengths[i])) n++;
  }
  return n;
}

// scanCandidates(topNames, dnaNames) -> {prefixes (sorted), fasta Set, blast Set, top Set, dna Set, topNames}
// Candidates: prefixes of dna/*.dna.toplevel.fa.gz that also have .fai and .gzi, union prefixes of
// *.dna.toplevel.nal|.nin (BLAST volume files such as X.dna.toplevel.00.nin are excluded).
function scanCandidates(topNames, dnaNames) {
  const top = new Set(topNames);
  const dna = new Set(dnaNames);
  const fasta = new Set();
  const blast = new Set();
  dnaNames.forEach(function (name) {
    const m = DNA_FASTA_RE.exec(name);
    if (m && dna.has(name + '.fai') && dna.has(name + '.gzi')) fasta.add(m[1]);
  });
  topNames.forEach(function (name) {
    if (VOLUME_NIN_RE.test(name)) return;
    const m = DNA_BLASTDB_RE.exec(name);
    if (m) blast.add(m[1]);
  });
  const prefixes = Array.from(new Set(Array.from(fasta).concat(Array.from(blast)))).sort();
  return { prefixes: prefixes, fasta: fasta, blast: blast, top: top, dna: dna, topNames: topNames.slice() };
}

// blastDbFiles(scan, prefix, 'dna.toplevel' | 'cdna.all') -> {base, alias, index, seq, volumes[]} | null
// A .nal alias is preferred over .nin; a .nin DB also needs its .nsq.
function blastDbFiles(scan, prefix, kind) {
  const base = prefix + '.' + kind;
  if (scan.top.has(base + '.nal')) {
    const volume = new RegExp('^' + escapeRegExp(base) + '\\.\\d{2,}\\.nsq$');
    return {
      base: base, alias: base + '.nal', index: null, seq: null,
      volumes: scan.topNames.filter(function (n) { return volume.test(n); }).sort()
    };
  }
  if (scan.top.has(base + '.nin') && scan.top.has(base + '.nsq')) {
    return { base: base, alias: null, index: base + '.nin', seq: base + '.nsq', volumes: [] };
  }
  return null;
}

// maskSampleWindows(fai[, sample]) -> [{name, start, end}] (0-based, end exclusive)
function maskSampleWindows(fai, sample) {
  sample = sample || MASK_SAMPLE;
  const longest = Array.from(fai.lengths.entries())
    .filter(function (e) { return e[1] > 0; })
    .sort(function (a, b) { return b[1] - a[1] || (a[0] < b[0] ? -1 : a[0] > b[0] ? 1 : 0); })
    .slice(0, sample.sequences);
  const out = [];
  const seen = new Set();
  sample.fractions.forEach(function (f) {
    longest.forEach(function (e) {
      const len = e[1];
      const w = Math.min(sample.window, len);
      const start = Math.max(0, Math.min(len - w, Math.floor(len * f - w / 2)));
      const key = e[0] + ':' + start;
      if (seen.has(key)) return;
      seen.add(key);
      out.push({ name: e[0], start: start, end: start + w });
    });
  });
  return out;
}

function deepFreeze(result) {
  Object.freeze(result.fasta);
  Object.freeze(result.blastdb);
  result.warnings.forEach(Object.freeze);
  Object.freeze(result.warnings);
  return Object.freeze(result);
}

function ambiguousError(genome, prefixes, why) {
  return new PrimerHttpError(422, 'AMBIGUOUS_ASSEMBLY',
    'cannot choose between ' + prefixes.length + ' assemblies for ' + genome.system_name + ' (' + why +
    '); set primers.assembly_overrides', { system_name: genome.system_name, candidates: prefixes.slice() });
}

// ---- soft-mask detection ------------------------------------------------------------------------

let maskCache = new Map();

function openBgzipFasta(file) {
  const { BgzipIndexedFasta } = require('@gmod/indexedfasta');
  return new BgzipIndexedFasta({ path: file, faiPath: file + '.fai', gziPath: file + '.gzi' });
}

// -> 'soft_masked' | 'unmasked_copy' | null (could not be sampled)
async function detectSoftMask(file, fai, deps) {
  const fsp = deps.fs || fs.promises;
  const log = deps.log || console;
  let key = null;
  try {
    const st = await statFile(fsp, file);
    if (!st) return null;
    key = file + '\0' + statToken(st);
    if (maskCache.has(key)) return maskCache.get(key);
    const handle = (deps.openFasta || openBgzipFasta)(file);
    let bases = 0;
    let lower = false;
    const windows = maskSampleWindows(fai, deps.mask_sample || MASK_SAMPLE);
    for (let i = 0; i < windows.length && !lower; i++) {
      const w = windows[i];
      const seq = await handle.getSequence(w.name, w.start, w.end);
      if (typeof seq !== 'string') continue;
      bases += seq.length;
      lower = /[a-z]/.test(seq);
    }
    if (bases === 0) throw new Error('no bases could be read');
    const value = lower ? 'soft_masked' : 'unmasked_copy';
    if (maskCache.size >= MASK_CACHE_MAX) maskCache.delete(maskCache.keys().next().value);
    maskCache.set(key, value);
    return value;
  } catch (err) {
    logSafe(log, 'warn', 'primers: could not sample ' + path.basename(file) + ' for soft-masking: ' +
      (err && err.message ? err.message : String(err)));
    return null;
  }
}

// ---- resolution ---------------------------------------------------------------------------------

function makeInspector(fsp, dir, scan) {
  const dnaDir = path.join(dir, 'dna');
  const memo = new Map();
  function once(key, fn) {
    if (!memo.has(key)) memo.set(key, fn());
    return memo.get(key);
  }
  return {
    dir: dir,
    dnaDir: dnaDir,
    fai: function (prefix) {
      return once('fai\0' + prefix, async function () {
        if (!scan.fasta.has(prefix)) return null;
        return parseFai(await fsp.readFile(path.join(dnaDir, prefix + '.dna.toplevel.fa.gz.fai'), 'utf8'));
      });
    },
    stat: function (file) {
      return once('stat\0' + file, function () { return statFile(fsp, file); });
    }
  };
}

// Bytes of sequence in the dna BLAST DB (sum of volume .nsq for a .nal) or null without a DB.
async function blastSeqBytes(scan, prefix, insp) {
  const db = blastDbFiles(scan, prefix, 'dna.toplevel');
  if (!db) return null;
  const files = db.alias ? db.volumes : [db.seq];
  if (!files.length) return null;
  let total = 0;
  for (let i = 0; i < files.length; i++) {
    const st = await insp.stat(path.join(insp.dir, files[i]));
    if (!st) return null;
    total += st.size;
  }
  return total;
}

async function identicalAssemblies(prefixes, scan, insp) {
  if (prefixes.length < 2) return false;
  const faiIds = await Promise.all(prefixes.map(async function (p) { return faiIdentity(await insp.fai(p)); }));
  const seqBytes = await Promise.all(prefixes.map(function (p) { return blastSeqBytes(scan, p, insp); }));
  return faiIds.every(function (id) { return id === faiIds[0]; }) &&
    seqBytes.every(function (b) { return b === seqBytes[0]; }) &&
    (faiIds[0] !== null || seqBytes[0] !== null);
}

// Newest by BLAST DB mtime (.nal, else .nin; FASTA mtime when there is no DB); ties -> last prefix.
async function newestPrefix(prefixes, scan, insp) {
  const stamped = await Promise.all(prefixes.map(async function (p) {
    const db = blastDbFiles(scan, p, 'dna.toplevel');
    const file = db ? path.join(insp.dir, db.alias || db.index)
      : path.join(insp.dnaDir, p + '.dna.toplevel.fa.gz');
    const st = await insp.stat(file);
    return { prefix: p, mtime: st ? st.mtimeMs : -Infinity };
  }));
  stamped.sort(function (a, b) {
    return b.mtime - a.mtime || (a.prefix < b.prefix ? 1 : a.prefix > b.prefix ? -1 : 0);
  });
  return stamped[0].prefix;
}

// Pick order: override -> unique '.'+map_id suffix -> unique best region score -> sole candidate ->
// identical tie (newest + warning) -> 422 AMBIGUOUS_ASSEMBLY.
async function pickPrefix(genome, scan, cfg, insp) {
  const prefixes = scan.prefixes;
  const warnings = [];
  if (prefixes.length === 0) return { prefix: null, rule: 'none', warnings: warnings };

  const overrides = cfg.assembly_overrides;
  if (hasOwn(overrides, genome.system_name) && overrides[genome.system_name] !== null &&
      overrides[genome.system_name] !== undefined && overrides[genome.system_name] !== '') {
    const wanted = overrides[genome.system_name];
    if (typeof wanted === 'string' && prefixes.indexOf(wanted) >= 0) {
      return { prefix: wanted, rule: 'override', warnings: warnings };
    }
    throw ambiguousError(genome, prefixes, 'the configured assembly override matches no assembly files');
  }

  const bySuffix = prefixes.filter(function (p) { return p.endsWith('.' + genome.map_id); });
  if (bySuffix.length === 1) return { prefix: bySuffix[0], rule: 'map_id', warnings: warnings };

  const scores = await Promise.all(prefixes.map(async function (p) {
    return regionMatches(await insp.fai(p), genome.regions);
  }));
  const best = Math.max.apply(null, scores);
  const tied = prefixes.filter(function (p, i) { return scores[i] === best; });
  if (best > 0 && tied.length === 1) return { prefix: tied[0], rule: 'regions', warnings: warnings };

  if (prefixes.length === 1) return { prefix: prefixes[0], rule: 'sole', warnings: warnings };

  if (await identicalAssemblies(tied, scan, insp)) {
    const prefix = await newestPrefix(tied, scan, insp);
    warnings.push({
      code: 'AMBIGUOUS_ASSEMBLY',
      message: tied.length + ' identical assemblies found; using the one with the newest BLAST database'
    });
    return { prefix: prefix, rule: 'newest', warnings: warnings };
  }
  throw ambiguousError(genome, prefixes, 'no unique match for the map');
}

async function fingerprintOf(prefix, insp, scan, dnaDb, cdnaDb) {
  const files = [];
  if (scan.fasta.has(prefix)) {
    const fa = path.join(insp.dnaDir, prefix + '.dna.toplevel.fa.gz');
    files.push(['dna.fa.gz', fa], ['dna.fai', fa + '.fai']);
  } else {
    files.push(['dna.fa.gz', null], ['dna.fai', null]);
  }
  [['dna', dnaDb], ['cdna', cdnaDb]].forEach(function (pair) {
    const label = pair[0];
    const db = pair[1];
    if (!db) {
      files.push([label + '.db', null]);
    } else if (db.alias) {
      files.push([label + '.nal', path.join(insp.dir, db.alias)]);
      db.volumes.forEach(function (v) { files.push([label + '.volume', path.join(insp.dir, v)]); });
    } else {
      files.push([label + '.nin', path.join(insp.dir, db.index)], [label + '.nsq', path.join(insp.dir, db.seq)]);
    }
  });
  const parts = [prefix];
  for (let i = 0; i < files.length; i++) {
    parts.push(files[i][0] + '=' + statToken(files[i][1] ? await insp.stat(files[i][1]) : null));
  }
  return crypto.createHash('sha1').update(parts.join('\n')).digest('hex');
}

async function resolveUncached(genome, cfg, deps) {
  const fsp = deps.fs || fs.promises;
  const log = deps.log || console;
  const root = path.resolve(cfg.fasta_root);
  const dir = path.join(root, genome.system_name);
  if (path.dirname(dir) !== root || path.basename(dir) !== genome.system_name) {
    throw genomes.unknownGenomeError(genome.system_name);
  }
  const names = await Promise.all([readdirNames(fsp, dir), readdirNames(fsp, path.join(dir, 'dna'))]);
  const scan = scanCandidates(names[0], names[1]);
  const insp = makeInspector(fsp, dir, scan);
  const pick = await pickPrefix(genome, scan, cfg, insp);
  const prefix = pick.prefix;
  const result = {
    system_name: genome.system_name,
    taxon_id: genome.taxon_id,
    display_name: genome.display_name,
    map_id: genome.map_id,
    prefix: prefix,
    dir: dir,
    fasta: { dna: null, dna_sm: null },
    blastdb: { dna: null, cdna: null },
    repeat_masking: 'absent',
    total_bases: null,
    num_sequences: null,
    fingerprint: null,
    warnings: pick.warnings.slice()
  };
  if (prefix === null) return deepFreeze(result);

  const fai = await insp.fai(prefix);
  let smFai = null;
  if (fai) {
    result.total_bases = fai.total_bases;
    result.num_sequences = fai.num_sequences;
    // A map with no real regions (only an UNANCHORED bin) cannot be checked: no warning.
    const total = realRegionCount(genome.regions);
    const matched = regionMatches(fai, genome.regions);
    if (total > 0 && matched < REGION_MATCH_MIN_FRACTION * total) {
      result.warnings.push({
        code: 'ASSEMBLY_MISMATCH',
        message: 'only ' + matched + ' of ' + total + ' map regions match the assembly index by name and length'
      });
    }
    result.fasta.dna = path.join(insp.dnaDir, prefix + '.dna.toplevel.fa.gz');
    const smName = prefix + '.dna_sm.toplevel.fa.gz';
    if (scan.dna.has(smName) && scan.dna.has(smName + '.fai') && scan.dna.has(smName + '.gzi')) {
      const smPath = path.join(insp.dnaDir, smName);
      smFai = parseFai(await fsp.readFile(smPath + '.fai', 'utf8'));
      if (faiIdentity(smFai) === faiIdentity(fai)) {
        result.fasta.dna_sm = smPath;
      } else {
        logSafe(log, 'warn', 'primers: ignoring the dna_sm FASTA of ' + genome.system_name + ': its index differs from dna');
      }
    }
  }

  const dnaDb = blastDbFiles(scan, prefix, 'dna.toplevel');
  const cdnaDb = blastDbFiles(scan, prefix, 'cdna.all');
  result.blastdb.dna = dnaDb ? path.join(dir, dnaDb.base) : null;
  result.blastdb.cdna = cdnaDb ? path.join(dir, cdnaDb.base) : null;

  const forced = hasOwn(cfg.repeat_masking_overrides, genome.system_name)
    ? cfg.repeat_masking_overrides[genome.system_name] : undefined;
  if (REPEAT_MASKING.indexOf(forced) >= 0) {
    // The override wins, but a mask cannot come from a dna_sm FASTA that is not there.
    result.repeat_masking = forced !== 'absent' && !result.fasta.dna_sm ? 'absent' : forced;
  } else if (result.fasta.dna_sm) {
    const detected = await detectSoftMask(result.fasta.dna_sm, smFai, deps);
    if (detected) {
      result.repeat_masking = detected;
    } else {
      result.fasta.dna_sm = null;
      result.repeat_masking = 'absent';
    }
  }

  result.fingerprint = await fingerprintOf(prefix, insp, scan, dnaDb, cdnaDb);
  return deepFreeze(result);
}

let resolveCache = new Map();

function resolveCacheKey(root, genome, cfg) {
  const override = hasOwn(cfg.assembly_overrides, genome.system_name) ? String(cfg.assembly_overrides[genome.system_name]) : '';
  const masking = hasOwn(cfg.repeat_masking_overrides, genome.system_name) ? String(cfg.repeat_masking_overrides[genome.system_name]) : '';
  return [root, genome.system_name, genome.map_id, override, masking].join('\0');
}

function pruneResolveCache(now, ttl) {
  resolveCache.forEach(function (entry, key) {
    if (now - entry.at >= ttl) resolveCache.delete(key);
  });
  while (resolveCache.size >= RESOLVE_CACHE_MAX) resolveCache.delete(resolveCache.keys().next().value);
}

// resolve(system_name, deps) -> frozen result (spec §A.4.1 item 5):
//   {system_name, taxon_id, display_name, map_id, prefix, dir, fasta {dna, dna_sm}, blastdb {dna, cdna},
//    repeat_masking, total_bases, num_sequences, fingerprint, warnings [{code, message}]}
// Errors: 404 UNKNOWN_GENOME, 422 AMBIGUOUS_ASSEMBLY, 503 MONGO_UNAVAILABLE (catalog), 500 INTERNAL (config).
// deps: {cfg, catalog | mongo, fs (promises API), openFasta(file) -> {getSequence}, now, resolve_ttl_ms,
//        cache (false disables the resolve cache), log}
async function resolve(systemName, deps) {
  deps = deps || {};
  if (!genomes.isValidSystemName(systemName)) throw genomes.unknownGenomeError(systemName);
  const catalog = await genomes.getCatalog(deps);
  const genome = catalog.bySystemName.get(systemName);
  if (!genome) throw genomes.unknownGenomeError(systemName);
  const cfg = getConfig(deps);
  if (!cfg || typeof cfg.fasta_root !== 'string' || !path.isAbsolute(cfg.fasta_root)) {
    throw new PrimerHttpError(500, 'INTERNAL', 'primers.fasta_root must be an absolute path', {});
  }
  if (deps.cache === false) return resolveUncached(genome, cfg, deps);

  const root = path.resolve(cfg.fasta_root);
  const key = resolveCacheKey(root, genome, cfg);
  const now = (deps.now || Date.now)();
  const ttl = deps.resolve_ttl_ms !== undefined ? deps.resolve_ttl_ms : RESOLVE_TTL_MS;
  const hit = resolveCache.get(key);
  if (hit && now - hit.at < ttl) return hit.promise;
  if (resolveCache.size >= RESOLVE_CACHE_MAX) pruneResolveCache(now, ttl);
  const entry = { at: now, promise: resolveUncached(genome, cfg, deps) };
  resolveCache.set(key, entry);
  entry.promise.catch(function (err) {
    // 4xx answers (AMBIGUOUS_ASSEMBLY) are kept for the TTL; transient failures are forgotten.
    if (!(isPrimerHttpError(err) && err.status < 500) && resolveCache.get(key) === entry) resolveCache.delete(key);
  });
  return entry.promise;
}

// resolveMany(names, deps) -> [{system_name, resolved, error}] in input order, at most 8 at a time.
// Per-genome failures are captured in `error` (resolved null); a catalog failure rejects the call.
async function resolveMany(names, deps) {
  deps = deps || {};
  const list = Array.from(names || []);
  if (list.length === 0) return [];
  const catalog = await genomes.getCatalog(deps);
  const shared = Object.assign({}, deps, { catalog: catalog });
  const limit = Math.max(1, Math.min(deps.concurrency || RESOLVE_CONCURRENCY, list.length));
  const out = new Array(list.length);
  let next = 0;
  async function worker() {
    while (next < list.length) {
      const i = next++;
      try {
        out[i] = { system_name: list[i], resolved: await resolve(list[i], shared), error: null };
      } catch (err) {
        out[i] = { system_name: list[i], resolved: null, error: err };
      }
    }
  }
  const workers = [];
  for (let k = 0; k < limit; k++) workers.push(worker());
  await Promise.all(workers);
  return out;
}

function clearCache() {
  resolveCache = new Map();
  maskCache = new Map();
}

module.exports = {
  RESOLVE_TTL_MS,
  RESOLVE_CONCURRENCY,
  REGION_MATCH_MIN_FRACTION,
  MASK_SAMPLE,
  REPEAT_MASKING,
  resolve,
  resolveMany,
  clearCache,
  _internal: {
    parseFai,
    faiIdentity,
    isSyntheticRegion,
    realRegionCount,
    regionMatches,
    scanCandidates,
    blastDbFiles,
    maskSampleWindows,
    detectSoftMask
  }
};
