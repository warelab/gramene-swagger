'use strict';

// check.normalize(body, deps) for POST /primers/check (spec §A.2.3 handler rules, plan Overrides).
//
// Returns {request, resolved, kind, warnings, estimate: {cpu_s, total}, dbs}
//   request   the canonical request that (with dbs) determines the job id:
//             {system_name, mode, [gene_id, transcript_id], checks (sorted, unique, incl. specificity),
//              genomes (sorted; [] unless pangenome), params (all 8 filled), pairs [{id, left, right, [expected]}]}
//             params.max_amplifying_mismatches must be < ignore_mismatches (400 INVALID_PARAMS); when it is
//             omitted, the default (3) is lowered to ignore_mismatches - 1 if needed.
//   resolved  {assemblies: {system_name: resolved assembly}, gene: {...} | null, species: {taxon_id, name} | null}
//             (stored in the job doc for the worker; never returned to clients)
//   kind      'pangenome' when checks include pangenome, else 'specificity'
//   warnings  [{code, message}]
//   estimate  {cpu_s, total} from ./cost.js (total = number of BLAST tasks)
//   dbs       {system_name: fingerprint} for the reference and every pan-genome genome
//
// deps: {cfg, catalog, mongo, genomes, assemblies, log, mongo_timeout_ms, ...anything getCatalog/resolve accept}

const { PrimerHttpError } = require('../errors');
const cost = require('./cost');

const MODES = Object.freeze(['gene', 'transcript', 'region', 'sequence']);
const CHECKS = Object.freeze(['specificity', 'pangenome']);
const DEFAULT_MODE = 'region';
const SYSTEM_NAME_RE = /^[a-z0-9_]+$/;
const SYSTEM_NAME_MAX = 128;
const PAIR_ID_RE = /^[A-Za-z0-9_.:-]+$/;
const PAIR_ID_MAX = 64;
const PRIMER_RE = /^[ACGTacgt]{15,36}$/;
const ID_MAX = 255;
const BODY_KEYS = Object.freeze(['system_name', 'mode', 'gene_id', 'transcript_id', 'checks', 'genomes', 'params', 'pairs']);
const PAIR_KEYS = Object.freeze(['id', 'left', 'right', 'expected']);
const EXPECTED_KEYS = Object.freeze(['region', 'start', 'end']);
const PARAM_RULES = Object.freeze({
  max_product_size: { type: 'integer', min: 50, max: 10000 },
  ignore_mismatches: { type: 'integer', min: 3, max: 6 },
  min_total_mismatches: { type: 'integer', min: 0, max: 6 },
  min_3p_mismatches: { type: 'integer', min: 1, max: 5 },
  three_prime_window: { type: 'integer', min: 3, max: 10 },
  include_unlikely: { type: 'boolean' },
  repeat_site_threshold: { type: 'integer', min: 1, max: 100 },
  max_amplifying_mismatches: { type: 'integer', min: 0, max: 5 }
});
const BUILTIN_DEFAULTS = Object.freeze({
  max_product_size: 4000, ignore_mismatches: 6, min_total_mismatches: 2, min_3p_mismatches: 2,
  three_prime_window: 5, include_unlikely: false, repeat_site_threshold: 5, max_amplifying_mismatches: 3
});
// Products the check can call at all (plan Overrides): expected sizes above this are refused.
const MAX_CHECKABLE_PRODUCT = 10000;
const PRODUCT_RAISE_FACTOR = 1.2;
const ORTHOLOG_KINDS = Object.freeze(['ortholog_one2one', 'ortholog_one2many', 'ortholog_many2many']);
const MONGO_TIMEOUT_MS = 10000;
const MONGO_RETRY_AFTER_S = 30;

function isPlainObject(v) {
  return v !== null && typeof v === 'object' && !Array.isArray(v) &&
    (Object.getPrototypeOf(v) === Object.prototype || Object.getPrototypeOf(v) === null);
}

function hasOwn(obj, key) {
  return Object.prototype.hasOwnProperty.call(obj, key);
}

function invalid(message, details) {
  return new PrimerHttpError(400, 'INVALID_REQUEST', message, details || {});
}

// ---- shape validation (mirrors the frozen PrimerCheckRequest; swagger normally rejects first) ----

function checkKeys(obj, allowed, where) {
  Object.keys(obj).forEach(function (k) {
    if (allowed.indexOf(k) < 0) throw invalid('unknown field ' + where + k, { field: where + k });
  });
}

function checkString(v, field, max, re) {
  if (typeof v !== 'string' || v.length === 0 || v.length > max || (re && !re.test(v))) {
    throw invalid(field + ' is not valid', { field: field });
  }
}

function validateShape(body, ccfg) {
  if (!isPlainObject(body)) throw invalid('the request body must be a JSON object');
  checkKeys(body, BODY_KEYS, '');
  checkString(body.system_name, 'system_name', SYSTEM_NAME_MAX, SYSTEM_NAME_RE);
  if (body.mode !== undefined && MODES.indexOf(body.mode) < 0) {
    throw invalid('mode must be one of ' + MODES.join(', '), { field: 'mode' });
  }
  if (body.gene_id !== undefined) checkString(body.gene_id, 'gene_id', ID_MAX);
  if (body.transcript_id !== undefined) checkString(body.transcript_id, 'transcript_id', ID_MAX);
  if (body.checks !== undefined) {
    if (!Array.isArray(body.checks) || body.checks.length < 1 || body.checks.length > 2 ||
        body.checks.some(function (c) { return CHECKS.indexOf(c) < 0; })) {
      throw invalid('checks must be 1-2 of ' + CHECKS.join(', '), { field: 'checks' });
    }
  }
  const maxGenomes = Number.isInteger(ccfg.max_genomes) ? ccfg.max_genomes : 150;
  if (body.genomes !== undefined) {
    if (!Array.isArray(body.genomes) || body.genomes.length > maxGenomes) {
      throw invalid('genomes must be an array of at most ' + maxGenomes + ' system names', { field: 'genomes', max: maxGenomes });
    }
    body.genomes.forEach(function (g, i) { checkString(g, 'genomes[' + i + ']', SYSTEM_NAME_MAX, SYSTEM_NAME_RE); });
  }
  if (body.params !== undefined) {
    if (!isPlainObject(body.params)) throw invalid('params must be an object', { field: 'params' });
    Object.keys(body.params).forEach(function (k) {
      const rule = PARAM_RULES[k];
      if (!rule) throw invalid('unknown field params.' + k, { field: 'params.' + k });
      const v = body.params[k];
      const ok = rule.type === 'boolean'
        ? typeof v === 'boolean'
        : Number.isInteger(v) && v >= rule.min && v <= rule.max;
      if (!ok) {
        throw invalid('params.' + k + (rule.type === 'boolean' ? ' must be a boolean' : ' must be an integer in ' + rule.min + '..' + rule.max),
          { field: 'params.' + k });
      }
    });
  }
  const maxPairs = Number.isInteger(ccfg.max_pairs) ? ccfg.max_pairs : 10;
  if (!Array.isArray(body.pairs) || body.pairs.length < 1 || body.pairs.length > maxPairs) {
    throw invalid('pairs must be an array of 1-' + maxPairs + ' primer pairs', { field: 'pairs', max: maxPairs });
  }
  body.pairs.forEach(function (p, i) {
    const where = 'pairs[' + i + ']';
    if (!isPlainObject(p)) throw invalid(where + ' must be an object', { field: where });
    checkKeys(p, PAIR_KEYS, where + '.');
    checkString(p.id, where + '.id', PAIR_ID_MAX, PAIR_ID_RE);
    checkString(p.left, where + '.left', 36, PRIMER_RE);
    checkString(p.right, where + '.right', 36, PRIMER_RE);
    if (p.expected !== undefined) {
      const e = p.expected;
      if (!isPlainObject(e)) throw invalid(where + '.expected must be an object', { field: where + '.expected' });
      checkKeys(e, EXPECTED_KEYS, where + '.expected.');
      checkString(e.region, where + '.expected.region', ID_MAX);
      ['start', 'end'].forEach(function (k) {
        if (!Number.isSafeInteger(e[k]) || e[k] < 1) {
          throw invalid(where + '.expected.' + k + ' must be a positive integer', { field: where + '.expected.' + k });
        }
      });
    }
  });
}

// ---- mongo helpers ---------------------------------------------------------------------------------

function mongoUnavailable() {
  return new PrimerHttpError(503, 'MONGO_UNAVAILABLE', 'the gene database is temporarily unavailable',
    { retry_after_s: MONGO_RETRY_AFTER_S });
}

function withTimeout(promise, ms) {
  let timer;
  const timeout = new Promise(function (resolve, reject) {
    timer = setTimeout(function () { reject(new Error('mongo query timed out after ' + ms + ' ms')); }, ms);
    if (timer.unref) timer.unref();
  });
  return Promise.race([promise, timeout]).finally(function () { clearTimeout(timer); });
}

async function mongoQuery(deps, fn) {
  const log = deps.log || console;
  const ms = deps.mongo_timeout_ms || MONGO_TIMEOUT_MS;
  try {
    return await withTimeout(Promise.resolve().then(fn), ms);
  } catch (err) {
    if (err && err.name === 'PrimerHttpError') throw err;
    try {
      (log.error || log.log).call(log, 'primers check: gene lookup failed: ' + (err && err.message ? err.message : String(err)));
    } catch (e) { /* ignore */ }
    throw mongoUnavailable();
  }
}

function uniqueStrings(values) {
  const seen = new Set();
  const out = [];
  (values || []).forEach(function (v) {
    if (typeof v === 'string' && v !== '' && !seen.has(v)) {
      seen.add(v);
      out.push(v);
    }
  });
  return out;
}

// Gene doc for gene/transcript modes -> resolved.gene
//   {id, system_name, taxon_id, location {region, start, end, strand, map}, transcripts [ids],
//    canonical_transcript, transcript_id, orthologs {system_name: [gene ids, sorted]}}
async function lookupGene(geneId, systemName, requestedTranscript, deps) {
  const mongo = deps.mongo || require('gramene-mongodb-config');
  const genes = await mongoQuery(deps, function () { return mongo.genes.mongoCollection(); });
  if (!genes) throw mongoUnavailable();
  const fields = {
    _id: 1, system_name: 1, taxon_id: 1, location: 1,
    'gene_structure.canonical_transcript': 1, 'gene_structure.transcripts.id': 1
  };
  ORTHOLOG_KINDS.forEach(function (k) { fields['homology.homologous_genes.' + k] = 1; });
  const docs = await mongoQuery(deps, function () { return genes.find({ _id: geneId }, { fields: fields }).toArray(); });
  const doc = Array.isArray(docs) ? docs.find(function (d) { return d && d._id === geneId; }) : null;
  if (!doc) throw new PrimerHttpError(404, 'UNKNOWN_GENE', 'unknown gene ' + geneId, { gene_id: geneId });
  if (doc.system_name !== systemName) {
    throw new PrimerHttpError(400, 'SYSTEM_NAME_MISMATCH',
      'gene ' + geneId + ' belongs to ' + doc.system_name + ', not ' + systemName,
      { gene_id: geneId, system_name: systemName, gene_system_name: doc.system_name === undefined ? null : doc.system_name });
  }
  const gs = doc.gene_structure || {};
  const transcripts = uniqueStrings((Array.isArray(gs.transcripts) ? gs.transcripts : []).map(function (t) { return t && t.id; }));
  const canonical = typeof gs.canonical_transcript === 'string' && gs.canonical_transcript ? gs.canonical_transcript : null;
  let transcriptId;
  if (requestedTranscript !== undefined) {
    if (transcripts.indexOf(requestedTranscript) < 0) {
      throw new PrimerHttpError(404, 'UNKNOWN_TRANSCRIPT', 'gene ' + geneId + ' has no transcript ' + requestedTranscript,
        { gene_id: geneId, transcript_id: requestedTranscript });
    }
    transcriptId = requestedTranscript;
  } else {
    transcriptId = canonical || transcripts[0] || null;
  }

  const hg = (doc.homology && doc.homology.homologous_genes) || {};
  const orthologIds = uniqueStrings([].concat.apply([], ORTHOLOG_KINDS.map(function (k) {
    return Array.isArray(hg[k]) ? hg[k] : [];
  })));
  const orthologs = {};
  if (orthologIds.length) {
    const odocs = await mongoQuery(deps, function () {
      return genes.find({ _id: { $in: orthologIds } }, { fields: { system_name: 1 } }).toArray();
    });
    (odocs || []).forEach(function (d) {
      if (!d || typeof d._id !== 'string' || typeof d.system_name !== 'string') return;
      if (!hasOwn(orthologs, d.system_name)) orthologs[d.system_name] = [];
      orthologs[d.system_name].push(d._id);
    });
    Object.keys(orthologs).forEach(function (k) { orthologs[k] = uniqueStrings(orthologs[k]).sort(); });
  }
  const loc = doc.location || {};
  return {
    id: doc._id,
    system_name: doc.system_name,
    taxon_id: doc.taxon_id === undefined ? null : doc.taxon_id,
    location: {
      region: loc.region === undefined ? null : String(loc.region),
      start: Number.isFinite(loc.start) ? loc.start : null,
      end: Number.isFinite(loc.end) ? loc.end : null,
      strand: loc.strand === -1 || loc.strand === '-1' ? -1 : 1,
      map: loc.map === undefined ? null : loc.map
    },
    transcripts: transcripts,
    canonical_transcript: canonical,
    transcript_id: transcriptId,
    orthologs: orthologs
  };
}

// ---- helpers -----------------------------------------------------------------------------------------

function plain(asm) {
  return JSON.parse(JSON.stringify(asm));
}

function normalizeChecks(checks) {
  const set = new Set(Array.isArray(checks) ? checks : []);
  set.add('specificity');
  return Array.from(set).sort();
}

function assemblyWarnings(systemName, asm) {
  return (asm.warnings || []).map(function (w) {
    return { code: w.code, message: systemName + ': ' + w.message };
  });
}

function notCheckable(names, reasons, message) {
  return new PrimerHttpError(400, 'GENOME_NOT_CHECKABLE', message, { genomes: names, reasons: reasons });
}

function outcomeReason(outcome, dbKind) {
  if (outcome.error) {
    return outcome.error.code === 'AMBIGUOUS_ASSEMBLY' ? 'ambiguous_assembly' : 'assembly_unavailable';
  }
  const asm = outcome.resolved;
  if (!asm || !asm.blastdb || !asm.blastdb[dbKind]) return dbKind === 'cdna' ? 'no_cdna_blastdb' : 'no_blastdb';
  return null;
}

// ---- normalize ---------------------------------------------------------------------------------------

async function normalize(body, deps) {
  deps = deps || {};
  const cfg = deps.cfg || require('../config').get();
  const ccfg = cfg.check || {};
  const genomesLib = deps.genomes || require('../genomes');
  const assembliesLib = deps.assemblies || require('../assemblies');
  const warnings = [];

  validateShape(body, ccfg);

  const systemName = body.system_name;
  const mode = body.mode === undefined ? DEFAULT_MODE : body.mode;
  const checks = normalizeChecks(body.checks);
  const pangenome = checks.indexOf('pangenome') >= 0;
  const kind = pangenome ? 'pangenome' : 'specificity';

  // Pairs: unique ids, uppercase primers, <= max_unique_primers distinct primers.
  const ids = new Set();
  const primers = new Set();
  const honoursExpected = mode === 'gene' || mode === 'region';
  let expectedDropped = 0;
  const pairs = body.pairs.map(function (p) {
    if (ids.has(p.id)) {
      throw new PrimerHttpError(400, 'DUPLICATE_PAIR_ID', 'pair id ' + p.id + ' is used more than once', { pair_id: p.id });
    }
    ids.add(p.id);
    const out = { id: p.id, left: p.left.toUpperCase(), right: p.right.toUpperCase() };
    primers.add(out.left);
    primers.add(out.right);
    if (p.expected !== undefined) {
      if (honoursExpected) out.expected = { region: p.expected.region, start: p.expected.start, end: p.expected.end };
      else expectedDropped++;
    }
    return out;
  });
  const maxPrimers = Number.isInteger(ccfg.max_unique_primers) ? ccfg.max_unique_primers : 20;
  if (primers.size > maxPrimers) {
    throw new PrimerHttpError(400, 'TOO_MANY_PRIMERS',
      'the pairs contain ' + primers.size + ' distinct primers; at most ' + maxPrimers + ' can be checked at once',
      { unique_primers: primers.size, max: maxPrimers });
  }
  if (expectedDropped > 0) {
    warnings.push({ code: 'EXPECTED_IGNORED', message: 'expected product locations are ignored in ' + mode + ' mode' });
  }

  // Expected products: start <= end, size <= 10 kb, auto-raise max_product_size.
  let largestExpected = 0;
  pairs.forEach(function (p) {
    if (!p.expected) return;
    if (p.expected.start > p.expected.end) {
      throw invalid('pair ' + p.id + ': expected.start must not be greater than expected.end', { pair_id: p.id });
    }
    const size = p.expected.end - p.expected.start + 1;
    if (size > MAX_CHECKABLE_PRODUCT) {
      throw new PrimerHttpError(400, 'PRODUCT_TOO_LONG_TO_CHECK',
        'pair ' + p.id + ': the expected product is ' + size + ' bp; products over ' + MAX_CHECKABLE_PRODUCT + ' bp cannot be checked',
        { pair_id: p.id, size: size, max: MAX_CHECKABLE_PRODUCT });
    }
    largestExpected = Math.max(largestExpected, size);
  });

  const defaults = Object.assign({}, BUILTIN_DEFAULTS, ccfg.defaults || {});
  const params = {};
  Object.keys(PARAM_RULES).forEach(function (k) {
    params[k] = body.params && body.params[k] !== undefined ? body.params[k] : defaults[k];
  });
  // A product amplifies only with <= max_amplifying_mismatches edits per primer, and products with a primer at
  // ignore_mismatches or more are dropped, so the cap must stay below ignore_mismatches. An omitted cap is
  // lowered to fit (ignore_mismatches 3 → cap 2); an explicit one that does not fit is refused.
  const capGiven = !!body.params && body.params.max_amplifying_mismatches !== undefined;
  if (!capGiven && Number.isInteger(params.max_amplifying_mismatches) && params.max_amplifying_mismatches >= params.ignore_mismatches) {
    params.max_amplifying_mismatches = params.ignore_mismatches - 1;
  }
  if (!(Number.isInteger(params.max_amplifying_mismatches) && params.max_amplifying_mismatches >= 0 &&
      params.max_amplifying_mismatches < params.ignore_mismatches)) {
    throw new PrimerHttpError(400, 'INVALID_PARAMS',
      'params.max_amplifying_mismatches (' + params.max_amplifying_mismatches + ') must be at least 0 and below ignore_mismatches (' +
      params.ignore_mismatches + ')',
      { field: 'params.max_amplifying_mismatches', max_amplifying_mismatches: params.max_amplifying_mismatches, ignore_mismatches: params.ignore_mismatches });
  }
  if (largestExpected > 0) {
    const needed = Math.ceil(PRODUCT_RAISE_FACTOR * largestExpected);
    const raised = Math.min(MAX_CHECKABLE_PRODUCT, Math.max(params.max_product_size, needed));
    if (raised > params.max_product_size) {
      warnings.push({
        code: 'MAX_PRODUCT_SIZE_RAISED',
        message: 'max_product_size raised from ' + params.max_product_size + ' to ' + raised +
          ' to cover the largest expected product (' + largestExpected + ' bp)'
      });
      params.max_product_size = raised;
    }
  }

  // Reference genome must be in the catalog (404) before any gene or filesystem work.
  const catalog = await genomesLib.getCatalog(deps);
  const shared = Object.assign({}, deps, { catalog: catalog, cfg: cfg });
  const refGenome = catalog.bySystemName.get(systemName);
  if (!refGenome) throw genomesLib.unknownGenomeError(systemName);

  // Gene / transcript modes: gene doc, system_name match, transcript default, orthologs.
  let gene = null;
  const request = { system_name: systemName, mode: mode };
  if (mode === 'gene' || mode === 'transcript') {
    if (body.gene_id === undefined) {
      throw invalid('gene_id is required in ' + mode + ' mode', { field: 'gene_id' });
    }
    gene = await lookupGene(body.gene_id, systemName, body.transcript_id, shared);
    if (mode === 'transcript' && !gene.transcript_id) {
      throw new PrimerHttpError(404, 'UNKNOWN_TRANSCRIPT', 'gene ' + body.gene_id + ' has no transcripts', { gene_id: body.gene_id });
    }
    request.gene_id = gene.id;
    if (gene.transcript_id) request.transcript_id = gene.transcript_id;
  }

  // Reference assembly: dna BLAST DB always, cdna DB in transcript mode.
  const ref = await assembliesLib.resolve(systemName, shared);
  if (!ref.blastdb || !ref.blastdb.dna) {
    throw new PrimerHttpError(422, 'NO_BLASTDB', systemName + ' has no genome BLAST database', { system_name: systemName, db: 'dna' });
  }
  if (mode === 'transcript' && !ref.blastdb.cdna) {
    throw new PrimerHttpError(422, 'NO_BLASTDB', systemName + ' has no cDNA BLAST database', { system_name: systemName, db: 'cdna' });
  }
  assemblyWarnings(systemName, ref).forEach(function (w) { warnings.push(w); });

  // Pan-genome genomes.
  const dbKind = mode === 'transcript' ? 'cdna' : 'dna';
  let genomeAsms = [];
  if (pangenome) {
    const members = await genomesLib.sameSpecies(systemName, shared);
    const speciesNames = new Set(members.map(function (g) { return g.system_name; }));
    if (body.genomes === undefined) {
      const names = members.map(function (g) { return g.system_name; }).filter(function (n) { return n !== systemName; });
      const outcomes = await assembliesLib.resolveMany(names, shared);
      genomeAsms = outcomes.filter(function (o) { return outcomeReason(o, dbKind) === null; })
        .map(function (o) { return o.resolved; });
      if (genomeAsms.length === 0) {
        throw notCheckable([], {}, 'no other genome of this species has a ' + (dbKind === 'cdna' ? 'cDNA ' : '') + 'BLAST database');
      }
    } else {
      const requested = uniqueStrings(body.genomes).filter(function (n) { return n !== systemName; });
      requested.forEach(function (n) {
        if (!catalog.bySystemName.has(n)) throw genomesLib.unknownGenomeError(n);
      });
      const reasons = {};
      const candidates = [];
      requested.forEach(function (n) {
        if (speciesNames.has(n)) candidates.push(n);
        else reasons[n] = 'other_species';
      });
      const outcomes = await assembliesLib.resolveMany(candidates, shared);
      outcomes.forEach(function (o) {
        const why = outcomeReason(o, dbKind);
        if (why) reasons[o.system_name] = why;
        else genomeAsms.push(o.resolved);
      });
      const bad = requested.filter(function (n) { return hasOwn(reasons, n); });
      if (bad.length) {
        throw notCheckable(bad, reasons, bad.length + ' of the requested genomes cannot be checked: ' + bad.join(', '));
      }
      if (genomeAsms.length === 0) throw notCheckable([], {}, 'no genomes to check besides ' + systemName);
    }
    genomeAsms.sort(function (a, b) { return a.system_name < b.system_name ? -1 : a.system_name > b.system_name ? 1 : 0; });
  } else if (Array.isArray(body.genomes) && body.genomes.length) {
    warnings.push({ code: 'GENOMES_IGNORED', message: 'genomes apply only when checks include pangenome' });
  }

  request.checks = checks;
  request.genomes = genomeAsms.map(function (a) { return a.system_name; });
  request.params = params;
  request.pairs = pairs;

  const estimate = cost.estimate({ unique_primers: primers.size, mode: mode, reference: ref, pangenome: genomeAsms, cfg: cfg });
  cost.assertWithinLimit(estimate, cfg);

  const assemblies = {};
  const dbs = {};
  [ref].concat(genomeAsms).forEach(function (a) {
    assemblies[a.system_name] = plain(a);
    dbs[a.system_name] = a.fingerprint === undefined ? null : a.fingerprint;
  });

  return {
    request: request,
    resolved: {
      assemblies: assemblies,
      gene: gene,
      species: refGenome.species ? { taxon_id: refGenome.species.taxon_id, name: refGenome.species.name } : null
    },
    kind: kind,
    warnings: warnings,
    estimate: { cpu_s: estimate.cpu_s, total: estimate.total },
    dbs: dbs
  };
}

module.exports = {
  normalize,
  validateShape,
  lookupGene,
  MODES,
  CHECKS,
  PARAM_RULES,
  MAX_CHECKABLE_PRODUCT,
  PRODUCT_RAISE_FACTOR,
  ORTHOLOG_KINDS
};
