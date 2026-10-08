'use strict';

// POST /primers/design (spec §A.2.1, §A.5, §A.6).
//
// design(body, deps): normalize and validate the request (presets merged before the client's params,
// cross-field checks), then under the design semaphore and ONE deadline shared by template, repeat
// mask and Primer3: build the template, mask it, run primer3_core, and annotate every primer with
// its genomic mapping and junction overlap. Returns {template, pairs, explain, settings, engine, warnings}.

const boulder = require('./boulder');
const coords = require('./coords');
const { PrimerHttpError } = require('./errors');

const MODES = Object.freeze(['gene', 'transcript', 'region', 'sequence']);
const MASK_MODES = Object.freeze(['n_mask', 'three_prime']);
const SYSTEM_NAME_RE = /^[a-z0-9_]+$/;
const MAX_ID_LENGTH = 255;
const MAX_SYSTEM_NAME_LENGTH = 128;
const MAX_SEQUENCE_CHARS = 60000;
const MAX_EXCLUDED = 50;
const MAX_JUNCTIONS = 200; // Primer3 PR_MAX_INTERVAL_ARRAY
const DEFAULT_MIN_5_PRIME_OVERLAP = 7;
const DEFAULT_MIN_3_PRIME_OVERLAP = 4;

const TOP_LEVEL_FIELDS = Object.freeze(['mode', 'gene_id', 'transcript_id', 'system_name', 'region', 'sequence',
  'flank_up', 'flank_down', 'target', 'included', 'excluded', 'avoid_repeats', 'repeat_mask_mode',
  'junction_spanning', 'template_only', 'params']);
const REGION_FIELDS = Object.freeze(['region', 'start', 'end', 'strand']);

// Closed set of client params (PrimerDesignParams). int: integer; ranges: product_size_ranges.
const PARAM_SPECS = deepFreeze({
  opt_size: { int: true, min: 15, max: 36 },
  min_size: { int: true, min: 15, max: 36 },
  max_size: { int: true, min: 15, max: 36 },
  opt_tm: { min: 30, max: 90 },
  min_tm: { min: 30, max: 90 },
  max_tm: { min: 30, max: 90 },
  opt_gc: { min: 0, max: 100 },
  min_gc: { min: 0, max: 100 },
  max_gc: { min: 0, max: 100 },
  max_tm_diff: { min: 0, max: 30 },
  max_poly_x: { int: true, min: 0, max: 10 },
  gc_clamp: { int: true, min: 0, max: 5 },
  max_end_stability: { min: 0, max: 100 },
  max_ns: { int: true, min: 0, max: 5 },
  salt_monovalent: { min: 0, max: 1000 },
  salt_divalent: { min: 0, max: 100 },
  dntp_conc: { min: 0, max: 100 },
  dna_conc: { min: 0, max: 10000 },
  num_return: { int: true, min: 1, max: 20 },
  min_3_prime_overlap_of_junction: { int: true, min: 1, max: 20 },
  min_5_prime_overlap_of_junction: { int: true, min: 1, max: 20 },
  product_size_ranges: { ranges: true, minItems: 1, maxItems: 10, min: 20, max: 50000 }
});

const PRESETS = deepFreeze({
  pcr: {
    opt_size: 20, min_size: 18, max_size: 25,
    opt_tm: 60, min_tm: 57, max_tm: 63,
    min_gc: 30, max_gc: 70,
    max_tm_diff: 3, max_poly_x: 4,
    product_size_ranges: [[100, 1000]],
    num_return: 5
  },
  qpcr: {
    opt_size: 20, min_size: 18, max_size: 24,
    opt_tm: 60, min_tm: 58, max_tm: 62,
    min_gc: 35, max_gc: 65,
    max_tm_diff: 2, max_poly_x: 4,
    product_size_ranges: [[70, 150]],
    num_return: 5,
    min_3_prime_overlap_of_junction: DEFAULT_MIN_3_PRIME_OVERLAP,
    min_5_prime_overlap_of_junction: DEFAULT_MIN_5_PRIME_OVERLAP
  }
});

// Server-owned tags, in record order (§A.5 "Always sent").
const ALWAYS_TAGS = Object.freeze([
  ['PRIMER_TASK', 'generic'],
  ['PRIMER_PICK_LEFT_PRIMER', 1],
  ['PRIMER_PICK_RIGHT_PRIMER', 1],
  ['PRIMER_PICK_INTERNAL_OLIGO', 0],
  ['PRIMER_FIRST_BASE_INDEX', 1],
  ['PRIMER_EXPLAIN_FLAG', 1],
  ['PRIMER_LIBERAL_BASE', 1],
  ['PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT', 1],
  ['PRIMER_THERMODYNAMIC_TEMPLATE_ALIGNMENT', 0],
  ['PRIMER_PRODUCT_MIN_TM', 0],
  ['PRIMER_PRODUCT_MAX_TM', 150],
  ['P3_FILE_FLAG', 0]
]);

// Client params -> Primer3 tags, in record order.
const PARAM_TAGS = Object.freeze([
  ['opt_size', 'PRIMER_OPT_SIZE'],
  ['min_size', 'PRIMER_MIN_SIZE'],
  ['max_size', 'PRIMER_MAX_SIZE'],
  ['opt_tm', 'PRIMER_OPT_TM'],
  ['min_tm', 'PRIMER_MIN_TM'],
  ['max_tm', 'PRIMER_MAX_TM'],
  ['opt_gc', 'PRIMER_OPT_GC_PERCENT'],
  ['min_gc', 'PRIMER_MIN_GC'],
  ['max_gc', 'PRIMER_MAX_GC'],
  ['max_tm_diff', 'PRIMER_PAIR_MAX_DIFF_TM'],
  ['max_poly_x', 'PRIMER_MAX_POLY_X'],
  ['gc_clamp', 'PRIMER_GC_CLAMP'],
  ['max_end_stability', 'PRIMER_MAX_END_STABILITY'],
  ['max_ns', 'PRIMER_MAX_NS_ACCEPTED'],
  ['salt_monovalent', 'PRIMER_SALT_MONOVALENT'],
  ['salt_divalent', 'PRIMER_SALT_DIVALENT'],
  ['dntp_conc', 'PRIMER_DNTP_CONC'],
  ['dna_conc', 'PRIMER_DNA_CONC']
]);

const MESSAGES = Object.freeze({
  NO_PAIRS: 'Primer3 found no acceptable primer pairs; see explain for the reasons candidates were rejected',
  SINGLE_EXON_TRANSCRIPT: 'the transcript has a single exon, so primers were designed without the junction-spanning constraint',
  JUNCTIONS_TRUNCATED: 'the transcript has more than ' + MAX_JUNCTIONS + ' exon junctions; only ' + MAX_JUNCTIONS +
    ' were used for the junction-spanning constraint'
});

function deepFreeze(o) {
  if (o !== null && typeof o === 'object' && !Object.isFrozen(o)) {
    Object.freeze(o);
    Object.keys(o).forEach(function (k) { deepFreeze(o[k]); });
  }
  return o;
}

function logSafe(log, level, msg) {
  try {
    const l = log || console;
    (l[level] || l.error || l.log).call(l, msg);
  } catch (e) { /* logging must never throw */ }
}

function warning(code, message) {
  return { code: code, message: message || MESSAGES[code] };
}

function invalidRequest(message, details) {
  return new PrimerHttpError(400, 'INVALID_REQUEST', message, details || {});
}

function invalidParams(message, details) {
  return new PrimerHttpError(400, 'INVALID_PARAMS', message, details || {});
}

function isPosInt(v) {
  return Number.isSafeInteger(v) && v >= 1;
}

function absent(v) {
  return v === undefined || v === null;
}

// ---- request normalization -------------------------------------------------------------------------

function optString(body, field, maxLength) {
  const v = body[field];
  if (absent(v)) return undefined;
  if (typeof v !== 'string' || v === '' || v.length > maxLength) {
    throw invalidRequest(field + ' must be a non-empty string of at most ' + maxLength + ' characters', { field: field });
  }
  return v;
}

function optBool(body, field) {
  const v = body[field];
  if (absent(v)) return undefined;
  if (typeof v !== 'boolean') throw invalidRequest(field + ' must be true or false', { field: field });
  return v;
}

function optInt(body, field, min, max) {
  const v = body[field];
  if (absent(v)) return undefined;
  if (!Number.isSafeInteger(v) || v < min || v > max) {
    throw invalidRequest(field + ' must be an integer from ' + min + ' to ' + max, { field: field, min: min, max: max });
  }
  return v;
}

function parseInterval(v, field) {
  if (!Array.isArray(v) || v.length !== 2 || !isPosInt(v[0]) || !isPosInt(v[1])) {
    throw invalidRequest(field + ' must be [start, length] with positive integers (1-based template coordinates)', { field: field });
  }
  return [v[0], v[1]];
}

function optRegion(body) {
  const r = body.region;
  if (absent(r)) return undefined;
  if (typeof r !== 'object' || Array.isArray(r)) throw invalidRequest('region must be an object {region, start, end, strand}', { field: 'region' });
  Object.keys(r).forEach(function (k) {
    if (REGION_FIELDS.indexOf(k) < 0) throw invalidRequest('unknown field region.' + k.slice(0, 64), { field: 'region.' + k.slice(0, 64) });
  });
  if (typeof r.region !== 'string' || r.region === '' || r.region.length > MAX_ID_LENGTH) {
    throw invalidRequest('region.region must be a non-empty string of at most ' + MAX_ID_LENGTH + ' characters', { field: 'region.region' });
  }
  ['start', 'end'].forEach(function (k) {
    if (!isPosInt(r[k])) throw invalidRequest('region.' + k + ' must be an integer >= 1', { field: 'region.' + k });
  });
  const strand = absent(r.strand) ? 1 : r.strand;
  if (strand !== 1 && strand !== -1) throw invalidRequest('region.strand must be 1 or -1', { field: 'region.strand' });
  return { region: r.region, start: r.start, end: r.end, strand: strand };
}

function checkParamValue(key, v) {
  const spec = PARAM_SPECS[key];
  if (spec.ranges) {
    if (!Array.isArray(v) || v.length < spec.minItems || v.length > spec.maxItems) {
      throw invalidParams('product_size_ranges must be a list of 1 to ' + spec.maxItems + ' [min, max] ranges', { param: key });
    }
    return v.map(function (r, i) {
      if (!Array.isArray(r) || r.length !== 2 || !r.every(function (x) { return Number.isSafeInteger(x) && x >= spec.min && x <= spec.max; })) {
        throw invalidParams('product_size_ranges[' + i + '] must be [min, max] with integers from ' + spec.min + ' to ' + spec.max,
          { param: key, index: i });
      }
      return [r[0], r[1]];
    });
  }
  const ok = typeof v === 'number' && Number.isFinite(v) && (!spec.int || Number.isInteger(v)) && v >= spec.min && v <= spec.max;
  if (!ok) {
    throw invalidParams(key + ' must be ' + (spec.int ? 'an integer' : 'a number') + ' from ' + spec.min + ' to ' + spec.max,
      { param: key, min: spec.min, max: spec.max });
  }
  return v;
}

function parseParams(params) {
  if (absent(params)) return {};
  if (typeof params !== 'object' || Array.isArray(params)) throw invalidParams('params must be an object', { param: null });
  const out = {};
  Object.keys(params).forEach(function (k) {
    if (!Object.prototype.hasOwnProperty.call(PARAM_SPECS, k)) {
      throw invalidParams('unknown parameter ' + JSON.stringify(k.slice(0, 64)), { param: k.slice(0, 64) });
    }
    if (params[k] === undefined) return;
    out[k] = checkParamValue(k, params[k]);
  });
  return out;
}

function presetParams(preset) {
  return JSON.parse(JSON.stringify(PRESETS[preset]));
}

function orderCheck(p, name) {
  const lo = p['min_' + name];
  const opt = p['opt_' + name];
  const hi = p['max_' + name];
  if (lo !== undefined && hi !== undefined && lo > hi) {
    throw invalidParams('min_' + name + ' (' + lo + ') must not exceed max_' + name + ' (' + hi + ')', { param: 'min_' + name });
  }
  if (opt !== undefined && lo !== undefined && opt < lo) {
    throw invalidParams('opt_' + name + ' (' + opt + ') must not be below min_' + name + ' (' + lo + ')', { param: 'opt_' + name });
  }
  if (opt !== undefined && hi !== undefined && opt > hi) {
    throw invalidParams('opt_' + name + ' (' + opt + ') must not exceed max_' + name + ' (' + hi + ')', { param: 'opt_' + name });
  }
}

// Cross-field checks that need no template (§A.5). ctx: {junctionSpanning}
function checkParams(p, ctx) {
  orderCheck(p, 'size');
  orderCheck(p, 'tm');
  orderCheck(p, 'gc');
  const ranges = p.product_size_ranges;
  ranges.forEach(function (r, i) {
    if (!(r[0] < r[1])) throw invalidParams('product_size_ranges[' + i + '] must have min < max', { param: 'product_size_ranges', index: i });
  });
  const smallest = Math.min.apply(null, ranges.map(function (r) { return r[0]; }));
  if (p.max_size > smallest) {
    throw invalidParams('max_size (' + p.max_size + ') must not exceed the smallest product size (' + smallest + ')', { param: 'max_size' });
  }
  if (p.gc_clamp !== undefined && p.gc_clamp > p.min_size) {
    throw invalidParams('gc_clamp must not exceed min_size', { param: 'gc_clamp' });
  }
  if (ctx && ctx.junctionSpanning) {
    // Primer3 aborts the whole record otherwise (libprimer3.cc: overlap > PRIMER_MAX_SIZE / 2, integer division).
    const half = Math.floor(p.max_size / 2);
    [['min_5_prime_overlap_of_junction', DEFAULT_MIN_5_PRIME_OVERLAP], ['min_3_prime_overlap_of_junction', DEFAULT_MIN_3_PRIME_OVERLAP]]
      .forEach(function (kd) {
        const v = p[kd[0]] === undefined ? kd[1] : p[kd[0]];
        if (v > half) {
          throw invalidParams(kd[0] + ' (' + v + ') must be at most floor(max_size / 2) = ' + half + ' for junction-spanning primers',
            { param: kd[0], max: half });
        }
      });
  }
}

// normalize(body, cfg) -> req. Throws 400 INVALID_REQUEST / INVALID_PARAMS / REGION_OUT_OF_BOUNDS /
// TEMPLATE_TOO_LONG before any I/O.
//   req: {mode, gene_id, transcript_id, system_name, region {region,start,end,strand} | null, sequence,
//         flank_up, flank_down, target, included, excluded[], avoid_repeats, repeat_mask_mode,
//         junction_spanning, template_only, preset, params (effective)}
function normalize(body, cfg) {
  cfg = cfg || require('./config').get();
  if (body === null || typeof body !== 'object' || Array.isArray(body)) throw invalidRequest('the request body must be a JSON object');
  Object.keys(body).forEach(function (k) {
    if (TOP_LEVEL_FIELDS.indexOf(k) < 0) throw invalidRequest('unknown field ' + JSON.stringify(k.slice(0, 64)), { field: k.slice(0, 64) });
  });
  const mode = body.mode;
  if (MODES.indexOf(mode) < 0) throw invalidRequest('mode must be one of ' + MODES.join(', '), { field: 'mode' });

  const geneId = optString(body, 'gene_id', MAX_ID_LENGTH);
  const transcriptId = optString(body, 'transcript_id', MAX_ID_LENGTH);
  const systemName = optString(body, 'system_name', MAX_SYSTEM_NAME_LENGTH);
  if (systemName !== undefined && !SYSTEM_NAME_RE.test(systemName)) {
    throw invalidRequest('system_name must match ^[a-z0-9_]+$', { field: 'system_name' });
  }
  const region = optRegion(body);
  let sequence;
  if (!absent(body.sequence)) {
    if (typeof body.sequence !== 'string' || body.sequence.length > MAX_SEQUENCE_CHARS) {
      throw invalidRequest('sequence must be a string of at most ' + MAX_SEQUENCE_CHARS + ' characters', { field: 'sequence' });
    }
    sequence = body.sequence;
  }
  const flankUp = optInt(body, 'flank_up', 0, cfg.design.max_flank);
  const flankDown = optInt(body, 'flank_down', 0, cfg.design.max_flank);
  const target = absent(body.target) ? null : parseInterval(body.target, 'target');
  const included = absent(body.included) ? null : parseInterval(body.included, 'included');
  let excluded = [];
  if (!absent(body.excluded)) {
    if (!Array.isArray(body.excluded) || body.excluded.length > MAX_EXCLUDED) {
      throw invalidRequest('excluded must be a list of at most ' + MAX_EXCLUDED + ' [start, length] intervals', { field: 'excluded' });
    }
    excluded = body.excluded.map(function (iv, i) { return parseInterval(iv, 'excluded[' + i + ']'); });
  }
  const avoidRepeats = optBool(body, 'avoid_repeats') === true;
  const junctionSpanning = optBool(body, 'junction_spanning');
  const templateOnly = optBool(body, 'template_only') === true;
  const maskMode = absent(body.repeat_mask_mode) ? undefined : body.repeat_mask_mode;
  if (maskMode !== undefined && MASK_MODES.indexOf(maskMode) < 0) {
    throw invalidRequest('repeat_mask_mode must be n_mask or three_prime', { field: 'repeat_mask_mode' });
  }
  const userParams = parseParams(body.params);

  const geneMode = mode === 'gene' || mode === 'transcript';
  if (geneMode && geneId === undefined) throw invalidRequest('gene_id is required in ' + mode + ' mode', { field: 'gene_id' });
  if (mode === 'region' && systemName === undefined) throw invalidRequest('system_name is required in region mode', { field: 'system_name' });
  if (mode === 'region' && region === undefined) throw invalidRequest('region is required in region mode', { field: 'region' });
  if (mode === 'sequence' && sequence === undefined) throw invalidRequest('sequence is required in sequence mode', { field: 'sequence' });
  if (mode === 'region') {
    if (region.end < region.start) {
      throw new PrimerHttpError(400, 'REGION_OUT_OF_BOUNDS', 'region start must not exceed end',
        { region: region.region, start: region.start, end: region.end });
    }
    const length = region.end - region.start + 1;
    if (length > cfg.design.max_template_length) {
      throw new PrimerHttpError(400, 'TEMPLATE_TOO_LONG', 'the template would be ' + length + ' bp; the limit is ' +
        cfg.design.max_template_length + ' bp', { length: length, max: cfg.design.max_template_length });
    }
  }

  const preset = mode === 'transcript' ? 'qpcr' : 'pcr';
  const params = Object.assign(presetParams(preset), userParams);
  const req = {
    mode: mode,
    gene_id: geneMode ? geneId : null,
    transcript_id: geneMode && transcriptId !== undefined ? transcriptId : null,
    system_name: systemName === undefined ? null : systemName,
    region: mode === 'region' ? region : null,
    sequence: mode === 'sequence' ? sequence : null,
    flank_up: mode === 'gene' ? flankUp || 0 : 0,
    flank_down: mode === 'gene' ? flankDown || 0 : 0,
    target: target,
    included: included,
    excluded: excluded,
    avoid_repeats: avoidRepeats,
    repeat_mask_mode: avoidRepeats ? maskMode || 'n_mask' : null,
    junction_spanning: mode === 'transcript' ? junctionSpanning !== false : false,
    template_only: templateOnly,
    preset: preset,
    params: params
  };
  // n_mask turns masked bases into N, so no N may be accepted in a primer (forced, echoed in settings).
  if (req.repeat_mask_mode === 'n_mask') params.max_ns = 0;
  checkParams(params, { junctionSpanning: req.junction_spanning });
  return req;
}

// Checks that need the template: intervals inside it and a product range that can fit. template_only runs no
// Primer3, so its product ranges are not checked against the template length (every other check still applies).
function validateAgainstTemplate(req, template) {
  const len = template.length;
  const check = function (field, iv) {
    if (iv[0] > len || iv[0] + iv[1] - 1 > len) {
      throw new PrimerHttpError(400, 'INTERVAL_OUT_OF_BOUNDS',
        field + ' [' + iv[0] + ', ' + iv[1] + '] extends beyond the template (' + len + ' bp)',
        { field: field, interval: [iv[0], iv[1]], template_length: len });
    }
  };
  if (req.target) check('target', req.target);
  if (req.included) check('included', req.included);
  (req.excluded || []).forEach(function (iv, i) { check('excluded[' + i + ']', iv); });
  if (req.template_only) return;
  const ranges = req.params.product_size_ranges;
  if (!ranges.some(function (r) { return r[0] <= len; })) {
    throw invalidParams('every product_size_ranges entry starts beyond the template length (' + len + ' bp)',
      { param: 'product_size_ranges', template_length: len });
  }
}

// Junctions sent in SEQUENCE_OVERLAP_JUNCTION_LIST (null = no constraint). Pushes warnings.
function junctionList(req, template, warnings) {
  if (req.mode !== 'transcript' || !req.junction_spanning) return null;
  const all = template.junctions || [];
  if (all.length === 0) {
    warnings.push(warning('SINGLE_EXON_TRANSCRIPT'));
    return null;
  }
  if (all.length <= MAX_JUNCTIONS) return all.slice();
  let keep = all;
  if (req.included) {
    const s = req.included[0];
    const e = s + req.included[1] - 1;
    const inside = all.filter(function (j) { return j >= s && j < e; });
    if (inside.length > 0) keep = inside;
  }
  warnings.push(warning('JUNCTIONS_TRUNCATED'));
  return keep.slice(0, MAX_JUNCTIONS);
}

function formatInterval(iv) {
  return iv[0] + ',' + iv[1];
}

// Pure Boulder-IO tag builder (§A.5). Returns an insertion-ordered tags object for boulder.serialize.
//   template: {id, seq} - seq is the exact Primer3 template (masked, uppercase unless three_prime)
//   settings: {params (effective), target, included, excluded, junctions (list to send | null),
//              avoid_repeats, repeat_mask_mode}
function buildRecord(template, settings) {
  const p = settings.params;
  const tags = {};
  tags.SEQUENCE_ID = boulder.sanitizeId(template.id);
  tags.SEQUENCE_TEMPLATE = template.seq;
  if (settings.target) tags.SEQUENCE_TARGET = formatInterval(settings.target);
  if (settings.included) tags.SEQUENCE_INCLUDED_REGION = formatInterval(settings.included);
  if (settings.excluded && settings.excluded.length > 0) tags.SEQUENCE_EXCLUDED_REGION = settings.excluded.map(formatInterval).join(' ');
  ALWAYS_TAGS.forEach(function (kv) { tags[kv[0]] = kv[1]; });
  tags.PRIMER_NUM_RETURN = p.num_return;
  tags.PRIMER_PRODUCT_SIZE_RANGE = p.product_size_ranges.map(function (r) { return r[0] + '-' + r[1]; }).join(' ');
  const masking = settings.avoid_repeats ? settings.repeat_mask_mode || 'n_mask' : null;
  PARAM_TAGS.forEach(function (kt) {
    let v = p[kt[0]];
    if (kt[0] === 'max_ns' && masking === 'n_mask') v = 0;
    if (!absent(v)) tags[kt[1]] = v;
  });
  if (masking === 'n_mask' && tags.PRIMER_MAX_NS_ACCEPTED === undefined) tags.PRIMER_MAX_NS_ACCEPTED = 0;
  if (settings.junctions && settings.junctions.length > 0) {
    const min3 = absent(p.min_3_prime_overlap_of_junction) ? DEFAULT_MIN_3_PRIME_OVERLAP : p.min_3_prime_overlap_of_junction;
    const min5 = absent(p.min_5_prime_overlap_of_junction) ? DEFAULT_MIN_5_PRIME_OVERLAP : p.min_5_prime_overlap_of_junction;
    tags.SEQUENCE_OVERLAP_JUNCTION_LIST = settings.junctions.join(' ');
    tags.PRIMER_MIN_3_PRIME_OVERLAP_OF_JUNCTION = min3;
    tags.PRIMER_MIN_5_PRIME_OVERLAP_OF_JUNCTION = min5;
    tags.PRIMER_INTERNAL_MIN_3_PRIME_OVERLAP_OF_JUNCTION = min3;
    tags.PRIMER_INTERNAL_MIN_5_PRIME_OVERLAP_OF_JUNCTION = min5;
  }
  if (masking === 'three_prime') tags.PRIMER_LOWERCASE_MASKING = 1;
  return tags;
}

// Adds junction overlap and genomic mappings to extracted pairs (in place).
function decoratePairs(pairs, template, params) {
  const mapper = template.mapper;
  const junctions = template.junctions || [];
  const jopts = {
    min5: absent(params.min_5_prime_overlap_of_junction) ? DEFAULT_MIN_5_PRIME_OVERLAP : params.min_5_prime_overlap_of_junction,
    min3: absent(params.min_3_prime_overlap_of_junction) ? DEFAULT_MIN_3_PRIME_OVERLAP : params.min_3_prime_overlap_of_junction
  };
  pairs.forEach(function (pair) {
    ['left', 'right'].forEach(function (side) {
      const o = pair[side];
      o.seq = o.seq.toUpperCase();
      o.junction = junctions.length > 0 ? coords.junctionOverlap(side, o.start, o.end, junctions, jopts) : null;
      o.genomic = mapper ? mapper.primer(side, o.start, o.end) : null;
    });
    pair.product.genomic = mapper ? mapper.product(pair.product.start, pair.product.end) : null;
    if (mapper && mapper.kind === 'spliced') pair.product.genomic_size = mapper.genomicSize(pair.product.start, pair.product.end);
  });
  return pairs;
}

// Primer3 2.6.1 defaults of the salt terms in the product Tm (libprimer3.cc pr_set_default_global_args_2).
const PRIMER3_SALT_DEFAULTS = Object.freeze({ salt_monovalent: 50, salt_divalent: 1.5, dntp_conc: 0.6 });

// Primer3's product Tm, oligotm.c long_seq_tm as printed in PRIMER_PAIR_i_PRODUCT_TM (no DMSO or formamide):
//   81.5 + 16.6 log10(Na / 1000) + 41 GC / L - 600 / L,  Na = monovalent + divalent_to_monovalent(divalent, dNTP)
//   divalent_to_monovalent = 120 sqrt(divalent - dNTP), or 0 when divalent <= dNTP.
// Only uppercase G and C count as GC. params: effective design params (Primer3 defaults for absent salts).
// Returns null when the value is not finite.
function longSeqTm(seq, params) {
  const p = params || {};
  const pick = function (k) { return absent(p[k]) ? PRIMER3_SALT_DEFAULTS[k] : p[k]; };
  const mono = pick('salt_monovalent');
  const div = pick('salt_divalent');
  const dntp = pick('dntp_conc');
  const na = mono + (div > dntp ? 120 * Math.sqrt(div - dntp) : 0);
  const len = seq.length;
  let gc = 0;
  for (let i = 0; i < len; i++) {
    const ch = seq.charCodeAt(i);
    if (ch === 71 || ch === 67) gc++; // 'G', 'C'
  }
  const tm = 81.5 + 16.6 * Math.log10(na / 1000) + 41 * gc / len - 600 / len;
  return Number.isFinite(tm) ? tm : null;
}

// n_mask gives Primer3 the masked bases as N, which long_seq_tm counts as A/T, so the product Tm of a product
// that covers masked bases comes out too low. Recompute those on the unmasked template (the value three_prime
// mode reports), rounded like Primer3's %.4f. Only Primer3-trusted values (not null) are replaced. In place.
function unmaskProductTm(pairs, seq, runs, params) {
  pairs.forEach(function (pair) {
    if (absent(pair.product_tm) || !coords.overlapsRuns(pair.product.start, pair.product.end, runs)) return;
    const tm = longSeqTm(String(seq).slice(pair.product.start - 1, pair.product.end).toUpperCase(), params);
    pair.product_tm = tm === null ? null : Math.round(tm * 10000) / 10000;
  });
  return pairs;
}

function publicTemplate(t, mask) {
  return {
    mode: t.mode,
    system_name: t.system_name,
    gene_id: t.gene_id,
    transcript_id: t.transcript_id,
    region: t.region,
    start: t.start,
    end: t.end,
    strand: t.strand,
    length: t.length,
    seq: t.seq.toUpperCase(),
    masked: mask.masked,
    mask_source: mask.mask_source,
    mask: mask.mask,
    masked_fraction: mask.masked_fraction,
    features: t.features
  };
}

// ---- deadline --------------------------------------------------------------------------------------

function deadlineError(ms) {
  return new PrimerHttpError(504, 'DEADLINE_EXCEEDED', 'the design did not finish within ' + ms + ' ms', { deadline_ms: ms });
}

// One deadline for template, mask and Primer3. -> {signal, deadlineAt, remaining(), within(promise), dispose()}
function createDeadline(ms, parentSignal, now) {
  now = now || Date.now;
  const controller = new AbortController();
  const deadlineAt = now() + ms;
  const timer = setTimeout(function () { controller.abort(deadlineError(ms)); }, Math.max(0, ms));
  let onParent = null;
  const parentReason = function () {
    const r = parentSignal.reason;
    return r instanceof PrimerHttpError ? r : new PrimerHttpError(504, 'DEADLINE_EXCEEDED', 'the request was aborted', {});
  };
  if (parentSignal) {
    if (parentSignal.aborted) controller.abort(parentReason());
    else {
      onParent = function () { controller.abort(parentReason()); };
      parentSignal.addEventListener('abort', onParent, { once: true });
    }
  }
  const aborted = new Promise(function (resolve, reject) {
    if (controller.signal.aborted) return reject(controller.signal.reason);
    controller.signal.addEventListener('abort', function () { reject(controller.signal.reason); }, { once: true });
  });
  aborted.catch(function () { /* observed through within() */ });
  return {
    signal: controller.signal,
    deadlineAt: deadlineAt,
    remaining: function () { return deadlineAt - now(); },
    within: function (promise) { return Promise.race([promise, aborted]); },
    dispose: function () {
      clearTimeout(timer);
      if (onParent) parentSignal.removeEventListener('abort', onParent);
    }
  };
}

function versionOrNull(primer3, log) {
  return Promise.resolve().then(function () { return primer3.version(); }).catch(function (err) {
    logSafe(log, 'warn', 'primers: primer3 version unavailable: ' + (err && err.code));
    return null;
  });
}

// ---- orchestration ---------------------------------------------------------------------------------

async function runDesign(req, cfg, deps, dl) {
  const log = deps.log || console;
  const primer3 = deps.primer3 || require('./primer3');
  const templates = deps.template || require('./template');
  const tdeps = Object.assign({}, deps, { cfg: cfg, signal: dl.signal, log: log });
  const warnings = [];

  const template = await dl.within(templates.buildTemplate(req, tdeps));
  template.warnings.forEach(function (w) { warnings.push(w); });
  validateAgainstTemplate(req, template);
  const junctions = junctionList(req, template, warnings);

  let mask = { masked: false, mask_source: null, mask: [], masked_fraction: 0 };
  let primer3Seq = template.seq.toUpperCase();
  if (req.avoid_repeats) {
    const masker = deps.repeatMask || require('./repeat_mask');
    const m = await dl.within(masker.repeatMask(template, { mode: req.repeat_mask_mode, signal: dl.signal, deadline: dl.deadlineAt }, tdeps));
    mask = { masked: m.masked, mask_source: m.mask_source, mask: m.mask, masked_fraction: m.masked_fraction };
    primer3Seq = m.seq;
    m.warnings.forEach(function (w) { warnings.push(w); });
  }

  const response = {
    template: publicTemplate(template, mask),
    pairs: [],
    explain: null,
    settings: {
      preset: req.preset,
      junction_spanning: req.junction_spanning,
      avoid_repeats: req.avoid_repeats,
      repeat_mask_mode: req.repeat_mask_mode,
      params: req.params
    },
    engine: { primer3: null },
    warnings: warnings
  };

  if (req.template_only) {
    response.engine.primer3 = await dl.within(versionOrNull(primer3, log));
    return response;
  }

  const tags = buildRecord({ id: template.id, seq: primer3Seq }, {
    params: req.params, target: req.target, included: req.included, excluded: req.excluded,
    junctions: junctions, avoid_repeats: req.avoid_repeats, repeat_mask_mode: req.repeat_mask_mode
  });
  const remaining = Math.floor(dl.remaining());
  if (remaining <= 0) throw deadlineError(cfg.design.deadline_ms);
  const timeoutMs = Math.min(cfg.design.primer3_timeout_ms, remaining);
  const results = await dl.within(Promise.all([
    primer3.run(tags, { timeoutMs: timeoutMs, signal: dl.signal, log: log, tmpDir: cfg.tmp_dir }),
    versionOrNull(primer3, log)
  ]));
  const result = results[0];
  response.engine.primer3 = results[1];
  if (result.error) {
    throw new PrimerHttpError(400, 'PRIMER3_INPUT_ERROR', 'Primer3 rejected the input: ' + result.error, { primer3_error: result.error });
  }
  if (result.warning) warnings.push(warning('PRIMER3_WARNING', result.warning));
  response.pairs = decoratePairs(boulder.extractPairs(result.tags), template, req.params);
  if (req.repeat_mask_mode === 'n_mask' && mask.mask.length > 0) unmaskProductTm(response.pairs, template.seq, mask.mask, req.params);
  response.explain = boulder.extractExplain(result.tags);
  if (response.pairs.length === 0) warnings.push(warning('NO_PAIRS'));
  return response;
}

// design(body, deps) -> Promise<response>
//   deps (all optional): {cfg, log, signal (client abort), semaphore {acquire}, now, primer3 {run, version},
//     template {buildTemplate}, repeatMask {repeatMask}, and what template/repeat_mask take:
//     findGene | mongo, resolve | catalog, sequence, runBlast, cache, span_fetch_limit}
async function design(body, deps) {
  deps = deps || {};
  const cfg = deps.cfg || require('./config').get();
  if (cfg.enabled === false) throw new PrimerHttpError(503, 'FEATURE_DISABLED', 'primer design is disabled on this server', {});
  const req = normalize(body, cfg);
  const semaphore = deps.semaphore || require('./semaphore').designSemaphore();
  const release = await semaphore.acquire({ signal: deps.signal });
  let dl = null;
  try {
    dl = createDeadline(cfg.design.deadline_ms, deps.signal, deps.now);
    return await runDesign(req, cfg, deps, dl);
  } finally {
    if (dl) dl.dispose();
    release();
  }
}

module.exports = {
  design,
  normalize,
  parseParams,
  checkParams,
  validateAgainstTemplate,
  junctionList,
  buildRecord,
  decoratePairs,
  longSeqTm,
  unmaskProductTm,
  publicTemplate,
  createDeadline,
  PRESETS,
  PARAM_SPECS,
  ALWAYS_TAGS,
  PARAM_TAGS,
  MODES,
  MASK_MODES,
  MAX_JUNCTIONS,
  MESSAGES
};
