'use strict';

// POST /primers/genotyping/design request normalization (spec §2.7 handler rules 1-6). Pure: no I/O, so it runs
// before the deadline-bound work. Rules 7-15 need the catalog, Ensembl or the FASTA and run later; an allele that
// is '*' or contains N is syntactically accepted here and rejected by rule 12 (UNSUPPORTED_ALLELE) after rule 7.

const design = require('../design');
const presets = require('./presets');
const { PrimerHttpError } = require('../errors');

const SYSTEM_NAME_RE = /^[a-z0-9_]+$/;
const MAX_SYSTEM_NAME_LENGTH = 128;
const MAX_REGION_LENGTH = 255;
const ID_RE = /^[A-Za-z0-9][A-Za-z0-9_.:,*-]{0,254}$/;
const LABEL_RE = /^[A-Za-z0-9_.-]{1,40}$/;
const DEFAULT_MAX_ALLELE_LENGTH = 50;
const DEFAULT_MAX_SETS = 10;
const MASK_MODES = Object.freeze(['n_mask', 'three_prime']);

const TOP_LEVEL_FIELDS = Object.freeze(['system_name', 'variant', 'assay', 'avoid_repeats', 'repeat_mask_mode',
  'template_only', 'label', 'params']);
const VARIANT_FIELDS = Object.freeze(['id', 'region', 'position', 'ref', 'alt']);
const ASSAY_FIELDS = Object.freeze(['type', 'orientation', 'tails', 'deliberate_mismatch', 'mismatch_position', 'num_sets',
  'max_relaxation', 'neighbour_policy']);

// PrimerGenotypingParams: the closed subset of PrimerDesignParams a genotyping design accepts.
const PARAM_KEYS = Object.freeze(['opt_size', 'min_size', 'max_size', 'opt_tm', 'min_tm', 'max_tm', 'opt_gc', 'min_gc',
  'max_gc', 'max_tm_diff', 'max_poly_x', 'gc_clamp', 'max_end_stability', 'salt_monovalent', 'salt_divalent', 'dntp_conc',
  'dna_conc', 'product_size_ranges']);
// Narrower than design.PARAM_SPECS.product_size_ranges (10 ranges up to 50,000 bp).
const PRODUCT_RANGES = Object.freeze({ maxItems: 4, min: 20, max: 1000 });

const ENUMS = Object.freeze({
  type: Object.freeze(['kasp', 'as_pcr']),
  orientation: Object.freeze(['both', 'forward', 'reverse']),
  tails: Object.freeze(['none', 'ref_fam_alt_hex', 'ref_hex_alt_fam']),
  deliberate_mismatch: Object.freeze(['none', 'auto']),
  mismatch_position: Object.freeze([2, 3]),
  neighbour_policy: Object.freeze(['avoid_3p', 'ignore'])
});

// Effective assay defaults per type (§2.7 PrimerGenotypingAssay); num_sets comes from genotyping.num_sets_default
// when configured.
const ASSAY_DEFAULTS = Object.freeze({
  kasp: Object.freeze({ type: 'kasp', orientation: 'both', tails: 'ref_fam_alt_hex', deliberate_mismatch: 'none',
    mismatch_position: 2, num_sets: 6, max_relaxation: 2, neighbour_policy: 'avoid_3p' }),
  as_pcr: Object.freeze({ type: 'as_pcr', orientation: 'both', tails: 'none', deliberate_mismatch: 'auto',
    mismatch_position: 2, num_sets: 6, max_relaxation: 2, neighbour_policy: 'avoid_3p' })
});
const DEFAULT_TYPE = 'kasp';

function absent(v) {
  return v === undefined || v === null;
}

function isObject(v) {
  return v !== null && typeof v === 'object' && !Array.isArray(v);
}

function invalidRequest(message, details) {
  return new PrimerHttpError(400, 'INVALID_REQUEST', message, details || {});
}

function invalidParams(message, details) {
  return new PrimerHttpError(400, 'INVALID_PARAMS', message, details || {});
}

function invalidVariant(reason, message) {
  return new PrimerHttpError(400, 'INVALID_VARIANT', message, { reason: reason });
}

// Rule 1 (defence in depth; swagger rejects unknown keys first).
function checkKeys(obj, allowed, prefix) {
  Object.keys(obj).forEach(function (k) {
    if (allowed.indexOf(k) < 0) {
      const field = prefix + k.slice(0, 64);
      throw invalidRequest('unknown field ' + JSON.stringify(field), { field: field });
    }
  });
}

function optBool(obj, field) {
  const v = obj[field];
  if (absent(v)) return undefined;
  if (typeof v !== 'boolean') throw invalidRequest(field + ' must be true or false', { field: field });
  return v;
}

function optEnum(obj, key, prefix) {
  const v = obj[key];
  if (absent(v)) return undefined;
  if (ENUMS[key].indexOf(v) < 0) {
    throw invalidRequest(prefix + key + ' must be one of ' + ENUMS[key].join(', '), { field: prefix + key });
  }
  return v;
}

function optInt(obj, key, min, max, prefix) {
  const v = obj[key];
  if (absent(v)) return undefined;
  if (!Number.isSafeInteger(v) || v < min || v > max) {
    throw invalidRequest(prefix + key + ' must be an integer from ' + min + ' to ' + max, { field: prefix + key, min: min, max: max });
  }
  return v;
}

// 'empty' ('-'), 'bases' (A/C/G/T), 'unsupported' ('*' or containing N: rule 12, after I/O), 'invalid'.
function alleleClass(a) {
  if (a === '-') return 'empty';
  if (/^[ACGT]+$/.test(a)) return 'bases';
  if (a === '*' || /^[ACGTN]+$/.test(a)) return 'unsupported';
  return 'invalid';
}

function alleleString(v, field) {
  if (typeof v !== 'string' || v === '') throw invalidRequest(field + ' must be a non-empty string', { field: field });
  return v.toUpperCase();
}

// Rules 2 and 3. -> {input: 'id', id, alt | null} | {input: 'manual', region, position, ref, alt} (alleles uppercase)
function parseVariant(v, maxAllele) {
  if (!isObject(v)) {
    throw invalidRequest('variant is required: {id [, alt]} or {region, position, ref, alt}', { field: 'variant' });
  }
  const has = function (k) { return !absent(v[k]); };
  if (has('id')) {
    if (has('region') || has('position') || has('ref')) {
      throw invalidVariant('id_or_manual', 'variant takes either id (with an optional alt) or region, position, ref and alt');
    }
    if (typeof v.id !== 'string' || !ID_RE.test(v.id)) {
      throw invalidRequest('variant.id must match ' + ID_RE.source, { field: 'variant.id' });
    }
    let alt = null;
    if (has('alt')) {
      alt = alleleString(v.alt, 'variant.alt');
      if (alleleClass(alt) === 'invalid') throw invalidRequest('variant.alt must be A/C/G/T bases or "-"', { field: 'variant.alt' });
      if (alt.length > maxAllele) throw invalidVariant('allele_too_long', 'alleles are limited to ' + maxAllele + ' nt');
    }
    return { input: 'id', id: v.id, alt: alt };
  }
  if (!(has('region') && has('position') && has('ref') && has('alt'))) {
    throw invalidVariant('id_or_manual', 'variant needs either id (with an optional alt) or all of region, position, ref and alt');
  }
  if (typeof v.region !== 'string' || v.region === '' || v.region.length > MAX_REGION_LENGTH) {
    throw invalidRequest('variant.region must be a non-empty string of at most ' + MAX_REGION_LENGTH + ' characters', { field: 'variant.region' });
  }
  if (!Number.isSafeInteger(v.position) || v.position < 1) {
    throw invalidRequest('variant.position must be an integer >= 1', { field: 'variant.position' });
  }
  const ref = alleleString(v.ref, 'variant.ref');
  const alt = alleleString(v.alt, 'variant.alt');
  const rc = alleleClass(ref);
  const ac = alleleClass(alt);
  if (rc === 'invalid' || ac === 'invalid' || ref === alt || (rc === 'empty' && ac !== 'bases') || (ac === 'empty' && rc !== 'bases')) {
    throw invalidVariant('alleles', 'ref and alt must be different A/C/G/T alleles, and at most one of them may be "-"');
  }
  if (ref.length > maxAllele || alt.length > maxAllele) {
    throw invalidVariant('allele_too_long', 'alleles are limited to ' + maxAllele + ' nt');
  }
  return { input: 'manual', region: v.region, position: v.position, ref: ref, alt: alt };
}

// Effective assay: the client's values over the defaults of its type (rule 6: any combination is allowed).
function parseAssay(a, g) {
  const src = isObject(a) ? a : {};
  const p = 'assay.';
  const type = optEnum(src, 'type', p) || DEFAULT_TYPE;
  const d = ASSAY_DEFAULTS[type];
  const maxSets = Number.isSafeInteger(g.max_sets) && g.max_sets >= 1 ? g.max_sets : DEFAULT_MAX_SETS;
  const defaultSets = Number.isSafeInteger(g.num_sets_default) ? g.num_sets_default : d.num_sets;
  const pick = function (v, dflt) { return v === undefined ? dflt : v; };
  const mismatchPosition = src.mismatch_position;
  if (!absent(mismatchPosition) && ENUMS.mismatch_position.indexOf(mismatchPosition) < 0) {
    throw invalidRequest('assay.mismatch_position must be 2 or 3', { field: 'assay.mismatch_position' });
  }
  return {
    type: type,
    orientation: pick(optEnum(src, 'orientation', p), d.orientation),
    tails: pick(optEnum(src, 'tails', p), d.tails),
    deliberate_mismatch: pick(optEnum(src, 'deliberate_mismatch', p), d.deliberate_mismatch),
    mismatch_position: absent(mismatchPosition) ? d.mismatch_position : mismatchPosition,
    num_sets: pick(optInt(src, 'num_sets', 1, maxSets, p), Math.min(defaultSets, maxSets)),
    max_relaxation: pick(optInt(src, 'max_relaxation', 0, presets.MAX_LEVEL, p), d.max_relaxation),
    neighbour_policy: pick(optEnum(src, 'neighbour_policy', p), d.neighbour_policy)
  };
}

function parseUserParams(params) {
  if (absent(params)) return {};
  const p = design.parseParams(params);
  if (p.product_size_ranges !== undefined) {
    const ranges = p.product_size_ranges;
    if (ranges.length > PRODUCT_RANGES.maxItems) {
      throw invalidParams('product_size_ranges must be a list of 1 to ' + PRODUCT_RANGES.maxItems + ' [min, max] ranges',
        { param: 'product_size_ranges' });
    }
    ranges.forEach(function (r, i) {
      if (!r.every(function (x) { return x >= PRODUCT_RANGES.min && x <= PRODUCT_RANGES.max; })) {
        throw invalidParams('product_size_ranges[' + i + '] must be [min, max] with integers from ' + PRODUCT_RANGES.min + ' to ' +
          PRODUCT_RANGES.max, { param: 'product_size_ranges', index: i });
      }
    });
  }
  return p;
}

// Rule 5 at every level the design may run, so a relaxation level never reaches Primer3 with crossed bounds.
function checkLevel(params, level) {
  try {
    design.checkParams(params);
  } catch (e) {
    if (level === 0 || !(e instanceof PrimerHttpError)) throw e;
    throw invalidParams(e.message + ' at relaxation level ' + level, Object.assign({}, e.details, { level: level }));
  }
}

function numberOr(v, dflt) {
  return typeof v === 'number' && Number.isFinite(v) ? v : dflt;
}

// normalize(body, cfg) -> req. Throws 400 INVALID_REQUEST / INVALID_VARIANT / INVALID_PARAMS.
//   req: {system_name, variant_input: 'id' | 'manual', variant, assay (effective, §2.8 key order),
//         orientations ['forward', 'reverse'] or one of them, avoid_repeats, repeat_mask_mode (null unless
//         avoid_repeats), template_only, label (null when absent), preset, user_params, pinned (PARAM_KEYS order),
//         levels [params of level 0 .. assay.max_relaxation], params (= levels[0]), ladder [{level, changes}],
//         floors {as_min_tm, as_min_gc}}
//   cfg: the primers config; reads genotyping.{num_sets_default, max_sets, as_min_tm, as_min_gc} and
//        variation.max_allele_length.
function normalize(body, cfg) {
  cfg = cfg || require('../config').get();
  const g = cfg.genotyping || {};
  const maxAllele = Number.isSafeInteger(cfg.variation && cfg.variation.max_allele_length)
    ? cfg.variation.max_allele_length : DEFAULT_MAX_ALLELE_LENGTH;
  if (!isObject(body)) throw invalidRequest('the request body must be a JSON object');
  checkKeys(body, TOP_LEVEL_FIELDS, '');
  if (isObject(body.variant)) checkKeys(body.variant, VARIANT_FIELDS, 'variant.');
  if (isObject(body.assay)) checkKeys(body.assay, ASSAY_FIELDS, 'assay.');
  if (isObject(body.params)) checkKeys(body.params, PARAM_KEYS, 'params.');
  if (!absent(body.assay) && !isObject(body.assay)) throw invalidRequest('assay must be an object', { field: 'assay' });

  const systemName = body.system_name;
  if (typeof systemName !== 'string' || systemName.length > MAX_SYSTEM_NAME_LENGTH || !SYSTEM_NAME_RE.test(systemName)) {
    throw invalidRequest('system_name is required and must match ^[a-z0-9_]+$ (at most ' + MAX_SYSTEM_NAME_LENGTH + ' characters)',
      { field: 'system_name' });
  }
  const variant = parseVariant(body.variant, maxAllele);
  const assay = parseAssay(body.assay, g);
  const avoidRepeats = optBool(body, 'avoid_repeats') === true;
  const maskMode = absent(body.repeat_mask_mode) ? undefined : body.repeat_mask_mode;
  if (maskMode !== undefined && MASK_MODES.indexOf(maskMode) < 0) {
    throw invalidRequest('repeat_mask_mode must be n_mask or three_prime', { field: 'repeat_mask_mode' });
  }
  const templateOnly = optBool(body, 'template_only') === true;
  let label = null;
  if (!absent(body.label)) {
    if (typeof body.label !== 'string' || !LABEL_RE.test(body.label)) {
      throw invalidRequest('label must match ' + LABEL_RE.source, { field: 'label' });
    }
    label = body.label;
  }

  const userParams = parseUserParams(body.params);
  const levels = [];
  for (let level = 0; level <= assay.max_relaxation; level++) {
    const p = presets.levelParams(assay.type, level, userParams);
    checkLevel(p, level);
    levels.push(p);
  }
  const variantOut = Object.assign({}, variant);
  delete variantOut.input;
  return {
    system_name: systemName,
    variant_input: variant.input,
    variant: variantOut,
    assay: assay,
    orientations: assay.orientation === 'both' ? ['forward', 'reverse'] : [assay.orientation],
    avoid_repeats: avoidRepeats,
    repeat_mask_mode: avoidRepeats ? maskMode || 'n_mask' : null,
    template_only: templateOnly,
    label: label,
    preset: assay.type,
    user_params: userParams,
    pinned: PARAM_KEYS.filter(function (k) { return userParams[k] !== undefined; }),
    levels: levels,
    params: levels[0],
    ladder: levels.slice(1).map(function (p, i) { return { level: i + 1, changes: presets.paramChanges(levels[i], p) }; }),
    floors: {
      as_min_tm: numberOr(g.as_min_tm, presets.FLOORS.as_min_tm),
      as_min_gc: numberOr(g.as_min_gc, presets.FLOORS.as_min_gc)
    }
  };
}

module.exports = {
  normalize,
  ASSAY_DEFAULTS,
  PARAM_KEYS,
  ENUMS,
  TOP_LEVEL_FIELDS,
  VARIANT_FIELDS,
  ASSAY_FIELDS,
  PRODUCT_RANGES
};
