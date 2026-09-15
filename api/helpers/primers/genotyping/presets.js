'use strict';

// Genotyping presets, the relaxation ladder and the hard floors (spec §4.6). Frozen tables and pure helpers.
//
// Level 0 is the preset. Each ladder level applies its changes on top of the previous level, except the params the
// client set: those are pinned and win at every level. At every level the effective product-range minimum is raised
// to 2 x max_size + 1, so a pair can never be shorter than its two primers.

function deepFreeze(o) {
  if (o !== null && typeof o === 'object' && !Object.isFrozen(o)) {
    Object.freeze(o);
    Object.keys(o).forEach(function (k) { deepFreeze(o[k]); });
  }
  return o;
}

function clone(v) {
  return JSON.parse(JSON.stringify(v));
}

// Level-0 params per assay type (§4.6 tables; sources there).
const PRESETS = deepFreeze({
  kasp: {
    opt_size: 22, min_size: 18, max_size: 30,
    opt_tm: 60, min_tm: 57, max_tm: 63,
    min_gc: 30, max_gc: 70,
    max_tm_diff: 3, max_poly_x: 5,
    product_size_ranges: [[50, 120]]
  },
  as_pcr: {
    opt_size: 24, min_size: 18, max_size: 30,
    opt_tm: 60, min_tm: 57, max_tm: 63,
    min_gc: 30, max_gc: 70,
    max_tm_diff: 3, max_poly_x: 4,
    product_size_ranges: [[150, 300]]
  }
});

// LADDER[type][level]: the changes that level applies on top of the previous one (level 0 has none).
const LADDER = deepFreeze({
  kasp: [
    {},
    { max_size: 32, min_tm: 55, max_tm: 65, min_gc: 20, max_gc: 80, product_size_ranges: [[50, 150]] },
    { min_tm: 52, max_tm_diff: 6 }
  ],
  as_pcr: [
    {},
    { max_size: 32, min_tm: 55, max_tm: 65, min_gc: 20, max_gc: 80 },
    { min_tm: 52, max_tm_diff: 6 }
  ]
});

const MAX_LEVEL = 2;

// Never relaxed (§3.1 genotyping.as_min_tm / as_min_gc): Tm of every allele-specific primer (matched_tm for a
// deliberate-mismatch primer), GC of the derived ones.
const FLOORS = deepFreeze({ as_min_tm: 52, as_min_gc: 15 });

// levelParams(type, level, userParams) -> effective params of one ladder level (a fresh object).
//   userParams: the client's validated params (genotyping/request.js); each one is pinned.
function levelParams(type, level, userParams) {
  if (!Object.prototype.hasOwnProperty.call(PRESETS, type)) throw new TypeError('unknown assay type ' + JSON.stringify(type));
  if (!Number.isInteger(level) || level < 0 || level > MAX_LEVEL) throw new RangeError('level must be an integer from 0 to ' + MAX_LEVEL);
  const p = clone(PRESETS[type]);
  for (let l = 1; l <= level; l++) Object.assign(p, clone(LADDER[type][l]));
  const user = userParams || {};
  Object.keys(user).forEach(function (k) {
    if (user[k] !== undefined) p[k] = clone(user[k]);
  });
  p.product_size_ranges = p.product_size_ranges.map(function (r) { return [Math.max(r[0], 2 * p.max_size + 1), r[1]]; });
  return p;
}

// The effective changes from one level's params to the next, in the order of `next` (§2.8 attempts[].changes,
// settings.ladder[].changes). Pinned params never differ, so they never appear.
function paramChanges(prev, next) {
  const out = {};
  Object.keys(next).forEach(function (k) {
    if (JSON.stringify(prev[k]) !== JSON.stringify(next[k])) out[k] = clone(next[k]);
  });
  return out;
}

module.exports = {
  PRESETS,
  LADDER,
  FLOORS,
  MAX_LEVEL,
  levelParams,
  paramChanges
};
