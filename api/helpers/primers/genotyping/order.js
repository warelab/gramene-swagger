'use strict';

// Order sheet of a genotyping design (spec §4.17): oligo names, the three order rows of a set and the KASP primer
// mix. Pure.

const decimal = require('./decimal');
const mismatch = require('./mismatch');

const LABEL_MAX = 40;
const ROLE_SUFFIX = Object.freeze({ as_ref: 'REF', as_alt: 'ALT', common: 'COM' });
const ROLES = Object.freeze(['as_ref', 'as_alt', 'common']);

// KASP primer mix: the 12:12:30 ratio at 100 uM is from Makhoul et al. 2020; the water to 100 uL is arithmetic.
const KASP_MIX = Object.freeze({ stock_uM: 100, as_ref_uL: 12, as_alt_uL: 12, common_uL: 30, water_uL: 46, total_uL: 100,
  source: 'Makhoul et al. 2020 (12:12:30 at 100 uM); water to 100 uL inferred' });

// [^A-Za-z0-9_.-] -> _, at most 40 characters.
function sanitizeLabel(s) {
  return String(s).replace(/[^A-Za-z0-9_.-]/g, '_').slice(0, LABEL_MAX);
}

// §2.1 label: the request label, else the first variant id, else the key; sanitized.
function labelFor(label, variant) {
  const raw = label || (variant && Array.isArray(variant.ids) && variant.ids.length ? variant.ids[0] : null) || (variant && variant.key);
  return sanitizeLabel(raw);
}

// {label}_{set id}_{REF|ALT|COM}, plus _FAM/_HEX for a tailed allele-specific primer.
function oligoName(label, setId, role, dye) {
  if (!Object.prototype.hasOwnProperty.call(ROLE_SUFFIX, role)) throw new TypeError('role must be as_ref, as_alt or common');
  return label + '_' + setId + '_' + ROLE_SUFFIX[role] + (role !== 'common' && dye ? '_' + dye : '');
}

// The row notes: the deliberate mismatch, and the tailed hairpin rounded once to 1 decimal whenever the oligo has a
// TAILED_STRUCTURE hairpin issue of either severity. issues: the set's internal issues (value_exact kept).
function notesFor(role, oligo, issues) {
  const parts = [];
  if (oligo.deliberate_mismatch) parts.push(mismatch.note(oligo.deliberate_mismatch));
  const hairpin = (issues || []).find(function (i) {
    return i.code === 'TAILED_STRUCTURE' && i.details.oligo === role && i.details.metric === 'hairpin_th';
  });
  if (hairpin) parts.push('tail hairpin ' + decimal.round(decimal.micro(hairpin.value_exact), 1).toFixed(1) + ' °C');
  return parts.join('; ');
}

// orderRows(set, {label, variant_key}) -> the three PrimerGenotypingOrderRow objects, REF, ALT, common.
//   set: {id, key, orientation, primers {as_ref, as_alt, common}, products {ref, alt}, issues (internal)}
function orderRows(set, opts) {
  return ROLES.map(function (role) {
    const o = set.primers[role];
    return {
      name: oligoName(opts.label, set.id, role, o.dye),
      set_id: set.id,
      set_key: set.key,
      role: role,
      allele: o.allele,
      dye: o.dye,
      order_seq: o.order_seq,
      target_seq: o.target_seq,
      tail_seq: o.tail_seq,
      length: o.order_len,
      tm: o.tm,
      gc: o.gc,
      orientation: set.orientation,
      product_size_ref: set.products.ref.size,
      product_size_alt: set.products.alt.size,
      variant_key: opts.variant_key,
      notes: notesFor(role, o, set.issues)
    };
  });
}

// assay.kasp_mix: the mix for type kasp, null otherwise.
function kaspMix(type) {
  return type === 'kasp' ? Object.assign({}, KASP_MIX) : null;
}

module.exports = {
  sanitizeLabel,
  labelFor,
  oligoName,
  orderRows,
  notesFor,
  kaspMix,
  KASP_MIX,
  LABEL_MAX
};
