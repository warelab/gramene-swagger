'use strict';

// The genotyping template (spec §4.3, §4.5): one plus-strand REF window around a resolved variant, its ALT haplotype,
// the template features, the per-orientation room and N checks, and the repeat-mask exemption of the
// allele-specific windows.
//
// Template coordinates are 1-based inclusive on the plus strand, t = g - template.start + 1. alt_seq indexes like seq
// up to the variant and is shifted by alt_offset = len(vcf.alt) - len(vcf.ref) after it (§2.1). Intervals use the
// PrimerInterval [start, length] encoding.

const templates = require('../template');
const { PrimerHttpError } = require('../errors');

// §4.3: flank = max(template_flank, largest product-range upper bound + FLANK_PRODUCT_PAD).
const FLANK_PRODUCT_PAD = 40;

function isPosInt(v) {
  return Number.isSafeInteger(v) && v >= 1;
}

function interval(start, end) {
  return [start, end - start + 1];
}

// Over the levels a design may run: the smallest product minimum, the largest product maximum and the largest
// max_size (product minimums already raised by presets.levelParams).
function levelBounds(levels) {
  if (!Array.isArray(levels) || levels.length === 0) throw new TypeError('deps.req.levels must list the effective params of each level');
  let productMin = Infinity;
  let productMax = 0;
  let maxSize = 0;
  levels.forEach(function (p) {
    p.product_size_ranges.forEach(function (r) {
      productMin = Math.min(productMin, r[0]);
      productMax = Math.max(productMax, r[1]);
    });
    maxSize = Math.max(maxSize, p.max_size);
  });
  return { productMin: productMin, productMax: productMax, maxSize: maxSize };
}

// §4.3 flank F = max(genotyping.template_flank, the largest product-range upper bound over the levels + 40).
function templateFlank(levels, templateFlankCfg) {
  return Math.max(templateFlankCfg, levelBounds(levels).productMax + FLANK_PRODUCT_PAD);
}

// §4.3 template window of a canonical entry: [max(1, min(fwd, rev) - F), min(region_length, max(fwd, rev) + F)]. The
// design fetches the neighbours of this window before the semaphore, so it is computed from the entry alone.
function templateWindow(variant, regionLength, flank) {
  const fwd = variant.discriminating.forward.position;
  const rev = variant.discriminating.reverse.position;
  return { start: Math.max(1, Math.min(fwd, rev) - flank), end: Math.min(regionLength, Math.max(fwd, rev) + flank) };
}

function checkVariant(variant, cfg) {
  if (!variant || typeof variant.region !== 'string' || !variant.vcf || !isPosInt(variant.vcf.position) ||
      typeof variant.vcf.ref !== 'string' || typeof variant.vcf.alt !== 'string') {
    throw new TypeError('variant must be a canonical entry with region and vcf {position, ref, alt}');
  }
  if (!variant.zone || !variant.discriminating) {
    // recordsToEntries leaves both null when the event slides more than variation.max_shift (REPEAT_TOO_LONG).
    const max = cfg.variation && Number.isSafeInteger(cfg.variation.max_shift) ? cfg.variation.max_shift : null;
    throw new PrimerHttpError(400, 'VARIANT_TOO_REPETITIVE',
      'the variant at ' + variant.region + ':' + variant.vcf.position + ' can slide too far to have a stable discriminating base',
      { region: variant.region, position: variant.vcf.position, shift: variant.shift === undefined ? null : variant.shift, max: max });
  }
}

// buildVariantTemplate(variant, resolved, deps) -> Promise<vt>
//   variant: a canonical entry (variation/normalize.js: region, vcf, zone, discriminating), REF already verified
//   resolved: the assembly (assemblies.resolve); fasta.dna is read
//   deps: {cfg, req (genotyping/request.js normalize; its levels), sequence {fetch, regionLength}, log}; passed on to
//         template.buildRegionTemplate, whose assembly lookup is answered with `resolved`
// vt: {id, system_name, region, start, end, strand: 1, length, alt_length, seq, alt_seq (uppercase),
//      region_length, features {variant, zone, discriminating {forward, reverse}, alt_offset, exempt [fwd, rev]},
//      product_min, skip {forward, reverse: null | 'too_close_to_end' | 'n_in_primer_window'},
//      region_template (the region-mode template, for repeat_mask.repeatMask), warnings}
// Throws 422 NO_SEQUENCE, 404 UNKNOWN_REGION, 400 VARIANT_TOO_REPETITIVE, 400 VARIANT_TOO_CLOSE_TO_END (neither
// orientation has room for the smallest product), 400 REF_MISMATCH (the template disagrees with vcf.ref).
async function buildVariantTemplate(variant, resolved, deps) {
  deps = deps || {};
  const cfg = deps.cfg || require('../config').get();
  const g = cfg.genotyping || {};
  ['template_flank', 'mask_exempt_pad'].forEach(function (k) {
    if (!Number.isSafeInteger(g[k]) || g[k] < 0) throw new TypeError('primers config: genotyping.' + k + ' must be a non-negative integer');
  });
  const sequence = deps.sequence || require('../sequence');
  checkVariant(variant, cfg);
  const fasta = resolved && resolved.fasta && resolved.fasta.dna;
  if (!fasta) {
    const name = resolved && resolved.system_name;
    throw new PrimerHttpError(422, 'NO_SEQUENCE', 'no genome sequence is available for ' + name, { system_name: name });
  }
  const bounds = levelBounds(deps.req && deps.req.levels);
  const vcf = variant.vcf;
  const fwd = variant.discriminating.forward.position;
  const rev = variant.discriminating.reverse.position;

  const regionLength = await templates.requireRegionLength(sequence, fasta, variant.region);
  const win = templateWindow(variant, regionLength, templateFlank(deps.req.levels, g.template_flank));
  const start = win.start;
  const end = win.end;

  // Room (§4.3), before any sequence is read.
  const skip = { forward: null, reverse: null };
  if (end - variant.zone.end + 1 < bounds.productMin) skip.forward = 'too_close_to_end';
  if (variant.zone.start - start + 1 < bounds.productMin) skip.reverse = 'too_close_to_end';
  if (skip.forward && skip.reverse) {
    throw new PrimerHttpError(400, 'VARIANT_TOO_CLOSE_TO_END',
      'the variant at ' + variant.region + ':' + vcf.position + ' is too close to the end of the sequence for a ' +
      bounds.productMin + ' bp product in either orientation',
      { region: variant.region, position: vcf.position, region_length: regionLength, needed: bounds.productMin });
  }

  const tdeps = Object.assign({}, deps, { cfg: cfg, sequence: sequence, resolve: async function () { return resolved; } });
  const tpl = await templates.buildRegionTemplate({
    system_name: resolved.system_name,
    region: { region: variant.region, start: start, end: end, strand: 1 }
  }, tdeps);
  const seq = tpl.seq.toUpperCase();
  const length = seq.length;
  const t = function (gpos) { return gpos - start + 1; };

  const tv = t(vcf.position);
  const genomeRef = seq.slice(tv - 1, tv - 1 + vcf.ref.length);
  if (genomeRef !== vcf.ref) {
    throw new PrimerHttpError(400, 'REF_MISMATCH', 'the reference allele ' + vcf.ref + ' does not match the genome bases ' +
      genomeRef + ' at ' + variant.region + ':' + vcf.position,
    { region: variant.region, position: vcf.position, given: vcf.ref, genome: genomeRef });
  }
  const altSeq = seq.slice(0, tv - 1) + vcf.alt + seq.slice(tv - 1 + vcf.ref.length);

  const tFwd = t(fwd);
  const tRev = t(rev);
  const zone = { start: t(variant.zone.start), end: t(variant.zone.end) };
  const pad = g.mask_exempt_pad;
  const exempt = [
    interval(Math.max(1, tFwd - pad), zone.end),
    interval(zone.start, Math.min(length, tRev + pad))
  ];

  // An N in the allele-specific window of the longest primer any level allows skips that orientation (§4.3).
  if (!skip.forward && seq.slice(Math.max(0, tFwd - bounds.maxSize), zone.end).indexOf('N') >= 0) skip.forward = 'n_in_primer_window';
  if (!skip.reverse && seq.slice(zone.start - 1, Math.min(length, tRev + bounds.maxSize - 1)).indexOf('N') >= 0) {
    skip.reverse = 'n_in_primer_window';
  }

  return {
    id: tpl.id,
    system_name: tpl.system_name,
    region: variant.region,
    start: start,
    end: end,
    strand: 1,
    length: length,
    alt_length: altSeq.length,
    seq: seq,
    alt_seq: altSeq,
    region_length: regionLength,
    features: {
      variant: { start: tv, end: tv + vcf.ref.length - 1 },
      zone: zone,
      discriminating: { forward: tFwd, reverse: tRev },
      alt_offset: vcf.alt.length - vcf.ref.length,
      exempt: exempt
    },
    product_min: bounds.productMin,
    skip: skip,
    region_template: tpl,
    warnings: tpl.warnings.slice()
  };
}

// The Primer3 template of one orientation (§4.5): the masked template with that orientation's allele-specific window
// (features.exempt) restored to the REF bases. maskedSeq null (no repeat mask) -> vt.seq.
// -> {seq, exempt [start, length], masked_bases (exempted bases the mask had changed)}
function orientationTemplate(vt, maskedSeq, orientation) {
  if (orientation !== 'forward' && orientation !== 'reverse') throw new TypeError('orientation must be forward or reverse');
  const span = vt.features.exempt[orientation === 'forward' ? 0 : 1];
  if (maskedSeq === null || maskedSeq === undefined) return { seq: vt.seq, exempt: span, masked_bases: 0 };
  if (typeof maskedSeq !== 'string' || maskedSeq.length !== vt.length) throw new RangeError('the masked template must have the template length');
  const s = span[0] - 1;
  const e = s + span[1];
  let masked = 0;
  for (let i = s; i < e; i++) if (maskedSeq[i] !== vt.seq[i]) masked++;
  return { seq: maskedSeq.slice(0, s) + vt.seq.slice(s, e) + maskedSeq.slice(e), exempt: span, masked_bases: masked };
}

// PrimerGenotypingTemplate. mask: {masked, mask_source, mask, masked_fraction} (repeat_mask.repeatMask, or unmasked).
function publicVariantTemplate(vt, mask) {
  const m = mask || { masked: false, mask_source: null, mask: [], masked_fraction: 0 };
  return {
    system_name: vt.system_name,
    region: vt.region,
    start: vt.start,
    end: vt.end,
    strand: 1,
    length: vt.length,
    alt_length: vt.alt_length,
    seq: vt.seq,
    alt_seq: vt.alt_seq,
    masked: m.masked,
    mask_source: m.mask_source,
    mask: m.mask,
    masked_fraction: m.masked_fraction,
    features: vt.features
  };
}

module.exports = {
  buildVariantTemplate,
  orientationTemplate,
  publicVariantTemplate,
  levelBounds,
  templateFlank,
  templateWindow,
  FLANK_PRODUCT_PAD
};
