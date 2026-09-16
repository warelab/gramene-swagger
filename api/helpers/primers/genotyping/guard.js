'use strict';

// The common-primer guard of a genotyping design run (spec §4.4). targetTag builds the SEQUENCE_TARGET that keeps
// Primer3's common primer at least guard_gap bases outside the variant zone, on the side where it must sit, while the
// forced allele-specific primer on the other side of the target is untouched. commonPrimerOk is the post-filter that
// re-checks each returned pair on both haplotypes. The zone already holds both discriminating positions and both
// ALT 3' anchor mappings (§3.5), so the guard also protects the ALT haplotype.
//
// Coordinates are genotyping template coordinates (genotyping/template.js): 1-based inclusive on the plus strand.

function assertInt(v, name) {
  if (!Number.isSafeInteger(v)) throw new TypeError(name + ' must be an integer');
}

function assertSpan(s, name) {
  if (!s || typeof s !== 'object') throw new TypeError(name + ' must be {start, end}');
  assertInt(s.start, name + '.start');
  assertInt(s.end, name + '.end');
  if (s.end < s.start) throw new RangeError(name + ' must have start <= end');
}

function assertOrientation(o) {
  if (o !== 'forward' && o !== 'reverse') throw new TypeError('orientation must be forward or reverse');
}

function intersects(a, b) {
  return a.start <= b.end && b.start <= a.end;
}

function shifted(s, delta) {
  return { start: s.start + delta, end: s.end + delta };
}

// targetTag(zone, orientation, template, gap) -> "start,length" for SEQUENCE_TARGET, or null when the template edge
// leaves no base on the common-primer side.
//   zone: template.features.zone; template: {length}; gap: genotyping.guard_gap
// forward: the gap bases right after the zone; reverse: the gap bases right before it. The gap is shortened (never
// below 1) when the template edge is closer.
function targetTag(zone, orientation, template, gap) {
  assertSpan(zone, 'zone');
  assertOrientation(orientation);
  assertInt(template && template.length, 'template.length');
  assertInt(gap, 'gap');
  if (gap < 1) throw new RangeError('gap must be at least 1');
  if (zone.start < 1 || zone.end > template.length) throw new RangeError('the zone must lie inside the template');
  if (orientation === 'forward') {
    const len = Math.min(gap, template.length - zone.end);
    return len < 1 ? null : (zone.end + 1) + ',' + len;
  }
  const len = Math.min(gap, zone.start - 1);
  return len < 1 ? null : (zone.start - len) + ',' + len;
}

// commonPrimerOk(set, zone, haplotypes) -> boolean
//   set: {orientation, as: {start, end}, common: {start, end}}: REF template footprints of the allele-specific and
//        the common primer of one design-run pair
//   zone: template.features.zone
//   haplotypes: {delta} with delta = len(vcf.alt) - len(vcf.ref), or the variant template itself ({length, alt_length})
// Rejects a common primer that
//   1. touches [zone.start - 1, zone.end + 1] on the REF haplotype, or
//   2. overlaps the allele-specific primer on either haplotype. On the ALT haplotype the bases after the variant move
//      by delta: a reverse allele-specific primer moves (§4.9), and so does a common primer downstream of the zone.
// The third post-filter rule (the ALT check_primers run must find the common primer) belongs to scoring.
function commonPrimerOk(set, zone, haplotypes) {
  if (!set || typeof set !== 'object') throw new TypeError('set must be {orientation, as, common}');
  assertOrientation(set.orientation);
  assertSpan(set.as, 'set.as');
  assertSpan(set.common, 'set.common');
  assertSpan(zone, 'zone');
  const h = haplotypes || {};
  const delta = Number.isSafeInteger(h.delta) ? h.delta : h.alt_length - h.length;
  assertInt(delta, 'haplotypes.delta');

  const common = set.common;
  if (intersects(common, { start: zone.start - 1, end: zone.end + 1 })) return false;
  if (intersects(common, set.as)) return false;
  const asAlt = set.orientation === 'reverse' ? shifted(set.as, delta) : set.as;
  const commonAlt = common.start > zone.end ? shifted(common, delta) : common;
  return !intersects(commonAlt, asAlt);
}

module.exports = {
  targetTag,
  commonPrimerOk
};
