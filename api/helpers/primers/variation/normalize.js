'use strict';

// Variant resolution, pure half (spec §3.4, §3.5, §3.9, §4.17). Ensembl records and manual input become
// canonical entries keyed by the left-aligned VCF triple, with the shift, the shift tract, both
// discriminating positions and the guarded zone.
//
// Pure: no network, no FASTA I/O and no config reads. Genome bases come from an injected accessor
// {bases(start, end) -> uppercase plus-strand string, regionLength}, normally one sequenceWindow() over the
// single sequence.fetch string of the operation (§3.6). Coordinates are 1-based inclusive on the reference
// plus strand. A normalized variant `v` is {position, ref, alt}: left-aligned VCF with an anchor base for
// indels (§2.1). Nothing here derives a fractional number, so the §2.1 rounding rule has no work to do.

const { PrimerHttpError } = require('../errors');

// The config defaults of primers.variation (§3.1); callers normally pass cfg.variation itself.
const DEFAULTS = Object.freeze({
  max_allele_length: 50,
  max_shift: 1000,
  ems_source_pattern: '^EMS_',
  submission_flank: 50
});

// §3.5 compares the haplotypes inside W bases of the VCF span, W doubling from 200 up to 1,600. Scanning
// outward from the variant finds the same first difference for every W that contains it, so W survives
// here only as the bound of the scan.
const MAX_W = 1600;

const ALLELE_RE = /^([ACGTN]+|-|\*)$/;
const MANUAL_ALLELE_RE = /^([ACGT]+|-)$/;
const ID_RE = /^[A-Za-z0-9][A-Za-z0-9_.:,*-]{0,254}$/;
const RS_RE = /^rs[0-9]+$/;

// §4.17 two-base IUPAC codes, keyed by the alphabetically sorted pair.
const IUPAC = Object.freeze({ AC: 'M', AG: 'R', AT: 'W', CG: 'S', CT: 'Y', GT: 'K' });

// Issue codes in PrimerVariantIssue enum order; only STAR_ALLELE leaves an entry designable (§3.4).
const ISSUE_ORDER = ['REF_MISMATCH', 'STAR_ALLELE', 'ALLELE_TOO_LONG', 'UNSUPPORTED_ALLELE', 'REPEAT_TOO_LONG'];

function assertInt(v, name) {
  if (!Number.isSafeInteger(v)) throw new TypeError(name + ' must be an integer');
}

function opt(opts, key) {
  return opts && opts[key] !== undefined && opts[key] !== null ? opts[key] : DEFAULTS[key];
}

function cmp(a, b) {
  return a < b ? -1 : a > b ? 1 : 0;
}

function unique(list) {
  return list.filter(function (x, i) { return list.indexOf(x) === i; });
}

function codedRangeError(code, message) {
  const e = new RangeError(message);
  e.code = code;
  return e;
}

function invalidVariant(reason, message) {
  return new PrimerHttpError(400, 'INVALID_VARIANT', message, { reason: reason });
}

// ---- genome access ------------------------------------------------------------------------------------

// A read outside the fetched window: the caller re-fetches a wider window and retries (§3.6).
class SequenceWindowError extends RangeError {
  constructor(start, end, win) {
    super('bases ' + start + '-' + end + ' are outside the fetched window ' + win.start + '-' + win.end);
    this.name = 'SequenceWindowError';
    this.code = 'SEQUENCE_WINDOW';
    this.needed = { start: start, end: end };
  }
}

// One fetched plus-strand string starting at `start` as a genome accessor. Reads outside the region
// throw a RangeError with code OUTSIDE_REGION; reads outside the window throw SequenceWindowError.
function sequenceWindow(seq, start, regionLength) {
  if (typeof seq !== 'string') throw new TypeError('seq must be a string');
  assertInt(start, 'start');
  if (regionLength !== undefined && regionLength !== null) assertInt(regionLength, 'regionLength');
  const s = seq.toUpperCase();
  const win = {
    start: start,
    end: start + s.length - 1,
    regionLength: regionLength === undefined || regionLength === null ? Infinity : regionLength,
    bases: function (a, b) {
      if (b < a) return '';
      if (a < 1 || b > win.regionLength) {
        throw codedRangeError('OUTSIDE_REGION', 'bases ' + a + '-' + b + ' are outside the region (1-' + win.regionLength + ')');
      }
      if (a < win.start || b > win.end) throw new SequenceWindowError(a, b, win);
      return s.slice(a - win.start, b - win.start + 1);
    }
  };
  return win;
}

function regionLengthOf(genome) {
  return Number.isSafeInteger(genome.regionLength) ? genome.regionLength : Infinity;
}

function base(genome, p) {
  return genome.bases(p, p);
}

// ---- one variant --------------------------------------------------------------------------------------

// §3.5, vt-style: trim a shared last base, re-anchor an empty allele one base left, repeat; then trim
// shared first bases while both alleles keep at least two.
function leftNormalize(vcf, genome) {
  assertInt(vcf.position, 'position');
  let pos = vcf.position;
  let ref = String(vcf.ref).toUpperCase();
  let alt = String(vcf.alt).toUpperCase();
  if (ref === alt) throw new TypeError('ref and alt must differ');
  for (;;) {
    let changed = false;
    if (ref.length && alt.length && ref[ref.length - 1] === alt[alt.length - 1]) {
      ref = ref.slice(0, -1);
      alt = alt.slice(0, -1);
      changed = true;
    }
    if (!ref.length || !alt.length) {
      pos -= 1;
      const b = base(genome, pos);
      ref = b + ref;
      alt = b + alt;
      changed = true;
    }
    if (!changed) break;
  }
  while (ref.length >= 2 && alt.length >= 2 && ref[0] === alt[0]) {
    ref = ref.slice(1);
    alt = alt.slice(1);
    pos += 1;
  }
  return { position: pos, ref: ref, alt: alt };
}

// §2.1 kinds, for normalized alleles.
function kindOf(ref, alt) {
  if (ref.length === alt.length) return ref.length === 1 ? 'snv' : 'mnv';
  if (ref[0] === alt[0] && alt.length === 1) return 'deletion';
  if (ref[0] === alt[0] && ref.length === 1) return 'insertion';
  return 'complex';
}

// How many one-base slides to the right describe the same haplotype. Only an anchored indel can slide:
// the deleted (or inserted) string rotates while its first base equals the next reference base.
// Counting stops at max_shift + 1, so a result above max_shift means VARIANT_TOO_REPETITIVE (§3.5).
function shiftOf(v, genome, opts) {
  const max = opt(opts, 'max_shift');
  const kind = kindOf(v.ref, v.alt);
  if (kind !== 'deletion' && kind !== 'insertion') return 0;
  let cur = kind === 'deletion' ? v.ref.slice(1) : v.alt.slice(1);
  const last = kind === 'deletion' ? v.position + v.ref.length - 1 : v.position;
  const regionLength = regionLengthOf(genome);
  let s = 0;
  while (s <= max) {
    const p = last + s + 1;
    if (p > regionLength || base(genome, p) !== cur[0]) break;
    cur = cur.slice(1) + cur[0];
    s += 1;
  }
  return s;
}

// §3.5 shift tract: the reference bases the indel can slide through, anchor base excluded. It serves only
// in_shift_tract (§4.12); the §5.6 core is a different span that keeps the anchor.
function tractOf(v, shift) {
  return shift > 0 ? { start: v.position + 1, end: v.position + v.ref.length - 1 + shift } : null;
}

// Reference coordinate of the ALT-haplotype base at q, where q numbers the ALT haplotype like the
// reference up to the variant (so an ALT template index t maps to q = template.start + t - 1).
// null for an inserted base.
function altToRef(v, q) {
  const lr = v.ref.length;
  const la = v.alt.length;
  if (q < v.position) return q;
  const offset = q - v.position;
  if (offset < la) {
    if (lr === la) return v.position + offset; // an SNV or MNV base
    if (offset === 0 && v.ref[0] === v.alt[0]) return v.position; // the anchor base of an indel
    return null;
  }
  return v.position + lr + (offset - la);
}

// §3.5: the 3' base of a forward allele-specific primer is the first base, left to right, where the
// haplotypes differ; a reverse primer's is the first one right to left. Both scans run outward from the
// variant, reading v.ref for the REF allele (for a verified entry that is the genome) and at most W = 1,600
// bases past the shorter allele, as a W-window comparison would.
function discriminating(v, genome) {
  const p = v.position;
  const R = v.ref;
  const A = v.alt;
  const lr = R.length;
  const la = A.length;
  const bound = Math.min(lr, la) + MAX_W;
  const regionLength = regionLengthOf(genome);
  let forward = null;
  let reverse = null;
  for (let k = 0; k < bound && !forward; k++) {
    const rp = p + k;
    const ap = p + lr + (k - la);
    if ((k >= lr && rp > regionLength) || (k >= la && ap > regionLength)) break;
    const rb = k < lr ? R[k] : base(genome, rp);
    const ab = k < la ? A[k] : base(genome, ap);
    if (rb !== ab) forward = { position: rp, ref_base: rb, alt_base: ab, alt_maps_to: altToRef(v, p + k) };
  }
  for (let m = 0; m < bound && !reverse; m++) {
    const rp = p + lr - 1 - m;
    const ap = p + la - 1 - m;
    if ((m >= lr && rp < 1) || (m >= la && ap < 1)) break;
    const rb = m < lr ? R[lr - 1 - m] : base(genome, rp);
    const ab = m < la ? A[la - 1 - m] : base(genome, ap);
    if (rb !== ab) reverse = { position: rp, ref_base: rb, alt_base: ab, alt_maps_to: altToRef(v, ap) };
  }
  if (!forward || !reverse) {
    throw codedRangeError('NO_DISCRIMINATING_BASE', 'the haplotypes of ' + p + ':' + R + ':' + A +
      ' do not differ within ' + MAX_W + ' bp or the region');
  }
  return { forward: forward, reverse: reverse };
}

// §3.5 zone: the VCF span, both discriminating positions and both ALT 3' anchor mappings.
function zoneOf(v, d) {
  const rev = d.reverse.alt_maps_to === null ? Infinity : d.reverse.alt_maps_to;
  const fwd = d.forward.alt_maps_to === null ? -Infinity : d.forward.alt_maps_to;
  return {
    start: Math.min(v.position, d.reverse.position, rev),
    end: Math.max(v.position + v.ref.length - 1, d.forward.position, fwd)
  };
}

// The reference span [start, end] read on each haplotype (§5.2, §5.6). The span must cover the VCF span;
// the two strings differ only inside it.
function haplotypes(v, genome, start, end) {
  const vEnd = v.position + v.ref.length - 1;
  if (start > v.position || end < vEnd) throw new RangeError('the span ' + start + '-' + end + ' must cover the variant');
  const left = genome.bases(start, v.position - 1);
  const right = genome.bases(vEnd + 1, end);
  return { start: start, end: end, ref: left + genome.bases(v.position, vEnd) + right, alt: left + v.alt + right };
}

// Ensembl style {start, end, ref, alt} ('-' = empty allele; an insertion has start = end + 1).
function minimalOf(v) {
  let start = v.position;
  let ref = v.ref;
  let alt = v.alt;
  while (ref.length && alt.length && ref[0] === alt[0]) {
    ref = ref.slice(1);
    alt = alt.slice(1);
    start += 1;
  }
  return { start: start, end: start + ref.length - 1, ref: ref || '-', alt: alt || '-' };
}

// '1:11109 C/A', '1:11283-11283 A/-', '1:11502^11503 -/GT'.
function labelOf(region, m) {
  if (m.ref === '-') return region + ':' + m.end + '^' + m.start + ' -/' + m.alt;
  if (m.start === m.end && m.alt.length === 1 && m.alt !== '-') return region + ':' + m.start + ' ' + m.ref + '/' + m.alt;
  return region + ':' + m.start + '-' + m.end + ' ' + m.ref + '/' + m.alt;
}

// ---- records to entries -------------------------------------------------------------------------------

// §3.4: both insertion conventions (EVA end = start, SAP/BAP end = start - 1) insert before `start`.
function recordVcf(start, refAllele, alt, genome) {
  if (refAllele === '-') {
    const p = start - 1;
    const b = base(genome, p);
    return { position: p, ref: b, alt: b + alt };
  }
  if (alt === '-') {
    const p = start - 1;
    const b = base(genome, p);
    return { position: p, ref: b + refAllele, alt: b };
  }
  return { position: start, ref: refAllele, alt: alt };
}

// The requested id first, then rs ids by number, then the rest by code unit.
function orderIds(ids, requested) {
  const rest = ids.filter(function (id) { return id !== requested; });
  const rs = rest.filter(function (id) { return RS_RE.test(id); })
    .sort(function (a, b) { return a.length - b.length || cmp(a, b); });
  const other = rest.filter(function (id) { return !RS_RE.test(id); }).sort(cmp);
  return (ids.indexOf(requested) >= 0 ? [requested] : []).concat(rs, other);
}

const emsCache = new Map();
function emsTester(pattern) {
  if (pattern instanceof RegExp) return pattern;
  if (!emsCache.has(pattern)) emsCache.set(pattern, new RegExp(pattern));
  return emsCache.get(pattern);
}

function addIssue(issues, code, message, details) {
  const json = JSON.stringify([code, details]);
  if (issues.some(function (i) { return JSON.stringify([i.code, i.details]) === json; })) return;
  issues.push({ code: code, message: message, details: details });
}

function buildEntry(group, genome, o) {
  const v = group.v;
  const region = group.region;
  const ids = orderIds(unique(group.hits.map(function (h) { return h.id; }).filter(function (id) { return id !== null; })),
    o.requested_id);
  const rank = new Map(ids.map(function (id, i) { return [id, i]; }));
  // Real records in id order, then manual pseudo-records (stable sort).
  const hits = group.hits.slice().sort(function (a, b) {
    return (a.id === null ? ids.length : rank.get(a.id)) - (b.id === null ? ids.length : rank.get(b.id));
  });

  const records = ids.map(function (id) {
    const h = hits.find(function (x) { return x.id === id; });
    const source = typeof h.record.source === 'string' ? h.record.source : null;
    return { id: id, source: source, ems: source !== null && o.ems.test(source) };
  });

  const minimal = minimalOf(v);
  const issues = [];
  hits.forEach(function (h) {
    if (h.refOk) return;
    addIssue(issues, 'REF_MISMATCH', 'the reference allele ' + h.alleles[0] + ' at ' + region + ':' + h.record.start +
      ' does not match the genome (' + h.genomeRef + ')',
    { region: region, position: h.record.start, given: h.alleles[0], genome: h.genomeRef });
  });
  const starIds = unique(hits.filter(function (h) { return h.star; }).map(function (h) { return h.id; }));
  if (starIds.length) {
    addIssue(issues, 'STAR_ALLELE', 'the site also has a spanning-deletion allele (*), which is never designed',
      { ids: starIds });
  }
  hits.forEach(function (h) {
    const long = [h.alleles[0], h.alt].filter(function (a) { return a.length > o.max_allele_length; });
    if (long.length) {
      addIssue(issues, 'ALLELE_TOO_LONG', 'an allele is longer than ' + o.max_allele_length + ' nt',
        { length: Math.max.apply(null, long.map(function (a) { return a.length; })), max: o.max_allele_length });
    }
    [h.alleles[0], h.alt].forEach(function (a) {
      if (a.indexOf('N') >= 0) addIssue(issues, 'UNSUPPORTED_ALLELE', 'allele ' + a + ' contains N', { allele: a });
    });
  });

  const shift = shiftOf(v, genome, o);
  let disc = null;
  let zone = null;
  if (shift > o.max_shift) {
    addIssue(issues, 'REPEAT_TOO_LONG', 'the event can slide more than ' + o.max_shift + ' bp, so no stable discriminating base exists',
      { region: region, position: v.position, max: o.max_shift });
  } else {
    try {
      disc = discriminating(v, genome);
      zone = zoneOf(v, disc);
    } catch (e) {
      if (e.code !== 'NO_DISCRIMINATING_BASE') throw e;
      addIssue(issues, 'REPEAT_TOO_LONG', 'the haplotypes do not differ within ' + MAX_W + ' bp or the region',
        { region: region, position: v.position, max: o.max_shift });
    }
  }
  issues.sort(function (a, b) { return ISSUE_ORDER.indexOf(a.code) - ISSUE_ORDER.indexOf(b.code); });

  // The site's allele list comes from the merged record with the most alleles (real records first).
  let site = hits[0];
  hits.forEach(function (h) { if (h.alleles.length > site.alleles.length) site = h; });
  const otherAlts = unique(site.alleles.slice(1).filter(function (a) { return a !== site.alt; }));
  const consequenceHit = hits.find(function (h) {
    return typeof h.record.consequence_type === 'string' && h.record.consequence_type !== '';
  });
  const synonyms = o.requested_id !== undefined && ids.indexOf(o.requested_id) >= 0
    ? unique((o.synonyms || []).filter(function (s) { return typeof s === 'string' && ID_RE.test(s) && ids.indexOf(s) < 0; }))
    : [];

  return {
    key: group.key,
    ids: ids,
    synonyms: synonyms,
    label: labelOf(region, minimal),
    kind: kindOf(v.ref, v.alt),
    region: region,
    vcf: { position: v.position, ref: v.ref, alt: v.alt },
    minimal: minimal,
    // A manual pseudo-record reports its alleles in the minimal form, so both input styles give one entry.
    alleles: site.id === null ? [minimal.ref, minimal.alt] : site.alleles.slice(),
    multiallelic: otherAlts.length ? { alleles: site.alleles.slice(), other_alts: otherAlts } : null,
    shift: shift,
    zone: zone,
    discriminating: disc,
    records: records,
    ems: records.length > 0 && records.every(function (r) { return r.ems; }),
    consequence: consequenceHit ? consequenceHit.record.consequence_type : null,
    ref_verified: hits.every(function (h) { return h.refOk; }),
    designable: issues.every(function (i) { return i.code === 'STAR_ALLELE'; }),
    issues: issues
  };
}

// §3.4. `records` are validated compact overlap records {id, seq_region_name, start, end, alleles, source,
// consequence_type} or parseManual pseudo-records (id null). All must lie on one region, which `genome`
// reads. opts: the primers.variation keys (max_allele_length, max_shift, ems_source_pattern), plus
// `region` (default: each record's seq_region_name), `requested_id` (ordered first in ids) and `synonyms`
// (the raw lookup list, filtered by the id pattern into the entries holding requested_id).
// Returns {entries sorted by vcf.position then key, skipped: {count, reasons}}. A read outside the fetched
// window throws SequenceWindowError; the caller widens the window and calls again.
function recordsToEntries(records, genome, opts) {
  const o = {
    max_allele_length: opt(opts, 'max_allele_length'),
    max_shift: opt(opts, 'max_shift'),
    ems: emsTester(opt(opts, 'ems_source_pattern')),
    region: opts && opts.region !== undefined ? String(opts.region) : undefined,
    requested_id: opts && typeof opts.requested_id === 'string' ? opts.requested_id : undefined,
    synonyms: opts && Array.isArray(opts.synonyms) ? opts.synonyms : []
  };
  const skipped = { count: 0, reasons: {} };
  function skip(reason) {
    skipped.count += 1;
    skipped.reasons[reason] = (skipped.reasons[reason] || 0) + 1;
  }
  function edge(e) {
    if (e.code !== 'OUTSIDE_REGION') throw e;
    skip('region_edge');
  }

  const groups = new Map();
  const seen = new Set();
  (records || []).forEach(function (r) {
    const alleles = r && Array.isArray(r.alleles) && r.alleles.every(function (a) { return typeof a === 'string'; })
      ? r.alleles.map(function (a) { return a.toUpperCase(); }) : null;
    if (!alleles || alleles.length < 2 || !alleles.every(function (a) { return ALLELE_RE.test(a); })) return skip('invalid_alleles');
    if (!Number.isSafeInteger(r.start) || !Number.isSafeInteger(r.end)) return skip('invalid_coordinates');
    const region = o.region !== undefined ? o.region : r.seq_region_name;
    if (typeof region !== 'string' || region === '') return skip('invalid_region');
    const id = typeof r.id === 'string' ? r.id : null;
    const signature = JSON.stringify([id, r.start, r.end, alleles]);
    if (seen.has(signature)) return; // the same record from two overlapping chunks
    seen.add(signature);
    const refAllele = alleles[0];
    if (refAllele === '*') return skip('star_reference');
    let genomeRef = null;
    try {
      genomeRef = refAllele === '-' ? null : genome.bases(r.start, r.end);
    } catch (e) {
      return edge(e);
    }
    const refOk = refAllele === '-' || genomeRef === refAllele;
    const star = alleles.indexOf('*', 1) > 0;
    unique(alleles.slice(1)).forEach(function (alt) {
      if (alt === '*') return;
      if (alt === refAllele) return skip('identical_alleles');
      let v;
      try {
        v = leftNormalize(recordVcf(r.start, refAllele, alt, genome), genome);
      } catch (e) {
        return edge(e);
      }
      const key = region + ':' + v.position + ':' + v.ref + ':' + v.alt;
      let group = groups.get(key);
      if (!group) {
        group = { key: key, region: region, v: v, hits: [] };
        groups.set(key, group);
      }
      group.hits.push({ record: r, id: id, alleles: alleles, alt: alt, refOk: refOk, genomeRef: genomeRef, star: star });
    });
  });

  const entries = Array.from(groups.values()).map(function (g) { return buildEntry(g, genome, o); });
  entries.sort(function (a, b) { return a.vcf.position - b.vcf.position || cmp(a.key, b.key); });
  return { entries: entries, skipped: skipped };
}

// §2.3: an entry is in [start, end] when its minimal span overlaps it; an insertion when either flanking
// coordinate is inside.
function inWindow(minimal, start, end) {
  if (minimal.start > minimal.end) {
    return (minimal.end >= start && minimal.end <= end) || (minimal.start >= start && minimal.start <= end);
  }
  return minimal.start <= end && minimal.end >= start;
}

// ---- manual input -------------------------------------------------------------------------------------

// §3.9 and §2.7 rule 3: a manual {region, position, ref, alt} becomes a pseudo-record for
// recordsToEntries. VCF style (neither allele '-'): position is the first ref base. Ensembl style
// (exactly one '-'): a deletion starts at the first deleted base; an insertion's position is the base after
// the insertion point, so end = start - 1.
function parseManual(input, opts) {
  const max = opt(opts, 'max_allele_length');
  const v = input || {};
  if (typeof v.region !== 'string' || v.region === '' || !Number.isSafeInteger(v.position) || v.position < 1 ||
      typeof v.ref !== 'string' || typeof v.alt !== 'string') {
    throw invalidVariant('id_or_manual', 'a manual variant needs region, position, ref and alt');
  }
  const ref = v.ref.toUpperCase();
  const alt = v.alt.toUpperCase();
  if (ref === alt) throw invalidVariant('alleles', 'ref and alt must be different alleles, and at most one may be "-"');
  [ref, alt].forEach(function (a) {
    if (a === '*' || /^[ACGTN]+$/.test(a) && a.indexOf('N') >= 0) {
      throw new PrimerHttpError(400, 'UNSUPPORTED_ALLELE', 'allele ' + a + ' cannot be designed', { allele: a });
    }
  });
  if (!MANUAL_ALLELE_RE.test(ref) || !MANUAL_ALLELE_RE.test(alt)) {
    throw invalidVariant('alleles', 'alleles must be A/C/G/T strings or "-"');
  }
  if (ref.length > max || alt.length > max) {
    throw invalidVariant('allele_too_long', 'alleles are limited to ' + max + ' nt');
  }
  const refLength = ref === '-' ? 0 : ref.length;
  return {
    id: null,
    source: null,
    seq_region_name: v.region,
    start: v.position,
    end: v.position + refLength - 1,
    alleles: [ref, alt]
  };
}

// ---- submission string --------------------------------------------------------------------------------

// §4.17: `submission_flank` (50) reference bases each side of the minimal alleles, '[REF/ALT]' with '-'
// written as ''. Biallelic non-EMS SNV neighbours become IUPAC codes; a '*' allele does not make a site
// multi-allelic. Other non-EMS neighbours touching a flank (indels, MNVs, multi-allelic sites, SNVs whose REF
// was not verified) stay reference bases and are listed in omitted_ids for SUBMISSION_NEIGHBOURS_OMITTED.
// EMS neighbours are ignored. `neighbours` are recordsToEntries entries; the target's own key is skipped.
function submissionSequence(entry, genome, neighbours, opts) {
  const flank = opt(opts, 'submission_flank');
  const m = entry.minimal;
  const left = { start: Math.max(1, m.start - flank), end: m.start - 1 };
  const right = { start: m.end + 1, end: Math.min(regionLengthOf(genome), m.end + flank) };
  const leftBases = genome.bases(left.start, left.end).split('');
  const rightBases = genome.bases(right.start, right.end).split('');
  const snvs = new Map();
  const omitted = [];

  (neighbours || []).forEach(function (n) {
    if (!n || n.key === entry.key || n.ems) return;
    const lo = Math.min(n.minimal.start, n.minimal.end);
    const hi = Math.max(n.minimal.start, n.minimal.end);
    if (!(lo <= left.end && hi >= left.start) && !(lo <= right.end && hi >= right.start)) return;
    const multi = n.multiallelic !== null && n.multiallelic !== undefined &&
      n.multiallelic.other_alts.some(function (a) { return a !== '*'; });
    if (n.kind === 'snv' && n.ref_verified !== false && !multi) {
      if (!snvs.has(lo)) snvs.set(lo, []);
      snvs.get(lo).push(n);
    } else {
      omitted.push(n);
    }
  });
  snvs.forEach(function (list, pos) {
    const alts = unique(list.map(function (n) { return n.vcf.alt; }));
    const code = alts.length === 1 ? IUPAC[[list[0].vcf.ref, alts[0]].sort().join('')] : undefined;
    if (!code) return Array.prototype.push.apply(omitted, list);
    if (pos <= left.end) leftBases[pos - left.start] = code;
    else rightBases[pos - right.start] = code;
  });

  omitted.sort(function (a, b) { return a.vcf.position - b.vcf.position || cmp(a.key, b.key); });
  const ids = unique([].concat.apply([], omitted.map(function (n) { return n.ids; })));
  const allele = function (a) { return a === '-' ? '' : a; };
  return {
    sequence: leftBases.join('') + '[' + allele(m.ref) + '/' + allele(m.alt) + ']' + rightBases.join(''),
    omitted_ids: ids
  };
}

module.exports = {
  recordsToEntries,
  leftNormalize,
  shiftOf,
  tractOf,
  discriminating,
  zoneOf,
  altToRef,
  haplotypes,
  kindOf,
  minimalOf,
  labelOf,
  parseManual,
  submissionSequence,
  IUPAC,
  inWindow,
  sequenceWindow,
  SequenceWindowError,
  DEFAULTS
};
