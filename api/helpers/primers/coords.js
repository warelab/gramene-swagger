'use strict';

// Pure coordinate helpers (spec §A.3). Template coordinates are 1-based inclusive;
// intervals are [start, length]; genomic ranges are {start, end} with start <= end.

function assertInt(v, name) {
  if (!Number.isSafeInteger(v)) throw new TypeError(name + ' must be an integer');
}

function assertStrand(strand) {
  if (strand !== 1 && strand !== -1) throw new TypeError('strand must be 1 or -1');
}

function assertSide(side) {
  if (side !== 'left' && side !== 'right') throw new TypeError('side must be "left" or "right"');
}

// ---- footprints ---------------------------------------------------------------------------

// PRIMER_LEFT_i=pos,len -> [pos, pos+len-1]
function leftFootprint(pos, len) {
  assertInt(pos, 'pos');
  assertInt(len, 'len');
  return { start: pos, end: pos + len - 1 };
}

// PRIMER_RIGHT_i=pos,len names the rightmost template base -> [pos-len+1, pos]
function rightFootprint(pos, len) {
  assertInt(pos, 'pos');
  assertInt(len, 'len');
  return { start: pos - len + 1, end: pos };
}

function footprint(side, pos, len) {
  assertSide(side);
  return side === 'left' ? leftFootprint(pos, len) : rightFootprint(pos, len);
}

function productSize(leftPos, rightPos) {
  return rightPos - leftPos + 1;
}

function intervalToRange(interval) {
  if (!Array.isArray(interval) || interval.length !== 2) throw new TypeError('interval must be [start, length]');
  assertInt(interval[0], 'interval start');
  assertInt(interval[1], 'interval length');
  return { start: interval[0], end: interval[0] + interval[1] - 1 };
}

function rangeToInterval(start, end) {
  assertInt(start, 'start');
  assertInt(end, 'end');
  return [start, end - start + 1];
}

// Genomic strand whose 5'->3' sequence equals the primer: left = template strand, right = opposite.
function primerStrand(side, templateStrand) {
  assertSide(side);
  assertStrand(templateStrand);
  return side === 'left' ? templateStrand : -templateStrand;
}

// ---- gene-relative positions --------------------------------------------------------------

// Gene-relative 1-based position p (transcription order) -> genomic coordinate.
function geneRelativeToGenomic(loc, p) {
  assertStrand(loc.strand);
  assertInt(p, 'p');
  return loc.strand === 1 ? loc.start + p - 1 : loc.end - p + 1;
}

// Gene-relative exon {start, end} -> ascending genomic {start, end}.
function exonGenomicRange(loc, exon) {
  assertStrand(loc.strand);
  return loc.strand === 1
    ? { start: loc.start + exon.start - 1, end: loc.start + exon.end - 1 }
    : { start: loc.end - exon.end + 1, end: loc.end - exon.start + 1 };
}

// Exons in transcript order -> contiguous cDNA segments {id, t_start, t_end, g_start, g_end}.
function buildSegments(loc, exons) {
  if (!Array.isArray(exons) || exons.length === 0) throw new TypeError('exons must be a non-empty array');
  let cum = 0;
  return exons.map(function (e) {
    assertInt(e.start, 'exon start');
    assertInt(e.end, 'exon end');
    if (e.end < e.start) throw new TypeError('exon end < start');
    const n = e.end - e.start + 1;
    const g = exonGenomicRange(loc, e);
    const seg = { id: e.id, t_start: cum + 1, t_end: cum + n, g_start: g.start, g_end: g.end };
    cum += n;
    return seg;
  });
}

// Junction p = boundary between cDNA bases p and p+1: cumulative segment ends except the last.
function junctionsFromSegments(segments) {
  return segments.slice(0, -1).map(function (s) { return s.t_end; });
}

// ---- mappers ------------------------------------------------------------------------------

// Genomic template (gene/region mode): template spans genomic [start, end] on `strand`.
// t -> strand == 1 ? start + t - 1 : end - t + 1
function genomicMapper(opts) {
  const region = opts.region;
  const start = opts.start;
  const end = opts.end;
  const strand = opts.strand;
  assertInt(start, 'start');
  assertInt(end, 'end');
  assertStrand(strand);
  if (end < start) throw new TypeError('end < start');
  const length = end - start + 1;

  function toGenomic(t) {
    assertInt(t, 't');
    if (t < 1 || t > length) throw new RangeError('template position ' + t + ' outside 1..' + length);
    return strand === 1 ? start + t - 1 : end - t + 1;
  }
  function toTemplate(g) {
    assertInt(g, 'g');
    if (g < start || g > end) return null;
    return strand === 1 ? g - start + 1 : end - g + 1;
  }
  function blocks(a, b) {
    const x = toGenomic(a);
    const y = toGenomic(b);
    return [{ start: Math.min(x, y), end: Math.max(x, y) }];
  }
  function primer(side, a, b) {
    const bl = blocks(a, b);
    return { region: region, start: bl[0].start, end: bl[0].end, strand: primerStrand(side, strand), blocks: bl };
  }
  function product(a, b) {
    const bl = blocks(a, b);
    return { region: region, start: bl[0].start, end: bl[0].end, strand: strand };
  }
  return { kind: 'genomic', region: region, start: start, end: end, strand: strand, length: length,
    toGenomic: toGenomic, toTemplate: toTemplate, blocks: blocks, primer: primer, product: product };
}

// Spliced (transcript) template: segments from buildSegments(), in cDNA order.
// A span [a, b] is clipped to each segment and each piece mapped to an ascending genomic block:
//   strand +1: [g_start + lo - t_start, g_start + hi - t_start]
//   strand -1: [g_end - (hi - t_start), g_end - (lo - t_start)]
function splicedMapper(opts) {
  const region = opts.region;
  const strand = opts.strand;
  const segments = opts.segments;
  assertStrand(strand);
  if (!Array.isArray(segments) || segments.length === 0) throw new TypeError('segments must be a non-empty array');
  segments.forEach(function (s, i) {
    const expectStart = i === 0 ? 1 : segments[i - 1].t_end + 1;
    if (s.t_start !== expectStart || s.t_end < s.t_start || s.g_end - s.g_start !== s.t_end - s.t_start) {
      throw new TypeError('segments must be contiguous in cDNA and match their genomic lengths (segment ' + i + ')');
    }
  });
  const length = segments[segments.length - 1].t_end;

  function blocks(a, b) {
    assertInt(a, 'a');
    assertInt(b, 'b');
    if (a < 1 || b > length || a > b) throw new RangeError('template span ' + a + '-' + b + ' outside 1..' + length);
    const out = [];
    for (const s of segments) {
      const lo = Math.max(a, s.t_start);
      const hi = Math.min(b, s.t_end);
      if (lo > hi) continue;
      out.push(strand === 1
        ? { start: s.g_start + lo - s.t_start, end: s.g_start + hi - s.t_start }
        : { start: s.g_end - (hi - s.t_start), end: s.g_end - (lo - s.t_start) });
    }
    return out.sort(function (x, y) { return x.start - y.start; });
  }
  function toGenomic(t) {
    return blocks(t, t)[0].start;
  }
  function toTemplate(g) {
    assertInt(g, 'g');
    for (const s of segments) {
      if (g >= s.g_start && g <= s.g_end) return strand === 1 ? s.t_start + (g - s.g_start) : s.t_start + (s.g_end - g);
    }
    return null;
  }
  function primer(side, a, b) {
    const bl = blocks(a, b);
    return { region: region, start: bl[0].start, end: bl[bl.length - 1].end, strand: primerStrand(side, strand), blocks: bl };
  }
  // Envelope of the product's genomic footprint (= envelope of both primers' blocks).
  function product(a, b) {
    const x = toGenomic(a);
    const y = toGenomic(b);
    return { region: region, start: Math.min(x, y), end: Math.max(x, y), strand: strand };
  }
  function genomicSize(a, b) {
    const p = product(a, b);
    return p.end - p.start + 1;
  }
  return { kind: 'spliced', region: region, strand: strand, length: length, segments: segments,
    junctions: junctionsFromSegments(segments),
    toGenomic: toGenomic, toTemplate: toTemplate, blocks: blocks, primer: primer, product: product, genomicSize: genomicSize };
}

// ---- junction overlap ---------------------------------------------------------------------

// Left primer [a,b] spans j when j-a+1 >= min5 and b-j >= min3.
// Right primer [a,b] (5' end at b) spans j when j-a+1 >= min3 and b-j >= min5.
function spansJunction(side, a, b, j, min5, min3) {
  assertSide(side);
  return side === 'left'
    ? (j - a + 1 >= min5 && b - j >= min3)
    : (j - a + 1 >= min3 && b - j >= min5);
}

// First junction spanned by the primer footprint [a,b], as
// {position, overlap_5p, overlap_3p}, or null. Defaults are Primer3's (min5 7, min3 4).
function junctionOverlap(side, a, b, junctions, opts) {
  opts = opts || {};
  const min5 = opts.min5 === undefined ? 7 : opts.min5;
  const min3 = opts.min3 === undefined ? 4 : opts.min3;
  assertSide(side);
  for (const j of junctions || []) {
    if (j < a || j >= b) continue;
    if (spansJunction(side, a, b, j, min5, min3)) {
      return side === 'left'
        ? { position: j, overlap_5p: j - a + 1, overlap_3p: b - j }
        : { position: j, overlap_5p: b - j, overlap_3p: j - a + 1 };
    }
  }
  return null;
}

// ---- mask runs ([start, length], 1-based) ------------------------------------------------

// Sort, drop empty runs, merge overlapping or adjacent runs; optionally clip to 1..length.
function mergeRuns(runs, opts) {
  const limit = opts && opts.length;
  const ranges = [];
  (runs || []).forEach(function (r) {
    let s = r[0];
    let e = r[0] + r[1] - 1;
    if (limit !== undefined) {
      s = Math.max(1, s);
      e = Math.min(limit, e);
    }
    if (Number.isFinite(s) && Number.isFinite(e) && e >= s) ranges.push([s, e]);
  });
  ranges.sort(function (x, y) { return x[0] - y[0] || x[1] - y[1]; });
  const out = [];
  for (const r of ranges) {
    const last = out[out.length - 1];
    if (last && r[0] <= last[1] + 1) last[1] = Math.max(last[1], r[1]);
    else out.push([r[0], r[1]]);
  }
  return out.map(function (r) { return [r[0], r[1] - r[0] + 1]; });
}

// Truthy flags (index 0 = base 1) -> merged runs.
function runsFromFlags(flags) {
  const out = [];
  let runStart = 0;
  const n = flags.length;
  for (let i = 0; i <= n; i++) {
    const on = i < n && !!flags[i];
    if (on && !runStart) runStart = i + 1;
    else if (!on && runStart) {
      out.push([runStart, i + 1 - runStart]);
      runStart = 0;
    }
  }
  return out;
}

// Runs of lowercase letters (soft-masked bases).
function lowercaseRuns(seq) {
  const out = [];
  const re = /[a-z]+/g;
  let m;
  while ((m = re.exec(seq)) !== null) out.push([m.index + 1, m[0].length]);
  return out;
}

function maskedBases(runs) {
  return mergeRuns(runs).reduce(function (sum, r) { return sum + r[1]; }, 0);
}

function maskedFraction(runs, length) {
  if (!length) return 0;
  return maskedBases(mergeRuns(runs, { length: length })) / length;
}

// True when [start, end] overlaps any run.
function overlapsRuns(start, end, runs) {
  for (const r of runs || []) {
    if (r[0] <= end && r[0] + r[1] - 1 >= start) return true;
  }
  return false;
}

module.exports = {
  leftFootprint,
  rightFootprint,
  footprint,
  productSize,
  intervalToRange,
  rangeToInterval,
  primerStrand,
  geneRelativeToGenomic,
  exonGenomicRange,
  buildSegments,
  junctionsFromSegments,
  genomicMapper,
  splicedMapper,
  spansJunction,
  junctionOverlap,
  mergeRuns,
  runsFromFlags,
  lowercaseRuns,
  maskedBases,
  maskedFraction,
  overlapsRuns
};
