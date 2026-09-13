'use strict';

// Template construction for POST /primers/design (spec §A.4.3-§A.4.6).
//
// buildTemplate(req, deps) turns a normalized design request (design.normalize) into the
// sequence Primer3 designs on, plus the features and the coordinate mapper used to annotate the
// returned primers. Template coordinates are 1-based inclusive in template orientation:
// transcription orientation for gene and transcript modes, the requested strand in region mode.

const coords = require('./coords');
const { PrimerHttpError } = require('./errors');

const GENE_FIELDS = Object.freeze({ _id: 1, name: 1, system_name: 1, taxon_id: 1, location: 1, gene_structure: 1 });
const MAX_ID_LENGTH = 255;
const MIN_SEQUENCE_LENGTH = 20;
// Exon (or mask) windows are read with one fetch of their genomic envelope when the envelope is
// at most this long; longer genes are read window by window.
const SPAN_FETCH_LIMIT = 1000000;
const SEQUENCE_ALPHABET_RE = /^[ACGTRYKMSWBDHVN]+$/i;
const INVALID_BASE_RE = /[^ACGTRYKMSWBDHVN]/i;
const IUPAC_RE = /[RYKMSWBDHV]/gi;
const MONGO_RETRY_AFTER_S = 30;

function logSafe(log, level, msg) {
  try {
    const l = log || console;
    (l[level] || l.error || l.log).call(l, msg);
  } catch (e) { /* logging must never throw */ }
}

function isPosInt(v) {
  return Number.isSafeInteger(v) && v >= 1;
}

function configOf(deps) {
  return deps.cfg || require('./config').get();
}

function sequenceOf(deps) {
  return deps.sequence || require('./sequence');
}

function warning(code, message) {
  return { code: code, message: message };
}

function invalidRequest(message, details) {
  return new PrimerHttpError(400, 'INVALID_REQUEST', message, details || {});
}

function tooLong(length, max, hint) {
  const details = { length: length, max: max };
  if (hint) details.hint = hint;
  return new PrimerHttpError(400, 'TEMPLATE_TOO_LONG',
    'the template would be ' + length + ' bp; the limit is ' + max + ' bp' + (hint ? ' (' + hint + ')' : ''), details);
}

function structureError(geneId, reason, log) {
  logSafe(log, 'error', 'primers: inconsistent gene structure for ' + geneId + ': ' + reason);
  return new PrimerHttpError(500, 'GENE_STRUCTURE_MISMATCH', 'the stored gene structure of ' + geneId + ' is inconsistent',
    { gene_id: String(geneId) });
}

// ---- gene documents ------------------------------------------------------------------------

// genes.findOne({_id: gene_id}) with a string-equality lookup only. deps.findGene(id) replaces mongo.
async function findGene(geneId, deps) {
  deps = deps || {};
  if (typeof geneId !== 'string' || geneId === '' || geneId.length > MAX_ID_LENGTH) {
    throw invalidRequest('gene_id must be a non-empty string of at most ' + MAX_ID_LENGTH + ' characters', { field: 'gene_id' });
  }
  if (typeof deps.findGene === 'function') return (await deps.findGene(geneId)) || null;
  const mongo = deps.mongo || require('gramene-mongodb-config');
  let doc;
  try {
    // gramene-mongodb-config resolves mongoCollection() to undefined when the connection failed.
    const genes = await mongo.genes.mongoCollection();
    if (!genes) throw new Error('the genes collection is unavailable');
    doc = await genes.findOne({ _id: geneId }, { fields: Object.assign({}, GENE_FIELDS) });
  } catch (err) {
    logSafe(deps.log, 'error', 'primers: gene lookup failed: ' + (err && err.message ? err.message : String(err)));
    throw new PrimerHttpError(503, 'MONGO_UNAVAILABLE', 'the gene database is temporarily unavailable',
      { retry_after_s: MONGO_RETRY_AFTER_S });
  }
  return doc || null;
}

async function loadGene(req, deps) {
  const doc = await findGene(req.gene_id, deps);
  if (!doc || typeof doc !== 'object') {
    throw new PrimerHttpError(404, 'UNKNOWN_GENE', 'gene ' + JSON.stringify(req.gene_id) + ' was not found', { gene_id: req.gene_id });
  }
  if (typeof doc.system_name !== 'string' || doc.system_name === '') {
    throw structureError(req.gene_id, 'system_name is missing', deps.log);
  }
  if (req.system_name !== undefined && req.system_name !== null && req.system_name !== doc.system_name) {
    throw new PrimerHttpError(400, 'SYSTEM_NAME_MISMATCH',
      'gene ' + req.gene_id + ' belongs to ' + doc.system_name + ', not ' + JSON.stringify(String(req.system_name)),
      { gene_id: req.gene_id, system_name: String(req.system_name), gene_system_name: doc.system_name });
  }
  return doc;
}

function geneLocation(doc, log) {
  const loc = doc.location;
  if (!loc || typeof loc.region !== 'string' || loc.region === '' || !isPosInt(loc.start) || !isPosInt(loc.end) ||
      loc.end < loc.start || (loc.strand !== 1 && loc.strand !== -1)) {
    throw structureError(doc._id, 'location is missing or malformed', log);
  }
  return { region: loc.region, start: loc.start, end: loc.end, strand: loc.strand };
}

// transcript_id given -> that transcript or 404 UNKNOWN_TRANSCRIPT; else the canonical transcript,
// else the first one, else null.
function pickTranscript(doc, transcriptId) {
  const gs = doc.gene_structure || {};
  const list = (Array.isArray(gs.transcripts) ? gs.transcripts : []).filter(function (t) {
    return t && typeof t.id === 'string';
  });
  if (transcriptId !== undefined && transcriptId !== null) {
    const found = typeof transcriptId === 'string' ? list.find(function (t) { return t.id === transcriptId; }) : undefined;
    if (!found) {
      throw new PrimerHttpError(404, 'UNKNOWN_TRANSCRIPT',
        'transcript ' + JSON.stringify(String(transcriptId)) + ' is not a transcript of gene ' + doc._id,
        { gene_id: String(doc._id), transcript_id: String(transcriptId) });
    }
    return found;
  }
  return list.find(function (t) { return t.id === gs.canonical_transcript; }) || list[0] || null;
}

// transcript.exons (exon ids, transcription order) -> [{id, start, end}] gene-relative, validated.
function transcriptExons(doc, transcript, geneLength, log) {
  const gs = doc.gene_structure || {};
  const byId = new Map();
  (Array.isArray(gs.exons) ? gs.exons : []).forEach(function (e) {
    if (e && typeof e.id === 'string') byId.set(e.id, e);
  });
  const refs = Array.isArray(transcript.exons) ? transcript.exons : [];
  if (refs.length === 0) throw structureError(doc._id, 'transcript ' + transcript.id + ' has no exons', log);
  return refs.map(function (ref) {
    const e = typeof ref === 'string' ? byId.get(ref) : ref;
    if (!e || !isPosInt(e.start) || !isPosInt(e.end) || e.end < e.start || e.end > geneLength) {
      throw structureError(doc._id, 'exon ' + JSON.stringify(typeof ref === 'string' ? ref : ref && ref.id) +
        ' of ' + transcript.id + ' is missing or outside the gene', log);
    }
    return { id: typeof e.id === 'string' ? e.id : String(ref), start: e.start, end: e.end };
  });
}

// cDNA position -> gene-relative position through exons in transcription order (null if outside).
function cdnaToGeneRelative(exons, c) {
  let cum = 0;
  for (const e of exons) {
    const n = e.end - e.start + 1;
    if (c <= cum + n) return e.start + (c - cum) - 1;
    cum += n;
  }
  return null;
}

function validCds(cds, cdnaLength) {
  return !!cds && isPosInt(cds.start) && isPosInt(cds.end) && cds.start <= cds.end && cds.end <= cdnaLength;
}

function exonsLength(exons) {
  return exons.reduce(function (sum, e) { return sum + e.end - e.start + 1; }, 0);
}

function sameNumbers(a, b) {
  return a.length === b.length && a.every(function (v, i) { return v === b[i]; });
}

// ---- assemblies and sequence access ---------------------------------------------------------

async function resolveAssembly(systemName, deps) {
  if (typeof deps.resolve === 'function') return deps.resolve(systemName);
  const opts = {};
  ['cfg', 'mongo', 'catalog', 'log'].forEach(function (k) {
    if (deps[k] !== undefined) opts[k] = deps[k];
  });
  return require('./assemblies').resolve(systemName, opts);
}

function dnaPath(resolved) {
  const p = resolved && resolved.fasta && resolved.fasta.dna;
  if (!p) {
    const name = resolved && resolved.system_name;
    throw new PrimerHttpError(422, 'NO_SEQUENCE', 'no genome sequence is available for ' + name, { system_name: name });
  }
  return p;
}

function assemblyWarnings(resolved) {
  return (resolved && Array.isArray(resolved.warnings) ? resolved.warnings : []).map(function (w) {
    return warning(w.code, w.message);
  });
}

async function requireRegionLength(sequence, fastaPath, region) {
  const len = await sequence.regionLength(fastaPath, region);
  if (!isPosInt(len)) {
    throw new PrimerHttpError(404, 'UNKNOWN_REGION', 'region ' + JSON.stringify(String(region)) + ' is not in this assembly',
      { region: String(region) });
  }
  return len;
}

// Sequences of genomic windows [{start, end}] (ascending coordinates) of one region, each oriented
// on `strand` and returned in the order of `windows`. Reads the envelope once when it is short.
// opts: {maxLength (sequence.fetch cap), spanFetchLimit}
async function fetchWindows(sequence, fastaPath, region, strand, windows, opts) {
  opts = opts || {};
  if (!Array.isArray(windows) || windows.length === 0) return [];
  const limit = Math.min(opts.spanFetchLimit || SPAN_FETCH_LIMIT, opts.maxLength || Infinity);
  const fetchOpts = opts.maxLength ? { maxLength: opts.maxLength } : undefined;
  let lo = Infinity;
  let hi = -Infinity;
  windows.forEach(function (w) {
    if (w.start < lo) lo = w.start;
    if (w.end > hi) hi = w.end;
  });
  if (hi - lo + 1 <= limit) {
    const env = await sequence.fetch(fastaPath, region, lo, hi, strand, fetchOpts);
    return windows.map(function (w) {
      return strand === 1 ? env.slice(w.start - lo, w.end - lo + 1) : env.slice(hi - w.end, hi - w.start + 1);
    });
  }
  const out = [];
  for (const w of windows) out.push(await sequence.fetch(fastaPath, region, w.start, w.end, strand, fetchOpts));
  return out;
}

// Gene span +- flanks, clamped to the region (spec §A.4.3 step 3).
function flankExtent(loc, regionLength, up, down) {
  if (loc.strand === 1) {
    const start = Math.max(1, loc.start - up);
    return { start: start, end: Math.min(regionLength, loc.end + down), effUp: loc.start - start };
  }
  const end = Math.min(regionLength, loc.end + up);
  return { start: Math.max(1, loc.start - down), end: end, effUp: end - loc.end };
}

// ---- sequence mode ------------------------------------------------------------------------------

// Drop FASTA header lines, remove whitespace and digits, validate IUPAC DNA, convert ambiguity codes
// to N. opts: {keepCase (lowercase = the user's own mask), maxLength}. -> {seq, iupac, name, records}
// records: how many FASTA records with sequence were joined (headers with no sequence after them do not count).
function cleanSequence(raw, opts) {
  opts = opts || {};
  if (typeof raw !== 'string') throw new PrimerHttpError(400, 'INVALID_SEQUENCE', 'sequence must be a string', {});
  let name = null;
  let records = 0;
  let recordHasSequence = false;
  const body = raw.split(/\r\n|\r|\n/).filter(function (line) {
    const m = /^\s*>\s*(\S*)/.exec(line);
    if (m) {
      if (name === null && m[1]) name = m[1];
      recordHasSequence = false;
      return false;
    }
    if (!recordHasSequence && /[^\s\d]/.test(line)) {
      recordHasSequence = true;
      records++;
    }
    return true;
  }).join('').replace(/[\s\d]/g, '');
  if (body === '') throw new PrimerHttpError(400, 'INVALID_SEQUENCE', 'the sequence is empty', {});
  if (!SEQUENCE_ALPHABET_RE.test(body)) {
    const m = INVALID_BASE_RE.exec(body);
    throw new PrimerHttpError(400, 'INVALID_SEQUENCE',
      'the sequence contains ' + JSON.stringify(m[0]) + ' at position ' + (m.index + 1) + '; only IUPAC DNA letters are allowed',
      { position: m.index + 1, character: m[0] });
  }
  if (body.length < MIN_SEQUENCE_LENGTH) {
    throw new PrimerHttpError(400, 'INVALID_SEQUENCE', 'the sequence must be at least ' + MIN_SEQUENCE_LENGTH + ' nt long',
      { length: body.length, min: MIN_SEQUENCE_LENGTH });
  }
  if (opts.maxLength && body.length > opts.maxLength) throw tooLong(body.length, opts.maxLength);
  let iupac = 0;
  const seq = (opts.keepCase ? body : body.toUpperCase()).replace(IUPAC_RE, function (c) {
    iupac++;
    return c === c.toLowerCase() ? 'n' : 'N';
  });
  return { seq: seq, iupac: iupac, name: name, records: records };
}

// ---- builders -----------------------------------------------------------------------------------

function finishTemplate(t) {
  return Object.assign({
    mode: null,
    id: 'template',
    system_name: null,
    gene_id: null,
    transcript_id: null,
    region: null,
    start: null,
    end: null,
    strand: null,
    seq: '',
    features: {},
    junctions: [],
    exon_count: null,
    mapper: null,
    resolved: null,
    location: null,
    segments: null,
    warnings: []
  }, t, { length: t.seq.length });
}

async function buildGeneTemplate(req, deps) {
  const cfg = configOf(deps);
  const log = deps.log || console;
  const sequence = sequenceOf(deps);
  const doc = await loadGene(req, deps);
  const loc = geneLocation(doc, log);
  const transcript = pickTranscript(doc, req.transcript_id);
  const geneLen = loc.end - loc.start + 1;
  let exons = [];
  if (transcript) {
    try {
      exons = transcriptExons(doc, transcript, geneLen, log);
    } catch (err) {
      // The overlay is decoration in gene mode: an inconsistent structure drops it instead of failing.
      if (req.transcript_id !== undefined && req.transcript_id !== null) throw err;
      exons = [];
    }
  }
  const resolved = await resolveAssembly(doc.system_name, deps);
  const fasta = dnaPath(resolved);
  const regionLen = await requireRegionLength(sequence, fasta, loc.region);
  if (loc.end > regionLen) throw structureError(doc._id, 'the gene ends beyond region ' + loc.region, log);
  const ext = flankExtent(loc, regionLen, req.flank_up || 0, req.flank_down || 0);
  const length = ext.end - ext.start + 1;
  if (length > cfg.design.max_template_length) {
    throw tooLong(length, cfg.design.max_template_length, 'use mode transcript, or mode region for part of the gene');
  }
  const seq = await sequence.fetch(fasta, loc.region, ext.start, ext.end, loc.strand, { maxLength: cfg.design.max_fetch_length });
  const overlay = exons.length > 0 ? transcript : null;
  let cds = null;
  if (overlay && validCds(overlay.cds, exonsLength(exons))) {
    const s = cdnaToGeneRelative(exons, overlay.cds.start);
    const e = cdnaToGeneRelative(exons, overlay.cds.end);
    cds = { start: s + ext.effUp, end: e + ext.effUp };
  }
  return finishTemplate({
    mode: 'gene',
    id: doc._id,
    system_name: doc.system_name,
    gene_id: doc._id,
    transcript_id: overlay ? overlay.id : null,
    region: loc.region,
    start: ext.start,
    end: ext.end,
    strand: loc.strand,
    seq: seq,
    features: {
      gene: { start: ext.effUp + 1, end: ext.effUp + geneLen },
      exons: exons.map(function (e) {
        return { id: e.id, start: e.start + ext.effUp, end: e.end + ext.effUp, genomic: coords.exonGenomicRange(loc, e) };
      }),
      cds: cds,
      junctions: []
    },
    exon_count: overlay ? exons.length : null,
    mapper: coords.genomicMapper({ region: loc.region, start: ext.start, end: ext.end, strand: loc.strand }),
    resolved: resolved,
    location: loc,
    warnings: assemblyWarnings(resolved)
  });
}

async function buildTranscriptTemplate(req, deps) {
  const cfg = configOf(deps);
  const log = deps.log || console;
  const sequence = sequenceOf(deps);
  const doc = await loadGene(req, deps);
  const loc = geneLocation(doc, log);
  const transcript = pickTranscript(doc, req.transcript_id);
  if (!transcript) {
    throw new PrimerHttpError(404, 'UNKNOWN_TRANSCRIPT', 'gene ' + doc._id + ' has no transcripts', { gene_id: String(doc._id) });
  }
  const geneLen = loc.end - loc.start + 1;
  const exons = transcriptExons(doc, transcript, geneLen, log);
  const segments = coords.buildSegments(loc, exons);
  const cdnaLen = segments[segments.length - 1].t_end;
  if (transcript.length !== undefined && transcript.length !== null && transcript.length !== cdnaLen) {
    throw structureError(doc._id, 'the exons of ' + transcript.id + ' sum to ' + cdnaLen + ' nt but transcript.length is ' +
      transcript.length, log);
  }
  if (cdnaLen > cfg.design.max_template_length) throw tooLong(cdnaLen, cfg.design.max_template_length);
  const resolved = await resolveAssembly(doc.system_name, deps);
  const fasta = dnaPath(resolved);
  const regionLen = await requireRegionLength(sequence, fasta, loc.region);
  if (loc.end > regionLen) throw structureError(doc._id, 'the gene ends beyond region ' + loc.region, log);
  const parts = await fetchWindows(sequence, fasta, loc.region, loc.strand,
    segments.map(function (s) { return { start: s.g_start, end: s.g_end }; }),
    { maxLength: cfg.design.max_fetch_length, spanFetchLimit: deps.span_fetch_limit });
  const seq = parts.join('');
  if (seq.length !== cdnaLen) {
    throw structureError(doc._id, 'the spliced cDNA of ' + transcript.id + ' is ' + seq.length + ' nt, expected ' + cdnaLen, log);
  }
  const junctions = coords.junctionsFromSegments(segments);
  if (Array.isArray(transcript.exon_junctions) && !sameNumbers(transcript.exon_junctions, junctions)) {
    logSafe(log, 'warn', 'primers: derived junctions of ' + transcript.id + ' differ from exon_junctions; using the derived ones');
  }
  let gStart = Infinity;
  let gEnd = -Infinity;
  segments.forEach(function (s) {
    gStart = Math.min(gStart, s.g_start);
    gEnd = Math.max(gEnd, s.g_end);
  });
  return finishTemplate({
    mode: 'transcript',
    id: transcript.id,
    system_name: doc.system_name,
    gene_id: doc._id,
    transcript_id: transcript.id,
    region: loc.region,
    start: gStart,
    end: gEnd,
    strand: loc.strand,
    seq: seq,
    features: {
      gene: null,
      exons: segments.map(function (s) {
        return { id: s.id, start: s.t_start, end: s.t_end, genomic: { start: s.g_start, end: s.g_end } };
      }),
      cds: validCds(transcript.cds, cdnaLen) ? { start: transcript.cds.start, end: transcript.cds.end } : null,
      junctions: junctions.slice()
    },
    junctions: junctions,
    exon_count: exons.length,
    mapper: coords.splicedMapper({ region: loc.region, strand: loc.strand, segments: segments }),
    resolved: resolved,
    location: loc,
    segments: segments,
    warnings: assemblyWarnings(resolved)
  });
}

async function buildRegionTemplate(req, deps) {
  const cfg = configOf(deps);
  const sequence = sequenceOf(deps);
  const r = req.region;
  if (!r || typeof r !== 'object' || typeof r.region !== 'string' || r.region === '' || r.region.length > MAX_ID_LENGTH) {
    throw invalidRequest('region must be {region, start, end, strand} with a non-empty region name', { field: 'region' });
  }
  const strand = r.strand === undefined || r.strand === null ? 1 : r.strand;
  if (strand !== 1 && strand !== -1) throw invalidRequest('region.strand must be 1 or -1', { field: 'region.strand' });
  if (!isPosInt(r.start) || !isPosInt(r.end) || r.end < r.start) {
    throw new PrimerHttpError(400, 'REGION_OUT_OF_BOUNDS', 'region start and end must be integers with 1 <= start <= end',
      { region: r.region, start: r.start, end: r.end });
  }
  const length = r.end - r.start + 1;
  if (length > cfg.design.max_template_length) throw tooLong(length, cfg.design.max_template_length);
  if (typeof req.system_name !== 'string') throw invalidRequest('system_name is required in region mode', { field: 'system_name' });
  const resolved = await resolveAssembly(req.system_name, deps);
  const fasta = dnaPath(resolved);
  const regionLen = await requireRegionLength(sequence, fasta, r.region);
  if (r.end > regionLen) {
    throw new PrimerHttpError(400, 'REGION_OUT_OF_BOUNDS',
      'end ' + r.end + ' is beyond the end of region ' + r.region + ' (' + regionLen + ' bp)',
      { region: r.region, start: r.start, end: r.end, length: regionLen });
  }
  const seq = await sequence.fetch(fasta, r.region, r.start, r.end, strand, { maxLength: cfg.design.max_fetch_length });
  return finishTemplate({
    mode: 'region',
    id: resolved.system_name + '_' + r.region + '_' + r.start + '-' + r.end,
    system_name: resolved.system_name,
    region: r.region,
    start: r.start,
    end: r.end,
    strand: strand,
    seq: seq,
    features: {},
    mapper: coords.genomicMapper({ region: r.region, start: r.start, end: r.end, strand: strand }),
    resolved: resolved,
    warnings: assemblyWarnings(resolved)
  });
}

async function buildSequenceTemplate(req, deps) {
  const cfg = configOf(deps);
  const cleaned = cleanSequence(req.sequence, { keepCase: req.avoid_repeats === true, maxLength: cfg.design.max_template_length });
  let resolved = null;
  if (req.system_name !== undefined && req.system_name !== null) resolved = await resolveAssembly(req.system_name, deps);
  const warnings = assemblyWarnings(resolved);
  if (cleaned.records > 1) {
    // Joined, not rejected: pasting exons one after another to build a cDNA is a valid use.
    warnings.push(warning('MULTIPLE_RECORDS', cleaned.records + ' FASTA records were joined into one template; primers may span the joins'));
  }
  if (cleaned.iupac > 0) {
    warnings.push(warning('IUPAC_CONVERTED', cleaned.iupac + ' IUPAC ambiguity code' + (cleaned.iupac === 1 ? ' was' : 's were') +
      ' converted to N'));
  }
  return finishTemplate({
    mode: 'sequence',
    id: cleaned.name || 'sequence',
    system_name: resolved ? resolved.system_name : null,
    seq: cleaned.seq,
    features: {},
    resolved: resolved,
    warnings: warnings
  });
}

// buildTemplate(req, deps) -> template (internal; design.js builds the public view)
//   {mode, id, system_name, gene_id, transcript_id, region, start, end, strand, length, seq (case as read),
//    features, junctions[], exon_count, mapper (coords mapper | null), resolved (assembly | null),
//    location (gene location | null), segments (cDNA segments | null), warnings[{code, message}]}
// deps: {cfg, log, findGene(id) | mongo, resolve(system_name) | (cfg, mongo, catalog), sequence {fetch, regionLength},
//        span_fetch_limit}
async function buildTemplate(req, deps) {
  deps = deps || {};
  switch (req && req.mode) {
    case 'gene': return buildGeneTemplate(req, deps);
    case 'transcript': return buildTranscriptTemplate(req, deps);
    case 'region': return buildRegionTemplate(req, deps);
    case 'sequence': return buildSequenceTemplate(req, deps);
    default: throw invalidRequest('mode must be one of gene, transcript, region, sequence', { field: 'mode' });
  }
}

module.exports = {
  buildTemplate,
  cleanSequence,
  fetchWindows,
  flankExtent,
  findGene,
  pickTranscript,
  transcriptExons,
  cdnaToGeneRelative,
  GENE_FIELDS,
  MIN_SEQUENCE_LENGTH,
  SPAN_FETCH_LIMIT
};
