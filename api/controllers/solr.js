'use strict';

var _ = require('lodash');
var solrHelper = require('../helpers/solr');
const axios = require('axios');

// /gene_lists/validate guards. The effective ceiling used to be body-parser's implicit 100kb
// default, which surfaced as a bare 413; MAX_IDS makes the documented swagger maxItems the real
// limit and returns a 400 that says so.
const MAX_IDS = 6000;
const MAX_ID_LENGTH = 255;
const CHUNK_SIZE = 200;
const MAX_ROWS_PER_CHUNK = 5000;   // one input can legitimately match many genes (up to ~15 seen)
const SEP = '|';                   // {!terms} separator; gene ids never contain '|'

module.exports = {
  genes: genes,
  suggestions: suggestions,
  validate: validate
};

function genes(req, res) {
  solrRequest(req, res, solrHelper.streamGenes, {facet: true, stats: true});
}

function suggestions(req, res) {
  solrRequest(req, res, solrHelper.streamSuggestions);
}

function solrRequest(req, res, streamMethod, additionalParams) {
  var params = _.mapValues(req.swagger.params, 'value');
  if(additionalParams) {
    _.assign(params, additionalParams);
  }
  var solrStream = streamMethod(params);

  solrStream.pipe(res);
}

async function validate(req, res) {
  // Lookup only. Resolves each submitted identifier against BOTH the stable `id` and the
  // alternate-id field, and reports what it found. It deliberately does NOT compute a hash and
  // does NOT touch the genes core — merely validating a list the user may never save used to
  // write saved_search tags across a 5.4M-doc core, fire-and-forget, after the response had
  // already been sent. Saving (hash + core tagging) now lives in POST /gene_lists.
  const solrURL = solrHelper.genesURL + '/select';

  const raw = req.body;
  if (!Array.isArray(raw)) {
    return res.status(400).send('body must be a JSON array of identifier strings');
  }
  if (raw.length > MAX_IDS) {
    return res.status(400).send('too many identifiers: ' + raw.length + ' (limit ' + MAX_IDS + ')');
  }
  for (const v of raw) {
    if (typeof v !== 'string') {
      return res.status(400).send('every identifier must be a string');
    }
    if (v.length > MAX_ID_LENGTH) {
      return res.status(400).send('identifier longer than ' + MAX_ID_LENGTH + ' characters');
    }
  }

  // The client is responsible for trimming/splitting; we only drop empties and duplicates.
  const submitted = [...new Set(raw.filter(s => s.length))];
  if (!submitted.length) {
    return res.json({resolved: [], ambiguous: [], unknown: []});
  }

  // submitted (as typed) -> Set of matching gene ids
  const matches = new Map(submitted.map(s => [s, new Set()]));
  // lowercased submitted -> the forms the user typed (alt_id is indexed lowercased, so a hit
  // has to be attributed back to every input that lowercases to it)
  const byLower = new Map();
  submitted.forEach(s => {
    const lc = s.toLowerCase();
    if (!byLower.has(lc)) byLower.set(lc, []);
    byLower.get(lc).push(s);
  });

  async function lookup(chunk) {
    // {!terms} takes literal values — no Lucene escaping needed and no way to inject query
    // syntax, unlike the old `id:(a OR b OR ...)` string interpolation. It also bypasses
    // analysis, which is why alt_id (a `lowercase` field) must be lowercased by hand here.
    const idTerms  = chunk.join(SEP);
    const altTerms = chunk.map(s => s.toLowerCase()).join(SEP);
    const q = `{!terms f=id separator=${SEP}}${idTerms}`;
    const fq = `{!terms f=alt_id separator=${SEP}}${altTerms}`;
    // one gene can match several submitted ids, and one submitted id several genes, so ask for
    // more rows than inputs rather than capping at chunk.length (which silently truncated).
    const common = {fl: 'id,alt_id', rows: MAX_ROWS_PER_CHUNK, wt: 'json'};
    const [byId, byAlt] = await Promise.all([
      axios.post(solrURL, new URLSearchParams({...common, q}).toString()),
      axios.post(solrURL, new URLSearchParams({...common, q: fq}).toString())
    ]);
    return [...byId.data.response.docs, ...byAlt.data.response.docs];
  }

  try {
    for (let i = 0; i < submitted.length; i += CHUNK_SIZE) {
      const chunk = submitted.slice(i, i + CHUNK_SIZE);
      const inChunk = new Set(chunk);
      const docs = await lookup(chunk);
      for (const doc of docs) {
        // exact stable-id hit
        if (inChunk.has(doc.id)) matches.get(doc.id).add(doc.id);
        // alternate-id hit: attribute to whichever submitted forms lowercase to it
        for (const a of (doc.alt_id || [])) {
          for (const typed of (byLower.get(String(a).toLowerCase()) || [])) {
            if (inChunk.has(typed)) matches.get(typed).add(doc.id);
          }
        }
      }
    }
  } catch (error) {
    console.error('validate: solr lookup failed:', error && error.message || error);
    return res.status(502).send('gene lookup failed');
  }

  // An input matching more than one gene is reported as ambiguous for the client to resolve —
  // including when one of those matches is an exact stable id.
  const resolved = [], ambiguous = [], unknown = [];
  for (const input of submitted) {
    const hits = [...matches.get(input)].sort();
    if (hits.length === 0)      unknown.push(input);
    else if (hits.length === 1) resolved.push({input, id: hits[0]});
    else                        ambiguous.push({input, matches: hits});
  }
  res.json({resolved, ambiguous, unknown});
}
