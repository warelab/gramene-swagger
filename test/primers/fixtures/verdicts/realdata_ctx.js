'use strict';

// Real-data helpers for run_realdata.test.js (read-only: mongo, /scratch/olson/fasta, blastn).
// Builds job.resolved the way check/index.js normalize is specified to (§A.2.3, §B.1) and a worker-like
// ctx (§A.8.5) whose spawnLines is blast.spawnLinesLocal at nice 10.

const fs = require('fs');
const os = require('os');
const path = require('path');

const ROOT = path.resolve(__dirname, '..', '..', '..', '..');
require(path.join(ROOT, 'api/helpers/primers/node_compat'));

const blast = require(path.join(ROOT, 'api/helpers/primers/check/blast'));
const classify = require(path.join(ROOT, 'api/helpers/primers/check/classify'));

// §10.4 primers (sorghum_bicolor reference).
const P = Object.freeze({
  P1_L: 'GATCGACAATCCGACGATAGAAG',
  P1_R: 'GTGAACATCATGCTGCCCGATG',
  P2_L: 'GGACAGCTCCACAACATATCAG',
  P2_R: 'GGACATTTGAAGCCCATGGCC',
  P3_L: 'GATATCAGTGGAATCATAAGACCG',
  P3_R: 'CATCGATATCAGGATCTGGCTT',
  P5_L: 'CGACGAGGAGGCGCTGCGATG',
  m20_7_14: 'GGACAGATCCACATCATATC',
  del18: 'GGACAGCTCCACAACATTCAG',
  J_L: 'CCAACAAAGTCATGGATGCACT',
  PA_L: 'TTACCTCTCCTCCTTCCCGATC',
  PA_R: 'GCATAGTCTCATTAGCTAGCAC'
});

function mongoConfig() {
  return require(path.join(ROOT, 'node_modules/gramene-mongodb-config'));
}

function primersConfig(overrides) {
  const config = require(path.join(ROOT, 'api/helpers/primers/config'));
  const env = {};
  if (process.env.PRIMERS_FASTA_ROOT) env.PRIMERS_FASTA_ROOT = process.env.PRIMERS_FASTA_ROOT;
  return config._build({ env, fileConfig: {}, overrides: overrides || undefined }).config;
}

// resolved.gene as §A.2.3 describes: location, transcript ids, ortholog ids grouped by system_name.
async function geneContext(mongo, geneId) {
  const genes = await mongo.genes.mongoCollection();
  if (!genes) throw new Error('mongo genes collection unavailable');
  const g = await genes.findOne({ _id: geneId }, {
    fields: { location: 1, system_name: 1, 'gene_structure.transcripts.id': 1, 'gene_structure.canonical_transcript': 1, 'homology.homologous_genes': 1 }
  });
  if (!g) throw new Error('unknown gene ' + geneId);
  const hg = (g.homology && g.homology.homologous_genes) || {};
  const ids = [].concat(hg.ortholog_one2one || [], hg.ortholog_one2many || [], hg.ortholog_many2many || []);
  const docs = ids.length ? await genes.find({ _id: { $in: ids } }, { fields: { system_name: 1 } }).toArray() : [];
  const orthologs = {};
  for (const d of docs) (orthologs[d.system_name] = orthologs[d.system_name] || []).push(String(d._id));
  for (const k of Object.keys(orthologs)) orthologs[k].sort();
  return {
    id: String(g._id),
    system_name: g.system_name,
    location: g.location,
    transcripts: g.gene_structure.transcripts.map((t) => t.id),
    canonical_transcript: g.gene_structure.canonical_transcript,
    orthologs
  };
}

async function buildResolved(opts) {
  const assemblies = require(path.join(ROOT, 'api/helpers/primers/assemblies'));
  const names = [opts.systemName].concat(opts.genomes || []);
  const out = await assemblies.resolveMany(names, { cfg: opts.cfg, mongo: opts.mongo });
  const resolved = { assemblies: {}, errors: {}, gene: null };
  for (const o of out) {
    if (o.error) resolved.errors[o.system_name] = { code: o.error.code || 'INTERNAL', message: o.error.message };
    else resolved.assemblies[o.system_name] = JSON.parse(JSON.stringify(o.resolved));
  }
  if (opts.geneId) resolved.gene = await geneContext(opts.mongo, opts.geneId);
  return resolved;
}

// Normalized request (§A.2.3): defaults filled, primers uppercased, checks sorted with specificity.
function normalizedRequest(body) {
  const checks = Array.from(new Set(['specificity'].concat(body.checks || []))).sort();
  return {
    system_name: body.system_name,
    mode: body.mode || 'region',
    gene_id: body.gene_id,
    transcript_id: body.transcript_id,
    checks,
    genomes: checks.indexOf('pangenome') >= 0 ? (body.genomes || []).slice().sort() : [],
    params: Object.assign({}, classify.DEFAULT_PARAMS, body.params || {}),
    pairs: body.pairs.map((p) => Object.assign({}, p, { left: p.left.toUpperCase(), right: p.right.toUpperCase() }))
  };
}

// ctx.tmpdir() is a private mkdtemp directory (like the worker's per-job directory); call cleanup() after the run.
function makeCtx(opts) {
  const events = { progress: [], partials: [] };
  const controller = new AbortController();
  let tmp = null;
  const ctx = {
    jobId: opts.jobId || 'realdata',
    siteKey: 'primers_test',
    signal: controller.signal,
    config: opts.cfg,
    resolved: opts.resolved,
    procs: opts.procs,
    progress: (p) => { events.progress.push(Object.assign({}, p, { running: (p.running || []).slice() })); },
    partial: (r) => { events.partials.push({ specificity: !!r.specificity, pangenome_genomes: r.pangenome ? r.pangenome.pairs.map((x) => x.genomes.length) : null }); },
    tmpdir: () => {
      if (!tmp) tmp = fs.mkdtempSync(path.join(os.tmpdir(), 'primers-realdata-'));
      return tmp;
    },
    log: opts.log || { info() {}, warn() {}, error() {} },
    spawnLines: (cmd, args, o) => blast.spawnLinesLocal(cmd, args, Object.assign({}, o, { nice: 10, signal: (o && o.signal) || controller.signal })),
    mongo: opts.mongo
  };
  const cleanup = () => {
    if (tmp) fs.rmSync(tmp, { recursive: true, force: true });
    tmp = null;
  };
  return { ctx, events, controller, cleanup };
}

module.exports = { ROOT, P, mongoConfig, primersConfig, geneContext, buildResolved, normalizedRequest, makeCtx };
