'use strict';

// Hermetic "species" for run() tests: synthetic chromosomes with planted primer sites written as real
// bgzip FASTA (catalog fsfixture.bgzipFasta, read back through sequence.js), a fake blastn that reports
// every full-length ungapped hit with <= 5 mismatches (outfmt 6 genome or cDNA columns, per query,
// subject-contiguous), a fake blastdbcmd -info, resolve-shaped assemblies and fake mongo gene docs.
//
// standardWorld() builds the world used by verdicts.test.js and pangenome.test.js; call cleanup().

const fs = require('fs');
const os = require('os');
const path = require('path');

const { bgzipFasta } = require('../catalog/fsfixture');
const { fakeMongo } = require('./fake_mongo');

const MAX_FAKE_MM = 5;
const COMP = { A: 'T', C: 'G', G: 'C', T: 'A', N: 'N' };
const SUB = { A: 'C', C: 'G', G: 'T', T: 'A' };

function revcomp(s) {
  let o = '';
  for (let i = s.length - 1; i >= 0; i--) o += COMP[s[i]] || 'N';
  return o;
}

function rng(seed) {
  let a = seed >>> 0;
  return () => {
    a = (a + 0x6D2B79F5) | 0;
    let t = Math.imul(a ^ (a >>> 15), 1 | a);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

function randomSeq(len, seed) {
  const r = rng(seed);
  let s = '';
  for (let i = 0; i < len; i++) s += 'ACGT'[Math.floor(r() * 4)];
  return s;
}

// Substitutions at distances from the primer 3' end (1 = terminal base).
function mutate(primer, distances) {
  const a = primer.split('');
  for (const d of distances || []) {
    const i = a.length - d;
    a[i] = SUB[a[i]];
  }
  return a.join('');
}

let seedCounter = 1;

class Seq {
  constructor(name, len, seed) {
    this.name = name;
    this.arr = randomSeq(len, seed != null ? seed : 7919 * seedCounter++).split('');
  }

  put(pos1, s) {
    if (pos1 < 1 || pos1 - 1 + s.length > this.arr.length) throw new RangeError('site outside ' + this.name);
    for (let i = 0; i < s.length; i++) this.arr[pos1 - 1 + i] = s[i];
    return this;
  }

  // forward-facing site: primer 5' end at p5, primer reads the plus strand
  forward(p5, primer, dist) {
    return this.put(p5, mutate(primer, dist));
  }

  // reverse-facing site: primer 5' end at p5, primer reads the minus strand
  reverse(p5, primer, dist) {
    const m = mutate(primer, dist);
    return this.put(p5 - m.length + 1, revcomp(m));
  }

  lr(start, end, L, R, o) {
    return this.forward(start, L, o && o.l).reverse(end, R, o && o.r);
  }

  rl(start, end, L, R, o) {
    return this.forward(start, R, o && o.r).reverse(end, L, o && o.l);
  }

  get seq() {
    return this.arr.join('');
  }
}

function hamming(a, s, offset, max) {
  let mm = 0;
  for (let k = 0; k < a.length; k++) {
    if (a.charCodeAt(k) !== s.charCodeAt(offset + k)) {
      mm++;
      if (mm > max) return mm;
    }
  }
  return mm;
}

function parseQueryFasta(text) {
  const out = [];
  const lines = String(text || '').split('\n');
  for (let i = 0; i < lines.length; i++) {
    const m = /^>q(\d+)$/.exec(lines[i]);
    if (m) out[Number(m[1])] = lines[i + 1];
  }
  return out;
}

// genomes spec: { sys: { taxon_id, map_id, display_name, chromosomes: [Seq], cdna: [Seq] | null,
//                         fasta = true, dnaDb = true, failures = 0 (failing blastn calls per DB),
//                         timeouts = 0 (blastn calls per DB that report timedOut) } }
// world.calls records {cmd, args, cwd} for every spawn.
function buildWorld(spec) {
  const root = fs.mkdtempSync(path.join(os.tmpdir(), 'primers-verdicts-'));
  const dbs = new Map();
  const assemblies = {};
  const calls = [];
  for (const sys of Object.keys(spec.genomes)) {
    const g = spec.genomes[sys];
    const dir = path.join(root, sys);
    fs.mkdirSync(dir, { recursive: true });
    const chroms = (g.chromosomes || []).map((c) => ({ name: c.name, seq: c.seq }));
    const cdna = g.cdna ? g.cdna.map((c) => ({ name: c.name, seq: c.seq })) : null;
    let fasta = null;
    if (g.fasta !== false && chroms.length) fasta = bgzipFasta(dir, sys, 'dna', chroms);
    const dnaDb = g.dnaDb !== false && chroms.length ? path.join(dir, sys + '.dna.toplevel') : null;
    const cdnaDb = cdna ? path.join(dir, sys + '.cdna.all') : null;
    if (dnaDb) dbs.set(dnaDb, { kind: 'genome', seqs: chroms, failuresLeft: g.failures || 0, timeoutsLeft: g.timeouts || 0 });
    if (cdnaDb) dbs.set(cdnaDb, { kind: 'cdna', seqs: cdna, failuresLeft: g.failures || 0, timeoutsLeft: g.timeouts || 0 });
    assemblies[sys] = {
      system_name: sys,
      taxon_id: g.taxon_id != null ? g.taxon_id : null,
      display_name: g.display_name || sys.toUpperCase(),
      map_id: g.map_id || 'MAP_' + sys,
      prefix: sys,
      dir,
      fasta: { dna: fasta, dna_sm: null },
      blastdb: { dna: dnaDb, cdna: cdnaDb },
      repeat_masking: 'absent',
      total_bases: chroms.reduce((a, c) => a + c.seq.length, 0) || null,
      num_sequences: chroms.length || null,
      fingerprint: 'f'.repeat(40),
      warnings: g.warnings || []
    };
  }

  async function spawnLines(cmd, args, opts) {
    const o = opts || {};
    calls.push({ cmd, args: args.slice(), cwd: o.cwd });
    await new Promise((resolve) => setImmediate(resolve));
    if (o.signal && o.signal.aborted) return { code: null, signal: 'SIGTERM', stderrTail: '', aborted: true };
    const emit = (line) => { if (o.onLine) o.onLine(line); };
    if (args[0] === '-version') {
      emit('blastn: 2.13.0+');
      emit(' Package: blast 2.13.0, build fake');
      return { code: 0, signal: null, stderrTail: '' };
    }
    const db = args[args.indexOf('-db') + 1];
    const entry = dbs.get(db);
    if (!entry) return { code: 2, signal: null, stderrTail: 'BLAST Database error: No alias or index file found for nucleotide database [' + db + ']\n' };
    if (args[0] === '-info') {
      const total = entry.seqs.reduce((a, s) => a + s.seq.length, 0);
      emit('Database: ' + path.basename(db));
      emit('\t' + entry.seqs.length.toLocaleString('en-US') + ' sequences; ' + total.toLocaleString('en-US') + ' total bases');
      return { code: 0, signal: null, stderrTail: '' };
    }
    if (entry.timeoutsLeft > 0) {
      entry.timeoutsLeft--;
      return { code: null, signal: 'SIGTERM', stderrTail: '', timedOut: true };
    }
    if (entry.failuresLeft > 0) {
      entry.failuresLeft--;
      return { code: 2, signal: null, stderrTail: 'BLAST Database error: Could not find volume or alias file (' + db + '.00) referenced in alias file (' + db + '.nal).\n' };
    }
    const fmt = args[args.indexOf('-outfmt') + 1] || '';
    const cdna = / qseq( |$)/.test(fmt);
    const primers = parseQueryFasta(o.stdin);
    primers.forEach((p, q) => {
      const rc = revcomp(p);
      const n = p.length;
      for (const subj of entry.seqs) {
        const s = subj.seq;
        for (let i = 0; i + n <= s.length; i++) {
          const mp = hamming(p, s, i, MAX_FAKE_MM);
          if (mp <= MAX_FAKE_MM) {
            const cols = ['q' + q, subj.name, n, 1, n, i + 1, i + n, 'plus', mp];
            if (cdna) cols.push(p, s.substr(i, n));
            emit(cols.join('\t'));
          }
          const mm = hamming(rc, s, i, MAX_FAKE_MM);
          if (mm <= MAX_FAKE_MM) {
            const cols = ['q' + q, subj.name, n, 1, n, i + n, i + 1, 'minus', mm];
            if (cdna) cols.push(p, revcomp(s.substr(i, n)));
            emit(cols.join('\t'));
          }
        }
      }
    });
    return { code: 0, signal: null, stderrTail: '' };
  }

  return {
    root,
    dbs,
    assemblies,
    calls,
    spawnLines,
    blastCalls: (sys) => calls.filter((c) => c.args[0] !== '-version' && c.args[0] !== '-info' && (!sys || String(c.args[c.args.indexOf('-db') + 1]).indexOf(path.join(root, sys) + path.sep) === 0)),
    cleanup: () => fs.rmSync(root, { recursive: true, force: true })
  };
}

const PRIMERS = Object.freeze({
  L: 'ACGTTGCATGCCAGTACGATCG',
  R: 'TTGACCGGATCACTGGAAGTC',
  Q_L: 'GGCATTCAGCTAGGCTTACCAG',
  Q_R: 'CAGTTGGACGCATTCGAGATC',
  S_L: 'TCGAAGGTCTCACAGGTTGCA',
  S_R: 'AGCCTTGATCCGTAAGCTGCT',
  T_L: 'GTTCACGGATAGCATCCTGAGC',
  T_R: 'ACTCGGTAAGTCCATGCGTTA'
});

function geneDoc(id, taxonId, map, region, start, end, strand, transcripts) {
  return {
    _id: id,
    name: id,
    biotype: 'protein_coding',
    taxon_id: taxonId,
    location: { map, region, start, end, strand },
    gene_structure: { transcripts: (transcripts || []).map((t) => ({ id: t })) }
  };
}

function standardWorld() {
  const { L, R, Q_L, Q_R, S_L, S_R, T_L, T_R } = PRIMERS;
  const chr = (name, len, seed) => new Seq(name, len, seed);
  const world = buildWorld({
    genomes: {
      ref: {
        taxon_id: 1001,
        map_id: 'MAP_REF',
        display_name: 'Reference',
        chromosomes: [
          chr('1', 3000, 11).lr(1001, 1300, L, R),
          chr('2', 3000, 12).lr(1001, 1350, L, R, { l: [12, 16] }),
          chr('3', 3000, 13).rl(1001, 1300, L, R, { l: [1] }),
          chr('4', 3000, 14).lr(1001, 1300, L, R, { r: [1, 2] }),
          chr('5', 3000, 15).lr(1001, 1300, L, R, { l: [3, 6, 9, 12, 15, 18] }),
          chr('6', 3000, 16).lr(101, 400, Q_L, Q_R),
          chr('7', 3000, 17).lr(101, 350, Q_L, Q_R),
          chr('8', 3000, 18).lr(501, 700, S_L, S_R),
          chr('9', 3600, 19).lr(101, 3101, T_L, T_R)
        ],
        cdna: [
          chr('REFG1.1', 800, 21).lr(101, 400, L, R),
          chr('REFG1.2', 700, 22).lr(51, 350, L, R),
          chr('OTHER_G.1', 600, 23).lr(101, 380, L, R, { l: [10] }),
          chr('odd_tx', 500, 24).lr(51, 330, L, R),
          chr('UNMAPPED_X.3', 500, 25).lr(21, 300, L, R, { r: [9] }),
          chr('REFG8.1', 400, 26).lr(51, 250, S_L, S_R)
        ]
      },
      pan_a: {
        taxon_id: 2001,
        map_id: 'MAP_A',
        chromosomes: [chr('1', 3000, 31).lr(1101, 1400, L, R), chr('3', 5000, 32).lr(101, 4600, T_L, T_R)],
        cdna: [chr('PANA_G1.1', 500, 33).lr(51, 350, L, R)]
      },
      pan_b: {
        taxon_id: 2002,
        map_id: 'MAP_B',
        chromosomes: [chr('1', 3000, 41).lr(501, 810, L, R, { r: [3] })],
        cdna: [chr('PANB_G7.1', 500, 42).lr(51, 350, L, R, { l: [2] })]
      },
      pan_c: {
        taxon_id: 2003,
        map_id: 'MAP_C',
        chromosomes: [chr('1', 3000, 51).lr(201, 500, L, R), chr('2', 3000, 52).lr(1001, 1280, L, R)]
      },
      pan_d: { taxon_id: 2004, map_id: 'MAP_D', chromosomes: [chr('1', 3000, 61).lr(1001, 1300, L, R, { r: [1, 2] })] },
      pan_e: { taxon_id: 2005, map_id: 'MAP_E', chromosomes: [chr('1', 3000, 71).lr(1001, 1300, L, R)], dnaDb: false },
      pan_f: { taxon_id: 2006, map_id: 'MAP_F', chromosomes: [chr('1', 3000, 81).lr(1001, 1300, L, R)], failures: Infinity },
      pan_g: { taxon_id: 2007, map_id: 'MAP_G', chromosomes: [chr('1', 3000, 91).lr(1001, 1300, L, R)], failures: 1 },
      pan_h: { taxon_id: 2008, map_id: 'MAP_H', chromosomes: [chr('1', 3000, 101).lr(1001, 1300, L, R)], fasta: false },
      pan_t: { taxon_id: 2009, map_id: 'MAP_T', chromosomes: [chr('1', 3000, 111).lr(1001, 1300, L, R)], timeouts: 1 }
    }
  });
  world.primers = PRIMERS;
  world.genes = [
    geneDoc('REFG1', 1001, 'MAP_REF', '1', 500, 2000, 1, ['REFG1.1', 'REFG1.2']),
    geneDoc('REFG8', 1001, 'MAP_REF', '8', 400, 900, -1, ['REFG8.1']),
    geneDoc('OTHER_G', 1001, 'MAP_REF', '2', 900, 1500, 1, ['OTHER_G.1']),
    geneDoc('ODDG', 1001, 'MAP_REF', '3', 2000, 2500, 1, ['odd_tx']),
    geneDoc('PANA_G1', 2001, 'MAP_A', '1', 1000, 1500, 1, ['PANA_G1.1']),
    geneDoc('PANB_G7', 2002, 'MAP_B', '1', 400, 900, 1, ['PANB_G7.1']),
    geneDoc('PANC_G2', 2003, 'MAP_C', '2', 900, 1400, 1, []),
    geneDoc('PANG_G1', 2007, 'MAP_G', '1', 900, 1400, 1, [])
  ];
  world.gene = {
    id: 'REFG1',
    system_name: 'ref',
    location: { map: 'MAP_REF', region: '1', start: 500, end: 2000, strand: 1 },
    transcripts: ['REFG1.1', 'REFG1.2'],
    orthologs: { pan_a: ['PANA_G1'], pan_c: ['PANC_G2'], pan_g: ['PANG_G1'], pan_h: ['PANH_G1'] }
  };
  world.mongo = (opts) => fakeMongo(world.genes, opts);
  return world;
}

// Worker-like ctx for the fake world. tmpdir() is a private mkdtemp directory inside world.root (like the
// worker's per-job directory); opts.tmpdir replaces it.
function worldCtx(world, opts) {
  const o = opts || {};
  const events = { progress: [], partials: [] };
  const controller = new AbortController();
  let jobTmp = null;
  const privateTmp = () => {
    if (!jobTmp) jobTmp = fs.mkdtempSync(path.join(world.root, 'job-'));
    return jobTmp;
  };
  const DEFAULTS = require('../../../../api/helpers/primers/config').DEFAULTS;
  const ctx = {
    jobId: 'fakejob',
    siteKey: 'primers_test',
    signal: controller.signal,
    config: { blastn: '/fake/bin/blastn', blastdbcmd: '/fake/bin/blastdbcmd', check: Object.assign({}, DEFAULTS.check, o.check || {}) },
    resolved: { assemblies: world.assemblies, errors: o.errors || {}, gene: o.gene === undefined ? world.gene : o.gene },
    procs: o.procs || 4,
    progress: (p) => {
      events.progress.push(Object.assign({}, p, { running: (p.running || []).slice() }));
      if (o.onProgress) o.onProgress(p, controller);
    },
    partial: (r) => {
      events.partials.push(JSON.parse(JSON.stringify(r)));
    },
    tmpdir: o.tmpdir || privateTmp,
    log: { info() {}, warn() {}, error() {} },
    spawnLines: world.spawnLines,
    mongo: o.mongo !== undefined ? o.mongo : world.mongo().mongo
  };
  return { ctx, events, controller };
}

function request(body) {
  const classify = require('../../../../api/helpers/primers/check/classify');
  const checks = Array.from(new Set(['specificity'].concat(body.checks || []))).sort();
  return {
    system_name: body.system_name || 'ref',
    mode: body.mode || 'region',
    gene_id: body.gene_id,
    transcript_id: body.transcript_id,
    checks,
    genomes: checks.indexOf('pangenome') >= 0 ? (body.genomes || []).slice().sort() : [],
    params: Object.assign({}, classify.DEFAULT_PARAMS, body.params || {}),
    pairs: body.pairs
  };
}

module.exports = { MAX_FAKE_MM, revcomp, randomSeq, mutate, Seq, buildWorld, standardWorld, worldCtx, request, PRIMERS };
