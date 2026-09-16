'use strict';

require('../../../api/helpers/primers/node_compat');

const test = require('node:test');
const should = require('should');
const fs = require('fs');
const os = require('os');
const path = require('path');
const { execFileSync } = require('child_process');

const rm = require('../../../api/helpers/primers/repeat_mask');
const template = require('../../../api/helpers/primers/template');
const coords = require('../../../api/helpers/primers/coords');
const { revcomp } = require('../../../api/helpers/primers/sequence');
const { PrimerHttpError } = require('../../../api/helpers/primers/errors');
const stubs = require('../fixtures/design/stubs');
const fx = require('../fixtures/catalog/fsfixture');

const DNA = '/fake/syn/dna/Syn.dna.toplevel.fa.gz';
const DNA_SM = '/fake/syn/dna/Syn.dna_sm.toplevel.fa.gz';
const BLASTDB = '/fake/syn/Syn.dna.toplevel';
const C_SEQ = fx.randomSeq(1000, 21);
const BIG_SEQ = fx.randomSeq(45000, 5);

function lowerAt(seq, ranges) {
  let out = seq;
  ranges.forEach(function (r) {
    out = out.slice(0, r[0] - 1) + out.slice(r[0] - 1, r[1]).toLowerCase() + out.slice(r[1]);
  });
  return out;
}

const SM_SEQ = lowerAt(C_SEQ, [[11, 20], [290, 309], [560, 569]]);

function files(extra) {
  const f = {};
  f[DNA] = { c: { length: 1000, windows: [{ start: 1, seq: C_SEQ }] }, big: { length: 45000, windows: [{ start: 1, seq: BIG_SEQ }] } };
  f[DNA_SM] = { c: { length: 1000, windows: [{ start: 1, seq: SM_SEQ }] } };
  return Object.assign(f, extra || {});
}

function asm(overrides) {
  return stubs.resolvedStub(Object.assign({
    system_name: 'syn',
    fasta: { dna: DNA, dna_sm: null },
    blastdb: { dna: BLASTDB, cdna: null },
    repeat_masking: 'unmasked_copy'
  }, overrides || {}));
}

const NC1 = {
  _id: 'NC1', system_name: 'syn', location: { region: 'c', start: 101, end: 700, strand: -1 },
  gene_structure: {
    exons: [{ id: 'e1', start: 1, end: 150 }, { id: 'e2', start: 401, end: 600 }],
    transcripts: [{ id: 'NC1.1', length: 350, exons: ['e1', 'e2'] }]
  }
};

async function buildTemplate(req, assembly, extra) {
  const deps = Object.assign({
    cfg: stubs.cfg(),
    log: stubs.silentLog,
    findGene: async function () { return JSON.parse(JSON.stringify(NC1)); },
    resolve: async function () { return assembly; },
    sequence: stubs.sequenceStub(files())
  }, extra || {});
  return template.buildTemplate(req, deps);
}

function maskDeps(extra) {
  return Object.assign({ cfg: stubs.cfg(), log: stubs.silentLog, sequence: stubs.sequenceStub(files()), cache: rm.createLru(10) }, extra || {});
}

function hsp(q, s, e, len) {
  return 'q' + q + '\t' + s + '\t' + e + '\t' + (len === undefined ? Math.abs(e - s) + 1 : len);
}

// runBlast stub: records its call and replays `lines`.
function blastStub(lines, record) {
  return async function (opts) {
    if (record) record.push(opts);
    lines.forEach(function (l) { opts.onLine(l); });
  };
}

// ---- primitives -----------------------------------------------------------------------------------

test('applyMask: n_mask -> N, three_prime -> lowercase; unmasked bases uppercase', function () {
  rm.applyMask('acgtACGTacgt', [[2, 3], [10, 1]], 'n_mask').should.equal('ANNNACGTANGT');
  rm.applyMask('acgtACGTacgt', [[2, 3], [10, 1]], 'three_prime').should.equal('AcgtACGTAcGT');
  rm.applyMask('acgt', [], 'three_prime').should.equal('ACGT');
  rm.applyMask('ACGTACGT', [[7, 10]], 'n_mask').should.equal('ACGTACNN');
  rm.maskTags('n_mask').should.eql({ PRIMER_MAX_NS_ACCEPTED: 0 });
  rm.maskTags('three_prime').should.eql({ PRIMER_LOWERCASE_MASKING: 1 });
});

test('megablastArgs: the §B.14 megablast argument array from the config', function () {
  rm.megablastArgs('/scratch/x/Sb.dna.toplevel', stubs.cfg().repeat_mask).should.eql([
    '-task', 'megablast', '-db', '/scratch/x/Sb.dna.toplevel', '-query', '-', '-dust', 'no', '-soft_masking', 'false',
    '-evalue', '1e-10', '-perc_identity', '85', '-max_target_seqs', '500', '-max_hsps', '200', '-num_threads', '2',
    '-outfmt', '6 qseqid qstart qend length']);
  rm.megablastArgs('/x/multi.dna.toplevel.nal', stubs.cfg().repeat_mask)[3].should.equal('/x/multi.dna.toplevel');
  (function () { rm.megablastArgs('/x/a b', stubs.cfg().repeat_mask); }).should.throw(TypeError);
  (function () { rm.megablastArgs('-remote', stubs.cfg().repeat_mask); }).should.throw(TypeError);
});

test('parseHspLine accepts outfmt "6 qseqid qstart qend length" lines only', function () {
  rm.parseHspLine('q3\t10\t250\t241').should.eql({ query: 3, qstart: 10, qend: 250, length: 241 });
  should(rm.parseHspLine('Warning: something')).be.null();
  should(rm.parseHspLine('x3\t10\t250\t241')).be.null();
  should(rm.parseHspLine('q3\t10\tabc\t241')).be.null();
  should(rm.parseHspLine('')).be.null();
});

test('chunkQueries: 20 kb chunks overlapping by 1 kb cover the whole template', function () {
  const q = rm.chunkQueries('a'.repeat(45000), 20000, rm.CHUNK_OVERLAP);
  q.map(function (x) { return [x.offset, x.seq.length, x.lo, x.hi]; }).should.eql([[0, 20000, 1, 20000], [19000, 20000, 1, 20000], [38000, 7000, 1, 7000]]);
  q[0].seq.should.equal('A'.repeat(20000));
  rm.chunkQueries('ACGT'.repeat(100), 20000, 1000).length.should.equal(1);
});

test('depthRuns: >= min_depth HSPs of >= min_hsp_len (self-hit counts); short HSPs and pads ignored', function () {
  const queries = [{ seq: 'A'.repeat(300), offset: -100, lo: 101, hi: 250 }];
  const opts = { minLen: 50, minDepth: 3 };
  // depth 3 over query 1..300, projected through lo..hi -> template 1..150
  rm.depthRuns([{ query: 0, qstart: 1, qend: 300, length: 300 }, { query: 0, qstart: 1, qend: 300, length: 300 },
    { query: 0, qstart: 300, qend: 1, length: 300 }], queries, 150, opts).should.eql([[1, 150]]);
  // depth 2 only
  rm.depthRuns([{ query: 0, qstart: 1, qend: 300, length: 300 }, { query: 0, qstart: 1, qend: 300, length: 300 }], queries, 150, opts).should.eql([]);
  // the third HSP is too short
  rm.depthRuns([{ query: 0, qstart: 1, qend: 300, length: 300 }, { query: 0, qstart: 1, qend: 300, length: 300 },
    { query: 0, qstart: 120, qend: 160, length: 41 }], queries, 150, opts).should.eql([]);
  // a deep pile only inside the pad projects nothing
  rm.depthRuns([1, 2, 3].map(function () { return { query: 0, qstart: 1, qend: 100, length: 100 }; }), queries, 150, opts).should.eql([]);
});

// ---- blast depth mask -------------------------------------------------------------------------------

test('genomic template: 3 chunk queries, HSP depth projected with the chunk offset, BLAST_DEPTH_MASK', async function () {
  const assembly = asm();
  const t = await buildTemplate({ mode: 'region', system_name: 'syn', region: { region: 'big', start: 1, end: 45000, strand: 1 } }, assembly);
  const calls = [];
  const lines = [hsp(1, 100, 400), hsp(1, 100, 400), hsp(1, 400, 100),
    hsp(0, 1, 40), hsp(0, 1, 40), hsp(0, 1, 40), hsp(0, 1, 40),
    hsp(2, 10, 100), hsp(2, 10, 100), 'garbage line'];
  const out = await rm.repeatMask(t, { mode: 'n_mask' }, maskDeps({ runBlast: blastStub(lines, calls) }));
  calls.length.should.equal(1);
  calls[0].bin.should.equal(stubs.cfg().blastn);
  calls[0].args.should.eql(rm.megablastArgs(BLASTDB, stubs.cfg().repeat_mask));
  calls[0].timeoutMs.should.equal(20000);
  calls[0].tmpDir.should.equal(stubs.cfg().tmp_dir);
  calls[0].input.split('>').filter(Boolean).map(function (r) { return r.split('\n')[0]; }).should.eql(['q0', 'q1', 'q2']);
  out.mask.should.eql([[19100, 301]]);
  out.mask_source.should.equal('blast_depth');
  out.masked.should.be.true();
  out.masked_fraction.should.equal(Math.round(301 / 45000 * 10000) / 10000);
  out.warnings.map(function (w) { return w.code; }).should.eql(['BLAST_DEPTH_MASK']);
  out.seq.slice(19099, 19400).should.equal('N'.repeat(301));
  out.seq.slice(19098, 19099).should.equal(BIG_SEQ.charAt(19098));
  out.seq.length.should.equal(45000);
  out.tags.should.eql({ PRIMER_MAX_NS_ACCEPTED: 0 });
});

test('transcript template: exon +- 100 bp queries (clipped at the region start) are projected onto the cDNA', async function () {
  const assembly = asm();
  const t = await buildTemplate({ mode: 'transcript', gene_id: 'NC1' }, assembly);
  const calls = [];
  const lines = [hsp(0, 1, 350), hsp(0, 1, 350), hsp(0, 1, 350),
    hsp(1, 150, 160, 60), hsp(1, 150, 160, 60), hsp(1, 160, 150, 60),
    hsp(1, 1, 100), hsp(1, 1, 100), hsp(1, 1, 100)];
  const out = await rm.repeatMask(t, { mode: 'three_prime' }, maskDeps({ runBlast: blastStub(lines, calls) }));
  const records = calls[0].input.split('>').filter(Boolean).map(function (r) { return r.split('\n')[1]; });
  records.should.eql([revcomp(C_SEQ.slice(450, 800)), revcomp(C_SEQ.slice(0, 400))]);
  out.mask.should.eql([[1, 150], [200, 11]]);
  out.mask_source.should.equal('blast_depth');
  out.seq.should.equal(t.seq.slice(0, 150).toLowerCase() + t.seq.slice(150, 199) + t.seq.slice(199, 210).toLowerCase() + t.seq.slice(210));
  out.tags.should.eql({ PRIMER_LOWERCASE_MASKING: 1 });
});

test('megablast timeout is min(repeat_mask.timeout_ms, time left before the design deadline)', async function () {
  const t = await buildTemplate({ mode: 'region', system_name: 'syn', region: { region: 'c', start: 1, end: 500 } }, asm());
  const calls = [];
  await rm.repeatMask(t, { deadline: 1000 + 5000 }, maskDeps({ now: function () { return 1000; }, runBlast: blastStub([], calls) }));
  calls[0].timeoutMs.should.equal(5000);
});

test('megablast failure or timeout -> REPEAT_MASK_FAILED and an unmasked template', async function () {
  const t = await buildTemplate({ mode: 'region', system_name: 'syn', region: { region: 'c', start: 1, end: 500 } }, asm());
  for (const code of ['BLAST_TIMEOUT', 'BLAST_FAILED']) {
    const out = await rm.repeatMask(t, {}, maskDeps({
      runBlast: async function () { const e = new Error('boom'); e.code = code; throw e; }
    }));
    out.mask.should.eql([]);
    should(out.mask_source).be.null();
    out.masked.should.be.false();
    out.seq.should.equal(t.seq.toUpperCase());
    out.warnings.map(function (w) { return w.code; }).should.eql(['REPEAT_MASK_FAILED']);
  }
});

test('an aborted design deadline propagates its error instead of REPEAT_MASK_FAILED', async function () {
  const t = await buildTemplate({ mode: 'region', system_name: 'syn', region: { region: 'c', start: 1, end: 500 } }, asm());
  const ac = new AbortController();
  const reason = new PrimerHttpError(504, 'DEADLINE_EXCEEDED', 'the design did not finish', {});
  let err = null;
  try {
    await rm.repeatMask(t, { signal: ac.signal }, maskDeps({
      runBlast: function (opts) {
        return new Promise(function (resolve, reject) {
          opts.signal.addEventListener('abort', function () { reject(opts.signal.reason); });
          setTimeout(function () { ac.abort(reason); }, 5);
        });
      }
    }));
  } catch (e) {
    err = e;
  }
  should(err).equal(reason);
});

test('LRU cache: a repeated template does not re-run megablast; different exon structures use different keys', async function () {
  const assembly = asm();
  const t = await buildTemplate({ mode: 'transcript', gene_id: 'NC1' }, assembly);
  const calls = [];
  const deps = maskDeps({ runBlast: blastStub([hsp(0, 1, 350), hsp(0, 1, 350), hsp(0, 1, 350)], calls) });
  const a = await rm.repeatMask(t, {}, deps);
  const b = await rm.repeatMask(t, {}, deps);
  calls.length.should.equal(1);
  b.mask.should.eql(a.mask);
  b.mask_source.should.equal('blast_depth');
  b.warnings.map(function (w) { return w.code; }).should.eql(['BLAST_DEPTH_MASK']);
  const other = Object.assign({}, t, { segments: t.segments.map(function (s) { return Object.assign({}, s); }) });
  other.segments[1].g_start += 1;
  rm.blastCacheKey(other, assembly, stubs.cfg().repeat_mask).should.not.equal(rm.blastCacheKey(t, assembly, stubs.cfg().repeat_mask));
  const lru = rm.createLru(2);
  lru.set('a', 1); lru.set('b', 2); lru.get('a'); lru.set('c', 3);
  lru.has('a').should.be.true();
  lru.has('b').should.be.false();
  lru.size.should.equal(2);
});

// ---- soft mask -------------------------------------------------------------------------------------

test('soft_masked genome: region template on - strand takes lowercase runs of dna_sm in template orientation', async function () {
  const assembly = asm({ fasta: { dna: DNA, dna_sm: DNA_SM }, repeat_masking: 'soft_masked' });
  const t = await buildTemplate({ mode: 'region', system_name: 'syn', region: { region: 'c', start: 1, end: 100, strand: -1 } }, assembly);
  const calls = [];
  const out = await rm.repeatMask(t, {}, maskDeps({ runBlast: blastStub([], calls) }));
  calls.length.should.equal(0);
  out.mask_source.should.equal('softmask');
  out.mask.should.eql([[81, 10]]);
  out.warnings.should.eql([]);
  out.seq.slice(80, 90).should.equal('NNNNNNNNNN');
});

test('soft_masked genome: transcript template splices dna_sm so the mask follows the exons', async function () {
  const assembly = asm({ fasta: { dna: DNA, dna_sm: DNA_SM }, repeat_masking: 'soft_masked' });
  const t = await buildTemplate({ mode: 'transcript', gene_id: 'NC1' }, assembly);
  const out = await rm.repeatMask(t, {}, maskDeps());
  out.mask_source.should.equal('softmask');
  out.mask.should.eql([[132, 10], [151, 11]]);
});

test('soft mask that does not match the template sequence -> REPEAT_MASK_FAILED', async function () {
  const bad = files();
  bad[DNA_SM] = { c: { length: 1000, windows: [{ start: 1, seq: fx.randomSeq(1000, 99).toLowerCase() }] } };
  const assembly = asm({ fasta: { dna: DNA, dna_sm: DNA_SM }, repeat_masking: 'soft_masked' });
  const t = await buildTemplate({ mode: 'region', system_name: 'syn', region: { region: 'c', start: 1, end: 100 } }, assembly);
  const out = await rm.repeatMask(t, {}, maskDeps({ sequence: stubs.sequenceStub(bad) }));
  should(out.mask_source).be.null();
  out.warnings.map(function (w) { return w.code; }).should.eql(['REPEAT_MASK_FAILED']);
});

// ---- method selection ------------------------------------------------------------------------------

test('sequence mode: the user\'s lowercase is the mask; otherwise megablast with system_name; else NO_REPEAT_MASK', async function () {
  const lower = await buildTemplate({ mode: 'sequence', sequence: 'ACGTACGTacgtacgtACGTACGT', avoid_repeats: true }, null);
  const a = await rm.repeatMask(lower, {}, maskDeps());
  a.mask_source.should.equal('user_lowercase');
  a.mask.should.eql([[9, 8]]);
  a.seq.should.equal('ACGTACGTNNNNNNNNACGTACGT');

  const plain = await buildTemplate({ mode: 'sequence', sequence: 'ACGTACGTACGTACGTACGTACGT', avoid_repeats: true }, null);
  const b = await rm.repeatMask(plain, {}, maskDeps());
  should(b.mask_source).be.null();
  b.warnings.map(function (w) { return w.code; }).should.eql(['NO_REPEAT_MASK']);

  const withGenome = await buildTemplate({ mode: 'sequence', sequence: 'ACGTACGTACGTACGTACGTACGT', avoid_repeats: true, system_name: 'syn' }, asm());
  const calls = [];
  const c = await rm.repeatMask(withGenome, {}, maskDeps({ runBlast: blastStub([hsp(0, 1, 24, 60), hsp(0, 1, 24, 60), hsp(0, 1, 24, 60)], calls) }));
  calls.length.should.equal(1);
  c.mask_source.should.equal('blast_depth');
  c.mask.should.eql([[1, 24]]);
  c.warnings.map(function (w) { return w.code; }).should.eql(['BLAST_DEPTH_MASK', 'MOSTLY_REPEAT']);
});

test('NO_REPEAT_MASK when the genome has no BLAST DB or megablast masking is disabled', async function () {
  const noDb = asm({ blastdb: { dna: null, cdna: null } });
  const t1 = await buildTemplate({ mode: 'region', system_name: 'syn', region: { region: 'c', start: 1, end: 100 } }, noDb);
  (await rm.repeatMask(t1, {}, maskDeps())).warnings.map(function (w) { return w.code; }).should.eql(['NO_REPEAT_MASK']);
  const t2 = await buildTemplate({ mode: 'region', system_name: 'syn', region: { region: 'c', start: 1, end: 100 } }, asm());
  const calls = [];
  const out = await rm.repeatMask(t2, {}, maskDeps({ cfg: stubs.cfg({ repeat_mask: { enabled: false } }), runBlast: blastStub([], calls) }));
  calls.length.should.equal(0);
  out.warnings.map(function (w) { return w.code; }).should.eql(['NO_REPEAT_MASK']);
});

test('MOSTLY_REPEAT above 80% masked; invalid mode is 400', async function () {
  const assembly = asm({ fasta: { dna: DNA, dna_sm: DNA_SM }, repeat_masking: 'soft_masked' });
  const allLower = files();
  allLower[DNA_SM] = { c: { length: 1000, windows: [{ start: 1, seq: C_SEQ.slice(0, 81).toLowerCase() + C_SEQ.slice(81) }] } };
  const t = await buildTemplate({ mode: 'region', system_name: 'syn', region: { region: 'c', start: 1, end: 100 } }, assembly);
  const out = await rm.repeatMask(t, {}, maskDeps({ sequence: stubs.sequenceStub(allLower) }));
  out.masked_fraction.should.equal(0.81);
  out.warnings.map(function (w) { return w.code; }).should.eql(['MOSTLY_REPEAT']);
  let err = null;
  try { await rm.repeatMask(t, { mode: 'hard' }, maskDeps()); } catch (e) { err = e; }
  err.code.should.equal('INVALID_REQUEST');
});

// ---- the real spawn path, with fake blastn scripts -----------------------------------------------

function script(dir, name, body) {
  const file = path.join(dir, name);
  fs.writeFileSync(file, '#!/bin/sh\n' + body + '\n');
  fs.chmodSync(file, 0o755);
  return file;
}

function pgrepCount(pattern) {
  try {
    return execFileSync('pgrep', ['-f', pattern], { encoding: 'utf8' }).trim().split('\n').filter(Boolean).length;
  } catch (e) {
    return e.status === 1 ? 0 : -1;
  }
}

test('runMegablast: niced spawn, stdin, streamed lines (including an unterminated last line), exit codes', async function (t) {
  const dir = fs.mkdtempSync(path.join(os.tmpdir(), 'primers-rm-'));
  t.after(function () { fs.rmSync(dir, { recursive: true, force: true }); });
  const ok = script(dir, 'ok.sh', 'cat > "$(dirname "$0")/stdin.txt"\necho "$*" > "$(dirname "$0")/args.txt"\n' +
    'printf \'q0\\t1\\t60\\t60\\nq0\\t1\\t60\\t60\\n\'\nprintf \'q0\\t5\\t64\\t60\'');
  const lines = [];
  await rm.runMegablast({ bin: ok, args: ['-task', 'megablast', '-num_threads', '2'], input: '>q0\nACGT\n', timeoutMs: 10000,
    onLine: function (l) { lines.push(l); } });
  lines.should.eql(['q0\t1\t60\t60', 'q0\t1\t60\t60', 'q0\t5\t64\t60']);
  fs.readFileSync(path.join(dir, 'stdin.txt'), 'utf8').should.equal('>q0\nACGT\n');
  fs.readFileSync(path.join(dir, 'args.txt'), 'utf8').trim().should.equal('-task megablast -num_threads 2');

  const fail = script(dir, 'fail.sh', 'cat > /dev/null\necho "BLAST Database error: No alias or index file found for /scratch/x/db" >&2\nexit 2');
  let err = null;
  try { await rm.runMegablast({ bin: fail, args: [], input: '', timeoutMs: 10000, onLine: function () {} }); } catch (e) { err = e; }
  err.code.should.equal('BLAST_FAILED');
  err.message.should.match(/code 2/);
  err.message.should.not.match(/\/scratch\/x/);

  err = null;
  try { await rm.runMegablast({ bin: path.join(dir, 'missing'), args: [], input: '', timeoutMs: 10000, onLine: function () {} }); } catch (e) { err = e; }
  err.code.should.equal('BLAST_FAILED');
});

test('runMegablast: NCBI_DONT_USE_NCBIRC=1; cwd a private 0700 directory under tmpDir (never tmpdir itself), removed afterwards', async function (t) {
  const dir = fs.mkdtempSync(path.join(os.tmpdir(), 'primers-rm-'));
  t.after(function () { fs.rmSync(dir, { recursive: true, force: true }); });
  rm.SPAWN_ENV.should.eql({ PATH: '/usr/bin:/bin', NCBI_DONT_USE_NCBIRC: '1' });
  const probe = script(dir, 'probe.sh', 'cat > /dev/null\npwd -P\necho "NCBIRC=$NCBI_DONT_USE_NCBIRC HOME=${HOME:-unset}"\nstat -c %a .');
  const base = fs.mkdtempSync(path.join(dir, 'base-'));
  const lines = [];
  await rm.runMegablast({ bin: probe, args: [], input: '>q0\nACGT\n', timeoutMs: 10000, tmpDir: base, onLine: function (l) { lines.push(l); } });
  path.dirname(lines[0]).should.equal(fs.realpathSync(base));
  path.basename(lines[0]).indexOf('primers-megablast-').should.equal(0);
  lines.slice(1).should.eql(['NCBIRC=1 HOME=unset', '700']);
  fs.existsSync(lines[0]).should.be.false();
  fs.readdirSync(base).should.eql([]);

  // default base: a private directory under os.tmpdir(), never os.tmpdir() itself; removed after a failure too
  const cwdFile = path.join(dir, 'failcwd.txt');
  const failing = script(dir, 'failing.sh', 'cat > /dev/null\npwd -P > ' + cwdFile + '\nexit 2');
  let err = null;
  try { await rm.runMegablast({ bin: failing, args: [], input: '', timeoutMs: 10000, onLine: function () {} }); } catch (e) { err = e; }
  err.code.should.equal('BLAST_FAILED');
  const failCwd = fs.readFileSync(cwdFile, 'utf8').trim();
  failCwd.should.not.equal(fs.realpathSync(os.tmpdir()));
  path.dirname(failCwd).should.equal(fs.realpathSync(os.tmpdir()));
  fs.existsSync(failCwd).should.be.false();

  // a working directory that cannot be created is BLAST_FAILED (-> REPEAT_MASK_FAILED), and nothing is spawned
  err = null;
  try {
    await rm.runMegablast({ bin: failing, args: [], input: '', timeoutMs: 10000, tmpDir: path.join(dir, 'missing'), onLine: function () {} });
  } catch (e) { err = e; }
  err.code.should.equal('BLAST_FAILED');
  err.message.should.match(/working directory/);
});

test('runMegablast: timeout and abort SIGKILL the child promptly', async function (t) {
  const dir = fs.mkdtempSync(path.join(os.tmpdir(), 'primers-rm-'));
  t.after(function () { fs.rmSync(dir, { recursive: true, force: true }); });
  const tag = String(30 + (process.pid % 1000) / 1000);
  const hang = script(dir, 'hang.sh', 'exec sleep ' + tag);
  let started = Date.now();
  let err = null;
  try { await rm.runMegablast({ bin: hang, args: [], input: '>q0\nACGT\n', timeoutMs: 300, onLine: function () {} }); } catch (e) { err = e; }
  err.code.should.equal('BLAST_TIMEOUT');
  (Date.now() - started).should.be.below(3000);
  await new Promise(function (r) { setTimeout(r, 200); });
  const n = pgrepCount('sleep ' + tag);
  if (n >= 0) n.should.equal(0);

  const ac = new AbortController();
  const reason = new PrimerHttpError(504, 'DEADLINE_EXCEEDED', 'design deadline', {});
  setTimeout(function () { ac.abort(reason); }, 100);
  started = Date.now();
  err = null;
  try { await rm.runMegablast({ bin: hang, args: [], input: '', timeoutMs: 20000, signal: ac.signal, onLine: function () {} }); } catch (e) { err = e; }
  should(err).equal(reason);
  (Date.now() - started).should.be.below(3000);
  await new Promise(function (r) { setTimeout(r, 200); });
  const m = pgrepCount('sleep ' + tag);
  if (m >= 0) m.should.equal(0);
});

test('repeatMask end to end through a fake blastn binary (real spawn, parse, depth, projection)', async function (t) {
  const dir = fs.mkdtempSync(path.join(os.tmpdir(), 'primers-rm-'));
  t.after(function () { fs.rmSync(dir, { recursive: true, force: true }); });
  const fake = script(dir, 'blastn', 'cat > /dev/null\nfor i in 1 2 3; do printf \'q0\\t21\\t80\\t60\\n\'; done');
  const tpl = await buildTemplate({ mode: 'region', system_name: 'syn', region: { region: 'c', start: 101, end: 300, strand: -1 } }, asm());
  const out = await rm.repeatMask(tpl, {}, maskDeps({ cfg: stubs.cfg({ blastn: fake }) }));
  out.mask_source.should.equal('blast_depth');
  out.mask.should.eql([[21, 60]]);
  coords.maskedBases(out.mask).should.equal(60);
});
