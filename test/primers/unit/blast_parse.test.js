'use strict';

const { describe, it } = require('node:test');
const assert = require('node:assert');
const fs = require('fs');
const os = require('os');
const path = require('path');
const should = require('should');

const blast = require('../../../api/helpers/primers/check/blast');

const FIX = path.join(__dirname, '..', 'fixtures', 'check_core');
const readLines = (name) => fs.readFileSync(path.join(FIX, name), 'utf8').split('\n').filter(Boolean);
const P2_L = 'GGACAGCTCCACAACATATCAG';
const P2_R = 'GGACATTTGAAGCCCATGGCC';
const QUERY = [P2_L, P2_R, 'GGACAGATCCACATCATATC', 'GGACAGCTCCACAACATTCAG'];
const NODE = process.execPath;

describe('blast.buildArgs', () => {
  const db = '/scratch/olson/fasta/sorghum_bicolor/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel';

  it('builds exactly the §B.3 genome argument array', () => {
    should(blast.buildArgs({ target: 'genome', wordSize: 5, threads: 4, db, maxTargetSeqs: 867 })).eql([
      '-task', 'blastn-short', '-reward', '1', '-penalty', '-1', '-word_size', '5', '-ungapped',
      '-evalue', '30000', '-searchsp', '15000000000', '-dust', 'no', '-soft_masking', 'false',
      '-max_target_seqs', '5000', '-max_hsps', '100000',
      '-num_threads', '4', '-db', db, '-query', '-',
      '-outfmt', '6 qseqid sseqid qlen qstart qend sstart send sstrand mismatch']);
  });

  it('builds the cDNA variant with qseq sseq and the DB count', () => {
    const cdna = '/scratch/olson/fasta/sorghum_bicolor/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.cdna.all';
    const args = blast.buildArgs({ target: 'cdna', wordSize: 6, threads: 1, db: cdna, maxTargetSeqs: 47110 });
    should(args[args.indexOf('-max_target_seqs') + 1]).equal('47110');
    should(args[args.indexOf('-word_size') + 1]).equal('6');
    should(args[args.indexOf('-num_threads') + 1]).equal('1');
    should(args[args.length - 1]).equal('6 qseqid sseqid qlen qstart qend sstart send sstrand mismatch qseq sseq');
    should(blast.targetFromArgs(args)).equal('cdna');
  });

  it('uses max(5000, count) and strips .nal/.nin', () => {
    const a = blast.buildArgs({ target: 'genome', wordSize: 6, threads: 1, db: db + '.nal' });
    should(a[a.indexOf('-max_target_seqs') + 1]).equal('5000');
    should(a[a.indexOf('-db') + 1]).equal(db);
    should(blast.targetFromArgs(a)).equal('genome');
  });

  it('rejects unsafe or invalid options', () => {
    should(() => blast.buildArgs({ target: 'protein', wordSize: 5, threads: 1, db })).throw(TypeError);
    should(() => blast.buildArgs({ target: 'genome', wordSize: 3, threads: 1, db })).throw(TypeError);
    should(() => blast.buildArgs({ target: 'genome', wordSize: 5, threads: 0, db })).throw(TypeError);
    should(() => blast.buildArgs({ target: 'genome', wordSize: 5, threads: 1, db: '' })).throw(TypeError);
    should(() => blast.buildArgs({ target: 'genome', wordSize: 5, threads: 1, db: '/a b' })).throw(TypeError);
    should(() => blast.buildArgs({ target: 'genome', wordSize: 5, threads: 1, db: '-remote' })).throw(TypeError);
  });

  it('describes the engine as in §B.12', () => {
    should(blast.engineDescription(5)).equal('blastn-short r1 p-1 ws5 ungapped e30000 searchsp1.5e10');
  });
});

describe('blast.queryFasta and uniquePrimers', () => {
  it('writes >q<i> records of uppercase primers', () => {
    should(blast.queryFasta(['acgtacgtacgtacgt', P2_R])).equal('>q0\nACGTACGTACGTACGT\n>q1\n' + P2_R + '\n');
  });

  it('rejects non-ACGT input, injection and duplicates', () => {
    should(() => blast.queryFasta([])).throw(TypeError);
    should(() => blast.queryFasta(['ACGTN'])).throw(TypeError);
    should(() => blast.queryFasta(['ACGT\n>q9\nAAAA'])).throw(TypeError);
    should(() => blast.queryFasta(['ACGT', 'acgt'])).throw(/duplicate/);
    should(() => blast.queryFasta([42])).throw(TypeError);
  });

  it('dedupes case-insensitively in first-seen order', () => {
    should(blast.uniquePrimers(['acgt', P2_L, 'ACGT', P2_R, P2_L.toLowerCase()])).eql(['ACGT', P2_L, P2_R]);
  });
});

describe('blast.parseLine', () => {
  it('parses a plus-strand genome hit', () => {
    should(blast.parseLine('q0\t4\t22\t1\t22\t7423537\t7423558\tplus\t0', 'genome')).eql({
      q: 0, qseqid: 'q0', sseqid: '4', qlen: 22, qstart: 1, qend: 22, sstart: 7423537, send: 7423558, strand: 1, mismatch: 0
    });
  });

  it('parses a minus-strand genome hit (sstart > send) and CRLF', () => {
    const h = blast.parseLine('q0\t5\t22\t1\t22\t66890824\t66890803\tminus\t2\r', 'genome');
    should(h.strand).equal(-1);
    should(h.sstart).equal(66890824);
    should(h.send).equal(66890803);
    should(h.mismatch).equal(2);
  });

  it('parses cDNA hits with qseq and sseq', () => {
    const line = readLines('p2_cdna.tsv').find((l) => l.startsWith('q1\tSORBI_3005G183900.1\t'));
    const h = blast.parseLine(line, 'cdna');
    should(h).match({ q: 1, sseqid: 'SORBI_3005G183900.1', strand: -1, mismatch: 2, qseq: P2_R, sseq: 'GGACATCTGAAACCCATGGCC' });
  });

  it('returns null for blank and comment lines', () => {
    should(blast.parseLine('', 'genome')).be.null();
    should(blast.parseLine('# BLASTN 2.13.0+', 'genome')).be.null();
    should(blast.parseLine(undefined, 'genome')).be.null();
  });

  it('throws BLAST_PARSE on malformed lines', () => {
    const bad = [
      'q0\t4\t22\t1\t22\t7423537\t7423558\tplus',
      'q0\t4\t22\t1\t22\t7423537\t7423558\tplus\t0\tEXTRA',
      'x0\t4\t22\t1\t22\t7423537\t7423558\tplus\t0',
      'q0\t\t22\t1\t22\t7423537\t7423558\tplus\t0',
      'q0\t4\t22\t1\t22\t7423537\t7423558\tboth\t0',
      'q0\t4\t22\t1\t2x\t7423537\t7423558\tplus\t0',
      'q0\t4\t22\t5\t3\t7423537\t7423558\tplus\t0',
      'q0\t4\t22\t1\t23\t7423537\t7423558\tplus\t0',
      'q0\t4\t22\t1\t22\t7423558\t7423537\tplus\t0',
      'q0\t4\t22\t1\t22\t7423537\t7423558\tminus\t0',
      'q0\t4\t22\t1\t22\t7423537\t7423558\tplus\t-1'
    ];
    for (const line of bad) should(() => blast.parseLine(line, 'genome')).throw({ code: 'BLAST_PARSE' });
    // A genome-format line is malformed for a cDNA run.
    should(() => blast.parseLine('q0\t4\t22\t1\t22\t7423537\t7423558\tplus\t0', 'cdna')).throw({ code: 'BLAST_PARSE' });
  });

  it('parses every captured real line', () => {
    const lines = readLines('p2_genome_sample.tsv').concat(readLines('p2_genome_loci.tsv'));
    for (const line of lines) {
      const h = blast.parseLine(line, 'genome');
      should(h.q).be.within(0, 3);
      should(h.qend - h.qstart).equal(Math.abs(h.send - h.sstart));
    }
  });
});

describe('blast.parseDbInfo', () => {
  it('reads sequence count and total bases from blastdbcmd -info', () => {
    const info = blast.parseDbInfo(fs.readFileSync(path.join(FIX, 'blastdbcmd_info_cdna.txt'), 'utf8'));
    should(info.num_sequences).equal(47110);
    should(info.total_bases).equal(103922012);
    should(info.title).equal('Sorghum bicolor ssp. bicolor BTx623 cdna');
    should(info.volumes).eql(['/scratch/olson/fasta/sorghum_bicolor/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.cdna.all']);
  });

  it('runs blastdbcmd through an injected spawnLines', async () => {
    const text = fs.readFileSync(path.join(FIX, 'blastdbcmd_info_cdna.txt'), 'utf8');
    let seen;
    const spawnLines = async (cmd, args, opts) => {
      seen = { cmd, args };
      text.split('\n').forEach((l) => opts.onLine(l));
      return { code: 0, signal: null, stderrTail: '' };
    };
    const info = await blast.dbInfo({ cmd: '/home/olson/bin/blastdbcmd', db: '/x/y.cdna.all.nal', spawnLines });
    should(seen).eql({ cmd: '/home/olson/bin/blastdbcmd', args: ['-info', '-db', '/x/y.cdna.all'] });
    should(info.num_sequences).equal(47110);
    await assert.rejects(blast.dbInfo({ cmd: 'x', db: 'y', spawnLines: async () => ({ code: 2, stderrTail: 'no db' }) }),
      { code: 'BLASTDBCMD_FAILED' });
  });
});

describe('blast.runBlast with an injected spawnLines', () => {
  const lines = readLines('p2_genome_loci.tsv');
  const args = blast.buildArgs({ target: 'genome', wordSize: 5, threads: 4, db: '/db/x.dna.toplevel', maxTargetSeqs: 867 });

  function fakeSpawn(outLines, result, capture) {
    return async (cmd, a, opts) => {
      if (capture) Object.assign(capture, { cmd, args: a, stdin: opts.stdin, timeoutMs: opts.timeoutMs, signal: opts.signal });
      for (const l of outLines) opts.onLine(l);
      return Object.assign({ code: 0, signal: null, stderrTail: '', timedOut: false, aborted: false }, result);
    };
  }

  it('passes cwd through to spawnLines for blastn and blastdbcmd', async () => {
    const seen = [];
    const spawnLines = async (cmd, a, opts) => {
      seen.push(opts.cwd);
      if (a[0] === '-info') opts.onLine('\t47,110 sequences; 103,922,012 total bases');
      return { code: 0, signal: null, stderrTail: '' };
    };
    await blast.runBlast({ cmd: 'b', args, primers: QUERY, spawnLines, cwd: '/private/job' });
    await blast.dbInfo({ cmd: 'c', db: '/x/y.cdna.all', spawnLines, cwd: '/private/job' });
    should(seen).eql(['/private/job', '/private/job']);
    should(blast.CHILD_ENV).eql({ PATH: '/usr/bin:/bin', NCBI_DONT_USE_NCBIRC: '1' });
  });

  it('streams parsed hits with the primer attached', async () => {
    const capture = {};
    const hits = [];
    const res = await blast.runBlast({ cmd: '/home/olson/bin/blastn', args, primers: QUERY, spawnLines: fakeSpawn(lines, {}, capture), onHit: (h) => hits.push(h) });
    should(res.hits).equal(lines.length);
    should(hits.length).equal(lines.length);
    should(capture.stdin).equal(blast.queryFasta(QUERY));
    should(capture.timeoutMs).equal(600000);
    should(capture.cmd).equal('/home/olson/bin/blastn');
    for (const h of hits) should(h.primer).equal(QUERY[h.q]);
    should(hits.filter((h) => h.sseqid === '5' && h.primer === P2_L)).match([{ strand: -1, mismatch: 2, sstart: 66890824 }]);
  });

  it('ignores blank lines but fails on unparseable output', async () => {
    const ok = await blast.runBlast({ args, cmd: 'b', primers: QUERY, spawnLines: fakeSpawn(['', lines[0], '']) });
    should(ok.hits).equal(1);
    await assert.rejects(blast.runBlast({ args, cmd: 'b', primers: QUERY, spawnLines: fakeSpawn([lines[0], 'garbage']) }),
      (e) => e.code === 'BLAST_OUTPUT_INVALID' && e.badLines === 1 && e.firstBadLine === 'garbage');
    await assert.rejects(blast.runBlast({ args, cmd: 'b', primers: [P2_L], spawnLines: fakeSpawn(['q3\t4\t21\t1\t17\t7423537\t7423553\tplus\t0']) }),
      { code: 'BLAST_OUTPUT_INVALID' });
  });

  it('maps process outcomes to error codes', async () => {
    await assert.rejects(blast.runBlast({ args, cmd: 'b', primers: QUERY, spawnLines: fakeSpawn([], { code: 2, stderrTail: 'BLAST Database error' }) }),
      (e) => e.code === 'BLAST_FAILED' && e.stderrTail === 'BLAST Database error' && e.exitCode === 2);
    await assert.rejects(blast.runBlast({ args, cmd: 'b', primers: QUERY, spawnLines: fakeSpawn([], { code: null, signal: 'SIGTERM', timedOut: true }) }),
      { code: 'BLAST_TIMEOUT' });
    await assert.rejects(blast.runBlast({ args, cmd: 'b', primers: QUERY, spawnLines: fakeSpawn([], { code: null, signal: 'SIGTERM', aborted: true }) }),
      { code: 'ABORTED' });
    const ac = new AbortController();
    ac.abort();
    await assert.rejects(blast.runBlast({ args, cmd: 'b', primers: QUERY, signal: ac.signal, spawnLines: fakeSpawn([], { code: null, signal: 'SIGTERM' }) }),
      { code: 'ABORTED' });
  });

  it('propagates onHit errors and cDNA parsing', async () => {
    await assert.rejects(blast.runBlast({ args, cmd: 'b', primers: QUERY, spawnLines: fakeSpawn(lines), onHit: () => { throw new Error('boom'); } }), /boom/);
    const cdnaArgs = blast.buildArgs({ target: 'cdna', wordSize: 5, threads: 1, db: '/db/c', maxTargetSeqs: 47110 });
    const hits = [];
    await blast.runBlast({ args: cdnaArgs, cmd: 'b', primers: [P2_L, P2_R], spawnLines: fakeSpawn(readLines('p2_cdna.tsv')), onHit: (h) => hits.push(h) });
    should(hits.length).equal(12);
    should(hits.every((h) => h.qseq.length === h.sseq.length)).be.true();
  });
});

describe('blast.spawnLinesLocal (real child processes)', () => {
  const script = (src) => ['-e', src];

  it('splits lines across chunk boundaries and flushes the final partial line', async () => {
    const lines = [];
    const res = await blast.spawnLinesLocal(NODE, script(
      "const w=(s,d)=>setTimeout(()=>process.stdout.write(s),d);w('ab',0);w('c\\nde',30);w('f\\r\\n\\nlast',60);"), {
      onLine: (l) => lines.push(l), nice: false
    });
    should(res.code).equal(0);
    should(lines).eql(['abc', 'def', '', 'last']);
  });

  it('feeds stdin and keeps only a 2 KB stderr tail', async () => {
    const lines = [];
    const res = await blast.spawnLinesLocal(NODE, script(
      "let d='';process.stdin.on('data',c=>d+=c);process.stdin.on('end',()=>{process.stderr.write('x'.repeat(5000)+'END');" +
      "d.split('\\n').filter(Boolean).forEach(l=>console.log(l.toLowerCase()));process.exit(3);});"), {
      stdin: '>q0\nACGT\n', onLine: (l) => lines.push(l), nice: false
    });
    should(lines).eql(['>q0', 'acgt']);
    should(res.code).equal(3);
    should(res.stderrTail.length).equal(2048);
    should(res.stderrTail.endsWith('xEND')).be.true();
    should(res.timedOut).be.false();
    should(res.aborted).be.false();
  });

  it('runs under nice -n 10 by default with a minimal environment that disables .ncbirc', async () => {
    const lines = [];
    const base = os.getPriority();
    await blast.spawnLinesLocal(NODE, script("console.log(require('os').getPriority());console.log(JSON.stringify(process.env))"), {
      onLine: (l) => lines.push(l)
    });
    should(Number(lines[0])).equal(Math.min(19, base + 10));
    should(JSON.parse(lines[1])).eql({ PATH: '/usr/bin:/bin', NCBI_DONT_USE_NCBIRC: '1' });
  });

  it('runs in the given cwd (sec-ncbirc-cwd-tmp)', async () => {
    const dir = fs.mkdtempSync(path.join(os.tmpdir(), 'primers-cwd-'));
    try {
      const lines = [];
      await blast.spawnLinesLocal(NODE, script('console.log(process.cwd())'), { cwd: dir, nice: false, onLine: (l) => lines.push(l) });
      should(lines).eql([fs.realpathSync(dir)]);
    } finally {
      fs.rmSync(dir, { recursive: true, force: true });
    }
  });

  it('times out with SIGTERM and resolves timedOut', async () => {
    const t0 = Date.now();
    const res = await blast.spawnLinesLocal(NODE, script('setTimeout(()=>{},10000)'), { timeoutMs: 150, nice: false });
    should(res.timedOut).be.true();
    should(res.signal).equal('SIGTERM');
    should(Date.now() - t0).be.below(3000);
  });

  it('escalates to SIGKILL when SIGTERM is ignored', async () => {
    const t0 = Date.now();
    const res = await blast.spawnLinesLocal(NODE, script("process.on('SIGTERM',()=>{});console.log('ready');setInterval(()=>{},1000)"), {
      timeoutMs: 300, killGraceMs: 300, nice: false, onLine: () => {}
    });
    should(res.timedOut).be.true();
    should(res.signal).equal('SIGKILL');
    should(Date.now() - t0).be.below(3000);
    should(blast.KILL_GRACE_MS).equal(3000);
  });

  it('aborts via AbortSignal, and does not spawn when already aborted', async () => {
    const ac = new AbortController();
    const p = blast.spawnLinesLocal(NODE, script("console.log('go');setTimeout(()=>{},10000)"), {
      signal: ac.signal, nice: false, onLine: () => ac.abort()
    });
    const res = await p;
    should(res.aborted).be.true();
    should(res.signal).equal('SIGTERM');
    const pre = new AbortController();
    pre.abort();
    const r2 = await blast.spawnLinesLocal('/nonexistent/binary', [], { signal: pre.signal });
    should(r2).match({ aborted: true, code: null });
  });

  it('kills the child and rejects when onLine throws', async () => {
    const t0 = Date.now();
    await assert.rejects(blast.spawnLinesLocal(NODE, script("setInterval(()=>console.log('line'),5)"), {
      nice: false, onLine: () => { throw new Error('consumer failed'); }
    }), /consumer failed/);
    should(Date.now() - t0).be.below(3000);
  });

  it('rejects when the executable cannot be spawned', async () => {
    await assert.rejects(blast.spawnLinesLocal('/nonexistent/binary', [], { nice: false }), { code: 'ENOENT' });
    // Through nice the exec failure is an exit code instead.
    const res = await blast.spawnLinesLocal('/nonexistent/binary', []);
    should(res.code).equal(127);
    should(res.stderrTail).match(/nonexistent/);
  });

  it('drives runBlast end to end with a stand-in blastn', async () => {
    const dir = fs.mkdtempSync(path.join(os.tmpdir(), 'primers-blast-'));
    try {
      const fake = path.join(dir, 'fake_blastn.js');
      const out = readLines('p2_genome_loci.tsv').join('\n') + '\n';
      fs.writeFileSync(fake, "#!" + NODE + "\n'use strict';let d='';process.stdin.on('data',c=>d+=c);process.stdin.on('end',()=>{" +
        "if(!d.startsWith('>q0\\n')){process.stderr.write('bad stdin');process.exit(1);}process.stdout.write(" + JSON.stringify(out) + ");});\n");
      fs.chmodSync(fake, 0o755);
      const hits = [];
      const res = await blast.runBlast({ cmd: fake, args: ['-outfmt', blast.GENOME_OUTFMT], primers: QUERY, onHit: (h) => hits.push(h) });
      should(res.hits).equal(12);
      should(hits.length).equal(12);
    } finally {
      fs.rmSync(dir, { recursive: true, force: true });
    }
  });
});
