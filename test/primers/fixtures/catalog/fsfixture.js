'use strict';

// Builders for throw-away <fasta_root>/<system_name>/ trees used by the resolver tests.
// Real bgzip FASTA (BGZF blocks + .fai + .gzi) is written in pure Node so soft-mask sampling runs
// through the real @gmod/indexedfasta reader without needing bgzip/samtools on the test host.

const fs = require('fs');
const os = require('os');
const path = require('path');
const zlib = require('zlib');

const BGZF_EOF = Buffer.from('1f8b08040000000000ff0600424302001b0003000000000000000000', 'hex');
const BGZF_BLOCK_INPUT = 0xff00;
const LINE_WIDTH = 60;

function makeRoot(tag) {
  return fs.mkdtempSync(path.join(os.tmpdir(), 'primers-' + (tag || 'asm') + '-'));
}

function removeRoot(root) {
  if (root && root.indexOf(os.tmpdir()) === 0) fs.rmSync(root, { recursive: true, force: true });
}

// mtime in epoch seconds
function touch(file, content, mtime) {
  fs.mkdirSync(path.dirname(file), { recursive: true });
  fs.writeFileSync(file, content === undefined ? '' : content);
  if (mtime !== undefined) fs.utimesSync(file, mtime, mtime);
  return file;
}

// entries: [[name, length], ...] -> samtools-style .fai text
function faiText(entries) {
  let offset = 0;
  return entries.map(function (e) {
    offset += e[0].length + 2;
    const line = [e[0], e[1], offset, LINE_WIDTH, LINE_WIDTH + 1].join('\t');
    offset += e[1] + Math.ceil(e[1] / LINE_WIDTH);
    return line;
  }).join('\n') + '\n';
}

// Index-only FASTA (placeholder .fa.gz + .fai + .gzi): enough for candidate scanning and region scoring.
// kind: 'dna' | 'dna_sm' | 'dna_rm'; opts: {fai:false, gzi:false, mtime}
function indexOnlyFasta(dir, prefix, kind, entries, opts) {
  opts = opts || {};
  const file = path.join(dir, 'dna', prefix + '.' + kind + '.toplevel.fa.gz');
  touch(file, 'placeholder', opts.mtime);
  if (opts.fai !== false) touch(file + '.fai', faiText(entries), opts.mtime);
  if (opts.gzi !== false) touch(file + '.gzi', Buffer.alloc(8), opts.mtime);
  return file;
}

function bgzfBlocks(buf) {
  const blocks = [];
  for (let off = 0; off < buf.length; off += BGZF_BLOCK_INPUT) {
    const chunk = buf.subarray(off, Math.min(buf.length, off + BGZF_BLOCK_INPUT));
    const cdata = zlib.deflateRawSync(chunk, { level: 6 });
    const header = Buffer.from([0x1f, 0x8b, 8, 4, 0, 0, 0, 0, 0, 0xff, 6, 0, 0x42, 0x43, 2, 0, 0, 0]);
    const size = header.length + cdata.length + 8;
    if (size > 0x10000) throw new Error('BGZF block too large');
    header.writeUInt16LE(size - 1, 16);
    const trailer = Buffer.alloc(8);
    trailer.writeUInt32LE(zlib.crc32(chunk) >>> 0, 0);
    trailer.writeUInt32LE(chunk.length, 4);
    blocks.push({ data: Buffer.concat([header, cdata, trailer]), usize: chunk.length });
  }
  return blocks;
}

// Real bgzip FASTA with .fai and .gzi. seqs: [{name, seq}]; opts: {fai:false, gzi:false, mtime}
function bgzipFasta(dir, prefix, kind, seqs, opts) {
  opts = opts || {};
  const file = path.join(dir, 'dna', prefix + '.' + kind + '.toplevel.fa.gz');
  let text = '';
  const fai = [];
  seqs.forEach(function (s) {
    text += '>' + s.name + '\n';
    const offset = Buffer.byteLength(text);
    for (let i = 0; i < s.seq.length; i += LINE_WIDTH) text += s.seq.slice(i, i + LINE_WIDTH) + '\n';
    fai.push([s.name, s.seq.length, offset, LINE_WIDTH, LINE_WIDTH + 1].join('\t'));
  });
  const blocks = bgzfBlocks(Buffer.from(text));
  const gziPairs = [];
  let coffset = 0;
  let uoffset = 0;
  blocks.forEach(function (b, i) {
    if (i > 0) gziPairs.push([coffset, uoffset]);
    coffset += b.data.length;
    uoffset += b.usize;
  });
  touch(file, Buffer.concat(blocks.map(function (b) { return b.data; }).concat([BGZF_EOF])), opts.mtime);
  if (opts.fai !== false) touch(file + '.fai', fai.join('\n') + '\n', opts.mtime);
  if (opts.gzi !== false) {
    const gzi = Buffer.alloc(8 + 16 * gziPairs.length);
    gzi.writeBigUInt64LE(BigInt(gziPairs.length), 0);
    gziPairs.forEach(function (p, i) {
      gzi.writeBigUInt64LE(BigInt(p[0]), 8 + 16 * i);
      gzi.writeBigUInt64LE(BigInt(p[1]), 16 + 16 * i);
    });
    touch(file + '.gzi', gzi, opts.mtime);
  }
  return file;
}

// BLAST DB files for base (e.g. prefix + '.dna.toplevel').
// opts: {volumes: n (writes .nal + base.NN.nin/.nsq/.nhr), nsqBytes, nsq:false, mtime}
function blastDb(dir, base, opts) {
  opts = opts || {};
  const bytes = opts.nsqBytes || 100;
  if (opts.volumes) {
    const vols = [];
    for (let i = 0; i < opts.volumes; i++) vols.push(base + '.' + String(i).padStart(2, '0'));
    touch(path.join(dir, base + '.nal'), '#\nTITLE test\nDBLIST ' + vols.join(' ') + '\n', opts.mtime);
    vols.forEach(function (v) {
      touch(path.join(dir, v + '.nin'), 'i', opts.mtime);
      touch(path.join(dir, v + '.nhr'), 'h', opts.mtime);
      touch(path.join(dir, v + '.nsq'), Buffer.alloc(bytes), opts.mtime);
    });
    return;
  }
  touch(path.join(dir, base + '.nin'), 'i', opts.mtime);
  touch(path.join(dir, base + '.nhr'), 'h', opts.mtime);
  if (opts.nsq !== false) touch(path.join(dir, base + '.nsq'), Buffer.alloc(bytes), opts.mtime);
}

// Deterministic ACGT sequence.
function randomSeq(len, seed) {
  let x = (seed || 1) >>> 0;
  const out = new Array(len);
  for (let i = 0; i < len; i++) {
    x = (Math.imul(x, 1664525) + 1013904223) >>> 0;
    out[i] = 'ACGT'[x >>> 30];
  }
  return out.join('');
}

// fs.promises wrapper that records calls; overrides: {readdir, stat, readFile} replacing the real ones.
function spyFs(overrides) {
  overrides = overrides || {};
  const base = fs.promises;
  const calls = [];
  function wrapped(name) {
    return function () {
      const args = Array.prototype.slice.call(arguments);
      calls.push([name, args[0]]);
      return (overrides[name] || base[name]).apply(base, args);
    };
  }
  return { calls: calls, readdir: wrapped('readdir'), stat: wrapped('stat'), readFile: wrapped('readFile') };
}

module.exports = {
  makeRoot,
  removeRoot,
  touch,
  faiText,
  indexOnlyFasta,
  bgzipFasta,
  blastDb,
  randomSeq,
  spyFs
};
