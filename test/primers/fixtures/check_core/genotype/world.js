'use strict';

// The genotyping fake world of run_genotyping.test.js (genotyping spec §7.3): the verdicts fake world (fake blastn over
// planted sequence, real bgzip FASTA re-alignment) with a reference whose chromosome 1 carries the recorded sorghum_bicolor
// 1:9000-15500 at its own coordinates, and pan-genome genomes built from the rs871475760 locus 1:10600-11600.
// Three sets: KASP S1 (reverse) and S2 (forward) of §2.9 and the gel AS-PCR set A1 of §2.10(d) (-2 mismatch on both
// allele-specific primers).
//
//   genotypingWorld(opts) -> {world, body(genomes, extra), genomes: {name: description}, cleanup}
//     opts.referenceAltCopy: the reference also carries the ALT locus on a chromosome 2 (its allele becomes ambiguous)
//   withMegablast(world, opts) -> spawnLines that answers megablast calls from the world's sequences (ungapped, >= 90 %
//     identity, both strands) and passes everything else to world.spawnLines; opts.fail(system_name) -> exit code 2
//   locus(subs), paralog(subs, opts), revcomp

const path = require('path');

const W = require('../../verdicts/fake_world');
const { fakeMongo } = require('../../verdicts/fake_mongo');
const stubs = require('./stubs');

const ROOT = path.resolve(__dirname, '..', '..', '..', '..', '..');
const realign = require(path.join(ROOT, 'api/helpers/primers/check/realign'));
const classify = require(path.join(ROOT, 'api/helpers/primers/check/classify'));

const LOCUS = Object.freeze({ start: 10600, end: 11600 });
const VARIANT = Object.freeze({ region: '1', position: 11109, ref: 'C', alt: 'A' });
// Every primer site of the three sets: A1 common 10880-10903; S1 common 11068-11090, allele-specific 11109-11132;
// S2 allele-specific 11081-11109, common 11149-11172.
const SITES = Object.freeze([[10870, 10910], [11060, 11180]]);
// The 3' and -2 bases of the three common primers (S1 and A1 face forward, S2 reverse).
const COMMON_3P = Object.freeze([11090, 11149, 10903]);
const COMMON_2 = Object.freeze([11089, 11150, 10902]);
const COMMONS = Object.freeze([
  { primer: 'AGCTTCTCTAAGTGGTTATCCGA', face: 'F', p5: 11068, p3: 11090 },
  { primer: 'TCTTTGTCTACTGAGAAATCCAGA', face: 'R', p5: 11172, p3: 11149 },
  { primer: 'TGCATCAACAAATGTGCTATGTGT', face: 'F', p5: 10880, p3: 10903 }
]);
const SUB = { A: 'C', C: 'G', G: 'T', T: 'A' };

function other(p) {
  return stubs.bases(p, p) === 'C' ? 'G' : 'C';
}

function subsAt(positions) {
  const s = {};
  positions.forEach(function (p) { s[p] = other(p); });
  return s;
}

// Substitutions at each common primer's -1 and -2 bases that realign.js cannot absorb with a gap, so that its site is
// 3'-blocked (classify.isBlocked); the first pair of bases, in ACGT order, that does it.
function blockedCommonSubs() {
  const out = {};
  COMMONS.forEach(function (c) {
    const d1 = c.p3;
    const d2 = c.face === 'F' ? c.p3 - 1 : c.p3 + 1;
    const lo = Math.min(c.p5, c.p3) - realign.PAD;
    const hi = Math.max(c.p5, c.p3) + realign.PAD;
    const found = 'ACGT'.split('').some(function (b1) {
      return 'ACGT'.split('').some(function (b2) {
        if (b1 === stubs.bases(d1, d1) || b2 === stubs.bases(d2, d2)) return false;
        const win = stubs.bases(lo, hi).split('');
        win[d1 - lo] = b1;
        win[d2 - lo] = b2;
        const aln = realign.realignSiteOnWindow({ key: 'c', primer: c.primer, face: c.face, p3: c.p3 }, win.join(''), lo, 5);
        if (!classify.isBlocked(aln)) return false;
        out[d1] = b1;
        out[d2] = b2;
        return true;
      });
    });
    if (!found) throw new Error('genotype world: no blocking substitution for ' + c.primer);
  });
  return out;
}

// The reference locus with substitutions {position: base}.
function locus(subs) {
  const a = stubs.bases(LOCUS.start, LOCUS.end).split('');
  Object.keys(subs || {}).forEach(function (p) { a[Number(p) - LOCUS.start] = subs[p]; });
  return a.join('');
}

// A diverged copy: every 4th base outside the primer sites substituted (about 86 % identity over A1's 353 bp segment),
// then subs applied.
function paralog(subs) {
  const a = stubs.bases(LOCUS.start, LOCUS.end).split('');
  for (let p = LOCUS.start, k = 0; p <= LOCUS.end; p++) {
    if (SITES.some(function (r) { return p >= r[0] && p <= r[1]; })) continue;
    if (k++ % 4 === 0) a[p - LOCUS.start] = SUB[a[p - LOCUS.start]];
  }
  Object.keys(subs || {}).forEach(function (p) { a[Number(p) - LOCUS.start] = subs[p]; });
  return a.join('');
}

function chr(name, parts, seed) {
  return { name: name, seq: W.randomSeq(2000, seed) + parts + W.randomSeq(2000, seed + 1) };
}

const PAIRS = Object.freeze([
  { id: 'S1_REF', left: 'AGCTTCTCTAAGTGGTTATCCGA', right: 'ATCTTTGACTAGCGAGAAATTCAG', expected: { region: '1', start: 11068, end: 11132 } },
  { id: 'S1_ALT', left: 'AGCTTCTCTAAGTGGTTATCCGA', right: 'ATCTTTGACTAGCGAGAAATTCAT', expected: { region: '1', start: 11068, end: 11132 } },
  { id: 'S2_REF', left: 'GGTTATCCGAATATAGTCATACTCTATTC', right: 'TCTTTGTCTACTGAGAAATCCAGA', expected: { region: '1', start: 11081, end: 11172 } },
  { id: 'S2_ALT', left: 'GGTTATCCGAATATAGTCATACTCTATTA', right: 'TCTTTGTCTACTGAGAAATCCAGA', expected: { region: '1', start: 11081, end: 11172 } },
  { id: 'A1_REF', left: 'TGCATCAACAAATGTGCTATGTGT', right: 'ATCTTTGACTAGCGAGAAATTCGG', expected: { region: '1', start: 10880, end: 11132 } },
  { id: 'A1_ALT', left: 'TGCATCAACAAATGTGCTATGTGT', right: 'ATCTTTGACTAGCGAGAAATTCGT', expected: { region: '1', start: 10880, end: 11132 } }
]);

const SETS = Object.freeze([
  { id: 'S1', ref_pair: 'S1_REF', alt_pair: 'S1_ALT' },
  { id: 'S2', ref_pair: 'S2_REF', alt_pair: 'S2_ALT' },
  { id: 'A1', ref_pair: 'A1_REF', alt_pair: 'A1_ALT' }
]);

// Pan-genome genomes of the §7.3 world, in the order of its expectations, plus the fallback cases.
const GENOMES = Object.freeze({
  g01_ref: 'REF locus',
  g02_alt: 'ALT locus',
  g03_other: 'third allele T',
  g04_missing: 'no locus',
  g05_dup_alt: 'ALT locus twice (chromosome 2 inverted)',
  g06_paralog85: 'REF locus + an 86 % paralog whose common sites are blocked',
  g07_alt_paralog: 'REF locus + an amplifying ALT paralog',
  g08_ref_paralog: 'REF locus + an amplifying REF paralog',
  g09_common3p: "REF locus with a private SNP under each common primer's 3' base",
  g10_common2: "REF locus with a private SNP at each common primer's -2 base",
  m1_mb_alt: 'ALT locus whose allele-specific sites carry 6 mismatches each (no product; megablast finds the locus)',
  m2_mb_fail: 'ALT locus without products whose megablast fails',
  m3_nodb: 'no BLAST database',
  m4_blast_error: 'every BLAST fails'
});

// 6 substitutions in each allele-specific site (S1/A1 11109-11132, S2 11081-11109), none in the core 11108-11110.
const NO_AS_SITES = Object.freeze([11112, 11115, 11118, 11121, 11124, 11127, 11083, 11086, 11092, 11095, 11098, 11101]);

function genotypingWorld(opts) {
  const o = opts || {};
  const alt = { 11109: 'A' };
  const noSites = Object.assign(subsAt(NO_AS_SITES), alt);
  const refChromosomes = [{ name: '1', seq: W.randomSeq(8999, 501) + stubs.bases(9000, 15500) }];
  if (o.referenceAltCopy) refChromosomes.push(chr('2', locus(alt), 503));
  const spec = {
    ref: { taxon_id: 1001, map_id: 'MAP_REF', display_name: 'Reference', chromosomes: refChromosomes },
    g01_ref: { chromosomes: [chr('1', locus(), 601)] },
    g02_alt: { chromosomes: [chr('1', locus(alt), 611)] },
    g03_other: { chromosomes: [chr('1', locus({ 11109: 'T' }), 621)] },
    g04_missing: { chromosomes: [{ name: '1', seq: W.randomSeq(5000, 631) }] },
    g05_dup_alt: { chromosomes: [chr('1', locus(alt), 641), chr('2', W.revcomp(locus(alt)), 643)] },
    g06_paralog85: { chromosomes: [chr('1', locus(), 651), chr('2', paralog(Object.assign(blockedCommonSubs(), alt)), 653)] },
    g07_alt_paralog: { chromosomes: [chr('1', locus(), 661), chr('2', paralog(alt), 663)] },
    g08_ref_paralog: { chromosomes: [chr('1', locus(), 671), chr('2', paralog(), 673)] },
    g09_common3p: { chromosomes: [chr('1', locus(subsAt(COMMON_3P)), 681)] },
    g10_common2: { chromosomes: [chr('1', locus(subsAt(COMMON_2)), 691)] },
    m1_mb_alt: { chromosomes: [chr('1', locus(noSites), 701)] },
    m2_mb_fail: { chromosomes: [chr('1', locus(noSites), 711)] },
    m3_nodb: { chromosomes: [chr('1', locus(), 721)], dnaDb: false },
    m4_blast_error: { chromosomes: [chr('1', locus(), 731)], failures: Infinity }
  };
  const world = W.buildWorld({ genomes: spec });
  // No gene annotation in region mode; W.worldCtx still wants a mongo handle.
  world.genes = [];
  world.mongo = function (opts) { return fakeMongo(world.genes, opts); };
  function body(genomes, extra) {
    return Object.assign({
      system_name: 'ref',
      mode: 'region',
      checks: ['pangenome'],
      genomes: genomes,
      pairs: PAIRS.map(function (p) { return Object.assign({}, p, { expected: Object.assign({}, p.expected) }); }),
      genotyping: { variant: Object.assign({}, VARIANT), sets: SETS.map(function (s) { return Object.assign({}, s); }) }
    }, extra || {});
  }
  return { world: world, body: body, genomes: GENOMES, cleanup: world.cleanup };
}

// Megablast answered from the world's own sequences (the verdicts fake blastn only knows primer queries).
function withMegablast(world, opts) {
  const o = opts || {};
  const calls = [];
  const spawnLines = async function (cmd, args, so) {
    if (args[0] !== '-task' || args[1] !== 'megablast') return world.spawnLines(cmd, args, so);
    const db = args[args.indexOf('-db') + 1];
    const sys = Object.keys(world.assemblies).find(function (k) { return world.assemblies[k].blastdb.dna === db; });
    calls.push({ system_name: sys, args: args.slice(), cwd: so && so.cwd });
    await new Promise(function (resolve) { setImmediate(resolve); });
    if (o.fail && o.fail(sys)) return { code: 2, signal: null, stderrTail: 'BLAST Database error: fake failure\n' };
    const entry = world.dbs.get(db);
    const query = String(so.stdin).split('\n')[1];
    const n = query.length;
    const max = Math.floor(n / 10);
    const emit = function (line) { if (so.onLine) so.onLine(line); };
    (entry ? entry.seqs : []).forEach(function (subj) {
      [1, -1].forEach(function (strand) {
        const q = strand === 1 ? query : W.revcomp(query);
        for (let i = 0; i + n <= subj.seq.length; i++) {
          let mm = 0;
          for (let k = 0; k < n && mm <= max; k++) if (q.charCodeAt(k) !== subj.seq.charCodeAt(i + k)) mm++;
          if (mm > max) continue;
          const piece = subj.seq.substr(i, n);
          emit([subj.name, strand === 1 ? i + 1 : i + n, strand === 1 ? i + n : i + 1, strand === 1 ? 'plus' : 'minus',
            (100 * (n - mm) / n).toFixed(3), n, (1.8 * (n - mm)).toFixed(0), 1, n, query, strand === 1 ? piece : W.revcomp(piece)].join('\t'));
        }
      });
    });
    return { code: 0, signal: null, stderrTail: '' };
  };
  spawnLines.calls = calls;
  return spawnLines;
}

module.exports = { LOCUS, VARIANT, SITES, COMMON_3P, COMMON_2, COMMONS, PAIRS, SETS, GENOMES, NO_AS_SITES, locus, paralog, blockedCommonSubs, genotypingWorld, withMegablast };
