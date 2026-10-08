'use strict';

const { describe, it } = require('node:test');
const assert = require('node:assert');
const fs = require('fs');
const path = require('path');
const should = require('should');

const blast = require('../../../api/helpers/primers/check/blast');
const sites = require('../../../api/helpers/primers/check/sites');
const realign = require('../../../api/helpers/primers/check/realign');

const FIX = path.join(__dirname, '..', 'fixtures', 'check_core');
const DOC = JSON.parse(fs.readFileSync(path.join(FIX, 'sorghum_windows.json'), 'utf8'));
const readLines = (name) => fs.readFileSync(path.join(FIX, name), 'utf8').split('\n').filter(Boolean);
const P = DOC.primers;
const QUERY = DOC.query_order.map((k) => P[k]);

function windowFetch(windows) {
  const calls = [];
  const fn = async (fastaPath, region, start1, end1, strand) => {
    calls.push({ fastaPath, region, start1, end1, strand });
    for (const w of windows) {
      if (w.subject === region && start1 >= w.start && end1 <= w.end) return w.seq.slice(start1 - w.start, end1 - w.start + 1);
    }
    const e = new Error('outside fixture windows');
    e.code = 'REGION_OUT_OF_BOUNDS';
    throw e;
  };
  fn.calls = calls;
  return fn;
}

const bicolorFetch = () => windowFetch(Object.values(DOC.windows));
const locusSite = (q, subject) => {
  const line = readLines('p2_genome_loci.tsv').find((l) => l.startsWith('q' + q + '\t' + subject + '\t'));
  const h = blast.parseLine(line, 'genome');
  h.primer = QUERY[h.q];
  return sites.hitToSite(h);
};
const plainSeq = (subject, start, end) => {
  for (const w of Object.values(DOC.windows)) {
    if (w.subject === subject && start >= w.start && end <= w.end) return w.seq.slice(start - w.start, end - w.start + 1);
  }
  throw new Error('no window');
};

describe('realign.revcomp', () => {
  it('complements IUPAC codes and preserves case', () => {
    should(realign.revcomp('ACGTacgt')).equal('acgtACGT');
    should(realign.revcomp('RYKMSWBDHVN')).equal('NBDHVWSKMRY');
    should(realign.revcomp('ryn')).equal('nry');
    should(realign.revcomp('A-X')).equal('N-T');
  });
});

describe('realign.alignPrimer (DP rules)', () => {
  const primer = 'ACGTACGTAC';

  it('aligns a perfect site with free window ends', () => {
    const a = realign.alignPrimer(primer, 'TTT' + primer + 'GGG', 13, 5);
    should(a).match({ mm: 0, mm_3p: 0, terminal_mm: false, gaps: 0, mm_pos: [], startIdx: 4, endIdx: 13, ops: 'MMMMMMMMMM' });
  });

  it('reports a terminal mismatch at distance 1', () => {
    const a = realign.alignPrimer(primer, 'TTTACGTACGTAAGGG', 13, 5);
    should(a).match({ mm: 1, mm_3p: 1, terminal_mm: true, gaps: 0, mm_pos: [1], endIdx: 13 });
  });

  it('never ends with a gap: the 3\' base is forced onto a window base', () => {
    // With a terminal gap allowed this would be ACGTACGTAC + T/- ; instead T pairs with G.
    const a = realign.alignPrimer('ACGTACGTACT', 'TTTACGTACGTACGGG', 14, 5);
    should(a).match({ mm: 1, terminal_mm: true, gaps: 0, mm_pos: [1], endIdx: 14 });
    should(a.ops.endsWith('X')).be.true();
  });

  it('counts a window base skipped inside the primer as one gap', () => {
    const a = realign.alignPrimer(primer, 'TTTACGTAGCGTACGGG', 14, 5);
    should(a).match({ mm: 1, gaps: 1, mm_pos: [6], startIdx: 4, endIdx: 14, ops: 'MMMMMDMMMMM' });
    should(a.terminal_mm).be.false();
  });

  it('counts an extra primer base as one gap at its own position', () => {
    const a = realign.alignPrimer('ACGTATCGTAC', 'TTTACGTACGTACGGG', 13, 5);
    should(a).match({ mm: 1, gaps: 1, mm_pos: [6], startIdx: 4, endIdx: 13, ops: 'MMMMMIMMMMM' });
  });

  it('places equivalent gaps toward the 5\' end', () => {
    const a = realign.alignPrimer('GACTTTTCAG', 'CCCGACTTTTTCAGCCC', 14, 5);
    should(a).match({ mm: 1, gaps: 1, mm_pos: [8], startIdx: 4, endIdx: 14, ops: 'MMMDMMMMMMM' });
  });

  it('breaks end ties toward the expected 3\' end, then the lower index', () => {
    should(realign.alignPrimer('AAAAA', 'AAAAAAAAAAA', 8, 5)).match({ mm: 0, startIdx: 4, endIdx: 8 });
    should(realign.alignPrimer('AAAAA', 'AAAAAAAAAAA', 10, 5)).match({ endIdx: 10 });
    should(realign.alignPrimer('AAAAA', 'AAAAAAAAAAA', 99, 5)).match({ endIdx: 11 });
    should(realign.alignPrimer('AAAAA', 'AAAAAAAAAAA', undefined, 5)).match({ endIdx: 8 }); // m-3
    // Two perfect ends (7 and 12) equally far from 9.5: the lower index wins.
    should(realign.alignPrimer('AAAAC', 'GGAAAACAAAACGG', 9.5, 5)).match({ mm: 0, endIdx: 7 });
  });

  it('is case-insensitive and treats N as a mismatch', () => {
    should(realign.alignPrimer(primer, 'ttt' + primer.toLowerCase() + 'ggg', 13, 5)).match({ mm: 0 });
    // N at window index 8 pairs with primer base 5 → distance 10-5+1 = 6.
    should(realign.alignPrimer(primer, 'TTTACGTNCGTACGGG', 13, 5)).match({ mm: 1, mm_pos: [6], ops: 'MMMMXMMMMM' });
  });

  it('derives mm_3p from the window parameter', () => {
    const a = realign.alignPrimer(primer, 'TTTACGTAGGTCCGGG', 13, 5); // mismatches at primer 9? and 6
    should(a.mm_pos.length).equal(a.mm);
    should(realign.alignPrimer(primer, 'TTTACGTAGGTCCGGG', 13, 3).mm_3p).equal(a.mm_pos.filter((d) => d <= 3).length);
  });
});

describe('realign on real sorghum windows', () => {
  it('P2_L on chr5 (R face) → [10,7]', () => {
    const s = locusSite(0, '5');
    should(s).match({ face: 'R', p5: 66890824, p3: 66890803 });
    should(realign.revcomp(plainSeq('5', 66890803, 66890824))).equal('GGACAGCTCCACGACCTATCAG');
    const seq = plainSeq('5', 66890800, 66890827);
    const r = realign.realignSiteOnWindow(s, seq, 66890800, 5);
    should(r).match({ mm: 2, mm_3p: 0, terminal_mm: false, gaps: 0, mm_pos: [10, 7], p5: 66890824, p3: 66890803, approx: false });
  });

  it('P2_R on chr5 (F face) → [15,10]', () => {
    const s = locusSite(1, '5');
    const r = realign.realignSiteOnWindow(s, plainSeq('5', 66890612, 66890638), 66890612, 5);
    should(r).match({ mm: 2, mm_pos: [15, 10], p5: 66890615, p3: 66890635 });
  });

  it('P3_L in the sorghum_353 footprint → [1]', () => {
    const w = DOC.sorghum_353;
    const site = { key: 'P3_L', primer: P.P3_L, len: 24, subject: '4', face: 'F', strand: 1, p5: 7499931, p3: 7499954, t5: 0, t3: 0, mm: 1 };
    const seq = w.seq.slice(7499928 - w.start, 7499957 - w.start + 1);
    const r = realign.realignSiteOnWindow(site, seq, 7499928, 5);
    should(r).match({ mm: 1, mm_3p: 1, terminal_mm: true, gaps: 0, mm_pos: [1], p5: 7499931, p3: 7499954 });
  });

  it('del18 → mm 1 with one gap near the 3\' end', () => {
    const s = locusSite(3, '4');
    should(s).match({ t3: 4, p5: 7423537, p3: 7423557 });
    const r = realign.realignSiteOnWindow(s, plainSeq('4', 7423534, 7423560), 7423534, 5);
    should(r).match({ mm: 1, gaps: 1, mm_3p: 1, terminal_mm: false, mm_pos: [5], p5: 7423537, p3: 7423558 });
  });

  it('requires the primer sequence on the site', () => {
    const s = Object.assign(locusSite(0, '4'), { primer: null });
    should(() => realign.realignSiteOnWindow(s, plainSeq('4', 7423534, 7423561), 7423534, 5)).throw(TypeError);
  });

  it('m20_7_14 → [14,7]; ins12 → one inserted base', () => {
    const m20 = locusSite(2, '4');
    should(realign.realignSiteOnWindow(m20, plainSeq('4', 7423534, 7423559), 7423534, 5)).match({ mm: 2, mm_pos: [14, 7], gaps: 0 });
    const [line] = readLines('ins12_genome_loci.tsv');
    const h = blast.parseLine(line, 'genome');
    h.primer = P.ins12_T;
    const s = sites.hitToSite(h);
    const r = realign.realignSiteOnWindow(s, plainSeq('4', 7423534, 7423562), 7423534, 5);
    should(r).match({ mm: 1, gaps: 1, mm_pos: [11], p5: 7423537, p3: 7423558 });
  });
});

describe('realign.realignGenomeSites', () => {
  const allSites = () => readLines('p2_genome_loci.tsv').map((l) => {
    const h = blast.parseLine(l, 'genome');
    h.primer = QUERY[h.q];
    return sites.hitToSite(h);
  });

  it('returns results parallel to the input and merges nearby fetches', async () => {
    const input = allSites();
    const fetchWindow = bicolorFetch();
    const res = await realign.realignGenomeSites(input, { fetchWindow, fastaPath: '/fa.gz', threePrimeWindow: 5 });
    should(res.length).equal(input.length);
    res.forEach((r, k) => {
      should(r.key).equal(input[k].key);
      should(r.approx).be.false();
    });
    // 12 sites in 3 loci (two primers facing each other ~210 bp apart) → 3 fetches, plus strand.
    should(fetchWindow.calls.length).equal(3);
    should(fetchWindow.calls.every((c) => c.strand === 1 && c.fastaPath === '/fa.gz')).be.true();
    const byKey = realign.indexAlignments(input, res);
    should(byKey.get(QUERY[0] + ':5:R:66890824')).match({ mm_pos: [10, 7] });
    should(byKey.get(QUERY[3] + ':4:F:7423537')).match({ mm: 1, gaps: 1 });
  });

  it('retries up to the aligned end when the padded window runs off the region', async () => {
    // Region "t" is 25 bp and ends 1 bp after the HSP; the 3' tail (t3 = 2) plus the
    // pad run past its end, so the first read fails like sequence.fetch does.
    const region = 'GGGG' + QUERY[0].slice(0, 20) + 'C';
    const calls = [];
    const regionFetch = async (fp, name, s, e) => {
      calls.push([s, e]);
      if (e > region.length) {
        const err = new Error('short read');
        err.code = 'REGION_OUT_OF_BOUNDS';
        throw err;
      }
      return region.slice(s - 1, e);
    };
    const site = sites.hitToSite({ q: 0, primer: QUERY[0], qlen: 22, qstart: 1, qend: 20, sseqid: 't', sstart: 5, send: 24, strand: 1, mismatch: 0 });
    should(site.p3).equal(26);
    const [r] = await realign.realignGenomeSites([site], { fetchWindow: regionFetch, fastaPath: 'x' });
    should(calls).eql([[2, 29], [2, 24]]);
    should(r.approx).be.false();
    should(r).not.have.property('fetch_error');
    should(r.mm).be.aboveOrEqual(2);
  });

  it('clips the window start at 1', async () => {
    const fetchWindow = async (fp, name, s, e) => {
      if (s < 1) throw new Error('regionStart cannot be less than 0');
      return (QUERY[1] + 'AAAAAAAA').slice(s - 1, e);
    };
    const site = sites.hitToSite({ q: 1, primer: QUERY[1], qlen: 21, qstart: 1, qend: 21, sseqid: 'c', sstart: 1, send: 21, strand: 1, mismatch: 0 });
    const [r] = await realign.realignGenomeSites([site], { fetchWindow, fastaPath: 'x' });
    should(r).match({ mm: 0, p5: 1, p3: 21, approx: false });
  });

  it('uses regionLength when provided', async () => {
    const calls = [];
    const fetchWindow = async (fp, name, s, e) => {
      calls.push([s, e]);
      return ('TTT' + QUERY[1]).slice(s - 1, e);
    };
    const site = sites.hitToSite({ q: 1, primer: QUERY[1], qlen: 21, qstart: 1, qend: 21, sseqid: 'c', sstart: 4, send: 24, strand: 1, mismatch: 0 });
    const [r] = await realign.realignGenomeSites([site], { fetchWindow, fastaPath: 'x', regionLength: async () => 24 });
    should(calls).eql([[1, 24]]);
    should(r).match({ mm: 0, p5: 4, p3: 24 });
  });

  it('falls back to approximate values when FASTA is unavailable', async () => {
    const site = sites.hitToSite({ q: 0, primer: QUERY[0], qlen: 22, qstart: 1, qend: 19, sseqid: '9', sstart: 100, send: 118, strand: 1, mismatch: 1 });
    const failing = async () => { const e = new Error('nope'); e.code = 'UNKNOWN_REGION'; throw e; };
    const [a] = await realign.realignGenomeSites([site], { fetchWindow: failing, fastaPath: 'x', threePrimeWindow: 5 });
    should(a).match({ approx: true, fetch_error: 'UNKNOWN_REGION', mm: 1 + 2, mm_pos: null, mm_3p: 2, terminal_mm: true, p5: 100 });
    const [b] = await realign.realignGenomeSites([site], {});
    should(b).match({ approx: true, fetch_error: 'NO_FASTA_FOR_REALIGN' });
    const [c] = await realign.realignGenomeSites([site], { fetchWindow: async () => undefined, fastaPath: 'x' });
    should(c).match({ approx: true, fetch_error: 'NO_SEQUENCE' });
  });

  it('honours an aborted signal', async () => {
    const ac = new AbortController();
    ac.abort();
    await assert.rejects(realign.realignGenomeSites(allSites(), { fetchWindow: bicolorFetch(), fastaPath: 'x', signal: ac.signal }), { code: 'ABORTED' });
  });
});

describe('realign.realignCdnaHit', () => {
  const cdna = readLines('p2_cdna.tsv').map((l) => blast.parseLine(l, 'cdna'));

  it('gives exact positions from qseq/sseq for full-length hits', () => {
    const l = cdna.find((h) => h.q === 0 && h.sseqid === 'SORBI_3005G183900.1');
    should(realign.realignCdnaHit(l, { threePrimeWindow: 5 })).match({ mm: 2, mm_3p: 0, terminal_mm: false, gaps: 0, mm_pos: [10, 7], p5: 1069, p3: 1090, approx: false });
    const r = cdna.find((h) => h.q === 1 && h.sseqid === 'SORBI_3005G183900.1');
    should(realign.realignCdnaHit(r)).match({ mm: 2, mm_pos: [15, 10], p5: 1278, p3: 1258, approx: false });
    const perfect = cdna.find((h) => h.q === 0 && h.sseqid === 'SORBI_3004G087700.1');
    should(realign.realignCdnaHit(sites.hitToSite(perfect))).match({ mm: 0, mm_pos: [], approx: false });
  });

  it('adds lbU tail edits and marks tails approximate', () => {
    // Real cDNA hit: q1-11 of P2_L on SORBI_3003G216300.3 (t3 = 11).
    const h = blast.parseLine('q0\tSORBI_3003G216300.3\t22\t1\t11\t2494\t2484\tminus\t0\tGGACAGCTCCA\tGGACAGCTCCA', 'cdna');
    should(realign.realignCdnaHit(h, { threePrimeWindow: 5 })).match({ mm: 6, mm_pos: [11, 9, 7, 5, 3, 1], mm_3p: 3, terminal_mm: true, approx: true, p5: 2494, p3: 2473 });
    // Synthetic 5' tail: t5 = 3 → primer positions 3 and 1 (distances 18, 20); internal
    // mismatch at primer position 4+8 = 12 → distance 9.
    const t5 = { qlen: 20, qstart: 4, qend: 20, qseq: 'AAAAAAAAAAAAAAAAA', sseq: 'AAAAAAAATAAAAAAAA', strand: 1, sstart: 50, send: 66 };
    should(realign.realignCdnaHit(t5)).match({ mm: 3, mm_pos: [20, 18, 9], mm_3p: 0, terminal_mm: false, approx: true, p5: 47, p3: 66 });
  });

  it('handles gap characters and rejects bad input', () => {
    const del = { qlen: 10, qstart: 1, qend: 10, qseq: 'ACGTA-CGTAC', sseq: 'ACGTAGCGTAC', strand: 1, sstart: 1, send: 11 };
    should(realign.realignCdnaHit(del)).match({ mm: 1, gaps: 1, mm_pos: [6] });
    const ins = { qlen: 11, qstart: 1, qend: 11, qseq: 'ACGTATCGTAC', sseq: 'ACGTA-CGTAC', strand: 1, sstart: 1, send: 10 };
    should(realign.realignCdnaHit(ins)).match({ mm: 1, gaps: 1, mm_pos: [6] });
    should(() => realign.realignCdnaHit({ qlen: 5, qstart: 1, qend: 5, qseq: 'ACGTA' })).throw(TypeError);
  });
});

describe('realign.realignFromHit', () => {
  it('uses BLAST mismatches plus tail lower bounds without positions', () => {
    const s = { key: 'k', len: 22, t5: 2, t3: 3, mm: 1, p5: 10, p3: 31 };
    should(realign.realignFromHit(s, 5)).eql({ key: 'k', mm: 1 + 1 + 2, mm_3p: 2, terminal_mm: true, gaps: 0, mm_pos: null, p5: 10, p3: 31, approx: true });
    should(realign.tailDistances(2, 3, 22)).eql([3, 1, 21]);
  });

  it('reports a known lower bound instead of the tail bound when given (bounded sites)', () => {
    const s = { key: 'k', len: 22, t5: 0, t3: 4, mm: 4, p5: 100, p3: 121 };
    should(realign.realignFromHit(s, 5).mm).equal(6);
    should(realign.realignFromHit(s, 5, { mm: 5 })).eql({ key: 'k', mm: 5, mm_3p: 2, terminal_mm: false, gaps: 0, mm_pos: null, p5: 100, p3: 121, approx: true });
    should(realign.realignFromHit({ key: 'k', len: 22, t5: 0, t3: 1, mm: 0, p5: 100, p3: 121 }, 5, { mm: 0 })).match({ mm: 0, mm_3p: 0, terminal_mm: false });
  });

  it('clamps coordinates to the sequence (chk-coords-below-one)', () => {
    should(realign.realignFromHit({ key: 'k', len: 22, t5: 3, t3: 0, mm: 0, p5: -2, p3: 19 }, 5)).match({ mm: 2, p5: 1, p3: 19 });
    should(realign.realignFromHit({ key: 'k', len: 22, t5: 0, t3: 0, mm: 0, p5: 90, p3: 111 }, 5, { regionLength: 100 })).match({ p5: 90, p3: 100 });
  });
});

describe('realigned coordinates stay inside the sequence (chk-coords-below-one)', () => {
  const primer = QUERY[1]; // P2_R, 21 nt

  it('a primer whose 5\' bases overhang the contig start keeps the overhang edits but starts at 1', async () => {
    const contig = primer.slice(3) + 'ACGTTGCAACGTTGCA';
    const site = sites.hitToSite({ q: 1, primer, qlen: 21, qstart: 4, qend: 21, sseqid: 'ctg', sstart: 1, send: 18, strand: 1, mismatch: 0 });
    should(site).match({ p5: -2, p3: 18 });
    should(sites.passesParseFilter(site)).be.true();
    const fetchWindow = async (fp, name, s, e) => contig.slice(s - 1, e);
    const [r] = await realign.realignGenomeSites([site], { fetchWindow, fastaPath: 'x', regionLength: async () => contig.length });
    should(r).match({ mm: 3, gaps: 3, mm_pos: [21, 20, 19], mm_3p: 0, p5: 1, p3: 18, approx: false, ops: 'III' + 'M'.repeat(18) });
    const [r2] = await realign.realignGenomeSites([site], { fetchWindow, fastaPath: 'x' });
    should(r2).match({ mm: 3, p5: 1 });
    should(realign.realignSiteOnWindow(site, contig.slice(0, 21), 1, 5)).match({ p5: 1, p3: 18 });
  });

  it('a reverse-facing primer overhanging the contig end ends at the contig length', async () => {
    const contig = 'ACGTTGCAACGTTGCA' + realign.revcomp(primer).slice(0, 18);
    const L = contig.length;
    const site = sites.hitToSite({ q: 1, primer, qlen: 21, qstart: 4, qend: 21, sseqid: 'ctg', sstart: L, send: L - 17, strand: -1, mismatch: 0 });
    should(site).match({ face: 'R', p5: L + 3, p3: L - 17 });
    const fetchWindow = async (fp, name, s, e) => {
      if (e > L) throw Object.assign(new Error('past the end'), { code: 'REGION_OUT_OF_BOUNDS' });
      return contig.slice(s - 1, e);
    };
    const [r] = await realign.realignGenomeSites([site], { fetchWindow, fastaPath: 'x', regionLength: async () => L });
    should(r).match({ mm: 3, mm_pos: [21, 20, 19], p5: L, p3: L - 17, approx: false });
  });
});
