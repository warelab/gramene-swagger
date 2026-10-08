'use strict';

const { describe, it } = require('node:test');
const fs = require('fs');
const path = require('path');
const should = require('should');

const blast = require('../../../api/helpers/primers/check/blast');
const sites = require('../../../api/helpers/primers/check/sites');
const realign = require('../../../api/helpers/primers/check/realign');

const FIX = path.join(__dirname, '..', 'fixtures', 'check_core');
// Real blastn-short output (make_indel_fixture.js) for P2_L / L26 on sorghum_bicolor 4:7423301-7424000 with a
// planted insertion (ins22, and its reverse complement ins22rc) or deletion (del26); see indel_split.json notes.
const IND = JSON.parse(fs.readFileSync(path.join(FIX, 'indel_split.json'), 'utf8'));
const IQ = IND.query_order.map((k) => IND.primers[k]);
const readLines = (name) => fs.readFileSync(path.join(FIX, name), 'utf8').split('\n').filter(Boolean);
const QUERY = ['GGACAGCTCCACAACATATCAG', 'GGACATTTGAAGCCCATGGCC', 'GGACAGATCCACATCATATC', 'GGACAGCTCCACAACATTCAG'];
const hit = (line, primers) => {
  const h = blast.parseLine(line, 'genome');
  if (primers) h.primer = primers[h.q];
  return h;
};

describe('sites.hitToSite', () => {
  it('converts a plus-strand hit to an F site', () => {
    const s = sites.hitToSite(hit('q0\t4\t22\t1\t22\t7423537\t7423558\tplus\t0', QUERY));
    should(s).match({ q: 0, primer: QUERY[0], len: 22, subject: '4', strand: 1, face: 'F', p5: 7423537, p3: 7423558, t5: 0, t3: 0, mm: 0, lbU: 0, lbI: 0 });
    should(s.key).equal(QUERY[0] + ':4:F:7423537');
  });

  it('converts a minus-strand hit to an R site', () => {
    const s = sites.hitToSite(hit('q0\t5\t22\t1\t22\t66890824\t66890803\tminus\t2', QUERY));
    should(s).match({ subject: '5', strand: -1, face: 'R', p5: 66890824, p3: 66890803, mm: 2, lbU: 2, lbI: 2 });
  });

  it('extrapolates tails on both strands', () => {
    // Real line from the ws5 run: q2-17 of 22 on the minus strand.
    const m = sites.hitToSite(hit('q0\t4\t22\t2\t17\t19468232\t19468217\tminus\t0', QUERY));
    should(m).match({ face: 'R', t5: 1, t3: 5, p5: 19468233, p3: 19468212, lbU: 4, lbI: 2 });
    should(Math.abs(m.p5 - m.p3) + 1).equal(22);
    const p = sites.hitToSite({ q: 1, qlen: 21, qstart: 3, qend: 19, sseqid: 'chr', sstart: 1000, send: 1016, strand: 1, mismatch: 1 });
    should(p).match({ face: 'F', t5: 2, t3: 2, p5: 998, p3: 1018, lbU: 3, lbI: 3 });
    should(p.key).equal('q1:chr:F:998');
  });

  it('keeps cDNA alignment strings', () => {
    const line = readLines('p2_cdna.tsv').find((l) => l.startsWith('q1\tSORBI_3005G183900.1\t'));
    const h = blast.parseLine(line, 'cdna');
    const s = sites.hitToSite(h);
    should(s).match({ face: 'R', p5: 1278, p3: 1258, qseq: QUERY[1], sseq: 'GGACATCTGAAACCCATGGCC', qstart: 1, qend: 21 });
  });
});

describe('sites filters', () => {
  it('computes lbU and lbI', () => {
    should(sites.lowerBoundU(1, 3, 4)).equal(1 + 2 + 2);
    should(sites.lowerBoundI(1, 3, 4)).equal(3);
    should(sites.lowerBoundI(2, 0, 0)).equal(2);
  });

  it('keeps del18 (q1-17, t3=4) through the parse filter', () => {
    const lines = readLines('p2_genome_loci.tsv').filter((l) => l.startsWith('q3\t4\t'));
    should(lines.length).equal(2);
    for (const l of lines) {
      const s = sites.hitToSite(hit(l, QUERY));
      should(s).match({ t5: 0, t3: 4, mm: 0, lbU: 2, lbI: 1 });
      should(sites.passesParseFilter(s)).be.true();
    }
  });

  it('drops ins12 (q1-12, t3=11) under the plan formula', () => {
    // Real ws5 hit of GGACAGCTCCACTAACATATCAG (23 nt). lbU = ceil(11/2) = 6 > 5 and
    // t5+t3 = 11 > 8, so the §B.4/plan parse filter drops it. §A.12 expects it to
    // survive; that expectation contradicts the formula (reported as a spec issue).
    const [line] = readLines('ins12_genome_loci.tsv');
    const s = sites.hitToSite(hit(line));
    should(s).match({ t5: 0, t3: 11, mm: 0, lbU: 6, lbI: 1 });
    should(sites.passesParseFilter(s)).be.false();
  });

  it('applies the parse filter boundaries', () => {
    const site = (mm, t5, t3) => ({ mm, t5, t3, lbU: sites.lowerBoundU(mm, t5, t3), lbI: sites.lowerBoundI(mm, t5, t3) });
    should(sites.passesParseFilter(site(5, 0, 0))).be.true();
    should(sites.passesParseFilter(site(6, 0, 0))).be.false();
    should(sites.passesParseFilter(site(0, 0, 10))).be.true(); // lbU 5
    should(sites.passesParseFilter(site(0, 0, 11))).be.false(); // lbU 6, tails 11
    should(sites.passesParseFilter(site(3, 4, 4))).be.true(); // lbU 7, lbI 5, tails 8
    should(sites.passesParseFilter(site(3, 4, 5))).be.false(); // tails 9
    should(sites.passesParseFilter(site(4, 4, 4))).be.false(); // lbI 6
  });

  it('applies the job filter with the tail-dependent bound', () => {
    const s1 = { mm: 4, t5: 0, t3: 4, lbU: 6, lbI: 5 };
    should(sites.jobLowerBound(s1)).equal(5);
    should(sites.passesJobFilter(s1, 6)).be.true();
    should(sites.passesJobFilter(s1, 5)).be.false();
    const s2 = { mm: 1, t5: 0, t3: 9, lbU: 6, lbI: 2 };
    should(sites.jobLowerBound(s2)).equal(6); // tails > 8: lbU only
    should(sites.passesJobFilter(s2, 6)).be.false();
    should(sites.passesJobFilter({ mm: 2, t5: 0, t3: 0, lbU: 2, lbI: 2 }, 3)).be.true();
  });

  it('identifies near-perfect sites and hits', () => {
    should(sites.isNearPerfect({ t5: 0, t3: 0, mm: 1 })).be.true();
    should(sites.isNearPerfect({ t5: 0, t3: 0, mm: 2 })).be.false();
    should(sites.isNearPerfect({ t5: 1, t3: 0, mm: 0 })).be.false();
    should(sites.isNearPerfectHit({ qstart: 1, qend: 22, qlen: 22, mismatch: 0 })).be.true();
    should(sites.isNearPerfectHit({ qstart: 1, qend: 21, qlen: 22, mismatch: 0 })).be.false();
  });

  it('computes the site cap', () => {
    should(sites.siteCap(708735318)).equal(50000);
    should(sites.siteCap(3e9)).equal(180000);
    should(sites.siteCap(1000000001)).equal(60001);
    should(sites.siteCap(undefined)).equal(50000);
  });
});

describe('sites.SiteStore', () => {
  it('stores parse-filtered hits in Int32 columns and counts near-perfect hits', () => {
    const store = new sites.SiteStore({ primers: QUERY.map((s) => s.toLowerCase()), totalBases: 708735318, initialCapacity: 2 });
    should(store.cap).equal(50000);
    const lines = readLines('p2_genome_sample.tsv').concat(readLines('p2_genome_loci.tsv'));
    const expected = QUERY.map(() => ({ hits: 0, sites: 0, near: 0 }));
    for (const l of lines) {
      const h = hit(l, QUERY);
      const s = sites.hitToSite(h);
      const e = expected[h.q];
      e.hits++;
      if (sites.isNearPerfect(s)) e.near++;
      const kept = sites.passesParseFilter(s);
      if (kept) e.sites++;
      should(store.addHit(h)).equal(kept);
    }
    QUERY.forEach((seq, q) => {
      should(store.stats(q)).eql({ seq, len: seq.length, hits: expected[q].hits, sites: expected[q].sites, near_perfect_sites: expected[q].near, truncated: false });
    });
    should(store.stats(QUERY[0]).near_perfect_sites).equal(2);
    should(store.stats(QUERY[1]).near_perfect_sites).equal(2);
    should(store.totalSites()).equal(expected.reduce((a, e) => a + e.sites, 0));
    should(store._cols[0].p5).be.instanceOf(Int32Array);
  });

  it('materializes sites identical to hitToSite', () => {
    const store = new sites.SiteStore({ primers: QUERY, totalBases: 1 });
    const lines = readLines('p2_genome_loci.tsv');
    const direct = [];
    for (const l of lines) {
      const h = hit(l, QUERY);
      if (store.addHit(h)) direct.push(sites.hitToSite(h));
    }
    const fromStore = [];
    QUERY.forEach((seq, q) => { for (let i = 0; i < store.size(q); i++) fromStore.push(store.site(q, i)); });
    const byKey = (a, b) => (a.key < b.key ? -1 : 1);
    should(fromStore.sort(byKey)).eql(direct.sort(byKey));
    should(store.subjects.slice().sort()).eql(['4', '5']);
  });

  it('marks a primer truncated and repetitive past the cap', () => {
    const store = new sites.SiteStore({ primers: ['ACGTACGTACGTACGTAC'], siteCap: 3 });
    for (let i = 0; i < 5; i++) {
      store.addHit({ q: 0, sseqid: 'c' + (i % 2), qlen: 18, qstart: 1, qend: 18, sstart: 100 * (i + 1), send: 100 * (i + 1) + 17, strand: 1, mismatch: 0 });
    }
    should(store.stats(0)).match({ hits: 5, sites: 3, near_perfect_sites: 5, truncated: true });
    should(store.primerSummary(0, { repeat_site_threshold: 100 })).eql({ len: 18, near_perfect_sites: 5, repetitive: true, truncated: true });
  });

  it('derives repetitive from near-perfect sites and the threshold', () => {
    const store = new sites.SiteStore({ primers: ['ACGTACGTACGTACGTAC'] });
    for (let i = 0; i < 6; i++) {
      store.addHit({ q: 0, sseqid: '1', qlen: 18, qstart: 1, qend: 18, sstart: 1000 * (i + 1), send: 1000 * (i + 1) + 17, strand: 1, mismatch: i % 2 });
    }
    store.addHit({ q: 0, sseqid: '1', qlen: 18, qstart: 1, qend: 18, sstart: 99999, send: 99982, strand: -1, mismatch: 2 });
    should(store.primerSummary('acgtacgtacgtacgtac', { repeat_site_threshold: 5 })).eql({ len: 18, near_perfect_sites: 6, repetitive: true, truncated: false });
    should(store.primerSummary(0, { repeat_site_threshold: 6 }).repetitive).be.false();
    should(store.primerSummary(0).repetitive).be.true(); // default threshold 5
  });

  it('builds job-filtered columnar views', () => {
    const store = new sites.SiteStore({ primers: ['ACGTACGTACGTACGTACGT'] });
    const add = (mm, t3, sstart) => store.addHit({ q: 0, sseqid: 'x', qlen: 20, qstart: 1, qend: 20 - t3, sstart, send: sstart + 19 - t3, strand: 1, mismatch: mm });
    add(0, 0, 10); // lb 0
    add(4, 4, 20); // lbU 6, lbI 5 → lb 5
    add(5, 0, 30); // lb 5
    add(3, 0, 40); // lb 3
    const all = store.view(0);
    should(all.n).equal(4);
    should(Array.from(all.lb)).eql([0, 5, 5, 3]);
    const v = store.view(0, { ignoreMismatches: 5 });
    should(v.n).equal(2);
    should(Array.from(v.idx)).eql([0, 3]);
    should(v.p5[v.idx[1]]).equal(40);
    should(v.subjects[v.subj[0]]).equal('x');
    should(store.sites(0, { ignoreMismatches: 5 }).map((s) => s.p5)).eql([10, 40]);
  });

  it('keeps cDNA alignments when asked', () => {
    const store = new sites.SiteStore({ primers: QUERY.slice(0, 2), keepAlignments: true });
    for (const l of readLines('p2_cdna.tsv')) {
      const h = blast.parseLine(l, 'cdna');
      h.primer = QUERY[h.q];
      store.addHit(h);
    }
    should(store.size(0)).equal(6);
    const s = store.sites(1).find((x) => x.subject === 'SORBI_3005G183900.1');
    should(s).match({ face: 'R', p5: 1278, qseq: QUERY[1], sseq: 'GGACATCTGAAACCCATGGCC', qstart: 1, qend: 21 });
  });

  it('exports the gapped-site constants', () => {
    should(sites.MIN_HALF_ALIGNED).equal(11);
    should(sites.MERGE_P5_TOLERANCE).equal(3);
  });

  it('rejects bad input', () => {
    should(() => new sites.SiteStore({ primers: [] })).throw(TypeError);
    const store = new sites.SiteStore({ primers: ['ACGTACGTACGTACGT'] });
    should(() => store.addHit({ q: 3, qlen: 16, qstart: 1, qend: 16, sstart: 1, send: 16, strand: 1, mismatch: 0, sseqid: 'a' })).throw(RangeError);
    should(() => store.stats('TTTT')).throw(RangeError);
    should(() => store.site(0, 0)).throw(RangeError);
    should(store.indexOf('acgtacgtacgtacgt')).equal(0);
    should(store.indexOf(7)).equal(-1);
  });
});

describe('sites.SiteStore gapped sites from split HSPs (chk-indel-parse-filter)', () => {
  const hits = (target) => (target === 'cdna' ? IND.cdna_hits : IND.genome_hits).map((l) => blast.parseLine(l, target || 'genome'));
  const pick = (list, q, subject) => list.filter((h) => h.q === q && h.sseqid === subject);

  it('holds both halves of a mid-primer insertion (each dropped by the parse filter) and joins them in finalize()', () => {
    const halves = pick(hits(), 0, 'ins22');
    should(halves.map((h) => [h.qstart, h.qend, h.sstart, h.send, h.mismatch])).eql([[1, 11, 237, 247, 0], [12, 22, 249, 259, 0]]);
    const store = new sites.SiteStore({ primers: IQ });
    for (const h of halves) {
      const s = sites.hitToSite(h);
      should(s.lbU).equal(6);
      should(sites.passesParseFilter(s)).be.false();
      should(store.addHit(h)).be.false();
    }
    should(store.size(0)).equal(0); // what the pipeline saw before the fix
    should(store.finalize()).equal(1);
    should(store.finalize()).equal(0);
    should(store.stats(0)).match({ hits: 2, sites: 1, near_perfect_sites: 0, truncated: false });
    const [site] = store.sites(0, { ignoreMismatches: 6 });
    should(site).match({ primer: IQ[0], subject: 'ins22', face: 'F', strand: 1, p5: 237, p3: 259, t5: 0, t3: 0, mm: 1 });
    should(sites.jobLowerBound(site)).equal(1);
    should(store.gappedSites).equal(1);
    // the joined site re-aligns to the one-edit product the DP finds on the real window
    const r = realign.realignSiteOnWindow(site, IND.subjects.ins22.slice(233, 263), 234, 5);
    should(r).match({ mm: 1, gaps: 1, mm_3p: 0, p5: 237, p3: 259, approx: false });
  });

  it('joins minus-strand halves and a deletion whose halves overlap by one base', () => {
    const store = new sites.SiteStore({ primers: IQ });
    for (const h of hits()) store.addHit(h);
    should(store.finalize()).equal(5); // P2_L on ins22 and ins22rc; L26 on ins22, ins22rc and del26 (P2_L on del26 is kept as-is)
    const byKey = (q) => store.sites(q).map((s) => [s.subject, s.face, s.p5, s.p3, s.t5, s.t3, s.mm]).sort();
    should(byKey(0)).eql([
      ['del26', 'F', 237, 249 + 9, 0, 9, 0], // q1-13, t3 9: lbU 5 passes the parse filter unchanged
      ['ins22', 'F', 237, 259, 0, 0, 1],
      ['ins22rc', 'R', 465, 443, 0, 0, 1]
    ]);
    should(byKey(2).filter((s) => s[0] === 'del26')).eql([['del26', 'F', 237, 261, 0, 0, 1]]);
    should(IND.subjects.ins22rc.length).equal(701);
  });

  it('view() joins pending halves; halves that are far apart or contained are not joined', () => {
    const P22 = IQ[0];
    const store = new sites.SiteStore({ primers: [P22] });
    const h = (qstart, qend, sstart, strand) => ({ q: 0, sseqid: 'x', qlen: 22, qstart, qend, sstart, send: strand === 1 ? sstart + (qend - qstart) : sstart - (qend - qstart), strand, mismatch: 0 });
    store.addHit(h(1, 11, 1000, 1));
    store.addHit(h(12, 22, 1020, 1)); // extrapolated p5 1009: 9 away
    store.addHit(h(1, 11, 5000, 1));
    store.addHit(h(2, 12, 5001, 1)); // same diagonal, overlaps q2-11: joined it would still leave a 10 nt 3' tail
    store.addHit(h(12, 22, 2011, 1)); // p5 2000 ...
    store.addHit(h(1, 11, 2000, -1)); // ... but other strand
    should(store.size(0)).equal(0);
    should(store.view(0).n).equal(0);
    should(store.gappedSites).equal(0);
    const joinable = new sites.SiteStore({ primers: [P22] });
    joinable.addHit(h(1, 11, 3000, 1));
    joinable.addHit(h(12, 22, 3012, 1));
    should(joinable.view(0).n).equal(1);
    should(joinable.site(0, 0)).match({ p5: 3000, p3: 3022, mm: 1 });
  });

  it('cDNA stores join qseq/sseq so realignCdnaHit counts the indel', () => {
    const store = new sites.SiteStore({ primers: IQ, keepAlignments: true });
    for (const h of hits('cdna')) store.addHit(h);
    store.finalize();
    const ins = store.sites(0).find((s) => s.subject === 'ins22');
    should(ins).match({ qseq: IQ[0].slice(0, 11) + '-' + IQ[0].slice(11), sseq: IQ[0].slice(0, 11) + 'N' + IQ[0].slice(11), qstart: 1, qend: 22 });
    should(realign.realignCdnaHit(ins)).match({ mm: 1, gaps: 1, mm_pos: [12], approx: false, p5: 237, p3: 259 });
    const rc = store.sites(0).find((s) => s.subject === 'ins22rc');
    should(realign.realignCdnaHit(rc)).match({ mm: 1, gaps: 1, mm_pos: [12], p5: 465, p3: 443 });
    const del = store.sites(2).find((s) => s.subject === 'del26');
    should(del.qseq.replace(/-/g, '')).equal(IQ[2]);
    should(realign.realignCdnaHit(del)).match({ mm: 1, gaps: 1, approx: false, p5: 237, p3: 261 });
  });
});
