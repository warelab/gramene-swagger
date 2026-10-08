'use strict';

const { describe, it } = require('node:test');
const fs = require('fs');
const path = require('path');
const should = require('should');

const blast = require('../../../api/helpers/primers/check/blast');
const sites = require('../../../api/helpers/primers/check/sites');
const amplicons = require('../../../api/helpers/primers/check/amplicons');
const realign = require('../../../api/helpers/primers/check/realign');
const classify = require('../../../api/helpers/primers/check/classify');

const FIX = path.join(__dirname, '..', 'fixtures', 'check_core');
const DOC = JSON.parse(fs.readFileSync(path.join(FIX, 'sorghum_windows.json'), 'utf8'));
const readLines = (name) => fs.readFileSync(path.join(FIX, name), 'utf8').split('\n').filter(Boolean);
const P = DOC.primers;

const L = 'ACGTACGTACGTACGTACGT';
const R = 'TTTTGGGGCCCCAAAATTTT';

function mk(primer, subject, face, p5, extra) {
  const len = primer.length;
  return Object.assign({
    primer, len, subject, face, strand: face === 'F' ? 1 : -1, p5,
    p3: face === 'F' ? p5 + len - 1 : p5 - len + 1, t5: 0, t3: 0, mm: 0, lbU: 0, lbI: 0
  }, extra);
}

function toHit(s, q) {
  return {
    q, sseqid: s.subject, qlen: s.len, qstart: 1, qend: s.len, strand: s.strand, mismatch: s.mm,
    sstart: s.p5, send: s.face === 'F' ? s.p5 + s.len - 1 : s.p5 - s.len + 1
  };
}

const label = (c) => c.subject + ':' + c.orientation + ':' + c.start + '-' + c.end;

describe('amplicons.candidates', () => {
  const layout = () => ({
    [L]: [
      mk(L, 'chr1', 'F', 1000),
      mk(L, 'chr1', 'R', 5300),
      mk(L, 'chr2', 'F', 100),
      mk(L, 'chr2', 'R', 400),
      mk(L, 'chr4', 'F', 1),
      mk(L, 'chr6', 'F', 1),
      mk(L, 'chr7', 'F', 1),
      mk(L, 'chr8', 'F', 100),
      mk(L, 'chr9', 'R', 150)
    ],
    [R]: [
      mk(R, 'chr1', 'R', 1200),
      mk(R, 'chr1', 'F', 5000),
      mk(R, 'chr3', 'F', 10),
      mk(R, 'chr3', 'R', 110),
      mk(R, 'chr5', 'R', 100),
      mk(R, 'chr6', 'R', 4001),
      mk(R, 'chr7', 'R', 4000),
      mk(R, 'chr8', 'R', 138),
      mk(R, 'chr8', 'R', 139),
      mk(R, 'chr9', 'F', 200)
    ]
  });
  const expected = [
    'chr1:LR:1000-1200', // L F → R R
    'chr1:RL:5000-5300', // R F → L R
    'chr2:LL:100-400', // single-primer products
    'chr3:RR:10-110',
    'chr7:LR:1-4000', // exactly maxSize
    'chr8:LR:100-139' // size 40 = fwd + rev primer length; size 39 would overlap the footprints
  ];

  it('pairs F then R sites on one subject into LR/RL/LL/RR within size bounds', () => {
    const res = amplicons.candidates({ id: 'x', left: L, right: R }, layout(), { maxSize: 4000 });
    should(res.truncated).be.false();
    should(res.candidates.map(label).sort()).eql(expected.slice().sort());
    const lr = res.candidates.find((c) => c.orientation === 'LR' && c.subject === 'chr1');
    should(lr).match({ size: 201, start: 1000, end: 1200 });
    should(lr.fwd.primer).equal(L);
    should(lr.rev.primer).equal(R);
    should(lr.fwd.key).equal(L + ':chr1:F:1000');
    should(res.sites_considered).equal(19);
  });

  it('gives the same candidates from a Map and from a SiteStore', () => {
    const lay = layout();
    const fromMap = amplicons.candidates({ left: L.toLowerCase(), right: R }, new Map(Object.entries(lay)), { maxSize: 4000 });
    const store = new sites.SiteStore({ primers: [L, R] });
    lay[L].forEach((s) => store.addHit(toHit(s, 0)));
    lay[R].forEach((s) => store.addHit(toHit(s, 1)));
    const fromStore = amplicons.candidates({ left: L, right: R }, store, { maxSize: 4000 });
    should(fromMap.candidates.map(label)).eql(fromStore.candidates.map(label));
    should(fromStore.candidates.map(label).sort()).eql(expected.slice().sort());
    const c = fromStore.candidates.find((x) => x.orientation === 'RL');
    should(c.fwd).match({ primer: R, subject: 'chr1', face: 'F', p5: 5000, key: R + ':chr1:F:5000' });
    // The same site object is shared between candidates of one call.
    const again = amplicons.candidates({ left: L, right: R }, store, { maxSize: 4000 });
    should(again.candidates.map((x) => x.fwd.key)).eql(fromStore.candidates.map((x) => x.fwd.key));
  });

  it('orders candidates by subject name then position', () => {
    const res = amplicons.candidates({ left: L, right: R }, layout(), { maxSize: 4000 });
    should(res.candidates.map(label)).eql(expected);
  });

  it('respects the candidate cap and flags truncation', () => {
    const lay = { [L]: [], [R]: [] };
    for (let i = 0; i < 10; i++) {
      lay[L].push(mk(L, 'c', 'F', 1000 + i * 10));
      lay[R].push(mk(R, 'c', 'R', 3000 + i * 10));
    }
    should(amplicons.candidates({ left: L, right: R }, lay, { maxSize: 4000 })).match({ truncated: false });
    should(amplicons.candidates({ left: L, right: R }, lay, { maxSize: 4000 }).candidates.length).equal(100);
    const capped = amplicons.candidates({ left: L, right: R }, lay, { maxSize: 4000, maxCandidates: 30 });
    should(capped.candidates.length).equal(30);
    should(capped.truncated).be.true();
    const exact = amplicons.candidates({ left: L, right: R }, lay, { maxSize: 4000, maxCandidates: 100 });
    should(exact.truncated).be.false();
    should(amplicons.DEFAULT_MAX_CANDIDATES).equal(5000);
  });

  it('labels the F×R products of an identical-primer pair LR (chk-identical-lr)', () => {
    const lay = { [L]: [mk(L, 'c', 'F', 100), mk(L, 'c', 'R', 300)] };
    const res = amplicons.candidates({ left: L, right: L }, lay, { maxSize: 4000 });
    should(res.candidates.map(label)).eql(['c:LR:100-300']);
    const [c] = res.candidates;
    const alns = new Map([[c.fwd.key, { p5: 100, p3: 119, mm: 0, mm_pos: [], mm_3p: 0 }], [c.rev.key, { p5: 300, p3: 281, mm: 0, mm_pos: [], mm_3p: 0 }]]);
    const [amp] = amplicons.finalize(res.candidates, alns, {}).amplicons;
    should(amp).match({ orientation: 'LR', strand: 1, start: 100, end: 300, likelihood: 'likely' });
    // so it can be the expected on-target product (it never could as LL)
    should(require('../../../api/helpers/primers/check/specificity').selectGenome({ amplicons: [amp], mode: 'region', expected: { region: 'c', start: 100, end: 300 } }))
      .match({ verdict: 'specific', onTarget: amp });
  });

  it('requires size >= fwd primer length + rev primer length (non-overlapping footprints)', () => {
    const R25 = 'CCCCAAAATTTTGGGGAAAACCCCA';
    const lay = { [L]: [mk(L, 'c', 'F', 100)], [R25]: [mk(R25, 'c', 'R', 143), mk(R25, 'c', 'R', 144)] };
    should(amplicons.candidates({ left: L, right: R25 }, lay, { maxSize: 4000 }).candidates.map(label)).eql(['c:LR:100-144']);
    const rl = { [R25]: [mk(R25, 'd', 'F', 100)], [L]: [mk(L, 'd', 'R', 143), mk(L, 'd', 'R', 144)] };
    should(amplicons.candidates({ left: L, right: R25 }, rl, { maxSize: 4000 }).candidates.map(label)).eql(['d:RL:100-144']);
  });

  it('applies the job filter and drops duplicate sites', () => {
    const lay = {
      [L]: [mk(L, 'c', 'F', 100), mk(L, 'c', 'F', 100, { t3: 2, lbU: 1, lbI: 1 }), mk(L, 'd', 'F', 100, { mm: 5, lbU: 5, lbI: 5 })],
      [R]: [mk(R, 'c', 'R', 300), mk(R, 'd', 'R', 300)]
    };
    const all = amplicons.candidates({ left: L, right: R }, lay, { maxSize: 4000 });
    should(all.candidates.map(label)).eql(['c:LR:100-300', 'd:LR:100-300']);
    should(all.candidates[0].fwd.t3).equal(0); // best lower bound kept
    const filtered = amplicons.candidates({ left: L, right: R }, lay, { maxSize: 4000, ignoreMismatches: 5 });
    should(filtered.candidates.map(label)).eql(['c:LR:100-300']);
  });

  it('handles missing primers and validates maxSize', () => {
    should(amplicons.candidates({ left: L, right: R }, {}, { maxSize: 100 }).candidates).eql([]);
    should(() => amplicons.candidates({ left: L, right: R }, {}, {})).throw(TypeError);
  });

  it('computes the pan-genome maxSize', () => {
    should(amplicons.maxSizeFor('reference', 4000, 210)).equal(4000);
    should(amplicons.maxSizeFor('pangenome', 4000, 580)).equal(4000);
    should(amplicons.maxSizeFor('pangenome', 4000, 3000)).equal(5000);
    should(amplicons.maxSizeFor('pangenome', 1000, 1001)).equal(2002);
    should(amplicons.maxSizeFor('pangenome', 4000, null)).equal(4000);
  });
});

describe('amplicons.uniqueSites', () => {
  it('collects unique sites across pairs and caps them', () => {
    const a = mk(L, 'c', 'F', 1);
    const b = mk(R, 'c', 'R', 200);
    const c = mk(R, 'c', 'F', 300);
    const d = mk(L, 'c', 'R', 500);
    [a, b, c, d].forEach((s) => { s.key = sites.siteKey(s.primer, s.subject, s.face, s.p5); });
    const lists = [[{ fwd: a, rev: b }, { fwd: a, rev: d }], [{ fwd: c, rev: d }]];
    should(amplicons.uniqueSites(lists)).match({ truncated: false });
    should(amplicons.uniqueSites(lists).sites.map((s) => s.key)).eql([a.key, b.key, d.key, c.key]);
    should(amplicons.uniqueSites(lists[0]).sites.length).equal(3);
    const capped = amplicons.uniqueSites(lists, { maxSites: 2 });
    should(capped.sites.length).equal(2);
    should(capped.truncated).be.true();
    should(amplicons.DEFAULT_MAX_REALIGN_SITES).equal(20000);
  });
});

describe('amplicons.realignPlan', () => {
  const s = (primer, face, p5, extra) => {
    const x = mk(primer, 'c', face, p5, extra);
    x.key = sites.siteKey(primer, 'c', face, p5);
    return x;
  };

  it('re-aligns sites whose gap-aware bound is within max_amplifying_mismatches; the realign cap applies to re-aligned sites only', () => {
    const good = s(L, 'F', 100); // bound 0
    const three = s(R, 'R', 400, { mm: 3 }); // bound 3
    const four = s(R, 'R', 500, { mm: 4 }); // bound 4 > 3
    const tails = s(L, 'F', 200, { mm: 2, t3: 4 }); // lbU 4, lbI 3 → bound 3
    // one side of an indel: q1-13 of 22 (t3 9) has jobLowerBound (lbU) 5 but may re-align with 1 edit: lbI 1
    const half = s(L, 'F', 300, { t3: 9 });
    should([sites.jobLowerBound(half), amplicons.realignBound(half)]).eql([5, 1]);
    const lists = [[{ fwd: good, rev: three }, { fwd: good, rev: four }], [{ fwd: tails, rev: four }, { fwd: half, rev: three }]];
    const plan = amplicons.realignPlan(lists, { maxAmplifyingMismatches: 3 });
    should(plan.realign.map((x) => x.key)).eql([good.key, three.key, tails.key, half.key]);
    should(plan.bounded.map((x) => x.key)).eql([four.key]);
    should(plan.truncated).be.false();
    should(amplicons.realignPlan(lists, {}).realign.length).equal(5); // no cap: every site, like uniqueSites
    should(amplicons.realignPlan(lists[0]).bounded).eql([]);
    const capped = amplicons.realignPlan(lists, { maxAmplifyingMismatches: 3, maxSites: 1 });
    should(capped.realign.map((x) => x.key)).eql([good.key]);
    should(capped.bounded.map((x) => x.key)).eql([four.key]);
    should(capped.truncated).be.true();
  });

  it('bounded sites (realignFromHit with mm = realignBound) give approximate unlikely products, never likely or ignored', () => {
    const f = s(L, 'F', 100);
    const fa = { p5: 100, p3: 119, mm: 0, mm_pos: [], mm_3p: 0, terminal_mm: false, approx: false };
    const r4 = s(R, 'R', 400, { mm: 4 });
    const out = amplicons.finalize([{ subject: 'c', orientation: 'LR', fwd: f, rev: r4 }],
      new Map([[f.key, fa], [r4.key, realign.realignFromHit(r4, 5, { mm: amplicons.realignBound(r4) })]]), {});
    should(out.amplicons).eql([]);
    should(out.unlikely).match([{ likelihood: 'unlikely', right_mm: 4, right_mm_pos: null, approx: true }]);
    // lbU 6 (the plain fallback would be ignored) but gap-aware bound 5: stays unlikely
    const r5 = s(R, 'R', 700, { mm: 4, t3: 4 });
    should(realign.realignFromHit(r5, 5).mm).equal(6);
    const out5 = amplicons.finalize([{ subject: 'c', orientation: 'LR', fwd: f, rev: r5 }],
      new Map([[f.key, fa], [r5.key, realign.realignFromHit(r5, 5, { mm: amplicons.realignBound(r5) })]]), {});
    should(out5.counts).match({ unlikely: 1, ignored: 0 });
    should(out5.unlikely[0]).match({ right_mm: 5, right_3p_mm: 2, approx: true });
  });
});

describe('classify rules (§B.7)', () => {
  const aln = (mmPos, extra) => Object.assign({
    mm: mmPos.length, mm_pos: mmPos, mm_3p: mmPos.filter((d) => d <= 5).length, terminal_mm: mmPos.includes(1), gaps: 0, approx: false
  }, extra);
  const perfect = aln([]);

  it('has the plan defaults', () => {
    should(classify.DEFAULT_PARAMS).eql({
      max_product_size: 4000, ignore_mismatches: 6, min_total_mismatches: 2, min_3p_mismatches: 2,
      three_prime_window: 5, include_unlikely: false, repeat_site_threshold: 5, max_amplifying_mismatches: 3
    });
  });

  it('classifies the rule matrix with defaults', () => {
    const cases = [
      [perfect, perfect, 'likely'],
      [aln([20, 18, 16, 14, 6]), perfect, 'unlikely'], // 5 mm > max_amplifying_mismatches 3, none in the 3' window
      [aln([20, 18, 16, 14]), perfect, 'unlikely'], // 4 mm
      [perfect, aln([19, 15, 11, 7]), 'unlikely'],
      [aln([20, 16, 12]), aln([14, 9, 6]), 'likely'], // 3 mm per primer, outside the window
      [aln([20, 16, 1]), perfect, 'likely_weak'], // 3 mm with a terminal mismatch
      [aln([20, 18, 16, 1]), perfect, 'unlikely'], // the cap beats the terminal flag
      [aln([10, 7]), aln([15, 10]), 'likely'], // P2 chr5 off-target
      [aln([3]), perfect, 'likely'], // 1 mm < min_total_mismatches
      [aln([9, 3]), perfect, 'likely'], // 2 mm, only 1 in the window
      [aln([1]), perfect, 'likely_weak'], // terminal mismatch
      [perfect, aln([12, 1]), 'likely_weak'],
      [aln([4, 2]), perfect, 'unlikely'], // blocked
      [perfect, aln([20, 5, 4]), 'unlikely'],
      [aln([2, 1]), perfect, 'unlikely'], // blocked beats terminal
      [aln([22, 20, 18, 16, 14, 12]), perfect, 'ignored'], // 6 mm
      [perfect, aln([22, 20, 18, 16, 14, 12, 10]), 'ignored'],
      [aln([22, 20, 18, 16, 3, 2]), aln([5, 4]), 'ignored'] // ignored beats unlikely
    ];
    for (const [f, r, want] of cases) should(classify.classifyAmplicon(f, r)).equal(want);
  });

  it('honours the parameters', () => {
    // max_amplifying_mismatches: <= cap amplifies, cap < mm < ignore_mismatches is unlikely
    should(classify.classifyAmplicon(aln([20, 18, 16, 14, 6]), perfect, { max_amplifying_mismatches: 5 })).equal('likely');
    should(classify.classifyAmplicon(aln([20, 18, 16, 14]), perfect, { max_amplifying_mismatches: 4 })).equal('likely');
    should(classify.classifyAmplicon(aln([12]), perfect, { max_amplifying_mismatches: 0 })).equal('unlikely');
    should(classify.classifyAmplicon(perfect, perfect, { max_amplifying_mismatches: 0 })).equal('likely');
    should(classify.isOverAmplifyingCap(aln([9, 8, 7, 6]))).be.true();
    should(classify.isOverAmplifyingCap(aln([9, 8, 7]))).be.false();
    should(classify.classifyAmplicon(aln([20, 19, 18]), perfect, { ignore_mismatches: 3 })).equal('ignored');
    should(classify.classifyAmplicon(aln([2]), perfect, { min_total_mismatches: 0, min_3p_mismatches: 1 })).equal('unlikely');
    should(classify.classifyAmplicon(aln([2]), perfect, { min_total_mismatches: 2, min_3p_mismatches: 1 })).equal('likely');
    // mm_3p is recomputed from mm_pos for the requested window.
    should(classify.classifyAmplicon(aln([5, 4]), perfect, { three_prime_window: 5 })).equal('unlikely');
    should(classify.classifyAmplicon(aln([5, 4]), perfect, { three_prime_window: 3 })).equal('likely');
    should(classify.classifyAmplicon(aln([8, 7]), perfect, { three_prime_window: 8 })).equal('unlikely');
    // Without positions (approximate alignments) the stored mm_3p is used.
    should(classify.classifyAmplicon({ mm: 3, mm_3p: 2, mm_pos: null, terminal_mm: false }, perfect, { three_prime_window: 3 })).equal('unlikely');
  });

  it('exposes helpers', () => {
    should(classify.isBlocked(aln([4, 2]))).be.true();
    should(classify.isIgnored(aln([1, 2, 3, 4, 5, 6]))).be.true();
    should(classify.mismatches3p(aln([6, 5, 1]), 5)).equal(2);
    should(classify.isRepetitive(6, { repeat_site_threshold: 5 })).be.true();
    should(classify.isRepetitive(5, { repeat_site_threshold: 5 })).be.false();
    should(classify.isRepetitive(6)).be.true();
    should(classify.countsAsAmplicon('likely')).be.true();
    should(classify.countsAsAmplicon('likely_weak')).be.true();
    should(classify.countsAsAmplicon('unlikely')).be.false();
    should(classify.countsAsAmplicon('ignored')).be.false();
    should(classify.bestLikelihood(['unlikely', 'likely_weak'])).equal('likely_weak');
    should(classify.bestLikelihood([])).be.null();
  });
});

describe('amplicons.finalize', () => {
  const siteOf = (primer, subject, face, p5) => {
    const s = mk(primer, subject, face, p5);
    s.key = sites.siteKey(primer, subject, face, p5);
    return s;
  };
  const a = (p5, mmPos, extra) => Object.assign({
    p5, mm: mmPos.length, mm_pos: mmPos, mm_3p: mmPos.filter((d) => d <= 5).length, terminal_mm: mmPos.includes(1), gaps: 0, approx: false
  }, extra);

  it('uses refined coordinates, maps left/right by orientation and classifies', () => {
    const lf = siteOf(L, '4', 'F', 1000);
    const rr = siteOf(R, '4', 'R', 1209);
    const rf = siteOf(R, '5', 'F', 66890615);
    const lr = siteOf(L, '5', 'R', 66890824);
    const cands = [
      { subject: '4', orientation: 'LR', start: 1000, end: 1209, size: 210, fwd: lf, rev: rr },
      { subject: '5', orientation: 'RL', start: 66890615, end: 66890824, size: 210, fwd: rf, rev: lr }
    ];
    const alns = new Map([
      [lf.key, a(1001, [])],
      [rr.key, a(1210, [])],
      [rf.key, a(66890615, [15, 10])],
      [lr.key, a(66890824, [10, 7])]
    ]);
    const out = amplicons.finalize(cands, alns, {});
    should(out.counts).match({ candidates: 2, likely: 2, likely_weak: 0, unlikely: 0, ignored: 0, unaligned: 0 });
    should(out.amplicons).eql([
      { region: '4', start: 1001, end: 1210, size: 210, strand: 1, orientation: 'LR', likelihood: 'likely',
        left_mm: 0, right_mm: 0, left_3p_mm: 0, right_3p_mm: 0, left_mm_pos: [], right_mm_pos: [], terminal_mismatch: false, approx: false },
      { region: '5', start: 66890615, end: 66890824, size: 210, strand: -1, orientation: 'RL', likelihood: 'likely',
        left_mm: 2, right_mm: 2, left_3p_mm: 0, right_3p_mm: 0, left_mm_pos: [10, 7], right_mm_pos: [15, 10], terminal_mismatch: false, approx: false }
    ]);
  });

  it('drops ignored, lists unlikely separately, counts unaligned and dedupes', () => {
    const s = (face, p5, primer) => siteOf(primer || L, 'c', face, p5);
    const f1 = s('F', 100);
    const f2 = s('F', 101); // a second HSP that refines to the same start
    const r1 = s('R', 400, R);
    const f3 = s('F', 1000);
    const r3 = s('R', 1300, R);
    const f4 = s('F', 2000);
    const r4 = s('R', 2300, R);
    const f5 = s('F', 3000);
    const r5 = s('R', 3300, R);
    const lost = s('F', 5000);
    const cands = [
      { subject: 'c', orientation: 'LR', fwd: f1, rev: r1 },
      { subject: 'c', orientation: 'LR', fwd: f2, rev: r1 },
      { subject: 'c', orientation: 'LR', fwd: f3, rev: r3 },
      { subject: 'c', orientation: 'LR', fwd: f4, rev: r4 },
      { subject: 'c', orientation: 'LR', fwd: f5, rev: r5 },
      { subject: 'c', orientation: 'LR', fwd: lost, rev: r5 }
    ];
    const alns = new Map([
      [f1.key, a(100, [1])], // terminal → likely_weak
      [f2.key, a(100, [])], // same coordinates, better → kept
      [r1.key, a(400, [])],
      [f3.key, a(1000, [4, 2])], // blocked → unlikely
      [r3.key, a(1300, [])],
      [f4.key, a(2000, [20, 18, 16, 14, 12, 10])], // ignored
      [r4.key, a(2300, [])],
      [f5.key, a(3000, [12], { approx: true })],
      [r5.key, a(3300, [1])]
    ]);
    const out = amplicons.finalize(cands, (site) => alns.get(site.key), {});
    should(out.counts).match({ candidates: 6, likely: 1, likely_weak: 1, unlikely: 1, ignored: 1, unaligned: 1, duplicates: 1 });
    should(out.amplicons.map((x) => [x.start, x.likelihood, x.approx])).eql([[100, 'likely', false], [3000, 'likely_weak', true]]);
    should(out.amplicons[1]).match({ terminal_mismatch: true, right_mm_pos: [1], right_3p_mm: 1 });
    should(out.unlikely.map((x) => [x.start, x.left_3p_mm])).eql([[1000, 2]]);
  });

  it('reports the forward-facing site as left for single-primer products', () => {
    const f = siteOf(R, 'c', 'F', 10);
    const r = siteOf(R, 'c', 'R', 110);
    const out = amplicons.finalize([{ subject: 'c', orientation: 'RR', fwd: f, rev: r }],
      new Map([[f.key, a(10, [9])], [r.key, a(110, [])]]), {});
    should(out.amplicons[0]).match({ orientation: 'RR', strand: null, left_mm_pos: [9], right_mm_pos: [], size: 101 });
  });

  it('drops candidates whose refined ends cross', () => {
    const f = siteOf(L, 'c', 'F', 100);
    const r = siteOf(R, 'c', 'R', 120);
    const out = amplicons.finalize([{ subject: 'c', orientation: 'LR', fwd: f, rev: r }],
      new Map([[f.key, a(122, [])], [r.key, a(120, [])]]), {});
    should(out.counts.invalid).equal(1);
    should(out.amplicons).eql([]);
  });

  it('drops products whose re-aligned primer footprints overlap (counted as invalid)', () => {
    const f = siteOf(L, 'c', 'F', 100);
    const r1 = siteOf(R, 'c', 'R', 139);
    const r2 = siteOf(R, 'c', 'R', 138);
    const one = (rev, fa, ra) => amplicons.finalize([{ subject: 'c', orientation: 'LR', fwd: f, rev }], new Map([[f.key, fa], [rev.key, ra]]), {});
    should(one(r1, a(100, [], { p3: 119 }), a(139, [], { p3: 120 })).amplicons).match([{ start: 100, end: 139, size: 40 }]);
    should(one(r2, a(100, [], { p3: 119 }), a(138, [], { p3: 119 })).counts).match({ invalid: 1, likely: 0 });
    // a gapped forward alignment one base longer now touches the reverse footprint
    should(one(r1, a(100, [12], { p3: 120, gaps: 1 }), a(139, [], { p3: 120 })).counts.invalid).equal(1);
    // without p3 the site lengths give the footprints
    should(one(r2, a(100, []), a(138, [])).counts.invalid).equal(1);
    should(one(r1, a(100, []), a(139, [])).counts.likely).equal(1);
  });
});

describe('offline pipeline on captured sorghum_bicolor hits', () => {
  const QUERY = DOC.query_order.map((k) => P[k]);
  const windows = Object.values(DOC.windows);
  const fetchWindow = async (fastaPath, region, start1, end1) => {
    const w = windows.find((x) => x.subject === region && start1 >= x.start && end1 <= x.end);
    if (!w) throw Object.assign(new Error('outside fixture'), { code: 'REGION_OUT_OF_BOUNDS' });
    return w.seq.slice(start1 - w.start, end1 - w.start + 1);
  };

  async function pipeline(pairs) {
    const store = new sites.SiteStore({ primers: QUERY, totalBases: 708735318 });
    for (const line of readLines('p2_genome_loci.tsv')) store.addHit(blast.parseLine(line, 'genome'));
    const params = classify.withDefaults({});
    const cands = pairs.map((p) => amplicons.candidates(p, store, { maxSize: params.max_product_size, ignoreMismatches: params.ignore_mismatches }));
    const uniq = amplicons.uniqueSites(cands.map((c) => c.candidates));
    const alns = await realign.realignGenomeSites(uniq.sites, { fetchWindow, fastaPath: 'fixture', threePrimeWindow: params.three_prime_window });
    const byKey = realign.indexAlignments(uniq.sites, alns);
    return { store, results: cands.map((c) => amplicons.finalize(c.candidates, byKey, params)) };
  }

  it('reproduces the P2 on-target and both off-targets', async () => {
    const { store, results: [p2] } = await pipeline([{ id: 'P2', left: P.P2_L, right: P.P2_R }]);
    should(store.primerSummary(P.P2_L)).eql({ len: 22, near_perfect_sites: 2, repetitive: false, truncated: false });
    const brief = p2.amplicons.map((x) => ({ at: x.region + ':' + x.start + '-' + x.end, o: x.orientation, l: x.left_mm_pos, r: x.right_mm_pos, c: x.likelihood }));
    should(brief).eql([
      { at: '4:7423537-7423746', o: 'LR', l: [], r: [], c: 'likely' },
      { at: '4:7437317-7437526', o: 'LR', l: [], r: [], c: 'likely' },
      { at: '5:66890615-66890824', o: 'RL', l: [10, 7], r: [15, 10], c: 'likely' }
    ]);
    should(p2.amplicons.every((x) => x.size === 210 && !x.approx)).be.true();
  });

  it('finds m20_7_14 with 2 mismatches and del18 with one gap', async () => {
    const { results: [m20, del18] } = await pipeline([
      { id: 'm20', left: P.m20_7_14, right: P.P2_R },
      { id: 'del18', left: P.del18, right: P.P2_R }
    ]);
    should(m20.amplicons.find((x) => x.region === '4' && x.start === 7423537)).match({ end: 7423746, orientation: 'LR', left_mm: 2, left_mm_pos: [14, 7], right_mm: 0, likelihood: 'likely' });
    should(del18.amplicons.find((x) => x.region === '4' && x.start === 7423537)).match({ end: 7423746, orientation: 'LR', left_mm: 1, left_mm_pos: [5], left_3p_mm: 1, likelihood: 'likely' });
  });
});

// ---------------------------------------------------------------------------
// Real data (opt-in): runs blastn on sorghum_bicolor and reads its bgzip FASTA.
// PRIMERS_REALDATA=1 node --test test/primers/unit/amplicons_classify.test.js
// ---------------------------------------------------------------------------
const REAL = process.env.PRIMERS_REALDATA === '1';

describe('real data: sorghum_bicolor BLAST pipeline', { skip: REAL ? false : 'set PRIMERS_REALDATA=1' }, () => {
  const BLASTN = process.env.PRIMERS_BLASTN || '/home/olson/bin/blastn';
  const BLASTDBCMD = process.env.PRIMERS_BLASTDBCMD || '/home/olson/bin/blastdbcmd';
  const DIR = '/scratch/olson/fasta/sorghum_bicolor';
  const DNA_DB = DIR + '/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel';
  const CDNA_DB = DIR + '/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.cdna.all';
  const FASTA = DIR + '/dna/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel.fa.gz';
  const THREADS = Math.max(1, Math.min(8, Number(process.env.PRIMERS_BLAST_THREADS) || 4));
  const QUERY = DOC.query_order.map((k) => P[k]);
  const params = classify.withDefaults({});
  const timings = {};
  const shared = {};

  function faiInfo() {
    const rows = fs.readFileSync(FASTA + '.fai', 'utf8').split('\n').filter(Boolean);
    return { num_sequences: rows.length, total_bases: rows.reduce((a, r) => a + Number(r.split('\t')[1]), 0) };
  }

  function fastaFetch() {
    const { BgzipIndexedFasta } = require('@gmod/indexedfasta');
    const fa = new BgzipIndexedFasta({ path: FASTA, faiPath: FASTA + '.fai', gziPath: FASTA + '.gzi' });
    return async (fastaPath, region, start1, end1) => fa.getSequence(region, start1 - 1, end1);
  }

  it('P2, m20_7_14 and del18 at word size 5: blast → sites → amplicons → realign → classify', { timeout: 600000 }, async (t) => {
    const fai = faiInfo();
    should(fai.total_bases).equal(708735318);
    const store = new sites.SiteStore({ primers: QUERY, totalBases: fai.total_bases });
    let t0 = Date.now();
    const run = await blast.runBlast({
      cmd: BLASTN,
      args: blast.buildArgs({ target: 'genome', wordSize: 5, threads: THREADS, db: DNA_DB, maxTargetSeqs: fai.num_sequences }),
      primers: QUERY,
      onHit: (h) => store.addHit(h)
    });
    timings.blast_ws5_4primers_ms = Date.now() - t0;
    timings.hits = run.hits;
    timings.sites = store.totalSites();
    should(run.hits).be.above(50000);

    // m20_7_14 has a site at 4:7423537 at word size 5.
    should(store.sites(P.m20_7_14).some((s) => s.subject === '4' && s.face === 'F' && s.p5 === 7423537 && s.mm === 2)).be.true();
    // del18 keeps its q1-17 site (t3 = 4).
    should(store.sites(P.del18).some((s) => s.subject === '4' && s.p5 === 7423537 && s.t3 === 4)).be.true();
    should(store.primerSummary(P.P2_L, params)).eql({ len: 22, near_perfect_sites: 2, repetitive: false, truncated: false });

    const pairs = [
      { id: 'P2', left: P.P2_L, right: P.P2_R },
      { id: 'm20_7_14', left: P.m20_7_14, right: P.P2_R },
      { id: 'del18', left: P.del18, right: P.P2_R }
    ];
    t0 = Date.now();
    const cands = pairs.map((p) => amplicons.candidates(p, store, { maxSize: params.max_product_size, ignoreMismatches: params.ignore_mismatches }));
    timings.candidates_ms = Date.now() - t0;
    timings.candidates = cands.map((c) => c.candidates.length);
    cands.forEach((c) => should(c.truncated).be.false());
    const uniq = amplicons.uniqueSites(cands.map((c) => c.candidates));
    should(uniq.truncated).be.false();
    t0 = Date.now();
    const alns = await realign.realignGenomeSites(uniq.sites, { fetchWindow: fastaFetch(), fastaPath: FASTA, threePrimeWindow: params.three_prime_window });
    timings.realign_ms = Date.now() - t0;
    timings.realigned_sites = uniq.sites.length;
    should(alns.filter((x) => x.approx).length).equal(0);
    const byKey = realign.indexAlignments(uniq.sites, alns);
    t0 = Date.now();
    const [p2, m20, del18] = cands.map((c) => amplicons.finalize(c.candidates, byKey, params));
    timings.finalize_ms = Date.now() - t0;
    timings.p2_counts = p2.counts;

    const at = (res, region, start, end, orientation) => res.amplicons.find((x) => x.region === region && x.start === start && x.end === end && x.orientation === orientation);
    should(at(p2, '4', 7423537, 7423746, 'LR')).match({ left_mm: 0, right_mm: 0, likelihood: 'likely' });
    should(at(p2, '4', 7437317, 7437526, 'LR')).match({ left_mm: 0, right_mm: 0, likelihood: 'likely' });
    should(at(p2, '5', 66890615, 66890824, 'RL')).match({ left_mm_pos: [10, 7], right_mm_pos: [15, 10], left_3p_mm: 0, right_3p_mm: 0, likelihood: 'likely' });
    should(at(m20, '4', 7423537, 7423746, 'LR')).match({ left_mm: 2, left_mm_pos: [14, 7], likelihood: 'likely' });
    should(at(del18, '4', 7423537, 7423746, 'LR')).match({ left_mm: 1, left_mm_pos: [5], likelihood: 'likely' });
    shared.sites = uniq.sites;
    shared.alns = alns;
    t.diagnostic('timings ' + JSON.stringify(timings));
  });

  it('re-alignment through sequence.js (fetch + regionLength) matches direct indexedfasta reads', { timeout: 600000 }, async (t) => {
    should(shared.sites).be.ok();
    require('../../../api/helpers/primers/node_compat');
    const sequence = require('../../../api/helpers/primers/sequence');
    const t0 = Date.now();
    const alns = await realign.realignGenomeSites(shared.sites, {
      fetchWindow: (fp, region, start1, end1, strand) => sequence.fetch(fp, region, start1, end1, strand, { maxLength: 2000000 }),
      regionLength: sequence.regionLength,
      fastaPath: FASTA,
      threePrimeWindow: params.three_prime_window
    });
    should(alns).eql(shared.alns);
    t.diagnostic('realign_via_sequence_js_ms ' + (Date.now() - t0) + ' sites ' + shared.sites.length);
  });

  it('m20_7_14 is missed at word size 7', { timeout: 600000 }, async (t) => {
    const fai = faiInfo();
    const store = new sites.SiteStore({ primers: [P.m20_7_14], totalBases: fai.total_bases });
    const t0 = Date.now();
    const run = await blast.runBlast({
      cmd: BLASTN,
      args: blast.buildArgs({ target: 'genome', wordSize: 7, threads: THREADS, db: DNA_DB, maxTargetSeqs: fai.num_sequences }),
      primers: [P.m20_7_14],
      onHit: (h) => store.addHit(h)
    });
    should(run.hits).be.above(0);
    should(store.stats(0).hits).equal(run.hits);
    should(store.sites(0).some((s) => s.subject === '4' && Math.abs(s.p5 - 7423537) <= 3)).be.false();
    t.diagnostic('blast_ws7_m20_ms ' + (Date.now() - t0) + ' hits ' + run.hits);
  });

  it('cDNA target: P2 on transcripts with qseq/sseq positions', { timeout: 600000 }, async (t) => {
    let t0 = Date.now();
    const info = await blast.dbInfo({ cmd: BLASTDBCMD, db: CDNA_DB });
    should(info.num_sequences).equal(47110);
    const primers = [P.P2_L, P.P2_R];
    const store = new sites.SiteStore({ primers, totalBases: info.total_bases, keepAlignments: true });
    await blast.runBlast({
      cmd: BLASTN,
      args: blast.buildArgs({ target: 'cdna', wordSize: 5, threads: THREADS, db: CDNA_DB, maxTargetSeqs: info.num_sequences }),
      primers,
      onHit: (h) => store.addHit(h)
    });
    const blastMs = Date.now() - t0;
    t0 = Date.now();
    const c = amplicons.candidates({ id: 'P2', left: P.P2_L, right: P.P2_R }, store, { maxSize: params.max_product_size, ignoreMismatches: params.ignore_mismatches });
    const uniq = amplicons.uniqueSites(c.candidates);
    const alns = uniq.sites.map((s) => realign.realignCdnaHit(s, { threePrimeWindow: params.three_prime_window }));
    const out = amplicons.finalize(c.candidates, realign.indexAlignments(uniq.sites, alns), params);
    const byTx = (id) => out.amplicons.find((x) => x.region === id && x.orientation === 'LR');
    for (const id of ['SORBI_3004G087700.1', 'SORBI_3004G087700.2', 'SORBI_3004G087700.3', 'SORBI_3004G087800.1']) {
      should(byTx(id)).match({ left_mm: 0, right_mm: 0, size: 210, likelihood: 'likely', approx: false });
    }
    for (const id of ['SORBI_3005G183900.1', 'SORBI_3005G183900.2']) {
      should(byTx(id)).match({ left_mm_pos: [10, 7], right_mm_pos: [15, 10], likelihood: 'likely', approx: false });
    }
    t.diagnostic('cdna blast_ms ' + blastMs + ' pipeline_ms ' + (Date.now() - t0) + ' candidates ' + c.candidates.length + ' amplicons ' + out.amplicons.length);
  });
});

describe('offline pipeline on split-HSP indel hits (chk-indel-parse-filter)', () => {
  // Real blastn-short output for P2_L / L26 against sorghum_bicolor 4:7423301-7424000 with a planted
  // insertion (ins22, reverse complement ins22rc) or deletion (del26); see indel_split.json.
  const IND = JSON.parse(fs.readFileSync(path.join(FIX, 'indel_split.json'), 'utf8'));
  const IQ = IND.query_order.map((k) => IND.primers[k]);
  const fetchWindow = async (fastaPath, region, start1, end1) => {
    const seq = IND.subjects[region];
    if (!seq || start1 < 1 || end1 > seq.length) throw Object.assign(new Error('outside fixture'), { code: 'REGION_OUT_OF_BOUNDS' });
    return seq.slice(start1 - 1, end1);
  };

  it('the joined gapped site pairs into a likely one-edit product on both strands', async () => {
    const store = new sites.SiteStore({ primers: IQ });
    for (const l of IND.genome_hits) store.addHit(blast.parseLine(l, 'genome'));
    const params = classify.withDefaults({});
    const pairs = [{ id: 'P2', left: IQ[0], right: IQ[1] }, { id: 'L26', left: IQ[2], right: IQ[1] }];
    const cands = pairs.map((p) => amplicons.candidates(p, store, { maxSize: params.max_product_size, ignoreMismatches: params.ignore_mismatches }));
    const plan = amplicons.realignPlan(cands.map((c) => c.candidates), { maxAmplifyingMismatches: params.max_amplifying_mismatches });
    const alns = await realign.realignGenomeSites(plan.realign, { fetchWindow, fastaPath: 'fixture', regionLength: async (fp, r) => IND.subjects[r].length });
    should(alns.every((x) => !x.approx)).be.true();
    const byKey = realign.indexAlignments(plan.realign, alns);
    for (const s of plan.bounded) byKey.set(s.key, realign.realignFromHit(s, params.three_prime_window, { mm: amplicons.realignBound(s) }));
    const [p2, l26] = cands.map((c) => amplicons.finalize(c.candidates, byKey, params));
    const at = (res, region, start, end, o) => res.amplicons.find((x) => x.region === region && x.start === start && x.end === end && x.orientation === o);
    should(at(p2, 'ins22', 237, 447, 'LR')).match({ size: 211, left_mm: 1, left_mm_pos: [12], right_mm: 0, likelihood: 'likely', approx: false });
    should(at(p2, 'ins22rc', 255, 465, 'RL')).match({ size: 211, left_mm: 1, right_mm: 0, likelihood: 'likely', approx: false });
    // P2_L q1-13 (t3 9) passed the parse filter already; its lbU 5 exceeds the cap but it must still be re-aligned
    should(at(p2, 'del26', 237, 445, 'LR')).match({ left_mm: 1, likelihood: 'likely', approx: false });
    should(at(l26, 'ins22', 237, 447, 'LR')).match({ left_mm: 1, likelihood: 'likely' });
    should(at(l26, 'ins22rc', 255, 465, 'RL')).match({ left_mm: 1, likelihood: 'likely' });
    should(at(l26, 'del26', 237, 445, 'LR')).match({ size: 209, left_mm: 1, likelihood: 'likely' });
    should(store.gappedSites).equal(5);
  });
});
