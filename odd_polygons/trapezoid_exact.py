"""[Deposit] Helper, imported by the scripts of this directory. Paper: "Perpendicular billiards in rational polygons of odd denominator". Research name: probes/s931_trapezoid_exact.py, copied verbatim except imports (deposit names), paths (data/ and figs/ under odd_polygons/) and this paragraph. Run through run_all.py; the usage lines below are the research repo's. Bracketed IDs and section marks (§sNNN) refer to the project's internal ledgers.

s931 keeper (queue `s930-a`): (M⁺) at EVERY side length for the two bounded N = 7 quadrilateral families `(1,2,5,6)/7`,
`(1,3,4,6)/7` — an EXACT map of the change set (every prong word, every breakpoint) on a range of the slice.

Both families are TRAPEZOIDS (sides 0 ∥ 2: a_0 + a_3 = a_1 + a_2 = N): the heptagonal triangles (1,2,4)/7 and (1,3,3)/7 cut parallel to
the side V0V1 at leg length t (slice L = (1, t, ·, ·), t ∈ (0, 1/(1 + 2cos(2π/7))) — the cut side → 0 at the right end, the strip
gets thin at the left end).

PRE-REGISTERED (s931, before any range was tracked).
  Claim under test (T): on a compact range [T0, T1] inside the slice, (i) every slot's prong has finitely many words, tracked exactly
  (`poly_exact.Fam.track`: each word valid on the open interval between the nearest roots of its vertex forms); (ii) no prong
  caps; (iii) at every t ∈ [T0, T1] — every open piece AND every breakpoint — the vertex-ending prongs pair by reversal into p ≤ 1
  swapped pairs, or into p ≥ 2 pairs with no signed relation among their lengths (exact). Then Theorem P ([POLY-THM-P]) gives
  (M⁺) ⊕ CP at every t of the range. Refuted by a breakpoint carrying two pairs whose lengths satisfy a signed relation, by a
  persistent pair set with an identical relation, or by a capped prong.
  Expected from s930 (measured, float): every t* carries one pair (`no coincident t*`), max prong steps 63 on the central 80%;
  so p ≤ 1 everywhere on the central range. The ends are NOT covered by a compact range: near t = 0 the strip transit costs
  ~1/t bounces (breakpoints expected to accumulate), near t_max the cut side shrinks to the apex of the triangle — handled
  separately (the ledger section).
  Could come out otherwise: a coincident breakpoint (two independent alignment equations sharing a root in ℚ(cos π/14)) — the
  float scan resolved t* only to 1e-20 and checked pair COUNTS, not exact coincidence; or a persistent vertex prong (an identically
  aligned pair) on a piece.
  Control: the tracker reproduces `Poly.prong` (float, s910) on all 14 slots at sample t (checked in `data/scratch/s931/t_check.py`),
  and every breakpoint is re-derived as an exact vertex configuration (`bp`), its vertex prong counted.
PRIOR ART: §s930 (2), (3) ⊕ [NCYL-721] (b) (the two families: no coincident t*, bounded prong length 63 on the central 80%, float);
  §s928 (2) [POLY-ZERO-ONE] (change set ⊆ Z₁, countable); §s928 (3) truncations of triangles (vertex cuts, not parallel cuts);
  `rulings.py --grep` *change set* → [NCYL-721] only; *1,2,5,6* → [NCYL-721]; *heptagonal* → none; *trapezoid* → (grepped s931).
DATA: `data_index.py --grep 1256` → only `s440_claim_store_summary` keys (unrelated); no exact polygon store.
Output: data/s931_trapezoid_<fam>_<T0>_<T1>.json.  Usage: python probes/s931_trapezoid_exact.py FAM T0 T1   (FAM 1256 | 1346; T rationals)
"""
import _setup  # noqa: F401  (deposit: this directory on the path, cwd = odd_polygons/)
import json
import sys
import time
from fractions import Fraction as Fr

from poly_exact import Fam

FAMS = {'1256': [1, 2, 5, 6], '1346': [1, 3, 4, 6]}


def run(tag, T0, T1, cap=20000):
    if tag in FAMS:
        a, NN = FAMS[tag], 7
    else:                                   # 'N:a0,a1,a2,a3'
        nn, aa = tag.split(':')
        NN, a = int(nn), [int(x) for x in aa.split(',')]
    E = Fam.closing(a, NN, [(1, 0), (0, 1)])
    if T1 is None:                          # just below t_max = 1/c (top side 1 − ct)
        T1 = Fr(999, 1000) / Fr(str(float(-E.L[2].c1)))
        T1 = Fr(T1.numerator // 10**6 if False else T1).limit_denominator(10**6)
    sg = E.sgn
    slots = E.slots()
    per = {}
    t_start = time.time()
    for (i, d) in slots:
        per[(i, d)] = E.track(i, d, T0, T1, cap)
        assert all(p['end'][0] != 'cap' for p in per[(i, d)]), ('CAP', i, d)
    # global breakpoints
    bps = []
    for key, pcs in per.items():
        for p in pcs[:-1]:
            bps.append((p['hi'], key))
    bps.sort(key=lambda x: x[0].mp(40))
    groups = []
    for x, key in bps:
        if groups and (x - groups[-1][0]).iszero():
            groups[-1][1].append(key)
        else:
            groups.append([x, [key]])

    def piece_at(key, x):
        """the piece of `key` whose open interval contains x, or ('bp', piece) if x is its breakpoint"""
        for p in per[key]:
            lo_ok = p['lo'] is None or sg(x - p['lo']) > 0
            hi_ok = p['hi'] is None or sg(p['hi'] - x) > 0
            if lo_ok and hi_ok:
                return p
            if p['hi'] is not None and (p['hi'] - x).iszero():
                return p['bp']
        raise RuntimeError

    def pairs_of(prs):
        """vertex prongs -> pairs by reversal; returns (pairs, unpaired)"""
        vp = {}
        for key, p in prs.items():
            if p['end'][0] == 'vertex':
                vp[key] = p
        pairs, seen, unp = [], set(), []
        for key, p in vp.items():
            if key in seen:
                continue
            rev = (p['end'][1], (p['arr'] + 2 * NN) % (4 * NN))
            if rev in vp and rev not in seen:
                pairs.append((key, rev, p['length'].at(T_ := None) if False else p['length']))
                seen |= {key, rev}
            else:
                unp.append(key)
        return pairs, unp

    def relation(pairs, x):
        """exact signed relation among pair lengths at t = x?"""
        import itertools
        ls = [p[2].at(x) for p in pairs]
        for eps in itertools.product((-1, 0, 1), repeat=len(ls)):
            if any(eps) and next(e for e in eps if e) == 1:
                s = sum((e * l for e, l in zip(eps, ls)), E.K.zero)
                if s.iszero():
                    return eps
        return None
    report = {'fam': a, 'N': NN, 'T0': str(T0), 'T1': str(T1), 'slots': len(slots),
              'pieces_per_slot': {f'{k[0]},{k[1]}': len(v) for k, v in per.items()},
              'max_steps': max(p['steps'] for v in per.values() for p in v),
              'n_breakpoints': len(groups), 'bp': [], 'intervals_p': {}, 'fail': []}
    # open intervals: persistent vertex prongs only
    cuts = [None] + [g[0] for g in groups] + [None]
    for j in range(len(cuts) - 1):
        lo, hi = cuts[j], cuts[j + 1]
        lo_v = lo if lo is not None else E.K.of(T0)
        hi_v = hi if hi is not None else E.K.of(T1)
        mid = (lo_v + hi_v) * Fr(1, 2)
        prs = {key: piece_at(key, mid) for key in slots}
        pairs, unp = pairs_of(prs)
        p = len(pairs)
        report['intervals_p'][p] = report['intervals_p'].get(p, 0) + 1
        if unp:
            report['fail'].append(('unpaired on interval', float(mid), unp))
        if p >= 2:
            # identical relation? lengths affine; check at two points
            r1 = relation(pairs, mid)
            if r1:
                report['fail'].append(('relation on interval', float(mid), r1))
    for x, keys in groups:
        prs = {key: piece_at(key, x) for key in slots}
        pairs, unp = pairs_of(prs)
        rel = relation(pairs, x) if len(pairs) >= 2 else None
        row = {'t': float(x), 't_exact': [str(c) for c in x.c], 'slots_changing': [list(k) for k in keys], 'p': len(pairs),
               'lengths': [float(pp[2].at(x)) for pp in pairs], 'relation': rel, 'unpaired': [list(k) for k in unp],
               'steps': max(prs[k]['steps'] for k in keys)}
        report['bp'].append(row)
        if unp or rel:
            report['fail'].append(('breakpoint', float(x), unp, rel))
    report['p_at_breakpoints'] = {}
    for r in report['bp']:
        report['p_at_breakpoints'][r['p']] = report['p_at_breakpoints'].get(r['p'], 0) + 1
    report['seconds'] = round(time.time() - t_start, 1)
    return report


if __name__ == '__main__':
    tag, T0 = sys.argv[1], Fr(sys.argv[2])
    T1 = None if sys.argv[3] == 'max' else Fr(sys.argv[3])
    rep = run(tag, T0, T1)
    out = f'data/s931_trapezoid_{tag.replace(":", "_").replace(",", "")}_{sys.argv[2].replace("/", "_")}_{sys.argv[3].replace("/", "_")}.json'
    json.dump(rep, open(out, 'w'), indent=1, default=str)
    print(tag, 'range', T0, T1, 'breakpoints', rep['n_breakpoints'], 'p at bps', rep['p_at_breakpoints'],
          'intervals p', rep['intervals_p'], 'max steps', rep['max_steps'], 'fail', rep['fail'][:5], rep['seconds'], 's')
