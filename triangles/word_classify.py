"""[Deposit] helper: strands and the affine classification of slot walks into word families. Paper: "Complete periodicity of perpendicular billiards in rational triangles of prime denominator". Research name: probes/s939_q5_families.py, copied verbatim except imports (deposit names), paths (data/ and figs/ under triangles/) and this paragraph. Run through run_all.py; the usage lines below are the research repo's. Bracketed IDs and section marks (§sNNN) refer to the project's internal ledgers.

s939 `s938-d` (i) survey: EVERY perpendicular-class slot walk of the two q = 5 families, classified into word families.

Families P1 = (3, 12m+5, 18m+7), P2 = (3, 12m+4, 18m+8), N = 30m + 15. For each m the class has N launch slots; each walk
(`s929_asym_witness.walk`, bounce cap 40m + 100 — the default 400 truncates the loops, 24m + 2k + 17 bounces, from m ≈ 16)
ends either at a vertex (a connection) or perpendicularly on a side (an ι-invariant connection, Theorem P type (2)).
MEASURED (m = 4 … 26, both families): exactly 2m slots end at a vertex — the m loops L_{m,k} and their reverses — so p = m;
all other slots end perpendicularly; no walk capped.

CLASSIFICATION. A word = runs of alternating 0/2 reflections (about V0) separated by single side-1 reflections. At each start
vertex the slots, in direction order, split into three STRANDS (slot residue mod 6); t = position in the strand. Each strand
cuts into maximal segments on which the run vector moves by a constant increment (same run count, same end kind). From m = 5
the segment count is constant — 30 (P1) / 32 (P2) — and from m = 6 every segment's runs are EXACTLY affine in (m, t)
(A m + B t + C per run, fitted on m = 6 … 19, every member checked) — EXCEPT that words also change with m: P2's V0 slot 3
is [0, 6m+2, 1, 2m+6, 3m+5] for m ≤ 19 and [0, 6m+2, m] for every m = 20 … 26 (same at 40 and 80 digits) — so the m ≤ 19
structure is NOT the all-m structure; `m_transitions` lists every family whose fitted word breaks on m = 20 … 26.
The end side alternates 0/2 with the parities.
BOUNDARIES (first/last t of each family as functions of m): each is fitted as ⌊(a m + b)/q⌋ (q ≤ 12) on m = 6 … 19 and
VALIDATED on held-out m = 20 … 26 (`boundaries_failing_holdout`): the fits with q ≥ 7 are overfits and break out of sample —
those boundaries are not rational floors but transcendental curves t*(m) (e.g. P1 4|5, slope ≈ 0.43; 21|22, ≈ 5.57). A family boundary is where the walk passes through a vertex; at a
rational boundary lattice points stay off it by exact arithmetic, but at a transcendental one a lattice point ON the boundary curve
for some m would be an extra vertex connection (an extra swapped pair, p > m). So (i) for every m needs THREE things: (1) the
asymptotic families (x → 0) certified on x ≤ x(m*), m* past the last transition in m (`s938-e`'s pipeline with
perpendicular endings and floor-type ranges); (2) a rigorous per-m check for m < m*; (3) an exact non-vanishing argument at
every transcendental boundary (a cyclotomic question, like [POLY-Q5-COS]).

PRIOR ART: §s911 (2) (p = m measured, m small), §s938 (2) (the loops' word), §s939 (1) ([POLY-Q5-COS]); `rulings.py --grep`
*word family*, *slot walk* → none. DATA: `data_index.py --grep loop` → none for this family.
Output: `data/s939_q5_walks_{P1,P2}.json` (every walk: slot, kind, end, runs, first letters), `data/s939_q5_families_{P1,P2}.json`
(the family table: runs affine in (m, t), t-ranges per m, floor fits). Runtime ~ 10 min (both families, 2 processes)."""
import collections
import itertools
import json
import sys
from fractions import Fraction as F

import mpmath as mp
import sympy as sp

import _setup  # noqa: F401  (deposit: lib/ on the path, cwd = triangles/)
from poly_slots import Poly
import slot_walk as W

mp.mp.dps = 40
ABC = {'P1': lambda m: (3, 12 * m + 5, 18 * m + 7), 'P2': lambda m: (3, 12 * m + 4, 18 * m + 8)}
MS_WALK, MS_FIT = range(4, 27), range(6, 20)


def walks(fam, m):
    N, abc = 30 * m + 15, ABC[fam](m)
    P = Poly(list(abc), N, lengths=[mp.sin(mp.pi * x / N) for x in (abc[2], abc[0], abc[1])])
    rows = []
    for k in sorted(P.slots(P.classes()[0])):
        r = W.walk(P, *k, mp.mpf(10) ** -20, cap=40 * m + 100)
        parts = ''.join(map(str, r[5])).split('1')
        rows.append([list(k), r[0], str(r[1]), [len(x) for x in parts], ''.join(x[:1] or '-' for x in parts)])
    return rows


def strands(rows):
    out = {}
    for v in (0, 1, 2):
        R = sorted([r for r in rows if r[0][0] == v], key=lambda r: r[0][1])
        ds = [r[0][1] for r in R]
        step = min(b - a for a, b in zip(ds, ds[1:])) if len(ds) > 1 else 1
        for res in sorted({(d - ds[0]) // step % 3 for d in ds}):
            out[(v, res)] = [r for r in R if (r[0][1] - ds[0]) // step % 3 == res]
    return out


def segments(rows):
    segs = []
    for (v, res), st in strands(rows).items():
        cur = []
        for t, r in enumerate(st):
            if cur:
                prev = cur[-1][1]
                ok = len(r[3]) == len(prev[3]) and r[1] == prev[1]
                if ok and len(cur) >= 2:
                    ok = [a - b for a, b in zip(r[3], prev[3])] == [a - b for a, b in zip(prev[3], cur[-2][1][3])]
                if not ok:
                    segs.append((v, res, cur))
                    cur = []
            cur.append((t, r))
        if cur:
            segs.append((v, res, cur))
    return segs


def affine_fit(P):
    """exact fit val = A m + B t + C through all points (m, t, val), or None"""
    for a, b, c in itertools.combinations(range(len(P)), 3):
        M = sp.Matrix([[P[i][0], P[i][1], 1] for i in (a, b, c)])
        if M.det() != 0:
            break
    else:
        vals = {v for _, _, v in P}
        ms = sorted({m for m, _, _ in P})
        if len(ms) < 2:
            return None
        (m0, _, v0), (m1, _, v1) = min(P), max(P)
        A = sp.Rational(v1 - v0, m1 - m0)
        C = v0 - A * m0
        return (A, 0, C) if all(A * m + C == v for m, _, v in P) else None
    A, B, C = M.inv() * sp.Matrix([P[i][2] for i in (a, b, c)])
    return (A, B, C) if all(A * m + B * t + C == v for m, t, v in P) else None


def floor_fit(vals):
    for q in range(1, 13):
        for a in range(-12 * q, 12 * q + 1):
            for b in range(-40 * q, 40 * q + 1):
                if all((a * m + b) // q == v for m, v in vals.items()):
                    return (a, b, q)
    return None


def classify(fam, W_):
    S = {m: segments(W_[m]) for m in MS_WALK}
    counts = {m: len(S[m]) for m in MS_WALK}
    n = counts[MS_FIT[0]]
    assert all(counts[m] == n for m in MS_FIT), counts
    table, allaff = [], True
    for i in range(n):
        pts = collections.defaultdict(list)
        rng = {}
        for m in MS_FIT:
            v, res, cur = S[m][i]
            rng[m] = (cur[0][0], cur[-1][0])
            for t, r in cur:
                for c, x in enumerate(r[3]):
                    pts[c].append((m, t, x))
        fits = {c: affine_fit(P) for c, P in pts.items()}
        aff = all(f is not None for f in fits.values())
        allaff &= aff
        v, res, cur = S[MS_FIT[0]][i]
        table.append({'i': i, 'vertex': v, 'strand': res, 'kind': cur[0][1][1], 'first_letters': cur[0][1][4],
                      'nruns': len(cur[0][1][3]), 'affine': aff,
                      'runs': {c: [str(x) for x in f] for c, f in fits.items() if f},
                      't_lo': floor_fit({m: a for m, (a, b) in rng.items()}), 't_hi': floor_fit({m: b for m, (a, b) in rng.items()}),
                      't_range_by_m': rng})
    return counts, allaff, table


def main(fam):
    W_ = {m: walks(fam, m) for m in MS_WALK}
    json.dump(W_, open(f'data/s939_q5_walks_{fam}.json', 'w'))
    kinds = {m: dict(collections.Counter(r[1] for r in W_[m])) for m in MS_WALK}
    assert all(kinds[m].get('vertex', 0) == 2 * m and 'cap' not in kinds[m] for m in MS_WALK), kinds
    counts, allaff, table = classify(fam, W_)
    irr = [t['i'] for t in table if t['t_lo'] is None or t['t_hi'] is None]
    json.dump({'family': fam, 'm_walk': [4, 20], 'm_fit': [6, 20], 'kinds': kinds, 'segment_counts': counts,
               'all_affine': allaff, 'irrational_boundary_families': irr, 'families': table},
              open(f'data/s939_q5_families_{fam}.json', 'w'), indent=1, default=str)
    # words that break the m <= 19 fit on m = 20..26
    trans = []
    for t in table:
        for m in range(20, 27):
            seg = segments(W_[m])
            if len(seg) != len(table) or len(seg[t['i']][2][0][1][3]) != t['nruns'] or any(
                    sp.Rational(t['runs'][c][0]) * m + sp.Rational(t['runs'][c][1]) * tt + sp.Rational(t['runs'][c][2]) != x
                    for tt, r in seg[t['i']][2] for c, x in enumerate(r[3]) if c in t['runs']):
                trans.append([t['i'], m])
                break
    out = json.load(open(f'data/s939_q5_families_{fam}.json'))
    out['m_transitions'] = trans
    hold = []                                        # every floor fit (m = 6..19) checked on HELD-OUT m = 20..26
    for t in table:
        for name, j in (('t_lo', 0), ('t_hi', -1)):
            fit = t[name]
            okh = fit is not None and all(len(segments(W_[m])) == len(table) and
                                          (fit[0] * m + fit[1]) // fit[2] == segments(W_[m])[t['i']][2][j][0] for m in range(20, 27))
            if not okh:
                hold.append([t['i'], name, str(fit)])
    out['boundaries_failing_holdout'] = hold
    json.dump(out, open(f'data/s939_q5_families_{fam}.json', 'w'), indent=1, default=str)
    print(fam, '| p = m, no caps (m = 4..26) | segments', counts, '| all affine (m 6..19)', allaff, '| unfitted boundary', irr,
          '| first break on m = 20..26 [family, m]', trans, '| boundaries failing holdout', len(hold), [h[:2] for h in hold])
    assert allaff


if __name__ == '__main__':
    main(sys.argv[1])
