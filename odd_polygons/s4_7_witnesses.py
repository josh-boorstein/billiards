"""[Deposit] Example 4.7 (a), (b): the two asymmetric free-cylinder witnesses. Paper: "Perpendicular billiards in rational polygons of odd denominator". Research name: probes/s929_asym_witness.py, copied verbatim except imports (deposit names), paths (data/ and figs/ under odd_polygons/) and this paragraph. Run through run_all.py; the usage lines below are the research repo's. Bracketed IDs and section marks (§sNNN) refer to the project's internal ledgers.

s929 keeper (queue `s928-d`): WITNESS — free area at odd N in a quadrilateral with NO polygon symmetry.

Input: angles a/N (all four distinct ⟹ no kite / parallelogram / isosceles-trapezoid symmetry) and an approximate t*
from `probes/s929_asym_free.py` (family L = (1, t, L2(t), L3(t))). Steps, all at DPS digits:
  (1) t* resolved CHEAPLY: a prong whose outcome flips across t* ('ends perpendicular within 12 bounces' vs not) is
      bisected with a 50-bounce cap — no capped walk in the loop (the walk-capped bisection of s929_asym_free is slow
      and lands where the leak time equals the cap, not at t*);
  (2) the class's prongs re-walked with a NEAR-VERTEX tolerance NEAR (the step-(1) bisection stops at the edge of
      the vertex-tolerance window, so connections miss their vertex by ~the tolerance and the instrument's own test misses them —
      the s929 p = 0 artefact; NOT a geometry floor):
      every vertex-ending prong, its reverse slot, word and LENGTH — the swapped pairs and their signed relation;
  (3) the area test R (`s928_poly_thmP.row`, c = 1);
  (4) an independent Monte Carlo of the free share (uniform point in the table x uniform class direction = uniform on S;
      free = no perpendicular collision within 400 bounces; reports the largest steps-to-perpendicular among the
      non-free samples, which must sit far below the cap).
  (5) R at t* ± 1e-6 (isolation).
PRIOR ART: §s928 (2), (4), [NCYL-719] (*a free region bounded by connections from two distinct swapped pairs (a slit,
or equal-length loops)* — the residual this witnesses), `probes/s929_asym_free.py` (the scan that found t*).
DATA: `data/s929_asym_free.json` (scan rows; t* to ~1e-12). Output: data/s929_asym_witness.json.
Usage: python probes/s929_asym_witness.py  (runs the WITNESSES list; ~3 min).
"""
import _setup  # noqa: F401  (deposit: this directory on the path, cwd = odd_polygons/)
import collections
import json
import random
import sys

import mpmath as mp

import poly_slots as S
import area_test as T
import s4_7_scan as M

DPS = 130
WITNESSES = [([2, 4, 5, 7], 9, '0.095459452979'), ([4, 5, 9, 12], 15, '0.54070669528')]


def resolve(a, N, tc, w=mp.mpf('1e-9')):
    lo, hi = tc - w, tc + w
    P0, P1 = M.make(a, N, lo), M.make(a, N, hi)
    for key in P0.slots(P0.classes()[0]):
        r0, r1 = P0.prong(*key, 50), P1.prong(*key, 50)
        q0, q1 = r0[0] == 'perp' and r0[3] <= 12, r1[0] == 'perp' and r1[3] <= 12
        if q0 != q1:
            break
    else:
        return None, None
    for _ in range(500):
        m = (lo + hi) / 2
        r = M.make(a, N, m).prong(*key, 50)
        if (r[0] == 'perp' and r[3] <= 12) == q0:
            lo = m
        else:
            hi = m
        if hi - lo < mp.mpf(10) ** -(DPS - 15):
            break
    return (lo + hi) / 2, key


def walk(P, i, d, near, cap=400):
    N, Q, L, word = P.N, P.V[i], mp.mpf(0), []
    for st in range(cap):
        h = P.first_hit(Q, d, None)
        L += T.dist(Q, h[2])
        if h[0] == 'vertex':
            return ('vertex', h[1], d, st + 1, L, word)
        _, j, X, s = h
        if s < near or s > 1 - near:
            return ('vertex', j if s < near else (j + 1) % P.k, d, st + 1, L, word)
        word.append(j)
        if (d - P.m[j] - N) % (2 * N) == 0:
            return ('perp', j, d, st + 1, L, word)
        d = (2 * P.m[j] - d) % (4 * N)
        Q = X
    return ('cap', None, None, cap, L, word)


def pairs(P, near):
    N = P.N
    C = P.classes()[0]
    res = {k: walk(P, *k, near) for k in P.slots(C)}
    out = []
    for k, r in sorted(res.items()):
        if r[0] != 'vertex':
            continue
        rev = (r[1], (r[2] + 2 * N) % (4 * N))
        if rev != k and rev in res and res[rev][0] == 'vertex' and k < rev:
            out.append({'slot': k, 'rev': rev, 'end_vertex': r[1], 'bounces': r[3], 'length': mp.nstr(r[4], 30),
                        'word': r[5], 'loop': k[0] == r[1]})
    return out, dict(collections.Counter(r[0] for r in res.values()))


def monte_carlo(P, n, seed=929, cap=400):
    N, C = P.N, sorted(P.classes()[0])
    xs, ys = [v[0] for v in P.V], [v[1] for v in P.V]

    def inside(x, y):
        c = False
        for i in range(P.k):
            (x1, y1), (x2, y2) = P.V[i], P.V[(i + 1) % P.k]
            if (y1 > y) != (y2 > y) and x < x1 + (y - y1) * (x2 - x1) / (y2 - y1):
                c = not c
        return c
    rnd, free, worst = random.Random(seed), 0, 0
    for _ in range(n):
        while True:
            x = mp.mpf(rnd.uniform(float(min(xs)), float(max(xs))))
            y = mp.mpf(rnd.uniform(float(min(ys)), float(max(ys))))
            if inside(x, y):
                break
        d, Q, hit = rnd.choice(C), (x, y), None
        for st in range(cap):
            h = P.first_hit(Q, d, None)
            if h[0] == 'vertex':
                hit = st
                break
            _, j, X, s = h
            if (d - P.m[j] - N) % (2 * N) == 0:
                hit = st
                break
            d, Q = (2 * P.m[j] - d) % (4 * N), X
        if hit is None:
            free += 1
        else:
            worst = max(worst, hit)
    p = free / n
    return {'n': n, 'free': free, 'share': p, 'sigma': (p * (1 - p) / n) ** 0.5, 'max_steps_nonfree': worst}


def main(nmc=3000):
    out = []
    for a, N, tc in WITNESSES:
        mp.mp.dps = DPS
        S.TOL = mp.mpf(10) ** -73
        ts, key = resolve(a, N, mp.mpf(tc))
        if ts is None:
            print(f'{a}/{N}: no flipping prong near {tc}')
            continue
        P = M.make(a, N, ts)
        pr, counts = pairs(P, mp.mpf(10) ** -45)
        lens = [mp.mpf(x['length']) for x in pr]
        rel = None
        if len(lens) >= 2:
            import itertools
            rel = min((abs(sum(e * l for e, l in zip(eps, lens))), eps)
                      for eps in itertools.product((-1, 0, 1), repeat=len(lens)) if any(eps) and
                      next(e for e in eps if e) == 1)
        poly = next((q for q in (mp.findpoly(ts, dg, maxcoeff=10 ** 8) for dg in range(1, 7)) if q), None)   # lowest degree (s929 audit W: deg 6 returned the square)
        mp.mp.dps = 60
        S.TOL = mp.mpf(10) ** -30
        P60 = M.make(a, N, ts)
        R = [r['one_minus_R'] for r in T.row('w', P60, cap=4000)]
        Rnb = {d: [r['one_minus_R'] for r in T.row('nb', M.make(a, N, ts + mp.mpf(d)), cap=4000)]
               for d in ('1e-6', '-1e-6')}
        mp.mp.dps = 40
        S.TOL = mp.mpf(10) ** -32
        mc = monte_carlo(M.make(a, N, ts), nmc)
        row = {'a': a, 'N': N, 't_star': mp.nstr(ts, 60), 'flip_prong': key, 'minpoly_deg<=6': poly,
               'prong_counts': counts, 'pairs': pr, 'min_signed_rel': None if rel is None else
               (mp.nstr(rel[0], 5), rel[1]), 'one_minus_R': R, 'one_minus_R_pm1e-6': Rnb, 'mc': mc}
        out.append(row)
        print(json.dumps(row, indent=1, default=str), flush=True)
    json.dump(out, open('data/s929_asym_witness.json', 'w'), indent=1, default=str)


if __name__ == '__main__':
    main(*(int(x) for x in sys.argv[1:2]))
