"""Theorem 4.2 (the slot identity), Proposition 5.1 (the band count) and Theorem 5.2 (the band bound), on every class
of every scalene rational triangle with odd N <= 21 or even N <= 34 (456 triangles, 820 classes) -- the polygon note,
Section 9, first item.

For each triangle (angles pi*(a, b, c)/N, a < b < c, gcd 1) and each class C of perpendicular directions:
  * the CLASSIFIER (`beams.classify`) gives every vertical side's branches, exactly: n_sigma, and its retracing,
    mirror and translation branches;
  * the TRACER (`slots.trace_class`) follows every slot's prong, exactly: f_sigma (prongs meeting sigma at a right
    angle first), the vertex-to-vertex prongs matched with their reverses (p pairs), and s' (prongs reaching no
    vertex within the cap).
Checked, side by side and class by class:
  (i)   Theorem 4.2 (ii), per vertical side: f_sigma = n_sigma - 1 (the singular points are the prongs that meet sigma
        at a right angle first), and S = N at odd N, S = N/2 - h at even N (Lemma 4.1);
  (ii)  Theorem 4.2: sum_{sigma vertical} n_sigma = v + S - 2p - s';
  (iii) Proposition 5.1: sum n_sigma = 2(R + M + T), R the retracing bands, M the mirror and T the translation pairs;
  (iv)  Theorem 5.2: at even N, on the class whose vertical sides meet at the vertex with even numerator,
        R + T = ceil(N/4) - M - p - s'/2; on the other class R = 0 and M + T = floor(N/4) - p - s'/2; at odd N, M = 0
        and R + T = (N + 3)/2 - p - s'/2.
The identities are not vacuous (the control): on the class of (7, 8, 15) perpendicular to the side opposite the
right angle's neighbour -- the class that is not completely periodic, found by Delecroix and Ruth [DR22] -- the tracer
finds s' = 6 and p = 1, and sum n_sigma = 15 - 2 + 1 - 2 - 6 = 6.

At even N the unfolding carries each slot in both orientations: the tracer's slot, pair and cap counts are twice the
note's (`slots.py`), and are halved here.  About five minutes.
"""
import math
import time
from math import gcd

import common
from common import check, done
from beams import classify
from slots import classes, trace_class, pair_reversals

NMAX_ODD, NMAX_EVEN = 21, 34


def scalene():
    for N in range(3, NMAX_EVEN + 1):
        if N % 2 and N > NMAX_ODD:
            continue
        for a in range(1, N):
            for b in range(a + 1, N):
                c = N - a - b
                if c > b and gcd(gcd(a, b), c) == 1:
                    yield (a, b, c)


def main():
    t0 = time.time()
    n_tri = n_cls = 0
    bad = {k: [] for k in ('classifier', 'per_side', 'S', 'identity', 'pairs', 'band_count', 'band_bound')}
    control = None
    for abc in scalene():
        N = sum(abc)
        n_tri += 1
        side = {s: classify(abc, s) for s in range(3)}
        if not all(r['ok'] for r in side.values()):
            bad['classifier'].append(abc)
            continue
        even_vertex = next((i for i in range(3) if abc[i] % 2 == 0), None) if N % 2 == 0 else None
        for c in classes(abc):
            n_cls += 1
            T, sl, geo = trace_class(abc, c['j'])
            mult = 1 if N % 2 else 2
            pairs, _ = pair_reversals(T, sl, geo)
            if pairs is None:
                bad['pairs'].append((abc, c['j']))
                continue
            p, sp = len(pairs) / mult, sum(r['end'] == 'cap' for r in sl) / mult
            S = len(sl) / mult
            v, h = len(c['vert']), len(c['hor'])
            f = {}
            for r in sl:
                if r['end'] == 'perp':
                    f[r['side']] = f.get(r['side'], 0) + 1
            n = {s: side[s]['n'] for s in c['vert']}
            if any(f.get(s, 0) != mult * (n[s] - 1) for s in c['vert']):
                bad['per_side'].append((abc, c['j']))
            if S != (N if N % 2 else N // 2 - h):
                bad['S'].append((abc, c['j']))
            sum_n = sum(n.values())
            if sum_n != v + S - 2 * p - sp:
                bad['identity'].append((abc, c['j']))
            R = sum(side[s]['R'] for s in c['vert']) / 2
            M = sum(side[s]['M'] for s in c['vert']) / 2
            Tp = sum(side[s]['T'] for s in c['vert']) / 2
            if sum_n != 2 * (R + M + Tp):
                bad['band_count'].append((abc, c['j']))
            if N % 2:
                ok = M == 0 and R + Tp == (N + 3) / 2 - p - sp / 2
            elif v == 2:
                at_even = {0: (0, 2), 1: (0, 1), 2: (1, 2)}[even_vertex]      # the two sides through V_k
                ok = set(c['vert']) == set(at_even) and R + Tp == math.ceil(N / 4) - M - p - sp / 2
            else:
                ok = R == 0 and M + Tp == N // 4 - p - sp / 2
            if not ok:
                bad['band_bound'].append((abc, c['j']))
            if abc == (7, 8, 15) and v == 1:
                control = (sp, p, sum_n, v, S)
    print(f'{n_tri} triangles, {n_cls} classes, {time.time() - t0:.0f} s')
    check(f'456 triangles and 820 classes, every side classified (exact tiling of every entry edge)',
          (n_tri, n_cls) == (456, 820) and not bad['classifier'])
    check('every vertex-to-vertex prong matched with its reverse', not bad['pairs'])
    check('Theorem 4.2 (ii): f_sigma = n_sigma - 1 on every vertical side of every class', not bad['per_side'])
    check('Lemma 4.1: S = N at odd N, S = N/2 - h at even N', not bad['S'])
    check('Theorem 4.2: sum n_sigma = v + S - 2p - s\' on all 820 classes', not bad['identity'])
    check('Proposition 5.1: sum n_sigma = 2(R + M + T) on all 820 classes', not bad['band_count'])
    check('Theorem 5.2: the band bound with its exact defect on all 820 classes', not bad['band_bound'])
    check(f'control (7, 8, 15), the class that is not completely periodic: s\' = 6, p = 1, '
          f'sum n_sigma = 15 - 2 + 1 - 2 - 6 = 6 (got s\', p, sum = {control[:3] if control else None})',
          control is not None and control[:3] == (6, 1, 6) and control[3:] == (1, 13))
    for k, v in bad.items():
        if v:
            print(f'  {k}: {v[:5]}')
    done()


if __name__ == '__main__':
    main()
