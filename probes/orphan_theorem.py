#!/usr/bin/env python3
"""
orphan_theorem.py — verifier for the Orphan / Degeneracy Theorem.

Claim: for alpha=(P/Q)(pi/2), gcd(P,Q)=1, 0<P<Q, the perpendicular partition of
(0,1) has, on a COMPLETE partition, exactly ONE self-paired ("orphan") branch if
P is odd and ZERO if P is even. The orphan is the DEGENERATE branch on which the
first-return map T is the identity (retracing orbits: they hit a wall head-on and
return to their own launch height).

Two independent checks per (P,Q):
  A. exact width-pairing parity (cyclotomic keys) on the certified-complete tree;
  B. the geometric signature: sample the orphan and confirm T(s)=s + an interior
     head-on wall hit, and a paired branch has neither.
"""
import math
from collections import Counter

import numpy as np

from cos_poly_fast import CosPoly
from n_exact import _canonical, _reduction_table, exact_key
from right_triangle_billiards import RightTriangleBilliard
from word_tree_nofrac import build_partition_tree


def branch_data(P, Q, depth):
    """Return (branch_edges, orphan_count, complete) using exact width-keys."""
    g = math.gcd(P, Q)
    P, Q = P // g, Q // g
    alpha = (P / Q) * (math.pi / 2)
    N = 2 * Q
    _, _, bp, trunc = build_partition_tree(
        alpha, max_depth=depth, closure_tol=1e-6, min_width=1e-11,
        return_raw=True)
    red, D = _reduction_table(N)
    zk = _canonical([0] * D, 0)
    ok = exact_key(np.array([1]), 0, P, N, red, D)
    seen = {}
    for poly, _vt in bp:
        k = exact_key(poly.coeffs, poly.shift, P, N, red, D)
        if k == zk or k == ok or k in seen:
            continue
        seen[k] = (poly.evaluate(alpha), poly)
    items = sorted(seen.values())
    polys = ([CosPoly(np.array([0]), 0)] + [p for _, p in items] +
             [CosPoly(np.array([1]), 0)])
    wk = [exact_key((polys[i + 1] - polys[i]).coeffs,
                    (polys[i + 1] - polys[i]).shift, P, N, red, D)
          for i in range(len(polys) - 1)]
    c = Counter(wk)
    orphan_keys = {k for k in c if c[k] % 2 == 1}
    edges = [0.0] + sorted(v for v, _ in items) + [1.0]
    orphan_idx = [i for i in range(len(edges) - 1)
                  if wk[i] in orphan_keys]
    return edges, orphan_idx, trunc == 0


def certified(P, Q, cap=40000):
    d = 16 * Q
    while d <= cap:
        edges, orph, comp = branch_data(P, Q, d)
        if comp:
            return edges, orph, True, d
        d *= 2
    return edges, orph, False, d // 2


def geometric_check(P, Q, lo, hi):
    """T(mid)=mid?  interior head-on wall hit?  (retrace signature)."""
    alpha = (P / Q) * (math.pi / 2)
    bil = RightTriangleBilliard(alpha)
    s = 0.5 * (lo + hi)
    tr = bil.trace(s, max_hits=200000, stop_at_perpendicular_return=True)
    Ts = tr.l1_returns[-1].s if (tr.returned_perpendicular and
                                 tr.l1_returns) else None
    nrm = {'L1': (1, 0), 'L2': (0, 1),
           'H': (math.sin(alpha), -math.cos(alpha))}
    headon = 1.0
    for h in tr.hits[:-1]:                       # exclude final L1 return
        vx, vy = math.cos(h.incoming_angle), math.sin(h.incoming_angle)
        nx, ny = nrm[h.side]
        headon = min(headon, abs(vx * ny - vy * nx))
    return Ts, s, headon


if __name__ == '__main__':
    ODD = [(3, 8), (5, 9), (3, 11), (5, 11), (7, 11), (3, 13), (5, 13),
           (7, 13), (9, 13), (11, 13), (3, 17), (5, 17), (3, 19), (5, 23),
           (9, 23), (7, 25), (11, 25), (13, 25), (9, 29), (11, 31)]
    EVEN = [(2, 7), (4, 9), (2, 11), (4, 11), (2, 13), (6, 13), (8, 13)]

    print("A. EXACT width-pairing parity on certified-complete partitions")
    print("   odd P  -> expect 1 orphan;  even P -> expect 0")
    okA = True
    for P, Q in ODD:
        edges, orph, comp, d = certified(P, Q)
        exp = 1
        good = comp and len(orph) == exp
        okA &= good
        print(f"   {P}/{Q}: n={len(edges)-1} orphans={len(orph)} "
              f"complete={comp} depth={d} {'OK' if good else 'CHECK'}")
    for P, Q in EVEN:
        edges, orph, comp, d = certified(P, Q, cap=80000)
        good = comp and len(orph) == 0
        okA &= good
        print(f"   {P}/{Q}: n={len(edges)-1} orphans={len(orph)} "
              f"complete={comp} depth={d} {'OK' if good else 'CHECK'}")
    print(f"   => A passed: {okA}\n")

    print("B. Geometric retrace signature (odd P): orphan T(s)=s + interior "
          "head-on; a paired branch has neither")
    for P, Q in [(3, 8), (5, 9), (3, 11), (7, 13)]:
        edges, orph, comp, d = certified(P, Q)
        i = orph[0]
        Ts, s, headon = geometric_check(P, Q, edges[i], edges[i + 1])
        # a paired (non-orphan) branch
        j = next(k for k in range(len(edges) - 1) if k not in orph)
        Ts2, s2, headon2 = geometric_check(P, Q, edges[j], edges[j + 1])
        print(f"   {P}/{Q}: ORPHAN s={s:.4f} T(s)-s={Ts-s:+.1e} "
              f"interior_headon={headon:.1e}  |  PAIRED s={s2:.4f} "
              f"T(s)-s={Ts2-s2:+.2e} interior_headon={headon2:.2e}")
