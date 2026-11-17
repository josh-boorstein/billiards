"""[Deposit] Helper, imported by the scripts of this directory. Paper: "Perpendicular billiards in rational polygons of odd denominator". Research name: probes/s933_regular_pieces.py, copied verbatim except imports (deposit names), paths (data/ and figs/ under odd_polygons/) and this paragraph. Run through run_all.py; the usage lines below are the research repo's. Bracketed IDs and section marks (§sNNN) refer to the project's internal ledgers.

s933 keeper (queue `s929-e` re-pointed → answered): EXPLICIT ODD-N NON-CP CONVEX QUADRILATERALS — regular-polygon pieces,
inscribed quadrilaterals, and an ASYMMETRIC infinite family that is neither a mirror nor a fold.

NOT PRE-REGISTERED: exploratory. Found by `probes/s933_diagonal_census.py` (the (H) census) — its N = 9 hit (2,2,7,7)/9 turned
out to be 4 consecutive vertices of the regular 9-gon, and its mirror-free N = 15, 17 hits (11,2,13,4)/15, (13,2,15,4)/17 one
family. Everything below is a re-run of what was then measured, with the instrument `engine/flow_poly.py` (exact; a
`no_periodic` component is a Boshernitzan/Gordan certificate, i.e. a computer proof of non-CP for that row) and an
independent tracer check (`s910_poly_slots.Poly.prong`, capped prongs) on single rows.
(a) REGULAR PIECES: k consecutive vertices of the regular n-gon, n odd ≤ 33, k = 3 … n−1: angles (k−2, k−2, n−2, …, n−2)/n,
    sides 2 sin((k−1)π/n) (the diagonal) and 2 sin(π/n). (b) the k = 4 trapezoid (2,2,n−2,n−2)/n, odd n ≤ 41, ⊕ isolation:
    the same angles with top/leg ≠ 1 (rational ratios 999/1000 … 10). (c) INSCRIBED quadrilaterals: vertices at n-th roots of
    unity, arcs (c0, c1, c2, c3) (Σ = n, up to dihedral relabelling), angles n − c_{m−1} − c_m, sides 2 sin(c_m π/n), odd n ≤ NMAX.
    (d) THE ASYMMETRIC FAMILY Π_N = (N−4, 4, N−2, 2)/N with sides (1, 2c−1, 2c, 2c), c = cos(2π/N) (= the isosceles triangle with
    apex 2π/N at V3 glued along its base to the triangle (N−6, 8, N−2)/2N), odd N = 7 … 45.
CONTROLS: the kite (41,15,67,15)/69 (positive, flow_poly gate); every CP row is a negative.
COULD HAVE COME OUT OTHERWISE: everything CP (as on s930's 16 float slices) — then (H) would stand.
PRIOR ART: §s929 (2) (non-CP kites, N ≥ 69 — the only odd-N non-CP polygons before s933), §s933 (1) (fold octagons), §s930 (2)
  ⊕ [NCYL-721] (no asymmetric non-CP on 16 slices); [LIT-007]/[LIT-011] (regular polygons as VEECH surfaces for right triangles —
  a different object); `rulings.py --grep` *regular polygon* → [LIT-*] only, *nonagon*, *inscribed* → none.
DATA: `data_index.py --grep inscribed` → none.
RESULT (s933): regular pieces non-CP only for k = 4, odd n ≡ 1 (mod 4), n ≥ 9 (all k to n = 33; k = 4 to 41), isolated in
top/leg; inscribed n ≤ 41: 152/6875 non-CP — 151 mirror (isosceles trapezoids), 1 none (arcs (2,3,5,11)/21); Π_N CP for N ≤ 13,
non-CP for N = 15 … 45 with 2⌊(N−7)/8⌋ no_periodic (extended to 101 in-thread, data/s933_asym_extend.json).
→ general_triangles.md §s933 (4), [POLY-ASYM-NONCP].
Output: data/s933_regular_pieces.json.  Usage: python probes/s933_regular_pieces.py [NMAX_INSCRIBED=41]   (≈ 6 min)
"""
import _setup  # noqa: F401  (deposit: this directory on the path, cwd = odd_polygons/)
import itertools
import json
import math
import sys
import time

import flow_poly as FP
from flow_decomp import madd, mmul


def bad_of(P):
    r = FP.decompose_poly(P, FP.perp_dirs_poly(P)[0], cap=50000)
    nb = r['counts'].get('no_periodic', 0) + r['counts'].get('undetermined', 0)
    return nb, r['counts'], sorted(round(c['mass'], 5) for c in r['components'] if c['kind'] != 'cylinder')


def regular_piece(n, k):
    a = [k - 2, k - 2] + [n - 2] * (k - 2)
    return FP.PolyT(a, n, [FP.sin_mono(k - 1, n)] + [FP.sin_mono(1, n)] * (k - 1))


def iso_trap(n, p, q):
    """(2,2,n−2,n−2)/n with top/leg = p/q"""
    N4 = 4 * n
    s1 = FP.sin_mono(1, n)
    leg, top = FP.mscale(s1, q), FP.mscale(s1, p)
    base = madd(top, mmul(leg, {4 % N4: 1, (-4) % N4: 1}, N4))   # + 2cos(2π/n)·leg
    return FP.PolyT([2, 2, n - 2, n - 2], n, [base, leg, top, leg])


def inscribed(c):
    n, k = sum(c), len(c)
    a = [n - c[m - 1] - c[m] for m in range(k)]
    return a, FP.PolyT(a, n, [FP.sin_mono(x, n) for x in c])


def canon(c):
    k = len(c)
    cs = [tuple(c[(i + r) % k] for i in range(k)) for r in range(k)]
    return min(cs + [tuple(reversed(x)) for x in cs])


def sym_type(c):
    k = len(c)
    if any(all(c[(r - i) % k] == c[i] for i in range(k)) for r in range(k)):
        return 'mirror'
    if any(all(c[(i + r) % k] == c[i] for i in range(k)) for r in range(1, k)):
        return 'rotation'
    return 'none'


def asym(N):
    """Π_N exact: sides ∝ (1, 2c−1, 2c, 2c) — built as T₁ ∪ T₂ (s933_diagonal_census.quad(N, N−6, 4, N−2, N−2, 2, N−2))"""
    from diagonal_census import quad
    return quad(N, N - 6, 4, N - 2, N - 2, 2, N - 2)


def main():
    NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 41
    out = {}
    t = time.time()
    out['regular_pieces'] = {}
    for n in range(5, 34, 2):
        row = {}
        for k in range(3, n):
            P = regular_piece(n, k)
            row[k] = bad_of(P)[0]
        out['regular_pieces'][n] = row
        print('regular', n, 'non-CP k:', [k for k, b in row.items() if b], flush=True)
    out['trap4'] = {n: bad_of(regular_piece(n, 4)) for n in range(5, 42, 2)}
    out['trap4_isolation'] = {n: {f'{p}/{q}': bad_of(iso_trap(n, p, q))[1] for p, q in
                                  [(1, 1), (999, 1000), (1001, 1000), (1, 2), (2, 1), (3, 5), (7, 5), (1, 10), (10, 1)]}
                              for n in (9, 13)}
    print('trap4', {n: v[0] for n, v in out['trap4'].items()}, flush=True)
    ins = {}
    for n in range(5, NMAX + 1, 2):
        seen, rows = set(), []
        for c in itertools.product(range(1, n), repeat=3):
            c4 = n - sum(c)
            if c4 < 1:
                continue
            cc = canon(list(c) + [c4])
            if cc in seen:
                continue
            seen.add(cc)
            a = [n - cc[m - 1] - cc[m] for m in range(4)]
            if math.gcd(math.gcd(*a), n) != 1 or max(a) >= n:
                continue
            nb, counts, masses = bad_of(inscribed(list(cc))[1])
            rows.append({'arcs': cc, 'angles': a, 'sym': sym_type(list(cc)), 'bad': nb, 'counts': counts})
        bad = [r for r in rows if r['bad']]
        ins[n] = {'total': len(rows), 'nonCP': len(bad), 'by_sym': {s: sum(r['sym'] == s for r in bad) for s in
                                                                     ('mirror', 'rotation', 'none')},
                  'nonCP_rows': bad}
        print('inscribed', n, ins[n]['total'], ins[n]['nonCP'], ins[n]['by_sym'], f'{time.time() - t:.0f}s', flush=True)
    out['inscribed'] = ins
    out['asym_family'] = {N: bad_of(asym(N)) for N in range(7, 46, 2)}
    print('asym', {N: v[0] for N, v in out['asym_family'].items()}, flush=True)
    json.dump(out, open('data/s933_regular_pieces.json', 'w'), indent=1, default=str)


if __name__ == '__main__':
    main()
