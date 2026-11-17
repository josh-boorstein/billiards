"""[Deposit] Helper, imported by the scripts of this directory. Paper: "Perpendicular billiards in rational polygons of odd denominator". Research name: probes/s928_poly_thmP.py, copied verbatim except imports (deposit names), paths (data/ and figs/ under odd_polygons/) and this paragraph. Run through run_all.py; the usage lines below are the research repo's. Bracketed IDs and section marks (§sNNN) refer to the project's internal ledgers.

s928 keeper (queue `s927-d`): Theorem P [TRI-THM-P] in rational POLYGONS — the area test of (M⁺).

PRE-REGISTERED (s928, before any row was run).
  Claim (hand port, `general_triangles.md` §s928): at odd N, in every rational k-gon (convex or not),
  no signed relation among the swapped-pair lengths ⟹ every regular orbit of the perpendicular class
  meets a side perpendicularly (M⁺); in particular p <= 1 ⟹ (M⁺).
  Instrument (cylinder-free): each invariant cylinder crosses Fix(ι) in exactly two fixed transversals,
  which are billiard BRANCHES (components of a vertical side minus its singular points); a branch b of
  width h_b whose perpendicular launch first collides perpendicularly after length λ_b lies on a closed
  leaf of length 2λ_b. Hence  area(K) = c · Σ_b h_b λ_b  with c = 1 at odd N (one vertical arc per side
  type on X = S, 2N copies) and c = 2 at even N (on S: two vertical arcs per side type).
  Ratio R := c Σ_b h_b λ_b / (2N · area(Π)) = area(K)/area(S); (M⁺) on S ⟺ R = 1.
  Predictions: (i) R = 1 (to 1e-25) on every odd-N row whose swapped-pair lengths have no signed
  relation — including the p = 1 truncations of [NCYL-709]; (ii) λ is constant on each branch
  (checked at 1/3 and 2/3); (iii) CONTROL, must fire: even-N rows with p >= 1 have R < 1 on S (the
  §s927 (2) even control: free cylinders exactly on the p >= 1 rows); even p = 0 rows give R = 1
  (normalisation check for c = 2).
  Could fail: a free cylinder or minimal component at odd N (R < 1); a wrong normalisation shows on
  the validation triangles (odd p = 0 rows, R must be 1) before any polygon is read.
  PRIOR ART: §s910 (2) ([POLY-CP-CERT], the port of Theorem C), [NCYL-709] (odd-N polygons with p = 1),
  §s927 (2)-(3), [NCYL-719]; `rulings.py --grep polygon` → [LIT-012], [NCYL-709]; no area test of
  invariant cylinders in the ledgers (grep "area(K)", "Σ h λ" → none).
  DATA: `data/s910_poly_sweep.jsonl` (p, s′ per class; no branch widths or lengths) — not reusable for R.
Output: data/s928_poly_thmP.json.  Run from the repo root with the project PYTHONPATH (~2.5 min).
CASES lists 14 odd-N polygons; 3 with an angle a_i = N ((2,3,4,9)/9 s1, s2; (3,4,5,6,9)/9 s1) are skipped by Poly's
singular solve (stdout only), leaving 11 (audit C9).
"""
import _setup  # noqa: F401  (deposit: this directory on the path, cwd = odd_polygons/)
import itertools
import json
import sys

import mpmath as mp

from poly_slots import Poly

mp.mp.dps = 60


def dist(A, B):
    return mp.sqrt((A[0] - B[0]) ** 2 + (A[1] - B[1]) ** 2)


def area(poly):
    V = poly.V
    return abs(sum(V[i][0] * V[(i + 1) % poly.k][1] - V[(i + 1) % poly.k][0] * V[i][1]
                   for i in range(poly.k))) / 2


def walk(poly, P, d, cap):
    """from P in direction d: ('perp'|'vertex'|'cap', length, last side)"""
    N, lam = poly.N, mp.mpf(0)
    for _ in range(cap):
        h = poly.first_hit(P, d, None)
        lam += dist(P, h[2])
        if h[0] == 'vertex':
            return ('vertex', lam, h[1])
        _, j, X, s = h
        if (d - poly.m[j] - N) % (2 * N) == 0:
            return ('perp', lam, j)
        d = (2 * poly.m[j] - d) % (4 * N)
        P = X
    return ('cap', lam, None)


def side_point(poly, j, s):
    A, B = poly.V[j], poly.V[(j + 1) % poly.k]
    return (A[0] + s * (B[0] - A[0]), A[1] + s * (B[1] - A[1]))


def row(label, poly, cap=20000):
    N, k = poly.N, poly.k
    out = []
    for C in poly.classes():
        vert = [j for j in range(k) if (poly.idx[j] + N) % (4 * N) in C]
        res = {(i, d): poly.prong(i, d, cap) for (i, d) in poly.slots(C)}
        perp = [(r[1], r[2]) for r in res.values() if r[0] == 'perp']
        sprime = sum(1 for r in res.values() if r[0] == 'cap')
        pairs = []
        for (i, d), r in res.items():
            if r[0] != 'vertex':
                continue
            rev = (r[1], (r[2] + 2 * N) % (4 * N))
            if rev != (i, d) and res.get(rev, ('',))[0] == 'vertex' and (i, d) < rev:
                w = walk(poly, poly.V[i], d, cap)
                pairs.append(w[1])
        p = len(pairs)
        # signed relations among pair lengths
        rel = None
        if p >= 2:
            best = None
            for eps in itertools.product((-1, 0, 1), repeat=p):
                if any(eps) and next(e for e in eps if e) == 1:
                    v = abs(sum(e * l for e, l in zip(eps, pairs)))
                    best = v if best is None or v < best else best
            rel = float(best)
        S, lam_bad, capped = mp.mpf(0), 0, 0
        for j in vert:
            Lj = dist(poly.V[j], poly.V[(j + 1) % k])
            pts = sorted([mp.mpf(0)] + [s for (jj, s) in perp if jj == j] + [mp.mpf(1)])
            d0 = (poly.idx[j] + N) % (4 * N)
            for a, b in zip(pts, pts[1:]):
                if b - a < mp.mpf(10) ** -25:
                    continue
                ws = [walk(poly, side_point(poly, j, a + (b - a) * t), d0, cap) for t in (mp.mpf(1) / 3, mp.mpf(2) / 3)]
                if any(w[0] != 'perp' for w in ws):
                    capped += 1
                    continue
                if abs(ws[0][1] - ws[1][1]) > mp.mpf(10) ** -25:
                    lam_bad += 1
                S += (b - a) * Lj * ws[0][1]
        c = 1 if N % 2 else 2
        R = c * S / (2 * N * area(poly)) if capped == 0 else mp.nan   # a skipped capped branch would fake R < 1 (audit s928b G 1-13)
        out.append({'label': label, 'a': poly.a, 'N': N, 'k': k, 'v': len(vert), 'p': p, 'sprime': sprime,
                    'pair_lengths': [float(x) for x in pairs], 'min_signed_rel': rel,
                    'R': mp.nstr(R, 30), 'one_minus_R': float(1 - R), 'lam_nonconst': lam_bad,
                    'branch_capped': capped})
    return out


def truncations():
    """[NCYL-709]: (2,5,8)/15 with V1 cut at distance t -> (2,10,10,8)/15, p = 1"""
    T = Poly([2, 5, 8], 15, seed=0)
    V0, V1, V2 = T.V

    def pt(A, B, t):
        d = dist(A, B)
        return (B[0] + (A[0] - B[0]) * t / d, B[1] + (A[1] - B[1]) * t / d)
    for t in ('0.01', '0.05', '0.15', '0.3'):
        t = mp.mpf(t)
        P, Pp = pt(V0, V1, t), pt(V2, V1, t)
        L = [dist(V0, P), dist(P, Pp), dist(Pp, V2), dist(V2, V0)]
        yield f'trunc (2,5,8)/15 t={mp.nstr(t, 3)}', Poly([2, 10, 10, 8], 15, lengths=L)


CASES = [
    # validation: odd-N triangles (R must be 1; (2,5,8)/15 has p = 1)
    ('tri', [2, 3, 4], 9, 0), ('tri', [2, 5, 8], 15, 0), ('tri', [1, 3, 5], 9, 0), ('tri', [3, 5, 7], 15, 0),
    # odd-N polygons (s910 cases ⊕ more seeds)
    ('quad', [2, 3, 4, 9], 9, 1), ('quad', [2, 3, 4, 9], 9, 2), ('quad', [3, 3, 4, 4], 7, 1),
    ('quad', [1, 2, 5, 6], 7, 3), ('quad', [2, 4, 5, 7], 9, 4), ('quad', [3, 4, 5, 6], 9, 5),
    ('quad', [3, 4, 5, 6], 9, 11), ('quad', [4, 5, 9, 12], 15, 1), ('quad', [2, 7, 8, 13], 15, 2),
    ('pent', [3, 4, 5, 6, 9], 9, 1), ('pent', [4, 5, 5, 6, 7], 9, 2), ('pent', [5, 6, 8, 12, 14], 15, 3),
    ('nonconvex', [1, 1, 11, 5], 9, 8), ('nonconvex', [2, 1, 13, 2], 9, 9),
    # CONTROL: even N (R < 1 expected exactly where p >= 1)
    ('tri', [1, 2, 3], 6, 0), ('tri', [3, 4, 5], 12, 0), ('tri', [2, 3, 5], 10, 0), ('tri', [3, 5, 8], 16, 0),
    ('tri', [5, 6, 7], 18, 0), ('tri', [3, 7, 10], 20, 0),
    ('quad', [1, 3, 5, 7], 8, 1), ('quad', [2, 3, 5, 6], 8, 2), ('quad', [2, 2, 3, 5], 6, 6),
    ('quad', [1, 2, 3, 10], 8, 1), ('quad', [2, 3, 5, 6], 8, 5),
]


def cut(a, N, i, x, t, seed=0):
    """truncate vertex i of the triangle a/N at distance t: angles (x, y), x + y = N + a_i"""
    a3 = [a[(i - 1) % 3], a[i], a[(i + 1) % 3]]
    T = Poly(a3, N, seed=seed)
    V0, V1, V2 = T.V
    y = N + a3[1] - x
    t = mp.mpf(t)
    sx, sy, sa = mp.sin(mp.pi * x / N), mp.sin(mp.pi * y / N), mp.sin(mp.pi * a3[1] / N)
    L = [dist(V0, V1) - t, t * sa / sy, dist(V1, V2) - t * sx / sy, dist(V2, V0)]
    return Poly([a3[0], x, y, a3[2]], N, lengths=L)


# PART 2 (pre-registered with part 1): truncations of the stored odd-N p >= 2 triangles (§s927 (4)), every
# vertex, two cut splits, two depths; rows kept whatever p comes out (a truncation may destroy pairs).
TRUNC2 = [([3, 13, 29], 45), ([3, 14, 28], 45), ([5, 21, 49], 75), ([3, 19, 41], 63)]


def main():
    rows = []
    for kind, a, N, seed in CASES:
        lab = f'{kind} {tuple(a)}/{N} s{seed}'
        try:
            P = Poly(a, N, seed=seed)
        except (ValueError, AssertionError, ZeroDivisionError) as e:
            print(f'{lab}: skipped ({e})')
            continue
        rows += row(lab, P)
    for lab, P in truncations():
        rows += row(lab, P)
    for a, N in TRUNC2:
        for i in range(3):
            for x in sorted({(N + a[i] + 1) // 2, a[i] + 1 + (N - a[i]) // 4}):
                if not (a[i] < x < N and a[i] < N + a[i] - x < N):
                    continue
                for t in ('0.02', '0.1'):
                    lab = f'cut {tuple(a)}/{N} v{i} x={x} t={t}'
                    try:
                        P = cut(a, N, i, x, t)
                    except (ValueError, AssertionError, ZeroDivisionError) as e:
                        print(f'{lab}: skipped ({e})')
                        continue
                    rows += row(lab, P)
    json.dump(rows, open('data/s928_poly_thmP.json', 'w'), indent=1)
    for r in rows:
        print(f"{r['label']:34s} cls v={r['v']} p={r['p']} s'={r['sprime']} 1-R={r['one_minus_R']:+.2e} "
              f"rel={r['min_signed_rel'] if r['min_signed_rel'] is None else round(r['min_signed_rel'], 4)} "
              f"lamNC={r['lam_nonconst']} cap={r['branch_capped']}")


if __name__ == '__main__':
    main()
