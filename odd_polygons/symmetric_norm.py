"""[Deposit] Helper, imported by the scripts of this directory. Paper: "Perpendicular billiards in rational polygons of odd denominator". Research name: probes/s625_symmetric_norm.py, copied verbatim except imports (deposit names), paths (data/ and figs/ under odd_polygons/) and this paragraph. Run through run_all.py; the usage lines below are the research repo's. Bracketed IDs and section marks (§sNNN) refer to the project's internal ledgers.

s625 -- [NCYL-633] DO-NOT (c): is each known minimal torus's loop relation the SHORTEST relation of
its triangle in an S3-SYMMETRIC norm, not only in the basis-dependent Euclidean (K, L)?

PRE-REGISTERED.  A relation K alpha + L beta in pi Z is (K, L, 0) in Z^3 / (1, 1, 1) (gamma = pi -
alpha - beta); permuting the angles permutes the coordinates.  Norms: the hexagonal form
H = K^2 - K L + L^2 (the projection perpendicular to (1,1,1), up to 2/3) and the coset norms
min_t |(K+t, L+t, t)|_p, p = 1, 2, inf.  The loop line's shortest LATTICE vector (not its primitive
integer vector: (7,8,15)'s loops are (2,2)) against the lattice minimum; unique up to sign, or TIE.
COULD COME OUT OTHERWISE: a strictly shorter non-parallel relation in some symmetric norm.
PRIOR ART: [NCYL-633] (3) (Euclidean only; DO-NOT (c)).  DATA: data/s624_formal_rotation.json
(the loops' (K, L)); (13,26,119)'s K = -2L from [NCYL-633] ↓s624.
RESULT (s625): the ten ⊕ (13,26,119): OK in all five norms, unique.  (5,7,11) (the side-parallel
class, not one of the ten): its (1,-4) loses to 5 alpha + 3 beta = 2 pi in H (21 vs 19) and L2 (14 vs
13), ties in L1 / Linf.  -> [NCYL-633] ↓s625.
"""
import _setup  # noqa: F401  (deposit: this directory on the path, cwd = odd_polygons/)
import json
from math import gcd


def coset(v, p):
    K, L = v
    best = None
    for t in range(-abs(K) - abs(L) - 1, abs(K) + abs(L) + 2):
        w = (K + t, L + t, t)
        n = max(map(abs, w)) if p == 'inf' else sum(abs(x) ** p for x in w)
        best = n if best is None else min(best, n)
    return best


NORMS = {'E(K,L)': lambda v: v[0] ** 2 + v[1] ** 2, 'H': lambda v: v[0] ** 2 - v[0] * v[1] + v[1] ** 2,
         'L1': lambda v: coset(v, 1), 'L2': lambda v: coset(v, 2), 'Linf': lambda v: coset(v, 'inf')}


def main():
    dirs = {}
    for r in json.load(open('data/s624_formal_rotation.json')):
        for comp in r['min_K']:
            for K, L in comp:
                if (K, L) == (0, 0):
                    continue
                g = gcd(K, L)
                k, l = K // g, L // g
                if k < 0 or (k == 0 and l < 0):
                    k, l = -k, -l
                dirs.setdefault(tuple(r['abc']), set()).add((k, l))
    dirs[(13, 26, 119)] = {(2, -1)}
    rows = []
    for abc, ds in sorted(dirs.items()):
        a, b, c = abc
        d = a + b + c
        lat = [(K, L) for K in range(-d, d + 1) for L in range(-d, d + 1)
               if (K, L) != (0, 0) and (K * a + L * b) % d == 0]
        (k0, l0), = ds
        s = next(s for s in range(1, d + 1) if (s * k0 * a + s * l0 * b) % d == 0)
        k0, l0 = s * k0, s * l0
        row = {'abc': abc, 'loop': (k0, l0)}
        for nm, f in NORMS.items():
            m = min(f(v) for v in lat)
            tied = sorted({v if v > (0, 0) else (-v[0], -v[1]) for v in lat if f(v) == m})
            nonpar = [v for v in tied if v[0] * l0 - v[1] * k0 != 0]
            row[nm] = ('OK' if f((k0, l0)) == m and not nonpar else 'TIE' if f((k0, l0)) == m else 'FAIL',
                       f((k0, l0)), m, nonpar[:2])
        rows.append(row)
        print(abc, d, (k0, l0), ' '.join(f"{nm}:{row[nm][0]}({row[nm][1]} vs {row[nm][2]})" for nm in NORMS))
    json.dump(rows, open('data/s625_symmetric_norm.json', 'w'), indent=1, default=str)


if __name__ == '__main__':
    main()
