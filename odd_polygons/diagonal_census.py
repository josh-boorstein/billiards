"""[Deposit] Helper, imported by the scripts of this directory. Paper: "Perpendicular billiards in rational polygons of odd denominator". Research name: probes/s933_diagonal_census.py, copied verbatim except imports (deposit names), paths (data/ and figs/ under odd_polygons/) and this paragraph. Run through run_all.py; the usage lines below are the research repo's. Bracketed IDs and section marks (§sNNN) refer to the project's internal ledgers.

s933 keeper (queue `s929-e`, re-pointed): is every odd-N non-CP polygon class a MIRROR or a FOLD? — an exact census of the
convex quadrilaterals cut by a rational diagonal into two rational triangles.

PRE-REGISTERED (s933, before the census; a 2-shape timing run was seen, both CP).
HYPOTHESIS (H) (coined s933, §s933 (1)): every non-CP perpendicular class of an odd-N rational polygon comes from a mirror symmetry
(kites, [POLY-KITE-NONCP]) or a reflection fold onto an even-N non-CP triangle ([POLY-FOLD-ODD]). Refuted by ONE non-CP
quadrilateral with no mirror and no fold.
OBJECTS: odd N, convex Π = T₁ ∪ T₂ glued along the diagonal V0V2, T₁ = (x, 2a1, x′) and T₂ = (y′, 2a3, y) in units π/(2N)
(T₁ on V0 V1 V2, T₂ on V2 V3 V0; parents of even N_e = 2N allowed), Π = ((x+y)/2, a1, (x′+y′)/2, a3)/N, every angle < N,
gcd(angles, N) = 1. Exact sides: T₁ scaled by sin(2a3·π/2N), T₂ by sin(2a1·π/2N) (law of sines; products of sines lie in ℤ[ζ_{4N}]/4).
These are the countable 'special' shapes where length coincidences live (every kite is one: T₂ = T₁'s mirror).
INSTRUMENT: `engine/flow_poly.py` — exact CP / non-CP of Π's single perpendicular class (one decomposition per shape).
PREDICTION: no non-CP shape among those with no mirror symmetry, at every N ≤ NMAX.
COULD HAVE COME OUT OTHERWISE: a `no_periodic` component at a shape with no mirror (then check for a fold by hand before
claiming anything).  POSITIVE CONTROL: the kite (41,15,67,15)/69 built by the same constructor (T₁ = T₂′ = (41,30,67) half-units)
must give `no_periodic` (4 of 70, §s929 (2)).  Mirror shapes (kites: T₂ = mirror of T₁; and any quadrilateral whose
(angle, side) cycle is a reversal of itself) are flagged and reported separately.
PRIOR ART: §s930 (2) ⊕ [NCYL-721] (16 float slices, no asymmetric non-CP); §s929 (2) (the only odd-N non-CP known before s933 —
kites at N = 69, 127, 135); §s933 (1) (the fold octagons, non-convex); `rulings.py --grep` *diagonal* → other objects, *non-CP
polygon* → [NCYL-721]. → no census of rational-diagonal quadrilaterals on record.
DATA: `data_index.py --grep quadrilateral` → none.
Output: data/s933_diagonal_census.json.  Usage: python probes/s933_diagonal_census.py NMAX
"""
import _setup  # noqa: F401  (deposit: this directory on the path, cwd = odd_polygons/)
import json
import math
import sys
import time

import flow_poly as FP
from flow_decomp import mmul


def quad(N, x, a1, xp, y, a3, yp):
    """x, xp, y, yp in units π/(2N); a1, a3 in units π/N"""
    d, n4 = N, 4 * N
    S = lambda t: FP.sin_half(t, d)
    L = [mmul(S(xp), S(2 * a3), n4), mmul(S(x), S(2 * a3), n4), mmul(S(y), S(2 * a1), n4), mmul(S(yp), S(2 * a1), n4)]
    return FP.PolyT([(x + y) // 2, a1, (xp + yp) // 2, a3], d, L)


def mirror_free(x, a1, xp, y, a3, yp):
    """Π's (angle, side) cycle, sides by exact law-of-sines factors; a mirror ⟺ the cycle equals a reversal of itself"""
    ang = [(x + y), 2 * a1, (xp + yp), 2 * a3]
    # side lengths (float, 12 digits): sin(xp)sin(2a3), sin(x)sin(2a3), sin(y)sin(2a1), sin(yp)sin(2a1)  (half-units)
    s = lambda t: math.sin(math.pi * t / (2 * N_CUR))
    L = [round(s(xp) * s(2 * a3), 10), round(s(x) * s(2 * a3), 10), round(s(y) * s(2 * a1), 10), round(s(yp) * s(2 * a1), 10)]
    n = 4
    for r in range(n):
        if all(ang[(r - i) % n] == ang[i] and L[(r - i - 1) % n] == L[i] for i in range(n)):
            return False
    return True


N_CUR = 0


def census(N):
    global N_CUR
    N_CUR = N
    rows, seen = [], set()
    for a1 in range(1, N):
        for a3 in range(1, N):
            for x in range(1, 2 * N - 2 * a1):
                xp = 2 * N - 2 * a1 - x
                for y in range(1, 2 * N - 2 * a3):
                    yp = 2 * N - 2 * a3 - y
                    if (x + y) % 2 or (xp + yp) % 2 or x + y >= 2 * N or xp + yp >= 2 * N:
                        continue
                    ang = [(x + y) // 2, a1, (xp + yp) // 2, a3]
                    if math.gcd(math.gcd(*ang), N) != 1:
                        continue
                    # dedupe the dihedral relabellings of the same (Π, diagonal): rotate by 2, reverse
                    key = min((x, a1, xp, y, a3, yp), (y, a3, yp, x, a1, xp), (xp, a1, x, yp, a3, y), (yp, a3, y, xp, a1, x))
                    if key in seen:
                        continue
                    seen.add(key)
                    P = quad(N, x, a1, xp, y, a3, yp)
                    r = FP.decompose_poly(P, FP.perp_dirs_poly(P)[0], cap=20000)
                    bad = r['counts'].get('no_periodic', 0) + r['counts'].get('undetermined', 0)
                    rows.append({'key': key, 'angles': ang, 'mirror_free': mirror_free(x, a1, xp, y, a3, yp),
                                 'counts': r['counts'], 'bad': bad})
    return rows


def main():
    NMAX = int(sys.argv[1])
    out = {}
    t = time.time()
    ctl = FP.decompose_poly(quad(69, 41, 15, 67, 41, 15, 67), FP.perp_dirs_poly(quad(69, 41, 15, 67, 41, 15, 67))[0])
    assert ctl['counts'].get('no_periodic') == 4, ctl['counts']
    out['control'] = {'shape': '(41,15,67,15)/69', 'counts': ctl['counts']}
    print('control fires:', ctl['counts'], flush=True)
    for N in range(5, NMAX + 1, 2):
        rows = census(N)
        nb = [r for r in rows if r['bad']]
        nbf = [r for r in nb if r['mirror_free']]
        out[str(N)] = {'shapes': len(rows), 'mirror_free': sum(r['mirror_free'] for r in rows), 'nonCP': len(nb),
                       'nonCP_mirror_free': len(nbf), 'nonCP_rows': nb, 'secs': round(time.time() - t)}
        print(f"N={N}: {len(rows)} shapes ({out[str(N)]['mirror_free']} mirror-free), non-CP {len(nb)}, "
              f"non-CP mirror-free {len(nbf)}  [{time.time() - t:.0f}s]", flush=True)
        json.dump(out, open('data/s933_diagonal_census.json', 'w'), indent=1)


if __name__ == '__main__':
    main()
