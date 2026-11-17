"""[Deposit] Theorem 5.4: the identities transferred to real x; Lemma 5.5 inputs. Paper: "Perpendicular billiards in rational polygons of odd denominator". Research name: probes/s952_s1_family.py, copied verbatim except imports (deposit names), paths (data/ and figs/ under odd_polygons/) and this paragraph. Run through run_all.py; the usage lines below are the research repo's. Bracketed IDs and section marks (§sNNN) refer to the project's internal ledgers.

s952 keeper (queue `s952-b`): S₁ as a family in the angle — as [POLY-PI-X] did for Π_N. S₁(m)'s polygon at x is
S₁(m; x) = (π−(3m+1)x, 4mx, π−(m+1)x, 2x), sides ∝ (2 sinx cos mx, 2 sinx cos 3mx, sin 4mx, sin 4mx), isosceles apex; at
x = aπ/N it is (N−3ma−a, 4ma, N−ma−a, 2a)/N (census quad(N, N−6ma, 4ma, N−2ma, N−2a, 2a, N−2a)).

WHY IT SHOULD TRANSFER (checked in stage B): s951's identities are rational identities in INDEPENDENT symbols z, w
(`s951_series_first_look.formal`, sympy) — true at every z = e^{ix/2}, w = e^{imx}; REF = DIR² monomial-wise; the copy rule is
reflection composition; so the walk at any real x is a chain of reflected copies of S₁(m; x) (the [POLY-PI-X] argument). The
inequality certificate `probes/s951_s1_allm_cert.py` is stated in (θ, ξ) = (mx, 1/m) for m ≥ 4 — region θ(14 + ξ) ≤ π, a
continuum in θ, no N — and per m on x ∈ (0, X], X ≥ π/(14m + 1), for m = 2, 3; m = 1 is [POLY-PI-X]. So it covers every real
x ≤ π/(14m + 1). Start copy (3m−1)/2: billiard direction d = i·e^{i(3m−2)x}; at x = aπ/N, N odd, d ∈ ±i·μ_N (perpendicular).

STAGE A — PRE-REGISTERED (s952, before any row ran): for (m, a) ∈ {(2,2), (2,3), (3,2), (4,3)}:
  (A1) at the first two odd N > a(14m+1) with gcd(a, N) = 1: S₁(m; aπ/N) is non-CP in its perpendicular class, with a non-periodic
       component at rotation number 4 sin²(maπ/N);
  (A2) at the largest such N < a(14m+1): CP (the a = 1 onset law, §s946 (3), extended);
  (A3) the census polygon IS S₁(m; x) (angles exact, side ratios to 1e-12);
  (A4) no-fold, MEASURED: s946's local brute force (`data/scratch/s946/audit2/lf.py`, M = 300) leaves 0 survivors on each A1 row
       (Lemma F′'s hand case analysis was done for a = 1 only — not claimed for a ≥ 2).
  ⚠ OUTCOME (s952): A1, A3, A4 hold; A2 FAILED — at the largest admissible N < a(14m+1) every row is NON-CP with the torus and 0
  containment violations, and the edge control x = π/(14m+1)·(1 + 10⁻³) shows 0. My error: A2 extrapolated s951's CERTIFIED region
  θ(14 + ξ) ≤ π, not the template's true edge (at a = 1 no odd N lies between them). AMENDMENT, PRE-REGISTERED before running (H′):
  the template's edge is θ = mx = π/14 for every m (s951's corner; Π(x)'s edge at m = 1): containment 0 at θ = (π/14)(1 − 10⁻³),
  > 0 at (π/14)(1 + 10⁻³), m = 2, 3, 4, 6, 10; and at the largest admissible N < 14ma the class is CP (`h_prime()`).
STAGE B: s951's formal identities re-run (independent z, w); s951's certificate verdict read; containment numerically at every A1
  row's x (0 violations) and at x = π/(14m + 1)·(1 + 10⁻³) (must be > 0 — the region's edge is sharp).
PRIOR ART: §s952 (1) ([POLY-PI-X], the transfer argument), §s951 (2)–(4) ([POLY-S1-ALL], its certificate), §s946 (3) (onset law,
  Lemma F′); queue s952-b; `rulings.py --grep` series → [NCYL-736] (S₂ only). DATA: `data/s951_s1_allm_cert.json`, `data/s946_apex_census.json`.
→ general_triangles.md §s952 (4).   Usage: python probes/s952_s1_family.py   → data/s952_s1_family.json
"""
import _setup  # noqa: F401  (deposit: this directory on the path, cwd = odd_polygons/)
import json
import math
import sys
import time

import flow_poly as FP                                  # noqa: E402
import thmC_identities as S1                     # noqa: E402
from thmC_m1 import components                   # noqa: E402
from diagonal_census import quad                   # noqa: E402
from fold_survivors import survivors              # noqa: E402  (s958: was data/scratch/s946/audit2/lf.py)

OUT = 'data/s952_s1_family.json'
CASES = [(2, 2), (2, 3), (3, 2), (4, 3)]


def poly(N, m, a):
    return quad(N, N - 6 * m * a, 4 * m * a, N - 2 * m * a, N - 2 * a, 2 * a, N - 2 * a)


def shape_err(P, N, m, a):
    import cmath
    x = a * math.pi / N
    L = [abs(sum(c * cmath.exp(2j * math.pi * e / P.N) for e, c in l.items())) for l in P.L]
    ref = [2 * math.sin(x) * math.cos(m * x), 2 * math.sin(x) * math.cos(3 * m * x), math.sin(4 * m * x), math.sin(4 * m * x)]
    return max(abs(L[i] / L[0] / (ref[i] / ref[0]) - 1) for i in range(4))


def row(N, m, a, fold=False):
    t = time.time()
    P = poly(N, m, a)
    ang = [N - 3 * m * a - a, 4 * m * a, N - m * a - a, 2 * a]
    assert list(P.a) == ang, (list(P.a), ang)
    se = shape_err(P, N, m, a)
    j = FP.perp_dirs_poly(P)[0]
    kinds, rhos = components(P, j)
    target = 4 * math.sin(m * a * math.pi / N) ** 2
    r = {'N': N, 'm': m, 'a': a, 'angles': ang, 'shape_err': se, 'kinds': kinds, 'rhos': rhos, 'target': target,
         'nonCP': bool(kinds.get('no_periodic') or kinds.get('undetermined')),
         'torus_seen': any(x is not None and abs(x - target) < 1e-9 for x in rhos),
         'containment_viol': S1.check(m, N / a)}
    if fold:
        r['fold_survivors'] = len(survivors(ang, N, M=300))
    r['secs'] = round(time.time() - t, 1)
    print({k: v for k, v in r.items() if k != 'rhos'}, flush=True)
    return r


def main():
    out = {'A1': [], 'A2': []}
    for m, a in CASES:
        on = a * (14 * m + 1)
        above = [N for N in range(on + 1, on + 60) if N % 2 and math.gcd(a, N) == 1][:2]
        below = [N for N in range(on - 1, 0, -1) if N % 2 and math.gcd(a, N) == 1][0]
        out['A1'] += [row(N, m, a, fold=True) for N in above]
        out['A2'].append(row(below, m, a))
    out['B'] = {'formal': S1.formal(),
                'cert_verdict': json.load(open('data/s951_s1_allm_cert.json')).get('verdict'),
                'edge_control_viol': {m: S1.check(m, (14 * m + 1) / (1 + 1e-3)) for m in (2, 3, 4, 6)}}
    out['verdict'] = {
        'A1': all(r['nonCP'] and r['torus_seen'] for r in out['A1']),
        'A2_CP': [(r['N'], r['m'], r['a'], not r['nonCP']) for r in out['A2']],
        'A3': max(r['shape_err'] for r in out['A1'] + out['A2']) < 1e-12,
        'A4_fold_survivors': [r['fold_survivors'] for r in out['A1']],
        'B_containment_A1': [r['containment_viol'] for r in out['A1']],
        'B_edge_fires': all(v > 0 for v in out['B']['edge_control_viol'].values())}
    print('VERDICT', out['verdict'], out['B'], flush=True)
    json.dump(out, open(OUT, 'w'), indent=1, default=str)


def h_prime():
    res = {'edge': {m: (S1.check(m, 14 * m / (1 - 1e-3)), S1.check(m, 14 * m / (1 + 1e-3))) for m in (2, 3, 4, 6, 10)}, 'rows': []}
    for m, a in CASES:
        below = [N for N in range(14 * m * a - 1, 0, -1) if N % 2 and math.gcd(a, N) == 1][0]
        res['rows'].append(row(below, m, a))
    res['verdict'] = {'edge_0_below': all(v[0] == 0 for v in res['edge'].values()),
                      'edge_fires_above': all(v[1] > 0 for v in res['edge'].values()),
                      'CP_below_14ma': [(r['N'], r['m'], r['a'], not r['nonCP']) for r in res['rows']]}
    print('H-PRIME', res['verdict'], res['edge'], flush=True)
    out = json.load(open(OUT))
    out['H_prime'] = res
    json.dump(out, open(OUT, 'w'), indent=1, default=str)


if __name__ == '__main__':
    h_prime() if 'hprime' in sys.argv else main()
