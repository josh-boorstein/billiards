"""[Deposit] helper: polygon families (imported for its row reader). Paper: "Complete periodicity of perpendicular billiards in rational triangles of prime denominator". Research name: probes/s929_asym_free.py, copied verbatim except imports (deposit names), paths (data/ and figs/ under triangles/) and this paragraph. Run through run_all.py; the usage lines below are the research repo's. Bracketed IDs and section marks (§sNNN) refer to the project's internal ledgers.

s929 keeper (queue `s928-d`): do odd-N free regions need a polygon SYMMETRY? — asymmetric quadrilaterals.

PRE-REGISTERED (s929, before any row was run).
  Target under test (`s928-d`): at odd N a free region (R < 1) needs a polygon symmetry (central symmetry `τ_h`
  or a kite doubling). Refuted by ONE asymmetric polygon with R < 1.
  Hand argument predicting REFUTATION (s929, to be written up only if the rows agree): for fixed angles the
  surfaces X_t are identified by a copy-preserving Φ commuting with ι (§s928 (2)), so ι_* on H_1 is
  t-independent. A cylinder C whose holonomy's flow-transverse component is a NON-identically-zero linear form
  in the side lengths crosses the flow direction at some t*; if ι(C) = C then Per_t(γ) = −Dι Per_t(γ) for all t,
  i.e. the transverse component vanishes identically — contradiction. So at t* ι(C) ≠ C, the two are distinct
  maximal parallel cylinders, hence disjoint, so C misses Fix(ι): a FREE cylinder, with no symmetry needed.
  Locating t*: a flow-aligned free cylinder's boundary is made of saddle connections in perpendicular-class
  directions realised only at t*, so t* is a point where some PRONG's walk changes word (it passes over a vertex).
  Method: for each asymmetric odd-N angle vector, the family L = (1, t, L2(t), L3(t)); prong outcomes on a t-grid;
  each change between neighbours is bisected (mpmath, 60 digits, to 1e-40) to the vertex-hit t*; the area test
  R (`s928_poly_thmP.row`, c = 1) is run AT t* and at t* ± 1e-6 (both must give R = 1 if the free area is
  isolated, as the zero-one law (§s928 (2)) predicts off an open bad set).
  Predictions: (i) R < 1 at some t* in every asymmetric family whose prong words change (refutes `s928-d`);
  (ii) R = 1 at t* ± 1e-6 on every such row (null part, Z₁); (iii) CONTROL that fires: the grid points
  themselves (generic t) give R = 1 — and the instrument sees free area at the known rhombus (3,2,3,2)/5.
  Could fail: R = 1 at every t* (then free regions at Z₁ are not automatic and the hand argument has a gap —
  e.g. the tilting cylinder does not survive to t*); R < 1 at t* ± 1e-6 (an open bad set: (M⁺) a.e. FAILS).
  Not every t* need carry free area: a connection realised at t* can be interior to a non-free configuration.
  AMENDED s929, after the first rows (the pre-registration above is kept as written):
  (A1) the hand argument has a GAP — a cylinder near the flow direction crosses the flow's invariant cylinders C_i in
  ONE transverse sense, so its transverse holonomy is Σ m_i h_i > 0; it cannot tilt through ω while the flow
  decomposition persists. Free area at t* needs the invariant cylinders to thin out (h_i -> 0, walks blowing up), so
  the argument predicts nothing. (A2) the parallelogram free region at ρ = 1.1 is NOT a closed cylinder: free orbits
  return with shift ±0.0753 (a rotation-like region bounded by two equal-length, 31.0, connections) — so the
  coincidence test is the VERTEX-PRONG COUNT at t* (>= 3: a second pair realised at the same t*), not a second
  bisection landing on the same t*. (A3) POSITIVE CONTROL (s929, `(5,2,5,2)/7`, window ρ ∈ [1.05, 1.15], 21 grid points):
  14 t* flagged, 13 with free area (ρ = 21/20, 19/18, 53/50, 15/14, 43/40, 27/25, 21/19, 10/9, 28/25, …;
  1 − R ∈ {0.0056, 0.1202, 0.2347}); data `data/scratch/s929/para_ctrl.json`. Some t* rows read p = 0 with R < 1:
  at a bisected t* (1e-40 off the special point) the leak branches are narrower than row()'s 1e-25 skip width — an
  instrument artefact, not a Theorem P counterexample.
  PRIOR ART: §s928 (2) (zero-one law; Φ commutes with ι), §s928 (4) (kites, parallelograms at isolated ratios,
  `s928-f`), [NCYL-719], [NCYL-709]; `rulings.py --grep` *free cylinder* → [NCYL-456], [NCYL-457] (right-triangle
  cone polarity, unrelated); *codimension*, *isolated* → nothing on polygon families; general_triangles.md
  grep *codim* → §s910's *forced on whole length families, not codimension-one accidents* (even N) only.
  DATA: `data/s928_poly_thmP.json`, `data/s928_kites.json` (R at fixed shapes only; no family scans) → none reusable.
Output: data/s929_asym_free.json.  Run from the repo root with the project PYTHONPATH.
"""
import _setup  # noqa: F401  (deposit: lib/ on the path, cwd = triangles/)
import collections
import json
import sys

import mpmath as mp

from poly_slots import Poly
from poly_thmP import row

mp.mp.dps = 60


def unit(d, N):
    th = mp.pi * d / (2 * N)
    return (mp.cos(th), mp.sin(th))


def closing(a, N, t):
    idx = [0]
    for j in range(1, 4):
        idx.append(idx[-1] + 2 * (N - a[j]))
    u = [unit(x % (4 * N), N) for x in idx]
    sx, sy = u[0][0] + t * u[1][0], u[0][1] + t * u[1][1]
    A = mp.matrix([[u[2][0], u[3][0]], [u[2][1], u[3][1]]])
    sol = mp.lu_solve(A, mp.matrix([-sx, -sy]))
    return [mp.mpf(1), t, sol[0], sol[1]]


def make(a, N, t):
    L = closing(a, N, t)
    if min(L) <= 0:
        return None
    try:
        return Poly(a, N, lengths=L)
    except (AssertionError, ValueError, ZeroDivisionError):
        return None


def prong_sig(poly, cap):
    out = {}
    for C in poly.classes():
        for (i, d) in poly.slots(C):
            r = poly.prong(i, d, cap)
            out[(i, d)] = (r[0], r[1], r[3])          # kind, end side/vertex, step count
    return out


def bisect(a, N, key, t0, t1, s0, cap, tol):
    """t0 < t1, prong `key` has signature s0 at t0 and something else at t1"""
    while t1 - t0 > tol:
        tm = (t0 + t1) / 2
        P = make(a, N, tm)
        if P is None:
            return None
        sm = P.prong(*key, cap)
        if (sm[0], sm[1], sm[3]) == s0:
            t0 = tm
        else:
            t1 = tm
    return (t0 + t1) / 2


FAMILIES = [([2, 7, 8, 13], 15), ([1, 2, 5, 6], 7), ([2, 4, 5, 7], 9), ([3, 4, 5, 6], 9), ([4, 5, 9, 12], 15)]


def is_symmetric(a, L):
    """kite (two adjacent-pair reflections), parallelogram, isosceles trapezoid — by angles AND lengths"""
    k = 4
    for s in range(k):     # reflection symmetry through vertex s ⟺ a[s+i] = a[s-i], L matching
        if all(a[(s + i) % k] == a[(s - i) % k] for i in range(k)) and \
           all(abs(L[(s + i) % k] - L[(s - 1 - i) % k]) < 1e-20 for i in range(k)):
            return 'kite'
    for s in range(k):     # reflection through midpoints of side s and side s+2
        if all(a[(s + 1 + i) % k] == a[(s - i) % k] for i in range(k)) and \
           all(abs(L[(s + i) % k] - L[(s - i) % k]) < 1e-20 for i in range(k)):
            return 'trapezoid'
    if a[0] == a[2] and a[1] == a[3] and abs(L[0] - L[2]) < 1e-20:
        return 'parallelogram'
    return None


def main(ngrid=120, cap=4000, tol=mp.mpf(10) ** -40, margin=mp.mpf('0.1'), nsingle=4):
    """grid on the central (1 - 2·margin) of the valid interval (the ends degenerate: a side -> 0, walks blow up);
    every prong word change bisected; t* COINCIDENT across >= 2 prongs (to 1e-30) are the candidate cylinder
    crossings (a flow-aligned free cylinder has >= 2 boundary connections, each seen by two prongs); the area
    test runs at every cluster of >= 2 distinct pairs (by prong count >= 3, i.e. more than one pair's two ends)
    and at the first `nsingle` singletons as a baseline."""
    rows = []
    for fam in FAMILIES:
        a, N = fam[0], fam[1]
        # L2, L3 are affine in t: the valid range is the open interval where all four sides are positive
        l0, l1 = closing(a, N, mp.mpf(0)), closing(a, N, mp.mpf(1))
        lo, hi = mp.mpf(0), mp.mpf(50)
        for j in (2, 3):
            sl = l1[j] - l0[j]
            if sl > 0:
                lo = max(lo, -l0[j] / sl)
            elif sl < 0:
                hi = min(hi, -l0[j] / sl)
        w = hi - lo
        lo, hi = lo + margin * w, hi - margin * w
        if len(fam) > 2:                                   # explicit window (controls)
            lo, hi = mp.mpf(fam[2]), mp.mpf(fam[3])
        ts = [lo + (hi - lo) * mp.mpf(i) / (ngrid - 1) for i in range(ngrid)]
        grid = [(t, make(a, N, t)) for t in ts]
        grid = [(t, P) for t, P in grid if P is not None]
        sigs = [(t, P, prong_sig(P, cap)) for t, P in grid]
        print(f'{tuple(a)}/{N}: {len(grid)} valid t in ({mp.nstr(grid[0][0], 4)}, {mp.nstr(grid[-1][0], 4)})', flush=True)
        for t, P, _ in sigs[:: max(1, len(sigs) // 3)]:          # control: generic grid points
            for r in row(f'grid {tuple(a)}/{N} t={mp.nstr(t, 4)}', P, cap=cap):
                r.update(kind='grid', t=mp.nstr(t, 45), sym=is_symmetric(a, P.L), nprong=0)
                rows.append(r)
                print(f"  grid t={mp.nstr(t, 4)} p={r['p']} s'={r['sprime']} 1-R={r['one_minus_R']:+.3e}", flush=True)
        hits = []                                                  # (t*, prong)
        for (t0, P0, s0), (t1, P1, s1) in zip(sigs, sigs[1:]):
            for key in s0:
                if key in s1 and s0[key] != s1[key] and 'cap' not in (s0[key][0], s1[key][0]):
                    ts_ = bisect(a, N, key, t0, t1, s0[key], cap, tol)
                    if ts_ is not None:
                        hits.append((ts_, key))
        hits.sort()
        clusters = []
        for ts_, key in hits:
            if clusters and ts_ - clusters[-1][0] < mp.mpf(10) ** -30:
                clusters[-1][1].append(key)
            else:
                clusters.append([ts_, [key]])
        multi = [c for c in clusters if len(c[1]) >= 3]
        print(f'  {len(hits)} prong changes, {len(clusters)} distinct t*, {len(multi)} with >= 3 prongs '
              f'(sizes {sorted(len(c[1]) for c in clusters)[-8:]})', flush=True)
        # at each t*, count prongs ending at a VERTEX (one swapped pair = its two ends): a homologous partner of
        # the new connection is realised at the same t* automatically, so >= 3 vertex prongs flags it
        nv = []
        for c in clusters:
            P = make(a, N, c[0])
            sg = prong_sig(P, cap) if P is not None else {}
            nv.append(sum(1 for v in sg.values() if v[0] == 'vertex'))
            c.append(nv[-1])
        print(f'  vertex-prong counts at t*: {sorted(collections.Counter(nv).items())}', flush=True)
        todo = [c for c in clusters if len(c[1]) >= 3 or c[2] >= 3] + [c for c in clusters if c[2] < 3][:nsingle]
        for ts_, keys, _nv in todo:
            for off, kind in ((0, 'tstar'), (mp.mpf(10) ** -6, 'plus'), (-mp.mpf(10) ** -6, 'minus')):
                P = make(a, N, ts_ + off)
                if P is None:
                    continue
                for r in row(f'{kind} {tuple(a)}/{N} t={mp.nstr(ts_, 12)}', P, cap=cap):
                    r.update(kind=kind, t=mp.nstr(ts_ + off, 45), sym=is_symmetric(a, P.L), nprong=len(keys))
                    rows.append(r)
                    if kind == 'tstar' or abs(r['one_minus_R']) > 1e-20 or r['branch_capped']:
                        print(f"  {kind:5s} t={mp.nstr(ts_ + off, 12)} nprong={len(keys)} sym={r['sym']} p={r['p']} "
                              f"s'={r['sprime']} cap={r['branch_capped']} 1-R={r['one_minus_R']:+.3e}", flush=True)
        json.dump(rows, open('data/s929_asym_free.json', 'w'), indent=1)


if __name__ == '__main__':
    main(*(int(x) for x in sys.argv[1:2]))
