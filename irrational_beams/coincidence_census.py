"""[Deposit] Section 5: the floating-point finiteness retest (test). Paper: "Finite perpendicular-beam partitions at irrational angles". Research name: probes/s973_coincidence_census.py, copied verbatim except imports (deposit names), paths (data/ and figs/ under irrational_beams/) and this paragraph. Run through run_all.py.

s973 — which irrational α give a FINITE L1 perpendicular-beam partition? Census over FIRST-LEG VERTEX COINCIDENCES.

PRE-REGISTERED (s973, queue s972-b). §s972 (5)/(6): both known irrational finite points are first-leg vertex coincidences in the
unfolding around O (images R_k at height h(R_k) = cot α sin 2kα, A_{2k+1} at h(A_{2k+1}) = D_k(α) = sin((2k+1)α)/sin α): the family
α_n is h(R_n) = ½, the point D₂ = ½ is h(A_5) = ½ (certified exactly s973, 6 branches, probes/exact_partition.py).
CENSUS: every root α ∈ (0.03, π/2 − 0.005) of
   (T1) h_v = ½,   (T2) h_v + h_w = 1 (v ≠ w),   (T3) h_v = h_w (v ≠ w),
v, w ∈ {R_k : 1 ≤ k ≤ K} ∪ {A_{2k+1} : 0 ≤ k ≤ K}, K = 12, rational α/π (x = 2α/π with denominator ≤ 400) dropped.
Each root: float trace (s959 trace_nb) at ε = 1e-6 (screen), then ε = 1e-9 for screen-finite ones; 'finite' := no capped piece and
same branch count at both ε. Hypotheses:
  (H1) the finite points are SPARSE among the coincidences (§s972 (6): D₃ = ½, D₄ = ½ infinite-looking) — most roots infinite.
  (H2) every finite point has σ₁ = 4cos²α a unit (Theorem I′, [BEAM-INF-UNIT]: finite ⇒ unit) — a float findpoly identification,
       to be certified by exact_partition. A non-unit finite point would contradict a PROVED theorem ⇒ instrument error.
  (H3) CONTROL: no finite point in (π/4, π/3) ([Tro05] Lemma 12 via [LIT-060]: n = ∞ there for irrational α) and none on the
       leg-swapped side where the two-legs corollary forbids it (4 − σ₁ a unit too — both legs finite is impossible).
Could come out otherwise: finite points abundant (H1 fails), or a finite point in (π/4, π/3) (H3 fires: instrument bug).
PRIOR ART: §s972 (1)–(6) (general_triangles.md); rulings --grep 'finite partition' / 'unit locus' → [CPRIME-025] (T_G, other object),
[NCYL-744], [LIT-060] (Troubetzkoy: n = ∞ on (π/4, π/3)); nothing on a coincidence census.
DATA: data_index --grep unit → data/s972_unit_family.json (the family; reproduced by the census as a control: α_n must appear).
Output data/coincidence_census.json.
Usage: python probes/coincidence_census.py [K] [NPROC]
"""
import _setup  # noqa: F401  (deposit: this directory on the path, cwd = irrational_beams/)
import json
import math
import os
import sys
import time
from fractions import Fraction
from multiprocessing import Pool

import numpy as np
from scipy.optimize import brentq

sys.path[:0] = [os.path.dirname(os.path.abspath(__file__))]
K = int(sys.argv[1]) if len(sys.argv) > 1 else 12
NPROC = int(sys.argv[2]) if len(sys.argv) > 2 else 3
OUT = 'data/s973_coincidence_census_K%d.json' % K


def verts(K):
    v = [('R%d' % k, lambda a, k=k: math.cos(a) / math.sin(a) * math.sin(2 * k * a)) for k in range(1, K + 1)]
    v += [('A%d' % (2 * k + 1), lambda a, k=k: math.sin((2 * k + 1) * a) / math.sin(a)) for k in range(0, K + 1)]
    return v


def roots(f, lo=0.03, hi=math.pi / 2 - 0.005, n=20000):
    g = np.linspace(lo, hi, n)
    vals = np.array([f(a) for a in g])
    out = []
    for i in range(n - 1):
        if vals[i] == 0:
            out.append(g[i])
        elif vals[i] * vals[i + 1] < 0:
            try:
                out.append(brentq(f, g[i], g[i + 1], xtol=1e-15, rtol=1e-15))
            except ValueError:
                pass
    return out


def candidates(K):
    V = verts(K)
    C = []
    for i, (n1, h1) in enumerate(V):
        for a in roots(lambda a: h1(a) - 0.5):
            C.append((a, 'T1', n1, ''))
        for n2, h2 in V[i + 1:]:
            for a in roots(lambda a: h1(a) + h2(a) - 1):
                C.append((a, 'T2', n1, n2))
            for a in roots(lambda a: h1(a) - h2(a)):
                C.append((a, 'T3', n1, n2))
    # drop rational x = 2α/π (small denominator) and dedupe
    keep = {}
    for a, t, n1, n2 in C:
        x = 2 * a / math.pi
        fr = Fraction(x).limit_denominator(400)
        if abs(float(fr) - x) < 1e-12:
            continue
        key = round(a, 11)
        keep.setdefault(key, (a, []))[1].append('%s:%s%s' % (t, n1, ('=' + n2) if n2 else ''))
    return sorted(keep.values())


def test(item):
    a, tags = item
    import launch_side as L, beam_families as B
    T = L.table('right', a, 0)
    r6 = B.trace_nb(T, eps_min=1e-6, budget=3 * 10 ** 6)
    row = dict(alpha=a, x=2 * a / math.pi, tags=tags, n6=len(r6['branches']), c6=len(r6['capped']))
    if r6['capped'] or not r6['branches']:
        row['finite'] = False
        return row
    r9 = B.trace_nb(T, eps_min=1e-9, budget=3 * 10 ** 7)
    row.update(n9=len(r9['branches']), c9=len(r9['capped']))
    row['finite'] = (not r9['capped']) and len(r9['branches']) == len(r6['branches'])
    if row['finite']:
        row['branches'] = [(float(b[0]), float(b[1]), int(b[2])) for b in sorted(r9['branches'])]
        row['sigma1'] = 4 * math.cos(a) ** 2
    return row


if __name__ == '__main__':
    t0 = time.time()
    C = candidates(K)
    print(f'{len(C)} candidate roots (K = {K}), {time.time() - t0:.0f}s', flush=True)
    with Pool(NPROC) as P:
        rows = P.map(test, C, chunksize=8)
    fin = [r for r in rows if r['finite']]
    json.dump(dict(K=K, n_candidates=len(C), rows=rows), open(OUT, 'w'), indent=0)
    print(f'{len(fin)} finite of {len(rows)}; {time.time() - t0:.0f}s')
    for r in fin:
        print(f"  x={r['x']:.10f} α={r['alpha']:.10f} σ1={r['sigma1']:.10f} n={r['n9']} nref={[b[2] for b in r['branches']][:10]} "
              f"tags={r['tags'][:4]}")
