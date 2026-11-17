"""[Deposit] Remark 4.5 (the kite census) and Example 4.6 (parallelograms, the four ratios). Paper: "Perpendicular billiards in rational polygons of odd denominator". Research name: probes/s928_kites.py, copied verbatim except imports (deposit names), paths (data/ and figs/ under odd_polygons/) and this paragraph. Run through run_all.py; the usage lines below are the research repo's. Bracketed IDs and section marks (§sNNN) refer to the project's internal ledgers.

s928 keeper (§s928 (4)): (M⁺) FAILS at odd N for SYMMETRIC polygons — the area test of `probes/s928_poly_thmP.py`
on (a) every kite (a, b/2, c, b/2)/(Ne/2) doubled from a primitive triangle (a,b,c)/Ne across the side opposite its
even angle b, Ne ≡ 2 (mod 4), Ne <= NEMAX (default 30): odd N; (b) odd-N parallelograms (a, b, a, b)/N at side
ratios rho (centrally symmetric for every rho).  Not pre-registered (exploration after the (3,41,41,41)/63 kite);
independent check: the audit's Monte Carlo `data/scratch/s928/audit/mc_area.py` (see §s928 (4)).
PRIOR ART: §s927 (2) (even control on S), §s920 (3') (doubling across a leg), [NCYL-709]; rulings --grep kite /
parallelogram / rhombus -> none.  DATA: none.  Output: data/s928_kites.json.  ~4 min.
every kite (a, b/2, c, b/2)/(Ne/2) from a primitive triangle (a,b,c)/Ne, Ne ≡ 2 mod 4, b even (apex V1),
area test on the kite's perpendicular class (odd N) — does any odd-N kite have a free region?"""
import _setup  # noqa: F401  (deposit: this directory on the path, cwd = odd_polygons/)
import json, sys
from math import gcd
import mpmath as mp
from poly_slots import Poly
from area_test import row
mp.mp.dps = 50
NeMax = int(sys.argv[1]) if len(sys.argv) > 1 else 30
out = []
for Ne in range(6, NeMax + 1, 4):
    N = Ne // 2
    seen = set()
    for b in range(2, Ne, 2):
        for a in range(1, Ne - b):
            c = Ne - a - b
            if c < 1 or gcd(gcd(a, b), gcd(c, Ne)) != 1:
                continue
            key = (min(a, c), b, max(a, c))
            if key in seen:
                continue
            seen.add(key)
            ka = [a, b // 2, c, b // 2]
            if gcd(gcd(*ka), N) != 1:
                continue
            sa, sc = mp.sin(mp.pi * a / Ne), mp.sin(mp.pi * c / Ne)
            # triangle V0 V1 V2, angles a (V0), b (V1), c (V2) in pi/Ne: |V0V1| ∝ sin c, |V1V2| ∝ sin a
            try:
                K = Poly(ka, N, lengths=[sc, sa, sa, sc])
            except Exception as e:
                print('skip', ka, N, e); continue
            for r in row(f'kite {tuple(ka)}/{N} <- {(a, b, c)}/{Ne}', K, cap=40000):
                out.append(r)
                flag = '  <== FREE AREA' if abs(r['one_minus_R']) > 1e-20 or r['sprime'] else ''
                print(f"{r['label']:40s} p={r['p']} s'={r['sprime']} 1-R={r['one_minus_R']:+.3e} cap={r['branch_capped']}{flag}", flush=True)
for (pa, N) in [([3, 2, 3, 2], 5), ([4, 3, 4, 3], 7), ([5, 2, 5, 2], 7)]:
    for rho in ('1', '1.1', '1.37', '0.6'):
        rho = mp.mpf(rho)
        P = Poly(pa, N, lengths=[mp.mpf(1), rho, mp.mpf(1), rho])
        for r in row(f'para {tuple(pa)}/{N} rho={mp.nstr(rho, 3)}', P):
            out.append(r)
            print(f"{r['label']:40s} p={r['p']} s'={r['sprime']} 1-R={r['one_minus_R']:+.3e} rel={r['min_signed_rel']}", flush=True)
json.dump(out, open('data/s928_kites.json', 'w'), indent=1)
