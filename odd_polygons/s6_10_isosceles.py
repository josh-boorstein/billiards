"""[Deposit] Remark 6.10(b): the isosceles control. Paper: "Perpendicular billiards in rational polygons of odd denominator". Research name: probes/s932_iso_test.py, copied verbatim except imports (deposit names), paths (data/ and figs/ under odd_polygons/) and this paragraph. Run through run_all.py; the usage lines below are the research repo's. Bracketed IDs and section marks (§sNNN) refer to the project's internal ledgers.

s932 keeper (queue `s931-c`): the FALSIFICATION TEST of [POLY-TRAP-ACYC] on the isosceles trapezoids (1,1,N−1,N−1)/N —
acyclic (chain 0/1), but their double walls carry two swapped pairs of EQUAL length (mirror diagonals, §s931 (3) control (a)),
so Theorem P is silent there; the theorem predicts no free area anyway. Every p ≥ 2 wall in D ∈ [1/20, 6] (exact member census)
is evaluated by the area test `s928_poly_thmP.row` (R = 1 ⟺ no free area). Could have come out otherwise: a free cylinder
bounded by the two mirror diagonals (the kite / parallelogram failure mode of §s928 (4)).
Result (s932): N = 5, 7, 9 — two p = 2 walls each (D = 1 and D ≈ 4.6–4.9), min signed relation 0, 1 − R = 0 to 50 digits.
PRIOR ART: §s931 (3) control (a) (*whether those shapes carry free area was not checked*).  DATA: none.
Usage: python probes/s932_iso_test.py 7:1,1,6,6 5:1,1,4,4 9:1,1,8,8
"""
import _setup  # noqa: F401  (deposit: this directory on the path, cwd = odd_polygons/)
import sys, mpmath as mp
import multistage_walls as M, free_scan as FS
from area_test import row
mp.mp.dps = 50
for tag in sys.argv[1:]:
    N=int(tag.split(':')[0]); a=[int(x) for x in tag.split(':')[1].split(',')]
    groups,fams,G=M.member_census(tag, __import__("fractions").Fraction(1,20), 6)
    c=G.c_slice; dbl=[(D,js) for D,js in groups if len(js)>=4]
    print(tag,'members',sum(len(js) for D,js in groups),'distinct D',len(groups),'p>=2 walls',len(dbl), flush=True)
    for D,js in dbl[:4]:
        t=1/(D.mp(50)+c.mp(50))
        P=FS.make((a,N,[1,0]),t); r=row('dbl',P,cap=20000)[0]
        print('  D=%.6f t=%.8f'%(float(D),float(t)),'p',r['p'],'min_rel',r['min_signed_rel'],'1-R',r['one_minus_R'],'capped',r['branch_capped'],'sprime',r['sprime'], flush=True)
