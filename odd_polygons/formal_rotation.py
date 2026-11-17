"""[Deposit] Helper, imported by the scripts of this directory. Paper: "Perpendicular billiards in rational polygons of odd denominator". Research name: probes/s624_formal_rotation.py, copied verbatim except imports (deposit names), paths (data/ and figs/ under odd_polygons/) and this paragraph. Run through run_all.py; the usage lines below are the research repo's. Bracketed IDs and section marks (§sNNN) refer to the project's internal ledgers.

s624 -- the CONTINUATION question made canonical: the formal net rotation of every closed loop of a
direction class (cylinder cores, minimal-component return loops), in the bad classes of the ten and
in CP classes for comparison.

PRE-REGISTERED.  A DIRECTION has no canonical continuation off a rational triangle (the class
{+-theta_0 + 2k alpha + 2l beta} collapses a whole coset of (k, l) onto it: the Letter's 14a - 3pi vs
-16a + 5pi); an object with a finite WORD does: crossing sides s_1..s_n rotates by
2 sum_i (-1)^(i+1) phi_{s_i}, a FORMAL form 2(K alpha + L beta) (+ multiples of pi), phi_0 = 0,
phi_1 = pi - beta, phi_2 = alpha (the side lines of flow_decomp's Triangle).  Per copy crossing
the contribution is (-1)^e(g) phi_{exit side} (e = the copy's orientation), additive along a path, so
it is carried through the induction exactly like the return time (class HJK).  A closed loop is
STRUCTURAL (persists under every small deformation of the triangle: the Galperin-Stepin-Vorobets /
Galperin-Zvonkine alternating-side-sum criterion) iff (K, L) = (0, 0); in the right-triangle family
(gamma = pi/2 fixed, beta = pi/2 - alpha) iff K = L.  Otherwise it closes only because the angles
are rational: a COINCIDENCE.
CHECK (must hold on every loop): K a + L b == 0 (mod d) -- it closes at the rational triangle.
Minimal components are induced into one edge first (probes/s624_loop_holonomy.py), so their return
loops are closed curves; a cylinder's single label returns to itself exactly.
MEASURED: per class, the cylinders split structural / coincidence (count and share of mass), and the
(K, L) of the minimal components' loops.
COULD COME OUT OTHERWISE / what would be informative: CP classes all-structural and bad classes
carrying coincidence cylinders (a mechanism shape: non-CP needs a coincidence); or coincidences
everywhere (then 'coincidence' is what a rational triangle is, and this invariant explains nothing --
Hooper 2007: every periodic path in a RIGHT triangle is unstable under general deformation, so the
2-parameter criterion is expected to fire on every right-triangle loop; the family criterion is the
one to read there).
PRIOR ART: paper_orphan.md (GZ03 alternating-side-sum criterion, GSV92, Hoo07 cited; not applied to
these components); Letter §3.1 ('every return word has net rotation zero', T0's Sigma_H return words,
in the right-triangle family); rulings.py --grep 'net rotation' / 'stable periodic' / 'formal
rotation' -> none.  DATA: data_index --grep rotation -> s461_selfsim (other object); none.
RESULT (s624; data/s624_formal_rotation.json, data/s624_formal_rotation_offgrid.json): the
pre-registered mechanism shape FAILS -- coincidence cylinders are everywhere, CP and bad classes
alike (2-parameter structural cylinder mass: bad 0.078-0.464, CP 0.009-0.431, no split); in T0
(right family) the BAD class is all-structural and the CP class carries the one coincidence pair
(+-8, -+7).  UNEXPECTED, and the finding: every minimal component's return loops have formal
rotations on ONE LINE (K, L) ~ (k0, l0) -- 13/13 rows (the ten, (5,7,11)'s two parallel classes)
-- and (k0, l0) is the SHORTEST vector of the closure lattice {K a + L b == 0 mod d} in all 11
triangles: the torus is structural along the one-parameter family k0 alpha + l0 beta = const.
CONTROL: off-grid scalene tori ((2,5,6)) rank 2 in 8/8; right / isosceles off-grid mostly rank 1
along their symmetry family ((4,5,9), T0 K = L; (7,7,16) K = -L; 3 rank-2 exceptions).
-> [NCYL-633].
"""
import _setup  # noqa: F401  (deposit: this directory on the path, cwd = odd_polygons/)
import json
import sys
import time
from collections import Counter

import flow_decomp as fd
import loop_holonomy as lh

PHI = {0: (0, 0), 1: (0, -1), 2: (1, 0)}                 # formal (alpha, beta) part of phi_s


class HJK(lh.HJ):
    __slots__ = ('K',)

    def __init__(self, h, J, K):
        super().__init__(h, J)
        self.K = K

    def __add__(self, o):
        r = lh.HJ.__add__(self, o)
        return HJK(r.h, r.J, (self.K[0] + o.K[0], self.K[1] + o.K[1]))

    def __rmul__(self, k):
        r = lh.HJ.__rmul__(self, k)
        return HJK(r.h, r.J, (k * self.K[0], k * self.K[1]))


def classify(abc, om, cap=20000, extra=200000):
    T = fd.Triangle(abc)
    F = T.F
    pieces, lengths, top, bot, meta = fd.build_iet(T, om)
    J0 = {p: T.red(T.g_apply(pieces[p][3][0], T.t[pieces[p][4][1]])) for p in range(len(pieces))}
    K0 = {}
    for p in range(len(pieces)):
        (g, _), (_, s2) = pieces[p][3], pieces[p][4]
        sg = -1 if g[1] else 1
        K0[p] = (sg * PHI[s2][0], sg * PHI[s2][1])
    segs = {}
    for dom, ln, img, key, key2 in pieces:
        end = F.add(dom, ln)
        if key not in segs:
            segs[key] = [dom, end]
        else:
            if F.cmp(dom, segs[key][0]) < 0:
                segs[key][0] = dom
            if F.cmp(end, segs[key][1]) > 0:
                segs[key][1] = end
    seglist = list(segs.values())

    def inside_one(a, ln):
        b = F.add(a, ln)
        return any(F.cmp(lo, a) <= 0 and F.cmp(b, hi) <= 0 for lo, hi in seglist)
    Dc = lh.PD(F, lengths, J0)
    Dc.h = {l: HJK(1, J0[l], K0[l]) for l in lengths}
    comps = Dc.run_tracked(top, bot, cap + extra, inside_one)
    L = meta['L']
    out = []
    a, b, c = abc
    d = a + b + c
    for comp, st in comps:
        mass = F.zero
        for l in comp.top:
            mass = F.add(mass, F.smul(Dc.h[l].h, Dc.lam[l]))
        Ks = [Dc.h[l].K for l in comp.top]
        for K in Ks:
            assert (K[0] * a + K[1] * b) % d == 0, ('loop does not close formally', abc, K)
        out.append({'kind': comp.kind, 'n': len(comp.top), 'mass': F.fval(mass) / F.fval(L),
                    'K': [list(K) for K in Ks],
                    'structural2': [K == (0, 0) for K in Ks],
                    'structural_right': [K[0] == K[1] for K in Ks]})
    return out


def main():
    LIST = [(7, 8, 15), (7, 7, 16), (4, 4, 7), (13, 20, 23), (15, 30, 47), (11, 45, 50),
            (23, 29, 58), (19, 38, 55), (11, 22, 103), (30, 41, 67)]
    CP_EXTRA = [(3, 4, 5), (2, 5, 6), (4, 5, 9), (3, 5, 7), (5, 7, 11), (2, 3, 4)]
    which = sys.argv[1:] or ['all']
    rows = []
    for abc in LIST + CP_EXTRA:
        dirs = [('perp', s, j) for s, j in enumerate(fd.perp_dirs(abc))] + \
               [('par', s, j) for s, j in enumerate(fd.side_dirs(abc))]
        seen = set()
        for name, s, j in dirs:
            # one direction per class: classes = j mod 4 up to sign (rotations by zeta^4, reflection)
            cls = min(j % 4, (-j) % 4)
            if cls in seen:
                continue
            seen.add(cls)
            t0 = time.time()
            try:
                comps = classify(abc, j)
            except (AssertionError, RuntimeError) as e:
                print(abc, name, s, 'ERROR', str(e)[:100], flush=True)
                continue
            cyl = [c for c in comps if c['kind'] == 'cylinder']
            mins = [c for c in comps if c['kind'].startswith('no_periodic')]
            und = [c for c in comps if c['kind'] == 'undetermined']
            s2 = sum(c['structural2'][0] for c in cyl)
            sr = sum(c['structural_right'][0] for c in cyl)
            ms2 = sum(c['mass'] for c in cyl if c['structural2'][0])
            msr = sum(c['mass'] for c in cyl if c['structural_right'][0])
            row = {'abc': abc, 'dir': f'{name}{s}', 'class': cls, 'ncyl': len(cyl), 'nmin': len(mins),
                   'nund': len(und), 'cyl_structural2': s2, 'cyl_structural_right': sr,
                   'mass_structural2': ms2, 'mass_structural_right': msr,
                   'cyl_K': Counter(tuple(c['K'][0]) for c in cyl).most_common(8),
                   'min_K': [c['K'] for c in mins]}
            rows.append(row)
            print(abc, f'{name}{s}', 'cyl', len(cyl), 'min', len(mins), 'und', len(und),
                  '| structural cyl (2-param) %d, mass %.3f | (right-family) %d, mass %.3f' % (s2, ms2, sr, msr),
                  '| cyl K', row['cyl_K'][:5], '| min K', [sorted(set(map(tuple, k))) for k in row['min_K']][:2],
                  f'{time.time() - t0:.0f}s', flush=True)
    json.dump(rows, open('data/s624_formal_rotation.json', 'w'), indent=1, default=str)


if __name__ == '__main__':
    main()
