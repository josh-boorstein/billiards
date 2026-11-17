"""[Deposit] Remark 6.10(a), (b): the wall families at a0 = 1, N ≤ 13, and the isosceles control. Paper: "Perpendicular billiards in rational polygons of odd denominator". Research name: probes/s931_twist_walls.py, copied verbatim except imports (deposit names), paths (data/ and figs/ under odd_polygons/) and this paragraph. Run through run_all.py; the usage lines below are the research repo's. Bracketed IDs and section marks (§sNNN) refer to the project's internal ledgers.

s931 keeper (queue `s930-a`): THE TWIST-WALL PROOF of (M⁺) at every side length of the parallel-cut trapezoids `(1,2,5,6)/7`,
`(1,3,4,6)/7` — exact cap automaton ⊕ wall families in the strip length D ⊕ every double-wall coincidence solved over ℤ.

SETTING (s931, hand). The slice L = (1, t, 1 − ct, ·) (c = 1 + 2cos(2π/7)) is the triangle (1,2,4)/7 resp. (1,3,3)/7 cut parallel to
V0V1. Rescaled by 1/t it is RIGHT CAP R (B_R, V1, V2) ∪ a STRIP [x_L, x_R] × [0, h] (h = sin(a_1π/7), length D = 1/t − c) ∪ LEFT CAP L
(V0, B_L, V3); the caps do not depend on t. A class ray crossing the strip in direction θ moves vertically by D·tanθ, folded into
[0, h]; so t enters ONLY through D, D ∈ (0, ∞) ⟷ t ∈ (0, t_max).
  The CAP MAPS (exact, here): a ray entering a cap through its virtual side x = x_{R/L} at height y ∈ (0, h) in a class direction ends
  perpendicularly, hits a vertex (finitely many y), or leaves through the virtual side at height c + λy. LEMMA (checked exactly
  here, asserted): L absorbs — every entering ray ends perpendicularly in L except at finitely many heights (vertex hits). Hence every
  prong is (stage 0: cap path from its vertex) → ≤ 1 transit → R → ≤ 1 transit → L, at most two transits, for EVERY D.
  WALLS: a prong reaches a vertex iff D satisfies W·D = P + Q·N, N ∈ ℤ^{1 or 2} (the fold counts), with explicit realisability
  inequalities. Two different connections at one D ⟺ an integer solution of one ℚ(cos π/14)-linear equation = 6 rational equations in
  ≤ 4 integer unknowns — solved exactly for every pair of families.

PRE-REGISTERED (s931, before the families were enumerated).
  Claim under test (T): (i) L absorbs (both families); (ii) no stage-0 path ends at a vertex (no persistent pair); (iii) every pair of
  wall families either has IDENTICAL D-sets and is a connection with its own reversal, or meets in finitely many D, and at every
  realised double-wall D the exact pair lengths satisfy no signed relation. Then Theorem P ([POLY-THM-P]) gives (M⁺) ⊕ CP at EVERY
  t ∈ (0, t_max): off the walls p = 0, on one wall p = 1, at the solved coincidences p ≥ 2 with no relation.
  Refuted by: an exit piece of L; a stage-0 vertex end; two non-reverse families with identical D-sets (a double wall on a whole
  progression); or a realised coincidence with a relation.
  Could come out otherwise: (iii) is a rank question — two families whose K-coefficients are ℚ-dependent meet in an infinite lattice
  of D; nothing forces full rank.
  Control: the cap automaton reproduces the s931 float twist walker and the real trapezoid (tw_validate, 6/6 rows, 0 mismatches);
  every realised coincidence is re-derived by the exact tracker `poly_exact.Fam.prong` on the real trapezoid at t = 1/(D + c); the wall
  families are checked against the exact breakpoints of `probes/s931_trapezoid_exact.py` on [1/20, 0.445] (every tracked breakpoint
  must lie on a family, and every family member in that range must be a tracked breakpoint).
PRIOR ART: §s930 (2), (3) ⊕ [NCYL-721] (b) (the families, float); `probes/s931_trapezoid_exact.py` (this session, compact range);
  `rulings.py --grep` *twist* → (s931: hits are right-triangle twist maps, not polygon strips); *cap map* → none; *strip length* → none.
DATA: none (`data_index.py --grep trapezoid` → none).
Output: data/s931_twist_walls_<fam>.json.  Usage: python probes/s931_twist_walls.py FAM   (FAM 1256 | 1346 | N:a0,a1,a2,a3 — a trapezoid, a0 + a3 = a1 + a2 = N, 2a1 < N)
  Series (s931, `data/scratch/s931/series.sh`): every (1, a, N − a, N − 1)/N, N = 5 … 13; crosscheck(tag, T0, T1, fn) is the control.
"""
import _setup  # noqa: F401  (deposit: this directory on the path, cwd = odd_polygons/)
import itertools
import json
import math
import sys
from fractions import Fraction as Fr

import mpmath

from poly_exact import field, _Sign, Aff, Fam

N = 7
FAMS = {'1256': [1, 2, 5, 6], '1346': [1, 3, 4, 6]}


class Geo:
    def __init__(self, a):
        self.a = a
        K = self.K = field(N)
        self.sg = _Sign(K)
        th = K.gen
        c = [K.of(2), th]
        for j in range(2, 4 * N + 1):
            c.append(th * c[-1] - c[-2])
        self._c2 = c
        idx = [0]
        for j in range(1, 4):
            idx.append(idx[-1] + 2 * (N - a[j]))
        self.idx = [x % (4 * N) for x in idx]
        self.m = [x % (2 * N) for x in self.idx]
        A0, A1 = 2 * a[0], 2 * a[1]
        self.h = h = self.sin(A1)
        Z = K.zero
        P = lambda x, y: (Aff(K.of(x), Z), Aff(K.of(y), Z))
        self.R = dict(verts=[P(Z, Z), P(self.cos(A1), Z), P(Z, h)], labels=['B', 'V1', 'V2'], sides=[0, 1, 'v'])
        cotA0 = self.cos(A0) / self.sin(A0)
        self.L = dict(verts=[P(-h * cotA0, Z), P(Z, Z), P(Z, h)], labels=['V0', 'B', 'V3'], sides=[0, 'v', 3])
        self._maps = {}

    def cos(self, j):
        return self._c2[j % (4 * N)] * Fr(1, 2)

    def sin(self, j):
        return self.cos(N - j)

    def tan_abs(self, d):
        return abs(self.sin(d) / self.cos(d))

    def trace(self, cap, d, start, y0, entry, cap_steps=200):
        K, sg = self.K, self.sg
        y0 = K.of(y0) if not hasattr(y0, 'c') else y0
        ux, uy = self.cos(d), self.sin(d)
        S = start
        copy = list(cap['verts'])
        labels, sides = cap['labels'], cap['sides']
        eps, r, ent = 1, 0, entry
        word, forms = [], []

        def form(W):
            return (W[1] - S[1]).scale(ux) - (W[0] - S[0]).scale(uy)
        for step in range(cap_steps):
            fs = []
            for j in range(3):
                if step == 0 and isinstance(ent, tuple) and ent[1] == j:
                    fs.append(None)
                    continue
                f = form(copy[j])
                fs.append(f)
                forms.append(f)
            sv = [None if f is None else sg(f.at(y0)) for f in fs]
            zero = [j for j in range(3) if sv[j] == 0]
            dcopy = (eps * (d - r)) % (4 * N)
            if zero:
                assert len(zero) == 1
                W = copy[zero[0]]
                ln = (W[0] - S[0]).scale(ux) + (W[1] - S[1]).scale(uy)
                out = dict(kind='vertex', vertex=labels[zero[0]], length=ln, arr=dcopy, steps=step)
                break
            entry_side = ent if not isinstance(ent, tuple) else None
            ex = None
            for j in range(3):
                if j == entry_side:
                    continue
                s1, s2 = sv[j], sv[(j + 1) % 3]
                if s1 is None or s2 is None:
                    continue
                if s1 != s2:
                    assert ex is None
                    ex = j
            assert ex is not None
            word.append(sides[ex])
            if sides[ex] == 'v':
                fa, fb = fs[ex], fs[(ex + 1) % 3]
                den = fa - fb
                assert den.const()
                lam = fa.scale(K.one / den.c0)
                A, B = cap['verts'][ex], cap['verts'][(ex + 1) % 3]
                hy = A[1] + lam.scale((B[1] - A[1]).c0)
                Xx = copy[ex][0] + lam.scale((copy[(ex + 1) % 3][0] - copy[ex][0]).c0)
                Xy = copy[ex][1] + lam.scale((copy[(ex + 1) % 3][1] - copy[ex][1]).c0)
                ln = (Xx - S[0]).scale(ux) + (Xy - S[1]).scale(uy)
                out = dict(kind='exit', y=hy, d=dcopy, length=ln, steps=step)
                break
            if (dcopy - self.m[sides[ex]] - N) % (2 * N) == 0:
                out = dict(kind='perp', side=sides[ex], steps=step, d=dcopy)
                break
            phi = (r + eps * self.idx[sides[ex]]) % (4 * N)
            A = copy[ex]
            c2, s2 = self.cos(2 * phi), self.sin(2 * phi)

            def refl(P):
                dx, dy = P[0] - A[0], P[1] - A[1]
                return (A[0] + dx.scale(c2) + dy.scale(s2), A[1] + dx.scale(s2) - dy.scale(c2))
            copy = [refl(P) for P in copy]
            eps, r = -eps, (2 * phi - r) % (4 * N)
            ent = ex
        else:
            raise RuntimeError('cap trace capped')
        lo = hi = None
        for f in forms:
            if f.c1.iszero():
                continue
            if sg(f.at(y0)) == 0:
                continue
            rt = f.root()
            s = sg(rt - y0)
            if s > 0 and (hi is None or sg(rt - hi) < 0):
                hi = rt
            elif s < 0 and (lo is None or sg(rt - lo) > 0):
                lo = rt
        out.update(word=tuple(word), lo=lo, hi=hi)
        return out

    def cap_map(self, which, d):
        if (which, d) in self._maps:
            return self._maps[(which, d)]
        K, sg = self.K, self.sg
        cap = self.R if which == 'R' else self.L
        vs = cap['sides'].index('v')
        start = (Aff(K.zero, K.zero), Aff(K.zero, K.one))
        pieces, y, left = [], Fr(1, 10**12), K.zero
        while True:
            o = self.trace(cap, d, start, y, vs)
            if o['lo'] is not None and sg(o['lo'] - left) > 0:
                y = Fr(mpmath.nstr((left.mp(30) + K.of(y).mp(30)) / 2, 28))
                continue
            o['lo'] = left
            pieces.append(o)
            hi = o['hi']
            if hi is None or sg(hi - self.h) >= 0:
                o['hi'] = self.h
                break
            o['bp'] = self.trace(cap, d, start, hi, vs)
            left = hi
            y = Fr(mpmath.nstr(hi.mp(30) + mpmath.mpf(10) ** -15, 28))
        self._maps[(which, d)] = pieces
        return pieces


def mov_right(G, d):
    return G.sg(G.cos(d)) > 0


def flip(d):
    return (-d) % (4 * N)


def up(d):
    return d if math.sin(math.pi * d / (2 * N)) > 0 else flip(d)


def run(tag):
    global N
    if tag in FAMS:
        a = FAMS[tag]
        N = 7
    else:                                   # 'N:a0,a1,a2,a3'
        nn, aa = tag.split(':')
        N, a = int(nn), [int(x) for x in aa.split(',')]
    assert a[0] + a[3] == N and a[1] + a[2] == N and 2 * a[1] < N, 'trapezoid with a right-triangle right cap expected'
    G = Geo(a)
    K, sg, h = G.K, G.sg, G.h
    E = Fam.closing(a, N, [(1, 0), (0, 1)])
    G.c_slice = -E.L[2].c1                  # top side = 1 − c t, so D = top/t = 1/t − c
    slots = E.slots()
    rep = {'fam': a, 'checks': {}, 'families': [], 'coincidences': [], 'fail': []}
    # (i) L absorbs
    absorbs = True
    for d in range(1, 4 * N, 2):
        if mov_right(G, d) or sg(G.cos(d)) == 0:
            continue
        for p in G.cap_map('L', d):
            if p['kind'] != 'perp':
                absorbs = False
    rep['checks']['L_absorbs'] = absorbs
    if not absorbs:
        rep['fail'].append('L does not absorb')
    # stage 0
    capR, capL = G.R, G.L
    stage0 = {}
    for (i, d) in slots:
        cosd = sg(G.cos(d))
        if i == 1:
            o = G.trace(capR, d, capR['verts'][1], 0, ('vtx', 1))
        elif i == 0:
            o = G.trace(capL, d, capL['verts'][0], 0, ('vtx', 0))
        elif i == 2:       # V2: into R if moving right, strip if moving left, down = perp at B_R
            if cosd > 0:
                o = G.trace(capR, d, capR['verts'][2], 0, ('vtx', 2))
            elif cosd == 0:
                o = dict(kind='perp', side=0)
            else:
                o = dict(kind='exit', y=Aff(h, K.zero), d=d, length=Aff(K.zero, K.zero))
        else:              # V3
            if cosd < 0:
                o = G.trace(capL, d, capL['verts'][2], 0, ('vtx', 2))
            elif cosd == 0:
                o = dict(kind='perp', side=0)
            else:
                o = dict(kind='exit', y=Aff(h, K.zero), d=d, length=Aff(K.zero, K.zero))
        stage0[(i, d)] = o
        if o['kind'] == 'vertex':
            rep['fail'].append(('stage-0 vertex end (persistent pair)', i, d, o['vertex']))
    rep['checks']['stage0'] = {f'{k[0]},{k[1]}': v['kind'] for k, v in stage0.items()}
    VLAB = {'V0': 0, 'V1': 1, 'V2': 2, 'V3': 3}
    fams = []
    for (i, d), o in stage0.items():
        if o['kind'] != 'exit':
            continue
        y0, d0 = o['y'].c0, o['d']
        X = 'L' if not mov_right(G, d0) else 'R'
        s0 = sg(G.sin(d0))
        T0 = G.tan_abs(d0)
        for e1 in (1, -1):
            d1 = d0 if e1 == 1 else flip(d0)
            targets = [(h, 'V2' if X == 'R' else 'V3', up(d1))] if e1 == 1 else []   # the top corner is reached moving UP
            for p in G.cap_map(X, d1):
                if 'bp' in p and p['bp']['kind'] == 'vertex' and p['bp']['vertex'] != 'B':
                    targets.append((p['hi'], p['bp']['vertex'], p['bp']['arr']))
            for yt, lab, arr in targets:
                fams.append(dict(slot=(i, d), X=X, transits=1, e=(e1,), target=lab, arr=arr, W=K.of(e1 * s0) * T0,
                                 P=yt - e1 * y0, Q=[-2 * h], cons=[('y1', e1, y0, s0 * e1, T0, None, None)]))
            if X == 'L':
                continue
            for p in G.cap_map('R', d1):
                if p['kind'] != 'exit':
                    continue
                cR, lam, d1x = p['y'].c0, p['y'].c1, p['d']
                s1, T1 = sg(G.sin(d1x)), G.tan_abs(d1x)
                for e2 in (1, -1):
                    d2 = d1x if e2 == 1 else flip(d1x)
                    targets = [(h, 'V3', up(d2))] if e2 == 1 else []
                    for q in G.cap_map('L', d2):
                        if 'bp' in q and q['bp']['kind'] == 'vertex' and q['bp']['vertex'] != 'B':
                            targets.append((q['hi'], q['bp']['vertex'], q['bp']['arr']))
                    for yt, lab, arr in targets:
                        W = K.of(e2) * (lam * (e1 * s0) * T0 + K.of(s1) * T1)
                        P = yt - e2 * cR - e2 * lam * e1 * y0
                        Q = [-2 * h * e2 * lam, -2 * h]
                        fams.append(dict(slot=(i, d), X='R', transits=2, e=(e1, e2), target=lab, arr=arr, W=W, P=P, Q=Q,
                                         piece=(p['lo'], p['hi']), y0=y0, s0e1=s0 * e1, T0=T0, d1=d1))
    rep['n_families'] = len(fams)
    assert all(not F['W'].iszero() for F in fams), 'a D-independent wall (W = 0)'   # audit s931 R3

    def D_of(F, Nv):
        return (F['P'] + sum((q * n for q, n in zip(F['Q'], Nv)), K.zero)) / F['W']

    def realised(F, Nv):
        D = D_of(F, Nv)
        if sg(D) <= 0:
            return False, D
        if F['transits'] == 1:
            return True, D            # y1 = target ∈ (0, h] by construction
        y1 = F['e'][0] * F['y0'] + K.of(F['s0e1']) * F['T0'] * D + 2 * h * Nv[0]
        lo, hi = F['piece']
        return (sg(y1 - lo) > 0 and sg(hi - y1) > 0), D

    def vec(x):
        c = list(x.c) + [Fr(0)] * (K.deg - len(x.c))
        return c

    def solve(cols, rhs):
        """exact solve of Σ x_j cols_j = rhs over ℚ → None (inconsistent) or (x0, kernel basis)"""
        m = len(cols)
        A = [[vec(cj)[r] for cj in cols] + [vec(rhs)[r]] for r in range(K.deg)]
        piv, row = [], 0
        for col in range(m):
            pr = next((r for r in range(row, K.deg) if A[r][col] != 0), None)
            if pr is None:
                continue
            A[row], A[pr] = A[pr], A[row]
            pv = A[row][col]
            A[row] = [x / pv for x in A[row]]
            for r in range(K.deg):
                if r != row and A[r][col] != 0:
                    f = A[r][col]
                    A[r] = [x - f * y for x, y in zip(A[r], A[row])]
            piv.append(col)
            row += 1
        for r in range(row, K.deg):
            if A[r][m] != 0:
                return None
        x0 = [Fr(0)] * m
        for r, col in enumerate(piv):
            x0[col] = A[r][m]
        ker = []
        for j in range(m):
            if j in piv:
                continue
            v = [Fr(0)] * m
            v[j] = Fr(1)
            for r, col in enumerate(piv):
                v[col] = -A[r][j]
            ker.append(v)
        return x0, ker

    def lin_constraints(F, Nv_aff):
        """constraints of family F as K-affine forms in the lattice parameter k: list of (a, b) meaning a + b k > 0.
        Nv_aff = list of (n0, n1) rationals: N_i = n0 + n1 k."""
        Pn = F['P'] + sum((q * n0 for q, (n0, n1) in zip(F['Q'], Nv_aff)), K.zero)
        Pk = sum((q * n1 for q, (n0, n1) in zip(F['Q'], Nv_aff)), K.zero)
        D0, D1 = Pn / F['W'], Pk / F['W']
        out = [(D0, D1)]
        if F['transits'] == 2:
            y10 = F['e'][0] * F['y0'] + K.of(F['s0e1']) * F['T0'] * D0 + 2 * h * Nv_aff[0][0]
            y11 = K.of(F['s0e1']) * F['T0'] * D1 + 2 * h * Nv_aff[0][1]
            lo, hi = F['piece']
            out.append((y10 - lo, y11))
            out.append((hi - y10, -y11))
        return out, (D0, D1)

    def k_range(cons):
        """integer k with a + b k > 0 for all constraints: (lo, hi) with None = unbounded, or 'empty'"""
        lo = hi = None
        for a_, b_ in cons:
            sb = sg(b_)
            if sb == 0:
                if sg(a_) <= 0:
                    return 'empty'
                continue
            r = -a_ / b_
            rf = r.mp(40)
            if sb > 0:       # k > r
                kk = int(mpmath.floor(rf)) + 1
                while sg(a_ + b_ * (kk - 1)) > 0:
                    kk -= 1
                while sg(a_ + b_ * kk) <= 0:
                    kk += 1
                lo = kk if lo is None else max(lo, kk)
            else:            # k < r
                kk = int(mpmath.ceil(rf)) - 1
                while sg(a_ + b_ * (kk + 1)) > 0:
                    kk += 1
                while sg(a_ + b_ * kk) <= 0:
                    kk -= 1
                hi = kk if hi is None else min(hi, kk)
        if lo is not None and hi is not None and lo > hi:
            return 'empty'
        return (lo, hi)

    from math import lcm
    n_none = n_unique = n_line = 0
    infinite, coinc = [], []
    for u, v in itertools.combinations(range(len(fams)), 2):
        F, F2 = fams[u], fams[v]
        if F['slot'] == F2['slot']:
            continue                       # one prong ends at one vertex
        cols = [q * F2['W'] for q in F['Q']] + [-(q * F['W']) for q in F2['Q']]
        rhs = F2['P'] * F['W'] - F['P'] * F2['W']
        sol = solve(cols, rhs)
        if sol is None:
            n_none += 1
            continue
        x0, ker = sol
        n1 = len(F['Q'])
        f = len(ker)
        if f == 0:
            n_unique += 1
        else:
            n_line += 1
        m = len(x0)
        # x = x0 + Σ λ_j ker_j, λ ∈ ℤ^f (the free coordinates); the other coordinates integral ⟺ λ in residue classes mod per
        per = 1
        for i in range(m):
            per = lcm(per, x0[i].denominator, *[kv[i].denominator for kv in ker])
        residues = []
        for r in itertools.product(range(per), repeat=f):
            if all((x0[i] + sum((r[j] * ker[j][i] for j in range(f)), Fr(0))).denominator == 1 for i in range(m)):
                residues.append(r)
        for r in residues:
            # N_i = n0_i + Σ_j per·ker_j[i]·k_j
            base = [x0[i] + sum((r[j] * ker[j][i] for j in range(f)), Fr(0)) for i in range(m)]
            gens = [[per * ker[j][i] for j in range(f)] for i in range(m)]

            def cons_of(Fx, off):
                """constraints a + b·k > 0, b ∈ K^f"""
                Pn = Fx['P'] + sum((q * base[off + i] for i, q in enumerate(Fx['Q'])), K.zero)
                Pk = [sum((q * gens[off + i][j] for i, q in enumerate(Fx['Q'])), K.zero) for j in range(f)]
                D0, D1 = Pn / Fx['W'], [p / Fx['W'] for p in Pk]
                out = [(D0, D1)]
                if Fx['transits'] == 2:
                    y10 = Fx['e'][0] * Fx['y0'] + K.of(Fx['s0e1']) * Fx['T0'] * D0 + 2 * h * base[off]
                    y11 = [K.of(Fx['s0e1']) * Fx['T0'] * D1[j] + 2 * h * gens[off][j] for j in range(f)]
                    lo, hi = Fx['piece']
                    out.append((y10 - lo, y11))
                    out.append((hi - y10, [-z for z in y11]))
                return out, (D0, D1)
            c1, Dk = cons_of(F, 0)
            c2, _ = cons_of(F2, n1)
            cons = c1 + c2
            if f == 0:
                if all(sg(a_) > 0 for a_, _ in cons):
                    coinc.append((u, v, Dk[0]))
                continue
            # numeric bounding box of {k : a + b·k > 0} (LP), then exact enumeration
            from scipy.optimize import linprog
            A_ub = [[-float(bj) for bj in b_] for a_, b_ in cons]
            b_ub = [float(a_) for a_, b_ in cons]
            box, unbounded = [], False
            for j in range(f):
                for sgn_ in (1, -1):
                    cvec = [0.0] * f
                    cvec[j] = sgn_
                    res = linprog(cvec, A_ub=A_ub, b_ub=b_ub, bounds=[(None, None)] * f, method='highs')
                    if res.status == 3:
                        unbounded = True
                    elif res.status == 2:
                        box = None
                        break
                    else:
                        box.append(sgn_ * res.fun)
                if box is None or unbounded:
                    break
            if box is None:
                continue
            if unbounded:
                infinite.append((u, v))
                continue
            rng = [range(int(math.floor(box[2 * j])) - 3, int(math.ceil(box[2 * j + 1])) + 4) for j in range(f)]   # audit s931 R5
            for kv in itertools.product(*rng):
                if all(sg(a_ + sum((b_[j] * kv[j] for j in range(f)), K.zero)) > 0 for a_, b_ in cons):
                    coinc.append((u, v, Dk[0] + sum((Dk[1][j] * kv[j] for j in range(f)), K.zero)))
    # audit s931 R4: REGULAR ends of an R-exit piece (y1 = 0: the base at B_R; or a 'B' breakpoint of the R map) — the trajectory
    # continues there as the limit of the adjacent pieces, so a lattice point of F with y1 on such an end would be a wall no family
    # records. Solve (y_b − e1 y0 − 2h N1) W = (P + Q·N)·A1 over ℤ² (A1 = s0 e1 T0); any integral solution with D > 0 fails.
    bw_rows, n_bend = [], 0
    for j, F in enumerate(fams):
        if F['transits'] != 2:
            continue
        pcs = G.cap_map('R', F['d1'])
        kk = next(i for i, pp in enumerate(pcs) if (pp['lo'] - F['piece'][0]).iszero() and (pp['hi'] - F['piece'][1]).iszero())
        ends = []
        if kk == 0:
            ends.append(K.zero)
        elif pcs[kk - 1]['bp']['kind'] == 'vertex' and pcs[kk - 1]['bp']['vertex'] == 'B':
            ends.append(F['piece'][0])
        if 'bp' in pcs[kk] and pcs[kk]['bp']['kind'] == 'vertex' and pcs[kk]['bp']['vertex'] == 'B':
            ends.append(F['piece'][1])
        for yb in ends:
            n_bend += 1
            A1 = K.of(F['s0e1']) * F['T0']
            q1, q2 = F['Q']
            sol = solve([-(2 * h) * F['W'] - q1 * A1, -(q2 * A1)], F['P'] * A1 - (yb - F['e'][0] * F['y0']) * F['W'])
            if sol is None:
                continue
            x0, ker = sol
            bad = bool(ker) or (all(x.denominator == 1 for x in x0) and sg(D_of(F, [int(x0[0]), int(x0[1])])) > 0)
            bw_rows.append(dict(fam=j, kernel=len(ker), integral_realised=bad))
            if bad:
                rep['fail'].append(('regular-end wall', j))
    rep['regular_ends'] = {'ends': n_bend, 'consistent': len(bw_rows), 'bad': sum(r['integral_realised'] for r in bw_rows)}
    rep['pair_stats'] = {'none': n_none, 'unique_rational': n_unique, 'line': n_line,
                         'infinite': len(infinite), 'finite_coincidences': len(coinc)}

    def is_reverse(F, F2):
        (vi, d), (wi, d2) = F['slot'], F2['slot']
        return (F['transits'] == F2['transits'] and VLAB[F['target']] == wi and VLAB[F2['target']] == vi
                and (F['arr'] + 2 * N) % (4 * N) == d2 and (F2['arr'] + 2 * N) % (4 * N) == d)
    rep['infinite'] = [dict(u=u, v=v, a=(fams[u]['slot'], fams[u]['target'], fams[u]['transits'], fams[u]['e']),
                            b=(fams[v]['slot'], fams[v]['target'], fams[v]['transits'], fams[v]['e']),
                            reverse_ends=is_reverse(fams[u], fams[v])) for u, v in infinite]
    rep['coincidences'] = [dict(u=u, v=v, D=float(D), D_exact=[str(c) for c in D.c]) for u, v, D in coinc]
    rep['_coinc'] = coinc
    for j, F in enumerate(fams):
        rep['families'].append(dict(j=j, slot=list(F['slot']), X=F['X'], transits=F['transits'], e=list(F['e']),
                                    target=F['target'], arr=F['arr'], W=float(F['W']), P=float(F['P']),
                                    Q=[float(q) for q in F['Q']]))
    return rep, G, E, fams, D_of, realised


if __name__ == '__main__':
    tag = sys.argv[1]
    rep, G, E, fams, D_of, realised = run(tag)
    rep.pop('_coinc', None)
    json.dump(rep, open(f'data/s931_twist_walls_{tag.replace(":", "_").replace(",", "")}.json', 'w'), indent=1, default=str)
    print(tag, 'checks', {k: v for k, v in rep['checks'].items() if k != 'stage0'})
    print('  stage0 kinds', rep['checks']['stage0'])
    print('  families', rep['n_families'], 'pair stats', rep['pair_stats'])
    print('  fail', rep['fail'], ' regular ends', rep['regular_ends'])
    for dg in rep['infinite']:
        if not dg['reverse_ends']:
            print('  INFINITE NON-REVERSE', dg)
    print('  infinite overlaps all reversals:', all(dg['reverse_ends'] for dg in rep['infinite']))
    for c in rep['coincidences'][:30]:
        print('  coincidence', c['u'], c['v'], c['D'])


def crosscheck(tag, T0=Fr(1, 20), T1=Fr(4450, 10000), fn=None):
    """control: every wall-family member with t in [T0, T1] is a breakpoint of the exact t-tracker
    (`probes/s931_trapezoid_exact.py`, its saved JSON), and every tracked breakpoint is a family member (exact equality)."""
    rep, G, E, fams, D_of, realised = run(tag)
    K, sg, h = G.K, G.sg, G.h
    c = G.c_slice
    Dlo, Dhi = K.one / K.of(T1) - c, K.one / K.of(T0) - c
    members = []
    for j, F in enumerate(fams):
        def ok(Nv):
            r, D = realised(F, Nv)
            return r and sg(D - Dlo) >= 0 and sg(Dhi - D) >= 0, D
        if F['transits'] == 1:
            q = F['Q'][0]
            # D = (P + q n)/W  monotone in n: scan a generous range
            n_ = [(Dlo * F['W'] - F['P']) / q, (Dhi * F['W'] - F['P']) / q]
            lo_, hi_ = sorted(int(mpmath.floor(x.mp(30))) for x in n_)
            for n in range(lo_ - 2, hi_ + 3):
                r, D = ok([n])
                if r:
                    members.append((D, j))
        else:
            # y1 = e1 y0 + s0e1 T0 D + 2h N1 ∈ (0, h): N1 range from the D range
            ends = [(-(F['e'][0] * F['y0']) - K.of(F['s0e1']) * F['T0'] * D) / (2 * h) for D in (Dlo, Dhi)]
            lo1, hi1 = sorted(int(mpmath.floor(x.mp(30))) for x in ends)
            for n1 in range(lo1 - 2, hi1 + 3):
                q1, q2 = F['Q']
                n_ = [(D * F['W'] - F['P'] - q1 * n1) / q2 for D in (Dlo, Dhi)]
                lo2, hi2 = sorted(int(mpmath.floor(x.mp(30))) for x in n_)
                for n2 in range(lo2 - 2, hi2 + 3):
                    r, D = ok([n1, n2])
                    if r:
                        members.append((D, j))
    # distinct D values
    distinct = []
    for D, j in sorted(members, key=lambda x: x[0].mp(30)):
        if distinct and (D - distinct[-1][0]).iszero():
            distinct[-1][1].append(j)
        else:
            distinct.append([D, [j]])
    fn = fn if fn else f"data/s931_trapezoid_{tag}_1_20_4450_10000.json" if (T0, T1) == (Fr(1, 20), Fr(4450, 10000)) else f"data/s931_trapezoid_{tag}_{str(T0).replace('/', '_')}_{str(T1).replace('/', '_')}.json"
    tr = json.load(open(fn))
    tracked = []
    for b in tr['bp']:
        t = K.of([Fr(x) for x in b['t_exact']])
        tracked.append(K.one / t - c)
    tracked.sort(key=lambda x: x.mp(30))
    tset = list(tracked)
    unmatched_fam = [float(D) for D, js in distinct if not any((D - x).iszero() for x in tset)]
    unmatched_tr = [float(x) for x in tset if not any((D - x).iszero() for D, js in distinct)]
    per_D = sorted(set(len(js) for D, js in distinct))
    return dict(tag=tag, family_D=len(distinct), tracked=len(tset), unmatched_family=unmatched_fam,
                unmatched_tracked=unmatched_tr, families_per_D=per_D)
