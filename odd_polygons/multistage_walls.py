"""[Deposit] Helper, imported by the scripts of this directory. Paper: "Perpendicular billiards in rational polygons of odd denominator". Research name: probes/s932_multistage_walls.py, copied verbatim except imports (deposit names), paths (data/ and figs/ under odd_polygons/) and this paragraph. Run through run_all.py; the usage lines below are the research repo's. Bracketed IDs and section marks (§sNNN) refer to the project's internal ledgers.

s932 keeper (queue `s931-c`): (M⁺) ⊕ CP at EVERY side length for trapezoids with an ACYCLIC cap transfer graph — the
MULTI-STAGE generalisation of `probes/s931_twist_walls.py` ([POLY-TRAP-ALL]'s proof) to any number of strip transits.

SETTING (§s931 (1), (3), (4)). A trapezoid (a0, a1, N − a1, N − a0)/N, both base angles acute, rescaled: right cap R (corner a1π/N),
strip of length D ∈ (0, ∞), left cap L (corner a0π/N); D is the only shape parameter. A prong = stage-0 cap path, then transits
alternating between the caps, each a fold y_{k} = e_k (y′ + s T D) + 2h n_k, e_k = ±1, n_k ∈ ℤ, through exact cap-map pieces
y′ = c + λ y. By induction y_k = α_k + β_k D + γ_k·n (α, β, γ exact in K = ℚ(cos π/2N)). A WALL FAMILY ends at a vertex height of
the k-th cap (a vertex breakpoint, or the top corner h, reached moving UP): β_k D = (y_t − α_k) − γ_k·n, n ∈ ℤ^k, realised iff
D > 0 and lo_j < y_j < hi_j at every intermediate stage j. ACYCLIC transfer graph ([NCYL-722] (b)) ⟹ finitely many families.
Coincidence of two families: one K-equation in (n, n′) = deg rational equations, solved exactly; the integer solutions inside the
realisation cone are enumerated by recursive LP bounding (relaxed by 1e-6, so a float INFEASIBLE is conservative) ⊕ exact checks
at the leaves. Regular ends of intermediate pieces (y_j = 0, or a `B` breakpoint): solved as an extra equation, as s931 audit R4.

PRE-REGISTERED (s932, before any acyclic family was enumerated).
  Claim under test (T): for (2,3,6,7)/9 and (3,4,5,6)/9 — (i) no stage-0 vertex end; (ii) every wall family has β ≠ 0; (iii) every
  pair of families from different slots either meets nowhere, or overlaps on an infinite set AND is a connection with its own
  reversal; no realised finite meeting of non-reversals; (iv) no regular-end lattice point. Then p ≤ 1 at every D and Theorem P
  ([POLY-THM-P]) gives (M⁺) ⊕ CP at every side length.
  Refuted by: a stage-0 vertex end; a non-reversal infinite overlap; a realised finite meeting (then check the signed relation of
  the two pair lengths before concluding anything); a regular-end point.
  Could come out otherwise: (iii) is a rank question in up to 2·(chain + 1) integer unknowns against deg = φ(4N)/2 = 6 equations —
  at chain 5 the kernel has dimension ≥ 6, nothing forces the realisation cone to miss every lattice point.
  Controls: (R) REGRESSION — on (1,2,5,6)/7, (1,3,4,6)/7 this must reproduce `s931_twist_walls`' family counts 34, 28, 0 realised
  finite meetings, infinite overlaps 25, 20 all reversals; (F) FIRES — the isosceles (1,1,6,6)/7 must report non-reversal
  infinite overlaps (s931: 4); (X) the wall members reproduce the exact t-tracker's breakpoints (`probes/s931_trapezoid_exact.py`)
  on [1/20, t_max) exactly, both ways (`crosscheck`).
PRIOR ART: §s931 (3) ⊕ `probes/s931_twist_walls.py` (two transits, absorbing L); §s931 (4) ⊕ [NCYL-722] (the transfer graph, chains
  3 and 5, float max transits 4 and 6); `rulings.py --grep` *multi-stage*, *acyclic* → [NCYL-722] only.
DATA: `data/s931_twist_survey.json` (float transits, read for the expected maxima); `data_index.py --grep 2367` → none.
Output: data/s932_multistage_<tag>.json.  Usage: python probes/s932_multistage_walls.py N:a0,a1,a2,a3 [crosscheck]
"""
import _setup  # noqa: F401  (deposit: this directory on the path, cwd = odd_polygons/)
import itertools
import json
import math
import sys
import time
from fractions import Fraction as Fr
from math import lcm

import mpmath
from scipy.optimize import linprog

import s6_10_walls as TW
from poly_exact import Aff, Fam


def setup(tag):
    nn, aa = tag.split(':')
    N, a = int(nn), [int(x) for x in aa.split(',')]
    assert a[0] + a[3] == N and a[1] + a[2] == N and 2 * a[1] < N and 2 * a[0] < N
    TW.N = N
    return N, a, TW.Geo(a)


def stage0(G, a, N, slots):
    K, sg, h = G.K, G.sg, G.h
    capR, capL = G.R, G.L
    out = {}
    for (i, d) in slots:
        cosd = sg(G.cos(d))
        if i == 1:
            o = G.trace(capR, d, capR['verts'][1], 0, ('vtx', 1))
        elif i == 0:
            o = G.trace(capL, d, capL['verts'][0], 0, ('vtx', 0))
        elif i == 2:
            if cosd > 0:
                o = G.trace(capR, d, capR['verts'][2], 0, ('vtx', 2))
            elif cosd == 0:
                o = dict(kind='perp', side=0)
            else:
                o = dict(kind='exit', y=Aff(h, K.zero), d=d, length=Aff(K.zero, K.zero))
        else:
            if cosd < 0:
                o = G.trace(capL, d, capL['verts'][2], 0, ('vtx', 2))
            elif cosd == 0:
                o = dict(kind='perp', side=0)
            else:
                o = dict(kind='exit', y=Aff(h, K.zero), d=d, length=Aff(K.zero, K.zero))
        out[(i, d)] = o
    return out


def families(G, st0, maxT=20):
    """DFS over transits. A family: slot, target, arr, transits k, W (=β_k), P (=y_t − α_k), Q (=−γ_k), stages: list of
    (α_j, β_j, γ_j, lo_j, hi_j, lo_regular, hi_regular) for the intermediate stages j < k."""
    K, sg, h = G.K, G.sg, G.h
    fams = []
    top = {'R': 'V2', 'L': 'V3'}

    def rec(slot, al, be, ga, dex, X, stages, path):
        if len(ga) >= maxT:
            raise RuntimeError('transit depth exceeded — transfer graph cyclic?')
        s, T = sg(G.sin(dex)), G.tan_abs(dex)
        for e in (1, -1):
            d1 = dex if e == 1 else TW.flip(dex)
            a1 = al * e
            b1 = (be + K.of(s) * T) * e
            g1 = [g * e for g in ga] + [K.of(2) * h]
            pieces = G.cap_map(X, d1)
            targets = [(h, top[X], TW.up(d1))] if e == 1 else []
            for p in pieces:
                if 'bp' in p and p['bp']['kind'] == 'vertex' and p['bp']['vertex'] != 'B':
                    targets.append((p['hi'], p['bp']['vertex'], p['bp']['arr']))
            for yt, lab, arr in targets:
                fams.append(dict(slot=slot, target=lab, arr=arr, transits=len(g1), W=b1, P=yt - a1, Q=[-g for g in g1],
                                 stages=list(stages), path=path + [(X, d1, e)]))
            for k_, p in enumerate(pieces):
                if p['kind'] != 'exit':
                    continue
                lo_reg = (k_ == 0) or (pieces[k_ - 1]['bp']['kind'] == 'vertex' and pieces[k_ - 1]['bp']['vertex'] == 'B')
                hi_reg = ('bp' in p and p['bp']['kind'] == 'vertex' and p['bp']['vertex'] == 'B')
                c, lam = p['y'].c0, p['y'].c1
                st = stages + [(a1, b1, g1, p['lo'], p['hi'], lo_reg, hi_reg)]
                Y = 'L' if X == 'R' else 'R'
                rec(slot, c + lam * a1, lam * b1, [lam * g for g in g1], p['d'], Y, st, path + [(X, d1, e, k_)])
    for slot, o in st0.items():
        if o['kind'] != 'exit':
            continue
        X = 'R' if TW.mov_right(G, o['d']) else 'L'
        rec(slot, o['y'].c0, K.zero, [], o['d'], X, [], [])
    return fams


class Solver:
    def __init__(self, G):
        self.G, self.K, self.sg = G, G.K, G.sg

    def vec(self, x):
        return list(x.c) + [Fr(0)] * (self.K.deg - len(x.c))

    def solve(self, cols, rhs):
        K = self.K
        m = len(cols)
        A = [[self.vec(cj)[r] for cj in cols] + [self.vec(rhs)[r]] for r in range(K.deg)]
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

    def lattice(self, cols, rhs, Dexpr, conlist):
        """integer n with Σ n_j cols_j = rhs (K), inside the cone. Dexpr = (P, Q, W): D = (P + Q·n[:len Q])/W.
        conlist: list of (off, al, be, ga, lo, hi) meaning lo < al + be D + ga·n[off:off+len ga] < hi (lo/hi None = none).
        Returns ('none'|'finite'|'infinite', points, kernel_dim) — points = list of D values (K) for finite."""
        K, sg = self.K, self.sg
        m = len(cols)
        sol = self.solve(cols, rhs) if cols else ([], [])
        if sol is None:
            return 'none', [], None
        x0, ker = sol
        f = len(ker)
        per = 1
        for i in range(m):
            per = lcm(per, x0[i].denominator, *[kv[i].denominator for kv in ker])
        residues = [r for r in itertools.product(range(per), repeat=f)
                    if all((x0[i] + sum((r[j] * ker[j][i] for j in range(f)), Fr(0))).denominator == 1 for i in range(m))]
        P, Q, W = Dexpr
        pts, infinite = [], False
        for r in residues:
            base = [x0[i] + sum((r[j] * ker[j][i] for j in range(f)), Fr(0)) for i in range(m)]
            gens = [[per * ker[j][i] for j in range(f)] for i in range(m)]
            Dn = (P + sum((q * base[i] for i, q in enumerate(Q)), K.zero)) / W
            Dk = [sum((q * gens[i][j] for i, q in enumerate(Q)), K.zero) / W for j in range(f)]
            cons = [(Dn, Dk)]
            for off, al, be, ga, lo, hi in conlist:
                yn = al + be * Dn + sum((g * base[off + i] for i, g in enumerate(ga)), K.zero)
                yk = [be * Dk[j] + sum((g * gens[off + i][j] for i, g in enumerate(ga)), K.zero) for j in range(f)]
                if lo is not None:
                    cons.append((yn - lo, yk))
                if hi is not None:
                    cons.append((hi - yn, [-z for z in yk]))
            if f == 0:
                if all(sg(a_) > 0 for a_, _ in cons):
                    pts.append(Dn)
                continue
            st, ks = self.enum(cons, f)
            if st == 'infinite':
                infinite = True
                continue
            for kv in ks:
                pts.append(Dn + sum((Dk[j] * kv[j] for j in range(f)), K.zero))
        if infinite:
            return 'infinite', pts, f
        return ('finite' if pts else 'none'), pts, f

    def enum(self, cons, f):
        """integer k ∈ ℤ^f with a + b·k > 0 for all (a, b): recursive LP bounding (relaxed), exact leaf checks."""
        sg = self.sg
        A = [[-float(bj) for bj in b_] for a_, b_ in cons]
        bb = [float(a_) + 1e-6 * (1 + abs(float(a_))) for a_, b_ in cons]
        out = []

        def rec(fixed):
            j = len(fixed)
            if j == f:
                kv = fixed
                if all(sg(a_ + sum((b_[i] * kv[i] for i in range(f)), self.K.zero)) > 0 for a_, b_ in cons):
                    out.append(tuple(kv))
                return 'ok'
            nf = f - j
            Aj = [row[j:] for row in A]
            bj = [bb[r] - sum(A[r][i] * fixed[i] for i in range(j)) for r in range(len(A))]
            lohi = []
            for sgn_ in (1, -1):
                cvec = [0.0] * nf
                cvec[0] = sgn_
                res = linprog(cvec, A_ub=Aj, b_ub=bj, bounds=[(None, None)] * nf, method='highs')
                if res.status == 3:
                    return 'infinite'
                if res.status == 2:
                    return 'ok'
                if res.status != 0:
                    raise RuntimeError(f'LP status {res.status}')
                lohi.append(sgn_ * res.fun)
            for kj in range(int(math.floor(lohi[0])) - 1, int(math.ceil(lohi[1])) + 2):
                if rec(fixed + [kj]) == 'infinite':
                    return 'infinite'
            return 'ok'
        if rec([]) == 'infinite':
            return 'infinite', None
        return 'finite', out


def fam_cons(F, off):
    return [(off, al, be, ga, lo, hi) for (al, be, ga, lo, hi, _, _) in F['stages']]


def is_reverse(F, F2, N):
    VLAB = {'V0': 0, 'V1': 1, 'V2': 2, 'V3': 3}
    (vi, d), (wi, d2) = F['slot'], F2['slot']
    return (F['transits'] == F2['transits'] and VLAB[F['target']] == wi and VLAB[F2['target']] == vi
            and (F['arr'] + 2 * N) % (4 * N) == d2 and (F2['arr'] + 2 * N) % (4 * N) == d)


def run(tag, verbose=True):
    t0 = time.time()
    N, a, G = setup(tag)
    K, sg, h = G.K, G.sg, G.h
    E = Fam.closing(a, N, [(1, 0), (0, 1)])
    G.c_slice = -E.L[2].c1
    slots = E.slots()
    rep = {'tag': tag, 'fail': []}
    st0 = stage0(G, a, N, slots)
    for k, o in st0.items():
        if o['kind'] == 'vertex':
            rep['fail'].append(('stage-0 vertex end', k, o['vertex']))
    fams = families(G, st0)
    rep['n_families'] = len(fams)
    rep['transit_hist'] = dict(sorted(collections_counter(F['transits'] for F in fams).items()))
    if any(F['W'].iszero() for F in fams):
        rep['fail'].append('W = 0 family')
    S = Solver(G)
    stats = collections_counter()
    infinite, coinc, kdims = [], [], collections_counter()
    for u, v in itertools.combinations(range(len(fams)), 2):
        F, F2 = fams[u], fams[v]
        if F['slot'] == F2['slot']:
            continue
        n1 = len(F['Q'])
        cols = [q * F2['W'] for q in F['Q']] + [-(q * F['W']) for q in F2['Q']]
        rhs = F2['P'] * F['W'] - F['P'] * F2['W']
        st, pts, fdim = S.lattice(cols, rhs, (F['P'], F['Q'], F['W']), fam_cons(F, 0) + fam_cons(F2, n1))
        stats[st] += 1
        kdims[fdim] += 1
        if st == 'infinite':
            infinite.append((u, v))
        for D in pts:
            coinc.append((u, v, D))
    # regular ends of intermediate pieces
    n_end, bad_end = 0, []
    for j, F in enumerate(fams):
        for si, (al, be, ga, lo, hi, lo_reg, hi_reg) in enumerate(F['stages']):
            for yb, reg in ((lo, lo_reg), (hi, hi_reg)):
                if not reg:
                    continue
                n_end += 1
                # wall: W D = P + Q·n ; end: al + be D + ga·n[:len ga] = yb  ⟹  (yb − al − ga·n)·W = be (P + Q·n)
                kq = len(F['Q'])
                gpad = list(ga) + [K.zero] * (kq - len(ga))
                cols = [-(g * F['W']) - be * q for g, q in zip(gpad, F['Q'])]
                rhs = be * F['P'] - (yb - al) * F['W']
                others = [(0, a2, b2, g2, l2, h2) for k2, (a2, b2, g2, l2, h2, _, _) in enumerate(F['stages']) if k2 != si]
                st, pts, fdim = S.lattice(cols, rhs, (F['P'], F['Q'], F['W']), others)
                if st != 'none':
                    bad_end.append(dict(fam=j, stage=si, status=st, n=len(pts)))
    rep['regular_ends'] = {'ends': n_end, 'bad': bad_end}
    if bad_end:
        rep['fail'].append(('regular-end lattice points', len(bad_end)))
    rep['pair_stats'] = dict(stats)
    rep['kernel_dims'] = {str(k): v for k, v in kdims.items()}
    rep['infinite'] = [dict(u=u, v=v, a=(list(fams[u]['slot']), fams[u]['target'], fams[u]['transits']),
                            b=(list(fams[v]['slot']), fams[v]['target'], fams[v]['transits']),
                            reverse=is_reverse(fams[u], fams[v], N)) for u, v in infinite]
    rep['n_infinite'] = len(infinite)
    rep['n_infinite_nonreverse'] = sum(1 for x in rep['infinite'] if not x['reverse'])
    rep['coincidences'] = [dict(u=u, v=v, D=float(D), D_exact=[str(c) for c in D.c],
                                reverse=is_reverse(fams[u], fams[v], N)) for u, v, D in coinc]
    rep['n_finite_meetings'] = len(coinc)
    rep['families'] = [dict(j=j, slot=list(F['slot']), target=F['target'], arr=F['arr'], transits=F['transits'],
                            path=[list(map(str, p)) for p in F['path']], W=float(F['W'])) for j, F in enumerate(fams)]
    rep['seconds'] = round(time.time() - t0, 1)
    rep['_fams'], rep['_G'] = fams, G
    if verbose:
        print(f"{tag}: families {len(fams)} transits {rep['transit_hist']} pairs {rep['pair_stats']} kernel dims "
              f"{rep['kernel_dims']} infinite {len(infinite)} (non-reverse {rep['n_infinite_nonreverse']}) finite meetings "
              f"{len(coinc)} regular ends {n_end} bad {len(bad_end)} fail {rep['fail'][:4]} {rep['seconds']}s", flush=True)
    return rep


def collections_counter(it=()):
    import collections
    return collections.Counter(it)


def members(rep, Dlo, Dhi):
    """every realised member D ∈ [Dlo, Dhi] of every family (enumerate the lattice of the wall equation alone)."""
    fams, G = rep['_fams'], rep['_G']
    K, sg = G.K, G.sg
    S = Solver(G)
    out = []
    for j, F in enumerate(fams):
        k = len(F['Q'])
        # n free in ℤ^k, constraints: stages, D ∈ [Dlo, Dhi]
        cons = fam_cons(F, 0)
        # D = (P + Q n)/W as a pseudo-stage with al = 0, be = 1, ga = 0
        cons = cons + [(0, K.zero, K.one, [], Dlo - K.of(Fr(1, 10**30)), Dhi + K.of(Fr(1, 10**30)))]
        # identity system: no equation, kernel = all of ℤ^k
        x0, ker = [Fr(0)] * k, [[Fr(int(i == j_)) for i in range(k)] for j_ in range(k)]
        S_ = S

        class _S(Solver):
            def solve(self, cols, rhs):
                return x0, ker
        st, pts, _ = _S(G).lattice([K.zero] * k, K.zero, (F['P'], F['Q'], F['W']), cons)
        assert st != 'infinite'
        out += [(D, j) for D in pts]
    return out


def crosscheck(tag, fn):
    from fractions import Fraction as Fr_
    N, a, G = setup(tag)
    E = Fam.closing(a, N, [(1, 0), (0, 1)])
    G.c_slice = -E.L[2].c1
    rep = {'_fams': families(G, stage0(G, a, N, E.slots())), '_G': G}
    K, sg = G.K, G.sg
    c = G.c_slice
    tr = json.load(open(fn))
    T0 = Fr_(tr['range'][0]) if 'range' in tr else None
    tracked = sorted((K.one / K.of([Fr_(x) for x in b['t_exact']]) - c for b in tr['bp']), key=lambda x: x.mp(30))
    Tlo, Thi = [K.of(Fr_(x)) if not hasattr(x, 'c') else x for x in (tr['T0'], tr['T1'])]
    Dlo, Dhi = K.one / Thi - c, K.one / Tlo - c
    mem = members(rep, Dlo, Dhi)
    distinct = []
    for D, j in sorted(mem, key=lambda x: x[0].mp(30)):
        if distinct and (D - distinct[-1][0]).iszero():
            distinct[-1][1].append(j)
        else:
            distinct.append([D, [j]])
    unf = [float(D) for D, js in distinct if not any((D - x).iszero() for x in tracked)]
    unt = [float(x) for x in tracked if not any((D - x).iszero() for D, js in distinct)]
    return dict(tag=tag, family_D=len(distinct), tracked=len(tracked), unmatched_family=unf, unmatched_tracked=unt,
                families_per_D=sorted(set(len(js) for D, js in distinct)))


if __name__ == '__main__':
    tag = sys.argv[1]
    if len(sys.argv) > 2 and sys.argv[2] == 'crosscheck':
        print(crosscheck(tag, sys.argv[3]))
        sys.exit()
    rep = run(tag)
    rep.pop('_fams')
    rep.pop('_G')
    json.dump(rep, open(f'data/s932_multistage_{tag.replace(":", "_").replace(",", "")}.json', 'w'), indent=1, default=str)


# ---------------------------------------------------------------------------------------------------------------------
# UNBOUNDED pairs (s932): the LP relaxation is unbounded iff the tube direction v lies in the real span of the kernel —
# the two walls are PARALLEL in the twist torus. Exact decision: R_v = smallest rational subspace containing v (in kernel
# coordinates k), integer coordinates z = (z_b, z_a) with z_a spanning ℤ^f ∩ R_v (column HNF of the functionals vanishing on
# R_v); the projection of the cone onto z_b is bounded (enumerated, recursive LP); each z_b-slice is OPEN and recedes along v,
# and v is dense in the torus R_v/(ℤ^f ∩ R_v) (minimality of R_v) ⟹ a nonempty slice holds INFINITELY many lattice points,
# an empty one none. So an unbounded pair meets nowhere or on an infinite set — never finitely often.

def dn_of(F):
    """exact tube direction dn/dD ∈ K^k: every stage height held constant, the wall defines D."""
    dn = []
    for (al, be, ga, lo, hi, _, _) in F['stages']:
        acc = be + sum((g * x for g, x in zip(ga[:-1], dn)), ga[-1] * 0)
        dn.append(-acc / ga[-1])
    Q, W = F['Q'], F['W']
    acc = W - sum((q * x for q, x in zip(Q[:-1], dn)), W * 0)
    dn.append(acc / Q[-1])
    return dn


def col_hnf(Wm):
    """integer matrix Wm (r × f, list of lists of int): unimodular V (f × f) with Wm·V = [H | 0]; returns (V, rank)."""
    r, f = len(Wm), len(Wm[0]) if Wm else 0
    A = [row[:] for row in Wm]
    V = [[int(i == j) for j in range(f)] for i in range(f)]

    def colop(i, j, q):          # col_j -= q col_i
        for row in A:
            row[j] -= q * row[i]
        for row in V:
            row[j] -= q * row[i]

    def swap(i, j):
        for row in A:
            row[i], row[j] = row[j], row[i]
        for row in V:
            row[i], row[j] = row[j], row[i]
    piv = 0
    for i in range(r):
        if piv >= f:
            break
        while True:
            nz = [j for j in range(piv, f) if A[i][j] != 0]
            if not nz:
                break
            jm = min(nz, key=lambda j: abs(A[i][j]))
            swap(piv, jm)
            done = True
            for j in range(piv + 1, f):
                if A[i][j] != 0:
                    colop(piv, j, A[i][j] // A[i][piv])
                    if A[i][j] != 0:
                        done = False
            if done:
                break
        if any(A[i][j] != 0 for j in range(piv, f)):
            piv += 1
    return V, piv


def refine_unbounded(S, F, F2):
    """exact verdict for a pair whose LP relaxation is unbounded: ('none' | 'infinite', detail)."""
    K, sg = S.K, S.sg
    n1 = len(F['Q'])
    cols = [q * F2['W'] for q in F['Q']] + [-(q * F['W']) for q in F2['Q']]
    rhs = F2['P'] * F['W'] - F['P'] * F2['W']
    sol = S.solve(cols, rhs)
    if sol is None:
        return 'none', 'inconsistent'
    x0, ker = sol
    f, m = len(ker), len(x0)
    v = dn_of(F) + dn_of(F2)                                  # ∈ K^m, per unit D
    per = 1
    for i in range(m):
        per = lcm(per, x0[i].denominator, *[kv[i].denominator for kv in ker])
    residues = [r for r in itertools.product(range(per), repeat=f)
                if all((x0[i] + sum((r[j] * ker[j][i] for j in range(f)), Fr(0))).denominator == 1 for i in range(m))]
    gens = [[per * ker[j][i] for j in range(f)] for i in range(m)]
    # λ_b ∈ Q^f with gens·λ_b = v_b (v_b = the Q-components of v): least-squares exact via the pivot rows of gens
    vb = [[S.vec(v[i])[b] for i in range(m)] for b in range(K.deg)]
    rows = []
    M_ = []
    for i in range(m):                                        # greedy independent rows of gens
        cand = M_ + [gens[i]]
        if rank_q(cand) > len(M_):
            M_.append(gens[i])
            rows.append(i)
        if len(M_) == f:
            break
    lam = []
    for b in range(K.deg):
        x = solve_q([gens[i] for i in rows], [vb[b][i] for i in rows])
        assert all(sum((gens[i][j] * x[j] for j in range(f)), Fr(0)) == vb[b][i] for i in range(m)), 'v not in kernel span'
        lam.append(x)
    # functionals vanishing on R_v = span(lam): integer basis of the left null space
    ann = nullspace_q(lam, f)                                 # list of rational vectors w with w·λ_b = 0
    Wint = [[int(x * lcm(*[y.denominator for y in w])) for x in w] for w in ann]
    V, rk = col_hnf(Wint) if Wint else ([[int(i == j) for j in range(f)] for i in range(f)], 0)
    s = f - rk
    out = []
    for r in residues:
        base = [x0[i] + sum((r[j] * ker[j][i] for j in range(f)), Fr(0)) for i in range(m)]
        # k = V z
        G2 = [[sum((gens[i][j] * V[j][c] for j in range(f)), Fr(0)) for c in range(f)] for i in range(m)]
        P_, Q_, W_ = F['P'], F['Q'], F['W']
        Dn = (P_ + sum((q * base[i] for i, q in enumerate(Q_)), K.zero)) / W_
        Dk = [sum((q * G2[i][c] for i, q in enumerate(Q_)), K.zero) / W_ for c in range(f)]
        cons = [(Dn, Dk)]
        for off, al, be, ga, lo, hi in fam_cons(F, 0) + fam_cons(F2, n1):
            yn = al + be * Dn + sum((g * base[off + i] for i, g in enumerate(ga)), K.zero)
            yk = [be * Dk[c] + sum((g * G2[off + i][c] for i, g in enumerate(ga)), K.zero) for c in range(f)]
            cons.append((yn - lo, yk))
            cons.append((hi - yn, [-z for z in yk]))
        A = [[-float(bj) for bj in b_] for a_, b_ in cons]
        bb = [float(a_) + 1e-6 * (1 + abs(float(a_))) for a_, b_ in cons]
        zbs = []

        def rec(fixed):
            j = len(fixed)
            if j == rk:
                zbs.append(list(fixed))
                return
            nf = f - j
            Aj = [row[j:] for row in A]
            bj = [bb[q] - sum(A[q][i] * fixed[i] for i in range(j)) for q in range(len(A))]
            lohi = []
            for sgn_ in (1, -1):
                cvec = [0.0] * nf
                cvec[0] = sgn_
                res = linprog(cvec, A_ub=Aj, b_ub=bj, bounds=[(None, None)] * nf, method='highs')
                if res.status == 2:
                    return
                assert res.status == 0, f'z_b projection not bounded (status {res.status})'
                lohi.append(sgn_ * res.fun)
            for kj in range(int(math.floor(lohi[0])) - 1, int(math.ceil(lohi[1])) + 2):
                rec(fixed + [kj])
        rec([])
        for zb in zbs:
            # slice feasibility (strict), Chebyshev centre in z_a, then an EXACT rational interior point
            Aa = [row[rk:] for row in A]
            ba = [float(a_) - sum(-float(b_[i]) * zb[i] for i in range(rk)) for a_, b_ in cons]
            nrm = [math.sqrt(sum(x * x for x in row)) or 1.0 for row in Aa]
            res = linprog([0.0] * s + [-1.0], A_ub=[row + [nr] for row, nr in zip(Aa, nrm)], b_ub=ba,
                          bounds=[(None, None)] * s + [(0, 1)], method='highs')
            if res.status != 0 or res.x[-1] <= 1e-9:
                # conservative guard: feasible once relaxed ⟹ a slice too thin to decide in floats — reported as MARGINAL
                bar = [x + 1e-6 * (1 + abs(x)) for x in ba]
                r3 = linprog([0.0] * s, A_ub=Aa, b_ub=bar, bounds=[(None, None)] * s, method='highs') if s else None
                if (r3 is not None and r3.status == 0) or (s == 0 and all(x > -1e-6 for x in bar)):
                    # EXACT emptiness certificate: a constraint with b = 0, a ≤ 0, or two OPPOSITE constraints
                    # b_j = −μ b_i (μ > 0) with a_j + μ a_i ≤ 0 (a zero-width slab) — exact in K
                    ex = [(a_ + sum((b_[c] * zb[c] for c in range(rk)), K.zero), b_[rk:]) for a_, b_ in cons]
                    if not empty_certificate(ex, sg):
                        out.append(dict(residue=list(r), zb=zb, marginal=True))
                continue
            za = [Fr(x).limit_denominator(10**9) for x in res.x[:s]]
            z = [Fr(x) for x in zb] + za
            if all(sg(a_ + sum((b_[c] * z[c] for c in range(f)), K.zero)) > 0 for a_, b_ in cons):
                out.append(dict(residue=list(r), zb=zb, radius=float(res.x[-1])))
            else:
                out.append(dict(residue=list(r), zb=zb, radius=float(res.x[-1]), exact_check='FAILED'))
    if any(o.get('marginal') for o in out):
        return 'marginal', dict(f=f, s=s, slices=out)
    return ('infinite' if out else 'none'), dict(f=f, s=s, slices=out)


def rank_q(rows):
    return len(_rref(rows)[1])


def _rref(rows):
    A = [[Fr(x) for x in r] for r in rows]
    piv, row = [], 0
    ncol = len(A[0]) if A else 0
    for col in range(ncol):
        pr = next((r for r in range(row, len(A)) if A[r][col] != 0), None)
        if pr is None:
            continue
        A[row], A[pr] = A[pr], A[row]
        pv = A[row][col]
        A[row] = [x / pv for x in A[row]]
        for r in range(len(A)):
            if r != row and A[r][col] != 0:
                fct = A[r][col]
                A[r] = [x - fct * y for x, y in zip(A[r], A[row])]
        piv.append(col)
        row += 1
    return A, piv


def solve_q(M, b):
    """square nonsingular rational solve."""
    n = len(M)
    A, piv = _rref([list(M[i]) + [b[i]] for i in range(n)])
    assert piv == list(range(n))
    return [A[i][n] for i in range(n)]


def nullspace_q(rows, f):
    """rational basis of {w ∈ Q^f : w·row = 0 for every row}."""
    if not rows:
        return [[Fr(int(i == j)) for i in range(f)] for j in range(f)]
    A, piv = _rref(rows)
    out = []
    for j in range(f):
        if j in piv:
            continue
        w = [Fr(0)] * f
        w[j] = Fr(1)
        for r, col in enumerate(piv):
            w[col] = -A[r][j]
        out.append(w)
    return out


# ---------------------------------------------------------------------------------------------------------------------
# THE PROOF DRIVER (s932). D ≤ D_tr is covered EXACTLY by the t-tracker (`probes/s931_trapezoid_exact.py` on [1/20, t_max):
# p at every breakpoint, the signed relation checked); the walls cover D > D_tr. One family per connection (start vertex <
# end vertex; loops keep both orientations — their reversal pairs are recognised by is_reverse ⊕ an infinite verdict).
# Per pair: exact solve → LP max D on the relaxed cone (relaxed by 1e-6, so the bound is conservative) → bounded ≤ D_tr:
# tracker's; bounded > D_tr: exact enumeration of the lattice points with D > D_tr; unbounded: refine_unbounded.

VIDX = {'V0': 0, 'V1': 1, 'V2': 2, 'V3': 3}


def keep_orientation(F):
    return F['slot'][0] <= VIDX[F['target']]


def pair_verdict(S, F, F2, Dtr, Dlo=None):
    K, sg = S.K, S.sg
    n1 = len(F['Q'])
    cols = [q * F2['W'] for q in F['Q']] + [-(q * F['W']) for q in F2['Q']]
    rhs = F2['P'] * F['W'] - F['P'] * F2['W']
    sol = S.solve(cols, rhs)
    if sol is None:
        return 'none', None
    x0, ker = sol
    f = len(ker)
    Dn = (F['P'] + sum((q * x0[i] for i, q in enumerate(F['Q'])), K.zero)) / F['W']
    Dk = [sum((q * ker[j][i] for i, q in enumerate(F['Q'])), K.zero) / F['W'] for j in range(f)]
    cons = [(Dn, Dk)]
    for off, al, be, ga, lo, hi in fam_cons(F, 0) + fam_cons(F2, n1):
        yn = al + be * Dn + sum((g * x0[off + i] for i, g in enumerate(ga)), K.zero)
        yk = [be * Dk[j] + sum((g * ker[j][off + i] for i, g in enumerate(ga)), K.zero) for j in range(f)]
        cons.append((yn - lo, yk))
        cons.append((hi - yn, [-z for z in yk]))
    if f == 0:
        if all(x.denominator == 1 for x in x0) and all(sg(a_) > 0 for a_, _ in cons):
            if Dlo is not None and sg(Dn - Dlo) < 0:
                return 'finite_thin', [Dn]
            return ('finite_low' if sg(Dn - Dtr) <= 0 else 'finite_high'), [Dn]
        return 'none', None
    A = [[-float(b) for b in bb] for aa, bb in cons]
    bub = [float(aa) + 1e-6 * (1 + abs(float(aa))) for aa, bb in cons]
    r = linprog([-float(x) for x in Dk], A_ub=A, b_ub=bub, bounds=[(None, None)] * f, method='highs')
    if r.status == 2:
        return 'none', None
    if r.status == 3:
        st, det = refine_unbounded(S, F, F2)
        return {'infinite': 'infinite', 'none': 'none_parallel', 'marginal': 'MARGINAL'}[st], det
    assert r.status == 0, r.status
    dmax = float(Dn) - r.fun
    if Dlo is not None:                     # the thin end D < D_lo, outside the tracker's range: exact enumeration
        r2 = linprog([float(x) for x in Dk], A_ub=A, b_ub=bub, bounds=[(None, None)] * f, method='highs')
        if r2.status == 0 and float(Dn) + r2.fun < float(Dlo) + 1e-6:
            conl = fam_cons(F, 0) + fam_cons(F2, n1) + [(0, K.zero, K.one, [], None, Dlo)]
            st, pts, _ = S.lattice(cols, rhs, (F['P'], F['Q'], F['W']), conl)
            assert st != 'infinite'
            if pts:
                return 'finite_thin', pts
    if dmax <= float(Dtr) - 1e-6:
        return 'bounded_low', dmax
    # exact enumeration of lattice points with D > D_tr
    conl = fam_cons(F, 0) + fam_cons(F2, n1) + [(0, K.zero, K.one, [], Dtr, None)]
    st, pts, _ = S.lattice(cols, rhs, (F['P'], F['Q'], F['W']), conl)
    assert st != 'infinite'
    return ('finite_high' if pts else 'bounded_high_empty'), [p for p in pts] if pts else dmax


def run2(tag, Dtr, shard=(0, 1), out=None, every=2000, Dlo=None):
    import collections
    import os
    t0 = time.time()
    N, a, G = setup(tag)
    K, sg = G.K, G.sg
    E = Fam.closing(a, N, [(1, 0), (0, 1)])
    G.c_slice = -E.L[2].c1
    st0 = stage0(G, a, N, E.slots())
    fail = [('stage-0 vertex end', k) for k, o in st0.items() if o['kind'] == 'vertex']
    fams = families(G, st0)
    keep = [j for j, F in enumerate(fams) if keep_orientation(F)]
    S = Solver(G)
    Dtr = K.of(Dtr) if not hasattr(Dtr, 'c') else Dtr
    Dlo = None if Dlo is None else (K.of(Dlo) if not hasattr(Dlo, 'c') else Dlo)
    pairs = [(u, v) for u, v in itertools.combinations(keep, 2) if fams[u]['slot'] != fams[v]['slot']]
    i0, ns = shard
    mine = pairs[i0::ns]
    state = {'done': 0, 'counts': {}, 'events': []}
    if out and os.path.exists(out):
        state = json.load(open(out))
    cnt = collections.Counter(state['counts'])
    for idx in range(state['done'], len(mine)):
        u, v = mine[idx]
        st, det = pair_verdict(S, fams[u], fams[v], Dtr, Dlo)
        if st == 'infinite':
            st = 'infinite_reverse' if is_reverse(fams[u], fams[v], N) else 'INFINITE_NONREVERSE'
        cnt[st] += 1
        if st in ('INFINITE_NONREVERSE', 'finite_high', 'finite_low', 'finite_thin', 'MARGINAL'):
            state['events'].append(dict(u=u, v=v, st=st, D=[float(x) for x in det] if isinstance(det, list) else None,
                                        D_exact=[[str(c) for c in x.c] for x in det] if isinstance(det, list) else None))
        if out and (idx + 1) % every == 0:
            state.update(done=idx + 1, counts=dict(cnt))
            json.dump(state, open(out + '.tmp', 'w'))
            os.replace(out + '.tmp', out)
    state.update(done=len(mine), counts=dict(cnt), finished=True, n_families=len(fams), n_kept=len(keep),
                 n_pairs=len(mine), Dtr=float(Dtr), Dlo=None if Dlo is None else float(Dlo), fail=fail, seconds=round(time.time() - t0, 1))
    if out:
        json.dump(state, open(out, 'w'), indent=1)
    return state, fams, G


def member_census(tag, Dlo, Dhi):
    """every realised member D ∈ (Dlo, Dhi) of every family (ALL families, both orientations), grouped by exact D:
    p(D) = (#families at D)/2 (each connection is found from both ends). Used for the THIN END (0, D_lo] that the t-tracker
    does not reach, and for control (X) on the tracker's range."""
    N, a, G = setup(tag)
    K = G.K
    E = Fam.closing(a, N, [(1, 0), (0, 1)])
    G.c_slice = -E.L[2].c1
    fams = families(G, stage0(G, a, N, E.slots()))
    rep = {'_fams': fams, '_G': G}
    Dlo = K.of(Dlo) if not hasattr(Dlo, 'c') else Dlo
    Dhi = K.of(Dhi) if not hasattr(Dhi, 'c') else Dhi
    mem = members(rep, Dlo, Dhi)
    groups = []
    for D, j in sorted(mem, key=lambda x: x[0].mp(40)):
        if groups and (D - groups[-1][0]).iszero():
            groups[-1][1].append(j)
        else:
            groups.append([D, [j]])
    return groups, fams, G


def regular_ends(tag):
    """audit-R4 analogue for any number of transits: a lattice point of family F with an intermediate height y_j EXACTLY on a
    regular end of its piece (y_j = 0, or a `B` breakpoint) is a wall no family records (both adjacent families exclude it).
    Solve wall ⊕ (y_j = y_b) over ℤ inside the remaining cone. Unbounded relaxations are decided by the same parallel-wall
    argument: the extra equation is one more K-linear condition, so 'infinite' here = the lattice-point set of an unbounded
    cone, decided by refine-style slicing — reported, and must be 'none' for the proof."""
    N, a, G = setup(tag)
    K = G.K
    E = Fam.closing(a, N, [(1, 0), (0, 1)])
    fams = families(G, stage0(G, a, N, E.slots()))
    assert all(not F['W'].iszero() for F in fams), 'W = 0 family'
    S = Solver(G)
    import collections
    cnt, bad = collections.Counter(), []
    for j, F in enumerate(fams):
        for si, (al, be, ga, lo, hi, lo_reg, hi_reg) in enumerate(F['stages']):
            for yb, reg in ((lo, lo_reg), (hi, hi_reg)):
                if not reg:
                    continue
                kq = len(F['Q'])
                gpad = list(ga) + [K.zero] * (kq - len(ga))
                cols = [-(g * F['W']) - be * q for g, q in zip(gpad, F['Q'])]
                rhs = be * F['P'] - (yb - al) * F['W']
                others = [(0, a2, b2, g2, l2, h2) for k2, (a2, b2, g2, l2, h2, _, _) in enumerate(F['stages']) if k2 != si]
                st, pts, fdim = S.lattice(cols, rhs, (F['P'], F['Q'], F['W']), others)
                cnt[(st, fdim)] += 1
                if st != 'none':
                    bad.append(dict(fam=j, stage=si, status=st, kernel=fdim, n=len(pts), D=[float(x) for x in pts][:5]))
    return dict(n_families=len(fams), verdicts={str(k): v for k, v in cnt.items()}, bad=bad)


def empty_certificate(ex, sg):
    """exact: is {z : a + b·z > 0 ∀(a, b) ∈ ex} empty by a one- or two-constraint certificate?"""
    for a_, b_ in ex:
        if all(x.iszero() for x in b_) and sg(a_) <= 0:
            return True
    for i in range(len(ex)):
        ai, bi = ex[i]
        t = next((q for q, x in enumerate(bi) if not x.iszero()), None)
        if t is None:
            continue
        for j in range(i + 1, len(ex)):
            aj, bj = ex[j]
            if bj[t].iszero():
                continue
            mu = -bj[t] / bi[t]
            if sg(mu) <= 0:
                continue
            if all((bj[q] + mu * bi[q]).iszero() for q in range(len(bi))) and sg(aj + mu * ai) <= 0:
                return True
    return False
