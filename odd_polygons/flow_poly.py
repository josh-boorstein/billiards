"""[Deposit] Helper, imported by the scripts of this directory. Paper: "Perpendicular billiards in rational polygons of odd denominator". Research name: engine/flow_poly.py, copied verbatim except imports (deposit names), paths (data/ and figs/ under odd_polygons/) and this paragraph. Run through run_all.py; the usage lines below are the research repo's. Bracketed IDs and section marks (§sNNN) refer to the project's internal ledgers.

flow_poly — `flow_decomp`'s exact flow decomposition, for a CONVEX rational POLYGON with cyclotomic vertices (s933).

`engine/flow_decomp.py` decides CP / non-CP of a direction on a rational TRIANGLE exactly (intervalxt's induction ⊕
Boshernitzan, every component certified, Kac balance asserted). Only its geometry is triangle-specific: `Triangle`
builds the 2d copies D_d·Δ in ℤ[ζ_{4d}] and `_build_iet` reads entry/exit sides per copy using the third vertex. This
module supplies the same interface for a convex k-gon whose side lengths are given EXACTLY as monomial sums in ζ_{4d}
(integer coefficients — rescale rational lengths first), and an entry/exit test by orientation (a convex CCW polygon
has its interior to the LEFT of every side; a reflected copy is CW). Everything downstream (`Decomposer`, the
certificates, Kac) is `flow_decomp`'s, unchanged.

  P = PolyT(a, d, L)        a = interior angles in units π/d (Σ = (k − 2)d, gcd(a, d) = 1), L = exact side lengths
                            (monomial dicts {exponent mod 4d: int}); closure asserted EXACTLY.
  decompose_poly(P, j, cap) as `flow_decomp.decompose`; j = direction exponent (ω = ζ^j) or a monomial dict.
  perp_dirs_poly(P)         j ⊥ side s, per side.
  PolyT(a, d, L, M)         (s934) field multiplier: angles in π/d, lengths in ζ_{4dM} (e.g. rational-diagonal shapes on a
                            finer grid π/(2dM)); `mlift(u, M, 4d)` lifts a ζ_{4d} length. Gate: the M = 2 lift reproduces M = 1.
  sin_mono(x, d)            2i·sin(xπ/d)·(−i) = 2 sin(xπ/d) as a monomial dict (exact).
  trapezoid_cut(a0, a1, d, p, q)   the triangle (a0, a1, d − a0 − a1)/d cut parallel to its side V0V1, keeping the
                            fraction μ = p/q of each leg (μ → 1 is the triangle): angles (a0, a1, d − a1, d − a0),
                            lengths 2q·(sin a2, μ sin a0, (1 − μ) sin a2, μ sin a1).
  mplus_poly(P, j)          (s936) the exact (M⁺) verdict: every component a cylinder AND every cylinder crosses a side copy
                            PERPENDICULAR to ω (`perp_keys`); a cylinder crossing none is FREE. `crossings_poly` = per-component
                            crossing vectors (`flow_decomp.crossings` for a PolyT). ⚠ `free_mass` is the TRANSVERSAL mass Σhλ/L
                            (as `decompose`'s 'mass'), NOT the area fraction 1 − R of §s928 (rhombus: 0.309 vs φ/5 = 0.324).
  marked_components(T, j, verts)   (s938) EXACT POLE MEMBERSHIP: for each marked vertex V_i and each surface point over it
                            (`vertex_points`: copies sharing the corner), the components whose closure holds it — its
                            one-sided germs (`marked_germs`, convex only) tracked through the induction by `MarkDecomposer`
                            (a germ in a removed left part is pushed by the pre-step map: same orbit). One component = the
                            point is interior; two = it sits on their common boundary. GATES (exact, every call): each
                            germ's T-image ends in the same component; a cylinder germ returns after exactly the period.
NON-CONVEX (s936): a reflex angle (d < a < 2d) is allowed. The convex IET pairs every entry side with every exit side whose
x-range overlaps (one segment per flow line); a non-convex copy is paired by a SWEEP (`_pair_sweep`): break at every vertex
x, order the covering sides by exact height (`_above`, orientation tests only), pair consecutive entry/exit sides — the
first exit along ω. Alternation and the two-sided tiling of [0, L) are asserted, so a wrong pairing cannot pass silently.
GATE (`_selftest`): a triangle through PolyT reproduces `flow_decomp.decompose` on T₀ = (7,8,15)/30 (every perp class)
and on (7,7,16)/30 (counts ⊕ masses identical); the sweep reproduces the convex pieces IDENTICALLY (28 convex rows,
s936); `mplus_poly` reproduces §s927 (2) (even control N ≤ 24, no a_i | N: 55/87 direction rows free; odd N 5–29: 356/356
M⁺ — `data/scratch/s936_g2.py`, the selftest runs N ≤ 16); free cylinders on the rhombus (3,2,3,2)/5, none at 1000:2001;
the non-convex mirror-free octagon (134,30,67,15,82,30,41,15)/69 non-CP, as [POLY-FOLD-ODD] proves.
"""
import _setup  # noqa: F401  (deposit: this directory on the path, cwd = odd_polygons/)
import functools
import math

import flow_decomp as FD
from flow_decomp import (Cyc, CycNP, Comp, Decomposer, madd, mneg, mshift, mconj, mmul)


def mlift(u, M, N):
    """a monomial dict in ζ_N re-expressed in ζ_{NM} (exponents × M)"""
    return {(e * M) % (N * M): c for e, c in u.items()}


def sin_mono(x, d):
    """2·sin(xπ/d) = −i(ζ^{2x} − ζ^{−2x}) with ζ = ζ_{4d}, −i = ζ^{3d}"""
    N = 4 * d
    return {(2 * x + 3 * d) % N: 1, (-2 * x + 3 * d) % N: -1} if (4 * x) % N else {}


def sin_half(y, d):
    """2·sin(yπ/(2d)) — half-integer angles (y odd allowed)"""
    N = 4 * d
    return {(y + 3 * d) % N: 1, (-y + 3 * d) % N: -1} if (2 * y) % N else {}


def mscale(u, n):
    return {e: n * c for e, c in u.items() if n * c}


class PolyT:
    """duck-typed `flow_decomp.Triangle` for a convex rational polygon"""

    def __init__(self, a, d, L, M=1):
        k = len(a)
        assert sum(a) == (k - 2) * d and math.gcd(math.gcd(*a), d) == 1
        assert all(0 < x < 2 * d for x in a), 'interior angles in (0, 2π)'
        self.convex = all(x <= d for x in a)                  # s936: a reflex vertex → the sweep pairing in `_build_iet_poly`
        self.abc, self.a, self.k = tuple(a), list(a), k
        self.d = d
        self.M = M                                            # s934: lengths in ζ_{4dM}, angles still in π/d
        self.N = N = 4 * d * M
        self.F = Cyc(N)
        # side s = V_s V_{s+1}, line direction m_s (units π/d): m_0 = 0, m_{s+1} = m_s + (d − a_{s+1})
        m = [0]
        for s in range(1, k):
            m.append(m[-1] + d - a[s])
        self.m = [x % d for x in m]
        self.mfull = [x % (2 * d) for x in m]                 # directed (mod 2π)
        V = [{}]
        for s in range(k - 1):
            V.append(madd(V[-1], mshift(L[s], 2 * M * self.mfull[s], N)))
        close = madd(V[-1], mshift(L[k - 1], 2 * M * self.mfull[k - 1], N))
        assert not any(self.red(close)), 'polygon does not close (exact)'
        self.V, self.L = V, L
        self.sides = [(s, (s + 1) % k) for s in range(k)]
        self.t = []
        for s, (i, _) in enumerate(self.sides):
            P = V[i]
            LP = mshift(mconj(P, N), 4 * M * self.m[s], N)
            self.t.append(madd(P, mneg(LP)))

    def g_apply(self, g, p):
        k, e = g
        q = mconj(p, self.N) if e else p
        return mshift(q, 4 * k * self.M, self.N)

    def cross2(self, u, v):
        N = self.N
        w = madd(mmul(mconj(u, N), v, N), mneg(mmul(u, mconj(v, N), N)))
        return mshift(w, 3 * self.d * self.M, N)

    g_across = FD.Triangle.g_across
    red = FD.Triangle.red


def perp_dirs_poly(P):
    return [((P.d + 2 * m) * P.M) % P.N for m in P.m]


def _above(T, F, A, B):
    """s936: +1 if segment B lies above segment A (further along ω) over their common open x-range, −1 if below.
    A, B = (p, q) monomial points, each directed with x increasing, non-crossing (sides of a simple polygon), overlapping
    in x on an open interval. In the frame x = 2ω×p, y = 2ω·p (orientation-REVERSING), p is above the directed line P→Q
    iff cross(Q − P, p − P) < 0 in the original plane. Standard sweep comparator: decide B against line A if B's endpoints
    lie weakly on one side of it, else A against line B."""
    def side(P, Q, p):                                        # +1 above, −1 below, 0 on the line
        return -F.sign(F.red(T.cross2(madd(Q, mneg(P)), madd(p, mneg(P)))))
    s0, s1 = side(*A, B[0]), side(*A, B[1])
    if s0 >= 0 and s1 >= 0 and (s0 or s1):
        return 1
    if s0 <= 0 and s1 <= 0 and (s0 or s1):
        return -1
    t0, t1 = side(*B, A[0]), side(*B, A[1])
    assert (t0 or t1) and not (t0 * t1 < 0), 'sides cross or are collinear'
    return -1 if (t0 > 0 or t1 > 0) else 1


def _pair_sweep(T, F, W, xs, ent, ext):
    """s936: for a (possibly non-convex) simple polygon copy, the first exit side hit along ω from each entry point.
    Break the x-axis at every vertex x; on each elementary interval the covering sides are totally ordered by height
    (`_above`) and alternate entry, exit, entry, … from below; pair consecutive ones. Returns {entry s: [(u, w, s2)]},
    maximal x-runs with one exit side (adjacent elementary intervals with the same pair merged)."""
    bks = []
    for x in xs:
        if not any(F.cmp(x, y) == 0 for y in bks):
            bks.append(x)
    bks.sort(key=functools.cmp_to_key(F.cmp))
    seg = {}
    for (s, lo, hi) in ent + ext:
        i, k2 = T.sides[s]
        p, q = W[i], W[k2]
        if F.cmp(xs[i], xs[k2]) > 0:
            p, q = q, p
        seg[s] = (p, q, lo, hi)
    kind = {s: 0 for s, _, _ in ent}
    kind.update({s: 1 for s, _, _ in ext})
    out = {s: [] for s, _, _ in ent}
    for u, w in zip(bks, bks[1:]):
        cov = [s for s in seg if F.cmp(seg[s][2], u) <= 0 and F.cmp(seg[s][3], w) >= 0]
        cov.sort(key=functools.cmp_to_key(lambda a, b: -_above(T, F, seg[a][:2], seg[b][:2])))   # ascending height
        assert len(cov) % 2 == 0 and all(kind[s] == (n % 2) for n, s in enumerate(cov)), 'entry/exit do not alternate'
        for n in range(0, len(cov), 2):
            s, s2 = cov[n], cov[n + 1]
            run = out[s]
            if run and run[-1][2] == s2 and F.cmp(run[-1][1], u) == 0:
                run[-1] = (run[-1][0], w, s2)
            else:
                run.append((u, w, s2))
    return out


def _build_iet_poly(T, j, sweep=None):
    """`flow_decomp._build_iet` with the inward test by orientation (CCW base copy). Convex: every entry side is paired
    with every exit side whose x-range overlaps (one segment per flow line). Non-convex, or `sweep=True` (the gate):
    the first exit side along ω, by `_pair_sweep`."""
    if sweep is None:
        sweep = not T.convex
    F, d, N = FD._fast(T), T.d, T.N
    om = {e % N: c for e, c in j.items()} if isinstance(j, dict) else {j % N: 1}
    X = lambda p: F.red(T.cross2(om, p))
    copies = [(k, e) for e in (0, 1) for k in range(d)]
    info = {}
    for g in copies:
        W = [T.g_apply(g, v) for v in T.V]
        xs = [X(w) for w in W]
        orient = 1 if g[1] == 0 else -1
        ent, ext = [], []
        for s, (i, k2) in enumerate(T.sides):
            e = madd(W[k2], mneg(W[i]))
            c_om = F.sign(F.red(T.cross2(e, om)))
            if c_om == 0:
                continue
            lo, hi = (xs[i], xs[k2]) if F.cmp(xs[i], xs[k2]) < 0 else (xs[k2], xs[i])
            (ent if c_om == orient else ext).append((s, lo, hi))
        info[g] = (xs, ent, ext, _pair_sweep(T, F, W, xs, ent, ext) if sweep else None)
    order = [(g, s, lo, hi) for g in copies for (s, lo, hi) in info[g][1]]
    off = {}
    pos = F.zero
    for g, s, lo, hi in order:
        off[(g, s)] = (pos, lo, hi)
        pos = F.add(pos, F.sub(hi, lo))
    L = pos
    pieces = []
    for g, s, lo, hi in order:
        O, _, _ = off[(g, s)]
        if sweep:
            runs = info[g][3][s]
        else:
            runs = []
            for (s2, lo2, hi2) in info[g][2]:
                u = lo if F.cmp(lo, lo2) >= 0 else lo2
                w = hi if F.cmp(hi, hi2) <= 0 else hi2
                if F.cmp(w, u) > 0:
                    runs.append((u, w, s2))
        for (u, w, s2) in runs:
            g2 = T.g_across(g, s2)
            key2 = (g2, s2)
            assert key2 in off, ('exit side is not the neighbour entry', g, s2, g2)
            O2, lo2e, _ = off[key2]
            shift = X(T.g_apply(g, T.t[s2]))
            dom = F.add(O, F.sub(u, lo))
            img = F.add(O2, F.sub(F.sub(u, shift), lo2e))
            pieces.append((dom, F.sub(w, u), img, (g, s), key2))
    n = len(pieces)
    top = sorted(range(n), key=functools.cmp_to_key(lambda p, q: F.cmp(pieces[p][0], pieces[q][0])))
    bot = sorted(range(n), key=functools.cmp_to_key(lambda p, q: F.cmp(pieces[p][2], pieces[q][2])))
    for seq, idx in ((top, 0), (bot, 2)):
        cur = F.zero
        for p in seq:
            assert F.cmp(pieces[p][idx], cur) == 0, 'intervals do not tile [0, L)'
            cur = F.add(cur, pieces[p][1])
        assert F.cmp(cur, L) == 0
    lengths = {p: pieces[p][1] for p in range(n)}
    return pieces, lengths, top, bot, {'L': L, 'ncopies': len(copies), 'nentry': len(order), 'off': off}


def decompose_poly(T, j, cap=5000, ctl=None):
    pieces, lengths, top, bot, meta = _build_iet_poly(T, j)
    F = FD._fast(T)
    Dc = Decomposer(F, lengths)
    comps = Dc.run(top, bot, cap, ctl=ctl) if ctl is not None else Dc.run(top, bot, cap)
    out, total = [], F.zero
    for c in comps:
        meas, mass = F.zero, F.zero
        for l in c.top:
            meas = F.add(meas, Dc.lam[l])
            mass = F.add(mass, F.smul(Dc.h[l], Dc.lam[l]))
        total = F.add(total, mass)
        out.append({'kind': c.kind, 'size': len(c.top), 'steps': c.steps,
                    'mass': F.fval(mass) / F.fval(meta['L']),
                    'period': Dc.h[c.top[0]] if len(c.top) == 1 else None, 'cert': c.cert})
    assert F.cmp(total, meta['L']) == 0, 'Kac mass balance failed'
    kinds = {}
    for o in out:
        kinds[o['kind']] = kinds.get(o['kind'], 0) + 1
    return {'a': T.abc, 'j': j, 'N': T.N, 'npieces': len(pieces), 'counts': kinds, 'components': out}


def perp_keys(T, j):
    """s936: the entry keys (g, s) whose side copy is PERPENDICULAR to ω — crossing one = hitting a side perpendicularly"""
    om = {e % T.N: c for e, c in j.items()} if isinstance(j, dict) else {j % T.N: 1}
    iom = mshift(om, T.d * T.M, T.N)                          # i·ω
    keys = []
    for g in [(k, e) for e in (0, 1) for k in range(T.d)]:
        W = [T.g_apply(g, v) for v in T.V]
        for s, (i, k2) in enumerate(T.sides):
            if not any(T.red(T.cross2(iom, madd(W[k2], mneg(W[i]))))):
                keys.append((g, s))
    return keys


def crossings_poly(T, j, keys=None, cap=5000, sweep=None):
    """s936: `flow_decomp.crossings` for a PolyT — every component, and for a cylinder its crossing count at each key
    (default `perp_keys`). Kac asserted exactly on the return times."""
    keys = perp_keys(T, j) if keys is None else keys
    pieces, lengths, top, bot, meta = _build_iet_poly(T, j, sweep)
    idx = {key: i for i, key in enumerate(keys)}
    F = FD._fast(T)
    Dc = Decomposer(F, lengths)
    zero = (0,) * len(keys)
    Dc.h = {}
    for p in lengths:
        k = list(zero)
        if pieces[p][3] in idx:
            k[idx[pieces[p][3]]] = 1
        Dc.h[p] = FD._HK(1, tuple(k))
    comps = Dc.run(top, bot, cap)
    total, out = F.zero, []
    for c in comps:
        for l in c.top:
            total = F.add(total, F.smul(Dc.h[l].h, Dc.lam[l]))
        one = len(c.top) == 1
        out.append({'kind': c.kind, 'size': len(c.top), 'period': Dc.h[c.top[0]].h if one else None,
                    'k': Dc.h[c.top[0]].k if one else None,
                    'mass': F.fval(F.add(F.zero, F.smul(Dc.h[c.top[0]].h, Dc.lam[c.top[0]]))) / F.fval(meta['L'])
                    if one else None})
    assert F.cmp(total, meta['L']) == 0, 'Kac mass balance failed'
    return out, keys


def mplus_poly(T, j, cap=5000, sweep=None):
    """s936: the exact (M⁺) verdict — every regular orbit of the class meets a side perpendicularly ⟺ every component is a
    cylinder AND every cylinder crosses a perpendicular side copy. A cylinder crossing none is FREE. Returns counts ⊕
    the free cylinders' masses (fractions of the total transversal) ⊕ 'Mplus'."""
    out, keys = crossings_poly(T, j, None, cap, sweep)
    cyl = [o for o in out if o['kind'] == 'cylinder']
    free = [o for o in cyl if not any(o['k'])]
    other = [o for o in out if o['kind'] != 'cylinder']
    return {'a': T.abc, 'd': T.d, 'j': j, 'nkeys': len(keys), 'ncyl': len(cyl), 'nfree': len(free),
            'free_mass': sorted(o['mass'] for o in free), 'nother': len(other),   # free_mass: transversal, NOT area
            'other_kinds': sorted({o['kind'] for o in other}), 'cp': not other,
            'Mplus': not other and not free}


def marked_germs(T, j, i, off):
    """s938: the one-sided GERMS of vertex V_i's leaf on the transversal. In copy g (convex), the leaf through W = g·V_i is
    the chord at x = X(W), so the points just left / right of it inside g enter through the entry sides whose x-range
    contains X(W): a germ (pos, −1) for each with lo < X(W), (pos, +1) for each with X(W) < hi. Returns
    [(pos, side, g, s)]; each germ lies in exactly ONE component (components are finite unions of intervals)."""
    assert T.convex, 'marked_germs: convex copies only'
    F, N = FD._fast(T), T.N
    om = {e % N: c for e, c in j.items()} if isinstance(j, dict) else {j % N: 1}
    out = []
    for (g, s), (O, lo, hi) in off.items():
        x = F.red(T.cross2(om, T.g_apply(g, T.V[i])))
        c1, c2 = F.cmp(lo, x), F.cmp(x, hi)
        if c1 <= 0 and c2 <= 0:
            pos = F.add(O, F.sub(x, lo))
            if c1 < 0:
                out.append((pos, -1, g, s))
            if c2 < 0:
                out.append((pos, 1, g, s))
    return out


def vertex_points(T, i):
    """s938: the copies sharing V_i's corner, grouped into the surface points over V_i (orbits of the two side
    reflections at V_i). Returns {g: point index}."""
    k = T.k
    s_in, s_out = (i - 1) % k, i
    copies = [(kk, e) for e in (0, 1) for kk in range(T.d)]
    pid, n = {}, 0
    for g in copies:
        if g in pid:
            continue
        stack = [g]
        pid[g] = n
        while stack:
            h = stack.pop()
            for s in (s_in, s_out):
                h2 = T.g_across(h, s)
                if h2 not in pid:
                    pid[h2] = n
                    stack.append(h2)
        n += 1
    return pid


class MarkDecomposer(Decomposer):
    """s938: `Decomposer` carrying TRACKED GERMS through the induction, to say which component holds a marked point.
    A germ is (y, σ): the point just left (σ = −1) / right (+1) of y on the current interval of its component (left end
    0). Every interval change removes a LEFT part [0, R) — Zorich step (top or bottom: R = the winner's length loss),
    connection merge (R = the merged length), split (`_reduce`: the tail [P, total) becomes a new component) — and a germ
    in a removed part is pushed by the PRE-STEP map until it re-enters: same orbit, so same component. Gates in
    `marked_components` (forward image in the same component; a cylinder germ returns after exactly its period)."""

    def __init__(self, F, lengths, germs):
        super().__init__(F, lengths)
        self.pts = [[None, y, s] for (y, s) in germs]   # [comp, y, σ]

    def _gcmp(self, y, s, b):
        c = self.F.cmp(y, b)
        return c if c else s

    def _prefix(self, lst, lam):
        F, cur, d = self.F, self.F.zero, {}
        for l in lst:
            d[l] = cur
            cur = F.add(cur, lam[l])
        return d

    def _advance(self, c, R, otop, obot, olam):
        """the germs of c after the left part [0, R) is removed; (otop, obot, olam) = the pre-step IET"""
        F, pt, pb = self.F, None, None
        for p in self.pts:
            if p[0] is not c:
                continue
            while self._gcmp(p[1], p[2], R) < 0:
                if pt is None:
                    pt, pb = self._prefix(otop, olam), self._prefix(obot, olam)
                y, s = p[1], p[2]
                l = next(l for l in otop if self._gcmp(y, s, pt[l]) > 0 and self._gcmp(y, s, F.add(pt[l], olam[l])) < 0)
                p[1] = self._jump(y, s, pt[l], F.add(pt[l], olam[l]), F.sub(pb[l], pt[l]), R)
            p[1] = F.sub(p[1], R)

    def _jump(self, y, s, lo, hi, D, R):
        """T_old^k(y) for y in label [lo, hi) with translation D: k = 1, or — a TWISTING label maps into itself, with
        k up to the Dehn-twist multiplicity (huge) — the closed-form k at which the orbit leaves [lo, min(hi, R))
        (exact floor division; germ ties broken by σ). Every intermediate point stays in the label, so each step is +D."""
        F = self.F
        sg = F.sign(D)
        if sg > 0:
            b = R if self._gcmp(R, 1, hi) < 0 else hi           # stop at R or at the label's right end
            q = self._floordiv(F.sub(b, y), D) if F.cmp(b, y) > 0 else 0
            k = q if q >= 1 and self._gcmp(F.add(y, F.smul(q, D)), s, b) > 0 else q + 1
            return F.add(y, F.smul(k, D))
        if sg < 0:
            nD = F.neg(D)
            q = self._floordiv(F.sub(y, lo), nD)
            k = q if q >= 1 and self._gcmp(F.sub(y, F.smul(q, nD)), s, lo) < 0 else q + 1
            return F.sub(y, F.smul(k, nD))
        return y

    def _snap(self, c):
        return list(c.top), list(c.bot), {l: self.lam[l] for l in c.top}

    def run(self, top, bot, cap):
        F = self.F
        c0 = Comp(top, bot)
        for p in self.pts:
            p[0] = c0
        work, done = [c0], []
        while work:
            c = work.pop()
            while True:
                rest = self._reduce(c)
                if rest is not None:
                    work.append(rest)
                    c.trk = c.saf0 = None
                    P = F.zero
                    for l in c.top:
                        P = F.add(P, self.lam[l])
                    for p in self.pts:
                        if p[0] is c and self._gcmp(p[1], p[2], P) > 0:
                            p[0], p[1] = rest, F.sub(p[1], P)
                if len(c.top) == 1:
                    c.kind = 'cylinder'
                    break
                t0, b0 = c.top[0], c.bot[0]
                if t0 != b0 and F.cmp(self.lam[t0], self.lam[b0]) == 0:
                    tracked = any(p[0] is c for p in self.pts)
                    snap = self._snap(c) if tracked else None
                    it = c.bot.index(t0)
                    c.bot[it] = b0
                    self.h[b0] += self.h[t0]
                    del c.bot[0]
                    del c.top[0]
                    c.trk = c.saf0 = None
                    if tracked:
                        self._advance(c, snap[2][t0], *snap)
                    continue
                if c.steps >= cap:
                    c.kind = 'undetermined'
                    break
                if self.autosimilar(c):
                    c.kind, c.cert = 'no_periodic', 'autosimilar'
                    break
                if c.steps % FD.BOSH_EVERY == 0 and c.steps > 0 and self.saf_nonzero(c):
                    y = self.boshernitzan(c)
                    if y is not None:
                        c.kind, c.cert = 'no_periodic', y
                        break
                c.steps += 1
                hit = False
                for first, second in ((c.top, c.bot), (c.bot, c.top)):
                    tracked = any(p[0] is c for p in self.pts)
                    snap = self._snap(c) if tracked else None
                    w = first[0]
                    lw = self.lam[w]
                    hit = self._zorich(first, second)
                    if tracked:
                        R = F.sub(lw, self.lam[w])
                        if F.sign(R) != 0:
                            self._advance(c, R, *snap)
                    if hit:
                        break
            done.append(c)
        return done


def marked_components(T, j, verts, cap=5000, gate=True):
    """s938: EXACT pole membership — for each marked vertex V_i (i in `verts`; a π/2 vertex is a regular point of the
    translation surface, a pole of S/τ) and each surface point over it (`vertex_points`), the components its germs lie
    in. Returns {'comps': [{kind, size, period, mass}], 'marks': {i: {point: sorted comp indices}}, 'germs': n}.
    GATES (exact, `gate=True`): each germ's T-image is tracked too and must end in the SAME component; a germ in a
    cylinder must return to itself under the original T after exactly that cylinder's period (no earlier)."""
    pieces, lengths, top, bot, meta = _build_iet_poly(T, j)
    F = FD._fast(T)
    n = len(pieces)

    starts = [pieces[p][0] for p in top]                  # `top` = pieces sorted by domain start

    def Tmap(y, s):                                       # binary search: the last piece whose start is left of the germ
        lo, hi = 0, n - 1
        while lo < hi:
            mid = (lo + hi + 1) // 2
            if (F.cmp(y, starts[mid]) or s) > 0:
                lo = mid
            else:
                hi = mid - 1
        p = top[lo]
        dom, ln, img = pieces[p][0], pieces[p][1], pieces[p][2]
        assert (F.cmp(y, dom) or s) > 0 and (F.cmp(y, F.add(dom, ln)) or s) < 0, 'germ not in its piece'
        return F.add(img, F.sub(y, dom))

    germs, gmeta = [], []
    for i in verts:
        pid = vertex_points(T, i)
        for (pos, s, g, sd) in marked_germs(T, j, i, meta['off']):
            germs.append((pos, s))
            gmeta.append((i, pid[g]))
    ng = len(germs)
    if gate:
        germs = germs + [(Tmap(y, s), s) for (y, s) in germs]
    Dc = MarkDecomposer(F, lengths, germs)
    comps = Dc.run(top, bot, cap)
    cidx = {id(c): k for k, c in enumerate(comps)}
    where = [cidx[id(p[0])] for p in Dc.pts]
    out, total = [], F.zero
    for c in comps:
        mass = F.zero
        for l in c.top:
            mass = F.add(mass, F.smul(Dc.h[l], Dc.lam[l]))
        total = F.add(total, mass)
        out.append({'kind': c.kind, 'size': len(c.top), 'period': Dc.h[c.top[0]] if len(c.top) == 1 else None,
                    'mass': F.fval(mass) / F.fval(meta['L'])})
    assert F.cmp(total, meta['L']) == 0, 'Kac mass balance failed'
    if gate:
        assert where[:ng] == where[ng:], 'gate: a germ and its T-image ended in different components'
        for k in range(ng):
            c = out[where[k]]
            if c['kind'] == 'cylinder':
                y, s = germs[k]
                z = y
                for step in range(1, c['period'] + 1):
                    z = Tmap(z, s)
                    if F.cmp(z, y) == 0:
                        break
                assert step == c['period'] and F.cmp(z, y) == 0, ('gate: cylinder return', step, c['period'])
    marks = {}
    for k, (i, q) in enumerate(gmeta):
        marks.setdefault(i, {}).setdefault(q, set()).add(where[k])
    marks = {i: {q: sorted(v) for q, v in d.items()} for i, d in marks.items()}
    return {'comps': out, 'marks': marks, 'germs': ng}


def triangle_poly(abc):
    """a triangle as a PolyT (law of sines: |V_jV_{j+1}| = 2 sin(a_{j+2}π/d))"""
    d = sum(abc)
    return PolyT(list(abc), d, [sin_mono(abc[(j + 2) % 3], d) for j in range(3)])


def trapezoid_cut(a0, a1, d, p, q):
    a2 = d - a0 - a1
    assert a2 > 0 and 0 < p < q
    L = [mscale(sin_mono(a2, d), q), mscale(sin_mono(a0, d), p), mscale(sin_mono(a2, d), q - p),
         mscale(sin_mono(a1, d), p)]
    return PolyT([a0, a1, d - a1, d - a0], d, L)


def kite_poly(abc, Ne, vtx):
    """§s929 (2)'s kite (a, b/2, c, b/2)/(Ne/2) doubled from T = abc/Ne across the side opposite the even vertex vtx"""
    b = abc[vtx]
    a, c = abc[(vtx + 1) % 3], abc[(vtx - 1) % 3]
    d = Ne // 2
    sa, sc = sin_half(a, d), sin_half(c, d)
    return PolyT([a, b // 2, c, b // 2], d, [sc, sa, sa, sc])


def _summ(r):
    return (r['counts'], sorted(round(c['mass'], 12) for c in r['components']))


def _selftest():
    for abc in [(7, 8, 15), (7, 7, 16), (3, 4, 5)]:
        P = triangle_poly(abc)
        for s, j in enumerate(perp_dirs_poly(P)):
            r1 = decompose_poly(P, j)
            r2 = FD.decompose(list(abc), j)
            # the triangle's own convention: compare by direction (the same ω), the copies may be ordered differently
            assert _summ(r1) == (r2['counts'], sorted(round(c['mass'], 12) for c in r2['components'])), (abc, j, r1, r2)
            print('ok', abc, 'side', s, 'j', j, r1['counts'])
    # a genuine polygon: the s929 kites have T's surface; their perp class must carry no_periodic iff T's pair class does
    for abc, Ne, vtx, want in [((30, 41, 67), 138, 0, True), ((7, 8, 15), 30, 1, False), ((23, 29, 58), 110, 2, False)]:
        P = kite_poly(abc, Ne, vtx)
        r = decompose_poly(P, perp_dirs_poly(P)[0])
        bad = 'no_periodic' in r['counts']
        assert bad == want, (abc, r['counts'])
        print('ok kite', P.abc, '/', P.d, r['counts'])
    # s934: the field multiplier M — the same polygon with its lengths lifted to ζ_{4dM} gives the identical decomposition
    for abc, Ne, vtx in [((30, 41, 67), 138, 0), ((7, 8, 15), 30, 1)]:
        P1 = kite_poly(abc, Ne, vtx)
        P2 = PolyT(P1.a, P1.d, [mlift(x, 2, 4 * P1.d) for x in P1.L], M=2)
        r1, r2 = (decompose_poly(P, perp_dirs_poly(P)[0]) for P in (P1, P2))
        assert _summ(r1) == _summ(r2), (abc, r1['counts'], r2['counts'])
        print('ok M=2 lift', P1.abc, '/', P1.d, r2['counts'])


    # s936: the sweep pairing reproduces the convex pieces; (M⁺) controls; a non-convex fold octagon
    for P in [triangle_poly((7, 8, 15)), kite_poly((30, 41, 67), 138, 0), trapezoid_cut(2, 3, 9, 1, 3)]:
        for j in perp_dirs_poly(P):
            p1, p2 = (_build_iet_poly(P, j, sweep=s)[0] for s in (False, True))
            assert sorted(map(str, p1)) == sorted(map(str, p2)), (P.abc, j)
    print('ok sweep = convex pairing')
    rows = free = 0
    for N in range(4, 17, 2):
        for a in range(1, N):
            for b in range(a, N - a):
                c = N - a - b
                if c >= b and math.gcd(math.gcd(a, b), c) == 1 and all(N % x for x in (a, b, c)):
                    P = triangle_poly((a, b, c))
                    for j in perp_dirs_poly(P):
                        rows += 1
                        free += mplus_poly(P, j)['nfree'] > 0
    assert (rows, free) == (21, 14), (rows, free)              # the N ≤ 16 part of §s927 (2)'s 55/87
    assert all(mplus_poly(triangle_poly(t), perp_dirs_poly(triangle_poly(t))[0])['Mplus'] for t in [(2, 3, 4), (3, 5, 7)])
    one = {0: 1}
    assert mplus_poly(PolyT([3, 2, 3, 2], 5, [one] * 4), 5)['nfree'] == 2
    assert mplus_poly(PolyT([3, 2, 3, 2], 5, [{0: 1000}, {0: 2001}] * 2), 5)['Mplus']
    print('ok M+ controls (even N <= 16 free rows 14/21; rhombus free; 1000:2001 M+)')
    T, Ne = (30, 41, 67), 138
    P = PolyT([134, 30, 67, 15, 82, 30, 41, 15], 69, [sin_mono(T[(t + 2) % 3], Ne) for t in [2, 2, 2, 0, 0, 0, 0, 2]], M=2)
    assert not P.convex
    r = mplus_poly(P, perp_dirs_poly(P)[0])
    assert not r['cp'], r
    print('ok non-convex octagon', P.abc, '/ 69', r['ncyl'], 'cyl', r['nother'], r['other_kinds'])
    # s938: marked components — the decomposition is the unmarked one; T₀'s ⊥hyp class puts its right angle in the
    # no-periodic components on 8 of 15 points (the gates run inside)
    P = triangle_poly((7, 8, 15))
    j = perp_dirs_poly(P)[2]
    r = marked_components(P, j, [2])
    assert sorted(round(c['mass'], 12) for c in r['comps']) == _summ(decompose_poly(P, j))[1]
    nonper = sum(any(r['comps'][k]['kind'] != 'cylinder' for k in v) for v in r['marks'][2].values())
    assert (len(r['marks'][2]), nonper) == (15, 8), r['marks']
    print('ok marked components T0 hyp class: right-angle points in no-periodic comps', nonper, '/ 15')


if __name__ == '__main__':
    _selftest()
