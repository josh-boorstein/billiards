"""[Deposit] Helper, imported by the scripts of this directory. Paper: "Perpendicular billiards in rational polygons of odd denominator". Research name: engine/poly_exact.py, copied verbatim except imports (deposit names), paths (data/ and figs/ under odd_polygons/) and this paragraph. Run through run_all.py; the usage lines below are the research repo's. Bracketed IDs and section marks (§sNNN) refer to the project's internal ledgers.

poly_exact — EXACT PRONG TRACKING ON A ONE-PARAMETER FAMILY OF CONVEX RATIONAL POLYGONS, sides AFFINE in `t`, over
K = ℚ(cos(π/(2N))). s931 (queue `s930-a`).

WHY. On a family `Π_t` with fixed angles `π a_i/N` and side lengths affine in `t`, every unfolded vertex is AFFINE in `t`
(reflections have constant linear part), so a prong's word is fixed by the SIGNS of the affine forms
`f_W(t) = cross(u, W(t) − V_i(t))` over the vertices `W` of the copies it visits (convexity: the chord through a convex copy is
fixed by the vertex sign pattern). So a prong's word is valid EXACTLY on the open interval between the nearest roots of those
forms, and the word changes at each root (the chord passes a vertex there). `track()` walks a slot's prong across `[T0, T1]`
breakpoint by breakpoint with every comparison an exact field operation (signs CERTIFIED: a 50-digit evaluation with an
explicit error bound, falling back to `numfield`'s interval sign).

OBJECTS.
  Fam(a, N, L)        a = angle numerators (Σ a = (k − 2)N, all a_i < N: CONVEX asserted); L = k side lengths, each a pair
                      (c0, c1) of K-elements meaning c0 + c1·t (closure checked EXACTLY). Fam.closing(a, N, free) builds L from
                      the first k − 2 sides (affine) by exact closure.
    .K  .cos(k)  .sin(k)  (k in units π/(2N))   .slots()  (the perpendicular class of side 0, as Poly.slots)
    .prong(i, d, t0, cap)  exact walk at t0 (int / Fraction / K-element) → dict: word (sides hit), end ('perp', j) /
                      ('vertex', label) / ('cap',), length (an affine form: the billiard length, dot(u, W − V_i) at a vertex
                      end; None otherwise), lo / hi (the open validity interval, K-elements, None = unbounded), persistent
                      (True iff the vertex end holds identically in t), arr (the arrival direction at a vertex end).
    .track(i, d, T0, T1, cap)  list of pieces (lo, hi, word, end, length) covering (T0, T1) with the breakpoints between them;
                      gaps are bisected until the exact left end of each piece equals the previous right end.
  Aff                 c0 + c1·t over K: + − scalar·, .at(t0), .root() (= −c0/c1), .const().

TRAPS.
  ⚠ CONVEX polygons only (the exit side is read off the sign pattern); asserted on the angles.
  ⚠ the sides must stay positive on [T0, T1] — the caller passes an interval inside the family's range.
  ⚠ a breakpoint is where the chord passes a vertex of a VISITED copy — the prong ends there AT that `t` (a vertex
    connection, or a vertex hit before the old end); `prong(i, d, t*)` at the exact breakpoint returns that configuration.

stdlib + numfield (sympy at construction only) + mpmath.
"""
import _setup  # noqa: F401  (deposit: this directory on the path, cwd = odd_polygons/)
from fractions import Fraction
from math import gcd

import mpmath
import sympy as sp

from numfield import RealField, Elt

_FIELDS = {}


def field(N):
    """K = ℚ(θ), θ = 2cos(π/(2N))."""
    if N not in _FIELDS:
        x = sp.Symbol('x')
        mp_ = sp.minimal_polynomial(2 * sp.cos(sp.pi / (2 * N)), x)
        co = [int(c) for c in reversed(sp.Poly(mp_, x).all_coeffs())]
        _FIELDS[N] = RealField(co, float(2 * mpmath.cos(mpmath.pi / (2 * N))))
    return _FIELDS[N]


class _Sign:
    """certified sign with a fast path: 60-digit evaluation with an explicit error bound."""
    def __init__(self, K):
        self.K = K
        with mpmath.workdps(80):
            self.th = K.approx(70)
            self.pw = [self.th ** i for i in range(K.deg + 1)]

    def __call__(self, e):
        if not e.c:
            return 0
        with mpmath.workdps(60):
            v = mpmath.mpf(0)
            s = mpmath.mpf(0)
            for i, x in enumerate(e.c):
                term = mpmath.mpf(x.numerator) / x.denominator * self.pw[i]
                v += term
                s += abs(term)
            if abs(v) > s * mpmath.mpf(10) ** -50 + mpmath.mpf(10) ** -55:
                return 1 if v > 0 else -1
        return e.sign()


class Aff:
    __slots__ = ('c0', 'c1')

    def __init__(self, c0, c1):
        self.c0, self.c1 = c0, c1

    def __add__(self, o):
        return Aff(self.c0 + o.c0, self.c1 + o.c1)

    def __sub__(self, o):
        return Aff(self.c0 - o.c0, self.c1 - o.c1)

    def scale(self, k):
        return Aff(self.c0 * k, self.c1 * k)

    def at(self, t0):
        return self.c0 + self.c1 * t0

    def const(self):
        return self.c1.iszero()

    def root(self):
        return -self.c0 / self.c1

    def __repr__(self):
        return f'Aff({float(self.c0):.6g} + {float(self.c1):.6g} t)'


class Fam:
    def __init__(self, a, N, L):
        k = len(a)
        assert sum(a) == (k - 2) * N and all(0 < x < N for x in a), 'convex rational polygon expected'
        self.a, self.N, self.k = list(a), N, k
        self.K = K = field(N)
        self.sgn = _Sign(K)
        th = K.gen
        c = [K.of(2), th]                       # c[j] = 2cos(jπ/(2N))
        for j in range(2, 4 * N + 1):
            c.append(th * c[-1] - c[-2])
        self._c2 = c
        idx = [0]
        for j in range(1, k):
            idx.append(idx[-1] + 2 * (N - a[j]))
        self.idx = [x % (4 * N) for x in idx]
        self.m = [x % (2 * N) for x in self.idx]
        self.L = [Aff(K.of(l[0]), K.of(l[1])) for l in L]
        Vx, Vy = [Aff(K.zero, K.zero)], [Aff(K.zero, K.zero)]
        for j in range(k):
            Vx.append(Vx[-1] + self.L[j].scale(self.cos(self.idx[j])))
            Vy.append(Vy[-1] + self.L[j].scale(self.sin(self.idx[j])))
        assert Vx[k].c0.iszero() and Vx[k].c1.iszero() and Vy[k].c0.iszero() and Vy[k].c1.iszero(), 'does not close'
        self.V = list(zip(Vx[:k], Vy[:k]))

    @classmethod
    def closing(cls, a, N, free):
        """free = k − 2 affine side lengths (pairs of int/Fraction/K-elt); the last two by exact closure."""
        k = len(a)
        K = field(N)
        tmp = cls.__new__(cls)
        tmp.N, tmp.K = N, K
        th = K.gen
        c = [K.of(2), th]
        for j in range(2, 4 * N + 1):
            c.append(th * c[-1] - c[-2])
        tmp._c2 = c
        idx = [0]
        for j in range(1, k):
            idx.append(idx[-1] + 2 * (N - a[j]))
        idx = [x % (4 * N) for x in idx]
        free = [(K.of(f[0]), K.of(f[1])) for f in free]
        out = list(free)
        sol = []
        for comp in (0, 1):
            sx = [sum((f[comp] * tmp.cos(idx[j]) for j, f in enumerate(free)), K.zero)]
            sy = [sum((f[comp] * tmp.sin(idx[j]) for j, f in enumerate(free)), K.zero)]
            a11, a12 = tmp.cos(idx[k - 2]), tmp.cos(idx[k - 1])
            a21, a22 = tmp.sin(idx[k - 2]), tmp.sin(idx[k - 1])
            det = a11 * a22 - a12 * a21
            x = (-sx[0] * a22 + a12 * sy[0]) / det
            y = (-a11 * sy[0] + a21 * sx[0]) / det
            sol.append((x, y))
        out.append((sol[0][0], sol[1][0]))
        out.append((sol[0][1], sol[1][1]))
        return cls(a, N, out)

    # trig on the π/(2N) grid
    def cos(self, j):
        j %= 4 * self.N
        return self._c2[j] * Fraction(1, 2)

    def sin(self, j):
        return self.cos(self.N - j)

    def slots(self):
        N, k = self.N, self.k
        d0 = (self.idx[0] + N) % (4 * N)
        orb, todo = {d0}, [d0]
        while todo:
            d = todo.pop()
            for mj in self.m:
                e = (2 * mj - d) % (4 * N)
                if e not in orb:
                    orb.add(e)
                    todo.append(e)
        self.C = orb
        out = []
        for i in range(k):
            for s in range(1, 2 * self.a[i]):
                d = (self.idx[i] + s) % (4 * N)
                if d in orb:
                    out.append((i, d))
        return out

    def _reflect(self, P, A, phi):
        """reflect point P (pair of Aff) in the line through A with direction index phi"""
        c2, s2 = self.cos(2 * phi), self.sin(2 * phi)
        dx, dy = P[0] - A[0], P[1] - A[1]
        return (A[0] + dx.scale(c2) + dy.scale(s2), A[1] + dx.scale(s2) - dy.scale(c2))

    def prong(self, i, d, t0, cap=4000, K_t=None):
        K, N, k, sg = self.K, self.N, self.k, self.sgn
        t0 = K.of(t0) if not isinstance(t0, Elt) else t0
        ux, uy = self.cos(d), self.sin(d)
        S = self.V[i]
        copy = list(self.V)
        eps, r = 1, 0                           # copy linear part: direction x -> r + eps*x
        entry = None
        word, forms = [], []

        def form(W):
            return (W[1] - S[1]).scale(ux) - (W[0] - S[0]).scale(uy)   # cross(u, W − S)
        out = None
        for step in range(cap):
            fs = []
            for j in range(k):
                if step == 0 and j == i:
                    fs.append(None)
                    continue
                f = form(copy[j])
                fs.append(f)
                forms.append(f)
            sv = [None if f is None else sg(f.at(t0)) for f in fs]
            zero = [j for j in range(k) if sv[j] == 0]
            if zero:
                assert len(zero) == 1, 'two vertices on the chord'
                W = copy[zero[0]]
                ln = (W[0] - S[0]).scale(ux) + (W[1] - S[1]).scale(uy)
                f = fs[zero[0]]
                darr = (eps * (d - r)) % (4 * N)
                out = dict(end=('vertex', zero[0]), length=ln, persistent=f.c0.iszero() and f.c1.iszero(),
                           arr=darr, steps=step)
                break
            ex = None
            for j in range(k):
                if j == entry:
                    continue
                s1, s2 = sv[j], sv[(j + 1) % k]
                if s1 is None or s2 is None:
                    continue
                if s1 != s2:
                    assert ex is None, 'two exit sides (non-convex?)'
                    ex = j
            assert ex is not None
            word.append(ex)
            dcopy = (eps * (d - r)) % (4 * N)  # billiard direction in the original frame
            if (dcopy - self.m[ex] - N) % (2 * N) == 0:
                out = dict(end=('perp', ex), length=None, persistent=None, arr=dcopy, steps=step)
                break
            phi = (r + eps * self.idx[ex]) % (4 * N)
            A = copy[ex]
            copy = [self._reflect(P, A, phi) for P in copy]
            eps, r = -eps, (2 * phi - r) % (4 * N)
            entry = ex
        else:
            out = dict(end=('cap',), length=None, persistent=None, arr=None, steps=cap)
        lo = hi = None
        for f in forms:
            if f.c1.iszero():
                continue
            v = sg(f.at(t0))
            if v == 0:
                continue                        # the vertex end itself, at an isolated breakpoint
            rt = f.root()
            s = sg(rt - t0)
            if s > 0 and (hi is None or sg(rt - hi) < 0):
                hi = rt
            elif s < 0 and (lo is None or sg(rt - lo) > 0):
                lo = rt
        if out['end'][0] == 'vertex' and not out['persistent']:
            lo = hi = t0                        # an isolated vertex connection: valid at t0 only
        out.update(word=tuple(word), lo=lo, hi=hi)
        return out

    def track(self, i, d, T0, T1, cap=4000, delta=Fraction(1, 10**9), log=None):
        """pieces covering (T0, T1); T0, T1 rationals (or K-elements). Returns list of dicts (prong output ⊕ the exact
        breakpoint configuration 'bp' at each right end)."""
        K, sg = self.K, self.sgn
        T0, T1 = K.of(T0), K.of(T1)

        def rat_after(x, gap):
            q = Fraction(mpmath.nstr(x.mp(40) + gap, 35))
            while sg(K.of(q) - x) <= 0:
                gap *= 2
                q = Fraction(mpmath.nstr(x.mp(40) + gap, 35))
            return q
        pieces = []
        t = rat_after(T0, mpmath.mpf(10) ** -30)
        left = T0
        while True:
            pr = self.prong(i, d, t, cap)
            if pr['lo'] is not None and sg(pr['lo'] - left) > 0:
                # gap between `left` and this piece: bisect toward `left`
                tt = Fraction(mpmath.nstr((left.mp(40) + K.of(t).mp(40)) / 2, 35))
                if sg(K.of(tt) - left) <= 0:
                    tt = rat_after(left, mpmath.mpf(10) ** -32)
                t = tt
                continue
            pieces.append(pr)
            if log:
                log(pr)
            hi = pr['hi']
            if hi is None or sg(hi - T1) >= 0:
                break
            pr['bp'] = self.prong(i, d, hi, cap)
            left = hi
            gap = mpmath.mpf(delta.numerator) / delta.denominator * max(mpmath.mpf(1), abs(hi.mp(20)))
            t = rat_after(hi, gap)
        return pieces
