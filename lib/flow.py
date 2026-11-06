"""EXACT flow decomposition of a rational-triangle billiard in a direction of the
`pi/(2d)` grid: every component of the unfolded translation surface is certified a CYLINDER or
WITHOUT PERIODIC TRAJECTORY (Boshernitzan), or left UNDETERMINED at the step cap.

A port of the decision procedure of flatsurf's intervalxt (Delecroix--Rueth;
libintervalxt/src/{interval_exchange_transformation,component}.cc), in pure Python over the exact
cyclotomic field below, so that it runs anywhere the rest of `lib/` does.

THE SURFACE.  Triangle `(a, b, c)` with angles `a pi/d, b pi/d, c pi/d` at V0, V1, V2, `d = a+b+c`,
`gcd(a,b,c) = 1`.  Field `Q(zeta_N)`, `N = 4d`, so `E(m) := exp(i m pi/d) = zeta^(2m)` and
`i = zeta^d`.  V0 = 0, V1 = 2 sin(c pi/d), V2 = 2 sin(b pi/d) E(a): all in Z[zeta].  The unfolding
`X` is the 2d copies `g Delta`, `g in D_d` linear (`g = (k, 0): z -> E(2k) z`,
`(k, 1): z -> E(2k) conj z`), side `s` of copy `g` glued to side `s` of copy `g o L_s`
(`L_s` the linear part of the reflection in side `s`).  This is the minimal translation cover
(flatsurf's `minimal_cover("translation")`), marked points INCLUDED: a vertex of angle `pi/2` is a
cone angle `2 pi` point that still cuts the transversal, so it can split a cylinder in two -- a
cylinder COUNT is convention-dependent, a minimal-component count is not.

THE MAP.  Flow direction `omega = zeta^j`.  Transverse coordinate `x_g(P) = 2 cross(omega, g P)`
(exact, real).  Transversal = every edge, taken as the ENTRY side of the copy the flow enters
through, laid end to end on `[0, L)`.  One crossing of one copy is `T`: exact translations, and
`T` is an IET (checked: its images tile `[0, L)` EXACTLY, `build_iet` asserts it).  The geometry
is carried as MONOMIAL sums in `zeta` (every vertex is two roots of unity, `D_d` permutes roots of
unity), reduced to the power basis only for comparison -- no dense multiply anywhere.

RETURN TIMES.  Each label carries `h` = its return time in copy crossings, updated through every
Zorich step, Dehn twist and merge, so a cylinder reports its exact PERIOD and a component its Kac
MASS `sum h*lam / L`; `sum h*lam = L` over all labels is asserted EXACTLY on every run.

THE INDUCTION (intervalxt, faithfully).  Left Zorich-accelerated Rauzy induction, top then bottom
(`_zorich`, with `subtractRepeated`'s Dehn twist); after each round: `_reduce` splits a reducible
permutation (a separating connection / invariant subset); equal first lengths = a saddle
connection, merged; a single interval = a CYLINDER.  Every `BOSH_EVERY` rounds, Boshernitzan
(gated, as in intervalxt, on SAF != 0, where alone it can succeed -- by an O(nD) exact random-
projection test, `saf_nonzero`, since the dense test is O(nD^2); the LP has a `LP_SECONDS` limit,
a stalled solve being 'no certificate yet'): a periodic orbit through a component visits interval `j` `a_j >= 0` times with
`sum a_j t_j = 0` (`t_j` the translations, as Q-vectors in the power basis) -- so a Q-linear
functional `y` with `y . t_j > 0` for ALL `j` (Gordan) PROVES there is no periodic trajectory.
`y` is found by an LP and VERIFIED IN EXACT INTEGERS; the certificate is returned.
SAF = 0 components (intervalxt's self-similarity loop is NOT ported) end UNDETERMINED.

SIGNS.  Exact zero test (reduced coefficient vector); a nonzero sign from a float dot product with
a rigorous error bound, else an exact fixed-point integer dot product against cosines rounded to
2^-256 (error <= sum |a_k|, precision quadrupled on a tie), else a certified `mpmath.iv` enclosure
(as `cyclofield.py`, which is the cross-check reference; `_selftest`).

ARITHMETIC.  The drivers (`decompose`, `crossings`, `n_flow`, `n_flow_pos`, `m_par_flow`,
`build_iet`) compute in `CycNP` -- the same field with elements as numpy int64 vectors under an
overflow guard (object arrays of Python ints past 2^62, demoted back when they fit) -- and hand out
`Cyc` tuples; the middle sign stage is a vectorised exact 128-bit fixed-point dot product.  The
arithmetic is an implementation choice only: the outputs are identical to the tuple arithmetic.

RUN CONTROL.  `RunControl` -> `ctl=` of `n_flow` / `crossings` / `Decomposer.run`: progress
lines, atomic checkpoints, a stop file, resume (for the long rows, e.g. `3/385`).

Entry: `decompose(abc, j, cap)` -> dict (components with kind / size / transverse measure /
certificate).  `perp_dirs(abc)`, `side_dirs(abc)` give `j` for the directions perpendicular /
parallel to each side.
"""
from fractions import Fraction
import math
from operator import add as _op_add, sub as _op_sub

import numpy as np
import sympy as sp
from mpmath import iv

BOSH_EVERY = 16
SAF_PRIMES = (8388593, 8388587, 8388581)   # primes < 2^23 for `saf_nonzero`
LP_SECONDS = 20
AUTOSIM = True          # intervalxt's SAF-0 self-similarity loop (SimilarityTracker), see `autosimilar`
AUTOSIM_AFTER = 256     # arm `autosimilar` only after this many steps with no merge/reduce (0 = s635 behaviour)


# ---------------------------------------------------------------------------------------------
# exact field Z[zeta_N], power basis, integer coefficient tuples
# ---------------------------------------------------------------------------------------------
_CYC_TABLES = {}        # N -> (D, tab, C, fixed-point tables): immutable, shared by every Cyc(N)


class Cyc:
    def __init__(self, N):
        self.N = N
        self.nsign_iv = 0
        self.nsign_fx = 0
        hit = _CYC_TABLES.get(N)                     # the tables are immutable -- shared per N
        if hit is not None:
            self.D, self.tab, self.C, self._fx = hit
            self.zero = (0,) * self.D
            self.one = self.tab[0]
            return
        z = sp.symbols('z')
        phi = [int(c) for c in reversed(sp.Poly(sp.cyclotomic_poly(N, z), z).all_coeffs())]
        self.D = D = len(phi) - 1
        # zeta^k in the power basis, k = 0..N-1
        tab = []
        cur = [0] * D
        cur[0] = 1
        for k in range(N):
            tab.append(tuple(cur))
            # multiply by zeta: shift, reduce zeta^D = -sum phi_i zeta^i
            top = cur[-1]
            cur = [0] + cur[:-1]
            if top:
                for i in range(D):
                    cur[i] -= top * phi[i]
        self.tab = tab
        self.C = [math.cos(2 * math.pi * k / N) for k in range(D)]
        self.zero = (0,) * D
        self.one = tab[0]
        self._fx = {}                                # bits -> fixed-point cosines
        _CYC_TABLES[N] = (self.D, self.tab, self.C, self._fx)

    def z(self, k):
        return self.tab[k % self.N]

    def add(self, u, v):
        return tuple(map(_op_add, u, v))

    def sub(self, u, v):
        return tuple(map(_op_sub, u, v))

    def maxabs(self, u):
        return max(map(abs, u))

    def neg(self, u):
        return tuple(-a for a in u)

    def smul(self, n, u):
        return tuple(n * a for a in u)

    def mul(self, u, v):
        D, N = self.D, self.N
        acc = [0] * D
        for i, a in enumerate(u):
            if not a:
                continue
            for j, b in enumerate(v):
                if b:
                    t = self.tab[(i + j) % N]
                    ab = a * b
                    for k in range(D):
                        if t[k]:
                            acc[k] += ab * t[k]
        return tuple(acc)

    def conj(self, u):
        acc = [0] * self.D
        for k, a in enumerate(u):
            if a:
                t = self.tab[(-k) % self.N]
                for i in range(self.D):
                    if t[i]:
                        acc[i] += a * t[i]
        return tuple(acc)

    def cross2(self, u, v):
        """2 * cross(u, v) = 2 Im(conj(u) v) = -i (conj(u) v - u conj(v)); real."""
        w = self.sub(self.mul(self.conj(u), v), self.mul(u, self.conj(v)))
        return self.mul(self.neg(self.z(self.N // 4)), w)

    def fval(self, u):
        return sum(a * c for a, c in zip(u, self.C) if a)

    def mpval(self, u, dps=None):
        """High-precision point value (an ESTIMATE; never used to decide)."""
        from mpmath import mp, mpf, cos, pi
        m = max(1, max(map(abs, u)))
        with mp.workdps(dps or (40 + len(str(m)))):
            return sum(mpf(int(a)) * cos(2 * pi * k / self.N) for k, a in enumerate(u) if a)

    def sign(self, u):
        """Exact sign of a REAL element."""
        s = 0.0
        m = 0
        for a, c in zip(u, self.C):
            if a:
                s += float(a) * c
                m += abs(a)
        if m == 0:
            return 0
        bound = float(m) * (self.D + 4) * 4.5e-16
        if s > bound:
            return 1
        if s < -bound:
            return -1
        return self._sign_fx(u, m)

    def _fx_table(self, bits):
        """round(cos(2 pi k / N) * 2^bits) for k < D, each within 1 of the true scaled value
        (computed at bits + 64 binary digits, then rounded)."""
        t = self._fx.get(bits)
        if t is None:
            from mpmath import mp, mpf, cos, pi, nint
            with mp.workprec(bits + 64):
                sc = mpf(2) ** bits
                t = [int(nint(cos(2 * pi * k / self.N) * sc)) for k in range(self.D)]
            self._fx[bits] = t
        return t

    def _sign_fx(self, u, m):
        """Exact sign by a fixed-point integer dot product (s761; replaced a per-call interval
        evaluation of every cosine, 70% of n_flow's time at Q = 151).  With C_k the table above,
        |2^bits * value - sum a_k C_k| <= sum |a_k| = m, so |sum| > m decides the sign.  Doubles the
        precision on a near-tie; the interval enclosure stays as the last resort."""
        self.nsign_fx += 1
        bits = 256
        while bits <= 4096:
            S = sum(a * c for a, c in zip(u, self._fx_table(bits)) if a)
            if S > m:
                return 1
            if S < -m:
                return -1
            bits *= 4
        return self._sign_iv(u)

    def _sign_iv(self, u):
        self.nsign_iv += 1
        dps = 60
        while True:
            old = iv.dps
            iv.dps = dps
            try:
                tot = iv.mpf(0)
                for k, a in enumerate(u):
                    if a:
                        tot += iv.mpf(int(a)) * iv.cos(2 * iv.pi * k / self.N)
            finally:
                iv.dps = old
            if tot.a > 0:
                return 1
            if tot.b < 0:
                return -1
            dps *= 2
            if dps > 8000:
                raise RuntimeError('enclosure did not separate a nonzero element from 0')

    def cmp(self, u, v):
        return self.sign(self.sub(u, v))


_I64 = np.dtype(np.int64)
_CYCNP_TABLES = {}      # N -> CycNP's numpy tables and caches (immutable / deterministic), shared
_LIM = 1 << 62          # CycNP invariant: an int64 element has every |coeff| < 2^62, so u +- v cannot wrap


class CycNP(Cyc):
    """The SAME field as `Cyc`, elements as numpy vectors: int64 while every |coeff| < 2^62
    (the common case -- 99.9% of the ops at 144/233), else an object array of Python ints, demoted
    back to int64 as soon as it fits again.  Linear ops are one numpy op plus an overflow guard
    (u +- v of two guarded int64 vectors cannot wrap; the result is re-guarded, and redone over Z on
    a miss; `smul` pre-checks |n| * max|u|), so every value is EXACT, as in `Cyc`; `sign` is the
    same three-stage decision (float with the same rigorous bound, fixed-point, interval) -- an exact
    sign, so every decision and output is identical to `Cyc`'s (`_selftest`, s869 regress).
    INTERNAL to the drivers of this module (`decompose`, `crossings`, `n_flow`, `n_flow_pos`,
    `m_par_flow`, `build_iet`), which hand out tuples exactly as before; `Cyc` keeps the tuple API
    for every probe that builds its own `Decomposer`.  Elements are never mutated in place (results
    may share an operand's array)."""

    def __init__(self, N, base=None):
        if base is None:
            base = Cyc(N)
        self.N, self.D, self.tab, self.C = base.N, base.D, base.tab, base.C
        self.nsign_iv = self.nsign_fx = 0
        self._fx = base._fx
        sh = _CYCNP_TABLES.get(self.N)               # immutable numpy tables + caches, shared per N
        if sh is None:
            taba = np.array(self.tab, dtype=np.int64)
            taba.flags.writeable = False
            assert int(np.abs(taba).max()) < 2 ** 20
            sh = _CYCNP_TABLES[self.N] = {'Ca': np.array(self.C, dtype=np.float64), 'taba': taba,
                                          'mpcos': {}, 'limbs': None}
        self._sh = sh
        self.Ca, self.taba = sh['Ca'], sh['taba']
        self.zero = self.taba[0] * 0
        self.zero.flags.writeable = False
        self.one = self.taba[0]
        self.bound_c = (self.D + 4) * 4.5e-16
        self._mpcos = sh['mpcos']                    # mp.prec -> {k: cos(2 pi k / N)} (mpval)

    # conversions ------------------------------------------------------------------------------
    def vec(self, t):
        """Any coefficient sequence (tuple / list / array) -> a guarded element."""
        if isinstance(t, np.ndarray) and t.dtype is _I64:
            return t if int(np.abs(t).max()) < _LIM else self._norm(t.astype(object))
        try:
            a = np.array(t, dtype=np.int64)
        except OverflowError:
            return self._norm(np.array([int(x) for x in t], dtype=object))
        if int(np.abs(a).max()) < _LIM:
            return a
        return np.array([int(x) for x in t], dtype=object)

    @staticmethod
    def tup(u):
        """Element -> the tuple of Python ints `Cyc` would hold."""
        return tuple(u.tolist()) if isinstance(u, np.ndarray) else tuple(u)

    @staticmethod
    def _norm(r):
        """Object array -> int64 when every |coeff| < 2^62."""
        if max(map(abs, r)) < _LIM:
            return r.astype(np.int64)
        return r

    def red(self, mono):
        """Monomial dict -> element (`Triangle.red`, vectorised: coefficients < 2^20 by construction,
        checked, so each sum of < 2^21 terms stays far inside int64)."""
        if not mono:
            return self.zero
        es = list(mono.keys())
        cs = list(mono.values())
        if len(cs) < 2 ** 20 and max(map(abs, cs)) < 2 ** 20:
            return self.vec(np.array(cs, dtype=np.int64) @ self.taba[es])
        acc = [0] * self.D
        for e, c in mono.items():
            for i, x in enumerate(self.tab[e]):
                if x:
                    acc[i] += c * x
        return self.vec(acc)

    # field ops --------------------------------------------------------------------------------
    def z(self, k):
        return self.taba[k % self.N]

    def add(self, u, v):
        r = u + v
        if r.dtype is _I64:
            if np.abs(r).max() < _LIM:
                return r
            return self._norm(u.astype(object) + v.astype(object))
        return self._norm(r)

    def sub(self, u, v):
        r = u - v
        if r.dtype is _I64:
            if np.abs(r).max() < _LIM:
                return r
            return self._norm(u.astype(object) - v.astype(object))
        return self._norm(r)

    def neg(self, u):
        return -u                                    # |coeff| < 2^62 is symmetric; object stays exact

    def smul(self, n, u):
        n = int(n)
        if u.dtype is _I64:
            if n == 0:
                return self.zero
            if n == 1:
                return u
            if abs(n) < _LIM and abs(n) * int(np.abs(u).max()) < _LIM:
                return u * n
            return self._norm(u.astype(object) * n)
        return self._norm(u * n)

    def maxabs(self, u):
        return int(np.abs(u).max())

    def mul(self, u, v):
        return self.vec(Cyc.mul(self, u.tolist(), v.tolist()))

    def conj(self, u):
        return self.vec(Cyc.conj(self, u.tolist()))

    def fval(self, u):
        return Cyc.fval(self, u.tolist() if isinstance(u, np.ndarray) else u)   # same sum order as Cyc

    def mpval(self, u, dps=None):
        """`Cyc.mpval`, bit-identical: the same terms summed in the same order at the same precision;
        only cos(2 pi k / N) is cached per (precision, k) -- a deterministic function of both."""
        from mpmath import mp, mpf, cos, pi
        from mpmath.libmp import mpf_add, mpf_mul, mpf_pos, from_int, fzero
        ul = u.tolist() if isinstance(u, np.ndarray) else list(u)
        m = max(1, max(map(abs, ul)))
        with mp.workdps(dps or (40 + len(str(m)))):
            prec, rnd = mp._prec_rounding
            cs = self._mpcos.get(prec)
            if cs is None:
                cs = self._mpcos[prec] = {}
            # the raw libmp calls that `sum(mpf(int(a)) * cos(...) ...)` makes: mpf(int) = mpf_pos(from_int),
            # mpf * mpf = mpf_mul, 0 + x = mpf_add(x, fzero), x + y = mpf_add(x, y) -- all at (prec, rnd)
            tot = None
            for k, a in enumerate(ul):
                if a:
                    c = cs.get(k)
                    if c is None:
                        c = cs[k] = cos(2 * pi * k / self.N)._mpf_
                    t = mpf_mul(mpf_pos(from_int(int(a)), prec, rnd), c, prec, rnd)
                    tot = mpf_add(t, fzero, prec, rnd) if tot is None else mpf_add(tot, t, prec, rnd)
            return 0 if tot is None else mpf(tot)

    def quot_est(self, a, b):
        """An ESTIMATE of floor(a / b) for `Decomposer._floordiv`, which corrects it EXACTLY (so the
        quotient, hence every output, does not depend on it): a float dot product when b's value
        clears its error by 2^20, else the fixed-point table (b's error <= sum |b_k| ulps, raised to
        clear by 2^40), else None (the caller's mpmath path)."""
        fa = float(a.astype(np.float64) @ self.Ca)
        fbv = b.astype(np.float64)
        fb = float(fbv @ self.Ca)
        if abs(fb) > float(np.abs(fbv).sum()) * self.bound_c * 2.0 ** 20:
            return max(int(fa // fb) - 1, 0)
        al, bl = a.tolist(), b.tolist()
        ma, mb = sum(map(abs, al)), sum(map(abs, bl))
        bits = 256
        while bits <= 4096:
            C = self._fx_table(bits)
            Sb = sum(x * c for x, c in zip(bl, C) if x)
            if abs(Sb) > (max(ma, mb) + 1) << 40:
                Sa = sum(x * c for x, c in zip(al, C) if x)
                return max(Sa // Sb - 1, 0)
            bits *= 4
        return None

    def sign(self, u):
        """`Cyc.sign`: float dot product (numpy; its rounding error is within the same bound
        m (D+4) 4.5e-16, m = sum |a_k|), else the exact fixed-point / interval stages on Python ints."""
        f = u.astype(np.float64)
        m = float(np.abs(f).sum())
        if m == 0:
            return 0
        s = float(f @ self.Ca)
        bound = m * self.bound_c
        if s > bound:
            return 1
        if s < -bound:
            return -1
        if u.dtype is _I64:
            # an upper bound of the exact sum |u_k|: m is the float sum of the (rounded) |u_k|, relative
            # error <= (log2 D + 1) 2^-53 << 2^-40
            r = self._sign_fx_np(u, int(m * (1.0 + 2.0 ** -40)) + 1)
            if r:
                return r
        ul = u.tolist()
        return self._sign_fx(ul, sum(map(abs, ul)))

    FXNP_BITS = 128     # the vectorised fixed-point stage; on record every sign the float
                        # stage leaves has |value| > 2^-64 sum |a_k| (89/144), so 128 bits decide them

    def _fx_limbs(self):
        """`_fx_table(FXNP_BITS)` as a D x L float64 matrix of 20-bit integer limbs
        (C_k = sum_j limb_kj 2^(20 j), limbs in [0, 2^20) except the signed top one, in {-1, 0, 1})."""
        t = self._sh['limbs']
        if t is None:
            C = self._fx_table(self.FXNP_BITS)
            L = (self.FXNP_BITS + 2) // 20 + 2
            rows = []
            for c in C:
                row = []
                for j in range(L - 1):
                    row.append(c & 0xFFFFF)
                    c >>= 20
                assert abs(c) <= 1
                row.append(c)                        # signed top limb
                rows.append(row)
            t = self._sh['limbs'] = np.array(rows, dtype=np.float64)
            t.flags.writeable = False
        return t

    _USH = np.array([[0], [21], [42]], dtype=np.int64)
    _UMASK = np.array([[0x1FFFFF], [0x1FFFFF], [-1]], dtype=np.int64)

    def _sign_fx_np(self, u, m):
        """A fixed-point stage like `_sign_fx`'s, at `FXNP_BITS`, vectorised and EXACT: u (int64,
        |u_k| < 2^63) is cut into three limbs (21, 21 bits unsigned, then the signed rest, |.| <= 2^21)
        and multiplied by the 20-bit limbs of the table as a float64 BLAS product: every product is an
        integer < 2^41 and every partial sum of the D < 2^11 terms an integer < 2^52, so the float
        arithmetic is exact in any summation order; recombined over Z it is S = sum_k u_k C_k exactly.
        Each C_k is within 1 of 2^bits cos(2 pi k / N), so |2^bits value - S| <= sum |u_k| <= `m` (any
        upper bound) and |S| > m decides the sign -- `_sign_fx`'s argument.  0 = undecided (the caller
        then runs `_sign_fx` itself: 256, 1024, 4096 bits, the enclosure)."""
        if self.D >= 2048:
            return 0
        U = ((u >> self._USH) & self._UMASK).astype(np.float64)     # 3 x D, u = sum_i U_i 2^(21 i)
        r0, r1, r2 = (U @ self._fx_limbs()).astype(np.int64).tolist()   # 3 x L, exact
        S = 0
        for j in range(len(r0) - 1, -1, -1):        # S = sum_ij r_i[j] 2^(21 i + 20 j), Horner in j
            S = (S << 20) + r0[j] + (r1[j] << 21) + (r2[j] << 42)
        self.nsign_fx += 1
        if S > m:
            return 1
        if S < -m:
            return -1
        return 0

    def _sign_iv(self, u):
        return Cyc._sign_iv(self, u.tolist() if isinstance(u, np.ndarray) else u)

    def cmp(self, u, v):
        if u.dtype is _I64 and v.dtype is _I64:
            return self.sign(u - v)                  # |u - v| < 2^63: exact in int64, only signed
        return self.sign(self.sub(u, v))


def _fast(T):
    """The `CycNP` twin of `T.F` (one per Triangle, so sign counters stay per run as with `T.F`)."""
    FV = getattr(T, '_fv', None)
    if FV is None:
        FV = T._fv = CycNP(T.N, T.F)
    return FV


# ---------------------------------------------------------------------------------------------
# the unfolded surface and its one-crossing IET
# ---------------------------------------------------------------------------------------------
class Triangle:
    """Geometry as MONOMIAL sums {exponent mod N: coeff} in zeta: every vertex is two roots of unity
    and every g in D_d maps a root of unity to one, so nothing here needs a dense multiply; an
    element is reduced to the power basis (`red`) only when compared or accumulated."""

    def __init__(self, abc):
        a, b, c = abc
        assert math.gcd(math.gcd(a, b), c) == 1
        self.abc = abc
        self.d = d = a + b + c
        self.N = N = 4 * d
        self.F = Cyc(N)
        mi = 3 * d                                   # -i = zeta^(3d)
        S = lambda m, sh: {(2 * m + sh) % N: 1, (-2 * m + sh) % N: -1}   # 2i sin(m pi/d) zeta^sh
        self.V = [{}, S(c, mi), S(b, mi + 2 * a)]
        # sides: s0 = V0V1, s1 = V1V2, s2 = V2V0; line direction E(m_s)
        self.sides = [(0, 1), (1, 2), (2, 0)]
        self.m = [0, (d - b) % d, a % d]
        # affine reflection in side s: z -> E(2 m_s) conj(z) + t_s, t_s = P - L_s P for P on line
        self.t = []
        for s, (i, _) in enumerate(self.sides):
            P = self.V[i]
            LP = mshift(mconj(P, N), 4 * self.m[s], N)
            self.t.append(madd(P, mneg(LP)))

    def g_apply(self, g, p):
        k, e = g
        q = mconj(p, self.N) if e else p
        return mshift(q, 4 * k, self.N)

    def g_across(self, g, s):
        k, e = g
        m = self.m[s]
        return ((k + m) % self.d, 1) if e == 0 else ((k - m) % self.d, 0)

    def red(self, mono):
        F = self.F
        acc = [0] * F.D
        for e, c in mono.items():
            if c:
                t = F.tab[e]
                for i in range(F.D):
                    if t[i]:
                        acc[i] += c * t[i]
        return tuple(acc)

    def cross2(self, u, v):
        """2 cross(u, v) = -i (conj(u) v - u conj(v)), monomial in, monomial out."""
        N = self.N
        w = madd(mmul(mconj(u, N), v, N), mneg(mmul(u, mconj(v, N), N)))
        return mshift(w, 3 * self.d, N)


def madd(u, v):
    w = dict(u)
    for e, c in v.items():
        w[e] = w.get(e, 0) + c
    return {e: c for e, c in w.items() if c}


def mneg(u):
    return {e: -c for e, c in u.items()}


def mshift(u, k, N):
    return {(e + k) % N: c for e, c in u.items()}


def mconj(u, N):
    return {(-e) % N: c for e, c in u.items()}


def mmul(u, v, N):
    w = {}
    for e1, c1 in u.items():
        for e2, c2 in v.items():
            e = (e1 + e2) % N
            w[e] = w.get(e, 0) + c1 * c2
    return {e: c for e, c in w.items() if c}


def perp_dirs(abc):
    """j with omega = zeta^j perpendicular to side s (s = 0, 1, 2): i E(m_s)."""
    T = Triangle(abc)
    return [(T.d + 2 * m) % (4 * T.d) for m in T.m]


def side_dirs(abc):
    T = Triangle(abc)
    return [(2 * m) % (4 * T.d) for m in T.m]


def build_iet(T, j):
    """The one-crossing IET of direction omega = zeta^j (or, `j` a monomial dict {exponent: coeff},
    omega = that element of Z[zeta] -- any direction with exact holonomy, e.g. a saddle connection's).
    Returns (labels data, lengths, top, bottom, meta).  Asserts the image intervals tile [0, L)
    exactly.  Field elements are `Cyc` tuples (computed in `CycNP`, `_build_iet`)."""
    pieces, lengths, top, bot, meta = _build_iet(T, j)
    tu = CycNP.tup
    pieces = [(tu(a), tu(w), tu(b), k1, k2) for a, w, b, k1, k2 in pieces]
    lengths = {p: pieces[p][1] for p in lengths}
    meta['L'] = tu(meta['L'])
    return pieces, lengths, top, bot, meta


def _build_iet(T, j):
    """`build_iet` with the elements native to `_fast(T)` (numpy vectors)."""
    F, d, N = _fast(T), T.d, T.N
    om = {e % N: c for e, c in j.items()} if isinstance(j, dict) else {j % N: 1}
    X = lambda p: F.red(T.cross2(om, p))
    copies = [(k, e) for e in (0, 1) for k in range(d)]
    # per copy: vertex x-coords (unfolded), classify sides entry / exit / parallel
    info = {}
    for g in copies:
        W = [T.g_apply(g, v) for v in T.V]
        xs = [X(w) for w in W]
        ent, ext = [], []
        for s, (i, k2) in enumerate(T.sides):
            kk = 3 - i - k2
            e = madd(W[k2], mneg(W[i]))
            c_om = F.sign(F.red(T.cross2(e, om)))
            if c_om == 0:
                continue                             # side parallel to the flow
            c_in = F.sign(F.red(T.cross2(e, madd(W[kk], mneg(W[i])))))
            lo, hi = (xs[i], xs[k2]) if F.cmp(xs[i], xs[k2]) < 0 else (xs[k2], xs[i])
            (ent if c_om == c_in else ext).append((s, lo, hi))
        info[g] = (xs, ent, ext)
    # entry segments laid end to end
    order = [(g, s, lo, hi) for g in copies for (s, lo, hi) in info[g][1]]
    off = {}
    pos = F.zero
    for g, s, lo, hi in order:
        off[(g, s)] = (pos, lo, hi)
        pos = F.add(pos, F.sub(hi, lo))
    L = pos
    pieces = []                                      # (dom_start, length, img_start, (g,s), (g',s'))
    for g, s, lo, hi in order:
        O, _, _ = off[(g, s)]
        for (s2, lo2, hi2) in info[g][2]:
            u = lo if F.cmp(lo, lo2) >= 0 else lo2
            w = hi if F.cmp(hi, hi2) <= 0 else hi2
            if F.cmp(w, u) <= 0:
                continue
            g2 = T.g_across(g, s2)
            key2 = (g2, s2)
            assert key2 in off, ('exit side is not the neighbour entry', g, s2, g2)
            O2, lo2e, _ = off[key2]
            shift = X(T.g_apply(g, T.t[s2]))         # x_{g'}(P) = x_g(P) - x_g(t_s)
            dom = F.add(O, F.sub(u, lo))
            img = F.add(O2, F.sub(F.sub(u, shift), lo2e))
            pieces.append((dom, F.sub(w, u), img, (g, s), key2))
    n = len(pieces)
    import functools
    top = sorted(range(n), key=functools.cmp_to_key(lambda p, q: F.cmp(pieces[p][0], pieces[q][0])))
    bot = sorted(range(n), key=functools.cmp_to_key(lambda p, q: F.cmp(pieces[p][2], pieces[q][2])))
    # exact tiling checks
    for seq, idx in ((top, 0), (bot, 2)):
        cur = F.zero
        for p in seq:
            assert F.cmp(pieces[p][idx], cur) == 0, 'intervals do not tile [0, L)'
            cur = F.add(cur, pieces[p][1])
        assert F.cmp(cur, L) == 0
    lengths = {p: pieces[p][1] for p in range(n)}
    return pieces, lengths, top, bot, {'L': L, 'ncopies': len(copies), 'nentry': len(order)}


# ---------------------------------------------------------------------------------------------
# intervalxt's induction
# ---------------------------------------------------------------------------------------------
class Comp:
    __slots__ = ('top', 'bot', 'kind', 'cert', 'steps', 'trk', 'saf0')

    def __init__(self, top, bot):
        self.top, self.bot = list(top), list(bot)
        self.kind, self.cert, self.steps = None, None, 0
        self.trk, self.saf0 = None, None             # autosimilar state; reset on merge / reduce


class Decomposer:
    def __init__(self, F, lengths):
        self.F = F
        self.lam = dict(lengths)
        self.h = {l: 1 for l in lengths}            # return time (copy crossings) of each label

    def _zorich(self, top, bot):
        """One left Zorich step with `top`'s first interval as the (possible) winner.  Mutates
        `bot` and lengths.  Returns True iff a connection is visible at the left end."""
        F, lam = self.F, self.lam
        t = top[0]
        if bot[0] == t:
            return True
        stack, S = [], F.zero
        k = 0
        while True:
            b = bot[k]
            if b == t:
                # Dehn twist: all of bot[:k] (total S < lam[t]) precedes t
                lt = lam[t]
                q = self._floordiv(lt, S)
                lt = F.sub(lt, F.smul(q, S))
                if F.sign(lt) == 0:
                    lt = F.add(lt, S)
                    q -= 1
                nmove = len(stack)                  # default stop = last stack label
                npart = 0                           # labels subtracted in the partial twist
                for idx, lab in enumerate(stack):
                    if F.cmp(lam[lab], lt) >= 0:
                        nmove = idx if idx > 0 else len(stack)
                        npart = idx
                        break
                    lt = F.sub(lt, lam[lab])
                else:
                    raise RuntimeError('floor division inconsistent')
                lam[t] = lt
                # return times: every stack label composes t q times, the partial ones once more
                ht = self.h[t]
                for idx, lab in enumerate(stack):
                    self.h[lab] += (q + (1 if idx < npart else 0)) * ht
                k = nmove
                break
            S2 = F.add(S, lam[b])
            if F.cmp(S2, lam[t]) >= 0:
                lam[t] = F.sub(lam[t], S)
                for lab in stack:
                    self.h[lab] += self.h[t]
                break
            stack.append(b)
            S = S2
            k += 1
        if k:
            it = bot.index(t)
            bot[:] = bot[k:it] + bot[:k] + bot[it:]
        return F.cmp(lam[top[0]], lam[bot[0]]) == 0

    def _floordiv(self, a, b):
        F = self.F
        qe = getattr(F, 'quot_est', None)            # CycNP's estimate (float / fixed point);
        q = qe(a, b) if qe is not None else None     # the loops below make q exact whatever it is
        if q is not None:
            pass
        elif abs(fb := F.fval(b)) > (m := max(max(map(abs, a)), max(map(abs, b)))) * F.D * 1e-12:
            q = max(int(F.fval(a) // fb) - 1, 0)
        else:                                       # cancellation: estimate at high precision,
            from mpmath import floor                # raised until b's estimate clears zero
            dps = 40 + len(str(m))
            while True:
                vb = F.mpval(b, dps)
                if abs(vb) > m * F.D * 10.0 ** (-dps + 5):
                    break
                dps *= 2
                if dps > 8000:
                    raise RuntimeError('floordiv: divisor indistinguishable from 0')
            from mpmath import mp
            with mp.workdps(dps):
                q = max(int(floor(F.mpval(a, dps) / vb)) - 1, 0)
        while F.sign(F.sub(a, F.smul(q + 1, b))) >= 0:
            q += 1
        while q > 0 and F.sign(F.sub(a, F.smul(q, b))) < 0:
            q -= 1
        return q

    @staticmethod
    def _reduce(c):
        """Shortest prefix of top and bot with equal label sets; split the rest off.
        O(n) -- `odd` = the labels seen in exactly one of the two prefixes (each label occurs
        once in top and once in bot), so the prefix sets are equal iff it is empty; the s623 loop
        compared the two sets at every k, O(n^2).  The whole prefix (k = n - 1) is not a split.
        For int labels (always, from `build_iet`) it is vectorised: top[:k+1] and bot[:k+1] (k+1
        distinct labels each) are the same set iff max_{i<=k} (position of top[i] in bot) == k."""
        top, bot = c.top, c.bot
        n = len(top)
        if n > 16 and type(top[0]) is int:
            t = np.array(top)
            inv = np.empty(int(t.max()) + 1, dtype=np.int64)
            inv[np.array(bot)] = np.arange(n)
            hit = np.flatnonzero(np.maximum.accumulate(inv[t[:-1]]) == np.arange(n - 1))
            if not len(hit):
                return None
            k = int(hit[0])
            rest = Comp(top[k + 1:], bot[k + 1:])
            c.top, c.bot = top[:k + 1], bot[:k + 1]
            return rest
        odd = set()
        for k in range(len(top) - 1):
            a, b = top[k], bot[k]
            if a != b:
                if a in odd:
                    odd.remove(a)
                else:
                    odd.add(a)
                if b in odd:
                    odd.remove(b)
                else:
                    odd.add(b)
            if not odd:
                rest = Comp(top[k + 1:], bot[k + 1:])
                c.top, c.bot = top[:k + 1], bot[:k + 1]
                return rest
        return None

    def translations(self, c):
        F, lam = self.F, self.lam
        pt, pb = {}, {}
        cur = F.zero
        for l in c.top:
            pt[l] = cur
            cur = F.add(cur, lam[l])
        cur = F.zero
        for l in c.bot:
            pb[l] = cur
            cur = F.add(cur, lam[l])
        return {l: F.sub(pb[l], pt[l]) for l in c.top}

    def saf_zero(self, c):
        tr = self.translations(c)
        lam = np.array([[int(x) for x in self.lam[l]] for l in c.top], dtype=object)   # Python ints:
        t = np.array([[int(x) for x in tr[l]] for l in c.top], dtype=object)           # no int64 wrap
        W = lam.T.dot(t)                            # sum_j lam_j (x) t_j
        A = W - W.T
        return all(x == 0 for x in A.flat)

    def saf_nonzero(self, c, trials=3):
        """SAF(c) != 0, PROVED when True: the projection SAF(u, v) = sum_j (u.lam_j)(v.t_j) - (v.lam_j)(u.t_j)
        onto random functionals u, v is computed EXACTLY MOD A PRIME p (int64 numpy, s761 -- the
        object-dtype version over Z was 90% of n_flow's time once signs were fast), and nonzero mod p
        means nonzero over Z.  False means every projection vanished -- SAF = 0 up to a miss
        (probability <= 2/p ~ 2e-7 per trial for SAF != 0 not divisible by p; the trials use
        different primes).  O(n D) per trial (the dense `saf_zero` is O(n D^2)).  Used only to
        decide whether to ATTEMPT Boshernitzan, which can succeed only when SAF != 0 (intervalxt's
        gate; 122/122 certified no-periodic components on record have SAF != 0, s626), and to arm
        `autosimilar`, which certifies by exact proportionality: a miss costs a certificate, never a
        wrong one."""
        pos = {l: i for i, l in enumerate(c.top)}
        bidx = np.array([pos[l] for l in c.bot])
        rng = np.random.default_rng(len(c.top) * 1000003 + c.steps)
        D = self.F.D
        assert D < 2 ** 17                          # products < 2^46, D-term sums < 2^63
        cached = D < 2048 and isinstance(self.F, CycNP) and type(c.top[0]) is int
        Lz = big = None
        for t in range(trials):
            p = SAF_PRIMES[t % len(SAF_PRIMES)]
            if cached:
                Lp, sel = self._saf_rows(c.top, p)  # float64 residues mod p: rows Lp[sel]
            else:
                sel = None
                if Lz is None and big is None:
                    rows = [self.lam[l] for l in c.top]
                    try:
                        Lz = np.array(rows, dtype=np.int64)
                    except OverflowError:
                        big = np.array(rows, dtype=object)
                Lp = (Lz - (Lz // p) * p) if big is None else (big % p).astype(np.int64)   # = Lz % p
            u = rng.integers(0, p, D)
            v = rng.integers(0, p, D)
            # PROJECT FIRST.  t_j = (bottom-order prefix sum) - (top-order prefix sum) of the
            # lambda rows is linear, so T @ u is the same prefix-sum difference of the n-vector Lp @ u:
            # the same residues mod p as the s761 dense version (prefix sums of the n x D matrix),
            # for 2 matvecs instead of 4 plus 2 n x D cumsums.  Entries < 2^23 * n: no wrap.
            if D < 2048:    # exact float64 BLAS: 12-bit halves of u, v; products < 2^35, D-sums < 2^46
                M = Lp.astype(np.float64, copy=False) @ np.stack(
                    (u & 4095, u >> 12, v & 4095, v >> 12), axis=1).astype(np.float64)
                M = M.astype(np.int64) if sel is None else M[sel].astype(np.int64)
                Luv = np.stack((M[:, 0] + (M[:, 1] << 12), M[:, 2] + (M[:, 3] << 12)), axis=1) % p
            else:
                Luv = (Lp @ np.stack((u, v), axis=1)) % p     # n x 2: (Lp @ u, Lp @ v) mod p
            pt = np.cumsum(Luv, axis=0) - Luv               # exclusive prefix sums, top order
            pbb = np.cumsum(Luv[bidx], axis=0) - Luv[bidx]
            pb = np.empty_like(pbb)
            pb[bidx] = pbb                                  # back to top order
            Tuv = (pb - pt) % p
            Lu, Lv, Tu, Tv = Luv[:, 0], Luv[:, 1], Tuv[:, 0], Tuv[:, 1]
            if (((Lu * Tv) % p - (Lv * Tu) % p).sum()) % p != 0:
                return True
        return False

    def _saf_rows(self, labels, p):
        """The rows lam[l] mod p (l in `labels`) as an n x D float64 matrix (exact: residues < 2^23), from a
        per-prime cache of every label's reduced row, refreshed only where lam[l] is a NEW object (a
        CycNP element is never mutated in place, and `_zorich` changes one label's length per step) --
        so a call costs a gather, not the n x D reduction mod p.  Returns (matrix, selector): the
        rows are matrix[selector] (selector None = all, in order); a large component projects the
        whole cache and selects after."""
        cache = self.__dict__.setdefault('_safc', {})
        ent = cache.get(p)
        if ent is None:
            n0 = max(self.lam) + 1
            ent = cache[p] = (np.zeros((n0, self.F.D)), [None] * n0)
        G, ref = ent
        lam = self.lam
        for l in labels:
            x = lam[l]
            if ref[l] is not x:
                if x.dtype is _I64:
                    G[l] = x - (x // p) * p
                else:
                    G[l] = (x % p).astype(np.int64)
                ref[l] = x
        if 2 * len(labels) >= len(ref):
            return G, labels                         # project every cached row, then select (no n x D copy)
        return G[labels], None

    def autosimilar(self, c):
        """intervalxt's SimilarityTracker: on a SAF-0 component, save the POSITIONAL state
        (permutation, length vector) at Brent checkpoints (ttl 1, 2, 4, ...) and report True when the
        current state is the saved one RESCALED -- exact proportionality in the field.  Zorich
        induction is scale-invariant and depends only on the positional state, so the induction then
        cycles forever without a connection: NO PERIODIC TRAJECTORY (intervalxt's
        WITHOUT_PERIODIC_TRAJECTORY_AUTO_SIMILAR).  Call at the loop head (same phase every step);
        callers reset c.trk / c.saf0 to None on a merge or a reduce.
        ARMED only after `AUTOSIM_AFTER` consecutive steps with no merge or reduce.  A component
        that is converging to cylinders merges every few steps, and each merge reset the cache, so
        the SAF test re-ran after nearly every merge and the exact proportionality test (`F.mul`,
        O(D^2) at degree D) fired whenever a permutation recurred -- ~80% of `n_flow`'s time on CP rows
        (13/105: 7.2 s -> 1.6 s with AUTOSIM off), for a certificate that has fired on NO row on record
        ([NCYL-641]).  A self-similar component never merges again once in its cycle, so delaying the
        arming delays its certificate by `AUTOSIM_AFTER` steps and loses nothing; `AUTOSIM_AFTER = 0`
        is exactly the s635 code path."""
        if not AUTOSIM:
            return False
        if c.trk is None:
            c.trk = {'age': 0, 'ttl': 0, 'bound': 1, 'perm': None, 'lens': None}
        k = c.trk
        k['age'] += 1
        if k['age'] <= AUTOSIM_AFTER:                # not yet armed -- see the docstring
            return False
        if c.saf0 is None:
            c.saf0 = not self.saf_nonzero(c)
        if not c.saf0:
            return False
        F, lam = self.F, self.lam
        pos = {l: i for i, l in enumerate(c.top)}
        perm = [pos[l] for l in c.bot]
        if k['ttl'] <= 0:
            k['perm'], k['lens'] = perm, [lam[l] for l in c.top]
            k['ttl'] = k['bound']
            k['bound'] *= 2
            return False
        k['ttl'] -= 1
        if perm != k['perm']:
            return False
        s0 = k['lens']
        a0 = lam[c.top[0]]
        return all(F.cmp(F.mul(lam[l], s0[0]), F.mul(a0, s0[i])) == 0 for i, l in enumerate(c.top))

    def boshernitzan(self, c):
        """A verified y in Z^D with y . t_j > 0 for all j, or None."""
        from scipy.optimize import linprog
        tr = self.translations(c)
        T = [[int(x) for x in tr[l]] for l in c.top]   # Python ints: y . t must not wrap
        D = self.F.D
        M = np.array([[float(x) for x in t] for t in T])
        sc = np.abs(M).max()
        if sc == 0:
            return None
        M = M / sc
        # max s st M y >= s, -1 <= y <= 1
        cost = np.zeros(D + 1)
        cost[-1] = -1.0
        A = np.hstack([-M, np.ones((len(T), 1))])
        res = linprog(cost, A_ub=A, b_ub=np.zeros(len(T)),
                      bounds=[(-1, 1)] * D + [(None, 1)], method='highs',
                      options={'time_limit': LP_SECONDS})   # a stalled LP is 'no certificate yet'
        if res.status != 0 or res.x[-1] <= 1e-12:
            return None
        y0 = res.x[:D]
        for scale in (10 ** 6, 10 ** 9, 10 ** 12, 10 ** 15):
            y = [int(round(v * scale)) for v in y0]
            if all(sum(a * b for a, b in zip(y, t)) > 0 for t in T):
                return y
        return None

    def run(self, top, bot, cap, ctl=None):
        """`ctl` (a `RunControl`, s761): called once per round with the full state -- progress lines,
        checkpoints, a stop file; if it holds a checkpoint of THIS run, the run resumes from it."""
        F = self.F
        work = [Comp(top, bot)]
        done = []
        if ctl is not None:
            work, done = ctl.start(self, work, done)
        while work:
            c = work.pop()
            while True:
                if ctl is not None:
                    ctl(self, work, done, c)
                rest = self._reduce(c)
                if rest is not None:
                    work.append(rest)
                    c.trk = c.saf0 = None
                if len(c.top) == 1:
                    c.kind = 'cylinder'
                    break
                # saddle connection at the left end: merge (intervalxt NON_SEPARATING)
                t0, b0 = c.top[0], c.bot[0]
                if t0 != b0 and F.cmp(self.lam[t0], self.lam[b0]) == 0:
                    it = c.bot.index(t0)
                    c.bot[it] = b0
                    self.h[b0] += self.h[t0]
                    del c.bot[0]
                    del c.top[0]
                    c.trk = c.saf0 = None
                    continue
                if c.steps >= cap:
                    c.kind = 'undetermined'
                    break
                if self.autosimilar(c):
                    c.kind, c.cert = 'no_periodic', 'autosimilar'
                    break
                # intervalxt skips Boshernitzan when SAF = 0, where it cannot succeed; the gate is the
                # O(nD) random-projection test (s626 -- the dense test cost 70% at D = 208, and an
                # ungated LP on an all-cylinder component can stall HiGHS for many minutes)
                if c.steps % BOSH_EVERY == 0 and c.steps > 0 and self.saf_nonzero(c):
                    y = self.boshernitzan(c)
                    if y is not None:
                        c.kind, c.cert = 'no_periodic', y
                        break
                c.steps += 1
                if self._zorich(c.top, c.bot):
                    continue
                if self._zorich(c.bot, c.top):
                    continue
            done.append(c)
        return done


def decompose(abc, j, cap=5000):
    T = Triangle(abc)
    pieces, lengths, top, bot, meta = _build_iet(T, j)
    F = _fast(T)                                     # numpy-vector elements, tuple-identical outputs
    Dc = Decomposer(F, lengths)
    comps = Dc.run(top, bot, cap)
    out = []
    total = F.zero
    for c in comps:
        meas, mass = F.zero, F.zero
        for l in c.top:
            meas = F.add(meas, Dc.lam[l])
            mass = F.add(mass, F.smul(Dc.h[l], Dc.lam[l]))
        total = F.add(total, mass)
        out.append({'kind': c.kind, 'size': len(c.top), 'steps': c.steps,
                    'induced_len': F.fval(meas), 'mass': F.fval(mass) / F.fval(meta['L']),
                    'period': Dc.h[c.top[0]] if len(c.top) == 1 else None, 'cert': c.cert,
                    'saf0': Dc.saf_zero(c) if c.kind != 'cylinder' else None})
    # Kac: every orbit of a component meets its induced interval, so sum h*lam over ALL labels is
    # the whole transversal -- EXACTLY.  A mis-tracked return time or a lost/duplicated piece fails.
    assert F.cmp(total, meta['L']) == 0, 'Kac mass balance failed'
    kinds = {}
    for o in out:
        kinds[o['kind']] = kinds.get(o['kind'], 0) + 1
    return {'abc': abc, 'j': j, 'N': 4 * T.d, 'npieces': len(pieces), 'meta_ncopies': meta['ncopies'],
            'L': F.fval(meta['L']), 'counts': kinds, 'components': out,
            'nsign_iv': F.nsign_iv + T.F.nsign_iv}


class StopRun(Exception):
    """Raised by `RunControl` when its stop file appears -- AFTER the checkpoint is written."""


class RunControl:
    """Progress, checkpoint and stop for one long `Decomposer.run`.  Pass as `ctl=` to
    `n_flow` / `crossings` / `Decomposer.run`.
      log(line)   every `every` s: rounds, components finished (cylinders / other), components
                  pending, and `found` = the fraction of the transversal's Kac mass already in
                  finished CYLINDERS (float, display only) -- tends to 1 iff the direction is CP;
                  a minimal component holds the rest at the end.
      ckpt        path; the state (lengths, return times, pending/finished components) is pickled
                  there atomically every `ckpt_every` s and on stop.  `start` resumes from it when
                  its `tag` matches (the IET is rebuilt deterministically by the caller).
      stop        path; when it exists, checkpoint and raise `StopRun`.
      status      optional callable(dict), called with each progress record (e.g. a status file)."""

    def __init__(self, tag, ckpt=None, stop=None, every=60, ckpt_every=300, log=print, status=None):
        import time
        self.tag, self.ckpt, self.stop = tag, ckpt, stop
        self.every, self.ckpt_every, self.log, self.status = every, ckpt_every, log, status
        self.L = None
        self.t0 = self.t_log = self.t_ck = time.monotonic()
        self.resumed = False
        self.prior = 0.0                             # seconds spent before the resumed checkpoint
        self._nd, self._mass = 0, 0.0

    def start(self, Dc, work, done):
        import os
        import pickle
        if self.ckpt and os.path.exists(self.ckpt):
            st = pickle.load(open(self.ckpt, 'rb'))
            if st['tag'] == self.tag:
                Dc.lam, Dc.h = st['lam'], st['h']
                if hasattr(Dc.F, 'vec'):                 # tuples on disk -> the run's elements
                    Dc.lam = {l: Dc.F.vec(x) for l, x in Dc.lam.items()}
                self.resumed, self.prior = True, st['secs']
                self.log(f"[{self.tag}] resumed from checkpoint: {len(st['done'])} components done, "
                         f"{len(st['work'])} pending, {st['secs']:.0f} s already spent")
                return st['work'], st['done']
            self.log(f"[{self.tag}] checkpoint tag {st['tag']} does not match -- ignored")
        return work, done

    def save(self, Dc, work, done, c):
        import os
        import pickle
        import time
        lam = {l: tuple(int(a) for a in x) for l, x in Dc.lam.items()}   # tuples, as before
        st = {'tag': self.tag, 'lam': lam, 'h': Dc.h, 'work': work + [c], 'done': done,
              'secs': self.prior + time.monotonic() - self.t0}
        tmp = self.ckpt + '.tmp'
        with open(tmp, 'wb') as fh:
            pickle.dump(st, fh, protocol=pickle.HIGHEST_PROTOCOL)
        os.replace(tmp, self.ckpt)

    def record(self, Dc, work, done, c):
        import time
        F = Dc.F
        for cc in done[self._nd:]:
            if cc.kind == 'cylinder':
                for l in cc.top:
                    h = Dc.h[l]
                    self._mass += float(getattr(h, 'h', h)) * F.fval(Dc.lam[l])
        self._nd = len(done)
        ncyl = sum(1 for cc in done if cc.kind == 'cylinder')
        rounds = sum(cc.steps for cc in done) + sum(cc.steps for cc in work) + c.steps
        return {'tag': self.tag, 'secs': round(self.prior + time.monotonic() - self.t0, 1),
                'rounds': rounds, 'done': len(done), 'cylinders': ncyl, 'other': len(done) - ncyl,
                'pending': len(work) + 1, 'labels_pending': len(c.top) + sum(len(w.top) for w in work),
                'found': (self._mass / F.fval(self.L)) if self.L is not None else None}

    def __call__(self, Dc, work, done, c):
        import os
        import time
        now = time.monotonic()
        stop = bool(self.stop) and os.path.exists(self.stop)
        if now - self.t_log >= self.every or stop:
            self.t_log = now
            r = self.record(Dc, work, done, c)
            fd = f"{r['found']:.6f}" if r['found'] is not None else '?'
            self.log(f"[{self.tag}] {r['secs']:.0f} s  rounds {r['rounds']}  done {r['done']} "
                     f"({r['cylinders']} cyl, {r['other']} other)  pending {r['pending']} "
                     f"({r['labels_pending']} labels)  found {fd}")
            if self.status:
                self.status(r)
        if self.ckpt and (now - self.t_ck >= self.ckpt_every or stop):
            self.t_ck = now
            self.save(Dc, work, done, c)
        if stop:
            raise StopRun(self.tag)


class _HK:
    """(return time, crossings of each tracked entry key) -- composes like the int h."""
    __slots__ = ('h', 'k')

    def __init__(self, h, k):
        self.h, self.k = h, k

    def __add__(self, o):
        return _HK(self.h + o.h, tuple(a + b for a, b in zip(self.k, o.k)))

    def __rmul__(self, q):
        return _HK(q * self.h, tuple(q * a for a in self.k))


def crossings(abc, j, keys, cap=20000, ctl=None):
    """`decompose`, with every component's labels also carrying how many times their return orbit
    crosses each entry key `(g, s)` in `keys` (a copy crossing entered through side s of copy g).
    A cylinder's count at a key is its crossing multiplicity k_c there: the number of intervals it
    cuts on that edge.  Kac is asserted exactly, on the return times, as in `decompose`.
    Field elements ('lam', the edge lengths) are `Cyc` tuples."""
    out, edge, _ = _crossings(abc, j, keys, cap, ctl)
    tu = CycNP.tup
    for o in out:
        if o['lam'] is not None:
            o['lam'] = tu(o['lam'])
    return out, {key: tu(x) for key, x in edge.items()}


def _crossings(abc, j, keys, cap=20000, ctl=None):
    """`crossings` with elements native to its `CycNP` field, returned third."""
    T = Triangle(abc)
    pieces, lengths, top, bot, meta = _build_iet(T, j)
    idx = {key: i for i, key in enumerate(keys)}
    Dc = Decomposer(_fast(T), lengths)
    zero = (0,) * len(keys)
    Dc.h = {}
    for p in lengths:
        k = list(zero)
        if pieces[p][3] in idx:
            k[idx[pieces[p][3]]] = 1
        Dc.h[p] = _HK(1, tuple(k))
    if ctl is not None:
        ctl.L = meta['L']
    comps = Dc.run(top, bot, cap, ctl)
    F = Dc.F
    total = F.zero
    out = []
    for c in comps:
        for l in c.top:
            total = F.add(total, F.smul(Dc.h[l].h, Dc.lam[l]))
        out.append({'kind': c.kind, 'size': len(c.top),
                    'period': Dc.h[c.top[0]].h if len(c.top) == 1 else None,
                    'k': Dc.h[c.top[0]].k if len(c.top) == 1 else None,
                    'lam': Dc.lam[c.top[0]] if len(c.top) == 1 else None})
    assert F.cmp(total, meta['L']) == 0, 'Kac mass balance failed'
    edge = {key: F.zero for key in keys}             # exact length of each tracked entry edge
    for p in lengths:
        if pieces[p][3] in idx:
            edge[pieces[p][3]] = F.add(edge[pieces[p][3]], lengths[p])
    return out, edge, F


def _dps(F, u):
    """Working digits for a float READ-OUT of u: 50, or 40 + the digits of u's largest coefficient if that is more.
    A fixed 50 cancelled catastrophically once coefficients passed ~1e30 (377/610: branch widths read NEGATIVE, -1.5e-23);
    rows whose coefficients fit keep exactly the old 50-digit value."""
    m = max(1, max(abs(int(a)) for a in (u.tolist() if hasattr(u, 'tolist') else u)))
    return max(50, 40 + len(str(m)))


def n_flow(P, Q, cap=20000, beam='1', ctl=None):
    """n(P/Q), the number of branches of the perpendicular partition of L1 (the orphan paper's
    object; `lib/partition.py`), from the flow decomposition: the right triangle
    `(P, Q-P, Q)` (V0 = O, V1 = A, V2 = R; side 1 = AR = L1), direction perpendicular to L1.  Every
    perpendicular beam retraces, so a branch is one crossing of an L1 edge by a cylinder and
    n = sum_c k_c(E) for any L1 edge E perpendicular to the flow -- computed at EVERY such edge and
    required equal.  Returns {'n', 'ncyl', 'ncomp', 'edges', 'tiled', 'all_cyl', 'widths',
    'prefix_len'}: one entry per branch, as a MULTISET (sorted by prefix length, then width -- NOT in
    order along L1: positions and boundary vertex types are not carried).  width = the cylinder's
    transverse width / |E| (L1 has length 1 in the paper's triangle); prefix_len = period / 2 (the
    closed orbit is the prefix and its reversal).  NOT every
    perpendicular class is CP (8/15 = T0's centre holds a minimal component), but every beam from L1 retraces, so no minimal component meets an L1 edge: the
    CERTIFICATE is exact -- sum_c k_c(E) * width_c = |E|, the cylinders tile the edge -- and n is
    returned only when it holds (then non-cylinder components are irrelevant to n).
    beam = 'h': the same for the beam perpendicular to the hypotenuse OA (side 0; `partition(...,
    beam='h')`); widths are then fractions of |OA|."""
    assert math.gcd(P, Q) == 1 and 0 < P < Q and beam in ('1', 'h')
    abc = (P, Q - P, Q)
    T = Triangle(abc)
    side = 1 if beam == '1' else 0
    v1, v2 = T.sides[side]
    j = perp_dirs(abc)[side]
    om = {j % T.N: 1}
    keys = []
    for g in [(k, e) for e in (0, 1) for k in range(T.d)]:
        W = [T.g_apply(g, v) for v in T.V]
        e_ = madd(W[v2], mneg(W[v1]))
        if not any(T.red(T.cross2(mshift(om, T.d, T.N), e_))):   # e parallel to i*omega
            keys.append((g, side))                  # the beam side of copy g, perpendicular to omega
    comps, edge, F = _crossings(abc, j, keys, cap, ctl)
    cyl = [c for c in comps if c['kind'] == 'cylinder']
    sums, tiled = [], []
    for i, key in enumerate(keys):
        if F.sign(edge[key]) == 0:
            continue                                # not an entry key (the flow exits there)
        cover = F.zero
        for c in cyl:
            cover = F.add(cover, F.smul(c['k'][i], c['lam']))
        sums.append(sum(c['k'][i] for c in cyl))
        tiled.append(F.cmp(cover, edge[key]) == 0)
    ok = bool(sums) and all(tiled) and len(set(sums)) == 1
    widths, plen = [], []
    if ok:
        i0 = next(i for i, key in enumerate(keys) if F.sign(edge[key]))
        E = F.mpval(edge[keys[i0]], _dps(F, edge[keys[i0]]))
        for c in cyl:
            wc = float(F.mpval(c['lam'], _dps(F, c['lam'])) / E)       # fval cancels on small elements
            for _ in range(c['k'][i0]):
                widths.append(wc)
                assert c['period'] % 2 == 0, 'a retracing orbit has even period'
                plen.append(c['period'] // 2)
        order = sorted(range(len(widths)), key=lambda t: (plen[t], widths[t]))
        widths, plen = [widths[t] for t in order], [plen[t] for t in order]
    return {'n': sums[0] if ok else None, 'ncyl': len(cyl), 'ncomp': len(comps), 'edges': sums,
            'tiled': tiled, 'all_cyl': len(cyl) == len(comps),
            'kinds': {k: sum(1 for c in comps if c['kind'] == k) for k in sorted({c['kind'] for c in comps})},
            'widths': widths, 'prefix_len': plen}


def m_par_flow(P, Q, cap=20000):
    """m_par(P/Q) of `necklace_ledger.md` §6.18 -- the number of distinct cylinders of H's class met by
    the L1 strip `L1:2Q+P` -- from the exact flow decomposition PARALLEL TO H (side 0), tracing no
    prong.  The unfolding contains the pi-rotation, so every billiard cylinder lifts to TWO
    translation cylinders, one per orientation: m_par = #(translation cylinders with positive
    crossing at either `L1` entry key of label 0) / 2, label := flow exponent in the copy's frame minus
    `4 m_1` (mod N), the copy frame of key `(k, e)` being `j - 4k` (`e = 0`) or `-(j - 4k)` (`e = 1`).
    CERTIFIED: both label-0 keys tile exactly (sum_c k_c lam_c = |E|) and the met count is even;
    `m` is None otherwise.  Gate: 16/16 against `prong_jit` rows (`Q <= 29`), and 11 rows
    where the prong tracer CAPPED (`P = 3, 5, 7, 9`, `Q <= 151`) -- all all-cylinder, all >= the
    tracer's lower bound.  Cost: 1.5-2.9 s at Q = 91, 7-11 s at Q = 151 (P = 6-10) before the
    numpy-vector arithmetic, 3-4x less after it."""
    assert math.gcd(P, Q) == 1 and 0 < P < Q
    abc = (P, Q - P, Q)
    T = Triangle(abc)
    j = side_dirs(abc)[0]
    keys = [((k, e), 1) for e in (0, 1) for k in range(T.d)]
    comps, edge, F = _crossings(abc, j, keys, cap)
    cyl = [c for c in comps if c['kind'] == 'cylinder']
    K0 = []
    for i, ((k, e), _) in enumerate(keys):
        jp = ((j - 4 * k) if e == 0 else -(j - 4 * k)) % T.N
        if (jp - 4 * T.m[1]) % T.N == 0 and F.sign(edge[keys[i]]) != 0:
            K0.append(i)
    tiled = []
    for i in K0:
        cover = F.zero
        for c in cyl:
            cover = F.add(cover, F.smul(c['k'][i], c['lam']))
        tiled.append(F.cmp(cover, edge[keys[i]]) == 0)
    met = sum(1 for c in cyl if any(c['k'][i] > 0 for i in K0))
    ok = bool(K0) and all(tiled) and met % 2 == 0
    return {'m': met // 2 if ok else None, 'met_lifts': met, 'ncyl': len(cyl), 'ncomp': len(comps),
            'all_cyl': len(cyl) == len(comps), 'tiled': tiled}


class PosDecomposer(Decomposer):
    """`Decomposer` that also carries WHERE each label's return orbit crosses the tracked entry edges.

    `self.pos[l]` = list of `(key index, offset)`: a point at distance `u` from the left end of label
    `l`'s domain crosses that key's edge at `offset + u` (edge coordinate from its `lo` end, exact).
    Every induction step is a union with a shift: a label that composes the winner `t` at
    relative position `off` (its place among the stacked labels) gains `t`'s crossings shifted by
    `off + i*S` for each of its `i < q'` passes (`S` the stack's total length), and the winner's
    domain loses its left part of length `R = lam_old - lam_new`, shifting its own offsets by `R`.
    The same formulas hold for the inverse step (`_zorich(bot, top)`) and the connection merge.
    Nothing here certifies itself: `n_flow_pos` gates every run by the exact tiling of each edge."""

    def __init__(self, F, lengths, pos0, maxlist=64):
        super().__init__(F, lengths)
        self.pos = {l: list(pos0.get(l, ())) for l in lengths}
        self.maxlist = maxlist

    def _absorb(self, lab, t, off, S, npass):
        if not self.pos[t]:
            return
        F = self.F
        add = []
        for i in range(npass):
            sh = F.add(off, F.smul(i, S)) if i else off
            add += [(kk, F.add(p, sh)) for kk, p in self.pos[t]]
        self.pos[lab] += add
        if len(self.pos[lab]) > self.maxlist:
            raise RuntimeError('crossing list overflow (a non-cylinder component?)')

    def _shift(self, t, R):
        F = self.F
        self.pos[t] = [(kk, F.add(p, R)) for kk, p in self.pos[t]]

    def _zorich(self, top, bot):
        F, lam = self.F, self.lam
        t = top[0]
        if bot[0] == t:
            return True
        stack, S = [], F.zero
        k = 0
        while True:
            b = bot[k]
            if b == t:
                lt0 = lt = lam[t]
                q = self._floordiv(lt, S)
                lt = F.sub(lt, F.smul(q, S))
                if F.sign(lt) == 0:
                    lt = F.add(lt, S)
                    q -= 1
                nmove = len(stack)
                npart = 0
                for idx, lab in enumerate(stack):
                    if F.cmp(lam[lab], lt) >= 0:
                        nmove = idx if idx > 0 else len(stack)
                        npart = idx
                        break
                    lt = F.sub(lt, lam[lab])
                else:
                    raise RuntimeError('floor division inconsistent')
                lam[t] = lt
                ht = self.h[t]
                off = F.zero
                for idx, lab in enumerate(stack):
                    npass = q + (1 if idx < npart else 0)
                    self.h[lab] += npass * ht
                    self._absorb(lab, t, off, S, npass)
                    off = F.add(off, lam[lab])
                self._shift(t, F.sub(lt0, lt))
                k = nmove
                break
            S2 = F.add(S, lam[b])
            if F.cmp(S2, lam[t]) >= 0:
                lam[t] = F.sub(lam[t], S)
                off = F.zero
                for lab in stack:
                    self.h[lab] += self.h[t]
                    self._absorb(lab, t, off, S, 1)
                    off = F.add(off, lam[lab])
                self._shift(t, S)
                break
            stack.append(b)
            S = S2
            k += 1
        if k:
            it = bot.index(t)
            bot[:] = bot[k:it] + bot[:k] + bot[it:]
        return F.cmp(lam[top[0]], lam[bot[0]]) == 0

    def run(self, top, bot, cap):
        # `Decomposer.run` with the connection merge carrying positions (b0 composes t0 at offset 0)
        F = self.F
        work = [Comp(top, bot)]
        done = []
        while work:
            c = work.pop()
            while True:
                rest = self._reduce(c)
                if rest is not None:
                    work.append(rest)
                if len(c.top) == 1:
                    c.kind = 'cylinder'
                    break
                t0, b0 = c.top[0], c.bot[0]
                if t0 != b0 and F.cmp(self.lam[t0], self.lam[b0]) == 0:
                    it = c.bot.index(t0)
                    c.bot[it] = b0
                    self.h[b0] += self.h[t0]
                    self._absorb(b0, t0, F.zero, F.zero, 1)
                    del c.bot[0]
                    del c.top[0]
                    continue
                if c.steps >= cap:
                    c.kind = 'undetermined'
                    break
                if c.steps % BOSH_EVERY == 0 and c.steps > 0 and self.saf_nonzero(c):
                    y = self.boshernitzan(c)
                    if y is not None:
                        c.kind, c.cert = 'no_periodic', y
                        break
                c.steps += 1
                if self._zorich(c.top, c.bot):
                    continue
                if self._zorich(c.bot, c.top):
                    continue
            done.append(c)
        return done


def n_flow_pos(P, Q, cap=20000, beam='1', exact=False):
    """`n_flow` with POSITIONS: each L1 branch as `(start, width, prefix_len, cyl)` on the edge,
    `start` and `width` as fractions of |E| (evaluated at `_dps` digits: 50, more for huge coefficients — s869), sorted by `start`, measured from
    the edge's `lo` end; `lo_vertex` names that end ('A' or 'R' for beam '1'; 'O' or 'A' for 'h').
    GATE, exact, on every entry edge: the cylinders' crossing intervals, sorted by exact start, are
    contiguous from 0 and end at |E| (a tiling -- a bookkeeping error in the position algebra breaks
    it), and each cylinder's crossing count equals its induction count `k`.  Returns None for n
    unless the gate holds; also returns n_flow's fields for cross-checking.  `exact=True` adds
    `out['exact']` = {'N', 'E', 'branches': [(start, width)]}: the same branches, exact power-basis
    tuples in Z[zeta_N] in transverse units (divide by E for the fraction of the edge).
    COST: ~60 s at 13/385, ~7x faster than the PRE-s767 `n_flow`. ⚠ against today's `n_flow` it is
    ~5-8% SLOWER (9.5-25.8 s vs 8.8-23.9 s at 50/151, 75/151, 100/201) and returns no `all_cyl` -- use `n_flow` for CP.
    s869 (CycNP arithmetic, outputs identical): 55/89 0.7 s (was 2.8), 89/144 1.3-1.6 s (6.4), 144/233 8.4-8.9 s
    (61), 233/377 19.6 s (160) -- timed on a loaded 4-core box."""
    assert math.gcd(P, Q) == 1 and 0 < P < Q and beam in ('1', 'h')
    abc = (P, Q - P, Q)
    T = Triangle(abc)
    side = 1 if beam == '1' else 0
    v1, v2 = T.sides[side]
    j = perp_dirs(abc)[side]
    om = {j % T.N: 1}
    keys = []
    for g in [(k, e) for e in (0, 1) for k in range(T.d)]:
        W = [T.g_apply(g, v) for v in T.V]
        e_ = madd(W[v2], mneg(W[v1]))
        if not any(T.red(T.cross2(mshift(om, T.d, T.N), e_))):
            keys.append((g, side))
    F = _fast(T)                                     # numpy-vector elements; tuples handed out
    pieces, lengths, top, bot, meta = _build_iet(T, j)
    idx = {key: i for i, key in enumerate(keys)}
    # entry-segment origin of each key: the least domain start among its pieces
    origin = {}
    for p in lengths:
        kk = pieces[p][3]
        if kk in idx and (kk not in origin or F.cmp(pieces[p][0], origin[kk]) < 0):
            origin[kk] = pieces[p][0]
    pos0 = {p: [(idx[pieces[p][3]], F.sub(pieces[p][0], origin[pieces[p][3]]))]
            for p in lengths if pieces[p][3] in idx}
    Dc = PosDecomposer(F, lengths, pos0)
    zero = (0,) * len(keys)
    Dc.h = {}
    for p in lengths:
        kv = list(zero)
        if pieces[p][3] in idx:
            kv[idx[pieces[p][3]]] = 1
        Dc.h[p] = _HK(1, tuple(kv))
    comps = Dc.run(top, bot, cap)
    total = F.zero
    for c in comps:
        for l in c.top:
            total = F.add(total, F.smul(Dc.h[l].h, Dc.lam[l]))
    assert F.cmp(total, meta['L']) == 0, 'Kac mass balance failed'
    edge = {key: F.zero for key in keys}
    for p in lengths:
        if pieces[p][3] in idx:
            edge[pieces[p][3]] = F.add(edge[pieces[p][3]], lengths[p])
    cyl = [c.top[0] for c in comps if c.kind == 'cylinder']
    ncomp = len(comps)
    import functools
    ok, per_key = True, []
    for i, key in enumerate(keys):
        if F.sign(edge[key]) == 0:
            continue
        iv = [(p, Dc.lam[l], l) for l in cyl for kk, p in Dc.pos[l] if kk == i]
        cnt = sum(Dc.h[l].k[i] for l in cyl)
        iv.sort(key=functools.cmp_to_key(lambda a, b: F.cmp(a[0], b[0])))
        cur = F.zero
        for p, w, _ in iv:
            if F.cmp(p, cur) != 0:
                ok = False
                break
            cur = F.add(cur, w)
        ok = ok and F.cmp(cur, edge[key]) == 0 and len(iv) == cnt
        per_key.append((i, iv))
    if not ok or not per_key:
        return {'n': None, 'ok': False, 'ncyl': len(cyl), 'ncomp': ncomp}
    i0, iv = per_key[0]
    E = F.mpval(edge[keys[i0]], _dps(F, edge[keys[i0]]))
    g0 = keys[i0][0]
    W = [T.g_apply(g0, v) for v in T.V]
    x1, x2 = F.red(T.cross2(om, W[v1])), F.red(T.cross2(om, W[v2]))
    names = {0: 'O', 1: 'A', 2: 'R'}
    lo_v = names[v1] if F.cmp(x1, x2) < 0 else names[v2]
    br = [(float(F.mpval(p, _dps(F, p)) / E), float(F.mpval(w, _dps(F, w)) / E), Dc.h[l].h // 2, cyl.index(l))
          for p, w, l in iv]
    out = {'n': len(iv), 'ok': True, 'ncyl': len(cyl), 'ncomp': ncomp, 'branches': br,
           'lo_vertex': lo_v, 'nkeys': len(per_key)}
    if exact:       # the same branches as exact power-basis tuples in Z[zeta_N] (N = 8Q), edge units
        out['exact'] = {'N': T.N, 'E': F.tup(edge[keys[i0]]),
                        'branches': [(F.tup(p), F.tup(w)) for p, w, _ in iv]}
    return out


def _selftest():
    """Field against lib/cyclofield.py on random elements."""
    import random
    import sys
    sys.path.insert(0, '/workarea/billiards_deposit')
    from lib.cyclofield import CycloField
    N = 60
    F, G = Cyc(N), CycloField(N)
    rng = random.Random(1)
    for _ in range(200):
        u = tuple(rng.randint(-5, 5) for _ in range(F.D))
        v = tuple(rng.randint(-5, 5) for _ in range(F.D))
        gu = sum((G.zeta(k) * c for k, c in enumerate(u) if c), G.zero)
        gv = sum((G.zeta(k) * c for k, c in enumerate(v) if c), G.zero)
        w = F.mul(u, v)
        gw = sum((G.zeta(k) * c for k, c in enumerate(w) if c), G.zero)
        assert (gw - gu * gv).is_zero()
        r = F.add(u, F.conj(u))
        gr = gu + gu.conj()
        assert F.sign(r) == gr.sign()
        if any(r):                                   # the fixed-point path vs the enclosure
            m = sum(map(abs, r))
            assert F._sign_fx(r, m) == F._sign_iv(r) == gr.sign()
    # near-ties: (2 cos(2 pi/N))^k against its integer-coefficient expansion minus a tiny shift
    G1 = Cyc(1540)
    for k in range(1, 6):
        u = list(G1.mul(G1.add(G1.z(1), G1.z(-1)), G1.one) if k == 1 else u)
        if k > 1:
            u = list(G1.mul(tuple(u), G1.add(G1.z(1), G1.z(-1))))
        for e in (10 ** 12, -10 ** 12):
            w = tuple(a * e for a in u)
            w = G1.sub(w, G1.smul(int(e * float(G1.fval(tuple(u)))), G1.one))
            assert G1._sign_fx(w, sum(map(abs, w))) == G1._sign_iv(w)
    # CycNP against Cyc -- every op, coefficients from 0 to 2^100 (int64 / object promotion and
    # demotion), and near-ties through the vectorised fixed-point stage (against the enclosure)
    for N, reps in ((60, 150), (1540, 50), (1864, 30)):
        F, V = Cyc(N), CycNP(N)
        tu = V.tup
        for it in range(reps):
            sz = rng.choice((3, 30, 61, 62, 63, 64, 100))
            u = tuple(rng.randint(-2 ** sz, 2 ** sz) if rng.random() < 0.7 else 0 for _ in range(F.D))
            v = tuple(rng.randint(-2 ** sz, 2 ** sz) if rng.random() < 0.7 else 0 for _ in range(F.D))
            U, W = V.vec(u), V.vec(v)
            assert tu(U) == u and (U.dtype == np.int64) == (max(map(abs, u)) < 2 ** 62)
            assert tu(V.add(U, W)) == F.add(u, v) and tu(V.sub(U, W)) == F.sub(u, v)
            assert tu(V.neg(U)) == F.neg(u)
            for n in (0, 1, -1, 3, 2 ** 20, -2 ** 40, 2 ** 70):
                assert tu(V.smul(n, U)) == F.smul(n, u)
            assert V.sign(U) == F.sign(u) and V.cmp(U, W) == F.cmp(u, v)
            assert V.fval(U) == F.fval(u) and V.mpval(U) == F.mpval(u) and V.mpval(U, 50) == F.mpval(u, 50)
            if N == 60 and it < 20:
                assert tu(V.mul(U, W)) == F.mul(u, v) and tu(V.conj(U)) == F.conj(u)
            r = F.add(u, F.conj(u))                      # real: a near-tie against a scaled integer
            x = F.sub(F.smul(2 ** 40, r), F.smul(int(2 ** 40 * F.fval(r)), F.one))
            if max(map(abs, x)) < 2 ** 62 and any(x):
                X = V.vec(x)
                d = V._sign_fx_np(X, sum(map(abs, x)))
                assert d in (0, F._sign_iv(x)) and V.sign(X) == F._sign_iv(x)
        D = V.vec([2 ** 61] * F.D)                       # overflow on add -> object, back to int64 on sub
        S2 = V.add(D, D)
        assert S2.dtype == object and V.sub(S2, D).dtype == np.int64 and tu(V.sub(S2, D)) == tu(D)
    print('selftest ok: mul / conj / sign agree with cyclofield on 200 elements; fixed-point sign = '
          'enclosure on those and on 10 near-ties at N = 1540; CycNP = Cyc on 230 elements at N = 60, 1540, 1864')


if __name__ == '__main__':
    _selftest()
