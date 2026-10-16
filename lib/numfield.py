"""Exact arithmetic in a real number field Q(theta), with certified signs, and in its CM
extension K(w), w^2 = t*w - 1 -- for unfolding billiard words in a right triangle whose angle
is an ALGEBRAIC irrational multiple of pi.

  RealField(poly, approx)  poly = integer coefficients, LOW to HIGH, irreducible over Q
                           (factored with sympy and asserted); theta := the real root nearest
                           `approx`, kept as an isolating rational interval refined on demand.
                           The constructor raises if that root is not within 1e-6 of `approx`.
  Elt                      + - * / ==, exact: an element is its coefficient vector in the
                           power basis, unique because `poly` is irreducible.  .sign() is
                           CERTIFIED -- a rigorous interval evaluation (mpmath.iv) on a
                           refining enclosure of theta; it terminates for every nonzero element
                           and never returns a wrong sign.  Comparisons go through .sign().
  CM(K, t)                 K(w), w^2 = t*w - 1, for t in K with |t| < 2: w = e^{i alpha},
                           t = 2 cos alpha, alpha in (0, pi).
    CElt a + b*w           + - * .conj() .inv() .norm() (= |z|^2, in K) .re() (= a + b*t/2)
                           .im_over_sin() (= b, since Im(a + b*w) = b*sin alpha)
                           .im_sign() (exact) .complex(dps) (for display only)

No decision in this module is made on a floating-point value.  stdlib + mpmath, and sympy at
construction only.
"""
from fractions import Fraction

import mpmath
import sympy as sp

_iv = mpmath.iv


def _trim(c):
    c = list(c)
    while c and c[-1] == 0:
        c.pop()
    return c


def _pmul(a, b):
    if not a or not b:
        return []
    r = [Fraction(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                r[i + j] += x * y
    return r


def _pdivmod(a, b):
    a = [Fraction(x) for x in _trim(a)]
    b = [Fraction(x) for x in _trim(b)]
    q = [Fraction(0)] * max(len(a) - len(b) + 1, 1)
    while len(a) >= len(b) and a:
        c = a[-1] / b[-1]
        k = len(a) - len(b)
        q[k] = c
        for i, y in enumerate(b):
            a[i + k] -= c * y
        a = _trim(a)
    return _trim(q), a


class RealField:
    def __init__(self, poly, approx, name='θ'):
        poly = [int(c) for c in poly]
        X = sp.Symbol('X')
        P = sp.Poly(list(reversed(poly)), X)
        fl = sp.factor_list(P)[1]
        if len(fl) != 1 or fl[0][1] != 1:
            raise ValueError(f'poly {poly} is not irreducible over Q: {fl}')
        self.poly = poly
        self.deg = len(poly) - 1
        self.name = name
        self._P = P
        # sympy's isolating intervals can be WIDE, so choosing by "nearest midpoint" can pick a conjugate root. Refine every
        # interval first, choose by distance to the INTERVAL, then assert the root is near approx.
        a = Fraction(str(mpmath.mpf(approx)))
        ivs = P.intervals(eps=Fraction(1, 10**12))
        dist = lambda lo, hi: max(lo - a, a - hi, Fraction(0))
        best = min(ivs, key=lambda iv: dist(Fraction(iv[0][0].p, iv[0][0].q), Fraction(iv[0][1].p, iv[0][1].q)))
        self._lo = Fraction(best[0][0].p, best[0][0].q)
        self._hi = Fraction(best[0][1].p, best[0][1].q)
        if abs((self._lo + self._hi) / 2 - a) > Fraction(1, 10**6) * max(1, abs(a)):
            raise ValueError(f'no root of {poly} within 1e-6 of approx {float(a)} (nearest {float((self._lo + self._hi) / 2)})')
        self._mod = [Fraction(c, poly[-1]) for c in poly]          # monic modulus
        self.zero = Elt(self, [])
        self.one = Elt(self, [1])
        self.gen = Elt(self, [0, 1]) if self.deg > 1 else Elt(self, [-Fraction(poly[0], poly[1])])
        self._enc = {}

    # --- root enclosure -------------------------------------------------------------------------------------------------
    def _pval(self, r):
        v = Fraction(0)
        for c in reversed(self.poly):
            v = v * r + c
        return v

    def refine(self, width):
        """bisect the isolating interval until hi − lo ≤ width (exact rational signs of poly)."""
        if self._lo == self._hi:
            return
        slo = self._pval(self._lo) > 0
        while self._hi - self._lo > width:
            mid = (self._lo + self._hi) / 2
            v = self._pval(mid)
            if v == 0:
                self._lo = self._hi = mid
                return
            if (v > 0) == slo:
                self._lo = mid
            else:
                self._hi = mid

    def _fast_refine(self, bits):
        """jump to width 2^-bits using an mpmath Newton root, verified by an exact sign change inside the current interval."""
        if self._lo == self._hi or self._hi - self._lo <= Fraction(1, 2**bits):
            return
        with mpmath.workprec(bits + 40):
            f = lambda X: sum(c * X**i for i, c in enumerate(self.poly))
            mid = mpmath.mpf(self._lo.numerator) / self._lo.denominator / 2 + mpmath.mpf(self._hi.numerator) / self._hi.denominator / 2
            try:
                r = mpmath.findroot(f, mid)
                rf = Fraction(mpmath.nstr(r, int(bits * 0.302) + 20, strip_zeros=False))
            except Exception:
                self.refine(Fraction(1, 2**bits)); return
        h = Fraction(1, 2**(bits + 1))
        lo, hi = max(rf - h, self._lo), min(rf + h, self._hi)
        if lo < hi and (self._pval(lo) > 0) != (self._pval(hi) > 0) and self._pval(lo) != 0:
            self._lo, self._hi = lo, hi
        else:
            self.refine(Fraction(1, 2**bits))

    def enclosure(self, bits):
        if bits not in self._enc:
            self._fast_refine(bits)
            with mpmath.workprec(bits + 20):
                _iv.prec = bits + 20
                self._enc[bits] = (_iv.mpf(self._lo.numerator) / self._lo.denominator,
                                   _iv.mpf(self._hi.numerator) / self._hi.denominator)
        return self._enc[bits]

    def approx(self, dps=30):
        self._fast_refine(int(dps * 3.33) + 10)
        with mpmath.workdps(dps + 10):
            return (mpmath.mpf(self._lo.numerator) / self._lo.denominator + mpmath.mpf(self._hi.numerator) / self._hi.denominator) / 2

    # --- constructors ---------------------------------------------------------------------------------------------------
    def of(self, c):
        if isinstance(c, Elt):
            assert c.K is self
            return c
        if isinstance(c, (int, Fraction)):
            return Elt(self, [c])
        return Elt(self, list(c))

    def from_poly(self, expr, sym):
        """a sympy polynomial (rational coefficients) in `sym`, evaluated at θ."""
        P = sp.Poly(sp.expand(expr), sym)
        cs = [Fraction(int(sp.Rational(c).p), int(sp.Rational(c).q)) for c in reversed(P.all_coeffs())]
        return Elt(self, cs)

    def __repr__(self):
        return f'RealField({self.poly}, θ≈{mpmath.nstr(self.approx(15), 15)})'


class Elt:
    __slots__ = ('K', 'c')

    def __init__(self, K, c):
        self.K = K
        c = [Fraction(x) for x in c]
        if len(c) > K.deg:
            c = _pdivmod(c, K._mod)[1]
        self.c = tuple(_trim(c))

    def _co(self, o):
        if isinstance(o, Elt):
            assert o.K is self.K, 'elements of different fields'
            return o
        return self.K.of(o)

    def __add__(self, o):
        o = self._co(o)
        n = max(len(self.c), len(o.c))
        return Elt(self.K, [(self.c[i] if i < len(self.c) else 0) + (o.c[i] if i < len(o.c) else 0) for i in range(n)])

    __radd__ = __add__

    def __neg__(self):
        return Elt(self.K, [-x for x in self.c])

    def __sub__(self, o):
        return self + (-self._co(o))

    def __rsub__(self, o):
        return self._co(o) - self

    def __mul__(self, o):
        if isinstance(o, (int, Fraction)):
            return Elt(self.K, [x * o for x in self.c])
        o = self._co(o)
        return Elt(self.K, _pmul(list(self.c), list(o.c)))

    __rmul__ = __mul__

    def inv(self):
        if not self.c:
            raise ZeroDivisionError
        # extended Euclid: s·a + u·m = 1
        r0, r1 = [Fraction(x) for x in self.K._mod], list(self.c)
        s0, s1 = [], [Fraction(1)]
        while len(_trim(r1)) > 1:
            q, r = _pdivmod(r0, r1)
            r0, r1 = r1, r
            qs = _pmul(q, s1)
            n = max(len(s0), len(qs))
            s0, s1 = s1, _trim([(s0[i] if i < len(s0) else 0) - (qs[i] if i < len(qs) else 0) for i in range(n)])
        c = _trim(r1)[0]
        return Elt(self.K, [x / c for x in s1])

    def __truediv__(self, o):
        if isinstance(o, (int, Fraction)):
            return Elt(self.K, [x / o for x in self.c])
        return self * self._co(o).inv()

    def __rtruediv__(self, o):
        return self._co(o) * self.inv()

    def __pow__(self, k):
        r, b = self.K.one, self
        if k < 0:
            b, k = self.inv(), -k
        while k:
            if k & 1:
                r = r * b
            b = b * b
            k >>= 1
        return r

    def __eq__(self, o):
        if isinstance(o, (int, Fraction, Elt)):
            return (self - o).c == ()
        return NotImplemented

    def __hash__(self):
        return hash(self.c)

    def iszero(self):
        return self.c == ()

    def sign(self):
        if not self.c:
            return 0
        if len(self.c) == 1:
            return 1 if self.c[0] > 0 else -1
        bits = 64 + 4 * max(max(abs(x.numerator).bit_length(), x.denominator.bit_length()) for x in self.c)
        while True:
            lo, hi = self.K.enclosure(bits)
            with mpmath.workprec(bits + 20):
                _iv.prec = bits + 20
                th = _iv.mpf([lo.a, hi.b])
                v = _iv.mpf(0)
                for x in reversed(self.c):
                    v = v * th + _iv.mpf(x.numerator) / x.denominator
            if v.a > 0:
                return 1
            if v.b < 0:
                return -1
            bits *= 2
            if bits > 1 << 20:
                raise RuntimeError('sign: precision runaway (element should be nonzero)')

    def __lt__(self, o):
        return (self - o).sign() < 0

    def __gt__(self, o):
        return (self - o).sign() > 0

    def __le__(self, o):
        return (self - o).sign() <= 0

    def __ge__(self, o):
        return (self - o).sign() >= 0

    def __abs__(self):
        return -self if self.sign() < 0 else self

    def mp(self, dps=30):
        th = self.K.approx(dps + 10)
        with mpmath.workdps(dps + 10):
            v = mpmath.mpf(0)
            for x in reversed(self.c):
                v = v * th + mpmath.mpf(x.numerator) / x.denominator
            return +v

    def __float__(self):
        return float(self.mp(20))

    def height(self):
        return max((max(abs(x.numerator), x.denominator) for x in self.c), default=0).bit_length()

    def __repr__(self):
        return f'Elt({mpmath.nstr(self.mp(15), 15)}; h={self.height()}b)'


class CM:
    """K(w), w² = t·w − 1, t = 2cos α ∈ K real with |t| < 2 ⇒ w = e^{iα}, α ∈ (0, π)."""

    def __init__(self, K, t):
        self.K = K
        self.t = K.of(t)
        assert self.t > -2 and self.t < 2
        self.w = CElt(self, K.zero, K.one)
        self.one = CElt(self, K.one, K.zero)
        self.zero = CElt(self, K.zero, K.zero)

    def of(self, a, b=0):
        return CElt(self, self.K.of(a), self.K.of(b))


class CElt:
    __slots__ = ('C', 'a', 'b')

    def __init__(self, C, a, b):
        self.C, self.a, self.b = C, a, b

    def _co(self, o):
        if isinstance(o, CElt):
            return o
        return CElt(self.C, self.C.K.of(o), self.C.K.zero)

    def __add__(self, o):
        o = self._co(o)
        return CElt(self.C, self.a + o.a, self.b + o.b)

    __radd__ = __add__

    def __neg__(self):
        return CElt(self.C, -self.a, -self.b)

    def __sub__(self, o):
        return self + (-self._co(o))

    def __rsub__(self, o):
        return self._co(o) - self

    def __mul__(self, o):
        if isinstance(o, (int, Fraction, Elt)):
            return CElt(self.C, self.a * o, self.b * o)
        o = self._co(o)
        # (a + bw)(c + dw) = ac + (ad + bc)w + bd·w², w² = tw − 1
        bd = self.b * o.b
        return CElt(self.C, self.a * o.a - bd, self.a * o.b + self.b * o.a + bd * self.C.t)

    __rmul__ = __mul__

    def conj(self):
        # w̄ = 1/w = t − w
        return CElt(self.C, self.a + self.b * self.C.t, -self.b)

    def norm(self):
        return self.a * self.a + self.a * self.b * self.C.t + self.b * self.b

    def inv(self):
        return self.conj() * self.norm().inv()

    def __truediv__(self, o):
        if isinstance(o, (int, Fraction, Elt)):
            return CElt(self.C, self.a / o, self.b / o)
        return self * self._co(o).inv()

    def __pow__(self, k):
        r, b = self.C.one, self
        if k < 0:
            b, k = self.inv(), -k
        while k:
            if k & 1:
                r = r * b
            b = b * b
            k >>= 1
        return r

    def re(self):
        return self.a + self.b * self.C.t / 2

    def im_over_sin(self):
        """Im(z)/sin α — an element of K with the sign of Im z."""
        return self.b

    def im_sign(self):
        return self.b.sign()

    def iszero(self):
        return self.a.iszero() and self.b.iszero()

    def __eq__(self, o):
        if isinstance(o, (int, Fraction, Elt, CElt)):
            return (self - o).iszero()
        return NotImplemented

    def __hash__(self):
        return hash((self.a.c, self.b.c))

    def complex(self, dps=30):
        tt = self.C.t.mp(dps)
        with mpmath.workdps(dps + 10):
            al = mpmath.acos(tt / 2)
            w = mpmath.expj(al)
            return self.a.mp(dps) + self.b.mp(dps) * w

    def __repr__(self):
        z = self.complex(12)
        return f'CElt({mpmath.nstr(z, 12)})'

