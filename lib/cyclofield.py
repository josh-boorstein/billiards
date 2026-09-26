"""Exact arithmetic in the cyclotomic field Q(zeta_N).

Elements are polynomials in z over Q reduced modulo the cyclotomic polynomial Phi_N(z);
z stands for zeta = exp(2 pi i / N).  Equality is decided on the reduced polynomial, so
it is exact.  The sign of a REAL element is read off a certified interval enclosure
(mpmath.iv), whose precision is raised until the enclosure excludes zero; a nonzero real
element always has a sign, so this terminates, and zero is decided before it is tried.
No decision anywhere in this module is made on a floating-point value.

For a right triangle with centre P/Q the natural choice is N = 4Q: then
zeta = exp(i pi / 2Q), and cos(m pi / 2Q), sin(m pi / 2Q) and i = zeta^Q all lie in the
field, which is what the unfolding in `triangle.py` needs.
"""
from fractions import Fraction

import sympy as sp
from mpmath import iv, mp, mpf

_z = sp.symbols('z')


class CycloField:
    def __init__(self, N):
        self.N = N
        self.phi = sp.Poly(sp.cyclotomic_poly(N, _z), _z, domain='QQ')
        self.degree = self.phi.degree()
        self.zero = self._wrap(sp.Poly(0, _z, domain='QQ'))
        self.one = self(1)
        self.i = self.zeta(N // 4) if N % 4 == 0 else None

    # -- construction ------------------------------------------------------------------
    def _wrap(self, poly):
        return Elem(self, poly.rem(self.phi))

    def __call__(self, q):
        """The rational number q (int, Fraction or sympy Rational) as a field element."""
        q = sp.Rational(q) if not isinstance(q, Fraction) else sp.Rational(q.numerator,
                                                                            q.denominator)
        return self._wrap(sp.Poly(q, _z, domain='QQ'))

    def zeta(self, k):
        return self._wrap(sp.Poly(_z ** (k % self.N), _z, domain='QQ'))

    def cos(self, k):
        """cos(2 pi k / N)."""
        return (self.zeta(k) + self.zeta(-k)) * sp.Rational(1, 2)

    def sin(self, k):
        """sin(2 pi k / N); needs 4 | N."""
        return (self.zeta(k) - self.zeta(-k)) * self.i * sp.Rational(-1, 2)


class Elem:
    __slots__ = ('F', 'p')

    def __init__(self, F, poly):
        self.F, self.p = F, poly

    def _coerce(self, o):
        if isinstance(o, Elem):
            return o
        return self.F(o)

    def __add__(self, o):
        return Elem(self.F, self.p.add(self._coerce(o).p))

    __radd__ = __add__

    def __sub__(self, o):
        return Elem(self.F, self.p.sub(self._coerce(o).p))

    def __rsub__(self, o):
        return self._coerce(o) - self

    def __neg__(self):
        return Elem(self.F, self.p.neg())

    def __mul__(self, o):
        if isinstance(o, (int, Fraction, sp.Rational)):
            return Elem(self.F, self.p.mul_ground(sp.Rational(o) if not isinstance(o, Fraction)
                                                  else sp.Rational(o.numerator, o.denominator)))
        return Elem(self.F, self.p.mul(o.p).rem(self.F.phi))

    __rmul__ = __mul__

    def inverse(self):
        if self.is_zero():
            raise ZeroDivisionError('inverse of 0 in Q(zeta_N)')
        inv = sp.invert(self.p.as_expr(), self.F.phi.as_expr(), _z)
        return Elem(self.F, sp.Poly(inv, _z, domain='QQ').rem(self.F.phi))

    def __truediv__(self, o):
        if isinstance(o, (int, Fraction, sp.Rational)):
            return self * (1 / sp.Rational(o) if not isinstance(o, Fraction)
                           else sp.Rational(o.denominator, o.numerator))
        return self * o.inverse()

    def __rtruediv__(self, o):
        return self._coerce(o) * self.inverse()

    # -- exact predicates --------------------------------------------------------------
    def is_zero(self):
        return self.p.is_zero

    def __eq__(self, o):
        return (self - o).is_zero()

    def __hash__(self):
        return hash(tuple(self.p.all_coeffs()))

    def coeffs(self):
        """Ascending coefficients c_k with self = sum c_k zeta^k (k < degree)."""
        return [sp.Rational(c) for c in reversed(self.p.all_coeffs())]

    def conj(self):
        acc = self.F.zero
        for k, c in enumerate(self.coeffs()):
            if c:
                acc = acc + self.F.zeta(-k) * c
        return acc

    def is_real(self):
        return (self - self.conj()).is_zero()

    # -- certified real value ----------------------------------------------------------
    def enclosure(self, dps=60):
        """A certified interval containing this (real) element."""
        assert self.is_real(), 'enclosure() of a non-real element'
        N = self.F.N
        old = iv.dps
        iv.dps = dps
        try:
            tot = iv.mpf(0)
            for k, c in enumerate(self.coeffs()):
                if c:
                    tot += iv.mpf(int(c.p)) / iv.mpf(int(c.q)) * iv.cos(2 * iv.pi * k / N)
        finally:
            iv.dps = old
        return tot

    def sign(self):
        """-1, 0 or +1, exactly: zero by the polynomial, a nonzero sign by an enclosure."""
        if self.is_zero():
            return 0
        dps = 60
        while True:
            e = self.enclosure(dps)
            if e.a > 0:
                return 1
            if e.b < 0:
                return -1
            dps *= 2
            if dps > 4000:
                raise RuntimeError('enclosure did not separate a nonzero element from 0')

    def __lt__(self, o):
        return (self._coerce(o) - self).sign() > 0

    def __gt__(self, o):
        return (self - self._coerce(o)).sign() > 0

    def __le__(self, o):
        return not self > o

    def __ge__(self, o):
        return not self < o

    def value(self, dps=40):
        """Point value for display only (never used to decide anything)."""
        with mp.workdps(dps + 10):
            N = self.F.N
            tot = mp.mpc(0)
            for k, c in enumerate(self.coeffs()):
                if c:
                    tot += mpf(int(c.p)) / int(c.q) * mp.expjpi(mpf(2 * k) / N)
            return +tot.real if self.is_real() else +tot

    def __repr__(self):
        v = self.value(25)
        return f'<Q(zeta_{self.F.N}) {mp.nstr(v, 20)}>'
