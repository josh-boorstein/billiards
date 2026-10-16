"""Section 5.3, Theorem 5.3: the infinite family of non-perpendicular cylinders in the class
perpendicular to OB, accumulating at y_c = (sqrt 17 - 1)/4.

Statements checked (the words X, Y = X reversed, P1 are the paper's):
  * (1) for every n >= 1 the cubic Q_n is irreducible and a_n/pi is irrational: the resultant
    behind the hand argument is -1, and, as a check, Q_n has no rational root and a/pi is
    irrational (lib/mirror.py) for n <= 60;
  * the blocks X and Y unfold to translations with b_Y = -b_X and
    b_X = 2 (w^2 - 1)(w^2 + 1)(w^2 + 2)/w^4 (frame O = 0, A = 1, B = i tan a, w = e^{ia});
  * (2) ALIGNMENT, symbolically in m: the image of O after X^m H (n = 2m), resp. P1 Y^m
    (n = 2m + 1), lies on the line from O in the direction 2 j a (j = -3, resp. -4) iff
    Q_n(y) = 0 -- the side function, linear in m, is proportional to Q_n;
  * (2) CLEARANCE, for every n: the 274 constraint functions (each crossed edge's endpoints on
    opposite sides of that line, at the first and last block copy -- the sign is linear in the
    block index), as rational functions of t = 2 cos a on the family, have no real root on the
    family's interval and the required sign at an interior point (Sturm sequences over Q);
  * (2)-(3) formally, for every n: at every hit on OB the direction relative to the copy is
    +-w^{2r} with r != 0, so no path meets a side at a right angle;
  * (3) as an independent check, the actual members n = 1..8 are certified exactly by
    lib/mirror.py: the cylinder word (the diagonal there and back, through the two corners) is
    a primitive non-perpendicular periodic cylinder of period 30n + 6 perpendicular to OB;
  * at y_c the block X is itself a non-perpendicular periodic cylinder of period 30, whose
    word is a 13-bounce corner diagonal there and back.
"""
import sys
from fractions import Fraction

import mpmath as mp
import sympy as sp

from common import check, done
import mirror as m

X = 'HABHBAHAHBAHAHAHBAHAHAHABHAHAB'
Y = X[::-1]
P1 = 'HAHAHABHAHABHBAH'

w, t, y = sp.symbols('w t y')
cj = lambda e: e.subs(w, 1 / w)
RHO = {'A': (sp.Integer(1), sp.Integer(0)), 'B': (sp.Integer(-1), sp.Integer(0)), 'H': (w**-2, 1 - w**-2)}
VERT = {'O': sp.Integer(0), 'A': sp.Integer(1), 'B': (w**2 - 1) / (w**2 + 1)}     # i·tan α = (w² − 1)/(w² + 1)
SIDE = {'A': ('O', 'A'), 'B': ('O', 'B'), 'H': ('A', 'B')}


class Iso:
    """z ↦ a·z + b (f = 0) or a·z̄ + b (f = 1)."""

    def __init__(self, a=sp.Integer(1), b=sp.Integer(0), f=0):
        self.a, self.b, self.f = a, b, f

    def __call__(self, z):
        return sp.simplify(self.a * (z if self.f == 0 else cj(z)) + self.b)

    def then(self, ch):                       # self ∘ ρ_ch
        c, d = RHO[ch]
        if self.f == 0:
            return Iso(sp.expand(self.a * c), sp.expand(self.a * d + self.b), 1)
        return Iso(sp.expand(self.a * cj(c)), sp.expand(self.a * cj(d) + self.b), 0)


def edges(word, g0=None):
    """crossed edges of `word` from copy g0: list of (letter, P-image, Q-image, iso before crossing)."""
    g = g0 or Iso()
    out = []
    for ch in word:
        P, Q = SIDE[ch]
        out.append((ch, g(VERT[P]), g(VERT[Q]), g))
        g = g.then(ch)
    return out, g


def im_t(p, q):
    """Im(p̄·q) = sin α · M(t)/D(t); returns the rational function M/D in t (D = |den|² > 0 off the unit-circle zeros)."""
    e = sp.together(sp.expand(cj(p) * q))
    num, den = sp.fraction(e)
    N = sp.expand(num * cj(den)); Dn = sp.expand(den * cj(den))
    anti = sp.expand(N - cj(N))                                    # = (w − 1/w)·Msym
    Msym = sp.cancel(anti / (w - 1 / w))
    return sp.cancel(sym_to_t(Msym) / sym_to_t(Dn))                # Im N = (N − N̄)/(2i) = sin α·Msym(t)


def sym_to_t(e):
    """a symmetric Laurent polynomial in w → polynomial in t = w + 1/w."""
    e = sp.expand(e)
    pw = sp.Poly(sp.expand(e * w**200), w)
    coeffs = {m - 200: c for (m,), c in pw.terms()}
    res = sp.Integer(0)
    T = {0: sp.Integer(2), 1: t}
    def Tk(k):
        if k not in T:
            T[k] = sp.expand(t * Tk(k - 1) - Tk(k - 2))            # w^k + w^{-k}
        return T[k]
    for k, c in coeffs.items():
        if k > 0:
            assert sp.simplify(c - coeffs.get(-k, 0)) == 0, 'not symmetric'
            res += c * Tk(k)
        elif k == 0:
            res += c
    return sp.expand(res)


def build(words):
    X, Y, P1 = words['X'], words['Y'], words['P1']
    ex, gX = edges(X)
    bX = gX.b; assert sp.simplify(gX.a - 1) == 0 and gX.f == 0
    eh, gH = edges('H', Iso(1, 0, 0))
    eP, gP = edges(P1)
    ey, gY = edges(Y)
    bY = gY.b; assert sp.simplify(gY.a - 1) == 0 and gY.f == 0
    aP, bP = gP.a, gP.b
    assert gP.f == 0
    cons = []
    # even, n = 2m: O′ = b_H + m·b_X; block vertices V₀ + k·b_X (k = 0…m−1); final edge = H of copy g_X^m: V + m·b_X
    E0, E1 = gH.b, bX
    for i, (ch, Pp, Qp, g) in enumerate(ex):
        for lab, V0 in (('P', Pp), ('Q', Qp)):
            cons.append(dict(par='even', where=f'X[{i}]{ch}.{lab}', V0=V0, D=bX, E0=E0, E1=E1, kind='block'))
    for lab, V in (('P', VERT['A']), ('Q', VERT['B'])):
        cons.append(dict(par='even', where=f'final-H.{lab}', V0=V, D=bX, E0=E0, E1=E1, kind='final'))
    # odd, n = 2m+1: O′ = b_P + m·a_P·b_Y; prefix vertices fixed; block vertices a_P(V₀ + k b_Y) + b_P
    E0o, E1o = bP, sp.expand(aP * bY)
    for i, (ch, Pp, Qp, g) in enumerate(eP):
        for lab, V0 in (('P', Pp), ('Q', Qp)):
            cons.append(dict(par='odd', where=f'P1[{i}]{ch}.{lab}', V0=V0, D=sp.Integer(0), E0=E0o, E1=E1o, kind='prefix'))
    for i, (ch, Pp, Qp, g) in enumerate(ey):
        for lab, V0 in (('P', Pp), ('Q', Qp)):
            V = sp.expand(aP * V0 + bP)
            cons.append(dict(par='odd', where=f'Y[{i}]{ch}.{lab}', V0=V, D=E1o, E0=E0o, E1=E1o, kind='block'))
    return cons


def m_of_t(par):
    S = 8 * y**3 + 8 * y**2 - 6 * y - 4; T = 4 * y**2 - 2 * y - 1
    n = -T / S
    m = n / 2 if par == 'even' else (n - 1) / 2
    return sp.cancel(m.subs(y, (t**2 - 2) / 2))


def constraint_funcs(c):
    """F at the relevant k as rational functions of t (k = 0 and k = m−1 for blocks; k = m for the final edge)."""
    A = im_t(c['E0'], c['V0']); B = im_t(c['E1'], c['V0'])
    C = im_t(c['E0'], c['D']) if c['D'] != 0 else sp.Integer(0)
    MK = im_t(c['E1'], c['D']) if c['D'] != 0 else sp.Integer(0)
    assert sp.simplify(MK) == 0
    m = m_of_t(c['par'])
    if c['kind'] == 'prefix':
        return {'k=-': sp.cancel(A + m * B)}
    if c['kind'] == 'final':
        return {'k=m': sp.cancel(A + m * B + m * C)}
    return {'k=0': sp.cancel(A + m * B), 'k=m-1': sp.cancel(A + m * B + (m - 1) * C)}



def t_of_n(n, dps=60):
    mp.mp.dps = dps
    yy = mp.findroot(lambda s: 8 * n * s**3 + (8 * n + 4) * s**2 - (6 * n + 2) * s - (4 * n + 1), (mp.sqrt(17) - 1) / 4)
    return mp.sqrt(2 + 2 * yy)


def rat_near(v, up):
    f = Fraction(str(mp.nstr(v, 40, min_fixed=-1, max_fixed=1)))
    return f + Fraction(1, 10**25) if up else f - Fraction(1, 10**25)


TC_POLY = sp.Poly(t**4 - 3 * t**2 - 2, t)


def certify(f, a, b, sign):
    """0 roots of numerator and denominator on [a, b] after removing the t_c factor (positive on (t_c, ∞)); returns dict."""
    num, den = sp.fraction(sp.cancel(f))
    out = dict(removed_tc=0)
    parts = {}
    for lab, e in (('num', num), ('den', den)):
        P = sp.Poly(sp.expand(e), t)
        while P.degree() > 0 and P.rem(TC_POLY).is_zero:
            P = P.quo(TC_POLY); out['removed_tc'] += (1 if lab == 'num' else -1)
        parts[lab] = P
        out[lab + '_roots'] = int(P.count_roots(sp.Rational(a.numerator, a.denominator), sp.Rational(b.numerator, b.denominator)))
    mid = (a + b) / 2
    val = sp.Rational(mid.numerator, mid.denominator)
    sv = sp.sign(parts['num'].eval(val)) * sp.sign(parts['den'].eval(val))
    out['sign_mid'] = int(sv)
    out['ok'] = out['num_roots'] == 0 and out['den_roots'] == 0 and int(sv) == sign
    return out




def Qn(n):
    return f'{n}*(8*y**3 + 8*y**2 - 6*y - 4) + (4*y**2 - 2*y - 1)'


def corner_cylinder_word(K, W):
    """The closed word of the cylinder around the corner diagonal W (Lemma 6.3's reversal):
    W, two bounces at one corner, W reversed, two at the other -- the corner letters are the
    combination that certifies (at most one does)."""
    for c1 in ('AB', 'BA'):
        for c2 in ('AB', 'BA'):
            word = W + c1 + W[::-1] + c2
            if m.cyclic_ok(word):
                c = m.certify(K, word)
                if c['ok']:
                    return word, c
    return None, None


def main():
    # (1)
    r = sp.symbols('r')
    R = sp.resultant(r**3 + 2 * r**2 - 3 * r - 4, r**2 - r - 1, r)
    check(f'(1) resultant of (r+1)(r^2+r-4) and r^2-r-1 is {R}', R == -1)
    ok = all(sp.Poly(sp.sympify(Qn(n)), y).is_irreducible and m.angle_is_irrational(Qn(n)) for n in range(1, 61))
    check('(1) Q_n irreducible and a_n/pi irrational, n = 1..60', ok)

    # blocks
    ex, gX = edges(X)
    ey, gY = edges(Y)
    bX = sp.expand(gX.b)
    okb = (sp.simplify(gX.a - 1) == 0 and gX.f == 0 and sp.simplify(gY.a - 1) == 0 and gY.f == 0
           and sp.simplify(gY.b + gX.b) == 0
           and sp.simplify(bX - 2 * (w**2 - 1) * (w**2 + 1) * (w**2 + 2) / w**4) == 0)
    check('X, Y unfold to translations, b_Y = -b_X = -2(w^2-1)(w^2+1)(w^2+2)/w^4', okb)

    # (2) alignment, symbolic in m
    S = 8 * y**3 + 8 * y**2 - 6 * y - 4
    Tq = 4 * y**2 - 2 * y - 1
    eh, gH = edges('H', Iso(1, 0, 0))
    eP, gP = edges(P1)
    for par, j, E0, E1, nfun in (('even', -3, gH.b, bX, lambda mm: 2 * mm),
                                 ('odd', -4, gP.b, sp.expand(gP.a * gY.b), lambda mm: 2 * mm + 1)):
        d = w ** (2 * j)
        A0 = im_t(d, E0)                       # Im(conj(delta) O') / sin a = A0 + m A1
        A1 = im_t(d, E1)
        mm = sp.Symbol('m')
        Qm = (nfun(mm) * S + Tq).subs(y, (t**2 - 2) / 2)
        # A0 + m A1 = 0 iff Q_n = 0, identically in m: A0 * dQ/dm - A1 * Q(m=0) == 0
        ident = sp.simplify(A0 * sp.diff(Qm, mm) - A1 * Qm.subs(mm, 0)) == 0 and sp.simplify(A1) != 0
        check(f'(2) alignment, {par} n: the side function of O\' is proportional to Q_n, for all m', ident)
        dw = w ** (2 * (j + 1))                # control: the neighbouring direction index must fail
        B0, B1 = im_t(dw, E0), im_t(dw, E1)
        check(f'    control: with j + 1 in place of j the identity fails',
              sp.simplify(B0 * sp.diff(Qm, mm) - B1 * Qm.subs(mm, 0)) != 0)

    # (2) clearance: the 274 Sturm certificates
    words = dict(X=X, Y=Y, P1=P1)
    cons = build(words)
    tc = mp.sqrt((3 + mp.sqrt(17)) / 2)
    a = rat_near(tc, up=False)
    tops = {('even', 'block'): 2, ('even', 'final'): 2, ('odd', 'prefix'): 1, ('odd', 'block'): 3}
    ref = {'even': t_of_n(4), 'odd': t_of_n(3)}
    nf, fails = 0, 0
    for c in cons:
        fs = constraint_funcs(c)
        n0 = tops[(c['par'], c['kind'])]
        b = rat_near(t_of_n(n0), up=True)
        for kk, f in fs.items():
            nf += 1
            sref = int(mp.sign(sp.lambdify(t, f, 'mpmath')(ref[c['par']])))
            if sref == 0 or not certify(f, a, b, sref)['ok']:
                fails += 1
    check(f'(2) clearance: {nf} constraint functions, none with a root on the interval, signs as required '
          f'({fails} failures)', nf == 274 and fails == 0)

    # (2)-(3) formal non-perpendicularity at OB hits
    viol = 0
    for wd, j in ((X + 'H', -3), (P1 + Y, -4)):
        g = Iso()
        v = w ** (2 * j)
        for ch in wd:
            u = sp.simplify(v / g.a) if g.f == 0 else sp.simplify(cj(v / g.a))
            if ch == 'B':
                pw = sp.Poly(sp.expand(u * w**100), w)
                if len(pw.terms()) != 1 or pw.terms()[0][0][0] == 100:
                    viol += 1
            g = g.then(ch)
    check('(2)-(3) at every OB hit the relative direction is +-w^{2r}, r != 0', viol == 0)

    # (3) the actual members, exactly
    for n in range(1, 9):
        mm = n // 2
        W = X * mm + 'H' if n % 2 == 0 else P1 + Y * mm
        yn = 2 * mp.mpf(1)
        with mp.workdps(40):
            yn = mp.findroot(lambda s: 8 * n * s**3 + (8 * n + 4) * s**2 - (6 * n + 2) * s - (4 * n + 1),
                             (mp.sqrt(17) - 1) / 4 + mp.mpf(1) / (10 * n))
        K = m.t_field(Qn(n), float(yn))
        word, cc = corner_cylinder_word(K, W)
        check(f'(3) n = {n}: diagonal of {len(W)} = 15n+1 bounces; cylinder of period {cc["period"] if cc else "-"} '
              f'= 30n+6, perp OB, non-perpendicular, primitive',
              cc is not None and len(W) == 15 * n + 1 and cc['period'] == 30 * n + 6 and cc['cls'] == 'perp OB'
              and not cc['perpendicular'] and cc['primitive'])

    # the limit point
    Kc = m.t_field('2*y**2 + y - 2', (17 ** 0.5 - 1) / 4)
    cx = m.certify(Kc, X)
    rot = [X[k:] + X[:k] for k in range(len(X))]
    core = any(r[15:28] == r[:13][::-1] and r[13:15] in ('AB', 'BA') and r[28:30] in ('AB', 'BA') for r in rot)
    check('y_c: X is a primitive non-perpendicular periodic cylinder of period 30, perp OB',
          cx['ok'] and cx['period'] == 30 and not cx['perpendicular'] and cx['primitive'] and cx['cls'] == 'perp OB')
    check('y_c: its word is a 13-bounce corner diagonal there and back (W, corner, W reversed, corner)', core)
    check('y_c: and Q_n has no root there (the family accumulates at y_c, not on it)',
          all(sp.sympify(Qn(n)).subs(y, (sp.sqrt(17) - 1) / 4).simplify() != 0 for n in (1, 2, 3)))
    done()


if __name__ == '__main__':
    main()
