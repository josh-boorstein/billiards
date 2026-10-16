"""Section 6.3, Proposition 6.6 and Table 2: five corner words whose tiles cover (0, pi/2).

Frame: O = 0, A = 1, B = i t, t = tan a.  A word W unfolds by g_i = rho_{w1} o ... o rho_{wi};
the i-th crossed edge is g_{i-1}(side w_i) and O_end = g_k(O).  By Lemma 6.2 (W starts and ends
with H, no doubled letter) a is in the tile R_W iff for every crossed edge the endpoints P, Q
satisfy f_P f_Q < 0, f_P := Im(conj(O_end) P).  Every f_P is (polynomial in t)/(1 + t^2)^j; its
sign-carrying factor is a polynomial in s = t^2.  Everything here is exact (sympy): the tiles
are unions of intervals between real roots, decided at interior rational points.

Statements checked:
  * Table 2: for each of the three words, the sign-carrying factor pair of every crossed edge
    is the printed one (up to a positive constant);
  * the tiles: R_HAHABH = (0, a1), R_HAHBHAHABH = (a2, a3), R_HAHBH = (a4, a5), each a single
    interval, with the printed end polynomials (5s = 3; the smaller positive roots of
    s^3 - 5s^2 - 93s + 41 and s^3 - 77s^2 + 51s + 1; s^2 + 10s - 7 and 7s^2 - 10s - 1) and
    values a1 = 0.65905.., a2 = 0.58131.., a3 = 0.69222.., a4 = 0.68108.., a5 = 0.88970..;
  * the leg swaps HBHBAH, HBHAHBHBAH have tiles pi/2 - R_W, computed directly;
  * the cover: a2 < a1, a4 < a3, pi/2 - a3 < a5, so the five tiles cover (0, pi/2);
  * none of the five words is a palindrome of odd length (Lemma 6.4);
  * all five are non-degenerate (Section 6.2): O_end = F_W(u), u = e^{-2ia}, with F_W the printed
    polynomials 1 + u - 2u^2, 1 - 2u^2 - u^3 + 2u^4, 1 - 2u^2 + u^3, 1 - 2u + 3u^3 - 2u^4 and
    1 - 3u + 2u^2, and no F_W(u) is c u^m F_W(1/u) (so R_W is not a monomial); control: a
    self-reciprocal polynomial is flagged;
  * Lemma 6.3/6.4 at one algebraic angle inside each tile: the corner path's cylinder, of
    period 2k + 4, is certified exactly by lib/mirror.py and is not perpendicular.
"""
import math

import sympy as sp

from common import check, done
import mirror as m

t = sp.symbols('t', positive=True)
s = sp.Symbol('s')
I = sp.I
SIDE = {'A': (sp.Integer(0), sp.Integer(1)), 'B': (sp.Integer(0), I * t), 'H': (sp.Integer(1), I * t)}
WORDS = ['HAHABH', 'HAHBHAHABH', 'HAHBH']


def _refl(p, q):
    d = q - p
    return lambda z: sp.expand(p + d / sp.conjugate(d) * sp.conjugate(z - p))


def unfold(W):
    """Edges (vertex-image pairs) and O_end of word W, exact in t."""
    g = lambda z: z
    edges = []
    for w in W:
        p, q = SIDE[w]
        edges.append((sp.simplify(g(p)), sp.simplify(g(q))))
        g = (lambda f, r: (lambda z: f(r(z))))(g, _refl(p, q))
    return edges, sp.simplify(g(sp.Integer(0)))


def numer(V, Oe):
    f = sp.together(sp.simplify(sp.im(sp.expand(sp.conjugate(Oe) * V))))
    return sp.fraction(sp.factor(f))[0]


def sign_core(num):
    """Sign-carrying factor of a numerator, as a polynomial in s = t² where possible."""
    c, facs = sp.factor_list(num)
    out = sp.sign(c)
    for f, e in facs:
        if f == t:
            continue
        g = sp.expand(f)
        if all(m % 2 == 0 for (m,) in sp.Poly(g, t).monoms()):
            g = sp.expand(g.subs(t, sp.sqrt(s)))
        out *= g ** e
    return sp.factor(sp.expand(out))


def tile(W):
    edges, Oe = unfold(W)
    pairs = [(numer(P, Oe), numer(Q, Oe)) for P, Q in edges]
    polys = {sp.expand(f) for pr in pairs for x in pr for f, _ in sp.factor_list(x)[1] if sp.Poly(f, t).degree() > 0}
    crit = sorted({r for f in polys for r in sp.real_roots(sp.Poly(f, t)) if r > 0}, key=float)
    pts = [sp.Integer(0)] + crit + [None]
    ivs = []
    for lo, hi in zip(pts, pts[1:]):
        mid = sp.Rational(sp.nsimplify(float(lo) + 1, rational=True)) if hi is None else \
            sp.nsimplify((sp.N(lo, 40) + sp.N(hi, 40)) / 2, rational=True, tolerance=1e-30)
        if all(sp.sign(p.subs(t, mid)) * sp.sign(q.subs(t, mid)) < 0 for p, q in pairs):
            ivs.append((lo, hi))
    rows = [(str(sign_core(p)), str(sign_core(q))) for p, q in pairs]
    return rows, ivs



TABLE2 = {
    'HAHABH': ['(5*s - 3, 7 - s)', '(1, 5*s - 3)', '(5*s - 3, s + 9)', '(1, 5*s - 3)', '(1, -1)', '(1, -1)'],
    'HAHBH': ['(7*s**2 - 10*s - 1, -(s**2 - 14*s + 1))', '(3*s - 1, 7*s**2 - 10*s - 1)',
              '(7*s**2 - 10*s - 1, s**2 + 10*s - 7)', '(s - 3, s**2 + 10*s - 7)',
              '(s**2 - 14*s + 1, s**2 + 10*s - 7)'],
    'HAHBHAHABH': ['(-7*s**3 + 67*s**2 - 53*s + 1, (s - 3)*(s**2 - 26*s + 5))',
                   '((7 - s)*(3*s - 1), -7*s**3 + 67*s**2 - 53*s + 1)',
                   '(-7*s**3 + 67*s**2 - 53*s + 1, -s**3 + 5*s**2 + 93*s - 41)',
                   '(-s**2 + 18*s - 13, -s**3 + 5*s**2 + 93*s - 41)',
                   '(-s**3 + 77*s**2 - 51*s - 1, -s**3 + 5*s**2 + 93*s - 41)',
                   '(3*s - 1, -s**3 + 77*s**2 - 51*s - 1)',
                   '(-s**3 + 77*s**2 - 51*s - 1, s**3 - 5*s**2 + 99*s - 23)',
                   '(3*s - 1, -s**3 + 77*s**2 - 51*s - 1)', '(3*s - 1, -(s - 3)**2)',
                   '(s**2 + 18*s + 1, -(s - 3)**2)'],
}
ENDS = {'HAHABH': (None, 5 * s - 3, 0.65905), 'HAHBHAHABH': (s**3 - 5 * s**2 - 93 * s + 41, s**3 - 77 * s**2 + 51 * s + 1, (0.58131, 0.69222)),
        'HAHBH': (s**2 + 10 * s - 7, 7 * s**2 - 10 * s - 1, (0.68108, 0.88970))}


u = sp.Symbol('u')
FW = {'HAHABH': 1 + u - 2 * u**2, 'HAHBHAHABH': 1 - 2 * u**2 - u**3 + 2 * u**4, 'HAHBH': 1 - 2 * u**2 + u**3,
      'HBHAHBHBAH': 1 - 2 * u + 3 * u**3 - 2 * u**4, 'HBHBAH': 1 - 3 * u + 2 * u**2}
REFL_U = {'A': (1, 0), 'B': (-1, 0), 'H': (u, 1 - u)}        # rho z = c conj(z) + d, coefficients in Z[u]


def F_word(W):
    """O_end = g_k(0) as a Laurent polynomial in u = e^{-2ia}; on |u| = 1, conj acts by u -> 1/u."""
    z = sp.Integer(0)
    for ch in reversed(W):
        c, d = REFL_U[ch]
        z = sp.expand(c * z.subs(u, 1 / u) + d)
    return z


def degenerate(Fw):
    """True iff Fw(u) = c u^m Fw(1/u) for a constant c and an integer m."""
    num, den = sp.fraction(sp.cancel(Fw / sp.expand(Fw.subs(u, 1 / u))))
    return len(sp.Poly(sp.expand(num), u).terms()) == 1 and len(sp.Poly(sp.expand(den), u).terms()) == 1


def same_up_to_positive(a, b):
    tt = sp.Symbol('t', positive=True)
    r = sp.cancel(sp.sympify(a).subs(s, tt**2).subs(sp.Symbol('t'), tt)
                  / sp.sympify(b).subs(s, tt**2).subs(sp.Symbol('t'), tt))
    return r.is_number and r > 0


def end_ok(r, poly, value):
    """tile end r (a root in t = tan a) is the smaller positive root of `poly` in s, with a = value.."""
    s2 = sp.nsimplify(r ** 2) if False else r ** 2
    roots = sorted(x for x in sp.real_roots(sp.Poly(poly, s)) if x > 0)
    return sp.simplify(s2 - roots[0]) == 0 and abs(math.atan(float(r)) - value) < 1e-5


def main():
    ivs = {}
    for W in ('HAHABH', 'HAHBHAHABH', 'HAHBH'):
        rows, iv = tile(W)
        ok = len(rows) == len(TABLE2[W]) and all(
            same_up_to_positive(sp.sympify(x)[0], sp.sympify(y)[0]) and same_up_to_positive(sp.sympify(x)[1], sp.sympify(y)[1])
            for x, y in zip(rows, TABLE2[W]))
        check(f'Table 2: {W}, all {len(rows)} sign-factor pairs as printed', ok)
        check(f'{W}: the tile is a single interval', len(iv) == 1)
        ivs[W] = iv[0]
    lo, hi = ivs['HAHABH']
    check('R_HAHABH = (0, a1), 5 tan^2 a1 = 3, a1 = 0.65905..',
          lo == 0 and sp.simplify(5 * hi**2 - 3) == 0 and abs(math.atan(float(hi)) - 0.65905) < 1e-5)
    lo, hi = ivs['HAHBHAHABH']
    check('R_HAHBHAHABH = (a2, a3): smaller positive roots of s^3-5s^2-93s+41 and s^3-77s^2+51s+1; '
          'a2 = 0.58131.., a3 = 0.69222..',
          end_ok(lo, ENDS['HAHBHAHABH'][0], 0.58131) and end_ok(hi, ENDS['HAHBHAHABH'][1], 0.69222))
    lo, hi = ivs['HAHBH']
    check('R_HAHBH = (a4, a5): s^2+10s-7 = 0, 7s^2-10s-1 = 0; a4 = 0.68108.., a5 = 0.88970..',
          end_ok(lo, ENDS['HAHBH'][0], 0.68108) and end_ok(hi, ENDS['HAHBH'][1], 0.88970))
    check('R_HAHBH is |cos 2a| < (sqrt 2 - 1)/2',
          all(abs(abs(float((1 - r**2) / (1 + r**2))) - (math.sqrt(2) - 1) / 2) < 1e-12 for r in ivs['HAHBH']))
    a1 = ivs['HAHABH'][1]
    a2, a3 = ivs['HAHBHAHABH']
    a4, a5 = ivs['HAHBH']
    check('cover: a2 < a1, a4 < a3, pi/2 - a3 < a5 (compared on tan a exactly)',
          bool(sp.simplify(a2 - a1) < 0) and bool(sp.simplify(a4 - a3) < 0) and bool(sp.simplify(1 / a3 - a5) < 0))
    for W, Wp in (('HAHABH', 'HBHBAH'), ('HAHBHAHABH', 'HBHAHBHBAH')):
        _, ivp = tile(Wp)
        lo, hi = ivs[W]
        exp_lo, exp_hi = 1 / hi, (None if lo == 0 else 1 / lo)
        ok = len(ivp) == 1 and sp.simplify(ivp[0][0] - exp_lo) == 0 and (
            (ivp[0][1] is None and exp_hi is None) or (exp_hi is not None and sp.simplify(ivp[0][1] - exp_hi) == 0))
        check(f'{Wp}: tile = pi/2 - R_{W}, computed directly', ok)
    words = ['HAHABH', 'HAHBHAHABH', 'HAHBH', 'HBHAHBHBAH', 'HBHBAH']
    check('none of the five words is an odd palindrome', not any(w == w[::-1] and len(w) % 2 for w in words))
    for W in words:
        Fw = F_word(W)
        check(f'{W}: F_W = {FW[W]}, non-degenerate', sp.expand(Fw - FW[W]) == 0 and not degenerate(Fw))
    check('control: 1 + 2u + u^2 (R = u^2) is flagged degenerate', degenerate(1 + 2 * u + u**2))
    # Lemma 6.3/6.4 at a rational cos 2a inside each tile
    for W, y in (('HAHABH', '1/2 + 1/5'), ('HAHBHAHABH', '19/50'), ('HAHBH', '1/10'),
                 ('HBHAHBHBAH', '-19/50'), ('HBHBAH', '-7/10')):
        yv = sp.Rational(sp.sympify(y))
        a = math.acos(float(yv)) / 2
        lo, hi = (float(sp.atan(x)) if x is not None else math.pi / 2 for x in tile(W)[1][0])
        K = m.t_field(f'{yv.q}*y - ({yv.p})', str(yv))
        found = None
        for c1 in ('AB', 'BA'):
            for c2 in ('AB', 'BA'):
                word = W + c1 + W[::-1] + c2
                if m.cyclic_ok(word):
                    cc = m.certify(K, word, mirror=False)
                    if cc['ok']:
                        found = cc
        check(f'Lemma 6.3/6.4: {W} at cos 2a = {yv} (a = {a:.4f} in its tile): cylinder of period 2k+4 = '
              f'{2 * len(W) + 4}, not perpendicular',
              lo < a < hi and found is not None and found['period'] == 2 * len(W) + 4 and not found['perpendicular'])
    done()


if __name__ == '__main__':
    main()
