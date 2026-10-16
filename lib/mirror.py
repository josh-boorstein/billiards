"""Mirror-direction cylinders of a right triangle with an algebraic, irrational angle: exact
certification of a periodic cylinder from its word alone.

The triangle.  `Tri(a)` has its right angle at O and the angle a at A.  We work in the field
K = Q(t), t = 2 cos a, and its CM extension K(w), w = e^{ia} (`numfield.CM`), with the triangle
scaled by t: O = 0, A = t, B = 2w - t (= 2i sin a), so that every vertex image lies in K(w).
y = cos 2a = (t^2 - 2)/2 is in K.  The sides are named by letters: `A` = OA, `B` = OB,
`H` = AB (the hypotenuse); a word is the cyclic sequence of sides a path hits.

The certificate (`certify`).  Let W = w_1 ... w_n be a cyclic word with no two equal letters
adjacent (cyclically).  Unfold it: g_0 = id, g_i = g_{i-1} o (reflection in the side w_i), copies
T_i = g_i(T), crossed edges e_i = g_{i-1}(side w_i).  W is the word of a periodic cylinder iff

  (1) g_n is a translation, z -> z + beta, beta != 0  (exact: linear part 1, orientation kept);
  (2) some line parallel to beta separates the two endpoints of every e_i strictly.

Writing F(z) := Im(conj(beta) z)/sin a (an element of K), (2) says the open intervals between
F(P) and F(Q), over the edges (P, Q), have a common point; their intersection (lo, hi) is the
set of offsets of the cylinder's lines, and the cylinder has width (hi - lo) sin a / |beta| in
the scaled frame.  Sufficiency is the corridor argument: consecutive edges are two distinct
sides of the same copy, a line through no vertex crossing both meets that copy in the chord
between them, and since T_n = T_0 + beta (the edge e_n is the side w_n of T_n, whose
translate by -beta is the side of T_0 through which the line entered) the crossings are met in
the order e_1, ..., e_n, e_1 + beta, ... -- the order is forced, because the crossing of e_n lies
|beta| beyond that of e_n - beta.  Every line with offset in (lo, hi) is therefore a periodic
billiard path following W, and no line outside is.

The direction of beta is identified exactly: it is parallel to w^e (direction e*a) or to
i w^e (direction pi/2 + e*a) for a unique integer e when a/pi is irrational, which names the
mirror class -- e even: perpendicular to OB resp. OA; e odd: parallel resp. perpendicular to H.
A path of the cylinder meets a side at a right angle iff beta is perpendicular to some crossed
edge, Re(conj(beta)(Q - P)) = 0, decided exactly; the paths near the corners O_end are not an
issue because no vertex lies on any line of the cylinder.

No decision here is made on a floating-point value.  Floats appear only to pick candidate
exponents e (each then checked exactly) and in the reported angles.
"""
import math
from fractions import Fraction

import mpmath
import sympy as sp

from numfield import CM, RealField

Y, T, Z = sp.symbols('y t z')
LETTERS = 'ABH'
CLASSES = ('perp OB', 'perp OA', 'par H', 'perp H')


# ------------------------------------------------------------------------------ the field

def t_field(minpoly_y, y_approx):
    """K = Q(t), t = 2 cos a, for the triangle with cos 2a = y.

    `minpoly_y` is an integer polynomial in y (a sympy expression or a string) vanishing at
    y; `y_approx` locates y among its roots (a float, or a Fraction / str for exactness).
    t is the positive root of t^2 = 2y + 2, and K is cut out by the factor of
    f((t^2 - 2)/2) that vanishes there.
    """
    f = sp.Poly(sp.sympify(minpoly_y), Y)
    yv = float(Fraction(y_approx)) if isinstance(y_approx, str) else float(y_approx)
    assert -1 < yv < 1
    tv = math.sqrt(2 * yv + 2)
    Pt = sp.Poly(sp.expand(f.as_expr().subs(Y, (T ** 2 - 2) / 2) * 2 ** f.degree()), T)
    facs = [g for g, _ in sp.factor_list(Pt)[1]]
    g = min(facs, key=lambda g: abs(complex(sp.Poly(g, T).eval(tv))))
    cs = [int(c) for c in reversed(sp.Poly(g, T).all_coeffs())]
    return RealField(cs, tv, name='t')


def minpoly_y_of(K):
    """The minimal polynomial of y = (t^2 - 2)/2 over Q, as a primitive integer sympy Poly."""
    tpoly = sum(int(c) * T ** i for i, c in enumerate(K.poly))
    # resultant route (exact, no root objects): y = (t^2 - 2)/2
    r = sp.resultant(sp.Poly(tpoly, T).as_expr(), 2 * Y - (T ** 2 - 2), T)
    fl = [g for g, _ in sp.factor_list(sp.Poly(r, Y))[1]]
    yv = (K.approx(30) ** 2 - 2) / 2
    g = min(fl, key=lambda g: abs(complex(sp.Poly(g, Y).eval(float(yv)))))
    return sp.Poly(g, Y)


def angle_is_irrational(minpoly_y):
    """a/pi is irrational, for cos 2a a root of the irreducible `minpoly_y`.

    a/pi in Q iff 2 cos 2a = 2 cos(2 pi k/n) for some n, i.e. iff the minimal polynomial of
    2y is psi_n, the minimal polynomial of 2 cos(2 pi/n).  That needs 2y to be an algebraic
    integer and deg = phi(n)/2 (or 1 for n <= 2), which leaves finitely many n to compare.
    """
    f = sp.Poly(sp.sympify(minpoly_y), Y)
    assert f.is_irreducible
    g = sp.Poly(sp.expand(f.as_expr().subs(Y, Z / 2) * 2 ** f.degree()), Z)
    g = sp.Poly(g.as_expr() / g.LC(), Z)
    if not all(c.is_integer for c in g.all_coeffs()):
        return True
    d = g.degree()
    for n in range(1, 8 * d * d + 7):
        phi = sp.totient(n)
        if (phi // 2 if n > 2 else 1) != d:
            continue
        if sp.Poly(sp.minimal_polynomial(2 * sp.cos(2 * sp.pi / n), Z), Z) == g:
            return False
    return True


# --------------------------------------------------------------------------- the triangle

class Triangle:
    """Tri(a) in the t-scaled frame; maps are (lin, trans, anti): z -> lin*(conj z if anti) + trans."""
    SIDE = {'A': ('O', 'A'), 'B': ('O', 'B'), 'H': ('A', 'B')}

    def __init__(self, K):
        self.K = K
        self.C = CM(K, K.gen)
        C, w, t = self.C, self.C.w, K.gen
        self.w, self.t = w, t
        self.V = {'O': C.zero, 'A': C.of(t), 'B': 2 * w - t}
        u = w.conj() * w.conj()                                  # e^{-2ia}
        self.RHO = {'A': (C.one, C.zero), 'B': (-C.one, C.zero), 'H': (u, (C.one - u) * t)}
        self.ident = (C.one, C.zero, 0)

    def then(self, g, ch):
        a, b, f = g
        c, d = self.RHO[ch]
        if f == 0:
            return (a * c, a * d + b, 1)
        return (a * c.conj(), a * d.conj() + b, 0)

    @staticmethod
    def apply(g, z):
        a, b, f = g
        return a * (z if f == 0 else z.conj()) + b

    def unfold(self, word):
        """[(g_{i-1}, edge e_i endpoints)] for i = 1..n, and g_n."""
        g, out = self.ident, []
        for ch in word:
            P, Q = self.SIDE[ch]
            out.append((self.apply(g, self.V[P]), self.apply(g, self.V[Q])))
            g = self.then(g, ch)
        return out, g

    def wpow(self, e):
        return self.w ** e if e >= 0 else self.w.conj() ** (-e)


def primitive(word):
    n = len(word)
    return all(word != word[k:] + word[:k] for k in range(1, n) if n % k == 0)


def cyclic_ok(word):
    return len(word) >= 2 and set(word) <= set(LETTERS) and all(
        word[i] != word[(i + 1) % len(word)] for i in range(len(word)))


# ---------------------------------------------------------------------------- certificate

def direction_class(tri, beta, nmax):
    """(class, e) with beta parallel to w^e or i w^e, |e| <= nmax, checked exactly; None if none."""
    K, C, w, t = tri.K, tri.C, tri.w, tri.t
    a = math.acos(float(t.mp(30)) / 2)
    z = beta.complex(30)
    phi = float(mpmath.arg(z))
    ivec = w - t / 2                                              # i sin a: a positive multiple of i
    for itype in (0, 1):
        base = phi - (math.pi / 2 if itype else 0.0)
        kmax = int((nmax + 2) * a / math.pi) + 2
        for k in range(-kmax, kmax + 1):
            ef = (base + k * math.pi) / a
            e = round(ef)
            if abs(ef - e) > 1e-6 or abs(e) > nmax:
                continue
            d = tri.wpow(e) * (ivec if itype else C.one)
            if (d.conj() * beta).im_over_sin().iszero():
                cls = ('perp OA' if itype else 'perp OB') if e % 2 == 0 else ('perp H' if itype else 'par H')
                return cls, e
    return None


def certify(K, word, mirror=True):
    """Certify `word` as the word of a periodic cylinder of Tri(a) (see the module docstring).
    With `mirror=False` the direction need not be a mirror direction (cls, e are then None
    unless it is one): the cylinder is still certified, e.g. a corner word's at a generic angle.

    Returns a dict: ok (the cylinder exists, in a mirror direction), period, cls, e,
    perpendicular (some path meets a side at a right angle), lo, hi (offsets, in K), width2
    (the squared width in the frame OA = 1, an element of K), and float summaries.
    """
    out = {'ok': False, 'period': len(word)}
    if not cyclic_ok(word):
        out['why'] = 'not a cyclic billiard word'
        return out
    out['primitive'] = primitive(word)
    tri = Triangle(K)
    edges, g = tri.unfold(word)
    lin, beta, anti = g
    if anti != 0 or not (lin - tri.C.one).iszero() or beta.iszero():
        out['why'] = 'does not unfold to a nonzero translation'
        return out
    dc = direction_class(tri, beta, len(word) + 4)
    if dc is None and mirror:
        out['why'] = 'translation not in a mirror direction'
        return out
    out['cls'], out['e'] = dc if dc else (None, None)
    F = lambda z: (beta.conj() * z).im_over_sin()
    lo = hi = None
    perp = False
    min_ang = math.inf
    bz = beta.complex(20)
    for P, Q in edges:
        fp, fq = F(P), F(Q)
        s = (fp - fq).sign()
        if s == 0:
            out['why'] = 'an edge parallel to the translation'
            return out
        a, b = (fq, fp) if s > 0 else (fp, fq)
        lo = a if lo is None or (a - lo).sign() > 0 else lo
        hi = b if hi is None or (b - hi).sign() < 0 else hi
        if (beta.conj() * (Q - P)).re().iszero():
            perp = True
        ez = (Q - P).complex(20)
        cosang = abs(float(mpmath.re(mpmath.conj(bz) * ez))) / float(abs(bz) * abs(ez))
        min_ang = min(min_ang, math.asin(min(1.0, cosang)))      # angle to the edge's normal
    if (hi - lo).sign() <= 0:
        out['why'] = 'no line separates every edge (empty strip)'
        return out
    t = tri.t
    sin2 = 1 - t * t / 4
    width2 = (hi - lo) * (hi - lo) * sin2 / (beta.norm() * t * t)
    out.update(ok=True, perpendicular=perp, lo=lo, hi=hi, beta=beta, width2=width2,
               width=math.sqrt(float(width2.mp(30))), min_angle_to_normal=min_ang, tri=tri)
    return out


def launch_interval(cert, side='B'):
    """The cylinder's trace on the side `side` of the copy T_0 it leaves (W must end in `side`):
    the interval of parameters lam in (0, 1) along that side from its first to its second
    endpoint (O -> B for `B`), as a pair of elements of K, clipped to [0, 1]."""
    tri, lo, hi, beta = cert['tri'], cert['lo'], cert['hi'], cert['beta']
    P, Q = (tri.V[v] for v in Triangle.SIDE[side])
    F = lambda z: (beta.conj() * z).im_over_sin()
    fP, fQ = F(P), F(Q)
    l1, l2 = (lo - fP) / (fQ - fP), (hi - fP) / (fQ - fP)
    if (l1 - l2).sign() > 0:
        l1, l2 = l2, l1
    K = tri.K
    l1 = K.zero if l1.sign() < 0 else l1
    l2 = K.one if (l2 - 1).sign() > 0 else l2
    return l1, l2


# ------------------------------------------------------------------------- word proposals

def trace_word(y, theta, side, lam, nmax, dps=60):
    """PROPOSE a word: the billiard path of Tri(a), cos 2a = y, leaving the point at parameter
    `lam` along side `side` in direction `theta`, traced at `dps` digits until it first returns
    to its starting point with its starting direction, within 10^(-dps/2), or `nmax` bounces.
    Returns the word (ending in `side`) or None.  Used only to propose words; nothing is
    concluded from it -- `certify` decides."""
    with mpmath.workdps(dps):
        y = mpmath.mpf(y) if not isinstance(y, Fraction) else mpmath.mpf(y.numerator) / y.denominator
        a = mpmath.acos(y) / 2
        V = {'O': mpmath.mpc(0), 'A': mpmath.mpc(1), 'B': mpmath.mpc(0, mpmath.tan(a))}
        P, Q = (V[v] for v in Triangle.SIDE[side])
        p = P + lam * (Q - P)
        p0 = p
        d = mpmath.expj(theta)
        d0 = d
        last, word = side, ''
        tol = mpmath.mpf(10) ** (-dps // 2)
        for _ in range(nmax):
            best = None
            for ch, (i, j) in Triangle.SIDE.items():
                if ch == last:
                    continue
                Pp, Qq = V[i], V[j]
                e = Qq - Pp
                den = mpmath.im(mpmath.conj(d) * e)
                if den == 0:
                    continue
                tt = mpmath.im(mpmath.conj(e) * (Pp - p)) / mpmath.im(mpmath.conj(e) * d)
                uu = mpmath.im(mpmath.conj(d) * (p - Pp)) / mpmath.im(mpmath.conj(d) * e)
                if tt > tol and -tol <= uu <= 1 + tol and (best is None or tt < best[0]):
                    best = (tt, ch, e)
            if best is None:
                return None
            tt, ch, e = best
            p = p + tt * d
            e = e / abs(e)
            d = e * e * mpmath.conj(d)
            word += ch
            last = ch
            if ch == side and abs(d - d0) < tol and abs(p - p0) < tol:
                return word
        return None
