"""Proposition 8.2's containment step -- the one computer-assisted step of the polygon note: in
Delta_1 = (3, 6m + 2, 12m + 4) and Delta_2 = (3, 6m + 1, 12m + 5), N = 18m + 9, x = pi/N, the k-th swapped prong
(k = 0, ..., m - 1) crosses exactly the side copies its word says, for every m >= 1 at once.  Exact: sympy over
Q(sqrt 3), no floating point.

THE MODEL (the note's proof of Proposition 8.2).  The k-th prong leaves V_2 in direction
psi_k = -5pi/6 + (8 + 3k)x (Delta_1), -5pi/6 + (7 + 3k)x (Delta_2), with word
w_{m,k} = (02)^{2m-1-k} 0 . 1 . (02)^{2m+2+k} . 1 . (02)^k 0 (sides 0 = V0V1, 1 = V1V2, 2 = V2V0).  Its development is
three stars of copies about images O_1, O_2, O_3 of V_0 (the second the reflection of the first in e_1, the third of
the second in e_2), and with 12mx = 2pi/3 - 6x and y := kx every angle is a fixed combination of pi/6, x and y, on the
region 0 < x <= pi/27, 0 <= y <= Y(x) := (pi - 27x)/18.  The prong's line meets each required side copy iff the
SIDE FUNCTIONS side(z) := Im((z - S) e^{-i psi}) of the copy's two ends have opposite signs; fourteen conditions
(each star's extreme crossed rays, the cuts of e_1 and e_2, and two used only at k >= 1).

WHAT IS CHECKED.
  (D) the direction identity: Im((E - S) e^{-i psi}) is the zero Laurent polynomial in q = e^{ix}, r = e^{iy}
      (the endpoint E is the image of V_2; the note proves this by hand, and it is re-checked here);
  (H) every side function used is exactly a(x) cos 3y + b(x) sin 3y (only r^{+-3} occur), so on an interval of
      3y-length < pi it has constant sign iff its two end values agree in sign;
  (S) each end function -- f(x, 0) and f(x, Y(x)) on x <= pi/27, and f(x, x), f(x, Y(x)) on x <= pi/45 for the two
      conditions used only at k >= 1 -- is (A(t) + sqrt3 B(t))/(1 + t^2)^d in t = tan(x/4); the norm A^2 - 3B^2 in
      Q[t] (degree at most 92), its factor t^j removed, has no root in (0, 3/100] (resp. (0, 7/400]) by an exact
      Sturm count, except side(O_1) of Delta_1 at y = 0, whose norm vanishes at the conjugate point x = pi/30 and
      which has no root there over Q(sqrt 3) (exact count); and the signs read exactly at t = 1/100 satisfy all
      fourteen conditions in both families.
(3/100 > tan(pi/108) and 7/400 > tan(pi/180) are asserted.)  A few minutes.
"""
import sys
import time

import common
from common import check, done

import sympy as sp

x, y = sp.symbols('x y', real=True)
q, r, w = sp.symbols('q r w')
t = sp.symbols('t', real=True)
pi, I, S3 = sp.pi, sp.I, sp.sqrt(3)


def model(fam):
    A = 3*x
    B, C = (pi/3 - x, 2*pi/3 - 2*x) if fam == 'D1' else (pi/3 - 2*x, 2*pi/3 - x)
    psi = -5*pi/6 + (8 if fam == 'D1' else 7)*x + 3*y
    rC, rB = sp.sin(C), sp.sin(B)
    h2 = 2*rB*rC

    def ang1(c0, c1, c2):        # star-1 angle of ray kappa = c0 + c1 m + c2 k  (12 m x = 2pi/3 - 6x, y = k x)
        return 3*x*c0 + c1*(2*pi/3 - 6*x)/4 + 3*c2*y
    rad = lambda c0: rC if c0 % 2 == 0 else rB       # c1 in 4Z, c2 even: parity of kappa = parity of c0
    phi1 = ang1(2, -4, 2) - (pi/2 - B)               # e1: V1 ray -(4m-2-2k), V2 ray -(4m-1-2k)
    sig = lambda a: 2*phi1 + pi - a                  # star 2 = reflection of star 1 in e1
    phi2 = sig(ang1(6, 0, 4)) + (pi/2 - B)           # e2: V1 ray 6+4k, V2 ray 5+4k
    tau = lambda a: 2*(phi2 - phi1) + a              # star 3 = reflection of star 2 in e2
    O = {1: sp.Integer(0), 2: h2*sp.exp(I*phi1)}
    O[3] = O[2] + h2*sp.exp(I*phi2)
    S = rB*sp.exp(I*A)
    E = O[3] + rB*sp.exp(I*tau(ang1(7, 0, 6)))
    rot = {1: lambda a: a, 2: sig, 3: tau}
    spec = {'s1:0': (1, (0, 0, 0)), 's1:-(4m-2-2k)': (1, (2, -4, 2)), 's1:-1': (1, (-1, 0, 0)),
            's1:-(4m-3-2k)': (1, (3, -4, 2)), 's1:-(4m-1-2k)': (1, (1, -4, 2)),
            's2:-(4m-2-2k)': (2, (2, -4, 2)), 's2:4+4k': (2, (4, 0, 4)), 's2:-(4m-3-2k)': (2, (3, -4, 2)),
            's2:5+4k': (2, (5, 0, 4)), 's2:6+4k': (2, (6, 0, 4)),
            's3:6+4k': (3, (6, 0, 4)), 's3:6+6k': (3, (6, 0, 6)), 's3:7+4k': (3, (7, 0, 4)),
            's3:5+6k': (3, (5, 0, 6))}
    P = {n: O[st] + rad(c[0])*sp.exp(I*rot[st](ang1(*c))) for n, (st, c) in spec.items()}
    P.update({'O1': O[1], 'O2': O[2], 'O3': O[3]})
    return psi, S, E, P


# (name, point a, point b, k>=1 only): the line crosses segment ab
CONDS = [('star1 ray 0', 'O1', 's1:0', 0), ('star1 ray -(4m-2-2k)', 'O1', 's1:-(4m-2-2k)', 0),
         ('star1 ray -1', 'O1', 's1:-1', 0), ('star1 ray -(4m-3-2k)', 'O1', 's1:-(4m-3-2k)', 0),
         ('e1 cut', 's1:-(4m-2-2k)', 's1:-(4m-1-2k)', 0),
         ('star2 ray -(4m-2-2k)', 'O2', 's2:-(4m-2-2k)', 0), ('star2 ray 4+4k', 'O2', 's2:4+4k', 0),
         ('star2 ray -(4m-3-2k)', 'O2', 's2:-(4m-3-2k)', 0), ('star2 ray 5+4k', 'O2', 's2:5+4k', 0),
         ('e2 cut', 's2:5+4k', 's2:6+4k', 0),
         ('star3 ray 6+4k', 'O3', 's3:6+4k', 0), ('star3 ray 6+6k', 'O3', 's3:6+6k', 0),
         ('star3 ray 7+4k', 'O3', 's3:7+4k', 1), ('star3 ray 5+6k', 'O3', 's3:5+6k', 1)]


def laurent_qr(e):
    e = sp.expand(sp.expand(e).rewrite(sp.exp), power_exp=True)
    e = sp.expand(e.subs({sp.exp(I*x): q, sp.exp(I*y): r}))
    e = sp.expand(e.subs({sp.exp(-I*x): 1/q, sp.exp(-I*y): 1/r}))
    assert not (e.free_symbols & {x, y}), 'x or y survived'
    return e


def harmonic(e):
    """e (Laurent in q, r) -> (c3, cm3) with e == c3 r^3 + cm3 r^-3 exactly, else AssertionError."""
    E20 = sp.Poly(sp.expand(e*r**20), r)
    pw = {mon[0] - 20: c for mon, c in E20.terms()}
    assert set(pw) <= {3, -3}, sorted(pw)
    return pw.get(3, 0), pw.get(-3, 0)


def end_fn(c3, cm3, which):
    """one-variable end function as a Laurent expression in w = e^{iu}, u = x/2."""
    r3 = {'y=0': sp.Integer(1), 'y=Y': sp.exp(I*pi/6)*w**-9, 'y=x': w**6}[which]
    e = sp.expand(c3.subs(q, w**2)*r3 + cm3.subs(q, w**2)/r3)
    return sp.expand(e.rewrite(sp.cos).subs(sp.exp(I*pi/6), (S3 + I)/2)) if e.has(sp.exp) else e


def to_t(e):
    """Laurent poly in w (coefficients in Q(sqrt3, i)), real on |w| = 1 -> (A, B) in Q[t]:
    (1+t^2)^D e = A + sqrt3 B, w = (1+it)/(1-it)."""
    e = sp.expand(e)
    terms = sp.Poly(sp.expand(e*w**200), w).terms()
    ks = [mon[0] - 200 for mon, _ in terms]
    D = max(abs(k) for k in ks)
    tot = 0
    for (mon, c) in terms:
        k = mon[0] - 200
        tot += c*(1 + I*t)**(D + k)*(1 - I*t)**(D - k)
    tot = sp.expand(tot)
    re, im = sp.expand(sp.re(tot)), sp.expand(sp.im(tot))
    assert sp.expand(im) == 0, 'not real'
    A, B = 0, 0
    for (k,), co in sp.Poly(re, t).terms():
        co = sp.expand(co)
        b = co.coeff(S3)
        a = sp.expand(co - b*S3)
        assert a.is_Rational and b.is_Rational, co
        A += a*t**k
        B += b*t**k
    return sp.Poly(A, t, domain='QQ'), sp.Poly(B, t, domain='QQ')


def strip_t(Pn):
    k = 0
    while not Pn.is_zero and Pn.eval(0) == 0:
        Pn = sp.Poly(sp.cancel(Pn.as_expr()/t), t, domain=Pn.domain)
        k += 1
    return Pn, k


def sign_exact(A, B, tv):
    a, b = A.eval(tv), B.eval(tv)
    if a >= 0 and b >= 0:
        return 0 if a == 0 and b == 0 else 1
    if a <= 0 and b <= 0:
        return -1
    return (1 if a > 0 else -1)*(1 if a*a > 3*b*b else -1)


def certify(e, T, tread):
    """sign of the end function on (0, T] in t, exact; None if it has a root there."""
    A, B = to_t(e)
    Nm, _ = strip_t(A**2 - 3*B**2)
    hits = Nm.count_roots(0, T)
    fb = None
    if hits:
        F = sp.Poly(A.as_expr() + S3*B.as_expr(), t, extension=S3)
        F, _ = strip_t(F)
        fb = F.count_roots(0, T)
    ok = hits == 0 or fb == 0
    return (sign_exact(A, B, tread) if ok else None), {'deg_N': Nm.degree(), 'norm_hits': int(hits),
                                                       'fallback': fb}


def main():
    T_ALL, T_K1, TREAD = sp.Rational(3, 100), sp.Rational(7, 400), sp.Rational(1, 100)
    check('3/100 > tan(pi/108), 7/400 > tan(pi/180), and x = 4 atan(1/100) < pi/45',
          sp.tan(pi/108).evalf(30) < T_ALL and sp.tan(pi/180).evalf(30) < T_K1
          and 4*sp.atan(TREAD).evalf(30) < (pi/45).evalf(30))
    for fam, name in (('D1', 'Delta_1 = (3, 6m+2, 12m+4)'), ('D2', 'Delta_2 = (3, 6m+1, 12m+5)')):
        t0 = time.time()
        psi, S, E, P = model(fam)
        Z = sp.expand((E - S)*sp.exp(-I*psi))
        dz = laurent_qr((Z - sp.conjugate(Z))/(2*I))
        check(f'{name}: (D) Im((E - S) e^(-i psi)) is the zero Laurent polynomial', sp.expand(sp.expand_complex(dz)) == 0)
        need = {}
        for _, a, b, k1 in CONDS:
            for n in (a, b):
                need.setdefault(n, set()).update({('y=0', 0), ('y=Y', 0)} if not k1 else {('y=x', 1), ('y=Y', 1)})
        sides, harm_ok, hits = {}, True, []
        for n, ends in need.items():
            z = P[n]
            wv = sp.expand((z - S)*sp.exp(-I*psi))
            try:
                c3, cm3 = harmonic(laurent_qr((wv - sp.conjugate(wv))/(2*I)))
            except AssertionError:
                harm_ok = False
                continue
            sides[n] = {}
            for which, k1 in sorted(ends):
                sg, info = certify(end_fn(c3, cm3, which), T_K1 if k1 else T_ALL, TREAD)
                key = f'{which}|k1' if k1 else which
                sides[n][key] = dict(info, sign=sg)
                if info['norm_hits']:
                    hits.append((n, key, info['fallback']))
        check(f'{name}: (H) every side function is exactly a(x) cos 3y + b(x) sin 3y', harm_ok)
        degs = [v['deg_N'] for s in sides.values() for v in s.values()]
        exc = ', except side(O_1) at y = 0, which has none over Q(sqrt 3)' if fam == 'D1' else ''
        check(f'{name}: (S) every norm has degree <= 92 (max {max(degs)}), and no root on its range{exc} '
              f'(norm hits: {hits})', max(degs) <= 92 and (hits == [] if fam == 'D2' else hits == [('O1', 'y=0', 0)]))
        conds = []
        for c, a, b, k1 in CONDS:
            ends = ('y=x|k1', 'y=Y|k1') if k1 else ('y=0', 'y=Y')
            sa = [sides[a][e]['sign'] for e in ends]
            sb = [sides[b][e]['sign'] for e in ends]
            const = None not in sa + sb and sa[0] == sa[1] and sb[0] == sb[1]
            conds.append((c, bool(const and sa[0]*sb[0] < 0)))
        check(f'{name}: (S) the signs at t = 1/100 satisfy all fourteen conditions '
              f'(failing: {[c for c, ok in conds if not ok]}) [{time.time() - t0:.0f} s]',
              len(conds) == 14 and all(ok for _, ok in conds))
    done()


if __name__ == '__main__':
    main()
