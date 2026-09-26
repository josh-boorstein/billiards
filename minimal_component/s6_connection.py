"""Section 6, "What the direction does carry": the saddle connection and its arithmetic.

Statements checked:
  * The billiard path leaving O at 12 degrees to L2 reaches the vertex A after 13
    reflections (it arrives at A on its 14th contact with the boundary): the word is
    proposed by the tracer and the arrival is verified EXACTLY in Q(zeta_60), every
    earlier contact strictly inside its side.
  * Leaves launched parallel to it within 1e-25 on either side show the minimal component
    on one side (they reach Sigma_Hbar, which lies in M) and a closed orbit, i.e. a
    cylinder, on the other.
  * Its arrival condition, derived by formal development of the word, is the vanishing
    sine sum 2 sin 12 + 2 sin 24 + sin 36 = sin 60 + sin 72 (degrees); that identity and
    the boundary-5 identity sin 12 + 2 sin 24 + sin 36 = sin 48 + sin 60 both hold exactly.
  * Neither identity follows from the 3-family (with the sign relations) alone, nor from
    the 5-family alone: a combination of the first kind still vanishes with 12 degrees
    replaced by 60, of the second kind with 12 replaced by 36, and both identities fail at
    both angles.  (Checked also algebraically: as elements of Z[x]/(x^30 + 1), x =
    e^{i pi/30}, neither vanishes at a primitive 12th or 20th root of unity.)
  * The connection persists under deformation of the triangle, but its direction moves
    at a rate that is not an even integer, so it leaves the directions +-2k alpha + m pi
    of the horizontal flow at once.
"""
import sympy as sp
from mpmath import mp, mpf, pi, sin, cos, diff, findroot, nstr, radians

from t0 import T, F, M_H, check, done
from develop import develop
from tracer import Billiard

mp.dps = 60


def arrival_sines(word, alpha_deg=48, phi_deg=12):
    """The arrival condition Im(e^{-i phi} G(A)) sin(alpha) = 0 of a path from O in
    direction phi through the copy G(A) of A, G the development of `word`, as an integer
    sine sum {angle in degrees: coefficient}, reduced to angles in (0, 90]."""
    a, k, b = develop(word)
    terms = {}

    def add(deg, c):
        deg %= 360
        sg = 1
        if deg >= 180:
            deg, sg = deg - 180, -1
        if deg > 90:
            deg = 180 - deg
        if deg % 180:
            terms[deg] = terms.get(deg, 0) + sg * c
    psi = lambda h, q: 2 * h * alpha_deg + 180 * q
    # A = e^{i alpha}/sin(alpha) and L = cos(alpha)/sin(alpha); 2 cos(alpha) sin(x) =
    # sin(x + alpha) + sin(x - alpha)
    add(psi(*a) + (alpha_deg if k == 0 else -alpha_deg) - phi_deg, 2)
    for (h, q), c in b.items():
        add(psi(h, q) - phi_deg + alpha_deg, c)
        add(psi(h, q) - phi_deg - alpha_deg, c)
    return {d: c for d, c in sorted(terms.items()) if c}


def leaf_class(B, px, py, vx, vy, cap=3000):
    """'M' if the leaf reaches Sigma_Hbar (L1, outgoing 11pi/15), 'closed' if it comes
    back to its launch point in its launch direction, else '?'."""
    tx, ty = cos(11 * pi / 15), sin(11 * pi / 15)
    x0, y0, u0, w0 = px, py, vx, vy
    tol = mpf(10) ** -30
    for j in range(1, cap + 1):
        st = B.step(px, py, vx, vy)
        if st is None:
            return 'lost', j
        t, side = st
        if abs(vx - u0) < tol and abs(vy - w0) < tol:
            dx, dy = x0 - px, y0 - py
            tt = dx * vx + dy * vy
            if 0 < tt < t and abs(dx * vy - dy * vx) < tol:
                return 'closed', j
        px, py = px + t * vx, py + t * vy
        vx, vy = B.reflect(side, vx, vy)
        if side == '1' and 0 < py < 1 and abs(vx - tx) < tol and abs(vy - ty) < tol:
            return 'M', j
    return '?', cap


def fam_value(S, unit_deg):
    """The sine sum S (degrees, multiples of 12) with 12 degrees replaced by unit_deg."""
    return sum(c * sin(radians(d // 12 * unit_deg)) for d, c in S.items())


def classify(S):
    """Does the sine sum vanish at a primitive 12th / 20th / 60th root (x^30 = -1)?"""
    x = sp.Symbol('x')
    acc = {}
    for d, c in S.items():
        e = d // 6                                     # x = e^{i pi/30}: 6 degrees
        for ee, s2 in ((e, 1), (-e, -1)):
            ee %= 60
            s3 = 1
            if ee >= 30:
                ee, s3 = ee - 30, -1
            acc[ee] = acc.get(ee, 0) + c * s2 * s3
    poly = sum(v * x ** e for e, v in acc.items())
    z = lambda n: sp.expand(sp.rem(poly, sp.cyclotomic_poly(n, x), x)) == 0
    return z(60), z(12), z(20)


def main():
    print('Section 6 -- the saddle connection and its arithmetic\n')
    B = Billiard(mp.pi * 8 / 30)
    orb = B.orbit(mpf(0), mpf(0), pi / 15, 40)
    hit = next(j for j, (s, x, y) in enumerate(orb, 1)
               if abs(x - B.L) < mpf(10) ** -40 and abs(y - 1) < mpf(10) ** -40)
    word = ''.join(s for s, _, _ in orb[:hit])
    print(f'  tracer: O at 12 degrees arrives at A on contact {hit}; word {word}')
    st = T.unfold(word, 2, start=T.launch_point(F.zero, F.zero))
    Ax, Ay = st[-1][4]
    check('exactly: the path arrives at A = (cot alpha, 1)', Ax == T.L and Ay == 1)
    interior = all(f1.is_zero() and f0.sign() > 0
                   for s_ in st[:-1] for (_, _, f0, f1) in T.constraints(s_))
    check('every earlier contact strictly inside its side', interior)
    check(f'so A is reached after {hit - 1} reflections (expect 13)', hit - 1 == 13)

    # leaves either side
    ux, uy = cos(pi / 15), sin(pi / 15)
    px, py = mpf(0), mpf(0)
    mids = []
    for s, x, y in orb[:hit]:
        mids.append(((px + x) / 2, (py + y) / 2, ux, uy))
        ux, uy = B.reflect(s, ux, uy)
        px, py = x, y
    sides_ok = True
    for kk in (0, 4, 8, 12):
        mx, my, vx, vy = mids[kk]
        nx, ny = -vy, vx
        e = mpf('1e-25')
        l = leaf_class(B, mx + e * nx, my + e * ny, vx, vy)[0]
        r = leaf_class(B, mx - e * nx, my - e * ny, vx, vy)[0]
        print(f'  segment {kk + 1:2d}: left {l}, right {r}')
        sides_ok &= {l, r} == {'M', 'closed'}
    check('within 1e-25: the minimal component on one side, a cylinder on the other', sides_ok)

    S_conn = arrival_sines(word[:-1])
    S_id = {12: 2, 24: 2, 36: 1, 60: -1, 72: -1}
    S_b5 = {12: 1, 24: 2, 36: 1, 48: -1, 60: -1}
    print(f'  arrival condition (degrees: coefficient) = {S_conn}')
    check('the arrival condition is twice 2 sin 12 + 2 sin 24 + sin 36 - sin 60 - sin 72',
          S_conn == {d: 2 * c for d, c in S_id.items()})
    ex = lambda S: sum((F.sin(d // 6) * c for d, c in S.items()), F.zero)
    check('2 sin 12 + 2 sin 24 + sin 36 = sin 60 + sin 72, exactly', ex(S_id).is_zero())
    check('sin 12 + 2 sin 24 + sin 36 = sin 48 + sin 60, exactly', ex(S_b5).is_zero())
    for name, S in (('connection', S_id), ('boundary 5', S_b5)):
        v60, v36 = fam_value(S, 60), fam_value(S, 36)
        z60, z12, z20 = classify(S)
        print(f'  {name}: at 60 degrees {nstr(v60, 6)}, at 36 degrees {nstr(v36, 6)}')
        check(f'{name}: fails at both substitutions -> needs both the 3- and 5-families',
              abs(v60) > 0.1 and abs(v36) > 0.1)
        check(f'{name}: vanishes at zeta_60, not at zeta_12 nor zeta_20', z60 and not z12 and not z20)

    # the deformed connection: F(alpha, phi) = 0 along the fixed word
    a, k, b = develop(word[:-1])

    def Fc(al, ph):
        tot = 2 * sin((2 * a[0] * al + a[1] * pi) + (al if k == 0 else -al) - ph)
        for (h, q), c in b.items():
            ps = 2 * h * al + q * pi
            tot += c * (sin(ps - ph + al) + sin(ps - ph - al))
        return tot
    a0, p0 = 4 * pi / 15, pi / 15
    check('F(4pi/15, pi/15) = 0', abs(Fc(a0, p0)) < mpf(10) ** -50)
    rate = -diff(lambda t: Fc(t, p0), a0) / diff(lambda t: Fc(a0, t), p0)
    print(f'  d phi / d alpha along the connection = {nstr(rate, 12)}')
    check('the rate is not an even integer', abs(rate / 2 - mp.nint(rate / 2)) > mpf('1e-3'))
    done()


if __name__ == '__main__':
    main()
