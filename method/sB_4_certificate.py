"""III Appendix B.4: the certificate behind Theorem 7.11 (the (4, 3) ladder) -- the fourteen inequalities (I1)-(I14)
of Appendix B.2 (c)-(e) at every rung m >= 3, and the prong clause at the two rungs m = 1, 2.

Setting (III B.1): P = 4, Q = 4m + 3, theta = pi/2Q, L = cot 4 theta. The return map to L1 sends a state (y, b) along
s(b, b'; y) = [L (sin b theta - sin b' theta) - y cos b theta] / cos b' theta. The prong R_4 (B.2 (c)) has heights
    y1 = -s(-4, -2; 0),   y2 = -s(2, 4; y1),     y3 = s(-4, -10; y2),
    y4 = s(-10, -12; y3), y5 = -s(-12, -10; y4), y6 = s(10, 4; y5),
and folds provided
    (I1)-(I12)  0 < y_i < 1 (i = 1..6)          (I13)  y6 cos 4 theta < theta_0 := L (3 sin 4 theta - 4 sin 2 theta);
the prong R_{Q-1} (B.2 (d)) folds at even m provided
    (I14)  2 cos 4 theta sin 3 theta + L [cos 7 theta - sin 4 theta - 4 sin theta (cos 4 theta - sin 7 theta)/sin 4 theta] < theta_0.
Each is the positivity of a quotient num/den, den a product of sin 4 theta and of cos k theta (k <= 12), all positive for
k theta < pi/2, i.e. k < Q -- true at m >= 3 (Q >= 15). num is a trigonometric polynomial sum_j a_j e^{i j theta} with
Gaussian-rational a_j, built here exactly (stdlib fractions; no floating point in the Taylor tier).

Statements checked, as B.2 (e) and B.4 state them:
  * the fourteen are built from the B.2 formulas above, and every denominator is a product of sin 4 theta and of
    cos k theta with k <= 12;
  * INTERVAL TIER: at every m with 3 <= m <= 3000 every one of the fourteen is positive in outward-rounded interval
    arithmetic (mpmath.iv, pi enclosed; PASS iff the lower end of every interval is > 0); (I14) at even m, as B.2 (d)
    needs it -- and, beyond III, at odd m too;
  * TAYLOR TIER: for 0 < theta <= 1/7600 each numerator is positive, from its exact Taylor coefficients at theta = 0,
    c_n = sum_j a_j (i j)^n / n!, and the Lagrange remainder sum_j |a_j| |j|^N theta^N / N!: with n0 the first nonzero
    order, num/theta^n0 >= c_n0 - sum_{n0<n<N} |c_n| theta0^(n-n0) - B_N theta0^(N-n0) > 0;
  * theta <= 1/7600 contains every m >= 2984 (rigorous: pi < 3.1415927), so the two tiers overlap;
  * each numerator vanishes at theta = 0 to order 1, 2 or 3, and the least relative margin of the Taylor bound is 40%
    (the fourteen numerators, their orders N and margins are written to sB_4_record.json -- the deposited certificate);
  * the closed form is the billiard: at m = 3..8 a rigorous interval trace of R_4 departs L1 along +2, -4, -10, -12,
    +10, +4 at heights that agree with y1..y6, then reaches L1 head-on;
  * m = 1, 2 (Q = 7, 11): a rigorous interval-arithmetic trace of every R-prong R_2, ..., R_{Q-1} reaches L1 head-on,
    every side choice interval-separated and every hit interval-interior to its side (L1 hits 3, 4, 1 at Q = 7 and
    3, 8, 3, 2, 5 at Q = 11, the hit at the head-on arrival included);
  * CONTROLS (must fire): (a) the Taylor test rejects sin^n theta (1 - 7601 sin theta), n = 1, 2, 3 -- leading
    coefficient positive, negative at theta = 1/7600 -- and accepts the same with 7599, positive on (0, 1/7600]: the
    remainder bound is neither vacuous nor slack at the edge of its range;
    (b) at m = 1, 2 the conditions are outside their range (B.2 needs m >= 3) and the interval tier reports some of
    them non-positive; (c) the same tracer does NOT certify the corner prongs O2, A2 of the class 0 at 8/15, which is
    not completely periodic (III Section 6.3: they are the two ends of one O-A connection).

About a minute.
"""
import json
import math
import os
from fractions import Fraction as F

from mpmath import iv, mp

from common import HERE, check, done

THETA0 = F(1, 7600)
M_MAX = 3000


# ---------- exact arithmetic: Gaussian rationals, Laurent polynomials in z = e^{i theta}, quotients ----------

class G:
    __slots__ = ('re', 'im')

    def __init__(self, re, im=F(0)):
        self.re, self.im = F(re), F(im)

    def __add__(s, o): return G(s.re + o.re, s.im + o.im)
    def __mul__(s, o): return G(s.re * o.re - s.im * o.im, s.re * o.im + s.im * o.re)
    def neg(s): return G(-s.re, -s.im)
    def zero(s): return s.re == 0 and s.im == 0
    def absub(s): return abs(s.re) + abs(s.im)          # an upper bound for |a|


class Lp:
    def __init__(self, d=None):
        self.d = {k: v for k, v in (d or {}).items() if not v.zero()}

    def __add__(s, o):
        d = dict(s.d)
        for k, v in o.d.items():
            d[k] = d[k] + v if k in d else v
        return Lp(d)

    def __neg__(s): return Lp({k: v.neg() for k, v in s.d.items()})
    def __sub__(s, o): return s + (-o)

    def __mul__(s, o):
        d = {}
        for a, x in s.d.items():
            for b, y in o.d.items():
                d[a + b] = d[a + b] + x * y if a + b in d else x * y
        return Lp(d)


def lS(k):   # sin k theta = (z^k - z^-k) / 2i
    return Lp({k: G(0, F(-1, 2)), -k: G(0, F(1, 2))}) if k else Lp()


def lC(k):
    return Lp({k: G(F(1, 2)), -k: G(F(1, 2))}) if k else Lp({0: G(1)})


class Qt:
    """num / den, den a tuple of positive primitives ('S', k) / ('C', k)."""
    prim = None

    def __init__(self, num, den=()):
        self.num, self.den = num, tuple(sorted(den))

    @staticmethod
    def lp(p): return lS(p[1]) if p[0] == 'S' else lC(p[1])

    def lift(s, den):
        rest = list(den)
        for p in s.den:
            rest.remove(p)
        n = s.num
        for p in rest:
            n = n * Qt.lp(p)
        return n

    def common(s, o):
        a, out = list(s.den), list(s.den)
        for p in o.den:
            if p in a:
                a.remove(p)
            else:
                out.append(p)
        return tuple(sorted(out))

    def __add__(s, o):
        o = o if isinstance(o, Qt) else const(o)
        d = s.common(o)
        return Qt(s.lift(d) + o.lift(d), d)

    __radd__ = __add__
    def __neg__(s): return Qt(-s.num, s.den)
    def __sub__(s, o): return s + (-(o if isinstance(o, Qt) else const(o)))
    def __rsub__(s, o): return const(o) - s

    def __mul__(s, o):
        o = o if isinstance(o, Qt) else const(o)
        return Qt(s.num * o.num, s.den + o.den)

    __rmul__ = __mul__

    def __truediv__(s, o):
        assert isinstance(o, Qt) and o.den == () and o.prim is not None, 'division by a primitive only'
        return Qt(s.num, s.den + (o.prim,))


def const(c): return Qt(Lp({0: G(F(c))}))


def qS(k):
    q = Qt(lS(k)) if k >= 0 else -Qt(lS(-k))
    q.prim = ('S', k) if k > 0 else None
    return q


def qC(k):
    q = Qt(lC(abs(k)))
    q.prim = ('C', abs(k))
    return q


# ---------- the fourteen conditions, generic in the backend ----------

def conditions(S, C, div):
    """[(name, value)] for (I1)-(I14) as B.2 writes them. S, C: sin k theta, cos k theta; div(x, prim) = x / prim."""
    L = lambda x: div(x * C(4), S(4))
    s = lambda b, bp, y: div(L(S(b) - S(bp)) - y * C(b), C(bp))
    th0 = L(3 * S(4) - 4 * S(2))
    y1 = -s(-4, -2, 0 * S(1))
    y2 = -s(2, 4, y1)
    y3 = s(-4, -10, y2)
    y4 = s(-10, -12, y3)
    y5 = -s(-12, -10, y4)
    y6 = s(10, 4, y5)
    out = []
    for i, y in enumerate((y1, y2, y3, y4, y5, y6), 1):
        out += [(f'y{i} > 0', y), (f'y{i} < 1', 1 - y)]
    out.append(('y6 cos4t < theta0', th0 - y6 * C(4)))
    out.append(('(I14) R_{Q-1} fold', th0 - (2 * C(4) * S(3)
                                            + L(C(7) - S(4) - div(4 * S(1) * (C(4) - S(7)), S(4))))))
    return out


def heights(S, C, div):
    L = lambda x: div(x * C(4), S(4))
    s = lambda b, bp, y: div(L(S(b) - S(bp)) - y * C(b), C(bp))
    ys = [-s(-4, -2, 0 * S(1))]
    for (b, bp, sg) in ((2, 4, -1), (-4, -10, 1), (-10, -12, 1), (-12, -10, -1), (10, 4, 1)):
        ys.append(sg * s(b, bp, ys[-1]))
    return ys


def exact_conditions():
    return conditions(qS, qC, lambda x, p: x / p)


def iv_conditions(m):
    th = iv.pi / (2 * (4 * m + 3))
    return conditions(lambda k: iv.sin(k * th), lambda k: iv.cos(k * th), lambda x, p: x / p)


# ---------- the Taylor tier ----------

def taylor(num, max_extra=8):
    """(ok, n0, c_n0, N, margin) for a Laurent numerator: positive on (0, THETA0] by the Lagrange bound."""
    a = num.d

    def coef(n):
        tot = G(0)
        for j, v in a.items():
            p = G(1)
            for _ in range(n):
                p = p * G(0, j)
            tot = tot + v * p
        assert tot.im == 0, 'non-real Taylor coefficient'
        return tot.re / math.factorial(n)
    c = []
    n0 = None
    for n in range(40):
        c.append(coef(n))
        if c[-1] != 0:
            n0 = n
            break
    if n0 is None or c[n0] < 0:
        return False, n0, c[-1] if c else None, None, None
    for N in range(n0 + 1, n0 + 1 + max_extra):
        while len(c) < N:
            c.append(coef(len(c)))
        BN = sum(v.absub() * F(abs(j)) ** N for j, v in a.items()) / math.factorial(N)
        low = c[n0] - sum(abs(c[n]) * THETA0 ** (n - n0) for n in range(n0 + 1, N)) - BN * THETA0 ** (N - n0)
        if low > 0:
            return True, n0, c[n0], N, low / c[n0]
    return False, n0, c[n0], None, None


def sin_power(n):
    p = Lp({0: G(1)})
    for _ in range(n):
        p = p * lS(1)
    return p


# ---------- the rigorous interval tracer (from a corner of T) ----------

def trace(P, Q, start, u0, cap=20000):
    """('headon', [(departure b, height y)] at each L1 hit) or ('cap', ...); raises RuntimeError on any unresolved
    step (a side choice not interval-separated, or a hit not interval-interior to its side)."""
    e = iv.pi / (2 * Q)
    al = P * e
    L = iv.cos(al) / iv.sin(al)
    sa, ca = iv.sin(al), iv.cos(al)
    x, y, excl = {'O': (iv.mpf(0), iv.mpf(0), {'L2', 'H'}), 'A': (L, iv.mpf(1), {'L1', 'H'}),
                  'R': (L, iv.mpf(0), {'L1', 'L2'})}[start]
    u, hits = u0 % (4 * Q), []
    for k in range(cap):
        dx, dy = iv.cos(u * e), iv.sin(u * e)
        cands = []
        if 'L2' not in excl and u % (2 * Q):
            if dy.b < 0:
                cands.append((-y / dy, 'L2'))
            elif not dy.a > 0:
                raise RuntimeError('ambiguous dy')
        if 'L1' not in excl and (u - Q) % (2 * Q):
            if dx.a > 0:
                cands.append(((L - x) / dx, 'L1'))
            elif not dx.b < 0:
                raise RuntimeError('ambiguous dx')
        h, dh = x * sa - y * ca, dx * sa - dy * ca
        if 'H' not in excl and (u - P) % (2 * Q):
            if dh.b < 0:
                cands.append((-h / dh, 'H'))
            elif not dh.a > 0:
                raise RuntimeError('ambiguous dh')
        if not cands:
            raise RuntimeError('no forward side')
        cands.sort(key=lambda c: c[0].a)
        t, side = cands[0]
        if not t.a > 0 or (len(cands) > 1 and not t.b < cands[1][0].a):
            raise RuntimeError(f'ambiguous side at bounce {k}')
        x, y = x + t * dx, y + t * dy
        if side == 'L1':
            if not (y.a > 0 and y.b < 1):
                raise RuntimeError(f'L1 hit not interior at bounce {k}')
            if u % (4 * Q) == 0:
                hits.append((0, y))
                return 'headon', hits
            u = (2 * Q - u) % (4 * Q)
            hits.append((u - 2 * Q, y))                  # the departure b = u - 2Q in (-Q, Q)
        else:
            if not (x.a > 0 and x.b < L.a):
                raise RuntimeError(f'{side} hit not interior at bounce {k}')
            u = (-u) % (4 * Q) if side == 'L2' else (2 * P - u) % (4 * Q)
        excl = {side}
    return 'cap', hits


def main():
    iv.dps = 50
    mp.dps = 40

    # the fourteen, exactly
    ex = exact_conditions()
    check(f'(I1)-(I14) built from B.2: {len(ex)} conditions', len(ex) == 14)
    prims = sorted({p for _, q in ex for p in q.den})
    check(f'every denominator is a product of sin 4t and cos kt, k <= 12: primitives {prims}',
          all(p == ('S', 4) or (p[0] == 'C' and p[1] <= 12) for p in prims))
    # transcription: the exact numerator over its denominator equals the interval value (m = 3, 100)
    worst = 0
    for m in (3, 100):
        th = mp.pi / (2 * (4 * m + 3))
        vals = iv_conditions(m)
        for (nm, q), (_, v) in zip(ex, vals):
            num = sum(mp.mpc(float(c.re), float(c.im)) * mp.expj(j * th) for j, c in q.num.d.items())
            den = mp.mpf(1)
            for p in q.den:
                den *= mp.sin(p[1] * th) if p[0] == 'S' else mp.cos(p[1] * th)
            worst = max(worst, abs(num.real / den - mp.mpf(v.mid)) / max(abs(mp.mpf(v.mid)), mp.mpf(1e-30)))
    check(f'exact numerators / denominators = the interval values at m = 3, 100 (max rel. diff {float(worst):.1e})',
          worst < 1e-12)

    # interval tier
    lows = [None] * 14
    odd14 = None
    for m in range(3, M_MAX + 1):
        for i, (_, v) in enumerate(iv_conditions(m)):
            lo = float(v.a)
            if i == 13 and m % 2:
                odd14 = lo if odd14 is None else min(odd14, lo)
                continue
            if lows[i] is None or lo < lows[i][0]:
                lows[i] = (lo, m)
    check(f'interval tier, every m in [3, {M_MAX}]: all fourteen positive (least lower end {min(l[0] for l in lows):.3e} '
          f'at m = {min(lows)[1]}; (I14) at even m)', all(l[0] > 0 for l in lows))
    check(f'(beyond III) (I14) is positive at odd m in [3, {M_MAX}] too (least lower end {odd14:.3e})', odd14 > 0)

    # Taylor tier
    PI_UP = F(31415927, 10 ** 7)
    m_lo = next(m for m in range(1, 10 ** 5) if PI_UP / (2 * (4 * m + 3)) <= THETA0)
    check(f'theta <= 1/7600 for every m >= {m_lo} (pi < 3.1415927): the tiers overlap on [{m_lo}, {M_MAX}]',
          m_lo == 2984 and m_lo <= M_MAX)
    rec, ok_all = [], True
    for nm, q in ex:
        ok, n0, cn0, N, marg = taylor(q.num)
        ok_all &= ok
        rec.append({'condition': nm, 'ok': ok, 'order': n0, 'leading_coefficient': str(cn0), 'N': N,
                    'relative_margin': str(marg), 'relative_margin_float': float(marg) if marg else None,
                    'denominator': [f'{p[0]}{p[1]}' for p in q.den],
                    'numerator': {str(j): [str(c.re), str(c.im)] for j, c in sorted(q.num.d.items())}})
    check('Taylor tier, 0 < theta <= 1/7600: all fourteen numerators positive', ok_all)
    orders = sorted({r['order'] for r in rec})
    check(f'orders of vanishing at theta = 0: {orders} (III: 1, 2 or 3)', orders == [1, 2, 3])
    mins = min(rec, key=lambda r: r['relative_margin_float'])
    check(f'least relative margin {mins["relative_margin_float"]:.4f} at {mins["condition"]} (III: 40%)',
          math.floor(100 * mins['relative_margin_float']) == 40)
    with open(os.path.join(HERE, 'sB_4_record.json'), 'w') as f:
        json.dump({'theta0': '1/7600', 'interval_range': [3, M_MAX], 'conditions': rec}, f, indent=1)
    print('  wrote sB_4_record.json (the fourteen numerators, N and the margins)')

    # the closed form is the billiard: R_4 at m = 3..8 (a trace widens about one digit per bounce: dps 150)
    iv.dps = 150
    agree = True
    for m in range(3, 9):
        Q = 4 * m + 3
        st, hits = trace(4, Q, 'R', 2 * Q - 4)
        th = iv.pi / (2 * Q)
        ys = heights(lambda k: iv.sin(k * th), lambda k: iv.cos(k * th), lambda x, p: x / p)
        deps = [b for b, _ in hits[:6]]
        close = all(abs(float(hy.mid) - float(y.mid)) < 1e-25 for (_, hy), y in zip(hits, ys))
        agree &= st == 'headon' and deps == [2, -4, -10, -12, 10, 4] and close
    check('R_4 at m = 3..8 (rigorous trace): departs +2, -4, -10, -12, +10, +4 at heights y1..y6, then head-on', agree)

    # m = 1, 2: every R-prong, rigorous trace
    want = {7: [3, 4, 1], 11: [3, 8, 3, 2, 5]}
    for Q in (7, 11):
        got = []
        for k in range(2, Q, 2):
            st, hits = trace(4, Q, 'R', 2 * Q - k)
            got.append(len(hits) if st == 'headon' else None)
        check(f'4/{Q}: every R-prong R_2..R_{Q - 1} reaches L1 head-on, L1 hits {got}', got == want[Q])

    # CONTROLS
    # (a) sin^n t (1 - k sin t), n = 1, 2, 3: leading coefficient 1 > 0; negative at t = 1/7600 when k = 7601 (it
    #     must be rejected), positive on (0, 1/7600] when k = 7599 (accepting it shows the bound is sharp there)
    rej = [not taylor(sin_power(n) - Lp({0: G(7601)}) * sin_power(n + 1))[0] for n in (1, 2, 3)]
    acc = [taylor(sin_power(n) - Lp({0: G(7599)}) * sin_power(n + 1))[0] for n in (1, 2, 3)]
    check(f'CONTROL (a): sin^n t (1 - k sin t), n = 1..3: rejected at k = 7601 {rej}, accepted at k = 7599 {acc}',
          all(rej) and all(acc))
    neg = []
    for m in (1, 2):
        bad = [nm for nm, v in iv_conditions(m) if not v.a > 0]
        neg.append(len(bad))
    check(f'CONTROL (b): at m = 1, 2 (outside B.2\'s range m >= 3) the interval tier reports {neg} conditions '
          f'not positive', all(n > 0 for n in neg))
    iv.dps = 150
    res = {}
    for st, u0, lab in (('O', 2, 'O2'), ('A', 2 * 15 + 8 + 2, 'A2')):
        try:
            res[lab] = trace(8, 15, st, u0)[0]
        except RuntimeError as ex_:
            res[lab] = 'unresolved'
    check(f'CONTROL (c): 8/15 class 0, the tracer does not certify O2, A2 head-on: {res}',
          all(v != 'headon' for v in res.values()))
    done()


if __name__ == '__main__':
    main()
