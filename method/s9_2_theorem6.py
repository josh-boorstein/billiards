"""III Section 9.2, Theorem 6: infinitely many counterexamples to Boshernitzan's Conjecture 2.2(a) on the line
(x, 2x, 1 - 3x) pi -- the identities and the interval arithmetic its proof rests on, on each of the three intervals
    I_1 = [3749379/23000000, 750069/4600000] (contains 15/92),  I_2 = [2099637/24200000, 10503509/121000000] (21/242),
    I_3 = [1437689/16875000, 23005103/270000000] (inside 23/270 + (1.1, 1.9) 10^-5).

The objects, as Section 9.2 defines them. The triangle Delta_x is placed with V0 = 0, V1 = 2 sin 3A, V2 = 2 sin 2A e^{iA},
A = pi x; sides s0 = V0V1, s1 = V1V2, s2 = V2V0. Every side direction is a multiple of A mod pi, so a direction class is
named by a FORMAL label (j, m mod 2), theta = j A + m pi, the same object for every x (reflection in s0: j -> -j, in s2:
j -> 2 - j, in s1: j -> -4 - j; m -> -m each time). A word is a sequence of (side, formal direction) strips; along a fixed
word every hit parameter is affine in the start parameter u, with coefficients Laurent polynomials in z = e^{i pi x}
over Z[i] divided by one common D0 -- exact (class LP below). For each interval: a transversal J (a sub-interval
`win` of one strip), and the words whose first returns tile J, PROPOSED by a 50-digit tracer at one point x0 of I and
then CERTIFIED on all of I: every strip constraint of every word, evaluated at both ends of the word's piece (each
constraint is affine in u, so its two ends decide it), is either identically 0 or positive on I by interval arithmetic
(mpmath.iv, I enclosed outward, bisection).

Statements checked, per interval:
  * the side-parallel class: j odd on I_1 and I_3, j even on I_2;
  * four words tile J, the tiling joins are exact identities, every constraint is identically 0 or certified positive on
    I, the start direction is inward on I, and every word returns with slope +1 (exactly);
  * the four words merge into three blocks of equal translation, of lengths l1, l2, l3, with translations l2 + l3,
    l3 - l1, -(l1 + l2) identically, and l1, l2, l3 > 0 on I: the first return is a symmetric exchange of three
    intervals, rotation number rho = (l2 + l3)/(|J| + l2) = Num/Den;
  * after dividing by i, Num and Den are integer Laurent polynomials in w = z^2; their gcd is a product of cyclotomic
    polynomials of orders 1, 2, 3, 6 (zero only at x in {0, 1/6, 1/3, 1/2}, none in I); Den != 0 on I; the reduced pair
    is coprime, the reduced Den is not a monomial, and their joint exponent span is at most 20;
  * phi(d) <= 20 only for d <= 66, and no rational of I has d <= 66 -- the least denominators are 92, 242, 223;
  * a periodic orbit in the same direction: a rational start on a side, return translation identically 0, slope +1,
    every crossing certified inside on I, start direction inward on I -- 282, 118 and 126 bounces;
  * 149/914 and 179/1098 lie in I_1 with even denominator (III's examples); 15/92 in I_1, 21/242 in I_2, and
    I_3 inside 23/270 + (1.1, 1.9) 10^-5.
  CONTROLS (must fire): (a) on [x0 - 3.5e-5, x0 - 3e-5], past I_1's chamber (which ends near x0 - 2.8e-5), the same
  words fail constraints with a genuine negative interval; (b) the periodic-orbit test does not certify a start on J
  (u = 3/10, inside the minimal component). (A start next to u* is no control: it lies in the same cylinder and is
  periodic too.)
  Cross-check (beyond III): at 15/92 the w-span before the gcd is 26 (the independent read's figure).

About five minutes.
"""
import math
from fractions import Fraction as F

import mpmath as mp
import sympy
from mpmath import iv

from common import check, done

mp.mp.dps = 50
iv.dps = 40
SIDES = [(0, 1), (1, 2), (2, 0)]

INTERVALS = [
    dict(name='I_1', x0=F(15, 92), strip=(2, 7, 1), win=(F(1, 5), F(9, 20)),
         I=(F(3749379, 23000000), F(750069, 4600000)), words=[220, 404, 404, 184], groups=[[0], [1, 2], [3]],
         periodic=(0, F(69, 200), 7, 1, 282), least_d=92),
    dict(name='I_2', x0=F(21, 242), strip=(2, 0, 0), win=(F(4, 25), F(7, 25)),
         I=(F(2099637, 24200000), F(10503509, 121000000)), words=[178, 290, 290, 112], groups=[[0], [1, 2], [3]],
         periodic=(2, F(11, 100), 0, 0, 118), least_d=242),
    dict(name='I_3', x0=F(460081, 5400000), strip=(2, -11, 1), win=(F(7, 50), F(8, 25)),
         I=(F(1437689, 16875000), F(23005103, 270000000)), words=[314, 1060, 746, 746], groups=[[0], [1], [2, 3]],
         periodic=(0, F(1, 400), -11, 1, 126), least_d=223),
]


def mpf(q): return mp.mpf(q.numerator) / q.denominator


def refl(s, j, m):
    return ((-j, (-m) % 2) if s == 0 else (2 - j, (-m) % 2) if s == 2 else (-4 - j, (-m) % 2))


# ---------- the 50-digit tracer (proposes words; nothing below trusts it) ----------

def mtri(x):
    A = mp.pi * mpf(x)
    return [mp.mpc(0), mp.mpc(2 * mp.sin(3 * A), 0), 2 * mp.sin(2 * A) * mp.expj(A)], A


def mcross(a, b): return mp.im(mp.conj(a) * b)


def mtrace_word(V, A, s, lam, j, m, strip, win, cap):
    """Free trace from (side s, parameter lam, direction (j, m)) to the first arrival on `strip` inside win."""
    i, k = SIDES[s]
    p = V[i] + lam * (V[k] - V[i])
    word = []
    for _ in range(cap):
        dv = mp.expj(j * A + m * mp.pi)
        best = None
        for s2 in range(3):
            if s2 == s:
                continue
            i2, k2 = SIDES[s2]
            e = V[k2] - V[i2]
            den = mcross(dv, e)
            if den == 0:
                continue
            l2 = mcross(V[i2] - p, dv) / den
            if 0 <= l2 <= 1:
                t = mcross(V[i2] - p, e) / den
                if t > 0 and (best is None or t < best[0]):
                    best = (t, s2, l2)
        if best is None:
            return None
        _, s, lam = best
        i, k = SIDES[s]
        p = V[i] + lam * (V[k] - V[i])
        j, m = refl(s, j, m)
        word.append((s, j, m))
        if (s, j, m) == strip and win[0] <= lam <= win[1]:
            return word
    return None


def affine_word(V, A, strip, word):
    """Per step, the hit parameter a + b u along the fixed word."""
    out = []
    for u in (mp.mpf(0), mp.mpf(1)):
        s, j, m = strip
        i, k = SIDES[s]
        p = V[i] + u * (V[k] - V[i])
        lams = []
        for (s2, j2, m2) in word:
            dv = mp.expj(j * A + m * mp.pi)
            i2, k2 = SIDES[s2]
            e = V[k2] - V[i2]
            l2 = mcross(V[i2] - p, dv) / mcross(dv, e)
            p = V[i2] + l2 * e
            lams.append(l2)
            s, j, m = s2, j2, m2
        out.append(lams)
    return [(a, b - a) for a, b in zip(*out)]


def tile(x0, strip, win, cap=3000):
    V, A = mtri(x0)
    a0, b0 = mpf(win[0]), mpf(win[1])
    u, words = a0, []
    while b0 - u > mp.mpf(10) ** -40:
        us = u + mp.mpf(10) ** -30
        word = mtrace_word(V, A, strip[0], us, strip[1], strip[2], strip, (a0, b0), cap)
        assert word is not None, 'tracer cap'
        aff = affine_word(V, A, strip, word)
        lo, hi, pattern = mp.mpf(-10), mp.mpf(10), []
        for w_idx, ((a, b), st) in enumerate(zip(aff, word)):
            cons = [(0, 1)]
            if st == strip and w_idx < len(word) - 1:      # an earlier arrival must lie outside win
                val = a + b * us
                cons.append((-10, a0) if val < a0 else (b0, 10))
                pattern.append(0 if val < a0 else 1)
            if w_idx == len(word) - 1:
                cons.append((a0, b0))
            for lower, upper in cons:
                if b == 0:
                    continue
                u1, u2 = (lower - a) / b, (upper - a) / b
                if b < 0:
                    u1, u2 = u2, u1
                lo, hi = max(lo, u1), min(hi, u2)
        lo, hi = max(lo, a0), min(hi, b0)
        words.append({'word': tuple(st[0] for st in word), 'pattern': tuple(pattern), 'N': len(word),
                      'lo': lo, 'hi': hi, 'tr': aff[-1][0]})
        assert hi > u, 'tiling stuck'
        u = hi
    return words


# ---------- exact: Laurent polynomials in z = e^{i pi x} over Z[i] ----------

class LP:
    __slots__ = ('c',)

    def __init__(self, c=None):
        self.c = {e: v for e, v in (c or {}).items() if v != (0, 0)}

    @staticmethod
    def mono(e, re=1, im=0): return LP({e: (re, im)})

    def __add__(self, o):
        c = dict(self.c)
        for e, (a, b) in o.c.items():
            x, y = c.get(e, (0, 0))
            c[e] = (x + a, y + b)
        return LP(c)

    def __neg__(self): return LP({e: (-a, -b) for e, (a, b) in self.c.items()})
    def __sub__(self, o): return self + (-o)

    def __mul__(self, o):
        if isinstance(o, int):
            return LP({e: (a * o, b * o) for e, (a, b) in self.c.items()})
        c = {}
        for e1, (a1, b1) in self.c.items():
            for e2, (a2, b2) in o.c.items():
                x, y = c.get(e1 + e2, (0, 0))
                c[e1 + e2] = (x + a1 * a2 - b1 * b2, y + a1 * b2 + b1 * a2)
        return LP(c)

    def conj(self): return LP({-e: (a, -b) for e, (a, b) in self.c.items()})
    def zero(self): return not self.c
    def span(self): return (min(self.c), max(self.c)) if self.c else (0, 0)
    def ev(self, z): return sum(mp.mpc(a, b) * z ** e for e, (a, b) in self.c.items())

    def civ(self, X):
        """(Re, Im) at z = e^{i pi x}, x in the interval X."""
        re_, im_ = iv.mpf(0), iv.mpf(0)
        for e, (a, b) in self.c.items():
            t = iv.pi * e * X
            cs, sn = iv.cos(t), iv.sin(t)
            re_ += a * cs - b * sn
            im_ += a * sn + b * cs
        return re_, im_


def lp_cross(a, b): return a.conj() * b - a * b.conj()        # = 2i cross(a, b)


def lp_tri():
    s = lambda k: (LP.mono(k) - LP.mono(-k)) * LP.mono(0, 0, -1)   # 2 sin kA = -i (z^k - z^-k)
    return [LP(), s(3), s(2) * LP.mono(1)]


SIDE_ROT = {0: LP.mono(0), 2: LP.mono(2), 1: LP.mono(-4)}           # e^{2 i sigma_s}


def lp_unfold(word):
    """The copies of Delta_x along the word, each as its three vertices (exact)."""
    base = lp_tri()
    R, eps, t = LP.mono(0), 0, LP()
    copies = []
    for s in word:
        P = base[SIDES[s][0]]
        rot = SIDE_ROT[s]
        shift = P - rot * P.conj()
        g = lambda w, R=R, eps=eps, t=t: R * (w.conj() if eps else w) + t
        R, eps, t = (R * rot, 1, g(shift)) if eps == 0 else (R * rot.conj(), 0, g(shift))
        copies.append((s, [R * (v.conj() if eps else v) + t for v in base]))
    return copies


def enclose(I): return mpf(I[0]) - mp.mpf(10) ** -40, mpf(I[1]) + mp.mpf(10) ** -40


def positive(num, den, lo, hi, sg=1, maxdep=18):
    """sg Re(num conj(den)) > 0 on [lo, hi] by bisection: (True, depth) / (False, 'negative' | 'depth')."""
    stack, depth = [(lo, hi, 0)], 0
    while stack:
        a, b, dep = stack.pop()
        X = iv.mpf([a, b])
        nr, ni = num.civ(X)
        dr, di = den.civ(X)
        v = sg * (nr * dr + ni * di)
        if v.a > 0:
            depth = max(depth, dep)
            continue
        if v.b < 0:
            return False, 'negative'
        if dep > maxdep:
            return False, 'depth'
        stack += [(a, (a + b) / 2, dep + 1), ((a + b) / 2, b, dep + 1)]
    return True, depth


def im_positive(p, lo, hi):
    stack = [(lo, hi, 0)]
    while stack:
        a, b, dep = stack.pop()
        _, im_ = p.civ(iv.mpf([a, b]))
        if im_.a > 0:
            continue
        if im_.b < 0 or dep > 18:
            return False
        stack += [(a, (a + b) / 2, dep + 1), ((a + b) / 2, b, dep + 1)]
    return True


# ---------- the certificate ----------

def certificate(c, I=None):
    I = I or c['I']
    x0, strip, win = c['x0'], c['strip'], c['win']
    words = tile(x0, strip, win)
    z0 = mp.expj(mp.pi * mpf(x0))
    base = lp_tri()
    s0, j0, m0 = strip
    i0, k0 = SIDES[s0]
    dv = LP.mono(j0) * (-1 if m0 else 1)
    D0 = lp_cross(base[k0] - base[i0], dv)
    W = []
    for w in words:
        prev, steps, jj, mm = base, [], j0, m0
        for (s, cp) in lp_unfold(w['word']):
            i, k = SIDES[s]
            jj, mm = refl(s, jj, mm)
            steps.append((lp_cross(prev[i] - base[i0], dv), lp_cross(dv, prev[k] - prev[i]), (s, jj, mm) == strip))
            prev = cp
        W.append(steps)

    def endpoint(bp):
        """the constraint binding at the numerical endpoint bp: u = (C q - E p)/(D0 q), c = p/q."""
        for steps in W:
            for (C, E, _) in steps:
                for cc in (F(0), F(1), win[0], win[1]):
                    num = C * cc.denominator - E * cc.numerator
                    if abs(num.ev(z0) / (D0.ev(z0) * cc.denominator) - bp) < mp.mpf(10) ** -35:
                        return (C, E, cc)
        raise AssertionError('no endpoint formula')
    wf = (mpf(win[0]), mpf(win[1]))
    ends = [(('const', win[0]) if abs(w['lo'] - wf[0]) < mp.mpf(10) ** -40 else ('f', endpoint(w['lo'])),
             ('const', win[1]) if abs(w['hi'] - wf[1]) < mp.mpf(10) ** -40 else ('f', endpoint(w['hi'])))
            for w in words]

    def U(e):                                          # u = U_num / U_den, exactly
        if e[0] == 'const':
            return D0 * e[1].numerator, D0 * e[1].denominator
        C, E, cc = e[1]
        return C * cc.denominator - E * cc.numerator, D0 * cc.denominator
    tiling = all(((U(ends[i][1])[0] * U(ends[i + 1][0])[1]) - (U(ends[i + 1][0])[0] * U(ends[i][1])[1])).zero()
                 for i in range(len(words) - 1))
    checks = []
    for wi, steps in enumerate(W):
        pattern, pi_ = list(words[wi]['pattern']), 0
        for kst, (C, E, arr) in enumerate(steps):
            last = kst == len(steps) - 1
            cons = [(F(0), +1), (F(1), -1)]
            if arr and not last:
                cons.append((win[0], -1) if pattern[pi_] == 0 else (win[1], +1))
                pi_ += 1
            if last:
                cons += [(win[0], +1), (win[1], -1)]
            for e in ends[wi]:
                un, ud = U(e)
                for cc, sg in cons:
                    checks.append(((C * ud - un * D0) * cc.denominator - E * ud * cc.numerator,
                                   E * ud * cc.denominator, sg))
    lo, hi = enclose(I)
    nzero, fails, depth = 0, [], 0
    for num, den, sg in checks:
        if num.zero():
            nzero += 1
            continue
        ok, info = positive(num, den, lo, hi, sg)
        if ok:
            depth = max(depth, info)
        else:
            fails.append(info)
    slope1 = all((D0 + steps[-1][1]).zero() for steps in W)
    inward = im_positive(D0, lo, hi)
    return dict(words=words, W=W, ends=ends, U=U, D0=D0, tiling=tiling, n_checks=len(checks), n_zero=nzero,
                fails=fails, depth=depth, slope1=slope1, inward=inward)


def sub(a, b): return (a[0] * b[1] - b[0] * a[1], a[1] * b[1])
def add(a, b): return (a[0] * b[1] + b[0] * a[1], a[1] * b[1])
def const(q): return (LP.mono(0) * q.numerator, LP.mono(0) * q.denominator)


def structure(c, R):
    """The three blocks, the symmetric 3-IET identities, l1, l2, l3 > 0 on I, and rho = Num/Den."""
    trs = [w['tr'] for w in R['words']]
    groups = [[0]]
    for i in range(1, len(trs)):
        (groups[-1].append(i) if abs(trs[i] - trs[groups[-1][-1]]) < mp.mpf(10) ** -35 else groups.append([i]))
    U, ends, win = R['U'], R['ends'], c['win']
    tau = [(steps[-1][0], steps[-1][1]) for steps in R['W']]      # final lam = C/E + u
    b1, b2 = U(ends[groups[0][-1]][1]), U(ends[groups[1][-1]][1])
    l1, l2, l3 = sub(b1, const(win[0])), sub(b2, b1), sub(const(win[1]), b2)
    z = lambda r: r[0].zero()
    ids = (all(z(sub(tau[g[0]], tau[k])) for g in groups for k in g)
           and z(sub(tau[groups[0][0]], add(l2, l3))) and z(sub(tau[groups[1][0]], sub(l3, l1)))
           and z(add(tau[groups[2][0]], add(l1, l2))))
    lo, hi = enclose(c['I'])
    lens_pos = all(positive(l[0], l[1], lo, hi)[0] for l in (l1, l2, l3))
    # rho = (l2 + l3)/(T + l2), T = |J|; over the common denominator D0 (the s630 form)
    (n1, d1), (n2, d2) = b1, b2

    def q_of(d):
        e = next(iter(R['D0'].c))
        a0, bb0 = R['D0'].c[e]
        a, bb = d.c[e]
        q = F(a, a0) if a0 else F(bb, bb0)
        assert q.denominator == 1 and (R['D0'] * int(q)).c == d.c
        return int(q)
    q1, q2 = q_of(d1), q_of(d2)
    w0, w1 = win
    M = w0.denominator * w1.denominator * q1 * q2
    Num = R['D0'] * int(w1 * M) - n1 * (M // q1)
    Den = R['D0'] * int((w1 - w0) * M) + n2 * (M // q2) - n1 * (M // q1)
    return groups, ids, lens_pos, Num, Den


def galois(Num, Den):
    """Reduce rho = Num/Den to coprime integer polynomials in w = z^2; return the facts Section 9.2 uses."""
    parity = {e % 2 for e in list(Num.c) + list(Den.c)}
    kinds = {('re' if b == 0 else 'im' if a == 0 else 'mixed') for p in (Num, Den) for (a, b) in p.c.values()}
    assert len(parity) == 1 and kinds == {'im'}, (parity, kinds)
    off = parity.pop()
    t = sympy.Symbol('t')
    lo = (min(Num.span()[0], Den.span()[0]) - off) // 2
    poly = lambda p: sympy.Poly(sum(sympy.Integer(b) * t ** ((e - off) // 2 - lo) for e, (a, b) in p.c.items()), t)
    Pn, Pd = poly(Num), poly(Den)
    G = sympy.gcd(Pn, Pd)
    span_before = max(Pn.degree(), Pd.degree())
    Rn, Rd = sympy.div(Pn, G)[0], sympy.div(Pd, G)[0]
    orders = []
    for f, _ in sympy.factor_list(G.as_expr())[1]:
        n = next((n for n in range(1, 13) if sympy.Poly(f, t) == sympy.Poly(sympy.cyclotomic_poly(n, t), t)
                  or sympy.Poly(-f, t) == sympy.Poly(sympy.cyclotomic_poly(n, t), t)), None)
        orders.append(n)

    def low(P):
        cs = P.all_coeffs()[::-1]
        return next(i for i, v in enumerate(cs) if v != 0)
    span = max(Rn.degree(), Rd.degree()) - min(low(Rn), low(Rd))
    coprime = sympy.gcd(Rn, Rd).degree() == 0
    den_monomial = len([v for v in Rd.all_coeffs() if v != 0]) == 1
    return dict(span_before=span_before, orders=orders, span=span, coprime=coprime, den_monomial=den_monomial)


def least_denominator(I):
    d = 1
    while True:
        a = math.ceil(I[0] * d)
        if F(a, d) <= I[1]:
            return d, F(a, d)
        d += 1


def periodic(c, side, u, j0, m0, nbounce=None):
    """The orbit from side `side` at parameter u in direction (j0, m0): its word to the first return to the same strip
    at the same parameter (50-digit tracer), then exactly: slope +1, translation identically 0, every crossing of the
    word inside its side on all of I, start direction inward on I."""
    x0, I = c['x0'], c['I']
    V, A = mtri(x0)
    strip = (side, j0, m0)
    word = mtrace_word(V, A, side, mpf(u), j0, m0, strip, (mpf(u) - mp.mpf(10) ** -30, mpf(u) + mp.mpf(10) ** -30),
                       4000)
    if word is None:
        return None
    base = lp_tri()
    i0, k0 = SIDES[side]
    dv = LP.mono(j0) * (-1 if m0 else 1)
    D0 = lp_cross(base[k0] - base[i0], dv)
    prev, cons = base, []
    for (sw, cp) in lp_unfold([w[0] for w in word]):
        ii, kk = SIDES[sw]
        cons.append((lp_cross(prev[ii] - base[i0], dv), lp_cross(dv, prev[kk] - prev[ii])))
        prev = cp
    C, E = cons[-1]
    slope1, trans0 = (D0 + E).zero(), C.zero()
    lo, hi = enclose(I)
    inside = all(positive((Cc * u.denominator - D0 * u.numerator) * cc.denominator - Ec * u.denominator * cc.numerator,
                          Ec * u.denominator * cc.denominator, lo, hi, sg)[0]
                 for (Cc, Ec) in cons[:-1] for cc, sg in ((F(0), 1), (F(1), -1)))
    return dict(N=len(word), slope1=slope1, trans0=trans0, inside=inside, inward=im_positive(D0, lo, hi))


def main():
    for c in INTERVALS:
        nm = c['name']
        print(f'--- {nm} = [{c["I"][0]}, {c["I"][1]}], strip (side, j, m) = {c["strip"]}, J = {c["win"]}')
        j = c['strip'][1]
        check(f'{nm}: side-parallel class, j = {j} ({"even" if j % 2 == 0 else "odd"})',
              (j % 2 == 0) == (nm == 'I_2'))
        R = certificate(c)
        Ns = [w['N'] for w in R['words']]
        check(f'{nm}: {len(Ns)} words tile J, bounces {Ns}; joins exact identities: {R["tiling"]}',
              Ns == c['words'] and R['tiling'])
        check(f'{nm}: {R["n_checks"]} constraints, {R["n_zero"]} identically 0, {len(R["fails"])} not certified '
              f'(bisection depth <= {R["depth"]})', not R['fails'])
        check(f'{nm}: every word returns with slope +1; start direction inward on I', R['slope1'] and R['inward'])
        groups, ids, lens_pos, Num, Den = structure(c, R)
        check(f'{nm}: blocks {groups}; translations l2 + l3, l3 - l1, -(l1 + l2) identically; l1, l2, l3 > 0 on I',
              groups == c['groups'] and ids and lens_pos)
        g = galois(Num, Den)
        if nm == 'I_1':
            check(f'(cross-check) 15/92: w-span before the gcd {g["span_before"]} (the independent read: 26)',
                  g['span_before'] == 26)
        check(f'{nm}: gcd(Num, Den) cyclotomic of orders {sorted(g["orders"])} (zeros only at x in '
              f'{{0, 1/6, 1/3, 1/2}}); reduced pair coprime, Den not a monomial; joint span {g["span"]} <= 20',
              all(o in (1, 2, 3, 6) for o in g['orders']) and g['coprime'] and not g['den_monomial']
              and g['span'] <= 20)
        dmin, rmin = least_denominator(c['I'])
        check(f'{nm}: least denominator in I is {dmin} ({rmin}) > 66', dmin == c['least_d'] and dmin > 66)
        side, u, j0, m0, nb = c['periodic']
        P = periodic(c, side, u, j0, m0)
        check(f'{nm}: periodic orbit from side {side} at u = {u}: {P["N"]} bounces, slope +1, translation identically 0, '
              f'crossings inside on I, inward on I',
              P['N'] == nb and P['slope1'] and P['trans0'] and P['inside'] and P['inward'])
        if nm == 'I_1':
            Rw = certificate(c, I=(c['x0'] - F(35, 10 ** 6), c['x0'] - F(30, 10 ** 6)))
            check(f'CONTROL (a): on [x0 - 3.5e-5, x0 - 3e-5], past the chamber, {len(Rw["fails"])} constraints fail, '
                  f'{Rw["fails"].count("negative")} with a negative interval', 'negative' in Rw['fails'])
            Pc = periodic(c, c['strip'][0], F(3, 10), c['strip'][1], c['strip'][2])
            check(f'CONTROL (b): a start on J (u = 3/10, in the minimal component) is not certified periodic: '
                  f'{"no return within 4000 bounces" if Pc is None else Pc}',
                  Pc is None or not (Pc['trans0'] and Pc['slope1'] and Pc['inside']))
    I1, I2, I3 = (c['I'] for c in INTERVALS)
    lo3, hi3 = I3[0] - F(23, 270), I3[1] - F(23, 270)
    check(f'15/92 in I_1, 21/242 in I_2; I_3 - 23/270 = [{float(lo3):.4e}, {float(hi3):.4e}] inside (1.1e-5, 1.9e-5)',
          I1[0] <= F(15, 92) <= I1[1] and I2[0] <= F(21, 242) <= I2[1]
          and F(11, 10 ** 6) < lo3 and hi3 < F(19, 10 ** 6))
    ex = [F(149, 914), F(179, 1098)]
    check(f'phi(d) <= 20 only for d <= {max(d for d in range(1, 1000) if sympy.totient(d) <= 20)}',
          max(d for d in range(1, 1000) if sympy.totient(d) <= 20) == 66)
    check('149/914 and 179/1098 lie in I_1, in lowest terms with d even',
          all(I1[0] <= q <= I1[1] and q.denominator % 2 == 0 and q in (F(149, 914), F(179, 1098)) for q in ex)
          and math.gcd(149, 914) == 1 and math.gcd(179, 1098) == 1)
    done()


if __name__ == '__main__':
    main()
