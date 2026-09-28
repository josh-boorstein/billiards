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
(intervalxt gates it on SAF != 0, where alone it can succeed; the SAF test cost 70% of the run at
`D = 208` and the certificate is verified exactly anyway, so here SAF is only reported): a periodic orbit through a component visits interval `j` `a_j >= 0` times with
`sum a_j t_j = 0` (`t_j` the translations, as Q-vectors in the power basis) -- so a Q-linear
functional `y` with `y . t_j > 0` for ALL `j` (Gordan) PROVES there is no periodic trajectory.
`y` is found by an LP and VERIFIED IN EXACT INTEGERS; the certificate is returned.
SAF = 0 components (intervalxt's self-similarity loop is NOT ported) end UNDETERMINED.

SIGNS.  Exact zero test (reduced coefficient vector); a nonzero sign from a float dot product with
a rigorous error bound, else a certified `mpmath.iv` enclosure (as `cyclofield.py`, which is
the cross-check reference; `_selftest`).

Entry: `decompose(abc, j, cap)` -> dict (components with kind / size / transverse measure /
certificate).  `perp_dirs(abc)`, `side_dirs(abc)` give `j` for the directions perpendicular /
parallel to each side.
"""
from fractions import Fraction
import math

import numpy as np
import sympy as sp
from mpmath import iv

BOSH_EVERY = 16


# ---------------------------------------------------------------------------------------------
# exact field Z[zeta_N], power basis, integer coefficient tuples
# ---------------------------------------------------------------------------------------------
class Cyc:
    def __init__(self, N):
        self.N = N
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
        self.nsign_iv = 0

    def z(self, k):
        return self.tab[k % self.N]

    def add(self, u, v):
        return tuple(a + b for a, b in zip(u, v))

    def sub(self, u, v):
        return tuple(a - b for a, b in zip(u, v))

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
    exactly."""
    F, d, N = T.F, T.d, T.N
    om = {e % N: c for e, c in j.items()} if isinstance(j, dict) else {j % N: 1}
    X = lambda p: T.red(T.cross2(om, p))
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
            c_om = F.sign(T.red(T.cross2(e, om)))
            if c_om == 0:
                continue                             # side parallel to the flow
            c_in = F.sign(T.red(T.cross2(e, madd(W[kk], mneg(W[i])))))
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
    __slots__ = ('top', 'bot', 'kind', 'cert', 'steps')

    def __init__(self, top, bot):
        self.top, self.bot = list(top), list(bot)
        self.kind, self.cert, self.steps = None, None, 0


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
        fa, fb = F.fval(a), F.fval(b)
        m = max(max(map(abs, a)), max(map(abs, b)))
        if abs(fb) > m * F.D * 1e-12:
            q = max(int(fa // fb) - 1, 0)
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
        """Shortest prefix of top and bot with equal label sets; split the rest off."""
        seen_t, seen_b = set(), set()
        n = len(c.top)
        for k in range(n):
            seen_t.add(c.top[k])
            seen_b.add(c.bot[k])
            if seen_t == seen_b:
                if k == n - 1:
                    return None
                rest = Comp(c.top[k + 1:], c.bot[k + 1:])
                c.top, c.bot = c.top[:k + 1], c.bot[:k + 1]
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
        lam = np.array([self.lam[l] for l in c.top], dtype=object)
        t = np.array([tr[l] for l in c.top], dtype=object)
        W = lam.T.dot(t)                            # sum_j lam_j (x) t_j
        A = W - W.T
        return all(x == 0 for x in A.flat)

    def boshernitzan(self, c):
        """A verified y in Z^D with y . t_j > 0 for all j, or None."""
        from scipy.optimize import linprog
        tr = self.translations(c)
        T = [tr[l] for l in c.top]
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
                      bounds=[(-1, 1)] * D + [(None, 1)], method='highs')
        if res.status != 0 or res.x[-1] <= 1e-12:
            return None
        y0 = res.x[:D]
        for scale in (10 ** 6, 10 ** 9, 10 ** 12, 10 ** 15):
            y = [int(round(v * scale)) for v in y0]
            if all(sum(a * b for a, b in zip(y, t)) > 0 for t in T):
                return y
        return None

    def run(self, top, bot, cap):
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
                # saddle connection at the left end: merge (intervalxt NON_SEPARATING)
                t0, b0 = c.top[0], c.bot[0]
                if t0 != b0 and F.cmp(self.lam[t0], self.lam[b0]) == 0:
                    it = c.bot.index(t0)
                    c.bot[it] = b0
                    self.h[b0] += self.h[t0]
                    del c.bot[0]
                    del c.top[0]
                    continue
                if c.steps >= cap:
                    c.kind = 'undetermined'
                    break
                # (intervalxt skips Boshernitzan when SAF = 0, where it cannot succeed; here the
                # certificate is verified exactly, so the costly SAF test is only reported, not gated)
                if c.steps % BOSH_EVERY == 0 and c.steps > 0:
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
    pieces, lengths, top, bot, meta = build_iet(T, j)
    Dc = Decomposer(T.F, lengths)
    comps = Dc.run(top, bot, cap)
    F = T.F
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
            'nsign_iv': F.nsign_iv}


class _HK:
    """(return time, crossings of each tracked entry key) -- composes like the int h."""
    __slots__ = ('h', 'k')

    def __init__(self, h, k):
        self.h, self.k = h, k

    def __add__(self, o):
        return _HK(self.h + o.h, tuple(a + b for a, b in zip(self.k, o.k)))

    def __rmul__(self, q):
        return _HK(q * self.h, tuple(q * a for a in self.k))


def crossings(abc, j, keys, cap=20000):
    """`decompose`, with every component's labels also carrying how many times their return orbit
    crosses each entry key `(g, s)` in `keys` (a copy crossing entered through side s of copy g).
    A cylinder's count at a key is its crossing multiplicity k_c there: the number of intervals it
    cuts on that edge.  Kac is asserted exactly, on the return times, as in `decompose`."""
    T = Triangle(abc)
    pieces, lengths, top, bot, meta = build_iet(T, j)
    idx = {key: i for i, key in enumerate(keys)}
    Dc = Decomposer(T.F, lengths)
    zero = (0,) * len(keys)
    Dc.h = {}
    for p in lengths:
        k = list(zero)
        if pieces[p][3] in idx:
            k[idx[pieces[p][3]]] = 1
        Dc.h[p] = _HK(1, tuple(k))
    comps = Dc.run(top, bot, cap)
    F = T.F
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
    return out, edge


def n_flow(P, Q, cap=20000, beam='1'):
    """n(P/Q), the number of branches of the perpendicular partition of L1 (the orphan paper's
    object; `partition.py`), from the flow decomposition: the right triangle
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
    comps, edge = crossings(abc, j, keys, cap)
    F = T.F
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
        E = F.mpval(edge[keys[i0]], 50)
        for c in cyl:
            wc = float(F.mpval(c['lam'], 50) / E)       # fval cancels on small elements
            for _ in range(c['k'][i0]):
                widths.append(wc)
                assert c['period'] % 2 == 0, 'a retracing orbit has even period'
                plen.append(c['period'] // 2)
        order = sorted(range(len(widths)), key=lambda t: (plen[t], widths[t]))
        widths, plen = [widths[t] for t in order], [plen[t] for t in order]
    return {'n': sums[0] if ok else None, 'ncyl': len(cyl), 'ncomp': len(comps), 'edges': sums,
            'tiled': tiled, 'all_cyl': len(cyl) == len(comps),
            'widths': widths, 'prefix_len': plen}


def _selftest():
    """Field against cyclofield.py on random elements."""
    import random
    import os
    import sys
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from cyclofield import CycloField
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
    print('selftest ok: mul / conj / sign agree with cyclofield on 200 elements')


if __name__ == '__main__':
    _selftest()
