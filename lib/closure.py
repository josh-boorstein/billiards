"""The exact closure certificate: is a perpendicular direction completely periodic?

THE OBJECT.  For a right triangle with centre P/Q, the unfolding S_alpha is the holonomy
double cover of a genus-zero quadratic differential B on the sphere, with cone points Z_O,
Z_A (the images of the acute vertices) and Q simple poles R_t (the images of the right
angle).  A direction on S_alpha is completely periodic iff it is on B, since a cylinder
contains no singularity (its preimage is one or two annuli) and the holonomy involution
permutes the cylinders upstairs.  On B a direction is completely periodic iff every
separatrix is a saddle connection, i.e. every leaf leaving a cone point or pole arrives
at one.  This module walks every separatrix and reports whether each arrives.

THE MODEL.  B is Q "kites" (each the rhombus of four triangle copies about R, modulo the
half-turn) glued in a ring along copies of H.  In a fixed direction the foliation crosses
each kite by one of two maps between its two interfaces:
  * the NEAR (shorter, in transverse measure) interface is mapped bijectively onto a
    sub-interval of the FAR one ("through band");
  * the rest of the far interface is folded onto itself by a half-turn ("fold band");
    the pole R_t is the fold's fixed point, and the fold's end is the sector corner,
    Z_O or Z_A.
With c[t] = (c0 + 2 P t) mod 2Q and ell[t] = sin(c[t] pi / 2Q) the transverse measure of
interface t, the fold lies on the larger of the kite's two interfaces, and which
sub-interval the through band maps onto is decided by whether the kite's angular sector
contains the direction (the O-sector rule) or not (the A-sector rule).  The two classes
of perpendicular direction at a centre are eps = 0 and eps = 1 (c0 = P, resp. P + Q at
odd Q and P + 1 at even Q).  At odd Q, eps = 0 is the direction perpendicular to L1 and
eps = 1 the one perpendicular to L2; at even Q, eps = 0 holds both legs'
perpendiculars and eps = 1 the hypotenuse's.

THE ARITHMETIC.  Every transverse coordinate is a Z-combination of the ell's with
half-integer coefficients, so twice a coordinate is exactly an element of Z[zeta_4Q]
(2i sin(c pi/2Q) = zeta^c - zeta^-c), carried as a reduced integer vector modulo
Phi_4Q.  An ARRIVAL (a leaf hitting a singularity) is decided by exact equality of these
vectors -- never by a float.  The ORDER tests (which band a point lies in) use a float
alongside the vector; whenever a float gap is within the float's own error bound
(times a safety factor) the sign is instead decided exactly, from the vector, at a
precision raised until it is proved (`certified_sign`).

THE VERDICT IS ONE-SIDED.  `capped == 0` (every separatrix arrived) PROVES complete
periodicity.  A separatrix that reaches the step cap returns NO verdict -- not a negative
one.  So this instrument can confirm complete periodicity and can never refute it.
"""
import math

import numpy as np
from sympy import Poly, cyclotomic_poly, symbols

TOL = 1e-6              # floor of the exact-decision gate
EPS = 2.2e-16           # double-precision unit roundoff
GATE_KAPPA = 1024.0     # safety factor between the gate and the float error bound
HP_PREC0 = 113          # first precision (bits) tried by certified_sign


class Ring:
    """Z[zeta_4Q] as reduced integer vectors modulo Phi_4Q (only addition is needed)."""

    def __init__(self, Q):
        self.Q = Q
        n = 4 * Q
        x = symbols('x')
        coeffs = [int(c) for c in Poly(cyclotomic_poly(n, x), x).all_coeffs()][::-1]
        d = len(coeffs) - 1
        self.d = d
        red = np.zeros((n, d), dtype=np.int64)           # red[k] = x^k mod Phi
        cur = np.zeros(d, dtype=np.int64)
        cur[0] = 1
        tail = -np.array(coeffs[:d], dtype=np.int64)
        for k in range(n):
            red[k] = cur
            lead = cur[d - 1]
            cur = np.roll(cur, 1)
            cur[0] = 0
            if lead:
                cur = cur + lead * tail
        # imvals[k] = Im(zeta^k)/4, so a coordinate's real value is dot(v, imvals)
        self.imvals = np.array([math.sin(k * math.pi / (2 * Q)) / 4.0 for k in range(d)])
        self._hp = {}
        self.gen = np.zeros((2 * Q, d), dtype=np.int64)  # gen[c] = zeta^c - zeta^-c
        for c in range(2 * Q):
            self.gen[c] = red[c] - red[(n - c) % n]

    def hp_imvals(self, bits):
        """round(imvals * 2**bits) as exact Python ints, cached."""
        tab = self._hp.get(bits)
        if tab is None:
            from mpmath import mp
            with mp.workprec(bits + 64):
                scale = mp.mpf(2) ** bits
                tab = [int(mp.nint(mp.sin(k * mp.pi / (2 * self.Q)) / 4 * scale))
                       for k in range(self.d)]
            self._hp[bits] = tab
        return tab


class Coord:
    """A transverse coordinate: exact vector for 2x, float for ordering."""
    __slots__ = ('v', 'f')

    def __init__(self, v, f):
        self.v = v
        self.f = f

    def __add__(self, o):
        return Coord(self.v + o.v, self.f + o.f)

    def __sub__(self, o):
        return Coord(self.v - o.v, self.f - o.f)

    def eq(self, o):
        return np.array_equal(self.v, o.v)

    def is_zero(self):
        return not self.v.any()


def half(a):
    """a/2, exact: it is only applied to sums of two ell's, whose vectors are even."""
    assert not (a.v & 1).any(), 'half() on an odd vector'
    return Coord(a.v // 2, a.f / 2.0)


_RING_CACHE = {}


def get_ring(Q):
    if Q not in _RING_CACHE:
        _RING_CACHE[Q] = Ring(Q)
    return _RING_CACHE[Q]


class Chain:
    def __init__(self, P, Q, eps, step_cap=3_000_000):
        assert math.gcd(P, Q) == 1 and 0 < P < Q and eps in (0, 1)
        self.P, self.Q, self.eps, self.step_cap = P, Q, eps, step_cap
        self.even = (Q % 2 == 0)
        twoQ = 2 * Q
        if self.even:
            c0 = P if eps == 0 else (P + 1) % twoQ
        else:
            c0 = P if eps == 0 else (P + Q) % twoQ
        self.c = [(c0 + 2 * P * t) % twoQ for t in range(Q)]
        R = get_ring(Q)
        self.ring = R
        self.ell = [Coord(2 * R.gen[self.c[t]], math.sin(self.c[t] * math.pi / twoQ))
                    for t in range(Q)]
        self.zero = Coord(np.zeros(R.d, dtype=np.int64), 0.0)
        self.m_idx = None if self.even else (Q - 1) // 2
        self.far, self.near, self.osec, self.tie = {}, {}, {}, {}
        for k in range(Q):
            lo, hi = (k - 1) % Q, k
            a, b = self.ell[lo], self.ell[hi]
            self.tie[k] = a.eq(b)
            if a.f >= b.f:
                self.far[k], self.near[k] = lo, hi
            else:
                self.far[k], self.near[k] = hi, lo
            self.osec[k] = 0 < (twoQ - self.c[lo]) % twoQ < 2 * P
        self.delta, self.flo, self.fhi, self.fsum, self.ffix = {}, {}, {}, {}, {}
        for k in range(Q):
            M, m = self.ell[self.far[k]], self.ell[self.near[k]]
            self.delta[k] = M - m
            lo, hi = (self.zero, M - m) if self.osec[k] else (m, M)
            self.flo[k], self.fhi[k] = lo, hi
            self.fsum[k] = lo + hi
            self.ffix[k] = half(self.fsum[k])
        self.exact_calls = 0
        self.sign_calls = 0
        self.maxcoef = 0
        self.min_gate_ratio = float('inf')
        self._gate = TOL

    # -- exact-when-close comparison ----------------------------------------------------
    def gate(self, mag):
        """Half-width inside which a comparison is decided exactly: GATE_KAPPA times the
        float error bound of a coordinate whose largest coefficient is `mag`."""
        return max(TOL, GATE_KAPPA * EPS * self.ring.d * max(mag, 1))

    def certified_sign(self, w):
        """Sign of dot(w, imvals), proved.  Every coordinate is purely imaginary in
        Z[zeta_4Q], and Im is injective there, so w != 0 has a nonzero value and the
        precision loop terminates."""
        mx = int(np.abs(w).max())
        if mx == 0:
            return 0
        nz = [int(k) for k in np.nonzero(w)[0]]
        wl = [int(w[k]) for k in nz]
        bits = HP_PREC0
        while bits <= 1 << 16:
            tab = self.ring.hp_imvals(bits)
            s = 0
            for k, c in zip(nz, wl):
                s += c * tab[k]
            # s is exact in units of 2**-bits; each table entry is off by <= 1/2 unit
            if abs(s) > (self.ring.d * mx) // 2 + 1:
                return 1 if s > 0 else -1
            bits *= 2
        raise RuntimeError('certified_sign did not converge')

    def cmp(self, a, b):
        d = a.f - b.f
        g = self._gate
        if -g < d < g:
            self.exact_calls += 1
            if a.eq(b):
                return 0
            self.sign_calls += 1
            return self.certified_sign(a.v - b.v)
        r = abs(d) / g
        if r < self.min_gate_ratio:
            self.min_gate_ratio = r
        return -1 if d < 0 else 1

    def eq_exact(self, a, b):
        d = a.f - b.f
        g = self._gate
        if not (-g < d < g):
            return False
        self.exact_calls += 1
        return a.eq(b)

    # -- the kite maps ------------------------------------------------------------------
    def other_edge(self, k, j):
        lo, hi = (k - 1) % self.Q, k
        return hi if j == lo else lo

    def next_kite(self, k, j):
        return (k - 1) % self.Q if j == (k - 1) % self.Q else (k + 1) % self.Q

    def cross(self, k, j, x):
        """Enter kite k through interface j at coordinate x.
        Returns ('pole', k) or (j', x', transited)."""
        if self.tie[k]:
            return self.other_edge(k, j), x, True
        if j == self.near[k]:
            x2 = x + self.delta[k] if self.osec[k] else x
            return self.far[k], x2, True
        if self.cmp(x, self.flo[k]) >= 0 and self.cmp(x, self.fhi[k]) <= 0:
            if self.cmp(x, self.ffix[k]) == 0:
                return 'pole', k, False
            return j, self.fsum[k] - x, False
        x2 = x - self.delta[k] if self.osec[k] else x
        return self.near[k], x2, True

    def trace(self, j, x, k):
        """Follow the leaf from interface point (j, x) into kite k until it arrives at a
        singularity or the step cap.  Returns (end, steps); end is ('R', k), ('ZO', k),
        ('ZA', k) or 'CAP'."""
        steps = 0
        imv = self.ring.imvals
        while True:
            steps += 1
            if steps > self.step_cap:
                return 'CAP', steps
            if not steps & 511:
                # refresh the float from the exact vector, and track the magnitude
                x.f = float(x.v @ imv)
                a = int(np.abs(x.v).max())
                self._gate = self.gate(a)
                if a > self.maxcoef:
                    self.maxcoef = a
            r = self.cross(k, j, x)
            if r[0] == 'pole':
                return ('R', r[1]), steps
            j2, x2, _ = r
            if -self._gate < x2.f < self._gate and x2.is_zero():
                return ('ZO', k), steps
            if self.eq_exact(x2, self.ell[j2]):
                return ('ZA', k), steps
            k = self.next_kite(k, j2)
            j, x = j2, x2

    def separatrices(self):
        """Every separatrix of B in this direction.  Those lying on the horizontal arcs of
        the fixed set of the reversing involution are saddle connections by construction
        and are listed with `in_fix` and not walked; every other one is walked, from each
        of its two ends (so each such saddle connection is found twice)."""
        Q = self.Q
        seps, starts = [], []
        if self.even:
            for k in range(Q):
                if self.tie[k]:
                    seps.append({'ends': [('R', k), ('ZO' if k == 0 else 'ZA', k)],
                                 'in_fix': True, 'steps': 0})
            for t in range(Q):
                if self.ell[t].is_zero():
                    seps.append({'ends': [('ZO', None), ('ZA', None)], 'in_fix': True,
                                 'steps': 0})
        else:
            seps.append({'ends': [('R', 0), ('ZO' if self.eps == 0 else 'ZA', 0)],
                         'in_fix': True, 'steps': 0})
            if self.ell[self.m_idx].is_zero():
                seps.append({'ends': [('ZO', None), ('ZA', None)], 'in_fix': True,
                             'steps': 0})
        for k in range(Q):
            if self.tie[k]:
                continue
            lo, hi, fix = self.flo[k], self.fhi[k], self.ffix[k]
            far = self.far[k]
            starts.append((far, fix, k, ('R', k)))                    # the pole's prong
            if self.ell[self.near[k]].is_zero():
                continue                    # the corner prongs here ARE the horizontal arc
            starts.append((far, hi if self.osec[k] else lo, k,
                           ('ZO' if self.osec[k] else 'ZA', k)))      # the corner's prong
        for j, x, k, tag in starts:
            end, steps = self.trace(j, x, self.next_kite(k, j))
            seps.append({'ends': [tag, end], 'in_fix': False, 'steps': steps})
        return seps


def certify(P, Q, eps, step_cap=3_000_000, with_seps=False):
    """The closure certificate of one perpendicular direction class.

    Returns a dict: `cp` (True iff every walked separatrix arrived -- PROVES complete
    periodicity; False means only "undecided at this cap"), `n_walked`, `capped`,
    `h` (separatrices on the fixed arcs, not walked), `max_steps` (longest walk that
    arrived), `tot_steps`, `E` (number of saddle connections, when cp), and the exact-path
    counters."""
    ch = Chain(P, Q, eps, step_cap)
    seps = ch.separatrices()
    walked = [s for s in seps if not s['in_fix']]
    fixed = [s for s in seps if s['in_fix']]
    capped = sum(1 for s in walked if s['ends'][1] == 'CAP')
    arrived = [s['steps'] for s in walked if s['ends'][1] != 'CAP']
    out = {'P': P, 'Q': Q, 'eps': eps, 'cap': step_cap, 'cp': capped == 0,
           'n_walked': len(walked), 'capped': capped, 'h': len(fixed),
           'max_steps': max(arrived, default=0),
           'tot_steps': sum(s['steps'] for s in walked),
           'E': (len(fixed) + len(walked) // 2) if capped == 0 else None,
           'exact_calls': ch.exact_calls, 'sign_calls': ch.sign_calls,
           'safety_certified': GATE_KAPPA * ch.min_gate_ratio}
    if with_seps:
        out['seps'] = seps
    return out
