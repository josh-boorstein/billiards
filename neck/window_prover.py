"""A mechanical prover for bounded interval-model walks (the necklace paper, Section 10).

For a family (P, r, side) -- Q = mP + r, class eps = 0, and `side` the NEAR arc (the L1 copy in
kite 0) or the FAR arc (the last vertical arc) -- every leaf from the arc's interval is followed
into the chain D symbolically, and the cell count is certified for ALL large m.

Method.  The walk is run at a concrete m0 to fix the discrete path; each coordinate is carried
as +-Y + (a Z-combination of s_j = sin(j x)), Y the start point.  Along a monotone run the
intervals nest (Proposition 7.2), so a stretch of crossings is decided by the interval at its
lowest point: the emitted constraints are y in I_t at every local minimum passed and at the last
interface before a fold, and y not in I_next at the fold.  Every interface that appears must have
c_t within JMAX of 0 or 2Q (so sin(pi c_t / 2Q) = s_j with j small and m-independent) -- otherwise
the family is reported UNBOUNDED and nothing is claimed.  Bands are the maximal Y-intervals with
one decision sequence; their endpoints are the tight constraints, solved symbolically.  Each
constraint is checked at the band's two endpoints (it is affine in Y), each check an odd
polynomial in S = sin x certified root-free on (0, sin(pi/2Q0)] by exact real-root isolation.
Cells = bands + interior poles.  The symbolic output is required to be IDENTICAL at m0, m0+1,
m0+2, m0+3 (the m-independence gate).  Fold-kite orientation is fixed in float at m0; it is
uniform because every listed interface index is <= 23 < Q.  `arm_S` checks the skeleton the
proof reads -- the local minima and the O-sector kites -- over a range of m; it is not the
m-independence argument, which is the paper's explicit interface list (Theorem 10.7's proof).
"""
from __future__ import annotations

import json
import math
import os
import sys

import mpmath as mp
import sympy as sp

mp.mp.dps = 60
JMAX = 80
X = sp.symbols("x", positive=True)
Y = sp.symbols("Y")


def s(j):
    return sp.sin(j * X)


class Fam:
    def __init__(self, P, r, m):
        self.P, self.r, self.m = P, r, m
        Q = self.Q = m * P + r
        assert math.gcd(P, Q) == 1
        tQ = 2 * Q
        self.even = Q % 2 == 0
        self.M = Q // 2 - 1 if self.even else (Q - 1) // 2
        self.c = [(P + 2 * P * t) % tQ for t in range(Q)]
        self.x = mp.pi / tQ
        self.j = [min(ci, tQ - ci) for ci in self.c]
        self.osec = {k: 0 < (tQ - self.c[(k - 1) % Q]) % tQ < 2 * P for k in range(1, self.M + 1)}
        # numeric and symbolic interval endpoints, with u_0 = 0
        self.lo_n = {0: mp.mpf(0)}
        self.lo_s = {0: sp.Integer(0)}
        for k in range(1, self.M + 1):
            a, b = k - 1, k                      # kite k between interfaces a and b
            ea, eb = self.eln(a), self.eln(b)
            if ea >= eb:                         # far = a (known), near = b
                if self.osec[k]:
                    self.lo_n[b] = self.lo_n[a] + (ea - eb)
                    self.lo_s[b] = self.lo_s[a] + (self.els(a) - self.els(b))
                else:
                    self.lo_n[b], self.lo_s[b] = self.lo_n[a], self.lo_s[a]
            else:                                # far = b, near = a known: u_a = u_b + d (O) or u_b
                if self.osec[k]:
                    self.lo_n[b] = self.lo_n[a] - (eb - ea)
                    self.lo_s[b] = self.lo_s[a] - (self.els(b) - self.els(a))
                else:
                    self.lo_n[b], self.lo_s[b] = self.lo_n[a], self.lo_s[a]

    def eln(self, t):
        return mp.sin(self.c[t] * self.x)

    def els(self, t):
        return s(self.j[t])

    def small(self, t):
        return self.j[t] <= JMAX

    def hi_n(self, t):
        return self.lo_n[t] + self.eln(t)

    def hi_s(self, t):
        return self.lo_s[t] + self.els(t)


def walk(F, y0, side, sym=True, maxsteps=200000):
    """Follow the leaf from start point `y0`.  Returns (fate, (constraints, folds)).  With
    `sym=False` only the fate and the fold kites are computed (fast; used to find the bands).
    Constraints are symbolic `expr > 0` in `Y`; folds are (kite, fixed point, y before)."""
    yn = y0
    ys = Y if sym else None
    t, d = (0, +1) if side == "near" else (F.M, -1)
    cons, folds, minima, run_down = [], [], [], None

    def member(tt):
        if sym:
            cons.append((ys - F.lo_s[tt], f"in I{tt} lo"))
            cons.append((F.hi_s(tt) - ys, f"in I{tt} hi"))

    for _ in range(maxsteps):
        k = t if d == -1 else t + 1
        if k == 0 or k == F.M + 1:
            fate = "N" if k == 0 else "F"
            break
        nxt = k - 1 if d == -1 else k
        if F.eln(nxt) >= F.eln(t):                 # uphill: forced
            if run_down is not None:
                minima.append(run_down); run_down = None
            t = nxt
            continue
        if F.lo_n[nxt] < yn < F.hi_n(nxt):          # downhill crossing
            run_down = nxt
            t = nxt
            continue
        # fold in kite k on the tail of interface t
        for tt in minima + [t, nxt]:
            if not F.small(tt):
                return "UNBOUNDED", (tt, F.j[tt])
        for tt in set(minima + [t]):
            member(tt)
        minima, run_down = [], None
        above = yn >= F.hi_n(nxt)
        if above:
            rlo_n, rhi_n = F.hi_n(nxt), F.hi_n(t)
        else:
            rlo_n, rhi_n = F.lo_n[t], F.lo_n[nxt]
        if sym:
            if above:
                cons.append((ys - F.hi_s(nxt), f"out I{nxt} above"))
                rlo_s, rhi_s = F.hi_s(nxt), F.hi_s(t)
            else:
                cons.append((F.lo_s[nxt] - ys, f"out I{nxt} below"))
                rlo_s, rhi_s = F.lo_s[t], F.lo_s[nxt]
            folds.append((k, (rlo_s + rhi_s) / 2, ys))
            ys = sp.expand(rlo_s + rhi_s - ys)
        else:
            folds.append((k, None, None))
        yn = rlo_n + rhi_n - yn
        d = -d
    else:
        return "CAP", None
    if run_down is not None:
        minima.append(run_down)
    for tt in minima:
        if not F.small(tt):
            return "UNBOUNDED", (tt, F.j[tt])
        member(tt)
    return fate, (cons, folds)


def start_interval(F, side):
    if side == "near":
        return 0, (F.lo_s[0], F.hi_s(0)), (F.lo_n[0], F.hi_n(0))
    # far: the window is the last local minimum's interval; the tails are handled analytically
    t = F.M
    while t > 0 and F.eln(t - 1) < F.eln(t):
        t -= 1
    return t, (F.lo_s[t], F.hi_s(t)), (F.lo_n[t], F.hi_n(t))


def bands(F, side, N=4000):
    tw, (Ls, Hs), (Ln, Hn) = start_interval(F, side)
    out, prev = [], None
    for i in range(1, N):
        y0 = Ln + (Hn - Ln) * i / N
        fate, data = walk(F, y0, side, sym=False)
        if fate in ("UNBOUNDED", "CAP"):
            return tw, None, (fate, data)
        sig = (fate, tuple(f[0] for f in data[1]))
        if prev is None or prev[0] != sig:
            prev = [sig, y0, y0]
            out.append(prev)
        else:
            prev[2] = y0
    for b in out:                                   # symbolic walk once per band
        ymid = (b[1] + b[2]) / 2
        fate, (cons, folds) = walk(F, ymid, side, sym=True)
        b.extend([cons, folds])
    return tw, out, (Ls, Hs)


def signature(F, side):
    tw, bl, extra = bands(F, side)
    if bl is None:
        return None, extra
    # m-independence is judged on the INTERFACE INDICES of each fold kite, (j_{k-1}, j_k): the
    # fold kites sit at the profile's minima, which move with m relative to both ends.
    rel = lambda k: (F.j[k - 1], F.j[k])
    return [(b[0][0], tuple(rel(k) for k in b[0][1])) for b in bl], (tw, bl, extra)


def to_poly(e, Ssym):
    p = sp.expand(sp.expand_trig(e))
    p = sp.expand(p.subs(sp.cos(X) ** 2, 1 - sp.sin(X) ** 2))
    p = sp.expand(p.replace(lambda b: b.is_Pow and b.base == sp.cos(X),
                            lambda b: (1 - sp.sin(X) ** 2) ** (b.exp // 2) * sp.cos(X) ** (b.exp % 2)))
    p = sp.expand(p.subs(sp.sin(X), Ssym))
    if p.has(sp.cos(X)):
        raise ValueError("not an odd polynomial in sin x")
    return sp.Poly(p, Ssym)


def q0_for(e):
    """least Q0 such that e > 0 on (0, sin(pi/2Q0)]; None if e <= 0 near 0; 0 if e == 0."""
    Ssym = sp.symbols("S")
    if sp.simplify(e) == 0:
        return 0
    q = to_poly(e, Ssym)
    if q.is_zero:
        return 0
    while q.eval(0) == 0:
        q = sp.Poly(sp.cancel(q.as_expr() / Ssym), Ssym)
    if sp.sign(q.eval(sp.Rational(1, 10**7))) <= 0:
        return None
    pos = sorted(float(r) for r in q.real_roots() if r > 0)
    if not pos:
        return 1
    v = math.pi / (2 * math.asin(min(pos[0], 1.0)))
    # the root may be EXACTLY sin(pi/2N) (an identity at Q = N): then Q = N fails and N + 1 is least
    return int(round(v)) + 1 if abs(v - round(v)) < 1e-9 else math.floor(v) + 1


def prove(P, r, side, m0=60, verbose=True):
    """Return (cells_in_start_interval, Q0, (window t, M), number of checks) or None."""
    sigs = []
    for dm in range(4):
        F = Fam(P, r, m0 + dm)
        sig, extra = signature(F, side)
        if sig is None:
            if verbose:
                print(f"  ({P},{r},{side}) m={m0 + dm}: {extra}")
            return None
        sigs.append(sig)
    if any(sg != sigs[0] for sg in sigs):
        if verbose:
            print(f"  ({P},{r},{side}): signature depends on m: {sigs}")
        return None
    F = Fam(P, r, m0)
    tw, bl, (Ls, Hs) = bands(F, side)
    # band endpoints: tight constraints; build symbolic endpoints by solving at the numeric split
    ends = [Ls]
    for b1, b2 in zip(bl, bl[1:]):
        # the split value: some constraint of b1 or b2 changes sign; find the constraint of b2
        # (or b1) whose zero lies between b1's last and b2's first sample
        ya, yb = b1[2], b2[1]
        found = None
        for cons in (b1[3], b2[3]):
            for e, tag in cons:
                sol = sp.solve(sp.Eq(e, 0), Y)
                if not sol:
                    continue
                z = sol[0]
                zn = mp.mpf(sp.N(z.subs(X, F.x), 50))
                if ya - mp.mpf(10) ** -30 <= zn <= yb + mp.mpf(10) ** -30:
                    found = sp.expand(z)
                    break
            if found is not None:
                break
        if found is None:
            if verbose:
                print("  cannot identify a split point"); return None
        ends.append(found)
    ends.append(Hs)
    worst = 0
    cells = 0
    checks = []
    for i, b in enumerate(bl):
        lo_e, hi_e = ends[i], ends[i + 1]
        checks.append((f"band {i} nonempty", hi_e - lo_e))
        for e, tag in b[3]:
            for E, nm in ((e.subs(Y, lo_e), "lo"), (e.subs(Y, hi_e), "hi")):
                checks.append((f"band {i} {tag} @ {nm}", E))
        cells += 1
        for (k, fix, ybefore) in b[4]:
            # does the image of the band under the pre-fold map contain the fixed point?
            a_, b_ = ybefore.subs(Y, lo_e), ybefore.subs(Y, hi_e)
            fa = float(sp.N((a_ - fix).subs(X, F.x)))
            fb = float(sp.N((b_ - fix).subs(X, F.x)))
            if fa * fb < 0:
                cells += 1
                lo_, hi_ = (a_, b_) if fa < 0 else (b_, a_)
                checks.append((f"band {i} pole of kite {k} strictly inside (lo)", fix - lo_))
                checks.append((f"band {i} pole of kite {k} strictly inside (hi)", hi_ - fix))
            else:
                far_end = a_ if abs(fa) < abs(fb) else b_
                sgn = 1 if (fa > 0) else -1
                checks.append((f"band {i} misses pole of kite {k}", sgn * (far_end - fix)))
    for name, E in checks:
        q0 = q0_for(sp.expand(E))
        if q0 is None:
            if verbose:
                print(f"    FAIL {name}: nonpositive near 0")
            return None
        if q0 == 0:
            # an identity equality: allowed only for a constraint evaluated at a band endpoint (the
            # tight constraint that DEFINES that endpoint); a zero "nonempty" or pole check fails
            if "@ lo" in name or "@ hi" in name:
                continue
            if verbose:
                print(f"    FAIL {name}: identically zero")
            return None
        worst = max(worst, q0)
    if verbose:
        print(f"  ({P},{r},{side}): bands {[(b[0][0], len(b[0][1])) for b in bl]}  cells in start "
              f"interval = {cells}; {len(checks)} checks certified for Q ≥ {worst}; window at "
              f"t = {tw} (M = {F.M})")
    return cells, worst, (tw, F.M), len(checks)


def skeleton(P, r, m):
    """The data the symbolic proof reads: the local minima of the path profile and the O-sector
    kites, each with its interface indices, in order; the far arc's tail length is M - t_last."""
    Q = m * P + r
    tQ = 2 * Q
    M = Q // 2 - 1 if Q % 2 == 0 else (Q - 1) // 2
    c = [(P + 2 * P * t) % tQ for t in range(Q)]
    j = [min(ci, tQ - ci) for ci in c]
    ell = [math.sin(math.pi * ci / tQ) for ci in c]
    mins = tuple((j[t - 1], j[t], j[t + 1]) for t in range(1, M) if ell[t] < ell[t - 1] and ell[t] < ell[t + 1])
    osec = tuple((j[k - 1], j[k]) for k in range(1, M + 1) if 0 < (tQ - c[k - 1]) % tQ < 2 * P)
    return mins, osec, j[0], j[M]


def arm_S(P, r, m_lo, m_hi=2000):
    ref = skeleton(P, r, m_lo)
    bad = [m for m in range(m_lo, m_hi + 1) if math.gcd(P, m * P + r) == 1
           and skeleton(P, r, m)[:2] != ref[:2]]
    print(f"  S: ({P},{r}) skeleton for {m_lo} ≤ m ≤ {m_hi}: {'identical' if not bad else 'DIFFERS at m=' + str(bad[:5])}"
          f"  minima {ref[0]}  O-sector {ref[1]}")
    return not bad


