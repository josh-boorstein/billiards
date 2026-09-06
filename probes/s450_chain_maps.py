#!/usr/bin/env python3.13
"""s450_chain_maps.py -- PRE-REGISTERED.  Queue item (s'): the necklace's INTERVAL MAPS.

PRIOR ART: `rulings.py --grep` on 'interval exchange' -> NOTHING; 'separatrix diagram' ->
[NCYL-280] (`E = Q`, the genus-0 base, the a/b/d/kappa edge bookkeeping), [OPS-214],
[OPS-218] (`E = h + f + 2p`); 'bigon' -> [NCYL-280], [NCYL-287] (the necklace model),
[NCYL-291] (the `+P` gluing order + the coordinates); 'transverse measure'/'rotation by
P/Q' -> [NCYL-291] only.  s446 computes `p` by TRACING billiard prongs (caps `Q <= 30`,
`data/s446_covering_q30.json`, `527` scored rows); s449 pins the necklace's COORDINATES
but implements no dynamics ("the per-kite flat geometry" is explicitly left undone).
No probe in the repo implements the foliation of `B` as interval maps.
=> what is NEW here is the per-kite FLAT GEOMETRY and the resulting tracer; the necklace,
its gluing order, its coordinates and `E = h + f + 2p` are all prior art and are USED.

------------------------------------------------------------------------------------------
THE PER-KITE FLAT GEOMETRY (derived here; this is the piece s449 left open).

Set-up is [NCYL-291]'s verbatim.  `B` = `Q` kites glued along `H`-edges; kite `t` sits
between interface `t-1` and interface `t`; interface `t` has direction `a + P(2t+1)`.  The
deck makes all `Q` directions of one parity equivalent, so WLOG `psi = eps`, and then
`a = eps` (`eps = 0`) or `eps + Q` (`eps = 1`), giving the offsets from `psi`

    c[t] = (c0 + 2*P*t) mod 2Q,   c0 = P (eps = 0) / P + Q (eps = 1),
    ell[t] = sin(c[t]*pi/(2Q))  >= 0   -- the interface's TRANSVERSE measure.

`c[t]` runs over ALL residues of one parity, so exactly one interface is horizontal
(`ell = 0`, lone class) or perpendicular (`ell = 1`, paired class), and it is `t = m`.

A kite is the rhombus of 4 triangle copies about `R` modulo the half-turn.  Put the two
side vectors at `-alpha` and `+alpha` (`2alpha = P*pi/Q` = the `Z_O` angle) and sweep by
the horizontal direction.  The four vertices `O, A, O', A'` split the sweep into THREE
bands, and the half-turn identifies the outer two.  What comes out, with `M`/`m` the
larger/smaller of the two interface measures and the coordinate on every interface being
the TRANSVERSE distance from `Z_O`:

  (i)   the fold lives on the LARGER ("far") edge; the smaller ("near") edge is entirely
        a THROUGH band, mapped bijectively onto a sub-interval of the far edge;
  (ii)  WHICH sub-interval is decided by which corner's angular sector contains `psi`.
        Kite `t` occupies the sector `[c[t-1], c[t]]` around `Z_O`, so
            O-SECTOR  (`0` in that arc): through = far`[M-m, M]`, `near x -> x + (M-m)`;
                                         fold = far`[0, M-m]`, `x -> (M-m) - x`;
            A-SECTOR  (otherwise)      : through = far`[0, m]`,   `near x -> x`;
                                         fold = far`[m, M]`,     `x -> (m+M) - x`;
  (iii) the pole `R_t`'s single prong is the FOLD'S FIXED POINT, `(M-m)/2` resp. `(M+m)/2`;
  (iv)  the fold interval's other endpoint IS the sector corner, so the `Z_O` (resp `Z_A`)
        separatrix into kite `t` is the fold image of that corner -- which is why `Z_O` has
        exactly one prong per O-sector kite, i.e. `P` of them, and `Z_A` has `Q-P`.

  ⇒ TWO DEGENERACIES, both forced and both self-consistent:
  (v)   `M = m` happens iff the kite's bisector is horizontal or vertical, i.e. iff
        `c[t-1] + P == 0 mod Q` -- exactly ONE kite, and it is always kite `0`.  There the
        fold is empty and `R_0`'s prong runs along the horizontal diagonal to `Z_O`
        (`eps = 0`) or `Z_A` (`eps = 1`): that IS [NCYL-283]'s horizontal `Fix(iota)` arm.
  (vi)  `m = 0` (lone class) at the two kites flanking the horizontal interface: the through
        band is empty, the fold is the whole far edge, and the O-/A-sector rules COINCIDE.
        The horizontal interface is itself a `Z_O`-`Z_A` saddle connection -- the second
        horizontal `Fix` arm, which is why `h = 2` exactly in the lone class.

  ⇒ `Fix(iota)` in these coordinates: kite `0`'s two diagonal arms (`Z_O`-`R_0` and
    `R_0`-`Z_A`, one horizontal one vertical) plus interface `m`.  The VERTICAL arm of
    kite `0` is crossed by EVERY leaf that transits kite `0` (a tie kite has no fold, and
    the vertical diagonal spans the whole sweep), and it is the `L1`-copy (`eps = 0`) resp.
    the `L2`-copy (`eps = 1`).  Interface `m` is the `H`-copy, vertical in the paired class.

  ⇒ therefore, with `G` the separatrix diagram, `f` = #separatrices meeting a vertical arc,
    `h` = #separatrices lying IN `Fix`, and [OPS-218]'s `Q = h + f + 2p`:
        p  =  (Q - h - f)/2,      n_sigma  =  (#separatrices meeting sigma's arc) + 1.

ARITHMETIC.  Every coordinate is a Z-combination of the `ell`'s with half-integer
coefficients, so `2x` is exactly an element of `Z[zeta_4Q]` (`2i*sin(c*pi/2Q) = z^c - z^-c`).
Coordinates are carried as REDUCED integer vectors mod `Phi_{4Q}` (dimension `2*phi(Q)`)
alongside a float: EQUALITY (= a leaf hitting a singularity, which is the generic event
here, not a rare one) is decided EXACTLY on the vector, and the float is used only for the
`<`/`>` band test, whose margin is monitored.  ⚠ [OPS-215]: no float arrival test anywhere.

HYPOTHESES.  ⚠ Read [OPS-041]/[OPS-222] first: an arm that cannot fail is not evidence.
  H1  ⇒⇒ THE CROSS-INSTRUMENT TEST, AND IT CAN FAIL ON EVERY ROW AND IN EVERY COLUMN.
      Against s446's `527` TRACED rows (`Q <= 30`): the per-side cell counts `n_sigma`
      (1 or 2 integers per row, values `2..30`), `f`, `h`, `v` and `n_nofold = 2p`.  s446
      traces billiard prongs with `right_triangle_billiards`; this file does pure kite
      arithmetic and shares NO code, NO instrument and NO input with it.  A wrong sector
      rule, a wrong fold end, a wrong gluing order or a wrong `M`/`m` choice all break
      `n_sigma` immediately.  ⚠ `n_nofold == 0` on all 527 is NOT the test (it is constant
      -- [OPS-041]); the discriminating columns are the `n_sigma` values.
  H2  CONTROL, must FAIL.  Re-run H1 with a WRONG gluing order (`+1`, `+2`) and with a
      WRONG sector partition (`shift`, `parity`).  ⚠⚠ AND ONE CONTROL IS VACUOUS BY
      CONSTRUCTION AND SCORES A PERFECT `347/347` -- `sector='complement'`; see
      `wrong_sector()`, which owns the argument, and [OPS-223].  It is kept, and kept
      labelled, because a control that is a SYMMETRY OF THE READOUT is the failure mode
      worth having on the page.
  H3  STRUCTURE, free.  Every separatrix meets the vertical arcs at most ONCE (s446 step
      (3)); `#separatrices == Q`; prong-ends `== 2Q`; `h` matches [NCYL-283]'s parity rule.
  H5  THE PAIRING, and it is a RESTATEMENT not a test.  Every traced separatrix joins the
      prong of kite `k` to the prong of kite `-k`, same type -- `158/158` rows, `Q <= 19`.
      ⚠⚠ DO NOT QUOTE IT AS EVIDENCE FOR `p = 0`: `p` counts iota-SWAPPED pairs, so for
      `k != 0` the ends `{k, sigma(k)}` being iota-stable forces `sigma(k) = -k` and
      conversely -- the pairing IS `p = 0` ([OPS-041]/[OPS-222], [NCYL-292]).  It is
      recorded because it is the form the open residual is now stated in.
  H4  COVERAGE, the point of the exercise.  `p == 0` above `Q = 30`.  ⚠ COVERAGE IS NOT A
      PROOF -- `p = 0` stays open whatever this returns ([NCYL-287]).

  python3.13 probes/s450_chain_maps.py [validate|sweep] [--qmax=N]
"""
from __future__ import annotations

import json
import math
import sys
import time

import numpy as np
from sympy import Poly, cyclotomic_poly, symbols

SIDES = ("L2", "L1", "H")
STEP_CAP = 20_000_000
TOL = 1e-6            # float gate for the exact vector tests; `margin` audits it


# ---------------------------------------------------------------- exact ring scaffolding

class Ring:
    """`Z[zeta_{4Q}]` as reduced integer vectors mod `Phi_{4Q}`; only ADDITION is needed."""

    def __init__(self, Q):
        n = 4 * Q
        x = symbols("x")
        coeffs = [int(c) for c in
                  Poly(cyclotomic_poly(n, x), x).all_coeffs()][::-1]     # ascending
        d = len(coeffs) - 1
        self.d = d
        red = np.zeros((n, d), dtype=np.int64)      # red[k] = x^k mod Phi
        cur = np.zeros(d, dtype=np.int64)
        cur[0] = 1
        tail = -np.array(coeffs[:d], dtype=np.int64)   # x^d == tail
        for k in range(n):
            red[k] = cur
            lead = cur[d - 1]
            cur = np.roll(cur, 1)
            cur[0] = 0
            if lead:
                cur = cur + lead * tail
        # imvals[k] = Im(zeta^k)/4, so a coordinate's real value is dot(v, imvals)
        self.imvals = np.array([math.sin(k * math.pi / (2 * Q)) / 4.0
                                for k in range(d)], dtype=float)
        # gen[c] = zeta^c - zeta^-c  == 2i*sin(c*pi/(2Q))
        self.gen = np.zeros((2 * Q, d), dtype=np.int64)
        for c in range(2 * Q):
            self.gen[c] = red[c] - red[(n - c) % n]


class Coord:
    """A transverse coordinate: exact vector for `2x`, float for ordering."""
    __slots__ = ("v", "f")

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
    """`a/2`.  Every `ell` enters with an EVEN vector (`2*gen`), so sums of `ell`s are
    even and this division is exact; it is only ever applied to `lo + hi`."""
    assert not (a.v & 1).any(), "half() on an odd vector -- would lose exactness"
    return Coord(a.v // 2, a.f / 2.0)


_RING_CACHE = {}


def get_ring(Q):
    if Q not in _RING_CACHE:
        _RING_CACHE[Q] = Ring(Q)
    return _RING_CACHE[Q]


# ---------------------------------------------------------------- the model

def wrong_sector(osec, Q, P, mode):
    """Deliberately wrong sector assignments, for H2.

    ⚠ `complement` IS VACUOUS AND IS KEPT ONLY TO RECORD THAT ([OPS-219]).  Flipping
    EVERY kite is the global relabelling `Z_O <-> Z_A`, i.e. the coordinate flip
    `x -> ell - x` on every interface simultaneously, so it conjugates the whole system
    and every count this file reads (`f`, `n_sigma`, `p`) is invariant.  It scored
    `160/162` -- exactly the true model -- which is what a control measuring nothing
    looks like.  The sector rule's content is the PARTITION of the kites, so a control
    must move kites BETWEEN the classes: `shift` and `parity` do.
    """
    if mode == "complement":
        return {k: not v for k, v in osec.items()}
    if mode == "shift":                       # same number of O-sector kites, moved by 1
        return {k: osec[(k - 1) % Q] for k in range(Q)}
    if mode == "parity":                      # right count is P; this is ~Q/2
        return {k: (k % 2 == 0) for k in range(Q)}
    raise ValueError(mode)


class Chain:
    def __init__(self, P, Q, eps, order=None, sector=None):
        assert Q % 2 == 1 and math.gcd(P, Q) == 1 and 0 < P < Q
        self.P, self.Q, self.eps = P, Q, eps
        self.lone = (eps == P % 2)
        twoQ = 2 * Q
        step = 2 * (P if order is None else order)
        c0 = P if eps == 0 else (P + Q) % twoQ
        self.c = [(c0 + step * t) % twoQ for t in range(Q)]
        R = get_ring(Q)
        self.ring = R
        # ell[t] as a Coord (vector carries 2*ell)
        self.ell = [Coord(2 * R.gen[self.c[t]], math.sin(self.c[t] * math.pi / twoQ))
                    for t in range(Q)]
        self.zero = Coord(np.zeros(R.d, dtype=np.int64), 0.0)
        self.m_idx = (Q - 1) // 2
        # per-kite data
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
        if sector is not None:
            self.osec = wrong_sector(self.osec, Q, P, sector)
        # per-kite CONSTANTS, precomputed once (the hot loop must not allocate these)
        self.delta, self.flo, self.fhi, self.fsum, self.ffix = {}, {}, {}, {}, {}
        for k in range(Q):
            M, m = self.ell[self.far[k]], self.ell[self.near[k]]
            self.delta[k] = M - m
            lo, hi = (self.zero, M - m) if self.osec[k] else (m, M)
            self.flo[k], self.fhi[k] = lo, hi
            self.fsum[k] = lo + hi
            self.ffix[k] = half(self.fsum[k])
        self.margin = 1.0
        self.exact_calls = 0
        self.maxcoef = 0

    def safety(self):
        """`margin / (worst-case float error)` -- THE AUDIT OF THE FLOAT GATE.

        Equality is always exact, but the `<`/`>` band test reads the float, whose error
        is bounded by `eps * d * max|coef|` (the coordinate's float is recomputed from the
        exact vector every 512 steps, so this is a bound on the DOT PRODUCT, not on an
        accumulation).  `margin` is the smallest gap ever resolved as non-equal.  ⚠ A row
        is only as trustworthy as this ratio: `< 1` means a band test could have taken the
        wrong branch.  Measured on the s450 sweep: `>= 20` on all `234` SCORED rows
        (worst `10/37` eps1, `1.85e-08` against `2.3e-10`), while every ratio below that
        sits on a CAPPED row, which is unscorable anyway.  ⚠ If this ever approaches `1`,
        the fix is a high-precision sign on the gated comparisons -- they are RARE (`80`
        exact calls in `2.09e6` steps), so it costs almost nothing; it is not written
        because it has never been needed and an unexercised code path is worse.
        """
        err = 2.2e-16 * self.ring.d * max(self.maxcoef, 1)
        return self.margin / err if err else float("inf")

    # -- geometry helpers -------------------------------------------------
    def fold_lo_hi(self, k):
        return self.flo[k], self.fhi[k], self.ffix[k]

    def cmp(self, a, b):
        """`-1/0/+1`.  FLOAT-GATED: the exact vector test runs only when the floats are
        within `TOL`; outside that the float sign is certain by many orders of magnitude.
        `self.margin` records the smallest float gap ever resolved as NON-equal, which is
        the number that says whether the gate was safe ([OPS-215] -- no float arrival
        test: an ARRIVAL is always decided on the exact vector)."""
        d = a.f - b.f
        if -TOL < d < TOL:
            self.exact_calls += 1
            if a.eq(b):
                return 0
        if d < 0:
            if -d < self.margin:
                self.margin = -d
            return -1
        if d < self.margin:
            self.margin = d
        return 1

    def eq_exact(self, a, b):
        """Exact equality, float-gated the same way."""
        d = a.f - b.f
        if not (-TOL < d < TOL):
            return False
        self.exact_calls += 1
        return a.eq(b)

    def other_edge(self, k, j):
        lo, hi = (k - 1) % self.Q, k
        return hi if j == lo else lo

    def next_kite(self, k, j):
        """the kite on the far side of edge `j` from kite `k`."""
        return (k - 1) % self.Q if j == (k - 1) % self.Q else (k + 1) % self.Q

    # -- one kite crossing ------------------------------------------------
    def cross(self, k, j, x):
        """Enter kite `k` through edge `j` at coord `x`.

        Returns `('pole', k)` or `(j', x', transited)`."""
        if self.tie[k]:
            return self.other_edge(k, j), x, True
        if j == self.near[k]:
            x2 = x + self.delta[k] if self.osec[k] else x
            return self.far[k], x2, True
        if self.cmp(x, self.flo[k]) >= 0 and self.cmp(x, self.fhi[k]) <= 0:
            if self.cmp(x, self.ffix[k]) == 0:
                return "pole", k, False
            return j, self.fsum[k] - x, False
        x2 = x - self.delta[k] if self.osec[k] else x
        return self.near[k], x2, True

    # -- separatrix tracing -----------------------------------------------
    def trace(self, j, x, k, resume=None):
        """Follow the leaf from the interface point `(j, x)` into kite `k`.

        Returns (end_label, transits_of_kite0, interior_points_on_interface_m, steps,
        state).

        ⇒ RESUMABLE (s458, user-directed standing rule -- `CLAUDE.md` Measurement
        discipline).  `state` is `None` unless the walk hit `STEP_CAP`, in which case it
        is the JSON-able walk state to hand straight back as `resume=` at a HIGHER cap.
        Before s458 the cap DISCARDED `(j, x, k)`, so `3e6 -> 3e7` restarted every capped
        separatrix from step `0` and an escalation ladder cost `Sigma caps` rather than the
        final cap.  The state is ~`phi(4Q)` ints against `3e7` steps of compute.
        ⚠ `state` is APPENDED to the tuple, never inserted: pre-s458 callers index it
        positionally (`s451_quotient_path.py` takes `[0]`)."""
        if resume is None:
            t0 = 0
            onm = 1 if (j == self.m_idx) else 0
            steps = 0
        else:
            # ⚠ `x.f` is RECOMPUTED from the exact vector rather than restored from JSON:
            # it is a DERIVED float that the tracer already refreshes every 512 steps, so
            # rebuilding it is strictly better than round-tripping it.
            j, k = resume["j"], resume["k"]
            x = Coord(np.array(resume["x"], dtype=np.int64), 0.0)
            x.f = float(x.v @ self.ring.imvals)
            t0, onm, steps = resume["t0"], resume["onm"], resume["steps"]
        imv = self.ring.imvals
        while True:
            steps += 1
            if steps > STEP_CAP:
                # `steps - 1` is the number of crossings actually COMPLETED -- resuming
                # with it re-does the crossing this iteration skipped, and no other. The
                # REPORTED count keeps its pre-s458 value (`STEP_CAP + 1`) so that
                # `max_steps`/`tot_steps` and every ruling quoting them are unchanged.
                return ("CAP", t0, onm, steps,
                        {"j": j, "k": k, "x": x.v.tolist(),
                         "t0": t0, "onm": onm, "steps": steps - 1})
            if not steps & 511:
                # kill float DRIFT: recompute the float from the exact vector, and audit
                # the integer magnitude (int64 overflow would silently break equality)
                x.f = float(x.v @ imv)
                a = int(np.abs(x.v).max())
                if a > self.maxcoef:
                    self.maxcoef = a
            r = self.cross(k, j, x)
            if r[0] == "pole":
                return (("R", r[1]), t0, onm, steps, None)
            j2, x2, transited = r
            if transited and k == 0:
                t0 += 1
            # the corner is reached from INSIDE kite `k`, so `k` names the prong-end
            if -TOL < x2.f < TOL and x2.is_zero():
                return (("ZO", k), t0, onm, steps, None)
            if self.eq_exact(x2, self.ell[j2]):
                return (("ZA", k), t0, onm, steps, None)
            if j2 == self.m_idx:
                onm += 1
            k = self.next_kite(k, j2)
            j, x = j2, x2

    # -- the separatrix diagram -------------------------------------------
    def separatrices(self, prior=None):
        """All `Q` separatrices, with their `Fix(iota)` incidences.

        ⇒ RESUMABLE (s458).  Pass a list previously returned by this method as `prior` and
        ONLY its `CAP` entries are re-run -- each one CONTINUED from its stored walk state
        rather than restarted.  Everything already resolved is taken verbatim, so raising
        `STEP_CAP` costs the capped region and nothing else.  (Before s458 a cap raise
        re-traced every separatrix in the row: `14/31 eps0` re-walked all `58` when `6`
        were capped.)  `prior` must come from the SAME `(P, Q, eps)` and the same model
        kwargs -- `starts` is rebuilt here and the traced entries are positional."""
        Q = self.Q
        seps, starts = [], []
        # (v) the tie kite's horizontal Fix arm: R_0 -- Z_O (eps=0) / Z_A (eps=1)
        seps.append({"ends": [("R", 0), ("ZO" if self.eps == 0 else "ZA", 0)],
                     "in_fix": True, "t0": 0, "onm": 0, "steps": 0})
        # (vi) the horizontal interface, lone class only
        if self.ell[self.m_idx].is_zero():
            seps.append({"ends": [("ZO", None), ("ZA", None)], "in_fix": True,
                         "t0": 0, "onm": 0, "steps": 0})
        for k in range(Q):
            if self.tie[k]:
                continue
            lo, hi, fix = self.fold_lo_hi(k)
            far = self.far[k]
            if self.ell[self.near[k]].is_zero():
                # (vi) degenerate flank: fold covers the whole far edge; only the pole
                # prong starts here (the corner leaves ARE the horizontal interface).
                starts.append((far, fix, k, ("R", k)))
                continue
            starts.append((far, fix, k, ("R", k)))                      # pole prong
            starts.append((far, hi if self.osec[k] else lo, k,
                           ("ZO" if self.osec[k] else "ZA", k)))        # corner prong
        # `prior`'s traced entries are positional against `starts`; the `in_fix` prefix is
        # rebuilt above and never traced, so split it off before indexing.
        old = None
        if prior is not None:
            old = [s for s in prior if not s["in_fix"]]
            if len(old) != len(starts):
                raise ValueError(f"resume mismatch: {len(old)} stored traces vs "
                                 f"{len(starts)} starts -- wrong (P,Q,eps) or model kwargs")
        for i, (j, x, k, tag) in enumerate(starts):
            if old is not None and old[i]["ends"][1] != "CAP":
                seps.append(dict(old[i]))            # already resolved -- do NOT re-walk
                continue
            rs = old[i].get("state") if old is not None else None
            end, t0, onm, steps, state = self.trace(
                j, x, self.next_kite(k, j), resume=rs)
            rec = {"ends": [tag, end], "in_fix": False,
                   "t0": t0, "onm": onm, "steps": steps}
            if state is not None:
                rec["state"] = state                 # the cap's resume point
            seps.append(rec)
        return seps

    def readout(self, prior=None, with_seps=False):
        """⇒ `prior`/`with_seps` are the s458 resume hooks; with neither, byte-identical
        to the pre-s458 readout.  `with_seps` returns the per-separatrix table, which is
        what a store must persist to be able to resume (and is also the only way to
        decompose `tot_steps` after the fact).

        `prior` is a previous readout dict taken `with_seps=True`, at a LOWER `STEP_CAP`.
        ⚠⚠ THE FLOAT-GATE AUDIT IS MERGED, NOT RECOMPUTED, and that is not cosmetic: a
        resumed row only re-walks the capped separatrices, so `margin`/`maxcoef`/
        `exact_calls` gathered this pass cover a SUBSET of the row and would report a
        rosier `safety` than the same row measured from scratch ([OPS-223] quotes `safety`
        with every table).  Merging min/max/sum reproduces the from-scratch values
        exactly -- `_resume_equivalence()` is the test."""
        Q, eps = self.Q, self.eps
        seps = self.separatrices(prior=(prior["seps"] if prior is not None else None))
        # each traced separatrix is found from BOTH ends -> dedupe by halving
        traced = [s for s in seps if not s["in_fix"]]
        fixed = [s for s in seps if s["in_fix"]]
        n_edges = len(fixed) + len(traced) // 2
        h = len(fixed)
        m_vertical = not self.ell[self.m_idx].is_zero()
        # incidences, counted once per EDGE (each traced edge appears twice in `traced`)
        f_arc1 = sum(s["t0"] > 0 for s in traced) // 2
        f_arcm = (sum(s["onm"] > 0 for s in traced) // 2) if m_vertical else 0
        multi = sum(1 for s in traced if s["t0"] + (s["onm"] if m_vertical else 0) > 1)
        capped = sum(1 for s in traced if s["ends"][1] == "CAP")
        # H5 -- the pairing.  ⚠ THIS IS `p = 0` WRITTEN OUT, NOT EVIDENCE FOR IT
        # ([NCYL-292]): `p` counts iota-SWAPPED pairs, so `sigma(k) = -k` and `p = 0` are
        # the same statement.  Recorded because it also fixes the edge TYPES, and because
        # it is the form the residual is now stated in.
        sym = types = 0
        for s in traced:
            a_, b_ = s["ends"]
            if b_ == "CAP":
                continue
            types += (a_[0] == b_[0])
            sym += (a_[0] == b_[0] and b_[1] == (-a_[1]) % Q)
        f = f_arc1 + f_arcm
        v = 2 if m_vertical else 1
        arc1_side = "L1" if eps == 0 else "L2"
        n_sigma = {arc1_side: f_arc1 + 1}
        if m_vertical:
            n_sigma["H"] = f_arcm + 1
        if prior is not None:                       # see the docstring -- NOT cosmetic
            self.margin = min(self.margin, prior["margin"])
            self.maxcoef = max(self.maxcoef, prior["maxcoef"])
            # ⚠ `.get(…, 0)`: `exact_calls` is a COUNTER, so an older record that never
            # stored it contributes nothing and the merge stays sound (unlike `margin`/
            # `maxcoef`, which are bounds and are required fields).
            self.exact_calls += prior.get("exact_calls", 0)
        out = {"P": self.P, "Q": Q, "eps": eps, "lone": self.lone,
               "E": n_edges, "h": h, "v": v, "f": f, "n_sigma": n_sigma,
               "p2": Q - h - f, "multi": multi, "capped": capped,
               "sym": sym, "types": types, "n_traced": len(traced),
               "max_steps": max([s["steps"] for s in traced], default=0),
               "tot_steps": sum(s["steps"] for s in traced),
               "exact_calls": self.exact_calls, "margin": self.margin,
               "maxcoef": self.maxcoef, "safety": self.safety()}
        if with_seps:
            out["seps"] = seps
        return out


def model_row(P, Q, eps, prior=None, with_seps=False, **kw):
    return Chain(P, Q, eps, **kw).readout(prior=prior, with_seps=with_seps)


# ------------------------------------------------- s458: the resume-equivalence control

# (P, Q, eps, low_cap, high_cap) -- chosen so the LOW cap actually bites.  The set covers
# the three regimes a resume must get right, and the middle one is the easy miss:
#   * caps low, RESOLVES high        -- the resume must finish the walk correctly;
#   * caps at BOTH                   -- the state must survive a SECOND round trip;
#   * caps at neither                -- the resume must re-walk nothing at all.
_RESUME_CASES = [(4, 21, 0, 1_000_000, 5_000_000),
                 (10, 23, 0, 1_000_000, 5_000_000),
                 (7, 15, 1, 100_000, 500_000),
                 (3, 11, 0, 1_000_000, 5_000_000)]


def _resume_equivalence(cases=None, verbose=True):
    """⇒⇒ THE CONTROL FOR THE s458 RESUME, AND IT IS THE ONE THAT CAN FAIL.

    Escalating in two stages (`low` then resume at `high`) must give a readout IDENTICAL
    to measuring once at `high` -- every scalar, including `max_steps`, `tot_steps` and
    the `safety` audit.  If the resume dropped a crossing, resumed off by one, or lost
    the float-gate accumulators, this is where it shows.

    ⚠ THE PRIOR IS ROUND-TRIPPED THROUGH JSON FIRST, deliberately: that is the real path
    (the state lives in a store), and it is what turns the `("R", k)` end tuples into
    lists and the coordinate vector into a list of ints.  Testing the in-memory object
    would skip the only serialization the mechanism actually depends on."""
    global STEP_CAP
    keep = STEP_CAP
    rows, ok = [], True
    try:
        for (P, Q, eps, lo, hi) in (cases or _RESUME_CASES):
            STEP_CAP = hi                                   # (A) one shot at the high cap
            t0 = time.time()
            one = model_row(P, Q, eps, with_seps=True)
            t_one = time.time() - t0

            STEP_CAP = lo                                   # (B) low, then resume at high
            first = json.loads(json.dumps(model_row(P, Q, eps, with_seps=True)))
            STEP_CAP = hi
            t1 = time.time()
            two = model_row(P, Q, eps, prior=first, with_seps=True)
            t_two = time.time() - t1

            fields = [k for k in one if k not in ("seps", "n_sigma")]
            bad = [k for k in fields if one[k] != two[k]]
            if one["n_sigma"] != two["n_sigma"]:
                bad.append("n_sigma")
            # steps the resume pass did NOT have to re-walk: everything already resolved
            # at `lo`, plus the prefix each capped separatrix had already covered.
            saved = sum(s["steps"] for s in first["seps"]
                        if not s["in_fix"] and s["ends"][1] != "CAP")
            saved += sum(s["state"]["steps"] for s in first["seps"]
                         if not s["in_fix"] and s["ends"][1] == "CAP")
            n_cap = sum(1 for s in first["seps"]
                        if not s["in_fix"] and s["ends"][1] == "CAP")
            rows.append({"P": P, "Q": Q, "eps": eps, "lo": lo, "hi": hi,
                         "capped_at_lo": n_cap, "n_traced": one["n_traced"],
                         "tot_steps": one["tot_steps"], "steps_not_rewalked": saved,
                         "identical": not bad, "mismatched": bad,
                         "secs_oneshot": round(t_one, 1), "secs_resume": round(t_two, 1)})
            ok &= not bad
            if verbose:
                print(f"  {P}/{Q} eps{eps}  cap {lo}->{hi}  "
                      f"{n_cap}/{one['n_traced']} capped at lo  "
                      f"{'IDENTICAL' if not bad else 'MISMATCH ' + ','.join(bad)}  "
                      f"| saved {saved} steps  ({t_one:.1f}s one-shot vs "
                      f"{t_two:.1f}s resume)", flush=True)
    finally:
        STEP_CAP = keep
    if verbose:
        tot = sum(r["tot_steps"] for r in rows)
        sav = sum(r["steps_not_rewalked"] for r in rows)
        print(f"\n  {sum(r['identical'] for r in rows)}/{len(rows)} identical; "
              f"resume re-walked {tot - sav} of {tot} steps "
              f"({100 * (tot - sav) / max(tot, 1):.1f}%)")
        print("  ⚠ a case where the low cap does NOT bite contributes 0 saving and is "
              "still a real test\n    (it checks the resume re-walks NOTHING).")
    return {"ok": ok, "rows": rows}


# ---------------------------------------------------------------- H1/H2: validation

def validate(path="data/s446_covering_q30.json", order=None, sector=None,
             qmax=None, verbose=True, label="model"):
    ref = json.load(open(path))
    # ⚠ SCOPE: the necklace has a single iota-fixed kite only at ODD Q ([NCYL-291]);
    # even-Q rows are OUT OF SCOPE for this model and are never scored.
    rows = [r for r in ref["rows"] if r.get("status") == "ok" and r["Q"] % 2 == 1]
    if qmax:
        rows = [r for r in rows if r["Q"] <= qmax]
    tot = {"rows": 0, "nsig": 0, "f": 0, "h": 0, "v": 0, "p": 0, "E": 0,
           "multi": 0, "capped": 0}
    misses = []
    worst = 1.0
    for r in rows:
        try:
            mo = model_row(r["P"], r["Q"], r["eps"], order=order, sector=sector)
        except Exception as exc:                                   # noqa: BLE001
            misses.append((r["P"], r["Q"], r["eps"], f"EXC {exc}"))
            tot["rows"] += 1
            continue
        tot["rows"] += 1
        worst = min(worst, mo["margin"])
        ok_ns = mo["n_sigma"] == r["n_sigma"]
        tot["nsig"] += ok_ns
        tot["f"] += mo["f"] == r["f"]
        tot["h"] += mo["h"] == r["h"]
        tot["v"] += mo["v"] == r["v"]
        tot["p"] += mo["p2"] == r["n_nofold"]
        tot["E"] += mo["E"] == r["Q"]
        tot["multi"] += mo["multi"] == 0
        tot["capped"] += mo["capped"] == 0
        if not ok_ns and len(misses) < 12:
            misses.append((r["P"], r["Q"], r["eps"], "model", mo["n_sigma"],
                           "s446", r["n_sigma"]))
    if verbose:
        n = tot["rows"]
        print(f"  [{label}]  {n} traced rows from {path}")
        for key, name in [("nsig", "n_sigma per side == s446 TRACED"),
                          ("f", "f == s446 f"), ("h", "h"), ("v", "v"),
                          ("p", "2p == s446 n_nofold"), ("E", "E == Q"),
                          ("multi", "<=1 vertical hit per separatrix"),
                          ("capped", "no capped trace")]:
            print(f"      {name:38s}: {tot[key]}/{n}")
        print(f"      min float margin                      : {worst:.3e}")
        for x in misses[:12]:
            print("       miss", x)
    return tot, misses, worst


# ---------------------------------------------------------------- H4: coverage sweep

def sweep(qmax=49, qmin=31, out_path="data/s450_chain.json", budget_s=None,
          cap=3_000_000, verbose=True):
    """H4 coverage.  ⚠ THREE OUTCOMES, and the middle one is NOT a verdict:
    OK = every separatrix closed and `p == 0` (an EXACT certificate, since termination
    is decided on the ring vector); UNSCORABLE = some trace exceeded `cap` (the class may
    be non-CP, or the saddle connection may simply be longer than the cap -- `3/23` eps=0
    closes at `1740348` steps, past BOTH prior tracers' caps, so a cap is never evidence
    of non-CP); BAD = everything closed and `p != 0`, which would be the real event."""
    global STEP_CAP
    STEP_CAP = cap
    t0 = time.time()
    rows, bad, unscored = [], [], []
    n = 0
    stop = False
    for Q in range(qmin if qmin % 2 else qmin + 1, qmax + 1, 2):
        if stop:
            break
        for P in range(1, Q):
            if math.gcd(P, Q) != 1:
                continue
            if budget_s and time.time() - t0 > budget_s:
                print(f"  [budget stop inside Q={Q}]", flush=True)
                stop = True
                break
            for eps in (0, 1):
                mo = model_row(P, Q, eps)
                n += 1
                key = (mo["P"], mo["Q"], mo["eps"], mo["lone"])
                if mo["capped"]:
                    unscored.append(key)
                elif mo["p2"] or mo["multi"] or mo["E"] != Q:
                    bad.append(mo)
                    print("   BAD", {k: mo[k] for k in
                                     ("P", "Q", "eps", "E", "h", "f", "p2", "multi")},
                          flush=True)
                rows.append({k: mo[k] for k in ("P", "Q", "eps", "lone", "E", "h", "v",
                                                "f", "p2", "capped", "max_steps",
                                                "margin", "maxcoef", "safety")})
        if verbose:
            print(f"    Q<={Q} done  [{round(time.time()-t0,1)}s, {n} rows, "
                  f"{len(unscored)} unscorable, {len(bad)} BAD]", flush=True)
    scored = n - len(unscored)
    if verbose:
        print(f"\nSWEEP {qmin}<=Q<={qmax}: {n} class rows visited, {scored} SCORED "
              f"(all closed), {len(unscored)} UNSCORABLE (trace > {cap} steps), "
              f"{len(bad)} BAD  [{round(time.time()-t0,1)}s]")
        print(f"  p == 0 on {scored - len(bad)}/{scored} scored rows")
        sc = [r for r in rows if not r["capped"]]
        print(f"  worst float margin {min([r['margin'] for r in rows], default=1):.3e}"
              f" | max |coef| {max([r['maxcoef'] for r in rows], default=0)}")
        print(f"  ⇒ FLOAT-GATE AUDIT over the SCORED rows: worst margin "
              f"{min([r['margin'] for r in sc], default=1):.3e}, worst safety ratio "
              f"{min([r['safety'] for r in sc], default=0):.1f} (must be >> 1)")
        print(f"  UNSCORABLE: {unscored}")
    json.dump({"qmin": qmin, "qmax": qmax, "visited": n, "scored": scored,
               "cap": cap, "unscorable": unscored,
               "bad": [{k: b[k] for k in ('P', 'Q', 'eps', 'E', 'h', 'f', 'p2')}
                       for b in bad],
               "rows": rows}, open(out_path, "w"), indent=1)
    print(f"wrote {out_path}")
    return n, scored, unscored, bad


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "validate"
    kw = {}
    for a in sys.argv[2:]:
        for flag, key, cast in (("--qmax=", "qmax", int), ("--qmin=", "qmin", int),
                                ("--budget=", "budget_s", float), ("--cap=", "cap", int)):
            if a.startswith(flag):
                kw[key] = cast(a.split("=")[1])
    if mode == "validate":
        global STEP_CAP
        print("H1  MODEL vs s446's TRACED rows (odd Q only -- see the SCOPE note)")
        validate(qmax=kw.get("qmax"))
        # ⚠ the controls run at a SMALL qmax and a SMALL cap on purpose: a broken model
        # does not terminate, so every one of its traces runs to the cap.
        cq, STEP_CAP = kw.get("qmax", 15), kw.get("cap", 20_000)
        print(f"\nH2  CONTROLS at Q <= {cq} -- these MUST fail")
        validate(qmax=cq, label="TRUE MODEL (reference)")
        for o in (1, 2):
            validate(qmax=cq, order=o, label=f"gluing order = +{o}")
        for md in ("complement", "shift", "parity"):
            validate(qmax=cq, sector=md, label=f"sector={md}")
    elif mode == "sweep":
        sweep(**kw)
    elif mode == "resume":
        print("s458  RESUME EQUIVALENCE -- two-stage escalation vs one shot\n")
        r = _resume_equivalence()
        print("VERDICT:", "PASS" if r["ok"] else "FAIL")
        sys.exit(0 if r["ok"] else 1)
    else:
        print(__doc__)


if __name__ == "__main__":
    main()
