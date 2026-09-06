#!/usr/bin/env python3
"""
exact_unfold.py -- EXACT affine unfolding of a right-triangle billiard orbit over
`Q(zeta_4Q)`, for ANY `(P,Q)` and any direction on the `pi/(2Q)` grid.

WHAT IT IS.  `probes/s319_cut_exact.py` built this machine and hard-wired it to
`8/15` / `Q(zeta_60)` / the `pi/Q` direction grid (its `REFL` has `8 - m` literal
and its `cosm/sinm` are `pi/15`).  This is the same lever, generalised, and
sharing the field layer with `engine/cyclo_heights.CosBasis` instead of carrying
a private sympy one.

THE LEVER ([NCYL-087], s319's WHY).  Along a FIXED side word the direction is a
constant, so the position after `j` steps is AFFINE in the launch coordinate `s`:
`p_j = A_j + B_j*s` with `A_j, B_j` exact field elements.  The reflections
    L1 (x = L):  theta -> pi - theta        u -> 2Q - u
    L2 (y = 0):  theta -> -theta            u -> -u
    H  (axis alpha): theta -> 2*alpha-theta u -> 2P - u          (all mod 4Q)
keep `u` on the `pi/(2Q)` grid, and they preserve its PARITY -- which is the
direction-class parity of the necklace strand.

WHY THE REAL SUBFIELD IS ENOUGH, and why this is cheap.  `sin(u*pi/(2Q)) =
cos((Q-u)*pi/(2Q))`, so every coordinate that appears is a `c_j` or a product /
quotient of them.  Nothing complex is ever built, and `CosBasis` handles it all.

⇒⇒ THE PAYOFF THAT MOTIVATED IT: THE CIRCUMFERENCE NEEDS NO LAUNCH POINT.
Each step's ray parameter is `t_j = t0_j + t1_j*s`, and the direction vectors are
UNIT, so the orbit's arc length is `sum_j t_j = (sum_j t0_j) + (sum_j t1_j)*s`.
For a genuine CYLINDER word the length is constant across the cylinder, i.e.
`sum_j t1_j = 0` -- so the circumference is `sum_j t0_j`, an exact field element
depending only on the WORD and the side geometry.  The launch coordinate `s` is
a float cell midpoint and never enters.  `sum_j t1_j == 0` is therefore not an
assumption but a CERTIFICATE, and it can fail (it does, on a non-closing word).

GEOMETRY (`probes/s315_cylinder_count.geometry`, exactly -- not re-derived):
O = (0,0), R = (L,0), A = (L,1), `L = cot alpha`, `alpha = P*pi/(2Q)`;
L1: x = L, base R, tangent (0,1), length 1;  L2: y = 0, base O, tangent (1,0),
length L;  H: x = L*y, base O, tangent (cos alpha, sin alpha), length csc alpha.

⚠ TRAPS.
 (i)  The word must be a genuine FIRST-return / full-period word for the launch
      strip.  Nothing in the affine algebra notices an earlier return -- s319's
      trap, carried.  `closure()` checks the return map is the identity, which
      is the strip-level form of that check.
 (ii) `ring_cyc.RingContext` is the WRONG field (`Q(zeta_2Q)`, short by a factor
      of 2) -- [OPS-031], and `cyclo_heights` trap (iii).  Use `CosBasis(Q)`,
      which is `4Q`.
 (iii) A correct-looking VALUE is no evidence the field was big enough, since
      answers often land in a small subfield (s319).  The certificates here are
      exact-zero identities, not value agreements.
"""
import math
from fractions import Fraction

from cyclo_heights import CosBasis

SIDES = ("L1", "L2", "H")


# ------------------------------------------------------------------ vector helpers
def vadd(a, b):
    return [x + y for x, y in zip(a, b)]


def vsub(a, b):
    return [x - y for x, y in zip(a, b)]


def vscale(c, a):
    c = Fraction(c)
    return [c * x for x in a]


def vzero(D):
    return [Fraction(0)] * D


def is_zero(a):
    return not any(a)


class ExactUnfold:
    """Exact affine unfolding for the right triangle with `alpha = P*pi/(2Q)`."""

    def __init__(self, P, Q):
        self.P, self.Q, self.N = P, Q, 4 * Q
        self.B = B = CosBasis(Q)
        self.D = B.D
        self.cos_a = B.cos(P)              # cos(alpha)
        self.sin_a = B.cos(Q - P)          # sin(alpha) = cos((Q-P)pi/2Q)
        self.L = B.mul(self.cos_a, B.inv(self.sin_a))          # cot(alpha)
        self.ZERO, self.ONE = vzero(B.D), B.cos(0)
        self._inv_cache = {}
        # side base point / unit tangent / length, exactly
        self.geo = {
            "L1": (self.ONE, (self.L, self.ZERO), (self.ZERO, self.ONE)),
            "L2": (self.L, (self.ZERO, self.ZERO), (self.ONE, self.ZERO)),
            "H": (B.inv(self.sin_a), (self.ZERO, self.ZERO),
                  (self.cos_a, self.sin_a)),
        }

    # -------------------------------------------------------------- primitives
    def dirvec(self, u):
        """Unit direction at grid index `u`: `(cos(u*pi/2Q), sin(u*pi/2Q))`."""
        return self.B.cos(u), self.B.cos(self.Q - u)

    def reflect(self, side, u):
        if side == "L1":
            return (2 * self.Q - u) % self.N
        if side == "L2":
            return (-u) % self.N
        if side == "H":
            return (2 * self.P - u) % self.N
        raise ValueError(side)

    def _inv(self, key, vec):
        v = self._inv_cache.get(key)
        if v is None:
            v = self.B.inv(vec)
            self._inv_cache[key] = v
        return v

    def _step_time(self, side, u, A, Bc):
        """Ray parameter to reach `side` from `A + B*s` in direction `u`, as the
        affine pair `(t0, t1)`.  Pure algebra -- no test that this side is the one
        actually hit first (that is the caller's word)."""
        B_, L = self.B, self.L
        dx, dy = self.dirvec(u)
        if side == "L1":                                   # x = L
            num0, num1 = vsub(L, A[0]), vscale(-1, Bc[0])
            den, dkey = dx, ("dx", u)
        elif side == "L2":                                 # y = 0
            num0, num1 = vscale(-1, A[1]), vscale(-1, Bc[1])
            den, dkey = dy, ("dy", u)
        else:                                              # H: x - L*y = 0
            num0 = vscale(-1, vsub(A[0], B_.mul(L, A[1])))
            num1 = vscale(-1, vsub(Bc[0], B_.mul(L, Bc[1])))
            den, dkey = vsub(dx, B_.mul(L, dy)), ("dh", u)
        iden = self._inv(dkey, den)
        return B_.mul(num0, iden), B_.mul(num1, iden)

    # -------------------------------------------------------------- the unfolding
    def unfold(self, side0, u0, word):
        """Force the side sequence `word` from the launch side `side0` at grid
        direction `u0`, launched at `base0 + s*tan0`.

        Returns a dict with the exact affine data:
          `circ` / `circ_s`  -- arc length `= circ + circ_s * s`
          `A`, `Bc`          -- terminal position `A + Bc*s` (pairs of field elts)
          `u`                -- terminal direction index
          `steps`            -- per-bounce `(side, t0, t1, A, Bc)` AFTER the move,
                                BEFORE the reflection
        """
        _, base, tan = self.geo[side0]
        A = [base[0], base[1]]
        Bc = [tan[0], tan[1]]
        u = u0 % self.N
        circ, circ_s = vzero(self.D), vzero(self.D)
        steps = []
        for S in word:
            t0, t1 = self._step_time(S, u, A, Bc)
            dx, dy = self.dirvec(u)
            A = [vadd(A[0], self.B.mul(t0, dx)), vadd(A[1], self.B.mul(t0, dy))]
            Bc = [vadd(Bc[0], self.B.mul(t1, dx)), vadd(Bc[1], self.B.mul(t1, dy))]
            circ, circ_s = vadd(circ, t0), vadd(circ_s, t1)
            steps.append((S, t0, t1, (A[0], A[1]), (Bc[0], Bc[1])))
            u = self.reflect(S, u)
        return {"circ": circ, "circ_s": circ_s, "A": A, "Bc": Bc, "u": u,
                "steps": steps}

    # -------------------------------------------------------------- certificates
    def closure(self, res, side0, u0):
        """Is the unfolded word a genuine CLOSED cylinder word?

        Three exact identities, each able to fail:
          `dir_ok`   the terminal direction is the launch direction;
          `ret_id`   the return map on the launch side is the IDENTITY, i.e. the
                     terminal position is `base0 + s*tan0` again -- checked as an
                     identity in `s`, not at a point;
          `len_const` the arc length does not depend on `s` (`circ_s == 0`).
        `len_const` is what licenses reading the circumference off `circ` alone.
        """
        _, base, tan = self.geo[side0]
        dir_ok = (res["u"] == u0 % self.N)
        ret_id = (is_zero(vsub(res["A"][0], base[0]))
                  and is_zero(vsub(res["A"][1], base[1]))
                  and is_zero(vsub(res["Bc"][0], tan[0]))
                  and is_zero(vsub(res["Bc"][1], tan[1])))
        len_const = is_zero(res["circ_s"])
        return {"dir_ok": dir_ok, "ret_id": ret_id, "len_const": len_const,
                "closed": dir_ok and ret_id and len_const}

    def value(self, vec):
        return self.B.value(vec)

    def cos_form(self, vec, tol=None):
        """`vec` in the `{c_j}` basis, as a dict `j -> Fraction`, or None if it is
        not in the real subfield's canonical basis span (which cannot happen for
        anything this module builds, but is checked rather than assumed)."""
        from cyclo_heights import cos_basis_form
        return cos_basis_form(self.B, vec)


def circumference(P, Q, side0, u0, word):
    """Convenience: the exact circumference of the cylinder carrying `word`, plus
    its certificate dict.  Returns `(circ_vec, certs, unfolder)`."""
    uf = ExactUnfold(P, Q)
    res = uf.unfold(side0, u0, word)
    return res["circ"], uf.closure(res, side0, u0), uf
