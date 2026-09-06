#!/usr/bin/env python3
"""
cyclo_heights.py -- EXACT cylinder heights of a `pi/(2Q)` direction class, as elements
of `Q(zeta_4Q)`, from the BANKED integer crossing matrix.  No tracing, no ring run, no
float anywhere in the arithmetic, and -- unlike the s476 route this replaces -- no PSLQ.

WHAT IT COMPUTES.  [NCYL-247]'s flux identity, one equation per strip `(sigma,u)`:

        Sum_c  k_c(sigma,u) * h_c  =  |sigma| * flux(sigma,u)                     (F)

`k_c` is the integer crossing multiplicity, banked per cylinder as `k_per_strip` in
`data/s439_exact_store/`.  Stacked over every strip this is `K h = b` with `K` an
INTEGER matrix of full column rank ([NCYL-330], `242/242`).  The right-hand side is
elementary trigonometry in `alpha = P*pi/(2Q)` -- from `s315_cylinder_count.geometry`,
with `theta_u = u*pi/(2Q)` and `c_j := cos(j*pi/(2Q))`:

        L1:  |sigma| = 1,       n = (-1, 0)         ->  b = -c_u
        L2:  |sigma| = cot a,   n = (0, 1)          ->  b = (c_P / s_P) * s_u
        H :  |sigma| = csc a,   n = (sin a, -cos a) ->  b = s_{P-u} / s_P

  ⇒⇒ **SOLVE FOR `y_c := sin(alpha) * h_c`, NOT FOR `h_c`, AND THE WHOLE COMPUTATION
     BECOMES INDEX ARITHMETIC.**  Multiplying (F) through by `s_P` clears the only
     denominator, and `sin(j pi/2Q) = c_{Q-j}` turns every sine into a cosine:

        L1:  Sum_c k_c y_c = -c_{Q-P} * c_u        L2:  ... = c_P * c_{Q-u}
        H :  Sum_c k_c y_c =  c_{Q-P+u}

     Every right-hand side is a single cosine or a PRODUCT OF TWO, and
     `c_a*c_b = (c_{a+b} + c_{a-b})/2` is index arithmetic -- so no polynomial
     multiplication is ever performed, and `K` stays INTEGER so the elimination is
     pure `Fraction` scalar work.  ⚠ Scaling by the fixed field element `sin(alpha)`
     is a `Q`-linear isomorphism of `Q(zeta_4Q)`, so `dim_Q span{h_c}` is UNCHANGED:
     rank questions may be answered on `y` directly, with no inversion.  Recovering
     `h_c` itself (for display) does need one inversion and is `heights()`.

WHY `zeta_4Q` AND NOT `zeta_2Q`.  `ring_cyc.RingContext` -- and hence [NCYL-288]'s
`probes/s448_width_rank.py` -- works in `Q(zeta_2Q)`, which is correct for the
`Sigma_L1` partition (`u = 2Q`, even) but CANNOT REPRESENT the other class: at `3/7`
the lone heights are combinations of `cos(pi/14)`, `cos(3pi/14)`, `cos(5pi/14)` --
ODD cosines, of degree 6, not in the degree-3 subfield.  The `4Q` grid is the whole
point of the oblique strand ([NCYL-240]), so the field has to be `Q(zeta_4Q)`.

WHAT THE OVER-DETERMINATION BUYS -- THE CERTIFICATE.  `C` unknowns against every strip
that meets the class (19 at `3/7` lone, hundreds at large `Q`).  `solve()` uses a
full-column-rank subsystem and then checks `K y == b` EXACTLY on every remaining row.
That residual is `0` or it is not; there is no tolerance to tune, and it is a genuine
consistency check between the COMBINATORIAL layer (`k_c`, crossing counts) and the
METRIC layer (`|sigma|*flux`, closed form).
  ⚠ It is NOT an independent check against the store's Kac `height` -- for that, compare
  `value()` against the record ([NCYL-330]'s caveat: the counting reading and the
  Kac-integrating reading are two readings of the SAME traced partition).

⚠⚠ TRAPS.
 (i)  **`parity` IS THE CLASS LABEL THAT MATTERS HERE, NOT `class_role`.**  Which class
      holds `Sigma_L1` (`u = 2Q`, always even) is the parity-0 one; `paired`/`lone`
      cross-cuts parity (in `data/s439_exact_store`: paired sits at parity 0 on 82 rows
      and parity 1 on 45).  [NCYL-288]'s scope is the `Sigma_L1` partition = parity 0.
 (ii) **THE FIELD OF `y` IS NOT THE FIELD OF `h`.**  `y = sin(alpha)*h` and `sin(alpha)`
      is itself a cosine of index `Q-P`, so the scaling mixes parity.  Ranks transfer;
      FIELD MEMBERSHIP DOES NOT.  Test membership with `in_half_field()`, which carries
      the `(-1)^(Q-P)` twist, or unscale first with `heights()`.
 (iii) **`_reduction_table(N)` IS CALLED WITH `N = 4Q`, NOT `2Q`.**  Passing `2Q` here
      silently gives the subfield and every odd-cosine height becomes unrepresentable
      (the solve then fails on the residual rather than returning a wrong answer, but
      the error message will point at the data, not at the basis).

Public API
    CosBasis(Q)                 .cos(j) .cos_prod(a,b) .value(v) .galois(v,k) .inv(v)
    strip_rhs(basis, side, u, P)              -> the scaled RHS of (F) for one strip
    solve(P, Q, cylinders)                    -> ExactHeights (y-form + certificate)
    ExactHeights.heights()                    -> h_c coordinate vectors (one inversion)
    ExactHeights.rank()                       -> dim_Q span{h_c}, exact
    ExactHeights.in_half_field()              -> is every h_c in Q(zeta_2Q)?
    qrank(rows)                               -> exact rank over Q of Fraction vectors
    cos_basis_form(basis, vec)                -> readable {c_j} coordinates (unique)
"""
import math
from fractions import Fraction

from ring_cyc import _reduction_table


# ------------------------------------------------------------------ the field
class CosBasis:
    """Exact arithmetic in `Q(zeta_N)`, `N = 4Q`, `zeta = exp(i*pi/(2Q))`, in the
    power basis `{1, x, ..., x^(D-1)}` mod `Phi_N`, `D = phi(N)`.  Vectors are lists
    of `Fraction` of length `D`.  `c_j := cos(j*pi/(2Q)) = (zeta^j + zeta^-j)/2`."""

    def __init__(self, Q):
        self.Q = Q
        self.N = N = 4 * Q
        self.red, self.D = _reduction_table(N)
        self._cos_cache = {}
        self._basis_js = None

    def cos(self, j):
        """`c_j` as a length-`D` Fraction vector.  Handles any integer `j`."""
        j %= self.N
        v = self._cos_cache.get(j)
        if v is None:
            a, b = self.red[j], self.red[(-j) % self.N]
            v = [Fraction(int(x) + int(y), 2) for x, y in zip(a, b)]
            self._cos_cache[j] = v
        return v

    def cos_prod(self, a, b):
        """`c_a * c_b = (c_{a+b} + c_{a-b}) / 2` -- no polynomial multiplication."""
        u, v = self.cos(a + b), self.cos(a - b)
        return [(x + y) / 2 for x, y in zip(u, v)]

    def value(self, vec):
        """Float value.  The power basis is `x^i = exp(2*pi*i*i/N)`; every element we
        build is real, so the imaginary parts cancel and the real part is the value."""
        return float(sum(float(c) * math.cos(2.0 * math.pi * i / self.N)
                         for i, c in enumerate(vec) if c))

    def value_mp(self, vec, dps=50):
        from mpmath import mp, mpf, cos, pi
        old, mp.dps = mp.dps, dps
        try:
            return sum(mpf(c.numerator) / c.denominator
                       * cos(2 * pi * i / self.N)
                       for i, c in enumerate(vec) if c)
        finally:
            mp.dps = old

    def galois(self, vec, k):
        """Apply `zeta -> zeta^k` (`gcd(k,N) = 1`).  `x^i` has coordinate vector
        `red[i*k mod N]`."""
        assert math.gcd(k, self.N) == 1, f"{k} not a unit mod {self.N}"
        out = [Fraction(0)] * self.D
        for i, c in enumerate(vec):
            if not c:
                continue
            for d, r in enumerate(self.red[(i * k) % self.N]):
                if r:
                    out[d] += c * int(r)
        return out

    def mul(self, a, b):
        """General product mod `Phi_N`.  Only needed by `inv`; the solve never calls it."""
        conv = [Fraction(0)] * (2 * self.D - 1)
        for i, x in enumerate(a):
            if not x:
                continue
            for j, y in enumerate(b):
                if y:
                    conv[i + j] += x * y
        out = [Fraction(0)] * self.D
        for m, c in enumerate(conv):
            if not c:
                continue
            for d, r in enumerate(self.red[m] if m < self.N else self.red[m % self.N]):
                if r:
                    out[d] += c * int(r)
        return out

    def inv(self, vec):
        """Multiplicative inverse mod `Phi_N`, by extended Euclid over `Q[x]`."""
        phi = self._phi_coeffs()
        r0, r1 = phi, _trim(list(vec))
        s0, s1 = [Fraction(0)], [Fraction(1)]
        while r1:
            q, r = _divmod_poly(r0, r1)
            r0, r1 = r1, r
            s0, s1 = s1, _trim(_sub_poly(s0, _mul_poly(q, s1)))
        assert len(r0) == 1 and r0[0] != 0, "not invertible (element is zero)"
        out = [x / r0[0] for x in s0] + [Fraction(0)] * self.D
        return _reduce_to(out, self)

    def _phi_coeffs(self):
        """`Phi_N` as ascending Fraction coefficients, read off the reduction table:
        `x^D = -(red[D] as a poly)` since `red[D] = x^D mod Phi_N`."""
        if getattr(self, "_phi", None) is None:
            self._phi = _trim([Fraction(-int(c)) for c in self.red[self.D]]
                              + [Fraction(1)])
        return self._phi

    def basis_indices(self):
        """A canonical `Q`-basis of the REAL subfield as cosines: `{c_j}` for the
        smallest `j = 0, 1, 2, ...` that keep increasing the rank, `phi(N)/2` of them."""
        if self._basis_js is None:
            js, rows = [], []
            for j in range(0, self.N):
                cand = rows + [self.cos(j)]
                if qrank(cand) > len(rows):
                    rows, js = cand, js + [j]
                    if len(js) == self.D // 2:
                        break
            self._basis_js = js
        return self._basis_js


def _trim(p):
    while p and p[-1] == 0:
        p.pop()
    return p


def _sub_poly(a, b):
    n = max(len(a), len(b))
    return [(a[i] if i < len(a) else Fraction(0))
            - (b[i] if i < len(b) else Fraction(0)) for i in range(n)]


def _mul_poly(a, b):
    if not a or not b:
        return []
    out = [Fraction(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                if y:
                    out[i + j] += x * y
    return out


def _divmod_poly(a, b):
    a, b = list(a), _trim(list(b))
    q = [Fraction(0)] * max(1, len(a) - len(b) + 1)
    while len(a) >= len(b) and _trim(a):
        d = len(a) - len(b)
        f = a[-1] / b[-1]
        q[d] = f
        for i, x in enumerate(b):
            a[i + d] -= f * x
        _trim(a)
    return _trim(q), a


def _reduce_to(coeffs, basis):
    out = [Fraction(0)] * basis.D
    for m, c in enumerate(coeffs):
        if not c:
            continue
        for d, r in enumerate(basis.red[m % basis.N]):
            if r:
                out[d] += c * int(r)
    return out


# ------------------------------------------------------------------ linear algebra
def qrank(rows):
    """Exact rank over `Q` of a list of Fraction vectors.  Same routine shape as
    `s448_width_rank._rank`, which [NCYL-288]'s H3 control validated against
    pseudo-random vectors of the same shape."""
    M = [list(r) for r in rows]
    if not M:
        return 0
    n, D, r = len(M), len(M[0]), 0
    for col in range(D):
        piv = next((i for i in range(r, n) if M[i][col] != 0), None)
        if piv is None:
            continue
        M[r], M[piv] = M[piv], M[r]
        pv = M[r][col]
        M[r] = [x / pv for x in M[r]]
        for i in range(n):
            if i != r and M[i][col] != 0:
                f = M[i][col]
                M[i] = [x - f * y for x, y in zip(M[i], M[r])]
        r += 1
        if r == n:
            break
    return r


# ------------------------------------------------------------------ the RHS of (F)
def strip_rhs(basis, side, u, P):
    """`sin(alpha) * |sigma| * flux(sigma,u)` for one strip, exactly.

    Derived from `s315_cylinder_count.geometry` (`|L1| = 1`, `|L2| = cot a`,
    `|H| = csc a`; inward normals `(-1,0)`, `(0,1)`, `(sin a, -cos a)`) with
    `alpha = P*pi/(2Q)` and `sin(j*pi/2Q) = c_{Q-j}`.  Matches the closed forms
    [OPS-259] quotes, multiplied through by `sin(alpha) = c_{Q-P}`."""
    Q = basis.Q
    if side == "L1":
        return [-x for x in basis.cos_prod(Q - P, u)]
    if side == "L2":
        return basis.cos_prod(P, Q - u)
    if side == "H":
        return basis.cos(Q - P + u)
    raise ValueError(f"unknown side {side!r}")


# ------------------------------------------------------------------ the solve
class ExactHeights:
    """Result of `solve`.  `Y[c]` is `sin(alpha)*h_c` as a Fraction vector."""

    def __init__(self, P, Q, basis, Y, strips, n_rows, residual_bad, rank_K):
        self.P, self.Q, self.basis = P, Q, basis
        self.Y, self.strips = Y, strips
        self.C, self.D = len(Y), basis.D
        self.n_rows, self.residual_bad, self.rank_K = n_rows, residual_bad, rank_K
        self._H = None

    # -- the certificate -----------------------------------------------------
    @property
    def certified(self):
        """Every strip equation -- including the ones the solve did not use -- is
        satisfied EXACTLY, and `K` had full column rank."""
        return self.residual_bad == 0 and self.rank_K == self.C

    def values(self):
        """Float `h_c`, for comparison against a stored Kac height."""
        s = math.sin(self.P * math.pi / (2 * self.Q))
        return [self.basis.value(y) / s for y in self.Y]

    # -- the objects ---------------------------------------------------------
    def heights(self):
        """`h_c` coordinate vectors -- one field inversion (`1/sin alpha`).  Not
        needed for `rank()`; needed for `cos_basis_form` and for field membership
        read directly rather than through the twist."""
        if self._H is None:
            inv = self.basis.inv(self.basis.cos(self.Q - self.P))
            self._H = [self.basis.mul(y, inv) for y in self.Y]
        return self._H

    def rank(self):
        """`dim_Q span{h_c}`, exact.  Computed on `Y`: scaling by the fixed field
        element `sin(alpha)` is a `Q`-linear isomorphism, so the rank is the same."""
        return qrank(self.Y)

    def rank_on(self, mask):
        """`dim_Q span{h_c : mask[c]}` -- e.g. the `Sigma_L1`-met sub-collection,
        which is [NCYL-288]'s object."""
        rows = [y for y, m in zip(self.Y, mask) if m]
        return qrank(rows)

    def eigen_sign(self):
        """`+1` / `-1` if EVERY `h_c` is an eigenvector of the Galois element
        `sigma : zeta -> zeta^(1+2Q)` with that eigenvalue; `None` otherwise.

        `sigma` fixes `zeta^2 = zeta_2Q` and sends `c_j -> (-1)^j c_j`, so its `+1`
        eigenspace on the real field is `Q(zeta_2Q)+` (the even cosines) and its `-1`
        eigenspace is the odd-cosine complement.  BOTH have dimension `phi(2Q)/2`,
        which is what caps the height span -- so being an eigenvector is the load-
        bearing property, and merely *not being fixed* is not enough.

        ⚠ Trap (ii): tested on `Y = sin(alpha)*h`.  `sigma(c_{Q-P}) = (-1)^(Q-P)
        c_{Q-P}`, so `sigma(h) = e*h` iff `sigma(y) = (-1)^(Q-P) * e * y`."""
        k = 1 + 2 * self.Q
        tw = 1 if (self.Q - self.P) % 2 == 0 else -1
        out = None
        for y in self.Y:
            g = self.basis.galois(y, k)
            if all(a == b for a, b in zip(g, y)):
                e = tw
            elif all(a == -b for a, b in zip(g, y)):
                e = -tw
            else:
                return None
            if out is None:
                out = e
            elif out != e:
                return None
        return out

    def in_half_field(self):
        """Is every `h_c` in `Q(zeta_2Q)+`?  I.e. `eigen_sign() == +1`."""
        return self.eigen_sign() == 1


def solve(P, Q, cylinders, basis=None):
    """Exact heights of one direction class.

    `cylinders` is the store's per-class list; each entry needs `k_per_strip`, a
    dict `"SIDE:u" -> k`.  Returns an `ExactHeights`; check `.certified`."""
    basis = basis or CosBasis(Q)
    C = len(cylinders)
    if C == 0:
        raise ValueError("no cylinders")
    strips = sorted({s for c in cylinders for s in c["k_per_strip"]})
    K = [[Fraction(c["k_per_strip"].get(s, 0)) for c in cylinders] for s in strips]
    B = []
    for s in strips:
        side, us = s.split(":")
        B.append(strip_rhs(basis, side, int(us), P))

    M = [K[i] + B[i] for i in range(len(strips))]
    r = 0
    for col in range(C):
        piv = next((i for i in range(r, len(M)) if M[i][col] != 0), None)
        if piv is None:
            continue
        M[r], M[piv] = M[piv], M[r]
        pv = M[r][col]
        M[r] = [x / pv for x in M[r]]
        for i in range(len(M)):
            if i != r and M[i][col] != 0:
                f = M[i][col]
                M[i] = [x - f * y for x, y in zip(M[i], M[r])]
        r += 1
    rank_K = r
    if r < C:
        return ExactHeights(P, Q, basis, [[Fraction(0)] * basis.D] * C,
                            strips, len(strips), -1, rank_K)
    Y = [M[i][C:] for i in range(C)]

    bad = 0
    for i in range(len(strips)):
        row = K[i]
        for d in range(basis.D):
            if sum(row[c] * Y[c][d] for c in range(C) if row[c]) != B[i][d]:
                bad += 1
                break
    return ExactHeights(P, Q, basis, Y, strips, len(strips), bad, rank_K)


def solve_record(record, role, basis=None):
    """`solve` straight off a `data/s439_exact_store/<P>_<Q>.json` record."""
    cl = record["classes"][role]
    return solve(record["P"], record["Q"], cl["cylinders"], basis=basis)


# ------------------------------------------------------------------ display
def cos_basis_form(basis, vec):
    """`vec` as `{c_j}` coordinates in `basis.basis_indices()` -- the canonical
    smallest-index `Q`-basis of the real subfield, so the representation is UNIQUE.
    Returns `{j: Fraction}` for the nonzero coordinates, or `None` if `vec` is not
    real (which would be a bug, not a datum)."""
    js = basis.basis_indices()
    cols = [basis.cos(j) for j in js]
    n = len(js)
    M = [[cols[k][d] for k in range(n)] + [vec[d]] for d in range(basis.D)]
    r, where = 0, {}
    for col in range(n):
        piv = next((i for i in range(r, len(M)) if M[i][col] != 0), None)
        if piv is None:
            continue
        M[r], M[piv] = M[piv], M[r]
        pv = M[r][col]
        M[r] = [x / pv for x in M[r]]
        for i in range(len(M)):
            if i != r and M[i][col] != 0:
                f = M[i][col]
                M[i] = [x - f * y for x, y in zip(M[i], M[r])]
        where[col] = r
        r += 1
    for i in range(r, len(M)):
        if M[i][n] != 0:
            return None
    return {js[col]: M[row][n] for col, row in where.items() if M[row][n] != 0}


def fmt_cos(form):
    if not form:
        return "0"
    out = []
    for j, a in sorted(form.items()):
        s = "-" if a < 0 else "+"
        m = abs(a)
        out.append(f" {s} " + (f"c{j}" if m == 1 else f"{m}*c{j}"))
    return "".join(out).strip().lstrip("+ ")
