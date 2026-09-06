#!/usr/bin/env python3
"""ring_cyc.py — fixed-angle exact arithmetic in the cyclotomic ring Z[ζ_N].

For α = (P/Q)·(π/2), every billiard coordinate is a CosPoly value
    V = (1/2^s) Σ_k c_k cos(2kα),
and cos(2kα)=cos(kPπ/Q)=(ζ^{kP}+ζ^{-kP})/2 with ζ = ζ_N = exp(iπ/Q), N=2Q.
So 2^{s+1}·V ∈ Z[ζ_N] = Z[x]/Φ_N(x): an INTEGER vector of length D=φ(N) in the
integral power basis {1,ζ,…,ζ^{D-1}} (a Z-basis of the ring of integers).

Why this exists
---------------
A CosPoly's Chebyshev DEGREE grows ~linearly with reflection depth, so its
coefficients blow up (O(depth²) work per reflect; int64 overflow by ~depth 200
in the redundant {cos mπ/Q} spanning set). This representation is CANONICAL and
BOUNDED: the integer entries track the true coordinate magnitude (~linear in
depth), every reflect step is O(D²) regardless of depth, and float evaluation
from the reduced vector is numerically STABLE (no growing cancellation). This is
the depth-independent fast path the golden-orphan measurement needs (§81, OPEN #2).

RingElem = (integer vector length D, power-of-2 shift); value = ring(vec)/2^shift.
Ops: +, -, * (ring or int), /2 (shift), and EXACT /(1±cos2α) via a precomputed
ring inverse. Division is exact because the true quotient is an algebraic integer
and {ζ^i} is a Z-basis, so its coordinates are integers.
"""

import math
from fractions import Fraction

import numpy as np

from n_exact import _reduction_table  # reuse the vetted cyclotomic table


def _rational_solve(M, b):
    """Solve M x = b exactly over Q. M: D×D int rows, b: length-D int list.
    Returns list of Fraction. Assumes M invertible (raises on singular)."""
    D = len(M)
    A = [[Fraction(M[i][j]) for j in range(D)] + [Fraction(b[i])] for i in range(D)]
    for col in range(D):
        piv = next((r for r in range(col, D) if A[r][col] != 0), None)
        if piv is None:
            raise ValueError("singular matrix in ring inverse")
        A[col], A[piv] = A[piv], A[col]
        pv = A[col][col]
        A[col] = [v / pv for v in A[col]]
        for r in range(D):
            if r != col and A[r][col] != 0:
                f = A[r][col]
                A[r] = [A[r][j] - f * A[col][j] for j in range(D + 1)]
    return [A[i][D] for i in range(D)]


# Three primes just below 2^31 — products of residues fit in int64, so the
# mod-p elimination stays in fast numpy int64; CRT of the three lifts to a
# ~93-bit modulus, ample for reconstructing the (tiny) ring-inverse rationals.
_SOLVE_PRIMES = (2147483647, 2147483477, 2147483423)


def _solve_mod_p(M, b, p):
    """Solve M x ≡ b (mod p) by Gauss-Jordan in numpy int64. M: D×D int64,
    b: length-D int64. Returns x mod p (int64 array) or None if singular mod p."""
    D = M.shape[0]
    A = M % p
    rhs = (b % p).astype(np.int64)
    for k in range(D):
        piv = k
        while piv < D and A[piv, k] == 0:
            piv += 1
        if piv == D:
            return None
        if piv != k:
            A[[k, piv]] = A[[piv, k]]
            rhs[k], rhs[piv] = rhs[piv].copy(), rhs[k].copy()
        try:
            inv = pow(int(A[k, k]), -1, p)
        except ValueError:
            return None
        A[k] = (A[k] * inv) % p
        rhs[k] = (rhs[k] * inv) % p
        col = A[:, k].copy()
        col[k] = 0
        A = (A - np.outer(col, A[k])) % p
        rhs = (rhs - col * rhs[k]) % p
    return rhs


def _rat_recon(a, m):
    """Rational reconstruction: (n, d), gcd(n,d)=1, d>0, n ≡ a·d (mod m),
    |n|,d ≤ √(m/2). Returns (n, d) or None."""
    a %= m
    bound = math.isqrt(m // 2)
    r0, r1 = m, a
    s0, s1 = 0, 1
    while r1 > bound:
        q = r0 // r1
        r0, r1 = r1, r0 - q * r1
        s0, s1 = s1, s0 - q * s1
    n, d = r1, s1
    if d < 0:
        n, d = -n, -d
    if d == 0 or d > bound or math.gcd(abs(n), d) != 1:
        return None
    return n, d


def _crt_solve(M, b, primes=_SOLVE_PRIMES):
    """Solve M x = b exactly over Q via multi-prime residues + CRT + rational
    reconstruction. M, b: numpy int64. Returns list of Fraction, or None if the
    system is singular mod a prime or reconstruction fails (caller falls back)."""
    residues, mods = [], []
    for p in primes:
        xp = _solve_mod_p(M, b, p)
        if xp is None:
            return None
        residues.append(xp.astype(object))
        mods.append(p)
    # CRT the residues coordinate-wise to x mod (∏ p)
    Mprod = 1
    for p in mods:
        Mprod *= p
    x = np.zeros(len(b), dtype=object)
    for p, r in zip(mods, residues):
        Mp = Mprod // p
        coeff = (Mp * pow(Mp % p, -1, p)) % Mprod
        x = (x + r * coeff) % Mprod
    out = []
    for xi in x.tolist():
        rec = _rat_recon(int(xi), Mprod)
        if rec is None:
            return None
        out.append(Fraction(rec[0], rec[1]))
    return out


class RingContext:
    """Precomputed fixed-angle context for α=(P/Q)(π/2): reduction table,
    ring constants (1±cos2α), and their exact ring inverses."""

    def __init__(self, P, Q):
        g = math.gcd(P, Q)
        P, Q = P // g, Q // g
        if not (0 < P < Q):
            raise ValueError("need 0 < P < Q")
        self.P, self.Q = P, Q
        self.N = N = 2 * Q
        self.red, self.D = _reduction_table(N)
        D = self.D
        self.alpha = (P / Q) * (math.pi / 2)
        # float basis {cos(iπ/Q)} for stable evaluation (Re of ζ^i)
        self._cosb = [math.cos(math.pi * i / Q) for i in range(D)]
        self._cosb_np = np.array(self._cosb, dtype=np.float64)
        # reduction matrix: row m = x^m mod Φ_N, for m in [0, 2D-1) (2D-2 ≤ N-1
        # for even N=2Q, so red has these rows). conv(a,b)·REDmat reduces a*b.
        self._REDmat = np.array([self.red[m] for m in range(2 * D - 1)],
                                dtype=np.int64)
        # The top D rows are the identity (x^m mod Φ = x^m for m<D), so only the
        # HIGH half reduces: reduce(conv) = conv[:D] + conv[D:] @ REDhi. Φ_N is
        # SPARSE ⇒ REDhi is ~0.3% dense at large D (1317/434940 at D=660), so the
        # reduce is done as a bincount over its nonzeros in float64 (numpy has no
        # int64 BLAS path) rather than a dense matmul — ~2.5× faster at D=660,
        # bit-for-bit identical. Exact for the tiny coeffs these orbits produce
        # (max_bits ~5-9, §82); the _mul guard enforces the 2^52 bound.
        _rr, _cc = np.nonzero(self._REDmat[D:])
        self._red_rows = _rr                       # index into conv[D:]
        self._red_cols = _cc                       # index into the output coord
        self._red_vals = self._REDmat[D:][_rr, _cc].astype(np.float64)
        # per-column |·| sum of the high block → bound on the reduction's partial
        # sums (+1 for the pass-through identity part) for the float-exactness guard
        self._red_colabs = 1 + int(np.abs(self._REDmat[D:]).sum(axis=0).max(initial=0))

        # 2(1±cos2α) as integer vectors: 2·e0 ± (ζ^P + ζ^{-P})
        two = 2 * np.asarray(self.red[0], dtype=np.int64)
        cP = np.asarray(self.red[P % N], dtype=np.int64)
        cN = np.asarray(self.red[(N - P) % N], dtype=np.int64)
        self.g_plus = two + cP + cN
        self.g_minus = two - cP - cN

        # exact ring inverses of the integer elements g_plus, g_minus
        self.invp = self._ring_inverse(self.g_plus)   # (V, c): g_plus^{-1}=V/c
        self.invm = self._ring_inverse(self.g_minus)

        # (1±cos2α) themselves as RingElems: g/2  → vector g, shift 1
        self.one_plus = RingElem(self.g_plus.copy(), 1, self)
        self.one_minus = RingElem(self.g_minus.copy(), 1, self)

    # -- low-level integer-vector ring ops (length-D numpy int64) -------------
    def _mul(self, a, b):
        """Multiply two integer ring vectors, reduce mod Φ_N. O(D²) via numpy.

        Inputs and output are numpy int64 arrays (length D). When coefficients
        are small enough to stay integer-exact in float64 (the usual case —
        max_bits ~5-9 over any orbit, §82) the polynomial product and mod-Φ
        reduction are done in float64 so `convolve` hits SIMD and the sparse
        reduce uses bincount (numpy has no fast int64 path); the result is
        rounded back. If the magnitudes could exceed the 2^52 integer-exact
        bound, it falls back to the exact int64 path (guarded at 2^60)."""
        am = int(np.abs(a).max(initial=0))
        bm = int(np.abs(b).max(initial=0))
        # convolve accumulates ≤ D products each ≤ am·bm; float is exact iff <2^52
        if am * bm * self.D < (1 << 52):
            conv = np.convolve(a.astype(np.float64), b.astype(np.float64))
            # reduction partial sums ≤ maxc·(Σ|REDhi col|+1); exact iff <2^52
            if float(np.abs(conv).max(initial=0)) * self._red_colabs < (1 << 52):
                hi = conv[self.D:]
                red = np.bincount(self._red_cols, weights=hi[self._red_rows] * self._red_vals,
                                  minlength=self.D)
                return np.rint(conv[:self.D] + red).astype(np.int64)
        # exact int64 fallback for rare large-coefficient elements
        res = np.convolve(a, b) @ self._REDmat
        assert int(np.abs(res).max(initial=0)) < (1 << 60), "int64 overflow risk"
        return res

    def _ring_inverse(self, g):
        """Return (V, c): int64 vector V and int c with g·(V/c)=1 in the ring.

        Solves the D×D system g·x = 1 for the inverse coordinates. The default
        path is a modular solve + CRT + rational reconstruction (O(D³) but in
        fast numpy int64 mod-p arithmetic — the ring inverse has tiny coords, so
        one 31-bit prime usually reconstructs it), ~40× faster than exact
        Fraction Gaussian elimination at D≈660. It escalates 1 prime → 3 primes
        → exact Fraction solve; each candidate is accepted only if it passes the
        g·V == c·1 verification, so an under-resolved modular solve escalates
        rather than corrupting the result."""
        D = self.D
        cols = []
        for j in range(D):
            ej = np.zeros(D, dtype=np.int64)
            ej[j] = 1
            cols.append(self._mul(g, ej))
        # M[i][j] = (g·x^j)_i = column j of the multiply-by-g matrix
        M = np.array(cols, dtype=np.int64).T
        e0 = np.zeros(D, dtype=np.int64)
        e0[0] = 1

        for solve in (lambda: _crt_solve(M, e0, _SOLVE_PRIMES[:1]),
                      lambda: _crt_solve(M, e0, _SOLVE_PRIMES),
                      lambda: _rational_solve(M.tolist(), e0.tolist())):
            x = solve()
            if x is None:
                continue
            c = 1
            for v in x:
                c = c * v.denominator // math.gcd(c, v.denominator)
            V = np.array([int(v * c) for v in x], dtype=np.int64)
            chk = self._mul(g, V)
            if int(chk[0]) == c and not np.any(chk[1:]):
                return (V, c)
        raise AssertionError("ring inverse verification failed")

    def from_cospoly(self, poly):
        """Convert a CosPoly (value at this α) into a reduced RingElem."""
        N, D, P, red = self.N, self.D, self.P, self.red
        W = np.zeros(D, dtype=np.int64)
        for k, ck in enumerate(poly.coeffs):
            ck = int(ck)
            if ck == 0:
                continue
            m = (k * P) % N
            rm = np.asarray(red[m], dtype=np.int64)
            rn = np.asarray(red[(N - m) % N], dtype=np.int64)
            W += ck * (rm + rn)
        e = RingElem(W, poly.shift + 1, self)
        e.reduce()
        return e

    def zero(self):
        return RingElem(np.zeros(self.D, dtype=np.int64), 0, self)


class RingElem:
    """value = ring(vec)/2^shift, vec an integer vector length D in Z[ζ_N]."""
    __slots__ = ('vec', 'shift', 'ctx')

    def __init__(self, vec, shift, ctx):
        # vec is a numpy int64 array (length D); normalize only if a caller
        # hands in a list. Internal ops already produce int64 arrays, so the
        # common path is a no-op.
        self.vec = vec if type(vec) is np.ndarray else np.asarray(vec, dtype=np.int64)
        self.shift = shift
        self.ctx = ctx

    def copy(self):
        return RingElem(self.vec.copy(), self.shift, self.ctx)

    def reduce(self):
        """Cancel common powers of 2 between vec and shift (bounds the shift)."""
        if self.shift == 0:
            return self
        vec = self.vec
        nz = vec[vec != 0]
        if nz.size == 0:
            self.shift = 0
            return self
        # trailing-zero count of the least-divisible entry: a & -a isolates the
        # lowest set bit (a power of 2); the min over entries is 2**(min tz).
        a = np.abs(nz)
        min_low = int((a & (-a)).min())
        tz = min_low.bit_length() - 1
        mtz = tz if tz < self.shift else self.shift
        if mtz > 0:
            self.vec = vec >> mtz   # exact: every entry divisible by 2**mtz
            self.shift -= mtz
        return self

    def _align(self, o):
        if self.shift == o.shift:
            return self.vec, o.vec, self.shift
        if self.shift > o.shift:
            d = self.shift - o.shift
            return self.vec, o.vec << d, self.shift
        d = o.shift - self.shift
        return self.vec << d, o.vec, o.shift

    def __add__(self, o):
        a, b, s = self._align(o)
        return RingElem(a + b, s, self.ctx)

    def __sub__(self, o):
        a, b, s = self._align(o)
        return RingElem(a - b, s, self.ctx)

    def __mul__(self, o):
        if isinstance(o, int):
            return RingElem(self.vec * o, self.shift, self.ctx)
        return RingElem(self.ctx._mul(self.vec, o.vec),
                        self.shift + o.shift, self.ctx)

    __rmul__ = __mul__

    def div2(self):
        return RingElem(self.vec.copy(), self.shift + 1, self.ctx)

    def _div_by(self, inv):
        V, c = inv
        prod = self.ctx._mul(self.vec, V)
        p2 = 2 * prod
        assert not np.any(p2 % c), "inexact division by (1±cos2α)"
        r = RingElem(p2 // c, self.shift, self.ctx)   # floor == exact here
        r.reduce()
        return r

    def div_one_plus(self):
        return self._div_by(self.ctx.invp)

    def div_one_minus(self):
        return self._div_by(self.ctx.invm)

    def is_zero(self):
        return not self.vec.any()

    def to_float(self):
        return float(self.vec @ self.ctx._cosb_np) / (2 ** self.shift)

    def to_mpf(self, dps=50):
        """Exact value to arbitrary precision: Σ vec[i]·cos(iπ/Q) / 2^shift."""
        from mpmath import mp, mpf, cos, pi
        mp.dps = dps
        Q = self.ctx.Q
        tot = mpf(0)
        for i, c in enumerate(self.vec):
            if c:
                tot += int(c) * cos(mpf(i) * pi / Q)
        return tot / (mpf(2) ** int(self.shift))

    def max_bits(self):
        return int(np.abs(self.vec).max(initial=0)).bit_length()


# ---------------------------------------------------------------------------
# Self-test: cross-check ring arithmetic against CosPoly.evaluate at the angle.
# ---------------------------------------------------------------------------
def _selftest():
    import random
    from cos_poly import CosPoly

    random.seed(1)
    cases = [(1, 5), (2, 7), (3, 7), (2, 11), (3, 11), (5, 13), (3, 21), (5, 34)]
    worst = 0.0
    for P, Q in cases:
        ctx = RingContext(P, Q)
        alpha = ctx.alpha

        def rand_poly():
            deg = random.randint(0, 6)
            coeffs = [random.randint(-5, 5) for _ in range(deg + 1)]
            if not any(coeffs):
                coeffs[0] = 1
            return CosPoly(coeffs, random.randint(0, 3))

        for _ in range(40):
            f = rand_poly()
            h = rand_poly()
            rf = ctx.from_cospoly(f)
            rh = ctx.from_cospoly(h)

            # value round-trip
            worst = max(worst, abs(rf.to_float() - f.evaluate(alpha)))
            # add / sub / mul
            worst = max(worst, abs((rf + rh).to_float()
                                   - (f + h).evaluate(alpha)))
            worst = max(worst, abs((rf - rh).to_float()
                                   - (f - h).evaluate(alpha)))
            worst = max(worst, abs((rf * rh).to_float()
                                   - (f * h).evaluate(alpha)))
            worst = max(worst, abs((rf * 3).to_float()
                                   - (f * 3).evaluate(alpha)))
            worst = max(worst, abs(rf.div2().to_float()
                                   - (f.evaluate(alpha) / 2)))

            # exact division by (1±cos2α): build g = f·(1±cos2α), divide back
            gp = f * CosPoly([1, 1], 0)
            gm = f * CosPoly([1, -1], 0)
            back_p = ctx.from_cospoly(gp).div_one_plus()
            back_m = ctx.from_cospoly(gm).div_one_minus()
            worst = max(worst, abs(back_p.to_float() - f.evaluate(alpha)))
            worst = max(worst, abs(back_m.to_float() - f.evaluate(alpha)))
            # and exactness: back_p should equal rf as a reduced ring element
            assert np.array_equal(back_p.reduce().vec, rf.copy().reduce().vec), \
                f"div_one_plus not exact for P/Q={P}/{Q}"
            assert np.array_equal(back_m.reduce().vec, rf.copy().reduce().vec), \
                f"div_one_minus not exact for P/Q={P}/{Q}"

        print(f"  P/Q={P}/{Q}: D={ctx.D:3d}  inv_c=({ctx.invp[1]},{ctx.invm[1]})  ok")
    print(f"ring_cyc self-test PASSED, worst value error = {worst:.2e}")


if __name__ == '__main__':
    _selftest()
