#!/usr/bin/env python3
"""
n_exact.py — EXACT branch counter n(P/Q) for the perpendicular rational
partition of 2alpha = (P/Q)*pi  (i.e. alpha = (P/Q)*(pi/2)).

Why this exists
---------------
`word_tree_nofrac.build_partition_tree` dedups boundaries by comparing FLOAT
positions with a fixed min_width. That is unreliable two ways (see paper_audit
and session_state):
  (a) OVER-count: the tree emits spurious ~1e-15-width boundaries at the s->0
      (and s->1) edges — their EXACT algebraic value is 0 (resp. 1) but float
      noise makes them look like distinct interior boundaries.
  (b) UNDER-count at large q: genuine distinct boundaries can sit closer than
      the float tolerance and get merged.

Fix: dedup boundaries by their EXACT algebraic value, not floats.

The mechanism
-------------
Every boundary is a CosPoly:  value = (1/2^s) * sum_k c_k cos(2k*alpha).
With alpha = (P/Q)*(pi/2):   2k*alpha = k*P*pi/Q = 2*pi*(kP)/N,  N = 2Q.
So cos(2k*alpha) = cos(2*pi*(kP)/N) = (zeta^{kP} + zeta^{-kP}) / 2,
where zeta = exp(2*pi*i/N) is a primitive N-th root of unity.

Thus 2^{s+1} * value = sum_k c_k (zeta^{kP mod N} + zeta^{(-kP) mod N})
is an element of the cyclotomic integer ring Z[zeta_N] = Z[x]/(Phi_N(x)).
Representing it as an integer vector of length deg(Phi_N)=phi(N) is CANONICAL:
two boundaries have equal algebraic value IFF their reduced (vector, shift)
keys are identical. This is exact — no floats in the comparison.

Zero-width branches: a boundary whose exact key equals that of 0 (the s->0
edge) or 1 (the s->1 edge) is dropped, as is any duplicate of an interior
boundary (same key). What survives are the genuine interior splits; the branch
count n = #distinct-interior-boundaries + 1.
"""

import math
from functools import lru_cache

import numpy as np

from word_tree_nofrac import build_partition_tree


# ----------------------------------------------------------------------------
# Cyclotomic polynomial Phi_N(x) as integer coeffs (low -> high), monic.
# ----------------------------------------------------------------------------
def _polydiv_exact(a, b):
    """Divide polynomial a by monic-ish b exactly; assumes exact division.

    Coeffs low->high. b's leading coeff must divide as we go (true for
    x^n-1 / product-of-cyclotomics: leading coeff 1)."""
    a = a[:]
    db = len(b) - 1
    lead = b[-1]
    q = [0] * (len(a) - db)
    for i in range(len(a) - 1, db - 1, -1):
        coeff = a[i]
        if coeff == 0:
            continue
        assert coeff % lead == 0, "non-exact cyclotomic division"
        c = coeff // lead
        q[i - db] = c
        for j in range(len(b)):
            a[i - db + j] -= c * b[j]
    # remainder a[:db] must be all zero
    assert all(v == 0 for v in a[:db]), "nonzero remainder in cyclotomic div"
    return q


@lru_cache(maxsize=None)
def cyclotomic_poly(n):
    """Phi_n(x), integer coeffs low->high, monic. Degree = phi(n)."""
    # x^n - 1
    num = [-1] + [0] * (n - 1) + [1]
    for d in range(1, n):
        if n % d == 0:
            num = _polydiv_exact(num, list(cyclotomic_poly(d)))
    return tuple(num)


def _reduction_table(N):
    """red[m] = x^m mod Phi_N(x) as an int list of length D=deg(Phi_N),
    for m in [0, N)."""
    phi = list(cyclotomic_poly(N))
    D = len(phi) - 1                      # degree = phi(N)
    # Phi_N is monic: leading coeff phi[D] == 1
    assert phi[D] == 1
    red = []
    cur = [0] * D
    cur[0] = 1                            # x^0 = 1
    red.append(cur[:])
    for _ in range(1, N):
        # multiply cur by x  (shift up), then reduce mod Phi_N if degree hits D
        nxt = [0] + cur[:]               # length D+1, top coeff at index D
        top = nxt[D]
        if top != 0:
            # x^D = -(phi[0..D-1]) since Phi_N monic:  x^D + sum phi[i] x^i = 0
            for i in range(D):
                nxt[i] -= top * phi[i]
            nxt[D] = 0
        cur = nxt[:D]
        red.append(cur[:])
    return red, D


# ----------------------------------------------------------------------------
# Exact value key of a boundary CosPoly at alpha = (P/Q)*(pi/2).
# ----------------------------------------------------------------------------
def _canonical(vec, shift):
    """Reduce (integer vector / 2^shift) by cancelling common powers of 2.
    Returns a hashable canonical key (tuple(vec), shift). Zero -> ((0,),0)."""
    if not any(vec):
        return ((0,) * len(vec), 0)
    min_tz = shift
    for c in vec:
        if c == 0:
            continue
        v = abs(int(c))
        tz = 0
        while (v & 1) == 0 and tz < min_tz:
            v >>= 1
            tz += 1
        if tz < min_tz:
            min_tz = tz
        if min_tz == 0:
            break
    if min_tz > 0:
        vec = [c >> min_tz for c in vec]
        shift -= min_tz
    return (tuple(vec), shift)


def exact_key(coeffs, shift, P, N, red, D):
    """Canonical exact key of value = (1/2^shift) sum_k coeffs[k] cos(2k alpha).

    2^{shift+1} * value = sum_k coeffs[k] (x^{kP mod N} + x^{-kP mod N})
    reduced mod Phi_N(x)."""
    W = [0] * D
    for k, ck in enumerate(coeffs):
        ck = int(ck)
        if ck == 0:
            continue
        m = (k * P) % N
        rm = red[m]
        rn = red[(N - m) % N]
        for i in range(D):
            W[i] += ck * (rm[i] + rn[i])
    return _canonical(W, shift + 1)


# ----------------------------------------------------------------------------
# The exact branch counter.
# ----------------------------------------------------------------------------
def n_pq_exact(P, Q, max_depth=None, closure_tol=1e-6, return_detail=False,
               min_width=1e-11):
    """Exact branch count n(P/Q) of the perpendicular rational partition of
    2alpha = (P/Q)*pi.  P/Q must be in lowest terms with P/Q < 1.

    max_depth defaults to a generous 12*Q (bump for even-P centers, which
    hide thin high-period branches — validate stability by raising it).

    min_width=1e-11 prunes the spurious ~1e-15-width degenerate corridors the
    tree emits at the s->0 edge (documented filter). Genuine boundaries at the
    q we study have width >> 1e-11, so this is safe and it keeps those spurious
    non-closing corridors from spoiling the completeness certificate."""
    g = math.gcd(P, Q)
    P, Q = P // g, Q // g
    if not (0 < P < Q):
        raise ValueError("need 0 < P < Q (P/Q in (0,1))")
    alpha = (P / Q) * (math.pi / 2)
    N = 2 * Q
    if max_depth is None:
        max_depth = 12 * Q + 60

    _, _, polys, n_trunc = build_partition_tree(
        alpha, max_depth=max_depth, closure_tol=closure_tol,
        min_width=min_width, return_raw=True)

    red, D = _reduction_table(N)
    zero_key = _canonical([0] * D, 0)
    one_key = exact_key(np.array([1], dtype=np.int64), 0, P, N, red, D)

    distinct = {}          # key -> float value (for reporting)
    n_zero = n_one = n_dup = 0
    for poly, _vtype in polys:
        key = exact_key(poly.coeffs, poly.shift, P, N, red, D)
        if key == zero_key:
            n_zero += 1
            continue
        if key == one_key:
            n_one += 1
            continue
        if key in distinct:
            n_dup += 1
            continue
        distinct[key] = poly.evaluate(alpha)

    n = len(distinct) + 1
    if return_detail:
        return {
            'P': P, 'Q': Q, 'n': n, 'raw_polys': len(polys),
            'distinct_boundaries': len(distinct),
            'dropped_zero': n_zero, 'dropped_one': n_one, 'dropped_dup': n_dup,
            'boundaries': sorted(distinct.values()), 'max_depth': max_depth,
            'depth_truncations': n_trunc, 'complete': n_trunc == 0,
        }
    return n


def n_pq_certified(P, Q, start_depth=None, cap=200000, verbose=False):
    """Exact, PROVABLY-COMPLETE n(P/Q).

    Increases the depth cap until the tree reports zero depth-truncations
    (every branch closed via orbit closure) — then the count is exact and
    certified complete, immune to the plateau trap (a count that is constant
    over several depths but still hiding thin high-period branches).

    Returns (n, certified, depth_used). certified=False means the cap was hit
    before completeness (the partition genuinely needs a deeper tree)."""
    g = math.gcd(P, Q)
    P, Q = P // g, Q // g
    depth = start_depth if start_depth is not None else 16 * Q
    det = None
    while depth <= cap:
        det = n_pq_exact(P, Q, max_depth=depth, return_detail=True)
        if verbose:
            print(f"    depth {depth}: n={det['n']} "
                  f"trunc={det['depth_truncations']}", flush=True)
        if det['complete']:
            return det['n'], True, depth
        depth *= 2
    if det is None:                       # start_depth already exceeded cap
        det = n_pq_exact(P, Q, max_depth=cap, return_detail=True)
        return det['n'], det['complete'], cap
    return det['n'], False, depth // 2


if __name__ == '__main__':
    print("SANITY: cyclotomic degrees")
    for N in (6, 8, 10, 16, 26):
        print(f"  deg Phi_{N} = {len(cyclotomic_poly(N)) - 1} "
              f"(phi={sum(1 for k in range(1, N) if math.gcd(k, N) == 1)})")

    print("\nEXACT-KEY vs float cross-check (values that ARE equal must key-eq)")
    # cos(4a) and 2cos(2a)^2-1 are the same value: keys must match.
    from cos_poly_fast import CosPoly
    P, Q = 3, 8
    N = 2 * Q
    red, D = _reduction_table(N)
    a = CosPoly(np.array([0, 0, 1], dtype=np.int64), 0)          # cos(4a)
    b = CosPoly(np.array([0, 1], dtype=np.int64), 0)
    b = b * b * 2 - 1                                            # 2cos(2a)^2-1
    ka = exact_key(a.coeffs, a.shift, P, N, red, D)
    kb = exact_key(b.coeffs, b.shift, P, N, red, D)
    alpha = (P / Q) * (math.pi / 2)
    print(f"  cos4a float={a.evaluate(alpha):.6f}  2c2^2-1 float="
          f"{b.evaluate(alpha):.6f}  keys_equal={ka == kb}")

    print("\nVALIDATION: Theorem C  n(2/q) = q-1  (even numerator P=2)")
    for q in [7, 9, 11, 13, 15, 19, 23, 25, 31, 35, 41]:
        d = n_pq_exact(2, q, return_detail=True)
        ok = 'OK' if d['n'] == q - 1 else f"!! expected {q-1}"
        print(f"  n(2/{q:<3}) = {d['n']:<3} (raw {d['raw_polys']}, "
              f"zero {d['dropped_zero']}, dup {d['dropped_dup']})  {ok}")
