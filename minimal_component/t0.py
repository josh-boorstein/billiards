"""Shared set-up for the scripts of this directory: the triangle T0 and its words.

T0 is the right triangle with angles (4/15, 7/30, 1/2) pi: centre P/Q = 8/15,
alpha = 4 pi / 15.  Directions are integers m mod 60 standing for m pi / 30
(`lib/triangle.py`); the horizontal direction pi is m = 30, and the transversal of
Section 3, Sigma_Hbar = L1 read in direction theta_Hbar = 11 pi / 15, is m = 22.

The words are PROPOSED by the 60-digit tracer from a scan of L1 and then VERIFIED
exactly by the scripts that use them (a word's exact validity interval, and the exact
tiling of [0, 1] by those intervals).  A word the tracer missed would leave a gap in the
tiling, and a wrong word would have an empty interval; either fails loudly.
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'lib'))

from mpmath import mp, mpf, nstr  # noqa: E402

from triangle import RightTriangle  # noqa: E402
from tracer import Billiard, scan_words  # noqa: E402

P, Q = 8, 15
M_H = 2 * Q - P            # theta_Hbar = (2Q - P) pi / 2Q = 11 pi / 15
M_PERP = 2 * Q             # the horizontal direction, pi
SCAN = 1000                # launch points of the proposing scan

T = RightTriangle(P, Q)
F = T.F
KAPPA = T.L                # cot(alpha), the length of L2


def c(k):
    """cos(k pi / 30) in F."""
    return F.cos(k)


def cpi15(j):
    """cos(j pi / 15) in F -- the basis {1, cos(pi/15), cos(2pi/15), cos(pi/5)} of K."""
    return F.cos(2 * j)


K_BASIS = None


def k_coords(x):
    """Coordinates of x in the Q-basis {1, cos(pi/15), cos(2pi/15), cos(pi/5)} of
    K = Q(zeta_30)^+, solved exactly over Q; raises if x is not in K."""
    import sympy as sp
    global K_BASIS
    if K_BASIS is None:
        K_BASIS = [F.one, cpi15(1), cpi15(2), cpi15(3)]
    n = F.degree
    cols = [b.coeffs() + [0] * (n - len(b.coeffs())) for b in K_BASIS]
    rhs = x.coeffs() + [0] * (n - len(x.coeffs()))
    A = sp.Matrix(n, 4, lambda i, j: cols[j][i])
    sol, params = A.gauss_jordan_solve(sp.Matrix(rhs))
    assert not params, 'basis is not independent'
    return tuple(sol)


def return_words():
    """The first-return words to Sigma_Hbar, in increasing s."""
    B = Billiard(mp.pi * P / (2 * Q))
    seen = scan_words(B, M_H * mp.pi / (2 * Q), SCAN, kind='return')
    return sorted(seen, key=seen.get)


def branch_words():
    """The perpendicular-beam branch words on L1 (launch horizontal, up to the
    head-on hit), in increasing s."""
    B = Billiard(mp.pi * P / (2 * Q))
    seen = scan_words(B, mp.pi, SCAN, kind='head_on')
    return sorted(seen, key=seen.get)


def trunc(x, places):
    """Truncate (not round) a positive real to `places` decimals, as the paper prints."""
    s = nstr(x, places + 10, strip_zeros=False)
    whole, frac = s.split('.')
    return whole + '.' + frac[:places]


def check(label, ok):
    print(f'  [{"ok" if ok else "FAIL"}] {label}')
    if not ok:
        check.failures += 1
    return ok


check.failures = 0


def done():
    if check.failures:
        print(f'\n{check.failures} CHECK(S) FAILED')
        sys.exit(1)
    print('\nall checks passed')
