#!/usr/bin/env python3.13
"""zlattice.py -- exact Z-lattice arithmetic on integer row vectors.

Hermite-style echelon basis, rank, membership, integer coordinates, and the class of a
vector in `Lambda / n*Lambda`.  Everything is exact: Python `int` only, no float, no
Fraction, no sympy.

WHY IT EXISTS.  Several strands need "is this exact algebraic number an integer
combination of those ones?" -- `ring_cyc` / `s450_chain_maps.Coord` already give an
element of `Z[zeta]` as an integer vector in a genuine Z-basis, so the question is
literally lattice membership, and the repo had no exact answer for it (`s448_width_rank`
computes a Q-RANK with Fractions, which is a different question: rank is about the
Q-span, membership is about the Z-span, and `x in span_Q` does not give `x in span_Z`).

    >>> L = Lattice([[2, 0], [0, 2]])
    >>> L.rank, L.contains([1, 1]), L.contains([2, 4]), L.coords([2, 4])
    (2, False, True, [1, 2])
    >>> L.klass([2, 4], 2)          # class in Lambda / 2 Lambda
    (1, 0)

⚠ TRAP -- CAST numpy int64 OUT FIRST.  `Coord.v` is `np.int64`; the gcd combinations
below multiply entries by cofactors and overflow SILENTLY at 2**63.  Pass
`[int(x) for x in v]`, or use `from_coords` which does it for you.

⚠ `klass(x, n)` is only defined for `x` IN the lattice; it returns `None` otherwise.  The
class is taken in the ECHELON basis this module computes, which is canonical for a given
generating set only up to the basis -- so compare classes from ONE `Lattice` object, never
across two built from different generators.
"""
from __future__ import annotations

from math import gcd


def _ext_gcd(a: int, b: int) -> tuple[int, int, int]:
    """(g, x, y) with `a*x + b*y == g == gcd(a, b)`, `g >= 0`."""
    old_r, r = a, b
    old_x, x = 1, 0
    old_y, y = 0, 1
    while r:
        q = old_r // r
        old_r, r = r, old_r - q * r
        old_x, x = x, old_x - q * x
        old_y, y = y, old_y - q * y
    if old_r < 0:
        return -old_r, -old_x, -old_y
    return old_r, old_x, old_y


class Lattice:
    """The Z-span of the given integer rows, as an echelon basis with distinct pivots."""

    __slots__ = ("dim", "basis")

    def __init__(self, rows):
        rows = [[int(x) for x in r] for r in rows]
        self.dim = len(rows[0]) if rows else 0
        self.basis: list[tuple[int, list[int]]] = []     # (pivot column, row), sorted
        for r in rows:
            self._insert(r)

    def _insert(self, r: list[int]) -> None:
        while True:
            p = next((i for i, x in enumerate(r) if x), None)
            if p is None:
                return
            hit = next((b for b in self.basis if b[0] == p), None)
            if hit is None:
                if r[p] < 0:
                    r = [-x for x in r]
                self.basis.append((p, r))
                self.basis.sort(key=lambda t: t[0])
                return
            _, br = hit
            a, b = r[p], br[p]
            g, x, y = _ext_gcd(a, b)
            nb = [x * u + y * v for u, v in zip(r, br)]          # pivot entry == g
            r = [(b // g) * u - (a // g) * v for u, v in zip(r, br)]   # r[p] == 0 now
            self.basis[self.basis.index(hit)] = (p, nb)

    @property
    def rank(self) -> int:
        return len(self.basis)

    def coords(self, x):
        """Integer coordinates of `x` in the echelon basis, or `None` if `x` is not in
        the lattice.  `sum(c[i] * basis[i] for i) == x`."""
        x = [int(v) for v in x]
        out = []
        for p, row in self.basis:
            if x[p] % row[p]:
                return None
            c = x[p] // row[p]
            out.append(c)
            if c:
                x = [u - c * v for u, v in zip(x, row)]
        return None if any(x) else out

    def contains(self, x) -> bool:
        return self.coords(x) is not None

    def klass(self, x, n: int = 2):
        """The class of `x` in `Lambda / n*Lambda`, as a tuple of residues; `None` if
        `x` is not in the lattice."""
        c = self.coords(x)
        return None if c is None else tuple(v % n for v in c)


def from_coords(coords):
    """`Lattice` from a list of `s450_chain_maps.Coord`-like objects (anything with a
    `.v` array): casts out of numpy, which this module requires (see the TRAP above)."""
    return Lattice([[int(t) for t in c.v] for c in coords])


def rational_rows(vals):
    """Integer rows of length 1 for a list of `Fraction`s, cleared by a common
    denominator -- so the same `Lattice` machinery covers the rank-1 rational case."""
    den = 1
    for v in vals:
        den = den * v.denominator // gcd(den, v.denominator)
    return [[int(v * den)] for v in vals]
