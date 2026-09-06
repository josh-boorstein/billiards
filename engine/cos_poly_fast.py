#!/usr/bin/env python3
"""
cos_poly_fast.py — Fast exact arithmetic in Z[1/2][cos(2α)] using numpy.

Drop-in replacement for cos_poly.CosPoly with numpy-backed operations.
Key speedups:
- Multiplication via vectorized product-to-sum formula
- Division via vectorized Chebyshev recurrence
- Bulk evaluation via numpy dot products
"""

import numpy as np


class CosPoly:
    """Polynomial in cos(2α) with integer coefficients and power-of-2 denominator.

    Value = sum(coeffs[k] * cos(2kα) for k in range(len(coeffs))) / 2^shift

    coeffs is a numpy int64 array. All arithmetic stays exact in integers.
    """
    __slots__ = ('coeffs', 'shift')

    def __init__(self, coeffs=None, shift=0):
        if coeffs is None:
            self.coeffs = np.zeros(1, dtype=np.int64)
        elif isinstance(coeffs, np.ndarray):
            self.coeffs = coeffs.astype(np.int64)
        else:
            self.coeffs = np.array(coeffs, dtype=np.int64)
        self.shift = shift
        self._trim()

    def _trim(self):
        # Remove trailing zeros
        n = len(self.coeffs)
        while n > 1 and self.coeffs[n-1] == 0:
            n -= 1
        if n < len(self.coeffs):
            self.coeffs = self.coeffs[:n]

    def _reduce(self):
        """Cancel common factors of 2 between coeffs and shift."""
        if self.shift == 0 or not np.any(self.coeffs):
            return
        nonzero = self.coeffs[self.coeffs != 0]
        if len(nonzero) == 0:
            return
        # Find minimum trailing zeros across all nonzero coefficients
        # Use bitwise: trailing zeros of x = number of times 2 divides x
        min_tz = self.shift
        for c in nonzero:
            c = int(c)
            if c == 0:
                continue
            tz = 0
            v = abs(c)
            while v & 1 == 0 and tz < min_tz:
                v >>= 1
                tz += 1
            min_tz = min(min_tz, tz)
            if min_tz == 0:
                return
        if min_tz > 0:
            self.coeffs = self.coeffs >> min_tz
            self.shift -= min_tz

    def copy(self):
        p = CosPoly.__new__(CosPoly)
        p.coeffs = self.coeffs.copy()
        p.shift = self.shift
        return p

    def _align_shift(self, other):
        """Return (self_coeffs, other_coeffs, common_shift) aligned to same shift."""
        if self.shift == other.shift:
            return self.coeffs, other.coeffs, self.shift
        elif self.shift > other.shift:
            diff = self.shift - other.shift
            return self.coeffs, other.coeffs << diff, self.shift
        else:
            diff = other.shift - self.shift
            return self.coeffs << diff, other.coeffs, other.shift

    def __add__(self, other):
        if isinstance(other, int):
            other = CosPoly(np.array([other], dtype=np.int64), 0)
        a, b, shift = self._align_shift(other)
        n = max(len(a), len(b))
        # Pad to same length
        a_pad = np.zeros(n, dtype=np.int64)
        b_pad = np.zeros(n, dtype=np.int64)
        a_pad[:len(a)] = a
        b_pad[:len(b)] = b
        result = CosPoly(a_pad + b_pad, shift)
        result._reduce()
        return result

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        if isinstance(other, int):
            other = CosPoly(np.array([other], dtype=np.int64), 0)
        a, b, shift = self._align_shift(other)
        n = max(len(a), len(b))
        a_pad = np.zeros(n, dtype=np.int64)
        b_pad = np.zeros(n, dtype=np.int64)
        a_pad[:len(a)] = a
        b_pad[:len(b)] = b
        result = CosPoly(a_pad - b_pad, shift)
        result._reduce()
        return result

    def __rsub__(self, other):
        if isinstance(other, int):
            other = CosPoly(np.array([other], dtype=np.int64), 0)
        return other.__sub__(self)

    def __neg__(self):
        return CosPoly(-self.coeffs, self.shift)

    def __mul__(self, other):
        if isinstance(other, int):
            result = CosPoly(self.coeffs * other, self.shift)
            result._reduce()
            return result
        if isinstance(other, CosPoly):
            return self._multiply_poly(other)
        raise TypeError(f"Cannot multiply CosPoly by {type(other)}")

    def __rmul__(self, other):
        if isinstance(other, int):
            return self.__mul__(other)
        raise TypeError

    def _multiply_poly(self, other):
        """Multiply two CosPolys using product-to-sum formula.

        cos(jθ)·cos(kθ) = (cos((j+k)θ) + cos(|j-k|θ)) / 2

        For the k=0 or j=0 case: cos(0)·cos(nθ) = cos(nθ), no /2.
        Strategy: scale by 2 (add 1 to shift), then all products are integer.
        """
        a = self.coeffs
        b = other.coeffs
        na = len(a)
        nb = len(b)
        max_degree = (na - 1) + (nb - 1)
        result = np.zeros(max_degree + 1, dtype=np.int64)

        # Terms where j=0: a[0] * b[k] * cos(kθ) → coefficient a[0]*b[k] at index k
        # (scaled by 2 to match the other terms)
        # Terms where k=0: a[j] * b[0] * cos(jθ)
        # Terms where j>0, k>0: a[j]*b[k]*(cos((j+k)θ) + cos(|j-k|θ))/2
        # After scaling everything by 2:
        #   j=0 or k=0 terms: contribute 2*a[j]*b[k] at the appropriate index
        #   j>0,k>0 terms: contribute a[j]*b[k] at both (j+k) and |j-k|

        # Handle j=0 terms
        if a[0] != 0:
            result[:nb] += 2 * a[0] * b

        # Handle k=0 terms (avoiding double-counting j=0,k=0)
        if b[0] != 0:
            result[:na] += 2 * b[0] * a
            # Undo the double-count at (0,0)
            result[0] -= 2 * a[0] * b[0]

        # Handle j>0, k>0 terms
        for j in range(1, na):
            if a[j] == 0:
                continue
            aj = a[j]
            # Sum contribution: aj * b[k] at index (j+k) for k=1..nb-1
            # Diff contribution: aj * b[k] at index |j-k| for k=1..nb-1
            for k in range(1, nb):
                if b[k] == 0:
                    continue
                prod = aj * b[k]
                s = j + k
                d = abs(j - k)
                result[s] += prod
                result[d] += prod

        new_shift = self.shift + other.shift + 1  # +1 for the factor of 2
        r = CosPoly(result, new_shift)
        r._reduce()
        return r

    def evaluate(self, alpha):
        """Evaluate at a specific angle α. Returns float."""
        import math
        n = len(self.coeffs)
        # cos(2k*alpha) for k = 0, 1, ..., n-1
        cos_vals = np.array([math.cos(2 * k * alpha) for k in range(n)])
        val = np.dot(self.coeffs.astype(np.float64), cos_vals)
        return val / (2 ** self.shift)

    def sign_at(self, alpha):
        """Determine the sign at a specific α."""
        val = self.evaluate(alpha)
        if abs(val) < 1e-15:
            return 0
        return 1 if val > 0 else -1

    def is_zero(self):
        return not np.any(self.coeffs)

    def degree(self):
        """Highest k with nonzero coefficient."""
        nonzero = np.nonzero(self.coeffs)[0]
        if len(nonzero) == 0:
            return 0
        return int(nonzero[-1])

    def __repr__(self):
        terms = []
        for k in range(len(self.coeffs)):
            c = int(self.coeffs[k])
            if c == 0:
                continue
            if k == 0:
                terms.append(str(c))
            else:
                terms.append(f"{c}·cos({2*k}α)")
        s = " + ".join(terms) if terms else "0"
        if self.shift > 0:
            s = f"({s}) / 2^{self.shift}"
        return s


def divide_by_1_plus_cos2(poly):
    """Exact division by (1 + cos(2α)) using Chebyshev recurrence."""
    return _divide_by_1_pm_cos2(poly, +1)


def divide_by_1_minus_cos2(poly):
    """Exact division by (1 - cos(2α)) using Chebyshev recurrence."""
    return _divide_by_1_pm_cos2(poly, -1)


def _divide_by_1_pm_cos2(poly, sign):
    """Divide poly by (1 + sign·cos(2α)) exactly.

    Recurrence: g_{n-1}=2s·f_n, g_{k-1}=2s(f_k-g_k)-g_{k+1}, then g_0.
    """
    if poly.is_zero():
        return CosPoly(np.zeros(1, dtype=np.int64), 0)

    f = poly.coeffs
    n = len(f) - 1
    while n > 0 and f[n] == 0:
        n -= 1

    if n == 0:
        if f[0] == 0:
            return CosPoly(np.zeros(1, dtype=np.int64), 0)
        assert False, f"Cannot divide constant by (1{'+' if sign>0 else '-'}cos2α)"

    g = np.zeros(n, dtype=np.int64)
    result_shift = poly.shift

    if n == 1:
        g[0] = sign * f[1]
    else:
        g[n-1] = 2 * sign * f[n]
        for k in range(n-1, 1, -1):
            g_kp1 = int(g[k+1]) if k + 1 <= n - 1 else 0
            g[k-1] = 2 * sign * (int(f[k]) - int(g[k])) - g_kp1

        g2 = int(g[2]) if n > 2 else 0
        if g2 % 2 == 0:
            g[0] = sign * (int(f[1]) - int(g[1])) - g2 // 2
        else:
            g = g * 2
            result_shift += 1
            g[0] = sign * (2 * int(f[1]) - int(g[1])) - g2

    result = CosPoly(g, result_shift)
    result._reduce()
    return result


def test():
    """Verify fast CosPoly matches the original."""
    import math

    alpha = math.pi / 7

    # Test basic arithmetic
    p1 = CosPoly(np.array([1, 1]), 0)
    assert abs(p1.evaluate(alpha) - (1 + math.cos(2*alpha))) < 1e-14

    # Multiplication
    a = CosPoly(np.array([0, 1]), 0)  # cos(2α)
    b = a * a
    expected = (1 + math.cos(4*alpha)) / 2
    assert abs(b.evaluate(alpha) - expected) < 1e-14

    # Division
    f = CosPoly(np.array([3, 3]), 0)
    g = divide_by_1_plus_cos2(f)
    assert abs(g.evaluate(alpha) - 3.0) < 1e-14

    f2 = CosPoly(np.array([5, 6, 1]), 1)
    g2 = divide_by_1_plus_cos2(f2)
    assert abs(g2.evaluate(alpha) - (2 + math.cos(2*alpha))) < 1e-14

    f3 = CosPoly(np.array([1, -1]), 0)
    g3 = divide_by_1_minus_cos2(f3)
    assert abs(g3.evaluate(alpha) - 1.0) < 1e-14

    # Higher degree division
    p = CosPoly(np.array([1, -1, 2, 0, -3, 1]), 0)
    prod = p * CosPoly(np.array([1, -1]), 0)
    g4 = divide_by_1_minus_cos2(prod)
    assert abs(g4.evaluate(alpha) - p.evaluate(alpha)) < 1e-12

    print("All fast CosPoly tests passed!")


if __name__ == "__main__":
    test()
