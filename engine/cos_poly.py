#!/usr/bin/env python3
"""
cos_poly.py — Exact arithmetic in Z[cos(2α)] for billiard computations.

A CosPoly represents an element of Z[1/2][cos(2α)] as an integer coefficient
vector [a₀, a₁, a₂, ...] meaning a₀ + a₁·cos(2α) + a₂·cos(4α) + ...
with a common power-of-2 denominator tracked separately.

Key operations:
- Addition/subtraction: componentwise
- Multiplication by cos(2kα): use cos(A)cos(B) = (cos(A+B) + cos(A-B))/2
- Multiplication by another CosPoly: distribute + product-to-sum
- Evaluation at a specific α: sum(a_k * cos(2kα)) / 2^shift
- Sign determination: evaluate at target α
"""

from fractions import Fraction


class CosPoly:
    """Polynomial in cos(2α) with rational coefficients.

    Stored as integer coefficients + a power-of-2 denominator (shift).
    Value = sum(coeffs[k] * cos(2kα) for k in range(len(coeffs))) / 2^shift
    """
    __slots__ = ('coeffs', 'shift')

    def __init__(self, coeffs=None, shift=0):
        if coeffs is None:
            self.coeffs = [0]
        else:
            self.coeffs = list(coeffs)
        self.shift = shift
        self._trim()

    def _trim(self):
        while len(self.coeffs) > 1 and self.coeffs[-1] == 0:
            self.coeffs.pop()

    def _reduce(self):
        """Cancel common factors of 2 between coeffs and shift."""
        if not any(self.coeffs) or self.shift == 0:
            return
        # Find min trailing zeros across all nonzero coefficients
        min_tz = self.shift
        for c in self.coeffs:
            if c != 0:
                tz = 0
                v = abs(c)
                while v & 1 == 0 and tz < min_tz:
                    v >>= 1
                    tz += 1
                min_tz = min(min_tz, tz)
                if min_tz == 0:
                    return
        if min_tz > 0:
            self.coeffs = [c >> min_tz for c in self.coeffs]
            self.shift -= min_tz

    def copy(self):
        p = CosPoly.__new__(CosPoly)
        p.coeffs = self.coeffs[:]
        p.shift = self.shift
        return p

    @staticmethod
    def constant(value):
        """Create a CosPoly from an integer or Fraction."""
        if isinstance(value, int):
            return CosPoly([value], 0)
        elif isinstance(value, Fraction):
            # Express as integer / power of 2 if possible
            d = value.denominator
            shift = 0
            while d > 1:
                if d % 2 == 0:
                    d //= 2
                    shift += 1
                else:
                    raise ValueError(f"Denominator {value.denominator} is not a power of 2")
            return CosPoly([value.numerator * (1 if shift == 0 else 1)], shift)
        else:
            raise TypeError(f"Cannot create CosPoly from {type(value)}")

    @staticmethod
    def cos_2k(k):
        """Create cos(2kα) as a CosPoly."""
        if k == 0:
            return CosPoly([1], 0)
        coeffs = [0] * (k + 1)
        coeffs[k] = 1
        return CosPoly(coeffs, 0)

    @staticmethod
    def from_cot_alpha():
        """cot(α) = cos(α)/sin(α). NOT in Z[cos(2α)] directly.

        Instead, represent cot²(α) = cos²α/sin²α = (1+cos2α)/(1-cos2α).
        This is a rational function, not a polynomial.

        For our purposes, we track coordinates SCALED by sin²α powers.
        """
        raise NotImplementedError("Use ExactTriangleState instead")

    def _align(self, other):
        """Return copies with matching shift (higher shift = smaller denominator)."""
        if self.shift == other.shift:
            return self.coeffs[:], other.coeffs[:], self.shift
        elif self.shift > other.shift:
            diff = self.shift - other.shift
            return self.coeffs[:], [c << diff for c in other.coeffs], self.shift
        else:
            diff = other.shift - self.shift
            return [c << diff for c in self.coeffs], other.coeffs[:], other.shift

    def __add__(self, other):
        if isinstance(other, int):
            other = CosPoly([other], 0)
        a, b, shift = self._align(other)
        n = max(len(a), len(b))
        a += [0] * (n - len(a))
        b += [0] * (n - len(b))
        result = CosPoly([a[i] + b[i] for i in range(n)], shift)
        result._reduce()
        return result

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        if isinstance(other, int):
            other = CosPoly([other], 0)
        a, b, shift = self._align(other)
        n = max(len(a), len(b))
        a += [0] * (n - len(a))
        b += [0] * (n - len(b))
        result = CosPoly([a[i] - b[i] for i in range(n)], shift)
        result._reduce()
        return result

    def __rsub__(self, other):
        if isinstance(other, int):
            other = CosPoly([other], 0)
        return other.__sub__(self)

    def __neg__(self):
        return CosPoly([-c for c in self.coeffs], self.shift)

    def __mul__(self, other):
        if isinstance(other, int):
            result = CosPoly([c * other for c in self.coeffs], self.shift)
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
        """Multiply two CosPolys using product-to-sum:
        cos(A)·cos(B) = (cos(A+B) + cos(A-B)) / 2
        """
        # Each term: self.coeffs[j] * cos(2jα) * other.coeffs[k] * cos(2kα)
        #          = self.coeffs[j] * other.coeffs[k] * (cos(2(j+k)α) + cos(2|j-k|α)) / 2
        max_degree = (len(self.coeffs) - 1) + (len(other.coeffs) - 1)
        result = [0] * (max_degree + 1)

        for j, cj in enumerate(self.coeffs):
            if cj == 0:
                continue
            for k, ck in enumerate(other.coeffs):
                if ck == 0:
                    continue
                prod = cj * ck
                # cos(2jα)·cos(2kα) = (cos(2(j+k)α) + cos(2|j-k|α)) / 2
                # But if j=0: cos(0)·cos(2kα) = cos(2kα), no division by 2
                # If k=0: cos(2jα)·cos(0) = cos(2jα), no division by 2
                if j == 0:
                    # 1 * cos(2kα) = cos(2kα)
                    if k >= len(result):
                        result.extend([0] * (k - len(result) + 1))
                    result[k] += prod
                elif k == 0:
                    # cos(2jα) * 1 = cos(2jα)
                    if j >= len(result):
                        result.extend([0] * (j - len(result) + 1))
                    result[j] += prod
                else:
                    # Need the /2: shift up by 1
                    # But we're accumulating before applying shift...
                    # Handle by pre-multiplying all terms by 2 and adding 1 to shift at the end
                    # Actually, let's handle this differently:
                    # We'll accumulate in a separate pass
                    pass

        # The above is getting complicated. Let's do it cleanly:
        # Represent everything with an extra shift for the products with j>0 and k>0
        max_degree = (len(self.coeffs) - 1) + (len(other.coeffs) - 1)
        result_coeffs = [0] * (max_degree + 1)
        extra_shift = 0  # How many extra /2 factors we need

        # Strategy: multiply all coefficients by 2 for each product pair,
        # then divide by 2 at the end. But nested products make this messy.
        #
        # Better: work with doubled coefficients throughout.
        # For terms involving cos(j)*cos(k) with j,k > 0:
        #   product = (cos(j+k) + cos(|j-k|)) / 2
        # So multiply the whole result by 2, then we get integer coefficients,
        # and add 1 to the shift.

        result_coeffs = [0] * (max_degree + 1)

        for j, cj in enumerate(self.coeffs):
            if cj == 0:
                continue
            for k, ck in enumerate(other.coeffs):
                if ck == 0:
                    continue
                prod = cj * ck
                if j == 0 or k == 0:
                    # No /2 needed: multiply by 2 to match the shifted representation
                    idx = j + k  # = max(j,k) since one is 0
                    if idx >= len(result_coeffs):
                        result_coeffs.extend([0] * (idx - len(result_coeffs) + 1))
                    result_coeffs[idx] += 2 * prod
                else:
                    # cos(2jα)·cos(2kα) = (cos(2(j+k)α) + cos(2|j-k|α)) / 2
                    # After multiplying by 2: just add prod to both indices
                    s = j + k
                    d = abs(j - k)
                    if s >= len(result_coeffs):
                        result_coeffs.extend([0] * (s - len(result_coeffs) + 1))
                    result_coeffs[s] += prod
                    result_coeffs[d] += prod

        new_shift = self.shift + other.shift + 1  # +1 for the factor of 2 we introduced
        result = CosPoly(result_coeffs, new_shift)
        result._reduce()
        return result

    def scale_by_cos2k(self, k):
        """Multiply by cos(2kα). Returns new CosPoly.

        cos(2kα)·cos(2jα) = (cos(2(k+j)α) + cos(2|k-j|α)) / 2
        """
        if k == 0:
            return self.copy()

        max_degree = len(self.coeffs) - 1 + k
        result = [0] * (max_degree + 1)

        for j, cj in enumerate(self.coeffs):
            if cj == 0:
                continue
            if j == 0:
                # 1 * cos(2kα) = cos(2kα), but in shifted form: 2*cj at index k
                result[k] += 2 * cj
            else:
                # cos(2jα) * cos(2kα) = (cos(2(j+k)α) + cos(2|j-k|α)) / 2
                # In shifted (+1) form: cj at each index
                s = j + k
                d = abs(j - k)
                if s >= len(result):
                    result.extend([0] * (s - len(result) + 1))
                result[s] += cj
                result[d] += cj

        r = CosPoly(result, self.shift + 1)
        r._reduce()
        return r

    def evaluate(self, alpha):
        """Evaluate at a specific angle α. Returns float."""
        import math
        val = sum(c * math.cos(2 * k * alpha) for k, c in enumerate(self.coeffs))
        return val / (2 ** self.shift)

    def sign_at(self, alpha):
        """Determine the sign at a specific α. Returns -1, 0, or 1."""
        val = self.evaluate(alpha)
        if abs(val) < 1e-15:
            return 0
        return 1 if val > 0 else -1

    def is_zero(self):
        return all(c == 0 for c in self.coeffs)

    def __repr__(self):
        terms = []
        for k, c in enumerate(self.coeffs):
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

    def degree(self):
        """Highest k with nonzero coefficient."""
        for k in range(len(self.coeffs) - 1, -1, -1):
            if self.coeffs[k] != 0:
                return k
        return 0


def test():
    """Verify basic arithmetic."""
    import math

    alpha = math.pi / 7
    c2a = math.cos(2 * alpha)

    # Test: 1 + cos(2α)
    p = CosPoly([1, 1], 0)
    assert abs(p.evaluate(alpha) - (1 + c2a)) < 1e-14

    # Test: cos(2α) * cos(2α) = (1 + cos(4α)) / 2
    a = CosPoly([0, 1], 0)  # cos(2α)
    b = a * a
    expected = (1 + math.cos(4 * alpha)) / 2
    assert abs(b.evaluate(alpha) - expected) < 1e-14, f"{b.evaluate(alpha)} vs {expected}"

    # Test: cos(2α) * cos(4α) = (cos(6α) + cos(2α)) / 2
    c = CosPoly([0, 0, 1], 0)  # cos(4α)
    d = a * c
    expected = (math.cos(6*alpha) + math.cos(2*alpha)) / 2
    assert abs(d.evaluate(alpha) - expected) < 1e-14

    # Test scale_by_cos2k
    p2 = CosPoly([3, -1, 2], 0)  # 3 - cos(2α) + 2cos(4α)
    scaled = p2.scale_by_cos2k(1)  # multiply by cos(2α)
    # Manually: 3cos(2α) + (-1)(1+cos(4α))/2 + 2(cos(6α)+cos(2α))/2
    #         = 3cos(2α) - 1/2 - cos(4α)/2 + cos(6α) + cos(2α)
    #         = -1/2 + 4cos(2α) - cos(4α)/2 + cos(6α)
    expected = p2.evaluate(alpha) * math.cos(2*alpha)
    assert abs(scaled.evaluate(alpha) - expected) < 1e-13, f"{scaled.evaluate(alpha)} vs {expected}"

    # Test: (1 - cos(2α)) / 2 = sin²α
    sin2 = CosPoly([1, -1], 1)  # (1 - cos(2α)) / 2
    assert abs(sin2.evaluate(alpha) - math.sin(alpha)**2) < 1e-14

    # Test polynomial multiplication
    # (1 + cos2α)(1 - cos2α) = 1 - cos²(2α) = sin²(2α) = (1-cos4α)/2
    p1 = CosPoly([1, 1], 0)
    p2 = CosPoly([1, -1], 0)
    prod = p1 * p2
    expected = (1 + c2a) * (1 - c2a)
    assert abs(prod.evaluate(alpha) - expected) < 1e-13, f"{prod.evaluate(alpha)} vs {expected}"

    print("All CosPoly tests passed!")
    print(f"  1 + cos(2α) = {CosPoly([1,1],0)} = {CosPoly([1,1],0).evaluate(alpha):.6f}")
    print(f"  cos²(2α) = {a*a}")
    print(f"  sin²α = {sin2} = {sin2.evaluate(alpha):.6f}")
    print(f"  (1+c2)(1-c2) = {prod}")


if __name__ == "__main__":
    test()
