"""Theorem 6.2 (1)-(2), the beams of an isosceles triangle at even N, on every I = (a, a, b) with even N = 2a + b <= 30
and gcd(a, b) = 1 (47 triangles) -- the polygon note, Section 9, fourth item.

I has angles pi*(a, a, b)/N at V_0, V_1 (the base vertices) and V_2 (the apex); side 0 = V0V1 is the base, sides 1, 2
the legs.  Its right half T has 2 alpha/pi = P/Q := (b/2)/(N/2) (in lowest terms, as gcd(a, b) = 1), written
(P, Q - P, Q): apex O = V_0, base vertex A = V_1, right angle R = V_2, so side 0 is the hypotenuse H, side 1 the half
base L1, side 2 the axis L2.  n_H, n(P/Q) are the numbers of branches of T's beams perpendicular to H and to L1, and
nu_H, nu their numbers of retracing branches.  Every count is the classifier's (`beams.classify`), exact.  Checked:
  (1) legs: every branch of each leg beam retraces, onto the other leg, and n_leg = (n_H + nu_H)/2;
  (2) base: no branch retraces; there is one mirror pair if nu = 1 (and then Q is even), none otherwise, and every
      other branch lies in a translation pair, (n(P/Q) - nu)/2 of them.
A minute or two.
"""
import time
from math import gcd

import common
from common import check, done
from beams import classify

NMAX = 30


def main():
    t0 = time.time()
    rows = []
    for N in range(4, NMAX + 1, 2):
        for b in range(2, N - 1, 2):
            a = (N - b) // 2
            if gcd(a, b) != 1:
                continue
            I = (a, a, b)
            P, Q = b // 2, N // 2
            T = (P, Q - P, Q)
            base, leg1, leg2 = (classify(I, s) for s in range(3))
            H, L1 = classify(T, 0), classify(T, 1)
            ok = all(r['ok'] for r in (base, leg1, leg2, H, L1))
            legs = ok and all(lg['R'] == lg['n'] and set(lg['retr_by_side']) <= {other}
                              and lg['n'] == (H['n'] + H['R']) / 2
                              for lg, other in ((leg1, 2), (leg2, 1)))
            nu = L1['R'] if ok else None
            base_ok = ok and base['R'] == 0 and base['M'] == (2 if nu == 1 else 0) \
                and (nu == 0 or Q % 2 == 0) and base['T'] / 2 == (L1['n'] - nu) / 2
            rows.append((I, (P, Q), ok, legs, base_ok, nu))
    print(f'{len(rows)} triangles, {time.time() - t0:.0f} s; nu = 1 at {[r[1] for r in rows if r[5] == 1]}')
    check('47 isosceles triangles with even N <= 30, every beam of I and of its half classified',
          len(rows) == 47 and all(r[2] for r in rows))
    check('Theorem 6.2 (1): every leg branch retraces onto the other leg, n_leg = (n_H + nu_H)/2',
          all(r[3] for r in rows))
    check('Theorem 6.2 (2): no base branch retraces; one mirror pair iff nu = 1 (Q even); (n(P/Q) - nu)/2 '
          'translation pairs', all(r[4] for r in rows))
    for r in rows:
        if not (r[2] and r[3] and r[4]):
            print('  failing:', r)
    done()


if __name__ == '__main__':
    main()
