"""Remark 8.5 at odd N <= 39 -- the polygon note, Section 9, fifth item: for every isosceles triangle I = (a, a, b) with
odd N = 2a + b <= 39 and gcd(a, b) = 1 (158 triangles), the one class of I is completely periodic with p = 0 --
157 read off an exact certificate for the right half, and (17, 17, 3), whose half needs separatrices of length above
10^9, by an exact decomposition of the direction; the exact decomposition of each triangle's own directions agrees
on all 158.

Two instruments, independent of each other:
  (a) THE HALF.  I's unfolding is the unfolding of its right half T, at 2 alpha/pi = b/N, and I's class is T's class
      in which the hypotenuse and the half base are vertical: eps = 0 at the odd denominator N (Remark 8.5;
      `closure.py`).  `necklace.readout(b, N, 0)` walks every separatrix of that class exactly: `capped == 0` PROVES
      complete periodicity, and on a certified class `swapped` is the number of swapped pairs of T, so
      p_T^{OA} = p_T^{OR} = 0 and p_I = 2 p_T^{OA} + p_T^{OR} = 0.  A capped walk is no verdict.
  (b) THE TRIANGLE ITSELF.  `flow.decompose` on I in the direction perpendicular to its base (at odd N every side's
      perpendicular is in the one class): every component certified a cylinder.
(17, 17, 3) is run in (a) with a cap of 10^5 steps and reported capped; (b) certifies it, and Remark 8.5's
equivalence (N = 37 is prime, so unconditionally) gives p = 0 there.  About ten minutes, most of it the half of
(10, 10, 19), whose longest separatrix has 12 733 910 steps.
"""
import time
from math import gcd

import common
from common import check, done
import flow
from necklace import readout

NMAX, CAP, HARD = 39, 20_000_000, (17, 17, 3)


def main():
    t0 = time.time()
    rows = []
    for N in range(3, NMAX + 1, 2):
        for b in range(1, N - 1, 2):
            a = (N - b) // 2
            if gcd(a, b) != 1:
                continue
            I = (a, a, b)
            h = readout(b, N, 0, step_cap=10**5 if I == HARD else CAP)
            d = flow.decompose(I, flow.perp_dirs(I)[0], cap=20000)
            rows.append((I, h['capped'], h['swapped'], h['max_steps'], d['counts']))
    print(f'{len(rows)} triangles, {time.time() - t0:.0f} s; longest separatrix of a certified half: '
          f'{max(r[3] for r in rows if r[1] == 0)} steps')
    cert = [r for r in rows if r[1] == 0]
    cyl = [r for r in rows if set(r[4]) == {'cylinder'}]
    check('158 isosceles triangles with odd N <= 39', len(rows) == 158)
    check(f'(a) the half\'s class certified completely periodic on {len(cert)} = 157, all but (17, 17, 3)',
          len(cert) == 157 and HARD not in [r[0] for r in cert])
    check('(a) on each of the 157, no swapped pair in the half: p_I = 0', all(r[2] == 0 for r in cert))
    check(f'(b) the decomposition of I\'s own direction is all cylinders on {len(cyl)} = 158, (17, 17, 3) included',
          len(cyl) == 158)
    hard = next(r for r in rows if r[0] == HARD)
    print(f'  (17, 17, 3): the half capped ({hard[1]} separatrices at 10^5 steps); its own decomposition {hard[4]}')
    done()


if __name__ == '__main__':
    main()
