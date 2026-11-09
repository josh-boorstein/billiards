"""Isosceles triangles at even N with b < 2a -- the polygon note, Section 9, fifth item (continued): for every
I = (a, a, b) with even N = 2a + b <= 60, b < 2a, b != a and gcd(a, b) = 1 (91 triangles, 31 of them with a < b < 2a),
both classes of perpendicular directions are completely periodic, by an exact decomposition -- including the members
outside the hypotheses of Proposition 8.6's second sentence (b = 2, or b >= 4 and either N has at most one odd prime
factor or a is not +-2 modulo any odd prime q | N), where the note proves complete periodicity.  On the same instrument the minimal components are found (the
control): in a perpendicular class of (7, 7, 16) (b > a; found by Delecroix and Ruth [DR22]) and, for (4, 4, 7), in
the side directions, outside its perpendicular class.

`flow.decompose(I, j)` certifies every component of the unfolded surface in direction j a CYLINDER or WITHOUT
PERIODIC TRAJECTORY ('no_periodic', a minimal component), or leaves it undetermined at the step cap.  The two
classes are those of the base's and of a leg's perpendicular (`flow.perp_dirs`).  Prop 8.6 shows p >= 1 on the class
containing the base direction, so on these triangles the converse of Theorem 7.8 fails.  Under a minute.
"""
import math
import time

import common
from common import check, done
import flow

NMAX = 60


def odd_primes(n):
    return [q for q in range(3, n + 1, 2) if n % q == 0 and all(q % r for r in range(3, int(q ** 0.5) + 1, 2))]


def covered(a, b, N):
    """Proposition 8.6's second sentence: both classes proved completely periodic."""
    qs = odd_primes(N)
    return b == 2 or (b >= 4 and (len(qs) <= 1 or all(a % q not in (2, q - 2) for q in qs)))


def main():
    t0 = time.time()
    rows = []
    for N in range(4, NMAX + 1, 2):
        for a in range(1, N // 2):
            b = N - 2 * a
            if b <= 0 or b >= 2 * a or b == a or math.gcd(a, b) != 1:
                continue
            I = (a, a, b)
            pd = flow.perp_dirs(I)
            rows.append((I, [flow.decompose(I, j, cap=20000)['counts'] for j in (pd[0], pd[1])]))
    cp = lambda c: set(c) == {'cylinder'}
    acute = [r for r in rows if r[0][0] < r[0][2]]
    print(f'{len(rows)} triangles ({len(acute)} with a < b < 2a), {time.time() - t0:.0f} s')
    check('91 triangles with even N <= 60, b < 2a, b != a; 31 with a < b < 2a', (len(rows), len(acute)) == (91, 31))
    check('both classes completely periodic on all 91 (exact decomposition, every component a cylinder)',
          all(cp(c) for _, cs in rows for c in cs))
    out = [r[0] for r in rows if not covered(r[0][0], r[0][2], sum(r[0]))]
    check(f'... including the {len(out)} members outside the hypotheses of Proposition 8.6\'s second sentence ({out})',
          len(out) > 0 and all(cp(c) for r in rows if r[0] in out for c in r[1]))
    c1 = flow.decompose((7, 7, 16), flow.perp_dirs((7, 7, 16))[0], cap=20000)['counts']
    c2 = flow.decompose((4, 4, 7), flow.side_dirs((4, 4, 7))[0], cap=20000)['counts']
    c3 = [flow.decompose((4, 4, 7), j, cap=20000)['counts'] for j in flow.perp_dirs((4, 4, 7))[:2]]
    check(f'control: minimal components in the class of (7, 7, 16) perpendicular to its base ({c1})', 'no_periodic' in c1)
    check(f'control: minimal components in a side direction of (4, 4, 7) ({c2}), none in its perpendicular classes '
          f'({c3})', 'no_periodic' in c2 and all(cp(c) for c in c3))
    done()


if __name__ == '__main__':
    main()
