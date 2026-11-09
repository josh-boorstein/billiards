"""The odd-N census of the polygon note, Section 9, second item: among the 1578 scalene rational triangles with odd
N <= 51, p >= 1 occurs on exactly 18, every one among the 626 with some gcd(a_i, N) >= 3, and every one completely
periodic (s' = 0); with N <= 21 they are (2, 5, 8), (3, 4, 8), (3, 5, 13) and (4, 6, 11), where sum n_sigma = N + 1;
at the prime values N = 23, 29, ..., 47, p = 0 on every triangle.

At odd N there is one class of perpendicular directions and every side is vertical (Section 2).  For each triangle
(a < b < c, gcd 1) the tracer (`slots.trace_class`) follows all N prongs exactly: p is the number of vertex-to-vertex
prongs matched with their reverses, s' the number reaching no vertex within the cap.  A prong capped at 30 000 steps
is re-traced at 2 000 000 (escalation of the capped triangles only); the verdict counts are read after it.  For the
four N <= 21 triangles with p >= 1 the classifier (`beams.classify`) gives sum n_sigma over the three sides.

A capped prong is no verdict: s' = 0 is read only where every prong ended.  The census claims s' = 0 only on the 18
rows with p >= 1; the script prints how many rows of the 1578 still have a capped prong after escalation.
About four minutes.
"""
import time
from math import gcd

import common
from common import check, done
from beams import classify
from slots import classes, trace_class, pair_reversals

NMAX = 51
PRIMES = [23, 29, 31, 37, 41, 43, 47]


def census_row(abc, cap):
    c = classes(abc)[0]
    T, sl, geo = trace_class(abc, c['j'], cap=cap)
    pairs, _ = pair_reversals(T, sl, geo)
    return {'p': None if pairs is None else len(pairs), 'sprime': sum(r['end'] == 'cap' for r in sl), 'S': len(sl)}


def main():
    t0 = time.time()
    rows = {}
    for N in range(3, NMAX + 1, 2):
        for a in range(1, N):
            for b in range(a + 1, N):
                c = N - a - b
                if c > b and gcd(gcd(a, b), c) == 1:
                    rows[(a, b, c)] = census_row((a, b, c), 30_000)
    capped = [abc for abc, r in rows.items() if r['sprime']]
    for abc in capped:
        rows[abc] = census_row(abc, 2_000_000)
    still = [abc for abc, r in rows.items() if r['sprime']]
    print(f'{len(rows)} triangles, {len(capped)} with a prong capped at 30 000, {len(still)} still capped at 2 000 000; '
          f'{time.time() - t0:.0f} s')
    big_gcd = {abc for abc in rows if max(gcd(x, sum(abc)) for x in abc) >= 3}
    hit = sorted((abc for abc, r in rows.items() if r['p']), key=lambda t: (sum(t), t))
    print('p >= 1 on:', ', '.join(f'{t} (p = {rows[t]["p"]})' for t in hit))
    check('1578 scalene triangles with odd N <= 51, every prong matched or capped, S = N on each',
          len(rows) == 1578 and all(r['p'] is not None and r['S'] == sum(t) for t, r in rows.items()))
    check(f'{len(big_gcd)} = 626 of them have some gcd(a_i, N) >= 3', len(big_gcd) == 626)
    check(f'p >= 1 on exactly {len(hit)} = 18', len(hit) == 18)
    check('every p >= 1 triangle has some gcd(a_i, N) >= 3', set(hit) <= big_gcd)
    check('every p >= 1 triangle is completely periodic: s\' = 0', all(rows[t]['sprime'] == 0 for t in hit))
    small = [t for t in hit if sum(t) <= 21]
    check(f'with N <= 21 they are (2, 5, 8), (3, 4, 8), (3, 5, 13), (4, 6, 11) (got {small})',
          set(small) == {(2, 5, 8), (3, 4, 8), (3, 5, 13), (4, 6, 11)})
    sums = {t: sum(classify(t, s)['n'] for s in range(3)) for t in small}
    check(f'... where sum n_sigma = N + 1 (got {[(t, sums[t]) for t in small]})',
          all(sums[t] == sum(t) + 1 for t in small))
    check(f'at prime N = 23, ..., 47, p = 0 on every triangle',
          all(not r['p'] for t, r in rows.items() if sum(t) in PRIMES))
    done()


if __name__ == '__main__':
    main()
