"""Theorem 1.4's families on N <= 303 -- the polygon note, Section 9, third item: every member is completely
periodic, and p is 1 on the residues of Proposition 8.1 and (b - 2)/6, resp. (b - 1)/6, on those of Proposition 8.2
-- Proposition 8.2's lower bound p >= m, attained.

The families (b >= 4, 3 not dividing b):  F1 = (3, b, 2b) for b even, N = 3(b + 1);  F2 = (3, b, 2b + 3) for b odd,
N = 3(b + 2).  Proposition 8.1 covers b = 4 (mod 6) in F1 and b = 5 (mod 6) in F2 (one swapped pair); Proposition 8.2
covers b = 6m + 2 in F1 and b = 6m + 1 in F2 (m swapped pairs).  For every member with N <= 303 (65 triangles) the
tracer (`slots.trace_class`) follows all N prongs of the one class exactly, with a cap of 200 000 steps: s' = 0 (every
prong ends, so the class is completely periodic by Proposition 7.6) and p, the number of vertex-to-vertex prongs
matched with their reverses, equals the value above.  About a minute and a half.
"""
import time

import common
from common import check, done
from slots import classes, trace_class, pair_reversals

NMAX, CAP = 303, 200_000


def members():
    for b in range(4, NMAX):
        if b % 3 == 0:
            continue
        if b % 2 == 0 and 3 * (b + 1) <= NMAX:
            yield 'F1', b, (3, b, 2 * b)
        if b % 2 == 1 and b >= 5 and 3 * (b + 2) <= NMAX:
            yield 'F2', b, (3, b, 2 * b + 3)


def predicted_p(fam, b):
    if fam == 'F1':
        return (b - 2) // 6 if b % 6 == 2 else 1
    return (b - 1) // 6 if b % 6 == 1 else 1


def main():
    t0 = time.time()
    rows = []
    for fam, b, abc in members():
        c = classes(abc)[0]
        T, sl, geo = trace_class(abc, c['j'], cap=CAP)
        pairs, _ = pair_reversals(T, sl, geo)
        rows.append((fam, b, abc, len(sl), sum(r['end'] == 'cap' for r in sl),
                     None if pairs is None else len(pairs), predicted_p(fam, b)))
    print(f'{len(rows)} members, {time.time() - t0:.0f} s; largest p: '
          f'{max(rows, key=lambda r: r[5] or 0)[2]} with p = {max(r[5] or 0 for r in rows)}')
    check('65 members with N <= 303, all N slots traced on each', len(rows) == 65 and all(r[3] == sum(r[2]) for r in rows))
    check('every member completely periodic: s\' = 0', all(r[4] == 0 for r in rows))
    p81 = [r for r in rows if (r[0], r[1] % 6) in (('F1', 4), ('F2', 5))]
    p82 = [r for r in rows if (r[0], r[1] % 6) in (('F1', 2), ('F2', 1))]
    check(f'p = 1 on the {len(p81)} members of Proposition 8.1', all(r[5] == 1 for r in p81))
    check(f'p = (b - 2)/6, resp. (b - 1)/6, on the {len(p82)} members of Proposition 8.2', all(r[5] == r[6] for r in p82))
    check('every member is in one of the two propositions', len(p81) + len(p82) == len(rows))
    done()


if __name__ == '__main__':
    main()
