"""Section 5, "What a census of perpendicular directions decides" (and the abstract):
every direction of S^1(2d') is an image, under the triangle's reflection group, of a
direction perpendicular to a side.

At the centre P/Q the right triangle has angles (P/2Q, (Q-P)/2Q, 1/2) pi, so d = d' = 2Q
and S^1(2d') is the set of multiples of pi/2Q -- the 4Q residues m mod 4Q in the units of
`lib/triangle.py`.  The reflection group acts on them by m -> 2Q - m, -m, 2P - m.  The
claim checked, at every coprime centre with Q <= 120 (default), is:
  * the orbit of the direction perpendicular to L1 (m = 0) is the even residues;
  * the odd residues form one further orbit, containing the direction perpendicular to L2
    (m = Q) at odd Q, and the direction perpendicular to H (m = P + Q) at even Q, where
    both legs' perpendiculars are even;
  * so the perpendicular direction classes exhaust S^1(2d').
Usage: python s5_coverage.py [QMAX]
"""
import sys
from math import gcd

from t0 import check, done


def orbit(m0, P, Q):
    M = 4 * Q
    seen, todo = {m0 % M}, [m0 % M]
    while todo:
        m = todo.pop()
        for n in ((2 * Q - m) % M, (-m) % M, (2 * P - m) % M):
            if n not in seen:
                seen.add(n)
                todo.append(n)
    return seen


def main(qmax=120):
    print(f'Section 5 -- perpendicular directions cover S^1(2d\'), coprime centres Q <= {qmax}\n')
    n = good = 0
    for Q in range(2, qmax + 1):
        for P in range(1, Q):
            if gcd(P, Q) != 1:
                continue
            n += 1
            M = 4 * Q
            e = orbit(0, P, Q)
            other = Q if Q % 2 else P + Q
            o = orbit(other, P, Q)
            ok = (e == set(range(0, M, 2)) and o == set(range(1, M, 2)))
            if Q % 2 == 0:
                ok = ok and (Q in e)                      # L2's perpendicular is even
            good += ok
            if not ok:
                print(f'  FAILS at {P}/{Q}')
    print(f'  {good} / {n} centres')
    check(f'every coprime centre with Q <= {qmax}: two orbits, even and odd, covering all 4Q',
          good == n)
    done()


if __name__ == '__main__':
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 120)
