"""Theorem 5.5, the three-sided refinement, at the 116 coprime centres with 4 <= Q <= 19.

Statements checked, for the beams perpendicular to each of the three sides:
  * a side carries a retracing cell (a branch whose orbit meets another side head-on and
    retraces, J = id) exactly when:  L1 iff P odd,  L2 iff P + Q odd,  H iff Q odd;
  * so exactly two sides carry one and the third none, and each carries exactly ONE;
  * the pair is L1-L2 (P odd, Q even), L1-H (P odd, Q odd), L2-H (P even, Q odd);
  * the two retracing cells are the two ends of one double normal: each meets its partner
    side head-on, and their widths (as lengths, not fractions of the side) are EQUAL, decided
    exactly in Z[zeta_4Q];
  * every beam's partition is complete, tiles its side exactly, and has J an involution.
The paper also records an enumeration of its developed model at Q <= 60; that congruence count is
elementary and is not repeated here.
"""
import math

from common import check, done
from partition import Partition, partition

LETTER = {'1': 'L1', '2': 'L2', 'h': 'H'}


def retracing(P, Q, beam):
    r = partition(P, Q, beam=beam)
    assert r['open'] == 0
    return r, [r['leaves'][i] for i in r['orphans']]


def main():
    print('Theorem 5.5 -- the three-sided refinement\n')
    rows = [(P, Q) for Q in range(4, 20) for P in range(1, Q) if math.gcd(P, Q) == 1]
    pattern_ok = one_each = partner_ok = width_ok = engine_ok = 0
    for P, Q in rows:
        got, cells = {}, {}
        for b in '12h':
            r, c = retracing(P, Q, b)
            engine_ok += r['tiles'] and r['J_ok']
            got[b] = len(c)
            cells[b] = c
        exp = {'1': P % 2, '2': (P + Q) % 2, 'h': Q % 2}
        pattern_ok += got == exp
        one_each += sorted(got.values()) == [0, 1, 1]
        sides = [b for b in '12h' if got[b]]
        if len(sides) == 2:
            a, b = sides
            ca, cb = cells[a][0], cells[b][0]
            partner_ok += (ca['head_on'] == b and cb['head_on'] == a)
            R = Partition(P, Q).ring
            wa, wb = ca['hi'] - ca['lo'], cb['hi'] - cb['lo']
            # equal lengths, up to the orientation each placement gives its side
            width_ok += R.im_sign_exact(wa - wb) == 0 or R.im_sign_exact(wa + wb) == 0
    n = len(rows)
    check(f'{n} coprime centres with 4 <= Q <= 19', n == 116)
    check(f'retracing sides are L1 iff P odd, L2 iff P+Q odd, H iff Q odd: {pattern_ok}/{n}',
          pattern_ok == n)
    check(f'exactly two sides carry one retracing cell each, the third none: {one_each}/{n}',
          one_each == n)
    check(f'each retracing cell turns round on the partner side: {partner_ok}/{n}',
          partner_ok == n)
    check(f'the two cells have equal length (one double normal), exactly: {width_ok}/{n}',
          width_ok == n)
    check(f'all {3 * n} beam partitions complete, tiling, J an involution: {engine_ok}/{3 * n}',
          engine_ok == 3 * n)
    done()


if __name__ == '__main__':
    main()
