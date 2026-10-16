"""Section 5.1: Examples 5.1 and 5.2, certified exactly.

Statements checked:
  * (Ex 5.1) cos 2a = 3/4 (legs 1 : sqrt 7; a/pi irrational): the word HAHABHABHBAHAHAB is a
    primitive periodic cylinder of period 16 in the class perpendicular to OB, direction 2a;
    it unfolds to the translation 4 (cos 2a, sin 2a) (OA = 1); its launch interval on OB is
    (0, sqrt7/42) and its perpendicular width sqrt7/56 -- all exact; its smallest angles to
    the normals of the sides are arcsin(1/8) (OA), 0.236.. (H) and 2a (OB); no path meets a
    side at a right angle.
  * (Ex 5.2) cos 2a = 1/4 (legs sqrt 3 : sqrt 5): the word HABHBAHAHBHAHABHAHBHAHAB is a
    primitive non-perpendicular periodic cylinder of period 24 in the class perpendicular to
    H, certified in Q(2 cos a) = Q(sqrt 10), with launch interval (7 sqrt15/129, 5 sqrt15/43)
    on OB and smallest angle to a normal 0.09898...
  * Controls: the same words do NOT certify at a nearby rational angle (cos 2a = 3/4 + 1/1000,
    1/4 + 1/1000), and a one-letter change of either word does not certify.
"""
import math
from fractions import Fraction

from common import check, done
import mirror as m


def side_angles(cert, word):
    """Smallest angle to the normal, per side letter, over one period (floats)."""
    tri = cert['tri']
    edges, _ = tri.unfold(word)
    bz = cert['beta'].complex(30)
    out = {}
    for (P, Q), ch in zip(edges, word):
        ez = (Q - P).complex(30)
        c = abs((bz.conjugate() * ez).real) / (abs(bz) * abs(ez))
        out[ch] = min(out.get(ch, 9), math.asin(min(1.0, float(c))))
    return out


def main():
    # Example 5.1
    w1 = 'HAHABHABHBAHAHAB'
    K = m.t_field('4*y - 3', '3/4')
    c = m.certify(K, w1)
    check('Ex 5.1: a/pi irrational (cos 2a = 3/4)', m.angle_is_irrational('4*y - 3'))
    check('Ex 5.1: periodic cylinder, primitive, period 16', c['ok'] and c['primitive'] and c['period'] == 16)
    check('Ex 5.1: direction 2a (class perpendicular to OB, e = 2)', c.get('cls') == 'perp OB' and abs(c['e']) == 2)
    tri = c['tri']
    # translation, in the frame OA = 1: beta / t; claim 4 e(2a) = 4 w^2
    check('Ex 5.1: unfolds to the translation 4 (cos 2a, sin 2a) exactly',
          (c['beta'] - tri.wpow(2) * (4 * tri.t)).iszero())
    l1, l2 = m.launch_interval(c, 'B')
    # y-coordinate on OB = lam * tan a; tan a = 1/sqrt 7, so sqrt7/42 <-> lam = 1/6
    check('Ex 5.1: launch interval on OB is (0, sqrt7/42), i.e. lam in (0, 1/6) exactly',
          l1.iszero() and (l2 - Fraction(1, 6)).iszero())
    check('Ex 5.1: width sqrt7/56 exactly (width^2 = 7/3136)', (c['width2'] - Fraction(7, 3136)).iszero())
    ang = side_angles(c, w1)
    a = math.acos(0.75) / 2
    check(f'Ex 5.1: smallest angles to the normals OA {ang["A"]:.6f} = arcsin(1/8), '
          f'H {ang["H"]:.4f}, OB {ang["B"]:.6f} = 2a',
          abs(ang['A'] - math.asin(1 / 8)) < 1e-12 and abs(ang['H'] - 0.236) < 5e-4
          and abs(ang['B'] - 2 * a) < 1e-12)
    check('Ex 5.1: no path meets a side at a right angle', not c['perpendicular'])

    # Example 5.2
    w2 = 'HABHBAHAHBHAHABHAHBHAHAB'
    K2 = m.t_field('4*y - 1', '1/4')
    c2 = m.certify(K2, w2)
    check('Ex 5.2: a/pi irrational (cos 2a = 1/4)', m.angle_is_irrational('4*y - 1'))
    check('Ex 5.2: periodic cylinder, primitive, period 24, class perpendicular to H, non-perpendicular',
          c2['ok'] and c2['primitive'] and c2['period'] == 24 and c2.get('cls') == 'perp H'
          and not c2['perpendicular'])
    l1, l2 = m.launch_interval(c2, 'B')
    # tan a = sqrt(3/5) = sqrt15/5: 7 sqrt15/129 <-> lam = 35/129, 5 sqrt15/43 <-> lam = 75/129
    check('Ex 5.2: launch interval on OB is (7 sqrt15/129, 5 sqrt15/43) exactly',
          (l1 - Fraction(35, 129)).iszero() and (l2 - Fraction(75, 129)).iszero())
    check(f'Ex 5.2: smallest angle to a normal {c2["min_angle_to_normal"]:.6f} = 0.09898..',
          0.09898 <= c2['min_angle_to_normal'] < 0.09899)

    # controls
    for word, y in ((w1, Fraction(751, 1000)), (w2, Fraction(251, 1000))):
        Kc = m.t_field(f'{y.denominator}*y - {y.numerator}', str(y))
        check(f'control: {word[:6]}.. does not close at cos 2a = {y}', not m.certify(Kc, word)['ok'])
    for word, K_ in ((w1, K), (w2, K2)):
        bad = 0
        for i in range(len(word)):
            for ch in 'ABH':
                v = word[:i] + ch + word[i + 1:]
                if v != word and m.cyclic_ok(v) and m.certify(K_, v)['ok']:
                    bad += 1
        check(f'control: no one-letter change of {word[:6]}.. certifies', bad == 0)
    done()


if __name__ == '__main__':
    main()
