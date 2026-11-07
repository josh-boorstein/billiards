"""Remark 6.10: a non-perpendicular periodic cylinder in Troubetzkoy's perpendicular class, inside his range.

[Tro05] states, for an irrational right triangle whose smaller angle a satisfies pi/6 < a < pi/4, that the periodic
orbits perpendicular to the leg adjacent to a fill the invariant surface of that direction (abstract, (iii)). In the
paper's coordinates (angle a at A) that leg is OA and the surface is the class perp OA. Theorem 6.1 makes the
exceptional set of perp OA dense, so the claim fails on a dense set; this script certifies one explicit instance.

Statements checked:
  * y = (sqrt 7 - 1)/4 is the positive root of 8y^2 + 4y - 3, and 0 < y < 1/2 exactly, i.e. pi/6 < a < pi/4, so a is
    the smaller angle and OA the leg adjacent to it; a/pi is irrational;
  * the word in `tro05_record.json` is a primitive, non-perpendicular periodic cylinder in a mirror direction,
    certified EXACTLY (`lib/mirror.py`), of period 944, in the class perp OA at y (its width printed);
  * controls: the word does not close at the nearby rational cos 2a = 2/5, and no one-letter change of it certifies
    (a fixed sample of positions).
"""
import math
import random
from fractions import Fraction

import sympy as sp

from common import check, done, load, normalise
import mirror as m


def main():
    r = load('tro05_record.json')['rows'][0]
    s7 = sp.sqrt(7)
    yv = (s7 - 1) / 4
    check('(sqrt 7 - 1)/4 is the positive root of 8y^2 + 4y - 3', sp.expand(8 * yv ** 2 + 4 * yv - 3) == 0)
    check('0 < (sqrt 7 - 1)/4 < 1/2 exactly, i.e. pi/6 < a < pi/4', bool(yv > 0) and bool(yv < sp.Rational(1, 2)))
    check('a/pi irrational at cos 2a = (sqrt 7 - 1)/4', m.angle_is_irrational('8*y**2 + 4*y - 3'))
    K = m.t_field(r['minpoly'], float(r['y']))
    c = m.certify(K, r['word'])
    ay, cls = normalise(float(r['y']), c.get('cls'))
    print(f"    certified at y = {float(r['y']):+.6f} in class {c.get('cls', '-')} -> at |y|: {cls}, period {c['period']}, "
          f"width {c.get('width', 0):.3e}, a = {math.degrees(math.acos(ay) / 2):.3f} deg", flush=True)
    check('primitive non-perpendicular periodic cylinder of period 944 in the class perp OA at y = (sqrt 7 - 1)/4',
          c['ok'] and c['primitive'] and not c['perpendicular'] and c['period'] == 944 and cls == 'perp OA'
          and abs(ay - float(yv)) < 1e-12)

    # controls
    yc = Fraction(2, 5)
    Kc = m.t_field(f'{yc.denominator}*y + {yc.numerator}', str(-float(yc)))
    check('control: the word does not close at cos 2a = -2/5 (the leg-swapped 2/5)', not m.certify(Kc, r['word'])['ok'])
    rng = random.Random(892)
    word, bad, tried = r['word'], 0, 0
    for i in sorted(rng.sample(range(len(word)), 20)):
        for ch in 'ABH':
            v = word[:i] + ch + word[i + 1:]
            if v != word and m.cyclic_ok(v):
                tried += 1
                bad += bool(m.certify(K, v)['ok'])
    check(f'control: none of {tried} one-letter changes of the word certifies', bad == 0)
    done()


if __name__ == '__main__':
    main()
