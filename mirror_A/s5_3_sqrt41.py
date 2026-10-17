"""Section 5.3: the triangle cos 2a = (1 + sqrt 41)/8 carries non-perpendicular periodic
cylinders in all four mirror classes.

Each of the four words in `sqrt41_record.json` (rows of the research census; the census
itself is left for later work) is certified EXACTLY at its algebraic y (`lib/mirror.py`); the
class is read off the certificate, then normalised to the paper's convention (a triangle is
{y, -y}, listed at y > 0 where OB is the shorter leg).

Statements checked:
  * y = (1 + sqrt 41)/8 is the positive root of 8y^2 - 2y - 5, and a/pi is irrational;
  * each word is a primitive, non-perpendicular periodic cylinder in a mirror direction (its
    period printed), and the four classes at |y| are perp OB, perp OA, perp H and par H;
  * controls: no word closes at the nearby rational cos 2a = 37/40, and no one-letter change
    of a word certifies (a fixed sample of positions per word).
"""
import math
import random
from fractions import Fraction

from common import check, done, load, normalise
import mirror as m


def main():
    rows = load('sqrt41_record.json')['rows']
    yv = (1 + math.sqrt(41)) / 8
    check('(1 + sqrt 41)/8 is the positive root of 8y^2 - 2y - 5',
          abs(8 * yv * yv - 2 * yv - 5) < 1e-12 and (1 - math.sqrt(41)) / 8 < 0)
    check('a/pi irrational at cos 2a = (1 + sqrt 41)/8', m.angle_is_irrational('8*y**2 - 2*y - 5'))
    classes = set()
    certs = []
    for r in rows:
        K = m.t_field(r['minpoly'], float(r['y']))
        c = m.certify(K, r['word'])
        ok = c['ok'] and c['primitive'] and not c['perpendicular']
        ay, cls = normalise(float(r['y']), c.get('cls'))
        print(f"    {r['src']:5s} y = {float(r['y']):+.6f}  certified class {c.get('cls', '-'):8s} -> at |y|: "
              f"{cls:8s} period {c['period']:5d}  width {c.get('width', 0):.3e}", flush=True)
        check(f'{cls}: primitive non-perpendicular periodic cylinder of period {c["period"]}, at |y| = (1 + sqrt 41)/8',
              ok and abs(ay - yv) < 1e-12)
        classes.add(cls)
        certs.append((r, K))
    check('all four classes: perp OB, perp OA, perp H, par H', classes == {'perp OB', 'perp OA', 'perp H', 'par H'})

    # controls
    yc = Fraction(37, 40)
    for r, _ in certs:
        Kc = m.t_field(f'{yc.denominator}*y - {yc.numerator}' if float(r['y']) > 0
                       else f'{yc.denominator}*y + {yc.numerator}', str(math.copysign(float(yc), float(r['y']))))
        check(f"control: the period-{r['period']} word does not close at cos 2a = {'' if float(r['y']) > 0 else '-'}{yc}",
              not m.certify(Kc, r['word'])['ok'])
    rng = random.Random(735)
    for r, K in certs:
        word, bad, tried = r['word'], 0, 0
        for i in sorted(rng.sample(range(len(word)), 20)):
            for ch in 'ABH':
                v = word[:i] + ch + word[i + 1:]
                if v != word and m.cyclic_ok(v):
                    tried += 1
                    bad += bool(m.certify(K, v)['ok'])
        check(f"control: none of {tried} one-letter changes of the period-{r['period']} word certifies", bad == 0)
    done()


if __name__ == '__main__':
    main()
