"""Section 5.2, Table 1: non-perpendicular periodic cylinders at rational cos 2a = y, q <= 20.

Each entry of the table is certified exactly from the word in `table1_record.json`
(`lib/mirror.py`); the class is read off the certificate.  The words were proposed by the
research searches' tracer; the certificate does not depend on how they were found.  The
table's blank cells are cells where the search found nothing, which is not a proof of absence,
and nothing here checks them.

Statements checked:
  * every rational y in the table gives an irrational right triangle (Niven: y not in
    {0, +-1/2}), with legs OB : OA as printed (OB the shorter, y > 0);
  * each of the 14 entries is a primitive, non-perpendicular periodic cylinder of the printed
    period, in the printed class;
  * the triangles 1 : sqrt(4m - 1), y = (2m - 1)/(2m), m = 2..5 carry cylinders in both the
    class perpendicular to OB and the class parallel to H.
"""
import math
from fractions import Fraction

from common import check, done, load
import mirror as m

TABLE = {  # y: (legs OB : OA as squares, {class: period})
    '1/4': ((3, 5), {'perp H': 24}),
    '1/3': ((1, 2), {'par H': 52}),
    '5/8': ((3, 13), {'perp H': 356}),
    '2/3': ((1, 5), {'perp H': 3600}),
    '7/10': ((3, 17), {'perp OB': 106, 'par H': 210}),
    '3/4': ((1, 7), {'perp OB': 16, 'par H': 30}),
    '5/6': ((1, 11), {'perp OB': 84, 'par H': 164}),
    '7/8': ((1, 15), {'perp OB': 338, 'par H': 670}),
    '9/10': ((1, 19), {'perp OB': 7356, 'par H': 14252}),
}


def main():
    rec = load('table1_record.json')['rows']
    got = {}
    for r in rec:
        y = Fraction(r['y'])
        K = m.t_field(r['minpoly'], r['y'])
        c = m.certify(K, r['word'])
        ok = c['ok'] and c['primitive'] and not c['perpendicular']
        got.setdefault(r['y'], {})[c.get('cls')] = c['period'] if ok else None
        print(f"    y = {r['y']:5s}  {c.get('cls', '-'):8s} period {c['period']:6d}  "
              f"width {c.get('width', 0):.3e}  certified {ok}", flush=True)
    for ys, ((s1, s2), cells) in TABLE.items():
        y = Fraction(ys)
        legs = (1 - y) / (1 + y)                       # (OB/OA)^2 = tan^2 a
        check(f'y = {ys}: a/pi irrational, legs OB : OA = sqrt{s1} : sqrt{s2}',
              y not in (0, Fraction(1, 2), Fraction(-1, 2)) and legs == Fraction(s1, s2)
              and m.angle_is_irrational(f'{y.denominator}*y - {y.numerator}'))
        check(f'y = {ys}: certified {cells}', got.get(ys) == cells)
    both = all({'perp OB', 'par H'} <= set(got.get(str(Fraction(2 * k - 1, 2 * k)), {}))
               for k in range(2, 6))
    check('y = (2m-1)/(2m), m = 2..5: both perp OB and par H', both)
    done()


if __name__ == '__main__':
    main()
