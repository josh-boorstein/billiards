"""The closed-form counts: Proposition 5.2 (P = 3), Proposition 5.3 (P = 2) and the table of
Proposition 9.4 (P = 1, Q - 1, Q - 2), at every centre of the partition record in their ranges.

Statements checked:
  * (Prop 5.2) for P = 3, Q >= 5, gcd(3, Q) = 1: the partition is [orphan] (A R)^m -- the orphan
    branch at one end of L1, then 2m boundaries strictly alternating between the two vertex
    types A and R -- with m = round(Q/3) - 1, so n(3/Q) = 2 round(Q/3) - 1;
  * (Prop 5.3) for P = 2, Q odd, Q >= 3: n(2/Q) = Q - 1;
  * (Prop 9.4) n = 1 at P = 1; n = Q - 1 at P = Q - 1; n = 2 floor((Q-1)/4) + 1 at P = Q - 2,
    Q odd.
The paper states ranges for these (Prop 5.2 certified to Q <= 26; Prop 9.4 for Q = 5..51, the
P = Q - 1 row to odd Q <= 29 and even Q <= 51); every centre of the record in the family is
checked, and the range actually covered is printed.
"""
from common import check, done, record


def main():
    rec = record()
    print('Closed-form counts\n')
    fam = {'P=3': [], 'P=2': [], 'P=1': [], 'P=Q-1': [], 'P=Q-2': []}
    for r in rec.values():
        P, Q = r['P'], r['Q']
        if P == 3 and Q >= 5 and Q % 3:
            fam['P=3'].append(r)
        if P == 2 and Q % 2 and Q >= 3:
            fam['P=2'].append(r)
        if P == 1:
            fam['P=1'].append(r)
        if P == Q - 1 and Q >= 3:
            fam['P=Q-1'].append(r)
        if P == Q - 2 and Q % 2 and Q >= 5:
            fam['P=Q-2'].append(r)

    def span(rows):
        qs = sorted(r['Q'] for r in rows)
        return f'{len(qs)} centres, Q = {qs[0]}..{qs[-1]}' if qs else 'none'

    rows = fam['P=3']
    cnt = sum(r['n'] == 2 * round(r['Q'] / 3) - 1 for r in rows)
    check(f'Prop 5.2: n(3/Q) = 2 round(Q/3) - 1 [{span(rows)}]: {cnt}/{len(rows)}',
          cnt == len(rows))

    def shape(r):
        m = round(r['Q'] / 3) - 1
        b, o = r['bounds'], r['orphans']
        alt = all(x in 'AR' for x in b) and all(b[i] != b[i + 1] for i in range(len(b) - 1))
        return len(b) == 2 * m and alt and o in ([0], [r['n'] - 1])
    cnt = sum(shape(r) for r in rows)
    check(f'Prop 5.2: the partition is [orphan](A R)^m, m = round(Q/3) - 1: {cnt}/{len(rows)}',
          cnt == len(rows))
    rows = fam['P=2']
    cnt = sum(r['n'] == r['Q'] - 1 for r in rows)
    check(f'Prop 5.3: n(2/Q) = Q - 1 at odd Q [{span(rows)}]: {cnt}/{len(rows)}', cnt == len(rows))
    for key, f in [('P=1', lambda Q: 1), ('P=Q-1', lambda Q: Q - 1),
                   ('P=Q-2', lambda Q: 2 * ((Q - 1) // 4) + 1)]:
        rows = fam[key]
        cnt = sum(r['n'] == f(r['Q']) for r in rows)
        check(f'Prop 9.4, row {key} [{span(rows)}]: {cnt}/{len(rows)}', cnt == len(rows))
    done()


if __name__ == '__main__':
    main()
