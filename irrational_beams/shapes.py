"""[Deposit] Section 5, second paragraph: the fan decomposition of the certified census angles with at most four beams. Paper:
"Finite perpendicular-beam partitions at irrational angles". Research name: data/scratch/s979/classify_n4.py (a scratch driver over
fan_trace.fans, the floating-point fan scan in the polygon P), rewritten for the deposit to count DISTINCT certified angles (the
research driver counted rows, so x = 2/5, found by two coincidences, twice, and it predates the certification of x = 0.1211654606;
the reconciliation is the research ledger's, and the counts below are the paper's).

A FLOAT trace, as the paper says: it proposes the fans, it certifies nothing. Classes:
  single   -- one fan;
  resolved -- exactly two fans, one of steps (m, m), the other palindromic, both closing perpendicularly (fan_trace 'headon');
  partial  -- the scan closes the (m, m) fan and resolves the second fan, of width below 0.004, only in part (a 'cap' piece).
Input: the certified points of data/census_summary_<mode>.json's source rows (data/census_exact_<mode>.json over the record).
Output: data/shapes.json. Paper: 192 angles = 6 single (five roots of (F) with J = 1, and x = 2/3) + 184 resolved + 2 partial.
The scan is run at each point's full double-precision x. (At x rounded to ten digits the research reconciliation counted 183 + 3:
at x = 0.1211654606 the rounding leaves a capped sliver inside the 0.0011-wide second fan, which at the record's x, equal to the
root of the certified minimal polynomial to 3e-18, is resolved as one palindromic fan closing perpendicularly.)
"""
import json
import sys

import _setup  # noqa: F401  (deposit: this directory on the path, cwd = irrational_beams/)
from fan_trace import fans


def certified_rows(mode):
    rows = json.load(open('data/census_record.json'))['rows']
    redo = {round(float(r['x']), 10): r for r in json.load(open(f'data/census_exact_{mode}.json'))}
    out = {}
    for r in rows:
        if int(r['n9']) > 4:
            continue
        o = redo.get(round(float(r['x']), 10), r)
        if str(o.get('CERTIFIED')) == 'True':
            out[round(float(r['x']), 9)] = float(r['x'])
    return sorted(out.values())


def classify(x):
    m = int(1 / x)
    F = fans(x)
    pal = [list(st) == list(st)[::-1] for _, _, st, _ in F]
    mm = [list(st) == [m, m] for _, _, st, _ in F]
    if len(F) == 1:
        return 'single', F
    if len(F) == 2 and sum(mm) == 1 and all(pal) and all(w == 'headon' for *_, w in F):
        return 'resolved', F
    capped = [(lo, hi) for lo, hi, _, w in F if w != 'headon']
    rest = [hi - lo for (lo, hi, st, w), is_mm in zip(F, mm) if not is_mm]
    if sum(mm) == 1 and capped and sum(rest) < 0.004:
        return 'partial', F
    return 'other', F


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else 'sample'
    xs = certified_rows(mode)
    res = {}
    for x in xs:
        c, F = classify(x)
        res.setdefault(c, []).append(dict(x=x, fans=[(lo, hi, list(st), w) for lo, hi, st, w in F]))
    counts = {k: len(v) for k, v in res.items()}
    singles = sorted(round(r['x'], 10) for r in res.get('single', []))
    json.dump(dict(counts=counts, angles=len(xs), classes=res), open('data/shapes.json', 'w'))
    print(f'{len(xs)} certified angles with n <= 4: {counts}; single-fan x = {singles}')
    ok = len(xs) == 192 and counts == {'single': 6, 'resolved': 184, 'partial': 2}
    print('ALL OK' if ok else 'FAIL')


if __name__ == '__main__':
    main()
