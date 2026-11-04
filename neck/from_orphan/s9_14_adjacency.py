"""Remark 9.14: at every certified centre, the interior branch boundaries grazing a copy of R are
exactly as many as the mirror pairs.

The proof of Proposition 9.2(ii) makes every R-graze boundary pair-internal, and distinct ones give
distinct pairs, so  #{R-graze boundaries} <= #{mirror pairs}  always, with equality exactly when
every pair is adjacent.  The count is read off the partition: a boundary's source vertex is
recorded with it, and #{mirror pairs} = (n - #orphans) / 2.  Every centre of the record, including
8/15, where complete periodicity fails; a slice (Q <= 12) is recomputed from scratch.
"""
from common import FULL, centre, check, done, record


def counts(r):
    return r['bounds'].count('R'), (r['n'] - len(r['orphans'])) // 2


def main():
    rec = record()
    rows = sorted(rec.values(), key=lambda r: (r['Q'], r['P']))
    print(f'Remark 9.14 -- adjacency count at {len(rows)} certified centres, '
          f'Q up to {max(r["Q"] for r in rows)}\n')
    eq = [r for r in rows if counts(r)[0] == counts(r)[1]]
    check(f'#R-graze boundaries = #mirror pairs: {len(eq)}/{len(rows)}', len(eq) == len(rows))
    check('the count includes 8/15 (not completely periodic)', '8/15' in rec)
    check(f'... and both P <= 2 families: {sum(r["P"] <= 2 for r in rows)} centres',
          any(r['P'] == 1 for r in rows) and any(r['P'] == 2 for r in rows))
    slice_ = [r for r in rows if r['Q'] <= (201 if FULL else 12)]
    same = sum(counts(centre(r['P'], r['Q'], fresh=True)) == counts(r) for r in slice_)
    check(f'recomputed from scratch, Q <= {201 if FULL else 12}: {same}/{len(slice_)}',
          same == len(slice_))
    done()


if __name__ == '__main__':
    main()
