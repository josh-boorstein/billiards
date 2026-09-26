"""Section 4, the second route to Theorem 2: the branch count n(8/15) = 6.

Statements checked:
  * Launched perpendicular to L1 (horizontally), the orbits from L1 fall into exactly
    six branches -- maximal intervals of launch heights sharing one reflection word up
    to the perpendicular (head-on) hit, after which the orbit retraces itself.  The six
    exact validity intervals TILE [0, 1] = L1, so there is no seventh branch.
  * Every launch other than the five cuts ends in a head-on hit, hence is periodic, so L1
    meets no minimal component (Section 3, "Why not the obvious strip"; Figure 1(a)).
  * The branches fall into three adjacent pairs of equal width, one pair per cylinder;
    the pair-sums total 1 = |L1|.  The half-words of a pair differ by one transposition
    of a '1' and a '2' (Appendix).
  * n(8/15) = 6 < Q - 1 = 14, so by [Orph, Thm 9.15] the direction is not completely
    periodic.

How: the branch words are proposed by the tracer and each is certified by its exact
validity interval in Q(zeta_60); the tiling is exact (`lib/triangle.py`).  This replaces
the partition tree the paper's Figure 1 caption refers to by an equivalent, smaller
certificate: a tiling of L1 by necessary-and-sufficient validity intervals.
"""
from t0 import T, Q, M_PERP, branch_words, check, done, trunc


def branches():
    """The exact tiling of L1 by the horizontal branch words, in increasing s."""
    words = branch_words()
    rows, tiles, defects = T.tile(words, M_PERP)
    return rows, tiles, defects


def cylinders(rows):
    """Pair the branches into cylinders: consecutive pairs of equal width."""
    cyl = []
    for a, b in zip(rows[0::2], rows[1::2]):
        cyl.append((a, b))
    return cyl


def main():
    print('Section 4 -- the branch count\n')
    rows, tiles, defects = branches()
    n = len(rows)
    for r in rows:
        print(f'  [{trunc(r["lo"].value(), 9)}, {trunc(r["hi"].value(), 9)})  '
              f'{r["word"]}  cut at top: {r["hi_src"][3] if r["hi_src"] else "end of L1"}')
    check(f'{n} branch words tile [0, 1] exactly', tiles and not defects)
    check('each ends at its first head-on hit (so every non-cut launch is periodic)',
          all(T.head_on_ok(r['word'], M_PERP) for r in rows))
    check('the five internal cuts are vertex grazes',
          all(r['hi_src'][3] in ('O', 'R', 'A') for r in rows[:-1]))
    cyl = cylinders(rows)
    eq_w = all((a['hi'] - a['lo']) == (b['hi'] - b['lo']) for a, b in cyl)
    check('three adjacent pairs of equal width', len(cyl) == 3 and eq_w)
    total = sum(((a['hi'] - a['lo']) + (b['hi'] - b['lo']) for a, b in cyl[1:]),
                (cyl[0][0]['hi'] - cyl[0][0]['lo']) + (cyl[0][1]['hi'] - cyl[0][1]['lo']))
    check('pair-sums total 1 = |L1|', total == 1)

    def one_swap(u, v):
        d = [i for i, (x, y) in enumerate(zip(u, v)) if x != y]
        return (len(u) == len(v) and len(d) == 2 and d[1] == d[0] + 1 and
                {u[d[0]], u[d[1]]} == {'1', '2'} and u[d[0]] == v[d[1]])
    check('the half-words of each pair differ by one transposition of 1 and 2',
          all(one_swap(a['word'], b['word']) for a, b in cyl))
    check(f'n(8/15) = {n} < Q - 1 = {Q - 1}', n == 6 and n < Q - 1)
    done()


if __name__ == '__main__':
    main()
