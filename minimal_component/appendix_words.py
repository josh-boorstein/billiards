"""Appendix: the words of the certificate.

Prints the two tables of the Appendix -- the seven return words of Section 3.1 with their
validity intervals, and the six branch half-words of the three cylinders of Section 3.3 --
from the exact computation, and checks them against the tables as printed in the paper.
Words list the sides met in order: '1' = L1, '2' = L2, 'h' = H.
"""
from t0 import T, M_H, M_PERP, return_words, branch_words, check, done, trunc

PAPER_RETURN = [
    (6, '0.000000', '0.581886', 'h12h21'),
    (52, '0.581886', '0.604227', 'h1h2h21h21h1h2h21h2h1h12h12h2h1h12h21h12h21h12h12h21'),
    (46, '0.604227', '0.720227', 'h1h2h21h21h1h2h21h2h1h12h12h2h1h12h21h12h21h21'),
    (46, '0.720227', '0.813284', 'h1h2h21h21h1h2h12h2h1h12h12h2h1h12h21h12h21h21'),
    (46, '0.813284', '0.858568', 'h1h2h21h21h1h2h12h2h1h12h12h2h1h12h21h12h12h21'),
    (66, '0.858568', '0.929284',
     'h1h2h21h21h1h2h12h2h1h12h1h2h21h2h1h12h12h2h1h12h21h12h21h12h12h21'),
    (66, '0.929284', '1.000000',
     'h1h2h21h21h1h2h12h2h1h21h1h2h21h2h1h12h12h2h1h12h21h12h21h12h12h21'),
]
PAPER_CYL = [
    ('C1', 25, 13, 'h2h1h12h1h2h1'), ('C1', 25, 13, 'h2h1h21h1h2h1'),
    ('C2', 53, 27, 'h2h1h2h1h12h12h21h1h2h1h2h1'), ('C2', 53, 27, 'h2h1h2h1h12h21h21h1h2h1h2h1'),
    ('C3', 9, 5, 'h12h1'), ('C3', 9, 5, 'h21h1'),
]


def six(x):
    """Six places, the paper's convention for this table (rounded)."""
    return f'{float(x.value(30)):.6f}'


def main():
    print('Appendix -- the words of the certificate\n')
    rows, tiles, _ = T.tile(return_words(), M_H)
    got = [(len(r['word']), six(r['lo']), six(r['hi']), r['word']) for r in rows]
    for g in got:
        print(f'  {g[0]:3d}  [{g[1]}, {g[2]})  {g[3]}')
    check('the seven return words and intervals, as printed', tiles and got == PAPER_RETURN)

    rows, tiles, _ = T.tile(branch_words(), M_PERP)
    by_word = {r['word']: r for r in rows}
    ok = tiles
    for name, N, hits, w in PAPER_CYL:
        ok &= w in by_word and len(w) == hits and 2 * hits - 1 == N
    # the pairs: equal widths
    for (n1, _, _, w1), (n2, _, _, w2) in zip(PAPER_CYL[0::2], PAPER_CYL[1::2]):
        a, b = by_word[w1], by_word[w2]
        ok &= n1 == n2 and (a['hi'] - a['lo']) == (b['hi'] - b['lo'])
    print()
    for name, N, hits, w in PAPER_CYL:
        print(f'  {name}  N = {N:2d}  hits {hits:2d}  {w}')
    check('the six cylinder half-words, as printed (pairs of equal width)', ok)
    done()


if __name__ == '__main__':
    main()
