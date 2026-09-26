"""Section 3.1: the first return map to Sigma_Hbar is an exchange of three intervals.

Statements checked:
  * The first return to Sigma_Hbar = {(cot alpha, s)} with outgoing theta_Hbar = 11pi/15
    has seven return-word classes with return counts N = 6 / 52 / 46,46,46 / 66,66.
  * Along each word the return map is s -> s + c EXACTLY (the returning point lies on L1
    and the s-coefficient is exactly 1), with only three distinct translations,
    +(l2+l3), (l3-l1), -(l1+l2).
  * The six internal boundaries (from the exact tiling, `s3_4_completeness.py`) are all
    vertex grazes, at the positions, vertices and segments of the paper's table; the
    translation changes at exactly two of them, s = l1 and s = l1 + l2.
  * l1 = 1 + 4cos(pi/5) - 4cos(2pi/15) = 1 - 4 sin(pi/30),
    l2 = 3 + 4cos(pi/15) - 4cos(2pi/15) - 4cos(pi/5),
    l3 = -3 - 4cos(pi/15) + 8cos(2pi/15);
    coordinates (1,0,-4,4), (3,4,-4,-4), (-3,-4,8,0) in the basis
    {1, cos(pi/15), cos(2pi/15), cos(pi/5)} of K, summing to (1,0,0,0).
  * The printed decimals (truncated, 19 places).

How: each word is unfolded exactly in Q(zeta_60) (`lib/triangle.py`); every cut is one
linear solve over that field.  No value is fitted or recognised from a decimal.  The
words are proposed by the tracer and certified by their exact validity intervals.
"""
from t0 import T, F, M_H, cpi15, c, k_coords, return_words, check, done, trunc

PAPER_TABLE = [  # boundary: s (19 places), vertex, segment, words, translation changes
    ('0.5818861469293861144', 'O', 3, (6, 52), True),
    ('0.6042267417944153876', 'R', 45, (52, 46), True),
    ('0.7202272239678215794', 'R', 16, (46, 46), False),
    ('0.8132836683297223304', 'R', 42, (46, 46), False),
    ('0.8585683010062570444', 'O', 27, (46, 66), False),
    ('0.9292841505031285222', 'R', 23, (66, 66), False),
]


def lambdas():
    l1 = 1 + 4 * cpi15(3) - 4 * cpi15(2)
    l2 = 3 + 4 * cpi15(1) - 4 * cpi15(2) - 4 * cpi15(3)
    l3 = -3 - 4 * cpi15(1) + 8 * cpi15(2)
    return l1, l2, l3


def main():
    print('Section 3.1 -- the first return map to Sigma_Hbar\n')
    words = return_words()
    Ns = [len(w) for w in words]
    check(f'seven return words, N = {Ns}', Ns == [6, 52, 46, 46, 46, 66, 66])

    l1, l2, l3 = lambdas()
    sin_pi30 = c(14)                                   # sin(pi/30) = cos(14 pi/30)
    check('l1 = 1 - 4 sin(pi/30)', l1 == 1 - 4 * sin_pi30)
    coords = [k_coords(x) for x in (l1, l2, l3)]
    check(f'coordinates {coords}',
          coords == [(1, 0, -4, 4), (3, 4, -4, -4), (-3, -4, 8, 0)])
    check('l1 + l2 + l3 = 1 exactly', l1 + l2 + l3 == 1)
    for name, x, want in (('l1', l1, '0.5818861469293861144'),
                          ('l2', l2, '0.0223405948650292732'),
                          ('l3', l3, '0.3957732582055846123')):
        check(f'{name} = {trunc(x.value(), 19)}...', trunc(x.value(), 19) == want)

    # each word: exact translation
    shifts = []
    for w in words:
        st = T.unfold(w, M_H)
        Ax, Ay = st[-1][4]
        Bx, By = st[-1][5]
        ok = Ax == T.L and Bx.is_zero() and By == 1 and st[-1][6] == M_H
        check(f'N={len(w):3d}: returns to L1 in direction theta_Hbar as s -> s + c, '
              f'c = {trunc(abs(Ay.value()), 12)}{"" if Ay > 0 else " (negative)"}', ok)
        shifts.append(Ay)
    want = [l2 + l3, l3 - l1] + [-(l1 + l2)] * 5
    check('translations: +(l2+l3), (l3-l1), -(l1+l2) x5', all(a == b for a, b in zip(shifts, want)))

    # the boundaries, from the exact tiling
    rows, tiles, defects = T.tile(words, M_H)
    check('the seven validity intervals tile [0,1] (see s3_4_completeness.py)', tiles)
    print('\n  boundary  s                        arrives at  segment  words     translation')
    changes = []
    for b, (a, nxt, paper) in enumerate(zip(rows, rows[1:], PAPER_TABLE), 1):
        s = a['hi']
        j, side, name, vtx = a['hi_src']
        j2, _, _, vtx2 = nxt['lo_src']
        same = (shifts[b - 1] == shifts[b])
        ch = not same
        if ch:
            changes.append(s)
        print(f'  {b}         {trunc(s.value(), 19)}  {vtx}           {j:3d}      '
              f'{len(a["word"])} -> {len(nxt["word"])}  {"changes" if ch else "unchanged"}')
        ok = (trunc(s.value(), 19) == paper[0] and vtx == paper[1] == vtx2 and
              j == paper[2] and (len(a['word']), len(nxt['word'])) == paper[3] and
              ch == paper[4])
        check(f'boundary {b} matches the paper\'s row', ok)
    check('exactly two boundaries change the translation, at l1 and l1 + l2',
          len(changes) == 2 and changes[0] == l1 and changes[1] == l1 + l2)
    done()


if __name__ == '__main__':
    main()
