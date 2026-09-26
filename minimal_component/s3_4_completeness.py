"""Section 3.4: why "three branches" is complete, and the mass budget.

Statements checked:
  * Each of the seven words is a FIRST return to Sigma_Hbar in direction theta_Hbar: no
    earlier letter reflects off L1 into theta_Hbar (direction arithmetic only, so exact).
  * Each word's region of validity is a single interval with exact endpoints, cut out by
    18 to 198 affine constraints; the seven intervals TILE [0, 1] -- consecutive
    endpoints are equal as field elements, no gap, no overlap -- including at s -> 0 and
    s -> 1, where no constraint binds.  Since the itinerary at each s is unique, this
    excludes an eighth word.
  * Of the six internal boundaries exactly two change the translation, at l1 and l1 + l2,
    and all six are vertex grazes (the binding constraint on each side names a vertex).
  * The mass budget: half of area(B) is 6.7530303322337...; mu(M) exceeds it by
    0.6177788...; each of the two iota-blocks carries 3.68540457604279214..., 54.57% of
    its disc; and the excess is not (sqrt 5 - 1)/2 = 0.6180339887....
"""
from mpmath import sqrt, mpf

from t0 import T, F, Q, KAPPA, M_H, return_words, check, done, trunc
from s3_1_return_map import lambdas
from s3_3_mass_identity import kac_side


def main():
    print('Section 3.4 -- completeness of the three branches, and the mass budget\n')
    words = return_words()
    check('all seven are first-return words', all(T.first_return_ok(w, M_H) for w in words))
    rows, tiles, defects = T.tile(words, M_H)
    ncon = [r['n_constraints'] for r in rows]
    print(f'  constraints per word: {ncon}')
    check('18 to 198 constraints each', min(ncon) == 18 and max(ncon) == 198)
    check('the seven validity intervals tile [0, 1] exactly', tiles and not defects)
    check('no constraint binds as s -> 0 or s -> 1',
          rows[0]['lo_src'] is None and rows[-1]['hi_src'] is None)
    grazes = all(a['hi_src'][3] is not None and b['lo_src'][3] is not None
                 for a, b in zip(rows, rows[1:]))
    check('all six internal boundaries are vertex grazes', grazes)
    l1, l2, _ = lambdas()
    shifts = [T.unfold(r['word'], M_H)[-1][4][1] for r in rows]
    changes = [a['hi'] for a, sa, sb in zip(rows, shifts, shifts[1:]) if not sa == sb]
    check('exactly two change the translation, at l1 and l1 + l2',
          len(changes) == 2 and changes[0] == l1 and changes[1] == l1 + l2)

    mu, _ = kac_side()
    half = KAPPA * Q / 2
    excess = mu - half
    block = mu / 2
    print(f'\n  area(B)/2 = {trunc(half.value(), 13)};  mu(M) - area(B)/2 = '
          f'{trunc(excess.value(), 7)};  each block {trunc(block.value(), 17)}')
    check('half of area(B) = 6.7530303322337...', trunc(half.value(), 13) == '6.7530303322337')
    check('mu(M) exceeds it by 0.6177788...', trunc(excess.value(), 7) == '0.6177788')
    check('each block 3.68540457604279214...', trunc(block.value(), 17) == '3.68540457604279214')
    pct = block.value() / half.value() * 100
    check(f'  = {trunc(pct, 2)}% of its disc', trunc(pct, 2) == '54.57')
    golden = (sqrt(5) - 1) / 2
    check('the excess is not (sqrt 5 - 1)/2 (they agree to three places only)',
          abs(excess.value() - golden) > mpf('1e-4') and trunc(golden, 10) == '0.6180339887')
    done()


if __name__ == '__main__':
    main()
