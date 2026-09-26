"""Section 3.2: minimality and unique ergodicity -- Lemma 3, Lemma 4, the roof.

Statements checked:
  * Lemma 3: the 3 x 4 integer matrix of coordinates of l1, l2, l3 (Section 3.1) has
    rank 3 over Q, so l1, l2, l3 are linearly independent over Q.
  * The symmetric permutation (3 2 1) is irreducible (Keane's condition then applies).
  * Lemma 4: tau is the first return to [0,1) of the rotation by l2 + l3 of the circle
    R/(1 + l2)Z, with return time 1 on [0, l1) u [l1 + l2, 1) and 2 on [l1, l1 + l2):
    the three cases are exact identities in K (with l1 + l2 + l3 = 1), and a numerical
    simulation on random points agrees.  (l2 + l3)/(1 + l2) = p/q forces
    p l1 + (2p - q) l2 + (p - q) l3 = 0.
  * The tower over ([0,1), tau) with that roof has total length 1 + l2, the circle's
    length (l1 + 2 l2 + l3 = 1 + l2 exactly).
  * The roof r (length of the first-return path to Sigma_Hbar) is constant on each of the
    seven words -- the s-coefficient of the path length is exactly zero -- and takes the
    four values 2.972579301909..., 19.136181234849..., 22.108760536759...,
    27.750141644774..., so 2.9 < r < 27.8.
"""
import random

import sympy as sp
from mpmath import mp, mpf, nstr

from t0 import T, M_H, k_coords, return_words, check, done, trunc
from s3_1_return_map import lambdas


def main():
    print('Section 3.2 -- minimality and unique ergodicity\n')
    l1, l2, l3 = lambdas()
    Mx = sp.Matrix([k_coords(x) for x in (l1, l2, l3)])
    print(f'  coordinate matrix = {Mx.tolist()}')
    check(f'Lemma 3: rank {Mx.rank()} = 3 over Q', Mx.rank() == 3)
    perm = (3, 2, 1)
    irreducible = all(set(perm[:k]) != set(range(1, k + 1)) for k in range(1, 3))
    check('(3 2 1) is irreducible', irreducible)

    # Lemma 4, exactly
    rot = l2 + l3
    circ = 1 + l2
    check('case [0, l1): x + (l2+l3) < 1 up to the endpoint, since l1 + l2 + l3 = 1',
          l1 + rot == 1)
    check('case [l1+l2, 1): x + (l2+l3) >= 1 + l2 from the endpoint, and the wrap adds '
          'l3 - 1 = -(l1 + l2)',
          (l1 + l2) + rot == circ and rot - circ == l3 - 1 and l3 - 1 == -(l1 + l2))
    check('case [l1, l1+l2): lands in [1, 1+l2), second step gives x + l3 - l1',
          l1 + rot == 1 and (l1 + l2) + rot == circ and 2 * rot - circ == l3 - l1)
    p_, q_ = sp.symbols('p q')
    L1, L2, L3 = sp.symbols('l1 l2 l3')
    rel = sp.expand(q_ * (L2 + L3) - p_ * (L1 + 2 * L2 + L3))
    check('(l2+l3)/(1+l2) = p/q  <=>  p l1 + (2p - q) l2 + (p - q) l3 = 0',
          sp.expand(rel + (p_ * L1 + (2 * p_ - q_) * L2 + (p_ - q_) * L3)) == 0)
    check('1 + l2 = l1 + 2 l2 + l3 exactly (the tower fills the circle)',
          circ == l1 + 2 * l2 + l3)

    # Lemma 4, numerically
    mp.dps = 50
    L1v, L2v, L3v = (x.value(50) for x in (l1, l2, l3))
    Rv, Cv = L2v + L3v, 1 + L2v

    def tau(x):
        if x < L1v:
            return x + L2v + L3v
        if x < L1v + L2v:
            return x + L3v - L1v
        return x - L1v - L2v

    def induced(x):
        y, n = x, 0
        while True:
            y = (y + Rv) % Cv
            n += 1
            if y < 1:
                return y, n
    random.seed(1)
    worst, times = mpf(0), {1: 0, 2: 0}
    for _ in range(20000):
        x = mpf(random.random())
        y, n = induced(x)
        worst = max(worst, abs(y - tau(x)))
        times[n] += 1
    print(f'  20000 random points: max |tau - induced| = {nstr(worst, 3)}, return times {times}')
    check('Lemma 4 numerically', worst < mpf(10) ** -45 and set(times) == {1, 2})

    # the roof
    words = return_words()
    roofs = []
    for w in words:
        T0_, T1_ = T.path_length(T.unfold(w, M_H))
        check(f'N={len(w):3d}: roof constant (s-coefficient exactly 0), '
              f'r = {trunc(T0_.value(), 12)}', T1_.is_zero())
        roofs.append(T0_)
    vals = sorted({trunc(r.value(), 12) for r in roofs})
    print(f'  distinct roof values: {vals}')
    check('four values, as printed',
          vals == ['19.136181234849', '2.972579301909', '22.108760536759', '27.750141644774'])
    check('2.9 < r < 27.8', all(2.9 < float(r.value()) < 27.8 for r in roofs))
    done()


if __name__ == '__main__':
    main()
