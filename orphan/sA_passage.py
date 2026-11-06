"""Appendix A.1 / Remark A.4: the coefficient growth along the longest passage at 2 alpha / pi
= 11/41, recomputed from scratch (about a minute).

Statements checked:
  * the partition at 11/41 is complete, tiles L1 and has J an involution, and its longest
    passage (prefix to the head-on hit) has 397,153 reflections;
  * every vertex is carried as an integer vector of length 4Q = 164 over the powers of zeta,
    reduced modulo Phi_164 at creation; the largest coefficient of the two vertices of the copy
    that ends that passage is printed, with its ratio to the passage length, and the ratio is
    checked to be below 1 -- the linear growth of Remark A.4 with a small constant.  The same
    bound is checked over the end copies of every returning branch.
The number the paper quotes for the end coefficient is a measurement, not a statement: this
script prints the value the engine gives.
"""
import numpy as np

from common import check, done
from partition import partition


def main():
    print('Appendix A.1 -- the longest passage at 11/41\n')
    r = partition(11, 41)
    check(f"11/41: complete, tiling, J an involution (n = {r['n']})",
          r['open'] == 0 and r['tiles'] and r['J_ok'])
    L = r['leaves']
    m = max(len(x['word']) for x in L)
    check(f'longest passage: {m:,} reflections (paper: 397,153)', m == 397153)
    longest = [x for x in L if len(x['word']) == m]
    end = max(int(np.abs(np.asarray(x[k])).max()) for x in longest for k in ('o_copy', 't_copy'))
    check(f'largest coefficient at its end: {end:,}, ratio to m {end / m:.3f} (< 1, linear growth)',
          end < m)
    worst = max(int(np.abs(np.asarray(x[k])).max()) / len(x['word'])
                for x in L if x['returns'] for k in ('o_copy', 't_copy'))
    check(f'every returning branch: largest end coefficient / passage length <= {worst:.3f} (< 1)',
          worst < 1)
    done()


if __name__ == '__main__':
    main()
