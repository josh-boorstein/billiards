"""Sections 10.2-10.3: the counts on the lattice locus -- Theorems 10.2, 10.3 and 10.5 -- at every
centre of their rows with Q <= 60 (the paper: "the certificate of Appendix A confirms all four at
every Q <= 60").

Statements checked:
  * (Theorem 10.2) n(1/Q) = 1 and n((Q-1)/Q) = Q - 1 at every 2 <= Q <= 60; at P = 1 the whole
    leg is the orphan branch (J the identity) -- by `lib/partition.py`, Q <= 30;
  * (Theorem 10.3) n((Q-2)/Q) = 2 floor((Q-1)/4) + 1 at every odd 3 <= Q <= 59;
  * (Theorem 10.5) n(2/Q) = Q - 1 at every odd 5 <= Q <= 59.
Every count is the exact flow decomposition's (`lib/flow.py`: the branches are the cylinders'
crossings of L1, certified by their exact tiling of it); the P = 1 orphan is the exact partition's.
"""
import os
import sys

from common import check, done

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'lib'))
from flow import n_flow  # noqa: E402
from partition import partition  # noqa: E402

QMAX = 60


def rows(name, pairs, f):
    bad, untiled = [], []
    for P, Q in pairs:
        r = n_flow(P, Q)
        if not all(r['tiled']):
            untiled.append(f'{P}/{Q}')
        if r['n'] != f(Q):
            bad.append((P, Q, r['n']))
    qs = [Q for _, Q in pairs]
    check(f'{name} at {len(pairs)} centres, Q = {min(qs)}..{max(qs)}: mismatches {bad}'
          + (f', untiled {untiled}' if untiled else ''), not bad and not untiled)


def main():
    print(f'Sections 10.2-10.3 -- the lattice locus, Q <= {QMAX}\n')
    rows('Theorem 10.2: n(1/Q) = 1', [(1, Q) for Q in range(2, QMAX + 1)], lambda Q: 1)
    rows('Theorem 10.2: n((Q-1)/Q) = Q - 1', [(Q - 1, Q) for Q in range(3, QMAX + 1)],
         lambda Q: Q - 1)
    rows('Theorem 10.3: n((Q-2)/Q) = 2 floor((Q-1)/4) + 1',
         [(Q - 2, Q) for Q in range(3, QMAX, 2)], lambda Q: 2 * ((Q - 1) // 4) + 1)
    rows('Theorem 10.5: n(2/Q) = Q - 1', [(2, Q) for Q in range(5, QMAX, 2)], lambda Q: Q - 1)
    ok = 0
    for Q in range(2, 31):
        r = partition(1, Q)
        ok += r['open'] == 0 and r['n'] == 1 and r['orphans'] == [0]
    check(f'Theorem 10.2: at P = 1 the whole leg is the orphan branch, Q = 2..30 (exact partition): '
          f'{ok}/29', ok == 29)
    done()


if __name__ == '__main__':
    main()
