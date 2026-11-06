"""The branch counts the paper quotes at named centres, each recomputed from scratch by two exact
engines: the perpendicular partition (`lib/partition.py`) and the flow decomposition
(`lib/flow.py`).

Statements checked:
  * n(8/15) = 6, n(5/22) = 5, n(5/24) = 17                       (Remark 10.12, Section 13 item 1);
  * n(17/22) = 17, n(19/24) = 7                                  (Section 13 item 1);
  * n(2/3) = 2                                                   (Section 5.1).
For the partition: complete (no interval open), its branches tile the leg exactly, and J is an
exact width-preserving involution on the returning branches.  For the flow: the cylinders tile
every L1 edge exactly.  The counts quoted with Theorems 10.7 and 10.11 and Corollary 10.8
(n(5/11), n(5/6), n(3/4), n_H(3/7)) are in `s10_windows.py`, and Proposition 11.1's in
`s11_cp_band.py`.
"""
import os
import sys

from common import check, done

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'lib'))
from flow import n_flow  # noqa: E402
from partition import partition  # noqa: E402

NAMED = [((8, 15), 6), ((5, 22), 5), ((5, 24), 17), ((17, 22), 17), ((19, 24), 7), ((2, 3), 2)]


def main():
    print('Named-centre branch counts\n')
    for (P, Q), n in NAMED:
        r, f = partition(P, Q), n_flow(P, Q)
        ok = r['open'] == 0 and r['tiles'] and r['J_ok'] and all(f['tiled']) and r['n'] == f['n'] == n
        check(f'n({P}/{Q}) = {r["n"]} (partition) = {f["n"]} (flow)  (paper: {n})', ok)
    done()


if __name__ == '__main__':
    main()
