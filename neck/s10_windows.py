"""Section 10: the bounded-walk families, certified for all large Q, and their finite rows.

Statements checked:
  * the prover reproduces the two hand proofs as its controls: P = 3 far arc, r = 2 -> 1 cell,
    r = 1 -> 3 cells for Q >= 10 (Theorem 10.2 at rho = 1); P = 5, r = 1, near arc -> 5 cells
    for Q >= 14 (Theorem 10.6);
  * (Theorem 10.7) P = 5, r = 4, far: 3 window cells, 114 inequalities, root-free for Q >= 28;
    P = 7, r = 6, far: 3 window cells, 90 inequalities, root-free for Q >= 20;
  * the prover refuses a family whose walk is unbounded (P = 5, r = 2): nothing is claimed;
  * the families' conclusions, n(5, 5m+1) = 5 (Q >= 14), n(5, 5m+4) = 4m + 1 (Q >= 28),
    n(7, 7m+6) = 6m + 3 (Q >= 20), hold at every such centre with Q <= 60 (--full: Q <= 100)
    by an independent engine, the exact flow decomposition (`lib/flow.py`: the branches are the
    cylinders' crossings of L1, certified by their exact tiling of it);
  * the finite rows the theorems do not cover, by exact computation: n(5/11) = 7, n(5/6) = 5,
    n_H(3/7) = 5 (the beam perpendicular to the hypotenuse), and the rows Q = 24 at P = 5 and
    Q = 13, 20 at P = 7, whose certified counts agree with 4m + 1 and 6m + 3 -- each by the
    flow decomposition AND by the exact perpendicular partition (`lib/partition.py`).
"""
from common import FULL, check, done
from flow import n_flow
from partition import partition
from window_prover import prove

QMAX = 100 if FULL else 60


def main():
    print('Section 10 -- bounded walks\n')
    want = {(3, 2, 'far'): (1, None, None), (3, 1, 'far'): (3, 10, None),
            (5, 1, 'near'): (5, 14, None), (5, 4, 'far'): (3, 28, 114),
            (7, 6, 'far'): (3, 20, 90)}
    for (P, r, side), (cells, q0, nchk) in want.items():
        out = prove(P, r, side, verbose=False)
        ok = out is not None and out[0] == cells and (q0 is None or out[1] == q0) \
            and (nchk is None or out[3] == nchk)
        got = None if out is None else f'{out[0]} cells, {out[3]} checks, Q >= {out[1]}'
        check(f'({P}, {r}, {side}): {got}', ok)
    check('(5, 2, near): unbounded, nothing claimed', prove(5, 2, 'near', verbose=False) is None)
    fam = [(5, 1, 14, lambda m: 5), (5, 4, 28, lambda m: 4 * m + 1), (7, 6, 20, lambda m: 6 * m + 3)]
    for P, r, q0, f in fam:
        qs = [Q for Q in range(q0, QMAX + 1) if Q % P == r]
        bad = [(Q, n) for Q in qs if (n := n_flow(P, Q)['n']) != f(Q // P)]
        check(f'n({P}, {P}m+{r}) at every Q in [{q0}, {QMAX}] ({len(qs)} centres), flow '
              f'decomposition: mismatches {bad}', not bad)
    fin = [((5, 11), 7), ((5, 6), 5), ((5, 24), 17), ((7, 13), 9), ((7, 20), 15)]
    for (P, Q), n in fin:
        a, b = n_flow(P, Q)['n'], partition(P, Q)['n']
        check(f'n({P}/{Q}) = {a} (flow) = {b} (partition)  (paper: {n})', a == b == n)
    a, b = n_flow(3, 7, beam='h')['n'], partition(3, 7, beam='h')['n']
    check(f'n_H(3/7) = {a} (flow) = {b} (partition)  (paper: 5)', a == b == 5)
    done()


if __name__ == '__main__':
    main()
