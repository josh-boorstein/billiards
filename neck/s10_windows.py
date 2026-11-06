"""Section 10.4: the bounded walks of Theorems 10.8 and 10.12, certified for all large Q, and the
counts they give at every centre of the band.

Statements checked:
  * the prover (`window_prover.py`) certifies the two walks: P = 3, far arc, r3 = 2 -> 1 window
    cell and r3 = 1 -> 3 window cells for Q >= 10 (Theorem 10.8); P = 5, r = 1, near arc -> 5
    cells for Q >= 14 (Theorem 10.12) -- the printed conditions are `sA3_conditions.py`;
  * the prover refuses a family whose walk is unbounded (P = 5, r = 2): nothing is claimed;
  * Theorem 10.8's counts, n_H(3/Q) = (Q + 4)/3 or (Q + 8)/3 at odd Q and
    n((Q-3)/Q) = n_L2(3/Q) = (Q + 1)/3 or (Q + 5)/3 at even Q (r3 = 2, resp. 1), at every
    3 does-not-divide Q in [5, 60] (--full: 100), and Corollary 10.9's count n(3/Q) = 2 round(Q/3) - 1
    (equality: every such class is completely periodic, Remark 10.10a) with n(3/4) = 3 -- by an
    independent engine, the exact flow decomposition (`lib/flow.py`: the branches are the cylinders'
    crossings of the beam's side, certified by their exact tiling of it);
  * Theorem 10.12's count n(5/Q) = 5 at every Q = 5m + 1, m >= 3, Q <= 60 (--full: 100), by the
    flow decomposition, and its exceptional rows n(5/11) = 7, n(5/6) = 5;
  * the row Q = 7 of Theorem 10.8, n_H(3/7) = 5 (the beam perpendicular to the hypotenuse);
    the finite rows by the flow decomposition AND by the exact perpendicular partition
    (`lib/partition.py`).
"""
from common import FULL, check, done
from flow import n_flow
from partition import partition
from window_prover import prove

QMAX = 100 if FULL else 60


def main():
    print('Section 10.4 -- bounded walks\n')
    want = {(3, 2, 'far'): (1, None), (3, 1, 'far'): (3, 10), (5, 1, 'near'): (5, 14)}
    for (P, r, side), (cells, q0) in want.items():
        out = prove(P, r, side, verbose=False)
        ok = out is not None and out[0] == cells and (q0 is None or out[1] == q0)
        got = None if out is None else f'{out[0]} cells, {out[3]} checks, Q >= {out[1]}'
        check(f'({P}, {r}, {side}): {got}', ok)
    check('(5, 2, near): unbounded, nothing claimed', prove(5, 2, 'near', verbose=False) is None)

    qs = [Q for Q in range(5, QMAX + 1) if Q % 3]
    bad = []
    for Q in qs:
        r3 = Q % 3
        if Q % 2:
            got, want_ = n_flow(3, Q, beam='h')['n'], (Q + 4) // 3 if r3 == 2 else (Q + 8) // 3
        else:
            got, want_ = n_flow(Q - 3, Q)['n'], (Q + 1) // 3 if r3 == 2 else (Q + 5) // 3
        if got != want_:
            bad.append((Q, got, want_))
    check(f'Theorem 10.8: n_H(3/Q) at odd Q, n((Q-3)/Q) at even Q, every Q in [5, {QMAX}], '
          f'3 not dividing Q ({len(qs)} centres): mismatches {bad}', not bad)
    bad = []
    for Q in qs:
        f = n_flow(3, Q)
        if not (f['all_cyl'] and f['n'] == 2 * round(Q / 3) - 1):
            bad.append((Q, f['n'], f['all_cyl']))
    check(f'Corollary 10.9: n(3/Q) = 2 round(Q/3) - 1, the class completely periodic, at the same '
          f'{len(qs)} centres: exceptions {bad}', not bad)
    a, b = n_flow(3, 4)['n'], partition(3, 4)['n']
    check(f'Corollary 10.9 at Q = 4: n(3/4) = {a} (flow) = {b} (partition), the formula reads 1',
          a == b == 3 and 2 * round(4 / 3) - 1 == 1)

    qs = [Q for Q in range(16, QMAX + 1) if Q % 5 == 1]
    bad = [(Q, n) for Q in qs if (n := n_flow(5, Q)['n']) != 5]
    check(f'Theorem 10.12: n(5/Q) = 5 at every Q = 5m + 1 in [16, {QMAX}] ({len(qs)} centres), '
          f'flow decomposition: mismatches {bad}', not bad)
    for (P, Q), n in [((5, 11), 7), ((5, 6), 5)]:
        a, b = n_flow(P, Q)['n'], partition(P, Q)['n']
        check(f'n({P}/{Q}) = {a} (flow) = {b} (partition)  (paper: {n})', a == b == n)
    a, b = n_flow(3, 7, beam='h')['n'], partition(3, 7, beam='h')['n']
    check(f'Theorem 10.8 at Q = 7: n_H(3/7) = {a} (flow) = {b} (partition)  (paper: 5)', a == b == 5)
    done()

if __name__ == '__main__':
    main()
