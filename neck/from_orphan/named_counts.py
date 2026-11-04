"""The branch counts at named centres that moved from the orphan paper to the necklace paper at
the s740 re-cut, each recomputed from scratch.  Statement numbers are the OLD orphan paper's
(until the necklace re-sync re-keys them).

Statements checked:
  * n(17/22) = 17, n(19/24) = 7                                      (old Section 1.2, Section 10 item 2);
  * n(2/3) = 2                                                       (old Section 9.1);
  * at Q = 4 Proposition 5.2's formula reads 1 against a true n(3/4) = 3   (old Section 5.1).
Each partition is complete (no interval open), its branches tile the leg exactly, and J is an
exact width-preserving involution on the returning branches.
"""
from common import centre, check, done

NAMED = [((17, 22), 17), ((19, 24), 7), ((2, 3), 2), ((3, 4), 3)]


def main():
    print('Named-centre branch counts (carried from the orphan paper)\n')
    for (P, Q), n in NAMED:
        r = centre(P, Q, fresh=True)
        ok = r['open'] == 0 and r['tiles'] and r['J_ok'] and r['n'] == n
        check(f'n({P}/{Q}) = {r.get("n")}  (paper: {n}); complete, tiling, J an involution', ok)
    check("Prop 5.2's formula 2 round(Q/3) - 1 reads 1 at Q = 4, against n(3/4) = 3",
          2 * round(4 / 3) - 1 == 1 and centre(3, 4, fresh=True)['n'] == 3)
    done()


if __name__ == '__main__':
    main()
