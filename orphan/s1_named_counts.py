"""The branch counts the paper states at named centres (Sections 1.2, 5, 5.2, 9.1, 10 and
Appendix A.3), each recomputed from scratch.

Statements checked:
  * n(5/22) = 5, n(5/24) = 17, n(17/22) = 17, n(19/24) = 7          (Section 1.2, Remark 5.4,
    Section 10 item 2);
  * n(8/15) = 6, not Q - 1 = 14                                     (Section 1.2, Remark 5.4);
  * n(5/9) = 7 (odd, P odd) and n(2/9) = 8 (even, P even)           (Figure 3);
  * n(2/3) = 2                                                      (Section 9.1);
  * at Q = 4 Proposition 5.2's formula reads 1 against a true n(3/4) = 3   (Section 5.1);
  * n(3/29) = 19, the plateau example of Appendix A.3.
Each partition is complete (no interval open), its branches tile the leg exactly, and J is an
exact width-preserving involution on the returning branches.
"""
from common import centre, check, done

NAMED = [((5, 22), 5), ((5, 24), 17), ((17, 22), 17), ((19, 24), 7), ((8, 15), 6),
         ((5, 9), 7), ((2, 9), 8), ((2, 3), 2), ((3, 4), 3), ((3, 29), 19)]


def main():
    print('Named-centre branch counts\n')
    for (P, Q), n in NAMED:
        r = centre(P, Q, fresh=True)
        ok = r['open'] == 0 and r['tiles'] and r['J_ok'] and r['n'] == n
        check(f'n({P}/{Q}) = {r.get("n")}  (paper: {n}); complete, tiling, J an involution', ok)
    check('n(8/15) = 6 < Q - 1 = 14', centre(8, 15, fresh=True)['n'] < 14)
    check("Prop 5.2's formula 2 round(Q/3) - 1 reads 1 at Q = 4, against n(3/4) = 3",
          2 * round(4 / 3) - 1 == 1 and centre(3, 4, fresh=True)['n'] == 3)
    done()


if __name__ == '__main__':
    main()
