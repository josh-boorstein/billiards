"""The branch counts and the one word length the paper states at named centres (Section 4.3,
Remark 5.2, Figure 3, Appendix A.3), each recomputed from scratch.

Statements checked:
  * n(8/15) = 6, against Q - 1 = 14                                  (Remark 5.2);
  * n(5/22) = 5 while n(5/24) = 17                                   (Remark 5.2);
  * n(5/9) = 7 (odd, P odd, one orphan) and n(2/9) = 8 (even, no orphan)   (Figure 3);
  * n(3/29) = 19, the plateau example                                (Appendix A.3);
  * the orphan's return word at 5/24 has 311 letters, against 2Q = 48  (Section 4.3): the
    engine stops at the head-on letter, and by Lemma 2.1a the return word is that prefix
    followed by its reversal about the head-on letter, so 2 x prefix - 1 letters.
Each partition is complete (no interval open), its branches tile the leg exactly, and J is an
exact width-preserving involution on the returning branches.

The named counts of the companion necklace paper moved with it (`neck/from_orphan/`).
"""
from common import centre, check, done

NAMED = [((8, 15), 6), ((5, 22), 5), ((5, 24), 17), ((5, 9), 7), ((2, 9), 8), ((3, 29), 19)]


def main():
    print('Named-centre branch counts\n')
    for (P, Q), n in NAMED:
        r = centre(P, Q, fresh=True)
        ok = r['open'] == 0 and r['tiles'] and r['J_ok'] and r['n'] == n
        check(f'n({P}/{Q}) = {r.get("n")}  (paper: {n}); complete, tiling, J an involution', ok)
    check('n(8/15) = 6 < Q - 1 = 14', centre(8, 15, fresh=True)['n'] < 14)
    for (P, Q), k in (((5, 9), 1), ((2, 9), 0)):
        check(f'{P}/{Q}: {k} orphan branch(es)', len(centre(P, Q, fresh=True)['orphans']) == k)
    r = centre(5, 24, fresh=True)
    (i,) = r['orphans']
    length = 2 * r['prefix_len'][i] - 1
    check(f"5/24: the orphan's return word has {length} letters (paper: 311), against 2Q = 48",
          length == 311)
    done()


if __name__ == '__main__':
    main()
