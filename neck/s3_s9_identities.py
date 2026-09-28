"""Sections 3 and 9.1: the counting identities on the separatrix diagram, at every class row of
the census (odd Q <= 39, even Q <= 40).

Statements checked (h, v, f, D, E, the swapped edges: `lib/necklace.py`):
  * (Lemma 3.2a) no resolved saddle connection meets the vertical arcs of Fix(iota) more than
    once -- at every odd-Q row;
  * (Theorem 3.1) the critical graph of a completely periodic class has E = Q edges;
  * (Theorem 3.3) on a completely periodic class Q = h + f + 2p, i.e. D = Q - h - f equals
    twice the number of iota-swapped edge pairs -- and it is 0 at every CP row here, so the
    covering law (Corollary 3.4) holds on all of them;
  * (Corollary 9.8) D = 0 certifies complete periodicity: every row with D = 0 is a CP row;
  * (after Corollary 9.8) the swapped-edge detector 2p = D - s/2, s the non-closing rays (the
    capped walks that never met a vertical arc -- one that crossed before capping is already in
    f, so is not counted again: `cap_conf`), is cutoff-free and fires on exactly the two rows 8/15 eps 0 and 7/15 eps 1,
    where D = 8 = 2p + s/2 = 2 + 6 at 8/15 eps 0;
  * (after Theorem 9.8e) at even Q every completely periodic eps = 0 row has C = Q/2 cylinders,
    C = (Q - b)/2 + kappa - 1 read off the critical graph (Theorem 3.1).
A CP row is walked at exactly the cap that closes it; a capped row at 3e6 (D is then an upper
bound, and the detector is cutoff-free).  The default run takes the rows needing <= 3e6 steps;
--full takes every row.
"""
from common import FULL, QUICK_CAP, census_rows, check, done, row_readout


def main():
    rows = [r for r in census_rows() if FULL or r['capped'] or r['max_steps'] <= QUICK_CAP]
    print(f'Sections 3 and 9.1 -- {len(rows)} class rows\n')
    odd = [r for r in rows if r['Q'] % 2]
    lemma = E_ok = ident = d0_cp = n_cp_odd = n_d0 = 0
    fires, even_cp0, c_ok, n_cp = [], 0, 0, 0
    at815 = None
    for r in rows:
        o = row_readout(r)
        cp = o['capped'] == 0
        if cp:
            n_cp += 1
            E_ok += o['E'] == r['Q']
        if r['Q'] % 2:
            lemma += o['multi'] == 0
            if cp:
                n_cp_odd += 1
                ident += o['D'] == 2 * o['swapped'] == 0
            if o['D'] == 0:
                n_d0 += 1
                d0_cp += cp
            det = o['D'] - o['cap_conf'] // 2
            if det != 0:
                fires.append((r['P'], r['Q'], r['eps'], det, 2 * o['swapped']))
            if (r['P'], r['Q'], r['eps']) == (8, 15, 0):
                at815 = (o['D'], o['cap_conf'] // 2, det)
        elif cp and r['eps'] == 0:
            even_cp0 += 1
            c_ok += o['C'] == r['Q'] // 2
    check(f'Lemma 3.2a: no resolved edge meets the vertical arcs twice: {lemma}/{len(odd)} odd rows',
          lemma == len(odd))
    check(f'Theorem 3.1: E = Q on every CP row: {E_ok}/{n_cp}', E_ok == n_cp)
    check(f'Theorem 3.3 / Cor 3.4: D = 2p = 0 on every odd CP row: {ident}/{n_cp_odd}',
          ident == n_cp_odd)
    check(f'Corollary 9.8: every row with D = 0 is CP: {d0_cp}/{n_d0}', d0_cp == n_d0)
    check(f'the detector D - s/2 fires only at 8/15 eps 0 and 7/15 eps 1, value 2 = twice the '
          f'resolved swapped pairs (p = 1): {fires}',
          sorted(fires) == [(7, 15, 1, 2, 2), (8, 15, 0, 2, 2)])
    check(f'8/15 eps 0: D = 8 = 2p + s/2 = 2 + 6: (D, s/2, 2p) = {at815}', at815 == (8, 6, 2))
    check(f'even Q: C = Q/2 on every CP eps = 0 row: {c_ok}/{even_cp0}', c_ok == even_cp0)
    done()


if __name__ == '__main__':
    main()
