"""Section 11 and Proposition 11.1: the complete-periodicity verdict on every class row of the band
2 <= Q <= 60, and the two rows where it fails; with Remarks 10.10 and 10.12's band statements.

Statements checked (`band_record.json`, from `lib/flow.py`'s exact flow decomposition; a slice is
recomputed from scratch on every run, all of it with --full):
  * (Section 11) all 2202 class rows -- both classes at each of the 1101 coprime centres, 1/2
    included -- are decided: the class eps = 0 of P/Q by the direction perpendicular to L1; at odd
    Q the class eps = 1 of P/Q as the class eps = 0 of (Q-P)/Q; at even Q the class eps = 1 by the
    direction perpendicular to H.  2200 are completely periodic (every component a cylinder, the
    cylinders tiling the beam's side exactly), and the 2 others are 8/15 eps = 0 and 7/15 eps = 1;
  * (Section 11) at 8/15 the direction perpendicular to L1 has 6 cylinders and 2 components
    without periodic trajectory, each with a Boshernitzan certificate and positive mass; the other
    class is completely periodic -- perpendicular to L2 and to OA (one direction at odd Q) every
    component is a cylinder, the same ones, and that class is the L1 class at 7/15;
  * (Section 11) perpendicular to H at every even Q <= 60 (373 centres) the class is completely
    periodic and n_H = Q, as Corollary 7.11 requires; at odd Q <= 31 the H decomposition agrees with
    the verdict of the paired class on all 212 rows;
  * (Section 11) the exact closure certificate (`neck/census_record.json`, the
    minimal-component paper's census, `lib/closure.py`) gives the same verdicts on the 424 class
    rows with odd Q <= 31: complete periodicity on all but the same two, which it leaves capped;
  * (Proposition 11.1) n(8/15) = 6 by the flow and by `lib/partition.py`; and the second witness,
    n_L2(8/15) + n_H(8/15) = 11 + 5 = 16 = Q + 1, so D = 0 on the other class (Corollary 7.8);
  * (Remark 10.10) every P = 3 centre with Q <= 60 is completely periodic, and so is 3/385, with
    n = 255 = 2 round(385/3) - 1, recomputed on every run (about a minute);
  * (Remark 10.12) among the even-P centres off the lattice locus with Q <= 60 (P not in
    {2, Q - 2, Q - 1}, i.e. min(P, Q - P) > 2) there are 307, and 8/15 is the only one whose count is not maximal;
  * control: at nine other odd-Q centres (Q = 13, 15, 17) every perpendicular direction is
    completely periodic, so the verdict at 8/15 is not the instrument's default.
"""
import json
import math
import os
import sys
import time

from common import FULL, check, done

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'lib'))
import flow  # noqa: E402
from necklace import interval_model  # noqa: E402
from partition import partition  # noqa: E402

RECORD = os.path.join(HERE, 'band_record.json')
CENSUS = os.path.join(HERE, 'census_record.json')
QMAX = 60
SLICE = 60 if FULL else 20
CONTROLS = [(2, 13), (4, 13), (6, 13), (8, 13), (2, 15), (4, 15), (6, 17), (8, 17), (10, 17)]


def cp(r):
    return r['ncyl'] == r['ncomp'] and r['tiled']


def perp(P, Q):
    """The decomposition perpendicular to each side of the triangle of P/Q (sides 0, 1, 2:
    OA = H, AR = L1, RO = L2)."""
    abc = (P, Q - P, Q)
    return [flow.decompose(abc, j, cap=20000) for j in flow.perp_dirs(abc)]


def cylinders_only(r):
    return set(r['counts']) == {'cylinder'}


def periods(r):
    return sorted(c['period'] for c in r['components'])


def main():
    rec = json.load(open(RECORD))
    L1, H = rec['L1'], rec['H']
    centres = [(P, Q) for Q in range(2, QMAX + 1) for P in range(1, Q) if math.gcd(P, Q) == 1]
    print(f'Section 11 -- {len(centres)} centres, Q <= {QMAX}\n')

    # the record against a fresh computation, on the slice
    sl = [(P, Q) for P, Q in centres if Q <= SLICE]
    same = 0
    for P, Q in sl:
        r = flow.n_flow(P, Q)
        o = L1[f'{P}/{Q}']
        same += (r['n'], r['all_cyl'], all(r['tiled'])) == (o['n'], o['ncyl'] == o['ncomp'], o['tiled'])
    slh = [k for k in H if int(k.split('/')[1]) <= SLICE]
    sameh = 0
    for k in slh:
        P, Q = map(int, k.split('/'))
        r = flow.n_flow(P, Q, beam='h')
        sameh += (r['n'], r['all_cyl'], all(r['tiled'])) == (H[k]['n'], cp(H[k]), H[k]['tiled'])
    check(f'record recomputed from scratch, Q <= {SLICE}: perpendicular to L1 {same}/{len(sl)}, '
          f'to H {sameh}/{len(slh)}', same == len(sl) and sameh == len(slh))

    # the class rows
    assert all(f'{P}/{Q}' in L1 for P, Q in centres) and len(L1) == len(centres)
    verdict = {}
    for P, Q in centres:
        verdict[(P, Q, 0)] = cp(L1[f'{P}/{Q}'])
        if Q % 2:
            verdict[(P, Q, 1)] = cp(L1[f'{Q - P}/{Q}'])
        else:
            verdict[(P, Q, 1)] = cp(H[f'{P}/{Q}'])
    non = sorted(k for k, v in verdict.items() if not v)
    check(f'{len(verdict)} class rows, {sum(verdict.values())} completely periodic, the others {non}',
          len(verdict) == 2202 and sum(verdict.values()) == 2200 and non == [(7, 15, 1), (8, 15, 0)])
    untiled = [k for k, v in L1.items() if not v['tiled']] + [k for k, v in H.items() if not v['tiled']]
    check(f'every decomposition tiles the beam\'s side exactly (branch counts certified): untiled '
          f'{untiled}', not untiled)

    oa, l1, l2 = perp(8, 15)
    free = [c for c in l1['components'] if c['kind'] == 'no_periodic']
    ncyl = l1['counts'].get('cylinder', 0)
    check(f'8/15 perpendicular to L1: {ncyl} cylinders and {len(free)} components without periodic '
          f'trajectory (paper: 6 and 2)', ncyl == 6 and len(free) == 2
          and set(l1['counts']) <= {'cylinder', 'no_periodic'})
    masses = ', '.join(f'{c["mass"]:.4f}' for c in free)
    check(f'  each with a Boshernitzan certificate and positive mass ({masses})',
          all(c['cert'] is not None and c['mass'] > 0 for c in free))
    check(f'8/15 perpendicular to L2: all {l2["counts"].get("cylinder", 0)} components cylinders; '
          f'perpendicular to OA the same ones', cylinders_only(l2) and cylinders_only(oa)
          and periods(oa) == periods(l2))
    _, l1_7, _ = perp(7, 15)
    check('7/15 perpendicular to L1 (= the L2 class at 8/15): all cylinders', cylinders_only(l1_7))

    ev = [k for k in H if int(k.split('/')[1]) % 2 == 0]
    good = [k for k in ev if cp(H[k]) and H[k]['n'] == int(k.split('/')[1])]
    check(f'perpendicular to H at even Q <= {QMAX}: {len(good)}/{len(ev)} completely periodic with '
          f'n_H = Q (paper: 373)', len(ev) == 373 and len(good) == 373)
    od = [k for k in H if int(k.split('/')[1]) % 2]
    agree = 0
    for k in od:
        P, Q = map(int, k.split('/'))
        lone0 = interval_model(P, Q, 0)[0]
        lone0 = lone0.ell[lone0.m_idx].is_zero()
        paired = 1 if lone0 else 0                 # the class whose interface m is vertical
        agree += cp(H[k]) == verdict[(P, Q, paired)]
    check(f'perpendicular to H at odd Q <= 31: agrees with the paired class on {agree}/{len(od)} '
          f'(paper: 212)', len(od) == 212 and agree == 212)

    cen = [r for r in json.load(open(CENSUS))['rows'] if r['Q'] % 2 and r['Q'] <= 31]
    same = sum((r['capped'] == 0) == verdict[(r['P'], r['Q'], r['eps'])] for r in cen)
    check(f'closure certificate, odd Q <= 31: the same verdict on {same}/{len(cen)} class rows '
          f'(paper: 424)', len(cen) == 424 and same == 424)

    nf = flow.n_flow(8, 15)
    check(f'Proposition 11.1: n(8/15) = {nf["n"]} by the flow, tiled, not CP; '
          f'= {partition(8, 15)["n"]} by the partition', nf['n'] == 6 == partition(8, 15)['n']
          and all(nf['tiled']) and not nf['all_cyl'])
    nl2, nh = flow.n_flow(7, 15)['n'], flow.n_flow(8, 15, beam='h')['n']
    check(f'Proposition 11.1, second witness: n_L2(8/15) + n_H(8/15) = {nl2} + {nh} = Q + 1',
          (nl2, nh) == (11, 5))

    p3 = [(P, Q) for P, Q in centres if P == 3]
    check(f'Remark 10.10: every P = 3 centre with Q <= {QMAX} completely periodic: '
          f'{sum(cp(L1[f"3/{Q}"]) for _, Q in p3)}/{len(p3)}', all(cp(L1[f'3/{Q}']) for _, Q in p3))
    t = time.time()
    r = flow.n_flow(3, 385)
    check(f'Remark 10.10: 3/385 completely periodic, n = {r["n"]} = 2 round(385/3) - 1, tiled '
          f'({time.time() - t:.0f} s)', r['all_cyl'] and all(r['tiled']) and r['n'] == 255)

    off = [(P, Q) for P, Q in centres if P % 2 == 0 and min(P, Q - P) > 2]
    low = [f'{P}/{Q}' for P, Q in off if L1[f'{P}/{Q}']['n'] != Q - 1]
    check(f'Remark 10.12: {len(off)} even-P centres off the lattice locus, Q <= {QMAX} (paper: 307); '
          f'count not maximal at {low}', len(off) == 307 and low == ['8/15'])

    bad = [f'{P}/{Q}' for P, Q in CONTROLS if not all(cylinders_only(r) for r in perp(P, Q))]
    check(f'control: every perpendicular direction completely periodic at {len(CONTROLS)} other '
          f'odd-Q centres' + (f' (fails at {", ".join(bad)})' if bad else ''), not bad)
    done()


if __name__ == '__main__':
    main()
