"""Section 5 (and the abstract, Section 6): the census of perpendicular directions.

Every class row (P, Q, eps) -- a centre and one of its two perpendicular direction
classes -- is decided by the exact closure certificate of `lib/closure.py`: every
separatrix of the base B is walked until it arrives at a singularity, arrival decided by
exact equality in Z[zeta_4Q].  capped == 0 PROVES complete periodicity; a capped row is
NO verdict.

Statements checked (odd Q <= 31, the table of Section 5):
  * 212 coprime centres -> 424 classes, all decided; 154 non-Veech centres
    (min(P, Q-P) >= 3) -> 308 classes -> 154 distinct directions on 77 triangles;
    153 certified completely periodic; exactly one direction is not certified, the
    direction of Theorem 1 (8/15 eps 0, and its leg-swap image 7/15 eps 1), which Section 3
    settles; the 58 Veech centres' 116 classes are all certified, none capped.
  * The non-Veech population is phi(Q) - 4 directions at each odd Q >= 7 and none at
    Q = 3, 5: 0 0 2 2 6 8 4 12 14 8 18 16 14 24 26.
  * The class rows of T0 (Section 5, "What the certificate is"; Section 6): the direction
    of Theorem 1 has h = 2 separatrices on the fixed arcs and 2Q - 2h = 26 walked; the
    other perpendicular direction has h = 1 and 28 walked, all of which close, the longest
    after 40 steps.
  * Beyond the table (from the record, see below): odd Q = 33..39 add 88 non-Veech
    directions, 84 certified and 4 undecided (at Q = 33 and 37); over odd Q <= 39,
    242 non-Veech directions, 237 certified, 1 not completely periodic, 4 undecided.
    Even Q = 4..40: 172 directions, 38 Veech (all certified) and 134 non-Veech, of which
    130 are certified and 4 -- 11/26, 11/32, 15/38, 9/40 -- undecided.  No row that
    resolves resolves otherwise than completely periodic.
  * So 2.2(e) holds at every rational right triangle with d = 2Q <= 50 other than T0,
    and the first undecided direction is 11/26 (d = 52); up to d <= 80 eight directions
    are undecided and none fails.
  * 15 | Q: at Q = 30 the 6 non-Veech directions are all certified.

MODES
  python s5_census.py            recompute odd Q <= 31 at cap 3e6 (about an hour) and
                                 check it against the record and the paper
  python s5_census.py --quick    the same for odd Q <= 21 (a few minutes)
  python s5_census.py --record   read every verdict from census_record.json and check
                                 all the figures above without recomputing
  python s5_census.py --rows Q [Q ...] [--cap N]
                                 recompute the given denominators at cap N (default: the
                                 record's cap per row) and compare with the record
census_record.json was produced by the same algorithm in the research repository
(gated identical to lib/closure.py on every class row with Q <= 25).  Its rows past the
3e6 base cap were escalated to caps up to 1e9; to reproduce one, pass that cap.
"""
import json
import os
import sys
from math import gcd

from t0 import check, done

import closure

HERE = os.path.dirname(os.path.abspath(__file__))
RECORD = os.path.join(HERE, 'census_record.json')


def phi(n):
    return sum(1 for k in range(1, n + 1) if gcd(k, n) == 1)


def veech(P, Q):
    return min(P, Q - P) <= 2


def centres(qs):
    return [(P, Q) for Q in qs for P in range(1, Q) if gcd(P, Q) == 1]


def partner(P, Q, e):
    """The leg swap: same direction on the same surface, legs relabelled."""
    return (Q - P, Q, 1 - e) if Q % 2 else (Q - P, Q, e)


def directions(rows, qs):
    """Group class rows into directions (leg-swap pairs).  Returns list of dicts."""
    out, seen = [], set()
    for P, Q in centres(qs):
        for e in (0, 1):
            k = (P, Q, e)
            if k in seen:
                continue
            k2 = partner(*k)
            seen |= {k, k2}
            a, b = rows[k], rows[k2]
            out.append({'key': k, 'pair': k2, 'veech': veech(P, Q),
                        'cp': a['capped'] == 0 or b['capped'] == 0,
                        'agree': (a['capped'] == 0) == (b['capped'] == 0)})
    return out


def load_record():
    d = json.load(open(RECORD))
    return {(r['P'], r['Q'], r['eps']): r for r in d['rows']}


def table31(rows):
    qs = list(range(3, 32, 2))
    C = centres(qs)
    check(f'{len(C)} coprime centres -> {2 * len(C)} classes (expect 212 -> 424)',
          len(C) == 212 and all((P, Q, e) in rows for P, Q in C for e in (0, 1)))
    nv = [(P, Q) for P, Q in C if not veech(P, Q)]
    check(f'{len(nv)} non-Veech centres -> {2 * len(nv)} classes (expect 154 -> 308)',
          len(nv) == 154)
    tri = {frozenset((P, Q - P)) | {Q} for P, Q in nv}
    check(f'{len(tri)} distinct non-Veech triangles (expect 77)', len(tri) == 77)
    D = directions(rows, qs)
    Dn = [d for d in D if not d['veech']]
    check(f'{len(Dn)} non-Veech directions (expect 154)', len(Dn) == 154)
    check('leg-swap partners agree on every direction', all(d['agree'] for d in D))
    cp = [d for d in Dn if d['cp']]
    bad = [d for d in Dn if not d['cp']]
    check(f'{len(cp)} certified completely periodic (expect 153)', len(cp) == 153)
    check(f'not certified: {[d["key"] for d in bad]} -- the direction of Theorem 1 only',
          len(bad) == 1 and {bad[0]['key'], bad[0]['pair']} == {(8, 15, 0), (7, 15, 1)})
    vc = [(P, Q, e) for P, Q in C if veech(P, Q) for e in (0, 1)]
    check(f'{len(vc)} Veech classes, all certified with no capped separatrix (expect 116)',
          len(vc) == 116 and all(rows[k]['capped'] == 0 for k in vc))
    per = [sum(1 for d in Dn if d['key'][1] == Q) for Q in qs]
    check(f'non-Veech directions per Q = {per}',
          per == [0, 0, 2, 2, 6, 8, 4, 12, 14, 8, 18, 16, 14, 24, 26] and sum(per) == 154)
    check('= phi(Q) - 4 at each odd Q >= 7', all(n == phi(Q) - 4 for n, Q in zip(per, qs) if Q >= 7))


def beyond(rows):
    D = directions(rows, range(33, 40, 2))
    Dn = [d for d in D if not d['veech']]
    und = sorted({min(d['key'], d['pair']) for d in Dn if not d['cp']})
    check(f'odd Q = 33..39: {len(Dn)} non-Veech directions, {len(Dn) - len(und)} certified, '
          f'{len(und)} undecided {und}',
          len(Dn) == 88 and len(und) == 4 and {k[1] for k in und} == {33, 37})
    Dall = [d for d in directions(rows, range(3, 40, 2)) if not d['veech']]
    ncp = sum(d['cp'] for d in Dall)
    check(f'odd Q <= 39: {len(Dall)} non-Veech, {ncp} certified, 1 not CP, '
          f'{len(Dall) - ncp - 1} undecided', len(Dall) == 242 and ncp == 237)
    E = directions(rows, range(4, 41, 2))
    Ev = [d for d in E if d['veech']]
    En = [d for d in E if not d['veech']]
    und_e = sorted(min(d['key'], d['pair']) for d in En if not d['cp'])
    names = sorted(f'{min(k[0], k[1] - k[0])}/{k[1]}' for k in und_e)
    check(f'even Q = 4..40: {len(E)} directions, {len(Ev)} Veech (all certified: '
          f'{all(d["cp"] for d in Ev)}), {len(En)} non-Veech, '
          f'{len(En) - len(und_e)} certified, undecided {names}',
          len(E) == 172 and len(Ev) == 38 and all(d['cp'] for d in Ev) and len(En) == 134 and
          sorted(names, key=lambda s: int(s.split('/')[1])) == ['11/26', '11/32', '15/38', '9/40'])
    # 2.2(e) at d = 2Q <= 50, i.e. Q <= 25: every direction decided except T0
    small = [d for d in directions(rows, range(3, 26)) if not d['cp']]
    check('d <= 50: the only undecided or failing direction is T0\'s',
          [sorted((d['key'], d['pair']))[0] for d in small] == [(7, 15, 1)])
    upto80 = [d for d in directions(rows, range(26, 41)) if not d['cp']]
    check(f'26 <= Q <= 40 (d <= 80): {len(upto80)} undecided (expect 8), none failing',
          len(upto80) == 8)
    first = min((d['key'][1] for d in upto80))
    check(f'first undecided at Q = {first} (d = {2 * first}), 11/26', first == 26)
    q30 = [d for d in directions(rows, [30]) if not d['veech']]
    check(f'Q = 30: {len(q30)} non-Veech directions, all certified',
          len(q30) == 6 and all(d['cp'] for d in q30))


def t0_rows():
    a = closure.certify(8, 15, 0, 200_000)
    b = closure.certify(8, 15, 1, 3_000_000)
    check(f'Theorem 1\'s direction (8/15 eps 0): h = {a["h"]}, {a["n_walked"]} walked',
          a['h'] == 2 and a['n_walked'] == 26)
    check(f'the other perpendicular direction (8/15 eps 1): h = {b["h"]}, {b["n_walked"]} '
          f'walked, all close: {b["cp"]}, longest {b["max_steps"]} steps',
          b['h'] == 1 and b['n_walked'] == 28 and b['cp'] and b['max_steps'] == 40)


def recompute(qs, cap=None, rec=None):
    rows = {}
    last = None
    for P, Q in centres(qs):
        if Q != last:
            print(f'  Q = {Q} ...', flush=True)
            last = Q
        for e in (0, 1):
            c = cap if cap is not None else (rec[(P, Q, e)]['cap'] if rec else 3_000_000)
            if (P, Q) in ((8, 15), (7, 15)) and cap is None:
                c = 3_000_000
            r = closure.certify(P, Q, e, c)
            rows[(P, Q, e)] = r
            if rec is not None and (P, Q, e) in rec:
                old = rec[(P, Q, e)]
                same = (r['capped'] == 0) == (old['capped'] == 0) and \
                    (r['capped'] or r['max_steps'] == old['max_steps'])
                if not same:
                    print(f'  DISAGREES with record: {P}/{Q} eps {e}: {r["capped"]} capped '
                          f'(record {old["capped"]} at {old["cap"]})')
                    check.failures += 1
    return rows


def main(argv):
    print('Section 5 -- the census\n')
    rec = load_record()
    t0_rows()
    if '--record' in argv:
        table31(rec)
        beyond(rec)
    elif '--rows' in argv:
        i = argv.index('--rows')
        qs = [int(a) for a in argv[i + 1:] if a.isdigit()]
        cap = int(float(argv[argv.index('--cap') + 1])) if '--cap' in argv else None
        if cap is not None:
            qs = [q for q in qs if q != cap]
        recompute(qs, cap, rec)
    else:
        qmax = 21 if '--quick' in argv else 31
        rows = recompute(range(3, qmax + 1, 2), 3_000_000, rec)
        if qmax == 31:
            table31(rows)
        else:
            D = [d for d in directions(rows, range(3, qmax + 1, 2)) if not d['veech']]
            check(f'odd Q <= {qmax}: every non-Veech direction but T0\'s certified',
                  [d['key'] for d in D if not d['cp']] == [(7, 15, 1)] or
                  sorted(k for d in D if not d['cp'] for k in (d['key'], d['pair'])) ==
                  [(7, 15, 1), (8, 15, 0)])
    done()


if __name__ == '__main__':
    main(sys.argv[1:])
