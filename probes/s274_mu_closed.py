#!/usr/bin/env python3
"""
s274_mu_closed.py — PRE-REGISTERED: the odd-`P` closed forms for the cross-section
mass `mu_d`, and a falsification test of the natural closed-form route to odd-`P`
`<N>`.

CONTEXT.  §s270 (2) defines mu_d := sum_sigma len_sigma * Adir_sigma and records
the PARITY RULE for Adir (Adir_L1 = csc h iff Q odd; Adir_H = csc h iff P odd;
Adir_L2 = cot h always), h = pi/2Q, alpha = pi P/2Q.  §s270 (4) assembles it for
EVEN P only, where cot a + csc a = cot(a/2) collapses two terms and gives
    <N> = 2 mu_d - 1 = 2 csc(pi/2Q) + 2 cot(pi P/4Q) cot(pi/2Q) - 1.
This file assembles the OTHER TWO parity classes and scores all three, then tests
the natural route to an odd-P analogue.

H1 (assembly).  With csch/coth = csc/cot(h), csca/cota = csc/cot(alpha):
    P even (forces Q odd):  mu_d = csch + cota*coth + csca*coth = csch + cot(a/2)*coth
    P odd,  Q odd:          mu_d = csch + cota*coth + csca*csch
    P odd,  Q even:         mu_d = coth + cota*coth + csca*csch
  Scored against the MEASURED `M` column of the two s270 stores.  Tolerance 1e-13
  relative (these are the same float sum, so machine precision is expected).

H2 (the falsification).  `<N> = 2 nu(Reach) - 1` is a theorem (§s273 (1)), so an
  odd-P closed form of the same SHAPE as H1 needs nu(Reach) to be the same
  three-term sum over a SUB-COLLECTION of the 2Q directions.  The natural
  candidate sub-collection is a coset of a sublattice: S(d, off) = {k < 2Q :
  k = off mod d}, d | 2Q.  Test every such S against the measured
  nu(Reach) = (<N>+1)/2 on every odd-P row.  A hit means rel err ~1e-13 (the
  tolerance the real identities in this strand hold to); anything at 1e-2 is noise.

PRIOR ART: grepped 'cross-section mass', 'Adir', 'csc(h)', 'cot(h)', 'mu_partial',
'parity lattice', 'sublattice', 'coset', '2Q-1' across all .md ->
  - §s270 (2) (`cf_width_laws.md` 1205-1213) ALREADY OWNS the parity rule and the
    definition; this file does NOT re-book it, only assembles + scores it.
  - §s270 (4) owns the even-P assembly and the cot a + csc a = cot(a/2) step.
  - `rulings.md` [NB] line 147 owns 'g_sw tracks n/Q to 0.644-1.460';
    [NCYL] line 230 owns 'g_sw = n/(Q-1) is DEAD'.  Both are ratio laws for g_sw,
    NOT direction-set descriptions, so H2 is not covered by either.
  - No hit for the P=1 word length 2Q-1, or for any sublattice/coset test.
"""
import json
import math

STORES = ('data/s270_kac_sweep.json', 'data/s270_mean_word_length.json')
SIDES = ('L1', 'L2', 'H')


def mu_closed(P, Q):
    """Closed-form billiard-map cross-section mass (all three parity classes)."""
    h = math.pi / (2 * Q)
    a = 0.5 * math.pi * P / Q
    csch, coth = 1 / math.sin(h), 1 / math.tan(h)
    csca, cota = 1 / math.sin(a), 1 / math.tan(a)
    if P % 2 == 0:                                   # gcd(P,Q)=1 forces Q odd
        return csch + cota * coth + csca * coth      # = csch + cot(a/2)*coth
    if Q % 2 == 1:
        return csch + cota * coth + csca * csch
    return coth + cota * coth + csca * csch


def sub_mass(P, Q, S):
    """Cross-section mass restricted to the direction subset S."""
    a = 0.5 * math.pi * P / Q
    ell = {'L1': 1.0, 'L2': 1 / math.tan(a), 'H': 1 / math.sin(a)}
    nu = {'L1': math.pi, 'L2': 0.5 * math.pi, 'H': a - 0.5 * math.pi}
    return sum(ell[s] * 0.5 * sum(abs(math.cos(k * math.pi / Q - nu[s])) for k in S)
               for s in SIDES)


def load():
    rows = []
    for path in STORES:
        for r in json.load(open(path)):
            rows.append({'P': r['P'], 'Q': r['Q'], 'A': r['A'], 'M': r['M'],
                         'src': path.split('/')[-1]})
    return rows


def h1(rows):
    out, worst = [], {}
    for r in rows:
        P, Q = r['P'], r['Q']
        e = abs(r['M'] - mu_closed(P, Q)) / r['M']
        cls = ('evenP' if P % 2 == 0 else
               ('oddP_oddQ' if Q % 2 else 'oddP_evenQ'))
        out.append({**r, 'mu_closed': mu_closed(P, Q), 'relerr': e, 'cls': cls})
        if e > worst.get(cls, (0, ''))[0]:
            worst[cls] = (e, f'{P}/{Q}')
    return out, worst


def h2(rows):
    out = []
    for r in rows:
        P, Q = r['P'], r['Q']
        if P % 2 == 0:
            continue
        nuR = (r['A'] + 1) / 2
        best = (float('inf'), None, None)
        for d in range(1, 2 * Q):
            if (2 * Q) % d:
                continue
            for off in range(d):
                v = sub_mass(P, Q, range(off, 2 * Q, d))
                e = abs(v - nuR) / nuR
                if e < best[0]:
                    best = (e, d, off)
        out.append({'P': P, 'Q': Q, 'nu_reach': nuR, 'mu': r['M'],
                    'frac': nuR / r['M'], 'best_relerr': best[0],
                    'best_d': best[1], 'best_off': best[2]})
    return out


def main():
    rows = load()
    scored, worst = h1(rows)
    print(f'H1 — mu_d closed form, {len(scored)} rows')
    for cls in sorted(worst):
        print(f'  {cls:>11}: worst rel err {worst[cls][0]:.2e}  ({worst[cls][1]})')
    print(f'  VERDICT: {"PASS" if max(s["relerr"] for s in scored) < 1e-13 else "FAIL"}')

    sub = h2({(r["P"], r["Q"]): r for r in rows}.values())
    hits = [s for s in sub if s['best_relerr'] < 1e-13]
    print(f'\nH2 — is nu(Reach) a direction-coset mass?  {len(sub)} odd-P rows')
    print(f'  best rel err over ALL (d|2Q, offset): '
          f'min {min(s["best_relerr"] for s in sub):.1e}, '
          f'median {sorted(s["best_relerr"] for s in sub)[len(sub)//2]:.1e}')
    print(f'  hits at 1e-13: {len(hits)}  -> VERDICT: '
          f'{"REFUTED (no coset reproduces the swept mass)" if not hits else "HIT"}')

    p1 = [s for s in scored if s['P'] == 1]
    print(f'\nP=1 control ({len(p1)} rows): <N> = 2Q-1 exactly? '
          f'{all(abs(s["A"] - (2 * s["Q"] - 1)) < 1e-9 for s in p1)}'
          f'   => nu(Reach) = Q, vs mu_d ~ (8/pi^2)Q^2')

    json.dump({'h1': scored, 'h2': sub}, open('data/s274_mu_closed.json', 'w'), indent=1)
    print('\nwrote data/s274_mu_closed.json')


if __name__ == '__main__':
    main()
