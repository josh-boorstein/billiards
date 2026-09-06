"""s377 PRE-REGISTERED -- is (V2) a RESIDUAL for [NCYL-216]'s converse, or is it the
converse restated?  And: what ARE the `R`-copies on the surface?

WHY.  s374 named (R★) as the converse's residual; s375 showed (R★) is FALSE at a regular
non-`R` passage and named (V2) instead -- *the limit orbit at a `J`-fold meets no CONE
point* -- measured `689/689`, and booked the converse as PROVED MODULO (V2).  Two
residuals in two sessions is the [OPS-039] structured-failure tell, so this probe audits
the INSTRUMENT and the LOGIC instead of measuring a third one.

  ⇒ PART A, THE AUDIT.  A cell boundary is BY DEFINITION a height whose ray grazes a
  vertex (Lemma A / `foundations.md` §1), so the limit orbit at `y` passes through its own
  graze vertex `V(y)`.  And for `P ≥ 3` and `r = Q−P ≥ 3` every vertex except `R` is a
  cone point (`O`: `2πP`; `A`: `2π(Q−P)` or `π(Q−P)` -- `veech_reframing.md` §1).  If both
  hold then

        (V2)  ⟹  V(y) is not a cone point  ⟹  V(y) = R,

  which is Proposition 2's CONCLUSION, reached without any of §B′.3 (1a)-(1e).  The
  hypothesis would then imply the conclusion and (V2) would be the converse restated, not
  a reduction of it -- and the `689/689` scoring it could not have come out otherwise
  ([OPS-041]), because it was scored on FOLDS and every fold in the sample is an
  `R`-boundary.  A1/A2/A3 score exactly that, off s375's own stored rows (no tracing).

  ⚠ NAMED EXPECTED-WRONG READING, and it is the one that would SAVE (V2): that the graze
  vertex is NOT always met by the limit orbit -- e.g. that `graze_hi_V` records a vertex
  the ray only approaches.  A1 is the test.  If A1 fails, the audit collapses.

  ⇒ PART B, THE STRUCTURE.  Weakening (V2) to what §B′.3 (1d) actually uses -- `γ⁻ = γ⁺`
  -- does not escape A3 either: at a cone point of angle `2πk` the two one-sided limits
  leave along the first and last of the `k` outgoing separatrices, distinct iff `k ≥ 2`.
  So the honest question is what a self-reverse first-return arc IS on the surface.  Its
  midpoint is fixed by the direction-reversing involution `τ` (`τ*ω = −ω`, the
  hyperelliptic involution -- `veech_reframing.md` §3e), i.e. it is a WEIERSTRASS point.
  B scores the resulting identification against the independent genus count: the `Q`
  copies of `R` are `τ`-fixed (they lift the `Q` simple poles of the genus-0 base
  `Q(P−2, Q−P−2, −1^Q)`, all of ODD order), `O` is `τ`-fixed iff `P−2` is odd and `A` iff
  `Q−P−2` is odd, and the total must be exactly `2g+2 = 2⌊Q/2⌋+2`.

PRIOR ART: grepped `rulings.md` + the ledgers for 'saddle loop', 'self-reverse',
'Weierstrass', 'hyperelliptic', 'involution', 'cone point', 'V2' ->
[NCYL-220] (s375's "PROVED modulo (V2)" -- the claim audited here), [NCYL-216]/[NCYL-217]
(the forward direction, PROVED, untouched by this), [NCYL-133] (`R` is the unique
reversing vertex -- Part B reframes it), [T-WEIER] + `veech_reframing.md` §3e (`O` is a
Weierstrass point ⟺ `P` odd; the deck involution IS hyperelliptic) and §1/§1a (the cone
table + `g = ⌊Q/2⌋`, PROVED), [RENORM-005] (`τ` is invisible on a one-orientation section
-- why the deck never showed up in the IET work), [OPS-041] (a denominator that cannot
fail), [OPS-071] (scoring a law against a column derived from it).
⚠ Nothing anywhere identifies the `R`-copies as Weierstrass points; `veech_reframing.md`
§1 records only that they are REGULAR ("the right angle is invisible").  Part B is new.

USAGE
    PYTHONPATH=.:engine:archive/scripts_2026-07:archive/scripts_2026-08:probes \
      .venv/bin/python3.13 probes/s377_residual_audit.py [--qmax 200]
"""
import argparse
import collections
import json
import math


SRC = 'data/s375_foldmid.json'


def regular(V, P, Q):
    """Is the vertex copy `V` a REGULAR point of the unfolding (cone angle 2π)?

    `veech_reframing.md` §1/§1a: cone angles are `2πP` at `O`, `2π` at each of the `Q`
    copies of `R`, and `2π(Q−P)` (Q even) / `π(Q−P)` (Q odd) at `A`.
    """
    if V == 'R':
        return True
    if V == 'O':
        return P in (1, 2)
    if V == 'A':
        return (Q - P) in (1, 2)
    raise ValueError(V)


def part_a():
    rows = json.load(open(SRC))['rows']
    c = collections.Counter()
    off = []
    for r in rows:
        P, Q, V = r['P'], r['Q'], r['V']
        clus = {x[0] for x in r['clus']}
        nocone = not r['conepts']
        c['rows'] += 1
        # A1  the limit orbit meets its own graze vertex
        c['A1'] += int(V in clus)
        # A2  (V2) holds  <=>  the graze vertex is REGULAR
        c['A2'] += int(nocone == regular(V, P, Q))
        if nocone != regular(V, P, Q):
            off.append((P, Q, V, r['clus'], r['conepts']))
        # A3  on the generic regime, (V2) <=> the graze is at R  [= the conclusion]
        if P >= 3 and Q - P >= 3:
            c['A3_n'] += 1
            c['A3'] += int(nocone == (V == 'R'))
        # A4  the [OPS-041] denominator: folds with a non-R graze
        if r['fold']:
            c['folds'] += 1
            c['folds_nonR'] += int(V != 'R')
    print('=== PART A -- is (V2) the converse restated? ===')
    print(f'  rows scored                                        {c["rows"]}')
    print(f'  A1  limit orbit meets its own graze vertex         {c["A1"]}/{c["rows"]}')
    print(f'  A2  (V2)  <=>  graze vertex REGULAR                {c["A2"]}/{c["rows"]}')
    print(f'  A3  (V2)  <=>  graze at R   [P>=3, r>=3]           {c["A3"]}/{c["A3_n"]}')
    print(f'  A4  folds in the sample                            {c["folds"]}')
    print(f'      of them with a NON-R graze  [the (V2) denominator]  '
          f'{c["folds_nonR"]}')
    if off:
        print(f'  ⚠ A2 exceptions ({len(off)}):')
        for t in off[:12]:
            print('     ', t)
    ok = (c['A1'] == c['rows'] and c['A2'] == c['rows'] and c['A3'] == c['A3_n'])
    print(f'  ⇒ (V2) is EQUIVALENT to the conclusion on P>=3,r>=3: {ok}; '
          f'and the {c["folds"]} rows scoring it contain {c["folds_nonR"]} '
          f'that could have failed.')
    # ASSERT THE DIAGNOSIS ([OPS-043]) -- if any of these ever stops holding, (V2) is
    # a genuine residual after all and this section must be re-opened.
    assert c['A1'] == c['rows'], 'A1: a graze vertex is NOT met -- audit collapses'
    assert c['A2'] == c['rows'], 'A2: (V2) is not coextensive with graze-regularity'
    assert c['A3'] == c['A3_n'], 'A3: (V2) is not the conclusion on the generic regime'
    assert c['folds_nonR'] == 0, 'A4: a fold with a non-R graze exists -- converse FALSE'
    return c


def part_b(qmax):
    """`τ`-fixed points (Weierstrass) = lifts of the ODD-order base singularities.

    Base (Apisa, `veech_reframing.md` §3a) is the genus-0 `Q(P−2, Q−P−2, −1^Q)`:
    orders `P−2` at O, `Q−P−2` at A, and `Q` simple poles (order `−1`, ODD) which lift to
    the `Q` REGULAR points -- the copies of `R`.  Weierstrass count must be `2g+2`.
    """
    bad, n, splits = [], 0, collections.Counter()
    for Q in range(3, qmax + 1):
        for P in range(1, Q):
            if math.gcd(P, Q) != 1:
                continue
            n += 1
            fixed_R = Q                       # the Q simple poles, order -1, odd
            fixed_O = int((P - 2) % 2 != 0)   # P odd
            fixed_A = int((Q - P - 2) % 2 != 0)
            g = Q // 2
            if fixed_R + fixed_O + fixed_A != 2 * g + 2:
                bad.append((P, Q, fixed_R, fixed_O, fixed_A, 2 * g + 2))
            splits[(P % 2, Q % 2, fixed_O, fixed_A)] += 1
    print('\n=== PART B -- the R-copies are Weierstrass points ===')
    print(f'  coprime (P,Q), 3 <= Q <= {qmax}:                     {n}')
    print(f'  Q + [P odd] + [Q−P odd]  ==  2g+2 = 2⌊Q/2⌋+2:      {n - len(bad)}/{n}')
    print('  (P%2, Q%2) -> O fixed, A fixed:')
    for k in sorted(splits):
        print(f'     P{"odd" if k[0] else "even":>4}  Q{"odd" if k[1] else "even":>5}'
              f'   O {bool(k[2])!s:>5}   A {bool(k[3])!s:>5}   ({splits[k]} centres)')
    if bad:
        print(f'  ⚠ {len(bad)} MISMATCHES: {bad[:8]}')
    print('  ⇒ the Q copies of R exhaust the Weierstrass points other than O and A.')
    assert not bad, 'B: Weierstrass count does not match 2g+2'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--qmax', type=int, default=200)
    a = ap.parse_args()
    part_a()
    part_b(a.qmax)
    print('\nno violations.')


if __name__ == '__main__':
    main()
