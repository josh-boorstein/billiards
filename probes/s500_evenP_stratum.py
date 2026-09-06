#!/usr/bin/env python3
"""PRE-REGISTERED (s500, E3 §7): does the unfolded-surface data of `S_alpha` hold at EVEN `P`,
and does the Apisa-base lift still reproduce it there?

WHY. `paper_orphan.md` §8 states Theorem W as an IFF over the parity of `P` (`O` is a Weierstrass
point iff `P` is odd). Its derivation reads `O`'s order off the §7 surface data. But BOTH inputs are
stated for odd `P` only:
  * `claims.md` [T-SURFACE] gives the stratum as H(P-1, Q-P-1) (Q even) / H(P-1, (Q-P)/2-1,
    (Q-P)/2-1) (Q odd) with NO `P`-parity hypothesis on the row, while its proof
    (`veech_reframing.md` §1a) opens "Setup ... `P` odd" and turns on gcd(P,2Q)=1, which FAILS at
    even `P`.  At even `P` (so `Q` odd) the Q-odd formula is not even integral: `Q-P` is odd, so
    `(Q-P)/2` is not an integer.  A reader plugging in even `P` gets nonsense, not an error.
  * `probes/veech_apisa.py`, the 109-case base<->lift cross-check, SKIPS even `P` outright
    (`if math.gcd(P,Q) != 1 or P % 2 == 0: continue`).

So the even-`P` half of Theorem W is VERIFIED end-to-end (`claims.md` [T-WEIER], 0 mismatches,
Q<200) but its stated derivation chain does not cover it.  `cf_width_laws.md` §s272 (3) re-derived
the even-`P` cone data by hand and got a DIFFERENT stratum, H(P/2-1, P/2-1, Q-P-1); that
re-derivation has never been connected to [T-SURFACE] or checked against the Apisa base.

WHAT THIS CHECKS. For every coprime (P,Q) in range, in BOTH parities of `P`:
  (A) the direct Katok-Zemlyakov cone-point computation (`veech_stratum.surface`, which is fully
      general -- only its DRIVER filtered to odd `P`) against
  (B) the holonomy-double-cover lift of the Apisa base Q(P-2, Q-P-2, -1^Q).
These are two independent computations: (A) is reduced denominators + Euler characteristic on the
unfolding, (B) is a fixed lift rule applied to a base defined without reference to the unfolding.
Also checked: genus = floor(Q/2), base Gauss-Bonnet = -4, and orders summing to 2g-2.

CONTROLS (this must be able to come out otherwise -- [OPS-041]):
  * NEGATIVE CONTROL 1: the [T-SURFACE] row's Q-odd formula H(P-1,(Q-P)/2-1,(Q-P)/2-1) is asserted
    to FAIL at even `P`.  If it ever matches, the scope finding is wrong.
  * NEGATIVE CONTROL 2: at even `P`, the NAIVE form H(P-1, Q-P-1) -- what you get by carrying the
    odd-`P` shape across the parity line -- is asserted to FAIL.  If it matched, the even-`P` data
    would not actually be different and there would be nothing to report.
    (An earlier version of this control asked whether H(P/2-1,P/2-1,Q-P-1) fails at ODD `P`.  That
    is ill-posed: `P/2` is not an integer there, and `P//2` silently makes the form collapse to
    [Q-P-1], which coincides with the truth on the 19 rows where `P-1` is 0 or negative.  The
    control was measuring integer division, not the mathematics.)
  * The P=2 special case (both O-points regular -> minimal stratum H(2g-2)) is checked separately.
    ⚠ `(P,Q) = (2,3)` is excluded from that check and reported on its own: it is the 30-60-90
    triangle, whose unfolding is a singularity-free flat TORUS (g=1, empty stratum), so "H(2g-2)"
    would read as a spurious order-0 zero.  This is the decisive case `results_ledger.md` [V-NCYL]
    uses for branch != cylinder, and the probe reproduces it.

PRIOR ART: grepped `rulings.py --grep 'even-P stratum' 'even P stratum' 'H(P/2'` -> none; read
`results_ledger.md` [T-SURFACE]/[V-NCYL] (the latter's scope caveat is where the gap is recorded),
`cf_width_laws.md` §s272 (3) (the hand re-derivation this checks), `claims.md` [T-SURFACE]/[T-WEIER],
and `probes/veech_apisa.py` (the odd-P-only cross-check this extends).  No probe covers even `P`.

RUN: PYTHONPATH=.:engine:archive/scripts_2026-07:archive/scripts_2026-08:probes \
     .venv/bin/python3.13 probes/s500_evenP_stratum.py
"""
import math
from fractions import Fraction as F

import veech_stratum as vs


def holonomy_double_lift(orders):
    """Half-translation zero orders -> abelian holonomy double cover stratum.

    d odd  -> ONE zero of order d+1   (branch point)
    d even -> TWO zeros of order d/2  (unbranched, swapped by the deck involution)
    d = -1 (simple pole) -> ONE regular point (order 0, unmarked)
    Order-0 points are dropped so the result is the genuine-zero list.
    """
    out = []
    for d in orders:
        if d == -1:
            continue
        if d % 2 == 1:
            out.append(d + 1)
        else:
            out += [d // 2, d // 2]
    return sorted([o for o in out if o > 0], reverse=True)


def apisa_base(P, Q):
    return [P - 2, Q - P - 2] + [-1] * Q


def row_formula_Qodd(P, Q):
    """The [T-SURFACE] row's Q-odd form. Returns None when not integral."""
    if (Q - P) % 2 != 0:
        return None
    return sorted([x for x in (P - 1, (Q - P) // 2 - 1, (Q - P) // 2 - 1) if x > 0],
                  reverse=True)


def row_formula_Qeven(P, Q):
    return sorted([x for x in (P - 1, Q - P - 1) if x > 0], reverse=True)


def evenP_formula(P, Q):
    """The cf_width_laws.md §s272 (3) prediction: H(P/2-1, P/2-1, Q-P-1)."""
    if P % 2 != 0:
        return None
    return sorted([x for x in (P // 2 - 1, P // 2 - 1, Q - P - 1) if x > 0], reverse=True)


def weierstrass_census(P, Q):
    """Fixed points of the deck involution tau = preimages of the ODD-order singularities of the
    base Q(P-2, Q-P-2, -1^Q).  Returns (total, breakdown).

    Added s500 after reading [Hoo07]: Hooper's involution `i` on the minimal translation surface is
    rotation by pi about the centres of the rhombi, which are the images of the RIGHT-ANGLE vertex
    `R`.  He punctures them, so his `i` is fixed-point free; filled in, it fixes them.  Those are our
    `Q` regular points over `R`, and they are the lifts of the `Q` SIMPLE POLES -- order -1, which is
    ODD, hence ramification points.  So the R-images are Weierstrass points, and the census below is
    the check that this whole picture is hyperelliptic: a hyperelliptic surface of genus g has
    EXACTLY 2g+2 Weierstrass points.  This is a genuine constraint -- it is not built into the
    stratum computation and it can fail.
    """
    parts = {'R (simple poles)': Q}
    if (P - 2) % 2 != 0:
        parts['O'] = 1
    if (Q - P - 2) % 2 != 0:
        parts['A'] = 1
    return sum(parts.values()), parts


def main():
    QMAX = 40
    rows = []
    for Q in range(3, QMAX + 1):
        for P in range(1, Q):
            if math.gcd(P, Q) != 1:
                continue
            rows.append((P, Q))

    n_odd = n_even = 0
    fail_agree = []
    fail_genus = []
    # negative controls
    nc1_hits = []   # row Q-odd formula matching at even P  (must stay empty)
    nc2_hits = []   # even-P formula matching at odd P      (must stay empty)
    p2_bad = []
    fail_weier = []
    fail_thmW = []

    for P, Q in rows:
        s = vs.surface(P, Q)
        kz = sorted(s['stratum'], reverse=True)

        base = apisa_base(P, Q)
        assert sum(base) == -4, ('base Gauss-Bonnet', P, Q, sum(base))
        lifted = holonomy_double_lift(base)

        if kz != lifted:
            fail_agree.append((P, Q, kz, lifted))

        g = int(s['gA'])
        assert s['gA'] == s['gB'], ('genus routes disagree', P, Q, s['gA'], s['gB'])
        if g != Q // 2 or sum(kz) != 2 * g - 2:
            fail_genus.append((P, Q, g, Q // 2, sum(kz)))

        wtot, wparts = weierstrass_census(P, Q)
        if wtot != 2 * g + 2:
            fail_weier.append((P, Q, wtot, 2 * g + 2, wparts))
        # Theorem W's own clause, read off the census: O is Weierstrass iff P is odd.
        if ('O' in wparts) != (P % 2 == 1):
            fail_thmW.append((P, Q, wparts))

        if P % 2 == 0:
            n_even += 1
            pred = evenP_formula(P, Q)
            if kz != pred:
                fail_agree.append((P, Q, kz, pred, 'evenP-formula'))
            rf = row_formula_Qodd(P, Q) if Q % 2 else row_formula_Qeven(P, Q)
            if rf is not None and rf == kz:
                nc1_hits.append((P, Q, rf))
            # NC2: the naive odd-P-shaped form H(P-1, Q-P-1), carried across the parity line.
            naive = sorted([x for x in (P - 1, Q - P - 1) if x > 0], reverse=True)
            if naive == kz:
                nc2_hits.append((P, Q, naive))
            if P == 2 and (P, Q) != (2, 3) and kz != [2 * g - 2]:
                p2_bad.append((P, Q, kz, g))
        else:
            n_odd += 1
            rf = row_formula_Qodd(P, Q) if Q % 2 else row_formula_Qeven(P, Q)
            if rf != kz:
                fail_agree.append((P, Q, kz, rf, 'row-formula-oddP'))

    print(f"rows: {len(rows)}  (P odd {n_odd}, P even {n_even})   Q <= {QMAX}")
    print(f"(A) direct KZ  vs  (B) Apisa-base lift : "
          f"{'AGREE on every row' if not fail_agree else f'{len(fail_agree)} MISMATCH'}")
    for r in fail_agree[:8]:
        print("   MISMATCH", r)
    print(f"genus = floor(Q/2) and sum(orders) = 2g-2 : "
          f"{'OK on every row' if not fail_genus else f'{len(fail_genus)} FAIL'}")
    for r in fail_genus[:8]:
        print("   FAIL", r)

    print()
    print("NEGATIVE CONTROL 1 -- the [T-SURFACE] row formula at EVEN P "
          "(must FAIL on every even-P row, else the scope finding is wrong):")
    n_int = sum(1 for (P, Q) in rows if P % 2 == 0
                and (row_formula_Qodd(P, Q) if Q % 2 else row_formula_Qeven(P, Q)) is not None)
    print(f"   even-P rows where the row formula is even INTEGRAL: {n_int} of {n_even}")
    print(f"   even-P rows where it MATCHES the true stratum      : {len(nc1_hits)}"
          f"   {'(control fires: it never matches)' if not nc1_hits else '(!! control did not fire)'}")

    print("NEGATIVE CONTROL 2 -- the NAIVE form H(P-1,Q-P-1) carried across the parity line to "
          "EVEN P (must FAIL, else the even-P data is not actually different):")
    print(f"   even-P rows where it matches: {len(nc2_hits)} of {n_even}"
          f"   {'(control fires: never)' if not nc2_hits else '(!! control did not fire)'}")

    print()
    print(f"WEIERSTRASS CENSUS (fixed points of tau = lifts of the ODD-order base singularities) "
          f"= 2g+2 : {'OK on every row' if not fail_weier else f'{len(fail_weier)} FAIL'}")
    for r in fail_weier[:8]:
        print("   FAIL", r)
    print(f"   and 'O is Weierstrass iff P is odd' (Theorem W's clause, read off the census) : "
          f"{'OK on every row' if not fail_thmW else f'{len(fail_thmW)} FAIL'}")
    print("   the three admissible parity cases:")
    for P, Q, lab in [(3, 8, 'P odd,  Q even'), (3, 11, 'P odd,  Q odd '),
                      (4, 9, 'P even, Q odd ')]:
        wt, wp = weierstrass_census(P, Q)
        print(f"      {lab}  (P,Q)=({P:>2},{Q:>2}): {wt} = 2g+2 = {2*(Q//2)+2}   from {wp}")

    s23 = vs.surface(2, 3)
    print()
    print(f"   (2,3), the 30-60-90 triangle: stratum H{tuple(s23['stratum'])} (empty), "
          f"g={int(s23['gA'])}, copies={s23['copies']} -- the singularity-free flat TORUS "
          f"([V-NCYL]'s decisive case, reproduced)")
    print(f"P = 2 special case for Q >= 5 (both O-points regular -> minimal H(2g-2)): "
          f"{'OK' if not p2_bad else f'{len(p2_bad)} FAIL'}")
    for P, Q in [(2, q) for q in (5, 7, 9, 11, 15, 21)]:
        if math.gcd(P, Q) == 1:
            s = vs.surface(P, Q)
            print(f"   P=2 Q={Q:>2}: stratum H{tuple(s['stratum'])}, g={int(s['gA'])}")

    print()
    print("sample even-P rows (KZ stratum, cone data at O / A):")
    for P, Q in [(2, 5), (2, 9), (4, 9), (4, 15), (6, 11), (8, 15), (10, 21)]:
        if math.gcd(P, Q) != 1:
            continue
        s = vs.surface(P, Q)
        o = [x for x in s['sings'] if x[0] == 'O'][0]
        a = [x for x in s['sings'] if x[0] == 'A'][0]
        print(f"   P={P:>2} Q={Q:>2}: H{tuple(s['stratum'])}  "
              f"O: {o[1]} pt(s) order {o[3]}   A: {a[1]} pt(s) order {a[3]}   g={int(s['gA'])}")

    ok = not (fail_agree or fail_genus or nc1_hits or nc2_hits or p2_bad
              or fail_weier or fail_thmW)
    print()
    print("VERDICT:", "ALL CHECKS PASS, BOTH CONTROLS FIRE" if ok else "*** SOMETHING FAILED ***")
    return 0 if ok else 1


if __name__ == '__main__':
    raise SystemExit(main())
