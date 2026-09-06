#!/usr/bin/env python3.13
"""s454_composite_q.py -- PRE-REGISTERED.  Queue item (delta): extend [NCYL-298]'s
criterion (*) from PRIME `Q` to COMPOSITE `Q`.  It does more than that -- the term-count
argument is replaced outright by an exact factorisation, so the prime hypothesis, the
`Q >= 11` threshold and the `Q = 5, 7` direct-verification exception all disappear.

PRIOR ART: `rulings.py --grep 'composite q' 'phi_q' 'symmetric difference'` -> NO rulings
match; `grep -n 'Phi_15\\|Phi_Q\\|Φ_Q' *.md` -> only [NCYL-298] itself, its `TOOLS.md` row
and the s453 handoff, i.e. the statement of the gap and nothing attacking it.  The object
(the criterion (*), the exponent map `alpha_i = -2 P i mod Q`, the four-index symmetric
difference, the admissible pair range `0 <= a < b <= m-1`) is entirely [NCYL-298] /
`probes/s453_pole_module.py::support_mod_Q`, and is USED verbatim, not redefined.
[NCYL-299]'s rank law is NOT used -- the whole point is that this route needs no lattice.
=> NEW here: the factorisation H2 and the exact zero-set characterisation H3.
⚠ [NCYL-296] bars routing this residual through CF depth / renormalization.  It does not
  bite: nothing below reads a CF digit.  ⚠ [OPS-228] (the content-division trap) does not
  bite either: nothing below touches a `Coord` vector -- the exponents are plain ints.

==========================================================================================
THE THEOREM (proved here; it REPLACES [NCYL-298]'s prime-`Q` term count).

Setting, all of it [NCYL-298] and unchanged.  `Q >= 5` odd, `gcd(P, Q) = 1`,
`m = (Q-1)/2`, `beta = -2P`.  Reducing (*) mod `2` and then mod `Phi_Q(x)` sends
`eta_c |-> x^c + x^{-c}` with node exponents `alpha_i = beta i mod Q`, so the criterion
for the pole pair `(a, b)` -- edges `a` and `b` of the quotient path, i.e. node indices
`{a, a+1}` and `{b, b+1}` -- is the nonvanishing of

        f_{a,b}  :=  sum_{i in {a, a+1, b, b+1}}  ( x^{beta i} + x^{-beta i} )

in `F_2[x] / Phi_Q(x)`.  ⚠ The reduction is a RING SURJECTION
`Z[zeta_{4Q}]/2 = F_2[x]/(Phi_Q^2)  ->>  F_2[x]/Phi_Q` (mod `2`, `Phi_{4Q} == Phi_Q^2`),
so nonzero downstairs implies nonzero upstairs -- the direction the criterion needs.

⇒⇒ THE FACTORISATION, which is the whole content.  Work in `F_2[x]/(x^Q - 1)` and put
`y = x^beta`.  Grouping the two consecutive indices of each edge,

        (y^a + y^{-a}) + (y^{a+1} + y^{-(a+1)})  =  (y + 1) * ( y^a + y^{-(a+1)} ),

and summing the two edges, then pulling `y^{-(a+b)}` out of the inverse pair,

    ⇒⇒   f_{a,b}  =  (y + 1) * (y^a + y^b) * (1 + y^{-(a+b+1)}).            [H2, exact]

⇒⇒ THE ZERO SET.  Let `zeta` be ANY primitive `Q`-th root of unity in `F_2bar` and put
`xi = zeta^beta`, again primitive because `gcd(beta, Q) = gcd(2P, Q) = 1` (`Q` odd).
Since `F_2[x]/Phi_Q` embeds in `F_2bar` factorwise, `f_{a,b} == 0 mod Phi_Q` iff
`f_{a,b}(zeta) = 0` for every such `zeta`, and

        f_{a,b}(zeta)  =  (xi + 1) * (xi^a + xi^b) * (1 + xi^{-(a+b+1)}).

`xi != 1`, so the three factors vanish iff respectively: never; `Q | a - b`; `Q | a+b+1`.

    ⇒⇒   f_{a,b} == 0  (mod Phi_Q)   <=>   Q | (a - b)   or   Q | (a + b + 1).   [H3]

⇒⇒ THE RANGE DISCHARGES BOTH.  The admissible pairs are `0 <= a < b <= m-1`, so
        `1 <= b - a <= m - 1 < Q`        and        `1 <= a + b + 1 <= 2m - 2 = Q - 3`,
neither of which is `0 mod Q`.  ⇒ `f_{a,b} != 0` for EVERY odd `Q >= 5`, every coprime
`P`, every admissible pair, and BOTH CLASSES.                                        ∎

⇒⇒ THE COROLLARY, AND IT IS WHAT SETTLES WHAT s453 ALREADY HAD.  Each of the two
  vanishing conditions makes the SUPPORT EMPTY outright: `Q | a - b` doubles every index,
  and `Q | a + b + 1` gives `beta b == -beta(a+1)` and `beta(b+1) == -beta a`, so edge
  `b`'s two exponent pairs are edge `a`'s, swapped.  Hence

    ⇒⇒   support EMPTY   <=>   `f == 0` in `F_2[x]/(x^Q - 1)`   <=>   `f == 0 mod Phi_Q`,

  all three equivalent, PROVED above and measured `785076/785076` (H3).  ⇒ PASSING TO THE
  CYCLOTOMIC QUOTIENT LOSES NO INFORMATION ON THIS FAMILY, so [NCYL-298]'s H8 arm
  (`support NONEMPTY`, `1650636/1650636`) was never the weaker statement it was booked as
  -- at composite `Q` just as much as at prime `Q`, it WAS the criterion.

⇒ NO PRIMALITY, NO TERM COUNT, NO SUPPORT SIZE.  [NCYL-298]'s argument needed
  `|support| <= 8 < Q` together with "at prime `Q` the only nonzero multiple of `Phi_Q`
  of degree `< Q` is `Phi_Q` itself".  The factorisation never counts terms.
⇒ AND IT SUBSUMES THE TWO EXCEPTIONS: `Q = 5, 7` were "below the term count, covered by
  direct verification" -- they are now inside the theorem.

⚠⚠ TWO CORRECTIONS TO [NCYL-298]'s SCOPE CAVEAT, both found here and both ALGEBRAIC.
  (i) `Phi_15 mod 2` has degree `8` and **`7`** terms, not `8` (`x^8+x^7+x^5+x^4+x^3+x+1`;
      checked against sympy, H6).  The caveat's "`Phi_15` has exactly `8` terms, which is
      the largest support this construction produces, so `Q = 15` turns on whether one
      `8`-set equals another" is therefore wrong on the figure AND on the reading.
  (ii) More substantively, the caveat OVERSTATED the gap.  The gap was real *in the
      ARGUMENT* -- the term count genuinely does not cover composite `Q`, and that much
      stands -- but by the corollary the CONCLUSION held anyway and H8's composite-`Q`
      rows were already measuring the criterion itself.  ⇒ What s454 supplies is the
      PROOF, not the coverage; the coverage was there.  ⚠ Do not re-quote the `8`-terms
      sentence, and do not describe s453's composite-`Q` rows as "weaker evidence".
⇒ ROBUST TO AN OFF-BY-ONE IN THE RANGE: even at `a, b in [0, m]` the bound reads
  `a + b + 1 <= 2m = Q - 1 < Q`, so the conclusion survives the wider index range too.
  It is SHARP one step further out: `a + b + 1 = Q` is realisable and does vanish (H3).

⚠⚠ WHAT THIS DOES AND DOES NOT MOVE.  It closes the composite-`Q` scope caveat on
  [NCYL-298]'s criterion and NOTHING ELSE.  Every other qualification stands verbatim:
  the lemma is SILENT ON `j' = j` (`s_j - s_j = 0 in 2 Lambda` always), so `p = 0` is
  STILL NOT PROVED ON EITHER CLASS; it is CLASS-BLIND; and it lives WITHIN [NCYL-294]'s
  model, whose global-coordinate closure is VERIFIED (`420/420`), not proved.  The links
  ABOVE the arithmetic -- H7's closed form `s_j == L_j + L_{j+1} (mod 2)`, derivation-
  checked `1288/1288` -- are untouched in status by anything here.

HYPOTHESES.  ⚠ [OPS-041]: an arm that cannot fail is not evidence.
  H1  DERIVATION CHECK, NOT EVIDENCE.  On real `Chain` rows: `ch.c[m-i] % Q == (-2Pi) % Q`
      and `m == (Q-1)//2`.  This is the bridge from the geometry to the arithmetic and it
      is a re-run of what `support_mod_Q` already computes -- it can only catch an
      indexing slip, which is exactly why it is here and why it is not evidence.
  H2  ⇒⇒ THE FACTORISATION, AS AN EXACT POLYNOMIAL IDENTITY in `F_2[x]/(x^Q - 1)`, over
      the FULL index range `0 <= a <= b <= Q-1` (not the admissible one) and every
      coprime `beta`.  It can fail on any single triple.
  H3  ⇒⇒ THE ZERO SET, AND THE ARM THAT COULD KILL THE THEOREM.  `f == 0 mod Phi_Q` is
      scored against the predicate `Q | (a-b) or Q | (a+b+1)` over the same full range.
      ⚠ NON-VACUOUS BY CONSTRUCTION: the predicted-zero set is LARGE and nonempty, so
      this is not a `0 == 0` tautology -- the criterion genuinely fails off the admissible
      range, and the theorem is the observation that the range excludes it.  A run
      reporting `zeros = 0` must be booked as a VACUOUS audit, not a pass.
  H4  THE RANGE ARITHMETIC: max `a+b+1` over admissible pairs, against `Q - 3`.
  H5  ⇒⇒ THE COMPOSITE ROWS, MEASURED DIRECTLY MOD `Phi_Q`.  s453's H8 checked only that
      the SUPPORT IS NONEMPTY, i.e. `f != 0` in `F_2[x]/(x^Q - 1)` -- strictly WEAKER than
      `f != 0 mod Phi_Q`, and the gap between them IS the composite-`Q` gap.  H5 runs the
      honest long division on every admissible pair at every odd composite `Q` in range.
  H6  INSTRUMENT AUDIT: the bitmask `GF(2)` reduction against sympy's `Poly(..., modulus=2)`
      remainder, different code and different library, on a sample.

RESULTS (`all --qmax=99 --full=61`, `35.4 s`, no tracing and no `Chain` outside H1).
  H2  `785076/785076` exact identities, `0` bad -- full index range, every coprime `beta`,
      odd `Q = 5..61`.
  H3  ⇒⇒ zero set matches the predicate `785076/785076`, `0` disagreements, with
      **`48864` observed zeros**: the arm is NON-VACUOUS, the criterion really does fail
      off the admissible range, and the range is what excludes it.  Third equivalence
      (support empty) also `785076/785076`.  ⇒ `0` pairs with support nonempty yet zero
      mod `Phi_Q`, at prime or composite `Q` -- the corollary above, measured.
  H4  `55/55` centres: max `a+b+1` over admissible pairs is exactly `Q - 3`.
  H5  ⇒⇒ `6274464/6274464` admissible pairs nonzero mod `Phi_Q`, honest long division,
      odd `Q = 5..99` plus `105, 165, 195, 225, 231, 255, 315` -- **`32` composite `Q` of
      `55`**, and `3.8x` s453's pair count at the strictly-stated mod-`Phi_Q` level.
  H1  `840/840` real rows (`420` paired, `420` lone, odd `Q <= 45`): `alpha_i == -2Pi mod Q`
      and `m == (Q-1)//2`.  DERIVATION CHECK -- it can only catch an indexing slip.
  H6  `480/480` against sympy `GF(2)`, including the `Phi_Q mod 2` masks themselves.

  python3.13 probes/s454_composite_q.py [all|arith|bridge] [--qmax=N] [--full=N]
"""
from __future__ import annotations

import json
import math
import sys
import time

OUT = "data/s454_composite_q.json"
STRESS = (105, 165, 195, 225, 231, 255, 315)


# ---------------------------------------------------------------- GF(2) polynomials
# A polynomial over F_2 is an int bitmask: bit `i` is the coefficient of `x^i`.

def _phi_int(Q):
    n, r = Q, Q
    p = 2
    while p * p <= n:
        if n % p == 0:
            while n % p == 0:
                n //= p
            r -= r // p
        p += 1
    if n > 1:
        r -= r // n
    return r


def phi_mask(Q):
    """`Phi_Q(x) mod 2` as a bitmask, by exact integer division of `x^Q - 1` by the
    lower cyclotomics -- stdlib only, no sympy in the hot path."""
    num = (1 << Q) | 1                      # x^Q - 1 == x^Q + 1 over F_2
    for d in range(1, Q):
        if Q % d == 0:
            num = _divexact(num, phi_mask(d)) if d > 1 else _divexact(num, 0b11)
    return num


def _divexact(f, g):
    """Exact quotient `f / g` over `F_2` (caller guarantees divisibility)."""
    dg = g.bit_length() - 1
    q = 0
    while f and f.bit_length() - 1 >= dg:
        sh = f.bit_length() - 1 - dg
        q ^= 1 << sh
        f ^= g << sh
    return q


def rem2(f, g):
    dg = g.bit_length() - 1
    while f and f.bit_length() - 1 >= dg:
        f ^= g << (f.bit_length() - 1 - dg)
    return f


def mulQ(u, v, Q):
    """XOR-convolution in `F_2[x]/(x^Q - 1)`."""
    out = 0
    for i in range(u.bit_length()):
        if (u >> i) & 1:
            for j in range(v.bit_length()):
                if (v >> j) & 1:
                    out ^= 1 << ((i + j) % Q)
    return out


def fpoly(Q, beta, a, b):
    """`f_{a,b} = sum_{i in {a,a+1,b,b+1}} (x^{beta i} + x^{-beta i})` in `F_2[x]/(x^Q-1)`.

    Verbatim the symmetric difference `s453_pole_module.support_mod_Q` forms, with the
    exponent map `alpha_i = beta i mod Q` (H1 checks that against real `Chain` rows)."""
    sup = 0
    for i in (a, a + 1, b, b + 1):
        e = (beta * i) % Q
        sup ^= 1 << e
        sup ^= 1 << ((-e) % Q)
    return sup


def fsupport(Q, beta, a, b):
    """The same object as an exponent SET -- `|support| <= 8`, so reducing mod `Phi_Q`
    costs `<= 8` table lookups rather than a scan of all `Q` bit positions."""
    sup = set()
    for i in (a, a + 1, b, b + 1):
        e = (beta * i) % Q
        sup ^= {e}
        sup ^= {(-e) % Q}
    return sup


def factored(Q, beta, a, b):
    """`(y+1)(y^a + y^b)(1 + y^{-(a+b+1)})`, `y = x^beta`, in `F_2[x]/(x^Q-1)`."""
    y = lambda k: 1 << ((beta * k) % Q)
    return mulQ(mulQ(y(0) ^ y(1), y(a) ^ y(b), Q), 1 ^ y(-(a + b + 1)), Q)


# ---------------------------------------------------------------- hypotheses

def arith(qmax, fullmax):
    """H2/H3/H4/H5 -- pure arithmetic, no dynamics, no tracing, no `Chain`."""
    res = {"h2": [0, 0], "h2_bad": [], "h3": [0, 0], "h3_bad": [], "h3_zeros": 0,
           "h3_pairs": 0, "h3_gap": 0, "h3_gap_composite": 0, "h3_gap_ex": [],
           "h3_supp": [0, 0], "h3_supp_ex": [],
           "h4": [0, 0], "h4_max": {}, "h5": {}, "h5_bad": [],
           "h5_pairs": 0, "h5_ok": 0, "by_Q": {}}
    qs = [q for q in range(5, qmax + 1, 2)] + [q for q in STRESS if q > qmax]
    for Q in qs:
        m = (Q - 1) // 2
        ph = phi_mask(Q)
        xr = [rem2(1 << e, ph) for e in range(Q)]      # x^e mod Phi_Q, precomputed
        betas = [b for b in range(1, Q) if math.gcd(b, Q) == 1]
        full = Q <= fullmax
        z_here = p_here = 0
        for beta in betas:
            # -- H2 / H3 over the FULL index range (the falsifiable, non-vacuous arm)
            if full:
                for a in range(Q):
                    for b in range(a, Q):
                        f = fpoly(Q, beta, a, b)
                        res["h2"][0 if f == factored(Q, beta, a, b) else 1] += 1
                        if f != factored(Q, beta, a, b) and len(res["h2_bad"]) < 8:
                            res["h2_bad"].append([Q, beta, a, b])
                        red = 0
                        for e in fsupport(Q, beta, a, b):
                            red ^= xr[e]
                        z = (red == 0)
                        pred = ((a - b) % Q == 0) or ((a + b + 1) % Q == 0)
                        res["h3"][0 if z == pred else 1] += 1
                        res["h3_pairs"] += 1
                        z_here += z
                        p_here += pred
                        # ⇒⇒ THE GAP s453's H8 COULD NOT SEE, EXHIBITED: support
                        # NONEMPTY (so H8 passes) yet ZERO mod Phi_Q (so the criterion
                        # FAILS).  Every one of these is a pair the term-count argument
                        # cannot rule out at composite `Q`.
                        # ⇒⇒ THE THIRD EQUIVALENCE: is SUPPORT EMPTY the same
                        # condition?  If yes, s453's H8 (`support nonempty`) was never
                        # weaker than the criterion, at ANY `Q` -- which retroactively
                        # settles what its composite-`Q` rows were worth.
                        res["h3_supp"][0 if (f == 0) == pred else 1] += 1
                        if (f == 0) != pred and len(res["h3_supp_ex"]) < 8:
                            res["h3_supp_ex"].append([Q, beta, a, b])
                        if z and f:
                            res["h3_gap"] += 1
                            if not _is_prime(Q):
                                res["h3_gap_composite"] += 1
                            if len(res["h3_gap_ex"]) < 8:
                                res["h3_gap_ex"].append(
                                    [Q, beta, a, b, len(fsupport(Q, beta, a, b))])
                        if z != pred and len(res["h3_bad"]) < 8:
                            res["h3_bad"].append([Q, beta, a, b, int(z), int(pred)])
            # -- H5: the ADMISSIBLE pairs, direct mod-Phi_Q, at every Q incl. composite
            for a in range(m):
                for b in range(a + 1, m):
                    red = 0
                    for e in fsupport(Q, beta, a, b):
                        red ^= xr[e]
                    res["h5_pairs"] += 1
                    if red:
                        res["h5_ok"] += 1
                    elif len(res["h5_bad"]) < 8:
                        res["h5_bad"].append([Q, beta, a, b])
        res["h3_zeros"] += z_here
        # -- H4: the range arithmetic that discharges the zero set
        mx = (m - 1) + (m - 2) + 1 if m >= 2 else None
        if mx is not None:
            res["h4"][0 if mx == Q - 3 else 1] += 1
            res["h4_max"][str(Q)] = [mx, Q - 3]
        res["by_Q"][str(Q)] = {"m": m, "prime": _is_prime(Q), "phi": _phi_int(Q),
                               "deg_phiQ": ph.bit_length() - 1,
                               "terms_phiQ": bin(ph).count("1"),
                               "betas": len(betas), "full_range": full,
                               "zeros_full": z_here if full else None}
    return res


def _is_prime(n):
    if n < 2:
        return False
    d = 2
    while d * d <= n:
        if n % d == 0:
            return False
        d += 1
    return True


def bridge(qmax):
    """H1 -- the geometry->arithmetic bridge, a DERIVATION CHECK and not evidence."""
    import s450_chain_maps as cm            # noqa: F401  (imported by s451_quotient_path)
    import s451_quotient_path as qp
    res = {"rows": 0, "ok": 0, "bad": 0, "bad_ex": [], "m_ok": 0, "m_bad": 0,
           "classes": {"paired": 0, "lone": 0}}
    for Q in range(5, qmax + 1, 2):
        for P in range(1, Q):
            if math.gcd(P, Q) != 1:
                continue
            for eps in (0, 1):
                ch, lo, hi, _closes = qp.build(P, Q, eps)
                m = ch.m_idx
                res["rows"] += 1
                res["classes"]["lone" if ch.lone else "paired"] += 1
                res["m_ok" if m == (Q - 1) // 2 else "m_bad"] += 1
                good = all(ch.c[m - i] % Q == (-2 * P * i) % Q for i in range(m + 1))
                res["ok" if good else "bad"] += 1
                if not good and len(res["bad_ex"]) < 5:
                    res["bad_ex"].append([P, Q, eps])
    return res


def audit(qmax):
    """H6 -- the bitmask reduction against sympy `GF(2)`, different code and library."""
    from sympy import Poly, cyclotomic_poly, symbols
    x = symbols("x")
    ok = bad = 0
    ex = []
    for Q in range(5, qmax + 1, 2):
        ph = phi_mask(Q)
        # the mask itself, against sympy's cyclotomic
        sy = Poly(cyclotomic_poly(Q, x), x, modulus=2)
        sy_mask = 0
        for i, v in enumerate(sy.all_coeffs()[::-1]):
            if int(v) % 2:
                sy_mask |= 1 << i
        if sy_mask != ph:
            bad += 1
            ex.append([Q, "PHI-MASK"])
            continue
        m = (Q - 1) // 2
        for beta in (1, 2, Q - 1):
            if math.gcd(beta, Q) != 1:
                continue
            for a in range(min(m, 4)):
                for b in range(a + 1, min(m, 6)):
                    f = fpoly(Q, beta, a, b)
                    mine = rem2(f, ph) == 0
                    fp = Poly([1 if (f >> e) & 1 else 0
                               for e in range(f.bit_length() - 1, -1, -1)] or [0],
                              x, modulus=2)
                    theirs = (fp % sy).is_zero
                    if mine == bool(theirs):
                        ok += 1
                    else:
                        bad += 1
                        if len(ex) < 8:
                            ex.append([Q, beta, a, b, mine, bool(theirs)])
    return {"ok": ok, "bad": bad, "ex": ex}


# ---------------------------------------------------------------- main

def main():
    mode = "all"
    qmax = 99
    fullmax = 61
    for a in sys.argv[1:]:
        if a.startswith("--qmax="):
            qmax = int(a.split("=")[1])
        elif a.startswith("--full="):
            fullmax = int(a.split("=")[1])
        elif not a.startswith("-"):
            mode = a
    t0 = time.time()
    out = {"mode": mode, "qmax": qmax, "fullmax": fullmax}
    if mode in ("all", "arith"):
        out["arith"] = arith(qmax, fullmax)
        r = out["arith"]
        print(f"H2 factorisation (full range, all beta): ok={r['h2'][0]} "
              f"BAD={r['h2'][1]}  {r['h2_bad'][:3]}")
        print(f"H3 zero set == (Q|a-b or Q|a+b+1): agree={r['h3'][0]} "
              f"DISAGREE={r['h3'][1]}  over {r['h3_pairs']} pairs; "
              f"observed zeros={r['h3_zeros']}  (NON-VACUITY: zeros must be > 0)")
        print(f"   of those, SUPPORT NONEMPTY yet zero mod Phi_Q = {r['h3_gap']} "
              f"({r['h3_gap_composite']} at composite Q) -- exactly the pairs s453's H8 "
              f"passes but the criterion fails; ex {r['h3_gap_ex'][:3]}")
        print(f"   3rd equivalence (support EMPTY == same predicate): "
              f"agree={r['h3_supp'][0]} DISAGREE={r['h3_supp'][1]} {r['h3_supp_ex'][:3]}"
              f"  ⇒ support-nonempty IS the criterion on this family")
        print(f"H4 max(a+b+1) over admissible pairs == Q-3: ok={r['h4'][0]} "
              f"BAD={r['h4'][1]}")
        nprime = sum(1 for v in r["by_Q"].values() if not v["prime"])
        print(f"H5 admissible pairs nonzero mod Phi_Q: {r['h5_ok']}/{r['h5_pairs']}  "
              f"BAD={len(r['h5_bad'])} {r['h5_bad'][:3]}   "
              f"({nprime} composite Q of {len(r['by_Q'])})")
        p15 = r["by_Q"].get("15")
        if p15:
            print(f"   Q=15 check: deg Phi_15 mod 2 = {p15['deg_phiQ']}, "
                  f"terms = {p15['terms_phiQ']}")
    if mode in ("all", "bridge"):
        out["bridge"] = bridge(min(qmax, 45))
        b = out["bridge"]
        print(f"H1 bridge alpha_i == -2Pi mod Q (DERIVATION, not evidence): "
              f"ok={b['ok']} BAD={b['bad']} {b['bad_ex'][:3]}; "
              f"m==(Q-1)//2 ok={b['m_ok']} BAD={b['m_bad']}; "
              f"rows={b['rows']} {b['classes']}")
    if mode == "all":
        out["audit"] = audit(min(qmax, 31))
        a = out["audit"]
        print(f"H6 bitmask vs sympy GF(2): ok={a['ok']} BAD={a['bad']} {a['ex'][:3]}")
    out["secs"] = round(time.time() - t0, 1)
    with open(OUT, "w") as fh:
        json.dump(out, fh)
    print(f"-> {OUT}  ({out['secs']} s)")


if __name__ == "__main__":
    main()
