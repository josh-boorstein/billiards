#!/usr/bin/env python3.13
"""s494 -- KEEPER: the arms behind `necklace_paper.md` Lemma 5.4 (the interiority invariant).

PRE-REGISTERED.  Written while promoting [NCYL-370] into `claims.md` `[NECK-SELF]` and
into the manuscript as Lemma 5.4 + Theorem 5.5, to check the FOLD clause in the form the
paper states it rather than the form the ruling states it.

PRIOR ART: grepped `rulings.py --grep 'interiority' 'nesting' 'fold sum'` and
`grep -n 'unshared\\|annulus' *.md` -> [NCYL-370] (Lemma I, in the ruling's own
`y in (L+a, L+b)` form), [NCYL-326] (Prop 7.2, nesting is definitional, `5032/5032`),
[NCYL-301] (Thm 5.2 Steps 1-4), [NCYL-308] (the bar s492 removed).  The ANNULUS form
below -- `s_e` = the sum of the two UNSHARED endpoints -- appears in NO ruling; it is
this file's restatement, and it is the one the manuscript prints, so it needs a check
of its own rather than an inherited one.

⇒⇒ WHAT IS CHECKED, AND WHY EACH LINE CAN FAIL.

  A1  Adjacent path interfaces share EXACTLY ONE endpoint.        (Prop 7.2 -- a regression.)
  A2  `s_e` == the sum of the two UNSHARED endpoints of the nested pair.
      ⇒ equivalently: `y |-> s_e - y` is the reflection of the OPEN ANNULUS
        `J_far \\ closure(J_near)` in its own midpoint, hence maps it onto itself.
      This is the fold clause of Lemma 5.4 and it is the ONE line that is new here.
  A3  `s_e/2` is interior to the LARGER interval and to that one only.  (Thm 5.2 Step 2.)
  A4  A CORNER launch coordinate is strictly interior to `J_far`.     (s492's finding.)

  C1  ⇐ THE CONTROL FOR A2, and it is what makes A2 evidence: substitute a SHARED
      endpoint for the near interface's unshared one and re-score.  A2 is an identity
      about which endpoints get summed, so a scorer that cannot distinguish the two
      would read `2240/2240` either way ([OPS-041]).

⚠ A4 IS NOT FULL AND THE SHORTFALL IS THE POINT.  It reads `2106/2240`; the misses are
EXACTLY the DEGENERATE FLANKS (a node's interval is a point, so the "free" endpoint IS
the shared one), which `s456_neckgb.corner_start` already returns `None` on.  So the
invariant's base holds on every corner prong that EXISTS -- `2106/2106` -- and the paper
now states the tie and degenerate-flank exclusions on the face of the definition, because
without them Lemma 5.4's base is false.  A5 below asserts that diagnosis in code
([OPS-043]) rather than leaving it to prose.

RESULTS (`Q <= 25`, `268` class rows, `2240` non-tie edges):
  A1 2240/2240 · A2 2240/2240 · A3 2240/2240 · A4 2106/2240 (misses == degenerate 134/134)
  A5 corner_start exists 2106, its `y` strictly interior to its own node 2106/2106
  C1 134/2240 -- fires only where the two endpoints coincide, i.e. it discriminates.

⚠ DO NOT quote A1/A3 as evidence for anything: A1 is Prop 7.2 (PROVED, [NCYL-326]) and A3
is Thm 5.2 Step 2 (PROVED, [NCYL-301]).  They are regressions on the tracer.  The evidence
in this file is A2 against C1, and A4's shortfall against A5.

Run: PYTHONPATH=.:engine:archive/scripts_2026-07:archive/scripts_2026-08:probes \\
     .venv/bin/python3.13 probes/s494_lemma_interiority.py [--qmax 25]
"""
from __future__ import annotations

import sys
from math import gcd

import s456_neckgb as gb


def run(qmax: int = 25, qmin: int = 5) -> dict:
    r = {k: 0 for k in ("rows", "edges", "A1", "A2", "A3", "A4", "C1",
                        "degen", "miss_degen", "cs_exists", "A5")}
    for Q in range(qmin, qmax + 1, 2):
        for P in range(1, Q):
            if gcd(P, Q) != 1:
                continue
            for eps in (0, 1):
                got = gb._chain(P, Q, eps)
                if got is None:
                    continue
                ch, Lo, Hi, S, ops = got
                r["rows"] += 1
                for j in range(len(S)):
                    a, b = j, j + 1
                    c = ops.cmp(ops.sub(Hi[a], Lo[a]), ops.sub(Hi[b], Lo[b]))
                    if c == 0:
                        continue                      # tie kite: no fold, no corner prong
                    far, near = (a, b) if c > 0 else (b, a)
                    r["edges"] += 1

                    lo_sh, hi_sh = ops.eq(Lo[a], Lo[b]), ops.eq(Hi[a], Hi[b])
                    if lo_sh ^ hi_sh:
                        r["A1"] += 1

                    # the two UNSHARED endpoints of the nested pair
                    u_near, u_far = ((Hi[near], Hi[far]) if lo_sh else (Lo[near], Lo[far]))
                    if ops.eq(S[j], u_near + u_far):
                        r["A2"] += 1

                    # C1: swap the near unshared endpoint for the SHARED one
                    shared = Lo[far] if lo_sh else Hi[far]
                    if ops.eq(S[j], shared + u_far):
                        r["C1"] += 1

                    h = ops.half(S[j])
                    in_far = ops.cmp(Lo[far], h) < 0 and ops.cmp(h, Hi[far]) < 0
                    in_near = ops.cmp(Lo[near], h) < 0 and ops.cmp(h, Hi[near]) < 0
                    if in_far and not in_near:
                        r["A3"] += 1

                    cpt = Hi[near] if lo_sh else Lo[near]
                    inside = ops.cmp(Lo[far], cpt) < 0 and ops.cmp(cpt, Hi[far]) < 0
                    if inside:
                        r["A4"] += 1
                    degen = ops.eq(Lo[near], Hi[near]) or ops.eq(Lo[far], Hi[far])
                    if degen:
                        r["degen"] += 1
                    if not inside and degen:
                        r["miss_degen"] += 1

                    st = gb.corner_start(Lo, Hi, S, j, ops)
                    if st is not None:
                        r["cs_exists"] += 1
                        (node, _d, y), _n, _f = st
                        if ops.cmp(Lo[node], y) < 0 and ops.cmp(y, Hi[node]) < 0:
                            r["A5"] += 1
    return r


def main() -> int:
    qmax = 25
    if "--qmax" in sys.argv:
        qmax = int(sys.argv[sys.argv.index("--qmax") + 1])
    r = run(qmax)
    n = r["edges"]
    print(f"s494 Lemma 5.4 -- Q <= {qmax}: {r['rows']} class rows, {n} non-tie edges")
    print(f"  A1 shares exactly one endpoint      : {r['A1']}/{n}")
    print(f"  A2 s_e == sum of UNSHARED endpoints : {r['A2']}/{n}   (the annulus reflects onto itself)")
    print(f"  A3 s_e/2 interior to the LARGER only: {r['A3']}/{n}")
    print(f"  A4 corner launch inside J_far       : {r['A4']}/{n}"
          f"   [misses {n - r['A4']}, degenerate flanks {r['degen']}, "
          f"misses that are degenerate {r['miss_degen']}]")
    print(f"  A5 corner_start exists {r['cs_exists']}, y strictly interior to its own node"
          f" : {r['A5']}/{r['cs_exists']}")
    print(f"  C1 CONTROL (shared endpoint swapped): {r['C1']}/{n}   <- must fire only on the degenerate flanks")

    ok = (r["A1"] == n and r["A2"] == n and r["A3"] == n
          and r["miss_degen"] == n - r["A4"]
          and r["A5"] == r["cs_exists"]
          and r["C1"] == r["degen"] and r["C1"] < n)
    print("  ⇒ ALL ARMS AS BANKED" if ok else "  ⇒ ⚠ DIVERGENCE FROM THE BANKED RESULT")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
