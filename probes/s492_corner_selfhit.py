#!/usr/bin/env python3
"""s492_corner_selfhit.py -- PRE-REGISTERED.  DOES THE (G0b) WITNESS JOIN *DISTINCT*
CORNERS, OR IS IT A SELF-HIT?

WHY IT MATTERS, AND IT IS THE WHOLE POINT OF THE FILE.  The pole half of `p = 0` is proved
by a SPLIT: Thm 5.1 (`j' != j`, transverse holonomy, arithmetic in `Lambda/2Lambda`) and
Thm 5.2 (`j' = j`, reversibility + palindrome, NO arithmetic at all).  (G0b) -- the corner
half -- has never been split that way; every route in [NCYL-308]/[NCYL-323] treats "corner
prong reaches a corner" as ONE statement.  If the known witness lies entirely in the
DISTINCT-corner half, then the SELF-HIT half at corners has **no counterexample**, and
[OPS-242]'s gate -- which is what kills a row-independent instrument -- does NOT bind on
it.  That would make the corner self-hit half a candidate for a CP-free proof on the Thm
5.2 template, whose ONE missing input is already named ([NCYL-308] (1): Step 2 needs
`s_e/2` interior to exactly one neighbour, and a corner is the shared endpoint of both).

  THE QUESTION, operationally: for each of the `4` census `Z` prongs, the corner prong
  launched at edge `j` arrives at a point `y`.  Which edge's corner is that?  Compare `y`
  against every `corner_start(j')` launch coordinate on the same row -- they all live in
  the one path coordinate, so the comparison is exact.

  PRE-REGISTERED PREDICTION (and it could come out either way, [OPS-041]): DISTINCT.
  [NCYL-324] (6) reports the two `7/15 eps1` prongs as `FFCCZ` / `CCFFZ`, "exact reverses",
  and [NCYL-317] (1) as "the `4 Z` are one connection walked from its two ends" -- one
  connection walked from both ends is `2 -> 5` and `5 -> 2` if the ends are DISTINCT, and
  `2 -> 2`, `5 -> 5` (two connections) if they are not.  Neither source states which.
  ⚠ NAMED WRONG-HYPOTHESIS IS DELIBERATELY ABSENT ([OPS-044]: self-applied it is 0 for 4).

  ⚠ WHAT A "DISTINCT" ANSWER DOES *NOT* BUY, stated before the run: it does not prove the
  corner self-hit half, it does not supply Step 2's arrival node, and it does not touch
  (G0b-form) or [NCYL-323].  It decides ONE thing -- whether the self-hit half is
  counterexample-free in the census -- and hence whether the Thm 5.2 template is even
  admissible there.  A "SELF" answer kills that route outright, which is why it is worth
  5 seconds.

CONTROLS.
  C1 EXHAUSTIVE-MATCH: the arrival `y` must match EXACTLY ONE `j'`.  Zero matches (the
     arrival is not a corner-launch point at all) or several would mean the identification
     is not well posed and the verdict must be withheld.
  C2 REVERSAL: if prong `j` arrives at `j'`, prong `j'` must arrive at `j`.  Fails => the
     "one connection walked from both ends" reading is wrong and so is this file's frame.
  C3 STEP/FOLD REPRODUCTION: `steps` and `folds` per prong must reproduce [NCYL-324] (6)'s
     `5` and `2`.  A bug check on this file, never an addition to it.
  ⚠ C3 is NOT evidence ([OPS-041]) -- it is predicted by the published record.

PRIOR ART: `rulings.py --grep` on 'self-hit' -> [NCYL-301] (the POLE theorem), [NCYL-298];
'distinct corner', 'corner to corner', 'same corner', 'Z-Z', 'zero to zero' -> nothing;
'arrival node' -> [NCYL-308] (1) only, which names it MISSING and does not ask which
corners the witness joins; 'split' + 'G0b' -> [NCYL-317] (the rel/form split, a DIFFERENT
axis -- arithmetic mechanism, not endpoint identity).  ⇒ ID-grep per [OPS-246] on
[NCYL-307]/[NCYL-308]/[NCYL-317]/[NCYL-319]/[NCYL-323]/[NCYL-324]: [NCYL-307] H5 owns the
`4 Z` census and records `(P,Q,eps,j)` only; [NCYL-324] (6) adds the words and directions
and classifies (U)/(D); NEITHER records the ARRIVAL corner.  ⇒ the question is unasked.
⚠ [OPS-242] does not bite: this is a measurement of the counterexample row, not an
instrument that must admit it.

RUN: PYTHONPATH=.:engine:archive/scripts_2026-07:archive/scripts_2026-08:probes \
     .venv/bin/python3.13 probes/s492_corner_selfhit.py
"""

from __future__ import annotations

import json
import math
import sys

import s451_quotient_path as qp
import s453_pole_module as pm
import s454_self_hit as sh
import s456_neckgb as nb

CAP = 4000
QMAX = 33


def _chain(P, Q, eps):
    ch, lo, hi, closes = qp.build(P, Q, eps)
    if not closes:
        return None
    sf = pm.fold_sums(ch, lo, hi)
    Lo, Hi, S = pm.path_view(ch, lo, hi, sf)
    return Lo, Hi, S, pm.CoordOps(ch)


def _launch_points(Lo, Hi, S, ops):
    """`{j: y_j}` -- every corner prong's launch coordinate on this row."""
    out = {}
    for j in range(len(Lo) - 1):
        cs = nb.corner_start(Lo, Hi, S, j, ops)
        if cs is not None:
            out[j] = cs[0][2]
    return out


def _walk(Lo, Hi, S, ops, st0, cap=CAP):
    """Returns (end_kind, final_state, steps, folds)."""
    st, steps, folds = st0, 0, 0
    for _ in range(cap):
        i, d, _y = st
        if not (0 <= i + d < len(Lo)):
            return "ESC", st, steps, folds
        nxt, kind, _e = sh.step(Lo, Hi, S, st, ops)
        steps += 1
        folds += int(kind == "fold")
        if nxt is None:
            return kind, st, steps, folds
        st = nxt
    return "CAP", st, steps, folds


def census(qmax=QMAX):
    """Every `Z` prong in the census, with the ARRIVAL corner identified."""
    res = {"rows": 0, "prongs": 0, "Z": 0, "z_rows": [],
           "c1_unique": 0, "c1_none": 0, "c1_multi": 0}
    for Q in range(5, qmax + 1, 2):
        for P in range(1, Q):
            if math.gcd(P, Q) != 1:
                continue
            for eps in (0, 1):
                c = _chain(P, Q, eps)
                if c is None:
                    continue
                Lo, Hi, S, ops = c
                res["rows"] += 1
                pts = _launch_points(Lo, Hi, S, ops)
                for j, _y0 in sorted(pts.items()):
                    cs = nb.corner_start(Lo, Hi, S, j, ops)
                    res["prongs"] += 1
                    kind, st, steps, folds = _walk(Lo, Hi, S, ops, cs[0])
                    if kind != "Z":
                        continue
                    res["Z"] += 1
                    yz = st[2]
                    hits = [jj for jj, yy in pts.items() if ops.eq(yy, yz)]
                    if len(hits) == 1:
                        res["c1_unique"] += 1
                    elif not hits:
                        res["c1_none"] += 1
                    else:
                        res["c1_multi"] += 1
                    res["z_rows"].append(
                        {"P": P, "Q": Q, "eps": eps, "from_j": j,
                         "to_j": hits[0] if len(hits) == 1 else None,
                         "n_match": len(hits), "steps": steps, "folds": folds,
                         "self": (len(hits) == 1 and hits[0] == j)})
    return res


def step2_hypothesis(qmax=QMAX):
    """ARM 2 -- DOES STEP 2's INTERIORITY HYPOTHESIS ACTUALLY FAIL AT A CORNER PRONG?

    [NCYL-308] (1) states the obstruction as *a corner is the SHARED endpoint, hence an
    endpoint of BOTH*.  But `corner_start` does NOT launch at the shared endpoint -- it
    launches at the FREE (non-shared) endpoint of the strictly SMALLER interface, the fold
    image of the corner.  And [NCYL-324] H5 classifies all `4` census `Z` as kind (D),
    i.e. arriving at exactly that kind of point, with `U = 0`.  So the point Step 2 must
    localise, on every witness the census has, is NOT the point [NCYL-308] prices.

    ⇒ per the s306 rule (when the residual is *does PROVED thing X cover case Y?*, read
    X's PROOF, not the prose citing it), score the hypothesis itself:

      H2a  the corner-prong LAUNCH point `y_j` is interior to exactly ONE of `J_j`,
           `J_{j+1}` -- and it is the LARGER.  (This is Step 2's hypothesis verbatim, with
           `s_e/2` replaced by `y_j`.)
      H2b  CONTROL, and it is the discriminating one: the SHARED endpoint of the same pair
           is interior to NEITHER.  This is [NCYL-308] (1)'s actual claim.  If H2a and H2b
           come out the SAME, the distinction drawn here is not real and the arm is void.
      H2c  CALIBRATION: the POLE point `S[j]/2` -- Step 2's original subject -- interior to
           exactly one, the larger.  Reproduces [NCYL-301] H3's `4392/4392`; a bug check on
           this file, NEVER evidence ([OPS-041]).
    """
    res = {"rows": 0, "edges": 0,
           "h2a_one_larger": 0, "h2a_one_smaller": 0, "h2a_none": 0, "h2a_both": 0,
           "h2b_none": 0, "h2b_one": 0, "h2b_both": 0,
           "h2c_one_larger": 0, "h2c_other": 0, "bad": []}
    for Q in range(5, qmax + 1, 2):
        for P in range(1, Q):
            if math.gcd(P, Q) != 1:
                continue
            for eps in (0, 1):
                c = _chain(P, Q, eps)
                if c is None:
                    continue
                Lo, Hi, S, ops = c
                res["rows"] += 1
                for j in range(len(Lo) - 1):
                    cs = nb.corner_start(Lo, Hi, S, j, ops)
                    if cs is None:
                        continue
                    res["edges"] += 1
                    _st, near, far = cs
                    y = cs[0][2]

                    def interior(pt, i):
                        return (ops.cmp(Lo[i], pt) < 0) and (ops.cmp(pt, Hi[i]) < 0)

                    inn, inf = interior(y, near), interior(y, far)
                    if inf and not inn:
                        res["h2a_one_larger"] += 1
                    elif inn and not inf:
                        res["h2a_one_smaller"] += 1
                    elif not inn and not inf:
                        res["h2a_none"] += 1
                        if len(res["bad"]) < 10:
                            res["bad"].append({"P": P, "Q": Q, "eps": eps, "j": j,
                                               "arm": "h2a_none"})
                    else:
                        res["h2a_both"] += 1

                    shared = Lo[j] if ops.eq(Lo[j], Lo[j + 1]) else Hi[j]
                    sn, sf_ = interior(shared, near), interior(shared, far)
                    res["h2b_" + ("none" if not (sn or sf_)
                                  else ("both" if (sn and sf_) else "one"))] += 1

                    half = ops.half(S[j]) if hasattr(ops, "half") else None
                    if half is not None:
                        pn, pf = interior(half, near), interior(half, far)
                        res["h2c_one_larger" if (pf and not pn) else "h2c_other"] += 1
    return res


def accumulation(P, Q, eps, src_j, cap=200000,
                 rungs=(1000, 10000, 100000, 200000)):
    """ARM 3 -- IS THE `Z`-`Z` CONNECTION IN THE CLOSURE OF THE MINIMAL COMPONENT?

    (user-raised: the (G0b) witness sits on the repo's anomaly row, and [NCYL-321] proves
    that row non-CP VIA its minimal component -- so ask whether the two objects touch.)

    A never-closing corner prong is dense in the minimal component `M` ([NCYL-002]'s `2` of
    `6`; [NCYL-210] `susp(Sigma_11)`, uniquely ergodic).  So for a corner point `y_j'`:

        `y_j'` in closure(M)  =>  min-distance over the orbit -> 0 as the orbit lengthens;
        `y_j'` in a CYLINDER  =>  min-distance PLATEAUS at a positive value.

    ⇒ reporting the min distance ACROSS RUNGS is the discriminant, not its value at one
    cap -- [NCYL-002]'s own *FLAT across four decades* methodology.  Distances are taken
    only at states sitting on `y_j'`'s OWN node, so this is not a bare 1-D coincidence.

    ⚠ [OPS-041] -- it can come out either way, and the internal contrast is the point: the
    same orbit is scored against the `Z` corners AND against the ESC corners on the same
    row.  If EVERY corner reads `-> 0` the arm is vacuous (the orbit would be dense in the
    whole cross-section, contradicting [NCYL-210]'s `3` cylinders of positive mass) and the
    verdict must be withheld.  ⚠ Float `.f` is used for the DISTANCE only; the exact `.v`
    decides node identity and every `eq` elsewhere in this file.
    """
    c = _chain(P, Q, eps)
    Lo, Hi, S, ops = c
    targets = {}
    for j in range(len(Lo) - 1):
        cs = nb.corner_start(Lo, Hi, S, j, ops)
        if cs is not None:
            targets[j] = (cs[0][0], cs[0][2])          # (node, y)
    src = nb.corner_start(Lo, Hi, S, src_j, ops)
    best = {j: [float("inf")] * len(rungs) for j in targets}
    st, n = src[0], 0
    for _ in range(cap):
        i, d, y = st
        for j, (node, ty) in targets.items():
            if i == node:
                dist = abs(y.f - ty.f)
                for r, R in enumerate(rungs):
                    if n < R and dist < best[j][r]:
                        best[j][r] = dist
        if not (0 <= i + d < len(Lo)):
            break
        nxt, kind, _e = sh.step(Lo, Hi, S, st, ops)
        if nxt is None:
            break
        st, n = nxt, n + 1
    return {"P": P, "Q": Q, "eps": eps, "src_j": src_j, "steps": n,
            "rungs": list(rungs),
            "min_dist": {str(j): [None if v == float("inf") else v for v in b]
                         for j, b in sorted(best.items())}}


def boundary_census(P, Q, eps, src_j, cap=200000, rungs=(1000, 10000, 100000, 200000)):
    """ARM 4 -- WHAT IS THE WHOLE OF `dM`?  An UNBIASED census, not a chosen shortlist.

    Arm 3 asked whether TWO chosen points sit in `closure(M)`.  This asks it of EVERY
    distinguished point of the cross-section at once -- every interface endpoint
    `Lo[i]`/`Hi[i]`, every POLE point `S[j]/2` ([NCYL-301]'s fold fixed point), and every
    CORNER launch point -- so the answer is a partition of the row's special points into
    `closure(M)` and the cylinders, with nothing pre-selected.

    ⇒⇒ THE QUESTION IT DECIDES.  If `dM` turns out to consist of the `Z`-`Z` connection
    and NOTHING ELSE, then in this family a minimal component is bounded by corner
    connections alone -- i.e. `M` exists <=> a `Z`-`Z` connection exists, which is the
    BICONDITIONAL that would convert the CP hypothesis into (G0b).  If instead POLE prongs
    also bound `M`, the connection is one boundary piece among several and no such
    conversion is available.  ⚠ EITHER ANSWER IS INFORMATIVE, and the second is the one
    that costs the route -- name it now so a null is not re-read as a delay ([OPS-041]).

    ⚠ NOT-EVIDENCE, stated in advance: `M` carries POSITIVE measure (`55.56%` of the
    boundary cross-section, [NCYL-082]), so a LOT of points will read `-> 0` and that is
    expected, not a finding.  The informative cells are the ones that PLATEAU, and the
    verdict is the PARTITION, never a single distance.  ⚠ A point that plateaus at a
    distance comparable to the orbit's own resolution is NOT decided -- report the rung
    trajectory, not the endpoint.
    """
    Lo, Hi, S, ops = _chain(P, Q, eps)
    tgt = {}
    for i in range(len(Lo)):
        tgt[("lo", i)] = (i, Lo[i])
        tgt[("hi", i)] = (i, Hi[i])
    for j in range(len(Lo) - 1):
        st = sh.pole_start(Lo, Hi, S, j, ops)
        tgt[("pole", j)] = (st[0], st[2])
        cs = nb.corner_start(Lo, Hi, S, j, ops)
        if cs is not None:
            tgt[("corner", j)] = (cs[0][0], cs[0][2])
    best = {k: [float("inf")] * len(rungs) for k in tgt}
    st, n = nb.corner_start(Lo, Hi, S, src_j, ops)[0], 0
    while n < cap:
        i, d, y = st
        for k, (node, ty) in tgt.items():
            if i == node:
                dist = abs(y.f - ty.f)
                for r, R in enumerate(rungs):
                    if n < R and dist < best[k][r]:
                        best[k][r] = dist
        if not (0 <= i + d < len(Lo)):
            break
        nxt, _kind, _e = sh.step(Lo, Hi, S, st, ops)
        if nxt is None:
            break
        st, n = nxt, n + 1
    return {"P": P, "Q": Q, "eps": eps, "src_j": src_j, "steps": n, "rungs": list(rungs),
            "min_dist": {f"{a}{b}": [None if v == float("inf") else v for v in vs]
                         for (a, b), vs in sorted(best.items())}}


def piece_table(rows=((7, 15, 0), (7, 15, 1), (4, 15, 0), (4, 15, 1)), cap=200000):
    """ARM 5 (user-requested) -- `7/15` vs `4/15`, PIECE BY PIECE.

    `4/15` is [NCYL-321]'s named negative control for `8/15`: same `Q`, same
    `class_role = lone`, same `perp_side = L1`, same `6` separatrices -- and `6/6` CLOSED,
    so "lone + L1 + 6 separatrices" does not force the minimal-component signature.  Both
    `min(P,Q-P)` values are >= 3, so BOTH rows are non-Veech ([NCYL-327]'s turn law:
    `turns = 6` at `7/15`, `turns = 3` at `4/15`) -- the comparison is inside the open
    locus, not across the Veech boundary.

    Tabulates, per class row: the role, the interface-length profile (the object whose
    extrema are the turns), and EVERY prong of both families -- the `Q-1)/2` CORNER prongs
    ([NCYL-002]'s separatrices) and the pole prongs -- with fate, steps and folds.

    ⚠ NOT-EVIDENCE.  The `4/15` side is a control, so its `all ESC` is EXPECTED and is not
    a finding ([OPS-041]); what is informative is any piece that differs in a way the
    CP/non-CP split does not already predict.  ⚠ A `CAP` verdict is not "never closes" --
    it is "not resolved at this cap" ([NCYL-262]: a cap is never a verdict); only
    [NCYL-002]'s four-decade flatness licenses the stronger reading, and only for the
    corner prongs it measured.
    """
    out = []
    for P, Q, eps in rows:
        c = _chain(P, Q, eps)
        if c is None:
            out.append({"P": P, "Q": Q, "eps": eps, "chain": "does not close"})
            continue
        Lo, Hi, S, ops = c
        prof = [float(_lenf(Lo, Hi, i)) for i in range(len(Lo))]
        rec = {"P": P, "Q": Q, "eps": eps,
               "role": "lone" if (eps % 2) == (P % 2) else "paired",
               "m": len(Lo) - 1, "nodes": len(Lo),
               "min_P": min(P, Q - P), "turns": min(P, Q - P) - 1,
               "profile": [round(x, 6) for x in prof],
               "corner": [], "pole": []}
        for j in range(len(Lo) - 1):
            cs = nb.corner_start(Lo, Hi, S, j, ops)
            if cs is None:
                rec["corner"].append({"j": j, "fate": "tie/degenerate"})
            else:
                k, _st, n, f = _walk(Lo, Hi, S, ops, cs[0], cap=cap)
                rec["corner"].append({"j": j, "fate": k, "steps": n, "folds": f})
            k, _st, n, f = _walk(Lo, Hi, S, ops,
                                 sh.pole_start(Lo, Hi, S, j, ops), cap=cap)
            rec["pole"].append({"j": j, "fate": k, "steps": n, "folds": f})
        for fam in ("corner", "pole"):
            tally = {}
            for e in rec[fam]:
                tally[e["fate"]] = tally.get(e["fate"], 0) + 1
            rec[fam + "_tally"] = tally
        out.append(rec)
    return out


def _lenf(Lo, Hi, i):
    return Hi[i].f - Lo[i].f


def cap_escalation(P=7, Q=15, eps=1, rungs=(200000, 2000000, 20000000),
                   fams=(("pole", 0), ("pole", 1), ("pole", 2), ("pole", 5),
                         ("corner", 1), ("corner", 3), ("corner", 2), ("corner", 5))):
    """ARM 6 -- THE CAP ESCALATION THAT MAKES ARM 4's BOUNDARY CLAIM QUOTABLE.

    A cap is never a verdict ([NCYL-262]).  [NCYL-002] licenses "never closes" at
    `1e3 -> 2e7` for the CORNER prongs it measured; the POLE prongs were capped at `2e5`
    by arms 3/4 and were undecided, and they are exactly the ones arm 4 first mis-read as
    BOUNDING `M`.  Re-runs them to [NCYL-002]'s own range.

    ⇒ Records the per-cap LADDER, not just the final rung: an escalation that reports only
    its last cap cannot supply the plateau shape, and the fold trajectory across decades is
    what distinguishes an orbit dense in `M` from one about to close.  ⚠ It RESTARTS rather
    than resuming (the walk state is one `(i, d, y)` triple and the whole ladder is ~20 s),
    so the `Sigma caps` cost is real but negligible here -- do NOT copy this shape to a run
    where a rung is expensive.
    """
    Lo, Hi, S, ops = _chain(P, Q, eps)
    out = []
    for fam, j in fams:
        st = (sh.pole_start(Lo, Hi, S, j, ops) if fam == "pole"
              else nb.corner_start(Lo, Hi, S, j, ops)[0])
        ladder, n, folds, fate = [], 0, 0, None
        for R in rungs:
            while n < R:
                i, d, _y = st
                if not (0 <= i + d < len(Lo)):
                    fate = "ESC"
                    break
                nxt, kind, _e = sh.step(Lo, Hi, S, st, ops)
                folds += int(kind == "fold")
                n += 1
                if nxt is None:
                    fate = kind
                    break
                st = nxt
            ladder.append({"cap": R, "steps": n, "folds": folds, "fate": fate or "CAP"})
            if fate:
                break
        out.append({"fam": fam, "j": j, "ladder": ladder, "fate": fate or "CAP"})
    return out


def verdict(res):
    """C2 reversal + the headline split."""
    by = {(r["P"], r["Q"], r["eps"], r["from_j"]): r for r in res["z_rows"]}
    rev_ok, rev_bad = 0, []
    for k, r in by.items():
        P, Q, eps, j = k
        t = r["to_j"]
        back = by.get((P, Q, eps, t))
        if back is not None and back["to_j"] == j:
            rev_ok += 1
        else:
            rev_bad.append({**r, "back_to": None if back is None else back["to_j"]})
    n_self = sum(int(r["self"]) for r in res["z_rows"])
    return {"c2_reversal_ok": rev_ok, "c2_bad": rev_bad,
            "SELF": n_self, "DISTINCT": len(res["z_rows"]) - n_self}


if __name__ == "__main__":
    qmax = int(sys.argv[1]) if len(sys.argv) > 1 else QMAX
    r = census(qmax)
    v = verdict(r)
    out = {"census": r, "verdict": v, "step2": step2_hypothesis(qmax),
           "accumulation": [accumulation(P, Q, eps, j)
                            for (P, Q, eps) in ((8, 15, 0), (7, 15, 1))
                            for j in (1, 3)],
           "boundary": [boundary_census(8, 15, 0, 1), boundary_census(8, 15, 0, 3)]}
    print(json.dumps(out, indent=1, sort_keys=True))
    with open("data/s492_corner_selfhit.json", "w") as fh:
        json.dump(out, fh, indent=1, sort_keys=True)
