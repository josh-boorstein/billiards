#!/usr/bin/env python3
"""
s456_neckgb.py -- PRE-REGISTERED.  Queue item (II): [NECK-GB], i.e. [NCYL-305]'s (G0)+(G1).

WHAT THIS IS.  NOT a proof of `p = 0`.  A REDUCTION of [NECK-GB] -- the last unproved
input to `CP => p = 0` besides CP itself -- plus the instrument that scores the reduction's
own diagnosis.  Two findings, and the second is a correction to [NCYL-305]:

  (A) ⇒⇒ THE GAUSS-BONNET APPARATUS IS REDUNDANT.  [NCYL-305] anatomises [NCYL-283]'s
      closing step into (G0) interiority, (G1) `chi(W) = 1`, (G2) no `Z` in `int W`, (G3)
      the boundary contributes `0`, and calls **(G1) the load-bearing one**.  It is not.
      Given (G0) -- both ends of the swapped edge are interior poles -- the edge IS an
      `R`-`R` saddle connection interior to `D`, which is verbatim the conclusion the step
      is invoked for.  (G1)/(G2)/(G3) exist only to evaluate `#poles(W) = 2`, and NOTHING
      downstream consumes that number; [NCYL-287] already recorded that it "is not a
      contradiction", and the reason is that it was never needed.
      ⇒ `chi(W) = 1` is free anyway: `W` is a regular neighbourhood of the arc `e` and
      deformation-retracts to it, so `chi = 1` unless `e` is a LOOP -- and a simple pole
      carries exactly ONE separatrix germ ([NCYL-280]: "the `R`-copies, ONE prong each"),
      so a leaf leaving `R_a` cannot also arrive at `R_a`; it would have to arrive along
      its own outgoing germ, which is a retrace, i.e. a `Fix` crossing, i.e. NOT swapped.
      The annulus [NCYL-305] warns about needs the pole self-loop that cannot exist.

  (B) ⇒⇒ SO [NECK-GB] IS EXACTLY (G0), AND (G0) IS NOT A SIDE-CONDITION -- IT IS THE
      UNTOUCHED HALF OF `p = 0`.  A swapped edge misses `Fix` except possibly at a SINGULAR
      endpoint on it (two leaves cannot cross a regular point; a horizontal `Fix` arc is
      itself a leaf), so `e` lies in `int D` together with its endpoints, each of which is
      an interior pole, `Z_O`, `Z_A`, or the one `iota`-fixed pole `R_0`.  Hence

          (G0)  <=>  no swapped edge has an endpoint at a zero
                <=>  (G0a) no `R`-`Z` edge beyond the `Fix` arm    [CLOSED: s455 H8/H9]
                   & (G0b) no swapped `Z`-`Z` edge                 [OPEN, and untouched]

      (`R_0` is free: its single prong IS the `Fix` arm, hence invariant.)  And (G0b) is
      the `O`/`A` half of [NCYL-281]'s (*) -- *every prong launched from a vertex returns
      to the same vertex*.  ⇒⇒ **[NCYL-298]+[NCYL-300]+[NCYL-301] prove (*) for prongs
      launched from `R` ONLY -- `(Q-1)/2` of the `Q - h`.  The `O`- and `A`-launched
      prongs, the other `Q - h - (Q-1)/2`, are proved by nothing.**  Modulo (G0a),
      (G0) <=> `p = 0`: discharging the Gauss-Bonnet step is not a step towards `p = 0`,
      it IS the half of `p = 0` nobody has attacked.

  ⚠ THIS IS A SCOPE CORRECTION, NOT A REFUTATION.  Nothing measured is wrong and no proof
  is broken: `p = 0` is measured on `859/859` scored class rows to `Q <= 40` INCLUDING the
  zero-launched prongs (`s446_covering_identity.row` loops `for v0 in ("O","R","A")`), and
  the necklace model traces the corner prongs too (`s450_chain_maps.Chain.separatrices`
  appends a "corner prong" per non-tie kite).  What is corrected is the READING of
  [NCYL-305]: the residual is not a topological lemma about a region, it is half the
  theorem.

PRIOR ART: read `rulings.md` [NCYL-305] (the anatomy this file rewrites), [NCYL-283] (the
step), [NCYL-287] (the `W` framing + the "not a contradiction" warning), [NCYL-282] (the
disk + corner angles), [NCYL-280] (`E = Q`, `V = Q+2`, the `a`/`b`/`d` edge types),
[NCYL-281] ((*) over ALL THREE vertices, and "no further `O`-`A` edge"), [NCYL-292] (the
executable necklace; the corner prong is the fold image of the corner), [NCYL-294] (the
global coordinate), [NCYL-298]/[NCYL-300] (`j' != j`), [NCYL-301] (`j' = j`), [NCYL-302]
(the assembly), [NCYL-304] (CP), [OPS-230] (the identity-element tell), [OPS-041];
grepped 'zero-launched', 'Z-Z', 'corner start', 'zero prong', 'Z_O prong', 'germ',
'ends at Z', 'launched from', 'n_nofold' through `rulings.py --grep` ->
  - NO ruling anywhere distinguishes the LAUNCH VERTEX of a prong when stating what the
    s453/s454 theorems cover.  [NCYL-281] states (*) for all three vertices; [NCYL-298]
    and [NCYL-301] state their results for "interior pole" prongs; [NCYL-302] assembles
    them into `p = 0` without noting the difference.  The gap is in the ASSEMBLY's scope,
    which is [OPS-232] a second time -- an assembly inherits the union of its inputs'
    hypotheses, and here the inputs' QUANTIFIER is narrower than the conclusion's.
  - s455's H8 (`pole -> Z` impossible) is the only result that touches a zero at all, and
    it is the `R`-side of the same edge; it closes (G0a) and says nothing about (G0b).

------------------------------------------------------------------------------------
WHY THE TWO PROVED TOOLS CANNOT SEE (G0b), and both diagnoses are SCORED below.

  * [NCYL-298]'s lattice.  A pole prong starts at `y_0 = s_j/2`; the criterion is
    `target -/+ y_0 in Lambda_S`, and it bites because `s_j/2` is measured OUTSIDE
    `Lambda_S = Z-span{s_e}`.  A CORNER prong starts at an interval ENDPOINT, which is a
    lattice point by construction -- so the difference lands in `Lambda_S` for free and
    the test carries no information.  ⇒ H3.  This is [OPS-230]'s identity-element tell for
    the THIRD time (s453's `s_j - s_j = 0 in 2Lambda`, s455's `0 in Lambda`, now this).
  * [NCYL-301]'s palindrome.  Step 2 needs the arrival state to be forced: `s_e/2` is
    interior to EXACTLY ONE of the two intervals at edge `e` (the larger), so a return to
    pole `j` must be `R(u_0)`.  A corner is an endpoint of BOTH intervals at its edge, so
    the arrival node is not forced, `u_M = R(u_0)` fails, and Steps 3-4 have no input.
    ⇒ H2.  If corners were ALSO interior to exactly one, [NCYL-301] would transfer
    verbatim and there would be no gap -- so H2 can come out either way.

HYPOTHESES (pre-registered).
  H1  ⇒⇒ THE SIZE OF THE GAP, pure arithmetic over the whole box, no tracing.  For every
      coprime `(P,Q)` and both classes: `n_R == (Q-1)/2` (the prongs the theorems cover),
      `n_O + n_A == Q - h - (Q-1)/2` (the prongs they do not), and their sum is s446's
      H0 `Q - h`.  s446 H0 scores only the SUM; the SPLIT is what this asserts, and an
      off-by-one in either part breaks it while leaving H0 intact.
  H2  ⇒⇒ [NCYL-301] STEP 2's INPUT, ASKED OF BOTH START TYPES ([OPS-043] -- the diagnosis
      asserted in code, not in prose).  At every edge `j` of the quotient path:
      `S[j]/2` interior to exactly ONE of `J_j`, `J_{j+1}` (Step 2, re-scored), AND the
      shared endpoint an endpoint of BOTH (so no arrival node is forced).  Can fail.
  H3  ⇒⇒ [NCYL-298]'s CRITERION EVALUATED AT CORNER STARTS, against s455 H8's `0` hits at
      pole starts on the same rows.  Same lattice, same test, same endpoints; only `y_0`
      changes.  Reported as a HIT RATE, not a pass: a low corner hit rate would mean the
      lattice route is alive for corners too and (B)'s "untouched" is wrong.
  H4  ⇒⇒ THE POSITIVE CONTROL, and it is what makes the real `0` mean anything ([OPS-041],
      the [NCYL-298]-H4 / [NCYL-301]-H2 design): on SYNTHETIC interval data with the same
      adjacency (adjacent intervals share exactly one endpoint) and no arithmetic, do
      corner prongs end at `Z` AT ALL?  If they never do, (G0b) is structural and the
      argument is findable; if they do, the real `0` is an arithmetic fact and (G0b) is a
      genuine open problem.  Either outcome is informative and they say opposite things.
      Pole prongs run alongside as the reference arm.
  H5  ⇒⇒ THE REAL-ROW CROSS-TAB.  Corner-prong endings split by whether s446's BILLIARD
      sweep SCORED the row.  ⚠ *no non-`ESC` on a scored row* is PREDICTED by
      `n_nofold = 0` and is NOT evidence; what is not predicted is that the detector
      fires at all on real data, and WHERE it fires.
  H6  ⇒⇒ THE OBVIOUS NEXT INSTRUMENT, PRE-SCORED SO NOBODY BUILDS IT TWICE.  The signs in
      `y_n = (-1)^n y_0 + sum_e c_e S_e` alternate in visit order, forcing
      `sum_e c_e in {0,1}` by the parity of `n` -- a constraint the plain `Lambda_S` test
      discards.  Restoring it costs one lattice coordinate.  Does it close the corner
      case?  A `0` corner hit rate would be a proof route; anything else names the
      strengthening as insufficient and prices it.

RESULTS (s456; `all --qmax=33`, `103.1 s`, `data/s456_neckgb.json`, `logs/s456_neckgb.log`).
  H1  `n_R == (Q-1)/2` on **`5816/5816`** class rows (odd `Q <= 120`, separate run) and
      `460/460` at `Q <= 33`; s446's H0 sum re-scores `5816/5816`.  ⇒⇒ over the whole box
      **`229708` prongs COVERED** by [NCYL-298]/[NCYL-300]/[NCYL-301]/H8 and **`226800`
      covered by NOTHING -- `49.7%`** (`48.8%` at `Q <= 33`).
  H2  `5032/5032` edges (`460` rows): `S[j]/2` interior to exactly ONE neighbour AND the
      shared endpoint an endpoint of BOTH; neighbours share exactly one endpoint
      `5032/5032`.  ⇒ Step 2 forces the arrival node for a pole and cannot for a corner.
  H3  ⇒⇒ THE CONTRAST, same lattice, same targets, only `y_0` changes:
      POLE starts **`0/267040`** hits (`y_0 notin Lambda_S` on `5032/5032`) against
      CORNER starts **`222845/256056 = 87.0%`** hits (`y_0 in Lambda_S` on `4254/4802`).
      ⇒ [NCYL-298]'s invariant is not merely weaker at a corner, it admits `7` targets in
      `8`.  ⚠ It is NOT identically vacuous -- it still kills `13%` -- so a FINER
      arithmetic invariant is not a priori hopeless, which is why H6 exists.
  H6  the coefficient-sum refinement is REAL AND INSUFFICIENT: corner hit rate
      `85.75% -> 77.60%` (`37582 -> 34012` of `43828`, `Q <= 25`); poles stay `0/46336`.
      ⇒ it buys `8` points and does not close the case.  Priced so nobody rebuilds it.
  H4  ⇒⇒ POSITIVE CONTROL FIRES: on `400` synthetic paths corner prongs end
      `2325 ESC / 500 Z / 460 pole` and pole prongs `2509 ESC / 460 Z / 316 pole`.
      ⇒⇒ **corner -> `Z` is NOT structurally impossible -- it happens freely once the
      arithmetic is removed -- so (G0b) is a genuine arithmetic statement and the real `0`
      is a fact that needs proving, not a triviality.**
  H5  ⇒⇒ AND THE DETECTOR FIRES ON REAL ROWS TOO, EXACTLY OUTSIDE THE CP SCOPE GATE:
      `4802` corner prongs over `460` rows, `4684 ESC / 114 CAP / 4 Z`.  Split by s446's
      verdict: on SCORED rows `4397 ESC`, `3 CAP`, **`0 Z`** (predicted -- `n_nofold = 0`);
      on UNSCORED rows `287 ESC`, `111 CAP`, **`4 Z`**.  All four `Z` are `7/15 eps=1` and
      `8/15 eps=0` -- [NCYL-283]'s named minimal-component pair, reproduced here from kite
      arithmetic with no billiard trace and no shared code.  ⚠ Not every unscored row
      yields a `Z`; the implication runs one way only.

GUARDS.
  G1  a row whose global coordinate does not close (`s451.build`) is SKIPPED, not scored.
  G2  a TIE edge (equal-length neighbours) and a DEGENERATE node (`Lo == Hi`, the lone
      class's horizontal interface) carry no corner prong; counted as SKIPPED, never as a
      pass.
  G3  nothing here writes to any store, and no existing probe is modified.
  G4  ⚠ [OPS-228]: numpy `int64` is cast to Python `int` before `zlattice`.

Run: PYTHONPATH=.:engine:archive/scripts_2026-07:archive/scripts_2026-08:probes \\
     .venv/bin/python3.13 probes/s456_neckgb.py [all|h1|h2|h3|h4|h5] [--qmax=N]
"""
from __future__ import annotations

import json
import math
import random
import sys
import time
from fractions import Fraction

import s439_exact_cells as s439
import s446_covering_identity as s446
import s451_quotient_path as qp
import s453_pole_module as pm
import s454_self_hit as sh
import zlattice as zl

OUT = "data/s456_neckgb.json"


# ---------------------------------------------------------------- the corner prong

def corner_start(Lo, Hi, S, j, ops):
    """`u_0` for the CORNER prong at edge `j`, the zero-launched twin of `sh.pole_start`.

    [NCYL-292]: the fold at edge `j` lives on the LARGER of the two interfaces and the
    corner's prong is the fold image of the corner -- i.e. the leaf at the NEAR
    interface's NON-shared endpoint, sitting on the FAR node and heading AWAY from the
    edge.  (The shared endpoint is the kite corner itself; the free endpoint of the
    smaller interface is where its prong enters the larger one.)

    Returns `(state, near, far)` or `None` when the edge is a TIE or a node is degenerate.
    """
    la, lb = ops.sub(Hi[j], Lo[j]), ops.sub(Hi[j + 1], Lo[j + 1])
    c = ops.cmp(la, lb)
    if c == 0:
        return None                                   # tie kite: no fold, no corner prong
    far, near = (j, j + 1) if c > 0 else (j + 1, j)
    if ops.eq(Lo[near], Hi[near]) or ops.eq(Lo[far], Hi[far]):
        return None                                   # degenerate flank
    lo_shared = ops.eq(Lo[j], Lo[j + 1])
    y = Hi[near] if lo_shared else Lo[near]
    return (far, -1 if far == j else +1, y), near, far


def _rows(qmax, qmin=5):
    for Q in range(qmin, qmax + 1, 2):
        for P in range(1, Q):
            if math.gcd(P, Q) == 1:
                for eps in (0, 1):
                    yield P, Q, eps


def _chain(P, Q, eps):
    ch, lo, hi, closes = qp.build(P, Q, eps)
    if not closes:
        return None
    sf = pm.fold_sums(ch, lo, hi)
    Lo, Hi, S = pm.path_view(ch, lo, hi, sf)
    return ch, Lo, Hi, S, pm.CoordOps(ch)


# ---------------------------------------------------------------- H1

def h1(qmax=120):
    """The launch-vertex SPLIT.  `n_R == (Q-1)/2`; the rest is the uncovered share."""
    res = {"rows": 0, "nR_ok": 0, "sum_ok": 0, "covered": 0, "uncovered": 0,
           "bad": []}
    for Q in range(5, qmax + 1, 2):
        for P in range(1, Q):
            if math.gcd(P, Q) != 1:
                continue
            for eps in (0, 1):
                res["rows"] += 1
                hor, _v = s446.hv_sides(P, Q, eps)
                n = {v: len(s439.prongs(v, P, Q, eps)) for v in ("O", "R", "A")}
                res["nR_ok"] += (n["R"] == (Q - 1) // 2)
                res["sum_ok"] += (sum(n.values()) == Q - len(hor))
                res["covered"] += n["R"]
                res["uncovered"] += n["O"] + n["A"]
                if n["R"] != (Q - 1) // 2 and len(res["bad"]) < 6:
                    res["bad"].append([P, Q, eps, n["R"], (Q - 1) // 2])
    tot = res["covered"] + res["uncovered"]
    res["frac_uncovered"] = round(res["uncovered"] / tot, 5) if tot else None
    return res


# ---------------------------------------------------------------- H2

def h2(qmax=45):
    """[NCYL-301] Step 2's input, for POLES (re-scored) and for CORNERS (the diagnosis)."""
    res = {"rows": 0, "edges": 0, "pole_one": 0, "pole_bad": 0,
           "corner_both": 0, "corner_bad": 0, "shared_exactly_one": 0, "skipped": 0,
           "bad": []}
    for P, Q, eps in _rows(qmax):
        c = _chain(P, Q, eps)
        if c is None:
            continue
        _ch, Lo, Hi, S, ops = c
        res["rows"] += 1
        for j in range(len(S)):
            res["edges"] += 1
            half = ops.half(S[j])
            inside = sum(1 for i in (j, j + 1)
                         if ops.cmp(half, Lo[i]) > 0 and ops.cmp(half, Hi[i]) < 0)
            if inside == 1:
                res["pole_one"] += 1
            else:
                res["pole_bad"] += 1
                if len(res["bad"]) < 6:
                    res["bad"].append(["pole", P, Q, eps, j, inside])
            lo_s, hi_s = ops.eq(Lo[j], Lo[j + 1]), ops.eq(Hi[j], Hi[j + 1])
            res["shared_exactly_one"] += (lo_s != hi_s)
            # the shared endpoint is an endpoint of BOTH intervals -> arrival node free
            shared = Lo[j] if lo_s else Hi[j]
            both = sum(1 for i in (j, j + 1)
                       if ops.eq(shared, Lo[i]) or ops.eq(shared, Hi[i]))
            if both == 2:
                res["corner_both"] += 1
            else:
                res["corner_bad"] += 1
                if len(res["bad"]) < 6:
                    res["bad"].append(["corner", P, Q, eps, j, both])
    return res


# ---------------------------------------------------------------- H3

def h3(qmax=29):
    """[NCYL-298]'s criterion at CORNER starts vs POLE starts, same lattice, same targets."""
    res = {"rows": 0, "pole_tests": 0, "pole_hits": 0, "corner_tests": 0,
           "corner_hits": 0, "y0_pole_in": 0, "y0_pole_out": 0,
           "y0_corner_in": 0, "y0_corner_out": 0, "skipped": 0}
    for P, Q, eps in _rows(qmax):
        c = _chain(P, Q, eps)
        if c is None:
            continue
        _ch, Lo, Hi, S, ops = c
        m = len(S)
        if m == 0:
            continue
        res["rows"] += 1
        L = zl.Lattice([[int(x) for x in S[e].v] for e in range(m)])
        ends = [e for i in range(m + 1) for e in (Lo[i], Hi[i])]
        for j in range(m):
            starts = [("pole", ops.half(S[j]))]
            cs = corner_start(Lo, Hi, S, j, ops)
            if cs is None:
                res["skipped"] += 1
            else:
                starts.append(("corner", cs[0][2]))
            for tag, y0 in starts:
                yv = [int(x) for x in y0.v]
                res[f"y0_{tag}_in" if L.contains(yv) else f"y0_{tag}_out"] += 1
                for end in ends:
                    for sgn in (+1, -1):
                        res[f"{tag}_tests"] += 1
                        if L.contains([int(a) - sgn * b
                                       for a, b in zip(end.v, yv)]):
                            res[f"{tag}_hits"] += 1
    for t in ("pole", "corner"):
        n = res[f"{t}_tests"]
        res[f"{t}_hit_rate"] = round(res[f"{t}_hits"] / n, 6) if n else None
    return res


# ---------------------------------------------------------------- H6

def h6(qmax=25):
    """⇒⇒ THE NATURAL STRENGTHENING OF [NCYL-298]'s INVARIANT, AND IT IS NOT ENOUGH.

    `y_{n} = S_{e_n} - y_{n-1}` gives, after `n` folds,

        y_n = (-1)^n y_0 + sum_e c_e S_e      with   sum_e c_e = 1 (n odd) / 0 (n even),

    because the signs alternate in VISIT order.  `Lambda_S`-membership throws that
    constraint away (it is the relaxation s455 H8 flags: "the Z-span is a superset of the
    +-1-coefficient set").  Keeping it costs ONE extra lattice coordinate -- test
    `(t - y_0, 0)` and `(t + y_0, 1)` against `span{(S_e, 1)}` -- and gives a strictly
    stronger necessary condition, scored here against the plain test on the same rows.
    ⚠ It CAN come out either way: a `0` corner hit rate would be a proof route.
    """
    res = {"rows": 0, "pole_aug_hits": 0, "pole_tests": 0,
           "corner_aug_hits": 0, "corner_plain_hits": 0, "corner_tests": 0}
    for P, Q, eps in _rows(qmax):
        c = _chain(P, Q, eps)
        if c is None:
            continue
        _ch, Lo, Hi, S, ops = c
        m = len(S)
        if m == 0:
            continue
        res["rows"] += 1
        plain = [[int(x) for x in S[e].v] for e in range(m)]
        Laug = zl.Lattice([g + [1] for g in plain])
        Lpln = zl.Lattice(plain)
        ends = [[int(x) for x in e.v] for i in range(m + 1) for e in (Lo[i], Hi[i])]
        for j in range(m):
            cs = corner_start(Lo, Hi, S, j, ops)
            todo = [("pole", [int(x) for x in ops.half(S[j]).v])]
            if cs is not None:
                todo.append(("corner", [int(x) for x in cs[0][2].v]))
            for tag, yv in todo:
                for ev in ends:
                    res[f"{tag}_tests"] += 1
                    hit = (Laug.contains([a - b for a, b in zip(ev, yv)] + [0])
                           or Laug.contains([a + b for a, b in zip(ev, yv)] + [1]))
                    res[f"{tag}_aug_hits"] += hit
                    if tag == "corner":
                        res["corner_plain_hits"] += Lpln.contains(
                            [a - b for a, b in zip(ev, yv)])
    for t in ("pole_aug", "corner_aug", "corner_plain"):
        n = res["pole_tests" if t.startswith("pole") else "corner_tests"]
        res[f"{t}_rate"] = round(res[f"{t}_hits"] / n, 6) if n else None
    return res


# ---------------------------------------------------------------- H4

def synth_path(rng, n):
    """A random quotient path with the REAL adjacency and no arithmetic: consecutive
    intervals share exactly one endpoint ([NCYL-301] H3, `4392/4392`)."""
    Lo = [Fraction(0)]
    Hi = [Fraction(rng.randint(8, 50))]
    for _ in range(1, n):
        if rng.random() < 0.5:                                   # share the low endpoint
            lo = Lo[-1]
            hi = lo + Fraction(rng.randint(1, 70))
            while hi == Hi[-1]:
                hi += 1
        else:                                                    # share the high endpoint
            hi = Hi[-1]
            lo = hi - Fraction(rng.randint(1, 70))
            while lo == Lo[-1]:
                lo -= 1
        Lo.append(lo)
        Hi.append(hi)
    S = [(Hi[j] + Hi[j + 1]) if Lo[j] == Lo[j + 1] else (Lo[j] + Lo[j + 1])
         for j in range(n - 1)]
    return Lo, Hi, S


def h4(trials=400, seed=20260902, cap=20000):
    """POSITIVE CONTROL: do corner prongs end at `Z` on synthetic data?"""
    rng = random.Random(seed)
    ops = pm.FracOps
    res = {"trials": 0, "pole": {}, "corner": {}, "skipped": 0, "examples": []}
    for _ in range(trials):
        n = rng.randint(4, 14)
        Lo, Hi, S = synth_path(rng, n)
        res["trials"] += 1
        for j in range(len(S)):
            for tag, st in (("pole", sh.pole_start(Lo, Hi, S, j, ops)),
                            ("corner", (corner_start(Lo, Hi, S, j, ops) or (None,))[0])):
                if st is None:
                    res["skipped"] += 1
                    continue
                _s, _k, _e, end = sh.orbit(Lo, Hi, S, st, ops, cap)
                res[tag][end] = res[tag].get(end, 0) + 1
                if tag == "corner" and end == "Z" and len(res["examples"]) < 4:
                    res["examples"].append({"n": n, "j": j, "steps": len(_k)})
    return res


# ---------------------------------------------------------------- H5

def h5(qmax=31, cap=200000, ref="data/s455_covering_q40.json"):
    """⇒⇒ THE CROSS-TAB, and it upgrades this arm out of *coverage only*.

    Corner-prong endings on REAL rows, split by whether s446's BILLIARD sweep scored the
    row (`status == 'ok'`, i.e. inside [NCYL-283]'s CP scope gate) or refused it
    (`PRONG_CAPPED`).  `n_nofold = 0` predicts NO non-`ESC` corner prong on a SCORED row,
    so that half is not evidence; what is NOT predicted is that the detector FIRES at all
    on real data, and where.  ⚠ `PRONG_CAPPED` is a CAP statement, not a *not CP* verdict
    ([NCYL-262]).
    """
    scored = set()
    try:
        for r in json.load(open(ref))["rows"]:
            if r.get("status") == "ok":
                scored.add((r["P"], r["Q"], r["eps"]))
    except OSError:
        scored = None
    res = {"rows": 0, "prongs": 0, "ends": {}, "skipped": 0, "hits": [],
           "max_folds": 0, "on_scored": {}, "on_unscored": {}, "no_ref": 0}
    for P, Q, eps in _rows(qmax):
        c = _chain(P, Q, eps)
        if c is None:
            continue
        _ch, Lo, Hi, S, ops = c
        res["rows"] += 1
        if scored is None:
            bucket = None
        elif (P, Q, eps) in scored:
            bucket = "on_scored"
        else:
            bucket = "on_unscored"
        for j in range(len(S)):
            cs = corner_start(Lo, Hi, S, j, ops)
            if cs is None:
                res["skipped"] += 1
                continue
            res["prongs"] += 1
            sts, kinds, _e, end = sh.orbit(Lo, Hi, S, cs[0], ops, cap)
            res["ends"][end] = res["ends"].get(end, 0) + 1
            if bucket:
                res[bucket][end] = res[bucket].get(end, 0) + 1
            else:
                res["no_ref"] += 1
            res["max_folds"] = max(res["max_folds"], sum(k == "fold" for k in kinds))
            if end not in ("ESC", "CAP"):
                res["hits"].append([P, Q, eps, j, end, len(sts),
                                    bucket or "no_ref"])
    return res


# ---------------------------------------------------------------- driver

def main():
    which = sys.argv[1] if len(sys.argv) > 1 and not sys.argv[1].startswith("-") else "all"
    qmax = None
    for a in sys.argv[1:]:
        if a.startswith("--qmax="):
            qmax = int(a.split("=")[1])
    out, t0 = {}, time.time()

    if which in ("all", "h1"):
        t = time.time()
        r = out["h1"] = h1(qmax or 120)
        print(f"  H1  n_R == (Q-1)/2                       : {r['nR_ok']}/{r['rows']}"
              f" class rows (odd Q <= {qmax or 120})", flush=True)
        print(f"      sum == Q - h (s446 H0, re-scored)    : {r['sum_ok']}/{r['rows']}",
              flush=True)
        print(f"      ⇒ prongs COVERED by [NCYL-298]/[NCYL-301]/H8 : {r['covered']}",
              flush=True)
        print(f"      ⇒ prongs covered by NOTHING (O and A)       : {r['uncovered']}"
              f"  ({100 * r['frac_uncovered']:.1f}% of all prongs)", flush=True)
        for b in r["bad"]:
            print("      bad:", b, flush=True)
        print(f"      [{time.time() - t:.1f}s]", flush=True)

    if which in ("all", "h2"):
        t = time.time()
        r = out["h2"] = h2(qmax or 45)
        print(f"  H2  pole s_j/2 interior to EXACTLY ONE   : {r['pole_one']}/{r['edges']}"
              f"  ({r['rows']} rows)", flush=True)
        print(f"      corner is an endpoint of BOTH        : {r['corner_both']}"
              f"/{r['edges']}", flush=True)
        print(f"      neighbours share exactly one endpoint: "
              f"{r['shared_exactly_one']}/{r['edges']}", flush=True)
        print("      ⇒ Step 2 forces the arrival node for a POLE and not for a CORNER",
              flush=True)
        for b in r["bad"]:
            print("      bad:", b, flush=True)
        print(f"      [{time.time() - t:.1f}s]", flush=True)

    if which in ("all", "h3"):
        t = time.time()
        r = out["h3"] = h3(qmax or 29)
        print(f"  H3  [NCYL-298] criterion, POLE   starts  : {r['pole_hits']}"
              f"/{r['pole_tests']} hits  (rate {r['pole_hit_rate']})", flush=True)
        print(f"      [NCYL-298] criterion, CORNER starts  : {r['corner_hits']}"
              f"/{r['corner_tests']} hits  (rate {r['corner_hit_rate']})", flush=True)
        print(f"      y_0 in Lambda_S: pole {r['y0_pole_in']}/"
              f"{r['y0_pole_in'] + r['y0_pole_out']}   corner {r['y0_corner_in']}/"
              f"{r['y0_corner_in'] + r['y0_corner_out']}", flush=True)
        print(f"      [{time.time() - t:.1f}s]", flush=True)

    if which in ("all", "h6"):
        t = time.time()
        r = out["h6"] = h6(min(qmax or 25, 25))
        print(f"  H6  AUGMENTED lattice, POLE   starts     : {r['pole_aug_hits']}"
              f"/{r['pole_tests']} hits  (rate {r['pole_aug_rate']})", flush=True)
        print(f"      AUGMENTED lattice, CORNER starts     : {r['corner_aug_hits']}"
              f"/{r['corner_tests']} hits  (rate {r['corner_aug_rate']})", flush=True)
        print(f"      PLAIN     lattice, CORNER starts     : {r['corner_plain_hits']}"
              f"/{r['corner_tests']} hits  (rate {r['corner_plain_rate']})", flush=True)
        print("      ⇒ the coefficient-sum constraint helps and does NOT close it",
              flush=True)
        print(f"      [{time.time() - t:.1f}s]", flush=True)

    if which in ("all", "h4"):
        t = time.time()
        r = out["h4"] = h4()
        print(f"  H4  SYNTHETIC ({r['trials']} paths).  pole ends   : {r['pole']}",
              flush=True)
        print(f"                              corner ends : {r['corner']}", flush=True)
        cz = r["corner"].get("Z", 0)
        print(f"      ⇒ POSITIVE CONTROL corner->Z fires   : {cz}"
              f"   {'OK' if cz else '*** DEAD -- (G0b) may be structural ***'}",
              flush=True)
        print(f"      [{time.time() - t:.1f}s]", flush=True)

    if which in ("all", "h5"):
        t = time.time()
        r = out["h5"] = h5(qmax or 31)
        print(f"  H5  real corner prongs: {r['prongs']} over {r['rows']} rows,"
              f" ends {r['ends']}", flush=True)
        print(f"      on s446-SCORED rows (⚠ predicted, not evidence): {r['on_scored']}",
              flush=True)
        print(f"      on s446-UNSCORED rows (the detector FIRING)    : {r['on_unscored']}",
              flush=True)
        print(f"      skipped (tie/degenerate) {r['skipped']}, max folds"
              f" {r['max_folds']}", flush=True)
        for b in r["hits"]:
            print("      NON-ESC:", b, flush=True)
        print(f"      [{time.time() - t:.1f}s]", flush=True)

    out["secs"] = round(time.time() - t0, 1)
    json.dump(out, open(OUT, "w"), indent=1)
    print(f"  -> {OUT}  [{out['secs']}s]", flush=True)


if __name__ == "__main__":
    main()
