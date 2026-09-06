#!/usr/bin/env python3
"""s469_profile_proof.py -- PRE-REGISTERED.  Queue items (s468-a) and (s468-c).

⇒⇒ WHAT THIS FILE IS.  A PROOF, plus the guards on its joints.  [NCYL-324] (1) records
*adjacent path interfaces nest and share exactly one endpoint* as MEASURED (`5032/5032`,
odd `Q <= 33`) and calls it "the ONE lemma between this line and a theorem".  It is not a
measurement at all: it is an IDENTITY, forced by the two lines that DEFINE `delta[k]` and
`u`, and once that is seen the whole extreme-`P` escape follows with no cap, no arithmetic
input beyond a sine comparison, and no bound on `Q`.

⚠⚠ EVERY SWEEP BELOW IS A GUARD ON A JOINT, NOT A GROUND ([NCYL-318] scope (i)/[OPS-041]).
The theorem does not rest on the row counts; the row counts catch a misread of the code.

--------------------------------------------------------------------------------------
THE PROOF.

L1 (THE NESTING IDENTITY -- definitional).  `s450_chain_maps.Chain.__init__` sets, for
   every kite `k`,   `far[k], near[k] = the two edges of k ordered by ell`  (line 232-235,
   `far` the LARGER), and then  `delta[k] = ell[far[k]] - ell[near[k]] >= 0`  (line 243).
   `s451_quotient_path.build` propagates the global left endpoint by, for each kite,
        `u[near] = u[far] + delta[k]`   if `osec[k]`,      `u[near] = u[far]`   if not
   (lines 156-161; the `else` branch writes the SAME relation solved for `u[far]`).  With
   `lo[t] = u[t]` and `hi[t] = u[t] + ell[t]` (lines 164-165):

        not osec[k]:  lo[near] = lo[far]                and  hi[near] = hi[far] - delta[k]
            osec[k]:  lo[near] = lo[far] + delta[k]     and  hi[near] = hi[far]

   ⇒ in BOTH cases  `J_near ⊆ J_far`, sharing the LO endpoint when `osec[k]` is false and
   the HI endpoint when it is true, with the containment STRICT iff `delta[k] > 0`.
   ⚠ This needs no closure verdict and no bound on `Q`: it is two substitutions.

L2 (THE TIE IS ALWAYS KITE `0`, AND THE PATH NEVER USES KITE `0`).  `delta[k] = 0` iff
   `ell[(k-1) % Q] = ell[k]`.  With `c[t] = (c0 + 2Pt) mod 2Q` and `ell[t] = sin(c[t]π/2Q)`,
   and `c[t] ∈ [0, 2Q)` so that `c[t]π/2Q ∈ [0, π)`, equality of sines forces
   `c[k] = c[k-1]` or `c[k] + c[k-1] = 2Q`.  The first is impossible (`c[k] - c[k-1] ≡ 2P`,
   `0 < 2P < 2Q`).  The second reads `c0 + P(2k-1) ≡ 0 (mod Q)`; since `c0 ≡ P (mod Q)` in
   both parities (`c0 = P` or `P + Q`), it becomes `2Pk ≡ 0 (mod Q)`, and `gcd(2P, Q) = 1`
   for odd `Q` gives `k ≡ 0`.  ⇒ the UNIQUE tie kite is `k = 0`.
   Path node `i` is interface `m - i` (`s453_pole_module.path_view`), so the adjacency
   `(i, i+1)` is kite `m - i` with `i ∈ [0, m-1]`, i.e. kites `1 .. m`.  Kite `0` is never
   a path adjacency.  ⇒ L1 applies with `delta > 0` at every path adjacency.

L3 (THE PROFILE IS A FOLDED ROTATION).  Write `w(x) = min(x, 2Q - x)`.  Then
   `|J_i| = ell[m-i] = sin(w(x_i)π/2Q)` with `x_i = (c0 + 2P(m-i)) mod 2Q`, so the profile
   is monotone-equivalent to `a_i := w(x_i)`, an arithmetic progression of step `2P` read
   through the fold at `{0, Q}`.  Three facts, each a one-liner:
     (a) `x_0 = c[m] = (c0 + P(Q-1)) mod 2Q ∈ {0, Q}` -- for `c0 = P`, `c[m] = PQ mod 2Q`,
         which is `Q` for odd `P` and `0` for even `P`; adding `Q` (the `eps = 1` case)
         swaps them.  ⇒ the walk STARTS on a fold point and `a_0` is the extreme value.
     (b) the value set `{a_i}` is `{1, 3, …, Q}` (odd `a_0 = Q`) or `{0, 2, …, Q-1}`
         (even, `a_0 = 0`) -- an AP of step `2` with `m+1` distinct terms.  Hence exactly
         ONE of the two fold points lies in the value set, so no INTERIOR `x_i` is a fold
         point.
     (c) `t ↦ Q-1-t` is the `w`-pairing of interfaces, and the path `t ∈ [0, m]` takes one
         representative of each pair -- which is why (b)'s terms are distinct.

L4 (THE TURN COUNT).  `turns` := the number of interior local extrema of `a`.  On each of
   the open arcs `(0, Q)` and `(Q, 2Q)`, `w` is strictly monotone, so a turn occurs exactly
   when a step crosses a fold point; the step length is `2p` with `p := min(P, Q-P) <= (Q-1)/2`,
   which is `< Q`, so a step crosses at most one.  Total travel is `2pm = p(Q-1)`, from a
   start ON a fold point, and the fold points are `Q` apart, so the number of crossings is
   `floor((pQ - p)/Q) = p - 1`.  Each crossing contributes its extremum to an INTERIOR
   index: it lands on whichever of the two samples is nearer the fold point, and
     - the FIRST crossing (travel `Q`) cannot land on `i = 0`, which sits a full `Q` away
       while `x_1` sits `Q - 2p` away;
     - the LAST crossing (travel `(p-1)Q`) lands on `i = m` only if `p(Q-2) >= (p-1)Q`,
       i.e. `2p >= Q`, which is false.
   ⇒⇒ `turns = min(P, Q-P) - 1`.  In particular `turns = 0` iff `P ∈ {1, Q-1}` and
   `turns = 1` iff `P ∈ {2, Q-2}`.

L5 (THE LAUNCH RULE).  `s456_neckgb.corner_start` puts the prong on the LARGER node of the
   edge, heading AWAY from the edge.  So its first step is uphill unless the launch node is
   a LOCAL MAXIMUM of the profile; and by induction the run stays uphill unless it meets a
   local maximum.  ⇒⇒ IF THE PROFILE HAS NO INTERIOR LOCAL MAXIMUM, every corner prong
   crosses (L1: uphill entry into a strictly larger interface is forced, and a cross leaves
   `y` unchanged) to the end of the path and ESCAPES, with ZERO folds, in exactly
   `j + 1` steps when the larger node is `j`, and `m - j` steps when it is `j + 1`.
   ⚠ The step count is keyed on WHICH NODE OF THE EDGE IS LARGER, not on a global
   direction: on a VALLEY both cases occur on the same row (the first draft of H5 used the
   row-level direction and mispredicted `1250` of `7546` prongs, every one of them a
   valley -- the walk outcome was right, the closed form was wrong).

⇒⇒⇒ THEOREM.  Let `Q >= 5` be odd, `gcd(P, Q) = 1`, `eps ∈ {0,1}`, and suppose the profile
has no interior local MAXIMUM.  Then every corner prong escapes with zero folds in at most
`m = (Q-1)/2` steps.  ⇒ this holds for `turns = 0`, i.e. `P ∈ {1, Q-1}` (both `eps`), and
for the VALLEY half of `turns = 1`, i.e. `P ∈ {2, Q-2}` with the interior extremum a
minimum.  ⇒ [NCYL-324] (3)'s `4998/4998` is PROVED, at every odd `Q`, not just `Q <= 101`.

⚠ WHAT IS NOT PROVED: the TENT half of `turns = 1` (one interior local MAXIMUM).  H6 scores
it -- it escapes too, but by a bounded fold cascade, and that bound is measured.

--------------------------------------------------------------------------------------
PRIOR ART.  `rulings.py --grep` on 'folded rotation', 'local maximum', 'three distance'
-> nothing; 'unimodal' -> [WFLOOR-013]/[WFLOOR-023], a DIFFERENT object (the width ladder
`w_m`), no contact; 'three-distance' -> [SEAM-006]/[WFLOOR-034]/[WFLOOR-037], all the graze
orbit, no contact with the chain model; 'turning point' -> [WFLOOR-151]'s ν-walk peaks,
again a different object.  ⇒ ID-grep per [OPS-246] on [NCYL-292]/[NCYL-297]/[NCYL-308]/
[NCYL-319]/[NCYL-324]: [NCYL-292] uses the nesting at the ONE extremal kite (`ell[m] = 1`
is the strict max) and [NCYL-324] generalises it to the whole profile -- BOTH state it as
measured; no line anywhere says it is definitional.  ⇒⇒ THE ONE REAL HIT, AND IT IS ON THE
(s468-c) HALF: [NCYL-012] already records `literature_review.md`'s Veech criterion for this
family as `min(P, Q-P) ∈ {1, 2}`.  That is L4's `turns <= 1` verbatim, so the census's
"Veech = extreme `P`" coincidence [NCYL-325] could not resolve is not a coincidence of two
properties -- it is ONE property with two names, and [NCYL-325] re-derived it from the
regular-`n`-gon match without citing [NCYL-012].  ⊕ [NCYL-062] + [NCYL-316] (10) settle the
DEFICIT candidate from stored data (see H7).  ⚠ [NCYL-315] bars renormalization /
self-similarity and [NCYL-296] bars routing through CF depth: neither bites -- nothing here
reads a CF digit or induces a return map.  ⚠ [OPS-242]'s gate does not bind: this is a
statement about a NAMED region, not a row-independent instrument, and the `7/15 eps1`
counterexample row has `min(P,Q-P) = 7`, i.e. `turns = 6`, so it lies outside the theorem's
hypothesis and must NOT be read as either support or refutation ([NCYL-321] s468 marker).

HYPOTHESES (pre-registered).  H0-H5 are GUARDS: given the proof they cannot fail, and a
failure means I have misread the source, not that the mathematics is wrong ([OPS-222]).
H6 and H7 are the two arms that can come back either way.

  H0  L1, on a SUPERSET of the path: for every kite `k = 1 .. Q-1` of every class row,
      `J_near ⊆ J_far`, and the shared endpoint is LO iff `not osec[k]`, HI iff `osec[k]`.
      ⚠ NEGATIVE CONTROL, and it is what makes H0 more than a tautology: the SAME test with
      the `osec` branch INVERTED must FAIL on essentially every row.  If both pass, the
      test is not reading `osec` at all.
      ⚠⚠ THE ONE FLOAT-GATED STEP IN THE WHOLE CHAIN, AUDITED HERE: `Chain.__init__` picks
      `far`/`near` by `a.f >= b.f`, a FLOAT comparison.  Score it against the exact `cmp`,
      and report the worst gap.  ⇒ the gap has a closed-form floor: the two `w`-values of a
      non-tie kite differ by at least `2`, and the flattest place on `sin` is the top, so
      the worst case is `1 - cos(π/Q)` -- `4.53e-03` at `Q = 33` against a double epsilon of
      `2.2e-16`.  ⇒ the label is safe until `Q ~ 1e8`, and L1's NESTING does not depend on
      it at all (a swapped label would give `delta < 0` and nest the other way, still
      nested); only the identification *near is the smaller* consumes it.

  H1  L2: the tie set is exactly `{0}` on every row, and kite `0` is never a path adjacency.
      ⚠ NEGATIVE CONTROL: the tie set must be NON-EMPTY (a vacuous lemma proves nothing) --
      i.e. kite `0` really is a tie, so L2 is a statement about WHERE the tie is.

  H2  L3 (a)(b)(c): `x_0 ∈ {0, Q}`; no interior `x_i` is a fold point; the value multiset is
      an AP of step `2` with distinct terms; the pairing is `t ↦ Q-1-t`.

  H3  L4: crossings of `{0, Q}` == interior extrema == `min(P, Q-P) - 1`, over a range far
      past the census box.  ⚠ NEGATIVE CONTROL: `P - 1` (the un-folded guess) must FAIL on
      every `P > Q/2` row -- otherwise the `min` is decoration.

  H4  L4's corollary: monotone iff `P ∈ {1, Q-1}`; one interior extremum iff `P ∈ {2, Q-2}`.

  H5  ⇒⇒ THE THEOREM.  On every row with no interior local MAXIMUM: every corner prong
      launches uphill, ends ESC, folds `0` times, and takes EXACTLY `j+1` or `m-j` steps.
      ⚠⚠ NEGATIVE CONTROL, AND IT IS THE ARM THAT SHOWS THE HYPOTHESIS DOES WORK: on TENT
      rows (`turns = 1` with the extremum a MAXIMUM) the zero-fold conclusion must FAIL.
      If tents also gave `0` folds, "no interior local max" would be decoration and the
      theorem would be proving something weaker than it claims.

  H6  ⇒⇒ THE TENT, WHICH IS THE PART THE PROOF DOES NOT REACH, AND I DO NOT PREDICT THE
      SHAPE OF THE ANSWER.  On `P ∈ {2, Q-2}` tent rows, odd `Q <= 101`: terminal kind, fold
      count, and the SEQUENCE OF FOLD EDGES per prong.  ⇒ the question is whether the fold
      edges are DISTINCT (which would bound the cascade by `m` and give termination) and
      whether they march outward from the peak.  ⚠ A capped prong here would REFUTE
      [NCYL-325]'s `0 CAP` on the Veech locus, so this is a live test of a published figure.

  H8  ⇒⇒ THE TENT INVARIANT, WHICH IS WHAT H6's TERMINATION ACTUALLY RESTS ON, AND IT IS
      SHARPER THAN H6 SUSPECTED.  H6 asks whether the fold edges are distinct; measure the
      quantity that makes them so.  ⇒ at each fold, record the SMALLER interface's profile
      value `min(a_k, a_{k+1})`.  PRE-REGISTERED PREDICTION, from the `Q = 27` cascade read
      by eye: the sequence is STRICTLY DECREASING.  ⚠ I do NOT predict the stronger form
      (that it steps down the value AP by exactly `2` each time), and if that holds it is
      the lemma to try to prove, not the weak one.
      ⚠⚠ **AND THE OBVIOUS STRUCTURAL EXPLANATION IS PRE-REGISTERED AS A CONTROL AND IS
      EXPECTED TO FAIL:** if the whole interface family were a TOTAL CHAIN under inclusion,
      `y`'s level would determine everything and the invariant would be immediate.  Score
      it.  If tents are not total chains, the invariant is real and its cause is subtler --
      which is the informative outcome, because it says the remaining lemma is not a
      one-liner about nesting.

  H7  ⇒⇒ (s468-c): THE DEFICIT CANDIDATE, SCORED AGAINST STORED DATA.  [NCYL-325] (4) names
      DEFICIT `d = Q-1-n(P/Q)` as the cheapest untested third cause and predicts the check
      needs `n(P/Q)`.  It does not: `data/deficit_certified.json` already carries `d` for
      `121` certified rows.  Cross-tabulate `d` against the cap set on the overlap, and
      cross-tabulate `P` PARITY against the cap set on all `460`.
      ⚠⚠ BOTH OUTCOMES ARE INFORMATIVE.  If `d = 0` tracks the cap-free rows, deficit is the
      better predictor and Veech-ness is a label on it.  If `d` is orthogonal to the cap,
      deficit is dead and the field narrows to ONE property under two names.
      ⚠ [OPS-071]: `d` is measured by a cylinder enumerator that shares NO code with the
      chain model, so this is an independent column, unlike `turns`, which is computed from
      the very profile that drives the walk -- `turns` agreeing with the cap boundary is
      therefore much weaker evidence than `d` disagreeing with it would be.

RUN: PYTHONPATH=.:engine:archive/scripts_2026-07:archive/scripts_2026-08:probes \
     .venv/bin/python3.13 probes/s469_profile_proof.py [all|h0..h8] [--qmax=33] [--big=201]
"""
from __future__ import annotations

import json
import math
import sys

import s450_chain_maps as cm          # noqa: F401  (Chain: the definitions L1 quotes)
import s451_quotient_path as qp
import s453_pole_module as pm
import s456_neckgb as nb

QMAX = 33            # the class-row census box (460 rows)
BIG = 201            # the pure-arithmetic box, for L3/L4 which need no chain build
TENTMAX = 101        # H6's box, matching [NCYL-324]'s own extreme-P range
OUT = "data/s469_profile_proof.json"


# ------------------------------------------------------------------ scaffolding

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
    return ch, lo, hi, Lo, Hi, S, pm.CoordOps(ch)


def profile(P, Q, eps):
    """L3: the folded rotation `(x_i, a_i)`, `i = 0..m`.  Pure arithmetic -- no ring."""
    twoQ, m = 2 * Q, (Q - 1) // 2
    c0 = P if eps == 0 else (P + Q) % twoQ
    x = [(c0 + 2 * P * (m - i)) % twoQ for i in range(m + 1)]
    return x, [min(v, twoQ - v) for v in x], m, twoQ


def turns_of(a):
    s = [1 if a[i + 1] > a[i] else -1 for i in range(len(a) - 1)]
    return sum(1 for i in range(len(s) - 1) if s[i] != s[i + 1])


def shape_of(a):
    """'mono' | 'valley' | 'tent' | 'multi' -- classified by the interior extrema."""
    t = turns_of(a)
    if t == 0:
        return "mono"
    if t > 1:
        return "multi"
    s = [1 if a[i + 1] > a[i] else -1 for i in range(len(a) - 1)]
    return "tent" if s[0] > 0 else "valley"


def has_interior_max(a):
    return any(a[i - 1] < a[i] > a[i + 1] for i in range(1, len(a) - 1))


# ------------------------------------------------------------------ H0

def h0(qmax=QMAX):
    """L1 on every kite `k = 1..Q-1`, plus the osec-inverted negative control."""
    r = {"rows": 0, "kites": 0, "nested": 0, "endpoint_ok": 0,
         "ctrl_inverted_ok": 0, "float_label_disagrees": 0,
         "worst_float_gap": None, "worst_float_gap_at": None,
         "predicted_floor": None, "bad": []}
    worst = float("inf")
    for P, Q, eps in _rows(qmax):
        c = _chain(P, Q, eps)
        if c is None:
            continue
        ch, lo, hi, _Lo, _Hi, _S, ops = c
        r["rows"] += 1
        for k in range(1, Q):
            if ch.tie[k]:
                continue
            f, n = ch.far[k], ch.near[k]
            r["kites"] += 1
            nested = ops.cmp(lo[f], lo[n]) <= 0 and ops.cmp(hi[n], hi[f]) <= 0
            share_lo, share_hi = lo[n].eq(lo[f]), hi[n].eq(hi[f])
            want_hi = bool(ch.osec[k])
            ok = (share_hi and not share_lo) if want_hi else (share_lo and not share_hi)
            ctrl = (share_lo and not share_hi) if want_hi else (share_hi and not share_lo)
            r["nested"] += int(nested)
            r["endpoint_ok"] += int(ok)
            r["ctrl_inverted_ok"] += int(ctrl)
            # the float-gated far/near label, audited against the exact comparison
            if ops.cmp(ch.ell[f], ch.ell[n]) <= 0:
                r["float_label_disagrees"] += 1
            gap = abs(ch.ell[f].f - ch.ell[n].f)
            if gap < worst:
                worst = gap
                r["worst_float_gap"] = gap
                r["worst_float_gap_at"] = {"P": P, "Q": Q, "eps": eps, "k": k}
                r["predicted_floor"] = 1.0 - math.cos(math.pi / Q)
            if not (nested and ok) and len(r["bad"]) < 20:
                r["bad"].append({"P": P, "Q": Q, "eps": eps, "k": k,
                                 "nested": nested, "share_lo": share_lo,
                                 "share_hi": share_hi, "osec": bool(ch.osec[k])})
    return r


# ------------------------------------------------------------------ H1

def h1(qmax=QMAX):
    """L2: the tie set is exactly {0}, and the path uses kites 1..m only."""
    r = {"rows": 0, "tieset_is_zero": 0, "tieset_empty": 0,
         "path_kites_exclude_0": 0, "bad": []}
    for P, Q, eps in _rows(qmax):
        c = _chain(P, Q, eps)
        if c is None:
            continue
        ch = c[0]
        r["rows"] += 1
        ts = sorted(k for k in range(Q) if ch.tie[k])
        r["tieset_is_zero"] += int(ts == [0])
        r["tieset_empty"] += int(ts == [])
        m = (Q - 1) // 2
        pk = {m - i for i in range(m)}          # kite of path adjacency (i, i+1)
        r["path_kites_exclude_0"] += int(0 not in pk and pk <= set(range(1, Q)))
        if ts != [0] and len(r["bad"]) < 20:
            r["bad"].append({"P": P, "Q": Q, "eps": eps, "tieset": ts})
    return r


# ------------------------------------------------------------------ H2

def h2(qmax=QMAX, big=BIG):
    """L3 (a)(b)(c), on the pure-arithmetic box AND against the real ring profile."""
    r = {"arith_rows": 0, "x0_is_foldpoint": 0, "no_interior_foldpoint": 0,
         "valueset_is_ap2": 0, "pairing_ok": 0, "ring_rows": 0, "ring_agrees": 0,
         "bad": []}
    for P, Q, eps in _rows(big):
        x, a, m, twoQ = profile(P, Q, eps)
        r["arith_rows"] += 1
        r["x0_is_foldpoint"] += int(x[0] in (0, Q))
        r["no_interior_foldpoint"] += int(all(v not in (0, Q) for v in x[1:]))
        s = sorted(a)
        r["valueset_is_ap2"] += int(len(set(a)) == len(a)
                                    and all(s[i + 1] - s[i] == 2 for i in range(len(s) - 1)))
        c0 = P if eps == 0 else (P + Q) % twoQ
        cc = [(c0 + 2 * P * t) % twoQ for t in range(Q)]
        r["pairing_ok"] += int(all(min(cc[t], twoQ - cc[t]) == min(cc[Q - 1 - t], twoQ - cc[Q - 1 - t])
                                   for t in range(Q)))
    # the arithmetic profile must reproduce the RING profile's ORDER (the thing L5 uses)
    for P, Q, eps in _rows(qmax):
        c = _chain(P, Q, eps)
        if c is None:
            continue
        _ch, _lo, _hi, Lo, Hi, _S, ops = c
        _x, a, m, _t = profile(P, Q, eps)
        r["ring_rows"] += 1
        L = [ops.sub(Hi[i], Lo[i]) for i in range(m + 1)]
        ok = all((ops.cmp(L[i + 1], L[i]) > 0) == (a[i + 1] > a[i]) for i in range(m))
        r["ring_agrees"] += int(ok)
        if not ok and len(r["bad"]) < 20:
            r["bad"].append({"P": P, "Q": Q, "eps": eps})
    return r


# ------------------------------------------------------------------ H3

def h3(big=BIG):
    """L4: crossings == interior extrema == min(P,Q-P)-1.  Control: P-1 must fail."""
    r = {"rows": 0, "turns_eq_minpq": 0, "cross_eq_turns": 0,
         "ctrl_pminus1_rows": 0, "ctrl_pminus1_ok": 0, "bad": []}
    for P, Q, eps in _rows(big):
        x, a, m, twoQ = profile(P, Q, eps)
        p = min(P, Q - P)
        t = turns_of(a)
        r["rows"] += 1
        r["turns_eq_minpq"] += int(t == p - 1)
        step = 2 * p
        cross = 0
        for i in range(m):
            d = (x[i + 1] - x[i]) % twoQ
            dirn = 1 if d == step else -1
            cross += sum(1 for s2 in range(1, step + 1)
                         if (x[i] + dirn * s2) % twoQ in (0, Q))
        r["cross_eq_turns"] += int(cross == t)
        if P > Q / 2:                                    # the control's target region
            r["ctrl_pminus1_rows"] += 1
            r["ctrl_pminus1_ok"] += int(t == P - 1)
        if t != p - 1 and len(r["bad"]) < 20:
            r["bad"].append({"P": P, "Q": Q, "eps": eps, "turns": t, "pred": p - 1})
    return r


# ------------------------------------------------------------------ H4

def h4(big=BIG):
    """monotone iff P in {1,Q-1}; one interior extremum iff P in {2,Q-2}."""
    r = {"rows": 0, "mono_iff_extreme": 0, "one_iff_second": 0,
         "n_mono": 0, "n_one": 0, "shapes": {}, "bad": []}
    for P, Q, eps in _rows(big):
        _x, a, _m, _t = profile(P, Q, eps)
        t = turns_of(a)
        sh = shape_of(a)
        r["rows"] += 1
        r["shapes"][sh] = r["shapes"].get(sh, 0) + 1
        r["n_mono"] += int(t == 0)
        r["n_one"] += int(t == 1)
        ok1 = (t == 0) == (P in (1, Q - 1))
        ok2 = (t == 1) == (P in (2, Q - 2))
        r["mono_iff_extreme"] += int(ok1)
        r["one_iff_second"] += int(ok2)
        if not (ok1 and ok2) and len(r["bad"]) < 20:
            r["bad"].append({"P": P, "Q": Q, "eps": eps, "turns": t})
    return r


# ------------------------------------------------------------------ H5

def _prong_walk(Lo, Hi, S, ops, st0, cap=400_000):
    """Return (kind, steps, folds, fold_edges).  Mirrors `pm.walk`'s step rule exactly."""
    n = len(Lo)
    i, d, y = st0
    folds, steps, fe = 0, 0, []
    for _ in range(cap):
        i2 = i + d
        steps += 1
        if i2 < 0 or i2 >= n:
            return "ESC", steps, folds, fe
        c0, c1 = ops.cmp(y, Lo[i2]), ops.cmp(y, Hi[i2])
        if c0 == 0 or c1 == 0:
            return "Z", steps, folds, fe
        if c0 > 0 and c1 < 0:
            i = i2
            continue
        k = min(i, i2)
        if ops.eq(y, ops.half(S[k])):
            return "R", steps, folds, fe
        y = ops.sub(S[k], y)
        d = -d
        folds += 1
        fe.append(k)
    return "CAP", steps, folds, fe


def h5(qmax=TENTMAX):
    """THE THEOREM, plus the tent negative control."""
    r = {"nomax_rows": 0, "nomax_prongs": 0, "nomax_esc": 0, "nomax_zerofold": 0,
         "nomax_steps_exact": 0, "by_shape": {},
         "ctrl_tent_rows": 0, "ctrl_tent_prongs": 0, "ctrl_tent_zerofold": 0,
         "bad": []}
    for P, Q, eps in _rows(qmax):
        _x, a, m, _t = profile(P, Q, eps)
        sh = shape_of(a)
        nomax = not has_interior_max(a)
        if not (nomax or sh == "tent"):
            continue                       # H5 is about the theorem's region + its control
        c = _chain(P, Q, eps)
        if c is None:
            continue
        _ch, _lo, _hi, Lo, Hi, S, ops = c
        L = [ops.sub(Hi[i], Lo[i]) for i in range(m + 1)]
        if nomax:
            r["nomax_rows"] += 1
        else:
            r["ctrl_tent_rows"] += 1
        b = r["by_shape"].setdefault(sh, {"rows": 0, "prongs": 0, "esc": 0,
                                          "zerofold": 0, "maxfold": 0})
        b["rows"] += 1
        for j in range(m):
            cs = nb.corner_start(Lo, Hi, S, j, ops)
            if cs is None:
                continue
            st0 = cs[0]
            kind, steps, folds, _fe = _prong_walk(Lo, Hi, S, ops, st0)
            b["prongs"] += 1
            b["esc"] += int(kind == "ESC")
            b["zerofold"] += int(folds == 0)
            b["maxfold"] = max(b["maxfold"], folds)
            if nomax:
                # launch must be uphill
                i0, d0, _y0 = st0
                i2 = i0 + d0
                up = (i2 < 0 or i2 > m) or ops.cmp(L[i2], L[i0]) > 0
                pred = (j + 1) if i0 == j else (m - j)     # keyed on the LARGER node
                r["nomax_prongs"] += 1
                r["nomax_esc"] += int(kind == "ESC")
                r["nomax_zerofold"] += int(folds == 0)
                r["nomax_steps_exact"] += int(steps == pred)
                if not (up and kind == "ESC" and folds == 0 and steps == pred) \
                        and len(r["bad"]) < 20:
                    r["bad"].append({"P": P, "Q": Q, "eps": eps, "j": j, "shape": sh,
                                     "up": up, "kind": kind, "folds": folds,
                                     "steps": steps, "pred": pred})
            else:
                r["ctrl_tent_prongs"] += 1
                r["ctrl_tent_zerofold"] += int(folds == 0)
    return r


# ------------------------------------------------------------------ H6

def h6(qmax=TENTMAX):
    """The TENT cascade: is it bounded, are the fold edges distinct, do they march out?"""
    r = {"rows": 0, "prongs": 0, "esc": 0, "cap": 0, "z": 0,
         "folds_le_m": 0, "edges_distinct": 0, "marches_outward": 0,
         "maxfold_eq_mminus1_rows": 0, "examples": [], "bad": []}
    for P, Q, eps in _rows(qmax):
        if min(P, Q - P) != 2:
            continue
        _x, a, m, _t = profile(P, Q, eps)
        if shape_of(a) != "tent":
            continue
        c = _chain(P, Q, eps)
        if c is None:
            continue
        _ch, _lo, _hi, Lo, Hi, S, ops = c
        peak = max(range(m + 1), key=lambda i: a[i])
        r["rows"] += 1
        rowmax = 0
        for j in range(m):
            cs = nb.corner_start(Lo, Hi, S, j, ops)
            if cs is None:
                continue
            kind, _st, folds, fe = _prong_walk(Lo, Hi, S, ops, cs[0])
            r["prongs"] += 1
            r["esc"] += int(kind == "ESC")
            r["cap"] += int(kind == "CAP")
            r["z"] += int(kind == "Z")
            r["folds_le_m"] += int(folds <= m)
            r["edges_distinct"] += int(len(set(fe)) == len(fe))
            # "marches outward": |edge - peak| strictly increases along the cascade
            dist = [abs(e - peak) for e in fe]
            r["marches_outward"] += int(all(dist[i + 1] > dist[i]
                                            for i in range(len(dist) - 1)))
            rowmax = max(rowmax, folds)
            if Q <= 27 and len(r["examples"]) < 12:
                r["examples"].append({"P": P, "Q": Q, "eps": eps, "j": j,
                                      "peak": peak, "kind": kind, "edges": fe})
            if kind != "ESC" and len(r["bad"]) < 20:
                r["bad"].append({"P": P, "Q": Q, "eps": eps, "j": j, "kind": kind})
        r["maxfold_eq_mminus1_rows"] += int(rowmax == m - 1)
    return r


# ------------------------------------------------------------------ H8

def _is_total_chain(Lo, Hi, ops):
    n = len(Lo)
    for x in range(n):
        for y in range(x + 1, n):
            if not ((ops.cmp(Lo[y], Lo[x]) <= 0 and ops.cmp(Hi[x], Hi[y]) <= 0)
                    or (ops.cmp(Lo[x], Lo[y]) <= 0 and ops.cmp(Hi[y], Hi[x]) <= 0)):
                return False
    return True


def h8(qmax=TENTMAX):
    """The TENT fold-level invariant, plus the total-chain control (expected to fail)."""
    r = {"tent_rows": 0, "tent_prongs": 0, "strictly_decreasing": 0, "step_is_minus2": 0,
         "ctrl_chain_rows_by_shape": {}, "bad": []}
    for P, Q, eps in _rows(qmax):
        _x, a, m, _t = profile(P, Q, eps)
        sh = shape_of(a)
        c = _chain(P, Q, eps)
        if c is None:
            continue
        _ch, _lo, _hi, Lo, Hi, S, ops = c
        if Q <= 51:                                   # the control, all shapes, cheap box
            b = r["ctrl_chain_rows_by_shape"].setdefault(sh, {"rows": 0, "chains": 0})
            b["rows"] += 1
            b["chains"] += int(_is_total_chain(Lo, Hi, ops))
        if sh != "tent":
            continue
        r["tent_rows"] += 1
        n = m + 1
        for j in range(m):
            cs = nb.corner_start(Lo, Hi, S, j, ops)
            if cs is None:
                continue
            i, d, y = cs[0]
            seq = []
            for _ in range(400_000):
                i2 = i + d
                if i2 < 0 or i2 >= n:
                    break
                c0, c1 = ops.cmp(y, Lo[i2]), ops.cmp(y, Hi[i2])
                if c0 == 0 or c1 == 0:
                    break
                if c0 > 0 and c1 < 0:
                    i = i2
                    continue
                k = min(i, i2)
                if ops.eq(y, ops.half(S[k])):
                    break
                seq.append(min(a[k], a[k + 1]))
                y = ops.sub(S[k], y)
                d = -d
            r["tent_prongs"] += 1
            dec = all(seq[t + 1] < seq[t] for t in range(len(seq) - 1))
            two = all(seq[t] - seq[t + 1] == 2 for t in range(len(seq) - 1))
            r["strictly_decreasing"] += int(dec)
            r["step_is_minus2"] += int(two)
            if not dec and len(r["bad"]) < 20:
                r["bad"].append({"P": P, "Q": Q, "eps": eps, "j": j, "seq": seq[:12]})
    return r


# ------------------------------------------------------------------ H7

def h7():
    """(s468-c): score DEFICIT and P-PARITY against the cap set.  Stored data only."""
    with open("data/deficit_certified.json") as fh:
        dc = [x for x in json.load(fh) if x.get("certified")]
    with open("data/s468_veech_split.json") as fh:
        vs = json.load(fh)
    dmap = {(x["P"], x["Q"]): x["d"] for x in dc}
    cen = {}
    for x in vs:
        cen.setdefault((x["P"], x["Q"]), []).append(x)
    ov = [k for k in dmap if k in cen]
    rows = []
    for (P, Q) in sorted(ov, key=lambda k: (k[1], k[0])):
        rr = cen[(P, Q)]
        rows.append({"P": P, "Q": Q, "minpq": min(P, Q - P), "d": dmap[(P, Q)],
                     "veech": rr[0]["veech"], "Podd": P % 2,
                     "cap": sum(z["cap"] for z in rr),
                     "maxfold": max(z["maxfold"] for z in rr)})
    par = {}
    for x in vs:
        b = par.setdefault(x["P"] % 2, {"rows": 0, "caprows": 0, "capprongs": 0})
        b["rows"] += 1
        b["caprows"] += int(x["cap"] > 0)
        b["capprongs"] += x["cap"]
    return {
        "certified_rows": len(dmap),
        "odd_Q_certified": sum(1 for (_P, Q) in dmap if Q % 2),
        "overlap_rows": len(rows),
        "overlap_all_odd_P": all(x["Podd"] == 1 for x in rows),
        "overlap_d_zero": sum(1 for x in rows if x["d"] == 0),
        "overlap_d_range": [min(x["d"] for x in rows), max(x["d"] for x in rows)],
        "veech_d_range": [min((x["d"] for x in rows if x["veech"]), default=None),
                          max((x["d"] for x in rows if x["veech"]), default=None)],
        "capped_d_values": sorted(x["d"] for x in rows if x["cap"] > 0),
        "by_P_parity": par,
        "rows": rows,
    }


# ------------------------------------------------------------------ main

ARMS = {"h0": h0, "h1": h1, "h2": h2, "h3": h3, "h4": h4, "h5": h5, "h6": h6,
        "h7": h7, "h8": h8}


def main(argv):
    which = [a for a in argv[1:] if not a.startswith("--")] or ["all"]
    kw = {}
    for a in argv[1:]:
        if a.startswith("--qmax="):
            kw["qmax"] = int(a.split("=")[1])
        if a.startswith("--big="):
            kw["big"] = int(a.split("=")[1])
    names = list(ARMS) if which == ["all"] else which
    out = {}
    for nm in names:
        fn = ARMS[nm]
        ok = {k: v for k, v in kw.items()
              if k in fn.__code__.co_varnames[:fn.__code__.co_argcount]}
        out[nm] = fn(**ok)
        r = dict(out[nm])
        r.pop("rows", None) if isinstance(r.get("rows"), list) else None
        print(f"== {nm} ==")
        print(json.dumps({k: v for k, v in r.items()
                          if k not in ("rows",) or not isinstance(v, list)},
                         default=str)[:1600])
    with open(OUT, "w") as fh:
        json.dump(out, fh, indent=1, default=str)
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main(sys.argv)
