#!/usr/bin/env python3
"""s468_ordering.py -- PRE-REGISTERED.  Queue item (s464-b): THE ORDERING TESTS.

[NCYL-319] (3) closed the arithmetic side of (G0b-rel) NEGATIVELY for the discriminating
purpose: a `+-1` fold word reaching a genuine relation exists at `Q = 15, 21, 33, 39`, and
three of those four are witness-free.  Its verdict: *what keeps the others witness-free is
the INTERVAL-ORDERING TESTS, not the arithmetic and not the word combinatorics*, so the
next instrument must be about the tests.  This file is that instrument.

WHAT THIS IS.  NOT a proof of (G0b) and NOT a claim about `Q = 15`.  It is a structural
reading of `s454_self_hit.step`'s two comparisons, and the claim it makes is that they are
NOT a free binary choice: on a step into the LARGER of two nested interfaces the outcome is
FORCED, so the test does real work on DOWNHILL steps only.

  ⇒⇒ THE OBSERVATION.  `step` at state `(i, d, y)` compares `y` against `Lo[i2]`, `Hi[i2]`,
  `i2 = i + d`.  The path interfaces are NESTED and adjacent ones share EXACTLY ONE
  endpoint (H0, and `s456_neckgb.corner_start` already relies on the shared-endpoint half
  via its `lo_shared` branch).  Hence with `J_i := [Lo[i], Hi[i]]`, adjacency gives
  `J_i subset J_{i2}` or `J_{i2} subset J_i`, decided by length.  So:

      UPHILL   `|J_{i2}| > |J_i|`, i.e. `J_i subset J_{i2}`:  every `y in J_i` lies in
               `J_{i2}`, so the step is CROSS -- FORCED -- unless `y` is the SHARED
               endpoint, which is an endpoint of BOTH and therefore reads `Z`.
               ⇒ on an uphill step the ordering test decides nothing; `y` is UNCHANGED by a
               cross, so an uphill run carries one fixed `y` and terminates only by
               running off the path (ESC) or by `y` coinciding with a shared endpoint (Z).

      DOWNHILL `|J_{i2}| < |J_i|`, i.e. `J_{i2} subset J_i`:  `y` may be interior to
               `J_{i2}` (CROSS), outside it (FOLD), or at either of its endpoints (Z).
               ⇒ this is the ONLY place the ordering test is two-sided, and the only place
               a FOLD can be born.

  ⇒ (G0b) is therefore a question about the DESCENT structure of the interface-length
  profile `|J_i|`, plus the coincidences `y = ` an interface endpoint.  ⚠ The profile is
  `|J_i| = ell` at the node's own `c`, and [NCYL-294]'s CF1 makes `c_i` an ARITHMETIC
  PROGRESSION mod `2Q` -- so the descent set is the descent set of a rotation, and rows
  with a MONOTONE profile have no downhill step in one of the two directions at all.

PRIOR ART: `rulings.py --grep` on 'ordering test' -> [NCYL-319] only (which named the
residual and did not attack it); 'feasible', 'polytope', 'inequality', 'y0 interval',
'piecewise' -> nothing on this; 'realisable' -> [NCYL-300]/[NCYL-319]/[NCYL-322], none
about the tests; 'corner prong' -> [NCYL-307]/[NCYL-308]/[NCYL-309]/[NCYL-312]/[NCYL-313]/
[NCYL-319]; 'periodic orbit', 'closes up', 'returns to its start', 'prong fate' -> nothing.
⇒ ID-grep per [OPS-246] on [NCYL-308]/[NCYL-319]/[NCYL-322]/[NCYL-292]/[NCYL-297]:
⚠⚠ THE ONE REAL HIT, AND IT IS THE THING THIS GENERALISES.  [NCYL-292]'s paired-only
one-kite partial IS this mechanism at a single kite -- *`ell[m] = 1` is the strict max, so
kite `m` passes everything through* -- and [NCYL-297] measures its consequence (kite `m` a
fold DESERT in paired, a fold HOTSPOT in lone).  Both are about INTERIOR POLE prongs at the
ONE extremal kite.  ⇒ what is unbuilt, and is this file, is the same nesting argument run
(a) over the WHOLE profile rather than its maximum and (b) on CORNER prongs, which are the
objects (G0b) quantifies over ([NCYL-307] H5) and which [NCYL-301]/[NCYL-298] explicitly do
NOT cover ([NCYL-308] (1): a corner is the shared endpoint, so no arrival node is forced).
⚠ [NCYL-315] bars renormalization / self-similarity: does not bite -- nothing below reads a
CF digit or induces a return map.  ⚠ [NCYL-296] bars routing through CF depth: same.
⚠ user-directed s460 bars EXTENDING THE RANGE as a way to characterise: does not bite --
the box is the existing census box and no bigger box is proposed as an answer.
⚠ [OPS-242]'s GATE: (G0b) is FALSE without its hypothesis (the `7/15 eps1` witness is exact
and terminating), so a ROW-INDEPENDENT instrument is dead a priori.  This one is
row-dependent throughout -- it consumes the row's own `Lo`/`Hi`/`S` and its own profile --
and H5 RUNS IT ON THE COUNTEREXAMPLE ROW, which is the gate's operative demand.

HYPOTHESES (pre-registered).

  H0  ⇒⇒ THE STRUCTURAL INPUT, and everything else is void if it fails.  For every adjacent
      pair on every class row, odd `Q <= QMAX`: the two interfaces share EXACTLY ONE
      endpoint, and one CONTAINS the other.  PREDICTION: 100% on both counts.
      ⚠ It can fail on either count and the failure modes differ: no shared endpoint would
      break `corner_start` too; sharing one endpoint WITHOUT nesting is impossible in 1-D
      only if both are non-degenerate, and the LONE class carries a DEGENERATE interface
      (`ell = 0`, [NCYL-319] H0), so the degenerate rows are the ones to watch.

  H1  ⇒⇒ THE FORCED-CROSS LAW, and ONLY ITS CONVERSE IS EVIDENCE.  Two parts.
      (a) On an UPHILL step the outcome is CROSS, or `Z` at the shared endpoint -- never
          FOLD.  ⚠⚠ GIVEN H0 THIS IS A DERIVATION AND CANNOT FAIL ([OPS-222]); it is scored
          only so a reader can see the algebra ran, and it is NOT evidence.
      (b) ⇒⇒ THE CONVERSE IS THE ARM: on a DOWNHILL step, is the test genuinely two-sided?
          Score CROSS-vs-FOLD over all downhill steps of all real corner-prong orbits.
          ⚠ If downhill ALWAYS folds, the walk is combinatorially determined by the profile
          alone, the ordering test carries no information at all, and (G0b) reduces to a
          statement about the descent word -- a MUCH stronger and quite different result.
          ⚠ If downhill folds only sometimes, the test is doing real work and (b) prices how
          much.  BOTH outcomes are informative and I do not predict which.

  H2  ⇒⇒ THE ARM THAT CAN CLOSE ROWS: MONOTONE-PROFILE ROWS ESCAPE, CAP-FREE.  On a row
      whose profile `|J_i|` is monotone along the path, one of the two directions has NO
      downhill step, so a prong heading that way crosses to the end of the path and ESCapes,
      and a prong heading the other way ESCapes as soon as one fold reverses it.  ⇒ on such
      a row every corner prong terminates in a BOUNDED number of steps with no cap, so
      (G0b) holds there by a terminating argument rather than by a capped census.
      PRE-REGISTERED PREDICTION: every corner prong of a monotone row ends ESC (or `Z` by an
      endpoint coincidence), never CAP, in `<= 2 * len(path)` steps.
      ⚠ It CAN fail: a fold changes `y`, and nothing above forbids a SECOND fold before the
      reversed prong reaches the end -- if downhill steps exist in the reversed direction
      too, the bound is wrong and the row is not closed.  That is exactly why H2 is scored
      and not asserted.
      ⚠⚠ AND THE COVERAGE IS THE POINT, NOT THE PASS RATE: report how many rows are
      monotone.  A law that closes 3 rows out of 460 is a curiosity ([OPS-041] -- quote the
      coverage or do not quote the law).

  H3  ⇒⇒ THE PROFILE'S DESCENT STRUCTURE, MEASURED, since H2's reach is set by it.  Per
      class row: `ndesc` = number of `i` with `|J_{i+1}| < |J_i|`, and the same for the
      reversed direction.  ⇒ `ndesc = 0` in one direction IS H2's monotone class.  Also
      report the distribution against `P`, since CF1 makes `c_i` an arithmetic progression
      and the descent count of a rotation is a three-distance quantity.
      ⚠ This is a MEASUREMENT of an arithmetic fact, not an arm; it cannot fail, it
      calibrates H2's coverage.

  H4  ⇒⇒ THE FATE CENSUS UNDER THE DESCENT LENS -- how far the mechanism reaches beyond
      monotone rows.  For every corner prong on every row: run the REAL walk (`sh.orbit`,
      the audited stepper) to `CAP4` and record terminal kind, step count, number of folds,
      and the maximum number of folds between consecutive ESC-ward runs.  ⇒ the quantity
      that matters is FOLD COUNT: a prong that folds `k` times and then escapes is decided;
      a prong still folding at the cap is not.
      ⚠ NOT EVIDENCE FOR (G0b): a capped null is not a verdict ([NCYL-262]/[OPS-041]), and
      the ESC counts here REPRODUCE [NCYL-307] H5 (`4802` prongs -> `4684 ESC / 114 CAP /
      4 Z` at `Q <= 33`) rather than adding to them.  ⇒ AGREEMENT WITH H5's PUBLISHED SPLIT
      IS A BUG-CHECK ON THIS FILE, and a DISAGREEMENT is a bug in this file, not a finding.

  H5  ⇒⇒ THE HIGHEST-VALUE OUTPUT, AND IT IS ABOUT THE COUNTEREXAMPLE: WHICH KIND OF `Z` IS
      THE `7/15 eps1` WITNESS?  The mechanism admits exactly two kinds --
        (U) an UPHILL `Z`: `y` equals the SHARED endpoint of the pair being entered;
        (D) a DOWNHILL `Z`: `y` equals an endpoint of the strictly SMALLER interface.
      Classify all `4` census `Z` prongs ([NCYL-307] H5, `7/15 eps1` and `8/15 eps0`).
      ⚠⚠ THIS IS THE ARM WITH A REAL STAKE, because the two kinds carry DIFFERENT residuals:
      (U) makes (G0b) a statement about `y_0` avoiding a set of `<= len(path)` shared
      endpoints along ONE uphill run -- a finite, explicit, cap-free avoidance condition;
      (D) makes it a statement about the fold orbit and keeps the dynamics in play.
      ⇒ I do NOT predict which, and the two outcomes name DIFFERENT next instruments.

  H7  ⇒⇒⇒ THE ARM THIS FILE EXISTS FOR, AND IT IS ONLY AVAILABLE BECAUSE OF H1(a): THE
      PROFILE-RESTRICTED `+-1` ENUMERATION.  H1(a) says a FOLD can be born only on a
      DOWNHILL step -- and the uphill/downhill label is a function of `(i, d)` ALONE, read
      off the interface lengths, with NO reference to `y`.  ⇒ it can therefore be imposed
      on [NCYL-319]'s relaxation WITHOUT putting the `y`-comparison back, i.e. without
      un-relaxing the ordering test.  Soundness is what a necessary condition needs and it
      holds by H1(a): every real orbit satisfies `no uphill fold` (`0` counterexamples).
      ⇒ re-run `s464_realisable.relaxed_reach` verbatim EXCEPT that FOLD is offered only
      when `|J_{i2}| < |J_i|`, and re-classify every admitted arrival `form`/`triv`/`rel`.
      ⇒⇒ THIS SITS STRICTLY BETWEEN TWO MEASURED NEIGHBOURS, exactly as s464's did: above
      it, s464's unrestricted enumeration (`1.31%` of pairs at depth `20`, `rel` at
      `Q = 15, 21, 33, 39`); below it, the real walk (`4 Z`, `Q = 15` only).
      ⚠⚠ BOTH OUTCOMES ARE INFORMATIVE AND I DO NOT PREDICT WHICH:
        `rel` SURVIVES at `Q = 21, 33, 39` ⇒ the descent structure is NOT the obstruction
             either; [NCYL-319]'s *dynamical* verdict stands unrefined and the residual is
             the genuinely `y`-DEPENDENT half of the test.
        `rel` DIES at `Q = 21, 33, 39` and survives at `Q = 15` ⇒ the obstruction is the
             DESCENT WORD -- a `y`-independent, purely combinatorial criterion -- and
             [NCYL-319] (3)'s *the obstruction is DYNAMICAL* is REFINED: it is dynamical
             only in the weak sense of depending on the profile, not on the orbit.
      ⚠ POSITIVE CONTROL, PREDICTED, HENCE A BUG-CHECK AND NEVER EVIDENCE ([OPS-242]): the
      real `7/15 eps1` witness MUST still be admitted -- it is a real orbit, so H1(a)
      guarantees it violates no uphill-fold constraint.  A miss is an implementation bug.
      ⚠ Depth- and state-capped exactly as s464 was, and `closed` is reported per
      [OPS-245]: a null with `closed = False` is depth-limited and says nothing.

  H6  ⇒⇒ INDEPENDENT VERIFICATION, because H1/H2 classify steps by a PROFILE comparison
      while the truth is `s454_self_hit.step`'s own two comparisons, and a shared bug
      between my classifier and the walker would be invisible ([WFLOOR-124]/[OPS-244]).
      Re-derive every step's kind from `sh.step` directly and require agreement with the
      uphill/downhill prediction on EVERY step of EVERY orbit traced.  A single mismatch
      fails the file.

⚠⚠ WHAT NONE OF THIS IS.  It does not close (G0b), it moves no `claims.md` tier and no
γ=1 DAG node, and it lives inside [NCYL-294]'s model, whose global-coordinate closure is
VERIFIED (`420/420`), not proved.  `CP` remains a hypothesis of [NECK-P0]/[NECK-COVER].

RUN: PYTHONPATH=.:engine:archive/scripts_2026-07:archive/scripts_2026-08:probes \
     .venv/bin/python3.13 probes/s468_ordering.py [all|h0|h1|h2|h3|h4|h5|h6] [--qmax=33]
"""
from __future__ import annotations

import json
import math
import sys
import time

import s450_chain_maps as cm
import s451_quotient_path as qp
import s453_pole_module as pm
import s454_self_hit as sh
import s456_neckgb as nb

QMAX = 33
CAP4 = 200_000        # H4 fate census; [NCYL-307] H5's own run capped comparably
OUT = "data/s468_ordering.json"


# ------------------------------------------------------------------ shared scaffolding

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


def _len(Lo, Hi, i, ops):
    return ops.sub(Hi[i], Lo[i])


def adjacency(Lo, Hi, ops, i, j):
    """Structure of the adjacent pair `(i, j)`.  Returns a dict, or None if not adjacent.

    `shared`  : 'lo' | 'hi' | None      -- which endpoint the two interfaces share
    `nest`    : 'i_in_j' | 'j_in_i' | 'equal' | None
    `cmp_len` : sign of `|J_j| - |J_i|` -- the UPHILL/DOWNHILL sign, computed from the
                LENGTHS and never from the float profile.
    """
    slo, shi = ops.eq(Lo[i], Lo[j]), ops.eq(Hi[i], Hi[j])
    shared = "lo" if (slo and not shi) else ("hi" if (shi and not slo) else None)
    a, b = ops.cmp(Lo[i], Lo[j]), ops.cmp(Hi[i], Hi[j])
    if a <= 0 and b >= 0:
        nest = "equal" if (a == 0 and b == 0) else "j_in_i"
    elif a >= 0 and b <= 0:
        nest = "i_in_j"
    else:
        nest = None
    return {"shared": shared, "nest": nest,
            "cmp_len": ops.cmp(_len(Lo, Hi, j, ops), _len(Lo, Hi, i, ops))}


# ------------------------------------------------------------------ H0

def h0(qmax=QMAX):
    """THE STRUCTURAL INPUT: adjacent interfaces share exactly one endpoint and nest."""
    res = {"rows": 0, "pairs": 0, "shared_one": 0, "nested": 0,
           "degenerate_rows": 0, "bad": []}
    for P, Q, eps in _rows(qmax):
        c = _chain(P, Q, eps)
        if c is None:
            continue
        _ch, _lo, _hi, Lo, Hi, S, ops = c
        res["rows"] += 1
        if any(ops.eq(Lo[i], Hi[i]) for i in range(len(Lo))):
            res["degenerate_rows"] += 1
        for i in range(len(Lo) - 1):
            res["pairs"] += 1
            ad = adjacency(Lo, Hi, ops, i, i + 1)
            ok_s = ad["shared"] is not None
            ok_n = ad["nest"] is not None
            res["shared_one"] += int(ok_s)
            res["nested"] += int(ok_n)
            if not (ok_s and ok_n) and len(res["bad"]) < 20:
                res["bad"].append({"P": P, "Q": Q, "eps": eps, "i": i, **ad})
    return res


# ------------------------------------------------------------------ H1 / H6

def _walk_steps(Lo, Hi, S, ops, st0, cap):
    """Real walk from `st0`.  Yields per-step records with BOTH classifications.

    `pred` is this file's profile-based prediction (uphill => forced), `kind` is the
    audited `sh.step` verdict.  H6 requires them consistent on every step.
    """
    st, out = st0, []
    for _ in range(cap):
        i, d, y = st
        i2 = i + d
        if i2 < 0 or i2 >= len(Lo):
            out.append({"kind": "ESC", "dir": None})
            return out, "ESC", st
        ad = adjacency(Lo, Hi, ops, i, i2)
        updown = "up" if ad["cmp_len"] > 0 else ("down" if ad["cmp_len"] < 0 else "flat")
        nxt, kind, _e = sh.step(Lo, Hi, S, st, ops)
        rec = {"kind": kind, "dir": updown, "shared": ad["shared"], "nest": ad["nest"]}
        if kind == "Z":
            # (U) vs (D) of H5: is the hit AT the shared endpoint, or at the free endpoint
            # of the strictly smaller interface?
            at_shared = ((ad["shared"] == "lo" and ops.eq(y, Lo[i2]))
                         or (ad["shared"] == "hi" and ops.eq(y, Hi[i2])))
            rec["zkind"] = "U" if at_shared else "D"
        out.append(rec)
        if nxt is None:
            return out, kind, st
        st = nxt
    return out, "CAP", st


def h1(qmax=QMAX, cap=CAP4):
    """(a) uphill never folds -- a DERIVATION CHECK.  (b) is downhill two-sided -- the ARM."""
    res = {"rows": 0, "prongs": 0, "steps": 0,
           "up_cross": 0, "up_fold": 0, "up_z": 0, "up_other": 0,
           "down_cross": 0, "down_fold": 0, "down_z": 0, "down_other": 0,
           "flat": 0, "h6_mismatch": [], "capped": 0}
    for P, Q, eps in _rows(qmax):
        c = _chain(P, Q, eps)
        if c is None:
            continue
        _ch, _lo, _hi, Lo, Hi, S, ops = c
        res["rows"] += 1
        for j in range(len(S)):
            cs = nb.corner_start(Lo, Hi, S, j, ops)
            if cs is None:
                continue
            res["prongs"] += 1
            steps, end, _st = _walk_steps(Lo, Hi, S, ops, cs[0], cap)
            res["capped"] += int(end == "CAP")
            for r in steps:
                if r["dir"] is None:
                    continue
                res["steps"] += 1
                if r["dir"] == "flat":
                    res["flat"] += 1
                    continue
                pre = "up_" if r["dir"] == "up" else "down_"
                k = r["kind"]
                key = pre + (k if k in ("cross", "fold") else
                             ("z" if k == "Z" else "other"))
                res[key] = res.get(key, 0) + 1
                # H6: an uphill FOLD contradicts H0+nesting outright.
                if r["dir"] == "up" and k == "fold" and len(res["h6_mismatch"]) < 20:
                    res["h6_mismatch"].append({"P": P, "Q": Q, "eps": eps, "j": j, **r})
    tot_down = res["down_cross"] + res["down_fold"]
    res["down_fold_rate"] = round(res["down_fold"] / tot_down, 6) if tot_down else None
    return res


# ------------------------------------------------------------------ H3 (profile descents)

def profile_descents(Lo, Hi, ops):
    """`(ndesc_fwd, ndesc_bwd, n)` for the interface-length profile along the path."""
    n = len(Lo)
    L = [_len(Lo, Hi, i, ops) for i in range(n)]
    fwd = sum(1 for i in range(n - 1) if ops.cmp(L[i + 1], L[i]) < 0)
    bwd = sum(1 for i in range(n - 1) if ops.cmp(L[i], L[i + 1]) < 0)
    return fwd, bwd, n


def h3(qmax=QMAX):
    """THE DESCENT STRUCTURE -- a measurement, and the calibration of H2's coverage."""
    res = {"rows": 0, "monotone": 0, "mono_rows": [], "by_Q": {}, "hist": {}}
    for P, Q, eps in _rows(qmax):
        c = _chain(P, Q, eps)
        if c is None:
            continue
        _ch, _lo, _hi, Lo, Hi, S, ops = c
        fwd, bwd, n = profile_descents(Lo, Hi, ops)
        res["rows"] += 1
        mono = (fwd == 0 or bwd == 0)
        res["monotone"] += int(mono)
        if mono and len(res["mono_rows"]) < 60:
            res["mono_rows"].append({"P": P, "Q": Q, "eps": eps, "fwd": fwd,
                                     "bwd": bwd, "n": n})
        k = str(Q)
        res["by_Q"].setdefault(k, {"rows": 0, "mono": 0})
        res["by_Q"][k]["rows"] += 1
        res["by_Q"][k]["mono"] += int(mono)
        m = str(min(fwd, bwd))
        res["hist"][m] = res["hist"].get(m, 0) + 1
    res["mono_rate"] = round(res["monotone"] / res["rows"], 6) if res["rows"] else None
    return res


# ------------------------------------------------------------------ H2

def h2(qmax=QMAX, cap=CAP4):
    """MONOTONE ROWS ESCAPE, CAP-FREE -- and the bound `<= 2*len(path)` is scored."""
    res = {"mono_rows": 0, "prongs": 0, "esc": 0, "z": 0, "pole": 0, "capped": 0,
           "within_bound": 0, "max_steps": 0, "max_folds": 0, "over_bound": []}
    for P, Q, eps in _rows(qmax):
        c = _chain(P, Q, eps)
        if c is None:
            continue
        _ch, _lo, _hi, Lo, Hi, S, ops = c
        fwd, bwd, n = profile_descents(Lo, Hi, ops)
        if not (fwd == 0 or bwd == 0):
            continue
        res["mono_rows"] += 1
        bound = 2 * n
        for j in range(len(S)):
            cs = nb.corner_start(Lo, Hi, S, j, ops)
            if cs is None:
                continue
            res["prongs"] += 1
            steps, end, _st = _walk_steps(Lo, Hi, S, ops, cs[0], cap)
            nst = len(steps)
            nfold = sum(1 for r in steps if r["kind"] == "fold")
            res["max_steps"] = max(res["max_steps"], nst)
            res["max_folds"] = max(res["max_folds"], nfold)
            res[{"ESC": "esc", "Z": "z", "pole": "pole", "CAP": "capped"}
                .get(end, "pole")] += 1
            if nst <= bound:
                res["within_bound"] += 1
            elif len(res["over_bound"]) < 20:
                res["over_bound"].append({"P": P, "Q": Q, "eps": eps, "j": j,
                                          "steps": nst, "bound": bound, "end": end,
                                          "folds": nfold})
    return res


# ------------------------------------------------------------------ H4

def h4(qmax=QMAX, cap=CAP4):
    """THE FATE CENSUS under the descent lens.  Reproduces [NCYL-307] H5 as a bug-check."""
    res = {"rows": 0, "prongs": 0, "ESC": 0, "Z": 0, "pole": 0, "CAP": 0,
           "fold_hist": {}, "z_rows": [], "max_folds_esc": 0}
    for P, Q, eps in _rows(qmax):
        c = _chain(P, Q, eps)
        if c is None:
            continue
        _ch, _lo, _hi, Lo, Hi, S, ops = c
        res["rows"] += 1
        for j in range(len(S)):
            cs = nb.corner_start(Lo, Hi, S, j, ops)
            if cs is None:
                continue
            res["prongs"] += 1
            steps, end, _st = _walk_steps(Lo, Hi, S, ops, cs[0], cap)
            res[end] = res.get(end, 0) + 1
            nfold = sum(1 for r in steps if r["kind"] == "fold")
            b = ("0" if nfold == 0 else "1" if nfold == 1 else "2" if nfold == 2
                 else "3-9" if nfold < 10 else "10-99" if nfold < 100 else "100+")
            res["fold_hist"][b] = res["fold_hist"].get(b, 0) + 1
            if end == "ESC":
                res["max_folds_esc"] = max(res["max_folds_esc"], nfold)
            if end == "Z":
                res["z_rows"].append({"P": P, "Q": Q, "eps": eps, "j": j,
                                      "steps": len(steps), "folds": nfold,
                                      "zkind": steps[-1].get("zkind")})
    return res


# ------------------------------------------------------------------ H5

def h5(qmax=QMAX, cap=CAP4):
    """WHICH KIND OF `Z` IS THE WITNESS -- (U) shared endpoint or (D) smaller interface."""
    res = {"hits": [], "U": 0, "D": 0}
    for P, Q, eps in _rows(qmax):
        c = _chain(P, Q, eps)
        if c is None:
            continue
        _ch, _lo, _hi, Lo, Hi, S, ops = c
        for j in range(len(S)):
            cs = nb.corner_start(Lo, Hi, S, j, ops)
            if cs is None:
                continue
            steps, end, _st = _walk_steps(Lo, Hi, S, ops, cs[0], cap)
            if end != "Z":
                continue
            last = steps[-1]
            res[last.get("zkind", "D")] += 1
            res["hits"].append({
                "P": P, "Q": Q, "eps": eps, "j": j, "steps": len(steps),
                "zkind": last.get("zkind"), "dir": last["dir"],
                "shared": last["shared"], "nest": last["nest"],
                "folds": sum(1 for r in steps if r["kind"] == "fold"),
                "word": "".join({"cross": "C", "fold": "F", "Z": "Z"}.get(r["kind"], "?")
                                for r in steps),
                "dirs": "".join({"up": "u", "down": "d", "flat": "-"}.get(r["dir"], "?")
                                for r in steps if r["dir"]),
            })
    return res


# ------------------------------------------------------------------ H6

def h6(qmax=QMAX, cap=CAP4):
    """INDEPENDENT VERIFICATION: profile prediction vs the audited `sh.step`, every step."""
    res = {"steps": 0, "agree": 0, "mismatch": []}
    for P, Q, eps in _rows(qmax):
        c = _chain(P, Q, eps)
        if c is None:
            continue
        _ch, _lo, _hi, Lo, Hi, S, ops = c
        for j in range(len(S)):
            cs = nb.corner_start(Lo, Hi, S, j, ops)
            if cs is None:
                continue
            steps, _end, _st = _walk_steps(Lo, Hi, S, ops, cs[0], cap)
            for r in steps:
                if r["dir"] in (None, "flat"):
                    continue
                res["steps"] += 1
                # the prediction: uphill => NOT fold.  downhill => anything.
                ok = not (r["dir"] == "up" and r["kind"] == "fold")
                res["agree"] += int(ok)
                if not ok and len(res["mismatch"]) < 20:
                    res["mismatch"].append({"P": P, "Q": Q, "eps": eps, "j": j, **r})
    return res


# ------------------------------------------------------------------ H7

def restricted_reach(Lo, Hi, S, ops, st0, sym0, sLo, sHi, sS,
                     depth, maxstates, widest, restrict=True):
    """`s464_realisable.relaxed_reach` with ONE change: FOLD only on a DOWNHILL step.

    Verbatim otherwise -- same `Z` test, same pole handling, same dedup key, same `closed`
    semantics.  `restrict=False` reproduces s464 exactly and is the calibration arm.
    """
    key0 = (st0[0], st0[1], st0[2].v.tobytes())
    seen = {key0}
    frontier = [(st0, tuple(sym0))]
    zhits, truncated, reached, closed = [], False, 0, False
    for dstep in range(depth):
        nxt = []
        for (i, d, y), ysym in frontier:
            i2 = i + d
            if i2 < 0 or i2 >= len(Lo):
                continue
            c0, c1 = ops.cmp(y, Lo[i2]), ops.cmp(y, Hi[i2])
            if c0 == 0 or c1 == 0:
                tgt = sLo[i2] if c0 == 0 else sHi[i2]
                zhits.append({"steps": dstep + 1, "node": i2,
                              "at": "lo" if c0 == 0 else "hi",
                              "a": [int(x - t) for x, t in zip(ysym, tgt)]})
                continue
            cands = [((i2, d, y), ysym)]                       # CROSS, always offered
            e = min(i, i2)
            downhill = ops.cmp(_len(Lo, Hi, i2, ops), _len(Lo, Hi, i, ops)) < 0
            if (not restrict or downhill) and not ops.eq(y, ops.half(S[e])):
                yf = ops.sub(S[e], y)
                if widest is None or (ops.cmp(yf, Lo[widest]) >= 0
                                      and ops.cmp(yf, Hi[widest]) <= 0):
                    cands.append(((i, -d, yf),
                                  tuple(s - x for s, x in zip(sS[e], ysym))))
            for st, sy in cands:
                k = (st[0], st[1], st[2].v.tobytes())
                if k in seen:
                    continue
                seen.add(k)
                nxt.append((st, sy))
            if len(seen) > maxstates:
                truncated = True
                break
        reached = dstep + 1
        if truncated:
            break
        if not nxt:
            closed = True
            break
        frontier = nxt
    return zhits, len(seen), reached, truncated, closed


def h7(qmax=None, depth=20, maxstates=40_000, qs=(9, 15, 21, 27, 33, 39),
       restrict=True):
    """THE PROFILE-RESTRICTED ENUMERATION.  `restrict=False` reproduces s464 as calibration."""
    import s464_realisable as s464
    from s462_zz_relation import check_sym, sym_build, sym_path

    res = {"qs": list(qs), "restrict": restrict, "depth": depth, "by_Q": {},
           "hits": [], "witness_admitted": None}
    for Q in qs:
        acc = {"rows": 0, "prongs": 0, "pairs": 0, "form": 0, "triv": 0, "rel": 0,
               "prongs_rel": 0, "prongs_any": 0, "states": 0, "truncated": 0,
               "closed": 0, "rel_rows": []}
        for P in range(1, Q):
            if math.gcd(P, Q) != 1:
                continue
            for eps in (0, 1):
                c = _chain(P, Q, eps)
                if c is None:
                    continue
                ch, lo, hi, Lo, Hi, S, ops = c
                LO, HI = sym_build(ch)
                good, _w = check_sym(ch, lo, hi, LO, HI)
                if not good:
                    continue
                sLo, sHi, sS = sym_path(ch, LO, HI, lo, hi)
                red, order = s464.reduce_map(ch)
                ndist = len(order)
                widest = s464.nested_ok(Lo, Hi, ops)
                acc["rows"] += 1
                row_rel = 0
                for j in range(len(S)):
                    cs = nb.corner_start(Lo, Hi, S, j, ops)
                    if cs is None:
                        continue
                    st0 = cs[0]
                    sym0 = s464.start_sym(Lo, Hi, sLo, sHi, ops, st0[2])
                    if sym0 is None:
                        continue
                    acc["prongs"] += 1
                    acc["pairs"] += 2 * len(Lo)
                    zh, ns, _dr, tr, cl = restricted_reach(
                        Lo, Hi, S, ops, st0, sym0, sLo, sHi, sS,
                        depth, maxstates, widest, restrict)
                    acc["states"] += ns
                    acc["truncated"] += int(tr)
                    acc["closed"] += int(cl)
                    acc["prongs_any"] += int(bool(zh))
                    hit_rel = False
                    for h in zh:
                        kind, ar = s464.classify(h["a"], red, ndist)
                        acc[kind] += 1
                        if kind == "rel":
                            hit_rel = True
                            if len(res["hits"]) < 40:
                                res["hits"].append({"P": P, "Q": Q, "eps": eps, "j": j,
                                                    "steps": h["steps"], "a_red": ar})
                    acc["prongs_rel"] += int(hit_rel)
                    row_rel += int(hit_rel)
                    if (P, Q, eps) == (7, 15, 1) and zh:
                        res["witness_admitted"] = True
                if row_rel:
                    acc["rel_rows"].append({"P": P, "eps": eps, "prongs_rel": row_rel})
        res["by_Q"][str(Q)] = acc
    return res


# ------------------------------------------------------------------ driver

def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    qmax = QMAX
    for a in sys.argv[1:]:
        if a.startswith("--qmax="):
            qmax = int(a.split("=")[1])
    which = args[0] if args else "all"
    arms = {"h0": h0, "h1": h1, "h2": h2, "h3": h3, "h4": h4, "h5": h5, "h6": h6,
            "h7": lambda q: h7(q, restrict=True),
            "h7cal": lambda q: h7(q, restrict=False)}
    run = list(arms) if which == "all" else [which]
    out, t0 = {"qmax": qmax}, time.time()
    for name in run:
        t = time.time()
        out[name] = arms[name](qmax)
        print(f"[{name}] {time.time() - t:.1f}s")
        print(json.dumps(out[name], indent=1, default=str)[:2600])
    out["elapsed"] = round(time.time() - t0, 1)
    with open(OUT, "w") as f:
        json.dump(out, f, indent=1, default=str)
    print(f"\nwrote {OUT}  ({out['elapsed']}s)")


if __name__ == "__main__":
    main()
