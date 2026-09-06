#!/usr/bin/env python3
"""
s459_corner_deep.py -- PRE-REGISTERED.  Queue head (III'), via the ONE move that does not
need P3's gate: run [NCYL-307] H5's CORNER-prong detector DEEP on the rows it left capped.

WHAT THIS IS.  NOT a CP characterisation and NOT a bigger-cap re-run of the separatrix
census.  A deep, exact, RESUMABLE re-walk of the `114` corner prongs that [NCYL-307] H5
left at `CAP` after `200000` steps -- with `14/31 eps=0` and `17/31 eps=1` first, because
those are `2` of the `4` rows [NCYL-312] leaves un-certified at `Q <= 31`.

WHY THIS AND NOT A BIGGER SEPARATRIX CAP.  [NCYL-312]'s `4` survivors are not one kind of
object.  `7/15 eps1` / `8/15 eps0` already carry an INDEPENDENT, TERMINATING,
tolerance-free witness of CP failure -- [NCYL-307] H5's corner prong reaching `Z` in `5`
exact steps, to which [NCYL-262]'s *a cap is not a verdict* does not apply.  `14/31 eps0`
/ `17/31 eps1` carry NO such witness and are ONLY capped.  The detector that produced the
witness was never run past `200000` steps on them: [NCYL-309] says so in as many words --
"[NCYL-307] H5 leaves `111` CORNER prongs `CAP` on these same unscored rows, so *fully
resolved* here means *fully resolved for the POLE half*".  ⇒ so the question is not
whether a bigger cap closes a separatrix; it is whether the terminating witness exists on
the other two rows.  Both outcomes are informative and they say different things:
  * a `Z` fires  ->  `14/31 eps0` / `17/31 eps1` are non-CP on evidence that is NOT a cap,
    and [NCYL-312]'s surviving set is CHARACTERISED at `Q <= 31`: all `4` non-CP rows
    carry a terminating witness, and P3's gate stops being the binding constraint on them;
  * no `Z` fires  ->  that is itself the split, and it is the first thing separating the
    two pairs by something other than trace length ([OPS-041] -- a pattern in the `4` at
    separation factor `1.30` is a pattern in length; a terminating witness is not).
⚠ The implication runs ONE way ([NCYL-307]): no `Z` is NOT evidence of CP.

PRIOR ART: read `rulings.md` [NCYL-307] (H5, the detector and the `4 Z`), [NCYL-308] (the
four routes priced as dead -- none of them is this one; this re-runs the EXISTING detector
deeper, it does not build a new arithmetic invariant), [NCYL-309] (the pole half, its
deep-cap run, its control, and the sentence naming this gap), [NCYL-312] (the `4`
survivors and why they are two different kinds), [NCYL-262] (a cap is not a verdict),
[NCYL-292] (the corner prong IS the fold image of the corner), [OPS-236] (resumable capped
runs -- MANDATORY here), [OPS-041], [OPS-215] (no float arrival test);
grepped 'corner prong', 'corner'+'deep', '111 CORNER', 'Z ending', 'terminating witness'
through `rulings.py --grep` ->
  - NOTHING has run a corner prong past `200000` steps.  [NCYL-309]'s deep-cap run is
    `pm.walk`, i.e. POLE prongs on the quotient path, a different start (`s_j/2` vs an
    interface endpoint) and a different object.
  - [NCYL-308]'s dead list is about ARITHMETIC INVARIANTS (`Lambda_S`, the coefficient-sum
    augmentation, [NCYL-301]'s palindrome, a topological attack).  A deeper cap on the
    existing exact tracer is not on it and is not priced anywhere.
  - ⚠ [NCYL-312]'s un-checked coincidence (`14/31`/`17/31` unresolved on BOTH objects) is
    NOT what this settles; C4 below only records the two objects side by side.

HYPOTHESES / CONTROLS (pre-registered; C1 and C2 are the ones that can fail hard).
  C1  ⇒⇒ WALKER EQUIVALENCE, and it is the arm that could invalidate everything else.
      `fast_orbit` is a flat re-implementation of `s454_self_hit.orbit` (no state history,
      in-place numpy, so a `1e8`-step walk fits in memory and in the session).  It must
      reproduce the REFERENCE walker's end label AND step count on ALL `4802` corner
      prongs at the s456 cap -- i.e. re-derive `4684 ESC / 114 CAP / 4 Z` exactly.
  C2  ⇒⇒ RESUME EQUIVALENCE ([OPS-236], mandatory).  A two-stage escalation `lo -> hi`
      must equal a ONE-SHOT run at `hi` in EVERY scalar -- end, steps, folds, and the
      MERGED audit (`margin`, `exact_calls`) -- with the intermediate state round-tripped
      through JSON, because that is the real path.  ⚠ The trap [OPS-236] names: `margin`
      and `exact_calls` are gathered on the traces you actually walk, so a resumed prong
      must MERGE them (min / sum), never recompute.  The suite must include a prong that
      caps at `lo`, one that caps at BOTH, and one that caps at NEITHER.
  C3  ⇒⇒ THE "IT IS NOT JUST CLOSING EVERYTHING" CONTROL, [NCYL-309]'s design.  Report how
      many prongs remain `CAP` at the deep cap.  If every one of the `114` resolves, the
      deep run has no discrimination left and a `Z` among them means correspondingly less;
      if some stay capped, the ones that DID resolve mean something.
  C4  ⇒ THE FLOAT-GATE AUDIT, quoted whatever it says ([OPS-215]).  TWO numbers, because
      at this depth they are different questions.  (i) `margin`, the smallest float gap
      ever resolved as NON-equal -- an ARRIVAL is always exact, so the gate can only
      mislabel a NON-arrival, and only if `margin` collapses.  (ii) `drift`, the largest
      discrepancy ever observed between the CARRIED float and the one recomputed from the
      exact vector.  ⚠ `drift` is a number no previous corner/pole run has: neither
      `s454_self_hit.orbit` nor `s453_pole_module.walk` refreshes, and [NCYL-309] ran the
      latter to `8141619` folds.  The band test is safe exactly while `drift << margin`.
  H   ⇒⇒ THE QUESTION.  Do the `5 + 5` capped corner prongs on `14/31 eps=0` and
      `17/31 eps=1` terminate at `Z`?

GUARDS.
  G1  the reference walker `sh.orbit` is NOT modified and is not imported for anything but
      C1; `s456_neckgb.corner_start` supplies the start, unmodified.
  G2  a TIE edge / degenerate flank carries no corner prong (`corner_start -> None`):
      counted as SKIPPED, never as a pass.
  G3  nothing here writes to `cell_store` or to any s45x store; this probe owns
      `data/s459_corner_deep/` only.  Per-prong atomic writes ([OPS-074]'s pattern).
  G4  ⚠ [OPS-228]: numpy `int64` never leaves this file as a JSON scalar -- state vectors
      are cast to Python `int` on write.
  G5  ⚠ a stored prong record is REUSED only if its cap is >= the one asked for, and a
      resume happens only on a STRICT escalation (re-running at the same cap is what
      `--force` is for, and resuming there would return the stored answer without walking
      -- i.e. silently not audit).  Same contract as `s457_cp_char.measure`.

RESULTS (s459; `data/s459_corner_deep/` = 32 class rows, `logs/s459_{corner,cover,hardened}.log`).
  C1  ⇒⇒ PASS, and BOTH arms.  `4802/4802` corner prongs agree with `s454_self_hit.orbit`
      on the end label AND on the step count, re-deriving [NCYL-307] H5's
      `4684 ESC / 114 CAP / 4 Z` exactly -- once bit-for-bit (`refresh=0, precise=False`)
      and once in the deep run's hardened configuration (`refresh=512, precise=True`), so
      the two departures below are SCORED at the reference cap, not assumed harmless.
  C2  ⇒⇒ PASS, `5/5` identical, all three regimes (caps at `lo`; caps at BOTH round trips;
      caps at NEITHER, which re-walks exactly `0` steps), prior round-tripped through JSON.
      ⇒⇒ **AND IT EARNED ITS KEEP: IT FAILED FIRST, `1/5`, ON A DEFECT [OPS-236] DOES NOT
      COVER** -- the float-refresh SCHEDULE was keyed on a LOCAL counter, which restarts at
      the resume boundary, so a resumed walk refreshed at different folds than a one-shot
      and reported a different `margin`/`drift`.  Fixed by keying the schedule on the
      CUMULATIVE fold count -- which REMOVES the state rather than carrying it.  → [OPS-238].
  H   ⇒⇒ THE ANSWER, AND IT IS THE NEGATIVE BRANCH -- WHICH IS THE INFORMATIVE ONE HERE.
      `14/31 eps=0` and `17/31 eps=1` have their corner half **EXHAUSTIVELY RESOLVED** at
      cap `1e8`: `14/14 ESC`, `0` capped, longest `34957215` steps (separation `2.86x`).
      **NO terminating witness, and NOT for lack of depth.**  Against `7/15 eps=1` /
      `8/15 eps=0`, which keep their `2 Z` (at `5` steps) AND still hold `2` capped corner
      prongs at `1e8` -- `500x` the s456 cap.  ⇒ the `4` [NCYL-312] survivors are split by
      something that is not a cap on either side, for the first time.
  C3  ⇒⇒ DISCRIMINATION RETAINED, and this is what makes H mean anything.  Of [NCYL-307]
      H5's `114` capped corner prongs, `106` resolve by `1e8` and **every one to `ESC`**;
      `8` stay capped (`2+2` on the witnessed pair, `1` each on `4/33 eps0`, `7/33 eps1`,
      `26/33 eps0`, `29/33 eps1`).  ⇒⇒ **AND THE `Z` COUNT NEVER MOVES: exactly `4`, on the
      same two rows, across a `500x` cap raise.**  The detector does not manufacture
      witnesses when you look harder.
  C4  ⚠⚠ THE AUDIT THAT FORCED THE HARDENING, quoted because it came out BADLY.  Worst
      measured float drift `8.015e-09` against a worst margin of `6.320e-10` -- a
      `margin/drift` ratio of **`0.08`** at `7/33 eps1 j=3`, i.e. a band test that COULD
      have taken the wrong branch.  ⇒ so the gate was hardened per `s450_chain_maps.safety`'s
      own prescription (`precise=True`: in-gate signs resolved in `80` digits), which is
      nearly free -- `186` gate hits in `5.7e7` steps.
  C5  ⇒ AND THE HARDENING CHANGED NOTHING, WHICH IS ALSO A RESULT: exact gate vs float gate,
      **`370` prongs identical, `0` differ** (`34` not comparable -- different caps between
      the two stores).  A separate arm re-walked the four survivor rows with the reference's
      UNREFRESHED float (`refresh=-512`, audit-only) : `40/40` SAME, and the unrefreshed
      drift is the same ORDER as the refreshed per-window drift (`3.7e-09` vs `4.6e-09` at
      `7/15 eps1 j=1`, `1e8` steps), i.e. the error does NOT accumulate -- the reflection
      `y -> S_e - y` resets it.  ⇒⇒ **so `s454_self_hit.orbit` / `s453_pole_module.walk` are
      vindicated at depth ON THIS EVIDENCE and do NOT need retrofitting** ([NCYL-309]'s
      `8141619`-fold run included).  → [OPS-239].
  ⚠⚠ NOT EVIDENCE.  `14/31 eps0` and `17/31 eps1` returning IDENTICAL numbers is FORCED by
      the leg swap ([OPS-235]), not a corroboration.  And *no `Z`* is NOT *CP*: the
      implication runs one way ([NCYL-307]) -- what H buys is the removal of a TRUNCATION
      explanation, not a verdict.  The verdict on those two rows came from the other
      instrument ([NCYL-313]).

Run: PYTHONPATH=.:engine:archive/scripts_2026-07:archive/scripts_2026-08:probes \\
     .venv/bin/python3.13 probes/s459_corner_deep.py \\
       [c1|c2|deep|summary|all] [--cap=N] [--qmax=N] [--budget=SECS] [--force]
"""
from __future__ import annotations

import glob
import json
import math
import os
import sys
import time

import mpmath as mp
import numpy as np

import s450_chain_maps as cm
import s453_pole_module as pm
import s454_self_hit as sh
import s456_neckgb as gb

OUTDIR = "data/s459_corner_deep"
SUMMARY = "data/s459_corner_deep.json"
REF_CAP = 200_000                      # the s456 H5 cap -- C1's reference point
TOL = cm.TOL

# ⇒ THE TWO ROWS THE WHOLE ITEM IS ABOUT ([NCYL-312]'s survivors with no witness).
TARGETS = [(14, 31, 0), (17, 31, 1)]
# ⇒ and the two that ALREADY carry one, as the positive control for the detector itself.
WITNESSED = [(7, 15, 1), (8, 15, 0)]


# ---------------------------------------------------------------- the exact sign
#
# ⇒⇒ [OPS-215]'s FIX, AND `s450_chain_maps.safety`'s OWN PRESCRIPTION, FINALLY EXERCISED.
# The band test `<`/`>` reads a float whose error C4 measures as `drift`; on the s459
# sweep the worst `margin/drift` came out at `0.08` -- i.e. a comparison that COULD have
# taken the wrong branch.  `safety`'s docstring says what to do and why it is cheap: the
# gated comparisons are RARE (here `186` exact calls in `5.7e7` steps), so resolving them
# in high precision costs nothing measurable.  It was never written because it had never
# been needed; C4 is what made it needed.
_MPIM = {}


def _mp_imvals(Q, d, dps=80):
    key = (Q, d, dps)
    if key not in _MPIM:
        with mp.workdps(dps):
            _MPIM[key] = [mp.sin(k * mp.pi / (2 * Q)) / 4 for k in range(d)]
    return _MPIM[key]


def exact_sign(dv, mpim, dps=80):
    """Sign of the real value of a NONZERO coordinate difference, in `dps` digits.

    Returns `-1`/`+1`, or `0` if the value is indistinguishable from zero at this
    precision -- which would mean the vector representation is not faithful on the
    reachable set, and is reported rather than swallowed.
    """
    with mp.workdps(dps):
        tot = mp.mpf(0)
        for k in np.nonzero(dv)[0]:
            tot += int(dv[k]) * mpim[k]
        if abs(tot) < mp.mpf(10) ** (-(dps - 20)):
            return 0
        return 1 if tot > 0 else -1


# ---------------------------------------------------------------- the flat walker

class Row:
    """Precomputed float/vector views of one class row's quotient path."""

    __slots__ = ("P", "Q", "eps", "n", "Lo_v", "Lo_f", "Hi_v", "Hi_f",
                 "S_v", "S_f", "hS_v", "Lo", "Hi", "S", "ops", "ch", "imvals", "mpim")

    def __init__(self, P, Q, eps, ch, Lo, Hi, S, ops):
        self.P, self.Q, self.eps, self.ch, self.ops = P, Q, eps, ch, ops
        self.Lo, self.Hi, self.S = Lo, Hi, S
        self.n = len(Lo)
        self.Lo_v = [np.asarray(c.v) for c in Lo]
        self.Hi_v = [np.asarray(c.v) for c in Hi]
        self.Lo_f = [float(c.f) for c in Lo]
        self.Hi_f = [float(c.f) for c in Hi]
        self.S_v = [np.asarray(c.v) for c in S]
        self.S_f = [float(c.f) for c in S]
        self.imvals = ch.ring.imvals          # `f == v @ imvals`; the refresh needs it
        self.mpim = _mp_imvals(Q, ch.ring.d)  # the same basis in 80 digits, for the gate
        # ⚠ `cm.half` ASSERTS an even vector.  The reference calls it only on the edges it
        # actually folds at, so precomputing eagerly would raise on an edge the reference
        # never reaches -- a divergence in the WRONG direction (a crash, not a verdict).
        # Hence lazy, memoised, and identical in reach to the reference.
        self.hS_v = [None] * len(S)

    def half_S(self, e):
        h = self.hS_v[e]
        if h is None:
            h = self.hS_v[e] = np.asarray(cm.half(self.S[e]).v)
        return h


def load_row(P, Q, eps):
    c = gb._chain(P, Q, eps)
    if c is None:
        return None
    ch, Lo, Hi, S, ops = c
    return Row(P, Q, eps, ch, Lo, Hi, S, ops)


def fast_orbit(row, st, cap, steps0=0, folds0=0, margin0=1.0, exact0=0, drift0=0.0,
               refresh=512, precise=True):
    """Flat, resumable forward orbit.  Semantics are `s454_self_hit.step` VERBATIM --
    C1 is what asserts that in code rather than in this docstring ([OPS-043]).

    `st` is `(i, d, yv, yf)` with `yv` a numpy int vector and `yf` its float.
    Returns `(end, steps, folds, state, margin, exact_calls, drift)`; `steps`/`folds` are
    CUMULATIVE (so a resumed prong reports the whole walk), and `margin`/`exact_calls`/
    `drift` are MERGED with the prior run's (min / sum / max), never recomputed
    ([OPS-236]'s audit-accumulator trap).  `state` is `None` once the prong has ended.

    ⚠⚠ `refresh` IS NOT COSMETIC AND IT IS A DEPARTURE FROM THE REFERENCE.  `Coord.f`
    is carried by float SUBTRACTION (`Coord.__sub__`), so over a deep walk it drifts away
    from the exact value `v @ imvals`, and the band test `<`/`>` reads that float --
    [OPS-215]'s failure mode exactly.  `s450_chain_maps` recomputes the float from the
    exact vector every `512` steps for precisely this reason; `s454_self_hit.orbit` and
    `s453_pole_module.walk` do NOT, and both have been run deep.  Here the float is
    refreshed every `refresh` folds AND the discrepancy is recorded as `drift`, so the
    audit is a measurement rather than an assumption.  `refresh=0` reproduces the
    reference bit-for-bit, which is what C1 uses as its second arm; `refresh<0` is
    AUDIT-ONLY -- it measures the reference's drift at `|refresh|`-fold intervals WITHOUT
    correcting it, which is the only way to price what an unrefreshed deep run is worth.
    """
    i, d, yv, yf = st
    n = row.n
    Lo_v, Lo_f, Hi_v, Hi_f = row.Lo_v, row.Lo_f, row.Hi_v, row.Hi_f
    S_v, S_f, half_S = row.S_v, row.S_f, row.half_S
    imv, mpim = row.imvals, row.mpim
    margin, exact, drift = margin0, exact0, drift0
    steps, folds = steps0, folds0
    audit_every = refresh if refresh > 0 else -refresh
    budget = cap - steps0
    array_equal = np.array_equal
    for _ in range(budget):
        i2 = i + d
        if i2 < 0 or i2 >= n:
            return "ESC", steps + 1, folds, None, margin, exact, drift
        steps += 1
        # --- cmp(y, Lo[i2]) ------------------------------------------------------
        g = yf - Lo_f[i2]
        c0 = -1 if g < 0 else 1
        if -TOL < g < TOL:
            exact += 1
            if array_equal(yv, Lo_v[i2]):
                return "Z", steps, folds, None, margin, exact, drift
            if precise:                         # ⇐ the float sign is NOT trusted in-gate
                sg = exact_sign(yv - Lo_v[i2], mpim)
                if sg == 0:                     # ⚠ loud, never swallowed: a nonzero
                    raise RuntimeError(         # vector with a zero value would mean the
                        f"non-faithful coordinate at {row.P}/{row.Q} eps{row.eps} "
                        f"step {steps}")        # representation is not faithful
                c0 = sg
        a = -g if g < 0 else g
        if a < margin:
            margin = a
        # --- cmp(y, Hi[i2]) ------------------------------------------------------
        g = yf - Hi_f[i2]
        c1 = -1 if g < 0 else 1
        if -TOL < g < TOL:
            exact += 1
            if array_equal(yv, Hi_v[i2]):
                return "Z", steps, folds, None, margin, exact, drift
            if precise:
                sg = exact_sign(yv - Hi_v[i2], mpim)
                if sg == 0:
                    raise RuntimeError(
                        f"non-faithful coordinate at {row.P}/{row.Q} eps{row.eps} "
                        f"step {steps}")
                c1 = sg
        a = -g if g < 0 else g
        if a < margin:
            margin = a
        # --- cross / fold --------------------------------------------------------
        if c0 > 0 and c1 < 0:
            i = i2
            continue
        e = i if i < i2 else i2
        if array_equal(yv, half_S(e)):
            return "pole", steps, folds, None, margin, exact, drift
        yv = S_v[e] - yv
        yf = S_f[e] - yf
        d = -d
        folds += 1
        # ⚠⚠ THE SCHEDULE IS KEYED ON THE CUMULATIVE FOLD COUNT, NOT ON A LOCAL COUNTER.
        # A local counter restarts at the resume boundary, so a resumed walk refreshes at
        # DIFFERENT folds than a one-shot and reports a different `margin`/`drift` --
        # [OPS-236]'s audit-accumulator trap, one level deeper: it is not only the
        # accumulators that must survive a resume, it is the SCHEDULE that feeds them.
        # C2 caught this; keying on `folds` removes the state instead of carrying it.
        if refresh and folds % audit_every == 0:
            exact_f = float(yv @ imv)          # the value the float is SUPPOSED to hold
            a = yf - exact_f
            if a < 0:
                a = -a
            if a > drift:
                drift = a
            if refresh > 0:
                yf = exact_f                   # kill the accumulation, do not just audit
    return "CAP", steps, folds, (i, d, yv, yf), margin, exact, drift


def start_state(row, j):
    """`corner_start` -> the flat state, or `None` (tie edge / degenerate flank)."""
    cs = gb.corner_start(row.Lo, row.Hi, row.S, j, row.ops)
    if cs is None:
        return None
    i, d, y = cs[0]
    return (i, d, np.asarray(y.v), float(y.f))


# ---------------------------------------------------------------- the store

def _path(P, Q, eps):
    return os.path.join(OUTDIR, f"{P}_{Q}_{eps}.json")


def _write(P, Q, eps, rec):
    os.makedirs(OUTDIR, exist_ok=True)
    tmp = _path(P, Q, eps) + ".tmp"
    json.dump(rec, open(tmp, "w"))
    os.replace(tmp, _path(P, Q, eps))            # atomic: a kill costs the in-flight row


def _ser(state):
    if state is None:
        return None
    i, d, yv, yf = state
    return [int(i), int(d), [int(x) for x in yv], float(yf)]     # G4


def _deser(s):
    if s is None:
        return None
    i, d, yv, yf = s
    return (i, d, np.asarray(yv, dtype=np.int64), yf)


def measure_row(P, Q, eps, cap, force=False, verbose=True):
    """Every corner prong of one class row, at `cap`, RESUMING the capped ones.

    ⇒⇒ [OPS-236], all three levels: PRONG-level (a prong already resolved at a lower cap
    is not re-walked at all), and STEP-level (a prong that capped resumes from its stored
    walk state, so the escalation costs `hi - lo`, not `hi`).  The audit accumulators are
    MERGED, not recomputed.
    """
    rec = {}
    if os.path.exists(_path(P, Q, eps)):
        rec = json.load(open(_path(P, Q, eps)))
        if not force and rec.get("cap", 0) >= cap:
            return rec
    row = load_row(P, Q, eps)
    if row is None:
        return None
    old = {str(p["j"]): p for p in rec.get("prongs", [])}
    prongs, skipped, t0 = [], 0, time.time()
    for j in range(len(row.S)):
        st0 = start_state(row, j)
        if st0 is None:
            skipped += 1                                            # G2
            continue
        prev = old.get(str(j))
        if prev is not None and prev["end"] != "CAP":
            prongs.append(prev)                                     # resolved: never re-walk
            continue
        if prev is not None and prev.get("state") is not None and prev["cap"] < cap:
            st = _deser(prev["state"])                              # STEP-level resume
            end, steps, folds, state, margin, exact, drift = fast_orbit(
                row, st, cap, prev["steps"], prev["folds"],
                prev["margin"], prev["exact_calls"], prev.get("drift", 0.0))
        else:
            end, steps, folds, state, margin, exact, drift = fast_orbit(row, st0, cap)
        prongs.append({"j": j, "end": end, "steps": steps, "folds": folds,
                       "state": _ser(state), "margin": margin,
                       "exact_calls": exact, "drift": drift, "cap": cap})
    out = {"P": P, "Q": Q, "eps": eps, "cap": cap, "skipped": skipped,
           # ⚠ the walker CONFIGURATION is provenance, not a measurement: it is recorded,
           # never defaulted or inferred ([OPS-237]).
           "walker": {"refresh": 512, "precise": True},
           "prongs": prongs, "secs": round(time.time() - t0, 1)}
    # ⚠ the per-cap trajectory, on EVERY path including the first ([OPS-236]).
    out["runs"] = dict(rec.get("runs", {}))
    out["runs"][str(cap)] = {
        "capped": sum(p["end"] == "CAP" for p in prongs),
        "prongs": len(prongs),
        "ends": _tally(prongs),
        "max_steps": max((p["steps"] for p in prongs), default=0)}
    _write(P, Q, eps, out)
    if verbose:
        print(f"  {P}/{Q} eps{eps}: {_tally(prongs)}  "
              f"capped={out['runs'][str(cap)]['capped']}/{len(prongs)}  "
              f"max_steps={out['runs'][str(cap)]['max_steps']}  [{out['secs']}s]",
              flush=True)
    return out


def _tally(prongs):
    t = {}
    for p in prongs:
        t[p["end"]] = t.get(p["end"], 0) + 1
    return t


# ---------------------------------------------------------------- C1

def c1(qmax=33, cap=REF_CAP, refresh=0, precise=False, verbose=True):
    """⇒⇒ WALKER EQUIVALENCE against `s454_self_hit.orbit` on every corner prong.

    Must re-derive [NCYL-307] H5's `4684 ESC / 114 CAP / 4 Z` and agree PRONG BY PRONG on
    the end label and the step count.  This is the arm that can invalidate everything
    below it, so it runs over the whole box, not a sample.

    ⚠ `refresh=0` is the bit-for-bit arm (the reference does not refresh its float).  The
    driver runs it AGAIN at `refresh=512` -- the setting the deep run uses -- so that the
    departure documented on `fast_orbit` is SCORED at the reference cap rather than
    assumed harmless.
    """
    res = {"prongs": 0, "end_ok": 0, "steps_ok": 0, "bad": [], "ends": {},
           "ref_ends": {}, "rows": 0, "skipped": 0}
    t0 = time.time()
    for P, Q, eps in gb._rows(qmax):
        row = load_row(P, Q, eps)
        if row is None:
            continue
        res["rows"] += 1
        for j in range(len(row.S)):
            cs = gb.corner_start(row.Lo, row.Hi, row.S, j, row.ops)
            if cs is None:
                res["skipped"] += 1
                continue
            res["prongs"] += 1
            sts, kinds, _e, ref_end = sh.orbit(row.Lo, row.Hi, row.S, cs[0],
                                               row.ops, cap)
            end, steps, _f, _st, _m, _x, _dr = fast_orbit(
                row, start_state(row, j), cap, refresh=refresh, precise=precise)
            res["ends"][end] = res["ends"].get(end, 0) + 1
            res["ref_ends"][ref_end] = res["ref_ends"].get(ref_end, 0) + 1
            ok_e = (end == ref_end)
            ok_s = (steps == len(kinds))
            res["end_ok"] += ok_e
            res["steps_ok"] += ok_s
            if not (ok_e and ok_s) and len(res["bad"]) < 12:
                res["bad"].append([P, Q, eps, j, end, ref_end, steps, len(kinds)])
    res["secs"] = round(time.time() - t0, 1)
    if verbose:
        print(f"  C1  fast vs reference on {res['prongs']} corner prongs "
              f"({res['rows']} rows, cap {cap}, refresh={refresh}, precise={precise})")
        print(f"      end label agrees : {res['end_ok']}/{res['prongs']}")
        print(f"      step count agrees: {res['steps_ok']}/{res['prongs']}")
        print(f"      fast ends {res['ends']}   reference ends {res['ref_ends']}")
        for b in res["bad"]:
            print("      MISMATCH:", b)
        print(f"      VERDICT: {'PASS' if res['end_ok'] == res['steps_ok'] == res['prongs'] else '*** FAIL ***'}"
              f"  [{res['secs']}s]", flush=True)
    return res


# ---------------------------------------------------------------- C2

def c2(cases=None, verbose=True):
    """⇒⇒ RESUME EQUIVALENCE ([OPS-236]).  Two-stage `lo -> hi` == one shot at `hi`, in
    end, steps, folds AND the merged audit -- with the intermediate state ROUND-TRIPPED
    THROUGH JSON, because that is the real path.

    The suite must contain a prong that caps at `lo` only, one that caps at BOTH (so the
    state survives two round trips) and one that caps at NEITHER (which checks the resume
    re-walks NOTHING).  ⚠ A case where the low cap does not bite contributes no saving and
    is still a real test.
    """
    if cases is None:
        cases = [((14, 31, 0), 2, 20_000, 400_000),      # caps at lo, resolves(?) at hi
                 ((14, 31, 0), 2, 20_000, 60_000),       # caps at BOTH round trips
                 ((7, 15, 1), 2, 20_000, 200_000),       # the `Z` prong: caps at NEITHER
                 ((14, 31, 0), 6, 20_000, 200_000),      # ESC at 77483: lo bites, hi does not
                 ((8, 15, 0), 3, 50_000, 500_000)]
    res = {"cases": [], "identical": 0, "n": 0, "one_shot_steps": 0, "escalation_steps": 0}
    for (P, Q, eps), j, lo, hi in cases:
        row = load_row(P, Q, eps)
        st0 = start_state(row, j)
        if st0 is None:
            continue
        res["n"] += 1
        a = fast_orbit(row, st0, hi)                                   # (A) one shot
        b1 = fast_orbit(row, st0, lo)                                  # (B) low ...
        ser = json.loads(json.dumps(_ser(b1[3])))                      # ... through JSON
        if ser is None:
            b = b1                                     # resolved at lo: nothing to resume
            walked = 0
        else:
            b = fast_orbit(row, _deser(ser), hi, b1[1], b1[2], b1[4], b1[5],
                           b1[6])
            walked = b[1] - b1[1]
        same = (a[0] == b[0] and a[1] == b[1] and a[2] == b[2]
                and a[4] == b[4] and a[5] == b[5] and a[6] == b[6])
        res["identical"] += same
        res["one_shot_steps"] += a[1]
        res["escalation_steps"] += walked
        res["cases"].append({"row": [P, Q, eps], "j": j, "lo": lo, "hi": hi,
                             "capped_at_lo": b1[0] == "CAP", "end": a[0],
                             "identical": bool(same),
                             "one_shot": [a[0], a[1], a[2], a[4], a[5], a[6]],
                             "resumed": [b[0], b[1], b[2], b[4], b[5], b[6]]})
        if verbose:
            tag = "IDENTICAL" if same else "*** DIFFERS ***"
            print(f"  {P}/{Q} eps{eps} j={j}  cap {lo}->{hi}  "
                  f"{'caps at lo' if b1[0]=='CAP' else 'resolves at lo'}  {tag}  "
                  f"| re-walked {walked} of {a[1]} steps", flush=True)
    if verbose:
        # ⚠ the honest metric is the ESCALATION's cost against a one-shot at `hi`, since
        # the low-cap run already exists in the store; a fresh two-stage ladder saves
        # nothing by construction and quoting `lo + tail` would hide that ([OPS-236]'s
        # `52.7%` caveat -- the ratio is a function of where the cap sits, not a law).
        print(f"  C2  {res['identical']}/{res['n']} identical; the ESCALATION walked "
              f"{res['escalation_steps']} of the {res['one_shot_steps']} steps a one-shot "
              f"at `hi` would cost")
        print(f"      VERDICT: {'PASS' if res['identical'] == res['n'] else '*** FAIL ***'}",
              flush=True)
    return res


# ---------------------------------------------------------------- the deep run

def deep(cap=100_000_000, qmax=33, budget_s=None, force=False, only=None,
         verbose=True):
    """Escalate every row that still holds a capped corner prong, TARGETS first.

    ⚠ The queue is rebuilt from the store on every call, so a killed run is resumed by
    re-issuing the same command; and a row already resolved at a lower cap is skipped
    entirely (`measure_row`'s prong-level guard).  STOP file: `<OUTDIR>/STOP`.
    """
    os.makedirs(OUTDIR, exist_ok=True)
    todo = list(only) if only else _queue(qmax)
    t0 = time.time()
    if verbose:
        print(f"deep  {len(todo)} rows with capped corner prongs -> cap {cap}"
              f"   STOP file: {OUTDIR}/STOP", flush=True)
    done = []
    for P, Q, eps in todo:
        if os.path.exists(os.path.join(OUTDIR, "STOP")):
            print("  [STOP file]", flush=True)
            break
        if budget_s and time.time() - t0 > budget_s:
            print("  [budget]", flush=True)
            break
        measure_row(P, Q, eps, cap, force=force, verbose=verbose)
        done.append((P, Q, eps))
    return done


def _queue(qmax):
    """Rows holding a `CAP` corner prong at the shallow cap -- TARGETS and the two
    WITNESSED control rows first, then the rest by `Q`."""
    seed = list(TARGETS) + list(WITNESSED)
    rest = []
    for P, Q, eps in gb._rows(qmax):
        if (P, Q, eps) in seed:
            continue
        f = _path(P, Q, eps)
        if os.path.exists(f):
            r = json.load(open(f))
            if any(p["end"] == "CAP" for p in r["prongs"]):
                rest.append((P, Q, eps))
            continue
        row = load_row(P, Q, eps)
        if row is None:
            continue
        for j in range(len(row.S)):
            st = start_state(row, j)
            if st is None:
                continue
            if fast_orbit(row, st, REF_CAP)[0] == "CAP":
                rest.append((P, Q, eps))
                break
    return seed + sorted(rest, key=lambda k: (k[1], k[0], k[2]))


# ---------------------------------------------------------------- summary

def summary(verbose=True):
    """C3 + C4 + H, read off the store."""
    rows, ends, capped, hits = 0, {}, [], []
    worst_margin, worst_at, maxsteps, worst_drift = 1.0, None, 0, 0.0
    worst_ratio, worst_ratio_at = float("inf"), None
    caps = {}
    for f in sorted(glob.glob(os.path.join(OUTDIR, "*.json"))):
        r = json.load(open(f))
        rows += 1
        caps[r["cap"]] = caps.get(r["cap"], 0) + 1
        for p in r["prongs"]:
            ends[p["end"]] = ends.get(p["end"], 0) + 1
            maxsteps = max(maxsteps, p["steps"])
            if p["margin"] < worst_margin:
                worst_margin, worst_at = p["margin"], [r["P"], r["Q"], r["eps"], p["j"]]
            worst_drift = max(worst_drift, p.get("drift", 0.0))
            # ⚠ the safety denominator is the MEASURED float error, not `TOL`.  `TOL` is
            # the width of the gate that decides WHEN to consult the exact vector; what
            # can flip a band test is the drift.  `s450_chain_maps.safety` uses a BOUND
            # (`eps * d * maxcoef`); `drift` is the same quantity measured, so it is the
            # honest denominator here.
            if p.get("drift", 0.0) > 0:
                rat = p["margin"] / p["drift"]
                if rat < worst_ratio:
                    worst_ratio = rat
                    worst_ratio_at = [r["P"], r["Q"], r["eps"], p["j"]]
            if p["end"] == "CAP":
                capped.append([r["P"], r["Q"], r["eps"], p["j"], p["steps"], p["cap"]])
            elif p["end"] != "ESC":
                hits.append([r["P"], r["Q"], r["eps"], p["j"], p["end"], p["steps"]])
    tgt = {}
    for P, Q, eps in TARGETS + WITNESSED:
        f = _path(P, Q, eps)
        if os.path.exists(f):
            r = json.load(open(f))
            tgt[f"{P}/{Q} eps{eps}"] = {"cap": r["cap"], "ends": _tally(r["prongs"]),
                                        "runs": sorted(r["runs"])}
    out = {"rows": rows, "ends": ends, "caps": caps, "capped": capped, "hits": hits,
           "max_steps": maxsteps, "worst_margin": worst_margin,
           "worst_margin_at": worst_at, "worst_drift": worst_drift,
           "safety": worst_ratio, "safety_at": worst_ratio_at,
           "targets": tgt}
    if verbose:
        print(f"  rows in store {rows}   caps {caps}")
        print(f"  C3  ends {ends};  still CAP: {len(capped)}"
              f"   {'<-- discrimination retained' if capped else '*** everything closed ***'}")
        print(f"  C4  worst float margin {worst_margin:.3e} at {worst_margin_str(worst_at)}"
              f";  worst measured drift {worst_drift:.3e}")
        print(f"      worst margin/drift ratio {worst_ratio:.2f}x at "
              f"{worst_margin_str(worst_ratio_at)}   (gate TOL {TOL})")
        print(f"      max steps walked {maxsteps}")
        print("  H   the four [NCYL-312] survivors:")
        for k, v in tgt.items():
            print(f"        {k:>14}  {v['ends']}   cap {v['cap']}  runs {v['runs']}")
        for h in hits:
            print("      NON-ESC:", h)
    return out


def worst_margin_str(at):
    return "n/a" if at is None else f"{at[0]}/{at[1]} eps{at[2]} j={at[3]}"


# ---------------------------------------------------------------- driver

def main():
    mode = sys.argv[1] if len(sys.argv) > 1 and not sys.argv[1].startswith("-") else "all"
    kw = {}
    for a in sys.argv[1:]:
        for flag, key, cast in (("--cap=", "cap", int), ("--qmax=", "qmax", int),
                                ("--budget=", "budget_s", float)):
            if a.startswith(flag):
                kw[key] = cast(a.split("=")[1])
    force = "--force" in sys.argv
    out = {}
    if mode in ("all", "c1"):
        out["c1"] = c1(kw.get("qmax", 33), refresh=0, precise=False)
        print("      -- again in the DEEP RUN's configuration (refresh + exact gate) --")
        out["c1_hardened"] = c1(kw.get("qmax", 33), refresh=512, precise=True)
    if mode in ("all", "c2"):
        out["c2"] = c2()
    if mode in ("all", "deep"):
        deep(cap=kw.get("cap", 100_000_000), qmax=kw.get("qmax", 33),
             budget_s=kw.get("budget_s"), force=force)
    if mode in ("all", "deep", "summary"):
        out["summary"] = summary()
    if mode in ("all", "c1", "c2", "summary"):
        json.dump(out, open(SUMMARY, "w"), indent=1)
        print(f"  -> {SUMMARY}")


if __name__ == "__main__":
    main()
