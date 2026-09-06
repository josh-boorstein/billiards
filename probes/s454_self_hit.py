#!/usr/bin/env python3.13
"""s454_self_hit.py -- PRE-REGISTERED.  Queue item (gamma): the SELF-HIT residual that
[NCYL-298] left open -- *no interior pole prong returns to its OWN pole*.  It closes, and
it closes on REVERSIBILITY ALONE: no lattice, no cyclotomics, no arithmetic of any kind.

PRIOR ART: `rulings.py --grep 'reversib'` -> NOTHING; `'palindrome'` -> [NCYL-294] (fold
counts EVEN on ESCAPING prongs, "the `rho` symmetry making the trajectory a palindrome
about its Fix crossing" -- a DIFFERENT involution and a different conclusion, and it is
listed there as NOT evidence), [NCYL-283]; `'self-hit'` / `'own pole'` -> [NCYL-298] only,
which states the residual and says the module argument CANNOT see it (`0 in 2*Lambda`
always).  `'symmetric point'` -> nothing.  The candidate technique named in the s453
handoff is *"the reversibility structure (`T^{-1} = RTR`, and the pole state is a symmetric
point, so a self-hit is a closed orbit through a symmetric point)"* -- so the ROUTE is
prior art (it is the queue item's own suggestion); the LEMMA and its proof are not.
The walker, the path and the fold rule are [NCYL-294]/[NCYL-292] and are USED verbatim
(`s453_pole_module.walk` is the reference implementation this file's tracer is scored
against).  `p`'s definition and the `p != 0` shape are [NCYL-283]/[NCYL-287].
⚠ [NCYL-296] bars routing this through CF depth / renormalization: does not bite, nothing
  below reads a CF digit.  ⚠ [NCYL-201] warns that a group-action REPHRASING of loneness
  is a restatement: does not bite either -- this is class-blind and rephrases nothing.
⚠⚠ [OPS-227] (the `Fix` pun): checked.  The involution below is `R(i,d,y) = (i,-d,y)`,
  TIME REVERSAL on the quotient path.  It is NOT [NCYL-294]'s `iota` (whose quotient the
  path already is) and NOT [NCYL-200]'s `Fix(r)`.  Three different involutions; this one
  is new here and is deliberately called `R`, not `Fix`.

==========================================================================================
THE THEOREM (proved here).   NO INTERIOR POLE PRONG RETURNS TO ITS OWN POLE.

Work in [NCYL-294]'s quotient path `J_0 ... J_m`, `m = (Q-1)/2`.  A STATE is `(i, d, y)`:
node `i in [0,m]`, direction `d = +-1`, transverse coordinate `y in int J_i`.  The step
map `T` (verbatim `s453_pole_module.walk`) is, with `i2 = i + d` and `e = min(i, i2)`:

    i2 out of range              -> ESCAPE (the prong reaches an end = a `Fix` component)
    y in int J_{i2}              -> CROSS:  (i2,  d, y)          [`i` moves, `d` and `y` fixed]
    y an endpoint of J_{i2}      -> Z-hit  (the prong runs into a cone point)
    otherwise, y == s_e/2        -> POLE HIT at edge `e`         [the prong ends on a pole]
    otherwise                    -> FOLD:   (i, -d, s_e - y)     [`i` fixed, `d` flips]

⇒⇒ STEP 1 -- REVERSIBILITY.  Let `R(i, d, y) = (i, -d, y)`.  Then `T^{-1} = R T R` on every
state of a Z-free orbit.  For a CROSS this is immediate.  For a FOLD: adjacent intervals
share exactly ONE endpoint (H3), so writing the shared one as a common `lo` (the `hi` case
is the mirror), blocking at `e` means `y in J_i \\ J_{i2}`, i.e. `hi_{i2} < y <= hi_i`, and
`s_e = hi_i + hi_{i2}` gives `s_e - y in [hi_{i2}, hi_i)` -- **still in `J_i`, still outside
`J_{i2}`.**  So the reflected point is blocked by the same edge, and applying the fold twice
returns `y`.  ⇒ the backward step of a fold is a fold at the same edge.                [H1]

⇒⇒ STEP 2 -- THE PRONG'S TWO ENDS ARE `R`-IMAGES OF EACH OTHER.  `s_e/2` lies in the
interior of EXACTLY ONE of the two intervals at edge `e`, namely the LARGER (H3).  Call it
`i_0`.  The prong of pole `j` starts at `u_0 = (i_0, d_0, s_j/2)` with `d_0` pointing AWAY
from edge `j`.  If it self-hits, the hit is detected at some state `u_M` with `y = s_j/2`
and blocked by edge `j` -- so `u_M`'s node also contains `s_j/2` in its interior, forcing
node `= i_0`, and its direction points TOWARD edge `j`, forcing `d = -d_0`.

    ⇒⇒   u_M  =  R(u_0).                                                              [H3]

⇒⇒ STEP 3 -- THE ORBIT IS A PALINDROME.  Put `w_k := R(u_{M-k})`.  Then `w_0 = R(u_M) =
R(R(u_0)) = u_0`, and `T(w_k) = T(R u_{M-k}) = R T^{-1}(u_{M-k}) = R(u_{M-k-1}) = w_{k+1}`
for `k <= M-1` by Step 1.  Two forward orbits with the same initial state coincide, so

    ⇒⇒   u_k  =  R(u_{M-k})        for every `k = 0 .. M`.

⇒⇒ STEP 4 -- THE MIDDLE IS A POLE, WHICH IS THE CONTRADICTION.
  * `M` EVEN (including `M = 0`): `k = M/2` gives `u_{M/2} = R(u_{M/2})`, i.e. `d = -d`.
    IMPOSSIBLE.
  * `M` ODD: `k = (M-1)/2` gives `u_{(M+1)/2} = R(u_{(M-1)/2})` -- same node, same `y`,
    OPPOSITE direction.  A CROSS changes the node and keeps `d`; only a FOLD flips `d`, and
    a fold also sends `y |-> s_e - y`.  Equality of the `y`'s forces `s_e - y = y`, i.e.
    **`y = s_e/2`: the state at index `(M-1)/2` is sitting on the pole of edge `e`.**  But
    then `T` reports a POLE HIT there and the prong ENDS at index `(M-1)/2 < M`, so it never
    reaches `u_M`.  CONTRADICTION.
                                                                                        ∎
⇒⇒ WHY IT NEEDS NO ARITHMETIC, AND WHY THAT IS THE POINT.  Nothing above uses `Q`, `P`, the
lengths, the lattice `Lambda` or any cyclotomic fact -- only that `T` is reversible and that
the pole is the fold's FIXED POINT.  ⇒ it holds for EVERY odd `Q`, BOTH CLASSES, and for
arbitrary (even non-arithmetic) interval data, which is exactly what H4 measures: `18640`
pole hits in a random rational family, **`0`** of them self-hits.  ⇒ It also RETROACTIVELY
EXPLAINS [NCYL-298] H6's `0 of 336`, which s453 correctly refused to read as evidence of
emptiness -- it was a theorem, not a small sample.

⇒⇒⇒ THE CONSEQUENCE, AND IT IS THE HEADLINE.  [NCYL-298] + [NCYL-300] exclude a connection
from interior pole `j` to a DIFFERENT interior pole `j'`; the above excludes `j' = j`.
⇒ **THERE IS NO INTERIOR POLE-TO-POLE SADDLE CONNECTION AT ALL.**  [NCYL-283]'s final
derivation runs Gauss-Bonnet on the region an `iota`-swapped pair cuts off inside the disk
and concludes it must contain *"exactly TWO simple poles, i.e. an `R`-`R` saddle connection
interior to the disk"* -- so `p != 0` IMPLIES such a connection, and the contrapositive is

    ⇒⇒⇒   p = 0,  BOTH CLASSES, EVERY ODD `Q`.

whence [NCYL-283](4)'s `Sum n_sigma = Q - h + v - 2p` reads `Q-1`/`Q`/`Q+1` -- [NCYL-093]'s
COVERING LAW, open since s320 -- and (M1), and [NCYL-287]'s `C = n_int + 1`.

⚠⚠⚠ SCOPE -- READ THIS BEFORE QUOTING ANY OF THE ABOVE.  `p = 0` IS PROVED *MODULO ITS
  INPUTS*, AND TWO OF THEM ARE NOT THEOREMS:
  (i)  ⇒⇒ [NCYL-294]'s MODEL is the weakest link.  That the global coordinate CLOSES is
       VERIFIED on `420/420` rows, NOT PROVED; everything here lives inside that model.
  (ii) [NCYL-283]'s `p != 0 => R-R connection` is DERIVED and NOT SCORED (its own words).
  (iii) [NCYL-298]'s closed form `s_j == L_j + L_{j+1} (mod 2)` is derivation-checked
       `1288/1288`; [NCYL-300] and the theorem above are unconditional given the model.
  ⇒ **The honest sentence is: `p = 0` now has a PROOF whose only unproved inputs are the
  model's closure and one derived Gauss-Bonnet step -- not `p = 0` is a theorem full stop.**
  ⚠ And this says NOTHING about `[C.10]`, `[S2]`, `[R-DICH-CONV]`, `[M0-CHARGED]` or the
  F-leg: `p = 0` is the NECKLACE strand, and no gamma=1 DAG node depends on it.

HYPOTHESES.  ⚠ [OPS-041]: an arm that cannot fail is not evidence.
  H1  ⇒⇒ REVERSIBILITY, STEP BY STEP: `R T R == T^{-1}` at every state of every traced
      orbit, on REAL `Coord` rows AND on synthetic `Fraction` paths.  This is the ONE
      structural input of the proof and it can fail at any single step -- most plausibly at
      a fold whose reflection lands on an endpoint.
  H2  ⇒⇒ THE PALINDROME LEMMA, TESTED POSITIVELY AND NON-VACUOUSLY.  Run a POLE-TRANSPARENT
      walker (a pole is treated as an ordinary fold, so orbits do not stop) from arbitrary
      start states.  Steps 3-4 predict: **every return to `R(u_0)` has ODD step count and a
      POLE-fold exactly at the midpoint.**  Both halves can fail on any single return, and
      returns are plentiful -- this is the arm that tests the mechanism rather than the
      conclusion.
  H3  THE OBJECT, which is what makes Step 2 true: adjacent intervals share EXACTLY one
      endpoint, and `s_j/2` is interior to EXACTLY the LARGER of the two neighbours.
      ⚠ If `s_j/2` were interior to BOTH, no pole hit at that edge is possible at all and
      the theorem is vacuous there; if to NEITHER, the start state is invalid.  Scored.
  H4  THE DETECTOR IS LIVE.  A random rational family produces thousands of pole hits.
      ⚠⚠ `self_hits = 0` is now PREDICTED BY THE THEOREM and is therefore NOT EVIDENCE
      ([OPS-041]).  What IS evidence is `cross_hits > 0`: the detector demonstrably fires,
      so `0` self-hits is not a broken instrument.  A run with `hits = 0` is VACUOUS.
  H5  THE REAL-ROW ENDING CENSUS, and it is reported for SCOPE, not as support: how many
      prongs are resolved by measurement at all.  ⚠ The CAPPED rows are exactly the ones
      measurement CANNOT settle and the theorem can.
  H6  THE TERMINAL-STATE CLAIM ON INSTANCES WHERE A HIT HAPPENS: for every observed pole
      hit at edge `e`, the hitting node is the LARGER node of edge `e` (Step 2's forcing,
      checked where it bites rather than where it is assumed).

RESULTS (`all --qmax=31 --cap=200000 --trials=12000`, `81.5 s`).
  H1  ⇒⇒ `2648643/2648643` steps, `0` bad -- `817187` synthetic (`Fraction`) ⊕ `1831456`
      REAL (`Coord`, exact cyclotomic).  The proof's one structural input, and it held at
      every step of every orbit traced in either arithmetic.
  H2  ⇒⇒ `7201` returns to `R(u_0)`, and ALL `7201` have ODD step count with a POLE-fold
      exactly at the midpoint -- `0` bad.  **This is the arm that tests the MECHANISM, and
      each of the `7201` could have failed on either half.**  ⚠ On REAL rows H2 finds `0`
      returns, and that is expected, not a gap: a real return WOULD BE a self-hit, so the
      positive test necessarily lives in the synthetic family.
  H3  `4392/4392` adjacent pairs share exactly ONE endpoint; `4392/4392` poles `s_j/2`
      interior to exactly the LARGER neighbour (odd `Q = 5..31`, both classes).  Step 2
      rests on the second one.
  H4  `6822` pole hits, `6822` CROSS, **`0` SELF**.  ⚠⚠ The `0` is PREDICTED and is NOT
      evidence; the evidence is `cross = 6822 > 0`, i.e. the detector demonstrably fires.
  H5  REAL deep ending census, `cap = 200000`: `4306` ESC, `86` CAP, `0` Z, `0` R.
      ⚠⚠ **THE `86` CAPPED PRONGS ARE THE POINT: measurement cannot settle them at any
      cap (fold counts reach `77194` and climb -- [NCYL-294]), and the theorem settles
      them for free, because it never looks at the length of the orbit.**
  H6  `6822/6822`: every observed pole hit is made from the LARGER node of its edge.

  python3.13 probes/s454_self_hit.py [all|real|synth] [--qmax=N] [--trials=N]
"""
from __future__ import annotations

import json
import math
import random
import sys
import time
from fractions import Fraction

import s450_chain_maps as cm                                        # noqa: F401
import s451_quotient_path as qp
import s453_pole_module as pm

OUT = "data/s454_self_hit.json"


# ---------------------------------------------------------------- the tracer

def step(Lo, Hi, S, st, ops, transparent=False):
    """One application of `T`.  Returns `(next_state, kind, edge)`.

    `kind` in {'cross','fold','pole','ESC','Z'}.  With `transparent=True` a POLE is taken
    as an ordinary fold (it IS one -- `y |-> s_e - y` fixes `s_e/2`), so the orbit runs on;
    that is what H2 needs in order to observe returns at all.
    """
    i, d, y = st
    i2 = i + d
    if i2 < 0 or i2 >= len(Lo):
        return None, "ESC", None
    c0, c1 = ops.cmp(y, Lo[i2]), ops.cmp(y, Hi[i2])
    if c0 == 0 or c1 == 0:
        return None, "Z", None
    if c0 > 0 and c1 < 0:
        return (i2, d, y), "cross", None
    e = min(i, i2)
    is_pole = ops.eq(y, ops.half(S[e]))
    if is_pole and not transparent:
        return None, "pole", e
    return (i, -d, ops.sub(S[e], y)), ("pole" if is_pole else "fold"), e


def orbit(Lo, Hi, S, st0, ops, cap, transparent=False, stop_at=None):
    """Forward orbit from `st0`.  Returns (states, kinds, edges, end)."""
    sts, kinds, edges = [st0], [], []
    st = st0
    for _ in range(cap):
        nxt, kind, e = step(Lo, Hi, S, st, ops, transparent)
        kinds.append(kind)
        edges.append(e)
        if nxt is None:
            return sts, kinds, edges, kind
        st = nxt
        sts.append(st)
        if stop_at is not None and _eqstate(st, stop_at, ops):
            return sts, kinds, edges, "RETURN"
    return sts, kinds, edges, "CAP"


def _eqstate(a, b, ops):
    return a[0] == b[0] and a[1] == b[1] and ops.eq(a[2], b[2])


def pole_start(Lo, Hi, S, j, ops):
    """`u_0` for the prong of the pole on edge `j`: the LARGER node, heading AWAY."""
    y = ops.half(S[j])
    la, lb = ops.sub(Hi[j], Lo[j]), ops.sub(Hi[j + 1], Lo[j + 1])
    i = j if ops.cmp(la, lb) > 0 else j + 1
    return (i, -1 if i == j else +1, y)


def check_reversible(Lo, Hi, S, sts, ops):
    """H1 -- `R T R == T^{-1}` at every state of the orbit.  Returns (ok, bad)."""
    ok = bad = 0
    for k in range(1, len(sts)):
        i, d, y = sts[k]
        back, _kind, _e = step(Lo, Hi, S, (i, -d, y), ops, transparent=True)
        if back is not None and _eqstate((back[0], -back[1], back[2]), sts[k - 1], ops):
            ok += 1
        else:
            bad += 1
    return ok, bad


# ---------------------------------------------------------------- real rows

def real(qmax, cap, shallow=4000):
    res = {"rows": 0, "share_ok": 0, "share_bad": 0, "geom_ok": 0, "geom_bad": 0,
           "geom_ex": [], "rev_ok": 0, "rev_bad": 0, "rev_ex": [], "ends": {},
           "prongs": 0, "h2_returns": 0, "h2_ok": 0, "h2_bad": 0, "h2_ex": [],
           "hits": 0, "self_hits": 0, "h6_ok": 0, "h6_bad": 0, "deep_ends": {}}
    for Q in range(5, qmax + 1, 2):
        for P in range(1, Q):
            if math.gcd(P, Q) != 1:
                continue
            for eps in (0, 1):
                ch, lo, hi, closes = qp.build(P, Q, eps)
                if not closes:
                    continue
                s = pm.fold_sums(ch, lo, hi)
                Lo, Hi, S = pm.path_view(ch, lo, hi, s)
                ops = pm.CoordOps(ch)
                m = ch.m_idx
                res["rows"] += 1
                # -- H3(a): adjacent intervals share exactly ONE endpoint
                for i in range(m):
                    sl, sh = ops.eq(Lo[i], Lo[i + 1]), ops.eq(Hi[i], Hi[i + 1])
                    res["share_ok" if sl != sh else "share_bad"] += 1
                # -- H3(b): s_j/2 interior to exactly the LARGER neighbour
                for j in range(m):
                    y = ops.half(S[j])
                    inside = [k for k in (j, j + 1)
                              if ops.cmp(y, Lo[k]) > 0 and ops.cmp(y, Hi[k]) < 0]
                    la, lb = ops.sub(Hi[j], Lo[j]), ops.sub(Hi[j + 1], Lo[j + 1])
                    larger = j if ops.cmp(la, lb) > 0 else j + 1
                    good = inside == [larger]
                    res["geom_ok" if good else "geom_bad"] += 1
                    if not good and len(res["geom_ex"]) < 5:
                        res["geom_ex"].append([P, Q, eps, j, inside, larger])
                # -- H1/H5/H6: trace each prong
                for j in range(m):
                    u0 = pole_start(Lo, Hi, S, j, ops)
                    sts, kinds, edges, end = orbit(Lo, Hi, S, u0, ops, shallow)
                    res["prongs"] += 1
                    res["ends"][end] = res["ends"].get(end, 0) + 1
                    o, b = check_reversible(Lo, Hi, S, sts, ops)
                    res["rev_ok"] += o
                    res["rev_bad"] += b
                    if b and len(res["rev_ex"]) < 5:
                        res["rev_ex"].append([P, Q, eps, j])
                    if end == "pole":
                        res["hits"] += 1
                        e = edges[-1]
                        if e == j:
                            res["self_hits"] += 1
                        node = sts[-1][0]
                        la, lb = ops.sub(Hi[e], Lo[e]), ops.sub(Hi[e + 1], Lo[e + 1])
                        larger = e if ops.cmp(la, lb) > 0 else e + 1
                        res["h6_ok" if node == larger else "h6_bad"] += 1
                # -- H5 (deep): the ENDING census at the full cap, via the reference
                # implementation `s453_pole_module.walk`, which stores no states.  Reported
                # for SCOPE: the CAPPED prongs are the ones measurement cannot settle.
                for j in range(m):
                    end, _folds = pm.walk(Lo, Hi, S, j, ops, cap=cap)
                    lab = end[0] if end[0] != "R" else ("R-self" if end[1] == j else "R-x")
                    res["deep_ends"][lab] = res["deep_ends"].get(lab, 0) + 1
                # -- H2: pole-transparent returns to R(u_0), from every pole start
                for j in range(m):
                    u0 = pole_start(Lo, Hi, S, j, ops)
                    tgt = (u0[0], -u0[1], u0[2])
                    sts, kinds, edges, end = orbit(Lo, Hi, S, u0, ops, min(cap, 4000),
                                                   transparent=True, stop_at=tgt)
                    if end != "RETURN":
                        continue
                    N = len(sts) - 1
                    res["h2_returns"] += 1
                    good = (N % 2 == 1) and kinds[(N - 1) // 2] == "pole"
                    res["h2_ok" if good else "h2_bad"] += 1
                    if not good and len(res["h2_ex"]) < 5:
                        res["h2_ex"].append([P, Q, eps, j, N, kinds[(N - 1) // 2]])
    return res


# ---------------------------------------------------------------- synthetic

def rand_path(rng, n):
    """A random rational path: `n` nested intervals sharing alternating endpoints.
    Same construction as `s453_pole_module.synth`'s `rational` family."""
    L = [Fraction(rng.randint(1, 6)) for _ in range(n)]
    sig = [rng.choice((+1, -1)) for _ in range(n - 1)]
    Lo, Hi = [Fraction(0)], [L[0]]
    for j in range(n - 1):
        nlo = Lo[j] if sig[j] > 0 else (Lo[j] + L[j]) - L[j + 1]
        Lo.append(nlo)
        Hi.append(nlo + L[j + 1])
    S = [(Hi[j] + Hi[j + 1]) if sig[j] > 0 else (Lo[j] + Lo[j + 1]) for j in range(n - 1)]
    if any(Lo[i] >= Hi[i] for i in range(n)):
        return None
    return Lo, Hi, S


def synth(trials, seed=20260901):
    """H1/H2/H4/H6 on arbitrary (non-arithmetic) interval data -- the family in which the
    theorem's indifference to `Q`, `P` and the lengths is actually exercised."""
    rng = random.Random(seed)
    ops = pm.FracOps
    res = {"paths": 0, "prongs": 0, "ends": {}, "hits": 0, "self_hits": 0,
           "cross_hits": 0, "rev_ok": 0, "rev_bad": 0, "rev_ex": [],
           "h2_returns": 0, "h2_ok": 0, "h2_bad": 0, "h2_ex": [],
           "h6_ok": 0, "h6_bad": 0, "steps": 0}
    for _ in range(trials):
        n = rng.randint(3, 8)
        got = rand_path(rng, n)
        if got is None:
            continue
        Lo, Hi, S = got
        res["paths"] += 1
        for j in range(n - 1):
            u0 = pole_start(Lo, Hi, S, j, ops)
            if not (Lo[u0[0]] < u0[2] < Hi[u0[0]]):
                continue
            sts, kinds, edges, end = orbit(Lo, Hi, S, u0, ops, 20000)
            res["prongs"] += 1
            res["steps"] += len(sts)
            res["ends"][end] = res["ends"].get(end, 0) + 1
            o, b = check_reversible(Lo, Hi, S, sts, ops)
            res["rev_ok"] += o
            res["rev_bad"] += b
            if b and len(res["rev_ex"]) < 5:
                res["rev_ex"].append([n, j])
            if end == "pole":
                e = edges[-1]
                res["hits"] += 1
                res["self_hits" if e == j else "cross_hits"] += 1
                node = sts[-1][0]
                larger = e if (Hi[e] - Lo[e]) > (Hi[e + 1] - Lo[e + 1]) else e + 1
                res["h6_ok" if node == larger else "h6_bad"] += 1
        # -- H2 from ARBITRARY start states (not just poles): the mechanism test
        for _ in range(6):
            i0 = rng.randrange(n)
            if rng.random() < 0.35 and n > 1:
                e = rng.randrange(n - 1)
                y0 = S[e] / 2
                if not (Lo[i0] < y0 < Hi[i0]):
                    continue
            else:
                y0 = Lo[i0] + Fraction(rng.randint(1, 11), 12) * (Hi[i0] - Lo[i0])
            d0 = rng.choice((+1, -1))
            u0 = (i0, d0, y0)
            tgt = (i0, -d0, y0)
            sts, kinds, edges, end = orbit(Lo, Hi, S, u0, ops, 400,
                                           transparent=True, stop_at=tgt)
            o, b = check_reversible(Lo, Hi, S, sts, ops)
            res["rev_ok"] += o
            res["rev_bad"] += b
            if end != "RETURN":
                continue
            N = len(sts) - 1
            res["h2_returns"] += 1
            good = (N % 2 == 1) and kinds[(N - 1) // 2] == "pole"
            res["h2_ok" if good else "h2_bad"] += 1
            if not good and len(res["h2_ex"]) < 5:
                res["h2_ex"].append([n, i0, d0, N, kinds[(N - 1) // 2]])
    return res


# ---------------------------------------------------------------- main

def main():
    mode, qmax, trials, cap = "all", 31, 6000, 200_000
    for a in sys.argv[1:]:
        if a.startswith("--qmax="):
            qmax = int(a.split("=")[1])
        elif a.startswith("--trials="):
            trials = int(a.split("=")[1])
        elif a.startswith("--cap="):
            cap = int(a.split("=")[1])
        elif not a.startswith("-"):
            mode = a
    t0 = time.time()
    out = {"mode": mode, "qmax": qmax, "trials": trials, "cap": cap}
    if mode in ("all", "synth"):
        out["synth"] = synth(trials)
        r = out["synth"]
        print(f"SYNTH paths={r['paths']} prongs={r['prongs']} steps={r['steps']} "
              f"ends={r['ends']}")
        print(f"  H1 reversibility R T R == T^-1: ok={r['rev_ok']} BAD={r['rev_bad']} "
              f"{r['rev_ex'][:3]}")
        print(f"  H2 returns to R(u_0)={r['h2_returns']}  N ODD and midpoint a POLE: "
              f"ok={r['h2_ok']} BAD={r['h2_bad']} {r['h2_ex'][:3]}   "
              f"(⚠ VACUOUS if returns == 0)")
        print(f"  H4 pole hits={r['hits']} cross={r['cross_hits']} "
              f"SELF={r['self_hits']}   (⚠ cross>0 is the evidence; SELF=0 is PREDICTED)")
        print(f"  H6 hitting node is the LARGER: ok={r['h6_ok']} BAD={r['h6_bad']}")
    if mode in ("all", "real"):
        out["real"] = real(qmax, cap)
        r = out["real"]
        print(f"REAL rows={r['rows']} prongs={r['prongs']} ends(shallow)={r['ends']}")
        print(f"  H5 DEEP ending census (cap={cap}, reference `pm.walk`): {r['deep_ends']}"
              f"   ⚠ CAP = unresolved by measurement; the theorem is what settles those")
        print(f"  H3a adjacent share exactly ONE endpoint: ok={r['share_ok']} "
              f"BAD={r['share_bad']}")
        print(f"  H3b s_j/2 interior to exactly the LARGER: ok={r['geom_ok']} "
              f"BAD={r['geom_bad']} {r['geom_ex'][:3]}")
        print(f"  H1 reversibility: ok={r['rev_ok']} BAD={r['rev_bad']} {r['rev_ex'][:3]}")
        print(f"  H2 returns={r['h2_returns']} ok={r['h2_ok']} BAD={r['h2_bad']}")
        print(f"  H4 pole hits={r['hits']} SELF={r['self_hits']}   "
              f"H6 ok={r['h6_ok']} BAD={r['h6_bad']}")
    out["secs"] = round(time.time() - t0, 1)
    with open(OUT, "w") as fh:
        json.dump(out, fh)
    print(f"-> {OUT}  ({out['secs']} s)")


if __name__ == "__main__":
    main()
