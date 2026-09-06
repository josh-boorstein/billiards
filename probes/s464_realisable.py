#!/usr/bin/env python3
"""s464_realisable.py -- PRE-REGISTERED.  Queue item (s462-b): REALISABILITY.  With the
arithmetic side of (G0b-rel) fully understood -- [NCYL-318] gives the relation lattice
with NAMED GENERATORS -- the residual is *which relation vectors are attainable as the
arrival vector `a` of an actual fold word*.  `Q = 9, 21, 27, 33, 39` carry the SAME
`p = 3` generator family as `Q = 15` and produce NO witness, so defect `> 0` is necessary
and not sufficient, and what is missing is a statement about the WALK.

WHAT THIS IS.  NOT a proof of (G0b), and NOT a claim about `Q = 15`.  Two things:

  (A) ⇒⇒ THE EXACT `+-1` ENUMERATION THAT [NCYL-308] (3) RELAXED.  Its H6 restored the
      alternating-sign constraint by ONE augmented lattice coordinate and moved the corner
      hit rate `85.75% -> 77.60%`; the note itself flags the residual gap ("the Z-span is a
      superset of the +-1-coefficient set", s455 H8).  This file does not augment: it
      ENUMERATES the +-1 words, exactly, as a nondeterministic relaxation of the walk in
      which the interval-ORDERING tests are dropped and the ARITHMETIC is untouched.
      ⇒ that places it strictly between two things already measured: [NCYL-308] (3)'s
      lattice surrogate (`77.6%` blind) and [NCYL-308] (4)'s synthetic control (arithmetic
      REMOVED -- `Z` happens freely, `500` of `3285`).  So the score has two known
      neighbours and cannot be read as an isolated number.

  (B) ⇒⇒ AND A COORDINATE DEFECT IN THE `rel`/`form` SPLIT ITSELF, which has to be fixed
      before (A) can be scored.  [NCYL-317] (4) classifies an arrival by `a != 0` in the
      RAW basis `(ell_0 .. ell_{Q-1})`.  But that basis is NOT the basis the rank law is
      about: `c_t = c_0 + 2Pt (mod 2Q)` runs over every residue of one parity, so on the
      LONE class exactly one `t` has `ell_t = 0` (the horizontal interface) and the other
      `Q-1` fall into `(Q-1)/2` pairs `c ~ 2Q-c` of EQUAL measure.  ⇒ raw-coordinate
      relations `e_t` (zero measure) and `e_t - e_t'` (equal measures) exist at EVERY odd
      `Q`, prime or composite.  H0 measures this; it is arithmetic, not a conjecture.
      ⚠⚠ SO THE HONEST SPLIT IS THREE-WAY, on the REDUCED vector `a_red` (sum `a_t` over
      each distinct nonzero measure, drop the zero one):
          form  `a_raw == 0`   -- a formal integer identity
          triv  `a_raw != 0`, `a_red == 0` -- vanishes by the zero interface / the +- pairing
                                             ALONE, no cyclotomic relation consumed
          rel   `a_red != 0`   -- a genuine relation among DISTINCT measures consumed
      and [NCYL-318]'s theorem bars `rel` on the LONE class at prime `Q`.  It does NOT bar
      `triv` there, and `triv` is available at every `Q`.

PRIOR ART: `rulings.py --grep` on 'realisab' -> [WFLOOR-087], [NCYL-177], [CHARGE-005],
[OPS-019], [NCYL-300], none about this; 'fold word' -> [NCYL-317] only (which coined it);
'arrival vector', 'attainable' -> nothing.  Re-run of s463's gate, and it still holds.
⚠⚠ THE GATE THAT MATTERS IS [NCYL-308], whose four priced-dead routes this must not be:
(1) the palindrome transfer and (4) the topological/`chi` route are ROW-INDEPENDENT and
dead a priori ([OPS-242]) -- this is row-dependent, it consumes the row's own `S_e` and
its own path.  (2) the plain `Lambda_S` test and (3) the coefficient-sum augmentation ask
Z-SPAN membership; this asks whether a `+-1` WORD exists, which is the subset (3) itself
names as the thing it relaxes.  ⇒ it is route (3) done exactly rather than surrogately,
and its own pre-registered failure mode is that the exact version buys nothing over the
surrogate -- in which case route (3) is closed rather than merely priced.
⚠ [NCYL-315] bars renormalization / self-similarity: does not bite, nothing below reads a
CF digit or induces a return map.  ⚠ [NCYL-296] bars routing through CF depth: same.
⚠ user-directed s460 bars EXTENDING THE RANGE as a way to characterise: does not bite --
the `Q` set here is FIXED by the question (`3 | Q`, `Q <= 39`), and nothing below proposes
a bigger box as an answer.

HYPOTHESES (pre-registered).
  H0  ⇒⇒ THE COORDINATE STRUCTURE, and it is a MEASUREMENT of an arithmetic fact, not an
      arm: per class, the number of `t` with `ell_t = 0` and the multiset of multiplicities
      of the distinct nonzero measures.  PREDICTION: LONE = `1` zero + `(Q-1)/2` values of
      multiplicity `2`; PAIRED = `0` zeros + `(Q+1)/2` values, one of multiplicity `1` and
      the rest `2`.  ⇒ raw-coordinate relations exist everywhere, so [NCYL-317] (4)'s
      raw-basis `a != 0` is the wrong classifier and (B)'s three-way split replaces it.
      ⚠ If the prediction FAILS the whole of (B) is wrong and H1's re-classification is
      meaningless; it is scored per class-row and can fail on either count.
  H1  ⇒⇒ THE GATE, and a PASS IS PREDICTED -- it is a bug-check, never evidence
      ([OPS-242]/[OPS-041]).  Two parts.  (a) The relaxation must ADMIT the real witness at
      `7/15 eps1`, at depth `<= 5`, with the recorded `a` of [NCYL-317] (2).  The real walk
      IS one of the enumerated words, so a miss is an enumerator bug.  (b) That `a` must
      re-classify as `rel` under (B) -- `a_red != 0`.  ⚠ (b) CAN FAIL, and if it does then
      [NCYL-317] (2)'s "the witness consumes an arithmetic coincidence" is itself wrong and
      the mechanism split collapses; this is the one place (B) could bite backwards.
  H2  ⇒⇒ THE DISCRIMINATION ARM, and BOTH outcomes are informative.  Run the relaxation on
      every class row at `Q = 9, 15, 21, 27, 33, 39` and classify every admitted arrival.
        no `rel` arrival at `Q != 15`  ⇒ the obstruction is COMBINATORIAL: the relation
             exists, and no `+-1` fold word reaches it.  That is a criterion, and the next
             instrument is about words.
        `rel` arrivals at `Q != 15`    ⇒ the obstruction is DYNAMICAL: the word exists and
             the ordering tests refuse it.  Then no arithmetic instrument can close
             (G0b-rel) and the next instrument must be about the interval tests.
      ⚠ It is a DEPTH-CAPPED search: a null at depth `D` is not "no word exists" ([OPS-041]),
      and the depth reached and the truncation flag are reported with every null.
  H3  ⇒⇒ THE PRICING, against the two known neighbours.  Admission rate over the corner
      census at odd `Q <= QCEN`, split `form`/`triv`/`rel` and by class/primality.  ⚠ The
      comparison is to [NCYL-308] H6's `77.60%` and H4's synthetic `Z`-happens-freely, NOT
      to zero: this instrument is a RELAXATION, so a high rate prices it as blind and a low
      one makes the `+-1` structure the content.  Either way it retires route (3).
  H4  ⇒⇒ THE CONTROL THAT CAN REFUTE SOMETHING.  On the LONE class at PRIME `Q`,
      [NCYL-318] PROVES the relation space is `0`, so the relaxation must return `0` `rel`
      arrivals there -- over the whole prime-`Q` LONE census, not a sample.  A single `rel`
      arrival on such a row means the enumerator, the reduction map or [NCYL-318] is wrong.
      ⚠ `triv` and `form` arrivals there are NOT a failure -- they are the point of (B),
      and their count is the measurement that qualifies [NCYL-317] (4).

GUARDS.
  G1  a row whose global coordinate does not close (`s451.build`) is SKIPPED, not scored.
  G2  a TIE edge / degenerate node carries no corner prong; SKIPPED, never a pass.
  G3  nothing here writes to any store, and NO existing probe is modified -- `sym_build`,
      `sym_path`, `check_sym`, `sym_value` are IMPORTED from `s462_zz_relation`.
  G4  ⚠ [OPS-228]: numpy `int64` is cast to Python `int` before sympy / any gcd path.
  G5  the symbolic layer is re-checked against the ring layer on every row before any `a`
      is reported (s462's derivation guard, re-run here rather than inherited).
  G6  ⇒⇒ THE RELAXATION MUST BE SOUND -- every real walk must be one of the enumerated
      words.  Scored directly: replay `sh.orbit` from each corner start and assert every
      state it visits is in the enumerated reachable set (to the shared depth).  A
      relaxation that loses the real walk cannot support a NULL result, which is exactly
      what H2 may return.
  G7  the `y in J_0` prune is used only after NESTEDNESS is verified on the row
      (`J_{i+1} subset J_i` or the reverse, exactly); on a row where it fails the prune is
      dropped rather than the row.

Run: PYTHONPATH=.:engine:archive/scripts_2026-07:archive/scripts_2026-08:probes \\
     .venv/bin/python3.13 probes/s464_realisable.py [all|h0|h1|h2|h3|h4] \\
     [--depth=N] [--states=N] [--qcen=N]
"""
from __future__ import annotations

import json
import math
import sys
import time

import numpy as np

import s450_chain_maps as cm
import s451_quotient_path as qp
import s453_pole_module as pm
import s454_self_hit as sh
import s456_neckgb as nb
from s462_zz_relation import (check_sym, isprime, phi, sym_build, sym_path,
                              sym_value)

OUT = "data/s464_realisable.json"
WITNESS = [(7, 15, 1), (8, 15, 0)]
TARGET_Q = [9, 15, 21, 27, 33, 39]          # the (s462-b) set: 3 | Q, one witness among them

DEPTH = 20          # the settled rung; the 14 -> 20 trajectory is the point (see H2)
MAXSTATES = 40_000
QCEN = 21


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


# ------------------------------------------------------------------ the reduction map

def reduce_map(ch):
    """`red[t]` = index of `ell_t` among the DISTINCT NONZERO measures, or `None`.

    Keyed by the exact ring vector, so two `t` are identified iff their measures are
    EQUAL as elements of `Z[zeta_{4Q}]` -- never by float and never by the `c` label
    (`c ~ 2Q-c` is the reason they coincide, but the identification is made on the value).
    """
    key, red, order = {}, [], []
    for t in range(ch.Q):
        v = ch.ell[t].v
        if not v.any():
            red.append(None)
            continue
        b = v.tobytes()
        if b not in key:
            key[b] = len(order)
            order.append(t)
        red.append(key[b])
    return red, order


def reduce_vec(a, red, ndist):
    out = [0] * ndist
    for t, c in enumerate(a):
        if c and red[t] is not None:
            out[red[t]] += int(c)
    return out


def classify(a, red, ndist):
    """The three-way split of (B).  `a` is the RAW arrival vector."""
    if not any(a):
        return "form", [0] * ndist
    ar = reduce_vec(a, red, ndist)
    return ("triv" if not any(ar) else "rel"), ar


# ------------------------------------------------------------------ the relaxation

def nested_ok(Lo, Hi, ops):
    """G7: the path intervals are nested, so `y in J_0` (the widest) is a sound prune."""
    widest = max(range(len(Lo)), key=lambda i: float(ops.sub(Hi[i], Lo[i]).f))
    for i in range(len(Lo)):
        if ops.cmp(Lo[widest], Lo[i]) > 0 or ops.cmp(Hi[widest], Hi[i]) < 0:
            return None
    return widest


def relaxed_reach(Lo, Hi, S, ops, st0, sym0, sLo, sHi, sS,
                  depth=DEPTH, maxstates=MAXSTATES, widest=None):
    """BFS over the NONDETERMINISTIC relaxation of the walk.

    At every state BOTH moves are offered -- CROSS to `i+d` and FOLD at `e = min(i,i+d)` --
    regardless of the interval-ordering test that would decide between them in
    `s454_self_hit.step`.  Everything else is verbatim: `Z` when `y` is an endpoint of the
    interval being entered, `pole` when a fold lands on `s_e/2`, `ESC` off the end.  ⇒ the
    real orbit is ONE of the enumerated words (G6 scores this), and the enumeration is
    exactly the set of `+-1` alternating fold words on the row's own path.

    Returns `(zhits, nstates, depth_reached, truncated, closed)`; each hit carries the RAW
    arrival vector `a = ysym - tgtsym` of the first word reaching it.

    ⇒⇒ `closed` IS THE DIFFERENCE BETWEEN A CAPPED NULL AND AN EXHAUSTIVE ONE, and it is
    the only reason a null here is quotable.  States are deduped on `(node, dir, y)`, so an
    EMPTY frontier means the reachable set is closed under both moves -- no word of ANY
    length reaches a `Z`, not merely no word of length `<= depth`.  A null with
    `closed = False` is depth-limited and says nothing ([OPS-041]).
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
                continue                                            # ESC, terminal
            c0, c1 = ops.cmp(y, Lo[i2]), ops.cmp(y, Hi[i2])
            if c0 == 0 or c1 == 0:                                  # Z, terminal
                tgt = sLo[i2] if c0 == 0 else sHi[i2]
                zhits.append({"steps": dstep + 1, "node": i2,
                              "at": "lo" if c0 == 0 else "hi",
                              "a": [int(x - t) for x, t in zip(ysym, tgt)]})
                continue
            cands = [((i2, d, y), ysym)]                            # CROSS, always offered
            e = min(i, i2)
            if not ops.eq(y, ops.half(S[e])):                       # else: pole, terminal
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


def start_sym(Lo, Hi, sLo, sHi, ops, y):
    for idx in range(len(Lo)):
        if ops.eq(y, Lo[idx]):
            return list(sLo[idx])
        if ops.eq(y, Hi[idx]):
            return list(sHi[idx])
    return None


def scan_row(P, Q, eps, depth=DEPTH, maxstates=MAXSTATES, prune=True):
    """Every corner prong of one class row, through the relaxation, classified."""
    c = _chain(P, Q, eps)
    if c is None:
        return None
    ch, lo, hi, Lo, Hi, S, ops = c
    LO, HI = sym_build(ch)
    good, where = check_sym(ch, lo, hi, LO, HI)                     # G5
    if not good:
        return {"guard_bad": list(where)}
    sLo, sHi, sS = sym_path(ch, LO, HI, lo, hi)
    red, order = reduce_map(ch)
    ndist = len(order)
    widest = nested_ok(Lo, Hi, ops) if prune else None              # G7
    out = {"P": P, "Q": Q, "eps": eps, "ndist": ndist,
           "nested": widest is not None, "prongs": 0, "skipped": 0,
           "form": 0, "triv": 0, "rel": 0, "states": 0,
           "truncated": 0, "depth_reached": 0, "closed": 0,
           # ⇒ the denominator [NCYL-308] H6's `77.60%` uses: (corner start, target)
           #   PAIRS, a target being an endpoint of a node.  Quoting an arrival COUNT
           #   against a RATE would be comparing two different things.
           "pairs": 0, "targets": 0, "prongs_any": 0, "prongs_rel": 0, "hits": []}
    for j in range(len(S)):
        cs = nb.corner_start(Lo, Hi, S, j, ops)
        if cs is None:
            out["skipped"] += 1                                     # G2
            continue
        st0 = cs[0]
        sym0 = start_sym(Lo, Hi, sLo, sHi, ops, st0[2])
        if sym0 is None:
            out["skipped"] += 1
            continue
        out["prongs"] += 1
        out["pairs"] += 2 * len(Lo)
        zh, ns, dr, tr, cl = relaxed_reach(Lo, Hi, S, ops, st0, sym0, sLo, sHi, sS,
                                           depth, maxstates, widest)
        out["states"] += ns
        out["truncated"] += int(tr)
        out["closed"] += int(cl)
        out["depth_reached"] = max(out["depth_reached"], dr)
        seen_tgt = set()
        for h in zh:
            seen_tgt.add((h["node"], h["at"]))
            kind, ar = classify(h["a"], red, ndist)
            # the vanishing is a THEOREM of the construction, but it is free to check
            assert not sym_value(ch, h["a"]).any(), (P, Q, eps, j)
            out[kind] += 1
            if kind == "rel":
                out["hits"].append({"j": j, "steps": h["steps"], "kind": kind,
                                    "a": h["a"], "a_red": ar})
        out["targets"] += len(seen_tgt)
        out["prongs_any"] += int(bool(zh))
        out["prongs_rel"] += int(any(classify(h["a"], red, ndist)[0] == "rel"
                                     for h in zh))
    return out


# ------------------------------------------------------------------ H0

def h0(qmax=45):
    """THE COORDINATE STRUCTURE: zeros and multiplicities, per class."""
    res = {"checked": 0, "lone_ok": 0, "paired_ok": 0, "bad": [], "sample": []}
    for Q in range(5, qmax + 1, 2):
        for lone in (True, False):
            ch = cm.Chain(1, Q, 1 if lone else 0)
            mult = {}
            nzero = 0
            for t in range(Q):
                v = ch.ell[t].v
                if not v.any():
                    nzero += 1
                else:
                    b = v.tobytes()
                    mult[b] = mult.get(b, 0) + 1
            ms = sorted(mult.values())
            res["checked"] += 1
            if lone:
                ok = (nzero == 1 and len(ms) == (Q - 1) // 2 and set(ms) == {2})
                res["lone_ok"] += ok
            else:
                ok = (nzero == 0 and len(ms) == (Q + 1) // 2
                      and sorted(ms) == [1] + [2] * ((Q - 1) // 2))
                res["paired_ok"] += ok
            if not ok:
                res["bad"].append([Q, "lone" if lone else "paired", nzero, ms[:8]])
            if Q in (9, 15):
                res["sample"].append({"Q": Q, "class": "lone" if lone else "paired",
                                      "zeros": nzero, "ndist": len(ms),
                                      "mult": ms})
    res["raw_relations_exist_everywhere"] = (res["lone_ok"] + res["paired_ok"]
                                             == res["checked"])
    return res


# ------------------------------------------------------------------ H1 (the gate) + G6

def h1(depth=DEPTH, maxstates=MAXSTATES):
    """THE GATE: the relaxation must admit the real witness, and it must re-classify."""
    res = {"cases": [], "soundness": {}, "admits": None, "witness_kind": None}
    for (P, Q, eps) in WITNESS:
        r = scan_row(P, Q, eps, depth, maxstates)
        rel = [h for h in r["hits"] if h["kind"] == "rel"]
        res["cases"].append({"P": P, "Q": Q, "eps": eps,
                             "form": r["form"], "triv": r["triv"], "rel": r["rel"],
                             "min_rel_steps": min([h["steps"] for h in rel], default=None),
                             "rel_hits": rel[:4]})
    a15 = None
    for h in res["cases"][0]["rel_hits"]:
        if h["steps"] <= 5:
            a15 = h
            break
    res["admits"] = a15 is not None
    res["witness_kind"] = a15["kind"] if a15 else None
    res["witness_a"] = a15["a"] if a15 else None
    res["witness_a_red"] = a15["a_red"] if a15 else None
    # [NCYL-317] (2)'s recorded vector, for identification.  ⚠ UP TO SIGN: (1) of that
    # ruling says the connection's two ENDS carry `a` and `-a`, so a sign-strict test
    # would fail on the other end and mean nothing.
    ref = [0, 1, -1, 2, -1, 2, 0, 0]
    res["matches_NCYL317_a"] = bool(
        a15 and (list(a15["a"][:8]) == ref
                 or list(a15["a"][:8]) == [-x for x in ref]))

    # G6: soundness -- every state of the REAL orbit is in the enumerated set.
    # ⚠⚠ SCORED ONLY ON PRONGS WHOSE BFS DID NOT HIT THE STATE CAP.  A truncated BFS can
    # drop the real walk itself, so on those prongs this is not a relaxation at all -- and
    # its null is worth even less than a depth-capped one.  One blended number hides that.
    ok = bad = 0
    trunc_prongs = 0
    for (P, Q, eps) in WITNESS + [(1, 9, 0), (2, 21, 1), (4, 27, 0), (5, 33, 1)]:
        c = _chain(P, Q, eps)
        if c is None:
            continue
        ch, lo, hi, Lo, Hi, S, ops = c
        LO, HI = sym_build(ch)
        sLo, sHi, sS = sym_path(ch, LO, HI, lo, hi)
        widest = nested_ok(Lo, Hi, ops)
        for j in range(len(S)):
            cs = nb.corner_start(Lo, Hi, S, j, ops)
            if cs is None:
                continue
            sym0 = start_sym(Lo, Hi, sLo, sHi, ops, cs[0][2])
            if sym0 is None:
                continue
            sts, kinds, _e, _end = sh.orbit(Lo, Hi, S, cs[0], ops, depth)
            # re-run the BFS keeping the visited KEY set, then test containment
            keys, ktr = _reach_keys(Lo, Hi, S, ops, cs[0], depth, maxstates, widest)
            if ktr:
                trunc_prongs += 1
                continue
            for st in sts[:depth + 1]:
                k = (st[0], st[1], st[2].v.tobytes())
                if k in keys:
                    ok += 1
                else:
                    bad += 1
    res["soundness"] = {"states_covered": ok, "states_missed": bad,
                        "sound_on_untruncated": bad == 0,
                        "prongs_state_capped": trunc_prongs}
    return res


def _reach_keys(Lo, Hi, S, ops, st0, depth, maxstates, widest):
    """The reachable KEY set of the relaxation (G6's half of `relaxed_reach`)."""
    key0 = (st0[0], st0[1], st0[2].v.tobytes())
    seen = {key0}
    frontier = [st0]
    for _ in range(depth):
        nxt = []
        for (i, d, y) in frontier:
            i2 = i + d
            if i2 < 0 or i2 >= len(Lo):
                continue
            c0, c1 = ops.cmp(y, Lo[i2]), ops.cmp(y, Hi[i2])
            if c0 == 0 or c1 == 0:
                continue
            cands = [(i2, d, y)]
            e = min(i, i2)
            if not ops.eq(y, ops.half(S[e])):
                yf = ops.sub(S[e], y)
                if widest is None or (ops.cmp(yf, Lo[widest]) >= 0
                                      and ops.cmp(yf, Hi[widest]) <= 0):
                    cands.append((i, -d, yf))
            for st in cands:
                k = (st[0], st[1], st[2].v.tobytes())
                if k in seen:
                    continue
                seen.add(k)
                nxt.append(st)
            if len(seen) > maxstates:
                return seen, True
        if not nxt:
            break
        frontier = nxt
    return seen, False


# ------------------------------------------------------------------ H2

def h2(depth=DEPTH, maxstates=MAXSTATES, qs=None):
    """THE DISCRIMINATION ARM: the `3 | Q` set, `Q = 15` among them."""
    out = {"by_Q": {}, "depth": depth, "maxstates": maxstates, "rel_rows": []}
    for Q in (qs or TARGET_Q):
        agg = {"rows": 0, "prongs": 0, "form": 0, "triv": 0, "rel": 0,
               "rows_with_rel": 0, "truncated_prongs": 0, "states": 0,
               "closed_prongs": 0}
        for P in range(1, Q):
            if math.gcd(P, Q) != 1:
                continue
            for eps in (0, 1):
                r = scan_row(P, Q, eps, depth, maxstates)
                if r is None or "guard_bad" in r:
                    continue
                agg["rows"] += 1
                for k in ("prongs", "form", "triv", "rel", "states"):
                    agg[k] += r[k]
                agg["truncated_prongs"] += r["truncated"]
                agg["closed_prongs"] += r["closed"]
                if r["rel"]:
                    agg["rows_with_rel"] += 1
                    out["rel_rows"].append({"P": P, "Q": Q, "eps": eps,
                                            "rel": r["rel"],
                                            "min_steps": min(h["steps"]
                                                             for h in r["hits"]),
                                            "hits": r["hits"][:3]})
        out["by_Q"][str(Q)] = agg
    q15 = out["by_Q"].get("15", {})
    others = [v for k, v in out["by_Q"].items() if k != "15"]
    out["only_15_has_rel"] = bool(q15.get("rel", 0) > 0
                                  and all(v["rel"] == 0 for v in others))
    out["obstruction"] = ("COMBINATORIAL" if out["only_15_has_rel"]
                          else "DYNAMICAL-or-DEEPER")
    return out


# ------------------------------------------------------------------ H3

def h3(qcen=QCEN, depth=DEPTH, maxstates=MAXSTATES):
    """THE PRICING, against [NCYL-308] H6's `77.60%` and H4's synthetic control."""
    out = {"qcen": qcen, "depth": depth, "rows": 0, "prongs": 0,
           "prongs_admitting_any": 0, "prongs_admitting_rel": 0,
           "form": 0, "triv": 0, "rel": 0, "truncated_prongs": 0,
           "closed_prongs": 0,
           "by_class": {"lone": {"rows": 0, "rel": 0, "triv": 0, "form": 0},
                        "paired": {"rows": 0, "rel": 0, "triv": 0, "form": 0}}}
    for P, Q, eps in _rows(qcen):
        r = scan_row(P, Q, eps, depth, maxstates)
        if r is None or "guard_bad" in r:
            continue
        lone = (eps == P % 2)
        cl = "lone" if lone else "paired"
        out["rows"] += 1
        out["prongs"] += r["prongs"]
        out["truncated_prongs"] += r["truncated"]
        out["closed_prongs"] += r["closed"]
        for k in ("pairs", "targets"):
            out[k] = out.get(k, 0) + r[k]
        out["prongs_admitting_any"] += r["prongs_any"]
        out["prongs_admitting_rel"] += r["prongs_rel"]
        for k in ("form", "triv", "rel"):
            out[k] += r[k]
            out["by_class"][cl][k] += r[k]
        out["by_class"][cl]["rows"] += 1
    tot = out["form"] + out["triv"] + out["rel"]
    out["arrivals_total"] = tot
    out["rel_share"] = round(out["rel"] / tot, 4) if tot else None
    out["pair_admission_rate"] = (round(out["targets"] / out["pairs"], 5)
                                  if out.get("pairs") else None)
    out["prong_admission_rate"] = (round(out["prongs_admitting_any"] / out["prongs"], 5)
                                   if out["prongs"] else None)
    return out


# ------------------------------------------------------------------ H4

def h4(qmax=31, depth=DEPTH, maxstates=MAXSTATES):
    """THE CONTROL: no `rel` arrival on the LONE class at PRIME `Q` ([NCYL-318])."""
    out = {"qmax": qmax, "rows": 0, "rel": 0, "triv": 0, "form": 0,
           "violations": [], "by_Q": {}}
    for Q in range(5, qmax + 1, 2):
        if not isprime(Q):
            continue
        agg = {"rows": 0, "rel": 0, "triv": 0, "form": 0, "prongs": 0,
               "closed_prongs": 0, "truncated_prongs": 0}
        for P in range(1, Q):
            if math.gcd(P, Q) != 1:
                continue
            for eps in (0, 1):
                if eps != P % 2:                       # LONE only
                    continue
                r = scan_row(P, Q, eps, depth, maxstates)
                if r is None or "guard_bad" in r:
                    continue
                agg["rows"] += 1
                for k in ("rel", "triv", "form", "prongs"):
                    agg[k] += r[k]
                agg["closed_prongs"] += r["closed"]
                agg["truncated_prongs"] += r["truncated"]
                if r["rel"]:
                    out["violations"].append({"P": P, "Q": Q, "eps": eps,
                                              "hits": r["hits"][:2]})
        out["by_Q"][str(Q)] = agg
        out["rows"] += agg["rows"]
        for k in ("rel", "triv", "form"):
            out[k] += agg[k]
    out["prongs"] = sum(a["prongs"] for a in out["by_Q"].values())
    out["closed_prongs"] = sum(a["closed_prongs"] for a in out["by_Q"].values())
    out["truncated_prongs"] = sum(a["truncated_prongs"] for a in out["by_Q"].values())
    out["control_holds"] = (out["rel"] == 0)
    out["null_is_exhaustive"] = (out["closed_prongs"] == out["prongs"])
    out["defect_zero_bars_rel_only"] = {"rel": out["rel"], "triv": out["triv"],
                                        "form": out["form"]}
    return out


# ------------------------------------------------------------------ H5

def h5(qs=None, depth=DEPTH, maxstates=MAXSTATES):
    """THE INDEPENDENT VERIFICATION OF EVERY `rel` VECTOR, EXACTLY, OFF THE RING.

    ⚠⚠ THIS IS THE ARM THAT MAKES H2 QUOTABLE, and it is [WFLOOR-124]'s hazard applied to
    an ARRIVAL VECTOR rather than to a rank: `a` is produced by `s462_zz_relation`'s
    symbolic layer and reduced by `reduce_map`, both of which read `ch.ell[t].v` -- the
    SAME integer-vector representation the vanishing test uses.  A bug shared by the
    producer and the checker is invisible.

    So every distinct reduced vector is re-tested against `sin(c*pi/(2Q))` in an
    INDEPENDENT representation: with `zeta = zeta_{4Q}`, `2i*sin(c*pi/2Q) = zeta^c -
    zeta^{-c}`, so the claim is that `sum_i a_i (x^{c_i} - x^{4Q-c_i})` reduces to `0`
    modulo `Phi_{4Q}(x)`.  That is `s463_rank_law.reduction`, sympy polynomial remainder,
    sharing no code with `engine/ring_cyc.py`.

    ⚠⚠ AND IT IS DELIBERATELY *NOT* `sympy.simplify` ON A SUM OF SINES.  The first form of
    this arm was, and it reported `45/115` -- every `Q = 15` vector verified and every
    `Q = 21/33/39` vector "failed".  That was a FALSE NEGATIVE of the checker, not a defect
    in the data: `simplify` is a heuristic, not a decision procedure for a vanishing sum of
    roots of unity, and it gives up as the degree grows.  Reduction mod `Phi_{4Q}` IS a
    decision procedure.  ⇒ a checker that fails exactly where the instrument is new looks
    identical to an instrument that is wrong ([OPS-043]).
    """
    from s463_rank_law import reduction
    out = {"checked": 0, "vanishes": 0, "nonzero": 0, "bad": [], "by_Q": {}}
    for Q in (qs or TARGET_Q):
        R, f = reduction(4 * Q)
        agg = {"checked": 0, "vanishes": 0, "nonzero": 0}
        for P in range(1, Q):
            if math.gcd(P, Q) != 1:
                continue
            for eps in (0, 1):
                r = scan_row(P, Q, eps, depth, maxstates)
                if r is None or "guard_bad" in r or not r["hits"]:
                    continue
                ch = cm.Chain(P, Q, eps)
                red, order = reduce_map(ch)
                cs = [int(ch.c[t]) % (4 * Q) for t in order]
                seen = set()
                for h in r["hits"]:
                    key = tuple(h["a_red"])
                    if key in seen or not any(key):
                        continue
                    seen.add(key)
                    acc = [0] * f
                    for co, c in zip(key, cs):
                        if not co:
                            continue
                        hi, lo_ = R[c % (4 * Q)], R[(-c) % (4 * Q)]
                        for i in range(f):
                            acc[i] += int(co) * (hi[i] - lo_[i])
                    z = not any(acc)
                    agg["checked"] += 1
                    agg["vanishes"] += int(z)
                    agg["nonzero"] += 1
                    if not z:
                        out["bad"].append({"P": P, "Q": Q, "eps": eps,
                                           "a_red": list(key), "c": cs})
        out["by_Q"][str(Q)] = agg
        for k in ("checked", "vanishes", "nonzero"):
            out[k] += agg[k]
    out["all_verified"] = (out["checked"] > 0 and out["checked"] == out["vanishes"]
                           == out["nonzero"])
    return out


# ------------------------------------------------------------------ driver

def main():
    which = sys.argv[1] if len(sys.argv) > 1 and not sys.argv[1].startswith("-") else "all"
    depth, states, qcen, outf, qs = DEPTH, MAXSTATES, QCEN, OUT, None
    for a in sys.argv[1:]:
        if a.startswith("--out="):
            outf = a.split("=")[1]
        if a.startswith("--depth="):
            depth = int(a.split("=")[1])
        if a.startswith("--states="):
            states = int(a.split("=")[1])
        if a.startswith("--qcen="):
            qcen = int(a.split("=")[1])
        if a.startswith("--qs="):
            qs = [int(x) for x in a.split("=")[1].split(",")]
    out, t0 = {}, time.time()

    if which in ("all", "h0"):
        t = time.time()
        r = out["h0"] = h0()
        print(f"  H0  LONE  = 1 zero + (Q-1)/2 values of multiplicity 2 : "
              f"{r['lone_ok']}/{r['checked'] // 2}", flush=True)
        print(f"      PAIRED= 0 zeros + one singleton + the rest 2      : "
              f"{r['paired_ok']}/{r['checked'] // 2}", flush=True)
        print(f"      ⇒ RAW-basis relations exist at EVERY odd Q        : "
              f"{r['raw_relations_exist_everywhere']}", flush=True)
        print(f"        (so [NCYL-317] (4)'s `a != 0` in the RAW basis is not the "
              f"rank law's classifier; the three-way split replaces it)", flush=True)
        for s in r["sample"]:
            print(f"        Q={s['Q']:2d} {s['class']:6s} zeros={s['zeros']} "
                  f"ndist={s['ndist']} mult={s['mult']}", flush=True)
        for b in r["bad"][:6]:
            print("      bad:", b, flush=True)
        print(f"      [{time.time() - t:.1f}s]", flush=True)

    if which in ("all", "h1"):
        t = time.time()
        r = out["h1"] = h1(depth, states)
        for c in r["cases"]:
            print(f"  H1  {c['P']}/{c['Q']} eps{c['eps']}: relaxation admits  "
                  f"form={c['form']} triv={c['triv']} rel={c['rel']}"
                  f"  (min rel depth {c['min_rel_steps']})", flush=True)
        print(f"      ⇒ the real witness is admitted at depth <= 5      : {r['admits']}"
              f"   (a PASS is PREDICTED -- bug-check, not evidence)", flush=True)
        print(f"      ⇒ and it RE-CLASSIFIES as                          : "
              f"{r['witness_kind']}"
              f"   (this one CAN fail: `triv` would sink [NCYL-317] (2))", flush=True)
        print(f"        a     = {r['witness_a']}", flush=True)
        print(f"        a_red = {r['witness_a_red']}", flush=True)
        print(f"        matches [NCYL-317] (2)'s recorded vector          : "
              f"{r['matches_NCYL317_a']}", flush=True)
        s = r["soundness"]
        print(f"      G6 SOUNDNESS (untruncated prongs only): real-orbit states inside "
              f"the relaxation {s['states_covered']}/"
              f"{s['states_covered'] + s['states_missed']}"
              f"  -> sound={s['sound_on_untruncated']}", flush=True)
        print(f"      ⚠ {s['prongs_state_capped']} prongs hit the STATE cap and are "
              f"EXCLUDED -- a truncated BFS can drop the real walk, so it is not a "
              f"relaxation there and its null is worth less than a depth-capped one",
              flush=True)
        print(f"      [{time.time() - t:.1f}s]", flush=True)

    if which in ("all", "h2"):
        t = time.time()
        r = out["h2"] = h2(depth, states, qs)
        for Q in (qs or TARGET_Q):
            a = r["by_Q"][str(Q)]
            print(f"  H2  Q={Q:3d}  rows {a['rows']:3d}  prongs {a['prongs']:4d}  "
                  f"arrivals form={a['form']:5d} triv={a['triv']:5d} rel={a['rel']:4d}"
                  f"   rows_with_rel={a['rows_with_rel']:3d}"
                  f"   EXHAUSTED {a['closed_prongs']}/{a['prongs']}"
                  f"  trunc={a['truncated_prongs']}", flush=True)
        print(f"      ⇒ `rel` arrivals ONLY at Q = 15                    : "
              f"{r['only_15_has_rel']}", flush=True)
        print(f"      ⇒ the obstruction at Q = 9,21,27,33,39 is          : "
              f"{r['obstruction']}", flush=True)
        print(f"      ⚠ depth-capped at {r['depth']}: a NULL is not "
              f"'no word exists' ([OPS-041])", flush=True)
        for rr in r["rel_rows"][:8]:
            print(f"        rel row {rr['P']}/{rr['Q']} eps{rr['eps']}: "
                  f"{rr['rel']} hits, min depth {rr['min_steps']}", flush=True)
        print(f"      [{time.time() - t:.1f}s]", flush=True)

    if which in ("all", "h3"):
        t = time.time()
        r = out["h3"] = h3(qcen, depth, states)
        print(f"  H3  census odd Q <= {r['qcen']}: {r['rows']} rows, "
              f"{r['prongs']} corner prongs "
              f"(EXHAUSTED {r['closed_prongs']}, truncated {r['truncated_prongs']}), "
              f"{r['arrivals_total']} admitted arrivals", flush=True)
        print(f"      form={r['form']}  triv={r['triv']}  rel={r['rel']}"
              f"   rel share {r['rel_share']}", flush=True)
        for cl in ("lone", "paired"):
            c = r["by_class"][cl]
            print(f"      {cl:6s} rows {c['rows']:3d}: form={c['form']:5d} "
                  f"triv={c['triv']:5d} rel={c['rel']:4d}", flush=True)
        print(f"      ⇒ (start,target) PAIR admission "
              f"{r['targets']}/{r['pairs']} = {r['pair_admission_rate']}"
              f"    -- [NCYL-308] H3 plain lattice 0.8700, H6 augmented 0.7760",
              flush=True)
        print(f"      ⇒ prongs admitting ANY arrival "
              f"{r['prongs_admitting_any']}/{r['prongs']} "
              f"= {r['prong_admission_rate']}"
              f"   (of which `rel`: {r['prongs_admitting_rel']})", flush=True)
        print(f"      ⚠ the comparison is to those two, NOT to zero: this is a "
              f"RELAXATION, so the question is how much the +-1 structure buys",
              flush=True)
        print(f"      [{time.time() - t:.1f}s]", flush=True)

    if which in ("all", "h4"):
        t = time.time()
        r = out["h4"] = h4(31, depth, states)
        print(f"  H4  prime-Q LONE census (Q <= {r['qmax']}): {r['rows']} rows", flush=True)
        print(f"      rel={r['rel']}  triv={r['triv']}  form={r['form']}"
              f"   over {r['prongs']} corner prongs", flush=True)
        print(f"      ⇒ the null is EXHAUSTIVE, not depth-capped         : "
              f"{r['null_is_exhaustive']}"
              f"   ({r['closed_prongs']}/{r['prongs']} prongs closed, "
              f"{r['truncated_prongs']} truncated)", flush=True)
        print(f"      ⇒ [NCYL-318]'s zero relation space bars `rel` there: "
              f"{r['control_holds']}"
              f"   ({len(r['violations'])} violations)", flush=True)
        print(f"      ⚠ `triv`/`form` there are NOT failures -- they are the measurement "
              f"that qualifies [NCYL-317] (4)", flush=True)
        print(f"      [{time.time() - t:.1f}s]", flush=True)

    if which in ("all", "h5"):
        t = time.time()
        r = out["h5"] = h5(qs, depth, states)
        print(f"  H5  distinct `rel` reduced vectors re-tested in sympy, off the ring: "
              f"{r['checked']}", flush=True)
        print(f"      vanish exactly mod Phi_4Q                : "
              f"{r['vanishes']}/{r['checked']}", flush=True)
        print(f"      and are NONZERO (a genuine relation)     : "
              f"{r['nonzero']}/{r['checked']}", flush=True)
        print(f"      by Q: " + ", ".join(f"{q}:{a['vanishes']}/{a['checked']}"
                                          for q, a in r["by_Q"].items()), flush=True)
        print(f"      ⇒ every `rel` arrival is a real relation : {r['all_verified']}",
              flush=True)
        for b in r["bad"][:4]:
            print("      bad:", b, flush=True)
        print(f"      [{time.time() - t:.1f}s]", flush=True)

    out["secs"] = round(time.time() - t0, 1)
    out["params"] = {"depth": depth, "maxstates": states, "qcen": qcen}
    json.dump(out, open(outf, "w"), indent=1)
    print(f"  -> {outf}  [{out['secs']}s]", flush=True)


if __name__ == "__main__":
    main()
