#!/usr/bin/env python3.13
"""s467_form_lattice.py -- PRE-REGISTERED.  Queue item (s466-a): FEED THE TRANSVERSAL
INTO `(G0b-form)`.

PRIOR ART ([OPS-241], run before the probe was designed; the steer was gated as well as
the experiment).  `rulings.py --grep`: 'coefficient space' -> [NCYL-317] (which COINS
`(G0b-form)` as "honest `Lambda_S` membership in COEFFICIENT space" and does not test
it); 'fold-sum lattice' -> [NCYL-299], [NCYL-322]; 'lattice membership' -> none;
'codimension' -> none; 'formal identity' / 'integer identity' -> [NCYL-317]; 'value
class' / 'representative' / 'form arrival' -> [NCYL-319], [NCYL-322]; 'coefficient sum'
-> [NCYL-159], [NCYL-298]; 'parity invariant' / 'even sum' -> none.  ⊕ THE ID-GREPS
[OPS-246] PRESCRIBES -- `grep -rn '\\[NCYL-299\\]|\\[NCYL-308\\]|\\[NCYL-318\\]|
\\[NCYL-322\\]'` over `*.md probes/*.py engine/*.py` -- return `s455_cp_hypothesis.py`
h8/h9, `s453_pole_module.py` H3a, `s466_triv_gap.py` h3/h4, `s462/s464`.

⇒ WHAT CAME BACK, AND IT IS A GAP RATHER THAN AN ANSWER.  Every `Lambda`-membership
test this strand has run is in RING coordinates: [NCYL-308] H3 (plain, corner starts,
`87.0%` blind), H6 (the alternating-sign augmentation, `77.60%`), `s455` h8 (poles,
`0/267040`).  [NCYL-308]'s OWN s462 forward marker says why they die uniformly -- they
"ask MEMBERSHIP of a target in a lattice built from the `S_e` INSIDE the ring, and so
silently inherit whatever `Z`-relations the ring supplies".  **The COEFFICIENT-space
version -- the one that does NOT inherit them -- has never been run.**  `s466_triv_gap`
h3/h4 built the symbolic lattice but asked a different question of it (`Z*e_{t*}`
membership; the coset-meets-a-line relaxation).  ⇒ NOTHING BARS THIS AND NOTHING
ANSWERS IT.

⚠ [OPS-242]'s GATE, AND THE POINT IS THAT IT DOES **NOT** BIND HERE -- which is a design
fact worth stating, not a licence.  The gate kills a ROW-INDEPENDENT instrument when the
statement is FALSE without its hypothesis, and (G0b) is: the `7/15 eps1` connection is an
exact terminating counterexample.  But that connection is `rel` (`a_red != 0`), and
`(G0b-form)` has NO known counterexample anywhere -- all `4` observed `Z` are `rel`, `0`
are `form` ([NCYL-317] (4)).  So there is no row to run the candidate on, and a
row-independent obstruction to `(G0b-form)` is not a priori dead the way one to (G0b) is.
⚠ That is NOT evidence `(G0b-form)` is true: the base is ONE real object plus s464's
capped searches ([NCYL-319] scope (ii)-(iv)).

------------------------------------------------------------------------------------
THE OBJECT.  In [NCYL-294]'s model the fold at edge `e` is `y -> S_e - y`, so after `n`
folds `y_n = (-1)^n y_0 + sum_i (-1)^{n-i} S_{a_i}` ([NCYL-298]'s lemma).  Carried in
SYMBOLIC coordinates (`sym_build`/`sym_path`, integer vectors over `(ell_0..ell_{Q-1})`,
no ring) the arrival vector at a corner target is `a = y_n^sym - tgt`, and

    (G0b-form)  is  `a == 0` -- a FORMAL integer identity, no ring relation consumed.

`a == 0` with `n` EVEN reads `tgt - y_0 in Lambda_d`, with `n` ODD `tgt + y_0 in
sS[e] + Lambda_d`, where the ALTERNATING coefficients sum to `0` resp. `1` -- so the
honest necessary condition is the AUGMENTED one, `(tgt -/+ y_0, 0/1) in
span_Z{(sS[e], 1)}`, and the plain `span_Z{sS[e]}` test is its relaxation.  ⚠ NAME THE
LATTICE ([NCYL-299] amendment / [OPS-227]): `Lambda_S := span_Z{sS[e]}` is NOT
[NCYL-299]'s `Lambda := span_Z{sS[a] - sS[b]}`; H5 scores that they differ here too.

⇒⇒ WHAT THE TRANSVERSAL ADDS, WHICH IS THE WHOLE POINT OF THE ITEM.  [NCYL-322] (2):
every symbolic vector is supported inside `{0..t*}`, `t* = (Q-1)/2`, exactly one
representative per value class, because `sym_build`/`sym_path` live on the `iota`-
QUOTIENT path.  So the coefficient space is a CANONICAL `Z^d`, `d = (Q+1)/2`, and two
things become askable that are meaningless in the ring:
  (i)  a RANK/CODIMENSION count -- `rank Lambda_S` against `d` (H0);
  (ii) the COEFFICIENT-SUM functional `Sigma(v) = sum_t v_t`, which the ring destroys
       (there `ell`-relations move mass between coordinates) but which is well defined
       here.  DERIVED: `sym_build` sets `U[0] = 0` and steps by `e_far - e_near`
       (`Sigma = 0`), so `Sigma(LO[t]) = 0` and `Sigma(HI[t]) = 1` at every `t`; hence
       `Sigma(sS[e]) in {0, 2}` -- ALWAYS EVEN -- and `Lambda_S subset ker(Sigma mod 2)`.
       ⇒ **a `form` arrival would need `Sigma(tgt -/+ y_0)` EVEN, i.e. start and target
       of the SAME endpoint type.**  Cap-free, row-independent, every odd `Q`.
  ⚠ THIS IS A DIFFERENT SUM FROM [NCYL-298]/[NCYL-308] (3)'s.  Theirs sums the FOLD
  coefficients over edges (`sum_e c_e in {0,1}`); this sums the vector's own coordinates
  over interface indices.  Do not conflate them -- both appear in H2, on different axes.

⇒⇒ AND THE ARM THAT CAN KILL (ii), PRE-REGISTERED AS THE EXPECTED FAILURE ([OPS-243],
which s466 flagged one level up).  A corner is the SHARED endpoint of two intervals, so
if one geometric corner is BOTH some `Lo[i]` and some `Hi[i']`, its `Sigma` is `0` in one
labelling and `1` in the other -- the parity would be a property of the COORDINATES, not
of the objects, and the obstruction would be vacuous.  `sym_walk` itself resolves the
ambiguity by taking whichever matches first, `Lo` before `Hi`, which is exactly the tell.
H1's control decides it, and I expect it to fire.

------------------------------------------------------------------------------------
THE ARMS.  H0/H1/H2/H3/H4/H5 all scored.

  H0  THE RANK.  Per row: `m`, `d = (Q+1)/2`, `rank Lambda_S`, `codim = d - rank`, and
      the index `[Z^d : Lambda_S]` when full rank.  ⚠ COULD COME OUT EITHER WAY AND THAT
      IS THE POINT: `codim = 0` with index `1` means `Lambda_S` is ALL of `Z^d` and the
      route is blind before it starts; `codim > 0` is a genuine linear obstruction.

  H1  THE COEFFICIENT-SUM PREMISES, scored exactly rather than asserted:
      `Sigma(sLo[i]) == 0`, `Sigma(sHi[i]) == 1`, `Sigma(sS[e]) in {0,2}`, all rows.
      ⚠ THE CONTROL, WHICH IS THE ARM THAT DECIDES WHETHER (ii) MEANS ANYTHING: does any
      geometric corner carry BOTH labels (`Lo[i].v == Hi[i'].v` as exact ring vectors)?
      If yes, the LO/HI type is a coordinate artefact and the parity obstruction is void.

  H2  THE TEST.  Over all ordered corner pairs and both parities: admission of the
      AUGMENTED condition, and of the plain `Lambda_S` relaxation for comparability with
      [NCYL-308] H3's `87.0%`.  ⚠ DENOMINATOR STATED, and the degenerate `tgt == y_0` at
      even `n` (`0 in Lambda_S` always -- [OPS-230]'s identity-element tell, three
      instances already) is EXCLUDED and counted separately.  Split LONE/PAIRED x
      prime/composite `Q`, since the target region is LONE prime-`Q` ([NCYL-318]).

  H3  [OPS-230] DIRECTLY: are the symbolic endpoints themselves in `Lambda_S`?  In the
      RING they are, BY CONSTRUCTION at a corner, and that is precisely why [NCYL-308]
      H3 is `87%` blind.  If the same holds symbolically the route degenerates the same
      way and this probe says so.

  H4  POSITIVE CONTROL, MANDATORY ([OPS-041]).  A NECESSARY condition that REJECTS a real
      configuration is a bug, not a theorem.  Walk real corner prongs and check the
      condition at every visited state (it must hold by construction).  ⊕ re-score the
      stored `s464` arrivals: any recorded `a == 0` would be a `form` counterexample.
      ⚠ BOTH ARE PREDICTED -- bug-checks, not evidence.

  H5  THE NAMING AUDIT ([OPS-227]/[NCYL-299]'s amendment).  `rank span{sS[e]}` vs
      `rank span{sS[e] - sS[0]}` -- if they differ, the two objects are distinct here and
      `s466_triv_gap` h3's docstring, which calls `span_Z{sS[e]}` "[NCYL-299]'s object",
      is mislabelled.  ⚠ Scores the LABEL; [NCYL-322] H3's COMPUTATION is unaffected.

SCOPE.  Odd `Q` ([NCYL-291]).  Everything lives inside [NCYL-294]'s model, whose
global-coordinate closure is VERIFIED (`420/420`), not proved ([NCYL-322] scope (v)).
⚠ A necessary condition can only ever REFUTE; admission decides nothing.

USAGE:  PYTHONPATH=.:engine:archive/scripts_2026-07:archive/scripts_2026-08:probes \\
        .venv/bin/python3.13 probes/s467_form_lattice.py [h0|h1|h2|h3|h4|h5|all] \\
        [--qmax=N]
"""
from __future__ import annotations

import glob
import json
import math
import sys

import s451_quotient_path as qp
import s453_pole_module as pm
import s454_self_hit as sh
import s456_neckgb as nb
import zlattice as zl
from s462_zz_relation import sym_build, sym_path, sym_walk

OUT = "data/s467_form_lattice.json"
QMAX = 33


def rows(qmax):
    for Q in range(5, qmax + 1, 2):
        for P in range(1, Q):
            if math.gcd(P, Q) == 1:
                for eps in (0, 1):
                    yield P, Q, eps


def build(P, Q, eps):
    """The full pack: chain, ring path view, symbolic path view."""
    ch, lo, hi, closes = qp.build(P, Q, eps)
    if not closes:
        return None
    sf = pm.fold_sums(ch, lo, hi)
    Lo, Hi, S = pm.path_view(ch, lo, hi, sf)
    ops = pm.CoordOps(ch)
    LO, HI = sym_build(ch)
    sLo, sHi, sS = sym_path(ch, LO, HI, lo, hi)
    iv = lambda A: [[int(x) for x in v] for v in A]          # noqa: E731
    return ch, lo, hi, Lo, Hi, S, ops, iv(sLo), iv(sHi), iv(sS)


def trunc(v, d):
    """Restrict to the transversal `{0..t*}` (asserted empty outside by [NCYL-322])."""
    return list(v[:d])


def is_prime(n):
    return n > 1 and all(n % k for k in range(2, int(n ** 0.5) + 1))


def sigma(v):
    return sum(v)


# ------------------------------------------------------------------ H0
def h0(qmax, verbose=True):
    n, bad_supp, law, tab = 0, [], {}, []
    for P, Q, eps in rows(qmax):
        b = build(P, Q, eps)
        if b is None:
            continue
        ch, _lo, _hi, _Lo, _Hi, _S, _ops, sLo, sHi, sS = b
        d, m = (Q + 1) // 2, ch.m_idx
        n += 1
        for v in sLo + sHi + sS:
            if any(v[d:]):
                bad_supp.append((P, Q, eps))
                break
        L = zl.Lattice([trunc(v, d) for v in sS])
        r = L.rank
        # index [Z^d : Lambda_S] = product of echelon pivot entries, when full rank
        idx = None
        if r == d:
            idx = 1
            for pv, row in L.basis:
                idx *= abs(row[pv])
        key = (d - r, "full" if r == d else "deficient")
        law[key] = law.get(key, 0) + 1
        tab.append([P, Q, eps, int(ch.lone), m, d, r, d - r, idx])
    codims = sorted({t[7] for t in tab})
    if verbose:
        print(f"H0  the transversal is the ambient: every symbolic vector supported "
              f"inside {{0..t*}} on {n - len(bad_supp)}/{n} rows "
              f"{'PASS' if not bad_supp else 'FAIL'}")
        print(f"    codim(Lambda_S in Z^d), d = (Q+1)/2:  observed {codims}")
        for c in codims:
            sub = [t for t in tab if t[7] == c]
            idxs = sorted({t[8] for t in sub if t[8] is not None})
            print(f"      codim {c}: {len(sub)}/{n} rows"
                  + (f"   index [Z^d : Lambda_S] in {idxs[:6]}"
                     f"{' ...' if len(idxs) > 6 else ''}" if idxs else ""))
        ms = sorted({(t[4], t[5]) for t in tab})
        print(f"    (m, d) pairs: {ms[:8]}{' ...' if len(ms) > 8 else ''}")
        if codims == [0] and all(t[8] == 1 for t in tab):
            print("    ⇒ ⚠ Lambda_S IS ALL OF Z^d -- the membership route is BLIND "
                  "before it starts.")
    return {"rows": n, "bad_support": bad_supp[:10], "codims": codims,
            "table": tab, "pass": not bad_supp}


# ------------------------------------------------------------------ H1
def h1(qmax, verbose=True):
    n, bad, both_lab, both_rows = 0, [], 0, []
    lone_clash, zero_iface = 0, 0
    for P, Q, eps in rows(qmax):
        b = build(P, Q, eps)
        if b is None:
            continue
        _ch, _lo, _hi, Lo, Hi, _S, _ops, sLo, sHi, sS = b
        n += 1
        if any(sigma(v) != 0 for v in sLo):
            bad.append((P, Q, eps, "sigma(LO)!=0"))
        if any(sigma(v) != 1 for v in sHi):
            bad.append((P, Q, eps, "sigma(HI)!=1"))
        if any(sigma(v) not in (0, 2) for v in sS):
            bad.append((P, Q, eps, "sigma(sS) odd", [sigma(v) for v in sS]))
        # ⚠ THE CONTROL: is the LO/HI label of a geometric corner canonical?
        clash = [(i, j) for i in range(len(Lo)) for j in range(len(Hi))
                 if (Lo[i].v == Hi[j].v).all()]
        if clash:
            both_lab += 1
            lone_clash += _ch.lone
            # is EVERY clash the degenerate (zero-length) interval, i.e. the zero
            # interface at t* that only the LONE class carries ([NCYL-292] (vi))?
            if all((Lo[j].v == Hi[j].v).all() for _i, j in clash):
                zero_iface += 1
            if len(both_rows) < 6:
                both_rows.append((P, Q, eps, int(_ch.lone), clash[:3]))
    if verbose:
        print(f"H1  Sigma(sLo) == 0, Sigma(sHi) == 1, Sigma(sS) in {{0,2}} (always EVEN) "
              f"-- {n - len(bad)}/{n} rows   "
              f"{'PASS' if not bad else 'FAIL ' + str(bad[:3])}")
        print(f"    ⇒ Lambda_S subset ker(Sigma mod 2), so a `form` arrival needs "
              f"Sigma(tgt -/+ y_0) EVEN.")
        print(f"    ⚠ CONTROL -- is the LO/HI label of a corner canonical?  a corner "
              f"value carries BOTH labels on {both_lab}/{n} rows, of which "
              f"{lone_clash} are LONE and {zero_iface} are the DEGENERATE interval "
              f"(Lo[j] == Hi[j]) alone")
        if both_lab and lone_clash == both_lab and zero_iface == both_lab:
            print(f"      ⇒⇒ AND THAT IS THE WHOLE STORY: the ambiguity is EXACTLY the "
                  f"zero interface at t*, which only LONE carries ([NCYL-292] (vi)), "
                  f"and the two labels differ by e_t* -- [NCYL-322]'s redundant "
                  f"coordinate, not a new one.")
            print(f"      ⇒ SO Sigma mod 2 IS WELL DEFINED ON PAIRED and the parity "
                  f"obstruction BITES there; on LONE it is broken by e_t* "
                  f"(Sigma(e_t*) = 1).")
        elif both_lab:
            print(f"      ⇒⇒ THE PARITY OBSTRUCTION IS VOID: Sigma is a property of the "
                  f"LABELLING, not of the corner ([OPS-243]).  e.g. {both_rows[:2]}")
        else:
            print(f"      ⇒⇒ THE LABEL IS CANONICAL, so the parity obstruction is about "
                  f"the OBJECTS and it BITES.")
    return {"rows": n, "n_bad": len(bad), "bad": bad[:10],
            "rows_with_both_labels": both_lab, "lone_clash": lone_clash,
            "zero_interface_only": zero_iface, "examples": both_rows,
            "parity_bites_on_paired": both_lab == lone_clash == zero_iface,
            "pass": not bad}


# ------------------------------------------------------------------ H2
def h2(qmax, verbose=True):
    agg = {}
    tot = {"tests": 0, "aug_hits": 0, "plain_hits": 0, "degenerate": 0}
    for P, Q, eps in rows(qmax):
        b = build(P, Q, eps)
        if b is None:
            continue
        ch, _lo, _hi, _Lo, _Hi, _S, _ops, sLo, sHi, sS = b
        d = (Q + 1) // 2
        gen = [trunc(v, d) for v in sS]
        L = zl.Lattice(gen)                       # Lambda_S
        La = zl.Lattice([g + [1] for g in gen])   # the augmented lattice
        ends = [trunc(v, d) for v in sLo] + [trunc(v, d) for v in sHi]
        cell = agg.setdefault((bool(ch.lone), is_prime(Q)),
                              {"rows": 0, "tests": 0, "aug": 0, "plain": 0, "deg": 0})
        cell["rows"] += 1
        for y0 in ends:
            for tgt in ends:
                for par in (0, 1):                # n even / n odd
                    v = ([t - x for t, x in zip(tgt, y0)] if par == 0
                         else [t + x for t, x in zip(tgt, y0)])
                    if par == 0 and not any(v):
                        cell["deg"] += 1
                        tot["degenerate"] += 1
                        continue
                    cell["tests"] += 1
                    tot["tests"] += 1
                    if La.contains(v + [par]):
                        cell["aug"] += 1
                        tot["aug_hits"] += 1
                    if L.contains(v):
                        cell["plain"] += 1
                        tot["plain_hits"] += 1
    if verbose:
        print(f"H2  necessary condition for `form`, over ORDERED corner pairs x both "
              f"parities.  denominator = {tot['tests']} (start,target,parity) triples; "
              f"{tot['degenerate']} degenerate (tgt == y_0, n even) EXCLUDED "
              f"([OPS-230]).")
        pa = 100.0 * tot["aug_hits"] / max(1, tot["tests"])
        pp = 100.0 * tot["plain_hits"] / max(1, tot["tests"])
        print(f"    AUGMENTED (the honest condition): {tot['aug_hits']}/{tot['tests']} "
              f"= {pa:.2f}% admitted")
        print(f"    plain Lambda_S (its relaxation):  {tot['plain_hits']}/"
              f"{tot['tests']} = {pp:.2f}%   "
              f"⚠ compare [NCYL-308] H3's 87.0% in RING coordinates -- DIFFERENT "
              f"denominator, same shape of question")
        print(f"    {'class':7s} {'Q':9s} {'rows':>5s} {'tests':>8s} {'aug':>8s} "
              f"{'aug%':>7s} {'plain%':>7s}")
        for (lone, pr) in sorted(agg, key=lambda k: (not k[0], not k[1])):
            c = agg[(lone, pr)]
            print(f"    {'LONE' if lone else 'PAIRED':7s} "
                  f"{'prime' if pr else 'composite':9s} {c['rows']:5d} "
                  f"{c['tests']:8d} {c['aug']:8d} "
                  f"{100.0 * c['aug'] / max(1, c['tests']):6.2f}% "
                  f"{100.0 * c['plain'] / max(1, c['tests']):6.2f}%")
    return {"total": tot,
            "cells": {f"{'LONE' if k[0] else 'PAIRED'}/"
                      f"{'prime' if k[1] else 'composite'}": v
                      for k, v in agg.items()},
            "pass": True}


# ------------------------------------------------------------------ H3
def h3(qmax, verbose=True):
    n, ends_in, ends_tot = 0, 0, 0
    allin = 0
    for P, Q, eps in rows(qmax):
        b = build(P, Q, eps)
        if b is None:
            continue
        _ch, _lo, _hi, _Lo, _Hi, _S, _ops, sLo, sHi, sS = b
        d = (Q + 1) // 2
        L = zl.Lattice([trunc(v, d) for v in sS])
        n += 1
        k = 0
        ends = [trunc(v, d) for v in sLo] + [trunc(v, d) for v in sHi]
        for e in ends:
            ends_tot += 1
            if L.contains(e):
                ends_in += 1
                k += 1
        allin += (k == len(ends))
    verdict = ("SAME DEGENERACY HERE." if allin == n
               else "NOT the same degeneracy here.")
    if verbose:
        print(f"H3  [OPS-230]: are the symbolic corner endpoints themselves in "
              f"Lambda_S?  {ends_in}/{ends_tot} endpoints "
              f"({100.0 * ends_in / max(1, ends_tot):.2f}%), and ALL of them on "
              f"{allin}/{n} rows")
        print(f"    ⚠ in the RING they are, by construction, which is why [NCYL-308] H3 "
              f"is 87% blind.  {verdict}")
    return {"rows": n, "endpoints_in": ends_in, "endpoints": ends_tot,
            "rows_all_in": allin, "pass": True}


# ------------------------------------------------------------------ H4
def h4(qmax, verbose=True):
    """POSITIVE CONTROL: the condition must ADMIT every REAL arrival.

    At a real `Z` the walk gives `a = y_n - tgt` with `y_n = (-1)^n y0 + sum(+-1) sS`,
    so `(a + tgt) - (-1)^n y0` lies in `Lambda_S` at augmented coordinate `n mod 2`,
    BY CONSTRUCTION.  `sym_walk` does not return `tgt`, so the checkable form is the
    EXISTENTIAL one over the endpoint set (the real target is in it): if NO endpoint
    works, the walk, the symbolic layer or the lattice is broken.
    """
    ok, dead, rowsn, zs = 0, 0, 0, 0
    for P, Q, eps in rows(qmax):
        if Q > 21:            # walking is the cost; this is a control, not a census
            continue
        b = build(P, Q, eps)
        if b is None:
            continue
        ch, _lo, _hi, Lo, Hi, S, ops, sLo, sHi, sS = b
        d = (Q + 1) // 2
        La = zl.Lattice([trunc(v, d) + [1] for v in sS])
        rowsn += 1
        for j in range(ch.m_idx):
            cs = nb.corner_start(Lo, Hi, S, j, ops)
            if not cs:
                continue
            end, a, folds, _steps = sym_walk(ch, Lo, Hi, S, sLo, sHi, sS, cs[0], ops,
                                             cap=20_000)
            if end != "Z" or a is None:
                continue
            zs += 1
            par, sgn = folds % 2, (-1) ** (folds % 2)
            _i0, _d0, y0 = cs[0]
            y0s = None
            for idx in range(len(Lo)):
                if ops.eq(y0, Lo[idx]):
                    y0s = trunc(sLo[idx], d)
                    break
                if ops.eq(y0, Hi[idx]):
                    y0s = trunc(sHi[idx], d)
                    break
            if y0s is None:
                continue
            av = trunc(a, d)
            hit = False
            for t in [trunc(v, d) for v in sLo] + [trunc(v, d) for v in sHi]:
                w = [ai + ti - sgn * yi for ai, ti, yi in zip(av, t, y0s)]
                if La.contains(w + [par]):
                    hit = True
                    break
            ok += hit
            dead += not hit
    if dead:
        verdict = "FAIL -- rejects a real arrival; the instrument is UNSOUND"
    elif ok:
        verdict = "PASS"
    else:
        verdict = "DEAD -- no real Z arrival in range; the control is VACUOUS"
    if verbose:
        print(f"H4  POSITIVE CONTROL on {rowsn} rows (Q <= 21), {zs} real `Z` arrivals: "
              f"admitted {ok}, rejected {dead}   {verdict}")
        print("    ⚠ PREDICTED if it passes -- a bug-check, not evidence ([OPS-041]).")
    # any stored arrival with a == 0 would BE a `form` counterexample
    seen, forms = 0, 0
    for fn in sorted(glob.glob("data/s46*.json")):
        try:
            blob = json.load(open(fn))
        except Exception:
            continue
        for a in _walk_a(blob):
            seen += 1
            forms += not any(a)
    note = ("none -- consistent with [NCYL-317] (4)" if not forms
            else "⚠⚠ FORM COUNTEREXAMPLE PRESENT")
    if verbose:
        print(f"    ⊕ re-scored {seen} stored arrival vectors: {forms} have a == 0   "
              f"{note}")
    return {"rows": rowsn, "z_arrivals": zs, "admitted": ok, "rejected": dead,
            "stored_scored": seen, "stored_form": forms,
            "pass": dead == 0 and ok > 0}


def _walk_a(o):
    if isinstance(o, dict):
        if isinstance(o.get("a"), list):
            yield o["a"]
        for v in o.values():
            yield from _walk_a(v)
    elif isinstance(o, list):
        for v in o:
            yield from _walk_a(v)


# ------------------------------------------------------------------ H5
def h5(qmax, verbose=True):
    n, differ = 0, 0
    ex = []
    for P, Q, eps in rows(qmax):
        b = build(P, Q, eps)
        if b is None:
            continue
        _ch, _lo, _hi, _Lo, _Hi, _S, _ops, _sLo, _sHi, sS = b
        d = (Q + 1) // 2
        g = [trunc(v, d) for v in sS]
        rs = zl.Lattice(g).rank
        rd = (zl.Lattice([[a - b0 for a, b0 in zip(v, g[0])] for v in g[1:]]).rank
              if len(g) > 1 else 0)
        n += 1
        if rs != rd:
            differ += 1
            if len(ex) < 4:
                ex.append((P, Q, eps, rs, rd))
    if verbose:
        print(f"H5  NAMING AUDIT ([OPS-227]): rank span{{sS[e]}} vs "
              f"rank span{{sS[e] - sS[0]}} differ on {differ}/{n} rows   e.g. "
              f"{ex[:3]}")
        print(f"    ⇒ they are DIFFERENT objects here too, so `s466_triv_gap` h3's "
              f"\"Lambda = span_Z{{sS[e]}} ([NCYL-299]'s object)\" is a MISLABEL; "
              f"[NCYL-299]'s object is the DIFFERENCE lattice.  Its COMPUTATION "
              f"(zl.Lattice(sS)) is unaffected.")
    return {"rows": n, "differ": differ, "examples": ex, "pass": True}


# ------------------------------------------------------------------ H6
def _kernel_functional(gen, d):
    """The primitive integer `lam` in Z^d with `lam . g == 0` for every generator.

    `codim == 1` (H0) makes it unique up to sign.  Integer nullspace by fraction-free
    elimination -- no sympy, so it shares no code with the rank path above.
    """
    from fractions import Fraction
    M = [[Fraction(x) for x in g] for g in gen]
    piv, r = [], 0
    for c in range(d):
        pr = next((i for i in range(r, len(M)) if M[i][c]), None)
        if pr is None:
            continue
        M[r], M[pr] = M[pr], M[r]
        pv = M[r][c]
        M[r] = [x / pv for x in M[r]]
        for i in range(len(M)):
            if i != r and M[i][c]:
                f = M[i][c]
                M[i] = [x - f * y for x, y in zip(M[i], M[r])]
        piv.append(c)
        r += 1
        if r == len(M):
            break
    free = [c for c in range(d) if c not in piv]
    if len(free) != 1:
        return None
    f0 = free[0]
    lam = [Fraction(0)] * d
    lam[f0] = Fraction(1)
    for i, c in enumerate(piv):
        lam[c] = -M[i][f0]
    den = 1
    for x in lam:
        den = den * x.denominator // math.gcd(den, x.denominator)
    lam = [int(x * den) for x in lam]
    g = 0
    for x in lam:
        g = math.gcd(g, abs(x))
    lam = [x // g for x in lam]
    return lam if next((x for x in lam if x), 0) > 0 else [-x for x in lam]


def h6(qmax, verbose=True):
    """WHAT IS THE codim-1 FUNCTIONAL?  If it has a closed form, the obstruction is
    explicit and cap-free."""
    n, shapes, bad, ex = 0, {}, 0, []
    pm1 = {'n': 0, 'all_pm1': 0, 'alternating': 0, 'constant': 0}
    for P, Q, eps in rows(qmax):
        b = build(P, Q, eps)
        if b is None:
            continue
        ch, _lo, _hi, _Lo, _Hi, _S, _ops, sLo, sHi, sS = b
        d, ts = (Q + 1) // 2, (Q - 1) // 2
        gen = [trunc(v, d) for v in sS]
        lam = _kernel_functional(gen, d)
        n += 1
        if lam is None or any(sum(l * g for l, g in zip(lam, gg)) for gg in gen):
            bad += 1
            continue
        # DESCRIBE it: is it +- e_{t*}?  is it constant?  what are its values?
        vals = sorted(set(lam))
        nz = [t for t, x in enumerate(lam) if x]
        tag = ("e_t*" if nz == [ts] and abs(lam[ts]) == 1
               else "single" if len(nz) == 1
               else "constant" if len(vals) == 1
               else f"vals{vals}")
        key = (tag, bool(ch.lone))
        shapes[key] = shapes.get(key, 0) + 1
        pm1["all_pm1"] += all(abs(x) == 1 for x in lam)
        pm1["alternating"] += (lam == [(-1) ** t for t in range(d)]
                               or lam == [-(-1) ** t for t in range(d)])
        pm1["constant"] += len(set(lam)) == 1
        pm1["n"] += 1
        if len(ex) < 6 and Q <= 15:
            ex.append((P, Q, eps, int(ch.lone), lam))
    if verbose:
        print(f"H6  the codim-1 kernel functional lam (lam . sS[e] == 0 for all e), "
              f"recomputed fraction-free: found on {n - bad}/{n} rows")
        for (tag, lone), c in sorted(shapes.items(), key=lambda kv: -kv[1]):
            print(f"      {'LONE' if lone else 'PAIRED':7s} lam = {tag:12s} {c} rows")
        print(f"      ⇒ lam is a SIGN VECTOR (|lam_t| == 1 for every t) on "
              f"{pm1['all_pm1']}/{pm1['n']} rows; alternating (-1)^t on "
              f"{pm1['alternating']}; constant (i.e. lam == Sigma) on "
              f"{pm1['constant']}")
        for e in ex[:4]:
            print(f"      P={e[0]} Q={e[1]} eps={e[2]} "
                  f"{'LONE' if e[3] else 'PAIRED'}: lam = {e[4]}")
    return {"rows": n, "no_functional": bad,
            "shapes": {f"{k[0]}/{'LONE' if k[1] else 'PAIRED'}": v
                       for k, v in shapes.items()},
            "examples": ex, "sign_vector": pm1, "pass": bad == 0}


# ------------------------------------------------------------------ H7
def h7(qmax, verbose=True):
    """⇒⇒ THE DECIDING ARM.  On LONE, `form` is `a in Z*e_{t*}`, NOT `a == 0`
    ([NCYL-322] (4)) -- so the honest test there is membership in `Lambda_S + Z e_{t*}`.
    [NCYL-322] H3 says `e_{t*}` is NOT in `Lambda_S`, so that sum has rank `d`: the
    codim-1 obstruction is CONSUMED exactly on the class where `(G0b-form)` is the whole
    residual.  Scored, not asserted."""
    n, full, tot, hits = 0, 0, 0, 0
    paired_tot, paired_hits = 0, 0
    for P, Q, eps in rows(qmax):
        b = build(P, Q, eps)
        if b is None:
            continue
        ch, _lo, _hi, _Lo, _Hi, _S, _ops, sLo, sHi, sS = b
        d, ts = (Q + 1) // 2, (Q - 1) // 2
        gen = [trunc(v, d) for v in sS]
        em = [0] * d
        em[ts] = 1
        ends = [trunc(v, d) for v in sLo] + [trunc(v, d) for v in sHi]
        if ch.lone:
            n += 1
            Lp = zl.Lattice([g + [1] for g in gen] + [em + [0]])
            full += (zl.Lattice(gen + [em]).rank == d)
            for y0 in ends:
                for tgt in ends:
                    for par in (0, 1):
                        v = ([t - x for t, x in zip(tgt, y0)] if par == 0
                             else [t + x for t, x in zip(tgt, y0)])
                        if par == 0 and not any(v):
                            continue
                        tot += 1
                        hits += Lp.contains(v + [par])
        else:
            La = zl.Lattice([g + [1] for g in gen])
            for y0 in ends:
                for tgt in ends:
                    for par in (0, 1):
                        v = ([t - x for t, x in zip(tgt, y0)] if par == 0
                             else [t + x for t, x in zip(tgt, y0)])
                        if par == 0 and not any(v):
                            continue
                        paired_tot += 1
                        paired_hits += La.contains(v + [par])
    if verbose:
        print(f"H7  ⇒⇒ [NCYL-322]'s LONE CORRECTION: on LONE, `form` means "
              f"`a in Z*e_t*`, so the test is `Lambda_S + Z e_t*`.")
        print(f"    rank(Lambda_S + Z e_t*) == d on {full}/{n} LONE rows   "
              f"{'⇒ THE codim-1 OBSTRUCTION IS CONSUMED ON LONE.' if full == n else ''}")
        print(f"    corrected LONE admission: {hits}/{tot} = "
              f"{100.0 * hits / max(1, tot):.2f}%   "
              f"(uncorrected LONE reads ~16.3% -- H2)")
        print(f"    PAIRED, unchanged (no zero interface): {paired_hits}/{paired_tot} = "
              f"{100.0 * paired_hits / max(1, paired_tot):.2f}%")
    return {"lone_rows": n, "lone_full_rank": full,
            "lone_tests": tot, "lone_hits": hits,
            "paired_tests": paired_tot, "paired_hits": paired_hits,
            "pass": True}


# ------------------------------------------------------------------
def main():
    arms, qmax = "all", QMAX
    for arg in sys.argv[1:]:
        if arg.startswith("--qmax="):
            qmax = int(arg.split("=")[1])
        else:
            arms = arg
    res = {"qmax": qmax}
    print("s467_form_lattice.py -- (s466-a): feed the transversal into `(G0b-form)`\n")
    for k, fn in (("h0", h0), ("h1", h1), ("h2", h2), ("h3", h3), ("h4", h4),
                  ("h5", h5), ("h6", h6), ("h7", h7)):
        if arms in ("all", k):
            res[k] = fn(qmax)
            print()
    if arms == "all":
        json.dump(res, open(OUT, "w"), indent=1)
        print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
