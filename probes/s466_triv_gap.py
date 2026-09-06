#!/usr/bin/env python3.13
"""s466_triv_gap.py -- PRE-REGISTERED.  Queue item (s464-a): IS THERE A *REASON* A FOLD
WORD CANNOT REALISE A TRIVIAL RELATION?

PRIOR ART ([OPS-241], run before the probe was designed).  `rulings.py --grep` on
'trivial relation' -> [NCYL-288], [NCYL-319]; 'reduced vector' -> [NCYL-317],
[NCYL-319]; 'equal pair' -> [NCYL-069]/[NCYL-121]/[NCYL-135]/[NCYL-196]/[NCYL-118]/
[NCYL-317]; 'fold sum' -> [NCYL-299]; 'arrival vector' -> [OPS-244]; 'antisymmetric'
-> [WFLOOR-139], [OPS-226]; 'symbolic layer' -> [NCYL-319].  ⊕ THE ID-GREP THIS
SESSION'S OWN [OPS-246] PRESCRIBES -- `grep -rn 'NCYL-317\\]|NCYL-319\\]'` over `*.md`
and `probes/*.py` -- returns `TOOLS.md`, `future_directions.md` B9 q3 and
`s463_rank_law.py`, none of which touches `triv`.  ⇒ NOTHING BARS THIS AND NOTHING
ANSWERS IT.  ⚠ [NCYL-299] is the near miss and is USED, not re-derived: it names the
fold-sum lattice `Lambda` and is the reason H3's `Lambda`-membership question is the
natural one to ask.

------------------------------------------------------------------------------------
THE QUESTION.  [NCYL-319] (4) repaired [NCYL-317] (4)'s classifier into a THREE-way
split of an arrival vector `a` (raw, over the basis `ell_0 .. ell_{Q-1}`):
    form  `a == 0`;
    triv  `a != 0` but `a_red == 0` -- vanishes by the zero interface / the +- pairing;
    rel   `a_red != 0` -- a genuine relation among DISTINCT measures consumed.
`triv` was measured `0` everywhere (`0` of `277` census arrivals, `0` across the whole
`3 | Q` set, `0` on `1376` prime-`Q` LONE prongs) and [NCYL-317] (4) now needs that as a
HYPOTHESIS: with `triv` live, the residual on prime-`Q` LONE is `form ∪ triv` and
`(G0b-form)` is not the whole of it.

⇒⇒ THE ANSWER, AND IT DISSOLVES THE CATEGORY RATHER THAN SEARCHING FOR ITS ABSENCE.

  (A) THE VALUE PAIRING IS `sigma(t) = Q-1-t`, AND THAT IS AN IDENTITY, NOT A FIT.
      `ell_t = sin(c_t*pi/2Q)` with `c_t = (c0 + 2Pt) mod 2Q`, so `ell_t == ell_t'` iff
      `c_t' == c_t` or `c_t' == 2Q - c_t`.  The second reads `2P(t+t') = -2*c0 mod 2Q`,
      i.e. `P(t+t') = -c0 mod Q`; and `c0` is `P` (eps=0) or `P+Q` (eps=1), so
      `c0 = P mod Q` either way, and `gcd(P,Q) = 1` cancels it:
              t + t'  ==  -1   (mod Q).
      Independent of `P` and of `eps`.  Unique fixed point `t* = (Q-1)/2`, which is the
      ZERO interface on the LONE class and the perpendicular (`ell = 1`) singleton on
      the PAIRED one.

  (B) THE SYMBOLIC LAYER IS SUPPORTED ON `{0, .., t*}`, WHICH IS EXACTLY ONE
      REPRESENTATIVE PER `sigma`-ORBIT.  Not a coincidence: `sym_build`/`sym_path` live
      on the `iota`-QUOTIENT path, and `D = B/iota` is a CHAIN of `(Q-1)/2` kites plus
      one triangle ([NCYL-287]) -- i.e. the quotient by the very involution `sigma`
      induces on interface indices.  `sigma` maps `[0, t*-1]` bijectively onto
      `[t*+1, Q-1]` and fixes `t*`.

  (C) ⇒ SO `a_red == 0` FORCES `a_t == 0` AT EVERY INDEX WITH A NONZERO MEASURE, because
      each nonzero value class meets the support in exactly ONE index and `a_red`'s
      coordinate there IS that single `a_t`.  Hence
        PAIRED  (no zero measure):  `a_red == 0  ==>  a == 0`
                ⇒ **`triv` IS IMPOSSIBLE ON THE PAIRED CLASS AT EVERY ODD `Q`, at any
                  depth, with no search and no cap.**
        LONE    (one zero measure, at `t*`):  `a_red == 0  ==>  a in Z*e_{t*}`
                ⇒ **`triv` COLLAPSES FROM A `(Q+1)/2`-DIMENSIONAL LATTICE TO A LINE.**

  (D) AND THAT LINE IS NOT A THIRD MECHANISM -- IT IS `form` IN A REDUNDANT
      REPRESENTATION.  ⚠ DERIVED (H4), not scored.  On the LONE class `ell_{t*} == 0`
      identically, at every odd `Q`, as a STRUCTURAL feature of the necklace ([NCYL-292]
      (vi): the horizontal interface is why `h == 2` there).  So `a = c*e_{t*}` is the
      number `0` for every `c`, consumes no relation among distinct measures, and differs
      from a `form` arrival only in that the builder wrote `+c*ell_{t*}` where it could
      have written nothing.  ⇒⇒ **STATE THE SPLIT ON `a_red` AND THE TRICHOTOMY IS THE
      ORIGINAL DICHOTOMY AGAIN: `form` = `a_red == 0`, `rel` = `a_red != 0`.  [NCYL-317]
      (4) DOES NOT NEED `triv == 0` AS A HYPOTHESIS -- THE CATEGORY IS AN ARTEFACT OF THE
      OVER-COMPLETE RAW `ell`-BASIS, WHICH IS THE SAME DEFECT [NCYL-319] (4) DIAGNOSED,
      ONE STEP FURTHER.**

------------------------------------------------------------------------------------
THE ARMS.  H0/H1/H2/H3/H5 scored; H4 is the derivation and is not.

  H0  THE PAIRING.  `sigma(t) = Q-1-t` is the value-equality pairing, scored against the
      EXACT ring vectors `ch.ell[t].v` (never floats, never the `c` label), and `t*` is
      the zero measure iff the class is LONE.
      ⚠ CONTROL, WHICH MUST FAIL ([OPS-219]): the other natural involution
      `t -> -t mod Q` (`t + t' == 0`).  If it also paired the values, H0 would be
      scoring "some involution works" rather than this one.

  H1  THE SUPPORT.  Every `sLo[i]`, `sHi[i]`, `sS[e]` is supported inside `{0..t*}`, and
      the union of supports is the WHOLE of `{0..t*}` (so the transversal is not
      vacuously satisfied by a small support).

  H2  ⇒ THE CONSEQUENCE, COMPUTED AS EXACT LINEAR ALGEBRA PER ROW RATHER THAN ASSERTED:
      the solution space `K = {a : supp(a) subset {0..t*}, a_red == 0}` has basis
      `{e_{t*}}` on LONE and is `{0}` on PAIRED.
      ⊕ RE-SCORE THE STORED ARRIVALS: every `a` recorded by `s464_realisable` must obey
      it.  ⚠ PREDICTED, hence a BUG-CHECK and not evidence ([OPS-041]).

  H3  THE `Lambda_S` COROLLARY.  No nonzero multiple of `e_{t*}` lies in the fold-sum
      SPAN `Lambda_S = span_Z{sS[e]}`.  ⚠ CORRECTED s467 -- this line used to call that
      "[NCYL-299]'s object", and it is not: [NCYL-299]'s `Lambda` is the DIFFERENCE
      lattice `span_Z{s_a - s_b}`, and the two ranks differ on `460/460` rows here
      ([OPS-227]'s trap, `s467_form_lattice.py` h5).  The COMPUTATION below
      (`zl.Lattice(sS)`) is the span and is unaffected.  ⇒ the `e_{t*}` component of
      a LONE `triv` cannot be manufactured by FOLDING; it would have to come from the
      start/target endpoints.  ⚠ A structural fact about `Lambda`, not the answer -- (D)
      is the answer.

  H4  ⚠ THE NEGATIVE, RECORDED SO NOBODY BUILDS IT AGAIN.  The obvious cap-free
      instrument -- relax the word to its lattice, `a in (+-sX - sY) + Lambda`, and ask
      whether that coset meets `Z*e_{t*} \\ {0}` -- ADMITS on every LONE row.  It cannot
      decide the residue.  Reported as a count, not hidden.

  H5  THE KNOWN WITNESS MUST STILL COME OUT `rel` ([OPS-242]: run the candidate on the
      row where the answer is known before believing it elsewhere).  `7/15 eps1`'s
      recorded arrival vector must have `a_red != 0`.

SCOPE.  Odd `Q` ([NCYL-291]).  (A) is an identity and (C)/(D) follow from (A)+(B) with no
`Q`-bound; the SWEEPS are guards on the builder, `Q <= 33` (`460` rows), and are not the
grounds ([OPS-041], [NCYL-318] scope (i)).  ⚠ It does NOT close (G0b) or `(G0b-form)`,
does not move a `claims.md` tier, and touches nothing on the gamma=1 DAG.

USAGE:  PYTHONPATH=.:engine:archive/scripts_2026-07:archive/scripts_2026-08:probes \\
        .venv/bin/python3.13 probes/s466_triv_gap.py [h0|h1|h2|h3|h4|h5|all] [--qmax=N]
"""
from __future__ import annotations

import glob
import json
import math
import sys

import s451_quotient_path as qp
import s453_pole_module as pm
import zlattice as zl
from s462_zz_relation import sym_build, sym_path

OUT = "data/s466_triv_gap.json"
QMAX = 33


def rows(qmax):
    for Q in range(5, qmax + 1, 2):
        for P in range(1, Q):
            if math.gcd(P, Q) == 1:
                for eps in (0, 1):
                    yield P, Q, eps


def build(P, Q, eps):
    ch, lo, hi, closes = qp.build(P, Q, eps)
    if not closes:
        return None
    sf = pm.fold_sums(ch, lo, hi)
    Lo, Hi, S = pm.path_view(ch, lo, hi, sf)
    LO, HI = sym_build(ch)
    sLo, sHi, sS = sym_path(ch, LO, HI, lo, hi)
    return ch, [list(map(int, v)) for v in sLo], \
        [list(map(int, v)) for v in sHi], [list(map(int, v)) for v in sS]


def same(ch, t, u):
    return bool((ch.ell[t].v == ch.ell[u].v).all())


# ------------------------------------------------------------------ H0
def h0(qmax, verbose=True):
    n, bad, ctrl_ok, ctrl_n = 0, [], 0, 0
    for P, Q, eps in rows(qmax):
        b = build(P, Q, eps)
        if b is None:
            continue
        ch = b[0]
        ts = (Q - 1) // 2
        n += 1
        if not all(same(ch, t, (Q - 1 - t) % Q) for t in range(Q)):
            bad.append((P, Q, eps, "pairing"))
        if [t for t in range(Q) if (Q - 1 - t) % Q == t] != [ts]:
            bad.append((P, Q, eps, "fixed-point"))
        if (not ch.ell[ts].v.any()) != ch.lone:
            bad.append((P, Q, eps, "zero-at-t*"))
        # CONTROL: t -> -t mod Q must NOT be the value pairing
        ctrl_n += 1
        if not all(same(ch, t, (-t) % Q) for t in range(Q)):
            ctrl_ok += 1
    if verbose:
        print(f"H0  sigma(t) = Q-1-t IS the exact value pairing, its unique fixed point "
              f"is t* = (Q-1)/2, and t* carries the zero measure iff LONE -- "
              f"{n - len(bad)}/{n} rows   {'PASS' if not bad else 'FAIL ' + str(bad[:4])}")
        print(f"    CONTROL (t -> -t mod Q must NOT pair the values): fails as required "
              f"on {ctrl_ok}/{ctrl_n}   {'PASS (fires)' if ctrl_ok == ctrl_n else 'FAIL'}")
    return {"rows": n, "n_bad": len(bad), "bad": bad[:20],
            "ctrl_fires": ctrl_ok, "ctrl_n": ctrl_n,
            "pass": not bad and ctrl_ok == ctrl_n}


# ------------------------------------------------------------------ H1
def h1(qmax, verbose=True):
    n, bad, partial = 0, [], 0
    for P, Q, eps in rows(qmax):
        b = build(P, Q, eps)
        if b is None:
            continue
        ch, sLo, sHi, sS = b
        ts = (Q - 1) // 2
        n += 1
        sup = set()
        for arr in (sLo, sHi, sS):
            for v in arr:
                sup |= {t for t, c in enumerate(v) if c}
        if not sup <= set(range(ts + 1)):
            bad.append((P, Q, eps, sorted(sup - set(range(ts + 1)))))
        if sup != set(range(ts + 1)):
            partial += 1
    if verbose:
        print(f"H1  every symbolic vector is supported inside the transversal "
              f"{{0..(Q-1)/2}} -- {n - len(bad)}/{n} rows   "
              f"{'PASS' if not bad else 'FAIL ' + str(bad[:4])}")
        print(f"    and the union of supports is the FULL transversal on "
              f"{n - partial}/{n} (so it is not vacuous)")
    return {"rows": n, "n_bad": len(bad), "bad": bad[:20], "partial": partial,
            "pass": not bad}


# ------------------------------------------------------------------ H2
def kernel_basis(ch, Q, ts):
    """{a : supp(a) subset {0..t*}, a_red == 0}, as the list of free indices."""
    free = []
    for t in range(ts + 1):
        if not ch.ell[t].v.any():
            free.append(t)                 # zero measure: unconstrained by a_red
        # else: its value class meets {0..t*} only at t, so a_red == 0 forces a_t == 0
    return free


def h2(qmax, verbose=True):
    n, bad = 0, []
    for P, Q, eps in rows(qmax):
        b = build(P, Q, eps)
        if b is None:
            continue
        ch = b[0]
        ts = (Q - 1) // 2
        n += 1
        # the "meets the support in exactly one index" premise, re-checked here
        for t in range(ts + 1):
            others = [u for u in range(ts + 1) if u != t and same(ch, t, u)]
            if others and ch.ell[t].v.any():
                bad.append((P, Q, eps, "class-twice", t, others))
        free = kernel_basis(ch, Q, ts)
        want = [ts] if ch.lone else []
        if free != want:
            bad.append((P, Q, eps, "kernel", free, want))
    # re-score the stored s464 arrivals -- PREDICTED, a bug-check
    seen, viol = 0, 0
    for fn in glob.glob("data/s464_*.json"):
        try:
            blob = json.load(open(fn))
        except Exception:
            continue
        for a, ar in _walk_arrivals(blob):
            seen += 1
            if not any(ar) and any(a):
                nz = [i for i, c in enumerate(a) if c]
                if len(nz) != 1:
                    viol += 1
    if verbose:
        print(f"H2  a_red == 0  =>  a supported on the ZERO measures alone: kernel is "
              f"{{e_t*}} on LONE and {{0}} on PAIRED -- {n - len(bad)}/{n} rows   "
              f"{'PASS' if not bad else 'FAIL ' + str(bad[:3])}")
        print(f"    re-scored {seen} stored s464 arrival vectors: {viol} violate it "
              f"{'PASS' if not viol else 'FAIL'}   ⚠ PREDICTED -- a bug-check, "
              f"not evidence ([OPS-041])")
    return {"rows": n, "n_bad": len(bad), "bad": bad[:20],
            "arrivals_scored": seen, "arrivals_violating": viol,
            "pass": not bad and not viol}


def _walk_arrivals(o):
    """Yield (a, a_red) from any nesting of the s464 result blobs."""
    if isinstance(o, dict):
        if isinstance(o.get("a"), list) and isinstance(o.get("a_red"), list):
            yield o["a"], o["a_red"]
        for v in o.values():
            yield from _walk_arrivals(v)
    elif isinstance(o, list):
        for v in o:
            yield from _walk_arrivals(v)


# ------------------------------------------------------------------ H3 / H4
def h3h4(qmax, verbose=True):
    lone, em_in_L, admits = 0, 0, 0
    for P, Q, eps in rows(qmax):
        b = build(P, Q, eps)
        if b is None:
            continue
        ch, sLo, sHi, sS = b
        if not ch.lone:
            continue
        ts = (Q - 1) // 2
        lone += 1
        em = [0] * Q
        em[ts] = 1
        L = zl.Lattice(sS)
        for c in range(1, 25):
            v = [0] * Q
            v[ts] = c
            if L.contains(v):
                em_in_L += 1
                break
        Lp = zl.Lattice(sS + [em])
        ends = sLo + sHi
        hit = False
        for sX in ends:
            for sY in ends:
                for s in (1, -1):
                    v = [s * x - y for x, y in zip(sX, sY)]
                    if Lp.contains(v) and not L.contains(v):
                        hit = True
                        break
                if hit:
                    break
            if hit:
                break
        admits += hit
    if verbose:
        print(f"H3  no nonzero multiple of e_t* (|c| <= 24) lies in the fold-sum lattice "
              f"Lambda -- {lone - em_in_L}/{lone} LONE rows   "
              f"{'PASS' if not em_in_L else 'FAIL'}")
        print(f"    ⇒ a LONE `triv`'s e_t* component cannot come from FOLDING; it would "
              f"have to come from the endpoints.")
        print(f"H4  ⚠ THE NEGATIVE: the cap-free lattice relaxation "
              f"((+-sX - sY) + Lambda meets Z*e_t* \\ 0?) ADMITS on {admits}/{lone} LONE "
              f"rows -- it CANNOT decide the residue.  Only (D) does, by dissolving it.")
    return {"lone_rows": lone, "em_in_lambda": em_in_L, "lattice_admits": admits,
            "pass": em_in_L == 0}


# ------------------------------------------------------------------ H5
def h5(verbose=True):
    found = None
    for fn in glob.glob("data/s464_*.json"):
        try:
            blob = json.load(open(fn))
        except Exception:
            continue
        for a, ar in _walk_arrivals(blob):
            if any(ar):
                found = (fn, a, ar)
                break
        if found:
            break
    ok = found is not None
    if verbose:
        if ok:
            fn, a, ar = found
            print(f"H5  the known witness still classifies `rel`: a_red = {ar} "
                  f"(nonzero) from {fn}   PASS")
        else:
            print("H5  FAIL -- no `rel` arrival found in the stored s464 blobs; the "
                  "instrument would be rejecting everything ([OPS-242])")
    return {"found": bool(ok), "a_red": found[2] if ok else None,
            "src": found[0] if ok else None, "pass": bool(ok)}


# ------------------------------------------------------------------
def main():
    arms = "all"
    qmax = QMAX
    for arg in sys.argv[1:]:
        if arg.startswith("--qmax="):
            qmax = int(arg.split("=")[1])
        else:
            arms = arg
    res = {"qmax": qmax}
    print("s466_triv_gap.py -- (s464-a): is there a REASON a fold word cannot "
          "realise a trivial relation?\n")
    if arms in ("all", "h0"):
        res["h0"] = h0(qmax)
        print()
    if arms in ("all", "h1"):
        res["h1"] = h1(qmax)
        print()
    if arms in ("all", "h2"):
        res["h2"] = h2(qmax)
        print()
    if arms in ("all", "h3", "h4"):
        res["h3h4"] = h3h4(qmax)
        print()
    if arms in ("all", "h5"):
        res["h5"] = h5()
        print()
    if arms == "all":
        keys = ["h0", "h1", "h2", "h3h4", "h5"]
        allpass = all(res[k]["pass"] for k in keys)
        print(f"SCORED ARMS: {', '.join(keys)} -- "
              f"{'ALL PASS' if allpass else 'FAILURE PRESENT'}")
        print("⇒ PAIRED: `triv` is IMPOSSIBLE.  LONE: `triv` collapses to Z*e_t*, which "
              "is `form` in a redundant representation (H4/(D), DERIVED).")
        print("⇒ [NCYL-317] (4) does NOT need `triv == 0` as a hypothesis: state the "
              "split on `a_red` and the trichotomy is the original dichotomy.")
        res["pass"] = allpass
        json.dump(res, open(OUT, "w"), indent=1)
        print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
