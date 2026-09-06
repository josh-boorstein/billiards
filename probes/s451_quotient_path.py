#!/usr/bin/env python3.13
"""s451_quotient_path.py -- PRE-REGISTERED.  Queue item (w): the `ell` ORDERING, and what
it turns `p = 0` into.

PRIOR ART: `rulings.py --grep` on 'nested interval', 'rotation coding', 'sturmian' ->
NOTHING; 'three-distance' -> [WFLOOR-034] (the mechanism on the graze orbit `{kP/Q}`,
RUNG cells, F-leg strand), [SEAM-006] (do NOT pose the seam half as EQUIDISTRIBUTION),
[EXPO-001] (the three-distance model for `lambda(P,r)` is REFUTED at `P=13, r=1`);
'palindrome' -> nothing on this object; 'quotient' -> [OPS-214] (push a topological count
to the genus-0 quotient), nothing on an `iota` quotient of `B`.  The necklace, its `+P`
gluing order and its coordinates are [NCYL-291]; the per-kite flat geometry, the tracer
and the `k <-> -k` pairing are [NCYL-292]; the bijection form of `p = 0` is [NCYL-287];
`E = h + f + 2p` is [OPS-218].  ⚠ NO probe in the repo puts a GLOBAL coordinate on the
chain or quotients it by `iota`; and the F-leg strand's three-distance machinery
([WFLOOR-034]) has never been pointed at the NCYL strand.
=> what is NEW here is the global coordinate, the `iota`-quotient PATH, the closed forms
for its data, and the two negatives in H4/H5.  Everything else is prior art and is USED.

------------------------------------------------------------------------------------------
THE GLOBAL COORDINATE (this is what makes the rest possible).

[NCYL-292]'s per-kite maps move the transverse coordinate by `+-(M-m)` on every THROUGH
crossing, which is why a trajectory reads as structureless.  But those maps say exactly
that kite `k`'s two interfaces are ALIGNED -- at their TOP ends if the kite is O-sector
(`x -> x + (M-m)`), at their BOTTOM ends if A-sector (`x -> x`).  So put every interface
on ONE line:

    I_t = [u_t, u_t + ell_t],     u_near = u_far + (M-m)  (O-sector) / u_far  (A-sector).

Then (i) consecutive interfaces are NESTED (they share an endpoint and differ in length --
automatic, not a finding); (ii) the through band is exactly `I_{t-1} ∩ I_t` = the smaller;
(iii) the fold is the reflection of the larger's protruding tail about its own midpoint;
and (iv) ⇒⇒ **THE LEAF'S COORDINATE `y` IS CONSTANT BETWEEN FOLDS.**  A leaf traverses the
maximal arc of `S(y) = {t : y ∈ I_t}` containing its start, then reflects and reverses.
⚠ `u` is only well defined if it CLOSES around the cycle; H2 checks that exactly.

THE `iota` QUOTIENT (the system halves, and the target becomes an ESCAPE statement).

`ell_{-1-t} = ell_t` is algebra (`c_{-1-t} = -c_t mod 2Q`), and H2 gives `u_{-1-t} = u_t`,
so `I_{rho(t)} = I_t` for `rho(t) = -1-t` -- the interval sequence around the cycle is a
PALINDROME.  `rho` is the reflection of `Z/Q` with axis points interface `m` and the kite-0
gap, i.e. EXACTLY [NCYL-283]'s two `Fix(iota)` components.  Quotienting:

    a PATH  J_0 -k_0- J_1 -k_1- ... -k_{m-1}- J_m,   J_i := I_{m-i} = I_{m+i},
    one interior pole per node `k_i`, and the two ENDS are the two Fix components.

  ⇒⇒ `p = 0`  <=>  EVERY interior pole of the path escapes to an END of it.
  Consistency: [NCYL-287] says `C = #escapes + 1`, and `C(paired) = ceil(Q/2) = m + 1`
  ([NCYL-241]), so all `m` quotient poles must escape.  The counts agree -- but that is a
  CONSISTENCY CHECK, not evidence: both sides are the same count ([OPS-041]).

THE CLOSED FORMS -- and this is the CF handle (w) was opened for.

  CF1  |J_i| = |cos(i*P*pi/Q)|  (paired) / |sin(i*P*pi/Q)|  (lone),  i = 0..m.
       Equivalently `|J_i| = sin(pi * ||i*beta - 1/2||)` with `beta = P/Q` and `||.||` the
       distance to the nearest integer (paired): the node is LARGE where the rotation
       orbit `{i*beta}` is near `1/2` and SMALL where it is near `0`.
  CF2  `k_i` is TOP-aligned  <=>  frac((i+1)*beta + 1/2) < beta   (paired)
                              <=>  frac((i+1)*beta)       < beta   (lone),
       i.e. `(i+1)*beta` lands in the length-`beta` window at `1/2`.
  ⇒ BOTH the node sizes and the alignment coding are codings of the SAME rotation orbit
    `{i*P/Q}` against the point `1/2`, which is the `iota` axis.  The ascent/descent coding
    (`|J_{i+1}| > |J_i|`) is then a Sturmian coding of slope `beta` by an interval of
    length `1/2`.  ⚠⚠ CF1 and CF2 are DERIVED (from `c_{m-i} = c_m - 2Pi`, `c_m = Q`/`0`,
    and `osec`'s definition) and then checked; H3 confirms MY ALGEBRA and is NOT evidence
    about the dynamics ([OPS-222]).

  ⇒ THE GRAZE-ORBIT LINK.  `{i*P/Q}` is the orbit [WFLOOR-034] hangs three-distance
    structure on in the F-leg strand.  It is a HANDLE, not a result.  ⚠⚠ AND NOT A NEW
    SIGHTING OF THE ORBIT: [NCYL-291]'s headline already says `Delta` IS A ROTATION BY
    `P/Q`, and `osec` was already an interval coding of it.  What is new is the
    RE-CENTRING on the `iota` axis `1/2` and the FIRST-ENTRY form in H6 -- not the CF
    connection itself ([OPS-225]).

HYPOTHESES.  ⚠ Read [OPS-041]/[OPS-222] first: an arm that cannot fail is not evidence.
  H1  ⇒⇒ THE FAITHFULNESS GATE.  The global-`u` tracer must reproduce `Chain.trace`'s end
      label on every interior pole prong.  ⚠⚠ IT IMPORTS `Chain` AND READS
      `far/near/osec/delta/ell` FROM IT, SO IT SHARES THE MODEL -- the tell is an import
      ([OPS-222]).  It tests the RECOORDINATIZATION, not the model, and that is worth
      running precisely because a recoordinatization is the kind of thing that breaks
      silently.  ⚠ Also NOT evidence: the ends coming out at the MIRROR pole `R_{-k}` --
      that is [NCYL-292] H5, which IS `p = 0` written out ([NCYL-292] DO-NOT (1)).
  H2  ⇒⇒ CLOSURE, and it is load-bearing: walking the alignment rule around the cycle must
      return `u_0`.  If it fails there is no global coordinate and everything above dies.
      Exact (`Coord.eq`), no tolerance.
  H3  THE CLOSED FORMS CF1 + CF2.  ⚠ A DERIVATION CHECK, NOT A TEST (see above).
  H4  ⇒⇒ THE FIRST GENUINE NEGATIVE, AND IT COULD HAVE COME OUT AT `100%`: does every pole
      escape on its FIRST run -- i.e. is the escape criterion STATIC (`S(y)`'s arc through
      the start is `rho`-invariant), with no dynamics at all?  A `100%` here would have
      reduced `p = 0` to a statement about one subset of `Z/Q`.
  H5  ⇒⇒ THE SECOND, AND IT COULD HAVE COME OUT BOUNDED: the distribution of FOLD COUNTS to
      escape.  A uniform bound (or a bound in `Q`) would license a finite-step induction of
      the kind [NCYL-292]'s one-kite partial starts.
  H6  ⇒⇒ THE SHARPENING, AND IT IS THE FORM THE RESIDUAL SHOULD BE ATTACKED IN: run the
      SAME dynamics from CF DATA ALONE -- state `(i, d)`, inputs only `L[i]` (CF1) and
      `align[i]` (CF2) -- and score it against `qtrace` under the HALF-TRAJECTORY rule
      (`q_folds == 2*r_folds` and the full trajectory ends at the mirror pole `R_{-k}`).
      ⚠ **THIS ONE CAN AND DID FAIL: a first draft that advanced `i` but not `d` on a
      through crossing read `534/638`.**  It also tests the quotient picture itself --
      the factor `2` is the claim that the full trajectory is the quotient one doubled
      about its Fix crossing.
  ⚠ NOT A HYPOTHESIS: "every fold count is EVEN".  Escape + the `rho` symmetry make the
    trajectory a palindrome about its Fix crossing, so evenness is forced by H1 and cannot
    fail ([OPS-041]).  It is printed as a consistency line only.

RESULTS (s451; `all --qmax=31`, `423.7 s`, `data/s451_quotient.json`).
  H1  `1900/1900` end labels agree, `0` disagree, `12` capped (`Q <= 19`).  ⚠ `1900/1900`
      also end at the mirror pole `R_{-k}` -- NOT evidence, see H1 above.
  H2  ⇒⇒ `420/420` rows close EXACTLY (`Q = 5..31`, all `P`, both `eps`, both classes).
      The global coordinate exists.
  H3  CF1 `4812/4812`, CF2 `4392/4392`.  Derivation check.
  H4  ⇒⇒ FIRST-RUN ESCAPE IS `1064/2196 = 48.45%` (paired) AND `450/1986 = 22.66%` (lone).
      ⇒ THERE IS NO STATIC ESCAPE CRITERION.  `p = 0` is irreducibly dynamical, and the
      lone class is the harder half by a factor of two.
  H5  ⇒⇒ FOLDS-TO-ESCAPE ARE UNBOUNDED: max `40842` (paired), `77194` (lone), at a
      `200000`-step cap, still climbing with `Q`.  `never-escaped = 0`; `12`/`76` capped.
      ⇒ NO FINITE-STEP INDUCTION CAN WORK, which closes off the natural continuation of
      [NCYL-292]'s one-kite partial.
  H6  ⇒⇒ `638/638` (paired, `Q <= 21`) -- THE DYNAMICS IS FULLY DETERMINED BY THE
      ROTATION ORBIT `{i*beta}`: `L[i] = sin(pi*||i*beta - 1/2||)` and the Sturmian
      alignment coding, nothing else.  ⇒ WITHIN a maximal run of equal alignment one
      endpoint is shared, so `d` is CONSTANT across the run and the leaf passes node
      after node while `L[i] > d` -- i.e. the run length is the FIRST-ENTRY TIME of
      `{i*beta}` into the window `||i*beta - 1/2|| <= arcsin(d)/pi`, which is exactly
      what three-distance / the Ostrowski expansion of `beta` controls.  ⚠ **AND THE GAP
      IS THE COMPOSITION, NOT THE STEP: each fold resets `d -> (L[i] +- L[nb]) - d`, so
      the next run is a first-entry with a DIFFERENT window.  `p = 0` is *this sequence
      of first-entry problems runs off an end*.**  Three-distance supplies one step.
  Instrument audit: worst `Chain.safety()` = `405` on scored poles, `239` on capped -- far
  above `1`, so the float band gate is sound here.  ⚠ Reset `ch.margin`/`ch.maxcoef` PER
  TRACE; aggregating them per ROW pairs one trace's margin with another's coefficient and
  reads `0.2` on this same data ([OPS-224]).

  python3.13 probes/s451_quotient_path.py [all|gate] [--qmax=N]
"""
from __future__ import annotations

import collections
import json
import math
import sys
import time

import s450_chain_maps as cm

TOL = 1e-9
OUT = "data/s451_quotient.json"


# ---------------------------------------------------------------- the global coordinate

def build(P, Q, eps):
    """`Chain` + the global interval endpoints `lo[t] <= hi[t]`, and the closure verdict."""
    ch = cm.Chain(P, Q, eps)
    u = {0: ch.zero}
    for k in range(1, Q):                     # kites 1..Q-1 determine every u_t from u_0
        far, near, d = ch.far[k], ch.near[k], ch.delta[k]
        if far in u:
            u[near] = u[far] + d if ch.osec[k] else u[far]
        else:
            u[far] = u[near] - d if ch.osec[k] else u[near]
    far, near, d = ch.far[0], ch.near[0], ch.delta[0]          # H2: kite 0 must agree
    closes = (u[far] + d if ch.osec[0] else u[far]).eq(u[near])
    lo = {t: u[t] for t in range(Q)}
    hi = {t: u[t] + ch.ell[t] for t in range(Q)}
    return ch, lo, hi, closes


def pole_y(ch, lo, hi, k):
    """The prong of `R_k`: the midpoint of the tail `I_far \\ I_near`, and the direction
    AWAY from kite `k` (the fold's fixed point emits one prong, on the far interface)."""
    far, near = ch.far[k], ch.near[k]
    bot = lo[far].eq(lo[near])                       # bottom-aligned -> tail sits on top
    y = cm.half((hi[near] + hi[far]) if bot else (lo[far] + lo[near]))
    return y, (-1 if far == (k - 1) % ch.Q else +1), far


def qtrace(ch, lo, hi, k, cap=200_000):
    """Follow `R_k`'s prong in the GLOBAL coordinate.  `y` is constant between folds.

    Returns (end_label, folds, crossed_Fix, steps).  `crossed_Fix` is True as soon as the
    leaf reaches interface `m` or transits the tie kite 0 -- the two `Fix(iota)` arcs.
    """
    Q, m = ch.Q, ch.m_idx
    imv = ch.ring.imvals
    y, d, t = pole_y(ch, lo, hi, k)
    folds, crossed, steps = 0, (t == m), 0
    while True:
        steps += 1
        if steps > cap:
            return ("CAP", None), folds, crossed, steps
        if not steps & 511:
            # ⚠ REQUIRED, not hygiene: `y` is rebuilt by `s - y` at every fold, so its
            # float drifts linearly in the step count.  Recompute it from the EXACT
            # vector and audit the integer magnitude, exactly as `Chain.trace` does --
            # without this the `<`/`>` band test below is unsound past ~1e5 steps.
            y.f = float(y.v @ imv)
            a = int(abs(y.v).max())
            if a > ch.maxcoef:
                ch.maxcoef = a
        t2 = (t + d) % Q
        kk = t2 if d == 1 else (t2 + 1) % Q
        # `Chain.cmp` is the AUDITED comparison: float-gated, exact on ties, and it
        # records `ch.margin`.  `ch.safety()` then prices the gate ([NCYL-292]).
        cl, chg = ch.cmp(y, lo[t2]), ch.cmp(y, hi[t2])
        if cl == 0:
            return ("ZO", kk), folds, crossed, steps
        if chg == 0:
            return ("ZA", kk), folds, crossed, steps
        if cl > 0 and chg < 0:                                       # through
            if (d == 1 and t == Q - 1) or (d == -1 and t == 0):
                crossed = True                                       # transited kite 0
            t = t2
            if t == m:
                crossed = True
            continue
        s = (hi[t2] + hi[t]) if lo[t].eq(lo[t2]) else (lo[t] + lo[t2])   # the fold's tail
        if y.eq(cm.half(s)):
            return ("R", kk), folds, crossed, steps
        y = s - y
        d = -d
        folds += 1


def first_run_escapes(ch, lo, hi, k):
    """H4: does `R_k`'s prong reach a `Fix` arc BEFORE its first fold?  Purely static --
    it asks whether `S(y)`'s arc through the start is `rho`-invariant."""
    Q, m = ch.Q, ch.m_idx
    y, d, t = pole_y(ch, lo, hi, k)
    crossed, steps = (t == m), 0
    while steps <= 2 * Q:
        steps += 1
        t2 = (t + d) % Q
        if not (ch.cmp(y, lo[t2]) > 0 and ch.cmp(y, hi[t2]) < 0):
            return crossed
        if (d == 1 and t == Q - 1) or (d == -1 and t == 0):
            crossed = True
        t = t2
        if t == m:
            crossed = True
    return crossed


def reduced(ch, m, L, align, i0, d0, dirn, cap=200_000):
    """H6 -- THE SAME DYNAMICS FROM CF DATA ALONE.

    State `(i, d)`: `i` = node index in the quotient path, `d` = height above `J_i`'s
    BOTTOM end, plus a direction.  Inputs are ONLY `L[i] = |J_i|` (CF1) and `align[i]`
    (CF2, computed from `P/Q` arithmetic) -- no `lo`/`hi`, no `far`/`near`, no `osec`.

    Crossing node `j` from `J_i` to the neighbour `J_nb`:
        BOTTOM-aligned  d' = d                     (bottoms agree)
        TOP-aligned     d' = d - L[i] + L[nb]      (tops agree)
      through iff `0 < d' < L[nb]`; otherwise FOLD inside `J_i`,
        d -> (L[i] + L[nb]) - d   (BOTTOM)   /   (L[i] - L[nb]) - d   (TOP),
      and the direction flips.  Poles are the fold fixed points; ESCAPE is running off
      either END of the path (`j < 0` or `j > m-1`), the two `Fix(iota)` components.

    ⚠ THE BUG THIS COST: on a THROUGH crossing you must advance `d` as well as `i`.
    Advancing only `i` reads `534/638` -- close enough to look like a near-miss model
    rather than a one-line bookkeeping error.
    """
    i, d, folds, steps, crossed = i0, d0, 0, 0, (i0 == 0)
    while True:
        steps += 1
        if steps > cap:
            return "CAP", folds, crossed
        j = i if dirn > 0 else i - 1
        if j < 0 or j > m - 1:
            return "ESC", folds, True
        nb = i + dirn
        d2 = (d - L[i] + L[nb]) if align[j] else d
        c0, cL = ch.cmp(d2, ch.zero), ch.cmp(d2, L[nb])
        if c0 == 0 or cL == 0:
            return "Z", folds, crossed
        if c0 > 0 and cL < 0:
            i, d = nb, d2
            if i == 0:
                crossed = True
            continue
        s = (L[i] - L[nb]) if align[j] else (L[i] + L[nb])
        if ch.cmp(d, cm.half(s)) == 0:
            return "R", folds, crossed
        d = s - d
        dirn = -dirn
        folds += 1


def h6(qmax):
    """H6 scored against `qtrace` under the HALF-TRAJECTORY rule."""
    ok = bad = 0
    ex = []
    for P, Q, eps in centres(qmax):
        ch, lo, hi, _ = build(P, Q, eps)
        if ch.lone:
            continue
        m = ch.m_idx
        L = [ch.ell[m - i] for i in range(m + 1)]
        align = [(((i + 1) * P / Q + 0.5) % 1.0) < P / Q - 1e-12 for i in range(m)]
        for k in range(1, m + 1):
            if ch.tie[k] or ch.ell[ch.near[k]].is_zero():
                continue
            j = m - k
            big, small = (j, j + 1) if ch.cmp(L[j], L[j + 1]) > 0 else (j + 1, j)
            s = (L[big] - L[small]) if align[j] else (L[big] + L[small])
            r_end, r_folds, _ = reduced(ch, m, L, align, big, cm.half(s),
                                        -1 if big == j else +1)
            q_end, q_folds, _, _ = qtrace(ch, lo, hi, k)
            if r_end == "CAP" or q_end[0] == "CAP":
                continue
            good = (r_end == "ESC" and q_end[0] == "R"
                    and q_end[1] == (-k) % ch.Q and q_folds == 2 * r_folds)
            if good:
                ok += 1
            else:
                bad += 1
                if len(ex) < 5:
                    ex.append([P, Q, eps, k, r_end, r_folds, str(q_end), q_folds])
    return {"ok": ok, "bad": bad, "examples": ex}


# ---------------------------------------------------------------- hypotheses

def centres(qmax):
    for Q in range(5, qmax + 1, 2):
        for P in range(1, Q):
            if math.gcd(P, Q) == 1:
                for eps in (0, 1):
                    yield P, Q, eps


def h2_h3_h4_h5(qmax):
    res = {"h2_rows": 0, "h2_fail": 0, "h2_examples": [],
           "h3_len": [0, 0], "h3_align": [0, 0], "h3_examples": [],
           "h4": {}, "h5": {}, "even": [0, 0]}
    for cls in ("paired", "lone"):
        res["h4"][cls] = {"poles": 0, "first_run": 0}
        res["h5"][cls] = {"poles": 0, "capped": 0, "never": 0,
                          "hist": collections.Counter(), "max": 0}
    res["safety_scored"] = float("inf")
    res["safety_capped"] = float("inf")
    for P, Q, eps in centres(qmax):
        ch, lo, hi, closes = build(P, Q, eps)
        m, cls = ch.m_idx, ("lone" if ch.lone else "paired")
        res["h2_rows"] += 1
        if not closes:
            res["h2_fail"] += 1
            if len(res["h2_examples"]) < 5:
                res["h2_examples"].append([P, Q, eps])
        # H3 -- CF1 (node lengths) and CF2 (alignment coding)
        for i in range(m + 1):
            want = (abs(math.cos(i * P * math.pi / Q)) if cls == "paired"
                    else abs(math.sin(i * P * math.pi / Q)))
            ok = abs(ch.ell[m - i].f - want) <= 1e-11
            res["h3_len"][0 if ok else 1] += 1
        for i in range(m):
            sh = 0.5 if cls == "paired" else 0.0
            want = ((i + 1) * P / Q + sh) % 1.0 < P / Q - 1e-12
            ok = bool(ch.osec[m - i]) == want
            res["h3_align"][0 if ok else 1] += 1
            if not ok and len(res["h3_examples"]) < 5:
                res["h3_examples"].append([P, Q, eps, i])
        # H4 / H5 -- one representative per iota-PAIR of poles (k and -k share y)
        for k in range(1, m + 1):
            if ch.tie[k] or ch.ell[ch.near[k]].is_zero():
                continue                                  # lone class degenerate flanks
            res["h4"][cls]["poles"] += 1
            res["h4"][cls]["first_run"] += first_run_escapes(ch, lo, hi, k)
            ch.margin, ch.maxcoef = 1.0, 0
            end, folds, crossed, _ = qtrace(ch, lo, hi, k)
            s5 = res["h5"][cls]
            s5["poles"] += 1
            if end[0] == "CAP":
                s5["capped"] += 1
                res["safety_capped"] = min(res["safety_capped"], ch.safety())
                continue
            res["safety_scored"] = min(res["safety_scored"], ch.safety())
            s5["hist"][folds] += 1
            s5["max"] = max(s5["max"], folds)
            s5["never"] += (not crossed)
            res["even"][0 if folds % 2 == 0 else 1] += 1
    for cls in ("paired", "lone"):
        res["h5"][cls]["hist"] = dict(sorted(res["h5"][cls]["hist"].items()))
    return res


def h1_gate(qmax):
    """H1 -- end labels against `Chain.trace`.  SLOW: it runs s450's tracer."""
    agree = disagree = capped = mirror = scored = 0
    ex = []
    for P, Q, eps in centres(qmax):
        ch, lo, hi, _ = build(P, Q, eps)
        for k in range(1, Q):
            if ch.tie[k] or ch.ell[ch.near[k]].is_zero():
                continue
            fx, j = ch.ffix[k], ch.far[k]
            ref = ch.trace(j, cm.Coord(fx.v.copy(), fx.f), ch.next_kite(k, j))[0]
            got, _, _, _ = qtrace(ch, lo, hi, k)
            if ref == "CAP" or got[0] == "CAP":
                capped += 1
                continue
            if ref == got:
                agree += 1
            else:
                disagree += 1
                if len(ex) < 5:
                    ex.append([P, Q, eps, k, str(ref), str(got)])
            scored += 1
            mirror += (got[0] == "R" and got[1] == (-k) % Q)
    return {"agree": agree, "disagree": disagree, "capped": capped,
            "mirror": mirror, "scored": scored, "examples": ex}


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "all"
    qmax = 31
    for a in sys.argv[2:]:
        if a.startswith("--qmax="):
            qmax = int(a.split("=")[1])
    t0 = time.time()
    out = {"qmax": qmax, "mode": mode}
    if mode in ("all", "gate"):
        out["h1"] = h1_gate(min(qmax, 19))
        g = out["h1"]
        print(f"H1 faithfulness  agree={g['agree']}/{g['scored']}  disagree={g['disagree']}"
              f"  capped={g['capped']}  mirror(NOT evidence)={g['mirror']}/{g['scored']}")
        if g["examples"]:
            print("   ", g["examples"])
    if mode in ("all", "h6"):
        out["h6"] = h6(min(qmax, 21))
        print(f"H6 CF-only dynamics vs qtrace (half-trajectory rule): "
              f"{out['h6']['ok']}/{out['h6']['ok']+out['h6']['bad']}   "
              f"{out['h6']['examples']}")
    if mode == "all":
        r = h2_h3_h4_h5(qmax)
        out.update(r)
        print(f"H2 closure       {r['h2_rows']-r['h2_fail']}/{r['h2_rows']} rows close "
              f"EXACTLY   fails={r['h2_examples']}")
        print(f"H3 CF1 lengths   {r['h3_len'][0]}/{sum(r['h3_len'])}   "
              f"CF2 alignment {r['h3_align'][0]}/{sum(r['h3_align'])}   "
              f"(DERIVATION CHECK, not evidence)  {r['h3_examples']}")
        for cls in ("paired", "lone"):
            a, b = r["h4"][cls]["first_run"], r["h4"][cls]["poles"]
            s5 = r["h5"][cls]
            print(f"H4 {cls:6s} first-run escape {a}/{b} "
                  f"({100*a/max(b,1):.2f}%)  -- a STATIC criterion would need 100%")
            print(f"H5 {cls:6s} folds-to-escape  max={s5['max']}  "
                  f"never-escaped={s5['never']}  capped={s5['capped']}/{s5['poles']}")
        print(f"   consistency (forced, not evidence): even folds "
              f"{r['even'][0]}/{sum(r['even'])}")
        print(f"   INSTRUMENT AUDIT  worst Chain.safety() on SCORED poles = "
              f"{r['safety_scored']:.3g}   (on capped: {r['safety_capped']:.3g});"
              f"  < 1 means a band test could have branched wrongly")
    out["elapsed_s"] = round(time.time() - t0, 1)
    with open(OUT, "w") as f:
        json.dump(out, f, indent=1, default=str)
    print(f"-> {OUT}   [{out['elapsed_s']}s]")


if __name__ == "__main__":
    main()
