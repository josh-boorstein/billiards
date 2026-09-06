#!/usr/bin/env python3.13
"""s455_closure_proof.py -- PRE-REGISTERED.  Queue item (i): PROVE that [NCYL-294]'s
global coordinate CLOSES.  It does, and the proof is four lines of trigonometry.

PRIOR ART: `rulings.py --grep` on 'global coordinate' -> [NCYL-294] (H2, THE TARGET: the
closure is VERIFIED `420/420`, `Q = 5..31`, NOT proved -- and s454's amendment names
proving it "the highest-value open item on this strand"), [NCYL-292] (the per-kite flat
geometry and the `+-(M-m)` through-shift the alignment rule is read off), [NCYL-291] (the
necklace, its `+P` gluing order, `Delta` = rotation by `P/Q`); 'osec' -> [NCYL-294],
[OPS-225]; 'sector rule' -> [OPS-223]; 'closure' / 'telescop' / 'sum-to-product' /
'odd symmetry' -> NOTHING on this object.  ⚠ [OPS-225] (grepping NEW vocabulary will not
find prior art in the OLD) applied against the PARENT [NCYL-294]: its H2 IS this
statement, it is open by its own words, and no ruling anywhere attempts it.  The control
design is [OPS-219]'s (wrong-parameter controls in the SAME pass, read FIRST) and its
warning is honoured -- `sector="complement"` is VACUOUS and is NOT used here.
=> what is NEW is the PROOF; the object, the model and the verification are prior art.

------------------------------------------------------------------------------------------
WHAT CLOSURE IS (reduced to one identity -- this is Step 0 and it is the whole trick).

`s451_quotient_path.build` lays every interface `t` on one line as `I_t = [u_t, u_t+ell_t]`
by walking kites `1..Q-1`, then asks whether kite `0` agrees.  Kite `k` joins interfaces
`lo = k-1` and `hi = k` and imposes

    u[near] - u[far] = delta[k]  (osec[k])   /   0  (otherwise),    delta[k] = ell[far]-ell[near].

⇒⇒ **BOTH BRANCHES ARE SYMMETRIC IN far/near**: the first says `u[hi]+ell[hi] = u[lo]+ell[lo]`
(align at the TOP ends), the second says `u[hi] = u[lo]` (align at the BOTTOM ends).  So the
far/near split -- the only place the ORDER of the two lengths enters -- CANCELS, and with
`g[k] := u[k]-u[k-1]` the cycle telescopes.  Closure is therefore exactly

    ⇒⇒  S := sum_{k in Z/Q}  osec[k] * ( ell[k-1] - ell[k] )  ==  0.

THE THEOREM.  `Q` odd, `0 < P < Q`, `gcd(P,Q) = 1`, `eps in {0,1}`; `c_t = (c_0 + 2Pt) mod 2Q`
with `c_0 = P` (eps=0) / `P+Q` (eps=1); `ell_t = sin(pi c_t / 2Q)`; `osec[k] = 1` iff
`(-c_{k-1}) mod 2Q` lies in `(0, 2P)`.  Then `S = 0`.  ⇒ THE GLOBAL COORDINATE EXISTS.

PROOF.
 (1) REINDEX by `t = k-1`:  `S = sum_t w_t (ell_t - ell_{t+1})`,  `w_t = [(-c_t) mod 2Q in (0,2P)]`.
 (2) WINDOW.  `c_t in [0,2Q)`, so `w_t = 1  <=>  2Q - 2P < c_t < 2Q`.
 (3) SUCCESSOR, INCLUDING THE WRAP.  `c_{t+1} = (c_t + 2P) mod 2Q` for EVERY `t in Z/Q`;
     at `t = Q-1` this is `c_0` because `2PQ = 0 mod 2Q`.
 (4) ⇒⇒ THE MOD-`2Q` REDUCTION BITES EXACTLY ON THE SUPPORT, AND IT FLIPS A SIGN.  If
     `w_t = 1` then `c_t + 2P in (2Q, 2Q+2P)`, so `c_{t+1} = c_t + 2P - 2Q` and
     `ell_{t+1} = sin(pi(c_t+2P)/2Q - pi) = -sin(pi(c_t+2P)/2Q)`.  Hence the summand is a
     SUM, not a difference, and sum-to-product gives
         `ell_t - ell_{t+1} = 2 cos(pi P/2Q) * sin(pi (c_t + P)/2Q)`.
     ⚠ This is the step that makes the identity true: the window half-width `2P` and the
     step `2P` are THE SAME NUMBER, so the support is precisely the set that wraps.
 (5) THE SUPPORT AS A SET.  `gcd(2P,2Q) = 2`, so `{c_t} = {x in [0,2Q) : x = c_0 mod 2}` and
     the support is `{2Q - r : 0 < r < 2P, r = c_0 mod 2}`; there
     `sin(pi(c_t+P)/2Q) = sin(pi + pi(P-r)/2Q) = -sin(pi(P-r)/2Q)`.
 (6) ⇒⇒ ODD SYMMETRY.  With `s = P - r` the index set is `{s : |s| <= P-1, s = P-c_0 mod 2}`,
     which is SYMMETRIC UNDER `s -> -s` (a parity class is), and `sin(pi s/2Q)` is ODD.
     The sum vanishes termwise-in-pairs.  ∎

⇒ SCOPE, and it is WIDER than the model needs: steps (1)-(6) never use `Q` ODD and never
  use which `eps` produced `c_0` -- only `gcd(P,Q)=1`, `0<P<Q`, and `c_0` fixed.  H5 checks
  that stated generality rather than leaving it as a claim.  ⚠ It is NOT vacuous width: the
  identity is FALSE under every perturbation of the sector rule or the window (H4).

HYPOTHESES.  ⚠ [OPS-041]/[OPS-219]: the controls are in the same pass and are read FIRST.
  H1  ⇒⇒ THE REDUCTION (Step 0): `build`'s `closes` verdict `<=>` `S == 0`, exact.
      ⚠⚠ **ON REAL ROWS THIS IS VACUOUS -- both sides are True on all `420`.**  It is
      therefore scored on the PERTURBED models too, where both sides can and do fail; only
      the both-False rows make the equivalence informative.
  H2  THE PROOF STEPS, each checked on its own: far/near symmetry of the kite constraint,
      the reindex, the successor-with-wrap, the window, the sum-to-product form, the
      s-set identification, and its negation symmetry.  ⚠ A DERIVATION CHECK of MY ALGEBRA,
      not evidence about the model ([OPS-222]).
  H3  COVERAGE: `S == 0` ring-exact at `Q` FAR beyond [NCYL-294]'s `31`.  ⚠ This is a
      PREDICTION OF THE PROOF, so it cannot be evidence FOR it; what it tests is that the
      proof describes the CODE (a mis-transcribed model would diverge at some `Q`).
  H4  ⇒⇒ THE ARM THAT CAN FAIL, AND THE ONLY ONE: the identity must BREAK under
      wrong-parameter controls -- `sector="shift"` (same O-count, kites moved between the
      classes), `sector="parity"` (wrong count), `order != P` (wrong gluing step), and the
      window half-width `2P -> 2(P+-1)`.  A control scoring `420/420` would mean `S = 0` is
      bookkeeping and the proof explains nothing.  ⚠ `sector="complement"` is EXCLUDED: it
      is the global relabelling `Z_O <-> Z_A` and measures nothing ([OPS-219]).
  H5  THE SCOPE CLAIM: `S == 0` also at EVEN `Q` and for arbitrary `c_0` of either parity.
      Out of the model's scope by [NCYL-291] (the necklace has one `iota`-fixed kite only
      at odd `Q`) -- this scores the ALGEBRA's stated generality, nothing about billiards.

RESULTS (s455; `all --qmax=201`, `26.4 s`, `data/s455_closure.json`, `logs/s455_closure.log`).
  H4  ⇒⇒ READ FIRST ([OPS-219]).  EVERY CONTROL BREAKS IT, `Q = 5..31`: `sector=shift`
      **`14/420`**, `sector=parity` **`0/420`**, `order=P+2` **`10/316`**, window `W=P+1`
      **`0/392`**, window `W=P-1` **`14/392`** -- against **`420/420`** for the true model.
      ⇒ `S = 0` IS SPECIFIC TO THE SECTOR RULE, NOT FORCED BY THE BOOKKEEPING.
  H1  ⇒⇒ `closes <=> S==0` on **`1260/1260`** rows across the true model and two perturbed
      ones, with **`826` BOTH-FALSE** rows -- the equivalence is informative in BOTH
      directions.  (True model alone: `420/420`, all both-True, i.e. VACUOUS on its own.)
  H2  every step `840/840` (`Q = 5..45`, all `P`, both `eps`): far/near symmetry, reindex,
      successor, window, sum-to-product, s-set, negation symmetry.
  H3  ⇒⇒ `S == 0` EXACTLY on **`16560/16560`** rows, `Q = 5..201` odd, all `P`, both `eps`
      -- `[NCYL-294]`'s `420` extended **`39x`** with zero failures, in `25 s`.
  H5  `S == 0` on **`1324/1324`** even-`Q` rows (`Q = 4..80`) and **`4748/4748`** rows over
      ALL `c_0 in [0,2Q)` at `Q = 5..25` -- the proof's stated generality holds.

  ⇒⇒ CONSEQUENCE: [NCYL-294]'s H2 moves VERIFIED -> PROVED, and `p = 0` ([NCYL-302]) loses
     its weakest input.  ⚠ It does NOT become unconditional: [NCYL-283]'s Gauss-Bonnet step
     is still DERIVED AND NOT SCORED.  See [NCYL-303].

  python3.13 probes/s455_closure_proof.py [all|h1|h2|h3|h4|h5] [--qmax=N]
"""
from __future__ import annotations

import json
import math
import sys
import time

import s450_chain_maps as cm

OUT = "data/s455_closure.json"


# ---------------------------------------------------------------- the reduced quantity

def S_zero(ch, osec=None):
    """`S = sum_k osec[k]*(ell[k-1]-ell[k])`, EXACT (ring `Coord`).  True iff it vanishes."""
    osec = ch.osec if osec is None else osec
    S = ch.zero
    for k in range(ch.Q):
        if osec[k]:
            S = S + ch.ell[(k - 1) % ch.Q] - ch.ell[k]
    return S.is_zero()


def closes(ch):
    """`s451_quotient_path.build`'s verdict, but on a SUPPLIED chain so the perturbed
    sector/order controls are honoured (the library `build` rebuilds its own Chain)."""
    Q = ch.Q
    u = {0: ch.zero}
    for k in range(1, Q):
        far, near, d = ch.far[k], ch.near[k], ch.delta[k]
        if far in u:
            u[near] = u[far] + d if ch.osec[k] else u[far]
        else:
            u[far] = u[near] - d if ch.osec[k] else u[near]
    far, near, d = ch.far[0], ch.near[0], ch.delta[0]
    return bool((u[far] + d if ch.osec[0] else u[far]).eq(u[near]))


def odd_rows(qmax, qmin=5):
    for Q in range(qmin | 1, qmax + 1, 2):
        for P in range(1, Q):
            if math.gcd(P, Q) == 1:
                for eps in (0, 1):
                    yield P, Q, eps


# ---------------------------------------------------------------- H1: the reduction

def h1(qmax=31):
    """closes <=> S==0.  ⚠ VACUOUS on the true model; the perturbed arms carry it."""
    out = {}
    for label, kw in [("true", {}), ("sector=shift", {"sector": "shift"}),
                      ("sector=parity", {"sector": "parity"})]:
        n = agree = bt = bf = 0
        for P, Q, eps in odd_rows(qmax):
            ch = cm.Chain(P, Q, eps, **kw)
            a, b = closes(ch), bool(S_zero(ch))
            n += 1
            agree += (a == b)
            bt += (a and b)
            bf += ((not a) and (not b))
        out[label] = {"rows": n, "agree": agree, "both_true": bt, "both_false": bf}
    tot = sum(v["rows"] for v in out.values())
    out["TOTAL"] = {"rows": tot, "agree": sum(v["agree"] for v in out.values() if "agree" in v),
                    "both_false": sum(v["both_false"] for v in out.values() if "both_false" in v)}
    return out


# ---------------------------------------------------------------- H2: the proof steps

def h2(qmax=45):
    """Each step of the PROOF, separately.  ⚠ derivation check, not evidence."""
    keys = ("farnear", "reindex", "successor", "window", "sumprod", "sset", "symm")
    ok = dict.fromkeys(keys, 0)
    n = 0
    for P, Q, eps in odd_rows(qmax):
        ch = cm.Chain(P, Q, eps)
        n += 1
        twoQ, c, c0 = 2 * Q, ch.c, ch.c[0]
        sin = lambda x: math.sin(x * math.pi / twoQ)                      # noqa: E731

        # (0) the kite constraint is SYMMETRIC in far/near -- the far/near split cancels
        ok["farnear"] += all(
            (ch.far[k], ch.near[k]) in (((k - 1) % Q, k), (k, (k - 1) % Q))
            and abs((sin(c[ch.far[k]]) - sin(c[ch.near[k]])) - ch.delta[k].f) < 1e-9
            for k in range(Q))

        w = [1 if 0 < (twoQ - c[t]) % twoQ < 2 * P else 0 for t in range(Q)]
        ok["reindex"] += all(w[(k - 1) % Q] == ch.osec[k] for k in range(Q))
        ok["successor"] += all(c[(t + 1) % Q] == (c[t] + 2 * P) % twoQ for t in range(Q))
        ok["window"] += all(w[t] == int(twoQ - 2 * P < c[t] < twoQ) for t in range(Q))

        S = sum(sin(c[(k - 1) % Q]) - sin(c[k]) for k in range(Q) if ch.osec[k])
        inner = sum(math.sin((c[t] + P) * math.pi / twoQ) for t in range(Q) if w[t])
        ok["sumprod"] += abs(S - 2 * math.cos(P * math.pi / twoQ) * inner) < 1e-9

        sset = sorted(P - (twoQ - c[t]) for t in range(Q) if w[t])
        pred = sorted(s for s in range(-(P - 1), P) if (s - (P - c0)) % 2 == 0)
        ok["sset"] += (sset == pred)
        ok["symm"] += (sset == sorted(-s for s in sset))
    return {"rows": n, **{k: ok[k] for k in keys}}


# ---------------------------------------------------------------- H3: coverage

def h3(qmax=81):
    """The proof's own PREDICTION at Q far beyond [NCYL-294]'s 31.  Tests that the proof
    describes the CODE; it cannot be evidence for the proof."""
    n = good = 0
    fails = []
    for P, Q, eps in odd_rows(qmax):
        n += 1
        if S_zero(cm.Chain(P, Q, eps)):
            good += 1
        elif len(fails) < 8:
            fails.append([P, Q, eps])
    return {"rows": n, "S_zero": good, "qmax": qmax, "fails": fails}


# ---------------------------------------------------------------- H4: THE CONTROLS

def h4(qmax=31):
    """⇒⇒ The only arm that can fail.  Every wrong-parameter model must BREAK closure."""
    res = {k: 0 for k in ("true", "sector=shift", "sector=parity", "order=P+2",
                          "window=P+1", "window=P-1")}
    applic = dict(res)
    n = 0
    for P, Q, eps in odd_rows(qmax):
        n += 1
        twoQ = 2 * Q
        for label, kw in [("true", {}), ("sector=shift", {"sector": "shift"}),
                          ("sector=parity", {"sector": "parity"})]:
            res[label] += bool(S_zero(cm.Chain(P, Q, eps, **kw)))
            applic[label] += 1
        o = P + 2
        if 0 < o < Q and math.gcd(o, Q) == 1:
            res["order=P+2"] += bool(S_zero(cm.Chain(P, Q, eps, order=o)))
            applic["order=P+2"] += 1
        ch = cm.Chain(P, Q, eps)
        for label, W in [("window=P+1", P + 1), ("window=P-1", P - 1)]:
            if not 0 < W < Q:
                continue
            osec = {k: int(0 < (twoQ - ch.c[(k - 1) % Q]) % twoQ < 2 * W) for k in range(Q)}
            res[label] += bool(S_zero(ch, osec))
            applic[label] += 1
    return {"rows": n, "S_zero": res, "applicable": applic}


# ---------------------------------------------------------------- H5: the scope claim

def h5(qmax_even=80, qmax_c0=25):
    """The proof never uses Q odd nor which `eps` made `c_0`.  Score that generality.
    ⚠ OUT OF THE MODEL'S SCOPE ([NCYL-291]) -- this is about the ALGEBRA only."""
    n_e = ok_e = 0
    for Q in range(4, qmax_even + 1, 2):
        for P in range(1, Q):
            if math.gcd(P, Q) != 1:
                continue
            for eps in (0, 1):
                try:
                    ch = cm.Chain.__new__(cm.Chain)
                    _init_raw(ch, P, Q, (P if eps == 0 else P + Q) % (2 * Q))
                except Exception:                                          # noqa: BLE001
                    continue
                n_e += 1
                ok_e += bool(S_zero(ch))
    n_c = ok_c = 0
    for Q in range(5, qmax_c0 + 1, 2):
        for P in range(1, Q):
            if math.gcd(P, Q) != 1:
                continue
            for c0 in range(2 * Q):
                ch = cm.Chain.__new__(cm.Chain)
                _init_raw(ch, P, Q, c0)
                n_c += 1
                ok_c += bool(S_zero(ch))
    return {"even_rows": n_e, "even_S_zero": ok_e,
            "c0_rows": n_c, "c0_S_zero": ok_c}


def _init_raw(ch, P, Q, c0):
    """A bare Chain carrying only what `S_zero` reads (`c`, `ell`, `osec`, `zero`, `Q`).
    ⚠ NOT a billiard model -- `Chain.__init__` asserts Q odd, and at even Q / arbitrary
    `c_0` the necklace does not exist ([NCYL-291]).  This exercises the ALGEBRA only."""
    twoQ = 2 * Q
    R = cm.get_ring(Q)
    ch.P, ch.Q, ch.ring = P, Q, R
    ch.c = [(c0 + 2 * P * t) % twoQ for t in range(Q)]
    ch.ell = [cm.Coord(2 * R.gen[ch.c[t]], math.sin(ch.c[t] * math.pi / twoQ))
              for t in range(Q)]
    ch.zero = cm.Coord(cm.np.zeros(R.d, dtype=cm.np.int64), 0.0)
    ch.osec = {k: int(0 < (twoQ - ch.c[(k - 1) % Q]) % twoQ < 2 * P) for k in range(Q)}


# ---------------------------------------------------------------- driver

def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "all"
    qmax = 81
    for a in sys.argv[2:]:
        if a.startswith("--qmax="):
            qmax = int(a.split("=")[1])
    t0 = time.time()
    out = {"mode": mode, "qmax": qmax}

    if mode in ("all", "h4"):                      # ⚠ [OPS-219]: read the controls FIRST
        out["h4"] = r = h4(min(qmax, 31))
        print("H4 ⇒⇒ CONTROLS (must BREAK; true model must not)")
        for k in r["S_zero"]:
            print(f"     S==0 on {k:14s}: {r['S_zero'][k]}/{r['applicable'][k]}")
    if mode in ("all", "h1"):
        out["h1"] = r = h1(min(qmax, 31))
        t = r["TOTAL"]
        print(f"H1 closes <=> S==0 : {t['agree']}/{t['rows']}   "
              f"BOTH-FALSE (the informative rows) = {t['both_false']}")
        for k in ("true", "sector=shift", "sector=parity"):
            v = r[k]
            print(f"     {k:14s}: {v['agree']}/{v['rows']}  "
                  f"[both True {v['both_true']}, both False {v['both_false']}]")
    if mode in ("all", "h2"):
        out["h2"] = r = h2(min(qmax, 45))
        print(f"H2 proof steps (DERIVATION CHECK, {r['rows']} rows): " +
              "  ".join(f"{k}={r[k]}" for k in r if k != "rows"))
    if mode in ("all", "h3"):
        out["h3"] = r = h3(qmax)
        print(f"H3 coverage S==0 EXACT: {r['S_zero']}/{r['rows']}  (Q<={r['qmax']})"
              f"   fails={r['fails']}")
    if mode in ("all", "h5"):
        out["h5"] = r = h5(min(qmax, 80), min(qmax, 25))
        print(f"H5 scope: even Q {r['even_S_zero']}/{r['even_rows']}   "
              f"all c_0 {r['c0_S_zero']}/{r['c0_rows']}")

    out["elapsed_s"] = round(time.time() - t0, 1)
    with open(OUT, "w") as f:
        json.dump(out, f, indent=1)
    print(f"-> {OUT}   [{out['elapsed_s']}s]")


if __name__ == "__main__":
    main()
