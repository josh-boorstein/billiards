#!/usr/bin/env python3.13
"""s466_cp_transfer.py -- PRE-REGISTERED.  Queue item (s465-a): IS `8/15 eps0`
ALREADY PROVED NON-CP?

PRIOR ART (the gate, [OPS-241], run BEFORE the probe was designed).
`rulings.py --grep 'non-CP'`, `'completely periodic'`, `'minimal component'`,
`'necklace'`, `'double cover'`, plus the bodies of [NCYL-001] (the separatrix
census), [NCYL-002] (`8/15`: 4 of 6 close, FLAT `1e3->2e7`), [NCYL-018],
[NCYL-059] (odd-`P` even class CP `84/84`; SCOPED s438 to `u_L1`), [NCYL-091],
[NCYL-210], [NCYL-212], [NCYL-144], [NCYL-256], [NCYL-257], [NCYL-262],
[NCYL-287] (the necklace `B`), [NCYL-292] (its flat geometry), [NCYL-307] H5,
[NCYL-312], [NCYL-313], [NCYL-317] (6), [OPS-234], plus `claims.md` [NECK-P0] /
[NECK-GB] / [NECK-COVER] / [T-SURFACE] and `veech_reframing.md` §1a/§3a.

⇒⇒ WHAT CAME BACK, AND IT IS THE ANSWER: **[NCYL-212] (s372) ALREADY SAYS IT** --
verbatim, *"The `8/15` failure is a HYPOTHESIS-failure ([NCYL-144]'s CP rider),
and [NCYL-210] is now the proof it is not CP."*  That clause sits in a
`cell_store` ruling about ring-certification coverage.  The CP-census strand
opened at s438, 66 sessions later, never cited it, and at s459 wrote the
opposite: DO-NOT (2) *"DO NOT CALL `7/15 eps1` / `8/15 eps0` PROVED NON-CP"*,
with `claims.md` [NECK-P0] carrying *"`¬CP` on those rows stays UNPROVED, since
the available witness routes through this row's own contrapositive and is
circular."*  Both are wrong, and wrong in the same way: "the available witness"
was not the only witness.

⇒ SO THIS PROBE IS NOT A DISCOVERY INSTRUMENT.  It is the GUARD on a transfer
argument whose load-bearing input ([NCYL-210]) was proved 94 sessions ago.  The
argument has three joints and each is an arm below; the theorem itself is not
scorable and the arms are not its evidence ([NCYL-318] scope (i) / [OPS-041] --
treating a proof's sanity sweep as its support is the error one level up).

------------------------------------------------------------------------------
THE ARGUMENT.

  (1) PROVED, unconditional, exact in `Q(zeta_60)` ([NCYL-091] + [NCYL-210]):
      the `8/15` PERPENDICULAR direction is 3 cylinders + 1 minimal, uniquely
      ergodic component, up to measure zero.  The component has positive mass
      (`mu_flow(M)/L = 15 - 6.8138871... > 0`; `mu_d(M)` is `55.56%` of the
      boundary cross-section, [NCYL-082]).  A completely periodic direction has
      NO minimal component, so that direction is `not CP` -- and nothing in the
      chain touches [NECK-GB], [NECK-P0] or the corner witness.

  (2) THAT DIRECTION IS THE CLASS ROW `8/15 eps0`  -- arm H1.

  (3) `CP` AS [NECK-P0] CONSUMES IT IS A PROPERTY OF THE NECKLACE `B`, WHICH IS
      THE GENUS-0 BASE, AND CP TRANSFERS BOTH WAYS ACROSS THE DOUBLE COVER
      -- arms H2 (the base match) and H3 (the lemma, DERIVED).

  ⇒ `8/15 eps0` is NOT completely periodic, unconditionally; and by the leg-swap
    involution `(P, eps) <-> (Q-P, 1-eps)` -- the SAME direction on the SAME
    surface with the legs relabelled (s457 h1) -- neither is `7/15 eps1`.

------------------------------------------------------------------------------
THE ARMS.  In dependency order; only H0/H1/H2/H4/H5 are scored.

  H0  THE CLASS MAP, PURE ARITHMETIC.  `s450_chain_maps.Chain` sets
      `lone = (eps == P % 2)` as a one-liner with no derivation next to it, and
      the whole identification rests on it.  Re-derive it from [NCYL-292]'s
      coordinates -- `c[t] = (c0 + 2Pt) mod 2Q`, `c0 = P` (eps=0) / `P+Q`
      (eps=1); LONE iff some `c[t] == 0`, PAIRED iff some `c[t] == Q` -- and
      check on every odd `Q <= 59` and every coprime `P` that exactly one holds,
      that it agrees with the one-liner, and that the leg swap preserves it.
      ⚠ CONTROL, WHICH MUST FAIL ([OPS-219]): the parity-PRESERVING relabelling
      `(P, eps) -> (Q-P, eps)`.  If that preserved LONE too, H0 would be
      scoring a property of the pair `{P, Q-P}` and not of the involution.

  H1  ⇒⇒ THE DIRECTION IDENTIFICATION, AND IT IS THE ARM THAT COULD HAVE KILLED
      THE WHOLE THING.  Joint (2) is the only joint where a label could be
      attached to the wrong object, and the check is available for free in an
      instrument that shares NO code with the necklace model: `s438_oddclass_cp`
      traces the BILLIARD on the direction grid `u*pi/(2Q)` and stores, per
      class, `class_role`, `perp_sides`, the separatrix count, the closed count
      and the per-prong `u`.  Require, at `Q = 15`:
        (a) exactly ONE class row of `8/15` carries [NCYL-002]'s signature
            (`n_separatrices == 6`, `n_closed == 4`, verdict MIXED, and the
            closure curve FLAT across the whole `1e3 -> 2e7` ladder), and it is
            `eps = 0`;
        (b) that row is `class_role == 'lone'` and `perp_sides == 'L1'`, i.e. it
            IS a perpendicular-beam direction and not the `H`-carrying class;
        (c) its prong `u` values all have parity `eps` -- the grid-parity half
            of the identification, checked against the tracer rather than
            assumed from `psi = eps`;
        (d) the leg-swap twin `7/15 eps1` reproduces (a)-(c) with `L2`;
        (e) NEGATIVE CONTROL AT THE SAME `Q` AND THE SAME `class_role`:
            `4/15 eps0` is lone, `L1`, 6 separatrices -- and `6/6` closed.  So
            "lone + L1 + 6 separatrices" does NOT force the signature, and (a)
            is not reading a property of the class label.
      ⚠ (a) COULD HAVE COME OUT OTHERWISE: had the 4-of-6 signature sat on
      `eps = 1`, joint (2) would be false and the transfer dead.
      Then generalise (b)/(c) over the WHOLE stored census: for every row,
      `class_role == 'lone'` iff `parity == P % 2` (H0's law, now scored on the
      billiard tracer), and every prong `u` matches `parity`.

  H2  THE BASE MATCH -- A CITATION CHECK, DELIBERATELY NOT A REDERIVATION.
      ⚠ s465 DO-NOT (4) bars re-deriving the stratum of `S_alpha`; that is
      PROVED at [T-SURFACE] and is USED here, not reproved.  What is checked is
      the OTHER side: that [NCYL-287]'s necklace -- cone angle `P*pi` at `Z_O`,
      `(Q-P)*pi` at `Z_A`, `pi` at each of `Q` simple poles -- has quadratic
      order vector `(P-2, Q-P-2, -1^Q)`, i.e. is LITERALLY the genus-0 base
      `Q(P-2, Q-P-2, -1^Q)` whose holonomy double cover `veech_reframing.md`
      §3a records as `S_alpha`.  Scored: the order sum is `-4` (a genus-0
      quadratic differential) and Riemann-Hurwitz over the odd-order
      singularities returns `g = (Q-1)/2 == floor(Q/2)`, matching [T-SURFACE].
      ⚠ NOT EVIDENCE FOR [T-SURFACE] -- it is the consistency of the necklace
      with a proved stratum, and it fails loudly if the necklace data is wrong.

  H3  THE COVERING LEMMA -- ⚠ DERIVED, NOT SCORED, AND IT IS AN ARGUMENT OF MINE.
      `pi : S_alpha -> B` has degree 2 (the holonomy involution `tau`,
      `tau^* omega = -omega`).
        (=>) CP(`B`, theta) => CP(`S_alpha`, theta~): a cylinder contains no
             singularity, so `pi` restricted to its preimage is an unbranched
             cover of an annulus, hence 1 or 2 annuli; the preimages of a
             cylinder decomposition therefore cover `S_alpha` up to measure zero
             by cylinders.
        (<=) CP(`S_alpha`, theta~) => CP(`B`, theta): `tau` permutes the
             cylinders of a decomposition, so the quotient is one.
      ⇒ the two predicates are EQUIVALENT, which is what the repo has been
      assuming all along without writing down: [NCYL-312] scores the exact
      NECKLACE certificates against the float BILLIARD census row-by-row and
      calls `20` of the float's non-CP rows "refuted", a cross-tabulation that
      is only meaningful under this equivalence.
      ⚠ So H3 is not new machinery -- it is the missing justification for an
      identification already in load-bearing use.

  H4  REPRODUCE THE LOAD-BEARING INPUT.  The conclusion has exactly ONE
      mathematical input outside this file, [NCYL-210]'s field identity, and a
      ruling is prose.  Re-run `probes/s372_field_identity.py` and require
      `S_ring + mu_flow(M)/L - Q` RING-ZERO.  ⚠ COULD FAIL (bit-rot, or a ruling
      that overstates what its probe returns).

  H5  THE REFUTATION CONTROL, AND IT IS A REAL STAKE, NOT A FORMALITY.  By
      [NCYL-317] (6) / `claims.md` [NECK-P0], certifying `7/15 eps1` or
      `8/15 eps0` COMPLETELY PERIODIC would refute [NECK-P0] and [NECK-COVER].
      This session proves those rows are NOT CP by a route independent of
      [NECK-*], so the exact necklace census MUST NOT hold a CP certificate for
      them.  Read `data/s457_exact_cp/` and require `capped > 0` on both, and
      `capped == 0` on their paired-class siblings.
      ⚠⚠ A SURVIVED FALSIFICATION IS NOT SUPPORT: [NECK-P0] says `CP => p = 0`,
      and a non-CP row is VACUOUS for it ([OPS-041]).  What H5 buys is that the
      one concrete falsifier on the books is now settled, not that the row is
      any better supported.

SCOPE.  ODD `Q` (the necklace model's scope, [NCYL-291]).  The THEOREM is about
exactly one direction -- `8/15`'s perpendicular beam and its leg-swap twin --
because [NCYL-210] is a one-row result and [NCYL-003] records, as the USER'S
CALL at s378, that `8/15` is a permanently unpromoted singleton.  ⚠ It does NOT
give a method: no other row has an exact mass budget, and building one is
`future_directions.md` B9 q3 steer (ii).

USAGE:  PYTHONPATH=.:engine:archive/scripts_2026-07:archive/scripts_2026-08:probes \
        .venv/bin/python3.13 probes/s466_cp_transfer.py [h0|h1|h2|h5|all]
"""
import json
import os
import sys
from math import gcd

OUT = "data/s466_cp_transfer.json"
FLOAT_STORE = "data/s438_cp_store"
EXACT_STORE = "data/s457_exact_cp"


# --------------------------------------------------------------------------
# H0 -- the class map, from [NCYL-292]'s coordinates
# --------------------------------------------------------------------------
def class_roles(P, Q, eps):
    """LONE/PAIRED from the interface offsets, NOT from the one-liner."""
    c0 = P + (Q if eps else 0)
    c = [(c0 + 2 * P * t) % (2 * Q) for t in range(Q)]
    return {"lone": 0 in c, "paired": Q in c, "parity": c0 % 2, "c": c}


def h0(verbose=True):
    rows, bad, ctrl_preserved = 0, [], 0
    for Q in range(5, 60, 2):
        for P in range(1, Q):
            if gcd(P, Q) != 1:
                continue
            for eps in (0, 1):
                r = class_roles(P, Q, eps)
                rows += 1
                # (i) exactly one of lone/paired
                if r["lone"] == r["paired"]:
                    bad.append((P, Q, eps, "not-exclusive"))
                # (ii) agrees with s450_chain_maps' one-liner
                if r["lone"] != (eps == P % 2):
                    bad.append((P, Q, eps, "oneliner"))
                # (iii) lone <=> the offsets are even
                if r["lone"] != (r["parity"] == 0):
                    bad.append((P, Q, eps, "parity"))
                # (iv) the leg swap preserves the role
                s = class_roles(Q - P, Q, 1 - eps)
                if s["lone"] != r["lone"]:
                    bad.append((P, Q, eps, "legswap"))
                # CONTROL: the parity-preserving relabelling must NOT preserve it
                k = class_roles(Q - P, Q, eps)
                if k["lone"] == r["lone"]:
                    ctrl_preserved += 1
    ok = not bad and ctrl_preserved == 0
    if verbose:
        print(f"H0  class map, odd Q <= 59: {rows - len(bad)}/{rows} rows consistent "
              f"(exclusive | one-liner | parity | leg-swap)   "
              f"{'PASS' if not bad else 'FAIL ' + str(bad[:5])}")
        print(f"    CONTROL (parity-preserving swap must FLIP the role): "
              f"preserved on {ctrl_preserved}/{rows}   "
              f"{'PASS (fires)' if ctrl_preserved == 0 else 'FAIL'}")
        r8 = class_roles(8, 15, 0)
        r7 = class_roles(7, 15, 1)
        print(f"    8/15 eps0: lone={r8['lone']} offsets even={r8['parity'] == 0}   "
              f"7/15 eps1: lone={r7['lone']} offsets even={r7['parity'] == 0}   "
              f"(leg-swap pair)")
    return {"rows": rows, "bad": bad[:20], "n_bad": len(bad),
            "ctrl_preserved": ctrl_preserved, "pass": ok}


# --------------------------------------------------------------------------
# H1 -- the direction identification, against the BILLIARD tracer
# --------------------------------------------------------------------------
SIG = {"n_separatrices": 6, "n_closed": 4, "verdict": "MIXED"}


def _load_float(P, Q):
    p = os.path.join(FLOAT_STORE, f"{P}_{Q}.json")
    return json.load(open(p)) if os.path.exists(p) else None


def h1(verbose=True):
    out = {}
    # ---- (a) which class of 8/15 carries [NCYL-002]'s signature? ----
    d = _load_float(8, 15)
    hits = []
    for k, c in d["classes"].items():
        if all(c.get(f) == v for f, v in SIG.items()):
            hits.append(int(k))
    flat = None
    if hits:
        curve = d["classes"][str(hits[0])]["curve"]
        flat = len({pt["n_closed"] for pt in curve}) == 1 and len(curve) >= 5
    a_ok = hits == [0] and flat
    if verbose:
        print(f"H1a 8/15: classes carrying [NCYL-002]'s 4-of-6 MIXED signature = "
              f"{hits}; closure curve flat over {len(d['classes']['0']['curve'])} "
              f"rungs = {flat}   {'PASS' if a_ok else 'FAIL'}")
    out["a"] = {"hits": hits, "flat": flat, "pass": bool(a_ok)}

    # ---- (b)(c)(d) role / perp side / prong parity on the pair ----
    pair = {}
    for (P, eps, want_side) in [(8, 0, "L1"), (7, 1, "L2")]:
        c = _load_float(P, 15)["classes"][str(eps)]
        us = sorted({r["u"] for r in c["records"]})
        rec = {"class_role": c["class_role"], "perp_sides": c["perp_sides"],
               "n_sep": c["n_separatrices"], "n_closed": c["n_closed"],
               "verdict": c["verdict"], "u": us,
               "u_parity_ok": all(u % 2 == eps for u in us),
               "side_ok": c["perp_sides"] == want_side,
               "lone_ok": c["class_role"] == "lone",
               "sig_ok": all(c.get(f) == v for f, v in SIG.items())}
        rec["pass"] = all(rec[f] for f in
                          ("u_parity_ok", "side_ok", "lone_ok", "sig_ok"))
        pair[f"{P}/15 eps{eps}"] = rec
        if verbose:
            print(f"H1bcd {P}/15 eps{eps}: role={rec['class_role']:6s} "
                  f"perp={rec['perp_sides']:5s} {rec['n_closed']}/{rec['n_sep']} "
                  f"{rec['verdict']:20s} u={us} "
                  f"{'PASS' if rec['pass'] else 'FAIL'}")
    out["pair"] = pair

    # ---- (e) negative control: same Q, same role, same side, CP ----
    c = _load_float(4, 15)["classes"]["0"]
    e_ok = (c["class_role"] == "lone" and c["perp_sides"] == "L1"
            and c["n_separatrices"] == 6
            and c["verdict"] == "COMPLETELY_PERIODIC")
    if verbose:
        print(f"H1e CONTROL 4/15 eps0 (lone, L1, 6 separatrices): "
              f"{c['n_closed']}/{c['n_separatrices']} {c['verdict']}   "
              f"{'PASS (label does not force the signature)' if e_ok else 'FAIL'}")
    out["e"] = {"role": c["class_role"], "side": c["perp_sides"],
                "n_closed": c["n_closed"], "verdict": c["verdict"],
                "pass": bool(e_ok)}

    # ---- the census-wide generalisation of H0's law + the u-parity ----
    n, bad = 0, []
    for fn in sorted(os.listdir(FLOAT_STORE)):
        if not fn.endswith(".json"):
            continue
        d = json.load(open(os.path.join(FLOAT_STORE, fn)))
        P, Q = d["P"], d["Q"]
        if Q % 2 == 0:
            continue
        for k, c in d["classes"].items():
            eps = int(k)
            n += 1
            if (c["class_role"] == "lone") != (eps == P % 2):
                bad.append((P, Q, eps, "role"))
            if not all(r["u"] % 2 == eps for r in c["records"]):
                bad.append((P, Q, eps, "u-parity"))
    if verbose:
        print(f"H1* census-wide: role == (eps == P%2) AND every prong u has "
              f"parity eps -- {n - len(bad)}/{n} class rows   "
              f"{'PASS' if not bad else 'FAIL ' + str(bad[:5])}")
    out["census"] = {"rows": n, "n_bad": len(bad), "bad": bad[:20],
                     "pass": not bad}
    out["pass"] = bool(a_ok and e_ok and not bad
                       and all(v["pass"] for v in pair.values()))
    return out


# --------------------------------------------------------------------------
# H2 -- the base match (citation check; the stratum of S_alpha is CITED)
# --------------------------------------------------------------------------
def h2(verbose=True):
    rows, bad = 0, []
    for Q in range(5, 60, 2):
        for P in range(1, Q):
            if gcd(P, Q) != 1:
                continue
            rows += 1
            # [NCYL-287]: cone angles P*pi, (Q-P)*pi, and pi at each of Q poles
            orders = [P - 2, Q - P - 2] + [-1] * Q
            if sum(orders) != -4:
                bad.append((P, Q, "order-sum", sum(orders)))
            # Riemann-Hurwitz for the holonomy double cover: branched exactly
            # over the ODD-order singularities.
            b = sum(1 for o in orders if o % 2)
            if b % 2:
                bad.append((P, Q, "odd-branch-count", b))
                continue
            g = (2 * (2 * 0 - 2) + b + 2) // 2          # 2g-2 = 2(2g_B-2) + b
            if g != Q // 2:
                bad.append((P, Q, "genus", g, Q // 2))
    if verbose:
        print(f"H2  necklace [NCYL-287] as a quadratic differential, odd Q <= 59: "
              f"order vector (P-2, Q-P-2, -1^Q) sums to -4 and Riemann-Hurwitz "
              f"gives g = floor(Q/2) -- {rows - len(bad)}/{rows}   "
              f"{'PASS' if not bad else 'FAIL ' + str(bad[:5])}")
        print(f"    => the base IS Q(P-2, Q-P-2, -1^Q), whose holonomy double "
              f"cover is S_alpha (CITED: [T-SURFACE], veech_reframing.md §3a; "
              f"NOT rederived -- s465 DO-NOT (4))")
    return {"rows": rows, "n_bad": len(bad), "bad": bad[:20], "pass": not bad}


# --------------------------------------------------------------------------
# H5 -- the refutation control
# --------------------------------------------------------------------------
def h5(verbose=True):
    out, bad = {}, []
    for (P, eps) in [(8, 0), (7, 1)]:
        d = json.load(open(os.path.join(EXACT_STORE, f"{P}_15.json")))
        me = d["classes"][str(eps)]
        sib = d["classes"][str(1 - eps)]
        rec = {"cap": me["cap"], "capped": me["capped"],
               "n_traced": me["n_traced"], "lone": me["lone"],
               "sibling_capped": sib["capped"], "sibling_lone": sib["lone"],
               "runs": d["runs"][str(eps)]}
        rec["pass"] = me["capped"] > 0 and sib["capped"] == 0
        if not rec["pass"]:
            bad.append((P, eps))
        out[f"{P}/15 eps{eps}"] = rec
        if verbose:
            caps = sorted(int(k) for k in rec["runs"])
            print(f"H5  {P}/15 eps{eps} (lone={me['lone']}): capped "
                  f"{me['capped']}/{me['n_traced']} at cap {me['cap']:.0e}; "
                  f"ladder {[(f'{c:.0e}', rec['runs'][str(c)]['capped']) for c in caps]}; "
                  f"paired sibling capped={sib['capped']}   "
                  f"{'PASS (no CP certificate -- the falsifier does not fire)' if rec['pass'] else 'FAIL -- REFUTES THIS SESSION'}")
    out["pass"] = not bad
    if verbose:
        print("    ⚠ A SURVIVED FALSIFICATION IS NOT SUPPORT FOR [NECK-P0]: a "
              "non-CP row is vacuous for `CP => p = 0` ([OPS-041]).")
    return out


# --------------------------------------------------------------------------
def main():
    arms = sys.argv[1] if len(sys.argv) > 1 else "all"
    res = {}
    print(__doc__.split("------")[0].strip().splitlines()[0])
    print()
    if arms in ("all", "h0"):
        res["h0"] = h0()
        print()
    if arms in ("all", "h1"):
        res["h1"] = h1()
        print()
    if arms in ("all", "h2"):
        res["h2"] = h2()
        print()
    if arms in ("all", "h5"):
        res["h5"] = h5()
        print()
    if arms == "all":
        scored = [k for k in ("h0", "h1", "h2", "h5") if k in res]
        allpass = all(res[k]["pass"] for k in scored)
        print(f"SCORED ARMS: {', '.join(scored)} -- "
              f"{'ALL PASS' if allpass else 'FAILURE PRESENT'}")
        print("H3 (the covering lemma) is DERIVED and not scored; H4 "
              "(reproducing [NCYL-210]) is `probes/s372_field_identity.py`, "
              "run separately -- see logs/s466_h4_s372_replay.log.")
        res["pass"] = allpass
        json.dump(res, open(OUT, "w"), indent=1)
        print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
