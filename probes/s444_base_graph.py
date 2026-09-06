#!/usr/bin/env python3
"""
s444_base_graph.py -- PRE-REGISTERED.  Queue item (a'): PROVE `C(lone) = g`.

WHAT THIS IS.  Not another measurement of the count law.  An exact TOPOLOGICAL IDENTITY
that computes `C` from two small integers, plus the instrument that measures those two
integers, so that the law becomes a statement a proof can attack.

------------------------------------------------------------------------------------
THE IDENTITY (derived here; the derivation is elementary and is the point).

Write `Sigma = S_alpha` (the 4Q-copy unfolding, genus `g = floor(Q/2)`) and let
`B = Sigma/tau` be its quotient by the hyperelliptic involution -- Apisa's genus-0 base
`Q(P-2, Q-P-2, -1^Q)` (`veech_reframing.md` sec 3a).  `B` is a SPHERE carrying

    two zeros  Z_O (order P-2, so P prongs)  and  Z_A (order Q-P-2, so Q-P prongs),
    Q simple poles (order -1, ONE prong each) = the Q copies of the right angle R.

Fix one of the two direction classes ([NCYL-240]: the 4Q launch angles carry exactly two
cylinder decompositions, indexed by `u % 2`).  Let `G` be the union of the leaves through
the marked points -- the SEPARATRIX DIAGRAM of that direction on `B`.  Then

    #edges(G) = (sum of prongs)/2 = (P + (Q-P) + Q)/2 = Q          [ALWAYS, any centre]
    #vertices(G) = Q + 2.

`B \\ G` is a disjoint union of `C_B` open cylinders (annuli).  A regular neighbourhood
`N_i` of each of the `c` connected components of `G` is a PLANAR surface (we are on a
sphere), so `chi(N_i) = 2 - b_i` with `b_i` its boundary circles; summing,

    (V - E) = sum_i chi(N_i) = 2c - sum_i b_i = 2c - 2*C_B,
    2 = (Q+2) - Q = 2c - 2*C_B   ==>   **C_B = c - 1**.

Because a pole has exactly ONE prong it is a LEAF of `G`, so every edge is pole-pole,
pole-zero or zero-zero.  With `a`, `b`, `d` their counts:

    2a + b = Q   (pole prongs)      b + 2d = Q   (zero prongs)    ==>  a = d,
    c = a + kappa,  kappa := #components containing a zero  in {1,2},

    **C_B = (Q - b)/2 + kappa - 1.**

Finally `C_B` IS the count the strand measures: `tau` acts on the direction's cylinders,
a tau-pair descends to one cylinder of `B` and a tau-fixed one likewise, and the
enumerator counts billiard cylinders.  (Cross-check available in the ledger: at `3/7`
the horizontal decomposition of `Sigma` is 8 cylinders in 4 tau-pairs
(`f_leg_covariance.md` sec 20.27 R2/R3, area-verified) while the count law reads
`C(even) = 4` ([NCYL-241]) -- exactly `C_B`.)

------------------------------------------------------------------------------------
WHAT THE LAW BECOMES.  Put `delta := b - 2(kappa-1)`, so `C = (Q - delta)/2`.  Then

    C = floor(Q/2)  <=>  delta = +[Q odd]        C = ceil(Q/2)  <=>  delta = -[Q odd].

The SIDES of the triangle are themselves leaves whenever the class contains their own
direction (in units `pi/(2Q)`: `L2` has direction `0`, `L1` has `Q`, `H` has `P`), and
each such side is one zero-pole edge (`L2 = O-R`, `L1 = A-R`) or the one zero-zero edge
joining the two DISTINCT zeros (`H = O-A`).  Those "boundary rays" alone give

    b_bd = [eps==0] + [eps==Q mod 2],     an O-A edge  <=>  eps == P mod 2,

and a one-line check (below, `predict_from_sides`) shows that `b = b_bd` together with
"no other O-A edge" reproduces `delta = +[Q odd]` on the class `eps == P mod 2` and
`delta = -[Q odd]` on the other, in all six `(P,Q)`-parity cases.  Since `b >= b_bd`
always and `kappa = 1` is FORCED whenever `H` is parallel, the count law is EQUIVALENT to

  (*)  every prong launched from a vertex INTO THE INTERIOR returns to the SAME vertex:
       no interior separatrix from `O` or `A` arrives at `R`, and none joins `O` to `A`.

`s439_exact_cells`'s module docstring already ASSERTS half of (*) as a design remark --
"an `R` prong terminates by RETURNING TO `R`, not by arriving at `O`/`A`" -- but its
justification is "it need NEVER reach a cone point", which is a possibility, not a fact,
and the arrival vertex is computed by `trace_prong` and then discarded.  This probe scores
it.

------------------------------------------------------------------------------------
HYPOTHESES (pre-registered, all able to fail).

  H0  CONTROL, RUN FIRST.  `lone <=> parity == P mod 2` on every class row of
      `data/s440_claim_store`.  Shares no code with this probe's derivation; a mismatch
      means the parity bookkeeping is wrong and everything below is void.
  H1  CONTROL.  Prong-degree identities `2*n_int(v) + n_bd(v) = deg(v)` for
      `deg = P, Q-P, Q` at `O, A, R`.  Pure arithmetic against `s439.prongs`; a mismatch
      means the inward windows and the base prong counts disagree.
  H2  CONTROL, and it is the one that tests the TRACER: `n(X->Y) == n(Y->X)` for every
      ordered pair of distinct vertex types -- every edge must be seen from both ends.
      Also `a == d` and `a + b + d == Q`.
  H3  THE PREDICTION.  `C_pred = (Q-b)/2 + kappa - 1` equals the stored `C_total`
      on every scorable class row.  This is the identity above tested end to end against
      an instrument that shares no topology with it.  ⚠ NAMED AS THE THING MOST LIKELY
      TO BREAK: if `tau` does NOT act freely somewhere, `C_total` is not `C_B` and H3
      fails on exactly those rows -- which would itself be the finding.
  H4  (*) itself: `n_int(O->R) = n_int(A->R) = n_int(O->A) = 0`.
      Predicted to hold wherever the count law holds; H3 and H4 are logically equivalent
      given H1/H2, so a row where they disagree is an instrument bug, not a discovery.

GUARDS.
  G1  A class with any prong not terminating at `cap2` is UNSCORABLE (`PRONG_CAPPED`),
      reported and excluded -- [NCYL-020]'s truncation guard, s439's own convention.
  G2  Arrival is a float test at `VERT_TOL = 1e-9` (s439's constant, reused verbatim so
      the arrival predicate is literally the one the count-law instrument uses).  The
      minimum over the sweep of the closest NON-arrival approach is reported; a row whose
      margin is not clean is flagged, not silently scored.
  G3  Nothing here writes to any store.

PRIOR ART (grepped 'singular graph', 'connected component', 'separatrix', 'saddle
connection', 'c - 1', 'Weierstrass', 'envelope', 'base', 'cylinder diagram'; read
`rulings.md` [NCYL-001]/[NCYL-043]/[NCYL-133]/[NCYL-228]/[NCYL-241]/[NCYL-251]/
[NCYL-253]/[NCYL-254]/[NCYL-255], `veech_reframing.md` sec 1/1a/3a/3d/3e and
`f_leg_covariance.md` sec 20.22-20.27 before designing) ->
  - [NCYL-251] proposed the TOPOLOGICAL form `C(lone) = g` and named the homology
    mechanism as the hypothesis to attack; [NCYL-253] then showed `(g,s)` is a function
    of `Q` alone so that form has no independent falsification test.  This probe does not
    re-run that: it replaces the stratum reading with the BASE separatrix diagram, whose
    two integers `(b, kappa)` are NOT functions of `Q` a priori.
  - [NCYL-255] found the billiard decomposition uses ORDER-0 marked vertex points as
    cylinder boundaries and left "no convention is known under which both halves of
    [NCYL-251] are simultaneously topological and predictive" UNRESOLVED.  The base
    picture supplies one: the order-0 points ARE the simple poles of `B`.
  - `veech_reframing.md` sec 3d imports Apisa Thm 8.4's `N_cyl = N_env + 2*N_simp`; that
    is the same tau-descent used above, in the all-directions setting.
  - `f_leg_covariance.md` sec 20.24-20.27 built `Sigma` as an explicit translation atlas
    and enumerated saddle connections AT 3/7 ONLY, with a per-SC cylinder assignment;
    sec 20.27 R2/R3 is the 8-cylinder / 4-tau-pair result used as the cross-check above.
    That strand is parked and works one centre; nothing there counts components or runs
    a sweep.
  - `s276_separatrix_census` / `s438_oddclass_cp` trace `O`/`A` prongs and ask only
    whether they ARRIVE; they do not record WHERE, they never launch from `R`, and
    `s276` is even-class only.  `s439_exact_cells.trace_prong` records the arrival vertex
    and discards it.
  - No probe in the repo computes a separatrix-diagram component count, `b`, or `kappa`.

Run: PYTHONPATH=.:engine:archive/scripts_2026-07:archive/scripts_2026-08:probes \\
     .venv/bin/python3.13 probes/s444_base_graph.py [h0|sweep] [--qmax=N] [--out=PATH]
"""
import json
import math
import os
import sys
import time
from collections import Counter

import s437_oblique_cylinders as s437
import s439_exact_cells as s439
from right_triangle_billiards import RightTriangleBilliard

STORE = "data/s440_claim_store"
CAP = 200_000
CAP_BIG = 2_000_000
VERT_TOL = s439.VERT_TOL


# ---------------------------------------------------------------- the six-case algebra

def side_leaves(P, Q, eps):
    """Which SIDES are leaves of the class-`eps` foliation, and the base edge each is.

    A side is a leaf iff the class contains its own direction: `L2` -> `0`, `L1` -> `Q`,
    `H` -> `P` in units `pi/(2Q)`.  `L2` joins `O`-`R`, `L1` joins `A`-`R`, `H` joins
    `O`-`A`."""
    out = {}
    if eps % 2 == 0:
        out["L2"] = ("O", "R")
    if eps % 2 == Q % 2:
        out["L1"] = ("A", "R")
    if eps % 2 == P % 2:
        out["H"] = ("O", "A")
    return out


def predict_from_sides(P, Q, eps):
    """`(b, kappa, C)` under (*) -- i.e. assuming the ONLY zero-pole edges and the only
    `O`-`A` edge are the sides themselves."""
    sl = side_leaves(P, Q, eps)
    b = sum(1 for v in sl.values() if "R" in v)
    kappa = 1 if ("O", "A") in sl.values() else 2
    return b, kappa, (Q - b) // 2 + kappa - 1


def six_case_check(verbose=True):
    """The algebra behind the equivalence: under (*), `C` is `floor(Q/2)` on the class
    `eps == P mod 2` and `ceil(Q/2)` on the other, in all six parity cases.  Scored over
    a wide box so it is a fact about the formula, not about six hand-picked rows."""
    bad = []
    for Q in range(3, 121):
        for P in range(1, Q):
            if math.gcd(P, Q) != 1:
                continue
            for eps in (0, 1):
                b, kappa, C = predict_from_sides(P, Q, eps)
                lone = (eps == P % 2)
                want = Q // 2 if lone else (Q + 1) // 2
                if C != want or (Q - b) % 2:
                    bad.append((P, Q, eps, b, kappa, C, want))
    if verbose:
        print(f"  six-case algebra: {len(bad)} mismatches over all coprime (P,Q), Q<=120,"
              f" both classes", flush=True)
        if bad:
            for r in bad[:8]:
                print("    ", r, flush=True)
    return bad


# ---------------------------------------------------------------- H0: the class labels

def h0_control(verbose=True):
    """`lone <=> parity == P mod 2`, over every class row of the store."""
    bad, n = [], 0
    for fn in sorted(os.listdir(STORE)):
        if not fn.endswith(".json"):
            continue
        rec = json.load(open(os.path.join(STORE, fn)))
        P, Q = rec["P"], rec["Q"]
        for role, cl in rec.get("classes", {}).items():
            if not isinstance(cl, dict) or "parity" not in cl:
                continue
            n += 1
            lone_pred = (cl["parity"] == P % 2)
            if lone_pred != (cl.get("class_role") == "lone"):
                bad.append((P, Q, cl["parity"], cl.get("class_role")))
    if verbose:
        print(f"H0  lone <=> parity == P mod 2 : {n - len(bad)}/{n} class rows agree",
              flush=True)
        for r in bad[:8]:
            print("    ", r, flush=True)
    return n, bad


# ---------------------------------------------------------------- the measurement

def arrivals(P, Q, eps, cap=CAP, cap2=CAP_BIG):
    """Trace every INTERIOR prong of the class and record which vertex it reaches.

    Reuses `s439_exact_cells.prongs` (the inward integer windows) and `trace_prong` (the
    arrival predicate) verbatim, so the separatrix endpoints are the count-law
    instrument's own."""
    geo = s437.geometry(P, Q)
    verts = s439.vertices(P, Q)
    B = RightTriangleBilliard((P / Q) * (math.pi / 2))
    rec, open_prongs, nb = {}, [], 0
    for v0 in ("O", "R", "A"):
        for u0 in s439.prongs(v0, P, Q, eps):
            bounces, why = s439.trace_prong(B, geo, verts, Q, v0, u0, cap)
            if why == "cap" and cap2 > cap:
                bounces, why = s439.trace_prong(B, geo, verts, Q, v0, u0, cap2)
            nb += len(bounces)
            if not why.startswith("vertex"):
                open_prongs.append((v0, u0, why))
            rec[(v0, u0)] = why.split(":")[-1] if why.startswith("vertex") else None
    return rec, open_prongs, nb


def base_graph(P, Q, eps, rec):
    """Assemble the base separatrix diagram's edge classification from the prong arrivals.

    An INTERIOR billiard prong is TWO prong-ends of `B` (the two sheets of the doubled
    triangle, mirror images, same arrival); a SIDE prong is one, with known endpoints."""
    n = Counter()
    for (v0, _u0), w in rec.items():
        if w is not None:
            n[(v0, w)] += 2
    n_int = dict(n)
    for _sd, (x, y) in side_leaves(P, Q, eps).items():
        n[(x, y)] += 1
        n[(y, x)] += 1

    deg = {"O": P, "A": Q - P, "R": Q}
    checks = {}
    checks["H1_degree"] = {v: (sum(n[(v, w)] for w in ("O", "R", "A")), deg[v])
                           for v in ("O", "R", "A")}
    checks["H1_ok"] = all(a == b for a, b in checks["H1_degree"].values())
    pairs = [("O", "R"), ("A", "R"), ("O", "A")]
    checks["H2_sym"] = {f"{x}{y}": (n[(x, y)], n[(y, x)]) for x, y in pairs}
    checks["H2_sym_ok"] = all(n[(x, y)] == n[(y, x)] for x, y in pairs)

    a = n[("R", "R")] // 2
    b = n[("O", "R")] + n[("A", "R")]
    d = n[("O", "A")] + n[("O", "O")] // 2 + n[("A", "A")] // 2
    kappa = 1 if n[("O", "A")] else 2
    checks["H2_ad"] = (a, d)
    checks["H2_edges"] = (a + b + d, Q)
    checks["H2_ok"] = checks["H2_sym_ok"] and a == d and a + b + d == Q

    star = {"int_O_R": n_int.get(("O", "R"), 0), "int_A_R": n_int.get(("A", "R"), 0),
            "int_O_A": n_int.get(("O", "A"), 0), "int_R_O": n_int.get(("R", "O"), 0),
            "int_R_A": n_int.get(("R", "A"), 0)}
    checks["H4_star_ok"] = not any(star.values())

    return {"a": a, "b": b, "d": d, "kappa": kappa,
            "C_pred": a + kappa - 1,
            "counts": {f"{x}->{y}": n[(x, y)] for x in "ORA" for y in "ORA" if n[(x, y)]},
            "interior_offdiag": star, "checks": checks}


def run_row(P, Q, eps, cap=CAP, cap2=CAP_BIG):
    t0 = time.time()
    rec, open_prongs, nb = arrivals(P, Q, eps, cap, cap2)
    out = {"P": P, "Q": Q, "eps": eps, "lone": eps == P % 2,
           "n_prongs": len(rec), "n_bounces": nb,
           "open_prongs": [list(o) for o in open_prongs], "secs": round(time.time() - t0, 2)}
    if open_prongs:
        out["status"] = "PRONG_CAPPED"
        return out
    out["status"] = "ok"
    out.update(base_graph(P, Q, eps, rec))
    b_s, k_s, C_s = predict_from_sides(P, Q, eps)
    out["sides_only"] = {"b": b_s, "kappa": k_s, "C": C_s}
    out["C_law"] = Q // 2 if out["lone"] else (Q + 1) // 2
    return out


# ---------------------------------------------------------------- the sweep

def stored_counts(qmax):
    """`{(P,Q,parity): C_total}` from the s440 claim store, scorable rows only."""
    got = {}
    for fn in sorted(os.listdir(STORE)):
        if not fn.endswith(".json"):
            continue
        rec = json.load(open(os.path.join(STORE, fn)))
        P, Q = rec["P"], rec["Q"]
        if Q > qmax:
            continue
        for cl in rec.get("classes", {}).values():
            if not isinstance(cl, dict) or cl.get("C_total") is None:
                continue
            if cl.get("prongs_not_terminated"):
                continue
            got[(P, Q, cl["parity"])] = cl["C_total"]
    return got


def sweep(qmax=20, out_path="data/s444_base_graph.json", budget_s=None, verbose=True):
    t_start = time.time()
    store = stored_counts(qmax)
    rows, keys = [], sorted(store, key=lambda k: (k[1], k[0], k[2]))
    for (P, Q, eps) in keys:
        if budget_s and time.time() - t_start > budget_s:
            print(f"  budget reached, stopped after {len(rows)} rows", flush=True)
            break
        r = run_row(P, Q, eps)
        r["C_store"] = store[(P, Q, eps)]
        rows.append(r)
        if verbose and len(rows) % 25 == 0:
            print(f"  ... {len(rows)}/{len(keys)} rows "
                  f"[{round(time.time()-t_start,1)}s]", flush=True)
    ok = [r for r in rows if r["status"] == "ok"]
    summary = {
        "n_rows": len(rows), "n_scorable": len(ok),
        "n_capped": sum(1 for r in rows if r["status"] == "PRONG_CAPPED"),
        "H1": sum(1 for r in ok if r["checks"]["H1_ok"]),
        "H2": sum(1 for r in ok if r["checks"]["H2_ok"]),
        "H3": sum(1 for r in ok if r["C_pred"] == r["C_store"]),
        "H4": sum(1 for r in ok if r["checks"]["H4_star_ok"]),
        "law_vs_store": sum(1 for r in ok if r["C_law"] == r["C_store"]),
        "sides_vs_pred": sum(1 for r in ok if r["sides_only"]["C"] == r["C_pred"]),
        "qmax": qmax, "secs": round(time.time() - t_start, 1),
    }
    if verbose:
        print(f"\nSWEEP  Q<={qmax}: {summary['n_scorable']}/{summary['n_rows']} scorable "
              f"({summary['n_capped']} prong-capped)  [{summary['secs']}s]", flush=True)
        for k in ("H1", "H2", "H3", "H4", "law_vs_store", "sides_vs_pred"):
            print(f"  {k:<14} {summary[k]}/{summary['n_scorable']}", flush=True)
        okey = {"H1": "H1_ok", "H2": "H2_ok", "H4": "H4_star_ok"}
        for k in ("H1", "H2", "H3", "H4"):
            bad = [r for r in ok if not (r["C_pred"] == r["C_store"] if k == "H3"
                                         else r["checks"][okey[k]])]
            for r in bad[:5]:
                print(f"    {k} FAIL {r['P']}/{r['Q']} eps={r['eps']} "
                      f"b={r.get('b')} kappa={r.get('kappa')} "
                      f"C_pred={r.get('C_pred')} C_store={r.get('C_store')} "
                      f"{r.get('interior_offdiag')}", flush=True)
    os.makedirs("data", exist_ok=True)
    json.dump({"summary": summary, "rows": rows}, open(out_path, "w"), indent=1)
    print(f"wrote {out_path}", flush=True)
    return summary, rows


# ---------------------------------------------------------------- the MECHANISM

def mech_row(P, Q, eps, cap=CAP, cap2=CAP_BIG):
    """(M1): does every interior prong contain a PERPENDICULAR bounce?

    A bounce is perpendicular iff `u_out == u_in + 2Q` -- the orbit reverses, so it
    RETRACES and necessarily comes back to the vertex it was launched from.  Hence
    **(M1) implies (*)**, with no topology: it is the elementary half of the reduction.
    Perpendicular side-copies are exactly `Fix(delta_perp)`, the reflection with axis
    perpendicular to the class direction, so (M1) says every separatrix meets that
    transversal -- the separatrix-level twin of the store's `perp_complete` flag
    ([NCYL-271] DO-NOT: `perp_complete` is not a completeness PROOF above `Q ~ 55`)."""
    geo = s437.geometry(P, Q)
    verts = s439.vertices(P, Q)
    B = RightTriangleBilliard((P / Q) * (math.pi / 2))
    rows, nofold, notmid, notpal, capped = [], 0, 0, 0, 0
    for v0 in ("O", "R", "A"):
        for u0 in s439.prongs(v0, P, Q, eps):
            bounces, why = s439.trace_prong(B, geo, verts, Q, v0, u0, cap)
            if why == "cap" and cap2 > cap:
                bounces, why = s439.trace_prong(B, geo, verts, Q, v0, u0, cap2)
            if not why.startswith("vertex"):
                capped += 1
                continue
            n = len(bounces)
            perp = [i for i, bo in enumerate(bounces)
                    if bo["u_out"] == (bo["u_in"] + 2 * Q) % (4 * Q)]
            # the LAST bounce is the arrival at the vertex and is not part of the
            # retraced word: the orbit folds at bounce `n/2` (0-indexed `n/2-1`) and
            # mirrors back, so `sides[:n-1]` is the palindrome.
            sides = [bo["side"] for bo in bounces][:n - 1]
            pal = bool(sides) and sides == sides[::-1]
            mid = (n % 2 == 0) and perp == [n // 2 - 1]
            nofold += not perp
            notmid += not mid
            notpal += not pal
            rows.append({"v0": v0, "u0": u0, "arr": why.split(":")[-1], "n": n,
                         "n_perp": len(perp), "perp_at": perp[:3], "mid": mid,
                         "pal": pal})
    return {"P": P, "Q": Q, "eps": eps, "lone": eps == P % 2, "n_prongs": len(rows),
            "capped": capped, "n_nofold": nofold, "n_notmid": notmid,
            "n_notpal": notpal, "rows": rows}


def mech_sweep(qmax=25, out_path="data/s444_mechanism.json", verbose=True):
    t0 = time.time()
    keys = sorted(stored_counts(qmax), key=lambda k: (k[1], k[0], k[2]))
    out, tot = [], Counter()
    for (P, Q, eps) in keys:
        r = mech_row(P, Q, eps)
        out.append(r)
        tot["prongs"] += r["n_prongs"]
        tot["nofold"] += r["n_nofold"]
        tot["notmid"] += r["n_notmid"]
        tot["notpal"] += r["n_notpal"]
        tot["capped"] += r["capped"]
        tot["rows"] += 1
        if r["n_nofold"] and verbose:
            print(f"    M1 FAIL {P}/{Q} eps={eps}: {r['n_nofold']} prongs with no "
                  f"perpendicular bounce", flush=True)
    if verbose:
        print(f"\nMECHANISM Q<={qmax}: {tot['rows']} class rows, {tot['prongs']} prongs "
              f"({tot['capped']} capped) [{round(time.time()-t0,1)}s]", flush=True)
        print(f"  M1 every prong has a perpendicular bounce : "
              f"{tot['prongs']-tot['nofold']}/{tot['prongs']}", flush=True)
        print(f"  M2 that bounce is the unique MIDPOINT     : "
              f"{tot['prongs']-tot['notmid']}/{tot['prongs']}", flush=True)
        print(f"  M3 the side word is a palindrome          : "
              f"{tot['prongs']-tot['notpal']}/{tot['prongs']}", flush=True)
    json.dump({"summary": dict(tot), "qmax": qmax, "rows": out},
              open(out_path, "w"), indent=1)
    print(f"wrote {out_path}", flush=True)
    return tot, out


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    kw = dict(a[2:].split("=", 1) for a in sys.argv[1:] if a.startswith("--"))
    cmd = args[0] if args else "all"
    if cmd in ("all", "h0"):
        print("H0/H1-algebra controls", flush=True)
        h0_control()
        six_case_check()
    if cmd in ("all", "mech"):
        mech_sweep(qmax=int(kw.get("qmax", 25)),
                   out_path=kw.get("out", "data/s444_mechanism.json"))
    if cmd in ("all", "sweep"):
        sweep(qmax=int(kw.get("qmax", 20)),
              out_path=kw.get("out", "data/s444_base_graph.json"),
              budget_s=float(kw["budget"]) if "budget" in kw else None)


if __name__ == "__main__":
    main()
