#!/usr/bin/env python3
"""
s315_oddP_census.py -- PRE-REGISTERED: the s276 separatrix census (COMPLETE / MIXED /
UNRESOLVED classification), extended to ODD `P` for the first time.

WHY ODD P IS NEW.  Every separatrix census this repo has ever run (`s276_separatrix_census.py`,
`s276_plateau.py`, `s276_census_sweep.py`) used EVEN-`P` rows (`8/15`, `4/15`, `2/15`, `4/13`,
`6/17`, ...).  Odd `P` is the corner where `g_len` (the perpendicular-beam swept fraction) is
measured deficient by far more: `0.299..0.735` on 9 rows `P=3..9` (`cf_width_laws.md` §s275 H_D
block, `data/s275_flow_exact.json` key `odd`), against a handful of percent at even `P`.  Nobody
has asked what KIND of set the missing area is at odd `P` -- cylinders, a minimal component, or
something else -- and this census is the instrument that can tell, exactly as it told at `8/15`
(even `P`, s276: 4/15/2/15/4/13/6/17 close 100%, `8/15` alone has a genuine minimal component).

THE DERIVATION THIS PROBE IS GATED ON.  `s276_separatrix_census.py`'s docstring and
`rulings.md` [NCYL-001]/[NCYL-006] describe its cone-point set ("two acute vertices O, A;
R regular; `P/2-1` prongs at O, `(Q-P-1)/2` at A") as the EVEN-`P` derivation
(`[NCYL-054]`: that stratum, `H(P/2-1, P/2-1, Q-P-1)`, was RE-DERIVED at even `P` because
`veech_reframing.md` `§1a`'s own "Setup" line opens "`P` odd" and turns on `gcd(P,2Q)=1`,
which FAILS for even `P`).  So the odd-`P` cone-point data is not a re-derivation at all --
it is `§1a`'s OWN PROVED table (s161, "`T-SURFACE` PROVED" in `claims.md`), read directly:

    P odd, gcd(P,Q)=1.  Q even: O order P-1 (1 pt), A order Q-P-1 (1 pt), R regular (Q pts).
                         Q odd:  O order P-1 (1 pt), A order (Q-P)/2-1 (2 pts), R regular (Q pts).
    Genus g = floor(Q/2) both parities; Sum(orders) = 2g-2 (Gauss-Bonnet, proved in §1a).

`cone_table()` below reproduces this and CHECKS Gauss-Bonnet in code before any trace runs
(`derive_and_validate`, run once per row at start of `main`) -- if it fails, the run stops
before spending any compute, per the brief's explicit stop condition.

THE GEOMETRIC ENGINE NEEDS NO CHANGE.  `s276_separatrix_census.py`'s `verts()`/`inward_dirs()`
never reference cone orders or Q-parity at all -- they compute the vertex O angle directly as
`alpha=(P/Q)(pi/2)` (the paper's `2*alpha/pi=P/Q` convention, identical in both parities) and
enumerate grid directions `k*pi/Q` strictly inside the true Euclidean vertex angle at O or A.
That is exactly the physical operation of tracing a separatrix from the ACTUAL triangle vertex,
and it is parity-agnostic by construction: the code was never given even-P-specific knowledge to
begin with, so it needs no odd-P-specific patch.  What IS new here is the INDEPENDENT check that
this purely-geometric enumeration reproduces the count implied by the (odd-P, §1a) cone-point
table -- i.e. a second, unrelated derivation of the same number, which is exactly what
"validate independently" should mean:

    interior separatrix count (from inward_dirs, geometric)  ==  ceil(Q/2) - 1

verified below for every row before it is traced (`derive_and_validate`; the real-interval
integer-counting derivation of the right-hand side is in that function's docstring -- an
earlier draft guessed `Sum(orders)/2` and that guess FAILS at `3/7`, caught by this same
check).  ⚠ This is a DIFFERENT family from [NCYL-001]'s worked total `(Q-3)/2`: that formula
is the EVEN-P family (P even forces Q odd, and there O splits instead of A), a different
stratum on the same Q, not directly comparable.  Within the odd-P family here, `ceil(Q/2)-1`
covers BOTH Q parities uniformly (`= Q/2-1` for Q even, `= (Q-1)/2` for Q odd) -- the Q-even
case is new: no prior even-P row could ever have Q even (even P + gcd(P,Q)=1 forces Q odd),
so this is the first time that branch of the geometry has been exercised at all.

R-REGULARITY AT ODD P.  `§1a`'s table gives `q_R=2`, `m_R=1` -- i.e. order 0 -- UNCONDITIONALLY
(not parity-split in its own derivation); this is confirmed by the Gauss-Bonnet check itself
(R contributes 0 to Sum(orders) in both the Q-even and Q-odd rows above, and the check passes).
So R is regular at odd P by the SAME proof that gives the rest of the table, not an assumption.

CLASSIFICATION (exactly three classes, no fourth option):
  COMPLETELY PERIODIC -- every separatrix closes (arrives at a cone point) within the cap.
  MIXED               -- >=1 separatrix open at the cap AND the closed-count ladder is FLAT
                          (equal) over the last 3 rungs -- s276_plateau.classify's PLATEAU rule.
  UNRESOLVED          -- >=1 separatrix open at the cap AND the ladder is still RISING --
                          truncation, not a deficit (s276_plateau's RISING rule, `rulings.md`
                          [NCYL-005]).

CALIBRATION GATE (must pass before any odd-P verdict is trusted): re-run, on this exact code
path, the four full-sweep ring-exact even-P controls (`4/15`,`2/15`,`4/13`,`6/17` -> must be
COMPLETELY PERIODIC) and the known deficit `8/15` (-> must be MIXED, 4/6 close, 2 never, flat
across `1e3..2e7`).  If any control fails, the run stops and reports the failure; no odd-P row
is trusted from a mis-calibrated instrument.

TOLERANCE-EROSION CHECK (`rulings.md` [NCYL-006]).  Any row that comes back COMPLETELY PERIODIC
at `tol=1e-9` is re-run at `tol=1e-11` (two orders tighter -- the floor set by the billiard
engine's own internal EPS=1e-11, `right_triangle_billiards.EPS`) on the SAME ladder.  If the
tightened run still closes 100%, the verdict is CONFIRMED; if a separatrix that closed at 1e-9
now fails to arrive at 1e-11, the original closure was a tolerance artifact and the row is
reclassified using the tightened data.  This is the hypothesis flagged in the brief as most
likely to be wrong: a false COMPLETELY PERIODIC manufactured by tolerance erosion.

PRIOR ART (grepped 'odd-P', 'odd P separatrix', 'q_A', 'q_O', '[NCYL-054]', '[NCYL-001]'..
'[NCYL-010]', 'separatrix census' against `rulings.md`, `cf_width_laws.md`, `veech_reframing.md`) ->
  - `rulings.md` [NCYL-001]..[NCYL-010]: the even-P census instrument, its calibration, its two
    named gaps (arrival tolerance, along-side prongs) -- all REUSED here, none re-derived.
  - `rulings.md` [NCYL-054]: the even-P stratum re-derivation notice -- confirms odd-P is the
    UN-swapped, already-PROVED case (`veech_reframing.md` §1a), not a case needing new work.
  - `veech_reframing.md` §1a: the odd-P cone-point/genus/stratum table, PROVED s161 -- the
    derivation this probe reads off and Gauss-Bonnet-checks in code.
  - `cf_width_laws.md` §s276, §s275 (6): the census methodology and the odd-P g_len deficit
    numbers (`0.299..0.735`) that motivate this run.
  - No prior probe runs a separatrix census at odd P.  (`probes/s186_atlas.py` /
    `f_leg_covariance.md` §20.23-20.25 work an odd-P surface at `3/7` via a translation-surface
    ATLAS, a different instrument for a different question -- not a census, and not reused here
    per `s276_separatrix_census.py`'s own note that it "stays in the billiard picture, so it
    needs no atlas".)

Run: PYTHONPATH=.:engine:archive/scripts_2026-07:archive/scripts_2026-08:probes \
     .venv/bin/python3.13 probes/s315_oddP_census.py 2>&1 | tee logs/s315_census.log
"""
import json
import math
import sys
import time

import s276_separatrix_census as sc  # reused verbatim: verts, inward_dirs, trace_separatrix

CAP = 20000000
LADDER = (1000, 10000, 100000, 1000000, 5000000, 20000000)
TOL_STD = 1e-9
TOL_TIGHT = 1e-11

CALIBRATION = [
    (4, 15, "COMPLETELY PERIODIC", "CONTROL full sweep, same Q as 8/15"),
    (2, 15, "COMPLETELY PERIODIC", "CONTROL Veech locus, same Q"),
    (4, 13, "COMPLETELY PERIODIC", "CONTROL ring-exact, non-Veech"),
    (6, 17, "COMPLETELY PERIODIC", "CONTROL ring-exact"),
    (8, 15, "MIXED", "known minimal component (s276)"),
]

PRIMARY = [
    (3, 7), (3, 8), (3, 10), (3, 11), (5, 12), (3, 13), (5, 13), (7, 16), (9, 20),
]

P1_ROWS = [(1, 5), (1, 7), (1, 9)]

# Secondary: EXHAUSTIVE odd-P, Q<=20 sweep (all P odd, 1<=P<Q, gcd(P,Q)=1), minus the rows
# already covered by PRIMARY/P1_ROWS above.  DECISION (revised in-session): an initial 43-row
# curated subset ran in well under the time budget (full run incl. calibration+primary+P1
# finished in 134s), so there is no reason to curate -- the exhaustive 72-row set is cheap
# enough to just run in full, giving a genuine "spread" deliverable instead of a sample of one.
SECONDARY = [
    (1, 4), (3, 4), (3, 5), (1, 6), (5, 6), (5, 7), (1, 8), (5, 8), (7, 8), (5, 9), (7, 9),
    (1, 10), (7, 10), (9, 10), (1, 11), (5, 11), (7, 11), (9, 11), (1, 12), (7, 12), (11, 12),
    (1, 13), (7, 13), (9, 13), (11, 13), (1, 14), (3, 14), (5, 14), (9, 14), (11, 14), (13, 14),
    (1, 15), (7, 15), (11, 15), (13, 15), (1, 16), (3, 16), (5, 16), (9, 16), (11, 16), (13, 16),
    (15, 16), (1, 17), (3, 17), (5, 17), (7, 17), (9, 17), (11, 17), (13, 17), (15, 17),
    (1, 18), (5, 18), (7, 18), (11, 18), (13, 18), (17, 18), (1, 19), (3, 19), (5, 19), (7, 19),
    (9, 19), (11, 19), (13, 19), (15, 19), (17, 19), (1, 20), (3, 20), (7, 20), (11, 20),
    (13, 20), (17, 20), (19, 20),
]

TIME_BUDGET_SECS = 4.5 * 3600  # leave margin inside the 5h hard stop
T_START = time.time()


def time_left():
    return TIME_BUDGET_SECS - (time.time() - T_START)


def cone_table(P, Q):
    """The cone-point table, read directly from `veech_reframing.md` §1a (PROVED s161).
    Returns (cones, r_regular, g) where cones lists {vertex, order, count} for the SINGULAR
    points only (R, order 0, is reported separately -- regular, never a separatrix endpoint).

    ⚠ GENERALISED TO EVEN `P` at s316 (user-raised: the whole s315 cylinder census is odd `P`
    because THIS function asserted it).  The odd-`P` branches below are byte-identical to the
    s315 version; the even-`P` branch is new.  §1a's rule: write each vertex angle as a
    fraction of pi in lowest terms `m/q`; the class splits into `N/q = 2Q/q` cone points, each
    of cone angle `2*pi*m`, i.e. order `m-1`.  Angles are `O = P/(2Q)`, `R = 1/2`,
    `A = (Q-P)/(2Q)`.
      * `P` ODD  ⟹ gcd(P,2Q)=1 ⟹ O stays `P/(2Q)`: ONE point of order `P-1`.
        `A`: gcd(Q-P,2Q)=gcd(Q-P,2)`, so Q even ⟹ one point of order `Q-P-1`;
        Q odd ⟹ TWO points of order `(Q-P)/2 - 1`.
      * `P` EVEN ⟹ `Q` is ODD (coprimality), gcd(P,2Q)=2, so `O = (P/2)/Q`: **TWO** points of
        order `P/2 - 1`; and `Q-P` is odd so gcd(Q-P,2Q)=1 and `A` is ONE point of order
        `Q-P-1`.  Exactly the odd-`P` odd-`Q` case with the roles of O and A EXCHANGED --
        which is the same leg asymmetry that makes the s316 leg-swap law even-`Q` only.
    Gauss-Bonnet holds in every branch (`sum orders == 2g-2`); `derive_and_validate` checks it
    and cross-checks the separatrix count against the independent geometric `inward_dirs`."""
    assert math.gcd(P, Q) == 1, f"P={P}, Q={Q} not coprime"
    g = Q // 2
    if P % 2 == 1:
        if Q % 2 == 0:
            cones = [
                {"vertex": "O", "order": P - 1, "count": 1},
                {"vertex": "A", "order": Q - P - 1, "count": 1},
            ]
        else:
            cones = [
                {"vertex": "O", "order": P - 1, "count": 1},
                {"vertex": "A", "order": (Q - P) // 2 - 1, "count": 2},
            ]
    else:
        assert Q % 2 == 1, f"P even forces Q odd by coprimality; got P={P}, Q={Q}"
        cones = [
            {"vertex": "O", "order": P // 2 - 1, "count": 2},
            {"vertex": "A", "order": Q - P - 1, "count": 1},
        ]
    r_regular = {"vertex": "R", "order": 0, "count": Q}
    return cones, r_regular, g


def derive_and_validate(P, Q):
    """Cross-check the §1a cone table (Gauss-Bonnet) against the INDEPENDENT geometric
    enumeration `inward_dirs` uses (no reference to cone orders at all).  Returns a dict;
    raises AssertionError (caught by the caller, which STOPS the run) on any mismatch.

    THE SECOND FORMULA (derived here, not asserted): the interior separatrix count is NOT
    simply `sum_orders/2` in general -- an early draft of this probe assumed that and it
    FAILED at `3/7` (predicted 2, actual 3).  Redone from the vertex-angle geometry directly:
    `inward_dirs` counts integers strictly inside an open real interval of width `w` starting
    at a point `h`.  At O, `h=0` (an exact grid point) and `w=P/2` (a half-integer, P odd),
    giving `floor((P-1)/2) = (P-1)/2` interior integers, independent of Q's parity.  At A,
    `h=Q+P/2` (always a half-integer, since P is odd) and `w=(Q-P)/2`:
      - Q even: Q-P is odd, so w is a half-integer too and h+w is an INTEGER -- the open
        interval (half-integer, integer) contains `floor(w) = (Q-P-1)/2` integers.
      - Q odd:  Q-P is even, so w is an integer and h+w is again a half-integer -- the open
        interval (half-integer, half-integer) contains exactly `w = (Q-P)/2` integers.
    Summing: total = (P-1)/2 + (Q-P-1)/2 = Q/2-1 = g-1 (Q even); total = (P-1)/2+(Q-P)/2
    = (Q-1)/2 = g (Q odd).  Unified: **total = ceil(Q/2) - 1** in both parities.  This is a
    genuinely SECOND derivation (real-interval integer-counting, not stratum theory) and it
    reproduces every one of the 12 hand/numeric spot-checks run while building this probe."""
    cones, r_regular, g = cone_table(P, Q)
    sum_orders = sum(c["order"] * c["count"] for c in cones)  # R contributes 0
    gauss_bonnet_ok = (sum_orders == 2 * g - 2)

    _, _, L, alpha = sc.verts(P, Q)
    geo_O = len(sc.inward_dirs(P, Q, "O", L))
    geo_A = len(sc.inward_dirs(P, Q, "A", L))
    geo_total = geo_O + geo_A
    # ⚠ s316: the closed form `ceil(Q/2) - 1` in the docstring above is the ODD-`P` value
    # (its derivation says "h = Q + P/2 (always a half-integer, since P is odd)").  Replaced
    # by the SAME interval arithmetic done exactly, which is parity-general: `inward_dirs`
    # keeps k with theta = k*pi/Q strictly inside the vertex cone, i.e. k strictly inside
    # (0, P/2) at O and strictly inside (Q + P/2, 3Q/2) at A.  Fractions, so no float edge.
    from fractions import Fraction as _F
    def _open_count(a, b):
        """# integers strictly between a and b (exact rationals)."""
        lo = math.floor(a) + 1
        hi = math.ceil(b) - 1
        if a == lo - 1: pass                 # a integral -> lo already excludes it
        if b == hi + 1: pass                 # b integral -> hi already excludes it
        return max(0, hi - lo + 1)
    pred_O = _open_count(_F(0), _F(P, 2))
    pred_A = _open_count(_F(Q) + _F(P, 2), _F(3 * Q, 2))
    predicted_total = pred_O + pred_A
    # closed forms, for the record: odd P -> ceil(Q/2)-1 ; even P (hence Q odd) -> (Q-3)/2
    closed_form = (-(-Q // 2) - 1) if P % 2 else ((Q - 3) // 2)
    geometric_ok = (geo_total == predicted_total)

    return {
        "P": P, "Q": Q, "genus": g,
        "cones": cones, "R": r_regular,
        "sum_orders": sum_orders,
        "gauss_bonnet_2g_minus_2": 2 * g - 2,
        "gauss_bonnet_ok": gauss_bonnet_ok,
        "geo_count_O": geo_O, "geo_count_A": geo_A, "geo_count_total": geo_total,
        "predicted_total_ceilQ2_minus_1": predicted_total,   # name kept for callers
        "predicted_O": pred_O, "predicted_A": pred_A,
        "closed_form_total": closed_form,
        "closed_form_ok": (predicted_total == closed_form),
        "geometric_cross_check_ok": geometric_ok,
        "valid": gauss_bonnet_ok and geometric_ok,
    }


def run_row_tol(P, Q, cap, ladder, tol):
    """Like `sc.run_row` but with an explicit arrival tolerance (sc.run_row/trace_separatrix's
    default `tol=ARRIVE_TOL` is bound at function-definition time, so it cannot be overridden
    by mutating the module global -- this wrapper passes `tol` explicitly instead)."""
    cone, reg, L, alpha = sc.verts(P, Q)
    B = sc.RightTriangleBilliard(alpha)  # sc imports this name at module scope
    recs = []
    for v0 in ("O", "A"):
        for k0 in sc.inward_dirs(P, Q, v0, L):
            r = sc.trace_separatrix(B, P, Q, v0, k0, cap, cone, reg, tol=tol)
            r.update({"from": v0, "k": k0})
            recs.append(r)
    closed = [r for r in recs if r["arrived"]]
    curve = [{"depth": d,
              "n_closed": sum(1 for r in closed if r["hits"] <= d),
              "frac_closed": sum(1 for r in closed if r["hits"] <= d) / max(1, len(recs))}
             for d in ladder if d <= cap]
    return {"P": P, "Q": Q, "cap": cap, "tol": tol, "n_separatrices": len(recs),
            "n_closed": len(closed), "frac_closed": len(closed) / max(1, len(recs)),
            "curve": curve,
            "worst_miss_of_open": max([r["best_miss"] for r in recs
                                       if not r["arrived"]], default=None),
            "min_miss_of_open": min([r["best_miss"] for r in recs
                                     if not r["arrived"]], default=None),
            "near_R_total": sum(r["near_R"] for r in recs),
            "records": recs}


def classify(row_result):
    """PLATEAU/RISING rule from `probes/s276_plateau.py::classify`, plus the ALL_CLOSED case."""
    if row_result["n_closed"] == row_result["n_separatrices"]:
        return "COMPLETELY PERIODIC"
    counts = [c["n_closed"] for c in row_result["curve"]]
    if len(counts) < 3:
        return "UNRESOLVED"  # not enough ladder rungs to call a plateau; conservative
    if counts[-1] == counts[-2] == counts[-3] and counts[-1] < row_result["n_separatrices"]:
        return "MIXED"
    return "UNRESOLVED"


def census_row(P, Q, tag, cap=CAP, ladder=LADDER, tol=TOL_STD):
    t0 = time.time()
    r = run_row_tol(P, Q, cap, ladder, tol)
    r["tag"] = tag
    r["secs"] = round(time.time() - t0, 1)
    r["verdict"] = classify(r)
    opens = [x for x in r["records"] if not x["arrived"]]
    r["open_dirs"] = [{"from": x["from"], "k": x["k"], "best_miss": x["best_miss"],
                        "hits_at_cap": x["hits"]} for x in opens]
    r.pop("records", None)
    return r


def maybe_tighten(P, Q, base_result, ladder=LADDER, cap=CAP):
    """[NCYL-006] guard: if a row reads COMPLETELY PERIODIC at TOL_STD, re-run at TOL_TIGHT
    (two orders tighter) to check the closure is not tolerance erosion faking a saddle
    connection.  Returns (tight_result_or_None, final_verdict, note)."""
    if base_result["verdict"] != "COMPLETELY PERIODIC":
        return None, base_result["verdict"], "no tightening needed (not COMPLETE at std tol)"
    tight = census_row(P, Q, base_result["tag"] + " [TIGHTENED 1e-11]", cap, ladder, TOL_TIGHT)
    if tight["verdict"] == "COMPLETELY PERIODIC":
        return tight, "COMPLETELY PERIODIC", "CONFIRMED at tol=1e-11: not tolerance erosion"
    return tight, tight["verdict"], (
        f"REVISED: closed 100%% at tol=1e-9 but only "
        f"{tight['n_closed']}/{tight['n_separatrices']} at tol=1e-11 -- "
        f"the 1e-9 closure(s) were tolerance artifacts")


def dump(out, path):
    with open(path, "w") as f:
        json.dump(out, f, indent=1)
    print(f"wrote {path}", flush=True)


def main():
    decisions = [
        "secondary row list is the FULL exhaustive 72-row Q<=20 odd-P sweep (revised "
        "in-session from an initial 43-row curated subset once the run time was measured "
        "at 134s total -- see SECONDARY comment in source).",
        "UNRESOLVED is the conservative default when a row has <3 ladder rungs available "
        "(cap smaller than the 3rd-from-last ladder rung) -- never called MIXED without a "
        "3-rung-flat check.",
        "tolerance-erosion recheck runs ONLY on COMPLETELY PERIODIC verdicts (the class the "
        "brief flags as most likely to be a false positive), not on MIXED/UNRESOLVED -- a "
        "tighter tolerance can only ever turn an arrival into a non-arrival, never the reverse, "
        "so MIXED/UNRESOLVED rows would not change under tightening and re-running them would "
        "just burn the CAP a second time for no new information.",
    ]
    out = {"decisions": decisions, "derivation": [], "calibration": [], "primary": [],
           "p1_rows": [], "secondary": [], "calibration_gate_passed": None,
           "stopped_early": False, "wall_clock_secs": None}

    print("=== DERIVATION + GAUSS-BONNET / GEOMETRIC CROSS-CHECK (odd-P rows only) ===",
          flush=True)
    derivation_ok = True
    for (P, Q) in PRIMARY + P1_ROWS + SECONDARY:
        d = derive_and_validate(P, Q)
        out["derivation"].append(d)
        status = "OK" if d["valid"] else "FAIL"
        print(f"  {P:>2}/{Q:<3} g={d['genus']:<2} cones={d['cones']} "
              f"sum_orders={d['sum_orders']:<3} 2g-2={d['gauss_bonnet_2g_minus_2']:<3} "
              f"geo_total={d['geo_count_total']:<3} predicted={d['predicted_total_ceilQ2_minus_1']:<3} "
              f"[{status}]", flush=True)
        if not d["valid"]:
            derivation_ok = False
    dump(out, "data/s315_oddP_census.json")
    if not derivation_ok:
        print("DERIVATION FAILED GAUSS-BONNET / GEOMETRIC CROSS-CHECK -- STOPPING, "
              "no traces run.", flush=True)
        out["stopped_early"] = True
        out["wall_clock_secs"] = round(time.time() - T_START, 1)
        dump(out, "data/s315_oddP_census.json")
        return

    print("\n=== CALIBRATION GATE (even-P controls, must reproduce known verdicts) ===",
          flush=True)
    gate_ok = True
    for (P, Q, expect, tag) in CALIBRATION:
        r = census_row(P, Q, tag)
        r["expected"] = expect
        r["gate_pass"] = (r["verdict"] == expect)
        gate_ok = gate_ok and r["gate_pass"]
        out["calibration"].append(r)
        print(f"  {P:>2}/{Q:<3} {tag:<32} verdict={r['verdict']:<20} expect={expect:<20} "
              f"{'PASS' if r['gate_pass'] else 'FAIL'} sep={r['n_separatrices']} "
              f"closed={r['n_closed']} [{r['secs']}s]", flush=True)
    out["calibration_gate_passed"] = gate_ok
    dump(out, "data/s315_oddP_census.json")
    if not gate_ok:
        print("CALIBRATION GATE FAILED -- STOPPING, no odd-P rows will be run.", flush=True)
        out["stopped_early"] = True
        out["wall_clock_secs"] = round(time.time() - T_START, 1)
        dump(out, "data/s315_oddP_census.json")
        return

    print("\n=== PRIMARY ODD-P ROWS ===", flush=True)
    for (P, Q) in PRIMARY:
        if time_left() < 60:
            out["stopped_early"] = True
            print("TIME BUDGET NEAR EXHAUSTION -- stopping before primary row "
                  f"{P}/{Q}.", flush=True)
            break
        r = census_row(P, Q, "primary")
        tight, final_verdict, note = maybe_tighten(P, Q, r)
        r["tolerance_check"] = {"tight_result": tight, "note": note}
        r["final_verdict"] = final_verdict
        out["primary"].append(r)
        dump(out, "data/s315_oddP_census.json")
        print(f"  {P:>2}/{Q:<3} verdict={r['verdict']:<20} final={final_verdict:<20} "
              f"sep={r['n_separatrices']} closed={r['n_closed']} "
              f"min_miss_open={r['min_miss_of_open']} [{r['secs']}s] :: {note}", flush=True)

    print("\n=== P=1 ROWS (extreme deficit corner, g_len=2/Q exactly) ===", flush=True)
    for (P, Q) in P1_ROWS:
        if time_left() < 60:
            out["stopped_early"] = True
            break
        r = census_row(P, Q, "P=1 control")
        tight, final_verdict, note = maybe_tighten(P, Q, r)
        r["tolerance_check"] = {"tight_result": tight, "note": note}
        r["final_verdict"] = final_verdict
        out["p1_rows"].append(r)
        dump(out, "data/s315_oddP_census.json")
        print(f"  {P:>2}/{Q:<3} verdict={r['verdict']:<20} final={final_verdict:<20} "
              f"sep={r['n_separatrices']} closed={r['n_closed']} [{r['secs']}s] :: {note}",
              flush=True)

    print("\n=== SECONDARY SPREAD (if time allows) ===", flush=True)
    for (P, Q) in SECONDARY:
        if time_left() < 30:
            out["stopped_early"] = True
            print(f"TIME BUDGET NEAR EXHAUSTION -- stopping secondary sweep before "
                  f"{P}/{Q}; {len(out['secondary'])}/{len(SECONDARY)} secondary rows done.",
                  flush=True)
            break
        r = census_row(P, Q, "secondary")
        tight, final_verdict, note = maybe_tighten(P, Q, r)
        r["tolerance_check"] = {"tight_result": tight, "note": note}
        r["final_verdict"] = final_verdict
        out["secondary"].append(r)
        dump(out, "data/s315_oddP_census.json")
        print(f"  {P:>2}/{Q:<3} verdict={r['verdict']:<20} final={final_verdict:<20} "
              f"sep={r['n_separatrices']} closed={r['n_closed']} [{r['secs']}s] :: {note}",
              flush=True)

    out["wall_clock_secs"] = round(time.time() - T_START, 1)
    dump(out, "data/s315_oddP_census.json")
    print(f"\nDONE. wall_clock={out['wall_clock_secs']}s stopped_early={out['stopped_early']}",
          flush=True)


if __name__ == "__main__":
    sys.exit(main())
