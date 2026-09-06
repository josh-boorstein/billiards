#!/usr/bin/env python3
"""
s276_separatrix_census.py — PRE-REGISTERED: the COMPLETE instrument for
`future_directions.md` B9 q5 (a).  Sampling can never prove "there is no other
cylinder"; a separatrix census can, because the objects it enumerates are finite.

THE ARGUMENT THIS IMPLEMENTS.  On a translation surface, every cylinder in a given
direction is bounded by saddle connections in that direction, and every boundary
component of a cylinder contains at least one SEPARATRIX -- a leaf emanating from a cone
point.  In a fixed direction the number of outgoing separatrices is FINITE and known from
the stratum: a singularity of order d has d+1 outgoing prongs, so at 8/15, where the
even-P stratum is H(3,3,6) (`cf_width_laws.md` §s272 (3), `rulings.md` [NCYL]), there are
4 + 4 + 7 = 15 of them.  Hence:

    every cylinder in the perpendicular direction is bounded by one of finitely many
    separatrices  =>  trace them all, and no cylinder can hide.

That is the completeness the flux-sampling probe (`probes/s276_periodic_mass.py`) cannot
have: sampling misses a cylinder that is thinner than the sample spacing, and the
topological bound only says an unmet cylinder would have area >= (1 - g_len)/6 = 9.1% of
the surface, which does NOT bound its width.  A census has no such gap.

IN THE BILLIARD PICTURE.  The cone points are the two ACUTE vertices: O (angle alpha) and
A (angle pi/2 - alpha).  The right-angle vertex R is REGULAR -- it unfolds to cone angle
2*pi (`rulings.md` [NCYL], "the right angle unfolds to cone angle 2pi, so it cuts the WORD
without bounding a cylinder") -- so a leaf through R is not a separatrix and must not be
counted as a saddle-connection endpoint.  A separatrix is therefore a trajectory launched
FROM O or A in a grid direction, and it is a SADDLE CONNECTION exactly when it arrives at
O or A.

  H2 (minimal component, s275 (6)):  some separatrices never arrive at a cone point, at
      any depth -- they are dense in the minimal component.
  H1 (missing cylinders):  every separatrix closes into a saddle connection, and the
      saddle connections bound more cylinders than the 3 the beam meets.

THE CONTROL IS WHAT MAKES IT SOUND.  On a completely periodic direction EVERY separatrix
is a saddle connection.  So the run is only worth reading if the full-sweep control rows
at the same Q (4/15, 2/15) return 100% closure at a depth where 8/15 does not.  Reported
side by side, same depth, same code path.

GUARDS.
 G1  A separatrix that has not arrived by the cap is NOT proved aperiodic -- the depth
     curve is reported, and only a PLATEAU against a saturating control is read as
     evidence (`rulings.md` [NCYL] truncation guard, s275 (7) / s274 (7)(ii)).
 G2  Arrival is a float test, so every claimed saddle connection reports its miss
     distance; a near-miss at 1e-6 is not an arrival.  Both tolerances are swept.
 G3  Passages near the REGULAR vertex R are counted separately and never treated as
     arrivals -- confusing them would manufacture saddle connections that bound nothing.
 G4  Launch directions are restricted to those pointing into the triangle from the
     vertex, checked against the two adjacent side directions (the `s276_drift_guard`
     first draft launched a perpendicular OUTWARD from L1 and its own control caught it;
     the same class of error is guarded here explicitly).
 G6  Drift horizon: measured at `probes/s276_drift_guard.py` -- every recurrence of a
     known periodic orbit is detected out to 1.28e6 hits with worst return error 6.2e-10,
     three orders inside the 1e-7 tolerance.  So non-arrival at these depths is not drift.

PRIOR ART (grepped 'separatrix', 'saddle connection', 'prong', 'cone point', 'minimal
component', 'cylinder diagram'; read `rulings.md` [NCYL] and `f_leg_covariance.md`
§20.23-20.25 before designing) ->
  - `f_leg_covariance.md` §20.23 G1 OWNS the prong-count identity used above ("a
    singularity of order d has d+2 horizontal separatrix prongs" in its half-translation
    convention, #SC = (sum prongs)/2) and §20.25 R3 owns an EXACT separatrix tracer with
    the two traps this probe inherits: dedup by SURFACE identity rather than developed
    coordinates, and DISCARD spurious saddle connections that assign to no cylinder
    (s187 R3 found pole-launch artifacts that way).  Those sections work the ODD-P
    surface at 3/7 via the `probes/s186_atlas.py` translation atlas; this probe stays in
    the billiard picture, so it needs no atlas and no even-P re-derivation of the gluing.
  - `cf_width_laws.md` §s275 (6) owns the finding under test; §s272 (1) owns the
    regular-vertex fact that G3 implements; §s273 (7) owns the stratum and genus.
  - `rulings.md` [NCYL]: moduli commensurability is DEAD as a discriminant -- not retried
    here; the cylinder-diagram BUILD is struck -- this probe does not assemble one, it
    only counts arrivals.
  - No prior probe runs a separatrix census at even P.

Run: PYTHONPATH=.:engine:archive/scripts_2026-07:archive/scripts_2026-08:probes \
     .venv/bin/python3.13 probes/s276_separatrix_census.py
"""
import json
import math
import time

from right_triangle_billiards import RightTriangleBilliard

CAP = 2000000
LADDER = (1000, 10000, 100000, 1000000, 2000000)
ARRIVE_TOL = 1e-9
ROWS = [
    (8, 15, "DEFICIT"),
    (4, 15, "CONTROL full sweep, same Q"),
    (2, 15, "CONTROL Veech locus, same Q"),
    (4, 13, "CONTROL ring-exact, non-Veech"),
    (6, 17, "CONTROL ring-exact"),
]


def verts(P, Q):
    alpha = (P / Q) * (math.pi / 2)
    L = 1.0 / math.tan(alpha)
    # O = acute vertex (angle alpha), A = acute vertex (angle pi/2 - alpha),
    # R = the RIGHT-ANGLE vertex, which is a REGULAR point of the surface (G3).
    return {"O": (0.0, 0.0), "A": (L, 1.0)}, {"R": (L, 0.0)}, L, alpha


def inward_dirs(P, Q, vname, L):
    """Grid directions strictly inside the triangle's cone at vertex `vname` (G4)."""
    if vname == "O":
        lo, hi = 0.0, math.atan2(1.0, L)          # between L2 (+x) and H
    else:                                          # A: between H (pointing back to O)
        lo, hi = math.pi + math.atan2(1.0, L), 1.5 * math.pi   # and L1 (-y)
    out = []
    for k in range(2 * Q):
        th = (k * math.pi / Q) % (2 * math.pi)
        if lo + 1e-12 < th < hi - 1e-12:
            out.append(k)
    return out


def trace_separatrix(B, P, Q, v0, k0, cap, cone, reg, tol=ARRIVE_TOL):
    """Launch from cone point `v0` in grid direction k0; return arrival data."""
    px, py = cone[v0]
    th = k0 * math.pi / Q
    vx, vy = math.cos(th), math.sin(th)
    near_R = 0
    best = (1e9, None, 0)
    for n in range(1, cap + 1):
        cand = B._candidate_times(px, py, vx, vy)
        if not cand:
            return {"arrived": None, "hits": n, "reason": "escaped",
                    "near_R": near_R, "best_miss": best[0]}
        t, side = min(cand, key=lambda z: z[0])
        px += t * vx
        py += t * vy
        for nm, (qx, qy) in cone.items():
            d = math.hypot(px - qx, py - qy)
            if d < best[0]:
                best = (d, nm, n)
            if d < tol:
                return {"arrived": nm, "hits": n, "reason": "saddle_connection",
                        "near_R": near_R, "best_miss": d}
        for _, (qx, qy) in reg.items():             # G3: R is regular, never an arrival
            if math.hypot(px - qx, py - qy) < tol:
                near_R += 1
        vx, vy = B.reflect_velocity(side, vx, vy)
    return {"arrived": None, "hits": cap, "reason": "cap", "near_R": near_R,
            "best_miss": best[0], "best_at": best[2]}


def run_row(P, Q, cap=CAP):
    cone, reg, L, alpha = verts(P, Q)
    B = RightTriangleBilliard(alpha)
    recs = []
    for v0 in ("O", "A"):
        for k0 in inward_dirs(P, Q, v0, L):
            r = trace_separatrix(B, P, Q, v0, k0, cap, cone, reg)
            r.update({"from": v0, "k": k0})
            recs.append(r)
    closed = [r for r in recs if r["arrived"]]
    curve = [{"depth": d,
              "n_closed": sum(1 for r in closed if r["hits"] <= d),
              "frac_closed": sum(1 for r in closed if r["hits"] <= d) / max(1, len(recs))}
             for d in LADDER if d <= cap]
    return {"P": P, "Q": Q, "cap": cap, "n_separatrices": len(recs),
            "n_closed": len(closed), "frac_closed": len(closed) / max(1, len(recs)),
            "curve": curve,
            "worst_miss_of_open": max([r["best_miss"] for r in recs
                                       if not r["arrived"]], default=None),
            "min_miss_of_open": min([r["best_miss"] for r in recs
                                     if not r["arrived"]], default=None),
            "near_R_total": sum(r["near_R"] for r in recs),
            "records": recs}


def main():
    out = []
    for (P, Q, tag) in ROWS:
        t0 = time.time()
        r = run_row(P, Q)
        r["tag"] = tag
        r["secs"] = round(time.time() - t0, 1)
        out.append(r)
        print(f"{P:>2}/{Q:<3} {tag:<28} separatrices={r['n_separatrices']:<3} "
              f"closed={r['n_closed']:<3} frac={r['frac_closed']:.3f} "
              f"min_miss_open={r['min_miss_of_open']} nearR={r['near_R_total']} "
              f"[{r['secs']}s]", flush=True)
        for c in r["curve"]:
            print(f"        d={c['depth']:<9} closed={c['n_closed']:<3} "
                  f"frac={c['frac_closed']:.3f}", flush=True)
    with open("data/s276_separatrix_census.json", "w") as f:
        json.dump(out, f, indent=1)
    print("wrote data/s276_separatrix_census.json", flush=True)


if __name__ == "__main__":
    main()
