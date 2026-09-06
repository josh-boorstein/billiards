"""PRE-REGISTERED — s332, item (a) first move: is [NCYL-123]'s vertex gap real?

CLAIM UNDER TEST (the READ's conclusion, asserted here in code rather than in prose):
  `orphan_theorem.md` Lemma B's step "a velocity reversal ... occurs ONLY at a
  perpendicular wall hit" is applied only at launch heights INTERIOR to a cell, and a
  cell interior is vertex-free BY THE DEFINITION of the partition (`foundations.md` §1:
  a cell is a maximal s-interval of constant word; a cell boundary is a zero of a
  vertex-graze function η_{k,V}).  Hence the vertex alternative is excluded by the
  DOMAIN of T, not by an unstated hypothesis.

FOUR ASSERTIONS, each of which could come out otherwise ([OPS-041]):
  A  orphan-interior orbits are vertex-free with a MACROSCOPIC margin
     (min over all reflections of dist-to-nearest-vertex >> 0).
  B  the turnaround (head-on) hit is in the INTERIOR of its side — the thing
     [NCYL-072]/[NCYL-080] flag as unexcluded.
  C  NEGATIVE CONTROL, the one that makes A/B non-vacuous: at the orphan cell's two
     ENDPOINTS the same margin collapses to O(δ), i.e. the boundary IS a vertex graze.
     Without C, A is just "a generic height misses finitely many points".
  D  WHICH vertex each orphan endpoint grazes.  R (right angle) is the interesting
     answer: R has cone angle 2π (REGULAR — `computational_findings.md` §101) and is a
     2D corner retroreflector (R_L1·R_L2 = rotation by π ⇒ v ↦ −v with NO perpendicular
     wall hit), so R is the one vertex that can actually produce the reversal Lemma B's
     "ONLY" would forbid.  If R-grazes sit on cell boundaries, the escape is real and
     exactly excluded.

PRIOR ART: grepped 'vertex gap', 'NCYL-072', 'NCYL-080', 'NCYL-123', 'retroreflector',
  'cone angle 2π', 'first-return' →
   [NCYL-123] (s328) flags the gap's textual home as Lemma B and leaves THIS read undone;
   [NCYL-072]/[NCYL-080] carry it as "degenerate turnaround at a VERTEX (saddle
   connection), unexcluded";  `computational_findings.md` §101 gives R cone angle 2π
   (regular, so NOT a saddle connection — the gap's own phrasing mis-names it);
   `cf_width_laws.md` §s218/§s219/§s220 own the R-retroreflector / J-fold-centre story
   in the α→0 strand ([T-ESELECT] gap-C, criterion num_V | Q) — a DIFFERENT strand, cited
   as corroboration only, not imported.  No probe anywhere measures dist-to-vertex along
   an orphan orbit.

FLOAT, deliberately: the assertions are macroscopic separations (O(0.1) vs O(1e-9)), not
fine measurements.  Explore with float, certify with the ring — nothing here needs the ring.
"""
import json
import math
from math import gcd

from right_triangle_billiards import RightTriangleBilliard

TOL_RETURN = 1e-9      # |T(s) - s| for "retraces"
TOL_HEADON = 1e-9      # |v . d_side| for "head-on"
DELTA = 1e-9           # how far inside/outside a boundary we probe


def _ctx(P, Q):
    return RightTriangleBilliard((P / Q) * math.pi / 2)


def _side_dirs(bil):
    return {"L1": (0.0, 1.0), "L2": (1.0, 0.0), "H": bil.axis_H}


def orbit_report(bil, s):
    """Trace the first-return orbit at height s; report vertex margins + turnaround."""
    tr = bil.trace(s, max_hits=20000)
    if not tr.hits:
        return None
    verts = {"O": bil.O, "R": bil.R, "A": bil.A}
    dirs = _side_dirs(bil)

    min_d, min_v, min_i = float("inf"), None, None
    headons = []
    for h in tr.hits:
        for name, (vx, vy) in verts.items():
            d = math.hypot(h.x - vx, h.y - vy)
            if d < min_d:
                min_d, min_v, min_i = d, name, h.index
        dx, dy = dirs[h.side]
        tang = abs(math.cos(h.incoming_angle) * dx + math.sin(h.incoming_angle) * dy)
        if tang < TOL_HEADON and h.index < len(tr.hits):   # interior head-on
            dv = min(math.hypot(h.x - v[0], h.y - v[1]) for v in verts.values())
            headons.append({"index": h.index, "side": h.side, "dist_to_vertex": dv})

    return {
        "s": s,
        "nhits": len(tr.hits),
        "returned": tr.returned_perpendicular,
        "T_s": tr.hits[-1].y if tr.returned_perpendicular else None,
        "min_dist_to_vertex": min_d,
        "argmin_vertex": min_v,
        "argmin_hit": min_i,
        "word_len": len(tr.full_word),
        "headons": headons,
    }


def is_orphan(bil, s):
    r = orbit_report(bil, s)
    return r is not None and r["returned"] and abs(r["T_s"] - s) < TOL_RETURN


def find_orphan_cell(bil, ngrid=4000):
    """Grid-scan for a retracing height, then bisect both cell endpoints on retrace."""
    seed = None
    for i in range(1, ngrid):
        s = i / ngrid
        if is_orphan(bil, s):
            seed = s
            break
    if seed is None:
        return None

    def push(direction):
        lo, hi = seed, seed + direction * 1.0
        hi = min(max(hi, 1e-12), 1 - 1e-12)
        # hi is outside (or the interval end); bisect on the retrace predicate
        if is_orphan(bil, hi):
            return hi
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            if is_orphan(bil, mid):
                lo = mid
            else:
                hi = mid
        return lo

    return push(-1), push(+1)


def run(centres):
    rows = []
    for P, Q in centres:
        bil = _ctx(P, Q)
        cell = find_orphan_cell(bil)
        if cell is None:
            rows.append({"P": P, "Q": Q, "status": "no-orphan-found"})
            continue
        lo, hi = cell
        mid = 0.5 * (lo + hi)

        interior = orbit_report(bil, mid)
        # C: probe just INSIDE each endpoint
        edge_lo = orbit_report(bil, lo + DELTA * 0.5)
        edge_hi = orbit_report(bil, hi - DELTA * 0.5)

        rows.append({
            "P": P, "Q": Q, "status": "ok",
            "cell": [lo, hi], "width": hi - lo,
            "A_interior_margin": interior["min_dist_to_vertex"],
            "A_interior_word_len": interior["word_len"],
            "B_headons": interior["headons"],
            "C_lo_margin": edge_lo["min_dist_to_vertex"],
            "C_hi_margin": edge_hi["min_dist_to_vertex"],
            "D_lo_vertex": edge_lo["argmin_vertex"],
            "D_hi_vertex": edge_hi["argmin_vertex"],
            "D_lo_hit": edge_lo["argmin_hit"], "D_lo_nhits": edge_lo["nhits"],
            "D_hi_hit": edge_hi["argmin_hit"], "D_hi_nhits": edge_hi["nhits"],
        })
    return rows


if __name__ == "__main__":
    centres = [(P, Q) for Q in range(5, 20) for P in range(1, Q)
               if P % 2 == 1 and gcd(P, Q) == 1]
    rows = run(centres)
    with open("data/s332_vertexfree.json", "w") as f:
        json.dump(rows, f, indent=1)

    ok = [r for r in rows if r["status"] == "ok"]
    print(f"centres: {len(rows)}  orphan located: {len(ok)}")
    print(f"{'P/Q':>7} {'width':>10} {'A margin':>10} {'B ho-dist':>10} "
          f"{'C lo':>9} {'C hi':>9}  D(lo,hi)")
    for r in ok:
        ho = min((h["dist_to_vertex"] for h in r["B_headons"]), default=float("nan"))
        print(f"{r['P']}/{r['Q']:>4} {r['width']:10.3e} {r['A_interior_margin']:10.3e} "
              f"{ho:10.3e} {r['C_lo_margin']:9.2e} {r['C_hi_margin']:9.2e}  "
              f"({r['D_lo_vertex']},{r['D_hi_vertex']})")

    import statistics
    print("\n--- assertions ---")
    A = min(r["A_interior_margin"] for r in ok)
    B = min(min((h["dist_to_vertex"] for h in r["B_headons"]), default=float("inf"))
            for r in ok)
    nho = [len(r["B_headons"]) for r in ok]
    C = max(max(r["C_lo_margin"], r["C_hi_margin"]) for r in ok)
    print(f"A  worst orphan-interior vertex margin      : {A:.4e}   (want >> 0)")
    print(f"B  worst turnaround dist-to-vertex          : {B:.4e}   (want >> 0)")
    print(f"   interior head-on count per orphan        : {sorted(set(nho))}")
    print(f"C  worst endpoint margin at delta={DELTA:.0e}     : {C:.4e}   (want ~ delta)")
    from collections import Counter
    print(f"D  grazed vertex at orphan endpoints        : "
          f"{dict(Counter([r['D_lo_vertex'] for r in ok] + [r['D_hi_vertex'] for r in ok]))}")
