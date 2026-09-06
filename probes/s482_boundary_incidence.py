#!/usr/bin/env python3
"""
s482_boundary_incidence.py -- PRE-REGISTERED: queue item (s481-a), PROVE [NCYL-340].

THE ITEM, from the s481 handoff: "(s481-a) PROVE [NCYL-340]: the edge identity is exact
on `240/240` and has no derivation; it is a length/edge bookkeeping statement about an
object whose two sides are both exact."  [NCYL-340]'s s481 layer already split the status
-- the COUNT half is derivable arithmetic (`24460/24460` to `Q = 200`), so what is open is
the LENGTH half

      sum_c circ_c  =  sum_{prongs} ell_p  +  sum_{sigma : u_side(sigma) == p mod 2} |sigma|

⇒⇒ THE DERIVATION, and it is an ASSEMBLY of things the repo already proved plus one
standard lemma, not a new measurement.  Let `B` be the class's genus-0 half-translation
base ([NCYL-280], `B = Sigma/tau`), on which the class direction is a foliation.
  (i)   [NCYL-280], PROVED: in the class direction `B`'s separatrix diagram has `E = Q`
        edges and `2Q` ray-ends -- `P` at `Z_O`, `Q-P` at `Z_A`, `1` at each of the `Q`
        simple poles (the `R`-copies).
  (ii)  [NCYL-281], PROVED (its "boundary rays" clause): those `2Q` ray-ends are exactly
        the `2` rays of each billiard prong SLOT (the forward ray in sheet `u_0` and the
        reverse ray in sheet `u_0 + 2Q`) together with ONE ray per (side-leaf, endpoint).
        Counting: `2*n_prongs + 2*S = 2Q`, i.e. the COUNT half.  Hence every edge of the
        diagram is the lift of a traced prong or of a side-leaf, and the lift is a local
        isometry so the lengths are `ell_p` and `|sigma|` unchanged.
  (iii) ⇒⇒ THE LEMMA THIS PROBE TESTS, and it is the only step with no repo precedent.
        For a COMPLETELY PERIODIC direction on a compact half-translation surface,
              sum_{cylinders} circ_c  =  sum_{edges} ell_e.
        Proof: `X` is the union of the closed cylinders, glued along the singular leaf
        set, which is the union of the edges.  Each cylinder `C` has two boundary
        components, each a closed curve of length `circ_C` concatenated out of edges.
        Each edge has exactly TWO sides, and each side lies on exactly one boundary
        component of exactly one cylinder.  Counting total boundary length two ways,
              sum_C (|d+C| + |d-C|) = 2 sum_C circ_C  and  = sum_e 2 ell_e.   QED
        (A fold at a simple pole puts BOTH sides of one edge on the SAME boundary
        component; the count is unaffected, and `2|sigma|` on one cylinder is exactly what
        that looks like.)
⇒ (i)+(ii)+(iii) give [NCYL-340]'s length half, conditional on complete periodicity and on
nothing else.  ⚠ IN PARTICULAR IT DOES NOT USE (M1)/(*) ([NCYL-282], the open residual):
the boundary-length count does not care WHICH vertex an edge lands on, only that the edge
set is the prongs and the side-leaves.  The `a`/`b`/`d` SPLIT needs (*); this does not.

WHAT THIS PROBE ADDS, because the sum alone is one number a row and (iii) is where an
assembly could be wrong: it computes the INCIDENCE -- which cylinder sits on each side of
each edge -- and scores the lemma PER CYLINDER.  That is `C` equations a row instead of
`1`, and it tests the mechanism (each edge used exactly twice, distributed correctly)
rather than the conclusion.

HYPOTHESES (all able to fail)
  H1  ASSIGNMENT.  Offsetting a point of each edge perpendicular by `+eps` and `-eps` and
      tracing to closure gives, on both sides and at BOTH of two eps scales, a closed
      orbit matching exactly one store cylinder on `(circumference, k_per_strip)`.
  H2  ⇒⇒ THE LEMMA, PER CYLINDER.  For every cylinder `c`,
            sum over the edge-sides incident to `c` of `ell_e`  ==  2 * circ_c.
  H2b ⇒ THE SAME LEMMA ONE LEVEL FINER, and it is free once the `leaf` id exists: each
      cylinder has EXACTLY TWO boundary components, and EACH has length `circ_c` -- H2 is
      only their sum, so a wrong split between the two would pass H2 and fail this.
  H3  THE FOLD PREDICTION.  A side-leaf is a boundary of a SINGLE cylinder on both of its
      sides (the two sheets `u_side` and `u_side + 2Q` are `tau`-identified) -- so every
      side-leaf folds.  An interior prong may or may not.
  H4  NEGATIVE CONTROLS, and they are what stop H2 being bookkeeping that cannot fail.
      (C1) score H2 with each incident edge counted ONCE instead of per side;
      (C2) score H2 against a rotated (deranged) edge->cylinder assignment on the same row.
      Both must MISS on rows where H2 passes.
  ⚠ NO BRANCH IS NOMINATED as the one I expect to fail ([OPS-044]).

⚠ WHAT THIS IS NOT.  It does not touch (M1)/(*) ([NCYL-282]), it does not count SURFACE
cylinders on `Sigma` ([NCYL-338]/[NCYL-339] -- everything here is on `B`, which IS the
billiard's own class phase space and whose cylinder count [NCYL-280] already matched to
the store), and it does not discharge [NCYL-340] DO-NOT (b): the right-hand side still
includes `R` prongs and `R` is not a cone point.

PRIOR ART: `rulings.py` on [NCYL-340], [NCYL-280], [NCYL-281], [NCYL-282], [NCYL-292],
[NCYL-338], [NCYL-339]; `rulings.py --grep 'separatrix diagram'`, `'top boundary'`,
`'boundary length'`.
  -> [NCYL-280] owns `E = Q` and the ray degrees (`P`, `Q-P`, `1`x`Q`) -- PROVED.
  -> [NCYL-281] owns the side-leaves as edges and `b_bd = [eps==0] + [eps==Q]`, and the
     COUNT half's algebra -- PROVED; (*) itself MEASURED.
  -> [NCYL-292] owns `B`'s separatrix diagram as the necklace-of-kites `iota`-pairing --
     the finest description of the diagram in the repo, and it carries NO edge LENGTHS.
  -> `--grep 'top boundary'`/`'boundary length'`: no rulings.  The length bookkeeping of
     the diagram is not in the repo; the cylinder<->edge INCIDENCE is not either.

Usage:
  PYTHONPATH=.:engine:archive/scripts_2026-07:archive/scripts_2026-08:probes \
    .venv/bin/python3.13 probes/s482_boundary_incidence.py [run [QMAX] | report]
"""
import glob
import json
import math
import os
import sys
import time

import s437_oblique_cylinders as s437
import s439_exact_cells as s439
from right_triangle_billiards import RightTriangleBilliard

OUT = "data/s482_boundary_incidence.json"
STORE = "data/s439_exact_store"
CLOSE_TOL = 1e-7          # s437/s439's POS_TOL; closure of the offset orbit
EPS = (1e-6, 2e-7)        # the two offset scales; H1 requires they agree
ORBIT_CAP = 60_000        # store max wordlen is 3394; this is 17x headroom


def u_side(P, Q):
    """The direction index of the SIDE ITSELF ([NCYL-340]): `L2 -> 0`, `H -> P`,
    `L1 -> Q`, in units of `pi/(2Q)`."""
    return {"L1": Q, "L2": 0, "H": P}


def _pt(geo, side, w):
    _, _, base, tan = geo[side]
    return (base[0] + w * tan[0], base[1] + w * tan[1])


def _inward_normal(geo, side, centroid):
    """Unit normal of `side` pointing into the triangle."""
    _, _, base, tan = geo[side]
    n = (-tan[1], tan[0])
    mid = (base[0] + 0.5 * tan[0], base[1] + 0.5 * tan[1])
    if (centroid[0] - mid[0]) * n[0] + (centroid[1] - mid[1]) * n[1] < 0:
        n = (-n[0], -n[1])
    return n


def edge_list(P, Q, par, cap):
    """Every edge of `B`'s separatrix diagram in the class direction, as a geometric
    object: its length and, per SIDE of the edge, a launch (point, direction).

    Prongs come from `s439.prongs` / `s439.trace_prong` (the same generator [NCYL-280]'s
    diagram is traced from); side-leaves are the sides whose own direction has the class
    parity ([NCYL-340] / [NCYL-281]'s boundary rays)."""
    geo = s437.geometry(P, Q)
    verts = s439.vertices(P, Q)
    B = RightTriangleBilliard((P / Q) * (math.pi / 2))
    cen = ((verts["O"][0] + verts["R"][0] + verts["A"][0]) / 3.0,
           (verts["O"][1] + verts["R"][1] + verts["A"][1]) / 3.0)
    out = []

    for v0 in ("O", "R", "A"):
        for u0 in s439.prongs(v0, P, Q, par):
            bs, why = s439.trace_prong(B, geo, verts, Q, v0, u0, cap)
            pts = [verts[v0]] + [_pt(geo, b["side"], b["w"]) for b in bs]
            dirs = [u0] + [b["u_out"] for b in bs]
            segs = [(math.dist(pts[i], pts[i + 1]), i) for i in range(len(pts) - 1)]
            ell = sum(s[0] for s in segs)
            # sample the LONGEST segment: maximal clearance from both endpoints.
            _, i = max(segs)
            mid = (0.5 * (pts[i][0] + pts[i + 1][0]), 0.5 * (pts[i][1] + pts[i + 1][1]))
            th = dirs[i] * math.pi / (2 * Q)
            nrm = (-math.sin(th), math.cos(th))
            out.append({"kind": "prong", "v0": v0, "u0": u0, "why": why,
                        "n": len(bs), "ell": ell,
                        "sides": [{"pt": mid, "n": nrm, "sgn": +1, "u": dirs[i]},
                                  {"pt": mid, "n": nrm, "sgn": -1, "u": dirs[i]}]})

    us = u_side(P, Q)
    for s in ("L2", "H", "L1"):
        if us[s] % 2 != par:
            continue
        L, _, base, tan = geo[s]
        mid = (base[0] + 0.5 * L * tan[0], base[1] + 0.5 * L * tan[1])
        nrm = _inward_normal(geo, s, cen)
        # A side-leaf has both of its sides INSIDE the triangle: they are the two sheets
        # `u_side` and `u_side + 2Q`, not two sides of a line in the plane.
        out.append({"kind": "side", "side": s, "ell": L,
                    "sides": [{"pt": mid, "n": nrm, "sgn": +1, "u": us[s]},
                              {"pt": mid, "n": nrm, "sgn": +1, "u": (us[s] + 2 * Q) % (4 * Q)}]})
    return out, geo, B


def orbit_at(B, geo, Q, pt, u0, cap=ORBIT_CAP):
    """One full period of the class orbit through an interior point, as
    `(circumference, k_per_strip, leaf_id)` -- the first two are the store's own cylinder
    signature (`s439._attach_kc`: the tag is `side:u_out`).

    ⇒⇒ `leaf_id` IS WHAT RESOLVES THE BOUNDARY COMPONENTS, and it costs nothing: it is the
    least crossing key of the whole orbit, so two points on the SAME leaf get the same id.
    Offset points sit at perpendicular distance `eps` from their edge, i.e. at transverse
    coordinate `eps` from the boundary component they came from — so all the edge-sides of
    ONE boundary component lie on ONE leaf and share an id, and a cylinder's edge-sides
    split into exactly TWO ids.  Subset-summing the lengths would not have been unique."""
    px, py = pt
    th = u0 * math.pi / (2 * Q)
    vx, vy = math.cos(th), math.sin(th)
    cand = B._candidate_times(px, py, vx, vy)
    if not cand:
        return None
    t, side = min(cand, key=lambda z: z[0])
    px += t * vx
    py += t * vy
    vx, vy = B.reflect_velocity(side, vx, vy)
    u = s439._u_of(vx, vy, Q)
    if u is None:
        return None
    _, _, base, tan = geo[side]
    w = (px - base[0]) * tan[0] + (py - base[1]) * tan[1]
    s0 = (side, u, w)
    sig, circ, leaf = {}, 0.0, None
    for _ in range(cap):
        sig[f"{side}:{u}"] = sig.get(f"{side}:{u}", 0) + 1
        k = (side, u, round(w, 9))
        leaf = k if leaf is None or k < leaf else leaf
        cand = B._candidate_times(px, py, vx, vy)
        if not cand:
            return None
        t, nside = min(cand, key=lambda z: z[0])
        circ += t
        px += t * vx
        py += t * vy
        vx, vy = B.reflect_velocity(nside, vx, vy)
        nu = s439._u_of(vx, vy, Q)
        if nu is None:
            return None
        _, _, base, tan = geo[nside]
        nw = (px - base[0]) * tan[0] + (py - base[1]) * tan[1]
        side, u, w = nside, nu, nw
        if side == s0[0] and u == s0[1] and abs(w - s0[2]) < CLOSE_TOL:
            return circ, sig, leaf
    return None


def match(res, cyls):
    """Index of the unique store cylinder with this signature, or `None`."""
    if res is None:
        return None
    circ, sig, _leaf = res
    hit = [i for i, c in enumerate(cyls) if c["k_per_strip"] == sig]
    if len(hit) != 1:
        return None
    if abs(circ - cyls[hit[0]]["circumference"]) > 1e-6 * max(1.0, circ):
        return None
    return hit[0]


def do_row(P, Q, role, rec, cl):
    par = cl["parity"] % 2
    cyls = cl["cylinders"]
    t0 = time.time()
    edges, geo, B = edge_list(P, Q, par, rec["prong_cap"])

    inc = []            # (edge index, side index) -> cylinder index
    n_unres = 0
    for ei, e in enumerate(edges):
        for si, sd in enumerate(e["sides"]):
            got, leaf = [], None
            for eps in EPS:
                p = (sd["pt"][0] + sd["sgn"] * eps * sd["n"][0],
                     sd["pt"][1] + sd["sgn"] * eps * sd["n"][1])
                res = orbit_at(B, geo, Q, p, sd["u"])
                got.append(match(res, cyls))
                if eps == EPS[0] and res is not None:
                    leaf = "%s:%d:%.9f" % res[2]
            ci = got[0] if (got[0] is not None and got[0] == got[1]) else None
            if ci is None:
                n_unres += 1
            inc.append({"e": ei, "s": si, "c": ci, "leaf": leaf})

    # ---- H2: per cylinder, sum of incident edge-side lengths == 2 * circumference
    tot = [0.0] * len(cyls)
    once = [0.0] * len(cyls)          # C1: each incident edge counted once
    seen = [set() for _ in cyls]
    for r in inc:
        if r["c"] is None:
            continue
        tot[r["c"]] += edges[r["e"]]["ell"]
        if r["e"] not in seen[r["c"]]:
            seen[r["c"]].add(r["e"])
            once[r["c"]] += edges[r["e"]]["ell"]
    resid = [abs(tot[i] - 2 * c["circumference"]) / (2 * c["circumference"])
             for i, c in enumerate(cyls)]
    resid_c1 = [abs(once[i] - 2 * c["circumference"]) / (2 * c["circumference"])
                for i, c in enumerate(cyls)]

    # ---- H2b: the BOUNDARY COMPONENTS, grouped by leaf id.  Each cylinder must have
    # exactly TWO, and each must have total length == circ_c (H2 is only their SUM).
    comp = {}
    for r in inc:
        if r["c"] is None or r["leaf"] is None:
            continue
        comp.setdefault((r["c"], r["leaf"]), []).append(r["e"])
    per_cyl = {}
    for (ci, lf), es in comp.items():
        per_cyl.setdefault(ci, []).append((lf, es))
    n_two = sum(1 for ci in per_cyl if len(per_cyl[ci]) == 2)
    comp_resid = []
    for ci, gs in per_cyl.items():
        for lf, es in gs:
            tot = sum(edges[e]["ell"] for e in es)
            comp_resid.append(abs(tot - cyls[ci]["circumference"])
                              / cyls[ci]["circumference"])

    # ---- C2: rotate the assignment by one cylinder (a derangement when C > 1)
    tot2 = [0.0] * len(cyls)
    for r in inc:
        if r["c"] is None:
            continue
        tot2[(r["c"] + 1) % len(cyls)] += edges[r["e"]]["ell"]
    resid_c2 = [abs(tot2[i] - 2 * c["circumference"]) / (2 * c["circumference"])
                for i, c in enumerate(cyls)]

    # ---- H3: does each edge fold (both sides in the same cylinder)?
    folds = {"side": [0, 0], "prong": [0, 0]}     # [fold, no-fold]
    for ei, e in enumerate(edges):
        cs = [r["c"] for r in inc if r["e"] == ei]
        if any(c is None for c in cs):
            continue
        folds[e["kind"]][0 if cs[0] == cs[1] else 1] += 1

    ok = (n_unres == 0 and len(cyls) > 0)
    # ⚠ Score the controls PER CYLINDER, not per row ([OPS-041]): a row-level "the
    # control missed somewhere" hides how many of the `C` equations could not fail.
    return {"P": P, "Q": Q, "role": role, "parity": par, "C": len(cyls),
            "n_edges": len(edges), "n_unresolved": n_unres,
            "h1_ok": ok,
            "h2_worst": max(resid, default=None),
            "h2_ok": ok and max(resid, default=1.0) < 1e-9,
            "c1_worst": min(resid_c1, default=None),
            "c1_pass": sum(x < 1e-9 for x in resid_c1),
            "c2_worst": min(resid_c2, default=None),
            "c2_pass": sum(x < 1e-9 for x in resid_c2),
            "n_cyl_two_comps": n_two,
            "h2b_worst": max(comp_resid, default=None),
            "h2b_ok": ok and n_two == len(cyls) and max(comp_resid, default=1.0) < 1e-9,
            "components": {f"{ci}|{lf}": es for (ci, lf), es in comp.items()},
            "folds": folds,
            "sum_edges": sum(e["ell"] for e in edges),
            "sum_circ": sum(c["circumference"] for c in cyls),
            "incidence": [[r["e"], r["s"], r["c"], r["leaf"]] for r in inc],
            "edge_kinds": [e["kind"] for e in edges],
            "edge_ells": [e["ell"] for e in edges],
            "sec": round(time.time() - t0, 2)}


def run(qmax=10 ** 9):
    out = json.load(open(OUT)) if os.path.exists(OUT) else {}
    for f in sorted(glob.glob(f"{STORE}/*.json"), key=lambda p: (
            int(os.path.basename(p)[:-5].split("_")[1]),
            int(os.path.basename(p)[:-5].split("_")[0]))):
        rec = json.load(open(f))
        if rec["Q"] > qmax:
            continue
        for role, cl in rec["classes"].items():
            key = f"{rec['P']}/{rec['Q']}/{role}"
            if key in out or not cl.get("cylinders"):
                continue
            r = do_row(rec["P"], rec["Q"], role, rec, cl)
            out[key] = r
            print(f"{key:16s} C={r['C']:3d} E={r['n_edges']:3d} unres={r['n_unresolved']:2d} "
                  f"h2={r['h2_worst']:.2e} c1={r['c1_worst']:.2e} c2={r['c2_worst']:.2e} "
                  f"{r['sec']:6.1f}s", flush=True)
            json.dump(out, open(OUT, "w"))
    report()


def report():
    d = json.load(open(OUT))
    rows = list(d.values())
    n = len(rows)
    ok = [r for r in rows if r["h1_ok"]]
    print(f"\nROWS {n}   H1 (assignment, both eps agree) {len(ok)}/{n}")
    if not ok:
        return
    print(f"H2 (per cylinder, sum edge-sides == 2*circ)  {sum(r['h2_ok'] for r in ok)}/{len(ok)}"
          f"   worst rel {max(r['h2_worst'] for r in ok):.3e}")
    print(f"   cylinders scored {sum(r['C'] for r in ok)}   edges {sum(r['n_edges'] for r in ok)}")
    ncyl = sum(r["C"] for r in ok)
    print(f"C1 (each edge counted ONCE)    misses {ncyl - sum(r['c1_pass'] for r in ok)}"
          f"/{ncyl} cylinders   smallest miss {min(r['c1_worst'] for r in ok):.3e}")
    m = [r for r in ok if r["C"] > 1]
    ncyl2 = sum(r["C"] for r in m)
    print(f"C2 (assignment rotated by one) misses {ncyl2 - sum(r['c2_pass'] for r in m)}"
          f"/{ncyl2} cylinders   smallest miss "
          f"{min((r['c2_worst'] for r in m), default=float('nan')):.3e}")
    fs = [0, 0]
    fp = [0, 0]
    for r in ok:
        fs[0] += r["folds"]["side"][0]
        fs[1] += r["folds"]["side"][1]
        fp[0] += r["folds"]["prong"][0]
        fp[1] += r["folds"]["prong"][1]
    print(f"H2b boundary components: cylinders with exactly TWO "
          f"{sum(r['n_cyl_two_comps'] for r in ok)}/{ncyl}   each of length == circ_c, "
          f"worst rel {max((r['h2b_worst'] for r in ok if r['h2b_worst'] is not None), default=0):.3e}"
          f"   rows clean {sum(r['h2b_ok'] for r in ok)}/{len(ok)}")
    print(f"H3 side-leaves fold {fs[0]}/{fs[0]+fs[1]}    prongs fold {fp[0]}/{fp[0]+fp[1]}")
    bad = [f"{r['P']}/{r['Q']}/{r['role']}" for r in rows if not (r["h1_ok"] and r["h2_ok"])]
    if bad:
        print("NOT CLEAN:", " ".join(bad[:20]))
    qs = sorted({r["Q"] for r in ok})
    print(f"COVERAGE Q = {qs[0]}..{qs[-1]}, {len(ok)} classes, full-mode store rows")


def degrees(qmax=300):
    """⇒⇒ STEP (ii) OF THE PROOF, PER VERTEX AND WITH NO TRACING.  [NCYL-340]'s H1 is the
    SUM `n_prongs + S == Q`; the proof needs it vertex by vertex, because that is what
    says the diagram's ray-ends at `Z_O` / `Z_A` / the poles are exactly the prong slots
    (2 rays each) and the side-leaf endpoints (1 ray each):

        2*N_O + b_O == P        2*N_A + b_A == Q - P        2*N_R + b_R == Q

    with `N_V` the class prong slots at `V` (`s439.prongs`) and `b_V` the sides incident to
    `V` whose own direction has the class parity.  Degrees `P`/`Q-P`/`Q` are [NCYL-280]'s
    (`Z_O`, `Z_A`, and `Q` poles of one prong each).  Pure arithmetic: `40x`+ the reach of
    any traced check, and independent of the store."""
    at = {"O": ("L2", "H"), "A": ("H", "L1"), "R": ("L2", "L1")}
    tot = bad = 0
    worst = []
    for Q in range(3, qmax + 1):
        for P in range(1, Q):
            if math.gcd(P, Q) != 1:
                continue
            us = u_side(P, Q)
            for par in (0, 1):
                deg = {"O": P, "A": Q - P, "R": Q}
                for v in ("O", "A", "R"):
                    N = len(s439.prongs(v, P, Q, par))
                    b = sum(1 for s in at[v] if us[s] % 2 == par)
                    tot += 1
                    if 2 * N + b != deg[v]:
                        bad += 1
                        if len(worst) < 10:
                            worst.append((P, Q, par, v, 2 * N + b, deg[v]))
    print(f"PER-VERTEX RAY IDENTITY  {tot - bad}/{tot} vertex-rows, "
          f"coprime (P,Q) with Q <= {qmax}, both parities")
    if worst:
        print("  MISSES:", worst)


if __name__ == "__main__":
    a = sys.argv[1:]
    if a and a[0] == "report":
        report()
    elif a and a[0] == "degrees":
        degrees(int(a[1]) if len(a) > 1 else 300)
    else:
        run(int(a[1]) if len(a) > 1 else 10 ** 9)
