#!/usr/bin/env python3
"""
s446_covering_identity.py -- PRE-REGISTERED.  Queue item (a''): PROVE (M1).

WHAT THIS IS.  NOT a proof of (M1).  A DERIVATION, from [NCYL-280]'s already-proved
`E = Q`, that **(M1) is not a new lemma**: it is the CLOSED FORM of [NCYL-093]'s
COVERING LAW -- open since s320 -- and s444's own `n_nofold` counter, which it reported
as `0` and then dropped, IS that law's DEFECT term.  Plus the instrument that scores the
resulting exact identity on one pass of traces.

------------------------------------------------------------------------------------
THE DERIVATION.  Fix a centre `P/Q`, a class `eps`, and assume complete periodicity.
Work on `B = Sigma/tau`, [NCYL-280]'s genus-0 base, and let `psi` be the class
direction (a line direction, i.e. an element of `Z/2Q` of parity `eps`).

(1) ONE COPY PER DIRECTION.  `Sigma` is the `4Q`-copy unfolding, so its copies are
    indexed by the dihedral group `D` of order `4Q`, and each side type has `2Q` edges
    on `Sigma`.  A copy of side `sigma` carried by `g in D` has direction `g(s)` where
    `s = s(sigma)` is the base direction (`L2 -> 0`, `L1 -> Q`, `H -> P`, units
    `pi/(2Q)`); rotations give `s + 2m` and reflections `2b - s`, so EVERY copy has
    direction `= s (mod 2)` and each of the `Q` line-directions of that parity is
    realised by exactly `2` edges.  `tau` (rotation by `pi`) acts trivially on line
    directions and pairs those two, so **on `B` each side type has exactly ONE copy of
    each line direction of its own parity**, `Q` copies, `3Q` edges in all -- which is
    `B`'s edge count (`2Q` triangles).

(2) `Fix(iota)`.  `iota` is the involution of `B` induced by the reflection
    `delta_psi` of `Sigma` with axis direction `psi` (equivalently by
    `delta_{psi+Q}`: the two differ by `tau`).  A reflection of `Sigma` fixes no
    triangle copy, so its fixed set is a union of EDGES, namely the copies whose
    direction is the axis.  Hence

        Fix(iota) = { the copy of direction psi  of each sigma with s(sigma) = eps }
                  U { the copy of direction psi+Q of each sigma with s(sigma) = eps+Q }

    -- `h` HORIZONTAL arcs (leaves of the foliation) and `v` VERTICAL arcs
    (perpendicular walls), with

        h = #{sigma : s(sigma) = eps (mod 2)},   v = #{sigma : s(sigma) = eps+Q (mod 2)}.

    At odd `Q`, `h + v = 3` (one copy of each side, [NCYL-282]'s disk boundary); at even
    `Q`, `h = v` and a side either supplies BOTH arcs or neither.  A side's vertical copy
    is exactly the transversal of the perpendicular beam `perp sigma`, so `v` counts the
    sides with `e_sigma = eps` in [NCYL-112]'s notation.

(3) THE EDGE BOOKKEEPING.  `G` = the class's separatrix diagram on `B`, `E = Q` edges
    ([NCYL-280], PROVED).  `iota` permutes them.  An edge inside `Fix(iota)` must be a
    leaf, hence one of the `h` horizontal arcs; an `iota`-invariant edge NOT inside
    `Fix(iota)` is reversed by `iota` and so meets `Fix(iota)` in exactly one point,
    necessarily on a VERTICAL arc; and two leaves cannot cross, so every point of
    `G` on a vertical arc lies on an `iota`-invariant edge.  Therefore

        **Q = h + f + 2p**,   f = #(G  ^  vertical arcs),  p = # iota-swapped edge PAIRS.

(4) BOTH TERMS ARE THINGS ALREADY MEASURED.
    * `f`.  The points of `G` on side `sigma`'s vertical arc are exactly the interior
      CELL BOUNDARIES of the `perp sigma` beam, so `f = sum_{sigma: e_sigma = eps}
      (n_sigma - 1)` and

        **sum_{sigma : e_sigma = eps} n_sigma  =  Q - h + v - 2p.**

      With `p = 0` this is `Q-1` / `Q` / `Q+1` in the three parity classes -- [NCYL-093]'s
      table verbatim, and [NCYL-114]'s pair law `n_sigma + n_sigma' = Q+1` (`Q` odd) /
      `Q` (`Q` even) verbatim.  So the covering law is not an independent law: it is
      `E = Q` plus `p = 0`.
    * `p`.  An INTERIOR billiard prong is two prong-ends of `B`, swapped by `iota`
      ([NCYL-280]'s doubling).  So an `iota`-invariant non-horizontal edge is ONE billiard
      prong and an `iota`-swapped pair is TWO, giving `Q - h = f + 2p` -- i.e.
      `#interior prongs = Q - h` (H0 below) -- and a billiard prong is `iota`-invariant
      exactly when its orbit meets a perpendicular wall.  **Hence `2p` IS s444's
      `n_nofold`**, the count of prongs with no perpendicular bounce, and

        **(M1)  <=>  p = 0  <=>  the closed-form covering law.**

    (M1) is therefore the SAME open statement as [NCYL-093]'s covering half, which
    s320 measured on 29 rows `Q = 5..19` and s404 re-measured per cylinder
    ([NCYL-236]) -- 125 sessions of prior work on a residual s444 posed as new.
    ⚠ It also SHARPENS both: the covering law's failure is not qualitative, it is
    `2p` separatrices that avoid every perpendicular wall, and `p >= 0` makes the
    covering law an UPPER bound `2*C_total <= Q - h + v` before any measurement.

------------------------------------------------------------------------------------
HYPOTHESES (pre-registered, all able to fail).

  H0  CONTROL, RUN FIRST, PURE ARITHMETIC OVER THE WHOLE BOX.  `#interior prongs = Q - h`
      for every coprime `(P,Q)`, `Q <= 120`, both classes -- s439's inward windows on one
      side, this file's `h` on the other, sharing no code.  A mismatch means the
      `Fix(iota)` accounting in (2) is wrong and (3)(4) are void.
  H1  CONTROL.  For each side, `e_sigma = perp_index(sigma) mod 2` (s437's GEOMETRIC
      inward-normal index) equals `s(sigma) + Q mod 2` (this file's algebra), and the
      strip `(sigma, perp_index(sigma))` is present in `boundary_map` exactly for the
      `v` sides.  Shares no code with (2).
  H2  THE PREDICTION.  `sum_{sigma: e_sigma = eps} n_sigma == Q - h + v - n_nofold`,
      with `n_sigma` from `s439.boundary_map`'s cut list and `n_nofold` from the SAME
      traces.  ⚠ FREE CONTROL: `Q - h - f` must be EVEN and `>= 0` -- an off-by-one in
      any of the three terms breaks parity, and no tuning can restore it.
  H3  COVERAGE.  Re-read (M1) (`n_nofold == 0`) above s444's `Q <= 20` mechanism sweep.
      By H2 every new row is also a new row for the covering law, whose own sample stops
      at `Q <= 19`.

GUARDS.
  G1  A class with any prong not terminating at `cap2` is UNSCORABLE (`PRONG_CAPPED`) --
      s439's convention, [NCYL-020]'s truncation guard.  This is also the CP scope gate:
      `8/15` has a minimal component and is excluded there, which is why [NCYL-115]'s
      refutation of the three-side total is not a refutation of this identity.
  G2  One pass of traces feeds BOTH readings, so `f` and `n_nofold` cannot disagree
      about which prongs were traced.
  G3  Nothing here writes to any store.

PRIOR ART (grepped 'covering law', 'met_L1', 'n_H', 'Sigma_H', 'transversal', 'Fix(',
'perp_complete', 'palindrome', 'retrace', 'time reversal', 'n = Q-1', '8/15'; read
`rulings.md` [NCYL-093]/[NCYL-098]/[NCYL-112]/[NCYL-113]/[NCYL-114]/[NCYL-115]/
[NCYL-132]/[NCYL-133]/[NCYL-216]/[NCYL-217]/[NCYL-220]/[NCYL-236]/[NCYL-237]/
[NCYL-241]/[NCYL-271]/[NCYL-280]/[NCYL-281]/[NCYL-282], `foundations.md` sec 2/3/7/8 and
`cf_width_laws.md` sec s320/s323b) ->
  - [NCYL-093] IS the target: "all three classes are the one line
    `sum_{sigma : nu_sigma in the grid} n_sigma = 2*C_total` (= `Q-1`, `Q`, `Q+1`)",
    with "orphan in overlap" PROVED ([NCYL-132]) and the COVERING half MEASURED.  This
    probe does not re-measure that law; it DERIVES its closed form from `E = Q` and
    identifies its defect with s444's counter.
  - [NCYL-112] supplies `e_L1 = 0`, `e_L2 = Q`, `e_H = P+Q` -- rederived here as
    `s(sigma) + Q` and used as H1's control.
  - [NCYL-115] records `sum_sigma n_sigma = 2Q` REFUTED at `8/15`; under this identity
    the shortfall `8` is `2p` summed over the two classes, and `8/15` is out of the CP
    scope anyway (G1).
  - [NCYL-282] states (M1) and its disk, and calls it "the separatrix-level twin of
    `perp_complete`" -- it does not connect it to [NCYL-093].  s444's `mech_row` computes
    `n_nofold` and reports it only as a pass/fail.
  - [NCYL-220]/[NCYL-217] use the same reversal at the CELL level (a `J`-fold retraces
    via `R`); this is the SEPARATRIX-level statement and neither implies the other.
  - No probe in the repo computes `f`, `p`, or `h`/`v`, or scores `n_sigma` against a
    closed form derived from the base graph.

Run: PYTHONPATH=.:engine:archive/scripts_2026-07:archive/scripts_2026-08:probes \\
     .venv/bin/python3.13 probes/s446_covering_identity.py [h0|sweep] [--qmax=N]
"""
import json
import math
import sys
import time
from collections import Counter

from right_triangle_billiards import RightTriangleBilliard
import s437_oblique_cylinders as s437
import s439_exact_cells as s439
import s444_base_graph as s444

SIDES = ("L2", "L1", "H")
CAP = 200_000
CAP_BIG = 2_000_000


def base_dir(side, P, Q):
    """The side's OWN direction, units `pi/(2Q)`: `L2 -> 0`, `L1 -> Q`, `H -> P`."""
    return {"L2": 0, "L1": Q, "H": P}[side]


def hv_sides(P, Q, eps):
    """`(horizontal, vertical)` side lists: a side supplies a horizontal arc of
    `Fix(iota)` iff its own direction is in the class, a vertical arc iff its
    perpendicular is."""
    hor = [s for s in SIDES if base_dir(s, P, Q) % 2 == eps % 2]
    ver = [s for s in SIDES if (base_dir(s, P, Q) + Q) % 2 == eps % 2]
    return hor, ver


def n_prongs(P, Q, eps):
    return sum(len(s439.prongs(v, P, Q, eps)) for v in ("O", "R", "A"))


# ---------------------------------------------------------------- H0: pure arithmetic

def h0_control(qmax=120, verbose=True):
    """`#interior prongs == Q - h`, whole box.  s439's inward windows against this
    file's `Fix(iota)` count; no geometry, no tracing, nothing to tune."""
    bad, rows = [], 0
    for Q in range(3, qmax + 1):
        for P in range(1, Q):
            if math.gcd(P, Q) != 1:
                continue
            for eps in (0, 1):
                rows += 1
                hor, _ver = hv_sides(P, Q, eps)
                if n_prongs(P, Q, eps) != Q - len(hor):
                    bad.append((P, Q, eps, n_prongs(P, Q, eps), Q - len(hor)))
    if verbose:
        print(f"  H0  #interior prongs == Q - h : {rows - len(bad)}/{rows} class rows"
              f" (coprime, Q<={qmax})", flush=True)
        for r in bad[:6]:
            print("     ", r, flush=True)
    return bad, rows


def h1_control(qmax=40, verbose=True):
    """`e_sigma` from s437's geometric inward normal == `s(sigma)+Q` from the algebra."""
    bad, rows = [], 0
    for Q in range(3, qmax + 1):
        for P in range(1, Q):
            if math.gcd(P, Q) != 1:
                continue
            for side in SIDES:
                rows += 1
                if s437.perp_index(side, P, Q) % 2 != (base_dir(side, P, Q) + Q) % 2:
                    bad.append((P, Q, side))
    if verbose:
        print(f"  H1  e_sigma geometric == algebraic : {rows - len(bad)}/{rows}"
              f" (side rows, Q<={qmax})", flush=True)
    return bad, rows


# ---------------------------------------------------------------- the measurement

def row(P, Q, eps, cap=CAP, cap2=CAP_BIG):
    """ONE pass of prong traces; BOTH readings come off it.

    `n_nofold` = prongs with no perpendicular bounce ( = `2p` );
    `n_sigma`  = cells of the `perp sigma` beam, from the same bounces via s439's own
                 `_cut_list` (so the dedup convention is the enumerator's)."""
    t0 = time.time()
    geo = s437.geometry(P, Q)
    verts = s439.vertices(P, Q)
    B = RightTriangleBilliard((P / Q) * (math.pi / 2))
    hor, ver = hv_sides(P, Q, eps)
    want = {side: s437.perp_index(side, P, Q) for side in ver}

    bnds = {side: [] for side in ver}
    nofold, capped, nprong = 0, 0, 0
    for v0 in ("O", "R", "A"):
        for u0 in s439.prongs(v0, P, Q, eps):
            nprong += 1
            bounces, why = s439.trace_prong(B, geo, verts, Q, v0, u0, cap)
            if why == "cap" and cap2 > cap:
                bounces, why = s439.trace_prong(B, geo, verts, Q, v0, u0, cap2)
            if not why.startswith("vertex"):
                capped += 1
                continue
            if not any(bo["u_out"] == (bo["u_in"] + 2 * Q) % (4 * Q) for bo in bounces):
                nofold += 1
            for bo in bounces:                       # s439.boundary_map's own rule
                side = bo["side"]
                if side in want:
                    if bo["u_out"] == want[side]:
                        bnds[side].append(bo["w"])
                    if (bo["u_in"] + 2 * Q) % (4 * Q) == want[side]:
                        bnds[side].append(bo["w"])

    out = {"P": P, "Q": Q, "eps": eps, "lone": eps == P % 2, "n_prongs": nprong,
           "capped": capped, "h": len(hor), "v": len(ver), "hor": hor, "ver": ver,
           "secs": round(time.time() - t0, 2)}
    if capped:
        out["status"] = "PRONG_CAPPED"
        return out
    n_sigma = {side: len(s439._cut_list(geo, side, bnds[side])) - 1 for side in ver}
    f = sum(n - 1 for n in n_sigma.values())
    out.update({"status": "ok", "n_sigma": n_sigma, "f": f, "n_nofold": nofold,
                "sum_n": sum(n_sigma.values()),
                "pred_sum": Q - len(hor) + len(ver) - nofold,
                "defect": Q - len(hor) - f,
                "H0_ok": nprong == Q - len(hor)})
    out["H2_ok"] = out["sum_n"] == out["pred_sum"]
    out["parity_ok"] = out["defect"] % 2 == 0 and out["defect"] >= 0
    out["defect_eq_nofold"] = out["defect"] == nofold
    return out


def sweep(qmax=20, out_path="data/s446_covering.json", budget_s=None, verbose=True):
    t0 = time.time()
    tot = Counter()
    rows = []
    keys = [(P, Q, eps) for Q in range(3, qmax + 1) for P in range(1, Q)
            for eps in (0, 1) if math.gcd(P, Q) == 1]
    keys.sort(key=lambda k: (k[1], k[0], k[2]))
    for (P, Q, eps) in keys:
        if budget_s and time.time() - t0 > budget_s:
            print(f"  [budget stop at Q={Q}]", flush=True)
            break
        r = row(P, Q, eps)
        rows.append(r)
        tot["rows"] += 1
        if r["status"] != "ok":
            tot["capped_rows"] += 1
            if verbose:
                print(f"    PRONG_CAPPED {P}/{Q} eps={eps}", flush=True)
            continue
        tot["scored"] += 1
        tot["H0"] += r["H0_ok"]
        tot["H2"] += r["H2_ok"]
        tot["parity"] += r["parity_ok"]
        tot["defeq"] += r["defect_eq_nofold"]
        tot["M1"] += (r["n_nofold"] == 0)
        tot["prongs"] += r["n_prongs"]
        if not r["H2_ok"] or not r["parity_ok"]:
            print(f"    H2 MISS {P}/{Q} eps={eps}: sum_n={r['sum_n']} "
                  f"pred={r['pred_sum']} n_sigma={r['n_sigma']} "
                  f"defect={r['defect']} nofold={r['n_nofold']}", flush=True)
    if verbose:
        s = tot["scored"]
        print(f"\nSWEEP Q<={qmax}: {tot['rows']} class rows, {s} scored, "
              f"{tot['capped_rows']} PRONG_CAPPED, {tot['prongs']} prongs "
              f"[{round(time.time()-t0,1)}s]", flush=True)
        print(f"  H0 #prongs == Q-h                       : {tot['H0']}/{s}", flush=True)
        print(f"  H2 sum n_sigma == Q-h+v-n_nofold        : {tot['H2']}/{s}", flush=True)
        print(f"  -- free control: Q-h-f even and >= 0    : {tot['parity']}/{s}",
              flush=True)
        print(f"  -- defect Q-h-f == n_nofold             : {tot['defeq']}/{s}",
              flush=True)
        print(f"  M1 n_nofold == 0 (=> p = 0)             : {tot['M1']}/{s}", flush=True)
    json.dump({"summary": dict(tot), "qmax": qmax, "rows": rows},
              open(out_path, "w"), indent=1)
    print(f"wrote {out_path}", flush=True)
    return tot, rows


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "h0"
    kw = {}
    for a in sys.argv[2:]:
        if a.startswith("--qmax="):
            kw["qmax"] = int(a.split("=")[1])
        if a.startswith("--budget="):
            kw["budget_s"] = float(a.split("=")[1])
        if a.startswith("--out="):
            kw["out"] = a.split("=")[1]
    if mode == "h0":
        print("CONTROLS")
        h0_control(kw.get("qmax", 120))
        h1_control(min(kw.get("qmax", 120), 40))
        print("\nsix-case algebra (s444, for reference):")
        s444.six_case_check()
    elif mode == "sweep":
        out = kw.pop("out", "data/s446_covering.json")
        sweep(out_path=out, **kw)
    else:
        print(__doc__)


if __name__ == "__main__":
    main()
