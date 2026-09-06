#!/usr/bin/env python3
"""
s437_oblique_cylinders.py -- PRE-REGISTERED: enumerate the cylinders of BOTH
pi/(2Q) direction classes, from all three sides, and compare their lengths and
widths.  User-raised (s437): "I want to see how the lengths and widths of
cylinders change across the full set of pi/(2Q) angles.  Some of those will
overlap the cylinders for those perpendicular to the sides."

WHAT IS ACTUALLY NEW.  `probes/s315_cylinder_count.py` already sweeps all three
sides -- but its direction loop is `for k in range(2*Q)` at `theta = k*pi/Q`,
i.e. only the EVEN absolute direction index.  Every cylinder census, met-set,
covering law and `chi` column in this repo lives in that one class.  By
[NCYL-112] the direction index parity is CONSERVED, so `theta = u*pi/(2Q)` with
`u` ODD is a second, disjoint class -- a second direction on the surface, never
enumerated.  This probe is `s315_cylinder_count` with the direction grid refined
`2Q -> 4Q` and every result split by `u % 2`.

WHY ONLY TWO CLASSES (the structural claim H3 tests).  Reachable directions from
`u0` are `{u0 + even} u {-u0 + even}` = the parity class of `u0` ([NCYL-112]);
`u` and `u+2` differ by the deck rotation `pi/Q`.  So the 4Q launch angles carry
only TWO cylinder decompositions, and what varies with `u` inside a class is the
TRANSVERSAL, not the cylinders.  The three perpendicular strips have indices
`u_L1 = 2Q`, `u_L2 = Q`, `u_H = 3Q + P`, whose parities are `0`, `Q mod 2`,
`(P+Q) mod 2` -- [NCYL-112]'s `e_sigma` re-derived geometrically, and the reason
the classes always split 2+1 ([NCYL-113]).

HYPOTHESES (pre-registered).
  H1  CONTROL, RUN FIRST.  Restricted to EVEN `u`, this probe reproduces
      `s315_cylinder_count.run_row` on the shared rows (3,7),(3,8),(3,10),
      (5,12),(5,13): identical `n_cylinders_total`, `total_area_ratio` to 1e-9,
      identical `n_cylinders_meeting_sigma1`.  If this fails the refinement
      broke the enumerator and nothing below is readable.
  H2  Each parity class is INDEPENDENTLY a complete decomposition of the same
      surface: `total_area_ratio = 1.0` per class (Kac, against
      `area_exact = Q*cot(alpha)`), and `C_total <= g + s - 1` per class.
  H3  The class, not the angle, carries the cylinders: within one parity every
      `(side,u)` strip's slices canonicalise into ONE cylinder set, with
      `n_bad_slices = 0`.  Refutable -- if two strips of the same parity produced
      disjoint cylinder sets the two-class collapse would be wrong.
  H4  THE USER'S QUESTION.  Report per class the multiset of cylinder
      circumferences (LENGTHS) and heights `area/circumference` (WIDTHS), the
      counts `C_total`, and how the two classes differ.  No prediction is
      pre-registered for the comparison itself; it is the measurement.
  H5  THE OVERLAP.  Per cylinder, the FULL set of `(side,u)` strips that meet it
      -- generalising `meets_sigma1`, which `s315` records for `(L1, 2Q)` alone
      and `s320_perp_H` extends to the second PERPENDICULAR strip only.  Read
      off: (a) which cylinders are met by NO perpendicular strip at all, (b)
      whether the oblique strips of one side cover a class.

  CONTROLS, each able to fail independently.
  C1 = H1.  C2 = H2's per-class Kac total.  C3 = H2's topological bound.
  C4  Perpendicular-index parity: the computed `u_perp(sigma)` parities must
      equal [NCYL-112]'s `e_L1 = 0`, `e_L2 = Q mod 2`, `e_H = (P+Q) mod 2` on
      every row -- a free re-derivation of a PROVED statement.
  C5  Class membership of the perpendicular strips must split 2+1 on every row
      ([NCYL-113]).
  C6  Even-`P` rows are included (2/5, 2/7, 4/9): `s315`'s own ROWS are all odd
      `P`, and [NCYL-075] found even `P` behaves differently (no orphan, full
      sweep).  A law that only holds at odd `P` will show here.

  ⚠ NOT pre-registered as expected-to-fail: [OPS-044] bars the self-applied
  named-wrong-hypothesis device.  The scoring discipline used instead is
  [OPS-041] -- for each reported number, could it have come out otherwise?

PRIOR ART: grepped 'oblique', 'non-perpendicular', 'transversal', 'crossing
multiplicity', 'chi', 'R-cycle', 'even coset', 'odd coset', 'direction index',
'2Q directions', 'first_return(', 'outgoing_angle', 'theta_k', 'Sigma_k',
'covering law', 'lone side' across rulings.md / cf_width_laws.md /
computational_findings.md / future_directions.md / TOOLS.md / notation.md /
foundations.md ->
  - [NCYL-112] PROVES the conserved direction-index parity and gives `e_sigma`;
    [NCYL-113] the 2+1 split and the unique double normal.  This probe consumes
    both and re-derives them as controls C4/C5.  Neither was ever used to
    ENUMERATE the second class.
  - `probes/s315_cylinder_count.py` ([NCYL-060], TOOLS.md) is the enumerator
    being generalised.  Its four measured TRAPS are inherited verbatim and are
    NOT re-derived: (i) area by `height*circumference` overcounts by the crossing
    multiplicity (11.9x at 3/7) -- use Kac's next-segment `width*flux*seg0`;
    (ii) `total_area_ratio = 1.0` is blind to DUPLICATION -- pair it with the
    topological bound; (iii) do not merge on (word length + circumference) alone
    -- require word-content agreement under best cyclic rotation; (iv) NSAMP=400
    under-resolves at Q ~ 18-19.  `discover_cylinder_intervals`'s raw-word
    (not canon) boundary signature and its two-endpoint seg0 trapezoid are also
    load-bearing bug fixes and are carried unchanged.
  - `probes/s318_kac_missing.py` sweeps strips `Sigma_k = (L1, k*pi/Q)` at 8/15
    -- again `pi/Q`, so even class, and L1 only.
  - `probes/s320_perp_H.py`'s `perp_index(side,P,Q)` returns `None` when a side
    normal is "off the grid" ([NCYL-032]: L1 always, L2 iff Q even, H iff
    P = Q mod 2).  "Off the grid" is exactly "in the ODD class" -- so that
    None is this probe's subject, not an absence.
  - [NCYL-114] records even-`Q` `n_H` as OPEN with a resolution-UNSTABLE float
    value ("do not quote it; it is the lone side there").  At even `Q`,
    `e_H = 1`: the lone side IS the odd class.  Not re-measured blind here --
    it comes out of the class enumeration.
  - [NCYL-093]/[NCYL-236] own the covering laws, both stated for the two
    PERPENDICULAR transversals of one class.  H5 is their oblique extension.
  - `computational_findings.md` sections 43 and 45 are the ONLY prior oblique
    numerics.  Both are at `alpha = 1.0 RADIANS` (irrational), where the
    direction set is infinite and no cylinder decomposition exists; section 43's
    "the involution holds for ALL directions pi + n*alpha" is tested there at
    `alpha = pi/5` and `pi/7` only (both `P = 2`).  Neither is evidence at a
    rational centre and neither is reused.
  - `future_directions.md` C1 "non-perpendicular launch" is the roster item;
    it has never been run.
  - No hit anywhere for a cylinder census at an odd direction index.

Run: PYTHONPATH=.:engine:archive/scripts_2026-07:archive/scripts_2026-08:probes \
     .venv/bin/python3.13 probes/s437_oblique_cylinders.py
"""
import hashlib
import json
import math
import sys
import time

from right_triangle_billiards import RightTriangleBilliard
import s315_oddP_census as census
import s315_cylinder_count as s315

POS_TOL = 1e-7
CAP = 20000
NSAMP_COARSE = 1200  # ⚠ NOT s315's 400.  s315's trap (iv) says 400 under-resolves at Q ~ 18-19;
                     # on the 4Q grid it under-resolves much earlier.  MEASURED this session
                     # (data/s437_convergence.json): at 400, three of eight rows carry a
                     # SPURIOUS extra cylinder -- 3/8 odd (C 5 -> 4, area_ratio 0.999870 -> 1.0),
                     # 4/9 even (5 -> 4, 0.9999995 -> 1.0), 3/10 odd (6 -> 5) -- and each spurious
                     # cylinder presents as "met by no perpendicular strip", i.e. exactly the
                     # headline H5 would report.  All three are stable at 1200 and again at 3000.
                     # ⚠ A missed cell BOUNDARY leaves `sum_width = |side|` intact, so the
                     # per-strip completeness check CANNOT see this failure (it read 0/402
                     # incomplete strips while three rows were wrong); the area ratio and the
                     # count are what move.  Override with argv[2].
REFINE_STEPS = 40

ROWS = [
    (3, 7), (2, 7), (3, 8), (2, 5), (4, 9), (3, 10), (5, 12), (5, 13),
]
H1_ROWS = [(3, 7), (3, 8), (3, 10), (5, 12), (5, 13)]   # rows s315 also runs


# ---------------------------------------------------------------- geometry
def geometry(P, Q):
    """Verbatim `s315_cylinder_count.geometry` -- not re-derived."""
    return s315.geometry(P, Q)


def perp_index(side, P, Q):
    """Direction index `u` (unit pi/(2Q), mod 4Q) of `side`'s INWARD normal.
    Geometric, not table-driven: L1 -> pi, L2 -> pi/2, H -> 3pi/2 + alpha."""
    geo = geometry(P, Q)
    nx, ny = geo[side][1]
    ang = math.atan2(ny, nx) % (2.0 * math.pi)
    uf = ang * 2.0 * Q / math.pi
    u = int(round(uf))
    assert abs(uf - u) < 1e-6 * Q, f"{side} normal off the 4Q grid: {uf}"
    return u % (4 * Q)


def flux_weight(side, u, Q, geo):
    th = u * math.pi / (2 * Q)
    _, nrm, _, _ = geo[side]
    return math.cos(th) * nrm[0] + math.sin(th) * nrm[1]


def trace_full_period(B, geo, Q, side0, x0, u0, cap):
    """`s315_cylinder_count.trace_full_period` with the direction grid refined
    2Q -> 4Q.  Returns (word, circumference, seg0, closed)."""
    _, _, base0, tan0 = geo[side0]
    px = base0[0] + x0 * tan0[0]
    py = base0[1] + x0 * tan0[1]
    th = u0 * math.pi / (2 * Q)
    vx, vy = math.cos(th), math.sin(th)

    word = []
    circ = 0.0
    seg0 = None
    for _ in range(1, cap + 1):
        cand = B._candidate_times(px, py, vx, vy)
        if not cand:
            return None, None, seg0, False
        t, side = min(cand, key=lambda z: z[0])
        if seg0 is None:
            seg0 = t
        circ += t
        px += t * vx
        py += t * vy
        word.append(side)
        vx, vy = B.reflect_velocity(side, vx, vy)
        uf = (math.atan2(vy, vx) % (2.0 * math.pi)) * 2.0 * Q / math.pi
        if abs(uf - round(uf)) > 1e-6 * 2 * Q:
            return None, None, seg0, False
        u = int(round(uf)) % (4 * Q)
        _, _, base, tan = geo[side]
        w = (px - base[0]) * tan[0] + (py - base[1]) * tan[1]
        if side == side0 and u == u0 and abs(w - x0) < POS_TOL:
            return tuple(word), circ, seg0, True
    return None, None, seg0, False


def discover_cylinder_intervals(B, geo, Q, side, u, cap):
    """`s315_cylinder_count.discover_cylinder_intervals`, 4Q grid.  The two
    load-bearing details are carried unchanged: the boundary signature is the
    RAW word (canonicalising here masks real boundaries), and `seg0` is the
    two-endpoint trapezoid (exact for the affine integrand; a midpoint sample
    overcounted 3/7's area by 8%)."""
    ell = geo[side][0]
    xs = [ell * (i + 0.5) / NSAMP_COARSE for i in range(NSAMP_COARSE)]
    sigs = []
    for x in xs:
        w, c, s0, ok = trace_full_period(B, geo, Q, side, x, u, cap)
        sigs.append(w if ok else None)

    boundaries = [0.0]
    for i in range(NSAMP_COARSE - 1):
        if sigs[i] != sigs[i + 1]:
            lo, hi = xs[i], xs[i + 1]
            sig_lo = sigs[i]
            for _ in range(REFINE_STEPS):
                mid = 0.5 * (lo + hi)
                w, c, s0, ok = trace_full_period(B, geo, Q, side, mid, u, cap)
                if (w if ok else None) == sig_lo:
                    lo = mid
                else:
                    hi = mid
            b = 0.5 * (lo + hi)
            if b - boundaries[-1] > 1e-9:
                boundaries.append(b)
    boundaries.append(ell)

    out = []
    for lo, hi in zip(boundaries[:-1], boundaries[1:]):
        width = hi - lo
        if width < 1e-10:
            continue
        mid = 0.5 * (lo + hi)
        w, c, s0_mid, ok = trace_full_period(B, geo, Q, side, mid, u, cap)
        eps = width * 1e-6
        _, _, s0_lo, ok_lo = trace_full_period(B, geo, Q, side, lo + eps, u, cap)
        _, _, s0_hi, ok_hi = trace_full_period(B, geo, Q, side, hi - eps, u, cap)
        seg0_avg = (s0_lo + s0_hi) / 2.0 if (ok_lo and ok_hi) else s0_mid
        canon = s315.canon(w) if ok else None
        out.append({"side": side, "u": u, "lo": lo, "hi": hi, "width": width,
                    "word": w, "canon": canon,
                    # stable short join key for the CSV export: the canonical word is the
                    # cylinder's identity, but it runs to 10^2..10^4 letters, so cells and
                    # cylinders are joined on its digest instead of on the word itself.
                    "key": (hashlib.sha1(repr(canon).encode()).hexdigest()[:10]
                            if canon is not None else None),
                    "circumference": c, "seg0": seg0_avg, "ok": ok})
    return out


# ---------------------------------------------------------------- per-class assembly
def assemble(slices, P, Q, area_exact, perp_u):
    """Group one parity class's slices into cylinders.  Kac area only
    (`width*flux*seg0`) -- `height*circumference` overcounts by the crossing
    multiplicity (s315 trap (i))."""
    groups = {}
    for s in slices:
        groups.setdefault(s["canon"], []).append(s)

    cylinders = []
    for cw, grp in groups.items():
        area = sum(s["area_contrib"] for s in grp)
        circs = [s["circumference"] for s in grp]
        circ_lo, circ_hi = min(circs), max(circs)
        circ_mid = 0.5 * (circ_lo + circ_hi)
        met = sorted({(s["side"], s["u"]) for s in grp})
        keys = sorted({s["key"] for s in grp})
        met_perp = sorted(sd for sd, uu in perp_u.items() if (sd, uu) in met)
        cylinders.append({
            "canon_word": cw,
            "keys": keys,
            "wordlen": len(cw),
            "n_slices": len(grp),
            "circumference_min": circ_lo, "circumference_max": circ_hi,
            "circumference_consistent": (circ_hi - circ_lo) < 1e-6 * max(1.0, circ_hi),
            "area": area,
            "height_derived": area / circ_mid if circ_mid else None,
            "sum_width_flux": sum(s["width"] * s["flux"] for s in grp),
            "met_strips": [list(m) for m in met],
            "n_met_strips": len(met),
            "met_perp_sides": met_perp,
            "meets_sigma1": ("L1", perp_u["L1"]) in met,
        })

    n_before = len(cylinders)
    # ⚠ `merge_duplicate_cylinders` rebuilds a merged entry from a FIXED key set and drops
    # every field it does not know about -- the same trap s320_perp_H records for
    # `meets_sigma1`.  A merged cylinder would therefore come back with NO met-strip set and
    # read as "met by no perpendicular", which is exactly the quantity H5 reports.  s315's
    # calibrated merge logic stays authoritative; the union is re-applied afterwards, matching
    # each merged entry back to its group through the `merged_areas` it carries.
    by_area = {c["area"]: c for c in cylinders}
    cylinders = s315.merge_duplicate_cylinders(cylinders)
    for c in cylinders:
        if "merged_from_n_groups" not in c:
            continue
        grp = [by_area[a] for a in c["merged_areas"] if a in by_area]
        met = sorted({tuple(m) for g in grp for m in g["met_strips"]})
        c["met_strips"] = [list(m) for m in met]
        c["n_met_strips"] = len(met)
        c["met_perp_sides"] = sorted({sd for sd, uu in perp_u.items() if (sd, uu) in met})
        c["keys"] = sorted({k for g in grp for k in g["keys"]})
        c["wordlen"] = len(c["canon_word"][0]) if c["canon_word"] else None
        c["sum_width_flux"] = sum(g["sum_width_flux"] for g in grp)
        c["merge_recovered_groups"] = len(grp) == c["merged_from_n_groups"]
        cm = 0.5 * (c["circumference_min"] + c["circumference_max"])
        c["height_derived"] = c["area"] / cm if cm else None   # merge() leaves this None
    total_area = sum(c["area"] for c in cylinders)
    bound, g, s = s315.topological_bound(P, Q)
    return {
        "C_total": len(cylinders),
        "C_before_merge": n_before,
        "total_area_ratio": total_area / area_exact,
        "topological_bound": bound, "genus": g, "n_singularities": s,
        "topological_bound_ok": len(cylinders) <= bound,
        "circumferences": sorted(round(0.5 * (c["circumference_min"]
                                              + c["circumference_max"]), 9)
                                 for c in cylinders),
        "heights": sorted(round(c["height_derived"], 9) for c in cylinders
                          if c["height_derived"] is not None),
        "cylinders": cylinders,
    }


def run_row(P, Q, cap=CAP):
    alpha = (P / Q) * (math.pi / 2)
    B = RightTriangleBilliard(alpha)
    geo = geometry(P, Q)
    area_exact = Q / math.tan(alpha)

    perp_u = {sd: perp_index(sd, P, Q) for sd in ("L1", "L2", "H")}
    e_pred = {"L1": 0, "L2": Q % 2, "H": (P + Q) % 2}          # C4: [NCYL-112]
    c4_ok = all(perp_u[sd] % 2 == e_pred[sd] for sd in perp_u)
    parities = [perp_u[sd] % 2 for sd in ("L1", "L2", "H")]
    c5_ok = sorted(parities.count(0) for _ in [0])[0] in (1, 2) and \
        {parities.count(0), parities.count(1)} == {1, 2}       # C5: 2+1 split

    by_parity = {0: [], 1: []}
    bad = {0: 0, 1: 0}
    strips = []
    t0 = time.time()
    for side in geo:
        for u in range(4 * Q):
            fl = flux_weight(side, u, Q, geo)
            if fl <= 1e-12:
                continue
            ivs = discover_cylinder_intervals(B, geo, Q, side, u, cap)
            good = []
            for iv in ivs:
                iv["flux"] = fl
                iv["area_contrib"] = iv["width"] * fl * iv["seg0"]
                if iv["ok"]:
                    by_parity[u % 2].append(iv)
                    good.append(iv)
                else:
                    bad[u % 2] += 1
            # PER-STRIP VIEW -- the beam launched from `side` at angle u*pi/(2Q).  `n_cells` is
            # the generalisation of [NCYL-093]/[NCYL-114]'s `n_sigma` off the perpendicular;
            # `sum_width_flux` is the flux identity `sum_c k_c h_c = |side| * sin(phi)`, which
            # is a self-check on the strip and equals `|side| * flux` if the strip is complete.
            strips.append({
                "side": side, "u": u, "parity": u % 2,
                "is_perp": u == perp_u[side],
                "flux": fl,
                "n_cells": len(good),
                "n_bad": len(ivs) - len(good),
                "n_distinct_cylinders": len({iv["canon"] for iv in good}),
                # ⚠ PER-CELL rows, in POSITION order along the side.  An earlier draft kept
                # two INDEPENDENTLY SORTED lists (`cell_widths`, `wordlens`); zipping those
                # pairs the i-th narrowest cell with the i-th shortest word, which are not
                # the same cell.  `key` joins a cell to its cylinder in the class table.
                "cells": [{"lo": round(iv["lo"], 12), "hi": round(iv["hi"], 12),
                           "width": iv["width"], "wordlen": len(iv["word"]),
                           "key": iv["key"]} for iv in good],
                "cell_widths": sorted(round(iv["width"], 9) for iv in good),
                "wordlens": sorted(len(iv["word"]) for iv in good),
                "sum_width": sum(iv["width"] for iv in good),
                "side_length": geo[side][0],
                "sum_width_flux": sum(iv["width"] * fl for iv in good),
            })

    out = {"P": P, "Q": Q, "area_exact": area_exact, "alpha": alpha,
           "perp_u": perp_u, "perp_parity": {k: v % 2 for k, v in perp_u.items()},
           "e_pred_NCYL112": e_pred, "C4_perp_parity_ok": c4_ok,
           "C5_split_2_1_ok": c5_ok,
           "n_bad_slices": bad, "seconds": None, "strips": strips, "classes": {}}
    for par in (0, 1):
        res = assemble(by_parity[par], P, Q, area_exact, perp_u)
        res["n_bad_slices"] = bad[par]
        res["n_slices"] = len(by_parity[par])
        res["perp_sides_in_class"] = [sd for sd in ("L1", "L2", "H")
                                      if perp_u[sd] % 2 == par]
        res["n_cyl_met_by_no_perp"] = sum(1 for c in res["cylinders"]
                                          if not c.get("met_perp_sides"))
        out["classes"][str(par)] = res
    out["seconds"] = round(time.time() - t0, 1)
    return out


# ---------------------------------------------------------------- H1 control
def h1_control(rows):
    """Restrict to EVEN u and compare against s315_cylinder_count.run_row."""
    report = []
    for (P, Q) in rows:
        mine = run_row(P, Q)
        theirs = s315.run_row(P, Q)
        c0 = mine["classes"]["0"]
        report.append({
            "P": P, "Q": Q,
            "C_total_mine": c0["C_total"], "C_total_s315": theirs["n_cylinders_total"],
            "C_match": c0["C_total"] == theirs["n_cylinders_total"],
            "area_mine": c0["total_area_ratio"], "area_s315": theirs["total_area_ratio"],
            "area_match": abs(c0["total_area_ratio"] - theirs["total_area_ratio"]) < 1e-9,
            "met_mine": sum(1 for c in c0["cylinders"] if c.get("meets_sigma1")),
            "met_s315": theirs["n_cylinders_meeting_sigma1"],
        })
        r = report[-1]
        r["met_match"] = r["met_mine"] == r["met_s315"]
        print(f"  H1 {P}/{Q}: C {r['C_total_mine']} vs {r['C_total_s315']} "
              f"{'OK' if r['C_match'] else 'MISMATCH'} | area "
              f"{r['area_mine']:.9f} vs {r['area_s315']:.9f} "
              f"{'OK' if r['area_match'] else 'MISMATCH'} | met "
              f"{r['met_mine']} vs {r['met_s315']} "
              f"{'OK' if r['met_match'] else 'MISMATCH'}")
        sys.stdout.flush()
    return report


# ---------------------------------------------------------------- CSV export
def _fmt(x, nd=9):
    return "" if x is None else (f"{x:.{nd}g}" if isinstance(x, float) else str(x))


def export_csv(src="data/s437_oblique_cylinders.json", stem="data/s437"):
    """Three flat tables off the run JSON.  Columns are derived from the records
    ([OPS-073]: no hand-maintained parallel column list), except for the few
    joined/derived keys named explicitly below.  The canonical WORD is deliberately
    NOT exported -- it runs to 146+ letters and a merged cylinder carries a tuple of
    them; `wordlen` is the usable handle."""
    import csv
    d = json.load(open(src))
    nsamp = d.get("nsamp_coarse")

    # (1) one row per CLASS -- the headline table
    cls_rows = []
    for r in d["rows"]:
        P, Q = r["P"], r["Q"]
        for par in ("0", "1"):
            c = r["classes"][par]
            perp = c["perp_sides_in_class"]
            lone = len(perp) == 1
            L, H = c["circumferences"], c["heights"]
            cls_rows.append({
                "P": P, "Q": Q, "P_parity": P % 2, "Q_parity": Q % 2,
                "class_parity": int(par), "class_role": "lone" if lone else "paired",
                "perp_sides_in_class": "|".join(perp),
                "C_total": c["C_total"], "C_before_merge": c["C_before_merge"],
                "C_pred": (Q // 2) if lone else -(-Q // 2),
                "total_area_ratio": _fmt(c["total_area_ratio"], 12),
                "area_exact": _fmt(r["area_exact"]),
                "topological_bound": c["topological_bound"],
                "topological_bound_ok": c["topological_bound_ok"],
                "genus": c["genus"], "n_singularities": c["n_singularities"],
                "n_slices": c["n_slices"], "n_bad_slices": c["n_bad_slices"],
                "n_cyl_met_by_no_perp": c["n_cyl_met_by_no_perp"],
                "sum_circumference": _fmt(sum(L)), "max_circumference": _fmt(max(L)),
                "min_circumference": _fmt(min(L)),
                "min_height": _fmt(min(H)) if H else "", "max_height": _fmt(max(H)) if H else "",
                "nsamp_coarse": nsamp,
            })

    # (2) one row per CYLINDER
    cyl_rows = []
    for r in d["rows"]:
        P, Q = r["P"], r["Q"]
        for par in ("0", "1"):
            c = r["classes"][par]
            perp = c["perp_sides_in_class"]
            for i, cy in enumerate(sorted(c["cylinders"],
                                          key=lambda z: z["circumference_min"])):
                circ = 0.5 * (cy["circumference_min"] + cy["circumference_max"])
                cyl_rows.append({
                    "P": P, "Q": Q, "class_parity": int(par),
                    "class_role": "lone" if len(perp) == 1 else "paired",
                    "cyl_rank_by_circ": i,
                    "circumference": _fmt(circ),
                    "height": _fmt(cy["height_derived"]),
                    "area": _fmt(cy["area"]),
                    "modulus": _fmt(cy["height_derived"] / circ if circ else None),
                    "wordlen": cy.get("wordlen"),
                    "n_slices": cy["n_slices"],
                    "n_met_strips": cy.get("n_met_strips"),
                    "met_perp_sides": "|".join(cy.get("met_perp_sides") or []),
                    "meets_sigma1": cy.get("meets_sigma1"),
                    "met_strips": "|".join(f"{s}:{u}" for s, u in cy.get("met_strips", [])),
                    "cyl_keys": "|".join(cy.get("keys") or []),
                    "circumference_consistent": cy["circumference_consistent"],
                    "merged_from_n_groups": cy.get("merged_from_n_groups", 1),
                    "nsamp_coarse": nsamp,
                })

    # (3) one row per STRIP (side, u) -- n(sigma,u) as a function of launch angle
    strip_rows = []
    for r in d["rows"]:
        P, Q = r["P"], r["Q"]
        for s in sorted(r["strips"], key=lambda z: (z["side"], z["u"])):
            wl = s["wordlens"]
            strip_rows.append({
                "P": P, "Q": Q, "side": s["side"], "u": s["u"],
                "j_from_side": (s["u"] - Q) % (4 * Q) if s["side"] == "L1" else "",
                "angle_deg": _fmt(180.0 * s["u"] / (2 * Q)),
                "parity": s["parity"], "is_perp": s["is_perp"],
                "class_role": ("lone" if len(r["classes"][str(s["parity"])]
                                          ["perp_sides_in_class"]) == 1 else "paired"),
                "flux": _fmt(s["flux"]), "n_cells": s["n_cells"], "n_bad": s["n_bad"],
                "n_distinct_cylinders": s["n_distinct_cylinders"],
                "side_length": _fmt(s["side_length"]),
                "sum_width": _fmt(s["sum_width"]),
                "sum_width_over_side": _fmt(s["sum_width"] / s["side_length"], 12),
                "min_cell_width": _fmt(min(s["cell_widths"])) if s["cell_widths"] else "",
                "max_cell_width": _fmt(max(s["cell_widths"])) if s["cell_widths"] else "",
                "min_wordlen": min(wl) if wl else "", "max_wordlen": max(wl) if wl else "",
                "nsamp_coarse": nsamp,
            })

    # (4) one row per CELL of constant word -- the full partition of every strip, in
    # POSITION order along the side.  Join to table (2) on `cyl_key` in `cyl_keys`.
    cell_rows = []
    for r in d["rows"]:
        P, Q = r["P"], r["Q"]
        for s in sorted(r["strips"], key=lambda z: (z["side"], z["u"])):
            for i, cell in enumerate(s["cells"]):
                cell_rows.append({"P": P, "Q": Q, "side": s["side"], "u": s["u"],
                                  "parity": s["parity"], "is_perp": s["is_perp"],
                                  "cell_index": i,
                                  "lo": _fmt(cell["lo"], 12), "hi": _fmt(cell["hi"], 12),
                                  "width": _fmt(cell["width"]),
                                  "wordlen": cell["wordlen"], "cyl_key": cell["key"],
                                  "nsamp_coarse": nsamp})

    for name, rows in (("classes", cls_rows), ("cylinders", cyl_rows),
                       ("strips", strip_rows), ("cells", cell_rows)):
        path = f"{stem}_{name}.csv"
        with open(path, "w", newline="") as f:
            wtr = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            wtr.writeheader()
            wtr.writerows(rows)
        print(f"  wrote {path}  ({len(rows)} rows x {len(rows[0])} cols)")
    return {"classes": len(cls_rows), "cylinders": len(cyl_rows),
            "strips": len(strip_rows), "cells": len(cell_rows)}


def main():
    global NSAMP_COARSE
    only = sys.argv[1] if len(sys.argv) > 1 else None
    if len(sys.argv) > 2:
        NSAMP_COARSE = int(sys.argv[2])
    out = {"nsamp_coarse": NSAMP_COARSE}
    print(f"NSAMP_COARSE = {NSAMP_COARSE}")

    if only in (None, "h1"):
        print("H1 CONTROL -- even-u restriction vs s315_cylinder_count")
        out["H1"] = h1_control(H1_ROWS)

    if only == "csv":
        print("CSV EXPORT")
        out["csv"] = export_csv()
        return

    if only in (None, "rows"):
        print("\nFULL 4Q SWEEP")
        out["rows"] = []
        for (P, Q) in ROWS:
            r = run_row(P, Q)
            out["rows"].append(r)
            print(f"\n{P}/{Q}  (alpha={r['alpha']:.6f})  {r['seconds']}s   "
                  f"perp u={r['perp_u']} parity={r['perp_parity']}  "
                  f"C4={'OK' if r['C4_perp_parity_ok'] else 'FAIL'} "
                  f"C5={'OK' if r['C5_split_2_1_ok'] else 'FAIL'}")
            for par in ("0", "1"):
                c = r["classes"][par]
                print(f"   class u%2={par}: C_total={c['C_total']} "
                      f"(pre-merge {c['C_before_merge']})  "
                      f"area_ratio={c['total_area_ratio']:.9f}  "
                      f"bound {c['C_total']}<={c['topological_bound']} "
                      f"{'OK' if c['topological_bound_ok'] else 'FAIL'}  "
                      f"bad={c['n_bad_slices']}  perp sides in class "
                      f"{c['perp_sides_in_class']}  "
                      f"cyl met by no perp = {c['n_cyl_met_by_no_perp']}")
                print(f"      lengths {[round(x,5) for x in c['circumferences']]}")
                print(f"      widths  {[round(x,7) for x in c['heights']]}")
            for side in ("L1", "L2", "H"):
                for par in (0, 1):
                    row = [s for s in r["strips"] if s["side"] == side and s["parity"] == par]
                    row.sort(key=lambda s: s["u"])
                    tag = "".join("*" if s["is_perp"] else " " for s in row)
                    print(f"      n(cells) {side} par{par}: "
                          + " ".join(f"{s['u']}:{s['n_cells']}"
                                     + ("*" if s["is_perp"] else "") for s in row))
            sys.stdout.flush()

    with open("data/s437_oblique_cylinders.json", "w") as f:
        json.dump(out, f, indent=1)
    print("\nwrote data/s437_oblique_cylinders.json")
    if only in (None, "rows"):
        print("CSV EXPORT")
        export_csv()


if __name__ == "__main__":
    main()
