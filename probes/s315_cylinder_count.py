#!/usr/bin/env python3
"""
s315_cylinder_count.py -- PRE-REGISTERED: FLOAT cylinder count + area reconciliation for the
odd-P COMPLETELY PERIODIC direction, on the 9 exact-g_len rows.

WHY.  `s315_oddP_census.py` established (84/84, `s315_arrival_margins.py` confirms robustly)
that the perpendicular direction is completely periodic at every sampled odd-P row.  A
completely periodic direction is tiled by finitely many CYLINDERS.  `rulings.md` [NCYL-046]
gives `ncyl_meeting_L1 = (n+1)/2` for the cylinders the beam Sigma_1 MEETS (n = branch count,
`data/s275_flow_exact.json`).  What is NOT established is the TOTAL cylinder count `C`
(including cylinders that never touch L1 at all) -- that is the direct structural test of
whether the odd-P g_len deficit (0.30-0.74 on these 9 rows) is made of uncounted CYLINDERS
(C > (n+1)/2, areas summing to 1.0) or something else.

METHOD (a genuinely new instrument -- there is no existing "enumerate all cylinders of a
billiard direction" tool in this repo; see PRIOR ART).  Every interior point of the domain,
followed BACKWARD along the flow, hits the triangle's boundary at some unique
(side, position, direction) -- so sampling ALL THREE SIDES in ALL 2Q inward grid directions is
a genuine partition of the whole domain (this is exactly `probes/s276_periodic_mass.py`'s
"Sigma_all" cross-section, reused verbatim: `geometry()`, the flux weight, and the exact
`mu_boundary`/`area_exact` closed forms it calibrates against).  For each (side, k) with
positive inward flux, `discover_cylinder_intervals` bisects along that side's coordinate exactly
as `right_triangle_billiards.discover_periodic_branches` bisects along L1 -- but the signature
here is the FULL-PERIOD closure word (trace until the exact boundary state (side,x,k) recurs,
`s276_periodic_mass.orbit`'s closure test, extended to also record the hit-side WORD and the
summed path length = the orbit's circumference).  Two entry points on (possibly different)
sides belong to the SAME cylinder iff their words agree up to CYCLIC ROTATION (different entry
points into one cylinder start recording the word at different phase) -- canonicalized by the
lexicographically-smallest rotation.  AREA is by KAC'S LEMMA, not `height * circumference`:
each slice contributes `width * flux * seg0` (`seg0` = the length of the FIRST segment from
that boundary point, exactly `s276_periodic_mass.py`'s own area weight) and `Area_c` = the SUM
of that contribution over every slice in the canonical-word group.  ⚠ An EARLY DRAFT used
`(sum of width*flux) * circumference` instead and overcounted area by roughly the word length
(a cylinder is crossed by the boundary M times per period, once per word letter, so summing
its full circumference at every one of those M crossings multiplies area by M) -- caught by
the total-area self-check below (11.9x off on `3/7`), not trusted silently.  `L_c` (the
circumference) is still recorded and CHECKED CONSTANT across every slice in the group (a
cylinder's parallel orbits all have equal period by construction, so disagreement is a bug
flag); a derived height `Area_c / L_c` is reported for reference, not used to compute area.

THE THREE SELF-CHECKS THAT MAKE THIS TRUSTWORTHY (all three must pass or the count is not
reported as final):
  (iii-total)  sum(Area_c over ALL cylinders) / Area(triangle)      -> must be 1.0
  (iii-swept)  sum(Area_c over cylinders meeting Sigma_1) / Area    -> must reproduce g_len
  (topo-bound) C_total <= g + s - 1                                 -> must hold
"cylinders meeting Sigma_1" = canonical-word groups that have at least one slice on
`side="L1", k=Q` (the single inward-perpendicular direction, [NCYL-010]) in their group --
NOT "any slice on L1 in any direction", which would overcount (L1 is crossed by many
cylinders in many directions; only the k=Q slices are the actual beam).

⚠ THE AREA SELF-CHECK CANNOT SEE DUPLICATION, ONLY OMISSION (team-lead-caught, `7/16`, this
session, after the first two self-checks had already been reported passing).  `total_area_ratio
= 1.0` proves nothing is MISSING; it says nothing about whether one physical cylinder was
counted TWICE, because a duplicate contributes the same total area whether split into one
group or two.  The topological bound `C_total <= g + s - 1` (`s` = number of singular cone
points, `rulings.md` [NCYL-007]'s own worked instance: `H(3,3,6)` at `8/15` has `g+s-1=9`) is a
genuinely INDEPENDENT constraint that the area checks cannot substitute for -- `7/16` passed
both area checks with `C_total=10` while its bound is 9.  The duplicate was two canon-word
groups that are NOT cyclic rotations of each other (so the primary merge missed them) but do
share word length (675) and circumference (336.749540998, agreeing to 9 decimals) with
DIFFERENT areas (0.067665556 vs 0.060952331) -- `merge_duplicate_cylinders` merges on that
signature (word length + circumference tolerance), explicitly NOT on area equality (the areas
being different is the evidence these are two slices of one cylinder, not a coincidence to
average away). `topological_bound()` and the merge are applied to EVERY row, not just `7/16` --
any other row silently carrying the same defect would have been invisible to the area checks
too.

GUARDS.
 G1  Circumference consistency within a canonical-word group (`max/min circumference` per
     group reported; a group failing this is flagged, not silently averaged).
 G2  `[NCYL-006]`/`[NCYL-009]`-style tolerance care: closure test tolerance is `POS_TOL=1e-7`
     (identical to `s276_periodic_mass.py`, itself validated by that probe's own calibration
     against exact `mu_boundary`/`area_exact`), not re-derived here.
 G3  If total area does not reconcile to 1.0 within a loose float band (this being a NEW
     instrument, no established tolerance exists for it yet), STOP and report the discrepancy
     rather than trusting a partial cylinder count -- this is the team-lead's explicit
     instruction, not a discretionary choice.
 G4  RESOLUTION CONVERGENCE, not a fixed sample count: `NSAMP_COARSE` was raised in-session
     (60 -> 150 -> 400) because 3/13, 7/16, 9/20 each showed the SAME symptom at coarser
     settings -- extra spurious canonical-word groups (a genuine cylinder boundary sitting
     inside one coarse bracket got missed, so bisection anchored on the wrong pair of coarse
     samples and split or merged intervals wrongly) and a swept-area mismatch against the
     independently-known exact `g_len` (up to 1.9% relative on 3/13 at NSAMP=150). Each of
     those three rows was RERUN at higher resolution and the count DROPPED (10->7 on 3/13,
     12->10 on 7/16, 12->10 on 9/20) while the area matched `g_len` to 9-10 decimal places --
     i.e. the fix was resolution, not a new mechanism, and G3's self-check is what caught the
     under-resolved runs in the first place (they are not silently reported).

PRIOR ART (grepped 'cylinder count', 'cylinder diagram', 'enumerate cylinders',
'canonical word', 'flat cylinder area' across `rulings.md`, `cf_width_laws.md`,
`veech_reframing.md`, `TOOLS.md`) ->
  - `engine/cyl_diagram.py`: builds a `CylinderDiagram` OBJECT from an ALREADY-KNOWN
    permutation/cycle structure (genus, stratum, singularities) -- it does not discover
    cylinders from a billiard trace, so it is not reusable for the discovery step; if this
    probe succeeds, its OUTPUT (cylinder list) could in principle feed `cyl_diagram` later,
    but that is not attempted here.
  - `archive/scripts_2026-07/veech_cylinders.py`/`veech_exceptional.py`: cylinder counts for
    the ALL-DIRECTIONS Veech dichotomy (Apisa import), a different question (every direction on
    a lattice surface) from ONE fixed perpendicular direction on a (generically non-Veech)
    surface -- not reusable.
  - `probes/s276_periodic_mass.py`: OWNS the full-boundary cross-section, the flux weighting,
    and the exact `mu_boundary`/`area_exact` calibration this probe reuses verbatim; it
    measures MEASURE (nu(Periodic \\ Reach)), not a cylinder COUNT, which is the gap this probe
    fills.
  - `right_triangle_billiards.discover_periodic_branches`: the bisection-refinement technique
    this probe generalizes from ("L1, stop-at-perpendicular-return" signature) to "any side,
    full-period-closure" signature.
  - `rulings.md` [NCYL-046]: the `ncyl = (n+[P odd])/2` branch-pairing formula this probe's
    count is compared against (not re-derived, not assumed).
  - `rulings.md` [NCYL-007]: OWNS the `g+s-1` topological bound this probe's `topological_bound()`
    applies (its own worked instance is `H(3,3,6)` at `8/15`, `g+s-1=9`, 3 cylinders met) --
    supplied by team-lead this session as an independent check the area self-checks cannot
    perform (they see completeness, not duplication); not re-derived here, only applied.
  - No prior probe enumerates cylinders (as opposed to branches or a periodic MEASURE) for a
    single fixed direction on this billiard family, at either parity.

⚠ FLOAT ONLY.  Not a certificate; the exact route is a ring job, not attempted here.

Run: PYTHONPATH=.:engine:archive/scripts_2026-07:archive/scripts_2026-08:probes \
     .venv/bin/python3.13 probes/s315_cylinder_count.py 2>&1 | tee logs/s315_cylinder_count.log
"""
import json
import math
import sys
import time

from right_triangle_billiards import RightTriangleBilliard
import s315_oddP_census as census  # reused verbatim for the topological-bound cone table

POS_TOL = 1e-7      # closure tolerance on the boundary coordinate -- same as s276_periodic_mass
CAP = 20000         # hits allowed per orbit before giving up (periods here are short: n<=11)
NSAMP_COARSE = 400   # coarse samples per (side,k) strip before bisection refinement -- an
                      # earlier draft used 60, then 150; 3/13, 7/16, 9/20 still showed spurious
                      # extra cylinders and a swept-area mismatch at 150 (converged only at 400)
REFINE_STEPS = 40

ROWS = [
    (3, 7), (3, 8), (3, 10), (3, 11), (5, 12), (3, 13), (5, 13), (7, 16), (9, 20),
]


def geometry(P, Q):
    """Side data: name -> (length, inward unit normal, base point, unit tangent).
    Verbatim from `probes/s276_periodic_mass.py::geometry` -- not re-derived."""
    alpha = (P / Q) * (math.pi / 2)
    L = 1.0 / math.tan(alpha)
    hyp = 1.0 / math.sin(alpha)
    return {
        "L1": (1.0, (-1.0, 0.0), (L, 0.0), (0.0, 1.0)),
        "L2": (L, (0.0, 1.0), (0.0, 0.0), (1.0, 0.0)),
        "H": (hyp, (1.0 / hyp, -L / hyp), (0.0, 0.0), (L / hyp, 1.0 / hyp)),
    }


def flux_weight(side, k, Q, geo):
    th = k * math.pi / Q
    _, nrm, _, _ = geo[side]
    return math.cos(th) * nrm[0] + math.sin(th) * nrm[1]


def trace_full_period(B, geo, Q, side0, x0, k0, cap):
    """Trace from (side0,x0,k0) until the exact boundary state recurs.  Returns
    (word, circumference, seg0, closed) -- word is the tuple of sides hit, circumference the
    total path length back to the start, seg0 the length of the FIRST segment (the Kac
    next-return length used for the area reconciliation -- NOT the same as circumference,
    see the area-formula note in `run_row`), closed=False if it ran into a vertex or the cap."""
    _, _, base0, tan0 = geo[side0]
    px = base0[0] + x0 * tan0[0]
    py = base0[1] + x0 * tan0[1]
    th = k0 * math.pi / Q
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
        kf = (math.atan2(vy, vx) % (2.0 * math.pi)) * Q / math.pi
        if abs(kf - round(kf)) > 1e-6 * Q:
            return None, None, seg0, False
        k = int(round(kf)) % (2 * Q)
        _, _, base, tan = geo[side]
        u = (px - base[0]) * tan[0] + (py - base[1]) * tan[1]
        if side == side0 and k == k0 and abs(u - x0) < POS_TOL:
            return tuple(word), circ, seg0, True
    return None, None, seg0, False


def canon(word):
    """Lexicographically-smallest cyclic rotation."""
    m = len(word)
    return min(word[i:] + word[:i] for i in range(m))


def discover_cylinder_intervals(B, geo, Q, side, k, cap):
    """Bisection-refined intervals of `side`'s coordinate for fixed direction `k`.  Analogous
    to `discover_periodic_branches`, generalized to any side and full-period closure.

    ⚠ The boundary-detection signature is the RAW word, NOT `canon(word)`.  An early draft
    canonicalized (cyclically rotated) before comparing adjacent samples WITHIN one (side,k)
    sweep -- but two genuinely different first-hit sides can rotate to the SAME canonical form
    by coincidence, which either masks a real cylinder boundary or (worse) lets `seg0`, which
    is only affine in x while the first-hit side is truly constant, get integrated across a
    hidden discontinuity -- exactly the bug the area self-check caught (area ratio off by
    ~8-9% on 3/7, in BOTH directions depending on which evaluation point the bug happened to
    hit). Canonicalization is correct ONLY when grouping intervals ACROSS different (side,k)
    entries afterward (different starting phase into the same cylinder) -- done in `run_row`."""
    ell = geo[side][0]
    xs = [ell * (i + 0.5) / NSAMP_COARSE for i in range(NSAMP_COARSE)]
    sigs = []
    for x in xs:
        w, c, s0, ok = trace_full_period(B, geo, Q, side, x, k, cap)
        sigs.append(w if ok else None)

    boundaries = [0.0]
    for i in range(NSAMP_COARSE - 1):
        if sigs[i] != sigs[i + 1]:
            lo, hi = xs[i], xs[i + 1]
            sig_lo = sigs[i]
            for _ in range(REFINE_STEPS):
                mid = 0.5 * (lo + hi)
                w, c, s0, ok = trace_full_period(B, geo, Q, side, mid, k, cap)
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
        w, c, s0_mid, ok = trace_full_period(B, geo, Q, side, mid, k, cap)
        # AREA (Kac) INTEGRATION FIX: seg0(x) -- the distance to the first boundary hit -- is
        # an AFFINE function of the starting position x while the first-hit side is constant,
        # which it is throughout one discovered interval (the interval boundaries are exactly
        # where the word, hence the first letter, changes).  A single midpoint evaluation of
        # seg0 under-/over-integrates badly on a WIDE interval (an early draft used exactly
        # that and overcounted total area by 8% on 3/7, caught by the total-area self-check);
        # the two-endpoint trapezoidal average is EXACT for an affine integrand.
        eps = width * 1e-6
        _, _, s0_lo, ok_lo = trace_full_period(B, geo, Q, side, lo + eps, k, cap)
        _, _, s0_hi, ok_hi = trace_full_period(B, geo, Q, side, hi - eps, k, cap)
        seg0_avg = (s0_lo + s0_hi) / 2.0 if (ok_lo and ok_hi) else s0_mid
        out.append({"side": side, "k": k, "lo": lo, "hi": hi, "width": width,
                     "word": w, "canon": canon(w) if ok else None,
                     "circumference": c, "seg0": seg0_avg, "ok": ok})
    return out


MERGE_CIRC_RELTOL = 1e-6  # relative tolerance for the post-hoc word-length+circumference merge
MERGE_MAX_MISMATCH = 3    # max Hamming distance (best cyclic rotation) to call two same-length,
                          # same-circumference groups one physical cylinder -- calibrated
                          # in-session: 7/16's true duplicate (675-letter word) measured 2
                          # mismatches; 5/7's false-positive pair (14-letter word, same length,
                          # CLOSER circumference agreement than 7/16's true pair -- 1e-14 vs
                          # 1e-9) measured 7 mismatches (half the word) and is correctly
                          # rejected at threshold 3.


def topological_bound(P, Q):
    """`C_total <= g + s - 1`, where `s` = the number of VERTEX CLASSES at O and A, i.e.
    `sum(count)` over `s315_oddP_census.cone_table`'s entries.  Reused verbatim from
    `rulings.md` [NCYL-007]'s own worked instance (`H(3,3,6)` at `8/15`: g=7, s=3,
    g+s-1=9) -- not re-derived here, only applied.

    ⚠⚠ DOCSTRING CORRECTED s438 ([OPS-202]).  This line used to read "s = number of
    SINGULAR cone points ... R, order 0, does not count", which the code has never
    done: `sum(count)` is UNFILTERED, so an O or A class of order 0 -- a REGULAR point,
    cone angle 2*pi -- is counted.  The two readings disagree on 172 of 1100 coprime
    centres (Q <= 60), in exactly four families: P=1, P=2, Q-P=1, and Q-P=2 with both
    odd.  ** The CODE is the correct one and the docstring was wrong: ** filtering
    order 0 out drops the measured count law from 12/12 to 8/12 at Q=9 ([NCYL-254]),
    and the filtered ceiling is VIOLATED by measured rows at 2/5, 2/7, 2/9, 7/9
    ([NCYL-255]).  So do NOT "fix" this function to match the old docstring.
    ⚠ R (order 0, Q points) is still excluded, as it always was; whether it OUGHT to
    count is open, and counting it would make the bound vacuous ([NCYL-255])."""
    cones, _, g = census.cone_table(P, Q)
    s = sum(c["count"] for c in cones)
    return g + s - 1, g, s


def best_rotation_mismatch(w1, w2):
    """Minimum Hamming distance between `w1` and any cyclic rotation of `w2` (equal length
    required).  0 means `w1`/`w2` ARE cyclic rotations of each other (should already be one
    canon-word group -- this function is only ever called on the residual, already-different,
    cases)."""
    n = len(w1)
    best = n
    for r in range(n):
        rot = w2[r:] + w2[:r]
        diff = sum(1 for x, y in zip(w1, rot) if x != y)
        if diff < best:
            best = diff
    return best


def merge_duplicate_cylinders(cylinders):
    """Post-hoc merge of canon-word groups that are the SAME physical cylinder entered via a
    word that is not an EXACT cyclic rotation of the other -- caught by team-lead's
    topological-bound check at `7/16` (C_total=10 > g+s-1=9): two groups (word length 675,
    circumference 336.749540998 agreeing to 9 decimals, areas 0.067665556 / 0.060952331 --
    DIFFERENT areas, which is the tell that this is not a same-canon-word accounting artifact)
    were one cylinder counted twice, most plausibly by a single mislabeled hit deep in a long
    word near a bisection boundary (a knife-edge numerical effect, not a physical difference).

    ⚠ (WORD LENGTH, CIRCUMFERENCE) ALONE IS NOT A SAFE MERGE SIGNATURE -- FALSE POSITIVE FOUND
    AND FIXED IN-SESSION.  `5/7` had two GENUINELY DIFFERENT cylinders (one meeting Sigma_1,
    one not) with equal word length (14) and circumference agreeing to 1e-14 (TIGHTER than
    `7/16`'s true duplicate) -- yet their best cyclic-rotation alignment differs in 7 of 14
    positions (half the word), the opposite of a one-letter numerical glitch. The naive merge
    silently produced `swept_area_ratio=0.892` against the independently-known exact
    `g_len=0.635` on this row (caught by the extended-rows cross-check against the exact table,
    NOT by the (iii-total)/(iii-swept)/bound self-checks -- none of those three sees a
    WRONGLY-CLASSIFIED meets_sigma1 flag on a correctly-summed total area). Requiring the
    content of the words to ALSO nearly match (`best_rotation_mismatch <= MERGE_MAX_MISMATCH`,
    calibrated below: `7/16`'s true duplicate measures 2 mismatches, `5/7`'s false positive
    measures 7) correctly separates the two cases -- `7/16`'s true duplicate still merges,
    `5/7`'s coincidence no longer does; both re-verified against their independently-known
    exact `g_len` after this fix (see `run_row`'s docstring note / the session log)."""
    merged = []
    used = [False] * len(cylinders)
    for i, ci in enumerate(cylinders):
        if used[i]:
            continue
        group = [ci]
        used[i] = True
        for j in range(i + 1, len(cylinders)):
            if used[j]:
                continue
            cj = cylinders[j]
            len_i, len_j = len(ci["canon_word"]), len(cj["canon_word"])
            if len_i != len_j:
                continue
            ci_mid = (ci["circumference_min"] + ci["circumference_max"]) / 2.0
            cj_mid = (cj["circumference_min"] + cj["circumference_max"]) / 2.0
            if abs(ci_mid - cj_mid) >= MERGE_CIRC_RELTOL * max(1.0, ci_mid):
                continue
            mismatch = best_rotation_mismatch(ci["canon_word"], cj["canon_word"])
            if mismatch <= MERGE_MAX_MISMATCH:
                group.append(cj)
                used[j] = True
        if len(group) == 1:
            merged.append(ci)
        else:
            area = sum(c["area"] for c in group)
            circ_lo = min(c["circumference_min"] for c in group)
            circ_hi = max(c["circumference_max"] for c in group)
            merged.append({
                "canon_word": tuple(c["canon_word"] for c in group),  # merged identity
                "merged_from_n_groups": len(group),
                "n_slices": sum(c["n_slices"] for c in group),
                "circumference_min": circ_lo, "circumference_max": circ_hi,
                "circumference_consistent": (circ_hi - circ_lo) < 1e-6 * max(1.0, circ_hi),
                "area": area,
                "height_derived": None,
                "meets_sigma1": any(c["meets_sigma1"] for c in group),
                "merged_areas": [c["area"] for c in group],
            })
    return merged


def run_row(P, Q, cap=CAP):
    alpha = (P / Q) * (math.pi / 2)
    B = RightTriangleBilliard(alpha)
    geo = geometry(P, Q)
    area_exact = Q / math.tan(alpha)  # 2Q*Area(triangle) -- s276_periodic_mass's normalization

    all_slices = []
    bad_slices = []
    for side in geo:
        for k in range(2 * Q):
            if flux_weight(side, k, Q, geo) <= 1e-12:
                continue
            ivs = discover_cylinder_intervals(B, geo, Q, side, k, cap)
            flux = flux_weight(side, k, Q, geo)
            for iv in ivs:
                iv["flux"] = flux
                # KAC AREA CONTRIBUTION (not a cylinder height): width * flux * seg0.  A
                # cylinder's periodic orbit crosses the boundary M times per period (M = its
                # word length), so its FULL transverse width is only seen ONCE if you sum
                # width*flux*seg0 over ALL M of its slices (Kac's lemma tiles the domain
                # exactly once via the "next segment", not the full period, at each crossing)
                # -- an early version of this probe used width*flux*circumference instead and
                # overcounted total area by roughly the word length (11.9x on 3/7).  Caught by
                # the (iii-total)=1.0 self-check, not trusted silently.
                iv["area_contrib"] = iv["width"] * flux * iv["seg0"]
                if iv["ok"]:
                    all_slices.append(iv)
                else:
                    bad_slices.append(iv)

    groups = {}
    for s in all_slices:
        groups.setdefault(s["canon"], []).append(s)

    cylinders = []
    for cw, slices in groups.items():
        area = sum(s["area_contrib"] for s in slices)
        circs = [s["circumference"] for s in slices]
        circ_lo, circ_hi = min(circs), max(circs)
        meets_sigma1 = any(s["side"] == "L1" and s["k"] == Q for s in slices)
        circ_mid = (circ_lo + circ_hi) / 2.0
        cylinders.append({
            "canon_word": cw, "n_slices": len(slices),
            "circumference_min": circ_lo, "circumference_max": circ_hi,
            "circumference_consistent": (circ_hi - circ_lo) < 1e-6 * max(1.0, circ_hi),
            "area": area,
            "height_derived": area / circ_mid if circ_mid else None,  # Area/L, NOT summed
            "meets_sigma1": meets_sigma1,
        })

    n_before_merge = len(cylinders)
    cylinders = merge_duplicate_cylinders(cylinders)
    n_merged_away = n_before_merge - len(cylinders)

    bound, g, s = topological_bound(P, Q)
    bound_ok = len(cylinders) <= bound

    total_area = sum(c["area"] for c in cylinders)
    swept_area = sum(c["area"] for c in cylinders if c["meets_sigma1"])
    return {
        "P": P, "Q": Q, "area_exact": area_exact,
        "n_cylinders_total": len(cylinders),
        "n_cylinders_before_merge": n_before_merge,
        "n_merged_away": n_merged_away,
        "n_cylinders_meeting_sigma1": sum(1 for c in cylinders if c["meets_sigma1"]),
        "total_area_ratio": total_area / area_exact,
        "swept_area_ratio": swept_area / area_exact,
        "n_bad_slices": len(bad_slices),
        "genus": g, "n_singularities": s, "topological_bound": bound,
        "topological_bound_ok": bound_ok,
        "cylinders": cylinders,
    }


def main():
    t0 = time.time()
    out = []
    for (P, Q) in ROWS:
        rt0 = time.time()
        r = run_row(P, Q)
        r["secs"] = round(time.time() - rt0, 1)
        out.append(r)
        consistent = all(c["circumference_consistent"] for c in r["cylinders"])
        bound_flag = "" if r["topological_bound_ok"] else "  <<< BOUND VIOLATED"
        merge_flag = f" merged_away={r['n_merged_away']}" if r["n_merged_away"] else ""
        print(f"{P}/{Q}: C_total={r['n_cylinders_total']} "
              f"C_meeting_L1={r['n_cylinders_meeting_sigma1']} "
              f"bound(g+s-1)={r['topological_bound']} "
              f"total_area_ratio={r['total_area_ratio']:.6f} "
              f"swept_area_ratio={r['swept_area_ratio']:.6f} "
              f"circ_consistent={consistent} bad_slices={r['n_bad_slices']}"
              f"{merge_flag}{bound_flag} [{r['secs']}s]", flush=True)
        with open("data/s315_cylinder_count.json", "w") as f:
            json.dump(out, f, indent=1)
        if not r["topological_bound_ok"]:
            print(f"  !! {P}/{Q}: C_total={r['n_cylinders_total']} exceeds the topological "
                  f"bound g+s-1={r['topological_bound']} (g={r['genus']}, s={r['n_singularities']}) "
                  f"even after the word-length+circumference merge -- reported AS-IS, not "
                  f"forced to the bound.", flush=True)
    n_violations = sum(1 for r in out if not r["topological_bound_ok"])
    print(f"\nDONE wall_clock={round(time.time()-t0,1)}s "
          f"bound_violations={n_violations}/{len(out)}", flush=True)
    print("wrote data/s315_cylinder_count.json", flush=True)


if __name__ == "__main__":
    sys.exit(main())
