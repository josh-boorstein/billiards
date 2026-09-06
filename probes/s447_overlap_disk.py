#!/usr/bin/env python3
"""
s447_overlap_disk.py -- PRE-REGISTERED.  Queue items (g') and (a'').

WHAT THIS IS.  A PROOF of [NCYL-093]'s SECOND measured half -- *overlap subset {orphan}*,
i.e. at most ONE cylinder is met by BOTH perpendicular transversals -- plus the instrument
that scores the geometric MODEL the proof rests on (not the proof's conclusion).
⚠ It does NOT prove the COVERING half (`p = 0` / (M1)); see "WHAT THIS DOES NOT DO".

PRIOR ART: grepped 'overlap', 'met_H', 'met_L2', 'non-crossing', 'nested', 'straddl',
'Gauss', 'maximal cylinder', 'double normal', 'at most one shared'; read `rulings.md`
[NCYL-071]/[NCYL-072]/[NCYL-076]/[NCYL-078]/[NCYL-080]/[NCYL-093]/[NCYL-100]/[NCYL-107]/
[NCYL-114]/[NCYL-123]/[NCYL-132]/[NCYL-236]/[NCYL-280]/[NCYL-281]/[NCYL-282]/[NCYL-283]) ->
  - [NCYL-093] IS the target: "overlap subset {orphan}" and the COVERING half are the two
    MEASURED halves; "orphan IN overlap" is PROVED ([NCYL-080] + [NCYL-132]).
  - [NCYL-123] is the audit that says the covering half has "no argument, no candidate and
    no named attack anywhere" -- and it says the same of §s316 (9)'s claim that the "at most
    one shared" half "has a candidate argument": *the argument is nowhere in the repo*.
    So this file is that argument, and it is a DIFFERENT argument from the covering one.
  - [NCYL-236] (s404) measured both set laws per cylinder on `58` rows, `Q <= 20`, and is
    the row set H3/H4 must agree with; `probes/s404_second_strip.py` owns that reading.
  - [NCYL-100]/[NCYL-107] ring-certify the IDENTITY `n + n_H = Q+1`, explicitly NOT the
    covering argument -- so no certification here is inherited.
  - No probe in the repo reads the ORDER of the cells along a strip against the order of
    the other strip, or tests any planarity/non-crossing property of the pairing.

------------------------------------------------------------------------------------
THE PROOF.  Everything is on [NCYL-283]'s picture, and nothing here re-derives it:
`B = Sigma/tau` is the genus-0 base ([NCYL-280]), `iota` the involution induced by the
class reflection, `Fix(iota)` a circle made of the `h` horizontal and `v` vertical
side-copies, and `B \\ Fix(iota) = D t D'`, two disks swapped by `iota` ([NCYL-282]).

(0) THE THREE CORNERS OF `D` ARE THE THREE VERTICES, SO `int D` CARRIES NO ZERO.
    `Fix(iota)`'s arcs are side-copies and its corners are vertex-copies, so `Z_O` and
    `Z_A` -- the only two zeros of `B` -- lie ON `Fix(iota)`.  Hence **every singularity
    in `int D` is a SIMPLE POLE** (an `R`-copy, cone angle `pi`).  Gauss-Bonnet closes on
    the pole count in all three parity classes, which is the arithmetic check below:
        odd `Q`  (h,v)=(1,2) or (2,1): corners `Z_O`,`Z_A`,`R_*`   -> `n_int = (Q-1)/2`
        even `Q`, eps even   (h=v=2) : corners `Z_O`,`Z_A`,`R_1`,`R_2` -> `n_int = Q/2 - 1`
        even `Q`, eps odd    (h=v=1) : corners `Z_O`,`Z_A`         -> `n_int = Q/2`
    and `2*n_int + #(poles on Fix)` is `Q` in every one of them.

(1) A CYLINDER THAT MEETS `Fix(iota)` MEETS IT IN EXACTLY TWO VERTICAL SEGMENTS.
    `iota` maps cylinders to cylinders.  It cannot preserve one freely: a free-invariant
    cylinder is connected and misses `Fix`, so it lies in `D`, but then `iota` sends it
    into `D'`.  It cannot fix one whose `Fix`-set is HORIZONTAL: that would put a closed
    leaf inside `Fix(iota)`, and every arc of `Fix(iota)` ends at a singularity.  So an
    `iota`-invariant cylinder is fixed with `iota|_Z : (x,t) -> (-x,t)` in its own flat
    coordinate: **two full vertical segments crossing it**, cutting it into two rectangles
    swapped by `iota`, one of which is `Z ^ D`.  (Cylinders meeting `Fix` at all is the
    `p = 0` question -- NOT assumed anywhere below.)

(2) SO "MET BY THE PERPENDICULAR BEAM OF `sigma`" = "HAS A `Fix`-SEGMENT ON `sigma`'s
    VERTICAL ARC" ([NCYL-283] (4): the points of `G` on that arc are the beam's interior
    cell boundaries), and **the OVERLAP is the set of cylinders whose two segments lie on
    DIFFERENT vertical arcs.**  Call these STRADDLERS and let `s = #straddlers`.

(3) THE STRADDLERS ARE NESTED.  Their rectangles are disjoint subsets of the disk `D`,
    each with its two vertical sides on the two different vertical arcs, so the chords
    joining their vertical sides are disjoint arcs in a disk with endpoints on the
    boundary: a NON-CROSSING family, and every straddler separates the two horizontal
    boundary arcs from each other.  Hence they are linearly ordered.

(4) ⇒⇒ **`s <= 1`.**  Suppose `Rect_1`, `Rect_2` are consecutive straddlers and let `U`
    be the region of `D` between them.
    * If `U` is empty, `top(Rect_1) = bottom(Rect_2)` is a horizontal arc bounding two
      DISTINCT maximal cylinders, so it carries a singularity.  Its endpoints lie in the
      INTERIORS of the two vertical arcs (they are strictly between the straddlers'
      other corners), which are regular; and an interior singularity of the arc is a
      simple pole by (0), whose total cone angle `pi` cannot supply the angle `pi` that
      EACH of the two cylinders needs there.  Contradiction.
    * Otherwise `U` is a disk with four corners of interior angle `pi/2` (a vertical arc
      of `Fix` meets a horizontal cylinder boundary orthogonally) and no other
      singularity on its two vertical sides.  Gauss-Bonnet,
          Sum_{int}(2pi - theta) + Sum_{bdry}(pi - theta) = 2pi*chi(U) = 2pi,
      spends the whole `2pi` on those four corners (`4 * pi/2`).  Every remaining term is
      >= 0: an interior singularity is a pole (`2pi - pi = +pi`), and a pole on one of the
      two horizontal sides has `U`-side angle `theta <= pi` (`+pi - theta >= 0`); no zero
      can appear, by (0).  So **`U` contains no singularity at all** -- and then
      `top(Rect_1)` is a singularity-free horizontal arc whose endpoints are regular, so
      the cylinder boundary component `top(Rect_1) u iota(top(Rect_1))` is a closed leaf
      with no singularity on it, which no MAXIMAL cylinder has.  Contradiction.
    ⇒ no two straddlers, `s <= 1`.

(5) AND `s >= 1` IS ALREADY PROVED: the orphan retraces, so it is perpendicular to `L1`
    at one end and to the second side at its temporal midpoint ([NCYL-072] + Lemma B,
    vertex gap DISCHARGED at [NCYL-132]), i.e. it has one segment on each vertical arc.
    ⇒⇒ **`overlap = {orphan}` EXACTLY, in every class that has two vertical arcs.**
    This closes [NCYL-093]'s "overlap subset {orphan}" half and [NCYL-236]'s even-`Q`
    sibling `met_L1 ^ met_L2 = {orphan}`.

WHAT THIS DOES NOT DO.  It does not touch `p = 0` (M1).  The same Gauss-Bonnet on the
disk `W` cut off by a cylinder interior to `D` returns `#poles(W) = 2` -- which is
[NCYL-283]'s own already-derived failure shape (an interior `R`-`R` saddle connection),
not a contradiction.  The two halves of [NCYL-093] are therefore NOT symmetric: the
overlap half is a nesting/maximality argument and closes; the covering half is an
existence statement about every pole's prong and does not.

------------------------------------------------------------------------------------
HYPOTHESES (pre-registered).  ⚠ H2 was pre-registered as the model test and turned out to be
VACUOUS -- it is kept, scored, and quoted as nothing ([OPS-219]); H6 replaces it and is the one
that can fail.

  H1  CONTROL, and it validates the reading before anything is concluded from it.  In the
      concatenated cell sequence along the two vertical arcs, EVERY canonical word occurs
      exactly TWICE.  This is (1) at the level of individual cylinders; the aggregate
      `sum n_sigma = 2*C` is [NCYL-093]'s law, but WHICH cells pair with which is not, and
      a mis-paired strip breaks this without breaking the total.
  H2  ⚠⚠ SCORED BUT **NOT EVIDENCE**, AND THE PROBE SAYS SO IN ITS OWN OUTPUT.  Read in the
      geometric order the matching is NON-CROSSING.  It looked like the model test -- a
      random `2C`-point matching is non-crossing with probability `~C_C/(2C-1)!!`, `1/135135`
      at `C = 8` -- but it is IMPLIED by structure this repo already asserts: s439
      `boundary_map`'s docstring records that a cylinder's TWO crossings of a strip are
      *"separated INSIDE by an `R` hit"*, i.e. they are ADJACENT cells (H7), and a matching
      of adjacent pairs plus one long chord is non-crossing in EVERY arc order.  Its
      negative control confirms exactly that: all three WRONG orders also score `100%`.
      ⇒ [OPS-041]: kept, reported, and quoted as nothing.
  H3  `s == 1`.  ⚠ THIS IS THE THEOREM'S CONCLUSION, so it is a CONSISTENCY CHECK, not
      evidence for the proof -- EXCEPT on the `{L2,H}` class (even `P`, odd `Q`), where the
      per-cylinder overlap law has never been measured at all: [NCYL-236]/`s404` excludes
      even `P` by scope (*"Even `P` has NO second transversal"* -- true of the `Sigma_L1`
      FLOW direction), and [NCYL-114] gets that column's COUNT law by transport from the
      partner row, never from a second per-cylinder flag.  Those rows can fail.
  H4  the straddler is the ORPHAN, identified INDEPENDENTLY by [NCYL-076]'s odd-letter
      signature ([OPS-071]): a retracing word is `hinge + u + hinge + reverse(u)`, so the
      orphan is the unique cylinder whose word has ODD counts of BOTH hinge sides and an
      EVEN count of the third.  ⇒ the hinges are predicted to be exactly the two VERTICAL
      sides -- `{L1,H}`, `{L1,L2}` are [NCYL-076]'s own two cases and `{L2,H}` is new.
  H5  the pole arithmetic of (0): `2*n_int + #poles_on_Fix == Q` in each parity class,
      pure arithmetic over the whole box.  An off-by-one in the corner set breaks it.
  H6  ⇒⇒ THE MODEL TEST THAT CAN FAIL, AND IT IS THE PROOF'S OWN BOOKKEEPING RUN ON A
      REGION IT DOES NOT NEED.  The straddler cuts `D` into two pieces; Gauss-Bonnet on
      each gives `#poles = sum_{corners on that piece}(theta_i/pi - 1) + 1`, and every
      non-straddling cylinder wraps exactly ONE interior pole (its two adjacent cells are
      the two sides of that pole's prong), so the two piece-populations are FORCED:
          m_in  = #cylinders strictly between the straddler and the arcs' junction,
          m_out = #cylinders on the other side,
      with `m_in`/`m_out` read off the corner angles alone (`Z_O -> P/2`, `Z_A -> (Q-P)/2`,
      `R -> 1/2`, in units of `pi`).  This PINS THE ORPHAN'S POSITION on both strips --
      `(m_in, m_out) = ((Q-P)/2, (P-1)/2)` for `{L1,H}`, `(P/2, (Q-P-1)/2)` for `{L2,H}`,
      `((Q-P-1)/2, (P-1)/2)` for `{L1,L2}` -- which no counting law implies (they all
      constrain only `m_in + m_out + 1 = C`).  A shift of the straddler by one cylinder
      breaks it and preserves every law in [NCYL-093]/[NCYL-114]/[NCYL-241].
  H7  CONTROL (prior art, s439): every NON-straddling cylinder's two cells are ADJACENT on
      their arc, and the straddler's are not.  Scored because H2 rests on it.

ARC ORDER (derived, not fitted -- getting it wrong is what the negative control detects).
`geo` puts `L1` at base `R` tangent towards `A`, `L2` at base `O` towards `R`, `H` at base
`O` towards `A`.  `Fix(iota)`'s cyclic order is the triangle's own, `L1 - A - H - O - L2 -
R - L1`, except at even `Q` where each side supplies BOTH a horizontal and a vertical copy
and the cycle is `L1^v - Z_A - L1^h - R_2 - L2^v - Z_O - L2^h - R_1` (at an `R` corner the
two sides are perpendicular, so exactly one of the two copies there is vertical).  Reading
the two vertical arcs head-to-tail along that cycle gives
    odd  `Q`, verticals {L1,H} : L1 increasing `w`,  then H decreasing `w`
    odd  `Q`, verticals {L2,H} : H  decreasing `w`,  then L2 increasing `w`
    even `Q`, verticals {L1,L2}: L1 increasing `w`,  then L2 decreasing `w`
(the whole sequence may be reversed; non-crossing and `s` are invariant under that).

GUARDS.
  G1  a class with any prong not terminating at `cap2` is UNSCORABLE (`PRONG_CAPPED`),
      s439's convention -- the CP scope gate ([NCYL-283], [NCYL-262]).
  G2  a cell whose closure trace did not close (`canon is None`) makes the row UNSCORABLE
      rather than silently dropping a chord endpoint.
  G3  nothing here writes to any store.

Run: PYTHONPATH=.:engine:archive/scripts_2026-07:archive/scripts_2026-08:probes \\
     .venv/bin/python3.13 probes/s447_overlap_disk.py [h5|sweep] [--qmax=N] [--budget=S]
"""
import json
import math
import sys
import time
from collections import Counter

from right_triangle_billiards import RightTriangleBilliard
import s315_cylinder_count as s315
import s437_oblique_cylinders as s437
import s439_exact_cells as s439
import s446_covering_identity as s446

CAP = 200_000
CAP_BIG = 2_000_000


# ---------------------------------------------------------------- H5: pole arithmetic

def geometry_of_class(P, Q, eps):
    """`(h_sides, v_sides, n_int, poles_on_fix)` -- the disk's singularity budget from
    Gauss-Bonnet, `n_int` = simple poles interior to `D`.

    A disk with corner angles `theta_i` and `n_int` interior cone points of angle `pi`
    satisfies `n_int*pi + sum(pi - theta_i) = 2pi`.  The corners are the vertex-copies
    where consecutive arcs of `Fix(iota)` meet; `Z_O` has angle `P*pi/2` in `D`, `Z_A`
    `(Q-P)*pi/2`, and an `R`-corner `pi/2`."""
    hor, ver = s446.hv_sides(P, Q, eps)
    arcs = len(hor) + len(ver)
    if arcs == 3:                      # odd Q: corners Z_O, Z_A, R_*
        r_corners = 1
    elif arcs == 4:                    # even Q, eps even: corners Z_O, Z_A, R_1, R_2
        r_corners = 2
    elif arcs == 2:                    # even Q, eps odd: corners Z_O, Z_A only
        r_corners = 0
    else:
        raise AssertionError(f"{P}/{Q} eps={eps}: {arcs} arcs")
    # n_int = 2 - sum(1 - theta_i/pi) over corners
    defect = (1 - P / 2.0) + (1 - (Q - P) / 2.0) + r_corners * 0.5
    n_int = 2.0 - defect
    return hor, ver, n_int, r_corners


def h5_control(qmax=120, verbose=True):
    """`n_int` is a non-negative integer and `2*n_int + r_corners == Q`, whole box."""
    bad, rows = [], 0
    for Q in range(3, qmax + 1):
        for P in range(1, Q):
            if math.gcd(P, Q) != 1:
                continue
            for eps in (0, 1):
                rows += 1
                _h, _v, n_int, rc = geometry_of_class(P, Q, eps)
                if abs(n_int - round(n_int)) > 1e-9 or round(n_int) < 0 \
                        or 2 * round(n_int) + rc != Q:
                    bad.append((P, Q, eps, n_int, rc))
    if verbose:
        print(f"  H5  2*n_int + #poles_on_Fix == Q : {rows - len(bad)}/{rows}"
              f" class rows (coprime, Q<={qmax})", flush=True)
        for r in bad[:6]:
            print("     ", r, flush=True)
    return bad, rows


# ---------------------------------------------------------------- the arc order

def arc_reading(ver, P, Q):
    """`[(side, reverse), (side, reverse)]` -- the two vertical arcs head-to-tail along
    `Fix(iota)`, `reverse=True` meaning decreasing `w`.  Derived in ARC ORDER above."""
    v = set(ver)
    if v == {"L1", "H"}:
        return [("L1", False), ("H", True)]
    if v == {"L2", "H"}:
        return [("H", True), ("L2", False)]
    if v == {"L1", "L2"}:
        return [("L1", False), ("L2", True)]
    raise AssertionError(f"unexpected vertical set {ver}")


WRONG_ORDERS = {
    "flip2": lambda rd: [rd[0], (rd[1][0], not rd[1][1])],
    "flip1": lambda rd: [(rd[0][0], not rd[0][1]), rd[1]],
    "swap": lambda rd: [rd[1], rd[0]],
}


def _noncrossing(seq):
    """`seq` = labels, each occurring twice.  True iff the induced perfect matching is
    non-crossing (a bracket-matching test on a stack)."""
    stack = []
    for x in seq:
        if stack and stack[-1] == x:
            stack.pop()
        else:
            stack.append(x)
    return not stack


def _straddlers(seq, n_first):
    """Labels with one occurrence in `seq[:n_first]` and one in `seq[n_first:]`."""
    a, b = Counter(seq[:n_first]), Counter(seq[n_first:])
    return sorted(k for k in a if b.get(k))


THETA = {"O": lambda P, Q: P / 2.0, "A": lambda P, Q: (Q - P) / 2.0,
         "R": lambda P, Q: 0.5}


def _hinges(canon_letters, ver):
    """[NCYL-076]'s signature: ODD counts on both hinge sides, EVEN on the third.
    Predicted hinges = the two VERTICAL sides."""
    c = Counter(canon_letters)
    return all(c[s] % 2 == 1 for s in ver) and \
        all(c[s] % 2 == 0 for s in ("L1", "L2", "H") if s not in ver)


def corner_split(ver, P, Q):
    """`(inner_corners, outer_corners)` -- the corners of `Fix(iota)` on each side of a
    straddler, along the arc order of `arc_reading`.  Derived in ARC ORDER + H6."""
    v = set(ver)
    if v == {"L1", "H"}:                 # Gamma = R_* -L1- A -H- Z_O
        return ["A"], ["R", "O"]
    if v == {"L2", "H"}:                 # Gamma = Z_A -H- O -L2- R_*
        return ["O"], ["R", "A"]
    if v == {"L1", "L2"}:                # Gamma = R_1 -L1^v- Z_A -L1^h- R_2 -L2^v- Z_O
        return ["A", "R"], ["O", "R"]
    raise AssertionError(ver)


def m_predicted(ver, P, Q):
    """`(m_in, m_out)` from the corner angles alone -- H6."""
    inn, out = corner_split(ver, P, Q)
    f = lambda cs: sum(THETA[c](P, Q) - 1 for c in cs) + 1
    return f(inn), f(out)


# ---------------------------------------------------------------- the measurement

def row(P, Q, cap=CAP, cap2=CAP_BIG):
    t0 = time.time()
    eps = next(e for e in (0, 1) if len(s446.hv_sides(P, Q, e)[1]) == 2)
    hor, ver, n_int, rc = geometry_of_class(P, Q, eps)
    out = {"P": P, "Q": Q, "eps": eps, "lone": eps == P % 2, "hor": hor, "ver": ver,
           "n_int_poles": int(round(n_int))}

    bmap, reasons, _nb, _n1, aborted = s439.boundary_map(
        P, Q, eps, cap=cap, cap2=cap2, abort_on_open=True)
    if aborted or any(not r["reason"].startswith("vertex") for r in reasons.values()):
        out["status"] = "PRONG_CAPPED"
        out["secs"] = round(time.time() - t0, 1)
        return out

    geo = s437.geometry(P, Q)
    B = RightTriangleBilliard((P / Q) * (math.pi / 2))
    cells = {}
    for side in ver:
        u = s437.perp_index(side, P, Q)
        cl, _mg = s439.strip_cells(B, geo, Q, side, u, bmap[(side, u)], cap)
        if any(c["canon"] is None for c in cl):
            out["status"] = "CELL_OPEN"
            out["secs"] = round(time.time() - t0, 1)
            return out
        cells[side] = cl                       # already sorted by `lo`

    reading = arc_reading(ver, P, Q)

    def build(rd):
        seq = []
        for side, rev in rd:
            lab = [str(c["canon"]) for c in cells[side]]
            seq += lab[::-1] if rev else lab
        return seq

    seq = build(reading)
    n_first = len(cells[reading[0][0]])
    out["n_sigma"] = {side: len(cells[side]) for side in ver}
    out["sum_n"] = len(seq)

    cnt = Counter(seq)
    out["H1_ok"] = all(v == 2 for v in cnt.values())
    out["C_met"] = len(cnt)
    if not out["H1_ok"]:
        out["status"] = "H1_FAIL"
        out["mult"] = sorted(Counter(cnt.values()).items())
        out["secs"] = round(time.time() - t0, 1)
        return out

    out["H2_ok"] = _noncrossing(seq)
    out["negctl"] = {k: _noncrossing(build(f(reading))) for k, f in WRONG_ORDERS.items()}
    strad = _straddlers(seq, n_first)
    out["s"] = len(strad)
    out["H3_ok"] = len(strad) == 1

    orph = sorted(c for c in cnt if _hinges(_canon_letters(c), ver))
    out["n_orphan_candidates"] = len(orph)
    out["H4_ok"] = len(orph) == 1 and orph == strad

    # H6: the straddler's POSITION, from the corner angles alone.
    if out["H3_ok"]:
        i, j = (k for k, x in enumerate(seq) if x == strad[0])
        out["m_in"], out["m_out"] = (j - i - 1) / 2.0, (len(seq) - (j - i + 1)) / 2.0
        mi, mo = m_predicted(ver, P, Q)
        out["m_pred"] = [mi, mo]
        out["H6_ok"] = (out["m_in"], out["m_out"]) == (mi, mo)
        # NEGATIVE CONTROLS for H6, both able to be right only by coincidence:
        # (a) the two pieces' predictions SWAPPED, (b) the straddler moved one cylinder.
        out["neg_H6_swap"] = (out["m_in"], out["m_out"]) == (mo, mi)
        out["neg_H6_shift"] = (out["m_in"], out["m_out"]) == (mi + 1, mo - 1)
        # H7: every non-straddler occupies two ADJACENT cells of one arc.
        blocks, k0 = [], 0
        for side, _rev in reading:
            blocks.append(seq[k0:k0 + len(cells[side])])
            k0 += len(cells[side])
        adj = all(b[t] != b[t + 1] or True for b in blocks for t in range(len(b) - 1))
        bad_adj = []
        for lab in cnt:
            if lab == strad[0]:
                continue
            pos = [k for k, x in enumerate(seq) if x == lab]
            if pos[1] - pos[0] != 1:
                bad_adj.append(lab)
        out["H7_ok"] = not bad_adj and adj
        out["n_nonadjacent"] = len(bad_adj)
    out["status"] = "ok"
    out["secs"] = round(time.time() - t0, 1)
    return out


def _canon_letters(canon_repr):
    """`s315.canon` returns a tuple of side names; `str()` of it is the label used as the
    chord key.  Recover the letters without re-tracing."""
    return [t.strip().strip("'\"") for t in
            canon_repr.strip("()").split(",") if t.strip()]


def sweep(qmax=20, out_path="data/s447_overlap.json", budget_s=None, verbose=True):
    t0 = time.time()
    tot, rows = Counter(), []
    keys = [(P, Q) for Q in range(3, qmax + 1) for P in range(1, Q)
            if math.gcd(P, Q) == 1]
    keys.sort(key=lambda k: (k[1], k[0]))
    for (P, Q) in keys:
        if budget_s and time.time() - t0 > budget_s:
            print(f"  [budget stop at Q={Q}]", flush=True)
            break
        r = row(P, Q)
        rows.append(r)
        tot["rows"] += 1
        if r["status"] != "ok":
            tot[r["status"]] += 1
            if verbose:
                print(f"    {r['status']} {P}/{Q} eps={r['eps']}", flush=True)
            continue
        tot["scored"] += 1
        tot["cls_" + "".join(sorted(r["ver"]))] += 1
        for k in ("H1", "H2", "H3", "H4", "H6", "H7"):
            tot[k] += bool(r.get(k + "_ok"))
        for k, ok in r["negctl"].items():
            tot["neg_" + k] += bool(ok)
        for k in ("neg_H6_swap", "neg_H6_shift"):
            tot[k] += bool(r.get(k))
        tot["mvals"] = tot.get("mvals") or 0
        if not (r["H3_ok"] and r["H4_ok"] and r.get("H6_ok") and r.get("H7_ok")):
            print(f"    MISS {P}/{Q} eps={r['eps']} ver={r['ver']} s={r['s']} "
                  f"H4={r['H4_ok']} H6={r.get('H6_ok')} "
                  f"m=({r.get('m_in')},{r.get('m_out')}) pred={r.get('m_pred')} "
                  f"H7={r.get('H7_ok')} n={r['n_sigma']}", flush=True)
    if verbose:
        s = tot["scored"]
        print(f"\nSWEEP Q<={qmax}: {tot['rows']} rows, {s} scored, "
              f"{tot['PRONG_CAPPED']} PRONG_CAPPED, {tot['CELL_OPEN']} CELL_OPEN "
              f"[{round(time.time()-t0,1)}s]", flush=True)
        print(f"  H1 every canon occurs exactly twice   : {tot['H1']}/{s}", flush=True)
        print(f"  H3 s == 1  (overlap = ONE cylinder)   : {tot['H3']}/{s}", flush=True)
        print(f"  H4 the straddler is the orphan        : {tot['H4']}/{s}", flush=True)
        print(f"  H6 straddler POSITION == corner angles: {tot['H6']}/{s}", flush=True)
        print(f"  H7 non-straddlers occupy ADJACENT cells: {tot['H7']}/{s}", flush=True)
        print(f"  by class: " + "  ".join(f"{k[4:]}={v}" for k, v in sorted(tot.items())
                                          if k.startswith("cls_")), flush=True)
        print(f"  -- H6 negative controls (must be RARE): swapped pieces "
              f"{tot['neg_H6_swap']}/{s}, straddler shifted one cylinder "
              f"{tot['neg_H6_shift']}/{s}", flush=True)
        dist = Counter((r["m_in"], r["m_out"]) for r in rows if r.get("H6_ok"))
        print(f"  -- distinct (m_in,m_out) values realised: {len(dist)}", flush=True)
        print(f"  ⚠ H2 non-crossing {tot['H2']}/{s} -- NOT EVIDENCE, implied by H7;"
              f" its negative control confirms that:", flush=True)
        for k in WRONG_ORDERS:
            print(f"       {k:6s} non-crossing          : {tot['neg_' + k]}/{s}",
                  flush=True)
    json.dump({"summary": dict(tot), "qmax": qmax, "rows": rows},
              open(out_path, "w"), indent=1)
    print(f"wrote {out_path}", flush=True)
    return tot, rows


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "h5"
    kw = {}
    for a in sys.argv[2:]:
        if a.startswith("--qmax="):
            kw["qmax"] = int(a.split("=")[1])
        if a.startswith("--budget="):
            kw["budget_s"] = float(a.split("=")[1])
        if a.startswith("--out="):
            kw["out_path"] = a.split("=")[1]
    if mode == "h5":
        print("CONTROLS")
        h5_control(kw.get("qmax", 120))
    elif mode == "sweep":
        sweep(**kw)
    else:
        P, Q = (int(x) for x in mode.split("/"))
        print(json.dumps(row(P, Q), indent=1, default=str))


if __name__ == "__main__":
    main()
