#!/usr/bin/env python3
"""
s503_paper_figures.py -- the four figures of `paper_orphan.md` (Paper I, item E5).

NOT a probe and NOT a claim: every object drawn here is already stated and proved (or,
where it is not, already labelled as such) in the draft.  Two figures are DATA -- drawn
from the certified partition tree and from the developed ray -- and two are SCHEMATIC, and
each says which on its own face, because a cartoon that reads as a measurement is its own
trap (the `probes/s447_figure.py` convention, followed here).

⚠ IT READS NO `data/` DIRECTORY AND NO STORE, deliberately: it ships in the public deposit,
and [OPS-269] found that a `data/` read is in nobody's dependency graph, so a script that
takes one dies on a clean clone after doing all its work.  The partitions come from
`n_exact.n_pq_exact` with the depth escalated until the tree certifies complete.

PRIOR ART: grepped `rulings.py --grep adjacent`, `--grep figure`, and `orphan_theorem.md`
for the pairing structure; listed `probes/*figure*`/`*plot*`.  Came back: (i) [NCYL-186] --
"the cells come in adjacent equal-width pairs plus the orphan" is PROVED as far as the
EQUAL-WIDTH PAIRING goes (Thm D + Thm A / Lemma B'), but **POSITIONAL ADJACENCY of a pair
is OPEN** (measured 198/198 centres; proved only for pairs separated by an `R`-graze,
[NCYL-217]/[NCYL-218], and ⚠ [NCYL-223] corrects the claim that it is the same statement
as [NCYL-216]'s converse -- it is a second open question).  ⇒ Figure 3 therefore draws the
COMPUTED pairing at its two centres and its caption asserts only what Theorem 2.2 gives;
it must not be read as claiming adjacency, which is not among §10's open items and is not
being added to them.  (ii) `probes/s444_figure.py` and `probes/s447_figure.py` already draw
a genus-zero base -- but for the NECKLACE strand (the separatrix diagram of one cylinder
class; `B` as a necklace of `Q` kites).  Figure 4 here draws a different reading of the
same surface: the singularity data of `Q(P-2, Q-P-2, -1^Q)` marked by RAMIFICATION, which
is Theorem W's picture and is drawn nowhere.  (iii) No existing script draws the triangle,
the development, or the transversal partition.

  s503_fig1_triangle.png     (a) SCHEMATIC.  `T` with its vertices, its three labelled
                             sides and its angles, and the launch parametrisation of §2.1:
                             the beam leaves `(cot a, s)`, `s` in (0,1), with velocity
                             (-1, 0).  (b) DATA.  The orphan orbit at `3/5`, FOLDED, from
                             the float tracer -- it strikes `H` head-on at its temporal
                             midpoint and retraces (Lemma 4.5 / Theorem 4.6).
  s503_fig2_development.png  DATA.  The same orbit UNFOLDED, and its even-`P` twin.  The
                             ray is the horizontal line `y = s` and the triangle is
                             reflected across each side the ray exits (Lemma 3.1).
                             (a) `3/5`: the ray meets a vertical copy of `H` -- the UNIQUE
                             interior vertical wall-copy of the return period (Prop 4.4),
                             hence the head-on hit.  (b) `2/5`: no interior copy is
                             vertical at all, and the first vertical copy the ray meets is
                             a copy of `L1` -- the perpendicular return (Prop 4.3).
  s503_fig3_transversal.png  DATA.  The exact branch partition of the launch leg, from the
                             certified tree, with each `J`-pair joined by an arc.  The
                             pairing is computed by TRACING each branch's midpoint to its
                             perpendicular return -- never inferred from equal widths, so
                             that the equal widths are an observation and not a tautology
                             ([OPS-071]).  (a) `5/9`, `P` odd: `n = 7`, three pairs and the
                             orphan, which is its own partner.  (b) `2/9`, `P` even:
                             `n = 8`, four pairs, no orphan (Theorem 4.6).
  s503_fig4_base.png         SCHEMATIC.  The genus-zero base `Q(P-2, Q-P-2, -1^Q)` of
                             Prop 8.1, with its `Q` simple poles, and its singularities
                             marked by parity: ODD order = a branch point of the holonomy
                             double cover = a Weierstrass point of `S_a` (Prop 8.2).
                             (a) `5/9`: `P` odd, so `O` is a branch point -- Theorem W's
                             positive case.  (b) `4/9`: `P` even, so `O` has order 2, is
                             not a branch point, and has two preimages exchanged by `tau`.

Every panel's numbers are checked before it is drawn (`verify()`), and the checks are the
figures' point rather than decoration: the turnaround index must be half the return's hit
count (Lemma 4.5), the turnaround type must follow the `Q`-parity rule (Theorem 4.6), the
partition must have `sum w = 1` and the traced pairing must be an involution with exactly
`[P odd]` fixed points (Theorem 2.2 + Theorem 4.6).

Run: PYTHONPATH=.:engine:archive/scripts_2026-07:archive/scripts_2026-08:probes \\
     .venv/bin/python3.13 probes/s503_paper_figures.py
"""
import math
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Arc, FancyArrowPatch, Polygon

from n_exact import n_pq_exact
from right_triangle_billiards import RightTriangleBilliard

OUT = "probes/figs"
INK, ACC, PERP, FADE = "#1f2933", "#2f6f9f", "#c0392b", "#b8c2cc"
SIDE_COL = {"L1": "#2f6f9f", "L2": "#4a8f6d", "H": "#c08a3e"}
SIDES = {"L1": ("R", "A"), "L2": ("O", "R"), "H": ("O", "A")}


# --------------------------------------------------------------------- geometry helpers

def triangle(P, Q):
    """`T_0` of §2.1: O = (0,0), R = (cot a, 0), A = (cot a, 1)."""
    alpha = (P / Q) * (math.pi / 2)
    c = 1.0 / math.tan(alpha)
    return alpha, c, {"O": (0.0, 0.0), "R": (c, 0.0), "A": (c, 1.0)}


def develop(P, Q, s, max_copies=400):
    """Unfold the ray `y = s` (Lemma 3.1): reflect `T` across each side the ray exits.

    Returns (copies, stop), where `copies` is a list of
    `(vertices, exit_side, exit_x, is_vertical)` in order and `stop` names why the walk
    ended.  The walk ends at the FIRST vertical wall-copy, which by Lemma 3.1(2) is where
    the folded orbit meets a side head-on -- an interior side for an orphan, and the copy
    of `L1` that is the perpendicular return otherwise.
    """
    _alpha, c, tri = triangle(P, Q)
    x = c
    copies = []
    for _ in range(max_copies):
        best = None
        for name, (u, v) in SIDES.items():
            (x1, y1), (x2, y2) = tri[u], tri[v]
            if abs(y2 - y1) < 1e-12:
                continue
            t = (s - y1) / (y2 - y1)
            if not (-1e-9 <= t <= 1.0 + 1e-9):
                continue
            xc = x1 + t * (x2 - x1)
            if xc < x - 1e-9 and (best is None or xc > best[1]):
                best = (name, xc, (x1, y1), (x2, y2))
        if best is None:
            return copies, "no-exit"
        name, xc, p1, p2 = best
        vertical = abs(p2[0] - p1[0]) < 1e-9
        copies.append((dict(tri), name, xc, vertical))
        if vertical:
            return copies, name
        dx, dy = p2[0] - p1[0], p2[1] - p1[1]
        norm = math.hypot(dx, dy)
        dx, dy = dx / norm, dy / norm

        def reflect(pt, p1=p1, dx=dx, dy=dy):
            wx, wy = pt[0] - p1[0], pt[1] - p1[1]
            d = wx * dx + wy * dy
            return (p1[0] + 2 * d * dx - wx, p1[1] + 2 * d * dy - wy)

        tri = {k: reflect(v) for k, v in tri.items()}
        x = xc
    return copies, "maxed"


_PART_CACHE = {}


def partition(P, Q):
    """The exact branch partition, with the traced `J`-pairing.

    ⚠ The boundaries come from the CERTIFIED partition tree of Appendix A
    (`n_exact.n_pq_exact`, depth escalated until it reports zero truncations), NOT from a
    stored measurement -- so this runs from a clean checkout with no data directory, which
    is the defect [OPS-269] found the hard way: a `data/` read is in nobody's dependency
    graph.

    The pairing is read off the RETURN MAP, never off equal widths: each branch's midpoint
    is traced to its perpendicular return and the branch containing that height is the
    partner.  Equal width, and the orphan being the unique fixed point, are then CHECKS on
    Theorem 2.2 and Theorem 4.6 rather than the way the picture was built ([OPS-071]).
    """
    if (P, Q) in _PART_CACHE:
        return _PART_CACHE[(P, Q)]
    depth = 16 * Q
    while True:
        det = n_pq_exact(P, Q, max_depth=depth, return_detail=True)
        if det["complete"]:
            break
        depth *= 2
        if depth > 200000:
            raise RuntimeError(f"partition at {P}/{Q} did not certify complete")
    edges = [0.0] + list(det["boundaries"]) + [1.0]
    cells = [{"lo": a, "hi": b, "width": b - a} for a, b in zip(edges[:-1], edges[1:])]
    B = RightTriangleBilliard((P / Q) * (math.pi / 2))
    partner = []
    for cell in cells:
        mid = 0.5 * (cell["lo"] + cell["hi"])
        tr = B.trace(mid, max_hits=200000)
        if not tr.returned_perpendicular:
            partner.append(None)
            continue
        js = tr.l1_returns[-1].s
        partner.append(next((j for j, y in enumerate(cells)
                             if y["lo"] - 1e-9 <= js <= y["hi"] + 1e-9), None))
    for i, cell in enumerate(cells):
        cell["is_orphan"] = (partner[i] == i)
    _PART_CACHE[(P, Q)] = (cells, partner)
    return cells, partner


def _hits(P, Q, cell):
    """Hit count of a branch's midpoint orbit -- the shortest one draws most legibly."""
    B = RightTriangleBilliard((P / Q) * (math.pi / 2))
    return len(B.trace(0.5 * (cell["lo"] + cell["hi"]), max_hits=200000).hits)


def folded_orbit(P, Q, s):
    """The folded orbit and the index of its INTERIOR head-on hit (None if there is none).

    ⚠ The perpendicular return is itself a head-on hit, on `L1`: the ray meets the vertical
    `L1`-copy at right angles and the reflection reverses the velocity there too.  It is
    not an interior retrace, and Prop. 4.3's proof excludes it in exactly these words, so
    the final hit of a completed return is skipped here.
    """
    B = RightTriangleBilliard((P / Q) * (math.pi / 2))
    tr = B.trace(s, max_hits=200000)
    pts = [B.start_point(s)] + [(h.x, h.y) for h in tr.hits]
    interior = tr.hits[:-1] if tr.returned_perpendicular else tr.hits
    head_on = None
    for h in interior:
        vx, vy = math.cos(h.incoming_angle), math.sin(h.incoming_angle)
        rx, ry = B.reflect_velocity(h.side, vx, vy)
        if abs(rx + vx) < 1e-7 and abs(ry + vy) < 1e-7:
            head_on = h.index
            break
    return tr, pts, head_on


# ------------------------------------------------------------------------ verification

def verify():
    """Every claim a caption makes, checked before anything is drawn."""
    out = []

    # Fig 1b / Fig 2a -- the orphan at 3/5.
    cells, _ = partition(3, 5)
    orphan = next(c for c in cells if c["is_orphan"])
    s = 0.5 * (orphan["lo"] + orphan["hi"])
    tr, _pts, head_on = folded_orbit(3, 5, s)
    copies, stop = develop(3, 5, s)
    assert head_on is not None, "3/5 orphan has no head-on hit"
    assert head_on * 2 == len(tr.hits), (head_on, len(tr.hits))       # Lemma 4.5
    assert stop == "H", stop                                          # Thm 4.6, Q odd
    assert len(copies) == head_on, (len(copies), head_on)             # Lemma 3.1(2)
    assert sum(1 for c in copies[:-1] if c[3]) == 0                   # Prop 4.4: unique
    out.append(f"3/5 orphan: {len(tr.hits)} hits, head-on at {head_on} = half, "
               f"turnaround type {stop}, {len(copies)} developed copies, "
               f"{sum(1 for c in copies if c[3])} vertical")

    # Fig 2b -- an even-`P` branch at 2/5.
    cells2, _ = partition(2, 5)
    assert not any(c["is_orphan"] for c in cells2)                    # Prop 4.3
    cell = min(cells2, key=lambda z: _hits(2, 5, z))
    s2 = 0.5 * (cell["lo"] + cell["hi"])
    tr2, _p2, head_on2 = folded_orbit(2, 5, s2)
    copies2, stop2 = develop(2, 5, s2)
    assert head_on2 is None, head_on2                                 # Prop 4.3
    assert stop2 == "L1", stop2                                       # the return itself
    assert sum(1 for c in copies2[:-1] if c[3]) == 0                  # no interior vertical
    out.append(f"2/5 branch: {len(tr2.hits)} hits, no head-on, first vertical copy is "
               f"{stop2} at developed copy {len(copies2)}")

    # Fig 3 -- the two partitions and their pairings.
    for P, Q in ((5, 9), (2, 9)):
        cells3, partner = partition(P, Q)
        n = len(cells3)
        assert abs(sum(c["width"] for c in cells3) - 1.0) < 1e-9      # CHK, sum w = 1
        assert None not in partner, (P, Q, partner)
        assert all(partner[partner[i]] == i for i in range(n))        # Thm 2.2: J^2 = id
        fixed = [i for i in range(n) if partner[i] == i]
        assert len(fixed) == (P % 2), (P, Q, fixed)                   # Thm 4.6
        assert all(c["is_orphan"] == (partner[i] == i)
                   for i, c in enumerate(cells3))
        gaps = [abs(cells3[i]["width"] - cells3[partner[i]]["width"]) for i in range(n)]
        assert max(gaps) < 1e-9, max(gaps)                            # equal widths
        assert n % 2 == P % 2, (n, P)                                 # Thm 4.6's parity
        # CONTROL, and it is the one that could come out otherwise ([OPS-041]): the
        # boundaries above are the exact CosPoly tree's; the plain float tracer is a
        # different engine, and a branch missed by either shows up as a count mismatch.
        B = RightTriangleBilliard((P / Q) * (math.pi / 2))
        flt = B.discover_periodic_branches(n_samples=4000)
        assert len(flt) == n, (P, Q, len(flt), n)
        drift = max(abs(a[0] - c["lo"]) for a, c in zip(flt, cells3))
        assert drift < 1e-6, (P, Q, drift)
        out.append(f"{P}/{Q}: n = {n}, sum w = 1, traced pairing is an involution with "
                   f"{len(fixed)} fixed point(s), worst in-pair width gap {max(gaps):.2e}; "
                   f"float tracer agrees on n and on every boundary to {drift:.1e}")

    # Fig 4 -- the base's singularity data.
    for P, Q in ((5, 9), (4, 9)):
        orders = [P - 2, Q - P - 2] + [-1] * Q
        g = Q // 2
        weier = Q + (P % 2) + ((Q - P) % 2)
        assert sum(orders) == -4, (P, Q, sum(orders))                 # Gauss-Bonnet
        assert weier == 2 * g + 2, (P, Q, weier, g)                   # Prop 8.2
        out.append(f"{P}/{Q} base Q({P-2}, {Q-P-2}, -1^{Q}): sum of orders = -4, "
                   f"g = {g}, Weierstrass points = {weier} = 2g+2")

    return out


# ------------------------------------------------------------------------------ fig 1

def fig1(path):
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(11.6, 4.9))

    # ---- (a) the labelled triangle, SCHEMATIC (a drawing angle, not a centre)
    alpha = math.radians(40.0)
    c = 1.0 / math.tan(alpha)
    O, R, A = (0.0, 0.0), (c, 0.0), (c, 1.0)
    axL.add_patch(Polygon([O, R, A], closed=True, fc="#f4f6f8", ec="none"))
    for name, (u, v) in SIDES.items():
        p, q = {"O": O, "R": R, "A": A}[u], {"O": O, "R": R, "A": A}[v]
        axL.plot([p[0], q[0]], [p[1], q[1]], color=SIDE_COL[name], lw=2.6, zorder=3)
    axL.text(c * 0.5, -0.11, r"$L_2$", color=SIDE_COL["L2"], ha="center", fontsize=13)
    axL.text(c + 0.07, 0.36, r"$L_1$", color=SIDE_COL["L1"], va="center", fontsize=13)
    axL.text(c * 0.30, 0.36, r"$H$", color=SIDE_COL["H"], ha="center", fontsize=13)

    axL.add_patch(Arc(O, 0.62, 0.62, theta1=0, theta2=40.0, color=INK, lw=1.2))
    axL.text(0.40, 0.10, r"$\alpha$", color=INK, fontsize=12)
    axL.plot([c - 0.075, c - 0.075, c], [0, 0.075, 0.075], color=INK, lw=1.1)
    axL.add_patch(Arc(A, 0.52, 0.52, theta1=230.0, theta2=270.0, color=INK, lw=1.2))
    axL.text(c - 0.31, 0.80, r"$\frac{\pi}{2}-\alpha$", color=INK, fontsize=11)

    s0 = 0.62
    axL.plot([0, c], [s0, s0], color=FADE, lw=1.0, ls=(0, (4, 3)), zorder=2)
    axL.add_patch(FancyArrowPatch((c, s0), (c - 0.42, s0), color=PERP, lw=2.0,
                                  arrowstyle="-|>", mutation_scale=15, zorder=4))
    axL.plot([c], [s0], "o", color=PERP, ms=6, zorder=5)
    axL.text(c + 0.07, s0, r"$(\cot\alpha,\, s)$", color=PERP, va="center", fontsize=11)
    axL.text(c - 0.30, s0 - 0.07, r"$(-1,0)$", color=PERP, fontsize=10,
             ha="center", va="top")
    for pt, lab, dx, dy in ((O, "O", -0.10, -0.09), (R, "R", 0.05, -0.11),
                            (A, "A", 0.05, 0.04)):
        axL.plot([pt[0]], [pt[1]], "o", color=INK, ms=5, zorder=5)
        axL.text(pt[0] + dx, pt[1] + dy, lab, color=INK, fontsize=12, fontweight="bold")
    axL.annotate("", xy=(c + 0.30, 1.0), xytext=(c + 0.30, 0.0),
                 arrowprops=dict(arrowstyle="<->", color=INK, lw=1.0))
    axL.text(c + 0.36, 0.50, r"$s\in(0,1)$", color=INK, rotation=90,
             va="center", fontsize=10)
    axL.set_title("(a)  SCHEMATIC — the triangle and the launch of §2.1\n ",
                  fontsize=11, color=INK, loc="left")
    axL.set_xlim(-0.30, c + 0.72)
    axL.set_ylim(-0.24, 1.18)

    # ---- (b) the orphan orbit at 3/5, folded, DATA
    P, Q = 3, 5
    _a, c2, tri = triangle(P, Q)
    cells, _ = partition(P, Q)
    orphan = next(z for z in cells if z["is_orphan"])
    s = 0.5 * (orphan["lo"] + orphan["hi"])
    tr, pts, head_on = folded_orbit(P, Q, s)

    axR.add_patch(Polygon([tri["O"], tri["R"], tri["A"]], closed=True,
                          fc="#f4f6f8", ec="none"))
    for name, (u, v) in SIDES.items():
        p, q = tri[u], tri[v]
        axR.plot([p[0], q[0]], [p[1], q[1]], color=SIDE_COL[name], lw=2.4, zorder=3)
    half = pts[:head_on + 1]
    axR.plot([p[0] for p in half], [p[1] for p in half],
             color=PERP, lw=1.7, zorder=4)
    axR.plot([p[0] for p in pts[head_on:]], [p[1] for p in pts[head_on:]],
             color=ACC, lw=1.7, ls=(0, (5, 2.5)), zorder=4)
    hx, hy = pts[head_on]
    axR.plot([hx], [hy], "o", color=PERP, ms=9, zorder=6)
    axR.plot([c2], [s], "o", color=INK, ms=6, zorder=6)
    axR.annotate("head-on hit on $H$\n(the turnaround)", xy=(hx, hy),
                 xytext=(hx - 0.62, hy + 0.16), color=PERP, fontsize=10, ha="center",
                 arrowprops=dict(arrowstyle="->", color=PERP, lw=1.1))
    axR.text(c2 + 0.05, s, r"launch $=$ return", color=INK, va="center", fontsize=10)
    axR.set_title(f"(b)  DATA — the orphan branch at $2\\alpha/\\pi = {P}/{Q}$:\n"
                  f"{len(tr.hits)} hits; out (solid) and back (dashed) coincide",
                  fontsize=11, color=INK, loc="left")
    axR.set_xlim(-0.30, c2 + 0.92)
    axR.set_ylim(-0.16, 1.24)

    for ax in (axL, axR):
        ax.set_aspect("equal")
        ax.axis("off")
    fig.tight_layout()
    _save(fig, path)


# ------------------------------------------------------------------------------ fig 2

def _draw_development(ax, P, Q, s, label):
    copies, stop = develop(P, Q, s)
    xs = [c[2] for c in copies]
    x0 = 1.0 / math.tan((P / Q) * (math.pi / 2))
    for verts, _name, _xc, _vert in copies:
        ax.add_patch(Polygon([verts["O"], verts["R"], verts["A"]], closed=True,
                             fc="#f6f7f9", ec=FADE, lw=0.6, zorder=1))
    prev_x = x0
    for i, (verts, exit_name, xc, vert) in enumerate(copies):
        for name, (u, v) in SIDES.items():
            p, q = verts[u], verts[v]
            is_the_wall = vert and name == exit_name
            ax.plot([p[0], q[0]], [p[1], q[1]],
                    color=PERP if is_the_wall else SIDE_COL[name],
                    lw=3.6 if is_the_wall else 1.5,
                    zorder=5 if is_the_wall else 3)
        ax.text(0.5 * (prev_x + xc), s - 0.115, str(i + 1), ha="center", va="top",
                fontsize=8.5, color="#7a848e", zorder=8)
        prev_x = xc
    ax.plot([min(xs), x0], [s, s], color=INK, lw=1.6, zorder=6)
    ax.add_patch(FancyArrowPatch((x0, s), (x0 - 0.30, s), color=INK, lw=1.6,
                                 arrowstyle="-|>", mutation_scale=13, zorder=7))
    ax.plot([x0], [s], "o", color=INK, ms=5, zorder=8)
    ax.text(x0 + 0.06, s, r"$(\cot\alpha, s)$", fontsize=9, color=INK, va="center")
    xw = copies[-1][2]
    ax.plot([xw], [s], "o", color=PERP, ms=9, zorder=9)
    ax.plot([xw, xw + 0.13], [s + 0.13, s + 0.13], color=PERP, lw=1.2, zorder=9)
    ax.plot([xw + 0.13, xw + 0.13], [s, s + 0.13], color=PERP, lw=1.2, zorder=9)
    ax.set_title(label, fontsize=10.5, color=INK, loc="left")
    ax.set_aspect("equal")
    ax.set_xlim(min(xs) - 0.22, x0 + 0.62)
    ax.axis("off")
    return copies, stop


TEX_SIDE = {"L1": "L_1", "L2": "L_2", "H": "H"}


def fig2(path):
    fig, (axA, axB) = plt.subplots(2, 1, figsize=(12.4, 6.4))

    cells, _ = partition(3, 5)
    orphan = next(z for z in cells if z["is_orphan"])
    s = 0.5 * (orphan["lo"] + orphan["hi"])
    tr, _pts, _head_on = folded_orbit(3, 5, s)
    copiesA, stopA = develop(3, 5, s)
    _draw_development(
        axA, 3, 5, s,
        f"(a)  DATA — $2\\alpha/\\pi = 3/5$, $P$ odd.  The ray meets the UNIQUE interior "
        f"vertical wall-copy,\n"
        f"      of type ${TEX_SIDE[stopA]}$ (Prop. 4.4), at developed copy "
        f"{len(copiesA)} of the {len(tr.hits)}-hit return — its temporal midpoint "
        f"(Lemma 4.5).")

    cells2, _ = partition(2, 5)
    cell = min(cells2, key=lambda z: _hits(2, 5, z))
    s2 = 0.5 * (cell["lo"] + cell["hi"])
    copiesB, stopB = develop(2, 5, s2)
    _draw_development(
        axB, 2, 5, s2,
        f"(b)  DATA — $2\\alpha/\\pi = 2/5$, $P$ even.  NO interior copy is vertical "
        f"(Prop. 4.3); the first vertical copy the ray\n"
        f"      meets is a copy of ${TEX_SIDE[stopB]}$ — the perpendicular return itself, "
        f"at copy {len(copiesB)}.")

    handles = [plt.Line2D([], [], color=SIDE_COL[k], lw=2.2,
                          label={"L1": "$L_1$", "L2": "$L_2$", "H": "$H$"}[k])
               for k in ("L1", "L2", "H")]
    handles.append(plt.Line2D([], [], color=PERP, lw=3.6,
                              label="the vertical wall-copy"))
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False,
               fontsize=9.5, bbox_to_anchor=(0.5, 0.005))
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    _save(fig, path)


# ------------------------------------------------------------------------------ fig 3

def _draw_partition(ax, P, Q, label):
    cells, partner = partition(P, Q)
    n = len(cells)
    pool = ["#7f8c9a", "#4a8f6d", "#8f6da8", "#c08a3e", "#3f7f8f", "#a05a72",
            "#6d7fa8", "#8f8f4a", "#5a8fa0", "#a86d4a"]
    y0, h = 0.0, 0.30
    seen = {}
    for i, cell in enumerate(cells):
        j = partner[i]
        key = tuple(sorted((i, j)))
        col = seen.setdefault(key, pool[len(seen) % len(pool)])
        is_orphan = (j == i)
        ax.add_patch(plt.Rectangle((cell["lo"], y0), cell["width"], h,
                                   fc=PERP if is_orphan else col,
                                   ec="white", lw=1.0,
                                   alpha=0.95 if is_orphan else 0.55, zorder=2))
        mid = 0.5 * (cell["lo"] + cell["hi"])
        if cell["width"] > 0.035:
            ax.text(mid, y0 + h / 2, str(i + 1), ha="center", va="center",
                    color="white" if is_orphan else INK, fontsize=9, zorder=3)
    for i, cell in enumerate(cells):
        j = partner[i]
        if j <= i:
            continue
        a = 0.5 * (cells[i]["lo"] + cells[i]["hi"])
        b = 0.5 * (cells[j]["lo"] + cells[j]["hi"])
        r = abs(b - a)
        ax.add_patch(Arc((0.5 * (a + b), y0 + h), r, 0.52,
                         theta1=0, theta2=180, color=ACC, lw=1.4, zorder=4))
    orph = [i for i in range(n) if partner[i] == i]
    if orph:
        i = orph[0]
        mid = 0.5 * (cells[i]["lo"] + cells[i]["hi"])
        ax.annotate("the orphan:\n$J=\\mathrm{id}$", xy=(mid, y0 + h),
                    xytext=(mid, y0 + h + 0.40), color=PERP, fontsize=10,
                    ha="center", arrowprops=dict(arrowstyle="->", color=PERP, lw=1.2))
    ax.plot([0, 1], [y0, y0], color=INK, lw=1.2, zorder=5)
    for t in (0.0, 1.0):
        ax.plot([t, t], [y0 - 0.035, y0 + 0.035], color=INK, lw=1.2, zorder=5)
    ax.text(0.0, y0 - 0.13, "0", ha="center", color=INK, fontsize=10)
    ax.text(1.0, y0 - 0.13, "1", ha="center", color=INK, fontsize=10)
    ax.set_title(label, fontsize=11, color=INK, loc="left")
    ax.set_xlim(-0.045, 1.045)
    ax.set_ylim(y0 - 0.24, y0 + h + 0.78)
    ax.axis("off")
    return cells, partner


def fig3(path):
    fig, (axA, axB) = plt.subplots(2, 1, figsize=(11.4, 5.4))
    cA, pA = _draw_partition(
        axA, 5, 9,
        "(a)  DATA — $2\\alpha/\\pi = 5/9$, $P$ odd:  $n = 7$, three mirror pairs "
        "(arcs) of equal width and one orphan.")
    cB, pB = _draw_partition(
        axB, 2, 9,
        "(b)  DATA — $2\\alpha/\\pi = 2/9$, $P$ even:  $n = 8$, four mirror pairs "
        "and no orphan.")
    fig.text(0.012, 0.012,
             "Branch endpoints are exact (the partition tree of Appendix A); the arcs are "
             "the first-return map $J$, computed by tracing each\nbranch to its "
             "perpendicular return, so the equality of paired widths is an observation "
             "here and not a definition.",
             fontsize=8.5, color="#5a6570")
    fig.tight_layout(rect=(0, 0.075, 1, 1))
    _save(fig, path)


# ------------------------------------------------------------------------------ fig 4

def _draw_base(ax, P, Q, label):
    ax.add_patch(plt.Circle((0, 0), 1.0, fc="#f4f6f8", ec=FADE, lw=1.2, zorder=1))
    ordO, ordA = P - 2, Q - P - 2
    g = Q // 2

    def mark(x, y, order, name, angle_lbl, side):
        odd = (order % 2 != 0)
        ax.plot([x], [y], marker="*" if odd else "o", ms=21 if odd else 13,
                color=PERP if odd else INK,
                markerfacecolor=PERP if odd else "white",
                markeredgecolor=PERP if odd else INK, mew=1.8, zorder=6)
        note = "one preimage,\nfixed by $\\tau$" if odd else "two preimages,\nswapped by $\\tau$"
        ax.text(x + side * 0.16, y + 0.10,
                f"{name}: order ${order}$\ncone ${angle_lbl}\\pi$\n{note}",
                ha="left" if side > 0 else "right", va="top", fontsize=9,
                color=PERP if odd else INK, zorder=7)

    mark(-0.30, 0.66, ordO, "$O$", P, -1)
    mark(0.30, 0.66, ordA, "$A$", Q - P, +1)
    for k in range(Q):
        th = math.pi + math.pi * (k + 0.5) / Q
        x, y = 0.64 * math.cos(th), -0.28 + 0.36 * math.sin(th)
        ax.plot([x], [y], marker="*", ms=12, color=PERP, zorder=6)
    ax.text(0.0, -0.86, f"${Q}$ simple poles over $R$: order $-1$, all odd",
            ha="center", va="top", fontsize=9.5, color=PERP)
    weier = Q + (P % 2) + ((Q - P) % 2)
    ax.text(0.0, 1.18,
            f"$Q({ordO},\\ {ordA},\\ -1^{{{Q}}})$    "
            f"$\\sum \\mathrm{{orders}} = {ordO}+{ordA}-{Q} = -4$",
            ha="center", fontsize=10, color=INK)
    ax.text(0.0, -1.16,
            f"Weierstrass points $= {Q} + {P % 2} + {(Q - P) % 2} = {weier} "
            f"= 2g+2$,   $g = {g}$",
            ha="center", fontsize=9.5, color=INK)
    ax.set_title(label, fontsize=10.5, color=INK, loc="left")
    ax.set_xlim(-1.42, 1.42)
    ax.set_ylim(-1.40, 1.44)
    ax.set_aspect("equal")
    ax.axis("off")


def fig4(path):
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(12.2, 6.0))
    _draw_base(axA, 5, 9,
               "(a)  SCHEMATIC — $5/9$, $P$ odd:  $O$ has ODD order,\n"
               "so it is a branch point and a Weierstrass point (Thm. W).")
    _draw_base(axB, 4, 9,
               "(b)  SCHEMATIC — $4/9$, $P$ even:  $O$ has EVEN order,\n"
               "so it is not, and has two preimages exchanged by $\\tau$.")
    fig.text(0.012, 0.012,
             "★ = odd order = a branch point of the holonomy double cover = a "
             "Weierstrass point of $S_\\alpha$ (Prop. 8.2).   ○ = even order.\n"
             "Positions on the sphere are schematic; the orders, the count and the "
             "ramification are not.",
             fontsize=8.5, color="#5a6570")
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    _save(fig, path)


# ----------------------------------------------------------------------------- driver

def _save(fig, path):
    fig.savefig(path + ".png", dpi=200)
    fig.savefig(path + ".pdf")
    plt.close(fig)
    print(f"  wrote {path}.png and {path}.pdf")


def main():
    os.makedirs(OUT, exist_ok=True)
    print("verifying every claim the captions make ...")
    for line in verify():
        print("  ok  " + line)
    print("drawing ...")
    fig1(f"{OUT}/s503_fig1_triangle")
    fig2(f"{OUT}/s503_fig2_development")
    fig3(f"{OUT}/s503_fig3_transversal")
    fig4(f"{OUT}/s503_fig4_base")


if __name__ == "__main__":
    main()
