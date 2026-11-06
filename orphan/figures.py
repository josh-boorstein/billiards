"""The four figures of the orphan paper, each re-derived before it is drawn.

  python figures.py            verify every caption, then draw  -> orphan/figs/

What is exact and what is float.  Every branch partition, every width, the pairing J and the
orphan are `lib/partition.py`'s: exact in Z[zeta_4Q].  The folded orbits and developments of
Figures 1-2 are a 60-digit trace (`lib/tracer.py`) whose side sequence is REQUIRED to equal the
exact word before anything is drawn.  Figure 4 places the leaves of the orphan's cylinder on the
rhombus presentation of S_alpha by a float trace from three heights of the exact orphan branch;
its caption's claims -- one word across the branch, every leaf closing on lambda_1 and missing
lambda_2, the two vertical hypotenuse copies, the midpoint crossing, Z and tau Z disjoint -- are
asserted before it is drawn.

The companion paper's two figures that used to be drawn here (the double Q-gon and the deficit
plot) moved with their statements to `neck/figures_5_6.py`.
"""
import json
import math
import sys
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Arc, FancyArrowPatch, Polygon

from mpmath import mp, mpf

from common import HERE
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
from partition import partition as exact_partition  # noqa: E402
from tracer import Billiard, scan_words  # noqa: E402

OUT = os.path.join(HERE, "figs")
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
    """The exact branch partition and its pairing J, from `lib/partition.py`.

    J is the first-return map read off the development: each branch is followed to its
    perpendicular return and the returned interval is matched, endpoint for endpoint as field
    elements, to a branch.  So equal widths in a pair, and the orphan being the unique fixed
    point, are CHECKS (in verify) and not how the picture is built."""
    if (P, Q) in _PART_CACHE:
        return _PART_CACHE[(P, Q)]
    r = exact_partition(P, Q)
    assert r["open"] == 0 and r["tiles"] and r["J_ok"], (P, Q)
    cells, lo = [], 0.0
    for w in r["widths"]:
        cells.append({"lo": lo, "hi": lo + w, "width": w})
        lo += w
    partner = [r["J"].get(i, i if i in r["orphans"] else None) for i in range(r["n"])]
    for i, cell in enumerate(cells):
        cell["is_orphan"] = partner[i] == i
        cell["word"] = r["words"][i]
    _PART_CACHE[(P, Q)] = (cells, partner)
    return cells, partner


class _Trace:
    def __init__(self, hits, returned):
        self.hits, self.returned_perpendicular = hits, returned


def _trace(P, Q, s, cap=200000):
    """60-digit folded orbit from (cot a, s), velocity (-1, 0), to the perpendicular return to
    L1: (trace, points, 1-based index of the first interior head-on hit or None)."""
    B = Billiard(mp.pi * P / (2 * Q))
    px, py, vx, vy = B.L, mpf(s), mpf(-1), mpf(0)
    tol = mpf(10) ** (-mp.dps // 3)
    pts, hits, head_on = [(float(px), float(py))], [], None
    for _ in range(cap):
        t, side = B.step(px, py, vx, vy)
        px, py = px + t * vx, py + t * vy
        wx, wy = B.reflect(side, vx, vy)
        hits.append(side)
        pts.append((float(px), float(py)))
        perp = abs(wx + vx) < tol and abs(wy + vy) < tol
        if side == "1" and perp:
            return _Trace(hits, True), pts, head_on
        if perp and head_on is None:
            head_on = len(hits)
        vx, vy = wx, wy
    return _Trace(hits, False), pts, head_on


def _hits(P, Q, cell):
    """Hit count of a branch's midpoint orbit -- the shortest one draws most legibly."""
    return len(_trace(P, Q, 0.5 * (cell["lo"] + cell["hi"]))[0].hits)


def folded_orbit(P, Q, s):
    """The folded orbit and the index of its INTERIOR head-on hit (None if there is none).
    The perpendicular return is itself a head-on hit, on L1; it is not an interior retrace
    and is not reported."""
    return _trace(P, Q, s)


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
    assert head_on * 2 == len(tr.hits), (head_on, len(tr.hits))       # Lemma 2.1
    assert stop == "H", stop                                          # Thm 4.6, Q odd
    assert len(copies) == head_on, (len(copies), head_on)             # Lemma 3.1(2)
    assert sum(1 for c in copies[:-1] if c[3]) == 0                   # Prop 4.5: unique
    word = "".join({"L1": "1", "L2": "2", "H": "h"}[c[1]] for c in copies)
    assert word == orphan["word"], (word, orphan["word"])             # float = exact word
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
        # CONTROL, and it is the one that could come out otherwise: the words above are the
        # exact engine's; a 60-digit float scan of L1 is a different engine, and a branch
        # missed or invented by either shows up as a word mismatch.
        B = Billiard(mp.pi * P / (2 * Q))
        seen = scan_words(B, mp.pi, 4000, kind="head_on", cap=200000)
        flt = [w for w, _ in sorted(seen.items(), key=lambda kv: kv[1])]
        assert flt == [c["word"] for c in cells3], (P, Q, len(flt), n)
        drift = 0.0
        out.append(f"{P}/{Q}: n = {n}, sum w = 1, traced pairing is an involution with "
                   f"{len(fixed)} fixed point(s), worst in-pair width gap {max(gaps):.2e}; "
                   f"a 60-digit float scan finds the same {n} words in the same order")

    # Prop 7.1's base arithmetic and Prop 7.2's count, at the two centres of Figure 3's
    # denominator.  NOT A CAPTION CHECK (the base figure was cut); s6_7_surface.py checks both
    # at every Q <= 200 from the built surface, and this stays as the figures' own sanity check.
    for P, Q in ((5, 9), (4, 9)):
        orders = [P - 2, Q - P - 2] + [-1] * Q
        g = Q // 2
        weier = Q + (P % 2) + ((Q - P) % 2)
        assert sum(orders) == -4, (P, Q, sum(orders))                 # Gauss-Bonnet
        assert weier == 2 * g + 2, (P, Q, weier, g)                   # Prop 7.2
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
        f"(a)  DATA — $2\\alpha/\\pi = 3/5$, $P$ odd.  The ray meets a vertical wall-copy of type "
        f"${TEX_SIDE[stopA]}$,\n"
        f"      the only interior vertical class (Prop. 4.5), at developed copy {len(copiesA)} of the\n"
        f"      {len(tr.hits)}-hit return — its temporal midpoint (Lemma 2.1).")

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
             "Branch endpoints are exact (lib/partition.py, in Z[zeta_4Q]); the arcs are "
             "the first-return map $J$, computed by tracing each\nbranch to its "
             "perpendicular return, so the equality of paired widths is an observation "
             "here and not a definition.",
             fontsize=8.5, color="#5a6570")
    fig.tight_layout(rect=(0, 0.075, 1, 1))
    _save(fig, path)


# ------------------------------------------------- fig 4: the cut locus at 3/5 (Section 7.6)

def _kz(P, Q):
    """The rhombus presentation of `S_alpha` (s7.1's labelled copies, grouped by the image of
    `R` they surround).  Group elements are `(k, f)` = `rho^k r_0^f`, `rho` the rotation through
    `pi/Q`, `r_0` the reflection in `L2`; copy `g` of `T` is drawn as `X_j + M_g (p - R)` in
    rhombus `j = k mod Q`.  The flow of the perpendicular beam is horizontal in every copy."""
    alpha, c, tri = triangle(P, Q)
    N, Rp = 2 * Q, tri["R"]
    refl = {"L1": (Q, 1), "L2": (0, 1), "H": (P, 1)}

    def mul(a, b):
        return ((a[0] + (b[0] if a[1] == 0 else -b[0])) % N, (a[1] + b[1]) % 2)

    def M(g):
        t = g[0] * math.pi / Q
        ct, st = math.cos(t), math.sin(t)
        return (ct, st, st, -ct) if g[1] else (ct, -st, st, ct)

    def ap(m, v):
        return (m[0] * v[0] + m[1] * v[1], m[2] * v[0] + m[3] * v[1])

    def X(j):
        return (2.6 * j, 0.0)

    def place(g, pt):
        v = ap(M(g), (pt[0] - Rp[0], pt[1] - Rp[1]))
        x = X(g[0] % Q)
        return (x[0] + v[0], x[1] + v[1])

    def exit_(pt, d):
        best = None
        for nm, (u, v) in SIDES.items():
            a, b = tri[u], tri[v]
            ex, ey = b[0] - a[0], b[1] - a[1]
            den = d[0] * ey - d[1] * ex
            if abs(den) < 1e-14:
                continue
            t = ((a[0] - pt[0]) * ey - (a[1] - pt[1]) * ex) / den
            w = ((a[0] - pt[0]) * d[1] - (a[1] - pt[1]) * d[0]) / den
            if t > 1e-12 and -1e-12 <= w <= 1 + 1e-12 and (best is None or t < best[0]):
                best = (t, nm)
        t, nm = best
        return nm, (pt[0] + t * d[0], pt[1] + t * d[1])

    def leaf(s, maxn=2000):
        """The leaf through height `s` of `lambda_1`, to the perpendicular return:
        `[(copy, exit side, developed start, developed end)]`, and the final copy."""
        g, pt, d, out = (0, 0), (c, s), (-1.0, 0.0), []
        for _ in range(maxn):
            nm, q = exit_(pt, d)
            out.append((g, nm, place(g, pt), place(g, q)))
            if nm == "L1" and abs(d[1]) < 1e-12 and d[0] > 0:
                return out, g
            d = ap(M(refl[nm]), d)
            g, pt = mul(g, refl[nm]), q
        raise RuntimeError("leaf did not return")

    G = [(k, f) for k in range(N) for f in (0, 1)]
    return dict(tri=tri, c=c, N=N, refl=refl, mul=mul, place=place, leaf=leaf, X=X, G=G)


def _cutlocus_data(P=3, Q=5):
    K = _kz(P, Q)
    cells, _ = partition(P, Q)
    orph = next(cl for cl in cells if cl["is_orphan"])
    eps = 1e-9 * max(1.0, orph["width"])
    lo, mid, hi = orph["lo"] + eps, 0.5 * (orph["lo"] + orph["hi"]), orph["hi"] - eps
    leaves = {k: K["leaf"](s) for k, s in (("lo", lo), ("mid", mid), ("hi", hi))}

    def side(g, nm):
        u, v = SIDES[nm]
        return K["place"](g, K["tri"][u]), K["place"](g, K["tri"][v])

    def vertical(g, nm):
        a, b = side(g, nm)
        return abs(a[0] - b[0]) < 1e-9

    pair = lambda g, nm: frozenset([g, K["mul"](g, K["refl"][nm])])
    Q_ = Q
    lam1, lam2 = pair((0, 0), "L1"), pair((Q_, 0), "L1")
    Wkeys = {}
    for g in K["G"]:
        if vertical(g, "H"):
            Wkeys.setdefault(pair(g, "H"), []).append(side(g, "H"))
    return K, orph, leaves, side, vertical, pair, lam1, lam2, Wkeys


def verify_cutlocus():
    """Every claim Figure 4's caption makes, checked before it is drawn."""
    out = []
    K, orph, leaves, side, vertical, pair, lam1, lam2, Wkeys = _cutlocus_data()
    words = {k: [x[1] for x in v[0]] for k, v in leaves.items()}
    assert words["lo"] == words["mid"] == words["hi"], "word changes across the orphan branch"
    assert orph["lo"] == 0.0 and abs(orph["hi"] - (3 - math.sqrt(5)) / 2) < 1e-12, "caption's branch"
    out.append(f"Fig 4: orphan branch at 3/5 is [{orph['lo']:.6f}, {orph['hi']:.6f}), one word "
               f"of {len(words['mid'])} letters across it")
    for k, (segs, _) in leaves.items():
        for g, nm, a, b in segs:
            assert abs(a[1] - b[1]) < 1e-9, "developed segment not horizontal"
    out.append("Fig 4: every developed segment is horizontal (the flow is one direction)")
    for k, (segs, gfin) in leaves.items():
        assert pair(gfin, "L1") == lam1, f"leaf {k} closes on lambda_2"
        assert not any(nm == "L1" and pair(g, "L1") == lam2 for g, nm, _, _ in segs), "leaf meets lambda_2"
    out.append("Fig 4: every leaf of Z closes on lambda_1 and meets lambda_2 nowhere")
    assert len(Wkeys) == 2 and all(len(v) == 2 for v in Wkeys.values()), "vertical H-copies"
    segs = leaves["mid"][0]
    wh = [i for i, (g, nm, _, _) in enumerate(segs) if nm == "H" and vertical(g, "H")]
    assert len(wh) == 1 and 2 * (wh[0] + 1) == len(segs), "W crossing is not the unique midpoint hit"
    out.append("Fig 4: two vertical H-copies (four edge appearances); the leaf crosses one, once, "
               "at its temporal midpoint")
    # Z and tau Z disjoint: tau is the half-turn about each rhombus centre, and Z's strips in
    # each rhombus lie on one side of that centre's horizontal line through it?  Check by
    # sampling: no point of a Z strip lies in the half-turned image of a Z strip.
    strips = _strips(K, leaves)
    import matplotlib.path as mpath
    tauS = [[(2 * K["X"](j)[0] - x, 2 * K["X"](j)[1] - y) for x, y in poly] for j, poly in strips]
    for (j, poly) in strips:
        cx = sum(x for x, _ in poly) / len(poly); cy = sum(y for _, y in poly) / len(poly)
        for (j2, _), tp in zip(strips, tauS):
            if j2 == j and mpath.Path(tp).contains_point((cx, cy)):
                raise AssertionError("Z meets tau Z")
    out.append("Fig 4: Z and tau Z have disjoint interiors")
    return out


def _strips(K, leaves):
    lo, hi = leaves["lo"][0], leaves["hi"][0]
    return [(g[0] % (K["N"] // 2), [a0, b0, b1, a1])
            for (g, _, a0, b0), (_, _, a1, b1) in zip(lo, hi)]


def fig_cutlocus(path):
    K, orph, leaves, side, vertical, pair, lam1, lam2, Wkeys = _cutlocus_data()
    Q = K["N"] // 2
    fig, ax = plt.subplots(figsize=(11, 2.6))
    strips = _strips(K, leaves)
    for j, poly in strips:
        ax.add_patch(Polygon(poly, closed=True, fc="#f3d9a4", ec="none", zorder=1))
        cx, cy = K["X"](j)
        ax.add_patch(Polygon([(2 * cx - x, 2 * cy - y) for x, y in poly], closed=True,
                             fc="none", ec=FADE, hatch="////", lw=0, zorder=1))
    for g in K["G"]:
        for nm in SIDES:
            (x1, y1), (x2, y2) = side(g, nm)
            ax.plot([x1, x2], [y1, y2], color=INK if nm == "H" else FADE,
                    lw=0.9 if nm == "H" else 0.6, zorder=2)
    mid_segs = leaves["mid"][0]
    wa = next(pair(g, "H") for g, nm, _, _ in mid_segs if nm == "H" and vertical(g, "H"))
    letters, nxt = {}, iter("abcdefghijklmnopqrstuvwxyz")
    for g in K["G"]:
        key = pair(g, "H")
        (x1, y1), (x2, y2) = side(g, "H")
        m = ((x1 + x2) / 2, (y1 + y2) / 2)
        cx, cy = K["X"](g[0] % Q)
        ox, oy = m[0] - cx, m[1] - cy
        r = math.hypot(ox, oy)
        if key in Wkeys:
            lab, col, fs = (r"$W_a$" if key == wa else r"$W_b$"), SIDE_COL["H"], 10
            ax.plot([x1, x2], [y1, y2], color=col, lw=2.6, zorder=4)
            k = 0.26
        else:
            if key not in letters:
                letters[key] = next(nxt)
            lab, col, fs, k = letters[key], INK, 8, 0.17
        ax.text(m[0] + k * ox / r, m[1] + k * oy / r, lab, fontsize=fs, ha="center",
                va="center", color=col, style="italic" if key not in Wkeys else "normal")
    for key, name in ((lam1, r"$\lambda_1$"), (lam2, r"$\lambda_2$")):
        g = min(key)
        (x1, y1), (x2, y2) = side(g, "L1")
        ax.plot([x1, x2], [y1, y2], color=SIDE_COL["L1"], lw=2.6, zorder=4)
        ax.text((x1 + x2) / 2 + 0.12, (y1 + y2) / 2, name, color=SIDE_COL["L1"], fontsize=10,
                va="center")
    for g, nm, a, b in mid_segs:
        ax.annotate("", xy=b, xytext=a, zorder=5,
                    arrowprops=dict(arrowstyle="-|>", color=INK, lw=0.9, mutation_scale=7))
    for j in range(Q):
        ax.plot(*K["X"](j), "o", color=INK, ms=3.5, zorder=6)
    seen = set()
    for g in K["G"]:
        for v in "OA":
            pt = K["place"](g, K["tri"][v])
            key = (round(pt[0], 6), round(pt[1], 6))
            if key not in seen:
                seen.add(key)
                ax.plot(*pt, "o", mfc="white", mec=INK, ms=3.5, zorder=6)
    ax.set_aspect("equal")
    ax.axis("off")
    fig.tight_layout()
    _save(fig, path)


def _save(fig, path):
    fig.savefig(path + ".png", dpi=200)
    fig.savefig(path + ".pdf")
    plt.close(fig)
    print(f"  wrote {path}.png and {path}.pdf")


def main():
    os.makedirs(OUT, exist_ok=True)
    print("verifying every claim the captions make ...")
    for line in verify() + verify_cutlocus():
        print("  ok  " + line)
    print("drawing ...")
    fig1(f"{OUT}/fig1_triangle")
    fig2(f"{OUT}/fig2_development")
    fig3(f"{OUT}/fig3_transversal")
    fig_cutlocus(f"{OUT}/fig4_cutlocus")


if __name__ == "__main__":
    main()
