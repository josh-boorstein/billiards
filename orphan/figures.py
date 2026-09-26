"""The five figures of the orphan paper, each re-derived before it is drawn.

  python figures.py            verify every caption, then draw  -> orphan/figs/
  python figures.py --full     also recompute Figure 5's counts from scratch (hours)

What is exact and what is float.  Every branch partition, every width, the pairing J and the
orphan are `lib/partition.py`'s: exact in Z[zeta_4Q].  Two things are drawn from floats and say
so: the folded orbits and developments (a 60-digit trace, `lib/tracer.py`, whose side sequence is
REQUIRED to equal the exact word before anything is drawn), and Figure 4's double Q-gon
(`lib/trans_surface.py`, a float polygon model cut on separatrices and flowed; its closed forms are
checked to 1e-9).  Figure 5 plots the partition record (`partition_record.json`), and its caption's
arithmetic -- parity, the even-P exception, the P = 1 and P = 3 laws, the coverage -- is checked
on every run; `--full` recomputes every count in it.
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
import trans_surface as TS  # noqa: E402
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
    assert head_on * 2 == len(tr.hits), (head_on, len(tr.hits))       # Lemma 4.5
    assert stop == "H", stop                                          # Thm 4.6, Q odd
    assert len(copies) == head_on, (len(copies), head_on)             # Lemma 3.1(2)
    assert sum(1 for c in copies[:-1] if c[3]) == 0                   # Prop 4.4: unique
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

    # s8.1-s8.3's base arithmetic.  ⚠ NO LONGER A CAPTION CHECK: the base figure was
    # CUT at s594 and its drawer with it, but the prose of s8.1, s8.2 and Theorem 8.3
    # still asserts every number below, so the check outlives the picture it came from.
    for P, Q in ((5, 9), (4, 9)):
        orders = [P - 2, Q - P - 2] + [-1] * Q
        g = Q // 2
        weier = Q + (P % 2) + ((Q - P) % 2)
        assert sum(orders) == -4, (P, Q, sum(orders))                 # Gauss-Bonnet
        assert weier == 2 * g + 2, (P, Q, weier, g)                   # Prop 8.2
        out.append(f"{P}/{Q} base Q({P-2}, {Q-P-2}, -1^{Q}): sum of orders = -4, "
                   f"g = {g}, Weierstrass points = {weier} = 2g+2")


    # Fig 5 -- the double `Q`-gon of s9.5.  Every closed form Theorem 9.8's proof states,
    # re-derived from the polygon model by separatrix cut and flow, at three odd `Q`.
    for Q in (5, 7, 11):
        D = qgon(Q)
        m, a = D["m"], D["a"]
        assert len(D["H"]) == m and len(D["V"]) == m, (Q, len(D["H"]), len(D["V"]))
        for i in range(1, m + 1):                                   # Step 1, Step 2
            h = 2 * math.sin(2 * math.pi * i / Q)
            assert abs(D["H"][i]["height"] - h) < 1e-9
            assert abs(D["H"][i]["circ"] - 4 / math.tan(math.pi / Q)
                       * math.sin(2 * math.pi * i / Q)) < 1e-8
            assert abs(D["V"][i]["height"] - h * math.tan(a)) < 1e-9
        mods = {round(D["H"][i]["height"] / D["H"][i]["circ"], 12) for i in D["H"]}
        assert mods == {round(math.tan(math.pi / Q) / 2, 12)}, mods   # one modulus
        hp = [D["V"][i]["height"] for i in D["V"]]
        assert len({round(z, 10) for z in hp}) == m                   # pairwise distinct
        assert abs(sum(hp) - 1.0) < 1e-9, sum(hp)                     # = |L1|
        # Steps 3 and 4: `S0` meets each perpendicular cylinder twice, `L1` exactly once,
        # the two halves exchanged by `sigma`, and the leg is partitioned.
        tot = 0.0
        for i, iv in D["cuts"].items():
            assert len(iv) == 2, (Q, i, iv)                           # Step 3
            pos = [z for z in iv if z[0] > -1e-12]
            assert len(pos) == 1, (Q, i, iv)                          # Step 4
            lo, hi_ = pos[0]
            assert abs((hi_ - lo) - D["V"][i]["height"]) < 1e-9
            assert any(abs(z[0] + hi_) < 1e-9 and abs(z[1] + lo) < 1e-9
                       for z in iv), (Q, i, iv)                       # sigma-mirror
            tot += hi_ - lo
        assert abs(tot - 1.0) < 1e-9, tot
        assert 2 * m == Q - 1                                         # n(2/Q) = Q - 1
        out.append(f"double {Q}-gon: {m} cylinders each direction, one modulus "
                   f"tan(pi/{Q})/2, perp heights distinct and summing to 1, "
                   f"S0 meets each perp cylinder twice and L1 once, n(2/{Q}) = {Q - 1}")

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
             "Branch endpoints are exact (lib/partition.py, in Z[zeta_4Q]); the arcs are "
             "the first-return map $J$, computed by tracing each\nbranch to its "
             "perpendicular return, so the equality of paired widths is an observation "
             "here and not a definition.",
             fontsize=8.5, color="#5a6570")
    fig.tight_layout(rect=(0, 0.075, 1, 1))
    _save(fig, path)


# ------------------------------------------------------------------------------ fig 4

# --------------------------------------------- Figure 4: the double `Q`-gon of s9.5

# The colour of a cylinder is its INDEX in the side-parallel decomposition, and the same
# index colours its image under `Phi` in the perpendicular one -- that pairing is the
# figure's whole point, so the two panels must never be coloured independently.
CYL_COL = ["#2f6f9f", "#c0392b", "#4a8f6d", "#c08a3e", "#7d5ba6",
           "#3d9aa1", "#a85a2f", "#6b7b8c"]


def qgon(Q):
    """Everything Figure 5 draws, and everything its caption asserts, at one odd `Q`.

    ⚠ DIRECTION INDEXING IS THE TRAP HERE (`s529_leg_surjectivity` trap (iv)): on
    `double_Qgon` `axis = 1` is the SIDE-parallel (horizontal) direction -- s9.5's Step 1,
    where every closed form is proved -- and `axis = 0` is the PERPENDICULAR (vertical)
    one.  Getting them the wrong way round draws a picture of the wrong theorem.

    ⚠ AND THE STRIPS ARE NOT THE CYLINDERS (`s520_polygon_staircase` trap (i)): that
    shortcut is the `2Q`-gon's.  Here the decomposition comes from `Surface.cylinders`,
    which cuts on separatrices and flows; nothing below assumes a strip is a cylinder.
    """
    S = TS.double_Qgon(Q)
    m = (Q - 1) // 2
    H = S.cylinders(1)                       # side-parallel  (s9.5 Step 1)
    V = S.cylinders(0)                       # perpendicular  (s9.5 Step 2)
    a = math.pi / (2 * Q)

    # index both decompositions by `i = 1..m` through the closed forms of Steps 1 and 2,
    # so the two panels share a colour exactly when `Phi` relates the two cylinders.
    def by_height(cyls, pred):
        out = {}
        for c in cyls:
            i = min(range(1, m + 1), key=lambda j: abs(c["height"] - pred(j)))
            assert abs(c["height"] - pred(i)) < 1e-9, (c["height"], pred(i))
            assert i not in out, f"two cylinders at index {i}"
            out[i] = c
        return out

    hi = by_height(H, lambda i: 2 * math.sin(2 * math.pi * i / Q))
    vi = by_height(V, lambda i: 2 * math.tan(a) * math.sin(2 * math.pi * i / Q))

    # `S0` is side 0 of `Pi`: horizontal, length 2, midpoint at `x = 0` an image of `R`.
    (x0, y0), (x1, y1) = S.side(0, 0)
    assert abs(y0 - y1) < 1e-12, "side 0 is not horizontal"
    assert abs(abs(x1 - x0) - 2.0) < 1e-9, abs(x1 - x0)
    assert abs(x0 + x1) < 1e-12, "side 0 is not centred on x = 0"

    # the partition of `S0` by vertical cylinder: a member of `V` in polygon 0 is an
    # x-interval, and `S0` spans `x in [-1, 1]`, so no flowing is needed here.
    cuts = {}
    for i, c in vi.items():
        iv = []
        for (pp, lo, hi_) in c["members"]:
            if pp != 0:
                continue
            lo_, hi__ = max(lo, -1.0), min(hi_, 1.0)
            if hi__ - lo_ > 1e-9:
                iv.append((lo_, hi__))
        cuts[i] = sorted(iv)
    return dict(S=S, Q=Q, m=m, a=a, H=hi, V=vi, S0=(x0, x1, y0), cuts=cuts)


def _draw_qgon_panel(ax, D, axis, title):
    """One polygon pair with the `axis`-direction cylinders shaded."""
    S, Q = D["S"], D["Q"]
    R = 1.0 / math.sin(math.pi / Q)
    dx = 2.35 * R                                     # display offset for `Pi'`
    cyls = D["H"] if axis == 1 else D["V"]
    for i, c in sorted(cyls.items()):
        col = CYL_COL[(i - 1) % len(CYL_COL)]
        for (pp, lo, hi_) in c["members"]:
            piece = TS._clip(S.polys[pp], axis, lo, hi_)
            if len(piece) < 3:
                continue
            off = dx if pp == 1 else 0.0
            ax.add_patch(Polygon([(x + off, y) for (x, y) in piece], closed=True,
                                 facecolor=col, edgecolor="none", alpha=0.42, zorder=1))
    for pp, name in ((0, "$\\Pi$"), (1, "$\\Pi\'$")):
        off = dx if pp == 1 else 0.0
        pts = [(x + off, y) for (x, y) in S.polys[pp]]
        ax.add_patch(Polygon(pts, closed=True, facecolor="none",
                             edgecolor=INK, lw=1.3, zorder=3))
        cx = sum(x for x, _ in pts) / len(pts)
        cy = sum(y for _, y in pts) / len(pts)
        ax.text(cx, cy, name, ha="center", va="center", fontsize=13,
                color=INK, zorder=4)
        # side numbers -- the gluing is side `r` of `Pi` to side `r` of `Pi'`.
        # ⚠ pushed OUT along the side midpoint's radius and set in ink, not grey: ref-3's
        # note on Figure 2 was that low-contrast numbers over the drawing do not read.
        for r in range(Q):
            if pp == 0 and r == 0:
                continue          # that is `S_0`, named in full below -- a digit there
            (ax0, ay0), (ax1, ay1) = S.side(pp, r)   # would land on the marked point `R`
            mx, my = (ax0 + ax1) / 2, (ay0 + ay1) / 2
            ax.text(mx * (1 + 0.13) + off, my * (1 + 0.13), str(r),
                    ha="center", va="center", fontsize=7.5, color=INK,
                    alpha=0.75, zorder=4)
    ax.set_title(title, fontsize=10.5, color=INK)
    ax.set_aspect("equal")
    ax.axis("off")
    return dx


def fig4(path):
    """The double regular `Q`-gon of s9.5 -- the model, both decompositions, and `S0`.

    ⇒ THIS IS FIGURE 4.  It REPLACED the genus-zero base picture at s594 (referee 3:
    that figure's "entire informational content is `3 is odd, 2 is even`", it was
    schematic in a paper whose other figures are not, and it hung on the result s592
    demoted, while s9.5 -- which carries Theorem 9.8 -- had no illustration at all).
    """
    Q = 7
    D = qgon(Q)
    S, m = D["S"], D["m"]
    x0, x1, y0 = D["S0"]
    R = 1.0 / math.sin(math.pi / Q)

    fig = plt.figure(figsize=(12.6, 5.6))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 0.34], hspace=0.05, wspace=0.04)
    axA, axB = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])
    axS = fig.add_subplot(gs[1, :])

    # ---- (a) the side-parallel direction, and `Gamma`
    dx = _draw_qgon_panel(axA, D, 1, "(a) DATA. The side-parallel direction (Step 1), "
                                     f"and $\\Gamma$ (Step 3)")
    for off in (0.0, dx):
        ys = [y for (_, y) in S.polys[0 if off == 0.0 else 1]]
        axA.plot([off, off], [min(ys), max(ys)], color=PERP, lw=1.6, zorder=5)
    axA.text(0.0, max(y for (_, y) in S.polys[0]) + 0.18, "$\\Gamma$", color=PERP,
             ha="center", fontsize=11)
    axA.plot([x0, x1], [y0, y0], color=INK, lw=3.2, solid_capstyle="butt", zorder=6)
    axA.plot([0.0], [y0], marker="o", ms=5.5, color=INK, zorder=7)
    axA.text(0.0, y0 - 0.78, "$R$", ha="center", fontsize=10, color=INK)
    axA.text(1.02, y0 - 0.34, "$\\lambda_1\\!=\\!L_1$", ha="left", fontsize=9, color=ACC)
    axA.text(-1.02, y0 - 0.34, "$\\lambda_2$", ha="right", fontsize=9, color=ACC)
    axA.text(0.0, y0 + 0.26, "$S_0 = \\Lambda$ (side 0)", ha="center", fontsize=9.5,
             color=INK)

    # ---- (b) the perpendicular direction
    _draw_qgon_panel(axB, D, 0, "(b) DATA. The perpendicular direction (Step 2), "
                                "the image of (a) under $\\Phi$")
    axB.plot([x0, x1], [y0, y0], color=INK, lw=3.2, solid_capstyle="butt", zorder=6)
    axB.plot([0.0], [y0], marker="o", ms=5.5, color=INK, zorder=7)
    axB.text(0.0, y0 - 0.62, "$S_0$ crosses each cylinder twice", ha="center",
             fontsize=9, color=INK)

    # ---- `S0` magnified: the partition of the leg, one interval per cylinder
    for i, iv in sorted(D["cuts"].items()):
        col = CYL_COL[(i - 1) % len(CYL_COL)]
        for (lo, hi_) in iv:
            axS.add_patch(plt.Rectangle((lo, 0.0), hi_ - lo, 1.0, facecolor=col,
                                        edgecolor="white", lw=1.2, alpha=0.62))
            axS.text((lo + hi_) / 2, 0.5, str(i), ha="center", va="center",
                     fontsize=9, color=INK)
    axS.plot([0, 0], [-0.18, 1.18], color=INK, lw=1.4)
    axS.plot([0.0], [0.0], marker="o", ms=6, color=INK, clip_on=False, zorder=6)
    axS.annotate("", xy=(1.0, -0.42), xytext=(0.0, -0.42),
                 arrowprops=dict(arrowstyle="<->", color=ACC, lw=1.2))
    axS.text(0.5, -0.78, "$\\lambda_1 = L_1$, length 1", ha="center", fontsize=9,
             color=ACC)
    axS.annotate("", xy=(0.0, -0.42), xytext=(-1.0, -0.42),
                 arrowprops=dict(arrowstyle="<->", color=FADE, lw=1.2))
    axS.text(-0.5, -0.78, "$\\lambda_2 = \\sigma(\\lambda_1)$", ha="center",
             fontsize=9, color=FADE)
    axS.text(0.06, 1.12, "$R$", ha="left", fontsize=10, color=INK)
    axS.set_title("(c) DATA. $S_0 = \\Lambda$ magnified: its partition by perpendicular "
                  "cylinder (Steps 3 and 4). Each index appears exactly once in "
                  "$\\lambda_1$ and once in $\\lambda_2$.", fontsize=10.5, color=INK)
    axS.set_xlim(-1.08, 1.08)
    axS.set_ylim(-0.95, 1.30)
    axS.axis("off")

    # ⚠ NOT `tight_layout`: it does not understand `set_aspect("equal")` axes and warns,
    # and the row gap it leaves is most of the figure.  Limits are set by hand instead.
    for ax in (axA, axB):
        ax.set_xlim(-R - 0.45, 2.35 * R + R + 0.45)
        ax.set_ylim(-R - 1.15, R + 0.75)
    fig.subplots_adjust(left=0.015, right=0.985, top=0.93, bottom=0.03)
    _save(fig, path)


# ----------------------------------------------------------------------------- driver

# --------------------------------------------------------------- figure 5: the deficit
#
# The counts come from the partition record (every centre of the window the census closed, on
# the exact engine), and `--full` recomputes each from scratch.  The table the figure was first
# drawn from -- 160 counts certified on the CosPoly tree -- is kept as a CONTROL: every one of
# them must equal the record's count.
LEGACY_COUNTS = {
    1: {2: 1, 3: 1, 4: 1, 5: 1, 7: 1, 8: 1, 11: 1, 16: 1, 25: 1},
    2: {3: 2, 5: 4, 7: 6, 9: 8, 11: 10, 13: 12, 15: 14, 17: 16, 19: 18, 21: 20, 23: 22,
        25: 24, 41: 40},
    3: {4: 3, 5: 3, 7: 3, 8: 5, 10: 5, 11: 7, 13: 7, 14: 9, 16: 9, 17: 11, 19: 11,
        20: 13, 22: 13, 25: 15, 26: 17, 28: 17, 29: 19, 35: 23},
    4: {5: 4, 9: 8, 11: 10, 13: 12, 15: 14, 17: 16, 19: 18, 23: 22, 27: 26, 31: 30,
        35: 34, 39: 38},
    5: {6: 5, 7: 3, 8: 3, 11: 7, 12: 5, 13: 9, 14: 7, 16: 5, 17: 5, 18: 11, 19: 11,
        21: 5, 22: 5, 23: 13, 24: 17, 26: 5, 27: 7, 28: 15, 29: 21, 31: 5, 32: 11,
        33: 17, 34: 25, 36: 5, 37: 13, 38: 19, 39: 29, 41: 5, 42: 15, 47: 17, 52: 19,
        57: 21},
    6: {7: 6, 13: 12, 17: 16, 19: 18, 23: 22, 25: 24, 29: 28, 31: 30, 35: 34, 37: 36,
        41: 40, 43: 42},
    7: {8: 7, 9: 5, 10: 5, 15: 11, 16: 3, 17: 9, 18: 13, 19: 11, 20: 15, 22: 17, 23: 3,
        24: 5, 25: 17, 26: 11, 27: 21, 29: 15, 30: 15, 31: 17, 32: 15, 33: 17, 34: 27,
        36: 23, 37: 19, 38: 7, 39: 15, 40: 23, 45: 7, 46: 17, 52: 7, 59: 13},
    8: {9: 8, 15: 6, 17: 16, 19: 18, 21: 20, 23: 22, 25: 24, 29: 28, 33: 32, 37: 36,
        39: 38, 45: 44},
    9: {10: 9, 11: 5, 19: 15, 20: 11, 22: 13, 23: 13, 25: 13, 26: 17, 28: 9, 29: 9,
        31: 13, 32: 13, 34: 19, 35: 21, 37: 15, 40: 17, 41: 13, 43: 23, 44: 27, 47: 17,
        49: 9, 58: 23},
}
QMAX5 = 60


def window_counts():
    from common import centre, record
    rec = record()
    flat = {}
    for P in range(1, 10):
        for Q in range(P + 1, QMAX5 + 1):
            if math.gcd(P, Q) != 1:
                continue
            r = centre(P, Q) if "--full" in sys.argv else rec.get(f"{P}/{Q}")
            if r is not None and r.get("n") is not None:
                flat[(P, Q)] = r["n"]
    return flat


def verify_counts_table():
    """Everything Figure 5's caption states, checked on the counts actually plotted."""
    out = []
    flat = window_counts()
    legacy = {(P, Q): n for P, row in LEGACY_COUNTS.items() for Q, n in row.items()}
    miss = [k for k in legacy if k not in flat]
    bad = [(k, legacy[k], flat[k]) for k in legacy if k in flat and flat[k] != legacy[k]]
    assert not bad, f"legacy CosPoly counts disagree with the exact engine: {bad}"
    out.append(f"figure 5: {len(legacy) - len(miss)} of the {len(legacy)} legacy CosPoly counts "
               f"are in the record and all agree" + (f" ({len(miss)} not yet recorded)"
                                                      if miss else ""))
    bad = [(P, Q) for (P, Q), n in flat.items() if n % 2 != P % 2]
    assert not bad, f"Corollary 5.1 (n = P mod 2) fails at {bad}"
    out.append(f"figure 5: n = P (mod 2) on all {len(flat)} plotted rows (Corollary 5.1)")
    bado = [(P, Q) for (P, Q), n in flat.items()
            if P % 2 == 1 and ((Q - 1 - n) - Q) % 2 != 0]
    assert not bado, f"d = Q (mod 2) fails at odd P: {bado}"
    out.append("figure 5: d = Q (mod 2) on every odd-P row (s10 item 2's congruence)")
    ex = sorted((P, Q, Q - 1 - n) for (P, Q), n in flat.items()
                if P % 2 == 0 and n != Q - 1)
    assert ex == [(8, 15, 8)], f"even-P exceptions: {ex}"
    out.append("figure 5: d = 0 at every even-P centre but 8/15, where d = 8")
    assert all(n == 1 for (P, Q), n in flat.items() if P == 1)
    b3 = [(Q, n) for (P, Q), n in flat.items()
          if P == 3 and Q >= 5 and n != 2 * round(Q / 3) - 1]
    assert not b3, f"Proposition 5.2 fails at {b3}"
    out.append("figure 5: n = 1 at P = 1; Proposition 5.2 holds on every P = 3 row")
    dmax = max(Q - 1 - n for (P, Q), n in flat.items() if P % 2 == 1 and P >= 5)
    out.append(f"figure 5: largest odd-P deficit in the window is d = {dmax}")
    cov, tot = {}, 0
    for P in range(1, 10):
        allq = [Q for Q in range(P + 1, QMAX5 + 1) if math.gcd(P, Q) == 1]
        tot += len(allq)
        cov[P] = round(100 * sum(1 for (p, _q) in flat if p == P) / len(allq))
    assert tot == 322, tot
    out.append(f"figure 5: coverage {len(flat)}/{tot} of the window; by P: "
               + ", ".join(f"{P}:{cov[P]}%" for P in range(1, 10)))
    return out, flat


def fig5(path):
    """DATA.  The deficit `d = Q-1-n` of every certified centre with `P <= 9`, `Q <= 60`.

    Deficit rather than the count itself, so that what is PROVED becomes the axis: every
    even-`P` centre but `8/15` has `d = 0`, and the whole plotting area is then the open
    odd-`P` range that s10 item 2 is about.
    """
    lines, flat = verify_counts_table()
    fig, ax = plt.subplots(figsize=(7.4, 4.5))
    cols = {5: "#c0392b", 7: "#7b4397", 9: "#1f7a5a"}

    ax.axhline(0, color=FADE, lw=1.2, zorder=1)
    ax.plot([2, QMAX5], [0, QMAX5 - 2], color=FADE, lw=1.0, ls=":", zorder=1)
    ax.text(QMAX5 - 1, QMAX5 - 4, r"$d=Q-2$  (i.e. $n=1$)", color="#8a949e",
            fontsize=8, ha="right", va="top")

    # the three families with a known count get their own light marks: lumping them
    # would put P = 1 on the d = Q-2 bound and P even on d = 0 in one colour, which
    # reads as one class behaving two ways.
    ev = sorted((Q, Q - 1 - n) for (P, Q), n in flat.items() if P % 2 == 0)
    ax.scatter(*zip(*ev), s=22, color=FADE, marker="s", lw=0, zorder=2,
               label=r"$P$ even (verified; $P=2$ is Theorem 9.8)")
    for P, mk, lab in ((1, "^", r"$P=1$ ($n=1$)"),
                       (3, "v", r"$P=3$ (Proposition 5.2)")):
        pts = sorted((Q, Q - 1 - n) for (p, Q), n in flat.items() if p == P)
        ax.plot(*zip(*pts), color="#8a949e", lw=.7, ls="-", alpha=.6, zorder=2)
        ax.scatter(*zip(*pts), s=22, color="#8a949e", marker=mk, lw=0, zorder=2,
                   label=lab)
    for P in (5, 7, 9):
        pts = [(Q, Q - 1 - n) for (p, Q), n in flat.items() if p == P]
        ax.scatter(*zip(*pts), s=30, color=cols[P], lw=0, alpha=.85, zorder=3,
                   label=rf"$P={P}$  (open)")

    ax.annotate(r"$8/15$", xy=(15, 8), xytext=(19, 26), fontsize=9, color=INK,
                arrowprops=dict(arrowstyle="-", color=INK, lw=.7, shrinkA=0, shrinkB=3))
    ax.text(QMAX5 - 1, 1.8, r"$d=0$: the count is maximal, $n=Q-1$",
            color="#8a949e", fontsize=8, ha="right")
    ax.set_xlabel(r"$Q$")
    ax.set_ylabel(r"deficit $\ d = Q-1-n(P/Q)$")
    ax.set_xlim(0, QMAX5 + 2)
    ax.set_ylim(-3, QMAX5 - 9)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.legend(frameon=False, fontsize=9, loc="upper left")
    fig.tight_layout()
    _save(fig, path)
    plt.close(fig)
    return lines


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
    fig1(f"{OUT}/fig1_triangle")
    fig2(f"{OUT}/fig2_development")
    fig3(f"{OUT}/fig3_transversal")
    fig4(f"{OUT}/fig4_qgon")
    for line in fig5(f"{OUT}/fig5_deficit"):
        print("  ok  " + line)


if __name__ == "__main__":
    main()
