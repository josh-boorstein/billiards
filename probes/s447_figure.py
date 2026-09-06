#!/usr/bin/env python3
"""
s447_figure.py -- illustrations for [NCYL-285]/[NCYL-286]/[NCYL-287].  NOT a probe: it
draws objects those rulings already measured, off the same instruments, so nothing here
is a new claim.  Two of the four figures are DATA and two are SCHEMATIC, and each says
which on its own face -- a cartoon that reads as a measurement is its own trap.

  s447_overlap.png    DATA.  One panel per parity class (`{L1,H}`, `{L1,L2}`, `{L2,H}`).
                      The triangle, the cells of BOTH perpendicular strips coloured by
                      cylinder off `s439.strip_cells`, and the unique cylinder met by
                      both -- the orphan -- with its actual traced double-normal core.
                      [NCYL-285]'s conclusion, in the billiard picture.
  s447_chords.png     DATA.  The same cells read in the geometric order along
                      `Fix(iota)`: every cylinder an ADJACENT pair except one chord that
                      spans the junction of the two vertical arcs.  `m_in`/`m_out` are
                      counted off the picture and compared with the corner-angle
                      prediction ([NCYL-286]).
  s447_disk.png       SCHEMATIC.  (a) the disk `D` -- `Fix(iota)` as its boundary, the
                      three corners, the interior poles, and the rectangle decomposition.
                      (b) the step that closes the proof: two straddlers would leave the
                      region `U`, whose Gauss-Bonnet budget is spent entirely on four
                      right angles, so `U` carries no singularity -- and then a maximal
                      cylinder has a singularity-free boundary component ([NCYL-285]).
  s447_necklace.png   SCHEMATIC.  `B` as a NECKLACE OF `Q` KITES ([NCYL-287]): `Q` radial
                      `H`-edges, inner circle collapsing to `Z_O` and outer to `Z_A`, one
                      simple pole per kite, `iota` the reflection fixing one kite and the
                      opposite edge, `D` one half.

Run: PYTHONPATH=.:engine:archive/scripts_2026-07:archive/scripts_2026-08:probes \\
     .venv/bin/python3.13 probes/s447_figure.py
"""
import math

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Arc, Polygon, Wedge

from right_triangle_billiards import RightTriangleBilliard
import s437_oblique_cylinders as s437
import s439_exact_cells as s439
import s446_covering_identity as s446
import s447_overlap_disk as s447

INK, ACC, PERP, FADE = "#1f2933", "#2f6f9f", "#c0392b", "#b8c2cc"
POOL = ["#7f8c9a", "#4a8f6d", "#8f6da8", "#c08a3e", "#3f7f8f", "#a05a72",
        "#6d7fa8", "#8f8f4a", "#5a8fa0", "#a86d4a", "#6da88f", "#8a5a8f"]


# ------------------------------------------------------------------ shared data build

def build(P, Q):
    """The pairing, off exactly the instruments `s447_overlap_disk.row` uses."""
    eps = next(e for e in (0, 1) if len(s446.hv_sides(P, Q, e)[1]) == 2)
    hor, ver = s446.hv_sides(P, Q, eps)
    bmap, _r, _n, _o, _a = s439.boundary_map(P, Q, eps, cap=200_000,
                                             cap2=2_000_000, abort_on_open=True)
    geo = s437.geometry(P, Q)
    B = RightTriangleBilliard((P / Q) * (math.pi / 2))
    cells, us = {}, {}
    for side in ver:
        us[side] = s437.perp_index(side, P, Q)
        cells[side], _ = s439.strip_cells(B, geo, Q, side, us[side],
                                          bmap[(side, us[side])], 200_000)
    reading = s447.arc_reading(ver, P, Q)
    seq = []
    for side, rev in reading:
        lab = [str(c["canon"]) for c in cells[side]]
        seq += lab[::-1] if rev else lab
    strad = s447._straddlers(seq, len(cells[reading[0][0]]))
    return dict(P=P, Q=Q, eps=eps, hor=hor, ver=ver, cells=cells, us=us,
                reading=reading, seq=seq, strad=strad[0] if strad else None,
                geo=geo, B=B)


def double_normal(d, side, x0, cap=200_000):
    """The orphan's CORE: launch from `(side, x0)` along the inward normal and stop at
    the PERPENDICULAR bounce (`v` exactly reversed), which is its temporal midpoint.
    That half-orbit IS the double normal of [NCYL-072]/[NCYL-080] -- the full period is
    this path traversed twice, so drawing it once is what shows the object."""
    geo, B, Q = d["geo"], d["B"], d["Q"]
    _, _, base, tan = geo[side]
    px, py = base[0] + x0 * tan[0], base[1] + x0 * tan[1]
    th = d["us"][side] * math.pi / (2 * Q)
    vx, vy = math.cos(th), math.sin(th)
    pts = [(px, py)]
    for _ in range(cap):
        cand = B._candidate_times(px, py, vx, vy)
        if not cand:
            return pts, None
        t, hit = min(cand, key=lambda z: z[0])
        px, py = px + t * vx, py + t * vy
        pts.append((px, py))
        nx, ny = geo[hit][1]
        dot = vx * nx + vy * ny
        wx, wy = vx - 2 * dot * nx, vy - 2 * dot * ny
        if wx * vx + wy * vy < -1 + 1e-9:            # exact reversal = head-on
            return pts, hit
        vx, vy = wx, wy
    return pts, None


def colours(d):
    order = []
    for lab in d["seq"]:
        if lab not in order:
            order.append(lab)
    return {lab: (PERP if lab == d["strad"] else POOL[i % len(POOL)])
            for i, lab in enumerate(order)}


# ------------------------------------------------------------------ F1: the billiard

def fig_overlap(path="data/s447_overlap.png", rows=((3, 7), (3, 8), (4, 9))):
    fig, axes = plt.subplots(1, len(rows), figsize=(4.6 * len(rows), 5.0))
    for ax, (P, Q) in zip(axes, rows):
        d = build(P, Q)
        geo, col = d["geo"], colours(d)
        L = 1.0 / math.tan((P / Q) * (math.pi / 2))
        tri = [(0, 0), (L, 0), (L, 1)]
        ax.add_patch(Polygon(tri, closed=True, fc="#f7f9fb", ec=INK, lw=1.6, zorder=1))
        for nm, (x, y), dx, dy, ha in ((r"$Z_O$", (0, 0), -7, -6, "right"),
                                       (r"$R_*$", (L, 0), 9, -6, "left"),
                                       (r"$Z_A$", (L, 1), 9, 5, "left")):
            ax.plot([x], [y], "o", ms=5, color=INK, zorder=6)
            ax.annotate(nm, (x, y), textcoords="offset points", xytext=(dx, dy),
                        fontsize=11, color=INK, va="center", ha=ha)
        # the DOUBLE NORMAL of the unique shared cylinder
        st = d["strad"]
        side0 = d["reading"][0][0]
        c0 = next(c for c in d["cells"][side0] if str(c["canon"]) == st)
        pts, foot = double_normal(d, side0, 0.5 * (c0["lo"] + c0["hi"]))
        ax.plot([p[0] for p in pts], [p[1] for p in pts], "-", lw=1.3,
                color=PERP, alpha=.9, zorder=3)
        for pt, s_ in ((pts[0], side0), (pts[-1], foot)):
            ax.plot([pt[0]], [pt[1]], "o", ms=7, mfc="white", mec=PERP, mew=2.0,
                    zorder=7)
        # the cells of both perpendicular strips, drawn just outside their side
        off = 0.036 * max(L, 1.0)
        for side in d["ver"]:
            ell, nrm, base, tan = geo[side]
            for c in d["cells"][side]:
                xs = [base[0] + w * tan[0] - off * nrm[0] for w in (c["lo"], c["hi"])]
                ys = [base[1] + w * tan[1] - off * nrm[1] for w in (c["lo"], c["hi"])]
                lab = str(c["canon"])
                ax.plot(xs, ys, "-", lw=6.0 if lab == st else 4.0,
                        color=col[lab], solid_capstyle="butt", zorder=5)
        legend = "   ".join(rf"$\perp${s}: $n={len(d['cells'][s])}$" for s in d["ver"])
        ax.annotate(legend + rf"    $\sum n={len(d['seq'])}=2C$", (.5, -.02),
                    xycoords="axes fraction", fontsize=10, color=ACC, ha="center")
        ax.annotate("red = the orphan's DOUBLE NORMAL,\nperpendicular at both open circles",
                    (.5, .045), xycoords="axes fraction", fontsize=8.5,
                    color=PERP, ha="center", va="bottom")
        ax.set_title(rf"${P}/{Q}$  class $\varepsilon={d['eps']}$   "
                     rf"verticals $\{{{','.join(d['ver'])}\}}$" "\n"
                     rf"$C={len(col)}$,  overlap $=\{{$orphan$\}}$ (red)",
                     fontsize=10.5, color=INK)
        pad = 0.14 * max(L, 1.0)
        ax.set_xlim(-pad - 0.16, L + pad + 0.10)
        ax.set_ylim(-pad - 0.30, 1 + pad + 0.10)
        ax.set_aspect("equal")
        ax.axis("off")
    fig.suptitle("DATA — the two perpendicular strips, cells coloured by cylinder off "
                 "`s439.strip_cells`; red = the one cylinder met by BOTH, with its "
                 "traced double-normal core  [NCYL-285]", fontsize=10, color=INK)
    fig.tight_layout(rect=(0, 0, 1, .93))
    fig.savefig(path, dpi=170)
    print("wrote", path)


# ------------------------------------------------------------------ F2: chord diagram

def fig_chords(path="data/s447_chords.png", rows=((3, 19), (7, 20), (4, 17))):
    fig, axes = plt.subplots(len(rows), 1, figsize=(11.4, 2.55 * len(rows)))
    for ax, (P, Q) in zip(axes, rows):
        d = build(P, Q)
        col = colours(d)
        seq, st = d["seq"], d["strad"]
        n1 = len(d["cells"][d["reading"][0][0]])
        gap = 1.3                                    # visual break at the junction
        xs = [i + (gap if i >= n1 else 0) for i in range(len(seq))]
        ax.plot([xs[0] - .5, xs[n1 - 1] + .35], [0, 0], "-", lw=3.0, color=ACC)
        ax.plot([xs[n1] - .35, xs[-1] + .5], [0, 0], "-", lw=3.0, color=ACC)
        pos = {}
        for i, lab in enumerate(seq):
            pos.setdefault(lab, []).append(i)
            ax.plot([xs[i]], [0], "|", ms=13, color=INK, zorder=4)
        HMAX = 2.0          # tall chords are FLATTENED so they stay in frame
        for lab, (i, j) in pos.items():
            r = (xs[j] - xs[i]) / 2.0
            ax.add_patch(Arc(((xs[i] + xs[j]) / 2.0, 0), 2 * r, 2 * min(r, HMAX),
                             theta1=0, theta2=180, lw=2.8 if lab == st else 1.5,
                             color=col[lab], zorder=3))
        mid = (xs[n1 - 1] + xs[n1]) / 2.0
        ax.plot([mid, mid], [-.6, 2.6], ":", lw=1.4, color=FADE, zorder=1)
        i, j = pos[st]
        m_in, m_out = (j - i - 1) / 2, (len(seq) - (j - i + 1)) / 2
        mi, mo = s447.m_predicted(d["ver"], P, Q)
        ax.annotate("junction of the two\nvertical arcs", (mid, -.55), fontsize=8,
                    color=FADE, ha="center", va="top")
        ax.annotate(rf"$m_{{in}}={m_in:g}$   (corner-angle prediction {mi:g})",
                    (mid, 2.75), fontsize=9.5, color=INK, ha="center", va="bottom")
        ax.annotate(rf"$m_{{out}}={m_out:g}$   (corner-angle prediction {mo:g})",
                    (xs[0] - .4, 2.75), fontsize=9.5, color=INK, ha="left",
                    va="bottom")
        ax.set_title(rf"${P}/{Q}$   $\{{{','.join(d['ver'])}\}}$   "
                     rf"$n=({','.join(str(len(d['cells'][s])) for s in d['ver'])})$, "
                     rf"$C={len(col)}$ — every cylinder an ADJACENT pair except the one "
                     rf"red chord spanning the junction", fontsize=9.5, color=INK)
        ax.set_xlim(xs[0] - 1.2, xs[-1] + 1.2)
        ax.set_ylim(-1.5, 3.5)
        ax.axis("off")
    fig.suptitle("DATA — the same cells read in the geometric order along Fix(iota); "
                 "the straddler's POSITION is forced by the corner angles  [NCYL-286]",
                 fontsize=10, color=INK)
    fig.tight_layout(rect=(0, 0, 1, .95))
    fig.savefig(path, dpi=170)
    print("wrote", path)


# ------------------------------------------------------------------ F3: the disk

def _pt(a, r=1.0):
    t = math.radians(a)
    return r * math.cos(t), r * math.sin(t)


def _arc(a1, a2, r=1.0, n=60):
    return [_pt(a1 + (a2 - a1) * i / n, r) for i in range(n + 1)]


def _disk_frame(ax):
    """The circle `Fix(iota)`, cut at the three vertex-copies: `L1` (top) and `H`
    (right) are the VERTICAL arcs, `L2` (left) is the horizontal one -- a leaf.
    Boundary angles: `Z_A` 30, `R_*` 150, `Z_O` 270."""
    ax.add_patch(plt.Circle((0, 0), 1.0, fc="#f7f9fb", ec="none", zorder=0))
    for a, b, c, lab, la in ((30, 150, ACC, r"$L_1$ (vertical arc)", 88),
                             (270, 390, ACC, r"$H$ (vertical arc)", 328),
                             (150, 270, INK, r"$L_2$ (horizontal — a leaf)", 200)):
        ax.add_patch(Arc((0, 0), 2, 2, theta1=a, theta2=b, lw=3.6, color=c, zorder=6))
        ax.annotate(lab, _pt(la, 1.33), fontsize=10, color=c, ha="center", va="center")
    for nm, a, dx, dy in ((r"$Z_A$", 30, 16, 7), (r"$R_*$", 150, -16, 7),
                          (r"$Z_O$", 270, 18, -7)):
        x, y = _pt(a)
        ax.plot([x], [y], "o", ms=8, color=INK, zorder=8)
        ax.annotate(nm, (x, y), textcoords="offset points", xytext=(dx, dy),
                    fontsize=12, color=INK, ha="center", va="center")
    ax.set_xlim(-1.72, 1.72)
    ax.set_aspect("equal")
    ax.axis("off")


def _straddler(ax, a1, a2, b1, b2, alpha=.24):
    """The rectangle whose vertical sides are the arcs `[a1,a2]` (on `L1`) and
    `[b1,b2]` (on `H`); its horizontal sides are the chords `a2->b1` and `b2->a1`."""
    ax.add_patch(Polygon(_arc(a1, a2) + _arc(b1, b2), closed=True, fc=PERP,
                         alpha=alpha, ec=PERP, lw=2.0, zorder=3))


def _wrap(ax, c1, c3, depth=.34, c=POOL[1]):
    """A pole-wrap cylinder: both feet on the SAME arc and ADJACENT, sharing the point
    where that pole's prong crosses.  Drawn with the pole and its prong."""
    ax.add_patch(Polygon(_arc(c1, c3) + _arc(c3, c1, r=1 - depth), closed=True, fc=c,
                         alpha=.16, ec=c, lw=1.6, zorder=2))
    cm = (c1 + c3) / 2.0
    px, py = _pt(cm, 1 - depth * .55)
    ax.plot([px, _pt(cm)[0]], [py, _pt(cm)[1]], "-", lw=1.4, color=c, zorder=4)
    ax.plot([px], [py], "o", ms=7, mfc="white", mec=c, mew=2.0, zorder=5)


def fig_disk(path="data/s447_disk.png"):
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(13.6, 8.0))

    # ---- (a) the decomposition ----
    _disk_frame(axL)
    _straddler(axL, 76, 94, 318, 336)
    axL.annotate("the ONE straddler = the orphan" "\n" "(one foot on each vertical arc)",
                 (-.14, .06), fontsize=9.5, color=PERP, ha="center")
    for c1, c3 in ((100, 122), (126, 146), (280, 300), (296, 314)):
        _wrap(axL, c1, c3)
    axL.annotate("every OTHER cylinder wraps exactly ONE interior pole:" "\n"
                 "its two feet are ADJACENT cells of one arc, separated" "\n"
                 "by the point where that pole's prong crosses",
                 (0, -1.46), fontsize=9.5, color=POOL[1], ha="center", va="center")
    axL.set_ylim(-1.78, 1.44)
    axL.set_title(r"(a) SCHEMATIC — the disk $D = B\setminus\mathrm{Fix}(\iota)$."
                  "\n"
                  r"Both zeros sit ON the boundary, so every singularity of "
                  r"$\mathrm{int}\,D$ is a SIMPLE POLE.", fontsize=10.5, color=INK)

    # ---- (b) the contradiction ----
    _disk_frame(axR)
    _straddler(axR, 52, 68, 344, 360, alpha=.30)
    _straddler(axR, 96, 112, 300, 316, alpha=.30)
    axR.add_patch(Polygon(_arc(68, 96) + _arc(316, 344), closed=True, fc="#c9a227",
                          alpha=.34, ec="#8a6d10", lw=1.8, zorder=2))
    for a in (68, 96, 316, 344):
        x, y = _pt(a)
        axR.plot([x], [y], "s", ms=7, mfc="white", mec=INK, mew=1.5, zorder=9)
    axR.annotate("$U$", (.60, .00), fontsize=17, color="#8a6d10", ha="center")
    axR.annotate(r"$Rect_1$", (.74, .34), fontsize=10, color=PERP, ha="center",
                 rotation=-64)
    axR.annotate(r"$Rect_2$", (.45, -.46), fontsize=10, color=PERP, ha="center",
                 rotation=-64)
    axR.annotate(r"the four white corners have interior angle $\pi/2$" "\n"
                 "(a vertical arc meets a horizontal cylinder boundary orthogonally)"
                 "\n"
                 r"$\Rightarrow$ Gauss–Bonnet spends its whole $2\pi$ on them" "\n"
                 r"$\Rightarrow$ every remaining term is $\geq 0$, so all vanish" "\n"
                 r"$\Rightarrow$ $U$ carries NO singularity" "\n"
                 r"$\Rightarrow$ $\mathrm{top}(Rect_1)\cup\iota\,\mathrm{top}(Rect_1)$"
                 r" is a closed leaf with no" "\n"
                 "      singularity on it — which no MAXIMAL cylinder has." "\n"
                 r"$\Rightarrow$ CONTRADICTION, so $s \leq 1$.",
                 (0, -1.52), fontsize=9.5, color=INK, ha="center", va="center")
    axR.set_ylim(-2.02, 1.44)
    axR.set_title("(b) SCHEMATIC — why two straddlers cannot both exist:" "\n"
                  "[NCYL-285] step (4), the region $U$ between them", fontsize=10.5,
                  color=INK)
    fig.tight_layout(rect=(0, 0, 1, .97))
    fig.savefig(path, dpi=170)
    print("wrote", path)


# ------------------------------------------------------------------ F4: the necklace

def fig_necklace(path="data/s447_necklace.png", Q=7, P=3):
    fig, ax = plt.subplots(figsize=(10.6, 7.8))
    ri, ro = 0.46, 1.20
    ax.add_patch(plt.Circle((0, 0), ro, fc="#f7f9fb", ec=INK, lw=2.0, zorder=1))
    ax.add_patch(plt.Circle((0, 0), ri, fc="white", ec=INK, lw=2.0, zorder=2))
    for k in range(Q):
        a = 2 * math.pi * k / Q + math.pi / 2
        ax.plot([ri * math.cos(a), ro * math.cos(a)],
                [ri * math.sin(a), ro * math.sin(a)], "-", lw=2.2,
                color=INK, zorder=3)                                  # H-edge
        am = a + math.pi / Q
        rp = 0.5 * (ri + ro)
        px, py = rp * math.cos(am), rp * math.sin(am)
        ax.plot([px], [py], "o", ms=8, mfc="white", mec=PERP, mew=2.0, zorder=5)
        ax.plot([ri * math.cos(am), px], [ri * math.sin(am), py], "-", lw=1.5,
                color=POOL[0], zorder=4)                              # L2-edge
        ax.plot([px, ro * math.cos(am)], [py, ro * math.sin(am)], "-", lw=1.5,
                color=POOL[3], zorder=4)                              # L1-edge
    axis = math.pi / 2 + math.pi / Q
    ax.plot([-1.52 * math.cos(axis), 1.52 * math.cos(axis)],
            [-1.52 * math.sin(axis), 1.52 * math.sin(axis)], "--", lw=1.7,
            color=ACC, zorder=6)
    ax.annotate(r"axis of $\iota$ — fixes one kite" "\n" r"and the opposite $H$-edge",
                (-0.86, 1.42), fontsize=9.5, color=ACC, ha="center")
    ax.annotate(r"$D$ = one half of the necklace" "\n"
                r"= a CHAIN of $(Q\!-\!1)/2$ kites" "\n"
                r"plus the half-kite on the axis",
                (1.96, -0.70), fontsize=9.5, color=ACC, ha="center")
    ax.annotate(r"$Z_O$", (0, 0), fontsize=15, color=INK, ha="center", va="center")
    ax.annotate(r"$Z_A$", (0, ro + 0.10), fontsize=15, color=INK, ha="center")
    ax.annotate(r"the INNER circle collapses to the single point $Z_O$ "
                r"($Q$ kite corners of angle $2\alpha$, total $P\pi$);" "\n"
                r"the OUTER circle collapses to $Z_A$ "
                r"($Q$ corners of angle $\pi-2\alpha$, total $(Q\!-\!P)\pi$)",
                (0, -1.50), fontsize=9.5, color=INK, ha="center", va="center")
    for lab, c, y in ((r"$H$-edge  $Z_O\!-\!Z_A$  (the kite gluings)", INK, 1.44),
                      (r"$L_2$-edge  $Z_O\!-\!R_i$", POOL[0], 1.30),
                      (r"$L_1$-edge  $R_i\!-\!Z_A$", POOL[3], 1.16),
                      (r"$R_i$ — one simple pole per kite, cone angle $\pi$",
                       PERP, 1.02)):
        ax.annotate(lab, (1.60, y), fontsize=9, color=c, ha="center")
    ax.set_xlim(-2.05, 2.95)
    ax.set_ylim(-1.78, 1.62)
    ax.set_aspect("equal")
    ax.axis("off")
    fig.suptitle(rf"SCHEMATIC — $B$ is a NECKLACE OF $Q$ KITES ($Q={Q}$ drawn):  "
                 rf"$V=Q+2$, $E=3Q$, $F=2Q$, $\chi=2$" "\n"
                 r"a kite = two triangle copies glued along BOTH legs = the rhombus "
                 r"about $R$ modulo rotation by $\pi$   [NCYL-287]",
                 fontsize=11, color=INK)
    fig.tight_layout(rect=(0, 0, 1, .92))
    fig.savefig(path, dpi=170)
    print("wrote", path)


if __name__ == "__main__":
    fig_overlap()
    fig_chords()
    fig_disk()
    fig_necklace()
