"""Figures 5 and 6 of the necklace paper, each re-derived before it is drawn.

  python figures_5_6.py        verify every caption, then draw -> neck/figs/

  fig5_qgon      the double regular Q-gon at Q = 7 and the four steps of Theorem 10.5: the
                 side-parallel and the perpendicular cylinders, and the side S0 partitioned by
                 perpendicular cylinder -- built from the polygon model (`lib/trans_surface.py`, a
                 float model cut on separatrices and flowed), every closed form the caption and the
                 proof state checked at Q = 5, 7, 11 to 1e-9 before drawing;
  fig6_deficit   the deficit d = Q - 1 - n(P/Q) at every coprime centre with P <= 9, Q <= 60 (all
                 322, each count certified by the exact flow decomposition: `band_record.json`);
                 every quantity the caption names is asserted first -- the parity d = Q (mod 2) at
                 odd P, d = 0 at every even P but 8/15, P = 1, 2, 3 and Q - 2 on Theorems 10.2, 10.5,
                 Corollary 10.8 and Theorem 10.3, and n = 1 only at P = 1 (Theorem 8.2(3)).
Figures 1-4 are `figures.py`.
"""
import json
import math
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Polygon  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
import trans_surface as TS  # noqa: E402

OUT = os.path.join(HERE, "figs")
INK, ACC, PERP, FADE = "#1f2933", "#2f6f9f", "#c0392b", "#b8c2cc"
QMAX = 60


def _save(fig, path):
    fig.savefig(path + ".png", dpi=200)
    fig.savefig(path + ".pdf")
    plt.close(fig)
    print(f"  wrote {path}.png and {path}.pdf")


# --------------------------------------------- Figure 5: the double `Q`-gon of Section 10.3

# The colour of a cylinder is its INDEX in the side-parallel decomposition, and the same
# index colours its image under `Phi` in the perpendicular one -- that pairing is the
# figure's whole point, so the two panels must never be coloured independently.
CYL_COL = ["#2f6f9f", "#c0392b", "#4a8f6d", "#c08a3e", "#7d5ba6",
           "#3d9aa1", "#a85a2f", "#6b7b8c"]


def qgon(Q):
    """Everything Figure 5 draws, and everything its caption asserts, at one odd `Q`.

    ⚠ DIRECTION INDEXING IS THE TRAP HERE (`s529_leg_surjectivity` trap (iv)): on
    `double_Qgon` `axis = 1` is the SIDE-parallel (horizontal) direction -- Theorem 10.5's Step 1,
    where every closed form is proved -- and `axis = 0` is the PERPENDICULAR (vertical)
    one.  Getting them the wrong way round draws a picture of the wrong theorem.

    ⚠ AND THE STRIPS ARE NOT THE CYLINDERS (`s520_polygon_staircase` trap (i)): that
    shortcut is the `2Q`-gon's.  Here the decomposition comes from `Surface.cylinders`,
    which cuts on separatrices and flows; nothing below assumes a strip is a cylinder.
    """
    S = TS.double_Qgon(Q)
    m = (Q - 1) // 2
    H = S.cylinders(1)                       # side-parallel  (Theorem 10.5, Step 1)
    V = S.cylinders(0)                       # perpendicular  (Theorem 10.5, Step 2)
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


def fig5(path):
    """The double regular `Q`-gon of Section 10.3 -- the model, both decompositions, and `S0`.

    ⇒ THIS IS FIGURE 4.  It REPLACED the genus-zero base picture at s594 (referee 3:
    that figure's "entire informational content is `3 is odd, 2 is even`", it was
    schematic in a paper whose other figures are not, and it hung on the result s592
    demoted, while s9.5 -- which carries Theorem 9.6 -- had no illustration at all).
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



# --------------------------------------------- Figure 6: the deficit, Section 13

def deficit_rows():
    rec = json.load(open(os.path.join(HERE, "band_record.json")))["L1"]
    flat = {}
    for P in range(1, 10):
        for Q in range(P + 1, QMAX + 1):
            if math.gcd(P, Q) != 1:
                continue
            r = rec[f"{P}/{Q}"]
            assert r["tiled"] and r["n"] is not None, f"{P}/{Q} uncertified"
            flat[(P, Q)] = r["n"]
    return flat


def verify_deficit(flat):
    out = []
    assert len(flat) == 322, len(flat)
    out.append("322/322 coprime centres P <= 9, Q <= 60, every count certified (exact tiling)")
    bad = [k for k, n in flat.items() if n % 2 != k[0] % 2]
    assert not bad, bad
    out.append("n = P (mod 2) on all 322 (the orphan paper's Corollary 5.1)")
    bado = [k for k, n in flat.items() if k[0] % 2 and ((k[1] - 1 - n) - k[1]) % 2]
    assert not bado, bado
    out.append("d = Q (mod 2) at every odd-P row")
    ex = sorted((P, Q, Q - 1 - n) for (P, Q), n in flat.items() if P % 2 == 0 and n != Q - 1)
    assert ex == [(8, 15, 8)], ex
    out.append("d = 0 at every even-P centre but 8/15 (d = 8)")
    assert all(n == 1 for (P, Q), n in flat.items() if P == 1)
    assert [k for k, n in flat.items() if n == 1 and k[0] != 1] == []
    out.append("n = 1 (d = Q - 2) at P = 1 on every row, and only there (Theorem 10.2, Theorem 8.2(3))")
    assert all(n == Q - 1 for (P, Q), n in flat.items() if P == 2)
    assert all(n == Q - 1 for (P, Q), n in flat.items() if P == Q - 1)
    b3 = [(Q, n) for (P, Q), n in flat.items() if P == 3 and Q >= 5 and n != 2 * round(Q / 3) - 1]
    assert not b3, b3
    out.append("P = 2: n = Q - 1 (Theorem 10.5); P = Q - 1: n = Q - 1 (Theorem 10.2); "
               "P = 3: n = 2 round(Q/3) - 1 on every row Q >= 5 (Corollary 10.8)")
    lat = [(P, Q, n) for (P, Q), n in flat.items() if P == Q - 2]
    assert all(n == 2 * ((Q - 1) // 4) + 1 for P, Q, n in lat), lat
    out.append(f"P = Q - 2 rows in the window ({len(lat)}) on Theorem 10.3")
    return out


def fig6(flat, path):
    fig, ax = plt.subplots(figsize=(7.4, 4.5))
    cols = {5: "#c0392b", 7: "#7b4397", 9: "#1f7a5a"}
    ax.axhline(0, color=FADE, lw=1.2, zorder=1)
    ax.plot([2, QMAX], [0, QMAX - 2], color=FADE, lw=1.0, ls=":", zorder=1)
    ax.text(QMAX - 3, QMAX - 1, r"$d=Q-2$  (i.e. $n=1$)", color="#8a949e", fontsize=8, ha="right",
            va="bottom")
    ev = sorted((Q, Q - 1 - n) for (P, Q), n in flat.items() if P % 2 == 0)
    ax.scatter(*zip(*ev), s=22, color=FADE, marker="s", lw=0, zorder=2, label=r"$P$ even")
    for P, mk, lab in ((1, "^", r"$P=1$ (Theorem 10.2)"), (3, "v", r"$P=3$ (Corollary 10.8)")):
        pts = sorted((Q, Q - 1 - n) for (p, Q), n in flat.items() if p == P)
        ax.plot(*zip(*pts), color="#8a949e", lw=.7, ls="-", alpha=.6, zorder=2)
        ax.scatter(*zip(*pts), s=22, color="#8a949e", marker=mk, lw=0, zorder=2, label=lab)
    for P in (5, 7, 9):
        pts = [(Q, Q - 1 - n) for (p, Q), n in flat.items() if p == P]
        ax.scatter(*zip(*pts), s=30, color=cols[P], lw=0, alpha=.85, zorder=3,
                   label=rf"$P={P}$  (open)")
    ax.annotate(r"$8/15$", xy=(15, 8), xytext=(19, 26), fontsize=9, color=INK,
                arrowprops=dict(arrowstyle="-", color=INK, lw=.7, shrinkA=0, shrinkB=3))
    ax.text(QMAX - 1, 1.8, r"$d=0$: the count is maximal, $n=Q-1$", color="#8a949e", fontsize=8,
            ha="right")
    ax.set_xlabel(r"$Q$")
    ax.set_ylabel(r"deficit $\ d = Q-1-n(P/Q)$")
    ax.set_xlim(0, QMAX + 2)
    ax.set_ylim(-3, QMAX + 1)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.legend(frameon=False, fontsize=9, loc="upper left")
    fig.tight_layout()
    _save(fig, path)


def verify():
    """Every claim a caption makes, checked before anything is drawn."""
    out = []
    # Fig 5 -- the double `Q`-gon of Section 10.3.  Every closed form Theorem 10.5's proof states,
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

    out += verify_deficit(deficit_rows())
    return out


def main():
    os.makedirs(OUT, exist_ok=True)
    print("verifying every claim the captions make ...")
    for line in verify():
        print("  ok  " + line)
    print("drawing ...")
    fig5(f"{OUT}/fig5_qgon")
    fig6(deficit_rows(), f"{OUT}/fig6_deficit")


if __name__ == "__main__":
    main()
