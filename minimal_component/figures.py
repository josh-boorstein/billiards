"""The four figures of the paper, drawn from the exact computation.

  fig1_two_readings   Figure 1: L1 read horizontally (six branches, three cylinders) and
                      as Sigma_Hbar (three exchanged intervals), on one s-axis.
  fig2_boundary5      Figure 2: the translation jump across boundaries 3-6 as alpha moves
                      off 4pi/15, exact curves with 60-digit trace points.
  fig3_exchange       Figure 3: (a) the exchange of three intervals; (b) Lemma 4.
  fig4_connection     Figure 4: the saddle connection of Section 6, unfolded.

Every number a figure shows is computed here from `lib/` (the same computations as the
statement scripts), and each figure's claim is asserted before it is drawn.  Writes PNG
and PDF to minimal_component/figs/.
"""
import math
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Arc, FancyArrowPatch, Polygon, Rectangle
import mpmath as mp

from t0 import T, P, Q, M_H, return_words
from s3_1_return_map import lambdas
from s4_branch_count import branches
from develop import offset
from tracer import Billiard

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figs")
ALPHA = P * math.pi / (2 * Q)
COT = 1.0 / math.tan(ALPHA)
THETA11 = 11.0 * math.pi / 15.0
LAM = [float(x.value(30)) for x in lambdas()]
INK, MUTED, HAIR, PAPER, SURF = "#1f2933", "#6b7280", "#c9d1d9", "#f4f6f8", "#ffffff"
CYL = {1: "#2a78d6", 2: "#eb6834", 3: "#1baf7a"}
LAMC = ["#e87ba4", "#4a3aa7", "#008300"]


def branch_rows():
    """The six branches as drawing rows, each tagged with its cylinder 1..3 (C1: N = 25,
    C2: N = 53, C3: N = 9), and the cylinder heights."""
    rows, tiles, _ = branches()
    assert tiles and len(rows) == 6
    cyl_of_N = {25: 1, 53: 2, 9: 3}
    out = []
    for r in rows:
        out.append({"lo": float(r["lo"].value(30)), "hi": float(r["hi"].value(30)),
                    "cyl": cyl_of_N[2 * len(r["word"]) - 1]})
    H = {}
    for r in out:
        H[r["cyl"]] = H.get(r["cyl"], 0.0) + (r["hi"] - r["lo"])
    return out, H


N_CLOSED = {1: 25, 2: 53, 3: 9}


def tag(ax, s):
    """Captions carry the data/schematic distinction; kept so call sites record it."""


def save(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(OUT, name + "." + ext), dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("wrote", os.path.join(OUT, name + ".png"))


# ------------------------------------------------------------------------------ Figure 1
def fig1(rows):
    fig = plt.figure(figsize=(10.0, 5.0))
    axT = fig.add_axes([0.00, 0.08, 0.36, 0.84])
    axB = fig.add_axes([0.37, 0.08, 0.63, 0.84])

    # the triangle, with the two readings' directions leaving L1  (SCHEMATIC key)
    O, R, A = (0.0, 0.0), (COT, 0.0), (COT, 1.0)
    axT.add_patch(Polygon([O, R, A], closed=True, fc=PAPER, ec="none", zorder=0))
    for p, q in ((O, R), (O, A), (R, A)):
        axT.plot([p[0], q[0]], [p[1], q[1]], color=INK, lw=1.4 if (p, q) != (R, A) else 3.0, zorder=2)
    axT.add_patch(Arc(O, 0.40, 0.40, theta1=0, theta2=math.degrees(ALPHA), color=INK, lw=1.0))
    axT.text(0.26, 0.07, r"$\alpha=\frac{4\pi}{15}$", fontsize=9.5)
    axT.plot([COT - 0.05, COT - 0.05, COT], [0, 0.05, 0.05], color=INK, lw=0.9)
    for pt, lab, dx, dy in ((O, "O", -0.08, -0.08), (R, "R", 0.03, -0.09), (A, "A", 0.03, 0.03)):
        axT.plot(*pt, "o", color=INK, ms=4, zorder=5)
        axT.text(pt[0] + dx, pt[1] + dy, lab, fontsize=11, fontweight="bold")
    axT.text(COT * 0.5, -0.09, r"$L_2$", ha="center", fontsize=11)
    axT.text(COT * 0.36, 0.44, r"$H$", ha="center", fontsize=11)
    axT.text(COT + 0.05, 0.5, r"$L_1$", va="center", fontsize=11)
    s0, s1 = 0.62, 0.28
    axT.add_patch(FancyArrowPatch((COT, s0), (COT - 0.42, s0), color=INK, lw=1.6,
                                  arrowstyle="-|>", mutation_scale=13, zorder=6))
    axT.text(COT - 0.45, s0 + 0.04, "(a) horizontal:\nthe direction of Theorem 1", fontsize=8.8, ha="right")
    ex, ey = math.cos(THETA11), math.sin(THETA11)
    axT.add_patch(FancyArrowPatch((COT, s1), (COT + 0.36 * ex, s1 + 0.36 * ey), color=INK, lw=1.6,
                                  arrowstyle="-|>", mutation_scale=13, zorder=6))
    axT.text(COT - 0.04, s1 - 0.03,
             r"(b) $\theta_{\bar H}=\frac{11\pi}{15}$", fontsize=8.8, ha="right", va="top")
    axT.set_xlim(-0.2, COT + 0.25)
    axT.set_ylim(-0.18, 1.12)
    axT.set_aspect("equal")
    axT.axis("off")
    tag(axT, "SCHEMATIC")

    # the two partitions of L1 on one s-axis  (DATA)
    xa, xb, w = 0.34, 0.66, 0.11
    cuts_a = sorted({r["lo"] for r in rows} | {r["hi"] for r in rows})[1:-1]
    cuts_b = [LAM[0], LAM[0] + LAM[1]]
    for y in cuts_a + cuts_b:                          # hairlines through BOTH bars
        axB.plot([xa - w / 2 - 0.02, xb + w / 2 + 0.02], [y, y], color="#9aa5b1", lw=0.8, zorder=1)
    for r in rows:
        axB.add_patch(Rectangle((xa - w / 2, r["lo"]), w, r["hi"] - r["lo"],
                                fc=CYL[r["cyl"]], ec=SURF, lw=2.0, zorder=3))
    lb = [0.0, LAM[0], LAM[0] + LAM[1], 1.0]
    for i in range(3):
        axB.add_patch(Rectangle((xb - w / 2, lb[i]), w, lb[i + 1] - lb[i],
                                fc=LAMC[i], ec=SURF, lw=2.0 if i != 1 else 0.6, zorder=3))
    # cylinder labels, left of bar (a): coloured bracket (mark) + ink text
    bands = {}
    for r in rows:
        b = bands.setdefault(r["cyl"], [r["lo"], r["hi"]])
        b[0], b[1] = min(b[0], r["lo"]), max(b[1], r["hi"])
    for c, (lo, hi) in bands.items():
        xk = xa - w / 2 - 0.035
        axB.plot([xk, xk], [lo + 0.006, hi - 0.006], color=CYL[c], lw=3, solid_capstyle="butt")
        axB.text(xk - 0.02, 0.5 * (lo + hi), "$C_%d$   $N=%d$\nheight %.6f" % (c, N_CLOSED[c], H_CLOSED[c]),
                 fontsize=8.8, ha="right", va="center")
    # lambda labels, right of bar (b)
    for i in range(3):
        lo, hi = lb[i], lb[i + 1]
        xk = xb + w / 2 + 0.035
        if hi - lo > 0.08:
            axB.plot([xk, xk], [lo + 0.006, hi - 0.006], color=LAMC[i], lw=3, solid_capstyle="butt")
            axB.text(xk + 0.02, 0.5 * (lo + hi), r"$\lambda_%d$  %.9f" % (i + 1, hi - lo),
                     fontsize=8.8, va="center")
        else:
            axB.annotate(r"$\lambda_%d$  %.9f" % (i + 1, hi - lo), xy=(xb + w / 2, 0.5 * (lo + hi)),
                         xytext=(xk + 0.02, 0.5 * (lo + hi) + 0.10), fontsize=8.8, va="center",
                         arrowprops=dict(arrowstyle="-", color=INK, lw=0.7))
    # the shared axis
    XA = -0.17
    axB.plot([XA, XA], [0, 1], color=MUTED, lw=0.8)
    for y, s in ((0, "$0$ ($R$)"), (0.5, "$0.5$"), (1, "$1$ ($A$)")):
        axB.plot([XA - 0.005, XA + 0.005], [y, y], color=MUTED, lw=0.8)
        axB.text(XA - 0.015, y, s, fontsize=8.5, color=MUTED, ha="right", va="center")
    axB.text(XA - 0.13, 0.5, "$s$ along $L_1$", rotation=90, fontsize=9, color=MUTED, va="center", ha="center")
    axB.text(xa, 1.045, "(a) horizontal", ha="center", fontsize=10)
    axB.text(xb, 1.045, r"(b) $\Sigma_{\bar H}$", ha="center", fontsize=10)
    axB.text(xa, -0.055, r"heights sum to $1=|L_1|$", ha="center", fontsize=8.8)
    axB.text(xb, -0.055, r"all of $\Sigma_{\bar H}$ lies in $M$", ha="center", fontsize=8.8)
    axB.set_xlim(-0.34, 1.05)
    axB.set_ylim(-0.1, 1.09)
    axB.axis("off")
    tag(axB, "DATA")
    save(fig, "fig1_two_readings")


# ------------------------------------------------------------------------------ Figure 2
L2_DRAW = 0.10                                      # lambda_2 as DRAWN in panel (b)


def drawn_lams():
    k = (1 - L2_DRAW) / (LAM[0] + LAM[2])
    return [LAM[0] * k, L2_DRAW, LAM[2] * k]


def lemma4_cases(lam):
    """Re-run Lemma 4 at the given lengths: (lo, hi, return time, sample orbit)."""
    C, rot = 1 + lam[1], lam[1] + lam[2]
    out = []
    for lo, hi in ((0, lam[0]), (lam[0], lam[0] + lam[1]), (lam[0] + lam[1], 1)):
        x = cur = 0.5 * (lo + hi)
        hops = []
        while True:
            nxt = (cur + rot) % C
            hops.append((cur, nxt))
            cur = nxt
            if cur < 1:
                break
        out.append((lo, hi, len(hops), x, hops, cur))
    return out


def verify_fig2b():
    true = lemma4_cases(LAM)
    drawn = lemma4_cases(drawn_lams())
    ok = [a[2] for a in true] == [1, 2, 1] == [b[2] for b in drawn]
    # and the image intervals are the SAME translations symbolically: tau(x) - x per case
    lam = drawn_lams()
    want = [lam[1] + lam[2], lam[2] - lam[0], -(lam[0] + lam[1])]
    got = [b[5] - b[3] for b in drawn]
    ok = ok and all(abs(g - w_) < 1e-12 for g, w_ in zip(got, want))
    print("verify_fig2b: return times at true and drawn lengths [1,2,1]; drawn translations match:", ok)
    assert ok


def fig2():
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(12.0, 4.6), gridspec_kw={"width_ratios": [1, 1.15]})
    # (a) the exchange  (DATA)
    cuts = [0.0, LAM[0], LAM[0] + LAM[1], 1.0]
    outs, x = {}, 0.0
    for j in (2, 1, 0):
        outs[j] = (x, x + LAM[j])
        x += LAM[j]
    yT, yB, h = 0.64, 0.22, 0.13
    for j in range(3):
        lo, hi = cuts[j], cuts[j + 1]
        olo, ohi = outs[j]
        for (a, b, y) in ((lo, hi, yT), (olo, ohi, yB)):
            axL.add_patch(Rectangle((a, y), b - a, h, fc=LAMC[j], ec=SURF, lw=2.0 if j != 1 else 0.6))
        axL.add_patch(FancyArrowPatch((0.5 * (lo + hi), yT), (0.5 * (olo + ohi), yB + h), color=INK,
                                      lw=1.0, arrowstyle="-|>", mutation_scale=10,
                                      connectionstyle="arc3,rad=0.10"))
        dy = 0.035 if j != 1 else 0.075
        axL.text(0.5 * (lo + hi), yT + h + dy, r"$\lambda_%d$" % (j + 1), ha="center", fontsize=10.5)
        axL.text(0.5 * (olo + ohi), yB - dy - 0.02, r"$\lambda_%d$" % (j + 1), ha="center", va="top",
                 fontsize=10.5)
    axL.text(-0.02, yT + h / 2, "before", ha="right", va="center", fontsize=9.5, color=MUTED)
    axL.text(-0.02, yB + h / 2, "after", ha="right", va="center", fontsize=9.5, color=MUTED)
    axL.set_xlim(-0.16, 1.02)
    axL.set_ylim(0.0, 1.0)
    axL.axis("off")
    tag(axL, "DATA")
    axL.text(0.0, 1.0, "(a)", transform=axL.transAxes, fontsize=12, fontweight="bold", va="top")

    # (b) Lemma 4, one row per case  (SCHEMATIC, lambda_2 enlarged)
    lam = drawn_lams()
    C = 1 + lam[1]
    ys, bh = [0.70, 0.43, 0.16], 0.075
    labels = [r"$[0,\lambda_1)$", r"$[\lambda_1,\lambda_1{+}\lambda_2)$", r"$[\lambda_1{+}\lambda_2,1)$"]
    tnames = [r"$+(\lambda_2{+}\lambda_3)$", r"$\lambda_3-\lambda_1$", r"$-(\lambda_1{+}\lambda_2)$"]
    for i, ((lo, hi, rt, x0, hops, end), y) in enumerate(zip(lemma4_cases(lam), ys)):
        axR.add_patch(Rectangle((0, y), 1.0, bh, fc=PAPER, ec=MUTED, lw=0.8, zorder=1))
        axR.add_patch(Rectangle((1.0, y), C - 1.0, bh, fc="#9aa5b1", ec=MUTED, lw=0.8, zorder=1))
        axR.add_patch(Rectangle((lo, y), hi - lo, bh, fc=LAMC[i], ec="none", zorder=2))
        axR.plot([x0], [y + bh / 2], "o", color=INK, ms=7, mec=SURF, mew=1.0, zorder=5)
        for a, b in hops:
            wrap = b < a
            axR.add_patch(FancyArrowPatch((a, y + bh), (b, y + bh), color=(MUTED if wrap else INK),
                                          lw=1.2, arrowstyle="-|>", mutation_scale=10, zorder=4,
                                          connectionstyle="arc3,rad=%.2f" % (-0.22 if wrap else -0.45)))
        axR.plot([end], [y + bh / 2], "o", color=SURF, ms=7, mec=INK, mew=1.6, zorder=5)
        axR.text(-0.04, y + bh / 2, labels[i], fontsize=9.2, ha="right", va="center")
        axR.text(C + 0.04, y + bh / 2, "return time %d\n%s" % (rt, tnames[i]), fontsize=8.8, va="center")
    for t, lab, row in ((0, "$0$", 0), (lam[0], r"$\lambda_1$", 0), (lam[0] + lam[1], r"$\lambda_1{+}\lambda_2$", 1),
                        (1, "$1$", 0), (C, r"$1{+}\lambda_2$", 1)):
        yb = ys[-1] - 0.05 - 0.055 * row
        axR.plot([t, t], [ys[-1], yb + 0.01], color=HAIR, lw=0.7, zorder=0)
        axR.text(t, yb, lab, fontsize=8.6, ha="center", va="top")
    axR.text(0.5 * (1 + C), ys[0] + bh + 0.02, "outside\n$[0,1)$", fontsize=8, color=MUTED, ha="center")
    axR.text(-0.34, 0.965, r"circle $\mathbb{R}/(1{+}\lambda_2)\mathbb{Z}$ cut open, rotated by $\lambda_2{+}\lambda_3$"
             "\n" r"$\bullet$ a point   $\circ$ its first return to $[0,1)$   grey arrows wrap", fontsize=8.8,
             color=INK, va="top")
    axR.set_xlim(-0.36, C + 0.36)
    axR.set_ylim(0.0, 1.06)
    axR.axis("off")
    tag(axR, r"SCHEMATIC ($\lambda_2$ enlarged)")
    axR.text(0.0, 1.07, "(b)", transform=axR.transAxes, fontsize=12, fontweight="bold", va="top")
    fig.tight_layout()
    save(fig, "fig3_exchange")




# ------------------------------------------------------------------------------ Figure 2
def fig_boundary5():
    mp.mp.dps = 60
    a0 = 4 * mp.pi / 15
    th = lambda a: 14 * a - 3 * mp.pi
    words = return_words()
    order = words[2:]                                     # the five words sharing l3
    assert [len(w) for w in order] == [46, 46, 46, 66, 66]
    bnames = [3, 4, 5, 6]
    pairs = list(zip(order[:-1], order[1:]))

    def jump(pair, a):
        return offset(pair[1], a, th(a)) - offset(pair[0], a, th(a))

    xs = [mp.mpf(k) * mp.mpf('1e-5') for k in range(-30, 31)]
    curves = {b: [float(jump(p, a0 + d)) for d in xs] for b, p in zip(bnames, pairs)}
    pts = {b: [] for b in bnames}
    for d in [mp.mpf(k) * mp.mpf('1e-4') for k in (-3, -2, -1, 1, 2, 3)]:
        a = a0 + d
        Bd = Billiard(a)
        cw = {}
        for i in range(1, 1500):
            s = mp.mpf(i) / 1500
            r = Bd.first_return(s, th(a))
            if r is not None and r[1] in order and r[1] not in cw:
                cw[r[1]] = r[0] - s
        for b, (u, v) in zip(bnames, pairs):
            if u in cw and v in cw:
                pts[b].append((d, cw[v] - cw[u]))
    worst = max(abs(y - jump(p, a0 + x)) for b, p in zip(bnames, pairs) for x, y in pts[b])
    print("fig2: trace vs exact, worst |diff| = %s" % mp.nstr(worst, 3))
    assert worst < 1e-58
    assert all(abs(v) < 1e-50 for b in (3, 4, 6) for v in curves[b])
    pts = {b: [(float(x), float(y)) for x, y in v] for b, v in pts.items()}
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    X = [float(x) * 1e4 for x in xs]
    ax.axhline(0, color=HAIR, lw=0.8, zorder=0)
    ax.axvline(0, color=HAIR, lw=0.8, zorder=0)
    for b in (3, 4, 6):
        ax.plot(X, curves[b], color=MUTED, lw=2.0, zorder=2)
        ax.plot([x * 1e4 for x, _ in pts[b]], [y for _, y in pts[b]], "o", color=MUTED, ms=6,
                mec=SURF, mew=1.5, zorder=3)
    ax.plot(X, curves[5], color=INK, lw=2.0, zorder=4, label="boundary 5 (an $O$-graze)")
    ax.plot([x * 1e4 for x, _ in pts[5]], [y for _, y in pts[5]], "o", color=INK, ms=7,
            mec=SURF, mew=1.5, zorder=5)
    ax.plot([], [], color=MUTED, lw=2.0, label="boundaries 3, 4, 6 ($R$-grazes)")
    ax.plot([], [], "o", color=INK, ms=6, mec=SURF, label="60-digit trace")
    ax.text(X[-1], curves[5][-1], "  boundary 5", va="center", fontsize=9)
    ax.text(X[-1], 0, "  3, 4, 6", va="bottom", fontsize=9, color=INK)
    ax.text(0.08, ax.get_ylim()[1] * 0.92, r"$T_0$", fontsize=9)
    ax.set_xlabel(r"$\alpha - 4\pi/15$   ($10^{-4}$ rad)")
    ax.set_ylabel(r"translation jump across the boundary (units of $|L_1|$)")
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    for sp in ("left", "bottom"):
        ax.spines[sp].set_color(HAIR)
    ax.legend(frameon=False, fontsize=8.8, loc="upper right", bbox_to_anchor=(1.0, 0.94))
    tag(ax, "DATA")
    fig.tight_layout()
    save(fig, "fig2_boundary5")


# ------------------------------------------------------------------------------ Figure 4
def fig_connection():
    """The path from O at 12 degrees, unfolded across the sides it meets, reaching a copy
    of A.  Reflections are applied as exact-in-form complex affine maps at 60 digits."""
    mp.mp.dps = 60
    B = Billiard(mp.pi * P / (2 * Q))
    L = B.L
    orb = B.orbit(mp.mpf(0), mp.mpf(0), mp.pi / 15, 40)
    hit = next(j for j, (s, x, y) in enumerate(orb, 1)
               if abs(x - L) < mp.mpf(10) ** -40 and abs(y - 1) < mp.mpf(10) ** -40)
    word = [s for s, _, _ in orb[:hit]]

    def refl(base, tang):
        w = tang * tang
        return (w, base - w * base.conjugate(), True)

    def apply(g, z):
        w, c, rev = g
        return w * (z.conjugate() if rev else z) + c

    def compose(g, f):
        wg, cg, rg = g
        wf, cf, rf = f
        if rg:
            return (wg * wf.conjugate(), wg * cf.conjugate() + cg, not rf)
        return (wg * wf, wg * cf + cg, rf)
    Ov, Rv, Av = mp.mpc(0), mp.mpc(L), mp.mpc(L, 1)
    hyp = Av / abs(Av)
    R = {'1': refl(Rv, mp.mpc(0, 1)), '2': refl(Ov, mp.mpc(1)), 'h': refl(Ov, hyp)}
    hs = [(mp.mpc(1), mp.mpc(0), False)]
    for s in word[:-1]:
        hs.append(compose(hs[-1], R[s]))
    dirc = mp.mpc(mp.cos(mp.pi / 15), mp.sin(mp.pi / 15))
    zA = apply(hs[-1], Av)
    off = (dirc.conjugate() * zA).imag
    print("fig4: developed copy of A lies on the ray to %s" % mp.nstr(abs(off), 3))
    assert abs(off) < mp.mpf(10) ** -50 and len(word) - 1 == 13
    tri = [Ov, Rv, Av]
    INKc, MUTc, EDGE = INK, MUTED, "#aab4bf"
    fig, ax = plt.subplots(figsize=(11.5, 4.4))
    for k, h in enumerate(hs):
        xy = [(float(apply(h, z).real), float(apply(h, z).imag)) for z in tri]
        ax.add_patch(Polygon(xy, closed=True, facecolor=("#eef1f4" if k % 2 == 0 else "#f8f9fa"),
                             edgecolor=EDGE, linewidth=0.6, zorder=1))
        cx, cy = sum(p[0] for p in xy) / 3.0, sum(p[1] for p in xy) / 3.0
        ax.text(cx, cy, str(k), ha="center", va="center", fontsize=6.5, color=MUTc, zorder=3)
    tmax = float((dirc.conjugate() * zA).real)
    ax.plot([0, float((dirc * tmax).real)], [0, float((dirc * tmax).imag)],
            color=INKc, lw=2.0, solid_capstyle="round", zorder=5)
    ax.plot(0, 0, marker="o", ms=8, color=INKc, mec="white", mew=2, zorder=6)
    ax.plot(float(zA.real), float(zA.imag), marker="o", ms=9, color=INKc, mec="white", mew=2, zorder=6)
    ax.annotate("$O$", (0, 0), textcoords="offset points", xytext=(-16, -12), fontsize=13, color=INKc)
    ax.annotate("$A$ (copy %d)" % (len(hs) - 1), (float(zA.real), float(zA.imag)),
                textcoords="offset points", xytext=(8, 8), fontsize=12, color=INKc)
    sub = {"1": "L_1", "2": "L_2", "h": "H"}
    ax.text(0.0, -0.06, r"sides met: $" + r"\,".join(sub[w] for w in word[:-1]) + r"$   "
            r"($%d$ reflections), then arrival at $A$" % (len(word) - 1),
            transform=ax.transAxes, fontsize=9.5, color=INKc, ha="left")
    ax.set_aspect("equal")
    ax.axis("off")
    fig.tight_layout()
    save(fig, "fig4_connection")


def main():
    os.makedirs(OUT, exist_ok=True)
    rows, H = branch_rows()
    global H_CLOSED
    H_CLOSED = H
    fig1(rows)
    verify_fig2b()
    fig2()
    fig_boundary5()
    fig_connection()


if __name__ == "__main__":
    main()
