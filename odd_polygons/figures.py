"""[Deposit] Figures 1-4 (content asserted against exact data first). Paper: "Perpendicular billiards in rational polygons of odd denominator". Research name: probes/s956_figures.py, copied verbatim except imports (deposit names), paths (data/ and figs/ under odd_polygons/) and this paragraph. Run through run_all.py; the usage lines below are the research repo's. Bracketed IDs and section marks (§sNNN) refer to the project's internal ledgers.

s956 (`s955-e`): the figures of paper 2 (`tex/oddpoly/oddpoly.tex`).

  Fig 1  (§5.2, Prop 5.5) the rhombus (3,2,3,2)π/5: (a) its three swapped pairs — the short diagonal V0V2 and the mirror couple,
         the loops at V1 and V3; (b) the free cylinder Z unfolded along σ1 σ3 σ2 σ0 (leaves from σ0 at 54°, s ∈ (φ⁻², 1)), a straight
         strip whose two edges are the loops.
  Fig 2  (§6.2, Theorem C) S₁(1; π/17) = Π_17: (a) the quadrilateral with J = J0 ∪ J1 on σ3 and the torus's strips; (b) the
         bands developed along W0 (from J0) and W1 (from J1), flow horizontal, the same start copy; each ends on a translate.
  Fig 3  (§7.2) the trapezoid (1, 2, 5, 6)π/7: (a) left cap ⊕ strip ⊕ right cap, a perpendicular launch from σ1; (b) the left
         cap (angle π/N at O) developed around O into the regular N-gon (Lemma 7.3; Lemma 7.5 with c = 1): the launch's left-cap
         part is a chord, and the foot F of the perpendicular from O lies on a copy of a real side.
  Fig 4  (App A) Π_31 schematic: the chains T_0–T_1–T_2 and h(T_m), glued along s_0, s_1 (Theorem A.3); s_2 facing the
         invariant cylinders (measured, App A). No computation.
CONTROLS (asserted): Fig 1 — every leaf at s ∈ (φ⁻², 1) has word 1,3,2,0, returns to its start point and direction, has length
  2φH; the band height (1 − φ⁻²) sin 54° = 1/2; the unfolding's last copy is a translate by 2φH along the flow; the band edges
  pass through copies of V3 and V1 (Prop 5.5). Fig 2 — the band words read off the trace equal the paper's W0, W1 (§6.2 (i));
  |J_i| = λ_i/(2 cos((6m−1)x)) (§6.2 (ii)(b) as corrected s956: the ledger's λ are in the transverse coordinate X = −2 Re);
  J0 returns onto [lo + |J1|, hi1), J1 onto [lo, lo + |J1|); the pins (lo → V0 after 6m−1 letters, hi0 → V2 after 12m+4,
  hi1 → V0 after 6m+1; float, to 1e-8); both developments end on a translate of the start copy. Fig 3 — the launch meets a real
  side of the left cap at a right angle where the chord's perpendicular foot develops (Lemma 7.3), the foot is at a multiple of
  π/N, the folded path and the chord have the same length; the 2N copies close up.
Paths are float-traced, for drawing only. Not pre-registered: figures. PRIOR ART: probes/s949_tri_figures.py (paper 1, same
style), probes/s941_fig_decomp.py (Π_31 decomposition graph, research version). DATA: none.
Output: probes/figs/s956_fig{1,2,3,4}.pdf (+ .png), copied to tex/oddpoly/fig{1,2,3,4}.pdf.
Run from the repo root:
PYTHONPATH=.:engine:archive/scripts_2026-07:archive/scripts_2026-08:probes .venv/bin/python3.13 probes/s956_figures.py [1 2 3 4]
"""
import _setup  # noqa: F401  (deposit: this directory on the path, cwd = odd_polygons/)
import math
import shutil
import sys

import matplotlib
import matplotlib.patches
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
plt.rcParams["pdf.fonttype"] = 42                  # TrueType: Hexagon safe-PDF (s900, s909)
plt.rcParams["mathtext.fontset"] = "cm"
import numpy as np

BLUE, ORANGE, AQUA, VIOLET = '#2a78d6', '#eb6834', '#1baf7a', '#6250d6'
INK, MUTED, FAINT, FILL = '#222222', '#8a8a85', '#d6d5cf', '#f4f3ef'
TOL = 1e-10
PHI = (1 + 5 ** 0.5) / 2


# ------------------------------------------------------------------------------------------------- geometry
def polygon(angles, lengths):
    """V0 at 0, side σ0 = V0V1 along +x, counter-clockwise; angles[i] = interior angle at V_i (radians)."""
    V, th, p = [], 0.0, np.zeros(2)
    for i in range(len(lengths)):
        V.append(p.copy())
        p = p + lengths[i] * np.array([math.cos(th), math.sin(th)])
        th += math.pi - angles[(i + 1) % len(angles)]
    assert np.linalg.norm(p) < 1e-9, ('polygon does not close', p)
    return V


def unit(a):
    return np.array([math.cos(a), math.sin(a)])


def hit(x, v, p, q):
    """(s, u): x + s v meets p + u (q − p), s > 0, u ∈ [0, 1]; or None."""
    d = q - p
    M = np.array([[v[0], -d[0]], [v[1], -d[1]]])
    if abs(np.linalg.det(M)) < 1e-15:
        return None
    s, u = np.linalg.solve(M, p - x)
    return (s, u) if s > TOL and -1e-12 <= u <= 1 + 1e-12 else None


def reflect_pt(P, a, b):
    d = (b - a) / np.linalg.norm(b - a)
    w = P - a
    return a + 2 * np.dot(w, d) * d - w


def step(V, x, v, last):
    """next side hit from x along v (not the side `last`): (point, side, u)."""
    k = len(V)
    best = None
    for j in range(k):
        if j == last:
            continue
        h = hit(x, v, V[j], V[(j + 1) % k])
        if h and (best is None or h[0] < best[0]):
            best = (h[0], h[1], j)
    s, u, j = best
    return x + s * v, j, u


def trace(V, x, v, side, stop_side=None, maxl=500):
    """billiard path from x (on side `side`) along v. Stops at a vertex ('vertex', i), a right angle ('perp', j), or on
    returning to `stop_side` ('back', j). Returns pts, word, end, final direction."""
    k = len(V)
    pts, word, last = [x.copy()], [], side
    for _ in range(maxl):
        x, j, u = step(V, x, v, last)
        pts.append(x.copy())
        word.append(j)
        if u < 1e-9 or u > 1 - 1e-9:
            return pts, word, ('vertex', j if u < 0.5 else (j + 1) % k), v
        t = (V[(j + 1) % k] - V[j]) / np.linalg.norm(V[(j + 1) % k] - V[j])
        if abs(np.dot(v, t)) < 1e-9:
            return pts, word, ('perp', j), v
        v = 2 * np.dot(v, t) * t - v
        if j == stop_side:
            return pts, word, ('back', j), v
        last = j
    raise RuntimeError('trace did not stop')


def right_angle_mark(ax, P, side_dir, path_dir, size, color=INK):
    a, b = side_dir / np.linalg.norm(side_dir) * size, path_dir / np.linalg.norm(path_dir) * size
    ax.plot([P[0] + a[0], P[0] + a[0] + b[0], P[0] + b[0]], [P[1] + a[1], P[1] + a[1] + b[1], P[1] + b[1]],
            color=color, lw=0.7, zorder=6)


def draw_poly(ax, V, fc=FILL, ec=INK, lw=1.0, z=1):
    ax.add_patch(Polygon(V, closed=True, fc=fc, ec=ec, lw=lw, zorder=z))


def vlabel(ax, V, i, text, off=0.07, fs=10):
    c = sum(V) / len(V)
    d = V[i] - c
    d = d / np.linalg.norm(d)
    ax.text(*(V[i] + off * d), text, ha='center', va='center', fontsize=fs)


# ------------------------------------------------------------------------------------------------- Fig 1 (§5): rhombus
def fig_rhombus():
    A, B = math.radians(108), math.radians(72)
    V = polygon([A, B, A, B], [1, 1, 1, 1])
    H = math.sin(B)
    ang = math.radians(54)
    v0 = unit(ang)

    def leaf(s):
        return trace(V, np.array([s, 0.0]), v0.copy(), 0, stop_side=0)

    lo, hi = PHI ** -2, 1.0
    for s in np.linspace(lo, hi, 41)[1:-1]:
        pts, word, end, vf = leaf(s)
        assert word == [1, 3, 2, 0] and end == ('back', 0), (s, word, end)
        assert abs(pts[-1][0] - s) < 1e-9 and np.linalg.norm(vf - v0) < 1e-9
        L = sum(np.linalg.norm(pts[i + 1] - pts[i]) for i in range(len(pts) - 1))
        assert abs(L - 2 * PHI * H) < 1e-9, (L, 2 * PHI * H)
    assert abs((1 - PHI ** -2) * math.sin(ang) - 0.5) < 1e-12
    # boundary leaves: s = φ⁻² passes through V3; s = 1 starts at V1
    eps = 1e-6
    pl, _, _, _ = leaf(lo + eps)
    ph, _, _, _ = leaf(hi - eps)
    assert min(np.linalg.norm(p - V[3]) for p in pl) < 1e-6 or \
        min(abs((pl[i + 1] - pl[i])[0] * (V[3] - pl[i])[1] - (pl[i + 1] - pl[i])[1] * (V[3] - pl[i])[0]) / np.linalg.norm(pl[i + 1] - pl[i]) for i in range(4)) < 1e-5
    assert np.linalg.norm(ph[0] - V[1]) < 1e-5

    # the unfolding along the word: the base copy and its reflections; the fourth is a translate of the base
    copies, T = [list(V)], list(V)
    for j in [1, 3, 2, 0]:
        a, b = T[j], T[(j + 1) % 4]
        T = [reflect_pt(P, a, b) for P in T]
        copies.append(T)
    shift = copies[-1][0] - V[0]
    assert all(np.linalg.norm(copies[-1][i] - V[i] - shift) < 1e-9 for i in range(4)), 'not a translate'
    L = 2 * PHI * H
    assert abs(np.linalg.norm(shift) - L) < 1e-9 and abs(np.dot(shift / L, v0) - 1) < 1e-9
    # the band edges pass through copies of V3 and V1
    P0, P1 = np.array([lo, 0.0]), np.array([hi, 0.0])
    def on_line(Q, P):
        w = Q - P
        return abs(v0[0] * w[1] - v0[1] * w[0]) < 1e-9
    v3s = [T[3] for T in copies if on_line(T[3], P0)]
    v1s = [T[1] for T in copies if on_line(T[1], P1)]
    assert v3s and v1s

    fig, (ax, bx) = plt.subplots(1, 2, figsize=(5.0, 2.0), gridspec_kw=dict(width_ratios=[1, 1.75]))
    # (a) the three swapped pairs
    draw_poly(ax, V, fc='white')
    ax.plot(*np.array(ph).T, color=ORANGE, lw=1.5, zorder=4)
    ax.plot(*np.array(pl).T, color=BLUE, lw=1.5, zorder=4)
    ax.plot(*np.array([V[0], V[2]]).T, color=VIOLET, lw=1.3, ls=(0, (4, 2)), zorder=4)
    for i in range(4):
        ax.plot(*V[i], 'o', color=INK, ms=3, zorder=6)
    for i, t in enumerate(['$V_0$', '$V_1$', '$V_2$', '$V_3$']):
        vlabel(ax, V, i, t, off=0.10)
    ax.plot([], [], color=ORANGE, lw=1.5, label=r'loop at $V_1$')
    ax.plot([], [], color=BLUE, lw=1.5, label=r'loop at $V_3$')
    ax.plot([], [], color=VIOLET, lw=1.3, ls=(0, (4, 2)), label=r'diagonal $V_0V_2$')
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, 0.02), ncol=1, fontsize=8, frameon=False, handlelength=2.2)
    ax.set_aspect('equal'); ax.axis('off')
    # (b) the cylinder unfolded, drawn with the flow horizontal
    rot = -ang
    Rm = np.array([[math.cos(rot), -math.sin(rot)], [math.sin(rot), math.cos(rot)]])
    R = lambda P: Rm @ P
    for i, T in enumerate(copies):
        draw_poly(bx, [R(P) for P in T], fc=FILL if 0 < i < 4 else '#e3e8f3', ec=MUTED if 0 < i < 4 else INK,
                  lw=0.7 if 0 < i < 4 else 0.9)
    for i, T in enumerate(copies):                  # s957: number the copies 1-5 (reviewer E row 41)
        cen = R(sum(T) / 4)
        mid_y = R(0.5 * (np.array([lo, 0.0]) + np.array([hi, 0.0])))[1]
        far = max((R(P) for P in T), key=lambda q: abs(q[1] - mid_y))   # the copy's tip away from the band
        bx.text(*(cen + 0.55 * (far - cen)), str(i + 1), ha='center', va='center', fontsize=7.5, color=INK, zorder=7)
    band = [R(P0), R(P1), R(P1 + shift), R(P0 + shift)]
    bx.add_patch(Polygon(band, closed=True, fc=AQUA, alpha=0.22, ec='none', zorder=2))
    bx.plot(*np.array([R(P1), R(P1 + shift)]).T, color=ORANGE, lw=1.5, zorder=4)
    bx.plot(*np.array([R(P0), R(P0 + shift)]).T, color=BLUE, lw=1.5, zorder=4)
    for P in (P0, P1):
        pass
    bx.plot(*np.array([R(V[0]), R(V[1])]).T, color=INK, lw=1.8, zorder=5)
    bx.plot(*np.array([R(copies[-1][0]), R(copies[-1][1])]).T, color=INK, lw=1.8, zorder=5)
    for Q in v3s:
        bx.plot(*R(Q), 'o', color=BLUE, ms=4, zorder=6)
    for Q in v1s:
        bx.plot(*R(Q), 'o', color=ORANGE, ms=4, zorder=6)
    c = 0.72 * band[0] + 0.28 * band[1] + 0.5 * (band[3] - band[0]) + 0.28 * (band[2] - band[3])
    bx.text(*c, '$Z$', ha='center', va='center', fontsize=11)
    bx.text(*(R(V[0]) + R(np.array([0.5, -0.12]))), r'$\sigma_0$', ha='center', va='center', fontsize=9)
    bx.set_aspect('equal'); bx.axis('off')
    fig.tight_layout()
    fig.text(0.01, 0.97, '(a)', fontsize=9, va='top'); fig.text(0.40, 0.97, '(b)', fontsize=9, va='top')
    return fig, dict(band=(lo, hi), leaf=L, v3_on_edge=len(v3s), v1_on_edge=len(v1s))


# ------------------------------------------------------------------------------------------------- Fig 3 (§7): trapezoid
def fig_trap():
    N, a0, a1, D = 7, 1, 2, 1.2
    x = math.pi / N
    l0, l1 = 1 / math.tan(a0 * x), 1 / math.tan(a1 * x)      # cap widths at height 1
    V = [np.array([0.0, 0.0]), np.array([l0 + D + l1, 0.0]), np.array([l0 + D, 1.0]), np.array([l0, 1.0])]
    BL, BR = np.array([l0, 0.0]), np.array([l0 + D, 0.0])
    # a perpendicular launch from σ1 (the right cap's slanted side) that crosses the strip and ends at a right angle in the
    # left cap; among launch points pick the one whose development has |OF| closest to 0.55 |OB_L| (a legible chord)
    nrm = np.array([-(V[2] - V[1])[1], (V[2] - V[1])[0]]); nrm /= np.linalg.norm(nrm)     # inward normal of σ1
    best = None
    for u in np.linspace(0.15, 0.85, 351):
        x0 = V[1] + u * (V[2] - V[1])
        pts, word, end, vf = trace(V, x0, nrm.copy(), 1)
        if end[0] != 'perp' or pts[-1][0] >= l0 - 1e-6 or len(pts) > 6:
            continue
        cross = sum(1 for k in range(len(pts) - 1) if (pts[k][0] - l0) * (pts[k + 1][0] - l0) < 0)
        if cross != 1:
            continue
        bounce = [P for P in pts[1:-1]]
        if any(min(np.linalg.norm(P - BL), np.linalg.norm(P - BR), np.linalg.norm(P - V[2]), np.linalg.norm(P - V[3])) < 0.15
               for P in bounce):
            continue
        k_in = max(i for i in range(len(pts) - 1) if pts[i][0] >= l0 - 1e-12 and pts[i + 1][0] < l0)
        dv = (pts[k_in + 1] - pts[k_in]) / np.linalg.norm(pts[k_in + 1] - pts[k_in])
        hy = (pts[k_in] + (l0 - pts[k_in][0]) / dv[0] * dv)[1]
        if not 0.2 < hy < 0.8:
            continue
        sc = abs(np.linalg.norm(pts[-1]) / l0 - 0.55)
        if best is None or sc < best[0]:
            best = (sc, u, pts, word, end, vf)
    assert best is not None
    _, u, pts, word, end, vf = best
    ang = math.atan2(nrm[1], nrm[0])
    k_in = max(i for i in range(len(pts) - 1) if pts[i][0] >= l0 - 1e-12 and pts[i + 1][0] < l0)
    seg = pts[k_in:]
    dirv = (seg[1] - seg[0]) / np.linalg.norm(seg[1] - seg[0])
    tt = (l0 - seg[0][0]) / dirv[0]
    entry = seg[0] + tt * dirv                                 # on the virtual side x = l0
    in_cap = [entry] + seg[1:]
    # the development: the line through entry in the entry direction, in the base copy (the cap itself)
    O = V[0]
    foot = entry + np.dot(O - entry, dirv) * dirv
    psiF = math.atan2(foot[1], foot[0])
    assert abs(psiF / x - round(psiF / x)) < 1e-9, psiF / x   # Lemma 7.3: the foot lies on a real-side copy
    Lpath = sum(np.linalg.norm(in_cap[i + 1] - in_cap[i]) for i in range(len(in_cap) - 1))
    assert abs(Lpath - np.linalg.norm(foot - entry)) < 1e-9   # the folded path has the chord's length up to F
    assert abs(np.linalg.norm(foot - O) - np.linalg.norm(in_cap[-1] - O)) < 1e-9

    fig, (ax, bx) = plt.subplots(1, 2, figsize=(5.0, 2.15), gridspec_kw=dict(width_ratios=[1.45, 1]))
    # (a)
    draw_poly(ax, V, fc='white')
    ax.add_patch(Polygon([V[0], BL, V[3]], closed=True, fc=BLUE, alpha=0.10, ec='none'))
    ax.add_patch(Polygon([BR, V[1], V[2]], closed=True, fc=ORANGE, alpha=0.10, ec='none'))
    for B, T in ((BL, V[3]), (BR, V[2])):
        ax.plot(*np.array([B, T]).T, color=MUTED, lw=0.8, ls=(0, (3, 2)))
    ax.plot(*np.array(pts).T, color=INK, lw=1.0, zorder=4)
    ax.plot(*np.array(in_cap).T, color=AQUA, lw=1.8, zorder=5)
    j = end[1]
    right_angle_mark(ax, pts[-1], V[(j + 1) % 4] - V[j], pts[-2] - pts[-1], 0.07)
    right_angle_mark(ax, pts[0], V[2] - V[1], pts[1] - pts[0], 0.07)
    for i, t in enumerate(['$V_0$', '$V_1$', '$V_2$', '$V_3$']):
        vlabel(ax, V, i, t, off=0.17)
    ax.text(*(BL + np.array([0, -0.1])), '$B_L$', ha='center', va='top', fontsize=9)
    ax.text(*(BR + np.array([0, -0.1])), '$B_R$', ha='center', va='top', fontsize=9)
    ax.text(l0 * 0.5, -0.12, 'left cap', ha='center', va='top', fontsize=8, color=BLUE)
    ax.text(l0 + D * 0.42, -0.12, 'strip', ha='center', va='top', fontsize=8, color=MUTED)
    ax.text(l0 + D + l1 * 0.62, -0.3, 'right cap', ha='center', va='top', fontsize=8, color=ORANGE)
    ax.annotate('', xy=(l0, 1.2), xytext=(l0 + D, 1.2), arrowprops=dict(arrowstyle='<->', lw=0.7, color=MUTED))
    ax.text(l0 + D / 2, 1.25, '$D$', ha='center', va='bottom', fontsize=9, color=MUTED)
    ax.set_title('(a)', fontsize=9, loc='left')
    ax.set_ylim(-0.45, 1.45)
    ax.set_aspect('equal'); ax.axis('off')
    # (b): 2N copies of the left cap around O
    cap = [O, BL, V[3]]
    copies, T = [], [p.copy() for p in cap]
    for i in range(2 * N):
        copies.append(T)
        # alternate reflections in the two real sides through O: side O-B (index 1) then O-V3 (index 2)
        k = 1 if i % 2 == 0 else 2
        T = [T[0], reflect_pt(T[1], T[0], T[k]) if k == 2 else T[1], reflect_pt(T[2], T[0], T[1]) if k == 1 else T[2]]
    assert np.allclose(T, cap), 'copies do not close up'
    for i, T in enumerate(copies):
        draw_poly(bx, T, fc=BLUE if i == 0 else FILL, ec=FAINT, lw=0.5)
        if i == 0:
            bx.patches[-1].set_alpha(0.18)
    for T in copies:
        for k in (1, 2):
            bx.plot(*np.array([T[0], T[k]]).T, color=MUTED, lw=0.6, zorder=2)
    # regular N-gon outline: virtual-side copies
    for T in copies:
        bx.plot(*np.array([T[1], T[2]]).T, color=INK, lw=0.9, zorder=3)
    # the chord: from entry to the far boundary; the actual path is entry -> F
    exits = [h for T in copies for h in [hit(entry + 1e-9 * dirv, dirv, T[1], T[2])] if h]
    far = entry + min(h[0] for h in exits) * dirv
    bx.plot(*np.array([entry, foot]).T, color=AQUA, lw=1.8, zorder=5)
    bx.plot(*np.array([foot, far]).T, color=AQUA, lw=0.9, ls=(0, (3, 2)), zorder=5)
    bx.plot(*np.array([O, foot]).T, color=INK, lw=1.2, zorder=4)
    right_angle_mark(bx, foot, O - foot, entry - foot, 0.06)
    bx.plot(*O, 'o', color=INK, ms=3, zorder=6)
    bx.text(O[0] - 0.06, O[1] - 0.06, '$O$', ha='right', va='top', fontsize=9)
    bx.text(foot[0] - 0.10, foot[1] + 0.04, '$F$', ha='right', va='bottom', fontsize=9)
    bx.set_title('(b)', fontsize=9, loc='left')
    bx.set_aspect('equal'); bx.axis('off')
    fig.tight_layout()
    return fig, dict(prong_dir_units=round(ang / (math.pi / (2 * N)), 6), word=word, end=end,
                     foot_angle_units=round(psiF / x, 6))


# ------------------------------------------------------------------------------------------------- Fig 2 (§6): S_1(1; pi/17)
def s1_template(m):
    W0 = '23' * (3 * m - 1) + '2' + '013' + '23' * (m - 1) + '2' + '02' + '32' * (2 * m - 1) + '3' + '102' + '32' * (2 * m - 1) + '3' + '1013'
    W1 = W0.replace('102' + '32' * (2 * m - 1) + '3', '', 1)
    assert len(W0) == 16 * m + 8 and len(W1) == 12 * m + 6
    return W0, W1


def fig_s1(m=1, N=17):
    x = math.pi / N
    V = polygon([math.pi - (3 * m + 1) * x, 4 * m * x, math.pi - (m + 1) * x, 2 * x],
                [2 * math.sin(x) * math.cos(m * x), 2 * math.sin(x) * math.cos(3 * m * x), math.sin(4 * m * x), math.sin(4 * m * x)])
    W0, W1 = s1_template(m)
    d = unit(math.pi / 2 + (3 * m - 2) * x)
    t3 = (V[0] - V[3]) / np.linalg.norm(V[0] - V[3]); S = np.linalg.norm(V[0] - V[3])

    P0 = lambda s: V[3] + s * t3

    def band(s, maxl=400):
        """leaf from the point of σ3 at distance s from V3: (points, word, end, return position)."""
        p, v, last, w, pts = V[3] + s * t3, d.copy(), 3, '', [V[3] + s * t3]
        for _ in range(maxl):
            try:
                p, j, u = step(V, p, v, last)
            except TypeError:
                return pts, w, 'lost', None
            pts.append(p); w += str(j)
            if u < 1e-11 or u > 1 - 1e-11:
                return pts, w, ('vertex', j if u < 0.5 else (j + 1) % 4), None
            tt = V[(j + 1) % 4] - V[j]; tt /= np.linalg.norm(tt)
            v = 2 * np.dot(v, tt) * tt - v
            if j == 3 and np.linalg.norm(v - d) < 1e-9:
                return pts, w, 'back', float(np.dot(p - V[3], t3))
            last = j
        return pts, w, 'long', None

    def cls(s):
        _, w, e, _ = band(s)
        return 'W0' if (e == 'back' and w == W0) else 'W1' if (e == 'back' and w == W1) else 'X'

    g = np.linspace(0, S, 3001)[1:-1]
    c = [cls(s) for s in g]
    i0 = c.index('W0'); k0 = max(i for i in range(len(c)) if c[i] == 'W0'); k1 = max(i for i in range(len(c)) if c[i] == 'W1')
    assert all(ci == 'W0' for ci in c[i0:k0 + 1]) and all(ci == 'W1' for ci in c[k0 + 1:k1 + 1])

    def bis(a, b):
        ca = cls(a)
        for _ in range(60):
            mid = (a + b) / 2
            if cls(mid) == ca:
                a = mid
            else:
                b = mid
        return (a + b) / 2
    lo, hi0, hi1 = bis(g[i0 - 1], g[i0]), bis(g[k0], g[k0 + 1]), bis(g[k1], g[k1 + 1])
    lam0 = 4 * math.sin(x) * math.sin(2 * m * x) * math.sin(4 * m * x) * (2 * math.cos(2 * m * x) - 1)
    lam1 = 16 * math.sin(x) * math.sin(m * x) ** 2 * math.sin(2 * m * x) * math.sin(4 * m * x)
    cf = 2 * math.cos((6 * m - 1) * x)
    J0, J1 = hi0 - lo, hi1 - hi0
    assert abs(J0 - lam0 / cf) < 1e-9 and abs(J1 - lam1 / cf) < 1e-9, (J0, lam0 / cf, J1, lam1 / cf)
    # the return positions (§6.2 (ii)(b)) and the pins (§6.2 (i))
    for s in np.linspace(lo, hi0, 9)[1:-1]:
        assert abs(band(s)[3] - (s + J1)) < 1e-9
    for s in np.linspace(hi0, hi1, 5)[1:-1]:
        assert abs(band(s)[3] - (s - J0)) < 1e-9
    for s0, wd, n, vtx in ((lo, W0, 6 * m - 1, 0), (hi0, W0, 12 * m + 4, 2), (hi1, W1, 6 * m + 1, 0)):
        p, v, last, w = P0(s0), d.copy(), 3, ''
        for _ in range(n + 1):                         # the (n+1)-th segment ends at the vertex (float: to 1e-8)
            p, j, u = step(V, p, v, last); w += str(j)
            tt = V[(j + 1) % 4] - V[j]; tt /= np.linalg.norm(tt)
            v = 2 * np.dot(v, tt) * tt - v; last = j
        assert w[:n] == wd[:n] and np.linalg.norm(p - V[vtx]) < 1e-8, (s0, w, wd[:n], np.linalg.norm(p - V[vtx]))

    # the development along W0 (from J0) and W1 (from J1)
    def develop(word):
        T, out = list(V), [list(V)]
        for ch in word:
            j = int(ch)
            a, b = T[j], T[(j + 1) % 4]
            T = [reflect_pt(P, a, b) for P in T]
            out.append(T)
        return out
    C0, C1 = develop(W0), develop(W1)
    sh0, sh1 = C0[-1][3] - V[3], C1[-1][3] - V[3]
    for C, sh in ((C0, sh0), (C1, sh1)):
        assert all(np.linalg.norm(C[-1][i] - V[i] - sh) < 1e-9 for i in range(4))     # (a): a translate
    P = lambda s: V[3] + s * t3
    # the endpoints of the bands land where (b) says
    def par(u):
        return abs(u[0] * d[1] - u[1] * d[0])
    assert par(P(lo + J1) + sh0 - P(lo)) < 1e-9 and par(P(hi1) + sh0 - P(hi0)) < 1e-9
    assert par(P(lo) + sh1 - P(hi0)) < 1e-9 and par(P(lo + J1) + sh1 - P(hi1)) < 1e-9

    rot = -math.atan2(d[1], d[0])
    Rm = np.array([[math.cos(rot), -math.sin(rot)], [math.sin(rot), math.cos(rot)]])
    R = lambda Q: Rm @ (Q - V[3])

    fig = plt.figure(figsize=(5.0, 3.5))           # s957: (a) enlarged on its own row (reviewer E row 39)
    ax = fig.add_axes([0.27, 0.6, 0.46, 0.4])
    bx = fig.add_axes([0.0, 0.0, 1.0, 0.6])
    # (a) the polygon, one leaf from J0 (solid) and one from J1 (dashed), J0 and J1 labelled
    draw_poly(ax, V, fc='white', z=1)
    for (s_, col, ls_) in ((lo + 0.5 * J0, BLUE, '-'), (hi0 + 0.5 * J1, ORANGE, (0, (3, 1.5)))):
        pts_, _, _, _ = band(s_)
        ax.plot(*np.array(pts_).T, color=col, lw=0.45, ls=ls_, zorder=3)
    ax.plot(*np.array([P(lo), P(hi0)]).T, color=BLUE, lw=3, solid_capstyle='butt', zorder=4)
    ax.plot(*np.array([P(hi0), P(hi1)]).T, color=ORANGE, lw=3, solid_capstyle='butt', zorder=4)
    for i, tx in enumerate(['$V_0$', '$V_1$', '$V_2$', '$V_3$']):
        vlabel(ax, V, i, tx, off=0.03, fs=9)
    nn = np.array([-t3[1], t3[0]])
    ax.text(*(P(lo + 0.5 * J0) - 0.05 * nn), '$J_0$', ha='center', va='center', fontsize=9, color=BLUE)
    ax.text(*(P(hi0 + 0.5 * J1) - 0.05 * nn), '$J_1$', ha='center', va='center', fontsize=9, color=ORANGE)
    ax.set_aspect('equal'); ax.axis('off')
    # (b) the developments, flow horizontal, in one frame: W0 on top, W1 below (offset), the same start copy
    gap = np.array([0.0, -1.75])
    for C, s_a, s_b, col, sh, wname, la, lb, off in ((C0, lo, hi0, BLUE, sh0, 'W_0', lo + J1, hi1, np.zeros(2)),
                                                    (C1, hi0, hi1, ORANGE, sh1, 'W_1', lo, lo + J1, gap)):
        RR = lambda Q: R(Q) + off
        for i, T in enumerate(C):
            first_last = i == 0 or i == len(C) - 1
            draw_poly(bx, [RR(Q) for Q in T], fc='#e3e8f3' if first_last else FILL, ec=INK if first_last else FAINT,
                      lw=0.8 if first_last else 0.4)
        A, B, A2, B2 = RR(P(s_a)), RR(P(s_b)), RR(P(la) + sh), RR(P(lb) + sh)
        bx.add_patch(Polygon([A, B, B2, A2], closed=True, fc=col, alpha=0.6, ec='none', zorder=3))
        if wname == 'W_1':                         # s957: the thin band also dashed-outlined (greyscale)
            bx.add_patch(Polygon([A, B, B2, A2], closed=True, fc='none', ec=INK, lw=0.5, ls=(0, (2, 1.5)), zorder=4))
        for o2 in (np.zeros(2), sh):
            for (u0, u1, cc) in ((lo, hi0, BLUE), (hi0, hi1, ORANGE)):
                Q0, Q1 = RR(P(u0) + o2), RR(P(u1) + o2)
                bx.plot([Q0[0], Q1[0]], [Q0[1], Q1[1]], color=cc, lw=2.6, solid_capstyle='butt', zorder=5)
        bx.text(*(RR(V[3]) + np.array([-0.25, 0.0])), f'${wname}$', ha='right', va='center', fontsize=10, color=col)
    for C, n, vtx, off in ((C0, 6 * m - 1, 0, np.zeros(2)), (C0, 12 * m + 4, 2, np.zeros(2)), (C1, 6 * m + 1, 0, gap)):
        Q = R(C[n][vtx]) + off              # the (n+1)-th segment runs in copy C[n] and ends at its vertex
        bx.plot(*Q, 'o', color=INK, ms=3, zorder=6)
    bx.set_aspect('equal'); bx.axis('off')
    fig.text(0.0, 0.99, '(a)', fontsize=9, va='top'); fig.text(0.0, 0.6, '(b)', fontsize=9, va='top')
    return fig, dict(lo=lo, hi0=hi0, hi1=hi1, J0=J0, J1=J1, cos_factor=cf, rot=J1 / (J0 + J1), shift0=float(np.dot(sh0, d)),
                     shift1=float(np.dot(sh1, d)))


# ------------------------------------------------------------------------------------------------- Fig 4 (App A): Pi_31's chains
def fig_chain(N=31):
    """schematic of App A at N = 31 (M = 2): T_m glued to T_{m+1} along s_m (Theorem A.3, proved); s_M facing the invariant
    cylinders, and the 2(M+1) tori being all the minimal components (both MEASURED, App A). No computation: the slit labels
    are Theorem A.3's, s_m = {l_m, c_{4m+5}, c_{4m+7}}."""
    M = (N - 15) // 8
    assert M == 2 and 2 * (M + 1) == 2 * ((N - 7) // 8)
    fig, ax = plt.subplots(figsize=(5.0, 1.2))
    xs = [-3.75, -2.55, -1.35, 0.0, 1.35, 2.55, 3.75]
    names = ['$T_0$', '$T_1$', '$T_2$', None, '$h(T_2)$', '$h(T_1)$', '$h(T_0)$']
    for i in range(6):
        measured = i in (2, 3)
        col = [BLUE, ORANGE, AQUA, AQUA, ORANGE, BLUE][i]
        ax.plot([xs[i], xs[i + 1]], [0, 0], color=col, lw=2.2, ls=(0, (3, 2)) if measured else '-', zorder=1)
        m = [0, 1, 2, 2, 1, 0][i]
        cx_ = -0.87 if i == 2 else 0.87 if i == 3 else (xs[i] + xs[i + 1]) / 2
        ax.text(cx_, 0.06, f'$s_{m}$' if i < 3 else f'$h(s_{m})$', ha='center', va='bottom', fontsize=8.5, color=col)
        if i < 3:
            ax.text(cx_, -0.08, rf'$\ell_{m},\,c_{{{4 * m + 5}}},\,c_{{{4 * m + 7}}}$', ha='center', va='top',
                    fontsize=8, color=INK)
    for x0, nm in zip(xs, names):
        if nm is None:
            ax.add_patch(matplotlib.patches.FancyBboxPatch((x0 - 0.38, -0.2), 0.76, 0.4, boxstyle='round,pad=0.02',
                                                           fc=FILL, ec=MUTED, lw=0.8, zorder=2))
            ax.text(x0, 0, 'invariant\ncylinders', ha='center', va='center', fontsize=6.5, color=INK, zorder=3)
        else:
            ax.plot(x0, 0, 'o', color=INK, ms=8, zorder=3)
            ax.text(x0, 0.17, nm, ha='center', va='bottom', fontsize=8.5)
    ax.set_xlim(-4.0, 4.0); ax.set_ylim(-0.4, 0.42)
    ax.set_aspect('equal'); ax.axis('off')
    fig.tight_layout()
    return fig, dict(N=N, M=M)


FIGS = {'1': fig_rhombus, '2': fig_s1, '3': fig_trap, '4': fig_chain}     # numbered in print order

if __name__ == '__main__':
    import os
    os.makedirs('figs', exist_ok=True)
    which = sys.argv[1:] or sorted(FIGS)
    for name in which:
        fig, info = FIGS[name]()
        for ext in ('pdf', 'png'):
            fig.savefig(f'figs/fig{name}.{ext}', dpi=200, bbox_inches='tight')
        print('fig' + name, info)
    print('ALL OK')
