"""[Deposit] Figure 1 (the polygon P, Section 3.1). Paper: "Finite perpendicular-beam partitions at irrational angles". Research name: probes/s987_beam2_figures.py (which draws the figures of the companion paper; its Figure 1 is this paper's Figure 1), copied verbatim except the output path (figs/ under irrational_beams/), save() no longer copying into the paper's source tree, and the main block drawing Figure 1 only. Run through run_all.py.

s987 (`s974-b`): the three figures of paper 2 (`tex/beam2/beam2.tex`).

  Fig 1  (§3.1) the polygon P: copies T_j of T around O, the sides S_k tangent to the circle of radius r = cot α,
         and one trajectory leaving S_0 horizontally at height s; the moment p is the distance of a chord from O.
  Fig 2  (§4.1) (a) one flight of the kicked particle (4.1) at G = 1/24, the top cell of Theorem A(F)'s J = 3 limit;
         (b) the same fan in P at m = 3 (the root of (F) with J = 3): steps S_3, S_0, S_3, S_6, S_3, S_6.
  Fig 3  (§6, Example 6.4) (a) the two fans of Theorem A(ST), j = 2, m = 4, in P; (b) the cells of the limit calculus
         in the (s, e) plane near e* = 3/4, coloured by closing hit (white: not closed within the cap).
  Fig 4  (§7.2) the self-similarity of Proposition 7.3: the leg of G = 45/116 (theta = 13/29 = [0; 2, 4, 2, 1]) -- two
         rungs, then in (0, eps1), eps1 = 3/29, a scaled copy of the leg of G2 = 5/12 (t = 1/3 = [0; 2, 1]); cells and
         closing hits exact (limit calculus); asserted: 5 and 3 cells, copy ends = eps1 x the ends at G2, hits as in (S).
CONTROLS (asserted): every traced hit satisfies the reflection identities of Lemma 3.1 (t = (p' − p)/(2 cos ψ));
  Fig 2(b) and Fig 3(a) reproduce the step sequences of Theorem A and return perpendicularly at height 2c − s;
  Fig 3(b) at e* = 3/4 has exactly the two cells (0, ½) (hit 6) and (½, 1) (hit 2) (Theorem 7.4 / Example 6.4).
Not pre-registered: figures. PRIOR ART: probes/s963_beaminf_figures.py (style). DATA: none.
Output: probes/figs/s987_fig{1,2,3,4}.pdf (+ .png), copied to tex/beam2/fig{1,2,3,4}.pdf.
Run from the repo root:
PYTHONPATH=.:engine:archive/scripts_2026-07:archive/scripts_2026-08:probes .venv/bin/python3.13 probes/s987_beam2_figures.py
"""
import math
import shutil

import _setup  # noqa: F401  (deposit: this directory on the path, cwd = irrational_beams/)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
plt.rcParams["pdf.fonttype"] = 42                  # TrueType: Hexagon safe-PDF (s900, s909)
plt.rcParams["mathtext.fontset"] = "cm"
import numpy as np
from scipy.optimize import brentq

BLUE, RED = '#2a78d6', '#c8553d'
INK, MUTED, FAINT, FILL = '#222222', '#8a8a85', '#d6d5cf', '#eef3fb'
OUT = 'figs/fig{}'


def save(fig, k):
    for ext in ('pdf', 'png'):
        fig.savefig(OUT.format(k) + '.' + ext, dpi=200, bbox_inches='tight')
    print('wrote', OUT.format(k) + '.pdf')


# ---------------------------------------------------------------- geometry of P
def nk(k, a):
    return np.array([math.cos(2 * k * a), math.sin(2 * k * a)])


def tk(k, a):
    return np.array([-math.sin(2 * k * a), math.cos(2 * k * a)])


def side(k, a):
    r = 1 / math.tan(a)
    return np.array([r * nk(k, a) - tk(k, a), r * nk(k, a) + tk(k, a)])


def trace(a, s, m, maxhits=60, partial=False):
    """billiard in P from (r, s) in direction π, until a perpendicular hit. Returns hit list
    [(k, t, point)] and the step list. Next side = first line of S_{k ± Δ}, 1 ≤ Δ ≤ m, in the sense of sgn p, met
    from inside at |t| < 1 (Lemma 3.2)."""
    r = 1 / math.tan(a)
    X = np.array([r, s])
    d = np.array([-1.0, 0.0])
    k = 0
    hits = [(0, s, X.copy())]
    for _ in range(maxhits):
        p = X[0] * d[1] - X[1] * d[0]
        g = 1 if p > 0 else -1
        best = None
        for D in range(1, m + 1):
            kk = k + g * D
            n = nk(kk, a)
            dn = d @ n
            if dn <= 0:
                continue
            lam = (r - X @ n) / dn
            if lam <= 1e-12:
                continue
            Y = X + lam * d
            t = Y @ tk(kk, a)
            if abs(t) < 1 and (best is None or lam < best[0]):
                best = (lam, kk, t, Y)
        assert best is not None, 'no next side'
        lam, kk, t, Y = best
        n = nk(kk, a)
        psi_in = math.atan2(d @ tk(kk, a), d @ n)
        dn = d - 2 * (d @ n) * n
        p_out = Y[0] * dn[1] - Y[1] * dn[0]
        assert abs(t - (p_out - p) / (2 * math.cos(psi_in))) < 1e-9          # Lemma 3.1
        X, d, k = Y, dn, kk
        hits.append((kk, t, Y.copy()))
        if abs(psi_in) < 1e-9:
            return hits
    if partial:
        return hits
    raise RuntimeError('no perpendicular return')


def draw_P(ax, a, ks, circle=True, sectors=True, side_labels=(), lw=1.2, faint=()):
    r = 1 / math.tan(a)
    if circle:
        th = np.linspace(0, 2 * math.pi, 400)
        ax.plot(r * np.cos(th), r * np.sin(th), color=FAINT, lw=0.8, zorder=1)
    for k in ks:
        S = side(k, a)
        if sectors:
            for P_ in S:
                ax.plot([0, P_[0]], [0, P_[1]], color=FAINT, lw=0.5, zorder=1)
            R = r * nk(k, a)
            ax.plot([0, R[0]], [0, R[1]], color=FAINT, lw=0.5, ls=':', zorder=1)
        ax.plot(S[:, 0], S[:, 1], color=MUTED if k in faint else INK, lw=lw,
                ls='--' if k in faint else '-', zorder=3)
    for k, txt, off in side_labels:
        R = (r + off) * nk(k, a)
        ax.text(R[0], R[1], txt, ha='center', va='center', fontsize=9, color=INK)
    ax.plot([0], [0], 'o', color=INK, ms=3, zorder=4)


def draw_traj(ax, hits, color, lw=1.0, alpha=1.0):
    pts = np.array([h[2] for h in hits])
    ax.plot(pts[:, 0], pts[:, 1], color=color, lw=lw, alpha=alpha, zorder=5)


# ---------------------------------------------------------------- Fig 1: the polygon P
def fig1():
    m = 3
    a = math.pi / (2 * m) - 0.035
    r = 1 / math.tan(a)
    fig, ax = plt.subplots(figsize=(5.2, 5.2))
    ks = range(0, 2 * m + 1)
    draw_P(ax, a, ks, side_labels=[(1, r'$S_1$', 0.32), (2, r'$S_2$', 0.32), (3, r'$S_3$', 0.32),
                                   (4, r'$S_4$', 0.32), (5, r'$S_5$', 0.32), (6, r'$S_6$', -0.32)],
           faint=(2 * m,))
    ax.text(r + 0.14, -0.9, r'$S_0$', fontsize=9, ha='left', va='center')
    # the triangle T_0 and the sector labels
    T0 = np.array([[0, 0], [r, 0], [r, 1]])
    ax.fill(T0[:, 0], T0[:, 1], color=FILL, zorder=0)
    ax.text(r * 0.82, 0.5, r'$T_0$', fontsize=9, ha='center')
    ax.text(r + 0.08, 0.5, r'$L_1$', fontsize=9, ha='left', va='center')
    ax.text(0.12, -0.17, r'$O$', fontsize=9)
    ax.text(r - 0.08, -0.13, r'$R_0$', fontsize=8, ha='right')
    ax.text(r + 0.06, 1.05, r'$A$', fontsize=8)
    # one trajectory
    s = 0.22
    hits = trace(a, s, m, maxhits=4, partial=True)
    draw_traj(ax, hits, BLUE, lw=1.1)
    # the moment of the first chord: the foot of the perpendicular from O
    Y0, Y1 = hits[0][2], hits[1][2]
    u = (Y1 - Y0) / np.linalg.norm(Y1 - Y0)
    foot = Y0 - (Y0 @ u) * u
    ax.plot([0, foot[0]], [0, foot[1]], color=RED, lw=0.9, zorder=6)
    ax.text(foot[0] * 0.5 + 0.06, foot[1] * 0.5, r'$p$', color=RED, fontsize=10, ha='left', va='center')
    ax.plot([Y0[0]], [Y0[1]], 'o', color=BLUE, ms=3.5, zorder=6)
    ax.text(Y0[0] + 0.07, Y0[1], r'$s$', color=BLUE, fontsize=9, va='center')
    ax.set_aspect('equal')
    lim = r + 0.75
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.axis('off')
    steps = [h[0] for h in hits]
    print('fig1: m =', m, 'alpha =', round(a, 4), 'sides', steps)
    save(fig, 1)


# ---------------------------------------------------------------- Fig 2: a flight and its fan
def tg_orbit(G, s, nmax=40):
    """T_G orbit of (0, s/2) (Definition 4.2) until the fold; returns positions y_0..y_{ℓ-1}, momenta, and ℓ."""
    mu, y = 0, s / 2
    ys, mus = [y], [mu]
    for nu in range(nmax):
        mu2 = mu + (1 if y > 0 else -1)
        n = math.floor(y - mu2 * G + 0.5)
        y2 = y - 2 * mu2 * G - 2 * n
        if mu2 == -mu - (1 if y > 0 else -1) and n == 0:          # T(X) = J(X): on the fold
            return ys, mus, nu + 1
        mu, y = mu2, y2
        ys.append(y)
        mus.append(mu)
    raise RuntimeError


def F_root(m, J):
    def f(d):
        a = math.pi / (2 * m) - d
        return (-1) ** J * math.cos(a) * math.sin((2 * J + 1) * m * a) - math.sin((m - 1) * a)
    dstar = math.pi / (2 * m * (2 * m * J + 1))
    return math.pi / (2 * m) - brentq(f, 1e-12, dstar * (1 - 1e-12))


def fig2():
    J, i = 3, 3
    G = 1 / (2 * J * (J + 1))
    s = 0.8
    ys, mus, ell = tg_orbit(G, s)
    assert ell == 2 * i, ell
    # full period: second half is the J-reflection (same positions, reversed)
    yfull = ys + ys[::-1]
    nus = np.arange(len(yfull))
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(10.4, 4.6), gridspec_kw={'width_ratios': [1.05, 1]})
    ax.axhline(0, color=MUTED, lw=0.7)
    ax.axhspan(-0.5, 0.5, color=FILL, zorder=0)
    ax.plot(nus, yfull, color=FAINT, lw=0.9, zorder=2)
    ax.plot(nus[:ell], yfull[:ell], 'o-', color=BLUE, ms=5, lw=1.2, zorder=3, label='out')
    ax.plot(nus[ell - 1:], yfull[ell - 1:], 'o--', color=BLUE, ms=4, lw=0.9, mfc='white', zorder=3,
            label='back (the reflection $J$)')
    for nu in range(ell):
        ax.text(nu, ys[nu] + (0.035 if ys[nu] > 0 else -0.06), rf'${mus[nu]}$', ha='center', fontsize=8,
                color=MUTED)
    ax.annotate('at rest', xy=(ell - 0.5, ys[-1]), xytext=(ell - 0.5, ys[-1] - 0.11), ha='center', fontsize=8,
                arrowprops=dict(arrowstyle='-', color=MUTED, lw=0.6))
    ax.annotate('launch at rest', xy=(0, ys[0]), xytext=(0.6, ys[0] + 0.07), fontsize=8,
                arrowprops=dict(arrowstyle='-', color=MUTED, lw=0.6))
    ax.set_xlabel(r'step $\nu$')
    ax.set_ylabel(r'position $y_\nu$')
    ax.set_ylim(-0.52, 0.52)
    ax.set_xticks(range(0, len(yfull), 2))
    ax.set_title(r'(a) a flight of the particle, $G = 1/24$', fontsize=10)
    ax.spines[['top', 'right']].set_visible(False)
    # (b) the fan in P at m = 3
    m = 3
    a = F_root(m, J)
    r = 1 / math.tan(a)
    d = math.pi / (2 * m) - a
    b = [2 * r * sum(math.sin(2 * m * l * d) for l in range(1, k + 1)) for k in range(J + 1)]
    assert abs(b[J] - 1) < 1e-9
    draw_P(bx, a, range(-1, 2 * m + 1), sectors=False, faint=(2 * m,),
           side_labels=[(0, r'$S_0$', 0.32), (3, r'$S_3$', 0.42), (6, r'$S_6$', -0.3)])
    c = (b[i - 1] + b[i]) / 2
    for s_ in np.linspace(b[i - 1], b[i], 9)[1:-1]:
        hits = trace(a, s_, m)
        steps = [h[0] for h in hits]
        assert steps == [0, 3, 0, 3, 6, 3, 6], steps
        assert abs(abs(hits[-1][1]) - (2 * c - s_)) < 1e-9                 # returns at height 2c − s
        if abs(s_ - 0.8) > 0.07:
            draw_traj(bx, hits, BLUE, lw=0.6, alpha=0.25)
    hits = trace(a, 0.8, m)
    draw_traj(bx, hits, BLUE, lw=1.3)
    for q, (_, _, Y) in enumerate(hits[1:], start=1):
        off = 0.16 if Y[0] > 0 else -0.16
        bx.text(Y[0] + off, Y[1], rf'${q}$', fontsize=8, color=BLUE, ha='center', va='center')
    bx.plot([r, r], [b[i - 1], b[i]], color=RED, lw=3, solid_capstyle='butt', zorder=6)
    bx.text(r + 0.32, c + 0.12, 'fan', color=RED, fontsize=9, va='center')
    bx.set_aspect('equal')
    lim = r + 0.6
    bx.set_xlim(-lim, lim)
    bx.set_ylim(-lim, lim)
    bx.axis('off')
    bx.set_title(r'(b) the same fan in $P$: Theorem A(F), $J = 3$, $m = 3$', fontsize=10)
    print('fig2: flight y =', [round(y, 4) for y in ys], 'mu =', mus, '; fan', (round(b[i - 1], 4), 1.0),
          'alpha/pi =', round(a / math.pi, 6))
    save(fig, 2)


# ---------------------------------------------------------------- Fig 3: Example 6.4
def ST_root(m, j):
    def T(a):
        return (2 * math.cos(a) * math.sin(j * a) * math.sin((2 * m - j + 1) * a) + math.sin(a) ** 2
                - 2 * math.sin(j * a) ** 2)
    return brentq(T, math.pi / (2 * m + 1) + 1e-12, math.pi / (2 * m) - 1e-12)


def calc_close(e, s, cap):
    """vectorised limit calculus (Definition 4.4) over an array of launches s at one e: closing hit, 0 if not
    closed within cap (a graze, which has measure zero, is ignored)."""
    k1 = np.zeros_like(s, dtype=np.int64)
    k2 = np.zeros_like(s, dtype=np.int64)
    p = s.copy()
    out = np.zeros(s.shape, dtype=np.int64)
    live = np.ones(s.shape, dtype=bool)
    for h in range(1, cap + 1):
        g = np.where(p > 0, 1, -1)
        k1n = -k1 + g
        c = k1n * e - k2                                   # t' = -p + k1n e + k2' with k2' = -k2 - 2gδ
        delta = np.floor((c - p) / (2 * g) + 0.5).astype(np.int64)
        k2n = -k2 - 2 * g * delta
        closed = live & (k1n == 0) & (k2n == 0)
        out[closed] = h
        live &= ~closed
        p = -p + 2 * (k1n * e + k2n)
        k1, k2 = k1n, k2n
        if not live.any():
            break
    return out


def fig3():
    # (a) the two fans of Theorem A(ST), j = 2, at m = 4
    m, j = 4, 2
    a = ST_root(m, j)
    r = 1 / math.tan(a)
    eps = math.pi - 2 * m * a
    K = 4 * math.cos(a) / math.sin(a) ** 2
    phi = eps / 2

    def S(nu):
        k = nu // 2
        if nu % 2 == 0:
            return -K * math.sin(k * a) * math.cos(k * a + phi) * math.sin(a - phi)
        return K * math.sin((k + 1) * a) * math.cos(k * a + phi) * math.sin(phi)
    eta1 = (S(2) + S(1)) / 2 + math.cos(2 * a)
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(10.4, 4.8), gridspec_kw={'width_ratios': [1, 1.1]})
    draw_P(ax, a, range(-1, 2 * m + 1), sectors=False, faint=(2 * m,),
           side_labels=[(0, r'$S_0$', 0.3), (m, r'$S_4$', 0.3), (2 * m - 1, r'$S_7$', 0.3)])
    for fan, (lo, hi, col) in {1: (eta1, 1, RED), 2: (0, eta1, BLUE)}.items():
        for s_ in np.linspace(lo, hi, 5)[1:-1]:
            hits = trace(a, s_, m)
            steps = [h[0] for h in hits]
            # steps D_1 = (m, m); D_2 = (m, m-1, m-2, m-2, m-1, m) as side indices
            D = [m, m] if fan == 1 else [m, m - 1, m - 2, m - 2, m - 1, m]
            diffs = [abs(steps[q + 1] - steps[q]) for q in range(len(steps) - 1)]
            assert diffs == D, (fan, steps)
            draw_traj(ax, hits, col, lw=0.8, alpha=0.55)
        ax.plot([r, r], [lo, hi], color=col, lw=3, solid_capstyle='butt', zorder=6)
    ax.text(r + 0.12, (eta1 + 1) / 2, r'$D_1$', color=RED, fontsize=10, va='center')
    ax.text(r + 0.12, eta1 / 2, r'$D_2$', color=BLUE, fontsize=10, va='center')
    ax.set_aspect('equal')
    lim = r + 0.6
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.axis('off')
    G = eps / (2 * a)
    ax.set_title(rf'(a) Theorem A(ST), $j = 2$, $m = 4$ ($G = {G:.4f}$)', fontsize=10)
    # (b) cells of the limit calculus in (s, e)
    es = np.linspace(0.70, 0.80, 401)
    ss = (np.arange(1200) + 0.5) / 1200
    cap = 400
    img = np.zeros((len(es), len(ss)), dtype=np.int64)
    for row, e in enumerate(es):
        img[row] = calc_close(e, ss, cap)
    # control at e* = 3/4: two cells, (0, 1/2) closing at 6, (1/2, 1) closing at 2
    ctrl = calc_close(0.75, ss, cap)
    assert set(ctrl[ss < 0.5].tolist()) == {6} and set(ctrl[ss > 0.5].tolist()) == {2}, set(ctrl.tolist())
    levels = [0, 1, 3, 7, 15, 31, 63, 127, cap + 1]
    cols = ['#ffffff', '#f4b183', '#9ecae1', '#6baed6', '#4292c6', '#2171b5', '#08519c', '#08306b']
    cmap = ListedColormap(cols)
    norm = BoundaryNorm(levels, cmap.N)
    im = bx.imshow(img, origin='lower', aspect='auto', extent=(0, 1, es[0], es[-1]), cmap=cmap, norm=norm,
                   interpolation='nearest')
    bx.axhline(0.75, color=RED, lw=0.8, ls='--')
    bx.text(1.01, 0.75, r'$e^* = \frac{3}{4}$', color=RED, fontsize=9, va='center', transform=bx.get_yaxis_transform())
    bx.set_xlabel(r'launch $s$')
    bx.set_ylabel(r'$e = 2G$')
    cb = fig.colorbar(im, ax=bx, fraction=0.05, pad=0.13, ticks=[2, 5, 11, 23, 47, 95, 263])
    cb.ax.set_yticklabels(['2', '4–7', '8–15', '16–31', '32–63', '64–127', r'$\geq 128$'], fontsize=7)
    cb.set_label('closing hit (white: open at %d hits)' % cap, fontsize=8)
    bx.set_title(r'(b) cells of the limit calculus near $e^* = \frac{3}{4}$', fontsize=10)
    unclosed = [(round(float(e), 3), round(float((img[row] == 0).mean()), 3)) for row, e in enumerate(es)
                if row % 50 == 0]
    print('fig3: alpha/pi =', round(a / math.pi, 6), 'G =', round(G, 5), 'eta1 =', round(eta1, 5),
          '; unclosed fraction by e:', unclosed)
    save(fig, 3)


# ---------------------------------------------------------------- Fig 4: the self-similarity (Prop renorm)
def leg_cells(e, M=2000, cap=4000):
    """exact limit-calculus cells of the leg at rational e, merged from M midpoint launches by closing hit."""
    import sys
    sys.path.insert(0, 'probes')
    from s983_tg_class import walk
    from fractions import Fraction as Fr
    pts = [Fr(2 * i + 1, 2 * M) for i in range(M)]
    cells = []
    for s_ in pts:
        r_ = walk(e, s_, cap)
        assert r_[0] == 'close', (s_, r_)
        if cells and cells[-1][2] == r_[1]:
            cells[-1][1] = s_
        else:
            cells.append([s_, s_, r_[1]])
    return [(float(c[0]), float(c[1]), c[2]) for c in cells]


def fig4():
    from fractions import Fraction as Fr
    th = Fr(13, 29)
    G = Fr(1, 2) - th / 4
    t = Fr(1, 3)
    G2 = Fr(1, 2) - t / 4
    eps1 = 1 - 2 * th
    assert G == Fr(45, 116) and G2 == Fr(5, 12) and eps1 == Fr(3, 29)
    C = leg_cells(2 * G, M=2900)
    C2 = leg_cells(2 * G2, M=1200)
    assert len(C) == 5 and len(C2) == 3, (C, C2)
    Q1, Q2 = 2, 2 * 4 + 1                                   # c1 = 2, c2 = 4
    for (lo, hi, h), (lo2, hi2, h2) in zip(C[:3], C2):
        assert abs(lo - float(eps1) * lo2) < 2e-3 and abs(hi - float(eps1) * hi2) < 2e-3
        N2 = h2 // 2 - 1
        s2 = (lo2 + hi2) / 2
        assert h == 2 * (Q2 * (N2 + 1) + Q1 * (math.floor(s2 + (N2 + 1) * float(t)) + 1)), (h, h2)
    fig, ax = plt.subplots(figsize=(10.4, 3.0))
    pal = ['#08519c', '#4292c6', '#9ecae1', '#f4b183', '#e6550d']
    ends = [0, 1 / 29, 2 / 29, 3 / 29, 16 / 29, 1]
    hits_ = [c[2] for c in C]
    y1, y2, hgt = 1.6, 0.0, 0.42
    for k in range(5):
        ax.add_patch(plt.Rectangle((ends[k], y1), ends[k + 1] - ends[k], hgt, color=pal[k], ec='white', lw=0.8))
    for k in (3, 4):
        ax.text((ends[k] + ends[k + 1]) / 2, y1 + hgt / 2, f'rung, hit {hits_[k]}', ha='center', va='center', fontsize=8,
                color=INK)
    ax.text(3 / 58, y1 + hgt + 0.08, f'hits {hits_[0]}, {hits_[1]}, {hits_[2]}', ha='center', fontsize=7, color=INK)
    ends2 = [0, 1 / 3, 2 / 3, 1]
    hits2 = [c[2] for c in C2]
    lab2 = ['bottom cell', 'rung', 'rung']
    for k in range(3):
        ax.add_patch(plt.Rectangle((ends2[k], y2), ends2[k + 1] - ends2[k], hgt, color=pal[k], ec='white', lw=0.8))
        ax.text((ends2[k] + ends2[k + 1]) / 2, y2 + hgt / 2, f'{lab2[k]}, hit {hits2[k]}', ha='center', va='center',
                fontsize=8, color='white' if k == 0 else INK)
    for x0, x1 in ((0, 0), (3 / 29, 1)):
        ax.plot([x0, x1], [y1, y2 + hgt], color=MUTED, lw=0.7, ls='--')
    ax.text(-0.02, y1 + hgt / 2, r'$G = \frac{45}{116}$', ha='right', va='center', fontsize=10)
    ax.text(-0.02, y2 + hgt / 2, r"$G'' = \frac{5}{12}$", ha='right', va='center', fontsize=10)
    for x, lab in ((0.0, r'$0$'), (3 / 29, r'$\varepsilon_1 = \frac{3}{29}$'), (16 / 29, r'$1 - \vartheta$'), (1.0, r'$1$')):
        ax.text(x, y1 - 0.12, lab, ha='center', va='top', fontsize=8)
    ax.text(0.5, y2 - 0.1, r"$s'' = s/\varepsilon_1$", ha='center', va='top', fontsize=8)
    ax.set_xlim(-0.16, 1.03)
    ax.set_ylim(-0.35, 2.25)
    ax.axis('off')
    print('fig4: cells G', [(round(a_, 4), round(b_, 4), h) for a_, b_, h in C],
          'G2', [(round(a_, 4), round(b_, 4), h) for a_, b_, h in C2])
    save(fig, 4)


if __name__ == '__main__':
    fig1()
    print('ALL OK')
