"""The three figures of the polygon note, drawn from the exact classifier -> polygon/figs/.

  Figure 1  branch partitions, coloured by kind, with each band's / pair's midpoint path:
            (a) (2, 3, 4), every side (one class, N = 9 odd); (b) Example 5.4's (1, 2, 9), the class of V0V1, V1V2.
  Figure 2  slots and prongs of (2, 3, 4): (a) the N = 9 slot directions; (b) each prong to its first perpendicular
            collision -- the singular points (Theorem 4.2's bijection).
  Figure 3  Proposition 8.1's fan at b = 10 (Delta = (3, 10, 20), N = 33): the chord from V_2 to an image of V_1 across
            n = 6 copies about V_0, and the folded prong in Delta.
Branches: `beams.classify` (exact).  Paths: traced in floating point here, for drawing only.
CONTROLS (asserted, so a figure cannot silently disagree with the statement it illustrates):
  Figure 1 -- each midpoint path's head-on side / return point matches the classifier's kind and partner, and each
              band / pair label lands on exactly two branches;
  Figure 2 -- the singular points from the prongs = the classifier's branch endpoints, side by side (Theorem 4.2);
  Figure 3 -- the folded prong reaches V_1 after exactly 6 reflections with no perpendicular collision.
"""
import math

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams["pdf.fonttype"] = 42                  # TrueType fonts in the PDFs
import numpy as np

import common  # noqa: F401  (puts lib/ on the path)
from beams import classify

SIDES = [(0, 1), (1, 2), (2, 0)]
NAMES = ['V0V1', 'V1V2', 'V2V0']
BLUE, ORANGE, AQUA, VIOLET = '#2a78d6', '#eb6834', '#1baf7a', '#6250d6'
INK, MUTED, FAINT = '#222222', '#8a8a85', '#d6d5cf'
KIND_COL = {'retr': BLUE, 'transl': ORANGE, 'mirror': AQUA}
KIND_LBL = {'retr': 'R', 'transl': 'T', 'mirror': 'M'}
TOL = 1e-9


def triangle(abc):
    N = sum(abc)
    A, B, C = (math.pi * a / N for a in abc)
    return [np.array([0.0, 0.0]), np.array([math.sin(C), 0.0]), math.sin(B) * np.array([math.cos(A), math.sin(A)])]


def inward(V, s):
    p, q = V[SIDES[s][0]], V[SIDES[s][1]]
    d = (q - p) / np.linalg.norm(q - p)
    n = np.array([-d[1], d[0]])
    return n if np.dot(sum(V) / 3 - p, n) > 0 else -n


def point(V, s, t):
    p, q = V[SIDES[s][0]], V[SIDES[s][1]]
    return p + t * (q - p)


def trace(V, x, v, skip, maxb=400):
    """billiard path from x in direction v, `skip` = sides not to hit first. Stops at the first
    perpendicular collision ('perp', side, t) or vertex ('vertex', index). Returns (points, end, word)."""
    pts, word = [x.copy()], []
    for _ in range(maxb):
        best = None
        for s, (i, k) in enumerate(SIDES):
            if s in skip:
                continue
            p, e = V[i], V[k] - V[i]
            M = np.array([v, -e]).T
            if abs(np.linalg.det(M)) < 1e-14:
                continue
            tt, u = np.linalg.solve(M, p - x)
            if tt > 1e-11 and -TOL <= u <= 1 + TOL and (best is None or tt < best[0]):
                best = (tt, s, u)
        tt, s, u = best
        x = x + tt * v
        pts.append(x.copy())
        if u < 1e-7 or u > 1 - 1e-7:
            return np.array(pts), ('vertex', SIDES[s][0] if u < 0.5 else SIDES[s][1]), word
        d = (V[SIDES[s][1]] - V[SIDES[s][0]])
        d = d / np.linalg.norm(d)
        if abs(np.dot(v, d)) < 1e-9:
            return np.array(pts), ('perp', s, u), word
        word.append(s)
        v = 2 * np.dot(v, d) * d - v
        skip = {s}
    raise RuntimeError('cap')


def launch(V, s, t):
    return trace(V, point(V, s, t), inward(V, s), {s})


def draw_tri(ax, V, fill=None, lw=1.1, labels=True):
    P = np.array(V + [V[0]])
    if fill:
        ax.fill(P[:, 0], P[:, 1], color=fill, lw=0, zorder=0)
    ax.plot(P[:, 0], P[:, 1], color=INK, lw=lw, solid_joinstyle='round', zorder=2)
    if labels:
        cen = sum(V) / 3
        for i, v in enumerate(V):
            d = (v - cen) / np.linalg.norm(v - cen)
            ax.text(*(v + 0.05 * d), f'$V_{i}$', color=MUTED, fontsize=9, ha='center', va='center')
    ax.set_aspect('equal')
    ax.axis('off')
    if labels:
        P = np.array(V)
        ax.set_xlim(P[:, 0].min() - 0.1, P[:, 0].max() + 0.1)
        ax.set_ylim(P[:, 1].min() - 0.1, P[:, 1].max() + 0.1)


def side_branch(ax, V, s, t0, t1, col, off=0.022, lw=4.5, gap=0.006):
    p, q = V[SIDES[s][0]], V[SIDES[s][1]]
    L = np.linalg.norm(q - p)
    a, b = point(V, s, t0 + gap / L), point(V, s, t1 - gap / L)
    o = -inward(V, s) * off
    ax.plot([a[0] + o[0], b[0] + o[0]], [a[1] + o[1], b[1] + o[1]], color=col, lw=lw,
            solid_capstyle='butt', zorder=3)
    return (a + b) / 2 + o


def label_branch(ax, V, s, mid, text, off=0.05):
    o = -inward(V, s) * off
    ax.text(*(mid + o), text, fontsize=7.5, color=INK, ha='center', va='center')


def fig1():
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.7), gridspec_kw={'width_ratios': [1, 1.15]})
    for ax, abc, sides, title in [
            (axes[0], (2, 3, 4), (0, 1, 2), '(a)  $(2\\pi/9,\\ \\pi/3,\\ 4\\pi/9)$, every side'),
            (axes[1], (1, 2, 9), (0, 1), '(b)  $(\\pi/12,\\ \\pi/6,\\ 3\\pi/4)$, the class of $V_0V_1$, $V_1V_2$')]:
        V = triangle(abc)
        draw_tri(ax, V)
        br = {}
        for s in sides:
            r = classify(list(abc), s)
            assert r['ok']
            br[s] = r['branches']

        def which(s, t):
            for k, (t0, w, *_rest) in enumerate(br[s]):
                if t0 - 1e-9 <= t <= t0 + w + 1e-9:
                    return k

        # number each band / pair; both of its branches carry the label, so the pairing reads at a glance
        ids, count = {}, {'retr': 0, 'transl': 0, 'mirror': 0}
        for s in sides:
            for k, (t0, w, kind, tau, eps) in enumerate(br[s]):
                pts, end, word = launch(V, s, t0 + w / 2)
                assert end[0] == 'perp' and end[1] == tau, (abc, s, k, end, tau)
                if kind != 'retr':
                    assert abs((end[2]) - (t0 + w / 2)) > 1e-6        # a pair: lands elsewhere on s
                if (s, k) in ids:
                    continue
                partner = (end[1], which(end[1], end[2]))
                assert partner[1] is not None and br[partner[0]][partner[1]][2] == kind
                count[kind] += 1
                ids[(s, k)] = ids[partner] = f'{KIND_LBL[kind]}$_{count[kind]}$'
                # one path per band / pair (a band's midpoint lands on its partner's midpoint)
                ax.plot(pts[:, 0], pts[:, 1], color=KIND_COL[kind], lw=0.9, alpha=0.75, zorder=1)
                if kind == 'retr':
                    ax.plot(*pts[-1], marker='o', ms=3.5, color=KIND_COL[kind], mec='white', mew=0.8, zorder=4)
        labels = list(ids.values())
        assert all(labels.count(x) == 2 for x in labels), labels
        for s in sides:
            for k, (t0, w, kind, tau, eps) in enumerate(br[s]):
                mid = side_branch(ax, V, s, t0, t0 + w, KIND_COL[kind])
                label_branch(ax, V, s, mid, ids[(s, k)])
        ax.set_title(title, fontsize=9, color=INK, loc='left', pad=16)
    fig.text(0.01, 0.02, 'R$_i$  the two ends of retracing band $i$ (blue; the dot marks the head-on collision)      '
             'T$_i$  the two branches of translation pair $i$ (orange)',
             fontsize=7.5, color=MUTED)
    fig.subplots_adjust(left=0.02, right=0.99, top=0.9, bottom=0.08, wspace=0.08)
    return fig


def slot_dirs(abc):
    """N odd: the class C = all half-integer multiples of pi/N (side V0V1 horizontal). Returns per vertex
    the unit directions strictly inside the angle."""
    N = sum(abc)
    assert N % 2
    A, B = (math.pi * a / N for a in abc[:2])
    lo = {0: 0.0, 1: math.pi - B, 2: math.pi + A}        # angle of the first side direction at V_i
    out = {}
    for i, a in enumerate(abc):
        out[i] = [lo[i] + (j + 0.5) * math.pi / N for j in range(a)]
    return out


def fig2():
    abc = (2, 3, 4)
    V = triangle(abc)
    fig, axes = plt.subplots(1, 2, figsize=(8.6, 3.7))
    dirs = slot_dirs(abc)
    assert sum(len(d) for d in dirs.values()) == sum(abc)
    vcol = [BLUE, ORANGE, AQUA]
    ax = axes[0]
    draw_tri(ax, V)
    for i, ds in dirs.items():
        for th in ds:
            u = np.array([math.cos(th), math.sin(th)])
            ax.annotate('', xy=V[i] + 0.16 * u, xytext=V[i],
                        arrowprops=dict(arrowstyle='-|>', color=vcol[i], lw=1.2, mutation_scale=7))
    ax.set_title('(a)  the $N = 9$ slots: directions of $C$ inside each angle', fontsize=9, color=INK, loc='left', pad=16)
    ax = axes[1]
    draw_tri(ax, V)
    sing = {0: [], 1: [], 2: []}
    nverts = 0
    for i, ds in dirs.items():
        skip = {s for s, e in enumerate(SIDES) if i in e}
        for th in ds:
            u = np.array([math.cos(th), math.sin(th)])
            pts, end, word = trace(V, V[i].copy(), u, skip)
            if end[0] == 'vertex':
                nverts += 1
                continue
            sing[end[1]].append(end[2])
            ax.plot(pts[:, 0], pts[:, 1], color=vcol[i], lw=0.9, alpha=0.85, zorder=1)
            ax.plot(*pts[-1], marker='o', ms=4.5, color=vcol[i], mec='white', mew=0.9, zorder=4)
    assert nverts == 0                                  # (2,3,4): p = s' = 0
    for s in range(3):                                   # CONTROL: Theorem 4.2's bijection, side by side
        r = classify(list(abc), s)
        ends = sorted(b[0] for b in r['branches'])[1:]
        got = sorted(sing[s])
        assert len(got) == len(ends) and all(abs(a - b) < 1e-7 for a, b in zip(got, ends)), (s, got, ends)
    ax.set_title('(b)  each prong to its perpendicular collision: the singular points', fontsize=9, color=INK,
                 loc='left', pad=16)
    fig.text(0.01, 0.02, 'colour = the vertex of the slot (blue $V_0$, orange $V_1$, aqua $V_2$); '
             '$2 + 3 + 4 = 9$ prongs, $9 = \\Sigma_\\sigma (n_\\sigma - 1)$ singular points', fontsize=7.5, color=MUTED)
    fig.subplots_adjust(left=0.02, right=0.99, top=0.9, bottom=0.08, wspace=0.08)
    return fig, {s: len(v) for s, v in sing.items()}


def fig3(b=10):
    abc = (3, b, 2 * b)
    N = sum(abc)
    A = math.pi * 3 / N
    n = 2 * (b - 1) // 3
    V = triangle(abc)
    rC, rB = np.linalg.norm(V[1]), np.linalg.norm(V[2])   # |V0V1|, |V0V2|
    ray = lambda k: (rC if k % 2 == 0 else rB) * np.array([math.cos(k * A), math.sin(k * A)])
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.9), gridspec_kw={'width_ratios': [0.8, 1]})
    ax = axes[0]
    for k in range(-n, 1):
        tri = [np.zeros(2), ray(k), ray(k + 1)]
        draw_tri(ax, tri, fill='#eeeeea' if k == 0 else None, lw=0.6, labels=False)
    P, Q = ray(1), ray(-n)
    ax.plot([P[0], Q[0]], [P[1], Q[1]], color=BLUE, lw=1.8, zorder=3)
    for z, lbl, dd in [(np.zeros(2), '$V_0$', (-0.05, 0)), (P, '$Y = V_2$', (0.0, 0.05)),
                       (Q, '$Z$', (0.03, -0.04)), (ray(0), '$V_1$', (0.05, 0.0))]:
        ax.text(z[0] + dd[0], z[1] + dd[1], lbl, fontsize=9, color=INK, ha='center', va='center')
    ax.set_title(f'(a)  the fan about $V_0$, $b = {b}$: $\\Delta = ({3}, {b}, {2 * b})$, $N = {N}$', fontsize=9,
                 color=INK, loc='left', pad=16)
    ax = axes[1]
    draw_tri(ax, V, fill='#eeeeea')
    u = (Q - P) / np.linalg.norm(Q - P)
    pts, end, word = trace(V, V[2].copy(), u, {1, 2})
    assert end == ('vertex', 1) and len(word) == n and set(word) <= {0, 2}, (end, word)
    ax.plot(pts[:, 0], pts[:, 1], color=BLUE, lw=1.4, zorder=3)
    ax.set_title(f'(b)  folded: the prong $V_2 \\to V_1$, {n} reflections in $V_0V_1$, $V_2V_0$', fontsize=9,
                 color=INK, loc='left', pad=16)
    fig.text(0.01, 0.02, 'Shaded: $\\Delta$. The chord $YZ$ crosses each ray of the fan inside its side copy '
             '(Proposition 8.1), so its fold is a vertex-to-vertex path meeting no side at a right angle.',
             fontsize=7.5, color=MUTED)
    fig.subplots_adjust(left=0.02, right=0.99, top=0.9, bottom=0.08, wspace=0.05)
    return fig, word


def main():
    import os
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'figs')
    os.makedirs(out, exist_ok=True)
    f = fig1()
    for ext in ('png', 'pdf'):
        f.savefig(os.path.join(out, f'fig1_partitions.{ext}'), dpi=200)
    f, sing = fig2()
    for ext in ('png', 'pdf'):
        f.savefig(os.path.join(out, f'fig2_prongs.{ext}'), dpi=200)
    f, word = fig3()
    for ext in ('png', 'pdf'):
        f.savefig(os.path.join(out, f'fig3_fan.{ext}'), dpi=200)
    print('fig2 singular points per side:', sing, '| fig3 word:', ''.join(map(str, word)))


if __name__ == '__main__':
    main()
