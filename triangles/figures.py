"""[Deposit] Figures 1-3. Paper: "Complete periodicity of perpendicular billiards in rational triangles of prime denominator". Research name: probes/s949_tri_figures.py, copied verbatim except imports (deposit names), paths (data/ and figs/ under triangles/) and this paragraph. Run through run_all.py; the usage lines below are the research repo's. Bracketed IDs and section marks (§sNNN) refer to the project's internal ledgers.

s949 (`s911-a`): the three figures of paper 1 (`tex/tri/tri.tex`).

  Fig 1  (§1) the two cases of Theorem A against a counterexample: (a) a prong of a prime-N triangle that meets a
         side at a right angle and retraces to its vertex; (b) the swapped prong of (3,4,8), N = 15 (the note's
         Prop 8.1, F1 at b = 4): V2 -> V1 in 2 reflections, no right angle.
  Fig 2  (§4, Prop 4.2 step (iii)) schematic: the four germs at P_u, the sectors in the components A, B, the two
         admissible cyclic orders O_1, O_2, exchanged by the orientation-reversing iota.
  Fig 3  (§6.1, App A) the band prong W_{3,2,0} of (3,14,28), N = 45 (Theorem C's Delta_1 at m = 2, t = 2/3),
         developed straight across its three stars (fans about the images of V0).
CONTROLS (asserted): Fig 1(a) -- the path meets a side perpendicularly before any vertex; (a) is chosen as the
  shortest such prong over all scalene triangles at N = 7, 11, 13 and EVERY slot there is checked to retrace
  (Theorem A at those N, a float check); Fig 1(b) -- reaches V1 after exactly 2 reflections, no perpendicular
  collision; Fig 2 -- the two orders are exactly those consistent with the sector labels and iota swaps them
  (brute force over the 6 cyclic orders); Fig 3 -- the traced word equals W_{3,2,0} = 0|(20)^6|(02)^3 0, it ends at
  an image of V2, and its length equals lambda_0 of Prop 6.2.
Paths are float-traced, for drawing only. Not pre-registered: figures. PRIOR ART: probes/s756_figures.py (the note's
figures, same style). DATA: none.
Output: figs/fig{1,2,3}.pdf (+ .png).
Run from the repo root:
PYTHONPATH=.:engine:archive/scripts_2026-07:archive/scripts_2026-08:probes .venv/bin/python3.13 probes/s949_tri_figures.py
"""
import _setup  # noqa: F401  (deposit: lib/ on the path, cwd = triangles/)
import itertools
import math
import shutil

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Wedge, Polygon, FancyArrowPatch
plt.rcParams["pdf.fonttype"] = 42                  # TrueType: Hexagon safe-PDF (s900, s909)
plt.rcParams["mathtext.fontset"] = "cm"
import numpy as np

BLUE, ORANGE, AQUA, VIOLET = '#2a78d6', '#eb6834', '#1baf7a', '#6250d6'
INK, MUTED, FAINT = '#222222', '#8a8a85', '#d6d5cf'
SIDES = [(0, 1), (1, 2), (2, 0)]                   # side j runs V_j -> V_{j+1}: 0 = V0V1, 1 = V1V2, 2 = V2V0
TOL = 1e-9


def triangle(abc):
    """circumdiameter 1: V0 = 0, V1 = sin C on the x-axis, V2 = sin B e^{iA}."""
    N = sum(abc)
    A, B, C = (math.pi * a / N for a in abc)
    return [np.array([0.0, 0.0]), np.array([math.sin(C), 0.0]), math.sin(B) * np.array([math.cos(A), math.sin(A)])]


def hit(x, v, p, q):
    """parameter s > 0 along x + s v where it meets segment pq, and the position u in [0,1] on pq (or None)."""
    d = q - p
    M = np.array([[v[0], -d[0]], [v[1], -d[1]]])
    if abs(np.linalg.det(M)) < 1e-15:
        return None
    s, u = np.linalg.solve(M, p - x)
    return (s, u) if s > TOL and -TOL <= u <= 1 + TOL else None


def trace(V, k, ang, maxb=200):
    """billiard prong from vertex V_k in direction ang. Returns (points, end, word):
    end = ('perp', side) at the first right-angle collision, or ('vertex', index)."""
    x, v = V[k].copy(), np.array([math.cos(ang), math.sin(ang)])
    pts, word, last = [x.copy()], [], None
    for _ in range(maxb):
        best = None
        for j, (a, b) in enumerate(SIDES):
            if j == last or (last is None and k in (a, b)):
                continue
            h = hit(x, v, V[a], V[b])
            if h and (best is None or h[0] < best[0]):
                best = (h[0], h[1], j)
        s, u, j = best
        x = x + s * v
        pts.append(x.copy())
        a, b = SIDES[j]
        if u < 1e-7 or u > 1 - 1e-7:
            return pts, ('vertex', a if u < 0.5 else b), word
        t = (V[b] - V[a]) / np.linalg.norm(V[b] - V[a])
        word.append(j)
        if abs(np.dot(v, t)) < 1e-9:
            return pts, ('perp', j), word
        v = 2 * np.dot(v, t) * t - v
        last = j
    return pts, ('capped', None), word


def slots(abc):
    """(vertex, direction) with direction an odd multiple of pi/2N strictly inside the vertex angle."""
    N = sum(abc)
    V = triangle(abc)
    g = sum(V) / 3
    out = []
    for k in range(3):
        for i in range(4 * N):
            if i % 2 == 0:
                continue
            ang = i * math.pi / (2 * N)
            u = np.array([math.cos(ang), math.sin(ang)])
            a, b = [V[j] - V[k] for j in range(3) if j != k]
            cr = lambda p, q: p[0] * q[1] - p[1] * q[0]
            sab = cr(a, b)
            if cr(a, u) * sab > 1e-12 and cr(u, b) * sab > 1e-12:   # u strictly inside the (convex) angle
                out.append((k, ang))
    assert len(out) == N, (abc, len(out))           # S = N slots at odd N ([Poly, Lemma 4.1])
    return out


def is_prime(n):
    return n > 1 and all(n % d for d in range(2, int(n ** .5) + 1))


def draw_tri(ax, V, labels=True, off=0.06):
    ax.add_patch(Polygon(V, closed=True, fill=True, fc='#f4f3ef', ec=INK, lw=1.0, zorder=1))
    if labels:
        g = sum(V) / 3
        for i, P in enumerate(V):
            d = (P - g) / np.linalg.norm(P - g)
            ax.text(*(P + off * d), f'$V_{i}$', ha='center', va='center', fontsize=10, color=INK)


def right_angle_mark(ax, P, side_dir, path_dir, size=0.035):
    t = side_dir / np.linalg.norm(side_dir)
    n = -path_dir / np.linalg.norm(path_dir)
    ax.plot(*np.array([P + size * t, P + size * (t + n), P + size * n]).T, color=INK, lw=0.8, zorder=5)


# ------------------------------------------------------------------------------------------------- Fig 1
def fig1():
    # (a) shortest retracing prong over prime N = 7, 11, 13; check every slot retraces (Theorem A, float)
    cands = []
    for N in (7, 11, 13):
        for a in range(1, N):
            for b in range(a + 1, N - a):
                c = N - a - b
                if c <= b or math.gcd(math.gcd(a, b), c) != 1:
                    continue
                abc = (a, b, c)
                V = triangle(abc)
                for k, ang in slots(abc):
                    pts, end, word = trace(V, k, ang)
                    assert end[0] == 'perp', (abc, k, ang, end)       # every prong retraces at prime N
                    cands.append((len(word), abs(math.pi * c / N - math.pi / 3) + abs(math.pi * a / N - math.pi / 4), abc, k, ang))
    cands = [c for c in cands if c[0] >= 3]
    cands.sort(key=lambda c: (c[0], c[1]))
    _, _, abcA, kA, angA = cands[0]
    VA = triangle(abcA)
    ptsA, endA, wordA = trace(VA, kA, angA)
    # (b) (3,4,8): the swapped prong V2 -> V1
    abcB = (3, 4, 8)
    VB = triangle(abcB)
    sw = []
    for k, ang in slots(abcB):
        pts, end, word = trace(VB, k, ang)
        if end[0] == 'vertex':
            sw.append((k, ang, pts, end, word))
    assert len(sw) == 2, len(sw)                    # p = 1: one swapped pair = 2 slots
    k, ang, ptsB, endB, wordB = next(s for s in sw if s[0] == 2)
    assert endB == ('vertex', 1) and len(wordB) == 2 and set(wordB) <= {0, 2}, (endB, wordB)

    fig, axs = plt.subplots(1, 2, figsize=(6.4, 3.0))
    a0, a1, a2 = abcA
    for ax, V, pts, abc, ttl in ((axs[0], VA, ptsA, abcA, rf'$\mathrm{{(a)}}\ \ ({a0}, {a1}, {a2}),\ N = {sum(abcA)}$'),
                                 (axs[1], VB, ptsB, abcB, r'$\mathrm{(b)}\ \ (3, 4, 8),\ N = 15$')):
        draw_tri(ax, V)
        P = np.array(pts)
        ax.plot(P[:, 0], P[:, 1], color=BLUE if ax is axs[0] else ORANGE, lw=1.6, zorder=3)
        ax.plot(*P[0], 'o', color=INK, ms=4, zorder=4)
        ax.text(0.5, -0.04, ttl, transform=ax.transAxes, ha='center', va='top', fontsize=10)
        ax.set_aspect('equal'); ax.axis('off')
    # mark the right angle in (a)
    j = endA[1]
    a, b = SIDES[j]
    right_angle_mark(axs[0], ptsA[-1], VA[b] - VA[a], ptsA[-1] - ptsA[-2])
    axs[0].plot(*ptsA[-1], 's', color=BLUE, ms=3, zorder=4)
    axs[1].plot(*ptsB[-1], 'o', color=INK, ms=4, zorder=4)
    fig.tight_layout()
    return fig, dict(a=(abcA, kA, len(wordA)), b=(abcB, wordB))


# ------------------------------------------------------------------------------------------------- Fig 2
CCW = {'o': 'A', 'i': 'B', "o'": 'B', "i'": 'A'}     # sector just counterclockwise of each germ
CW = {g: ('B' if c == 'A' else 'A') for g, c in CCW.items()}
IOTA = {'o': "i'", 'i': "o'", "o'": 'i', "i'": 'o'}


def admissible():
    res = []
    for perm in itertools.permutations(['i', "o'", "i'"]):
        order = ('o',) + perm
        if all(CCW[order[n]] == CW[order[(n + 1) % 4]] for n in range(4)):
            res.append(order)
    return res


def iota_image(order):
    img = [IOTA[g] for g in order][::-1]             # orientation-reversing: ccw order of images is reversed
    r = img.index('o')
    return tuple(img[r:] + img[:r])


def fig2():
    adm = admissible()
    O1, O2 = ('o', 'i', "i'", "o'"), ('o', "o'", "i'", 'i')
    assert sorted(adm) == sorted([O1, O2]), adm
    assert iota_image(O1) == O2 and iota_image(O2) == O1
    lab = {'o': r"$\gamma_{\mathrm{out}}$", 'i': r"$\gamma_{\mathrm{in}}$",
           "o'": r"$\gamma'_{\mathrm{out}}$", "i'": r"$\gamma'_{\mathrm{in}}$"}
    col = {'A': BLUE, 'B': ORANGE}
    fig, axs = plt.subplots(1, 3, figsize=(6.6, 2.6), gridspec_kw=dict(width_ratios=[1, 0.42, 1]))
    for ax, order, name in ((axs[0], O1, r'$\mathcal{O}_1$'), (axs[2], O2, r'$\mathcal{O}_2$')):
        angs = [45 + 90 * n for n in range(4)]
        for n, g in enumerate(order):
            a0, a1 = angs[n], angs[n] + 90
            ax.add_patch(Wedge((0, 0), 1.0, a0, a1, fc=col[CCW[g]], alpha=0.16, ec='none'))
            mid = math.radians(a0 + 45)
            ax.text(0.62 * math.cos(mid), 0.62 * math.sin(mid), rf'$\mathbf{{{CCW[g]}}}$', color=col[CCW[g]],
                    ha='center', va='center', fontsize=12)
        for n, g in enumerate(order):
            t = math.radians(angs[n])
            e = np.array([math.cos(t), math.sin(t)])
            c = INK if g in ('o', 'i') else VIOLET
            if g in ('o', "o'"):                     # out-germ: arrow away from P_u
                ax.add_patch(FancyArrowPatch((0, 0), 1.0 * e, arrowstyle='-|>', mutation_scale=10, color=c, lw=1.4))
            else:                                    # in-germ: arrow towards P_u
                ax.add_patch(FancyArrowPatch(1.0 * e, 0.06 * e, arrowstyle='-|>', mutation_scale=10, color=c, lw=1.4))
            ax.text(*(1.25 * e), lab[g], ha='center', va='center', fontsize=10, color=c)
        ax.add_patch(FancyArrowPatch((1.32 * math.cos(math.radians(-20)), 1.32 * math.sin(math.radians(-20))),
                                     (1.32 * math.cos(math.radians(20)), 1.32 * math.sin(math.radians(20))),
                                     connectionstyle='arc3,rad=0.35', arrowstyle='-|>', mutation_scale=9, color=MUTED, lw=0.9))
        ax.plot(0, 0, 'o', color=INK, ms=4)
        ax.text(0.13, 0.0, '$P_u$', fontsize=9, color=INK, ha='left', va='center')
        ax.set_title(name, fontsize=11)
        ax.set_xlim(-1.55, 1.55); ax.set_ylim(-1.45, 1.45)
        ax.set_aspect('equal'); ax.axis('off')
    ax = axs[1]
    ax.add_patch(FancyArrowPatch((0.05, 0.12), (0.95, 0.12), arrowstyle='<|-|>', mutation_scale=12, color=INK, lw=1.2))
    ax.text(0.5, 0.25, r'$\iota$', ha='center', va='bottom', fontsize=13)
    ax.set_xlim(0, 1); ax.set_ylim(-0.6, 0.8); ax.axis('off')
    fig.tight_layout()
    return fig, dict(orders=adm)


# ------------------------------------------------------------------------------------------------- Fig 3
def reflect(P, a, b):
    d = (b - a) / np.linalg.norm(b - a)
    w = P - a
    return a + 2 * np.dot(w, d) * d - w


def develop(V, k, ang, maxl=400):
    """straight development of the prong from V_k: list of copies (3 points each) and the word."""
    T = [p.copy() for p in V]
    x, v = V[k].copy(), np.array([math.cos(ang), math.sin(ang)])
    copies, word, last = [list(T)], [], None
    for _ in range(maxl):
        best = None
        for j, (a, b) in enumerate(SIDES):
            if j == last or (last is None and k in (a, b)):
                continue
            h = hit(x, v, T[a], T[b])
            if h and (best is None or h[0] < best[0]):
                best = (h[0], h[1], j)
        s, u, j = best
        a, b = SIDES[j]
        if u < 1e-7 or u > 1 - 1e-7:
            end = a if u < 0.5 else b
            return copies, word, x + s * v, end
        word.append(j)
        c = 3 - a - b
        T = list(T)
        T[c] = reflect(T[c], T[a], T[b])
        copies.append(T)
        last = j
    raise RuntimeError('not closed')


def fig3():
    q, m, j, t = 3, 2, 0, 2 / 3
    abc = (3, 14, 28)
    N = sum(abc); n = q * (2 * m + 1)
    A = math.pi / n; B = math.pi * abc[1] / N; beta = math.pi / q
    assert abs(B - (2 * m + t) * A) < 1e-12           # K = 2m at q = 3
    V = triangle(abc)
    psi = -math.pi / 2 + (t - j) * A
    copies, word, E, end = develop(V, 2, psi)
    W = [0, 1] + [2, 0] * 6 + [1] + [0, 2] * 3 + [0]
    assert word == W, word
    assert end == 2                                   # ends at an image of V2
    lam = 2 * math.sin(B) * math.cos(t * A) / math.sin(beta / 2) * math.cos((2 * m - 1 - 2 * j) * A / 2)
    assert abs(np.linalg.norm(E - V[2]) - lam) < 1e-9, (np.linalg.norm(E - V[2]), lam)
    rot = -psi
    Rm = np.array([[math.cos(rot), -math.sin(rot)], [math.sin(rot), math.cos(rot)]])
    R = lambda P: Rm @ (P - V[2])
    copies = [[R(P) for P in T] for T in copies]
    E, S0 = R(E), R(V[2])
    fig, ax = plt.subplots(figsize=(6.6, 2.6))
    for i, T in enumerate(copies):
        ax.add_patch(Polygon(T, closed=True, fc='#f4f3ef' if i else '#e3e8f3', ec=FAINT if i else INK,
                             lw=0.6 if i else 1.0, zorder=1))
    # star centres: images of V0, in order of appearance
    cents = []
    for T in copies:
        if not cents or np.linalg.norm(T[0] - cents[-1]) > 1e-9:
            cents.append(T[0])
    assert len(cents) == q, len(cents)
    for i, c in enumerate(cents):
        ax.plot(*c, 'o', color=AQUA, ms=5, zorder=4)
        dy = 0.07 if c[1] > 0 else -0.07
        ax.text(c[0], c[1] + dy, f'$O_{i}$', ha='center', va='bottom' if c[1] > 0 else 'top', fontsize=9, color=AQUA)
    ax.plot([S0[0], E[0]], [S0[1], E[1]], color=ORANGE, lw=1.6, zorder=3)
    ax.plot(*S0, 'o', color=INK, ms=4, zorder=5)
    ax.plot(*E, 'o', color=INK, ms=4, zorder=5)
    ax.text(S0[0] - 0.05, S0[1], '$V_2$', ha='right', va='center', fontsize=10)
    ax.text(E[0] + 0.05, E[1], '$E$', ha='left', va='center', fontsize=10)
    ax.set_aspect('equal'); ax.axis('off')
    fig.tight_layout()
    return fig, dict(word_len=len(word), copies=len(copies), lam=lam)


if __name__ == '__main__':
    import os
    os.makedirs('figs', exist_ok=True)
    for name, f in (('fig1', fig1), ('fig2', fig2), ('fig3', fig3)):
        fig, info = f()
        for ext in ('pdf', 'png'):
            fig.savefig(f'figs/{name}.{ext}', dpi=200, bbox_inches='tight')
        print(name, info)
    print('ALL OK')
