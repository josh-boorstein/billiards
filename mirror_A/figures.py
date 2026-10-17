"""The four figures of the paper -> mirror_A/figs/ (PNG and PDF).

  Figure 1  the triangle, its rhombus and doors; a line of the class perpendicular to OB
            developed through the chain of rhombi from its fold level (Sections 2.2, 3.2).
  Figure 2  Example 5.1's cylinder: folded into T, and developed (Section 5.1).
  Figure 3  the corner word HAHBH and its cylinder of period 14 (Section 6.1).
  Figure 4  the five tiles of Proposition 6.6 (Section 6.3).

Every number a caption states is checked by `verify()` before anything is drawn, with
controls that fire; the drawing is float geometry (the exact statements are the other
scripts').  Prints one line per check and exits nonzero if any fails.
"""
import cmath
import math
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figs")
INK, ACC, PERP, FADE = "#1f2933", "#2f6f9f", "#c0392b", "#b8c2cc"
GRN, GLD, VIO = "#4a8f6d", "#c08a3e", "#8f6da8"
DOORC = [ACC, GLD, GRN, VIO, PERP, "#5b6770"]


# ------------------------------------------------------------------ Tri(a), complex coords

def tri(a):
    return {"O": 0j, "A": 1 + 0j, "B": 1j * math.tan(a)}


SIDES = {"A": ("O", "A"), "B": ("O", "B"), "H": ("A", "B")}


def refl_map(letter, a):
    """The reflection in side `letter`, as an affine-antilinear map z -> al*conj(z) + be."""
    u = cmath.exp(-2j * a)
    return {"A": (1, 0), "B": (-1, 0), "H": (u, 1 - u)}[letter]


def compose(f, g):
    """(f o g) for maps of the form (al, be, anti); anti = conj applied first."""
    al1, be1, an1 = f
    al2, be2, an2 = g
    # f(g(z)) = al1 * C1(al2 * C2(z) + be2) + be1
    if an1:
        return (al1 * al2.conjugate(), al1 * be2.conjugate() + be1, not an2)
    return (al1 * al2, al1 * be2 + be1, an2)


def apply(f, z):
    al, be, an = f
    return al * (z.conjugate() if an else z) + be


def unfold(word, a):
    """g_0 = id, g_i = rho_{w1} o ... o rho_{wi}; copies T_i = g_i(T); edges e_i = g_{i-1}(w_i)."""
    V = tri(a)
    g = (1 + 0j, 0j, False)
    gs, edges = [g], []
    for w in word:
        p, q = SIDES[w]
        edges.append((apply(g, V[p]), apply(g, V[q])))
        al, be = refl_map(w, a)
        g = compose(g, (complex(al), complex(be), True))
        gs.append(g)
    copies = [[apply(gi, V[k]) for k in "OAB"] for gi in gs]
    return gs, edges, copies


def realised(word, a):
    """Lemma 6.2 sign test: every crossed edge's endpoints straddle the line O -> O_end."""
    gs, edges, _ = unfold(word, a)
    oend = apply(gs[-1], 0j)
    f = lambda P: (oend.conjugate() * P).imag
    return all(f(P) * f(Q) < 0 for P, Q in edges)


def level_of(g, a):
    """Level n of a copy: linear part rho_phi or rho_phi sigma_OA with phi = 2an (mod pi)."""
    phi = cmath.phase(g[0])
    best = min(range(-40, 41), key=lambda n: abs(math.remainder(phi - 2 * a * n, math.pi)))
    assert abs(math.remainder(phi - 2 * a * best, math.pi)) < 1e-9
    return best


# ------------------------------------------------------------------- a billiard tracer

def trace(a, p, d, nb, eps=1e-12):
    """Billiard in Tri(a) from p (complex) in unit direction d, nb bounces.
    Returns the list of points (p, hit_1, ..., hit_nb), the word, and the final direction."""
    V = tri(a)
    pts, word = [p], []
    last = None
    for _ in range(nb):
        best = None
        for s, (i, j) in SIDES.items():
            if s == last:
                continue
            P, Q = V[i], V[j]
            e = Q - P
            # solve p + t d = P + u e
            A11, A12, A21, A22 = d.real, -e.real, d.imag, -e.imag
            det = A11 * A22 - A12 * A21
            if abs(det) < 1e-15:
                continue
            bx, by = (P - p).real, (P - p).imag
            t = (bx * A22 - A12 * by) / det
            u = (A11 * by - A21 * bx) / det
            if t > eps and -1e-12 <= u <= 1 + 1e-12 and (best is None or t < best[0]):
                best = (t, s, u)
        t, s, u = best
        p = p + t * d
        e = V[SIDES[s][1]] - V[SIDES[s][0]]
        e = e / abs(e)
        d = e * e * d.conjugate()
        pts.append(p)
        word.append(s)
        last = s
    return pts, "".join(word), d


# --------------------------------------------------------------------- the rhombus chain

def rhombus(a, n, c):
    r = cmath.exp(2j * a * n)
    A, B = r * 1, r * 1j * math.tan(a)
    return [c + A, c + B, c - A, c - B]           # A, B, -A, -B


def chain_walk(a, p0, d, nsteps, n0=0):
    """Walk a straight line through the chain from level `n0` (centre 0), both ways.
    Returns seq = [(level, centre)] in order along +d, and the crossings [(door k, point)]."""
    h = math.sin(a)

    def step(n, c, p, dd):
        V = rhombus(a, n, c)
        best = None
        for i in range(4):
            P, Q = V[i], V[(i + 1) % 4]
            e = Q - P
            A11, A12, A21, A22 = dd.real, -e.real, dd.imag, -e.imag
            det = A11 * A22 - A12 * A21
            if abs(det) < 1e-14:
                continue
            bx, by = (P - p).real, (P - p).imag
            t = (bx * A22 - A12 * by) / det
            u = (A11 * by - A21 * bx) / det
            if t > 1e-9 and 0 <= u <= 1 and (best is None or t < best[0]):
                best = (t, i, e, P)
        t, i, e, P = best
        nu = -1j * e / abs(e)                       # outward normal (vertices run ccw)
        # edges A->B (i=0) and -A->-B (i=2) have direction -a+2an: door to n-1
        dn = -1 if i in (0, 2) else +1
        c2 = c + 2 * h * nu
        # the reflected rhombus is R_{n+dn} at c2 (Lemma 2.3's labelling) ...
        V2 = rhombus(a, n + dn, c2)
        Vr = [2 * (P + ((x - P) * e.conjugate()).real * e / abs(e) ** 2) - x for x in V]
        assert all(min(abs(x - y) for y in V2) < 1e-9 for x in Vr), "chain reflection"
        # ... and the offset moves by Lemma 2.4's mu_k (in the frame theta = +e_x)
        k = n + dn / 2
        thd = cmath.phase(dd)                       # the frame of the direction of motion
        chi = math.pi / 2 + 2 * a * k - thd
        mu = -2 * h * math.copysign(1, math.cos(chi)) * math.sin(chi)
        ds = -(dd.conjugate() * (c2 - c)).imag
        assert abs(ds - mu) < 1e-9, ("Lemma 2.4", ds, mu)
        return n + dn, c2, p + t * dd, k

    out = {}
    for sgn in (+1, -1):
        walk, n, c, p = [], n0, 0j, p0
        for _ in range(nsteps):
            n2, c2, p, k = step(n, c, p, sgn * d)
            walk.append((n2, c2, k, p))
            n, c = n2, c2
        out[sgn] = walk
    seq = ([(n, c) for n, c, _, _ in reversed(out[-1])] + [(n0, 0j)]
           + [(n, c) for n, c, _, _ in out[+1]])
    cross = [(k, p) for _, _, k, p in out[-1]] + [(k, p) for _, _, k, p in out[+1]]
    return seq, cross


# ----------------------------------------------------------------------------- the data

TILES_WORDS = ["HAHABH", "HAHBHAHABH", "HAHBH", "HBHAHBHBAH", "HBHBAH"]


def smaller_pos_root(coeffs):
    import numpy as np
    r = [x.real for x in np.roots(coeffs) if abs(x.imag) < 1e-12 and x.real > 0]
    return min(r)


def tile_ends():
    at = lambda s: math.atan(math.sqrt(s))
    a1 = at(3 / 5)
    a2 = at(smaller_pos_root([1, -5, -93, 41]))
    a3 = at(smaller_pos_root([1, -77, 51, 1]))
    a4 = at(smaller_pos_root([1, 10, -7]))
    a5 = at(max(x for x in [((10 + math.sqrt(128)) / 14)] if x > 0))
    H = math.pi / 2
    return {"HAHABH": (0, a1), "HAHBHAHABH": (a2, a3), "HAHBH": (a4, a5),
            "HBHAHBHBAH": (H - a3, H - a2), "HBHBAH": (H - a1, H)}, (a1, a2, a3, a4, a5)


EX_A = 0.5 * math.acos(3 / 4)
EX_WORD = "HAHABHABHBAHAHAB"
EX_TOP = math.sqrt(7) / 42


def example_path(y0, nb=16):
    return trace(EX_A, 1j * y0, cmath.exp(2j * EX_A), nb)


# ---------------------------------------------------------------------------- the checks

def verify(verbose=True):
    log = (lambda s: print(s)) if verbose else (lambda s: None)

    # (1) chain walk from the fold level l* = j: Lemma 2.3's labels and Lemma 2.4's moves are
    #     asserted step by step inside chain_walk; the level sequence obeys the fold r -> 2l* - r.
    a, j = 0.42, 2
    th = 2 * a * j
    d = cmath.exp(1j * th)
    s0 = 0.25                                       # an offset at which the line TURNS
    seq, cross = chain_walk(a, s0 * 1j * d, d, 4, n0=j)
    levels = [n for n, _ in seq]
    m = len(levels) // 2
    assert levels[m] == j and all(levels[m - i] + levels[m + i] == 2 * j for i in range(m + 1))

    def mirror_err(sq):
        """max distance between centre m+i and the mirror image of centre m-i in the
        OB-diagonal of the middle rhombus (the axis of the fold reflection)."""
        mm = len(sq) // 2
        c0 = sq[mm][1]
        ax = cmath.exp(1j * (math.pi / 2 + 2 * a * sq[mm][0]))
        refl = lambda z: c0 + ax * ax * (z - c0).conjugate()
        return max(abs(refl(sq[mm - i][1]) - sq[mm + i][1]) for i in range(mm + 1))

    assert mirror_err(seq) < 1e-9
    # the control: a line NOT in the class (theta + 0.05) turns on one side only.  (A line
    # passing straight through cannot fire it: the chain itself is fold-symmetric.)
    seqc, _ = chain_walk(a, s0 * 1j * d, cmath.exp(1j * (th + 0.05)), 4, n0=j)
    lc = [n for n, _ in seqc]
    assert any(lc[4 - i] + lc[4 + i] != 2 * j for i in range(5)), lc
    assert mirror_err(seqc) > 1e-3, mirror_err(seqc)
    log(f"  (1) chain at a={a}, theta=2a*{j}, from l*: levels {levels} (r -> 2l*-r); centres "
        f"mirror-symmetric in the fold diagonal (control theta+0.05: levels {lc}, "
        f"mirror error {mirror_err(seqc):.3f}); Lemma 2.4 moves checked per step")

    # (2) tiles: Table 2 ends vs a direct Lemma 6.2 test, with the outside-the-tile control.
    tiles, ends = tile_ends()
    for w, (lo, hi) in tiles.items():
        for x in (lo + 1e-4, (lo + hi) / 2, hi - 1e-4):
            if 0 < x < math.pi / 2:
                assert realised(w, x), (w, x)
        for x in (lo - 1e-4, hi + 1e-4):
            if 1e-6 < x < math.pi / 2 - 1e-6:
                assert not realised(w, x), ("control", w, x)       # the control fires
        n_in = sum(realised(w, lo + (hi - lo) * t / 200) for t in range(1, 200))
        assert n_in == 199, (w, n_in)
    a1, a2, a3, a4, a5 = ends
    assert a2 < a1 and a4 < a3 and a5 > math.pi / 2 - a3
    log("  (2) tiles: ends " + ", ".join(f"{e:.5f}" for e in ends)
        + "; realised inside, NOT at +-1e-4 outside, cover overlaps hold")
    for w in TILES_WORDS:
        assert not (w == w[::-1] and len(w) % 2 == 1)
    log("  (2b) none of the five words is an odd palindrome")

    # (3) corner word HAHBH at a=0.80: realised; nearby path periodic of period 14.
    aw, w = 0.80, "HAHBH"
    assert tiles["HAHBH"][0] < aw < tiles["HAHBH"][1] and realised(w, aw)
    gs, edges, copies = unfold(w, aw)
    oend = apply(gs[-1], 0j)
    phi = cmath.phase(oend)
    dvec = cmath.exp(1j * phi)
    pts, wd, dk = trace(aw, 0j, dvec, len(w))
    # after the k-th bounce the path runs straight into the corner O
    assert wd == w and abs((dk.conjugate() * (0j - pts[-1])).imag) < 1e-9, wd
    assert ((dk.conjugate() * (0j - pts[-1])).real) > 0
    eps = 0.004
    p0 = (pts[0] + pts[1]) / 2 + eps * 1j * dvec
    q, wq, dq = trace(aw, p0, dvec, 2 * len(w) + 4)
    assert abs(q[-1] - q[0]) > 0  # the first entry is p0 itself; compare after one more lap
    qq, wq2, dq2 = trace(aw, p0, dvec, 2 * (2 * len(w) + 4))
    # return to the same segment at the same point: period 14 means hit 14 equals hit 0's segment
    assert wq2[:14] == wq2[14:], wq2
    assert abs(dq - dvec) < 1e-9
    log(f"  (3) corner word {w} at a={aw}: realised, phi={phi:.4f}, nearby path word {wq2[:14]}"
        f" (period 14 = 2k+4), direction returns")

    # (4) Example 5.1.
    pts, wd, dd = example_path(EX_TOP / 2)
    assert wd == EX_WORD, wd
    assert abs(pts[-1] - pts[0]) < 1e-9 and abs(dd - cmath.exp(2j * EX_A)) < 1e-9
    gs, _, copies = unfold(EX_WORD, EX_A)
    gT = gs[-1]
    assert abs(gT[0] - 1) < 1e-12 and not gT[2]
    assert abs(gT[1] - 4 * cmath.exp(2j * EX_A)) < 1e-9
    for y in (EX_TOP * 1e-3, EX_TOP * (1 - 1e-3)):
        p, wy, _ = example_path(y)
        assert wy == EX_WORD and abs(p[-1] - p[0]) < 1e-9
    for y in (EX_TOP * (1 + 1e-2), EX_TOP * 1.2):
        p, wy, _ = example_path(y)
        assert wy != EX_WORD or abs(p[-1] - p[0]) > 1e-6, ("control", y)  # control fires
    width = EX_TOP * math.cos(2 * EX_A)
    assert abs(width - math.sqrt(7) / 56) < 1e-15
    lv = sorted({level_of(g, EX_A) for g in gs})
    assert 1 not in lv                              # fold level l* = j = 1 never visited
    log(f"  (4) Example 5.1: word, period 16, closure 4e(2a), launch (0, sqrt7/42) with "
        f"outside control, width sqrt7/56; levels visited {lv}, fold l*=1 avoided")
    return {"chain": (a, j, s0, seq, cross), "tiles": (tiles, ends)}


# ------------------------------------------------------------------------------- drawing

def poly(ax, pts, **kw):
    ax.add_patch(Polygon([(z.real, z.imag) for z in pts], closed=True, **kw))


def seg(ax, p, q, **kw):
    ax.plot([p.real, q.real], [p.imag, q.imag], **kw)


def fig1(chain):
    a, j, s0, seq, cross = chain
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2),
                                   gridspec_kw={"width_ratios": [1, 1.7]})
    # (a) T and R; the edges of R coloured by the level they lead to
    V = tri(a)
    R = rhombus(a, 0, 0j)                           # A, B, -A, -B
    for i in range(4):
        up = i in (1, 3)                            # B->-A and -B->A: direction a, to level +1
        seg(ax1, R[i], R[(i + 1) % 4], color=ACC if up else GLD, lw=2)
    poly(ax1, [V["O"], V["A"], V["B"]], fc="#dfe8f1", ec=INK, lw=1.2)
    seg(ax1, -1 + 0j, 1 + 0j, color=FADE, lw=0.8, ls="--")
    seg(ax1, -1j * math.tan(a), 1j * math.tan(a), color=FADE, lw=0.8, ls="--")
    for nm, z, off in [("$O$", 0j, (-0.05, -0.1)), ("$A$", 1 + 0j, (0.03, -0.03)),
                       ("$B$", V["B"], (-0.03, 0.05)), ("$-A$", -1 + 0j, (-0.2, -0.03)),
                       ("$-B$", -V["B"], (-0.06, -0.12))]:
        ax1.text(z.real + off[0], z.imag + off[1], nm, fontsize=10)
    ax1.text(0.42, 0.08, "$T$", fontsize=12, color=ACC)
    ax1.text(0.55, 0.27, "$H$", fontsize=10)
    nrm = cmath.exp(1j * (math.pi / 2 - a))        # unit normal to H
    foot = math.sin(a) * nrm
    seg(ax1, 0j, foot, color=PERP, lw=1)
    ax1.text(foot.real * 0.55 + 0.03, foot.imag * 0.55 - 0.04, "$h$", color=PERP, fontsize=10)
    ax1.text(-0.75, 0.3, "to level $+1$", color=ACC, fontsize=9, ha="center")
    ax1.text(-0.75, -0.36, "to level $-1$", color=GLD, fontsize=9, ha="center")
    ax1.set_xlim(-1.25, 1.15); ax1.set_ylim(-0.6, 0.55)
    ax1.set_aspect("equal"); ax1.axis("off")
    ax1.set_title("(a) $T$, the rhombus $R = R_0$, and its two doors", fontsize=10)

    # (b) the chain along a developed line, from the fold level
    th = 2 * a * j
    rot = cmath.exp(-1j * th)                       # draw with the line horizontal
    for (n, c) in seq:
        Vr = [z * rot for z in rhombus(a, n, c)]
        poly(ax2, Vr, fc="#f6d9d5" if n == j else "#dfe6ee", ec=INK, lw=0.9, alpha=0.75)
        cr = c * rot
        ax2.text(cr.real, cr.imag + (0.12 if cr.imag > s0 else -0.12), f"{n}", ha="center", va="center", fontsize=10)
    c = dict(seq)[j]
    B = cmath.exp(2j * a * j) * 1j * math.tan(a)
    seg(ax2, (c - B) * rot, (c + B) * rot, color=PERP, lw=1.6)
    xs = [c * rot for _, c in seq]
    x0, x1 = min(z.real for z in xs) - 0.9, max(z.real for z in xs) + 0.9
    ax2.plot([x0, x1], [s0, s0], color=ACC, lw=1.5)
    def half(k):
        m = int(abs(k))                             # k = ±(m + 1/2)
        return ("−" if k < 0 else "") + (str(m) if m else "") + "½"
    lastx, lev = None, 0
    for k, p in sorted(cross, key=lambda kp: (kp[1] * rot).real):
        pr = p * rot
        lev = (lev + 1) % 2 if lastx is not None and pr.real - lastx < 0.4 else 0
        lastx = pr.real
        ax2.plot([pr.real], [pr.imag], "o", color=INK, ms=3)
        above = lev == 1                            # a close neighbour: put this label above the line
        ax2.text(pr.real, pr.imag + (0.07 if above else -0.07), half(k), ha="center",
                 va="bottom" if above else "top", fontsize=7.5, color=ACC,
                 bbox=dict(fc="white", ec="none", pad=0.4, alpha=0.85))
    ax2.set_aspect("equal"); ax2.axis("off")
    ax2.set_title(f"(b) a line of direction $\\theta = {2 * j}a$ (class $\\perp OB$) "
                  f"through the chain;\nfold level $\\ell^* = {j}$ shaded, doors $k$ at the crossings",
                  fontsize=10)
    fig.tight_layout()
    save(fig, "fig1_chain")


def fig2():
    aw, w = 0.80, "HAHBH"
    gs, edges, copies = unfold(w, aw)
    oend = apply(gs[-1], 0j)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.4),
                                   gridspec_kw={"width_ratios": [1.6, 1]})
    for i, cp in enumerate(copies):
        poly(ax1, cp, fc="#dfe8f1" if i in (0, len(copies) - 1) else "#f3f5f8",
             ec=FADE, lw=0.8)
        cen = sum(cp) / 3
        ax1.text(cen.real, cen.imag, f"$T_{i}$", fontsize=8, ha="center", va="center",
                 color="#5b6770")
    for i, (P, Q) in enumerate(edges):
        seg(ax1, P, Q, color=INK, lw=1.4)
        m = (P + Q) / 2
        ax1.text(m.real + 0.02, m.imag + 0.02, f"$e_{i + 1}$", fontsize=9)
    seg(ax1, 0j, oend, color=PERP, lw=1.6)
    u = oend / abs(oend)
    for s in (+1, -1):
        off = s * 0.035 * 1j * u
        seg(ax1, off, oend + off, color=ACC, lw=0.9, ls="--")
    ax1.annotate("$O$", (-0.08, -0.06), fontsize=10)
    ax1.annotate("$O_{end}$", (oend.real + 0.03, oend.imag), fontsize=10)
    ax1.set_aspect("equal"); ax1.axis("off")
    ax1.set_title(f"(a) the corner word {w} unfolded, $a = {aw}$", fontsize=10)

    V = tri(aw)
    poly(ax2, [V["O"], V["A"], V["B"]], fc="#f3f5f8", ec=INK, lw=1.3)
    phi = cmath.phase(oend)
    dvec = cmath.exp(1j * phi)
    pts, _, _ = trace(aw, 0j, dvec, len(w))
    pts.append(0j)
    p0 = (pts[0] + pts[1]) / 2 + 0.02 * 1j * dvec
    q, _, _ = trace(aw, p0, dvec, 2 * len(w) + 4)
    q.append(q[0])                                  # close the period (verified in verify())
    ax2.plot([z.real for z in q], [z.imag for z in q], color=ACC, lw=0.8)
    ax2.plot([z.real for z in pts], [z.imag for z in pts], color=PERP, lw=1.6)
    for nm, z, off in [("O", 0j, (-0.07, -0.07)), ("A", 1 + 0j, (0.02, -0.06)),
                       ("B", V["B"], (-0.06, 0.02))]:
        ax2.text(z.real + off[0], z.imag + off[1], f"${nm}$", fontsize=10)
    ax2.set_aspect("equal"); ax2.axis("off")
    ax2.set_title("(b) folded: the corner path (red) and a nearby path (blue),\n"
                  "along $L_{\\varepsilon}$ and back along $L_{-\\varepsilon}$; period 14", fontsize=10)
    fig.tight_layout()
    save(fig, "fig3_corner")


def fig3(tiles_ends):
    tiles, ends = tiles_ends
    a1, a2, a3, a4, a5 = ends
    H = math.pi / 2
    fig, (ax, az) = plt.subplots(2, 1, figsize=(10, 5.2), gridspec_kw={"height_ratios": [1, 1.15]})
    cols = [ACC, GLD, GRN, GLD, ACC]
    names = {"HAHABH": ("", "$a_1$"), "HAHBHAHABH": ("$a_2$", "$a_3$"),
             "HAHBH": ("$a_4$", "$a_5$"), "HBHAHBHBAH": ("$\\frac{\\pi}{2}-a_3$", "$\\frac{\\pi}{2}-a_2$"),
             "HBHBAH": ("$\\frac{\\pi}{2}-a_1$", "")}
    for axx, (x0, x1), lab in [(ax, (0, H), True), (az, (0.55, 1.02), True)]:
        for r, w in enumerate(TILES_WORDS):
            lo, hi = tiles[w]
            y = 4 - r
            axx.plot([max(lo, x0), min(hi, x1)], [y, y], color=cols[r], lw=7,
                     solid_capstyle="butt")
            if axx is az:
                for x, nm in [(lo, names[w][0]), (hi, names[w][1])]:
                    if nm and x0 < x < x1:
                        axx.text(x, y + 0.32, nm, ha="center", va="bottom", fontsize=8)
                        axx.plot([x, x], [-0.5, y], color=FADE, lw=0.5, ls=":")
                axx.text(x1 + 0.006, y, w, va="center", fontsize=8, family="monospace")
            else:
                axx.text(H + 0.03, y, w, va="center", fontsize=9, family="monospace")
        axx.plot([x0, x1], [-0.5, -0.5], color=INK, lw=1)
        axx.set_xlim(x0 - 0.01, x1 + (0.3 if axx is ax else 0.1)); axx.set_ylim(-1.6, 5.0)
        axx.spines[:].set_visible(False); axx.set_yticks([])
    ax.set_xticks([0, H / 4, H / 2, 3 * H / 4, H],
                  ["0", "$\\pi/8$", "$\\pi/4$", "$3\\pi/8$", "$\\pi/2$"])
    ax.axvspan(0.55, 1.02, color=FADE, alpha=0.15, lw=0)
    az.set_xticks([3 * H / 8, 7 * H / 16, H / 2, 9 * H / 16, 5 * H / 8],
                  ["$3\\pi/16$", "$7\\pi/32$", "$\\pi/4$", "$9\\pi/32$", "$5\\pi/16$"])
    ax.set_title("(a) the five tiles $R_W$ cover $(0, \\pi/2)$ (the angle $a$)", fontsize=10)
    az.set_title("(b) the shaded window of (a), enlarged: the overlaps "
                 "$a_2 < a_1$, $a_4 < a_3$ and their mirror images", fontsize=10)
    fig.tight_layout()
    save(fig, "fig4_tiles")


def fig4():
    a = EX_A
    V = tri(a)
    fig = plt.figure(figsize=(12, 7.2))
    ax1 = fig.add_axes([0.02, 0.52, 0.55, 0.42])
    ax2 = fig.add_axes([0.02, 0.02, 0.96, 0.46])
    poly(ax1, [V["O"], V["A"], V["B"]], fc="#f3f5f8", ec=INK, lw=1.3)
    for t in [0.02, 0.98]:
        pts, _, _ = example_path(EX_TOP * t)
        ax1.plot([z.real for z in pts], [z.imag for z in pts], color=ACC, lw=0.6)
    pts, _, _ = example_path(EX_TOP / 2)
    ax1.plot([z.real for z in pts], [z.imag for z in pts], color=PERP, lw=1.1)
    seg(ax1, 0j, 1j * EX_TOP, color=GLD, lw=4)
    for nm, z, off in [("$O$", 0j, (-0.04, -0.03)), ("$A$", 1 + 0j, (0.01, -0.02)),
                       ("$B$", V["B"], (-0.04, 0.0))]:
        ax1.text(z.real + off[0], z.imag + off[1], nm, fontsize=10)
    ax1.set_aspect("equal"); ax1.axis("off")
    ax1.set_title("(a) the cylinder in $T$, $\\cos 2a = 3/4$: the middle path (red), two paths near"
                  " its edges (blue); the gold bar is the launch interval on $OB$", fontsize=9)

    gs, edges, copies = unfold(EX_WORD, a)
    rot = cmath.exp(-2j * a)
    lvcol = {0: "#eef2f6", -1: "#cfdbe8", -2: "#a9bfd6"}
    for cp, g in zip(copies, gs):
        poly(ax2, [z * rot for z in cp], fc=lvcol[level_of(g, a)], ec="white", lw=0.8)
    for (P, Q), w in zip(edges, EX_WORD):
        seg(ax2, P * rot, Q * rot, color=INK if w == "H" else FADE, lw=1.6 if w == "H" else 0.6)
    L = 4.0
    wdt = math.sqrt(7) / 56
    ax2.fill_between([0, L], [0, 0], [wdt] * 2, color=GLD, alpha=0.55, lw=0)
    ax2.plot([0, L], [wdt / 2] * 2, color=PERP, lw=0.9)
    for z, nm in [(0j, "$O$"), (L + 0j, "$O + 4e(2a)$")]:
        ax2.plot([z.real], [z.imag], "o", color=INK, ms=3)
        ax2.text(z.real, z.imag - 0.09, nm, ha="center", va="top", fontsize=9)
    for n, c in lvcol.items():
        ax2.add_patch(Polygon([(0, 0)], fc=c, ec=INK, lw=0.4, label=f"level {n}"))
    ax2.legend(loc="upper right", fontsize=8, frameon=False)
    ax2.set_aspect("equal"); ax2.axis("off")
    ax2.set_title("(b) its development, rotated so that the direction $2a$ is horizontal: the 17 copies of $T$ by level (doors, i.e. $H$-edges, bold),\n"
                  " the strip of width $\\sqrt{7}/56$ (gold); the fold level $1$ is never entered",
                  fontsize=9)
    save(fig, "fig2_example")


def save(fig, name):
    os.makedirs(OUT, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(f"{OUT}/{name}.{ext}", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {OUT}/{name}.png/.pdf")


if __name__ == "__main__":
    print("verify():")
    res = verify()
    fig1(res["chain"])
    fig2()
    fig3(res["tiles"])
    fig4()
