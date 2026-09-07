#!/usr/bin/env python3
"""
s506_necklace_figures.py -- the four figures of `paper_necklace.md` (item D3).

NOT a probe and NOT a claim: every object drawn here is already stated and proved in the
draft.  Two figures are SCHEMATIC and two are DATA, and each says which on its own face,
because a cartoon that reads as a measurement is its own trap (the `probes/s447_figure.py`
convention, kept by `probes/s503_paper_figures.py` and kept here).

⚠ IT READS NO `data/` DIRECTORY AND NO STORE, deliberately -- [OPS-269]'s trap (i): a
`data/` read is in nobody's dependency graph, so a script that takes one dies on a clean
clone after doing all its work.  Everything comes from `s450_chain_maps.Chain` (exact
`Z[zeta_4Q]` kite arithmetic, no tracing) and `s451_quotient_path.build` (the global
coordinate), which are the instruments the draft's own results were derived on.

PRIOR ART: `rulings.py --grep kite`, `--grep interface`, `--grep necklace figure`; listed
`probes/*figure*`.  Came back: (i) `probes/s447_figure.py` ALREADY draws a necklace
schematic (`s447_necklace.png`, [NCYL-287]) -- Figure 1(b) here is the same object and is
deliberately a redraw for the paper (labels, the `iota` axis, `D` shaded), while 1(a) (the
kite BUILT from two copies of `T`), Figure 2 (the `+P` gluing cycle), Figure 3 (`D` with
its per-kite fold data) and Figure 4 (the interval model) are drawn nowhere.  (ii) The
`+P` gluing is [NCYL-291] and the per-kite flat geometry [NCYL-292]; the global coordinate,
its closure and the quotient path are [NCYL-294] with the proof at [NCYL-303].  (iii) ⚠
[NCYL-291] DO-NOT: `+1` is WRONG on `2900/3020` centres -- Figure 2 exists to say exactly
that, and its `P = 1` panel is the agreeing case, not the general one.  (iv) ⚠ [OPS-222]/
[NCYL-292] DO-NOT (1): the `k <-> -k` pairing IS `p = 0` restated, so no figure here draws
it as evidence; `p = 0` is drawn only as the OPEN escape statement of Corollary 4.2.

  s506_fig1_kite.png     SCHEMATIC.  (a) The kite of Proposition 2.1, built: two copies of
                         `T` glued along BOTH legs -- developed as the isoceles triangle
                         `O O'' A`, whose base is then folded at `R` by `x -> 2c - x`.  The
                         fold makes `R` an interior point of cone angle `pi` (a simple
                         pole), identifies `O` with `O''` into a corner of angle `2a`, and
                         leaves the two hypotenuse copies free as the bigon's two edges.
                         (b) `B` as the cyclic chain of `Q` kites: all `O`-corners at one
                         cone point `Z_O` of angle `P*pi`, all `A`-corners at `Z_A` of
                         angle `(Q-P)*pi`, one simple pole per kite, the `iota` axis, and
                         `D = B/iota` shaded.
  s506_fig2_gluing.png   DATA.  Remark 2.3, which is the single most common way to get the
                         model wrong.  Kite `K_j` carries `H`-edges at directions
                         `{P+2j, -P+2j} (mod 2Q)`, so a shared direction forces
                         `j' = j +- P`: the gluing cycle is the star polygon `{Q/P}` on the
                         kite INDICES, and geometric position `t` holds `K_{tP}`.
                         (a) `P = 1`, `Q = 7`: the cycle is the heptagon and the naive `+1`
                         reading is right.  (b) `P = 3`, `Q = 7`: the cycle is the `{7/3}`
                         star and the naive reading is wrong at every step but the first.
                         The matched direction is printed on each interface.
  s506_fig3_chain.png    DATA.  `D = B/iota` of Proposition 2.4 at `3/11`, both classes:
                         a fan of the half-kite `0` and kites `1..m`, all sharing `Z_O` and
                         `Z_A`.  Three boundary corners (`Z_O` at `P*pi/2`, `Z_A` at
                         `(Q-P)*pi/2`, `R_0` at `pi/2`) and `m` interior poles.  Each kite
                         is marked with its FAR (larger) and NEAR interface, its sector,
                         and its transverse measures.  (a) PAIRED `eps = 0`: `ell_m = 1` is
                         the strict maximum.  (b) LONE `eps = 1`: `ell_m = 0`, the two
                         kites flanking it have an empty through band, and that degeneracy
                         is `h = 2` rather than an imposed convention.
  s506_fig4_intervals.png  DATA.  The interval model of Theorem 4.1 and Corollary 4.2 at
                         `3/11` PAIRED.  (a) The `Q` interfaces laid on one line as
                         `I_t = [u_t, u_t + ell_t]` by the alignment rule: consecutive
                         intervals share an endpoint and NEST, the through band is their
                         intersection, and the sequence is a palindrome about `t = m` and
                         the kite-`0` gap.  The offsets close around the cycle -- that is
                         Theorem 4.1, and it is what makes the drawing well posed at all.
                         (b) The `iota` quotient: the path `J_0 .. J_m` with
                         `|J_i| = |cos(i*P*pi/Q)|`, one interior pole per edge, and the
                         first run of a pole prong -- constant in `y`, so it is a
                         horizontal segment across the maximal arc of `S(y)`.

Every panel's numbers are checked before it is drawn (`verify()`), and the checks ARE the
figures' point rather than decoration: the Euler count `V - E + F = 2`, the cone angles,
the gluing law `j' = j +- P` from the direction arithmetic alone, the closed form for
`ell_t`, the palindrome, the closure `S = 0`, the nesting of consecutive intervals, and
`|J_i| = |cos(i*P*pi/Q)|`.  Three of them are scored against a FIRING control -- the naive
`+1` gluing (which fails), the `+2` gluing (which fails), and the closure under a
perturbed sector rule (which fails) -- because a check no model can fail is not a check
([OPS-041]).

Run: PYTHONPATH=.:engine:archive/scripts_2026-07:archive/scripts_2026-08:probes \\
     .venv/bin/python3.13 probes/s506_necklace_figures.py
"""
import math
import os
from math import gcd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Arc, FancyArrowPatch, Polygon, Wedge

import s450_chain_maps as cm
import s451_quotient_path as qp

OUT = "probes/figs"
INK, ACC, PERP, FADE = "#1f2933", "#2f6f9f", "#c0392b", "#b8c2cc"
GRN, GLD, VIO = "#4a8f6d", "#c08a3e", "#8f6da8"


# ------------------------------------------------------------------------ the model, bare

def cvals(P, Q, eps, step=None):
    """`c_t` of Proposition 2.2, from the direction arithmetic and nothing else.

    `step` overrides the gluing step `P` -- that is how the `+1` and `+2` controls are
    built, and they must FAIL the checks the true step passes.
    """
    s = 2 * (P if step is None else step)
    c0 = P if eps == 0 else (P + Q) % (2 * Q)
    return [(c0 + s * t) % (2 * Q) for t in range(Q)]


def hedges(P, Q, j):
    """Kite `K_j = rho^j K_0` has its two `H`-edges at these directions (mod `2Q`)."""
    return ((P + 2 * j) % (2 * Q), (-P + 2 * j) % (2 * Q))


# --------------------------------------------------------------------------- the checks

def verify(verbose=True):
    """Everything a caption below asserts, checked first.  Raises on any failure."""
    log = (lambda s: print(s)) if verbose else (lambda s: None)
    rows = [(P, Q, eps) for Q in range(5, 32, 2) for P in range(1, Q)
            if gcd(P, Q) == 1 for eps in (0, 1)]

    # (1) Proposition 2.1: the 1-skeleton and the cone angles.
    for P, Q, _e in rows:
        V, E, F = Q + 2, 3 * Q, 2 * Q
        assert V - E + F == 2, (P, Q)
        alpha = P * math.pi / (2 * Q)                       # 2a/pi = P/Q
        assert abs(Q * 2 * alpha - P * math.pi) < 1e-12
        assert abs(Q * (math.pi - 2 * alpha) - (Q - P) * math.pi) < 1e-9
    log(f"  (1) V-E+F=2 and the two cone angles P*pi/(Q-P)*pi   {len(rows)}/{len(rows)}")

    # (2) Proposition 2.2: a shared H-direction forces j' = j +- P, and the interface
    #     between geometric positions t and t+1 is the H-copy at a + P(2t+1).
    ok = 0
    for P, Q, eps in rows:
        for j in range(Q):
            hi_j, lo_j = hedges(P, Q, j)
            partners = {jj for jj in range(Q) if set(hedges(P, Q, jj)) & {hi_j, lo_j}}
            assert partners == {j, (j + P) % Q, (j - P) % Q}, (P, Q, j, partners)
        a = eps if (P % 2 == eps % 2) else eps            # offset fixed by the class
        for t in range(Q):
            shared = set(hedges(P, Q, (t * P) % Q)) & set(hedges(P, Q, ((t + 1) * P) % Q))
            assert len(shared) == 1, (P, Q, t)
            assert shared.pop() == (P * (2 * t + 1)) % (2 * Q), (P, Q, t)
        ok += 1
    log(f"  (2) j' = j +- P, and interface(t,t+1) = P(2t+1) mod 2Q   {ok}/{len(rows)}")

    # (2b) CONTROL -- the naive +1 reading and the +2 reading must FAIL (2)'s second half.
    bad1 = bad2 = tot = 0
    for P, Q, eps in rows:
        tot += 1
        for step, hit in ((1, "1"), (2, "2")):
            good = all(set(hedges(P, Q, (t * step) % Q))
                       & set(hedges(P, Q, ((t + 1) * step) % Q)) for t in range(Q))
            if not good:
                if step == 1:
                    bad1 += 1
                else:
                    bad2 += 1
    assert bad1 > 0 and bad2 > 0
    log(f"  (2b) CONTROL: +1 gluing fails on {bad1}/{tot} rows, +2 on {bad2}/{tot} "
        f"-- the checks can fail")

    # (3) The closed form for ell_t -- and ⚠ IT CARRIES THE CLASS OFFSET `a`, which the
    #     draft's abstract, §1.2 and Prop. 2.2 all dropped (s506).  The interface between
    #     positions t and t+1 is the H-copy at direction `a + P(2t+1) mod 2Q` with
    #     `a = 0` (eps=0) / `Q` (eps=1), so its transverse measure is
    #     `|sin((a + P(2t+1))pi/2Q)|` -- the BARE `|sin(P(2t+1)pi/2Q)|` is the eps=0 case
    #     only, and at eps=1 it is `|cos(P(2t+1)pi/2Q)|` instead.  Both are `sin(pi c_t/2Q)`
    #     with Prop. 2.2's own `c_0`, which is why the error is invisible downstream.
    bare_fail = 0
    for P, Q, eps in rows:
        ch = cm.Chain(P, Q, eps)
        c = cvals(P, Q, eps)
        assert ch.c == c, (P, Q, eps)
        a = 0 if eps == 0 else Q
        for t in range(Q):
            want = abs(math.sin((a + P * (2 * t + 1)) * math.pi / (2 * Q)))
            alt = (abs(math.sin(P * (2 * t + 1) * math.pi / (2 * Q))) if eps == 0 else
                   abs(math.cos(P * (2 * t + 1) * math.pi / (2 * Q))))
            assert abs(ch.ell[t].f - want) < 1e-11, (P, Q, eps, t)
            assert abs(ch.ell[t].f - alt) < 1e-11, (P, Q, eps, t)
            assert abs(ch.ell[t].f - math.sin(math.pi * c[t] / (2 * Q))) < 1e-11
            if abs(ch.ell[t].f - abs(math.sin(P * (2 * t + 1)
                                              * math.pi / (2 * Q)))) >= 1e-11:
                bare_fail += 1
    log(f"  (3) ell_t = |sin((a + P(2t+1))pi/2Q)| = sin(pi c_t/2Q), a = 0 / Q   "
        f"{len(rows)}/{len(rows)}")
    assert bare_fail > 0
    log(f"  (3b) CONTROL: the OFFSETLESS form |sin(P(2t+1)pi/2Q)| fails on {bare_fail} "
        f"interfaces — every eps=1 interface in range; the offset is not cosmetic")

    # (4) Corollary 4.2's palindrome, and the tie: kite 0 and no other.
    for P, Q, eps in rows:
        ch = cm.Chain(P, Q, eps)
        for t in range(Q):
            assert ch.ell[(-1 - t) % Q].eq(ch.ell[t]), (P, Q, eps, t)
        assert [k for k in range(Q) if ch.tie[k]] == [0], (P, Q, eps)
    log(f"  (4) palindrome ell_{{-1-t}} = ell_t, and the tie is kite 0 alone   "
        f"{len(rows)}/{len(rows)}")

    # (5) Theorem 4.1: the global coordinate CLOSES -- and a perturbed sector rule does not.
    small = [r for r in rows if r[1] <= 19]
    for P, Q, eps in small:
        _ch, lo, hi, closes = qp.build(P, Q, eps)
        assert closes, (P, Q, eps)
        for t in range(Q):                              # I_t = [u_t, u_t + ell_t]
            assert hi[t].f >= lo[t].f - 1e-12
    log(f"  (5) the global coordinate closes   {len(small)}/{len(small)}")

    fails = 0
    for P, Q, eps in small:
        ch = cm.Chain(P, Q, eps, sector="shift")        # a FIRING control ([OPS-219])
        u = {0: ch.zero}
        for k in range(1, Q):
            far, near, d = ch.far[k], ch.near[k], ch.delta[k]
            if far in u:
                u[near] = u[far] + d if ch.osec[k] else u[far]
            else:
                u[far] = u[near] - d if ch.osec[k] else u[near]
        far, near, d = ch.far[0], ch.near[0], ch.delta[0]
        if not (u[far] + d if ch.osec[0] else u[far]).eq(u[near]):
            fails += 1
    assert fails > 0
    log(f"  (5b) CONTROL: the shifted sector rule fails to close on {fails}/{len(small)} "
        f"rows")

    # (6) The nesting Figure 4(a) draws: consecutive intervals share an endpoint, and one
    #     contains the other.  Definitional, not a finding ([NCYL-292]'s s469 marker) --
    #     asserted because the DRAWING is only legible if it is true.
    for P, Q, eps in small:
        _ch, lo, hi, _c = qp.build(P, Q, eps)
        for t in range(Q):
            u = (t + 1) % Q
            shares = abs(lo[t].f - lo[u].f) < 1e-11 or abs(hi[t].f - hi[u].f) < 1e-11
            nests = ((lo[t].f <= lo[u].f + 1e-11 and hi[u].f <= hi[t].f + 1e-11)
                     or (lo[u].f <= lo[t].f + 1e-11 and hi[t].f <= hi[u].f + 1e-11))
            assert shares and nests, (P, Q, eps, t)
    log(f"  (6) consecutive interfaces share an endpoint and nest   "
        f"{len(small)}/{len(small)}")

    # (7) Corollary 4.2's closed form for the quotient path.
    for P, Q, eps in small:
        ch = cm.Chain(P, Q, eps)
        m = ch.m_idx
        for i in range(m + 1):
            want = (abs(math.sin(i * P * math.pi / Q)) if ch.lone
                    else abs(math.cos(i * P * math.pi / Q)))
            assert abs(ch.ell[(m - i) % Q].f - want) < 1e-11, (P, Q, eps, i)
            assert ch.ell[(m + i) % Q].eq(ch.ell[(m - i) % Q])
    log(f"  (7) |J_i| = |cos(iPpi/Q)| PAIRED / |sin(iPpi/Q)| LONE   "
        f"{len(small)}/{len(small)}")

    # (8) The two degeneracies of Proposition 2.4, both FORCED: the tie is kite 0, and the
    #     empty through band occurs exactly at the two kites flanking a zero interface,
    #     which happens exactly in the LONE class.
    zc = {True: set(), False: set()}
    for P, Q, eps in rows:
        ch = cm.Chain(P, Q, eps)
        zero = [t for t in range(Q) if ch.ell[t].is_zero()]
        zc[ch.lone].add(len(zero))
        if zero:
            assert zero == [ch.m_idx], (P, Q, eps, zero)
    assert zc[True] == {1} and zc[False] == {0}, zc
    log(f"  (8) a zero interface occurs in the LONE class, at t = m, and only there   "
        f"{len(rows)}/{len(rows)}")
    return True


# ---------------------------------------------------------------- F1: the kite, and B

def fig1_kite(P=3, Q=7, path=None):
    path = path or f"{OUT}/s506_fig1_kite.png"
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(13.4, 6.2))
    alpha = P * math.pi / (2 * Q)
    c = 1.0 / math.tan(alpha)

    # ---- (a) the kite, built ----
    O, R, A, O2 = (0.0, 0.0), (c, 0.0), (c, 1.0), (2 * c, 0.0)
    axL.add_patch(Polygon([O, R, A], closed=True, fc="#eef4f9", ec="none", zorder=1))
    axL.add_patch(Polygon([R, O2, A], closed=True, fc="#f6f1e8", ec="none", zorder=1))
    axL.add_patch(Polygon([O, O2, A], closed=True, fc="none", ec=INK, lw=2.0, zorder=4))
    axL.plot([c, c], [0, 1], "--", lw=1.4, color=FADE, zorder=3)
    axL.annotate("the two copies of $T$,\nglued along $L_1$ (dashed)", (0.42 * c, 0.92),
                 fontsize=9, color=FADE, ha="center", va="bottom")
    for x0, x1, col in ((0, c, GRN), (c, 2 * c, GRN)):
        axL.plot([x0, x1], [0, 0], "-", lw=4.0, color=col, solid_capstyle="butt",
                 zorder=5)
    axL.add_patch(FancyArrowPatch((0.30 * c, -0.10), (1.70 * c, -0.10),
                                  arrowstyle="<->", mutation_scale=13, lw=1.5,
                                  color=GRN, zorder=6))
    axL.annotate(r"the base is FOLDED at $R$:  $x \mapsto 2c-x$" "\n"
                 r"(this is the gluing of the two $L_2$ legs)", (c, -0.155),
                 fontsize=9.5, color=GRN, ha="center", va="top")
    for nm, pt, dx, dy, ha in ((r"$O$", O, -12, -3, "right"), (r"$O''$", O2, 12, -3,
                                                               "left"),
                               (r"$A$", A, -13, 4, "right"), (r"$R$", R, 15, -13,
                                                              "left")):
        axL.plot([pt[0]], [pt[1]], "o", ms=7, color=INK, zorder=8)
        axL.annotate(nm, pt, textcoords="offset points", xytext=(dx, dy), fontsize=13,
                     color=INK, ha=ha, va="center")
    axL.plot([c], [0], "o", ms=11, mfc="white", mec=PERP, mew=2.2, zorder=9)
    for pt, th1, th2, lab, lx, ly in ((O, 0, math.degrees(alpha), r"$\alpha$", .32, .062),
                                      (O2, 180 - math.degrees(alpha), 180, r"$\alpha$",
                                       2 * c - .32, .062)):
        axL.add_patch(Arc(pt, .50, .50, theta1=th1, theta2=th2, lw=1.3, color=ACC,
                          zorder=6))
        axL.annotate(lab, (lx, ly), fontsize=10, color=ACC, ha="center")
    axL.annotate(r"$O \sim O''$ — one corner of angle $2\alpha$" "\n"
                 r"$A$ — one corner of angle $2(\pi/2-\alpha) = \pi - 2\alpha$",
                 (c, 1.94), fontsize=10, color=ACC, ha="center", va="top")
    axL.annotate(r"$R$ — the two right angles meet: cone angle $\pi$," "\n"
                 r"an interior SIMPLE POLE", (c, 1.56), fontsize=10, color=PERP,
                 ha="center", va="top")
    axL.annotate(r"$H$", (0.5 * c - .06, 0.5 + .04), fontsize=12, color=GLD, ha="center")
    axL.annotate(r"$H$", (1.5 * c + .06, 0.5 + .04), fontsize=12, color=GLD, ha="center")
    for a, b in ((O, A), (O2, A)):
        axL.plot([a[0], b[0]], [a[1], b[1]], "-", lw=3.2, color=GLD, zorder=5)
    axL.annotate("the two free edges are the copies of $H$ —\n"
                 "the bigon's two sides, and what the necklace glues along",
                 (c, -0.40), fontsize=9.5, color=GLD, ha="center", va="top")
    axL.set_title(rf"(a) SCHEMATIC — the KITE (Prop. 2.1), drawn at $2\alpha/\pi = {P}/{Q}$."
                  "\n"
                  r"Two copies of $T$ glued along BOTH legs = a bigon with one interior pole.",
                  fontsize=10.5, color=INK)
    axL.set_xlim(-0.42, 2 * c + 0.42)
    axL.set_ylim(-0.92, 2.06)
    axL.set_aspect("equal")
    axL.axis("off")

    # ---- (b) the necklace ----
    ri, ro = 0.50, 1.18
    axR.add_patch(plt.Circle((0, 0), ro, fc="#f7f9fb", ec=INK, lw=1.8, zorder=1))
    axR.add_patch(plt.Circle((0, 0), ri, fc="white", ec=INK, lw=1.8, zorder=2))
    # `iota` is the reflection in the vertical axis: it fixes the H-edge at 90 degrees
    # and, Q being odd, the kite straddling 270.  D is exactly one side of that axis.
    axR.add_patch(Wedge((0, 0), ro, 90, 270, width=ro - ri, fc=ACC, alpha=.13,
                        ec="none", zorder=1))
    for k in range(Q):
        a0 = 2 * math.pi * k / Q + math.pi / 2
        axR.plot([ri * math.cos(a0), ro * math.cos(a0)],
                 [ri * math.sin(a0), ro * math.sin(a0)], "-", lw=2.0, color=GLD, zorder=4)
        am = a0 - math.pi / Q
        rp = 0.5 * (ri + ro)
        px, py = rp * math.cos(am), rp * math.sin(am)
        axR.plot([ri * math.cos(am), px], [ri * math.sin(am), py], "-", lw=1.3,
                 color=GRN, zorder=3)
        axR.plot([px, ro * math.cos(am)], [py, ro * math.sin(am)], "-", lw=1.3,
                 color=ACC, zorder=3)
        axR.plot([px], [py], "o", ms=7, mfc="white", mec=PERP, mew=1.9, zorder=6)
    axR.plot([0, 0], [-1.34, 1.34], "--", lw=1.6, color=PERP, zorder=7)
    axR.annotate(r"axis of $\iota$", (0.05, 1.40), fontsize=10, color=PERP, ha="left")
    axR.annotate(r"$Z_O$", (0, 0), fontsize=13, color=INK, ha="center", va="center")
    axR.annotate(r"$Z_A$", (0, ro + 0.09), fontsize=13, color=INK, ha="center")
    axR.annotate(r"$D = B/\iota$ (shaded):" "\n"
                 r"the half-kite on the axis" "\n" r"and $m = (Q\!-\!1)/2$ kites",
                 (-1.42, -1.06), fontsize=10, color=ACC, ha="center")
    for lab, col, y in ((r"$H$-edge $Z_O\!-\!Z_A$ (the gluings)", GLD, 1.30),
                        (r"$L_2$-edge $Z_O\!-\!R_i$", GRN, 1.16),
                        (r"$L_1$-edge $R_i\!-\!Z_A$", ACC, 1.02),
                        (r"$R_i$ — one simple pole per kite", PERP, 0.88)):
        axR.annotate(lab, (1.44, y), fontsize=9, color=col, ha="center")
    axR.annotate(r"the inner circle collapses to $Z_O$: $Q$ corners of angle "
                 r"$2\alpha$, total $P\pi$" "\n"
                 r"the outer circle collapses to $Z_A$: $Q$ corners of angle "
                 r"$\pi-2\alpha$, total $(Q\!-\!P)\pi$" "\n"
                 rf"$V = Q+2$,  $E = 3Q$,  $F = 2Q$,  $\chi = 2$   (drawn: $Q = {Q}$)",
                 (0, -1.46), fontsize=9.5, color=INK, ha="center", va="center")
    axR.set_title(r"(b) SCHEMATIC — $B$ is a cyclic chain of $Q$ kites (Prop. 2.1)."
                  "\n"
                  r"$\iota$ fixes one kite and the opposite $H$-edge, so $D$ is a chain.",
                  fontsize=10.5, color=INK)
    axR.set_xlim(-1.72, 2.34)
    axR.set_ylim(-1.78, 1.52)
    axR.set_aspect("equal")
    axR.axis("off")

    fig.tight_layout()
    fig.savefig(path, dpi=200)
    fig.savefig(path.replace(".png", ".pdf"))
    print("wrote", path)


# ------------------------------------------------------------------ F2: the +P gluing

def _gluing_panel(ax, P, Q, tag):
    """The gluing cycle on the kite INDICES, drawn as the star polygon `{Q/P}`."""
    pos = {j: (math.cos(math.pi / 2 - 2 * math.pi * j / Q),
               math.sin(math.pi / 2 - 2 * math.pi * j / Q)) for j in range(Q)}
    ax.add_patch(plt.Circle((0, 0), 1.0, fc="none", ec=FADE, lw=1.0, ls=":", zorder=1))
    for t in range(Q):
        j, j2 = (t * P) % Q, ((t + 1) * P) % Q
        (x1, y1), (x2, y2) = pos[j], pos[j2]
        ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>",
                                     mutation_scale=12, lw=1.8,
                                     color=PERP if t == 0 else ACC,
                                     alpha=1.0 if t == 0 else .70,
                                     shrinkA=13, shrinkB=13, zorder=3))
        # the shared H-direction, placed a THIRD of the way along the arrow so the labels
        # stay apart when the cycle is a star and every chord runs near the centre
        d = (P * (2 * t + 1)) % (2 * Q)
        mx, my = x1 + 0.33 * (x2 - x1), y1 + 0.33 * (y2 - y1)
        ax.annotate(str(d), (mx, my), fontsize=8, color=INK, ha="center", va="center",
                    bbox=dict(fc="white", ec=FADE, lw=.5, pad=0.9, alpha=.92), zorder=4)
    for j in range(Q):
        x, y = pos[j]
        t = next(tt for tt in range(Q) if (tt * P) % Q == j)
        ax.plot([x], [y], "o", ms=21, mfc="white", mec=INK, mew=1.5, zorder=5)
        ax.annotate(rf"$K_{{{j}}}$", (x, y), fontsize=9.5, color=INK, ha="center",
                    va="center", zorder=6)
        ax.annotate(rf"$t\!=\!{t}$", (1.30 * x, 1.30 * y), fontsize=8.5, color=GRN,
                    ha="center", va="center", zorder=6)
    ax.set_title(tag, fontsize=10.5, color=INK)
    ax.set_xlim(-1.62, 1.62)
    ax.set_ylim(-1.62, 1.72)
    ax.set_aspect("equal")
    ax.axis("off")


def fig2_gluing(Q=7, path=None):
    path = path or f"{OUT}/s506_fig2_gluing.png"
    fig, axes = plt.subplots(1, 2, figsize=(12.6, 6.4))
    for ax, P in zip(axes, (1, 3)):
        naive = all(set(hedges(P, Q, t)) & set(hedges(P, Q, (t + 1) % Q))
                    for t in range(Q))
        _gluing_panel(
            ax, P, Q,
            rf"({'a' if P == 1 else 'b'}) $P = {P}$, $Q = {Q}$ — the gluing cycle is "
            rf"$\{{{Q}/{P}\}}$." "\n"
            + (r"Here $j' = j+1$ too, so the naive reading happens to agree."
               if naive else
               r"$j' = j+1$ shares NO direction: the naive reading is wrong."))
    fig.text(0.5, 0.045,
             r"green $t$ = geometric position, and position $t$ holds $K_{tP}$   ·   "
             r"boxed number = the shared $H$-direction $P(2t+1)$ mod $2Q$   ·   "
             r"red arrow = the step out of $K_0$",
             fontsize=9.5, color=INK, ha="center")
    fig.suptitle(r"DATA — Remark 2.3: kite $K_j$ has $H$-edges at $\{P+2j,\,-P+2j\}$, so a "
                 r"shared direction forces $j' = j \pm P$." "\n"
                 r"The naive $+1$ reading agrees with the truth only at $P = 1$ and "
                 r"$P = Q-1$.", fontsize=11, color=INK)
    fig.tight_layout(rect=(0, 0.055, 1, .90))
    fig.savefig(path, dpi=200)
    fig.savefig(path.replace(".png", ".pdf"))
    print("wrote", path)


# ------------------------------------------------------------------------- F3: the chain

def _chain_panel(ax, P, Q, eps):
    """`D` drawn in the GLOBAL coordinate: position `t` occupies `[t, t+1]` and is bounded
    by interfaces `t-1` and `t`, each a vertical segment `[u, u + ell]`.  That is the frame
    in which Proposition 2.4's content is visible rather than asserted -- the through band
    is the overlap of the two bounding segments and the fold is the protruding tail."""
    ch, lo, hi, closes = qp.build(P, Q, eps)
    assert closes
    m = ch.m_idx
    for t in range(1, m + 1):                            # kite at position t
        a, b = (t - 1) % Q, t
        far, near = ch.far[t], ch.near[t]
        band = (max(lo[a].f, lo[b].f), min(hi[a].f, hi[b].f))
        ax.add_patch(Polygon([(t, lo[a].f), (t, hi[a].f), (t + 1, hi[b].f),
                              (t + 1, lo[b].f)], closed=True, fc="#eef4f9", ec=FADE,
                             lw=.8, zorder=1))
        if band[1] > band[0]:                            # the THROUGH band
            ax.add_patch(Polygon([(t, band[0]), (t, band[1]), (t + 1, band[1]),
                                  (t + 1, band[0])], closed=True, fc=GRN, alpha=.22,
                                 ec="none", zorder=2))
        fx = t if far == a else t + 1                    # the FOLD, on the far interface
        flo, fhi = ch.flo[t].f + lo[far].f, ch.fhi[t].f + lo[far].f
        ax.plot([fx, fx], [flo, fhi], "-", lw=6.5, color=GLD, solid_capstyle="butt",
                zorder=4)
        ax.plot([fx], [ch.ffix[t].f + lo[far].f], "o", ms=8, mfc="white", mec=PERP,
                mew=2.0, zorder=6)
        ax.annotate(rf"$R_{{{t}}}$", (fx, ch.ffix[t].f + lo[far].f),
                    textcoords="offset points",
                    xytext=(11 if fx == t else -11, 9 if t % 2 else -9),
                    fontsize=9, color=PERP, ha="left" if fx == t else "right",
                    va="center", zorder=7)
        ax.annotate(rf"$K_{{{(t * P) % Q}}}$" "\n"
                    rf"{'$O$' if ch.osec[t] else '$A$'}-sec",
                    (t + .5, max(hi[a].f, hi[b].f)), textcoords="offset points",
                    xytext=(0, 5), fontsize=8.2, color=INK, ha="center", va="bottom")
        if band[1] <= band[0] + 1e-11:
            ax.annotate("EMPTY\nthrough band", (t + .5, (band[0] + band[1]) / 2),
                        fontsize=8, color=PERP, ha="center", va="center", zorder=8,
                        bbox=dict(fc="white", ec="none", pad=1.0, alpha=.85))
    for t in range(m + 1):                               # the interfaces themselves
        x = t + 1
        col = PERP if t == m else ACC
        ax.plot([x, x], [lo[t].f, hi[t].f], "-", lw=2.4, color=col, zorder=5)
        ax.annotate(rf"$\ell_{{{t}}}\!=\!{ch.ell[t].f:.3f}$", (x, lo[t].f),
                    textcoords="offset points", xytext=(0, -9), fontsize=8,
                    color=col, ha="center", va="top", rotation=90)
    # the left boundary: the kite-0 axis, a Fix arc, carrying the corner R_0.  Kite 0 is
    # the TIE (its two interfaces are equal), so the half of it inside D is a rectangle.
    ax.add_patch(Polygon([(0.5, lo[0].f), (0.5, hi[0].f), (1, hi[0].f), (1, lo[0].f)],
                         closed=True, fc="#eef4f9", ec=FADE, lw=.8, zorder=1))
    ax.plot([0.5, 0.5], [lo[0].f - .04, hi[0].f + .04], "-", lw=2.6, color=PERP,
            zorder=5)
    ax.plot([0.5], [(lo[0].f + hi[0].f) / 2], "s", ms=9, mfc="white", mec=PERP, mew=2.0,
            zorder=7)
    ax.annotate(r"$R_0$" "\n" r"($\pi/2$, on $\partial D$)", (0.5, hi[0].f),
                textcoords="offset points", xytext=(0, 9), fontsize=8.5, color=PERP,
                ha="center", va="bottom")
    ax.annotate("half-\nkite 0", (0.75, (lo[0].f + hi[0].f) / 2), fontsize=8,
                color=FADE, ha="center", va="center")
    cls = "LONE" if ch.lone else "PAIRED"
    note = (rf"$\ell_m = 1$ is the strict MAXIMUM: interface $m\!-\!1$ is kite $m$'s NEAR"
            "\n" r"edge, so every leaf entering kite $m$ passes through."
            if not ch.lone else
            rf"$\ell_m = 0$ — the fixed interface DEGENERATES. Its two flanking kites"
            "\n" r"have an empty through band, and that is $h = 2$: forced, not imposed.")
    ax.annotate(note, (0.5, -0.30), xycoords=("axes fraction", "axes fraction"),
                fontsize=9, color=PERP if ch.lone else ACC, ha="center", va="center")
    ax.set_title(rf"({'a' if not ch.lone else 'b'}) ${P}/{Q}$, $\varepsilon = {eps}$ — "
                 rf"{cls}.   $m = (Q\!-\!1)/2 = {m}$ interior poles $R_1,\dots,R_m$.",
                 fontsize=10.5, color=INK)
    ax.set_xlim(0.1, m + 1.75)
    ax.set_ylim(-0.62, 1.62)
    ax.axis("off")


def fig3_chain(P=3, Q=11, path=None):
    path = path or f"{OUT}/s506_fig3_chain.png"
    fig, axes = plt.subplots(2, 1, figsize=(11.6, 10.4))
    paired = 0 if (0 != P % 2) else 1                   # `lone` iff eps == P % 2
    for ax, eps in zip(axes, (paired, 1 - paired)):
        _chain_panel(ax, P, Q, eps)
    fig.text(0.5, 0.030,
             "Horizontal = the chain; vertical = the transverse coordinate of Thm. 4.1.  "
             "Every kite touches BOTH cone points, so the strip is a drawing convention:\n"
             r"the whole lower boundary is $Z_O$ (angle $P\pi/2$ in $D$), the whole upper "
             r"boundary is $Z_A$ (angle $(Q-P)\pi/2$)." "\n"
             "Green = the through band, gold = the fold, and the pole is the fold's "
             "FIXED POINT.", fontsize=9.5, color=INK, ha="center")
    fig.suptitle(r"DATA — $D = B/\iota$ (Prop. 2.4): the half-kite $0$ followed by "
                 r"$m$ kites, bounded by the two arcs of $\mathrm{Fix}(\iota)$." "\n"
                 r"On each kite the fold lives on the FAR (larger) interface and the "
                 r"NEAR one is entirely a through band.", fontsize=11, color=INK)
    fig.tight_layout(rect=(0, 0.075, 1, .93))
    fig.savefig(path, dpi=200)
    fig.savefig(path.replace(".png", ".pdf"))
    print("wrote", path)


# --------------------------------------------------------------------- F4: the intervals

def first_run(ch, lo, hi, k):
    """The maximal arc of `S(y) = {t : y in I_t}` through the prong of `R_k`, with `y`
    the prong's (constant) global coordinate.  No dynamics: the paper's point is that the
    leaf is constant on this arc, so the arc IS the run.

    ⚠ The band test is `Chain.cmp`, NOT a float comparison, because that is what
    `s451_quotient_path.first_run_escapes` uses -- a drawn arc that disagreed with the
    instrument's escape verdict would be a picture of a different object.
    """
    y, d, t0 = qp.pole_y(ch, lo, hi, k)
    run, t = [t0], t0
    while len(run) <= ch.Q:
        t2 = (t + d) % ch.Q
        if not (ch.cmp(y, lo[t2]) > 0 and ch.cmp(y, hi[t2]) < 0):
            break
        run.append(t2)
        t = t2
    return y.f, run, d


def fig4_intervals(P=3, Q=11, path=None):
    path = path or f"{OUT}/s506_fig4_intervals.png"
    eps = 0 if (0 != P % 2) else 1                      # the PAIRED class
    ch, lo, hi, closes = qp.build(P, Q, eps)
    assert closes
    m = ch.m_idx
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(13.8, 6.6))

    # ---- (a) the Q interfaces on one line ----
    for t in range(Q):
        y = Q - 1 - t
        col = PERP if t == m else (VIO if t == 0 else ACC)
        axL.plot([lo[t].f, hi[t].f], [y, y], "-", lw=6.0, color=col,
                 solid_capstyle="butt", zorder=3)
        for x in (lo[t].f, hi[t].f):
            axL.plot([x], [y], "|", ms=13, color=INK, zorder=4)
        axL.annotate(rf"$I_{{{t}}}$", (lo[t].f - 0.035, y), fontsize=9, color=col,
                     ha="right", va="center")
        axL.annotate(rf"${ch.ell[t].f:.3f}$", (hi[t].f + 0.035, y), fontsize=8.5,
                     color=FADE, ha="left", va="center")
        if t + 1 < Q:                                   # the shared endpoint
            u = t + 1
            xs = lo[t].f if abs(lo[t].f - lo[u].f) < 1e-11 else hi[t].f
            axL.plot([xs, xs], [y, y - 1], ":", lw=1.2, color=GRN, zorder=2)
    axL.annotate("green dots: the shared endpoint of consecutive interfaces\n"
                 "(their intersection is the THROUGH band; the protruding\n"
                 "tail of the larger is where the FOLD lives)",
                 (0.02, -1.55), fontsize=9, color=GRN, ha="left", va="center")
    axL.annotate(rf"the sequence is a PALINDROME about $t = m = {m}$"
                 "\n" r"and the kite-$0$ gap — that is $\iota$",
                 (0.02, -2.75), fontsize=9, color=PERP, ha="left", va="center")
    axL.set_title(r"(a) DATA — the $Q$ interfaces in the global coordinate," "\n"
                  r"$I_t = [u_t,\, u_t + \ell_t]$.  The offsets CLOSE around the cycle "
                  r"(Thm. 4.1)." "\n"
                  r"Without that this picture would not exist.", fontsize=10.5,
                  color=INK)
    axL.set_ylim(-3.6, Q + 0.4)
    axL.axis("off")

    # ---- (b) the quotient path, in the SAME global coordinate, with every prong's run --
    qidx = lambda t: min((m - t) % Q, (t - m) % Q)       # noqa: E731  position -> node
    ylo = min(lo[(m - i) % Q].f for i in range(m + 1))
    yhi = max(hi[(m - i) % Q].f for i in range(m + 1))
    pad = 0.10 * (yhi - ylo)
    for i in range(m + 1):
        t = (m - i) % Q
        axR.plot([i, i], [lo[t].f, hi[t].f], "-", lw=6.0, color=ACC,
                 solid_capstyle="butt", zorder=3)
        axR.annotate(rf"$J_{{{i}}}$", (i, yhi + 0.35 * pad), fontsize=10, color=ACC,
                     ha="center", va="bottom")
        axR.annotate(rf"${ch.ell[t].f:.3f}$", (i, ylo - 0.35 * pad), fontsize=8.5,
                     color=FADE, ha="center", va="top")
        assert abs(ch.ell[t].f - abs(math.cos(i * P * math.pi / Q))) < 1e-11
    esc = 0
    for i in range(m):                     # the pole on edge J_i--J_{i+1} is kite m-i's,
        k = (m - i) % Q                    # and it is INTERIOR to that edge: x = i + 0.5
        y, _d, _far = qp.pole_y(ch, lo, hi, k)
        run_y, run, _dd = first_run(ch, lo, hi, k)
        away = qp.first_run_escapes(ch, lo, hi, k)
        esc += bool(away)
        col = GRN if away else GLD
        qs = [qidx(tt) for tt in run] + [i + 0.5]
        axR.plot([min(qs), max(qs)], [run_y, run_y], "-", lw=1.8, color=col, alpha=.9,
                 zorder=6)
        axR.plot([i + 0.5], [y.f], "o", ms=7, mfc="white", mec=PERP, mew=1.8, zorder=7)
    axR.annotate(rf"one interior pole per EDGE of the path ({m} of them), each at its "
                 rf"fold's FIXED POINT" "\n"
                 r"horizontal bars = the first run, at CONSTANT $y$:   green = reaches "
                 rf"an end of the path ({esc}/{m} here),   gold = folds first",
                 (m / 2, ylo - 1.5 * pad), fontsize=9, color=INK, ha="center", va="top")
    axR.annotate(r"$|J_i| = |\cos(iP\pi/Q)|$   (PAIRED),   "
                 r"$|\sin(iP\pi/Q)|$   (LONE)", (m / 2, yhi + 2.6 * pad), fontsize=10.5,
                 color=ACC, ha="center", va="bottom")
    axR.annotate(r"$p = 0$  $\Longleftrightarrow$  every interior pole of the path "
                 r"escapes to an END of it" "\n"
                 r"— a statement about $(Q+1)/2$ nested intervals, and OPEN."
                 "\n"
                 r"There is no static criterion: escaping on the FIRST run is common and "
                 r"not universal (Rem. 4.3).",
                 (m / 2, ylo - 3.3 * pad), fontsize=9.5, color=INK, ha="center",
                 va="top")
    axR.set_title(r"(b) DATA — the $\iota$ quotient (Cor. 4.2): the path "
                  r"$J_0,\dots,J_m$, $J_i := I_{m-i} = I_{m+i}$," "\n"
                  rf"drawn in the same coordinate as (a).  ${P}/{Q}$, PAIRED.",
                  fontsize=10.5, color=INK)
    axR.set_xlim(-0.7, m + 0.7)
    axR.set_ylim(ylo - 7.5 * pad, yhi + 4.4 * pad)
    axR.axis("off")

    fig.tight_layout()
    fig.savefig(path, dpi=200)
    fig.savefig(path.replace(".png", ".pdf"))
    print("wrote", path)


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    print("verify():")
    verify()
    print("all checks pass — drawing\n")
    fig1_kite()
    fig2_gluing()
    fig3_chain()
    fig4_intervals()
