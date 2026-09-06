#!/usr/bin/env python3.13
"""s449_necklace.py -- PRE-REGISTERED.  The necklace's GLUING ORDER is the `+P` cycle, so
reading the chain geometrically the direction advances by `2P` per kite: a rotation by `P/Q`.

PRIOR ART: ran `rulings.py --grep` on 'necklace', 'kite', 'cyclic cover', 'branched', 'Z/Q'
-> [NCYL-287] (the necklace model, s447), [NCYL-282] (the disk), [NCYL-280] (`E = Q`, the
genus-0 base), [NCYL-283] (`Q = h + f + 2p`); 'cyclic cover'/'branched' return NOTHING in
`rulings.md` but the framing IS in the ledgers -- `veech_reframing.md` §3a / `results_ledger.md`
(Apisa's `y^{2Q} = z^P(z-1)^Q`, "the pillowcase double pulled back under `z -> z^Q`"), so
`B -> B/rho` = the Z/Q cover of the doubled triangle is PRIOR ART and is used here, not claimed.
What is NOT in the repo is the gluing ORDER (below) -- [NCYL-287] says the kites are "glued
cyclically" and that the deck rotates them, which is true of BOTH orders and pins neither.

------------------------------------------------------------------------------------------
THE MODEL, derived (no tracing anywhere in this file; it is pure arithmetic).

`B` is `Q` kites.  A kite = the rhombus of 4 triangle copies about `R`, mod rotation by `pi`
= a BIGON: two unit-length `H`-edges, corners `Z_O` (angle `2a`) and `Z_A` (`pi - 2a`), one
interior simple pole `R` (angle `pi`, ONE prong), and -- because the `+-` quotient folds the
rhombus's two diagonals onto the segments `R-Z_O` and `R-Z_A` -- an interior CROSS made of one
`L2`-copy and one `L1`-copy.  Directions are indexed mod `2Q` (angle `k*pi/(2Q)`), so
`s(L2), s(L1), s(H) = 0, Q, P`.

(1) In the kite's own frame its two `H`-edges sit at directions `+-P` and its cross at `0`
    and `Q`.  The deck `rho` rotates by `pi/Q` = `2` direction units, so kite `K_j = rho^j K_0`
    has edges at `P + 2j` and `-P + 2j`.
(2) ⇒⇒ TWO KITES ARE GLUED ONLY IF THEY SHARE AN EDGE DIRECTION, AND `-P + 2j' = P + 2j`
    FORCES `j' = j + P (mod Q)`.  So the necklace's geometric cycle is `j -> j + P`, NOT
    `j -> j + 1`: the naive reading of "glued cyclically" is CONSISTENT ONLY AT `P = 1`.
    (H1 scores exactly this, against the `+1` and `+k` orders as controls.)
(3) Index kites by geometric position `t` (`kite at position t` = `K_{tP}`).  Then kite `t`
    has edges at `a + 2Pt +- P` and the INTERFACE between positions `t` and `t+1` is the
    `H`-copy at direction `a + P(2t+1) (mod 2Q)` -- an arithmetic progression of step `2P`,
    i.e. the rotation by `P/Q` this project runs on.  `a` is fixed by the class (below).
(4) `iota` reverses the necklace.  With the fixed kite at `t = 0`, `iota: t -> -t`, so for
    odd `Q` the unique fixed INTERFACE is at `2t + 1 = Q`, i.e. `t = m := (Q-1)/2` -- pure
    combinatorics, independent of the class.  `Fix(iota)` = the cross of kite `0` (one `L2` +
    one `L1`) ⊕ that interface `h_*` (one `H`), which is [NCYL-283]'s `h + v = 3` at odd `Q`.
(5) ⇒ `Delta` = half of kite `0` (a single triangle, its two legs on the boundary) then the
    kites at `t = 1..m` in a row, with `h_*` as the far boundary -- [NCYL-287]'s chain, now
    with coordinates.

WHAT THE MODEL SAYS `p = 0` IS.  `Delta` is a disk with ONE horizontal boundary arc and TWO
vertical ones; its interior prong-ends number `m` (the interior poles, one prong each)
`+ (P-1)/2` (`Z_O`) `+ (Q-P)/2` (`Z_A`) `= Q - 1`.  `p = 0` says every one of them escapes to
a vertical arc, i.e. NO saddle connection lies in `Delta`.  ⚠ THIS FILE DOES NOT PROVE THAT.
It pins the model the statement lives in, and scores the model.

HYPOTHESES.  ⚠⚠ READ THE SCORES WITH THIS PARAGRAPH OR YOU WILL MISREAD THEM ([OPS-041],
[OPS-219]).  H1, H2 and H3 are DERIVATIONS THAT CANNOT FAIL -- their `356/356` and `712/712`
are the computer agreeing with algebra, not evidence.  This was checked, not assumed:
  * H1  `e(j) = {P+2j, -P+2j}`, and `e(j) & e(j') != {}` forces `j' = j +- P` by four one-line
        cases.  So the score is a tautology.  What H1 is FOR is the derivation itself and the
        one number in it that is not: on `316/356` centres (`Q <= 41`) the `+-P` adjacency
        DIFFERS from the naive `+-1` reading of "glued cyclically", which agree only at
        `P = 1, Q-1`.  ⚠ The first draft of H1 scored a `_closes` predicate that was true for
        EVERY step order by construction (`0/18608` on its own wrong-order control) -- the
        control fired, which is the only reason this docstring is not claiming a test.
  * H2  `offset()` CALLS `hv_sides`, so the comparison is not independent; expanding the four
        `(eps mod 2, P mod 2)` cases shows `d_star = a + Q*(P mod 2)` agrees with `hv_sides`
        identically.  ⚠ BOTH ARMS ARE STILL WORTH RECORDING -- `h_*` is VERTICAL on the `356`
        paired rows and HORIZONTAL on the `356` lone ones -- because that is the model
        REPRODUCING [NCYL-283]'s `h`/`v` split from the kite necklace rather than from side
        parities, in both arms.  A reproduction, not a measurement.
  * H3  `n_int = 2 - [(1-P/2) + (1-(Q-P)/2) + 1/2] = (Q-1)/2 = m` identically once `arcs = 3`
        and `r_corners = 1`, both of which hold for every odd `Q`.  Same status as H2.
  * H4  ⇒⇒ THE ONE GENUINE CROSS-INSTRUMENT TEST IN THIS FILE, AND IT COULD HAVE FAILED ON
        EVERY ROW: `#R-prongs into Delta == m == (Q-1)/2`.  `s439.prongs` computes inward
        prong windows from the billiard side geometry and knows nothing about kites; the
        necklace says `Delta` is a half-kite plus `m` whole ones, each with a single interior
        pole.  `712/712`.  Second leg (`total == Q - h`) is [NCYL-283]'s H0 re-run -- a
        REGRESSION check, not new evidence.
  ⚠⚠ NONE OF THIS BEARS ON `p = 0`.  Not one hypothesis here can fail in a way that tells you
    whether a saddle connection sits inside `Delta`.  This file pins COORDINATES on the
    residual; it does not move it.  Against REDUCTION-vs-RESTATEMENT it is a RESTATEMENT.

  python3.13 probes/s449_necklace.py [QMAX]
"""
from __future__ import annotations
import json
import math
import sys

sys.path.insert(0, "probes")
sys.path.insert(0, "engine")

import s439_exact_cells as s439                # noqa: E402
import s446_covering_identity as s446          # noqa: E402
import s447_overlap_disk as s447               # noqa: E402

SIDES = ("L2", "L1", "H")


def base_dir(side, P, Q):
    return {"L2": 0, "L1": Q, "H": P}[side]


# ---------------------------------------------------------------- the model

def offset(P, Q, eps):
    """`a` -- the direction of kite 0's `L2`-copy.

    Kite `t`'s cross is at directions `a + 2Pt` (`L2`) and `a + 2Pt + Q` (`L1`).  Kite 0's
    cross is `Fix(iota)`'s two non-`H` arcs, so `{a, a+Q} == {eps, eps+Q}` -- `a = eps` when
    `L2` is the horizontal arc, `a = eps + Q` when `L1` is.
    """
    hor, _ver = s446.hv_sides(P, Q, eps)
    if "L2" in hor:
        return eps % (2 * Q)
    if "L1" in hor:
        return (eps + Q) % (2 * Q)
    return None                                 # even Q -- no single fixed kite


def interfaces(P, Q, eps, order=None):
    """Directions of the `Q` interfaces, indexed by geometric position `t` (interface `t`
    separates positions `t` and `t+1`).  `order` overrides the `+P` gluing step (control)."""
    a, k = offset(P, Q, eps), P if order is None else order
    if a is None:
        return None
    return [(a + k * (2 * t + 1)) % (2 * Q) for t in range(Q)]


def kite_edges(P, Q, eps, t, order=None):
    a, k = offset(P, Q, eps), P if order is None else order
    return ((a + 2 * k * t - k) % (2 * Q), (a + 2 * k * t + k) % (2 * Q))


# ---------------------------------------------------------------- hypotheses

def deck_neighbours(P, Q, j):
    """Kites sharing an edge DIRECTION with `K_j`, from the deck alone.

    `rho` is a rotation by 2 direction units and `K_j = rho^j K_0`, so `K_j`'s edges are at
    `{P + 2j, -P + 2j}` -- that much is forced, with no choice of labelling.  Two kites can be
    glued only along a common direction, so this set IS the adjacency of the necklace.  ⚠ This
    is the whole content of H1: it is not a test of my parametrization (which is circular --
    see the s449 DO-NOT), it is a computation of who is adjacent to whom.
    """
    e = {(P + 2 * j) % (2 * Q), (-P + 2 * j) % (2 * Q)}
    return sorted(jp for jp in range(Q) if jp != j
                  and e & {(P + 2 * jp) % (2 * Q), (-P + 2 * jp) % (2 * Q)})


def h1_gluing(qmax, verbose=True):
    """Necklace adjacency is `j -> j +- P`, and it DIFFERS from the naive `j -> j +- 1`."""
    ok = bad = 0
    discriminating = agree_p1 = 0
    misses = []
    for Q in range(3, qmax + 1, 2):
        for P in range(1, Q):
            if math.gcd(P, Q) != 1:
                continue
            want = {(P) % Q, (-P) % Q}
            naive = {1 % Q, (-1) % Q}
            row_ok = all(set(deck_neighbours(P, Q, j)) == {(j + P) % Q, (j - P) % Q}
                         for j in range(Q))
            # the adjacency must also be a single Q-cycle (a necklace, not several loops)
            row_ok = row_ok and math.gcd(P, Q) == 1
            ok, bad = (ok + 1, bad) if row_ok else (ok, bad + 1)
            if not row_ok:
                misses.append((P, Q))
            if want == naive:
                agree_p1 += 1
            else:
                discriminating += 1
    if verbose:
        print(f"  H1  adjacency is `j -> j +- P`              : {ok}/{ok+bad} centres")
        print(f"      DISCRIMINATING (where `+-P` != `+-1`)   : {discriminating}/"
              f"{discriminating+agree_p1} centres -- the naive necklace is WRONG on these")
        print(f"      (`+-P` == `+-1` only at P = 1, Q-1      : {agree_p1} centres)")
        for r in misses[:6]:
            print("       miss", r)
    return ok, bad, discriminating, agree_p1


def h2_fixed_interface(qmax, verbose=True):
    """`h_*` sits at `t = m`, and the model's DIRECTION for it decides horizontal-vs-vertical
    in agreement with `hv_sides` -- in all FOUR (eps parity x P parity) cases.

    Model: `h_* = a + P*Q (mod 2Q)` with `a = eps` or `eps + Q`; `P*Q = Q*(P mod 2)`.  So the
    verdict is a two-parity computation, and `hv_sides` reaches it from side parities by a
    different route ([NCYL-283] step (2)).  Both arms occur, which is what makes it a test:
    `h_*` is VERTICAL on the paired rows and HORIZONTAL on the lone ones.
    """
    ok = bad = 0
    arms = {"vertical": 0, "horizontal": 0}
    misses = []
    for Q in range(3, qmax + 1, 2):
        for P in range(1, Q):
            if math.gcd(P, Q) != 1:
                continue
            for eps in (0, 1):
                a = offset(P, Q, eps)
                if a is None:
                    continue
                m = (Q - 1) // 2
                d_star = interfaces(P, Q, eps)[m]
                hor, ver = s446.hv_sides(P, Q, eps)
                model_vert = (d_star == (eps + Q) % (2 * Q))
                model_hor = (d_star == eps % (2 * Q))
                good = (model_vert and "H" in ver) or (model_hor and "H" in hor)
                if good:
                    ok += 1
                    arms["vertical" if model_vert else "horizontal"] += 1
                else:
                    bad += 1
                    misses.append((P, Q, eps, d_star, hor, ver))
    if verbose:
        print(f"  H2  `h_*` at t=m, its direction decides h/v : {ok}/{ok+bad} class rows")
        print(f"      BOTH ARMS POPULATED: vertical {arms['vertical']}"
              f" (paired) / horizontal {arms['horizontal']} (lone)")
        for r in misses[:6]:
            print("       miss", r)
    return ok, bad, arms, misses


def h3_budget(qmax, verbose=True):
    """Model `Fix(iota)` + interior-pole count vs `s447.geometry_of_class`."""
    ok = bad = 0
    misses = []
    for Q in range(3, qmax + 1, 2):
        for P in range(1, Q):
            if math.gcd(P, Q) != 1:
                continue
            for eps in (0, 1):
                if offset(P, Q, eps) is None:
                    continue
                hor, ver, n_int, rc = s447.geometry_of_class(P, Q, eps)
                m = (Q - 1) // 2
                # model: Fix(iota) = kite 0's cross (one L2 + one L1, one of each of h/v)
                # + h_* (one H); interior kites of Delta = m; the R-corner is R_0.
                if len(hor) + len(ver) == 3 and sorted(hor + ver) == sorted(SIDES) \
                        and round(n_int) == m and rc == 1 \
                        and 1 <= len(ver) <= 2:
                    ok += 1
                else:
                    bad += 1
                    misses.append((P, Q, eps, hor, ver, n_int, rc, m))
    if verbose:
        print(f"  H3  model budget == Gauss-Bonnet budget     : {ok}/{ok+bad} class rows")
        for r in misses[:6]:
            print("       miss", r)
    return ok, bad, misses


def h4_prongs(qmax, verbose=True):
    """The chain length IS the pole-prong count: `#R-prongs into Delta == m == (Q-1)/2`.

    The necklace says `Delta` is a half-kite plus `m` whole kites, each whole kite carrying
    exactly one interior pole with one prong.  `s439.prongs` computes the inward prong windows
    from the billiard side geometry and knows nothing about kites, so this can fail on every
    row.  Second leg: the total over `O`/`R`/`A` is `Q - h` -- that is [NCYL-283]'s H0
    (`8768/8768`), so it is a REGRESSION check here, not new evidence.
    """
    ok = bad = 0
    tot_ok = tot_bad = 0
    misses = []
    for Q in range(3, qmax + 1, 2):
        for P in range(1, Q):
            if math.gcd(P, Q) != 1:
                continue
            for eps in (0, 1):
                if offset(P, Q, eps) is None:
                    continue
                m = (Q - 1) // 2
                hor, _ver = s446.hv_sides(P, Q, eps)
                nr = len(s439.prongs("R", P, Q, eps))
                tot = s446.n_prongs(P, Q, eps)
                if nr == m:
                    ok += 1
                else:
                    bad += 1
                    misses.append((P, Q, eps, nr, m))
                tot_ok, tot_bad = ((tot_ok + 1, tot_bad) if tot == Q - len(hor)
                                   else (tot_ok, tot_bad + 1))
    if verbose:
        print(f"  H4  #R-prongs into Delta == m == (Q-1)/2    : {ok}/{ok+bad} class rows")
        print(f"      regression ([NCYL-283] H0, total == Q-h): {tot_ok}/{tot_ok+tot_bad}")
        for r in misses[:6]:
            print("       miss", r)
    return ok, bad, tot_ok, tot_bad, misses


def main():
    qmax = int(sys.argv[1]) if len(sys.argv) > 1 else 121
    print(f"s449 necklace model -- odd Q <= {qmax}, both classes, pure arithmetic\n")
    o1, b1, disc, agree = h1_gluing(qmax)
    o2, b2, arms, m2 = h2_fixed_interface(qmax)
    o3, b3, m3 = h3_budget(qmax)
    o4, b4, t_ok, t_bad, m4 = h4_prongs(qmax)
    out = {
        "qmax": qmax,
        "H1": {"ok": o1, "bad": b1, "discriminating": disc, "plusP_eq_plus1": agree},
        "H2": {"ok": o2, "bad": b2, "arms": arms, "misses": m2[:20]},
        "H3": {"ok": o3, "bad": b3, "misses": [list(map(str, r)) for r in m3[:20]]},
        "H4": {"ok": o4, "bad": b4, "regression_ok": t_ok, "regression_bad": t_bad,
               "misses": m4[:20]},
    }
    with open("data/s449_necklace.json", "w") as fh:
        json.dump(out, fh, indent=1)
    print("\n  -> data/s449_necklace.json")


if __name__ == "__main__":
    main()
