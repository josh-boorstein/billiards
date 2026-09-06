#!/usr/bin/env python3
"""
s462_zz_relation.py -- PRE-REGISTERED.  Queue item (II'): [NECK-GB] = (G0b), the swapped
`Z`-`Z` edge.  What the ONE known witness actually consumes, and what that costs a proof.

WHAT THIS IS.  NOT a proof of (G0b).  An ANATOMY of the only witness the strand has, and
the mechanism split it forces.  Three findings:

  (A) ⇒⇒ THE `4 Z` OF [NCYL-307] H5 ARE **ONE SADDLE CONNECTION**, not four.  The two ROWS
      (`7/15 eps1`, `8/15 eps0`) carry byte-identical quotient-path data -- forced by the
      leg swap ([OPS-235]) and checked here rather than assumed -- and WITHIN a row the two
      `Z` prongs (edges `j = 2` and `j = 5`) are the SAME connection walked from its two
      ends, exact time-reverses under [NCYL-301]'s `R`.  ⇒ `4 = 1 object x 2 ends x 2
      forced rows`, and every statement in the strand about "when (G0b) fails" rests on
      `n = 1`.

  (B) ⇒⇒ THAT CONNECTION IS EXACTLY ONE LINEAR RELATION AMONG THE INTERFACE MEASURES.
      Carrying the walk symbolically in the `ell_t` (integer coefficients, no ring), the
      arrival condition at `7/15 eps1` reduces to

          ell_1 - ell_2 + 2 ell_3 - ell_4 + 2 ell_5 = 0,   i.e.
          2 sin(12deg) + 2 sin(24deg) + sin(36deg) = sin(60deg) + sin(72deg),

      an exact vanishing sum of sines of rational angles (`sympy.simplify -> 0`).  So the
      witness does NOT arrive by a formal identity: it CONSUMES an arithmetic coincidence.

  (C) ⇒⇒ AND THE COINCIDENCE IS A FUNCTION OF `Q` ALONE.  The `Q`-rank of the `Q`
      interface measures of a class is `phi(Q)/2` at every odd `Q` and in BOTH classes, so
      the space of integer relations among them has dimension

          LONE   (even `c`, the class with the horizontal interface):  (Q-1)/2 - phi(Q)/2
          PAIRED (odd  `c`, the class with the perpendicular one)    :  (Q+1)/2 - phi(Q)/2

      which is `0` **exactly on the LONE class at PRIME `Q`**.  ⇒ on those rows there is
      no relation to consume and mechanism (B) cannot occur at all; the witness lives at
      `Q = 15` because `3 | 15` supplies the three-term relation
      `sin(t) + sin(60deg - t) = sin(60deg + t)`.

  ⇒⇒ SO (G0b) SPLITS INTO TWO DISJOINT MECHANISMS, and this is the point of the file:
        (G0b-rel)   the arrival vector `a` is NONZERO -- a ring relation is consumed.
                    IMPOSSIBLE on the LONE class at prime `Q` (H3).  The only known
                    witness is of this kind (H2/H5).
        (G0b-form)  `a == 0` -- the arrival is a FORMAL integer identity, i.e. genuine
                    `Lambda_S` membership in COEFFICIENT space.  This is the half
                    [NCYL-298]'s criterion addresses and [NCYL-308] H3/H6 price as `87%` /
                    `77.6%` blind at corner starts.  UNTOUCHED here.
      ⚠⚠ THE SPLIT IS NOT A PROOF AND MUST NOT BE QUOTED AS ONE.  It closes (G0b-rel) on
      a named region and leaves (G0b-form) exactly where [NCYL-308] left it.

  ⚠⚠ AND THE DESIGN CONSTRAINT THE WITNESS IMPOSES ON ANY CANDIDATE, which is why the
  four routes [NCYL-308] priced all died for ONE reason: (G0b) is FALSE without a
  hypothesis -- the witness is a real, exact, terminating `Z`-`Z` connection on a real
  flat surface.  So no ROW-INDEPENDENT instrument can ever prove it (that kills the
  palindrome transfer and the topological/`chi` route a priori, before any measurement),
  and every candidate must be run on `7/15 eps1` FIRST and must ADMIT it.  H6 is that
  gate, applied to the three instruments the strand already has.

PRIOR ART: `rulings.py --grep` on 'vanishing sum', 'roots of unity', 'Conway', 'Mann',
'cyclotomic coincidence' -> NOTHING, in `rulings.md` or in any root `.md` including
`literature_review.md`: no probe or ruling in this repo has asked whether the interface
measures are linearly INDEPENDENT.  'composite Q' -> [NCYL-298] (whose own term-count
argument degenerates at `Q = 15`, "`Phi_15 mod 2` has degree 8 and exactly 8 terms", the
largest support that construction produces) and [NCYL-300] (which CLOSES that gap for
POLES by an exact `F_2[x]/Phi_Q` factorisation, no primality anywhere) -- so the
composite-`Q` worry is dead on the pole side and this is its corner-side twin.
[NCYL-308] prices four routes as not worth re-running: the plain `Lambda_S` test at corner
starts, the coefficient-sum augmentation, transferring [NCYL-301]'s palindrome without an
arrival-node argument, and a topological/`chi` attack.  ⚠ THIS IS NONE OF THEM.  Those
three arithmetic tests all ask MEMBERSHIP of a target in a lattice built from the `S_e`
INSIDE the ring, and therefore silently inherit whatever `Z`-relations the ring supplies;
this asks whether those relations EXIST.  The two questions coincide only where the
defect is `0`, which is the region (C) names and no prior arm isolated.

HYPOTHESES (pre-registered).
  H1  ⇒⇒ MULTIPLICITY, and it can fail on either half.  (a) the two witness rows have
      identical `(Lo, Hi, S)` as exact ring vectors; (b) within a row, the two `Z` prongs
      are one connection -- reverse the first orbit state-by-state under `R(i,d,y) =
      (i,-d,y)` and it must reproduce the second orbit EXACTLY, states and kinds.  If (b)
      failed they would be two connections and `n = 2`, not `n = 1`.
  H2  ⇒⇒ THE RELATION, exactly, with its own derivation guard.  G1 first: the symbolic
      endpoints (integer vectors over `ell_0..ell_{Q-1}`, built by mirroring `s451.build`)
      must reproduce `lo`/`hi` as exact ring vectors on every row in range -- a derivation
      check, not evidence.  Then for each corner prong ending at `Z`, extract the integer
      vector `a` with `sum_t a_t ell_t = 0`, verify the vanishing EXACTLY in
      `Z[zeta_{4Q}]`, and report whether `a == 0`.  ⚠ IT CAN COME OUT EITHER WAY, and the
      two answers say opposite things: `a != 0` puts the witness in (G0b-rel) and makes
      (C) load-bearing; `a == 0` puts it in (G0b-form), kills (C) as a mechanism, and says
      the whole residual is the lattice question [NCYL-308] already priced.
  H3  ⇒⇒ THE RANK LAW, over the whole box and both classes: `rank_Q{ell_t} == phi(Q)/2`,
      hence defect `(Q-1)/2 - phi(Q)/2` (LONE) / `(Q+1)/2 - phi(Q)/2` (PAIRED), hence
      `defect == 0 <=> (Q prime and LONE)`.  Three separate claims, each scoreable per
      row; a wrong rank, a wrong distinct-count or a single prime `Q` with a lone-class
      relation breaks it.
  H4  ⇒⇒ THE SHORT-RELATION CENSUS, and it is the DISCRIMINATION arm.  Shortest integer
      relation by `L1` at each composite `Q`.  ⚠ If every composite `Q` looked like
      `Q = 15` the mechanism would be necessary-and-insufficient in the least informative
      way; the census is what prices that, and the answer must be quoted with it.
  H5  ⇒⇒ THE SPLIT, SCORED ON REAL ROWS.  Over the corner census: every prong ending at
      `Z` classified `rel` (a != 0) or `form` (a == 0), and every prong ending at `Z`
      must sit on a row whose class has defect `> 0`.  ⚠ A `Z` on a prime-`Q` LONE row
      would refute (C) outright.  ⚠ `0 Z` on prime-`Q` LONE rows is also PREDICTED by
      `p = 0` being measured there, so it is COVERAGE, not independent support
      ([OPS-041]); what is not predicted is the `rel`/`form` classification.
  H6  ⇒⇒ THE PRE-FLIGHT GATE, applied to the instruments that exist.  On the actual
      witness edge, [NCYL-298]'s lattice criterion, [NCYL-308] H6's augmented criterion
      and [NCYL-301]'s Step-2 test must each ADMIT the connection.  ⚠ A pass is PREDICTED
      -- a real connection that a necessary condition refutes is a BUG, not a discovery --
      so this is a bug-check and a reusable gate, never evidence.  It is kept because it
      is the cheapest thing a future candidate can be run through before it is built.

GUARDS.
  G1  a row whose global coordinate does not close (`s451.build`) is SKIPPED, not scored.
  G2  a TIE edge / degenerate node carries no corner prong; SKIPPED, never a pass.
  G3  nothing here writes to any store, and no existing probe is modified.
  G4  ⚠ [OPS-228]: numpy `int64` is cast to Python `int` before `zlattice`/`sympy`.
  G5  the symbolic layer is checked against the ring layer (H2's derivation guard) before
      any `a` is reported; a mismatch aborts the row rather than reporting a relation.

Run: PYTHONPATH=.:engine:archive/scripts_2026-07:archive/scripts_2026-08:probes \\
     .venv/bin/python3.13 probes/s462_zz_relation.py [all|h1|h2|h3|h4|h5|h6] [--qmax=N]
"""
from __future__ import annotations

import itertools
import json
import math
import sys
import time

import numpy as np
from sympy import Matrix

import s450_chain_maps as cm
import s451_quotient_path as qp
import s453_pole_module as pm
import s454_self_hit as sh
import s456_neckgb as nb
import zlattice as zl

OUT = "data/s462_zz_relation.json"
WITNESS = [(7, 15, 1), (8, 15, 0)]


# ---------------------------------------------------------------- number theory helpers

def phi(n):
    r, m, p = n, n, 2
    while p * p <= m:
        if m % p == 0:
            while m % p == 0:
                m //= p
            r -= r // p
        p += 1
    if m > 1:
        r -= r // m
    return r


def isprime(n):
    return n > 1 and all(n % d for d in range(2, int(n ** 0.5) + 1))


def _rows(qmax, qmin=5):
    for Q in range(qmin, qmax + 1, 2):
        for P in range(1, Q):
            if math.gcd(P, Q) == 1:
                for eps in (0, 1):
                    yield P, Q, eps


def _chain(P, Q, eps):
    ch, lo, hi, closes = qp.build(P, Q, eps)
    if not closes:
        return None
    sf = pm.fold_sums(ch, lo, hi)
    Lo, Hi, S = pm.path_view(ch, lo, hi, sf)
    return ch, lo, hi, Lo, Hi, S, pm.CoordOps(ch)


# ---------------------------------------------------------------- the symbolic layer

def sym_build(ch):
    """`lo`/`hi` as INTEGER VECTORS over the basis `(ell_0, .., ell_{Q-1})`.

    Mirrors `s451.build` line for line, with `Coord` arithmetic replaced by `Z^Q` vectors
    and `ch.delta[k] = ell_far - ell_near` replaced by `e_far - e_near`.  Nothing here
    touches the ring, which is what makes it an independent expression of the SAME
    quantity -- `check_sym` is the derivation guard that they agree.
    """
    Q = ch.Q
    z = lambda: [0] * Q                                        # noqa: E731
    U = {0: z()}
    for k in range(1, Q):
        far, near = ch.far[k], ch.near[k]
        d = z()
        d[far] += 1
        d[near] -= 1
        if far in U:
            U[near] = [a + b for a, b in zip(U[far], d)] if ch.osec[k] else list(U[far])
        else:
            U[far] = [a - b for a, b in zip(U[near], d)] if ch.osec[k] else list(U[near])
    LO = {t: U[t] for t in range(Q)}
    HI = {}
    for t in range(Q):
        v = list(U[t])
        v[t] += 1
        HI[t] = v
    return LO, HI


def sym_value(ch, a):
    """`sum_t a_t * ell_t` as an exact ring vector (the `Coord.v` convention: `2*x`)."""
    out = np.zeros(ch.ring.d, dtype=np.int64)
    for t, c in enumerate(a):
        if c:
            out = out + c * ch.ell[t].v
    return out


def check_sym(ch, lo, hi, LO, HI):
    """G5 / H2's derivation guard: the symbolic endpoints reproduce the ring ones."""
    for t in range(ch.Q):
        if not np.array_equal(sym_value(ch, LO[t]), lo[t].v):
            return False, ("lo", t)
        if not np.array_equal(sym_value(ch, HI[t]), hi[t].v):
            return False, ("hi", t)
    return True, None


def sym_path(ch, LO, HI, lo, hi):
    """The path view of the symbolic endpoints and fold sums (mirrors `pm.path_view`
    and `pm.fold_sums`; the lo/hi ring values decide WHICH endpoint is shared)."""
    Q, m = ch.Q, ch.m_idx
    sLo = [LO[m - i] for i in range(m + 1)]
    sHi = [HI[m - i] for i in range(m + 1)]
    sS = []
    for j in range(m):
        k = m - j
        a, b = (k - 1) % Q, k
        if lo[a].eq(lo[b]):
            sS.append([x + y for x, y in zip(HI[a], HI[b])])
        else:
            sS.append([x + y for x, y in zip(LO[a], LO[b])])
    return sLo, sHi, sS


def sym_walk(ch, Lo, Hi, S, sLo, sHi, sS, st0, ops, cap=200_000):
    """Walk `sh.orbit` and carry the symbolic `y` alongside.

    Returns `(end, a, folds, steps)`.  `a` is the integer relation vector at a `Z`
    arrival (`sym(y) - sym(endpoint)`), else `None`.
    """
    i, d, y = st0
    # locate the symbolic twin of the starting coordinate
    ysym = None
    for idx in range(len(Lo)):
        if ops.eq(y, Lo[idx]):
            ysym = list(sLo[idx])
            break
        if ops.eq(y, Hi[idx]):
            ysym = list(sHi[idx])
            break
    if ysym is None:
        return "NOSYM", None, 0, 0
    folds = 0
    for step in range(cap):
        i2 = i + d
        if i2 < 0 or i2 >= len(Lo):
            return "ESC", None, folds, step
        c0, c1 = ops.cmp(y, Lo[i2]), ops.cmp(y, Hi[i2])
        if c0 == 0 or c1 == 0:
            tgt = sLo[i2] if c0 == 0 else sHi[i2]
            return "Z", [x - t for x, t in zip(ysym, tgt)], folds, step
        if c0 > 0 and c1 < 0:
            i = i2
            continue
        e = min(i, i2)
        if ops.eq(y, ops.half(S[e])):
            return "R", None, folds, step
        y = ops.sub(S[e], y)
        ysym = [s - x for s, x in zip(sS[e], ysym)]
        d = -d
        folds += 1
    return "CAP", None, folds, cap


# ---------------------------------------------------------------- H1

def h1():
    """MULTIPLICITY: the `4 Z` are one saddle connection."""
    res = {"rows": [], "rows_identical": None, "reversal_ok": None, "objects": None}
    packs = []
    for (P, Q, eps) in WITNESS:
        c = _chain(P, Q, eps)
        assert c is not None, (P, Q, eps)
        ch, lo, hi, Lo, Hi, S, ops = c
        packs.append((P, Q, eps, ch, Lo, Hi, S, ops))
    # (a) the two rows are the same system.  ⚠ NOT as identical `Lo`/`Hi`: the leg swap
    # acts by the `Z_O <-> Z_A` relabelling, i.e. the per-interface flip `x -> ell - x`,
    # which is NOT a global affine map of the path coordinate (checked: neither
    # `y -> y + c` nor `y -> -y + c` conjugates them).  What IS identical is the node
    # LENGTHS as exact ring vectors, and the whole per-edge corner-prong readout.
    (_p0, _q0, _e0, _c0, L0, H0, S0, o0) = packs[0]
    (_p1, _q1, _e1, _c1, L1, H1v, S1, o1) = packs[1]
    lens_same = (len(L0) == len(L1) and all(
        np.array_equal(o0.sub(b, a).v, o1.sub(d, c).v)
        for a, b, c, d in zip(L0, H0, L1, H1v)))

    def readout(L, H, Sx, ops):
        tab = []
        for j in range(len(Sx)):
            cs = nb.corner_start(L, H, Sx, j, ops)
            if cs is None:
                tab.append((j, "SKIP", 0))
                continue
            _s, kinds, _e, end = sh.orbit(L, H, Sx, cs[0], ops, 200_000)
            tab.append((j, end, len(kinds)))
        return tab

    t0_, t1_ = readout(L0, H0, S0, o0), readout(L1, H1v, S1, o1)
    res["node_lengths_identical"] = bool(lens_same)
    res["readouts_identical"] = bool(t0_ == t1_)
    same = lens_same and t0_ == t1_
    res["rows_identical"] = bool(same)

    # (b) within row 0, the Z prongs are one connection: reverse one, get the other
    zprongs = []
    for j in range(len(S0)):
        cs = nb.corner_start(L0, H0, S0, j, o0)
        if cs is None:
            continue
        sts, kinds, _e, end = sh.orbit(L0, H0, S0, cs[0], o0, 200_000)
        res["rows"].append({"j": j, "end": end, "steps": len(kinds)})
        if end == "Z":
            zprongs.append((j, sts, kinds))
    ok = None
    if len(zprongs) == 2:
        (_ja, sa, ka), (_jb, sb, kb) = zprongs
        # R(u_{M-k}) must be the k-th state of the other orbit
        rev = [(i, -d, y) for (i, d, y) in reversed(sa)]
        ok = (len(rev) == len(sb)
              and all(sh._eqstate(x, y, o0) for x, y in zip(rev, sb)))
        # and the kind words must be reverses of each other (fold/cross, minus the end)
        ok = bool(ok and ka[:-1][::-1] == kb[:-1])
    res["n_zprongs"] = len(zprongs)
    res["reversal_ok"] = ok
    res["objects"] = 1 if (same and ok) else None
    return res


# ---------------------------------------------------------------- H2

def h2(qmax=33):
    """THE RELATION.  Derivation guard first, then the arrival vector at every `Z`."""
    res = {"rows": 0, "guard_ok": 0, "guard_bad": [], "z_hits": [],
           "n_rel": 0, "n_form": 0}
    for P, Q, eps in _rows(qmax):
        c = _chain(P, Q, eps)
        if c is None:
            continue
        ch, lo, hi, Lo, Hi, S, ops = c
        LO, HI = sym_build(ch)
        good, where = check_sym(ch, lo, hi, LO, HI)
        res["rows"] += 1
        if not good:
            res["guard_bad"].append([P, Q, eps, list(where)])
            continue
        res["guard_ok"] += 1
        sLo, sHi, sS = sym_path(ch, LO, HI, lo, hi)
        for j in range(len(S)):
            cs = nb.corner_start(Lo, Hi, S, j, ops)
            if cs is None:
                continue
            end, a, folds, steps = sym_walk(ch, Lo, Hi, S, sLo, sHi, sS, cs[0], ops)
            if end != "Z":
                continue
            vanishes = not sym_value(ch, a).any()
            nonzero = any(a)
            res["n_rel" if nonzero else "n_form"] += 1
            res["z_hits"].append({"P": P, "Q": Q, "eps": eps, "j": j,
                                  "folds": folds, "steps": steps,
                                  "a": [int(x) for x in a],
                                  "vanishes_exactly": bool(vanishes),
                                  "kind": "rel" if nonzero else "form"})
    return res


# ---------------------------------------------------------------- H3

def ell_values(Q, lone):
    """The DISTINCT nonzero interface measures of a class, keyed by `min(c, 2Q-c)`.

    `lone` <=> the class whose `c` are EVEN <=> it contains the horizontal interface
    (`ell = 0`); `paired` <=> `c` odd <=> it contains the perpendicular one (`ell = 1`).
    Using `P = 1`, since `c_t = c_0 + 2Pt` covers every residue of one parity for any
    coprime `P` -- the SET depends only on that parity ([NCYL-291]).
    """
    ch = cm.Chain(1, Q, 1 if lone else 0)        # P=1: eps=1 -> c even, eps=0 -> c odd
    out = {}
    for t in range(Q):
        c, v = int(ch.c[t]), ch.ell[t].v
        if v.any():
            out.setdefault(min(c, 2 * Q - c), v)
    return out


def h3(qmax=59):
    """THE RANK LAW and the defect formula."""
    res = {"checked": 0, "rank_ok": 0, "count_ok": 0, "defect_ok": 0,
           "zero_defect_rows": [], "bad": [], "table": []}
    for Q in range(5, qmax + 1, 2):
        for lone in (True, False):
            d = ell_values(Q, lone)
            ks = sorted(d)
            M = Matrix([[int(x) for x in d[k]] for k in ks])
            rk = int(M.rank())
            want_rk = phi(Q) // 2
            want_n = (Q - 1) // 2 if lone else (Q + 1) // 2
            defect = len(ks) - rk
            want_def = want_n - want_rk
            res["checked"] += 1
            res["rank_ok"] += (rk == want_rk)
            res["count_ok"] += (len(ks) == want_n)
            res["defect_ok"] += (defect == want_def)
            if rk != want_rk or len(ks) != want_n:
                res["bad"].append([Q, "lone" if lone else "paired", rk, want_rk,
                                   len(ks), want_n])
            if defect == 0:
                res["zero_defect_rows"].append([Q, "lone" if lone else "paired"])
            res["table"].append({"Q": Q, "class": "lone" if lone else "paired",
                                 "n_distinct": len(ks), "rank": rk, "defect": defect,
                                 "prime": isprime(Q)})
    zd = res["zero_defect_rows"]
    res["zero_defect_iff_prime_lone"] = all(isprime(q) and c == "lone" for q, c in zd) \
        and all(any(q == Q and c == "lone" for q, c in zd)
                for Q in range(5, qmax + 1, 2) if isprime(Q))
    return res


# ---------------------------------------------------------------- H4

def shortest_relation(Q, lone, span=2):
    d = ell_values(Q, lone)
    ks = sorted(d)
    M = Matrix([[int(x) for x in d[k]] for k in ks])
    ns = M.T.nullspace()
    if not ns:
        return 0, None, ks
    B = []
    for v in ns:
        den = 1
        for x in v:
            den = math.lcm(den, int(x.q))
        B.append([int(x * den) for x in v])
    B = Matrix(B)
    r = min(B.rows, 6)
    best = None
    for co in itertools.product(range(-span, span + 1), repeat=r):
        if not any(co):
            continue
        w = [int(sum(co[i] * B[i, j] for i in range(r))) for j in range(B.cols)]
        g = 0
        for x in w:
            g = math.gcd(g, x)
        if g > 1:
            w = [x // g for x in w]
        l1 = sum(abs(x) for x in w)
        if l1 and (best is None or l1 < best[0]):
            best = (l1, w)
    return len(ns), best, ks


def h4(qmax=45):
    """THE SHORT-RELATION CENSUS -- and its non-sufficiency, stated with it."""
    out = {"rows": []}
    for Q in range(5, qmax + 1, 2):
        for lone in (True, False):
            defect, best, ks = shortest_relation(Q, lone)
            rec = {"Q": Q, "class": "lone" if lone else "paired",
                   "prime": isprime(Q), "defect": defect,
                   "three_divides_Q": Q % 3 == 0}
            if best:
                l1, w = best
                rec["minL1"] = l1
                rec["relation"] = [[int(c), int(k)] for c, k in zip(w, ks) if c]
            out["rows"].append(rec)
    lone3 = [r for r in out["rows"] if r["class"] == "lone" and r["defect"]]
    out["lone_minL1_is_3_iff_3divQ"] = all(
        (r.get("minL1") == 3) == r["three_divides_Q"] for r in lone3)
    return out


# ---------------------------------------------------------------- H5

def h5(qmax=33):
    """THE SPLIT, on real rows: where the `Z` prongs are and which mechanism they use."""
    r2 = h2(qmax)
    res = {"z_hits": r2["z_hits"], "n_rel": r2["n_rel"], "n_form": r2["n_form"],
           "rows_scored": r2["guard_ok"], "on_zero_defect": [], "by_Q": {}}
    for hit in r2["z_hits"]:
        Q, P, eps = hit["Q"], hit["P"], hit["eps"]
        lone = (eps == P % 2)
        defect = ((Q - 1) // 2 if lone else (Q + 1) // 2) - phi(Q) // 2
        hit["class"] = "lone" if lone else "paired"
        hit["defect"] = defect
        res["by_Q"][str(Q)] = res["by_Q"].get(str(Q), 0) + 1
        if defect == 0:
            res["on_zero_defect"].append(hit)
    res["all_Z_have_a_relation_available"] = not res["on_zero_defect"]
    res["all_Z_are_relation_consuming"] = (r2["n_form"] == 0 and r2["n_rel"] > 0)
    return res


# ---------------------------------------------------------------- H6

def h6():
    """THE PRE-FLIGHT GATE: the three existing instruments must ADMIT the real witness."""
    out = {"cases": []}
    for (P, Q, eps) in WITNESS:
        c = _chain(P, Q, eps)
        ch, lo, hi, Lo, Hi, S, ops = c
        m = len(S)
        plain = [[int(x) for x in S[e].v] for e in range(m)]
        Lpln = zl.Lattice(plain)
        Laug = zl.Lattice([g + [1] for g in plain])
        for j in range(m):
            cs = nb.corner_start(Lo, Hi, S, j, ops)
            if cs is None:
                continue
            sts, kinds, _e, end = sh.orbit(Lo, Hi, S, cs[0], ops, 200_000)
            if end != "Z":
                continue
            last = sts[-1]
            i2 = last[0] + last[1]
            tgt = Lo[i2] if ops.eq(last[2], Lo[i2]) else Hi[i2]
            y0 = [int(x) for x in cs[0][2].v]
            ev = [int(x) for x in tgt.v]
            # [NCYL-298]: target -/+ y0 in Lambda_S
            adm_plain = (Lpln.contains([a - b for a, b in zip(ev, y0)])
                         or Lpln.contains([a + b for a, b in zip(ev, y0)]))
            # [NCYL-308] H6: the coefficient-sum augmentation
            adm_aug = (Laug.contains([a - b for a, b in zip(ev, y0)] + [0])
                       or Laug.contains([a + b for a, b in zip(ev, y0)] + [1]))
            # [NCYL-301] Step 2: is the corner an endpoint of BOTH intervals at its edge?
            lo_s = ops.eq(Lo[j], Lo[j + 1])
            shared = Lo[j] if lo_s else Hi[j]
            both = sum(1 for i in (j, j + 1)
                       if ops.eq(shared, Lo[i]) or ops.eq(shared, Hi[i]))
            out["cases"].append({"P": P, "Q": Q, "eps": eps, "j": j,
                                 "NCYL298_admits": bool(adm_plain),
                                 "NCYL308H6_admits": bool(adm_aug),
                                 "NCYL301_step2_arrival_forced": both != 2})
    out["all_admit"] = all(c["NCYL298_admits"] and c["NCYL308H6_admits"]
                           and not c["NCYL301_step2_arrival_forced"]
                           for c in out["cases"])
    return out


# ---------------------------------------------------------------- driver

def main():
    which = sys.argv[1] if len(sys.argv) > 1 and not sys.argv[1].startswith("-") else "all"
    qmax = None
    for a in sys.argv[1:]:
        if a.startswith("--qmax="):
            qmax = int(a.split("=")[1])
    out, t0 = {}, time.time()

    if which in ("all", "h1"):
        t = time.time()
        r = out["h1"] = h1()
        print(f"  H1  the two witness ROWS: node lengths identical  : "
              f"{r['node_lengths_identical']}", flush=True)
        print(f"      ... and per-edge corner readout identical      : "
              f"{r['readouts_identical']}"
              f"   (the leg swap, [OPS-235] -- checked, not assumed)", flush=True)
        print(f"      Z prongs in 7/15 eps1                          : {r['n_zprongs']}",
              flush=True)
        print(f"      they are exact time-reverses under R           : {r['reversal_ok']}",
              flush=True)
        print(f"      ⇒ distinct saddle connections behind the `4 Z` : {r['objects']}",
              flush=True)
        print(f"      per-prong endings: "
              + ", ".join(f"j{d['j']}:{d['end']}" for d in r["rows"]), flush=True)
        print(f"      [{time.time() - t:.1f}s]", flush=True)

    if which in ("all", "h3"):
        t = time.time()
        r = out["h3"] = h3(qmax or 59)
        n = r["checked"]
        print(f"  H3  rank == phi(Q)/2                     : {r['rank_ok']}/{n}"
              f"  (odd Q <= {qmax or 59}, both classes)", flush=True)
        print(f"      #distinct nonzero == (Q-+1)/2        : {r['count_ok']}/{n}",
              flush=True)
        print(f"      defect == the formula                : {r['defect_ok']}/{n}",
              flush=True)
        print(f"      ⇒ defect == 0  <=>  Q PRIME and LONE : "
              f"{r['zero_defect_iff_prime_lone']}"
              f"   ({len(r['zero_defect_rows'])} zero-defect classes)", flush=True)
        for b in r["bad"][:6]:
            print("      bad:", b, flush=True)
        print(f"      [{time.time() - t:.1f}s]", flush=True)

    if which in ("all", "h4"):
        t = time.time()
        r = out["h4"] = h4(min(qmax or 45, 45))
        print("  H4  shortest integer relation, LONE class "
              "(s(k) := sin(k*pi/(2Q))):", flush=True)
        for rec in r["rows"]:
            if rec["class"] != "lone":
                continue
            if not rec["defect"]:
                print(f"      Q={rec['Q']:3d} {'prime':9s} defect 0  -- NO RELATION",
                      flush=True)
            else:
                trm = " ".join(f"{c:+d}s({k})" for c, k in rec["relation"])
                print(f"      Q={rec['Q']:3d} {'COMPOSITE':9s} defect {rec['defect']}"
                      f"  minL1={rec['minL1']:3d}   {trm} = 0", flush=True)
        print(f"      ⇒ minL1 == 3 exactly when 3 | Q      : "
              f"{r['lone_minL1_is_3_iff_3divQ']}"
              f"   (the sin(t)+sin(60-t)=sin(60+t) family)", flush=True)
        print("      ⚠ NECESSARY, NOT SUFFICIENT: Q = 9,21,27,33,39 carry the same "
              "three-term relation and NO witness", flush=True)
        print(f"      [{time.time() - t:.1f}s]", flush=True)

    if which in ("all", "h2", "h5"):
        t = time.time()
        r = out["h5"] = h5(qmax or 33)
        print(f"  H2/H5  symbolic derivation guard passed on {r['rows_scored']} rows",
              flush=True)
        print(f"      corner prongs ending at Z: {r['n_rel']} RELATION-CONSUMING (a != 0),"
              f" {r['n_form']} FORMAL (a == 0)", flush=True)
        print(f"      by Q: {r['by_Q']}", flush=True)
        for h in r["z_hits"]:
            print(f"      Z: {h['P']}/{h['Q']} eps{h['eps']} j={h['j']} "
                  f"{h['class']} defect={h['defect']} folds={h['folds']} "
                  f"kind={h['kind']} vanishes_exactly={h['vanishes_exactly']}",
                  flush=True)
            print(f"         a = {h['a']}", flush=True)
        print(f"      ⇒ every Z sits where a relation EXISTS : "
              f"{r['all_Z_have_a_relation_available']}", flush=True)
        print(f"      ⇒ every Z is relation-consuming        : "
              f"{r['all_Z_are_relation_consuming']}", flush=True)
        print(f"      [{time.time() - t:.1f}s]", flush=True)

    if which in ("all", "h6"):
        t = time.time()
        r = out["h6"] = h6()
        for c in r["cases"]:
            print(f"  H6  {c['P']}/{c['Q']} eps{c['eps']} j={c['j']}:"
                  f"  [NCYL-298] admits {c['NCYL298_admits']},"
                  f"  [NCYL-308]H6 admits {c['NCYL308H6_admits']},"
                  f"  [NCYL-301] Step2 forces arrival"
                  f" {c['NCYL301_step2_arrival_forced']}", flush=True)
        print(f"      ⇒ all three ADMIT the real witness (a PASS is predicted; a refusal "
              f"would be a bug): {r['all_admit']}", flush=True)
        print(f"      [{time.time() - t:.1f}s]", flush=True)

    out["secs"] = round(time.time() - t0, 1)
    json.dump(out, open(OUT, "w"), indent=1)
    print(f"  -> {OUT}  [{out['secs']}s]", flush=True)


if __name__ == "__main__":
    main()
