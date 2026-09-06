"""s340_glen_enum.py -- PRE-REGISTERED: a `g_len`-CAPABLE INTERVAL ENUMERATOR.

THE QUESTION (the handoff's NEXT (a), since s338, PROMOTED s339).  Every remaining
alpha->0 ladder measurement is depth-walled: `(3,1)`/`(3,2)`/`(5,3)`/`(5,4)` guard-fail on
depth by Q ~ 30-44, and sec.s339 (8) showed the wall is set by the FAN's max T/chi -- the
same object the general law is about.  So the depth wall and the general law are ONE
obstruction, and the enumerator is the critical path rather than an engineering
convenience.  [OPS-050] states the build exactly: `s312_leaf_enum` is the right TECHNIQUE
(exact interval propagation, no grid, no ring, no depth wall, 234x the ring per partition)
but is scoped to the LEAF CLASS -- its cells sum to Sum w = 0.9036 at [5,2,3], and its
records carry no Rx_c and no head-on flag, so no T_c and no chi_c.  A g_len-capable
version "must stop pruning at Lemma-B violations (to reach seam cells) and record Rx --
AND THE PRUNE IS WHAT MAKES IT CHEAP".

WHAT THIS PROBE ADDS TO `s312_leaf_enum.enum_leaf`, and nothing else.
  (i)   THE DEVELOPED x.  The escape is the POINT reflection through the R vertex, so the
        O-corner's developed x obeys the exact twin of the height rule already in the
        model: alongside c' = 2 y_R - c, carry x' = 2 x_R - x with
        x_R = x + cot*cos(Omega).  Launch at x = 0, whose i = 0 R vertex sits at
        x = cot(alpha) -- which is sec.s275's launch abscissa.
  (ii)  LEMMA B's OTHER BRANCH instead of the prune.  At a subinterval whose chain
        violates the opposite-sides hypothesis, [WFLOOR-113]'s branch (Phi' = Omega + pi,
        eps' = +eps, j0' = 2 - (j* mod 2)) continues it rather than dropping it, so the
        seam cells are reached and the cells can sum to 1.
  (iii) A GROWABLE EXIT WINDOW.  The other branch REINFORCES the winding ([WFLOOR-113]),
        so an exit index may run past `s312_leaf_enum`'s fixed imax = 4a1+6; the exit scan
        grows on demand and splits only at vertices it actually consults.

THE MEASUREMENT THAT MAKES chi_c UNNECESSARY (recon, this session, and it is why there is
no orphan detector below).  Against sec.s275's ring cells, the model's FORWARD developed
displacement equals the ring's FIRST-RETURN length cell for cell:
        T_c / chi_c  ==  1 - x_{j*}/cot(alpha),
paired cells (chi = 2) and the orphan (chi = 1) alike -- 54/54 on the odd-P rows of the
recon.  The ring's T_c is the FULL PERIOD, and the model's word is the palindrome
fwd . reverse(fwd); on a paired cell the ring depth is the model's whole closure 2n+1 and
T_c is twice the forward displacement, while on the ORPHAN the ring depth is the model's
forward n and T_c is the forward displacement itself.  The factor 2 that distinguishes
them is exactly chi_c, so it CANCELS: g_len needs the forward half and never needs to know
which cell is the orphan.  chi is still reported (from closure vs the ring depth) but is
not on the path to S.

HYPOTHESES, FIXED BEFORE THE RUN.
  G1  T_c/chi_c == 1 - x_{j*}/cot on every cell of every ring-certified ODD-P row.
                                        PREDICT: holds, ~1e-12 relative.  (Recon: 54/54.)
  G2  With the other branch replacing the prune, the enumerator's cells sum to
      Sum w = 1 and reproduce the ring partition (count and boundaries).
                                        PREDICT: THIS IS THE ONE AT RISK.  sec.s313 (9)
                                        measured `full_trace` at 124/160 on out-of-domain
                                        closure lengths, "cause undiagnosed", and that is
                                        the same branch.  A shortfall here is the
                                        expected outcome and is the thing to diagnose;
                                        Sum w != 1 is the detector.
  G3  g_len from the enumerator matches the ring's on certified rows.
                                        PREDICT: holds wherever G2 holds, and NOWHERE
                                        else -- G2 is the gate, G3 is not independent.
  G4  Cost tracks total letters, not a depth cap, so rows the ring guard-fails are
      reachable.                        PREDICT: reachable; NO PREDICTION on the constant.
  G5  SCOPE.  The head-on test is the congruence (phi + eps*j*P) = 0 mod 2Q, so an
      IN-FAN head-on exists as soon as 2Q/gcd(P,2Q) < imax -- at even P, gcd = 2 and that
      index is Q.  PREDICT: the model is degenerate exactly on those rows and correct off
      them, so the criterion is NOT "P even": 2/7, 4/15, 2/9, 4/17 fail and 6/19, 6/23,
      8/17 pass.  (Recon: the 4 in-range rows are exactly the 4 that failed, 13/13 vs
      54/54.)  Every alpha->0 ladder of sec.s338/sec.s339 is odd P.

PRIOR ART: grepped rulings.md (all tags), TOOLS.md, cf_width_laws.md, session_state.md for
'g_len', 'enum_leaf', 'interval propagation', 'seam cell', 'Rx', 'other branch', 'prune',
'partition reader' ->
  - [OPS-050] (s338) OWNS "no fast float route to g_len exists" and STATES this build's
    two requirements verbatim; it is the reason this probe exists, not a bar on it.
  - [WFLOOR-111] (s312) OWNS the interval-propagation enumerator and its LEAF-CLASS scope;
    [WFLOOR-113] (s313) OWNS Lemma B's other branch; sec.s313 (9) OWNS the 124/160
    out-of-domain closure defect and records the cause as UNDIAGNOSED -- the open item
    this probe walks into, so G2 is pre-registered as the one at risk.
  - sec.s303 (7b) OWNS the C1 refutation ("the model reproduces the ring on every cell"),
    204/444, and is NOT overturned by anything here.
  - [WFLOOR-107]/[WFLOOR-109]/[WFLOOR-112] own the opposite-sides domain and its PROVED
    identification with {G_B <= a1}; s275 OWNS the ring quantities T_c, chi_c, the flow
    identity and the two completeness guards ([NCYL-019]).
  - [NCYL-020] owns "a truncated enumeration is not a missing orbit"; the Sum w = 1 guard
    below is its analogue for this instrument.
  - Nothing found that carries the DEVELOPED x through the fan chain: the model tracks
    heights only, and no probe reads a length off it.

WHAT THIS CANNOT DO (up front).  It is mpmath FLOAT against the ring's exact cyclotomic
arithmetic -- an EXPLORATORY instrument, never a certificate ([NCYL-019] is the
certificate, and TOOLS.md's FAN MODEL trap (ii) applies verbatim).  It moves no DAG node
and touches neither [C.10] nor [O-DEFICIT] nor gamma=1: every row it can reach has
g_len < 1, so it bears on the DEFICIT and never on `lone + CP => full sweep`.

VERDICT (s340; the hypotheses above are left exactly as registered).
  G1  HELD.  T_c/chi_c == 1 - x_{j*}/cot on 51/51 ring-certified odd-P rows, worst
      relative error 3.5e-16.
  G2  HELD, AND IT WAS THE ONE PREDICTED TO FAIL.  Sum_c w_c - 1 == 0.0 exactly on all
      51 rows, cell counts equal 51/51, boundaries to 2.7e-13.  sec.s313 (9)'s 124/160
      out-of-domain shortfall did NOT reappear.  ⚠ NOT DIAGNOSED -- three differences
      from `full_trace` are candidates (a growable exit window vs its fixed 4a1+8; the
      interval decision read off ENDPOINTS rather than a midpoint; the graze-sliver
      floor), and this session did not separate them.  Do NOT read G2 as overturning
      sec.s303 (7b) or as repairing `full_trace`.
  G3  HELD.  worst |g_len_model - g_len_ring| = 7.3e-14 over the 51 rows.
  G4  HELD IN KIND, MUCH WEAKER IN DEGREE THAN THE QUEUE ASSUMED -- see sec.s340.  The
      ring's 200000-letter cap is gone (complete partitions returned at 869233, 563681,
      532949, 368501, 286553 and 251789 letters), and a head-to-head at that boundary is
      1853.6 s for the ring at 3/31 (cap raised to 1.6e6) against 43.6 s for the
      enumerator at 3/34, one rung further.  But cost is linear in TOTAL LETTERS and the
      letters are GEOMETRIC in m (ratio ~3.5 on (3,1)), so a budget buys rungs only as
      log(budget): 18 new rungs over 7 ladders (0 to +5 each) and not a regime.
  G5  REFUTED, by its own pre-registered discriminators.  6/19, 6/23 and 8/17 are even P
      with 2Q/gcd(P,2Q) >= imax -- predicted CLEAN -- and all three break (10/18, 10/22,
      6/16 cells).  The criterion is not the head-on congruence; it is simply P EVEN, and
      the mechanism is UNKNOWN.  Scope bar: ODD P only.  Every alpha->0 ladder is odd P,
      so the strand is unaffected, but 10/23 and 20/23 are NOT reachable by this tool.
      ⚠⚠ SUPERSEDED s341 ([WFLOOR-115]) -- THE SCOPE BAR IS LIFTED AND THE MECHANISM IS
      KNOWN.  The head-on test was at the EXIT vertex where retrace needs it at the R
      VERTEX (see `enum_glen` below); the two agree exactly when j* is even, and at odd P
      the phase subgroup forces that, which is why 51/51 rows could not see it.  Fixed
      here: 40/40 ring-certified even-P rows now exact, odd P unchanged at 51/51 figure
      for figure, and 10/23 IS now reachable (true closure 12474905).  20/23 remains out
      of reach on the SECOND, independent bar a1 >= 2.  G5's own congruence criterion
      stays REFUTED -- it is the right object asked the wrong question (the PARITY of the
      index that fires, not whether the least positive solution falls inside the window).
      Regression harness: `probes/s341_evenp_headon.py`.

Run: PYTHONPATH=.:engine:archive/scripts_2026-07:archive/scripts_2026-08:probes \
     .venv/bin/python3.13 -u probes/s340_glen_enum.py [check|scope|reach|dps|ladder P r m...]
  check   G1/G2/G3 against every ring-certified row of a small box
  scope   G5, the head-on degeneracy criterion
  reach   G4, how far each alpha->0 ladder goes under a per-row wall-clock budget
  dps     the float self-check: re-run the deepest rows at dps 80 against dps 50
  ladder  run one (P,r) ladder past the ring's depth wall
"""
import json
import os
import sys
import time
from math import gcd

import mpmath as mp

from s302_model import base

mp.mp.dps = 50
OUT = 'data/s340_glen_enum.json'
TINY = mp.mpf(10) ** -25          # graze-ray width floor, as s312_leaf_enum


def headon_index(P, Q):
    """Least j > 0 with j*P = 0 (mod 2Q) from phi = 0 -- an IN-FAN head-on if it is
    smaller than the exit window.  G5's criterion."""
    return (2 * Q) // gcd(P, 2 * Q)


# --------------------------------------------------------------------------
# RESUMABLE ENUMERATION (s354).  A node is a COMPLETE state vector -- nothing
# outside it is consulted except `base(P, r, a1)`, which is a pure function of
# the centre -- so a node drained by a cap can be stored and restarted later and
# the continuation is BIT-IDENTICAL to never having stopped.  Before s354 the
# drain path was `capped += 1; continue`, a counter: the interval was thrown
# away, so the only way to push a capped centre further was to redo the whole
# transversal from the root.
#
# ⚠ Serialisation goes through mpmath's `_mpf_` limb tuple, NOT `str`/`float`.
# A decimal round-trip at dps 50 perturbs the interval endpoints, and an endpoint
# is a DECISION boundary here -- the enumerator splits on `lo < y < hi`.
def pack_node(nd):
    """(lo, hi, c, x, phi, eps, j0, n, nfan, nopp, v, gmax, hmax, wok) -> JSON-able.

    ⚠ THE FOUR WINDING FIELDS ARE s383 AND THEY ARE PURE PASSENGERS: no enumeration
    DECISION reads them, so a record packed before s383 resumes to a bit-identical
    partition.  What it cannot resume is the winding HISTORY -- `unpack_node` marks
    such a node `wok=False` and every cell descending from it reports `gb`/`hmax`/
    `e_end` as None rather than a walk restarted from 0 mid-chain.
    """
    lo, hi, c, x, phi, eps, j0, n, nfan, nopp, v, gmax, hmax, wok = nd
    return dict(lo=list(lo._mpf_), hi=list(hi._mpf_), c=list(c._mpf_),
                x=list(x._mpf_), phi=int(phi), eps=int(eps), j0=int(j0),
                n=int(n), nfan=int(nfan), nopp=int(nopp),
                v=int(v), gmax=int(gmax), hmax=int(hmax), wok=bool(wok))


def unpack_node(d):
    """Inverse of `pack_node`.  A PRE-s383 packed node carries no winding history,
    so it comes back `wok=False` -- see `pack_node`."""
    return (mp.make_mpf(tuple(d['lo'])), mp.make_mpf(tuple(d['hi'])),
            mp.make_mpf(tuple(d['c'])), mp.make_mpf(tuple(d['x'])),
            d['phi'], d['eps'], d['j0'], d['n'], d['nfan'], d['nopp'],
            d.get('v', 0), d.get('gmax', 0), d.get('hmax', 0),
            bool(d.get('wok', False)))


def node_width(d):
    """Width of a packed node's interval, as a float."""
    return float(mp.make_mpf(tuple(d['hi'])) - mp.make_mpf(tuple(d['lo'])))


def enum_glen(P, r, a1, max_fans=2000000, max_nodes=400000, max_letters=40000000,
              scan_cap=None, min_width=mp.mpf(10) ** -30, budget=None,
              resume=None, collect_pending=False):
    """Interval enumeration of the WHOLE transversal s in (0,1), with the developed x.

    Returns (cells, stats).  Each cell carries lo, hi, w, `fwd` (= T_c/chi_c, the forward
    developed displacement in units of cot(alpha)), `closure` (= 2n+1), `nopp` (how
    many fans took Lemma B's OTHER branch -- 0 exactly on the leaf class) and, s383,
    the WINDING SUMMARY `gb`/`hmax`/`e_end` (see THE WINDING RIDE-ALONG below).

    A node is (lo, hi, c, x, phi, eps, j0, n, nfan, nopp, v, gmax, hmax, wok) with
    `v`/`gmax`/`hmax`/`wok` the winding passengers (s383) and `phi` the EXACT INTEGER
    index of the fan angle in units of th = pi/(2Q), so head-on is the congruence
    (phi + eps*j*P) = 0 (mod 2Q) and not a float tolerance (as s312_leaf_enum).

    `resume`: a list of PACKED nodes (see `pack_node`) to seed the stack with,
    instead of the root node -- so a capped run continues where it stopped and
    never redoes resolved work.  The returned cells are then only the ones found
    in THIS pass; they are disjoint from any earlier pass by construction, so the
    caller concatenates (⚠ that is not the (closure, fwd) merge the NO MERGE STEP
    note below forbids -- these are disjoint intervals, not duplicate readings).

    `collect_pending`: return the nodes that were DRAINED rather than only a count,
    under `pending` (cap/budget) and `noexit_nodes` ([WFLOOR-113]'s dead end), with
    their total widths under `wpending` / `wnoexit`.  Those two widths plus
    `wsliver` and the cell widths account for the whole transversal -- the
    partition identity a PARTIAL record can be checked against.
    """
    B = base(P, r, a1)
    Q, th, cot, csc = B['Q'], B['th'], B['cot'], B['csc']
    scan_cap = scan_cap or 64 * (a1 + 2)
    rho = lambda i: cot if i % 2 == 0 else csc
    two_q = 2 * Q

    if resume:
        stack = [unpack_node(d) for d in resume]
    else:
        stack = [(mp.mpf(0), mp.mpf(1), mp.mpf(0), mp.mpf(0), 0, 1, 2, 0, 0, 0,
                  0, 0, 0, True)]
    cells = []
    # ⇒⇒ THE WINDING RIDE-ALONG (s383, [OPS-114]'s "the same argument applies to
    # `fill_fanword`").  `gb`/`hmax`/`e_end` are three statistics of the SAME per-fan
    # walk this loop already performs: at each fan the winding advances by
    # `eps*(jR//2)` -- `s357_fanword.trace_fanword` is this node loop at a POINT and
    # accumulates exactly that.  So the store needed a fourth writer only because the
    # loop computed the step and dropped it.  Carried as three INTEGERS in the node
    # (`v` the walk position, `gmax` the running max |v| over every position visited,
    # `hmax` the running max |step|), which is O(1) per fan -- the pointwise tracer's
    # `vs`/`steps` LISTS are what made it cost ~100 B/fan out to nfan ~ 1e6.
    # ⚠ These are PASSENGERS: no split, exit, head-on or branch decision reads them,
    # so the partition is bit-identical with or without them (gated in
    # `probes/s383_fanword_source.py`, which also scores the values themselves).
    # ⇒⇒ THE GRAZE MAP ([OPS-113]/sec.s380).  Every cut this enumerator makes is a ray
    # grazing a VERTEX -- that is what a cell boundary IS (`foundations.md` sec.1) -- and
    # the loop already knows which vertex, because it chose the radius: `rho(i)` is
    # `cot` at an R copy and `csc` at an A copy, and the two centre-splits (`c`, `cp`)
    # are O copies.  It was computing the pair and dropping it, which is why the store
    # needed a separate third-engine `fill_graze` pass at all.  Keyed by the split
    # HEIGHT so no node-tuple changes and `pack_node`/`recap` bit-identity is untouched.
    # ⚠ A RESUMED run only sees ITS OWN cuts; boundaries made before the resume are not
    # in the map (the caller gets None there, never a wrong vertex).
    enum_glen.last_grazes = grazes = {}

    def _graze(y, V, k):
        """Record the graze at height `y`, SHALLOWEST WINS.

        ⚠⚠ The stack is LIFO, so write order is not orbit order, and the same height
        is reached again by deeper copies -- last-write-wins flipped `(5,31)`'s `b2`
        from its true `('O', 12)` to a `('R', 24)` seen later.  The canonical graze is
        the FIRST one along the orbit ([OPS-097]), i.e. the SMALLEST `k`."""
        y = float(y)
        old = grazes.get(y)
        if old is None or k < old[1]:
            grazes[y] = (V, k)
    nodes = capped = noexit = maxscan = nsliver = 0
    nfalse = nheadon_odd = 0          # s341: false head-ons rejected / parity violations
    wsliver = mp.mpf(0)
    pending, noexit_nodes = [], []
    wpending = wnoexit = mp.mpf(0)
    stop = (time.time() + budget) if budget else None
    while stack:
        nd = stack.pop()
        lo, hi, c, x, phi, eps, j0, n, nfan, nopp, v, gmax, hmx, wok = nd
        nodes += 1
        if nodes > max_nodes or nfan > max_fans or n > max_letters:
            capped += 1
            wpending += hi - lo
            if collect_pending:
                pending.append(pack_node(nd))
            continue
        if stop and nodes % 4096 == 0 and time.time() > stop:
            max_nodes = 0                      # drain the stack into `capped`
            capped += 1
            wpending += hi - lo
            if collect_pending:
                pending.append(pack_node(nd))
            continue
        # GRAZE RAYS.  At a cell boundary the ray grazes a vertex and the chain need not
        # close at all ([OPS-050]); un-pruned, such a sliver runs forever -- one at
        # s = 1 of width 3.5e-49 was this build's first failure.  Drop below a floor far
        # under any real cell and REPORT the total width dropped, so the Sum w = 1 guard
        # still sees any genuine loss ([NCYL-020]: a truncation is not a blank).
        if hi - lo < min_width:
            nsliver += 1
            wsliver += (hi - lo)
            continue
        # --- the exit scan.  Consult vertices ONE AT A TIME and split only at the ones
        # that actually fall inside the interval, so the window grows on demand and the
        # enumerator never cuts at a vertex it did not need ([WFLOOR-113] reinforcement
        # pushes exit indices past s312_leaf_enum's fixed imax).  Below, `c` and every
        # consulted `y` are OUTSIDE (lo,hi), so for every s in (lo,hi) the comparisons
        # `y < s` and `c < s` are constant and read off the endpoints.
        jstar = None
        split = c if lo < c < hi else None
        if split is not None:
            # the fan CENTRE is an O copy; it was established by the crossing that
            # created this node, i.e. at the node's own accumulated depth.
            _graze(split, 'O', n)
        i = j0
        while split is None:
            if i - j0 > scan_cap:
                break
            y = c + rho(i) * mp.sin((phi + eps * i * P) * th)
            if lo < y < hi:
                split = y
                _graze(y, 'R' if i % 2 == 0 else 'A', n + (i - j0) + 1)
                break
            if (y <= lo) if (c <= lo) else (y >= hi):       # vertex on the centre's side
                jstar = i
                maxscan = max(maxscan, i - j0)
                break
            i += 1
        if split is not None:
            # an interval CUT consumes no fan, so the winding passengers ride across
            # unchanged -- the same reason `n`/`nfan`/`nopp` do.
            if split - lo > 0:
                stack.append((lo, split, c, x, phi, eps, j0, n, nfan, nopp,
                              v, gmax, hmx, wok))
            if hi - split > 0:
                stack.append((split, hi, c, x, phi, eps, j0, n, nfan, nopp,
                              v, gmax, hmx, wok))
            continue
        if jstar is None:
            noexit += 1
            wnoexit += hi - lo
            if collect_pending:
                noexit_nodes.append(pack_node(nd))
            continue
        n2 = n + (jstar - j0 + 1)
        ang = phi + eps * jstar * P
        jR = jstar - (jstar % 2)
        Om = phi + eps * jR * P
        # HEAD-ON -- ⚠ AMENDED s341: THE TEST IS AT THE **R VERTEX**, NOT THE EXIT VERTEX
        # ([WFLOOR-115]).  Retrace happens when the ray meets L1 PERPENDICULARLY, and the
        # triangle's right angle sits at R so OR ⊥ L1 -- hence the condition is that
        # O–v_{jR} be HORIZONTAL, sin(Om*th) = 0.  Testing sin at the EXIT vertex instead
        # is the SAME test whenever j* is even, and at odd P it always is, which is why
        # this was invisible on 51/51 rows: phi stays in 2P·Z (mod 2Q), so a head-on needs
        # P(2k + eps*j*) = 0 (mod 2Q), i.e. (2k + eps*j*) = 0 (mod 2Q/gcd(P,2Q)) -- at odd
        # P that modulus is 2Q, EVEN, forcing j* even; at even P (hence Q odd) it is Q,
        # ODD, and admits j* of BOTH parities.  An odd solution puts an A vertex at the
        # centre height WITHOUT making L1 vertical: a false head-on that closes the cell
        # far too early.  `jstar` is returned so the caller can assert the parity.
        if ang % two_q == 0 and jstar % 2:
            nfalse += 1               # the PRE-s341 test would have closed the cell here
        # THE FAN'S WINDING STEP, taken with the eps in force BEFORE the branch below
        # -- the order `trace_fanword` uses (`steps.append`, then `v +=`, then the
        # branch flips eps).  `jR >= 0`, so |step| == jR//2.
        step = eps * (jR // 2)
        v2, g2, h2 = v + step, max(gmax, abs(v + step)), max(hmx, jR // 2)
        if Om % two_q == 0:
            nheadon_odd += jstar % 2  # must stay 0: a head-on exits AT its R vertex
            xR = x + cot * ((-1) ** ((Om // two_q) % 2))
            cells.append(dict(lo=lo, hi=hi, fwd=1 - xR / cot, closure=2 * n2 + 1,
                              nfan=nfan + 1, nopp=nopp, jstar=jstar,
                              gb=(g2 if wok else None),
                              hmax=(h2 if wok else None),
                              e_end=(v2 if wok else None)))
            continue
        cp = 2 * (c + cot * mp.sin(Om * th)) - c
        if lo < cp < hi:                                    # the opp test's breakpoint
            # `cp` IS sec.s380's crossing lemma `y(O') = 2 y(R) - y(O)` -- the reflection
            # of the fan centre in the R vertex, the angle at R being right.  So this
            # boundary's graze vertex is the next copy's O.
            _graze(cp, 'O', n2 + 1)
            stack.append((lo, cp, c, x, phi, eps, j0, n, nfan, nopp,
                          v, gmax, hmx, wok))
            stack.append((cp, hi, c, x, phi, eps, j0, n, nfan, nopp,
                          v, gmax, hmx, wok))
            continue
        xp = 2 * (x + cot * mp.cos(Om * th)) - x
        yj = c + rho(jstar) * mp.sin(ang * th)
        s = (lo + hi) / 2
        opp = (cp - s) * (yj - s) < 0
        if opp:                                             # Lemma B
            stack.append((lo, hi, cp, xp, Om + two_q, -eps, 1 + (jstar % 2),
                          n2 + 1, nfan + 1, nopp, v2, g2, h2, wok))
        else:                                               # [WFLOOR-113], eps UNCHANGED
            stack.append((lo, hi, cp, xp, Om + two_q, eps, 2 - (jstar % 2),
                          n2 + 1, nfan + 1, nopp + 1, v2, g2, h2, wok))
    # NO MERGE STEP.  `s312_leaf_enum` needs one because it cuts at EVERY vertex of the
    # window, including ones past the exit, so one cell arrives as several nodes.  The
    # incremental scan above consults only vertices it needs, and every cut it makes is a
    # genuine decision boundary -- so merging would be wrong, not merely redundant: at
    # 3/7 the two J-PAIRED cells are adjacent with identical closure AND identical
    # forward displacement, and a (closure, fwd) merge silently fuses them (n = 2 for a
    # 3-cell row, with S still exactly right -- a count error invisible to the S guard).
    cells.sort(key=lambda z: z['lo'])
    graze = [z for z in cells if z['hi'] - z['lo'] < TINY]
    merged = [dict(z) for z in cells if z['hi'] - z['lo'] >= TINY]
    for z in merged:
        z['w'] = z['hi'] - z['lo']
    W = sum((z['w'] for z in merged), mp.mpf(0))
    S = sum((z['w'] * z['fwd'] for z in merged), mp.mpf(0))
    return merged, dict(nodes=nodes, capped=capped, noexit=noexit, graze=len(graze),
                        sliver=nsliver, wsliver=wsliver,
                        nfalse=nfalse, nheadon_odd=nheadon_odd,
                        maxscan=maxscan, raw=len(cells), n=len(merged), sumw=W, S=S,
                        g_len=S / Q, ok=(capped == 0 and noexit == 0),
                        maxclos=max((z['closure'] for z in merged), default=0),
                        nopp=sum(z['nopp'] for z in merged),
                        # cells whose winding history was lost to a PRE-s383 resume
                        # (`wok=False`).  A loss must COUNT, not read as a blank.
                        nowind=sum(1 for z in merged if z.get('gb') is None),
                        pending=pending, noexit_nodes=noexit_nodes,
                        wpending=wpending, wnoexit=wnoexit)


# ------------------------------------------------------------------ G1/G2/G3: the ring

def check_row(P, Q, max_depth=200000, verbose=True):
    """One row against the ring: cell count, boundaries, per-cell T_c/chi_c, and g_len."""
    from s275_flow_exact import cells_exact, swept_sum
    a1, r = divmod(Q, P)
    ctx, rows, guard = cells_exact(P, Q, max_depth=max_depth)
    if guard:
        if verbose:
            print(f'  {P}/{Q}: ring SKIPPED ({guard})')
        return None
    Sr = swept_sum(ctx, rows).to_float()
    t0 = time.time()
    cells, st = enum_glen(P, r, a1)
    secs = time.time() - t0
    # pair by MIDPOINT, never by index (a missed cell would shift every later pairing)
    worst_b = worst_t = 0.0
    matched = set()
    for z in cells:
        mid = float(z['lo'] + z['hi']) / 2
        j = min(range(len(rows)), key=lambda i: abs((rows[i]['lo'] + rows[i]['hi']) / 2 - mid))
        matched.add(j)
        d = rows[j]
        worst_b = max(worst_b, abs(float(z['lo']) - d['lo']), abs(float(z['hi']) - d['hi']))
        tc = d['T'].to_float() / d['chi']
        worst_t = max(worst_t, abs(float(z['fwd']) - tc) / max(1.0, abs(tc)))
    rec = dict(P=P, Q=Q, a1=a1, r=r, n_ring=len(rows), n_model=st['n'],
               n_matched=len(matched), sumw=float(st['sumw']),
               g_ring=Sr / Q, g_model=float(st['g_len']),
               dg=abs(Sr / Q - float(st['g_len'])), max_bnd_err=worst_b,
               max_T_err=worst_t, nopp=st['nopp'], nodes=st['nodes'],
               capped=st['capped'], noexit=st['noexit'], maxscan=st['maxscan'],
               ring_depth=max(d['depth'] for d in rows), model_clos=st['maxclos'],
               secs=round(secs, 2))
    rec['ok'] = (rec['n_ring'] == rec['n_model'] == rec['n_matched']
                 and abs(rec['sumw'] - 1) < 1e-12 and rec['dg'] < 1e-11
                 and rec['max_T_err'] < 1e-11 and st['ok'])
    if verbose:
        ok = rec['ok']
        print(f"  {P:2d}/{Q:<4d} a1={a1} r={r} | ring n={rec['n_ring']:3d} "
              f"model n={rec['n_model']:3d} match={rec['n_matched']:3d} "
              f"Sw-1={rec['sumw']-1:+.1e} bnd={worst_b:.1e} T={worst_t:.1e} "
              f"g={rec['g_model']:.9f} dg={rec['dg']:.1e} nopp={rec['nopp']:4d} "
              f"{secs:6.2f}s  {'OK' if ok else '** FAIL **'}", flush=True)
        rec['ok'] = ok
    return rec


def mode_check(box=None):
    box = box or [(P, Q) for Q in range(5, 28) for P in range(3, Q, 2)
                  if gcd(P, Q) == 1 and Q // P >= 2]
    print(f'=== G1/G2/G3 — enumerator vs ring, {len(box)} candidate rows (odd P) ===')
    rows = [r for r in (check_row(P, Q) for (P, Q) in box) if r]
    good = [r for r in rows if r['ok']]
    print(f'\n  {len(good)}/{len(rows)} certified rows reproduced exactly; '
          f'{len(box) - len(rows)} rows the RING could not certify')
    if rows:
        print(f"  worst |Sum w - 1| = {max(abs(r['sumw']-1) for r in rows):.1e}, "
              f"worst boundary {max(r['max_bnd_err'] for r in rows):.1e}, "
              f"worst T_c/chi_c {max(r['max_T_err'] for r in rows):.1e}, "
              f"worst |dg_len| = {max(r['dg'] for r in rows):.1e}")
    bad = [r for r in rows if not r['ok']]
    if bad:
        print('  FAILING ROWS: ' + ', '.join(f"{r['P']}/{r['Q']}" for r in bad))
    return rows


def mode_scope():
    """G5: the head-on congruence is degenerate iff 2Q/gcd(P,2Q) < the exit window."""
    print('=== G5 — head-on degeneracy: is the criterion `2Q/gcd(P,2Q) < imax`, not `P even`? ===')
    rows = []
    for (P, Q) in [(2, 7), (4, 15), (2, 9), (4, 17), (6, 19), (6, 23), (8, 17),
                   (3, 7), (5, 16), (3, 13), (5, 21), (7, 16), (9, 20)]:
        a1, r = divmod(Q, P)
        j0 = headon_index(P, Q)
        pred_bad = j0 < 4 * a1 + 6
        rec = check_row(P, Q, verbose=False)
        if rec is None:
            print(f'  {P:2d}/{Q:<4d} predicted {"DEGENERATE" if pred_bad else "clean":11s} '
                  f'(j={j0} vs imax={4*a1+6})  ring SKIPPED')
            continue
        rec.update(headon_j=j0, imax=4 * a1 + 6, predicted_degenerate=pred_bad)
        rows.append(rec)
        print(f'  {P:2d}/{Q:<4d} P{"even" if P % 2 == 0 else " odd"} '
              f'predicted {"DEGENERATE" if pred_bad else "clean":11s} '
              f'(j={j0} vs imax={4*a1+6})  ->  observed '
              f'{"BROKEN" if not rec["ok"] else "exact"}   '
              f'{"PREDICTION HELD" if pred_bad != rec["ok"] else "** PREDICTION FAILED **"}')
    held = sum(1 for r in rows if r['predicted_degenerate'] != r['ok'])
    print(f'  {held}/{len(rows)} rows match the criterion')
    return rows


def mode_ladder(P, r, ms):
    """G4: walk a ladder past the ring's depth wall."""
    print(f'=== G4 — ladder (P,r)=({P},{r}), Q = {P}m+{r} ===')
    rows = []
    for m in ms:
        Q = P * m + r
        if gcd(P, Q) != 1:
            continue
        t0 = time.time()
        cells, st = enum_glen(P, r, m)
        secs = time.time() - t0
        rec = dict(P=P, r=r, m=m, Q=Q, n=st['n'], sumw=float(st['sumw']),
                   S=float(st['S']), g_len=float(st['g_len']), nopp=st['nopp'],
                   maxclos=st['maxclos'], nodes=st['nodes'], ok=st['ok'],
                   capped=st['capped'], noexit=st['noexit'], secs=round(secs, 2))
        rows.append(rec)
        print(f"  m={m:4d} Q={Q:5d} n={rec['n']:4d} Sw-1={rec['sumw']-1:+.1e} "
              f"S={rec['S']:12.6f} g_len={rec['g_len']:.9f} S-n={rec['S']-rec['n']:9.5f} "
              f"maxclos={rec['maxclos']:9d} nopp={rec['nopp']:6d} "
              f"{'' if rec['ok'] else 'CAPPED '}{secs:7.2f}s", flush=True)
    return rows


LADDERS = [(3, 1), (3, 2), (5, 2), (5, 3), (5, 4), (7, 3), (7, 4)]


def mode_reach(budget=45.0, ladders=None):
    """G4: how far each ladder goes under a per-row wall-clock budget.

    Stops on TWO CONSECUTIVE budget failures, never one -- cost along a ladder is
    ERRATIC, not monotone in m (`s338_ladder_limit` trap (i), and it recurs here:
    (3,1) misses the budget at m=10 and then certifies m=11 in 43.6 s).
    """
    out = {}
    for (P, r) in (ladders or LADDERS):
        print(f'--- ladder ({P},{r}), budget {budget:.0f}s/row ---', flush=True)
        rows, fails = [], 0
        for m in range(2, 200):
            Q = P * m + r
            if gcd(P, Q) != 1:
                continue
            t0 = time.time()
            cells, st = enum_glen(P, r, m, budget=budget)
            dt = time.time() - t0
            if not st['ok']:
                fails += 1
                print(f'    m={m:3d} Q={Q:5d}  BUDGET ({dt:.0f}s)', flush=True)
                if fails >= 2:
                    break
                continue
            fails = 0
            rec = dict(m=m, Q=Q, n=st['n'], S=float(st['S']), g_len=float(st['g_len']),
                       sumw_err=float(st['sumw']) - 1, maxclos=st['maxclos'],
                       nopp=st['nopp'], nodes=st['nodes'], secs=round(dt, 2))
            rows.append(rec)
            print(f"    m={m:3d} Q={Q:5d} n={st['n']:5d} Sw-1={rec['sumw_err']:+.0e} "
                  f"S={rec['S']:11.5f} g_len={rec['g_len']:.9f} "
                  f"S-n={rec['S']-rec['n']:8.4f} maxclos={rec['maxclos']:10d} "
                  f"{dt:6.2f}s", flush=True)
        out[f'{P}_{r}'] = rows
    return out


def mode_dps(rows=((3, 1, 11), (5, 2, 14), (7, 4, 11), (5, 3, 12), (3, 1, 9))):
    """The float self-check.  This instrument is mpmath, NOT the ring: past the ring's
    reach there is no certificate, so the only available self-check is precision
    stability (TOOLS.md's FAN MODEL trap (ii)).  A row whose n or g_len moves between
    dps 50 and dps 80 is not to be quoted."""
    res, out = {}, []
    for dps in (50, 80):
        mp.mp.dps = dps
        for (P, r, m) in rows:
            cells, st = enum_glen(P, r, m)
            res[(dps, P, r, m)] = (st['n'], float(st['g_len']))
    mp.mp.dps = 50
    for (P, r, m) in rows:
        a, b = res[(50, P, r, m)], res[(80, P, r, m)]
        rec = dict(P=P, r=r, m=m, Q=P * m + r, n50=a[0], n80=b[0],
                   dg=abs(a[1] - b[1]), stable=(a[0] == b[0] and abs(a[1] - b[1]) < 1e-12))
        out.append(rec)
        print(f"  ({P},{r}) m={m:3d} Q={rec['Q']:4d}: n {a[0]} vs {b[0]}, "
              f"|dg_len| = {rec['dg']:.2e}  {'stable' if rec['stable'] else '** MOVED **'}")
    return out


def main():
    """⚠ ONE FILE PER MODE, deliberately ([OPS-013]).  A single shared store is read at
    process start and written at process end, so a short mode run silently CLOBBERS a
    long one that started earlier -- this probe's `scope` run did exactly that to its
    own `reach` and `dps` blocks while both were on disk."""
    mode = sys.argv[1] if len(sys.argv) > 1 else 'check'
    os.makedirs('data', exist_ok=True)
    if mode == 'check':
        out, res = OUT, mode_check()
    elif mode == 'scope':
        out, res = 'data/s340_scope.json', mode_scope()
    elif mode == 'reach':
        out = 'data/s340_reach.json'
        res = mode_reach(float(sys.argv[2]) if len(sys.argv) > 2 else 45.0)
    elif mode == 'dps':
        out, res = 'data/s340_dps.json', mode_dps()
    elif mode == 'ladder':
        P, r = int(sys.argv[2]), int(sys.argv[3])
        ms = [int(v) for v in sys.argv[4:]] or list(range(2, 12))
        out, res = f'data/s340_ladder_{P}_{r}.json', mode_ladder(P, r, ms)
    else:
        raise SystemExit(f'unknown mode {mode}')
    json.dump(res, open(out, 'w'), indent=1, default=str)
    print(f'wrote {out}')


if __name__ == '__main__':
    main()
