#!/usr/bin/env python3
"""
s439_exact_cells.py -- PRE-REGISTERED: enumerate a strip's cells of constant word by
GENERATING their boundaries from the vertices, instead of SEARCHING for them by
sampling + bisection.  User-directed (s439): *"I don't like the sampling approach ... If
there is a better approach here we use it."*

THE IDENTITY THE INSTRUMENT IS BUILT ON, and it is already the repo's own.
`foundations.md` §1: **a cell boundary IS a vertex graze** -- `x` on side `σ` is a cell
boundary of the strip `(σ, u)` exactly when the billiard orbit through `x` in direction
`u·π/(2Q)` passes through a vertex.  Every cylinder census in this repo
(`s315_cylinder_count`, `s318_kac_missing`, `s437_oblique_cylinders`) reads that
backwards -- it samples `NSAMP` launch points, watches the closure word change, and
bisects.  ⇒⇒ **RUN IT FORWARDS: launch FROM each vertex along each grid direction of the
class, and every bounce of that orbit is, by the same identity, a cell boundary of the
strip it lands in.**  The boundary set is then CONSTRUCTED, and there is no resolution
parameter anywhere in the instrument.

WHY THIS IS NOT A SPEED OPTIMISATION.  `probes/s439_cylinder_census.py`'s own design note
records the reason: at `NSAMP = 400` the `3/10` odd class reads `C = 6` with
`total_area_ratio = 1.0` to 10 dp (`data/s437_convergence.json`), so the Kac completeness
certificate is BLIND to that row's under-resolution and only a second resolution catches
it.  A sampled enumerator therefore cannot certify its own output -- it can only agree
with itself at two grid sizes.  Here `Σ width = |σ|` is again automatic and again proves
nothing ([OPS-200]), but the boundary set is exhaustive BY CONSTRUCTION rather than by
convergence, so the failure mode it guards is absent rather than untested.  The remaining
error is float drift along one trace, which is bounded and MEASURED below (H3), not a
sampling gap.

  ⚠ THE ONE PLACE THIS CAN STILL BE INCOMPLETE, stated up front.  A prong that does not
  TERMINATE at a vertex within the cap contributes only a truncated boundary list.  In a
  completely periodic direction every prong terminates ([NCYL-001]); in a direction with
  a minimal component some do not, and the detector is the same one as always -- Kac area
  short of the surface.  Every prong's termination reason is stored, and a row with any
  `cap` prong is reported `PRONG_CAPPED` and is not scored.  This is [NCYL-020]'s "a
  truncated enumeration is not a missing orbit" for this instrument.

WHY THE THREE VERTICES AND NOT TWO.  `R` is REGULAR on the surface ([NCYL-001] G3) and is
never a separatrix ARRIVAL, which is why the CP census (`s276_separatrix_census`,
`s438_oddclass_cp`) launches from `O` and `A` only.  But regular-on-the-surface is not
irrelevant-to-the-WORD: a billiard orbit that meets the right angle exactly reflects off
two perpendicular walls and RETRACES, and the orbits either side of it hit those two
walls in opposite orders -- so an `R` graze splits a cell of constant WORD while leaving
the CYLINDER intact.  Omitting `R` would give the cylinder-level partition directly, but
then a cell's raw word would not be constant across it and the comparison against
`s437_oblique_cylinders` would no longer be like-for-like.  ⇒ all three vertices are
launched; the `R`-split cells are re-merged downstream by `s315.merge_duplicate_cylinders`
exactly as they are today (its `≤3` word mismatches under best cyclic rotation is
calibrated for precisely this, s315 trap (iii)).

  ⚠ AND `R`'s PRONGS TERMINATE DIFFERENTLY.  A geodesic through a REGULAR point in a
  completely periodic direction is an ordinary closed leaf: it need never reach a cone
  point, so an `R` prong terminates by RETURNING TO `R`, not by arriving at `O`/`A`.  The
  tracer therefore stops at any of the three vertices, not at the two cone points.

THE INWARD WINDOWS ARE EXACT INTEGER INTERVALS, which removes the last float predicate
from the enumeration side.  In units of `π/(2Q)`: at `O` the interior cone spans
`u ∈ (0, P)`, at `A` `u ∈ (2Q+P, 3Q)`, at `R` `u ∈ (Q, 2Q)` -- open.  The first two are
`s438_oddclass_cp.inward_count_exact` verbatim, which is a free cross-instrument check
(H2); the third is new and is the `R` window.

HYPOTHESES (pre-registered).
  H1  CONTROL, RUN FIRST, ABLE TO FAIL, AND IT IS THE WHOLE POINT.  On s437's 8 centres
      this module reproduces `s437_oblique_cylinders.run_row` EXACTLY: per class the same
      `C_total` and `total_area_ratio` to `1e-9`, and per strip the same `n_cells`.  The
      two instruments share `assemble`/`canon`/`merge` and the tracer, and share NO
      boundary-finding code whatever -- one bisects a sampled sign change, the other
      launches from a vertex.  A disagreement is real and says which one is wrong.
  H2  The `O`/`A` prong counts equal `s438_oddclass_cp.inward_count_exact`, per centre,
      per parity.  Shares no code; a mismatch means the window derivation is wrong and
      the enumeration is not exhaustive -- a STOP condition, not a warning.
  H3  THE ERROR BUDGET, which is what replaces the resolution ladder.  Report per strip
      the minimum gap between distinct boundaries and the maximum disagreement against
      s437's bisected boundaries.  ⚠ NO tolerance is pre-registered as acceptable; the
      number is the measurement, and if boundary drift over a long prong ever approaches
      the minimum cell gap this instrument needs higher precision, not a bigger cap.
  H4  THE PAYOFF, and it is a prediction that can fail: cost is
      `O(#prongs + #cells)` traces per centre against the sampled enumerator's
      `O(#strips × NSAMP)`, i.e. a speedup of order `NSAMP`.  PREDICT: two to three
      orders of magnitude on s437's rows.  NO prediction on `Q`-reach -- prong LENGTH
      still grows with `Q` and that is not addressed here.

  ⚠ NOT pre-registered as expected-to-fail: [OPS-044] bars the self-applied device.  The
  scoring discipline is [OPS-041]: H1 and H2 can both come out otherwise, and `Σ width =
  |σ|` cannot, so it is reported and not counted ([OPS-200]).

PRIOR ART: grepped 'cell boundary', 'vertex graze', 'graze', 'separatrix', 'saddle
connection', 'exact boundar', 'interval propagation', 'enum_glen', 'enum_leaf',
'NSAMP', 'bisect' across rulings.md / TOOLS.md / foundations.md / cf_width_laws.md /
computational_findings.md / probes/ / engine/ ->
  - `foundations.md` §1 OWNS "a cell boundary is a zero of a vertex-graze function
    `η_{k,V}`"; this module is that statement used constructively.  It coins nothing.
  - `s340_glen_enum.enum_glen` (+ `enum_glen.last_grazes`, TOOLS.md) IS the repo's
    existing non-sampling cell enumerator and IS `cell_store`'s source -- exact interval
    propagation, no grid, no depth cap, and it already reports the boundary GRAZE PAIR
    `(vertex, k)`.  ⚠ It is scoped to the PERPENDICULAR beam and is a bespoke model of
    that beam's fan structure (CF parameters `(P, r, a₁)`, Lemma-B's other branch,
    [WFLOOR-111]/[WFLOOR-113]); it has no notion of a launch SIDE or a launch ANGLE, so
    it does not extend to the `4Q` oblique strips on three sides.  What transfers is its
    PRINCIPLE -- boundaries are grazes, so generate them -- which is what this does.
  - `s312_leaf_enum` is the same technique scoped to the leaf class ([OPS-050]).
  - `s186_boundaries`/`s187_boundaries2`/`s187_separatrix`/`s188_stripword` extract
    saddle connections EXACTLY, but from the `3/7` translation ATLAS of the parked Veech
    thread; they are one-centre tools and are not general enumerators.
  - `s276_separatrix_census`/`s438_oddclass_cp` trace prongs from `O`/`A` and ask only
    whether each ARRIVES; they discard the bounce sequence, which is exactly the data
    this module keeps.  The `O`/`A` window arithmetic is theirs (H2).
  - `s315_cylinder_count`/`s437_oblique_cylinders` OWN `assemble`/`canon`/
    `merge_duplicate_cylinders`/the Kac area and all four traps; reused verbatim so that
    H1 isolates the boundary finder as the only difference.
  - NO ruling anywhere states a general vertex-launched cell enumerator; the closest is
    [OPS-050]'s "no fast float route to `g_len` exists", which is about `g_len` on the
    perpendicular beam and is not contradicted here.

================================================================================
s440 -- QUEUE ITEM (i): `mode="claim"`, THE COUNT READER.  User-directed at s439 as a
BUILD, not a question.  The instrument above is correct and fast at small `Q` and simply
does not FINISH at large `Q`: `520` class-rows and `293` centres of the s439 sweeps came
back `PRONG_CAPPED` or `BUDGET`, and the drop is biased -- the LONE class ran `100%`
scored at `Q <= 13` and `19%` at `Q = 41-53` ([NCYL-259]), so `C(lone) = floor(Q/2)`, the
half of the law that survives ([NCYL-253]), is the half being selected on.

WHAT [OPS-205] ESTABLISHED, and it is the whole design brief.  The prong CAP was never
the binding cost -- at 7/20 all 19 prongs terminate at `1e6` in 4.8 s where 4 cap at
`2e5`.  What the longer prong buys you is `2e5` CELLS in a single strip, and the readers
above walk every one of them, canonicalising and hashing an `O(period)` word each time.
THAT is the hour.

THE TWO CHANGES.
  (1) A TWO-STAGE PRONG CAP: everything at `CAP`, then only the prongs that did not
      terminate re-traced at `CAP_BIG = 2e6`, and the class ABANDONED at the first prong
      still open there (it is `PRONG_CAPPED`, hence unscorable, hence pure waste).
  (2) `class_claim`: trace ONE cell per cylinder to closure and let the orbit CLAIM every
      cell it lands in, in every strip at once.  A cylinder's cells in a strip are one
      first-return cycle, so a single closed orbit visits all of them; the trace count is
      the CYLINDER count, not the cell count.  `k_c` and the met-sets come off the same
      trace (`visits`), so they stay exact.  The stop is [NCYL-247]'s certificate
      `Σ_c circ_c·h_c = Q·cot α`, not the strip list ([OPS-204]).

  ⚠ THIS IS NOT THE WIDTH-GROUPING READER [OPS-205] SKETCHED, and the difference is a
  tolerance.  Grouping a strip's cells by width first needs a threshold that [OPS-205]
  itself measured to be a PLATEAU (`1e-14 -> 6655` groups, `1e-9 -> 6`, the answer `8`
  only over `1e-12..1e-11`), re-established per row.  Claiming needs none: the orbit says
  which cells are its own.  Width enters only as a REJECTION test on one claim.

  ⚠⚠ THREE THINGS MEASURED HERE THAT CONTRADICT THE OBVIOUS DESIGN, each of which cost a
  run and each of which is a DO-NOT (details at the code):
    (a) THE WIDTH GUARD AT `1e-6` IS WRONG.  [NCYL-247]'s own figure for `width = h_c/flux`
        is a max relative deviation of `3.8e-05`; at `1e-6` the guard refused 445280 of
        ONE orbit's own crossings and the reader re-traced them as new cylinders.  `1e-3`.
    (b) AVERAGING THE CROSSING WIDTHS MAKES THE HEIGHT WORSE, NOT BETTER, by 3000x --
        `k_c ~ 10^4` crossings look like independent estimates whose mean should be
        tighter, but the deviations are ONE-SIDED (a crossing cell can be cut short by a
        further graze and never runs long).  Use the LAUNCH cell.
    (c) LAUNCH FROM A PERPENDICULAR STRIP; ORDERING BY COST INSTEAD BREAKS THE ROW.  The
        launch cell calibrates `h_c`, so it must be a FULL crossing, and on a strip cut
        finer than the cylinder's height it is a fragment.  7/20 lone: 10 traces and a
        `5e-10` certificate from `H:67`, against 12 word-families and `area_ratio_hc =
        0.73` from the cheapest strip.
  ⚠ And the word primitives had to be rewritten: `s315.canon` and
  `best_rotation_mismatch` are both QUADRATIC in the word length, which is `6.6e5` here.

  ⚠⚠ WHAT CLAIM MODE GIVES UP, stated up front: the per-CELL Kac integrand.  Its
  `total_area_ratio` is a +-10% ESTIMATE (`seg0` sampled per (cylinder, strip), not
  integrated per cell) -- kept only because it is INDEPENDENT of the certificate and so
  catches a gross error, never as a verdict.  The verdict column is `area_ratio_hc`.

CONTROL (pre-registered at s439, in `s439_exact_census.control_claim`): claim mode must
reproduce the FULL-strip reader on every centre already banked in
`data/s439_exact_store` -- `C_total`, `C_met_sigma1` AND the `k_c` multiset, per class.
The two readers share the boundary generator, `assemble`, `canon` and the tracer and
differ in exactly one thing, so a disagreement is a claim error and nothing else.
================================================================================

Run: PYTHONPATH=.:engine:archive/scripts_2026-07:archive/scripts_2026-08:probes \
     .venv/bin/python3.13 probes/s439_exact_cells.py control
"""
import bisect
import hashlib
import json
import math
import random
import sys
import time
from array import array
from collections import defaultdict

from right_triangle_billiards import RightTriangleBilliard
import s315_cylinder_count as s315
import s437_oblique_cylinders as s437

POS_TOL = 1e-7          # closure tolerance, inherited from s437/s315.  ⚠ A CEILING ONLY
                        # since s443 -- see `SCALED_CLOSURE` and `trace_closure_claim`.
SCALED_CLOSURE = True   # ⇒⇒ s443 / [OPS-211]: SCALE THE CLOSURE TEST WITH THE CYLINDER.
                        # A flat `POS_TOL` falsely closes any cylinder whose whole
                        # cross-section fits inside it (`h < 1e-7`), at HALF the true
                        # period -- 19 cylinders of `data/s440_claim_store`, first at
                        # Q = 37, tell = an ODD word length.  Set False to reproduce the
                        # pre-s443 instrument (the s443 A/B control does exactly that);
                        # it is not a tuning knob and there is no reason to flip it in
                        # new work.
VERT_TOL = 1e-9         # prong termination: distance to a vertex
DEDUP_TOL = 1e-10       # two boundary positions closer than this are ONE boundary
MIN_WIDTH = 1e-10
CAP = 200_000           # stage-1 prong cap / the s439 setting (see CAP_BIG)
CAP_BIG = 2_000_000     # ⇒ s440, [OPS-205]: stage-2 prong cap and the CLOSURE cap in claim
                        # mode.  Measured free: at 7/20's lone class all 19 prongs terminate
                        # at 1e6 (longest 663139 bounces, 4.8 s) where 4 cap at 2e5.
WIDTH_REL_TOL = 1e-3    # claim mode: a claimed cell's width against `h_c/flux`
                        # ([NCYL-247]).  ⚠ CALIBRATED, and 1e-6 IS WRONG -- [NCYL-247]'s
                        # own figure for that identity is a max relative deviation of
                        # 3.77e-05 (median 2.9e-06), i.e. the PIPELINE FLOOR, and at 7/20
                        # a cell of width 6e-06 cut by an 8.9e5-bounce prong carries
                        # exactly that.  At 1e-6 the guard refused 445280 of one orbit's
                        # own crossings and the reader re-traced them as new cylinders.
                        # The realised deviation is reported per class (`width_dev_max`)
                        # so the margin is a measurement, not an assumption.
TRACE_CAP = 200         # claim mode: closure traces per class before the row is abandoned.
                        # ⚠ the real guard is the census per-centre BUDGET; this only
                        # bounds the waste on a row whose certificate never closes.
CERT_TOL = 1e-3         # claim mode: the [NCYL-247] completeness certificate, used BOTH
                        # as the stop and (in the census) as the verdict, so the reader
                        # never abandons a row its own verdict would then have passed.
                        # ⚠ CALIBRATED AGAINST THE OBJECT IT COULD SKIP: the smallest
                        # cylinder area RATIO anywhere in `data/s439_exact_store` (1849
                        # cylinders, Q <= 29) is 6.3e-04, five orders of magnitude above
                        # this -- BUT the realised residual GROWS WITH Q as the prong
                        # drift does (4.7e-10 at 7/20, 9.1e-07 by Q = 27, 2.3e-06 at
                        # 13/24), so a tolerance tight enough to be a real completeness
                        # test today is one that stops passing tomorrow, and at 13/24 a
                        # 1e-06 gate did exactly that: the class was complete at 12
                        # traces, read `1.0000023`, and the reader went hunting.
                        # ⇒⇒ SO THE AREA IS NOT THE PRIMARY TEST ANY MORE.  The primary
                        # test is EXACT and integer -- every cell of the class's
                        # PERPENDICULAR strips is claimed, i.e. every crossing of the one
                        # transversal that cannot be fragmented ([NCYL-248]: `k ∈ {1,2}`
                        # there) is attributed to a traced cylinder.  This tolerance is
                        # now the GROSS-error gate beside it, and it is reported per row
                        # (`area_ratio_hc`) beside the smallest cylinder's own share
                        # (`min_area_ratio`) so the margin stays visible.


# ---------------------------------------------------------------- geometry

def vertices(P, Q):
    L = 1.0 / math.tan((P / Q) * (math.pi / 2))
    return {"O": (0.0, 0.0), "R": (L, 0.0), "A": (L, 1.0)}


def inward_window(vname, P, Q):
    """The OPEN integer window of grid directions pointing into the triangle at `vname`,
    in units of `π/(2Q)`.  `O` and `A` are `s438_oddclass_cp.inward_count_exact`'s (H2);
    `R` is the right angle, spanning direction `π/2` (towards `A`) to `π` (towards `O`)."""
    return {"O": (0, P), "A": (2 * Q + P, 3 * Q), "R": (Q, 2 * Q)}[vname]


def prongs(vname, P, Q, parity):
    a, b = inward_window(vname, P, Q)
    return [u for u in range(a + 1, b) if u % 2 == parity]


def _u_of(vx, vy, Q):
    """Grid index of a direction, or None if it has left the `4Q` grid."""
    uf = (math.atan2(vy, vx) % (2.0 * math.pi)) * 2.0 * Q / math.pi
    u = round(uf)
    if abs(uf - u) > 1e-6 * 2 * Q:
        return None
    return int(u) % (4 * Q)


# ------------------------------------------------------- word ops at 10^5 letters
# ⇒⇒ s440.  A cylinder word here reaches 6.6e5 letters ([OPS-205]), and BOTH of s315's
# word primitives are QUADRATIC in that length: `canon` builds all `n` rotations
# (`min(word[i:]+word[:i] ...)`), `best_rotation_mismatch` compares all `n` of them.  At
# `n = 6.6e5` that is 4e11 character operations -- not slow, IMPOSSIBLE.  Both are
# replaced below by O(n) / output-equivalent forms, and both replacements are CONTROLLED
# against s315's originals (`word_ops_control`), which stay the reference.
#   ⚠ The words are also carried as STRINGS here, one character per side letter, with the
#   map chosen ORDER-PRESERVING (`H < L1 < L2` as tuples, `a < b < c` as characters) so
#   that `canon_fast(w) == word_str(s315.canon(w))` is an equality and not a convention.
#   A 6.6e5-letter tuple of interned side strings costs ~5 MB in pointers alone; the
#   string costs 0.66 MB, and every downstream consumer (`assemble`'s grouping, `sha1`,
#   `len`, `merge`'s slicing/`zip`) is sequence-generic and does not notice.

_LET = {"H": "a", "L1": "b", "L2": "c"}


def word_str(word):
    return word if isinstance(word, str) else "".join(_LET[c] for c in word)


def least_rotation(s):
    """Booth's algorithm: index of the lexicographically least rotation of `s`, O(n)."""
    s = s + s
    n = len(s)
    f = [-1] * n
    k = 0
    for j in range(1, n):
        sj = s[j]
        i = f[j - k - 1]
        while i != -1 and sj != s[k + i + 1]:
            if sj < s[k + i + 1]:
                k = j - i - 1
            i = f[i]
        if sj != s[k + i + 1]:
            if sj < s[k]:
                k = j
            f[j - k] = -1
        else:
            f[j - k] = i + 1
    return k


def canon_fast(word):
    """`s315.canon` in O(n), returning the compact STRING form (see the note above)."""
    s = word_str(word)
    if not s:
        return s
    k = least_rotation(s)
    return s[k:] + s[:k]


EXACT_MISMATCH_MAX = 3000
ANCHOR = 24
_S315_MISMATCH = s315.best_rotation_mismatch     # bound BEFORE any patch: the reference


def mismatch_le3(w1, w2):
    """`s315.best_rotation_mismatch` in the ONLY regime its caller uses -- the caller
    (`merge_duplicate_cylinders`) tests `<= MERGE_MAX_MISMATCH = 3` and nothing else.  So
    this returns the exact minimum when it is `<= 3`, and `4` ("more than 3") otherwise.

    ⚠ THE ANCHOR ARGUMENT IS EXACT, NOT PROBABILISTIC, and that is the whole point: if
    some rotation of `w2` differs from `w1` in at most 3 positions, then at most 3 of the
    4 DISJOINT anchor windows taken from `w1` can contain a mismatch, so at least one
    anchor occurs VERBATIM in the doubled `w2` at the offset that realises that rotation.
    Searching for it therefore cannot miss the alignment; what is probabilistic is only
    how many SPURIOUS occurrences a 24-letter window has in a 3-letter alphabet
    (expectation `n/3^24`, i.e. none), and that is a cost, not an error.
    ⚠ Short words keep s315's exact quadratic routine -- there is nothing to gain and the
    anchor argument needs `n >= 4*ANCHOR`."""
    n = len(w1)
    if n != len(w2):
        return max(n, len(w2))
    if n <= EXACT_MISMATCH_MAX or n < 4 * ANCHOR:
        return _S315_MISMATCH(w1, w2)
    s1, s2 = word_str(w1), word_str(w2)
    dbl = s2 + s2[:n - 1]
    cands = set()
    for a in range(4):
        pos = (a * n) // 4
        pat = s1[pos:pos + ANCHOR]
        start = 0
        while True:
            occ = dbl.find(pat, start)
            if occ < 0 or occ >= n:
                break
            cands.add((occ - pos) % n)
            start = occ + 1
    best = 4
    for r in cands:
        rot = s2[r:] + s2[:r]
        d = 0
        for x, y in zip(s1, rot):
            if x != y:
                d += 1
                if d >= best:
                    break
        if d < best:
            best = d
    return best


def install_fast_word_ops():
    """⚠ DELIBERATE MONKEYPATCH, and it is the least-bad option.  `s315`'s merge logic is
    CALIBRATED (its `MERGE_MAX_MISMATCH = 3` separates 7/16's true duplicate from 5/7's
    false positive, s315 trap (iii)) and must stay the single implementation -- copying it
    here to swap one call would fork a load-bearing routine.  So the QUADRATIC primitive
    is replaced under it instead, by a function that agrees with it exactly wherever the
    merge criterion looks.  Idempotent; controlled by `word_ops_control`."""
    if getattr(s315.best_rotation_mismatch, "_s440_fast", False):
        return
    mismatch_le3._s440_fast = True
    s315.best_rotation_mismatch = mismatch_le3


def word_ops_control(n_random=300, verbose=True):
    """CONTROL, able to fail: `canon_fast` and `mismatch_le3` against s315's originals.
    Random words over the real 3-letter alphabet, plus the near-duplicate case the merge
    actually has to decide (a word with `d` planted mismatches at a random rotation)."""
    rng = random.Random(4402)
    letters = ["H", "L1", "L2"]
    bad_canon = bad_mm = 0
    for _ in range(n_random):
        m = rng.randint(1, 40)
        w = tuple(rng.choice(letters) for _ in range(m))
        if canon_fast(w) != word_str(s315.canon(w)):
            bad_canon += 1
    # long words: canon against a direct (but O(n^2)-free) min over rotations is not
    # available, so check the DEFINING property instead -- the result is a rotation of the
    # word and is <= every rotation -- on a sample of rotations.
    for _ in range(20):
        m = rng.randint(4000, 9000)
        s = "".join(rng.choice("abc") for _ in range(m))
        c = canon_fast(s)
        if len(c) != m or (s + s).find(c) < 0:
            bad_canon += 1
        for _ in range(40):
            r = rng.randrange(m)
            if c > s[r:] + s[:r]:
                bad_canon += 1
    # mismatch: planted-difference pairs, exact regime and anchor regime
    for _ in range(120):
        m = rng.choice([200, 1500, 4000, 20000])
        s = "".join(rng.choice("abc") for _ in range(m))
        r = rng.randrange(m)
        t = list(s[r:] + s[:r])
        d = rng.randint(0, 6)
        for p in rng.sample(range(m), d):
            t[p] = rng.choice([c for c in "abc" if c != t[p]])
        t = "".join(t)
        got = mismatch_le3(s, t)
        want = (_S315_MISMATCH(s, t) if m <= EXACT_MISMATCH_MAX
                else min(d, 4))
        # the planted count is an UPPER bound on the true minimum (another rotation may
        # do better), so only the <=3 verdict is compared, which is all the merge uses.
        if (got <= 3) != (want <= 3) or (got <= 3 and want <= 3 and got > want):
            bad_mm += 1
    if verbose:
        print(f"  word ops: canon_fast {bad_canon} mismatches, "
              f"mismatch_le3 {bad_mm} disagreements with s315", flush=True)
    return bad_canon == 0 and bad_mm == 0


# ---------------------------------------------------------------- the prong tracer

def trace_prong(B, geo, verts, Q, v0, u0, cap=CAP):
    """Launch from vertex `v0` along grid direction `u0` and record every bounce until
    the orbit reaches a vertex again.

    ⚠ Terminates at ANY of `O`/`R`/`A`, not at the two CONE points: an `R` prong is a
    regular closed leaf and comes back to `R` without ever meeting a cone point (see the
    module docstring).  Returns `(bounces, reason)`; `reason` is `vertex:<name>`,
    `cap`, `escaped` or `offgrid`."""
    px, py = verts[v0]
    th = u0 * math.pi / (2 * Q)
    vx, vy = math.cos(th), math.sin(th)
    u_in = u0
    out = []
    for _ in range(cap):
        cand = B._candidate_times(px, py, vx, vy)
        if not cand:
            return out, "escaped"
        t, side = min(cand, key=lambda z: z[0])
        px += t * vx
        py += t * vy
        _, _, base, tan = geo[side]
        w = (px - base[0]) * tan[0] + (py - base[1]) * tan[1]
        hit = None
        for nm, (qx, qy) in verts.items():
            if math.hypot(px - qx, py - qy) < VERT_TOL:
                hit = nm
                break
        vx, vy = B.reflect_velocity(side, vx, vy)
        u_out = _u_of(vx, vy, Q)
        if u_out is None:
            return out, "offgrid"
        out.append({"side": side, "w": w, "u_in": u_in, "u_out": u_out})
        if hit:
            return out, "vertex:" + hit
        u_in = u_out
    return out, "cap"


def boundary_map(P, Q, parity, cap=CAP, cap2=None, abort_on_open=False):
    """⇒⇒ THE INSTRUMENT.  Every cell boundary of every strip of one direction class,
    generated from the vertices.

    A bounce at `w` on side `σ` with incoming index `u_in` and outgoing `u_out` is a
    boundary of TWO strips: `(σ, u_out)`, whose backward orbit from `w` reaches the
    launch vertex, and `(σ, u_in + 2Q)`, the time-reverse, whose forward orbit does.
    Both are grazes and both are recorded; a boundary found twice dedups.

    ⚠⚠ A PRONG DOES NOT STOP AT THE OPPOSITE SIDE, AND THE FAN'S OWN GAPS ARE NOT
    CELLS.  It bounces until it reaches a vertex -- 885594 times at 7/20 -- and EVERY
    bounce is a boundary, each in whichever strip that bounce's outgoing direction picks
    out.  The fan's FIRST bounces therefore land in as many DIFFERENT strips as there are
    prongs (measured at 5/12: 2 of 2 from `O`, 6 of 6 from `R`, 3 of 3 from `A`, all
    different), so no two of them bound a common cell.  What bounds a cell is a pair of
    hits on the SAME side with the SAME outgoing direction, and those come from different
    prongs at very different DEPTHS: 5/12's perpendicular strip `H:41` has 12 cells cut
    by `R:17@0`, `A:31@1`, `R:13@4`, `A:35@5`, `R:19@10`, `A:33@14`, `R:23@54`, `O:1@89`,
    `R:21@209`, `R:15@255`, `O:3@418` -- one strip, depths spanning 0 to 418.
    ⇒ THE FINE STRUCTURE IS MADE BY DEEP BOUNCES, WHICH IS WHY CELL WIDTH TRACKS PRONG
    LENGTH AND NOT `Q`: at 5/12 the longest prong is 838 bounces and the narrowest cell
    is 4.6e-03; at 7/20 the longest is 885594 and the narrowest is 6.0e-06.
    ⚠ And the alternation in that list is not decoration: the two crossings of one
    cylinder are bounded OUTSIDE by cone-point (`O`/`A`) hits -- genuine saddle
    connections -- and separated INSIDE by an `R` hit, the regular closed leaf.  Drop `R`
    and the two crossings fuse into one interval of twice the width, which is why the
    heights would come out 2x too big.

    ⇒ TWO-STAGE CAP (s440).  `cap2 > cap` re-traces ONLY the prongs that did not
    terminate at `cap`, which is how the cap can be raised 10x for free ([OPS-205]:
    the raise costs nothing on a completely periodic class and 10x the wasted work on a
    class with a genuine minimal component, where nothing terminates at any cap).
    `abort_on_open` then stops the whole class at the FIRST prong still open at `cap2` --
    such a class is `PRONG_CAPPED` and unscorable, so the remaining prongs are ~2e6
    bounces each of certified-useless work.  ⚠ It makes `prongs_not_terminated` a LOWER
    BOUND on that row rather than the full list; `n_open_at_cap1` is the complete count
    at the stage-1 cap and is not truncated."""
    geo = s437.geometry(P, Q)
    verts = vertices(P, Q)
    B = RightTriangleBilliard((P / Q) * (math.pi / 2))
    bmap = defaultdict(lambda: array("d"))
    reasons, n_bounces = {}, 0

    def absorb(bounces):
        nonlocal n_bounces
        n_bounces += len(bounces)
        for b in bounces:
            bmap[(b["side"], b["u_out"])].append(b["w"])
            bmap[(b["side"], (b["u_in"] + 2 * Q) % (4 * Q))].append(b["w"])

    retry = []
    for v0 in ("O", "R", "A"):
        for u0 in prongs(v0, P, Q, parity):
            bounces, why = trace_prong(B, geo, verts, Q, v0, u0, cap)
            if why == "cap" and cap2 and cap2 > cap:
                retry.append((v0, u0))            # bounces discarded; re-traced below
                continue
            reasons[f"{v0}:{u0}"] = {"reason": why, "n": len(bounces)}
            absorb(bounces)
    n_open_at_cap1 = len(retry) + sum(1 for v in reasons.values()
                                      if not v["reason"].startswith("vertex"))
    for (v0, u0) in retry:
        bounces, why = trace_prong(B, geo, verts, Q, v0, u0, cap2)
        reasons[f"{v0}:{u0}"] = {"reason": why, "n": len(bounces), "stage": 2}
        absorb(bounces)
        if abort_on_open and not why.startswith("vertex"):
            return bmap, reasons, n_bounces, n_open_at_cap1, True
    return bmap, reasons, n_bounces, n_open_at_cap1, False


# ---------------------------------------------------------------- cells from boundaries

def _seg0(B, geo, side, x, u, Q):
    """First-segment length from `x` on `side` in direction `u` -- ONE step, not a full
    trace.  (`s437.discover_cylinder_intervals` spends a whole closure trace per endpoint
    to read this; it is the first `_candidate_times` and nothing more.)"""
    _, _, base, tan = geo[side]
    px, py = base[0] + x * tan[0], base[1] + x * tan[1]
    th = u * math.pi / (2 * Q)
    cand = B._candidate_times(px, py, math.cos(th), math.sin(th))
    return min(c[0] for c in cand) if cand else None


def strip_cells(B, geo, Q, side, u, bnds, cap):
    """Cells of the strip `(side,u)` from its EXACT boundary list.  Same record shape as
    `s437.discover_cylinder_intervals` so `s437.assemble` consumes it unchanged."""
    ell = geo[side][0]
    pts = sorted([0.0, ell] + [w for w in bnds if 0.0 < w < ell])
    cuts = [pts[0]]
    for w in pts[1:]:
        if w - cuts[-1] > DEDUP_TOL:
            cuts.append(w)
    cuts[-1] = ell
    out, min_gap = [], float("inf")
    for lo, hi in zip(cuts[:-1], cuts[1:]):
        width = hi - lo
        if width < MIN_WIDTH:
            continue
        min_gap = min(min_gap, width)
        mid = 0.5 * (lo + hi)
        w_, c, _, ok = s437.trace_full_period(B, geo, Q, side, mid, u, cap)
        eps = width * 1e-6
        s_lo = _seg0(B, geo, side, lo + eps, u, Q)
        s_hi = _seg0(B, geo, side, hi - eps, u, Q)
        seg0 = (s_lo + s_hi) / 2.0 if (s_lo and s_hi) else None
        canon = s315.canon(w_) if ok else None
        out.append({"side": side, "u": u, "lo": lo, "hi": hi, "width": width,
                    "word": w_, "canon": canon,
                    "key": (hashlib.sha1(repr(canon).encode()).hexdigest()[:10]
                            if canon is not None else None),
                    "circumference": c, "seg0": seg0, "ok": ok and seg0 is not None})
    return out, min_gap


# ---------------------------------------------------------------- first-return form

def first_return(B, geo, Q, side0, x0, u0, cap):
    """Trace from `x0` on `(side0,u0)` to its FIRST RETURN to that same strip.
    Returns `(w, word, length, ok)` -- landing position, the side letters traversed, the
    arc length.  This is `s437.trace_full_period` with the position test dropped: that
    function runs to the FULL PERIOD, i.e. it re-walks the whole cylinder once per cell."""
    _, _, base0, tan0 = geo[side0]
    px = base0[0] + x0 * tan0[0]
    py = base0[1] + x0 * tan0[1]
    th = u0 * math.pi / (2 * Q)
    vx, vy = math.cos(th), math.sin(th)
    word, arc = [], 0.0
    for _ in range(cap):
        cand = B._candidate_times(px, py, vx, vy)
        if not cand:
            return None, None, arc, False
        t, side = min(cand, key=lambda z: z[0])
        arc += t
        px += t * vx
        py += t * vy
        word.append(side)
        vx, vy = B.reflect_velocity(side, vx, vy)
        u = _u_of(vx, vy, Q)
        if u is None:
            return None, None, arc, False
        if side == side0 and u == u0:
            _, _, base, tan = geo[side]
            w = (px - base[0]) * tan[0] + (py - base[1]) * tan[1]
            return w, tuple(word), arc, True
    return None, None, arc, False


def strip_cells_fast(B, geo, Q, side, u, bnds, cap):
    """⇒⇒ THE COST FIX, and it is a change of FACTORISATION, not a tolerance.
    `strip_cells` traces every cell to closure, so a cylinder crossing the strip `k_c`
    times is re-walked `k_c` times -- total work `Σ_c k_c · period_c`.  Here each cell is
    traced only to its FIRST RETURN; the induced first-return map on the cells is a
    permutation, its CYCLES are exactly the cylinder traces ([NCYL-247]: a cell of the
    strip is one crossing of a cylinder), the cycle LENGTH is `k_c`, and the closure word
    is the concatenation of the return words around the cycle.  Total work
    `Σ_c period_c`, i.e. a saving of the mean `k_c` -- which is `~10` at `Q = 14` and
    GROWS with `Q`, since `Σ_c k_c h_c = |σ|·flux` is fixed while `h_min` shrinks.
    Nothing is approximated: the words produced are the same words."""
    ell = geo[side][0]
    pts = sorted([0.0, ell] + [w for w in bnds if 0.0 < w < ell])
    cuts = [pts[0]]
    for w in pts[1:]:
        if w - cuts[-1] > DEDUP_TOL:
            cuts.append(w)
    cuts[-1] = ell
    cells, min_gap = [], float("inf")
    for lo, hi in zip(cuts[:-1], cuts[1:]):
        if hi - lo < MIN_WIDTH:
            continue
        min_gap = min(min_gap, hi - lo)
        cells.append({"lo": lo, "hi": hi, "width": hi - lo})
    if not cells:
        return [], min_gap

    # (1) one first-return trace per cell
    for c in cells:
        mid = 0.5 * (c["lo"] + c["hi"])
        w, word, arc, ok = first_return(B, geo, Q, side, mid, u, cap)
        c.update({"ret_w": w, "ret_word": word, "ret_arc": arc, "ok": ok})

    # (2) the induced permutation: which cell does the return land in?
    los = [c["lo"] for c in cells]
    def locate(w):
        if w is None:
            return None
        i = bisect.bisect_right(los, w) - 1
        if 0 <= i < len(cells) and cells[i]["lo"] - DEDUP_TOL <= w <= cells[i]["hi"] + DEDUP_TOL:
            return i
        return None
    nxt = [locate(c["ret_w"]) for c in cells]

    # (3) cycles = cylinders; the closure word is the concatenation around the cycle
    seen, out = [False] * len(cells), []
    for i in range(len(cells)):
        if seen[i]:
            continue
        cyc, j = [], i
        while j is not None and not seen[j]:
            seen[j] = True
            cyc.append(j)
            j = nxt[j]
        closed = (j == i) and all(cells[t]["ok"] for t in cyc)
        # ⚠ a cycle that does not RETURN to its start is not a cylinder -- mark every
        # cell of it bad rather than inventing one ([NCYL-020]).
        words = tuple(x for t in cyc for x in (cells[t]["ret_word"] or ()))
        circ = sum(cells[t]["ret_arc"] for t in cyc)
        for pos, t in enumerate(cyc):
            c = cells[t]
            # rotate the concatenated word so it starts at THIS cell
            off = sum(len(cells[cyc[q]]["ret_word"] or ()) for q in range(pos))
            wrd = words[off:] + words[:off] if closed else None
            canon = s315.canon(wrd) if closed else None
            eps = c["width"] * 1e-6
            s_lo = _seg0(B, geo, side, c["lo"] + eps, u, Q)
            s_hi = _seg0(B, geo, side, c["hi"] - eps, u, Q)
            seg0 = (s_lo + s_hi) / 2.0 if (s_lo and s_hi) else None
            out.append({"side": side, "u": u, "lo": c["lo"], "hi": c["hi"],
                        "width": c["width"], "word": wrd, "canon": canon,
                        "key": (hashlib.sha1(repr(canon).encode()).hexdigest()[:10]
                                if canon is not None else None),
                        "circumference": circ if closed else None, "seg0": seg0,
                        "ok": closed and seg0 is not None})
    out.sort(key=lambda z: z["lo"])
    return out, min_gap


# ------------------------------------------------- claim mode: O(#cylinders) traces
# ⇒⇒ s440, QUEUE ITEM (i).  THE COST OF A ROW WAS NEVER THE PRONG CAP AND WAS NEVER THE
# DEPTH -- it is that the readers above traverse EVERY CELL, and a 6.6e5-bounce saddle
# connection deposits ~2e5 cells in a single strip ([OPS-205]).  Every one of those cells
# is then canonicalised and hashed, at O(period) each, which is the hour.
#
# ⇒⇒ THE FACTORISATION THAT REMOVES IT.  A cylinder's cells in a strip are ONE first-
# return cycle ([NCYL-247], `strip_cells_fast`'s docstring), so a SINGLE closed orbit of
# that cylinder visits ALL of them -- and, tracing it once to closure, it visits every
# cell of that cylinder in EVERY strip of the class at the same time.  So: trace one cell
# to closure, CLAIM every cell its orbit lands in, move to the next unclaimed cell.  The
# number of traces is then the number of cylinders (plus one per `R`-split word variant),
# not the number of cells, and `k_c` and the met-sets come out EXACT as a by-product --
# they are the claim counts.
#
#   ⚠ WHY NOT THE WIDTH-GROUPING ROUTE [OPS-205] SKETCHED.  Grouping a strip's cells by
#   width first (`141002 -> 8` groups there) needs a tolerance, and [OPS-205] measured
#   that tolerance to be a PLATEAU (`1e-14 -> 6655` groups, `1e-9 -> 6`, the answer `8`
#   only over `1e-12..1e-11`) -- a plateau that would have to be re-established per row.
#   Claiming needs no tolerance at all: the orbit says which cells are its own.  The width
#   is used only as a REJECTION test on an individual claim (below), where it is a
#   guard rather than the classifier.
#
#   ⚠ AND THE TWO FAILURE DIRECTIONS ARE NOT SYMMETRIC, which is why this is safe.
#   UNDER-claiming (float drift moves a landing into the neighbouring cell, the width test
#   rejects it) costs one extra trace and is then absorbed: the re-trace produces the same
#   canonical word and `assemble` groups the two together.  OVER-claiming (a cell of a
#   DIFFERENT cylinder claimed) would hide a cylinder -- and that is exactly what the
#   completeness certificate `Σ_c circ_c·h_c = Q·cot α` ([NCYL-247]/[OPS-204]) is for.  So
#   one direction self-heals and the other is certified.


def _cut_list(geo, side, bnds):
    """The strip's exact cut positions -- factored verbatim out of `strip_cells` /
    `strip_cells_fast` so all three readers cut a strip identically."""
    ell = geo[side][0]
    pts = sorted([0.0, ell] + [w for w in bnds if 0.0 < w < ell])
    cuts = [pts[0]]
    for w in pts[1:]:
        if w - cuts[-1] > DEDUP_TOL:
            cuts.append(w)
    cuts[-1] = ell
    return cuts


def trace_closure_claim(B, geo, Q, side0, x0, u0, cap, strips, cid, h_c):
    """Trace the closed orbit through ONE cell to full closure.  Two things come back
    from the one trace, and they are deliberately DIFFERENT quantities:

      * `visits[tag]` -- how many times the orbit CROSSES strip `tag`.  This is
        `k_c(σ,u)` by definition ([NCYL-247]: the crossing multiplicity, = the number of
        cells of that strip belonging to `c`), it is counted off the orbit itself, and it
        does not depend on the cell bookkeeping at all.
      * the CLAIMS -- cells marked as belonging to this cylinder so that the reader does
        not trace them again.  A claim is refused when the landing cell's width disagrees
        with `h_c/flux` beyond `WIDTH_REL_TOL`, i.e. when the landing is not safely inside
        a cell of this cylinder's height; refusing costs at most a repeated trace, and
        `rep[tag]` then falls back to the derived width.

    ⚠ THE TWO MUST NOT BE CONFLATED: `k_c` off the claims would undercount every strip
    whose partition is finer than this cylinder's crossings, and that is exactly where
    the reader was spending its time before ([OPS-205] measured the cost, not this).

    ⇒⇒ THE CLOSURE TEST SCALES WITH THE CYLINDER (s443, [OPS-211]), AND THE SCALE IS NOT
    A TUNED CONSTANT -- IT IS THE LAUNCH CELL ITSELF.  The launch is the MIDPOINT of a
    cell of width `h_c/flux0` ([NCYL-247]), so `|w - x0| < h_c/(2*flux0)` says exactly
    *the return landed in the cell it was launched from*, and nothing else.  ⚠ THE FLAT
    `POS_TOL = 1e-7` DOES NOT SAY THAT: a cylinder with `h < POS_TOL` has its ENTIRE
    cross-section inside the tolerance, so the test fires at the first return that is
    merely NEARBY -- measured at `10/37` par 1, that return sits `0.99976` of a full cell
    width BELOW `x0`, i.e. outside the launch cell altogether, and the true closure is at
    exactly twice the step with an offset of `3.1e-11` (pure drift).  Hence the halved
    period, the halved `circ`, the perpendicular `k = 1` for `2`, and the odd word length
    that is the tell.  ⚠ `POS_TOL` is KEPT as a ceiling, so this is a strict tightening
    and a NO-OP wherever the launch cell is wider than `2*POS_TOL`.
    ⚠ WHAT IT DOES NOT FIX: if the launch cell is itself a MERGED pair of crossings the
    partition failed to separate ([NCYL-271]'s `DEDUP_TOL` floor), `h_c` is wrong at the
    source and so is the window.  That is the other float floor, not this one."""
    _, _, base0, tan0 = geo[side0]
    px = base0[0] + x0 * tan0[0]
    py = base0[1] + x0 * tan0[1]
    th = u0 * math.pi / (2 * Q)
    vx, vy = math.cos(th), math.sin(th)
    word, circ, seg0 = [], 0.0, None
    visits, rep, sw, cnt = {}, {}, {}, {}
    n_refused, dev_max = 0, 0.0
    pos_tol = POS_TOL
    if SCALED_CLOSURE:
        st0 = strips.get((side0, u0))
        if st0 is not None and st0["flux"] > 0.0:
            pos_tol = min(POS_TOL, 0.5 * h_c / st0["flux"])
    for _ in range(cap):
        cand = B._candidate_times(px, py, vx, vy)
        if not cand:
            return None, circ, seg0, visits, rep, sw, cnt, False, n_refused, dev_max
        t, side = min(cand, key=lambda z: z[0])
        if seg0 is None:
            seg0 = t
        circ += t
        px += t * vx
        py += t * vy
        word.append(side)
        vx, vy = B.reflect_velocity(side, vx, vy)
        u = _u_of(vx, vy, Q)
        if u is None:
            return None, circ, seg0, visits, rep, sw, cnt, False, n_refused, dev_max
        _, _, base, tan = geo[side]
        w = (px - base[0]) * tan[0] + (py - base[1]) * tan[1]
        tag = (side, u)
        st = strips.get(tag)
        if st is not None:
            visits[tag] = visits.get(tag, 0) + 1
            cuts = st["cuts"]
            j = bisect.bisect_right(cuts, w) - 1
            if 0 <= j < st["n"] and st["claim"][j] == -1:
                exp = h_c / st["flux"]
                dev = abs((cuts[j + 1] - cuts[j]) - exp) / exp
                if dev <= WIDTH_REL_TOL:
                    st["claim"][j] = cid
                    if dev > dev_max:
                        dev_max = dev
                    sw[tag] = sw.get(tag, 0.0) + (cuts[j + 1] - cuts[j])
                    cnt[tag] = cnt.get(tag, 0) + 1
                    if tag not in rep:
                        rep[tag] = j
                else:
                    n_refused += 1
        if side == side0 and u == u0 and abs(w - x0) < pos_tol:
            return (tuple(word), circ, seg0, visits, rep, sw, cnt, True, n_refused,
                    dev_max)
    return None, circ, seg0, visits, rep, sw, cnt, False, n_refused, dev_max


def class_claim(B, geo, Q, par, bmap, perp_u, area_exact, cap=CAP_BIG,
                trace_cap=TRACE_CAP, verbose=False, floor=None, out_strips=None,
                out_traces=None):
    """⇒⇒ THE COUNT READER.  One parity class, `O(#cylinders)` closure traces.  Returns
    `(slices, strip_rows, diag)` with `slices` in `s437.assemble`'s record shape, one per
    (traced cylinder, strip) and carrying the crossing count `k`.

    ⇒ THE STOP IS THE CERTIFICATE, NOT THE STRIP LIST ([OPS-204]/[OPS-205]).  After every
    trace the class's accounted area `Σ_c circ_c·h_c` is compared with `Q·cot α`: a
    cylinder contributes its circumference (from the closure) times its height (`width·
    flux` from one crossing, [NCYL-247]), so when the sum reaches the surface area there
    is nothing left to find and the remaining cells are further crossings of cylinders
    already counted.  ⚠ Kac's total CANNOT be used for this -- it is an integral over the
    whole 3-side cross-section and only reaches 1 once every cell of every strip has been
    materialised, which is the cost being removed.
    ⚠ AND THE STOP IS SOUND IN BOTH DIRECTIONS: a MISSED cylinder leaves the sum short, a
    cylinder counted TWICE pushes it past 1, and either way the row fails its own verdict
    check and is reported `INCOMPLETE` rather than scored.

    ⚠ EVERY strip of the class is cut (that is a sort, not a trace); what is dropped
    relative to `strip_cells*` is the per-CELL word, seg0 and Kac contribution.  The Kac
    total coming out of this reader is therefore an ESTIMATE and must never be the
    verdict column.

    ⇒ `out_strips` / `out_traces` (s441): containers handed in to receive the internal
    `strips` table (`cuts` / `claim` / `flux` per tag) and the `traces` list, indexed by
    the `cid` the `claim` array stores.  READ-ONLY for the caller and purely additive --
    the `claim` array is the full (cylinder -> every cell it owns) map, which `slices`
    collapses to one representative per (cylinder, strip) via `rep`.  `s441_graze_depth`
    needs all of them, since its object is a MINIMUM over a cylinder's crossings."""
    t_phase = time.time()
    strips, min_gap = {}, float("inf")
    for side in geo:
        for u in range(par, 4 * Q, 2):
            fl = s437.flux_weight(side, u, Q, geo)
            if fl <= 1e-12:
                continue
            cuts = _cut_list(geo, side, bmap.get((side, u), []))
            n = len(cuts) - 1
            if n < 1:
                continue
            strips[(side, u)] = {"cuts": cuts, "claim": array("i", [-1]) * n,
                                 "flux": fl, "n": n, "side": side, "u": u}
            for i in range(n):
                g = cuts[i + 1] - cuts[i]
                if MIN_WIDTH <= g < min_gap:
                    min_gap = g
    # ⇒⇒ PERPENDICULAR STRIPS FIRST, THEN FEWEST CELLS -- AND THIS IS A CORRECTNESS
    # ORDER, NOT A SPEED ONE.  MEASURED, both ways, at 7/20's lone class: launched from
    # the perpendicular strip `H:67` (20 cells = the 10 cylinders crossed twice each,
    # [NCYL-248]) the reader takes 10 traces, refuses NO claim and closes the certificate
    # to 5e-10; launched from the cheapest strip instead it takes 17 traces, finds 12
    # word-families, accepts claims at the very edge of `WIDTH_REL_TOL` and assembles to
    # `area_ratio_hc = 0.73`.
    # ⇒ THE REASON, and it is the thing to know before touching this: the launch cell's
    # width is what the whole trace's height `h_c` is calibrated from, so the launch cell
    # must be a FULL crossing of the cylinder.  On a strip cut finer than the cylinder's
    # own height (`R`-grazes split a cell of constant WORD while leaving the cylinder
    # intact, [NCYL-258]) a cell is a FRAGMENT, `h_c` comes out too small, and the
    # crossings it then claims are the wrong ones.  A perpendicular strip cannot be
    # fragmented that way -- `k ∈ {1,2}` there ([NCYL-248]).
    # ⚠ The safety net if this heuristic ever fails is the certificate, and it FIRED on
    # the bad ordering above: the row read `INCOMPLETE` and was not scored ([NCYL-020]).
    if out_strips is not None:
        out_strips.update(strips)
    order = sorted(strips, key=lambda t: (0 if t[1] == perp_u[t[0]] else 1,
                                          strips[t]["n"], t[0], t[1]))
    if verbose:
        print(f"    cut {sum(st['n'] for st in strips.values())} cells over "
              f"{len(strips)} strips [{time.time() - t_phase:.1f}s]", flush=True)
    perp_tags = [t for t in strips if t[1] == perp_u[t[0]]]

    def perp_complete():
        """⇒⇒ THE EXACT COMPLETENESS TEST, and it is an integer one.  Every cell of
        every PERPENDICULAR strip of the class has been attributed to a traced cylinder.
        A cylinder can then be missing only if it crosses no perpendicular strip at all
        -- which is [NCYL-243]'s measured negative (no cylinder is oblique-only) -- and
        the area certificate is the second, independent guard against that."""
        for t in perp_tags:
            cl = strips[t]["claim"]
            for i in range(strips[t]["n"]):
                if cl[i] == -1:
                    return False
        return True

    traces, by_canon = [], {}
    n_bad = n_refused = n_retrace = 0
    dev_max, acc, closed, capped = 0.0, 0.0, False, False
    for tag in order:
        if (closed and floor is None) or capped:
            break
        st = strips[tag]
        cuts, claim = st["cuts"], st["claim"]
        for i in range(st["n"]):
            if claim[i] != -1:
                continue
            lo, hi = cuts[i], cuts[i + 1]
            if hi - lo < (floor or MIN_WIDTH):
                claim[i] = -2
                continue
            if len(traces) >= trace_cap:
                capped = True
                break
            cid = len(traces)
            claim[i] = cid
            t_tr = time.time()
            h_c = (hi - lo) * st["flux"]
            # ⇒ THE CLOSURE CAP MUST EXCEED THE *TRUE* PERIOD, WHICH THE SCALED CLOSURE
            # TEST MAKES UP TO TWICE THE OLD ONE (s443, [OPS-211]).  ⚠ MEASURED, and it
            # is the only way the fix can LOSE a cylinder: at `26/55` par 1 two thin
            # cylinders have a true period of `2156254` against `CAP_BIG = 2000000`, so
            # both traces ran out and the row went `29 -> 27` with the certificate short
            # by exactly their combined area share.  At `cap = 5e6` both come back at
            # `|w| = 2156254`, `circ` doubled, perpendicular `k = 1 -> 2`.
            # ⚠ IT COSTS NOTHING ON A ROW THAT CLOSES: the loop exits at closure, so the
            # cap is a ceiling on failure, not a workload.
            (word, circ, seg0, visits, rep, sw, cnt, ok, nref,
             dev) = trace_closure_claim(
                B, geo, Q, st["side"], 0.5 * (lo + hi), st["u"],
                2 * cap if SCALED_CLOSURE else cap, strips, cid, h_c)
            n_refused += nref
            dev_max = max(dev_max, dev)
            if not ok:
                claim[i] = -3
                n_bad += 1
                traces.append(None)
                continue
            rep.setdefault(tag, i)                  # the launch cell, claimed above
            sw.setdefault(tag, 0.0)
            sw[tag] += hi - lo
            cnt[tag] = cnt.get(tag, 0) + 1
            n_cl = sum(cnt.values())
            # ⇒⇒ THE HEIGHT IS THE LAUNCH CELL'S, AND AVERAGING THE OTHER CROSSINGS
            # MAKES IT WORSE -- MEASURED, and it is the opposite of what one would
            # assume.  The `k_c ~ 10^4` crossings look like independent estimates of
            # `h_c/flux` whose mean should be ~100x tighter; at 7/20's lone class the
            # mean moves the certificate residual from 4.7e-10 to 1.6e-06, i.e. 3000x
            # WORSE, and past `CERT_TOL`, so the reader then hunts for a cylinder that
            # is not missing.  ⇒ the deviations are ONE-SIDED: a crossing cell can be
            # cut short by a further graze but never runs long, so the sample mean is
            # biased down while the launch cell -- a full crossing of a PERPENDICULAR
            # strip -- is not.  `h_mean_rel` records the gap; it is a diagnostic.
            h_mean = (sum(strips[t]["flux"] * v for t, v in sw.items()) / n_cl
                      if n_cl else h_c)
            canon = canon_fast(word)
            tr = {"circ": circ, "seg0": seg0, "visits": visits, "rep": rep,
                  "canon": canon, "h_c": h_c, "n_claimed": n_cl,
                  "h_mean_rel": (h_mean - h_c) / h_c if h_c else None,
                  "wordlen": len(word), "word": word}
            traces.append(tr)
            if canon in by_canon:
                n_retrace += 1                      # same family reached twice: no area
            else:
                by_canon[canon] = tr
                acc += circ * tr["h_c"]
                closed = (abs(acc / area_exact - 1.0) < CERT_TOL
                          and perp_complete())
            if verbose:
                print(f"    trace {cid} {tag} cell {i}: |w|={len(word)} "
                      f"strips={len(visits)} k={sum(visits.values())} refused={nref} "
                      f"acc={acc / area_exact:.12f} [{time.time() - t_tr:.1f}s]",
                      flush=True)
            if closed and floor is None:
                break

    if out_traces is not None:
        out_traces.extend(traces)
    slices = []
    for tr in traces:
        if tr is None:
            continue
        key = hashlib.sha1(repr(tr["canon"]).encode()).hexdigest()[:10]
        for tag, k in tr["visits"].items():
            st = strips[tag]
            j = tr["rep"].get(tag)
            # ⚠ the width recorded is the cylinder's AVERAGED crossing width
            # `h_c/flux`, not the representative cell's own; `lo`/`hi` still locate that
            # cell.  When every landing on this strip was refused (the strip cuts this
            # cylinder finer than its own height) there is no representative at all --
            # the crossing is still REAL, so it is recorded and FLAGGED rather than
            # dropped from the met-set.
            width = tr["h_c"] / st["flux"]
            derived = j is None
            lo = hi = None
            if j is not None:
                lo, hi = st["cuts"][j], st["cuts"][j + 1]
            seg0 = None
            if lo is not None:
                eps = (hi - lo) * 1e-6
                a = _seg0(B, geo, st["side"], lo + eps, st["u"], Q)
                b = _seg0(B, geo, st["side"], hi - eps, st["u"], Q)
                seg0 = 0.5 * (a + b) if (a and b) else None
            slices.append({"side": st["side"], "u": st["u"], "lo": lo, "hi": hi,
                           "width": width, "word": tr["word"], "canon": tr["canon"],
                           "key": key, "circumference": tr["circ"],
                           "seg0": seg0, "ok": True, "flux": st["flux"],
                           "k": k, "width_derived": derived,
                           "area_contrib": (k * width * st["flux"] * seg0
                                            if seg0 else 0.0)})
    strip_rows, n_unclaimed = [], 0
    for tag, st in strips.items():
        cl, nb2, nu = st["claim"], 0, 0
        for i in range(st["n"]):
            v = cl[i]
            if v == -1:
                nu += 1
            elif v < -1:
                nb2 += 1
        n_unclaimed += nu
        strip_rows.append({"side": st["side"], "u": st["u"], "parity": par,
                           "is_perp": st["u"] == perp_u[st["side"]],
                           "n_cells": st["n"], "n_bad": nb2, "n_unclaimed": nu,
                           "sum_width": st["cuts"][st["n"]] - st["cuts"][0],
                           "side_length": geo[st["side"]][0]})
    diag = {"n_traces": len(traces), "n_trace_failed": n_bad, "n_retrace": n_retrace,
            "exhaustive_floor": floor,
            "h_mean_rel_worst": max((abs(t["h_mean_rel"]) for t in traces
                                     if t and t["h_mean_rel"] is not None), default=0.0),
            "n_distinct_words": len(by_canon), "width_dev_max": dev_max,
            "n_claims_refused": n_refused, "n_cells_unclaimed": n_unclaimed,
            "certificate_closed": closed, "area_from_hc_direct": acc / area_exact,
            "perp_complete": bool(perp_tags) and perp_complete(),
            "n_perp_strips": len(perp_tags),
            "trace_capped": capped, "n_cells": sum(st["n"] for st in strips.values()),
            "min_cell_gap": None if min_gap == float("inf") else min_gap,
            "n_strips": len(strips), "seconds": round(time.time() - t_phase, 2)}
    return slices, strip_rows, diag



# ---------------------------------------------------------------- k_c

def _attach_kc(res, slices, geo, Q):
    """⇒ `k_c(σ,u)`, THE OBLIQUE `chi` (queue item (e)) -- the CROSSING MULTIPLICITY of
    cylinder `c` in strip `(σ,u)`, which by [NCYL-247] is just how many CELLS of that
    strip belong to `c`.

    ⚠ It cannot be read off `met_strips`: `assemble` stores that as a SET of `(side,u)`,
    so a cylinder crossed 44 times and one crossed once look identical there.  It is
    recomputed from the slices through the cell-key digests, which `assemble` re-unions
    across a merge -- so a merged cylinder gets its crossings summed, not dropped (the
    same trap s437 records for `met_strips`/`height_derived`).

    Also attaches the C4 flux identity per strip: `Σ_c k_c·h_c == |σ|·flux`
    ([NCYL-247], 402/402 at s437).  Written as a RESIDUAL, not a boolean."""
    by_key = {}
    for s in slices:
        by_key.setdefault(s["key"], []).append(s)
    hgt = {}
    for c in res["cylinders"]:
        per = {}
        for k in (c.get("keys") or []):
            for s in by_key.get(k, []):
                tag = f"{s['side']}:{s['u']}"
                # ⚠ `s["k"]` (s440 claim mode) is the crossing count the slice STANDS
                # FOR; in the per-cell readers each slice is one crossing and it is 1.
                per[tag] = per.get(tag, 0) + s.get("k", 1)
        c["k_per_strip"] = per
        c["k_multiset"] = sorted(per.values())
        c["k_max"] = max(per.values(), default=0)
        for tag in per:
            hgt.setdefault(tag, 0.0)
    # C4: per strip, sum_c k_c h_c against |sigma| * flux
    chk = {}
    for c in res["cylinders"]:
        h = c["height_derived"]
        for tag, k in c["k_per_strip"].items():
            chk[tag] = chk.get(tag, 0.0) + k * (h or 0.0)
    res["flux_identity"] = {}
    for tag, tot in chk.items():
        side, u = tag.split(":")
        want = geo[side][0] * s437.flux_weight(side, int(u), Q, geo)
        res["flux_identity"][tag] = abs(tot - want) / max(want, 1e-30)
    res["flux_identity_worst"] = max(res["flux_identity"].values(), default=None)
    res["k_max_overall"] = max((c["k_max"] for c in res["cylinders"]), default=0)

    # ⇒ THE STRIP-LOCAL COMPLETENESS CERTIFICATE, and it is what lets a count be
    # certified from a FEW strips instead of the whole boundary.  Kac's total
    # (`Σ width·flux·seg0`) is an integral over the ENTIRE 3-side cross-section, so it
    # reads 1 only when every strip of the class has been materialised -- which is why
    # dropping strips cannot be certified by it.  [NCYL-247] gives the alternative:
    # a cell's width IS `h_c/flux`, so `width·flux` recovers the cylinder HEIGHT from ANY
    # SINGLE crossing, and `area_c = circumference_c · h_c`.
    # ⚠ This is NOT s315 trap (i)'s `height × circumference` overcount: that one uses the
    # strip's TOTAL width for the cylinder, i.e. `k_c·h_c`, and so overcounts by exactly
    # the crossing multiplicity (11.9x at 3/7).  One crossing, not the sum.
    # ⚠ The spread of `width·flux` across a cylinder's cells is [NCYL-247] re-measured on
    # every row here and is reported, not assumed.
    tot, spread = 0.0, 0.0
    for c in res["cylinders"]:
        # ⚠ SUM OVER THE PRE-MERGE GROUPS, MEDIAN WITHIN ONE (fixed s440).  A key is one
        # canonical word, i.e. one first-return cycle; a MERGED cylinder ([NCYL-258]: an
        # `R` graze splits a cell of constant word while leaving the cylinder intact)
        # carries two cycles whose heights ADD to the cylinder's height, so a median over
        # the pooled cells returns one of the two halves and the certificate reads short.
        # With one key this is exactly the previous expression.
        h, seen_any = 0.0, False
        for k in (c.get("keys") or []):
            hs = [s["width"] * s["flux"] for s in by_key.get(k, [])]
            if not hs:
                continue
            seen_any = True
            hk = sorted(hs)[len(hs) // 2]
            h += hk
            # the [NCYL-247] spread is WITHIN one cycle, where every crossing must give
            # the same height; across cycles a difference is the R-split, not an error.
            spread = max(spread, (max(hs) - min(hs)) / max(hk, 1e-30))
        if not seen_any:
            continue
        c["h_from_width"] = h
        cm = 0.5 * (c["circumference_min"] + c["circumference_max"])
        c["area_from_hc"] = cm * h
        tot += c["area_from_hc"]
    res["area_from_hc_total"] = tot
    res["h_from_width_spread"] = spread


# ---------------------------------------------------------------- per row

def run_row(P, Q, cap=CAP, strip_budget=None, fast=True, mode="full", cap2=None,
            floor=None):
    """Both classes, exactly.  `strip_budget=None` materialises every strip (the
    like-for-like setting for H1); an int stops adding strips to a class once its Kac
    area closes, which is sound HERE and would not be under sampling -- the boundary set
    is exhaustive by construction, so `area_ratio = 1` cannot be hiding a fused cell.

    ⇒ `mode="claim"` (s440) reads the class with `class_claim` instead: every strip, but
    `O(#cylinders)` closure traces instead of one per cell.  It is the only mode that
    finishes a row whose prongs run to `10^5`+ bounces, and it keeps `k_c`, the met-sets
    and `C_met` exact; what it gives up is the per-cell Kac integrand (so its
    `total_area_ratio` is an estimate -- the verdict column is `area_ratio_hc`).
    `cap2 > cap` re-traces only the prongs that did not terminate at `cap`."""
    t0 = time.time()
    alpha = (P / Q) * (math.pi / 2)
    B = RightTriangleBilliard(alpha)
    geo = s437.geometry(P, Q)
    area_exact = Q / math.tan(alpha)
    perp_u = {sd: s437.perp_index(sd, P, Q) for sd in ("L1", "L2", "H")}
    if mode == "claim":
        install_fast_word_ops()

    out = {"P": P, "Q": Q, "alpha": alpha, "area_exact": area_exact,
           "perp_u": perp_u, "cap": cap, "cap2": cap2, "mode": mode, "fast": fast,
           "classes": {}, "strips": []}
    for par in (0, 1):
        bmap, reasons, nb, n_open1, aborted = boundary_map(
            P, Q, par, cap, cap2=cap2, abort_on_open=(mode == "claim"))
        capped = [k for k, v in reasons.items()
                  if not v["reason"].startswith("vertex")]
        if mode == "claim":
            slices, strip_rows, diag = ([], [], {}) if capped else class_claim(
                B, geo, Q, par, bmap, perp_u, area_exact, cap=(cap2 or cap),
                floor=floor)
            del bmap
            out["strips"] += strip_rows
            res = s437.assemble(slices, P, Q, area_exact, perp_u)
            _attach_kc(res, slices, geo, Q)
            res["area_ratio_hc"] = res["area_from_hc_total"] / area_exact
            res.update({"n_slices": len(slices), "n_bad_slices": 0,
                        "perp_sides_in_class": [sd for sd in ("L1", "L2", "H")
                                                if perp_u[sd] % 2 == par],
                        "n_prongs": len(reasons), "n_prong_bounces": nb,
                        "prongs_not_terminated": capped,
                        "n_open_at_cap1": n_open1, "prong_scan_aborted": aborted,
                        "min_cell_gap": diag.get("min_cell_gap"),
                        "n_strips_used": diag.get("n_strips", 0), "claim": diag})
            res["n_cyl_met_by_no_perp"] = sum(1 for c in res["cylinders"]
                                              if not c.get("met_perp_sides"))
            out["classes"][str(par)] = res
            continue
        # strips of this class, cheapest (fewest cells) first, perpendicular first
        cand = []
        for side in geo:
            for u in range(par, 4 * Q, 2):
                if s437.flux_weight(side, u, Q, geo) <= 1e-12:
                    continue
                cand.append((0 if u == perp_u[side] else 1,
                             len(set(bmap.get((side, u), []))), side, u))
        cand.sort()
        slices, used, min_gap = [], [], float("inf")
        # ⇒⇒ SKIP THE CELL PHASE WHEN A PRONG DID NOT TERMINATE.  A class with an open
        # prong is `PRONG_CAPPED` and is never scored, so materialising its strips buys
        # nothing -- and it is precisely the expensive case: if the direction is not
        # completely periodic then cells do not return either, so EVERY cell also runs to
        # the cap.  Two count-mode shards sat >8 min on one such centre before this.
        # ⚠ This makes a non-CP centre CHEAP rather than absent; the row is still stored,
        # still listed by `analyse()`, and still not scored ([NCYL-020]).
        if capped:
            cand = []
        for _, _, side, u in cand:
            reader = strip_cells_fast if fast else strip_cells
            cells, mg = reader(B, geo, Q, side, u, bmap.get((side, u), []), cap)
            min_gap = min(min_gap, mg)
            fl = s437.flux_weight(side, u, Q, geo)
            good = []
            for iv in cells:
                iv["flux"] = fl
                iv["area_contrib"] = iv["width"] * fl * iv["seg0"] if iv["ok"] else 0.0
                if iv["ok"]:
                    slices.append(iv)
                    good.append(iv)
            used.append((side, u))
            out["strips"].append({"side": side, "u": u, "parity": par,
                                  "is_perp": u == perp_u[side],
                                  "n_cells": len(good),
                                  "n_bad": len(cells) - len(good),
                                  "sum_width": sum(iv["width"] for iv in cells),
                                  "side_length": geo[side][0]})
            if strip_budget is not None and len(used) >= strip_budget:
                # ⚠ the KAC total cannot be used here -- it is an integral over the whole
                # 3-side cross-section and only reaches 1 with every strip present.  The
                # strip-local certificate is [NCYL-247]'s `Σ_c circ_c·h_c`.
                res0 = s437.assemble(slices, P, Q, area_exact, perp_u)
                _attach_kc(res0, slices, geo, Q)
                if abs(res0["area_from_hc_total"] / area_exact - 1.0) < 1e-9:
                    break
        res = s437.assemble(slices, P, Q, area_exact, perp_u)
        _attach_kc(res, slices, geo, Q)
        res["area_ratio_hc"] = res["area_from_hc_total"] / area_exact
        res.update({"n_slices": len(slices), "n_bad_slices": 0,
                    "perp_sides_in_class": [sd for sd in ("L1", "L2", "H")
                                            if perp_u[sd] % 2 == par],
                    "n_prongs": len(reasons), "n_prong_bounces": nb,
                    "prongs_not_terminated": capped,
                    "n_open_at_cap1": n_open1, "prong_scan_aborted": aborted,
                    "min_cell_gap": None if min_gap == float("inf") else min_gap,
                    "n_strips_used": len(used)})
        res["n_cyl_met_by_no_perp"] = sum(1 for c in res["cylinders"]
                                          if not c.get("met_perp_sides"))
        out["classes"][str(par)] = res
    out["seconds"] = round(time.time() - t0, 2)
    return out


# ---------------------------------------------------------------- controls

def h2_control(Qmax=40, verbose=True):
    """H2: the `O`/`A` prong counts against `s438_oddclass_cp.inward_count_exact`, which
    derives them independently.  A STOP condition."""
    import s438_oddclass_cp as cp
    bad = []
    for Q in range(3, Qmax + 1):
        for P in range(1, Q):
            if math.gcd(P, Q) != 1:
                continue
            for par in (0, 1):
                for v in ("O", "A"):
                    if len(prongs(v, P, Q, par)) != cp.inward_count_exact(P, Q, v, par):
                        bad.append((P, Q, par, v))
    if verbose:
        print(f"  H2 prong-window vs s438 inward_count_exact: {len(bad)} mismatches "
              f"over Q<={Qmax}", flush=True)
    return not bad


def h1_control(rows=None, verbose=True):
    """H1: reproduce `s437_oblique_cylinders.run_row` exactly -- per class `C_total` and
    Kac ratio, per strip `n_cells`, plus the boundary-position agreement (H3)."""
    rows = rows or s437.ROWS
    rep, ok = [], True
    for (P, Q) in rows:
        t0 = time.time()
        mine = run_row(P, Q)
        t_mine = time.time() - t0
        t0 = time.time()
        theirs = s437.run_row(P, Q)
        t_theirs = time.time() - t0
        mstr = {(s["side"], s["u"]): s["n_cells"] for s in mine["strips"]}
        tstr = {(s["side"], s["u"]): s["n_cells"] for s in theirs["strips"]}
        cell_bad = sorted(k for k in set(mstr) | set(tstr)
                          if mstr.get(k) != tstr.get(k))
        row = {"P": P, "Q": Q, "sec_exact": round(t_mine, 2),
               "sec_sampled": round(t_theirs, 2),
               "speedup": round(t_theirs / max(t_mine, 1e-9), 1),
               "n_strip_cell_mismatch": len(cell_bad),
               "strip_cell_mismatch": cell_bad[:5], "classes": {}}
        for par in ("0", "1"):
            a, b = mine["classes"][par], theirs["classes"][par]
            same = (a["C_total"] == b["C_total"]
                    and abs(a["total_area_ratio"] - b["total_area_ratio"]) < 1e-9)
            ok &= same and not cell_bad
            row["classes"][par] = {
                "C_exact": a["C_total"], "C_sampled": b["C_total"],
                "area_exact": round(a["total_area_ratio"], 12),
                "area_sampled": round(b["total_area_ratio"], 12),
                "min_cell_gap": a["min_cell_gap"],
                "n_prongs": a["n_prongs"], "n_prong_bounces": a["n_prong_bounces"],
                "prongs_not_terminated": a["prongs_not_terminated"],
                "match": same}
        rep.append(row)
        if verbose:
            for par in ("0", "1"):
                c = row["classes"][par]
                print(f"  H1 {P}/{Q:<3} par{par}: C {c['C_exact']} vs "
                      f"{c['C_sampled']}  area {c['area_exact']:.12f} vs "
                      f"{c['area_sampled']:.12f}  gap={c['min_cell_gap']:.3e}  "
                      f"prongs={c['n_prongs']}/{c['n_prong_bounces']}b  "
                      f"{'OK' if c['match'] else '** MISMATCH **'}", flush=True)
            print(f"     strips: {len(cell_bad)} n_cells mismatches "
                  f"{cell_bad[:5]} | {row['sec_exact']}s vs {row['sec_sampled']}s "
                  f"= {row['speedup']}x", flush=True)
    return ok, rep


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "control"
    if mode == "control":
        print("H2  prong windows (STOP condition)", flush=True)
        ok2 = h2_control()
        print("H1  exact vs sampled, s437's rows", flush=True)
        rows = None
        if len(sys.argv) > 2:
            rows = [tuple(int(v) for v in a.split("/")) for a in sys.argv[2:]]
        ok1, rep = h1_control(rows)
        with open("data/s439_exact_cells_control.json", "w") as f:
            json.dump({"h2_ok": ok2, "h1_ok": ok1, "rows": rep}, f, indent=1)
        print(f"H1={'PASS' if ok1 else 'FAIL'}  H2={'PASS' if ok2 else 'FAIL'}",
              flush=True)
        print("wrote data/s439_exact_cells_control.json", flush=True)
        sys.exit(0 if (ok1 and ok2) else 1)


if __name__ == "__main__":
    main()
