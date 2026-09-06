#!/usr/bin/env python3.13
"""s453_pole_module.py -- PRE-REGISTERED.  Queue item (alpha): close the PAIRED half of
`p = 0`.  What it actually delivers is a LEMMA that is class-blind and a criterion.

PRIOR ART: `rulings.py --grep` on 'module' -> [NCYL-194] (the FIBRE lattice
`Lambda = 2cot a * sum_m Z sin(mP pi/Q)` saturates `phi(2Q)/2`; F-leg strand, a rank
question about `nfan`, NOT about separatrices), [NCYL-173], [NCYL-280]; 'lattice' ->
[NB-071], [WFLOOR-037], nothing on the necklace; 'linearly independent' -> [NCYL-288]
ONLY (the Q-rank of the CELL-WIDTH span, `r_w = min(C_met, phi(2Q)/2)`; a different
object, and a Q-rank, not a Z-module); 'saddle connection' -> [NCYL-006], [NCYL-268],
[NCYL-259] -- all about FINDING/TRACING them, none about an arithmetic obstruction to
one existing; 'holonomy' -> [NCYL-018], [NCYL-099], [NCYL-228], nothing transverse;
'trace field' -> NOTHING; 'mod 2' -> nothing on this object; 'parity obstruction' ->
NOTHING.  The quotient path, the global coordinate and the fold rule are [NCYL-294];
the per-kite geometry and the tracer are [NCYL-292]; `p = 0` as a bijection is
[NCYL-287]; the kite-`m` localisation is [NCYL-297].
=> NEW here: the transverse-holonomy LEMMA below, its `Lambda / 2 Lambda` criterion, and
   `engine/zlattice.py`.  Everything else is prior art and is USED.
⚠ [NCYL-066] is the standing warning next door -- "an exhaustive negative is only as
  strong as its ansatz".  It does not bite: nothing here searches an ansatz.  The lemma
  is a proof, and the criterion is an EXACT module membership, not a bounded search.
⚠ [NCYL-296] bars routing this residual through CF depth / renormalization.  It does not
  bite either: the criterion is arithmetic in `Z[zeta_{4Q}]` and never reads a CF digit.

==========================================================================================
THE LEMMA (proved here; this is the session's content, and it is class-blind).

Work in [NCYL-294]'s GLOBAL coordinate on the `iota`-quotient path `J_0 ... J_m`
(`J_i = I_{m-i} = I_{m+i}`, `m = (Q-1)/2`).  Two facts from there, both already
established: the leaf's transverse coordinate `y` is CONSTANT between folds, and the
fold at edge `j` (between `J_j` and `J_{j+1}`) is the REFLECTION

        y  |->  s_j - y,        s_j := lo_j + lo_{j+1}   (tops shared)
                                    or hi_j + hi_{j+1}   (bottoms shared),

whose fixed point `s_j / 2` is exactly the interior pole on that edge ([NCYL-292] (iii)).
So after `n` folds at edges `a_1, ..., a_n`,

        y_n = (-1)^n y_0 + sum_{i=1..n} (-1)^{n-i} s_{a_i}.

Put `Lambda := sum_{a,b} Z (s_a - s_b)`, a Z-module of finite rank in `R`.  Suppose the
prong of pole `j` (so `y_0 = s_j/2`) reaches pole `j'` (so `y_n = s_{j'}/2`).

  * `n` EVEN: the coefficient sum of `sum (-1)^{n-i} s_{a_i}` is `0`, so that sum lies in
    `Lambda`, and `(s_{j'} - s_j)/2 in Lambda`.
  * `n` ODD: the coefficient sum is `1`, so the sum lies in `s_{a_n} + Lambda`, giving
    `(s_{j'} + s_j)/2 in s_{a_n} + Lambda`; and `s_{a_n} - s_j in Lambda`, so again
    `(s_{j'} - s_j)/2 in Lambda`.

  ⇒⇒ EITHER PARITY GIVES THE SAME CONDITION:   `s_{j'} - s_j  in  2*Lambda`.

  ⇒⇒ CONTRAPOSITIVE (THE CRITERION).  Give pole `j` the class
        `v_j := [s_j - s_0]  in  Lambda / 2 Lambda  =  (Z/2)^r`.
     IF THE `v_j` ARE PAIRWISE DISTINCT, NO INTERIOR POLE PRONG CAN REACH A *DIFFERENT*
     INTERIOR POLE.                                                                   [Q.E.D.]

  ⇒ THE CLOSED FORM, AND IT IS WHAT MAKES THE CRITERION PROVABLE.  With
    `sig_j = +1` bottoms-shared / `-1` tops-shared and `u_j` the (integer-combination)
    offset of node `j`,  `s_j = 2 u_j + L_j + sig_j L_{j+1}`, so mod `2 Z[zeta_{4Q}]`
    the offset dies and the sign is invisible:

        ⇒⇒   s_j  ==  L_j + L_{j+1}   (mod 2),

    and since `Lambda <= Z[zeta_{4Q}]`, the criterion follows from the strictly stronger,
    strictly arithmetic

        (*)   eta_{c_j} + eta_{c_{j+1}} + eta_{c_{j'}} + eta_{c_{j'+1}}  !=  0
              in  Z[zeta_{4Q}] / 2,        `eta_c = zeta^c - zeta^-c = 2i sin(c pi/2Q)`.

    No lattice, no dynamics, no CF: FOUR POINTS OF THE ROTATION ORBIT `{c_t}`.

  ⇒⇒ AND (*) IS A THEOREM AT PRIME `Q`.  Mod `2`, `-1 == +1` and `zeta^{2Q} = -1 == 1`,
    so `eta_c == x^c + x^-c`; reduce further mod `Phi_Q(x)`, where `x^Q == 1`, and the
    node exponents are `alpha_i = c_{m-i} mod Q = -2 P i mod Q` (`c_m = Q` paired / `0`
    lone, so `alpha_0 = 0` and `L_0 == 0` outright).  For `a, b in [0, m]`:
        `beta a == beta b`  =>  `a == b`   (`gcd(2P, Q) = 1`, `|a-b| <= m < Q`), and
        `beta a == -beta b` =>  `a + b == 0 mod Q`, impossible for `1 <= a+b <= 2m = Q-1`.
    So the eight exponents `+-alpha` over the four indices are DISTINCT mod `Q` (the only
    coincidence is `alpha_0 = -alpha_0`, which cancels `L_0` and is why `L_0 == 0`), and
    `s_{j'} - s_j` reduces to a sum of at most `8` distinct monomials.  A nonzero multiple
    of `Phi_Q` of degree `< Q` at PRIME `Q` is `Phi_Q` itself, which has all `Q`
    coefficients `1` -- so at `Q >= 11` prime a support of size `<= 8` cannot vanish.
    ⇒ (*) HOLDS, HENCE THE CRITERION HOLDS, HENCE NO POLE REACHES A DIFFERENT POLE.  ∎
    ⚠ SCOPE: `Q = 5, 7` are below the term count and are covered by H7/H8 directly.
    ⚠⚠ COMPOSITE `Q` IS *NOT* PROVED AND THE GAP IS REAL, NOT COSMETIC: `Phi_15 mod 2`
      has degree `8` and exactly `8` terms, which is the largest support this construction
      produces -- so at `Q = 15` the argument comes down to whether one specific `8`-set
      coincides with one other `8`-set.  It never does in range, and that is MEASURED.

  ⇒ WHY THE OBVIOUS SUFFICIENT CONDITION IS *NOT* THE ROUTE.  If `s_{j+1} - s_j` were
    Z-independent they would be a basis of `Lambda`, the `v_j` would be `0` and the
    standard basis vectors, and distinctness would be free.  They are not: with
        s_{j+1} - s_j  =  -sig_j L_j  +  sig_{j+1} L_{j+2}   (derived here),
    a two-term expression triangular in the highest index, independence would follow from
    Q-independence of `{L_0..L_m}` -- but `L_i = |cos(iP pi/Q)|` lies in `Q(zeta_{2Q})^+`,
    of degree only `phi(Q)/2`, and H3a measures `rank Lambda = min(m-1, phi(Q)/2)`
    EXACTLY.  So the generators are dependent on `748` of `1288` rows and the mod-2
    refinement is doing real work on all of them.

⚠⚠ WHAT THE LEMMA DOES NOT DO, AND IT MUST BE QUOTED WITH THIS.  It is silent on
   `j' = j`: `s_j - s_j = 0` is in `2 Lambda` always.  A prong returning to ITS OWN pole
   before running off an end is a `p != 0` configuration too -- in the full surface it is
   `R_k -> R_{+-k}` without crossing `Fix`, which is not `iota`-invariant (an invariant
   connection is reversed by `iota` and so contains a fixed point), hence a swapped pair.
   ⇒ THIS CLOSES THE `j' != j` HALF OF THE QUANTIFIER AND LEAVES THE `j' = j` HALF OPEN.
   Do NOT write it up as `p = 0` for the paired class.

HYPOTHESES.  ⚠ [OPS-041]: an arm that cannot fail is not evidence.
  H1  DERIVATION CHECK, NOT EVIDENCE.  `half(s_k) == s451.pole_y(...)` for every kite,
      and the palindrome `s_k == s_{Q-k}`.  Both are re-runs of the expression `qtrace`
      itself uses, so they can only catch an indexing slip -- which is the point.
  H2  ⇒⇒ THE CRITERION, AND IT CAN FAIL ON EVERY ROW: are the `v_j` pairwise distinct?
      Registered on odd `Q = 5..31` (s451's range), then run OUT OF SAMPLE to `Q = 61`
      and on the three STRESS rows `Q = 105, 165, 195` where the sufficient condition
      above is arithmetically unavailable.  BOTH CLASSES.
  H3/H3a  `rank Lambda` against `m-1` (the free-pass condition) and against
      `min(m-1, phi(Q)/2)` (the SATURATION LAW -- the analogue of [NCYL-288]'s
      `r_w = min(C_met, phi(2Q)/2)` on a different object).  H2b records the smallest
      Hamming distance between two classes: the TIGHTNESS of each row's pass.
  H7  THE CLOSED FORM: `s_j == L_j + L_{j+1} mod 2` (a derivation check of the `2u_j`
      step), and then (*) itself, scored per PAIR.
  H8  ⇒⇒ THE PROOF STEP, CHECKED RATHER THAN ASSERTED: the mod-`Phi_Q` exponent support
      is nonempty and has size `<= 8`.  This is the arm that says the prime-`Q` proof
      above is not hand-waving.
  H4  ⇒⇒ THE INSTRUMENT AUDIT, AND IT IS THE ARM THAT CAN KILL THE LEMMA: run the SAME
      walk-and-fold dynamics on SYNTHETIC paths (rational lengths, where `Lambda` has
      rank 1 and classes collide wholesale; and scrambled-cyclotomic lengths) and record
      every observed pole->pole hit.  EVERY ONE must have `v_j == v_{j'}`.  A single hit
      with distinct classes refutes the lemma.  ⚠ If no hit is ever observed the audit is
      VACUOUS and must be reported as such, not as a pass.
  H5  THE WALKER CROSS-CHECK: on real rows, my interval-based path walker must agree with
      `s451.reduced`'s `(i, d)` height model on end label and fold count.  Different
      state variable, different code; shares only `Chain`.
  H6  WHAT THE OPEN HALF COSTS: among synthetic paths that DO exhibit a pole->pole hit,
      how many are self-hits (`j' = j`)?  That is the fraction of the failure space the
      lemma cannot see, measured rather than guessed.

RESULTS (s453; `all --qmax=45`, `86.8 s`, `data/s453_pole_module.json`,
         `logs/s453_pole_module.log`).  Coverage: odd `Q = 5..45`, every coprime `P`,
         both `eps`, BOTH CLASSES, plus the three stress rows `Q = 105, 165, 195`
         -- `1288` class-rows, `1650636` pole pairs.
  H0  `zlattice` rank vs sympy `40/40`; coords round-trip `40/40`.
  H1  `98424/98424` (derivation check).
  H2  ⇒⇒ CLASSES PAIRWISE DISTINCT `644/644` PAIRED AND `644/644` LONE, INCLUDING ALL
      THREE STRESS ROWS.  On `748` of the `1288` rows `rank < m-1`, so the mod-2
      refinement -- not independence -- is what passes them.
  H2b MIN HAMMING DISTANCE IS `1` ON `1288/1288` ROWS: every row has a pair of poles one
      bit from colliding.  The birthday null over the rows is `10^-57`, and `260` rows sit
      below `P = 0.9` individually.  ⇒ THIS COULD HAVE FAILED, ROW BY ROW.
  H3a ⇒⇒ `rank Lambda = min(m-1, phi(Q)/2)` EXACTLY, `1288/1288` -- a saturation law.
      (H3, `rank = m-1`, holds on only `540/1288`, which is the same statement.)
  H4  ⇒⇒ INSTRUMENT AUDIT PASSES AND IT WAS NOT VACUOUS.  `400` rational paths (all with
      colliding classes) produced `262` cross-hits; `300` scrambled-cyclotomic paths
      (only `93` colliding) produced `74`.  ⇒ `0` LEMMA VIOLATIONS out of `336` observed
      pole->pole hits -- every one landed in a colliding-class path, as the lemma
      requires, and the `207` non-colliding scrambled paths produced none.
  H5  `638/638` against `s451.reduced` (end label AND fold count), `0` capped.
  H6  ⇒⇒ `0` OF THE `336` HITS WERE SELF-HITS.  The half of the failure space the lemma
      cannot see was never realised in the control family.  ⚠ That is a reason to expect
      the residual is small; it is NOT evidence that it is empty, and it is SYNTHETIC.
  H7  Derivation `1288/1288`; (*) holds on `1650636/1650636` pairs.
  H8  Support nonempty `1650636/1650636`; max `|support| = 8`; no prime-`Q` row reaches
      `|support| >= Q`.  ⇒ THE PRIME-`Q` PROOF ABOVE IS COMPLETE.
  ⚠⚠ THE CRUDE-TEST TRAP, WHICH I WALKED INTO ([OPS-228]): the first version of
    `crude_separated` asked "does the raw `Coord.v` carry an odd coordinate?" and read
    `0/1650636`, i.e. it looked like a clean NEGATIVE saying no cheap test exists.  It is
    an artefact: `Chain` stores `ell` as `2 * gen[c]`, so EVERY vector in this file is
    even by construction and the test was measuring the encoding.  Divide out the content
    first and it reads `1650636/1650636` -- the exact opposite, and the route to (*).

  python3.13 probes/s453_pole_module.py [all|arith|synth|walk] [--qmax=N]
"""
from __future__ import annotations

import json
import math
import random
import sys
import time
from fractions import Fraction

import s450_chain_maps as cm
import s451_quotient_path as qp
import zlattice

OUT = "data/s453_pole_module.json"
STRESS = (105, 165, 195)


# ---------------------------------------------------------------- the objects

class CoordOps:
    """Exact-comparison adapter for `s450_chain_maps.Coord` (uses the AUDITED `ch.cmp`)."""

    def __init__(self, ch):
        self.ch = ch

    def cmp(self, a, b):
        return self.ch.cmp(a, b)

    def eq(self, a, b):
        return a.eq(b)

    def half(self, a):
        return cm.half(a)

    def sub(self, a, b):
        return a - b


class FracOps:
    """The same protocol over `Fraction` -- used only by the synthetic control."""

    @staticmethod
    def cmp(a, b):
        return (a > b) - (a < b)

    @staticmethod
    def eq(a, b):
        return a == b

    @staticmethod
    def half(a):
        return a / 2

    @staticmethod
    def sub(a, b):
        return a - b


def fold_sums(ch, lo, hi):
    """`s_k` for every kite `k`, in the GLOBAL coordinate.

    Kite `k` joins interfaces `k-1` and `k`; the fold there is the reflection about
    `s_k / 2`.  This is verbatim the expression `s451.qtrace` builds in its loop, which
    is what makes H1 a derivation check rather than a test.
    """
    Q = ch.Q
    s = {}
    for k in range(Q):
        a, b = (k - 1) % Q, k
        s[k] = (hi[a] + hi[b]) if lo[a].eq(lo[b]) else (lo[a] + lo[b])
    return s


def path_view(ch, lo, hi, s):
    """The quotient PATH: nodes `i = 0..m` are interfaces `m-i`, edge `j` is kite `m-j`."""
    m = ch.m_idx
    Lo = [lo[m - i] for i in range(m + 1)]
    Hi = [hi[m - i] for i in range(m + 1)]
    S = [s[m - j] for j in range(m)]
    return Lo, Hi, S


def walk(Lo, Hi, S, j0, ops, cap=200_000):
    """The prong of the pole on edge `j0`, along the path.  Returns (end, folds).

    `end` is `('ESC', side)` / `('R', j)` / `('Z', i)` / `('CAP', None)`.  The pole sits
    on the LARGER of the two nodes and the prong heads AWAY from its edge ([NCYL-292]
    (iii)); the coordinate is constant between folds ([NCYL-294]).
    """
    n = len(Lo)
    y = ops.half(S[j0])
    la = ops.sub(Hi[j0], Lo[j0])
    lb = ops.sub(Hi[j0 + 1], Lo[j0 + 1])
    i = j0 if ops.cmp(la, lb) > 0 else j0 + 1
    d = -1 if i == j0 else +1
    folds = 0
    for _ in range(cap):
        i2 = i + d
        if i2 < 0 or i2 >= n:
            return ("ESC", 0 if i2 < 0 else 1), folds
        c0, c1 = ops.cmp(y, Lo[i2]), ops.cmp(y, Hi[i2])
        if c0 == 0 or c1 == 0:
            return ("Z", i2), folds
        if c0 > 0 and c1 < 0:
            i = i2
            continue
        j = min(i, i2)
        if ops.eq(y, ops.half(S[j])):
            return ("R", j), folds
        y = ops.sub(S[j], y)
        d = -d
        folds += 1
    return ("CAP", None), folds


def classes(S):
    """`(rank, [v_j], pairwise_distinct)` for a list of `Coord` fold sums `S`."""
    base = S[0]
    diffs = [[int(t) for t in (x - base).v] for x in S]
    lat = zlattice.Lattice([[a - b for a, b in zip(diffs[i + 1], diffs[i])]
                            for i in range(len(S) - 1)])
    vs = [lat.klass(dv, 2) for dv in diffs]
    ok = all(v is not None for v in vs) and len(set(vs)) == len(vs)
    return lat.rank, vs, ok


def crude_separated(S):
    """The CRUDE test: after dividing out the CONTENT, does `s_{j'} - s_j` carry an odd
    coordinate in the power basis?  If so it is outside `2 * (g * Z^d) >= 2 * Lambda`,
    so it separates the pair with no lattice at all -- and if it fired on every pair the
    criterion would collapse to a statement mod 2 in the ring, which is the shape a proof
    for all `Q` would want.  Returns (pairs, separated_by_crude).

    ⚠ THE CONTENT DIVISION IS NOT HYGIENE.  `Chain` stores `ell` as `2 * gen[c]`, so
    EVERY vector in this file is even by construction and the undivided test fires on
    `0/1650636` pairs -- a measurement of the encoding, not of the arithmetic.
    """
    diffs = [S[b].v - S[a].v for a in range(len(S)) for b in range(a + 1, len(S))]
    g = 0
    for w in diffs:
        for t in w:
            g = math.gcd(g, int(abs(t)))
    if g == 0:
        return len(diffs), 0
    return len(diffs), sum(bool(((abs(w) // g) & 1).any()) for w in diffs)


def mod2_form(ch, S):
    """H7 -- THE CLOSED FORM OF THE CRITERION.

    `s_j = 2 u_j + L_j + sig_j L_{j+1}` with `u_j` an INTEGER combination of the `L`'s,
    so mod `2` the offset dies and `sig_j = +-1 == 1`:

        s_j  ==  L_j + L_{j+1}   (mod 2 Z[zeta_{4Q}]),

    and the criterion `s_{j'} - s_j not in 2 Lambda` follows from the purely arithmetic

        eta_{c_j} + eta_{c_{j+1}} + eta_{c_{j'}} + eta_{c_{j'+1}}  !=  0  in  Z[zeta]/2,

    a statement about FOUR points of the rotation orbit `{c_t}` with no dynamics, no
    lattice and no CF in it.  Returns (derivation_ok, pairs, separated).
    """
    m = ch.m_idx
    lv = [(ch.ell[m - i].v // 2) & 1 for i in range(m + 1)]
    w = [(lv[j] ^ lv[j + 1]) for j in range(m)]
    deriv = all(((S[j].v // 2) & 1 == w[j]).all() for j in range(m))
    tot = sep = 0
    for a in range(m):
        for b in range(a + 1, m):
            tot += 1
            sep += bool((w[a] ^ w[b]).any())
    return deriv, tot, sep


def support_mod_Q(ch):
    """H8 -- THE STEP A PROOF NEEDS, checked rather than asserted.

    Reduce the closed form of H7 one step further, mod `Phi_Q(x)` (where `x^Q == 1`):
    `eta_c == x^c + x^{-c}` and the node exponents are `alpha_i = c_{m-i} mod Q`, i.e.
    `alpha_i = -2 P i mod Q` (`c_m = Q` paired / `0` lone, so `alpha_0 = 0` and
    `L_0 == 0`).  Then `s_{j'} - s_j` reduces to a sum of DISTINCT monomials, whose
    exponent SUPPORT is the symmetric difference of `{+-alpha}` over the four indices.

      * support NONEMPTY  <=>  the reduction is nonzero in `F_2[x]/(x^Q - 1)`;
      * `|support| < Q` and `Q` PRIME  =>  nonzero mod `Phi_Q` too, because the only
        multiple of `Phi_Q` of degree `< Q` other than `0` is `Phi_Q` itself, which has
        all `Q` coefficients `1`.  ⇒ AT PRIME `Q` THAT IS A PROOF OF THE CRITERION.

    Returns (pairs, nonempty, max_support).
    """
    Q, m = ch.Q, ch.m_idx
    al = [ch.c[m - i] % Q for i in range(m + 1)]
    tot = ne = mx = 0
    for a in range(m):
        for b in range(a + 1, m):
            sup = set()
            for i in (a, a + 1, b, b + 1):
                for e in (al[i] % Q, (-al[i]) % Q):
                    sup ^= {e}
            tot += 1
            ne += bool(sup)
            mx = max(mx, len(sup))
    return tot, ne, mx


def min_hamming(vs):
    """Smallest Hamming distance between two of the classes -- the TIGHTNESS of the pass
    (a `1` means one bit stood between this row and a collision)."""
    best = 10 ** 9
    for a in range(len(vs)):
        for b in range(a + 1, len(vs)):
            best = min(best, sum(x != y for x, y in zip(vs[a], vs[b])))
    return best


def birthday(k, r):
    """P[`k` uniform classes in `(Z/2)^r` are pairwise DISTINCT] -- the null model this
    result is quoted against.  ⚠ A HEURISTIC NULL, not a p-value on data: it prices how
    much room `Lambda / 2 Lambda` leaves, nothing more."""
    n = 1 << r
    if k > n:
        return 0.0
    p = 1.0
    for i in range(k):
        p *= (n - i) / n
    return p


def lattice_audit(rows, vec):
    """H0 -- `zlattice` against an independent rank, plus a coords round-trip."""
    from sympy import Matrix
    lat = zlattice.Lattice(rows)
    ok_rank = (lat.rank == Matrix([list(r) for r in rows]).rank())
    c = lat.coords(vec)
    if c is None:
        return ok_rank, False
    recon = [0] * lat.dim
    for ci, (_, brow) in zip(c, lat.basis):
        recon = [u + ci * v for u, v in zip(recon, brow)]
    return ok_rank, recon == [int(v) for v in vec]


def classes_frac(S):
    """The same over `Fraction` (rank-1 rational lattice)."""
    rows = zlattice.rational_rows([x - S[0] for x in S])
    lat = zlattice.Lattice([[rows[i + 1][0] - rows[i][0]] for i in range(len(S) - 1)])
    vs = [lat.klass(r, 2) for r in rows]
    return lat.rank, vs, all(v is not None for v in vs) and len(set(vs)) == len(vs)


# ---------------------------------------------------------------- hypotheses

def centres(qmax, extra=()):
    for Q in list(range(5, qmax + 1, 2)) + list(extra):
        for P in range(1, Q):
            if math.gcd(P, Q) == 1:
                for eps in (0, 1):
                    yield P, Q, eps


def arith(qmax):
    res = {"h1_ok": 0, "h1_bad": 0, "h1_examples": [],
           "h2": {}, "h3_full": [0, 0], "h3_examples": [], "rows": 0,
           "by_Q": {}, "rank_by_Q": {}, "h0": [0, 0, 0],
           "worst_birthday": [1.0, None], "tightest": [1e9, None],
           "at_risk": 0, "logprod": 0.0, "h3a": [0, 0], "h3a_ex": [],
           "hamming": {}, "crude": [0, 0], "crude_ex": [],
           "h7": [0, 0], "h7_deriv": [0, 0], "h8": [0, 0, 0], "h8_bad": []}
    for cls in ("paired", "lone"):
        res["h2"][cls] = {"rows": 0, "distinct": 0, "fail": [],
                          "rank_lt_m1": 0, "rank_lt_m1_distinct": 0}
    for P, Q, eps in centres(qmax, STRESS):
        ch, lo, hi, closes = qp.build(P, Q, eps)
        m, cls = ch.m_idx, ("lone" if ch.lone else "paired")
        s = fold_sums(ch, lo, hi)
        res["rows"] += 1
        # H1 -- derivation check
        for k in range(1, Q):
            if ch.tie[k]:
                continue
            good = cm.half(s[k]).eq(qp.pole_y(ch, lo, hi, k)[0]) and s[k].eq(s[Q - k])
            res["h1_ok" if good else "h1_bad"] += 1
            if not good and len(res["h1_examples"]) < 5:
                res["h1_examples"].append([P, Q, eps, k])
        if not closes:                       # H2 of s451; the coordinate must exist
            res["h2"][cls]["fail"].append([P, Q, eps, "NO-CLOSURE"])
            continue
        # H2 / H3 -- the criterion, on the m poles of the quotient path
        _, _, S = path_view(ch, lo, hi, s)
        rank, vs, ok = classes(S)
        # H3a -- the saturation law, and H2b -- how close the pass came
        res["h3a"][0 if rank == min(m - 1, _phi(Q) // 2) else 1] += 1
        if rank != min(m - 1, _phi(Q) // 2) and len(res["h3a_ex"]) < 8:
            res["h3a_ex"].append([P, Q, eps, rank, m - 1, _phi(Q) // 2])
        if ok:
            mh = min_hamming(vs)
            res["hamming"][str(mh)] = res["hamming"].get(str(mh), 0) + 1
        t8, n8, x8 = support_mod_Q(ch)
        res["h8"][0] += t8
        res["h8"][1] += n8
        res["h8"][2] = max(res["h8"][2], x8)
        if _isprime(Q) and x8 >= Q:
            res["h8_bad"].append([P, Q, eps, x8])
        d7, t7, s7 = mod2_form(ch, S)
        res["h7"][0] += t7
        res["h7"][1] += s7
        res["h7_deriv"][0 if d7 else 1] += 1
        ct, cs = crude_separated(S)
        res["crude"][0] += ct
        res["crude"][1] += cs
        if cs < ct and len(res["crude_ex"]) < 8:
            res["crude_ex"].append([P, Q, eps, ct - cs, ct])
        d = res["h2"][cls]
        d["rows"] += 1
        d["distinct"] += ok
        if not ok and len(d["fail"]) < 8:
            d["fail"].append([P, Q, eps, rank, m])
        res["h3_full"][0 if rank == m - 1 else 1] += 1
        if rank != m - 1 and len(res["h3_examples"]) < 8:
            res["h3_examples"].append([P, Q, eps, rank, m - 1, _phi(Q)])
        if rank < m - 1:
            d["rank_lt_m1"] += 1
            d["rank_lt_m1_distinct"] += ok
        # how much ROOM the criterion had: `m` classes in `2**rank` slots
        pb = birthday(m, rank)
        if pb < res["worst_birthday"][0]:
            res["worst_birthday"] = [pb, [P, Q, eps, m, rank]]
        slack = (1 << rank) / max(m, 1)
        if slack < res["tightest"][0]:
            res["tightest"] = [slack, [P, Q, eps, m, rank]]
        # H0 -- audit the lattice code itself on this row's real generators
        if res["h0"][0] < 40:
            rows = [[int(a) - int(b) for a, b in zip((S[i + 1] - S[0]).v,
                                                     (S[i] - S[0]).v)]
                    for i in range(len(S) - 1)]
            r_ok, c_ok = lattice_audit(rows, [int(t) for t in (S[-1] - S[0]).v])
            res["h0"][0] += 1
            res["h0"][1] += r_ok
            res["h0"][2] += c_ok
        res["at_risk"] += (pb < 0.9)
        res["logprod"] += math.log(max(pb, 1e-300))
        b = res["by_Q"].setdefault(str(Q), [0, 0])
        b[0] += 1
        b[1] += ok
        rk = res["rank_by_Q"].setdefault(str(Q), [m, 10 ** 9, 0])
        rk[1] = min(rk[1], rank)
        rk[2] = max(rk[2], rank)
    return res


def _isprime(n):
    if n < 2:
        return False
    d = 2
    while d * d <= n:
        if n % d == 0:
            return False
        d += 1
    return True


def _phi(n):
    r, p = n, 2
    out = n
    while p * p <= r:
        if r % p == 0:
            while r % p == 0:
                r //= p
            out -= out // p
        p += 1
    if r > 1:
        out -= out // r
    return out


def walk_check(qmax):
    """H5 -- the interval walker against `s451.reduced`'s `(i, d)` height model."""
    ok = bad = skipped = 0
    ex = []
    for P, Q, eps in centres(qmax):
        ch, lo, hi, _ = qp.build(P, Q, eps)
        if ch.lone:
            continue
        m = ch.m_idx
        s = fold_sums(ch, lo, hi)
        Lo, Hi, S = path_view(ch, lo, hi, s)
        L = [ch.ell[m - i] for i in range(m + 1)]
        align = [(((i + 1) * P / Q + 0.5) % 1.0) < P / Q - 1e-12 for i in range(m)]
        ops = CoordOps(ch)
        for j in range(m):
            big = j if ch.cmp(L[j], L[j + 1]) > 0 else j + 1
            sred = (L[big] - L[2 * j + 1 - big]) if align[j] else (L[j] + L[j + 1])
            r_end, r_folds, _ = qp.reduced(ch, m, L, align, big, cm.half(sred),
                                           -1 if big == j else +1)
            w_end, w_folds = walk(Lo, Hi, S, j, ops)
            if r_end == "CAP" or w_end[0] == "CAP":
                skipped += 1
                continue
            same = ((r_end == "ESC") == (w_end[0] == "ESC")) and r_folds == w_folds
            if same:
                ok += 1
            else:
                bad += 1
                if len(ex) < 5:
                    ex.append([P, Q, eps, j, r_end, r_folds, str(w_end), w_folds])
    return {"ok": ok, "bad": bad, "skipped": skipped, "examples": ex}


def synth(trials, seed=20260901):
    """H4/H6 -- the instrument audit.  Synthetic paths; every pole->pole hit must have
    colliding classes.  Family A: small rational lengths (rank-1 `Lambda`, so classes
    collide wholesale and hits are expected).  Family B: a random permutation of a real
    row's cyclotomic lengths with a random alignment word."""
    rng = random.Random(seed)
    out = {}
    for fam in ("rational", "scrambled"):
        rec = {"paths": 0, "hits": 0, "self_hits": 0, "cross_hits": 0,
               "violations": [], "esc": 0, "z": 0, "cap": 0, "collide_paths": 0}
        for _ in range(trials):
            if fam == "rational":
                n = rng.randint(4, 9)
                L = [Fraction(rng.randint(1, 12)) for _ in range(n)]
                ops = FracOps
                ZERO = Fraction(0)
            else:
                Q = rng.choice([11, 13, 17, 19, 23])
                P = rng.choice([p for p in range(1, Q) if math.gcd(p, Q) == 1])
                ch = cm.Chain(P, Q, 0)
                n = ch.m_idx + 1
                L = [ch.ell[t] for t in rng.sample(range(Q), n)]
                ops = CoordOps(ch)
                ZERO = ch.zero
            sig = [rng.choice((+1, -1)) for _ in range(n - 1)]
            Lo, Hi, S = [ZERO], [L[0]], []
            for j in range(n - 1):
                nlo = Lo[j] if sig[j] > 0 else (Lo[j] + L[j]) - L[j + 1]
                Lo.append(nlo)
                Hi.append(nlo + L[j + 1])
            for j in range(n - 1):
                S.append((Hi[j] + Hi[j + 1]) if sig[j] > 0 else (Lo[j] + Lo[j + 1]))
            if any(ops.cmp(Lo[i], Hi[i]) >= 0 for i in range(n)):
                continue                              # a zero-length node: skip
            rank, vs, distinct = (classes_frac(S) if fam == "rational" else classes(S))
            if vs is None or any(v is None for v in vs):
                continue
            rec["paths"] += 1
            rec["collide_paths"] += (not distinct)
            for j in range(n - 1):
                end, _ = walk(Lo, Hi, S, j, ops, cap=40_000)
                if end[0] == "ESC":
                    rec["esc"] += 1
                elif end[0] == "Z":
                    rec["z"] += 1
                elif end[0] == "CAP":
                    rec["cap"] += 1
                else:                                              # ('R', j2)
                    j2 = end[1]
                    rec["hits"] += 1
                    rec["self_hits" if j2 == j else "cross_hits"] += 1
                    if j2 != j and vs[j] != vs[j2] and len(rec["violations"]) < 5:
                        rec["violations"].append([fam, n, j, j2, str(vs[j]), str(vs[j2])])
        out[fam] = rec
    return out


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "all"
    qmax = 45
    for a in sys.argv[2:]:
        if a.startswith("--qmax="):
            qmax = int(a.split("=")[1])
    t0 = time.time()
    out = {"mode": mode, "qmax": qmax, "stress": list(STRESS)}
    if mode in ("all", "arith"):
        r = arith(qmax)
        out["arith"] = r
        print(f"H1 derivation check  {r['h1_ok']}/{r['h1_ok']+r['h1_bad']}"
              f"   {r['h1_examples']}   (CANNOT FAIL unless indexing is wrong)")
        for cls in ("paired", "lone"):
            d = r["h2"][cls]
            print(f"H2 {cls:6s} classes pairwise DISTINCT {d['distinct']}/{d['rows']}"
                  f"   fails={d['fail']}")
            print(f"   of which rank < m-1 (mod-2 doing real work): "
                  f"{d['rank_lt_m1_distinct']}/{d['rank_lt_m1']}")
        print(f"H0 zlattice audit  rank vs sympy {r['h0'][1]}/{r['h0'][0]}   "
              f"coords round-trip {r['h0'][2]}/{r['h0'][0]}")
        print(f"H3 rank Lambda == m-1 (the SUFFICIENT condition)  "
              f"{r['h3_full'][0]}/{sum(r['h3_full'])}   [P,Q,eps,rank,m-1,phi]"
              f" {r['h3_examples']}")
        print(f"   ROOM: worst birthday-null P[distinct] = {r['worst_birthday'][0]:.3g}"
              f" at {r['worst_birthday'][1]}   tightest 2^rank/m ="
              f" {r['tightest'][0]:.3g} at {r['tightest'][1]}")
        print(f"   rows where distinctness was AT RISK under that null (P<0.9): "
              f"{r['at_risk']}/{r['rows']};  log10 prod P[all distinct] = "
              f"{r['logprod']/math.log(10):.1f}")
        print(f"H3a rank Lambda == min(m-1, phi(Q)/2)  {r['h3a'][0]}/{sum(r['h3a'])}"
              f"   {r['h3a_ex']}")
        print(f"H2b min Hamming distance between classes, distribution: "
              f"{dict(sorted(r['hamming'].items(), key=lambda t: int(t[0])))}")
        print(f"H7 CLOSED FORM  s_j == L_j + L_{{j+1}} mod 2: derivation "
              f"{r['h7_deriv'][0]}/{sum(r['h7_deriv'])}   |   four-term sums NONZERO "
              f"in Z[zeta]/2: {r['h7'][1]}/{r['h7'][0]} pairs")
        print(f"H8 mod Phi_Q support NONEMPTY {r['h8'][1]}/{r['h8'][0]} pairs;  max "
              f"|support| = {r['h8'][2]}  (proof needs < Q, i.e. Q >= 11 at prime Q); "
              f" prime-Q rows where it fails: {r['h8_bad'] or 'NONE'}")
        print(f"⇒ CRUDE TEST (raw vector has an ODD coordinate, no lattice needed): "
              f"{r['crude'][1]}/{r['crude'][0]} pairs   {r['crude_ex']}")
        for Q in STRESS:
            b, rk = r["by_Q"].get(str(Q)), r["rank_by_Q"].get(str(Q))
            if b:
                print(f"   STRESS Q={Q}: distinct {b[1]}/{b[0]}   m={rk[0]}"
                      f"  rank Lambda in [{rk[1]},{rk[2]}]  (phi={_phi(Q)})")
    if mode in ("all", "walk"):
        out["h5"] = walk_check(min(qmax, 21))
        g = out["h5"]
        print(f"H5 walker vs s451.reduced  {g['ok']}/{g['ok']+g['bad']}"
              f"  skipped(cap)={g['skipped']}  {g['examples']}")
    if mode in ("all", "synth"):
        out["h4"] = synth(400)
        for fam, rec in out["h4"].items():
            print(f"H4 {fam:10s} paths={rec['paths']} (colliding {rec['collide_paths']})"
                  f"  esc={rec['esc']} Z={rec['z']} cap={rec['cap']}"
                  f"  HITS={rec['hits']} (self {rec['self_hits']}, "
                  f"cross {rec['cross_hits']})")
            print(f"   ⇒ LEMMA VIOLATIONS (cross-hit with DISTINCT classes): "
                  f"{rec['violations'] or 'NONE'}"
                  + ("   ⚠ AUDIT VACUOUS: no cross-hit observed"
                     if rec["cross_hits"] == 0 else ""))
    out["elapsed_s"] = round(time.time() - t0, 1)
    with open(OUT, "w") as f:
        json.dump(out, f, indent=1, default=str)
    print(f"-> {OUT}   [{out['elapsed_s']}s]")


if __name__ == "__main__":
    main()
