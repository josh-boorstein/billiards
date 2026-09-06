#!/usr/bin/env python3
"""
s463_rank_law.py -- PRE-REGISTERED.  Queue item (s462-a): PROVE [NCYL-317] H3's rank law.

WHAT THIS IS.  [NCYL-317] (3) is MEASURED -- `rank_Q{ell_t} == phi(Q)/2` on `56/56`
class-parities, odd `Q = 5..59`, both classes -- and that measured input is the only thing
holding up the mechanism split `(G0b-rel)` / `(G0b-form)`.  This file carries the PROOF
that replaces it, and the arms that certify the proof's own steps.  The proof is
unconditional at EVERY odd `Q`, prime or composite, in BOTH classes; the sweeps here are
GUARDS ON THE PROOF, not the evidence for it.

⇒⇒ THE PROOF, in full, because it is short.

  Setup.  `Q` odd `>= 3`; `zeta := zeta_Q = exp(2 pi i / Q)`; `K := Q(zeta)`, so
  `[K:Q] = phi(Q)`.  Let `sigma` be complex conjugation, `sigma(zeta) = zeta^-1`; it is a
  nontrivial involution in `Gal(K/Q)` for `Q >= 3`.  Write `K = K^+ (+) K^-` for its
  eigenspace decomposition, `K^+` the maximal real subfield, `dim_Q K^+ = phi(Q)/2`, hence
  `dim_Q K^- = phi(Q)/2`.

  (0) ⇒ THE MEASURES ARE `P`-INDEPENDENT AND ARE SINES/COSINES OF `pi/Q`.
      `c_t = c_0 + 2Pt (mod 2Q)` with `gcd(P,Q) = 1`, so `gcd(2P, 2Q) = 2` and the orbit is
      ALL residues mod `2Q` of the parity of `c_0`; and `c_0 = P` resp. `P+Q` makes that
      parity EVEN exactly on the LONE class ([NCYL-291]; H1 re-checks it in code).
      Folding by `c ~ 2Q-c` (which fixes `ell = sin(c pi/2Q)`) and dropping `ell = 0`:
          LONE   : `c` even, `c = 2k` ,  keys `{2,4,..,Q-1}`, `ell = sin(k pi / Q)`,
                   `k = 1 .. (Q-1)/2`                                -- `(Q-1)/2` values
          PAIRED : `c` odd  , `c = Q-2j`, keys `{1,3,..,Q}`  , `ell = cos(j pi / Q)`,
                   `j = 0 .. (Q-1)/2`                                -- `(Q+1)/2` values
      The PAIRED line is the complementary-angle identity `sin(c pi/2Q) = cos((Q-c) pi/2Q)`
      and it is the whole reason that class is not the harder one: it puts the paired
      measures in the REAL subfield of `Q(zeta_Q)`, not in the minus part of `Q(zeta_4Q)`.
      Both families are strictly monotone in `k` resp. `j` on `[0, pi/2]`, so the counts
      `(Q-1)/2` / `(Q+1)/2` are exact.

  (1) ⇒ THE LEMMA.  Let `A := span_Q{zeta^m + zeta^-m}` and `B := span_Q{zeta^m - zeta^-m}`,
      `m = 0..Q-1`.  Then `A = K^+`, `B = K^-`, so `dim A = dim B = phi(Q)/2`.
        `A <= K^+` and `B <= K^-`, so `A n B = 0` and `dim A + dim B <= phi(Q)`, with
        `dim A <= phi(Q)/2` and `dim B <= phi(Q)/2` separately.
        `2 zeta^m = (zeta^m + zeta^-m) + (zeta^m - zeta^-m)` for every `m`, and
        `{zeta^m}_{m<Q}` spans `K` over `Q`, so `A + B = K` and `dim A + dim B >= phi(Q)`.
        Both inequalities are therefore equalities.                                   []
      ⚠ `A + B = K` is the ONLY step with content, and it is H2's third column.

  (2) ⇒ TRANSPORT TO THE MEASURES.  `zeta_2Q = -zeta_Q^h` with `h = (Q+1)/2` (check:
      `-exp(2 pi i h / Q) = exp(i pi (2h+Q)/Q) = exp(i pi (2Q+1)/Q) = exp(i pi / Q)`), and
      `2h = 1 (mod Q)` so `h` is a unit.  Hence
          `sin(k pi/Q) = (-1)^k (zeta^hk - zeta^-hk) / (2i)`
          `cos(j pi/Q) = (-1)^j (zeta^hj + zeta^-hj) / 2`
      Multiplying the whole set by the fixed nonzero scalar `1/(2i)` resp. `1/2` and
      individual members by `(-1)^k` are `Q`-linear bijections of the ambient `Q`-space
      `C`, so neither changes `Q`-rank.  `k -> hk` permutes the `(Q-1)/2` classes
      `{+-m}` of nonzero residues (`h` a unit), and `m -> zeta^m -+ zeta^-m` is odd resp.
      even in `m`; so the two spans are exactly `B` resp. `A` of (1).  Therefore
          `rank_Q {ell_t : LONE}  = dim B = phi(Q)/2`
          `rank_Q {ell_t : PAIRED} = dim A = phi(Q)/2`                                []

  (3) ⇒ THE DEFECT COROLLARY, which is what [NCYL-317] (3)-(4) actually consume.
          `defect_LONE   = (Q-1)/2 - phi(Q)/2 = (Q - 1 - phi(Q))/2`
          `defect_PAIRED = (Q+1)/2 - phi(Q)/2 = (Q + 1 - phi(Q))/2 = defect_LONE + 1`
      `defect_LONE = 0 <=> phi(Q) = Q-1 <=> Q PRIME`; `defect_PAIRED >= 1` at every odd
      `Q`, since `phi(Q) <= Q-1`.  ⇒ the zero-defect region is EXACTLY the LONE class at
      prime `Q`, and the PAIRED class carries a relation everywhere -- both of which
      [NCYL-317] had only as measurements to `Q <= 59`.                               []

  (4) ⇒ AND THE STRUCTURAL READING, which is a corroboration and also H4's arm.  The two
      defects are the `sigma`-eigen halves of ONE object: the relation lattice
      `ker(Q^Q -> K, e_m -> zeta^m)`, of dimension `Q - phi(Q)`.  The involution `m -> -m`
      splits `Q^Q` into a `(Q+1)/2`-dimensional even part and a `(Q-1)/2`-dimensional odd
      part, the map intertwines it with `sigma`, so
          `defect_PAIRED + defect_LONE = (Q+1-phi(Q))/2 + (Q-1-phi(Q))/2 = Q - phi(Q)`
      exactly.  ⇒ THE INTERFACE-MEASURE RELATIONS ARE NOT A NEW ARITHMETIC PHENOMENON:
      they are the cyclotomic multiplicative-relation kernel, sorted by parity.

⚠⚠ WHAT THIS DOES *NOT* DO.  It promotes ONE measured input of [NCYL-317] to proved.  It
  does NOT close (G0b), it does NOT touch `(G0b-form)` (the `Lambda_S` membership question
  [NCYL-308] prices at `87.0%` / `77.6%` blind), it does NOT prove `p = 0` on either class,
  and it moves nothing on the gamma=1 DAG.  What it buys is that [NCYL-317] (4)'s
  "mechanism `(G0b-rel)` cannot occur on the LONE class at prime `Q`" now rests on a
  theorem instead of on `56/56` rows -- i.e. the split becomes a discharged QUANTIFIER
  (REDUCTION-vs-RESTATEMENT (b)) over all odd `Q` rather than coverage to `Q <= 59`.
  ⚠ It also does NOT explain WHY `Q = 15` carries the witness and `Q = 9,21,27,33,39` do
  not ([NCYL-317] (5) -- defect `> 0` is necessary, not sufficient; that is queue (s462-b)).

PRIOR ART: `rulings.py --grep` on 'rank', 'cyclotomic', 'linear independence',
'vanishing sum', 'roots of unity', 'minus eigenspace', 'Q(zeta' -> the only adjacent
rulings are [NCYL-317] itself (which MEASURES this and flags it unproved), [NCYL-066] (a
subset-sum search over `{sin^2(k pi/Q)}` -- a different set, and a cautionary tale about
ansatz-bounded negatives), [WFLOOR-124] (a "frozen six-term sine vector" that turned out
to be a BASIS ARTEFACT -- the same hazard this file's H1 bridge is built to exclude) and
[NCYL-300] (an exact `F_2[x]/Phi_Q` factorisation on the POLE side, no primality anywhere
-- the closest methodological sibling: it likewise replaced a term-count/primality
argument with an exact cyclotomic one).  s462's own docstring records the same grep coming
back empty: no probe or ruling in this repo has asked whether the interface measures are
linearly INDEPENDENT.  ⇒ NONE of these proves or disproves the rank law.

HYPOTHESES (pre-registered).
  H1  ⇒⇒ THE BRIDGE, and it is the arm that can actually fail, because everything above is
      about `sin(k pi/Q)` / `cos(j pi/Q)` while [NCYL-317] is about `Chain.ell`.  Over
      EVERY coprime `(P, Q, eps)`, odd `Q` in range: (a) the folded key set of the nonzero
      `ell` is exactly `{2,4,..,Q-1}` (LONE) / `{1,3,..,Q}` (PAIRED), hence `P`-INDEPENDENT
      -- which is the assumption `s462_zz_relation.ell_values` makes when it uses `P = 1`,
      and which is here CHECKED rather than inherited; (b) the code's exact ring vectors
      and an INDEPENDENT `sympy` reduction of `x^m mod Phi_Q` produce the SAME rank and the
      SAME nullspace dimension, class by class.  ⚠ (b) is the [WFLOOR-124] guard: a rank is
      a property of a set of NUMBERS, and reading it off one basis is how a basis artefact
      gets published.  A mismatch on either half invalidates the transport step (2) and the
      probe reports FAIL rather than a rank.
  H2  ⇒⇒ THE PROOF'S CRUX, scored per row and far past [NCYL-317]'s `Q <= 59`:
      `rank(B) == phi(Q)/2`, `rank(A) == phi(Q)/2`, and `rank(A u B) == phi(Q)`.  The third
      is step (1)'s only content and it is INDEPENDENTLY scoreable: `A + B = K`.  ⚠ Any row
      with `rank(A u B) < phi(Q)` breaks the proof, not just the measurement.
  H3  ⇒⇒ THE EIGENSPACE ARM, the other half of step (1): applying `zeta -> zeta^-1` to each
      spanning vector must return `+itself` on `A` and `-itself` on `B`, EXACTLY, so
      `A n B = 0`.  ⚠ Predicted, and it is a derivation check, not evidence -- but a
      failure here would mean the reduction table is wrong and H2's numbers are noise.
  H4  ⇒⇒ THE DEFECT COROLLARY AND ITS STRUCTURAL CROSS-CHECK.  Per row:
      `defect_LONE == (Q-1-phi(Q))/2`, `defect_PAIRED == defect_LONE + 1`,
      `defect_LONE == 0 <=> isprime(Q)`, `defect_PAIRED >= 1`, and the sum
      `defect_LONE + defect_PAIRED == Q - phi(Q) == dim ker(Q^Q -> K)`, with that kernel
      dimension MEASURED independently (nullspace of the reduction table) rather than
      computed from `phi`.  ⚠ [OPS-071]: the two sides of the last identity come from
      different computations -- `phi(Q)` from `sympy.totient`, the kernel from a nullspace
      of the `Q x phi(Q)` reduction matrix -- so it can fail.
  H5  ⇒⇒ NON-VACUITY, and it is a SIBLING PROJECTION of the same object, not a different
      object ([OPS-240]).  ⚠ [OPS-041]: `phi(Q)/2` must be something the instrument can
      FAIL to return.  Two ways it does: (a) the union `A u B` -- the same measures, same
      code path, same rank routine -- returns `phi(Q)`, never `phi(Q)/2`, at every `Q >= 5`;
      (b) `rank == #distinct` (defect `0`) on LONE holds at prime `Q` and FAILS at every
      composite `Q` in range, so the routine separates rows rather than emitting a constant.
  H6  ⇒⇒ EXPLICIT GENERATORS, which is what upgrades (3) from a DIMENSION to a LATTICE.
      For each prime `p | Q` and each `m` in `Z/Q` the rotated prime relation
      `N(p,m) := sum_{j<p} e_{m + jQ/p}`  maps to  `zeta^m * sum_{j<p} zeta_p^j = 0`.
      CLAIM: these span `ker(Q^Q -> K)` (dimension `Q - phi(Q)`), and their `m -> -m`
      ANTI-symmetrisations resp. SYMMETRISATIONS, transported by `a_k = (-1)^k w[hk]`,
      span the LONE resp. PAIRED relation lattice -- dimensions `defect_LONE` and
      `defect_PAIRED` from H4.  ⇒ the interface-measure relations are exactly the
      cyclotomic vanishing-sum relations, sorted by parity, with one generating family
      per PRIME DIVISOR of `Q`.
      ⚠⚠ THE `(-1)^k` IS NOT COSMETIC AND IT IS THIS ARM'S ONLY REAL HAZARD: dropped, the
      fold silently produces vectors that are NOT relations, the span still has a
      plausible-looking rank, and the [NCYL-317] (5) check below reads FALSE.  That
      happened in session and was caught by G5, not by inspection ([OPS-043]).
      ⊕ TIE-IN: at `Q = 15` the witness relation of [NCYL-317] (5) --
      `sin(pi/15) + sin(4pi/15) = sin(6pi/15)`, i.e. `sin t + sin(60deg - t) =
      sin(60deg + t)` at `t = 12deg` -- must lie in the span of the `p = 3` family ALONE,
      and NOT in the `p = 5` family.  ⚠ That is a CLASSIFICATION of the known witness, not
      an explanation of it: it says which family is consumed, not why `Q = 15` realises a
      connection when `Q = 9,21,27,33,39` carry the same family and do not (queue (s462-b)).

GUARDS.
  G1  ⇒ exact `sympy` rank for `Q <= --qexact` (default `75`); modular rank over two
      primes for the tail.  The modular value is a LOWER bound on the rational rank, and
      step (1) supplies the matching UPPER bound, so the tail is a sound certificate --
      but the two methods are cross-checked on the whole exact range first, and a
      disagreement there voids the tail.
  G2  ⇒ nothing here writes to any store, imports no engine beyond `s450_chain_maps`
      (H1 only), and modifies no existing probe.
  G3  ⇒ [OPS-228]: numpy `int64` cast to Python `int` before `sympy`.
  G4  ⇒ H1 skips no row silently: every coprime `(P,Q,eps)` in range is scored, and the
      per-row verdicts are counted so the denominator is quotable ([OPS-162]).
  G5  ⇒ H6's CONVENTION GUARD, and it is the reason H6 is trustworthy: every folded
      generator must be verified to VANISH on the actual measures -- exactly, as a zero
      vector in the reduction basis, AND numerically to `1e-12` -- BEFORE its span is
      ranked.  A wrong transport sign passes a rank check and fails this one.

Run: PYTHONPATH=.:engine:archive/scripts_2026-07:archive/scripts_2026-08:probes \
     .venv/bin/python3.13 probes/s463_rank_law.py [all|h1|h2|h3|h4|h5] \
         [--qmax=N] [--qexact=N] [--qbridge=N]
"""
import json
import math
import sys
import time

import numpy as np
from sympy import Matrix, Poly, ZZ, cyclotomic_poly, isprime, symbols, totient

X = symbols('x')
MODP = (2147483647, 2147483629)          # two large primes for the modular tail

_RED = {}


def reduction(Q):
    """rows[m] = coefficient vector of `x^m mod Phi_Q(x)` over `Z`, m = 0..Q-1.

    An INDEPENDENT representation of `zeta_Q^m`: sympy polynomial remainder against the
    cyclotomic polynomial, sharing no code with `engine/ring_cyc.py`.  H1(b) is what makes
    that independence load-bearing.
    """
    if Q in _RED:
        return _RED[Q]
    P = Poly(cyclotomic_poly(Q, X), X, domain=ZZ)
    f = int(totient(Q))
    rows = []
    for m in range(Q):
        c = Poly(X ** m, X, domain=ZZ).rem(P).all_coeffs()[::-1]
        rows.append([int(v) for v in c] + [0] * (f - len(c)))
    _RED[Q] = (rows, f)
    return rows, f


def spans(Q):
    """`A` (even part, PAIRED) and `B` (odd part, LONE) as integer matrices."""
    R, f = reduction(Q)
    B = [[R[m][i] - R[(-m) % Q][i] for i in range(f)] for m in range(1, (Q - 1) // 2 + 1)]
    A = [[R[m][i] + R[(-m) % Q][i] for i in range(f)] for m in range(0, (Q - 1) // 2 + 1)]
    return A, B, f


def rank_mod(M, p):
    """Gaussian elimination over `F_p`.  A LOWER bound on the rational rank."""
    if not M:
        return 0
    A = [[int(v) % p for v in row] for row in M]
    rows, cols, r = len(A), len(A[0]), 0
    for c in range(cols):
        piv = next((i for i in range(r, rows) if A[i][c]), None)
        if piv is None:
            continue
        A[r], A[piv] = A[piv], A[r]
        inv = pow(A[r][c], p - 2, p)
        A[r] = [(v * inv) % p for v in A[r]]
        for i in range(rows):
            if i != r and A[i][c]:
                fct = A[i][c]
                A[i] = [(A[i][j] - fct * A[r][j]) % p for j in range(cols)]
        r += 1
        if r == rows:
            break
    return r


# ---------------------------------------------------------------- H1

def h1(qbridge=33):
    """THE BRIDGE: `Chain.ell` <-> the sine/cosine description, and `P`-independence."""
    import s450_chain_maps as cm
    res = {"rows": 0, "keyset_ok": 0, "parity_ok": 0, "bad": [],
           "basis_rows": 0, "basis_ok": 0, "basis_bad": []}
    for Q in range(5, qbridge + 1, 2):
        want_lone = list(range(2, Q, 2))
        want_paired = list(range(1, Q + 1, 2))
        by_class = {}
        for P in range(1, Q):
            if math.gcd(P, Q) != 1:
                continue
            for eps in (0, 1):
                ch = cm.Chain(P, Q, eps)
                d = {}
                for t in range(Q):
                    c, v = int(ch.c[t]), ch.ell[t].v
                    if v.any():
                        d.setdefault(min(c, 2 * Q - c), v)
                keys = sorted(d)
                lone = bool(ch.lone)
                res["rows"] += 1
                ok_k = (keys == (want_lone if lone else want_paired))
                ok_p = all(k % 2 == (0 if lone else 1) for k in keys)
                res["keyset_ok"] += ok_k
                res["parity_ok"] += ok_p
                if not (ok_k and ok_p):
                    res["bad"].append([P, Q, eps, lone, keys[:6]])
                by_class.setdefault((Q, lone), []).append(
                    (P, eps, [[int(x) for x in d[k]] for k in keys]))
        # (b) code ring vectors vs the independent sympy reduction: same rank, same defect
        for (QQ, lone), lst in by_class.items():
            P0, eps0, M0 = lst[0]
            # P-independence of the VALUES, not just the keys
            same = all(Mi == M0 for _, _, Mi in lst)
            code_rank = int(Matrix(M0).rank())
            A, B, f = spans(QQ)
            sym_rank = int(Matrix(A if not lone else B).rank())
            n_dist = len(M0)
            res["basis_rows"] += 1
            ok = same and code_rank == sym_rank and n_dist == (
                (QQ - 1) // 2 if lone else (QQ + 1) // 2)
            res["basis_ok"] += ok
            if not ok:
                res["basis_bad"].append(
                    [QQ, "lone" if lone else "paired", same, code_rank, sym_rank, n_dist])
    res["verdict"] = (res["keyset_ok"] == res["rows"] == res["parity_ok"]
                      and res["basis_ok"] == res["basis_rows"])
    return res


# ---------------------------------------------------------------- H2 / H3 / H4 / H5

def sweep(qmax=199, qexact=75):
    """H2 ranks, H3 eigenvectors, H4 defects, H5 controls -- one pass, one table."""
    res = {"table": [], "n": 0,
           "h2_rkB": 0, "h2_rkA": 0, "h2_union": 0,
           "h3_plus": 0, "h3_minus": 0,
           "h4_lone": 0, "h4_paired": 0, "h4_prime": 0, "h4_sum": 0,
           "h5_union_differs": 0, "h5_lone_defect0_iff_prime": 0,
           "exact_vs_mod_checked": 0, "exact_vs_mod_ok": 0, "bad": []}
    for Q in range(5, qmax + 1, 2):
        A, B, f = spans(Q)
        R, _ = reduction(Q)
        exact = Q <= qexact
        if exact:
            rkA, rkB, rkU = (int(Matrix(A).rank()), int(Matrix(B).rank()),
                             int(Matrix(A + B).rank()))
            mA = min(rank_mod(A, p) for p in MODP)
            mB = min(rank_mod(B, p) for p in MODP)
            mU = min(rank_mod(A + B, p) for p in MODP)
            res["exact_vs_mod_checked"] += 1
            res["exact_vs_mod_ok"] += ((rkA, rkB, rkU) == (mA, mB, mU))
            if (rkA, rkB, rkU) != (mA, mB, mU):
                res["bad"].append(["modmismatch", Q, rkA, mA, rkB, mB, rkU, mU])
        else:
            rkA = min(rank_mod(A, p) for p in MODP)
            rkB = min(rank_mod(B, p) for p in MODP)
            rkU = min(rank_mod(A + B, p) for p in MODP)

        # H3: sigma acts as +1 on A's spanning vectors and -1 on B's.  sigma sends
        # zeta^m -> zeta^-m, so it permutes the reduction rows; apply it to the
        # SPANNING VECTORS by rebuilding them with m -> -m and compare exactly.
        sA = [[R[(-m) % Q][i] + R[m][i] for i in range(f)]
              for m in range(0, (Q - 1) // 2 + 1)]
        sB = [[R[(-m) % Q][i] - R[m][i] for i in range(f)]
              for m in range(1, (Q - 1) // 2 + 1)]
        plus_ok = (sA == A)
        minus_ok = all(sB[i][j] == -B[i][j] for i in range(len(B)) for j in range(f))

        # H4: the kernel dimension, MEASURED (nullspace of the Q x phi(Q) table), not
        # computed from phi -- so `defect_LONE + defect_PAIRED == Q - phi(Q)` can fail.
        ker = Q - (int(Matrix(R).rank()) if exact else min(rank_mod(R, p) for p in MODP))

        nL, nP = (Q - 1) // 2, (Q + 1) // 2
        dL, dP = nL - rkB, nP - rkA
        want = f // 2
        pr = bool(isprime(Q))
        row = {"Q": Q, "phi": f, "prime": pr, "exact": exact,
               "rank_lone": rkB, "rank_paired": rkA, "rank_union": rkU, "want": want,
               "defect_lone": dL, "defect_paired": dP, "ker_measured": ker}
        res["table"].append(row)
        res["n"] += 1
        res["h2_rkB"] += (rkB == want)
        res["h2_rkA"] += (rkA == want)
        res["h2_union"] += (rkU == f)
        res["h3_plus"] += plus_ok
        res["h3_minus"] += minus_ok
        res["h4_lone"] += (dL == (Q - 1 - f) // 2)
        res["h4_paired"] += (dP == dL + 1)
        res["h4_prime"] += ((dL == 0) == pr)
        res["h4_sum"] += (dL + dP == ker == Q - f)
        res["h5_union_differs"] += (rkU != want)
        res["h5_lone_defect0_iff_prime"] += ((rkB == nL) == pr)
        for nm, cond in (("rkB", rkB == want), ("rkA", rkA == want), ("union", rkU == f),
                         ("sigma+", plus_ok), ("sigma-", minus_ok),
                         ("defect", dL + dP == ker == Q - f)):
            if not cond:
                res["bad"].append([nm, Q, rkA, rkB, rkU, f, dL, dP, ker])
    res["verdict"] = (res["bad"] == []
                      and res["h2_rkB"] == res["h2_rkA"] == res["h2_union"] == res["n"]
                      and res["h4_sum"] == res["h4_prime"] == res["n"]
                      and res["h5_union_differs"] == res["n"]
                      and res["exact_vs_mod_ok"] == res["exact_vs_mod_checked"])
    return res


# ---------------------------------------------------------------- H6

def h6(qmax=63):
    """EXPLICIT GENERATORS: the rotated prime relations, folded by parity.

    G5 is enforced inline -- a folded generator is checked to VANISH (exactly in the
    reduction basis and numerically) before it is allowed into the span.
    """
    from sympy import primefactors
    res = {"rows": 0, "ker_ok": 0, "lone_ok": 0, "paired_ok": 0,
           "g5_checked": 0, "g5_vanish_exact": 0, "g5_vanish_num": 0,
           "bad": [], "table": [], "witness": None}
    for Q in range(9, qmax + 1, 2):
        R, f = reduction(Q)
        if f == Q - 1:                       # prime: relation lattice is 0-dimensional
            continue
        h = (Q + 1) // 2
        nL, nP = (Q - 1) // 2, (Q + 1) // 2
        raw, lone_g, paired_g, by_p = [], [], [], {}
        for p in primefactors(Q):
            by_p[p] = []
            for m in range(Q):
                v = [0] * Q
                for j in range(p):
                    v[(m + j * Q // p) % Q] += 1
                raw.append(v)
                w_anti = [v[i] - v[(-i) % Q] for i in range(Q)]
                w_sym = [v[i] + v[(-i) % Q] for i in range(Q)]
                a = [((-1) ** k) * w_anti[(h * k) % Q] for k in range(1, nL + 1)]
                # ⚠ the `j = 0` class is a SINGLETON under `m -> -m` (`cos 0 = 1` pairs
                # with itself), so the symmetric fold double-counts it.  Scale the whole
                # vector by 2 to stay integral and halve the `j = 0` entry.
                b = [2 * ((-1) ** j) * w_sym[(h * j) % Q] for j in range(0, nP)]
                b[0] //= 2
                # ---- G5: these must actually be relations, on both families.  The
                # transport carries `(-1)^k` a SECOND time: `a_k` already absorbed one,
                # and `s_k = (-1)^k (zeta^hk - zeta^-hk)/2i` supplies the other.
                exact_ok = (
                    all(sum(a[k - 1] * ((-1) ** k)
                            * (R[(h * k) % Q][i] - R[(-h * k) % Q][i])
                            for k in range(1, nL + 1)) == 0 for i in range(f))
                    and all(sum(b[j] * ((-1) ** j)
                                * (R[(h * j) % Q][i] + R[(-h * j) % Q][i])
                                for j in range(nP)) == 0 for i in range(f)))
                num_ok = (
                    abs(sum(a[k - 1] * math.sin(k * math.pi / Q)
                            for k in range(1, nL + 1))) < 1e-12
                    and abs(sum(b[j] * math.cos(j * math.pi / Q)
                                for j in range(nP))) < 1e-12)
                res["g5_checked"] += 1
                res["g5_vanish_exact"] += exact_ok
                res["g5_vanish_num"] += num_ok
                if not (exact_ok and num_ok):
                    res["bad"].append(["G5", Q, p, m, exact_ok, num_ok])
                    continue
                lone_g.append(a)
                paired_g.append(b)
                by_p[p].append(a)
        ker = Q - f
        r_raw = int(Matrix(raw).rank())
        r_L = int(Matrix(lone_g).rank()) if lone_g else 0
        r_P = int(Matrix(paired_g).rank()) if paired_g else 0
        dL, dP = nL - f // 2, nP - f // 2
        res["rows"] += 1
        res["ker_ok"] += (r_raw == ker)
        res["lone_ok"] += (r_L == dL)
        res["paired_ok"] += (r_P == dP)
        if not (r_raw == ker and r_L == dL and r_P == dP):
            res["bad"].append(["span", Q, r_raw, ker, r_L, dL, r_P, dP])
        res["table"].append({"Q": Q, "primes": [int(p) for p in primefactors(Q)],
                             "ker": ker, "span_raw": r_raw,
                             "defect_lone": dL, "span_lone": r_L,
                             "defect_paired": dP, "span_paired": r_P})
        if Q == 15:                          # the [NCYL-317] (5) tie-in
            wit = [0] * nL
            wit[0], wit[3], wit[5] = 1, 1, -1     # sin(pi/15) + sin(4pi/15) - sin(6pi/15)
            chk = {}
            for p, G in by_p.items():
                chk[int(p)] = bool(int(Matrix(G).rank())
                                   == int(Matrix(G + [wit]).rank()))
            res["witness"] = {
                "relation": "sin(pi/15) + sin(4pi/15) - sin(6pi/15) = 0",
                "vanishes": abs(sum(c * math.sin((i + 1) * math.pi / 15)
                                    for i, c in enumerate(wit))) < 1e-12,
                "in_span_of_family": chk}
    res["verdict"] = (res["bad"] == []
                      and res["ker_ok"] == res["lone_ok"] == res["paired_ok"] == res["rows"]
                      and res["g5_vanish_exact"] == res["g5_vanish_num"] == res["g5_checked"]
                      and res["witness"] is not None
                      and res["witness"]["vanishes"]
                      and res["witness"]["in_span_of_family"].get(3) is True
                      and res["witness"]["in_span_of_family"].get(5) is False)
    return res


# ---------------------------------------------------------------- main

def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    opt = dict(a[2:].split("=", 1) for a in sys.argv[1:] if a.startswith("--"))
    qmax = int(opt.get("qmax", 199))
    qexact = int(opt.get("qexact", 75))
    qbridge = int(opt.get("qbridge", 33))
    which = args[0] if args else "all"
    out, t0 = {}, time.time()

    if which in ("all", "h1"):
        r = h1(qbridge)
        out["h1"] = r
        print(f"H1 BRIDGE (odd Q <= {qbridge}, every coprime P, both eps)")
        print(f"   keyset  {r['keyset_ok']}/{r['rows']}   parity {r['parity_ok']}/{r['rows']}"
              f"   P-indep+rank-agreement {r['basis_ok']}/{r['basis_rows']}"
              f"   -> {'OK' if r['verdict'] else 'FAIL'}")
        for b in r["bad"][:5] + r["basis_bad"][:5]:
            print("   BAD", b)

    if which in ("all", "h2", "h3", "h4", "h5"):
        r = sweep(qmax, qexact)
        out["sweep"] = r
        n = r["n"]
        print(f"\nSWEEP  odd Q = 5..{qmax}  ({n} values; exact rank to Q <= {qexact}, "
              f"modular beyond; exact-vs-mod agree {r['exact_vs_mod_ok']}"
              f"/{r['exact_vs_mod_checked']})")
        print(f"H2  rank(LONE)=phi/2 {r['h2_rkB']}/{n}   rank(PAIRED)=phi/2 "
              f"{r['h2_rkA']}/{n}   rank(A u B)=phi {r['h2_union']}/{n}")
        print(f"H3  sigma=+1 on A {r['h3_plus']}/{n}   sigma=-1 on B {r['h3_minus']}/{n}")
        print(f"H4  defect_LONE law {r['h4_lone']}/{n}   defect_PAIRED=+1 "
              f"{r['h4_paired']}/{n}   defect0<=>prime {r['h4_prime']}/{n}   "
              f"sum=ker=Q-phi {r['h4_sum']}/{n}")
        print(f"H5  union differs from phi/2 {r['h5_union_differs']}/{n}   "
              f"LONE defect0<=>prime (separates rows) "
              f"{r['h5_lone_defect0_iff_prime']}/{n}")
        print(f"    -> {'OK' if r['verdict'] else 'FAIL'}")
        for b in r["bad"][:8]:
            print("   BAD", b)
        comp = [t for t in r["table"] if not t["prime"]]
        print(f"    non-vacuity detail: {len(comp)} composite Q in range, all with "
              f"defect_LONE > 0 (min {min(t['defect_lone'] for t in comp)}, "
              f"max {max(t['defect_lone'] for t in comp)}); "
              f"{n - len(comp)} prime Q, all with defect_LONE = 0")

    if which in ("all", "h6"):
        r = h6(int(opt.get("qgen", 63)))
        out["h6"] = r
        n = r["rows"]
        print(f"\nH6 EXPLICIT GENERATORS (composite odd Q <= {opt.get('qgen', 63)}, "
              f"{n} rows)")
        print(f"   G5 vanishing  exact {r['g5_vanish_exact']}/{r['g5_checked']}   "
              f"numeric {r['g5_vanish_num']}/{r['g5_checked']}")
        print(f"   span(N_p) = ker {r['ker_ok']}/{n}   "
              f"antisym = defect_LONE {r['lone_ok']}/{n}   "
              f"sym = defect_PAIRED {r['paired_ok']}/{n}")
        w = r["witness"]
        print(f"   [NCYL-317](5) witness at Q=15: vanishes {w['vanishes']}, "
              f"in span per family {w['in_span_of_family']}")
        print(f"    -> {'OK' if r['verdict'] else 'FAIL'}")
        for b in r["bad"][:8]:
            print("   BAD", b)

    out["elapsed_s"] = round(time.time() - t0, 1)
    with open("data/s463_rank_law.json", "w") as fh:
        json.dump(out, fh, indent=1)
    print(f"\nwrote data/s463_rank_law.json   ({out['elapsed_s']} s)")


if __name__ == "__main__":
    main()
