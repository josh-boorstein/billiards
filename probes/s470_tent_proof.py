#!/usr/bin/env python3
"""s470_tent_proof.py -- PRE-REGISTERED.  Queue item (s469-a): the TENT step-down-by-`2`.

⇒⇒ WHAT THIS FILE IS.  A PROOF, plus the guards on its joints.  [NCYL-326] (4) records the
tent half of `turns = 1` as the part s469's theorem does NOT reach: `98` rows, `2450`
prongs, `0 CAP`, measured to odd `Q <= 101`, with the named next lemma *at successive folds
the smaller interface's profile value steps down the value AP by exactly `2`* and the
explicit warning that it is NOT a one-liner about nesting (tents are not total chains,
`24/48`).  ⇒⇒ **THE WARNING IS RIGHT AND THE LEMMA IS STILL SHORT: THE MISSING INGREDIENT
IS NOT AN ORDER FACT ABOUT THE INTERVALS, IT IS THE STRICT CONCAVITY OF `sin` ON
`[0, π/2]`.**  With it the fold levels, the fold COUNT and the escape all follow by one
induction, at every odd `Q`, with no cap.

⇒⇒⇒ CONSEQUENCE.  s469 proved the `turns = 0` rows (`P ∈ {1, Q-1}`) and the VALLEY half of
`turns = 1`.  This file proves the TENT half, which is the whole remainder of `turns <= 1`.
⇒ **(G0b) IS NOW A THEOREM ON THE ENTIRE `turns <= 1` LOCUS = `P ∈ {1, 2, Q-2, Q-1}` = the
VEECH rows ([NCYL-012]'s criterion `min(P,Q-P) ∈ {1,2}`, which [NCYL-327] identifies with
`turns <= 1`), at every odd `Q`.**
⚠⚠ AND IT CLOSES NOTHING WHERE (G0b) IS OPEN: [NCYL-325]'s `114` capped prongs are all on
NON-Veech rows, `turns >= 2`, and this argument does not reach them -- with two interior
extrema the two-ramp decomposition below is simply false.  No `claims.md` tier moves and no
γ=1 DAG node moves.

--------------------------------------------------------------------------------------
PRIOR ART: run with `rulings.py` (s469 left this gate explicitly NOT RUN for (s469-a), and
DO-NOT (8) bars delegating it).  `--grep` on 'step-down', 'fold edge', 'value ladder' ->
[NCYL-326] alone,
i.e. the statement of the item and nothing else.  'tent' -> noise only ([L3a-009], [NB-*]),
which is s469 DO-NOT (7)'s documented trap (`content`/`extent`/`consistent`) firing a third
time.  'unimodal' -> [WFLOOR-013]/[WFLOOR-023]/[WFLOOR-055]: a DIFFERENT object (the F-leg
width ladder `w_m`, resp. the residue ladder `A(n)`) -- ⚠ but [WFLOOR-055] is worth naming
because it is the same SHAPE of argument in the other strand (a sine family unimodal along
a run, with WHICH endpoint wins settled by an integer comparison); no shared object, no
shared code, and nothing there is imported here.  'cascade' -> [SELECTOR-002]/[NB-040]/
[NB-043], the seam cascade, unrelated.  'arithmetic progression' -> [NCYL-279]/[WFLOOR-106]/
[OPS-116], none about interface lengths.  ⇒ ID-grep per [OPS-246] on [NCYL-292]/[NCYL-294]/
[NCYL-297]/[NCYL-324]/[NCYL-326]: [NCYL-292] has the nesting at the ONE extremal kite,
[NCYL-324] the whole profile, [NCYL-326] the definitional proof of both -- **no line
anywhere states T2's CONCENTRICITY, computes `osec` in closed form on any row family, or
uses concavity of `sin` on this object.**  ⇒ **THE ONE HIT THAT CHANGED THE DESIGN:
[NCYL-326] (4)'s total-chain control (`24/48`) reads as *the interval order is not usable*;
T2 shows the `24` failures are the DEGENERATE node `a_0 = 0` alone, so the usable order was
there the whole time.  Recorded as a forward marker on [NCYL-326], not as a refutation --
the control's VERDICT (not a total chain) is correct, only its reading was.**
⚠ [OPS-242]'s gate does not bind: named-region statement, and `7/15 eps1` has `turns = 6`.
⚠ [NCYL-315] (renormalization) and [NCYL-296] (CF depth) do not bite -- no CF digit is read
and no return map is induced.

--------------------------------------------------------------------------------------
THE PROOF.

Notation is s469's: `m = (Q-1)/2`; path node `i` is interface `m-i`; `a_i = w(x_i)` with
`w(x) = min(x, 2Q-x)` and `x_i = (c0 + 2P(m-i)) mod 2Q`; `h_a := sin(aπ/2Q)`, so the node
lengths are `|J_i| = h_{a_i}` and comparisons of interfaces are comparisons of `a`.
`ν(k) := min(a_k, a_{k+1})` is the SMALLER value at edge `k` -- the quantity [NCYL-326] (4)
names.  A TENT is a row whose profile has exactly one interior extremum and it is a MAXIMUM.

T1 (WHICH ROWS ARE TENTS -- exact, from s469's L3/L4).  `turns = min(P,Q-P) - 1`, so
   `turns = 1` iff `P ∈ {2, Q-2}`; and s469 L3(a) gives `x_0 ∈ {0, Q}`, with `x_0 = 0` iff
   the profile starts at its MINIMUM, i.e. iff the row is a tent.  For `c0 = P` (`eps = 0`)
   `x_0 = PQ mod 2Q` is `0` for even `P` and `Q` for odd `P`, and `eps = 1` swaps them.
   `P = 2` is even and `P = Q-2` is odd.  ⇒⇒ **the tent rows are EXACTLY `(P,eps) = (2,0)`
   and `(Q-2,1)`, two per odd `Q`** -- and the other two `turns = 1` rows are the valleys
   s469 already proved.

T2 (CONCENTRICITY -- exact arithmetic on `osec`, NOT a measurement, and it is what the
   `24/48` total-chain control was really seeing).  `Chain.__init__` sets
   `osec[k] = 0 < (2Q - c[k-1]) mod 2Q < 2P`, i.e. `osec[k]` iff `c[k-1] > 2Q - 2P`;
   `build` then gives `lo[near] = lo[far]` when `osec[k]` is false and `hi[near] = hi[far]`
   when it is true.  Evaluate it on the two tent families (`c[t] = (c0 + 2Pt) mod 2Q`, all
   `c[t]` even, path kites are `1..m`):
     (2,0):    `c0 = 2`, `c[t] = (2+4t) mod 2Q`.  `osec[k]` needs `c[k-1] ∈ (2Q-4, 2Q)`,
               whose only even member is `2Q-2`; `2+4t ≡ 2Q-2` gives `2t ≡ -2 (mod Q)`,
               i.e. `t = Q-1`, i.e. `k = 0`.  ⇒ `osec` is FALSE on every path kite.
     (Q-2,1):  `c0 = 2Q-2`, `c[t] = (2Q-2-4t) mod 2Q`.  `osec[k]` needs `c[k-1] > 4`, so it
               fails only at `c[k-1] ∈ {0,2,4}`, i.e. `t = m, Q-1, m-1`, i.e. `k ∈ {m+1, 0,
               m}`.  ⇒ on the path kites `1..m`, `osec` is TRUE except at kite `m`.
   Path edge `j` is kite `m-j`, so kite `m` is EDGE `0`.  ⇒⇒ **every path edge except edge
   `0` shares the SAME endpoint (LO in the first family, HI in the second), so the nodes are
   CONCENTRIC: there is one anchor `c` with `J_i = [c, c+h_{a_i}]` for all `i` (family 1) or
   `J_i = [c-h_{a_i}, c]` for `i >= 1` (family 2).**
   ⚠⚠ **AND THE EXCEPTION IS NOT A GAP.  On a tent `a_0 = 0` (T1), so node `0` is a single
   POINT.**  In family 2 it sits at `Lo[1] = c - h_{a_1}`, off the anchor -- which is the
   sole reason `_is_total_chain` reported `24/48` -- but it costs nothing: `Lo[0] = Hi[0]`
   makes the walk's strict-containment test `Lo[0] < y < Hi[0]` VACUOUS, so node `0` is
   never crossed, exactly as `r < h_0 = 0` predicts; and edge `0`'s fold obeys the same
   recursion as every other (write `S_0 = Hi[0] + Hi[1] = Lo[0] + c` and substitute:
   `r' = c - (S_0 - y) = y - Lo[0] = h_{a_1} - r = h_0 + h_{a_1} - r`).

T3 (THE LEVEL COORDINATE).  Put `r := |y - c|`.  By T2, for every node `i >= 1`
   `y ∈ int J_i ⟺ r < h_{a_i}`, and node `0` is never crossed.  A fold at edge `k` sends
   `y ↦ S_k - y`; substituting T2's two forms gives, in BOTH families and at edge `0` too,
        **`r ↦ h_near + h_far - r`**,
   which maps the annulus `(h_near, h_far)` to itself.  ⇒ the whole walk is a 1-D system in
   the single scalar `r`, compared against the `m+1` numbers `h_{a_i}`.

T4 (THE TWO RAMPS INTERLEAVE, mod 4).  A tent is strictly increasing on `0..q` and strictly
   decreasing on `q..m` (`q` = the peak index, `a_q = Q-1`).  Take family 1 (`P = 2`,
   `c0 = 2`): `x_i = (2 + 4(m-i)) mod 2Q`, so `x` steps by `4` and `a = w(x)` is affine on
   each of the arcs `(0,Q)`, `(Q,2Q)`.  On the DOWN ramp `x_i < Q` and `a_i = x_i ≡ 2
   (mod 4)`; on the UP ramp `x_i > Q` and `a_i = 2Q - x_i ≡ 0 (mod 4)`, since `2Q ≡ 2
   (mod 4)` for odd `Q`.  Together with s469 L3(b) (the value set is `{0,2,…,Q-1}`, each
   value once) this says:
     ⇒⇒ **the up-ramp values are exactly the residues `≡ 0 (mod 4)` and the down-ramp values
     exactly those `≡ 2 (mod 4)`, so SORTED BY VALUE THE TWO RAMPS ALTERNATE, step `2`.**
   Corollaries used below: for a non-peak value `ν`, `ν+2` lies on the OTHER ramp and `ν+4`
   on the SAME one; and `ν`'s same-ramp neighbour toward the peak has value `ν+4`, EXCEPT
   for `ν = Q-3` (the top of the ramp not carrying the peak), whose neighbour is the peak
   `Q-1 = ν+2`.  Family 2 is the same computation with `c0 = 2Q-2`, `P = Q-2`.

T5 (CONCAVITY -- THE ONE ANALYTIC INPUT, AND THE WHOLE CONTENT).  Every value satisfies
   `0 <= a <= Q`, so `aπ/2Q ∈ [0, π/2]`, where `sin` is STRICTLY concave.  Hence for every
   value `ν` with `ν-2, ν+2 ∈ [0,Q]`,
        **`h_{ν-2} + h_{ν+2} < 2 h_ν`.**

T6 (THE INVARIANT, AND THE INDUCTION).  Let the corner prong of edge `j` have `ν_j :=
   ν(j)`.  `corner_start` places it at the NEAR interface's non-shared endpoint, i.e. at
   `r_0 = h_{ν_j}` exactly, on the FAR node, heading away from edge `j` -- which by T4 is
   toward the peak.  Claim:
        **`I(ν)`:  after the fold at level `ν`, `h_ν < r < h_{ν+2}`,**
   with the base case `I(ν_j)` read in the closed form `h_{ν_j} <= r_0 < h_{ν_j+2}` (both
   sides are then trivial monotonicity of `h`).
   ⇒ *Given `I(ν)`, the next fold is at level `ν-2`.*  The prong runs toward the peak,
   crossing every node it meets on its own ramp (all have value `> ν+4 > ...`, hence
   `h > r`), crosses the peak, and descends the OTHER ramp, whose values in decreasing order
   are `…, ν+6, ν+2, ν-2, …` (T4).  It CROSSES the `ν+2` node because `r < h_{ν+2}`, and it
   FOLDS at the `ν-2` node because `r > h_ν > h_{ν-2}`.  That edge has near value `ν-2` and
   far value `ν+2` (T4), so by T3
        `r' = h_{ν-2} + h_{ν+2} - r`.
   ⇒ *`I(ν-2)` holds.*   `r' > h_{ν-2} ⟺ r < h_{ν+2}` ✓ (that is `I(ν)`);  and
        `r' < h_ν  ⟺  r > h_{ν-2} + h_{ν+2} - h_ν`,
   whose right side is `< h_ν` **by T5**, while `r > h_ν` by `I(ν)`. ✓
   ⚠ The single exception `ν = Q-3` (far `= ν+2`, T4) only makes the second inequality
   `r > h_{ν-2}`, which is weaker still.

⇒⇒⇒ THEOREM (the tent).  Let `Q >= 5` be odd and let the row be a tent -- by T1, `(P,eps) ∈
{(2,0), (Q-2,1)}`.  Then for every edge `j` carrying a corner prong (`j = 1 .. m-1`), that
prong's fold levels are EXACTLY
        `ν_j - 2, ν_j - 4, …, 4, 2, 0`   (step exactly `-2`),
it folds exactly `ν_j / 2` times, and it then ESCAPES.
*Proof.*  T6 by induction from `ν_j` down.  The induction cannot stall above `0`, since each
step lowers `ν` by `2` and `I(ν)` supplies the next fold whenever `ν > 0`.  At `ν = 0` the
near node is node `0`, the far value is `4`, and `I(0)` reads `0 < r < h_2`; `h_2` is the
SMALLEST positive node length, so from there the prong crosses every node and runs off the
end of the path.  ∎
⇒ COROLLARIES, both of which are figures [NCYL-326] (4) measured: `ν_j` runs over
`{2,4,…,Q-3}` as `j` runs over `1..m-1` (each edge carries a distinct `ν`, T4), so the
maximum fold count on a row is `(Q-3)/2 = m-1`, and the number of tent prongs over odd
`Q <= 101` is `Σ_Q 2(m-1) = Σ_Q (Q-3) = 2+4+…+98 = 2450`.
⇒ AND [NCYL-326] (4)'s own named consequences follow: the fold EDGES are distinct (their
`ν` are), the cascade is bounded by `m`, and the walk terminates.

--------------------------------------------------------------------------------------
HYPOTHESES (pre-registered).  G0-G5 are GUARDS: given the proof they cannot fail, and a
failure means I have misread the source ([OPS-222]).  G6 is the one arm whose OUTCOME I do
not know in advance.

  G0  T1: over odd `Q <= QBIG`, the tents are exactly `(2,0)` and `(Q-2,1)`.
      ⚠ NEGATIVE CONTROL: `(2,1)` and `(Q-2,0)` must ALL be valleys -- if they were tents
      too, T1 would be selecting on `P` alone and `eps` would be decoration.

  G1  T2: on every tent row, `osec` on path kites matches the closed form, and the BUILT
      endpoints are concentric off node `0`.  ⚠ NEGATIVE CONTROL: the `osec`-INVERTED
      prediction must FAIL on every row -- otherwise the test is not reading `osec`.
      ⚠ Also report how many rows `_is_total_chain` rejects and how many of those become
      chains once node `0` is dropped: the claim is `24/48` -> `0` rejections.

  G2  T4: up-ramp values `≡ 0 (mod 4)`, down-ramp `≡ 2 (mod 4)` (family 1; the complementary
      pair in family 2), sorted values alternate ramps with step `2`, and the same-ramp
      neighbour of `ν` toward the peak has value `ν+4` for every `ν` except `Q-3`.
      ⚠ NEGATIVE CONTROL: the flat guess *far `= ν+2` always* must fail on every row but at
      the single value `Q-3`.

  G3  T5: the concavity margin `min_ν (2h_ν - h_{ν-2} - h_{ν+2})` is `> 0`, with its
      closed-form value.  ⚠⚠ NEGATIVE CONTROL, AND IT IS THE ONE THAT LOCATES THE CONTENT:
      re-run the margin with the LINEAR surrogate `h̃_a = a`.  It must come out EXACTLY `0`
      -- i.e. the induction has no slack without strict concavity, so T5 is load-bearing
      and not decoration.  Report the surrogate's failure count for `I(ν)` as well.

  G4  T6: on the REAL exact walk, after every fold at level `ν`, `h_ν < r < h_{ν+2}` -- and
      `r` is recomputed from the anchor with the audited exact `cmp`, not from the model.
      ⚠ NEGATIVE CONTROL: the parity-shifted claim `h_{ν-2} < r < h_ν` must fail on
      essentially every fold.

  G5  ⇒⇒ THE THEOREM.  Every corner prong of every tent row, odd `Q <= QMAX`: end `ESC`,
      fold-level sequence exactly `[ν_j-2, …, 0]`, fold count exactly `ν_j/2`, per-row max
      exactly `m-1`, and the prong total reproducing [NCYL-326] (4)'s `2450` at `Q <= 101`.
      ⚠ The predicted sequence is computed from the PROFILE alone -- pure arithmetic, no
      chain, no ring -- and scored against the walk, which shares no code with it.
      ⚠ NEGATIVE CONTROL: the same prediction applied to the VALLEY rows must fail
      everywhere (valleys fold `0` times, [NCYL-326] (2)).

  G6  ⇒⇒ THE ARM.  Does the argument reach ANY `turns >= 2` row?  T4's two-ramp
      decomposition is what fails there, so the pre-registered expectation is NO -- but the
      informative question is WHICH clause dies first, and I do not predict that.  Score, on
      `turns = 2` rows: (i) is the value set still an AP of step `2`?  (ii) do the ramps
      still partition the values into residue classes?  (iii) does `I(ν)` survive on the
      real walk?  ⚠⚠ A `YES` on (iii) would be the interesting outcome and would say the
      invariant is not tied to the two-ramp picture; I expect NO and I expect (ii) to be the
      first casualty.

--------------------------------------------------------------------------------------
RESULTS (`data/s470_tent_proof.json`, `logs/s470_tent_proof.log`; `--qmax=101 --qbig=401`).

  G0  `198/198` tents predicted, `0` unpredicted, `0` missed, odd `Q <= 401`.
      ⇒ CONTROL FIRES: `0/198` of the other-`eps` rows `(2,1)`/`(Q-2,0)` are tents.
  G1  `osec` matches T2's closed form on `2548/2548` path kites; the BUILT endpoints are
      concentric off node `0` on `98/98` rows; node `0` degenerate on `98/98`.
      ⇒ CONTROL FIRES: the `osec`-inverted prediction holds on `0/98`.
      ⇒⇒ **AND THE READING OF [NCYL-326] (4)'s CONTROL IS CORRECTED: `_is_total_chain`
      accepts `49/98` rows, but `98/98` once node `0` is dropped.** So the `24/48` that made
      the invariant look like it needed a subtler cause is the DEGENERATE node alone -- the
      verdict was right, the reading was not.  The usable order was there the whole time.
  G2  residue split `398/398`, alternation `398/398`, AP-step-`2` `398/398`;
      `far = ν+4` on `39800` values and `far = ν+2` on exactly `398` -- one per row, every
      one at `ν = Q-3`, as T4 says.  ⇒ CONTROL FIRES: *far `= ν+2` always* holds on
      `398/40198`.
  G3  concavity on `39800/39800` triples; worst margin `4.808502e-07` at `Q = 401, ν = 2`,
      matching the closed form `2 h_ν (1 - cos(π/Q))` to 12 digits.
      ⇒⇒ CONTROL FIRES AND IT LOCATES THE CONTENT: the LINEAR surrogate `h̃_a = a` gives
      margin EXACTLY `0` on `39800/39800`.  The induction has zero slack without strict
      concavity, so T5 is load-bearing and not decoration.
  G4  ⇒ `41650/41650` folds satisfy `h_ν < r` AND `41650/41650` satisfy `r < h_{ν+2}`, with
      `r` recomputed from the anchor by the audited exact `cmp`.
      ⇒ CONTROL FIRES: the parity-shifted claim holds on `0/41650`.
  G5  ⇒⇒⇒ THE THEOREM.  `98` rows, **`2450/2450` prongs ESC**, `0` R / `0` Z / `0` CAP;
      fold-level sequence equals the profile-only prediction on `2450/2450`; fold count
      equals `ν_j/2` on `2450/2450`; per-row max `= m-1` on `98/98`.
      ⚠ The prong total `2450` and the fold total `41650` are DERIVED, not counted:
      `Σ_Q 2(m-1) = Σ_Q (Q-3) = 2+4+…+98 = 2450` and `Σ_Q 2·Σ_{j}ν_j/2 = Σ_{m=2}^{50}
      m(m-1) = 2·C(51,3) = 41650`.  `2450` is [NCYL-326] (4)'s own measured figure.
      ⇒ CONTROL FIRES: on the VALLEY rows the tent prediction holds on `196/2548`, and
      `2548/2548` fold ZERO times ([NCYL-326] (2)).  ⚠ THE `196` ARE VACUOUS AGREEMENTS,
      NOT PARTIAL SUPPORT -- they are the two edges per row with `ν = 1`, where the
      prediction is the EMPTY sequence and zero folds meets it trivially.
  G6  ⇒⇒ THE ARM, AND MY NAMED PREDICTION WAS ONLY HALF RIGHT.  On `turns = 2` (`52` rows,
      `610` prongs, `Q <= 45`):
        (i)   the AP-of-step-`2` value set SURVIVES -- `52/52`.
        (ii)  the residue split DIES -- `0/52`.  ⇐ this is the casualty I named.
        (T2)  ⚠⚠ **CONCENTRICITY ALSO DIES -- `0/52` -- AND I DID NOT PREDICT THAT.**
              `osec` is genuinely mixed on the path kites there (e.g. `[0,0,0,1,0]` at
              `3/11 eps0`), so there is no anchor and the scalar `r` of T3 DOES NOT EXIST.
              That is UPSTREAM of (ii): I predicted the failure in the combinatorial layer
              and the first failure is in the geometric one.
        (iv)  the coordinate-free conclusion dies OUTRIGHT: `0/362` prong sequences are
              strictly decreasing and `0/362` step by `-2`.  Not a graceful degradation.
      ⊕ the size of the gap, on rows ONE step outside the proved region: `1709600` folds on
      `610` prongs (`2803` each) against `41650` on `2450` (`17` each), and `510/610` ESC
      within a `60000`-step cap.
      ⇒⇒ **WHAT THIS BUYS (s468-b): the open locus does not merely lack the interleaving --
      it lacks the single scalar coordinate.  An instrument there cannot be a 1-D
      comparison, which is what every clause of this proof is.**

⚠⚠ SCOPE.  This closes the tent half of `turns = 1` and nothing else.  With s469 it makes
(G0b) a theorem on `turns <= 1` = `P ∈ {1,2,Q-2,Q-1}` = the Veech rows, at every odd `Q`.
[NCYL-325]'s `114` capped prongs are all at `turns >= 2` and G6 says why none of this
reaches them.  No `claims.md` tier moves; no γ=1 DAG node moves; `(G0b-form)` untouched.

RUN: PYTHONPATH=.:engine:archive/scripts_2026-07:archive/scripts_2026-08:probes \
     .venv/bin/python3.13 probes/s470_tent_proof.py [all|g0..g6] [--qmax=101] [--qbig=401]
"""
from __future__ import annotations

import json
import math
import sys

import s450_chain_maps as cm          # noqa: F401  (Chain: the definitions T2 quotes)
import s451_quotient_path as qp
import s453_pole_module as pm
import s456_neckgb as nb

QMAX = 101           # the exact-walk box (2 tent rows per odd Q)
QBIG = 401           # the pure-arithmetic box for T1/T4
OUT = "data/s470_tent_proof.json"


# ------------------------------------------------------------------ shared

def profile(P, Q, eps):
    """s469 L3's folded rotation `(x_i, a_i)`, verbatim -- pure arithmetic."""
    twoQ, m = 2 * Q, (Q - 1) // 2
    c0 = P if eps == 0 else (P + Q) % twoQ
    x = [(c0 + 2 * P * (m - i)) % twoQ for i in range(m + 1)]
    return x, [min(v, twoQ - v) for v in x], m


def turns_of(a):
    s = [1 if a[i + 1] > a[i] else -1 for i in range(len(a) - 1)]
    return sum(1 for i in range(len(s) - 1) if s[i] != s[i + 1])


def shape_of(a):
    t = turns_of(a)
    if t == 0:
        return "mono"
    if t > 1:
        return "multi"
    return "tent" if a[1] > a[0] else "valley"


def tent_rows(qmax, qmin=5):
    for Q in range(qmin, qmax + 1, 2):
        for P, eps in ((2, 0), (Q - 2, 1)):
            if Q >= 5:
                yield P, Q, eps


def all_rows(qmax, qmin=5):
    for Q in range(qmin, qmax + 1, 2):
        for P in range(1, Q):
            if math.gcd(P, Q) == 1:
                for eps in (0, 1):
                    yield P, Q, eps


def chain_of(P, Q, eps):
    ch, lo, hi, closes = qp.build(P, Q, eps)
    if not closes:
        return None
    Lo, Hi, S = pm.path_view(ch, lo, hi, pm.fold_sums(ch, lo, hi))
    return ch, Lo, Hi, S, pm.CoordOps(ch)


def h(a, Q):
    return math.sin(a * math.pi / (2 * Q))


def nu_of(a, k):
    return min(a[k], a[k + 1])


def predicted(a, j):
    """T6's conclusion, from the PROFILE alone: the fold-level sequence of prong `j`."""
    v = nu_of(a, j)
    return list(range(v - 2, -1, -2))


# ------------------------------------------------------------------ G0

def g0(qbig=QBIG):
    """T1: the tents are exactly (2,0) and (Q-2,1)."""
    r = {"Q_range": qbig, "tent_rows": 0, "predicted_hit": 0, "predicted_miss": 0,
         "unpredicted_tent": 0, "ctrl_other_eps_tents": 0, "ctrl_other_eps_rows": 0,
         "bad": []}
    for P, Q, eps in all_rows(qbig):
        _x, a, _m = profile(P, Q, eps)
        sh = shape_of(a)
        want = (P, eps) in ((2, 0), (Q - 2, 1)) and Q >= 5
        if sh == "tent":
            r["tent_rows"] += 1
            if want:
                r["predicted_hit"] += 1
            else:
                r["unpredicted_tent"] += 1
                if len(r["bad"]) < 10:
                    r["bad"].append({"P": P, "Q": Q, "eps": eps, "shape": sh})
        elif want:
            r["predicted_miss"] += 1
            if len(r["bad"]) < 10:
                r["bad"].append({"P": P, "Q": Q, "eps": eps, "shape": sh})
        if (P, eps) in ((2, 1), (Q - 2, 0)):          # the negative control
            r["ctrl_other_eps_rows"] += 1
            r["ctrl_other_eps_tents"] += int(sh == "tent")
    return r


# ------------------------------------------------------------------ G1

def g1(qmax=QMAX):
    """T2: osec closed form, concentricity, and what the total-chain control was seeing."""
    r = {"rows": 0, "osec_kites": 0, "osec_matches_closed_form": 0,
         "concentric_off_node0": 0, "ctrl_inverted_ok": 0,
         "node0_degenerate": 0, "total_chain": 0, "chain_after_dropping_node0": 0,
         "bad": []}
    for P, Q, eps in tent_rows(qmax):
        c = chain_of(P, Q, eps)
        if c is None:
            continue
        ch, Lo, Hi, _S, ops = c
        m = ch.m_idx
        r["rows"] += 1
        # -- osec against the closed form of T2, on the path kites 1..m
        ok = True
        for k in range(1, m + 1):
            r["osec_kites"] += 1
            if P == 2:
                want = False
            else:                                     # P = Q-2
                want = (k != m)
            ok_k = (bool(ch.osec[k]) == want)
            r["osec_matches_closed_form"] += int(ok_k)
            ok &= ok_k
        # -- the BUILT endpoints: concentric off node 0
        if P == 2:
            conc = all(ops.eq(Lo[i], Lo[0]) for i in range(m + 1))
        else:
            conc = all(ops.eq(Hi[i], Hi[m]) for i in range(1, m + 1))
        r["concentric_off_node0"] += int(conc)
        # -- negative control: the INVERTED osec prediction must fail
        inv_ok = True
        for k in range(1, m + 1):
            want_inv = (P == 2) if P == 2 else (k == m)
            inv_ok &= (bool(ch.osec[k]) == want_inv)
        r["ctrl_inverted_ok"] += int(inv_ok)
        # -- node 0 degenerate, and the total-chain reading
        r["node0_degenerate"] += int(ops.eq(Lo[0], Hi[0]))
        r["total_chain"] += int(_is_chain(Lo, Hi, ops, range(m + 1)))
        r["chain_after_dropping_node0"] += int(_is_chain(Lo, Hi, ops, range(1, m + 1)))
        if not (ok and conc) and len(r["bad"]) < 10:
            r["bad"].append({"P": P, "Q": Q, "eps": eps, "osec_ok": ok, "conc": conc})
    return r


def _is_chain(Lo, Hi, ops, idx):
    idx = list(idx)
    for x in idx:
        for y in idx:
            if y <= x:
                continue
            if not ((ops.cmp(Lo[y], Lo[x]) <= 0 and ops.cmp(Hi[x], Hi[y]) <= 0)
                    or (ops.cmp(Lo[x], Lo[y]) <= 0 and ops.cmp(Hi[y], Hi[x]) <= 0)):
                return False
    return True


# ------------------------------------------------------------------ G2

def g2(qbig=QBIG):
    """T4: the ramps are residue classes mod 4, they alternate, and far = nu+4."""
    r = {"rows": 0, "residue_split_ok": 0, "alternates_ok": 0, "ap_step2_ok": 0,
         "far_is_nu_plus4": 0, "far_is_nu_plus2": 0, "far_other": 0,
         "ctrl_far_always_plus2": 0, "values_checked": 0, "bad": []}
    for P, Q, eps in tent_rows(qbig):
        _x, a, m = profile(P, Q, eps)
        if shape_of(a) != "tent":
            continue
        r["rows"] += 1
        q = a.index(max(a))
        up, dn = a[:q], a[q + 1:]                     # both EXCLUDING the peak
        ru = {v % 4 for v in up}
        rd = {v % 4 for v in dn}
        split = (len(ru) <= 1 and len(rd) <= 1 and ru != rd)
        r["residue_split_ok"] += int(split)
        sv = sorted(a)
        r["ap_step2_ok"] += int(all(sv[i + 1] - sv[i] == 2 for i in range(len(sv) - 1)))
        # alternation: consecutive sorted values lie on different ramps (peak excluded)
        side = {}
        for v in up:
            side[v] = "U"
        for v in dn:
            side[v] = "D"
        inner = [v for v in sv if v in side]
        r["alternates_ok"] += int(all(side[inner[i]] != side[inner[i + 1]]
                                      for i in range(len(inner) - 1)))
        # far value of each edge, against nu+4 / nu+2
        for k in range(m):
            v, f = nu_of(a, k), max(a[k], a[k + 1])
            r["values_checked"] += 1
            if f == v + 4:
                r["far_is_nu_plus4"] += 1
            elif f == v + 2:
                r["far_is_nu_plus2"] += 1
                if v != Q - 3 and len(r["bad"]) < 10:
                    r["bad"].append({"P": P, "Q": Q, "nu": v, "far": f, "why": "plus2 off Q-3"})
            else:
                r["far_other"] += 1
                if len(r["bad"]) < 10:
                    r["bad"].append({"P": P, "Q": Q, "nu": v, "far": f})
            r["ctrl_far_always_plus2"] += int(f == v + 2)
        if not split and len(r["bad"]) < 10:
            r["bad"].append({"P": P, "Q": Q, "up_res": sorted(ru), "dn_res": sorted(rd)})
    return r


# ------------------------------------------------------------------ G3

def g3(qbig=QBIG):
    """T5: the concavity margin, and the LINEAR surrogate control (must be exactly 0)."""
    r = {"rows": 0, "triples": 0, "concave": 0, "worst_margin": None, "worst_at": None,
         "closed_form_at_worst": None, "ctrl_linear_zero": 0, "ctrl_linear_nonzero": 0}
    worst = float("inf")
    for P, Q, eps in tent_rows(qbig):
        _x, a, _m = profile(P, Q, eps)
        if shape_of(a) != "tent":
            continue
        r["rows"] += 1
        for v in range(2, Q - 1, 2) if a[0] == 0 else range(3, Q, 2):
            if v - 2 < 0 or v + 2 > Q:
                continue
            r["triples"] += 1
            mar = 2 * h(v, Q) - h(v - 2, Q) - h(v + 2, Q)
            r["concave"] += int(mar > 0)
            if mar < worst:
                worst, r["worst_at"] = mar, {"P": P, "Q": Q, "eps": eps, "nu": v}
                # 2 sin t - sin(t-d) - sin(t+d) = 2 sin t (1 - cos d), d = pi/Q
                r["closed_form_at_worst"] = 2 * h(v, Q) * (1 - math.cos(math.pi / Q))
            lin = 2 * v - (v - 2) - (v + 2)           # the LINEAR surrogate
            r["ctrl_linear_zero"] += int(lin == 0)
            r["ctrl_linear_nonzero"] += int(lin != 0)
    r["worst_margin"] = worst
    return r


# ------------------------------------------------------------------ G4 + G5

def _anchor(P, Lo, Hi, m):
    """T3: (anchor, sign) so that r = sign*(y - anchor) >= 0."""
    return (Lo[0], +1) if P == 2 else (Hi[m], -1)


def walk_levels(a, Lo, Hi, S, ops, j, anchor, sign, m, hcoord, cap=400_000):
    """The exact walk of prong `j`; returns (end, [(nu, r_ok_lo, r_ok_hi, r_ok_shift)])."""
    cs = nb.corner_start(Lo, Hi, S, j, ops)
    if cs is None:
        return None, []
    (i, d, y), _near, _far = cs
    n = m + 1
    out = []
    for _ in range(cap):
        i2 = i + d
        if i2 < 0 or i2 >= n:
            return "ESC", out
        c0, c1 = ops.cmp(y, Lo[i2]), ops.cmp(y, Hi[i2])
        if c0 == 0 or c1 == 0:
            return "Z", out
        if c0 > 0 and c1 < 0:
            i = i2
            continue
        k = min(i, i2)
        if ops.eq(y, ops.half(S[k])):
            return "R", out
        v = nu_of(a, k)
        y = ops.sub(S[k], y)
        d = -d
        # G4: recompute r from the ANCHOR with the audited exact cmp
        rr = ops.sub(y, anchor) if sign > 0 else ops.sub(anchor, y)
        lo_ok = ops.cmp(rr, hcoord[v]) > 0
        hi_ok = (v + 2 not in hcoord) or ops.cmp(rr, hcoord[v + 2]) < 0
        sh_ok = (v - 2 in hcoord and ops.cmp(rr, hcoord[v - 2]) > 0
                 and ops.cmp(rr, hcoord[v]) < 0)          # the parity-shifted control
        out.append((v, lo_ok, hi_ok, sh_ok))
    return "CAP", out


def _hcoord(a, Lo, Hi, ops):
    """value -> the exact interface length carrying it."""
    return {a[i]: ops.sub(Hi[i], Lo[i]) for i in range(len(a))}


def g45(qmax=QMAX):
    """G4 (the invariant, exact) and G5 (the theorem), in one walk each."""
    r = {"rows": 0, "prongs": 0, "esc": 0, "other_end": {},
         "folds": 0, "inv_lo_ok": 0, "inv_hi_ok": 0, "ctrl_shift_ok": 0,
         "seq_matches_prediction": 0, "count_matches": 0,
         "rowmax_is_m_minus_1": 0, "prongs_to_101": 0, "bad": []}
    for P, Q, eps in tent_rows(qmax):
        _x, a, m = profile(P, Q, eps)
        if shape_of(a) != "tent":
            continue
        c = chain_of(P, Q, eps)
        if c is None:
            continue
        _ch, Lo, Hi, S, ops = c
        anchor, sign = _anchor(P, Lo, Hi, m)
        hc = _hcoord(a, Lo, Hi, ops)
        r["rows"] += 1
        rowmax = -1
        for j in range(m):
            end, seq = walk_levels(a, Lo, Hi, S, ops, j, anchor, sign, m, hc)
            if end is None:
                continue
            r["prongs"] += 1
            if Q <= 101:
                r["prongs_to_101"] += 1
            if end == "ESC":
                r["esc"] += 1
            else:
                r["other_end"][end] = r["other_end"].get(end, 0) + 1
            levels = [t[0] for t in seq]
            for _v, lo_ok, hi_ok, sh_ok in seq:
                r["folds"] += 1
                r["inv_lo_ok"] += int(lo_ok)
                r["inv_hi_ok"] += int(hi_ok)
                r["ctrl_shift_ok"] += int(sh_ok)
            pred = predicted(a, j)
            ok = (levels == pred)
            r["seq_matches_prediction"] += int(ok)
            r["count_matches"] += int(len(levels) == nu_of(a, j) // 2)
            rowmax = max(rowmax, len(levels))
            if not ok and len(r["bad"]) < 10:
                r["bad"].append({"P": P, "Q": Q, "eps": eps, "j": j,
                                 "got": levels[:12], "want": pred[:12], "end": end})
        r["rowmax_is_m_minus_1"] += int(rowmax == m - 1)
    return r


def g5_control(qmax=QMAX):
    """G5's negative control: the tent prediction must FAIL on the VALLEY rows."""
    r = {"valley_rows": 0, "valley_prongs": 0, "prediction_holds": 0, "zero_fold": 0}
    for Q in range(5, qmax + 1, 2):
        for P, eps in ((2, 1), (Q - 2, 0)):
            _x, a, m = profile(P, Q, eps)
            if shape_of(a) != "valley":
                continue
            c = chain_of(P, Q, eps)
            if c is None:
                continue
            _ch, Lo, Hi, S, ops = c
            hc = _hcoord(a, Lo, Hi, ops)
            anchor, sign = (Lo[0], +1) if ops.eq(Lo[0], Lo[m]) else (Hi[m], -1)
            r["valley_rows"] += 1
            for j in range(m):
                end, seq = walk_levels(a, Lo, Hi, S, ops, j, anchor, sign, m, hc)
                if end is None:
                    continue
                r["valley_prongs"] += 1
                r["zero_fold"] += int(len(seq) == 0)
                r["prediction_holds"] += int([t[0] for t in seq] == predicted(a, j))
    return r


# ------------------------------------------------------------------ G6

def g6(qmax=45):
    """THE ARM: does anything survive at turns >= 2?  Which clause dies first?"""
    r = {"rows": 0, "ap_step2_ok": 0, "residue_split_ok": 0, "two_ramps": 0,
         "prongs": 0, "inv_holds": 0, "inv_fails": 0, "folds": 0, "esc": 0,
         "seq_strictly_decreasing": 0, "seq_step_minus2": 0, "seq_scored": 0,
         "note": "turns==2 rows only; (i) AP, (ii) ramp residue split, (iii) I(nu), "
                 "(iv) the COORDINATE-FREE half: is the level sequence still -2?"}
    for P, Q, eps in all_rows(qmax):
        _x, a, m = profile(P, Q, eps)
        if min(P, Q - P) - 1 != 2:                    # turns == 2 exactly
            continue
        c = chain_of(P, Q, eps)
        if c is None:
            continue
        _ch, Lo, Hi, S, ops = c
        r["rows"] += 1
        sv = sorted(a)
        r["ap_step2_ok"] += int(all(sv[i + 1] - sv[i] == 2 for i in range(len(sv) - 1)))
        sgn = [1 if a[i + 1] > a[i] else -1 for i in range(m)]
        runs = 1 + sum(1 for i in range(m - 1) if sgn[i] != sgn[i + 1])
        r["two_ramps"] += int(runs == 2)
        # (ii) do the monotone runs still split the values by residue mod 4?
        segs, cur = [], [a[0]]
        for i in range(m):
            cur.append(a[i + 1])
            if i + 1 < m and sgn[i] != sgn[i + 1]:
                segs.append(cur)
                cur = [a[i + 1]]
        segs.append(cur)
        res = [{v % 4 for v in s} for s in segs]
        r["residue_split_ok"] += int(all(len(x) <= 1 for x in res))
        # (iii) the invariant on the real walk
        hc = _hcoord(a, Lo, Hi, ops)
        anchor, sign = (Lo[0], +1) if ops.eq(Lo[0], Lo[m]) else (Hi[m], -1)
        conc = (all(ops.eq(Lo[i], Lo[0]) for i in range(m + 1))
                or all(ops.eq(Hi[i], Hi[m]) for i in range(1, m + 1)))
        for j in range(m):
            end, seq = walk_levels(a, Lo, Hi, S, ops, j, anchor, sign, m, hc, cap=60_000)
            if end is None:
                continue
            r["prongs"] += 1
            r["folds"] += len(seq)
            r["esc"] += int(end == "ESC")
            if conc:
                for _v, lo_ok, hi_ok, _s in seq:
                    if lo_ok and hi_ok:
                        r["inv_holds"] += 1
                    else:
                        r["inv_fails"] += 1
            # (iv) the COORDINATE-FREE half -- needs no anchor, so it is scored on every row
            lv = [t[0] for t in seq]
            if len(lv) >= 2:
                r["seq_scored"] += 1
                r["seq_strictly_decreasing"] += int(all(lv[i + 1] < lv[i]
                                                        for i in range(len(lv) - 1)))
                r["seq_step_minus2"] += int(all(lv[i] - lv[i + 1] == 2
                                                for i in range(len(lv) - 1)))
        r.setdefault("concentric_rows", 0)
        r["concentric_rows"] += int(conc)
    return r


# ------------------------------------------------------------------ main

def main():
    args = [x for x in sys.argv[1:] if not x.startswith("--")]
    kw = dict(x[2:].split("=", 1) for x in sys.argv[1:] if x.startswith("--"))
    qmax = int(kw.get("qmax", QMAX))
    qbig = int(kw.get("qbig", QBIG))
    which = args[0] if args else "all"
    out = {}
    if which in ("all", "g0"):
        out["g0"] = g0(qbig)
    if which in ("all", "g1"):
        out["g1"] = g1(qmax)
    if which in ("all", "g2"):
        out["g2"] = g2(qbig)
    if which in ("all", "g3"):
        out["g3"] = g3(qbig)
    if which in ("all", "g4", "g5"):
        out["g45"] = g45(qmax)
        out["g5_control"] = g5_control(qmax)
    if which in ("all", "g6"):
        out["g6"] = g6(45)
    print(json.dumps(out, indent=1, default=str))
    if which == "all":
        with open(OUT, "w") as fh:
            json.dump(out, fh, indent=1, default=str)
        print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
