#!/usr/bin/env python3
"""
s472_legswap_existence.py -- PRE-REGISTERED.  Queue item (s447-a).

WHAT THIS IS.  A PROOF of the ONE scope caveat on `necklace_paper.md` Cor. 8.2:
`s >= 1` on the EVEN-`P` `{L2,H}` class, which [NCYL-285] left MEASURED (`104/104`,
`Q <= 32`) because the orphan theorems it quotes reach only `P` ODD.  Plus the
instrument that scores the proof's CLASS BOOKKEEPING and its transported OBJECT --
not its conclusion.

⚠ The proof needs NO new machinery.  It is a COMPOSITION of five already-PROVED
results, and the reason it sat open is recorded at the bottom ("WHY THIS WAS MISSED").

------------------------------------------------------------------------------------
THE PROOF (route 1 -- transport the OBJECT).  Let `Q` be odd, `P` even, `gcd(P,Q)=1`.

(E1) THE LEG SWAP IS A RELABELLING, NOT A TRANSPORT.  `T(P/Q)` and `T((Q-P)/Q)` are
     the SAME triangle with the two legs interchanged and `H` fixed ([NCYL-068]);
     hence the `⊥L2` beam of `T(P/Q)` IS the `⊥L1` beam of `T((Q-P)/Q)` and the two
     `⊥H` beams coincide -- "pure relabelling, no grid condition, valid at irrational
     `α`" ([NCYL-098], PROVED).
(E2) `Q` odd and `P` even give `Q-P` ODD, so Theorem D ([T-ORPHAN], PROVED all `Q`)
     gives `T((Q-P)/Q)`'s `⊥L1` beam a unique ORPHAN.
(E3) `Q` is odd, so [NCYL-072] (PROVED; the vertex gap DISCHARGED at s332 by
     [NCYL-132]/[NCYL-133]) puts that orphan's far turnaround perpendicular hit on `H`.
(E4) So its core geodesic is a DOUBLE NORMAL joining `L1` and `H`, and by [NCYL-080]'s
     time-reversal reading (PROVED) the cylinder carrying it is met by BOTH `⊥L1` and
     `⊥H` -- a STRADDLER of the `{L1,H}` class at `(Q-P)/Q`.
(E5) By (E1) that cylinder is met by `⊥L2` and `⊥H` at `P/Q`.  The classes match:
     the paired class of `P/Q` is `eps=1` with verticals `{L2,H}` and that of
     `(Q-P)/Q` is `eps=0` with verticals `{L1,H}` (H0 below, `4141/4141`).  ⇒ `s >= 1`.
(E6) With Theorem 8.1's `s <= 1` ([NCYL-285], PROVED in all three classes): **`s = 1`**,
     and `overlap = {the double normal}` on the even-`P` class too.  ∎

⇒ `s >= 1` at even `P` is UNCONDITIONAL -- it needs no CP, exactly like the odd-`P`
half.  CP enters only through `s <= 1`, which is where it already was.

ROUTE 2 (independent of (E2)-(E4); shares only (E1)).  A free lemma:

  PARITY LEMMA.  In a completely periodic class, for each vertical arc `sigma`,
  `n_sigma = 2*a_sigma + s` where `a_sigma` = #cylinders with BOTH `Fix`-segments on
  `sigma`'s arc.  (Every cylinder meeting `Fix` meets it in exactly two vertical
  segments -- Thm 8.1 step (1); the segments on `sigma`'s arc are exactly the `⊥sigma`
  beam's cells -- [NCYL-283] (4).)  Hence `n_sigma == s (mod 2)`.

  Theorem D's PARITY corollary (`n(P/Q)` odd iff `P` odd -- PROVED all `Q`, and the
  clause SURVIVED the s272 correction, which struck only the even-`P` VALUE `n = Q-1`)
  plus (E1) give `n_L2(P/Q) = n_L1((Q-P)/Q)` ODD.  So `s` is odd, and `s <= 1` gives
  `s = 1`.  ∎

------------------------------------------------------------------------------------
WHAT IS PRE-REGISTERED.  The proof above is not what this file scores -- a proof's
sweeps are guards on its joints, not its grounds ([NCYL-318] scope (i) / [OPS-041]).
The joint I can actually get wrong is the CLASS BOOKKEEPING of (E5) and the claim in
(E1) that the two beams are the SAME beam rather than two beams of equal size.

  H0  ARITHMETIC, whole box.  For every coprime `(P,Q)`, `Q` odd, `P` even: the paired
      class of `P/Q` is `eps=1` with `ver = {L2,H}`, and that of `(Q-P)/Q` is `eps=0`
      with `ver = {L1,H}`.  (Fails ⇒ (E5) is wrong.)
  H1  SET LEVEL, and this is the real content.  The `⊥L2` cell decomposition at `P/Q`
      and the `⊥L1` decomposition at `(Q-P)/Q` are the SAME cylinders: their canonical
      words agree as multisets after the relabelling `L1 <-> L2`.  Same for the two
      `⊥H` decompositions.  ⚠ Compared as CYCLIC words up to reversal (the traces have
      no shared orientation convention); the forward-only rate is reported separately
      so the reversal is visible, not hidden.
  H2  THE THEOREM.  The straddler at `P/Q` relabels to the straddler at `(Q-P)/Q`, and
      the latter carries [NCYL-076]'s odd-letter ORPHAN signature.  This is (E5)'s
      conclusion read as an identification of a NAMED cylinder, not a count.
  H3  THE PARITY LEMMA, EXACTLY (`n_sigma == 2*a_sigma + s`), not mod 2, on all three
      classes.  ⚠⚠ See NEGATIVE CONTROLS: the mod-2 form is nearly vacuous here.
  H4  ROUTE 2's prediction: `n_L2` is ODD on every even-`P` odd-`Q` row.

NEGATIVE CONTROLS -- every one must FAIL, and each kills a different way of being
right by accident:
  N1  the SAME-eps partner `(Q-P, Q, eps=1)`: kills "any partner row would match".
  N2  a NON-partner row `(P', Q)`, `P' not in {P, Q-P}`, same `Q`: kills "the word
      multiset is coarse at fixed `Q`".
  N3  the relabelling OMITTED (raw canonical words): kills "the match is insensitive
      to the leg swap", i.e. that H1 tests the swap at all.
  ⚠⚠ AND THE ONE THAT IS NOT A CONTROL BUT A WARNING, stated because it would
  otherwise read as support: `s == 1` on all `312` scored rows of `data/s447_overlap.json`,
  so ANY test of the form "`s == f(...) (mod 2)`" has a CONSTANT left-hand side and
  could only fail by some `n_sigma` coming back EVEN.  That is why H3 is scored in its
  EXACT form with the `a_sigma` distribution printed: `n_sigma = 2*a_sigma + s` can
  fail on any row where `a_sigma` varies, and it does vary.

PRIOR ART: read `rulings.md` [NCYL-068]/[NCYL-072]/[NCYL-076]/[NCYL-080]/[NCYL-092]/
[NCYL-093]/[NCYL-098]/[NCYL-112]/[NCYL-113]/[NCYL-114]/[NCYL-132]/[NCYL-136]/
[NCYL-237]/[NCYL-250]/[NCYL-282]/[NCYL-283]/[NCYL-285]/[NCYL-286]/[NCYL-287];
`rulings.py --grep 'leg swap' / 'double normal' / 'parity of n' / 'n_L2 odd' / 'even P';
grepped `claims.md` [T-ORPHAN], `orphan_theorem.md`, `paper.md` Thm D ->
  - [NCYL-113] IS the measured statement being proved: "EVERY CENTRE HAS EXACTLY ONE
    DOUBLE NORMAL ... `L2` <-> `P+Q` odd, `H` <-> `Q` odd", census `60/60`, `Q = 4..14`.
    It is the target, not an input.
  - [NCYL-112] is the near miss and it says so in its own words: the direction-index
    parity decides WHICH pair can carry a double normal, and is "**NECESSARY, NOT
    SUFFICIENT** -- so existence is measured ([NCYL-113]), not deduced".  True of the
    parity route; the leg swap is a different route and it does deduce it.
  - [NCYL-098] states (E1) verbatim and draws the `N_L2`-is-a-corollary consequence for
    the COUNTING function.  Nobody applied it to the OVERLAP.
  - [NCYL-237] upgrades the leg swap to an identity of SETS -- but at EVEN `Q` only, and
    it records exactly why the even-`P` rows were never tested: "at odd `Q` the leg-swap
    partner has even `P` and is not in the census at all".  H1 is its odd-`Q` twin.
  - [NCYL-285] is the `s <= 1` half (E6) and the source of the caveat; its own text
    already names the missing piece -- "Its proof, if wanted, is the missing piece there".
  - NO probe in the repo composes the leg swap with Theorem D, and no ruling states
    `s >= 1` at even `P` other than as MEASURED.

WHY THIS WAS MISSED (for the record, and it is [OPS-246]'s shape again).  Every
ingredient predates s447 by 100+ sessions.  Three things kept them apart: (i) the
even-`P` rows are excluded from the odd-`Q` censuses by construction ([NCYL-237]), so
the leg swap was never exercised in the direction that needed it; (ii) [NCYL-112]'s
correct "existence is measured, not deduced" is about the PARITY route and reads as a
statement about existence in general; and (iii) the queue item itself named the wrong
technique -- "the technique is in hand, [NCYL-285]'s four-right-angle disk argument" --
which is the `s <= 1` machine and cannot prove an existence statement at all, exactly
as [NCYL-287] says ("the overlap half is a FINITENESS statement ... the covering half
is an EXISTENCE statement").  The disk never enters this proof.

WHAT THIS DOES NOT DO.  It does not touch (G0b), `p = 0`, the covering law or anything
on the gamma=1 DAG; it does not give a method for non-CP; and it does not prove
UNIQUENESS of the double normal independently -- uniqueness is (E6), i.e. [NCYL-285].

Usage:
  PYTHONPATH=.:engine:archive/scripts_2026-07:archive/scripts_2026-08:probes \\
    .venv/bin/python3.13 probes/s472_legswap_existence.py h0 [--big=201]
  ... s472_legswap_existence.py sweep [--qmax=20] [--budget=1800]
  ... s472_legswap_existence.py all
"""
import json
import math
import os
import sys
import time
from collections import Counter

from right_triangle_billiards import RightTriangleBilliard
import s315_cylinder_count as s315
import s437_oblique_cylinders as s437
import s439_exact_cells as s439
import s446_covering_identity as s446
import s447_overlap_disk as s447

CAP = 20000
CAP_BIG = 200000
OUT = "data/s472_legswap.json"

SWAP = {"L1": "L2", "L2": "L1", "H": "H"}


def relabel(word, swap=True):
    """The canonical CYCLIC word after the leg relabelling `L1 <-> L2` (H fixed).

    ⚠ `s315.canon` is the lexicographically smallest rotation, and relabelling changes
    the lex order -- so the word must be RE-canonicalised after the swap, never before.
    """
    w = tuple(SWAP[x] for x in word) if swap else tuple(word)
    return s315.canon(w)


def cyc_key(word, swap=True):
    """`relabel` up to REVERSAL as well: the traces carry no shared orientation."""
    f = relabel(word, swap)
    r = s315.canon(tuple(reversed(f)))
    return min(f, r)


# ------------------------------------------------------------------ H0 (arithmetic)

def h0(big=201, verbose=True):
    """Whole box: paired(P,Q) = (eps 1, {L2,H}) and paired(Q-P,Q) = (eps 0, {L1,H})."""
    bad, n = [], 0
    for Q in range(3, big + 1, 2):
        for P in range(2, Q, 2):
            if math.gcd(P, Q) != 1:
                continue
            n += 1
            e1 = [e for e in (0, 1) if len(s446.hv_sides(P, Q, e)[1]) == 2]
            e2 = [e for e in (0, 1) if len(s446.hv_sides(Q - P, Q, e)[1]) == 2]
            if len(e1) != 1 or len(e2) != 1:
                bad.append((P, Q, "not-unique-paired-class"))
                continue
            v1 = sorted(s446.hv_sides(P, Q, e1[0])[1])
            v2 = sorted(s446.hv_sides(Q - P, Q, e2[0])[1])
            if not (e1[0] == 1 and v1 == ["H", "L2"] and e2[0] == 0 and v2 == ["H", "L1"]):
                bad.append((P, Q, e1[0], v1, e2[0], v2))
    if verbose:
        print(f"  H0  paired(P,Q)=(eps1,{{L2,H}}) & paired(Q-P,Q)=(eps0,{{L1,H}}) : "
              f"{n - len(bad)}/{n} rows (Q odd <= {big}, P even coprime)", flush=True)
        for b in bad[:6]:
            print("     ", b, flush=True)
    return bad, n


# ------------------------------------------------------------------ the cells

def cells_of(P, Q, cap=CAP, cap2=CAP_BIG):
    """`(eps, ver, {side: [cell,...]})` for the PAIRED class of `P/Q`, or `(status,)`.

    Same construction as `s447_overlap_disk.row` -- deliberately, so H1/H2 read the
    same objects [NCYL-285] scored and no second cell convention enters the repo.
    """
    eps = next(e for e in (0, 1) if len(s446.hv_sides(P, Q, e)[1]) == 2)
    hor, ver = s446.hv_sides(P, Q, eps)
    bmap, reasons, _nb, _n1, aborted = s439.boundary_map(
        P, Q, eps, cap=cap, cap2=cap2, abort_on_open=True)
    if aborted or any(not r["reason"].startswith("vertex") for r in reasons.values()):
        return "PRONG_CAPPED", eps, ver, None
    geo = s437.geometry(P, Q)
    B = RightTriangleBilliard((P / Q) * (math.pi / 2))
    cells = {}
    for side in ver:
        u = s437.perp_index(side, P, Q)
        cl, _mg = s439.strip_cells(B, geo, Q, side, u, bmap[(side, u)], cap)
        if any(c["canon"] is None for c in cl):
            return "CELL_OPEN", eps, ver, None
        cells[side] = cl
    return "ok", eps, ver, cells


def straddler_and_counts(ver, cells, P, Q):
    """`(straddler_label, s, {side: (n_sigma, a_sigma)})` on `s447`'s arc reading."""
    reading = s447.arc_reading(ver, P, Q)
    seq = []
    for side, rev in reading:
        lab = [str(c["canon"]) for c in cells[side]]
        seq += lab[::-1] if rev else lab
    n_first = len(cells[reading[0][0]])
    strad = s447._straddlers(seq, n_first)
    per = {}
    for side in ver:
        labs = [str(c["canon"]) for c in cells[side]]
        cnt = Counter(labs)
        per[side] = (len(labs), sum(1 for k, v in cnt.items() if v == 2))
    return (strad[0] if len(strad) == 1 else None), len(strad), per, seq


# ------------------------------------------------------------------ the sweep

def sweep(qmax=20, out_path=OUT, budget_s=None, verbose=True):
    """H1-H4 + N1-N3 on every even-`P` odd-`Q` row and its leg-swap partner.

    ROW-LEVEL RESUME (the s458 rule, level 1): completed rows are keyed `P/Q` in
    `out_path` and skipped, so raising `--qmax` costs only the new rows.
    """
    t0 = time.time()
    store = {}
    if os.path.exists(out_path):
        store = json.load(open(out_path)).get("rows", {})
    cache = {}

    def get(p, q):
        if (p, q) not in cache:
            cache[(p, q)] = cells_of(p, q)
        return cache[(p, q)]

    targets = [(P, Q) for Q in range(3, qmax + 1, 2)
               for P in range(2, Q, 2) if math.gcd(P, Q) == 1]
    for P, Q in targets:
        key = f"{P}/{Q}"
        if key in store:
            continue
        if budget_s and time.time() - t0 > budget_s:
            if verbose:
                print(f"  [budget {budget_s}s reached at {key}; "
                      f"{len(store)}/{len(targets)} rows done -- rerun to resume]",
                      flush=True)
            break
        st_a, eps_a, ver_a, cl_a = get(P, Q)
        st_b, eps_b, ver_b, cl_b = get(Q - P, Q)
        rec = {"P": P, "Q": Q, "eps": eps_a, "eps_partner": eps_b,
               "ver": sorted(ver_a), "ver_partner": sorted(ver_b),
               "status": st_a if st_a != "ok" else st_b}
        if st_a != "ok" or st_b != "ok":
            store[key] = rec
            continue

        # ---- H1: the two decompositions are the SAME cylinders, after L1<->L2.
        def multiset(cells, side, swap, rev_ok):
            f = cyc_key if rev_ok else relabel
            return sorted(f(c["canon"], swap) for c in cells[side])

        rec["H1_L"] = (multiset(cl_a, "L2", True, True)
                       == multiset(cl_b, "L1", False, True))
        rec["H1_H"] = (multiset(cl_a, "H", True, True)
                       == multiset(cl_b, "H", False, True))
        rec["H1_L_fwd"] = (multiset(cl_a, "L2", True, False)
                           == multiset(cl_b, "L1", False, False))
        rec["H1_H_fwd"] = (multiset(cl_a, "H", True, False)
                           == multiset(cl_b, "H", False, False))
        rec["H1_ok"] = rec["H1_L"] and rec["H1_H"]

        # ---- H2: the straddler transports, and the partner's is the orphan.
        sa, s_a, per_a, _ = straddler_and_counts(ver_a, cl_a, P, Q)
        sb, s_b, per_b, _ = straddler_and_counts(ver_b, cl_b, Q - P, Q)
        rec["s"], rec["s_partner"] = s_a, s_b
        if sa is not None and sb is not None:
            wa = cyc_key(tuple(s447._canon_letters(sa)), True)
            wb = cyc_key(tuple(s447._canon_letters(sb)), False)
            rec["H2_same_cylinder"] = (wa == wb)
            rec["H2_partner_is_orphan"] = s447._hinges(
                s447._canon_letters(sb), ver_b)
            rec["H2_ok"] = rec["H2_same_cylinder"] and rec["H2_partner_is_orphan"]
        else:
            rec["H2_ok"] = False

        # ---- H3: the parity lemma, EXACTLY.
        rec["n_sigma"] = {k: v[0] for k, v in per_a.items()}
        rec["a_sigma"] = {k: v[1] for k, v in per_a.items()}
        rec["H3_ok"] = all(n == 2 * a + s_a for n, a in per_a.values())
        rec["H3_ok_partner"] = all(n == 2 * a + s_b for n, a in per_b.values())

        # ---- H4: route 2's prediction.
        rec["H4_ok"] = per_a["L2"][0] % 2 == 1

        # ---- NOTE (not a control -- there is nothing to score).  The eps -> 1-eps
        #      pairing of (E5) is FORCED, not chosen: the same-eps class of the partner
        #      is the LONE one (v == 1) and carries no {L1,H} reading at all.
        e_same = 1 - eps_b
        rec["same_eps_ver"] = sorted(s446.hv_sides(Q - P, Q, e_same)[1])

        # ---- N2: a NON-partner row at the same Q.
        rec["N2"] = None
        for Pp in range(1, Q):
            if Pp in (P, Q - P) or math.gcd(Pp, Q) != 1:
                continue
            st_c, _e, ver_c, cl_c = get(Pp, Q)
            if st_c != "ok" or "L1" not in ver_c:
                continue
            rec["N2"] = (multiset(cl_a, "L2", True, True)
                         == multiset(cl_c, "L1", False, True))
            rec["N2_row"] = f"{Pp}/{Q}"
            break

        # ---- N3: the relabelling omitted.
        rec["N3"] = (sorted(cyc_key(c["canon"], False) for c in cl_a["L2"])
                     == sorted(cyc_key(c["canon"], False) for c in cl_b["L1"]))
        store[key] = rec
        if verbose:
            print(f"  {key:>8}  H1={rec['H1_ok']} H2={rec['H2_ok']} H3={rec['H3_ok']} "
                  f"H4={rec['H4_ok']}  n_L2={per_a['L2'][0]} s={s_a}  "
                  f"N2={rec['N2']} N3={rec['N3']}", flush=True)

    summary = report(store, verbose=verbose)
    json.dump({"qmax": qmax, "summary": summary, "rows": store},
              open(out_path, "w"), indent=1)
    if verbose:
        print(f"  [{time.time() - t0:.0f}s -> {out_path}]", flush=True)
    return store, summary


def report(store, verbose=True):
    ok = [r for r in store.values() if r.get("status") == "ok" or "H1_ok" in r]
    scored = [r for r in ok if "H1_ok" in r]
    n = len(scored)
    s = {
        "rows": len(store), "scored": n,
        "H1": sum(r["H1_ok"] for r in scored),
        "H1_fwd_only": sum(r["H1_L_fwd"] and r["H1_H_fwd"] for r in scored),
        "H2": sum(r["H2_ok"] for r in scored),
        "H3": sum(r["H3_ok"] for r in scored),
        "H3_partner": sum(r["H3_ok_partner"] for r in scored),
        "H4": sum(r["H4_ok"] for r in scored),
        "N2_fired": sum(1 for r in scored if r.get("N2") is False),
        "N2_scored": sum(1 for r in scored if r.get("N2") is not None),
        "N3_fired": sum(1 for r in scored if r["N3"] is False),
        "s_values": dict(Counter(r["s"] for r in scored)),
        "a_sigma_spread": dict(Counter(
            a for r in scored for a in r["a_sigma"].values())),
        "capped": sum(1 for r in store.values() if r.get("status") not in ("ok", None)),
    }
    if verbose:
        print(f"\n  SCORED {n} even-P odd-Q rows (+ partners)")
        print(f"    H1  same cylinders after L1<->L2 (cyclic, up to reversal) : "
              f"{s['H1']}/{n}      [forward-only: {s['H1_fwd_only']}/{n}]")
        print(f"    H2  straddler transports AND partner's is the orphan     : {s['H2']}/{n}")
        print(f"    H3  parity lemma n_sigma == 2*a_sigma + s (exact)        : {s['H3']}/{n}"
              f"   (partner rows {s['H3_partner']}/{n})")
        print(f"    H4  n_L2 ODD  (route 2's prediction)                     : {s['H4']}/{n}")
        print(f"    N2  non-partner row at same Q          MUST FAIL, fired  : "
              f"{s['N2_fired']}/{s['N2_scored']}")
        print(f"    N3  relabelling omitted                MUST FAIL, fired  : "
              f"{s['N3_fired']}/{n}")
        print(f"    a_sigma spread (H3 is not vacuous iff this varies): {s['a_sigma_spread']}")
        print(f"    s values: {s['s_values']}   capped/open rows: {s['capped']}")
    return s


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    kw = dict(a[2:].split("=", 1) for a in sys.argv[1:] if a.startswith("--"))
    mode = args[0] if args else "all"
    if mode in ("h0", "all"):
        h0(big=int(kw.get("big", 201)))
    if mode in ("sweep", "all"):
        sweep(qmax=int(kw.get("qmax", 20)),
              budget_s=float(kw["budget"]) if "budget" in kw else None)


if __name__ == "__main__":
    main()
