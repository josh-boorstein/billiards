#!/usr/bin/env python3
"""
s482_regular_leaves.py -- PRE-REGISTERED: queue item (s481-b), THE SURFACE-SIDE VERSION
OF [NCYL-340].

THE ITEM, from the s481 handoff: "(s481-b) the prong sum includes `R` prongs, which
[NCYL-338]/[NCYL-339] establish are NOT saddle connections -- so is there a cone-point-only
identity, and does it need [NCYL-339]'s merged cylinder list?  That is the item's original
quadratic-differential risk, now attackable."

⇒⇒ THE ANSWER, AND IT IS ONE RULE APPLIED THREE TIMES.  [NCYL-342] proved
`sum_c circ_c = sum_e ell_e` over the `Q` edges of `G`.  Which of those edges is a SADDLE
CONNECTION depends on WHICH SURFACE you are on, and an edge drops out of the sum exactly
when THE WHOLE CLOSED LEAF IT LIES ON is regular there -- i.e. when every endpoint visited
by that leaf is non-singular.  Endpoint orders on `B` are
`ord(Z_O) = P-2`, `ord(Z_A) = Q-P-2`, `ord(R-copy) = -1` ([NCYL-280]'s stratum
`Q(P-2, Q-P-2, -1^Q)`), so:
  * `B` **with its marked points** -- the diagram as [NCYL-280] builds it.  Nothing is
    regular; this is the (s481-a) level and `sum_c circ_c = sum_e ell_e`.
  * `B` **as a quadratic differential** -- an order-`0` endpoint is a MARKED point, not a
    cone point, so an edge with order-`0` at BOTH ends is a regular closed leaf.  Drop it,
    and the two cylinders it separated MERGE.
  * `Sigma` -- the orientation double cover.  A simple pole lifts to a REGULAR point, and
    so does an order-`0` point, so an edge with both ends in `{pole, order 0}` is regular
    upstairs.  Every non-regular edge has total preimage length `2·ell_e`.
⇒ Hence the cone-point-only identity the item asked for:

      sum_{cyl(Sigma)} circ  =  2 * ( sum_e ell_e  -  sum_{e regular on Sigma} ell_e )

and, one level down, `sum_{cyl(B)} circ = sum_e ell_e - sum_{e regular on B} ell_e`.
⚠ A regular closed leaf IS a boundary component of its cylinder, so its length is that
cylinder's circumference -- which is why dropping it subtracts exactly its length.

⚠⚠ THE RULE IS PER LEAF, NOT PER EDGE, AND THE PER-EDGE VERSION IS REFUTED -- I ran it
first and it MISSED on `21/242` rows (worst `1.18e-01`, every one a LONE row at
`P in {2, Q-2}`).  The counterexample is `2/5 lone`: the side-leaf `L2` runs `O`(order 0)
to `R`(pole), so BOTH its endpoints are regular on `Sigma` and the per-edge rule deletes
it -- but upstairs the leaf through it CONTINUES past both regular points into `H` and the
`A`-prong, which end on a genuine zero.  Its length is not deleted; it is absorbed into a
LONGER saddle connection.  An edge only leaves the sum when its whole leaf closes up
without meeting a singularity, which for a pole-pole fold it does and here it does not.
⇒ The boundary components are resolved by `s482_boundary_incidence`'s `leaf` id.

⇒⇒ AND THE ANSWER TO THE ITEM'S SECOND HALF IS *NO*, WITH A MECHANISM ATTACHED.  It does
not need [NCYL-339]'s merged list; it EXPLAINS it.  `ord(Z_O) = 0` iff `P = 2` and
`ord(Z_A) = 0` iff `P = Q-2`, and an `O`-prong exists only when the class has an inward
direction in the width-`2` window -- which is the PAIRED class.  So a regular edge on `B`
exists exactly on the `P in {2, Q-2}` PAIRED rows, which is [NCYL-336]'s family named from
the STRATUM rather than from a coincidence of circumferences.  ⇒ And it decides the arm
[NCYL-339] could not: `P in {1, Q-1}` has `ord(Z_O) = -1`, a pole, so there is NO regular
edge there and the store's list IS the surface list -- [NCYL-339] logged that arm as
"consistent, not forced".

HYPOTHESES (all able to fail)
  H1  A `B`-regular LEAF separates two DISTINCT cylinders, of EQUAL circumference, and its
      own length equals that circumference.  (It is a closed leaf and a whole boundary
      component, so all three follow -- and all three can be measured.)
  H2  The rows carrying a `B`-regular edge are EXACTLY [NCYL-336]'s `21`, and merging
      across those edges makes the moduli COMMENSURABLE on every one, where the unmerged
      list is commensurable on none.  A stratum selector reproducing an
      equal-circumference selector's family is the content.
  H3  On the `14` OTHER equal-circumference rows (`P in {1, Q-1}`) there is NO regular
      edge, hence no merge.
  H4  `sum_{cyl(Sigma)} circ == 2 * (sum_e ell_e - sum_{e on a Sigma-regular leaf} ell_e)`,
      with the `Sigma` list built by the lift rule: merge across `Sigma`-regular leaves
      joining DISTINCT cylinders; then a group carrying a `Sigma`-regular SELF-loop leaf
      lifts to ONE cylinder `(L, 2h)`, otherwise to TWO of `(L, h)`.  ⚠ A group carrying
      TWO self-loops would glue its lifts twice and close up into a torus component --
      impossible at genus `>= 2`, so it is counted and reported as a violation.
  H5  ⇒ THE CHECK THAT CAN FAIL AND IS NOT BOOKKEEPING: `#cyl(Sigma) <= g + s - 1`, the
      store's own `topological_bound` -- computed for the COVER and until now compared
      against a `B` count, which is [OPS-214]'s "wrong side of a cover".  Plus the area
      `sum circ*h == 2 * Q cot(alpha)`.
  C1  NEGATIVE CONTROL: merge by EQUAL CIRCUMFERENCE instead of by the stratum.  It picks
      up the `14` `P in {1,Q-1}` rows, where H3 says there is nothing to merge -- so the
      two selectors are NOT the same rule and the stratum one is the one with a reason.
  C2  NEGATIVE CONTROL: skip the merge.  Moduli must be incommensurable on the `21`.
  ⚠ NO BRANCH IS NOMINATED as the one I expect to fail ([OPS-044]).

⚠ SCOPE, and it is the honest limit.  `sum_{cyl(Sigma)} circ`'s left-hand side is
CONSTRUCTED from `B`'s data by the lift rule, not measured by a `Sigma`-side instrument --
the repo has none (`engine/cyl_diagram.py` scores an already-assembled diagram; it does not
build one from a centre).  So H4 is a DERIVATION whose inputs are measured, and what
carries falsifiable weight is H1/H2/H3/H5 and the controls.  Do not quote H4 as scored.

PRIOR ART: `rulings.py` on [NCYL-342], [NCYL-340], [NCYL-280], [NCYL-336], [NCYL-338],
[NCYL-339], [OPS-214]; `--grep 'marked point'`, `'orientation cover'`, `'stratum'`.
  -> [NCYL-336] owns the split's FAMILY (`P in {2,Q-2}`, paired, leg-swap-closed) and its
     arithmetic; [NCYL-338] owns that the store's cell list is nonetheless right;
     [NCYL-339] owns the parabolic exhibition and states the multiple `2` is NOT explained
     and that the `{1,Q-1}` arm is not decided.  None of the three reads the STRATUM.
  -> [NCYL-280] owns the stratum `Q(P-2, Q-P-2, -1^Q)` -- the input used here -- and
     [OPS-214] owns the "push to the genus-0 quotient" reading this inverts.
  -> [NCYL-342] owns the length identity this refines.

Usage:
  PYTHONPATH=.:engine:archive/scripts_2026-07:archive/scripts_2026-08:probes \
    .venv/bin/python3.13 probes/s482_regular_leaves.py [run | report]
"""
import collections
import glob
import json
import math
import sys
from fractions import Fraction

import s482_boundary_incidence as s482

OUT = "data/s482_regular_leaves.json"
STORE = "data/s439_exact_store"
INC = "data/s482_boundary_incidence.json"
SIDE_ENDS = {"L2": ("O", "R"), "L1": ("R", "A"), "H": ("O", "A")}


def orders(P, Q):
    """Singularity orders of `B = Q(P-2, Q-P-2, -1^Q)` ([NCYL-280]), by vertex."""
    return {"O": P - 2, "A": Q - P - 2, "R": -1}


def commensurable(mus, maxden=64, tol=1e-9):
    """Are all moduli rational multiples of a common one?  Practical form: every ratio
    `mu/mu_min` is within `tol` (relative) of a fraction with denominator `<= maxden`.
    ⚠ Not vacuous: `1/cos(pi/11) = 1.042217...` is not within `1e-9` of any such fraction,
    which is the case this has to reject ([NCYL-336])."""
    m0 = min(mus)
    for m in mus:
        r = m / m0
        f = Fraction(r).limit_denominator(maxden)
        if abs(r - float(f)) > tol * max(1.0, r):
            return False
    return True


def do_row(P, Q, role, rec, cl, row):
    cyls = cl["cylinders"]
    edges, _, _ = s482.edge_list(P, Q, row["parity"], rec["prong_cap"])
    ordv = orders(P, Q)

    # ---- endpoints, MEASURED (trace_prong's arrival vertex), not assumed via (*)
    meta, n_offdiag = [], 0
    for e in edges:
        if e["kind"] == "prong":
            arr = e["why"].split(":")[1] if e["why"].startswith("vertex:") else None
            ends = (e["v0"], arr)
            if arr != e["v0"]:
                n_offdiag += 1
        else:
            ends = SIDE_ENDS[e["side"]]
        o = [ordv[v] if v else None for v in ends]
        meta.append({"ends": ends, "ell": e["ell"],
                     "regB": all(x == 0 for x in o),
                     "regS": all(x in (0, -1) for x in o)})

    # ---- boundary components, and the LEAVES (a leaf = the two components an edge joins)
    side_of = {}                       # (edge, side) -> component key
    comp = {}                          # component key -> [edge, ...]
    for e, si, c, lf in row["incidence"]:
        k = (c, lf)
        side_of[(e, si)] = k
        comp.setdefault(k, []).append(e)
    leaf = {}                          # component key -> the paired component key
    for ei in range(len(edges)):
        k0, k1 = side_of[(ei, 0)], side_of[(ei, 1)]
        leaf[k0], leaf[k1] = k1, k0

    def leaf_regular(k, field):
        return all(meta[e][field] for e in comp[k]) and all(
            meta[e][field] for e in comp[leaf[k]])

    regB_leaf = {k for k in comp if leaf_regular(k, "regB")}
    regS_leaf = {k for k in comp if leaf_regular(k, "regS")}
    deadB = {e for k in regB_leaf for e in comp[k]}
    deadS = {e for k in regS_leaf for e in comp[k]}

    # ---- H1 + the B-merge, across B-regular leaves
    par = list(range(len(cyls)))

    def find(x):
        while par[x] != x:
            par[x] = par[par[x]]
            x = par[x]
        return x

    h1, seen = [], set()
    for k in regB_leaf:
        k2 = leaf[k]
        if (k2, k) in seen:
            continue
        seen.add((k, k2))
        a, b = k[0], k2[0]
        La, Lb = cyls[a]["circumference"], cyls[b]["circumference"]
        h1.append({"distinct": a != b, "eq_circ": abs(La - Lb) / La,
                   "len_vs_circ": abs(sum(meta[e]["ell"] for e in comp[k]) - La) / La})
        par[find(a)] = find(b)
    grp = collections.defaultdict(list)
    for i in range(len(cyls)):
        grp[find(i)].append(i)
    merged = {g: {"L": cyls[v[0]]["circumference"],
                  "h": sum(cyls[i]["height"] for i in v), "members": v}
              for g, v in grp.items()}

    # ---- the Sigma list: merge across Sigma-regular leaves joining DISTINCT cylinders,
    # then glue the two lifts of any group carrying a Sigma-regular SELF-loop.
    parS = list(range(len(cyls)))

    def findS(x):
        while parS[x] != x:
            parS[x] = parS[parS[x]]
            x = parS[x]
        return x
    selfloop = collections.Counter()
    seen = set()
    for k in regS_leaf:
        k2 = leaf[k]
        if (k2, k) in seen:
            continue
        seen.add((k, k2))
        if k[0] != k2[0]:
            parS[findS(k[0])] = findS(k2[0])
        else:
            selfloop[k[0]] += 1
    grpS = collections.defaultdict(list)
    for i in range(len(cyls)):
        grpS[findS(i)].append(i)
    sig, n_double_glue = [], 0
    for g, v in grpS.items():
        L = cyls[v[0]]["circumference"]
        h = sum(cyls[i]["height"] for i in v)
        nsl = sum(selfloop[i] for i in v)
        if nsl > 1:
            n_double_glue += 1
        if nsl:
            sig.append({"L": L, "h": 2 * h})
        else:
            sig += [{"L": L, "h": h}, {"L": L, "h": h}]

    sum_all = sum(m["ell"] for m in meta)
    circ_sig = sum(c["L"] for c in sig)
    pred_sig = 2 * (sum_all - sum(meta[e]["ell"] for e in deadS))
    circ_B = sum(c["L"] for c in merged.values())
    pred_B = sum_all - sum(meta[e]["ell"] for e in deadB)
    area_sig = sum(c["L"] * c["h"] for c in sig)
    area_want = 2 * Q / math.tan((P / Q) * (math.pi / 2))
    cc = collections.Counter(round(c["circumference"], 9) for c in cyls)

    return {
        "P": P, "Q": Q, "role": role, "C": len(cyls),
        "n_offdiag_prongs": n_offdiag,
        "n_regB_leaves": len(regB_leaf) // 2 + sum(1 for k in regB_leaf if leaf[k] == k),
        "n_deadB": len(deadB), "n_deadS": len(deadS),
        "h1": h1,
        "has_eq_circ": any(v > 1 for v in cc.values()),
        "C_B": len(merged), "C_sig": len(sig), "n_double_glue": n_double_glue,
        "mu_before_ok": commensurable([c["height"] / c["circumference"] for c in cyls]),
        "mu_after_ok": commensurable([c["h"] / c["L"] for c in merged.values()]),
        "h3_resid": abs(circ_B - pred_B) / pred_B,
        "h4_resid": abs(circ_sig - pred_sig) / pred_sig,
        "h5_bound": cl.get("topological_bound"),
        "h5_ok": (cl.get("topological_bound") is None
                  or len(sig) <= cl["topological_bound"]),
        "area_resid": abs(area_sig - area_want) / area_want,
        "sum_all": sum_all, "circ_sig": circ_sig, "circ_B": circ_B,
    }


def run():
    inc = json.load(open(INC))
    out = {}
    for key, row in inc.items():
        rec = json.load(open(f"{STORE}/{row['P']}_{row['Q']}.json"))
        out[key] = do_row(row["P"], row["Q"], row["role"], rec,
                          rec["classes"][row["role"]], row)
    json.dump(out, open(OUT, "w"))
    report()


def report():
    d = json.load(open(OUT))
    rows = list(d.values())
    n = len(rows)
    regB = [r for r in rows if r["n_deadB"]]
    eq = [r for r in rows if r["has_eq_circ"]]
    print(f"ROWS {n}    off-diagonal prongs (a prong not returning to its own vertex): "
          f"{sum(r['n_offdiag_prongs'] for r in rows)}")
    h1 = [x for r in rows for x in r["h1"]]
    print(f"\nH1  B-regular LEAVES {len(h1)} on {len(regB)} rows:  separate DISTINCT "
          f"cylinders {sum(x['distinct'] for x in h1)}/{len(h1)}   equal circumference worst "
          f"{max((x['eq_circ'] for x in h1), default=0):.2e}   leaf length == circ worst "
          f"{max((x['len_vs_circ'] for x in h1), default=0):.2e}")
    print(f"\nH2  rows with a B-regular edge {len(regB)}   "
          f"P in {{2,Q-2}} & paired: "
          f"{sum(1 for r in regB if r['P'] in (2, r['Q'] - 2) and r['role'] == 'paired')}"
          f"   moduli commensurable BEFORE merge {sum(r['mu_before_ok'] for r in regB)}"
          f"/{len(regB)}  AFTER {sum(r['mu_after_ok'] for r in regB)}/{len(regB)}"
          "   [C2 is the BEFORE column]")
    only_eq = [r for r in eq if not r["n_deadB"]]
    print(f"\nH3  equal-circumference rows {len(eq)};  of them WITHOUT a B-regular edge "
          f"{len(only_eq)}  -> P in {{1,Q-1}}: "
          f"{sum(1 for r in only_eq if r['P'] in (1, r['Q'] - 1))}/{len(only_eq)} "
          "(no merge; the store list IS the surface list there)")
    print(f"C1  equal-circumference SELECTOR would merge on {len(eq)} rows, the stratum "
          f"selector on {len(regB)} -- they differ on {len(only_eq)}")
    print(f"\nH3b `sum_cyl(B) circ == sum_e ell - sum_(B-dead) ell`  worst rel "
          f"{max(r['h3_resid'] for r in rows):.2e}   [B unmarked, list CONSTRUCTED]")
    print(f"H4  `sum_cyl(Sigma) circ == 2(sum_e ell - sum_Sigma-reg ell)`  worst rel "
          f"{max(r['h4_resid'] for r in rows):.2e}  ⚠ DERIVED, not independently measured")
    print(f"    double-glued groups (would be a torus component; must be 0): "
          f"{sum(r['n_double_glue'] for r in rows)}")
    print(f"H5  #cyl(Sigma) <= store topological_bound  {sum(r['h5_ok'] for r in rows)}/{n}"
          f"   (equality on {sum(1 for r in rows if r['C_sig'] == r['h5_bound'])})"
          f"   area worst {max(r['area_resid'] for r in rows):.2e}")
    qs = sorted({r["Q"] for r in rows})
    print(f"\nCOVERAGE Q = {qs[0]}..{qs[-1]}, {n} certified full-mode classes")


if __name__ == "__main__":
    if sys.argv[1:] and sys.argv[1] == "report":
        report()
    else:
        run()
