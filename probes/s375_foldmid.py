"""s375 PRE-REGISTERED -- does the limit orbit at a `J`-fold boundary REVERSE AT THE
`R` CORNER, at its own temporal midpoint?  The converse of [NCYL-216], as a prediction.

WHY.  s374 (`orphan_theorem.md` §B′.2) proved [NCYL-216]'s FORWARD direction and reduced
the converse to (R★) = "`T` is discontinuous at every non-`R` boundary", then measured
(R★) and could not explain it.  Working the time-reversal instead of the discontinuity
turns the converse into a statement with a falsifiable geometric prediction:

  Let `y` be an interior boundary whose flanking cells form a `J`-pair.  Then both are
  case (i) with the SAME `b`, and `b = 2y` (the left cell maps onto the right one only if
  `b − y = y`).  So `T(y⁻) = T(y⁺) = y`: both one-sided limit orbits run from the launch
  point `p` back to `p`, arriving perpendicular to `L1`.  Time reversal sends the orbit at
  `s` to the orbit at `T(s) = 2y − s`, so `rev(γ⁻) = γ⁺`.  If the two limits are the SAME
  path then that path is self-reverse, so its velocity reverses at its temporal midpoint
  `τ/2`; by Lemma 2 case (i) a case-(i) orbit has NO interior head-on, and the limit adds
  none except at a wall ENDPOINT, so the midpoint reversal is not a wall hit but a VERTEX
  -- and by [NCYL-133] the only vertex that reverses a velocity is `R`.  ⇒ `y` is an
  `R`-graze.  ⇒ (R★) was never the residual the converse needs.

  ⇒ PREDICTION H1, which is what this probe scores.  At an `R`-boundary the limit orbit is
  PALINDROMIC IN POSITION about a single centre, and that centre is a near-`R` corner
  cluster.  The orbit bounces `p → R → p`: it retraces via the CORNER, where the orphan
  retraces via a head-on WALL.  H1 can fail outright -- if the reversal turned out to be a
  head-on on `L1` at the palindrome centre, the proof above would be wrong.

  ⇒ H2 is the one REMAINING residual, (V2).  The step "if the two limits are the same
  path" is discharged exactly when the limit orbit meets no CONE point (a regular vertex
  does not split a geodesic; a cone point does).  So: does the limit orbit at a fold meet
  only REGULAR vertices?  Only `R` is regular in general -- `A` is regular iff `r ∈ {1,2}`
  and `O` iff `P ∈ {1,2}` -- so H2 reads "every near-vertex cluster is at `R`".

⚠ THE `r = 2` ROWS ARE NOT WHAT [NCYL-219] TOOK THEM FOR, and H4 re-reads them.  s374's
named risk set (a graze at a REGULAR non-`R` vertex) is realised only by 6 `A`-boundaries
at `r = 2`, and ALL SIX ARE ORPHAN ENDPOINTS (left cell `is_orphan`, `chi == 1`; stored
`T_lo = y` exactly, `T_hi = 1.0` exactly).  There the graze is not a passage BEFORE the
reversal -- the reversal happens INSIDE the corner cluster at `A`.  Traced at `(3,5)`:
both sides enter the same cluster at `A` (hits 3..7 below, 3..5 above) and turn around in
it.  So the risk set was probed where the "regular vertex cannot branch" hypothesis does
not apply, and [NCYL-219]'s refutation does not bear on it.

  H4 scores that re-reading: at a non-`R` boundary, is the reversal INSIDE the graze
  cluster (endpoint-of-the-head-on-wall type) or strictly after it (a genuine passage)?

⚠ CLUSTERS, NOT HITS ([OPS-096]).  Near a vertex the orbit runs along the corner for
3..21 consecutive bounces, so an index comparison across a cluster is meaningless -- which
is why s374's (1a) margin check (`graze_k < depth/chi`, "tightest row 1.5 bounces clear")
reports clearance at `(3,5)` while the graze and the reversal are the SAME corner event.
Every count here is over CLUSTERS, and each row reports the GAP to the nearest hit outside
the cluster so a near-tie is visible rather than silently resolved.

PRIOR ART: grepped `rulings.md` + the ledgers for 'temporal midpoint', 'cone point',
'regular vertex', 'retroreflect', 'wall-termination', 'palindromic in position' ->
[NCYL-133] (`R` is the unique reversing vertex and is REGULAR -- the input this proof
turns on), [NCYL-215] (a `J`-pair's words differ in exactly the two central letters: the
WORD-side shadow of H1, and its ⚠ says do not extend it to the ORPHAN's fold),
[NCYL-216]/[NCYL-217]/[NCYL-218]/[NCYL-219] (s374's forward proof and its two refuted
mechanisms), [NCYL-080] (Lemma B's temporal-midpoint head-on -- for the ORPHAN),
`cf_width_laws.md` §s219 (the ORPHAN's fold is a single-wall `L1` head-on and the
`R`-graze sits AWAY from it -- a different object from the `J`-fold scored here, and the
distinction is [NCYL-215]'s own warning), [OPS-096] (argmin over near-ties at a graze).
Nothing anywhere traces the LIMIT orbit at a boundary or scores position-palindromy.

USAGE
    PYTHONPATH=.:engine:archive/scripts_2026-07:archive/scripts_2026-08:probes \
      .venv/bin/python3.13 probes/s375_foldmid.py [--qmax 30]
"""
import argparse
import collections
import glob
import json
import math
import os

import cell_store as cs
from right_triangle_billiards import RightTriangleBilliard

OUT = 'data/s375_foldmid.json'


def clusters_of(hits, verts, delta):
    """Group consecutive hits lying within `delta` of the SAME vertex.

    Returns (list of (vertex, i0, i1, min_dist), gap) where `gap` is the smallest
    vertex-distance among hits NOT in any cluster -- the [OPS-096] margin.  A cluster is
    only a measurement if `gap` is far above `delta`.

    ⚠⚠ `delta` MUST scale with the launch offset `eps`, not be a fixed distance.  A
    GENUINE vertex passage sits at distance `O(eps)` from the vertex -- halve `eps` and it
    halves.  A NEAR-MISS sits at a distance fixed by the geometry and does not move.  With
    a fixed `delta = 1e-4` this probe booked 5 `R`-folds as meeting a CONE POINT at `A`;
    re-run at `eps/100` the `R` distance fell by exactly 100x (1.007e-7 -> 1.007e-9) while
    the `A` distance did not move at all (4.318e-5 -> 4.329e-5) and the hit count was
    identical, so all 5 were near-misses.  `min_dist/eps` is reported per cluster and its
    distribution printed, so the separation is auditable rather than asserted.
    """
    lab, dist = [], []
    for h in hits:
        d = {k: math.hypot(h.x - v[0], h.y - v[1]) for k, v in verts.items()}
        n = min(d, key=d.get)
        lab.append(n if d[n] < delta else None)
        dist.append(d[n])
    out, i = [], 0
    while i < len(lab):
        if lab[i] is None:
            i += 1
            continue
        j = i
        while j + 1 < len(lab) and lab[j + 1] == lab[i]:
            j += 1
        out.append((lab[i], i, j, min(dist[i:j + 1])))
        i = j + 1
    gap = min([dist[i] for i in range(len(lab)) if lab[i] is None], default=float('inf'))
    return out, gap


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--qmax', type=int, default=30)
    ap.add_argument('--maxhits', type=int, default=60000)
    a = ap.parse_args()

    tot = collections.Counter()
    by_v = collections.defaultdict(collections.Counter)
    rows, bad, ratio, nearmiss = [], [], [], []

    for p in sorted(glob.glob(os.path.join(cs.STORE_DIR, '*.json'))):
        rec = json.load(open(p))
        if rec.get('guard') or not rec.get('complete') or not rec.get('cells'):
            continue
        P, Q = rec['P'], rec['Q']
        if Q > a.qmax:
            continue
        cells = sorted(rec['cells'], key=lambda d: d['lo'])
        if any(c.get('graze_hi_V') is None for c in cells[:-1]):
            continue
        if any(c.get('chi') is None for c in cells):
            continue
        al = math.pi * P / (2.0 * Q)
        bil = RightTriangleBilliard(al)
        ct = 1.0 / math.tan(al)
        verts = {'O': (0.0, 0.0), 'R': (ct, 0.0), 'A': (ct, 1.0)}
        tot['centres'] += 1

        # ---- per-cell branch of T, from INTERIOR samples -----------------------
        # ⚠⚠ [OPS-100].  A fold is `b_lo = b_hi = 2y`, and it CANNOT be read off
        # `T(y − ε) ≈ y + ε`: on the orphan `T = id`, so `T(y − ε) = y − ε`, which is
        # within `2ε` of `y + ε` and registers as a fold for free.  The first draft of
        # this probe did exactly that and booked all 29 orphan-endpoint `A`-rows as
        # `J`-folds.  `b` must come from samples well inside each cell.
        br = []
        for c in cells:
            w = c['hi'] - c['lo']
            bb = []
            for f in (0.27, 0.63):
                s = c['lo'] + f * w
                try:
                    t = bil.trace(s, max_hits=a.maxhits)
                except Exception:                                      # noqa: BLE001
                    bb.append(None)
                    continue
                bb.append(s + t.l1_returns[-1].s
                          if (t.returned_perpendicular and t.l1_returns) else None)
            if None in bb:
                br.append(None)
            elif abs(bb[0] - bb[1]) < 1e-7:
                br.append(('refl', 0.5 * (bb[0] + bb[1])))
            else:
                br.append(('fix', None))          # T = id (the orphan)

        for i in range(len(cells) - 1):
            c, d = cells[i], cells[i + 1]
            V, y = c['graze_hi_V'], c['hi']
            eps = min(1e-7, 0.01 * min(c['hi'] - c['lo'], d['hi'] - d['lo']))
            if eps < 1e-12:
                by_v[V]['skip_thin'] += 1
                continue
            try:
                tlo = bil.trace(y - eps, max_hits=a.maxhits)
                thi = bil.trace(y + eps, max_hits=a.maxhits)
            except Exception:                                          # noqa: BLE001
                by_v[V]['skip_trace'] += 1
                continue
            if not (tlo.returned_perpendicular and thi.returned_perpendicular):
                by_v[V]['skip_noreturn'] += 1
                continue
            by_v[V]['n'] += 1

            # a fold is: BOTH flanking cells case (i), with the same b, and b = 2y
            fold = (br[i] is not None and br[i + 1] is not None
                    and br[i][0] == 'refl' and br[i + 1][0] == 'refl'
                    and abs(br[i][1] - br[i + 1][1]) < 1e-7
                    and abs(br[i][1] - 2 * y) < 1e-7)
            by_v[V]['fold'] += int(fold)
            if br[i] is None or br[i + 1] is None:
                by_v[V]['branch_unknown'] += 1

            hl = tlo.hits
            n = len(hl)
            # H1a  position-palindromy of the limit orbit about its own midpoint
            tolp = max(1e-6, 200 * eps)
            pal = n >= 2 and all(
                math.hypot(hl[k].x - hl[n - 2 - k].x, hl[k].y - hl[n - 2 - k].y) < tolp
                for k in range((n - 1) // 2))
            # H1b  the centre of that palindrome is a near-R cluster
            cl, gap = clusters_of(hl, verts, 50 * eps)
            for _v, _i0, _i1, _md in cl:
                ratio.append(_md / eps)
            mid = (n - 2) / 2.0
            centre_R = any(v == 'R' and i0 - 0.5 <= mid <= i1 + 0.5 for v, i0, i1, _ in cl)
            conepts = [v for v, _, _, _ in cl
                       if not ((v == 'R') or (v == 'A' and Q - P in (1, 2))
                               or (v == 'O' and P in (1, 2)))]
            # H4  is the reversal INSIDE the graze cluster?  the reversal is the
            # palindrome centre; the graze cluster is the FIRST cluster.
            rev_in_first = bool(cl) and cl[0][1] - 0.5 <= mid <= cl[0][2] + 0.5

            by_v[V]['pal'] += int(pal)
            by_v[V]['centre_R'] += int(centre_R)
            by_v[V]['no_conept'] += int(not conepts)
            by_v[V]['rev_in_graze_cluster'] += int(rev_in_first)
            if fold:
                by_v[V]['fold_pal'] += int(pal)
                by_v[V]['fold_centre_R'] += int(centre_R)
                by_v[V]['fold_no_conept'] += int(not conepts)
                if not (pal and centre_R):
                    bad.append(('H1', P, Q, i, V, y, pal, centre_R, n))
                if conepts:
                    bad.append(('H2', P, Q, i, V, y, conepts))
            nearmiss.append(gap / eps)
            rows.append(dict(P=P, Q=Q, i=i, V=V, y=y, fold=fold, pal=pal,
                             centre_R=centre_R, nclus=len(cl),
                             clus=[[v, i0, i1] for v, i0, i1, _ in cl],
                             conepts=conepts, gap=gap, nhits=n,
                             rev_in_graze_cluster=rev_in_first,
                             orphan_side=bool(c['chi'] == 1 or d['chi'] == 1)))

    print(f'COVERAGE: {tot["centres"]} centres, Q <= {a.qmax}; '
          f'{sum(by_v[v]["n"] for v in by_v)} interior boundaries traced')
    print(f'\n{"V":>3} {"n":>6} {"fold":>6} | {"pal":>6} {"centreR":>8} {"noCone":>7}'
          f' | {"foldPal":>8} {"foldCtrR":>9} {"foldNoCone":>11} {"revInGraze":>11}')
    for v in ('R', 'O', 'A'):
        c = by_v[v]
        print(f'{v:>3} {c["n"]:>6} {c["fold"]:>6} | {c["pal"]:>6} {c["centre_R"]:>8} '
              f'{c["no_conept"]:>7} | {c["fold_pal"]:>8} {c["fold_centre_R"]:>9} '
              f'{c["fold_no_conept"]:>11} {c["rev_in_graze_cluster"]:>11}'
              + (''.join(f'  ({k} {c[k]})' for k in
                         ('skip_thin', 'skip_trace', 'skip_noreturn') if c[k])))
    print('\nH1 WANT: every R row is a fold, palindromic, centre at R.')
    print('H2 WANT (= the residual V2): every fold meets NO cone point.')
    print('H3 CONTROL: O/A rows are NOT folds and NOT palindromic.')

    orp = [r for r in rows if r['V'] != 'R' and r['orphan_side']]
    print(f'\nH4  non-R boundaries with the ORPHAN on one side: {len(orp)}; of these, '
          f'reversal inside the graze cluster: '
          f'{sum(r["rev_in_graze_cluster"] for r in orp)}')
    reg = [r for r in rows if r['V'] == 'A' and r['Q'] - r['P'] in (1, 2)]
    print(f'    the [NCYL-219] risk set (A-boundary at r in {{1,2}}, A REGULAR): '
          f'{len(reg)} rows; orphan on a side: {sum(r["orphan_side"] for r in reg)}; '
          f'reversal inside the graze cluster: '
          f'{sum(r["rev_in_graze_cluster"] for r in reg)}')
    for r in reg:
        print(f'      ({r["P"]:3d},{r["Q"]:3d}) i={r["i"]} clusters={r["clus"]} '
              f'nhits={r["nhits"]} gap={r["gap"]:.3e} orphan_side={r["orphan_side"]}')

    # [OPS-096] done properly: the threshold is justified only if min_dist/eps is
    # BIMODAL with a wide empty band around 50.  Print both sides of it.
    ratio.sort(); nearmiss = sorted(x for x in nearmiss if x < float('inf'))
    if ratio:
        inb = [x for x in ratio if x < 50]
        print(f'\n[OPS-096] SEPARATION.  {len(ratio)} clusters accepted; '
              f'min_dist/eps in [{ratio[0]:.2f}, {ratio[-1]:.2f}] '
              f'(all < 50 by construction; {len(inb)} of them).')
    if nearmiss:
        print(f'    nearest REJECTED hit, gap/eps: min {nearmiss[0]:.3g}, '
              f'10th pct {nearmiss[len(nearmiss)//10]:.3g} -- the empty band is '
              f'[{ratio[-1] if ratio else 0:.2f}, {nearmiss[0]:.3g}].')

    if bad:
        print(f'\n⚠ VIOLATIONS ({len(bad)}):')
        for t in bad[:20]:
            print('   ', t)
    else:
        print('\nno violations of H1/H2.')

    json.dump(dict(by_vertex={k: dict(v) for k, v in by_v.items()},
                   rows=rows, bad=[list(map(str, t)) for t in bad[:400]]),
              open(OUT, 'w'), indent=1)
    print(f'wrote {OUT}')


if __name__ == '__main__':
    main()
