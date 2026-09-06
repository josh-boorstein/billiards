"""s374 PRE-REGISTERED -- is the first-return involution `T` CONTINUOUS at a cell
boundary, and is continuity equivalent to the boundary being an `R`-graze?

WHY THIS IS THE DECIDING MEASUREMENT FOR ITEM (f).  Working the geometry reduces
[NCYL-216] to two statements, of which the forward half is provable outright:

  FORWARD (proved, see the session write-up).  `R` is the unique vertex whose two side
  reflections compose to `−id`, and it is REGULAR ([NCYL-133]).  A boundary graze happens
  strictly before the perpendicular return, so the ray launched at an `R`-graze height `y`
  retraces and comes back to `y`: `T(y) = y`, and `T` is continuous there because the
  development does not branch at a regular point.  With `T(s) = b − s` on each cell
  (Thm A / Lemma B'), continuity forces `b = 2y` on BOTH flanking cells, so `T` reflects
  about `y` -- the neighbour IS the partner, of equal width, and neither is the orphan.

  RESIDUAL (R★), the whole remaining content.  At a NON-`R` boundary, is `T`
  discontinuous?  If yes, then `b` changes at every non-`R` boundary, each `b`-class is a
  run of cells joined by `R`-grazes, two consecutive `R`-grazes are impossible
  (`b = 2y₁ = 2y₂`), so every run has ≤ 2 cells -- which is exactly [NCYL-216]'s converse
  AND [NCYL-186]'s unproved "positionally adjacent" detail.

So this probe measures `b = s + T(s)` on both sides of every interior boundary and asks
whether `b` jumps.  `T` is read from the float tracer's PERPENDICULAR return
(`right_triangle_billiards.trace`, `stop_at_perpendicular_return=True`) -- note the first
`L1` hit is NOT generally the return; only a perpendicular one is ([NCYL-046], and
[NCYL-026]'s "a perpendicular `L1` hit is a VERTICAL developed copy").

  H4  Lemma 2's DICHOTOMY, per cell: either `T(s) = b − s` with `b` constant (case (i),
      slope −1, the cell maps to its partner) or `T(s) = s` identically (case (ii), the
      whole branch RETRACES).  ⚠ INSTRUMENT CHECK, not a finding.  ⚠⚠ The first draft of
      this probe scored only "`b = s + T(s)` is constant" and took 83 hits, ALL of them
      case (ii): on the orphan `T` is the identity, so `b = 2s` varies BY CONSTRUCTION.
      Lemma B says the orphan branch retraces -- the whole BRANCH, not one point -- and
      `P = 1` is entirely of this kind (one cell, `T = id`).  Scoring the slope-−1 form on
      a slope-+1 cell is a wrong-object error, not a tracer error.
  H4b `T = id` on a cell  <=>  `chi == 1`.  A genuine CROSS-ENGINE check: `chi` is
      RING-measured, the retrace is FLOAT-traced, and neither knows about the other.
  H5  `T` continuous at a boundary  <=>  the boundary is an `R`-graze.       [= (R★)]
  H6  at an `R` boundary the common limit is `T(y) = y` (a genuine FOLD, i.e. the ray
      launched at the boundary comes back to its own launch height).  H6 is the forward
      proof's own prediction, so a failure refutes the proof.

⚠ THE DISCRIMINATING ROWS ARE THE NON-`R` BOUNDARIES, and specifically the `A`-grazes at
`r = Q − P ∈ {1,2}`, where `A` is a REGULAR point (cone angle `2π·m/gcd(m,2Q)` with
interior angle `mπ/2Q`; `gcd(r,2Q) = gcd(r,2)`).  Naive sheet reasoning says a regular
vertex cannot branch the development, hence cannot move the return wall, hence H5 should
FAIL exactly there.  `probes/s374_graze_structure.py` says it does not fail -- the `r = 2`
family carries 6 such `A`-boundaries and still tiles perfectly -- so either the naive
reasoning is wrong or the effect hides somewhere else.  ⇒ These rows are reported
SEPARATELY; the pooled score would drown 6 rows in ~1500.

PRIOR ART: grepped `rulings.md` for 'continuous', 'return map', 'first return', 'T(s)',
'involution', 'NCYL-046', 'NCYL-026', 'NCYL-186', and `TOOLS.md` for the tracer ->
[NCYL-046] (Thm A double crossing; the perpendicular return is the object, and the
regular-vertex-graze route was left "POSSIBLE but NOT ESTABLISHED"), [NCYL-026] (a
perpendicular `L1` hit is a VERTICAL developed copy at an `x` independent of `s`, so `b`
is a cell invariant -- this is H4's source), [NCYL-186] (adjacency measured `198/198`,
unproved).  Nothing scores `T`-continuity per boundary, and nothing anywhere in the repo
classifies boundaries by whether `b` jumps.  TOOLS.md row: `right_triangle_billiards`
(float tracer, cheap localisation), trap T9 noted (not used here).

USAGE
    PYTHONPATH=.:engine:archive/scripts_2026-07:archive/scripts_2026-08:probes \
      .venv/bin/python3.13 probes/s374_tcont.py [--qmax 30]
"""
import argparse
import collections
import glob
import json
import math
import os

import cell_store as cs
from right_triangle_billiards import RightTriangleBilliard

OUT = 'data/s374_tcont.json'


def b_at(bil, s):
    """`b = s + T(s)` for the perpendicular launch at height `s`, or None."""
    try:
        tr = bil.trace(s, max_hits=200000)
    except Exception:                                          # noqa: BLE001
        return None
    if not tr.returned_perpendicular or not tr.l1_returns:
        return None
    return s + tr.l1_returns[-1].s


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--qmax', type=int, default=30)
    a = ap.parse_args()

    tot = collections.Counter()
    by_v = collections.defaultdict(collections.Counter)
    risk = []          # the regular-A rows, reported separately
    bad = []
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
        bil = RightTriangleBilliard(math.pi * P / (2.0 * Q))

        # --- classify each cell by Lemma 2's dichotomy (H4 = instrument check) --
        # branch[i] is None, or ('fix', None) for T = id, or ('refl', b) for T = b - s.
        branch = []
        for c in cells:
            lo, hi = c['lo'], c['hi']
            w = hi - lo
            s1, s2 = lo + 0.27 * w, lo + 0.63 * w
            b1, b2 = b_at(bil, s1), b_at(bil, s2)
            if b1 is None or b2 is None:
                branch.append(None)
                tot['cell_untraced'] += 1
                continue
            tot['cell_traced'] += 1
            if abs(b1 - b2) < 1e-7:                       # T(s) = b - s
                branch.append(('refl', 0.5 * (b1 + b2)))
            elif abs((b1 - 2 * s1)) < 1e-7 and abs(b2 - 2 * s2) < 1e-7:
                branch.append(('fix', None))              # T(s) = s
            else:
                branch.append(None)
                bad.append(('H4', P, Q, lo, hi, b1, b2))
                continue
            tot['H4_ok'] += 1
            is_fix = branch[-1][0] == 'fix'
            if is_fix == (c['chi'] == 1):
                tot['H4b_ok'] += 1
            else:
                bad.append(('H4b', P, Q, lo, hi, branch[-1][0], c['chi']))
        bs = branch

        # --- per interior boundary --------------------------------------------
        for i in range(len(cells) - 1):
            v = cells[i]['graze_hi_V']
            y = cells[i]['hi']
            if bs[i] is None or bs[i + 1] is None:
                by_v[v]['untraced'] += 1
                continue

            def limit(br, at):
                """One-sided limit of T at the boundary, from a cell's own branch."""
                return at if br[0] == 'fix' else br[1] - at

            tl, tr = limit(bs[i], y), limit(bs[i + 1], y)
            cont = abs(tl - tr) < 1e-7
            fold = cont and abs(tl - y) < 1e-7
            by_v[v]['n'] += 1
            by_v[v]['cont'] += int(cont)
            by_v[v]['fold'] += int(fold)
            if (v == 'R') != cont:
                bad.append(('H5', P, Q, i, v, y, tl, tr))
            if v == 'R' and not fold:
                bad.append(('H6', P, Q, i, v, y, tl, tr))
            # --- MECHANISM TEST (s374).  Conjecture: `T` jumps at `A` because a
            # VERTICAL developed `L1` wall TERMINATES at `A` -- `L1` is incident to
            # `A`, and the interior angle there is not `π/2`, so the reflected copy
            # is not collinear and the wall simply ends.  The ray on the side that
            # still meets that wall must then return AT its endpoint, i.e. at the
            # `A`-end of `L1`, which is `s = 1`.  `L1` is NOT incident to `O`, so no
            # wall ends there and this prediction should FAIL at `O` -- that
            # asymmetry is the test ([OPS-041]: it must be able to come out
            # otherwise).
            ends = min(abs(t - e) for t in (tl, tr) for e in (0.0, 1.0))
            by_v[v]['endpoint_return'] += int(ends < 1e-9)
            r = Q - P
            if v != 'R' and r in (1, 2):
                risk.append(dict(P=P, Q=Q, r=r, V=v, y=y, T_lo=tl, T_hi=tr,
                                 cont=cont, fold=fold, jump=tr - tl))
        tot['centres'] += 1

    print(f'COVERAGE: {tot["centres"]} centres with Q <= {a.qmax}; cells traced '
          f'{tot["cell_traced"]} (untraced {tot["cell_untraced"]})')
    print(f'H4  Lemma-2 dichotomy holds per cell (INSTRUMENT CHECK): '
          f'{tot["H4_ok"]}/{tot["cell_traced"]}')
    print(f'H4b T = id  <=>  chi == 1 (CROSS-ENGINE: float retrace vs ring chi): '
          f'{tot["H4b_ok"]}/{tot["H4_ok"]}')
    print(f'\nH5/H6 per boundary vertex:')
    print(f'{"V":>3} {"n":>6} {"T continuous":>13} {"fold T(y)=y":>13} '
          f'{"a side returns at s=0/1":>24}')
    for v in ('R', 'O', 'A'):
        c = by_v[v]
        print(f'{v:>3} {c["n"]:>6} {c["cont"]:>13} {c["fold"]:>13} '
              f'{c["endpoint_return"]:>24}'
              + (f'   (untraced {c["untraced"]})' if c['untraced'] else ''))
    print('\nWANT: R -> continuous AND fold on every row;  O/A -> 0 continuous.')

    print(f'\n⇒ THE DISCRIMINATING SUBSET -- non-R boundaries at r in {{1,2}} '
          f'(A is a REGULAR point there): {len(risk)} rows')
    for d in risk[:25]:
        print(f'   ({d["P"]:3d},{d["Q"]:3d}) r={d["r"]} V={d["V"]} y={d["y"]:.6f} '
              f'T(y-)={d["T_lo"]:.6f} T(y+)={d["T_hi"]:.6f} '
              f'jump={d["jump"]:+.3e} cont={d["cont"]}')

    if bad:
        print(f'\n⚠ VIOLATIONS ({len(bad)}):')
        for t in bad[:20]:
            print('   ', t)
    else:
        print('\nno violations of H4/H5/H6.')

    json.dump(dict(totals=dict(tot),
                   by_vertex={k: dict(v) for k, v in by_v.items()},
                   risk=risk, bad=[list(map(str, t)) for t in bad[:400]]),
              open(OUT, 'w'), indent=1)
    print(f'wrote {OUT}')


if __name__ == '__main__':
    main()
