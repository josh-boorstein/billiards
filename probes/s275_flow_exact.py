#!/usr/bin/env python3
"""
s275_flow_exact.py — PRE-REGISTERED: RING-CERTIFY the flow identity of s274 (7).

THE BUILD ASKED FOR BY `future_directions.md` B9 q4.  s274 (7) measured
    Sum_c w_c len_c = 2*Area*Q = Q*cot(pi P/2Q)      (even P, full sweep)
with the FLOAT tracer at cell midpoints, 52/53 rows to 1e-9.  A float agreement
at 1e-11 is not an identity.  This probe redoes it in EXACT cyclotomic ring
arithmetic, so the verdict is a ring-vector zero and not a tolerance.

THE EXACT FORM (derived here, and it removes cot alpha entirely).  In the
unfolded picture the beam is the horizontal line y = s and the developed L1
copy the orbit closes on is VERTICAL, at x = Rx_c * cot(alpha) where Rx_c is the
x-coefficient RingElem of the state find_orphan_ring's closure test stops in
(engine/orphan_region.py, leaf_detail=True, reason == 'closed').  The launch is
at x = cot(alpha), so the developed horizontal displacement is
        T_c = cot(alpha) * (1 - Rx_c),
and dividing the identity through by cot(alpha) leaves a statement with NO
transcendental factor at all:
        Sum_c w_c * (1 - Rx_c) / chi_c  ==  Q         exactly in Z[zeta_2Q].
Both sides are ring elements; the test is `.is_zero()` on the difference.

chi_c: MEASURED, not assumed (H_A below).  T_c above is the FULL PERIOD of the
closed geodesic, not the first-return length to Sigma_1 = {L1, perpendicular} --
the closure test demands the developed copy be a translate / x-mirror of the
original, which is a period closure.  A PAIRED cell's geodesic meets Sigma_1
TWICE per period (at s and at J(s) != s, s272 (1)'s double crossing), so its
first-return length is T_c/2; the ORPHAN's meets it once, so its return length
is T_c.  Hence chi_c = 2 - [c is the orphan], read from the leaf's EXACT headon
flag (_exact_vertical, a ring-zero test).  At even P there is no orphan and
chi == 2 identically.

HYPOTHESES
  H_A (the geometry).  cot(alpha)*(1 - Rx_c) / chi_c equals the float tracer's
      first-return length at the cell midpoint.  Tolerance 1e-9 relative.  This
      is the bridge between the exact quantity and the object s274 (7) measured;
      without it the ring test below could be certifying a different number.
  H_B (the identity, THE DISCHARGE CONDITION).  Sum_c w_c (1-Rx_c)/chi_c - Q is
      ring-ZERO on every guarded even-P row, EXCEPT the deficit row 8/15, which
      must come out non-zero (it has no orphan and g_len = 0.454, s274 (7)).
  H_C (P = 1).  Sum_c w_c (1-Rx_c)/chi_c == 2 exactly, i.e. g_len = 2/Q as an
      exact rational.  s274 (7) had this at float 0.0e+00 on 8 rows.
  H_D (odd P >= 3).  Reported, not asserted: the exact swept fraction
      g_len = (that sum)/Q.

GUARDS (a truncated or merged enumeration would fake a deficit).  A row is
CERTIFIED iff BOTH:
  (i)  every leaf's reason == 'closed' -- no 'depth', 'width' or 'degenerate'
       leaf (the s273 (8) certificate hole, now surfaced per leaf by the engine);
  (ii) Sum_c w_c == 1 EXACTLY as a ring element.  This is strictly stronger than
       s274's float 1e-9 version: find_orphan_ring drops leaves narrower than
       1e-13, and any such drop breaks the exact sum.  (i)+(ii) is also a
       completeness certificate that does NOT beg the question the way
       `n == Q-1` does -- s272 (4)'s harness blind spot.
A row failing a guard is REPORTED AS SKIPPED, never as a deviation.  The n-even
parity detector of s274 (7)(i) is computed and reported as a cross-check but is
subsumed by (i)+(ii).

PRIOR ART: grepped 'flow identity', '2Q*Area', 'Q cot', 'ring-certif',
'len_c', 'g_len', 'Rx_c', 'period', 'first-return', 'leaf_detail',
'completeness certificate' across all .md and TOOLS.md ->
  - `cf_width_laws.md` s274 (7) OWNS the float version and states the open build;
    `future_directions.md` B9 q4 states the discharge condition satisfied here.
  - `rulings.md` [NCYL] 230-232 own "length is constant on a cell", the flow
    identity as a FLOAT regularity, and P=1's g_len = 2/Q.
  - [NCYL] 246 owns the certificate hole (`unclosed == 0` alone is not enough);
    [OPS] (s272) owns "a completeness certificate must not be an instance of the
    hypothesis under test".  Guards (i)+(ii) are built to satisfy both.
  - [NCYL] 250 owns the double-crossing mechanism -- chi_c here IS that factor,
    measured per cell rather than assumed from the count.
  - No prior exact/ring computation of a PATH LENGTH anywhere in the repo: the
    ring is used for widths (orphan_width_ring) and words, never for lengths.

Run: PYTHONPATH=.:engine:archive/scripts_2026-07:archive/scripts_2026-08:probes \
     .venv/bin/python3.13 probes/s275_flow_exact.py
"""
import json
import math
import time

from cos_poly import CosPoly
from orphan_region import find_orphan_ring
from right_triangle_billiards import RightTriangleBilliard
from ring_cyc import RingContext

MAX_DEPTH = 200000
QMAX = 17          # even-P sweep range (Q odd); 20 rows, well past the 10 required


def _leaf_endpoints(d, lo_end, hi_end, lo_val, hi_val):
    """Exact (lo, hi) RingElems for one leaf -> (lo, hi, err).

    A leaf's `lo_poly`/`hi_poly` is None exactly at an OUTER endpoint of the
    descent interval, whose exact value the caller supplies (`lo_end`/`hi_end`,
    with `lo_val`/`hi_val` their floats).  For a full run that is 0 and 1; for a
    resumed WINDOW it is the pending interval's own stored endpoints.
    ⚠ Check order (hi, then lo) is `cells_exact`'s and is preserved so the two
    paths cannot return different guard strings on the same input.
    """
    hi = d['hi_poly']
    if hi is None:
        if abs(d['hi'] - hi_val) > 1e-9:
            return None, None, 'unmatched upper boundary'
        hi = hi_end
    lo = d['lo_poly']
    if lo is None:
        if abs(d['lo'] - lo_val) > 1e-9:
            return None, None, 'unmatched lower boundary'
        lo = lo_end
    return lo, hi, None


def _row_from_leaf(d, lo, hi, one):
    """The per-cell record. `T` is T_c / cot(alpha), an exact ring element."""
    return {'lo': d['lo'], 'hi': d['hi'], 'w': hi - lo,
            'T': one - d['state'].Rx_c,
            'chi': 1 if d['headon'] else 2,
            'vertical': (d['state'].Rx_c - d['state'].Ax_c).is_zero(),
            'depth': d['depth']}


def cells_exact(P, Q, max_depth=MAX_DEPTH):
    """Exact per-cell record. Returns (ctx, rows, guard) where guard is None if
    the row is certified, else the reason it was skipped."""
    ctx = RingContext(P, Q)
    one = ctx.from_cospoly(CosPoly([1], 0))
    zero = ctx.zero()
    L = find_orphan_ring(P, Q, 0.0, 1.0, max_depth=max_depth,
                         return_leaves=True, leaf_detail=True)
    bad = sorted({d['reason'] for d in L if d['reason'] != 'closed'})
    if bad:
        return ctx, None, f"uncertified leaves: {','.join(bad)}"
    rows, W = [], ctx.zero()
    for d in L:
        lo, hi, err = _leaf_endpoints(d, zero, one, 0.0, 1.0)
        if err:
            return ctx, None, err
        W = W + (hi - lo)
        rows.append(_row_from_leaf(d, lo, hi, one))
    W.reduce()
    if not (W - one).reduce().is_zero():
        return ctx, None, f'Sum w != 1 exactly (float {W.to_float():.3e})'
    if not all(r['vertical'] for r in rows):
        return ctx, None, 'closure state has a non-vertical L1 copy'
    return ctx, rows, None


def cells_partial(P, Q, max_depth=MAX_DEPTH, window=None, window_polys=None):
    """PARTIAL certification: keep the CLOSED cells, hand back the rest as
    resumable intervals.  Returns `(ctx, rows, pending, guard)`.

    WHY (s378, item (g)).  `cells_exact` is all-or-nothing: one leaf that does
    not close discards every certified cell in the centre, and `ring_cache`
    then stores `cells: []`.  Measured on five guarded centres, 45 of 51 leaves
    were individually certified and thrown away that way -- [OPS-070](b),
    "discarding state also discards the audit".  Each `closed` leaf carries its
    OWN certificate; what a guard destroys is only the COMPLETENESS certificate.

    `window=(lo, hi)` restricts the descent to that interval via
    `find_orphan_ring`'s existing `s_lo`/`s_hi` -- the engine is already
    region-restricted (`overlaps()` prunes the stack), which is what makes a
    resume cheap: the deep cost is paid on the sliver, not the transversal.
    s278 named exactly this ("a window-pruned descent, not a bigger max_steps",
    `cf_width_laws.md`) and [SEAM-001] records one done by hand to 400000.
    `window_polys=(A, B)` supplies the window's EXACT endpoints; without them a
    resumed window could not close its own invariant.

    INVARIANT, ASSERTED: `sum(rows.w) + sum(pending.w) == hi_end - lo_end` in the
    ring.  This is the analogue of `cell_store.verify()`'s partial invariant for
    glen, and it is also what checks the resume assumption -- that a windowed
    re-descent emits only leaves INSIDE the window.

    ⚠ `guard` keeps `cells_exact`'s string format, but the reason set separates
    `width` from `depth` rather than lumping them: [WFLOOR-065] is the precedent
    that a silent `min_width` prune once hid cells and a session mis-read it.
    """
    ctx = RingContext(P, Q)
    one = ctx.from_cospoly(CosPoly([1], 0))
    zero = ctx.zero()
    if window is None:
        s_lo, s_hi, lo_end, hi_end = 0.0, 1.0, zero, one
    else:
        s_lo, s_hi = window
        A, B = window_polys if window_polys else (None, None)
        lo_end = zero if A is None else A
        hi_end = one if B is None else B
    L = find_orphan_ring(P, Q, s_lo, s_hi, max_depth=max_depth,
                         return_leaves=True, leaf_detail=True)
    rows, pending, W = [], [], ctx.zero()
    for d in L:
        lo, hi, err = _leaf_endpoints(d, lo_end, hi_end, s_lo, s_hi)
        if err:
            return ctx, None, None, err
        w = hi - lo
        W = W + w
        if d['reason'] == 'closed':
            rows.append(_row_from_leaf(d, lo, hi, one))
        else:
            pending.append({'lo': d['lo'], 'hi': d['hi'], 'w': w,
                            'lo_poly': d['lo_poly'], 'hi_poly': d['hi_poly'],
                            'reason': d['reason'], 'depth': d['depth']})
    W.reduce()
    span = (hi_end - lo_end).reduce()
    if not (W - span).reduce().is_zero():
        return ctx, None, None, f'Sum w != span exactly (float {W.to_float():.3e})'
    if not all(r['vertical'] for r in rows):
        return ctx, None, None, 'closure state has a non-vertical L1 copy'
    if pending:
        bad = sorted({p['reason'] for p in pending})
        return ctx, rows, pending, f"uncertified leaves: {','.join(bad)}"
    return ctx, rows, [], None


def swept_sum(ctx, rows):
    """Sum_c w_c (1 - Rx_c)/chi_c as an exact RingElem."""
    S = ctx.zero()
    for r in rows:
        t = r['w'] * r['T']
        if r['chi'] == 2:
            t = t.div2()
        S = S + t
    return S.reduce()


def main():
    out = {'even': [], 'p1': [], 'odd': [], 'HA': []}

    # ---- H_A: the exact length IS the tracer's first-return length * chi ----
    print('H_A — cot(a)(1-Rx_c)/chi == float first-return length?')
    worst, ncell = 0.0, 0
    for (P, Q) in [(2, 7), (4, 15), (2, 3), (3, 11), (1, 5), (5, 13)]:
        ctx, rows, guard = cells_exact(P, Q)
        if guard:
            print(f'    {P}/{Q}: SKIPPED ({guard})')
            continue
        cot = 1.0 / math.tan(ctx.alpha)
        B = RightTriangleBilliard(ctx.alpha)
        for r in rows:
            s = 0.5 * (r['lo'] + r['hi'])
            tr = B.trace(s, max_hits=400000)
            if not tr.returned_perpendicular:
                continue
            pts = [B.start_point(s)] + [(h.x, h.y) for h in tr.hits]
            flt = sum(math.hypot(pts[i + 1][0] - pts[i][0], pts[i + 1][1] - pts[i][1])
                      for i in range(len(pts) - 1))
            ex = cot * r['T'].to_float() / r['chi']
            worst = max(worst, abs(ex - flt) / flt)
            ncell += 1
        out['HA'].append({'P': P, 'Q': Q, 'cells': len(rows),
                          'chi1': sum(1 for r in rows if r['chi'] == 1)})
    print(f'    {ncell} cells over {len(out["HA"])} rows, worst relative '
          f'error {worst:.1e} -> {"PASS" if worst < 1e-9 else "FAIL"}\n')
    out['HA_worst'] = worst

    # ---- H_B: the identity, exactly, on even P -----------------------------
    print(f'H_B — Sum_c w_c(1-Rx_c)/chi_c == Q EXACTLY in Z[zeta_2Q]  (even P)')
    even = [(P, Q) for Q in range(3, 26, 2) for P in range(2, Q, 2)
            if math.gcd(P, Q) == 1]
    ok = dev = skip = 0
    for (P, Q) in even:
        t0 = time.time()
        try:
            ctx, rows, guard = cells_exact(P, Q)
        except Exception as e:                                   # noqa: BLE001
            print(f'    skipped  {P}/{Q}: {type(e).__name__}: {e}')
            out['even'].append({'P': P, 'Q': Q, 'skip': str(e)})
            skip += 1
            continue
        if guard:
            print(f'    skipped  {P}/{Q}: {guard}')
            out['even'].append({'P': P, 'Q': Q, 'skip': guard})
            skip += 1
            continue
        S = swept_sum(ctx, rows)
        exact = (S - ctx.from_cospoly(CosPoly([Q], 0))).reduce().is_zero()
        rec = {'P': P, 'Q': Q, 'n': len(rows), 'exact': exact,
               'sum_float': S.to_float(), 'g_len': S.to_float() / Q,
               'n_even': len(rows) % 2 == 0, 'secs': round(time.time() - t0, 1),
               'max_depth_used': max(r['depth'] for r in rows)}
        out['even'].append(rec)
        tag = 'EXACT' if exact else 'NOT EXACT'
        print(f'    {tag:10s} {P}/{Q}: n={rec["n"]:3d} sum={rec["sum_float"]:.9f} '
              f'Q={Q} g_len={rec["g_len"]:.6f} ({rec["secs"]}s)', flush=True)
        if exact:
            ok += 1
        else:
            dev += 1
    print(f'    {ok} rows EXACT (ring-zero), {dev} not exact, {skip} skipped by guards')
    d815 = next((r for r in out['even'] if (r['P'], r['Q']) == (8, 15)), None)
    if d815 is not None:
        verdict = ('SKIPPED: ' + d815['skip']) if 'skip' in d815 else (
            'WRONG — came out exact' if d815['exact'] else
            f'as required: NOT exact, g_len={d815["g_len"]:.5f}')
        print(f'    8/15 (the deficit row): {verdict}')

    # ---- H_C: P = 1 --------------------------------------------------------
    print('\nH_C — P=1: Sum == 2 exactly (g_len = 2/Q)')
    for Q in (5, 7, 9, 11, 13, 17, 21, 25):
        ctx, rows, guard = cells_exact(1, Q)
        if guard:
            print(f'    1/{Q}: SKIPPED ({guard})')
            continue
        S = swept_sum(ctx, rows)
        two = (S - ctx.from_cospoly(CosPoly([2], 0))).reduce().is_zero()
        out['p1'].append({'Q': Q, 'n': len(rows), 'sum_is_2': two,
                          'sum_float': S.to_float()})
    print(f'    {sum(r["sum_is_2"] for r in out["p1"])}/{len(out["p1"])} rows have '
          f'the sum EXACTLY 2 -> g_len = 2/Q exactly')

    # ---- H_D: odd P >= 3, reported ----------------------------------------
    print('\nH_D — odd P >= 3 (reported, not asserted): exact swept fraction')
    for (P, Q) in [(3, 7), (3, 8), (3, 10), (3, 11), (5, 12), (3, 13), (5, 13),
                   (7, 16), (9, 20)]:
        ctx, rows, guard = cells_exact(P, Q)
        if guard:
            print(f'    {P}/{Q}: SKIPPED ({guard})')
            continue
        S = swept_sum(ctx, rows)
        g = S.to_float() / Q
        norphan = sum(1 for r in rows if r['chi'] == 1)
        out['odd'].append({'P': P, 'Q': Q, 'n': len(rows), 'orphans': norphan,
                           'sum_float': S.to_float(), 'g_len': g})
        print(f'    {P}/{Q}: n={len(rows)} orphans={norphan} g_len={g:.6f}')

    def _clean(o):
        return o
    json.dump(out, open('data/s275_flow_exact.json', 'w'), indent=1, default=str)
    print('\nwrote data/s275_flow_exact.json')


if __name__ == '__main__':
    main()
