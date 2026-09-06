#!/usr/bin/env python3
"""
ring_cache.py -- a PERSISTENT on-disk cache for `s275_flow_exact.cells_exact`.

WHY THIS EXISTS (s344, user-raised mid-session: "you have to rerun everything?
nothing is saved?").  `cells_exact(P, Q)` is the expensive step in every exact ring
question -- minutes per row -- and it is a PURE FUNCTION of `(P, Q, max_depth)`.  Session 344
ran eight separate exact-ring probes and recomputed the same rows in each: `(3,1) Q=13`
five times, `(7,3) Q=17/24/31` four times each, ~90 calls covering ~40 distinct rows.  That
is [OPS-052]'s lesson -- the same session had already applied it on the FLOAT side (the
`s344_position` store) and simply did not think to apply it to the ring.

So: call `cells_cached` instead of `cells_exact` and the second question about a row is free.

    from ring_cache import cells_cached
    ctx, rows, guard = cells_cached(P, Q, max_depth=400000)

Same return contract as `cells_exact`: `(ctx, rows, guard)`, `guard is None` iff certified.
`rows` entries carry `lo`, `hi`, `w`, `T` (RingElem), `chi`, `vertical`, `depth` exactly as
`cells_exact` produces them.

WHAT IS STORED.  One JSON per row at `data/ring_cells/<P>_<Q>.json`: each cell's `lo`/`hi`
floats plus the INTEGER vectors and 2-power shifts of `w` and `T`, the `chi` flag, and the
`max_depth` the row was computed at.  `RingContext(P, Q)` is rebuilt on load (cheap next to
the enumeration) and the vectors are wrapped back into `RingElem`.

⚠ TRAPS.
 (i)  **A guard-fail is cached WITH its `max_depth`.** A later call at a HIGHER `max_depth`
      recomputes; a later call at the same or lower depth returns the cached guard.  A
      CERTIFIED row is returned for any requested depth -- certification does not weaken.
 (ii) The cache is keyed by `(P, Q)` ONLY, not by the enumerator version.  If
      `s275_flow_exact` / `orphan_region` change semantics (as `s340`/`s341` changed the
      head-on predicate), **delete `data/ring_cells/` or bump `CACHE_TAG`** -- a stale cache
      is worse than no cache.  `CACHE_TAG` is written into every file and a mismatch
      silently recomputes, so bumping it is the safe move.
 (iii) It caches EXACT ring data only.  It is not a float store and it does not memoise
      `find_orphan_ring` itself.
"""
import json
import os

import numpy as np

from ring_cyc import RingContext, RingElem

CACHE_DIR = os.environ.get('RING_CACHE_DIR', 'data/ring_cells')
CACHE_TAG = 's341-headon-fix'          # bump when the enumerator's semantics change


def _path(P, Q):
    return os.path.join(CACHE_DIR, f'{P}_{Q}.json')


def _dump_elem(e):
    return {'vec': [int(x) for x in e.vec], 'shift': int(e.shift)}


def _load_elem(d, ctx):
    return RingElem(np.asarray(d['vec'], dtype=np.int64), d['shift'], ctx)


def _read_rec(P, Q):
    """The raw on-disk record, or None if absent/unreadable/stale-tagged."""
    p = _path(P, Q)
    if not os.path.exists(p):
        return None
    try:
        with open(p) as f:
            rec = json.load(f)
    except Exception:                                        # noqa: BLE001
        return None
    return rec if rec.get('tag') == CACHE_TAG else None


def _load_rows(cells, ctx):
    return [{'lo': c['lo'], 'hi': c['hi'],
             'w': _load_elem(c['w'], ctx), 'T': _load_elem(c['T'], ctx),
             'chi': c['chi'], 'vertical': c['vertical'], 'depth': c['depth']}
            for c in cells]


def _load_pending(pend, ctx):
    return [{'lo': p['lo'], 'hi': p['hi'], 'w': _load_elem(p['w'], ctx),
             'lo_poly': None if p['lo_poly'] is None else _load_elem(p['lo_poly'], ctx),
             'hi_poly': None if p['hi_poly'] is None else _load_elem(p['hi_poly'], ctx),
             'reason': p['reason'], 'depth': p['depth']} for p in pend]


def load(P, Q, max_depth, partial=False):
    """Return (ctx, rows, guard) from disk, or None if absent/stale/too shallow.

    ⚠ The DEFAULT return is unchanged and must stay that way -- `(ctx, rows, guard)`
    with `rows is None` whenever a guard is set.  46 call sites across ~20 probes
    consume this contract.

    `partial=True` (s378) returns `(ctx, rows, pending, guard)` instead, so a
    guarded row can hand back the cells it DID certify.  A LEGACY guarded record
    (written before s378: `cells: []`, no `pending` key) carries no information
    about what was resolved, so it is reported ABSENT rather than as "nothing
    pending" -- reading it the other way would silently claim a centre has no
    recoverable cells.
    """
    rec = _read_rec(P, Q)
    if rec is None:
        return None
    # a guard-fail is only reusable if it was computed at least as deep as we now want
    if rec['guard'] is not None and rec['max_depth'] < max_depth:
        return None
    ctx = RingContext(P, Q)
    if not partial:
        if rec['guard'] is not None:
            return ctx, None, rec['guard']
        return ctx, _load_rows(rec['cells'], ctx), None
    if rec['guard'] is not None and rec.get('pending') is None:
        return None                                  # legacy guard: no anatomy
    return (ctx, _load_rows(rec['cells'], ctx),
            _load_pending(rec.get('pending') or [], ctx), rec['guard'])


def store(P, Q, max_depth, rows, guard, pending=None):
    """⚠ s378: cells are written EVEN WHEN GUARDED.

    Until s378 this wrote `'cells': []` on any guard, which is the single line
    that destroyed the certified cells of 38 centres -- [OPS-070](b).  A guarded
    record now carries what it resolved plus the unresolved intervals (with their
    EXACT ring endpoints, without which a resumed window cannot close its own
    invariant).  No `CACHE_TAG` bump: per-cell semantics are unchanged and no
    existing record changes meaning -- this is a schema addition.
    """
    os.makedirs(CACHE_DIR, exist_ok=True)
    rec = {'tag': CACHE_TAG, 'P': P, 'Q': Q, 'max_depth': max_depth, 'guard': guard,
           'cells': [] if rows is None else
           [{'lo': float(c['lo']), 'hi': float(c['hi']),
             'w': _dump_elem(c['w']), 'T': _dump_elem(c['T']),
             'chi': int(c['chi']), 'vertical': bool(c['vertical']),
             'depth': int(c['depth'])} for c in rows]}
    if pending is not None:
        rec['pending'] = [
            {'lo': float(p['lo']), 'hi': float(p['hi']), 'w': _dump_elem(p['w']),
             'lo_poly': None if p['lo_poly'] is None else _dump_elem(p['lo_poly']),
             'hi_poly': None if p['hi_poly'] is None else _dump_elem(p['hi_poly']),
             'reason': p['reason'], 'depth': int(p['depth'])} for p in pending]
    tmp = _path(P, Q) + '.tmp'
    with open(tmp, 'w') as f:
        json.dump(rec, f)
    os.replace(tmp, _path(P, Q))       # atomic: a killed job never leaves a half file


def cells_cached(P, Q, max_depth=200000, verbose=False):
    """`s275_flow_exact.cells_exact` with a persistent on-disk cache.  Same contract."""
    hit = load(P, Q, max_depth)
    if hit is not None:
        if verbose:
            print(f'    [cache hit {P}/{Q}]', flush=True)
        return hit
    from s275_flow_exact import cells_exact          # imported late: heavy module
    ctx, rows, guard = cells_exact(P, Q, max_depth=max_depth)
    store(P, Q, max_depth, rows, guard)
    if verbose:
        print(f'    [computed + cached {P}/{Q}]', flush=True)
    return ctx, rows, guard


def cells_resumed(P, Q, max_depth=200000, verbose=False):
    """`cells_partial` with a cache that RESUMES instead of restarting (s378).

    Returns `(ctx, rows, pending, guard)`.  The point of the function: when a
    cached PARTIAL exists at a shallower depth, only its unresolved windows are
    re-run at the deeper cap -- `find_orphan_ring` is region-restricted, so the
    deep cost is paid on the sliver rather than on the whole transversal.  That
    is why "raise max_depth" cost 9x for 3x depth at (5,44) ([OPS-099]) and why
    this need not.

    ⚠ It does NOT make a hard centre certifiable -- narrow is not cheap
    ([OPS-070](c)).  It makes the attempt payable and the partial result keepable.
    """
    from cos_poly import CosPoly
    from s275_flow_exact import cells_partial            # imported late: heavy
    ctx = RingContext(P, Q)
    rec = _read_rec(P, Q)

    # certified already: certification does not weaken with depth
    if rec is not None and rec['guard'] is None and rec.get('cells'):
        return ctx, _load_rows(rec['cells'], ctx), [], None

    resumable = (rec is not None and rec.get('pending') is not None
                 and rec.get('cells') is not None
                 and rec.get('max_depth', 0) < max_depth)
    if resumable:
        rows = _load_rows(rec['cells'], ctx)
        new_pending = []
        for p in _load_pending(rec['pending'], ctx):
            _c, r, p2, g2 = cells_partial(P, Q, max_depth=max_depth,
                                          window=(p['lo'], p['hi']),
                                          window_polys=(p['lo_poly'], p['hi_poly']))
            if r is None:
                return ctx, None, None, g2        # a real error, not a cap stop
            rows.extend(r)
            new_pending.extend(p2)
        rows.sort(key=lambda d: d['lo'])
        if verbose:
            print(f'    [resumed {P}/{Q}: {len(rec["pending"])} windows -> '
                  f'{len(new_pending)} still pending]', flush=True)
    else:
        _c, rows, new_pending, _g = cells_partial(P, Q, max_depth=max_depth)
        if rows is None:
            return ctx, None, None, _g

    # THE MERGE INVARIANT, asserted: a resumed record must still tile [0,1].
    # This is what catches a bad merge, and it is the same check `cells_partial`
    # applies per window ([OPS-043]: assert the diagnosis in code).
    one = ctx.from_cospoly(CosPoly([1], 0))
    W = ctx.zero()
    for c in rows:
        W = W + c['w']
    for p in new_pending:
        W = W + p['w']
    if not (W.reduce() - one).reduce().is_zero():
        return ctx, None, None, f'merge broke Sum w = 1 (float {W.to_float():.6f})'

    guard = None
    if new_pending:
        guard = 'uncertified leaves: ' + ','.join(
            sorted({q['reason'] for q in new_pending}))
    store(P, Q, max_depth, rows, guard, new_pending)
    return ctx, rows, new_pending, guard


def stats():
    """(rows cached, certified, guard-failed) -- for a quick 'what do I already have'."""
    if not os.path.isdir(CACHE_DIR):
        return 0, 0, 0
    n = c = 0
    for fn in os.listdir(CACHE_DIR):
        if not fn.endswith('.json'):
            continue
        try:
            with open(os.path.join(CACHE_DIR, fn)) as f:
                rec = json.load(f)
        except Exception:                                    # noqa: BLE001
            continue
        n += 1
        c += rec.get('guard') is None
    return n, c, n - c


if __name__ == '__main__':
    n, c, gf = stats()
    print(f'{CACHE_DIR}: {n} rows cached ({c} certified, {gf} guard-failed), tag {CACHE_TAG}')
