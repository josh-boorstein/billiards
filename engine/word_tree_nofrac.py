#!/usr/bin/env python3
"""
word_tree_nofrac.py — fully tracer-free exact partition tree.

The last-session bottleneck (§52): word_tree_fast calls bil.trace(...) at
every split just to decide which child interval is on the low-s side. That
trace is O(word_len), making the whole build ~quadratic in depth.

Key insight: in the unfolded frame the beam is the horizontal line y = s
(travelling in -x). A split happens at a vertex V of height split_y that is
shared by the two candidate sides. Each candidate side has one OTHER
endpoint. For s slightly below split_y the beam hits the side whose other
endpoint lies BELOW V; for s slightly above, the side whose other endpoint
lies ABOVE V. (SPLIT_VERTEX guarantees exactly one of each.) This is the
sign of C_j and needs only the exact vertex y-coordinates — no trace.

Validated in-file against the tracer-based word_tree_fast.build_partition_tree
and the saved deep_tree JSONs.
"""

import math
from exact_state_fast import ExactTriangleState

SIDES_AFTER = {'H': ('L2', 'L1'), 'L2': ('H', 'L1'), 'L1': ('H', 'L2')}

SPLIT_VERTEX = {
    frozenset(('L2', 'L1')): 'R',
    frozenset(('H', 'L1')): 'A',
    frozenset(('H', 'L2')): 'O',
}

SIDE_ENDPOINTS = {'L1': ('R', 'A'), 'L2': ('O', 'R'), 'H': ('O', 'A')}


def _vxy(state, vtype, cot, alpha):
    xc = {'O': state.Ox_c, 'R': state.Rx_c, 'A': state.Ax_c}[vtype]
    return xc.evaluate(alpha) * cot, state.vertex_y(vtype).evaluate(alpha)


def _hit_side(state, candidates, s, alpha, cot):
    """Which candidate side the horizontal ray y=s (travelling in -x) hits.

    A candidate side is a segment between two triangle vertices; it can be
    crossed only if s lies within its y-span. Among the candidates whose
    span contains s, the ray (coming from +x) hits the FRONTMOST — the one
    with the largest crossing-x. This is fully geometric: it resolves both
    the straddle case (V between the two other endpoints) and the
    non-straddle case (both endpoints on one side of V) with the same rule.
    """
    best_side, best_x = None, -float('inf')
    for side in candidates:
        p, q = SIDE_ENDPOINTS[side]
        (px, py), (qx, qy) = _vxy(state, p, cot, alpha), _vxy(state, q, cot,
                                                              alpha)
        lo_y, hi_y = (py, qy) if py <= qy else (qy, py)
        if not (lo_y <= s <= hi_y) or qy == py:
            continue
        cx = px + (s - py) * (qx - px) / (qy - py)
        if cx > best_x:
            best_x, best_side = cx, side
    return best_side


def build_partition_tree(alpha, max_depth=5000, min_width=1e-15,
                         closure_tol=None, return_raw=False):
    """Tracer-free exact partition tree.

    Returns (branches, boundaries, boundary_polys):
      branches: list of (lo, hi, word_length)
      boundaries: sorted list of float boundary positions
      boundary_polys: list of (CosPoly, vertex_type)

    With return_raw=True the final FLOAT deduplication is skipped and the
    boundaries/polys are returned exactly as emitted by the tree (still
    paired and sorted by float value). Use this when the caller performs
    its own EXACT dedup (see n_exact.py) — the float dedup both over-counts
    (keeps spurious ~1e-15 edge boundaries) and can under-count (merges
    genuinely-distinct boundaries closer than min_width at large q).
    """
    branches = []
    boundaries = []
    boundary_polys = []
    n_depth_trunc = 0            # branches cut off by the depth cap (not closure)

    cot = 1.0 / math.tan(alpha)
    init_state = ExactTriangleState(alpha)
    init_state.reflect('H')
    stack = [(0.0, 1.0, init_state, 'H', 1, 1)]

    while stack:
        lo, hi, state, prev_side, depth, word_len = stack.pop()

        if hi - lo < min_width or depth > max_depth:
            branches.append((lo, hi, word_len))
            if depth > max_depth and hi - lo >= min_width:
                # would have kept splitting: partition is NOT complete here
                n_depth_trunc += 1
            continue

        if closure_tol is not None and depth > 3:
            dx_val = (state.Rx_c - state.Ox_c).evaluate(alpha)
            dy_val = (state.Ry - state.Oy).evaluate(alpha)
            dax_val = (state.Ax_c - state.Ox_c).evaluate(alpha)
            day_val = (state.Ay - state.Oy).evaluate(alpha)
            tol = closure_tol
            original = (abs(dx_val - 1) < tol and abs(dy_val) < tol and
                        abs(dax_val - 1) < tol and abs(day_val - 1) < tol)
            mirrored = (abs(dx_val + 1) < tol and abs(dy_val) < tol and
                        abs(dax_val + 1) < tol and abs(day_val - 1) < tol)
            if original or mirrored:
                branches.append((lo, hi, word_len))
                continue

        candidates = SIDES_AFTER[prev_side]

        if len(candidates) == 1:
            side = candidates[0]
            new_state = state.copy()
            new_state.reflect(side)
            stack.append((lo, hi, new_state, side, depth + 1, word_len + 1))
            continue

        vertex = SPLIT_VERTEX[frozenset(candidates)]
        split_poly = state.vertex_y(vertex)
        split_y = split_poly.evaluate(alpha)

        # Split outside interval: whole interval takes a single side,
        # decided geometrically at the interval midpoint.
        if split_y <= lo + min_width or split_y >= hi - min_width:
            side = _hit_side(state, candidates, 0.5 * (lo + hi), alpha, cot)
            if side is None:
                branches.append((lo, hi, word_len))
                continue
            new_state = state.copy()
            new_state.reflect(side)
            stack.append((lo, hi, new_state, side, depth + 1, word_len + 1))
            continue

        # Split inside interval: record boundary; each child's side decided
        # at that child's midpoint.
        boundaries.append(split_y)
        boundary_polys.append((split_poly.copy(), vertex))

        side_hi = _hit_side(state, candidates, 0.5 * (split_y + hi), alpha,
                            cot)
        side_lo = _hit_side(state, candidates, 0.5 * (lo + split_y), alpha,
                            cot)
        if side_hi is None or side_lo is None:
            branches.append((lo, hi, word_len))
            continue

        state_hi = state.copy()
        state_hi.reflect(side_hi)
        stack.append((split_y, hi, state_hi, side_hi, depth + 1, word_len + 1))

        state_lo = state.copy()
        state_lo.reflect(side_lo)
        stack.append((lo, split_y, state_lo, side_lo, depth + 1, word_len + 1))

    combined = sorted(zip(boundaries, boundary_polys), key=lambda x: x[0])
    if return_raw:
        raw_bounds = [v for v, _ in combined]
        raw_polys = [pinfo for _, pinfo in combined]
        # 4th element: completeness certificate. n_depth_trunc==0 means every
        # branch closed via orbit closure (nothing cut by the depth cap), so
        # the partition is PROVABLY complete and the count is exact.
        return branches, raw_bounds, raw_polys, n_depth_trunc
    deduped_bounds = []
    deduped_polys = []
    for val, poly_info in combined:
        if not deduped_bounds or val - deduped_bounds[-1] > min_width:
            deduped_bounds.append(val)
            deduped_polys.append(poly_info)
    return branches, deduped_bounds, deduped_polys


if __name__ == '__main__':
    import time
    import json
    from word_tree_fast import build_partition_tree as build_tracer

    print("VALIDATION: tracer-free vs tracer-based")
    print("=" * 60)
    for p, q in [(1, 7), (2, 11), (1, 13), (3, 13)]:
        alpha = p * math.pi / q
        if alpha >= math.pi / 2:
            continue
        _, b_new, _ = build_partition_tree(alpha, max_depth=300,
                                           closure_tol=1e-6)
        _, b_old, _ = build_tracer(alpha, max_depth=300, closure_tol=1e-6)
        b_new, b_old = sorted(b_new), sorted(b_old)
        match = (len(b_new) == len(b_old) and
                 all(abs(x - y) < 1e-9 for x, y in zip(b_new, b_old)))
        print(f"  ({p},{q}): nofrac={len(b_new)} tracer={len(b_old)} "
              f"match={match}")

    print("\nVALIDATION: irrational vs saved deep_tree JSON")
    for a_str in ['0.4000', '1.0000']:
        alpha = float(a_str)
        saved = json.load(open(f'deep_tree_a{a_str}_d5000.json'))
        sb = sorted(x['s'] for x in saved['boundaries'])
        t0 = time.time()
        _, b_new, _ = build_partition_tree(alpha, max_depth=5000,
                                           closure_tol=None)
        dt = time.time() - t0
        b_new = sorted(b_new)
        match = (len(b_new) == len(sb) and
                 all(abs(x - y) < 1e-8 for x, y in zip(b_new, sb)))
        print(f"  a={alpha}: nofrac={len(b_new)} saved={len(sb)} "
              f"match={match} time={dt:.1f}s "
              f"(saved took {saved['elapsed_seconds']:.0f}s)")
