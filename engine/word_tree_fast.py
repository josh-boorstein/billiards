#!/usr/bin/env python3
"""
word_tree_fast.py — Iterative partition tree with numpy-accelerated exact arithmetic.

Key improvements over word_tree.build_partition_tree_exact:
- Iterative (stack-based) instead of recursive: no recursion limit
- Uses cos_poly_fast (numpy backend): ~10-20x faster per reflection
- Same algorithm and results as the recursive version
"""

import math
from right_triangle_billiards import RightTriangleBilliard
from exact_state_fast import ExactTriangleState


SIDES_AFTER = {'H': ('L2', 'L1'), 'L2': ('H', 'L1'), 'L1': ('H', 'L2')}

SPLIT_VERTEX = {
    frozenset(('L2', 'L1')): 'R',
    frozenset(('H', 'L1')): 'A',
    frozenset(('H', 'L2')): 'O',
}


def build_partition_tree(alpha, max_depth=5000, min_width=1e-15, closure_tol=None):
    """Build partition using iterative exact tree.

    Returns (branches, boundaries, boundary_polys) where:
      branches: list of (lo, hi, word_length)
      boundaries: list of float boundary positions
      boundary_polys: list of (CosPoly, vertex_type) for each boundary
    """
    bil = RightTriangleBilliard(alpha)

    branches = []
    boundaries = []
    boundary_polys = []

    # Stack: (lo, hi, state, prev_side, depth, word_len)
    init_state = ExactTriangleState(alpha)
    init_state.reflect('H')
    stack = [(0.0, 1.0, init_state, 'H', 1, 1)]

    while stack:
        lo, hi, state, prev_side, depth, word_len = stack.pop()

        if hi - lo < min_width or depth > max_depth:
            branches.append((lo, hi, word_len))
            continue

        # Orbit closure check
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

        # Two candidates: find split point
        next_a, next_b = candidates
        vertex = SPLIT_VERTEX[frozenset(candidates)]
        split_poly = state.vertex_y(vertex)
        split_y = split_poly.evaluate(alpha)

        # Split outside interval — determine which side using cached ordering
        if split_y <= lo + min_width or split_y >= hi - min_width:
            # Need to know which side the interval is on.
            # Use a single trace at the midpoint to determine.
            s_mid = 0.5 * (lo + hi)
            try:
                tr_mid = bil.trace(s_mid, max_hits=word_len + 3,
                                   stop_at_perpendicular_return=False)
            except Exception:
                branches.append((lo, hi, word_len))
                continue
            if word_len >= len(tr_mid.hits):
                branches.append((lo, hi, word_len))
                continue
            the_side = tr_mid.hits[word_len].side
            if the_side not in (next_a, next_b):
                branches.append((lo, hi, word_len))
                continue
            new_state = state.copy()
            new_state.reflect(the_side)
            stack.append((lo, hi, new_state, the_side, depth + 1, word_len + 1))
            continue

        # Split inside interval — record boundary and determine ordering
        boundaries.append(split_y)
        boundary_polys.append((split_poly.copy(), vertex))

        # Determine ordering with one trace call
        offset = max(0.01 * (hi - lo), min_width)
        s_test = split_y - offset if split_y - offset > lo else split_y + offset
        try:
            tr_test = bil.trace(s_test, max_hits=word_len + 3,
                                stop_at_perpendicular_return=False)
        except Exception:
            branches.append((lo, hi, word_len))
            continue
        if word_len >= len(tr_test.hits):
            branches.append((lo, hi, word_len))
            continue

        test_side = tr_test.hits[word_len].side
        if test_side == next_a:
            side_lo, side_hi = (next_a, next_b) if s_test < split_y else (next_b, next_a)
        elif test_side == next_b:
            side_lo, side_hi = (next_b, next_a) if s_test < split_y else (next_a, next_b)
        else:
            branches.append((lo, hi, word_len))
            continue

        # Push UPPER first (so lower is processed first from stack)
        state_hi = state.copy()
        state_hi.reflect(side_hi)
        stack.append((split_y, hi, state_hi, side_hi, depth + 1, word_len + 1))

        state_lo = state.copy()
        state_lo.reflect(side_lo)
        stack.append((lo, split_y, state_lo, side_lo, depth + 1, word_len + 1))

    # Sort and deduplicate boundaries
    combined = sorted(zip(boundaries, boundary_polys), key=lambda x: x[0])
    deduped_bounds = []
    deduped_polys = []
    for val, poly_info in combined:
        if not deduped_bounds or val - deduped_bounds[-1] > min_width:
            deduped_bounds.append(val)
            deduped_polys.append(poly_info)

    return branches, deduped_bounds, deduped_polys


def main():
    """Benchmark against the recursive version."""
    import time

    print("FAST ITERATIVE TREE vs RECURSIVE TREE")
    print("=" * 60)

    # Import the old version for comparison
    from word_tree import build_partition_tree_exact

    for p, q in [(1, 7), (1, 11), (2, 11), (1, 13)]:
        alpha = p * math.pi / q
        if alpha >= math.pi / 2:
            continue

        # Old recursive version
        t0 = time.time()
        _, bounds_old, _ = build_partition_tree_exact(alpha, max_depth=200)
        t_old = time.time() - t0

        # New iterative version
        t0 = time.time()
        _, bounds_new, _ = build_partition_tree(alpha, max_depth=200)
        t_new = time.time() - t0

        # Verify
        match = len(bounds_old) == len(bounds_new)
        if match:
            for bo, bn in zip(sorted(bounds_old), sorted(bounds_new)):
                if abs(bo - bn) > 1e-8:
                    match = False
                    break

        speedup = t_old / t_new if t_new > 0 else float('inf')
        print(f"  ({p},{q}): old={t_old*1000:.0f}ms  new={t_new*1000:.0f}ms  "
              f"speedup={speedup:.1f}x  bounds={len(bounds_new)}  match={'YES' if match else 'NO'}")

    # Test at irrational angle with deep tree
    print()
    print("Deep irrational (α=1):")
    for depth in [500, 1000, 2000, 5000]:
        t0 = time.time()
        branches, bounds, polys = build_partition_tree(
            1.0, max_depth=depth, min_width=1e-15, closure_tol=None)
        t = time.time() - t0
        max_deg = max((p.degree() for p, _ in polys), default=0)
        print(f"  depth={depth}: {len(bounds)} boundaries, "
              f"max_deg={max_deg}, {t:.1f}s")


if __name__ == "__main__":
    main()
