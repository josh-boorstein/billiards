#!/usr/bin/env python3
"""
word_tree.py — Build the partition by enumerating reflection words.

Binary tree of reflection sequences. At each node, the beam can hit
one of two sides. The split point (boundary) is where the accumulated
position polynomial crosses 0 (vertex R) or 1 (vertex A) or the
x-coordinate crosses 0 (vertex O).

No orbit tracing needed — purely algebraic computation from the
reflection geometry.
"""

import math


def reflect_point_across_line(px, py, ax, ay, bx, by):
    """Reflect (px,py) across line through (ax,ay)-(bx,by)."""
    dx, dy = bx - ax, by - ay
    d2 = dx * dx + dy * dy
    if d2 < 1e-30:
        return px, py
    t = ((px - ax) * dx + (py - ay) * dy) / d2
    return 2.0 * (ax + t * dx) - px, 2.0 * (ay + t * dy) - py


class TriangleState:
    """Track the three vertex positions through reflections."""

    def __init__(self, alpha):
        self.alpha = alpha
        cot_a = math.cos(alpha) / math.sin(alpha)
        self.Ox, self.Oy = 0.0, 0.0
        self.Rx, self.Ry = cot_a, 0.0
        self.Ax, self.Ay = cot_a, 1.0

    def copy(self):
        s = TriangleState.__new__(TriangleState)
        s.alpha = self.alpha
        s.Ox, s.Oy = self.Ox, self.Oy
        s.Rx, s.Ry = self.Rx, self.Ry
        s.Ax, s.Ay = self.Ax, self.Ay
        return s

    def reflect(self, side):
        """Reflect the appropriate vertex across the current side."""
        if side == 'H':
            ax, ay, bx, by = self.Ox, self.Oy, self.Ax, self.Ay
            self.Rx, self.Ry = reflect_point_across_line(
                self.Rx, self.Ry, ax, ay, bx, by)
        elif side == 'L2':
            ax, ay, bx, by = self.Ox, self.Oy, self.Rx, self.Ry
            self.Ax, self.Ay = reflect_point_across_line(
                self.Ax, self.Ay, ax, ay, bx, by)
        elif side == 'L1':
            ax, ay, bx, by = self.Rx, self.Ry, self.Ax, self.Ay
            self.Ox, self.Oy = reflect_point_across_line(
                self.Ox, self.Oy, ax, ay, bx, by)

    def vertex_y(self, vtype):
        if vtype == 'O': return self.Oy
        if vtype == 'R': return self.Ry
        if vtype == 'A': return self.Ay

    def vertex_x(self, vtype):
        if vtype == 'O': return self.Ox
        if vtype == 'R': return self.Rx
        if vtype == 'A': return self.Ax


def next_sides_after(side):
    """Which sides can the beam hit next after reflecting off `side`?

    After each reflection, the beam can't hit the same side again.
    It can hit either of the other two sides.
    """
    if side == 'H':
        return ['L2', 'L1']
    elif side == 'L2':
        return ['H', 'L1']
    elif side == 'L1':
        return ['H', 'L2']
    return []


def find_split_point(state, next_a, next_b):
    """Find the s value where the beam switches from hitting next_a to next_b.

    The split occurs at a vertex encounter. The vertex between the two
    candidate sides determines the split type:

    Between L₂ and L₁ after H: vertex R (where L₂ meets... no, L₁ and L₂
      don't share a vertex. The split is at vertex R for L₂→L₁ or at
      vertex A for the other direction.)

    Actually, the split between two sides is at the vertex NOT shared
    by either candidate side with the reflection side.

    After H (= O-A): candidates L₂ (= O-R) and L₁ (= R-A).
      L₂ shares O with H. L₁ shares A with H.
      The vertex between L₂ and L₁ is R. Split at R: y_R of current state.

    After L₂ (= O-R): candidates H (= O-A) and L₁ (= R-A).
      H shares O with L₂. L₁ shares R with L₂.
      The vertex between H and L₁ is A. Split at A: y_A of current state.

    After L₁ (= R-A): candidates H (= O-A).
      Only one candidate, no split needed.
    """
    # The split vertex is the one shared by the two candidate sides
    # but not by the previous side.
    #
    # More precisely: candidates next_a and next_b. The beam goes from
    # the current triangle through one side or the other. The vertex
    # that separates the two regions is the one that both candidates
    # "compete" for.
    #
    # After H: L₂ takes beams below R, L₁ takes beams above R.
    #   Split at y = R.y in the unfolded frame.
    #   Beam height s < R.y → hits L₂ next. s > R.y → hits L₁ next.
    #   (Or vice versa, depending on the geometry.)
    #
    # After L₂: H takes beams below A, L₁ takes beams above A.
    #   Split at y = A.y.

    sides = {next_a, next_b}
    if sides == {'L2', 'L1'}:
        return state.Ry, 'R'  # R is between L₂ and L₁
    elif sides == {'H', 'L1'}:
        return state.Ay, 'A'  # A is between H and L₁
    elif sides == {'H', 'L2'}:
        return state.Oy, 'O'  # O is between H and L₂
    return None, None


def build_partition_tree(alpha, max_depth=200, min_width=1e-10,
                         closure_tol=1e-6):
    """Build the partition by recursive word enumeration.

    Returns list of (lo, hi, word) for each branch.

    closure_tol: tolerance for orbit closure detection. At rational
        α=pπ/q the triangle returns exactly, so 1e-6 (default) works.
        At irrational α the triangle never returns exactly — pass None
        to disable closure detection (tree becomes purely depth-limited)
        or use a very tight value like 1e-12.
    """
    from right_triangle_billiards import RightTriangleBilliard
    bil_ref = [RightTriangleBilliard(alpha)]
    alpha_ref = [alpha]

    branches = []
    boundaries = set()

    def zoom_bisect(lo, hi):
        """For narrow intervals, find boundaries by tracing + bisection."""
        bil = bil_ref[0]
        mh = max_depth + 10
        s_mid = 0.5 * (lo + hi)
        P_lo, P_mid, P_hi = None, None, None
        for s, label in [(lo + (hi-lo)*0.1, 'lo'), (s_mid, 'mid'), (hi - (hi-lo)*0.1, 'hi')]:
            try:
                tr = bil.trace(s, max_hits=mh, stop_at_perpendicular_return=True)
                if tr.returned_perpendicular:
                    if label == 'lo': P_lo = len(tr.hits)
                    elif label == 'mid': P_mid = len(tr.hits)
                    else: P_hi = len(tr.hits)
            except Exception:
                pass
        # If periods differ within the interval, there's a boundary — bisect
        if P_lo and P_hi and P_lo != P_hi:
            left, right = lo + (hi-lo)*0.1, hi - (hi-lo)*0.1
            for _ in range(50):
                mid = 0.5 * (left + right)
                try:
                    tr = bil.trace(mid, max_hits=mh, stop_at_perpendicular_return=True)
                    P = len(tr.hits) if tr.returned_perpendicular else 0
                except Exception:
                    P = 0
                if P == P_lo:
                    left = mid
                else:
                    right = mid
            boundaries.add(0.5 * (left + right))

    def recurse(lo, hi, word, state, prev_side, depth):
        precision_limit = max(min_width, depth * 2e-14)
        if hi - lo < precision_limit or depth > max_depth:
            if hi - lo > min_width and hi - lo < 1e-6:
                # Narrow but above min_width: try zoom bisection
                zoom_bisect(lo, hi)
            branches.append((lo, hi, list(word)))
            return

        # Check for orbit closure: triangle returns to original SHAPE
        # (relative positions preserved, up to horizontal translation
        # AND possible mirror reflection for odd-period orbits)
        if closure_tol is not None:
            cot_a = math.cos(alpha_ref[0]) / math.sin(alpha_ref[0])
            dx_RO = state.Rx - state.Ox
            dy_RO = state.Ry - state.Oy
            dx_AO = state.Ax - state.Ox
            dy_AO = state.Ay - state.Oy
            tol = closure_tol
            original = (abs(dx_RO - cot_a) < tol and abs(dy_RO) < tol and
                        abs(dx_AO - cot_a) < tol and abs(dy_AO - 1.0) < tol)
            mirrored = (abs(dx_RO + cot_a) < tol and abs(dy_RO) < tol and
                        abs(dx_AO + cot_a) < tol and abs(dy_AO - 1.0) < tol)
        else:
            original = False
            mirrored = False
        if depth > 3 and (original or mirrored):
            branches.append((lo, hi, list(word)))
            return

        candidates = next_sides_after(prev_side)

        if len(candidates) == 1:
            side = candidates[0]
            new_state = state.copy()
            new_state.reflect(side)
            word.append(side)

            # Check for perpendicular return: after L₁, the beam goes
            # left again. If the accumulated rotation brings it back to
            # perpendicular, the orbit closes.
            # Heuristic: check if the beam has returned to near-perpendicular
            # by checking if the triangle is close to its original orientation.
            # For now, just limit depth.

            recurse(lo, hi, word, new_state, side, depth + 1)
            word.pop()
            return

        # Two candidates: find the split point
        next_a, next_b = candidates
        split_y, split_vertex = find_split_point(state, next_a, next_b)

        if split_y is None:
            branches.append((lo, hi, list(word)))
            return

        # CASE 1: Split outside interval — no branching needed
        if split_y <= lo + 1e-12 or split_y >= hi - 1e-12:
            # All beams go to one side. Determine which by testing
            # a point in the middle of the interval.
            s_mid = 0.5 * (lo + hi)
            hit_idx = len(word)
            try:
                tr_mid = bil_ref[0].trace(
                    s_mid, max_hits=hit_idx + 3,
                    stop_at_perpendicular_return=False)
            except Exception:
                branches.append((lo, hi, list(word)))
                return
            if hit_idx >= len(tr_mid.hits):
                branches.append((lo, hi, list(word)))
                return
            the_side = tr_mid.hits[hit_idx].side
            if the_side not in (next_a, next_b):
                branches.append((lo, hi, list(word)))
                return
            new_state = state.copy()
            new_state.reflect(the_side)
            word.append(the_side)
            recurse(lo, hi, word, new_state, the_side, depth + 1)
            word.pop()
            return

        # CASE 2: Split inside interval — determine ordering with tracer
        boundaries.add(split_y)
        hit_idx = len(word)
        # Use offset proportional to interval width (not fixed 1e-8)
        offset = max(0.01 * (hi - lo), 1e-13)
        if split_y - offset > lo:
            s_test = split_y - offset
        elif split_y + offset < hi:
            s_test = split_y + offset
        else:
            # Interval too narrow to place test point reliably
            branches.append((lo, hi, list(word)))
            return
        try:
            tr_test = bil_ref[0].trace(
                s_test, max_hits=hit_idx + 3,
                stop_at_perpendicular_return=False)
        except Exception:
            branches.append((lo, hi, list(word)))
            return
        if hit_idx >= len(tr_test.hits):
            branches.append((lo, hi, list(word)))
            return

        test_side = tr_test.hits[hit_idx].side
        if test_side == next_a:
            side_lo, side_hi = (next_a, next_b) if s_test < split_y else (next_b, next_a)
        elif test_side == next_b:
            side_lo, side_hi = (next_b, next_a) if s_test < split_y else (next_a, next_b)
        else:
            branches.append((lo, hi, list(word)))
            return

        # Lower interval [lo, split_y]
        state_lo = state.copy()
        state_lo.reflect(side_lo)
        word.append(side_lo)
        recurse(lo, split_y, word, state_lo, side_lo, depth + 1)
        word.pop()

        # Upper interval [split_y, hi]
        state_hi = state.copy()
        state_hi.reflect(side_hi)
        word.append(side_hi)
        recurse(split_y, hi, word, state_hi, side_hi, depth + 1)
        word.pop()

    # Start: beam from L₁ going left, first hit is H
    init_state = TriangleState(alpha)
    init_state.reflect('H')
    recurse(0.0, 1.0, ['H'], init_state, 'H', 1)

    # Deduplicate near-identical boundaries (floating-point jitter)
    sorted_bounds = sorted(boundaries)
    deduped = []
    for b in sorted_bounds:
        if not deduped or b - deduped[-1] > min_width:
            deduped.append(b)

    return branches, deduped


def build_partition_tree_exact(alpha, max_depth=500, min_width=1e-14,
                               closure_tol=1e-6):
    """Build partition using exact polynomial arithmetic.

    Uses ExactTriangleState (CosPoly in Z[1/2][cos(2α)]) for split-point
    computation. No floating-point precision floor — can drill arbitrarily
    deep into fine structure.

    Returns list of (lo, hi, word, split_poly) for each branch, plus
    boundary list. split_poly is the CosPoly giving the exact boundary position.
    """
    from right_triangle_billiards import RightTriangleBilliard
    from exact_state import ExactTriangleState

    bil = RightTriangleBilliard(alpha)
    branches = []
    boundaries = []  # list of (float_val, CosPoly, vertex_type)

    def recurse(lo, hi, word, state, prev_side, depth):
        if hi - lo < min_width or depth > max_depth:
            branches.append((lo, hi, list(word)))
            return

        # Orbit closure check via exact polynomial comparison
        if closure_tol is not None and depth > 3:
            cot_a = math.cos(alpha) / math.sin(alpha)
            # Check if relative vertex positions match original
            # In exact repr: Rx_c - Ox_c should be 1, Ry-Oy should be 0, etc.
            dx_c = state.Rx_c - state.Ox_c
            dy = state.Ry - state.Oy
            dax_c = state.Ax_c - state.Ox_c
            day = state.Ay - state.Oy

            # Original: dx_c=1, dy=0, dax_c=1, day=1
            # Check by evaluating (cheap — just polynomial evaluation)
            dx_val = dx_c.evaluate(alpha)
            dy_val = dy.evaluate(alpha)
            dax_val = dax_c.evaluate(alpha)
            day_val = day.evaluate(alpha)

            tol = closure_tol
            original = (abs(dx_val - 1) < tol and abs(dy_val) < tol and
                        abs(dax_val - 1) < tol and abs(day_val - 1) < tol)
            mirrored = (abs(dx_val + 1) < tol and abs(dy_val) < tol and
                        abs(dax_val + 1) < tol and abs(day_val - 1) < tol)

            if original or mirrored:
                branches.append((lo, hi, list(word)))
                return

        candidates = next_sides_after(prev_side)

        if len(candidates) == 1:
            side = candidates[0]
            new_state = state.copy()
            new_state.reflect(side)
            word.append(side)
            recurse(lo, hi, word, new_state, side, depth + 1)
            word.pop()
            return

        # Two candidates: find the split point
        next_a, next_b = candidates
        sides = {next_a, next_b}
        if sides == {'L2', 'L1'}:
            split_poly = state.vertex_y('R')
            split_vertex = 'R'
        elif sides == {'H', 'L1'}:
            split_poly = state.vertex_y('A')
            split_vertex = 'A'
        elif sides == {'H', 'L2'}:
            split_poly = state.vertex_y('O')
            split_vertex = 'O'
        else:
            branches.append((lo, hi, list(word)))
            return

        split_y = split_poly.evaluate(alpha)

        # CASE 1: Split outside interval
        if split_y <= lo + min_width or split_y >= hi - min_width:
            # Determine which side by testing midpoint
            s_mid = 0.5 * (lo + hi)
            hit_idx = len(word)
            try:
                tr_mid = bil.trace(s_mid, max_hits=hit_idx + 3,
                                   stop_at_perpendicular_return=False)
            except Exception:
                branches.append((lo, hi, list(word)))
                return
            if hit_idx >= len(tr_mid.hits):
                branches.append((lo, hi, list(word)))
                return
            the_side = tr_mid.hits[hit_idx].side
            if the_side not in (next_a, next_b):
                branches.append((lo, hi, list(word)))
                return
            new_state = state.copy()
            new_state.reflect(the_side)
            word.append(the_side)
            recurse(lo, hi, word, new_state, the_side, depth + 1)
            word.pop()
            return

        # CASE 2: Split inside interval
        boundaries.append((split_y, split_poly.copy(), split_vertex))
        hit_idx = len(word)

        # Determine ordering: which side is below split, which above
        offset = max(0.01 * (hi - lo), min_width)
        if split_y - offset > lo:
            s_test = split_y - offset
        elif split_y + offset < hi:
            s_test = split_y + offset
        else:
            branches.append((lo, hi, list(word)))
            return
        try:
            tr_test = bil.trace(s_test, max_hits=hit_idx + 3,
                                stop_at_perpendicular_return=False)
        except Exception:
            branches.append((lo, hi, list(word)))
            return
        if hit_idx >= len(tr_test.hits):
            branches.append((lo, hi, list(word)))
            return

        test_side = tr_test.hits[hit_idx].side
        if test_side == next_a:
            side_lo, side_hi = (next_a, next_b) if s_test < split_y else (next_b, next_a)
        elif test_side == next_b:
            side_lo, side_hi = (next_b, next_a) if s_test < split_y else (next_a, next_b)
        else:
            branches.append((lo, hi, list(word)))
            return

        # Lower interval [lo, split_y]
        state_lo = state.copy()
        state_lo.reflect(side_lo)
        word.append(side_lo)
        recurse(lo, split_y, word, state_lo, side_lo, depth + 1)
        word.pop()

        # Upper interval [split_y, hi]
        state_hi = state.copy()
        state_hi.reflect(side_hi)
        word.append(side_hi)
        recurse(split_y, hi, word, state_hi, side_hi, depth + 1)
        word.pop()

    # Start: beam from L₁ going left, first hit is H
    init_state = ExactTriangleState(alpha)
    init_state.reflect('H')
    recurse(0.0, 1.0, ['H'], init_state, 'H', 1)

    # Sort and deduplicate boundaries
    boundaries.sort(key=lambda x: x[0])
    deduped = []
    deduped_polys = []
    for val, poly, vtype in boundaries:
        if not deduped or val - deduped[-1] > min_width:
            deduped.append(val)
            deduped_polys.append((poly, vtype))

    return branches, deduped, deduped_polys


def main():
    print("WORD TREE: ALGEBRAIC PARTITION COMPUTATION")
    print("=" * 70)

    from right_triangle_billiards import RightTriangleBilliard

    for p, q in [(1, 5), (2, 5), (1, 7), (3, 7), (1, 11), (5, 11),
                 (2, 11), (3, 11)]:
        alpha = p * math.pi / q
        if alpha >= math.pi / 2:
            continue

        branches, tree_bounds = build_partition_tree(alpha, max_depth=80)

        # Compare with tracer
        bil = RightTriangleBilliard(alpha)
        traced = bil.discover_periodic_branches(n_samples=600, max_hits=15000)
        actual_bounds = sorted([b[1] for b in traced[:-1]])

        # Match
        matched = 0
        for ab in actual_bounds:
            if any(abs(ab - tb) < 1e-5 for tb in tree_bounds):
                matched += 1

        print(f"\np={p} q={q}: tree found {len(tree_bounds)} boundaries, "
              f"tracer found {len(actual_bounds)}, "
              f"matched {matched}/{len(actual_bounds)}")

        if len(actual_bounds) <= 12:
            print(f"  Tree:   {[f'{b:.5f}' for b in tree_bounds[:15]]}")
            print(f"  Actual: {[f'{b:.5f}' for b in actual_bounds]}")

        # Show tree branches (first few)
        if len(branches) <= 20:
            print(f"  Branches ({len(branches)}):")
            for lo, hi, word in branches[:10]:
                w = ''.join(s[0] for s in word)
                if len(w) > 25:
                    w = w[:25] + '..'
                print(f"    [{lo:.6f}, {hi:.6f}] word={w}")


if __name__ == "__main__":
    main()
