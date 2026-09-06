#!/usr/bin/env python3
"""
orphan_region.py — region-restricted, TRACER-FREE orphan finder + exact width.

Combines word_tree_nofrac (fully algebraic split-ordering, no float tracer, so
no orbit-length error accumulation) with region pruning (only explore a window
[s_lo, s_hi] around the orphan). This is the exact tool the §81 golden-orphan
measurement needed: the archived region_tree.py still called bil.trace() for
side decisions, the source of the head-on float degradation at 55/89, and the
full certified() tree is too slow at golden angles (marginal orphan ⇒ deep
completeness over the WHOLE interval).

Orphan identification is INTRINSIC and exact: during development every branch
carries the minimum "verticality gap" of any interior (non-L1) side it hits.
An interior head-on hit is a vertical wall-copy (Lemma D2), so the orphan
branch has min_gap ≈ 0 (machine zero) while paired branches have min_gap ~ 0.05
- 0.3 — a ~12-order-of-magnitude separation. (L1 is excluded: an L1 head-on is
the trivial launch/return, not a retrace, Lemma 3.)

Boundaries are exact CosPolys; the width is the cyclotomic number
sum_k c_k cos(k P pi / Q)/2^shift evaluated to arbitrary precision (mpmath),
immune to floating error.
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
    best_side, best_x = None, -float('inf')
    for side in candidates:
        p, q = SIDE_ENDPOINTS[side]
        (px, py), (qx, qy) = _vxy(state, p, cot, alpha), _vxy(state, q, cot, alpha)
        lo_y, hi_y = (py, qy) if py <= qy else (qy, py)
        if not (lo_y <= s <= hi_y) or qy == py:
            continue
        cx = px + (s - py) * (qx - px) / (qy - py)
        if cx > best_x:
            best_x, best_side = cx, side
    return best_side


def _vgap(state, side, cot, alpha):
    """Verticality gap of the copy of `side`: |x1 - x2| of its endpoints.
    L1 excluded (its head-on is the trivial launch/return)."""
    if side == 'L1':
        return 9.9
    p, q = SIDE_ENDPOINTS[side]
    return abs(_vxy(state, p, cot, alpha)[0] - _vxy(state, q, cot, alpha)[0])


def find_orphan(P, Q, s_lo, s_hi, max_depth=40000, min_width=1e-16,
                closure_tol=1e-6):
    """Region-restricted tracer-free orphan finder.

    Returns dict with:
      lo, hi          float bounds of the orphan branch
      lo_poly,hi_poly exact CosPoly split points bounding it
      min_gap         interior-vertical gap of the orphan (≈ machine zero)
      runner_gap      smallest min_gap among the OTHER in-window branches
                      (sanity: should be ~0.01-0.3, i.e. clearly separated)
      n_branches      branches overlapping the window
    or None if no branch has a near-zero gap (orphan not in window).
    """
    alpha = (P / Q) * (math.pi / 2)
    cot = 1.0 / math.tan(alpha)

    def overlaps(lo, hi):
        return hi > s_lo and lo < s_hi

    boundaries = {}            # float -> (CosPoly, vertex)
    leaves = []                # (lo, hi, min_gap)

    init = ExactTriangleState(alpha)
    init.reflect('H')
    stack = [(0.0, 1.0, init, 'H', 1, 9.9)]

    while stack:
        lo, hi, state, prev, depth, mg = stack.pop()
        if not overlaps(lo, hi):
            continue
        if hi - lo < min_width or depth > max_depth:
            leaves.append((lo, hi, mg))
            continue
        if depth > 3:
            dx = (state.Rx_c - state.Ox_c).evaluate(alpha)
            dy = (state.Ry - state.Oy).evaluate(alpha)
            dax = (state.Ax_c - state.Ox_c).evaluate(alpha)
            day = (state.Ay - state.Oy).evaluate(alpha)
            t = closure_tol
            if ((abs(dx - 1) < t and abs(dy) < t and abs(dax - 1) < t and abs(day - 1) < t)
                    or (abs(dx + 1) < t and abs(dy) < t and abs(dax + 1) < t and abs(day - 1) < t)):
                leaves.append((lo, hi, mg))
                continue

        cand = SIDES_AFTER[prev]
        if len(cand) == 1:
            side = cand[0]
            g = min(mg, _vgap(state, side, cot, alpha))
            ns = state.copy(); ns.reflect(side)
            stack.append((lo, hi, ns, side, depth + 1, g))
            continue

        vtx = SPLIT_VERTEX[frozenset(cand)]
        sp = state.vertex_y(vtx)
        sy = sp.evaluate(alpha)

        if sy <= lo + min_width or sy >= hi - min_width:
            side = _hit_side(state, cand, 0.5 * (lo + hi), alpha, cot)
            if side is None:
                leaves.append((lo, hi, mg)); continue
            g = min(mg, _vgap(state, side, cot, alpha))
            ns = state.copy(); ns.reflect(side)
            stack.append((lo, hi, ns, side, depth + 1, g))
            continue

        if s_lo <= sy <= s_hi:
            boundaries[sy] = (sp.copy(), vtx)

        shi = _hit_side(state, cand, 0.5 * (sy + hi), alpha, cot)
        slo = _hit_side(state, cand, 0.5 * (lo + sy), alpha, cot)
        if shi is None or slo is None:
            leaves.append((lo, hi, mg)); continue
        if overlaps(sy, hi):
            g = min(mg, _vgap(state, shi, cot, alpha))
            s2 = state.copy(); s2.reflect(shi)
            stack.append((sy, hi, s2, shi, depth + 1, g))
        if overlaps(lo, sy):
            g = min(mg, _vgap(state, slo, cot, alpha))
            s2 = state.copy(); s2.reflect(slo)
            stack.append((lo, sy, s2, slo, depth + 1, g))

    leaves = [lf for lf in leaves if lf[1] - lf[0] > 1e-13]
    if not leaves:
        return None
    # gap≈0 leaves = the orphan PLUS its head-on-boundary slivers; the orphan
    # is the widest of them. Paired branches have gap ~ 0.01-0.3.
    headon = [lf for lf in leaves if lf[2] < 1e-8]
    if not headon:
        return None
    olo, ohi, omg = max(headon, key=lambda x: x[1] - x[0])
    # runner = closest paired branch (gap >= 1e-8): confirms separation
    paired = [lf[2] for lf in leaves if lf[2] >= 1e-8]
    runner = min(paired) if paired else 9.9

    # match exact polys to the orphan's float bounds
    def match(val):
        best, bd = None, 1e-9
        for b, poly in boundaries.items():
            if abs(b - val) < bd:
                bd = abs(b - val); best = poly
        return best[0] if best else None

    return {
        'lo': olo, 'hi': ohi, 'min_gap': omg, 'runner_gap': runner,
        'lo_poly': match(olo), 'hi_poly': match(ohi),
        'n_branches': len(leaves), 'P': P, 'Q': Q,
    }


def poly_hp(poly, P, Q, dps=50):
    """High-precision cyclotomic value: sum_k c_k cos(k P pi / Q)/2^shift."""
    from mpmath import mp, mpf, cos, pi
    mp.dps = dps
    tot = mpf(0)
    for k, c in enumerate(poly.coeffs):
        if c:
            tot += int(c) * cos(mpf(k) * P * pi / Q)
    return tot / (mpf(2) ** int(poly.shift))


def orphan_width_hp(P, Q, s_lo, s_hi, dps=50, **kw):
    """Exact orphan width (mpmath) + diagnostics."""
    r = find_orphan(P, Q, s_lo, s_hi, **kw)
    if r is None or r['lo_poly'] is None or r['hi_poly'] is None:
        return None, r
    w = poly_hp(r['hi_poly'], P, Q, dps) - poly_hp(r['lo_poly'], P, Q, dps)
    return w, r


# ===========================================================================
# RING BACKEND — depth-independent fast path for DEEP (golden) orphans.
#
# Same region-restricted descent, but vertex coordinates live in Z[ζ_{2Q}]
# (exact_state_ring.RingTriangleState): O(D²)/reflect regardless of depth, so
# depth 20000+ (55/89, 89/144 head-on) is feasible. Head-on detection is EXACT
# — an interior side copy is vertical iff its two x-coefficient RingElems are
# equal (ring-zero difference), no float threshold. Float evaluation from the
# reduced ring vector is numerically stable (bounded entries), so ordering
# decisions no longer degrade with depth.
# ===========================================================================
def _vxy_ring(state, vtype, cot):
    return state.vertex_xc(vtype).to_float() * cot, state.vertex_y(vtype).to_float()


def _hit_side_ring(state, candidates, s, cot):
    best_side, best_x = None, -float('inf')
    for side in candidates:
        p, q = SIDE_ENDPOINTS[side]
        (px, py), (qx, qy) = _vxy_ring(state, p, cot), _vxy_ring(state, q, cot)
        lo_y, hi_y = (py, qy) if py <= qy else (qy, py)
        if not (lo_y <= s <= hi_y) or qy == py:
            continue
        cx = px + (s - py) * (qx - px) / (qy - py)
        if cx > best_x:
            best_x, best_side = cx, side
    return best_side


def _vgap_ring(state, side, cot):
    if side == 'L1':
        return 9.9
    p, q = SIDE_ENDPOINTS[side]
    return abs(_vxy_ring(state, p, cot)[0] - _vxy_ring(state, q, cot)[0])


def _exact_vertical(state, side):
    """True iff the copy of `side` is EXACTLY vertical (head-on wall-copy).
    L1 excluded (its head-on is the trivial launch/return, Lemma 3)."""
    if side == 'L1':
        return False
    p, q = SIDE_ENDPOINTS[side]
    return (state.vertex_xc(p) - state.vertex_xc(q)).is_zero()


def find_orphan_ring(P, Q, s_lo, s_hi, max_depth=60000, min_width=1e-16,
                     closure_tol=1e-6, return_leaves=False, leaf_polys=False,
                     leaf_detail=False):
    """Depth-independent region-restricted orphan finder (ring backend).

    Returns the same dict shape as find_orphan, plus:
      exact_headon  True if the orphan branch hit an exactly-vertical interior
                    side copy (rigorous head-on certificate, no float threshold)
      lo_poly/hi_poly are RingElems (exact; use .to_mpf for the width).

    leaf_detail=True (s275) returns one dict per leaf instead of a tuple:
      lo, hi, min_gap, headon, depth, lo_poly, hi_poly (exact boundary RingElems,
      None at the interval ends 0/1), state (the RingTriangleState the descent
      STOPPED in) and reason, the termination cause:
        'closed'     the developed copy is a translate / x-mirror of the original
                     (perpendicular L1 return — the only CERTIFIED outcome),
        'depth'      max_depth reached, 'width' min_width reached,
        'degenerate' _hit_side_ring found no side (s273's certificate hole).
    A run is certified iff every leaf's reason == 'closed'. In a 'closed' leaf the
    perpendicular L1 hit sits at x = Rx_c·cotα, so the cell's first-return LENGTH
    is cotα·(1 − Rx_c) exactly — see probes/s275_flow_exact.py.
    """
    from exact_state_ring import RingTriangleState
    from ring_cyc import RingContext

    ctx = RingContext(P, Q)
    P, Q = ctx.P, ctx.Q
    alpha = ctx.alpha
    cot = 1.0 / math.tan(alpha)

    def overlaps(lo, hi):
        return hi > s_lo and lo < s_hi

    boundaries = {}      # float -> RingElem split point
    graze = {}           # float -> (depth, vertex) of the split that CREATED it (s349)
    leaves = []          # (lo, hi, min_gap_float, exact_headon_bool, depth, reason, state)

    def leaf(lo, hi, mg, ho, depth, reason, state):
        leaves.append((lo, hi, mg, ho, depth, reason,
                       state if leaf_detail else None))

    init = RingTriangleState(ctx)
    init.reflect('H')
    # stack entries: lo, hi, state, prev_side, depth, min_gap_float, headon_bool
    stack = [(0.0, 1.0, init, 'H', 1, 9.9, False)]

    while stack:
        lo, hi, state, prev, depth, mg, ho = stack.pop()
        if not overlaps(lo, hi):
            continue
        if hi - lo < min_width or depth > max_depth:
            leaf(lo, hi, mg, ho, depth,
                 'width' if hi - lo < min_width else 'depth', state)
            continue
        if depth > 3:
            dx = (state.Rx_c - state.Ox_c).to_float()
            dy = (state.Ry - state.Oy).to_float()
            dax = (state.Ax_c - state.Ox_c).to_float()
            day = (state.Ay - state.Oy).to_float()
            t = closure_tol
            if ((abs(dx - 1) < t and abs(dy) < t and abs(dax - 1) < t and abs(day - 1) < t)
                    or (abs(dx + 1) < t and abs(dy) < t and abs(dax + 1) < t and abs(day - 1) < t)):
                leaf(lo, hi, mg, ho, depth, 'closed', state)
                continue

        cand = SIDES_AFTER[prev]
        if len(cand) == 1:
            side = cand[0]
            g = min(mg, _vgap_ring(state, side, cot))
            h2 = ho or _exact_vertical(state, side)
            ns = state.copy(); ns.reflect(side)
            stack.append((lo, hi, ns, side, depth + 1, g, h2))
            continue

        vtx = SPLIT_VERTEX[frozenset(cand)]
        sp = state.vertex_y(vtx)
        sy = sp.to_float()

        if sy <= lo + min_width or sy >= hi - min_width:
            side = _hit_side_ring(state, cand, 0.5 * (lo + hi), cot)
            if side is None:
                leaf(lo, hi, mg, ho, depth, 'degenerate', state); continue
            g = min(mg, _vgap_ring(state, side, cot))
            h2 = ho or _exact_vertical(state, side)
            ns = state.copy(); ns.reflect(side)
            stack.append((lo, hi, ns, side, depth + 1, g, h2))
            continue

        if s_lo <= sy <= s_hi:
            boundaries[sy] = sp.copy()
            graze[sy] = (depth, vtx)

        shi = _hit_side_ring(state, cand, 0.5 * (sy + hi), cot)
        slo = _hit_side_ring(state, cand, 0.5 * (lo + sy), cot)
        if shi is None or slo is None:
            leaf(lo, hi, mg, ho, depth, 'degenerate', state); continue
        if overlaps(sy, hi):
            g = min(mg, _vgap_ring(state, shi, cot))
            h2 = ho or _exact_vertical(state, shi)
            s2 = state.copy(); s2.reflect(shi)
            stack.append((sy, hi, s2, shi, depth + 1, g, h2))
        if overlaps(lo, sy):
            g = min(mg, _vgap_ring(state, slo, cot))
            h2 = ho or _exact_vertical(state, slo)
            s2 = state.copy(); s2.reflect(slo)
            stack.append((lo, sy, s2, slo, depth + 1, g, h2))

    leaves = [lf for lf in leaves if lf[1] - lf[0] > 1e-13]

    # EXACT boundary RingElems, matched from the split points; the 0.0/1.0 outer
    # endpoints have no split poly -> None, value is 0 or 1.
    def _match_poly(val):
        best, bd = None, 1e-9
        for b, poly in boundaries.items():
            if abs(b - val) < bd:
                bd = abs(b - val); best = poly
        return best

    def _match_graze(val):
        """(depth, vertex) of the split that created this boundary; None at 0/1.

        s349: the GRAZE DEPTH is the reflection count of the developed copy whose
        vertex the boundary ray grazes -- the shallow/deep variable of s169 item 9
        and s170.  Measured from the descent, never from a numerator solve."""
        best, bd = None, 1e-9
        for b, g in graze.items():
            if abs(b - val) < bd:
                bd = abs(b - val); best = g
        return best

    if return_leaves:
        if leaf_detail:
            return [{'lo': lo, 'hi': hi, 'min_gap': mg, 'headon': ho,
                     'depth': dep, 'reason': rsn, 'state': st,
                     'lo_poly': _match_poly(lo), 'hi_poly': _match_poly(hi),
                     'lo_graze': _match_graze(lo), 'hi_graze': _match_graze(hi)}
                    for (lo, hi, mg, ho, dep, rsn, st) in leaves]
        if leaf_polys:
            return [(lo, hi, mg, ho, dep, _match_poly(lo), _match_poly(hi))
                    for (lo, hi, mg, ho, dep, _rsn, _st) in leaves]
        return [lf[:5] for lf in leaves]
    if not leaves:
        return None
    # Orphan = widest leaf with an EXACT head-on (ring-zero) interior wall-copy.
    headon = [lf for lf in leaves if lf[3]]
    if not headon:
        return None
    olo, ohi, omg = max(headon, key=lambda x: x[1] - x[0])[:3]
    paired = [lf[2] for lf in leaves if not lf[3]]
    runner = min(paired) if paired else 9.9

    def match(val):
        best, bd = None, 1e-9
        for b, poly in boundaries.items():
            if abs(b - val) < bd:
                bd = abs(b - val); best = poly
        return best

    return {
        'lo': olo, 'hi': ohi, 'min_gap': omg, 'runner_gap': runner,
        'exact_headon': True, 'lo_poly': match(olo), 'hi_poly': match(ohi),
        'n_branches': len(leaves), 'P': P, 'Q': Q,
    }


def trace_headon(ctx, s, max_steps=200000, closure_tol=1e-6):
    """Trace ONE perpendicular beam at launch height s (single path, no tree).

    Because the developed ray is horizontal (Lemma D1) and a branch's word fixes
    its development (Lemma A), a launch is in the (unique) orphan branch IFF its
    ray meets an EXACTLY-vertical interior wall-copy — true for the whole branch,
    not just its center. So `headon` is the orphan-interval indicator function.

    Returns (headon, min_gap, steps, closed):
      headon   True iff an interior wall-copy was hit exactly vertical (ring-zero)
      min_gap  smallest interior-wall verticality along the orbit (float; →0 at
               the orphan) — a continuous proximity signal for bracketing
      steps    reflections until perpendicular closure (or max_steps)
      closed   whether the orbit closed (perpendicular return) before max_steps
    """
    from exact_state_ring import RingTriangleState
    cot = 1.0 / math.tan(ctx.alpha)
    st = RingTriangleState(ctx)
    st.reflect('H')
    prev = 'H'
    ho = False
    mg = 9.9
    for step in range(1, max_steps):
        if step > 3:
            dx = (st.Rx_c - st.Ox_c).to_float()
            dy = (st.Ry - st.Oy).to_float()
            dax = (st.Ax_c - st.Ox_c).to_float()
            day = (st.Ay - st.Oy).to_float()
            t = closure_tol
            if ((abs(dx - 1) < t and abs(dy) < t and abs(dax - 1) < t and abs(day - 1) < t)
                    or (abs(dx + 1) < t and abs(dy) < t and abs(dax + 1) < t and abs(day - 1) < t)):
                return ho, mg, step, True
        cand = SIDES_AFTER[prev]
        side = cand[0] if len(cand) == 1 else _hit_side_ring(st, cand, s, cot)
        if side is None:
            return ho, mg, step, False
        mg = min(mg, _vgap_ring(st, side, cot))
        if _exact_vertical(st, side):
            ho = True
        st.reflect(side)
        prev = side
    return ho, mg, max_steps, False


def locate_orphan(P, Q, s_seed, max_steps=200000, bracket=1e-2, tol=1e-12):
    """Given a launch height s_seed KNOWN inside the orphan branch, bisect both
    boundaries to a float bracket. Returns (lo, hi, steps). Single-ray only."""
    from ring_cyc import RingContext
    ctx = RingContext(P, Q)
    ho, _, steps, _ = trace_headon(ctx, s_seed, max_steps)
    if not ho:
        raise ValueError(f"s_seed={s_seed} is not inside the orphan branch")

    def inside(s):
        return trace_headon(ctx, s, max_steps)[0]

    def bisect(inside_pt, outside_pt):
        while abs(inside_pt - outside_pt) > tol:
            mid = 0.5 * (inside_pt + outside_pt)
            if inside(mid):
                inside_pt = mid
            else:
                outside_pt = mid
        return 0.5 * (inside_pt + outside_pt)

    d = bracket
    while inside(s_seed - d):
        d *= 2
    lo = bisect(s_seed, s_seed - d)
    d = bracket
    while inside(s_seed + d):
        d *= 2
    hi = bisect(s_seed, s_seed + d)
    return lo, hi, steps


def orphan_width_ring(P, Q, s_lo, s_hi, dps=60, edge_tol=1e-9, **kw):
    """Exact orphan width via the ring backend (mpmath) + diagnostics.

    An orphan may abut the interval edge (lo→0 or hi→1); those endpoints have
    exact value 0 / 1 and are not recorded as split-point boundaries, so we
    supply them directly."""
    from mpmath import mp, mpf
    r = find_orphan_ring(P, Q, s_lo, s_hi, **kw)
    if r is None:
        return None, r
    mp.dps = dps
    lo_v = mpf(0) if r['lo'] < edge_tol else (
        r['lo_poly'].to_mpf(dps) if r['lo_poly'] is not None else None)
    hi_v = mpf(1) if r['hi'] > 1 - edge_tol else (
        r['hi_poly'].to_mpf(dps) if r['hi_poly'] is not None else None)
    if lo_v is None or hi_v is None:
        return None, r
    return hi_v - lo_v, r


if __name__ == '__main__':
    # CosPoly backend (shallow reference)
    for (P, Q, lo, hi) in [(13, 21, 0.255, 0.275), (21, 34, 0.138, 0.166)]:
        w, r = orphan_width_hp(P, Q, lo, hi)
        print(f"[cospoly] {P}/{Q}: orphan=[{r['lo']:.6f},{r['hi']:.6f}] "
              f"width={float(w):.8e} min_gap={r['min_gap']:.1e} "
              f"runner_gap={r['runner_gap']:.1e}")
    # Ring backend (must agree) — cross-check
    for (P, Q, lo, hi) in [(13, 21, 0.255, 0.275), (21, 34, 0.138, 0.166)]:
        w, r = orphan_width_ring(P, Q, lo, hi)
        print(f"[ring]    {P}/{Q}: orphan=[{r['lo']:.6f},{r['hi']:.6f}] "
              f"width={float(w):.8e} exact_headon={r['exact_headon']} "
              f"runner_gap={r['runner_gap']:.1e}")
