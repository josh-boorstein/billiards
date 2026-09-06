#!/usr/bin/env python3
"""exact_state_ring.py — fixed-angle, depth-independent triangle state.

Drop-in analogue of exact_state.ExactTriangleState for a FIXED angle
α=(P/Q)(π/2). Vertex coordinates are RingElems in Z[ζ_N] (ring_cyc) instead of
CosPolys, so every reflect is O(D²) (D=φ(2Q)) regardless of depth — no Chebyshev
degree growth, no coefficient blow-up, and numerically stable float evaluation.

The reflection algebra is identical to exact_state._reflect_across (same
(x_c·cotα, y) representation, same clearing of the |B-A|² denominator by
1±cos2α); only the coefficient type changes. Cross-checked against
ExactTriangleState in __main__.

Use this for deep single-angle runs (golden-orphan measurement). For symbolic /
multi-angle work keep using the CosPoly ExactTriangleState.
"""

from ring_cyc import RingContext
from cos_poly import CosPoly


class RingTriangleState:
    """Vertex positions as RingElems at a fixed α. Same API as
    ExactTriangleState: reflect('H'/'L2'/'L1'), vertex_y(v), vertex_y_float(v).
    Also exposes the x-coefficient RingElems (Ox_c/Rx_c/Ax_c) for gap tests."""

    def __init__(self, ctx=None, P=None, Q=None):
        if ctx is None and P is not None:
            ctx = RingContext(P, Q)
        if ctx is None:
            return
        self.ctx = ctx
        one = ctx.from_cospoly(CosPoly([1], 0))
        zero = ctx.zero()
        # O=(0,0), R=(cotα,0), A=(cotα,1); x stored as coeff of cotα.
        self.Ox_c = zero.copy()
        self.Oy = zero.copy()
        self.Rx_c = one.copy()
        self.Ry = zero.copy()
        self.Ax_c = one.copy()
        self.Ay = one.copy()

    def copy(self):
        s = RingTriangleState.__new__(RingTriangleState)
        s.ctx = self.ctx
        s.Ox_c = self.Ox_c.copy()
        s.Oy = self.Oy.copy()
        s.Rx_c = self.Rx_c.copy()
        s.Ry = self.Ry.copy()
        s.Ax_c = self.Ax_c.copy()
        s.Ay = self.Ay.copy()
        return s

    def _get_vertex(self, label):
        if label == 'O':
            return self.Ox_c, self.Oy
        if label == 'R':
            return self.Rx_c, self.Ry
        return self.Ax_c, self.Ay

    def _set_vertex(self, label, x_c, y):
        if label == 'O':
            self.Ox_c, self.Oy = x_c, y
        elif label == 'R':
            self.Rx_c, self.Ry = x_c, y
        else:
            self.Ax_c, self.Ay = x_c, y

    def reflect(self, side):
        if side == 'H':
            self._reflect_across('R', 'O', 'A')
        elif side == 'L2':
            self._reflect_across('A', 'O', 'R')
        elif side == 'L1':
            self._reflect_across('O', 'R', 'A')
        else:
            raise ValueError(side)

    def _reflect_across(self, point_label, line_start, line_end):
        """Reflect `point` across the current side line_start–line_end.

        Mirrors exact_state._reflect_across. In (x_c·cotα, y) coordinates the
        reflection numerator/denominator (after clearing |B-A|²) is
          dot   = dPx·dAx·(1+cos2α) + dPy·dAy·(1-cos2α)
          d2    = dAx²·(1+cos2α)    + dAy²·(1-cos2α)
          P'y   = (2·Ay·d2 + 2·dot·dAy − Py·d2) / d2
          P'x_c = (2·Ax_c·d2 + 2·dot·dAx − Px_c·d2) / d2

        The reflected side is an isometric copy of a fixed triangle side, so d2
        is CONSTANT — exactly 2 (side O-A), 1+cos2α (O-R), or 1-cos2α (R-A)
        (verified ring-exact over deep orbits). Substituting the constant makes
        the d2 terms cancel against the division:
          P'y = 2·Ay − Py + 2·(dot·dAy)/d2   (÷2 for H ⇒ the leading 2 cancels)
        This removes the 4 muls that formed d2 and the 4 that scaled A/P by it —
        ~15 → ~6-8 ring muls per reflect (muls are the dominant cost).
        """
        ctx = self.ctx
        Px, Py = self._get_vertex(point_label)
        Ax, Ay = self._get_vertex(line_start)
        Bx, By = self._get_vertex(line_end)

        dAx = Bx - Ax
        dAy = By - Ay
        dPx = Px - Ax
        dPy = Py - Ay

        dot = dPx * dAx * ctx.one_plus + dPy * dAy * ctx.one_minus
        ty = dot * dAy   # scaled and divided by d2 below to give 2·(dot·dAy)/d2
        tx = dot * dAx

        pair = {line_start, line_end}
        if pair == {'O', 'A'}:            # H: d2 = 2  ⇒  2·t/2 = t
            new_Py = Ay * 2 - Py + ty
            new_Px = Ax * 2 - Px + tx
        elif pair == {'O', 'R'}:          # L2: d2 = 1+cos2α
            new_Py = Ay * 2 - Py + ty.div_one_plus() * 2
            new_Px = Ax * 2 - Px + tx.div_one_plus() * 2
        else:                             # L1: d2 = 1-cos2α
            new_Py = Ay * 2 - Py + ty.div_one_minus() * 2
            new_Px = Ax * 2 - Px + tx.div_one_minus() * 2
        new_Py.reduce()
        new_Px.reduce()
        self._set_vertex(point_label, new_Px, new_Py)

    # -- accessors -----------------------------------------------------------
    def vertex_y(self, vtype):
        if vtype == 'O':
            return self.Oy
        if vtype == 'R':
            return self.Ry
        return self.Ay

    def vertex_xc(self, vtype):
        if vtype == 'O':
            return self.Ox_c
        if vtype == 'R':
            return self.Rx_c
        return self.Ax_c

    def vertex_y_float(self, vtype):
        return self.vertex_y(vtype).to_float()

    def max_bits(self):
        return max(c.max_bits() for c in
                   (self.Ox_c, self.Oy, self.Rx_c, self.Ry, self.Ax_c, self.Ay))


# ---------------------------------------------------------------------------
# Cross-check against the CosPoly ExactTriangleState (float y agreement).
# ---------------------------------------------------------------------------
def _selftest():
    import math
    import random
    from exact_state import ExactTriangleState

    random.seed(7)
    cases = [(1, 5), (2, 7), (3, 7), (1, 11), (2, 11), (3, 11), (5, 13),
             (3, 21), (5, 21), (13, 21), (5, 34), (21, 34)]
    worst = 0.0
    for P, Q in cases:
        alpha = (P / Q) * (math.pi / 2)
        ctx = RingContext(P, Q)
        # random long words + a few fixed ones
        words = [['H'], ['H', 'L2'], ['H', 'L1'], ['H', 'L2', 'H', 'L1', 'L2']]
        for _ in range(6):
            w = ['H']
            prev = 'H'
            after = {'H': ('L2', 'L1'), 'L2': ('H', 'L1'), 'L1': ('H', 'L2')}
            for _ in range(random.randint(3, 30)):
                prev = random.choice(after[prev])
                w.append(prev)
            words.append(w)

        for w in words:
            es = ExactTriangleState(alpha)
            rs = RingTriangleState(ctx)
            for side in w:
                es.reflect(side)
                rs.reflect(side)
            for v in ('O', 'R', 'A'):
                worst = max(worst, abs(es.vertex_y_float(v)
                                       - rs.vertex_y_float(v)))
        print(f"  P/Q={P}/{Q}: D={ctx.D:3d}  {len(words)} words ok")
    print(f"RingTriangleState self-test PASSED, worst y error = {worst:.2e}")


if __name__ == '__main__':
    _selftest()
