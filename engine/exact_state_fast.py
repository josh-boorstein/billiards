#!/usr/bin/env python3
"""
exact_state_fast.py — Fast exact triangle state using numpy-backed CosPoly.

Drop-in replacement for exact_state.ExactTriangleState using cos_poly_fast.
"""

import numpy as np
from cos_poly_fast import CosPoly, divide_by_1_plus_cos2, divide_by_1_minus_cos2


class ExactTriangleState:
    """Track vertex positions as exact polynomials in cos(2α).

    Representation: (x_c·cotα, y) where x_c and y are CosPolys.
    Reflection formulas are purely polynomial in this basis.
    """
    __slots__ = ('alpha', 'Ox_c', 'Oy', 'Rx_c', 'Ry', 'Ax_c', 'Ay')

    def __init__(self, alpha=None):
        if alpha is None:
            return
        self.alpha = alpha
        self.Ox_c = CosPoly(np.array([0], dtype=np.int64), 0)
        self.Oy = CosPoly(np.array([0], dtype=np.int64), 0)
        self.Rx_c = CosPoly(np.array([1], dtype=np.int64), 0)
        self.Ry = CosPoly(np.array([0], dtype=np.int64), 0)
        self.Ax_c = CosPoly(np.array([1], dtype=np.int64), 0)
        self.Ay = CosPoly(np.array([1], dtype=np.int64), 0)

    def copy(self):
        s = ExactTriangleState.__new__(ExactTriangleState)
        s.alpha = self.alpha
        s.Ox_c = self.Ox_c.copy()
        s.Oy = self.Oy.copy()
        s.Rx_c = self.Rx_c.copy()
        s.Ry = self.Ry.copy()
        s.Ax_c = self.Ax_c.copy()
        s.Ay = self.Ay.copy()
        return s

    def reflect(self, side):
        if side == 'H':
            self._reflect_across('R', 'O', 'A')
        elif side == 'L2':
            self._reflect_across('A', 'O', 'R')
        elif side == 'L1':
            self._reflect_across('O', 'R', 'A')

    def _reflect_across(self, point_label, line_start, line_end):
        Px_c, Py = self._get_vertex(point_label)
        Ax_c, Ay = self._get_vertex(line_start)
        Bx_c, By = self._get_vertex(line_end)

        dAx = Bx_c - Ax_c
        dAy = By - Ay
        dPx = Px_c - Ax_c
        dPy = Py - Ay

        one_plus_c2 = CosPoly(np.array([1, 1], dtype=np.int64), 0)
        one_minus_c2 = CosPoly(np.array([1, -1], dtype=np.int64), 0)

        dot_cleared = dPx * dAx * one_plus_c2 + dPy * dAy * one_minus_c2
        d2_cleared = dAx * dAx * one_plus_c2 + dAy * dAy * one_minus_c2

        new_Py_num = 2 * Ay * d2_cleared + 2 * dot_cleared * dAy - Py * d2_cleared
        new_Px_c_num = 2 * Ax_c * d2_cleared + 2 * dot_cleared * dAx - Px_c * d2_cleared

        side_pair = {line_start, line_end}
        if side_pair == {'O', 'A'}:
            new_Py_num.shift += 1
            new_Py_num._reduce()
            new_Px_c_num.shift += 1
            new_Px_c_num._reduce()
            self._set_vertex(point_label, new_Px_c_num, new_Py_num)
        elif side_pair == {'O', 'R'}:
            new_Py = divide_by_1_plus_cos2(new_Py_num)
            new_Px_c = divide_by_1_plus_cos2(new_Px_c_num)
            self._set_vertex(point_label, new_Px_c, new_Py)
        elif side_pair == {'R', 'A'}:
            new_Py = divide_by_1_minus_cos2(new_Py_num)
            new_Px_c = divide_by_1_minus_cos2(new_Px_c_num)
            self._set_vertex(point_label, new_Px_c, new_Py)

    def _get_vertex(self, label):
        if label == 'O': return self.Ox_c, self.Oy
        elif label == 'R': return self.Rx_c, self.Ry
        elif label == 'A': return self.Ax_c, self.Ay

    def _set_vertex(self, label, x_c, y):
        if label == 'O': self.Ox_c, self.Oy = x_c, y
        elif label == 'R': self.Rx_c, self.Ry = x_c, y
        elif label == 'A': self.Ax_c, self.Ay = x_c, y

    def vertex_y(self, vtype):
        if vtype == 'O': return self.Oy
        elif vtype == 'R': return self.Ry
        elif vtype == 'A': return self.Ay

    def vertex_y_float(self, vtype):
        return self.vertex_y(vtype).evaluate(self.alpha)


def test():
    """Verify against the original exact_state."""
    import math
    import sys
    sys.path.insert(0, '.')
    from word_tree import TriangleState

    for p, q in [(1, 5), (2, 7), (1, 11), (2, 11)]:
        alpha = p * math.pi / q
        if alpha >= math.pi / 2:
            continue

        words = [
            ['H', 'L2', 'H', 'L1', 'H', 'L2', 'L1'],
            ['H', 'L1', 'H', 'L2', 'H', 'L1', 'H', 'L2', 'L1', 'H'],
        ]

        for word in words:
            fs = TriangleState(alpha)
            es = ExactTriangleState(alpha)
            for s in word:
                fs.reflect(s)
                es.reflect(s)
            for v in 'ORA':
                err = abs(fs.vertex_y(v) - es.vertex_y_float(v))
                assert err < 1e-10, f"p={p} q={q} word={''.join(s[0] for s in word)} {v}: err={err}"

    print("ExactTriangleState (fast) matches float TriangleState — all tests passed!")


if __name__ == "__main__":
    test()
