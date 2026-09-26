"""A high-precision floating-point billiard tracer, used ONLY to propose words.

Every statement in this deposit is decided by exact arithmetic (`triangle.py`,
`closure.py`).  The tracer's one job is to suggest the combinatorial input -- which
sequence of sides an orbit meets -- that the exact code then verifies; a wrong proposal
is rejected there, never silently accepted.  It also measures the deformed-angle flow of
Section 3.1, where the angle is not a rational multiple of pi and the comparison with an
exact formula is the point.

The triangle is placed as in `triangle.py`: O = (0,0), R = (L,0), A = (L,1), L = cot(alpha).
Coordinates are never snapped to a side; a hit is accepted within 10^(12 - dps) of a
segment, which at the default 60 digits is far below any quantity this deposit prints.
Words use the letters '1' (L1), '2' (L2), 'h' (H).
"""
from mpmath import mp, mpf, sqrt, cos, sin

mp.dps = 60


class Billiard:
    def __init__(self, alpha):
        self.alpha = alpha
        self.L = 1 / mp.tan(alpha)
        n = sqrt(self.L ** 2 + 1)
        self.ex, self.ey = self.L / n, 1 / n              # unit vector along H, O -> A

    def step(self, px, py, vx, vy):
        """(t, side) of the next hit from (px, py) in direction (vx, vy), or None."""
        L, tiny = self.L, mpf(10) ** (-mp.dps + 12)
        best = None
        if vx > 0:                                        # L1: x = L
            t = (L - px) / vx
            y = py + t * vy
            if t > tiny and -tiny <= y <= 1 + tiny:
                best = (t, '1')
        if vy < 0:                                        # L2: y = 0
            t = -py / vy
            x = px + t * vx
            if t > tiny and -tiny <= x <= L + tiny and (best is None or t < best[0]):
                best = (t, '2')
        den = vx - L * vy                                 # H: x - L y = 0
        if den != 0:
            t = -(px - L * py) / den
            if t > tiny:
                x, y = px + t * vx, py + t * vy
                u = (x * self.ex + y * self.ey) / sqrt(L ** 2 + 1)
                if -tiny <= u <= 1 + tiny and (best is None or t < best[0]):
                    best = (t, 'h')
        return best

    def reflect(self, side, vx, vy):
        if side == '1':
            return -vx, vy
        if side == '2':
            return vx, -vy
        d = vx * self.ex + vy * self.ey
        return 2 * d * self.ex - vx, 2 * d * self.ey - vy

    def orbit(self, px, py, theta, n):
        """The first n hits from (px, py) in direction theta: list of (side, x, y)."""
        vx, vy = cos(theta), sin(theta)
        out = []
        for _ in range(n):
            st = self.step(px, py, vx, vy)
            if st is None:
                break
            t, side = st
            px, py = px + t * vx, py + t * vy
            out.append((side, px, py))
            vx, vy = self.reflect(side, vx, vy)
        return out

    def first_return(self, s, theta, cap=400):
        """Launch from (L, s) in direction theta; stop at the first return to L1 in
        direction theta.  Returns (image s, word) or None if no return within `cap`."""
        px, py = self.L, mpf(s)
        vx, vy = cos(theta), sin(theta)
        tx, ty = vx, vy
        tol = mpf(10) ** (-mp.dps // 3)
        word = []
        for _ in range(cap):
            st = self.step(px, py, vx, vy)
            if st is None:
                return None
            t, side = st
            px, py = px + t * vx, py + t * vy
            word.append(side)
            vx, vy = self.reflect(side, vx, vy)
            if side == '1' and abs(vx - tx) < tol and abs(vy - ty) < tol:
                return py, ''.join(word)
        return None

    def head_on_word(self, s, theta, cap=400):
        """Launch from (L, s) in direction theta; the word up to the first perpendicular
        hit (after which the orbit retraces itself), or None."""
        px, py = self.L, mpf(s)
        vx, vy = cos(theta), sin(theta)
        tol = mpf(10) ** (-mp.dps // 3)
        word = []
        for _ in range(cap):
            st = self.step(px, py, vx, vy)
            if st is None:
                return None
            t, side = st
            px, py = px + t * vx, py + t * vy
            word.append(side)
            wx, wy = self.reflect(side, vx, vy)
            if abs(wx + vx) < tol and abs(wy + vy) < tol:
                return ''.join(word)
            vx, vy = wx, wy
        return None


def scan_words(billiard, theta, n, kind='return', cap=400):
    """Propose the words seen from a uniform scan of n - 1 launch points on L1:
    {word: first s at which it was seen}.  `kind` is 'return' or 'head_on'."""
    seen = {}
    for i in range(1, n):
        s = mpf(i) / n
        if kind == 'return':
            r = billiard.first_return(s, theta, cap)
            w = r[1] if r else None
        else:
            w = billiard.head_on_word(s, theta, cap)
        if w is not None and w not in seen:
            seen[w] = s
    return seen
