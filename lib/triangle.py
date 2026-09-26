"""The right triangle with centre P/Q, and the exact unfolding of an orbit along a fixed word.

THE TRIANGLE.  alpha = P pi / 2Q is the acute angle at O = (0,0); the right angle is at
R = (L, 0) and the third vertex is A = (L, 1), with L = cot(alpha).  The sides are

    L1 = RA : x = L,  0 < y < 1       (vertical leg, opposite alpha)
    L2 = OR : y = 0,  0 < x < L       (horizontal leg)
    H  = OA : x = L y, 0 < y < 1      (hypotenuse)

This is the placement of the paper.  A word is a string over {'1', '2', 'h'}.

DIRECTIONS are integers m modulo 4Q, standing for the angle m pi / 2Q.  Reflection in a
side acts on the angle as phi -> 2(side angle) - phi, so in these units

    L1: m -> 2Q - m,     L2: m -> -m,     H: m -> 2P - m.

WHY AN ORBIT ALONG A FIXED WORD IS A FINITE EXACT COMPUTATION.  The reflections act on
directions by integer arithmetic, so along a fixed word the direction at every step is a
known constant, independent of where the orbit started.  If the starting point depends
affinely on a parameter s, so does every later point:  p_j = A_j + B_j s,  with A_j, B_j
exact elements of Q(zeta_4Q).  Each hitting time t_j = t0_j + t1_j s is affine too.
Consequently
  * a condition "the orbit arrives exactly at a vertex" is one linear equation in s;
  * the set of s for which the word IS the orbit (every t_j > 0, every hit strictly inside
    its side) is an intersection of half-lines: a single interval with exact endpoints;
  * the path length, sum of t_j (the direction vectors are unit vectors), is affine in s.
Nothing is subdivided and nothing is approximated; comparisons are exact
(`cyclofield.Elem.sign`).
"""
from cyclofield import CycloField

SIDES = ('1', '2', 'h')


class RightTriangle:
    def __init__(self, P, Q):
        from math import gcd
        assert 0 < P < Q and gcd(P, Q) == 1
        self.P, self.Q = P, Q
        self.F = CycloField(4 * Q)
        F = self.F
        self.L = F.cos(P) / F.sin(P)                      # cot(alpha), alpha = P units
        self.M = 4 * Q                                    # directions live mod 4Q

    # -- directions ---------------------------------------------------------------------
    def reflect(self, side, m):
        Q, P = self.Q, self.P
        return {'1': 2 * Q - m, '2': -m, 'h': 2 * P - m}[side] % self.M

    def dirvec(self, m):
        return self.F.cos(m), self.F.sin(m)

    def is_head_on(self, side, m):
        """Does direction m meet `side` perpendicularly (so that the orbit retraces)?"""
        return self.reflect(side, m) == (m + 2 * self.Q) % self.M

    # -- the unfolding ------------------------------------------------------------------
    def launch_L1(self):
        """The starting point (L, s) on L1, as (Ax, Ay, Bx, By)."""
        F = self.F
        return (self.L, F.zero, F.zero, F.one)

    def launch_point(self, x, y):
        F = self.F
        return (x, y, F.zero, F.zero)

    def unfold(self, word, m0, start=None):
        """Force the orbit through the sides of `word`, starting in direction m0.

        Returns a list of steps, one per letter:
            (side, m_in, t0, t1, (Ax, Ay), (Bx, By), m_out)
        where the hit point is (Ax + Bx s, Ay + By s), reached after time t0 + t1 s in
        direction m_in, and m_out is the direction after reflection.
        """
        Ax, Ay, Bx, By = start if start is not None else self.launch_L1()
        L = self.L
        m = m0 % self.M
        steps = []
        for S in word:
            dx, dy = self.dirvec(m)
            if S == '1':                                  # x = L
                t0, t1 = (L - Ax) / dx, -Bx / dx
            elif S == '2':                                # y = 0
                t0, t1 = -Ay / dy, -By / dy
            else:                                         # x - L y = 0
                den = dx - L * dy
                t0, t1 = -(Ax - L * Ay) / den, -(Bx - L * By) / den
            Ax, Ay = Ax + t0 * dx, Ay + t0 * dy
            Bx, By = Bx + t1 * dx, By + t1 * dy
            m2 = self.reflect(S, m)
            steps.append((S, m, t0, t1, (Ax, Ay), (Bx, By), m2))
            m = m2
        return steps

    # -- the constraints that make a forced word the actual orbit ------------------------
    def constraints(self, step):
        """The affine constraints f0 + f1 s > 0 at one step, as (name, vertex, f0, f1).

        `vertex` names the corner the orbit reaches when the constraint becomes tight
        (None for the time constraint)."""
        F, L = self.F, self.L
        S, m, t0, t1, (Ax, Ay), (Bx, By), _ = step
        out = [('t>0', None, t0, t1)]
        if S == '1':
            out += [('y>0', 'R', Ay, By), ('y<1', 'A', F.one - Ay, -By)]
        elif S == '2':
            out += [('x>0', 'O', Ax, Bx), ('x<L', 'R', L - Ax, -Bx)]
        else:
            out += [('y>0', 'O', Ay, By), ('y<1', 'A', F.one - Ay, -By)]
        return out

    def validity_interval(self, word, m0, lo=None, hi=None, start=None):
        """The exact interval (lo, hi) of s on which `word` is the orbit.

        Returns a dict with lo, hi, the binding constraint at each end (step index,
        side, name, vertex) and the number of constraints, or None if the set is empty.
        """
        F = self.F
        lo = F.zero if lo is None else lo
        hi = F.one if hi is None else hi
        lo_src = hi_src = None
        n = 0
        steps = self.unfold(word, m0, start)
        for j, st in enumerate(steps, 1):
            for name, vtx, f0, f1 in self.constraints(st):
                n += 1
                sg = f1.sign()
                if sg == 0:
                    if f0.sign() <= 0:
                        return None
                    continue
                b = -f0 / f1
                if sg > 0 and b > lo:
                    lo, lo_src = b, (j, st[0], name, vtx)
                elif sg < 0 and b < hi:
                    hi, hi_src = b, (j, st[0], name, vtx)
        if not lo < hi:
            return None
        return {'lo': lo, 'hi': hi, 'lo_src': lo_src, 'hi_src': hi_src,
                'n_constraints': n, 'steps': steps}

    # -- combinatorial conditions on a word (direction arithmetic only, hence exact) -----
    def first_return_ok(self, word, m0, side='1'):
        """Is `word` a FIRST return to (side, outgoing m0)?  True iff its last letter is
        the first at which the orbit reflects off `side` into direction m0."""
        m = m0 % self.M
        for j, S in enumerate(word, 1):
            m2 = self.reflect(S, m)
            if S == side and m2 == m0 % self.M:
                return j == len(word)
            m = m2
        return False

    def head_on_ok(self, word, m0):
        """Does `word` end at its FIRST perpendicular hit?"""
        m = m0 % self.M
        for j, S in enumerate(word, 1):
            if self.is_head_on(S, m):
                return j == len(word)
            m = self.reflect(S, m)
        return False

    @staticmethod
    def path_length(steps):
        """(T0, T1): the path length along the steps is T0 + T1 s."""
        T0 = sum((st[2] for st in steps[1:]), steps[0][2])
        T1 = sum((st[3] for st in steps[1:]), steps[0][3])
        return T0, T1

    def tile(self, words, m0, lo=None, hi=None):
        """Validity intervals of `words`, sorted, and whether they tile [lo, hi] exactly.

        Returns (rows, tiles, defects).  Consecutive endpoints are compared as field
        elements, so "abut" means equal, not close."""
        F = self.F
        lo = F.zero if lo is None else lo
        hi = F.one if hi is None else hi
        rows = []
        for w in words:
            v = self.validity_interval(w, m0)
            if v is None:
                return None, False, [('empty', w)]
            v['word'] = w
            rows.append(v)
        rows.sort(key=lambda r: r['lo'].value(30))
        defects = []
        if not rows[0]['lo'] == lo:
            defects.append(('start', rows[0]['word']))
        for a, b in zip(rows, rows[1:]):
            if not a['hi'] == b['lo']:
                defects.append(('gap' if a['hi'] < b['lo'] else 'overlap',
                                a['word'], b['word']))
        if not rows[-1]['hi'] == hi:
            defects.append(('end', rows[-1]['word']))
        return rows, not defects, defects
