"""The perpendicular partition of L1, exactly, at a fixed centre P/Q.

THE OBJECT (the orphan paper, Section 2.1).  In the triangle O = (0,0), R = (cot a, 0),
A = (cot a, 1), a = P pi / 2Q, the beam leaves (cot a, s), s in (0, 1), with velocity (-1, 0).
The WORD of s is the sequence of sides struck until the first perpendicular return to L1; a
BRANCH is a maximal open interval of s on which the word is constant, and n(P/Q) is their
number.  A trajectory that meets a side head-on retraces itself, so the word is determined
by its prefix up to the FIRST head-on hit on any side: if that side is L1 the prefix is the
whole word, and otherwise the word is the prefix followed by its own reversal.  This module
computes the prefixes.

THE DEVELOPMENT.  Unfold: reflect the triangle across each side the beam crosses, so the
beam becomes the horizontal line Im z = y moving in the -x direction.  Scale the triangle
by 2 sin a, so that with zeta = exp(i pi / 2Q)

    O = 0,   R = zeta^P + zeta^-P,   A = 2 zeta^P,          y = 2 s sin a in (0, Im A).

Every developed copy of the triangle is the image of this one under an isometry
z -> t + zeta^r z or z -> t + zeta^r conj(z), and reflection in a side of direction
zeta^k is z -> X + zeta^(2k) conj(z - X).  So every developed vertex is an INTEGER
combination of powers of zeta, and one reflection is an index shift, an index reversal
and two additions on its coefficient vector: no division, no rounding, nothing
approximated.

ONE STEP.  The beam enters a copy through a side [U, V] with Im U < y < Im V.  Let W be
the third vertex.  The beam leaves through [U, W] if y < Im W and through [V, W] if
y > Im W; y = Im W is a VERTEX GRAZE, which is a branch boundary.  So a whole interval of
heights is carried along at once, and split at Im W only when Im W falls strictly inside
it.  The copy's exit side is vertical exactly when the beam meets it head-on, which ends
the prefix.  Every branch boundary is therefore Im of some developed vertex, and every
decision is the sign of Im of an element of Z[zeta_4Q].

EXACTNESS.  The sign of Im z is decided on the element z - conj(z) (purely imaginary,
reduced modulo Phi_4Q): zero exactly when that vector is zero, and otherwise by a float
when the float's own error bound is far below the gap, and by a certified high-precision
evaluation when it is not -- the same gate as `closure.py`.  No decision is made on an
unchecked float, and there is no width floor: a branch of any width is found, because
nothing is sampled.

THE CERTIFICATE.  Proposition 3.5 of the orphan paper bounds n(P/Q) <= 2Q + 1, and every
branch ends at a head-on hit, so the recursion terminates; `partition` walks until every
interval has closed, with a step cap only as a guard.  A result with `open == 0` is the
complete partition: the leaves tile (0, 1) by construction, each with its exact
endpoints.
"""
import math

import numpy as np
from sympy import Poly, cyclotomic_poly, symbols

EPS = 2.2e-16
GATE_KAPPA = 1024.0
HP_PREC0 = 113
COEF_GUARD = 1 << 40        # far below int64; a sum of N such products still cannot overflow

LETTER = {frozenset('RA'): '1', frozenset('OR'): '2', frozenset('OA'): 'h'}


class Ring:
    """Z[zeta_N], N = 4Q.  Elements are carried REDUNDANTLY as integer vectors of length N
    (coefficients of zeta^0 .. zeta^(N-1)), so that multiplying by zeta^j is a rotation of
    the vector and complex conjugation is k -> -k.  `canon` reduces to the basis of
    Z[x]/Phi_N, where equality is decided."""

    def __init__(self, Q):
        self.Q = Q
        N = self.N = 4 * Q
        x = symbols('x')
        coeffs = [int(c) for c in Poly(cyclotomic_poly(N, x), x).all_coeffs()][::-1]
        d = self.d = len(coeffs) - 1
        red = np.zeros((N, d), dtype=np.int64)
        cur = np.zeros(d, dtype=np.int64)
        cur[0] = 1
        tail = -np.array(coeffs[:d], dtype=np.int64)
        for k in range(N):
            red[k] = cur
            lead = cur[d - 1]
            cur = np.roll(cur, 1)
            cur[0] = 0
            if lead:
                cur = cur + lead * tail
        self.red = red
        self.sin_red = np.array([math.sin(k * math.pi / (2 * Q)) for k in range(N)])
        self.neg = (-np.arange(N)) % N
        self._hp = {}
        self._zc = {}

    def zeta(self, k):
        v = np.zeros(self.N, dtype=np.int64)
        v[k % self.N] = 1
        return v

    def conj(self, v):
        return v[self.neg]

    def mulz(self, v, j):
        return np.roll(v, j % self.N)

    def canon(self, v):
        return v @ self.red

    def zconj(self, v, j):
        """zeta^j conj(v), in one gather: its zeta^m coefficient is v[(j - m) mod N]."""
        idx = self._zc.get(j)
        if idx is None:
            idx = self._zc[j] = (j - np.arange(self.N)) % self.N
        return v[idx]

    def reduce(self, v):
        """v reduced modulo Phi_N and carried back in the redundant form.  Applied to every
        new vertex: without it nothing ever cancels (zeta^(k+2Q) = -zeta^k is never used),
        the coefficients double at each reflection and int64 overflows within a few
        hundred steps.  Reduced, they grow linearly in the depth (the paper's Remark 2.4)."""
        out = np.zeros(self.N, dtype=np.int64)
        out[:self.d] = v @ self.red
        if int(np.abs(out).max()) > COEF_GUARD:
            raise OverflowError('coefficient guard exceeded -- the int64 representation '
                                'is no longer safely exact')
        return out

    def im_float(self, v):
        return float(v @ self.sin_red)

    def _hp_sin(self, bits):
        tab = self._hp.get(bits)
        if tab is None:
            from mpmath import mp
            with mp.workprec(bits + 64):
                scale = mp.mpf(2) ** bits
                tab = [int(mp.nint(mp.sin(k * mp.pi / (2 * self.Q)) * scale))
                       for k in range(self.d)]
            self._hp[bits] = tab
        return tab

    def im_sign_exact(self, v):
        """sign(Im v), proved.  u = canon(v - conj v) represents 2i Im v; Im is injective on
        the purely imaginary elements, so u == 0 iff Im v == 0, and otherwise the sum
        below is nonzero and the precision loop ends."""
        u = self.canon(v - self.conj(v))
        mx = int(np.abs(u).max())
        if mx == 0:
            return 0
        nz = [int(k) for k in np.nonzero(u)[0]]
        cl = [int(u[k]) for k in nz]
        bits = HP_PREC0
        while bits <= 1 << 16:
            tab = self._hp_sin(bits)
            s = sum(c * tab[k] for k, c in zip(nz, cl))
            if abs(s) > (self.d * mx) // 2 + 1:
                return 1 if s > 0 else -1
            bits *= 2
        raise RuntimeError('im_sign_exact did not converge')


class Counters:
    def __init__(self):
        self.steps = 0
        self.float_decisions = 0
        self.exact_decisions = 0
        self.min_gate_ratio = float('inf')


def im_cmp(ring, a, b, cnt):
    """sign(Im a - Im b), exact."""
    v = a - b
    f = ring.im_float(v)
    gate = GATE_KAPPA * EPS * ring.N * max(1, int(np.abs(v).max()))
    if abs(f) > gate:
        cnt.float_decisions += 1
        r = abs(f) / gate
        if r < cnt.min_gate_ratio:
            cnt.min_gate_ratio = r
        return 1 if f > 0 else -1
    cnt.exact_decisions += 1
    return ring.im_sign_exact(v)


def _v(x):
    return None if x is None else [int(c) for c in x]


def _src(x):
    return list(x)


def _leaf_out(L):
    return {'lo': _v(L['lo']), 'hi': _v(L['hi']), 'lo_src': _src(L['lo_src']),
            'hi_src': _src(L['hi_src']), 'word': L['word'], 'head_on': L['head_on'],
            'returns': L['returns'], 'o_copy': _v(L['o_copy']), 't_copy': _v(L['t_copy'])}


def _arr(x):
    return None if x is None else np.array(x, dtype=np.int64)


def _leaf_in(d):
    return dict(d, lo=_arr(d['lo']), hi=_arr(d['hi']), o_copy=_arr(d['o_copy']),
                t_copy=_arr(d['t_copy']), lo_src=tuple(d['lo_src']),
                hi_src=tuple(d['hi_src']))


def _state_out(st):
    lo, hi, lo_src, hi_src, vs, a, b, r, flip, word = st
    return {'lo': _v(lo[0]), 'hi': _v(hi[0]), 'lo_src': _src(lo_src), 'hi_src': _src(hi_src),
            'vs': [_v(x[0]) for x in vs], 'a': a, 'b': b, 'r': int(r), 'flip': bool(flip),
            'word': ''.join(word)}


def _state_in(d, mk):
    return (mk(_arr(d['lo'])), mk(_arr(d['hi'])), tuple(d['lo_src']), tuple(d['hi_src']),
            [mk(_arr(x)) for x in d['vs']], d['a'], d['b'], d['r'], d['flip'],
            list(d['word']))


class Partition:
    # The beam can be launched perpendicular to any side.  Its placement is the base triangle
    # moved by z -> zeta^r0 z (or zeta^r0 conj z) so that the launch side is vertical with
    # the triangle on its left; the launch side's ORIGIN vertex has Im 0, and the beam runs
    # in -x.  BEAMS[side] = (r0 as a function of P, Q; flip0; origin index; other index).
    BEAMS = {'1': (lambda P, Q: 0, False, 1, 2),           # L1 = RA, already vertical
             '2': (lambda P, Q: Q, False, 0, 1),           # L2 = OR, rotated by i
             'h': (lambda P, Q: Q + P, True, 0, 2)}        # H = OA, z -> zeta^(Q+P) conj z

    def __init__(self, P, Q, beam='1'):
        assert 0 < P < Q and math.gcd(P, Q) == 1 and beam in self.BEAMS
        self.P, self.Q, self.beam = P, Q, beam
        self.ring = R = Ring(Q)
        self.base_dir = {frozenset('OR'): 0, frozenset('RA'): Q, frozenset('OA'): P}
        r0f, flip0, oi, ti = self.BEAMS[beam]
        self.r0, self.flip0 = r0f(P, Q) % R.N, flip0
        base = [np.zeros(R.N, dtype=np.int64), R.zeta(P) + R.zeta(-P), 2 * R.zeta(P)]
        place = [R.reduce(R.mulz(R.conj(v) if flip0 else v, self.r0)) for v in base]
        self.O, self.R, self.A = place
        self.origin_idx, self.other_idx = oi, ti
        self.origin, self.top = place[oi], place[ti]
        self.beam_side = oi + ti                # side id = sum of its vertex indices
        assert R.im_sign_exact(self.origin) == 0 and R.im_sign_exact(self.top) > 0

    def side_dir(self, side, r, flip):
        b = self.base_dir[side]
        return (r - b if flip else r + b) % self.ring.N

    def is_vertical(self, k):
        return k % (2 * self.Q) == self.Q

    def run(self, step_cap=10**8, resume=None):
        """Walk every interval to its head-on hit.  Returns (leaves, open, counters).

        RESUMABLE.  An interval still walking when the total step count reaches `step_cap`
        is returned in `open` with its whole walk state; pass the `state` of a capped
        `partition` back as `resume` at a higher cap and ONLY those intervals are walked
        on, from where they stopped, so an escalation costs the increment and not the sum
        of the caps.  The result is identical to a single run at the higher cap.

        Vertices are indexed 0 = O, 1 = R, 2 = A, and a side by the SUM of its two vertex
        indices (1 = OR = L2, 2 = OA = H, 3 = RA = L1).  Each vertex is carried as
        (vector, float Im, max |coefficient|), the float and the magnitude computed once
        when the vertex is created.  A word is a list extended in place and copied only
        at a split, of which there are at most 2Q + 1."""
        R, cnt, N, Q = self.ring, Counters(), self.ring.N, self.Q
        letter = {1: '2', 2: 'h', 3: '1'}
        base = {1: 0, 2: self.P, 3: Q}
        gk = GATE_KAPPA * EPS * N

        def mk(v):
            return (v, R.im_float(v), int(np.abs(v).max()))

        def cmp(a, b):
            f = a[1] - b[1]
            gate = gk * max(1, a[2] + b[2])
            if f > gate or f < -gate:
                cnt.float_decisions += 1
                rr = abs(f) / gate
                if rr < cnt.min_gate_ratio:
                    cnt.min_gate_ratio = rr
                return 1 if f > 0 else -1
            cnt.exact_decisions += 1
            return R.im_sign_exact(a[0] - b[0])

        if resume is None:
            V0 = (mk(self.O), mk(self.R), mk(self.A))
            end = ('end', None)
            stack = [(V0[self.origin_idx], V0[self.other_idx], end, end, V0,
                      self.origin_idx, self.other_idx, self.r0, self.flip0, [])]
            leaves = []
        else:
            leaves = [_leaf_in(L) for L in resume['leaves']]
            stack = [_state_in(st, mk) for st in reversed(resume['open'])]
            cnt.steps = resume['steps']
        opened = []
        while stack:
            lo, hi, lo_src, hi_src, vs, a, b, r, flip, word = stack.pop()
            while True:
                if cnt.steps >= step_cap:
                    opened.append((lo, hi, lo_src, hi_src, vs, a, b, r, flip, word))
                    break
                cnt.steps += 1
                u, v = (b, a) if cmp(vs[a], vs[b]) > 0 else (a, b)      # Im U < Im V
                w = 3 - a - b
                W = vs[w]
                if cmp(W, lo) <= 0:                                     # y > Im W throughout
                    pieces = [(lo, hi, lo_src, hi_src, v)]
                elif cmp(W, hi) >= 0:                                   # y < Im W throughout
                    pieces = [(lo, hi, lo_src, hi_src, u)]
                else:                                                   # a vertex graze
                    src = ('graze', 'ORA'[w], len(word))
                    pieces = [(lo, W, lo_src, src, u), (W, hi, src, hi_src, v)]
                nxt = []
                for i, (plo, phi, ps_lo, ps_hi, keep) in enumerate(pieces):
                    side = keep + w                     # the exit side [keep, w]
                    k = ((r - base[side]) if flip else (r + base[side])) % N
                    wd = word if i == len(pieces) - 1 else list(word)
                    wd.append(letter[side])
                    if k % (2 * Q) == Q:                                # head-on
                        leaves.append({'lo': plo[0], 'hi': phi[0], 'lo_src': ps_lo,
                                       'hi_src': ps_hi, 'word': ''.join(wd),
                                       'head_on': letter[side],
                                       'returns': side == self.beam_side,
                                       'o_copy': (vs[self.origin_idx][0]
                                                  if side == self.beam_side else None),
                                       't_copy': (vs[self.other_idx][0]
                                                  if side == self.beam_side else None)})
                        continue
                    X = vs[keep][0]
                    z = 3 - side                         # the vertex off the exit side
                    Zn = mk(R.reduce(X + R.zconj(vs[z][0] - X, (2 * k) % N)))
                    nv = list(vs)
                    nv[z] = Zn
                    nxt.append((plo, phi, ps_lo, ps_hi, nv, keep, w, (2 * k - r) % N,
                                not flip, wd))
                if not nxt:
                    break
                stack.extend(nxt[:-1])
                lo, hi, lo_src, hi_src, vs, a, b, r, flip, word = nxt[-1]
        leaves.sort(key=lambda L: R.im_float(L['lo']))
        return leaves, opened, cnt

    # -- readouts -----------------------------------------------------------------------
    def s_value(self, y):
        """s = Im y / Im A as a float (display only)."""
        return self.ring.im_float(y) / self.ring.im_float(self.top)

    def eq_im(self, a, b):
        return self.ring.im_sign_exact(a - b) == 0

    def return_interval(self, leaf):
        """J on a leaf returning to the launch side: the image interval, as a pair of
        vectors whose Im are its endpoints (the return height is measured along the
        launch-side copy from its origin vertex).  None for a retracing (J = id) leaf."""
        if not leaf['returns']:
            return None
        Rc, Ac = leaf['o_copy'], leaf['t_copy']
        # the copy runs from Rc (origin) to Ac, vertical; height above Rc along it is
        # |Im y - Im Rc|, oriented so that Ac is at height Im top
        if im_cmp(self.ring, Ac, Rc, Counters()) > 0:
            return (leaf['lo'] - Rc, leaf['hi'] - Rc)
        return (Rc - leaf['hi'], Rc - leaf['lo'])


def partition(P, Q, step_cap=10**8, beam='1', resume=None):
    """The exact perpendicular partition at P/Q, with its own consistency checks.

    Returns a dict: n (number of branches), open (intervals not closed by the cap; the
    partition is complete iff 0), orphans (branches retracing on themselves, J = id),
    J_ok (J maps the returning branches onto returning branches, endpoints equal as field
    elements, and is an involution), tiles (consecutive leaves abut exactly), words, and
    the exactness counters."""
    T = Partition(P, Q, beam)
    leaves, opened, cnt = T.run(step_cap, resume)
    R = T.ring
    tiles = bool(leaves) and T.eq_im(leaves[0]['lo'], T.origin) and T.eq_im(leaves[-1]['hi'], T.top)
    for a, b in zip(leaves, leaves[1:]):
        tiles = tiles and T.eq_im(a['hi'], b['lo'])
    orph = [i for i, L in enumerate(leaves) if not L['returns']]
    J_ok = not opened
    image = {}
    if J_ok:
        for i, L in enumerate(leaves):
            ri = T.return_interval(L)
            if ri is None:
                continue
            hits = [j for j, M in enumerate(leaves)
                    if T.eq_im(M['lo'], ri[0]) and T.eq_im(M['hi'], ri[1])]
            if len(hits) != 1:
                J_ok = False
                break
            image[i] = hits[0]
        J_ok = J_ok and all(image.get(image[i]) == i for i in image)
    return {'P': P, 'Q': Q, 'beam': beam, 'n': len(leaves), 'open': len(opened), 'orphans': orph,
            'J': image, 'J_ok': J_ok, 'tiles': tiles,
            'words': [L['word'] for L in leaves],
            'widths': [T.s_value(L['hi']) - T.s_value(L['lo']) for L in leaves],
            'steps': cnt.steps, 'exact_decisions': cnt.exact_decisions,
            'float_decisions': cnt.float_decisions,
            'safety': GATE_KAPPA * cnt.min_gate_ratio, 'leaves': leaves,
            # the walk state, JSON-able, to hand back as `resume=` at a higher cap
            'state': ({'leaves': [_leaf_out(L) for L in leaves],
                       'open': [_state_out(st) for st in opened], 'steps': cnt.steps}
                      if opened else None)}
