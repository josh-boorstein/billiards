"""The necklace readout: the incidences of the separatrix diagram with Fix(iota).

`closure.py` walks every separatrix of the genus-zero base B and decides whether each arrives
at a singularity.  This module walks the same separatrices, by the same kite maps and the same
exact arithmetic, and records in addition where each one meets the fixed set of the reversing
involution iota.  At odd Q that fixed set is the horizontal arc(s) (walked by no separatrix:
they ARE separatrices, listed `in_fix`), the vertical arm of Fix(iota) in kite 0, and -- in the
class where it is vertical -- the middle interface m = (Q - 1)/2.  So along a walk

  * `t0`  counts the transits of kite 0, each of which crosses its vertical arm;
  * `onm` counts the crossings of interface m.

From those, per class (the necklace paper's notation):

  h   the separatrices lying in Fix(iota) (its horizontal side-copies);
  v   the vertical side-copies: 2 when interface m is vertical (the PAIRED class), else 1;
  f   the points of the critical graph on the vertical arcs, counted once per EDGE;
  D   = Q - h - f, Corollary 9.8's covering defect; D = 0 certifies complete periodicity
      (a capped walk can only undercount f, so a capped row's D is an upper bound);
  E   the number of edges of the critical graph (Theorem 3.1: E = Q on a CP class);
  multi   RESOLVED edges meeting the vertical arcs more than once -- Lemma 3.2a says none;
  swapped the iota-swapped PAIRS among the resolved edges: iota sends the prong germ in kite k to
      the one in kite -k, so an edge with germs (a, k), (b, l) has image (a, -k), (b, -l); it is
      iota-INVARIANT when that is the same edge and swapped otherwise.  (An edge joining germs
      k and -k of one kind is invariant -- the type (2) edges of Lemma 3.2a -- not swapped.)
      On a CP class this is p, so D = 2 * swapped (Theorem 3.3).
At even Q, iota's bookkeeping is not derived; the readout gives the CP certificate, E, and
Theorem 3.1's cylinder count C = (Q - b)/2 + kappa - 1 (b the pole-zero edges, kappa = 1 iff
some edge joins Z_O to Z_A, else 2).
"""
from closure import Chain


class NecklaceChain(Chain):
    def trace_fix(self, j, x, k):
        """`closure.Chain.trace`, also counting kite-0 transits and interface-m crossings.
        Returns (end, steps, t0, onm)."""
        steps, t0 = 0, 0
        onm = 1 if (not self.even and j == self.m_idx) else 0
        imv = self.ring.imvals
        while True:
            steps += 1
            if steps > self.step_cap:
                return 'CAP', steps, t0, onm
            if not steps & 511:
                x.f = float(x.v @ imv)
                a = int(abs(x.v).max())
                self._gate = self.gate(a)
                if a > self.maxcoef:
                    self.maxcoef = a
            r = self.cross(k, j, x)
            if r[0] == 'pole':
                return ('R', r[1]), steps, t0, onm
            j2, x2, transited = r
            if transited and k == 0:
                t0 += 1
            if -self._gate < x2.f < self._gate and x2.is_zero():
                return ('ZO', k), steps, t0, onm
            if self.eq_exact(x2, self.ell[j2]):
                return ('ZA', k), steps, t0, onm
            if not self.even and j2 == self.m_idx:
                onm += 1
            k = self.next_kite(k, j2)
            j, x = j2, x2

    def separatrices_fix(self):
        """`closure.Chain.separatrices`, each walked separatrix carrying t0 and onm."""
        seps = []
        # the start list is built exactly as closure.Chain.separatrices builds it
        Q = self.Q
        if self.even:
            for k in range(Q):
                if self.tie[k]:
                    seps.append({'ends': [('R', k), ('ZO' if k == 0 else 'ZA', k)],
                                 'in_fix': True, 'steps': 0, 't0': 0, 'onm': 0})
            for t in range(Q):
                if self.ell[t].is_zero():
                    seps.append({'ends': [('ZO', None), ('ZA', None)], 'in_fix': True,
                                 'steps': 0, 't0': 0, 'onm': 0})
        else:
            seps.append({'ends': [('R', 0), ('ZO' if self.eps == 0 else 'ZA', 0)],
                         'in_fix': True, 'steps': 0, 't0': 0, 'onm': 0})
            if self.ell[self.m_idx].is_zero():
                seps.append({'ends': [('ZO', None), ('ZA', None)], 'in_fix': True,
                             'steps': 0, 't0': 0, 'onm': 0})
        for k in range(Q):
            if self.tie[k]:
                continue
            far = self.far[k]
            starts = [(far, self.ffix[k], k, ('R', k))]
            if not self.ell[self.near[k]].is_zero():
                starts.append((far, self.fhi[k] if self.osec[k] else self.flo[k], k,
                               ('ZO' if self.osec[k] else 'ZA', k)))
            for j, x, kk, tag in starts:
                end, steps, t0, onm = self.trace_fix(j, type(x)(x.v.copy(), x.f),
                                                     self.next_kite(kk, j))
                seps.append({'ends': [tag, end], 'in_fix': False, 'steps': steps,
                             't0': t0, 'onm': onm})
        return seps


def readout(P, Q, eps, step_cap=3_000_000, with_seps=False):
    """The necklace readout of one class; see the module docstring for the fields."""
    ch = NecklaceChain(P, Q, eps, step_cap)
    seps = ch.separatrices_fix()
    traced = [s for s in seps if not s['in_fix']]
    fixed = [s for s in seps if s['in_fix']]
    capped = sum(1 for s in traced if s['ends'][1] == 'CAP')
    out = {'P': P, 'Q': Q, 'eps': eps, 'cap': step_cap, 'capped': capped,
           'cp': capped == 0, 'n_traced': len(traced), 'h': len(fixed),
           'E': len(fixed) + len(traced) // 2,
           'max_steps': max((s['steps'] for s in traced if s['ends'][1] != 'CAP'), default=0),
           'sign_calls': ch.sign_calls}
    resolved = [s for s in traced if s['ends'][1] != 'CAP']
    out['swapped'] = None
    if not ch.even:                          # iota's action on kites is odd-Q bookkeeping
        edges = {frozenset((tuple(s['ends'][0]), tuple(s['ends'][1]))) for s in resolved}
        image = lambda e: frozenset((a, (-k) % Q) for a, k in e)
        out['swapped'] = sum(1 for e in edges if image(e) != e) // 2
    if not ch.even:
        m_vertical = not ch.ell[ch.m_idx].is_zero()
        out['v'] = 2 if m_vertical else 1
        out['lone'] = not m_vertical
        f1 = sum(s['t0'] > 0 for s in traced) // 2
        fm = (sum(s['onm'] > 0 for s in traced) // 2) if m_vertical else 0
        out['f'] = f1 + fm
        out['n_sigma'] = {('L1' if eps == 0 else 'L2'): f1 + 1}
        if m_vertical:
            out['n_sigma']['H'] = fm + 1
        out['D'] = Q - out['h'] - out['f']
        out['multi'] = sum(1 for s in resolved
                           if s['t0'] + (s['onm'] if m_vertical else 0) > 1)
    else:
        out['C'] = None
        if capped == 0:
            ends = [(s['ends'][0][0], s['ends'][1][0]) for s in traced]
            zero = ('ZO', 'ZA')
            rz = sum((a == 'R') != (b == 'R') for a, b in ends) // 2
            rz += sum(1 for s in fixed if s['ends'][0][0] == 'R')
            oa = sum({a, b} == {'ZO', 'ZA'} for a, b in ends) // 2
            oa += sum(1 for s in fixed if tuple(e[0] for e in s['ends']) == zero)
            kappa = 1 if oa else 2
            if (Q - rz) % 2 == 0:
                out['C'] = (Q - rz) // 2 + kappa - 1
            out['b'], out['kappa'] = rz, kappa
    if with_seps:
        out['seps'] = seps
    return out


# -- the interval model (the necklace paper, Section 4) ---------------------------------------

def interval_model(P, Q, eps):
    """The global transverse coordinate of Theorem 4.1: interface t is the interval
    I_t = [u_t, u_t + ell_t], with u_0 = 0 and each kite fixing its far interface's position
    from its near one's (through band: equal, or shifted by the fold width in the O-sector).
    Returns (chain, lo, hi, closes): `closes` is the consistency of kite 0 with the other
    Q - 1 kites, decided exactly -- the statement that the coordinate is global."""
    from closure import Chain
    ch = Chain(P, Q, eps)
    u = {0: ch.zero}
    for k in range(1, Q):
        far, near, d = ch.far[k], ch.near[k], ch.delta[k]
        if far in u:
            u[near] = u[far] + d if ch.osec[k] else u[far]
        else:
            u[far] = u[near] - d if ch.osec[k] else u[near]
    far, near, d = ch.far[0], ch.near[0], ch.delta[0]
    closes = (u[far] + d if ch.osec[0] else u[far]).eq(u[near])
    lo = {t: u[t] for t in range(Q)}
    hi = {t: u[t] + ch.ell[t] for t in range(Q)}
    return ch, lo, hi, closes


def pole_y(ch, lo, hi, k):
    """The prong of the pole R_k in the global coordinate: the midpoint of the tail
    I_far minus I_near, the direction away from kite k, and the far interface."""
    from closure import half
    far, near = ch.far[k], ch.near[k]
    bot = lo[far].eq(lo[near])
    y = half((hi[near] + hi[far]) if bot else (lo[far] + lo[near]))
    return y, (-1 if far == (k - 1) % ch.Q else +1), far


def first_run_escapes(ch, lo, hi, k):
    """Does R_k's prong reach a Fix(iota) arc (interface m, or a transit of kite 0) before
    its first fold?"""
    Q, m = ch.Q, ch.m_idx
    y, d, t = pole_y(ch, lo, hi, k)
    crossed, steps = (t == m), 0
    while steps <= 2 * Q:
        steps += 1
        t2 = (t + d) % Q
        if not (ch.cmp(y, lo[t2]) > 0 and ch.cmp(y, hi[t2]) < 0):
            return crossed
        if (d == 1 and t == Q - 1) or (d == -1 and t == 0):
            crossed = True
        t = t2
        if t == m:
            crossed = True
    return crossed


def quotient_trace(ch, lo, hi, k, cap=200_000):
    """Follow R_k's prong in the global coordinate, where the leaf's height is constant
    between folds and a fold at kite k' sends y to (the fold's tail sum) - y.  Returns
    (end, folds, crossed_fix, steps): end is ('R'|'ZO'|'ZA', kite) or ('CAP', None), and
    crossed_fix whether the leaf reached interface m or transited kite 0 (a Fix(iota) arc).
    Every comparison is `closure.Chain.cmp`, exact when close; y's float is refreshed from
    its exact vector every 512 steps (it is rebuilt by s - y at each fold, so it drifts)."""
    from closure import half
    Q, m = ch.Q, ch.m_idx
    imv = ch.ring.imvals
    y, d, t = pole_y(ch, lo, hi, k)
    folds, crossed, steps = 0, (t == m), 0
    while True:
        steps += 1
        if steps > cap:
            return ('CAP', None), folds, crossed, steps
        if not steps & 511:
            y.f = float(y.v @ imv)
            a = int(abs(y.v).max())
            ch._gate = ch.gate(a)
        t2 = (t + d) % Q
        kk = t2 if d == 1 else (t2 + 1) % Q
        cl, chg = ch.cmp(y, lo[t2]), ch.cmp(y, hi[t2])
        if cl == 0:
            return ('ZO', kk), folds, crossed, steps
        if chg == 0:
            return ('ZA', kk), folds, crossed, steps
        if cl > 0 and chg < 0:                                   # through
            if (d == 1 and t == Q - 1) or (d == -1 and t == 0):
                crossed = True
            t = t2
            if t == m:
                crossed = True
            continue
        s = (hi[t2] + hi[t]) if lo[t].eq(lo[t2]) else (lo[t] + lo[t2])
        if y.eq(half(s)):
            return ('R', kk), folds, crossed, steps
        y = s - y
        d = -d
        folds += 1
