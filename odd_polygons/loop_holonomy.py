"""[Deposit] Helper, imported by the scripts of this directory. Paper: "Perpendicular billiards in rational polygons of odd denominator". Research name: probes/s624_loop_holonomy.py, copied verbatim except imports (deposit names), paths (data/ and figs/ under odd_polygons/) and this paragraph. Run through run_all.py; the usage lines below are the research repo's. Bracketed IDs and section marks (§sNNN) refer to the project's internal ledgers.

s624 -- s623-a: the HOLONOMY of the loops of each known minimal component, read on a transversal
inside ONE edge -- is every known minimal component a translation cover of a flat torus?

PRE-REGISTERED.  engine/flow_decomp.py's transversal is every entry edge laid end to end, so a
component's final interval can straddle a junction; there its IET translations carry segment-offset
terms and are NOT transverse periods of loops (s624_transverse_rank.py's r(M) is exposed to this).
Here: the induction is re-run with (i) the component's absolute start tracked through every Zorich
step, Dehn twist, merge and split (left induction removes a LEFT piece of the interval, of length
old - new of the winner; a split's suffix starts after the prefix), and (ii) each label's planar
GLUING SUM J = sum of g t_s over the copy crossings of its return path (composed exactly as the
return time h).  After Boshernitzan certifies a no-periodic component, left induction CONTINUES on it
until its interval lies inside one entry edge (the right end is fixed, lengths -> 0).  Then each
label's return loop -- flow x -> T(x), back along that edge -- is a closed curve IN M with holonomy
J_l, and its transverse period is -X(J_l).
CHECK (must hold exactly, every label): t_l = -X(J_l), t_l the IET translation.  A tracking error in
start, J or the segment test fails it.
MEASURED: rank_Q{t_l} (transverse rank of the loops) and rank_Q{J_l} (the holonomy group of the
loops, as vectors of Z[zeta_4d]).  rank J = 2 <=> every loop's holonomy lies in a rank-2 lattice
Z v1 + Z v2, i.e. M develops onto the flat torus C / (Z v1 + Z v2) (a translation cover of it).
COULD COME OUT OTHERWISE: the check fails (the instrument); rank t > 2 on some of the ten (then
s624_transverse_rank's r = 2 was a junction artefact); rank J > 2 (the component is not a torus
cover, and the 'rotation' is only a first-return fact).  CONTROL: an off-grid T0 direction whose
s624 inline run gave a no-periodic component of translation rank 4 must give rank t >= 3 here.
PRIOR ART: as s624_transverse_rank.py; ⊕ rulings.py --grep "torus cover" / "holonomy" -> none on the
ten.  DATA: -> none (J / start were never tracked; data/s623_rotation_degree.json has final perms).
ADDED IN-SESSION (after the first A run, before any write-up): AREA CHECK -- area(M) from the flow
boxes vs covol |Im(J1 conj J2)| for n = 2 (area_check); equal iff M is a flat torus minus a
zero-area, flow-parallel hole.
RESULT (s624; data/s624_loop_holonomy.json (A), data/s624_loop_holonomy_C.json (C)): every
component of the ten lies in one edge at certification; t = -X(J) exactly on every label; t-rank 2
on all; J-rank 2 (the seven n = 2) and 3 (T0, (7,7,16), (4,4,7)); area = covol to 30 digits on
28/28 n = 2 components; T0's two areas sum to 2 mu(M) s^2 with the Letter's mu(M) (s = 2 sin 8pi/30)
to 20 digits.  CONTROL C fires: off-grid T0 components with t-rank 4, J-rank 4 and 8, check exact.
Reading (topology, via rank H_1 = 2g + b + p - 1 and 'return loops generate H_1'): -> [NCYL-631].
"""
import _setup  # noqa: F401  (deposit: this directory on the path, cwd = odd_polygons/)
import json
import sys
import time

import sympy as sp

import flow_decomp as fd

LIST = [(7, 8, 15), (7, 7, 16), (4, 4, 7), (13, 20, 23), (15, 30, 47), (11, 45, 50),
        (23, 29, 58), (19, 38, 55), (11, 22, 103), (30, 41, 67)]


class HJ:
    """(return time, gluing sum) -- composes like the int h, so Decomposer's code runs unchanged."""
    __slots__ = ('h', 'J')

    def __init__(self, h, J):
        self.h, self.J = h, J

    def __add__(self, o):
        return HJ(self.h + o.h, tuple(a + b for a, b in zip(self.J, o.J)))

    def __rmul__(self, k):
        return HJ(k * self.h, tuple(k * a for a in self.J))


class PD(fd.Decomposer):
    def __init__(self, F, lengths, J0):
        super().__init__(F, lengths)
        self.h = {l: HJ(1, J0[l]) for l in lengths}
        self.start = None                                  # start of the component being worked

    def _zorich(self, top, bot):
        t = top[0]
        old = self.lam[t]
        r = super()._zorich(top, bot)
        self.start = self.F.add(self.start, self.F.sub(old, self.lam[t]))
        return r

    def _merge_or_none(self, c):
        F = self.F
        t0, b0 = c.top[0], c.bot[0]
        if t0 != b0 and F.cmp(self.lam[t0], self.lam[b0]) == 0:
            self.start = F.add(self.start, self.lam[t0])
            it = c.bot.index(t0)
            c.bot[it] = b0
            self.h[b0] = self.h[b0] + self.h[t0]
            del c.bot[0]
            del c.top[0]
            return True
        return False

    def run_tracked(self, top, bot, cap, inside_one):
        F = self.F
        work = [(fd.Comp(top, bot), F.zero)]
        done = []
        while work:
            c, self.start = work.pop()
            while True:
                rest = self._reduce(c)
                if rest is not None:
                    off = F.zero
                    for l in c.top:
                        off = F.add(off, self.lam[l])
                    work.append((rest, F.add(self.start, off)))
                    c.trk = c.saf0 = None
                if len(c.top) == 1:
                    c.kind = 'cylinder'
                    break
                if self._merge_or_none(c):
                    c.trk = c.saf0 = None
                    continue
                if c.steps >= cap:
                    c.kind = 'no_periodic_straddling' if c.kind == 'no_periodic' else 'undetermined'
                    break
                if c.kind != 'no_periodic' and self.autosimilar(c):     # s635 port (SAF-0 components)
                    c.kind, c.cert = 'no_periodic', 'autosimilar'
                if c.kind != 'no_periodic' and c.steps % fd.BOSH_EVERY == 0 and c.steps > 0:
                    y = self.boshernitzan(c)
                    if y is not None:
                        c.kind, c.cert = 'no_periodic', y
                if c.kind == 'no_periodic' and inside_one(self.start, self._len(c)):
                    break
                c.steps += 1
                if self._zorich(c.top, c.bot):
                    continue
                if self._zorich(c.bot, c.top):
                    continue
            done.append((c, self.start))
        return done

    def _len(self, c):
        s = self.F.zero
        for l in c.top:
            s = self.F.add(s, self.lam[l])
        return s


def rank(vecs):
    return int(sp.Matrix([list(v) for v in vecs]).rank()) if vecs else 0


def area_check(T, om, e, labels, Dc, tr, dps=300):
    """AREA of M from its flow boxes on the one-edge transversal (edge vector e), against the
    covolume |Im(J1 conj J2)| of the loop lattice (n = 2).  Each label's box has transverse width
    lam/2 (X is 2 cross) and constant along-flow height r = Re(conj(om)(J + t e / X(e)))/|om|:
    the loop x -> T(x) -> back along e has holonomy J, the edge part is t e / X(e).  A one-holed
    torus whose boundary is flow-parallel with zero net holonomy has area = covol exactly (the hole
    is bounded by parallel segments, so it has zero area).  High-precision evaluation (dps 300)."""
    from mpmath import mp, exp, pi, conj, re, im, fabs, sqrt
    if e is None:
        return None, None
    N = T.N
    with mp.workdps(dps):
        val = lambda m: sum(c * exp(2j * pi * k / N) for k, c in m.items() if c)
        vec = lambda v: {k: a for k, a in enumerate(v) if a}
        w = val(om)
        ez = val(e)
        Xe = T.F.mpval(T.red(T.cross2(om, e)), dps)
        area = 0
        for l in labels:
            J = val(vec(Dc.h[l].J))
            t = T.F.mpval(tr[l], dps)
            lam = T.F.mpval(Dc.lam[l], dps)
            r = re(conj(w) * (J + t * ez / Xe)) / fabs(w)
            area += lam / 2 * r / fabs(w)
        covol = None
        if len(labels) == 2:
            z1, z2 = (val(vec(Dc.h[l].J)) for l in labels)
            covol = fabs(im(z1 * conj(z2)))
        return mp.nstr(fabs(area), 30), (mp.nstr(covol, 30) if covol is not None else None)


def analyse(abc, om, cap=20000, extra=200000):
    T = fd.Triangle(abc)
    F = T.F
    pieces, lengths, top, bot, meta = fd.build_iet(T, om)
    J0 = {p: T.red(T.g_apply(pieces[p][3][0], T.t[pieces[p][4][1]])) for p in range(len(pieces))}
    # entry segments (by domain key), exact [lo, hi)
    segs = {}
    for dom, ln, img, key, key2 in pieces:
        end = F.add(dom, ln)
        if key not in segs:
            segs[key] = [dom, end]
        else:
            if F.cmp(dom, segs[key][0]) < 0:
                segs[key][0] = dom
            if F.cmp(end, segs[key][1]) > 0:
                segs[key][1] = end
    seglist = list(segs.values())
    segkeys = list(segs.keys())

    def edge_of(a, ln):
        b = F.add(a, ln)
        for key, (lo, hi) in zip(segkeys, seglist):
            if F.cmp(lo, a) <= 0 and F.cmp(b, hi) <= 0:
                g, s = key
                i, k2 = T.sides[s]
                return fd.madd(T.g_apply(g, T.V[k2]), fd.mneg(T.g_apply(g, T.V[i])))
        return None

    def inside_one(a, ln):
        b = F.add(a, ln)
        return any(F.cmp(lo, a) <= 0 and F.cmp(b, hi) <= 0 for lo, hi in seglist)

    om_mono = {om % T.N: 1} if isinstance(om, int) else om
    Dc = PD(F, lengths, J0)
    comps = Dc.run_tracked(top, bot, cap + extra, inside_one)
    out = []
    for c, st in comps:
        if c.kind != 'no_periodic':
            continue
        tr = Dc.translations(c)
        Js = [Dc.h[l].J for l in c.top]
        # X of a power-basis vector: convert to monomials {k: coeff}
        ok = True
        for l in c.top:
            Jm = {k: a for k, a in enumerate(Dc.h[l].J) if a}
            xJ = T.red(T.cross2(om_mono, Jm))
            if tuple(-a for a in xJ) != tuple(tr[l]):
                ok = False
        area, covol = area_check(T, om_mono, edge_of(st, Dc._len(c)), c.top, Dc, tr)
        out.append({'area_M': area, 'covol_J12': covol, 'n': len(c.top), 'steps': c.steps, 'inside_one_edge': inside_one(st, Dc._len(c)),
                    'check_t_eq_minus_XJ': ok, 't_rank': rank([tr[l] for l in c.top]),
                    'J_rank': rank(Js), 'h': [Dc.h[l].h for l in c.top],
                    'J': [list(j) for j in Js], 'lam': [list(Dc.lam[l]) for l in c.top]})
    return out


def main():
    which = sys.argv[1:] or ['A', 'C']
    res = {}
    if 'C' in which:
        # control: an off-grid T0 direction (vertex difference g V_i - V_k) with a rank-4 component
        T = fd.Triangle((7, 8, 15))
        rows = []
        for g in [(k, e) for e in (0, 1) for k in range(T.d)]:
            v = fd.madd(T.g_apply(g, T.V[1]), fd.mneg(T.V[2]))
            if not v or any(not any(T.red(T.cross2({j: 1}, v))) for j in range(T.N // 2)):
                continue
            t0 = time.time()
            o = analyse((7, 8, 15), v, cap=300, extra=3000)
            rows.append({'g': g, 'om': {int(k): int(c) for k, c in v.items()}, 'bad': o})
            print('C', g, [(r['n'], r['inside_one_edge'], r['check_t_eq_minus_XJ'], r['t_rank'], r['J_rank'])
                           for r in o], f'{time.time() - t0:.0f}s', flush=True)
            if any(r['t_rank'] >= 3 for r in o) and len(rows) >= 3:
                break
        json.dump(rows, open('data/s624_loop_holonomy_C.json', 'w'), indent=1, default=str)
    if 'A' in which:
        rows = []
        for abc in LIST:
            t0 = time.time()
            for name, js in (('perp', fd.perp_dirs(abc)), ('par', fd.side_dirs(abc))):
                hit = False
                for s, j in enumerate(js):
                    o = analyse(abc, j)
                    if o:
                        rows.append({'abc': abc, 'dir': f'{name}{s}', 'bad': o})
                        print('A', abc, f'{name}{s}', [(r['n'], r['inside_one_edge'], r['check_t_eq_minus_XJ'],
                                                        r['t_rank'], r['J_rank']) for r in o],
                              f'{time.time() - t0:.0f}s', flush=True)
                        hit = True
                        break
                if hit:
                    break
        res['A'] = rows
    if res:
        json.dump(res, open('data/s624_loop_holonomy.json', 'w'), indent=1, default=str)


if __name__ == '__main__':
    main()
