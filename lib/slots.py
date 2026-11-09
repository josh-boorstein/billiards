"""Slots and prongs of a direction class of ANY rational triangle, traced exactly.

THE OBJECT (the polygon note, Section 4).  For a class C of perpendicular directions, a SLOT is a vertex V with a
direction u in C strictly inside the angle at V; its PRONG is the billiard path leaving V along u, followed until it
meets a side at a right angle ('perp': a singular point, Theorem 4.2 (ii)), reaches a vertex ('vertex'), or neither
within the step cap ('cap': a prong that never reaches a vertex is counted by s').

THE INSTRUMENT.  The flow in the class direction omega on the unfolded surface (`flow.Triangle`), in the transverse
coordinate x = cross(omega, p), which the flow preserves.  Each forward slot is the middle vertex (in x) of a copy at
which omega points inward; its prong is pushed edge to edge until x equals a copy's exit-exit vertex (a vertex) or it
crosses a side copy perpendicular to omega.  Exact throughout: every comparison is in the cyclotomic field.

    classes(abc) -> [{'j', 'vert', 'hor'}]           one representative direction j per class, its vertical and
                                                    horizontal sides
    trace_class(abc, j, cap) -> (T, slots, geo)     every slot's prong
    pair_reversals(T, slots, geo) -> (pairs, closed)  the vertex-to-vertex prongs matched with their reverses

At even N the unfolding carries both orientations of each slot, so `len(slots) = 2S` and the counts of pairs and of
capped prongs are doubled; the scripts divide by `mult = 2`.
"""
from flow import Triangle, perp_dirs, madd, mneg


def trace_class(abc, j, cap=200000):
    T = Triangle(abc)
    F, d = T.F, T.d
    om = {j % T.N: 1}
    X = lambda p: T.red(T.cross2(om, p))
    copies = [(k, e) for e in (0, 1) for k in range(d)]
    geo = {}
    for g in copies:
        W = [T.g_apply(g, v) for v in T.V]
        xs = [X(w) for w in W]
        sides = {}
        for s, (i, k2) in enumerate(T.sides):
            kk = 3 - i - k2
            e = madd(W[k2], mneg(W[i]))
            c_om = F.sign(T.red(T.cross2(e, om)))
            if c_om == 0:
                sides[s] = ('par', None, None)
                continue
            c_in = F.sign(T.red(T.cross2(e, madd(W[kk], mneg(W[i])))))
            lo, hi = (xs[i], xs[k2]) if F.cmp(xs[i], xs[k2]) < 0 else (xs[k2], xs[i])
            # perpendicular side copy: e parallel to i*omega <=> cross(i*omega, e) = 0
            iom = {(j + d) % T.N: 1}
            perp = not any(T.red(T.cross2(iom, e)))
            sides[s] = ('ent' if c_om == c_in else 'ext', (lo, hi), perp)
        # middle vertex: x strictly between the other two
        mid = None
        for v in range(3):
            o = [w for w in range(3) if w != v]
            if F.cmp(xs[o[0]], xs[v]) * F.cmp(xs[v], xs[o[1]]) > 0:
                mid = v
        geo[g] = (W, xs, sides, mid)

    def side_opp(v):
        return [s for s, (i, k2) in enumerate(T.sides) if v not in (i, k2)][0]

    def step_exit(g, s2, x):
        """leaf at transverse x leaves copy g through side s2: return (g', s2, x) in g' coords."""
        g2 = T.g_across(g, s2)
        shift = X(T.g_apply(g, T.t[s2]))
        return g2, s2, F.sub(x, shift)

    def run(g, v):
        """forward prong from middle vertex v of copy g."""
        W, xs, sides, mid = geo[g]
        x = xs[v]
        s_ex = side_opp(v)
        assert sides[s_ex][0] == 'ext'
        cur = step_exit(g, s_ex, x)
        for it in range(cap):
            g, s, x = cur
            Wg, xg, sg, midg = geo[g]
            kind, (lo, hi), perp = sg[s]
            assert kind == 'ent' and F.cmp(lo, x) < 0 and F.cmp(x, hi) < 0, (kind, it)
            if perp:
                return {'end': 'perp', 'side': s, 'steps': it}
            exits = [t for t in range(3) if sg[t][0] == 'ext']
            if len(exits) == 1:
                cur = step_exit(g, exits[0], x)
                continue
            w = midg                                   # exit-exit vertex
            c = F.cmp(x, xg[w])
            if c == 0:
                return {'end': 'vertex', 'copy': g, 'vertex': w, 'steps': it}
            for t in exits:
                lo2, hi2 = sg[t][1]
                if F.cmp(lo2, x) < 0 and F.cmp(x, hi2) < 0:
                    cur = step_exit(g, t, x)
                    break
            else:
                raise AssertionError('no exit')
        return {'end': 'cap', 'steps': cap}

    slots = []
    for g in copies:
        W, xs, sides, mid = geo[g]
        if mid is None:
            continue
        s_opp = side_opp(mid)
        if sides[s_opp][0] != 'ext':
            continue                                   # omega points outward at mid: backward ray
        r = run(g, mid)
        r['start'] = (g, mid)
        slots.append(r)
    return T, slots, geo


def pair_reversals(T, slots, geo):
    """A closed non-palindromic prong ends at vertex w of copy g arriving along omega; its reverse is
    the forward prong from w in the copy reached by the reflection fixing omega's perpendicular...
    Rather than reconstruct it, match by billiard data: vertex index of the end must be the start
    vertex index of a distinct closed non-palindromic slot, with equal step count."""
    closed = [i for i, r in enumerate(slots) if r['end'] == 'vertex']
    used, pairs = set(), []
    for i in closed:
        if i in used:
            continue
        ri = slots[i]
        cand = [k for k in closed if k not in used and k != i
                and slots[k]['start'][1] == ri['vertex'] and ri['start'][1] == slots[k]['vertex']
                and slots[k]['steps'] == ri['steps']]
        if not cand:
            return None, closed
        used |= {i, cand[0]}
        pairs.append((i, cand[0]))
    return pairs, closed


def classes(abc):
    """beam direction classes: one representative j per class, its vertical sides (a copy
    perpendicular to omega) and horizontal sides (a copy parallel to omega)."""
    T = Triangle(abc)
    pj = perp_dirs(abc)
    out = []

    def having(dirmono):
        res = []
        for t, (i, k2) in enumerate(T.sides):
            for g in [(k, e) for e in (0, 1) for k in range(T.d)]:
                e_ = madd(T.g_apply(g, T.V[k2]), mneg(T.g_apply(g, T.V[i])))
                if not any(T.red(T.cross2(dirmono, e_))):
                    res.append(t)
                    break
        return res
    for s, j in enumerate(pj):
        vert = having({(j + T.d) % T.N: 1})
        hor = having({j % T.N: 1})
        if not any(set(vert) == set(c['vert']) for c in out):
            out.append({'j': j, 'vert': vert, 'hor': hor})
    return out
