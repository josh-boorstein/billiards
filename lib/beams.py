"""Every branch of a perpendicular beam of ANY rational triangle, exactly classified.

THE OBJECT (the polygon note, Sections 2-3).  The triangle has angles pi*(a_0, a_1, a_2)/N at V_0, V_1, V_2; side 0
is V0V1, side 1 is V1V2, side 2 is V2V0 (`flow.Triangle`'s labelling).  The perpendicular beam of a side sigma is cut
by its singular points into branches, and each branch is classified by the FIRST perpendicular (head-on) collision of
its paths: on a side tau != sigma the branch is RETRACING ('retr'); on sigma itself it is paired by the return map J,
and the slope of J is eps = (orientation of the launch copy of sigma) x (orientation of the return copy): eps = -1 a
MIRROR pair, eps = +1 a TRANSLATION pair (Theorem 3.1).

THE INSTRUMENT.  `flow.build_iet` is the exact one-crossing interval exchange of the perpendicular direction on the
unfolded surface; a position induction (`flow.PosDecomposer`) gives every branch's exact start and width on every
entry edge of sigma, gated by an exact tiling of EVERY such edge; each branch's exact midpoint (doubled
coordinates, no division) is then pushed through the exchange until it enters an edge perpendicular to the flow.
Counts and kinds are exact; the returned positions along the side are 40-digit evaluations, for display.

    classify(abc, side) -> {'ok', 'n', 'R', 'M', 'T', 'E', 'retr_by_side', 'branches': [(t_lo, t_w, kind, end_side, eps)]}

`R`, `M`, `T` count BRANCHES of the side (a retracing band has one branch on each of two sides, a pair two branches
on this one); `t` is measured from the side's first vertex in units of its length.
"""
import bisect
import functools
import math
from math import gcd

from flow import Triangle, build_iet, PosDecomposer, _HK, perp_dirs, madd, mneg, mshift


def gprime(a):
    """The genus g' of the quotient S/tau of an even-N triangle (Remark 5.3): 2 - 2g' = -N/2 + sum gcd(a_i, N/2)."""
    N = sum(a)
    if N % 2:
        return None
    h = N // 2
    chi = -h + sum(gcd(x % h if x % h else h, h) for x in a)
    return (2 - chi) // 2


def _perp_keys(T, om, side=None):
    """Entry-capable keys (g, s) whose edge is perpendicular to omega (s == side if given)."""
    out = []
    for g in [(k, e) for e in (0, 1) for k in range(T.d)]:
        W = [T.g_apply(g, v) for v in T.V]
        for s, (v1, v2) in enumerate(T.sides):
            if side is not None and s != side:
                continue
            e_ = madd(W[v2], mneg(W[v1]))
            if not any(T.red(T.cross2(mshift(om, T.d, T.N), e_))):
                out.append((g, s))
    return out


def _kinds(comps):
    k = {}
    for c in comps:
        k[c.kind] = k.get(c.kind, 0) + 1
    return k


def classify(abc, side, cap=20000, maxsteps=10**9):
    """Every branch of side `side`'s perpendicular beam, exactly classified. Returns a dict with
    'branches': [(t_lo, t_w, kind, end_side, eps)] sorted along the side from its first vertex
    (t in units of the side length), and the counts."""
    T = Triangle(abc)
    F = T.F
    j = perp_dirs(abc)[side]
    om = {j % T.N: 1}
    X = lambda p: T.red(T.cross2(om, p))
    keys = _perp_keys(T, om, side)
    pieces, lengths, top, bot, meta = build_iet(T, j)
    idx = {key: i for i, key in enumerate(keys)}
    origin = {}
    for p in lengths:
        kk = pieces[p][3]
        if kk not in origin or F.cmp(pieces[p][0], origin[kk]) < 0:
            origin[kk] = pieces[p][0]
    pos0 = {p: [(idx[pieces[p][3]], F.sub(pieces[p][0], origin[pieces[p][3]]))]
            for p in lengths if pieces[p][3] in idx}
    Dc = PosDecomposer(F, lengths, pos0)
    zero = (0,) * len(keys)
    Dc.h = {}
    for p in lengths:
        kv = list(zero)
        if pieces[p][3] in idx:
            kv[idx[pieces[p][3]]] = 1
        Dc.h[p] = _HK(1, tuple(kv))
    comps = Dc.run(top, bot, cap)
    total = F.zero
    for c in comps:
        for l in c.top:
            total = F.add(total, F.smul(Dc.h[l].h, Dc.lam[l]))
    assert F.cmp(total, meta['L']) == 0, 'Kac mass balance failed'
    edge = {key: F.zero for key in keys}
    for p in lengths:
        if pieces[p][3] in idx:
            edge[pieces[p][3]] = F.add(edge[pieces[p][3]], lengths[p])
    cyl = [c.top[0] for c in comps if c.kind == 'cylinder']
    per_key = []
    for i, key in enumerate(keys):
        if F.sign(edge[key]) == 0:
            continue
        iv = [(p, Dc.lam[l], l) for l in cyl for kk, p in Dc.pos[l] if kk == i]
        cnt = sum(Dc.h[l].k[i] for l in cyl)
        iv.sort(key=functools.cmp_to_key(lambda a, b: F.cmp(a[0], b[0])))
        cur = F.zero
        for p, w, _ in iv:
            if F.cmp(p, cur) != 0:
                return {'ok': False, 'why': 'gap', 'abc': abc, 'side': side}
            cur = F.add(cur, w)
        if F.cmp(cur, edge[key]) != 0 or len(iv) != cnt:
            return {'ok': False, 'why': 'tiling', 'abc': abc, 'side': side}
        per_key.append((i, iv))
    if not per_key:
        return {'ok': False, 'why': 'no entry edge', 'abc': abc, 'side': side}
    ns = {len(iv) for _, iv in per_key}
    if len(ns) != 1:
        return {'ok': False, 'why': 'edges disagree', 'abc': abc, 'side': side}

    # ---- exact midpoint trace, doubled coordinates ----
    allperp = set(_perp_keys(T, om))
    two = lambda u: F.smul(2, u)
    starts2 = []
    order = sorted(range(len(pieces)), key=functools.cmp_to_key(
        lambda p, q: F.cmp(pieces[p][0], pieces[q][0])))
    starts2 = [two(pieces[p][0]) for p in order]
    ends2 = [two(F.add(pieces[p][0], pieces[p][1])) for p in order]
    shift2 = [two(F.sub(pieces[p][2], pieces[p][0])) for p in order]
    fv = [float(F.mpval(s, 30)) for s in starts2]

    def locate(Y):
        y = float(F.mpval(Y, 30))
        k = bisect.bisect_right(fv, y) - 1
        for kk in (k, k - 1, k + 1, k - 2, k + 2):
            if 0 <= kk < len(order):
                if F.cmp(starts2[kk], Y) < 0 and F.cmp(Y, ends2[kk]) < 0:
                    return kk
                if F.cmp(starts2[kk], Y) == 0:
                    raise RuntimeError('midpoint hit a piece boundary (vertex)')
        raise RuntimeError('locate failed')

    def orient(g):
        v1, v2 = T.sides[side]
        W = [T.g_apply(g, v) for v in T.V]
        return 1 if F.cmp(X(W[v1]), X(W[v2])) < 0 else -1, X(W[v1]), X(W[v2])

    i0, iv0 = per_key[0]
    key0 = keys[i0]
    o0, x1, x2 = orient(key0[0])
    E0 = F.mpval(edge[key0], 40)
    res = []
    for p, w, l in iv0:
        Y = F.add(two(F.add(origin[key0], p)), w)        # 2 * (midpoint), concatenated coordinate
        steps = 0
        while True:
            kk = locate(Y)
            Y = F.add(Y, shift2[kk])
            kto = pieces[order[kk]][4]
            steps += 1
            if kto in allperp:
                break
            if steps > min(maxsteps, Dc.h[l].h + 2):       # one period always returns to the launch edge
                raise RuntimeError('trace cap')
        tau = kto[1]
        if tau != side:
            kind, eps = 'retr', 0
        else:
            o1, _, _ = orient(kto[0])
            eps = o0 * o1
            kind = 'mirror' if eps == -1 else 'transl'
        tl = float(F.mpval(p, 40) / E0)
        tw = float(F.mpval(w, 40) / E0)
        if o0 == -1:                                     # measure t from the side's first vertex
            tl = 1.0 - tl - tw
        res.append((tl, tw, kind, tau, eps, steps, Dc.h[l].h))
    res.sort()
    cnt = {'retr': 0, 'mirror': 0, 'transl': 0}
    rby = {}
    for r in res:
        cnt[r[2]] += 1
        if r[2] == 'retr':
            rby[r[3]] = rby.get(r[3], 0) + 1
    E = cnt['retr'] + cnt['transl'] / 2
    return {'ok': True, 'abc': abc, 'side': side, 'n': len(res), 'R': cnt['retr'], 'M': cnt['mirror'],
            'T': cnt['transl'], 'E': E, 'retr_by_side': rby, 'nkeys': len(per_key),
            'ncomp': len(comps), 'ncyl': len(cyl), 'kinds': _kinds(comps),
            'branches': [(round(a, 12), round(b, 12), k, t, e) for a, b, k, t, e, _, _ in res]}
