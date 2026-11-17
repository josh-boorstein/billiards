"""[Deposit] Helper, imported by the scripts of this directory. Paper: "Perpendicular billiards in rational polygons of odd denominator". Research name: engine/poly_seps.py, copied verbatim except imports (deposit names), paths (data/ and figs/ under odd_polygons/) and this paragraph. Run through run_all.py; the usage lines below are the research repo's. Bracketed IDs and section marks (§sNNN) refer to the project's internal ledgers.

poly_seps — EXACT saddle connections of a convex rational polygon in a direction, with the component on EACH SIDE (s940).

`s624_slit_pairs.trace_all` traces every forward ω-separatrix of a TRIANGLE's unfolding exactly through `flow_decomp`'s
one-crossing IET; `flow_poly` decomposes a convex polygon (`PolyT`) but traces nothing. This module ports the tracer to
`PolyT` and adds what Theorem P's proof reads (§s927 (3), §s928): for every saddle connection, which component lies on its
LEFT and on its RIGHT (read along ω), from two tracked germs at its start through `flow_poly.MarkDecomposer`; the ι-pairing
(ι = the reflection whose axis is ⊥ ω, on germs (g, v) ↦ (r g, v), r = (k_r, 1), k_r = (j + d)/2 — §s927's `swapped_pairs`);
and, per component, the frontier sides and the exact HOLONOMY BALANCE Σ_{left sides} hol − Σ_{right sides} hol = 0.

  seps(T, j, cap)          → {'conn': [...], 'capped': [...], 'comps': [...], 'pairs': [...]}:
                             conn[i] = {start (g, v), arr (g, v), hol (tuple, ℤ[ζ] reduced), length, steps, left, right (component
                             indices), partner (index of ι·γ; = i for type (2)), pair (index into 'pairs' or None)};
                             capped = separatrices not ending by `cap` (start, left, right); comps = flow_poly components;
                             pairs[k] = (i, i′) the swapped pairs (type (3)), i < i′.
  balance(R, c)            the frontier of component c: [(conn index, side ±1)] (+1 = c on its LEFT), and the exact balance
                             (asserted zero) ⊕ α(c) = per swapped pair L − R (Theorem P's relation, §s935 (3) Lemma S).
CONVENTIONS: x = 2 Im(conj(ω)·P) on the DEVELOPED copy, so a larger x is LEFT of ω; a germ (pos, +1) sits at x just above
its position — left of a leaf through it. A separatrix starts at a lower-chain vertex of a copy (two entry sides meet) and
ends when its orbit lands exactly on a vertex position ('lo' end of an entry segment, or a 'mid' breakpoint = the vertex
between two consecutive exit sides). Sides parallel to ω are refused (none at odd N in the perpendicular class, Lemma H).
GATES (`_selftest`): on triangles, the saddle connections (start, end, holonomy) are IDENTICAL to `s624_slit_pairs.trace_all`
and the swapped pairs reproduce `s927_free_region.swapped_pairs` (holonomies, p); every component's balance is exactly zero
(asserted on every call); Kac via `flow_poly`.
"""
import _setup  # noqa: F401  (deposit: this directory on the path, cwd = odd_polygons/)
import functools

import flow_decomp as FD
import flow_poly as FP
from flow_decomp import madd, mneg, mshift


def _tup(u):
    return FD.CycNP.tup(u)


def seps(T, j, cap=30000, mark=True):
    assert T.convex and T.M == 1, 'poly_seps: convex PolyT, M = 1'
    pieces, lengths, top, bot, meta = FP._build_iet_poly(T, j)
    F = FD._fast(T)
    N, d = T.N, T.d
    om = {e % N: c for e, c in j.items()} if isinstance(j, dict) else {j % N: 1}
    iom = mshift(om, d, N)
    X = lambda p: F.red(T.cross2(om, p))
    n = len(pieces)
    off = meta['off']
    copies = [(k, e) for e in (0, 1) for k in range(d)]
    # refuse sides parallel to ω
    for g in copies:
        W = [T.g_apply(g, v) for v in T.V]
        for s, (i, k2) in enumerate(T.sides):
            assert any(T.red(T.cross2(madd(W[k2], mneg(W[i])), om))), 'a side is parallel to ω'
    dom = [_tup(pieces[p][0]) for p in range(n)]
    ln = [_tup(pieces[p][1]) for p in range(n)]
    shift = [_tup(F.sub(pieces[p][2], pieces[p][0])) for p in range(n)]
    Jp = [_tup(T.red(T.g_apply(pieces[p][3][0], T.t[pieces[p][4][1]]))) for p in range(n)]
    Fc = T.F                                               # tuple arithmetic
    add = lambda a, b: tuple(x + y for x, y in zip(a, b))
    # vertex positions on the transversal
    brk = {}
    seg_end = {}
    for (g, s), (O, lo, hi) in off.items():
        i, k2 = T.sides[s]
        W = [T.g_apply(g, v) for v in T.V]
        vlo, vhi = (i, k2) if F.cmp(X(W[i]), X(W[k2])) < 0 else (k2, i)
        plo, phi_ = _tup(O), _tup(F.add(O, F.sub(hi, lo)))
        brk.setdefault(plo, []).append((g, vlo, 'lo'))
        brk.setdefault(phi_, []).append((g, vhi, 'hi'))
        seg_end[(g, s)] = (plo, phi_, vlo, vhi)
    byseg = {}
    for p in range(n):
        byseg.setdefault(pieces[p][3], []).append(p)
    for key, ps in byseg.items():
        ps.sort(key=functools.cmp_to_key(lambda a, b: Fc.cmp(dom[a], dom[b])))
        for q, p in zip(ps, ps[1:]):
            sa, sb = set(T.sides[pieces[q][4][1]]), set(T.sides[pieces[p][4][1]])
            v = (sa & sb)
            assert len(v) == 1, 'consecutive exit sides share no vertex'
            brk.setdefault(dom[p], []).append((key[0], v.pop(), 'mid'))
    order = list(top)                                       # pieces sorted by domain start (exact, from the builder)
    piece_at = {dom[p]: p for p in range(n)}
    # separatrix starts: lower-chain vertices (two entry sides of one copy meet)
    starts = []
    entry_by_g = {}
    for (g, s) in off:
        entry_by_g.setdefault(g, []).append(s)
    for g, ss in entry_by_g.items():
        for a in ss:
            for b in ss:
                if a == b:
                    continue
                pa, pha, vloa, vhia = seg_end[(g, a)]
                pb, phb, vlob, vhib = seg_end[(g, b)]
                if vloa == vhib:                           # a's lo vertex is b's hi vertex: v, a to the left of ω's leaf
                    starts.append(((g, vloa), pa, (pa, 1), (phb, -1)))
    # s940 FAST WALK: exact positions as int64 vectors (CycNP), decisions by float with a rigorous error bound, exact only
    # on a near tie — the same decisions as the tuple walk (gate: `_selftest` vs s624, identical connections)
    import numpy as np
    Ca = F.Ca
    domv = [F.vec(dom[p]) for p in range(n)]
    shv = [F.vec(shift[p]) for p in range(n)]
    Jv = [F.vec(Jp[p]) for p in range(n)]
    dom_f = [float(np.dot(Ca, domv[p])) for p in range(n)]
    dom_m = [int(np.abs(domv[p]).sum()) for p in range(n)]
    sf = np.array([dom_f[p] for p in order])
    bkeys = [k for k, ents in brk.items() if any(e[2] in ('lo', 'mid') for e in ents)]
    bvec = [F.vec(k) for k in bkeys]
    bf = np.array([float(np.dot(Ca, v)) for v in bvec])
    bm = [int(np.abs(v).sum()) for v in bvec]
    bord = np.argsort(bf)
    bfs = bf[bord]
    eps = F.bound_c

    def is_brk(xv, xf, xm):
        k = int(np.searchsorted(bfs, xf))
        for kk in (k - 1, k, k + 1):
            if 0 <= kk < len(bfs):
                b = int(bord[kk])
                if abs(bfs[kk] - xf) <= eps * (xm + bm[b]) + 1e-300 and np.array_equal(xv, bvec[b]):
                    return bkeys[b]
        return None

    def locate_f(xv, xf, xm):
        k = int(np.searchsorted(sf, xf, side='right')) - 1
        k = max(0, min(k, n - 1))
        while k > 0 and (xf - sf[k] < eps * (xm + dom_m[order[k]]) and F.sign(xv - domv[order[k]]) < 0
                         or xf - sf[k] < -eps * (xm + dom_m[order[k]])):
            k -= 1
        while k < n - 1 and (sf[k + 1] - xf < eps * (xm + dom_m[order[k + 1]]) and F.sign(xv - domv[order[k + 1]]) >= 0
                             or sf[k + 1] - xf < -eps * (xm + dom_m[order[k + 1]])):
            k += 1
        return order[k]
    conn, capped = [], []
    for (g, v), x0, gl, gr in starts:
        p = piece_at[x0]
        xv = F.vec(x0) + shv[p]
        Jv_ = Jv[p].copy()
        steps, hit, qprev = 1, None, p
        while steps < cap:
            xf = float(np.dot(Ca, xv))
            xm = int(np.abs(xv).sum())
            key = is_brk(xv, xf, xm)
            if key is not None:
                hit = [e for e in brk[key] if e[2] in ('lo', 'mid')]
                break
            q = locate_f(xv, xf, xm)
            xv = xv + shv[q]
            Jv_ = Jv_ + Jv[q]
            qprev = q
            steps += 1
        if not hit:
            capped.append({'start': (g, v), 'gl': gl, 'gr': gr, 'steps': steps})
            continue
        assert len(hit) == 1, hit
        J = tuple(int(a) for a in Jv_)
        g2, v2, kind = hit[0]
        arr = (g2, v2) if kind == 'mid' else (pieces[qprev][3][0], v2)
        Ws = T.red(T.g_apply(g, T.V[v]))
        We = T.red(T.g_apply(g2, T.V[v2]))
        hol = tuple(int(a - b + c) for a, b, c in zip(We, Ws, J))
        hd = {k: c for k, c in enumerate(hol) if c}
        assert not any(T.red(T.cross2(om, hd))), 'holonomy not parallel to ω'
        import cmath, math
        omabs = abs(sum(c * cmath.exp(2j * math.pi * k / N) for k, c in om.items()))
        length = -Fc.fval(T.red(T.cross2(iom, hd))) / 2 / omabs
        assert length > 0
        conn.append({'start': (g, v), 'arr': arr, 'hol': hol, 'length': length, 'steps': steps, 'gl': gl, 'gr': gr})
    # ι-pairing
    kr = (j + d) // 2 if not isinstance(j, dict) else None
    assert kr is not None and (j + d) % 2 == 0, 'ι needs j + d even (a perpendicular-class direction)'
    rg = lambda g: ((kr - g[0]) % d, 1 - g[1])
    bystart = {c['start']: i for i, c in enumerate(conn)}
    pairs = []
    for i, c in enumerate(conn):
        k = bystart.get((rg(c['arr'][0]), c['arr'][1]))
        c['partner'] = k
        if k is None:
            continue                                       # ι-image is a capped separatrix
        assert conn[k]['hol'] == c['hol'], 'ι-partners with different holonomy'
        if k > i:
            pairs.append((i, k))
    for c in conn:
        c['pair'] = None
    for t, (i, k) in enumerate(pairs):
        conn[i]['pair'] = conn[k]['pair'] = t
    R = {'conn': conn, 'capped': capped, 'pairs': pairs, 'comps': None, 'N': N, 'j': j}
    if mark:
        germs = [c['gl'] for c in conn] + [c['gr'] for c in conn] + [c['gl'] for c in capped] + [c['gr'] for c in capped]
        germs = [(F.vec(y) if hasattr(F, 'vec') else y, s) for (y, s) in germs]
        Dc = FP.MarkDecomposer(F, lengths, germs)
        comps = Dc.run(top, bot, 5000)
        cidx = {id(c): k for k, c in enumerate(comps)}
        where = [cidx[id(p[0])] for p in Dc.pts]
        nc, ncp = len(conn), len(capped)
        for i, c in enumerate(conn):
            c['left'], c['right'] = where[i], where[nc + i]
        for i, c in enumerate(capped):
            c['left'], c['right'] = where[2 * nc + i], where[2 * nc + ncp + i]
        out, total = [], F.zero
        for c in comps:
            mass = F.zero
            for l in c.top:
                mass = F.add(mass, F.smul(Dc.h[l], Dc.lam[l]))
            total = F.add(total, mass)
            out.append({'kind': c.kind, 'size': len(c.top), 'period': Dc.h[c.top[0]] if len(c.top) == 1 else None,
                        'mass': F.fval(mass) / F.fval(meta['L'])})
        assert F.cmp(total, meta['L']) == 0, 'Kac mass balance failed'
        R['comps'] = out
        for k in range(len(out)):
            balance(R, k)
    for c in conn:
        del c['gl'], c['gr']
    return R


def balance(R, k):
    """frontier sides of component k and the exact balance; α = per swapped pair (L − R)"""
    conn = R['conn']
    fr = []
    for i, c in enumerate(conn):
        if c['left'] == k and c['right'] != k:
            fr.append((i, 1))
        elif c['right'] == k and c['left'] != k:
            fr.append((i, -1))
    D = len(conn[0]['hol']) if conn else 0
    tot = [0] * D
    for i, s in fr:
        for t, h in enumerate(conn[i]['hol']):
            tot[t] += s * h
    # s940 audit fix: excuse the balance ONLY if a capped separatrix lies ON k's frontier (left ≠ right, one side k); a capped
    # separatrix interior to k (left = right = k: a minimal component's dense separatrix) does not touch the frontier
    open_fr = any((c['left'] == k) != (c['right'] == k) for c in R['capped'])
    assert not any(tot) or open_fr, ('balance fails', k, fr)
    R.setdefault('balance_checked', {})[k] = not open_fr
    alpha = {}
    for i, s in fr:
        t = conn[i]['pair']
        key = ('pair', t) if t is not None else ('self', i)
        alpha[key] = alpha.get(key, 0) + s
    return fr, {k_: v for k_, v in alpha.items() if v}


def _selftest():
    import slit_pairs as S6
    import free_region as FR
    for abc in [(2, 5, 8), (3, 17, 37), (3, 41, 61), (5, 13, 27)]:
        d = sum(abc)
        T = FP.triangle_poly(abc)
        j = sorted(set(FD.perp_dirs(abc)))[0]
        R = seps(T, j)
        ref = S6.trace_all(abc, j, 30000)[0]
        A = sorted((c['start'], c['arr'], c['hol']) for c in R['conn'])
        B = sorted(((tuple(o['start'][0]), o['start'][1]), (tuple(o['arr'][0]), o['arr'][1]), tuple(o['hol']))
                   for o in ref if o['end'] and not o.get('edge'))
        assert A == B, ('trace mismatch', abc, len(A), len(B))
        prs = FR.swapped_pairs(abc)[0]
        assert sorted(R['conn'][i]['hol'] for i, _ in R['pairs']) == sorted(h for h, _, _ in prs), abc
        kinds = [c['kind'] for c in R['comps']]
        print(abc, 'conn', len(R['conn']), 'pairs', len(R['pairs']), 'comps', len(kinds), 'all balances exact')
    print('poly_seps selftest OK')


if __name__ == '__main__':
    _selftest()
