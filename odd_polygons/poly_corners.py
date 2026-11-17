"""[Deposit] Helper, imported by the scripts of this directory. Paper: "Perpendicular billiards in rational polygons of odd denominator". Research name: engine/poly_corners.py, copied verbatim except imports (deposit names), paths (data/ and figs/ under odd_polygons/) and this paragraph. Run through run_all.py; the usage lines below are the research repo's. Bracketed IDs and section marks (§sNNN) refer to the project's internal ledgers.

poly_corners — FRONTIER CORNER ANGLES and Gauss–Bonnet of every component / 𝓕-component of a convex rational polygon's
flow in a direction (s945, `s944-c`), read off `poly_seps` germs placed in each vertex point's link.

Lemma GB (§s944 (7), [TRI-FREE-GB]): for a region W with completion W̄ along its frontier, corners are vertex points on ∂W̄ of
interior angle k_q·π, interior vertex points have angle 2π m_p, and Σ_p (1 − m_p) + ½ Σ_q (1 − k_q) = χ(W̄). This module
computes the k_q, the interior points and the boundary-cycle count b directly, so χ and the genus g = (2 − b − χ)/2 follow.

  links(T, j)                    → per vertex point, the copies around it in walk order (alternating s_out, s_in reflections
                                   from an e = 0 copy) with the ω-ray ('out') or −ω-ray ('in') each copy's wedge contains
  place(T, j, R)                 → per vertex point {'vertex', 'm', 'rays', 'comp'}: rays matched to `R['conn']` starts /
                                   arrivals and `R['capped']` starts; 'comp' = the component of each π-sector (sector t lies
                                   between ray t and ray t+1). Returns (pts, unmatched_in, orient_fails, reverse_fails).
  analyse(T, j, R, pts, S, cut_interior=False)
                                 → for a SET S of component indices (one component, or an 𝓕-component): corners
                                   [{'v','k','r1','r2'}], interior vertex points, χ (GB), b and genus (when not cut).
                                   cut_interior=True also cuts W open along connections with S on BOTH sides.
  free_components(R)             → (𝓕-components as sorted index lists, ι on components): free = ι(c) ≠ c, grouped
                                   through connections with free components on both sides.

CONVENTIONS: a sector right AFTER an out-ray in walk order is that ray's RIGHT side (orientation +1). `place` checks this
rule against every (out, in) neighbour pair and also counts the failures of the reversed rule, which is the control:
0 fails vs many (Π_15: 0 vs 44).
GATES (`_selftest`, ≈ 10 s): every non-stacked cylinder has b = 2, χ = 0, all k = 1; Π_N's minimal components cut along
their interior connections have χ = −1 (slit tori: genus 1, agreeing with the IET genus — rotations); Π_15's 𝓕-component
corners are (V0, 2), (V0, 2), (V2, 5).
⚠ TRAP — STACKED CYLINDERS: `flow_poly` reports two equal-period cylinders separated by a loop that has angle π on both
sides as ONE cylinder (the IET does not see the loop). `analyse` then reports χ = −1, two k = 2 corners at the loop's
vertex (Π_N, N ≡ 1 mod 8, at V1). Cut along the loop (cut_interior=True): χ = 0, every k = 1.
⚠ convex polygons only (inherited from `poly_seps`). ⚠ All vertex points are found, but the corner count assumes every
frontier ray is a saddle connection; `analyse` asserts that no capped separatrix lies on a frontier.
"""
import _setup  # noqa: F401  (deposit: this directory on the path, cwd = odd_polygons/)
import cmath
import collections
import math


def _cval(T, u):
    return sum(c * cmath.exp(2j * math.pi * e / T.N) for e, c in u.items())


def links(T, j):
    om = cmath.exp(2j * math.pi * j / T.N)
    k = T.k
    pts = []
    for i in range(k):
        s_in, s_out = (i - 1) % k, i
        seen = set()
        for g0 in [(kk, 0) for kk in range(T.d)]:
            if g0 in seen:
                continue
            seq, g, s = [], g0, s_out
            while True:
                seen.add(g)
                W = [_cval(T, T.g_apply(g, v)) for v in T.V]
                a, b = W[(i + 1) % k] - W[i], W[(i - 1) % k] - W[i]
                ab = (a.conjugate() * b).imag

                def inside(z):
                    return (a.conjugate() * z).imag * ab > 0 and (z.conjugate() * b).imag * ab > 0
                ray = 'out' if inside(om) else ('in' if inside(-om) else None)
                seq.append((g, i, ray))
                g = T.g_across(g, s)
                s = s_in if s == s_out else s_out
                if g == g0 and s == s_out:
                    break
            pts.append(seq)
    return pts


def place(T, j, R):
    conn, capped = R['conn'], R['capped']
    out_of, in_of = {}, collections.defaultdict(list)
    for idx, c in enumerate(conn):
        out_of[c['start']] = ('c', idx)
        in_of[c['arr']].append(('c', idx))
    for idx, c in enumerate(capped):
        out_of[c['start']] = ('x', idx)
    side = lambda r: (conn if r[0] == 'c' else capped)[r[1]]
    pts, unmatched = [], 0
    for seq in links(T, j):
        rays = []
        for g, v, ray in seq:
            if ray == 'out':
                rays.append(('out', out_of[(g, v)]))
            elif ray == 'in':
                lst = in_of.get((g, v), [])
                assert len(lst) <= 1
                if not lst:
                    unmatched += 1
                rays.append(('in', lst[0] if lst else None))
        pts.append({'vertex': seq[0][1], 'm': sum(r[0] == 'out' for r in rays), 'rays': rays})
    fails = {}
    for orient in (+1, -1):
        f = 0
        for r in pts:
            rays, n, comp = r['rays'], len(r['rays']), []
            for t in range(n):
                r1, r2 = rays[t], rays[(t + 1) % n]
                o, i_ = (r1, r2) if r1[0] == 'out' else (r2, r1)
                use_right = ((r1[0] == 'out') == (orient > 0))
                val = side(o[1])['right' if use_right else 'left']
                if i_[1] is not None and side(i_[1])['right' if use_right else 'left'] != val:
                    f += 1
                comp.append(val)
            if orient > 0:
                r['comp'] = comp
        fails[orient] = f
    return pts, unmatched, fails[1], fails[-1]


def analyse(T, j, R, pts, S, cut_interior=False):
    conn, capped = R['conn'], R['capped']
    S = set(S)

    def bdry(ref):
        if ref is None:
            return False
        if ref[0] == 'x':
            c = capped[ref[1]]
            assert not ((c['left'] in S) != (c['right'] in S)), 'capped separatrix on the frontier'
            return False
        c = conn[ref[1]]
        a, b = c['left'] in S, c['right'] in S
        return (a != b) or (cut_interior and a and b)
    corners, interior = [], []
    for pid, r in enumerate(pts):
        comp, rays = r['comp'], r['rays']
        n = len(rays)
        mine = [c in S for c in comp]
        if not any(mine):
            continue
        cuts = [t for t in range(n) if bdry(rays[t][1])]
        if not cuts:
            assert all(mine), 'sector run mixes S and non-S without a boundary ray'
            interior.append((pid, r['vertex'], r['m']))
            continue
        for a_i, t in enumerate(cuts):
            t2 = cuts[(a_i + 1) % len(cuts)]
            span = [(t + s) % n for s in range((t2 - t) % n or n)]
            if all(mine[x] for x in span):
                corners.append({'pid': pid, 'v': r['vertex'], 'k': len(span),
                                'r1': (rays[t][1][0], rays[t][1][1], rays[t][0]),
                                'r2': (rays[t2][1][0], rays[t2][1][1], rays[t2][0])})
            else:
                assert not any(mine[x] for x in span), 'sector run mixes S and non-S without a boundary ray'
    edges = []
    for i, c in enumerate(conn):
        a, b = c['left'] in S, c['right'] in S
        if a != b:
            edges.append(i)
        elif cut_interior and a and b:
            edges += [i, i]
    assert len(corners) == len(edges), (len(corners), len(edges))
    chi = sum(1 - m for _, _, m in interior) + sum(1 - c['k'] for c in corners) / 2
    out = {'corners': corners, 'interior': interior, 'chi': chi, 'n_edges': len(edges)}
    if not cut_interior:
        par = {}

        def f(x):
            par.setdefault(x, x)
            while par[x] != x:
                par[x] = par[par[x]]
                x = par[x]
            return x
        for i in edges:
            par[f(('c', i, 'out'))] = f(('c', i, 'in'))
        for c in corners:
            par[f(c['r1'])] = f(c['r2'])
        out['b'] = len({f(('c', i, e)) for i in edges for e in ('out', 'in')})
        out['genus'] = (2 - out['b'] - chi) / 2
    return out


def free_components(R):
    conn = R['conn']
    good = {}
    for rule in ('same', 'swap'):
        m, ok = {}, True
        for c in conn:
            k = c['partner']
            if k is None:
                continue
            a, b = (conn[k]['left'], conn[k]['right']) if rule == 'same' else (conn[k]['right'], conn[k]['left'])
            for x, y in ((c['left'], a), (c['right'], b)):
                if m.setdefault(x, y) != y:
                    ok = False
        if ok:
            good[rule] = m
    assert list(good) == ['same'], ('ι on components: expected exactly the same-side rule', list(good))
    io = good['same']
    free = {c for c in range(len(R['comps'])) if c in io and io[c] != c}
    adj = collections.defaultdict(set)
    for c in conn:
        if c['left'] in free and c['right'] in free:
            adj[c['left']].add(c['right'])
            adj[c['right']].add(c['left'])
    seen, out = set(), []
    for k in sorted(free):
        if k in seen:
            continue
        st, grp = [k], set()
        while st:
            x = st.pop()
            if x not in grp:
                grp.add(x)
                st += list(adj[x])
        seen |= grp
        out.append(sorted(grp))
    return out, io


def _selftest():
    import flow_poly as FP
    import poly_seps as PS
    from regular_pieces import asym
    for N in (15, 17, 23):
        P = asym(N)
        j = FP.perp_dirs_poly(P)[0]
        R = PS.seps(P, j)
        pts, um, fo, fr = place(P, j, R)
        assert fo == 0 and fr > 0, (fo, fr)
        for k, c in enumerate(R['comps']):
            A = analyse(P, j, R, pts, [k])
            if c['kind'] == 'cylinder':
                ok = A['b'] == 2 and A['chi'] == 0 and all(x['k'] == 1 for x in A['corners'])
                if not ok:      # stacked pair (trap): cutting along the loop must restore a cylinder pair
                    B = analyse(P, j, R, pts, [k], cut_interior=True)
                    assert B['chi'] == 0 and all(x['k'] == 1 for x in B['corners']), (N, k)
            else:
                B = analyse(P, j, R, pts, [k], cut_interior=True)
                assert B['chi'] == -1, (N, k, B['chi'])
        if N == 15:
            fc, _ = free_components(R)
            for g in fc:
                A = analyse(P, j, R, pts, g)
                assert sorted((x['v'], x['k']) for x in A['corners']) == [(0, 2), (0, 2), (2, 5)]
        print('N', N, 'orient fails 0 /', fr, 'gates OK')
    print('poly_corners selftest OK')


if __name__ == '__main__':
    _selftest()
