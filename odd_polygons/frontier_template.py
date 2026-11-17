"""[Deposit] Helper, imported by the scripts of this directory. Paper: "Perpendicular billiards in rational polygons of odd denominator". Research name: probes/s941_frontier_template.py, copied verbatim except imports (deposit names), paths (data/ and figs/ under odd_polygons/) and this paragraph. Run through run_all.py; the usage lines below are the research repo's. Bracketed IDs and section marks (§sNNN) refer to the project's internal ledgers.

[promoted s958 from data/scratch/s941/frontier_template.py, verbatim but this line and the sys.path line; imported by
probes/s941_chain_frontier.py (per_n)] s941 scratch (s940-a step 1): the frontier of the m-th torus of Π_N read off §s935 (4)'s TEMPLATE alone (no tracer).
Exact Laurent arithmetic (s935_pi_n_torus.L, N-free exponents) for numeric n = m + 1.
Strip k of J0: copy g_k, slab X ∈ [lo0 − S_k, hi0 − S_k]; J1: [hi0 − S'_k, hi1 − S'_k]. Region R = union of slabs.
Exposed edge = a slab endpoint in copy g not shared by another slab of R in g. Paths: A = J1.lo ⊕ J0.lo, B = J0.hi ⊕ J1.hi.
Vertex touch = edge X equals X(g V_i) for an endpoint V_i of the strip's exit side (exit end) — or entry side at step 0."""
import _setup  # noqa: F401  (deposit: this directory on the path, cwd = odd_polygons/)
import sys, math, json, collections
import pi_n_torus as KP
L, X, g_apply, g_across, V, SIDES = KP.L, KP.X, KP.g_apply, KP.g_across, KP.V, KP.SIDES


def word(n, which):
    rep = lambda w, k: w * k
    Q = rep('23', 2 * n) + rep('10', n)
    w = rep('23', 2 * n) + '20' + rep('10', n - 1) + '13' + rep('23', 2 * n - 2) + '20' + rep('10', n - 1)
    return w + Q * (2 if which == 0 else 1) + '13'


def segs(n, which):
    o, e, s = 2 * n - 1, 1, 3
    out = []
    for ch in word(n, which):
        s2 = int(ch)
        out.append((o, e, s, s2))
        (o, e) = g_across((o, e), s2)
        s = s2
    assert (o, e, s) == (2 * n - 1, 1, 3), 'no return'
    return out


def walk(seg):
    S, out = L(), []
    for (o, e, s, s2) in seg:
        out.append(S)
        S = S + KP.shift((o, e), s2)
    return out, S


def ev(f, N):
    return f.evx(math.pi / N).real


def frontier(n, Nev=None):
    Nev = Nev or 8 * n + 41
    J0, J1 = segs(n, 0), segs(n, 1)
    S0, T0 = walk(J0); S1, T1 = walk(J1)
    xv = lambda seg, k, i: X(g_apply(seg[k][:2], V[i]))
    lo0 = xv(J0, 4 * n + 1, 0) + S0[4 * n + 1]
    hi0 = xv(J0, 18 * n - 2, 2) + S0[18 * n - 2]
    hi1 = xv(J1, 6 * n + 1, 0) + S1[6 * n + 1]
    lam1 = hi1 - hi0
    assert ((lo0 - T0) - (lo0 + lam1)).iszero() and ((hi0 - T1) - lo0).iszero(), 'returns'
    edges = {('J0', 'lo'): [lo0 - s for s in S0], ('J0', 'hi'): [hi0 - s for s in S0],
             ('J1', 'lo'): [hi0 - s for s in S1], ('J1', 'hi'): [hi1 - s for s in S1]}
    SEG = {'J0': J0, 'J1': J1}
    # slab endpoints per copy, keyed by float value (exactness checked on matches)
    per = collections.defaultdict(list)
    for (J, side), xs in edges.items():
        for k, x in enumerate(xs):
            per[SEG[J][k][:2]].append((ev(x, Nev), x, J, side, k))
    exposed = {}
    for g, lst in per.items():
        lst.sort(key=lambda r: r[0])
        for r in lst:
            mates = [q for q in lst if q is not r and abs(q[0] - r[0]) < 1e-9 and q[3] != r[3]]
            for q in mates:
                assert (q[1] - r[1]).iszero(), 'float match not exact'
            # overlap check: slabs in one copy must be disjoint
            exposed[(r[2], r[3], r[4])] = not mates
    # disjointness of slabs per copy
    for g, lst in per.items():
        sl = sorted([(ev(edges[(J, 'lo')][k], Nev), ev(edges[(J, 'hi')][k], Nev)) for (_, _, J, sd, k) in lst if sd == 'lo'])
        for (a, b), (c, d) in zip(sl, sl[1:]):
            assert b <= c + 1e-12, ('overlap', g)

    def touches(J, side, k):
        o, e, s, s2 = SEG[J][k]
        x = edges[(J, side)][k]
        t = []
        for which, ss in (('in', s), ('out', s2)):
            for i in SIDES[ss]:
                if (x - xv(SEG[J], k, i)).iszero():
                    t.append((which, i))
        return t

    paths = {'A': [('J1', 'lo', k) for k in range(len(J1))] + [('J0', 'lo', k) for k in range(len(J0))],
             'B': [('J0', 'hi', k) for k in range(len(J0))] + [('J1', 'hi', k) for k in range(len(J1))]}
    res = {}
    for name, P in paths.items():
        rows = []
        for idx, (J, sd, k) in enumerate(P):
            rows.append(dict(i=idx, J=J, side=sd, k=k, copy=SEG[J][k][:2], exp=exposed[(J, sd, k)], touch=touches(J, sd, k)))
        res[name] = rows
    return res


def connections(rows):
    """split a path at vertex touches; report each piece: endpoints, #crossings, exposed pattern"""
    cuts = []   # positions (between steps) where a vertex is touched: position p = boundary between step p-1 and p
    for r in rows:
        for which, i in r['touch']:
            p = r['i'] + (1 if which == 'out' else 0)
            cuts.append((p, i))
    cuts = sorted(set(cuts))
    pieces = []
    for (p, i), (q, j) in zip(cuts, cuts[1:]):
        seg = rows[p:q]
        ex = ''.join('E' if r['exp'] else '.' for r in seg)
        pieces.append(dict(frm=(p, f'V{i}'), to=(q, f'V{j}'), steps=q - p, bounces=q - p - 1, exposed=ex))
    return cuts, pieces


if __name__ == '__main__':
    for n in [int(a) for a in sys.argv[1:]]:
        res = frontier(n)
        print(f'=== n = {n} (m = {n - 1}): chords c_(4n-3..4n+3) = c{4*n-3}, c{4*n-1}, c{4*n+1}, c{4*n+3}; loops l_(m-1), l_m bounces {8*n-4}, {8*n+4}')
        for name, rows in res.items():
            cuts, pieces = connections(rows)
            ex = ''.join('E' if r['exp'] else '.' for r in rows)
            print(f' path {name} (len {len(rows)}): exposed map {ex}')
            print(f'   touches {cuts}')
            for pc in pieces:
                if 'E' in pc['exposed']:
                    print('   piece', pc['frm'], '->', pc['to'], 'bounces', pc['bounces'], pc['exposed'])


def holonomy_check(n):
    """exact developed holonomy of each exposed piece vs ±c2(z^{2b} − z^{−2b}) (chord c_b) or the loop sum; piece ID = (copy, vertex)."""
    J0, J1 = segs(n, 0), segs(n, 1)
    res = frontier(n)
    SEG = {'J0': J0, 'J1': J1}
    def dev(seg):
        C, out = L(), []
        for (o, e, s, s2) in seg:
            out.append(C)
            C = C + g_apply((o, e), KP.T[s2])
        return out
    D = {'J0': dev(J0), 'J1': dev(J1)}
    c2 = KP.c2
    chord = lambda b: c2 * (L.mono(2 * b) - L.mono(-2 * b))
    out = []
    for name, rows in res.items():
        cuts, pieces = connections(rows)
        for pc in pieces:
            if 'E' not in pc['exposed']:
                continue
            assert set(pc['exposed']) == {'E'}, 'mixed exposure inside a piece'
            p, q = pc['frm'][0], pc['to'][0]
            r0, r1 = rows[p], rows[q - 1]
            vi, vj = int(pc['frm'][1][1]), int(pc['to'][1][1])
            # path crosses J between J-blocks: developing restarts; sum holonomy blockwise
            hol = L(); blk = []
            for r in rows[p:q]:
                if not blk or blk[-1][0] != r['J'] or r['k'] != blk[-1][2] + 1:
                    blk.append([r['J'], r['k'], r['k']])
                else:
                    blk[-1][2] = r['k']
            for bi, (J, k0, k1) in enumerate(blk):
                seg = SEG[J]
                st = g_apply(seg[k0][:2], V[vi]) + D[J][k0] if bi == 0 else None
                en = g_apply(seg[k1][:2], V[vj]) + D[J][k1] if bi == len(blk) - 1 else None
                # block contribution: from start point (vertex or J-crossing) to end point; J-crossing positions cancel
                # across blocks because J is one line: use developed coordinates of the block's first/last side crossings
                first = st if st is not None else D[J][k0]        # crossing of side 3 of start copy at the block start: C = 0 frame
                last = en if en is not None else D[J][k1 + 1] if k1 + 1 < len(D[J]) else None
                if last is None:      # block runs to the end of its segment: end at the return crossing
                    C = D[J][k1] + g_apply(seg[k1][:2], KP.T[seg[k1][3]])
                    last = C
                hol = hol + (last - first)
            b = pc['steps']
            vert = (hol + hol.conj()).iszero()
            if vi != vj:
                tag = [f'c{b}' for sgn in (1, -1) if (hol - sgn * chord(b)).iszero()]
            else:
                m = (b - 12) // 8
                lp = chord(4 * m + 5) + chord(4 * m + 7)
                tag = [f'l{m}' for sgn in (1, -1) if (hol - sgn * lp).iszero()]
            out.append(dict(path=name, steps=b, frm=pc['frm'], to=pc['to'], vertical=vert, ident=tag,
                            id=(rows[p]['copy'], vi)))
    return out


def structure(n):
    """touch positions and, per non-exposed piece, its gluing partner at the piece's first step."""
    J0, J1 = segs(n, 0), segs(n, 1)
    res = frontier(n)
    out = {}
    for name, rows in res.items():
        cuts, pieces = connections(rows)
        out[name] = dict(cuts=cuts, pieces=[(pc['frm'], pc['to'], pc['exposed'][0]) for pc in pieces])
    return out
