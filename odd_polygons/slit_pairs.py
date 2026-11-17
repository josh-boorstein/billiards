"""[Deposit] Helper, imported by the scripts of this directory. Paper: "Perpendicular billiards in rational polygons of odd denominator". Research name: probes/s624_slit_pairs.py, copied verbatim except imports (deposit names), paths (data/ and figs/ under odd_polygons/) and this paragraph. Run through run_all.py; the usage lines below are the research repo's. Bracketed IDs and section marks (§sNNN) refer to the project's internal ledgers.

s624 -- s623-a (2): the BOUNDARY of each known minimal torus -- is it a pair of homologous saddle
connections (a slit), and what does it join?

PRE-REGISTERED.  [NCYL-631]: in the seven rotation-type counterexamples the minimal component M is a
flat torus C/Lambda minus a zero-area, flow-parallel hole.  A flow-parallel slit lies on one flow line
of the torus, so it is invisible to the first-return map; in X it must appear as TWO saddle
connections parallel to omega with the SAME holonomy l*omega and the same endpoint cone points
(the slit's two sides).  Here: every forward omega-separatrix of X is traced EXACTLY through the
one-crossing IET of engine/flow_decomp.py (a separatrix from vertex i of copy g, where omega points
into g, is the T-orbit of the common endpoint of g's two entry segments; it ends at a cone point
exactly when the orbit lands on a breakpoint -- every breakpoint is an edge endpoint), to a cap.
For each saddle connection: start / end cone point (vertex classes under the side gluings), number of
crossings, holonomy l (exact: l*omega = end - start - sum of gluing jumps).  Then: groups of saddle
connections with identical holonomy and endpoints.
COULD COME OUT OTHERWISE: no equal-holonomy pair (then the hole is not a simple slit -- a tree of
segments, or the model of [NCYL-631] (3) is wrong); a pair whose length is incompatible with M.
CHECK: separatrix count = d (omega points into exactly d of the 2d copies at one vertex each: the
total angle 2d*pi gives d outgoing omega-directions); every infinite (capped) separatrix should be one entering a minimal component.
PRIOR ART: [NCYL-631] (the torus), [NCYL-629]; Letter §6 (T0's O -> A connection, 13 reflections, on
M's boundary -- the T0 family is two-holed and is NOT in scope here); rulings.py --grep "homologous" /
"slit" / "saddle connection pair" -> none on the ten.  DATA: data_index --grep saddle -> none for these
triangles.
ADDED IN-SESSION (after the first (15,30,47) run): copy edges parallel to omega (the d - #sources
missing separatrices -- 90 of 92 at (15,30,47)); the junction fix (a concatenated junction value is
the lo of the next segment, never the hi of the last); the SLIT TEST (angle_table): at the start Q
of each capped separatrix, the incoming saddle connections at angle +-pi; a third argument runs an
off-grid direction of s624_transverse_rank (B) for the CP-triangle comparison.
RESULT (s624; data/s624_slit_pairs_<abc>_<dir>.json): all seven rotation-type counterexamples, 4
capped separatrices each (one per torus), every other separatrix a saddle connection.  At each Q the
two arrivals at +-pi bound ONE flow-parallel slit of length l: (15,30,47) 157.833703 and (11,45,50)
78.727685 are homologous loop pairs; the other five have one side a single saddle connection and the
other broken by cone points ON the slit, summing exactly: (13,20,23) 668.124605 = 617.58272 +
49.209192 + 1.332693; (23,29,58) 56.112591 = 54.891295 + 1.221296; (19,38,55) 308.809955 =
306.810742 + 1.999213; (11,22,103) 331.253539 = 283.779087 + 47.474452; (30,41,67) 59.29811 =
25.125185 + 20.097805 + 12.812943 + 1.262176 -- the short pieces 1.33 / 1.22 / 2.00 / 1.26 are COPY
EDGES of the triangle (in 4 of 7).  COMPARISON: off-grid tori of the CP triangles (4,5,9), (2,5,6)
show the same homologous-loop slits (7921.752243, 12707.42854).  -> [NCYL-632].
ADDED s625: each traced separatrix also records its 'word' (the side crossed at every copy), for
probes/s625_slit_family.py.
"""
import _setup  # noqa: F401  (deposit: this directory on the path, cwd = odd_polygons/)
import json
import sys
import time

import numpy as np

import flow_decomp as fd


def cone_classes(T):
    """Union-find: copy g vertex i ~ copy g_across(g, s) vertex i for each side s through i."""
    par = {}

    def f(x):
        while par.setdefault(x, x) != x:
            par[x] = par[par[x]]
            x = par[x]
        return x
    copies = [(k, e) for e in (0, 1) for k in range(T.d)]
    for g in copies:
        for s, (i, k2) in enumerate(T.sides):
            g2 = T.g_across(g, s)
            for v in (i, k2):
                a, b = f((g, v)), f((g2, v))
                if a != b:
                    par[a] = b
    cls = {}
    for g in copies:
        for v in range(3):
            cls.setdefault(f((g, v)), []).append((g, v))
    return {m: r for r, ms in cls.items() for m in ms}, cls


def trace_all(abc, j, cap):
    T = fd.Triangle(abc)
    F, N = T.F, T.N
    om = {e % N: c for e, c in j.items()} if isinstance(j, dict) else {j % N: 1}
    iom = fd.mshift(om, T.d, N)                      # i * omega (i = zeta^d)
    import cmath, math
    omabs = abs(sum(c * cmath.exp(2j * math.pi * k / N) for k, c in om.items()))
    X = lambda p: T.red(T.cross2(om, p))
    pieces, lengths, top, bot, meta = fd.build_iet(T, j)
    n = len(pieces)
    D = F.D
    dom = [np.array(pieces[p][0], dtype=object) for p in range(n)]
    img = [np.array(pieces[p][2], dtype=object) for p in range(n)]
    tr = [img[p] - dom[p] for p in range(n)]
    Jp = [np.array(T.red(T.g_apply(pieces[p][3][0], T.t[pieces[p][4][1]])), dtype=object) for p in range(n)]
    order = top                                     # pieces sorted by domain start
    starts_f = np.array([F.fval(pieces[p][0]) for p in order])
    # breakpoints (exact tuples) -> the (copy, vertex) sitting there, from the entry segments
    # rebuild segment table: key (g, s) -> [lo, hi] in concatenated coordinates
    segs = {}
    for p in range(n):
        key = pieces[p][3]
        a = pieces[p][0]
        b = F.add(a, pieces[p][1])
        if key not in segs:
            segs[key] = [a, b]
        else:
            if F.cmp(a, segs[key][0]) < 0:
                segs[key][0] = a
            if F.cmp(b, segs[key][1]) > 0:
                segs[key][1] = b
    # which vertex of copy g is at a segment endpoint: entry side s = (i, k2); the endpoint with the
    # smaller transverse coordinate x_g is lo
    brk = {}                                        # exact tuple -> list of (g, vertex, key, end)
    for (g, s), (lo, hi) in segs.items():
        i, k2 = T.sides[s]
        W = [T.g_apply(g, v) for v in T.V]
        xi, xk = X(W[i]), X(W[k2])
        vlo, vhi = (i, k2) if F.cmp(xi, xk) < 0 else (k2, i)
        brk.setdefault(tuple(lo), []).append((g, vlo, (g, s), 'lo'))
        brk.setdefault(tuple(hi), []).append((g, vhi, (g, s), 'hi'))
    # interior piece boundaries = far-vertex preimages (the vertex opposite the entry side)
    for p in range(n):
        a = tuple(pieces[p][0])
        if a not in brk:
            g, s = pieces[p][3]
            i, k2 = T.sides[s]
            brk.setdefault(a, []).append((g, 3 - i - k2, (g, s), 'mid'))
    Lt = tuple(meta['L'])
    vclass, classes = cone_classes(T)

    def locate(x):
        """index into `order` of the piece whose [dom, dom+len) contains x (exact)."""
        xf = float(sum(float(a) * c for a, c in zip(x, F.C) if a))
        k = int(np.searchsorted(starts_f, xf, side='right')) - 1
        k = max(0, min(k, n - 1))
        # exact adjust
        while k > 0 and F.cmp(tuple(x), pieces[order[k]][0]) < 0:
            k -= 1
        while k < n - 1 and F.cmp(tuple(x), pieces[order[k + 1]][0]) >= 0:
            k += 1
        return order[k]

    # separatrix starts: for each copy g, the vertex where omega points inward is the common endpoint
    # of its two entry segments.  Start point x0 = that endpoint (as the lo of one segment or hi of
    # the other); first crossing: the piece of the segment with that endpoint, translated.
    seps = []
    bykey = {}
    for (g, s), (lo, hi) in segs.items():
        bykey.setdefault(g, []).append((s, lo, hi))
    for g, lst in bykey.items():
        if len(lst) != 2:
            continue                                # one entry side: omega parallel to a side here
        W = [T.g_apply(g, v) for v in T.V]
        (s1, lo1, hi1), (s2, lo2, hi2) = lst
        vs1, vs2 = set(T.sides[s1]), set(T.sides[s2])
        v = (vs1 & vs2).pop()
        # the endpoint of segment 1 at vertex v
        i, k2 = T.sides[s1]
        vlo = i if F.cmp(X(W[i]), X(W[k2])) < 0 else k2
        if v == vlo:
            x0 = lo1                               # piece starting at lo1
            p = locate(x0)
        else:
            x0 = hi1                               # limit from the left: the piece ending at hi1
            p = None
            for q in range(n):
                if pieces[q][3] == (g, s1) and F.cmp(F.add(pieces[q][0], pieces[q][1]), hi1) == 0:
                    p = q
        seps.append(((g, v), np.array(x0, dtype=object), p))
    out = []
    t0 = time.time()
    for (g, v), x0, p in seps:
        x = x0 + tr[p]
        J = Jp[p].copy()
        steps = 1
        hit = None
        qprev = p
        word = [pieces[p][4][1]]                    # exit side of every copy crossed (s625)
        while steps < cap:
            key = tuple(x)
            if key in brk:
                ents = [e for e in brk[key] if e[3] in ('lo', 'mid')]
                hit = ents
                break
            q = locate(x)
            x = x + tr[q]
            J = J + Jp[q]
            word.append(pieces[q][4][1])
            qprev = q
            steps += 1
        rec = {'start': [list(g), v], 'steps': steps, 'end': None, 'word': word}
        rec['start_cone'] = str(vclass[(g, v)])
        if hit:
            assert len(hit) == 1, hit
            g2, v2, key2, end = hit[0]
            arr = (g2, v2) if end == 'mid' else (pieces[qprev][3][0], v2)
            # holonomy: l*omega = W_end (in copy g2 frame) - W_start (copy g frame) - sum of jumps;
            # jumps = -g t_s per crossing, so + J
            Ws = T.red(T.g_apply(g, T.V[v]))
            We = T.red(T.g_apply(g2, T.V[v2]))
            hol = tuple(a - b + c for a, b, c in zip(We, Ws, J))
            # along-flow length: 2 dot(omega, hol) = cross2(i omega, hol) with i = zeta^d
            Yh = T.red(T.cross2(iom, {k: c for k, c in enumerate(hol) if c}))
            Xh = T.red(T.cross2(om, {k: c for k, c in enumerate(hol) if c}))
            rec.update(arr=[list(arr[0]), arr[1]], end=[list(g2), v2], end_cone=str(vclass[(g2, v2)]), hol=[int(a) for a in hol],
                       transverse_zero=not any(Xh), length=-F.fval(Yh) / 2 / omabs)
        out.append(rec)
    # saddle connections along copy edges parallel to omega (one per glued edge)
    seen = set()
    for g in [(k, e) for e in (0, 1) for k in range(T.d)]:
        W = [T.g_apply(g, v) for v in T.V]
        for s_, (i, k2) in enumerate(T.sides):
            e = fd.madd(W[k2], fd.mneg(W[i]))
            if any(T.red(T.cross2(e, om))):
                continue
            g2 = T.g_across(g, s_)
            if (g2, s_) in seen:
                continue
            seen.add((g, s_))
            Y = F.fval(T.red(T.cross2(iom, e))) / omabs   # = -2 dot(omega, e)/|omega|
            tail, head = (i, k2) if Y < 0 else (k2, i)
            hol = T.red(e) if Y < 0 else tuple(-a for a in T.red(e))
            out.append({'start': [list(g), tail], 'start_cone': str(vclass[(g, tail)]), 'steps': 0,
                        'edge': True, 'arr': [list(g), head], 'end': [list(g), head],
                        'end_cone': str(vclass[(g, head)]), 'hol': [int(a) for a in hol],
                        'transverse_zero': True, 'length': abs(Y) / 2})
    return out, time.time() - t0, len(classes), angle_table(T, j)


def angle_table(T, j):
    """Cumulative angle of every (copy, vertex) sector around its cone point, walking ccw, and the
    angular position of +omega / -omega inside each sector (None if outside).  Floats: directions
    are multiples of pi/(2d), so +-pi matches are exact to rounding."""
    import cmath
    import math
    N = T.N
    z = lambda m: complex(sum(c * cmath.exp(2j * math.pi * k / N) for k, c in m.items()))
    om = (sum(c * cmath.exp(2j * math.pi * k / N) for k, c in j.items()) if isinstance(j, dict)
          else cmath.exp(2j * math.pi * j / N))
    om = om / abs(om)
    copies = [(k, e) for e in (0, 1) for k in range(T.d)]
    done, table = set(), {}
    for g0 in copies:
        for v in range(3):
            if (g0, v) in done:
                continue
            g, A, cone = g0, 0.0, []
            while True:
                W = [z(T.g_apply(g, u)) for u in T.V]
                sides = [s for s, (i, k2) in enumerate(T.sides) if v in (i, k2)]
                us = []
                for s in sides:
                    i, k2 = T.sides[s]
                    o = k2 if i == v else i
                    us.append((s, W[o] - W[v]))
                (sa, ua), (sb, ub) = us
                if (ua.conjugate() * ub).imag < 0:
                    (sa, ua), (sb, ub) = (sb, ub), (sa, ua)
                alpha = cmath.phase(ub / ua)
                def pos(w):
                    t = cmath.phase(w / ua)
                    return t if -1e-12 <= t <= alpha + 1e-12 else None
                cone.append(((g, v), A, alpha, pos(om), pos(-om)))
                done.add((g, v))
                A += alpha
                g = T.g_across(g, sb)
                if g == g0:
                    break
            tot = A
            for (gv, A0, alpha, po, pm) in cone:
                table[gv] = {'A': A0, 'alpha': alpha, 'out': None if po is None else A0 + po,
                             'in': None if pm is None else A0 + pm, 'total': tot}
    return table


def main():
    abc = tuple(int(a) for a in sys.argv[1].split(',')) if len(sys.argv) > 1 else (15, 30, 47)
    cap = int(float(sys.argv[2])) if len(sys.argv) > 2 else 10 ** 6
    if len(sys.argv) > 3:
        # COMPARISON: the idx-th off-grid direction of s624_transverse_rank (B) on abc that carries a
        # certified minimal component (a CP triangle's torus, for the anatomy comparison)
        idx = int(sys.argv[3])
        B = [r for r in json.load(open('data/s624_transverse_rank.json'))['B']
             if tuple(r['abc']) == abc and any(c['kind'] == 'no_periodic' for c in r.get('bad', []))]
        j = {int(k): int(c) for k, c in B[idx]['om'].items()}
        name = f'off{idx}'
    else:
        rows = json.load(open('data/s624_loop_holonomy.json'))['A']
        r = next(r for r in rows if tuple(r['abc']) == abc)
        name = r['dir']
        s = int(name[-1])
        j = (fd.perp_dirs(abc) if name.startswith('perp') else fd.side_dirs(abc))[s]
    out, dt, ncone, ang = trace_all(abc, j, cap)
    sc = [o for o in out if o['end']]
    print(abc, name, 'separatrices', len(out), 'saddle connections', len(sc), 'capped', len(out) - len(sc),
          'cone points', ncone, f'{dt:.0f}s')
    groups = {}
    for o in sc:
        groups.setdefault((tuple(o['hol']), o['start_cone'], o['end_cone']), []).append(o)
    for k, g in sorted(groups.items(), key=lambda kv: kv[1][0]['length']):
        print('  len %.6f  steps %s  x%d  %s -> %s  transverse0=%s' % (
            g[0]['length'], [o['steps'] for o in g], len(g), k[1], k[2], g[0]['transverse_zero']))
    # SLIT TEST: at the start Q of each capped (M-entering) separatrix, the incoming saddle
    # connections arriving at angle out +- pi
    import math
    slits = []
    for o in out:
        if o['end']:
            continue
        g, v = tuple(o['start'][0]), o['start'][1]
        a = ang[(g, v)]
        tot, th = a['total'], a['out']
        arr = []
        for o2 in out:
            if not o2['end']:
                continue
            b = ang[(tuple(o2['arr'][0]), o2['arr'][1])]
            if b['total'] != tot or o2['end_cone'] != o['start_cone']:
                continue
            if b['in'] is None and not o2.get('edge'):
                continue
            thin = b['in'] if b['in'] is not None else None
            if thin is None:
                continue
            dth = (thin - th) % tot
            if abs(dth - math.pi) < 1e-7 or abs(dth - (tot - math.pi)) < 1e-7:
                arr.append(o2)
        slits.append({'Q': o['start_cone'], 'Q_angle_over_2pi': round(tot / (2 * math.pi), 6),
                      'sides': [{'from': x['start_cone'], 'length': x['length'], 'steps': x['steps'],
                                 'hol': x['hol'], 'edge': x.get('edge', False)} for x in arr]})
        print('  SLIT? Q =', o['start_cone'], 'angle/2pi', round(tot / (2 * math.pi), 3), '->',
              [(x['start_cone'], round(x['length'], 6), x['steps'], x.get('edge', False)) for x in arr],
              'equal hol:', len(arr) == 2 and arr[0]['hol'] == arr[1]['hol'])
    json.dump({'abc': abc, 'dir': name, 'j': j, 'cap': cap, 'seps': out, 'slits': slits},
              open(f'data/s624_slit_pairs_{abc[0]}_{abc[1]}_{abc[2]}_{name}.json', 'w'), indent=1, default=str)


if __name__ == '__main__':
    main()
