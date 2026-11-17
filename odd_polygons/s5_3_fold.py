"""[Deposit] Proposition 5.3(b) and Table 1: the fold obstruction on the catalogue, the octagons. Paper: "Perpendicular billiards in rational polygons of odd denominator". Research name: probes/s933_reflection_tilings.py, copied verbatim except imports (deposit names), paths (data/ and figs/ under odd_polygons/) and this paragraph. Run through run_all.py; the usage lines below are the research repo's. Bracketed IDs and section marks (§sNNN) refer to the project's internal ledgers.

s933 keeper (queue `s932-b`): odd-N polygons tiled by reflected copies of ONE catalogued non-CP triangle.

PRE-REGISTERED (s933, before the tiling run).
QUESTION. A polygon Π tiled edge-to-edge by copies of a triangle T, every interior edge a MIRROR between its two tiles
(a *reflection tiling*), folds onto T: each tile maps isometrically to T and the fold carries Π's billiard onto T's,
so X_Π is a finite branched translation cover of X_T and a direction class of Π is CP iff the T class it folds onto is
(a periodic orbit lifts to finitely many (tile, point, direction) states, hence periodic; a non-periodic one lifts to
non-periodic ones) — the [NCYL-627] (1) / §s929 (2) fold argument, which there had one mirror (isosceles; kite).
Is there such a Π with ODD N_Π, NO mirror symmetry (no isometry at all), whose perpendicular class folds onto a
non-CP class of a catalogued T (`data/s626_noncp_catalogue.json`, 28 rows)?
THE LOCAL OBSTRUCTION (hand, checked here as `feasibility`). Work in units π/N_e (N_e = T's N, even for 27 rows). The
line directions of T's sides are d0 = 0, d1 = −a1, d2 = a0 (mod N_e) and the reflection group acts on line directions
by x ↦ ±x + 2ℤ, so there are exactly TWO direction classes, the even and the odd directions. (i) At odd N_Π every
Π angle k·a_i is ≡ 0 mod 2^e (e = v₂(N_e)), so all Π sides have one parity σ, and Π's perpendicular class folds onto
the class of parity σ + N_e/2: σ is forced by the bad parity β, σ = β + N_e/2. (ii) Sides of T with parity ≠ σ can
never be boundary. (iii) A vertex point of type i is interior only with a full fan of 2N_e/a_i tiles, which closes
under reflections iff a_i | N_e; on the boundary its fan of k tiles alternates T's two sides at V_i, so its two end
edges are boundary-type: k even if only one side at V_i has parity σ, any k if both; and k·a_i < 2N_e, k·a_i ≡ 0
mod 2^e. A vertex type with no admissible boundary fan and a_i ∤ N_e kills every tiling.
PREDICTION (from the obstruction, computed before the search): 25/28 rows obstructed; feasible exactly the three
kite parents of §s929 (2) — (30,41,67)/138, (23,46,201)/270, (45,90,119)/254 — with B = the two sides at the even
vertex, the third side always interior (tiles pair into s929's kites), no interior vertex possible.
SEARCH. Kite trees (one kite = T ∪ its mirror across the non-B side) glued along B-edges, up to `--kites` kites;
reject overlaps, non-mirror contacts (a shared edge whose endpoint labels do not match), T-junctions, pinched
boundaries. For each embedded tiling: N_Π, the isometry group of the boundary (cyclic (length, angle) sequence vs
its rotations and reversals), convexity; and a `Poly` rebuild (exact closure ⊕ simplicity at 60 digits).
COULD HAVE COME OUT OTHERWISE: every multi-kite tree overlapping or pinching (the B-vertex angles are near-flat or
reflex: 2·67 = 134 of 138 units, 2·201 = 402 of 270), or every embedded one mirror-symmetric — then the answer is
NO at small tile counts and the kite would be the only odd-N tiled carrier.
PRIOR ART: §s929 (2) (the kite = the 2-tile case, degree-1 fold, three parents selected), [NCYL-627] (1) (isosceles
fold), §s928 (4)(b); `rulings.py --grep` *tiled* → other objects, *reflection tiling* → none, *fold* → other objects;
`general_triangles.md` grep *tiled by* → §s756 (isosceles only). → no multi-tile reflection tiling recorded.
DATA: `data_index.py --grep tiling` → s293/s374 cells (other objects); *kite* → s451/s461/s469 (other objects) → none.
RESULT (s933): prediction held — 25/28 obstructed, the three kite parents feasible. Up to 4 kites: 15 / 4 / 16 congruence
classes, 7 / 0 / 7 with NO isometry; smallest mirror-free: octagons at 3 kites, (134,30,67,15,82,30,41,15)/69,
(82,15,134,15,41,45,67,15)/69, (90,238,45,45,90,90,45,119)/127, (45,90,45,119,135,45,45,238)/127 — all non-convex (no
tiling with ≥ 3 kites is convex; the symmetric 2-kite hexagon of (45,90,119)/254 is). Consistency (`--classify 2`): both N = 69 octagons s' = 12 at caps 3000/20000, 18/17 of 200
free orbits open, 99 pairs, 34/24 slits. → general_triangles.md §s933 (1), [POLY-FOLD-ODD], [NCYL-721] ↓s933.
Output: data/s933_reflection_tilings.json. Run from the repo root with the project PYTHONPATH (`--kites K`, `--classify M`:
s930's classifier on the first M mirror-free tilings per parent — minutes at N = 69, much longer at N = 127).
"""
import _setup  # noqa: F401  (deposit: this directory on the path, cwd = odd_polygons/)
import itertools
import json
import math
import sys
from math import gcd

CAT = 'data/s626_noncp_catalogue.json'
EPS = 1e-9


def v2(n):
    e = 0
    while n % 2 == 0:
        n //= 2
        e += 1
    return e


def side_dirs(a, Ne):
    return [0, (-a[1]) % Ne, a[0] % Ne]


def bad_parity(v):
    a, Ne = v['abc'], v['d']
    d = side_dirs(a, Ne)
    bad, good = set(), set()
    for x in v['dirs']:
        s = x['side']
        par = (d[s] + (Ne // 2 if x['dir'].startswith('perp') else 0)) % 2
        (good if all(c == 'cylinder' for c in x['counts']) else bad).add(par)
    assert not bad & good
    return bad


def feasibility(v):
    a, Ne = v['abc'], v['d']
    if Ne % 2:
        return [{'skip': 'odd N_e'}]
    d, e, out = side_dirs(a, Ne), v2(Ne), []
    for beta in sorted(bad_parity(v)):
        sig = (beta + Ne // 2) % 2
        B = [j for j in range(3) if d[j] % 2 == sig]
        ok, notes = True, []
        for i in range(3):
            nb = sum(s in B for s in (i, (i - 1) % 3))
            interior = Ne % a[i] == 0
            ks = [k for k in range(1, 2 * Ne) if k * a[i] < 2 * Ne and (nb == 2 or k % 2 == 0)
                  and (k * a[i]) % (2 ** e) == 0] if nb else []
            if not ks and not interior:
                ok = False
                notes.append(f'V{i}({a[i]}): ' + ('must be interior, ' if nb == 0 else 'no boundary fan, ')
                             + f'{a[i]}∤{Ne}')
            else:
                notes.append(f'V{i}({a[i]}): fans {ks[:8]}' + (' or interior' if interior else ''))
        out.append({'beta': beta, 'sigma': sig, 'B': B, 'feasible': ok, 'notes': notes})
    return out


# ---------------------------------------------------------------- geometry (float search; Poly re-check at 60 digits)
def reflect(p, q, x):
    dx, dy = q[0] - p[0], q[1] - p[1]
    t = ((x[0] - p[0]) * dx + (x[1] - p[1]) * dy) / (dx * dx + dy * dy)
    fx, fy = p[0] + t * dx, p[1] + t * dy
    return (2 * fx - x[0], 2 * fy - x[1])


def tri_reflect(tile, j):
    """mirror of a labelled tile across its side j = V_j V_{j+1}"""
    V = list(tile)
    V[(j + 2) % 3] = reflect(V[j], V[(j + 1) % 3], V[(j + 2) % 3])
    return tuple(V)


def key_pt(p):
    return (round(p[0] / EPS / 1000), round(p[1] / EPS / 1000))


def key_tile(t):
    return tuple(key_pt(p) for p in t)


def area2(a, b, c):
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])


def tris_overlap(t, s):
    """interiors intersect (separating-axis on 2D triangles, with tolerance)"""
    for P, Q in ((t, s), (s, t)):
        for i in range(3):
            p, q = P[i], P[(i + 1) % 3]
            nx, ny = q[1] - p[1], p[0] - q[0]
            pa = [nx * x[0] + ny * x[1] for x in P]
            qa = [nx * x[0] + ny * x[1] for x in Q]
            L = math.hypot(nx, ny)
            if max(pa) <= min(qa) + 1e-7 * L or max(qa) <= min(pa) + 1e-7 * L:
                return False
    return True


def on_open_segment(x, p, q):
    if abs(area2(p, q, x)) > 1e-8 * math.dist(p, q):
        return False
    t = ((x[0] - p[0]) * (q[0] - p[0]) + (x[1] - p[1]) * (q[1] - p[1])) / math.dist(p, q) ** 2
    return 1e-7 < t < 1 - 1e-7


def triangle(a, Ne):
    A = [math.pi * x / Ne for x in a]
    # V0 at origin, side 0 along +x, CCW; |V0V1| = sin(a2)
    L0 = math.sin(A[2])
    V0, V1 = (0.0, 0.0), (L0, 0.0)
    L2 = math.sin(A[1])
    V2 = (L2 * math.cos(A[0]), L2 * math.sin(A[0]))
    return (V0, V1, V2)


def analyse(tiles, a, Ne, B):
    """validity of a tile set as a reflection tiling with boundary types ⊆ B; returns the boundary or a reason"""
    tl = list(tiles.values())
    for t, s in itertools.combinations(tl, 2):
        if tris_overlap(t, s):
            return None, 'overlap'
    edges = {}
    for ti, t in enumerate(tl):
        for j in range(3):
            p, q = t[j], t[(j + 1) % 3]
            k = frozenset((key_pt(p), key_pt(q)))
            edges.setdefault(k, []).append((ti, j, key_pt(p), key_pt(q)))
    for k, lst in edges.items():
        if len(lst) > 2:
            return None, 'edge>2'
        if len(lst) == 2:
            (t1, j1, p1, q1), (t2, j2, p2, q2) = lst
            if j1 != j2 or p1 != p2 or q1 != q2:
                return None, 'non-mirror contact'
        elif lst[0][1] not in B:
            return None, 'non-B side on boundary'
    pts = {key_pt(p): p for t in tl for p in t}
    for t in tl:
        for j in range(3):
            for x in pts.values():
                if on_open_segment(x, t[j], t[(j + 1) % 3]):
                    return None, 'T-junction'
    # orient consistently (a reflected tile is CW): the boundary keeps the tiling on its left
    succ = {}
    for lst in edges.values():
        if len(lst) == 1:
            ti, j, p, q = lst[0]
            t = tl[ti]
            if area2(*t) < 0:
                p, q = q, p
            if p in succ:
                return None, 'pinched'
            succ[p] = (q, j)
    start = next(iter(succ))
    cyc, cur = [], start
    for _ in range(len(succ) + 1):
        q, j = succ[cur]
        cyc.append((cur, q, j))
        cur = q
        if cur == start:
            break
    if len(cyc) != len(succ):
        return None, 'boundary not one cycle'
    return cyc, pts


def boundary_data(cyc, pts, Ne):
    """interior angles (units π/Ne, rounded exactly) and lengths, merging collinear edges"""
    P = [pts[c[0]] for c in cyc]
    n = len(P)
    ang, L = [], []
    for i in range(n):
        p0, p1, p2 = P[i - 1], P[i], P[(i + 1) % n]
        a1 = math.atan2(p1[1] - p0[1], p1[0] - p0[0])
        a2 = math.atan2(p2[1] - p1[1], p2[0] - p1[0])
        turn = (a2 - a1 + math.pi) % (2 * math.pi) - math.pi      # left turn > 0 (CCW boundary)
        interior = math.pi - turn
        u = interior * Ne / math.pi
        assert abs(u - round(u)) < 1e-6, u
        ang.append(round(u))
        L.append(math.dist(p1, p2))
    # merge straight points
    keep = [i for i in range(n) if ang[i] != Ne]
    if len(keep) < n:
        raise NotImplementedError('straight boundary point')
    return ang, L


def symmetries(ang, L):
    """isometries of the polygon: rotations / reflections of the cyclic (angle, length) sequence"""
    n = len(ang)
    seq = [(ang[i], round(L[i], 7)) for i in range(n)]          # angle at V_i, then side V_i V_{i+1}
    out = []
    for r in range(n):
        if r and all(seq[(i + r) % n] == seq[i] for i in range(n)):
            out.append(('rot', r))
    # reversal: angle at V_i ↔ angle at V_{-i}, side V_i V_{i+1} ↔ side V_{-i-1} V_{-i}
    for r in range(n):
        if all(ang[(r - i) % n] == ang[i] and round(L[(r - i - 1) % n], 7) == round(L[i], 7) for i in range(n)):
            out.append(('mirror', r))
    return out


def canon(ang, types):
    """congruence-class key: min over rotations ⊕ reversals of the cyclic (angle at V_i, type of V_i V_{i+1})"""
    n = len(ang)
    cands = []
    for r in range(n):
        cands.append(tuple((ang[(i + r) % n], types[(i + r) % n]) for i in range(n)))
        cands.append(tuple((ang[(r - i) % n], types[(r - i - 1) % n]) for i in range(n)))
    return min(cands)


def search(a, Ne, B, kmax):
    nonB = [j for j in range(3) if j not in B]
    T = triangle(a, Ne)
    if area2(*T) < 0:
        raise RuntimeError
    def close_kite(tiles):
        # add the mirror across every non-B side (the forced partner)
        changed = True
        while changed:
            changed = False
            for t in list(tiles.values()):
                for j in nonB:
                    s = tri_reflect(t, j)
                    if key_tile(s) not in tiles:
                        tiles[key_tile(s)] = s
                        changed = True
                        if len(tiles) > 2 * kmax + 4:
                            return tiles
        return tiles
    start = close_kite({key_tile(T): T})
    seen = {frozenset(start)}
    frontier = [start]
    results = []
    for size in range(1, kmax + 1):
        new = []
        for tiles in frontier:
            cyc, info = analyse(tiles, a, Ne, B)
            if cyc is not None:
                ang, L = boundary_data(cyc, info, Ne)
                types = [c[2] for c in cyc]
                g = gcd(Ne, *ang)
                results.append({'kites': size, 'tiles': len(tiles), 'k': len(ang), 'angles_Ne': ang,
                                'types': types, 'canon': canon(ang, types),
                                'N_Pi': Ne // g, 'angles': [x // g for x in ang], 'lengths': L,
                                'convex': all(x < Ne for x in ang), 'symmetries': symmetries(ang, L),
                                'tile_list': [list(map(list, t)) for t in tiles.values()]})
            else:
                if info in ('overlap', 'edge>2', 'non-mirror contact', 'T-junction'):
                    continue            # adding tiles never repairs these
            if size == kmax:
                continue
            for t in list(tiles.values()):
                for j in B:
                    s = tri_reflect(t, j)
                    if key_tile(s) in tiles:
                        continue
                    nt = dict(tiles)
                    nt[key_tile(s)] = s
                    nt = close_kite(nt)
                    fk = frozenset(nt)
                    if fk in seen:
                        continue
                    seen.add(fk)
                    new.append(nt)
        frontier = new
    return results


def make_poly(a, Ne, r):
    """Π as a `Poly` at 60 digits: exact side lengths sin(π a_{j+2}/Ne) by type (law of sines, |V_jV_{j+1}| ∝ sin a_{j+2})"""
    import mpmath as mp
    from poly_slots import Poly
    L = [mp.sin(mp.pi * a[(j + 2) % 3] / Ne) for j in r['types']]
    return Poly(r['angles'], r['N_Pi'], lengths=L)          # asserts closure at 1e-30 and simplicity


def classify(P, nfree):
    import free_scan as F
    cl = F.classify(P)
    clo = F.closure(P, nfree)
    return {'npairs': len(cl['pairs']), 'pairs': cl['pairs'], 'equal_len': cl['equal_len'], 'slits': cl['slits'],
            'sprime': cl['sprime'], 'kinds': cl['kinds'], 'closure': clo, 'open_free': clo['free'] - clo['closed']}


def main():
    kmax = int(sys.argv[sys.argv.index('--kites') + 1]) if '--kites' in sys.argv else 4
    ncl = int(sys.argv[sys.argv.index('--classify') + 1]) if '--classify' in sys.argv else 0
    cat = json.load(open(CAT))
    feas = {}
    for v in cat.values():
        feas[f"{tuple(v['abc'])}/{v['d']}"] = feasibility(v)
    nfeas = [k for k, f in feas.items() if any(x.get('feasible') for x in f)]
    print(f'feasible rows: {len(nfeas)}/{len(feas)}: {nfeas}')
    out = {'kmax': kmax, 'feasibility': feas, 'tilings': {}, 'classified': {}}
    for v in cat.values():
        k = f"{tuple(v['abc'])}/{v['d']}"
        for f in feas[k]:
            if not f.get('feasible'):
                continue
            res = search(v['abc'], v['d'], f['B'], kmax)
            byc = {}
            for r in res:
                byc.setdefault(r['canon'], r)
            reps = sorted(byc.values(), key=lambda r: (r['kites'], r['canon']))
            out['tilings'][k] = [{x: r[x] for x in r if x not in ('canon', 'tile_list', 'lengths')} |
                                 {'tile_list': r['tile_list']} for r in reps]
            for r in reps:
                print(f"  {k} kites={r['kites']} {r['k']}-gon N_Pi={r['N_Pi']} convex={r['convex']} "
                      f"sym={r['symmetries'] or 'NONE'} angles={r['angles']} types={r['types']}")
            free = [r for r in reps if not r['symmetries']]
            print(f"  {k}: {len(reps)} congruence classes, {len(free)} with no isometry")
            for r in free[:ncl]:
                P = make_poly(v['abc'], v['d'], r)
                c = classify(P, 200)
                out['classified'][f"{k}:{tuple(r['angles'])}/{r['N_Pi']}"] = c
                print(f"    classify {tuple(r['angles'])}/{r['N_Pi']}: s'={c['sprime']} pairs={c['npairs']} "
                      f"slits={c['slits']} closure={c['closure']}", flush=True)
                json.dump(out, open('data/s933_reflection_tilings.json', 'w'), indent=1, default=str)
    json.dump(out, open('data/s933_reflection_tilings.json', 'w'), indent=1, default=str)


if __name__ == '__main__':
    main()
