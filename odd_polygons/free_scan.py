"""[Deposit] Helper, imported by the scripts of this directory. Paper: "Perpendicular billiards in rational polygons of odd denominator". Research name: probes/s930_free_scan.py, copied verbatim except imports (deposit names), paths (data/ and figs/ under odd_polygons/) and this paragraph. Run through run_all.py; the usage lines below are the research repo's. Bracketed IDs and section marks (§sNNN) refer to the project's internal ledgers.

s930 keeper (queue `s929-e`, `s929-d`): free regions at isolated shapes of asymmetric odd-N polygons — which
TOPOLOGY do they have, and is any of them a MINIMAL component (an asymmetric non-CP polygon)?

PRE-REGISTERED (s930, before any row was run).
  Mechanism (hand, s930): a free region at an isolated shape t* is bounded by swapped connections aligned at t*; with
  e_L − e_R = ∂W in H_1(X, Σ) the alignment costs ONE condition when W's frontier is (i) two LOOPS (a cylinder, or a
  genus >= 1 surface with two boundary circles) or (ii) a SLIT: two connections P -> Q with equal holonomy (W of genus
  >= 1, Gauss–Bonnet: angle 2π at each slit end, no interior cone point ⟹ a slit torus). Case (i) at genus 0 is periodic
  automatically (s929's witnesses); in case (ii) the inner torus is periodic iff its two transverse periods are
  commensurable, which nothing forces at an algebraic t* — so a slit-bounded free region at an isolated t* is expected
  to be MINIMAL. The even-N counterexample tori are exactly slit tori ([NCYL-632] (2): partner τγ); at odd N homology
  must supply the partner. A single slit-bounded free region with a capped prong (s' > 0) at an asymmetric t* answers
  `s929-e` POSITIVELY.
  Target under test (T): at odd N, every free region at a flagged t* of an asymmetric polygon family is a union of
  CYLINDERS (every sampled free orbit closes) and s' = 0 — i.e. no asymmetric non-CP polygon on these families.
  Refuted by one flagged t* with s' > 0 that survives the near-vertex re-walk at a higher cap.
  Secondary (`s929-d`): on s929's three families with NO free t* at a 40-point grid — (2,7,8,13)/15, (1,2,5,6)/7,
  (3,4,5,6)/9 — does a finer grid (ngrid as run) find free t*? Count per family; a count of 0 is NOT evidence of
  absence (one change per prong per grid interval is found).
  Method: family = a 1-parameter SLICE of the closing lengths: L_free = base + t·e_1 (k − 2 free sides, the last two
  by closing); prong signatures (kind, end, steps) on an ngrid-point grid over the central 80% of the valid interval;
  each signature change bisected with cap = max(steps) + 2 (a short walk: the prong passes a vertex within that many
  bounces — no walk-capped bisection, s929 Don't); t* clustered to 1e-20 (AMENDED s930: 1e-30 split one t* into 2–3 clusters, the bisections stopping at different
  edges of the 1e-27 vertex window; the 'd' run used 1e-30 — duplicates removed at harvest); vertex prongs at t* counted by the
  NEAR-VERTEX walk (s929_asym_witness.walk, near = 1e-20; AMENDED s930 after the first control: 1e-28 sat on the
  edge of Poly's 1e-27 vertex window, where the bisection actually stops); flagged = >= 3 vertex
  prongs; at flagged t*: area test R (s928_poly_thmP.row, c = 1); if 1 − R > 1e-15: pairs (near walk: endpoints,
  loop / P -> Q, lengths; SLIT = two non-loop pairs with the same endpoint pair and lengths equal to 1e-25), s' at
  cap CAPS (near walk, every slot), and free-orbit closure: `nfree` points sampled uniformly in the table × class
  direction, walked <= 400 bounces; free ones re-walked <= 4000 bounces for a return to (side, direction, s ± 1e-20).
  CONTROLS: (C1, fires) `(5,2,5,2)/7` ρ ∈ [1.05, 1.15] — free regions at rational ρ (s929 A3), rotation-like (A2): the
  closure test must report them; (C2) the s929 witnesses (2,4,5,7)/9, (4,5,9,12)/15 re-found as CYLINDERS with
  loops / loop ⊕ chain frontier (the classifier's positive case for (i)).
  Could come out otherwise: a free t* whose free orbits do not close within 4000 bounces and with s' > 0 — the
  minimal case; or a slit frontier with all orbits closing (then the commensurability heuristic is wrong).
PRIOR ART: §s929 (1) ⊕ [NCYL-720] (scan design, the two witnesses, instrument notes (i)–(iii)), [NCYL-719] (*a slit, or
  equal-length loops* — Theorem P's residual), [NCYL-632] (2) (even-N minimal tori attached along one flow-parallel
  slit), §s928 (2) (zero-one law); `rulings.py --grep` *slit torus* → [NCYL-631], [NCYL-652], [NCYL-674] (even-N /
  right-triangle objects); *pentagon* → [LIT-002], [NCYL-709], [NCYL-720] (no pentagon family scan).
DATA: `data/s929_asym_free.json` (flagged rows only, quadrilaterals, no topology) → not reusable; `data_index.py --grep
  pentagon` → none.
Output: data/s930_free_scan_<tag>.json (one file per run tag).  Usage:
  python probes/s930_free_scan.py TAG [ngrid]     TAG in FAMILY_SETS.
"""
import _setup  # noqa: F401  (deposit: this directory on the path, cwd = odd_polygons/)
import collections
import json
import random
import sys
import time

import mpmath as mp

from poly_slots import Poly
from area_test import row
import s4_7_witnesses as W

mp.mp.dps = 50
NEAR = mp.mpf(10) ** -20   # s930: the bisection stops at the edge of Poly's vertex window (TOL·1e3 = 1e-27), not at TOL below;
                            # 1e-28 sat ON that edge and missed pairs at the parallelogram control (p = 0 at 19/18 + 5e-29)
TOL = mp.mpf(10) ** -38

# (angles, N, base free lengths (k-2 of them; index 1 is the slice parameter), optional explicit t-window)
FAMILY_SETS = {
    'ctrl': [([5, 2, 5, 2], 7, [1, 0], ('1.05', '1.15')),
             ([2, 4, 5, 7], 9, [1, 0], ('0.0950', '0.0960')),
             ([4, 5, 9, 12], 15, [1, 0], ('0.5400', '0.5415'))],
    'd': [([2, 7, 8, 13], 15, [1, 0]), ([1, 2, 5, 6], 7, [1, 0]), ([3, 4, 5, 6], 9, [1, 0])],
    'quad': [([1, 3, 4, 6], 7, [1, 0]), ([2, 3, 6, 7], 9, [1, 0]), ([2, 5, 7, 8], 11, [1, 0]),
             ([3, 5, 9, 13], 15, [1, 0]), ([1, 4, 6, 7], 9, [1, 0])],
    'pilot': [([3, 4, 5, 7, 8], 9, 'seed:1')],
    'pent': [([3, 4, 5, 7, 8], 9, 'seed:1'), ([3, 4, 5, 7, 8], 9, 'seed:2'),
             ([2, 4, 6, 7, 8], 9, 'seed:1'), ([4, 5, 7, 8, 9], 11, 'seed:1'),
             ([3, 7, 5, 8, 4], 9, 'seed:1'), ([3, 6, 7, 8, 9], 11, 'seed:1')],
}


def unit(d, N):
    th = mp.pi * d / (2 * N)
    return (mp.cos(th), mp.sin(th))


def closing(a, N, free):
    """free = the first k-2 side lengths; the last two from closure"""
    k = len(a)
    idx = [0]
    for j in range(1, k):
        idx.append(idx[-1] + 2 * (N - a[j]))
    u = [unit(x % (4 * N), N) for x in idx]
    sx = sum(free[j] * u[j][0] for j in range(k - 2))
    sy = sum(free[j] * u[j][1] for j in range(k - 2))
    A = mp.matrix([[u[k - 2][0], u[k - 1][0]], [u[k - 2][1], u[k - 1][1]]])
    sol = mp.lu_solve(A, mp.matrix([-sx, -sy]))
    return list(free) + [sol[0], sol[1]]


def base_of(fam):
    """fam[2] = explicit free lengths, or 'seed:n' = a valid random shape of Poly (normalised L0 = 1), sliced through
    with the parameter OFFSET from its own L1 (so t = 0 is the seed shape)"""
    if isinstance(fam[2], str):
        P = Poly(fam[0], fam[1], seed=int(fam[2].split(':')[1]))
        L = [x / P.L[0] for x in P.L]
        return L[:len(fam[0]) - 2]
    return [mp.mpf(x) for x in fam[2]]


def lengths(fam, t):
    base = base_of(fam)
    base[1] = base[1] + t
    return closing(fam[0], fam[1], base)


def make(fam, t):
    L = lengths(fam, t)
    if min(L) <= 0:
        return None
    try:
        return Poly(fam[0], fam[1], lengths=L)
    except (AssertionError, ValueError, ZeroDivisionError, TypeError):
        return None


def window(fam, margin=mp.mpf('0.1')):
    if len(fam) > 3:
        return mp.mpf(fam[3][0]), mp.mpf(fam[3][1])
    l0, l1 = lengths(fam, mp.mpf(0)), lengths(fam, mp.mpf(1))
    lo, hi = mp.mpf(-50), mp.mpf(50)
    for j in range(len(l0)):
        sl = l1[j] - l0[j]
        if sl > 0:
            lo = max(lo, -l0[j] / sl)
        elif sl < 0:
            hi = min(hi, -l0[j] / sl)
    # simple-polygon range: shrink to the valid sub-interval containing the midpoint
    w = hi - lo
    return lo + margin * w, hi - margin * w


def sig(P, cap):
    return {key: P.prong(*key, cap) for key in P.slots(P.classes()[0])}


def short(r):
    return (r[0], r[1], r[3])


def bisect(fam, key, t0, t1, s0, cap):
    lo, hi = t0, t1
    for _ in range(400):
        m = (lo + hi) / 2
        P = make(fam, m)
        if P is None:
            return None
        if short(P.prong(*key, cap)) == s0:
            lo = m
        else:
            hi = m
        if hi - lo < TOL:
            return (lo + hi) / 2
    return (lo + hi) / 2


def closure(P, nfree, seed=930, cap0=400, cap1=4000):
    """sample uniform (point, class direction); free = no perpendicular collision within cap0; for free ones,
    re-walk <= cap1 bounces looking for a return to (side, direction, s +- 1e-20)"""
    N, C = P.N, sorted(P.classes()[0])
    xs, ys = [v[0] for v in P.V], [v[1] for v in P.V]

    def inside(x, y):
        c = False
        for i in range(P.k):
            (x1, y1), (x2, y2) = P.V[i], P.V[(i + 1) % P.k]
            if (y1 > y) != (y2 > y) and x < x1 + (y - y1) * (x2 - x1) / (y2 - y1):
                c = not c
        return c
    rnd = random.Random(seed)
    tried = free = closed = late = 0
    periods = collections.Counter()
    while free < nfree and tried < 60 * nfree:
        tried += 1
        while True:
            x = mp.mpf(rnd.uniform(float(min(xs)), float(max(xs))))
            y = mp.mpf(rnd.uniform(float(min(ys)), float(max(ys))))
            if inside(x, y):
                break
        d, Q, first, outcome = rnd.choice(C), (x, y), None, 'open'
        for st in range(cap1):
            h = P.first_hit(Q, d, None)
            if h[0] == 'vertex':
                outcome = ('hit', st)
                break
            _, j, X, s = h
            if (d - P.m[j] - N) % (2 * N) == 0:
                outcome = ('hit', st)
                break
            if first is None:
                first = (j, d, s)
            elif j == first[0] and d == first[1] and abs(s - first[2]) < mp.mpf(10) ** -20:
                outcome = ('closed', st)
                break
            d, Q = (2 * P.m[j] - d) % (4 * N), X
        if outcome[0] == 'hit' and outcome[1] < cap0:
            continue                       # not free
        if outcome[0] == 'hit':
            late += 1                      # perpendicular only after cap0 bounces: reported, not counted free
            continue
        free += 1
        if outcome[0] == 'closed':
            closed += 1
            periods[outcome[1]] += 1
    return {'tried': tried, 'free': free, 'closed': closed, 'late_hit': late, 'periods': dict(periods.most_common(6))}


def classify(P, caps=(3000, 20000)):
    prs, kinds = W.pairs(P, NEAR)
    lens = [mp.mpf(q['length']) for q in prs]
    slits = []
    for i in range(len(prs)):
        for j in range(i + 1, len(prs)):
            a, b = prs[i], prs[j]
            if a['loop'] or b['loop']:
                continue
            ea = {a['slot'][0], a['end_vertex']}
            eb = {b['slot'][0], b['end_vertex']}
            if ea == eb and abs(lens[i] - lens[j]) < mp.mpf(10) ** -25:
                slits.append((i, j))
    eq = [(i, j) for i in range(len(prs)) for j in range(i + 1, len(prs)) if abs(lens[i] - lens[j]) < mp.mpf(10) ** -25]
    sp = {}
    for cap in caps:
        res = {key: W.walk(P, *key, NEAR, cap=cap) for key in P.slots(P.classes()[0])}
        sp[cap] = sum(1 for r in res.values() if r[0] == 'cap')
        if sp[cap] == 0:
            break
    return {'pairs': [{'from': q['slot'][0], 'to': q['end_vertex'], 'loop': q['loop'], 'len': q['length'][:20],
                       'bounces': q['bounces']} for q in prs],
            'equal_len': eq, 'slits': slits, 'sprime': {str(k): v for k, v in sp.items()}, 'kinds': kinds}


def scan(fam, ngrid, cap=3000, nfree=40):
    a, N = fam[0], fam[1]
    lo, hi = window(fam)
    ts = [lo + (hi - lo) * mp.mpf(i) / (ngrid - 1) for i in range(ngrid)]
    grid = [(t, make(fam, t)) for t in ts]
    grid = [(t, P) for t, P in grid if P is not None]
    sigs = [(t, {k: short(v) for k, v in sig(P, cap).items()}) for t, P in grid]
    label = f'{tuple(a)}/{N} base={fam[2]}'
    print(f'{label}: {len(grid)} valid grid t in ({mp.nstr(lo, 5)}, {mp.nstr(hi, 5)})', flush=True)
    hits = []
    for (t0, s0), (t1, s1) in zip(sigs, sigs[1:]):
        for key in s0:
            if key in s1 and s0[key] != s1[key] and 'cap' not in (s0[key][0], s1[key][0]):
                c = max(s0[key][2], s1[key][2]) + 2
                tt = bisect(fam, key, t0, t1, s0[key], c)
                if tt is not None:
                    hits.append((tt, key))
    hits.sort()
    clusters = []
    for tt, key in hits:
        if clusters and tt - clusters[-1][0] < mp.mpf(10) ** -20:   # s930: bisections stop at the ~1e-27 window edge, so 1e-30 split one t* into 2-3
            clusters[-1][1].append(key)
        else:
            clusters.append([tt, [key]])
    out = {'family': label, 'a': a, 'N': N, 'base': [str(x) for x in fam[2]], 'ngrid': ngrid,
           'window': [mp.nstr(lo, 20), mp.nstr(hi, 20)], 'n_tstar': len(clusters), 'flagged': [], 'grid_ctrl': []}
    for t, P in grid[:: max(1, len(grid) // 3)]:
        r = row('grid', P, cap=cap)[0]
        out['grid_ctrl'].append({'t': mp.nstr(t, 20), 'p': r['p'], 'sprime': r['sprime'], 'one_minus_R': r['one_minus_R']})
    nv_hist = collections.Counter()
    for tt, keys in clusters:
        P = make(fam, tt)
        if P is None:
            continue
        prs, kinds = W.pairs(P, NEAR)
        nv = kinds.get('vertex', 0)
        nv_hist[nv] += 1
        if nv < 3:
            continue
        r = row('tstar', P, cap=cap)[0]
        rec = {'t': mp.nstr(tt, 40), 'nprong_changed': len(keys), 'n_vertex_prongs': nv, 'p_near': len(prs),
               'one_minus_R': r['one_minus_R'], 'branch_capped': r['branch_capped']}
        free = (r['one_minus_R'] == r['one_minus_R']) and abs(r['one_minus_R']) > 1e-15
        if free or r['branch_capped']:
            rec.update(classify(P))
            rec['closure'] = closure(P, nfree)
        out['flagged'].append(rec)
        print(f"  t*={mp.nstr(tt, 14)} nv={nv} p={len(prs)} 1-R={r['one_minus_R']:+.3e} capbr={r['branch_capped']}"
              + (f" slits={rec.get('slits')} s'={rec.get('sprime')} closure={rec.get('closure')}" if 'closure' in rec else ''),
              flush=True)
    out['vertex_prong_hist'] = dict(sorted(nv_hist.items()))
    return out


def main(tag, ngrid=60):
    t0 = time.time()
    res = []
    for fam in FAMILY_SETS[tag]:
        res.append(scan(fam, ngrid))
        res[-1]['elapsed_s'] = round(time.time() - t0)
        print(f'  ... {res[-1]["n_tstar"]} t*, {len(res[-1]["flagged"])} flagged, {res[-1]["elapsed_s"]} s', flush=True)
        json.dump(res, open(f'data/s930_free_scan_{tag}.json', 'w'), indent=1)


if __name__ == '__main__':
    main(sys.argv[1], *(int(x) for x in sys.argv[2:3]))
