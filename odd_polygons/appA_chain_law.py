"""[Deposit] Appendix A: the chain law and corners, measured (N ≤ 63). Paper: "Perpendicular billiards in rational polygons of odd denominator". Research name: probes/s940_pin_frontier.py, copied verbatim except imports (deposit names), paths (data/ and figs/ under odd_polygons/) and this paragraph. Run through run_all.py; the usage lines below are the research repo's. Bracketed IDs and section marks (§sNNN) refer to the project's internal ledgers.

s940 `s937-a`: WHAT BOUNDS Π_N's TORI — the frontier of every minimal component of Π_N = (N−4, 4, N−2, 2)/N, side by side,
and the {−1,0,1} relation α(T) it carries (Theorem P's balance, §s927 (3); Lemma S, §s935 (3)).

PRE-REGISTERED (s940, before the first run). Instrument: `engine/poly_seps.py` (exact separatrix trace of the convex polygon,
the component on each side of every saddle connection via tracked germs, ι-pairing, the holonomy balance asserted exactly per
component). Direction ω ⊥ V0V1 (`flow_poly.perp_dirs_poly(P)[0]`, §s935's geometry).
PREDICTION (from [POLY-PI-LOOP] ⊕ §s935 (2)): the 2⌊(N − 7)/8⌋ minimal components come in pairs, one pair per loop ℓ_m; each
component's frontier is ONE member of the loop pair ℓ_m on one side and ONE member each of the chord pairs c_{4m+5}, c_{4m+7} on
the other, so α(T) = ±(ℓ_m − c_{4m+5} − c_{4m+7}) — the proved loop relation is exactly the balance of T_m. The torus is then
a flat torus minus a slit whose two sides are the loop and the two-chord chain (the [NCYL-632] shape, broken at one cone point).
COULD COME OUT OTHERWISE: further connections on a frontier (α with more than three pairs, e.g. at N = 105 where the chords carry
their own 13-term relation, §s935 (3)); both members of one pair on one component's frontier (allowed for minimal components
inside one 𝓕-component); a frontier with the loop on the SAME side as a chord; a capped separatrix on a frontier.
PRIOR ART: §s935 (2) (ladder, rotation 4 sin²(π/N), "pairing presumably the deck reflection in V0V1's line … not checked"),
§s935 (3) (Lemma S: every 𝓕-component carries a nonzero relation; `#𝓕 ≤ 2|S|`), §s937 (the m-th torus for every m), §s938 (3)
([POLY-PI-LOOP]: ℓ_m = c_{4m+5} + c_{4m+7} PROVED); [NCYL-632] (slit tori: one side a single connection, the other broken by cone
points); `rulings.py --grep` *frontier*, *slit* → [NCYL-632], [NCYL-631]; no frontier-by-side measurement of Π_N.
DATA: `data_index.py --grep frontier` → none; `data/s935_tori_63.json` (component itineraries, no sides).
RESULT (s940): the prediction FAILS from N = 23 on (holds at N ≤ 21, one loop): the tori touch each other. CHAIN LAW (read off N ≤ 31,
held on all 20 N = 15 … 41 odd, 47, 55, 63, 71, 79, 105; `data/scratch/s940/chain_law.py`): two chains T_0 – … – T_M; slit m (loop ℓ_m |
chords c_{4m+5}, c_{4m+7}) between T_m and T_{m+1}; α(T_m) = s_{m−1} − s_m, the chain's α = −s_M; masses = §s937's tori to 1e-9; the
N = 105 chord relation in no α. → general_triangles.md §s940 (3).
Output `data/s940_pin_frontier.json` (N ≤ 63 ⊕ 71, 79, 105 merged). USAGE: N N N …
"""
import _setup  # noqa: F401  (deposit: this directory on the path, cwd = odd_polygons/)
import json
import math
import sys

sys.path[:0] = ['.', 'engine', 'probes']
import flow_poly as FP
import poly_seps as PS
from regular_pieces import asym


def classify(R, N, scale):
    c2 = 2 * math.cos(2 * math.pi / N)
    chord = {b: 2 * c2 * math.sin(b * math.pi / N) * scale for b in range(1, N, 2)}
    loop = {m: chord[4 * m + 5] + chord[4 * m + 7] for m in range((N - 7) // 8 + 1) if 4 * m + 7 < N / 2}
    lab = []
    for c in R['conn']:
        l = c['length']
        v0, v1 = c['start'][1], c['arr'][1]
        tag = None
        if v0 == v1:
            hit = [m for m, x in loop.items() if abs(x - l) < 1e-9 * l]
            tag = f'loop{hit[0]}' if len(hit) == 1 else None
        else:
            hit = [b for b, x in chord.items() if abs(x - l) < 1e-9 * l]
            tag = f'c{hit[0]}' if len(hit) == 1 else None
        lab.append(tag or f'V{v0}->V{v1}:{l:.6f}')
    return lab


def run(N):
    P = asym(N)
    j = FP.perp_dirs_poly(P)[0]
    R = PS.seps(P, j)
    L0 = abs(sum(c * complex(math.cos(2 * math.pi * e / P.N), math.sin(2 * math.pi * e / P.N)) for e, c in P.L[0].items()))
    lab = classify(R, N, L0)
    pair_lab = [lab[i] for i, _ in R['pairs']]
    comps = R['comps']
    nonper = [k for k, c in enumerate(comps) if c['kind'] != 'cylinder']
    rows = []
    for k in nonper:
        fr, alpha = PS.balance(R, k)
        rows.append({'comp': k, 'kind': comps[k]['kind'], 'mass': comps[k]['mass'],
                     'frontier': [(lab[i], 'L' if s > 0 else 'R', R['conn'][i]['right' if s > 0 else 'left'],
                                   comps[R['conn'][i]['right' if s > 0 else 'left']]['kind'][:3]) for i, s in fr],
                     'alpha': {(pair_lab[t[1]] if t[0] == 'pair' else f'self:{lab[t[1]]}'): v for t, v in alpha.items()}})
    out = {'balance_asserted_all': all(R['balance_checked'].values()), 'N': N, 'p': len(R['pairs']), 'pairs': sorted(pair_lab), 'n_nonperiodic': len(nonper),
           'expected': 2 * ((N - 7) // 8), 'capped': len(R['capped']), 'rows': rows}
    return out


def verdict(o, N):
    ok = True
    for r in o['rows']:
        a = r['alpha']
        loops = [k for k in a if k.startswith('loop')]
        if len(loops) != 1 or len(a) != 3:
            ok = False
            continue
        m = int(loops[0][4:])
        want = {f'loop{m}': 1, f'c{4 * m + 5}': -1, f'c{4 * m + 7}': -1}
        s = a[loops[0]]
        ok &= all(a.get(k) == s * v for k, v in want.items()) and len(r['frontier']) == 3
    return ok


if __name__ == '__main__':
    res = []
    for N in [int(x) for x in sys.argv[1:]]:
        o = run(N)
        o['prediction_holds'] = verdict(o, N)
        res.append(o)
        print(f"N = {N}: p = {o['p']}, non-periodic {o['n_nonperiodic']} (expected {o['expected']}), capped {o['capped']}, "
              f"prediction {'HOLDS' if o['prediction_holds'] else 'FAILS'}", flush=True)
        for r in o['rows']:
            print('   ', r['comp'], r['kind'][:4], f"mass {r['mass']:.5f}", 'α', r['alpha'], '| frontier', r['frontier'], flush=True)
    json.dump(res, open('data/s940_pin_frontier.json', 'w'), indent=1, default=str)
