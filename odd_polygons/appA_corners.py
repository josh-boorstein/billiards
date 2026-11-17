"""[Deposit] Appendix A, measured: for odd 15 ≤ N ≤ 63 each free component has one boundary cycle, three corners and genus 2M + 2. Paper: "Perpendicular billiards in rational polygons of odd denominator". Research name: probes/s945_corner_census.py, copied verbatim except imports (deposit names), paths (data/ and figs/ under odd_polygons/) and this paragraph. Run through run_all.py; the usage lines below are the research repo's. Bracketed IDs and section marks (§sNNN) refer to the project's internal ledgers.

s945 keeper (queue `s944-c` ⊕ `s944-d`): frontier corner angles k_q of every 𝓕-component — Π_N = (N−4, 4, N−2, 2)/N for odd
15 ≤ N ≤ 63, and the 19 non-CP rows of the s933 diagonal census (odd N ≤ 17, all generic: every gcd(aᵢ, N) = 1).

NOT PRE-REGISTERED as a scorer; the two questions were fixed before the runs by the queue rows: `s944-c` — Lemma GB checked
directly on a real minimal component, and where Π_N's excess corners sit (§s944 (8) expected Σ(k − 1) = 2(2g − 2 + b) at the
two large-angle vertices); `s944-d` — the candidate `k = 3` input *excess corners need obtuse vertices* (a triangle has ≤ 1).
COULD COME OUT OTHERWISE: a cylinder with a corner k ≥ 2 or χ ≠ 0 (tool error); GB's genus ≠ the IET genus after cutting
interior connections; excess at V1/V3 of Π_N; a census row with excess only at acute vertices, or at ONE vertex point.
INSTRUMENT: `engine/poly_corners.py` (gated; orientation control 0 vs 44–156 fails per Π_N).
PRIOR ART: §s944 (7)–(8) (Lemma GB, Corollary, the Π_N control *not computed*); §s941 ([POLY-PI-FRONTIER]: T_m's frontier,
c₁, c₃ interior to T_0); §s940 (3) ([POLY-PI-CHAIN]); §s933 (4) ⊕ `data/s933_diagonal_census.json` (the 19 rows);
`rulings.py --grep` *corner* → other objects (WFLOOR/NB), *obtuse* → [NCYL-725] (orthic map, unrelated); no corner angle of
any minimal component on record.
DATA: `data_index.py --grep corner` → none; `data/s933_diagonal_census.json` (keys, counts — no corners).
RESULT → `general_triangles.md` §s945; [NCYL-732] (the obtuse input refuted).
Output `data/s945_corner_angles.json`.  Runtime ≈ 1 min.
"""
import _setup  # noqa: F401  (deposit: this directory on the path, cwd = odd_polygons/)
import json
import sys
import time

sys.path[:0] = ['.', 'engine', 'probes']
import flow_poly as FP
import poly_corners as PC
import poly_seps as PS
from diagonal_census import quad
from regular_pieces import asym


def run(P, cap=30000):
    j = FP.perp_dirs_poly(P)[0]
    R = PS.seps(P, j, cap)
    while any((c['left'] != c['right']) for c in R['capped']) and cap < 10 ** 7:   # a capped ray on SOME frontier: escalate
        cap *= 10
        R = PS.seps(P, j, cap)
    pts, um, fo, fr = PC.place(P, j, R)
    assert fo == 0 and fr > 0
    cyl_bad = stacked = 0
    for k, c in enumerate(R['comps']):
        if c['kind'] == 'cylinder':
            A = PC.analyse(P, j, R, pts, [k])
            if not (A['b'] == 2 and A['chi'] == 0 and all(x['k'] == 1 for x in A['corners'])):
                B = PC.analyse(P, j, R, pts, [k], cut_interior=True)
                if B['chi'] == 0 and all(x['k'] == 1 for x in B['corners']):
                    stacked += 1
                else:
                    cyl_bad += 1
    per = []
    for k, c in enumerate(R['comps']):
        if c['kind'] != 'cylinder':
            A = PC.analyse(P, j, R, pts, [k])
            B = PC.analyse(P, j, R, pts, [k], cut_interior=True)
            per.append({'comp': k, 'mass': c['mass'], 'chi': A['chi'], 'b': A['b'], 'g': A['genus'], 'chi_cut': B['chi'],
                        'n_interior_conn': sum(1 for x in R['conn'] if x['left'] == k == x['right']),
                        'corners': sorted((x['v'], x['k']) for x in A['corners'])})
    fc, _ = PC.free_components(R)
    F = []
    for g in fc:
        A = PC.analyse(P, j, R, pts, g)
        cr = sorted((x['v'], x['k']) for x in A['corners'])
        exc = [(v, k) for v, k in cr if k >= 2]
        F.append({'comps': g, 'kinds': sorted(R['comps'][c]['kind'][:3] for c in g), 'chi': A['chi'], 'b': A['b'],
                  'g': A['genus'], 'interior_pts': A['interior'], 'corners': cr,
                  'excess_vertices': sorted({v for v, _ in exc}),
                  'excess_at_acute': sorted({v for v, _ in exc if 2 * P.a[v] < P.d}),
                  'k_le_half_a': all(k <= P.a[v] // 2 for v, k in cr)})
    return {'cap': cap, 'a': P.a, 'N': P.d, 'orient_fails': fo, 'reverse_fails': fr, 'cyl_anomalies': cyl_bad, 'stacked_cyl': stacked,
            'minimal': per, 'F': F}


def main():
    t = time.time()
    out = {'pi': [], 'census': []}
    for N in range(15, 64, 2):
        o = run(asym(N))
        M = len(o['minimal']) // 2 - 1
        pred = sorted([(0, 2 * M + 2), (0, 2 * M + 2), (2, 4 * M + 5)])
        o['M'], o['law_ok'] = M, all(f['corners'] == pred and f['b'] == 1 and f['g'] == 2 * M + 2 for f in o['F'])
        o['torus_ok'] = all(x['chi_cut'] == -1 for x in o['minimal'])
        out['pi'].append(o)
        print('Pi', N, 'M', M, 'law', o['law_ok'], 'cut tori χ=−1', o['torus_ok'], 'cyl anomalies', o['cyl_anomalies'],
              'stacked', o['stacked_cyl'], f'{time.time() - t:.0f}s', flush=True)
    d = json.load(open('data/s933_diagonal_census.json'))
    for N in ['9', '13', '15', '17']:
        for r in d[N]['nonCP_rows']:
            o = run(quad(int(N), *r['key']))
            o['key'] = r['key']
            out['census'].append(o)
            print('census', o['a'], N, [(f['corners'], f['excess_at_acute']) for f in o['F'] if f['excess_vertices']][:1],
                  'cyl anomalies', o['cyl_anomalies'], flush=True)
    c = out['census']
    out['summary'] = {
        'pi_law_all': all(o['law_ok'] for o in out['pi']), 'pi_tori_all': all(o['torus_ok'] for o in out['pi']),
        'pi_N': [o['N'] for o in out['pi']],
        'census_rows': len(c),
        'rows_excess_at_acute': [(o['a'], o['N']) for o in c if any(f['excess_at_acute'] for f in o['F'])],
        'rows_single_excess_vertex': [(o['a'], o['N']) for o in c if any(len(f['excess_vertices']) == 1 for f in o['F'])],
        'rows_k_gt_half_a': [(o['a'], o['N']) for o in c + out['pi'] if not all(f['k_le_half_a'] for f in o['F'])],
        'cyl_anomalies': sum(o['cyl_anomalies'] for o in c + out['pi']),
        'secs': round(time.time() - t)}
    print(out['summary'])
    json.dump(out, open('data/s945_corner_angles.json', 'w'), indent=1)


if __name__ == '__main__':
    main()
