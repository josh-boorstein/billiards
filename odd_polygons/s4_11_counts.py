"""[Deposit] Remark 4.11: the free-component counts of the kite (41,15,67,15)π/69 and of (5,6,5,2)π/9. Paper: "Perpendicular billiards in rational polygons of odd denominator". Research name: probes/s958_free_counts.py, copied verbatim except imports (deposit names), paths (data/ and figs/ under odd_polygons/) and this paragraph. Run through run_all.py; the usage lines below are the research repo's. Bracketed IDs and section marks (§sNNN) refer to the project's internal ledgers.

s958 keeper (queue `s955-f`): the two exact counts printed in paper 2's Remark 4.11 (`rem:count`), with a data file.

(1) The non-periodic kite (41,15,67,15)π/69 (Example 5.2): 8 free components, each with exactly 2 frontier connections, p = 19.
(2) (5,6,5,2)π/9 at the shape key [3,6,3,7,2,7] of `s933_diagonal_census.quad`: one ι-invariant free component made of four free
    cylinders.
Both shapes are built by `s933_diagonal_census.quad` (exact sides in ℤ[ζ_{4N}]); the instrument is `engine/poly_seps.seps` (exact saddle
connections on `engine/flow_poly`) and `engine/poly_corners.free_components` (ι on components, free = ι(c) ≠ c, grouped through
connections with free components on both sides). A frontier connection of a free component G is a connection with exactly one side in G.
Previously computed in scratch only (§s945 (5) *Checks*, `data/scratch/s945/zcount.py`, `vprime.py`); this keeper re-derives them.
CONTROL: Π_15 (`s933_regular_pieces.asym(15)`): 2 free components, p = 5 (§s945 (5) *Checks*); the count can come out otherwise there.
PRIOR ART: §s945 (5)–(6) (the counts, scratch); `rulings.py --grep` *free component* → none on these counts. DATA: `data_index.py --grep
free_comp` → none.
Usage: python probes/s958_free_counts.py → data/s958_free_counts.json (seconds).
"""
import _setup  # noqa: F401  (deposit: this directory on the path, cwd = odd_polygons/)
import json

import flow_poly as FP
import poly_corners as PC
import poly_seps as PS
from diagonal_census import quad
from regular_pieces import asym

CASES = [('kite (41,15,67,15)/69', lambda: quad(69, 41, 15, 67, 41, 15, 67)),
         ('(5,6,5,2)/9 key [3,6,3,7,2,7]', lambda: quad(9, 3, 6, 3, 7, 2, 7)),
         ('control Pi_15', lambda: asym(15))]


def count(P):
    j = FP.perp_dirs_poly(P)[0]
    R = PS.seps(P, j, 300000)
    fc, io = PC.free_components(R)
    rows = []
    for G in fc:
        g = set(G)
        fr = [i for i, c in enumerate(R['conn']) if (c['left'] in g) != (c['right'] in g)]
        rows.append({'comps': list(G), 'kinds': [R['comps'][c]['kind'] for c in G],
                     'iota_invariant': {io[c] for c in G} == g, 'frontier': len(fr)})
    return {'angles': list(P.a) if hasattr(P, 'a') else None, 'p': len(R['pairs']), 'n_free_components': len(fc),
            'free_components': rows, 'capped': len(R['capped'])}


def main():
    out = {lab: count(mk()) for lab, mk in CASES}
    k, q, c = (out[lab] for lab, _ in CASES)
    out['verdict'] = {
        'kite: 8 free components, 2 frontier connections each, p = 19':
            k['n_free_components'] == 8 and all(r['frontier'] == 2 for r in k['free_components']) and k['p'] == 19,
        '(5,6,5,2)/9: one iota-invariant free component of four free cylinders':
            q['n_free_components'] == 1 and q['free_components'][0]['iota_invariant']
            and q['free_components'][0]['kinds'].count('cylinder') == 4 == len(q['free_components'][0]['kinds']),
        'control Pi_15: 2 free components, p = 5': c['n_free_components'] == 2 and c['p'] == 5}
    json.dump(out, open('data/s958_free_counts.json', 'w'), indent=1, default=str)
    for kk, v in out['verdict'].items():
        print('OK ' if v else 'FAIL', kk)
    print({lab: (out[lab]['n_free_components'], out[lab]['p'], [r['frontier'] for r in out[lab]['free_components']]) for lab, _ in CASES})


if __name__ == '__main__':
    main()
