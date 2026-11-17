"""[Deposit] Example 5.2: the three non-periodic kites (selection; --trace re-traces them). Paper: "Perpendicular billiards in rational polygons of odd denominator". Research name: probes/s929_kite_noncp.py, copied verbatim except imports (deposit names), paths (data/ and figs/ under odd_polygons/) and this paragraph. Run through run_all.py; the usage lines below are the research repo's. Bracketed IDs and section marks (§sNNN) refer to the project's internal ledgers.

s929 keeper (queue `s928-c`): a NON-CP perpendicular class at ODD N in a rational polygon — via the kite transfer.

Transfer (§s928 (4)(b), hand; audit check K `158/158`): the kite K = (a, b/2, c, b/2)/(Ne/2) doubled from T = (a, b, c)/Ne
across the side opposite its even angle b (Ne ≡ 2 mod 4, so N = Ne/2 is odd) has the SAME translation surface as T
(2N kite copies = 2Ne triangle copies, D_Ne = D_N ⊔ D_N μ), and K's perpendicular class is T's PAIR class (the normals
of the two sides at the b-vertex). So T's pair class non-CP (a minimal component) ⟺ K's perpendicular class non-CP.
Selection (complete over `data/s626_noncp_catalogue.json`, the 28 known non-CP triangles, flow_decomp-certified): Ne ≡ 2
(mod 4), and a bad perpendicular direction `perp s` with an EVEN angle at an endpoint vertex of side s (side s = V_s V_{s+1},
the catalogue's convention — fixed by the two controls below: the opposite reading predicts the reverse outcome on both).
Consistency check (not the certificate — the certificate is flow_decomp's on T): the kite's prongs in its class traced
with `s910_poly_slots.Poly.prong`; a minimal component shows as prongs surviving a long cap. CONTROLS: T₀ = (7,8,15)/30
(bad class = third side, not the pair class) and (23,29,58)/110 (bad side 0 not at the even vertex) — their kites must
show NO capped prong.
PRIOR ART: §s928 (4)(b) (the transfer; *whether an odd-N kite can inherit a minimal component needs an even triangle whose
PAIR class is non-CP with Ne ≡ 2 (mod 4)*), queue `s928-c`; `rulings.py --grep` *kite* → [NCYL-200]/[NCYL-282]/… (other
objects), *non-CP polygon* → none.  DATA: `data/s626_noncp_catalogue.json` (the certificates, read here).
AUDIT (s929, `audits/s929_audit_kite.md`): the session's 60000-bounce figure for the N = 69 kite came from an ad-hoc run; the
audit re-traced all three (`data/scratch/s929/audit_kite/`). Run `--trace` to populate `trace`/`controls` in the JSON.
Output: data/s929_kite_noncp.json.  `--trace` re-traces the kites (~15 min, the N = 127, 135 kites dominate); default
lists the selection only.  Run from the repo root with the project PYTHONPATH.
"""
import _setup  # noqa: F401  (deposit: this directory on the path, cwd = odd_polygons/)
import collections
import json
import sys

import mpmath as mp

from poly_slots import Poly

mp.mp.dps = 40


def select():
    cat = json.load(open('data/s626_noncp_catalogue.json'))
    out = []
    for key, v in cat.items():
        abc, Ne = v['abc'], v['d']
        for x in v['dirs']:
            if not x['dir'].startswith('perp') or all(c == 'cylinder' for c in x['counts']):
                continue
            s = x['side']
            for vtx in (s, (s + 1) % 3):
                if Ne % 4 == 2 and abc[vtx] % 2 == 0:
                    out.append({'T': abc, 'Ne': Ne, 'bad_side': s, 'b_vertex': vtx, 'counts': x['counts'],
                                'bad_mass': x['bad_mass']})
    # one row per (T, b-vertex): both pair-class sides are bad together (same class)
    rows = {}
    for r in out:
        k = (tuple(r['T']), r['b_vertex'])
        rows.setdefault(k, dict(r, bad_sides=[]))['bad_sides'].append(r['bad_side'])
    return list(rows.values())


def kite_of(abc, Ne, vtx):
    b = abc[vtx]
    a, c = abc[(vtx + 1) % 3], abc[(vtx - 1) % 3]
    N = Ne // 2
    ka = [a, b // 2, c, b // 2]
    sa, sc = mp.sin(mp.pi * a / Ne), mp.sin(mp.pi * c / Ne)
    return ka, N, Poly(ka, N, lengths=[sc, sa, sa, sc])


def trace(K, caps=(3000, 30000)):
    C = K.classes()
    assert len(C) == 1
    res = {k: K.prong(*k, caps[0])[0] for k in K.slots(C[0])}
    capped = [k for k, r in res.items() if r == 'cap']
    esc = {k: K.prong(*k, caps[1])[0] for k in capped}
    return {'slots': len(res), 'counts': dict(collections.Counter(res.values())), 'capped_at': caps,
            'still_capped': sorted(k for k, r in esc.items() if r == 'cap'),
            'escalated': dict(collections.Counter(esc.values()))}


CONTROLS = [((7, 8, 15), 30, 1), ((23, 29, 58), 110, 2)]   # (T, Ne, even vertex): bad class is NOT the pair class


def main():
    do_trace = '--trace' in sys.argv
    rows = select()
    for r in rows:
        ka, N, K = kite_of(r['T'], r['Ne'], r['b_vertex'])
        r.update(kite=ka, N=N, convex=all(x < N for x in ka))
        if do_trace:
            r['trace'] = trace(K)
        print(f"T={tuple(r['T'])}/{r['Ne']} b-vertex {r['b_vertex']} bad sides {r['bad_sides']} mass {r['bad_mass']:.4f} "
              f"-> kite {tuple(ka)}/{N} convex={r['convex']} {r.get('trace', '')}", flush=True)
    ctrl = []
    if do_trace:
        for abc, Ne, vtx in CONTROLS:
            ka, N, K = kite_of(list(abc), Ne, vtx)
            t = trace(K)
            ctrl.append({'T': abc, 'Ne': Ne, 'kite': ka, 'N': N, 'trace': t})
            print(f'control T={abc}/{Ne} -> kite {tuple(ka)}/{N}: {t}', flush=True)
    json.dump({'selection': rows, 'controls': ctrl}, open('data/s929_kite_noncp.json', 'w'), indent=1)


if __name__ == '__main__':
    main()
