"""[Deposit] Example 5.2 and Proposition 5.3(b): the catalogue of non-periodic triangle classes. Paper: "Perpendicular billiards in rational polygons of odd denominator". Research name: probes/s626_noncp_catalogue.py, copied verbatim except imports (deposit names), paths (data/ and figs/ under odd_polygons/) and this paragraph. Run through run_all.py; the usage lines below are the research repo's. Bracketed IDs and section marks (§sNNN) refer to the project's internal ledgers.

s626 -- s623-a: ONE CATALOGUE of every known Conjecture 2.2(a) counterexample, on one instrument.

WHAT.  The known non-CP triangles in a class of S^1(2d') -- the ten of [NCYL-626]/[NCYL-627] (T0,
flatsurvey's nine, (4,4,7)), the s624-a census's eight, the s625 windows' seven, the s626 line scan's
three -- sit in five data files and a dozen rulings, measured by three generations of the engine.
Here each is re-certified in ALL SIX grid directions (the three side-perpendicular and three
side-parallel classes) on the s626 flow_decomp (SAF gate), and every bad class gets the s624/s625
anatomy (probes/s624_loop_holonomy.analyse, probes/s624_formal_rotation.classify): per minimal
component the t = -X(J) check, t-rank, J-rank, area vs covolume; the loop-rotation lines; the
triangle's shortest relation in the hexagonal norm.  Per triangle also: d, which angle relations
hold with small coefficients (one angle twice another; isosceles; right), and provenance.
The side-parallel non-CP triangles at ODD d ([NCYL-634], Conj 2.2(b)'s excluded set, a different
population) are NOT catalogued here.
USE.  The object for the sort below and for the deferred family paper ([NCYL-626] ↓s623); a
catalogue, not a test.  THE SORT (scored in-thread on the output, against NAMED classes):
(S1) loop line = the triangle's 2:1 relation / a symmetry-family relation (isosceles K = -L, right
K = L) / other; (S2) bad class perpendicular vs parallel, and to which side (opposite the smallest,
middle, largest angle); (S3) minimal components per bad class and their total mass share.
COULD COME OUT OTHERWISE: a re-certification disagreeing with the row's source (an engine
regression, or a source error); a component failing the anatomy (rank, area, one-line loops).
PRIOR ART: [NCYL-626], [NCYL-627], [NCYL-629], [NCYL-631]..[NCYL-637]; data_index --grep
counterexample / nonCP -> no catalogue (other objects).
DATA: data/s623_flow_decomp.json (C), data/s624_line_census.json, data/s625_ward_windows.json,
data/s626_line_scan.json (the sources, read for provenance only); anatomy of 28 rows exists in
s624/s625/s626 outputs from the same probe code and is RECOMPUTED here for uniformity.
Output data/s626_noncp_catalogue.json (resumable, one row per triangle).
"""
import _setup  # noqa: F401  (deposit: this directory on the path, cwd = odd_polygons/)
import json
import os
import sys
import time
from math import gcd

import flow_decomp as fd
import formal_rotation as frot
import loop_holonomy as lh
import symmetric_norm as sn

OUT = 'data/s626_noncp_catalogue.json'
SOURCES = {
    'ten ([NCYL-626]/[NCYL-627])': [(7, 8, 15), (7, 7, 16), (4, 4, 7), (13, 20, 23), (15, 30, 47),
                                    (11, 45, 50), (23, 29, 58), (19, 38, 55), (11, 22, 103), (30, 41, 67)],
    's624 line census': [(13, 26, 119), (15, 30, 131), (17, 34, 145), (17, 34, 149), (17, 34, 169),
                         (31, 74, 99), (43, 86, 89), (49, 67, 98)],
    's625 windows': [(21, 42, 179), (23, 46, 201), (25, 50, 213), (31, 62, 159), (47, 94, 143),
                     (51, 102, 155), (59, 118, 165)],
    's626 line scan': [(17, 34, 217), (45, 90, 119), (51, 71, 102)],
}


def relations(abc):
    a, b, c = abc
    d = a + b + c
    out = []
    for (x, nx), (y, ny) in [((a, 'a'), (b, 'b')), ((a, 'a'), (c, 'c')), ((b, 'b'), (c, 'c'))]:
        if 2 * x == y:
            out.append(f'{ny} = 2{nx}')
        if 2 * y == x:
            out.append(f'{nx} = 2{ny}')
        if x == y:
            out.append(f'{nx} = {ny} (isosceles)')
    if 2 * max(abc) == d:
        out.append('right')
    return out


def anatomy(abc, name, j):
    a, b, c = abc
    d = a + b + c
    from mpmath import mpf
    hol = []
    for h in lh.analyse(abc, j):
        try:                                        # covol is absent when J has no rank-2 lattice
            A_, C_ = mpf(str(h['area_M'])), mpf(str(h['covol_J12']))
            eq = bool(abs(A_ - C_) <= mpf('1e-20') * abs(C_))
        except (ValueError, TypeError, KeyError):
            eq = None
        hol.append({'n': h['n'], 'check': h.get('check_t_eq_minus_XJ'), 't_rank': h.get('t_rank'),
                    'J_rank': h.get('J_rank'), 'area': str(h.get('area_M')),
                    'covol': str(h.get('covol_J12')), 'area_eq_covol': eq})
    lines = []
    for cc in frot.classify(abc, j):
        if cc['kind'].startswith('no_periodic'):
            ds = set()
            for K, L in cc['K']:
                if (K, L) != (0, 0):
                    g = gcd(K, L)
                    k, l = K // g, L // g
                    if k < 0 or (k == 0 and l < 0):
                        k, l = -k, -l
                    ds.add((k, l))
            lines.append(sorted(ds))
    lat = [(K, L) for K in range(-d, d + 1) for L in range(-d, d + 1)
           if (K, L) != (0, 0) and (K * a + L * b) % d == 0]
    m = min(sn.NORMS['H'](v) for v in lat)
    short = sorted({v if v > (0, 0) else (-v[0], -v[1]) for v in lat if sn.NORMS['H'](v) == m})
    return {'hol': hol, 'loop_lines': lines, 'shortest_H': short}


def main():
    rows = json.load(open(OUT)) if os.path.exists(OUT) else {}
    only = [tuple(int(v) for v in x.split(',')) for x in sys.argv[1:]]
    for src, lst in SOURCES.items():
        for abc in lst:
            key = '%d,%d,%d' % abc
            if key in rows or (only and abc not in only):
                continue
            t0 = time.time()
            d = sum(abc)
            dirs = []
            for kind, js in (('perp', fd.perp_dirs(abc)), ('par', fd.side_dirs(abc))):
                for s, j in enumerate(js):
                    r = fd.decompose(abc, j, 20000)
                    bad = [c for c in r['components'] if c['kind'] != 'cylinder']
                    row = {'dir': f'{kind}{s}', 'side': s, 'j': j, 'counts': r['counts'],
                           'bad_mass': sum(c['mass'] for c in bad if c.get('mass') is not None),
                           'max_period': max((c['period'] or 0) for c in r['components']
                                             if c['kind'] == 'cylinder')}
                    if 'no_periodic' in r['counts']:
                        row['anatomy'] = anatomy(abc, row['dir'], j)
                    dirs.append(row)
            rows[key] = {'abc': list(abc), 'd': d, 'source': src, 'relations': relations(abc),
                         'dirs': dirs, 'secs': round(time.time() - t0)}
            print(key, src, rows[key]['relations'],
                  [(x['dir'], x['counts'].get('no_periodic', 0)) for x in dirs if 'no_periodic' in x['counts']],
                  f"{rows[key]['secs']}s", flush=True)
            json.dump(rows, open(OUT, 'w'), indent=1, default=str)


if __name__ == '__main__':
    main()
