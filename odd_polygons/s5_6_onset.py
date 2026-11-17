"""[Deposit] Remark 5.6: the onset census. Paper: "Perpendicular billiards in rational polygons of odd denominator". Research name: probes/s954_onset_converse.py, copied verbatim except imports (deposit names), paths (data/ and figs/ under odd_polygons/) and this paragraph. Run through run_all.py; the usage lines below are the research repo's. Bracketed IDs and section marks (§sNNN) refer to the project's internal ledgers.

s954 keeper (queue `s953-c`, first look): the ONSET CONVERSE of [POLY-S1-X] — is S₁(m; aπ/N) CP in its perpendicular class
for EVERY admissible odd N below the edge, not just the largest one?

OBJECT. S₁(m; x) = (π−(3m+1)x, 4mx, π−(m+1)x, 2x), sides ∝ (2 sinx cos mx, 2 sinx cos 3mx, sin 4mx, sin 4mx); a genuine
quadrilateral iff cos 3mx > 0, i.e. x < π/(6m). At x = aπ/N (gcd(a, N) = 1, N odd) it is (N−3ma−a, 4ma, N−ma−a, 2a)/N =
census quad(N, N−6ma, 4ma, N−2ma, N−2a, 2a, N−2a) (needs N > 6ma). [POLY-S1-X]: non-CP for N > 14ma (the torus). The converse
claim is CP for EVERY admissible N in (6ma, 14ma) — finitely many N per (m, a), but the range reaches x ≈ π/(6m), far from the
torus edge θ = mx = π/14, where nothing about the torus template says anything.

PRE-REGISTERED (s954, before any row ran):
  (C) for every (m, a) with m ≤ 5, a ≤ 3, ma ≤ 9 and every admissible odd N in (6ma, 14ma) with gcd(a, N) = 1:
      the perpendicular class of S₁(m; aπ/N) is CP (flow_poly components all 'cylinder').
  Could it come out otherwise? Yes — any row with a 'no_periodic' / 'undetermined' component refutes (C) as stated; the
  control is that the same instrument returns non-CP with the torus at the first admissible N > 14ma per (m, a) (row 'ctrl').
  If (C) fails, the failing rows are recorded with their rotation numbers, and the honest converse is local (CP just below
  the edge), not global.
PRIOR ART: §s952 (4) Scope (CP below 14ma MEASURED at one N per (m, a), four pairs: N = 55, 83, 83, 167), §s952 (1) A2
  (a = 1…5 at the largest N < 14a, m = 1), §s946 (3) (onset 14m + 1 at a = 1); `rulings.py --grep onset` → nothing on S₁;
  grep `14ma`, `onset converse` → only the s953-c row and §s952.
DATA: `data/s952_s1_family.json` (A2/H′ rows), `data/s946_apex_census.json`; `data_index.py --grep s1_family` → miss.
→ general_triangles.md §s954.   Usage: python probes/s954_onset_converse.py [maxma]   → data/s954_onset_converse.json
"""
import _setup  # noqa: F401  (deposit: this directory on the path, cwd = odd_polygons/)
import json
import math
import sys
import time

import flow_poly as FP
from thmC_m1 import components
from thmC_family import poly

OUT = 'data/s954_onset_converse.json'


def row(N, m, a):
    t = time.time()
    P = poly(N, m, a)
    assert list(P.a) == [N - 3 * m * a - a, 4 * m * a, N - m * a - a, 2 * a]
    j = FP.perp_dirs_poly(P)[0]
    try:
        kinds, rhos = components(P, j)
        err = None
    except Exception as e:                       # noqa: BLE001 — record, do not hide
        kinds, rhos, err = {}, [], repr(e)[:200]
    r = {'N': N, 'm': m, 'a': a, 'x_over_edge': 14 * m * a / N, 'kinds': kinds,
         'rhos': [None if v is None else round(v, 12) for v in rhos],
         'torus_rho': 4 * math.sin(m * a * math.pi / N) ** 2 % 1,
         'nonCP': bool(kinds.get('no_periodic') or kinds.get('undetermined')) or err is not None,
         'err': err, 'secs': round(time.time() - t, 1)}
    print({k: v for k, v in r.items() if k != 'rhos'}, flush=True)
    return r


def main():
    maxma = int(sys.argv[1]) if len(sys.argv) > 1 else 9
    out = {'rows': [], 'ctrl': []}
    for m in range(1, 6):
        for a in range(1, 4):
            if m * a > maxma:
                continue
            for N in range(6 * m * a + 1, 14 * m * a):
                if N % 2 and math.gcd(a, N) == 1:
                    out['rows'].append(row(N, m, a))
            Nc = next(N for N in range(14 * m * a + 1, 14 * m * a + 40) if N % 2 and math.gcd(a, N) == 1)
            out['ctrl'].append(row(Nc, m, a))
            json.dump(out, open(OUT, 'w'), indent=1)
    bad = [(r['N'], r['m'], r['a']) for r in out['rows'] if r['nonCP']]
    out['verdict'] = {'rows': len(out['rows']), 'nonCP_below': bad, 'C_holds': not bad,
                      'ctrl_fires': all(r['nonCP'] for r in out['ctrl'])}
    print('VERDICT', out['verdict'], flush=True)
    json.dump(out, open(OUT, 'w'), indent=1)


if __name__ == '__main__':
    main()
