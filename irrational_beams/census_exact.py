"""[Deposit] Section 5 (the census table) and Appendix B (its certification). Paper: "Finite perpendicular-beam partitions at
irrational angles". New for the deposit: it replaces the research run's three small drivers (exact stage at n <= 4, at n = 6, and a
retry of the root-isolation failures), which looped census_sweep.exact_row over the sweep's output; the certificate itself is
census_sweep.exact_row -> exact_partition (research keepers, copied).

Input: the float-finite census points with at most six beams, either
  data/census_record.json     -- the research run's 804 points (sweep N = 5000, dx = 2.5e-5 on x in [0.1, 0.9]), or
  data/s977_connection_census_*.json -- a fresh sweep (census_sweep.py sweep ... then coincidence_census retest; run_all.py --sweep).
Each point is re-certified from scratch: the coincidence polynomial E(sigma) is rebuilt exactly from the two connection words, factored,
the factor vanishing at the point is isolated, and the beams are traced exactly in Q(sigma); a point is CERTIFIED only if the exact
beam intervals tile (0, 1) with exactly equal shared ends (exact_partition.run_poly).

  python census_exact.py sample        the 5 points the research record leaves uncertified, the 8 rational ones, and every 60th
                                       other point (about 25 points, a few minutes on two cores)
  python census_exact.py all           all 804 points (about four hours on two cores)
  python census_exact.py all --sweep   the points of a fresh sweep instead of the record
Output: data/census_exact_<mode>.json (one row per re-certified point) and data/census_summary_<mode>.json (the paper's figures:
re-certified rows replace the record's, the rest are read from the record in 'sample' mode).
"""
import glob
import json
import sys
from fractions import Fraction
from multiprocessing import Pool

import _setup  # noqa: F401  (deposit: this directory on the path, cwd = irrational_beams/)
import census_sweep as C

KEYS = ('minpoly', 'degree', 'unit', 'irrational', 'CERTIFIED', 'n_exact', 'sign_margin', 'exact')


def load_rows(sweep):
    if not sweep:
        return json.load(open('data/census_record.json'))['rows']
    rows = []
    for p in sorted(glob.glob('data/s977_connection_census_*.json')):
        rows += [r for r in json.load(open(p))['rows'] if r.get('finite') and (r.get('n9') or 99) <= 6]
    return rows


def is_true(v):
    return v is True or v == 'True'


def pick_sample(rows):
    idx = [i for i, r in enumerate(rows) if not is_true(r.get('CERTIFIED'))]          # the record's uncertified rows
    idx += [i for i, r in enumerate(rows) if str(r.get('irrational')) == 'False']      # the rational points
    rest = [i for i in range(len(rows)) if i not in idx]
    idx += rest[::60]
    return sorted(set(idx))


def work(r):
    o = C.exact_row(r)
    return {k: o.get(k) for k in KEYS}


def summary(rows):
    """The figures Section 5 prints, from the rows (each with its exact fields)."""
    out = {}
    for name, ns in (('n2_4', (2, 4)), ('n6', (6,))):
        R = [r for r in rows if int(r.get('n9') or 99) in ns]
        cert = [r for r in R if is_true(r.get('CERTIFIED')) and int(r.get('n_exact') or -1) == int(r['n9'])]
        rat = sorted({str(Fraction(float(r['x'])).limit_denominator(100)) for r in cert if str(r.get('irrational')) == 'False'})
        irr = [r for r in cert if str(r.get('irrational')) != 'False']
        degs = [int(r['degree']) for r in cert]
        out[name] = dict(found=len(R), certified=len(cert), rational=rat, irrational_all_units=all(is_true(r.get('unit')) for r in irr),
                         distinct_angles=len({round(float(r['x']), 9) for r in R}), deg_min=min(degs), deg_max=max(degs),
                         uncertified=sorted(round(float(r['x']), 10) for r in R if r not in cert),
                         all_even=all(int(r['n9']) % 2 == 0 for r in R))
    return out


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else 'sample'
    sweep = '--sweep' in sys.argv
    rows = load_rows(sweep)
    idx = pick_sample(rows) if mode == 'sample' else list(range(len(rows)))
    print(f'{len(rows)} census points with n <= 6; re-certifying {len(idx)} ({mode}{", fresh sweep" if sweep else ""})', flush=True)
    res = {}
    with Pool(2) as P:
        for i, o in zip(idx, P.imap(work, [rows[i] for i in idx], chunksize=1)):
            res[i] = o
            r = rows[i]
            agree = ''
            if is_true(r.get('CERTIFIED')):
                agree = 'agrees with record' if (o['minpoly'] == json.loads(str(r['minpoly'])) and is_true(o['CERTIFIED'])
                                                 and int(o['n_exact'] or -1) == int(r['n_exact'])) else 'DIFFERS FROM RECORD'
            print(f"  x={float(r['x']):.10f} n={r['n9']} deg={o['degree']} unit={o['unit']} CERT={o['CERTIFIED']} n_exact={o['n_exact']} "
                  f"{o['exact'] or ''} {agree}", flush=True)
    merged = []
    for i, r in enumerate(rows):
        m = dict(r)
        if i in res:
            m.update(res[i])
        merged.append(m)
    json.dump([dict(x=rows[i]['x'], n9=rows[i]['n9'], **res[i]) for i in idx], open(f'data/census_exact_{mode}.json', 'w'), default=str)
    S = summary(merged)
    S['recertified'] = len(idx)
    S['disagreements'] = sum(1 for i in idx if is_true(rows[i].get('CERTIFIED')) and not (
        res[i]['minpoly'] == json.loads(str(rows[i]['minpoly'])) and is_true(res[i]['CERTIFIED'])))
    json.dump(S, open(f'data/census_summary_{mode}.json', 'w'), indent=1)
    print(json.dumps(S))
    print('ALL OK' if S['disagreements'] == 0 else 'FAIL: re-certification disagrees with the record')


if __name__ == '__main__':
    main()
