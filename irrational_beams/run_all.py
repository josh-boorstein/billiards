"""Run every script of this directory (paper A: finite perpendicular-beam partitions at irrational angles) in order, and check each
one's verdict line and recorded output against the figures the paper prints; exit nonzero if anything fails.

  python run_all.py            the census re-certified on a sample (the record's 5 uncertified points, the 8 rational ones, every 60th
                               other point), the rest of its figures read from data/census_record.json
  python run_all.py --full     every one of the 804 census points re-certified from scratch (about four hours on two cores)
  python run_all.py --sweep    also re-run the census sweep itself first (x in [0.1, 0.9], dx = 2.5e-5, N = 5000; about an hour more),
                               and certify ITS points instead of the record's

The scripts are the research keepers, copied (see README.md); each writes its record to data/ and the figures to figs/.
"""
import json
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(HERE, 'data')


def load(name):
    with open(os.path.join(D, name)) as f:
        return json.load(f)


def grab(out, after, upto):
    """The float printed right after the substring `after` (up to the next `upto`)."""
    i = out.index(after) + len(after)
    return float(out[i:].split(upto)[0].strip())


def chk_firr_even(out):
    return [('Lemma 4.1, even J: 201 admissible patterns, 191 with a constant <= 120, 10 linear in m', '201 structures: 191' in out),
            ('even J: the only vanishing linear form is 4(m - 2) at m = 2', 'vanishing only 4(m−2) at m = 2' in out),
            ('even J: 319 rationals p\'/Q, Q <= 120, excluded, min |E| > 2.5e-5', '319 rationals' in out and grab(out, 'min |E| ≥', '(') > 2.5e-5)]


def chk_firr_odd(out):
    return [('Lemma 4.1, odd J: 201 patterns; linear forms 2(m-1), 4(m-2), 4(m-1), 4m', '201 structures' in out
             and '[(2, -2), (4, -8), (4, -4), (4, 0)]' in out and 'vanishes at m = [2]' in out),
            ('odd J >= 3: 254 rationals excluded, min |E| > 1.7e-4', '254 rationals' in out and grab(out, 'min |E| ≥', '(') > 1.7e-4)]


PAPER = {'n2_4': dict(found=193, certified=193, rational=['2/3', '2/5', '4/5'], distinct_angles=192, uncertified=[]),
         'n6': dict(found=611, certified=609, rational=['2/7', '4/7', '6/7', '8/15'], uncertified=[0.2446654957, 0.6729409541])}


def chk_census(mode):
    def f(out):
        S = load(f'census_summary_{mode}.json')
        res = [(f'{S["recertified"]} points re-certified from scratch, every one agreeing with the record', S['disagreements'] == 0)]
        for k, want in PAPER.items():
            got = S[k]
            row = 'n = 2 and 4' if k == 'n2_4' else 'n = 6'
            res.append((f'Section 5 table, {row}: found {got["found"]}, certified {got["certified"]}, rational x = {", ".join(got["rational"])}',
                        all(got[a] == b for a, b in want.items() if a in ('found', 'certified', 'rational'))))
            res.append((f'{row}: irrational points all units, every branch count even, uncertified {got["uncertified"]}',
                        got['irrational_all_units'] and got['all_even'] and got['uncertified'] == want['uncertified']))
        res.append(('first row: 192 distinct angles (x = 2/5 found by two coincidences)', S['n2_4']['distinct_angles'] == 192))
        res.append((f'minimal polynomials of sigma of degrees {min(S["n2_4"]["deg_min"], S["n6"]["deg_min"])} to '
                    f'{max(S["n2_4"]["deg_max"], S["n6"]["deg_max"])} (paper: 1 to 33)',
                    min(S['n2_4']['deg_min'], S['n6']['deg_min']) == 1 and max(S['n2_4']['deg_max'], S['n6']['deg_max']) == 33))
        return res
    return f


def chk_shapes(out):
    d = load('shapes.json')
    return [('Section 5: of the 192 certified angles with n <= 4, 6 single fans, 184 = (m, m) + one palindromic fan, 2 resolved in part',
             d['angles'] == 192 and d['counts'] == {'single': 6, 'resolved': 184, 'partial': 2})]


def steps(full, sweep):
    mode = 'all' if (full or sweep) else 'sample'
    S = [('firr_J_even.py', [], None, chk_firr_even),
         ('firr_J_odd.py', [], None, chk_firr_odd),
         ('st_count.py', [], 'C5 (control fires): HOLD', None),
         ('prop_w_m3.py', [], 'ALL OK', None)]
    if sweep:
        for lo in ('0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8'):
            hi = str(round(float(lo) + 0.1, 1))
            S.append(('census_sweep.py', ['sweep', lo, hi, '0.000025', '5000', '2'], None, None))
            S.append(('census_retest.py', [f'data/s977_connection_census_{lo}_{hi}.json'], 'finite (C0.test)', None))
    S.append(('census_exact.py', [mode] + (['--sweep'] if sweep else []), 'ALL OK', chk_census(mode)))
    if not sweep:
        S.append(('shapes.py', [mode], 'ALL OK', chk_shapes))
    S.append(('figures.py', [], 'ALL OK', None))
    return S


def main():
    full, sweep = '--full' in sys.argv, '--sweep' in sys.argv
    bad = []
    for script, args, verdict, post in steps(full, sweep):
        name = ' '.join([script] + args)
        print(f'=== {name}', flush=True)
        t0 = time.time()
        p = subprocess.run([sys.executable, os.path.join(HERE, script)] + args, cwd=HERE, capture_output=True, text=True)
        tail = (p.stdout.strip().split('\n') or [''])[-1]
        ok = p.returncode == 0 and (verdict is None or verdict in p.stdout)
        print(f'  [{"ok" if ok else "FAIL"}] exit {p.returncode}; last line: {tail[:160]}  ({time.time() - t0:.0f} s)', flush=True)
        if not ok:
            print(p.stderr[-2000:])
            bad.append(name)
            continue
        for label, c in (post(p.stdout) if post else []):
            print(f'  [{"ok" if c else "FAIL"}] {label}', flush=True)
            if not c:
                bad.append(label)
    print('ALL OK' if not bad else f'FAILED: {bad}')
    sys.exit(1 if bad else 0)


if __name__ == '__main__':
    main()
