"""Run every script of this directory (the beam note: infinitely many perpendicular beams in right triangles at all but countably
many irrational angles) in order, and check each one's verdict line and recorded output against the figures the note prints; exit
nonzero if anything fails.

  python run_all.py            under a minute on one core

The note's theorems are proved by hand. The scripts check its computer-assisted statements: Proposition A.1 (finite partitions with
8 and 10 beams at two irrational angles inside the range of Troubetzkoy's Lemma 12, the beams listed in Appendix A) and the period-944
example of Remark 5.3. The scripts are the research keepers, copied (see README.md); each writes its record to data/.
"""
import json
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(HERE, 'data')


def chk_cert(out):
    with open(os.path.join(D, 'certificates_lemma12.txt')) as f:
        t = f.read()
    return [('Proposition A.1, first angle (sigma^3 - 5 sigma^2 + 5 sigma + 1): 8 beams tiling [0, 1] exactly',
             'x^3 - 5*x^2 + 5*x + 1' in t and 'number of beams: 8;  tiling of [0, 1] checked exactly: True' in t),
            ('Proposition A.1, second angle (sigma^4 - 7 sigma^3 + 15 sigma^2 - 9 sigma - 1): 10 beams tiling [0, 1] exactly',
             'x^4 - 7*x^3 + 15*x^2 - 9*x - 1' in t and 'number of beams: 10;  tiling of [0, 1] checked exactly: True' in t)]


def chk_codes(out):
    return [('Appendix A, first angle: 8 beams; exact chain, fan pairing b_{2k-1} = (b_{2k-2} + b_{2k})/2, printed decimals, '
             '50-digit re-trace of every word', 'triangle 1: 8 beams, chain/fans/floats/trace all OK' in out),
            ('Appendix A, second angle: 10 beams; the same checks', 'triangle 2: 10 beams, chain/fans/floats/trace all OK' in out)]


def chk_c944(out):
    with open(os.path.join(D, 'c944_boundary.json')) as f:
        R = json.load(f)
    lo, hi = (float(v) for v in R['launch'])
    adj = R['adjacency'].get('l1', {})
    a = adj if isinstance(adj, dict) else {}
    pl = a.get('launch', [0, 0])
    return [('Remark 5.3: the non-perpendicular cylinder C of period 944 certified exactly; its strip leaves L_1 between the heights '
             '0.1007262 and 0.1017763 at the angle 0.848062 from the normal',
             'C certified: period 944' in out and round(lo, 7) == 0.1007262 and round(hi, 7) == 0.1017763
             and round(float(R['theta']), 6) == 0.848062),
            ('control: the midpoint of C\'s strip, traced at 50 digits, closes with C\'s word after 944 bounces',
             R['control']['status'] == 'CLOSED' and R['control']['n'] == 944 and R['control']['word_matches']),
            ('Remark 5.3: the adjacent strip, between the heights 0.0982543 and 0.1007262, belongs to a PERPENDICULAR periodic cylinder of '
             'period 7338, certified exactly; it shares C\'s lower boundary line exactly (in Q(2 cos alpha))',
             a.get('period') == 7338 and a.get('ok') and a.get('perpendicular') and a.get('shares_C_lo')
             and round(pl[0], 7) == 0.0982543 and round(pl[1], 7) == 0.1007262)]


STEPS = [('cert_lemma12.py', [], None, chk_cert),
         ('appendix_codes.py', [], 'wrote data/appendix_codes.tex', chk_codes),
         ('c944_boundary.py', [], 'wrote data/c944_boundary.json', chk_c944)]


def main():
    bad = []
    for script, args, verdict, post in STEPS:
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
