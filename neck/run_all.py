"""Run every statement script of this directory in order; exit nonzero if any fails.

  python run_all.py            records re-verified on a slice; census rows needing <= 3e6 steps
  python run_all.py --full     every record recomputed, every census row at its closing cap (hours)
s10_windows.py and sA3_conditions.py are symbolic (exact real-root isolation) and take about ten
and five minutes.
"""
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = ['s3_s7_identities.py', 's4_interval_model.py', 's5_1_folds.py', 's9_6_adjacency.py',
           's10_2_lattice_counts.py', 's10_windows.py', 'sA3_conditions.py', 's11_cp_band.py',
           'named_counts.py', 'figures.py', 'figures_5_6.py']


def main():
    extra = ['--full'] if '--full' in sys.argv else []
    bad = []
    for s in SCRIPTS:
        print(f'=== {s}', flush=True)
        if subprocess.run([sys.executable, os.path.join(HERE, s)] + extra).returncode:
            bad.append(s)
    print('\nFAILED: ' + ', '.join(bad) if bad else f'\nall {len(SCRIPTS)} scripts passed')
    sys.exit(1 if bad else 0)


if __name__ == '__main__':
    main()
