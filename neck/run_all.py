"""Run every statement script of this directory in order; exit nonzero if any fails.

  python run_all.py            the census rows needing <= 3e6 steps (tens of minutes)
  python run_all.py --full     every census row, each at the cap that closes it (hours)
s10_windows.py is symbolic (exact real-root isolation) and takes about an hour on its own.
"""
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = ['s3_s9_identities.py', 's4_interval_model.py', 's9_2_witness.py', 's10_windows.py',
           'figures.py']


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
