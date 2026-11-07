"""Run every statement script of this directory (the method paper, III) in order; exit nonzero if any fails.

  python run_all.py            about ten minutes; Proposition A.23's census recomputed to Q <= 40, read from its record above
  python run_all.py --full     the census recomputed to Q <= 80 (about half an hour more)
"""
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = ['s6_3_arrival.py', 's6_3_boundary.py',
           's6_4_field.py', 's6_4_return_map.py', 's6_4_tiling.py', 's6_4_minimality.py',
           's6_4_cylinders.py', 's6_4_mass.py', 's6_4_words.py',
           's9_2_theorem6.py',
           'sA_12_certificate.py', 'sA_15_finite_check.py', 'sA_23_census.py',
           'sB_4_certificate.py']
TAKES_FULL = {'sA_23_census.py'}


def main():
    full = '--full' in sys.argv
    bad = []
    for s in SCRIPTS:
        print(f'=== {s}', flush=True)
        extra = ['--full'] if full and s in TAKES_FULL else []
        if subprocess.run([sys.executable, os.path.join(HERE, s)] + extra).returncode:
            bad.append(s)
    print('\nFAILED: ' + ', '.join(bad) if bad else f'\nall {len(SCRIPTS)} scripts passed')
    sys.exit(1 if bad else 0)


if __name__ == '__main__':
    main()
