"""Run every statement script of this directory (the method paper, III) in order; exit nonzero if any fails.

  python run_all.py
"""
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = ['s6_3_arrival.py', 's6_3_boundary.py',
           's6_4_field.py', 's6_4_return_map.py', 's6_4_tiling.py', 's6_4_minimality.py',
           's6_4_cylinders.py', 's6_4_mass.py', 's6_4_words.py',
           'sA_12_certificate.py', 'sA_15_finite_check.py', 'sA_23_census.py']


def main():
    bad = []
    for s in SCRIPTS:
        print(f'=== {s}', flush=True)
        if subprocess.run([sys.executable, os.path.join(HERE, s)]).returncode:
            bad.append(s)
    print('\nFAILED: ' + ', '.join(bad) if bad else f'\nall {len(SCRIPTS)} scripts passed')
    sys.exit(1 if bad else 0)


if __name__ == '__main__':
    main()
