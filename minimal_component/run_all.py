"""Run every statement script of this directory in order; exit nonzero if any fails.

  python run_all.py            the census from its record (seconds)
  python run_all.py --full     the census recomputed for odd Q <= 31 (about an hour)
"""
import subprocess
import sys
import os

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = ['s2_phase_space.py', 's3_1_return_map.py', 's3_1_boundary5.py',
           's3_2_minimality.py', 's3_3_mass_identity.py', 's3_4_completeness.py',
           's4_branch_count.py', 's5_coverage.py', 's5_census.py', 's6_connection.py',
           'appendix_words.py']


def main():
    full = '--full' in sys.argv
    bad = []
    for s in SCRIPTS:
        args = [sys.executable, os.path.join(HERE, s)]
        if s == 's5_census.py' and not full:
            args.append('--record')
        print(f'=== {s}', flush=True)
        r = subprocess.run(args)
        if r.returncode:
            bad.append(s)
    print('\nFAILED: ' + ', '.join(bad) if bad else f'\nall {len(SCRIPTS)} scripts passed')
    sys.exit(1 if bad else 0)


if __name__ == '__main__':
    main()
