"""Run every statement script of this directory in order; exit nonzero if any fails (about a minute).

  python run_all.py
"""
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = ['s5_1_examples.py', 's5_2_rational_table.py', 's5_3_sqrt41.py', 's5_3_family.py',
           's6_3_tiles.py', 'figures.py']


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
