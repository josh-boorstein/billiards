"""Run every statement script of this directory in order; exit nonzero if any fails.

  python run_all.py            the census from partition_record.json (minutes)
  python run_all.py --full     every count recomputed from scratch (hours)
"""
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = ['s1_named_counts.py', 's4_theorem_a.py', 's5_3_three_sided.py',
           's6_7_surface.py', 'sA_passage.py', 'figures.py']


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
