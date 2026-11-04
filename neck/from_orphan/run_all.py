"""Run every script of this directory in order; exit nonzero if any fails.

  python run_all.py            the census from orphan/partition_record.json (minutes)
  python run_all.py --full     every count recomputed from scratch (hours)
"""
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = ['named_counts.py', 's5_closed_forms.py', 's9_1_folds.py', 's9_4_remark_9_6.py',
           's9_14_adjacency.py', 'figures.py']


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
