"""Run every statement script of this directory (the polygon note) in order; exit nonzero if any fails.

  python run_all.py            about twenty-five minutes (the slot census nine, the odd-N census six, the
                             isosceles halves five)
"""
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = ['s4_s5_slots_and_bands.py', 's6_isosceles_beams.py', 's8_2_containment.py', 's8_families.py',
           's8_5_isosceles_odd.py', 's8_6_isosceles_even.py', 's9_odd_census.py']


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
