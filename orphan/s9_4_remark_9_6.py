"""Remark 9.6: at 8/15 the perpendicular direction is not completely periodic -- three cylinders
and a minimal component whose masses sum exactly to Q cot(alpha) -- and n(8/15) = 6 < Q - 1.

That decomposition is the subject of the minimal-component paper and is checked, exactly in
Q(zeta_60), by its scripts; this one runs them rather than repeating them:
  * minimal_component/s3_2_minimality.py   the three-interval exchange and its minimality;
  * minimal_component/s3_3_mass_identity.py  Sigma_cyl + mu_flow(M) = Q cot(alpha), exactly;
  * minimal_component/s4_branch_count.py   the six branches tile L1, n(8/15) = 6.
"""
import os
import subprocess
import sys

from common import HERE, check, done

MC = os.path.join(HERE, '..', 'minimal_component')


def main():
    print('Remark 9.6 -- 8/15, via the minimal-component scripts\n')
    for s in ('s3_2_minimality.py', 's3_3_mass_identity.py', 's4_branch_count.py'):
        r = subprocess.run([sys.executable, os.path.join(MC, s)], capture_output=True, text=True)
        check(f'minimal_component/{s}', r.returncode == 0)
    done()


if __name__ == '__main__':
    main()
