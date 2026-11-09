"""Shared helpers for the scripts of this directory (the polygon note): `lib/` on the path, and the
one-line-per-statement reporting every script uses.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'lib'))


def check(label, ok):
    print(f'  [{"ok" if ok else "FAIL"}] {label}')
    if not ok:
        check.failures += 1
    return ok


check.failures = 0


def done():
    if check.failures:
        print(f'\n{check.failures} CHECK(S) FAILED')
        sys.exit(1)
    print('\nall checks passed')
