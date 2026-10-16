"""Shared set-up for the scripts of this directory (the mirror-direction paper, "Paper A":
*Non-perpendicular periodic orbits in mirror directions of right triangles*).

Every cylinder is decided by `lib/mirror.py`'s exact certificate: the word unfolds to a
translation exactly parallel to a mirror direction, and an open strip of parallel lines
separates every crossed edge strictly -- decided in Q(2 cos a) with certified signs.  How a
word was found is irrelevant to the certificate.  Each script prints one line per check and
exits nonzero if any fails.
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'lib'))


def load(name):
    return json.load(open(os.path.join(HERE, name)))


def check(label, ok):
    print(f'  [{"ok" if ok else "FAIL"}] {label}', flush=True)
    if not ok:
        check.failures += 1
    return ok


check.failures = 0


def done():
    if check.failures:
        print(f'\n{check.failures} CHECK(S) FAILED')
        sys.exit(1)
    print('\nall checks passed')


# The paper's convention: a triangle is the pair {y, -y} (the legs exchanged), listed at y > 0,
# where OB is the SHORTER leg.  A cylinder certified at y < 0 in the class perpendicular to OB
# is, at -y, perpendicular to OA, and vice versa; the two hypotenuse classes are unchanged.
SWAP = {'perp OB': 'perp OA', 'perp OA': 'perp OB', 'perp H': 'perp H', 'par H': 'par H'}


def normalise(y, cls):
    """(|y| as a float key, the class at |y|) for a cylinder certified at y in class cls."""
    return (abs(y), cls if y > 0 else SWAP[cls])
