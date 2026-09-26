"""Shared set-up for the scripts of this directory (the orphan paper).

Every branch count is computed by `lib/partition.py`: the perpendicular partition of a leg at a
fixed centre P/Q, developed exactly in Z[zeta_4Q] with every comparison decided exactly, and
complete when no interval is left open.  The scripts check statements of the paper against it,
one statement or one family per script, and print one line per check.

`partition_record.json` holds the census of every centre the paper's ranges need (n, the orphan
cells and their head-on side, the grazed vertex types in L1 order, the widths), written by the
same engine; a script reads it by default and recomputes it with `--full`.  Small ranges are
always recomputed.
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'lib'))

from partition import partition  # noqa: E402

RECORD = os.path.join(HERE, 'partition_record.json')
FULL = '--full' in sys.argv


def record():
    return json.load(open(RECORD))


def centre(P, Q, fresh=False, beam='1'):
    """One centre's partition summary: from the record, or recomputed (always when `fresh`,
    under --full, or for a beam other than L1)."""
    key = f'{P}/{Q}'
    if not (fresh or FULL or beam != '1'):
        rec = record().get(key)
        if rec is not None:
            return rec
    r = partition(P, Q, beam=beam)
    L = r['leaves']
    out = {'P': P, 'Q': Q, 'open': r['open']}
    if r['open'] == 0:
        out.update({'n': r['n'], 'orphans': r['orphans'],
                    'orphan_head_on': [L[i]['head_on'] for i in r['orphans']],
                    'J_ok': r['J_ok'], 'tiles': r['tiles'],
                    'bounds': ''.join(x['hi_src'][1] for x in L[:-1]),
                    'widths': r['widths'], 'prefix_len': [len(x['word']) for x in L]})
    return out


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
