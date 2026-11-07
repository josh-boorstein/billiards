"""Shared set-up for the scripts of this directory (the necklace paper).

The separatrix diagram of a perpendicular direction class is walked by `lib/closure.py`'s
exact kite maps; `lib/necklace.py` adds the incidences with Fix(iota) (h, v, f, D, E, the
swapped edges).  The class rows and the cap each needs are the census record (`census_record.json`,
checked by `census.py`), which decides every class
with odd Q <= 39 and even Q <= 40 by the same certificate.  Branch partitions come from the orphan
paper's record (`orphan/partition_record.json`, `lib/partition.py`) through `record` and `centre`;
flow decompositions from `band_record.json` (`lib/flow.py`), read by the scripts that use it.
"""
import importlib.util
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'lib'))

from necklace import readout  # noqa: E402

# the orphan paper's partition record and its per-centre loader (`orphan/common.py`), for the
# scripts that read branch partitions
_spec = importlib.util.spec_from_file_location(
    'orphan_common', os.path.join(HERE, '..', 'orphan', 'common.py'))
_oc = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_oc)
record, centre = _oc.record, _oc.centre

CENSUS = os.path.join(HERE, 'census_record.json')
FULL = '--full' in sys.argv
QUICK_CAP = 3_000_000


def census_rows():
    return json.load(open(CENSUS))['rows']


def row_readout(r, cap=None):
    """The readout of a census row: a CP row is walked at exactly the cap that closes it."""
    if cap is None:
        cap = r['max_steps'] + 1 if r['capped'] == 0 else QUICK_CAP
    return readout(r['P'], r['Q'], r['eps'], step_cap=cap, with_seps=False)


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
