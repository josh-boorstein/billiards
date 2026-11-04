"""Shared set-up for the scripts carried over from the orphan paper (`neck/from_orphan/`).

These scripts check statements that the s740 re-cut moved from the orphan paper into the necklace
paper.  They read the same partition record and use the same engine as `orphan/`, so this module
loads `orphan/common.py` and re-exports it; `HERE` is this directory.
"""
import importlib.util
import os

HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location(
    'orphan_common', os.path.join(HERE, '..', '..', 'orphan', 'common.py'))
_oc = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_oc)

FULL, RECORD = _oc.FULL, _oc.RECORD
record, centre, check, done = _oc.record, _oc.centre, _oc.check, _oc.done
