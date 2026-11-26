"""Imported first by every script of this directory (the beam note): puts this directory and ../lib on the path and makes it
the working directory, so the scripts' `data/...` paths land under `beam_note/`.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
for p in (HERE, os.path.join(HERE, '..', 'lib')):
    if p not in sys.path:
        sys.path.insert(0, p)
os.chdir(HERE)
os.makedirs('data', exist_ok=True)
