"""Imported first by every script of this directory (paper 1): puts `lib/` and this directory on the path and makes this
directory the working directory, so the scripts' `data/...` and `figs/...` paths land under `triangles/`.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
for p in (HERE, os.path.join(HERE, '..', 'lib')):
    if p not in sys.path:
        sys.path.insert(0, p)
os.chdir(HERE)
os.makedirs('data', exist_ok=True)
os.makedirs('figs', exist_ok=True)
