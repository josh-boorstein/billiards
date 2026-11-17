"""Imported first by every script of this directory (paper 2): puts this directory on the path and makes it the working
directory, so the scripts' `data/...` and `figs/...` paths land under `odd_polygons/`. The directory is self-contained: it uses
nothing from `lib/`.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
os.chdir(HERE)
os.makedirs('data', exist_ok=True)
os.makedirs('figs', exist_ok=True)
