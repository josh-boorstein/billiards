"""Put engine/ and probes/ on sys.path.

Modules in this deposit import one another by BARE NAME (there is no package and
no __init__.py) — the same layout the research code runs under, preserved here so
that the deposited scripts are byte-identical to the ones that produced the
results rather than a repackaged edit of them.

Import this first from any driver, or set the environment variable directly:

    PYTHONPATH=engine:probes python3 <script>
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
for sub in ("engine", "probes"):
    p = str(ROOT / sub)
    if p not in sys.path:
        sys.path.insert(0, p)
