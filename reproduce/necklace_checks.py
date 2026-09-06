#!/usr/bin/env python3
"""Reproduce two of the necklace paper's computational checks.

This is a smoke test over the fast end of the necklace strand, not a full
reproduction: it runs the two scripts that verify the counting identity of
Section 3 and the necklace prong census of Sections 2-3, both of which finish
in seconds.  The remaining per-result scripts are slower and are listed
individually in MANIFEST.md, mapped to the statement each one supports.

  s446_covering_identity   Section 3, the counting identity and covering law.
                           Controls: interior prong count == Q - h on 8768/8768
                           class rows (coprime, Q <= 120); e_sigma geometric ==
                           algebraic on 1464/1464 side rows (Q <= 40).

  s449_necklace            Sections 2-3, the necklace itself.  Checks the
                           R-prong count into Delta equals m = (Q-1)/2 on
                           6040/6040 class rows.

Both write their full output to data/.  Everything these scripts report is
computational evidence for the stated ranges; where the paper marks a statement
conditional on (G0b) or on complete periodicity, running these does not
discharge that hypothesis.

Usage:
    python3 reproduce/necklace_checks.py
"""
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

CHECKS = [
    ("s446_covering_identity", "Section 3: counting identity and covering law"),
    ("s449_necklace", "Sections 2-3: the necklace prong census"),
]


def run(name, blurb):
    print(f"\n{'=' * 72}\n{name}  --  {blurb}\n{'=' * 72}")
    env = {
        "PYTHONPATH": f"{ROOT / 'engine'}:{ROOT / 'probes'}",
        "PATH": "/usr/bin:/bin",
    }
    t0 = time.time()
    proc = subprocess.run(
        [sys.executable, str(ROOT / "probes" / f"{name}.py")],
        cwd=ROOT, env=env, capture_output=True, text=True,
    )
    print(proc.stdout, end="")
    if proc.returncode != 0:
        print(proc.stderr, end="")
    print(f"[{name}: exit {proc.returncode}, {time.time() - t0:.1f}s]")
    return proc.returncode == 0


def main():
    (ROOT / "data").mkdir(exist_ok=True)
    ok = all([run(name, blurb) for name, blurb in CHECKS])
    print("\nALL CHECKS RAN CLEANLY" if ok else "\nSOME CHECKS FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
