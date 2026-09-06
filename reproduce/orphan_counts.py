#!/usr/bin/env python3
"""Reproduce the certified branch counts n(P/Q) of the orphan paper.

This is the end-to-end check of Appendix A: it exercises all three of the
ingredients described there — the exact partition tree (A.1), the cyclotomic
deduplication that compares boundaries by VALUE rather than by coefficient
vector (A.2), and the completeness certificate that gates the count (A.3).

Two things are checked:

  1. n(2/Q) = Q - 1 on odd Q.  This is Theorem C, which is VERIFIED rather than
     proved: exact on every odd Q <= 45, with no proof of the general case.
     Running it here reproduces the evidence, not a proof.

  2. The plateau trap.  n(3/29) sits at 4 well past depth 10^3 before climbing
     to 19.  Any "stable over two depths" heuristic accepts the wrong answer;
     only the truncation certificate of A.3 is reliable.  Check (2) is the one
     that shows why the certificate is essential rather than cosmetic, and it
     is slow for that reason -- pass --slow to include it.

Usage:
    python3 reproduce/orphan_counts.py [--max-q 21] [--slow]
"""
import argparse
import time

import _path  # noqa: F401  (side effect: puts engine/ and probes/ on sys.path)

from n_exact import n_pq_certified


def theorem_c(max_q):
    """n(2/Q) = Q - 1 on odd Q.  Returns True if every row matched."""
    print(f"Theorem C:  n(2/Q) = Q - 1  on odd Q, 5 <= Q <= {max_q}")
    print(f"  {'Q':>4}  {'n(2/Q)':>7}  {'Q-1':>5}  {'depth':>7}  {'cert':>5}  match")
    ok = True
    for q in range(5, max_q + 1, 2):
        t0 = time.time()
        n, certified, depth = n_pq_certified(2, q)
        match = (n == q - 1) and certified
        ok &= match
        flag = "ok" if match else "MISMATCH"
        print(
            f"  {q:>4}  {n:>7}  {q - 1:>5}  {depth:>7}  {str(certified):>5}  "
            f"{flag}  ({time.time() - t0:.1f}s)"
        )
    return ok


def plateau():
    """n(3/29) = 19, not the 4 that a stability heuristic would report."""
    print("\nThe plateau trap (Appendix A.3):  n(3/29)")
    t0 = time.time()
    n, certified, depth = n_pq_certified(3, 29)
    print(f"  certified n(3/29) = {n} at depth {depth} (certified={certified})")
    print(f"  a 'stable over two depths' test would have accepted 4")
    print(f"  ({time.time() - t0:.1f}s)")
    return n == 19 and certified


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-q", type=int, default=21,
                    help="largest odd Q for the Theorem C sweep (default 21)")
    ap.add_argument("--slow", action="store_true",
                    help="also run the n(3/29) plateau check (minutes)")
    args = ap.parse_args()

    ok = theorem_c(args.max_q)
    if args.slow:
        ok &= plateau()
    print("\nALL CHECKS PASSED" if ok else "\nSOME CHECKS FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
