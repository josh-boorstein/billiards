#!/usr/bin/env python3
"""
right_triangle_billiards.py

Reusable numerical tools for studying billiards in a right triangle.

Geometry
--------
O = (0,0) is the angle alpha opposite leg 1.
R = (cot(alpha),0) is the right-angle vertex.
A = (cot(alpha),1) is the top of leg 1.
Leg 1 is the vertical segment R--A and has length 1.

The coordinate s is height on leg 1:
    start point = (cot(alpha), s), 0 < s < 1.

A perpendicular launch from leg 1 into the triangle has global velocity (-1,0),
i.e. outgoing_angle = pi.
"""

from __future__ import annotations

from dataclasses import dataclass
import argparse
import csv
import math
from typing import List, Optional, Sequence, Tuple


EPS = 1e-11


@dataclass(frozen=True)
class Hit:
    index: int
    side: str
    x: float
    y: float
    incoming_angle: float


@dataclass(frozen=True)
class L1Return:
    index: int
    hit_index: int
    s: float
    outgoing_angle: float
    word_since_last_l1: Tuple[str, ...]


@dataclass(frozen=True)
class Trace:
    alpha: float
    start_s: float
    hits: Tuple[Hit, ...]
    l1_returns: Tuple[L1Return, ...]
    full_word: Tuple[str, ...]
    returned_perpendicular: bool


@dataclass(frozen=True)
class Branch:
    lo: float
    hi: float
    signature: Tuple[str, ...]
    sample_s: float
    l1_maps: Tuple[Tuple[float, float], ...]
    max_fit_error: float


def _cross(ax: float, ay: float, bx: float, by: float) -> float:
    return ax * by - ay * bx


def _normalize(x: float, y: float) -> Tuple[float, float]:
    n = math.hypot(x, y)
    if n == 0:
        raise ValueError("zero vector")
    return x / n, y / n


def _angle(x: float, y: float) -> float:
    return math.atan2(y, x)


def _reflect_across_unit_axis(
    vx: float, vy: float, ex: float, ey: float
) -> Tuple[float, float]:
    dot = vx * ex + vy * ey
    return 2.0 * dot * ex - vx, 2.0 * dot * ey - vy


class RightTriangleBilliard:
    def __init__(self, alpha: float):
        if not (0.0 < alpha < math.pi / 2):
            raise ValueError("alpha must be between 0 and pi/2")
        self.alpha = float(alpha)
        self.leg2 = 1.0 / math.tan(alpha)
        self.hyp = 1.0 / math.sin(alpha)

        self.O = (0.0, 0.0)
        self.R = (self.leg2, 0.0)
        self.A = (self.leg2, 1.0)
        self.axis_H = _normalize(self.A[0], self.A[1])

    def start_point(self, s: float) -> Tuple[float, float]:
        return self.leg2, float(s)

    def _candidate_times(
        self, px: float, py: float, vx: float, vy: float
    ) -> List[Tuple[float, str]]:
        out: List[Tuple[float, str]] = []
        L = self.leg2

        # L1: x = L
        if vx > EPS:
            t = (L - px) / vx
            y = py + t * vy
            if t > EPS and -1e-9 <= y <= 1.0 + 1e-9:
                out.append((t, "L1"))

        # L2: y = 0
        if vy < -EPS:
            t = -py / vy
            x = px + t * vx
            if t > EPS and -1e-9 <= x <= L + 1e-9:
                out.append((t, "L2"))

        # H: line through O and A
        Ax, Ay = self.A
        denom = _cross(vx, vy, Ax, Ay)
        numer = -_cross(px, py, Ax, Ay)
        if abs(denom) > EPS:
            t = numer / denom
            if t > EPS:
                x = px + t * vx
                y = py + t * vy
                param = (x * Ax + y * Ay) / (Ax * Ax + Ay * Ay)
                if -1e-9 <= param <= 1.0 + 1e-9:
                    out.append((t, "H"))

        return out

    def reflect_velocity(self, side: str, vx: float, vy: float) -> Tuple[float, float]:
        if side == "L1":
            return -vx, vy
        if side == "L2":
            return vx, -vy
        if side == "H":
            vx2, vy2 = _reflect_across_unit_axis(vx, vy, *self.axis_H)
            return _normalize(vx2, vy2)
        raise ValueError(f"unknown side {side!r}")

    def trace(
        self,
        s: float,
        outgoing_angle: float = math.pi,
        max_hits: int = 5000,
        max_l1_returns: Optional[int] = None,
        stop_at_perpendicular_return: bool = True,
        perpendicular_tol: float = 1e-8,
    ) -> Trace:
        if not (0.0 < s < 1.0):
            raise ValueError("s must be in (0,1)")

        px, py = self.start_point(s)
        vx, vy = math.cos(outgoing_angle), math.sin(outgoing_angle)

        hits: List[Hit] = []
        returns: List[L1Return] = []
        word_since_l1: List[str] = []

        for hit_index in range(1, max_hits + 1):
            cand = self._candidate_times(px, py, vx, vy)
            if not cand:
                break

            t, side = min(cand, key=lambda z: z[0])
            px += t * vx
            py += t * vy

            if abs(px - self.leg2) < 1e-9:
                px = self.leg2
            if abs(py) < 1e-9:
                py = 0.0
            if abs(py - 1.0) < 1e-9:
                py = 1.0

            hits.append(Hit(hit_index, side, px, py, _angle(vx, vy)))
            word_since_l1.append(side)

            vx_ref, vy_ref = self.reflect_velocity(side, vx, vy)

            if side == "L1":
                returns.append(
                    L1Return(
                        index=len(returns) + 1,
                        hit_index=hit_index,
                        s=py,
                        outgoing_angle=_angle(vx_ref, vy_ref),
                        word_since_last_l1=tuple(word_since_l1),
                    )
                )
                word_since_l1 = []

                if stop_at_perpendicular_return:
                    if abs(vx_ref + 1.0) < perpendicular_tol and abs(vy_ref) < perpendicular_tol:
                        return Trace(self.alpha, s, tuple(hits), tuple(returns), tuple(h.side for h in hits), True)

                if max_l1_returns is not None and len(returns) >= max_l1_returns:
                    return Trace(self.alpha, s, tuple(hits), tuple(returns), tuple(h.side for h in hits), False)

            vx, vy = vx_ref, vy_ref

        return Trace(self.alpha, s, tuple(hits), tuple(returns), tuple(h.side for h in hits), False)

    def branch_signature(self, s: float, max_hits: int = 5000) -> Tuple[str, ...]:
        return self.trace(s, max_hits=max_hits, stop_at_perpendicular_return=True).full_word

    def discover_periodic_branches(
        self,
        n_samples: int = 1000,
        max_hits: int = 5000,
        refine_steps: int = 35,
    ) -> List[Tuple[float, float, Tuple[str, ...]]]:
        eps = 1e-7
        ss = [eps + (1.0 - 2.0 * eps) * i / (n_samples - 1) for i in range(n_samples)]
        sigs = [self.branch_signature(s, max_hits=max_hits) for s in ss]

        boundaries: List[float] = [0.0]
        for i in range(n_samples - 1):
            if sigs[i] != sigs[i + 1]:
                lo, hi = ss[i], ss[i + 1]
                sig_lo = sigs[i]
                for _ in range(refine_steps):
                    mid = 0.5 * (lo + hi)
                    if self.branch_signature(mid, max_hits=max_hits) == sig_lo:
                        lo = mid
                    else:
                        hi = mid
                b = 0.5 * (lo + hi)
                if b - boundaries[-1] > 1e-8:
                    boundaries.append(b)
        boundaries.append(1.0)

        intervals = []
        for lo, hi in zip(boundaries[:-1], boundaries[1:]):
            if hi - lo > 1e-10:
                mid = 0.5 * (lo + hi)
                intervals.append((lo, hi, self.branch_signature(mid, max_hits=max_hits)))
        return intervals

    def fit_l1_return_maps(
        self,
        lo: float,
        hi: float,
        max_hits: int = 5000,
    ) -> Branch:
        s1 = lo + 0.25 * (hi - lo)
        s2 = lo + 0.50 * (hi - lo)
        s3 = lo + 0.75 * (hi - lo)

        traces = [
            self.trace(s1, max_hits=max_hits, stop_at_perpendicular_return=True),
            self.trace(s2, max_hits=max_hits, stop_at_perpendicular_return=True),
            self.trace(s3, max_hits=max_hits, stop_at_perpendicular_return=True),
        ]

        n_returns = min(len(t.l1_returns) for t in traces)
        maps: List[Tuple[float, float]] = []
        max_err = 0.0

        for j in range(n_returns):
            y1 = traces[0].l1_returns[j].s
            y2 = traces[1].l1_returns[j].s
            y3 = traces[2].l1_returns[j].s
            a = (y2 - y1) / (s2 - s1)
            b = y1 - a * s1
            max_err = max(max_err, abs((a * s3 + b) - y3))
            maps.append((a, b))

        return Branch(
            lo=lo,
            hi=hi,
            signature=traces[1].full_word,
            sample_s=s2,
            l1_maps=tuple(maps),
            max_fit_error=max_err,
        )


def format_word(word: Sequence[str]) -> str:
    return " ".join(word)


def parse_alpha(expr: str) -> float:
    return float(eval(expr, {"pi": math.pi, "math": math, "__builtins__": {}}))


def demo(alpha_expr: str, n_samples: int, out_csv: Optional[str]) -> None:
    alpha = parse_alpha(alpha_expr)
    billiard = RightTriangleBilliard(alpha)

    intervals = billiard.discover_periodic_branches(n_samples=n_samples)
    branches = [billiard.fit_l1_return_maps(lo, hi) for lo, hi, _ in intervals]

    print(f"alpha = {alpha:.15g}")
    print(f"leg2 = cot(alpha) = {billiard.leg2:.15g}")
    print(f"found {len(branches)} intervals")
    print()

    for i, br in enumerate(branches, start=1):
        print(f"interval {i}: ({br.lo:.12f}, {br.hi:.12f})")
        print(f"  word: {format_word(br.signature)}")
        print("  L1 return maps s_j = a*s + b:")
        for j, (a, b) in enumerate(br.l1_maps, start=1):
            print(f"    {j}: {a:+.12f} * s {b:+.12f}")
        print(f"  max affine fit error: {br.max_fit_error:.3e}")
        print()

    if out_csv:
        with open(out_csv, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["branch", "lo", "hi", "word", "return_index", "a", "b", "max_fit_error"])
            for i, br in enumerate(branches, start=1):
                for j, (a, b) in enumerate(br.l1_maps, start=1):
                    w.writerow([i, br.lo, br.hi, format_word(br.signature), j, a, b, br.max_fit_error])
        print(f"wrote {out_csv}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--alpha", default="pi/7", help="e.g. 'pi/7', '2*pi/11', or '0.449'")
    parser.add_argument("--samples", type=int, default=1000)
    parser.add_argument("--csv", default=None)
    args = parser.parse_args()
    demo(args.alpha, args.samples, args.csv)


if __name__ == "__main__":
    main()
