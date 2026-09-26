#!/usr/bin/env python3
"""trans_surface.py -- translation surfaces from POLYGONS + TRANSLATION GLUINGS, and the
CYLINDER DECOMPOSITION of a coordinate direction by SEPARATRIX CUT.

⇒⇒ WHY THIS EXISTS.  `engine/cyl_diagram.py` builds a diagram from an ALREADY-KNOWN cycle
structure and does not discover one; `probes/s315_cylinder_count.py` discovers cylinders by
tracing the BILLIARD.  Neither takes a polygon model and returns its cylinders.
`probes/s520_polygon_staircase.py` did that for the regular `2Q`-gon, but only under the
assumption that a strip between consecutive VERTEX levels is a union of closed leaves --
which is false in general (it raises `StripsNotCylinders`), and is exactly what blocked the
`P in {2, Q-2}` arm ([NCYL-425] (9)).

⇒ THE FIX, AND IT IS THE ONLY IDEA HERE.  A strip fails to be a union of closed leaves when a
SEPARATRIX crosses its interior.  So mark not just the vertex levels but their CLOSURE UNDER
THE FLOW: push each marked level forward and backward through the gluing until it dies on a
cone point.  Between consecutive marked levels no separatrix can enter, so every sub-strip IS
a union of closed leaves, the first-return map permutes sub-strips of equal width, and the
cylinders are its cycles.  On a surface whose direction is periodic the closure is finite.

CONVENTIONS.
  * A `Surface` is a list of convex polygons (vertices CCW) plus a gluing
    `glue[(p, r)] = (p', r')`: side `r` of polygon `p` (from vertex `r` to `r+1`) is
    identified with side `r'` of polygon `p'` by a TRANSLATION, orientation-reversing on the
    boundary (`[a,b] -> [d,c]`).  The gluing must be an involution.
  * `axis = 0` means STRIPS ARE LEVEL SETS OF `x` and the flow is along `+y` (a "vertical"
    direction); `axis = 1` is the transpose.  Both directions must be non-parallel to some
    side for the strips to be two-dimensional; sides PARALLEL to the flow are fine (they are
    saddle connections and their level is always marked).

⚠ FLOAT.  Everything is `float` with an explicit `tol`.  Nothing here is ring-certified.  The
one structural assertion -- that a cycle's sub-strips all have the same width -- is CHECKED
and raises, so a silent wrong answer needs the check itself to be wrong.

TRAPS (measured, in this order).
  (i)  Read the crossing count as an AREA RATIO `area(H_k n V_t)/(h_k v_t)`, or as an
       arclength ratio; NEVER as a run count -- adjacent components merge ([OPS-290] (1)).
       `crossing_boxes` and `crossing_walk` are the two independent routes and disagreeing is
       a real failure, not a tolerance issue.
  (ii) The marked-level closure needs BOTH flow directions.  Forward only loses the levels
       whose backward separatrix dies on a different cone point.
  (iii) A level whose exit point IS a vertex terminates -- that is a saddle connection
       arriving at a cone point, not a level to propagate.
"""
from __future__ import annotations

import math
from typing import Dict, List, Sequence, Tuple

Pt = Tuple[float, float]

TOL = 1e-9


# --------------------------------------------------------------------------- surface

class Surface:
    """Convex polygons glued in pairs by translation along sides."""

    def __init__(self, polys: Sequence[Sequence[Pt]],
                 glue: Dict[Tuple[int, int], Tuple[int, int]]):
        self.polys = [list(p) for p in polys]
        self.glue = dict(glue)
        for (p, r), (q, s) in self.glue.items():
            assert self.glue[(q, s)] == (p, r), "gluing is not an involution"

    # -- basics ----------------------------------------------------------
    def nsides(self, p: int) -> int:
        return len(self.polys[p])

    def side(self, p: int, r: int) -> Tuple[Pt, Pt]:
        P = self.polys[p]
        n = len(P)
        return P[r % n], P[(r + 1) % n]

    def shift(self, p: int, r: int) -> Pt:
        """Translation carrying side `(p,r)` onto its partner: `[a,b] -> [d,c]`."""
        q, s = self.glue[(p, r % self.nsides(p))]
        a, _b = self.side(p, r)
        _c, d = self.side(q, s)
        return (d[0] - a[0], d[1] - a[1])

    def area(self) -> float:
        tot = 0.0
        for P in self.polys:
            s = 0.0
            for i in range(len(P)):
                a, b = P[i], P[(i + 1) % len(P)]
                s += a[0] * b[1] - b[0] * a[1]
            tot += abs(s) / 2
        return tot

    # -- cross-section ---------------------------------------------------
    def chord(self, p: int, axis: int, c: float, tol: float = TOL):
        """The chord of polygon `p` on the line {coord[axis] == c}.

        Returns `(lo, hi, side_lo, side_hi, vert_lo, vert_hi)` where `lo`/`hi` are the
        extreme values of the OTHER coordinate, `side_*` the side index the endpoint lies on
        (`None` if the endpoint is a vertex, i.e. the separatrix dies there), and `vert_*`
        whether the endpoint is a vertex.  `None` if the line misses the polygon."""
        P = self.polys[p]
        n = len(P)
        o = 1 - axis
        hits = []                       # (other-coord, side index or None, is_vertex)
        for r in range(n):
            a, b = P[r], P[(r + 1) % n]
            da = a[axis] - c
            db = b[axis] - c
            if abs(da) <= tol and abs(db) <= tol:
                # side parallel to the flow and ON the line: both endpoints are vertices
                hits.append((a[o], None, True))
                hits.append((b[o], None, True))
                continue
            if abs(da) <= tol:
                hits.append((a[o], None, True))
                continue
            if abs(db) <= tol:
                # `b` is ON the line: it is side `r+1`'s own `a` and is recorded there, as a
                # VERTEX.  Falling through would let `da*db < 0` fire on a `1e-16` residue
                # and append an interpolated hit at `t ~ 1` flagged `is_vertex = False`,
                # level with the genuine vertex hit -- so a saddle connection arriving at a
                # cone point reads as an ordinary side crossing and the walk runs past it.
                # Measured s528 at `Q = 7`, `axis = 0`, `double_Qgon`.  `marked_levels` and
                # `cylinders` are unaffected (checked, both builders, `Q = 5..31`).
                continue
            if da * db < 0:
                t = da / (da - db)
                hits.append((a[o] + t * (b[o] - a[o]), r, False))
        if not hits:
            return None
        hits.sort(key=lambda h: h[0])
        lo, hi = hits[0], hits[-1]
        if hi[0] - lo[0] <= tol:
            return None                 # the line only grazes a vertex
        return (lo[0], hi[0], lo[1], hi[1], lo[2], hi[2])

    # -- the separatrix-cut levels --------------------------------------
    def marked_levels(self, axis: int, tol: float = TOL,
                      cap: int = 200000) -> List[List[float]]:
        """Per polygon, the sorted CLOSURE of the vertex coordinates under the flow.

        Seeds every vertex coordinate of every polygon, then pushes each level up and down
        through the gluing until it dies on a cone point.  Raises if the closure exceeds
        `cap` (the direction is not periodic, or the model is wrong)."""
        npoly = len(self.polys)
        seen: List[Dict[int, float]] = [dict() for _ in range(npoly)]
        work: List[Tuple[int, float]] = []

        def key(c: float) -> int:
            return int(round(c / tol))

        def add(p: int, c: float) -> None:
            k = key(c)
            for kk in (k - 1, k, k + 1):
                if kk in seen[p]:
                    return
            seen[p][k] = c
            work.append((p, c))

        for p, P in enumerate(self.polys):
            for v in P:
                add(p, v[axis])

        while work:
            if sum(len(s) for s in seen) > cap:
                raise RuntimeError("marked-level closure exceeded cap -- direction not periodic?")
            p, c = work.pop()
            ch = self.chord(p, axis, c, tol)
            if ch is None:
                continue
            lo, hi, side_lo, side_hi, _vl, _vh = ch
            for (side, end) in ((side_hi, hi), (side_lo, lo)):
                if side is None:
                    continue            # the separatrix dies on a cone point
                q, _s = self.glue[(p, side)]
                sh = self.shift(p, side)
                pt = [0.0, 0.0]
                pt[axis] = c
                pt[1 - axis] = end
                add(q, pt[axis] + sh[axis])

        out = []
        for p in range(npoly):
            vals = sorted(seen[p].values())
            merged = [vals[0]]
            for v in vals[1:]:
                if v - merged[-1] > tol:
                    merged.append(v)
            out.append(merged)
        return out

    # -- the decomposition ----------------------------------------------
    def cylinders(self, axis: int, tol: float = TOL) -> List[Dict]:
        """Cylinder decomposition of the `axis`-direction, by separatrix cut.

        Each entry: `members` (list of `(poly, lo, hi)` sub-strips), `height` (their common
        width), `circ` (circumference), `legs` (per sub-strip `(poly, c, lo, hi)` -- the
        chord at the sub-strip's midpoint), `nstrip`.  Sorted by increasing height."""
        levels = self.marked_levels(axis, tol)
        strips: List[Tuple[int, float, float]] = []
        for p, lv in enumerate(levels):
            for i in range(len(lv) - 1):
                strips.append((p, lv[i], lv[i + 1]))
        index: Dict[int, List[Tuple[int, float, float]]] = {}
        for i, (p, a, b) in enumerate(strips):
            index.setdefault(p, []).append((i, a, b))

        def locate(p: int, c: float) -> int:
            for (i, a, b) in index.get(p, ()):
                if a - tol <= c <= b + tol:
                    return i
            raise RuntimeError("landed outside every sub-strip")

        # first-return permutation on sub-strips (flow in +direction)
        nxt: List[int] = [-1] * len(strips)
        leg: List[Tuple[int, float, float, float]] = [(0, 0.0, 0.0, 0.0)] * len(strips)
        for i, (p, a, b) in enumerate(strips):
            c = 0.5 * (a + b)
            ch = self.chord(p, axis, c, tol)
            assert ch is not None, "empty sub-strip"
            lo, hi, _sl, side_hi, _vl, _vh = ch
            leg[i] = (p, c, lo, hi)
            assert side_hi is not None, "sub-strip midpoint exits at a vertex"
            q, _s = self.glue[(p, side_hi)]
            sh = self.shift(p, side_hi)
            nxt[i] = locate(q, c + sh[axis])

        seen = [False] * len(strips)
        out: List[Dict] = []
        for i0 in range(len(strips)):
            if seen[i0]:
                continue
            cyc, i = [], i0
            while not seen[i]:
                seen[i] = True
                cyc.append(i)
                i = nxt[i]
            assert i == i0, "first-return map is not a permutation"
            ws = [strips[j][2] - strips[j][1] for j in cyc]
            if max(ws) - min(ws) > 1e-7 * max(1.0, max(ws)):
                raise RuntimeError(f"cycle mixes widths {sorted(set(round(w, 9) for w in ws))}")
            circ = sum(leg[j][3] - leg[j][2] for j in cyc)
            out.append({"members": [strips[j] for j in cyc],
                        "idx": cyc,
                        "height": sum(ws) / len(ws),
                        "circ": circ,
                        "legs": [leg[j] for j in cyc],
                        "nstrip": len(cyc)})
        out.sort(key=lambda c: c["height"])
        return out


# --------------------------------------------------------------- crossing matrices

def _clip(poly: List[Pt], axis: int, lo: float, hi: float) -> List[Pt]:
    """Sutherland-Hodgman clip of a convex polygon to the slab `lo <= coord[axis] <= hi`."""
    for sgn, bound in ((1, lo), (-1, hi)):
        out: List[Pt] = []
        n = len(poly)
        if n == 0:
            return out
        for i in range(n):
            a, b = poly[i], poly[(i + 1) % n]
            fa = sgn * (a[axis] - bound)
            fb = sgn * (b[axis] - bound)
            if fa >= 0:
                out.append(a)
            if (fa > 0 > fb) or (fa < 0 < fb):
                t = fa / (fa - fb)
                out.append((a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1])))
        poly = out
    return poly


def _area(poly: Sequence[Pt]) -> float:
    s = 0.0
    for i in range(len(poly)):
        a, b = poly[i], poly[(i + 1) % len(poly)]
        s += a[0] * b[1] - b[0] * a[1]
    return abs(s) / 2


def crossing_boxes(S: Surface, H: List[Dict], V: List[Dict]
                   ) -> Tuple[List[List[float]], List[List[float]]]:
    """`C[k][t] = area(H_k n V_t)/(h_k v_t)` by clipping each sub-strip BOX against its
    polygon.  Also returns the multiset of per-box area FRACTIONS, which is what a
    box-lemma-style proof has to pin down."""
    C = [[0.0] * len(V) for _ in H]
    fracs = [[[] for _ in V] for _ in H]
    for k, hc in enumerate(H):
        for (p, a, b) in hc["members"]:
            for t, vc in enumerate(V):
                for (p2, a2, b2) in vc["members"]:
                    if p2 != p:
                        continue
                    piece = _clip(_clip(S.polys[p], 0, a, b), 1, a2, b2)
                    if len(piece) < 3:
                        continue
                    ar = _area(piece)
                    box = (b - a) * (b2 - a2)
                    if ar > 1e-13 * max(1.0, box):
                        C[k][t] += ar / box
                        fracs[k][t].append(ar / box)
    return C, fracs


def crossing_walk(S: Surface, H: List[Dict], V: List[Dict]) -> List[List[float]]:
    """`C[k][t]` as an ARCLENGTH ratio along a closed leaf of `H_k`: independent of
    `crossing_boxes` (it uses the leg chords and the transverse sub-strip bands, never a
    polygon clip)."""
    band: Dict[Tuple[int, float, float], int] = {}
    for t, vc in enumerate(V):
        for m in vc["members"]:
            band[m] = t
    C = [[0.0] * len(V) for _ in H]
    for k, hc in enumerate(H):
        for (p, _c, lo, hi) in hc["legs"]:
            for (p2, a2, b2), t in band.items():
                if p2 != p:
                    continue
                ov = min(hi, b2) - max(lo, a2)
                if ov > 1e-13:
                    C[k][t] += ov
    for k in range(len(H)):
        for t in range(len(V)):
            C[k][t] /= V[t]["height"]
    return C


# ------------------------------------------------------------------ model builders

def regular_2Qgon(Q: int, R: float | None = None) -> Surface:
    """The regular `2Q`-gon with opposite sides identified -- the `P in {1, Q-1}` arm
    ([NCYL-325] (1)).  Default `R = 1/sin(pi/2Q)`, the repo's `L1 = 1` normalisation."""
    R = (1.0 / math.sin(math.pi / (2 * Q))) if R is None else R
    M = 2 * Q
    V = [(R * math.cos(r * math.pi / Q), R * math.sin(r * math.pi / Q)) for r in range(M)]
    glue = {(0, r): (0, (r + Q) % M) for r in range(M)}
    return Surface([V], glue)


def double_Qgon(Q: int, R: float | None = None) -> Surface:
    """Two regular `Q`-gons related by `z -> -z`, side `r` glued to side `r` -- the
    `P in {2, Q-2}` arm ([NCYL-325] (1), [NCYL-425] (9)).  Oriented so that side `0` is
    HORIZONTAL, hence a vertical symmetry axis through the opposite vertex."""
    R = (1.0 / math.sin(math.pi / Q)) if R is None else R
    th0 = -math.pi / 2 - math.pi / Q
    A = [(R * math.cos(th0 + 2 * math.pi * r / Q),
          R * math.sin(th0 + 2 * math.pi * r / Q)) for r in range(Q)]
    B = [(-x, -y) for (x, y) in A]
    glue: Dict[Tuple[int, int], Tuple[int, int]] = {}
    for r in range(Q):
        glue[(0, r)] = (1, r)
        glue[(1, r)] = (0, r)
    return Surface([A, B], glue)
