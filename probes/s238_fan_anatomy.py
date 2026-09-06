"""s238_fan_anatomy.py -- WHERE the landmark depths (2m-1, 4m-1, 6m-5) come from.

s237 (4) found, on every passing row of both s236/s237 sweeps (all r2, both parities of P,
every m), that the three top landmark boundaries are created at ring depths exactly

    dep(H, S, B1) = (2m-1, 4m-1, 6m-5)     independent of (P,r).

This probe exposes the mechanism: it traces the single ray at a height inside the beat cell
(B1, S) and prints, at every reflection, the developed geometry -- the O-image height, the
O->R tilt in units of pi/Q, the split vertex and its height.  The expected structure is a
decomposition into FANS: long runs of alternating H/L2 reflections that rotate the frame by
2*alpha per step around a fixed O-image (Oy constant, tilt stepping by +/-P), separated by
short transitions.  The landmark depths should then be readable as exact fan lengths.

Run: PYTHONPATH=.:engine:archive/scripts_2026-07:archive/scripts_2026-08:probes \
     .venv/bin/python3.13 -u probes/s238_fan_anatomy.py > logs/s238_fan_anatomy.log 2>&1
"""
import math
import sys

from orphan_region import SIDES_AFTER, SPLIT_VERTEX, _hit_side_ring
from exact_state_ring import RingTriangleState
from ring_cyc import RingContext


def trace(P, r, m, s, maxdepth=200):
    """Follow the ray at height s; return the per-depth record."""
    Q = P * m + r
    ctx = RingContext(P, Q)
    cot = 1.0 / math.tan(ctx.alpha)
    st = RingTriangleState(ctx)
    st.reflect('H')
    prev, depth = 'H', 1
    rec = []
    while depth <= maxdepth:
        cand = SIDES_AFTER[prev]
        vtx = SPLIT_VERTEX[frozenset(cand)]
        sy = st.vertex_y(vtx).to_float()
        Oy = st.Oy.to_float()
        dx_c = (st.Rx_c - st.Ox_c).to_float()
        dy = (st.Ry - st.Oy).to_float()
        tilt = math.atan2(dy, dx_c * cot) * Q / math.pi
        side = _hit_side_ring(st, cand, s, cot)
        rec.append(dict(depth=depth, prev=prev, vtx=vtx, sy=sy, Oy=Oy,
                        Ox_c=st.Ox_c.to_float(), tilt=tilt, side=side))
        if side is None:
            break
        st.reflect(side)
        prev = side
        depth += 1
    return ctx, rec


def landmarks(P, r, m):
    """Closed-form H, S, B1 (s236 (2))."""
    Q = P * m + r
    th = math.pi / (2.0 * Q)
    cot = 1.0 / math.tan(P * th)
    s_ = lambda k: cot * math.sin(math.pi * k / Q)
    H = s_(r)
    S = 2 * H - 1
    B1 = 2 * H - 2 * s_(P) + s_(P + r)
    return H, S, B1


def report(P, r, m, maxdepth=None):
    Q = P * m + r
    H, S, B1 = landmarks(P, r, m)
    s = 0.5 * (S + B1)
    if maxdepth is None:
        maxdepth = 8 * m + 20
    ctx, rec = trace(P, r, m, s, maxdepth)
    print(f"\n=== (P,r)=({P},{r}) r2={P-2*r} m={m} Q={Q}  "
          f"H={H:.9f} S={S:.9f} B1={B1:.9f}  ray s={s:.9f}")
    print(f"    predicted depths: H {2*m-1}  S {4*m-1}  B1 {6*m-5}")
    print(f"  {'d':>4} {'prev':>4} {'vtx':>3} {'split y':>13} {'Oy':>13} "
          f"{'Ox_c':>10} {'tiltQ/pi':>10} {'->':>4}  mark")
    for e in rec:
        mark = ''
        for nm, val in (('H', H), ('S', S), ('B1', B1)):
            if abs(e['sy'] - val) < 1e-9:
                mark = f'<== {nm}'
        onseg = 0.0 < e['sy'] < 1.0
        print(f"  {e['depth']:4d} {e['prev']:>4} {e['vtx']:>3} {e['sy']:13.8f} "
              f"{e['Oy']:13.8f} {e['Ox_c']:10.4f} {e['tilt']:10.4f} {str(e['side']):>4}"
              f"  {'*' if onseg else ' '} {mark}")


if __name__ == '__main__':
    if len(sys.argv) > 1:
        P, r, m = (int(x) for x in sys.argv[1:4])
        report(P, r, m)
    else:
        for (P, r, m) in [(5, 2, 5), (5, 2, 6), (7, 3, 5), (19, 8, 4)]:
            report(P, r, m)
