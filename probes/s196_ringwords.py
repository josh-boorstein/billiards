"""Exact ring-backed enumerator of ALL perpendicular cylinders WITH their frozen
words — the reliable-at-depth replacement for exp_words.build/veech_cylinders.partition,
which are float-depth-limited and UNDERCOUNT at large Q (s196: 10/33 true n=32, exp_words
gives 13). Reuses orphan_region's EXACT ring tree loop (RingTriangleState, O(D^2)/reflect,
depth-independent) and just records the reflection word at each leaf.

Word alphabet {H,1,2} = universal across all P/Q (H=head side, 1=L1, 2=L2).

Run: PYTHONPATH=.:archive/scripts_2026-07:probes .venv/bin/python3.13 probes/s196_ringwords.py
"""
import math, json
from math import gcd
from orphan_region import (SIDES_AFTER, SPLIT_VERTEX, _hit_side_ring)
from exact_state_ring import RingTriangleState
from ring_cyc import RingContext

WMAP = {'H': 'H', 'L1': '1', 'L2': '2'}
def wstr(word): return ''.join(WMAP[s] for s in word)


import mpmath as mp


def _closed(state, closure_tol, exact):
    """Branch closes iff the developed frame returns to identity (dev O->R = +-(1,0),
    O->A = +-(1,1)). Float pre-filter (fast) then HIGH-PRECISION confirm (tol 1e-25) so
    near-retraces do NOT close early (the s196 10/33 premature-closure trap)."""
    dx = (state.Rx_c - state.Ox_c).to_float(); dy = (state.Ry - state.Oy).to_float()
    dax = (state.Ax_c - state.Ox_c).to_float(); day = (state.Ay - state.Oy).to_float()
    t = closure_tol
    near = ((abs(dx-1) < t and abs(dy) < t and abs(dax-1) < t and abs(day-1) < t) or
            (abs(dx+1) < t and abs(dy) < t and abs(dax+1) < t and abs(day-1) < t))
    if not near or not exact:
        return near
    dxe = (state.Rx_c - state.Ox_c).to_mpf(dps=40); dye = (state.Ry - state.Oy).to_mpf(dps=40)
    daxe = (state.Ax_c - state.Ox_c).to_mpf(dps=40); daye = (state.Ay - state.Oy).to_mpf(dps=40)
    T = mp.mpf(10) ** -25
    return ((abs(dxe-1) < T and abs(dye) < T and abs(daxe-1) < T and abs(daye-1) < T) or
            (abs(dxe+1) < T and abs(dye) < T and abs(daxe+1) < T and abs(daye-1) < T))


def enumerate_branches(P, Q, max_depth=400000, min_width=1e-15, closure_tol=1e-6, exact=True,
                       keep_words=True):
    """All cylinders of S_{P/Q}: sorted list of (lo, hi, word_str). Exact ring state,
    high-precision closure (exact=True). Returns (leaves, n_unclosed) if diagnostics wanted.

    keep_words=False (s273) returns the word LENGTH (int) in slot 2 instead of the word
    string. Materializing the word costs `word + (side,)` per reflection, i.e. O(depth^2)
    copying per branch, which DOMINATES beyond depth ~1e5: measured on 8/19,
    t ~ 2.34e-4*D + 9.45e-9*D^2, so D=2e6 is 10.6 h with words and ~8 min without.
    The length is free and exact -- len(word) == depth identically, since the root is
    (depth 1, word ('H',)) and every push increments both.

    ⇒⇒ `enumerate_branches.last_grazes` (s373) -- THE VERTEX GRAZE THAT DEFINES EACH
    INTERIOR BOUNDARY, which this traversal has always computed and always discarded.
    A dict `{round(y, 10): (vertex, k)}` over the interior cell boundaries, where
    `vertex` is 'O' / 'R' / 'A' and `k` is the number of REFLECTIONS ALREADY MADE when
    the split happened -- so the ray launched at height `y` grazes `vertex` on the step
    AFTER its `k`-th bounce, i.e. `k` = bounces before the graze.  Free: the split point
    IS `state.vertex_y(vtx)` at a known depth (`foundations.md` §1 -- a cell boundary is
    a zero of a vertex-graze function), so this records what the split already knew.
    Set on every call; no kwarg, no signature change, same idiom as `last_unclosed`.
    ⚠ The two TRANSVERSAL ENDS (`s = 0`, `s = 1`) are NOT in the dict: they are the `L1`
    leg's own endpoints `R` and `A`, not grazes.  So `len(last_grazes) == n_cells - 1`.
    ⚠⚠ `k` is the FIRST/DEFINING graze, and that is a real choice: at a boundary the
    orbit typically runs along the corner for MANY consecutive near-vertex bounces
    (median 8, max 21 measured over 60 orphan endpoints), so "the closest approach" is an
    argmin over a cluster of near-ties and is NOT well defined -- that is exactly what
    made `probes/s332_vertexfree.py`'s `D_*_hit`/`D_*_vertex` noise-selected ([OPS-096]).
    Validated against the independent float tracer, 60/60 on both vertex and index
    (`k + 1` == the first near-vertex hit index)."""
    ctx = RingContext(P, Q)
    P, Q = ctx.P, ctx.Q
    grazes = {}
    degenerate = 0          # s273: leaves terminated by _hit_side_ring None (NOT closure)
    cot = 1.0 / math.tan(ctx.alpha)
    init = RingTriangleState(ctx); init.reflect('H')
    ext = (lambda w, s: w + (s,)) if keep_words else (lambda w, s: w + 1)
    stack = [(0.0, 1.0, init, 'H', 1, ('H',) if keep_words else 1)]
    leaves = []; unclosed = 0
    while stack:
        lo, hi, state, prev, depth, word = stack.pop()
        if hi - lo < min_width or depth > max_depth:
            leaves.append((lo, hi, word)); unclosed += (depth > max_depth); continue
        if depth > 3 and _closed(state, closure_tol, exact):
            leaves.append((lo, hi, word)); continue
        cand = SIDES_AFTER[prev]
        if len(cand) == 1:
            side = cand[0]; ns = state.copy(); ns.reflect(side)
            stack.append((lo, hi, ns, side, depth + 1, ext(word, side))); continue
        vtx = SPLIT_VERTEX[frozenset(cand)]
        sy = state.vertex_y(vtx).to_float()
        if sy <= lo + min_width or sy >= hi - min_width:
            side = _hit_side_ring(state, cand, 0.5 * (lo + hi), cot)
            if side is None: leaves.append((lo, hi, word)); degenerate += 1; continue
            ns = state.copy(); ns.reflect(side)
            stack.append((lo, hi, ns, side, depth + 1, ext(word, side))); continue
        shi = _hit_side_ring(state, cand, 0.5 * (sy + hi), cot)
        slo = _hit_side_ring(state, cand, 0.5 * (lo + sy), cot)
        if shi is None or slo is None:
            leaves.append((lo, hi, word)); degenerate += 1; continue
        grazes[round(sy, 10)] = (vtx, depth)      # s373: the boundary's DEFINING graze
        s2 = state.copy(); s2.reflect(shi); stack.append((sy, hi, s2, shi, depth + 1, ext(word, shi)))
        s2 = state.copy(); s2.reflect(slo); stack.append((lo, sy, s2, slo, depth + 1, ext(word, slo)))
    leaves = [lf for lf in leaves if lf[1] - lf[0] > 1e-13]
    leaves.sort()
    out = [(lo, hi, wstr(w) if keep_words else w) for (lo, hi, w) in leaves]
    enumerate_branches.last_unclosed = unclosed
    enumerate_branches.last_degenerate = degenerate
    enumerate_branches.last_grazes = grazes
    return out


def ab(word):
    return (word.count('H'), word.count('1'), word.count('2'))


def from_cf(digits):
    p, q = 0, 1
    for a in reversed(digits):
        p, q = q, a * q + p
    g = gcd(p, q)
    return p // g, q // g


def can_factor(word, vocab):
    """True if `word` segments into a concatenation of elements of `vocab` (DP)."""
    lens = sorted({len(v) for v in vocab})
    vset = set(vocab)
    n = len(word)
    ok = [False] * (n + 1); ok[0] = True
    for i in range(n):
        if not ok[i]:
            continue
        for L in lens:
            if i + L <= n and word[i:i + L] in vset:
                ok[i + L] = True
    return ok[n]


def family_view(name, convergents, save=None):
    """Enumerate each convergent, show compact per-cylinder fingerprints + verbatim
    persistence vs previous convergent + factorization over previous vocabulary."""
    import time
    print(f"\n############  {name}  ############")
    prev_words = None
    rec = []
    for (P, Q) in convergents:
        t = time.time()
        br = enumerate_branches(P, Q)
        dt = time.time() - t
        words = [w for _, _, w in br]
        widths = [hi - lo for lo, hi, _ in br]
        unc = enumerate_branches.last_unclosed
        prevset = set(prev_words) if prev_words else set()
        print(f"\n  {P}/{Q}:  n={len(br)}  unclosed={unc}  ({dt:.1f}s)")
        persist = 0
        for (lo, hi, w), wd in sorted(zip(br, widths), key=lambda x: -x[1]):
            a = ab(w)
            tag = "PERSIST" if w in prevset else ("factors" if prev_words and can_factor(w, prev_words) else "new")
            persist += (w in prevset)
            print(f"     ab={str(a):14s} wlen={len(w):6d} width={wd:.6e}  {tag}")
        if prev_words:
            nfac = sum(1 for w in words if can_factor(w, prev_words))
            print(f"     -> {persist}/{len(words)} persist verbatim; {nfac}/{len(words)} factor over previous vocabulary")
        rec.append(dict(P=P, Q=Q, n=len(br), words=words, widths=widths, unclosed=unc))
        prev_words = words
    if save:
        with open(save, "w") as f:
            json.dump(rec, f, default=str)
        print(f"\n  wrote {save}")
    return rec


if __name__ == '__main__':
    import time
    # VALIDATION: n vs known small values.
    # WARNING (s364): this line used to read "n must match the even-P law (n=Q-1)".
    # That law is Thm B', REFUTED s272 -- n(8/15) = 6, not 14 (claims.md [T-N2]).
    # Using it as a gate REJECTS A GENUINE EVEN-P DEFICIT ROW AS A BUG, which is the
    # practice rulings.md [OPS-001] bars. The rows below happen to be Q-1-conforming,
    # so the gate was never triggered -- do NOT extend this list on the strength of
    # the law. Proved even-P content: n is EVEN. n(2/Q) = Q-1 (Thm C) still holds.
    checks = [(2, 3, 2), (3, 5, 3), (5, 8, 3), (8, 13, 12), (2, 5, 4),
              (5, 12, 5), (12, 29, 28), (10, 33, 32), (3, 10, 5)]
    print("VALIDATION (n vs expected):")
    for (P, Q, exp) in checks:
        t = time.time()
        br = enumerate_branches(P, Q)
        print(f"  {P}/{Q}: n={len(br)}  expected={exp}  {'OK' if len(br)==exp else 'MISMATCH'}  {time.time()-t:.1f}s")
