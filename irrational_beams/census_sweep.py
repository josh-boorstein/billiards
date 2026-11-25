"""[Deposit] Section 5, the census: the grid sweep and the exact stage (exact_row). Paper: "Finite perpendicular-beam partitions at irrational angles". Research name: probes/s977_connection_census.py, copied verbatim except imports (deposit names), paths (data/ and figs/ under irrational_beams/) and this paragraph. Run through run_all.py.

s977 — CONNECTION CENSUS: finite L1 partitions found from their A–O saddle connection (user s977: *we can build the census using
connections*; queue s976-c, re-scoped).

WHY. §s977 (1): every finite point carries ONE equation, a horizontal A–O saddle connection E = y(gA) − y(g′O) (E(0) = ±1, so every root is
a unit), realised at exactly one root. So the connections are a SUPERSET of the finite points that reaches any degree — unlike the unit
census (complete in a degree × height box, §s976 (5)) and the coincidence censuses (first-leg heights only).
METHOD. For α on a grid, a depth-limited beam splitter carries pieces (an s-interval of parallel rays, one word) WITH their words; a piece
entering from side L1 / L2 / H can split only at O / A / R. Each piece knows the events bounding it (an earlier split, or s = 0 — the leaf
past O — or s = 1 — the leaf under A). An O-split of a piece bounded by an A-event (or an A-split of one bounded by an O-event) at gap < GAP
is a CANDIDATE: its two fixed-word heights h_A(α), h_O(α) are analytic, and their difference is solved for a root near the grid point
(secant). A root is kept iff (i) both words are REALISED there (every earlier side parameter in [0, 1], the leaf height in [0, 1]); then
(ii) a quick float finiteness screen (trace_nb, ε = 1e-6, small budget); then (iii) the exact pipeline (E by exact interpolation in σ,
minimal polynomial, unit / irrational test, exact_partition certificate).
SCOPE (what "complete" means here): every finite point in the x-range whose connection pair is a piece-end/split pair of the depth-≤ N,
width-≥ EPS splitting at some grid point with gap < GAP. Thin one-sided windows (§s976 (4): wedge tips ~1e-5) can fall between grid points.
PRE-REGISTERED.
  (C1) CONTROL: every known finite point in range whose connection steps are ≤ N is recovered — the 12 unit-census points, the [WIND-FAM]
       points, the F(m, J) points of the range. A miss is an instrument failure (grid, depth or GAP), reported with which.
  (C2) every recovered finite point carries a unit minimal polynomial (Theorem I′) — cannot fail if the certificate is right.
  (C3) THE QUESTION: how many finite points does the census find that NO previous census found, and how many realised A–O connections
       are NOT finite (the realisability-vs-finiteness gap)? Could come out otherwise: no new points (the known families exhaust the
       range at this depth), or most realised connections finite.
RESULT (s977, window x ∈ [0.28, 0.32], DX = 1e-4, N = 5000 → data/s977_connection_census_0.28_0.32.json; scaling runs on
[0.300, 0.304] → data/scratch/s977/scal_*.json):
  (C1) HOLDS on the window: 0.2923 (WF(3,2)), 0.3022 (WF(3,4)), 0.3054 (α₃) recovered (4, 4, 2 branches).
  (C2) 8/8 sampled CERTIFIED exactly (units, irrational; degrees 6–21, e.g. σ⁶ − 8σ⁵ + 15σ⁴ + 17σ³ − 56σ² + 18σ − 1 at x = 0.2805982960).
  (C3) 506 realised A–O connections, 227 finite by C0.test (ε = 1e-6 ⊕ 1e-9, equal counts) — 279 realised connections NOT finite; every
       finite branch count EVEN. NO SATURATION: on [0.300, 0.304] the finite count is 19 / 107 / 31 / 142 at (DX, N) = (1e-4, 5000) /
       (2.5e-5, 5000) / (1e-4, 2e4) / (2.5e-5, 2e4); by branch count n = 4 nearly saturates (2 → 3), n = 6 grows 2 → 21 and those
       cluster on both sides of WF(3,4) = 0.302222 (finite points accumulate at finite points). A total count is therefore scope-relative;
       a fixed-n count is the census quantity. Exact stage ~1 min / point (2 cores).
PRIOR ART: §s977 (1) (one connection per point; realisability), §s976 (5) (unit census), §s973 (2)/§s976 (4) (coincidence censuses);
rulings --grep 'connection census', 'separatrix census' → none; data_index --grep connection → none.
DATA: data/s976_unit_census_D3_H6_certified.json, data/s977_winding_family.json, data/s973_coincidence_census_K20_certified_hardened.json,
data/s976_census_ext_K20.json (the control lists).
Usage: python probes/census_sweep.py sweep X0 X1 DX [N] [NPROC]  → data/s977_connection_census_<X0>_<X1>.json
       python probes/census_sweep.py exact <json>                → certifies the screened survivors in place
"""
import _setup  # noqa: F401  (deposit: this directory on the path, cwd = irrational_beams/)
import json
import math
import os
import sys
import time
from multiprocessing import Pool

sys.path[:0] = [os.path.dirname(os.path.abspath(__file__))]

EPS = 1e-7          # pieces narrower than this are dropped
GAP = 0.03          # candidate pair gap (in s)
SIDES_VX = dict(L1=('R', 'A'), L2=('R', 'O'), H=('O', 'A'))
OPP = dict(L1='O', L2='A', H='R')     # the vertex a piece entering from this side can split at


def geom(a):
    cot = math.cos(a) / math.sin(a)
    R, A, O = (0.0, 0.0), (0.0, 1.0), (-cot, 0.0)
    return dict(L1=(R, A), L2=(R, O), H=(O, A))


def split_sweep(a, N, eps=EPS):
    """Depth-≤ N beam splitting at angle a. Returns candidate pairs: (gap, (wordA, kA), (wordO, kO), hA, hO)."""
    S = geom(a)
    # piece: lo, hi, px0, px1, py0, py1, dx, dy, last, k, word(list), ev_lo, ev_hi ; event = (type, word, k, h)
    ev0 = ('O', ('H', None), 1, 0.0)        # s = 0: the step-1 exit side H passes O (the leaf along L2 reaches O)
    ev1 = ('A', None, 0, 1.0)       # s = 1: the launch corner A (height identically 1)
    REFL = dict(L1=lambda A, B: (1 - A, -B), L2=lambda A, B: (-A, -B), H=lambda A, B: (-A, 2 - B))   # φ = Aπ + Bα ↦ 2ψ − φ
    stack = [(0.0, 1.0, 0.0, 0.0, 0.0, 1.0, (1, 0), 'L1', 0, None, ev0, ev1)]      # words are cons cells (side, prev)
    cands = []
    split_sweep.steps = split_sweep.splits = split_sweep.capped = 0
    while stack:
        lo, hi, px0, px1, py0, py1, ph, last, k, word, elo, ehi = stack.pop()
        split_sweep.steps += 1
        phi = ph[0] * math.pi + ph[1] * a
        dx, dy = math.cos(phi), math.sin(phi)
        if hi - lo < eps:
            continue
        if k >= N:
            split_sweep.capped += 1
            continue
        # exit side at the two ends: try each side ≠ last
        best = {}
        for name, (P, Q) in S.items():
            if name == last:
                continue
            ex, ey = Q[0] - P[0], Q[1] - P[1]
            den = dx * ey - dy * ex
            if abs(den) < 1e-300:
                continue
            wx0, wx1, wy0, wy1 = P[0] - px0, -px1, P[1] - py0, -py1
            u0 = (wx0 * dy - wy0 * dx) / den; u1 = (wx1 * dy - wy1 * dx) / den
            t0 = (wx0 * ey - wy0 * ex) / den; t1 = (wx1 * ey - wy1 * ex) / den
            best[name] = (u0, u1, t0, t1, ex, ey)
        names = list(best)
        # which side owns s: u ∈ [0, 1] and t > 0 ; at most one split (the vertex opposite `last`)
        V = OPP[last]
        nm1 = names[0]
        u0, u1 = best[nm1][0], best[nm1][1]
        # the split height: where side nm1's parameter hits its V-end
        vend = 0.0 if SIDES_VX[nm1][0] == V else 1.0
        sv = (vend - u0) / u1 if u1 != 0 else None
        pieces = []
        if sv is not None and lo < sv < hi:
            ev = (V, (nm1, word), k + 1, sv)
            split_sweep.splits += 1
            pieces = [(lo, sv, elo, ev), (sv, hi, ev, ehi)]
            # candidate: an A/O split adjacent to an O/A end
            if V in 'AO':
                for end in (elo, ehi):
                    if end[0] in 'AO' and end[0] != V and abs(end[3] - sv) < GAP:
                        cands.append((abs(end[3] - sv), end, ev))
        else:
            pieces = [(lo, hi, elo, ehi)]
        for (l, h, el, eh) in pieces:
            m = 0.5 * (l + h)
            name = None
            for nm, (u0, u1, t0, t1, ex, ey) in best.items():
                u = u0 + u1 * m; t = t0 + t1 * m
                if t > 0 and -1e-12 <= u <= 1 + 1e-12:
                    name = nm
                    break
            if name is None:
                continue
            u0, u1, t0, t1, ex, ey = best[name]
            nx0, nx1 = px0 + t0 * dx, px1 + t1 * dx
            ny0, ny1 = py0 + t0 * dy, py1 + t1 * dy
            nph = REFL[name](*ph)
            if name == 'L1' and nph[1] == 0 and nph[0] % 2 == 1:      # head-on return (symbolic, as trace_nb): done
                continue
            stack.append((l, h, nx0, nx1, ny0, ny1, nph, name, k + 1, (name, word), el, eh))
    return cands


def unc(w):
    out = []
    while w is not None:
        out.append(w[0]); w = w[1]
    return tuple(reversed(out))


def word_height(word, k, V, a):
    """fixed-word height (float) at which step k's side passes vertex V; k = 0 means the launch corner (height 1)."""
    if k == 0:
        return 1.0, 1.0
    S = geom(a)
    px0, px1, py0, py1, dx, dy = 0.0, 0.0, 0.0, 1.0, -1.0, 0.0
    worst = 1.0
    for j, name in enumerate(word[:k], start=1):
        P, Q = S[name]
        ex, ey = Q[0] - P[0], Q[1] - P[1]
        den = dx * ey - dy * ex
        wx0, wx1, wy0, wy1 = P[0] - px0, -px1, P[1] - py0, -py1
        u0 = (wx0 * dy - wy0 * dx) / den; u1 = (wx1 * dy - wy1 * dx) / den
        t0 = (wx0 * ey - wy0 * ex) / den; t1 = (wx1 * ey - wy1 * ex) / den
        if j == k:
            vend = 0.0 if SIDES_VX[name][0] == V else 1.0
            return (vend - u0) / u1, (u0, u1, t0, t1, px0, px1)
        px0, px1, py0, py1 = px0 + t0 * dx, px1 + t1 * dx, py0 + t0 * dy, py1 + t1 * dy
        ee = ex * ex + ey * ey; dot = dx * ex + dy * ey
        dx, dy = dot * ex * 2 / ee - dx, dot * ey * 2 / ee - dy
    raise ValueError


def realised(word, k, a, s):
    """min over steps j < k of the side parameter margins at height s (≥ 0 ⇔ the leaf follows the word up to step k)."""
    if k == 0:
        return 1.0
    S = geom(a)
    px, py, dx, dy = 0.0, s, -1.0, 0.0
    worst = 1.0
    for j, name in enumerate(word[:k - 1], start=1):
        P, Q = S[name]
        ex, ey = Q[0] - P[0], Q[1] - P[1]
        den = dx * ey - dy * ex
        wx, wy = P[0] - px, P[1] - py
        u = (wx * dy - wy * dx) / den; t = (wx * ey - wy * ex) / den
        worst = min(worst, u, 1 - u)
        if t <= 0:
            return -1.0
        px, py = px + t * dx, py + t * dy
        ee = ex * ex + ey * ey; dot = dx * ex + dy * ey
        dx, dy = dot * ex * 2 / ee - dx, dot * ey * 2 / ee - dy
    return worst


def solve_pair(eA, eO, a0, da):
    f = lambda a: word_height(eO[1], eO[2], 'O', a)[0] - word_height(eA[1], eA[2], 'A', a)[0]
    x0, x1 = a0, a0 + da * 0.5
    f0, f1 = f(x0), f(x1)
    for _ in range(60):
        if f1 == f0:
            break
        x2 = x1 - f1 * (x1 - x0) / (f1 - f0)
        if not (0.01 < x2 < 1.56):
            return None
        x0, f0, x1, f1 = x1, f1, x2, f(x2)
        if abs(f1) < 1e-14:
            break
    if abs(f1) > 1e-11 or abs(x1 - a0) > 3 * da:
        return None
    return x1


def grid_point(args):
    a, N, da = args
    out = []
    for gap, end, ev in split_sweep(a, N):
        end = (end[0], unc(end[1]), end[2], end[3]); ev = (ev[0], unc(ev[1]), ev[2], ev[3])
        eA, eO = (end, ev) if end[0] == 'A' else (ev, end)
        r = solve_pair(eA, eO, a, da)
        if r is None:
            continue
        s = word_height(eA[1], eA[2], 'A', r)[0]
        if not (-1e-9 <= s <= 1 + 1e-9):
            continue
        mA, mO = realised(eA[1], eA[2], r, s), realised(eO[1], eO[2], r, s)
        if min(mA, mO) < -1e-7:
            continue
        split_sweep(r, N)                                # pre-screen: at a finite point no piece wider than EPS reaches depth N
        out.append(dict(alpha=r, x=2 * r / math.pi, s=s, A=(' '.join(eA[1][:eA[2]]), eA[2]), O=(' '.join(eO[1][:eO[2]]), eO[2]),
                        margin=min(mA, mO), capped_at_root=split_sweep.capped))
    return out


def screen(row):
    import beam_families as B, launch_side as L
    r = B.trace_nb(L.table('right', row['alpha'], 0), eps_min=1e-6, budget=3 * 10 ** 6, cap=10 ** 5)
    row['screen_n'] = len(r['branches']); row['screen_capped'] = len(r['capped'])
    row['screen_finite'] = not r['capped']
    return row


def sweep(x0, x1, dx, N=100, nproc=2):
    t0 = time.time()
    da = dx * math.pi / 2
    off = (math.sqrt(5) - 1) / 2                     # irrational offset: no grid point at a small-denominator rational
    grid = [((x0 + (i + off) * dx) * math.pi / 2, N, da) for i in range(int(round((x1 - x0) / dx)))]
    found = {}
    with Pool(nproc) as P:
        for i, rows in enumerate(P.imap(grid_point, grid, chunksize=8)):
            for r in rows:
                key = round(r['alpha'], 10)
                if key not in found:
                    found[key] = r
            if i % 500 == 0:
                print(f'  {i}/{len(grid)} {time.time() - t0:.0f}s  connections {len(found)}', flush=True)
    rows = sorted(found.values(), key=lambda r: r['alpha'])
    pre = [r for r in rows if r['capped_at_root'] == 0]
    print(f'{len(rows)} realised A–O connections, {len(pre)} with no capped piece at the root ({time.time() - t0:.0f}s); screening', flush=True)
    with Pool(nproc) as P:
        scr = {round(r['alpha'], 10): r for r in P.imap(screen, pre, chunksize=4)}
    for r in rows:
        r.update(scr.get(round(r['alpha'], 10), dict(screen_finite=False)))
    fin = [r for r in rows if r['screen_finite']]
    print(f'{len(fin)} pass the float screen ({time.time() - t0:.0f}s)', flush=True)
    for r in rows:                                   # keep words only where they will be used (full-range runs: ~10⁴–10⁵ rows)
        if not r['screen_finite']:
            r['A'] = ('', r['A'][1]); r['O'] = ('', r['O'][1])
    out = 'data/s977_connection_census_%s_%s.json' % (x0, x1)
    json.dump(dict(x0=x0, x1=x1, dx=dx, N=N, EPS=EPS, GAP=GAP, rows=rows), open(out, 'w'), default=str)
    for r in fin:
        print(f"   x={r['x']:.10f} s={r['s']:.6f} n={r['screen_n']} A@{r['A'][1]} O@{r['O'][1]}")
    return out


if __name__ == '__main__' and sys.argv[1] == 'sweep':
    x0, x1, dx = float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4])
    N = int(sys.argv[5]) if len(sys.argv) > 5 else 100
    npr = int(sys.argv[6]) if len(sys.argv) > 6 else 2
    sweep(x0, x1, dx, N, npr)


# ---------- exact stage: E by exact interpolation (adaptive degree), minimal polynomial, certificate ----------
def exact_height_poly(word, k, V, dmax=400):
    """y of the vertex image as an exact polynomial in σ (rational part; the √D part must vanish) — degree found adaptively."""
    import sympy as sp
    from fractions import Fraction as Fr
    import fan_identity as FI
    S = sp.Symbol('S')
    if k == 0:
        return sp.Poly(1, S)
    w = list(word)
    d = 4
    while d <= dmax:
        nodes = [Fr(1, 7) + Fr(j, 13 * (d + 6)) * Fr(7, 2) for j in range(3 * d + 20)]
        nodes = [z for z in nodes if 0 < z < 4 and not FI._is_sq(z * (4 - z) / 4)][:d + 4]
        vals = [FI.word_graze(w, k, V, z) for z in nodes]
        if any(v.b != 0 for v in vals):
            raise ValueError('√D part ≠ 0')
        pts = [(sp.Rational(z.numerator, z.denominator), sp.Rational(v.a.numerator, v.a.denominator)) for z, v in zip(nodes, vals)]
        P = sp.Poly(sp.interpolate(pts[:d + 1], S), S)
        if all(P.eval(x) == y for x, y in pts[d + 1:]):
            return P
        d *= 2
    raise ValueError('degree > %d' % dmax)


def exact_row(row):
    import sympy as sp
    import exact_partition as X
    S = sp.Symbol('S')
    try:
        wa, ka = row['A'][0].split(), row['A'][1]
        wo, ko = row['O'][0].split(), row['O'][1]
        E = exact_height_poly(wa, ka, 'A') - exact_height_poly(wo, ko, 'O')
        sig = 4 * math.cos(row['alpha']) ** 2
        best = None
        for fac, _ in sp.factor_list(E.as_expr(), S)[1]:
            fp = sp.Poly(fac, S)
            try:
                if any(r.is_real and abs(float(r) - sig) < 1e-8 for r in fp.nroots(n=30, maxsteps=500)):
                    best = fp
            except Exception:                                # s979: nroots non-convergence (s978's x = 0.1211654606) →
                lo, hi = sp.Rational(sig - 1e-8), sp.Rational(sig + 1e-8)   # exact sign change on the bracket instead
                if fp.eval(lo) * fp.eval(hi) <= 0:
                    best = fp
        if best is None:
            return dict(exact='no factor at σ', E_deg=E.degree())
        co = [int(c) for c in best.all_coeffs()]
        g = math.gcd(*co); co = [c // g for c in co]
        if co[0] < 0:
            co = [-c for c in co]
        out = dict(minpoly=co, E_deg=E.degree(), E_const=str(E.all_coeffs()[-1]), **X.arith(co, sig))
        c = X.run_poly(co, sig - 1e-9, sig + 1e-9)
        out.update(CERTIFIED=bool(c['CERTIFIED']), n_exact=c.get('n'), sign_margin=c.get('sign_margin'),
                   maxnref=max([b['nref'] for b in c.get('branches', [])] or [0]))
        return out
    except Exception as ex:                                  # report, never hide
        return dict(exact='error: %s' % str(ex)[:120])


def exact_stage(path, nproc=2, limit=None):
    d = json.load(open(path))
    todo = [r for r in d['rows'] if r.get('finite') and 'CERTIFIED' not in r and 'exact' not in r][:limit]
    with Pool(nproc) as P:
        for r, o in zip(todo, P.imap(exact_row, todo, chunksize=1)):
            r.update(o)
            print(f"  x={r['x']:.10f} n={r.get('n9')} deg={o.get('degree')} unit={o.get('unit')} CERT={o.get('CERTIFIED')} "
                  f"{o.get('exact', '')} minpoly={o.get('minpoly')}", flush=True)
            json.dump(d, open(path, 'w'), default=str)


if __name__ == '__main__' and sys.argv[1] == 'exact':
    exact_stage(sys.argv[2], limit=int(sys.argv[3]) if len(sys.argv) > 3 else None)
