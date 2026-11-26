"""[Deposit] Proposition A.1: the floating-point beam tracer (proposes; the exact code decides). Paper: "Infinitely many perpendicular beams in right triangles at all but countably many irrational angles". Research name: probes/s959_beam_families.py, copied verbatim except imports (deposit names), paths (data/ under beam_note/) and this paragraph. Run through run_all.py.

s959 keeper (user: *the natural question is whether non-right triangles and other polygons show similar [γ = 1] behaviour —
do we have the tools?* → *yes* to building the tracer): N(ε) of the PERPENDICULAR-BEAM partition, for triangle families beyond
the right triangle, all through ONE instrument.

OBJECT (the right-triangle object of `foundations.md` §1, transplanted): launch perpendicular into the table from every point `s`
of the launch side (side 0 = V0V1, length normalised to 1); a BRANCH is a maximal `s`-interval with a constant reflection word
up to the FIRST PERPENDICULAR RETURN to side 0 (head-on arrival); N(ε) = #branches of width > ε.  For the right triangle (right
angle at V0 = R, V1 = A, V2 = O with angle α) side 0 is the launch leg RA, the paper's object.

INSTRUMENT.  Beam tracing in the folded table: a piece = an `s`-interval of parallel rays with ONE direction, carried as the affine
transverse map `t = σ s + t0` (`t = X·d⊥`; a reflection is an isometry, so `σ = ±1` and widths in `s` never change).  A piece splits
only at the vertex opposite its entry side.  The direction angle is tracked SYMBOLICALLY, `φ = ε·π/2 + rπ/D + Aα + Bβ`
(reflection off a side of line angle ψ: `φ ↦ 2ψ − φ`), so the head-on test on side 0 (`φ ≡ −π/2 mod 2π`) is an INTEGER test at
generic (ℚ-independent) angles, and the float direction never drifts.  Pieces narrower than `EPS_MIN` are dropped (their mass
recorded); a piece still running after `CAP` reflections is recorded as CAPPED (mass ⊕ width).  Float64 positions.

FAMILIES (one free angle unless stated; V2 is the vertex opposite the launch side, as O is in the right triangle):
  right  : (V0, V1, V2) = (π/2, π/2 − α, α)                                 — the control (the paper's object)
  pi3    : (π/3, 2π/3 − α, α)       pi5 : (π/5, 4π/5 − α, α)                — one rational angle, the level walk is 1-D
  x2x    : (2x, π − 3x, x)          x3x : (3x, π − 4x, x)                    — every angle in ℤx + ℚπ, no rational angle
  gen    : (β, π − α − β, α), α, β independent                              — TWO free angles, the direction walk is 2-D

PRE-REGISTERED (s959, before any family ran):
  (V0) VALIDATION — right triangle at rational centres (α = Pπ/2Q): the instrument's branch count equals the exact
       `flow_decomp.n_flow(P, Q)` on every centre tested, with `0` capped and `0` dropped mass at `EPS_MIN` below the narrowest branch.
       If (V0) fails nothing else is read.
  (H1) the one-free-angle families (pi3, pi5, x2x, x3x) track the right-triangle control: at matched ε, N(ε)/log(1/ε) stays within
       the spread of the right-triangle parameter sample.
  (H2) the two-free-angle family (gen) departs from the control: N(ε) larger at small ε, or a larger capped mass at the same CAP.
  ⚠ KNOWN TRAP ([§96-NEG], `computational_findings.md` §96 / line ~5215): N(ε) is convex over reachable decades — so NO absolute
       slope is read; only the cross-family comparison at matched depth, same instrument, same parameter count.
PRIOR ART: `future_directions.md` C1 (γ = 1 off the perpendicular — untouched), C2 (obtuse / non-right — untouched), C3 (two free
  angles ⇒ surfaces, PARKED); `rulings.py --grep` non-right, general triangle, obtuse, acute triangle, gamma → nothing on N(ε)
  beyond right triangles ([CPRIME-026]: rooms/doors do not exist off the right triangle); `tools_index.py --grep` irrational,
  polygon, general triangle → right-triangle irrational instruments only (fan IET, s660_deep, T_G) ⊕ rational-angle exact tools
  (flow_decomp, flow_poly, poly_exact).
DATA: `data_index.py --grep` beam / N_eps → right-triangle stores only.
RESULT (s959; data/beam_families.json, data/s959_beam_families_strat.json; numba kernel `trace_nb` = widest-first heap):
  (V0) HOLDS — 8/8 rational right centres n = n_flow, widths to 5e-14; heap kernel = python tracer on 3 irrational tables;
  independent single-ray tracer agrees on 303 (pi3) ⊕ 596 (x2x) random s ⊕ 290 gen branch interiors, 0 mismatches.
  CONTROL (added s959): the ISOSCELES (π/2−α, π/2−α, 2α) from its base — the non-right code path on a table that folds to the
  right triangle — gives N_iso(ε) = N_right(ε) at every ε on all 4 irrational α (exact partition equality at 4/9, 8/21; at odd
  P the orphan splits into two mirror halves).
  ⚠ A FIRST (depth-first) kernel left wide sibling pieces unprocessed when the step budget ran out and booked them CAPPED (pi3
  α = 1.761: 96 % 'capped', widest 0.83); its gen 'first look' (96 % not returned) was the same artefact — withdrawn.
  (H1) REFUTED at this depth: the one-free-angle families do NOT track the right control — N(1e-5) right 18–46; pi3 206–2743,
  pi5 72–1606, x3x 257–533 (one row 1), x2x 3–69; and the growth is multiplicative (×2–3 per decade) where the right
  triangle's is additive.  Exceptions are FINITE partitions at irrational parameter (x2x α = 0.064, 0.22, 0.226: N ≤ 7 for
  every ε; x3x α = 0.372: ONE branch, the whole side returns in 8 reflections).
  (H2) HOLDS in its mass form: within 10⁹ steps the gen tables return 23–88 % of the launch side, every one-free-angle table
  ≥ 95.5 %; N(ε ≥ 1e-4) of gen is NOT larger than pi3/pi5's.  Strata (acute / obtuse): see the strat file and §s959b.
Usage: python probes/beam_families.py validate | run [EPS_MIN] [CAP] [NPAR] → data/beam_families.json
"""
import _setup  # noqa: F401  (deposit: this directory on the path, cwd = beam_note/)
import json
import math
import random
import sys
import time
from fractions import Fraction

OUT = 'data/beam_families.json'
TWO_PI = 2 * math.pi


class Table:
    """Triangle V0 V1 V2, side 0 = V0V1 on the x-axis (V0 at 0, V1 at 1), V2 above.  Angles given SYMBOLICALLY as
    (Fraction q, m, n) meaning qπ + mα + nβ; `al`, `be` floats."""

    def __init__(self, angs, al, be=0.0, D=None):
        self.al, self.be = al, be
        self.sym = angs
        val = [float(q) * math.pi + m * al + n * be for (q, m, n) in angs]
        assert all(v > 0 for v in val) and abs(sum(val) - math.pi) < 1e-12, val
        a0, a1, a2 = val
        # V2 from the law of sines, side 0 = 1
        L2 = math.sin(a1) / math.sin(a2)                    # |V0V2|
        self.V = [(0.0, 0.0), (1.0, 0.0), (L2 * math.cos(a0), L2 * math.sin(a0))]
        self.sides = [(0, 1), (1, 2), (2, 0)]
        # line angles ψ_j mod π, symbolic: side0 = 0, side1 = V1V2 = π − a1, side2 = V2V0 = a0
        (q0, m0, n0), (q1, m1, n1), _ = angs
        self.psi = [(Fraction(0), 0, 0), (Fraction(1) - q1, -m1, -n1), (q0, m0, n0)]
        dens = [Fraction(p[0]).denominator for p in self.psi] + [2]
        self.D = D or math.lcm(*dens)
        self.psiI = [(int(p[0] * self.D), p[1], p[2]) for p in self.psi]   # rπ/D, integer r
        self.psiF = [float(p[0]) * math.pi + p[1] * al + p[2] * be for p in self.psi]

    def phi_float(self, st):
        e, r, A, B = st
        return e * math.pi / 2 + r * math.pi / self.D + A * self.al + B * self.be

    def headon0(self, st):
        e, r, A, B = st
        if A or B:
            return False
        return (r % (2 * self.D)) == (self.D if e == 1 else 0)

    def reflect(self, st, j):
        e, r, A, B = st
        pr, pA, pB = self.psiI[j]
        return (-e, (2 * pr - r) % (2 * self.D), 2 * pA - A, 2 * pB - B)


def trace(T, eps_min=1e-6, cap=10 ** 6, budget=3 * 10 ** 8):
    """All branches of the perpendicular beam from side 0 with width >= eps_min. Returns dict."""
    V = T.V
    branches, capped = [], []
    dropped = 0.0
    steps = 0
    # piece: (s0, s1, sigma, t0, entry side j, symbolic direction, nref)
    st0 = (1, 0, 0, 0)                                  # φ = π/2
    stack = [(0.0, 1.0, 1, 0.0, 0, st0, 0)]
    # at launch d = (0,1), d⊥ := (1,0) (so t = x on side 0, and t = s): σ = 1, t0 = 0
    while stack:
        s0, s1, sg, t0, j, st, nref = stack.pop()
        while True:
            if s1 - s0 < eps_min:
                dropped += s1 - s0
                break
            if nref >= cap or steps >= budget:
                capped.append((s0, s1, nref))
                break
            phi = T.phi_float(st)
            dx, dy = math.cos(phi), math.sin(phi)
            px, py = dy, -dx                            # d⊥ (fixed convention: rotate d by −90°)
            a, b = T.sides[j]
            o = 3 - a - b
            tv = lambda k: V[k][0] * px + V[k][1] * py
            ta, tb, to = tv(a), tv(b), tv(o)
            so = (to - t0) * sg                         # s of the ray through the opposite vertex
            # rays with t between ta and to exit through side (a, o); between to and tb through (o, b)
            sa, sb = (ta - t0) * sg, (tb - t0) * sg
            parts = []
            for (u, w, k) in ((min(sa, so), max(sa, so), (a, o)), (min(so, sb), max(so, sb), (o, b))):
                lo, hi = max(s0, u), min(s1, w)
                if hi > lo:
                    parts.append((lo, hi, k))
            steps += 1
            if not parts:                               # numerical sliver
                dropped += s1 - s0
                break
            # continue with the first part, push the rest
            for (lo, hi, k) in parts[1:]:
                pc = _exit(T, lo, hi, sg, t0, k, st, nref, px, py, branches)
                if pc is not None:
                    stack.append(pc)
            nxt = _exit(T, parts[0][0], parts[0][1], sg, t0, parts[0][2], st, nref, px, py, branches)
            if nxt is None:
                break
            s0, s1, sg, t0, j, st, nref = nxt
    return {'branches': branches, 'capped': capped, 'dropped': dropped, 'steps': steps}


def _exit(T, lo, hi, sg, t0, k, st, nref, px, py, branches):
    """The rays s ∈ [lo, hi] leave through side k = (a, b) (vertex indices); reflect, or record a branch on head-on side 0."""
    a, b = k
    j = T.sides.index((a, b)) if (a, b) in T.sides else T.sides.index((b, a))
    if j == 0 and T.headon0(st):
        branches.append((lo, hi, nref + 1))
        return None
    st2 = T.reflect(st, j)
    phi2 = T.phi_float(st2)
    dx2, dy2 = math.cos(phi2), math.sin(phi2)
    px2, py2 = dy2, -dx2
    # transverse map after reflection: for X on side j, t'(X) affine in t(X) with slope ±1; fit on the side's endpoints
    A_, B_ = T.V[T.sides[j][0]], T.V[T.sides[j][1]]
    tA = A_[0] * px + A_[1] * py
    tB = B_[0] * px + B_[1] * py
    uA = A_[0] * px2 + A_[1] * py2
    uB = B_[0] * px2 + B_[1] * py2
    slope = (uB - uA) / (tB - tA)
    s2 = 1 if slope > 0 else -1
    # t = sg*s + t0 ; t' = s2*(t - tA) + uA = s2*sg*s + (s2*(t0 - tA) + uA)
    return (lo, hi, s2 * sg, s2 * (t0 - tA) + uA, j, st2, nref + 1)


FR = Fraction
FAMILIES = {
    'right': lambda: ([(FR(1, 2), 0, 0), (FR(1, 2), -1, 0), (FR(0), 1, 0)], (0.05, 1.5)),
    'pi3': lambda: ([(FR(1, 3), 0, 0), (FR(2, 3), -1, 0), (FR(0), 1, 0)], (0.05, 2.0)),
    'pi5': lambda: ([(FR(1, 5), 0, 0), (FR(4, 5), -1, 0), (FR(0), 1, 0)], (0.05, 2.4)),
    'x2x': lambda: ([(FR(0), 2, 0), (FR(1), -3, 0), (FR(0), 1, 0)], (0.05, 1.0)),
    'x3x': lambda: ([(FR(0), 3, 0), (FR(1), -4, 0), (FR(0), 1, 0)], (0.05, 0.75)),
    'gen': lambda: ([(FR(0), 0, 1), (FR(1), -1, -1), (FR(0), 1, 0)], None),
}


def n_of_eps(br, grid):
    w = sorted((b[1] - b[0] for b in br), reverse=True)
    out, i = [], 0
    for e in grid:
        while i < len(w) and w[i] > e:
            i += 1
        out.append(i)
    return out


def validate():
    import flow_decomp as FD
    res = []
    for P, Q in [(1, 3), (2, 5), (3, 7), (4, 9), (5, 13), (3, 11), (7, 17), (8, 21)]:
        al = P * math.pi / (2 * Q)
        # rational α: fold into the rational part (exact head-on test)
        angs = [(FR(1, 2), 0, 0), (FR(1, 2) - FR(P, 2 * Q), 0, 0), (FR(P, 2 * Q), 0, 0)]
        T = Table(angs, 0.0)
        t = time.time()
        r = trace(T, eps_min=1e-12, cap=10 ** 7)
        nf = FD.n_flow(P, Q)
        wt = sorted(b[1] - b[0] for b in r['branches'])
        wf = sorted(float(x) for x in nf['widths'])
        wdiff = max(abs(x - y) for x, y in zip(wt, wf)) if len(wt) == len(wf) else None
        row = dict(P=P, Q=Q, n_tracer=len(r['branches']), n_flow=nf['n'], width_maxdiff=wdiff, capped=len(r['capped']), dropped=r['dropped'],
                   minw=min(b[1] - b[0] for b in r['branches']), mass=sum(b[1] - b[0] for b in r['branches']),
                   secs=round(time.time() - t, 2))
        print(row, flush=True)
        res.append(row)
    return res


def run(eps_min, cap, npar, seed=959, budget=10 ** 9):
    rnd = random.Random(seed)
    grid = [10 ** (-x / 4) for x in range(4, int(-4 * math.log10(eps_min)) + 1)]
    out = {'eps_min': eps_min, 'cap': cap, 'budget': budget, 'grid': grid, 'rows': []}
    for fam, mk in FAMILIES.items():
        angs, rng = mk()
        for i in range(npar):
            if fam == 'gen':
                while True:
                    al, be = rnd.uniform(0.05, 1.5), rnd.uniform(0.05, 1.5)
                    if al + be < math.pi - 0.05:
                        break
            else:
                al, be = rnd.uniform(*rng), 0.0
            T = Table(angs, al, be)
            t = time.time()
            r = trace_nb(T, eps_min=eps_min, cap=cap, budget=budget)
            row = dict(fam=fam, al=al, be=be, N=n_of_eps(r['branches'], grid), nbr=len(r['branches']),
                       mass=sum(b[1] - b[0] for b in r['branches']), capped_mass=sum(c[1] - c[0] for c in r['capped']),
                       ncapped=len(r['capped']), maxcapw=max((c[1] - c[0] for c in r['capped']), default=0.0), dropped=r['dropped'], steps=r['steps'], secs=round(time.time() - t, 1),
                       maxref=max((b[2] for b in r['branches']), default=0))
            out['rows'].append(row)
            print(fam, round(al, 4), round(be, 4), 'N(1e-2,1e-3,..)=', [row['N'][k] for k in range(4, len(grid), 4)],
                  'capped', round(row['capped_mass'], 4), 'maxcapw', '%.2g' % row['maxcapw'], 'dropped', '%.2g' % row['dropped'], 'steps', row['steps'], row['secs'], 's', flush=True)
            json.dump(out, open(OUT, 'w'))
    return out


def ray_check(T, s, maxref=10 ** 6):
    """Naive single-ray tracer (independent code): reflections until head-on arrival on side 0, by the symbolic direction."""
    x, y = s, 0.0
    st = (1, 0, 0, 0)
    j = 0
    for n in range(1, maxref + 1):
        phi = T.phi_float(st)
        dx, dy = math.cos(phi), math.sin(phi)
        best = None
        for k, (a, b) in enumerate(T.sides):
            if k == j:
                continue
            (ax, ay), (bx, by) = T.V[a], T.V[b]
            ex, ey = bx - ax, by - ay
            den = dx * ey - dy * ex
            if abs(den) < 1e-300:
                continue
            tt = ((ax - x) * ey - (ay - y) * ex) / den
            uu = ((ax - x) * dy - (ay - y) * dx) / den
            if tt > 1e-12 and -1e-12 <= uu <= 1 + 1e-12 and (best is None or tt < best[0]):
                best = (tt, k)
        tt, k = best
        x, y = x + tt * dx, y + tt * dy
        if k == 0 and T.headon0(st):
            return n
        st = T.reflect(st, k)
        j = k
    return None


def check_rays(nsample=300, seed=1):
    rnd = random.Random(seed)
    bad = tot = 0
    for fam in ('gen', 'pi3', 'x2x'):
        angs, rng = FAMILIES[fam]()
        al, be = (0.7123, 0.9311) if fam == 'gen' else (rnd.uniform(*rng), 0.0)
        T = Table(angs, al, be)
        r = trace_nb(T, eps_min=1e-5, cap=10 ** 7, budget=3 * 10 ** 8)
        br = sorted(r['branches'])
        for _ in range(nsample):
            s = rnd.random()
            hit = [b for b in br if b[0] < s < b[1] and b[1] - b[0] > 1e-4]
            if not hit:
                continue
            n = ray_check(T, s)
            tot += 1
            if n != hit[0][2]:
                bad += 1
                print('MISMATCH', fam, s, n, hit[0])
        print(fam, 'checked', tot, 'bad', bad, flush=True)
    return bad, tot


# ------------------------------------------------------------------------------------------------ numba kernel (same algorithm)
def _nb_kernel():
    import numba as nb
    import numpy as np

    @nb.njit(cache=True)
    def kern(V, sides_a, sides_b, psir, psiA, psiB, D, al, be, eps_min, cap, budget, maxst):
        # stack of pieces
        S0 = np.empty(maxst); S1 = np.empty(maxst); SG = np.empty(maxst, np.int64); T0 = np.empty(maxst)
        J = np.empty(maxst, np.int64); E = np.empty(maxst, np.int64); R = np.empty(maxst, np.int64)
        AA = np.empty(maxst, np.int64); BB = np.empty(maxst, np.int64); NR = np.empty(maxst, np.int64)
        nbmax = 2000000
        bl = np.empty(nbmax); bh = np.empty(nbmax); bn = np.empty(nbmax, np.int64)
        cl = np.empty(nbmax); ch = np.empty(nbmax)
        nbr = 0; ncap = 0; dropped = 0.0; steps = 0; top = 0
        S0[0] = 0.0; S1[0] = 1.0; SG[0] = 1; T0[0] = 0.0; J[0] = 0; E[0] = 1; R[0] = 0; AA[0] = 0; BB[0] = 0; NR[0] = 0
        top = 1
        twoD = 2 * D
        overflow = 0
        while top > 0:
            top -= 1
            s0 = S0[top]; s1 = S1[top]; sg = SG[top]; t0 = T0[top]; j = J[top]
            e = E[top]; r = R[top]; A = AA[top]; B = BB[top]; nref = NR[top]
            while True:
                if s1 - s0 < eps_min:
                    dropped += s1 - s0
                    break
                if nref >= cap or steps >= budget:
                    if ncap < nbmax:
                        cl[ncap] = s0; ch[ncap] = s1
                    ncap += 1
                    break
                phi = e * np.pi / 2 + r * np.pi / D + A * al + B * be
                dx = np.cos(phi); dy = np.sin(phi)
                px = dy; py = -dx
                a = sides_a[j]; b = sides_b[j]; o = 3 - a - b
                ta = V[a, 0] * px + V[a, 1] * py
                tb = V[b, 0] * px + V[b, 1] * py
                to = V[o, 0] * px + V[o, 1] * py
                so = (to - t0) * sg; sa = (ta - t0) * sg; sb = (tb - t0) * sg
                steps += 1
                np_ = 0
                plo = np.empty(2); phi_ = np.empty(2); pk1 = np.empty(2, np.int64); pk2 = np.empty(2, np.int64)
                for q in range(2):
                    if q == 0:
                        u = min(sa, so); w = max(sa, so); k1 = a; k2 = o
                    else:
                        u = min(so, sb); w = max(so, sb); k1 = o; k2 = b
                    lo = max(s0, u); hi = min(s1, w)
                    if hi > lo:
                        plo[np_] = lo; phi_[np_] = hi; pk1[np_] = k1; pk2[np_] = k2; np_ += 1
                if np_ == 0:
                    dropped += s1 - s0
                    break
                cont = False
                for q in range(np_ - 1, -1, -1):
                    lo = plo[q]; hi = phi_[q]; k1 = pk1[q]; k2 = pk2[q]
                    # side index of (k1, k2)
                    jj = 0
                    for z in range(3):
                        if (sides_a[z] == k1 and sides_b[z] == k2) or (sides_a[z] == k2 and sides_b[z] == k1):
                            jj = z
                    if jj == 0 and A == 0 and B == 0 and ((r % twoD + twoD) % twoD) == (D if e == 1 else 0):
                        if nbr < nbmax:
                            bl[nbr] = lo; bh[nbr] = hi; bn[nbr] = nref + 1
                        nbr += 1
                        continue
                    e2 = -e; r2 = ((2 * psir[jj] - r) % twoD + twoD) % twoD; A2 = 2 * psiA[jj] - A; B2 = 2 * psiB[jj] - B
                    phi2 = e2 * np.pi / 2 + r2 * np.pi / D + A2 * al + B2 * be
                    px2 = np.sin(phi2); py2 = -np.cos(phi2)
                    ia = sides_a[jj]; ib = sides_b[jj]
                    tA = V[ia, 0] * px + V[ia, 1] * py; tB = V[ib, 0] * px + V[ib, 1] * py
                    uA = V[ia, 0] * px2 + V[ia, 1] * py2; uB = V[ib, 0] * px2 + V[ib, 1] * py2
                    s2 = 1 if (uB - uA) / (tB - tA) > 0 else -1
                    nsg = s2 * sg; nt0 = s2 * (t0 - tA) + uA
                    if q == 0:
                        s0 = lo; s1 = hi; sg = nsg; t0 = nt0; j = jj; e = e2; r = r2; A = A2; B = B2; nref = nref + 1
                        cont = True
                    else:
                        if top >= maxst:
                            overflow += 1
                            dropped += hi - lo
                        else:
                            S0[top] = lo; S1[top] = hi; SG[top] = nsg; T0[top] = nt0; J[top] = jj; E[top] = e2
                            R[top] = r2; AA[top] = A2; BB[top] = B2; NR[top] = nref + 1
                            top += 1
                if not cont:
                    break
        return bl[:min(nbr, nbmax)].copy(), bh[:min(nbr, nbmax)].copy(), bn[:min(nbr, nbmax)].copy(), \
            cl[:min(ncap, nbmax)].copy(), ch[:min(ncap, nbmax)].copy(), dropped, steps, overflow

    @nb.njit(cache=True)
    def kern_heap(V, sides_a, sides_b, psir, psiA, psiB, D, al, be, eps_min, cap, budget, maxst):
        """WIDEST-FIRST: pieces live in slots; a binary max-heap on width picks the next piece; a piece runs until it splits
        (both children go back to the heap) or ends.  A budget cut therefore leaves only pieces narrower than every
        processed one: N(eps) is exact for eps >= the widest capped piece."""
        S0 = np.empty(maxst); S1 = np.empty(maxst); SG = np.empty(maxst, np.int64); T0 = np.empty(maxst)
        J = np.empty(maxst, np.int64); E = np.empty(maxst, np.int64); R = np.empty(maxst, np.int64)
        AA = np.empty(maxst, np.int64); BB = np.empty(maxst, np.int64); NR = np.empty(maxst, np.int64)
        heap = np.empty(maxst, np.int64); hn = 0
        free = np.empty(maxst, np.int64); nfree = 0; nslot = 0
        nbmax = 4000000
        bl = np.empty(nbmax); bh = np.empty(nbmax); bn = np.empty(nbmax, np.int64)
        cl = np.empty(nbmax); ch = np.empty(nbmax); cn = np.empty(nbmax, np.int64)
        nbr = 0; ncap = 0; dropped = 0.0; steps = 0; overflow = 0
        twoD = 2 * D
        # slot 0 = the launch
        S0[0] = 0.0; S1[0] = 1.0; SG[0] = 1; T0[0] = 0.0; J[0] = 0; E[0] = 1; R[0] = 0; AA[0] = 0; BB[0] = 0; NR[0] = 0
        nslot = 1; heap[0] = 0; hn = 1
        plo = np.empty(2); phi_ = np.empty(2); pk1 = np.empty(2, np.int64); pk2 = np.empty(2, np.int64)
        while hn > 0:
            # pop max
            sl = heap[0]
            hn -= 1
            if hn > 0:
                heap[0] = heap[hn]
                i = 0
                while True:
                    l = 2 * i + 1; rr = l + 1; m = i
                    if l < hn and S1[heap[l]] - S0[heap[l]] > S1[heap[m]] - S0[heap[m]]:
                        m = l
                    if rr < hn and S1[heap[rr]] - S0[heap[rr]] > S1[heap[m]] - S0[heap[m]]:
                        m = rr
                    if m == i:
                        break
                    tmp = heap[i]; heap[i] = heap[m]; heap[m] = tmp; i = m
            s0 = S0[sl]; s1 = S1[sl]; sg = SG[sl]; t0 = T0[sl]; j = J[sl]
            e = E[sl]; r = R[sl]; A = AA[sl]; B = BB[sl]; nref = NR[sl]
            free[nfree] = sl; nfree += 1
            if steps >= budget:
                # budget gone: everything left is capped, untouched
                if ncap < nbmax:
                    cl[ncap] = s0; ch[ncap] = s1; cn[ncap] = nref
                ncap += 1
                continue
            while True:
                if s1 - s0 < eps_min:
                    dropped += s1 - s0
                    break
                if nref >= cap or steps >= budget:
                    if ncap < nbmax:
                        cl[ncap] = s0; ch[ncap] = s1; cn[ncap] = nref
                    ncap += 1
                    break
                phi = e * np.pi / 2 + r * np.pi / D + A * al + B * be
                dx = np.cos(phi); dy = np.sin(phi)
                px = dy; py = -dx
                a = sides_a[j]; b = sides_b[j]; o = 3 - a - b
                ta = V[a, 0] * px + V[a, 1] * py
                tb = V[b, 0] * px + V[b, 1] * py
                to = V[o, 0] * px + V[o, 1] * py
                so = (to - t0) * sg; sa = (ta - t0) * sg; sb = (tb - t0) * sg
                steps += 1
                np_ = 0
                for q in range(2):
                    if q == 0:
                        u = min(sa, so); w = max(sa, so); k1 = a; k2 = o
                    else:
                        u = min(so, sb); w = max(so, sb); k1 = o; k2 = b
                    lo = max(s0, u); hi = min(s1, w)
                    if hi > lo:
                        plo[np_] = lo; phi_[np_] = hi; pk1[np_] = k1; pk2[np_] = k2; np_ += 1
                if np_ == 0:
                    dropped += s1 - s0
                    break
                cont = False
                for q in range(np_):
                    lo = plo[q]; hi = phi_[q]; k1 = pk1[q]; k2 = pk2[q]
                    jj = 0
                    for z in range(3):
                        if (sides_a[z] == k1 and sides_b[z] == k2) or (sides_a[z] == k2 and sides_b[z] == k1):
                            jj = z
                    if jj == 0 and A == 0 and B == 0 and ((r % twoD + twoD) % twoD) == (D if e == 1 else 0):
                        if nbr < nbmax:
                            bl[nbr] = lo; bh[nbr] = hi; bn[nbr] = nref + 1
                        nbr += 1
                        continue
                    e2 = -e; r2 = ((2 * psir[jj] - r) % twoD + twoD) % twoD; A2 = 2 * psiA[jj] - A; B2 = 2 * psiB[jj] - B
                    phi2 = e2 * np.pi / 2 + r2 * np.pi / D + A2 * al + B2 * be
                    px2 = np.sin(phi2); py2 = -np.cos(phi2)
                    ia = sides_a[jj]; ib = sides_b[jj]
                    tA = V[ia, 0] * px + V[ia, 1] * py; tB = V[ib, 0] * px + V[ib, 1] * py
                    uA = V[ia, 0] * px2 + V[ia, 1] * py2; uB = V[ib, 0] * px2 + V[ib, 1] * py2
                    s2 = 1 if (uB - uA) / (tB - tA) > 0 else -1
                    nsg = s2 * sg; nt0 = s2 * (t0 - tA) + uA
                    if np_ == 1:
                        s0 = lo; s1 = hi; sg = nsg; t0 = nt0; j = jj; e = e2; r = r2; A = A2; B = B2; nref = nref + 1
                        cont = True
                    else:
                        if hi - lo < eps_min:
                            dropped += hi - lo
                            continue
                        if nfree > 0:
                            nfree -= 1; ns = free[nfree]
                        elif nslot < maxst:
                            ns = nslot; nslot += 1
                        else:
                            overflow += 1; dropped += hi - lo
                            continue
                        S0[ns] = lo; S1[ns] = hi; SG[ns] = nsg; T0[ns] = nt0; J[ns] = jj; E[ns] = e2
                        R[ns] = r2; AA[ns] = A2; BB[ns] = B2; NR[ns] = nref + 1
                        # push
                        i = hn; heap[hn] = ns; hn += 1
                        while i > 0:
                            pa = (i - 1) // 2
                            if S1[heap[pa]] - S0[heap[pa]] >= S1[heap[i]] - S0[heap[i]]:
                                break
                            tmp = heap[i]; heap[i] = heap[pa]; heap[pa] = tmp; i = pa
                if not cont:
                    break
        return bl[:min(nbr, nbmax)].copy(), bh[:min(nbr, nbmax)].copy(), bn[:min(nbr, nbmax)].copy(), \
            cl[:min(ncap, nbmax)].copy(), ch[:min(ncap, nbmax)].copy(), dropped, steps, overflow, cn[:min(ncap, nbmax)].copy()

    return kern, kern_heap


_KERN = None


def trace_nb(T, eps_min=1e-6, cap=10 ** 7, budget=10 ** 10, maxst=1 << 22):
    global _KERN
    import numpy as np
    if _KERN is None:
        _KERN = _nb_kernel()
    V = np.array(T.V, dtype=np.float64)
    sa = np.array([a for a, b in T.sides], np.int64); sb = np.array([b for a, b in T.sides], np.int64)
    pr = np.array([p[0] for p in T.psiI], np.int64); pA = np.array([p[1] for p in T.psiI], np.int64)
    pB = np.array([p[2] for p in T.psiI], np.int64)
    bl, bh, bn, cl, ch, dropped, steps, ovf, cn = _KERN[1](V, sa, sb, pr, pA, pB, T.D, T.al, T.be, eps_min, cap, budget, maxst)
    assert ovf == 0, 'slot overflow'
    return {'branches': list(zip(bl.tolist(), bh.tolist(), bn.tolist())), 'capped': list(zip(cl.tolist(), ch.tolist(), cn.tolist())),
            'dropped': float(dropped), 'steps': int(steps)}

def run_strat(eps_min, cap, npar, budget=10 ** 9, seed=9592, out_path='data/s959_beam_families_strat.json'):
    """Follow-up: ACUTE vs OBTUSE strata per family (npar each), same instrument and depth."""
    rnd = random.Random(seed)
    grid = [10 ** (-x / 4) for x in range(4, int(-4 * math.log10(eps_min)) + 1)]
    out = {'eps_min': eps_min, 'cap': cap, 'budget': budget, 'grid': grid, 'rows': []}
    for fam in ('pi3', 'pi5', 'x2x', 'x3x', 'gen'):
        angs, rng = FAMILIES[fam]()
        for stratum in ('acute', 'obtuse'):
            got = 0
            while got < npar:
                if fam == 'gen':
                    al, be = rnd.uniform(0.05, 2.5), rnd.uniform(0.05, 2.5)
                else:
                    al, be = rnd.uniform(*rng), 0.0
                val = [float(q) * math.pi + m * al + n * be for (q, m, n) in angs]
                if min(val) <= 0.05 or abs(sum(val) - math.pi) > 1e-9:
                    continue
                ob = max(val) > math.pi / 2 + 0.02
                ac = max(val) < math.pi / 2 - 0.02
                if (stratum == 'acute' and not ac) or (stratum == 'obtuse' and not ob):
                    continue
                got += 1
                T = Table(angs, al, be)
                t = time.time()
                r = trace_nb(T, eps_min=eps_min, cap=cap, budget=budget)
                row = dict(fam=fam, stratum=stratum, al=al, be=be, angles=val, N=n_of_eps(r['branches'], grid),
                           mass=sum(b[1] - b[0] for b in r['branches']), capped_mass=sum(c[1] - c[0] for c in r['capped']),
                           maxcapw=max((c[1] - c[0] for c in r['capped']), default=0.0), dropped=r['dropped'],
                           steps=r['steps'], secs=round(time.time() - t, 1))
                out['rows'].append(row)
                print(fam, stratum, round(al, 4), round(be, 4), [row['N'][k] for k in range(4, len(grid), 4)],
                      'ret %.3f maxcapw %.1g' % (row['mass'], row['maxcapw']), row['secs'], 's', flush=True)
                json.dump(out, open(out_path, 'w'))
    return out



if __name__ == '__main__':
    if sys.argv[1] == 'validate':
        validate()
        check_rays()
    elif sys.argv[1] == 'rays':
        check_rays()
    elif sys.argv[1] == 'strat':
        run_strat(float(sys.argv[2]), int(float(sys.argv[3])), int(sys.argv[4]), budget=int(float(sys.argv[5])) if len(sys.argv) > 5 else 2 * 10 ** 8)
    else:
        eps = float(sys.argv[2]) if len(sys.argv) > 2 else 1e-5
        cap = int(float(sys.argv[3])) if len(sys.argv) > 3 else 10 ** 6
        npar = int(sys.argv[4]) if len(sys.argv) > 4 else 4
        run(eps, cap, npar)
