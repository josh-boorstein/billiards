"""[Deposit] Theorem 5.4 at m = 1 (Π_N): the per-x certificate on (0, 0.2095]; the Laurent ring. Paper: "Perpendicular billiards in rational polygons of odd denominator". Research name: probes/s935_pi_n_torus.py, copied verbatim except imports (deposit names), paths (data/ and figs/ under odd_polygons/) and this paragraph. Run through run_all.py; the usage lines below are the research repo's. Bracketed IDs and section marks (§sNNN) refer to the project's internal ledgers.

s935 keeper (queue `s933-b`): Π_N = (N−4, 4, N−2, 2)/N IS NON-CP FOR EVERY ODD N ≥ 15 — a FORMAL certificate, uniform in N.

NOT PRE-REGISTERED: exploratory. Found by dumping `flow_poly`'s non-periodic components of Π_N (N = 15 … 41): every one ends
as a 2-interval exchange; the first time the leading one has two labels its return times are (24, 18) at EVERY N, through the
same 42 strips (copy offsets relative to K = (N−1)/2). This file re-does that torus with no `flow_decomp` at all:

  FORMAL ARITHMETIC. Elements are Laurent polynomials in z = ζ_{4N} with Z[i] coefficients, using ζ^N = i, ζ^{2N} = −1; copy
  (o, e) means k = K + o, so its rotation ζ^{4k} = −z^{4o−2}. Every exponent that occurs is then N-free, so a formal identity
  holds for EVERY odd N, and a formal inequality is a real trig polynomial f(x), x = π/N.
  GEOMETRY = `flow_poly._build_iet_poly`'s: flow ω = ζ^N = i (⊥ side V0V1), transverse coordinate X(p) = −2 Re p, a side is an
  entry side of copy g iff sign(2 Re(side vector)) = orientation of g, crossing exit side s2 moves to g_across(g, s2) and
  shifts X by −X(g·t_{s2}).
  CERTIFICATE. J = [lo0, hi1) in copy (1, 1) on entry side 3, split at hi0; lo0, hi0, hi1 pinned by vertex coincidences
  (step 5 / 16 of J0, step 7 of J1). Checked: (1) Π_N closes formally; (2) the itinerary is combinatorially consistent and
  returns on side 3 of copy (1, 1); (3) T^24 J0 = [lo0 + λ1, hi1) and T^18 J1 = [lo0, lo0 + λ1) FORMALLY — J's first return is the
  ROTATION by λ1 mod λ0 + λ1; (4) λ0 = 4 sin4x sinx (2cos2x − 1), λ1 = 16 sin4x sin³x formally, so the rotation number is
  λ1/(λ0+λ1) = 4 sin²(π/N), irrational for every N ∉ {1, 2, 3, 4, 6}; (5) every strip lies in its entry AND exit side's
  X-range and those sides are entry/exit (254 conditions = 236 certified + 18 formally tight): f(x) ≥ 0 on x ∈ (0, π/15] RIGOROUSLY —
  exact rational Taylor coefficients ENCLOSED in iv, f/x^{n0} by interval Horner on 400 subintervals, Lagrange tail bound.
  First return: no intermediate strip enters copy (1, 1) through side 3 (asserted); only a return map is needed anyway —
  T-periodic ⟹ the return map periodic, impossible for an irrational rotation.
  ⟹ a point of J off the countable set of separatrix points has an infinite, vertex-free, non-periodic orbit: the
  perpendicular class of Π_N is not completely periodic, every odd N ≥ 15.
CONTROLS: N = 13 (Π_13 is CP, `data/s933_regular_pieces.json`): 19 inequalities FAIL (36, 64, 77 at N = 11, 9, 7); the
  certificate's x-range ends exactly at x = π/14. `--xcheck N`: `flow_poly` at a held-out N, every non-periodic component's
  rotation number at its 2-label snapshot equals 4 sin²(π/N) (independent of the itinerary; same flow_poly geometry —
  the geometry-free check is the audit's plain billiard simulation, `audits/s935_audit_torus.md`).
COULD HAVE COME OUT OTHERWISE: the return identities could hold only modulo Φ_{4N} for the measured N (they hold formally);
  an inequality could change sign inside (π/101, π/15] or as x → 0 (none does; leading orders 0, 2, 4, all positive).
PRIOR ART: §s933 (4) ([POLY-ASYM-NONCP], per-row certificates 15 ≤ N ≤ 101, "a proof for every odd N ≥ 15" listed as not done);
  [NCYL-721] ↓s933; §s924 (the formal-identity ⊕ Taylor ⊕ interval template, for triangles). Grep of the closed form
  `4sin²(π/N)` / `sin²(π/N)` / `2 − 2cos(2π/N)` → only `w(3/Q) = 1/(3 − 4sin²(π/2Q))` (computational_findings, orphan widths —
  another object); `rulings.py --grep` *rotation number* → none on polygons.
DATA: `data_index.py --grep torus` → none. Output `data/s935_pi_n_torus.json`.
→ general_triangles.md §s935.   Usage: python probes/s935_pi_n_torus.py [--xcheck N]   (≈ 1 min; --xcheck N=61 ≈ minutes)
"""
import _setup  # noqa: F401  (deposit: this directory on the path, cwd = odd_polygons/)
import cmath
import json
import math
import sys
from fractions import Fraction as Fr

import mpmath as mp
from mpmath import iv


class L:
    """Laurent polynomial in z = ζ_{4N}, coefficients Gaussian integers (stored as complex of small ints — exact)"""
    __slots__ = ('c',)

    def __init__(self, c=None):
        self.c = {k: v for k, v in (c or {}).items() if v != 0}

    @staticmethod
    def mono(b, coef=1):
        return L({b: complex(coef)})

    def __add__(s, o):
        o = o if isinstance(o, L) else L.mono(0, o)
        c = dict(s.c)
        for k, v in o.c.items():
            c[k] = c.get(k, 0) + v
        return L(c)
    __radd__ = __add__

    def __neg__(s):
        return L({k: -v for k, v in s.c.items()})

    def __sub__(s, o):
        return s + (-o)

    def __mul__(s, o):
        o = o if isinstance(o, L) else L.mono(0, o)
        c = {}
        for a, u in s.c.items():
            for b, v in o.c.items():
                c[a + b] = c.get(a + b, 0) + u * v
        return L(c)
    __rmul__ = __mul__

    def conj(s):
        return L({-k: v.conjugate() for k, v in s.c.items()})

    def evx(s, x):
        return sum(v * cmath.exp(1j * k * x / 2) for k, v in s.c.items())

    def iszero(s):
        return not s.c


z = L.mono
c2 = z(4) + z(-4)                                   # 2cos(2π/N)
LEN = [z(0), c2 - 1, c2, c2]                        # |V0V1| : |V1V2| : |V2V3| : |V3V0| = 1 : 2c−1 : 2c : 2c
DIR = [z(0), -z(-8), -z(-4), z(-8)]                 # ζ^{2 m_full}: m_full = 0, N−4, N−2, 2N−4 (units π/N)
REF = [z(0), z(-16), z(-8), z(-16)]                 # ζ^{4m}, m = m_full mod N
DM = [0, -4, -2, -4]                                # m mod N as an offset of K (g_across)
V = [L()]
for _s in range(3):
    V.append(V[-1] + LEN[_s] * DIR[_s])
SIDES = [(s, (s + 1) % 4) for s in range(4)]
T = [V[i] - REF[s] * V[i].conj() for s, (i, _) in enumerate(SIDES)]


def g_apply(g, p):
    o, e = g
    return -z(4 * o - 2) * (p.conj() if e else p)


def g_across(g, s):
    o, e = g
    return (o + DM[s], 1) if e == 0 else (o - DM[s], 0)


def X(p):
    return -(p + p.conj())


def side_range(g, s):
    i, k = SIDES[s]
    return X(g_apply(g, V[i])), X(g_apply(g, V[k]))


def shift(g, s2):
    return X(g_apply(g, T[s2]))


# the itinerary (o, e, entry side, exit side), read off flow_poly at N = 17 and identical for N = 17 … 41 (s935)
J0 = [(1, 1, 3, 2), (3, 0, 2, 3), (-1, 1, 3, 2), (1, 0, 2, 3), (-3, 1, 3, 2), (-1, 0, 2, 0), (-1, 1, 0, 1), (3, 0, 1, 3),
      (-1, 1, 3, 2), (1, 0, 2, 0), (1, 1, 0, 2), (3, 0, 2, 3), (-1, 1, 3, 2), (1, 0, 2, 3), (-3, 1, 3, 1), (1, 0, 1, 0),
      (1, 1, 0, 2), (3, 0, 2, 3), (-1, 1, 3, 2), (1, 0, 2, 3), (-3, 1, 3, 1), (1, 0, 1, 0), (1, 1, 0, 1), (5, 0, 1, 3)]
J1 = J0[:16] + J0[22:]
PINS = {'lo0': (J0, 5, 0), 'hi0': (J0, 16, 2), 'hi1': (J1, 7, 0)}   # (segment, step, vertex): endpoint = X(g V) + S_step
XMAX = Fr(2095, 10000)                                               # > π/15


def walk(seg):
    S, out = L(), []
    for (o, e, s, s2) in seg:
        out.append(S)
        S = S + shift((o, e), s2)
    return out, S


def pinned(name):
    seg, k, vi = PINS[name]
    Ss, _ = walk(seg)
    o, e, _, _ = seg[k]
    return X(g_apply((o, e), V[vi])) + Ss[k]


def sinL(k):                                           # sin(kx) = (z^{2k} − z^{−2k})/(2i)
    return (z(2 * k) - z(-2 * k)) * (-0.5j)


def taylor(f, M):
    out = []
    for n in range(M + 1):
        re = im = Fr(0)
        for b, v in f.c.items():
            cr, ci = Fr(int(v.real)), Fr(int(v.imag))
            assert cr == v.real and ci == v.imag
            t = Fr(b, 2) ** n / math.factorial(n)
            ir, ii = [(1, 0), (0, 1), (-1, 0), (0, -1)][n % 4]
            re += (cr * ir - ci * ii) * t
            im += (cr * ii + ci * ir) * t
        assert im == 0, 'f is not real'
        out.append(re)
    return out


def certify_nonneg(f, M=70, K=400):
    """f(x) > 0 on (0, XMAX], rigorously; returns (n0, certified lower bound of f/x^n0).  Sound (s935 audit G1): rational
    coefficients ENCLOSED as iv.mpf(num)/iv.mpf(den), subinterval endpoints built in iv; tail Σ_{n>M} |f_n| x^{n−n0} ≤
    Σ_b |c_b| y_b^{M+1} e^{y_b} / ((M+1)! XMAX^{n0}), y_b = |b| XMAX/2 (Lagrange), |c_b| ≤ |Re| + |Im|."""
    Tc = taylor(f, M)
    n0 = next(n for n, v in enumerate(Tc) if v != 0)
    assert Tc[n0] > 0, 'leading Taylor coefficient negative'
    with mp.workdps(50):
        Xm = mp.mpf(XMAX.numerator) / XMAX.denominator
        tail = mp.mpf(0)
        for b, v in f.c.items():
            y = mp.mpf(abs(b)) / 2 * Xm
            tail += (abs(int(v.real)) + abs(int(v.imag))) * y ** (M + 1) / mp.factorial(M + 1) * mp.e ** y
        tail = tail / Xm ** n0 * (1 + mp.mpf(10) ** -20)
    co = [iv.mpf(Tc[n].numerator) / iv.mpf(Tc[n].denominator) for n in range(n0, M + 1)]
    Xiv = iv.mpf(XMAX.numerator) / iv.mpf(XMAX.denominator)
    lowest = None
    for j in range(K):
        Xi = iv.mpf([(Xiv * j / K).a if j else iv.mpf(0).a, (Xiv * (j + 1) / K).b])
        acc = iv.mpf(0)
        for c in reversed(co):
            acc = acc * Xi + c
        lowest = acc.a if lowest is None or acc.a < lowest else lowest
    assert lowest > tail, 'interval bound not positive'
    return n0, float(lowest)


def main():
    iv.dps = 40
    out = {}
    assert (V[-1] + LEN[3] * DIR[3]).iszero(), 'Π_N does not close formally'
    for seg in (J0, J1):
        for a, b in zip(seg, seg[1:]):
            assert g_across(a[:2], a[3]) == b[:2] and a[3] == b[2], 'itinerary inconsistent'
        assert g_across(seg[-1][:2], seg[-1][3]) == (1, 1) and seg[-1][3] == 3, 'no return to side 3 of copy (1, 1)'
    for seg in (J0, J1):        # FIRST return (audit G2): no intermediate strip enters copy (1, 1) through side 3
        assert all(not ((o, e) == (1, 1) and s == 3) for (o, e, s, _) in seg[1:]), 'an intermediate strip re-enters J'
    lo0, hi0, hi1 = pinned('lo0'), pinned('hi0'), pinned('hi1')
    lam0, lam1 = hi0 - lo0, hi1 - hi0
    _, S0 = walk(J0)
    Ss1, S1 = walk(J1)
    assert ((lo0 - S0) - (lo0 + lam1)).iszero(), 'T^24 J0 is not [lo0 + λ1, hi1)'
    assert ((hi0 - S1) - lo0).iszero(), 'T^18 J1 is not [lo0, lo0 + λ1)'
    assert any(((hi0 - Ss1[16]) - v).iszero() for v in side_range(J1[16][:2], 1)), 'J1 lower end not on its separatrix'
    assert (lam0 - 4 * sinL(4) * sinL(1) * (c2 - 1)).iszero(), 'λ0 closed form'
    assert (lam1 - 16 * sinL(4) * sinL(1) * sinL(1) * sinL(1)).iszero(), 'λ1 closed form'
    out['identities'] = 'closure, itinerary, both returns, λ0 = 4 sin4x sinx (2cos2x−1), λ1 = 16 sin4x sin³x: all FORMAL'
    ineqs = []
    for tag, seg, lo, lam in (('J0', J0, lo0, lam0), ('J1', J1, hi0, lam1)):
        Ss, _ = walk(seg)
        for k, (o, e, s, s2) in enumerate(seg):
            g = (o, e)
            a, b = lo - Ss[k], lo - Ss[k] + lam
            for ss in (s, s2):
                p, q = side_range(g, ss)
                if p.evx(math.pi / 17).real > q.evx(math.pi / 17).real:
                    p, q = q, p
                ineqs += [(f'{tag}.{k}.s{ss}.lo', a - p), (f'{tag}.{k}.s{ss}.hi', q - b)]
            for ss, want in ((s, 1), (s2, -1)):
                i, kk = SIDES[ss]
                ev = g_apply(g, V[kk]) - g_apply(g, V[i])
                ineqs.append((f'{tag}.{k}.s{ss}.{"entry" if want == 1 else "exit"}', (ev + ev.conj()) * ((1 if e == 0 else -1) * want)))
    ineqs += [('lam0', lam0), ('lam1', lam1)]
    tight = [n for n, f in ineqs if f.iszero()]
    strict = [(n, f) for n, f in ineqs if not f.iszero()]
    cert = {n: certify_nonneg(f) for n, f in strict}
    out['inequalities'] = {'total': len(ineqs), 'formally_tight': len(tight), 'certified_on_(0,pi/15]': len(cert),
                           'leading_orders': sorted({v[0] for v in cert.values()}),
                           'min_lower_bound_f_over_x^n0': min(v[1] for v in cert.values())}
    ctrl = {N: sum(1 for _, f in strict if f.evx(math.pi / N).real <= 0) for N in (15, 13, 11, 9, 7)}
    assert ctrl[15] == 0 and ctrl[13] > 0, 'control did not fire'
    out['control_violations_at_N'] = ctrl
    out['rotation_number'] = {N: lam1.evx(math.pi / N).real / (lam0 + lam1).evx(math.pi / N).real for N in (15, 17, 101)}
    for N, r in out['rotation_number'].items():
        assert abs(r - 4 * math.sin(math.pi / N) ** 2) < 1e-12
    print(json.dumps(out, indent=1))
    if '--xcheck' in sys.argv:
        out['xcheck'] = xcheck(int(sys.argv[sys.argv.index('--xcheck') + 1]))
    json.dump(out, open('data/s935_pi_n_torus.json', 'w'), indent=1, default=str)


import flow_poly as FP          # noqa: E402  (xcheck only)
import flow_decomp as FD        # noqa: E402
from flow_decomp import Decomposer, Comp   # noqa: E402
from regular_pieces import asym       # noqa: E402


class PosDec(Decomposer):
    def tot(self, c):
        s = self.F.zero
        for l in c.top: s = self.F.add(s, self.lam[l])
        return s
    def run(self, top, bot, cap, L):
        F = self.F
        c0 = Comp(top, bot); self.OFF = {id(c0): F.zero}; self.SNAP = {}
        work, done, snaps = [c0], [], {}
        while work:
            c = work.pop()
            while True:
                T0 = self.tot(c)
                rest = self._reduce(c)
                if rest is not None:
                    self.OFF[id(rest)] = F.add(self.OFF[id(c)], self.tot(c)); work.append(rest); c.trk = c.saf0 = None
                if len(c.top) == 1:
                    c.kind = 'cylinder'; break
                if len(c.top) == 2 and id(c) not in snaps:
                    snaps[id(c)] = dict(off=self.OFF[id(c)], top=list(c.top), bot=list(c.bot), lam={l: self.lam[l] for l in c.top},
                                        h={l: self.h[l] for l in c.top})
                t0, b0 = c.top[0], c.bot[0]
                if t0 != b0 and F.cmp(self.lam[t0], self.lam[b0]) == 0:
                    T1 = self.tot(c)
                    it = c.bot.index(t0); c.bot[it] = b0; self.h[b0] += self.h[t0]
                    del c.bot[0]; del c.top[0]; c.trk = c.saf0 = None
                    self.OFF[id(c)] = F.add(self.OFF[id(c)], F.sub(T1, self.tot(c)))
                    continue
                if c.steps >= cap:
                    c.kind = 'undetermined'; break
                if c.steps % FD.BOSH_EVERY == 0 and c.steps > 0 and self.saf_nonzero(c):
                    y = self.boshernitzan(c)
                    if y is not None:
                        c.kind, c.cert = 'no_periodic', y; break
                c.steps += 1
                T1 = self.tot(c)
                r = self._zorich(c.top, c.bot) or self._zorich(c.bot, c.top)
                self.OFF[id(c)] = F.add(self.OFF[id(c)], F.sub(T1, self.tot(c)))
                if r: continue
            self.SNAP[id(c)] = snaps.get(id(c))
            done.append(c)
        return done

def run(N, P=None):
    P = P or asym(N)
    j = FP.perp_dirs_poly(P)[0]
    pieces, lengths, top, bot, meta = FP._build_iet_poly(P, j)
    F = FD._fast(P)
    Dc = PosDec(F, lengths)
    comps = Dc.run(top, bot, 50000, meta['L'])
    return P, F, Dc, comps, pieces, meta


def xcheck(N):
    """independent instrument: flow_poly at N; every non-periodic component's rotation number at its 2-label snapshot"""
    P, F, Dc, comps, pieces, meta = run(N)
    rows = []
    for c in comps:
        if c.kind == 'cylinder':
            continue
        s = Dc.SNAP[id(c)]
        lam = [F.fval(s['lam'][t]) for t in s['top']]
        rows.append({'h': [s['h'][t] for t in s['top']], 'rho': min(lam) / sum(lam)})
    target = 4 * math.sin(math.pi / N) ** 2
    res = {'N': N, 'n_nonperiodic': len(rows), 'expected_2floor': 2 * ((N - 7) // 8),
           'all_rho_eq_4sin2': all(abs(r['rho'] - target) < 1e-9 for r in rows), 'rows': rows}
    print('xcheck', {k: v for k, v in res.items() if k != 'rows'})
    return res


if __name__ == '__main__':
    main()
