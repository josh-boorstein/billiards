"""[Deposit] Section 5 and Appendix B: the exact certificate of a finite partition. Paper: "Finite perpendicular-beam partitions at irrational angles". Research name: probes/s973_exact_partition.py, copied verbatim except imports (deposit names), paths (data/ and figs/ under irrational_beams/) and this paragraph. Run through run_all.py.

s973 — EXACT certificate that a right triangle's L1 perpendicular-beam partition is FINITE, at an algebraic α given by the
minimal polynomial of σ₁ = 4cos²α (general form of probes/s972_unit_family.py, which certifies one fixed two-word family).

OBJECT (foundations / [BEAM-INF-UNIT]): R = (0,0), A = (0,1), O = (−cot α, 0); launch from (0, s), s ∈ (0,1), direction (−1, 0);
a branch = maximal s-interval with one side word up to the first head-on return to L1.
METHOD. Exact arithmetic in K = ℚ(σ)(t), t² = 1 − (σ/2 − 1)² (t = sin 2α); zero tests exact, signs from a 120-digit evaluation
(s976, queue s974-f: now RIGOROUS — σ in a rational bracket certified by an exact sign change, every sign from an mpmath.iv
enclosure that must exclude 0, the smallest margin reported as sign_margin; binding ends chosen by `pick`, exact on every overlap,
with ALL exactly-tied labels recorded (lo_ties/hi_ties); `locate` finds the non-automatic coincidences by a rational perturbation of σ).
Before s976 (asserted |v| > 1e-80): sign margins were not bounded a priori: at the two Lemma-12 angles the smallest |v| is 0.011 against
coefficient sums ≤ 1.3e3 (margin ≥ 1e111, audits/s974_audit_tro_l12.md row 11); for [UNIT-FINITE] n ≤ 8, min |v| = 0.024 (s974 audit).
Near-ties between candidate exits are frequent (52/62 at the Lemma-12 angles, all EXACTLY equal, s974); the recorded lo/hi vertex is ONE
of the tied constraints, so it does not say where the non-automatic coincidence sits (row 19 of that audit). Each launch is traced SYMBOLICALLY in s: positions are affine forms a + b·s over K (the direction does not
depend on s), so each side parameter u_k(s) is affine and the word's realisation interval ∩_k {0 ≤ u_k ≤ 1} is computed exactly
(convexity: a forward ray from inside meets the boundary once, so u ∈ [0,1] on a side makes it THE exit). The side at each step is
chosen at one rational interior launch s₀ (a rational near the float branch's midpoint). The partition is CERTIFIED finite iff the
exact intervals, taken left to right, tile [0, 1] with exactly-equal shared ends.
Seeds: the float tracer's branches (beam_families.trace_nb at ε = 1e-9); a capped/missing piece ⇒ no certificate attempted.
Not PRE-REGISTERED (an instrument). PRIOR ART: probes/s972_unit_family.py (fixed family, two launches per half); none general.
DATA: data/s972_unit_family.json (the family, reproduced as a control by `control`).
Usage: python probes/exact_partition.py control | poly c_d,...,c_0 lo hi   (σ root bracket)
"""
import _setup  # noqa: F401  (deposit: this directory on the path, cwd = irrational_beams/)
import math
import os
import sys
from fractions import Fraction

import sympy as sp
from mpmath import mp, mpf, iv

sys.path[:0] = [os.path.dirname(os.path.abspath(__file__))]
mp.dps = 120
iv.dps = 120
S = sp.Symbol('S')


class Field:
    def __init__(self, m, sig):
        self.m = sp.Poly(m, S, domain='QQ')
        self.sig = mpf(sig)
        self.D = sp.Poly(1 - (S / 2 - 1) ** 2, S, domain='QQ').rem(self.m)
        self.tnum = mp.sqrt(1 - (self.sig / 2 - 1) ** 2)
        # s976 (s974-f): RIGOROUS evaluation — σ in a rational bracket certified by an exact sign change of the min-poly
        if self.m.degree() == 1:
            c1, c0 = self.m.all_coeffs()
            r = -sp.Rational(c0) / sp.Rational(c1)
            self.sigI = iv.mpf(int(r.p)) / int(r.q)
        else:
            eps = mpf(10) ** -110
            lo, hi = (sp.Rational(mp.nstr(self.sig + d, 125, strip_zeros=False)) for d in (-eps, eps))
            assert self.m.eval(lo) * self.m.eval(hi) < 0, 'sigma bracket has no sign change'
            self.sigI = iv.mpf([iv.mpf(int(lo.p)) / int(lo.q), iv.mpf(int(hi.p)) / int(hi.q)])
        self.tI = iv.sqrt(1 - (self.sigI / 2 - 1) ** 2)
        self.minabs = None                      # smallest |v| lower bound over all sign decisions (the reported margin)
        self.nsign = 0

    def evI(self, a):
        if a.is_zero:
            return iv.mpf(0)
        return sum(iv.mpf(int(c.p)) / int(c.q) * self.sigI ** k for (k,), c in a.terms())

    def red(self, p): return p.rem(self.m)
    def q(self, a): return sp.Poly(sp.nsimplify(a) if not isinstance(a, Fraction) else sp.Rational(a.numerator, a.denominator), S, domain='QQ')
    def inv1(self, a): return sp.Poly(sp.invert(a.as_expr(), self.m.as_expr(), S), S, domain='QQ')

    def ev1(self, a):
        if a.is_zero:
            return mpf(0)
        return sum(mpf(int(c.p)) / mpf(int(c.q)) * self.sig ** k for (k,), c in a.terms())


class E:
    F = None

    def __init__(self, a, b=None):
        F = E.F
        self.a = F.red(a if isinstance(a, sp.Poly) else F.q(a))
        self.b = F.red(b if isinstance(b, sp.Poly) else F.q(0 if b is None else b))

    def __add__(s, o): o = o if isinstance(o, E) else E(o); return E(s.a + o.a, s.b + o.b)
    def __sub__(s, o): o = o if isinstance(o, E) else E(o); return E(s.a - o.a, s.b - o.b)
    def __neg__(s): return E(-s.a, -s.b)

    def __mul__(s, o):
        o = o if isinstance(o, E) else E(o); F = E.F
        return E(s.a * o.a + F.red(s.b * o.b * F.D), s.a * o.b + s.b * o.a)

    def inv(s):
        F = E.F
        n = F.red(s.a * s.a - F.red(s.b * s.b * F.D))
        ni = F.inv1(n)
        return E(F.red(s.a * ni), F.red(-s.b * ni))

    def __truediv__(s, o): o = o if isinstance(o, E) else E(o); return s * o.inv()
    def iszero(s): return s.a.is_zero and s.b.is_zero
    def num(s): return E.F.ev1(s.a) + E.F.ev1(s.b) * E.F.tnum

    def numI(s): return E.F.evI(s.a) + E.F.evI(s.b) * E.F.tI

    def sign(s):
        """Exact zero test, else the sign of an interval enclosure that must exclude 0 (s976 hardening, s974-f)."""
        if s.iszero():
            return 0
        v = s.numI()
        assert 0 not in v, 'sign undecidable'
        F = E.F
        low = min(abs(float(v.a)), abs(float(v.b)))
        F.minabs = low if F.minabs is None else min(F.minabs, low)
        F.nsign += 1
        return 1 if float(v.mid) > 0 else -1


class Aff:
    """c0 + c1·s over K."""

    def __init__(self, c0, c1=None):
        self.c0 = c0 if isinstance(c0, E) else E(c0)
        self.c1 = c1 if isinstance(c1, E) else E(0 if c1 is None else c1)

    def __add__(s, o): o = o if isinstance(o, Aff) else Aff(o); return Aff(s.c0 + o.c0, s.c1 + o.c1)
    def __sub__(s, o): o = o if isinstance(o, Aff) else Aff(o); return Aff(s.c0 - o.c0, s.c1 - o.c1)
    def scale(s, k): return Aff(s.c0 * k, s.c1 * k)
    def at(s, x): return s.c0 + s.c1 * x


def setup(m, sig):
    E.F = Field(m, sig)
    ONE, ZERO = E(1), E(0)
    cos2 = E(sp.Poly(S / 2 - 1, S, domain='QQ'))
    t = E(0, 1)
    cot = (ONE + cos2) / t
    R, A, O = (ZERO, ZERO), (ZERO, ONE), (-cot, ZERO)
    return dict(L1=(R, A), L2=(R, O), H=(O, A))


def cr(a, b): return a[0] * b[1] - a[1] * b[0]


def pick(cands, want_max):
    """Rigorous max (min) of exact candidates: interval enclosures, exact resolution of every overlap; returns (best, all labels
    EXACTLY equal to it). s976 hardening (s974-f; audits/s974_audit_tro_l12.md rows 12, 19)."""
    I = [c[0].numI() for c in cands]
    key = (lambda v: float(v.a)) if want_max else (lambda v: -float(v.b))
    b = max(range(len(cands)), key=lambda i: key(I[i]))
    for i in range(len(cands)):                    # anything overlapping the provisional best is compared exactly
        if i != b and not (I[i].b < I[b].a or I[i].a > I[b].b):
            sg = (cands[i][0] - cands[b][0]).sign()
            if sg != 0 and (sg > 0) == want_max:
                b = i
    ties = [c[1:] for i, c in enumerate(cands)
            if not (I[i].b < I[b].a or I[i].a > I[b].b) and (c[0] - cands[b][0]).iszero()]
    pick.near = getattr(pick, 'near', 0) + sum(1 for i in range(len(cands)) if i != b and not (I[i].b < I[b].a or I[i].a > I[b].b))
    return cands[b], ties


def trace_branch(SIDES, s0, maxk=10 ** 4):
    """Trace launch s0 (E) symbolically in s; return (word, lo, hi, end_heights) with [lo, hi] the exact realisation interval."""
    ONE, ZERO = E(1), E(0)
    px, py = Aff(ZERO), Aff(ZERO, ONE)
    d = (-ONE, ZERO)
    word, last = [], 'L1'
    lows, highs = [(ZERO, 0, 'launch-R')], [(ONE, 0, 'launch-A')]                      # constraints on s: (root, step, vertex)
    VX = dict(L1=('R', 'A'), L2=('R', 'O'), H=('O', 'A'))
    for k in range(maxk):
        best = None
        for name, (Pp, Qp) in SIDES.items():
            if name == last:
                continue
            e = (Qp[0] - Pp[0], Qp[1] - Pp[1])
            den = cr(d, e)
            if den.iszero():
                continue
            wx, wy = Aff(Pp[0]) - px, Aff(Pp[1]) - py   # w = P − p
            # u = cr(w, d)/den, tt = cr(w, e)/den  (affine in s)
            u = Aff(wx.c0 * d[1] - wy.c0 * d[0], wx.c1 * d[1] - wy.c1 * d[0]).scale(den.inv())
            tt = Aff(wx.c0 * e[1] - wy.c0 * e[0], wx.c1 * e[1] - wy.c1 * e[0]).scale(den.inv())
            u0, t0 = u.at(s0), tt.at(s0)
            if t0.sign() > 0 and u0.sign() > 0 and (u0 - ONE).sign() < 0:
                best = (name, tt, u, e)
        assert best is not None, 'no exit side (launch on a vertex graze?)'
        name, tt, u, e = best
        word.append(name)
        for target in (ZERO, ONE):                    # constraint u(s) ≥ 0 and u(s) ≤ 1
            if u.c1.iszero():
                continue
            root = (target - u.c0) / u.c1
            # u − target has sign of c1·(s − root); we need u ≥ 0 (target 0) / u ≤ 1 (target 1)
            want_pos = target is ZERO
            vx = VX[name][0 if target is ZERO else 1]
            if (u.c1.sign() > 0) == want_pos:
                lows.append((root, k + 1, vx))
            else:
                highs.append((root, k + 1, vx))
        px, py = px + tt.scale(d[0]), py + tt.scale(d[1])
        ee = e[0] * e[0] + e[1] * e[1]; dot = d[0] * e[0] + d[1] * e[1]
        d = (dot * e[0] * 2 / ee - d[0], dot * e[1] * 2 / ee - d[1])
        last = name
        if name == 'L1' and d[1].iszero() and (d[0] + ONE).iszero():
            lo, lo_ties = pick(lows, True); hi, hi_ties = pick(highs, False)
            trace_branch.ends = (lo[1:], hi[1:])
            trace_branch.ties = (lo_ties, hi_ties)
            return word, lo[0], hi[0], py
    raise AssertionError('no head-on return within %d' % maxk)


def certify(m, sig, seeds, verbose=True):
    """seeds: float launch points (one per float branch, sorted). Returns dict with exact tiling verdict."""
    SIDES = setup(m, sig)
    cur = E(0)
    rows, ok = [], True
    for s in sorted(seeds):
        s0 = E(Fraction(s).limit_denominator(10 ** 12))
        w, lo, hi, py = trace_branch(SIDES, s0)
        shared = (lo - cur).iszero()
        ok &= shared
        rows.append(dict(s0=s, word=' '.join(w), nref=len(w), lo=float(lo.num()), hi=float(hi.num()), shares_left_end=shared,
                         lo_vertex=trace_branch.ends[0], hi_vertex=trace_branch.ends[1],
                         lo_ties=trace_branch.ties[0], hi_ties=trace_branch.ties[1],
                         lo_exact=(str(lo.a.as_expr()), str(lo.b.as_expr()))))
        if verbose:
            print(f"  [{float(lo.num()):.12f}, {float(hi.num()):.12f}] nref={len(w)} shared={shared}  {' '.join(w)[:120]}", flush=True)
        cur = hi
    ok &= (cur - E(1)).iszero()
    return dict(CERTIFIED=bool(ok), branches=rows, sign_margin=E.F.minabs, n_sign=E.F.nsign)


def graze_events(labels):
    """Distinct physical graze events among tied (step, vertex) labels: (k, v) and (k+1, v) are one hit (end of a step = start of the
    next)."""
    ev = []
    for k, v in sorted(set(tuple(x) for x in labels)):
        if ev and ev[-1][1] == v and k - ev[-1][2] <= 1:
            ev[-1] = (ev[-1][0], v, k)
        else:
            ev.append((k, v, k))
    return [(k0, v) for k0, v, _ in ev]


def coincidence_ends(c):
    """EXACT classification of the shared ends of a certified partition (s976, s974-f): an end is a COINCIDENCE iff either side is
    bound by ≥ 2 distinct graze events at exactly the same s (e.g. an O-graze at step 44 and an A-graze at step 57 — a saddle
    connection); otherwise automatic (one graze bounds both sides). Validated on the two Lemma-12 points: flags audits/
    s974_audit_tro_l12.md row 19's ends (0.546171…, 0.157884…) and two more per point. A rational perturbation of σ opens a coincidence
    only on ONE side (`locate`), which is why row 19 saw one."""
    bs = c['branches']; out = []
    for i in range(1, len(bs)):
        evL, evR = graze_events(bs[i - 1]['hi_ties']), graze_events(bs[i]['lo_ties'])
        if max(len(evL), len(evR)) > 1:
            out.append(dict(end=i, s=bs[i]['lo'], left=evL, right=evR))
    return out


def locate(m, sig, seeds, words=None, deltas=(mpf('1e-4'), mpf('-7e-5'))):
    """One-sided heuristic (superseded by coincidence_ends): re-run the same seeds in ℚ(t) at rational σ' = σ + δ and list the ends that
    open. A coincidence opens only for the sign of δ that moves its later graze past the earlier one. s974-f."""
    out = []
    for d in deltas:
        q = sp.Rational(mp.nstr(sig + d, 12))
        c = certify(sp.Poly(q.q * S - q.p, S, domain='QQ').as_expr(), mpf(q.p) / q.q, seeds, verbose=False)
        bs = c['branches']
        same = words is None or [b['nref'] for b in bs] == list(words)
        out.append(dict(delta=float(d), words_same=bool(same), top_hi_vertex=bs[-1]['hi_vertex'],
                        top_at_1=bool(abs(bs[-1]['hi'] - 1) < 1e-30),
                        open_ends=[(i, bs[i]['lo'], bs[i]['lo_vertex']) for i in range(1, len(bs)) if not bs[i]['shares_left_end']]))
    return out


def float_seeds(x, eps=1e-9, budget=10 ** 8):
    import launch_side as L, beam_families as B
    r = B.trace_nb(L.table('right', math.pi * x / 2, 0), eps_min=eps, budget=budget)
    if r['capped']:
        return None, r
    return [0.5 * (b[0] + b[1]) for b in sorted(r['branches'])], r


def run_poly(coeffs, lo, hi, eps=1e-9):
    """coeffs: integer min-poly of σ (highest first); [lo, hi] isolates the real root."""
    m = sum(sp.Integer(c) * S ** (len(coeffs) - 1 - i) for i, c in enumerate(coeffs))
    try:
        rts = [r for r in sp.Poly(m, S).nroots(n=130, maxsteps=1000) if r.is_real and lo < r < hi]
    except Exception:                                  # s979: nroots non-convergence at high degree → real-root isolation instead
        rts = [sp.Float((a + b) / 2, 130) for (a, b), _ in sp.Poly(m, S).intervals(eps=sp.Rational(1, 10 ** 60)) if lo < a and b < hi]
    assert len(rts) == 1, rts
    sig = mpf(str(rts[0]))
    al = mp.acos(mp.sqrt(sig) / 2)
    x = float(2 * al / mp.pi)
    seeds, r = float_seeds(x, eps)
    if seeds is None:
        return dict(x=x, CERTIFIED=False, reason='float capped', n_float=len(r['branches']), n_capped=len(r['capped']))
    out = certify(m, sig, seeds)
    out.update(x=x, sigma=float(sig), minpoly=coeffs, n=len(seeds))
    return out


if __name__ == '__main__':
    if sys.argv[1] == 'control':
        for coeffs, lo, hi in (([1, -2, -1], 2, 3), ([1, -4, 3, -1], 3, 3.3)):
            r = run_poly(coeffs, lo, hi)
            print(coeffs, 'n =', r.get('n'), 'CERTIFIED', r['CERTIFIED'])
    elif sys.argv[1] == 'poly':
        coeffs = [int(c) for c in sys.argv[2].split(',')]
        r = run_poly(coeffs, float(sys.argv[3]), float(sys.argv[4]))
        print({k: v for k, v in r.items() if k != 'branches'})


# ---- census mode: build each finite row's coincidence equation as a polynomial in σ, factor, certify ----
def hpoly(name):
    """vertex height in the first-leg unfolding as a polynomial in σ = 4cos²α:
    A_{2k+1}: D_k = sin((2k+1)α)/sin α = U_{2k}(c);  R_k: cot α sin 2kα = c·U_{2k−1}(c)  (c = cos α, both even in c)."""
    c = sp.Symbol('c')
    k = int(name[1:])
    if name[0] == 'A':
        e = sp.chebyshevu(k - 1, c)
    else:
        e = sp.expand(c * sp.chebyshevu(2 * k - 1, c))
    return sp.expand(e.subs(c, sp.sqrt(S) / 2))


def tag_poly(tag):
    t, rest = tag.split(':')
    if t == 'T1':
        return sp.expand(hpoly(rest) - sp.Rational(1, 2))
    v, w = rest.split('=')
    if t == 'T2':
        return sp.expand(hpoly(v) + hpoly(w) - 1)
    if t == 'T3':
        return sp.expand(hpoly(v) - hpoly(w))
    if t == 'T4':                                   # s976 census extension
        return sp.expand(hpoly(v) + hpoly(w))
    if t in ('T5+1', 'T5-1'):
        return sp.expand(hpoly(v) - hpoly(w) - (1 if t == 'T5+1' else -1))
    raise ValueError(tag)


def minpoly_at(tags, sig):
    for tag in tags:
        P = tag_poly(tag)
        for fac, _ in sp.factor_list(P, S)[1]:
            fp = sp.Poly(fac, S)
            if fp.degree() > 0 and abs(fp.eval(sp.Float(sig, 30))) < 1e-8 * max(1, max(abs(float(x)) for x in fp.all_coeffs())):
                den = sp.ilcm(*[sp.fraction(x)[1] for x in fp.all_coeffs()])
                co = [int(x * den) for x in fp.all_coeffs()]
                g = sp.igcd(*co); co = [x // g for x in co]
                if co[0] < 0:
                    co = [-x for x in co]
                return co, tag
    return None, None


def arith(co, sig):
    m = sum(sp.Integer(c) * S ** (len(co) - 1 - i) for i, c in enumerate(co))
    rts = sp.Poly(m, S).nroots(n=40)
    unit = abs(co[-1]) == 1                     # constant term of the primitive min-poly: unit of ℤ[σ] needs ±1 there
    alg_int = abs(co[0]) == 1
    outside = any(abs(complex(r).imag) > 1e-20 or not (0 <= float(sp.re(r)) <= 4) for r in rts)
    # α/π rational ⇒ σ − 2 = 2cos 2α is an algebraic integer with all conjugates in [−2, 2] (Kronecker): fails ⇒ irrational
    return dict(unit=unit, alg_int=alg_int, irrational=(not alg_int) or outside, degree=len(co) - 1)


def census_certify(path):
    import json
    D = json.load(open(path))
    out = []
    for r in D['rows']:
        if not r['finite']:
            continue
        co, tag = minpoly_at(r['tags'], r['sigma1'])
        row = dict(x=r['x'], alpha=r['alpha'], sigma1=r['sigma1'], tags=r['tags'], n_float=r['n9'], minpoly=co, from_tag=tag)
        if co is None:
            row['CERTIFIED'] = False; row['reason'] = 'no factor'
        else:
            row.update(arith(co, r['sigma1']))
            m = sum(sp.Integer(c) * S ** (len(co) - 1 - i) for i, c in enumerate(co))
            rts = [z for z in sp.Poly(m, S).nroots(n=130) if z.is_real]
            sig = min(rts, key=lambda z: abs(float(z) - r['sigma1']))
            print(f"x={r['x']:.10f} α={r['alpha']:.6f} ({math.degrees(r['alpha']):.2f}°) minpoly={co} {tag}", flush=True)
            c = certify(m, mpf(str(sig)), [0.5 * (b[0] + b[1]) for b in r['branches']], verbose=False)
            row['CERTIFIED'] = c['CERTIFIED']
            row['branches_exact'] = [(b['lo'], b['hi'], b['nref'], b['lo_vertex'], b['hi_vertex'], b['word']) for b in c['branches']]
            row['ties'] = [(b['lo_ties'], b['hi_ties']) for b in c['branches']]          # s976: every exactly-tied label
            row['coincidences'] = coincidence_ends(c)
            row['sign_margin'] = c['sign_margin']; row['n_sign'] = c['n_sign']
            print(f"    n={len(c['branches'])} CERTIFIED={c['CERTIFIED']} unit={row['unit']} irrational={row['irrational']} deg={row['degree']} "
                  f"margin={c['sign_margin']:.3g} coincidence_ends={[round(z['s'], 6) for z in row['coincidences']]}",
                  flush=True)
        out.append(row)
    json.dump(out, open(path.replace('.json', '_certified%s.json' % ('_hardened' if os.environ.get('HARDENED_TAG') else '')), 'w'),
              indent=1, default=str)
    return out


if __name__ == '__main__' and sys.argv[1] == 'census':
    census_certify(sys.argv[2])
