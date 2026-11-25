"""[Deposit] Section 5: exact height polynomials of a word. Paper: "Finite perpendicular-beam partitions at irrational angles". Research name: probes/s977_fan_identity.py, copied verbatim except imports (deposit names), paths (data/ and figs/ under irrational_beams/) and this paragraph. Run through run_all.py.

s977 — THE FAN IDENTITY OF EACH FINITE POINT, AS POLYNOMIALS: which R-images are the fan centres, and which equations does
the geometry impose? (queue s976-d; user s977: *realisability of the conditions imposed by geometry must play a role in the units
that are allowed*.)

OBJECT. Theorem I′'s proof (general_triangles.md §s961 (2)–(4), [BEAM-INF-UNIT]): a finite L1 partition is tiled by m fans
[e_{k−1}, e_k], each two mirror beams split at its centre c_k, the landing height of a pole's prong; 2c_k = 2y(g_k(R)) ∈ σℤ[σ]
(Lemma (3)), and e_k = 2c_k − e_{k−1}, e_0 = 0, e_m = 1, so 1 = Σ_k (−1)^{m−k} 2c_k (THE FAN IDENTITY). Theorem I′ reads only the
constant term of that identity (⇒ σ a unit). Here we read the WHOLE identity: the polynomials f_k with 2c_k = σ f_k(σ).
INSTRUMENT. The certified exact partition (probes/exact_partition.py) gives every end's exactly-tied graze labels (step, vertex).
For each graze label we take its branch's word up to that step and compute the graze height as a function of σ by tracing the FIXED word
at rational σ′ (exact, ℚ(√D), D = σ′(4 − σ′)/4) and interpolating; verified at extra nodes and at σ₁ (to 1e-25). A centre end must
carry an R graze; a boundary end an O or A graze.
PRE-REGISTERED.
  (F1) ends alternate boundary / centre, n = 2m (audit P5 of s963: proved — a CHECK of the instrument, cannot fail if it works).
  (F2) every centre height is σ·f(σ), f ∈ ℤ[σ] (Lemma (3): proved — instrument check).
  (F3) every fan boundary e_k that is NOT an exact multi-graze end is a SINGLE O/A graze whose height polynomial equals
       2c_k − e_{k−1} IDENTICALLY in σ (mirror symmetry is structural: no equation spent there).
  (F4) THE QUESTION: how many INDEPENDENT polynomial conditions does each point carry — the top end (e_m = 1), plus every coincidence
       end (≥ 2 graze events at one s) whose difference polynomial is NOT identically 0? Each is divisible by the minimal polynomial.
       Could come out otherwise: one condition everywhere (the top), or several (an over-determined system that a unit happens to solve).
  (F5) at every OTHER real conjugate σ′ ∈ (0, 4), the same polynomial identity holds (algebra), so realisability alone decides: report
       which geometric inequality of the fan data (0 < e_{k−1} < c_k < e_k, ordering) fails there.
RESULT (s977, the 12 certified unit-census points → data/fan_identity.json, …_conjugates.json):
  (F1) HOLDS: n = 2m everywhere; all 29 odd ends are pure R grazes (fan centres), all 17 even ends carry an A or O graze.
  (F2) HOLDS: 67 R heights are σ·f(σ), f ∈ ½ℤ[σ] (denominator 2: 2y ∈ σℤ[σ] as Lemma (3)). NEW, same instrument: every A height
       (34) has INTEGER coefficients and constant term ±1, every O height (24) integer coefficients and constant term 0 — so an A–O
       height coincidence E = y(gA) − y(g′O) has E(0) = ±1 and ALL its roots are units (MEASURED; no proof written).
  (F3) HOLDS (by the count: one equation per point, below).
  (F4) ONE CONDITION PER POINT, 12/12: every non-identical difference at an end, the top or the bottom gives the SAME polynomial E (up to
       sign) — an A–O horizontal saddle connection; E = (minimal polynomial) × (companion: σ − 1, σ² − 4σ + 1, σ³ − 5σ² + 4σ + 1, or 1).
       The several coincidence ends of a point (three at L12_a) are one equation, not several.
  (F5) at every OTHER real root of E in (0, 4) — 8 conjugates, 13 companion roots — the SAME words are NOT realised (leaf height outside
       [0, 1] or an intermediate side parameter out of (0, 1), margins −1 … −10²; degenerate at the rational companions); at the 12 points
       the margins are 0 to 10⁻⁴⁸. The two conjugate pairs finite together are finite by DIFFERENT words (ℚ(ζ₇)⁺: the same E up to sign).
PRIOR ART: §s961 (2)–(4) (the fan identity, Lemma (3)); §s961 *Residual* (outer chain O ⇝ A, angle π); §s976 (4) (coincidence ends,
three per Lemma-12 point); [NCYL-748] (conjugates not finite together: different branch counts at the finite pairs); rulings --grep
'fan identity' → NCYL-617/652/666/672/673 (the RATIONAL necklace fan identity — not at irrational α); 'pole height', 'outer chain' →
none beyond [LIT-060]; 'realisab' → nothing on beam partitions.
DATA: data/s976_unit_census_D3_H6.json (float branches, seeds), data/s976_unit_census_D3_H6_certified.json; data_index → none.
Usage: python probes/fan_identity.py [x-prefix …] → data/fan_identity.json
"""
import _setup  # noqa: F401  (deposit: this directory on the path, cwd = irrational_beams/)
import json
import math
import os
import sys
from fractions import Fraction as Fr

import sympy as sp
from mpmath import mp, mpf

sys.path[:0] = [os.path.dirname(os.path.abspath(__file__))]
import exact_partition as X

S = sp.Symbol('S')
OUT = 'data/fan_identity.json'


# ---------- exact ℚ(√D) arithmetic at rational σ′ ----------
class Q2:
    D = None
    __slots__ = ('a', 'b')

    def __init__(self, a, b=0): self.a, self.b = Fr(a), Fr(b)
    def __add__(s, o): o = o if isinstance(o, Q2) else Q2(o); return Q2(s.a + o.a, s.b + o.b)
    def __sub__(s, o): o = o if isinstance(o, Q2) else Q2(o); return Q2(s.a - o.a, s.b - o.b)
    def __neg__(s): return Q2(-s.a, -s.b)
    def __mul__(s, o): o = o if isinstance(o, Q2) else Q2(o); return Q2(s.a * o.a + s.b * o.b * Q2.D, s.a * o.b + s.b * o.a)
    def inv(s):
        n = s.a * s.a - s.b * s.b * Q2.D
        return Q2(s.a / n, -s.b / n)
    def __truediv__(s, o): o = o if isinstance(o, Q2) else Q2(o); return s * o.inv()
    def iszero(s): return s.a == 0 and s.b == 0


def word_graze(word, step, vertex, sig):
    """height s at which the FIXED word's step-`step` exit side passes through `vertex` (u = 0 or 1), at rational σ′."""
    sig = Fr(sig)
    Q2.D = sig * (4 - sig) / 4
    t = Q2(0, 1)
    cot = (Q2(1) + Q2(sig / 2 - 1)) / t
    R, A, O = (Q2(0), Q2(0)), (Q2(0), Q2(1)), (-cot, Q2(0))
    SIDES = dict(L1=(R, A), L2=(R, O), H=(O, A))
    VX = dict(L1=('R', 'A'), L2=('R', 'O'), H=('O', 'A'))
    # position affine in s: p = (px0 + px1 s, py0 + py1 s)
    px0, px1, py0, py1 = Q2(0), Q2(0), Q2(0), Q2(1)
    d = (Q2(-1), Q2(0))
    for k, name in enumerate(word[:step], start=1):
        P, Qp = SIDES[name]
        e = (Qp[0] - P[0], Qp[1] - P[1])
        den = d[0] * e[1] - d[1] * e[0]
        wx0, wx1, wy0, wy1 = P[0] - px0, -px1, P[1] - py0, -py1
        u0 = (wx0 * d[1] - wy0 * d[0]) / den; u1 = (wx1 * d[1] - wy1 * d[0]) / den
        t0 = (wx0 * e[1] - wy0 * e[0]) / den; t1 = (wx1 * e[1] - wy1 * e[0]) / den
        if k == step:
            if vertex not in VX[name]:
                raise ValueError('vertex %s not on side %s' % (vertex, name))
            target = 0 if vertex == VX[name][0] else 1
            return (Q2(target) - u0) / u1
        px0, px1 = px0 + t0 * d[0], px1 + t1 * d[0]
        py0, py1 = py0 + t0 * d[1], py1 + t1 * d[1]
        ee = e[0] * e[0] + e[1] * e[1]; dot = d[0] * e[0] + d[1] * e[1]
        d = (dot * e[0] * 2 / ee - d[0], dot * e[1] * 2 / ee - d[1])
    raise ValueError('step beyond word')


def _is_sq(q):
    from math import isqrt
    return q >= 0 and isqrt(q.numerator) ** 2 == q.numerator and isqrt(q.denominator) ** 2 == q.denominator


def interp_poly(word, step, vertex, sig1, degmax):
    """graze height as (polynomial in σ for the rational part, polynomial for the √D part) — returns sympy Polys or None."""
    nodes = [Fr(1, 7) + Fr(k, 23 * (1 + degmax // 60)) for k in range(3 * degmax + 12)]           # rational σ′ in (0.14, ~3.7)
    nodes = [z for z in nodes if 0 < z < 4 and not _is_sq(z * (4 - z) / 4)][:degmax + 4]
    vals = [word_graze(word, step, vertex, z) for z in nodes]
    fit, chk = nodes[:degmax + 1], nodes[degmax + 1:]
    pa = sp.interpolate([(sp.Rational(z.numerator, z.denominator), sp.Rational(v.a.numerator, v.a.denominator)) for z, v in zip(fit, vals)], S)
    # √D part: b·√D with D = σ(4−σ)/4 — keep b as a rational function check; for R/O/A heights expected 0 or polynomial
    pb = sp.interpolate([(sp.Rational(z.numerator, z.denominator), sp.Rational(v.b.numerator, v.b.denominator)) for z, v in zip(fit, vals)], S)
    ok = all(sp.Rational(v.a.numerator, v.a.denominator) == pa.subs(S, sp.Rational(z.numerator, z.denominator)) and
             sp.Rational(v.b.numerator, v.b.denominator) == pb.subs(S, sp.Rational(z.numerator, z.denominator)) for z, v in zip(chk, vals[degmax + 1:]))
    return sp.Poly(pa, S), sp.Poly(pb, S), ok


def analyse(row):
    co = row['minpoly']
    m = sum(sp.Integer(c) * S ** (len(co) - 1 - i) for i, c in enumerate(co))
    rts = [z for z in sp.Poly(m, S).nroots(n=130) if z.is_real]
    sig = min(rts, key=lambda z: abs(float(z) - row['sigma1']))
    seeds = [0.5 * (b[0] + b[1]) for b in sorted(row['branches'])]
    c = X.certify(m, mpf(str(sig)), seeds, verbose=False)
    assert c['CERTIFIED']
    bs = c['branches']
    n = len(bs)
    ends = []
    for i in range(1, n):
        evL = X.graze_events(bs[i - 1]['hi_ties']); evR = X.graze_events(bs[i]['lo_ties'])
        ev = [('L', k, v) for k, v in evL] + [('R', k, v) for k, v in evR]
        ends.append(dict(i=i, s=bs[i]['lo'], ev=ev))
    out = dict(x=row['x'], minpoly=co, sigma=float(sig), n=n, words=[b['word'] for b in bs], ends=[])
    sigf = mp.mpf(str(sig))
    for e in ends:
        polys = []
        for side, k, v in e['ev']:
            if v.startswith('launch'):
                continue
            w = bs[e['i'] - 1 if side == 'L' else e['i']]['word'].split()
            pa, pb, ok = interp_poly(w, k, v, sig, degmax=k + 2)
            val = sum(mpf(int(cc.p)) / int(cc.q) * sigf ** j for (j,), cc in pa.terms()) if not pa.is_zero else mpf(0)
            valb = sum(mpf(int(cc.p)) / int(cc.q) * sigf ** j for (j,), cc in pb.terms()) if not pb.is_zero else mpf(0)
            polys.append(dict(side=side, step=k, vertex=v, a=[str(x) for x in pa.all_coeffs()], b=[str(x) for x in pb.all_coeffs()],
                              b_zero=pb.is_zero, interp_ok=ok, at_sigma=float(val), err=float(abs(val - e['s']) if pb.is_zero else -1)))
        kind = 'centre' if any(p['vertex'] == 'R' for p in polys) else 'boundary'
        # distinct height polynomials at this end
        distinct = []
        for p in polys:
            key = (tuple(p['a']), tuple(p['b']))
            if key not in [d[0] for d in distinct]:
                distinct.append((key, p))
        eqs = []
        for (ka, _kb), _p in distinct[1:]:
            diff = sp.Poly([sp.Rational(z) for z in ka], S) - sp.Poly([sp.Rational(z) for z in distinct[0][0][0]], S)
            r = diff.rem(sp.Poly(m, S))
            eqs.append(dict(identically_zero=diff.is_zero, divisible_by_minpoly=r.is_zero,
                            factors=str(sp.factor_list(diff.as_expr(), S)) if not diff.is_zero else '0'))
        out['ends'].append(dict(i=e['i'], s=e['s'], kind=kind, n_events=len(polys), n_distinct_polys=len(distinct), polys=polys, eqs=eqs))
    # top (s = 1, against the constant 1) and bottom (s = 0, against 0): each non-launch graze tied there is a polynomial condition
    for key, b, labels, const in (('top', bs[-1], bs[-1]['hi_ties'], 1), ('bottom', bs[0], bs[0]['lo_ties'], 0)):
        ev = X.graze_events(labels)
        tp = []
        for k, v in ev:
            if v.startswith('launch') or (k == 1):          # step 1 is the launch side's own end (R/A on L1, O via L2 at s = 0)
                continue
            pa, pb, ok = interp_poly(b['word'].split(), k, v, sig, degmax=k + 2)
            diff = pa - sp.Poly(const, S)
            tp.append(dict(step=k, vertex=v, a=[str(x) for x in pa.all_coeffs()], b_zero=pb.is_zero, interp_ok=ok,
                           identically=diff.is_zero, factors=str(sp.factor_list(diff.as_expr(), S)),
                           divisible=diff.rem(sp.Poly(m, S)).is_zero))
        out[key] = tp
        out[key + '_events'] = ev
    return out


if __name__ == '__main__' and not os.environ.get('CONJ'):
    D = json.load(open('data/s976_unit_census_D3_H6.json'))
    fin = sorted([r for r in D['rows'] if r['finite']], key=lambda r: r['x'])
    pre = sys.argv[1:]
    res = []
    for r in fin:
        if pre and not any(('%.10f' % r['x']).startswith(p) for p in pre):
            continue
        o = analyse(r)
        res.append(o)
        print(f"x={o['x']:.10f} minpoly={o['minpoly']} n={o['n']}", flush=True)
        for e in o['ends']:
            ps = ' | '.join(f"{p['side']}{p['step']}{p['vertex']} deg{len(p['a']) - 1}{'' if p['b_zero'] else '+√D'}{'' if p['interp_ok'] else ' !INTERP'}"
                            for p in e['polys'])
            eq = ' '.join('ID0' if q['identically_zero'] else ('EQ' if q['divisible_by_minpoly'] else 'NOTDIV') for q in e['eqs'])
            print(f"  end {e['i']} s={e['s']:.6f} {e['kind']:8s} distinct={e['n_distinct_polys']} [{eq}]  {ps}")
        for key in ('top', 'bottom'):
            print(f'  {key}:', [(t['step'], t['vertex'], 'ID' if t['identically'] else ('EQ' if t['divisible'] else 'NOTDIV'), t['factors'][:90]) for t in o[key]])
        for e in o['ends']:
            for q in e['eqs']:
                if not q['identically_zero']:
                    print(f"  end {e['i']} eq factors: {q['factors'][:110]}")
    json.dump(res, open(OUT, 'w'), indent=1, default=str)


# ---------- part 2 (F5): realisability of the SAME connection at every real root of its equation E ----------
def word_check(word, step, vertex, sig, dps=60):
    """float (mpmath) trace of the FIXED word at σ′: returns (s′, worst) — s′ the height where step `step` passes `vertex`, and worst =
    the most negative margin min(u_j, 1 − u_j, t_j) over the steps j < step at s′ (≥ 0 ⇔ the leaf really follows the word: realised)."""
    mp.dps = dps
    sig = mpf(sig)
    t = mp.sqrt(sig * (4 - sig) / 4)
    cot = (sig / 2) / t
    R, A, O = (mpf(0), mpf(0)), (mpf(0), mpf(1)), (-cot, mpf(0))
    SIDES = dict(L1=(R, A), L2=(R, O), H=(O, A))
    VX = dict(L1=('R', 'A'), L2=('R', 'O'), H=('O', 'A'))
    # first pass: affine coefficients to get s′
    px0, px1, py0, py1 = mpf(0), mpf(0), mpf(0), mpf(1)
    d = (mpf(-1), mpf(0)); steps = []
    for k, name in enumerate(word[:step], start=1):
        P, Qp = SIDES[name]
        e = (Qp[0] - P[0], Qp[1] - P[1]); den = d[0] * e[1] - d[1] * e[0]
        wx0, wx1, wy0, wy1 = P[0] - px0, -px1, P[1] - py0, -py1
        u0 = (wx0 * d[1] - wy0 * d[0]) / den; u1 = (wx1 * d[1] - wy1 * d[0]) / den
        t0 = (wx0 * e[1] - wy0 * e[0]) / den; t1 = (wx1 * e[1] - wy1 * e[0]) / den
        steps.append((u0, u1, t0, t1))
        px0, px1 = px0 + t0 * d[0], px1 + t1 * d[0]; py0, py1 = py0 + t0 * d[1], py1 + t1 * d[1]
        ee = e[0] * e[0] + e[1] * e[1]; dot = d[0] * e[0] + d[1] * e[1]
        d = (dot * e[0] * 2 / ee - d[0], dot * e[1] * 2 / ee - d[1])
    u0, u1, _, _ = steps[-1]
    target = 0 if vertex == VX[word[step - 1]][0] else 1
    s1 = (target - u0) / u1
    worst = mpf(1)
    for (u0, u1, t0, t1) in steps[:-1]:
        u = u0 + u1 * s1; tt = t0 + t1 * s1
        worst = min(worst, u, 1 - u, tt)
    return s1, worst


def conjugate_report(path=OUT):
    R = json.load(open(path))
    out = []
    for o in R:
        # the equation polynomials and their events (word, step, vertex) — one representative event pair per distinct equation
        evs = []
        bsw = o['words']
        for e in o['ends']:
            for p in e['polys']:
                w = bsw[e['i'] - 1 if p['side'] == 'L' else e['i']]
                evs.append((e['i'], w, p['step'], p['vertex'], tuple(p['a'])))
        for key, bi in (('top', len(bsw) - 1), ('bottom', 0)):
            for p in o[key]:
                evs.append((key, bsw[bi], p['step'], p['vertex'], tuple(p['a'])))
        # the equation(s)
        eqs = set()
        for e in o['ends']:
            for q in e['eqs']:
                if not q['identically_zero']:
                    eqs.add(q['factors'])
        for key in ('top', 'bottom'):
            for p in o[key]:
                if not p['identically']:
                    eqs.add(p['factors'])
        co = o['minpoly']
        mpoly = sp.Poly(sum(sp.Integer(c) * S ** (len(co) - 1 - i) for i, c in enumerate(co)), S)
        # full equation polynomial: rebuild from the first non-trivial pair
        Epolys = set()
        for e in o['ends']:
            cls = {}
            for p in e['polys']:
                cls.setdefault(tuple(p['a']), p)
            ks = list(cls)
            for k2 in ks[1:]:
                Epolys.add(sp.Poly([sp.Rational(z) for z in k2], S) - sp.Poly([sp.Rational(z) for z in ks[0]], S))
        for key, const in (('top', 1), ('bottom', 0)):
            for p in o[key]:
                Epolys.add(sp.Poly([sp.Rational(z) for z in p['a']], S) - sp.Poly(const, S))
        Epolys = [P for P in Epolys if not P.is_zero]
        rows = []
        for P in Epolys:
            for r in P.nroots(n=50):
                if not r.is_real or not (0 < r < 4):
                    continue
                rr = mpf(str(r))
                own = abs(float(mpoly.eval(sp.Float(str(r), 50)))) < 1e-20
                checks = []
                for (where, w, k, v, a) in evs:
                    if not (P - sp.Poly([sp.Rational(z) for z in a], S)).is_zero and \
                            not any((P + sp.Poly(c, S) - sp.Poly([sp.Rational(z) for z in a], S)).is_zero for c in (0,)):
                        pass
                    try:
                        s1, worst = word_check(w.split(), k, v, rr)
                    except ZeroDivisionError:
                        s1, worst = mpf('nan'), mpf(-9)          # degenerate (parallel side / vertex at the root)
                    checks.append(dict(where=where, step=k, vertex=v, s=float(s1), margin=float(worst)))
                x = float(2 * mp.acos(mp.sqrt(rr) / 2) / mp.pi)
                rows.append(dict(E=str(sp.factor_list(P.as_expr(), S)), root=float(rr), x=x, own_minpoly=own, is_point=abs(rr - o['sigma']) < 1e-12,
                                 min_margin=min(c['margin'] for c in checks), s_range=(min(c['s'] for c in checks), max(c['s'] for c in checks)),
                                 failing=[c for c in checks if c['margin'] < -1e-30 or not (-1e-30 < c['s'] < 1 + 1e-30)][:4]))
        out.append(dict(x=o['x'], minpoly=co, roots=rows))
        print(f"x={o['x']:.10f} minpoly={co}")
        seen = set()
        for r in rows:
            key = round(r['root'], 12)
            if key in seen:
                continue
            seen.add(key)
            tag = 'POINT' if r['is_point'] else ('conj ' if r['own_minpoly'] else 'compn')
            print(f"   {tag} σ′={r['root']:.6f} x′={r['x']:.6f}  min margin {r['min_margin']:+.3e}  s∈[{r['s_range'][0]:.4f},{r['s_range'][1]:.4f}]  "
                  f"fails: {[(f['where'], f['step'], f['vertex'], round(f['margin'], 4), round(f['s'], 4)) for f in r['failing']][:3]}")
    json.dump(out, open(OUT.replace('.json', '_conjugates.json'), 'w'), indent=1, default=str)


if __name__ == '__main__' and os.environ.get('CONJ'):
    conjugate_report()
