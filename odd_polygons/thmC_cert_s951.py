"""[Deposit] Theorem 5.4, step (iii), Appendix B.1: the certificate machinery (s951 region). Paper: "Perpendicular billiards in rational polygons of odd denominator". Research name: probes/s951_s1_allm_cert.py, copied verbatim except imports (deposit names), paths (data/ and figs/ under odd_polygons/) and this paragraph. Run through run_all.py; the usage lines below are the research repo's. Bracketed IDs and section marks (§sNNN) refer to the project's internal ledgers.

s951 keeper (queue `s946-b`): S₁(m) = (N−3m−1, 4m, N−m−1, 2)/N IS NON-CP FOR EVERY m ≥ 1 AND EVERY ODD N ≥ 14m + 1 — the
containment inequalities of the S₁(m) torus template (`probes/s951_series_first_look.py`: formal identities symbolic in m,
rotation number 4 sin²(mπ/N)) certified for all (m, N) at once. m = 1 is Π_N ([POLY-ASYM-ALL], §s935 (1)).

NOT PRE-REGISTERED: a certificate, not a hypothesis test. Inputs: `data/s951_s1_fams.pkl` (written by s951_series_first_look.py:
Laurent expressions in z = ζ_{4N} = e^{ix/2}, v = e^{imx}, u = z^{4dj} = e^{iφ}; x = π/N, θ = mx, φ = 2djx).
METHOD = `probes/s937_ladder_allm_cert.py`'s (exact binomial discharge ⊕ sinusoid lemma in φ ⊕ rigorous 2D interval certificate),
ported: region N ≥ 14m + 1 (θ(14 + ξ) ≤ π, ξ = x/θ = 1/m), coordinates θ = (1 − s)π/(14 + ξ), (s, ξ) ∈ [0, 1] × [0, 1/NMIN];
corner (s, ξ) = (0, 0) is θ = π/14 (m = ∞ on N = 14m + 1). Changes from s937, each needed here:
 (i) factoring tries sympy's default ℤ[i] factorisation first (fast), the gaussian=True split only if pairing fails;
 (ii) an EVEN power of a ≥ 3-term core is listed twice — certified of constant sign (so the square is > 0), contributing (±1)²;
 (iii) the corner-zero test is EXACT by reduction mod Φ_56 (one-sided — Φ_56 splits over ℚ(i): a True is a proof, a False only
      skips the gradient test, so soundness is unaffected; audit W4) (s937's sp.simplify cannot prove the heptagon identity
      cos π/7 − cos 2π/7 + cos 3π/7 = 1/2, which makes the core 2cos6θ − 2cos4θ + 2cos2θ − 1 vanish at the corner);
 (iv) NMIN = 4: at m = 3 the run (23)^{m−1}'s φ-derivative factor V − Z⁷ (sign of m/4 − 7/8) has the other sign, and at m = 2 a
      θ-divided core vanishes on ξ = 1/2 as θ → 0 — so m = 2, 3 are certified PER m: every family at that m, every repetition j,
      is a Laurent polynomial in z, certified ≥ 0 on x ∈ (0, X], X ≥ π/(14m + 1) rational, by s935's rigorous Taylor ⊕ interval
      Horner (`s935_pi_n_torus.certify_nonneg`).
CONTROLS (could it have come out otherwise?): (a) the 2D pipeline with the region widened to N ≥ 14m − 1 must fail on EXACTLY the
families numerically negative at N = 14m − 1 (m = 4, 5, 7), none negative at N = 14m + 1 — this exercises only the 2D-core path
(all 6 sit in J0's last word); so, per audit G6: (b) the per-m certificates with x widened to π/(14M − 1) must fail exactly on the
families negative at N = 14M − 1 (M = 2, 3); (c) every run one repetition longer (c = (γ, δ + 1), repetitions the walk never makes):
the unchanged 2D pipeline must flag EXACTLY the families that become negative at N = 14m + 1 — exercising the binomial j-closure
and the φ-end (sinusoid) atoms.
SINUSOID LEMMA (audit G5f): it needs only φ-interval length < 2π; here |φ| ≤ 2|d|·j·x with |d| = 2, j ≤ 3m − 2, so the length is
≤ 4(3m − 2)x < 12θ ≤ 12π/14 < π (s937's "≤ 4θ" bound does not transfer, and is not needed).
λ0 = 4 sinx sin2mx sin4mx (2cos2mx − 1), λ1 = 16 sinx sin²mx sin2mx sin4mx > 0 there (4mx < π, 2mx < π/3): J's two pieces are
non-empty, and its return is the rotation by 4 sin²(mπ/N) ∉ ℚ. ⟹ a point of J off the countable separatrix set has an infinite,
vertex-free, non-periodic orbit. Even m: the template's start copy has a half-integer offset — the formal copy is −(a real copy),
i.e. the unfolding of the reversed flow, which lies in the same perpendicular class (reflection in V0V1 reverses ω); orbits and
their periodicity are unchanged.
PRIOR ART: §s946 (3) ([POLY-PI-SERIES] MEASURED, onset 14m + 1 pre-registered 54/54); §s935 (1)/(4), §s937 (the method);
`rulings.py --grep` series, torus → nothing on S₁. DATA: `data/s951_s1_fams.pkl` (this session).
→ general_triangles.md §s951.   Usage: PYTHONPATH=… python probes/s951_s1_allm_cert.py [--no-control]   (≈ 6 min)"""
import _setup  # noqa: F401  (deposit: this directory on the path, cwd = odd_polygons/)
import os, sys, json, math, cmath, pickle, time, collections
import sympy as sp
import mpmath as mp
from mpmath import iv, mpf
from fractions import Fraction as Fr

ROOT = '.'   # deposit: cwd = odd_polygons/ (_setup)
OUTF = f'{ROOT}/data/s951_s1_allm_cert.json'
NMIN = 4
ASL = 14            # region N >= ASL*m + C
z, v, t, u = sp.symbols('z v t u')
Z, V = sp.symbols('Z V')
F0, F1, *_ = pickle.load(open(f'{ROOT}/data/s951_s1_fams.pkl', 'rb'))
FAMS = F0 + F1
mp.mp.dps = 50
iv.dps = 25
UNITS = {(1, 0): 0, (0, 1): 1, (-1, 0): 2, (0, -1): 3}          # ω = i^k

def gi(c):
    c = sp.nsimplify(c); re, im = sp.re(c), sp.im(c)
    return (Fr(int(sp.numer(re)), int(sp.denom(re))), Fr(int(sp.numer(im)), int(sp.denom(im))))

def expr_of(f, tag):
    return {'sign': lambda: f[3], 'a-p': lambda: f[3] - f[5], 'q-b': lambda: f[6] - f[4],
            'a-q': lambda: f[3] - f[6], 'p-b': lambda: f[5] - f[4]}[tag]()

def reps_of(f, n):
    c, kind = f[1], f[0].split('.')[3]
    r = range(1) if c == (0, 1) else range(c[0] * n + c[1])
    return [0] if kind == 'first' else ([j for j in r if j >= 1] if kind == 'rest' else list(r))

# ---------------------------------------------------------------- (0) orientation
def orient():
    rel = []
    for i, f in enumerate(FAMS):
        if f[4] is None:
            rel.append((i, 'sign')); continue
        fn = sp.lambdify((z, v, u), f[6] - f[5], 'mpmath'); sg = set()
        for n, N in ((4, 57), (5, 73), (7, 1001)):
            zz = mp.exp(1j * mp.pi / (2 * N))
            for j in reps_of(f, n):
                sg.add(mp.re(fn(zz, zz ** (2 * n), zz ** (4 * f[2] * j))) > 0)
        assert len(sg) == 1, f[0]
        rel += [(i, 'a-p'), (i, 'q-b')] if sg.pop() else [(i, 'a-q'), (i, 'p-b')]
    return rel

# ---------------------------------------------------------------- (1) factor, normalize, binomials
def _pair(KF, gens):
    K, fl0 = KF
    def recip(g):
        Pg = sp.Poly(g, *gens); sub = {sp.conjugate(x_): 1 / x_ for x_ in gens}
        mono = sp.Mul(*[x_ ** Pg.degree(x_) for x_ in gens])
        return sp.expand(sp.conjugate(g).subs(sub) * mono)
    fl, pend = [], list(fl0)
    while pend:
        g, m = pend.pop(0)
        if len(sp.Poly(g, *gens).terms()) == 1:
            fl.append((g, m)); continue
        gr = recip(g)
        if sp.cancel(gr / g).is_number:
            fl.append((g, m)); continue
        for k_, (h, mh) in enumerate(pend):
            if sp.cancel(gr / h).is_number:
                assert mh == m; pend.pop(k_); fl.append((sp.expand(g * h), m)); break
        else:
            raise ValueError(f'unpaired factor {g}')
    return K, fl

def factor_real(P, gens):
    PP = sp.Poly(P, *gens)
    if str(PP.domain) == 'ZZ':
        return sp.factor_list(P)
    try:
        return _pair(sp.factor_list(P), gens)
    except ValueError:
        pass
    K, fl0 = sp.factor_list(P, gaussian=True)
    def recip(g):
        Pg = sp.Poly(g, *gens); sub = {sp.conjugate(x_): 1 / x_ for x_ in gens}
        mono = sp.Mul(*[x_ ** Pg.degree(x_) for x_ in gens])
        return sp.expand(sp.conjugate(g).subs(sub) * mono)
    fl, pend = [], list(fl0)
    while pend:
        g, m = pend.pop(0)
        if len(sp.Poly(g, *gens).terms()) == 1:
            fl.append((g, m)); continue
        gr = recip(g)
        if sp.cancel(gr / g).is_number:
            fl.append((g, m)); continue
        for k_, (h, mh) in enumerate(pend):
            if sp.cancel(gr / h).is_number:
                assert mh == m; pend.pop(k_); fl.append((sp.expand(g * h), m)); break
        else:
            raise ValueError(f'unpaired factor {g}')
    return K, fl

def normT(T):
    """T {e: c} (e tuples of Fractions): (cen, k8, H) with c_{2cen−e} = ω conj(c_e), ω = i^k8 (μ = k8 π/4); h = e^{−iμ}Σ c e^{i(e−cen)·ang}"""
    nv = len(next(iter(T)))
    cen = tuple((min(e[k] for e in T) + max(e[k] for e in T)) / 2 for k in range(nv))
    om = None
    for e, c in T.items():
        e2 = tuple(2 * cen[k] - e[k] for k in range(nv))
        if e2 not in T: return None
        c2 = T[e2]; cc = (c[0], -c[1]); den = cc[0] ** 2 + cc[1] ** 2
        w = ((c2[0] * cc[0] + c2[1] * cc[1]) / den, (c2[1] * cc[0] - c2[0] * cc[1]) / den)
        if om is None: om = w
        elif om != w: return None
    if om not in UNITS: return None
    return cen, UNITS[om], {tuple(e[k] - cen[k] for k in range(nv)): c for e, c in T.items()}

def binom(H, k8, wlin, ends, nmin, C7):
    """binomial h = 2|c1| cos(ψ + β − μ), ψ = π w/N; wlin(Δ) -> (A, B, C) with w = A + B n + C j; ends = j-ends (j0, j1): j = j0 + j1 n.
    Returns (sign, edge-zero info) or ('FAIL', msg)."""
    (e1, c1), _ = sorted(H.items())
    re, im = c1
    if re != 0 and im != 0: raise ValueError('non-axis coefficient')
    beta = UNITS[(1 if re > 0 else -1 if re < 0 else 0, 1 if im > 0 else -1 if im < 0 else 0)]
    s8 = (2 * beta - k8) % 8
    A, B, C = wlin(e1)
    pq = [(A + C * j0, B + C * j1) for j0, j1 in ends]
    if all(p == 0 and q == 0 for p, q in pq): raise ValueError('constant binomial')
    tv = [Fr(0)] + [x_ for p, q in pq for x_ in ((p + q * nmin) / (ASL * nmin + C7), Fr(q, ASL))]
    tmin, tmax = min(tv), max(tv)
    zs = [Fr(1, 2) - Fr(s8, 4) + k for k in range(-4, 5)]
    if any(tmin < x_ < tmax for x_ in zs):
        return 'FAIL', f'zero inside: w-ends {pq}, s8 {s8}, t ∈ [{tmin}, {tmax}]'
    edge = [x_ for x_ in zs if x_ in (tmin, tmax)]
    attained = []
    for x_ in edge:
        if x_ == 0:
            for p, q in pq:
                if (p == 0 and q == 0) or (q != 0 and (-p / q).denominator == 1 and -p / q >= nmin):
                    attained.append('w = 0 at a j-end')
        elif any(x_ == (p + q * nmin) / (ASL * nmin + C7) for p, q in pq):
            attained.append(f't = {x_} at m = mmin, N = {ASL}m + {C7}')
    sg = 1 if math.cos(math.pi * (float((tmin + tmax) / 2) + s8 / 4)) > 0 else -1
    return sg, attained

def signed_factors(P, gens, wlin, ends, nmin, C7, src, denP=None):
    """sign bookkeeping of P/denP; returns dict(sign, tight, cores=[(H, k8, g)], fails=[...])"""
    K, fl = factor_real(P, gens)
    lists = [(fl, 1, K)]
    if denP is not None:
        Kd, fld = factor_real(denP, gens); lists.append((fld, -1, Kd))
    phase = [Fr(0)] * len(gens); k8tot = 0; sign = 1; tight = []; cores = []; fails = []
    Kc = complex(K) / (complex(lists[1][2]) if denP is not None else 1)
    for fac_list, ms, _ in lists:
        for g, m in fac_list:
            m = int(m)
            Pg = sp.Poly(g, *gens); T = {tuple(Fr(int(a)) for a in e): gi(c) for e, c in Pg.terms()}
            if len(T) == 1:
                (e1, c1), = T.items(); Kc *= complex(float(c1[0]), float(c1[1])) ** (ms * m)
                for k in range(len(gens)): phase[k] += ms * m * e1[k]
                continue
            nm = normT(T)
            if nm is None: raise ValueError(f'{src}: non-self-reciprocal factor {g}')
            cen, k8, H = nm
            for k in range(len(gens)): phase[k] += ms * m * cen[k]
            k8tot += ms * m * k8
            if len(H) == 2:
                r = binom(H, k8, wlin, ends, nmin, C7)
                if r[0] == 'FAIL': fails.append((str(g), r[1])); continue
                sg, att = r
                if att:
                    if ms < 0: raise ValueError(f'{src}: denominator vanishes at a lattice point')
                    tight.append((str(g), att))
                if m % 2: sign *= sg
            else:
                if ms < 0: raise ValueError(f'{src}: core in denominator')
                cores += [(H, k8, str(g))] * (2 if m % 2 == 0 else 1)   # even power: certified nonvanishing, contributes (+-1)^2
    assert all(p == 0 for p in phase), (src, phase)
    Kc *= cmath.exp(1j * math.pi / 4 * k8tot)
    assert abs(Kc.imag) < 1e-9 * abs(Kc), (src, Kc)
    sign *= 1 if Kc.real > 0 else -1
    return dict(sign=sign, tight=tight, cores=cores, fails=fails)

def stepA(rel, C7):
    out = []
    for (i, tag) in rel:
        f = FAMS[i]; name, c, d = f[0], f[1], f[2]; kind = name.split('.')[3]
        const = c == (0, 1) or kind == 'first'
        ends = [(Fr(0), Fr(0))] if const else [(Fr(1 if kind == 'rest' else 0), Fr(0)), (Fr(c[1] - 1), Fr(c[0]))]
        num, den = sp.fraction(sp.together(sp.expand(expr_of(f, tag))))
        if const: num, den = num.subs(u, 1), den.subs(u, 1)
        num = sp.expand(num)
        rec = dict(name=name, tag=tag, c=c, d=d, kind=kind)
        if num == 0:
            rec['status'] = 'ZERO'; out.append(rec); continue
        wlin = lambda e: (e[0] / 2, e[1], 2 * d * e[2])           # ψ/x for z^a v^b u^c offsets
        r = signed_factors(num, (z, v, u), wlin, ends, NMIN, C7, f'{name}:{tag}', denP=sp.expand(den))
        rec.update(sign=r['sign'], tight=r['tight'], fails=r['fails'], cores=r['cores'])
        out.append(rec)
    return out

# ---------------------------------------------------------------- (2) φ-cores -> (x, θ) functions
def to_ZV(H2):
    es = {}
    for (ez, ev), c in H2.items():
        assert (2 * ez).denominator == 1 and (2 * ev).denominator == 1
        es[(int(2 * ez), int(2 * ev))] = c
    mz = min(a for a, b in es); mv = min(b for a, b in es)
    P = sum((sp.Rational(c[0].numerator, c[0].denominator) + sp.I * sp.Rational(c[1].numerator, c[1].denominator)) * Z ** (a - mz) * V ** (b - mv)
            for (a, b), c in es.items() if c != (0, 0))
    return sp.expand(P), (mz, mv)

def reduce2d(H2, need, src, C7):
    """H2 {(ez, ev): c}: Σ c e^{i(ez x/2 + ev θ)} real. Returns atom dict or failure."""
    P, sh = to_ZV(H2)
    if P == 0: raise ValueError(f'{src}: identically zero')
    # exponents: Z^a V^b with ψ = x(a/4 + b n); monomial shift sh is folded into the phase check by multiplying back
    # the monomial shift Z^{mz} V^{mv} must cancel the factors' centres
    K, fl = factor_real(P, (Z, V))
    phase = [Fr(sh[0]), Fr(sh[1])]; k8tot = 0; sign = 1; tight = []; cores = []; fails = []; Kc = complex(K)
    for g, m in fl:
        m = int(m)
        Pg = sp.Poly(g, Z, V); T = {tuple(Fr(int(a)) for a in e): gi(c) for e, c in Pg.terms()}
        if len(T) == 1:
            (e1, c1), = T.items(); Kc *= complex(float(c1[0]), float(c1[1])) ** m
            phase[0] += m * e1[0]; phase[1] += m * e1[1]; continue
        nm = normT(T)
        if nm is None: raise ValueError(f'{src}: non-self-reciprocal factor {g}')
        cen, k8, H = nm
        phase[0] += m * cen[0]; phase[1] += m * cen[1]; k8tot += m * k8
        if len(H) == 2:
            rr = binom(H, k8, lambda e: (e[0] / 4, e[1] / 2, Fr(0)), [(Fr(0), Fr(0))], NMIN, C7)
            if rr[0] == 'FAIL': fails.append((str(g), rr[1])); continue
            sg, att = rr
            if att: tight.append((str(g), att))
            if m % 2: sign *= sg
        else:
            cores += [(H, k8, str(g))] * (2 if m % 2 == 0 else 1)
    assert phase == [0, 0], (src, phase)
    Kc *= cmath.exp(1j * math.pi / 4 * k8tot)
    assert abs(Kc.imag) < 1e-9 * abs(Kc), (src, Kc)
    sign *= 1 if Kc.real > 0 else -1
    return dict(src=src, req=need * sign, cores=cores, fails=fails, tight=tight)

ROT = {0: (1, 0), 2: (0, -1), 4: (-1, 0), 6: (0, 1)}            # e^{-iπ k8/4}, exact for even k8
def rotate(H, k8):
    r = ROT[k8 % 8]
    return {e: (cv[0] * r[0] - cv[1] * r[1], cv[0] * r[1] + cv[1] * r[0]) for e, cv in H.items()}

def core_sample_sign(Hr, d, c, kind):
    """sign of the real core at the lattice point m = 4, N = 59, j = middle of the run"""
    n, N = 4, 59; x = math.pi / N      # sample point m = 4, N = 59 (inside the 2D region)
    rr = range(1) if c == (0, 1) else range(c[0] * n + c[1])
    js = [0] if kind == 'first' else ([j for j in rr if j >= 1] if kind == 'rest' else list(rr))
    j = js[len(js) // 2]; th, ph = n * x, 2 * d * j * x
    val = sum(complex(float(cv[0]), float(cv[1])) * cmath.exp(1j * (float(e[0]) * x / 2 + float(e[1]) * th + float(e[2]) * ph))
              for e, cv in Hr.items())
    assert abs(val.imag) < 1e-9 * (1 + abs(val.real)) and val.real != 0
    return 1 if val.real > 0 else -1

def stepB(A, C7):
    """each core gets its own target sign σ_i (sampled); certified separately; the record needs Π σ_i = its binomial/constant sign."""
    atoms = []; rec_checks = []
    for rec in A:
        if rec.get('status') == 'ZERO' or not rec['cores']: continue
        c, d, kind = rec['c'], rec['d'], rec['kind']
        src0 = f"{rec['name']}:{rec['tag']}"; sig = []
        for ci, (H, k8, g) in enumerate(rec['cores']):
            src = src0 + (f':core{ci}' if len(rec['cores']) > 1 else '')
            assert k8 % 2 == 0, (src, k8)
            Hr = rotate(H, k8)
            need = core_sample_sign(Hr, d, c, kind); sig.append(need)
            us = set(e[2] for e in Hr)
            if us == {0}:
                atoms.append(reduce2d({(e[0], e[1]): cv for e, cv in Hr.items()}, need, src, C7)); continue
            assert us == {-1, 0, 1}, (src, us)
            jlo = 1 if kind == 'rest' else 0
            ends = {'lo': (Fr(2 * d * jlo), Fr(0)), 'hi': (Fr(2 * d * (c[1] - 1)), Fr(2 * d * c[0]))}    # φ = αx + βθ
            left = 'lo' if d > 0 else 'hi'
            for en, (al, be) in ends.items():
                Hv, Hd = {}, {}
                for (ez, ev, eu), cv in Hr.items():
                    key = (ez + 2 * eu * al, ev + eu * be)
                    a0 = Hv.get(key, (Fr(0), Fr(0))); Hv[key] = (a0[0] + cv[0], a0[1] + cv[1])
                    dc = (-eu * cv[1], eu * cv[0]); a1 = Hd.get(key, (Fr(0), Fr(0))); Hd[key] = (a1[0] + dc[0], a1[1] + dc[1])
                Hv = {k_: x_ for k_, x_ in Hv.items() if x_ != (0, 0)}; Hd = {k_: x_ for k_, x_ in Hd.items() if x_ != (0, 0)}
                atoms.append(reduce2d(Hv, need, src + f':val@{en}', C7))
                atoms.append(reduce2d(Hd, need if en == left else -need, src + f':dphi@{en}', C7))
        rec_checks.append((src0, math.prod(sig) == rec['sign']))
    return atoms, rec_checks

# ---------------------------------------------------------------- (3) 2D certificate
def ivq(q): return iv.mpf(q.numerator) / q.denominator

def sinc_iv(tt):
    a, b = tt.a, tt.b; mm = max(abs(a), abs(b))
    if mm <= 2:
        M = 12; acc = iv.mpf(0); t2 = tt * tt
        for j in range(M, -1, -1):
            acc = acc * t2 + iv.mpf((-1) ** j) / math.factorial(2 * j + 1)
        R = iv.mpf(mm) ** (2 * M + 2) / math.factorial(2 * M + 3)
        return acc + iv.mpf([-R.b, R.b])
    if a > 0 or b < 0: return iv.sin(tt) / tt
    raise ValueError('sinc straddle')

class Core:
    def __init__(s, H, k8, C7):
        rot = {0: (1, 0), 2: (0, -1), 4: (-1, 0), 6: (0, 1)}[k8 % 8]
        s.C7 = C7
        s.c = []
        for (ez, ev), cv in H.items():
            cr = (cv[0] * rot[0] - cv[1] * rot[1], cv[0] * rot[1] + cv[1] * rot[0])
            s.c.append((cr, ez / 4, ev / 2))        # (Z, V) exponents: angle = ez·x/4 + ev·θ/2 = θ(aξ + b), a = ez/4, b = ev/2
        s.cre = [(ivq(cr[0]), ivq(cr[1]), ivq(a), ivq(b)) for cr, a, b in s.c]
        s.k = 1 if sum(cr[0] for cr, a, b in s.c) == 0 else 0
    def theta(s, S, X): return (1 - S) * iv.pi / (ASL + s.C7 * X)
    def G(s, S, X):
        th = s.theta(S, X); acc = iv.mpf(0)
        for cr, ci, a, b in s.cre:
            L = a * X + b; y = th * L
            if s.k == 0: acc += cr * iv.cos(y) - ci * iv.sin(y)
            else: acc += L * sinc_iv(y / 2) * (-cr * iv.sin(y / 2) - ci * iv.cos(y / 2))
        return acc
    def dG(s, S, X):
        th = s.theta(S, X); ht = iv.mpf(0); hx = iv.mpf(0)
        for cr, ci, a, b in s.cre:
            L = a * X + b; y = th * L; rie = -cr * iv.sin(y) - ci * iv.cos(y)
            ht += L * rie; hx += th * a * rie
        ts = -iv.pi / (ASL + s.C7 * X); tx = -(1 - S) * s.C7 * iv.pi / (ASL + s.C7 * X) ** 2
        return ht * ts, ht * tx + hx
    def corner_zero(s):
        """EXACT: Re Σ c e^{iπ b/ASL} == 0, reduced mod Φ_{4·ASL}(t), t = e^{iπ/(2·ASL)} (b ∈ ½ℤ)"""
        t = sp.symbols('t'); K = 4 * ASL
        terms = []
        for cr, a, b in s.c:
            cc = sp.Rational(cr[0].numerator, cr[0].denominator) + sp.I * sp.Rational(cr[1].numerator, cr[1].denominator)
            e2 = 2 * b; assert e2.denominator == 1
            terms += [(cc, int(e2)), (sp.conjugate(cc), -int(e2))]
        sh = -min(e for _, e in terms)
        P = sp.expand(sum(c_ * t ** (e + sh) for c_, e in terms))
        return sp.rem(P, sp.cyclotomic_poly(K, t), t) == 0

def certify(core, sgn, maxdep=26):
    corner0 = core.k == 0 and core.corner_zero()
    stack = [(mpf(0), mpf(1), mpf(0), mpf(1) / NMIN, 0)]; nb = 0; lowest = None; cinfo = None
    while stack:
        s1, s2, x1, x2, dep = stack.pop(); nb += 1
        S = iv.mpf([s1, s2]); X = iv.mpf([x1, x2])
        if corner0 and s1 == 0 and x1 == 0:
            ds, dx = core.dG(S, X); Aa = (sgn * ds).a; Bb = (sgn * dx).a
            if Aa > 0 and Bb > 0:
                cinfo = dict(s2=float(s2), x2=float(x2), A=float(Aa), B=float(Bb)); continue
        else:
            lo = (sgn * core.G(S, X)).a
            if core.k == 0 and lo <= 0:
                Sc = iv.mpf((s1 + s2) / 2); Xc = iv.mpf((x1 + x2) / 2); ds, dx = core.dG(S, X)
                lo = max(lo, (sgn * (core.G(Sc, Xc) + ds * (S - Sc) + dx * (X - Xc))).a)
            if lo > 0:
                lowest = lo if lowest is None or lo < lowest else lowest; continue
        if dep > maxdep: return dict(ok=False, boxes=nb, where=[float(s1), float(s2), float(x1), float(x2)])
        if (s2 - s1) >= NMIN * (x2 - x1):
            m = (s1 + s2) / 2; stack += [(s1, m, x1, x2, dep + 1), (m, s2, x1, x2, dep + 1)]
        else:
            m = (x1 + x2) / 2; stack += [(s1, s2, x1, m, dep + 1), (s1, s2, m, x2, dep + 1)]
    return dict(ok=True, boxes=nb, lb=float(lowest) if lowest is not None else None, corner=cinfo, corner_zero=corner0)

def run(C7, rel):
    t0 = time.time()
    A = stepA(rel, C7)
    atoms, rec_checks = stepB(A, C7)
    distinct = {}
    for a in atoms:
        for (H, k8, g) in a['cores']:
            distinct.setdefault((g, k8), (H, k8))
    cert = {}
    for key, (H, k8) in distinct.items():
        core = Core(H, k8, C7)
        v0 = core.G(iv.mpf(0.5), iv.mpf(0.03)); sgn = 1 if v0.a > 0 else -1 if v0.b < 0 else 0
        assert sgn != 0
        r = certify(core, sgn); r.update(sign=sgn, k=core.k); cert[key] = r
    # (4) assembly -> flagged expressions
    flagged = set(); why = collections.defaultdict(list)
    for rec in A:
        if rec.get('status') == 'ZERO': continue
        key = f"{rec['name']}:{rec['tag']}"
        if rec['fails']: flagged.add(key); why[key].append('binomial (A)')
        if not rec['cores'] and rec['sign'] != 1 and not rec['fails']: flagged.add(key); why[key].append('sign (A)')
    for src0, okp in rec_checks:
        if not okp: flagged.add(src0); why[src0].append('core sign product')
    for a in atoms:
        key = ':'.join(a['src'].split(':')[:2])
        if a['fails']: flagged.add(key); why[key].append('binomial (B)'); continue
        prod = 1; okc = True
        for (H, k8, g) in a['cores']:
            r = cert[(g, k8)]
            if not r['ok']: okc = False
            prod *= r['sign']
        if not okc: flagged.add(key); why[key].append('2D core'); continue
        if prod != a['req']: flagged.add(key); why[key].append('sign (B)')
    return dict(C7=C7, secs=round(time.time() - t0), n_expr=len(A), zero=sum(1 for r in A if r.get('status') == 'ZERO'),
                core_free=sum(1 for r in A if r.get('status') != 'ZERO' and not r['cores']), atoms=len(atoms),
                distinct_cores=len(cert), cores_ok=sum(1 for r in cert.values() if r['ok']),
                corner_cores=sum(1 for r in cert.values() if r.get('corner')),
                min_lb=min((r['lb'] for r in cert.values() if r.get('lb')), default=None),
                max_boxes=max(r['boxes'] for r in cert.values()),
                corner_tests=[r['corner'] for r in cert.values() if r.get('corner')],
                tight=sorted({f"{r['name']}:{r['tag']}" for r in A if r.get('tight')}),
                zero_names=sorted(f"{r['name']}:{r['tag']}" for r in A if r.get('status') == 'ZERO'),
                flagged=sorted(flagged), why={k_: v_ for k_, v_ in why.items()})

def numeric_violations(rel, off):
    bad = set()
    for (i, tag) in rel:
        f = FAMS[i]; fn = sp.lambdify((z, v, u), expr_of(f, tag), 'mpmath')
        for n in (4, 5, 7):
            N = ASL * n + off; zz = mp.exp(1j * mp.pi / (2 * N))
            for j in reps_of(f, n):
                if mp.re(fn(zz, zz ** (2 * n), zz ** (4 * f[2] * j))) < -1e-25: bad.add(f'{f[0]}:{tag}')
    return bad

def per_m(M, off=1):
    """every family at m = M, every repetition j: Laurent in z, certified >= 0 on x ∈ (0, X], X = ⌈π/(14M+1)⌉ (7 decimals)"""
    import pi_n_torus as KP
    KP.XMAX = Fr(math.ceil(math.pi / (ASL * M + off) * 10 ** 7), 10 ** 7)
    assert float(KP.XMAX) >= math.pi / (ASL * M + off)
    xn = math.pi / (ASL * M + off)
    def toL(e):
        e = sp.cancel(sp.together(e)); num, den = sp.fraction(e)
        dP = sp.Poly(den, z); assert len(dP.terms()) == 1, den
        (dd,), dc = dP.terms()[0]
        co = {int(k[0]) - int(dd): sp.nsimplify(c / dc) for k, c in sp.Poly(sp.expand(num), z).terms()}
        dl = sp.ilcm(*[sp.fraction(sp.nsimplify(sp.re(c)))[1] for c in co.values()],
                     *[sp.fraction(sp.nsimplify(sp.im(c)))[1] for c in co.values()], 1)
        return KP.L({k: complex(int(sp.re(c) * dl), int(sp.im(c) * dl)) for k, c in co.items()})
    res = dict(M=M, off=off, XMAX=str(KP.XMAX), tight=0, certified=0, orders=set(), minlb=None, fails=[], numneg=[])
    x0 = math.pi / (ASL * M + 3)
    for f in FAMS:
        c, d = f[1], f[2]
        for j in (range(1) if c == (0, 1) else range(c[0] * M + c[1])):
            sub = {v: z ** (2 * M), u: z ** (4 * d * j)}
            if f[4] is None:
                exprs = [('sign', f[3])]
            else:
                p, q = toL(f[5].subs(sub)), toL(f[6].subs(sub))
                exprs = ([('a-p', f[3] - f[5]), ('q-b', f[6] - f[4])] if p.evx(x0).real <= q.evx(x0).real
                         else [('a-q', f[3] - f[6]), ('p-b', f[5] - f[4])])
            for tag, e in exprs:
                L = toL(e.subs(sub))
                if L.iszero():
                    res['tight'] += 1; continue
                if L.evx(xn).real < -1e-12:
                    res['numneg'].append((f[0], tag, j))
                try:
                    n0, lb = KP.certify_nonneg(L)
                    res['certified'] += 1; res['orders'].add(n0); res['minlb'] = lb if res['minlb'] is None else min(res['minlb'], lb)
                except AssertionError as ex:
                    res['fails'].append((f[0], tag, j, str(ex)))
    res['orders'] = sorted(res['orders'])
    res['fails_eq_numneg'] = sorted(a[:3] for a in res['fails']) == sorted(res['numneg'])
    return res

if __name__ == '__main__':
    t0 = time.time()
    rel = orient()
    main = run(1, rel)
    print('MAIN', {k_: v_ for k_, v_ in main.items() if k_ not in ('why', 'tight', 'corner_tests')}, flush=True)
    out = dict(main=main)
    if '--no-control' not in sys.argv:
        ctrl = run(-1, rel)
        violm = numeric_violations(rel, -1); violp = numeric_violations(rel, 1)
        out['control'] = dict(flagged=ctrl['flagged'], numeric_neg_14m_minus1=sorted(violm), numeric_neg_14m_plus1=sorted(violp),
                              exact_match=set(ctrl['flagged']) == violm, why=ctrl['why'])
        print('CONTROL flagged', len(ctrl['flagged']), 'numeric 14m-1', len(violm), 'match', set(ctrl['flagged']) == violm,
              ' numeric 14m+1', len(violp), flush=True)
    out['per_m'] = {M: per_m(M) for M in (2, 3)}
    print('PER-m', {M: {k_: r[k_] for k_ in ('certified', 'tight', 'fails')} for M, r in out['per_m'].items()}, flush=True)
    if '--no-control' not in sys.argv:
        # audit G6: controls for the per-m and the run (sinusoid ⊕ binomial j-closure) paths
        pc = {M: per_m(M, off=-1) for M in (2, 3)}
        out['control_per_m'] = {M: dict(fails=len(r['fails']), numneg=len(r['numneg']), exact=r['fails_eq_numneg']) for M, r in pc.items()}
        orig = list(FAMS)
        FAMS[:] = [(f[0], f[1] if f[1] == (0, 1) else (f[1][0], f[1][1] + 1)) + tuple(f[2:]) for f in orig]   # every run one repetition longer
        negx = numeric_violations(rel, 1); rx = run(1, rel)
        FAMS[:] = orig
        out['control_run_extension'] = dict(numeric_neg=sorted(negx), flagged=rx['flagged'], why=rx['why'],
                                            exact=set(rx['flagged']) == negx and bool(negx))
        print('CONTROL per-m (N >= 14M-1):', out['control_per_m'], ' run-extension K=1: numeric', len(negx), 'flagged',
              len(rx['flagged']), 'exact', out['control_run_extension']['exact'], flush=True)
    out['verdict'] = bool(not main['flagged'] and main['cores_ok'] == main['distinct_cores']
                          and all(not r['fails'] for r in out['per_m'].values())
                          and (out.get('control_per_m') is None or all(c['exact'] and c['fails'] for c in out['control_per_m'].values()))
                          and (out.get('control_run_extension') is None or out['control_run_extension']['exact'])
                          and (out.get('control') is None or (out['control']['exact_match'] and out['control']['flagged']
                                                              and not out['control']['numeric_neg_14m_plus1'])))
    out['secs'] = round(time.time() - t0)
    json.dump(out, open(OUTF, 'w'), indent=1, default=str)
    print('VERDICT', out['verdict'], out['secs'], 's →', OUTF)
