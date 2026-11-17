"""[Deposit] Proposition A.1: the two-variable certificate for m ≥ 7 and distinctness. Paper: "Perpendicular billiards in rational polygons of odd denominator". Research name: probes/s937_ladder_allm_cert.py, copied verbatim except imports (deposit names), paths (data/ and figs/ under odd_polygons/) and this paragraph. Run through run_all.py; the usage lines below are the research repo's. Bracketed IDs and section marks (§sNNN) refer to the project's internal ledgers.

s937 keeper (queue `s935-a`, the LOWER bound): THE m-TH TORUS OF Π_N EXISTS FOR EVERY m ≥ 7 AND EVERY ODD N ≥ 8m + 15 —
the containment inequalities of §s935 (4)'s template, certified for all (m, N) at once; with the per-m certificates (m ≤ 6) and the
formal identities (all m) this gives 2⌊(N − 7)/8⌋ distinct minimal components of Π_N's perpendicular class.

NOT PRE-REGISTERED: a certificate, not a hypothesis test. Inputs: the symbolic inequality families of `probes/s935_torus_ladder.py`
(`data/s935_ladder_fams.pkl`: Laurent expressions in z = ζ_{4N} = e^{ix/2}, v = z^{4n} = e^{iθ}, u = z^{4dj} = e^{iφ}; x = π/N,
θ = 2nx, φ = 2djx, j the repetition inside a run). n = m + 1 ≥ 8, N odd ≥ 8n + 7.
METHOD (every step exact except (3), which is rigorous interval arithmetic):
 (0) orientation of each side-range pair (p, q) read at lattice points — soundness does not depend on it (a − p ≥ 0, q − b ≥ 0 and
     a < b force p ≤ a < b ≤ q).
 (1) each expression = K·Π g^m / den, factored over ℤ (ℤ[i] with conjugate-reciprocal pairing when needed); every factor is
     self-reciprocal up to a unit ω ∈ {±1, ±i}, so g = monomial·e^{iμ}·h with h real; the total phase must cancel and K e^{iΣμ} be
     real (asserted). A BINOMIAL h is 2|c|cos(π(w/N + s/4)) with w = A + Bn + Cj linear: its sign over the lattice is decided
     exactly from the closure of {w/N} (no zero of cos inside; zeros on the boundary are recorded — attained ones are tight strips).
 (2) a CORE (≥ 3 terms) depends on (x), (x, θ) or (x, θ, φ); in φ it is a + 2Re(b e^{iφ}), a sinusoid on the φ-interval of the
     run (length ≤ 4θ < π). Lemma: positive at both ends, non-decreasing at the left end, non-increasing at the right end ⟹ positive
     on the interval. So each φ-core becomes four (x, θ) functions (value ⊕ φ-derivative at both ends), each re-factored as in (1).
 (3) the remaining 56 distinct residual cores G(θ, ξ) = h(θξ, θ)/θ^k (k ∈ {0, 1}, ξ = x/θ = 1/(2n)) are certified to have constant
     sign on Ω = {ξ ∈ [0, 1/16], θ ∈ (0, π/(4 + 7ξ)]} in coordinates θ = (1 − s)π/(4 + 7ξ), (s, ξ) ∈ [0, 1] × [0, 1/16] — Ω exactly,
     no overshoot of the N = 8n + 7 edge. k = 1 via h/θ = Re Σ c·iL·e^{iy/2}·sinc(y/2) (Re Σ c = 0 exact). The only zeros of a
     core on the closed Ω are at the corner (s, ξ) = (0, 0) (n = ∞ on N = 8n + 7; checked exactly), discharged by a gradient test:
     sgn·∂_sG > 0 and sgn·∂_ξG > 0 on the corner box ⟹ sgn·G ≥ sA + ξB > 0. Elsewhere: interval boxes (naive ∩ mean-value),
     adaptive bisection, mpmath iv.
 (4) assembly: for every expression the product of the certified signs equals the required sign.
CONTROLS: (a) the whole pipeline with the region widened to N ≥ 8n + 5 must fail, and must fail on EXACTLY the expressions that
are numerically negative at N = 8n + 5 (n = 8, 9, 12) — the s935 per-m controls (19 + 12m violations) extended; (b) the families
reproduce the per-m inequality lists of `data/s935_tori_63.json` (flow_poly-read itineraries, m = 1 … 6) as value multisets.
DISTINCTNESS: the copies (o, e) met by the m-th torus have o odd (start 2n − 1, every side shift DM even) in [−2n − 1, 2n + 3]
(max of the symbolic walk's linear forms); the time-reversing deck map h(p) = p̄, (o, e) ↦ (1 − o, 1 − e), sends them to even o
— disjoint mod N (window 4n + 5 < N); max offset 2n + 3 separates different m. Copy sets are component invariants.
PRIOR ART: §s935 (1) (m = 0, [POLY-ASYM-ALL]), §s935 (4) (identities all m, per-m certificates m ≤ 6, "all-m rigor needs a 3D
certificate in (x, 2nx, 2djx)"), §s935 (2) (the pairing "presumably the deck reflection in V0V1's line … not checked"); rulings
grep `torus ladder`, `Π_N` → [LIT-058] only (lit, no certificate). DATA: `data/s935_ladder_fams.pkl`, `data/s935_tori_63.json`.
→ general_triangles.md §s937.   Usage: PYTHONPATH=… python probes/s937_ladder_allm_cert.py [--no-control | --distinct-only]   (≈ 17 min; 8 without control)"""
import _setup  # noqa: F401  (deposit: this directory on the path, cwd = odd_polygons/)
import os, sys, json, math, cmath, pickle, time, collections
import sympy as sp
import mpmath as mp
from mpmath import iv, mpf
from fractions import Fraction as Fr

ROOT = '.'   # deposit: cwd = odd_polygons/ (_setup)
OUTF = f'{ROOT}/data/s937_ladder_allm_cert.json'
NMIN = 8
z, v, t, u = sp.symbols('z v t u')
Z, V = sp.symbols('Z V')
F0, F1, *_ = pickle.load(open(f'{ROOT}/data/s935_ladder_fams.pkl', 'rb'))
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
        for n, N in ((8, 71), (9, 81), (12, 1001)):
            zz = mp.exp(1j * mp.pi / (2 * N))
            for j in reps_of(f, n):
                sg.add(mp.re(fn(zz, zz ** (4 * n), zz ** (4 * f[2] * j))) > 0)
        assert len(sg) == 1, f[0]
        rel += [(i, 'a-p'), (i, 'q-b')] if sg.pop() else [(i, 'a-q'), (i, 'p-b')]
    return rel

# ---------------------------------------------------------------- (1) factor, normalize, binomials
def factor_real(P, gens):
    PP = sp.Poly(P, *gens)
    if str(PP.domain) == 'ZZ':
        return sp.factor_list(P)
    K, fl0 = sp.factor_list(P, gaussian=True)
    def recip(g):
        Pg = sp.Poly(g, *gens); sub = {sp.conjugate(x_): 1 / x_ for x_ in gens}
        mono = sp.Mul(*[x_ ** Pg.degree(x_) for x_ in gens])
        return sp.expand(sp.conjugate(g).subs(sub) * mono)
    fl, pend = [], list(fl0)
    while pend:
        g, m = pend.pop(0); gr = recip(g)
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
    tv = [Fr(0)] + [x_ for p, q in pq for x_ in ((p + q * nmin) / (8 * nmin + C7), q / 8)]
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
        elif any(x_ == (p + q * nmin) / (8 * nmin + C7) for p, q in pq):
            attained.append(f't = {x_} at n = nmin, N = 8n + {C7}')
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
                if m % 2 == 0: raise ValueError(f'{src}: even core')
                cores.append((H, k8, str(g)))
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
        wlin = lambda e: (e[0] / 2, 2 * e[1], 2 * d * e[2])           # ψ/x for z^a v^b u^c offsets
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
            rr = binom(H, k8, lambda e: (e[0] / 4, e[1], Fr(0)), [(Fr(0), Fr(0))], NMIN, C7)
            if rr[0] == 'FAIL': fails.append((str(g), rr[1])); continue
            sg, att = rr
            if att: tight.append((str(g), att))
            if m % 2: sign *= sg
        else:
            if m % 2 == 0: raise ValueError('even core')
            cores.append((H, k8, str(g)))
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
    """sign of the real core at the lattice point n = 8, N = 71, j = middle of the run"""
    n, N = 8, 71; x = math.pi / N
    rr = range(1) if c == (0, 1) else range(c[0] * n + c[1])
    js = [0] if kind == 'first' else ([j for j in rr if j >= 1] if kind == 'rest' else list(rr))
    j = js[len(js) // 2]; th, ph = 2 * n * x, 2 * d * j * x
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
            ends = {'lo': (Fr(2 * d * jlo), Fr(0)), 'hi': (Fr(2 * d * (c[1] - 1)), Fr(d * c[0]))}    # φ = αx + βθ
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
    def theta(s, S, X): return (1 - S) * iv.pi / (4 + s.C7 * X)
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
        ts = -iv.pi / (4 + s.C7 * X); tx = -(1 - S) * s.C7 * iv.pi / (4 + s.C7 * X) ** 2
        return ht * ts, ht * tx + hx
    def corner_zero(s):
        val = sum((sp.Rational(cr[0].numerator, cr[0].denominator) + sp.I * sp.Rational(cr[1].numerator, cr[1].denominator))
                  * sp.exp(sp.I * sp.pi * sp.Rational(b.numerator, b.denominator) / 4) for cr, a, b in s.c)
        return sp.simplify(sp.re(sp.expand_complex(val))) == 0

def certify(core, sgn, maxdep=26):
    corner0 = core.k == 0 and core.corner_zero()
    stack = [(mpf(0), mpf(1), mpf(0), mpf(1) / 16, 0)]; nb = 0; lowest = None; cinfo = None
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
        if (s2 - s1) >= 16 * (x2 - x1):
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
        for n in (8, 9, 12):
            N = 8 * n + off; zz = mp.exp(1j * mp.pi / (2 * N))
            for j in reps_of(f, n):
                if mp.re(fn(zz, zz ** (4 * n), zz ** (4 * f[2] * j))) < -1e-25: bad.add(f'{f[0]}:{tag}')
    return bad

def generator_validation():
    sys.path.insert(0, f'{ROOT}/probes')
    import pi_n_torus as KP
    from pi_n_torus import L, X, g_apply, side_range, shift, V as VV, SIDES
    D = json.load(open(f'{ROOT}/data/s935_tori_63.json')); sc = D['scale']
    lam = [(f, [sp.lambdify((z, v, u), e, 'mpmath') if e is not None else None for e in f[3:7]]) for f in FAMS]
    N1 = 211; x1 = math.pi / N1; res = {}
    for m, tt in enumerate(D['tori']):
        if m == 0: continue
        n = m + 1; h0 = tt['h'][0]; rows = tt['rows']
        st = [(r['o'], r['e'], r['s'], r['s2']) for r in rows]; J0, J1 = st[:h0], st[h0:]
        def walk(seg):
            S, out = L(), []
            for (o, e, s, s2) in seg: out.append(S); S = S + shift((o, e), s2)
            return out, S
        def pin(rs, seg, hi):
            Ss, _ = walk(seg)
            for k, r in enumerate(rs):
                uu = (r['u'] + (r['w'] if hi else 0)) / sc
                for vi in range(4):
                    if abs(r['vx'][vi] / sc - uu) < 1e-9: return X(g_apply(seg[k][:2], VV[vi])) + Ss[k]
        lo0 = pin(rows[:h0], J0, False); hi0 = pin(rows[:h0], J0, True); hi1 = pin(rows[h0:], J1, True)
        Avals = []
        for seg, lo, la in ((J0, lo0, hi0 - lo0), (J1, hi0, hi1 - hi0)):
            Ss, _ = walk(seg)
            for k, (o, e, s, s2) in enumerate(seg):
                g = (o, e); a, b = lo - Ss[k], lo - Ss[k] + la
                for ss in (s, s2):
                    p, q = side_range(g, ss)
                    if p.evx(x1).real > q.evx(x1).real: p, q = q, p
                    Avals += [(a - p).evx(x1).real, (q - b).evx(x1).real]
                for ss, want in ((s, 1), (s2, -1)):
                    i_, kk = SIDES[ss]; ev = g_apply(g, VV[kk]) - g_apply(g, VV[i_])
                    Avals.append(((ev + ev.conj()) * ((1 if e == 0 else -1) * want)).evx(x1).real)
        Bvals = []; zz = mp.exp(1j * mp.pi / (2 * N1))
        for f, fn in lam:
            for j in reps_of(f, n):
                ev_ = [float(mp.re(g_(zz, zz ** (4 * n), zz ** (4 * f[2] * j)))) if g_ else None for g_ in fn]
                if f[4] is None: Bvals.append(ev_[0])
                else:
                    a, b, p, q = ev_
                    if p > q: p, q = q, p
                    Bvals += [a - p, q - b]
        Avals.sort(); Bvals.sort()
        res[m] = dict(n_perm=len(Avals), n_fam=len(Bvals),
                      maxdev=max(abs(a - b) for a, b in zip(Avals, Bvals)) if len(Avals) == len(Bvals) else None)
    return res

def distinctness():
    import pi_n_torus as KP
    src = open(f'{ROOT}/appA_1_identities.py').read().split('st0 = walk(J0)')[0]
    ns = {}; exec(src, ns)
    J0, J1, START, g_across, DM = ns['J0'], ns['J1'], ns['START'], ns['g_across'], ns['DM']
    assert all(d_ % 2 == 0 for d_ in DM) and (START[0][0] * 1 + START[0][1]) % 2 == 1 and START[0][0] % 2 == 0
    # symbolic extremes: offsets visited are o0 + partial + j·d (j < γn + δ), all linear in n -> list of (slope, intercept)
    forms = []
    for segs in (J0, J1):
        o = list(START[0]); e = START[1]
        for w, c in segs:
            # unit: partial offsets (relative) and net change d
            part = []; oo = 0; ee = e
            for ch in w:
                part.append(oo); (x_, ee) = g_across((0, oo), ee, int(ch)); oo = x_[1]
            dnet = oo
            assert ee == e or c == (0, 1)
            for pp in part:
                forms.append((o[0], o[1] + pp))                                              # j = 0
                if c != (0, 1): forms.append((o[0] + dnet * c[0], o[1] + pp + dnet * (c[1] - 1)))   # j = γn + δ − 1
            if c == (0, 1): o = [o[0], o[1] + dnet]; e = ee
            else: o = [o[0] + dnet * c[0], o[1] + dnet * c[1]]
    # max = 2n + 3, min = −2n − 1 for every n ≥ NMIN (forms at j = γn + δ − 1 are meaningless where a run is empty, n small)
    okmax = all(a <= 2 and a * NMIN + b <= 2 * NMIN + 3 for a, b in forms) and (2, 3) in forms
    okmin = all(a >= -2 and a * NMIN + b >= -2 * NMIN - 1 for a, b in forms) and (-2, -1) in forms
    # small n (the per-m tori m ≤ 6, n ≤ 7) by explicit grammar expansion: offsets odd, range exactly [−2n − 1, 2n + 3]  (audit G11)
    small = {}
    for n in range(1, NMIN):
        S = set()
        for segs in (J0, J1):
            o = START[0][0] * n + START[0][1]; e = START[1]
            for w, c in segs:
                for _ in range(c[0] * n + c[1]):
                    for ch in w:
                        S.add((o, e)); (x_, e) = g_across((0, o), e, int(ch)); o = x_[1]
        os_ = [o for o, e in S]
        small[n] = all(o % 2 == 1 for o in os_) and min(os_) == -2 * n - 1 and max(os_) == 2 * n + 3
    return dict(DM=DM, start='2n − 1 (odd)', max_is_2n_plus_3=okmax and all(small.values()), min_is_minus_2n_minus_1=okmin,
                small_n_expansion=small, n_forms=len(forms))

def lambdas_positive():
    """audit G13: λ0 = 4 sin 2θ sin x (2cos 2x − 1), λ1 = 16 sin 2θ sin³x (formal, §s935 (4)) are > 0: 2θ ≤ 2π/(4 + 7ξ) < π/2, x ≤ π/71 < π/6."""
    return (2 * math.pi / 4 < math.pi) and (math.pi / 71 < math.pi / 6)

def verdict(out):
    main = out['main']
    return bool(not main['flagged'] and main['cores_ok'] == main['distinct_cores']
                and all(r['maxdev'] is not None and r['maxdev'] < 1e-12 for r in out['generator_validation'].values())
                and out['distinctness']['max_is_2n_plus_3'] and out['distinctness']['min_is_minus_2n_minus_1']
                and (out.get('control') is None or (out['control']['exact_match'] and not out['control']['numeric_neg_8n7']
                                                    and out['control']['flagged'])))

if __name__ == '__main__' and '--distinct-only' in sys.argv:       # recompute distinctness ⊕ verdict into an existing output
    out = json.load(open(OUTF)); out['distinctness'] = distinctness(); out['verdict'] = verdict(out)
    json.dump(out, open(OUTF, 'w'), indent=1, default=str); print('DISTINCT', out['distinctness'], 'VERDICT', out['verdict']); sys.exit()

if __name__ == '__main__':
    t0 = time.time()
    rel = orient()
    main = run(7, rel)
    print('MAIN', {k_: v_ for k_, v_ in main.items() if k_ not in ('why', 'tight', 'corner_tests')}, flush=True)
    out = dict(main=main)
    if '--no-control' not in sys.argv:
        ctrl = run(5, rel)
        viol5 = numeric_violations(rel, 5); viol7 = numeric_violations(rel, 7)
        out['control'] = dict(flagged=ctrl['flagged'], numeric_neg_8n5=sorted(viol5), numeric_neg_8n7=sorted(viol7),
                              exact_match=set(ctrl['flagged']) == viol5, why=ctrl['why'])
        print('CONTROL flagged', len(ctrl['flagged']), 'numeric 8n+5', len(viol5), 'match', set(ctrl['flagged']) == viol5,
              ' numeric 8n+7', len(viol7), flush=True)
    out['generator_validation'] = generator_validation(); print('GENVAL', out['generator_validation'], flush=True)
    out['distinctness'] = distinctness(); print('DISTINCT', out['distinctness'], flush=True)
    out['lambdas_positive'] = lambdas_positive()
    out['verdict'] = verdict(out)
    out['secs'] = round(time.time() - t0)
    json.dump(out, open(OUTF, 'w'), indent=1, default=str)
    print('VERDICT', out['verdict'], f"{out['secs']}s →", OUTF)
