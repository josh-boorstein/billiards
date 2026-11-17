"""[Deposit] Helper, imported by the scripts of this directory. Paper: "Perpendicular billiards in rational polygons of odd denominator". Research name: probes/s938_q5_loops.py, copied verbatim except imports (deposit names), paths (data/ and figs/ under odd_polygons/) and this paragraph. Run through run_all.py; the usage lines below are the research repo's. Bracketed IDs and section marks (§sNNN) refer to the project's internal ledgers.

s938 `s936-d`: the V2-LOOPS of §s911's q = 5 family T_m = (3, 12m+5, 18m+7)/(30m+15) — for EVERY m ≥ 2, the m loops
L_{m,k} (k = 0 … m−1) exist with an explicit word, their lengths are a 40-term closed form, and they are STRICTLY increasing
in k with gap ℓ(m,k+1) − ℓ(m,k) = −½ sin(3x)·K(x, kx) > 0, x = π/(60m+30) — so no two loops ever coincide (no exact zero);
every gap is ≥ 57.375·x·sin 3x and the SMALLEST (k = m − 2) is ~ 0.6075/m² (small-k gaps are order 1/m).
AUDITED `audits/s938_audit_q5loops.md` (21 V / 3 G / 0 F / 3 W, repaired here: outward region bounds G22/G23, span assert G7). Computer-assisted: exact formal identities ⊕ rigorous interval certificates.

WORD (measured m ≤ 8 against the walk, then PROVED legal for all m): from V2, launch slot 78m + 55 + 6k (units π/2N), runs of
alternating reflections in sides 0, 2 (both through V0 = 0, so each run is linear) of lengths
[4m−1−2k, 8m+2−2k, 8m+6+2k, 4m+4+2k, 1+2k], first letters 0, 2, 0, 0, 0, separated by four reflections in side 1.
LENGTH: 4ℓ = Σ_40 ±cos(π(αm + βk + γ)/2N) (exact, `S_formal`); with y = kx each term is cos(απ/60 + βy + (γ − α/2)x),
α ∈ 12ℤ — the π/5 of "q = 5".
LEGALITY (the developed line realises the word, for all (m, k)): in polygon i (the V0-star after i side-1 reflections) the
crossed rays of each parity form an arc — L cos ψ > d is concave in ψ and every run spans < π (`span_ok`, audit G7) — so it suffices that the EXTREME crossed ray ends of each parity lie beyond the line
(opposite the star centre c_i), the ray after the last lies on the centre side (exit through side 1), the launch points into
T, and run 5 ends exactly at V2 (an IDENTITY: every (β, δ)-frequency group vanishes in ℚ(ζ_240)). Each condition f(x, y)
is certified on x ∈ (0, π/150], x ≤ y ≤ π/60 − 1.5x (k ≥ 1) and y = 0 (k = 0): a_0 ≡ 0 checked exactly → f/x via
sinc/versine enclosures; degenerate corners (a_1 = 0 at y = 0, or a_0 = 0 at y = π/60 for c_3) by ratio charts.
GAP: G = −½ sin(3x) K, K = Σ cσ sin(φ + λu) with Σ cσ sin φ = 0 EXACTLY (ℚ(ζ_120)), H = K/u < 0 on ρ = x/u ∈ [0, 0.4],
u = π/60 − y ∈ (0, π/60] (sup H ≤ −45.9); CONTROL: ρ ≤ 0.8 does NOT certify (H(ρ, 0) = −147.7 + 221.6ρ).

SCOPE. Proved: the m loops exist and are pairwise distinct with the stated gap, every m ≥ 2. NOT proved: that they are ALL the
swapped connections of T_m (p = m is the walk's count, m ≤ 8 here, ≤ 5 in §s936 (4)); signed relations Σ ε_i ℓ_i = 0 with
more than two terms are not addressed; the twin family (3, 12m+4, 18m+8) is not treated.
PRIOR ART: §s936 (4) (the near-coincidences 1.6·10⁻² … 3.3·10⁻³, "never 0 at these m"), §s911 (the family), §s935 (2) / §s937
(Π_N's loop ℓ_m = c_{4m+5} + c_{4m+7} and its all-m certificate — the template); `rulings.py --grep` *q = 5 family*, *loop gap*,
*near-coincidence* → none. DATA: `data_index.py --grep loop` → s935/s937 Π_N tori only → none for this family.
Output: `data/s938_q5_loops.json`. Runtime ~ 2 min."""
import _setup  # noqa: F401  (deposit: this directory on the path, cwd = odd_polygons/)
import json
import sys
from fractions import Fraction

import mpmath as mp
from mpmath import iv
import sympy as sp

mp.mp.dps = 50
iv.dps = 30
OUT = 'data/s938_q5_loops.json'


# ---------------- formal monomials ζ^E, E = (α, β, γ) ↦ αm + βk + γ, ζ = ζ_{4N}, N = 30m + 15 ----------------
def add(u, v, s=1):
    r = dict(u)
    for e, c in v.items():
        r[e] = r.get(e, 0) + s * c
    return {e: c for e, c in r.items() if c}


def sh(u, E):
    return {tuple(a + b for a, b in zip(e, E)): c for e, c in u.items()}


def cj(u):
    return {tuple(-a for a in e): c for e, c in u.items()}


def neg(E):
    return tuple(-a for a in E)


def lin(E1, E2):
    return tuple(a + b for a, b in zip(E1, E2))


def mulc(E, n):
    return tuple(n * a for a in E)


def canon(f):                                       # exponents mod 4N = (120, 0, 60)
    out = {}
    for (a, b, g), c in f.items():
        q = a // 120
        e = (a - 120 * q, b, g - 60 * q)
        out[e] = out.get(e, 0) + c
    return {e: c for e, c in out.items() if c}


def fold(f):                                        # a real cosine sum: normalise each term to β ≥ 0
    out = {}
    for (a, b, g), c in f.items():
        e = (a, b, g) if b > 0 or (b == 0 and a >= 0) else (-a, -b, -g)
        out[e] = out.get(e, 0) + c
    return {e: c for e, c in out.items() if c}


TN = (90, 0, 45)                                    # −i = ζ^{3N}
def sin2(x):                                        # 2 sin(xπ/N) = −i(ζ^{2x} − ζ^{−2x})
    return {lin(mulc(x, 2), TN): 1, lin(mulc(x, -2), TN): -1}


A_, C_ = (0, 0, 3), (18, 0, 7)                      # angles a = 3, c = 18m + 7 (units π/N)
U1, IDX1 = (72, 0, 40), (36, 0, 20)                 # 2·idx1, idx1 = 2(N − b): side 1's direction
V1 = sin2(C_)
V2 = add(V1, sh(sin2(A_), IDX1))
C1 = add(V1, sh(cj(V1), U1), -1)                    # R1(z) = ζ^{U1} z̄ + C1;  R0(z) = z̄;  R2(z) = ζ^{12} z̄
D0 = (78, 6, 55)                                    # launch slot
RUNS = [(4, -2, -1), (8, -2, 2), (8, 2, 6), (4, 2, 4), (0, 2, 1)]
FIRST = [0, 2, 0, 0, 0]
ODD = [1, 0, 0, 0, 1]


def frun(r, first, odd):
    """the linear map of a run of r alternating letters from `first` (R0R2 = rotation ζ^{−12})"""
    if odd:
        return (1, mulc(lin(r, (0, 0, -1)), -6)) if first == 0 else (1, lin(mulc(r, 6), (0, 0, 6)))
    return (0, mulc(r, -6)) if first == 0 else (0, mulc(r, 6))


def compose_lin(f, g):
    (e1, E1), (e2, E2) = f, g
    return (e1 ^ e2, lin(E1, neg(E2)) if e1 else lin(E1, E2))


def aff_compose(F, Lm):
    e, E, b = F
    e2, E2 = Lm
    return (e ^ e2, lin(E, neg(E2)) if e else lin(E, E2), b)


def aff_R1(F):
    e, E, b = F
    return (e ^ 1, lin(E, neg(U1)) if e else lin(E, U1), add(sh(cj(C1) if e else C1, E), b))


def aff_apply(F, u):
    e, E, b = F
    return add(sh(cj(u) if e else u, E), b)


def S_formal():
    """4ℓ = w + w̄, w = conj(ω0)(F(V2) − V2), F = A1 R1 A2 R1 A3 R1 A4 R1 A5"""
    F = (0, (0, 0, 0), {})
    for i in range(5):
        F = aff_compose(F, frun(RUNS[i], FIRST[i], ODD[i]))
        if i < 4:
            F = aff_R1(F)
    w = sh(add(aff_apply(F, V2), V2, -1), neg(D0))
    return canon(add(w, cj(w)))


def evf(f, m, k, scale=1):
    N = 30 * m + 15
    return sum(c * mp.cos(mp.pi * (e[0] * m + e[1] * k + e[2]) / (2 * N)) for e, c in f.items()) / scale


# ---------------- (1) the formula against the walk, m ≤ 8 ----------------
def check_walk(ms=range(2, 9)):
    sys.path.insert(0, 'probes')
    from poly_slots import Poly
    import s4_7_witnesses as W
    S = S_formal()
    rows, bad = [], 0
    for m in ms:
        N = 30 * m + 15
        abc = (3, 12 * m + 5, 18 * m + 7)
        P = Poly(list(abc), N, lengths=[mp.sin(mp.pi * x / N) for x in (abc[2], abc[0], abc[1])])
        prs, kinds = W.pairs(P, mp.mpf(10) ** -20)
        assert len(prs) == m and all(q['loop'] for q in prs), (m, len(prs))
        for k, q in enumerate(prs):
            word = ''.join(map(str, q['word']))
            runs = [len(x) for x in word.split('1')]
            want = [r[0] * m + r[1] * k + r[2] for r in RUNS]
            ok_word = runs == want and q['slot'] == (2, 78 * m + 55 + 6 * k) and all(
                x[0] == str(f) and all(x[i] != x[i + 1] for i in range(len(x) - 1)) for x, f in zip(word.split('1'), FIRST))
            ell = evf(S, m, k, 4)
            err = abs(ell - mp.mpf(q['length']))
            bad += (not ok_word) or err > mp.mpf(10) ** -25
            rows.append({'m': m, 'k': k, 'len': mp.nstr(ell, 25), 'walk_word_affine': ok_word, 'err': float(err)})
    return rows, bad


# ---------------- (2) legality conditions ----------------
def fval(p):                                        # 2 Im(conj ω0 (p − V2))
    w = sh(add(p, V2, -1), neg(D0))
    return canon(sh(add(w, cj(w), -1), TN))


def re2(p):                                         # 2 Re(conj ω0 (p − V2))
    w = sh(add(p, V2, -1), neg(D0))
    return canon(add(w, cj(w)))


def conditions(k0):
    conds = []
    G = (0, (0, 0, 0), {})
    for i in range(5):
        r, s, odd = RUNS[i], FIRST[i], ODD[i]
        c = fval(aff_apply(G, {}))

        def ray(t, par, G=G, s=s):
            w = s if par == 0 else 2 - s
            return fval(aff_apply(aff_compose(G, frun(t, s, par)), V1 if w == 0 else V2))
        rm1, rm2 = lin(r, (0, 0, -1)), lin(r, (0, 0, -2))
        single = i == 4 and k0
        conds.append((f'run{i+1} t=0', 'opp', ray((0, 0, 0), 0), c))
        if not single:
            conds.append((f'run{i+1} t=1', 'opp', ray((0, 0, 1), 1), c))
            conds.append((f'run{i+1} last-even', 'opp', ray(rm1 if odd else rm2, 0), c))
            conds.append((f'run{i+1} last-odd', 'opp', ray(rm2 if odd else rm1, 1), c))
        if i < 4:
            conds.append((f'run{i+1} t=r (exit side 1)', 'same', ray(r, 1 if odd else 0), c))
        G = aff_compose(G, frun(r, s, odd))
        if i < 4:
            G = aff_R1(G)
    conds.append(('end at V2', 'zero', fval(aff_apply(G, V2)), None))
    conds.append(('launch: V0 ahead', 'pos', re2({}), None))
    conds.append(('launch: V1 ahead', 'pos', re2(V1), None))
    return conds


Z = sp.symbols('z')
PHI240 = sp.cyclotomic_poly(240, Z)


def an_zero(f, n):
    """a_n(y) = Σ c δ^n cos(θ + βy + nπ/2) ≡ 0 ⟸ (⟺ for β ≠ 0; only ⟸ is used — audit W15) per β the amplitude Σ c δ^n e^{i(θ + nπ/2)} = 0 in ℚ(ζ_240), ζ_240 = e^{iπ/120}"""
    groups = {}
    for (a, b, g), c in fold(f).items():
        groups.setdefault(b, []).append((c, a, Fraction(g) - Fraction(a, 2)))
    return all(sp.rem(sp.expand(sum(sp.Rational(d.numerator, d.denominator) ** n * c * Z ** ((2 * a + 60 * n) % 240)
                                    for c, a, d in ts)), PHI240, Z) == 0 for ts in groups.values())


def identity_zero(f):
    """f ≡ 0 as a function of (x, y): every (β, δ) group's amplitude vanishes"""
    grp = {}
    for (a, b, g), c in fold(f).items():
        grp.setdefault((b, Fraction(g) - Fraction(a, 2)), []).append((c, a))
    return all(sp.rem(sp.expand(sum(c * Z ** ((2 * a) % 240) for c, a in ts)), PHI240, Z) == 0 for ts in grp.values()), len(grp)


def Sc(t):
    t2 = t * t
    assert t2.b < 16
    return 1 - t2 / 6 + iv.mpf([0, 1]) * (t2 * t2 / 120)


def Vs(t):
    t2 = t * t
    assert t2.b < 16
    return iv.mpf(0.5) - t2 / 24 + iv.mpf([0, 1]) * (t2 * t2 / 720)


def Ws(t):
    t2 = t * t
    assert t2.b < 16
    return iv.mpf(-1) / 6 + iv.mpf([0, 1]) * (t2 / 120)


def prep(f):
    return [(c, iv.pi * a / 60, b, iv.mpf(g) - iv.mpf(a) / 2) for (a, b, g), c in fold(f).items()]


def F_direct(T, X, Y):
    return sum(c * iv.cos(th + b * Y + d * X) for c, th, b, d in T)


def F_H(T, X, Y):                                   # f/x when a_0 ≡ 0
    return sum(c * (-iv.sin(th + b * Y) * d * Sc(d * X) - iv.cos(th + b * Y) * d * d * X * Vs(d * X)) for c, th, b, d in T)


def F_corner0(T, R, Y):                             # f/(x y) near (0,0), ρ = x/y, when a_0 ≡ 0 and a_1(0) = 0
    tot = iv.mpf(0)
    for c, th, b, d in T:
        st, ct = iv.sin(th), iv.cos(th)
        dsin = -st * b * b * Y * Vs(b * Y) + ct * b * Sc(b * Y)
        tot += c * (-d * (dsin * Sc(d * R * Y) + st * d * d * R * R * Y * Ws(d * R * Y)) - iv.cos(th + b * Y) * d * d * R * Vs(d * R * Y))
    return tot


def F_cornerTop(T, R, U):                           # f/u near (0, π/60), ρ = x/u, when a_0(π/60) = 0
    tot = iv.mpf(0)
    for c, th, b, d in T:
        phi = th + b * iv.pi / 60
        lam = d * R - b
        tot += c * (-iv.cos(phi) * lam * lam * U * Vs(lam * U) - iv.sin(phi) * lam * Sc(lam * U))
    return tot


def F_x2(T, X, Y):                                  # f/x² at y = 0 when a_0 ≡ 0 and a_1(0) = 0
    return sum(c * (-iv.sin(th) * d * d * d * X * Ws(d * X) - iv.cos(th) * d * d * Vs(d * X)) for c, th, b, d in T)


def certify(fun, box, inside, sign, minw=1e-7, maxn=200000):
    stack, n = [box], 0
    while stack:
        x0, x1, y0, y1 = stack.pop()
        if not inside(x0, x1, y0, y1):
            continue
        v = fun(iv.mpf([x0, x1]), iv.mpf([y0, y1]))
        n += 1
        assert n < maxn
        if (v.a > 0 if sign > 0 else v.b < 0):
            continue
        if x1 - x0 < minw and y1 - y0 < minw:
            return False, n
        xm, ym = (x0 + x1) / 2, (y0 + y1) / 2
        stack += ([(x0, xm, y0, ym), (xm, x1, y0, ym), (x0, xm, ym, y1), (xm, x1, ym, y1)] if y1 > y0
                  else [(x0, xm, y0, y1), (xm, x1, y0, y1)])
    return True, n


PI60 = mp.pi / 60 + mp.mpf(10) ** -40                # OUTWARD (audit G23: mp.pi/60 rounds down, missing u = π/60, i.e. k = 0)
XMAX = mp.pi / 150 + mp.mpf(10) ** -40             # OUTWARD (m = 2 sits exactly at x = π/150)
RHO_TOP = mp.mpf(2) / 3 + mp.mpf(10) ** -30         # OUTWARD (audit G22: k = m − 1 sits exactly at ρ = 2/3)
IN_R1 = lambda x0, x1, y0, y1: y1 >= x0 and y0 <= PI60 - 1.5 * x0
ANY = lambda *a: True


def exact_zero_at(f, y_units):                      # a_0(y) = 0 at y = y_units·π/60, exactly
    e = sum(c * (Z ** ((2 * (a + b * y_units)) % 240) + Z ** ((-2 * (a + b * y_units)) % 240)) for (a, b, g), c in fold(f).items())
    return sp.rem(sp.expand(e), PHI240, Z) == 0


def a1_zero_at0(f):
    e = sum(c * sp.Rational(Fraction(g) - Fraction(a, 2)) * (Z ** ((2 * a) % 240) - Z ** ((-2 * a) % 240))
            for (a, b, g), c in fold(f).items())
    return sp.rem(sp.expand(e), PHI240, Z) == 0


def cert_sign(g, sg, k0):
    """sg·g > 0 on the k ≥ 1 region (k0 False) or on y = 0 (k0 True); returns (ok, method, boxes)"""
    T = prep(g)
    z0 = an_zero(g, 0)
    base = (lambda X, Y: F_H(T, X, Y)) if z0 else (lambda X, Y: F_direct(T, X, Y))
    if k0:
        if z0 and a1_zero_at0(g):
            ok, n = certify(lambda X, Y: F_x2(T, X, Y), (mp.mpf(0), XMAX, mp.mpf(0), mp.mpf(0)), ANY, sg)
            return ok, 'f/x² at y=0', n
        ok, n = certify(base, (mp.mpf(0), XMAX, mp.mpf(0), mp.mpf(0)), ANY, sg)
        return ok, 'f/x' if z0 else 'f', n
    Y0 = mp.pi / 600
    if z0 and a1_zero_at0(g):
        ok1, n1 = certify(lambda R, Y: F_corner0(T, R, Y), (mp.mpf(0), mp.mpf(1), mp.mpf(0), Y0), ANY, sg)
        ok2, n2 = certify(base, (mp.mpf(0), XMAX, Y0, PI60), IN_R1, sg)
        return ok1 and ok2, 'corner (0,0): f/(xy) in ρ = x/y ⊕ f/x', n1 + n2
    if not z0 and exact_zero_at(g, 1):
        ok1, n1 = certify(lambda R, U: F_cornerTop(T, R, U), (mp.mpf(0), RHO_TOP, mp.mpf(0), Y0), ANY, sg)
        ok2, n2 = certify(base, (mp.mpf(0), XMAX, mp.mpf(0), PI60 - Y0), IN_R1, sg)
        return ok1 and ok2, 'corner (0,π/60): f/u in ρ = x/u ⊕ f', n1 + n2
    ok, n = certify(base, (mp.mpf(0), XMAX, mp.mpf(0), PI60), IN_R1, sg)
    return ok, 'f/x' if z0 else 'f', n


def span_ok():
    """audit G7: the arc step needs each run's angular span (r − 1)·3π/N < π, i.e. 3(r − 1) < N = 30m + 15 for all
    0 ≤ k ≤ m − 1, m ≥ 2. r is affine in (m, k): check at both k-ends as an affine inequality in m (slope ⊕ m = 2, 3)."""
    for (a, b, g) in RUNS:
        for kk in ((0, 0), (1, -1)):                    # k = 0 and k = m − 1 as (coef of m, const)
            ra, rc = a + b * kk[0], g + b * kk[1]        # r = ra·m + rc
            assert 3 * ra <= 30 and all(3 * (ra * m + rc - 1) < 30 * m + 15 for m in (2, 3)), ('span', (a, b, g))
    return True


def legality():
    assert span_ok()
    out, allok = [], True
    for k0 in (False, True):
        for name, kind, f, c in conditions(k0):
            if kind == 'zero':
                ok, ng = identity_zero(f)
                out.append({'k0': k0, 'cond': name, 'ok': ok, 'method': f'identity: {ng} (β,δ) groups vanish in ℚ(ζ_240)'})
            elif kind == 'pos':
                ok, meth, n = cert_sign(f, 1, k0)
                out.append({'k0': k0, 'cond': name, 'ok': ok, 'method': meth, 'boxes': n})
            else:
                cs = 1 if evf(c, 50, 0 if k0 else 25) > 0 else -1
                okc, mc, nc = cert_sign(c, cs, k0)
                okf, mf, nf = cert_sign(f, -cs if kind == 'opp' else cs, k0)
                ok = okc and okf
                out.append({'k0': k0, 'cond': name, 'ok': ok, 'centre_sign': cs, 'method_c': mc, 'method_f': mf,
                            'boxes': nc + nf})
            allok &= ok
            print(('k=0 ' if k0 else 'k≥1 ') + f"{name:28s} {'OK' if ok else 'FAIL'}", flush=True)
    return out, allok


# ---------------- (3) the gap ----------------
def gap_cert(rho_max='0.4'):
    S = S_formal()
    terms = [(c, b // 6, a + 6 * (b // 6), Fraction(g) - Fraction(a, 2)) for (a, b, g), c in fold(S).items()]
    assert all(abs(b) == 6 for (a, b, g) in fold(S)), 'sin(βx/2) = σ sin 3x needs |β| = 6 (audit W17)'
    z = sp.symbols('z')
    ident = sp.rem(sp.expand(sum(c * s * (z ** (n % 120) - z ** ((-n) % 120)) for c, s, n, d in terms)),
                   sp.cyclotomic_poly(120, z), z) == 0
    SPH = [(c * s, iv.sin(iv.pi * n / 60), iv.cos(iv.pi * n / 60), s, iv.mpf(d.numerator) / d.denominator) for c, s, n, d in terms]

    def Hbox(R, U):
        tot = iv.mpf(0)
        for cs, sp_, cp, s, d in SPH:
            lam = (d + 3 * s) * R - 6 * s
            tot += cs * (-sp_ * lam * lam * U * Vs(lam * U) + cp * lam * Sc(lam * U))
        return tot
    ok, n = certify(Hbox, (mp.mpf(0), mp.mpf(rho_max), mp.mpf(0), PI60), ANY, -1, minw=1e-5)
    sup = Hbox(iv.mpf([0, mp.mpf(rho_max)]), iv.mpf([0, PI60])).b if ok else None
    return ident, ok, n, (float(sup) if sup is not None else None), terms


def main():
    res = {}
    S = S_formal()
    res['formula'] = {'nterms': len(S), 'terms': [[list(e), c] for e, c in sorted(S.items())],
                      'reading': '4ℓ(m,k) = Σ c·cos(π(αm + βk + γ)/(60m + 30))'}
    rows, bad = check_walk()
    res['walk_check'] = {'m': '2..8', 'rows': rows, 'mismatches': bad}
    print('walk check (word affine ⊕ exact length = walk length, m ≤ 8): mismatches', bad, flush=True)
    leg, allok = legality()
    res['legality'] = {'conditions': leg, 'all_certified': allok}
    print('LEGALITY all certified:', allok, flush=True)
    ident, ok, n, sup, terms = gap_cert(mp.mpf('0.4') + mp.mpf(10) ** -30)   # OUTWARD: k = m − 2 sits at ρ = 0.4
    ctl = gap_cert('0.8')
    res['gap'] = {'identity_sum_csin_phi_zero': ident, 'H_negative_certified': ok, 'boxes': n, 'sup_H': sup,
                  'control_rho_0.8_certifies': ctl[1],
                  'lower_bound': 'G ≥ ½ sin(3x)·2.5x·|sup H| = 57.375·x·sin 3x (sup H ≤ −45.9)', 'asymptote': 'min over k (k = m − 2): 0.6075/m²; small k: order 1/m'}
    print('GAP: identity', ident, '| H < 0 certified', ok, n, 'boxes, sup', sup, '| control (ρ ≤ 0.8) certifies:', ctl[1], flush=True)
    res['gaps_m_le_12'] = [[m, [mp.nstr(evf(S, m, k + 1, 4) - evf(S, m, k, 4), 8) for k in range(m - 1)]] for m in range(2, 13)]
    json.dump(res, open(OUT, 'w'), indent=1)
    assert bad == 0 and allok and ident and ok and not ctl[1]
    print('ALL OK →', OUT)


if __name__ == '__main__':
    main()
