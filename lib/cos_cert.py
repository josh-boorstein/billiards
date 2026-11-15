"""cos_cert — RIGOROUS SIGN CERTIFICATES for real cosine sums f(x, y) = Σ c·cos(πθ + βy + δx) (c integer, θ, β, δ rational) on
regions touching x = 0 — the interval half that `word_family` leaves to callers (s948; generalises the hand-written charts of
`probes/s938_q5_loops.py`).

COORDINATES. In a one-parameter triangle family with N = N_m·m + N_0, a word-family condition Σ c cos(π(αm + βt + γ)/2N) is
f(x, y) with x = π/2N, y = t·x, θ = α/(2N_m), δ = γ − α·N_0/N_m (`from_formal`). m → ∞ is x → 0, so a region "every m ≥ M0, t
between two affine bounds" is a polygon {0 < x ≤ X0, L(x) ≤ y ≤ U(x)} with L, U lines y = Y0 + s·x.

TAYLOR CHARTS. Every chart writes f = Σ c cos(φ(p) + λ(p)·v), v ∈ [0, V] the small variable, p a bounded parameter, φ, λ affine:
 - x-chart  (p = y, v = x):          φ = πθ + βy, λ = δ;
 - s-chart  at a corner (0, y_c):    y = y_c + σ s x  (p = s, v = x),  φ = πθ + βy_c, λ = δ + σβs;
 - ρ-chart  at a corner (0, y_c):    x = ρu, y = y_c + σu  (p = ρ, v = u), φ = πθ + βy_c, λ = σβ + δρ.
f/v^k = Σ_{j=k}^{K−1} v^{j−k} C_j(p) + R,  C_j(p) = Σ c λ^j cos(φ + jπ/2)/j!,  |R| ≤ Σ|c||λ|^K v^{K−k}/K!  (Lagrange; |cos⁽ᴷ⁾| ≤ 1),
where k is the EXACT vanishing order: C_j ≡ 0 for j < k is decided in ℚ(ζ_L) (`xorder`, `corner_order`), never numerically.
A sign certificate is bisection of the (p, v) box until every sub-box's enclosure has the sign (`certify`).

REGIONS (`certify_region`). Bounds are `Line(Y0, s, L0)`: y = Y0 + s·x, with L0 = Y0/π a Fraction when the corner (0, Y0) is
rational (then it may be a ZERO of the leading coefficient, decided exactly → s- and ρ-charts cover {σ(y − y_c) ∈ [σ s x, U0]},
the x-chart the rest), or L0 = None (a numerical bound, e.g. a band line off a transcendental boundary — no chart there; the
leading coefficient must be nonzero). Extra half-planes h = (A, b, c): πA + b·y + c·x ≥ 0 (a `needs r ≥ 2` constraint); used
CONSERVATIVELY — a box is skipped only when it lies entirely outside (an over-covered box can only make a certificate fail,
never pass wrongly). In corner charts a half-plane through the corner restricts s; one strictly positive at the corner is
ignored (over-coverage); one strictly negative removes the chart (the condition is not required there).
A THIN region (both bounds through the same rational corner) is a pure s-chart, p = s ∈ [s_L, s_U] (s_L = s_U: a single-t family).

API: `from_formal(f, Nm, N0)` → terms [(c, θ, β, δ)] · `xorder(T)` · `certify_region(T, sign, lo, hi, X0, halfplanes=())`
→ (ok, info) · `evalf(T, x, y)` (mpmath). PRECISION: `iv.dps` 30; boxes stop at width `minw`.
"""
from fractions import Fraction as Fr
from math import comb, factorial
import mpmath as mp
from mpmath import iv
import sympy as sp

iv.dps = 30
_Z = sp.symbols('z')


def from_formal(f, Nm, N0):
    """word_family formal sum {(α, β, γ): c} → [(c, θ, β, δ)] (Fractions), f = Σ c cos(πθ + βy + δx)"""
    out = {}
    for (a, b, g), c in f.items():
        a, b, g = Fr(a), Fr(b), Fr(g)
        key = (a / (2 * Nm), b, g - a * Fr(N0) / Nm)
        if key[1] < 0 or (key[1] == 0 and key[2] < 0) or (key[1] == 0 and key[2] == 0 and key[0] < 0):
            key = (-key[0], -key[1], -key[2])
        key = (key[0] % 2, key[1], key[2])
        out[key] = out.get(key, 0) + c
    return [(c, th, be, de) for (th, be, de), c in sorted(out.items()) if c]


# ---------------- exact tests in ℚ(ζ_L) ----------------
def _poly(pairs, real):
    """Σ w e^{iπr} (real=False) or Σ w cos(πr) (real=True) as a polynomial in z = e^{iπ/D}, D = lcm of denominators"""
    D = 1
    for _, r in pairs:
        D = sp.ilcm(D, Fr(r).denominator)
    L = 2 * D
    e = 0
    for w, r in pairs:
        k = int((Fr(r) * D) % L)
        wq = sp.Rational(Fr(w).numerator, Fr(w).denominator)
        e += wq * (_Z ** k + _Z ** ((-k) % L)) if real else wq * _Z ** k
    return e, L


def is_zero(pairs, real):
    if not pairs:
        return True
    e, L = _poly(pairs, real)
    return sp.rem(sp.expand(e), sp.cyclotomic_poly(L, _Z), _Z) == 0


def xcoef_zero(T, k):
    """C_k(y) = Σ c δ^k cos(πθ + βy + kπ/2)/k! ≡ 0 in y (exact): per |β| group the complex amplitude vanishes (β = 0: real part)"""
    grp = {}
    for c, th, be, de in T:
        w = c * Fr(de) ** k
        if w == 0:
            continue
        r = th + Fr(k, 2)
        if be < 0:
            be, r = -be, -r
        grp.setdefault(be, []).append((w, r))
    return all(is_zero(ps, real=(be == 0)) for be, ps in grp.items())


def xorder(T, kmax=6):
    for k in range(kmax + 1):
        if not xcoef_zero(T, k):
            return k
    return None                                          # identically zero to order kmax: an identity


def coef_at_zero(T, k, L0):
    """C_k(πL0) = 0 exactly?"""
    return is_zero([(c * Fr(de) ** k, th + be * L0 + Fr(k, 2)) for c, th, be, de in T], real=True)


def corner_order(T, L0, l0f, l1f, kmax=8):
    """vanishing order in v of f = Σ c cos(π(θ + βL0) + λ v), λ = l0f(term) + l1f(term)·p (Fractions): C_j(p) ≡ 0 iff every
    p-power coefficient Σ c binom(j,i) l1^i l0^{j−i} cos(π(θ + βL0 + j/2)) vanishes (exact)"""
    for j in range(kmax + 1):
        ok = True
        for i in range(j + 1):
            ps = [(c * comb(j, i) * Fr(l1f(t)) ** i * Fr(l0f(t)) ** (j - i), t[1] + t[2] * L0 + Fr(j, 2))
                  for t in T for c in (t[0],)]
            ps = [(w, r) for w, r in ps if w != 0]
            if not is_zero(ps, real=True):
                ok = False
                break
        if not ok:
            return j
    return None


# ---------------- enclosures ----------------
def _ivf(q):
    q = Fr(q)
    return iv.mpf(q.numerator) / q.denominator


def evalf(T, x, y):
    return sum(c * mp.cos(mp.pi * mp.mpf(th.numerator) / th.denominator + mp.mpf(be.numerator) / be.denominator * y
                          + mp.mpf(de.numerator) / de.denominator * x) for c, th, be, de in T)


def _taylor(terms, k, K, P, V):
    """terms [(c, φ(P) interval, λ(P) interval)]; enclosure of f/v^k on V = [0, v1] (or [v0, v1], v0 ≥ 0)"""
    tot = iv.mpf(0)
    for j in range(k, K):
        Cj = iv.mpf(0)
        for c, ph, lam in terms:
            Cj += c * lam ** j * iv.cos(ph + j * iv.pi / 2)
        tot += Cj / factorial(j) * (V ** (j - k) if j > k else 1)
    rem = sum(abs(c) * iv.mpf(max(abs(lam.a), abs(lam.b))) ** K for c, ph, lam in terms) / factorial(K)
    rem = rem * (iv.mpf(V.b) ** (K - k))
    return tot + iv.mpf([-1, 1]) * rem


def certify(fun, box, inside, sign, minw=1e-9, maxn=400000):
    """bisect box = (p0, p1, v0, v1) until sign·fun(P, V) > 0 on every box that `inside` keeps"""
    stack, n = [box], 0
    while stack:
        p0, p1, v0, v1 = stack.pop()
        if not inside(p0, p1, v0, v1):
            continue
        val = fun(iv.mpf([p0, p1]), iv.mpf([v0, v1]))
        n += 1
        if n >= maxn:
            return False, n, (p0, p1, v0, v1)
        if (val.a > 0) if sign > 0 else (val.b < 0):
            continue
        wp, wv = p1 - p0, v1 - v0
        if wp < minw and wv < minw:
            return False, n, (p0, p1, v0, v1)
        pm, vm = (p0 + p1) / 2, (v0 + v1) / 2
        if wp > 0 and (wp >= wv or wv == 0):
            stack += [(p0, pm, v0, v1), (pm, p1, v0, v1)]
        else:
            stack += [(p0, p1, v0, vm), (p0, p1, vm, v1)]
    return True, n, None


class Line:
    """y = Y0 + s·x; L0 = Y0/π as a Fraction when exact (a chartable corner), else None"""

    def __init__(self, Y0, s, L0=None):
        self.Y0, self.s, self.L0 = mp.mpf(Y0), s, L0

    @staticmethod
    def exact(L0, s):
        return Line(mp.pi * mp.mpf(Fr(L0).numerator) / Fr(L0).denominator, Fr(s), Fr(L0))

    def at(self, x):
        return self.Y0 + mp.mpf(Fr(self.s).numerator) / Fr(self.s).denominator * x if isinstance(self.s, Fr) else self.Y0 + self.s * x


def _hp_val(h, x, y):
    A, b, c = h
    return mp.pi * mp.mpf(Fr(A).numerator) / Fr(A).denominator + mp.mpf(Fr(b).numerator) / Fr(b).denominator * y \
        + mp.mpf(Fr(c).numerator) / Fr(c).denominator * x


def _box_outside(x0, x1, y0, y1, lo, hi, hps, eps=mp.mpf(10) ** -25):
    cs = [(x, y) for x in (x0, x1) for y in (y0, y1)]
    if all(y < lo.at(x) - eps for x, y in cs) or all(y > hi.at(x) + eps for x, y in cs):
        return True
    return any(all(_hp_val(h, x, y) < -eps for x, y in cs) for h in hps)


def _chart_terms_x(T, Y, k_unused=None):
    return [(c, iv.pi * _ivf(th) + _ivf(be) * Y, _ivf(de)) for c, th, be, de in T]


def _corner_terms(T, L0, lam):                           # lam(term) → interval λ
    return [(t[0], iv.pi * (_ivf(t[1]) + _ivf(t[2]) * _ivf(L0)), lam(t)) for t in T]


def rho_chart(T, L0, sg, k, extra=6):
    """ρ-chart at the corner (0, πL0), x = ρu, y = πL0 + sg·u, for f of x-order k: f = Σ_{a,b} F_ab x^a (sg·u)^b … ⇒
    f/x^k = Σ_n u^{n−k} P_n(ρ), P_n(ρ) = Σ_{a=k}^{n} F_{a,n−a} ρ^{a−k}, F_ab = Σ c δ^a (sgβ)^b cos(φ + nπ/2)/(a! b!), φ = π(θ + βL0)
    (F_ab = 0 for a < k: the x-order). n0 = the first n with P_n ≢ 0 (exact); returns (n0, enclosure(P, V) of f/(x^k u^{n0−k})).
    Tail: |P_n| ≤ Σ|c| |δ|^k/k! · Λ^{n−k}/(n−k)!, Λ = |δ|ρ + |β|, so Σ_{n≥K} u^{n−n0}|P_n| ≤ Σ|c||δ|^k/k! Λ^{K−k} e^{Λu} u^{K−n0}/(K−k)!."""
    def F_exact(a, b):
        n = a + b
        return [(Fr(c) * Fr(de) ** a * (sg * Fr(be)) ** b / (factorial(a) * factorial(b)), th + be * L0 + Fr(n, 2))
                for c, th, be, de in T if not (a and de == 0) and not (b and be == 0)]
    n0 = None
    for n in range(k, k + 12):
        if any(not is_zero(F_exact(a, n - a), real=True) for a in range(k, n + 1)):
            n0 = n
            break
    if n0 is None:
        return None, None
    K = n0 + extra

    def F_iv(a, b):
        tot = iv.mpf(0)
        for w, r in F_exact(a, b):
            tot += _ivf(w) * iv.cos(iv.pi * _ivf(r))
        return tot
    coef = {n: [F_iv(a, n - a) for a in range(k, n + 1)] for n in range(n0, K)}
    split = is_zero(F_exact(k, n0 - k), real=True)       # P_{n0}(0) = 0: f/(x^k u^{n0−k}) = ρ·A(ρ) + u·B(ρ, u) near ρ = 0
    tails = [(abs(c) * abs(_ivf(de)) ** k / factorial(k), abs(_ivf(de)), abs(_ivf(be))) for c, th, be, de in T]

    def tail(P, V, drop):
        tl = iv.mpf(0)
        for w, d, b in tails:
            Lam = d * iv.mpf(P.b) + b
            tl += w * Lam ** (K - k) * iv.exp(Lam * iv.mpf(V.b)) * iv.mpf(V.b) ** (K - n0 - drop) / factorial(K - k)
        return iv.mpf([-1, 1]) * tl

    def enc(P, V):
        if split and P.a == 0:                            # hull of A and B: certified iff both carry the sign (ρ, u ≥ 0)
            A = sum((cf * P ** (i - 1) for i, cf in enumerate(coef[n0]) if i >= 1), iv.mpf(0))
            B = sum((sum((cf * P ** i for i, cf in enumerate(coef[n])), iv.mpf(0)) * (V ** (n - n0 - 1) if n > n0 + 1 else 1)
                     for n in range(n0 + 1, K)), iv.mpf(0)) + tail(P, V, 1)
            return iv.mpf([min(A.a, B.a), max(A.b, B.b)])
        tot = iv.mpf(0)
        for n in range(n0, K):
            Pn = sum((cf * P ** i for i, cf in enumerate(coef[n])), iv.mpf(0))
            tot += Pn * (V ** (n - n0) if n > n0 else 1)
        tl = iv.mpf(0)
        for w, d, b in tails:
            Lam = d * iv.mpf(P.b) + b
            tl += w * Lam ** (K - k) * iv.exp(Lam * iv.mpf(V.b)) * iv.mpf(V.b) ** (K - n0) / factorial(K - k)
        return tot + iv.mpf([-1, 1]) * tl
    return n0, enc


def certify_region(T, sign, lo, hi, X0, halfplanes=(), S=8, U0=mp.mpf('0.01'), extra=6, minw=1e-9, maxn=400000):
    """sign·f > 0 on {0 < x ≤ X0, lo(x) ≤ y ≤ hi(x)} ∩ halfplanes. Returns (ok, info)."""
    info = {'charts': []}
    X0 = mp.mpf(X0)
    # ---- thin region: one rational corner
    if lo.L0 is not None and hi.L0 is not None and lo.L0 == hi.L0:
        L0, sl, su = lo.L0, Fr(lo.s), Fr(hi.s)
        for h in halfplanes:                                    # restrict s by half-planes through the corner
            A, b, c = h
            at = Fr(A) + Fr(b) * L0
            if at < 0:                                          # violated near x = 0; removed only if violated on all x ≤ X0
                smax = max(abs(Fr(lo.s)), abs(Fr(hi.s)))
                if mp.pi * mp.mpf(at.numerator) / at.denominator + (abs(Fr(b)) * smax + abs(Fr(c))) * X0 < 0:
                    info['charts'].append('thin: removed by half-plane')
                    return True, info
                continue
            if at == 0:                                         # b·s + c ≥ 0 on the chart (y = y_c + s x)
                if b > 0:
                    sl = max(sl, Fr(-c) / b)
                elif b < 0:
                    su = min(su, Fr(-c) / b)
                elif c < 0:
                    return True, info
        if sl > su:
            return True, info
        if sl == su:                                            # a single line: the order ON it (λ constant)
            k = corner_order(T, L0, lambda t: t[3] + t[2] * sl, lambda t: 0)
        else:
            k = corner_order(T, L0, lambda t: t[3], lambda t: t[2])
        if k is None:
            return False, {'error': 'identically zero on the chart'}
        K = k + extra
        fun = lambda P, V: _taylor(_corner_terms(T, L0, lambda t: _ivf(t[3]) + _ivf(t[2]) * P), k, K, P, V)
        ok, n, bad = certify(fun, (mp.mpf(sl.numerator) / sl.denominator, mp.mpf(su.numerator) / su.denominator, mp.mpf(0), X0),
                             lambda *a: True, sign, minw=minw, maxn=maxn)
        info['charts'].append({'thin s-chart': [str(sl), str(su)], 'order': k, 'boxes': n, 'ok': ok, 'bad': bad})
        return ok, info
    # ---- wide region
    k = xorder(T)
    if k is None:
        return False, {'error': 'identically zero'}
    K = k + extra
    info['xorder'] = k
    excl = []                                                   # (y_lo, y_hi) bands the corner charts cover
    allok = True
    for bnd, sg in ((lo, 1), (hi, -1)):
        if bnd.L0 is None or not coef_at_zero(T, k, bnd.L0):
            continue
        L0, s0 = bnd.L0, Fr(bnd.s)
        # s-chart: y = y_c + sg·s·x, s ∈ [sg·s0, S]  (the region side of the bound)
        s_lo = sg * s0
        hp_ok = True
        s_hi = Fr(S)
        for h in halfplanes:
            A, b, c = h
            at = Fr(A) + Fr(b) * L0
            if at < 0:
                hp_ok = False
            elif at == 0 and b != 0:
                # b·(y_c + sg s x) + c x ≥ 0 ⟺ (b sg s + c) ≥ 0
                bound = Fr(-c) / (b * sg)
                if b * sg > 0:
                    s_lo = max(s_lo, bound)
                else:
                    s_hi = min(s_hi, bound)
        if not hp_ok:                                           # the condition is not required at this corner: no chart, no
            info['charts'].append({'corner': str(L0), 'half-plane negative': True})   # exclusion (x-chart boxes there are skipped
            continue                                                                  # by `inside` once wholly outside it)
        kc = corner_order(T, L0, lambda t: t[3], lambda t, sg=sg: sg * t[2])
        Kc = kc + extra
        if s_lo <= s_hi:
            fun_s = lambda P, V, L0=L0, sg=sg, kc=kc, Kc=Kc: _taylor(
                _corner_terms(T, L0, lambda t: _ivf(t[3]) + sg * _ivf(t[2]) * P), kc, Kc, P, V)
            ok1, n1, b1 = certify(fun_s, (mp.mpf(s_lo.numerator) / s_lo.denominator, mp.mpf(s_hi.numerator) / s_hi.denominator,
                                          mp.mpf(0), X0), lambda *a: True, sign, minw=minw, maxn=maxn)
        else:
            ok1, n1, b1 = True, 0, None
        kr, fun_r = rho_chart(T, L0, sg, k, extra)
        if fun_r is None:
            return False, {'error': 'rho chart identically zero', 'corner': str(L0)}
        ok2, n2, b2 = certify(fun_r, (mp.mpf(0), mp.mpf(1) / S, mp.mpf(0), U0), lambda *a: True, sign, minw=minw, maxn=maxn)
        info['charts'].append({'corner': str(L0), 'side': sg, 's': [str(s_lo), str(s_hi)], 'order_s': kc, 'order_rho': kr,
                               'ok': ok1 and ok2, 'boxes': n1 + n2, 'bad': b1 or b2})
        allok &= ok1 and ok2
        excl.append((bnd.Y0 - U0, bnd.Y0 + U0))
    ylo = min(lo.at(0), lo.at(X0)) - mp.mpf(10) ** -20
    yhi = max(hi.at(0), hi.at(X0)) + mp.mpf(10) ** -20

    def inside(y0, y1, x0, x1):
        if any(y0 >= a and y1 <= b for a, b in excl):
            return False
        return not _box_outside(x0, x1, y0, y1, lo, hi, halfplanes)
    fun_x = lambda Y, X: _taylor(_chart_terms_x(T, Y), k, K, Y, X)
    ok, n, bad = certify(fun_x, (ylo, yhi, mp.mpf(0), X0), inside, sign, minw=minw, maxn=maxn)
    info['charts'].append({'x-chart': True, 'order': k, 'boxes': n, 'ok': ok, 'bad': bad})
    return allok and ok, info
