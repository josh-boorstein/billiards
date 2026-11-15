"""word_family — formal (symbolic in (m, t)) billiard words in a one-parameter family of rational triangles, and their
legality conditions as finite cosine sums. The reusable core of the s938 pipeline (`s938-e`), generalised (s947).

SETTING. Triangle family with angles a_i(m)·π/N(m), each affine in m (`Fam`), frame of `s910_poly_slots.Poly`: V0 = 0,
side j runs V_j → V_{j+1} with direction index idx_j (units π/2N; idx_0 = 0, idx_j = idx_{j−1} + 2(N − a_j)), side lengths
2 sin(a_{j+2}π/N) (side j is opposite vertex j+2). ζ := ζ_{4N} = e^{iπ/2N}. A formal sum is {E: c}, E = (α, β, γ) a TRIPLE
of Fractions meaning ζ^{αm + βt + γ}, c an integer. Every point of the unfolding (vertex images) is a formal sum; every
condition is a real cosine sum f(m, t) = Σ c cos(π(αm + βt + γ)/2N).

WORD FAMILY. A launch from vertex v at direction d(m, t) (formal exponent), and a word = runs of alternating letters in
the two sides through V0 (sides 0 and 2), separated by single side-1 letters: `runs` = list of (length E, first letter).
A run's composite map depends on its length's PARITY, so a family is evaluated on a parity class (m mod 2, t mod 2)
(`parity` argument); within it every run length must have a fixed parity (asserted).

CONDITIONS (`conditions`). With g(p) := 2 Im(conj ζ^d (p − S)) (side of the line), in each run (star centre c = image
of V0): the extreme crossed rays of each parity lie OPPOSITE c (g·g(c) < 0); the ray after the last lies on c's side (the
line leaves through side 1); on an arc of < π the crossed rays of one parity form an arc (ρ cos(ψ − ψ_n) > dist is an arc
of half-width < π/2), so the extremes suffice — the span condition (r − 1)·a_0 < N is returned as an integer check.
Endings: 'perp' on the last letter (a ray: Re(conj ζ^d (ray end − c)) ≡ 0; or side 1 after an empty last run:
Re(conj ζ^d (G(V2) − G(V1))) ≡ 0), 'vertex' (the end vertex image ON the line, g ≡ 0). Each condition is
(name, kind, f, c) with kind ∈ {'opp', 'same', 'zero'}.

CHECKS: `evalf(f, m, t)` numerical value; `identity_zero(f)` exact: f ≡ 0 as a function of (x, y), x = π/2N, y = t·x
(every (β, δ)-group's amplitude vanishes in ℚ(ζ_L)). The interval certification of 'opp'/'same' is NOT here (callers).
"""
from fractions import Fraction as Fr
import mpmath as mp
import sympy as sp

E0 = (Fr(0), Fr(0), Fr(0))


def E(a=0, b=0, c=0):
    return (Fr(a), Fr(b), Fr(c))


def eadd(*es):
    return tuple(sum(x[i] for x in es) for i in range(3))


def eneg(e):
    return tuple(-x for x in e)


def emul(e, k):
    return tuple(Fr(k) * x for x in e)


def add(u, v, s=1):
    r = dict(u)
    for e, c in v.items():
        r[e] = r.get(e, 0) + s * c
    return {e: c for e, c in r.items() if c}


def sh(u, e):
    return {eadd(k, e): c for k, c in u.items()}


def cj(u):
    return {eneg(k): c for k, c in u.items()}


class Fam:
    """a = [(am, a0)] * 3 (angle numerators a_i = am·m + a0), N = am-sum (must equal Σ a_i)"""

    def __init__(self, a):
        self.a = [E(am, 0, a0) for am, a0 in a]
        self.N = eadd(*self.a)
        N = self.N
        self.idx = [E0]
        for j in (1, 2):
            self.idx.append(eadd(self.idx[-1], emul(eadd(N, eneg(self.a[j])), 2)))
        self.mI = eneg(N)                                   # −i = ζ^{−N}

        def two_sin(e):                                     # 2 sin(eπ/N) = −i(ζ^{2e} − ζ^{−2e})
            return {eadd(self.mI, emul(e, 2)): 1, eadd(self.mI, emul(e, -2)): -1}
        L = [two_sin(self.a[(j + 2) % 3]) for j in range(3)]
        self.V = [{}]
        self.V.append(sh(L[0], self.idx[0]))
        self.V.append(add(self.V[1], sh(L[1], self.idx[1])))
        self.R = []                                          # R_j(z) = ζ^{2 idx_j} z̄ + c_j
        for j in range(3):
            F = self.red(emul(self.idx[j], 2))
            Vj = self.V[j]
            self.R.append((F, add(Vj, sh(cj(Vj), F), -1)))

    def red(self, e):
        """subtract a multiple of 4N that kills the m-coefficient, when one exists (exponents are mod 4N)"""
        q = e[0] / (4 * self.N[0])
        if q.denominator == 1:
            return eadd(e, emul(self.N, -4 * q))
        return e

    def Nval(self, m):
        return int(self.N[0] * m + self.N[2])


# G(z) = ζ^E·(z̄ if e else z) + b, as (e, E, b)
ID = (0, E0, {})


def compose_refl(G, Rj):
    e, Ex, b = G
    F, c = Rj
    if e == 0:
        return (1, eadd(Ex, F), add(b, sh(c, Ex)))
    return (0, eadd(Ex, eneg(F)), add(b, sh(cj(c), Ex)))


def apply(G, u):
    e, Ex, b = G
    return add(sh(cj(u) if e else u, Ex), b)


def compose_run(fam, G, r, first, parity_r):
    """G ∘ (R_first R_other)^{k} (∘ R_first if r odd), r = 2k + parity_r; letters alternate in sides 0, 2 (both through V0 = 0,
    so R_0 R_2 is the rotation ζ^{2 idx_0 − 2 idx_2} and R_2 R_0 its inverse)"""
    other = 2 - first
    Fa, Fb = fam.R[first][0], fam.R[other][0]
    k = emul(eadd(r, E(0, 0, -parity_r)), Fr(1, 2))
    # rotation per pair: R_a R_b (z) = ζ^{Fa} conj(ζ^{Fb} z̄) = ζ^{Fa − Fb} z ; k pairs: ζ^{k (Fa − Fb)}
    Rexp = emul_e(k, fam.red(eadd(Fa, eneg(Fb))))
    e, Ex, b = G
    # G ∘ rot: G(ζ^{X} z) = ζ^{Ex}·conj^e(ζ^X z) + b
    G2 = (e, eadd(Ex, eneg(Rexp)) if e else eadd(Ex, Rexp), b)
    if parity_r:
        G2 = compose_refl(G2, fam.R[first])
    return G2


def emul_e(k, F):
    """exponent product of an affine k = (km, kt, k0) with a CONSTANT-in-t exponent F = (Fm, 0, F0) — must stay affine:
    requires km·Fm = 0 and kt·Fm = 0 (F constant) or k constant"""
    km, kt, k0 = k
    Fm, Ft, F0 = F
    assert Ft == 0
    if Fm != 0:
        assert km == 0 and kt == 0, ('non-affine exponent', k, F)
        return (k0 * Fm, Fr(0), k0 * F0)
    return (km * F0, kt * F0, k0 * F0)


def val(e, m, t):
    return e[0] * m + e[1] * t + e[2]


class Word:
    """vertex v, slot exponent d (E triple), runs [(len E, first letter)], end in {'perp', 'vertex'}"""

    def __init__(self, fam, v, d, runs, end):
        self.fam, self.v, self.d, self.runs, self.end = fam, v, d, runs, end

    def g(self, p):                                     # 2 Im(conj ζ^d (p − S))
        w = sh(add(p, self.fam.V[self.v], -1), eneg(self.d))
        return add(sh(w, self.fam.mI), sh(cj(w), self.fam.mI), -1)

    def fwd(self, p):                                   # 2 Re(conj ζ^d p) (a direction test, p a difference)
        w = sh(p, eneg(self.d))
        return add(w, cj(w))

    def conditions(self, parity, tail_checks=True):
        """parity = (m mod 2, t mod 2); returns (conds, span_checks)"""
        fam, conds, spans = self.fam, [], []
        G = ID
        nr = len(self.runs)
        for i, (r, first) in enumerate(self.runs):
            pr = int(val(r, parity[0], parity[1])) % 2
            # parity fixed on the class: r(m + 2, t) − r(m, t) and r(m, t + 2) − r(m, t) must be even
            assert r[0].denominator == 1 and r[1].denominator == 1, r
            c = self.g(apply(G, {}))
            last = i == nr - 1

            def ray(u_e, par, G=G, first=first):        # the u-th crossed ray (u = 2j + par), as the vertex image it ends at
                letter = first if par == 0 else 2 - first
                Gu = compose_run(fam, G, u_e, first, par)
                return apply(Gu, fam.V[1] if letter == 0 else fam.V[2])
            rlen_is_zero = (r == E0)
            if not rlen_is_zero:
                spans.append(r)
                conds.append((f'run{i} u=0', 'opp', self.g(ray(E0, 0)), c))
                conds.append((f'run{i} u=1', 'opp', self.g(ray(E(0, 0, 1), 1)), c, ('needs r>=2', r)))
                lastE = eadd(r, E(0, 0, -1))                # u = r − 1, parity 1 − pr ... (u parity = (r−1) mod 2)
                prevE = eadd(r, E(0, 0, -2))
                conds.append((f'run{i} u=r-1', 'opp', self.g(ray(lastE, 1 - pr)), c))
                conds.append((f'run{i} u=r-2', 'opp', self.g(ray(prevE, pr)), c, ('needs r>=2', r)))
            if not last:
                if not (rlen_is_zero and i == 0):          # an empty first run (launch from V0): slot validity, no ray
                    conds.append((f'run{i} u=r exit', 'same', self.g(ray(r, pr)), c))
                G = compose_run(fam, G, r, first, pr)
                G = compose_refl(G, fam.R[1])
            else:
                if self.end == 'perp':
                    if rlen_is_zero:                           # perpendicular on side 1 of the current copy
                        conds.append(('perp side 1', 'zero', self.fwd(add(apply(G, fam.V[2]), apply(G, fam.V[1]), -1)), None))
                    else:
                        lastE = eadd(r, E(0, 0, -1))
                        conds.append(('perp last ray', 'zero', self.fwd(add(ray(lastE, 1 - pr), apply(G, {}), -1)), None))
                else:                                          # vertex endings are certified elsewhere (e.g. [TRI-BAND-QM])
                    conds.append(('end at vertex', 'vertex', None, None))
        return conds, spans


def evalf(f, fam, m, t):
    N = fam.Nval(m)
    return sum(c * mp.cos(mp.pi * val(e, m, t) / (2 * N)) for e, c in f.items())


def identity_zero(f, fam, L=None):
    """f ≡ 0 as a function of (m, t) over the reals: write π(αm + βt + γ)/2N = θ + β y + δ x with x = π/2N, y = t x,
    θ = απ/(2 N_m), δ = γ − α N_0/N_m; group by (β, δ) — f ≡ 0 iff every group's Σ c e^{iθ} has zero real... we test the
    stronger: Σ c ζ_L^{...} (as a cosine sum, folded to β ≥ 0) vanishes for every group."""
    Nm, N0 = fam.N[0], fam.N[2]
    grp = {}
    for (a, b, g), c in f.items():
        de = Fr(g) - Fr(a) * N0 / Nm
        if b < 0 or (b == 0 and de < 0):            # fold cos(πθ + βy + δx) = cos(−πθ − βy − δx) on the KEY (β, δ)
            a, b, g, de = -a, -b, -g, -de
        th = Fr(a) / (2 * Nm)                       # θ/π
        grp.setdefault((b, de), []).append((c, th))
    z = sp.symbols('z')
    ok = True
    def key_is_const(ts, grp):
        return any(k == (0, 0) and v is ts for k, v in grp.items())
    for ts in grp.values():
        den = 1
        for _, th in ts:
            den = sp.ilcm(den, th.denominator)
        Lm = 2 * den
        expr = sum(c * (z ** int((th * den) % Lm)) for c, th in ts)   # e^{iπ th} = ζ_{2den}^{th·den}
        # the group is Σ c cos(πth + φ) with common φ = βy + δx; vanishing ⟺ Σ c e^{iπ th} = 0 (φ non-constant);
        # the constant group (β, δ) = (0, 0) needs only its real part: Σ c (z^k + z^{−k}) = 0
        if ts and key_is_const(ts, grp):
            expr = sum(c * (z ** int((th * den) % Lm) + z ** int((-th * den) % Lm)) for c, th in ts)
        if sp.rem(sp.expand(expr), sp.cyclotomic_poly(Lm, z), z) != 0:
            ok = False
    return ok, len(grp)
