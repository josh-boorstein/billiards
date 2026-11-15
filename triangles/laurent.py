"""[Deposit] helper: exact Laurent polynomials in y = e^{iB} over Q(zeta_2n). Paper: "Complete periodicity of perpendicular billiards in rational triangles of prime denominator". Research name: probes/s944_laurent.py, copied verbatim except imports (deposit names), paths (data/ and figs/ under triangles/) and this paragraph. Run through run_all.py; the usage lines below are the research repo's. Bracketed IDs and section marks (§sNNN) refer to the project's internal ledgers.

s944 — exact holonomy of a triangle word as a Laurent polynomial in (x, y) = (e^{iA}, e^{iB}), integer coefficients.

Library for `s944_relation_torsion.py` (session-lineage, not engine). Convention = `s922_band_shooter.develop`:
V0 = 0 (angle A), V1 = sin C (angle B), V2 = sin B e^{iA} (angle C); letter 0 = side V0V1, 1 = V1V2, 2 = V2V0.
Everything is scaled by 2i so that 2i sin A = x - 1/x etc. have integer coefficients; reflections are carried on the
unit side DIRECTIONS (signed monomials), so no conjugation of a position is ever needed:
    V_k' = V_i + len(i,k) * D'(i,k),   D'(i,k) = D(i,j)^2 * conj D(i,k).
`hol(w)` = 2i * (V2_end - V2_start) as {(p, q): int}. Connection length ell = 2|H| (the shooter's `2|E - S|`), so
ell^2 = 4 |H|^2 = hol * hol~ with hol~(x, y) = hol(1/x, 1/y)  (|2i H|^2 = 4|H|^2; hol~ = conj on |x| = |y| = 1).
"""
from collections import defaultdict

SIDES = [(0, 1), (1, 2), (2, 0)]


def _mono_mul(m1, m2):  # signed monomial (s, p, q)
    return (m1[0] * m2[0], m1[1] + m2[1], m1[2] + m2[2])


def _mono_conj(m):
    return (m[0], -m[1], -m[2])


# scaled side lengths 2i*sin: side {1,2} opposite V0 -> A; {0,2} opposite V1 -> B; {0,1} opposite V2 -> C = pi - A - B
LEN = {frozenset((1, 2)): {(1, 0): 1, (-1, 0): -1},
       frozenset((0, 2)): {(0, 1): 1, (0, -1): -1},
       frozenset((0, 1)): {(1, 1): 1, (-1, -1): -1}}   # 2i sin(A+B) = xy - 1/(xy)


def padd(P, Q, c=1):
    R = defaultdict(int, P)
    for k, v in Q.items():
        R[k] += c * v
    return {k: v for k, v in R.items() if v}


def pmul_mono(P, m):
    s, p, q = m
    return {(a + p, b + q): s * v for (a, b), v in P.items()}


def pmul(P, Q):
    R = defaultdict(int)
    for (a, b), u in P.items():
        for (c, d), v in Q.items():
            R[(a + c, b + d)] += u * v
    return {k: v for k, v in R.items() if v}


def ptilde(P):
    return {(-a, -b): v for (a, b), v in P.items()}


def hol(word):
    # D[(i,j)] unit direction V_i -> V_j as a signed monomial
    D = {(0, 1): (1, 0, 0), (1, 2): (-1, 0, -1), (0, 2): (1, 1, 0)}
    for (i, j) in list(D):
        D[(j, i)] = (-D[(i, j)][0], D[(i, j)][1], D[(i, j)][2])
    V = [{}, dict(LEN[frozenset((0, 1))]), pmul_mono(LEN[frozenset((0, 2))], (1, 1, 0))]
    S = V[2]
    for ch in word:
        i, j = SIDES[int(ch)]
        k = 3 - i - j
        u2 = _mono_mul(D[(i, j)], D[(i, j)])
        nik = _mono_mul(u2, _mono_conj(D[(i, k)]))
        njk = _mono_mul(u2, _mono_conj(D[(j, k)]))
        V = list(V)
        V[k] = padd(V[i], pmul_mono(LEN[frozenset((i, k))], nik))
        D[(i, k)] = nik; D[(k, i)] = (-nik[0], nik[1], nik[2])
        D[(j, k)] = njk; D[(k, j)] = (-njk[0], njk[1], njk[2])
    return padd(V[2], S, -1)


def len2(word):
    """hol * hol~ = |2iH|^2 = 4|H|^2 = ell^2 (ell = 2|H|): integer Laurent poly in (x, y)."""
    h = hol(word)
    return pmul(h, ptilde(h))


def cyclo(m):
    """Integer coefficients (low->high) of the m-th cyclotomic polynomial."""
    import sympy as sp
    t = sp.symbols('t')
    return [int(c) for c in reversed(sp.Poly(sp.cyclotomic_poly(m, t), t).all_coeffs())]


def _reduce(coeffs, phi):
    c = list(coeffs); d = len(phi) - 1
    for e in range(len(c) - 1, d - 1, -1):
        f = c[e]
        if f:
            for r in range(d + 1):
                c[e - d + r] -= f * phi[r]
    return tuple(c[:d])


def on_line(P, n, phi=None):
    """Restrict to A = pi/n: x = zeta_{2n}. Returns {q: coeff tuple in Q(zeta_2n) basis} with zeros dropped (exact)."""
    phi = phi or cyclo(2 * n)
    byq = defaultdict(dict)
    for (p, q), v in P.items():
        e = p % (2 * n)
        byq[q][e] = byq[q].get(e, 0) + v
    out = {}
    for q, dct in byq.items():
        c = [0] * (2 * n)
        for e, v in dct.items():
            c[e] += v
        r = _reduce(c, phi)
        if any(r):
            out[q] = r
    return out


def span(dq):
    return (max(dq) - min(dq)) if dq else -1


if __name__ == '__main__':   # validation against the float/mp developer
    import sys, random, mpmath as mp
    import _setup  # noqa: F401  (deposit: lib/ on the path, cwd = triangles/)
    import band_shooter as SH
    mp.mp.dps = 30
    random.seed(1)
    worst = 0
    for trial in range(40):
        L = random.randint(1, 40)
        w = [random.choice('012')]
        while len(w) < L:
            w.append(random.choice([c for c in '012' if c != w[-1]]))
        w = ''.join(w)
        A, B = mp.mpf(random.random()), mp.mpf(random.random())
        S, cps = SH.develop(A, B, w)
        ref = cps[-1][2] - S
        x, y = mp.expj(A), mp.expj(B)
        h = sum(v * x ** p * y ** q for (p, q), v in hol(w).items()) / (2j)
        l2 = sum(v * x ** p * y ** q for (p, q), v in len2(w).items())
        worst = max(worst, abs(h - ref)); worst2 = max(globals().get("worst2", 0), abs(l2 - 4 * abs(ref) ** 2)); globals()["worst2"] = worst2
    print('max |hol - develop|, |len2 - 4|H|^2| over 40 random words:', mp.nstr(worst, 3), mp.nstr(worst2, 3))
