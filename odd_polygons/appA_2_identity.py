"""[Deposit] Proposition A.2: the loop relation as a formal identity. Paper: "Perpendicular billiards in rational polygons of odd denominator". Research name: probes/s938_pin_loop.py, copied verbatim except imports (deposit names), paths (data/ and figs/ under odd_polygons/) and this paragraph. Run through run_all.py; the usage lines below are the research repo's. Bracketed IDs and section marks (§sNNN) refer to the project's internal ledgers.

s938, queue `s937-a` (first half): Π_N's loop relation ℓ_m = c_{4m+5} + c_{4m+7} as a FORMAL identity, every odd N and m.
Π_N = (N−4, 4, N−2, 2)/N, sides [1, 2c−1, 2c, 2c] (c = cos 2π/N), closure exact. The loop word (measured N = 15 … 47 via the
walk, all rows; §s935 (2) to N = 103): W_m = (10)^{m+1} 1 (32)^{2m+2} 3 2 (01)^{m+1} from V0 at slot N − 14 − 8m — runs are
rotations about V1 (sides 0, 1) and V3 (sides 2, 3), so the developed endpoint is a formal ζ_{4N}-sum with exponents AN + Bm + C.
Result: 2ℓ − 2(c_{4m+5} + c_{4m+7}) reduces to the EMPTY sum (A ∈ {0,1} after ζ^{2N} = −1; distinct (A,B,C) are independent),
i.e. ℓ_m = 2[sin((4m+3)π/N) + sin((4m+5)π/N) + sin((4m+7)π/N) + sin((4m+9)π/N)] wherever W_m is realised. Legality of W_m for every odd N ≥ 8m + 15:
`probes/s938_pin_legal.py`. AUDITED `audits/s938_audit_pinloop.md` (25 V / 0 G / 0 F / 4 W).
PRIOR ART: §s935 (2) (the 10⁻¹⁵ float identity, "not proved"), §s937 (the tori); rulings grep *loop relation* → [NCYL-719]/[NCYL-720]
(homology partner, not the identity). DATA: none."""
import _setup  # noqa: F401  (deposit: this directory on the path, cwd = odd_polygons/)
import sys
import mpmath as mp
mp.mp.dps = 40


def add(u, v, s=1):
    r = dict(u)
    for e, c in v.items(): r[e] = r.get(e, 0) + s * c
    return {e: c for e, c in r.items() if c}
def sh(u, E): return {tuple(a + b for a, b in zip(e, E)): c for e, c in u.items()}
def cj(u): return {tuple(-a for a in e): c for e, c in u.items()}
def lin(a, b): return tuple(x + y for x, y in zip(a, b))
def mulc(E, n): return tuple(n * x for x in E)
def neg(E): return tuple(-x for x in E)
O = (0, 0, 0)
one = {O: 1}
twoc = {(0, 0, 4): 1, (0, 0, -4): 1}
L = [one, add(twoc, one, -1), twoc, twoc]
IDX = [O, (2, 0, -8), (2, 0, -4), (0, 0, -8)]          # side directions (units π/2N)
V = [{}]
for j in range(3): V.append(add(V[-1], sh(L[j], IDX[j])))
close = add(V[3], sh(L[3], IDX[3]))
# affine map (e, E, b): z ↦ ζ^E op_e(z) + b ; reflection in side j: z ↦ ζ^{2 idx_j} z̄ + (V_j − ζ^{2 idx_j} conj V_j)
def refl(j): return (1, mulc(IDX[j], 2), add(V[j], sh(cj(V[j]), mulc(IDX[j], 2)), -1))
def comp(F, G):            # F ∘ G
    e, E, b = F; e2, E2, b2 = G
    return (e ^ e2, lin(E, neg(E2)) if e else lin(E, E2), add(sh(cj(b2) if e else b2, E), b))
def rot_about(P, E):       # z ↦ ζ^E (z − P) + P
    return (0, E, add(P, sh(P, E), -1))
def rot_pow(j1, j2, n):    # n affine tuple
    common = set([j1, (j1 + 1) % 4]) & set([j2, (j2 + 1) % 4])
    P = V[common.pop()]
    E1 = lin(mulc(IDX[j1], 2), neg(mulc(IDX[j2], 2)))
    # E_total = n · E1 with n = (0, nm, nc), E1 = (A, 0, C): product only valid if E1 has no N part or n is constant
    A, B, C = E1
    A = A % 4                   # ζ^{4N} = 1
    nN, nm, nc = n
    assert nN == 0 and B == 0
    if A != 0:             # ζ^{A N} = i^A ... (ζ^{4N} = 1): keep A·n as N-coefficient only when n is constant
        assert nm == 0, ('N-part times symbolic m', E1, n)
    Et = (A * nc, C * nm if A == 0 else 0, C * nc)
    if A == 0: Et = (0, C * nm, C * nc)
    return rot_about(P, Et)


def red(u):                # ζ^{2N} = −1: A mod 2 with sign (A odd kept as ζ^N = i)
    out = {}
    for (A, B, C), c in u.items():
        q, r = divmod(A, 2)
        e = (r, B, C)
        out[e] = out.get(e, 0) + c * (1 if q % 2 == 0 else -1)
    return {e: c for e, c in out.items() if c}


def loop_len():
    """4ℓ... : w = conj(ω0)(F(V0) − V0), S = w + w̄ = 2ℓ (vertices unscaled here: side 1 = 1)"""
    m = (0, 1, 0)
    F = (0, O, {})
    F = comp(F, rot_pow(1, 0, lin(m, (0, 0, 1))))                  # (10)^{m+1}
    F = comp(F, refl(1))                                           # 1
    F = comp(F, rot_pow(3, 2, lin(mulc(m, 2), (0, 0, 2))))         # (32)^{2m+2}
    F = comp(F, refl(3)); F = comp(F, refl(2))                     # 3 2
    F = comp(F, rot_pow(0, 1, lin(m, (0, 0, 1))))                  # (01)^{m+1}
    e, E, b = F
    end = add(sh(cj(V[0]) if e else V[0], E), b)
    d0 = (1, -8, -14)                                              # slot N − 14 − 8m
    w = sh(add(end, V[0], -1), neg(d0))
    return red(add(w, cj(w)))


def main():
    S = loop_len()
    mi = (3, 0, 0)
    def two_sin(b):
        return red(sh({mulc(b, 2): 1, mulc(b, -2): -1}, mi))
    def mul(u, v):
        out = {}
        for e1, c1 in u.items():
            for e2, c2 in v.items():
                e = lin(e1, e2)
                out[e] = out.get(e, 0) + c1 * c2
        return red({e: c for e, c in out.items() if c})
    rhs = add(mul(twoc, two_sin((0, 4, 5))), mul(twoc, two_sin((0, 4, 7))))
    assert red(close) == {}, 'Π_N does not close formally'
    diff = red(add(S, add(rhs, rhs), -1))
    assert diff == {}, diff
    sys.path.insert(0, 'probes')
    import s4_7_witnesses as W
    from poly_slots import Poly
    nchk = 0
    for N in (15, 17, 23, 31, 33, 39, 41, 47):
        cc = mp.cos(2 * mp.pi / N)
        Pn = Poly([N - 4, 4, N - 2, 2], N, lengths=[1, 2 * cc - 1, 2 * cc, 2 * cc])
        for q in W.pairs(Pn, mp.mpf(10) ** -20)[0]:
            if q['loop']:
                m = (N - 14 - q['slot'][1]) // 8
                want = ('10' * (m + 1)) + '1' + ('32' * (2 * m + 2)) + '32' + ('01' * (m + 1))
                assert ''.join(map(str, q['word'])) == want and q['slot'][1] == N - 14 - 8 * m, (N, q['slot'])
                nchk += 1
    print('formal identity ℓ_m = c_{4m+5} + c_{4m+7}: PROVED (empty difference);', nchk, 'walk loops match W_m (N ≤ 47)')


if __name__ == '__main__':
    main()
