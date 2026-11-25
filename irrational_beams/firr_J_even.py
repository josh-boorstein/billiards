"""[Deposit] Lemma 4.1 (Appendix B), even J. Paper: "Finite perpendicular-beam partitions at irrational angles". Research name: probes/s975_fmj_irrational.py, copied verbatim except imports (deposit names), paths (data/ and figs/ under irrational_beams/) and this paragraph. Run through run_all.py.

s975 — queue s973-k: the root α of Theorem F's equation E(m, j) (§s974 (1)) is an IRRATIONAL multiple of π, for every m ≥ 2, j ≥ 1.

STATEMENT.  E(m, j): cos α · sin((4j+1)mα) = sin((m−1)α); Theorem F: exactly one root α = π/(2m) − δ, 0 < δ < δ* = π/(2m(4mj+1)).
Claim: α/π ∉ ℚ.
PROOF (Mann ⊕ a finite case split; this file checks every finite step by assertion).
 (0) Location.  If α = 2πp'/Q (gcd(p', Q) = 1), then p'/Q ∈ I := (1/(4m) − 1/(4m(4mj+1)), 1/(4m)); writing 1/(4m) − p'/Q = r/(4mQ),
     r = Q − 4mp' ≥ 1, the bound r(4mj+1) < Q gives Q ≥ 4mj + 2.
 (1) The relation.  With z = e^{iα} (order Q) and N = (4j+1)m, 2i·[sin(N+1)α + sin(N−1)α − 2 sin(m−1)α] = 0 reads
       z^{N+1} + z^{N−1} − z^{−N−1} − z^{−N+1} − 2z^{m−1} + 2z^{1−m} = 0            (six terms, exponents E_i, coefficients a_i).
 (2) Structures.  Group equal terms into classes (z^{E_i} = z^{E_k} ⇒ Q | E_i − E_k); drop classes with zero total coefficient; split
     the rest into MINIMAL vanishing subsums.  Size 1: impossible.  Size 2: c₁x + c₂y = 0 ⇒ |c₁| = |c₂|, y/x = ±1 ⇒ Q | 2(E_i − E_k).
     Size 3: a minimal 3-term relation with rational coefficients is a rotation of 1 + ω + ω² up to signs (Mann: ratios in μ₆; the
     non-real ones force equal |coefficients|) ⇒ equal |c| and Q | 6(E_i − E_k).  Size 4–6: Mann's theorem (H. B. Mann, Mathematika
     12 (1965) 107–117): ratios of a minimal relation with k terms are roots of unity of order dividing ∏_{p ≤ k} p | 30 ⇒ Q | 30(E_i − E_k).
     Enumerated: 201 admissible structures; each gives a set of linear forms u·N + v·m + w that Q divides.
 (3) Bounds.  For each structure the N-free part of the lattice of its forms (u = 0) contains either a nonzero constant c ≤ 120
     (191 structures ⇒ Q ≤ 120) or a form ∈ {4(m−1), 2(m−1), 4m, 4(m−2)} (10 structures), each < 4mj + 2 ≤ Q unless it vanishes —
     only 4(m−2) at m = 2, where the structure has classes {E1, E5}, {E3, E6} ⇒ Q | N = 8j + 2 ⇒ Q = 8j + 2, p' ≤ j, r = Q − 8p' ≥ 2,
     r(8j+1) ≥ 2(8j+1) > Q: contradicts (0).
 (4) Finite check.  Q ≤ 120 and Q ≥ 4mj + 2 ⇒ mj ≤ 29; every p'/Q ∈ I with Q ≤ 120, mj ≤ 29 has |E| ≥ 2.5e-5 (interval arithmetic, 60 digits; 319 rationals): not a root. ∎
 AUDITED s976: audits/s976_audit_fmj_irr.md.
Not PRE-REGISTERED (a proof's computer steps).  PRIOR ART: §s973 (6) (the α_n Galois argument; does not transfer), queue s973-k;
rulings --grep 'Mann' / 'roots of unity' → none on this family.  DATA: none (deterministic).
Usage: python probes/s975_fmj_irrational.py
"""
import _setup  # noqa: F401  (deposit: this directory on the path, cwd = irrational_beams/)
import itertools
import math
from fractions import Fraction as F

import sympy as sp
from mpmath import mp, iv
from sympy.matrices.normalforms import hermite_normal_form

E = {1: (1, 0, 1), 2: (1, 0, -1), 3: (-1, 0, -1), 4: (-1, 0, 1), 5: (0, 1, -1), 6: (0, -1, 1)}      # u·N + v·m + w
A = {1: 1, 2: 1, 3: -1, 4: -1, 5: -2, 6: 2}


def setparts(s):
    if not s:
        yield []
        return
    f, r = s[0], s[1:]
    for p in setparts(r):
        yield [[f]] + p
        for i in range(len(p)):
            yield p[:i] + [[f] + p[i]] + p[i + 1:]


def structures():
    for classes in setparts(list(E)):
        cons = [(i, k, 1) for C in classes for i, k in itertools.combinations(C, 2)]
        live = [C for C in classes if sum(A[i] for i in C) != 0]
        for sub in setparts(list(range(len(live)))):
            c2, ok = list(cons), True
            for S in sub:
                coefs = [abs(sum(A[i] for i in live[s])) for s in S]
                if len(S) == 1 or (len(S) in (2, 3) and len(set(coefs)) > 1):
                    ok = False
                    break
                mult = {2: 2, 3: 6}.get(len(S), 30)
                c2 += [(live[s][0], live[t][0], mult) for s, t in itertools.combinations(S, 2)]
            if ok:
                yield classes, sub, c2


def nfree(cons):
    """Generators of (a sublattice of) L ∩ {u = 0}; returns ('const', c) or ('lin', (v, w))."""
    vecs = [[mult * (a - b) for a, b in zip(E[i], E[k])] for i, k, mult in cons]
    B = [list(r) for r in hermite_normal_form(sp.Matrix(vecs).T).T.tolist() if any(r)]
    L0 = [r[1:] for r in B if r[0] == 0]
    for r1, r2 in itertools.combinations([r for r in B if r[0] != 0], 2):
        l = sp.ilcm(r1[0], r2[0])
        L0.append([l // r1[0] * x - l // r2[0] * y for x, y in zip(r1[1:], r2[1:])])
    L0 = [v for v in L0 if any(v)]
    M0 = sp.Matrix(L0)
    if M0.rank() == 2:
        g0 = [list(r) for r in hermite_normal_form(M0.T).T.tolist() if any(r)]
        return 'const', abs(sp.Matrix(g0).det()) // sp.igcd(*[r[0] for r in g0])
    v0 = [r for r in L0 if r[0] == 0]
    if v0:
        return 'const', abs(v0[0][1])
    g = L0[0]
    return 'lin', (g[0], g[1])


def finite_check():
    """Interval arithmetic (mpmath iv, 60 digits): every enclosure of E must exclude 0 (s976 audit row 13)."""
    iv.dps = 60
    n, mn = 0, None
    for m in range(2, 30):
        for j in range(1, 30):
            if m * j > 29:
                continue
            N = (4 * j + 1) * m
            lo, hi = F(1, 4 * m) - F(1, 4 * m * (4 * m * j + 1)), F(1, 4 * m)
            for Q in range(1, 121):
                for p in range(1, Q):
                    if math.gcd(p, Q) == 1 and lo < F(p, Q) < hi:
                        a = 2 * iv.pi * p / Q
                        v = iv.cos(a) * iv.sin(N * a) - iv.sin((m - 1) * a)
                        assert not (v.a <= 0 <= v.b), (m, j, p, Q)
                        low = min(abs(float(v.a)), abs(float(v.b)))
                        n += 1
                        mn = low if mn is None else min(mn, low)
    return n, mn


if __name__ == '__main__':
    S = list(structures())
    assert len(S) == 201
    kinds = {'const': [], 'lin': []}
    for classes, sub, cons in S:
        k, val = nfree(cons)
        kinds[k].append((classes, sub, val))
    assert len(kinds['const']) == 191 and max(kinds['const'], key=lambda x: x[2])[2] <= 120
    norm = lambda v: tuple(int(x) for x in (v if v[0] > 0 else (-v[0], -v[1])))
    lins = sorted({norm(v) for _, _, v in kinds['lin']})
    assert set(lins) <= {(4, -4), (4, -8), (4, 0), (2, -2)}, lins          # 4(m−1), 4(m−2), 4m, 2(m−1)
    zero_at = [(c, s) for c, s, v in kinds['lin'] if norm(v) == (4, -8)]
    assert len(zero_at) == 1 and zero_at[0][0] == [[2], [4], [1, 5], [3, 6]]
    for m in range(2, 200):                                                 # the linear forms stay below 4mj+2 ≤ Q
        for v, w in lins:
            assert v * m + w == 0 or abs(v * m + w) < 4 * m + 2
    n, mn = finite_check()
    assert n == 319 and mn > 1e-6
    print(f"201 structures: 191 with a constant ≤ 120, 10 linear in m (vanishing only 4(m−2) at m = 2, excluded by hand); "
          f"finite check: {n} rationals p'/Q (Q ≤ 120, mj ≤ 29) in the root intervals, min |E| ≥ {mn:.4g} (interval) ⇒ α/π ∉ ℚ for all m, j.")
