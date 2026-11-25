"""[Deposit] Lemma 4.1 (Appendix B), odd J. Paper: "Finite perpendicular-beam partitions at irrational angles". Research name: probes/s976_fmJ_irrational.py, copied verbatim except imports (deposit names), paths (data/ and figs/ under irrational_beams/) and this paragraph. Run through run_all.py.

s976 — queue s974-g: the root of F(m, J) at ODD J ≥ 3 (b_J = 1, §s976 (1)) is an IRRATIONAL multiple of π, every m ≥ 2.
Same proof as probes/s975_fmj_irrational.py (even J = 2j, [FMJ-IRR], AUDITED s976) with two changes:
 (1) relation: (−1)^J cos α sin Nα = sin((m−1)α), N = (2J+1)m; odd J ⇒ z^{N+1} + z^{N−1} − z^{−N−1} − z^{−N+1} + 2z^{m−1} − 2z^{1−m} = 0;
 (0) location 0 < π/(2m) − α < π/(2m(2mJ+1)) ⇒ r(2mJ+1) < Q, r = Q − 4mp' ≥ 1 ⇒ Q ≥ 2mJ + 2 ≥ 6m + 2 (J ≥ 3).
 (3) the m-linear N-free forms have |vm + w| ≤ 4m < 6m + 2 ≤ Q unless zero; the m = 2 exception (if present): Q | N = 4J + 2 ⇒ Q = 4J + 2,
     r even ≥ 2 ⇒ r(4J+1) > Q — contradiction (by hand, as in s975).
 (4) finite check: Q ≤ c_max and Q ≥ 2mJ + 2 over odd J ≥ 3, interval arithmetic.
J = 1 (α_n) is §s973 (6)'s hand proof (the bound 4m < 2m + 2 fails there).
Not PRE-REGISTERED (a proof's computer steps). PRIOR ART: §s975b, audits/s976_audit_fmj_irr.md. DATA: none.
AUDITED s976: audits/s976_audit_fmj_all.md.
Usage: python probes/s976_fmJ_irrational.py
"""
import _setup  # noqa: F401  (deposit: this directory on the path, cwd = irrational_beams/)
import itertools
import math
from fractions import Fraction as F

import sympy as sp
from mpmath import mp, iv
from sympy.matrices.normalforms import hermite_normal_form

E = {1: (1, 0, 1), 2: (1, 0, -1), 3: (-1, 0, -1), 4: (-1, 0, 1), 5: (0, 1, -1), 6: (0, -1, 1)}      # u·N + v·m + w
A = {1: 1, 2: 1, 3: -1, 4: -1, 5: 2, 6: -2}      # odd J: cos α sin Nα + sin((m−1)α) = 0


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


def finite_check(cmax):
    iv.dps = 60
    n, mn = 0, None
    for m in range(2, cmax):
        for J in range(3, cmax, 2):
            if 2 * m * J + 2 > cmax:
                continue
            N = (2 * J + 1) * m
            lo, hi = F(1, 4 * m) - F(1, 4 * m * (2 * m * J + 1)), F(1, 4 * m)
            for Q in range(1, cmax + 1):
                for p in range(1, Q):
                    if math.gcd(p, Q) == 1 and lo < F(p, Q) < hi:
                        a = 2 * iv.pi * p / Q
                        v = iv.cos(a) * iv.sin(N * a) + iv.sin((m - 1) * a)
                        assert not (v.a <= 0 <= v.b), (m, J, p, Q)
                        low = min(abs(float(v.a)), abs(float(v.b)))
                        n += 1
                        mn = low if mn is None else min(mn, low)
    return n, mn


if __name__ == '__main__':
    S = list(structures())
    kinds = {'const': [], 'lin': []}
    for classes, sub, cons in S:
        k, val = nfree(cons)
        kinds[k].append((classes, sub, val))
    cmax = int(max(v for _, _, v in kinds['const']))
    assert len(S) == 201 and len(kinds['const']) == 191 and cmax <= 120          # s976 audit row 19
    norm = lambda v: tuple(int(x) for x in (v if v[0] > 0 else (-v[0], -v[1])))
    lins = sorted({norm(v) for _, _, v in kinds['lin']})
    print(f"{len(S)} structures: {len(kinds['const'])} constant (max {cmax}), {len(kinds['lin'])} linear, forms {lins}")
    for v, w in lins:
        assert abs(v) <= 4 and abs(w) <= 8, (v, w)                      # |vm + w| ≤ 4m < 6m + 2 for m ≥ 2 unless zero
        zeros = [m for m in range(2, 50) if v * m + w == 0]
        assert zeros in ([], [2]) and (not zeros or (v, w) == (4, -8)), (v, w, zeros)
        if zeros:
            print('  form', (v, w), 'vanishes at m =', zeros, 'in structures',
                  [c for c, s, vv in kinds['lin'] if norm(vv) == (v, w)])
    exc = [c for c, s_, vv in kinds['lin'] if norm(vv) == (4, -8)]
    assert exc == [[[1], [3], [4, 5], [2, 6]]], exc                             # E4 − E5 = −N at m = 2: Q | N (hand step)
    n, mn = finite_check(cmax)
    assert n == 254 and mn > 1e-6
    print(f"finite check (odd J ≥ 3, Q ≤ {cmax}): {n} rationals, min |E| ≥ {mn:.4g} (interval)")
