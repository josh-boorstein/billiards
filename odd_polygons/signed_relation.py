"""[Deposit] Helper, imported by the scripts of this directory. Paper: "Perpendicular billiards in rational polygons of odd denominator". Research name: engine/signed_relation.py, copied verbatim except imports (deposit names), paths (data/ and figs/ under odd_polygons/) and this paragraph. Run through run_all.py; the usage lines below are the research repo's. Bracketed IDs and section marks (§sNNN) refer to the project's internal ledgers.

signed_relation — EXACT test of Theorem P's hypothesis: is there ε ∈ {−1,0,1}^p ∖ 0 with Σ εᵢ vᵢ = 0, for integer vectors vᵢ
(the swapped-pair holonomies as reduced coefficient vectors in ℤ[ζ_{4N}], or any ℤ^D coordinates of the lengths)? (s940, `s938-c`)

The relations live in the ℚ-kernel K of the D × p matrix [v₁ … v_p]. Compute K exactly (python-flint `fmpz_mat.nullspace`), put a
basis in reduced row-echelon form, so a kernel vector is determined by its values on the k FREE coordinates (the non-pivots of the
echelon form of Kᵀ, where Kᵀ's rows span K), and every other coordinate is a fixed ℚ-linear form in them. A signed relation = an
assignment of {−1,0,1} to the free coordinates, not all 0, making every dependent coordinate an integer in {−1,0,1}. Searched by
depth-first branch-and-bound over the free coordinates (prune when a dependent form can no longer reach [−1, 1]); `k = 0` ⇒ no
relation at once, with no search. Replaces `s927_free_region.sign_relation`'s 3^p brute force / float meet-in-the-middle, which
does not scale past p ≈ 20.

  qkernel(V)                 V = list of p integer vectors (length D) → (k, B): B = k × p list of Fraction rows spanning the
                             ℚ-kernel, RREF (pivot columns carry the identity). k = p − rank.
  signed_relation(V, limit)  → a {−1,0,1} list ε ≠ 0 with Σ εᵢ vᵢ = 0 (verified exactly), or None; raises if more than `limit`
                             search nodes are needed (no silent truncation).
  rank_q(V)                  the ℚ-rank of the vᵢ.
GATE (`_selftest`): agrees with brute force on random small instances (with planted relations and without); the s939 tightness
control cos(jπ/210), j odd ≤ 49 — rank 24, unique relation in {−1,0,1} — reproduced from the exact ζ₈₄₀ coordinates.
"""
import _setup  # noqa: F401  (deposit: this directory on the path, cwd = odd_polygons/)
from fractions import Fraction

import flint


def qkernel(V):
    p = len(V)
    if p == 0:
        return 0, []
    D = len(V[0])
    M = flint.fmpz_mat(D, p, [int(V[j][i]) for i in range(D) for j in range(p)])
    X, k = M.nullspace()
    if k == 0:
        return 0, []
    Kt = flint.fmpq_mat(k, p, [X[i, c] for c in range(k) for i in range(p)])     # rows span the kernel
    R, r = Kt.rref()
    assert r == k
    B = [[Fraction(int(R[i, j].p), int(R[i, j].q)) for j in range(p)] for i in range(k)]
    return k, B


def rank_q(V):
    if not V:
        return 0
    D = len(V[0])
    return flint.fmpz_mat(D, len(V), [int(V[j][i]) for i in range(D) for j in range(len(V))]).rank()


def _verify(V, eps):
    D = len(V[0])
    return any(eps) and all(sum(e * v[i] for e, v in zip(eps, V)) == 0 for i in range(D))


def signed_relation(V, limit=10 ** 7):
    k, B = qkernel(V)
    if k == 0:
        return None
    p = len(V)
    piv = []
    for row in B:
        piv.append(next(j for j in range(p) if row[j] != 0))
    # x = Σ_i a_i B_i with a_i = x[piv_i] (RREF: B_i[piv_l] = δ_il), so every coordinate j is Σ_i a_i B_i[j]
    dep = [j for j in range(p) if j not in piv]
    coef = {j: [B[i][j] for i in range(k)] for j in dep}
    # order free coordinates by how many dependent forms they touch (most constraining first)
    order = sorted(range(k), key=lambda i: -sum(1 for j in dep if coef[j][i] != 0))
    rem = {j: [sum(abs(coef[j][order[t]]) for t in range(s, k)) for s in range(k + 1)] for j in dep}
    a = [0] * k
    nodes = [0]
    found = []

    def dfs(s, part):
        nodes[0] += 1
        if nodes[0] > limit:
            raise RuntimeError(f'signed_relation: search exceeded {limit} nodes (k = {k}, p = {p})')
        for j in dep:
            r = rem[j][s]
            if part[j] - r > 1 or part[j] + r < -1:
                return False
        if s == k:
            if not any(a):
                return False
            if all(part[j].denominator == 1 for j in dep):
                found.append(list(a))
                return True
            return False
        i = order[s]
        for v in (0, 1, -1):
            if s == 0 and v == -1:
                continue                           # ε and −ε: fix the sign of the first nonzero free value... (see below)
            a[i] = v
            np_ = {j: part[j] + v * coef[j][i] for j in dep} if v else part
            if dfs(s + 1, np_):
                return True
        a[i] = 0
        return False

    # sign symmetry: any relation can be negated, so the first free coordinate (in `order`) may be taken in {0, 1}; if it is 0
    # the recursion continues freely — the −1 branch at s = 0 is redundant, at s > 0 it is not
    dfs(0, {j: Fraction(0) for j in dep})
    if not found:
        return None
    a = found[0]
    eps = [0] * p
    for i in range(k):
        eps[piv[i]] = a[i]
    for j in dep:
        eps[j] = int(sum(a[i] * coef[j][i] for i in range(k)))
    assert _verify(V, eps), 'signed_relation: candidate failed exact verification'
    return eps


def _brute(V):
    import itertools
    p = len(V)
    for eps in itertools.product((-1, 0, 1), repeat=p):
        if any(eps) and next(e for e in eps if e) > 0 and _verify(V, eps):
            return list(eps)
    return None


def _selftest():
    import random
    random.seed(940)
    n_rel = 0
    for trial in range(300):
        p = random.randint(1, 7)
        D = random.randint(1, 6)
        V = [[random.randint(-3, 3) for _ in range(D)] for _ in range(p)]
        if random.random() < 0.5 and p >= 2:          # plant a relation
            eps = [random.choice((-1, 0, 1)) for _ in range(p - 1)]
            V[-1] = [-sum(e * v[i] for e, v in zip(eps, V[:-1])) for i in range(D)]
        a, b = signed_relation(V), _brute(V)
        assert (a is None) == (b is None), (V, a, b)
        n_rel += a is not None
    # the s939 tightness control: cos(jπ/210), j = 1, 3, …, 2M − 1, as exact vectors in ℚ(ζ_840) (power basis mod Φ_840)
    import sympy as sp
    z = sp.symbols('z')
    Phi = sp.Poly(sp.cyclotomic_poly(840, z), z)
    D = Phi.degree()

    def vec(j):                                     # 2cos(jπ/210) = ζ^{2j} + ζ^{−2j}, ζ = e^{iπ/420}... use ζ_840 = e^{2πi/840}
        e1, e2 = (2 * j) % 840, (-2 * j) % 840
        r = sp.Poly(z ** e1 + z ** e2, z).rem(Phi)
        c = r.all_coeffs()[::-1]
        return [int(x) for x in c] + [0] * (D - len(c))
    for M_, want_rank, want_rel in ((24, 24, False), (25, 24, True)):
        V = [vec(j) for j in range(1, 2 * M_, 2)]
        assert rank_q(V) == want_rank
        eps = signed_relation(V)
        assert (eps is not None) == want_rel, (M_, eps)
    print(f'signed_relation selftest OK ({n_rel}/300 random with relations; r = 105 control: rank 24, relation at M = 25: {eps})')


if __name__ == '__main__':
    _selftest()
