"""[Deposit] Helper, imported by the scripts of this directory. Paper: "Perpendicular billiards in rational polygons of odd denominator". Research name: probes/s946_fold_survivors.py, copied verbatim except imports (deposit names), paths (data/ and figs/ under odd_polygons/) and this paragraph. Run through run_all.py; the usage lines below are the research repo's. Bracketed IDs and section marks (§sNNN) refer to the project's internal ledgers.

[promoted s958 from data/scratch/s946/audit2/lf.py, verbatim but this line; imported by probes/s952_s1_family.py] audit item 5: brute-force the LOCAL constraints behind Lemma F (independent of its case analysis).
A reflection tiling of convex Pi (corner angles B, units pi/N) by T = (t0, t1, t2): every vertex point carries ONE label;
a point of label i and angle alpha (corner B, flat boundary N, interior 2N) has alpha/t_i tiles (integer; interior: even,
i.e. N/t_i integer). Each tile has one vertex of each label => n = sum_{P of label i} alpha_P / t_i for every i.
So: each t_i is B/m or N/m; corners labelled with B/t ∈ Z>=1; n_i = c_i + (f_i + 2 g_i) N/t_i equal for i = 0,1,2, n >= 2,
every label used. Reports every (T, corner labelling) that survives, with the forced n (or 'free' if all three can grow)."""
import _setup  # noqa: F401  (deposit: this directory on the path, cwd = odd_polygons/)
from fractions import Fraction as Fr
import itertools, sys, json

def survivors(corners, N, M=400):
    S = set()
    for B in corners:
        for m in range(1, M):
            S.add(Fr(B, m))
    for m in range(2, M):
        S.add(Fr(N, m))
    Sl = sorted(S, reverse=True)
    out = []
    seen = set()
    for t0 in Sl:
        for t1 in Sl:
            if t1 > t0:
                continue
            t2 = N - t0 - t1
            if t2 <= 0 or t2 > t1 or t2 not in S:
                continue
            T = (t0, t1, t2)
            if T in seen:
                continue
            seen.add(T)
            for lab in itertools.product(range(3), repeat=len(corners)):
                if any((Fr(B) / T[l]).denominator != 1 for B, l in zip(corners, lab)):
                    continue
                c = [sum(Fr(B) / T[i] for B, l in zip(corners, lab) if l == i) for i in range(3)]
                q = [Fr(N) / T[i] for i in range(3)]
                free = [qq.denominator == 1 for qq in q]
                # n must be >= max c_i, n ≡ c_i mod q_i for free labels (q_i integer, step q_i: f flat adds q_i, interior adds 2q_i),
                # n == c_i for non-free; each label used
                fixed = [c[i] for i in range(3) if not free[i]]
                if fixed and len(set(fixed)) > 1:
                    continue
                def ok(n):
                    if n < 2:
                        return False
                    for i in range(3):
                        if free[i]:
                            r = n - c[i]
                            if r < 0 or r % q[i] != 0:
                                return False
                            if c[i] == 0 and r == 0:
                                return False
                        else:
                            if n != c[i] or c[i] == 0:
                                return False
                    return True
                if fixed:
                    n = fixed[0]
                    if ok(n):
                        out.append((T, lab, int(n)))
                else:
                    ns = [n for n in range(2, 2000) if ok(n)]
                    if ns:
                        out.append((T, lab, 'free:%s...' % ns[:5]))
    return out

