"""[Deposit] Theorem A(ST): checks of the uniqueness and beam-count proofs (Section 4.2). Paper: "Finite perpendicular-beam partitions at irrational angles". Research name: probes/s1009_st_count.py, copied verbatim except imports (deposit names), paths (data/ and figs/ under irrational_beams/) and this paragraph. Run through run_all.py.

s1009 — Theorem ST loose ends (queue s979-c; referee FA R7): the window root of T_j is UNIQUE and the branch count is EXACTLY 2j,
for every j >= 1, m >= max(2, 2j-1) ([ST-COUNT]; proof in general_triangles.md §s1009). This keeper re-checks every step of the hand
proof against the raw recursion (eq:recursion of paper A); it proves nothing by itself.

PRE-REGISTERED checks (all must HOLD; C5 is the control and must FIRE):
  C1  closed form  2 M_{2k} = kappa sin(k a) [sin((k-1)a + eps) - cos a sin(k a)]  vs the recursion (Sigma_nu summed from psi_nu), random
      (j, m, a) in the window, 1 <= k <= j-1: max error < 1e-30.
  C2  sign identity at the root:  2 M_{2j-2} = -(kappa sin((j-1)a) sin a / (2 sin(j a))) * Br,
      Br := sin((2j-1)a + eps) - sin((2j-1)a) + sin(a - eps), and Br > 0, at every root (j <= 40, six m >= 2j-1 each).
  C3  R(a) = [sin ja + sin((j+1)a) sin((j-1)a)/sin ja]/(2 cos a) equals (2 sin^2 ja - sin^2 a)/(2 cos a sin ja), and is increasing
      on (0, pi/(4j-2)] (j < 60, sampled); T_j = 2 cos a sin ja (sin y - R), y = eps + (j-1)a (random points).
  C4  exactly one sign change of T_j on a fine grid of the window (j <= 40, same m's) — the uniqueness the proof gives.
  C5  CONTROL: at m = 2j - 2 (j >= 3), where (4j-2)a + eps > pi, M_{2j-2} > 0 at the root (it is the centre of D'_j) — so the sign
      test of C2 could have come out the other way.
PRIOR ART: general_triangles.md §s979 (2) ("Not proved here (MEASURED 326/326): the window root is UNIQUE; the branch count is 2j —
  ... M_{2k} <= 0 ... reduces to k a <= |theta|"); rulings.py --grep "theorem ST" → nothing on ST; queue s979-c. The reduction to
  M_{2k} <= 0 and the theta form are §s979's; the root substitution (C2) and the monotone R (C3) are new.
DATA: none (data_index --grep ST_ → no store field; the 326-root check of §s979 lives in probes/s979_* outputs, not reused).
Usage: python probes/s1009_st_count.py
"""
import _setup  # noqa: F401  (deposit: this directory on the path, cwd = irrational_beams/)
import random
from mpmath import mp, mpf, sin, cos, cot, pi

mp.dps = 40


def descent(a, m, n):
    """psi_nu, Sigma_nu of the universal descent (steps m - nu + 1) from the recursion; p_nu = (-1)^nu (s - Sigma_nu)."""
    r = cot(a)
    psi, Sig = [mpf(0)], [mpf(0)]
    for nu in range(1, n + 1):
        ps = pi - 2 * (m - nu + 1) * a - psi[-1]
        psi.append(ps)
        Sig.append(Sig[-1] + (-1) ** (nu + 1) * 2 * r * sin(ps))
    return psi, Sig


def M(Sig, k):
    return (Sig[2 * k] + Sig[2 * k - 1]) / 2


def Tj(a, j, m):
    return 2 * cos(a) * sin(j * a) * sin((2 * m - j + 1) * a) + sin(a) ** 2 - 2 * sin(j * a) ** 2


def R(a, j):
    return (sin(j * a) + sin((j + 1) * a) * sin((j - 1) * a) / sin(j * a)) / (2 * cos(a))


def root(j, m):
    lo, hi = pi / (2 * m + 1), pi / (2 * m)
    flo = Tj(lo, j, m)
    for _ in range(160):
        mid = (lo + hi) / 2
        fm = Tj(mid, j, m)
        if (fm > 0) == (flo > 0):
            lo, flo = mid, fm
        else:
            hi = mid
    return (lo + hi) / 2


def main():
    random.seed(1009)
    ok = {}
    # C1
    e1 = mpf(0)
    for _ in range(300):
        j = random.randint(2, 25)
        m = random.randint(2 * j - 1, 6 * j)
        a = pi / (2 * m + 1) + (pi / (2 * m) - pi / (2 * m + 1)) * mpf(random.random())
        eps = pi - 2 * m * a
        kap = 4 * cos(a) / sin(a) ** 2
        _, Sig = descent(a, m, 2 * j)
        for k in range(1, j):
            e1 = max(e1, abs(2 * M(Sig, k) - kap * sin(k * a) * (sin((k - 1) * a + eps) - cos(a) * sin(k * a))))
    ok['C1'] = e1 < mpf(10) ** -30
    # C2, C4
    e2, minBr, n, bad4 = mpf(0), mpf(1), 0, 0
    for j in range(2, 41):
        for m in sorted({2 * j - 1, 2 * j, 2 * j + 1, 3 * j, 4 * j, 10 * j}):
            a = root(j, m)
            eps = pi - 2 * m * a
            kap = 4 * cos(a) / sin(a) ** 2
            _, Sig = descent(a, m, 2 * j)
            Br = sin((2 * j - 1) * a + eps) - sin((2 * j - 1) * a) + sin(a - eps)
            e2 = max(e2, abs(2 * M(Sig, j - 1) + kap * sin((j - 1) * a) * sin(a) / (2 * sin(j * a)) * Br))
            minBr = min(minBr, Br)
            n += 1
            lo, hi = pi / (2 * m + 1), pi / (2 * m)
            g = [Tj(lo + (hi - lo) * t / 3000, j, m) for t in range(1, 3000)]
            bad4 += (sum(1 for i in range(len(g) - 1) if g[i] * g[i + 1] < 0) != 1)
    ok['C2'] = e2 < mpf(10) ** -25 and minBr > 0
    ok['C4'] = bad4 == 0
    # C3
    e3, dec = mpf(0), 0
    for _ in range(200):
        j = random.randint(1, 30)
        m = random.randint(max(2, 2 * j - 1), 6 * j + 2)
        a = pi / (2 * m + 1) + (pi / (2 * m) - pi / (2 * m + 1)) * mpf(random.random())
        y = pi - 2 * m * a + (j - 1) * a
        e3 = max(e3, abs(R(a, j) - (2 * sin(j * a) ** 2 - sin(a) ** 2) / (2 * cos(a) * sin(j * a))),
                 abs(Tj(a, j, m) - 2 * cos(a) * sin(j * a) * (sin(y) - R(a, j))))
    for j in range(2, 60):
        v = [R(pi / (4 * j - 2) * t / 400, j) for t in range(1, 401)]
        dec += sum(1 for i in range(399) if v[i + 1] <= v[i])
    ok['C3'] = e3 < mpf(10) ** -30 and dec == 0
    # C5 control
    fired = all(M(descent(root(j, 2 * j - 2), 2 * j - 2, 2 * j)[1], j - 1) > 0 for j in range(3, 41))
    ok['C5 (control fires)'] = fired
    print(f'C1 max err {mp.nstr(e1, 3)} | C2 {n} roots, identity err {mp.nstr(e2, 3)}, min Br {mp.nstr(minBr, 4)} | '
          f'C3 err {mp.nstr(e3, 3)}, R non-increasing steps {dec} | C4 roots != 1: {bad4}')
    for k, v in ok.items():
        print(f'  {k}: {"HOLD" if v else "FAIL"}')
    assert all(ok.values())


if __name__ == '__main__':
    main()
