"""[Deposit] Section 5: the floating-point fan scan in the polygon P (shapes). Paper: "Finite perpendicular-beam partitions at irrational angles". Research name: data/scratch/s979/ftrace.py, copied verbatim except imports (deposit names), paths (data/ and figs/ under irrational_beams/) and this paragraph. Run through run_all.py."""
import _setup  # noqa: F401  (deposit: this directory on the path, cwd = irrational_beams/)
# s979 scratch: FLOAT P-tracer (fast) ⊕ adaptive fan scan. Same geometry as ptrace.py.
import math

def trace(al, s, maxstep=3000):
    r = 1 / math.tan(al)
    X0, X1 = r, s; Phi = math.atan2(s, r)
    a, b = 0, 1
    kcur = 0
    m = int(math.pi / (2 * al)) + 1
    ks = []
    for _ in range(maxstep):
        th = a * al + b * math.pi
        u0, u1 = math.cos(th), math.sin(th)
        best = None
        for k in range(kcur - m - 2, kcur + m + 3):
            if k == kcur: continue
            c, sn = math.cos(2 * k * al), math.sin(2 * k * al)
            nu = c * u0 + sn * u1
            if nu <= 0: continue
            lam = (r - (c * X0 + sn * X1)) / nu
            if lam <= 1e-13: continue
            Y0, Y1 = X0 + lam * u0, X1 + lam * u1
            PhY = Phi + math.atan2(X0 * Y1 - X1 * Y0, X0 * Y0 + X1 * Y1)
            if abs(PhY - 2 * k * al) <= al:
                if best is None or lam < best[0]: best = (lam, k, Y0, Y1, PhY)
        if best is None: return tuple(ks), 'lost'
        lam, k, X0, X1, Phi = best
        ks.append(k)
        if a - 2 * k == 0 and b % 2 == 0: return tuple(ks), 'headon'
        a, b = 4 * k - a, 1 - b
        kcur = k
    return tuple(ks), 'cap'

def steps(ks):
    return tuple([ks[0]] + [ks[j] - ks[j - 1] for j in range(1, len(ks))]) if ks else ()

def fans(x, N=200, tol=1e-9, maxstep=3000):
    al = x * math.pi / 2
    memo = {}
    def it(s):
        if s not in memo:
            ks, why = trace(al, s, maxstep)
            memo[s] = (steps(ks), why)
        return memo[s]
    pts = sorted(set([10.0 ** -e for e in range(12, 2, -1)] + [i / N for i in range(1, N)] + [1 - 10.0 ** -e for e in range(3, 13)]))
    out = []
    def refine(s0, s1):
        a, b = it(s0), it(s1)
        if a == b or s1 - s0 < tol: return
        mid = (s0 + s1) / 2
        refine(s0, mid); refine(mid, s1)
    for i in range(len(pts) - 1):
        refine(pts[i], pts[i + 1])
    ss = sorted(memo)
    for s in ss:
        st, why = memo[s]
        if out and out[-1][2] == st and out[-1][3] == why: out[-1][1] = s
        else: out.append([s, s, st, why])
    return out
