"""[Deposit] Proposition A.2: legality of the loop word. Paper: "Perpendicular billiards in rational polygons of odd denominator". Research name: probes/s938_pin_legal.py, copied verbatim except imports (deposit names), paths (data/ and figs/ under odd_polygons/) and this paragraph. Run through run_all.py; the usage lines below are the research repo's. Bracketed IDs and section marks (§sNNN) refer to the project's internal ledgers.

s938, queue `s937-a`: LEGALITY of Π_N's loop word W_m for EVERY odd N ≥ 8m + 15, m ≥ 0 — so with `s938_pin_loop` (the
formal length identity) ℓ_m = c_{4m+5} + c_{4m+7} is PROVED for all such (N, m). Computer-assisted.
W_m = [1010…1 (2m+3) about V1][3232…32 (4m+6) about V3][0101…01 (2m+2) about V1], from V0 at slot N − 14 − 8m. Π_N is convex and
every exit side is a ray of the run's centre P, so a step is legal iff f(P)·f(ray end) < 0 (f = 2 Im(conj ω0 (· − V0))); per run
the crossed rays of one parity form an arc (L cos ψ > d concave; spans (2m+2)·4π/N, (4m+5)·2π/N, (2m+1)·4π/N < π ⟸ N ≥ 8m + 15, `span_ok`), so the
extremes t = 0, 1, last-even, last-odd suffice; run-to-run transitions are the next run's t = 0; launch: V1, V2 ahead; the end at
V0 is a formal identity (empty sum). Each condition is certified on x = π/2N ∈ (0, π/30], y = mx: m ≥ 1, x ≤ y ≤ π/16 − 1.875x;
m = 0, y = 0 — exact Taylor-coefficient checks (a_0 ≡ 0, a_1 = 0 at the corners) ⊕ interval enclosures, corner charts ρ = x/y,
ρ = x/u (u = π/16 − y); all region bounds rounded OUTWARD (s938 audit G22/G23).
PRIOR ART: §s935 (2) (the word, measured to N = 103), §s938 (2) (the q = 5 legality template ⊕ its audit), §s938 (3) (the identity).
DATA: none."""
import _setup  # noqa: F401  (deposit: this directory on the path, cwd = odd_polygons/)
import sys
import mpmath as mp
from mpmath import iv
import appA_2_identity as P
import q5_loops as K
from appA_2_identity import add, sh, cj, lin, mulc, neg, O, V, refl, comp, red
mp.mp.dps = 40
iv.dps = 30
D0 = (1, -8, -14)
def prefix(s1, s2, t, par):
    half = tuple(c // 2 for c in (t if par == 0 else lin(t, (0, 0, -1))))
    assert all(2 * h == c for h, c in zip(half, (t if par == 0 else lin(t, (0, 0, -1)))))
    Rm = P.rot_pow(s1, s2, half)
    return Rm if par == 0 else comp(Rm, refl(s1))
def fval(p):
    w = sh(add(p, V[0], -1), neg(D0))
    return red(sh(add(w, cj(w), -1), (3, 0, 0)))
def re2(p):
    w = sh(add(p, V[0], -1), neg(D0))
    return red(add(w, cj(w)))
def other_end(side, Pidx):
    a, b = side, (side + 1) % 4
    return V[b] if a == Pidx else V[a]
RUNS = [(1, 0, 1, (0, 2, 3)), (3, 2, 3, (0, 4, 6)), (0, 1, 1, (0, 2, 2))]    # (first, second, centre index, length)
def conditions():
    conds = []
    G = (0, O, {})
    for i, (s1, s2, pc, r) in enumerate(RUNS):
        c = fval(P.add(sh(cj(V[pc]) if G[0] else V[pc], G[1]), G[2]))
        def ray(t, par):
            w = s1 if par == 0 else s2
            F = comp(G, prefix(s1, s2, t, par))
            end = other_end(w, pc)
            return fval(add(sh(cj(end) if F[0] else end, F[1]), F[2]))
        rodd = r[2] % 2 == 1 and r[1] % 2 == 0
        rm1, rm2 = lin(r, (0, 0, -1)), lin(r, (0, 0, -2))
        conds.append((f'run{i+1} t=0', ray(O, 0), c))
        conds.append((f'run{i+1} t=1', ray((0, 0, 1), 1), c))
        conds.append((f'run{i+1} last-even', ray(rm1 if rodd else rm2, 0), c))
        conds.append((f'run{i+1} last-odd', ray(rm2 if rodd else rm1, 1), c))
        G = comp(G, prefix(s1, s2, r, 1 if rodd else 0))
    end = fval(add(sh(cj(V[0]) if G[0] else V[0], G[1]), G[2]))
    return conds, end, [('V1 ahead', re2(V[1])), ('V2 ahead', re2(V[2]))]
def evf(f, N, m):
    return sum(c * mp.cos(mp.pi * (e[0] * N + e[1] * m + e[2]) / (2 * N)) for e, c in f.items())

YT = mp.pi / 16 + mp.mpf(10) ** -35           # outward
XM = mp.pi / 30 + mp.mpf(10) ** -35           # outward (N = 15 sits at x = π/30)
Y0 = mp.pi / 400
def fold(f):
    out = {}
    for (A, B, C), c in f.items():
        e = (A, B, C) if B > 0 or (B == 0 and C >= 0) else (-A, -B, -C)
        out[e] = out.get(e, 0) + c
    return {e: c for e, c in out.items() if c}
def prep(f): return [(c, iv.pi * A / 2, B, iv.mpf(C)) for (A, B, C), c in fold(f).items()]
def amp_zero(f, n):              # exact: per β group Σ c δ^n i^{A+n} = 0 (Gaussian integers)
    g = {}
    for (A, B, C), c in fold(f).items():
        g.setdefault(B, []).append(c * C ** n * (1j) ** ((A + n) % 4))
    return all(abs(sum(v)) == 0 for v in g.values())
def a1_zero_at0(f): return sum(c * C * round(mp.sin(mp.pi * A / 2)) for (A, B, C), c in fold(f).items()) == 0
def a0_zero_top(f):
    import sympy as sp
    z = sp.symbols('z')                       # cos(Aπ/2 + Bπ/16): ζ_64 = e^{iπ/32}: exponent 16A + 2B
    e = sum(c * (z ** ((16 * A + 2 * B) % 64) + z ** ((-16 * A - 2 * B) % 64)) for (A, B, C), c in fold(f).items())
    return sp.rem(sp.expand(e), sp.cyclotomic_poly(64, z), z) == 0
def Fd(T, X, Y): return sum(c * iv.cos(th + b * Y + d * X) for c, th, b, d in T)
def FH(T, X, Y): return sum(c * (-iv.sin(th + b * Y) * d * K.Sc(d * X) - iv.cos(th + b * Y) * d * d * X * K.Vs(d * X)) for c, th, b, d in T)
def Fc0(T, R, Y):
    tot = iv.mpf(0)
    for c, th, b, d in T:
        st, ct = iv.sin(th), iv.cos(th)
        dsin = -st * b * b * Y * K.Vs(b * Y) + ct * b * K.Sc(b * Y)
        tot += c * (-d * (dsin * K.Sc(d * R * Y) + st * d * d * R * R * Y * K.Ws(d * R * Y)) - iv.cos(th + b * Y) * d * d * R * K.Vs(d * R * Y))
    return tot
def FcT(T, R, U):
    tot = iv.mpf(0)
    for c, th, b, d in T:
        phi = th + b * iv.pi / 16; lam = d * R - b
        tot += c * (-iv.cos(phi) * lam * lam * U * K.Vs(lam * U) - iv.sin(phi) * lam * K.Sc(lam * U))
    return tot
INR = lambda x0, x1, y0, y1: y1 >= x0 and y0 <= YT - mp.mpf('1.875') * x0
ANY = lambda *a: True
def f00_zero(f): return sum(c * round(mp.cos(mp.pi * A / 2)) for (A, B, C), c in fold(f).items()) == 0
def a1_zero_top(f):
    import sympy as sp
    z = sp.symbols('z')
    e = sum(c * C * (z ** ((16 * A + 2 * B) % 64) - z ** ((-16 * A - 2 * B) % 64)) for (A, B, C), c in fold(f).items())
    return sp.rem(sp.expand(e), sp.cyclotomic_poly(64, z), z) == 0
def Fx2(T, X, Y): return sum(c * (-iv.sin(th) * d * d * d * X * K.Ws(d * X) - iv.cos(th) * d * d * K.Vs(d * X)) for c, th, b, d in T)
def Fc0d(T, R, Y):            # f(ρy, y)/y, direct f with f(0,0) = 0
    tot = iv.mpf(0)
    for c, th, b, d in T:
        lam = b + d * R
        tot += c * (-iv.cos(th) * lam * lam * Y * K.Vs(lam * Y) - iv.sin(th) * lam * K.Sc(lam * Y))
    return tot
def FcTH(T, R, U):            # (f/x)/u near (0, π/16), ρ = x/u, when a_0 ≡ 0 and a_1(π/16) = 0
    tot = iv.mpf(0)
    for c, th, b, d in T:
        phi = th + b * iv.pi / 16; sp_, cp = iv.sin(phi), iv.cos(phi)
        dsin = -sp_ * b * b * U * K.Vs(b * U) - cp * b * K.Sc(b * U)          # (sin(φ − βu) − sin φ)/u
        tot += c * (-d * K.Sc(d * R * U) * dsin - d * sp_ * d * d * R * R * U * K.Ws(d * R * U) - iv.cos(phi - b * U) * d * d * R * K.Vs(d * R * U))
    return tot
def cert(g, sg, m0):
    T = prep(g); z0 = amp_zero(g, 0)
    base = (lambda X, Y: FH(T, X, Y)) if z0 else (lambda X, Y: Fd(T, X, Y))
    if m0:
        if f00_zero(g):
            if a1_zero_at0(g):
                return K.certify(lambda X, Y: Fx2(T, X, Y), (mp.mpf(0), XM, mp.mpf(0), mp.mpf(0)), ANY, sg)[0], 'm=0 f/x²'
            return K.certify(lambda X, Y: FH(T, X, Y), (mp.mpf(0), XM, mp.mpf(0), mp.mpf(0)), ANY, sg)[0], 'm=0 f/x'
        return K.certify(lambda X, Y: Fd(T, X, Y), (mp.mpf(0), XM, mp.mpf(0), mp.mpf(0)), ANY, sg)[0], 'm=0 f'
    parts, ok, ylo, yhi = [], True, mp.mpf(0), YT
    if z0 and a1_zero_at0(g):
        ok &= K.certify(lambda R, Y: Fc0(T, R, Y), (mp.mpf(0), mp.mpf(1) + mp.mpf(10) ** -30, mp.mpf(0), Y0), ANY, sg)[0]; parts.append('corner0'); ylo = Y0
    if (not z0) and f00_zero(g):
        ok &= K.certify(lambda R, Y: Fc0d(T, R, Y), (mp.mpf(0), mp.mpf(1) + mp.mpf(10) ** -30, mp.mpf(0), Y0), ANY, sg)[0]; parts.append('corner0d'); ylo = Y0
    if (not z0) and a0_zero_top(g):
        ok &= K.certify(lambda R, U: FcT(T, R, U), (mp.mpf(0), mp.mpf(8) / 15 + mp.mpf(10) ** -30, mp.mpf(0), Y0), ANY, sg)[0]; parts.append('cornerTop'); yhi = YT - Y0
    if z0 and a1_zero_top(g):
        ok &= K.certify(lambda R, U: FcTH(T, R, U), (mp.mpf(0), mp.mpf(8) / 15 + mp.mpf(10) ** -30, mp.mpf(0), Y0), ANY, sg)[0]; parts.append('cornerTopH'); yhi = YT - Y0
    ok &= K.certify(base, (mp.mpf(0), XM, ylo, yhi), INR, sg)[0]
    return ok, '+'.join(parts + ['f/x' if z0 else 'f'])


def span_ok():
    """the arc step needs each run's span < π (audit W13): run 1 (2m+2)·4π/N, run 2 (4m+5)·2π/N, run 3 (2m+1)·4π/N — i.e.
    8m+8 < N, 8m+10 < N, 8m+4 < N — all implied by N ≥ 8m + 15 (affine in m: check the slopes and m = 0)"""
    for slope, const in ((8, 8), (8, 10), (8, 4)):
        assert slope <= 8 and const < 15
    return True


def main():
    assert span_ok()
    conds, end, ahead = conditions()
    assert end == {}, 'end identity'
    allok = True
    for m0 in (False, True):
        for n, f, c in conds:
            cs = 1 if evf(c, 201, 0 if m0 else 12) > 0 else -1
            okc, _ = cert(c, cs, m0)
            okf, mf = cert(f, -cs, m0)
            allok &= okc and okf
            print('m=0' if m0 else 'm≥1', f'{n:16s}', 'OK' if okc and okf else 'FAIL', mf, flush=True)
        for n, g in ahead:
            ok, mm = cert(g, 1, m0)
            allok &= ok
            print('m=0' if m0 else 'm≥1', f'{n:16s}', 'OK' if ok else 'FAIL', mm)
    # control: the flipped sign must fail on every ray condition (m ≥ 1)
    fired = sum(not cert(f, 1 if evf(c, 201, 12) > 0 else -1, False)[0] for n, f, c in conds)
    print('control: flipped-sign ray conditions that fail to certify:', fired, '/', len(conds))
    assert allok and fired == len(conds)
    print('ALL Π_N LEGALITY CERTIFIED: W_m is realised for every odd N ≥ 8m + 15 — ℓ_m = c_{4m+5} + c_{4m+7} PROVED')


if __name__ == '__main__':
    main()
