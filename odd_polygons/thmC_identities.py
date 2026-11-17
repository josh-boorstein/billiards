"""[Deposit] Theorem 5.4, step (i) and (ii), Appendix B.0: the template and the identities. Paper: "Perpendicular billiards in rational polygons of odd denominator". Research name: probes/s951_series_first_look.py, copied verbatim except imports (deposit names), paths (data/ and figs/ under odd_polygons/) and this paragraph. Run through run_all.py; the usage lines below are the research repo's. Bracketed IDs and section marks (§sNNN) refer to the project's internal ledgers.

s951 keeper (queue `s946-b`, first look): DOES Π_N's ALL-N PROOF GO SYMBOLIC IN m FOR [POLY-PI-SERIES]?
S₁(m) = (N−3m−1, 4m, N−m−1, 2)/N with the isosceles apex (census `quad(N, N−6m, 4m, N−2m, N−2, 2, N−2)`): YES — S₂(m) = (N−2m−2, 4m, N−2m, 2)/N: NOT BY THIS ROUTE.

NOT PRE-REGISTERED: exploratory (a first look). Found by reading flow_poly's non-periodic components of the series rows:
  S₁(m): an equal pair of TORI, rotation number 4 sin²(mπ/N) (exact match m = 1 … 5), return times (16m+8, 12m+6); the return
    words (data/scratch/s951/s1_itin.json, m = 1, 2, 3, 5, 6; m = 4 read at another cut of the same torus) are ONE template,
    exit sides from entry 3 of copy o₀ = (3m−1)/2, e = 1:
      J0 = (23)^{3m−1} 2 · 013 · (23)^{m−1} 2 · 02 · (32)^{2m−1} 3 · 102 · (32)^{2m−1} 3 · 1013          (16m + 8 strips)
      J1 = (23)^{3m−1} 2 · 013 · (23)^{m−1} 2 · 02 · (32)^{2m−1} 3 · 1013                                (12m + 6 strips)
    pins lo0 = J0 step 6m−1 (V0), hi0 = J0 step 12m+4 (V2), hi1 = J1 step 6m+1 (V0). At m = 1 this IS s935's Π_N torus
    (`s935_pi_n_torus.J0`, exit for exit): S₁(m) is Π_N's torus stretched, NOT a rung of Π_N's ladder (§s935 (4)).
  S₂(m): the non-periodic components are 2m-interval exchanges when certified (Boshernitzan), permutation of genus m with
    ONE cone point (minimal stratum), length module of full rank 2m — genus-m minimal components, no rotation template.
PROVED here (formal, symbolic in m; Laurent in z = ζ_{4N} = e^{ix/2}, w = z^{2m} = e^{imx}; runs as geometric sums):
  S₁(m) closes; both walks return to the start copy; T^{16m+8} J0 = [lo0 + λ1, hi1), T^{12m+6} J1 = [lo0, lo0 + λ1);
  λ0 = 4 sinx sin2mx sin4mx (2cos2mx − 1), λ1 = 16 sinx sin²mx sin2mx sin4mx (sides = law of sines / cos x, a positive scale)
  ⟹ rotation number 4 sin²(mπ/N), irrational for N ≥ 14m + 1 (2cos(2πm/N) has reduced denominator ≥ 14 > 6).
  Even m: o₀ is a half-integer — the formal copy is −(a real copy), i.e. the unfolding of the REVERSED flow, which is the same
  perpendicular class (the reflection in V0V1 reverses ω); the N/2-shifted start fails every row (asserted below).
MEASURED here: the containment inequalities (each strip inside its entry and exit side's X-range, sides entry/exit) at
  m = 1 … 8, N ∈ {onset, onset+2, onset+30, ≈4·onset}: 0 violations; CONTROL N = 14m − 1: 6 violations every m ≥ 2 (19 at m = 1,
  = §s935's Π₁₃ count) — the template's own onset is EXACTLY §s946 (3)'s pre-registered 14m + 1.
THE RIGOROUS all-(m, N) inequality certificate: `probes/s951_s1_allm_cert.py` (reads this file's `families()` output).
  NOT DONE: S₂(m) for all m (needs a genus-m certificate, queue `s951-b`).
COULD HAVE COME OUT OTHERWISE: the identities could hold only at the measured m (they hold formally in w); the template's
  inequalities could fail at some m or onset could differ from 14m + 1 (they do not); S₂'s components could have been tori
  (2-label snapshots, as for S₁) — they are not.
PRIOR ART: §s946 (3) ([POLY-PI-SERIES], MEASURED; onset law 54/54); §s935 (1), (4) (Π_N torus, ladder template), §s937 (the
  all-m certificate method), §s941; `rulings.py --grep` series, torus → nothing on these series. Closed-form grep `4sin²(mπ/N)` →
  only Π_N's `4 sin²(π/N)` (its m = 1 case). DATA: data/s946_apex_census.json (series rows); `data_index.py` → none on S₁/S₂ tori.
→ general_triangles.md §s951.   Usage: PYTHONPATH=… python probes/s951_series_first_look.py [--no-s2]   (≈ 2 min)
"""
import _setup  # noqa: F401  (deposit: this directory on the path, cwd = odd_polygons/)
import cmath
import math
import sys
from fractions import Fraction as Fr

import sympy as sp

# ---------------------------------------------------------------- formal frame, symbolic in m
z, w = sp.symbols('z w')
I = sp.I


def conj(e):
    return e.subs({z: 1 / z, w: 1 / w, I: -I}, simultaneous=True)


def X(p):
    return -(p + conj(p))


def zo(o):                                     # z^{4o}, o = (α, β) meaning αm + β
    a, b = Fr(o[0]), Fr(o[1])
    assert (2 * a).denominator == 1 and (4 * b).denominator == 1
    return w ** int(2 * a) * z ** int(4 * b)


sinx = (z ** 2 - z ** -2) / (2 * I)
LEN = [2 * sinx * (w + 1 / w) / 2, 2 * sinx * (w ** 3 + w ** -3) / 2, (w ** 4 - w ** -4) / (2 * I), (w ** 4 - w ** -4) / (2 * I)]
DIR = [sp.Integer(1), -w ** -4, -w ** -3 * z ** 2, w ** -3 * z ** -2]      # z^{2 m_full}, m_full = 0, N−4m, N−3m+1, 2N−3m−1
REF = [sp.Integer(1), w ** -8, w ** -6 * z ** 4, w ** -6 * z ** -4]       # z^{4 m_full}
DM = [(0, 0), (-4, 0), (-3, 1), (-3, -1)]                                 # m_full mod N, as (α, β)
SIDES = [(s, (s + 1) % 4) for s in range(4)]
V = [sp.Integer(0)]
for _s in range(3):
    V.append(V[-1] + LEN[_s] * DIR[_s])
T = [V[i] - REF[s] * conj(V[i]) for s, (i, _) in enumerate(SIDES)]


def rot(o):
    return -zo(o) * z ** -2


def across(o, e, s):
    (a, b), (da, db) = o, DM[s]
    return ((a + da, b + db), 1) if e == 0 else ((a - da, b - db), 0)


def unitA(e, word):
    oo, ee, A = (Fr(0), Fr(0)), e, 0
    for ch in word:
        s2 = int(ch)
        A += zo(oo) * (conj(T[s2]) if ee else T[s2])
        oo, ee = across(oo, ee, s2)
    return A, oo, ee


def run(state, word, c):
    S, o, e = state
    A, d, e1 = unitA(e, word)
    if c == (0, 1):
        return S + X(rot(o) * A), (o[0] + d[0], o[1] + d[1]), e1
    g, dl = c
    assert e1 == e and d[0] == 0, (word, d)
    r, rc = z ** int(4 * d[1]), w ** int(2 * d[1] * g) * z ** int(4 * d[1] * dl)
    return S + X(rot(o) * A * (rc - 1) / (r - 1)), (o[0] + d[1] * g, o[1] + d[1] * dl), e


START = (sp.Integer(0), (Fr(3, 2), Fr(-1, 2)), 1)
J0 = [('23', (3, -1)), ('2', (0, 1)), ('013', (0, 1)), ('23', (1, -1)), ('2', (0, 1)), ('02', (0, 1)), ('32', (2, -1)),
      ('3', (0, 1)), ('10', (0, 1)), ('2', (0, 1)), ('32', (2, -1)), ('3', (0, 1)), ('1013', (0, 1))]
J1 = [('23', (3, -1)), ('2', (0, 1)), ('01', (0, 1)), ('3', (0, 1)), ('23', (1, -1)), ('2', (0, 1)), ('02', (0, 1)),
      ('32', (2, -1)), ('3', (0, 1)), ('1013', (0, 1))]


def walk(segs):
    st, out = START, [START]
    for wd, c in segs:
        st = run(st, wd, c)
        out.append(st)
    return out


def pin(st, vi):
    S, o, e = st
    return X(rot(o) * (conj(V[vi]) if e else V[vi])) + S


def formal():
    assert sp.simplify(sp.expand(V[-1] + LEN[3] * DIR[3])) == 0, 'S1(m) does not close formally'
    s0, s1 = walk(J0), walk(J1)
    assert s0[-1][1:] == START[1:] and s1[-1][1:] == START[1:], 'a walk does not return to the start copy'
    lo0, hi0, hi1 = pin(s0[2], 0), pin(s0[9], 2), pin(s1[3], 0)
    lam0, lam1 = sp.together(hi0 - lo0), sp.together(hi1 - hi0)
    assert sp.simplify(sp.together(-s0[-1][0] - lam1)) == 0, 'T^{16m+8} J0 is not [lo0 + λ1, hi1)'
    assert sp.simplify(sp.together(s1[-1][0] - lam0)) == 0, 'T^{12m+6} J1 is not [lo0, lo0 + λ1)'
    sn = lambda a, b: (z ** (2 * a) * w ** b - z ** (-2 * a) * w ** (-b)) / (2 * I)   # sin(ax + bmx), w = e^{imx}
    c2m = w ** 2 + w ** -2
    # with sides scaled by 1/cos x: λ = (law-of-sines λ)/cos x
    L0 = 4 * sn(1, 0) * sn(0, 2) * sn(0, 4) * (c2m - 1)
    L1 = 16 * sn(1, 0) * sn(0, 1) ** 2 * sn(0, 2) * sn(0, 4)
    assert sp.simplify(sp.together(lam0 - L0)) == 0, 'λ0 closed form'
    assert sp.simplify(sp.together(lam1 - L1)) == 0, 'λ1 closed form'
    return 'closure, returns, T^{16m+8}J0, T^{12m+6}J1, λ0 = 4 sinx sin2mx sin4mx (2cos2mx−1), λ1 = 16 sinx sin²mx sin2mx sin4mx (sides/cos x): FORMAL in w'


FAMS_OUT = 'data/s951_s1_fams.pkl'


def families():
    """the containment inequalities as symbolic families, in `s937_ladder_allm_cert`'s input format: (name, c, d, lo_k, lo_k + λ,
    p, q) or (name, c, d, sign-expression, None, None, None); Laurent in z, v = w = e^{imx}, u = z^{4dj} (j = repetition inside
    a run of count c = (γ, δ): γm + δ). Consumed by `probes/s951_s1_allm_cert.py`."""
    import pickle
    u, v = sp.symbols('u v')
    cj = lambda e: e.subs({z: 1 / z, w: 1 / w, u: 1 / u, I: -I}, simultaneous=True)
    Xu = lambda p: -(p + cj(p))

    def fam(segs, lo, lam, tag):
        out, sts = [], walk(segs)
        for idx, (wd, c) in enumerate(segs):
            S0, o0, e0 = sts[idx]
            A, d, _ = unitA(e0, wd)
            if c == (0, 1):
                Sk, rot0, dd = S0, rot(o0), 0
            else:
                dd = int(d[1])
                Sk = S0 + Xu(rot(o0) * A * (u - 1) / (z ** (4 * dd) - 1))      # Σ_{i<j} unit sums, u = z^{4dd·j}
                rot0 = rot(o0) * u
            oo, ee = (Fr(0), Fr(0)), e0
            s_in = int(segs[idx - 1][0][-1]) if idx else 3
            for pos, ch in enumerate(wd):
                s2 = int(ch)
                R = rot0 * zo(oo)
                ga = lambda p_: R * (cj(p_) if ee else p_)
                lo_k = lo - Sk
                for ss in (s_in, s2):
                    i_, k_ = SIDES[ss]
                    out.append((f'{tag}.seg{idx}.{wd}.all.p{pos}.s{ss}', c, dd, lo_k, lo_k + lam, Xu(ga(V[i_])), Xu(ga(V[k_]))))
                sgn = 1 if ee == 0 else -1
                i_, k_ = SIDES[s_in]
                ev = ga(V[k_]) - ga(V[i_])
                out.append((f'{tag}.seg{idx}.{wd}.all.p{pos}.entry{s_in}', c, dd, (ev + cj(ev)) * sgn, None, None, None))
                i_, k_ = SIDES[s2]
                ev = ga(V[k_]) - ga(V[i_])
                out.append((f'{tag}.seg{idx}.{wd}.all.p{pos}.exit{s2}', c, dd, -(ev + cj(ev)) * sgn, None, None, None))
                Sk = Sk + Xu(R * (cj(T[s2]) if ee else T[s2]))
                oo, ee = across(oo, ee, s2)
                s_in = s2
        return out
    s0, s1 = walk(J0), walk(J1)
    lo0, hi0, hi1 = pin(s0[2], 0), pin(s0[9], 2), pin(s1[3], 0)
    F0, F1 = fam(J0, lo0, hi0 - lo0, 'J0'), fam(J1, hi0, hi1 - hi0, 'J1')
    sw = lambda f: tuple(sp.expand(x).subs(w, v) if isinstance(x, sp.Basic) else x for x in f)
    F0, F1 = [sw(f) for f in F0], [sw(f) for f in F1]
    pickle.dump((F0, F1), open(FAMS_OUT, 'wb'))
    return len(F0), len(F1)


# ---------------------------------------------------------------- numeric containment at concrete (m, N)
def expand(m, which):
    head = '23' * (3 * m - 1) + '2' + '013' + '23' * (m - 1) + '2' + '02' + '32' * (2 * m - 1) + '3'
    return head + ('102' + '32' * (2 * m - 1) + '3' if which == 0 else '') + '1013'


def check(m, N, shift=0.0, tol=1e-12):
    x = math.pi / N
    sub = {z: cmath.exp(1j * x / 2), w: cmath.exp(1j * m * x)}
    Vn = [complex(sp.N(vv.subs(sub))) for vv in V]
    Tn = [complex(sp.N(t.subs(sub))) for t in T]
    zz = sub[z]
    R = lambda o: -zz ** (4 * o - 2)
    ga = lambda o, e, p: R(o) * (p.conjugate() if e else p)
    Xn = lambda p: -2 * p.real
    DMn = [a * m + b for a, b in DM]
    acr = lambda o, e, s: (o + DMn[s], 1) if e == 0 else (o - DMn[s], 0)

    def wk(word):
        o, e, S, s_in, out = (3 * m - 1) / 2 + shift, 1, 0.0, 3, []
        for ch in word:
            s2 = int(ch)
            out.append((o, e, s_in, s2, S))
            S += Xn(ga(o, e, Tn[s2]))
            o, e = acr(o, e, s2)
            s_in = s2
        return out
    w0, w1 = wk(expand(m, 0)), wk(expand(m, 1))
    pv = lambda st, vi: Xn(ga(st[0], st[1], Vn[vi])) + st[4]
    lo0, hi0, hi1 = pv(w0[6 * m - 1], 0), pv(w0[12 * m + 4], 2), pv(w1[6 * m + 1], 0)
    viol = 0
    for wkk, lo, lam in ((w0, lo0, hi0 - lo0), (w1, hi0, hi1 - hi0)):
        for (o, e, s, s2, S) in wkk:
            a, b = lo - S, lo - S + lam
            for ss in (s, s2):
                i, k = SIDES[ss]
                p, q = sorted((Xn(ga(o, e, Vn[i])), Xn(ga(o, e, Vn[k]))))
                viol += min(a - p, q - b) < -tol
            for ss, want in ((s, 1), (s2, -1)):
                i, k = SIDES[ss]
                ev = ga(o, e, Vn[k]) - ga(o, e, Vn[i])
                viol += 2 * ev.real * (1 if e == 0 else -1) * want < -tol
    return viol


# ---------------------------------------------------------------- S₂: genus of the certified components
def s2_genus(rows):
    import flow_decomp as FD
    import flow_poly as FP
    import diagonal_census as DC
    from pi_n_torus import PosDec

    class Track(PosDec):
        def _reduce(self, c):
            self.LAST[id(c)] = (list(c.top), list(c.bot))
            return super()._reduce(c)

    def genus(top, bot):
        d = len(top)
        pi = {i + 1: bot.index(top[i]) + 1 for i in range(d)}
        pi[0], pi[d + 1] = 0, d + 1
        inv = {v_: k_ for k_, v_ in pi.items()}
        sig = {j: inv[pi[j] + 1] - 1 for j in range(d + 1)}
        seen, cyc = set(), 0
        for j in range(d + 1):
            if j not in seen:
                cyc += 1
                k_ = j
                while k_ not in seen:
                    seen.add(k_)
                    k_ = sig[k_]
        return (d - cyc + 1) // 2, cyc
    out = []
    for m, N in rows:
        DC.N_CUR = N
        P = DC.quad(N, N - 4 * m - 2, 4 * m, N - 4 * m + 2, N - 2, 2, N - 2)
        pieces, lengths, top, bot, meta = FP._build_iet_poly(P, FP.perp_dirs_poly(P)[0])
        Dc = Track(FD._fast(P), lengths)
        Dc.LAST = {}
        comps = Dc.run(top, bot, 50000, meta['L'])
        gs = [(len(Dc.LAST[id(c)][0]), genus(*Dc.LAST[id(c)])) for c in comps if c.kind == 'no_periodic']
        out.append((m, N, gs))
        assert len(gs) == 2 and all(d == 2 * m and g == (m, 1) for d, g in gs), ('S2 genus', m, N, gs)
    return out


def polygon_check(rows):
    """the formal S1(m) polygon is the census polygon `s933_diagonal_census.quad` (the one flow_poly decomposed), up to scale"""
    import diagonal_census as DC
    out = []
    for m, N in rows:
        DC.N_CUR = N
        P = DC.quad(N, N - 6 * m, 4 * m, N - 2 * m, N - 2, 2, N - 2)
        cen = [sum(complex(c) * cmath.exp(2j * math.pi * e / P.N) for e, c in v_.items()) for v_ in P.V]
        x = math.pi / N
        sub = {z: cmath.exp(1j * x / 2), w: cmath.exp(1j * m * x)}
        frm = [complex(sp.N(vv.subs(sub))) for vv in V]
        dev = max(abs(a / cen[1] - b / frm[1]) for a, b in zip(cen, frm))
        assert dev < 1e-12 and (cen[1] / frm[1]).real > 0 and abs((cen[1] / frm[1]).imag) < 1e-12, ('polygon mismatch', m, N, dev)
        out.append((m, N, dev))
    return out


if __name__ == '__main__':
    print('POLYGON = census quad, up to positive scale:', polygon_check([(1, 17), (2, 31), (3, 45), (6, 87)]), flush=True)
    print('FORMAL:', formal(), flush=True)
    print('FAMILIES (J0, J1) →', FAMS_OUT, families(), flush=True)
    for m in range(1, 9):
        on = 14 * m + 1
        res = {N: check(m, N) for N in (on, on + 2, on + 30, 4 * on + 1)}
        ctrl = check(m, on - 2)
        shifted = check(m, on, shift=on / 2)
        assert all(v_ == 0 for v_ in res.values()) and ctrl > 0 and shifted > 0, (m, res, ctrl, shifted)
        print(f'm={m}: violations {res}  control N={on - 2}: {ctrl}  N/2-shifted start: {shifted}', flush=True)
    if '--no-s2' not in sys.argv:
        print('S2 genus (m, N, [(labels, (genus, cone points))]):', s2_genus([(2, 31), (2, 35), (3, 47), (4, 63)]))
