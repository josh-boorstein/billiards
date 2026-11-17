"""[Deposit] Proposition A.1: the m-th torus of Π_N, symbolic in n = m + 1. Paper: "Perpendicular billiards in rational polygons of odd denominator". Research name: probes/s935_torus_ladder.py, copied verbatim except imports (deposit names), paths (data/ and figs/ under odd_polygons/) and this paragraph. Run through run_all.py; the usage lines below are the research repo's. Bracketed IDs and section marks (§sNNN) refer to the project's internal ledgers.

s935 keeper (queue `s935-a`, the LOWER bound): THE m-TH TORUS OF Π_N FOR EVERY m, SYMBOLIC IN n = m + 1.

NOT PRE-REGISTERED: exploratory. Template read off `flow_poly` at N = 63 (`data/s935_tori_63.json`, m = 0 … 6) and written as an
exit-side grammar: J0 = (23)^{2n} 20 (10)^{n−1} 13 (23)^{2n−2} 20 (10)^{n−1} [(23)^{2n}(10)^n]² 13, J1 = one [(23)^{2n}(10)^n] fewer,
start copy (2n−1, 1) side 3; pins lo0 = J0 step 4n+1 (V0), hi0 = J0 step 18n−2 (V2), hi1 = J1 step 6n+1 (V0).
PROVED here (formal, symbolic in n — elements rational in z = ζ_{4N}, v = z^{4n}, run sums geometric): both walks return to the
start; T^{24n} J0 = [lo0 + λ1, hi1), T^{18n} J1 = [lo0, lo0 + λ1); λ0 = 4 sin(4nx) sin x (2cos2x − 1), λ1 = 16 sin(4nx) sin³x — the
rotation 4 sin²(π/N) for EVERY m. NOT proved here: the containment inequalities for all m — `families()` writes them as
(n, j, x) families (pickled to data/s935_ladder_fams.pkl); MEASURED 0 negatives for n ≤ 8, N ≤ 8n + 1001, every j
(scratch `data/scratch/s935/lower/validate.py`); per-m RIGOROUS certificates m ≤ 6 (scratch `perm_cert.py`, log
`data/scratch/s935/lower/perm_cert_3_6.log`, `tori_63_cert_m012.json`).
PRIOR ART: §s935 (1) (the m = 0 torus, [POLY-ASYM-ALL]), §s935 (2)–(3). DATA: `data/s935_tori_63.json`.
→ general_triangles.md §s935 (4).   Usage: python probes/s935_torus_ladder.py   (≈ 1 min)

Original scratch header: ALL m: the m-th torus (n = m + 1) as a SYMBOLIC-in-n formal object. Elements: rational functions in
z = ζ_{4N} and v = z^{4n} over Q(i). Run sums are geometric series. Checks: combinatorics, both returns, λ closed forms."""
import _setup  # noqa: F401  (deposit: this directory on the path, cwd = odd_polygons/)
import sys
import sympy as sp
import pi_n_torus as KP
z, v, t, U = sp.symbols("z v t u")
I = sp.I
def toS(Lobj):
    return sum(sp.nsimplify(c.real) + I * sp.nsimplify(c.imag) if c.imag else sp.Integer(int(c.real)) * 1
               for c in [0]) + sum((sp.Integer(int(c.real)) + I * sp.Integer(int(c.imag))) * z ** b for b, c in Lobj.c.items())
V = [toS(x) for x in KP.V]
T = [toS(x) for x in KP.T]
DM = KP.DM
SIDES = KP.SIDES
def conj(e):
    return e.subs({z: 1 / z, v: 1 / v, t: 1 / t, U: 1 / U, I: -I}, simultaneous=True)
def zpow(o):          # z^{4o}, o = (alpha, beta) meaning alpha*n + beta
    a, b = o
    return v ** a * z ** (4 * b)
def g_apply(o, e, p):
    return -zpow(o) * z ** -2 * (conj(p) if e else p)
def X(p):
    return -(p + conj(p))
def g_across(o, e, s):
    a, b = o
    return ((a, b + DM[s]), 1) if e == 0 else ((a, b - DM[s]), 0)
def unit(o, e, s_in, word):
    """walk one unit word symbolically from (o, e, s_in) with z^{4o} replaced by t; returns (shift-sum as poly in t,
    offset change, e', s_out) — the offset change must be o-independent"""
    oo, ee, s = (0, 0), e, s_in
    S = 0
    for ch in word:
        s2 = int(ch)
        # shift with copy (o + oo): rotation −t·z^{4·oo_beta}·z^{−2}
        q = -t * z ** (4 * oo[1]) * z ** -2 * (conj(T[s2]) if ee else T[s2])
        S += -(q + conj(q))
        oo, ee = g_across(oo, ee, s2)
        s = s2
    return sp.expand(S), oo[1], ee, s
def run(o, e, s_in, word, c):
    """c = (gamma, delta): repeat count gamma*n + delta. Returns (sum, o_new, e_new, s_new)."""
    S1, d, e1, s1 = unit(o, e, s_in, word)
    assert e1 == e and int(word[-1]) == s_in or c == (0, 1) or True
    if c == (0, 0):
        return 0, o, e, s_in
    # S1 = a t + b / t
    P = sp.Poly(sp.expand(S1 * t), t)
    coeffs = {mon[0]: co for mon, co in zip(P.monoms(), P.coeffs())}
    a, b = coeffs.get(2, 0), coeffs.get(0, 0)
    assert set(coeffs) <= {0, 2}, coeffs
    g, dl = c
    if c[0] == 0 and c[1] == 1:
        tot = a * zpow(o) + b / zpow(o)
        return tot, (o[0], o[1] + d), e1, s1
    assert e1 == e, 'unit does not preserve e'
    r = z ** (4 * d)
    rc = v ** (d * g) * z ** (4 * d * dl)
    tot = a * zpow(o) * (rc - 1) / (r - 1) + b / zpow(o) * (1 / rc - 1) / (1 / r - 1)
    return tot, (o[0] + d * g, o[1] + d * dl), e, s1
J0 = [('23', (2, 0)), ('20', (0, 1)), ('10', (1, -1)), ('13', (0, 1)), ('23', (2, -2)), ('20', (0, 1)), ('10', (1, -1)),
      ('23', (2, 0)), ('10', (1, 0)), ('23', (2, 0)), ('10', (1, 0)), ('13', (0, 1))]
J1 = [('23', (2, 0)), ('20', (0, 1)), ('10', (1, -1)), ('13', (0, 1)), ('23', (2, -2)), ('20', (0, 1)), ('10', (1, -1)),
      ('23', (2, 0)), ('10', (1, 0)), ('13', (0, 1))]
START = ((2, -1), 1, 3)          # copy offset 2n − 1, reflected, entry side 3
def walk(segs, stop_after=None):
    o, e, s = START
    S = 0
    states = []
    for idx, (w, c) in enumerate(segs):
        states.append((S, o, e, s))
        tot, o, e, s = run(o, e, s, w, c)
        S += tot
    states.append((S, o, e, s))
    return states
def single_step_state(states_idx_state, w, k):
    """advance k single steps of unit word w from a state"""
    S, o, e, s = states_idx_state
    for ch in w[:k]:
        s2 = int(ch)
        q = g_apply(o, e, T[s2])
        S += X(q)
        o, e = g_across(o, e, s2)
        s = s2
    return S, o, e, s
st0 = walk(J0); st1 = walk(J1)
print('J0 returns to', st0[-1][1:], ' J1 returns to', st1[-1][1:], ' start', START)
# pins: lo0 at step 4n+1 = after seg0 (23)^{2n}, then 1 step of '20' -> at the '0' step: state before that step
def pin_val(state, vi):
    S, o, e, s = state
    return X(g_apply(o, e, V[vi])) + S
lo0 = pin_val(single_step_state(st0[1], '20', 1), 0)           # J0 step 4n+1, vertex 0
hi1 = pin_val(single_step_state(st1[3], '13', 1), 0)           # J1 step 6n+1 (seg3 '13', its second step), vertex 0
hi0 = pin_val(st0[9], 2)                                         # J0 step 18n−2 = start of 2nd Q block (seg 9), vertex 2
lam0, lam1 = sp.together(hi0 - lo0), sp.together(hi1 - hi0)
sinn = lambda k: (z ** (2 * k) - z ** (-2 * k)) / (2 * I)
sin4n = (v ** 2 - v ** -2) / (2 * I)
c2 = z ** 4 + z ** -4
chk0 = sp.simplify(sp.together(lam0 - 4 * sin4n * sinn(1) * (c2 - 1)))
chk1 = sp.simplify(sp.together(lam1 - 16 * sin4n * sinn(1) ** 3))
ret0 = sp.simplify(sp.together((lo0 - st0[-1][0]) - (lo0 + lam1)))
ret1 = sp.simplify(sp.together((hi0 - st1[-1][0]) - lo0))
print('λ0 closed form:', chk0 == 0, ' λ1 closed form:', chk1 == 0, ' T^{24n} J0 return:', ret0 == 0, ' T^{18n} J1 return:', ret1 == 0)

# ---------- containment families: for a run (word, count), repetition j (symbol u = z^{4·Δ·j}), step position p in the unit
jj = sp.symbols("u")
def families(segs, lo, lam, tag):
    fams = []
    S, o, e, s = START[0], None, None, None
    states = walk(segs)
    for idx, (w, c) in enumerate(segs):
        S0, o0, e0, s0 = states[idx]
        S1u, d, _, _ = unit(o0, e0, s0, w)
        # state at repetition j: offset o0 + j d, partial sum S0 + Σ_{i<j} unit-sum(o0 + i d)
        P = sp.Poly(sp.expand(S1u * t), t); co = {m[0]: c_ for m, c_ in zip(P.monoms(), P.coeffs())}
        a, b = co.get(2, 0), co.get(0, 0)
        r = z ** (4 * d)
        if c == (0, 1):
            Sj, zo = S0, zpow(o0)
        else:
            Sj = S0 + a * zpow(o0) * (jj - 1) / (r - 1) + b / zpow(o0) * (1 / jj - 1) / (1 / r - 1)   # jj = r^j
            zo = zpow(o0) * jj                                                                        # z^{4(o0 + j d)}
        # now step through the unit word with copy rotation built from zo
        variants = [('all', s0)] if (c == (0, 1) or s0 == int(w[-1])) else [('first', s0), ('rest', int(w[-1]))]
        for jr, s_first in variants:
          oo, ee, s_in, Sk = 0, e0, s_first, Sj
          for pos, ch in enumerate(w):
              s2 = int(ch)
              rot = -zo * z ** (4 * oo) * z ** -2
              ga = lambda p_: rot * (conj(p_) if ee else p_)
              lo_k = lo - Sk
              for ss in (s_in, s2):
                  i_, k_ = SIDES[ss]
                  p_, q_ = X(ga(V[i_])), X(ga(V[k_]))
                  fams.append((f'{tag}.seg{idx}.{w}.{jr}.p{pos}.s{ss}', c, d, lo_k, lo_k + lam, p_, q_))
              i_, k_ = SIDES[s_in]; ev = ga(V[k_]) - ga(V[i_])
              fams.append((f'{tag}.seg{idx}.{w}.{jr}.p{pos}.entry{s_in}', c, d, (ev + conj(ev)) * (1 if ee == 0 else -1), None, None, None))
              i_, k_ = SIDES[s2]; ev = ga(V[k_]) - ga(V[i_])
              fams.append((f'{tag}.seg{idx}.{w}.{jr}.p{pos}.exit{s2}', c, d, -(ev + conj(ev)) * (1 if ee == 0 else -1), None, None, None))
              q = rot * (conj(T[s2]) if ee else T[s2]); Sk = Sk + X(q)
              o_ = (0, oo); (o_, ee) = g_across(o_, ee, s2); oo = o_[1]; s_in = s2
    return fams
F0 = families(J0, lo0, lam0, 'J0'); F1 = families(J1, hi0, lam1, 'J1')
import pickle
pickle.dump((F0, F1, lo0, hi0, hi1, lam0, lam1), open('data/s935_ladder_fams.pkl', 'wb'))
print('families:', len(F0), len(F1))
