"""[Deposit] Proposition A.3: the frontier of the tori. Paper: "Perpendicular billiards in rational polygons of odd denominator". Research name: probes/s941_chain_frontier.py, copied verbatim except imports (deposit names), paths (data/ and figs/ under odd_polygons/) and this paragraph. Run through run_all.py; the usage lines below are the research repo's. Bracketed IDs and section marks (§sNNN) refer to the project's internal ledgers.

s941 keeper (queue `s940-a`): PROOF OF THE CHAIN LAW's FRONTIER HALF — for every odd N ≥ 15 and every m with 8m + 15 ≤ N,
the frontier of Π_N's m-th torus T_m (§s937, [POLY-PI-LOWER]) is EXACTLY slit m − 1 ⊕ slit m (slit 0 alone for m = 0):
  curve A (T_m on the +X side): c_{4n+3} · ℓ_{m−1} · c_{4n+1}      curve B (T_m on the −X side): c_{4n−1} · ℓ_m · c_{4n−3}
(n = m + 1; chord c_b = 2R sin(bπ/N), b strips V0↔V2; loop ℓ_m = c_{4m+5} + c_{4m+7}, 8m + 12 strips V0→V0), and slit m is
the SAME saddle connections in T_m's and T_{m+1}'s frontiers, on opposite sides. → general_triangles.md §s941.

NOT PRE-REGISTERED: a proof (formal identities ⊕ topology on top of the s937 certificate), not a hypothesis test.
SET-UP (§s935 (4) template; s935_torus_ladder grammar): R_m = closure of the saturation of J = [lo0, hi1) = the union of the
42n strips (vertical slabs [a, b] in convex copies; disjoint — J0, J1 are FIRST returns, s935). Every slab edge lies on one of
four edge-trajectories (lo0, hi0⁻, hi0⁺, hi1 over one pass); they form two paths A = J1.lo ⊕ J0.lo, B = J0.hi ⊕ J1.hi.
STEP 1 (touches). A path edge meets a vertex iff a containment condition of s937 (n ≥ 8) / §s935 (4)'s per-m certificates (n ≤ 7) is
ZERO. Those zeros are exactly 16 per n ≥ 2 (8 touches × entry/exit; s937's 12 formal ⊕ 4 attained, mapped one-to-one by audit s941
row 1; 18 at n = 1; every other condition certified strictly positive), and here
the 8 touches are located for ALL n ≥ 2 (symbolic, formal in z = ζ_{4N}, v = z^{4n}): on A at J1 position 18n−1 (V2), J0
positions 4n+2 (V0, pin lo0), 12n−2 (V0), 16n−1 (V2); on B at J0 18n−1 (V2, pin hi0), 22n−2 (V0), J1 6n+2 (V0, pin hi1),
10n−1 (V2). [check_touches] — and per n ≤ NMAX in exact Laurent arithmetic, no other formal touch [per_n].
STEP 2 (the glued parts). Before the first touch both paths run along hi0's trajectory (J0's slab on one side, J1's on the
other); after the last touch both arrive at lo0 + λ1 (T^{24n} lo0 = T^{18n} hi1 — s935's return identity), J0.lo and J1.hi
the same leaf with J0 on its +X and J1 on its −X side. Both parts are interior to R_m.
STEP 3 (the six pieces). Between touches: A1, A2, A3 on A, B1, B2, B3 on B; each is a saddle connection (vertex to vertex, no
vertex inside). Their developed holonomies equal ±c_{4n+3}, ±ℓ_{m−1}, ±c_{4n+1}, ±c_{4n−1}, ±ℓ_m, ±c_{4n−3} FORMALLY
(all n ≥ 2, [check_holonomy]); each is identified by its start corner (copy (o, e), vertex) — the forward vertical ray from a
corner into its copy is unique — with o linear in n [ids].
STEP 4 (exposure). R_m is flow-invariant, so a piece's exposure is constant between touches.
 (a) A1 ∋ lo0 and B2 ∋ hi1, and R_m meets J's line only in J (first return), so the sides of lo0 / hi1 facing away from J are
     outside R_m: A1, B2 are frontier.
 (b) A2, B1, B3 have the start corners of T_{m−1}'s B2, A1, A3 ([ids]: same linear forms at n − 1; for m = 1 read off the exact n = 1 data), which are T_{m−1}-slab edges on
     the OPPOSITE side; distinct minimal components have disjoint interiors (s937 distinctness), so these are frontier of T_m (m ≥ 1).
 (c) A3: Stokes on R_m — ∫_{∂R_m} dy = 0 (dy closed on the translation surface, ∂R_m vertical), i.e. Σ_{+X side} len =
     Σ_{−X side} len over the frontier. With A1, A2, B1, B2, B3 frontier and A3 not: c_{4n+3} + ℓ_{m−1} − c_{4n−1} − ℓ_m − c_{4n−3}
     = −c_{4n+1} ≠ 0 (ℓ_k = c_{4k+5} + c_{4k+7}); so A3 is frontier.
 m = 0 (n = 1): exact per-n data (9 touches, 18 zero conditions) — A: c7 (∋ lo0), c1, c3, c5; B: c3, ℓ_0 (∋ hi1), c1. A's and B's
 c1, c3 start at the SAME corners ((−1,1)V0, (3,0)V2): the same saddle connections with R_0 on both sides — interior. Stokes:
 c7 + [c5] = ℓ_0 = c5 + c7 ⟹ c5 frontier. (Corrected by audit s941 row 20; the earlier "B = ℓ_0" was wrong.) [n1_case]
CONCLUSION: frontier(T_m) = the six (three for m = 0); slit m ⊂ frontier(T_m) ∩ frontier(T_{m+1}) on opposite sides; so
α(T_m) = s_{m−1} − s_m and the chain T_0 – … – T_M carries α = −s_M; by the deck map h the same for h(T_m).
NOT PROVED HERE: that slit M faces invariant cylinders, and that the 2(M + 1) tori are ALL the minimal components (= `s935-a`).
PRIOR ART: §s940 (3) (the measured law, `probes/s940_pin_frontier.py`), §s937 ([POLY-PI-LOWER] certificate, its 16 zeros), §s935 (4)
(template, pins), §s938 (3) ([POLY-PI-LOOP]); rulings --grep frontier, slit → [NCYL-632], [NCYL-631] (slit tori, other polygons).
DATA: data/s940_pin_frontier.json (the measured frontiers it reproduces). Output data/s941_chain_frontier.json.
USAGE: PYTHONPATH=… python probes/s941_chain_frontier.py [NMAX]   (symbolic part ≈ 1–2 min)
"""
import _setup  # noqa: F401  (deposit: this directory on the path, cwd = odd_polygons/)
import sys, json, math
import sympy as sp

sys.path[:0] = ['probes']   # s958: frontier_template promoted to probes/s941_frontier_template.py
import pi_n_torus as KP

z, v, t = sp.symbols('z v t')
I = sp.I


def toS(Lobj):
    return sp.sympify(sum((sp.Integer(int(c.real)) + I * sp.Integer(int(c.imag))) * z ** b for b, c in Lobj.c.items()))


VS = [toS(x) for x in KP.V]
TS = [toS(x) for x in KP.T]
DM = KP.DM


def conj(e):
    return e.subs({z: 1 / z, v: 1 / v, t: 1 / t, I: -I}, simultaneous=True)


def zpow(o):                      # o = (alpha, beta): z^{4(alpha n + beta)}
    return v ** o[0] * z ** (4 * o[1])


def gap(o, e, p):
    return -zpow(o) * z ** -2 * (conj(p) if e else p)


def across(o, e, s):
    a, b = o
    return ((a, b + DM[s]), 1) if e == 0 else ((a, b - DM[s]), 0)


def Xs(p):
    return -(p + conj(p))


# grammar (s935_torus_ladder): (unit word of exit sides, repeat count gamma*n + delta)
G0 = [('23', (2, 0)), ('20', (0, 1)), ('10', (1, -1)), ('13', (0, 1)), ('23', (2, -2)), ('20', (0, 1)), ('10', (1, -1)),
      ('23', (2, 0)), ('10', (1, 0)), ('23', (2, 0)), ('10', (1, 0)), ('13', (0, 1))]
G1 = G0[:9] + G0[11:]
START = ((2, -1), 1, 3)


def unit_c(o, e, s_in, word):
    """complex developed shift of one unit from copy o (offset as t·z^{4 oo}), returns (A with sum = t·A, d, e', s')"""
    oo, ee, s, A = 0, e, s_in, 0
    for ch in word:
        s2 = int(ch)
        A += -z ** (4 * oo) * z ** -2 * (conj(TS[s2]) if ee else TS[s2])
        (_, oo), ee = across((0, oo), ee, s2)
        s = s2
    return sp.expand(A), oo, ee, s


def run_c(st, word, c):
    C, o, e, s = st
    A, d, e1, s1 = unit_c(o, e, s, word)
    if c == (0, 1):
        return (C + zpow(o) * A, (o[0], o[1] + d), e1, s1)
    assert e1 == e, 'unit flips e'
    g, dl = c
    r = z ** (4 * d)
    rc = v ** (d * g) * z ** (4 * d * dl)          # r^{gamma n + delta}
    return (C + zpow(o) * A * (rc - 1) / (r - 1), (o[0] + d * g, o[1] + d * dl), e, s1)


def states(G):
    st = [(sp.Integer(0),) + START]
    for w, c in G:
        st.append(run_c(st[-1], w, c))
    return st


def step(st, word, k):
    C, o, e, s = st
    for ch in word[:k]:
        s2 = int(ch)
        C = C + gap(o, e, TS[s2])
        o, e = across(o, e, s2)
        s = s2
    return (C, o, e, s)


def zero(expr):
    return sp.simplify(sp.together(sp.expand(expr))) == 0


def chord(b):          # i·c_b in developed units: c2 (z^{2b} − z^{−2b}); b = beta_n n + beta_0 via v = z^{4n}
    bn, b0 = b
    c2 = z ** 4 + z ** -4
    w = v ** sp.Rational(bn, 2) * z ** (2 * b0)
    return c2 * (w - 1 / w)


def symbolic():
    s0, s1 = states(G0), states(G1)
    # positions (state at the position = copy of the next step, its entry side)
    P = {'J0_4n+2': s0[2], 'J0_12n-2': s0[7], 'J0_16n-1': step(s0[8], '10', 1), 'J0_18n-1': step(s0[9], '23', 1),
         'J0_22n-2': s0[10], 'J0_end': s0[12], 'J1_18n-1': step(s1[9], '13', 1), 'J1_end': s1[10],
         'J1_6n+2': s1[4], 'J1_10n-1': step(s1[5], '20', 1)}
    Dv = lambda key, i: gap(P[key][1], P[key][2], VS[i]) + P[key][0]           # developed position of vertex V_i there
    Xedge = lambda base, key: base - Xs(P[key][0])
    Xv = lambda key, i: Xs(gap(P[key][1], P[key][2], VS[i]))
    # returns and pins (s935)
    assert P['J0_end'][1:] == START and P['J1_end'][1:] == START
    lo0 = Xv('J0_4n+2', 0) + Xs(P['J0_4n+2'][0])
    # hi0: J0 step 18n-2 exits side 2 at V2 -> position 18n-1 entry side 2 holds V2
    hi0 = Xv('J0_18n-1', 2) + Xs(P['J0_18n-1'][0])
    hi1 = Xv('J1_6n+2', 0) + Xs(P['J1_6n+2'][0])
    out = {}
    out['returns'] = [zero((lo0 - Xs(P['J0_end'][0])) - (lo0 + hi1 - hi0)), zero((hi0 - Xs(P['J1_end'][0])) - lo0)]
    # the five non-pin touches (the pins hold by definition): edge X − vertex X == 0
    out['touches'] = {'A J1.lo 18n-1 V2': zero(Xedge(hi0, 'J1_18n-1') - Xv('J1_18n-1', 2)),
                      'A J0.lo 12n-2 V0': zero(Xedge(lo0, 'J0_12n-2') - Xv('J0_12n-2', 0)),
                      'A J0.lo 16n-1 V2': zero(Xedge(lo0, 'J0_16n-1') - Xv('J0_16n-1', 2)),
                      'B J0.hi 22n-2 V0': zero(Xedge(hi0, 'J0_22n-2') - Xv('J0_22n-2', 0)),
                      'B J1.hi 10n-1 V2': zero(Xedge(hi1, 'J1_10n-1') - Xv('J1_10n-1', 2))}
    # holonomies (developed; across J the J-crossing point cancels: + C_end of the first block)
    hol = {'A1': P['J1_end'][0] - Dv('J1_18n-1', 2) + Dv('J0_4n+2', 0),
           'A2': Dv('J0_12n-2', 0) - Dv('J0_4n+2', 0),
           'A3': Dv('J0_16n-1', 2) - Dv('J0_12n-2', 0),
           'B1': Dv('J0_22n-2', 0) - Dv('J0_18n-1', 2),
           'B2': P['J0_end'][0] - Dv('J0_22n-2', 0) + Dv('J1_6n+2', 0),
           'B3': Dv('J1_10n-1', 2) - Dv('J1_6n+2', 0)}
    want = {'A1': chord((4, 3)), 'A2': chord((4, 1)) + chord((4, -1)),         # ℓ_{m-1} = c_{4n+1}... no: ℓ_{m-1} = c_{4m+1} + c_{4m+3} = c_{4n-3} + c_{4n-1}
            'A3': chord((4, 1)), 'B1': chord((4, -1)), 'B2': chord((4, 1)) + chord((4, 3)), 'B3': chord((4, -3))}
    want['A2'] = chord((4, -3)) + chord((4, -1))
    out['holonomy'] = {k: ('+' if zero(hol[k] - want[k]) else '-' if zero(hol[k] + want[k]) else 'FAIL') for k in hol}
    out['vertical'] = {k: zero(hol[k] + conj(hol[k])) for k in hol}
    ids = {'A1': ('J1_18n-1', 2), 'A2': ('J0_4n+2', 0), 'A3': ('J0_12n-2', 0), 'B1': ('J0_18n-1', 2), 'B2': ('J0_22n-2', 0),
           'B3': ('J1_6n+2', 0)}
    out['ids'] = {k: (P[key][1], P[key][2], vi) for k, (key, vi) in ids.items()}
    sh = lambda o: (o[0], o[1] - o[0])                                           # the same linear form evaluated at n − 1
    I_ = out['ids']
    out['neighbour_match'] = {'A2 = T_{m-1}.B2': (I_['A2'][0], I_['A2'][1:]) == (sh(I_['B2'][0]), I_['B2'][1:]),
                              'B1 = T_{m-1}.A1': (I_['B1'][0], I_['B1'][1:]) == (sh(I_['A1'][0]), I_['A1'][1:]),
                              'B3 = T_{m-1}.A3': (I_['B3'][0], I_['B3'][1:]) == (sh(I_['A3'][0]), I_['A3'][1:])}
    return out


def per_n(nmax):
    """exact Laurent check per n (formal in N): touch set, piece list ⊕ holonomy identification, exposure map (float
    adjacency — a cross-check of STEP 4, not used by the proof)."""
    import frontier_template as F
    rows = {}
    for n in range(1, nmax + 1):
        hc = F.holonomy_check(n)
        st = F.structure(n)
        rows[n] = dict(pieces=[(r['path'], r['steps'], r['ident'], r['vertical'], [list(r['id'][0]), r['id'][1]]) for r in hc],
                       touches={P: [list(c) for c in st[P]['cuts']] for P in 'AB'})
    return rows


def symbolic_ids():
    s0, s1 = states(G0), states(G1)
    P = {'A1': step(s1[9], '13', 1), 'A2': s0[2], 'A3': s0[7], 'B1': step(s0[9], '23', 1), 'B2': s0[10], 'B3': s1[4]}
    vi = {'A1': 2, 'A2': 0, 'A3': 0, 'B1': 2, 'B2': 0, 'B3': 0}
    return {k: (P[k][1], P[k][2], vi[k]) for k in P}


def control():
    """could it have come out otherwise: wrong chord index, wrong vertex, wrong neighbour shift must all FAIL."""
    global VS
    s0, s1 = states(G0), states(G1)
    a = step(s0[8], '10', 1); b = s0[7]
    D = lambda st, i: gap(st[1], st[2], VS[i]) + st[0]
    h = D(a, 2) - D(b, 0)                                   # A3, true value +c_{4n+1}
    r = {'A3 vs c_{4n+3}': zero(h - chord((4, 3))) or zero(h + chord((4, 3))),
         'A3 vs c_{4n-1}': zero(h - chord((4, -1))) or zero(h + chord((4, -1))),
         'A3 ending at V1 instead of V2': zero((D(a, 1) - D(b, 0)) - chord((4, 1))) or zero((D(a, 1) - D(b, 0)) + chord((4, 1)))}
    lo0 = Xs(gap(s0[2][1], s0[2][2], VS[0])) + Xs(s0[2][0])
    r['touch J0.lo 12n-2 at V3'] = zero(lo0 - Xs(b[0]) - Xs(gap(b[1], b[2], VS[3])))
    I_ = symbolic_ids()
    sh2 = lambda o: (o[0], o[1] - 2 * o[0])                                     # the linear form at n − 2 (wrong neighbour)
    r['neighbour shift n-2'] = all((I_[a][0], I_[a][1:]) == (sh2(I_[b][0]), I_[b][1:]) for a, b in (('A2', 'B2'), ('B1', 'A1'), ('B3', 'A3')))
    return {'cases': r, 'all_fail': not any(r.values())}


def n1_case():
    # (superseded by the corner argument, audit row 20; kept as a remark)
    f = lambda x: math.sin(5 * x) - math.sin(3 * x) - math.sin(x)
    # f = sin5x − sin3x − sinx = 2cos4x sin x − sin x = sin x (2cos 4x − 1) > 0 iff cos 4x > 1/2 iff x < π/12 ⊃ (0, π/15]
    return {'identity': 'sin5x − sin3x − sinx = sin x (2cos4x − 1) > 0 on (0, π/12)', 'f(pi/15)': f(math.pi / 15)}


def first_return_check():
    """no strip other than step 0 is in the start copy (2n − 1, 1) entering through side 3, for EVERY n ≥ 1 — ANALYTIC (audit
    s941 row 14; `data/scratch/s941/audit/first_return2.py`): every (run, unit position, repetition class j = 0 / j ≥ 1) whose strip
    has e = 1 and entry side 3 — the entry side of a repetition j ≥ 1 is the unit's LAST exit side — gives (A n + B) + j d = 0;
    the verdict is read off A, B, d, never by a scan. Returns the classes with a solution n ≥ 1 other than step 0."""
    from fractions import Fraction as Fr
    def acr(oo, ee, s):
        return (oo + DM[s], 1) if ee == 0 else (oo - DM[s], 0)
    bad = []
    for name, G in (('J0', G0), ('J1', G1)):
        o, e, s = START
        for idx, (w, c) in enumerate(G):
            oo, ee, pos = 0, e, []
            for ch in w:
                pos.append((oo, ee)); oo, ee = acr(oo, ee, int(ch))
            d = oo; g, dl = c
            for p, (rel, ee_p) in enumerate(pos):
                cls = [('j=0', s if p == 0 else int(w[p - 1]))]
                if c != (0, 1):
                    cls.append(('j>=1', int(w[-1]) if p == 0 else int(w[p - 1])))
                for cn, ent in cls:
                    if ee_p != 1 or ent != 3:
                        continue
                    A, B = o[0] - 2, o[1] + rel + 1
                    if cn == 'j=0':
                        if A == 0 and B == 0:
                            if not (idx == 0 and p == 0):
                                bad.append((name, idx, p, cn, 'all n'))
                        elif A and Fr(-B, A) >= 1 and Fr(-B, A).denominator == 1:
                            bad.append((name, idx, p, cn, f'n={Fr(-B, A)}'))
                    else:
                        if d == 0:
                            if A == 0 and B == 0 or (A and Fr(-B, A) >= 1 and Fr(-B, A).denominator == 1):
                                bad.append((name, idx, p, cn, 'd=0'))
                        else:
                            sl, t0 = Fr(-A, d), Fr(-B, d)
                            if sl != 0 or (t0 >= 1 and t0.denominator == 1):
                                bad.append((name, idx, p, cn, f'CHECK slope {sl} const {t0}'))
            if c == (0, 1):
                o, e = (o[0], o[1] + d), ee
            else:
                o = (o[0] + d * g, o[1] + d * dl)
            s = int(w[-1])
    return bad

if __name__ == '__main__':
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    out = {'symbolic_n_ge_2': symbolic()}
    print(json.dumps(out['symbolic_n_ge_2'], default=str, indent=1))
    out['per_n'] = per_n(nmax)
    for n, r in out['per_n'].items():
        print(n, [(p[0], p[1], p[2]) for p in r['pieces']])
    out['n1'] = n1_case()
    s = out['symbolic_n_ge_2']
    ok = all(s['returns']) and all(s['touches'].values()) and all(v != 'FAIL' for v in s['holonomy'].values()) \
        and all(s['vertical'].values()) and all(s['neighbour_match'].values())
    out['first_return_violations'] = first_return_check()
    ok = ok and not out['first_return_violations']
    out['control'] = control()
    ok = ok and out['control']['all_fail']
    out['verdict'] = ok
    print('VERDICT', ok)
    json.dump(out, open('data/s941_chain_frontier.json', 'w'), default=str, indent=1)
