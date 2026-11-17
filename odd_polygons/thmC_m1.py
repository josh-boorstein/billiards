"""[Deposit] Theorem 5.4 at m = 1 on the open interval to π/14 (Appendix B.1). Paper: "Perpendicular billiards in rational polygons of odd denominator". Research name: probes/s952_pi_family.py, copied verbatim except imports (deposit names), paths (data/ and figs/ under odd_polygons/) and this paragraph. Run through run_all.py; the usage lines below are the research repo's. Bracketed IDs and section marks (§sNNN) refer to the project's internal ledgers.

s952 keeper (queue `s951-d`): the apex-4 row (N−8, 8, N−4, 4)/N is Π_N's polygon at x → 2x — and more generally the
quadrilateral Π(x) = (π−4x, 4x, π−2x, 2x), sides 1 : 2c−1 : 2c : 2c (c = cos 2x), isosceles apex |V2V3| = |V3V0|, carries
§s935's torus for every x in s935's certified range, because s935's certificate never uses ζ^N = i except to NAME copies.

STAGE A — PRE-REGISTERED (s952, written before any row ran):
  (A1) for a ∈ {2, 3, 4, 5}, odd N > 14a with gcd(a, N) = 1, Π(aπ/N) = (N−4a, 4a, N−2a, 2a)/N (census
       quad(N, N−6a, 4a, N−2a, N−2a, 2a, N−2a)) is NON-CP in its perpendicular class, and at least one non-periodic component
       has 2-label rotation number 4 sin²(aπ/N). Rows: the first two admissible N above 14a.
  (A2) at the largest admissible N < 14a (x > π/14, s935's torus inequalities fail) the perpendicular class is CP —
       the s946 onset law at a = 2 (CP at N = 27) extended; a weaker outcome (non-CP but no component at 4 sin²(aπ/N))
       would still separate the torus from the row and is recorded as such.
  (A3) CONTROL OF THE CLASS ARGUMENT: at x = π/30 (a = 1, "N = 30" → the polygon (13, 2, 14, 1)/15) the torus direction
       i·e^{ix} = ζ_60^16 is NOT perpendicular (perp exponents 15 + 2m are odd): predict ζ_60^16 non-CP with a component
       at 4 sin²(π/30), and the perpendicular class CP (s946 census: nothing at apex 1). Could fail if the torus is not
       a local object (then the s935 transfer argument is wrong).
STAGE B — the proof (formal ⊕ rigorous): s935's identities are identities in the FREE Laurent ring ℤ[i][z, 1/z] (its class
  `L` never reduces), so they hold at z = e^{ix/2} for EVERY real x; the copy rule g_across is reflection composition, a
  free identity too, so the formal walk evaluated at x is a chain of reflected copies of Π(x) — a billiard trajectory.
  The 236 inequalities are re-certified on the whole open interval (0, π/14): s935's Taylor ⊕ interval Horner on (0, XA],
  direct interval subdivision on [XA, XB], and at the right end, for each f with f(π/14) = 0 (exact, mod Φ_56), f' < 0
  on [XB, π/14] by interval arithmetic. ⟹ for every x ∈ (0, π/14) the first return to J is the rotation by 4 sin²x.
  Rational x = aπ/N in lowest terms: N odd ⟹ i·e^{ix} = ±i·ζ_N^k (perpendicular class); 4 sin²x irrational for N ∉ {1, 2, 3, 4, 6} (here N ≥ 15).
PRIOR ART: §s935 (1) ([POLY-ASYM-ALL], its range note: inequalities fail exactly at x = π/14), §s946 (3) ([POLY-PI-SERIES]:
  apex-4 onset 29; census apexes ≤ 4 only), §s951 (2) (half-integer o₀ = reversed flow); queue s951-d; `rulings.py --grep`
  apex → nothing on Π(x); grep `x → 2x`, `N−4a`, `Π(x)` → only the s951-d row itself. DATA: `data/s946_apex_census.json`
  (apex-4 rows N ≤ 51), `data/s935_pi_n_torus.json`.
→ general_triangles.md §s952.   Usage: python probes/s952_pi_family.py [A|B|AB]   → data/s952_pi_family.json
"""
import _setup  # noqa: F401  (deposit: this directory on the path, cwd = odd_polygons/)
import json
import math
import sys
import time

import flow_poly as FP
import flow_decomp as FD
import pi_n_torus as S
from diagonal_census import quad

OUT = 'data/s952_pi_family.json'


def pi_poly(N, a2):
    """Π(x) at x = (a2/2)·π/N: a2 = 2a (a integer) gives (N−4a, 4a, N−2a, 2a)/N; a2 odd only for the (13,2,14,1)/15 control"""
    return quad(N, N - 3 * a2, 2 * a2, N - a2, N - a2, a2, N - a2)


def components(P, j, cap=50000):
    pieces, lengths, top, bot, meta = FP._build_iet_poly(P, j)
    F = FD._fast(P)
    Dc = S.PosDec(F, lengths)
    comps = Dc.run(top, bot, cap, meta['L'])
    kinds, rhos = {}, []
    for c in comps:
        kinds[c.kind] = kinds.get(c.kind, 0) + 1
        if c.kind != 'cylinder':
            s = Dc.SNAP.get(id(c))
            if s is None:
                rhos.append(None)
                continue
            lam = [F.fval(s['lam'][t]) for t in s['top']]
            rhos.append(min(lam) / sum(lam))
    return kinds, rhos


def shape_err(P, N, a2):
    """the census quad IS Π(x): side ratios 1 : 2c−1 : 2c : 2c (c = cos 2x) — angles fixed, so congruent up to scale"""
    import cmath
    L = [abs(sum(c * cmath.exp(2j * math.pi * e / P.N) for e, c in l.items())) for l in P.L]
    c = math.cos(a2 * math.pi / N)
    return max(abs(L[i] / L[0] / r - 1) for i, r in enumerate([1, 2 * c - 1, 2 * c, 2 * c]))


def row(N, a2, j=None):
    t = time.time()
    P = pi_poly(N, a2)
    assert list(P.a) == [N - 2 * a2, 2 * a2, N - a2, a2] and shape_err(P, N, a2) < 1e-12, 'not Π(x)'
    js = [j] if j is not None else FP.perp_dirs_poly(P)[:1]       # odd N: all four perpendiculars are ONE class
    target = 4 * math.sin(a2 * math.pi / (2 * N)) ** 2
    out = {'N': N, 'a2': a2, 'angles': list(P.a), 'target_rho': target, 'dirs': {}}
    for jj in js:
        kinds, rhos = components(P, jj)
        out['dirs'][str(jj)] = {'kinds': kinds, 'rhos': rhos,
                                'n_at_target': sum(1 for r in rhos if r is not None and abs(r - target) < 1e-9)}
    out['nonCP'] = any(v['kinds'].get('no_periodic', 0) or v['kinds'].get('undetermined', 0) for v in out['dirs'].values())
    out['torus_seen'] = any(v['n_at_target'] for v in out['dirs'].values())
    out['secs'] = round(time.time() - t, 1)
    print(f"N={N} a2={a2} angles={out['angles']} nonCP={out['nonCP']} torus={out['torus_seen']} "
          f"{ {k: (v['kinds'], v['n_at_target']) for k, v in out['dirs'].items()} } [{out['secs']}s]", flush=True)
    return out


def stage_a():
    res = {'A1': [], 'A2': [], 'A3': {}}
    for a in (2, 3, 4, 5):
        above = [N for N in range(14 * a + 1, 14 * a + 40, 2) if math.gcd(a, N) == 1][:2]
        below = [N for N in range(14 * a - 1, 0, -2) if math.gcd(a, N) == 1][0]
        for N in above:
            res['A1'].append(row(N, 2 * a))
        res['A2'].append(row(below, 2 * a))
    res['A3']['torus_dir'] = row(15, 1, j=16)
    res['A3']['perp'] = row(15, 1)
    res['verdict'] = {
        'A1': all(r['nonCP'] and r['torus_seen'] for r in res['A1']),
        'A2_CP': [(r['N'], r['a2'] // 2, not r['nonCP']) for r in res['A2']],
        'A2_no_torus': all(not r['torus_seen'] for r in res['A2']),
        'A3': res['A3']['torus_dir']['nonCP'] and res['A3']['torus_dir']['torus_seen'] and not res['A3']['perp']['nonCP']}
    print('VERDICT', res['verdict'], flush=True)
    return res


if __name__ == '__main__':
    which = sys.argv[1] if len(sys.argv) > 1 else 'AB'
    try:
        out = json.load(open(OUT))
    except FileNotFoundError:
        out = {}
    if 'A' in which:
        out['stage_A'] = stage_a()
        json.dump(out, open(OUT, 'w'), indent=1)


# ---------------------------------------------------------------- STAGE B: the certificate on (0, π/14)
XA = S.XMAX                                    # s935's Taylor ⊕ Horner range (0, 0.2095]


def s935_ineqs():
    """s935.main's identity asserts and inequality list, verbatim in content (the formal objects are s935's)"""
    from pi_n_torus import J0, J1, V, LEN, DIR, SIDES, g_across, g_apply, side_range, walk, pinned, sinL, c2
    assert (V[-1] + LEN[3] * DIR[3]).iszero()
    assert all((S.REF[k] - DIR[k] * DIR[k]).iszero() for k in range(4))   # reflection = DIR², no ζ^{4N} = 1 used
    for seg in (J0, J1):
        for p, q in zip(seg, seg[1:]):
            assert g_across(p[:2], p[3]) == q[:2] and p[3] == q[2]
        assert g_across(seg[-1][:2], seg[-1][3]) == (1, 1) and seg[-1][3] == 3
    for k in range(4):                         # the copy rule is reflection composition in the FREE ring (any z)
        for g in [(o, e) for o in range(-6, 7) for e in (0, 1)]:
            p = S.L({3: 1 + 2j, -5: -1})
            refl = S.REF[k] * p.conj()         # linear part of the reflection in side k
            assert (g_apply(g, refl) - g_apply(g_across(g, k), p)).iszero()
    lo0, hi0, hi1 = pinned('lo0'), pinned('hi0'), pinned('hi1')
    lam0, lam1 = hi0 - lo0, hi1 - hi0
    _, S0 = walk(J0)
    _, S1 = walk(J1)
    assert ((lo0 - S0) - (lo0 + lam1)).iszero() and ((hi0 - S1) - lo0).iszero()
    assert (lam0 - 4 * sinL(4) * sinL(1) * (c2 - 1)).iszero()
    assert (lam1 - 16 * sinL(4) * sinL(1) * sinL(1) * sinL(1)).iszero()
    ineqs = []
    for tag, seg, lo, lam in (('J0', J0, lo0, lam0), ('J1', J1, hi0, lam1)):
        Ss, _ = walk(seg)
        for k, (o, e, s, s2) in enumerate(seg):
            g = (o, e)
            a, b = lo - Ss[k], lo - Ss[k] + lam
            for ss in (s, s2):
                p, q = side_range(g, ss)
                if p.evx(math.pi / 17).real > q.evx(math.pi / 17).real:
                    p, q = q, p
                ineqs += [(f'{tag}.{k}.s{ss}.lo', a - p), (f'{tag}.{k}.s{ss}.hi', q - b)]
            for ss, want in ((s, 1), (s2, -1)):
                i, kk = SIDES[ss]
                ev = g_apply(g, V[kk]) - g_apply(g, V[i])
                ineqs.append((f'{tag}.{k}.s{ss}.{"entry" if want == 1 else "exit"}', (ev + ev.conj()) * ((1 if e == 0 else -1) * want)))
    ineqs += [('lam0', lam0), ('lam1', lam1)]
    return ineqs, lam0, lam1


def iv_eval(f, X, deriv=0):
    """rigorous enclosure of f^(deriv)(x), x ∈ X (iv), f = Σ v_b e^{i b x/2} real: Σ Re(v (ib/2)^deriv e^{ibx/2})"""
    from mpmath import iv
    acc = iv.mpf(0)
    for b, v in f.c.items():
        w = complex(v) * (1j * b / 2) ** deriv
        wr, wi = round(w.real * 2 ** deriv), round(w.imag * 2 ** deriv)      # exact: b/2 powers → scale by 2^deriv
        assert wr == w.real * 2 ** deriv and wi == w.imag * 2 ** deriv
        t = iv.mpf(b) * X / 2
        acc += (iv.mpf(wr) * iv.cos(t) - iv.mpf(wi) * iv.sin(t)) / iv.mpf(2 ** deriv)
    return acc


def pos_on(f, lo, hi, deriv=0, sign=1, depth=0):
    """sign·f^(deriv) > 0 on [lo, hi] (iv endpoints), by bisection; returns the number of leaf boxes or raises"""
    from mpmath import iv
    X = iv.mpf([lo.a, hi.b])
    v = iv_eval(f, X, deriv) * sign
    if v.a > 0:
        return 1
    assert depth < 40, 'subdivision did not certify'
    mid = (lo + hi) / 2
    return pos_on(f, lo, mid, deriv, sign, depth + 1) + pos_on(f, mid, hi, deriv, sign, depth + 1)


def zero_at_pi14(f):
    """EXACT: f(π/14) = 0 iff Φ_56 divides f(z) at z = ζ_56 (x = π/14 ⟺ z = e^{iπ/28}); i = z^14"""
    import sympy as sp
    zz = sp.Symbol('z')
    c = [0] * 56
    for b, v in f.c.items():
        c[b % 56] += int(v.real)
        c[(b + 14) % 56] += int(v.imag)
    poly = sum(ci * zz ** k for k, ci in enumerate(c) if ci)
    return sp.rem(sp.Poly(poly, zz), sp.Poly(sp.cyclotomic_poly(56, zz), zz)).is_zero


def stage_b():
    from mpmath import iv
    iv.dps = 40
    t = time.time()
    ineqs, lam0, lam1 = s935_ineqs()
    tight = [n for n, f in ineqs if f.iszero()]
    strict = [(n, f) for n, f in ineqs if not f.iszero()]
    # (i) (0, XA]: s935's Taylor ⊕ interval Horner, unchanged
    cert0 = {n: S.certify_nonneg(f) for n, f in strict}
    # (ii) [XA, XB] direct interval subdivision; (iii) [XB, π/14]: f > 0, or f(π/14) = 0 exactly and f' < 0
    PI14 = iv.pi / 14
    XAiv = iv.mpf(XA.numerator) / iv.mpf(XA.denominator)
    XB = PI14 - iv.mpf('0.004')
    zeros = [n for n, f in strict if zero_at_pi14(f)]
    boxes, endtest = 0, {}
    for n, f in strict:
        boxes += pos_on(f, XAiv, XB)
        if n in zeros:
            endtest[n] = ('f(π/14)=0 exact, f\'<0', pos_on(f, XB, PI14, deriv=1, sign=-1))
        else:
            endtest[n] = ('f>0', pos_on(f, XB, PI14))
    # control: the certificate must FAIL past π/14 (Π_13's x, the s935 control) — the open end is sharp
    ctl = {N: sum(1 for _, f in strict if f.evx(math.pi / N).real <= 0) for N in (14.5, 14.01, 13.99, 13.5, 13)}
    assert ctl[14.5] == 0 and ctl[14.01] == 0 and ctl[13.99] > 0, 'control did not fire at π/14'
    res = {'identities': 'closure, itinerary, copy rule = reflection composition (free ring), both returns, λ0, λ1: FORMAL (free Laurent)',
           'inequalities': len(ineqs), 'formally_tight': len(tight), 'strict': len(strict),
           'taylor_(0,XA]': {'XA': float(XA), 'min_lower_bound': min(v[1] for v in cert0.values()),
                             'orders': sorted({v[0] for v in cert0.values()})},
           'interval_[XA,XB]_boxes': boxes, 'XB': float(XB.a),
           'vanish_at_pi14': zeros, 'end_tests': {k: v[0] for k, v in endtest.items() if k in zeros},
           'control_violations_at_x=pi/N': {str(k): v for k, v in ctl.items()},
           'rho_check': {f'{a}/{N}': lam1.evx(a * math.pi / N).real / (lam0 + lam1).evx(a * math.pi / N).real
                         - 4 * math.sin(a * math.pi / N) ** 2 for a, N in ((2, 29), (3, 43), (5, 71), (7, 99))},
           'secs': round(time.time() - t)}
    print(json.dumps(res, indent=1), flush=True)
    return res


if __name__ == '__main__' and 'B' in (sys.argv[1] if len(sys.argv) > 1 else 'AB'):
    out = json.load(open(OUT)) if __import__('os').path.exists(OUT) else {}
    out['stage_B'] = stage_b()
    json.dump(out, open(OUT, 'w'), indent=1, default=str)
