"""[Deposit] Appendix B, Proposition B.1: the eight transcendental boundaries never pass through a lattice point. Paper: "Complete periodicity of perpendicular billiards in rational triangles of prime denominator". Research name: probes/s947_q3_graze.py, copied verbatim except imports (deposit names), paths (data/ and figs/ under triangles/) and this paragraph. Run through run_all.py; the usage lines below are the research repo's. Bracketed IDs and section marks (§sNNN) refer to the project's internal ledgers.

s947 `s947-a` part (3): THE q = 3 FAMILIES' TRANSCENDENTAL WORD-FAMILY BOUNDARIES NEVER PASS THROUGH A VERTEX IMAGE — each
boundary function is a single-frequency sum of 4 cosines, nonzero at EVERY integer point (m ≥ 1, t ∈ ℤ).

SETTING: `probes/s947_q3_families.py` (F1 = (3, 6m+2, 12m+4), F2 = (3, 6m+1, 12m+5), N = 18m + 9): 4 boundaries per family fail the
held-out floor fit. At such a boundary the walk switches family because the developed line passes a vertex image; the boundary
function f(m, t) = 2 Im(conj ζ^d (p − S)) is the lower family's EARLIEST condition (along the path) that fails at t = hi(m) + 1
(`engine/word_family.py` conditions; found automatically, the same at every m = 7 … 30, and equal up to sign to the upper family's
earliest failure at lo(m) − 1 — one vertex image seen from both sides), and f·(centre sign) changes sign
between t = hi and hi + 1 at every m = 7 … 30 (asserted). A lattice point with f = 0 would be an extra vertex connection.
RESULT (x = π/2N, ξ = e^{ix}): exactly, f = Σ_k A_k cos(φ_k + 6tx + δ_k x), 4 groups, A_k ∈ {2, 2√3}, φ_k ∈ (π/6)ℤ (the same formal
function in all four parity classes, asserted). As roots: P_k = A_k e^{iφ_k} ξ^{6t + δ_k} (φ_k = k·(6m + 3)x since π/6 = (N/3)x),
a 2√3 term = two roots at ratio e^{iπ/3}; f = ½ Σ (P_k + conj P_k), ≤ 10 roots of order dividing 8N.
- MANN (1965): a minimal vanishing sub-sum of ≤ 10 roots of unity has all ratios in μ_210. For m off E (`exceptional`) no ratio
  between roots of DIFFERENT terms on one side (t-free exponent αm + β) lies in μ_210, so a minimal piece holds roots of ≤ 1 term
  per side; one side alone cannot vanish (2ζ, 2ζ(1 + e^{iπ/3}) ≠ 0), and a 2√3 term cannot cancel a weight-1 root (|·| differ;
  its two roots split across pieces would need a same-side ratio e^{±iπ/3} between two other terms' roots — excluded off E). So f = 0
  needs a perfect matching P_i = −conj P_j of equal-weight terms (i = j allowed): 12t + o_i + o_j ≡ 2N (mod 4N) each; pairwise
  differences are t-free ⇒ finitely many m, each scanned over t mod 4N (`matchings`): NONE solvable.
- m ∈ E (= {2, 3, 17} on all eight — audit s948 F42 corrected the 2√3 roots, which had given ∅ on two): rigorous interval evaluation (mpmath iv) at every
  t of the period 2N/3: min |f| ≥ 7.6·10⁻⁴.
SCOPE: the 8 boundary functions of the families MEASURED on m ≤ 30 (no word transitions there; the all-m family structure is part (1)
of `s947-a`, not proved here). PROVED: those functions never vanish at an integer point. Template: `probes/s939_q5_graze.py`
([POLY-Q5-GRAZE]).
PRIOR ART: §s939 (3) ([POLY-Q5-GRAZE]: the same mechanism, q = 5, 3/4-cosine forms), [NCYL-734]-free (no holonomy); `rulings.py
--grep` *vanishing sum*, *Mann*, *graze* → none on this object. DATA: `data/s947_q3_walks_*.json`, `data/s947_q3_families_*.json`
(this session). Output: `data/s947_q3_graze.json`. Runtime ~ 3 min."""
import cmath, itertools, json, sys
import _setup  # noqa: F401  (deposit: lib/ on the path, cwd = triangles/)
from fractions import Fraction as Fr
import mpmath as mp
from mpmath import iv
import sympy as sp
import word_family as wf
import appB_1_families as QF
mp.mp.dps = 40
iv.dps = 30
FAMS = {'F1': [(0, 3), (6, 2), (12, 4)], 'F2': [(0, 3), (6, 1), (12, 5)]}
PARS = ((0, 0), (0, 1), (1, 0), (1, 1))


def words_for(famname):
    fam = wf.Fam(FAMS[famname])
    cls = json.load(open(f'data/s947_q3_families_{famname}.json'))
    out = []
    for T in cls['families']:
        v, res = T['vertex'], T['strand']
        runs = [wf.E(*[Fr(x) for x in T['runs'][str(c)]]) for c in range(T['nruns'])]
        rr = [(r, 0 if f == '-' else int(f)) for r, f in zip(runs, T['first_letters'])]
        d = wf.eadd(fam.idx[v], wf.E(0, 6, 1 + 2 * res))
        out.append((T, wf.Word(fam, v, d, rr, T['kind'])))
    return fam, cls, out


def failing(fam, w, m, t):
    conds, _ = w.conditions((m % 2, t % 2))
    out = []
    for cd in conds:
        name, kind, f, c = cd[:4]
        if len(cd) > 4 and wf.val(cd[4][1], m, t) < 2:
            continue
        if kind in ('opp', 'same'):
            fv, cv = wf.evalf(f, fam, m, t), wf.evalf(c, fam, m, t)
            if not ((fv * cv < 0) if kind == 'opp' else (fv * cv > 0)):
                out.append(name)
    return out


def reduce(f, fam):
    Nm, N0 = fam.N[0], fam.N[2]
    grp, z = {}, sp.symbols('z')
    for (a, b, g), c in f.items():
        de = Fr(g) - Fr(a) * N0 / Nm
        if b < 0 or (b == 0 and de < 0):
            a, b, g, de = -a, -b, -g, -de
        th = (Fr(a) / (2 * Nm)) % 2
        grp.setdefault((b, de), {})
        grp[(b, de)][th] = grp[(b, de)].get(th, 0) + c
    out = {}
    for k, d in grp.items():
        d = {th: c for th, c in d.items() if c}
        if not d:
            continue
        den = 1
        for th in d:
            den = sp.ilcm(den, th.denominator)
        L = 2 * den
        e = sum(c * z ** int((th * den) % L) for th, c in d.items())
        if k == (0, 0):
            e = sum(c * (z ** int((th * den) % L) + z ** int((-th * den) % L)) for th, c in d.items())
        if sp.rem(sp.expand(e), sp.cyclotomic_poly(L, z), z) == 0:
            continue
        out[k] = complex(sum(c * sp.exp(sp.I * sp.pi * th) for th, c in d.items()).evalf(30))
    return out


def terms_of(red):
    out = []
    for (b, de), amp in sorted(red.items()):
        assert b == 6
        A, ph = abs(amp), cmath.phase(amp) / cmath.pi
        k6 = round(ph * 6)
        assert abs(ph * 6 - k6) < 1e-12
        w = 'w1' if abs(A - 2) < 1e-12 else ('w3' if abs(A - 2 * 3 ** 0.5) < 1e-12 else None)
        assert w, A
        out.append((w, (Fr(6 * k6), 3 * k6 + de)))       # offset o = φ/x + δ = k6(6m + 3) + δ
    return out


def roots_of(terms):
    R = []
    for idx, (w, (am, a0)) in enumerate(terms):
        if w == 'w1':
            R.append((idx, (am, a0)))
        else:                                             # φ ∓ π/6, π/6 = (6m + 3) x  (audit s948 F42: was π/12)
            R += [(idx, (am - 6, a0 - 3)), (idx, (am + 6, a0 + 3))]
    return R


def exceptional(terms, M=210):
    E = set()
    for (i, (a1, b1)), (j, (a2, b2)) in itertools.combinations(roots_of(terms), 2):
        if i == j:
            continue
        al, be = a1 - a2, b1 - b2                         # ratio ξ^{al m + be}, order | 8N; ∈ μ_M ⟺ M(al m + be) ∈ 4N ℤ = (72m + 36)ℤ
        assert abs(al * M) < 72 * 399
        for k in range(-400, 401):
            den = al * M - 72 * k
            if den == 0:
                assert be * M - 36 * k != 0, 'ratio in μ_M for every m'
                continue
            mm = Fr(36 * k - be * M) / den
            if mm.denominator == 1 and mm >= 1:
                E.add(int(mm))
    return sorted(E)


def matchings(terms):
    n, sols = len(terms), []
    for perm in itertools.permutations(range(n)):
        if any(terms[i][0] != terms[perm[i]][0] for i in range(n)):
            continue
        eqs = [(terms[i][1][0] + terms[j][1][0], terms[i][1][1] + terms[j][1][1]) for i, j in enumerate(perm)]
        cand = None
        for (a1, b1), (a2, b2) in itertools.combinations(eqs, 2):
            al, be = a1 - a2, b1 - b2
            ms, allm = set(), False
            for k in range(-60, 61):
                den = al - 72 * k
                if den == 0:
                    allm |= (be - 36 * k == 0)
                    continue
                mm = Fr(36 * k - be) / den
                if mm.denominator == 1 and mm >= 1:
                    ms.add(int(mm))
            assert abs(al) < 72 * 59
            if not allm:
                cand = ms if cand is None else cand & ms
        if cand is None:
            sols.append((perm, 'no t-free consequence'))
            continue
        for mm in sorted(cand):
            N = 18 * mm + 9
            for tt in range(4 * N):
                if all(Fr(12 * tt + a * mm + b - 2 * N) % (4 * N) == 0 for a, b in eqs):
                    sols.append((perm, mm, tt))
                    break
    return sols


def f_iv(terms, m, t):
    N = 18 * m + 9
    x, tot = iv.pi / (2 * N), iv.mpf(0)
    for w, (am, a0) in terms:
        e = Fr(6 * t) + am * m + a0
        tot += (iv.mpf(2) if w == 'w1' else 2 * iv.sqrt(3)) * iv.cos(iv.mpf(e.numerator) / e.denominator * x)
    return tot


def main():
    out, allok = {}, True
    for famname in FAMS:
        fam, cls, W = words_for(famname)
        walks = json.load(open(f'data/s947_q3_walks_{famname}.json'))
        lows = sorted({h[0] for h in cls['holdout_fail'] if h[1] == 't_hi'})
        for i in lows:
            (T, w), (T2, w2) = W[i], W[i + 1]
            names, names2 = set(), set()
            for m in range(7, 31):
                seg = QF.segments(walks[str(m)]['rows'])
                thi = seg[i][2][-1][0]
                assert seg[i + 1][2][0][0] == thi + 1 and not failing(fam, w, m, thi) and not failing(fam, w2, m, thi + 1)
                fl, fl2 = failing(fam, w, m, thi + 1), failing(fam, w2, m, thi)
                assert fl and fl2
                names.add(fl[0]); names2.add(fl2[0])               # the EARLIEST failing condition along the path (later ones are downstream)
            assert len(names) == 1 and len(names2) == 1
            nm, nm2 = names.pop(), names2.pop()
            forms = []
            for par in PARS:
                f = [cd[2] for cd in w.conditions(par)[0] if cd[0] == nm][0]
                forms.append(terms_of(reduce(f, fam)))
            assert all(F == forms[0] for F in forms), 'parity-dependent boundary function'
            terms = forms[0]
            f2 = [cd[2] for cd in w2.conditions((0, 0))[0] if cd[0] == nm2][0]
            t2 = terms_of(reduce(f2, fam))
            neg = [(w_, (a, b + 0)) for w_, (a, b) in t2]
            same = sorted(map(str, t2)) == sorted(map(str, terms)) or all(            # same vertex image: equal up to sign
                abs(f_iv(t2, m, t).mid + s_ * f_iv(terms, m, t).mid) < 1e-20 for m, t in ((7, 3), (12, 40)) for s_ in (1,)) or all(
                abs(f_iv(t2, m, t).mid - f_iv(terms, m, t).mid) < 1e-20 for m, t in ((7, 3), (12, 40)))
            assert same, ('upper family sees a different boundary function', famname, i, nm, nm2)
            for m, t in ((7, 3), (11, 20), (30, 101), (8, 50)):          # the reduced form = the formal sum (exact grouping; sanity)
                f = [cd[2] for cd in w.conditions((m % 2, t % 2))[0] if cd[0] == nm][0]
                assert abs(wf.evalf(f, fam, m, t) - mp.mpf(f_iv(terms, m, t).mid)) < mp.mpf(10) ** -20
            # EXACT: f ≡ Σ over roots_of(terms) of 2 cos((6t + o)x)  (2√3 cos θ = 2cos(θ − π/6) + 2cos(θ + π/6)) — audit s948 G39
            for par in PARS:
                fp = [cd[2] for cd in w.conditions(par)[0] if cd[0] == nm][0]
                rec = {}
                for _, (am, a0) in roots_of(terms):
                    e = (Fr(am), Fr(6), Fr(a0))
                    rec[e] = rec.get(e, 0) + 2
                assert wf.identity_zero(wf.add(fp, rec, -1), fam)[0], ('4-term form not exact', famname, i, par)
            E = exceptional(terms)
            sol = matchings(terms)
            solv = [s for s in sol if len(s) == 3]
            nofree = [s for s in sol if len(s) == 2]
            mins = {}
            for m in E:
                N = 18 * m + 9
                vals = [f_iv(terms, m, t) for t in range(2 * N // 3)]
                mins[m] = (all(v.a > 0 or v.b < 0 for v in vals), float(min(min(abs(v.a), abs(v.b)) for v in vals)))
            ok = not solv and not nofree and all(v[0] for v in mins.values())
            allok &= ok
            key = f'{famname} {i}|{i + 1}'
            out[key] = {'condition': nm, 'terms': [(w_, f'{a}m + {b}') for w_, (a, b) in terms], 'E': E,
                        'solvable_matchings': [str(s) for s in solv], 'E_min_abs_f': {str(k): v[1] for k, v in mins.items()}, 'ok': ok}
            print(key, nm, '| f = Σ A cos((6t + o)x):', [(w_, f'{a}m+{b}') for w_, (a, b) in terms], '| E', E,
                  '| solvable matchings', len(solv), '| min|f| on E', min((v[1] for v in mins.values()), default='-'),
                  'OK' if ok else 'FAIL', flush=True)
    json.dump(out, open('data/s947_q3_graze.json', 'w'), indent=1)
    assert allok
    print('ALL OK: no lattice point on any transcendental boundary, both families, every m ≥ 1')


if __name__ == '__main__':
    main()
