"""[Deposit] Appendix C: the Theorem D certificate (ball arithmetic). Paper: "Complete periodicity of perpendicular billiards in rational triangles of prime denominator". Research name: probes/s946_family_cert.py, copied verbatim except imports (deposit names), paths (data/ and figs/ under triangles/) and this paragraph. Run through run_all.py; the usage lines below are the research repo's. Bracketed IDs and section marks (§sNNN) refer to the project's internal ledgers.

s946 keeper (queue s943-b): CERTIFY the four free cylinder families of [TRI-MPLUS-FREE] -- on an explicit closed interval
J of B/A on the line A = pi/n, EVERY rational triangle (A, B) (A = pi/n, odd N) carries an iota-swapped pair of FREE
perpendicular-class cylinders (core c and its cyclic reverse).

Certificate per family (core c from data/s944_window_search.jsonl; frame/letters = s922_band_shooter.develop):
 (i)  EXACT: c has even length and its linear part is the identity on the line, identically in y = e^{iB}
      (direction monomials of all three sides after c equal those before, x = zeta_2n) -- so c closes by a translation T(B).
 (ii) EXACT: T = hol(c)/(2i) has direction psi = pi/2 + s B + ct A with INTEGERS (s, ct): U := hol(c) x^-ct y^-s is real on the
      line (U = U~ exactly in Q(zeta_2n)[y^+-1]) -- so psi is a perpendicular-class direction at every rational point of the line.
 (iii) RIGOROUS (Lipschitz): the transverse strip width W(B) = min_k max(t_ik, t_jk) - max_k min(t_ik, t_jk) over the crossed side
      copies (t_v = transverse coordinate, a trigonometric polynomial in B with coefficients computed exactly in Q(zeta_2n), evaluated
      in float64 with a generous rounding bound err) is > 0 on J: Lip(W) <= 2 max_v Lip(t_v), and
      min(W(B_j), W(B_j+1)) - Lip h / 2 - 2 err > 0 on every cell; this float pass FINDS J (rounded inward to 8 decimals).
 (vi) RIGOROUS: the cell certificate re-run on the stated J in arb ball arithmetic (200 bits, exact integer Laurent
      coefficients, rigorous Lipschitz constant) -- this is the proof; (iii) and (v) are guards.
 (o)  c is a cyclic billiard itinerary (no two cyclically adjacent letters equal).
 (iv) c is NOT a cyclic rotation of its reverse.
 (v)  independent instrument: at every certified cell endpoint the float width agrees with s922_band_shooter.develop at 30 digits
      to within the rounding bound err; |T| > 0 on J (U bounded away from 0, same Lipschitz scheme).
Then at every point of J the lines with itinerary c in direction psi fill an open strip -> a cylinder C (closed leaves avoiding
every vertex point); iota reverses the flow and keeps side labels, so iota C has itinerary reverse(c); a maximal cylinder has one
cyclic itinerary (no vertex point inside it -- every vertex point is marked, incl. the regular V0 points of the line),
so (iv) => C != iota C (C = iota C would give its leaves both itineraries): two free cylinders.
CONTROLS: W < 0 just outside J (the certificate can fail); the 13 exact flow_decomp free points of [TRI-MPLUS-FREE] (§s943 (5),
§s944 (4)) must lie in the certified J of their family, or be listed as outside.
PRE-REGISTERED: (i), (ii), (iv) hold exactly for all four cores; (iii) certifies a J containing each family's exact points.
PRIOR ART: `rulings.py --grep` 'free cylinder' -> [NCYL-719] (amended s943: the free pairs exist), [NCYL-729], [NCYL-730];
'certify' family -> none; §s943 (5) gives the statement to prove and float/mp width scans only (no certificate); §s944 (4) the n=27 family.
DATA: data/s944_window_search.jsonl (cores = pre + w of the cylinder bands); no stored certificate (data_index --grep family -> none).
Output: data/s946_family_cert.json.
Run: PYTHONPATH=.:engine:archive/scripts_2026-07:archive/scripts_2026-08:probes .venv/bin/python3.13 probes/s946_family_cert.py
"""
import _setup  # noqa: F401  (deposit: lib/ on the path, cwd = triangles/)
import json
import cmath
import math
from collections import defaultdict
from fractions import Fraction

import numpy as np
import mpmath as mp

import band_shooter as SH

import laurent as LA

SIDES = [(0, 1), (1, 2), (2, 0)]

EXACT_POINTS = [  # [TRI-MPLUS-FREE] exact free points (flow_decomp), §s943 (5), §s944 (4)
    (137, 277, 819), (183, 370, 1094), (103, 414, 1028), (185, 374, 1106), (259, 1041, 2585), (365, 1467, 3643),
    (167, 673, 1665), (115, 923, 2067),
    (229, 463, 1369), (275, 556, 1644), (321, 649, 1919), (323, 653, 1931), (155, 623, 1547)]  # 13 (audit row 36)


def develop_states(word):
    """Laurent development (s944_laurent conventions, everything x 2i): list of vertex-position states and final directions."""
    D = {(0, 1): (1, 0, 0), (1, 2): (-1, 0, -1), (0, 2): (1, 1, 0)}
    for (i, j) in list(D):
        D[(j, i)] = (-D[(i, j)][0], D[(i, j)][1], D[(i, j)][2])
    D0 = dict(D)
    V = [{}, dict(LA.LEN[frozenset((0, 1))]), LA.pmul_mono(LA.LEN[frozenset((0, 2))], (1, 1, 0))]
    states = [V]
    for ch in word:
        i, j = SIDES[int(ch)]
        k = 3 - i - j
        u2 = LA._mono_mul(D[(i, j)], D[(i, j)])
        nik = LA._mono_mul(u2, LA._mono_conj(D[(i, k)]))
        njk = LA._mono_mul(u2, LA._mono_conj(D[(j, k)]))
        V = list(V)
        V[k] = LA.padd(V[i], LA.pmul_mono(LA.LEN[frozenset((i, k))], nik))
        D[(i, k)] = nik; D[(k, i)] = (-nik[0], nik[1], nik[2])
        D[(j, k)] = njk; D[(k, j)] = (-njk[0], njk[1], njk[2])
        states.append(V)
    return states, D0, D


def norm_mono(m, n):
    s, p, q = m
    f, r = divmod(p, n)
    return (s * (-1) ** f, r, q)


def linear_identity(D0, D, n):
    return all(norm_mono(D0[key], n) == norm_mono(D[key], n) for key in D0)


def on_line_exact(P, n):
    """{q: coeff tuple in Q(zeta_2n) basis} (exact)."""
    return LA.on_line(P, n)


def direction_class(hol, n):
    """Find integers (s, ct) with hol * x^-ct * y^-s real on the line; exact confirmation."""
    x = cmath.exp(1j * math.pi / n)
    # numeric: psi(B) = arg(hol/(2i)) ; s = d psi/dB (integer), ct from psi - pi/2 - sB
    def H(B):
        y = cmath.exp(1j * B)
        return sum(v * x ** p * y ** q for (p, q), v in hol.items()) / 2j
    B1, B2 = 0.3, 0.3 + 1e-6
    d = cmath.phase(H(B2) / H(B1)) / 1e-6
    s = round(d)
    rest = (cmath.phase(H(B1)) - math.pi / 2 - s * B1) / (math.pi / n)
    ct = round(rest)
    for ct_try in (ct, ct - 1, ct + 1):
        for s_try in (s,):
            U = LA.pmul_mono(hol, (1, -ct_try, -s_try))
            diff = LA.padd(U, LA.ptilde(U), -1)
            if not LA.on_line(diff, n):
                return s_try, ct_try, U
    return None


def coeffs_complex(P, n, shift_p, shift_q):
    """t-coordinate trig poly: alpha_q (complex) with value = Im(sum alpha_q e^{iqB})."""
    x = np.exp(1j * np.pi / n)
    out = defaultdict(complex)
    for (p, q), v in P.items():
        out[q + shift_q] += v * x ** (p + shift_p)
    return {q: a for q, a in out.items() if abs(a) > 0}


def ball_certify(n, c, s, ct, J0, J1, states):
    """(vi) RIGOROUS: the same cell certificate in arb ball arithmetic (python-flint, 200 bits), exact integer Laurent
    coefficients, rigorous Lipschitz constant (audit s946 row 27, `data/scratch/s946/audit/a3_ball.py`)."""
    import flint
    from flint import acb, arb, acb_mat
    flint.ctx.prec = 200
    A = arb.pi() / n
    xp = [acb(0, A * e).exp() for e in range(2 * n)]
    seen, polys, edges = {}, [], []
    for kk, ch in enumerate(c):
        i, j = SIDES[int(ch)]
        pr = []
        for v in (i, j):
            key = tuple(sorted(states[kk][v].items()))
            if key not in seen:
                seen[key] = len(polys); polys.append(states[kk][v])
            pr.append(seen[key])
        edges.append(pr)
    qs = sorted({q - s for P in polys for (p, q) in P})
    qi = {q: k for k, q in enumerate(qs)}
    rows = []
    for P in polys:
        row = [acb(0)] * len(qs)
        for (p, q), v in P.items():
            row[qi[q - s]] += v * xp[(p - ct) % (2 * n)]
        rows.append(row)
    Mb = acb_mat(rows)
    lipv = max(sum(abs(q) * rows[k][qi[q]].abs_upper() for q in qs) for k in range(len(rows)))
    LipW = lipv * A                     # = 2 * (sum |q||alpha_q| / 2) * A, per d(B/A)

    def Wlo(tB):
        B = arb(tB) * A
        T = Mb * acb_mat([[acb(0, B * q).exp()] for q in qs])
        t = [T[k, 0].imag / 2 for k in range(len(rows))]
        lo = [x.lower() for x in t]; hi = [x.upper() for x in t]
        return min(max(lo[a], lo[b]) for a, b in edges) - max(min(hi[a], hi[b]) for a, b in edges)

    stack = [(J0, J1, Wlo(J0), Wlo(J1))]
    cells, worst = 0, None
    while stack:
        u, v, wu, wv = stack.pop()
        mg = (min(wu, wv) - LipW * (arb(v) - arb(u)) / 2).lower()
        if mg > 0:
            cells += 1; worst = mg if worst is None or mg < worst else worst; continue
        assert v - u > 1e-12, ('ball certificate failed', n, u, v)
        m = (u + v) / 2; wm = Wlo(m)
        stack += [(u, m, wu, wm), (m, v, wm, wv)]
    return cells, float(worst)


def cert_family(n, c, lo, hi):
    states, D0, D = develop_states(c)
    assert len(c) % 2 == 0
    assert all(c[k] != c[(k + 1) % len(c)] for k in range(len(c))), 'not a cyclic billiard itinerary'  # audit row 19
    lin = linear_identity(D0, D, n)
    hol = LA.padd(states[-1][2], states[0][2], -1)
    trans_all = all(not LA.on_line(LA.padd(LA.padd(states[-1][v], states[0][v], -1), hol, -1), n) for v in range(3))
    dc = direction_class(hol, n)
    assert dc is not None, 'direction class not found'
    s, ct, U = dc
    rev = c[::-1]
    reversible = rev in c + c
    # transverse coordinate of a position V (x 2i): hol = U x^ct y^s, U real, so the true direction is -i x^ct y^s and
    # t = Im(x^-ct y^-s V) / 2 (up to a global sign, width-invariant)
    vid = {}
    polys = []
    edges = []
    for kk, ch in enumerate(c):
        i, j = SIDES[int(ch)]
        pair = []
        for v in (i, j):
            key = tuple(sorted(states[kk][v].items()))
            if key not in vid:
                vid[key] = len(polys)
                polys.append(coeffs_complex(states[kk][v], n, -ct, -s))
            pair.append(vid[key])
        edges.append(pair)
    qs = sorted({q for P in polys for q in P})
    qi = {q: k for k, q in enumerate(qs)}
    M = np.zeros((len(polys), len(qs)), dtype=complex)
    for r, P in enumerate(polys):
        for q, a in P.items():
            M[r, qi[q]] = a
    qv = np.array(qs, dtype=float)
    lip = float(np.max(np.abs(M) @ np.abs(qv))) / 2          # bound on |dt/dB| (t = Im(...)/2)
    # float64 rounding bound per t evaluation (generous): alpha_q built from integer Laurent coefficients (sum |v| = absV),
    # phases q*B with |q| <= qmax, B <= pi: error <= 1e-13 * absV * (1 + qmax)
    absV = max(sum(abs(v) for v in states[kk][v_].values()) for kk in range(len(states)) for v_ in range(3)) / 2
    qmax = max(abs(q) for q in qs)
    err = 1e-13 * absV * (1 + qmax)
    E = np.array(edges)
    A = math.pi / n

    def W(tB):
        B = np.asarray(tB) * A
        T = (np.exp(1j * np.outer(B, qv)) @ M.T).imag / 2       # shape (nB, nverts)
        ti, tj = T[:, E[:, 0]], T[:, E[:, 1]]
        return np.min(np.maximum(ti, tj), axis=1) - np.max(np.minimum(ti, tj), axis=1)

    # scan, then certify the maximal positive run
    grid = np.linspace(lo, hi, 4001)
    w = W(grid)
    pos = np.where(w > 0)[0]
    runs = np.split(pos, np.where(np.diff(pos) != 1)[0] + 1)
    run = max(runs, key=len)
    a_, b_ = math.ceil(grid[run[0]] * 1e8) / 1e8, math.floor(grid[run[-1]] * 1e8) / 1e8   # stated J, rounded INWARD
    LipW = 2 * lip * A                                          # d/d(B/A)
    # certify [a_, b_] adaptively: cell ok if min endpoint W - LipW*h/2 - 2 err > 0
    stack = [(a_, b_, float(W([a_])[0]), float(W([b_])[0]))]
    cells, worst, ends = 0, math.inf, set()
    while stack:
        u, v, wu, wv = stack.pop()
        h = v - u
        margin = min(wu, wv) - LipW * h / 2 - 2 * err
        if margin > 0:
            cells += 1; worst = min(worst, margin); ends.update([(u, wu), (v, wv)]); continue
        assert h > 1e-13, ('certificate failed', n, u, v, wu, wv)
        m = (u + v) / 2; wm = float(W([m])[0])
        stack += [(u, m, wu, wm), (m, v, wm, wv)]
    # (v) independent instrument: 30-digit s922_band_shooter.develop at every certified cell endpoint
    mp.mp.dps = 30
    dmax = 0.0
    for (u, wu) in sorted(ends):
        S_, cps = SH.develop(mp.pi / n, mp.mpf(u) * mp.pi / n, c)
        Tt = cps[-1][2] - cps[0][2]; d = Tt / abs(Tt)
        lo_, hi_ = -mp.inf, mp.inf
        for kk, ch in enumerate(c):
            i, j = SIDES[int(ch)]
            xi, xj = mp.im(mp.conj(d) * cps[kk][i]), mp.im(mp.conj(d) * cps[kk][j])
            lo_ = max(lo_, min(xi, xj)); hi_ = min(hi_, max(xi, xj))
        dmax = max(dmax, abs(float(hi_ - lo_) - wu))
    assert dmax < err, ('float width disagrees with the 30-digit developer', dmax, err)
    outside = [float(W([a_ - (b_ - a_) * 0.01])[0]), float(W([b_ + (b_ - a_) * 0.01])[0])]
    # |T| > 0 on J: U (real on the line) bounded away from 0, same Lipschitz scheme
    Uc = coeffs_complex(U, n, 0, 0)
    uq = np.array(sorted(Uc), dtype=float); ua = np.array([Uc[q] for q in sorted(Uc)])
    Uf = lambda tB: np.abs((np.exp(1j * np.outer(np.asarray(tB) * A, uq)) @ ua).real)
    LipU = float(np.abs(ua) @ np.abs(uq)) * A
    errU = 1e-13 * sum(abs(v) for v in U.values()) * (1 + float(np.max(np.abs(uq))))
    gU = np.linspace(a_, b_, 2001); vU = Uf(gU)
    U_ok = bool(np.min(vU) - LipU * (gU[1] - gU[0]) / 2 - errU > 0)
    ball_cells, ball_worst = ball_certify(n, c, s, ct, a_, b_, states)
    return dict(ball_cells=ball_cells, ball_worst_margin=ball_worst, n=n, L=len(c), linear_identity=lin, translation_all_vertices=trans_all, s=s, ct=ct, reversible=reversible,
                n_vertices=len(polys), q_span=[min(qs), max(qs)], lip_t=lip, err=err, J=[float(a_), float(b_)],
                width_max=float(w.max()), cells=cells, worst_margin=worst, width_just_outside=outside,
                T_nonzero=U_ok, min_abs_U=float(np.min(vU)), mp_endpoints=len(ends), mp_max_dev=dmax)


def main():
    cores = {}
    for l in open('data/s944_window_search.jsonl'):
        r = json.loads(l)
        for b in r['bands']:
            if b['cyl']:
                c = b['hits'][0]['pre'] + b['w']
                key = (r['n'], len(c))
                cores.setdefault(key, []).append((c, b['lo'], b['hi']))
    out = []
    for (n, L), lst in sorted(cores.items()):
        c = lst[0][0]
        canon = lambda z: min(min(z[i:] + z[:i] for i in range(len(z))), min(z[::-1][i:] + z[::-1][:i] for i in range(len(z))))
        assert len({canon(x[0]) for x in lst}) == 1, 'bands of one family give different cores'
        lo = min(x[1] for x in lst); hi = max(x[2] for x in lst)
        pad = 0.0008 * n / 9
        res = cert_family(n, c, lo - pad, hi + pad)
        pts = []
        for abc in EXACT_POINTS:
            N = sum(abc)
            if N % abc[0] == 0 and N // abc[0] == n:
                t = abc[1] / abc[0]
                if lo - 0.02 < t < hi + 0.02:
                    pts.append((list(abc), t, res['J'][0] <= t <= res['J'][1]))
        res['exact_points'] = pts
        # primitive odd-N points of the line in J: N = n a0 (a0 odd), gcd(a0, a1) = 1, a1 / a0 in J
        cnt = {}
        for Nmax in (10 ** 4, 10 ** 5, 10 ** 6):
            k = 0
            for a0 in range(1, Nmax // n + 1, 2):
                for a1 in range(math.ceil(res['J'][0] * a0), math.floor(res['J'][1] * a0) + 1):
                    if math.gcd(a0, a1) == 1:
                        k += 1
            cnt[Nmax] = k
        res['points_in_J'] = cnt
        res['core'] = c
        out.append(res)
        print(f"n={n} |c|={L}: linear id {res['linear_identity']}, translation {res['translation_all_vertices']}, psi = pi/2 + {res['s']}B + {res['ct']}A,"
              f" reversible {res['reversible']}; J = [{res['J'][0]:.7f}, {res['J'][1]:.7f}] certified ({res['cells']} cells,"
              f" worst margin {res['worst_margin']:.2e}; BALL {res['ball_cells']} cells, worst {res['ball_worst_margin']:.2e}, max width {res['width_max']:.2e}); just outside W = "
              f"{res['width_just_outside'][0]:+.1e}, {res['width_just_outside'][1]:+.1e}; exact pts in J: "
              f"{sum(p[2] for p in pts)}/{len(pts)} {[(p[0], round(p[1], 6), p[2]) for p in pts]}; |T|>0 {res['T_nonzero']}; err {res['err']:.1e}; "
              f"30-digit developer at {res['mp_endpoints']} endpoints: max dev {res['mp_max_dev']:.1e}; "
              f"odd-N primitive points in J (N <= 1e4/1e5/1e6): {list(res['points_in_J'].values())}")
    json.dump(out, open('data/s946_family_cert.json', 'w'), indent=1)


if __name__ == '__main__':
    main()
