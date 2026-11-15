"""[Deposit] Proposition 6.1: the exact symbolic identities (and numeric controls). Paper: "Complete periodicity of perpendicular billiards in rational triangles of prime denominator". Research name: probes/s921_band_qm.py, copied verbatim except imports (deposit names), paths (data/ and figs/ under triangles/) and this paragraph. Run through run_all.py; the usage lines below are the research repo's. Bracketed IDs and section marks (§sNNN) refer to the project's internal ledgers.

s921: the general (q, m) band -- Proposition F' (q = 3) and Prop P (q = 5) are the q = 3, 5 members of ONE word family.

Keeper for `general_triangles.md` §s921 ([TRI-BAND-QM]). For odd q >= 3, m >= 1 put n = q(2m+1), A = pi/n,
beta = pi/q, d = 4m+2, h = (q-1)/2, K = (n-2m-3)/2, gamma0 = (3n-2m-3)/2, and the word (sides 0 = V0V1, 1 = V1V2, 2 = V2V0)
    W_{q,m} = 0 | (20)^{r_1/2} | ... | (20)^{r_h/2} | (02)^{f_{h+1}/2} | ... | (02)^{f_{q-2}/2} | (02)^{2m-1}0,   '|' = the letter 1,
    r_p = p d + 2 (rising, p = 1..h),  f_p = (q-p) d - 2 (falling, p = h+1..q-2).
Claim: in the triangle with angle A at V0 and B = (K + t)A at V1, W_{q,m} is a V2 -> V2 prong in direction
psi = pi/2 + B - gamma0 A (a perpendicular-class direction at odd N) for EVERY t in (0, t_up), t_up in (2/3, 1) the root of
    (*)  sin(3beta/2 - A) sin((1-t)A) cos(beta/2 - tA) = sin(tA) cos(beta/2 + A - tA) sin(beta/2 - A),
and for no t in (-1, 0] or [t_up, 1). W is not a palindrome, so the prong is swapped (Prop S (iii)).
Part (b): the m words W_{q,m,j} (runs 2j+1 | r_p + 2j | f_p - 2j | 4m-1-2j, j = 0..m-1; j = 0 is W_{q,m}) are swapped prongs
in direction psi_j = psi - jA on the same band, so p >= m there (q = 3: Prop F'').

Checks (all must pass; prints ALL OK):
 (S) SYMBOLIC, exact (sympy, Laurent polynomials in x = e^{i beta/2}, y = e^{iA/2}, s = e^{itA/2}, z = e^{ik beta/2}):
     the closed forms of the signed distances of every load-bearing vertex, for p = 2k and p = 2k+1 with k FREE, for the first
     V0-star, and for the last V0-star (z^4 = -x^-2, i.e. h beta = pi/2 - beta/2) -- including the closure F(E) = 0;
     the t = 2/3 identity for (*); and the same with the prong index j free (v = e^{ijA/2}).
 (N) NUMERIC (mpmath, 40 digits): the closed-form V0-star data against a direct unfolding of the word; the unfolding is proper at
     t = 0.01, 1/3, 2/3, t_up - 1e-6 and NOT proper at t = -1e-6, t_up + 1e-6 (controls fire) on a (q, m) grid incl. composite q;
     closing direction psi = pi/2 + B - gamma0 A; all m prongs proper on the band, 2m distinct words, controls fire.
Run from the repo root (~15 s):
PYTHONPATH=.:engine:archive/scripts_2026-07:archive/scripts_2026-08:probes .venv/bin/python3.13 probes/s921_band_qm.py
"""
import _setup  # noqa: F401  (deposit: lib/ on the path, cwd = triangles/)
import mpmath as mp
import sympy as sp

SIDES = [(0, 1), (1, 2), (2, 0)]


def alt(start, L):
    other = '2' if start == '0' else '0'
    return ''.join(start if i % 2 == 0 else other for i in range(L))


def word(q, m):
    h = (q - 1)//2
    parts = [alt('0', 1)]
    parts += [alt('2', 4*m + 4 + i*(4*m + 2)) for i in range(h)]
    parts += [alt('0', 8*m + 2 + j*(4*m + 2)) for j in range(h - 2, -1, -1)]
    parts += [alt('0', 4*m - 1)]
    return '1'.join(parts)


def consts(q, m):
    n = q*(2*m + 1)
    return n, (n - 2*m - 3)//2, (3*n - 2*m - 3)//2


# ---------------------------------------------------------------- numeric: unfolding
def unfold(A, B, w):
    C = mp.pi - A - B
    V = [mp.mpc(0), mp.mpc(mp.sin(C)), mp.sin(B)*mp.expj(A)]
    cps = [V]
    for ch in w:
        i, j = SIDES[int(ch)]
        k = 3 - i - j
        u = (V[j] - V[i]) / abs(V[j] - V[i])
        V = list(V)
        V[k] = V[i] + u*u*mp.conj(V[k] - V[i])
        cps.append(V)
    return cps


def proper(q, m, t):
    """True iff the straight segment S -> E crosses every listed side copy strictly inside; also returns psi offset."""
    n, K, g0 = consts(q, m)
    A = mp.pi/n
    B = (K + t)*A
    w = word(q, m)
    cps = unfold(A, B, w)
    S, E = cps[0][2], cps[-1][2]
    dvec = E - S
    f = lambda z: mp.im((z - S)*mp.conj(dvec))
    ok = all(f(cps[k][SIDES[int(ch)][0]])*f(cps[k][SIDES[int(ch)][1]]) < 0 for k, ch in enumerate(w))
    c = (mp.pi/2 + B - mp.arg(dvec))/A
    return ok, c


def t_up(q, m):
    n, K, g0 = consts(q, m)
    A = mp.pi/n
    u = mp.pi/(2*q)
    g = lambda t: mp.sin(3*u - A)*mp.sin((1 - t)*A)*mp.cos(u - t*A) - mp.sin(t*A)*mp.cos(u + A - t*A)*mp.sin(u - A)
    return mp.findroot(g, (mp.mpf('0.5'), mp.mpf('0.9999')), solver='bisect')


def closed_F(q, m, t, p, ray, pos=False):
    """signed distance (left of travel +) of V0-star p's vertex on base ray `ray` (None = the centre), closed form.
    pos=True: return the position in the unfolding instead."""
    n, K, g0 = consts(q, m)
    A = mp.pi/n
    b = mp.pi/q
    B = (K + t)*A
    sB, sC = mp.sin(B), mp.sin(B + A)
    if p == 0:
        G = lambda zz: zz
    else:
        tau = lambda j: (-j*b) if j % 2 == 0 else ((j + 1)*b + 2*A)     # tau_j * A, mod 2 pi
        Op = -2j*sB*sC*mp.expj(B)*sum((-1)**j*mp.expj(tau(j)) for j in range(p))
        if p % 2 == 1:
            G = lambda zz: Op + mp.expj(2*B + tau(p - 1))*zz
        else:
            G = lambda zz: Op + mp.expj(tau(p - 1))*mp.conj(zz)
    zz = mp.mpc(0) if ray is None else (sC if ray % 2 == 0 else sB)*mp.expj(ray*A)
    S = sB*mp.expj(A)
    if pos:
        return G(zz)
    return -mp.re((G(zz) - S)*mp.expj(-B + g0*A))


def load_bearing(q, m):
    """(p, ray) of every vertex the properness argument uses (base-frame ray index, mod 2n)."""
    h = (q - 1)//2
    d = 4*m + 2
    out = [(0, 0), (0, -1)]
    for p in range(1, q - 1):
        if p <= h:
            r = p*d + 2
            out += [(p, 1), (p, r - 1), (p, r), (p, r + 1)]
        else:
            f = (q - p)*d - 2
            out += [(p, 0), (p, -1), (p, -(f - 1)), (p, -f)]
    e = 4*m - 1
    out += [(q - 1, 0), (q - 1, -1), (q - 1, -(e - 1)), (q - 1, -(e - 2)), (q - 1, -e)]
    return out


def numeric_checks():
    mp.mp.dps = 40
    bad = []
    for q, m in [(3, 1), (3, 4), (5, 1), (5, 3), (7, 1), (7, 2), (9, 1), (9, 3), (11, 2), (15, 1), (21, 1), (25, 2)]:
        n, K, g0 = consts(q, m)
        tu = t_up(q, m)
        if not (mp.mpf(2)/3 < tu < 1):
            bad.append(('t_up range', q, m, tu))
        for t in [mp.mpf('0.01'), mp.mpf(1)/3, mp.mpf(2)/3, tu - mp.mpf('1e-6')]:
            ok, c = proper(q, m, t)
            if not ok or abs(c - g0) > mp.mpf('1e-25'):
                bad.append(('inside', q, m, t, ok, c))
        for t in [mp.mpf('-1e-6'), tu + mp.mpf('1e-6')]:
            ok, c = proper(q, m, t)
            if ok:
                bad.append(('control did not fire', q, m, t))
        # closed-form V0-star centres vs the unfolding
        t = mp.mpf('0.41')
        A = mp.pi/n
        w = word(q, m)
        cps = unfold(A, (K + t)*A, w)
        S, E = cps[0][2], cps[-1][2]
        dd = (E - S)/abs(E - S)
        starts = [0] + [k + 1 for k, ch in enumerate(w) if ch == '1']
        for p, k0 in enumerate(starts):
            direct = mp.im((cps[k0][0] - S)*mp.conj(dd))
            if abs(direct - closed_F(q, m, t, p, None)) > mp.mpf('1e-30'):
                bad.append(('centre', q, m, p))
        if abs(closed_F(q, m, t, q - 1, -(4*m - 1))) > mp.mpf('1e-30'):
            bad.append(('closure', q, m))
        verts = [v for cp in cps for v in cp]
        for (p, ray) in load_bearing(q, m):
            P = closed_F(q, m, t, p, ray, pos=True)
            if min(abs(P - v) for v in verts) > mp.mpf('1e-30'):
                bad.append(('vertex not in unfolding', q, m, p, ray))
    return bad


# ---------------------------------------------------------------- symbolic
x, y, s, z = sp.symbols('x y s z', nonzero=True)
I = sp.I


def ee(v):
    return x**v[0]*y**v[1]*s**v[2]*z**v[3]


def cos_(v):
    return (ee(v) + ee(tuple(-a for a in v)))/2


def sin_(v):
    return (ee(v) - ee(tuple(-a for a in v)))/(2*I)


def add(*vs):
    return tuple(sum(c) for c in zip(*vs))


# angle vectors in units (beta/2, A/2, tA/2, k beta/2)
sB, sC, cb2, cU = cos_((1, 2, -2, 0)), cos_((1, 0, -2, 0)), cos_((1, 0, 0, 0)), cos_((1, 2, -2, 0))
F_even0 = 2*sB*sC*sin_((0, 0, 0, 2))*sin_((0, 2, 0, 2))/cb2 - sB*cos_((0, 2, -2, 0))
F_odd0 = F_even0 + 2*sB*sC*cos_((1, 2, 0, 4))


def tipF(parity, R, phi):
    if parity == 0:     # F_{2k}(R e^{i phi}) = F_{2k}(0) + R cos(2k beta + (2-t)A - phi)
        return F_even0 + R*cos_(add((0, 4, -2, 4), tuple(-a for a in phi)))
    return F_odd0 - R*cos_(add(phi, (-2, -4, 2, -4)))     # F_{2k+1}(0) - R cos(phi - (2k+1)beta - (2-t)A)


def F_first(R, phi):    # V0-star 0, frame = identity: R cos(phi - tA) - sB cos((1-t)A)
    return R*cos_(add(phi, (0, 0, -2, 0))) - sB*cos_((0, 2, -2, 0))


def is_zero(expr, last=False):
    num, den = sp.fraction(sp.together(sp.expand(expr)))
    num = sp.expand(num)
    if last:            # z^4 -> -x^-2  (k = h)
        P = sp.Poly(num, z)
        red = 0
        for (deg,), co in P.terms():
            assert deg % 4 == 0
            red += co*(-x**-2)**(deg//4)
        num = sp.expand(sp.together(red))
        num = sp.fraction(num)[0]
        num = sp.expand(num)
    return num == 0


def symbolic_checks():
    P_ = sin_((0, 2, -2, 0))*cos_((1, 0, -2, 0))                   # sin((1-t)A) cos(beta/2 - tA)
    Q_ = sin_((0, 0, 2, 0))*cos_((1, 2, -2, 0))*sin_((1, -2, 0, 0))  # sin(tA) cos(beta/2+A-tA) sin(beta/2-A)
    bad = []
    for par in (0, 1):
        pb = (0, 0, 0, 8) if par == 0 else (4, 0, 0, 8)              # 2 p beta
        pbh = (0, 0, 0, 4) if par == 0 else (2, 0, 0, 4)             # p beta
        sgn = 1 if par == 0 else -1                                  # (-1)^p
        Sp = sin_(add(pbh, (1, 2, 0, 0)))                            # sin(p beta + beta/2 + A)
        # J1: V1 vertex of the exit junction, at 2 p beta + 2A
        if not is_zero(cb2*tipF(par, sC, add(pb, (0, 4, 0, 0))) - (sgn*P_*Sp + Q_)):
            bad.append(('J1', par))
        # J2: V2 vertex of the exit junction, at 2 p beta + 3A
        if par == 0:
            cf = cos_(add((0, 0, 0, 2), (1, 0, 0, 0)))*sin_(add((0, 0, 0, 2), (0, 2, 0, 0)))    # cos((p+1)b/2) sin(pb/2 + A)
        else:
            cf = sin_((2, 0, 0, 2))*cos_((1, 2, 0, 2))                                        # sin((p+1)b/2) cos(pb/2 + A)
        if not is_zero(cb2*tipF(par, sB, add(pb, (0, 6, 0, 0))) - (-2*sgn)*cU*sin_((0, 0, 2, 0))*cf):
            bad.append(('J2', par))
        # R2: interior V2 of a rising V0-star, at 2 p beta + A
        Ip = sin_((0, 2, -2, 0))*sin_(add(pbh, (1, 0, 0, 0))) + sin_((0, 2, 0, 0))*sin_(add(pbh, (-1, 0, 2, 0))) \
            + sgn*sin_((0, 0, 2, 0))*sin_((1, -2, 0, 0))
        if not is_zero(cb2*tipF(par, sB, add(pb, (0, 2, 0, 0))) - sgn*cU*Ip):
            bad.append(('R2', par))
        # F2: interior V2 of a falling V0-star, at -A
        Fp = sin_((0, 2, -2, 0))*sin_(add(pbh, (-1, 4, 0, 0))) + sin_((0, 2, 0, 0))*sin_(add(pbh, (1, 4, -2, 0))) \
            - sgn*sin_((0, 0, 2, 0))*sin_((1, -2, 0, 0))
        if not is_zero(cb2*tipF(par, sB, (0, -2, 0, 0)) + sgn*cU*Fp):
            bad.append(('F2', par))
        # centre: F_p(0) = sB/cos(b/2) [sin(b/2 - A) sin tA + (-1)^{p+1} cos(b/2 - tA) cos(p beta + A)]
        cen = F_even0 if par == 0 else F_odd0
        if not is_zero(cb2*cen - sB*(sin_((1, -2, 0, 0))*sin_((0, 0, 2, 0)) - sgn*cos_((1, 0, -2, 0))*cos_(add(pbh, (0, 2, 0, 0))))):
            bad.append(('centre', par))
    # first V0-star: J1(0) (tip 0, V1) and J2(0) (tip -A, V2), even-p closed forms at p = 0
    if not is_zero(cb2*F_first(sC, (0, 0, 0, 0)) - (P_*sin_((1, 2, 0, 0)) + Q_)):
        bad.append(('J1 first',))
    if not is_zero(cb2*F_first(sB, (0, -2, 0, 0)) - (-2)*cU*sin_((0, 0, 2, 0))*cos_((1, 0, 0, 0))*sin_((0, 2, 0, 0))):
        bad.append(('J2 first',))
    if not is_zero(F_first(sB, (0, 2, 0, 0))):          # S = tip A on the line
        bad.append(('S on line',))
    # last V0-star p = q-1 = 2h (even), z^4 = -x^-2
    if not is_zero(tipF(0, sB, (-4, 6, 0, 0)), last=True):                    # E = tip -(4m-1)A = -2beta + 3A: closure
        bad.append(('closure',))
    if not is_zero(tipF(0, sC, (-4, 8, 0, 0)) + sin_((0, 2, 0, 0))*sin_((3, -2, -4, 0)), last=True):       # L1
        bad.append(('L1',))
    if not is_zero(tipF(0, sB, (-4, 10, 0, 0)) + 2*sin_((0, 2, 0, 0))*cU*sin_((2, -4, -2, 0)), last=True):  # L2
        bad.append(('L2',))
    # last V0-star V2 at -A (F2 at p = q-1): factored, all factors positive
    lastF2 = sin_((0, 2, -2, 0))*cU*sin_((1, -2, 0, 0))*cos_((1, -2, 0, 0))
    # exact constant: evaluate the ratio at a generic point satisfying z^4 = -x^-2, then test the identity exactly
    vals = {x: sp.exp(I*sp.Rational(3, 10)), y: sp.exp(I*sp.Rational(1, 13)), s: sp.exp(I*sp.Rational(1, 29))}
    vals[z] = sp.exp(I*(sp.pi - sp.Rational(3, 5))/4)
    r = sp.N((tipF(0, sB, (0, -2, 0, 0))/lastF2).subs(vals), 30)
    c = sp.nsimplify(sp.re(r), tolerance=1e-20, rational=True)
    if abs(sp.im(r)) > 1e-20 or not is_zero(tipF(0, sB, (0, -2, 0, 0)) - c*lastF2, last=True):
        bad.append(('last F2 factor', r))
    else:
        print('  last-V0-star V2(-A): F = (%s) * lastF2' % c)
    # (*) at t = 2/3:  LHS - RHS = sin(A/3) [2 cos(3u - 4A/3) sin(u - A/3) + sin(5A/3) + sin A]/2, u = beta/2
    u, a = sp.symbols('u a', positive=True)
    t = sp.Rational(2, 3)
    lhs = sp.sin(3*u - a)*sp.sin((1 - t)*a)*sp.cos(u - t*a) - sp.sin(t*a)*sp.cos(u + a - t*a)*sp.sin(u - a)
    rhs = sp.sin(a/3)*(2*sp.cos(3*u - 4*a/3)*sp.sin(u - a/3) + sp.sin(5*a/3) + sp.sin(a))/2
    if sp.simplify(sp.expand_trig(sp.expand((lhs - rhs).rewrite(sp.exp)))) != 0:
        bad.append(('star 2/3',))
    return bad


# ---------------------------------------------------------------- part (b): the m prongs W_{q,m,j}, j = 0..m-1  (p >= m)
def word_j(q, m, j):
    """runs: first 2j+1, rising r_p + 2j, falling f_p - 2j, last 4m-1-2j (j = 0 is W_{q,m})."""
    h = (q - 1)//2
    parts = [alt('0', 2*j + 1)]
    parts += [alt('2', 4*m + 4 + i*(4*m + 2) + 2*j) for i in range(h)]
    parts += [alt('0', 8*m + 2 + i*(4*m + 2) - 2*j) for i in range(h - 2, -1, -1)]
    parts += [alt('0', 4*m - 1 - 2*j)]
    return '1'.join(parts)


def numeric_checks_j():
    mp.mp.dps = 40
    bad = []
    for q, m in [(3, 3), (5, 2), (5, 4), (7, 3), (9, 2), (11, 3)]:
        n, K, g0 = consts(q, m)
        A = mp.pi/n
        tu = t_up(q, m)
        words = [word_j(q, m, j) for j in range(m)]
        allw = set(words) | set(w[::-1] for w in words)
        if len(allw) != 2*m:
            bad.append(('not 2m distinct words', q, m))
        for j, w in enumerate(words):
            for t in [mp.mpf('0.01'), mp.mpf(1)/3, mp.mpf(2)/3, tu - mp.mpf('1e-6')]:
                B = (K + t)*A
                cps = unfold(A, B, w)
                S, E = cps[0][2], cps[-1][2]
                dvec = E - S
                f = lambda zz: mp.im((zz - S)*mp.conj(dvec))
                ok = all(f(cps[k][SIDES[int(ch)][0]])*f(cps[k][SIDES[int(ch)][1]]) < 0 for k, ch in enumerate(w))
                c = (mp.pi/2 + B - mp.arg(dvec))/A
                if not ok or abs(c - (g0 + j)) > mp.mpf('1e-25'):
                    bad.append(('prong j', q, m, j, t, ok, c))
            B = (K - mp.mpf('1e-6'))*A        # control: below the band
            cps = unfold(A, B, w)
            S, E = cps[0][2], cps[-1][2]
            dvec = E - S
            f = lambda zz: mp.im((zz - S)*mp.conj(dvec))
            if all(f(cps[k][SIDES[int(ch)][0]])*f(cps[k][SIDES[int(ch)][1]]) < 0 for k, ch in enumerate(w)):
                bad.append(('control j did not fire', q, m, j))
    return bad


def symbolic_checks_j():
    """the same closed forms with the prong index j free: v = e^{ijA/2} (angle vectors get a 5th entry)."""
    v = sp.Symbol('v', nonzero=True)

    def E5(w):
        return x**w[0]*y**w[1]*s**w[2]*z**w[3]*v**w[4]

    def c5(w):
        return (E5(w) + E5(tuple(-a for a in w)))/2

    def s5(w):
        return (E5(w) - E5(tuple(-a for a in w)))/(2*I)

    sB5, sC5, cb5, cU5 = c5((1, 2, -2, 0, 0)), c5((1, 0, -2, 0, 0)), c5((1, 0, 0, 0, 0)), c5((1, 2, -2, 0, 0))
    Fe = 2*sB5*sC5*s5((0, 0, 0, 2, 0))*s5((0, 2, 0, 2, 2))/cb5 - sB5*c5((0, 2, -2, 0, 2))
    Fo = Fe + 2*sB5*sC5*c5((1, 2, 0, 4, 2))

    def tip5(par, R, phi):
        if par == 0:
            return Fe + R*c5(add((0, 4, -2, 4, 2), tuple(-a for a in phi)))
        return Fo - R*c5(add(phi, (-2, -4, 2, -4, -2)))

    def first5(R, phi):
        return R*c5(add(phi, (0, 0, -2, 0, 2))) - sB5*c5((0, 2, -2, 0, 2))

    P5 = s5((0, 2, -2, 0, 0))*c5((1, 0, -2, 0, 0))
    Q5 = s5((0, 0, 2, 0, 0))*cU5*s5((1, -2, 0, 0, -2))
    tests = []
    for par in (0, 1):
        pb = (0, 0, 0, 8, 0) if par == 0 else (4, 0, 0, 8, 0)
        pbh = (0, 0, 0, 4, 0) if par == 0 else (2, 0, 0, 4, 0)
        j2 = (0, 0, 0, 0, 4)
        sg = 1 if par == 0 else -1
        Sj = s5(add(pbh, (1, 2, 0, 0, 2)))
        tests.append(('J1', par, cb5*tip5(par, sC5, add(pb, (0, 4, 0, 0, 0), j2)), sg*P5*Sj + Q5, False))
        cf = c5(add((0, 0, 0, 2, 0), (1, 0, 0, 0, 0)))*s5(add((0, 0, 0, 2, 0), (0, 2, 0, 0, 2))) if par == 0 \
            else s5((2, 0, 0, 2, 0))*c5((1, 2, 0, 2, 2))
        tests.append(('J2', par, cb5*tip5(par, sB5, add(pb, (0, 6, 0, 0, 0), j2)), -2*sg*cU5*s5((0, 0, 2, 0, 0))*cf, False))
        Ij = s5((0, 2, -2, 0, 0))*s5(add(pbh, (1, 0, 0, 0, 2))) + s5((0, 2, 0, 0, 0))*s5(add(pbh, (-1, 0, 2, 0, 2))) \
            + sg*s5((0, 0, 2, 0, 0))*s5((1, -2, 0, 0, -2))
        tests.append(('R2', par, cb5*tip5(par, sB5, add(pb, (0, 2, 0, 0, 0), j2)), sg*cU5*Ij, False))
        Fj = s5((0, 2, -2, 0, 0))*s5(add(pbh, (-1, 4, 0, 0, 2))) + s5((0, 2, 0, 0, 0))*s5(add(pbh, (1, 4, -2, 0, 2))) \
            - sg*s5((0, 0, 2, 0, 0))*s5((1, -2, 0, 0, -2))
        tests.append(('F2', par, cb5*tip5(par, sB5, (0, -2, 0, 0, 0)), -sg*cU5*Fj, False))
    tests += [
        ('first S', '-', first5(sB5, (0, 2, 0, 0, 0)), 0, False),
        ('first J1(0) @-2jA', '-', cb5*first5(sC5, (0, 0, 0, 0, -4)), P5*s5((1, 2, 0, 0, 2)) + Q5, False),
        ('first J2(0) @-(2j+1)A', '-', cb5*first5(sB5, (0, -2, 0, 0, -4)), -2*cU5*s5((0, 0, 2, 0, 0))*c5((1, 0, 0, 0, 0))*s5((0, 2, 0, 0, 2)), False),
        ('first V1 0', '-', first5(sC5, (0, 0, 0, 0, 0)), s5((0, 2, 0, 0, 0))*s5((1, 2, -4, 0, 2)), False),
        ('first V2 -A', '-', first5(sB5, (0, -2, 0, 0, 0)), 2*s5((0, 2, 0, 0, 0))*cU5*s5((0, 0, -2, 0, 2)), False),
        ('first V2 -(2j-1)A', '-', first5(sB5, (0, 2, 0, 0, -4)), 2*cU5*s5((0, 2, -2, 0, 0))*s5((0, 0, 0, 0, 2)), False),
        ('closure', '-', tip5(0, sB5, (-4, 6, 0, 0, 4)), 0, True),
        ('last L1', '-', tip5(0, sC5, (-4, 8, 0, 0, 4)), -s5((0, 2, 0, 0, 0))*s5((3, -2, -4, 0, -2)), True),
        ('last L2', '-', tip5(0, sB5, (-4, 10, 0, 0, 4)), -2*s5((0, 2, 0, 0, 0))*cU5*s5((2, -4, -2, 0, -2)), True),
        ('last V2 -A', '-', tip5(0, sB5, (0, -2, 0, 0, 0)), -2*s5((0, 2, -2, 0, 0))*cU5*s5((2, -4, 0, 0, -2)), True)]
    bad = []
    for name, par, lhs, rhs, last in tests:
        if not is_zero(lhs - rhs, last=last):
            bad.append((name, par))
    return bad


def bridge_check_j():
    """audit W4: the j-general frame formulas against a direct unfolding of W_{q,m,j} (vertex signed distances), plus an
    upper control at t_up(j) + 1e-6 for the binding prong j = 0 handled in numeric_checks."""
    mp.mp.dps = 30
    bad = []
    for q, m, j in [(5, 3, 2), (7, 2, 1), (9, 3, 1)]:
        n, K, g0 = consts(q, m)
        A, b = mp.pi/n, mp.pi/q
        t = mp.mpf('0.43')
        B = (K + t)*A
        sB_, sC_ = mp.sin(B), mp.sin(B + A)
        psi = -mp.pi/2 + (t - j)*A
        w = word_j(q, m, j)
        cps = unfold(A, B, w)
        S = cps[0][2]
        Fdir = lambda zz: mp.im((zz - S)*mp.expj(-psi))
        starts = [0] + [k + 1 for k, ch in enumerate(w) if ch == '1']
        x_ = (1 + j)*A
        for p, k0 in enumerate(starts):
            O = cps[k0][0]
            if p == 0:
                cf = -sB_*mp.cos((1 - t + j)*A)
            else:
                cf = sB_/mp.cos(b/2)*(mp.sin(b/2 - x_)*mp.sin(t*A) + (-1)**(p + 1)*mp.cos(b/2 - t*A)*mp.cos(p*b + x_))
            if abs(Fdir(O) - cf) > mp.mpf('1e-25'):
                bad.append(('centre j', q, m, j, p))
        # V0-star 0 tips: F_0(R e^{i phi}) = R cos(phi - (t - j)A) - sB cos((1 - t + j)A)
        for ray in (0, -1, -2*j, -(2*j + 1)):
            R = sC_ if ray % 2 == 0 else sB_
            P = R*mp.expj(ray*A)
            if abs(Fdir(P) - (R*mp.cos(ray*A - (t - j)*A) - sB_*mp.cos((1 - t + j)*A))) > mp.mpf('1e-25'):
                bad.append(('V0-star-0 tip j', q, m, j, ray))
        if abs(Fdir(cps[-1][2])) > mp.mpf('1e-25'):
            bad.append(('closure j', q, m, j))
    return bad


def bridge_check():
    """the symbolic tip formula (S) and the numeric closed form (N) are the same function: compare at sample points."""
    mp.mp.dps = 30
    bad = []
    for q, m, t in [(9, 2, mp.mpf('0.41')), (7, 1, mp.mpf('0.77'))]:
        n, K, g0 = consts(q, m)
        A, b = mp.pi/n, mp.pi/q
        for p in range(1, q):
            k = p//2
            vals = {x: mp.expj(b/2), y: mp.expj(A/2), s: mp.expj(t*A/2), z: mp.expj(k*b/2)}
            for ray in (0, -1, 3, 7):
                R = sC if ray % 2 == 0 else sB
                symb = complex(sp.N(tipF(p % 2, R, (0, 2*ray, 0, 0)).subs({kk: sp.Float(str(mp.re(v)), 30) + I*sp.Float(str(mp.im(v)), 30) for kk, v in vals.items()}), 25))
                num = closed_F(q, m, t, p, ray)
                if abs(symb.imag) > 1e-20 or abs(symb.real - float(num)) > 1e-12:
                    bad.append((q, m, p, ray, symb, num))
    return bad


if __name__ == '__main__':
    bb = bridge_check()
    print('bridge (S) = (N):', 'OK' if not bb else bb, flush=True)
    bs = symbolic_checks()
    print('symbolic:', 'OK' if not bs else bs, flush=True)
    bn = numeric_checks()
    print('numeric:', 'OK' if not bn else bn, flush=True)
    bsj = symbolic_checks_j()
    print('symbolic, prong index j free:', 'OK' if not bsj else bsj, flush=True)
    bnj = numeric_checks_j()
    print('numeric, all m prongs:', 'OK' if not bnj else bnj, flush=True)
    bbj = bridge_check_j()
    print('bridge, prong j (frame = closed form):', 'OK' if not bbj else bbj, flush=True)
    print('ALL OK' if not (bs or bn or bb or bsj or bnj or bbj) else 'FAIL')
