"""Section 9.2, Proposition 9.10: the corner witness on the leg-swap pair 8/15 eps 0, 7/15 eps 1.

Statements checked:
  * (the witness, on the chain) at each row exactly four walks run from a zero of the base to a
    zero without meeting Fix(iota) (no kite-0 transit, no interface-m crossing), each in 5 exact
    steps; they are the two ends of two edges, the edges are iota-images of each other (the one
    swapped pair, p = 1), and in one half of the necklace (kites 1..(Q-1)/2) they are two prongs
    that are exact time-reverses;
  * the four corner arrivals the paper counts are ONE saddle connection -- 1 object x 2 ends x 2
    leg-swap rows: the leg swap (Z_O <-> Z_A, kites fixed) carries one row's edges onto the other's;
  * (the witness, unfolded) at 8/15 the prong from O in direction 12 degrees -- word
    L1 H L2 H L1 H L1 L2 H L1 L2 H L2 -- arrives EXACTLY at the vertex A at its 14th hit, every
    earlier hit strictly inside its side (exact signs); the other two prongs of the class from O,
    at 24 and 36 degrees, reach no vertex within 60 hits (the controls);
  * the arrival consumes exactly one relation: the arrival identity (Section 6.4) along that word,
    F_A = sin((M + eps P) theta) + sum over the L1 letters of [sin((m_i + P) theta) + sin((m_i - P) theta)],
    is, as a formal sum of sines, 2 sin 12 + 2 sin 24 + sin 36 - sin 60 - sin 72, and that relation
    holds exactly in Q(zeta_60);
  * transported to Z[C_15] (zeta_30 = -zeta_15^8) the relation is nonzero at a primitive cube root
    and at a primitive fifth root of unity, and lies in the span of the 3 | Q and 5 | Q generator
    families (Remark 6.14) together and of neither alone (exact ranks over Q).
The chain walks are `lib/necklace.py`'s; the word is proposed by the 60-digit tracer and decided by
`lib/triangle.py`'s exact unfolding.  Seconds.
"""
from mpmath import mp, mpf, pi
from sympy import Matrix, Poly, cyclotomic_poly, symbols

from common import check, done
from necklace import readout
from tracer import Billiard
from triangle import RightTriangle


def corner_walks(P, Q, eps):
    """The walks from a zero to a zero that meet Fix(iota) nowhere."""
    o = readout(P, Q, eps, step_cap=200_000, with_seps=True)
    out = []
    for s in o['seps']:
        a, b = s['ends']
        if s['in_fix'] or b == 'CAP' or a[0] == 'R' or b[0] == 'R':
            continue
        if s['t0'] == 0 and s['onm'] == 0:
            out.append(s)
    return o, out


def fold(k, Q):
    """sin(k pi / 2Q) as (sign, j) with 0 <= j <= Q: sin is odd, 4Q-periodic and sin(x) = sin(pi - x)."""
    k %= 4 * Q
    sg = 1
    if k > 2 * Q:
        k, sg = 4 * Q - k, -1
    if k > Q:
        k = 2 * Q - k
    return sg, k


def sine_sum(terms, Q):
    out = {}
    for c, k in terms:
        sg, j = fold(k, Q)
        if j:
            out[j] = out.get(j, 0) + sg * c
    return {j: c for j, c in out.items() if c}


def main():
    print('Section 9.2 -- the corner witness (Proposition 9.10)\n')
    rows = {}
    for P, Q, eps in [(7, 15, 1), (8, 15, 0)]:
        o, W = corner_walks(P, Q, eps)
        edges = {frozenset(map(tuple, s['ends'])) for s in W}
        iota = lambda e: frozenset((a, (-k) % Q) for a, k in e)
        half = [s for s in W if 1 <= s['ends'][0][1] <= (Q - 1) // 2]
        rev = len(half) == 2 and half[0]['ends'] == half[1]['ends'][::-1]
        ok = (len(W) == 4 and all(s['steps'] == 5 for s in W) and len(edges) == 2
              and iota(next(iter(edges))) in edges and iota(next(iter(edges))) != next(iter(edges))
              and o['swapped'] == 1 and rev)
        check(f'{P}/{Q} eps {eps}: {len(W)} corner walks off Fix(iota), steps '
              f'{sorted(s["steps"] for s in W)}, edges {sorted(sorted(e) for e in edges)}, '
              f'iota-swapped pair, p = {o["swapped"]}; one half: 2 prongs, time-reverses: {rev}', ok)
        rows[(P, Q, eps)] = edges
    swap = lambda e: frozenset(({'ZO': 'ZA', 'ZA': 'ZO'}[a], k) for a, k in e)
    check('leg swap Z_O <-> Z_A carries 7/15 eps 1\'s edges onto 8/15 eps 0\'s: 4 arrivals = '
          '1 connection x 2 ends x 2 rows',
          {swap(e) for e in rows[(7, 15, 1)]} == rows[(8, 15, 0)])

    # -- the unfolded prong at 8/15 ---------------------------------------------------------
    P, Q = 8, 15
    mp.dps = 60
    B, T = Billiard(P * pi / (2 * Q)), RightTriangle(P, Q)
    F = T.F
    verts = {'A': (T.L, F.one), 'R': (T.L, F.zero), 'O': (F.zero, F.zero)}
    arrivals = {}
    for u0 in (2, 4, 6):                        # the class's directions strictly inside the angle at O
        word = ''.join(s for s, _, _ in B.orbit(mpf(0), mpf(0), u0 * pi / (2 * Q), 60))
        steps = T.unfold(word, u0, start=T.launch_point(F.zero, F.zero))
        hit, inside = None, True
        for i, st in enumerate(steps):
            Ax, Ay = st[4]
            hit = next((v for v, (x, y) in verts.items() if (Ax - x).is_zero() and (Ay - y).is_zero()),
                       None)
            if hit:
                hit = (i + 1, hit)
                break
            inside &= all(f0.sign() > 0 for _, _, f0, _ in T.constraints(st))
        arrivals[u0] = (word[:hit[0]] if hit else word, hit, inside)
    w2, hit2, in2 = arrivals[2]
    check(f'8/15, prong from O at 12 deg: word {" ".join(w2[:-1])} + corner, exact arrival at '
          f'{hit2}, every earlier hit strictly inside its side: {in2}',
          hit2 == (14, 'A') and in2 and w2[:13] == '1h2h1h12h12h2')
    check(f'controls: prongs from O at 24, 36 deg reach no vertex within 60 hits: '
          f'{arrivals[4][1]}, {arrivals[6][1]}', arrivals[4][1] is None and arrivals[6][1] is None)

    # -- the relation it consumes ------------------------------------------------------------
    d = {'2': 0, 'h': P, '1': Q}
    word, u0 = w2[:13], 2
    D, Dk = [0], 0
    for j, s in enumerate(word, 1):
        Dk += (-1) ** (j - 1) * d[s]
        D.append(Dk)
    n = len(word)
    M, e = 2 * D[n] - u0, (-1) ** n
    terms = [(1, M + e * P)]
    for i, s in enumerate(word, 1):
        if s == '1':
            m = 2 * D[i - 1] - u0
            terms += [(1, m + P), (1, m - P)]
    FA = sine_sum(terms, Q)
    R = sine_sum([(2, 2), (2, 4), (1, 6), (-1, 10), (-1, 12)], Q)      # units of 6 degrees
    check(f'F_A along the word = {FA} (sin(6j deg): coeff) = 2 sin 12 + 2 sin 24 + sin 36 - sin 60 '
          f'- sin 72', FA == R)
    val = sum((F.sin(k) * c for k, c in R.items()), F.zero)
    check('the relation holds exactly in Q(zeta_60)', val.is_zero())

    # -- its place in Remark 6.14's generator families ---------------------------------------
    t = symbols('t')
    r = [0] * 15                                  # zeta_30^j -> (-1)^j t^(8j mod 15)
    for k, c in R.items():                        # sin(6k deg) = (zeta_30^(k/2) - zeta_30^(-k/2)) / 2i
        for jj, sg in ((k // 2, 1), (-(k // 2), -1)):
            r[(8 * jj) % 15] += sg * c * (-1) ** (jj % 2)
    rp = Poly(sum(c * t ** i for i, c in enumerate(r)), t)
    z15 = Poly(cyclotomic_poly(15, t), t)
    at = {q: not rp.rem(Poly(cyclotomic_poly(q, t), t)).is_zero for q in (3, 5)}
    check(f'transported to Z[C_15]: vanishes at a primitive 15th root ({rp.rem(z15).is_zero}), '
          f'nonzero at a primitive cube root ({at[3]}) and fifth root ({at[5]})',
          rp.rem(z15).is_zero and at[3] and at[5])
    vec = lambda g, m: [g[(i - m) % 15] for i in range(15)]
    g3 = [1 if i % 5 == 0 else 0 for i in range(15)]           # 1 + t^5 + t^10
    g5 = [1 if i % 3 == 0 else 0 for i in range(15)]           # 1 + t^3 + ... + t^12
    fam3 = [vec(g3, m) for m in range(15)]
    fam5 = [vec(g5, m) for m in range(15)]
    rk = lambda rows_: Matrix(rows_).rank()
    alone3, alone5 = rk(fam3 + [r]) > rk(fam3), rk(fam5 + [r]) > rk(fam5)
    both = rk(fam3 + fam5 + [r]) == rk(fam3 + fam5)
    check(f'in the span of the 3 | Q and 5 | Q families together ({both}), of neither alone '
          f'({alone3}, {alone5})', both and alone3 and alone5)
    done()


if __name__ == '__main__':
    main()
