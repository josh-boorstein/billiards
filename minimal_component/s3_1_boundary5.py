"""Section 3.1, "Why boundary 5 is different", and Figure 2's data.

Statements checked:
  * Carry the transversal direction along as theta = 14 alpha - 3 pi (equal to 11pi/15 at
    T0).  Every return word is a translation at every alpha (its linear part is the
    identity), so each word's translation is an explicit function of alpha.
  * Across boundaries 3, 4 and 6 (the R-grazes inside a return path) the two translations
    agree IDENTICALLY in alpha.
  * Across boundary 5 their difference is -(cot alpha / cos theta) times an integer sine
    sum equal, up to orientation, to
        V(alpha) = 2(sin 14a - sin 12a - sin 10a + 2 sin 8a - sin 4a),
    which vanishes at alpha = 4pi/15 -- equivalently sin 12 + 2 sin 24 + sin 36 =
    sin 48 + sin 60 (degrees), checked exactly in Q(zeta_60) -- with nonzero derivative.
  * (cot alpha / |cos theta|) V'(4pi/15) = 1.3456 x 85.21 = 114.7 per radian; and a
    60-digit trace of the deformed flow agrees with the exact jumps to 1e-59 and measures
    the jump growing at 114.7 per radian.
  * The other continuation, theta = -16 alpha + 5 pi, also equals 11pi/15 at T0, and its
    boundary-5 sine sum also vanishes at T0 with nonzero derivative.
  * Off T0 the exchange has more than four intervals, on both sides, but not in the same
    way.  For alpha > 4pi/15 boundary 5 opens into a gap exactly as wide as the jump,
    filled by new, longer return words with other translations.  For alpha < 4pi/15 the
    two words either side of boundary 5 still abut -- the break is a bare jump -- and the
    new words appear elsewhere on the transversal (at its right end, and near boundary 2).
    The paper's sentence "it opens into a gap, as wide as the jump" describes the first
    side only.

How: the translations as functions of alpha come from the formal development of each
word (`lib/develop.py`), which composes the three reflections symbolically and never
substitutes an angle; the trace is the independent 60-digit tracer (`lib/tracer.py`).
"""
from mpmath import mp, mpf, pi, sin, cos, cot, nstr

from t0 import F, return_words, check, done
from develop import sine_sum, offset, is_translation
from tracer import Billiard

mp.dps = 60
A0 = 4 * pi / 15
V_PAPER = {14: 2, 12: -2, 10: -2, 8: 4, 4: -2}          # V = sum C_n sin(n alpha)


def V(C, a):
    return sum(k * sin(n * a) for n, k in C.items())


def dV(C, a):
    return sum(k * n * cos(n * a) for n, k in C.items())


def main():
    print('Section 3.1 -- why boundary 5 is different\n')
    words = return_words()
    check('every return word is a translation at every alpha',
          all(is_translation(w) for w in words))
    I3 = words[2:]                               # boundaries 3..6 lie between these
    sums = [sine_sum(a, b, 14, -3) for a, b in zip(I3, I3[1:])]
    for b, S in zip((3, 4, 6), (sums[0], sums[1], sums[3])):
        check(f'boundary {b}: the two translations agree identically in alpha', S == {})
    S5 = sums[2]
    neg = {n: -k for n, k in S5.items()}
    check(f'boundary 5: sine sum {S5} = +/- V(alpha)', S5 == V_PAPER or neg == V_PAPER)

    # V(4pi/15) = 0 exactly: with alpha = 48 degrees, V/2 reduces to
    # sin 12 + 2 sin 24 + sin 36 - sin 48 - sin 60 (degrees); F.sin(k) = sin(6k degrees)
    ident = F.sin(2) + 2 * F.sin(4) + F.sin(6) - F.sin(8) - F.sin(10)
    check('sin 12 + 2 sin 24 + sin 36 = sin 48 + sin 60, exactly in Q(zeta_60)', ident.is_zero())
    check('and V(4pi/15) = 0 numerically too', abs(V(V_PAPER, A0)) < mpf(10) ** -55)
    d0 = dV(V_PAPER, A0)
    fac = cot(A0) / abs(cos(11 * pi / 15))
    print(f'\n  V\'(4pi/15) = {nstr(d0, 8)};  cot(alpha)/|cos theta| = {nstr(fac, 8)};'
          f'  product = {nstr(fac * d0, 8)} per radian')
    check('V\'(4pi/15) = 85.21..., nonzero', nstr(abs(d0), 4) == '85.21')
    check('cot(alpha)/|cos theta| = 1.3456', nstr(fac, 5) == '1.3456')
    check('predicted jump rate 114.7 per radian', nstr(abs(fac * d0), 4) == '114.7')

    # the other continuation
    check('14 alpha - 3 pi = -16 alpha + 5 pi = 11 pi/15 at T0',
          abs(14 * A0 - 3 * pi - 11 * pi / 15) < mpf(10) ** -55 and
          abs(-16 * A0 + 5 * pi - 11 * pi / 15) < mpf(10) ** -55)
    S5b = sine_sum(I3[2], I3[3], -16, 5)
    print(f'  theta = -16 alpha + 5 pi: boundary-5 sine sum {S5b}, '
          f'value at T0 {nstr(V(S5b, A0), 3)}, derivative {nstr(dV(S5b, A0), 6)}')
    check('other continuation: vanishes at T0, nonzero derivative',
          abs(V(S5b, A0)) < mpf(10) ** -55 and abs(dV(S5b, A0)) > 1)
    others_b = [sine_sum(a, b, -16, 5) for a, b in zip(I3, I3[1:])]
    check('other continuation: boundaries 3, 4, 6 identically zero there too',
          others_b[0] == others_b[1] == others_b[3] == {})

    # the 60-digit trace of the deformed flow
    th = lambda a: 14 * a - 3 * pi
    worst, meas = mpf(0), []
    for k in (-3, -2, -1, 1, 2, 3):
        d = k * mpf('1e-4')
        a = A0 + d
        B = Billiard(a)
        cw = {}
        for i in range(1, 1500):
            s = mpf(i) / 1500
            r = B.first_return(s, th(a))
            if r is not None and r[1] in I3 and r[1] not in cw:
                cw[r[1]] = r[0] - s
        for u, v in zip(I3, I3[1:]):
            exact = offset(v, a, th(a)) - offset(u, a, th(a))
            worst = max(worst, abs((cw[v] - cw[u]) - exact))
        meas.append((d, cw[I3[3]] - cw[I3[2]]))
    rate = (meas[-1][1] - meas[0][1]) / (meas[-1][0] - meas[0][0])
    print(f'\n  trace vs exact jumps, 24 points: worst |diff| = {nstr(worst, 3)}; '
          f'measured boundary-5 rate {nstr(abs(rate), 5)} per radian')
    check('the trace agrees with the exact jumps to 1e-59', worst < mpf(10) ** -58)
    check('measured rate 114.7 per radian', nstr(abs(rate), 4) == '114.7')

    # off T0: count the exchange's intervals by a scan, and look at boundary 5 itself
    for d in (mpf('1e-4'), mpf('-1e-4')):
        a = A0 + d
        B = Billiard(a)
        word = lambda s: (B.first_return(s, th(a), cap=3000) or (None, None))[1]
        grid = [mpf(i) / 2000 for i in range(1, 2000)]
        tags = [word(s) for s in grid]
        # maximal runs of one translation: known words of equal translation merge
        tr = {w: offset(w, a, th(a)) for w in set(tags) if w is not None and is_translation(w)}
        runs = []
        for t in tags:
            key = ('?', t) if t not in tr else ('c', nstr(tr[t], 25))
            if not runs or runs[-1] != key:
                runs.append(key)
        n_int = len(runs)
        new = sorted({len(t) for t in tags if t is not None and t not in words})
        lo_in = max(s for s, t in zip(grid, tags) if t == I3[2])
        hi_in = min(s for s, t in zip(grid, tags) if t == I3[3] and s > lo_in)
        between = [t for s, t in zip(grid, tags) if lo_in < s < hi_in]

        def edge(s_in, s_out, w_in):
            for _ in range(150):
                m = (s_in + s_out) / 2
                if word(m) == w_in:
                    s_in = m
                else:
                    s_out = m
            return s_in
        jump = abs(offset(I3[3], a, th(a)) - offset(I3[2], a, th(a)))
        print(f'  alpha - 4pi/15 = {nstr(d, 2)}: {n_int} intervals in the scan; new words of '
              f'lengths {new}')
        check('  more than four intervals', n_int > 4 and bool(new))
        if d > 0:
            mid = [s for s, t in zip(grid, tags) if lo_in < s < hi_in][0]
            e1 = edge(lo_in, mid, I3[2])
            e2 = edge(hi_in, mid, I3[3])
            print(f'    boundary 5 opens into a gap of width {nstr(e2 - e1, 12)}; '
                  f'|jump| = {nstr(jump, 12)}')
            check('  the gap at boundary 5 is as wide as the jump (to 1e-30)',
                  bool(between) and abs((e2 - e1) - jump) < mpf(10) ** -30)
            got = [t for t in between if t is not None]
            check(f'  and is filled by longer words than any of the seven '
                  f'({len(between) - len(got)} of {len(between)} scan points do not return '
                  f'within 3000 bounces)',
                  bool(got) and all(t not in words and len(t) > 66 for t in got))
        else:
            e1 = edge(lo_in, hi_in, I3[2])
            e2 = edge(hi_in, lo_in, I3[3])
            print(f'    words 5 and 6 abut: boundary {nstr(e1, 15)} / {nstr(e2, 15)}; '
                  f'|jump| = {nstr(jump, 12)}')
            check('  on this side boundary 5 is a bare jump (no gap)',
                  not between and abs(e2 - e1) < mpf(10) ** -40 and jump > mpf(10) ** -6)
    done()


if __name__ == '__main__':
    main()
