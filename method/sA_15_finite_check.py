"""Appendix A.5, Lemma A.15 (bad gates): the finite check its proof leaves.

Lemma A.15: let Q = 2^a 3^b >= 81, P odd and coprime to Q with 2 <= P <= Q - 2, and varpi in {0, 1}. Every
good class K(e) with e = varpi (mod 2) has a bad gate, except the class of e = N - 1 when Q = 3^(2k),
P = (Q + 1)/2 and varpi = 0. Its proof argues by regimes of P and leaves finitely many cases, all at the
eight Q = 2^a 3^b (b >= 1) with 80 < Q < 280 -- 81, 96, 108, 144, 162, 192, 216, 243 -- where they are
checked directly. This script is that direct check.

The objects, as A.5 defines them (exact rational arithmetic throughout, no floating point):
  * N = Q/3, D = 2Q/P, P = 2 eps (mod 3); the class of 0 < e < N is K(e) = {e, 2N - e, 2N + e}.
  * K(e) is GOOD iff frac(e/D) = frac(eP/2Q) lies in (1/2, 5/6) (eps = +1), resp. (1/6, 1/2) (eps = -1).
  * For 0 < n < P/2 the GATE z_n of the zero nD is the integer = varpi (mod 2) with |z_n - nD| < 1. A frame is
    good if its class is; the frames N and 2N are not good. A BAD GATE of K(e) is a gate with nD in
    (e, 2N + e) that is not good.
  * The cases the proof leaves (its regimes (iv) and (v)): Q/3 < P < 2Q/3, varpi not = N (mod 2), and either
    P outside the window [7Q/15, 5Q/9] with Q < 102 (regime (iv) needs Q >= 102), or P in the window with
    j = |Q - 2P| >= 5 and P < 126 (regime (v) covers P >= 126).

Statements checked:
  * the eight Q are exactly the Q = 2^a 3^b, b >= 1, with 80 < Q < 280;
  * the goodness rule is Lemma 5.2's at p = 3: at every class of every centre checked, the chain signs
    chi = (c(2N - e), c(e), -c(2N + e)), c(d) = sgn cos(d alpha) read from dP mod 4Q, are constant exactly
    when frac(eP/2Q) lies in the window above (an independent computation of the same predicate);
  * the cases the proof leaves are 258 good classes, 81: 37, 96: 42, 108: 6, 144: 15, 162: 29, 192: 30,
    216: 59, 243: 40, and each has a bad gate;
  * the superset with every window centre j >= 5 at Q < 280 (P >= 126 included: P = 127, 131, 133 at Q = 243,
    13 good classes each) is 297 good classes, and each has a bad gate;
  * the exception occurs at the eight Q only at 41/81 (81 = 3^4; 243 = 3^5 is an odd power), and that class
    has no bad gate -- so the exception is needed;
  * (beyond the finite check) every other good class at the eight Q has a bad gate;
  * CONTROLS (must fire): the bad-gate detector says no at the exception (above); with eps flipped in the
    window, the goodness rule disagrees with the chain signs. (Flipping eps everywhere does NOT break the
    conclusion at these eight Q -- the flipped window still leaves a bad gate in every class -- so that
    variant is no control here.)

Seconds.
"""
import math
from fractions import Fraction as F

from common import check, done

QS = [81, 96, 108, 144, 162, 192, 216, 243]


def three_smooth(lo, hi):
    out = []
    for Q in range(lo + 1, hi):
        m, b = Q, 0
        while m % 2 == 0:
            m //= 2
        while m % 3 == 0:
            m //= 3
            b += 1
        if m == 1 and b >= 1:
            out.append(Q)
    return out


def frac(x):
    return x - math.floor(x)


def gates(Q, P, varpi):
    """(n, z_n, nD) for 0 < n < P/2: z_n the integer = varpi (mod 2) within 1 of nD."""
    out = []
    for n in range(1, (P + 1) // 2):
        nD = F(2 * Q * n, P)
        assert nD.denominator != 1
        z = [c for c in range(math.floor(nD) - 1, math.floor(nD) + 3) if c % 2 == varpi and abs(c - nD) < 1]
        assert len(z) == 1
        out.append((n, z[0], nD))
    return out


def centre(Q, P, flip=False):
    """Per class of P/Q: (varpi, e, good, has_bad_gate, left_by_proof, exception); and the goodness cross-check."""
    N, D = Q // 3, F(2 * Q, P)
    eps = 1 if P % 3 == 2 else -1
    w = -eps if flip else eps

    def c(d):
        r = (d * P) % (4 * Q)
        assert r not in (Q, 3 * Q)
        return 1 if (r < Q or r > 3 * Q) else -1

    good, agree = {}, True
    for e in range(1, N):
        fe = frac(F(e * P, 2 * Q))
        rule = (F(1, 2) < fe < F(5, 6)) if eps == 1 else (F(1, 6) < fe < F(1, 2))
        chain = len({c(2 * N - e), c(e), -c(2 * N + e)}) == 1
        agree &= rule == chain
        good[e] = ((F(1, 2) < fe < F(5, 6)) if w == 1 else (F(1, 6) < fe < F(1, 2)))

    def frame_good(z):
        if z % N == 0:
            return False
        return good[z] if z < N else good[2 * N - z] if z < 2 * N else good[z - 2 * N]

    rows = []
    middle = F(Q, 3) < P < F(2 * Q, 3)
    inwin = F(7 * Q, 15) <= P <= F(5 * Q, 9)
    j = abs(Q - 2 * P)
    for varpi in (0, 1):
        G = gates(Q, P, varpi)
        left = (middle and varpi != N % 2
                and ((not inwin and Q < 102) or (inwin and j >= 5 and P < 126)))
        superset = middle and varpi != N % 2 and ((not inwin and Q < 102) or (inwin and j >= 5))
        for e in range(varpi or 2, N, 2):
            if not good[e]:
                continue
            bad = any(not frame_good(z) for n, z, nD in G if e < nD < 2 * N + e)
            exc = (Q % 2 == 1 and round(math.log(Q, 3)) % 2 == 0 and 3 ** round(math.log(Q, 3)) == Q
                   and P == (Q + 1) // 2 and varpi == 0 and e == N - 1)
            rows.append((varpi, e, bad, left, superset, exc))
    return rows, agree


def main():
    print('III Lemma A.15 -- the finite check its proof leaves\n')
    check(f'the eight Q: {QS}', three_smooth(80, 280) == QS)

    per_q, per_q_sup, all_bad, sup_bad, agree_all, exc_seen, others_ok = {}, {}, True, True, True, [], True
    for Q in QS:
        for P in range(3, Q - 1, 2):
            if math.gcd(P, Q) != 1:
                continue
            rows, agree = centre(Q, P)
            agree_all &= agree
            for varpi, e, bad, left, sup, exc in rows:
                if exc:
                    exc_seen.append((Q, P, varpi, e, bad))
                    continue
                others_ok &= bad
                if left:
                    per_q[Q] = per_q.get(Q, 0) + 1
                    all_bad &= bad
                if sup:
                    per_q_sup[Q] = per_q_sup.get(Q, 0) + 1
                    sup_bad &= bad
    check('goodness: the rule frac(eP/2Q) in the window agrees with constant chain signs at every class',
          agree_all)
    want = {81: 37, 96: 42, 108: 6, 144: 15, 162: 29, 192: 30, 216: 59, 243: 40}
    check(f'the proof leaves {sum(per_q.values())} good classes, per Q {per_q} (expect 258, {want})',
          per_q == want)
    check('each of them has a bad gate', all_bad)
    check(f'with every window centre j >= 5 included: {sum(per_q_sup.values())} good classes (expect 297), '
          f'each with a bad gate', sum(per_q_sup.values()) == 297 and sup_bad)
    check(f'the exception (Q = 3^(2k), P = (Q + 1)/2, varpi = 0, e = N - 1) occurs at {[x[:2] for x in exc_seen]} '
          f'and has no bad gate there', [x[:2] for x in exc_seen] == [(81, 41)] and not any(x[4] for x in exc_seen))
    check('(beyond the finite check) every other good class at the eight Q has a bad gate', others_ok)

    # CONTROLS. The bad-gate detector can say no: it does at the exception (asserted above). And the goodness
    # cross-check can fail: with eps flipped in the window it disagrees with the chain signs.
    flip_agree = True
    for Q in QS:
        for P in range(3, Q - 1, 2):
            if math.gcd(P, Q) != 1:
                continue
            N, eps = Q // 3, (1 if P % 3 == 2 else -1)
            for e in range(1, N):
                fe = frac(F(e * P, 2 * Q))
                flipped = (F(1, 2) < fe < F(5, 6)) if eps == -1 else (F(1, 6) < fe < F(1, 2))
                r = [(d * P) % (4 * Q) for d in (2 * N - e, e, 2 * N + e)]
                cs = [1 if (x < Q or x > 3 * Q) else -1 for x in r]
                flip_agree &= flipped == (len({cs[0], cs[1], -cs[2]}) == 1)
    check('CONTROL: with eps flipped, the goodness rule disagrees with the chain signs (fires)', not flip_agree)
    done()


if __name__ == '__main__':
    main()
