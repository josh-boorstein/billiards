"""Appendix A.3: the conditions C1-C12c (Theorem 10.8) and D1-D13 (Theorem 10.12), exactly as printed.

Statements checked:
  * each printed condition, written out below from the paper's list, is an odd polynomial in
    S = sin x, with no root in (0, sin(pi/20)] (C) resp. (0, sin(pi/28)] (D), and positive there --
    by exact real-root isolation (sympy) on the polynomial with the trivial factor S^k removed;
  * the remarks of Theorem 10.8's proof: C6 and C9a vanish to first order in x and are positive at
    third order, 16 x^3 and 32 x^3; at Q = 7 (x = pi/14) C3, C4 and C12b fail and C2 and C6 are
    equalities, C6 through its factor 8S^3 - 4S^2 - 4S + 1, which vanishes at S = sin(pi/14);
    C12a repeats C4 and C9b repeats C8;
  * the remark of Theorem 10.12's proof: the tightest condition is D9, the only one failing at
    Q = 11 (x = pi/22);
  * the control of A.3: the window prover (`window_prover.py`), which derives the conditions of a
    walk mechanically, generates every printed C-condition on the far arc at P = 3, r = 1, and
    every printed D-condition but D11 on the near arc at P = 5, r = 1 -- D11 (s_13 > s_3) is a step
    of the hand proof's chain p''' < a + s_3 < a + s_13, which the prover tests directly.
Exact zero tests at x = pi/14, pi/22 are remainders modulo the minimal polynomial of sin x; signs
are read at 50 digits.
"""
import sympy as sp

from common import check, done
import window_prover as W

X = W.X
S = sp.symbols('S')


def s(k):
    return sp.sin(k * X)


C = {'C1': s(5) - s(1) - s(3), 'C2': s(9) - s(5), 'C3': s(11) - (s(1) + s(3) + s(9) - s(5)),
     'C4': s(3) + s(9) - 2 * s(5), 'C6': 2 * s(9) + s(3) - 2 * s(5) - s(11),
     'C7': 2 * s(5) + s(11) - s(1) - 2 * s(3) - s(9),
     'C8': 3 * s(5) + s(11) - s(1) - 2 * s(3) - 2 * s(9),
     'C9a': 2 * s(1) + 2 * s(3) + 2 * s(9) - 3 * s(5) - s(11),
     'C9b': 3 * s(5) + s(11) - s(1) - 2 * s(3) - 2 * s(9),
     'C10': 3 * s(3) + s(1) + 2 * s(9) - 3 * s(5) - s(11),
     'C11': 2 * s(3) + 2 * s(9) - 2 * s(5) - s(11), 'C12a': s(3) + s(9) - 2 * s(5),
     'C12b': 4 * s(5) + 2 * s(11) - 2 * s(1) - 3 * s(3) - 3 * s(9),
     'C12c': 4 * s(3) + 4 * s(9) + s(1) - 5 * s(5) - 2 * s(11)}
a = s(7) - s(3)
p_star = a + 2 * s(9) - 2 * s(1) - s(3) - s(13)
D = {'D1': s(5) - a, 'D2': s(7) - s(5), 'D3': a + s(9) - s(5), 'D4': a + s(9) - s(1) - s(5),
     'D5': s(7) + s(9) - s(1) - 2 * s(3) - s(5), 'D6': s(1) + s(13) - s(9),
     'D7': s(7) + 2 * s(9) - s(1) - 2 * s(3) - s(5) - s(13), 'D8': p_star - a, 'D9': s(5) - p_star,
     'D10': 2 * s(1) + 2 * s(3) + s(13) - 2 * s(9), 'D11': s(13) - s(3),
     'D12': a + (s(9) - s(1)) / 2 - s(5), 'D13': s(1) + s(3) + s(13) - s(9) - (s(9) - s(1)) / 2}
# D8 is printed as p* > a, "that is 2 s_9 > 2 s_1 + s_3 + s_13": the two forms must agree
D8_PRINTED_FORM = 2 * s(9) - 2 * s(1) - s(3) - s(13)


def poly(e):
    return W.to_poly(sp.expand(e), S)


def strip(q):
    """(k, q / S^k) with q / S^k nonzero at 0."""
    k = 0
    while q.eval(0) == 0:
        q, k = sp.Poly(sp.cancel(q.as_expr() / S), S), k + 1
    return k, q


def positive_on(e, Q0):
    """e > 0 on (0, sin(pi/2Q0)]: no root of q = e / S^k in (0, b], b rational >= sin(pi/2Q0), and
    q(0) > 0."""
    k, q = strip(poly(e))
    b = sp.Rational(int(sp.N(sp.sin(sp.pi / (2 * Q0)), 30) * 10**12) + 1, 10**12)
    return q.eval(0) > 0 and q.count_roots(0, b) == 0, k, q


def at(e, Q):
    """(exactly zero?, sign at 50 digits) of e at x = pi/2Q."""
    mp = sp.minimal_polynomial(sp.sin(sp.pi / (2 * Q)), S)
    q = poly(e)
    zero = sp.rem(q.as_expr(), mp, S) == 0
    v = sp.N(q.as_expr().subs(S, sp.sin(sp.pi / (2 * Q))), 50)
    return zero, (0 if zero else (1 if v > 0 else -1))


def main():
    print('Appendix A.3 -- the printed conditions of Theorems 10.8 and 10.12\n')
    for name, lst, Q0 in (('C1-C12c', C, 10), ('D1-D13', D, 14)):
        bad = [k for k, e in lst.items() if not positive_on(e, Q0)[0]]
        check(f'{name}: {len(lst)} odd polynomials in sin x, root-free and positive on '
              f'(0, sin(pi/{2 * Q0})]' + (f' -- fail: {bad}' if bad else ''), not bad)
    check('D8 printed two ways (p* > a, 2 s_9 > 2 s_1 + s_3 + s_13): the same polynomial',
          poly(D['D8']) == poly(D8_PRINTED_FORM))
    check('C12a = C4 and C9b = C8, as printed', poly(C['C12a']) == poly(C['C4'])
          and poly(C['C9b']) == poly(C['C8']))
    lead = {}
    for k in ('C6', 'C9a'):
        _, kk, q = positive_on(C[k], 10)
        lead[k] = (kk, q.eval(0))
    check(f'C6, C9a vanish to first order and are 16 x^3, 32 x^3 at third order: {lead}',
          lead == {'C6': (3, 16), 'C9a': (3, 32)})
    sig7 = {k: at(e, 7) for k, e in C.items()}
    fails = sorted(k for k, (z, sg) in sig7.items() if sg < 0 and k not in ('C12a', 'C9b'))
    eqs = sorted(k for k, (z, sg) in sig7.items() if z and k not in ('C12a', 'C9b'))
    check(f'Q = 7: fail {fails}, equalities {eqs} (paper: C3, C4, C12b fail; C2, C6 equal)',
          fails == ['C12b', 'C3', 'C4'] and eqs == ['C2', 'C6'])
    f6 = 8 * S**3 - 4 * S**2 - 4 * S + 1
    check("C6's polynomial has the factor 8S^3 - 4S^2 - 4S + 1, zero at S = sin(pi/14)",
          sp.rem(poly(C['C6']).as_expr(), f6, S) == 0
          and sp.rem(f6, sp.minimal_polynomial(sp.sin(sp.pi / 14), S), S) == 0)
    sig11 = {k: at(e, 11)[1] for k, e in D.items()}
    neg = sorted(k for k, sg in sig11.items() if sg <= 0)
    check(f'Q = 11: the only failing condition is D9: {neg}', neg == ['D9'])

    def gen(P, r, side):
        seen, orig = [], W.q0_for
        W.q0_for = lambda e: (seen.append(sp.expand(e)), orig(e))[1]
        try:
            out = W.prove(P, r, side, verbose=False)
        finally:
            W.q0_for = orig
        norm = set()
        for e in seen:
            if e != 0:
                q = poly(e)
                norm.add(sp.Poly(q.as_expr() / abs(q.LC()), S).as_expr())
        return out, norm

    for name, lst, fam, expect_missing in (('C', C, (3, 1, 'far'), []),
                                           ('D', D, (5, 1, 'near'), ['D11'])):
        out, G = gen(*fam)
        miss = [k for k, e in lst.items()
                if sp.Poly(poly(e).as_expr() / abs(poly(e).LC()), S).as_expr() not in G]
        check(f'control: the prover on {fam} ({out[0]} cells, Q >= {out[1]}) generates every printed '
              f'{name}-condition except {miss or "none"}', miss == expect_missing)
    done()


if __name__ == '__main__':
    main()
