"""[Deposit] Proposition 4.3 at m = 3 (Appendix A, "The inequalities"): a symbolic check of the hand proof. Paper: "Finite
perpendicular-beam partitions at irrational angles". New for the deposit (the research check is the session's sympy run recorded in
the research ledger; it replaced an interval-arithmetic script). It checks the identities the proof uses -- at m = 3,
  B = 2 cot(a) sin(6a) = sigma^3 - 4 sigma^2 + 3 sigma = -2 sigma^2 + 7 sigma - 1  mod  sigma^3 - 2 sigma^2 - 4 sigma + 1,
  4 cos(eps + 2a) = -4 cos(4a) = 4 - 2 (sigma - 2)^2,
so that (b) is 6 sigma^2 - 20 sigma + 4 > 0 and (c) is -8 sigma^2 + 28 sigma - 8 > 0 -- and the exact rational values on the bracket
3.16 < sigma < 3.17, and compares with the numerical root.
"""
import sympy as sp

import _setup  # noqa: F401  (deposit: this directory on the path, cwd = irrational_beams/)


def main():
    c, S = sp.symbols('c sigma')
    ok = []
    B_c = sp.expand(2 * c * sp.chebyshevu(5, c))                          # 2 cos a * sin(6a)/sin(a)
    Bs = sp.expand(B_c.subs(c, sp.sqrt(S) / 2))
    cub = S ** 3 - 2 * S ** 2 - 4 * S + 1
    ok.append(('B = sigma^3 - 4 sigma^2 + 3 sigma', sp.expand(Bs - (S ** 3 - 4 * S ** 2 + 3 * S)) == 0))
    ok.append(('B = -2 sigma^2 + 7 sigma - 1 modulo the cubic', sp.rem(Bs, cub, S) == -2 * S ** 2 + 7 * S - 1))
    C4 = sp.expand((-4 * sp.chebyshevt(4, c)).subs(c, sp.sqrt(S) / 2))  # 4 cos(eps + 2a) = -4 cos 4a at m = 3
    ok.append(('4 cos(eps + 2a) = 4 - 2 (sigma - 2)^2', sp.expand(C4 - (4 - 2 * (S - 2) ** 2)) == 0))
    b = sp.rem(sp.expand(S + 1 - 3 * Bs), cub, S)
    cc = sp.rem(sp.expand(C4 - (S + 1 - 3 * Bs)), cub, S)
    ok.append(('(b) is 6 sigma^2 - 20 sigma + 4', sp.expand(b - (6 * S ** 2 - 20 * S + 4)) == 0))
    ok.append(('(c) is -8 sigma^2 + 28 sigma - 8', sp.expand(cc - (-8 * S ** 2 + 28 * S - 8)) == 0))
    lo, hi = sp.Rational(316, 100), sp.Rational(317, 100)
    ok.append(('the cubic changes sign on [3.16, 3.17] (and increases on sigma >= 3)',
               cub.subs(S, lo) < 0 < cub.subs(S, hi) and all(sp.diff(cub, S).subs(S, t) > 0 for t in (3, lo, hi))))
    ok.append(('(b) increases for sigma > 5/3 and is 446/625 > 0 at 3.16', sp.solve(sp.diff(b, S), S) == [sp.Rational(5, 3)]
               and b.subs(S, lo) == sp.Rational(446, 625)))
    ok.append(('(c) decreases for sigma > 7/4 and is 461/1250 > 0 at 3.17', sp.solve(sp.diff(cc, S), S) == [sp.Rational(7, 4)]
               and cc.subs(S, hi) == sp.Rational(461, 1250)))
    r = [x for x in sp.Poly(cub, S).nroots(n=30) if 3 < x < 3.25][0]
    al = sp.acos(sp.sqrt(r) / 2)
    e = sp.pi - 6 * al
    num = (abs(sp.N(2 * sp.cot(al) * sp.sin(e) - Bs.subs(S, r), 25)) < 1e-20
           and abs(sp.N(4 * sp.cos(e + 2 * al) - C4.subs(S, r), 25)) < 1e-20)
    ok.append((f'numerical root sigma = {sp.N(r, 8)}: B and 4 cos(eps + 2a) agree with the polynomials; margins '
               f'{sp.N(b.subs(S, r), 3)}, {sp.N(cc.subs(S, r), 3)}', num))
    for label, v in ok:
        print(f'  [{"ok" if v else "FAIL"}] {label}')
    print('ALL OK' if all(v for _, v in ok) else 'FAIL')


if __name__ == '__main__':
    main()
