"""Sections 6 and 7: the unfolded surface S_alpha, its hyperelliptic involution and its base,
BUILT from the gluing and checked against Theorem 6.1, Corollary 6.2, Propositions 7.1 and 7.2,
Theorem 7.3 and Remark 7.4 at every coprime (P, Q) with Q <= 200, and against McMullen's
Weierstrass criterion on the slice min(P, Q - P) = 1 (Section 7.5, Appendix A.4).

How.  Nothing below uses the formulas it checks.  A copy of T in the Katok-Zemlyakov
unfolding is labelled by its linear part g, an element z -> zeta^k z or z -> zeta^k conj z of
the dihedral group (zeta = exp(i pi / 2Q)); the copy across side s of copy g is g rho_s, rho_s
the reflection in the base side s.  The copies are generated from the identity, and
  * a cone point over a vertex v is an orbit of copies under the two reflections in the sides
    at v; its cone angle is (orbit size) x (angle at v), so its order is read off the orbit;
  * the genus is 1 - chi/2 with chi = V - E + F, V the number of those points, E = 3F/2;
  * tau is the rotation by pi, acting on copies on the left; its fixed points on S_alpha are
    exactly the tau-invariant vertex orbits (it fixes no copy and no edge, since no reflection
    is -1);
  * the base S_alpha / tau has, at a fixed point of order d, cone angle pi (d + 1), i.e. order
    d - 1, and at a swapped pair of order d one point of order 2d.
Statements checked:
  * (Thm 6.1) 4Q copies; genus g = floor(Q/2); the Q points over R are regular; the cone
    orders over O and A are the theorem's row, and sum to 2g - 2;
  * (Cor 6.2) at P = 1, P = 2, P = Q - 2 and P = Q - 1 the stratum is H(2g - 2) or
    H(g - 1, g - 1), as the corollary lists them;
  * (Props 7.1, 7.2) tau fixes exactly the Q points over R, the point over O iff P odd and the
    point over A iff Q - P odd -- 2g + 2 points, so tau is the hyperelliptic involution -- and
    the base's singular orders are {P - 2, Q - P - 2, -1 (Q times)}, Gauss-Bonnet sum -4;
  * (Thm 7.3) the cone point over O is tau-fixed iff P is odd; when P is even it has two
    preimages and tau swaps them;
  * (Rem 7.4) 5/8, 3/5 and 3/7 lie in H(4, 2), H(2, 0, 0) and H(2, 1, 1);
  * (Section 7.5) on min(P, Q - P) = 1, Q = 4..15: Kenyon-Smillie's model is the regular
    n-gon, n = 2Q, of genus g, and McMullen's criterion -- its cone point is a Weierstrass point
    iff n != 4g + 2 -- agrees with tau-fixedness of the built surface's cone points over O at
    P = Q - 1 and over A at P = 1, in both parities of Q: one fixed point when n = 4g, two
    swapped ones when n = 4g + 2 (the n-gon's two vertex classes).  g is the built genus.
"""
import math
from collections import Counter

from common import check, done

QMAX = 200


class Surface:
    def __init__(self, P, Q):
        self.P, self.Q, self.M = P, Q, 4 * Q
        self.base = {'L2': 0, 'L1': Q, 'H': P}                 # side directions, units pi/2Q
        self.rho = {s: (2 * d % self.M, 1) for s, d in self.base.items()}
        self.copies = self._generate()

    def mul(self, a, b):
        (k1, f1), (k2, f2) = a, b
        return ((k1 + k2) % self.M, f2) if f1 == 0 else ((k1 - k2) % self.M, 1 - f2)

    def _generate(self):
        seen, todo = {(0, 0)}, [(0, 0)]
        while todo:
            g = todo.pop()
            for r in self.rho.values():
                h = self.mul(g, r)
                if h not in seen:
                    seen.add(h)
                    todo.append(h)
        return seen

    def vertex_orbits(self, s1, s2):
        left, orbits = set(self.copies), []
        while left:
            g = left.pop()
            orb, todo = {g}, [g]
            while todo:
                x = todo.pop()
                for s in (s1, s2):
                    y = self.mul(x, self.rho[s])
                    if y not in orb:
                        orb.add(y)
                        todo.append(y)
            left -= orb
            orbits.append(frozenset(orb))
        return orbits

    def analyse(self):
        """(copies, genus, {vertex: [(order, tau-fixed, tau-image index)]})."""
        P, Q, M = self.P, self.Q, self.M
        angle = {'O': P, 'A': Q - P, 'R': Q}                   # vertex angles, units pi/2Q
        sides = {'O': ('L2', 'H'), 'A': ('L1', 'H'), 'R': ('L1', 'L2')}
        tau = (2 * Q % M, 0)
        pts = {}
        for v in 'ORA':
            orbs = self.vertex_orbits(*sides[v])
            pts[v] = []
            for orb in orbs:
                num = len(orb) * angle[v]                      # cone angle, units pi/2Q
                assert num % M == 0
                d = num // M - 1                               # order
                img = frozenset(self.mul(tau, g) for g in orb)
                pts[v].append((d, img == orb, orbs.index(img)))
        F = len(self.copies)
        V = sum(len(x) for x in pts.values())
        chi = V - 3 * F // 2 + F
        return F, 1 - chi // 2, pts


def expected_row(P, Q):
    if P % 2 and Q % 2 == 0:
        return [P - 1], [Q - P - 1]
    if P % 2 and Q % 2:
        return [P - 1], [(Q - P) // 2 - 1] * 2
    return [P // 2 - 1] * 2, [Q - P - 1]


def stratum(pts):
    return sorted((d for v in 'OA' for d, _, _ in pts[v]), reverse=True)


def main():
    print(f'Sections 6-7 -- S_alpha built from the gluing, every coprime (P, Q), Q <= {QMAX}\n')
    rows = [(P, Q) for Q in range(2, QMAX + 1) for P in range(1, Q) if math.gcd(P, Q) == 1]
    c = Counter()
    ext = Counter()
    for P, Q in rows:
        F, g, pts = Surface(P, Q).analyse()
        c['copies'] += F == 4 * Q
        c['genus'] += g == Q // 2
        eo, ea = expected_row(P, Q)
        dO = sorted(d for d, _, _ in pts['O'])
        dA = sorted(d for d, _, _ in pts['A'])
        c['row'] += dO == sorted(eo) and dA == sorted(ea)
        c['sum'] += sum(dO) + sum(dA) == 2 * g - 2
        c['R_regular'] += len(pts['R']) == Q and all(d == 0 for d, _, _ in pts['R'])
        fix = {v: sum(f for _, f, _ in pts[v]) for v in 'ORA'}
        c['W_R'] += fix['R'] == Q
        c['W_O'] += fix['O'] == (P % 2)
        c['W_A'] += fix['A'] == ((Q - P) % 2)
        c['W_count'] += sum(fix.values()) == 2 * g + 2
        c['T73'] += (len(pts['O']) == 1 and pts['O'][0][1]) if P % 2 else (
            len(pts['O']) == 2 and pts['O'][0][2] == 1 and pts['O'][1][2] == 0)
        base = []
        for v in 'ORA':
            done_pairs = Counter()
            for d, f, _ in pts[v]:
                if f:
                    base.append(d - 1)
                else:
                    done_pairs[d] += 1
            for d, m in done_pairs.items():
                base += [2 * d] * (m // 2)
        sing = sorted(e for e in base if e != 0)
        want = sorted([e for e in (P - 2, Q - P - 2) if e != 0] + [-1] * Q)
        c['base'] += sing == want and sum(base) == -4
        # Cor 6.2: the extremes, as the corollary lists them (zero orders are regular points)
        st = [d for d in stratum(pts) if d]
        cases = []
        if P == 1 or P == Q - 1:
            cases.append([2 * g - 2] if Q % 2 == 0 else [g - 1, g - 1])
        if P == 2 or (P == Q - 2 and Q % 2):
            cases.append([2 * g - 2])
        for want_st in cases:
            ext['n'] += 1
            ext['ok'] += st == [d for d in want_st if d]
    n = len(rows)
    n3 = sum(1 for P, Q in rows if Q >= 3)
    check(f'{n} coprime (P, Q) with Q <= {QMAX} ({n3} with Q >= 3; the extra one is 1/2)', n > 0)
    for key, label in [('copies', '4Q copies of T'),
                       ('genus', 'genus g = floor(Q/2), from V - E + F'),
                       ('row', "cone orders over O and over A = Theorem 6.1's row"),
                       ('sum', 'the orders sum to 2g - 2'),
                       ('R_regular', 'the Q points over R are regular'),
                       ('W_R', 'tau fixes all Q points over R'),
                       ('W_O', 'tau fixes the point over O iff P is odd'),
                       ('W_A', 'tau fixes the point over A iff Q - P is odd'),
                       ('W_count', 'tau has exactly 2g + 2 fixed points (hyperelliptic; Prop 7.2)'),
                       ('T73', 'Thm 7.3: one tau-fixed point over O at odd P, two swapped at even P'),
                       ('base', 'the base has orders {P-2, Q-P-2, -1^Q}, sum -4 (Prop 7.1)')]:
        check(f'{label}: {c[key]}/{n}', c[key] == n)
    check(f'Cor 6.2: the stratum at P = 1, 2, Q-2, Q-1 is the one listed: {ext["ok"]}/{ext["n"]}',
          ext['n'] > 0 and ext['ok'] == ext['n'])
    rem = {(5, 8): [4, 2], (3, 5): [2, 0, 0], (3, 7): [2, 1, 1]}
    for (P, Q), want_st in rem.items():
        st = stratum(Surface(P, Q).analyse()[2])
        check(f'Rem 7.4: {P}/{Q} lies in H{tuple(st)}', st == want_st)
    # Section 7.5: McMullen's criterion on min(P, Q - P) = 1
    agree = total = 0
    for Q in range(4, 16):
        for P, v in ((Q - 1, 'O'), (1, 'A')):
            if math.gcd(P, Q) != 1:
                continue
            _F, g, pts = Surface(P, Q).analyse()
            cone = [x for x in pts[v] if x[0] > 0]
            n_gon = 2 * Q
            mcm = n_gon != 4 * g + 2
            total += 1
            # n = 4g: one vertex class; n = 4g + 2: two, which tau must swap (so neither is fixed)
            agree += len(cone) == (1 if mcm else 2) and all(f == mcm for _, f, _ in cone)
    check(f'Section 7.5: McMullen\'s criterion (n = 2Q, cone point Weierstrass iff n != 4g + 2) '
          f'= tau-fixedness of the built surface, P in {{1, Q-1}}, Q = 4..15: {agree}/{total}',
          total > 0 and agree == total)
    done()


if __name__ == '__main__':
    main()
