"""Appendix A.4, Proposition A.12 (the certificate; computer-assisted).

Proposition A.12: let p in {5, 7}; take one of the five sub-cases -- u0 odd with K = 4 (`K4`) or K >= 6 (`K6`),
u0 even with kappa = 1 (`k1`), kappa = 3 (`k3`) or kappa >= 5 (`k5`); a sign eps; b in {1/2 - eps/p, 3/2 - eps/p}
(mod 2), with zeta as in Lemma A.10 (3); an r in (0, 1/2) outside the finite set
B = U_{1<=e<=6} (1/ep)Z  U  U_{1<=m<=5} (1/m)Z; and a sign lambda admissible at r. Then for some e in {1, 3, 5}
(u0 odd), resp. {2, 4, 6} (u0 even), and some sign, the chain (e, +-) at theta = e(1 - r)/2, with the statuses of
Lemma A.11 and the cell form of Lemma A.10, has a member y that is a carrier with g = +1 and, across omega_1, a member
that is a non-carrier or a carrier with g = +1.

Its proof: the data change only at points of B, so the assertion is constant on each open interval between
consecutive points of B in [0, 1/2], and it is checked in exact rational arithmetic at each midpoint: 1 672 cases
(p, sub-case, interval, eps, b, lambda), every one certified. This script is that check.

`status` is Lemma A.11's table, `cell_offset` and `cell_g` Lemma A.10's cell form (the position of omega_1 among the
members, and g on member i), `lams` the admissibility rule of Proposition A.12 -- each transcribed from the paper.
Signs of sin and cos at rational multiples of pi are exact (`ssin`, `scos`: t in units of pi).

Statements checked:
  * CERT: 1 672 cases, 0 uncertified.
  * CONTROLS (must fire): C1 withholding Lemma A.11's forced non-carrier statuses leaves cases uncertified;
    C2 at p = 5, u0 even, kappa >= 5, with e in {2, 4} only, cases are uncertified (e = 6 is needed);
    C3 on r in (1/2, 1) (outside the proposition) cases are uncertified.
  * (beyond the certificate) END TO END at every dense lead of every even Q = 2^a p^b <= QMAX (p in {5, 7},
    N = Q/p >= 8 even, P odd = +-2 (mod p), P < Q/2): the lead's own lambda is admissible at r = P/Q; on the chains
    (e, +-) with e in {1, 3, 5} resp. {2, 4, 6}, the exact rows reproduce Lemma A.10's cell form (offset, origins,
    g) and every forced status of Lemma A.11 is the true carrier status; and the chain the certificate names is
    dead -- its sign conditions (m = chi w - g >= 0 on carriers, chi w >= 0 off them, rho = +-1) have no integer
    solution w. (Each condition bounds w by -1, 0 or 1, so searching |w| <= 3 is exhaustive.)

  python sA_12_certificate.py            QMAX = 1500 (about two minutes)
  python sA_12_certificate.py --qmax N
"""
import math
import sys
from fractions import Fraction as F
from itertools import product

from common import check, done

def oddprime_part(Q):
    m = Q
    while m % 2 == 0: m //= 2
    if m == 1: return None
    for p in range(3, m + 1, 2):
        if m % p == 0:
            mm = m
            while mm % p == 0: mm //= p
            return p if mm == 1 else None

def mk(P, Q):
    delta = Q - P
    def c(d):
        r = (d * P) % (4 * Q); return 1 if (r < Q or r > 3 * Q) else -1
    def sinb(z):
        r = (z * delta) % (4 * Q); assert r % (2 * Q) != 0; return 1 if r < 2 * Q else -1
    def cosb(z):
        r = (z * delta + Q) % (4 * Q); assert r % (2 * Q) != 0; return 1 if r < 2 * Q else -1
    return c, sinb, cosb

def chain_rows(P, Q, N, p, kap, r):
    c, _, _ = mk(P, Q); cs = c(kap - 1); odd_u0 = (kap % 2 == 0)
    def tau(d):
        return (-1) ** (((kap - 1 - d) // 2) % 2) if kap % 2 == 0 else (-1) ** (((d - kap - 1) // 2) % 2)
    d0 = r % (2 * N)
    while d0 - 2 * N > -Q: d0 -= 2 * N
    out = []
    for j in range(p):
        d = d0 + 2 * N * j; s = 1 if d > 0 else -1; a = abs(d)
        out.append(dict(d=d, a=a, x=d if d > 0 else d + Q, chi=s * (-1) ** j * c(a), rho=s * (-1) ** j * cs * tau(a),
                        g=cs * c(a) * tau(a), car=(a < kap) if odd_u0 else (a > kap)))
    return sorted(out, key=lambda m: m['x'])

def feasible(rows):
    for w in range(-3, 4):
        if all((m['chi'] * (w - m['rho']) >= 0) if m['car'] else (m['chi'] * w >= 0) for m in rows): return True
    return False

def omega1(P, Q, eps, rows):
    a = rows[0]['d']; return F(2 * a + eps * (Q - P * a), 2) % Q

def ssin(t):
    t = t % 2; return 0 if t in (0, 1) else (1 if t < 1 else -1)
def scos(t): return ssin(t + F(1, 2))

ES = {1: (1, 3, 5), 0: (2, 4, 6)}           # key: u0 odd -> 1, u0 even -> 0
SUBS = {1: ('K4', 'K6'), 0: ('k1', 'k3', 'k5')}

def status(p, odd_u0, sub, neg_class, i, e, hide_fn=False):
    """Lemma A.11's table.  Member i of the class (e, POS/NEG); it is of negative origin iff i even XOR POS... :
    POS: even i positive; NEG: even i negative."""
    neg = (i % 2 == 0) == neg_class
    FN = 'U' if hide_fn else 'FN'
    if odd_u0:
        if neg and i == 0:
            if sub == 'K4': return FN if e <= 3 else 'FC'
            return FN
        if p == 5:
            if not neg and i == 4: return 'FC' if (sub == 'K4' and e <= 3) else 'U'
            return 'FC'
        if not neg: return 'FC' if (i <= 3 or (sub == 'K4' and i == 4)) else 'U'
        return 'FC' if (i >= 3 or sub == 'K4') else 'U'
    if sub == 'k1': return 'FC'
    if sub == 'k3':
        if not neg and i == 0 and e == 2: return FN
        if neg and i == p - 1 and e >= 6: return 'U'
        return 'FC'
    if not neg and i == 0: return FN if e <= 4 else 'U'
    if p == 5: return 'U' if (neg and i == 4) else 'FC'
    if not neg: return 'FC' if i >= 3 else 'U'
    return 'FC' if i <= 3 else 'U'

def cell_offset(odd_u0, eps, neg_class, theta):
    return (eps * (theta + (0 if odd_u0 else F(1, 2))) + (F(1, 2) if neg_class else 0)) % 1

def cell_g(odd_u0, lam, zeta, neg_class, i, theta, b):
    t = theta + i * b; posorig = not ((i % 2 == 0) == neg_class)
    if odd_u0: return lam * (ssin(t) if posorig else zeta * scos(t))
    return lam * (scos(t) if posorig else zeta * ssin(t))

def cert_class(p, odd_u0, sub, eps, lam, zeta, b, neg_class, theta, e, hide_fn=False):
    om = cell_offset(odd_u0, eps, neg_class, theta) * p        # member i below omega_1 iff i < om
    mem = [(i < om, status(p, odd_u0, sub, neg_class, i, e, hide_fn), cell_g(odd_u0, lam, zeta, neg_class, i, theta, b))
           for i in range(p)]
    for iy, y in enumerate(mem):
        if y[1] == 'FC' and y[2] == 1:
            for iw, w in enumerate(mem):
                if iw != iy and w[0] != y[0] and (w[1] == 'FN' or (w[1] == 'FC' and w[2] == 1)): return True
    return False

def zeta_of(p, b):
    # delta = 2 p b (b taken mod 2 determines delta mod 4p); zeta = (-1)^((delta-1)/2)
    dl = (2 * p * b) % (4 * p); assert dl.denominator == 1; dl = int(dl)
    return (-1) ** (((dl - 1) // 2) % 2)

def lams(odd_u0, sub, r, zeta):
    if odd_u0 and sub == 'K4':
        if r <= F(1, 4) or ssin(3 * r / 2) != ssin(5 * r / 2): return []
        return [zeta * scos(5 * (1 - r) / 2)]
    if not odd_u0 and sub == 'k1': return [scos(1 - r)]
    if not odd_u0 and sub == 'k3': return [scos(2 * (1 - r))] if r < F(1, 4) else []
    return [1, -1]

def Bset(p, lo, hi):
    s = {lo, hi}
    for m in [e * p for e in range(1, 7)] + [1, 2, 3, 4, 5]:
        for j in range(0, 2 * m + 1):
            x = F(j, m)
            if lo < x < hi: s.add(x)
    return sorted(s)

def certificate(p, odd_u0, sub, es, hide_fn=False, lo=F(0), hi=F(1, 2)):
    n, fails = 0, []
    bps = Bset(p, lo, hi)
    for eps, half in product((1, -1), (F(1, 2), F(3, 2))):
        b = (half - F(eps, p)) % 2; zeta = zeta_of(p, b)
        for a_, b_ in zip(bps, bps[1:]):
            rm = (a_ + b_) / 2
            for lam in lams(odd_u0, sub, rm, zeta):
                n += 1
                if not any(cert_class(p, odd_u0, sub, eps, lam, zeta, b, ng, (e * (1 - rm) / 2) % 2, e, hide_fn)
                           for e in es for ng in (False, True)):
                    fails.append((eps, half, lam, a_, b_))
    return n, fails

def choose(p, odd_u0, sub, eps, lam, zeta, b, r):
    for e in ES[1 if odd_u0 else 0]:
        for ng in (False, True):
            if cert_class(p, odd_u0, sub, eps, lam, zeta, b, ng, (e * (1 - r) / 2) % 2, e): return e, ng
    return None



def end_to_end(qmax):
    """The dense leads at every even Q = 2^a p^b <= qmax: admissible lambda, cell form, forced statuses, dead chain."""
    st = dict(D=0, L=0, L_fail=0, F=0, F_fail=0, adm_fail=0, nochoice=0, alive=0)
    for Q in range(16, qmax + 1, 2):
        p = oddprime_part(Q)
        if p not in (5, 7):
            continue
        N = Q // p
        if N % 2 or N < 8:
            continue
        for P in range(3, Q // 2 + 1, 2):
            if math.gcd(P, Q) != 1 or P % p not in (2, p - 2):
                continue
            eps = 1 if P % p == 2 else -1
            c, sinb, cosb = mk(P, Q)
            delta = Q - P
            b = F(delta, 2 * p) % 2
            zeta = (-1) ** (((delta - 1) // 2) % 2)
            for kap in range(1, Q - 1):
                if c(kap - 1) != c(kap + 1) or kap * P > Q * P - Q:     # not a lead, or Lemma 5.12's range
                    continue
                odd_u0 = (kap % 2 == 0)
                if (kap <= 4 * N) if odd_u0 else (kap >= Q - 4 * N):     # sparse: Proposition A.9's other half
                    continue
                st['D'] += 1
                if odd_u0:
                    sub, lam = ('K4' if Q - kap == 4 else 'K6'), sinb(kap - 1)
                else:
                    sub, lam = ('k1' if kap == 1 else 'k3' if kap == 3 else 'k5'), cosb(kap + 1)
                r = F(P, Q)
                if lam not in lams(odd_u0, sub, r, zeta):
                    st['adm_fail'] += 1
                for e in ES[1 if odd_u0 else 0]:
                    theta = F(e * delta, 2 * Q)
                    for ng in (False, True):
                        rows = chain_rows(P, Q, N, p, kap, e + N if ng else e)
                        if any(m['a'] % N == 0 for m in rows):
                            continue
                        st['L'] += 1
                        om = omega1(P, Q, eps, rows)
                        ok = rows[0]['x'] == e and (om - e) % Q == cell_offset(odd_u0, eps, ng, theta) * Q
                        for i, m in enumerate(rows):
                            ok &= (m['d'] > 0) == (not ((i % 2 == 0) == ng))
                            ok &= m['g'] == cell_g(odd_u0, lam, zeta, ng, i, theta, b)
                            s_ = status(p, odd_u0, sub, ng, i, e)
                            if s_ != 'U':
                                st['F'] += 1
                                st['F_fail'] += (s_ == 'FC') != m['car']
                        st['L_fail'] += not ok
                ch = choose(p, odd_u0, sub, eps, lam, zeta, b, r)
                if ch is None:
                    st['nochoice'] += 1
                    continue
                e, ng = ch
                st['alive'] += feasible(chain_rows(P, Q, N, p, kap, e + N if ng else e))
    return st


def main():
    qmax = int(sys.argv[sys.argv.index('--qmax') + 1]) if '--qmax' in sys.argv else 1500
    print('III Proposition A.12 -- the certificate\n')
    tot, unc, c1, c3 = 0, 0, 0, 0
    for p in (5, 7):
        for ou in (1, 0):
            for sub in SUBS[ou]:
                n, f = certificate(p, bool(ou), sub, ES[ou])
                tot, unc = tot + n, unc + len(f)
                c1 += len(certificate(p, bool(ou), sub, ES[ou], hide_fn=True)[1])
                if sub in ('K6', 'k5'):
                    c3 += len(certificate(p, bool(ou), sub, ES[ou], lo=F(1, 2), hi=F(1))[1])
    c2 = len(certificate(5, False, 'k5', (2, 4))[1])
    check(f'CERT: {tot} cases (expect 1672), {unc} uncertified', tot == 1672 and unc == 0)
    check(f'CONTROL C1: without the forced non-carrier statuses, {c1} cases uncertified (fires)', c1 > 0)
    check(f'CONTROL C2: p = 5, u0 even, kappa >= 5, e in {{2, 4}} only: {c2} cases uncertified (fires)', c2 > 0)
    check(f'CONTROL C3: on r in (1/2, 1), {c3} cases uncertified (fires)', c3 > 0)
    st = end_to_end(qmax)
    check(f'end to end, even Q <= {qmax}: {st["D"]} dense leads, every lambda admissible ({st["adm_fail"]} not), '
          f'a certified chain at each ({st["nochoice"]} without)', st['D'] > 0 and not st['adm_fail'] and not st['nochoice'])
    check(f'  cell form on {st["L"]} chains ({st["L_fail"]} mismatches); {st["F"]} forced statuses '
          f'({st["F_fail"]} wrong)', st['L'] > 0 and not st['L_fail'] and st['F'] > 0 and not st['F_fail'])
    check(f'  the named chain is dead at every lead ({st["alive"]} alive)', not st['alive'])
    done()


if __name__ == '__main__':
    main()
