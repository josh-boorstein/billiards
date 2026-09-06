#!/usr/bin/env python3
"""Veech reframing, step-1 follow-on: IMPORT Apisa 2110.07540 and reconcile with our
§101 surface identification (veech_stratum.py).

Apisa's picture for the right triangle Q_{a,n}, angles (a/2n, b/2n, 1/2)*pi, a+b=n,
gcd(a,2n)=1  [our dictionary:  n=Q,  {a,b}={P,Q-P},  the right angle = pi/2]:

  * The pillowcase double is a (2n)-differential on the sphere with cone angles
    (a/n, (n-a)/n, 1)*pi at 0, inf, 1.
  * Pull back under z -> z^n  =>  a QUADRATIC differential on the sphere
        (S,q) = Q(a-2, b-2, -1^n) = Q(P-2, Q-P-2, -1^Q),
    cone angles a*pi, (n-a)*pi, and n simple poles of angle pi.
  * Our unfolding S_alpha = (X,omega) is the HOLONOMY (orientation) DOUBLE COVER of (S,q).

This script (A) reconstructs our §101 stratum from Apisa's base by the standard
holonomy-double-cover lift rule and checks it AGREES with veech_stratum.surface(),
and (B) prints Apisa's Veech/non-Veech dichotomy (Thm 7.6 + Lemma 7.4) and the
all-directions weak-asymptotic constants (Thm 1.1, ccyl of Cor 8.3).

Lift rule (half-translation zero of order d -> abelian holonomy double cover):
  d odd  -> ONE zero of order d+1        (cover branched there)
  d even -> TWO zeros of order d/2       (unbranched, swapped by the involution)
  d = -1 (simple pole), unmarked -> ONE regular point (order 0)
"""
import math
from fractions import Fraction as F
import veech_stratum as vs


def holonomy_double_lift(orders):
    """orders: list of quadratic-differential zero/pole orders (poles = -1).
    Return the abelian stratum (list of zero orders, poles unmarked -> dropped)."""
    out = []
    for d in orders:
        if d == -1:
            continue                      # simple pole -> regular point (unmarked)
        if d % 2 == 1:                    # odd -> single zero, order d+1
            out.append(d + 1)
        else:                             # even -> two zeros, order d/2 each
            out += [d // 2, d // 2]
    # drop order-0 (regular MARKED points) to match veech_stratum's genuine-zero list
    return sorted([o for o in out if o > 0], reverse=True)


def apisa_base_stratum(P, Q):
    """(S,q) = Q(P-2, Q-P-2, -1^Q)."""
    return [P - 2, Q - P - 2] + [-1] * Q


def check_all():
    ok = True
    n_checked = 0
    for Q in range(5, 24):
        for P in range(1, Q):
            if math.gcd(P, Q) != 1 or P % 2 == 0:
                continue
            base = apisa_base_stratum(P, Q)
            # Gauss-Bonnet on the base: sum of orders = -4 (genus-0 quad diff)
            assert sum(base) == -4, (P, Q, sum(base))
            lifted = holonomy_double_lift(base)
            ours = vs.surface(P, Q)['stratum']
            if lifted != ours:
                print(f"  MISMATCH P={P} Q={Q}: Apisa-lift {lifted}  vs  §101 {ours}")
                ok = False
            n_checked += 1
    print(f"Holonomy-double-cover of Apisa base  ==  §101 stratum "
          f"for all {n_checked} cases (P odd, coprime, 5<=Q<24)?  {ok}")
    return ok


def dichotomy_table():
    print("\nApisa Veech/non-Veech dichotomy (Thm 7.6 + Lemma 7.4):")
    print("  rank one (Veech / dbl-cover-of-Veech)  <=>  min(P,Q-P) in {1,2}")
    print(f"  {'P':>3} {'Q':>3} {'min(P,Q-P)':>11} {'class':>26} {'base stratum'}")
    for Q in range(5, 14):
        for P in range(1, Q):
            if math.gcd(P, Q) != 1 or P % 2 == 0:
                continue
            mn = min(P, Q - P)
            if mn <= 2:
                cls = "VEECH (rank 1)"
            else:
                cls = "full quad-double (rank>=2)"
            base = f"Q({P-2},{Q-P-2},-1^{Q})"
            print(f"  {P:>3} {Q:>3} {mn:>11} {cls:>26} {base}")


def constants():
    print("\nAll-directions weak asymptotics (companion to our single-direction n(P/Q)):")
    print("  Thm 1.1:  N_P(L) ~ (1/16pi)(1-1/Q)(1+2/(P(Q-P))) L^2/area,  min(P,Q-P)>2")
    print("  Cor 8.3:  ccyl(base) = Q(Q-1)/(4pi^2) * (1 + 2/(P(Q-P)))   [k1=P-2,k2=Q-P-2]")
    print(f"  {'P':>3} {'Q':>3} {'(1-1/Q)(1+2/P(Q-P))':>22}")
    for (P, Q) in [(3, 7), (5, 12), (3, 11), (5, 8), (7, 12), (3, 8), (5, 14)]:
        if math.gcd(P, Q) != 1:
            continue
        val = F(Q - 1, Q) * (1 + F(2, P * (Q - P)))
        print(f"  {P:>3} {Q:>3} {str(val):>22}  = {float(val):.5f}")


def orphan_weierstrass():
    """Orphan = Weierstrass point (reframing step 2), as a branch-parity statement.

    S_alpha = holonomy double cover of (S,q); the deck involution tau (tau*omega=-omega)
    IS the hyperelliptic involution, quotient = sphere, Weierstrass points = fixed points
    of tau = ramification points = preimages of ODD-order singularities of (S,q).
    The O cone point has base order P-2.  So:
        O is a Weierstrass point  <=>  P-2 odd  <=>  P odd  <=>  orphan exists (Thm D).
    Check the two 'P odd' characterisations coincide over BOTH parities of P."""
    print("\nOrphan = Weierstrass point (Thm D reframed):  "
          "O is a branch pt of the holonomy cover <=> P odd <=> orphan exists")
    ok = True
    for Q in range(4, 16):
        for P in range(1, Q):
            if math.gcd(P, Q) != 1:
                continue
            o_order = P - 2                      # base order at O (cone angle P*pi)
            O_is_weierstrass = (o_order % 2 == 1)  # odd order => holonomy-cover branch pt
            orphan_exists = (P % 2 == 1)           # Thm D
            if O_is_weierstrass != orphan_exists:
                print(f"  MISMATCH P={P} Q={Q}")
                ok = False
    print(f"  'O Weierstrass' == 'orphan exists' for all P (both parities), Q<16?  {ok}")


if __name__ == '__main__':
    check_all()
    orphan_weierstrass()
    dichotomy_table()
    constants()
