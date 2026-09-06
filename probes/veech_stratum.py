#!/usr/bin/env python3
"""Step 1 of the Veech reframing: identify + VERIFY the translation surface S_alpha
for the right triangle with angles (at O,R,A) = (P, Q, Q-P)*pi/(2Q), P odd, gcd(P,Q)=1.

Two INDEPENDENT computations of the genus/stratum, cross-checked:
  (A) cone-point / Gauss-Bonnet:  unfold with N=2Q -> 2N=4Q copies; at a vertex of
      reduced angle pi*m/n there are (2N)/(2n) = N/n cone points each of cone angle
      2*pi*m (order m-1).  Sum of orders must be 2g-2.
  (B) Zemlyakov-Katok genus formula:  g = 1 + (N/2) * sum_i (m_i - 1)/n_i.
Also cross-checks even-P: n(P/Q)=Q-1 should equal 2g (a checkable prediction).
"""
import sys, math
from fractions import Fraction as F


def reduce_angle(k, twoQ):
    """angle = k*pi/twoQ; return (m,n) reduced with angle = pi*m/n."""
    g = math.gcd(k, twoQ)
    return k // g, twoQ // g


def surface(P, Q):
    twoQ = 2 * Q
    verts = {'O': P, 'R': Q, 'A': Q - P}          # numerators over 2Q
    red = {v: reduce_angle(k, twoQ) for v, k in verts.items()}
    N = 1
    for (m, n) in red.values():
        N = N * n // math.gcd(N, n)                # lcm of denominators
    copies = 2 * N
    # (A) cone points
    sings = []  # (vertex, #points, cone_angle_over_pi, order)
    for v, (m, n) in red.items():
        npts = N // n
        order = m - 1
        sings.append((v, npts, 2 * m, order))      # cone angle = 2*pi*m
    sum_orders = sum(npts * order for (_, npts, _, order) in sings)
    gA = F(sum_orders, 2) + 1                       # 2g-2 = sum orders
    # (B) ZK formula
    gB = 1 + F(N, 2) * sum(F(m - 1, n) for (m, n) in red.values())
    # nonzero-order singularities => stratum
    strat = []
    for (v, npts, _, order) in sings:
        if order > 0:
            strat += [order] * npts
    strat.sort(reverse=True)
    return dict(N=N, copies=copies, red=red, sings=sings, gA=gA, gB=gB, stratum=strat)


def main():
    print(f"{'P':>3} {'Q':>3} {'Qpar':>4} {'N':>4} {'copies':>6} {'gA':>4} {'gB':>4} {'floorQ/2':>8} {'stratum':>20} {'2g=?n(evenP)':>12}")
    ok = True
    grid = []
    for Q in range(3, 24):
        for P in range(1, Q):
            if math.gcd(P, Q) != 1 or P % 2 == 0:   # P odd, coprime
                continue
            grid.append((P, Q))
    for P, Q in grid:
        s = surface(P, Q)
        gA, gB = s['gA'], s['gB']
        assert gA.denominator == 1 and gB.denominator == 1, (P, Q, gA, gB)
        gA, gB = int(gA), int(gB)
        fq = Q // 2
        agree = (gA == gB == fq)
        ok &= agree
        even_note = ""
        # even-P prediction n=Q-1=2g ; here P odd so instead note 2g vs Q
        flag = "" if agree else "  <-- MISMATCH"
        if Q <= 13 or not agree:
            print(f"{P:>3} {Q:>3} {'odd' if Q%2 else 'even':>4} {s['N']:>4} {s['copies']:>6} {gA:>4} {gB:>4} {fq:>8} {str(s['stratum']):>20} {2*gA:>12}{flag}")
    print(f"\nAll {len(grid)} cases (P odd, coprime, Q<24): gA==gB==floor(Q/2) ? {ok}")
    # even-P sanity: build a few even-P surfaces and check 2g == Q-1 (= n for even P, Thm B')
    print("\nEven-P check (Thm B': n=Q-1 should = 2g):")
    for (P, Q) in [(2, 5), (2, 7), (4, 9), (2, 9), (4, 11), (6, 13), (2, 11)]:
        if math.gcd(P, Q) != 1:
            continue
        s = surface(P, Q)
        g = int(s['gA'])
        print(f"  P={P} Q={Q}: g={g}, 2g={2*g}, Q-1={Q-1}  {'OK' if 2*g == Q-1 else 'DIFF'}  stratum {s['stratum']}")


if __name__ == '__main__':
    main()
