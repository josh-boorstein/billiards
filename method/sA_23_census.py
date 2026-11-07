"""Appendix A.6, Proposition A.23 (computer-assisted): the census at Q <= 80.

Proposition A.23: let Q <= 80 have at most one odd prime factor, and let P be coprime to Q with
min(P, Q - P) >= 2. For every launch 0 < u0 < P no head-on-free O-launched connection of P/Q ends at O, A or R.

The method is A.6's, transcribed:
  * THE DIRECTION AUTOMATON of a launch: states (s_j, t_j, lambda_j mod 2, j mod 2, v_j) -- the side just met, the
    frame, the parity of the number of reflections in L1, the parity of j, the direction in Z/4Q; initial state
    (L1, 0, 1, 1, 2Q - u0); an edge for each next side the ray can meet (the open arcs of A.6, `adm`) and does not
    meet head-on (`is_headon`), the frame changing by (-1)^j * 2 at a letter H (here `h` counts it in units of 2).
    A final state is one whose side is opposite the target (L1 for O, L2 for A, H for R).
  * LINEARITY: an edge into L1 from frame t, parity lambda weighs 2(-1)^lambda (x^(t+1) + x^(t-1)), the launch
    2(x + x^-1), a final state its lead; the map L(a) = a(x) - zeta^(2 u0) a(x^-1) is carried as an integer vector
    modulo Phi_4Q (`zvec`), and each final state is joined to one sink by an edge carrying L of its lead.
  * THE CERTIFICATE: a functional F, integer on every edge vector, such that (1) no cycle has negative weight and
    (2) every path from the initial state to the sink has weight > -F(L(launch)). Both are checked by an EXACT
    integer Bellman-Ford computation (`verify`). F is FOUND by a floating-point LP (scipy's HiGHS, the dual of the
    flow relaxation), rounded to rationals and scaled to integers; the LP's word is never used -- a rounding that
    fails exact verification is retried at a finer denominator and otherwise reported. A case with no final state
    reachable is void.

Statements checked:
  * at every one-odd-prime Q <= QMAX (powers of 2 included), every 2 <= P <= Q - 2 coprime to Q, every
    0 < u0 < P and target O, A, R with a final state: the case is certified (QMAX = 80: 95 394 cases, all certified);
  * the frame signs (Lemma 3.3): on every edge into L1 of every automaton, (-1)^lambda equals the sign of
    cos(t alpha - phi) at the frame, exactly (0 violations);
  * CONTROLS (must fire): at 8/15 u0 = 2 and 7/15 u0 = 5, target A, the flow relaxation is FEASIBLE (no certificate
    exists), and T0's connection is a walk of the automaton whose arrival vector is exactly 0; at every Q <= QMAX
    with two odd prime factors some case is uncertified.

  python sA_23_census.py           every Q <= 40 recomputed, and Q <= 80 read from sA_23_record.json (a minute or two)
  python sA_23_census.py --full    every Q <= 80 recomputed (about half an hour); rewrites sA_23_record.json
"""
import json
import os
import sys
from fractions import Fraction
from math import gcd

import numpy as np
from scipy.optimize import linprog
from sympy import Poly, cyclotomic_poly, symbols

from common import check, done, HERE

RECORD = os.path.join(HERE, 'sA_23_record.json')
_z = symbols('z')
_RED = {}
LAST = {'A': 'L2', 'O': 'L1', 'R': 'H'}
LAUNCH = [(1, 2), (-1, 2)]


# ---------------------------------------------------------------- the open arcs of A.6
def off(s, P, Q):
    """the direction of a side, in units of theta = pi/2Q: L2 = 0, L1 = Q, H = P."""
    return {'L2': 0, 'L1': Q, 'H': P}[s]


def inarc(u, lo, hi, N):
    """u in the OPEN arc (lo, hi) of Z/N, read positively."""
    return (u - lo) % N != 0 and (u - lo) % N < (hi - lo) % N


def adm(side, u, P, Q):
    """the sides a ray leaving `side` in the direction u can meet."""
    N, out = 4 * Q, []
    if side == 'L2':
        if inarc(u, P, 2 * Q, N):
            out.append('H')
        if inarc(u, 0, Q, N):
            out.append('L1')
    elif side == 'L1':
        if inarc(u, Q, 2 * Q + P, N):
            out.append('H')
        if inarc(u, 2 * Q, 3 * Q, N):
            out.append('L2')
    else:
        if inarc(u, 2 * Q + P, 4 * Q, N):
            out.append('L2')
        if inarc(u, 3 * Q, 4 * Q + P, N):
            out.append('L1')
    return out


def is_headon(u, s2, P, Q):
    """the reflection in s2 of the direction u is head-on: 2 o(s2) - u = u + 2Q."""
    return (2 * off(s2, P, Q) - 2 * u - 2 * Q) % (4 * Q) == 0


# ---------------------------------------------------------------- the automaton and its weights
def odd_primes(Q):
    m, ps, f = Q, [], 3
    while m % 2 == 0:
        m //= 2
    while m > 1:
        if m % f == 0:
            ps.append(f)
            while m % f == 0:
                m //= f
        f += 2
    return ps


def reducer(M):
    """integer matrix: the coefficient vector of z^k (k mod M) -> its residue mod Phi_M."""
    if M not in _RED:
        phi = [int(c) for c in Poly(cyclotomic_poly(M, _z), _z).all_coeffs()[::-1]]
        d = len(phi) - 1
        R = np.zeros((d, M), dtype=np.int64)
        for k in range(M):
            v = np.zeros(max(k + 1, d), dtype=np.int64)
            v[k] = 1
            for top in range(k, d - 1, -1):
                if v[top]:
                    v[top - d: top + 1] -= v[top] * np.array(phi, dtype=np.int64)
            R[:, k] = v[:d]
        _RED[M] = R
    return _RED[M]


def zvec(terms, P, Q, u0):
    """L(sum e x^a) = sum e (zeta^(Pa) - zeta^(2 u0 - Pa)) as an integer vector mod Phi_4Q (x = zeta^P)."""
    M = 4 * Q
    c = np.zeros(M, dtype=np.int64)
    for a, e in terms:
        c[(P * a) % M] += e
        c[(2 * u0 - P * a) % M] -= e
    return reducer(M) @ c


def graph(P, Q, u0):
    st0 = ('L1', 0, 1, 1, (2 * Q - u0) % (4 * Q))
    seen, order, edges, i = {st0: 0}, [st0], [], 0
    while i < len(order):
        s, h, qc, jp, u = order[i]
        for c in adm(s, u, P, Q):
            if is_headon(u, c, P, Q):
                continue
            j = jp + 1
            h2 = h + ((1 if j % 2 else -1) if c == 'H' else 0)
            t = (c, h2, (qc + (c == 'L1')) % 2, j % 2, (2 * off(c, P, Q) - u) % (4 * Q))
            assert abs(h2) <= 2 * Q + 4                     # finite: the frames stay in the window
            if t not in seen:
                seen[t] = len(order)
                order.append(t)
            w = [(2 * h + 1, 2 * (-1) ** qc), (2 * h - 1, 2 * (-1) ** qc)] if c == 'L1' else []
            edges.append((i, seen[t], w))
        i += 1
    return order, edges


def lead(state, tg):
    s, h, qc, jp, u = state
    eps, sn = (1 if jp == 0 else -1), (-1) ** qc
    return {'A': [(2 * h + eps, 2 * sn)], 'O': [], 'R': [(2 * h + 1, sn), (2 * h - 1, sn)]}[tg]


def problem(P, Q, u0, tg, G=None):
    order, edges = G or graph(P, Q, u0)
    nS = len(order)
    fin = [f for f, st in enumerate(order) if st[0] == LAST[tg]]
    if not fin:
        return None
    E = list(edges) + [(f, nS, lead(order[f], tg)) for f in fin]
    Z = np.array([zvec(w, P, Q, u0) for _, _, w in E], dtype=np.int64)
    return nS, E, Z, zvec(LAUNCH, P, Q, u0)


# ---------------------------------------------------------------- the certificate
def farkas(prob):
    """the LP that FINDS F (floating point); None iff the flow relaxation is feasible."""
    nS, E, Z, z0 = prob
    phi = Z.shape[1]
    nv = nS + 1 + phi
    A = np.zeros((len(E) + 1, nv))
    b = np.zeros(len(E) + 1)
    for k, (x, yv, _) in enumerate(E):          # -(pi_y - pi_x + F.Z_e) <= 0
        A[k, yv] -= 1
        A[k, x] += 1
        A[k, nS + 1:] = -Z[k]
    A[-1, nS] += 1                               # pi_sink - pi_launch - F.z0 <= -1
    A[-1, 0] -= 1
    A[-1, nS + 1:] = -z0
    b[-1] = -1
    bounds = [(None, None)] * (nS + 1) + [(-100, 100)] * phi
    res = linprog(np.zeros(nv), A_ub=A, b_ub=b, bounds=bounds, method='highs')
    return res.x[nS + 1:] if res.status == 0 else None


def verify(prob, y, den):
    """EXACT: the integer functional from y rounded at denominator `den`; Bellman-Ford with integer weights."""
    nS, E, Z, z0 = prob
    yr = [Fraction(float(v)).limit_denominator(den) for v in y]
    D = 1
    for f in yr:
        D = D * f.denominator // gcd(D, f.denominator)
    yi = np.array([int(f * D) for f in yr], dtype=object)
    w = [int(np.dot(yi, Z[k].astype(object))) for k in range(len(E))]
    c0 = int(np.dot(yi, z0.astype(object)))
    dist = [None] * (nS + 1)
    dist[0] = 0
    changed = False
    for _ in range(nS + 2):
        changed = False
        for k, (x, yv, _) in enumerate(E):
            if dist[x] is not None and (dist[yv] is None or dist[x] + w[k] < dist[yv]):
                dist[yv] = dist[x] + w[k]
                changed = True
        if not changed:
            break
    if changed:
        return False                              # a negative cycle
    return dist[nS] is None or dist[nS] + c0 > 0


def certify(P, Q, u0, tg, G=None):
    prob = problem(P, Q, u0, tg, G)
    if prob is None:
        return 'void'
    y = farkas(prob)
    if y is None:
        return 'feasible'
    if any(verify(prob, y, den) for den in (10 ** 6, 10 ** 9, 10 ** 12)):
        return 'certified'
    return 'unverified'


def frame_sign_violations(P, Q, u0, G):
    order, edges = G
    bad = 0
    for a, _, w in edges:
        if not w:
            continue
        s, h, qc, jp, u = order[a]
        r = (u0 - 2 * h * P) % (4 * Q)           # cos(r pi / 2Q) > 0  <=>  r in (-Q, Q) mod 4Q
        assert r not in (Q, 3 * Q)
        bad += (r < Q or r > 3 * Q) != ((-1) ** qc > 0)
    return bad


def walk_vector(word, P, Q, u0, tg):
    """T0's connection as a walk of the automaton: its total arrival vector."""
    order, edges = graph(P, Q, u0)
    emap = {(a, order[b][0]): (b, w) for a, b, w in edges}
    cur, tot = 0, zvec(LAUNCH, P, Q, u0).copy()
    for c in word[1:]:
        b, w = emap[(cur, c)]
        tot = tot + zvec(w, P, Q, u0)
        cur = b
    assert order[cur][0] == LAST[tg]
    return tot + zvec(lead(order[cur], tg), P, Q, u0)


def census_q(Q):
    row = dict(n_odd_primes=len(odd_primes(Q)), cases=0, certified=0, uncertified=[], frame_sign_violations=0)
    for P in range(2, Q - 1):
        if gcd(P, Q) != 1:
            continue
        for u0 in range(1, P):
            G = graph(P, Q, u0)
            row['frame_sign_violations'] += frame_sign_violations(P, Q, u0, G)
            for tg in 'AOR':
                r = certify(P, Q, u0, tg, G)
                if r == 'void':
                    continue
                row['cases'] += 1
                if r == 'certified':
                    row['certified'] += 1
                else:
                    row['uncertified'].append([P, u0, tg, r])
    return row


def main():
    full = '--full' in sys.argv
    qrun = 80 if full else 40
    print(f'III Proposition A.23 -- the census (Q <= {qrun} recomputed)\n')
    w815 = 'L1 H L2 H L1 H L1 L2 H L1 L2 H L2'.split()
    w715 = 'L1 H L1 L2 H L1 L2 H L2 H L1 H L2'.split()
    for P, Q, u0, w in ((8, 15, 2, w815), (7, 15, 5, w715)):
        z, c = walk_vector(w, P, Q, u0, 'A'), certify(P, Q, u0, 'A')
        check(f'CONTROL: {P}/{Q} u0 = {u0}, target A: the relaxation is {c}, and T0\'s connection is a walk with '
              f'arrival vector {"0" if not z.any() else "NONZERO"} (fires)', c == 'feasible' and not z.any())

    rows = {}
    for Q in range(3, qrun + 1):
        rows[str(Q)] = census_q(Q)
        r = rows[str(Q)]
        print(f'  Q = {Q:2d} ({r["n_odd_primes"]} odd primes): {r["certified"]}/{r["cases"]} certified', flush=True)
    if full:
        json.dump({'QMAX': 80, 'rows': rows}, open(RECORD, 'w'), indent=1)
    rec = json.load(open(RECORD))['rows']
    if not full:
        same = all(rows[q]['cases'] == rec[q]['cases'] and rows[q]['certified'] == rec[q]['certified']
                   for q in rows)
        check(f'Q <= {qrun}: the recomputation agrees with the record, Q by Q', same)
        rows = {**rec, **rows}

    one = [q for q, r in rows.items() if r['n_odd_primes'] <= 1]
    two = [q for q, r in rows.items() if r['n_odd_primes'] >= 2]
    n1 = sum(rows[q]['cases'] for q in one)
    u1 = sum(rows[q]['cases'] - rows[q]['certified'] for q in one)
    check(f'one-odd-prime Q <= 80: {n1} cases with a final state (expect 95394), {u1} uncertified',
          n1 == 95394 and u1 == 0)
    fs = sum(r['frame_sign_violations'] for q, r in rows.items() if int(q) <= qrun)
    check(f'frame signs (Lemma 3.3) on every L1 edge, Q <= {qrun}: {fs} violations', fs == 0)
    check(f'CONTROL: every two-odd-prime Q <= 80 ({", ".join(two)}) has an uncertified case (fires)',
          bool(two) and all(rows[q]['cases'] > rows[q]['certified'] for q in two))
    done()


if __name__ == '__main__':
    main()
