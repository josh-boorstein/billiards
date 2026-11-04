"""Section 4: Theorem D (Theorems 4.6 and 4.7), Proposition 4.4 and Corollary 5.1, at every
centre of the partition record, and the record's own tallies (Appendix A.3).

Statements checked:
  * (Thm 4.6) an orphan branch -- one whose orbit meets a side head-on and retraces, J = id --
    exists iff P is odd, and it turns round on a side of type H when Q is odd and of type L2
    when Q is even;
  * (Thm 4.7) when P is odd the orphan branch is unique;
  * (Cor 5.1) n(P/Q) = P (mod 2) wherever P is even or the orphan exists (every record centre);
  * the partition is complete, tiles L1 exactly, and J is an exact width-preserving involution
    on the other branches;
  * (Prop 4.4) over one return period k = 0, ..., Q-1 the interior vertical wall-copies number
    exactly 1 when P is odd (type H at odd Q, L2 at even Q) and 0 when P is even -- the
    congruence count, at every coprime pair with 5 <= Q < 60;
  * (Appendix A.3) the record: 598 centres -- the stated ranges (P <= 9 and Q <= 60; the rows
    P = 1, 2, Q - 2, Q - 1 to Q = 51; P = 1, 2 to Q = 201), which hold 595, and the three named
    centres 11/41, 17/22, 19/24 outside them; 572 certified; the 26 open ones all have
    P in {3, 4, 6, 8} and 33 <= Q <= 59; of the certified, the 394 with P odd each have exactly
    one palindromic branch (245 of them at a Q that is not a prime power) and the 178 with P
    even have none.
That the orphan turns round at its temporal MIDPOINT is how the engine records it (the walk stops
at the first head-on hit and the orbit retraces), so it is a property of the computation, not a
check of it.  The record is re-verified: every centre it holds with fewer than 10^5 steps (461 of
the 572, a few minutes) is recomputed from scratch here and must agree; with --full, every centre
(hours).
"""
import math

from common import FULL, centre, check, done, record


def wall_copies(P, Q):
    """Interior vertical wall-copies over one return period, by type (Prop 4.4's congruences)."""
    out = {'L1': 0, 'L2': 0, 'H': 0}
    for k in range(Q):
        if k and (k * P) % Q == 0:
            out['L1'] += 1
        if (2 * k * P - Q) % (2 * Q) == 0:
            out['L2'] += 1
        if ((2 * k + 1) * P - Q) % (2 * Q) == 0:
            out['H'] += 1
    return out


def main():
    rec = record()
    rows = sorted(rec.values(), key=lambda r: (r['Q'], r['P']))
    print(f'Theorem D -- {len(rows)} centres of the partition record\n')
    exist = unique = side = parity = engine = 0
    for r in rows:
        P, Q = r['P'], r['Q']
        k = len(r['orphans'])
        exist += (k > 0) == (P % 2 == 1)
        unique += k <= 1
        want = 'h' if Q % 2 else '2'
        side += all(s == want for s in r['orphan_head_on'])
        parity += r['n'] % 2 == P % 2
        engine += r['tiles'] and r['J_ok']
    n = len(rows)
    check(f'an orphan exists iff P is odd: {exist}/{n}', exist == n)
    check(f'at most one orphan (unique when P odd): {unique}/{n}', unique == n)
    check(f'the orphan turns round on H at odd Q, on L2 at even Q: {side}/{n}', side == n)
    check(f'n = P (mod 2): {parity}/{n}', parity == n)
    check(f'complete, tiling, J an involution: {engine}/{n}', engine == n)
    cheap = rows if FULL else [r for r in rows if r.get('steps', 0) < 10 ** 5]
    agree = sum(1 for r in cheap
                if (lambda f: f['n'] == r['n'] and f['orphans'] == r['orphans'])(
                    centre(r['P'], r['Q'], fresh=True)))
    check(f'record re-verified from scratch at every centre {"" if FULL else "under 1e5 steps"}: '
          f'{agree}/{len(cheap)}',
          agree == len(cheap))
    pairs = [(P, Q) for Q in range(5, 60) for P in range(1, Q) if math.gcd(P, Q) == 1]
    ok = 0
    for P, Q in pairs:
        w = wall_copies(P, Q)
        if P % 2:
            ok += w == ({'L1': 0, 'L2': 0, 'H': 1} if Q % 2 else {'L1': 0, 'L2': 1, 'H': 0})
        else:
            ok += w == {'L1': 0, 'L2': 0, 'H': 0}
    check(f'Prop 4.4: one vertical wall-copy per return iff P odd, of the stated type, at every '
          f'coprime pair with 5 <= Q < 60: {ok}/{len(pairs)}', ok == len(pairs))

    ranges = set()
    for Q in range(2, 202):
        for P in range(1, Q):
            if math.gcd(P, Q) == 1 and ((P <= 9 and Q <= 60) or (P in (1, 2, Q - 2, Q - 1) and Q <= 51)
                                        or (P in (1, 2) and Q <= 201)):
                ranges.add((P, Q))
    named = {(11, 41), (17, 22), (19, 24)}
    have = {(r['P'], r['Q']) for r in rows}
    open_ = sorted((ranges | named) - have)
    check(f'A.3: the stated ranges hold {len(ranges)} centres, with the three named ones '
          f'{len(ranges | named)} (paper: 598)', len(ranges) == 595 and len(ranges | named) == 598
          and not named & ranges)
    check(f'A.3: {len(have)} certified (paper: 572), all inside ranges + named', len(have) == 572
          and have <= ranges | named)
    check(f'A.3: {len(open_)} open (paper: 26), all with P in {{3,4,6,8}} and 33 <= Q <= 59',
          len(open_) == 26 and all(P in (3, 4, 6, 8) and 33 <= Q <= 59 for P, Q in open_))
    odd = [r for r in rows if r['P'] % 2]
    even = [r for r in rows if r['P'] % 2 == 0]

    def prime_power(q):
        ps = {p for p in range(2, q + 1) if q % p == 0 and all(p % d for d in range(2, int(p ** .5) + 1))}
        return len(ps) == 1
    one = sum(len(r['orphans']) == 1 for r in odd)
    npp = sum(1 for r in odd if not prime_power(r['Q']))
    none = sum(len(r['orphans']) == 0 for r in even)
    check(f'A.3: odd P -- {one}/{len(odd)} with exactly one palindromic branch (paper: 394), '
          f'{npp} at a Q not a prime power (paper: 245)', one == len(odd) == 394 and npp == 245)
    check(f'A.3: even P -- {none}/{len(even)} with none (paper: 178)', none == len(even) == 178)
    done()


if __name__ == '__main__':
    main()
