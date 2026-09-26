"""Section 4: Theorem D (Theorems 4.6 and 4.7), Proposition 4.4 and Corollary 5.1, at every
centre of the partition record.

Statements checked:
  * (Thm 4.6) an orphan branch -- one whose orbit meets a side head-on and retraces, J = id --
    exists iff P is odd, and it turns round on a side of type H when Q is odd and of type L2
    when Q is even;
  * (Thm 4.7) when P is odd the orphan branch is unique;
  * (Cor 5.1) n(P/Q) = P (mod 2), and at odd P the deficit d = Q - 1 - n satisfies d = Q (mod 2);
  * the partition is complete, tiles L1 exactly, and J is an exact width-preserving involution
    on the other branches;
  * (Prop 4.4) over one return period k = 0, ..., Q-1 the interior vertical wall-copies number
    exactly 1 when P is odd (type H at odd Q, L2 at even Q) and 0 when P is even -- the
    congruence count, at every coprime pair with 5 <= Q < 60.
That the orphan turns round at its temporal MIDPOINT is how the engine records it (the walk stops
at the first head-on hit and the orbit retraces), so it is a property of the computation, not a
check of it.  The record is re-verified: every centre it holds with fewer than 10^6 steps is
recomputed from scratch here and must agree.
"""
import math

from common import centre, check, done, record


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
    exist = unique = side = parity = deficit = engine = 0
    for r in rows:
        P, Q = r['P'], r['Q']
        k = len(r['orphans'])
        exist += (k > 0) == (P % 2 == 1)
        unique += k <= 1
        want = 'h' if Q % 2 else '2'
        side += all(s == want for s in r['orphan_head_on'])
        parity += r['n'] % 2 == P % 2
        deficit += P % 2 == 0 or (Q - 1 - r['n']) % 2 == Q % 2
        engine += r['tiles'] and r['J_ok']
    n = len(rows)
    check(f'an orphan exists iff P is odd: {exist}/{n}', exist == n)
    check(f'at most one orphan (unique when P odd): {unique}/{n}', unique == n)
    check(f'the orphan turns round on H at odd Q, on L2 at even Q: {side}/{n}', side == n)
    check(f'n = P (mod 2): {parity}/{n}', parity == n)
    check(f'odd P: d = Q - 1 - n = Q (mod 2): {deficit}/{n}', deficit == n)
    check(f'complete, tiling, J an involution: {engine}/{n}', engine == n)
    cheap = [r for r in rows if r.get('steps', 0) < 10 ** 6]
    agree = sum(1 for r in cheap
                if (lambda f: f['n'] == r['n'] and f['orphans'] == r['orphans'])(
                    centre(r['P'], r['Q'], fresh=True)))
    check(f'record re-verified from scratch at every centre under 1e6 steps: {agree}/{len(cheap)}',
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
    done()


if __name__ == '__main__':
    main()
