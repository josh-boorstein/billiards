"""Section 4: the interval model -- Theorem 4.1 and Corollary 4.2 -- at every class row with odd
Q <= 31.

Statements checked:
  * (Theorem 4.1) the transverse coordinate is global: laying interface t down as
    I_t = [u_t, u_t + ell_t] kite by kite, the offsets close around the cycle -- decided
    exactly; consecutive intervals share an endpoint and nest;
  * (Corollary 4.2) the palindrome ell_{-1-t} = ell_t with kite 0 the only tie, and the quotient
    path |J_i| = |cos(i P pi / Q)| (paired) or |sin(i P pi / Q)| (lone).
The paper proves both at every Q >= 3; this is a finite check of the proved statements.
"""
import math

from common import check, done
from necklace import interval_model


def main():
    rows = [(P, Q, e) for Q in range(5, 32, 2) for P in range(1, Q)
            if math.gcd(P, Q) == 1 for e in (0, 1)]
    print(f'Section 4 -- the interval model, {len(rows)} class rows (odd Q <= 31)\n')
    closes = nests = pal = path = 0
    for P, Q, e in rows:
        ch, lo, hi, cl = interval_model(P, Q, e)
        lone = ch.ell[ch.m_idx].is_zero()
        closes += cl
        ok = True
        for t in range(Q):
            u = (t + 1) % Q
            shares = lo[t].eq(lo[u]) or hi[t].eq(hi[u])
            inside = ((ch.cmp(lo[t], lo[u]) <= 0 and ch.cmp(hi[u], hi[t]) <= 0)
                      or (ch.cmp(lo[u], lo[t]) <= 0 and ch.cmp(hi[t], hi[u]) <= 0))
            ok = ok and shares and inside
        nests += ok
        pal += all(ch.ell[(-1 - t) % Q].eq(ch.ell[t]) for t in range(Q)) and \
            [k for k in range(Q) if ch.tie[k]] == [0]
        m = ch.m_idx
        f = (lambda i: abs(math.sin(i * P * math.pi / Q))) if lone else \
            (lambda i: abs(math.cos(i * P * math.pi / Q)))
        path += all(abs(ch.ell[(m - i) % Q].f - f(i)) < 1e-11 for i in range(m + 1))
    n = len(rows)
    check(f'Theorem 4.1: the global coordinate closes, exactly: {closes}/{n}', closes == n)
    check(f'consecutive intervals share an endpoint and nest: {nests}/{n}', nests == n)
    check(f'Corollary 4.2: palindrome, kite 0 the only tie: {pal}/{n}', pal == n)
    check(f'Corollary 4.2: |J_i| = |cos(iP pi/Q)| paired, |sin| lone: {path}/{n}', path == n)
    done()


if __name__ == '__main__':
    main()
