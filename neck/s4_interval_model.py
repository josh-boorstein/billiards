"""Section 4: the interval model -- Theorem 4.1, Corollary 4.2 and Remark 4.3 -- at every class
row with odd Q <= 31 (Remark 4.3's statistics on the same population).

Statements checked:
  * (Theorem 4.1) the transverse coordinate is global: laying interface t down as
    I_t = [u_t, u_t + ell_t] kite by kite, the offsets close around the cycle -- decided
    exactly; consecutive intervals share an endpoint and nest;
  * (Corollary 4.2) the palindrome ell_{-1-t} = ell_t with kite 0 the only tie, and the quotient
    path |J_i| = |cos(i P pi / Q)| (paired) or |sin(i P pi / Q)| (lone);
  * (Remark 4.3(1)) a pole prong escapes to Fix(iota) on its FIRST run on about 48% of the
    paired prongs and 23% of the lone ones -- so there is no static escape criterion;
  * (Remark 4.3(2)) the folds before escape are unbounded: the maximum exceeds 4 * 10^4 on
    paired rows and 7 * 10^4 on lone rows within the range (a 200000-step cap per prong).
"""
import math

from common import check, done
from necklace import first_run_escapes, interval_model, quotient_trace


def main():
    rows = [(P, Q, e) for Q in range(5, 32, 2) for P in range(1, Q)
            if math.gcd(P, Q) == 1 for e in (0, 1)]
    print(f'Section 4 -- the interval model, {len(rows)} class rows (odd Q <= 31)\n')
    closes = nests = pal = path = 0
    esc = {True: [0, 0], False: [0, 0]}
    fmax = {True: 0, False: 0}
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
        for k in range(Q):
            if ch.tie[k]:
                continue
            esc[lone][0] += 1
            esc[lone][1] += first_run_escapes(ch, lo, hi, k)
            end, folds, crossed, _ = quotient_trace(ch, lo, hi, k)
            if crossed and end[0] != 'CAP':
                fmax[lone] = max(fmax[lone], folds)
    n = len(rows)
    check(f'Theorem 4.1: the global coordinate closes, exactly: {closes}/{n}', closes == n)
    check(f'consecutive intervals share an endpoint and nest: {nests}/{n}', nests == n)
    check(f'Corollary 4.2: palindrome, kite 0 the only tie: {pal}/{n}', pal == n)
    check(f'Corollary 4.2: |J_i| = |cos(iP pi/Q)| paired, |sin| lone: {path}/{n}', path == n)
    pp = esc[False][1] / esc[False][0]
    pl = esc[True][1] / esc[True][0]
    check(f'Remark 4.3(1): first-run escape {esc[False][1]}/{esc[False][0]} = {100 * pp:.1f}% '
          f'paired, {esc[True][1]}/{esc[True][0]} = {100 * pl:.1f}% lone (about 48% and 23%)',
          abs(100 * pp - 48) < 1.5 and abs(100 * pl - 23) < 1.5)
    check(f'Remark 4.3(2): most folds before escape {fmax[False]} paired (> 4e4), '
          f'{fmax[True]} lone (> 7e4)', fmax[False] > 4e4 and fmax[True] > 7e4)
    done()


if __name__ == '__main__':
    main()
