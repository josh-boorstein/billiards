"""s253 part 3 -- the INEQUALITY LIST the general-`m` fan table still owes, scanned over the
`a2=2` cone in CLOSED FORM: no ring, no partition, no stored data.  The `probes/
s239_proof_ineqs.py` analogue, with `m` free.

WHAT IS BEING SCANNED.  `probes/s253_interval_law.py` states the table as: fan `i` carries
sign `eps_i = (-1)^(i+1)`, an interval `[lo_i, hi_i]` of the lattice `{nP}`, and heights
`h_i(n) = c_i + eps_i*lam(n)*A(n)` with `A(n) = sin(nP*theta)/sin(P*theta)`, `lam(n) =
cos(alpha)` for even `n` (an R) and `1` for odd `n` (an A).  §s239 (3)'s principle -- THE
EXIT RULE IS THE AVOIDANCE RULE -- says the table holds iff, for every fan,

    (I)  every NON-exit vertex of the run is on the FAR side of the cell from `c_i`:
             lam(n)*A(n) >= D_far(i)  for all n in the run except the exit
    (II) the exit vertex is on the centre's OWN side:
             lam(n_ex)*A(n_ex) <= D_near(i)

with `D_near = |c_i - (near cell edge)|` and `D_far = D_near + w_m`.  Because `A` is unimodal
in `n` (peak at `n = Q/P`) and each run's exit is its endpoint FARTHEST from `Q/P`, the
interior minimum of each parity class sits at an ENDPOINT -- so (I) is two inequalities per
fan, not `O(a1)`.  This probe locates those endpoints and scores both.

⚠ WHAT IS AND IS NOT PROVED.  (II) and the two-endpoint reduction of (I) are what a general-
`m` induction has to establish; this probe VERIFIES them, it does not derive them.  What IS
proved this session, and does not depend on this scan, is the consequence: given the table,
`w_m = min(V_m, X_m)` with `V_m*sin(P*theta) = sin((2m+1)P*theta) - cos(alpha)*sin(2(mP+r)
*theta)`, hence `V_m = (R/sin alpha)*sin((2m+1)*alpha + phi)` with `R, phi` free of `m` --
a SINUSOID sampled on an arithmetic progression lying inside `(0, pi)`, so `V_m` is unimodal.
(§s239 (7)'s `W_s*sigma_P = sigma_{3P} - (sigma_{3P+2r}+sigma_{P+2r})/2` is its `m=1` case.)

PRIOR ART: as `probes/s253_interval_law.py`.  Additionally grepped `cf_width_laws.md` for
'proof_ineqs'/'cone scan'/'(star)'; §s239 (6) owns the `m=1` cone scan (749 rows, 80
families) and its two RECORDED MISSES -- G2 ("every split outside the seven named is
off-segment") is FALSE, and G3's palindrome clause is FALSE; only the SET forms survive.
This scan is written in the surviving form: it asks `>= U_m or <= X_m`, never "off-segment".

Run:
  PYTHONPATH=.:engine:archive/scripts_2026-07:archive/scripts_2026-08:probes \
      .venv/bin/python3.13 -u probes/s253_proof_ineqs.py > logs/s253_proof_ineqs.log 2>&1
"""
import json
import math
from math import atan2, cos, gcd, pi, sin, sqrt

OUT = 'data/s253_proof_ineqs.json'


def geom(P, Q, r):
    th = pi / (2 * Q)
    return dict(th=th, al=P * th, sa=sin(P * th), r2=P - 2 * r)


def intervals(m, a1):
    out = []
    for k in range(1, m + 1):
        out.append((2 * k, 2 * a1 + 3 - 2 * k, +1, +1))
        out.append((2 * k - 1, 2 * a1 - 2 * k, -1, -1))
    out.append((2 * m + 2, 2 * a1 - 2 * m, +1, +1))
    for l in range(1, 2 * m + 2):
        if l % 2:
            out.append((2 * m + 1 - l, 2 * a1 - 2 * m - 2 + l, -1, -1))
        else:
            out.append((2 * m + 3 - l, 2 * a1 - 2 * m + l, +1, +1))
    return out


def V_closed(m, P, Q, r):
    """`V_m` in the sinusoid form -- the width when the rung is non-degenerate."""
    g = geom(P, Q, r)
    C = 1 - cos(g['al']) * cos(g['r2'] * g['th'])
    D = cos(g['al']) * sin(g['r2'] * g['th'])
    return sqrt(C * C + D * D) * sin((2 * m + 1) * g['al'] + atan2(D, C)) / g['sa']


def chain_data(m, P, Q, r, a1):
    """Centres, cell edges and the per-fan inequality data, all from the interval law."""
    g = geom(P, Q, r)
    A = lambda n: sin(n * P * g['th']) / g['sa']
    lam = lambda n: cos(g['al']) if n % 2 == 0 else 1.0
    fans, c = [], 0.0
    for (lo, hi, eps, delta) in intervals(m, a1):
        ns = list(range(lo, hi + 1)) if delta > 0 else list(range(hi, lo - 1, -1))
        n_ex = ns[-1]
        n_R = n_ex if n_ex % 2 == 0 else n_ex - delta
        fans.append(dict(lo=lo, hi=hi, eps=eps, delta=delta, ns=ns, c=c, n_ex=n_ex))
        c = c + 2 * eps * cos(g['al']) * A(n_R)
    X = fans[2 * m]['c'] + fans[2 * m]['eps'] * lam(fans[2 * m]['n_ex']) \
        * A(fans[2 * m]['n_ex'])                       # apex exit R  = X_m
    U = X + min(V_closed(m, P, Q, r), X)               # min(V_m, X_m) -> U_m = X + w_m
    return g, A, lam, fans, X, U


def scan_row(m, P, Q, r, a1):
    g, A, lam, fans, X, U = chain_data(m, P, Q, r, a1)
    out = dict(P=P, Q=Q, r=r, a1=a1, m=m, X=X, U=U, I=True, II=True,
               ENDPT=True, marginI=math.inf, marginII=math.inf)
    for f in fans:
        c, eps = f['c'], f['eps']
        Dnear = (X - c) if eps > 0 else (c - U)
        Dfar = (U - c) if eps > 0 else (c - X)
        interior = f['ns'][:-1]
        for n in interior:                              # (I)
            slack = lam(n) * A(n) - Dfar
            out['marginI'] = min(out['marginI'], slack)
            if slack < -1e-9:
                out['I'] = False
        slack2 = Dnear - lam(f['n_ex']) * A(f['n_ex'])   # (II)
        out['marginII'] = min(out['marginII'], slack2)
        if slack2 < -1e-9:
            out['II'] = False
        # the two-endpoint reduction: per parity class, is the interior min at an endpoint?
        for par in (0, 1):
            cls = [n for n in interior if n % 2 == par]
            if len(cls) > 2:
                vals = [lam(n) * A(n) for n in cls]
                if min(vals) < min(vals[0], vals[-1]) - 1e-12:
                    out['ENDPT'] = False
    return out


def main():
    rows = []
    for P in range(5, 60):
        for r in range(1, P):
            if gcd(P, r) != 1 or not (2 * r < P < 3 * r):
                continue
            for a1 in range(3, 26):
                Q = a1 * P + r
                for m in range(0, (a1 - 1) // 2 + 1):
                    rows.append(scan_row(m, P, Q, r, a1))
    with open(OUT, 'w') as fh:
        json.dump([{k: v for k, v in w.items() if k != 'ns'} for w in rows], fh)
    n = len(rows)
    print(f'=== cone scan, CLOSED FORM (no ring, no partition, no stored data)')
    print(f'    P=5..59 with gcd(P,r)=1 and 2r<P<3r; a1=3..25; m=0..floor((a1-1)/2)')
    print(f'    {n} (P,r,a1,m) rows\n')
    for k, lab in (('I', '(I)  every non-exit vertex is on the FAR side of the cell'),
                   ('II', '(II) the exit vertex is on the centre\'s own side'),
                   ('ENDPT', 'two-endpoint reduction of (I) (per parity class)')):
        print(f'  {lab:<58} {sum(1 for w in rows if w[k]):>6}/{n}')
    print(f'\n  min margin on (I)  over all rows: {min(w["marginI"] for w in rows):+.3e}')
    print(f'  min margin on (II) over all rows: {min(w["marginII"] for w in rows):+.3e}')
    for k in ('I', 'II', 'ENDPT'):
        bad = [w for w in rows if not w[k]]
        if bad:
            print(f'\n  --- {len(bad)} rows failing {k} (listed, not hidden) ---')
            for w in bad[:15]:
                print(f'    P={w["P"]} r={w["r"]} a1={w["a1"]} m={w["m"]} '
                      f'mI={w["marginI"]:+.2e} mII={w["marginII"]:+.2e}')
            bym = {}
            for w in bad:
                bym[w['m']] = bym.get(w['m'], 0) + 1
            print(f'    by m: {dict(sorted(bym.items()))}')


if __name__ == '__main__':
    main()
