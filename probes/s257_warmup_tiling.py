"""s257 -- the WARM-UP tiling: the rung stack below `a1 = 2m+1`, and the leftover bound made
UNCONDITIONAL.

WHAT §s256 LEFT.  Its tiling (`Y_m := X_m+V_m`, `Y_0 = 1`, `Y_m - Y_{m+1} = 2V_m`) and its
leftover bound `|Y_{m*}| <= t := sin(be)/sin(al)` are stated only for stacks that TERMINATE IN
SCOPE.  On 1273 of 6916 cone rows the five §s254 conditions still hold at the scope edge and the
stack continues into WARM-UP rungs (`a1 <= 2m`), where §s252 (3) knew one level (`a1 = 2m`:
`4m` fans, closure `law-4`) and §s251 (4) recorded a bare `-12` at `a1 = 2m-1` with no chain.
The §s256 proof also carried `a1 >= 4`, for the side condition `al <= be+ph`.

WHAT THIS PROBE CHECKS -- the four things that close both gaps.

 (1) THE WARM CHAIN IS THE IN-SCOPE CHAIN REFLECTED IN `m -> a1-m`.  Write `u := a1-m`.  Rung
     `m`'s chain is the §s255 hand-off recursion with TURN INDEX `tau = 2*min(m,u)+1`, and its
     turn fan -- the `(min(m,u)+1)`-th ascending odd fan, whose standard `+`-exit would be the A
     at `hi = 2a1+3-lo` -- exits instead at the ADJACENT R: one step BELOW when `m < u`
     (in-scope: the apex, §s253 (1)), one step ABOVE when `m >= u` (warm).  Hence
         #fans = min(4m+2, 4u),   DEV = dev_model(m) or dev_model_warm(u),
         closure = min( 8(2m+1)a1 + 1 - 16m^2 ,  16(a1-m)(m+1) - 3 ).
     The second branch is §s251 (4)'s warm-up deficit in closed form: at level
     `j := 2m+1-a1 >= 1` the deficit is `8j-4`, i.e. `-4, -12, -20, ...` -- §s251 declined to fit
     it and §s252 (3) derived only `j=1`.

 (2) THE CLOSED FORMS EXTEND VERBATIM.  Reading the warm chain's edges off the model gives
         X^w = c_{2u} - q(u)      (fan `2u`'s penultimate split, an R at n = 2u)
         U^w = min( c_{2u} - G(u-1),  S(u) + F(u),  2X^w )
     and three identities collapse this onto §s253/§s254:
         (I-a)  c_{2u} - q(u) = X_m          -- the REFLECTION identity (proved, one line)
         (I-b)  S(u) + F(u) - X^w = V_m      -- the turn fan's A-candidate IS the width sinusoid
         (I-c)  c_{2u} - G(u-1) - X^w = cos(2u*al) = -cos(2m*al + 2rho) =: V^A,  and
                (V^A - V_m) sin(al) = -2 cos(2m*al + rho) sin(al - rho),
                so `V^A >= V_m  <=>  2m*al + rho >= pi/2  <=>  a1 <= 2m` -- EXACTLY the warm
                regime, with equality iff `a1 = 2m`.
     ⇒ the third candidate NEVER binds, `w_m = min(V_m, X_m)` verbatim, and §s256's tiling
     identity holds at EVERY `m`, in scope or not.

 (3) THE STACK CANNOT PASS `m = a1`.  `Th_{a1} = (2a1+1)al + ph = pi + be + ph`, so
         X_{a1} = 0   and   V_{a1} = -sin(be)/sin(al) = -t,   hence   Y_{a1} = -t  EXACTLY.
     So `POS` fails at `m = a1` at the latest and the stopping index obeys `m* <= a1`.

 (4) THE SIDE CONDITION IS UNNECESSARY.  §s256's lower half needs
     `Th_M + 2al <= pi + be + ph`, i.e. `Th_{m*} <= Th_{a1}`, i.e. `m* <= a1` -- which is (3),
     automatic.  §s256 instead bounded `Th_{m*}` through in-scope termination plus `a1 >= 4`,
     which forces `al <= be+ph`; that is equivalent to `sin(2 al) >= 2 sin(2 rho)` and FAILS on
     138 of the 6916 cone rows (all `a1 <= 8`).  ⇒ `|Y_{m*}| <= t` with NO scope condition, and
     the ratio is `1` exactly iff `m* = a1`.

⚠ WHAT IS NOT PROVED.  The warm chain's own avoidance list -- the §s254 analogue -- is
VERIFIED here (T4), not derived: the families are the same ones with the turn fan's guards
changed, `B-Ie(u)` becomes the EQUALITY that defines `X_m` (the role the apex's (II) plays in
scope), and the residual is again `POS`.  So the ceiling's residual is the five §s254
inequalities PLUS the warm list; what is gone is the SCOPE, not the conditionality.

PRIOR ART: grepped `rulings.md` [WFLOOR] (all 24 lines), [SELECTOR], [NB], [OPS] for
'warm-up', '2m+1', '-12', 'apex fan pair', 'scope', 'side condition', 'be+ph';
`cf_width_laws.md` for 'warm', 'a1 = 2m', 'declined', 'deficit', 'reflection', 'a1-m';
`f_leg_covariance.md` for 'scope-capped', '1273', 'al <= be+ph', 'in-scope termination'.
What came back: [WFLOOR] "Do NOT extend the fan table below `a1 = 2m+1`: the warm-up chains lose
the apex fan pair (`4m` fans, closure `law-4`) and the s255 induction assumes the untruncated
arms" (s255 DO-NOT v) -- that ruling bars extending the table by extrapolating the IN-SCOPE
closed forms; this session extends it by READING THE WARM CHAIN first (§s256 NEXT (a) asks for
exactly this, "the same computation with the apex pair removed"), so the ruling is respected,
not overridden, and the s255 induction is re-run with the turn fan's exit flipped rather than
assumed.  [WFLOOR] "do not quote the s256 leftover theorem outside its scope (1273 of 6916)" --
the item being closed.  §s251 (4) and §s252 (3) both explicitly DECLINE to fit the deficit
sequence: (1) is not a fit, it is the fan count times `2a1` minus the deficit sum, falsifiable
on every row.  No prior art on `m -> a1-m` reflection, on `X_{a1} = 0`, or on any warm-up chain
past `a1 = 2m`.

PRE-REGISTERED:
  T1  the identities (I-a),(I-b),(I-c), `X_{a1}=0`, `V_{a1}=-t`, and §s256's `Y_0=1`,
      `Y_m-Y_{m+1}=2V_m` -- now over the FULL range `m = 0..a1`, not just `m <= (a1-1)/2`.
  T2  the unified chain model against the RING reads of `data/s257_warmup_discover.json`:
      #fans, DEV, per-fan split counts, head-on, closure, and both cell edges.
  T3  the closure law `min(law, 16(a1-m)(m+1)-3)` and the level-`j` deficit `8j-4`.
  T4  the warm avoidance list: brute-force (I)+(II) over the warm chain; the two-endpoint
      reduction; which families can fail (predicted: POS only, as in scope).
  T5  against the REAL partitions (20 rows of `data/s252_unimodal.json` + 15 of
      `data/s257_warmup_discover.json`): the full stack in position order INCLUDING warm rungs,
      the leftover band `= |Y_{m*}|`, every non-rung cell under `t`.
  T6  the cone scan: `|Y_{m*}| <= t` on all 6916 rows (s256 could score only 5643), the max
      ratio, and how many rows the dropped side condition would have excluded.
  T7  the s256 corollaries with the scope removed: `1/t < P/r2`, `w_0 > t`, `M <= a1-1`.
  T8  the `r2 = 0` EDGE, `P = 2`: `t = 0` (no gap), `2m* = Q-1` (= Thm C, `n(2/Q) = Q-1`),
      `V_m = tan(pi/2Q) sin((2m+1)pi/Q)` = the archived half-pair width, `sum 2V_m = 1`.
      A LIMIT CHECK against an independently PROVED theorem, not an extension of scope.

Run:
  PYTHONPATH=.:engine:archive/scripts_2026-07:archive/scripts_2026-08:probes \
      .venv/bin/python3.13 -u probes/s257_warmup_tiling.py > logs/s257_warmup_tiling.log 2>&1
"""
import json
from math import cos, gcd, inf, pi, sin

from s254_avoidance_proof import Cell
from s252_chain_model import dev_model, dev_model_warm

OUT = 'data/s257_warmup_tiling.json'
TOL = 1e-9


# --------------------------------------------------------------- the unified chain
def turn_data(a1, m):
    """(turn index, warm?, #fans) for rung `m`; `u = a1-m`."""
    u = a1 - m
    warm = m >= u                                   # a1 <= 2m
    tau = 2 * min(m, u) + 1
    return tau, warm, min(4 * m + 2, 4 * u)


def chain_intervals(a1, m):
    """§s255 (3)'s hand-off recursion, with the turn fan's exit R taken below/above."""
    tau, warm, _ = turn_data(a1, m)
    fans, n_lead, delta, eps, i = [], 2, +1, +1, 1
    while True:
        S = (2 * a1 + 3) if eps > 0 else (2 * a1 - 1)
        other = (S + 1 - n_lead) if (i == tau and warm) else \
                (S - 1 - n_lead) if i == tau else (S - n_lead)
        lo, hi = min(n_lead, other), max(n_lead, other)
        if i == tau:
            delta = 1 if other >= n_lead else -1
        n_ex = hi if delta > 0 else lo
        fans.append(dict(i=i, lo=lo, hi=hi, eps=eps, delta=delta, n_ex=n_ex,
                         ns=(list(range(lo, hi + 1)) if delta > 0
                             else list(range(hi, lo - 1, -1)))))
        if n_ex == 0:
            return fans
        n_R = n_ex if n_ex % 2 == 0 else n_ex - delta
        n_lead = n_R - delta * (0 if n_ex % 2 == 0 else 1) - delta
        eps = -eps
        S2 = (2 * a1 + 3) if eps > 0 else (2 * a1 - 1)
        delta = 1 if n_lead <= S2 - n_lead else -1
        i += 1
        if i > 8 * a1 + 40:
            return None


def chain_geom(z):
    """Centres and split heights of rung `z.m`'s chain, from the interval list."""
    al, sa = z.al, z.sa
    A = lambda n: sin(n * al) / sa
    lam = lambda n: cos(al) if n % 2 == 0 else 1.0
    fans = chain_intervals(z.a1, z.m)
    if fans is None:
        return None, A, lam
    c = 0.0
    for f in fans:
        f['c'] = c
        n_ex = f['n_ex']
        n_R = n_ex if n_ex % 2 == 0 else n_ex - f['delta']
        c += 2 * f['eps'] * cos(al) * A(n_R)
    return fans, A, lam


def dev_list(z):
    tau, warm, _ = turn_data(z.a1, z.m)
    return dev_model_warm(z.a1 - z.m, z.P, z.r) if warm else dev_model(z.m, z.P, z.r)


def closure_pred(a1, m):
    return min(8 * (2 * m + 1) * a1 + 1 - 16 * m * m, 16 * (a1 - m) * (m + 1) - 3)


def cone(pmax=60, a1lo=3, a1hi=40):
    for P in range(5, pmax + 1):
        for r in range(1, P):
            if not (2 * r < P < 3 * r) or gcd(P, r) != 1:
                continue
            for a1 in range(a1lo, a1hi + 1):
                yield P, r, a1


def stack_top(P, r, a1):
    """`M` = last rung (POS holds), `m*` = M+1.  `X_{a1} = 0` caps the loop identically."""
    Q = a1 * P + r
    m = 0
    while m <= a1 and Cell(P, Q, r, a1, m).X > 1e-14:
        m += 1
    return m - 1, m


# --------------------------------------------------------------------------- T1
def t1_identities():
    w = {k: 0.0 for k in ('refl', 'VB', 'VA', 'Xa1', 'Va1', 'Y0', 'tel')}
    sign_bad, eq_bad, n = 0, 0, 0
    for P, r, a1 in cone(pmax=40, a1hi=25):
        Q = a1 * P + r
        z0 = Cell(P, Q, r, a1, 0)
        w['Y0'] = max(w['Y0'], abs(z0.X + z0.V - 1.0))
        za = Cell(P, Q, r, a1, a1)
        t = sin(za.be) / za.sa
        w['Xa1'] = max(w['Xa1'], abs(za.X))
        w['Va1'] = max(w['Va1'], abs(za.V + t))
        for m in range(0, a1 + 1):
            z, zn = Cell(P, Q, r, a1, m), Cell(P, Q, r, a1, m + 1)
            u = a1 - m
            al, ro, be, ph, R, sa = z.al, z.rho, z.be, z.ph, z.R, z.sa
            w['tel'] = max(w['tel'], abs((z.X + z.V) - (zn.X + zn.V) - 2 * z.V))
            if u >= 0:
                w['refl'] = max(w['refl'], abs(z.c2(u) - z.q(u) - z.X))
                w['VB'] = max(w['VB'], abs(R * sin(2 * u * al - be - ph) / sa - z.V))
                VA = cos(2 * u * al)
                w['VA'] = max(w['VA'], abs((VA - z.V) * sa
                                           + 2 * cos(2 * m * al + ro) * sin(al - ro)))
                if (VA >= z.V - 1e-13) != (2 * m >= a1):
                    sign_bad += 1
                if 2 * m == a1 and abs(VA - z.V) > 1e-12:
                    eq_bad += 1
            n += 1
    return dict(worst=w, rows=n, sign_bad=sign_bad, eq_bad=eq_bad)


# --------------------------------------------------------------------------- T2/T3
def t2_t3_chain():
    rows = json.load(open('data/s257_warmup_discover.json'))
    sc = {k: [0, 0] for k in ('fans', 'dev', 'nspl', 'ho', 'cl', 'X', 'U', 'w')}
    warm_sc = {k: [0, 0] for k in sc}
    bad = []
    for row in rows:
        P, Q, r, a1 = row['P'], row['Q'], row['r'], row['a1']
        for c in row['cells']:
            nf = c.get('n_fans')
            if not nf or c.get('headon') is None:
                continue
            # intrinsic index: in-scope (nf-2)/4, warm a1 - nf/4 -- disambiguated by CLOSURE
            cands = []
            if (nf - 2) % 4 == 0 and a1 >= 2 * ((nf - 2) // 4) + 1:
                cands.append((nf - 2) // 4)
            if nf % 4 == 0 and nf <= 4 * a1:
                mm = a1 - nf // 4
                if a1 < 2 * mm + 1:
                    cands.append(mm)
            m = next((k for k in cands if closure_pred(a1, k) == c['closure']), None)
            if m is None:
                continue
            z = Cell(P, Q, r, a1, m)
            fans, A, lam = chain_geom(z)
            tgt = warm_sc if a1 < 2 * m + 1 else sc
            nspl = [len(f['ns']) + (0 if k == 0 else 1) for k, f in enumerate(fans)]
            ho = sum(nspl)
            # the J-PARTNER carries the reversed word, so its per-fan split counts differ
            # (already at m=0: [7,9] vs [8,8]); score `nspl` on the rung cell only.
            rung_cell = abs(z.X - c['lo']) < 2e-8
            checks = dict(
                fans=(len(fans) == nf), dev=([int(round(x)) for x in c['devs']] == dev_list(z)),
                ho=(ho == c['headon']),
                cl=(2 * ho + 1 == c['closure'] == closure_pred(a1, m)),
                X=(abs(z.X - c['lo']) < 2e-8 or abs(z.X - c['hi']) < 2e-8),
                U=(abs(z.U - c['hi']) < 2e-8 or abs(z.U - (2 * c['hi'] - c['lo'])) < 2e-8),
                w=(abs(z.w - c['w']) < 2e-8))
            if rung_cell:
                checks['nspl'] = (nspl == c['nsplits'])
            for k, v in checks.items():
                tgt[k][0] += bool(v)
                tgt[k][1] += 1
            if not all(checks.values()):
                bad.append(dict(digits=row['digits'], m=m, nf=nf,
                                fail=[k for k, v in checks.items() if not v]))
    return dict(inscope={k: v for k, v in sc.items()},
                warm={k: v for k, v in warm_sc.items()}, bad=bad[:20], n_bad=len(bad))


def t3_deficits():
    """closure = min(law, 16(a1-m)(m+1)-3); at level j = 2m+1-a1 >= 1 the deficit is 8j-4."""
    bad = 0
    for a1 in range(3, 41):
        for m in range(0, a1 + 1):
            law = 8 * (2 * m + 1) * a1 + 1 - 16 * m * m
            warm = 16 * (a1 - m) * (m + 1) - 3
            j = 2 * m + 1 - a1
            if j <= 0:
                bad += (closure_pred(a1, m) != law)
            else:
                bad += (closure_pred(a1, m) != warm) or (law - warm != 8 * j - 4)
    return bad


# --------------------------------------------------------------------------- T4
def t4_avoidance(pmax=30, a1hi=14, brute=True):
    """(I)+(II) over the WARM chain, from the closed-form cell (X_m, U_m)."""
    fam = {}
    endpt = [0, 0]
    conf, n = {'TT': 0, 'TF': 0, 'FT': 0, 'FF': 0}, 0
    binds = {}
    for P, r, a1 in cone(pmax=pmax, a1lo=3, a1hi=a1hi):
        Q = a1 * P + r
        M, _ = stack_top(P, r, a1)
        for m in range(0, a1 + 1):
            if a1 >= 2 * m + 1:
                continue                                  # in-scope: §s254 owns it
            z = Cell(P, Q, r, a1, m)
            fans, A, lam = chain_geom(z)
            if fans is None:
                continue
            n += 1
            X, U = z.X, z.U
            worst, tag = X, 'POS'                          # POS and CEIL, as in scope
            if 1.0 - U < worst:
                worst, tag = 1.0 - U, 'CEIL'
            tight = []
            for f in fans:
                c, eps = f['c'], f['eps']
                dn = (X - c) if eps > 0 else (c - U)
                df = (U - c) if eps > 0 else (c - X)
                for k, nn in enumerate(f['ns'][:-1]):
                    s = lam(nn) * A(nn) - df
                    if s < worst:
                        worst, tag = s, f"I/fan{f['i']}"
                    if abs(s) < 1e-12:
                        tight.append(f"I/fan{f['i']}" + ('=2u' if f['i'] == 2 * (a1 - m)
                                                         else '=turn' if f['i'] == 2 * (a1 - m) + 1
                                                         else ''))
                s2 = dn - lam(f['n_ex']) * A(f['n_ex'])
                if s2 < worst:
                    worst, tag = s2, f"II/fan{f['i']}"
                if abs(s2) < 1e-12:
                    tight.append(f"II/fan{f['i']}" + ('=2u' if f['i'] == 2 * (a1 - m) else ''))
                # two-endpoint reduction, per parity class
                for par in (0, 1):
                    cls = [nn for nn in f['ns'][:-1] if nn % 2 == par]
                    if len(cls) > 2:
                        vals = [lam(nn) * A(nn) for nn in cls]
                        endpt[0] += 1
                        endpt[1] += min(vals) >= min(vals[0], vals[-1]) - 1e-12
            allok = worst >= -1e-9
            exists = m <= M
            conf[('T' if allok else 'F') + ('T' if exists else 'F')] += 1
            if not allok:
                binds[tag if tag in ('POS', 'CEIL') else tag.split('fan')[0] + 'fan?'] = \
                    binds.get(tag if tag in ('POS', 'CEIL')
                              else tag.split('fan')[0] + 'fan?', 0) + 1
            for k in set(tight):
                fam[k] = fam.get(k, 0) + 1
    return dict(rows=n, endpoint=endpt, existence=conf,
                binding=dict(sorted(binds.items(), key=lambda kv: -kv[1])),
                tight=dict(sorted(fam.items(), key=lambda kv: -kv[1])[:8]))


# --------------------------------------------------------------------------- T5
def t5_partitions():
    out = []
    src = [('data/s252_unimodal.json', None), ('data/s257_warmup_discover.json', 'disc')]
    for path, kind in src:
        for row in json.load(open(path)):
            P, Q, r, a1 = row['P'], row['Q'], row['r'], row['a1']
            thr = (P - 2 * r) / P
            M, mstar = stack_top(P, r, a1)
            zs = Cell(P, Q, r, a1, mstar)
            t = sin(zs.be) / zs.sa
            ystar = zs.X + zs.V
            zM = Cell(P, Q, r, a1, max(M, 0))
            cells = []
            hi = 1.0
            for c in row['cells']:
                w = c['w'] if kind is None else c['w']
                cells.append(dict(lo=hi - w, hi=hi, w=w))
                hi -= w
            # rung widths predicted by the closed forms, each twice (the J-pair)
            want = sorted([Cell(P, Q, r, a1, k).w for k in range(M + 1)] * 2)
            got, rest = [], []
            pool = [c['w'] for c in cells]
            for x in want:
                j = min(range(len(pool)), key=lambda i: abs(pool[i] - x)) if pool else None
                if j is not None and abs(pool[j] - x) < 2e-6:
                    got.append(pool.pop(j))
            rest = pool
            band = ((2 * zM.X, zM.X + zM.V) if zM.V > zM.X else (0.0, ystar))
            out.append(dict(digits=row['digits'], P=P, a1=a1, M=M, thr=thr, t=t,
                            ystar=ystar, pairs_ok=(len(got) == len(want)),
                            leftover=sum(rest), band_meas=abs(ystar),
                            agree=abs(sum(rest) - abs(ystar)) < 3e-3,
                            n_left=len(rest),
                            worst_left_over_t=(max(rest, default=0.0) / t),
                            deg=bool(zM.V > zM.X), warm=bool(a1 < 2 * M + 1)))
    return out


# --------------------------------------------------------------------------- T6/T7
def t6_scan():
    n, fail, worst, worst_row = 0, 0, 0.0, None
    warm_rows, side_fail, eq_rows = 0, 0, 0
    worst_side = (0.0, None)
    for P, r, a1 in cone():
        Q = a1 * P + r
        M, mstar = stack_top(P, r, a1)
        if M < 0:
            continue
        z = Cell(P, Q, r, a1, mstar)
        t = sin(z.be) / z.sa
        rt = abs(z.X + z.V) / t
        n += 1
        warm_rows += (a1 < 2 * M + 1)
        side = sin(2 * z.al) >= 2 * sin(2 * z.rho)          # <=> al <= be+ph
        side_fail += (not side)
        eq_rows += (mstar == a1)
        if rt > 1 + 1e-12:
            fail += 1
        if rt > worst:
            worst, worst_row = rt, (P, r, a1, M, mstar == a1)
        if side and rt > worst_side[0]:
            worst_side = (rt, (P, r, a1, M))
    return dict(rows=n, fail=fail, max_ratio=worst, worst_row=worst_row,
                warm_rows=warm_rows, side_fail=side_fail, mstar_eq_a1=eq_rows,
                max_ratio_side_ok=worst_side[0], worst_side_row=worst_side[1])


def t7_corollaries():
    from s256_leftover import stack as s256_stack
    bad = {'pigeonhole': 0, 'w0': 0, 'Mcap': 0, 'Mfloor': 0, 'POS_vs_s254list': 0,
           's256_rows': 0}
    for P, r, a1 in cone():
        Q = a1 * P + r
        M256, _, capped = s256_stack(P, r, a1)
        if not capped:                       # s256's own scored set: the full 16-family stack
            bad['s256_rows'] += 1
            bad['POS_vs_s254list'] += (stack_top(P, r, a1)[0] != M256)
        z0 = Cell(P, Q, r, a1, 0)
        t = sin(z0.be) / z0.sa
        M, _ = stack_top(P, r, a1)
        bad['pigeonhole'] += 1.0 / t >= P / (P - 2 * r)
        bad['w0'] += z0.w < t - 1e-15
        bad['Mcap'] += M > a1 - 1
        n_class = sum(1 for m in range(M + 1) if Cell(P, Q, r, a1, m).w >= t)
        bad['Mfloor'] += 2 * n_class > P / (P - 2 * r) + 1
    return bad


# --------------------------------------------------------------------------- T8
def t8_p2_limit(Qs=(7, 9, 11, 13, 15, 21, 31, 41, 61, 81)):
    """The `r2 = 0` EDGE of the cone, `P = 2`: does the machinery land on Theorem C?

    `P=2, r=1` has `r2 = P-2r = 0`, so `be = 0`, `D = 0`, `ph = 0`, `R = 1-cos(al)` and
    `t = sin(be)/sin(al) = 0` -- the leftover bound degenerates to "NO GAP".  Predictions:
        V_m = tan(pi/2Q) sin((2m+1) pi/Q)      (the archived half-pair width, at odd index)
        m*  = a1 = (Q-1)/2,  so 2m* = Q-1      (= Thm C / [T-N2], n(2/Q) = Q-1)
        sum 2V_m = 1,  Y_{m*} = 0,  and the bottom rung is the EXACT tie V = X.
    ⚠ This is a LIMIT CHECK, not an extension: the cone needs `2r < P` strictly, and the
    s254 warm-up `(star1) sin rho <= cos al sin(al-rho)` reads `1 <= cos al` here, i.e. FALSE.
    """
    from math import tan
    out = {'rows': 0, 'thmC': 0, 't': 0.0, 'Yend': 0.0, 'sum': 0.0, 'V': 0.0,
           'multiset': 0.0, 'tie': 0.0}
    for Q in Qs:
        P, r, a1 = 2, 1, (Q - 1) // 2
        M, ms = stack_top(P, r, a1)
        zs = Cell(P, Q, r, a1, ms)
        cs = [Cell(P, Q, r, a1, m) for m in range(ms)]
        note = sorted(tan(pi / (2 * Q)) * sin(k * pi / Q) for k in range(1, a1 + 1))
        out['rows'] += 1
        out['thmC'] += (2 * ms == Q - 1)
        out['t'] = max(out['t'], abs(sin(zs.be) / zs.sa))
        out['Yend'] = max(out['Yend'], abs(zs.X + zs.V))
        out['sum'] = max(out['sum'], abs(sum(2 * z.V for z in cs) - 1.0))
        out['V'] = max(out['V'], max(abs(z.V - tan(pi / (2 * Q)) * sin((2 * m + 1) * pi / Q))
                                     for m, z in enumerate(cs)))
        out['multiset'] = max(out['multiset'],
                              max(abs(a - b) for a, b in zip(note, sorted(z.V for z in cs))))
        out['tie'] = max(out['tie'], abs(cs[-1].X - cs[-1].V))
    return out


def main():
    res = {}

    res['T1'] = t1_identities()
    w = res['T1']['worst']
    print(f"T1  identities over {res['T1']['rows']} (P,r,a1,m) rows, m = 0..a1")
    print(f"    (I-a) X_m == c2(u)-q(u)          {w['refl']:.2e}")
    print(f"    (I-b) V_m == R sin(2u al-be-ph)  {w['VB']:.2e}")
    print(f"    (I-c) (VA-V) sa == -2cos(2m al+rho)sin(al-rho)  {w['VA']:.2e}")
    print(f"          VA >= V_m  <=>  a1 <= 2m : {res['T1']['sign_bad']} mismatches; "
          f"equality at a1 = 2m: {res['T1']['eq_bad']} failures")
    print(f"    X_{{a1}} == 0  {w['Xa1']:.2e}   V_{{a1}} == -t  {w['Va1']:.2e}")
    print(f"    s256:  Y_0 == 1  {w['Y0']:.2e}   telescope  {w['tel']:.2e}  (now to m = a1)")

    res['T2'] = t2_t3_chain()
    print('\nT2  unified chain model vs the RING reads (data/s257_warmup_discover.json)')
    for lab, d in (('in-scope', res['T2']['inscope']), ('WARM    ', res['T2']['warm'])):
        print(f"    {lab}: " + '  '.join(f'{k} {v[0]}/{v[1]}' for k, v in d.items()))
    if res['T2']['bad']:
        print(f"    {res['T2']['n_bad']} failing cells: {res['T2']['bad'][:6]}")

    res['T3'] = t3_deficits()
    print(f"\nT3  closure == min(law, 16(a1-m)(m+1)-3), level-j deficit 8j-4: "
          f"{res['T3']} mismatches over a1=3..40, m=0..a1")

    res['T4'] = t4_avoidance()
    print(f"\nT4  the WARM avoidance list, brute force over {res['T4']['rows']} warm rungs")
    print(f"    two-endpoint reduction of (I): {res['T4']['endpoint'][1]}/{res['T4']['endpoint'][0]}")
    print(f"    all-hold vs ring-index existence (m <= M_POS): {res['T4']['existence']}")
    print(f"    what binds when the list fails: {res['T4']['binding']}")
    print(f"    which inequalities are TIGHT (slack 0): {res['T4']['tight']}")

    res['T5'] = t5_partitions()
    print('\nT5  against the REAL partitions (20 s252 rows + 15 s257 rows)')
    for p in res['T5']:
        print(f"    {str(p['digits']):<13} a1={p['a1']:<2} M={p['M']} "
              f"{'WARM' if p['warm'] else '    '} deg={p['deg']:d} "
              f"pairs={p['pairs_ok']!s:<5} |Y*|={p['band_meas']:.5f} "
              f"leftover={p['leftover']:.5f} agree={p['agree']:d} "
              f"n_left={p['n_left']:<2} worst/t={p['worst_left_over_t']:.3f}")
    tot = len(res['T5'])
    print(f"    pairs {sum(p['pairs_ok'] for p in res['T5'])}/{tot} · "
          f"leftover == |Y*| {sum(p['agree'] for p in res['T5'])}/{tot} · "
          f"every leftover cell under t {sum(p['worst_left_over_t'] < 1 for p in res['T5'])}/{tot}")

    res['T6'] = t6_scan()
    d = res['T6']
    print(f"\nT6  cone scan, {d['rows']} rows (s256 could score 5643 of them)")
    print(f"    |Y_m*| <= t : {d['fail']} failures, max ratio {d['max_ratio']:.6f} "
          f"at P={d['worst_row'][0]} r={d['worst_row'][1]} a1={d['worst_row'][2]} "
          f"(m*=a1: {d['worst_row'][4]})")
    print(f"    stacks continuing into warm-up: {d['warm_rows']} · m* = a1 on {d['mstar_eq_a1']}")
    print(f"    rows the DROPPED side condition al<=be+ph would exclude: {d['side_fail']}; "
          f"max ratio on the rest {d['max_ratio_side_ok']:.6f} at {d['worst_side_row']}")

    res['T7'] = t7_corollaries()
    print(f"\nT7  s256 corollaries with the scope removed (failures): {res['T7']}")

    res['T8'] = t8_p2_limit()
    d = res['T8']
    print(f"\nT8  the r2 = 0 edge, P = 2 ({d['rows']} values of Q) -- a LIMIT CHECK, not an extension")
    print(f"    2m* == Q-1 (Thm C, n(2/Q) = Q-1): {d['thmC']}/{d['rows']}   "
          f"t == 0: {d['t']:.1e}   Y_m* == 0: {d['Yend']:.1e}")
    print(f"    V_m == tan(pi/2Q) sin((2m+1)pi/Q): {d['V']:.1e}   sum 2V_m == 1: {d['sum']:.1e}")
    print(f"    width multiset == the archived note's: {d['multiset']:.1e}   "
          f"bottom rung is the exact tie V = X: {d['tie']:.1e}")

    with open(OUT, 'w') as fh:
        json.dump(res, fh, indent=1)
    print(f'\n-> {OUT}')


if __name__ == '__main__':
    main()
