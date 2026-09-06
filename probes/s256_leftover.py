"""s256 -- the SECOND conjunct of §21.28 (ii), "no non-rung cell clears the floor", PROVED
from a TILING identity, at the threshold `sin(be)/sin(al)` rather than `r2/P`.

THE ROUTE.  §s254 (1) gives `V_m = R sin(Th_m)/sin(al)` and
`X_m = cos(al) R (cos Th_m + cos(be+ph))/sin^2(al)`.  Put

    Y_m := X_m + V_m        (= the rung's upper edge `U_m` when the rung is NON-degenerate;
                              in general `U_m = X_m + w_m = min(Y_m, 2 X_m)`)

Then the two §s254 identities `R cos(be+ph) = cos be - cos al`, `R sin(be+ph) = sin be` give

    (T1a)  Y_0 = 1                                (the ceiling -- §s242's capped cell)
    (T1b)  Y_m - Y_{m+1} = 2 V_m                  (<=> X_m - V_m = X_{m+1} + V_{m+1})
    (T1c)  Y_m = 1 - 2 R sin(m al) sin(m al + ph)/sin^2(al)

so rung `m`'s cell `(X_m,U_m)` plus its J-partner `(X_m-w_m, X_m)` -- equal width, §s241 (1)'s
J-pair = §s230's Mirror lemma -- tile EXACTLY `(Y_{m+1}, Y_m)`, with NO gap between consecutive
rungs.  Descending from the ceiling, the stack therefore tiles `(Y_M, 1)` for the last rung `M`,
and everything that is not a rung cell is ONE band of measure `|Y_{M+1}|`:

    non-degenerate bottom rung : the band is `(0, Y_{M+1})`, BELOW the stack;
    degenerate bottom rung     : the pair is `(0, 2X_M)` and the band is `(2X_M, Y_M)`, a MIDDLE
                                 band of measure `V_M - X_M = |Y_{M+1}|` -- where the unpaired
                                 Thm-D orphan lives (§s242 (8b)).

so the whole conjunct is ONE inequality on `|Y_{m*}|`, `m* = M+1`.  The sharp bound is NOT
`r2/P` (15 rows of 5643 exceed it, to 1.0202) but

    (H-A)  |Y_{m*}| <= sin(be)/sin(al) = (r2/P)(1 + O(a1^-2)),

proved by: (P1) POS fails at `m*` so `X_{m*} <= 0`, and `Th_{m*} >= pi-be-ph > pi/2` puts
`sin Th_{m*} <= sin(be+ph)`, whence `Y_{m*} <= V_{m*} <= R sin(be+ph)/sin al = sin be/sin al`;
(P2) `Y_{m*} = X_M - V_M >= -sin be/sin al` reduces to `cos(Th_M+al) + cos(be+ph-al) >= 0`
`<=> 2 cos((Th_M+be+ph)/2) cos((Th_M+2al-be-ph)/2) >= 0`, whose first factor is POS at `M` and
whose second needs the side condition `al <= be+ph` -- which is implied by in-scope termination
plus `a1 >= 4`: `Th_{m*} <= (a1+2)al + ph = pi/2 - rho + 2al + ph` against `Th_{m*} >= pi-be-ph`
gives `ph >= pi/4 + (rho-be)/2 - al`, so `be+ph-al >= pi/4 + (be+rho)/2 - 2al >= 0` once
`2al <= pi/a1 <= pi/4`.

WHY THIS IS NOT §21.29 (2)'s REFUTED LEFTOVER PIGEONHOLE.  That one summed to `M_floor` = the
last rung ABOVE THE FLOOR and asked `1 - 2 sum w_m < r2/P`: FALSE on 13/25 rows, `L/thr` to
2.30.  This sums to the last rung that EXISTS, subtracting the below-floor rungs too, and it is
an IDENTITY (T1b), not a pigeonhole.  §21.29 (2)'s "it holds precisely on the SATURATED rows",
i.e. when the two indices coincide, is the tell.

⚠ THE THRESHOLD MOVES, AND THAT TOUCHES A DO-NOT.  §21.28 (8) says "do not re-cut the class at
any threshold other than `r2/P`".  Its rationale (s249) is a cut BELOW the proved floor (a free
`0.5/P`), which INFLATES the class with thin seam cells and is what made the max look undecided.
`t := sin(be)/sin(al) >= r2/P` cuts ABOVE it, SHRINKING the class, and the two things the floor
was chosen for are re-checked here: T7 verifies `1/t < P/r2` (the pigeonhole is at least as
strong, so `M <= y - 1/2` is unchanged) and `w_0 >= t` (the class still starts at rung 0), and
reports `t/(r2/P)` against §21.28 (4)'s H-D1 margin 1.0401 for the charged cell.

PRIOR ART: grepped `rulings.md` [WFLOOR] (14 lines), [SELECTOR], [SEAM], [NB], [OPS];
`f_leg_covariance.md` for 'leftover', 'U_m', 'bottom region', 'telescop', 'tiling';
`cf_width_laws.md` for 'leftover', 'tile', 'J-pair', 'orphan slack', 'X_m - V_m'.  What came
back: the REFUTED leftover pigeonhole (§21.29 (2), distinguished above); `S = 2H-1` as a PROVED
endpoint ([M0-FAN], `f_leg` §21.22 (1b)) -- which IS `Y_1` at `m=0`, so (T1b)'s first case is
already proved elsewhere; §s241 (1)'s J-pair object (equal width, reversed word, H-C5 41/41).
NO prior art on the stack as a TILING, on `Y_m`, or on the leftover measured to the LAST
EXISTING rung -- new object.

⚠ SCOPE.  The closed forms are the in-scope ones, `a1 >= 2m+1` (§s252 (3)); below that the chain
loses its apex fan pair and the table is the warm-up one (`rulings.md` [WFLOOR], s255 DO-NOT v).
Rows where the five §s254 conditions still hold at the scope edge ("scope-capped") are reported
SEPARATELY: there the stack continues into warm-up rungs and this argument says nothing.

PRE-REGISTERED:
  T1  the three identities, to 1e-12, over the cone grid.
  T2  against the 20 REAL partitions of `data/s252_unimodal.json`: measured non-rung total
      measure `= |Y_{m*}|`; the cells above the band are the rung pairs in order with widths
      `w_m`; the band is where predicted; every non-rung cell is under `t`.
  T3  (H-A) over the cone grid, both against `t` (predicted 0 failures) and against `r2/P`
      (predicted to FAIL on a thin set -- record it, do not hide it).
  T4  the proof steps P1, P2 and the side condition `al <= be+ph`.
  T5  which of the 16 families binds at `m*` (predicted: POS, always).
  T6  `M_top` from the closed forms vs RING-verified presence, `data/s253_wide_scan.json`.
  T7  the threshold's consequences: `1/t < P/r2`, `w_0 >= t`, `max t/(r2/P)` vs 1.0401.

Run:
  PYTHONPATH=.:engine:archive/scripts_2026-07:archive/scripts_2026-08:probes \
      .venv/bin/python3.13 -u probes/s256_leftover.py > logs/s256_leftover.log 2>&1
"""
import json
from math import cos, gcd, pi, sin

from s254_avoidance_proof import Cell, families

OUT = 'data/s256_leftover.json'
TOL = 1e-9


def cone_centres(pmax=60, a1lo=3, a1hi=40):
    for P in range(5, pmax + 1):
        for r in range(1, P):
            if not (2 * r < P < 3 * r) or gcd(P, r) != 1:
                continue
            for a1 in range(a1lo, a1hi + 1):
                yield P, r, a1


def holds(z):
    """All 16 §s254 families at this rung (= ring existence, §s254 (7), 4000/4000)."""
    return min((s for _, _, s in families(z)), default=1.0) >= -TOL


def stack(P, r, a1):
    """(M, m*, scope_capped) -- last in-scope existing rung, first failure, cap flag."""
    Q = a1 * P + r
    mmax = (a1 - 1) // 2                      # well-formedness a1 >= 2m+1
    m = 0
    while m <= mmax and holds(Cell(P, Q, r, a1, m)):
        m += 1
    return m - 1, m, m > mmax


def clean_rows(**kw):
    """Cone rows whose stack terminates IN SCOPE (POS fails before the warm-up regime)."""
    for P, r, a1 in cone_centres(**kw):
        M, mstar, scope = stack(P, r, a1)
        if M < 0 or scope:
            continue
        yield P, r, a1, M, mstar


# --------------------------------------------------------------------------- T1
def t1_identities():
    worst, n = {'y0': 0.0, 'tel': 0.0, 'cf': 0.0}, 0
    for P, r, a1 in cone_centres(pmax=40, a1hi=25):
        Q = a1 * P + r
        z0 = Cell(P, Q, r, a1, 0)
        worst['y0'] = max(worst['y0'], abs(z0.X + z0.V - 1.0))
        for m in range((a1 - 1) // 2 + 1):
            z, zn = Cell(P, Q, r, a1, m), Cell(P, Q, r, a1, m + 1)
            worst['tel'] = max(worst['tel'],
                               abs((z.X + z.V) - (zn.X + zn.V) - 2 * z.V))
            closed = 1 - 2 * z.R * sin(m * z.al) * sin(m * z.al + z.ph) / z.sa ** 2
            worst['cf'] = max(worst['cf'], abs(z.X + z.V - closed))
            n += 1
    return worst, n


# --------------------------------------------------------------------------- T2
def t2_partitions():
    rows, out = json.load(open('data/s252_unimodal.json')), []
    for row in rows:
        P, Q, r, a1, thr = row['P'], row['Q'], row['r'], row['a1'], row['thr']
        M, mstar, scope = stack(P, r, a1)
        z = Cell(P, Q, r, a1, mstar)
        zM = Cell(P, Q, r, a1, max(M, 0))
        t = sin(z.be) / z.sa
        ystar = z.X + z.V
        cells, hi = [], 1.0
        for c in row['cells']:                       # data is sorted by descending lo
            cells.append(dict(lo=hi - c['w'], hi=hi, w=c['w'], m=c['m']))
            hi -= c['w']
        nonrung = [c for c in cells if c['m'] is None]
        meas = sum(c['w'] for c in nonrung)
        # the pairs above the band, in order, with the predicted widths
        pairs_ok, band_ok = None, None
        if not scope and M >= 0:
            want = [w for m in range(M + 1)
                    for w in (Cell(P, Q, r, a1, m).w,) * 2]
            got = [c['w'] for c in cells if c['m'] is not None]
            pairs_ok = (len(got) == len(want)
                        and max(abs(x - y) for x, y in zip(sorted(want), sorted(got))) < 2e-6)
            band = ((2 * zM.X, zM.X + zM.V) if zM.V > zM.X else (0.0, ystar))
            band_ok = all(band[0] - 1e-6 <= c['lo'] and c['hi'] <= band[1] + 1e-6
                          for c in nonrung)
        out.append(dict(digits=row['digits'], P=P, Q=Q, r=r, a1=a1, thr=thr, t=t,
                        M=M, scope=scope, ystar=ystar, meas_nonrung=meas,
                        agree=(not scope and abs(abs(ystar) - meas) < 2e-3),
                        pairs_ok=pairs_ok, band_ok=band_ok,
                        deg=bool(zM.V > zM.X), n_nonrung=len(nonrung),
                        worst_over_thr=(max([c['w'] for c in nonrung], default=0.0) / thr),
                        worst_over_t=(max([c['w'] for c in nonrung], default=0.0) / t)))
    return out


# ------------------------------------------------------------------------ T3-T5
def t3_t4_t5():
    rows, steps, binds = [], {'P1x': 0, 'P1s': 0, 'P2': 0, 'side': 0}, {}
    for P, r, a1, M, mstar in clean_rows():
        Q = a1 * P + r
        z, zM = Cell(P, Q, r, a1, mstar), Cell(P, Q, r, a1, M)
        al, be, ph = z.al, z.be, z.ph
        t, thr = sin(be) / z.sa, (P - 2 * r) / P
        ystar = z.X + z.V
        rows.append(dict(P=P, r=r, a1=a1, M=M, deg=bool(zM.V > zM.X),
                         y=ystar, over_t=abs(ystar) / t, over_thr=abs(ystar) / thr))
        tag = min(families(z), key=lambda f: f[2])[0]
        binds[tag] = binds.get(tag, 0) + 1
        steps['P1x'] += z.X > 1e-15                                   # POS fails as X <= 0
        steps['P1s'] += sin(z.Th) > sin(be + ph) + 1e-12
        steps['P2'] += cos(zM.Th + al) + cos(be + ph - al) < -1e-12
        steps['side'] += al > be + ph + 1e-12
    return rows, steps, binds


# --------------------------------------------------------------------------- T6
def t6_ring():
    seen = {}
    for row in json.load(open('data/s253_wide_scan.json')):
        seen.setdefault((row['P'], row['r'], row['a1']), set()).add(row['m'])
    bad = []
    for (P, r, a1), ms in sorted(seen.items()):
        M, _, _ = stack(P, r, a1)
        if max(ms) > M:
            bad.append(dict(P=P, r=r, a1=a1, ring_max=max(ms), M=M))
    return dict(centres=len(seen), bad=bad)


# --------------------------------------------------------------------------- T7
def t7_threshold():
    worst_ex, bad_pig, bad_w0 = {}, 0, 0
    for P, r, a1 in cone_centres():
        Q = a1 * P + r
        z0 = Cell(P, Q, r, a1, 0)
        t, thr = sin(z0.be) / z0.sa, (P - 2 * r) / P
        ex = t / thr
        for lo in (3, 4, 5):
            if a1 >= lo and ex > worst_ex.get(lo, (0,))[0]:
                worst_ex[lo] = (ex, P, r, a1)
        bad_pig += 1.0 / t >= P / (P - 2 * r)
        bad_w0 += z0.w < t - 1e-15
    return dict(worst_excess={k: list(v) for k, v in worst_ex.items()},
                bad_pigeonhole=bad_pig, bad_w0=bad_w0)


# --------------------------------------------------------------------------- T8
def t8_class_floor():
    """`w_0 >= t` is PROVED, not measured: `V_0 >= t <=> POS at m=0`, and `X_0 = H > 1/2 > t`.
    Bonus: §21.16 H-D1's measured `1-H < W_s` is the closed form `4 al + 2 ph <= pi`."""
    bad = {'V0': 0, 'H': 0, 'tbig': 0, 'hd1': 0, 'hd1_rows': 0}
    for P, r, a1 in cone_centres():
        Q = a1 * P + r
        z0, z1 = Cell(P, Q, r, a1, 0), Cell(P, Q, r, a1, 1)
        t = sin(z0.be) / z0.sa
        if not holds(z0):
            continue
        bad['V0'] += (z0.V < t - 1e-15)                       # <=> POS at m=0
        bad['H'] += (z0.X <= 0.5)                             # §s242's 2H > 1, automatic
        bad['tbig'] += (t >= 0.5)
        if holds(z1):                                         # H-D1 as a closed form
            bad['hd1_rows'] += 1
            bad['hd1'] += ((z1.V >= z0.V) != (4 * z0.al + 2 * z0.ph <= pi + 1e-12))
    return bad


def main():
    res = {}

    w, n = t1_identities()
    res['T1'] = dict(rows=n, worst=w)
    print(f"T1  identities over {n} (P,r,a1,m):  Y0-1 {w['y0']:.2e}   "
          f"telescope {w['tel']:.2e}   closed form {w['cf']:.2e}")

    part = t2_partitions()
    res['T2'] = part
    print('\nT2  against the 20 REAL partitions (data/s252_unimodal.json)')
    for p in part:
        print(f"  {str(p['digits']):<13} thr={p['thr']:.4f} M={p['M']}"
              f"{'*' if p['scope'] else ' '} deg={p['deg']:d} "
              f"|Y*|={abs(p['ystar']):.5f} meas={p['meas_nonrung']:.5f} "
              f"agree={p['agree']:d} pairs={p['pairs_ok']} band={p['band_ok']} "
              f"nonrung={p['n_nonrung']:<2} worst/thr={p['worst_over_thr']:.3f} "
              f"worst/t={p['worst_over_t']:.3f}")
    cl = [p for p in part if not p['scope']]
    print(f"  clean rows {len(cl)}/20:  |Y*| == measured non-rung total "
          f"{sum(p['agree'] for p in cl)}/{len(cl)} · pairs {sum(bool(p['pairs_ok']) for p in cl)}"
          f"/{len(cl)} · band {sum(bool(p['band_ok']) for p in cl)}/{len(cl)}")
    print(f"  every non-rung cell under t, ALL 20 rows: "
          f"{sum(p['worst_over_t'] < 1 for p in part)}/20   (under r2/P: "
          f"{sum(p['worst_over_thr'] < 1 for p in part)}/20)")

    rows, steps, binds = t3_t4_t5()
    res['T3'] = dict(rows=len(rows),
                     fail_t=sum(1 for d in rows if d['over_t'] > 1 + 1e-12),
                     fail_thr=sum(1 for d in rows if d['over_thr'] > 1),
                     max_over_t=max(d['over_t'] for d in rows),
                     max_over_thr=max(d['over_thr'] for d in rows),
                     worst_thr=sorted(rows, key=lambda d: -d['over_thr'])[:8])
    res['T4'], res['T5'] = steps, binds
    print(f"\nT3  (H-A) on {len(rows)} in-scope-terminating cone centres")
    print(f"  |Y*| <= t = sin(be)/sin(al):  {res['T3']['fail_t']} failures, "
          f"max ratio {res['T3']['max_over_t']:.6f}")
    print(f"  |Y*| <  r2/P               :  {res['T3']['fail_thr']} FAILURES, "
          f"max ratio {res['T3']['max_over_thr']:.4f}  <- the r2/P form is FALSE")
    for d in res['T3']['worst_thr'][:5]:
        print(f"     P={d['P']:<3} r={d['r']:<3} a1={d['a1']:<3} M={d['M']} "
              f"deg={d['deg']:d} |Y*|/thr={d['over_thr']:.4f}")
    print(f"\nT4  proof steps (failures):  X_m* <= 0 {steps['P1x']} · "
          f"sin(Th_m*) <= sin(be+ph) {steps['P1s']} · "
          f"cos(Th_M+al)+cos(be+ph-al) >= 0 {steps['P2']} · side al <= be+ph {steps['side']}")
    print(f"T5  binding family at m*: {binds}")

    res['T6'] = t6_ring()
    print(f"\nT6  ring cross-check, {res['T6']['centres']} scanned centres: "
          f"{len(res['T6']['bad'])} with a present rung above M")

    res['T7'] = t7_threshold()
    print(f"\nT7  threshold: pigeonhole 1/t < P/r2 failures {res['T7']['bad_pigeonhole']} · "
          f"w_0 >= t failures {res['T7']['bad_w0']}")
    for k, v in sorted(res['T7']['worst_excess'].items()):
        print(f"     max t/(r2/P) for a1>={k}: {v[0]:.5f} at P={v[1]} r={v[2]} a1={v[3]}"
              f"{'   (vs H-D1 margin 1.0401)' if k == 4 else ''}")

    res['T8'] = t8_class_floor()
    print(f"\nT8  class floor at t (failures):  V_0 >= t {res['T8']['V0']} · "
          f"X_0 = H > 1/2 {res['T8']['H']} · t < 1/2 {res['T8']['tbig']} · "
          f"H-D1 (w_1>=w_0) == (4al+2ph <= pi) {res['T8']['hd1']}/{res['T8']['hd1_rows']}")

    with open(OUT, 'w') as fh:
        json.dump(res, fh, indent=1)
    print(f'\n-> {OUT}')


if __name__ == '__main__':
    main()
