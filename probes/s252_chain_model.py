"""s252 part 3 -- THE GENERAL-`m` FAN TABLE, stated as a model and scored against the
`probes/s252_fan_chain.py` reads.

THE MODEL (this is the new content; `m >= 1`, `a1 >= 2m+1`).  Rung `m` is a chain of
`4m+2` O-fans ending in a head-on retrace.  Write the ASCENDING ARM

    A(m) = [ r, -P, P+r, -2P, 2P+r, -3P, ..., -mP, mP+r ]        (2m+1 entries)
         (odd  position j:  ((j-1)/2)P + r      even position j:  -(j/2)P)

then the exit-deviation list of the whole chain is the PALINDROME

    DEV(m) = A(m) ++ reverse(A(m)[0 : 2m]) ++ [0],               (4m+2 entries)

the trailing `0` being the head-on (the grazed `R` sits AT the fan centre).  The per-fan
split counts are `n_i = 2*a1 - e_i` with the DEFICIT list

    B(m) = [0, 1, 3, 5, ..., 4m-3]                               (2m entries)
    E(m) = B(m) ++ [4m, 4m-1] ++ reverse(B(m))[0 : 2m-1] ++ [-1] (4m+2 entries)

and `E(0) = [1, -1]` (§s242's two-fan table, which the general form degenerates on).

TWO CONSEQUENCES, both exact arithmetic once the table is granted:
  (C-A) `sum E(m) = 2*(2m-1)^2 + (8m-1) - 1 = 8m^2` identically, so
        `headon = (4m+2)*2*a1 - sum E(m) = 4(2m+1)*a1 - 8m^2`
        and `closure = 2*headon + 1 = 8(2m+1)*a1 + 1 - 16m^2`
        -- i.e. **the §s241 LADDER LAW IS A COROLLARY OF THE TABLE**, and its mysterious
        `-16m^2` is just twice the sum of the fan deficits.
  (C-B) the innermost fan count is `2*a1 - 4m`, so the table is well-formed exactly when
        `2*a1 - 4m >= 1`, i.e. **`a1 >= 2m+1`** -- which is precisely §s251 (4)'s measured
        warm-up boundary ("rung `m` is EXACT for `a1 >= 2m+1`").  The `-4` deficit at
        `a1 = 2m` is the chain losing that fan PAIR: 12 fans instead of 14 and depth `-2`.

PRIOR ART: as `probes/s252_fan_chain.py` -- §s238 (5) owns the m=1 deviation list
`r,P,P+r,P,r` (PROVED), §s239 (2) the six-fan table, §s242 (2) the two-fan table; no
general-`m` table exists anywhere in the ledger.  ⚠ §s251 (4) explicitly DECLINED to fit a
closed form to the `-4`/`-12` warm-up deficits (citing the s248/s249 misses); (C-B) is not
a fit -- it is a well-formedness condition read off the table, and it is falsifiable here.

⚠ INDEXING CORRECTION (found mid-session, applies to `probes/s252_fan_chain.py`'s output):
that probe labels a rung by its J-PAIR RANK from the ceiling, which is wrong whenever a
THIN pair sits inside the stack -- at `[6,2,5]` the rung it prints as `m=4` is the `m=3`
rung (12 fans, closure 189 = the `a1=2m` warm-up row §s251 (4) records).  §s242 (8b) warned
about exactly this for the unpaired ORPHAN; thin PAIRS do it too.  The intrinsic index is
`m = (n_fans - 2)/4` read off the chain, and that is what this script uses.

Run:
  PYTHONPATH=.:engine:archive/scripts_2026-07:archive/scripts_2026-08:probes \
      .venv/bin/python3.13 -u probes/s252_chain_model.py
"""
import json


def dev_model(m, P, r):
    A = [(((j - 1) // 2) * P + r) if j % 2 else (-(j // 2) * P)
         for j in range(1, 2 * m + 2)]
    return A + A[2 * m - 1::-1] + [0] if m else [r, 0]


def dev_model_warm(m, P, r):
    """The `a1 = 2m` WARM-UP chain: the palindrome with its APEX FAN PAIR missing.

    The arm turns around one step early, at `-mP`, never reaching `mP+r`, so the chain has
    `4m` fans instead of `4m+2`, the head-on lands 2 steps early and the closure is `law-4`.
    That IS §s251 (4)'s `-4` deficit, which §s251 declined to fit a form to.  Confirmed 5/5
    at `[4,2,3]` m=2, `[6,2,5]` m=3, `[6,2,7]` m=3, `[8,2,7]` m=4, `[10,2,9]` m=5.
    """
    A = [(((j - 1) // 2) * P + r) if j % 2 else (-(j // 2) * P)
         for j in range(1, 2 * m + 2)]
    return A[:2 * m] + A[2 * m - 2::-1] + [0]


def x_model(m, P, Q, r):
    """The rung's LOWER BOUNDARY `X_m`, in closed form, as a COROLLARY of the table.

    The fan centres obey `c_{i+1} = 2*rho_i - c_i` with `rho_i = c_i + s_{d_i}` (the grazed
    `R`'s height), i.e. `c_{i+1} = c_i + 2*s_{d_i}`, and `c_1 = 0`.  With F5 (the last fan's
    centre is `2*X_m`, 98/98) that gives `X_m = sum_j s_{d_j}` over the 4m+1 exit deviations.
    The palindrome `DEV(m)` repeats `A(m)[0:2m]` twice and `A(m)[2m]` once, and `s_{-n}=-s_n`:

        X_m = 2*[ sum_{k<m} s_{kP+r} - sum_{1<=k<=m} s_{kP} ] + s_{mP+r}.

    At `m=0` this is `s_r = H` and at `m=1` it is `2H - 2s_P + s_{P+r} = B1` -- the two heights
    §s242 and §s239 derived by hand.  ⚠ This is the payoff of proving the table: it is the
    closed form the UNIMODALITY residual (`f_leg_covariance.md` §21.29 (8)(b)) needs, and
    without it unimodality is a statement about measured numbers with nothing to induct on.
    """
    import mpmath as mp
    mp.mp.dps = 40
    th = mp.pi / (2 * Q)
    s_ = lambda k: mp.cot(P * th) * mp.sin(mp.pi * k / Q)
    return float(2 * (sum(s_(k * P + r) for k in range(m))
                      - sum(s_(k * P) for k in range(1, m + 1))) + s_(m * P + r))


def def_model(m):
    if m == 0:
        return [1, -1]
    B = [0] + [2 * i - 3 for i in range(2, 2 * m + 1)]
    return B + [4 * m, 4 * m - 1] + B[::-1][:2 * m - 1] + [-1]


def main():
    rows = json.load(open('data/s252_fan_chain.json'))
    tally, xtally = {}, [0, 0]
    warm, clean = [], []
    print(f'{"centre":<11} {"m_pos":>5} {"m_fan":>5} {"a1":>3} {"2m+1":>4} '
          f'{"DEV":>5} {"DEF":>5} {"ho":>6} {"pred":>6} {"cl":>6} {"pred":>6}  note')
    for w in rows:
        n = w['n_fans']
        if w['headon'] is None or (n - 2) % 4:
            warm.append((w, 'chain did not close as 4m+2 fans'))
            continue
        m = (n - 2) // 4
        P, r, a1 = w['P'], w['r'], w['a1']
        dev_ok = [int(round(x)) for x in w['devs']] == dev_model(m, P, r)
        def_ok = [2 * a1 - c for c in w['nsplits']] == def_model(m)
        ho_ok = w['headon'] == 4 * (2 * m + 1) * a1 - 8 * m * m
        cl_ok = w['closure'] == 8 * (2 * m + 1) * a1 + 1 - 16 * m * m
        inscope = a1 >= 2 * m + 1
        if inscope:                       # score X_m on the same in-scope set as DEV/DEF
            xtally[0] += abs(x_model(m, P, w['Q'], r) - w['X']) < 2e-8
            xtally[1] += 1
        rec = (dev_ok, def_ok, ho_ok, cl_ok)
        (clean if inscope else warm).append((w, rec))
        for k, v in zip(('DEV', 'DEF', 'HO', 'CL'), rec):
            key = (k, inscope)
            tally[key] = tally.get(key, [0, 0])
            tally[key][0] += bool(v)
            tally[key][1] += 1
        print(f'{str(w["digits"]):<11} {w["m"]:>5} {m:>5} {a1:>3} {2*m+1:>4} '
              f'{str(dev_ok):>5} {str(def_ok):>5} {w["headon"]:>6} '
              f'{4*(2*m+1)*a1-8*m*m:>6} {str(w["closure"]):>6} '
              f'{8*(2*m+1)*a1+1-16*m*m:>6}  '
              f'{"" if inscope else "WARM-UP a1<2m+1"}'
              f'{" m_pos MISLABEL" if m != w["m"] else ""}')

    print('\n--- score (IN SCOPE: a1 >= 2m+1) ---')
    print(f'  X_m closed form (corollary of the table): {xtally[0]}/{xtally[1]}')
    for k in ('DEV', 'DEF', 'HO', 'CL'):
        a = tally.get((k, True), [0, 0])
        b = tally.get((k, False), [0, 0])
        print(f'  {k}: in-scope {a[0]}/{a[1]}   out-of-scope(warm-up) {b[0]}/{b[1]}')

    print('\n--- (C-A) sum of deficits == 8m^2, symbolically checked m=1..40 ---')
    bad = [m for m in range(1, 41) if sum(def_model(m)) != 8 * m * m]
    print('  mismatches:', bad if bad else 'none')
    print('  => headon = (4m+2)*2a1 - 8m^2 = 4(2m+1)a1 - 8m^2, closure = 2*headon+1')

    print('\n--- (C-B) innermost fan count = 2a1 - 4m; well-formed iff a1 >= 2m+1 ---')
    for m in range(1, 6):
        E = def_model(m)
        print(f'  m={m}: E={E}  innermost deficit {max(E)} = 4m ✓' if max(E) == 4 * m
              else f'  m={m}: E={E}  innermost {max(E)} != 4m')

    print('\n--- WARM-UP / out-of-scope rows (the falsifiers, listed not hidden) ---')
    for w, rec in warm:
        m_lab = w['m'] if isinstance(rec, str) else (w['n_fans'] - 2) // 4
        print(f'  {w["digits"]} m_pos={w["m"]} a1={w["a1"]} fans={w["n_fans"]} '
              f'ho={w["headon"]} cl={w["closure"]} -> '
              f'{rec if isinstance(rec, str) else f"m_fan={m_lab} DEV/DEF/HO/CL={rec}"}')


if __name__ == '__main__':
    main()
