"""[Deposit] Section 5: launch-side helper of the tracer. Paper: "Finite perpendicular-beam partitions at irrational angles". Research name: probes/s960_launch_side.py, copied verbatim except imports (deposit names), paths (data/ and figs/ under irrational_beams/) and this paragraph. Run through run_all.py.

s960 keeper (queue s959-a: *is the right triangle's slow N(ε) the right ANGLE or the launch geometry?*): the s959 perpendicular-
beam N(ε) of one-free-angle triangle families, launched from EACH of the three sides, same instrument (`beam_families.trace_nb`).

OBJECT: `beam_families` docstring (branch = maximal launch-side interval with one word up to the first head-on return to the
launch side; widths normalised by the LAUNCH side's length).  Launching from side k = cyclically relabelling the angle tuple so that
side k becomes V0V1: side 1 (V1V2) → (a1, a2, a0), side 2 (V2V0) → (a2, a0, a1); orientation is preserved.

FAMILIES (angles at V0, V1, V2 for the s959 launch side 0; α the free angle):
  right : (π/2, π/2−α, α)   side 0 = leg RA (the paper's object), side 1 = HYPOTENUSE, side 2 = leg OR (= side 0 at π/2−α)
  pi3   : (π/3, 2π/3−α, α)   pi5 : (π/5, 4π/5−α, α)        — s959's 'fast' families, π/2 NOT in the angle lattice ℚ-part
  pi4   : (π/4, 3π/4−α, α)   pi6 : (π/6, 5π/6−α, α)        — NEW: π/2 IS a mirror direction of the unfolding group, no right angle

PRE-REGISTERED (s960, before any row ran):
  (C0) CONTROLS — (a) right side 2 at α equals right side 0 at π/2−α (branch-width multiset, to 1e-9), and differs at the same α;
       (b) rational right centres, odd P / odd Q: n(side 0) + n(side 1) = Q + 1 ([NCYL-100], an INDEPENDENT instrument there);
       (c) single-ray tracer `ray_check` agrees on random interior points of rotated tables.  If (C0) fails nothing is read.
  (L1) LAUNCH GEOMETRY — if the right triangle's slowness is the launch geometry (beam parallel to a side), the HYPOTENUSE launch
       (beam parallel to no side) is FAST like pi3/pi5 side 0, and pi3/pi5 from some other side are slow.
  (L2) ANGLE — if it is the right angle, the right triangle is slow from all three sides, pi3/pi5 fast from all three.
  (L3) LATTICE — pi4/pi6 slow (from side 0 at least) would locate it in π/2 ∈ lattice, not in the right angle.
  Score: N(1e-3), N(1e-4), N(1e-5) and the decade ratio N(1e-5)/N(1e-3), matched ε, matched budget; 'slow' = right-control range.
  ⚠ KNOWN TRAP (§96-NEG): convex over reachable decades — compare at matched depth only, no slope read as γ.
PRIOR ART: `rulings.py --grep hypotenuse / 'launch side' / c_H` → the right triangle's ⊥H beam is a studied object: [NCYL-099] the
  slope law and Diophantine driver are side-independent; [NCYL-100] n + n_H = Q + 1 (odd/odd, ring-exact); [NCYL-102]/[NCYL-130] the
  γ=1 proof does not transport to H (off-grazing / chart), γ_H open; [CPRIME-011] c_H = c′ under (H-GATE), the H GEOMETRY unmeasured
  at irrational α; [NCYL-098] side 2 = side 0 at π/2−α (relabelling).  Nothing on pi3/pi5/pi4/pi6 beyond §s959b (side 0 only).
DATA: `data_index.py --grep hypotenuse` → none; `data/beam_families*.json` (side 0, same α where reused).
PRE-REGISTERED (s960, after the main run, before the census ran) — (F) the FINITE regime and −I:
  in s960's main run finite rows occur in pi3/pi5 (6/36 family×side rows) and never in right/pi4/pi6 (0/54).  The splitting
  property: π ∈ the RATIONAL rotation group (−I available) — pi_n (V0 = π/n) has it iff n is even, right has it, x2x/x3x do not.
  Algebra: without −I a head-on return word has ODD length (an even word is a rotation, and its angle π would have to lie in the
  group).  PREDICTION: finite tables occur in odd-n pi_n and x2x/x3x, NEVER in even-n pi_n or right; every branch of a finite
  table has odd word length.  Census: `census` mode, 12 families × 24 α × 3 sides, ε = 1e-4, budget 1e7.  'finite' := no capped piece,
  dropped mass < 1e-12 (float slivers), returned mass ≥ 1 − 1e-9, N(1e-2) = N(1e-4).  Could come out otherwise: one finite even-n row refutes it.
RESULT (s960; data/launch_side.json 90 rows ε = 1e-5 budget 1e8; data/s960_launch_side_controls.json;
  data/s960_finite_census.json; data/s960_finite_scan.json 12 957 rows, α-step 0.004, ε = 1e-4, budget 1e6):
  (C0) HOLDS — (a) leg swap: counts equal 3/3, widths to 1e-9 except two 2.3e-7 branches at 3e-9 (deep-word float error), and
  the same-α comparison differs 3/3; (b) n + n_H = Q + 1 on 12/12 odd/odd centres; (c) single-ray 0/1047 mismatches.
  (L1) REFUTED — the right triangle is slow from ALL THREE sides: N(1e-5) medians 43 / 30 / 40 (leg RA / HYPOTENUSE / leg OR),
  decade ratio N(1e-5)/N(1e-3) 2.3 / 2.1 / 3.4, all 18 rows exact; pi3/pi5 fast from all three (ratios 5.6–10.1 per side).
  (L2) HOLDS at this depth.  (L3) PARTIAL — pi4/pi6 lie BETWEEN: ratio medians right 2.7 < pi4 4.2 < pi6 6.1, pi5 6.7, pi3 8.0;
  N(1e-5) 38 / 156 / 237 / 421 / 242.  The even-n ladder rises with n, the right triangle (n = 2) is its slow end.
  (F) HOLDS, control fires — finite tables on the α-scan: pi3 44, pi5 6, pi7 1, x2x 73, x3x 64 grid points (all −I-free);
  right / pi4 / pi6 / pi8: 0 of 6 558.  Every finite table returns on ODD words only.  Unregistered: pi_n never finite from
  side 1 (opposite the π/n vertex); x2x/x3x only from side 0.  Finite windows are short open intervals (pi3: 18 windows on side 0,
  longest 4 grid steps; x2x/x3x up to 32–34 steps).
  ⚠ KILLED IN THREAD (census): *mirror-pair (even-word) mass slows growth* — across families the even-word mass orders like
  speed (right 1.00, pi4 0.82 … pi12 0.54, odd 0), but WITHIN each even family Spearman(even mass, N(1e-4)) = +0.15 … +0.54,
  the wrong sign; the across-family order is confounded with n.  Right = 1.00 is Lemma T ⊕ no orphan at irrational α, not new.
  RESULT of the 24-α census: UNDER-POWERED — the pi3 finite set is a union of short open windows (width 0.003–0.006 in α near
  α ≈ 0.93–1.0), which 24 uniform samples miss; 0 finite in every pi_n, so the control could not fire.  Superseded by `scan`.
Usage: python probes/launch_side.py controls | census [NPAR] [NPROC] | scan [STEP] [NPROC] | run [EPS_MIN] [BUDGET] [NPAR] [NPROC] → data/launch_side.json
"""
import _setup  # noqa: F401  (deposit: this directory on the path, cwd = irrational_beams/)
import json
import math
import os
import random
import sys
import time
from fractions import Fraction as FR

import beam_families as B

OUT = 'data/launch_side.json'
FAM = {
    'right': ([(FR(1, 2), 0, 0), (FR(1, 2), -1, 0), (FR(0), 1, 0)], (0.05, 1.5)),
    'pi3': ([(FR(1, 3), 0, 0), (FR(2, 3), -1, 0), (FR(0), 1, 0)], (0.05, 2.0)),
    'pi5': ([(FR(1, 5), 0, 0), (FR(4, 5), -1, 0), (FR(0), 1, 0)], (0.05, 2.4)),
    'pi4': ([(FR(1, 4), 0, 0), (FR(3, 4), -1, 0), (FR(0), 1, 0)], (0.05, 2.25)),
    'pi6': ([(FR(1, 6), 0, 0), (FR(5, 6), -1, 0), (FR(0), 1, 0)], (0.05, 2.5)),
}
GRID = [10 ** (-x / 4) for x in range(4, 21)]          # 1e-1 … 1e-5


def rot(angs, k):
    return angs[k:] + angs[:k]


def table(fam, al, side):
    return B.Table(rot(FAM[fam][0], side), al)


def widths(r):
    return sorted(b[1] - b[0] for b in r['branches'])


def controls():
    out = {}
    # (a) leg swap
    rows = []
    for al in (0.3123, 0.6211, 1.0457):
        w2 = widths(B.trace_nb(table('right', al, 2), eps_min=1e-7, budget=3 * 10 ** 8))
        w0s = widths(B.trace_nb(table('right', math.pi / 2 - al, 0), eps_min=1e-7, budget=3 * 10 ** 8))
        w0 = widths(B.trace_nb(table('right', al, 0), eps_min=1e-7, budget=3 * 10 ** 8))
        eq = len(w2) == len(w0s) and max((abs(x - y) for x, y in zip(w2, w0s)), default=0) < 1e-9
        same = len(w2) == len(w0) and max((abs(x - y) for x, y in zip(w2, w0)), default=0) < 1e-9
        rows.append(dict(al=al, n2=len(w2), n0_swap=len(w0s), n0=len(w0), equal_swap=eq, equal_same_alpha=same))
        print('(a)', rows[-1], flush=True)
    out['leg_swap'] = rows
    # (b) n + n_H = Q + 1, odd P / odd Q, rational table (exact head-on test on the rational part)
    rows = []
    for P, Q in [(1, 3), (3, 5), (1, 5), (3, 7), (5, 7), (5, 9), (7, 9), (3, 11), (5, 13), (9, 13), (7, 15), (11, 17)]:
        a = FR(P, 2 * Q)
        angs = [(FR(1, 2), 0, 0), (FR(1, 2) - a, 0, 0), (a, 0, 0)]
        r0 = B.trace(B.Table(angs, 0.0), eps_min=1e-12, cap=10 ** 7)
        rH = B.trace(B.Table(rot(angs, 1), 0.0), eps_min=1e-12, cap=10 ** 7)
        n0, nH = len(r0['branches']), len(rH['branches'])
        mH = sum(b[1] - b[0] for b in rH['branches'])
        rows.append(dict(P=P, Q=Q, n0=n0, nH=nH, sum_ok=(n0 + nH == Q + 1), massH=mH, cappedH=len(rH['capped'])))
        print('(b)', rows[-1], flush=True)
    out['nH_law'] = rows
    # (c) single-ray tracer on rotated tables
    rnd = random.Random(960)
    bad = tot = 0
    for fam, side in (('right', 1), ('pi3', 1), ('pi3', 2), ('pi5', 2), ('pi4', 0), ('pi4', 1), ('pi6', 2)):
        al = rnd.uniform(*FAM[fam][1])
        T = table(fam, al, side)
        r = B.trace_nb(T, eps_min=1e-5, budget=2 * 10 ** 8)
        br = sorted(r['branches'])
        for _ in range(150):
            s = rnd.random()
            hit = [b for b in br if b[0] < s < b[1] and b[1] - b[0] > 1e-4]
            if not hit:
                continue
            n = B.ray_check(T, s)
            tot += 1
            if n != hit[0][2]:
                bad += 1
                print('MISMATCH', fam, side, al, s, n, hit[0], flush=True)
    out['rays'] = dict(checked=tot, bad=bad)
    print('(c) rays checked', tot, 'bad', bad, flush=True)
    json.dump(out, open('data/s960_launch_side_controls.json', 'w'), indent=1)
    return out


def _famdef(fam):
    if fam in FAM:
        return FAM[fam]
    if fam.startswith('pi'):
        n = int(fam[2:])
        return ([(FR(1, n), 0, 0), (FR(n - 1, n), -1, 0), (FR(0), 1, 0)], (0.05, math.pi * (1 - 1 / n) - 0.05))
    return dict(x2x=([(FR(0), 2, 0), (FR(1), -3, 0), (FR(0), 1, 0)], (0.05, 1.0)),
                x3x=([(FR(0), 3, 0), (FR(1), -4, 0), (FR(0), 1, 0)], (0.05, 0.75)))[fam]


CENSUS = ['right', 'pi3', 'pi4', 'pi5', 'pi6', 'pi7', 'pi8', 'pi9', 'pi10', 'pi12', 'x2x', 'x3x']


def _cen(job):
    fam, al, side = job[:3]
    budget = job[3] if len(job) > 3 else 10 ** 7
    angs = _famdef(fam)[0]
    T = B.Table(rot(angs, side), al)
    r = B.trace_nb(T, eps_min=1e-4, budget=budget)
    N = B.n_of_eps(r['branches'], GRID)
    mass = sum(b[1] - b[0] for b in r['branches'])
    fin = (not r['capped']) and r['dropped'] < 1e-12 and mass >= 1 - 1e-9 and N[4] == N[12]
    ks = [b[2] - 1 for b in r['branches']]                 # reflections before the head-on hit
    return dict(fam=fam, al=al, side=side, finite=bool(fin), N=N, mass=mass, steps=r['steps'],
                n_even=sum(1 for k in ks if k % 2 == 0), n_odd=sum(1 for k in ks if k % 2 == 1),
                even_mass=sum(b[1] - b[0] for b in r['branches'] if (b[2] - 1) % 2 == 0))


def census(npar, nproc, seed=9601):
    from multiprocessing import Pool
    rnd = random.Random(seed)
    jobs = []
    for fam in CENSUS:
        rng = _famdef(fam)[1]
        for _ in range(npar):
            al = rnd.uniform(*rng)
            jobs += [(fam, al, k) for k in range(3)]
    with Pool(nproc, initializer=os.nice, initargs=(10,)) as pool:
        rows = pool.map(_cen, jobs, chunksize=4)
    json.dump({'grid': GRID, 'rows': rows}, open('data/s960_finite_census.json', 'w'))
    return rows


SCAN = ['right', 'pi3', 'pi4', 'pi5', 'pi6', 'pi7', 'pi8', 'x2x', 'x3x']


def scan(step, nproc, budget=10 ** 6):
    """(F) at resolution: a uniform α-grid of step `step` over each family's range, 3 sides, budget 1e6 (finite tables end in
    ≤ 1e4 steps; a capped row is booked NOT finite).  Rows stream to data/s960_finite_scan.json (lean: no N vector)."""
    from multiprocessing import Pool
    jobs = []
    for fam in SCAN:
        lo, hi = _famdef(fam)[1]
        m = int((hi - lo) / step)
        jobs += [(fam, lo + step * (i + 0.5), k, budget) for i in range(m) for k in range(3)]
    rows = []
    with Pool(nproc, initializer=os.nice, initargs=(10,)) as pool:
        for r in pool.imap(_cen, jobs, chunksize=16):
            r.pop('N')
            rows.append(r)
            if len(rows) % 500 == 0:
                print(len(rows), '/', len(jobs), flush=True)
                json.dump({'step': step, 'budget': budget, 'rows': rows}, open('data/s960_finite_scan.json', 'w'))
    json.dump({'step': step, 'budget': budget, 'rows': rows}, open('data/s960_finite_scan.json', 'w'))
    print('DONE', len(rows), flush=True)
    return rows


def _one(job):
    fam, al, side, eps, budget = job
    os.nice(0)
    T = table(fam, al, side)
    t = time.time()
    r = B.trace_nb(T, eps_min=eps, budget=budget)
    ang = [float(q) * math.pi + m * al for (q, m, n) in rot(FAM[fam][0], side)]
    return dict(fam=fam, side=side, al=al, angles=ang, N=B.n_of_eps(r['branches'], GRID), nbr=len(r['branches']),
                mass=sum(b[1] - b[0] for b in r['branches']), capped_mass=sum(c[1] - c[0] for c in r['capped']),
                maxcapw=max((c[1] - c[0] for c in r['capped']), default=0.0), dropped=r['dropped'], steps=r['steps'],
                maxref=max((b[2] for b in r['branches']), default=0), secs=round(time.time() - t, 1))


def run(eps, budget, npar, nproc, seed=960):
    from multiprocessing import Pool
    rnd = random.Random(seed)
    jobs = []
    for fam, (_, rng) in FAM.items():
        for _ in range(npar):
            al = rnd.uniform(*rng)
            for side in range(3):
                jobs.append((fam, al, side, eps, budget))
    out = {'eps_min': eps, 'budget': budget, 'grid': GRID, 'rows': []}
    with Pool(nproc, initializer=os.nice, initargs=(10,)) as pool:
        for row in pool.imap_unordered(_one, jobs):
            out['rows'].append(row)
            print(row['fam'], row['side'], round(row['al'], 4), [row['N'][k] for k in (4, 8, 12, 16)],
                  'ret %.3f maxcapw %.1g' % (row['mass'], row['maxcapw']), row['secs'], 's', flush=True)
            json.dump(out, open(OUT, 'w'))
    return out


if __name__ == '__main__':
    if sys.argv[1] == 'controls':
        controls()
    elif sys.argv[1] == 'scan':
        scan(float(sys.argv[2]) if len(sys.argv) > 2 else 0.004, int(sys.argv[3]) if len(sys.argv) > 3 else 2)
    elif sys.argv[1] == 'census':
        census(int(sys.argv[2]) if len(sys.argv) > 2 else 24, int(sys.argv[3]) if len(sys.argv) > 3 else 2)
    else:
        eps = float(sys.argv[2]) if len(sys.argv) > 2 else 1e-5
        budget = int(float(sys.argv[3])) if len(sys.argv) > 3 else 10 ** 8
        npar = int(sys.argv[4]) if len(sys.argv) > 4 else 6
        nproc = int(sys.argv[5]) if len(sys.argv) > 5 else 2
        run(eps, budget, npar, nproc)
