"""[Deposit] Appendix B, tiling and the case 7 <= m <= 100. Paper: "Complete periodicity of perpendicular billiards in rational triangles of prime denominator". Research name: probes/s948_q3_lattice.py, copied verbatim except imports (deposit names), paths (data/ and figs/ under triangles/) and this paragraph. Run through run_all.py; the usage lines below are the research repo's. Bracketed IDs and section marks (§sNNN) refer to the project's internal ledgers.

s948 `s947-a`: the FINITE half and the COMBINATORIAL half of the q = 3 legality theorem (the m ≥ M0 continuum half is
`probes/s948_q3_legality.py`).

(T) TILING, every m ≥ 7 (exact affine identities in m): in each strand (v, res) — slots d = idx_v + 6t + 1 + 2res, t = 0 … t_end,
t_end = ⌊(2a_v − 2 − 2res)/6⌋ — the families' ranges are consecutive: lo_first = 0, hi_i + 1 = lo_{i+1}, hi_last = t_end, every
affine range non-empty for m ≥ 7; at a transcendental junction i|i+1 the joint range [lo_i, hi_{i+1}] is split by the sign of the
boundary function. Vertex families: exactly m slots each (so 2m in all).
(R) RUNS AND SPANS, every m ≥ 7: every run length is ≥ 1 (or ≡ 0) and satisfies the arc condition 3(r − 1) < N on its family's range
(affine in (m, t): checked at the range's ends for m = 7 and on the slopes; at a transcendental end the band t ≤ (y_c + η)/x).
(L) LATTICE, 7 ≤ m ≤ M0 (inclusive; the continuum half's X0 is outward too): at EVERY slot of every perpendicular family, every sign condition of its parity class holds (float64 with
an a-priori error bound Σ|c|·10⁻¹³ — arguments reduced EXACTLY mod 8N before the cosine — and an mpmath interval fallback when the
margin is smaller), r ≥ 2 gating as in `word_family`; transcendental junctions assigned by the sign of f.
PRIOR ART / DATA: as `probes/s948_q3_legality.py`. Usage: `probes/s948_q3_lattice.py F1|F2 [M0]` → `data/s948_q3_lattice_{fam}.json`."""
import json, sys, time
import _setup  # noqa: F401  (deposit: lib/ on the path, cwd = triangles/)
from fractions import Fraction as Fr
from math import lcm
import numpy as np
import mpmath as mp
from mpmath import iv
import word_family as wf
import appB_2_boundaries as G
import appB_3_legality as LG

mp.mp.dps = 30
iv.dps = 30
ABC = {'F1': [(0, 3), (6, 2), (12, 4)], 'F2': [(0, 3), (6, 1), (12, 5)]}
PARS = LG.PARS


def val(e, m, t):
    return Fr(e[0]) * m + Fr(e[1]) * t + Fr(e[2])


def tiling(famname, fams, bands):
    """exact affine checks, every m ≥ 7; returns list of problems"""
    bad = []
    av = ABC[famname]
    strands = {}
    for F in fams:
        strands.setdefault((F['v'], F['res']), []).append(F)
    for (v, res), L in sorted(strands.items()):
        am, a0 = av[v]
        # t_end = floor((2 a_v − 2 − 2 res)/6), a_v = am m + a0, am ∈ {0, 6, 12}: slope 2am/6 integer
        assert (2 * am) % 6 == 0
        tend = (Fr(2 * am, 6), Fr((2 * a0 - 2 - 2 * res) // 6))
        if L[0]['lo'] != (0, 0):
            bad.append(('first lo', v, res))
        last = L[-1]['hi']
        if last is None or last != tend:
            bad.append(('last hi', v, res, last, tend))
        for A, B in zip(L, L[1:]):
            if A['hi'] is None:
                assert B['lo'] is None and A['i'] in bands, (A['i'], B['i'])
                continue
            if B['lo'] is None or (A['hi'][0], A['hi'][1] + 1) != B['lo']:
                bad.append(('gap', A['i'], B['i']))
        for F in L:
            if F['lo'] is not None and F['hi'] is not None:
                d = (F['hi'][0] - F['lo'][0], F['hi'][1] - F['lo'][1])
                if d[0] < 0 or d[0] * 7 + d[1] < 0:
                    bad.append(('empty for some m >= 7', F['i']))
                if F['kind'] == 'vertex' and d != (1, -1):
                    bad.append(('vertex family not m slots', F['i'], d))
    return bad


def runs_spans(fams, bands, eta):
    """run lengths ≥ 1 (or ≡ 0) and 3(r − 1) < N = 18m + 9, on each perp family's range, all m ≥ 7"""
    bad = []
    for F in fams:
        if F['kind'] != 'perp':
            continue
        lo = F['lo']
        hi = F['hi']
        for r, first in F['w'].runs:
            if r == wf.E0:
                continue
            rm, rt, r0 = (Fr(x) for x in r)

            def at(tb):                                   # r at t = tb(m) as (slope, const) in m, tb = (a, b)
                return (rm + rt * tb[0], r0 + rt * tb[1])
            ends = []
            for tb, side in ((lo, 'lo'), (hi, 'hi')):
                if tb is not None:
                    ends.append(at(tb))
                else:                                       # transcendental end: t within the band, t·x ∈ [y_c − η, y_c + η]
                    B = bands[F['i'] - 1] if side == 'lo' else bands[F['i']]
                    for y in (B['yc'] - eta, B['yc'] + eta):   # t = y (36 m + 18)/π (real slope)
                        ends.append((rm + rt * y * 36 / mp.pi, r0 + rt * y * 18 / mp.pi))
            for s, c in ends:
                # r ≥ 1 for all m ≥ 7 (affine in m on the range's end; r affine in t, so the range ends suffice)
                if not (s >= 0 and s * 7 + c >= 1):
                    bad.append(('run < 1', F['i'], str(r)))
                # 3(r − 1) < 18 m + 9  ⟺  (18 − 3s) m + 9 − 3(c − 1) > 0
                if not (18 - 3 * s >= 0 and (18 - 3 * s) * 7 + 9 - 3 * (c - 1) > 0):
                    bad.append(('span', F['i'], str(r)))
    return bad


class FastF:
    """a formal sum evaluated at many (m, t) in float64 with exact argument reduction"""

    def __init__(self, f):
        items = list(f.items())
        D = 1
        for e, c in items:
            for q in e:
                D = lcm(D, Fr(q).denominator)
        self.D = D
        self.A = np.array([[int(Fr(e[0]) * D), int(Fr(e[1]) * D), int(Fr(e[2]) * D)] for e, c in items], dtype=object)
        self.c = np.array([c for e, c in items], dtype=float)
        self.sabs = float(np.abs(self.c).sum())
        self.f = f

    def __call__(self, m, t):
        """m, t int arrays (same length) → float values, error bound"""
        N = 18 * m + 9
        mod = 8 * N * self.D                                       # cos(π n/(2N D)) has period 4N·D in n·… (2π ↔ 4ND)
        tot = np.zeros(len(m))
        for (a, b, g), c in zip(self.A, self.c):
            n = (int(a) * m + int(b) * t + int(g)) % mod             # exact integer reduction (python ints in object arrays)
            tot += c * np.cos(np.pi * n.astype(np.float64) / (2.0 * N.astype(np.float64) * self.D))
        return tot, self.sabs * 1e-13


def iv_eval(f, m, t):
    N = 18 * m + 9
    return sum(c * iv.cos(iv.pi * iv.mpf(val(e, m, t).numerator) / (val(e, m, t).denominator * 2 * N)) for e, c in f.items())


def lattice(famname, fams, bands, M0):
    """every slot, 7 ≤ m < M0"""
    stats = {'slots': 0, 'conds': 0, 'fallback': 0, 'fail': []}
    cache = {}
    for F in fams:
        if F['kind'] != 'perp':
            continue
        for par in PARS:
            conds, _ = LG.cond_funcs(F['w'], par)
            cache[(F['i'], par)] = [(nm, kd, FastF(f), FastF(c), r) for nm, kd, f, c, r in conds if kd in ('opp', 'same')]
    bfun = {i: {par: FastF([cd[2] for cd in fams_by_i[i]['w'].conditions(par)[0] if cd[0] == B['name']][0]) for par in PARS}
            for i, B in bands.items()} if bands else {}
    strands = {}
    for F in fams:
        strands.setdefault((F['v'], F['res']), []).append(F)
    for m in range(7, M0 + 1):                                       # M0 inclusive (audit G17)
        # ASSIGNMENT: walk each strand upward; slot t goes to the first family F (in strand order) with t ≤ hi_F (affine hi) or,
        # at a transcendental hi, σ_F·f_F(m, t) > 0 — a partition of the strand (the last family's hi is the strand end, (T)).
        assign = {}
        for key, L in strands.items():
            tend = int(L[-1]['hi'][0] * m + L[-1]['hi'][1])
            for t in range(0, tend + 1):
                for F in L:
                    if F['hi'] is not None:
                        if t <= int(F['hi'][0] * m + F['hi'][1]):
                            break
                    else:
                        par = (m % 2, t % 2)
                        v, err = bfun[F['i']][par](np.array([m], dtype=object), np.array([t], dtype=object))
                        assert abs(v[0]) > err, ('boundary f too small', m, t)
                        if bands[F['i']]['sig_lo'][par] * v[0] > 0:
                            break
                assign.setdefault(F['i'], []).append(t)
        for F in fams:
            if F['kind'] != 'perp':
                if F['kind'] == 'vertex':                            # the vertex families receive exactly their m slots
                    lo, hi = (int(F[s][0] * m + F[s][1]) for s in ('lo', 'hi'))
                    if assign.get(F['i'], []) != list(range(lo, hi + 1)):
                        stats['fail'].append((F['i'], m, 'vertex family slots'))
                continue
            if F['lo'] is not None and assign.get(F['i']) and assign[F['i']][0] != int(F['lo'][0] * m + F['lo'][1]):
                stats['fail'].append((F['i'], m, 'assignment does not start at lo'))
            ts = np.array(assign.get(F['i'], []), dtype=object)
            if len(ts) == 0:
                continue
            for t in assign[F['i']]:                                 # runs ≥ 1 (or ≡ 0) and span, at every assigned slot (audit G34)
                for r, _ in F['w'].runs:
                    if r == wf.E0:
                        continue
                    rv = val(r, m, t)
                    stats['run_checks'] = stats.get('run_checks', 0) + 1
                    if rv < 1 or 3 * (rv - 1) >= 18 * m + 9:
                        stats['fail'].append((F['i'], m, t, 'run/span'))
            sel = np.ones(len(ts), dtype=bool)
            for par in PARS:
                idx = np.array([(m % 2, int(t) % 2) == par for t in ts]) & sel
                if not idx.any():
                    continue
                tt = ts[idx]
                mm = np.full(len(tt), m, dtype=object)
                stats['slots'] += len(tt)
                for nm, kd, ff, fc, r in cache[(F['i'], par)]:
                    keep = np.ones(len(tt), dtype=bool)
                    if r is not None:
                        keep = np.array([val(r, m, int(t)) >= 2 for t in tt])
                    if not keep.any():
                        continue
                    t2, m2 = tt[keep], mm[keep]
                    vf, ef = ff(m2, t2)
                    vc, ec = fc(m2, t2)
                    stats['conds'] += len(t2)
                    good = (vf * vc < 0) if kd == 'opp' else (vf * vc > 0)
                    safe = (np.abs(vf) > ef) & (np.abs(vc) > ec)
                    for k in np.where(~(good & safe))[0]:
                        stats['fallback'] += 1
                        a, b = iv_eval(ff.f, m, int(t2[k])), iv_eval(fc.f, m, int(t2[k]))
                        pr = a * b
                        okk = (pr.b < 0) if kd == 'opp' else (pr.a > 0)
                        if not okk:
                            stats['fail'].append((F['i'], m, int(t2[k]), nm))
        if m % 10 == 0:
            print(f'  m = {m}: slots {stats["slots"]}, conditions {stats["conds"]}, fallbacks {stats["fallback"]}, '
                  f'fails {len(stats["fail"])}', flush=True)
    return stats


def main(famname, M0):
    global fams_by_i
    leg = json.load(open(f'data/s948_q3_legality_{famname}.json'))
    eta = mp.mpf(leg['eta'])
    fam, W, trans = LG.setup(famname)
    fams = []
    for T, w in W:
        fams.append({'i': T['i'], 'v': T['vertex'], 'res': T['strand'], 'kind': T['kind'], 'w': w,
                     'lo': LG.aff(T['t_lo']), 'hi': LG.aff(T['t_hi'])})
    fams_by_i = {F['i']: F for F in fams}
    bands = {}
    for k, B in leg['bands'].items():
        bands[int(k)] = {'name': B['name'], 'yc': mp.mpf(B['yc']), 'sig_lo': {eval(p): v for p, v in B['sig_lo'].items()}}
    out = {'fam': famname, 'M0': M0}
    out['tiling_problems'] = tiling(famname, fams, bands)
    print('tiling problems:', out['tiling_problems'], flush=True)
    out['run_span_problems'] = runs_spans(fams, bands, eta)
    print('run/span problems:', out['run_span_problems'], flush=True)
    t0 = time.time()
    st = lattice(famname, fams, bands, M0)
    out['lattice'] = {k: (v if k != 'fail' else v[:50]) for k, v in st.items()}
    out['lattice']['n_fail'] = len(st['fail'])
    out['all_ok'] = not out['tiling_problems'] and not out['run_span_problems'] and not st['fail']
    print(f'lattice 7 ≤ m ≤ {M0}: {st["slots"]} slots, {st["conds"]} condition checks, {st["fallback"]} interval fallbacks, '
          f'{len(st["fail"])} failures ({time.time() - t0:.0f}s)', flush=True)
    json.dump(out, open(f'data/s948_q3_lattice_{famname}.json', 'w'), indent=1, default=str)
    print('ALL OK' if out['all_ok'] else 'FAILURES')


if __name__ == '__main__':
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 100)
