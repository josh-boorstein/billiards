"""[Deposit] Theorem 5.4, step (iii), Appendix B.1: the certificate to the edge mx < π/14, m ≥ 2, controls. Paper: "Perpendicular billiards in rational polygons of odd denominator". Research name: probes/s952_s1_edge_cert.py, copied verbatim except imports (deposit names), paths (data/ and figs/ under odd_polygons/) and this paragraph. Run through run_all.py; the usage lines below are the research repo's. Bracketed IDs and section marks (§sNNN) refer to the project's internal ledgers.

s952 keeper (queue `s952-b`): s951's S₁ certificate extended to the TRUE edge θ = mx < π/14 — so the torus template of
[POLY-S1-ALL] holds on the whole region x < π/(14m), every m ≥ 2 (m = 1 is [POLY-PI-X]), and S₁(m; aπ/N) is non-CP for coprime a,
odd N > 14ma (the onset law measured in `probes/s952_s1_family.py`, pre-registered amendment H′).

NOT PRE-REGISTERED: a certificate (its hypothesis, H′, was registered and tested in s952_s1_family.py). Method = s951's
(`probes/s951_s1_allm_cert.py`, imported, unchanged) with ONE change each:
 (2D, m ≥ 4) region offset C7 = 0: θ = (1 − s)π/14, (s, ξ) ∈ [0, 1] × [0, 1/4] — the edge s = 0 is θ = π/14 for EVERY ξ. A core that
   depends on θ alone (every angle coefficient a = 0) and vanishes at θ = π/14 (exact, mod Φ_56 — s951's `corner_zero`) vanishes on the
   whole edge; a box touching s = 0 is certified by sgn·∂G/∂s > 0 (interval) — G increases from 0 into the region. Cores that
   vanish at the corner only keep s951's corner test.
 (per m, m = 2, 3) x ∈ (0, π/(14M)): s935's Taylor ⊕ Horner on (0, XA] (XA = 0.97·π/(14M), rational), then on [XA, π/(14M)] either
   interval positivity or — for an expression vanishing at π/(14M) EXACTLY (mod Φ_{56M}, z = ζ_{56M}, i = z^{14M}) — f′ < 0 by
   interval (`s952_pi_family.pos_on`).
CONTROLS: (a) the 2D pipeline at C7 = −1 (beyond the edge) flags exactly the families numerically negative at N = 14m − 1 (s951's
  control, re-run); (b) s951's own per-m (`per_m`, not `per_m_open`) at x ≤ π/(14M − 1) fails exactly where numerically negative;
  (a′) the EDGE branch disabled (s951's certify at C7 = 0) must flag; (b′) `per_m_open` with a wrong cyclotomic order K = 56M + 4 must
  fail on exactly the end-zero expressions (audit s952_audit_s1x row 14: the new branches had no control that could fail); (c) the edge is SHARP
  numerically (s952_s1_family.h_prime: 0 violations at θ = (π/14)(1 − 10⁻³), 6 at (1 + 10⁻³), m = 2, 3, 4, 6, 10).
PRIOR ART: §s951 (4) (the certificate and its four port changes), §s952 (1) (the open-end technique for Π(x)), §s952 (4) (H′).
DATA: `data/s951_s1_fams.pkl` (families), `data/s952_s1_family.json` (H′).
→ general_triangles.md §s952 (4).   Usage: python probes/s952_s1_edge_cert.py [--no-control]   → data/s952_s1_edge_cert.json
"""
import _setup  # noqa: F401  (deposit: this directory on the path, cwd = odd_polygons/)
import json
import math
import sys
import time
from fractions import Fraction as Fr

import sympy as sp
from mpmath import iv, mpf

import thmC_cert_s951 as C1
from thmC_m1 import pos_on

OUT = 'data/s952_s1_edge_cert.json'
NMIN, ASL = C1.NMIN, C1.ASL


def certify_edge(core, sgn, maxdep=26):
    corner0 = core.k == 0 and core.corner_zero()
    edge0 = corner0 and all(a == 0 for _, a, _ in core.c) and core.C7 == 0
    stack = [(mpf(0), mpf(1), mpf(0), mpf(1) / NMIN, 0)]
    nb, lowest, cinfo, nedge = 0, None, None, 0
    while stack:
        s1, s2, x1, x2, dep = stack.pop()
        nb += 1
        S, X = iv.mpf([s1, s2]), iv.mpf([x1, x2])
        if edge0 and s1 == 0:
            ds, _ = core.dG(S, X)
            if (sgn * ds).a > 0:
                nedge += 1
                continue
        elif corner0 and s1 == 0 and x1 == 0:
            ds, dx = core.dG(S, X)
            Aa, Bb = (sgn * ds).a, (sgn * dx).a
            if Aa > 0 and Bb > 0:
                cinfo = dict(s2=float(s2), x2=float(x2), A=float(Aa), B=float(Bb))
                continue
        else:
            lo = (sgn * core.G(S, X)).a
            if core.k == 0 and lo <= 0:
                Sc, Xc = iv.mpf((s1 + s2) / 2), iv.mpf((x1 + x2) / 2)
                ds, dx = core.dG(S, X)
                lo = max(lo, (sgn * (core.G(Sc, Xc) + ds * (S - Sc) + dx * (X - Xc))).a)
            if lo > 0:
                lowest = lo if lowest is None or lo < lowest else lowest
                continue
        if dep > maxdep:
            return dict(ok=False, boxes=nb, where=[float(s1), float(s2), float(x1), float(x2)])
        if (s2 - s1) >= NMIN * (x2 - x1):
            m = (s1 + s2) / 2
            stack += [(s1, m, x1, x2, dep + 1), (m, s2, x1, x2, dep + 1)]
        else:
            m = (x1 + x2) / 2
            stack += [(s1, s2, x1, m, dep + 1), (s1, s2, m, x2, dep + 1)]
    return dict(ok=True, boxes=nb, lb=float(lowest) if lowest is not None else None, corner=cinfo, corner_zero=corner0,
                edge_zero=edge0, edge_boxes=nedge)


def zero_at(f, K):
    """EXACT: f(z) = 0 at z = ζ_K (x = 4π/K... i.e. z = e^{2πi/K}), i = z^{K/4}"""
    assert K % 4 == 0
    t = sp.Symbol('t')
    c = [0] * K
    for b, v in f.c.items():
        c[b % K] += int(v.real)
        c[(b + K // 4) % K] += int(v.imag)
    P = sum(ci * t ** k for k, ci in enumerate(c) if ci)
    return sp.rem(sp.Poly(P, t), sp.Poly(sp.cyclotomic_poly(K, t), t)).is_zero


def per_m_open(M, K=None):
    """every family at m = M, every repetition j, certified > 0 on x ∈ (0, π/(14M)) (open end)"""
    import pi_n_torus as KP
    z, v, u = C1.z, C1.v, C1.u
    XE = math.pi / (ASL * M)
    KP.XMAX = Fr(math.floor(0.97 * XE * 10 ** 7), 10 ** 7)
    XA = iv.mpf(KP.XMAX.numerator) / iv.mpf(KP.XMAX.denominator)
    XEiv = iv.pi / (ASL * M)
    K = K or 56 * M                                  # x = π/(14M) ⟺ z = e^{iπ/(28M)} = ζ_{56M} (K override = control)

    def toL(e):
        e = sp.cancel(sp.together(e))
        num, den = sp.fraction(e)
        dP = sp.Poly(den, z)
        assert len(dP.terms()) == 1, den
        (dd,), dc = dP.terms()[0]
        co = {int(k[0]) - int(dd): sp.nsimplify(c / dc) for k, c in sp.Poly(sp.expand(num), z).terms()}
        dl = sp.ilcm(*[sp.fraction(sp.nsimplify(sp.re(c)))[1] for c in co.values()],
                     *[sp.fraction(sp.nsimplify(sp.im(c)))[1] for c in co.values()], 1)
        return KP.L({k: complex(int(sp.re(c) * dl), int(sp.im(c) * dl)) for k, c in co.items()})

    res = dict(M=M, XE=XE, XA=str(KP.XMAX), tight=0, certified=0, end_zero=0, fails=[], boxes=0)
    x0 = math.pi / (ASL * M + 3)
    for f in C1.FAMS:
        c, d = f[1], f[2]
        for j in (range(1) if c == (0, 1) else range(c[0] * M + c[1])):
            sub = {v: z ** (2 * M), u: z ** (4 * d * j)}
            if f[4] is None:
                exprs = [('sign', f[3])]
            else:
                p, q = toL(f[5].subs(sub)), toL(f[6].subs(sub))
                exprs = ([('a-p', f[3] - f[5]), ('q-b', f[6] - f[4])] if p.evx(x0).real <= q.evx(x0).real
                         else [('a-q', f[3] - f[6]), ('p-b', f[5] - f[4])])
            for tag, e in exprs:
                L = toL(e.subs(sub))
                if L.iszero():
                    res['tight'] += 1
                    continue
                try:
                    KP.certify_nonneg(L)
                    if zero_at(L, K):
                        res['end_zero'] += 1
                        res['boxes'] += pos_on(L, XA, XEiv, deriv=1, sign=-1)
                    else:
                        res['boxes'] += pos_on(L, XA, XEiv)
                    res['certified'] += 1
                except AssertionError as ex:
                    res['fails'].append((f[0], tag, j, str(ex)[:60]))
    return res


def main():
    t0 = time.time()
    iv.dps = 40
    orig_certify = C1.certify
    C1.certify = certify_edge                         # run() calls the module-level name
    rel = C1.orient()
    main_ = C1.run(0, rel)
    edge = [r for r in main_.get('corner_tests', [])]
    print('MAIN C7=0', {k: v for k, v in main_.items() if k not in ('why', 'tight', 'corner_tests', 'zero_names')}, flush=True)
    out = {'main_C7_0': main_}
    out['per_m'] = {M: per_m_open(M) for M in (2, 3)}
    print('PER-m open', {M: {k: r[k] for k in ('certified', 'tight', 'end_zero', 'fails')} for M, r in out['per_m'].items()}, flush=True)
    if '--no-control' not in sys.argv:
        ctrl = C1.run(-1, rel)
        violm = C1.numeric_violations(rel, -1)
        out['control_2d'] = dict(flagged=ctrl['flagged'], numeric_neg_14m_minus1=sorted(violm), exact=set(ctrl['flagged']) == violm and bool(violm))
        pc = {M: C1.per_m(M, off=-1) for M in (2, 3)}
        out['control_per_m'] = {M: dict(fails=len(r['fails']), exact=r['fails_eq_numneg']) for M, r in pc.items()}
        # audit s952_audit_s1x row 14 (G): controls of the two NEW branches — each must FAIL when disabled
        C1.certify = orig_certify
        off = C1.run(0, rel)
        C1.certify = certify_edge
        out['control_edge_off'] = dict(flagged=off['flagged'], fires=bool(off['flagged']))
        out['control_end_wrongK'] = {M: dict(fails=len(r['fails']), end_zero=r['end_zero'],
                                             fires=len(r['fails']) == out['per_m'][M]['end_zero'] and r['end_zero'] == 0)
                                     for M, r in ((M, per_m_open(M, K=56 * M + 4)) for M in (2, 3))}
        print('CONTROLS', out['control_2d']['exact'], out['control_per_m'], out['control_edge_off']['fires'],
              out['control_end_wrongK'], flush=True)
    out['verdict'] = bool(not main_['flagged'] and main_['cores_ok'] == main_['distinct_cores']
                          and all(not r['fails'] for r in out['per_m'].values())
                          and (out.get('control_2d') is None or out['control_2d']['exact'])
                          and (out.get('control_per_m') is None or all(c['exact'] and c['fails'] for c in out['control_per_m'].values()))
                          and (out.get('control_edge_off') is None or out['control_edge_off']['fires'])
                          and (out.get('control_end_wrongK') is None or all(c['fires'] for c in out['control_end_wrongK'].values())))
    out['secs'] = round(time.time() - t0)
    json.dump(out, open(OUT, 'w'), indent=1, default=str)
    print('VERDICT', out['verdict'], out['secs'], 's')


if __name__ == '__main__':
    main()
