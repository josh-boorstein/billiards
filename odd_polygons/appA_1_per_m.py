"""[Deposit] Proposition A.1: the per-m certificates, m ≤ 6. Paper: "Perpendicular billiards in rational polygons of odd denominator". Research name: probes/s935_ladder_perm_cert.py, copied verbatim except imports (deposit names), paths (data/ and figs/ under odd_polygons/) and this paragraph. Run through run_all.py; the usage lines below are the research repo's. Bracketed IDs and section marks (§sNNN) refer to the project's internal ledgers.

s935 keeper, promoted from scratch in s958 (queue `s955-f`): the RIGOROUS per-m certificate of Π_N's m-th torus, m = 0 … 6
(paper 2, Prop A.1 `thm:ladder`, "per m for m ≤ 6"; the m ≥ 7 half is `probes/s937_ladder_allm_cert.py`).

For each torus m of the N = 63 template (`data/s935_tori_63.json`, written by `probes/s935_torus_ladder.py`; byte-identical to the
scratch `data/scratch/s935/lower/tori_63.json` this was run on in s935), valid for every odd N with x = π/N ≤ π/(8m + 15):
combinatorics of the two words (consecutive copies glue, both return to the start copy), the two returns EXACT (`iszero` in the
ring of `s935_pi_n_torus`), the rotation number 4 sin²(π/N) (float, at N = 8m + 15 and 101; the symbolic proof is s935_torus_ladder's),
and every containment inequality certified non-negative on (0, XMAX], XMAX = π/(8m + 15) rounded up (`s935_pi_n_torus.certify_nonneg`).
CONTROL: at N = 8m + 13, just outside the range, some inequality must be violated (`control_viol > 0`).
s935 record (scratch log `perm_cert_3_6.log`, `tori_63_cert_m012.json`): m = 0 … 6 all certified; tight 18 (m = 0), 16 (m ≥ 1);
min lower bound 0.207 (m = 0) … 0.0476 (m = 6); control violations 19, …, 91.
Copied verbatim from `data/scratch/s935/lower/perm_cert.py` except: this docstring, no absolute sys.path, the default input path, the
output path. PRIOR ART: §s935 (4) (the scratch run); DATA: the scratch records above.
Usage: python probes/s935_ladder_perm_cert.py [data/s935_tori_63.json [m,m,...]] → data/s935_ladder_perm_cert.json
"""
import _setup  # noqa: F401  (deposit: this directory on the path, cwd = odd_polygons/)
import sys, json, math
import pi_n_torus as KP
from pi_n_torus import L, X, g_apply, g_across, side_range, shift, V, SIDES, sinL, c2
from fractions import Fraction as Fr
from mpmath import iv
iv.dps = 40
IN = sys.argv[1] if len(sys.argv) > 1 else 'data/s935_tori_63.json'
D = json.load(open(IN)); N0 = D['N']; sc = D['scale']
res = {}
for m, t in enumerate(D['tori']):
    if len(sys.argv) > 2 and m not in [int(a) for a in sys.argv[2].split(',')]: continue
    h0, h1 = t['h']; rows = t['rows']
    st = [(r['o'], r['e'], r['s'], r['s2']) for r in rows]
    J0, J1 = st[:h0], st[h0:]
    def walk(seg):
        S, out = L(), []
        for (o, e, s, s2) in seg:
            out.append(S); S = S + shift((o, e), s2)
        return out, S
    def pin(rs, seg, hi):
        Ss, _ = walk(seg)
        for k, r in enumerate(rs):
            u = (r['u'] + (r['w'] if hi else 0)) / sc
            for vi in range(4):
                if abs(r['vx'][vi] / sc - u) < 1e-9:
                    return X(g_apply(seg[k][:2], V[vi])) + Ss[k], (k, vi)
        raise ValueError('no pin')
    lo0, p0 = pin(rows[:h0], J0, False); hi0, p1 = pin(rows[:h0], J0, True); hi1, p2 = pin(rows[h0:], J1, True)
    lam0, lam1 = hi0 - lo0, hi1 - hi0
    g0 = J0[0][:2]; s0 = J0[0][2]
    ok = True
    for seg in (J0, J1):
        for a, b in zip(seg, seg[1:]):
            ok &= g_across(a[:2], a[3]) == b[:2] and a[3] == b[2]
        ok &= g_across(seg[-1][:2], seg[-1][3]) == g0 and seg[-1][3] == s0
        ok &= all(not ((o, e) == g0 and s == s0) for (o, e, s, _) in seg[1:])
    ok &= J1[0][:3] == J0[0][:3]
    _, S0 = walk(J0); _, S1 = walk(J1)
    ret = ((lo0 - S0) - (lo0 + lam1)).iszero() and ((hi0 - S1) - lo0).iszero()
    rho = [lam1.evx(math.pi / n).real / (lam0 + lam1).evx(math.pi / n).real for n in (8 * m + 15, 101)]
    rho_ok = all(abs(r - 4 * math.sin(math.pi / n) ** 2) < 1e-12 for r, n in zip(rho, (8 * m + 15, 101)))
    ineqs = []
    for tag, seg, lo, lam in (('J0', J0, lo0, lam0), ('J1', J1, hi0, lam1)):
        Ss, _ = walk(seg)
        for k, (o, e, s, s2) in enumerate(seg):
            g = (o, e); a, b = lo - Ss[k], lo - Ss[k] + lam
            for ss in (s, s2):
                p, q = side_range(g, ss)
                if p.evx(math.pi / N0).real > q.evx(math.pi / N0).real: p, q = q, p
                ineqs += [a - p, q - b]
            for ss, want in ((s, 1), (s2, -1)):
                i, kk = SIDES[ss]; ev = g_apply(g, V[kk]) - g_apply(g, V[i])
                ineqs.append((ev + ev.conj()) * ((1 if e == 0 else -1) * want))
    ineqs += [lam0, lam1]
    strict = [f for f in ineqs if not f.iszero()]
    KP.XMAX = Fr(1, 10 ** 6) * int(math.pi / (8 * m + 15) * 10 ** 6 + 1)
    try:
        cert = [KP.certify_nonneg(f) for f in strict]; cok = True; mn = min(c[1] for c in cert)
    except AssertionError as ex:
        cok, mn = False, str(ex)
    # control: first N below the range should break
    Nc = 8 * m + 13
    viol = sum(1 for f in strict if f.evx(math.pi / Nc).real <= 0)
    res[m] = dict(h=t['h'], start=(g0, s0), comb=ok, returns=ret, rho_4sin2=rho_ok, n_ineq=len(ineqs), tight=len(ineqs) - len(strict),
                  certified=cok, min_lb=mn, XMAX=float(KP.XMAX), control_N=Nc, control_viol=viol,
                  lam0=lam0.c and {k: [v.real, v.imag] for k, v in lam0.c.items()}, lam1={k: [v.real, v.imag] for k, v in lam1.c.items()})
    print(m, {k: v for k, v in res[m].items() if k not in ('lam0', 'lam1')}, flush=True)
json.dump({m: {k: (str(v) if k in ('start',) else v) for k, v in r.items()} for m, r in res.items()},
          open('data/s935_ladder_perm_cert.json', 'w'), indent=0, default=str)
