"""[Deposit] Appendix B, the case m >= 100 (and the extension control). Paper: "Complete periodicity of perpendicular billiards in rational triangles of prime denominator". Research name: probes/s948_q3_legality.py, copied verbatim except imports (deposit names), paths (data/ and figs/ under triangles/) and this paragraph. Run through run_all.py; the usage lines below are the research repo's. Bracketed IDs and section marks (§sNNN) refer to the project's internal ledgers.

s948 `s947-a` part (1): THE q = 3 FAMILIES' PERPENDICULAR WORD FAMILIES ARE LEGAL FOR EVERY m ≥ M0 — interval certificates of
every crossing condition on each family's (x, y) region, x = π/2N → 0 (`engine/cos_cert.py` Taylor charts).

SETTING. `probes/s947_q3_families.py`: F1 = (3, 6m+2, 12m+4), F2 = (3, 6m+1, 12m+5), N = 18m + 9; each strand (v, res) of slots
d = idx_v + 6t + 1 + 2res is tiled by word families with t ∈ [lo(m), hi(m)], lo, hi affine in m except 4 TRANSCENDENTAL boundaries
per family ([TRI-Q3-GRAZE]). `engine/word_family.py` turns a family into crossing conditions (per parity class of (m, t)):
sign conditions 'opp'/'same' (a vertex image f against the star centre c) and the perpendicular ending (an exact identity).
REGION of a family for m ≥ M0: {0 < x ≤ X0, lo ≤ y ≤ hi}, y = t·x, an affine bound t = am + b being the line y = aπ/36 + (b − a/2)x.
A 'needs r ≥ 2' condition is required only on its half-plane. At a transcendental boundary i|i+1 with boundary function f (the
lower family's earliest failing condition at hi + 1; every function of families i, i+1 EXACTLY ± f is a band MEMBER) the zero set of f in the joint region is
pinned in a band Y1 < y < Y2 (f certified of the lower sign on y ≤ Y1, the upper on y ≥ Y2), every OTHER condition of family i is
certified up to Y2 and of family i+1 from Y1; inside the band f ≠ 0 at every lattice point ([TRI-Q3-GRAZE]) picks the family.
Vertex families (F1 16, 17; F2 20, 21) are not certified: they hold exactly 2m slots (tiling), the m band prongs of [TRI-BAND-QM]
(t = 2/3, 1/3 < t_up) and their reverses, so they are those prongs by counting.
PRIOR ART: §s947 (1)–(3) and its "Not done" list (this is item (1)); `probes/s938_q5_loops.py` (the hand charts this generalises);
`rulings.py --grep` *legality*, *interval certificate*, *word family* → nothing on this object. DATA: `data/s947_q3_*.json` (s947).
Usage: `probes/s948_q3_legality.py F1|F2 [M0]` → `data/s948_q3_legality_{fam}.json`."""
import json, sys, time
import _setup  # noqa: F401  (deposit: lib/ on the path, cwd = triangles/)
from fractions import Fraction as Fr
import mpmath as mp
import cos_cert as cc
import word_family as wf
import appB_2_boundaries as G

mp.mp.dps = 30
NM, N0 = 18, 9
PARS = ((0, 0), (0, 1), (1, 0), (1, 1))


def lkey(L):
    return (mp.nstr(L.Y0, 40), str(L.s), L.L0)


def aff(fit):
    """t_lo / t_hi fit [a, b, q] → (a, b) when q = 1 (affine), else None (transcendental)"""
    a, b, q = (int(x) for x in fit)
    return (Fr(a), Fr(b)) if q == 1 else None


def line_of(ab):
    a, b = ab
    return cc.Line.exact(a / 36, b - a / 2)


def hp_of(r, par):
    """r = (rm, rt, r0) ≥ 2 on the parity class `par` as a half-plane πA + b·y + c·x ≥ 0 (m x = π/36 − x/2): r is an integer of
    FIXED parity on the class, so r ≥ 2 ⟺ r ≥ 3 when r is odd there (the lattice-tight threshold)"""
    rm, rt, r0 = (Fr(v) for v in r)
    thr = 2 + int(rm * par[0] + rt * par[1] + r0) % 2
    return (rm / 36, rt, r0 - thr - rm / 2)


def subst_t(f, ab):
    """a single-t family: substitute t = a·m + b in a formal sum"""
    a, b = ab
    out = {}
    for (al, be, ga), c in f.items():
        e = (Fr(al) + Fr(be) * a, Fr(0), Fr(ga) + Fr(be) * b)
        out[e] = out.get(e, 0) + c
    return {e: c for e, c in out.items() if c}


def cond_funcs(w, par):
    """[(name, kind, f, c, r_or_None)] for the sign conditions and the 'zero' identities"""
    conds, spans = w.conditions(par)
    return [(cd[0], cd[1], cd[2], cd[3], cd[4][1] if len(cd) > 4 else None) for cd in conds], spans


def setup(famname):
    fam, cls, W = G.words_for(famname)
    graze = json.load(open('data/s947_q3_graze.json'))
    trans = {}                                            # i (lower family) → boundary condition names (lower, upper)
    for key, v in graze.items():
        if key.startswith(famname + ' '):
            i = int(key.split()[1].split('|')[0])
            trans[i] = v['condition']
    return fam, W, trans


def run(famname, M0=100, eta=mp.mpf('0.002')):
    global fam_name_global
    fam_name_global = famname
    fam, W, trans = setup(famname)
    X0 = mp.pi / (2 * (NM * M0 + N0)) + mp.mpf(10) ** -25          # OUTWARD (audit G17)
    t0 = time.time()
    fams = []
    for T, w in W:
        fams.append({'i': T['i'], 'kind': T['kind'], 'w': w, 'lo': aff(T['t_lo']), 'hi': aff(T['t_hi'])})
    import appB_1_families as QF
    walks = json.load(open(f'data/s947_q3_walks_{famname}.json'))
    seg30 = QF.segments(walks['30']['rows'])
    x30 = mp.pi / (2 * (NM * 30 + N0))
    # transcendental bands: the lower family's boundary function f_i (per parity), its x = 0 zero y_c
    bands = {}
    for i, nm in trans.items():
        fi = {par: [cd[2] for cd in W[i][1].conditions(par)[0] if cd[0] == nm][0] for par in PARS}
        T0 = cc.from_formal(fi[(0, 0)], NM, N0)
        k = cc.xorder(T0)
        Ck = lambda y, T0=T0, k=k: sum(c * (mp.mpf(de.numerator) / de.denominator) ** k * mp.cos(
            mp.pi * mp.mpf(th.numerator) / th.denominator + be * y + k * mp.pi / 2) for c, th, be, de in T0)
        thi = seg30[i][2][-1][0]
        yc = mp.findroot(Ck, (thi + mp.mpf(1) / 2) * x30)
        bands[i] = {'name': nm, 'f': fi, 'k': k, 'yc': yc, 'Y1': cc.Line(yc - eta, 0), 'Y2': cc.Line(yc + eta, 0),
                    'members': {i: [], i + 1: []}}
        print(f'band {i}|{i+1}: f = {nm}  k={k}  y_c = {mp.nstr(yc, 12)} (= {mp.nstr(yc * 36 / mp.pi, 12)}·π/36)', flush=True)

    def band_eps(i, par, g):
        """ε if g ≡ ε·f_i (exact), else 0"""
        f = bands[i]['f'][par]
        if wf.identity_zero(wf.add(f, g, -1), fam)[0]:
            return 1
        if wf.identity_zero(wf.add(f, g, 1), fam)[0]:
            return -1
        return 0

    def t_ref(F, mm):                                     # an interior lattice t of family F at m = mm (exact for affine bounds)
        xm = mp.pi / (2 * (NM * mm + N0))
        lo = int(F['lo'][0] * mm + F['lo'][1]) if F['lo'] is not None else int(mp.ceil((bands[F['i'] - 1]['yc'] + eta) / xm))
        hi = int(F['hi'][0] * mm + F['hi'][1]) if F['hi'] is not None else int(mp.floor((bands[F['i']]['yc'] - eta) / xm))
        assert lo <= hi, (F['i'], lo, hi)
        return (lo + hi) // 2, xm
    cache, results, allok = {}, [], True
    for F in fams:
        if F['kind'] != 'perp':
            continue
        i = F['i']
        lo_o = line_of(F['lo']) if F['lo'] is not None else bands[i - 1]['Y1']
        lo_b = None if F['lo'] is not None else bands[i - 1]['Y2']
        hi_o = line_of(F['hi']) if F['hi'] is not None else bands[i]['Y2']
        hi_b = None if F['hi'] is not None else bands[i]['Y1']
        tm, xm = t_ref(F, 200)
        for par in PARS:
            if F['lo'] is not None and F['lo'] == F['hi'] and int(F['lo'][0] * par[0] + F['lo'][1]) % 2 != par[1]:
                continue                                          # single-t family: t's parity is fixed by m's
            conds, spans = cond_funcs(F['w'], par)
            for name, kind, f, c, r in conds:
                if kind == 'zero':
                    single = F['lo'] is not None and F['lo'] == F['hi']
                    ok = wf.identity_zero(subst_t(f, F['lo']) if single else f, fam)[0]
                    results.append({'fam': i, 'par': par, 'cond': name, 'kind': 'zero', 'ok': ok})
                    allok &= ok
                    continue
                if kind not in ('opp', 'same'):
                    continue
                if r is not None and Fr(r[0]) == 0 and Fr(r[1]) == 0 and Fr(r[2]) < 2:
                    continue                                          # r < 2 constant: the condition is never required
                hps = (hp_of(r, par),) if r is not None else ()
                Tc, Tf = cc.from_formal(c, NM, N0), cc.from_formal(f, NM, N0)
                sc = 1 if cc.evalf(Tc, xm, tm * xm) > 0 else -1
                sf = -sc if kind == 'opp' else sc
                for lab, g, Tg, sg in (('c', c, Tc, sc), ('f', f, Tf, sf)):
                    lo_r, hi_r, memb = lo_o, hi_o, None
                    if hi_b is not None:
                        e = band_eps(i, par, g)
                        if e:
                            hi_r, memb = hi_b, ('lower', i, e)
                            bands[i]['members'][i].append((par, name, lab, sg, e))
                    if lo_b is not None:
                        e = band_eps(i - 1, par, g)
                        if e:
                            lo_r, memb = lo_b, ('upper', i - 1, e)
                            bands[i - 1]['members'][i].append((par, name, lab, sg, e))
                    key = (tuple(Tg), sg, lkey(lo_r), lkey(hi_r), hps)        # VALUES, not id() (audit F20: ids are reused)
                    if key not in cache:
                        cache[key] = cc.certify_region(Tg, sg, lo_r, hi_r, X0, hps)
                    ok, inf = cache[key]
                    allok &= ok
                    results.append({'fam': i, 'par': par, 'cond': name, 'func': lab, 'sign': sg, 'band': memb, 'ok': ok,
                                    'info': json.loads(json.dumps(inf, default=str))})
                    if not ok:
                        print(f'  FAIL fam {i} par {par} {name} [{lab}] {json.dumps(inf, default=str)[:300]}', flush=True)
        print(f'fam {i}: done ({len(cache)} distinct certificates, {time.time() - t0:.0f}s)', flush=True)
    # band consistency: at a band lattice point f_i ≠ 0 ([TRI-Q3-GRAZE]); σ := the sign family i requires of f_i. Every member
    # g = ε f_i of family i must require σ_g ε = σ, of family i+1 σ_g ε = −σ — then sign(f_i) selects a family all of whose band
    # members hold, its other functions being certified across the band.
    for i, B in bands.items():
        sig = {}
        for par, name, lab, sg, e in B['members'][i]:
            if name == B['name'] and lab == 'f':
                sig[par] = sg * e
        okb = len(sig) == 4 and all(sg * e == sig[par] for par, name, lab, sg, e in B['members'][i]) and \
            all(sg * e == -sig[par] for par, name, lab, sg, e in B['members'][i + 1]) and len(B['members'][i + 1]) > 0
        B['sign_consistent'] = okb
        B['sig_lo'] = sig
        allok &= okb
        print(f'band {i}|{i+1}: members lower {len(B["members"][i])}, upper {len(B["members"][i + 1])}, sign consistency', okb,
              flush=True)
    return {'fam': famname, 'M0': M0, 'X0': str(X0), 'eta': str(eta), 'all_ok': allok,
            'bands': {str(k): {'name': v['name'], 'k': v['k'], 'yc': str(v['yc']), 'sign_consistent': v['sign_consistent'],
                               'sig_lo': {str(p): s for p, s in v['sig_lo'].items()},
                               'members': {str(j): v['members'][j] for j in v['members']}} for k, v in bands.items()},
            'results': results}


def control(famname, M0=100):
    """CONTROL: each perp family whose affine upper bound is NOT its strand's end, with that bound moved one slot up
    (t ≤ hi + 1), must lose at least one certificate (the next slot belongs to another word)"""
    fam, W, trans = setup(famname)
    X0 = mp.pi / (2 * (NM * M0 + N0)) + mp.mpf(10) ** -25          # OUTWARD (audit G17)
    fams = [{'i': T['i'], 'kind': T['kind'], 'w': w, 'lo': aff(T['t_lo']), 'hi': aff(T['t_hi']), 'v': T['vertex'],
             'res': T['strand']} for T, w in W]
    out = {}
    for F, G2 in zip(fams, fams[1:] + [None]):
        if F['kind'] != 'perp' or F['lo'] is None or F['hi'] is None or G2 is None or (G2['v'], G2['res']) != (F['v'], F['res']):
            continue
        hi2 = (F['hi'][0], F['hi'][1] + 1)
        lo_r, hi_r = line_of(F['lo']), line_of(hi2)
        fired = False
        for par in PARS:
            conds, _ = cond_funcs(F['w'], par)
            for name, kind, f, c, r in conds:
                if kind not in ('opp', 'same'):
                    continue
                if r is not None and Fr(r[0]) == 0 and Fr(r[1]) == 0 and Fr(r[2]) < 2:
                    continue
                hps = (hp_of(r, par),) if r is not None else ()
                mm = 200
                xm = mp.pi / (2 * (NM * mm + N0))
                tm = int(F['lo'][0] * mm + F['lo'][1])
                Tc, Tf = cc.from_formal(c, NM, N0), cc.from_formal(f, NM, N0)
                sc = 1 if cc.evalf(Tc, xm, tm * xm) > 0 else -1
                sf = -sc if kind == 'opp' else sc
                for Tg, sg in ((Tc, sc), (Tf, sf)):
                    if not cc.certify_region(Tg, sg, lo_r, hi_r, X0, hps, maxn=20000)[0]:
                        fired = True
                        break
                if fired:
                    break
            if fired:
                break
        out[F['i']] = fired
        print(f'control fam {F["i"]}: extended to hi+1 — certificate lost: {fired}', flush=True)
    return out


if __name__ == '__main__':
    if sys.argv[1] == 'control':
        res = {f: control(f) for f in ('F1', 'F2')}
        json.dump(res, open('data/s948_q3_legality_control.json', 'w'), indent=1)
        print('CONTROL fired on', sum(v for r in res.values() for v in r.values()), 'of', sum(len(r) for r in res.values()))
        sys.exit(0)
    famname = sys.argv[1]
    M0 = int(sys.argv[2]) if len(sys.argv) > 2 else 100
    out = run(famname, M0)
    json.dump(out, open(f'data/s948_q3_legality_{famname}.json', 'w'), indent=1, default=str)
    print('ALL OK' if out['all_ok'] else 'FAILURES', flush=True)
