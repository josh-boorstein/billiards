"""[Deposit] Appendix B, the case m <= 6 (and the perturbed-word control). Paper: "Complete periodicity of perpendicular billiards in rational triangles of prime denominator". Research name: probes/s948_q3_small.py, copied verbatim except imports (deposit names), paths (data/ and figs/ under triangles/) and this paragraph. Run through run_all.py; the usage lines below are the research repo's. Bracketed IDs and section marks (§sNNN) refer to the project's internal ledgers.

s948 `s947-a` part (2): p = m AT SMALL m (the word families start at m = 7) — for F1 = (3, 6m+2, 12m+4), F2 = (3, 6m+1, 12m+5),
m = 1 … MMAX, EVERY perpendicular-class slot whose walk ends perpendicularly is CERTIFIED: its walk word (constant runs) is realised
— every crossing condition of `engine/word_family.py` holds, by mpmath interval arithmetic at the point — and its perpendicular
ending is an EXACT identity in ℚ(ζ_{4N}) (`cos_cert.is_zero`). With the 2m vertex-ending slots (counted) this gives p = m by the
counting argument of `probes/s948_q3_legality.py` ([TRI-BAND-QM] supplies m distinct swapped prongs).
Usage: `probes/s948_q3_small.py [MMAX]` → `data/s948_q3_small.json`. PRIOR ART / DATA: as `probes/s948_q3_legality.py`."""
import json, sys
import _setup  # noqa: F401  (deposit: lib/ on the path, cwd = triangles/)
from fractions import Fraction as Fr
import mpmath as mp
from mpmath import iv
import cos_cert as cc
import word_family as wf
import appB_1_families as QF

mp.mp.dps = 40
iv.dps = 40
FAMS = {'F1': [(0, 3), (6, 2), (12, 4)], 'F2': [(0, 3), (6, 1), (12, 5)]}


def at_point(f, m):
    """formal sum at fixed m (t = 0 since runs/direction are constant): [(c, n/(2N))] with cos(π·n/(2N))"""
    N = 18 * m + 9
    out = {}
    for (a, b, g), c in f.items():
        r = (Fr(a) * m + Fr(g)) / (2 * N)
        out[r] = out.get(r, 0) + c
    return [(c, r) for r, c in out.items() if c]


def iv_sum(pairs):
    return sum((c * iv.cos(iv.pi * iv.mpf(r.numerator) / r.denominator) for c, r in pairs), iv.mpf(0))


def certify_word(fam, m, v, d, runs, firsts):
    rr = [(wf.E(0, 0, r), 0 if f == '-' else int(f)) for r, f in zip(runs, firsts)]
    w = wf.Word(fam, v, wf.E(0, 0, d), rr, 'perp')
    conds, spans = w.conditions((m % 2, 0))
    ok = all(3 * (int(wf.val(r, m, 0)) - 1) < 18 * m + 9 for r in spans)
    for cd in conds:
        name, kd, f, c = cd[:4]
        if len(cd) > 4 and wf.val(cd[4][1], m, 0) < 2:
            continue
        if kd == 'zero':
            ok &= cc.is_zero(at_point(f, m), real=True)
        elif kd in ('opp', 'same'):
            pr = iv_sum(at_point(f, m)) * iv_sum(at_point(c, m))
            ok &= bool(pr.b < 0) if kd == 'opp' else bool(pr.a > 0)
    return ok


def control(famname, m):
    """CONTROL: every perp slot's word with one run length changed by ±1 (each run, each sign, keeping lengths ≥ 0) must FAIL"""
    fam = wf.Fam(FAMS[famname])
    tried = passed = 0
    for row in QF.walks(famname, m):
        (v, d), kind, _, runs, firsts = row
        if kind != 'perp':
            continue
        for j in range(len(runs)):
            for e in (-1, 1):
                r2 = list(runs); r2[j] += e
                if r2[j] < 0 or (r2[j] == 0) != (runs[j] == 0):
                    continue
                tried += 1
                passed += certify_word(fam, m, v, d, r2, firsts)
    return tried, passed


def check(famname, m):
    fam = wf.Fam(FAMS[famname])
    rows = QF.walks(famname, m)
    st = {'slots': len(rows), 'perp': 0, 'vertex': 0, 'certified': 0, 'fail': []}
    for row in rows:
        (v, d), kind, _, runs, firsts = row
        if kind == 'vertex':
            st['vertex'] += 1
            continue
        assert kind == 'perp', kind
        st['perp'] += 1
        rr = [(wf.E(0, 0, r), 0 if f == '-' else int(f)) for r, f in zip(runs, firsts)]
        w = wf.Word(fam, v, wf.E(0, 0, d), rr, 'perp')          # at fixed m every exponent is evaluated at m: d is the slot index
        conds, spans = w.conditions((m % 2, 0))
        ok = all(3 * (int(wf.val(r, m, 0)) - 1) < 18 * m + 9 for r in spans)
        for cd in conds:
            name, kd, f, c = cd[:4]
            if len(cd) > 4 and wf.val(cd[4][1], m, 0) < 2:
                continue
            if kd == 'zero':
                ok &= cc.is_zero(at_point(f, m), real=True)
            elif kd in ('opp', 'same'):
                pr = iv_sum(at_point(f, m)) * iv_sum(at_point(c, m))
                ok &= bool(pr.b < 0) if kd == 'opp' else bool(pr.a > 0)
        if ok:
            st['certified'] += 1
        else:
            st['fail'].append([v, d])
    return st


def main(MMAX):
    out, allok = {}, True
    for famname in FAMS:
        for m in range(1, MMAX + 1):
            st = check(famname, m)
            ok = st['vertex'] == 2 * m and st['certified'] == st['perp'] and not st['fail']
            allok &= ok
            out[f'{famname} m={m}'] = dict(st, ok=ok)
            print(famname, m, {k: v for k, v in st.items() if k != 'fail'}, 'fails', st['fail'][:5], 'OK' if ok else 'FAIL', flush=True)
    for famname in FAMS:
        for m in (3, 5):
            tried, passed = control(famname, m)
            out[f'control {famname} m={m}'] = {'perturbed_words': tried, 'wrongly_certified': passed}
            print('CONTROL', famname, m, 'perturbed words', tried, 'wrongly certified', passed, flush=True)
            allok &= passed == 0 and tried > 0
    json.dump(out, open('data/s948_q3_small.json', 'w'), indent=1)
    assert allok
    print('ALL OK')


if __name__ == '__main__':
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 6)
