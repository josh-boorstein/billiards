"""Run every script of this directory (paper 1: perpendicular billiards at prime denominator) in order, and check each
one's verdict line and recorded output against the figures the paper prints; exit nonzero if anything fails.

  python run_all.py            about fifteen minutes on one core
  python run_all.py --quick    skip step 1 and the steps that read its output (about two minutes)

The scripts are the research keepers, copied (see README.md); each writes its record to data/ and the figures to figs/.
"""
import json
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(HERE, 'data')


def load(name):
    with open(os.path.join(D, name)) as f:
        return json.load(f)


def chk_families(fam, n_perp):
    d = load(f's947_q3_families_{fam}.json')
    return [(f'{fam}: {n_perp} word families ending at a right angle (+ 2 band-prong families)', len(d['families']) == n_perp + 2),
            (f'{fam}: no word transitions in m on the held-out m = 20..30', d['transitions'] == []),
            (f'{fam}: exactly 4 transcendental boundaries', len(d['holdout_fail']) == 8)]


def chk_legality(fam, n_sign):
    d = load(f's948_q3_legality_{fam}.json')
    sign = sum(1 for r in d['results'] if r.get('kind') is None)
    return [(f'{fam}, m >= 100: all {n_sign} sign conditions and every perpendicular ending certified',
             d['all_ok'] and d['M0'] == 100 and sign == n_sign and all(r['ok'] for r in d['results']))]


def chk_control():
    d = load('s948_q3_legality_control.json')
    vals = [v for r in d.values() for v in r.values()]
    return [('control: each of the 19 families extended one slot loses its certificate', len(vals) == 19 and all(vals))]


def chk_lattice(fam, n_cond):
    d = load(f's948_q3_lattice_{fam}.json')
    L = d['lattice']
    return [(f'{fam}, tiling every m >= 7 and runs/spans', d['tiling_problems'] == [] and d['run_span_problems'] == []),
            (f'{fam}, 7 <= m <= 100: {L["slots"]} slots, {n_cond} conditions, no failure, no interval fallback',
             d['all_ok'] and L['slots'] == 81310 and L['conds'] == n_cond and L['n_fail'] == 0 and L['fallback'] == 0)]


def chk_small():
    d = load('s948_q3_small.json')
    rows = {k: v for k, v in d.items() if not k.startswith('control')}
    ctl = [v for k, v in d.items() if k.startswith('control')]
    ok_rows = all(v['ok'] and v['vertex'] == 2 * int(k.split('m=')[1]) and v['certified'] == v['perp'] for k, v in rows.items())
    return [('m <= 6: every perpendicular slot certified, exactly 2m vertex slots (both triangles)', len(rows) == 12 and ok_rows),
            ('control: 1159 perturbed words at m = 3, 5, none certified',
             sum(v['perturbed_words'] for v in ctl) == 1159 and all(v['wrongly_certified'] == 0 for v in ctl))]


J_PAPER = {(9, 576): (2.02159990, 2.02191910), (15, 712): (4.01917575, 4.01943080),
           (15, 756): (4.02988894, 4.03003593), (27, 1408): (8.02604251, 8.02619503)}


def chk_thmD():
    d = load('s946_family_cert.json')
    out = []
    got = {(r['n'], r['L']): r for r in d}
    out.append(('Theorem D: the four families (n, |c|) of the paper', set(got) == set(J_PAPER)))
    for k, r in sorted(got.items()):
        ok = (r['linear_identity'] and r['translation_all_vertices'] and not r['reversible'] and r['T_nonzero']
              and r['ball_cells'] > 0 and r['ball_worst_margin'] > 0 and max(r['width_just_outside']) < 0
              and abs(r['J'][0] - J_PAPER.get(k, (0, 0))[0]) < 1e-8
              and abs(r['J'][1] - J_PAPER.get(k, (0, 0))[1]) < 1e-8 and len(r['core']) == k[1])
        out.append((f'n = {k[0]}, |c| = {k[1]}: (i)-(iv) exact, J = [{r["J"][0]:.8f}, {r["J"][1]:.8f}] certified in ball '
                    f'arithmetic, W < 0 just outside', ok))
    pts = {}
    for r in d:
        for abc, t, inside in r['exact_points']:
            pts[tuple(abc)] = pts.get(tuple(abc), False) or inside
    out.append((f'the {len(pts)} exact free-cylinder points (flow decompositions) each lie in a certified J', len(pts) == 13 and all(pts.values())))
    return out


# (script, args, verdict substring or None (the script asserts), post-check or None, needs step 1)
STEPS = [
    ('s6_1_band_identities.py', [], 'ALL OK', None, False),
    ('s6_2_length_law.py', [], 'ALL OK', None, False),
    ('appB_1_families.py', ['F1'], None, lambda: chk_families('F1', 28), True),
    ('appB_1_families.py', ['F2'], None, lambda: chk_families('F2', 24), True),
    ('appB_2_boundaries.py', [], 'ALL OK', None, True),
    ('appB_3_legality.py', ['F1'], 'ALL OK', lambda: chk_legality('F1', 2084), True),
    ('appB_3_legality.py', ['F2'], 'ALL OK', lambda: chk_legality('F2', 1936), True),
    ('appB_3_legality.py', ['control'], 'CONTROL fired on 19 of 19', chk_control, True),
    ('appB_4_lattice.py', ['F1'], 'ALL OK', lambda: chk_lattice('F1', 787290), True),
    ('appB_4_lattice.py', ['F2'], 'ALL OK', lambda: chk_lattice('F2', 737279), True),
    ('appB_5_small.py', [], 'ALL OK', chk_small, False),
    ('appC_certificate.py', [], None, chk_thmD, False),
    ('figures.py', [], 'ALL OK', None, False),
]


def main():
    quick = '--quick' in sys.argv
    bad = []
    for script, args, verdict, post, needs1 in STEPS:
        if quick and needs1:
            continue
        name = ' '.join([script] + args)
        print(f'=== {name}', flush=True)
        t0 = time.time()
        p = subprocess.run([sys.executable, os.path.join(HERE, script)] + args, cwd=HERE, capture_output=True, text=True)
        tail = (p.stdout.strip().split('\n') or [''])[-1]
        ok = p.returncode == 0 and (verdict is None or verdict in p.stdout)
        print(f'  [{"ok" if ok else "FAIL"}] exit {p.returncode}; last line: {tail[:160]}  ({time.time() - t0:.0f} s)')
        if not ok:
            print(p.stderr[-2000:])
        checks = post() if (ok and post) else []
        for label, c in checks:
            print(f'  [{"ok" if c else "FAIL"}] {label}')
            ok &= bool(c)
        if not ok:
            bad.append(name)
    n = sum(1 for s in STEPS if not (quick and s[4]))
    print('\nFAILED: ' + ', '.join(bad) if bad else f'\nall {n} steps passed')
    sys.exit(1 if bad else 0)


if __name__ == '__main__':
    main()
