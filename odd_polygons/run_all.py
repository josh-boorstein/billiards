"""Run every script of this directory (paper 2: perpendicular billiards in rational polygons of odd denominator) in order, then
check each one's record against the figures the paper prints; exit nonzero if anything fails.

  python run_all.py              about two and a half hours on one core
  python run_all.py --quick      skip the long certificates, the catalogue and the steps that read their output (about 15 minutes;
                                 the catalogue is then the shipped research record)
  python run_all.py --check      only the checks, on the records already in data/

The scripts are the research keepers, copied (see README.md); each writes its record to data/ and the figures to figs/.
The catalogue of non-periodic triangle classes is recomputed from scratch (28 triangles, six directions each, about 20 minutes)
and compared with the research record, shipped as data/s626_noncp_catalogue_shipped.json.
"""
import json
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.environ.get('RUN_ALL_DATA', os.path.join(HERE, 'data'))
CAT = os.path.join(D, 's626_noncp_catalogue.json')
SHIPPED = os.path.join(D, 's626_noncp_catalogue_shipped.json')   # the research record, for comparison

# (script, arguments, about how many minutes, long?)   long steps and their readers are skipped by --quick
STEPS = [
    ('area_test.py', [], 3, False),
    ('s4_5_kite_census.py', [], 4, False),
    ('s4_4_rhombus.py', [], 1, False),
    ('s4_7_scan.py', [], 35, True),
    ('s4_7_witnesses.py', [], 3, True),
    ('s4_11_counts.py', [], 1, False),
    ('noncp_catalogue.py', [], 20, True),
    ('s5_2_kites.py', ['--trace'], 15, True),
    ('s5_3_fold.py', [], 2, False),
    ('pi_n_torus.py', [], 1, False),
    ('thmC_identities.py', [], 2, False),
    ('thmC_cert_s951.py', [], 11, True),
    ('thmC_m1.py', ['AB'], 1, False),
    ('thmC_cert.py', [], 15, True),
    ('thmC_family.py', [], 2, True),
    ('s5_6_onset.py', [], 4, False),
    ('s6_9_census.py', ['15'], 3, False),
    ('s6_9_five.py', [], 1, False),
] + [('s6_10_walls.py', [t], 2, False) for t in
     ['5:1,2,3,4', '1256', '1346', '9:1,2,7,8', '9:1,3,6,8', '9:1,4,5,8', '11:1,2,9,10', '11:1,3,8,10', '11:1,4,7,10',
      '11:1,5,6,10', '13:1,2,11,12', '13:1,3,10,12', '13:1,4,9,12', '13:1,5,8,12', '13:1,6,7,12',
      '5:1,1,4,4', '7:1,1,6,6', '9:1,1,8,8', '11:1,1,10,10']] + [
    ('s6_10_isosceles.py', ['7:1,1,6,6', '5:1,1,4,4', '9:1,1,8,8'], 3, False),
    ('appA_1_identities.py', [], 1, False),
    ('appA_1_per_m.py', [], 21, True),
    ('appA_1_all_m.py', [], 17, True),
    ('appA_2_identity.py', [], 1, False),
    ('appA_2_legality.py', [], 1, False),
    ('appA_3_frontier.py', [], 1, False),
    ('appA_chain_law.py', [str(n) for n in range(15, 64, 2)], 5, False),
    ('appA_corners.py', [], 3, False),
    ('figures.py', [], 1, False),
]


def prepare_catalogue():
    """Start the catalogue empty, so noncp_catalogue.py recomputes every row; the research record stays as SHIPPED."""
    if os.path.exists(CAT):
        os.remove(CAT)


def run_steps(quick):
    log = []
    if quick and not os.path.exists(CAT):            # --quick skips the catalogue step: read the research record instead
        import shutil
        shutil.copy(SHIPPED, CAT)
    for script, args, minutes, long_ in STEPS:
        if quick and long_:
            print(f'--- skip {script} {" ".join(args)} (long)', flush=True)
            continue
        if script == 'noncp_catalogue.py':
            prepare_catalogue()
        t = time.time()
        print(f'--- {script} {" ".join(args)} (about {minutes} min)', flush=True)
        os.makedirs(os.path.join(D, 'logs'), exist_ok=True)
        tag = '_'.join([script[:-3]] + [a.replace(':', '_').replace(',', '') for a in args[:1]])
        with open(os.path.join(D, 'logs', tag + '.txt'), 'w') as lf:
            r = subprocess.run([sys.executable, script] + args, cwd=HERE, stdout=lf, stderr=subprocess.STDOUT)
        log.append({'script': script, 'args': args, 'rc': r.returncode, 'secs': round(time.time() - t)})
        print(f'    exit {r.returncode}, {log[-1]["secs"]} s', flush=True)
        json.dump(log, open(os.path.join(D, 'run_all_log.json'), 'w'), indent=1)
    return log


def load(name):
    with open(os.path.join(D, name)) as f:
        return json.load(f)


def log_text(tag):
    p = os.path.join(D, 'logs', tag + '.txt')
    return open(p).read() if os.path.exists(p) else ''


def f0(x):
    return float(x[0] if isinstance(x, list) else x)


def chk_kites():
    k = [r for r in load('s928_kites.json') if r['label'].startswith('kite') and r['N'] not in r['a']]   # no flat vertex
    p1, p3 = [r for r in k if r['p'] == 1], [r for r in k if r['p'] >= 3]
    free = [r for r in p3 if f0(r['one_minus_R']) > 1e-9]
    best, bp = max(free, key=lambda r: f0(r['one_minus_R'])), max(k, key=lambda r: r['p'])
    return [('Remark 4.5: 134 kites, 89 with p = 1, 45 with p >= 3, every p odd', (len(k), len(p1), len(p3)) == (134, 89, 45)
             and all(r['p'] % 2 for r in k)),
            ('Remark 4.5: 44 free (42 convex, 2 darts); the one p >= 3 kite without, (9,2,17,2)/15',
             len(free) == 44 and sum(max(r['a']) < r['N'] for r in free) == 42
             and [(r['a'], r['N']) for r in p3 if r not in free] == [([9, 2, 17, 2], 15)]),
            ('Remark 4.5: every prong of all 134 reaches a vertex', all(r['sprime'] == 0 for r in k)),
            ('Remark 4.5: largest free fraction 0.735 at (7,6,7,6)/13 (p = 9); largest p 11 at (7,8,7,8)/15 (0.675)',
             (best['a'], best['N'], best['p'], round(f0(best['one_minus_R']), 3)) == ([7, 6, 7, 6], 13, 9, 0.735)
             and (bp['a'], bp['N'], bp['p'], round(f0(bp['one_minus_R']), 3)) == ([7, 8, 7, 8], 15, 11, 0.675))]


def chk_rhombus():
    t = log_text('s4_4_rhombus')
    return [('Proposition 4.4: the rhombus (3,2,3,2)/5 — the script ran and printed the exact free fraction phi/5',
             'φ/5' in t or 'phi/5' in t)]


def chk_asym():
    w = load('s929_asym_witness.json')
    ok = [round(float(r['t_star']), 12) for r in w] == [0.095459452979, 0.540706695280] and all(f0(r['one_minus_R']) > 0.01 for r in w)
    mc = all(abs(f0(r['one_minus_R']) - r['mc']['share']) < 3 * r['mc']['sigma'] for r in w)
    fam = {}
    for r in load('s929_asym_free.json'):
        if r['kind'] == 'tstar':
            fam.setdefault(r['label'].split('t=')[0].replace('tstar', '').strip(), []).append(f0(r['one_minus_R']))
    free = sorted(f for f, v in fam.items() if any(x > 1e-6 for x in v))
    return [('Example 4.7: the witnesses (2,4,5,7)/9 at t* = 0.0954594529790…, (4,5,9,12)/15 at t* = 0.5407066952799…, both free', ok),
            ('Example 4.7: the area agrees with the independent Monte Carlo estimate (within 3 sigma)', mc),
            ('Example 4.7: of five families scanned, two show free regions', len(fam) == 5 and free == ['(2, 4, 5, 7)/9', '(4, 5, 9, 12)/15'])]


def chk_counts():
    v = load('s958_free_counts.json')['verdict']
    return [(f'Remark 4.11: {k}', ok) for k, ok in v.items()]


def chk_catalogue():
    old, new = load('s626_noncp_catalogue_shipped.json'), load('s626_noncp_catalogue.json')
    same = sorted(old) == sorted(new) and all(
        [(x['dir'], x['counts']) for x in new[k]['dirs']] == [(x['dir'], x['counts']) for x in old[k]['dirs']] for k in old)
    return [('Example 5.2, Proposition 5.3(b): the 28 catalogued triangles re-certified in all six directions; every '
             'decomposition equals the research record', same)]


def chk_kites_noncp():
    s = load('s929_kite_noncp.json')['selection']
    got = [(x['kite'], x['N'], x['counts'].get('no_periodic')) for x in s]
    return [('Example 5.2: the kites (41,15,67,15)/69, (201,23,23,23)/135, (119,45,45,45)/127, each with components without '
             'periodic orbit', got == [([41, 15, 67, 15], 69, 4), ([201, 23, 23, 23], 135, 4), ([119, 45, 45, 45], 127, 4)]),
            ('Appendix B.2: the kite (41,15,67,15)/69 has 66 cylinders and 4 components without periodic orbit',
             s[0]['counts'] == {'cylinder': 66, 'no_periodic': 4})]


def chk_fold():
    t = load('s933_reflection_tilings.json')
    passing = sorted(k for k, v in t['feasibility'].items() if any(x.get('feasible') for x in v))
    return [('Proposition 5.3(b): of the 28 triangles exactly the three parents of Example 5.2 pass (a)',
             len(t['feasibility']) == 28 and passing == sorted(['(30, 41, 67)/138', '(23, 46, 201)/270', '(45, 90, 119)/254'])),
            ('Table 1: the four mirror-free octagons, their angles and side types',
             sorted((x['N_Pi'], x['angles'], x['types']) for v in t['tilings'].values() for x in v if x['k'] == 8 and not x['symmetries'])
             == sorted([(69, [134, 30, 67, 15, 82, 30, 41, 15], [2, 2, 2, 0, 0, 0, 0, 2]), (69, [82, 15, 134, 15, 41, 45, 67, 15], [0, 2, 2, 0, 0, 2, 2, 0]),
                        (127, [90, 238, 45, 45, 90, 90, 45, 119], [1, 1, 0, 0, 0, 0, 1, 1]),
                        (127, [45, 90, 45, 119, 135, 45, 45, 238], [0, 0, 1, 1, 0, 0, 1, 1])]))]


def chk_thmC():
    p = load('s935_pi_n_torus.json')
    iq = p['inequalities']
    e = load('s952_s1_edge_cert.json')
    mc = e['main_C7_0']
    pm = {int(m): r for m, r in e['per_m'].items()}
    a = load('s951_s1_allm_cert.json')
    f = load('s952_pi_family.json')['stage_B']
    return [('Theorem 5.4, m = 1: 254 conditions, 18 formally zero, 236 certified on (0, 0.2095], orders 0, 2, 4',
             (iq['total'], iq['formally_tight'], iq['certified_on_(0,pi/15]'], sorted(iq['leading_orders'])) == (254, 18, 236, [0, 2, 4])),
            ('Theorem 5.4, m = 1: the certificate continued to the edge pi/14 (Appendix B.1)', bool(f) and 'FORMAL' in str(f.get('identities', ''))),
            ('Theorem 5.4 (iii): 252 conditions, 16 formally zero, 2 more vanishing once w = z^{2m}, 30 cores, at most 7 boxes, '
             'smallest lower bound 0.0324', (mc['n_expr'], mc['zero'], len(mc['tight']), mc['distinct_cores'], mc['cores_ok'], mc['max_boxes'],
                                             round(mc['min_lb'], 4)) == (252, 16, 2, 30, 30, 7, 0.0324)),
            ('Theorem 5.4 (iii), m = 2, 3: 402 and 570 pairs certified, 18 formally zero, 6 end-zeros, no failure',
             [(pm[m]['certified'], pm[m]['tight'], pm[m]['end_zero'], pm[m]['fails']) for m in (2, 3)] == [(402, 18, 6, []), (570, 18, 6, [])]),
            ('Theorem 5.4: the certificates\' verdicts (s951 region and edge certificate, with their controls)', a['verdict'] is True and e['verdict'] is True),
            ('Appendix B.1 controls: past the edge, a lengthened run and an edge zero off its cyclotomic place each fail',
             all(bool(e[k]) for k in ('control_2d', 'control_per_m', 'control_edge_off', 'control_end_wrongK')))]


def chk_family():
    v = load('s952_s1_family.json')['verdict']
    return [('Theorem 5.4: the identities transferred to real x; no fold survivor (Lemma 5.5)', v['A1'] is True and v['A3'] is True
             and not any(v['A4_fold_survivors']))]


def chk_onset():
    v = load('s954_onset_converse.json')['verdict']
    return [('Remark 5.6: below the edge every S1(m; a pi/N) measured is periodic (188 rows); the control fires',
             v['rows'] == 188 and v['nonCP_below'] == [] and v['C_holds'] and v['ctrl_fires'])]


def chk_trap():
    c = [r for r in load('s932_transfer_census.json') if r['a'][0] <= r['a'][1] and 2 * r['a'][1] < r['N'] <= 15]
    acyc = [r for r in c if not r['cyclic']]
    five = load('s932_five_acyclic.json')
    out = [('Remark 6.9: 78 angle vectors, 32 acyclic, 27 of them with a0 = 1', (len(c), len(acyc), sum(r['a'][0] == 1 for r in acyc)) == (78, 32, 27)),
           ('Remark 6.9: (B) and (C) for the five without a pi/N corner (846, 3102, 3078, 230, 1760 wall families, none with W = 0)',
            five['all_ok'] and [r['fams'] for r in five['rows']] == [846, 3102, 3078, 230, 1760])]
    tags = ['5_1234', '1256', '1346', '9_1278', '9_1368', '9_1458', '11_12910', '11_13810', '11_14710', '11_15610', '13_121112',
            '13_131012', '13_14912', '13_15812', '13_16712']
    ok = True
    for t in tags:
        w = load(f's931_twist_walls_{t}.json')
        ok &= not w['fail'] and all(d['reverse_ends'] for d in w['infinite']) and w['pair_stats']['finite_coincidences'] == 0
    iso = [load(f's931_twist_walls_{t}.json') for t in ('5_1144', '7_1166', '9_1188', '11_111010')]
    out += [('Remark 6.10(a): for the 15 vectors with a0 = 1, N <= 13, two wall families meet only when one is the reverse of the other', ok),
            ('Remark 6.10(b): the isosceles controls — the computation detects meetings that are not reversals',
             all(not all(d['reverse_ends'] for d in w['infinite']) for w in iso))]
    return out


def chk_appA():
    out = []
    pm = load('s935_ladder_perm_cert.json')
    out.append(('Proposition A.1, m <= 6: every torus certified per m; the control at N = 8m + 13 is violated',
                sorted(int(m) for m in pm) == list(range(7)) and all(r['certified'] and r['returns'] and r['comb'] and r['control_viol'] > 0
                                                                    for r in pm.values())))
    g = load('s937_ladder_allm_cert.json')
    out.append(('Proposition A.1, m >= 7: 336 symbolic conditions certified; distinctness', g['verdict'] is True and g['main']['n_expr'] == 336))
    out.append(('Proposition A.3: the frontier of the tori', load('s941_chain_frontier.json')['verdict'] is True))
    q = load('s940_pin_frontier.json')
    out.append(('Appendix A: for odd 15 <= N <= 63 the 2(M + 1) tori are all the minimal components (measured)',
                sorted(r['N'] for r in q) == list(range(15, 64, 2)) and all(int(r['n_nonperiodic']) == 2 * ((r['N'] - 7) // 8) for r in q)))
    c = load('s945_corner_angles.json')['pi']
    out.append(('Appendix A, measured: for odd 15 <= N <= 63 each free component has one boundary cycle, three corners, genus 2M + 2',
                [o['N'] for o in c] == list(range(15, 64, 2)) and all(o['law_ok'] for o in c)
                and all(len(o['F']) == 2 and all(f['b'] == 1 and len(f['corners']) == 3 and f['g'] == 2 * ((o['N'] - 15) // 8) + 2
                                                 for f in o['F']) for o in c)))
    return out


def chk_steps():
    p = os.path.join(D, 'run_all_log.json')
    if not os.path.exists(p):
        return []
    log = load('run_all_log.json')
    return [(f'every script ran without error ({len(log)} steps; Propositions A.2 and the isosceles control assert their results)',
             all(x['rc'] == 0 for x in log))]


def chk_figs():
    return [('Figures 1-4 written', all(os.path.exists(os.path.join(HERE, 'figs', f'fig{i}.pdf')) for i in range(1, 5)))]


CHECKS = [(chk_kites, False), (chk_rhombus, False), (chk_asym, True), (chk_counts, False), (chk_catalogue, True),
          (chk_kites_noncp, True), (chk_fold, False), (chk_thmC, True), (chk_family, True), (chk_onset, False), (chk_trap, False),
          (chk_appA, True), (chk_figs, False), (chk_steps, False)]


def checks(quick):
    """(statement, passed) pairs: each script's verdict line and recorded output against the paper's figures."""
    res = []
    for f, long_ in CHECKS:
        if quick and long_:
            continue
        try:
            res += f()
        except Exception as ex:          # a missing or malformed record is a failure, not a crash
            res.append((f'{f.__name__}: {type(ex).__name__}: {ex}', False))
    return res


def main():
    quick = '--quick' in sys.argv
    if '--check' not in sys.argv:
        log = run_steps(quick)
        bad = [x for x in log if x['rc'] != 0]
        for x in bad:
            print('FAILED', x['script'], ' '.join(x['args']), 'exit', x['rc'])
    else:
        bad = []
    res = checks(quick)
    for s, ok in res:
        print('OK  ' if ok else 'FAIL', s)
    n_ok = sum(ok for _, ok in res)
    print(f'{n_ok}/{len(res)} checks pass; {len(bad)} step(s) failed')
    sys.exit(0 if not bad and n_ok == len(res) else 1)


if __name__ == '__main__':
    main()
