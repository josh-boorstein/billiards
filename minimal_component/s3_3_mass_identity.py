"""Section 3.3: the cylinder table, the Kac integral, and Proposition 5.

Statements checked:
  * The three cylinders (from the six branches of `s4_branch_count.py`): return count
    N = 25 / 53 / 9 (core orbit of N + 1 = 26 / 54 / 10 reflections), height = the pair
    sum of the two branch widths, return length l = the length of the first-return path
    to L1 in the perpendicular direction (the half-orbit up to the head-on hit, constant
    on the branch), area mass = height x l.  Every entry as printed (truncated).
  * Heights exactly 2cos(2pi/5), 1 - 2cos(2pi/5) - 2sin(pi/30), 2sin(pi/30); the three
    l / cot(alpha) exactly 1 + 4cos(pi/15) + 2cos(pi/5), 4 + 4cos(pi/15) + 6cos(2pi/15),
    2 + 2cos(2pi/15) - 2cos(pi/5); the ratios l / N = 0.2352..., 0.2275..., 0.2210....
  * Sigma_cyl / cot(alpha) = 15 g_len = 3 - 4cos(pi/15) + 12cos(2pi/15) - 4cos(pi/5)
    = 6.81388711127619849790209...; Sigma_cyl = 6.135251512382014...
  * The minimal component's mass from the interval exchange alone (Kac on Sigma_Hbar,
    flux |cos theta_Hbar| ds = cos(4pi/15) ds, roof constant on each of the seven words):
    mu(M) = 7.37080915208558428507..., mu(M)/cot(alpha) =
    12 + 4cos(pi/15) - 12cos(2pi/15) + 4cos(pi/5); the widths w_j and roofs r_j of the
    seven words as tabulated, with w1 = l1, w2 = l2, w3 + ... + w7 = l3.
  * Proposition 5: Sigma_cyl / cot(alpha) + mu(M) / cot(alpha) = Q = 15 exactly;
    Sigma_cyl / cot(alpha) = 6 l1 + 7 l2 + 8 l3 and mu(M) / cot(alpha) = 9 l1 + 8 l2 + 7 l3.

INDEPENDENCE.  The cylinder side is computed from the six horizontal branch words on L1
and nothing else; the Kac side from the seven return words on Sigma_Hbar and nothing
else.  No quantity of one enters the other, so their sum being the integer 15 is a
check that could have failed.
"""
from mpmath import nstr

from t0 import T, F, Q, KAPPA, M_H, M_PERP, cpi15, c, return_words, check, done, trunc
from s3_1_return_map import lambdas
from s4_branch_count import branches, cylinders


def cylinder_side():
    rows, tiles, _ = branches()
    assert tiles
    out = []
    for a, b in cylinders(rows):
        height = (a['hi'] - a['lo']) + (b['hi'] - b['lo'])
        la, ta = T.path_length(a['steps'])
        lb, tb = T.path_length(b['steps'])
        assert ta.is_zero() and tb.is_zero() and la == lb
        out.append({'N': 2 * len(a['word']) - 1, 'height': height, 'ell': la,
                    'area': height * la})
    return out


def kac_side():
    """mu(M) = cos(4pi/15) * sum_j int_{lo_j}^{hi_j} r_j(s) ds over the seven words."""
    words = return_words()
    rows, tiles, _ = T.tile(words, M_H)
    assert tiles
    acc = F.zero
    table = []
    for r in rows:
        T0_, T1_ = T.path_length(r['steps'])
        w = r['hi'] - r['lo']
        acc = acc + T0_ * w + T1_ * (r['hi'] * r['hi'] - r['lo'] * r['lo']) * (F.one / 2)
        table.append({'N': len(r['word']), 'w': w, 'r': T0_, 'r_const': T1_.is_zero()})
    flux = c(8)                                          # |cos(11pi/15)| = cos(4pi/15)
    return flux * acc, table


def main():
    print('Section 3.3 -- the mass identity\n')
    cyl = sorted(cylinder_side(), key=lambda d: {25: 0, 53: 1, 9: 2}[d['N']])
    paper = [(25, '0.618033988749', '5.880200614063', '3.634163840159'),
             (53, '0.172909084714', '12.059888547943', '2.085264290587'),
             (9, '0.209056926535', '1.989043790736', '0.415823381635')]
    print('  cylinder  N    height          l                area')
    for k, (d, p) in enumerate(zip(cyl, paper), 1):
        row = (d['N'], trunc(d['height'].value(), 12), trunc(d['ell'].value(), 12),
               trunc(d['area'].value(), 12))
        print(f'  C{k}        {row[0]:<4d} {row[1]}  {row[2]:<16s} {row[3]}')
        check(f'C{k} row as printed', row == p)
    check('core orbits make N + 1 = 26 / 54 / 10 reflections',
          [d['N'] + 1 for d in cyl] == [26, 54, 10])
    ratios = [trunc((d['ell'].value() / d['N']), 4) for d in cyl]
    check(f'l / N = {ratios}', ratios == ['0.2352', '0.2275', '0.2210'])
    sin30 = c(14)
    cos2pi5 = c(12)
    check('heights 2cos(2pi/5), 1 - 2cos(2pi/5) - 2sin(pi/30), 2sin(pi/30)',
          [d['height'] for d in cyl] == [2 * cos2pi5, 1 - 2 * cos2pi5 - 2 * sin30, 2 * sin30])
    want_l = [1 + 4 * cpi15(1) + 2 * cpi15(3), 4 + 4 * cpi15(1) + 6 * cpi15(2),
              2 + 2 * cpi15(2) - 2 * cpi15(3)]
    check('l / cot(alpha) as printed', all(d['ell'] / KAPPA == wl for d, wl in zip(cyl, want_l)))
    sigma = cyl[0]['area'] + cyl[1]['area'] + cyl[2]['area']
    g15 = sigma / KAPPA
    check(f'Sigma_cyl = {trunc(sigma.value(), 15)}...', trunc(sigma.value(), 15) == '6.135251512382014')
    check('15 g_len = Sigma_cyl / cot(alpha) = 3 - 4cos(pi/15) + 12cos(2pi/15) - 4cos(pi/5)',
          g15 == 3 - 4 * cpi15(1) + 12 * cpi15(2) - 4 * cpi15(3))
    check(f'  = {trunc(g15.value(), 23)}...', trunc(g15.value(), 23) == '6.81388711127619849790209')

    mu, table = kac_side()
    l1, l2, l3 = lambdas()
    print('\n  word  N    width w_j         roof r_j')
    pw = ['0.581886146929', '0.022340594865', '0.116000482173', '0.093056444361',
          '0.045284632676', '0.070715849496', '0.070715849496']
    pr = ['2.972579301909', '22.108760536759', '19.136181234849', '19.136181234849',
          '19.136181234849', '27.750141644774', '27.750141644774']
    for j, (t, w_, r_) in enumerate(zip(table, pw, pr), 1):
        got = (trunc(t['w'].value(), 12), trunc(t['r'].value(), 12))
        print(f'  {j}     {t["N"]:<4d} {got[0]}    {got[1]}')
        check(f'word {j} row as printed', got == (w_, r_) and t['r_const'])
    w = [t['w'] for t in table]
    check('w1 = l1, w2 = l2, w3 + ... + w7 = l3',
          w[0] == l1 and w[1] == l2 and w[2] + w[3] + w[4] + w[5] + w[6] == l3)
    muk = mu / KAPPA
    # the paper prints this one ROUNDED (7.37080915208558428506976...), unlike the tables
    check(f'mu(M) = {nstr(mu.value(), 21)}... (rounded)',
          nstr(mu.value(), 21) == '7.37080915208558428507')
    check('mu(M)/cot(alpha) = 12 + 4cos(pi/15) - 12cos(2pi/15) + 4cos(pi/5)',
          muk == 12 + 4 * cpi15(1) - 12 * cpi15(2) + 4 * cpi15(3))
    check('PROPOSITION 5: Sigma_cyl/cot(alpha) + mu(M)/cot(alpha) = 15 exactly', g15 + muk == Q)
    check('Sigma_cyl/cot(alpha) = 6 l1 + 7 l2 + 8 l3', g15 == 6 * l1 + 7 * l2 + 8 * l3)
    check('mu(M)/cot(alpha) = 9 l1 + 8 l2 + 7 l3', muk == 9 * l1 + 8 * l2 + 7 * l3)
    done()


if __name__ == '__main__':
    main()
