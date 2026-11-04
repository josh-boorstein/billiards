"""Section 9.1: at a branch boundary grazing the right angle the limit orbit folds back on
itself, and at a boundary grazing an acute vertex it never does.

A boundary y between branches c and c' is a FOLD when the first-return map J, extended to y
from each side, returns y to itself: the limit orbit at y runs to the grazed vertex and back.
J is affine of slope +-1 on each branch, so its one-sided limits at y are endpoints of the
image intervals, and "= y" is decided exactly in Z[zeta_4Q].
Statements checked, at every recorded centre with Q <= 30 (the paper's population is the 126
centres of that range carrying its cylinder data; this is every one the census closed):
  * every boundary at an image of R is a fold from each returning side (the corner retroreflector:
    at R the two leg reflections compose to -id), and its two flanking branches are J-partners
    or one of them is the orphan -- so they are two different crossings, Prop 9.2(ii) at R;
  * no boundary at an image of O or A is a fold from either returning side.
  * on the row P = Q - 2 (Q odd), a graze at A -- the corner of angle pi/Q -- flanks the
    orphan exactly when Q = 1 (mod 4), and then exactly once; otherwise the orphan is the
    topmost branch (a measured pattern Section 9.1 states and does not prove).
The orphan's J is the identity, so on the orphan's side of a boundary "J returns y to itself" is
vacuous; folds are tested on the returning sides.  (Counted naively, the 130 acute boundaries
flanking an orphan at these centres would read as one-sided folds.)
"""
from common import FULL, check, done, record
from partition import Partition, im_cmp, Counters


def limits(T, L):
    """(J at the lower end, J at the upper end) of a branch, as vectors whose Im is the
    height; the orphan's J is the identity."""
    if not L['returns']:
        return L['lo'], L['hi']
    a, b = T.return_interval(L)
    if im_cmp(T.ring, L['t_copy'], L['o_copy'], Counters()) > 0:
        return a, b                        # orientation preserved: lo -> a, hi -> b
    return b, a                            # reversed: lo -> Rc - lo = b, hi -> a


def main():
    rec = record()
    rows = sorted((r for r in rec.values() if r['Q'] <= 30 and (FULL or r['steps'] < 10 ** 6)),
                  key=lambda r: (r['Q'], r['P']))
    print(f'Section 9.1 -- folds at {len(rows)} centres with Q <= 30\n')
    nR = nOA = foldR = partner = foldOA = 0
    orow = opat = 0
    for r in rows:
        T = Partition(r['P'], r['Q'])
        leaves, opened, _ = T.run()
        assert not opened
        R = T.ring
        lim = [limits(T, L) for L in leaves]
        J = {}
        for i, L in enumerate(leaves):
            if L['returns']:
                img = T.return_interval(L)
                J[i] = next(j for j, M in enumerate(leaves)
                            if T.eq_im(M['lo'], img[0]) and T.eq_im(M['hi'], img[1]))
        if r['P'] == r['Q'] - 2 and r['Q'] % 2:
            o = next(i for i, L in enumerate(leaves) if not L['returns'])
            flank = [src[1] for src in (leaves[o]['lo_src'], leaves[o]['hi_src'])
                     if src[0] == 'graze'].count('A')
            orow += 1
            opat += (flank == 1) if r['Q'] % 4 == 1 else (flank == 0 and o == len(leaves) - 1)
        for i in range(len(leaves) - 1):
            y = leaves[i]['hi']
            v = leaves[i]['hi_src'][1]
            # the orphan's J is the identity, so its side of a boundary returns y to itself
            # vacuously: a fold is tested on RETURNING sides only (None = the orphan side)
            fl = [T.eq_im(lim[j][e], y) if leaves[j]['returns'] else None
                  for j, e in ((i, 1), (i + 1, 0))]
            if v == 'R':
                nR += 1
                foldR += all(f is not False for f in fl)
                partner += (J.get(i) == i + 1 or not leaves[i]['returns']
                            or not leaves[i + 1]['returns'])
            else:
                nOA += 1
                foldOA += any(f is True for f in fl)
    check(f'P = Q - 2: an A-graze flanks the orphan iff Q = 1 (mod 4), then once; else the '
          f'orphan is topmost: {opat}/{orow}', opat == orow and orow > 0)
    check(f'boundaries at R that fold (every returning side): {foldR}/{nR}', foldR == nR)
    check(f'... whose flanking branches are J-partners or include the orphan: {partner}/{nR}',
          partner == nR)
    check(f'boundaries at O or A that fold from a returning side: {foldOA} of {nOA}', foldOA == 0)
    done()


if __name__ == '__main__':
    main()
