"""s254 -- PROVING (I)+(II): the general-`m` avoidance list, in closed form and case-free.

WHAT §s253 LEFT.  The interval law (`cf_width_laws.md` §s253 (1)) is verified 162/162 but not
proved; what it owes is §s239 (3)'s avoidance principle at general `m` -- per fan,
   (I)  every NON-exit vertex of the run is on the FAR side of the cell from the fan centre,
   (II) the exit vertex is on the centre's OWN side,
scored 2544/2544 by `probes/s253_wide_scan.py` but derived nowhere.  This probe verifies the
DERIVATION: closed forms for every fan centre, the reduction of the whole list to 16 named
families, a normal form for each, and the side conditions each proof step needs.

THE THREE THINGS THAT MAKE IT GO.
 1. THE CENTRES ARE SINUSOIDS.  With `al = P*th`, `rho = r*th`, `th = pi/2Q`,
    `Lam = 2 cos(al) sin(al-rho)/sin^2(al)`, `Gam = 2 cos(al) sin(rho)/sin^2(al)`:
        c_{2k+1} = Lam*(sin rho - sin(2k al + rho)),   c_{2k} = Gam*(sin(al-rho) + sin((2k-1)al+rho))
    and the DESCENDING arm is the mirror of the ascending one in the cell's lower edge:
        T_{2j+1} = 2X_m - c_{2(m-j)+1},   T_{2j} = 2X_m - c_{2(m-j+1)}.
 2. `X_m` AND `V_m` ARE THE COSINE AND SINE OF ONE PHASE.  With `be = al-2rho`,
    `C = 1-cos(al)cos(be)`, `D = cos(al)sin(be)`, `R = sqrt(C^2+D^2)`, `ph = atan2(D,C)` and the
    MASTER PHASE `Th_m = (2m+1)al + ph` (§s253 (4)'s phase):
        V_m = R sin(Th)/sin(al)          [§s253 (4)]
        X_m = cos(al) R (cos Th + cos(be+ph))/sin^2(al)
    -- because `R cos(be+ph) = cos be - cos al` and `R sin(be+ph) = sin be`.  So the cell's
    LOWER EDGE and its WIDTH are the two components of one rotating vector.
 3. TWO ELEMENTARY LEMMAS DISCHARGE 12 OF THE 16 FAMILIES.  For `0 <= A <= B <= pi`:
        L0  cos A - cos B >= 0
        L1a A <= B-al  =>  cos A - cos(al) cos B >= sin(al) sin B     >= 0
        L1b A <= B-al  =>  cos(al) cos A - cos B >= sin(al) sin(B-al) >= 0
    Every family is `R` times a two-cosine expression of exactly this shape.

THE RESULT.  Writing `Th = Th_m`, only FOUR of the sixteen families can fail, and each has a
closed form:
        POS  (X_m > 0)      <=>  Th + be + ph < pi
        A-II(k)             <=>  cos(al) cos Th + cos((2k-2)al + 2rho - ph) >= 0
        B-Ie(k)             <=>  cos Th + cos((2k-1)al + 2rho - ph) >= 0
        B-Io(k)              =   B-Ie(k) + sin^2(al) cos(2k al)/R   [so B-Ie => B-Io]
and the min over `k` of A-II/B-Ie sits at `k = 1` or `k = m`.  Hence
    THEOREM.  cone `2rho < al < 3rho`, `a1 >= 2m+1`, `(2m+2)al + ph <= pi/2`, and the warm-up
    `sin rho <= cos(al) sin(al-rho)`  ==>  the whole (I)+(II) list holds at rung `m`.
Past `Th + al = pi/2` the four residual conditions bite; they are what terminates the stack.

⚠ WHAT IS AND IS NOT PROVED.  (a) The step "(I)+(II) ==> the fan table" is §s239 (3)'s exit-rule
geometry, asserted at general `m`, NOT proved here.  (b) The equivalence "full list <=> the rung
EXISTS" is VERIFIED here (4000/4000 against `data/s253_wide_scan.json`), not proved -- and it is
NOT a route to the rung-COUNT law, which is the `n(P/Q)` §74-77 MATH WALL (`cf_width_laws.md`
§1291); this is one centre's stack, a smaller object.  (c) The three residual families are
verified, not proved, outside the theorem's range.

PRIOR ART: grepped `rulings.md` for 'fan'/'table'/'induction'/'centre'/'warm-up'/'avoidance'
-> [NB] s253 (no blind cone scan -- respected: the scan here scores families, and the
existence claim is checked only against ring-verified rows), [WFLOOR] s252/s253 (unimodality
modulo the table -- this probe is exactly the "modulo"), [SELECTOR] s252 (index by the fan
chain, done).  Grepped `cf_width_laws.md` for 'closed form.*centre'/'(star)'/'Dirichlet'/
'sinusoid' -> §s239 (4) owns the `m=1` six-fan case analysis with its hypotheses (*) `s_P>2H`
and (***) `2B1>S`, §s242 the `m=0` one, §s238 (3) records that both hypotheses are "finite-`Q`
forms of the cone itself, each with an O(1/m^2) correction" -- the warm-up `sin rho <= cos(al)
sin(al-rho)` found here is the same phenomenon at general `m`.  No closed form for a fan CENTRE,
and no general-`m` case analysis, exists anywhere in the ledger.  Grepped `f_leg_covariance.md`
for '4m+2'/'H-D1'.

Run:
  PYTHONPATH=.:engine:archive/scripts_2026-07:archive/scripts_2026-08:probes \
      .venv/bin/python3.13 -u probes/s254_avoidance_proof.py > logs/s254_avoidance_proof.log 2>&1
"""
import json
import math
from math import atan2, cos, gcd, inf, pi, sin, sqrt

from s253_proof_ineqs import chain_data

OUT = 'data/s254_avoidance_proof.json'
TOL = 1e-9


# ----------------------------------------------------------------- closed forms
class Cell:
    """Every closed form of §s254 for one `(P, r, a1, m)`."""

    def __init__(self, P, Q, r, a1, m):
        self.P, self.Q, self.r, self.a1, self.m = P, Q, r, a1, m
        th = pi / (2 * Q)
        self.al = al = P * th
        self.rho = rho = r * th
        self.be = be = al - 2 * rho
        self.k = k = cos(al)
        self.sa = sa = sin(al)
        self.C = C = 1 - k * cos(be)
        self.D = D = k * sin(be)
        self.R = R = sqrt(C * C + D * D)
        self.ph = ph = atan2(D, C)
        self.Th = Th = (2 * m + 1) * al + ph
        self.Lam = 2 * k * sin(al - rho) / sa ** 2
        self.Gam = 2 * k * sin(rho) / sa ** 2
        self.V = R * sin(Th) / sa
        self.X = k * R * (cos(Th) + cos(be + ph)) / sa ** 2
        self.w = min(self.V, self.X)
        self.U = self.X + self.w

    # ascending centres (odd / even fan index) and the descending mirror
    def S(self, j):                       # c_{2j+1}
        return self.Lam * (sin(self.rho) - sin(2 * j * self.al + self.rho))

    def c2(self, j):                      # c_{2j}
        return self.Gam * (sin(self.al - self.rho) + sin((2 * j - 1) * self.al + self.rho))

    def T(self, l):                       # descending fan `l` = chain fan 2m+1+l
        j, rem = divmod(l, 2)
        return 2 * self.X - (self.S(self.m - j) if rem else self.c2(self.m - j + 1))

    # the split heights, by lattice index
    def p(self, j):
        return self.k * sin(2 * j * self.al + 2 * self.rho) / self.sa

    def q(self, j):
        return self.k * sin(2 * j * self.al) / self.sa

    def F(self, j):
        return sin((2 * j - 1) * self.al + 2 * self.rho) / self.sa

    def G(self, j):
        return sin((2 * j + 1) * self.al) / self.sa

    def Dk(self, j):                      # X - c_{2j-1}
        return self.X - self.S(j - 1)

    def Ek(self, j):                      # c_{2j} - X
        return 2 * self.p(j - 1) - self.Dk(j)


# ------------------------------------------------- the 16 families (guarded list)
def families(z):
    """(tag, k, slack) for every avoidance inequality of the rung-`m` interval law."""
    m, a1, w, X, U = z.m, z.a1, z.w, z.X, z.U
    out = [('POS', 0, X), ('CEIL', 0, 1.0 - U)]
    for k in range(1, m + 1):                                  # asc odd fan 2k-1 (eps=+1)
        if 2 * a1 + 2 - 2 * k >= 2 * k:
            out.append(('A-Ie', k, z.p(k - 1) - z.Dk(k) - w))
        if 2 * a1 + 1 - 2 * k >= 2 * k + 1:
            out.append(('A-Io', k, z.F(k) - z.Dk(k) - w))
        out.append(('A-II', k, z.Dk(k) - z.F(k - 1)))
        if 2 * k <= 2 * a1 - 2 * k:                            # asc even fan 2k (eps=-1)
            out.append(('B-Ie', k, z.q(k) - z.Ek(k)))
        if 2 * k + 1 <= 2 * a1 - 2 * k - 1:
            out.append(('B-Io', k, z.G(k) - z.Ek(k)))
        out.append(('B-II', k, (z.Ek(k) - w) - z.G(k - 1)))
    if 2 * m + 2 <= 2 * a1 - 2 * m - 2:                        # APEX fan 2m+1 (eps=+1)
        out.append(('C-Ie', m + 1, z.q(m + 1) - z.Dk(m + 1) - w))
    if 2 * m + 3 <= 2 * a1 - 2 * m - 1:
        out.append(('C-Io', m + 1, z.F(m + 1) - z.Dk(m + 1) - w))
    for k in range(0, m + 1):                                  # desc l odd (eps=-1)
        if 2 * k + 2 <= 2 * a1 - 2 * k - 2:
            out.append(('D-Ie', k, z.q(k + 1) - z.Dk(k + 1)))
        if 2 * k + 1 <= 2 * a1 - 2 * k - 1:
            out.append(('D-Io', k, z.G(k) - z.Dk(k + 1)))
        out.append(('D-II', k, (z.Dk(k + 1) - w) - z.q(k)))
    for k in range(1, m + 1):                                  # desc l even (eps=+1)
        if 2 * a1 - 2 * k >= 2 * k + 2:
            out.append(('E-Ie', k, z.p(k) - z.Ek(k) - w))
        if 2 * a1 - 2 * k + 1 >= 2 * k + 1:
            out.append(('E-Io', k, z.F(k) - z.Ek(k) - w))
        out.append(('E-II', k, z.Ek(k) - z.p(k - 1)))
    return out


# ------------------------------------------------------------- the normal forms
def normal_forms(z, k):
    """(tag, direct slack, normal form) for every family, at fan index `k`.

    `w` is relaxed to `V` where that is sound (it always is for a LOWER bound, `w <= V`);
    the three families where the relaxation is LOSSY carry both readings.
    """
    al, rho, be, kk, sa, R, ph, Th, m = z.al, z.rho, z.be, z.k, z.sa, z.R, z.ph, z.Th, z.m
    s2 = sa * sa
    V, X = z.V, z.X
    chi = (2 * k - 2) * al + 2 * rho
    out = []
    a_ie = R * (kk * cos((2 * k - 1) * al + ph) - cos(2 * m * al + ph))
    out.append(('A-Ie', (z.p(k - 1) - z.Dk(k) - V) * s2, a_ie))
    out.append(('A-Io', (z.F(k) - z.Dk(k) - V) * s2, a_ie + s2 * cos(chi)))
    b_ie = kk * R * (cos(Th) + cos(chi + al - ph))
    out.append(('A-II', (z.Dk(k) - z.F(k - 1)) * s2, R * (kk * cos(Th) + cos(chi - ph))))
    out.append(('B-Ie', (z.q(k) - z.Ek(k)) * s2, b_ie))
    out.append(('B-Io', (z.G(k) - z.Ek(k)) * s2, b_ie + s2 * cos(2 * k * al)))
    out.append(('B-II', (z.Ek(k) - V - z.G(k - 1)) * s2,
                R * (cos(2 * k * al + ph) - cos(2 * m * al + ph))))
    d_io = R * (cos((2 * k - 2) * al + ph) - kk * cos(Th))
    out.append(('D-Io', (z.G(k - 1) - z.Dk(k)) * s2, d_io))
    out.append(('D-Ie', (z.q(k) - z.Dk(k)) * s2, d_io + s2 * cos(2 * (k - 1) * al + 2 * al)))
    out.append(('D-II', (z.Dk(k) - V - z.q(k - 1)) * s2,
                R * (cos(Th + al) + kk * cos((2 * k - 2) * al - be - ph))))
    e_io = R * (cos(chi - ph) + cos(Th + al))
    out.append(('E-Io', (z.F(k) - z.Ek(k) - V) * s2, e_io))
    out.append(('E-Ie', (z.p(k) - z.Ek(k) - V) * s2, e_io + s2 * cos(2 * k * al + 2 * rho)))
    out.append(('E-II', (z.Ek(k) - z.p(k - 1)) * s2, kk * R * (cos((2 * k - 1) * al + ph) - cos(Th))))
    # the three cell-level forms
    out.append(('POS', X * s2, kk * R * (cos(Th) + cos(be + ph))))
    out.append(('CEIL', (1.0 - X - V) * s2, R * (cos(ph) - cos(2 * m * al + ph))))
    # apex, and the two degenerate (w -> X) readings that the V-relaxation loses
    out.append(('C-Io', (z.F(m + 1) - z.Dk(m + 1) - z.w) * s2, s2 * cos(2 * m * al + 2 * rho) - z.w * s2))
    out.append(('D-IIx', (z.Dk(k) - X - z.q(k - 1)) * s2,
                kk * R * (cos((2 * k - 2) * al - be - ph) - cos(be + ph))))
    out.append(('E-Iox', (z.F(k) - z.Ek(k) - X) * s2, R * (cos(chi - ph) - kk * cos(be + ph))))
    return out


# --------------------------------------------------------------- the proof steps
def proof_certificate(z):
    """Every side condition the s254 proof actually invokes, at the indices it invokes it.

    Ranges and guards are those of `families()`; the three two-branch families carry the
    branch SELECTOR as well, so "the selected branch is valid" is what gets scored.
    """
    al, rho, be, ph, Th, m, a1 = z.al, z.rho, z.be, z.ph, z.Th, z.m, z.a1
    s2 = z.sa ** 2
    out = [('standing/cone-lo', al - 2 * rho), ('standing/cone-hi', 3 * rho - al),
           ('standing/a1-fix', 1e-12 - abs(a1 * al + rho - pi / 2)),
           ('standing/Th<pi', pi - Th),
           ('standing/be+ph<pi/2', pi / 2 - (be + ph)),          # <= R cos(be+ph)=cos be-cos al>0
           ('standing/(4m+2)al+2rho<=pi', pi - ((4 * m + 2) * al + 2 * rho)),
           ('warmup/star1  sin rho <= cos al sin(al-rho)', cos(al) * sin(al - rho) - sin(rho)),
           ('CEIL/L0', 2 * m * al)]
    for k in range(1, m + 1):
        chi = (2 * k - 2) * al + 2 * rho
        if 2 * a1 + 2 - 2 * k >= 2 * k:
            out.append(('A-Ie/L1b', (2 * m * al + ph - al) - ((2 * k - 1) * al + ph)))
        if 2 * a1 + 1 - 2 * k >= 2 * k + 1:
            out.append(('A-Io/cos>=0', pi / 2 - chi))
        if 2 * k + 1 <= 2 * a1 - 2 * k - 1:
            out.append(('B-Io/cos>=0', pi / 2 - 2 * k * al))
        out.append(('B-II/L0', 2 * (m - k) * al))
        out.append(('E-II/L0', Th - ((2 * k - 1) * al + ph)))
        if 2 * a1 - 2 * k >= 2 * k + 2:
            out.append(('E-Ie/cos>=0', pi / 2 - (2 * k * al + 2 * rho)))
        if 2 * a1 - 2 * k + 1 >= 2 * k + 1:                      # E-Io, two branches
            if chi >= ph:                                        # branch V
                out.append(('E-Io/V:|chi-ph|+Th+al<=pi', pi - ((chi - ph) + Th + al)))
            else:                                                # branch X
                out.append(('E-Io/X:|chi-ph|<=ph-2rho', (ph - 2 * rho) - (ph - chi)))
    for k in range(0, m + 1):
        t = 2 * k * al
        if 2 * k + 2 <= 2 * a1 - 2 * k - 2:
            out.append(('D-Ie/cos>=0', pi / 2 - rho - (t + 2 * al)))
        if 2 * k + 1 <= 2 * a1 - 2 * k - 1:
            out.append(('D-Io/L1a', (Th - al) - (t + ph)))
        if t <= 2 * (be + ph):                                   # D-II branch X
            out.append(('D-II/X:t<=2(be+ph)', 2 * (be + ph) - t))
        else:                                                    # D-II branch V, via L1b
            out.append(('D-II/V:t+(2m+3)al<=pi+be+ph', (pi + be + ph) - (t + (2 * m + 3) * al)))
    if 2 * m + 2 <= 2 * a1 - 2 * m - 2:
        out.append(('C-Ie/cos>=0', pi / 2 - rho - (2 * m + 2) * al))
    if 2 * m + 3 <= 2 * a1 - 2 * m - 1:                          # C-Io <=> w <= cos(Th-be-ph)
        # V-branch: a sin Th + b cos Th >= 0 with a = cos(al)sin(al-rho)-sin(rho) > 0 (star1),
        # b = sin(al)sin(al-rho) > 0  <=>  Th + psi0 <= pi,  psi0 = atan2(b,a) in (0,pi/2).
        psi0 = atan2(sin(al) * sin(al - rho), cos(al) * sin(al - rho) - sin(rho))
        out.append(('C-Io/V: Th + psi0 <= pi', pi - (Th + psi0)))
    return out


# ------------------------------------------------------------------------- scan
def main():
    ex, seen = set(), set()
    for row in json.load(open('data/s253_wide_scan.json')):
        ex.add((row['P'], row['r'], row['a1'], row['m']))
        seen.add((row['P'], row['r'], row['a1']))

    fam_ag, nf_max, pc_ag = {}, {}, {}
    cf_max = {'S': 0.0, 'c2': 0.0, 'T': 0.0, 'X': 0.0, 'V': 0.0, 'U': 0.0}
    fam_vs_brute, nrow = 0, 0
    thm = {'in_range': 0, 'in_range_ok': 0, 'ok': 0, 'tot': 0, 'star1_fail': 0}
    conf = {'TT': 0, 'TF': 0, 'FT': 0, 'FF': 0}
    endpt = [0, 0]
    star1_bad, interval = [], [0, 0]
    okset = {}
    for P in range(5, 41):
        for r in range(1, P):
            if gcd(P, r) != 1 or not (2 * r < P < 3 * r):
                continue
            for a1 in range(3, 20):
                Q = a1 * P + r
                for m in range(0, (a1 - 1) // 2 + 1):
                    z = Cell(P, Q, r, a1, m)
                    nrow += 1
                    # T1 -- closed forms against s253's hand-off recursion
                    g, A, lam, fans, X, U = chain_data(m, P, Q, r, a1)
                    for j in range(0, m + 1):
                        cf_max['S'] = max(cf_max['S'], abs(fans[2 * j]['c'] - z.S(j)))
                    for j in range(1, m + 1):
                        cf_max['c2'] = max(cf_max['c2'], abs(fans[2 * j - 1]['c'] - z.c2(j)))
                    for l in range(1, 2 * m + 2):
                        cf_max['T'] = max(cf_max['T'], abs(fans[2 * m + l]['c'] - z.T(l)))
                    cf_max['X'] = max(cf_max['X'], abs(X - z.X))
                    cf_max['U'] = max(cf_max['U'], abs(U - z.U))
                    # T2 -- the family list reproduces the brute-force per-fan (I)/(II)
                    bI, bII = inf, inf
                    for f in fans:
                        c, eps = f['c'], f['eps']
                        dn = (X - c) if eps > 0 else (c - U)
                        df = (U - c) if eps > 0 else (c - X)
                        for n in f['ns'][:-1]:
                            bI = min(bI, lam(n) * A(n) - df)
                        bII = min(bII, dn - lam(f['n_ex']) * A(f['n_ex']))
                    fl = families(z)
                    fI = min([s for t, _, s in fl if t[-2:] in ('Ie', 'Io')] + [inf])
                    fII = min([s for t, _, s in fl if t.endswith('II')] + [0.0])
                    if abs(bI - fI) > 1e-9 or abs(bII - fII) > 1e-9:
                        fam_vs_brute += 1
                    # T3/T4 -- normal forms and proof-step conditions
                    for k in range(1, max(m, 1) + 2):
                        for tag, direct, nf in normal_forms(z, k):
                            nf_max[tag] = max(nf_max.get(tag, 0.0), abs(direct - nf))
                    for tag, s in proof_certificate(z):
                        a = pc_ag.setdefault(tag, [0, 0, inf])
                        a[0] += 1
                        a[2] = min(a[2], s)
                        if s < -1e-12:
                            a[1] += 1
                    # T5/T6 -- families, the theorem's range, the residual
                    allok = True
                    for tag, k, s in fl:
                        a = fam_ag.setdefault(tag, [0, 0, inf])
                        a[0] += 1
                        a[2] = min(a[2], s)
                        if s < -TOL:
                            a[1] += 1
                            allok = False
                    thm['tot'] += 1
                    thm['ok'] += allok
                    okset.setdefault((P, r, a1), []).append((m, allok))
                    star1 = cos(z.al) * sin(z.al - z.rho) - sin(z.rho)
                    thm['star1_fail'] += star1 < 0
                    if star1 < 0 and (P, r, a1) not in star1_bad:
                        star1_bad.append((P, r, a1))
                    if z.Th <= pi / 2:
                        thm['in_range'] += 1
                        thm['in_range_ok'] += allok
                    # endpoint reduction of the two residual families
                    for tag in ('A-II', 'B-Ie'):
                        vals = [s for t, _, s in fl if t == tag]
                        if len(vals) > 2:
                            endpt[0] += 1
                            endpt[1] += min(vals) >= min(vals[0], vals[-1]) - 1e-12
                    # T7 -- against ring existence
                    if (P, r, a1) in seen and 5 <= a1 <= 14 and m <= 6:
                        e = (P, r, a1, m) in ex
                        conf[('T' if allok else 'F') + ('T' if e else 'F')] += 1

    print(f'=== s254: the general-m avoidance list, closed form.  {nrow} (P,r,a1,m) rows')
    print(f'    P=5..40 in the a2=2 cone, a1=3..19, m=0..floor((a1-1)/2)\n')
    print('T1  closed forms vs the s253 hand-off recursion (max |err|):')
    print('    ' + '  '.join(f'{k}={v:.1e}' for k, v in cf_max.items()))
    print(f'\nT2  family list == brute-force per-fan (I)/(II):  mismatches = {fam_vs_brute}/{nrow}')
    print('\nT3  normal-form identities (max |direct - normal form|):')
    for t in sorted(nf_max):
        print(f'    {t:6s} {nf_max[t]:.2e}')
    print('\nT4  proof-step side conditions (each must never fail):')
    for t in sorted(pc_ag):
        n, f, mn = pc_ag[t]
        print(f'    {t:12s} {f:6d}/{n:<7d} fail   min slack {mn:+.3e}')
    print('\nT5  families over the blind cone (a1 >= 2m+1 only):')
    for t in sorted(fam_ag):
        n, f, mn = fam_ag[t]
        print(f'    {t:6s} {f:7d}/{n:<7d} fail   min slack {mn:+.3e}')
    resid = {'POS', 'A-II', 'B-Ie', 'B-Io'}
    uncond = sum(v[1] for t, v in fam_ag.items() if t not in resid)
    print(f'\nT6  UNCONDITIONAL families (all but POS/A-II/B-Ie/B-Io): {uncond} failures')
    print(f'    Th <= pi/2 sub-range: {thm["in_range_ok"]}/{thm["in_range"]} rows hold the WHOLE list; '
          f'{thm["ok"]}/{thm["tot"]} rows hold it overall; star1 fails on {thm["star1_fail"]} rows')
    print(f'    endpoint reduction of A-II/B-Ie (min over k at k=1 or k=m): {endpt[1]}/{endpt[0]}')
    for key, lst in okset.items():
        lst.sort()
        flags = [b for _, b in lst]
        interval[0] += 1
        interval[1] += flags == sorted(flags, reverse=True)     # an INITIAL interval of m
    print(f'    the ok-set is an INITIAL interval of m: {interval[1]}/{interval[0]} centres')
    print(f'    star1 (the only warm-up) fails on {len(star1_bad)} (P,r,a1): {star1_bad[:12]}')
    print(f'\nT7  full list vs RING existence (s253 wide-scan rows): {conf}')
    with open(OUT, 'w') as fh:
        json.dump({'rows': nrow, 'closed_forms': cf_max, 'fam_vs_brute': fam_vs_brute,
                   'normal_forms': nf_max,
                   'proof_conditions': {k: v[:2] + [v[2]] for k, v in pc_ag.items()},
                   'families': {k: v[:2] + [v[2]] for k, v in fam_ag.items()},
                   'theorem': thm, 'endpoint': endpt, 'existence': conf}, fh, indent=1)
    print(f'\n-> {OUT}')


if __name__ == '__main__':
    main()
