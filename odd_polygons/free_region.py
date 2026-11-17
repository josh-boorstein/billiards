"""[Deposit] Helper, imported by the scripts of this directory. Paper: "Perpendicular billiards in rational polygons of odd denominator". Research name: probes/s927_free_region.py, copied verbatim except imports (deposit names), paths (data/ and figs/ under odd_polygons/) and this paragraph. Run through run_all.py; the usage lines below are the research repo's. Bracketed IDs and section marks (§sNNN) refer to the project's internal ledgers.

s927 -- the FREE REGION of the perpendicular class at odd N: is every cylinder iota-invariant (P1),
and do the swapped-pair lengths satisfy Theorem P's sign-independence hypothesis?

NOT PRE-REGISTERED (exploratory, s927; the control and the counts were read after the first run).
Notation of general_triangles.md §s753 (2): X = S (2N copies, every vertex image marked), omega a
perpendicular-class direction, iota the reflection with axis perp omega, Fix(iota) = the vertical arcs
(at odd N exactly three side copies, one per side type).  An iota-invariant cylinder crosses Fix(iota)
(exactly twice); a FREE cylinder misses it.  The free region of §s927 = X minus the closure of the
saturation of Fix(iota).

MODES
  census LO HI OUT   every primitive triangle with odd N in [LO, HI]; flow_decomp.crossings with the
                     vertical arcs as keys; per row ncyl, nmin, inv (= cylinders crossing a vertical
                     arc), free, kmult (crossing counts).  P1 <=> free = 0; the cylinder count is
                     compared with R + T = (N + 3)/2 - p (§s753 Props S, B).
  even NMAX OUT      the control: every primitive even-N triangle (no a_i | N), every perpendicular
                     class.  Free cylinders here are expected (horizontal arcs, Lemma H of §s927).
  families OUT       Prop F'', the q = 5 family, W_{q,m} (q = 7, 9, 11), N <= 300: P1 and the sign test.
  formal MMAX OUT    Prop F'' bands m <= MMAX: every signed sum of the m prong lengths w_{m,k} (unfolded in floats,
                     A = pi/n held, B free) -- zero crossings where all its prongs are realised, and min |sum|.
  band OUT           every rational triangle (t, b, tn-t-b), t = 5, 7, 11, 13, N <= 300, on the m = 1, 2, 3 bands:
                     p vs the number of realised family prongs; P1; the sign test.
  sign IN OUT        the swapped pairs of each p >= 1 row of IN (a .jsonl census), traced exactly
                     (s624_slit_pairs.trace_all; pair = saddle connections exchanged by iota on germs),
                     and the exact test of Theorem P's hypothesis: no eps in {-1,0,1}^p \\ 0 with
                     sum eps_i * hol_i = 0 (all hol_i are positive multiples of omega).
COULD COME OUT OTHERWISE: a free cylinder or a minimal component at odd N (would refute (M+), resp.
Conj 2.2(b)); a signed relation among swapped-pair lengths (Theorem P silent there); the even-N
control showing no free cylinders (then P1 would not be parity-sensitive).
PRIOR ART: general_triangles.md §s753 (2) (Theorem C, Lemma L), §s754 (1) (CP <=> s' = 0), §s791 (1)/(3)
(formal side-parallel direction, every rational of I), §s919 (4) (the loop lemma); [NCYL-632] (mirror
pairs, slits), [NCYL-634] (side-parallel non-CP at odd d; Delecroix-Ruth #141); rulings.py --grep
side-parallel / invariant cylinder / free region / retrac -> nothing on iota-invariance of cylinders or
on a p = 1 certificate.
DATA: data/s912_census_53_99.jsonl (p, sprime; the 50 p >= 1 rows); data/s624_* (the even-N
catalogue); data_index --grep invariant -> none.
RESULT (s927; data/s927_free_region_*.json): census odd N <= 61, 3162 rows: every row CP, free = 0,
kmult = 2 everywhere, ncyl = (N + 3)/2 - p (p >= 1 on 19 rows); the 50 p >= 1 rows N = 57 ... 99: the
same.  Control even N <= 24 (87 CP class rows): 55 have free cylinders -- fires.  Sign test: 68 distinct
p >= 1 rows (19 + 50, (3,17,37)/57 in both; p up to 5) have NO signed relation, so Theorem P certifies (M+) on each
(CP itself is already decided by the same complete trace, CP <=> s' = 0).
Families (74 rows, N <= 297, p = m up to 16): free = 0, CP, p traced = p counted, NO signed relation.
Formal (float): on every Prop F'' band m <= 11 no signed sum of the family prong lengths crosses zero where its
prongs are realised; the minimum sits at the band's lower end B = 2mA and is > 0 (0.0021 at m = 6, lengths ~3).
Band points (96 rows, m <= 3, t = 5, 7, 11, 13): p = number of realised family prongs on 96/96, free = 0, CP, no relation.
-> general_triangles.md §s927.
"""
import _setup  # noqa: F401  (deposit: this directory on the path, cwd = odd_polygons/)
import itertools
import json
import math
import sys
import time

import flow_decomp as fd
import slit_pairs as sp


def vertical_keys(T, j):
    """Entry keys (g, s) of the one-crossing IET whose copy side is perpendicular to omega."""
    N, d = T.N, T.d
    iom = fd.mshift({j % N: 1}, d, N)
    pieces = fd._build_iet(T, j)[0]
    entry = {pc[3] for pc in pieces}
    keys = []
    for g in [(k, e) for e in (0, 1) for k in range(d)]:
        W = [T.g_apply(g, v) for v in T.V]
        for s, (i, k2) in enumerate(T.sides):
            e = fd.madd(W[k2], fd.mneg(W[i]))
            if not any(T.red(T.cross2(e, iom))) and (g, s) in entry:
                keys.append((g, s))
    return keys


def p1_row(abc, j, cap=20000):
    T = fd.Triangle(abc)
    keys = vertical_keys(T, j)
    out, _ = fd.crossings(abc, j, keys, cap)
    cyl = [o for o in out if o['kind'] == 'cylinder']
    inv = [o for o in cyl if any(o['k'])]
    return dict(abc=list(abc), j=j, nkeys=len(keys), ncyl=len(cyl), nmin=len(out) - len(cyl),
                inv=len(inv), free=len(cyl) - len(inv), kmult=sorted(sum(o['k']) for o in inv))


def triangles(N, skip_regular=False):
    for a in range(1, N):
        for b in range(a, N):
            c = N - a - b
            if c < b:
                break
            if math.gcd(math.gcd(a, b), c) != 1:
                continue
            if skip_regular and any(N % x == 0 for x in (a, b, c)):
                continue
            yield (a, b, c)


def swapped_pairs(abc, cap=30000):
    """Exact swapped pairs: (hol tuple, float length, steps) per pair; plus the number of capped seps."""
    T = fd.Triangle(abc)
    d = T.d
    j = sorted(set(fd.perp_dirs(abc)))[0]
    out = sp.trace_all(abc, j, cap)[0]
    assert (2 * j + 2 * d) % 4 == 0
    kr = ((2 * j + 2 * d) // 4) % d
    rg = lambda g: ((kr - g[0]) % d, 1 - g[1])     # iota = left multiplication by r
    sc = [o for o in out if o['end'] and not o.get('edge')]
    ncap = sum(1 for o in out if not o['end'])
    bystart = {(tuple(o['start'][0]), o['start'][1]): i for i, o in enumerate(sc)}
    pairs, seen = [], set()
    for i, o in enumerate(sc):
        # iota(gamma), read along omega, starts at iota(arrival germ) = (r g_arr, v_arr)
        k = bystart.get((rg(tuple(o['arr'][0])), o['arr'][1]))
        assert k is not None, 'iota image not traced'
        if k == i or (k, i) in seen:
            continue                                # type (2) or already paired
        seen.add((i, k))
        assert (tuple(sc[k]['arr'][0]), sc[k]['arr'][1]) == (rg(tuple(o['start'][0])), o['start'][1])
        assert tuple(sc[k]['hol']) == tuple(o['hol'])
        pairs.append((tuple(o['hol']), o['length'], o['steps']))
    return pairs, ncap


def _exact_zero(eps, hols):
    return not any(sum(e * h[k] for e, h in zip(eps, hols)) for k in range(len(hols[0])))


def sign_relation(hols, lengths=None, tol=1e-6):
    """A nonzero eps in {-1,0,1}^p with sum eps_i hol_i = 0 EXACTLY, or None.  p <= 9: brute force,
    exact.  Larger p (needs `lengths`): meet in the middle on the float lengths (every candidate with
    |sum eps_i l_i| < tol is confirmed or rejected exactly; a miss of a true relation would need a
    float error > tol on a sum of <= p lengths ~10, i.e. never at double precision)."""
    p = len(hols)
    if p <= 9 or lengths is None:
        for eps in itertools.product((-1, 0, 1), repeat=p):
            if any(eps) and next(e for e in eps if e) > 0 and _exact_zero(eps, hols):
                return list(eps)
        return None
    import bisect
    h1, h2 = p // 2, p - p // 2
    A = sorted((sum(e * l for e, l in zip(ea, lengths[:h1])), ea)
               for ea in itertools.product((-1, 0, 1), repeat=h1))
    keys = [a[0] for a in A]
    for eb in itertools.product((-1, 0, 1), repeat=h2):
        sb = sum(e * l for e, l in zip(eb, lengths[h1:]))
        k = bisect.bisect_left(keys, -sb - tol)
        while k < len(A) and keys[k] <= -sb + tol:
            eps = A[k][1] + eb
            if any(eps) and _exact_zero(eps, hols):
                return list(eps)
            k += 1
    return None


def main():
    mode = sys.argv[1]
    if mode == 'census':
        lo, hi, outp = int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
        rows = []
        for N in range(lo | 1, hi + 1, 2):
            for abc in triangles(N):
                r = p1_row(abc, sorted(set(fd.perp_dirs(abc)))[0])
                r['N'] = N
                rows.append(r)
        bad = [r for r in rows if r['free'] or r['nmin']]
        print(f'census odd N {lo}..{hi}: {len(rows)} rows, free/non-CP rows {len(bad)}, '
              f'kmult != 2: {sum(1 for r in rows if any(k != 2 for k in r["kmult"]))}, '
              f'p >= 1 rows: {sum(1 for r in rows if r["ncyl"] != (r["N"] + 3) // 2)}')
        json.dump(rows, open(outp, 'w'))
    elif mode == 'even':
        nmax, outp = int(sys.argv[2]), sys.argv[3]
        rows = []
        for N in range(6, nmax + 1, 2):
            for abc in triangles(N, skip_regular=True):
                for j in sorted(set(fd.perp_dirs(abc))):
                    r = p1_row(abc, j)
                    r['N'] = N
                    rows.append(r)
        cp = [r for r in rows if r['nmin'] == 0]
        print(f'even control N <= {nmax}: {len(rows)} class rows, CP {len(cp)}, '
              f'CP with free cylinders {sum(1 for r in cp if r["free"])}')
        json.dump(rows, open(outp, 'w'))
    elif mode == 'sign':
        inp, outp = sys.argv[2], sys.argv[3]
        src = [json.loads(l) for l in open(inp)] if inp.endswith('.jsonl') else json.load(open(inp))
        res = []
        for r in src:
            N = r.get('N', sum(r['abc']))
            p = r.get('p')
            if p is None and 'ncyl' in r:
                p = (N + 3) // 2 - r['ncyl']
            if not p:
                continue
            t = time.time()
            pairs, ncap = swapped_pairs(tuple(r['abc']))
            res.append(dict(abc=r['abc'], N=N, p=p, p_traced=len(pairs), ncapped=ncap,
                            lengths=[round(q[1], 6) for q in pairs],
                            relation=sign_relation([q[0] for q in pairs], [q[1] for q in pairs]), sec=round(time.time() - t, 1)))
            print(json.dumps(res[-1]), flush=True)
        print(f'sign test: {len(res)} rows, p mismatches {sum(1 for x in res if x["p"] != x["p_traced"])}, '
              f'relations {sum(1 for x in res if x["relation"])}')
        json.dump(res, open(outp, 'w'))


# ---------------------------------------------------------------------------------------------
# formal: Prop F'' prong lengths along their band (A = pi/n held, B free), float unfolding
# ---------------------------------------------------------------------------------------------
_SIDES = [(0, 1), (1, 2), (2, 0)]          # flow_decomp: s0 = V0V1, s1 = V1V2, s2 = V2V0


def f2_word(m, k):
    """Prop F'' bounce word w_{m,k} (general_triangles.md §s754 (3)), V2 -> V2."""
    return [0, 2] * (2*m - 1 - k) + [0, 1] + [0, 2] * (2*m + 2 + k) + [1] + [0, 2] * k + [0]


def _develop(A, B, w):
    import cmath
    C = math.pi - A - B
    T = [0j, complex(math.sin(C), 0), math.sin(B) * cmath.exp(1j * A)]
    start, Ts = T[2], [T]
    for s in w:
        a, b = T[_SIDES[s][0]], T[_SIDES[s][1]]
        d = (b - a) / abs(b - a)
        T = [a + d * ((P - a) / d).conjugate() for P in T]
        Ts.append(T)
    return start, T[2], Ts


def word_length(A, B, w):
    s, e, _ = _develop(A, B, w)
    return abs(e - s)                      # half flow_decomp's scale (sides sin, not 2 sin)


def word_realised(A, B, w, eps=1e-12):
    """the straight segment V2 -> V2' crosses every developed side of the word in its interior, in order."""
    start, end, Ts = _develop(A, B, w)
    d = end - start
    for T, s in zip(Ts[:-1], w):
        a, b = T[_SIDES[s][0]], T[_SIDES[s][1]]
        e = b - a
        den = d.real * (-e.imag) - d.imag * (-e.real)
        if abs(den) < 1e-300:
            return False
        r = a - start
        t = (r.real * (-e.imag) - r.imag * (-e.real)) / den
        u = (d.real * r.imag - d.imag * r.real) / den
        if not (eps < u < 1 - eps and 0 < t < 1):
            return False
    return True


def formal_scan(mmax, G=1500):
    """per m: zero crossings of every signed sum of the m prong lengths where all are realised, and the
    minimum |sum| (m <= 8).  Float evidence, not interval-certified."""
    import numpy as np
    out = []
    for m in range(1, mmax + 1):
        n = 6*m + 3
        A = math.pi / n
        u = np.linspace(0, 1, G + 2)[1:-1] ** 2
        Bs = 2*m*A + u * A
        L = np.array([[word_length(A, B, f2_word(m, k)) for B in Bs] for k in range(m)])
        R = np.array([[word_realised(A, B, f2_word(m, k)) for B in Bs] for k in range(m)])
        E_all = np.array([e for e in itertools.product((-1, 0, 1), repeat=m)
                          if sum(1 for x in e if x) >= 2 and next(x for x in e if x) > 0], dtype=float)
        nzero, best = 0, (math.inf, None, None)
        for c0 in range(0, len(E_all), 20000):
            E = E_all[c0:c0 + 20000]
            F = E @ L
            Rm = np.ones(F.shape, bool)
            for k in range(m):
                Rm &= (E[:, [k]] == 0) | R[[k]]
            nzero += int(((np.sign(F[:, :-1]) * np.sign(F[:, 1:]) < 0) & Rm[:, :-1] & Rm[:, 1:]).sum())
            if m <= 8:
                Fa = np.where(Rm, np.abs(F), np.inf)
                i, j = np.unravel_index(Fa.argmin(), Fa.shape)
                if Fa[i, j] < best[0]:
                    best = (float(Fa[i, j]), E[i].astype(int).tolist(), float(Bs[j] / A))
        upper = [float(Bs[R[k]].max() / A) if R[k].any() else None for k in range(m)]
        out.append(dict(m=m, n=n, band_upper_B_over_A=upper, zero_crossings=nzero,
                        min_abs=(best[0] if best[1] is not None else None) if m <= 8 else None, min_eps=best[1], min_at_B_over_A=best[2] if m <= 8 else None))
        print(json.dumps(out[-1])[:220], flush=True)
    return out


def band_points(ms=(1, 2, 3), ts=(5, 7, 11, 13), Nmax=300):
    """every rational triangle (t, b, tn - t - b) on the Prop F'' band of m (A = pi/n, B/A in (2m, 2m+1))."""
    res = []
    for m in ms:
        n = 6*m + 3
        for t in ts:
            N = t * n
            if N > Nmax:
                continue
            for b in range(2*m*t + 1, (2*m + 1)*t):
                c = N - t - b
                if math.gcd(math.gcd(t, b), c) != 1:
                    continue
                A, B = math.pi / n, b * math.pi / N
                nreal = sum(word_realised(A, B, f2_word(m, k)) for k in range(m))
                r = p1_row((t, b, c), sorted(set(fd.perp_dirs((t, b, c))))[0])
                pairs, ncap = swapped_pairs((t, b, c))
                res.append(dict(m=m, t=t, abc=[t, b, c], N=N, family_prongs_realised=nreal,
                                p=(N + 3)//2 - r['ncyl'], p_traced=len(pairs), ncapped=ncap, free=r['free'],
                                nmin=r['nmin'], relation=sign_relation([q[0] for q in pairs], [q[1] for q in pairs]) if pairs else None))
                print(json.dumps(res[-1])[:200], flush=True)
    return res


def family_rows():
    """Prop F'' (N = 18m+9), the q = 5 family (N = 30m+15) and the W_{q,m} bands q = 7, 9, 11 (N = 3q(2m+1)),
    every member with N <= 300 -- all have p = m (TRI-SWAP-FAM3, TRI-PENT-FAM, TRI-BAND-QM give p >= m)."""
    fams = []
    for m in range(1, 17):
        fams += [('F2a', m, (3, 6*m+2, 12*m+4)), ('F2b', m, (3, 6*m+1, 12*m+5))]
    for m in range(1, 10):
        fams += [('P1', m, (3, 12*m+5, 18*m+7)), ('P2', m, (3, 12*m+4, 18*m+8))]
    for q in (7, 9, 11):
        for m in range(1, 5):
            n = q*(2*m+1)
            K = (n - 2*m - 3)//2
            if 3*n <= 300:
                fams += [(f'Wq{q}', m, (3, b, 3*n - 3 - b)) for b in (3*K+1, 3*K+2)]
    return fams


if __name__ == '__main__':
    if sys.argv[1] == 'formal':
        json.dump(formal_scan(int(sys.argv[2])), open(sys.argv[3], 'w'))
    elif sys.argv[1] == 'band':
        res = band_points()
        print(f'band points: {len(res)} rows, p == realised family prongs {sum(1 for x in res if x["p"] == x["family_prongs_realised"])}, '
              f'free/non-CP {sum(1 for x in res if x["free"] or x["nmin"])}, relations {sum(1 for x in res if x["relation"])}')
        json.dump(res, open(sys.argv[2], 'w'))
    elif sys.argv[1] == 'families':
        res = []
        for name, m, abc in family_rows():
            t = time.time()
            N = sum(abc)
            r = p1_row(abc, sorted(set(fd.perp_dirs(abc)))[0])
            pairs, ncap = swapped_pairs(abc)
            res.append(dict(fam=name, m=m, abc=list(abc), N=N, ncyl=r['ncyl'], free=r['free'], nmin=r['nmin'],
                            p=(N + 3)//2 - r['ncyl'], p_traced=len(pairs), ncapped=ncap,
                            lengths=sorted(round(q[1], 6) for q in pairs),
                            relation=sign_relation([q[0] for q in pairs], [q[1] for q in pairs]) if pairs else None,
                            sec=round(time.time() - t, 1)))
            print(json.dumps(res[-1])[:200], flush=True)
        print(f'families: {len(res)} rows, free/non-CP {sum(1 for x in res if x["free"] or x["nmin"])}, '
              f'p mismatches {sum(1 for x in res if x["p"] != x["p_traced"])}, relations {sum(1 for x in res if x["relation"])}')
        json.dump(res, open(sys.argv[2], 'w'))
    else:
        main()
