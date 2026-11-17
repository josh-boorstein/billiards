"""[Deposit] Helper, imported by the scripts of this directory. Paper: "Perpendicular billiards in rational polygons of odd denominator". Research name: probes/s910_poly_slots.py, copied verbatim except imports (deposit names), paths (data/ and figs/ under odd_polygons/) and this paragraph. Run through run_all.py; the usage lines below are the research repo's. Bracketed IDs and section marks (§sNNN) refer to the project's internal ledgers.

s910 (`s755-i`): the slot identity and `p = 0 ⟹ CP` beyond triangles -- a prong census in
rational POLYGONS (quadrilaterals, pentagons, one non-convex), directions exact, positions mpmath.

PRE-REGISTERED (s910, before any polygon was traced).
  Claim under test (hand reading of the polygon note, `tex/polygon/polygon.tex` §2–§4, §7): the
  proofs of Thm 4.2 (slot identity) and Thm 7.8 (`p = 0` ⟹ CP) use the triangle only through the
  slot count. For a rational k-gon with angles pi*a_i/N (gcd(a, N) = 1, Σ a_i = (k − 2)N) and a
  class C of perpendicular directions:
    (S) S = (k − 2)N at odd N;  S = (k − 2)N/2 − h at even N  (h = #sides whose direction is in C);
    (B) Σ_{σ ⊥ C} n_σ = v + S − 2p − s'   (each singular point of a vertical side = one slot whose
        prong meets a side perpendicularly before any vertex, on that side);
    (P) the other prongs pair by reversal (p pairs) or never reach a vertex (s');
    (T) p = 0 ⟹ s' = 0.
  Predictions: (i) the slot census equals (S) on every row; (ii) per vertical side, the prong
  perpendicular-hit points are exactly the discontinuities of the beam's side word, checked by
  sampling the beam (every sampled word change brackets >= 1 prong point; every prong point that is
  alone in its sampling gap separates different words); (iii) the reversal pairing closes on every
  row; (iv) on every row with p = 0 the census finds s' = 0 at the cap.
  Could fail: (ii) compares two instruments (vertex prongs vs. side launches); a triangle-only step
  in the bijection shows as a mismatch. (iv) can fire: a generic quadrilateral with fixed angles is
  a one-parameter family whose perpendicular class need not be CP; any capped prong on a p = 0 row
  refutes (T) (or flags the cap). Triangles (2,3,4)/9 (p = s' = 0, Σ n = 12) and (2,5,8)/15
  (Σ n = N + 1, p = 1) are the validation rows (`data/s753_slot_count.json`).
  PRIOR ART: `rulings.py --grep polygon` → only [LIT-012] (EMMW quadrilaterals, no cylinder count);
  `paper_polygon.md` port row P: "Prop S, Theorem C beyond triangles — not checked"; grep of
  `orphan_theorem.md`, `necklace_ledger.md`, `literature_review.md` for quadrilateral / k-gon /
  (k − 2)N → none. The literature half (counts of perpendicular orbits in polygons) is the s910
  literature agent's question (B), `audits/s910_lit_open_problems.md`.
  DATA: `data_index.py` has no polygon (k >= 4) census; `data/s753_slot_count.json` is the
  triangle reference.
Instrument: directions are exact integers mod 4N (unit pi/(2N)); reflection in side j: d -> 2m_j − d;
perpendicular to side j iff d ≡ m_j + N (mod 2N). Positions in mpmath (60 digits), vertex hit =
within 1e-30 of a vertex. Output: data/s910_poly_slots.json.
Run from the repo root: PYTHONPATH=.:engine:archive/scripts_2026-07:archive/scripts_2026-08:probes
"""
import _setup  # noqa: F401  (deposit: this directory on the path, cwd = odd_polygons/)
import json
import random
import sys
from math import gcd

import mpmath as mp

mp.mp.dps = 60
TOL = mp.mpf(10) ** -30


class Poly:
    def __init__(self, a, N, lengths=None, seed=0):
        k = len(a)
        assert sum(a) == (k - 2) * N and gcd(gcd(*a), N) == 1 if k > 2 else True
        self.a, self.N, self.k = a, N, k
        # side j runs V_j -> V_{j+1}; direction index (unit pi/(2N)), CCW polygon
        idx = [0]
        for j in range(1, k):
            idx.append(idx[-1] + 2 * (N - a[j]))
        self.idx = [x % (4 * N) for x in idx]
        self.m = [x % (2 * N) for x in self.idx]          # line direction mod pi
        u = [self.unit(x) for x in self.idx]
        rnd = random.Random(seed)
        if lengths is None:
            for _try in range(3000):
                L = [mp.mpf(rnd.uniform(0.5, 1.5)) for _ in range(k - 2)]
                sx = sum(L[i] * u[i][0] for i in range(k - 2))
                sy = sum(L[i] * u[i][1] for i in range(k - 2))
                # solve x*u[k-2] + y*u[k-1] = -(sx, sy)
                A = mp.matrix([[u[k - 2][0], u[k - 1][0]], [u[k - 2][1], u[k - 1][1]]])
                sol = mp.lu_solve(A, mp.matrix([-sx, -sy]))
                if sol[0] > 0.05 and sol[1] > 0.05:
                    L += [sol[0], sol[1]]
                    break
            else:
                raise ValueError("no closing lengths found (angles may force a non-positive side)")
            lengths = L
        self.L = lengths
        V = [(mp.mpf(0), mp.mpf(0))]
        for j in range(k - 1):
            V.append((V[-1][0] + lengths[j] * u[j][0], V[-1][1] + lengths[j] * u[j][1]))
        self.V = V
        closex = V[-1][0] + lengths[k - 1] * u[k - 1][0]
        closey = V[-1][1] + lengths[k - 1] * u[k - 1][1]
        assert abs(closex) < TOL and abs(closey) < TOL, "polygon does not close"
        self.check_simple()

    def unit(self, d):
        ang = mp.pi * d / (2 * self.N)
        return (mp.cos(ang), mp.sin(ang))

    def check_simple(self):
        k = self.k
        for i in range(k):
            for j in range(i + 2, k):
                if i == 0 and j == k - 1:
                    continue
                if seg_cross(self.V[i], self.V[(i + 1) % k], self.V[j], self.V[(j + 1) % k]):
                    raise ValueError("not simple")

    def classes(self):
        """the perpendicular classes, as sets of direction indices mod 4N"""
        N = self.N
        out = []
        for j in range(self.k):
            d0 = (self.idx[j] + N) % (4 * N)
            orb = {d0}
            todo = [d0]
            while todo:
                d = todo.pop()
                for mj in self.m:
                    e = (2 * mj - d) % (4 * N)
                    if e not in orb:
                        orb.add(e)
                        todo.append(e)
            if orb not in out:
                out.append(orb)
        return out

    def slots(self, C):
        N, k = self.N, self.k
        out = []
        for i in range(k):
            lo = self.idx[i]
            for t in range(1, 2 * self.a[i]):
                d = (lo + t) % (4 * N)
                if d in C:
                    out.append((i, d))
        return out

    def first_hit(self, P, d, skip):
        """first side hit by the ray from P in direction d; skip = set of sides/vertex to ignore at t~0"""
        ux, uy = self.unit(d)
        best = None
        for j in range(self.k):
            A, B = self.V[j], self.V[(j + 1) % self.k]
            ex, ey = B[0] - A[0], B[1] - A[1]
            den = ux * ey - uy * ex
            if abs(den) < TOL:
                continue
            wx, wy = A[0] - P[0], A[1] - P[1]
            t = (wx * ey - wy * ex) / den
            s = (wx * uy - wy * ux) / den
            if t <= TOL * 10 or s < -TOL or s > 1 + TOL:
                continue
            if best is None or t < best[0]:
                best = (t, j, s)
        t, j, s = best
        X = (P[0] + t * ux, P[1] + t * uy)
        if s < TOL * 1e3:
            return ('vertex', j, X)
        if s > 1 - TOL * 1e3:
            return ('vertex', (j + 1) % self.k, X)
        return ('side', j, X, s)

    def prong(self, i, d, cap):
        P = self.V[i]
        N = self.N
        for step in range(cap):
            h = self.first_hit(P, d, None)
            if h[0] == 'vertex':
                return ('vertex', h[1], d, step)
            _, j, X, s = h
            if (d - self.m[j] - N) % (2 * N) == 0:
                return ('perp', j, s, step)
            d = (2 * self.m[j] - d) % (4 * N)
            P = X
        return ('cap', None, None, cap)

    def beam_word(self, j, s, cap):
        A, B = self.V[j], self.V[(j + 1) % self.k]
        P = (A[0] + s * (B[0] - A[0]), A[1] + s * (B[1] - A[1]))
        N = self.N
        d = (self.idx[j] + N) % (4 * N)
        word = []
        for step in range(cap):
            h = self.first_hit(P, d, None)
            if h[0] == 'vertex':
                return ('sing', tuple(word))
            _, jj, X, ss = h
            word.append(jj)
            if (d - self.m[jj] - N) % (2 * N) == 0:
                return ('ok', tuple(word))
            d = (2 * self.m[jj] - d) % (4 * N)
            P = X
        return ('cap', tuple(word))


def seg_cross(p1, p2, p3, p4):
    def orient(a, b, c):
        return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
    o1, o2 = orient(p1, p2, p3), orient(p1, p2, p4)
    o3, o4 = orient(p3, p4, p1), orient(p3, p4, p2)
    return o1 * o2 < 0 and o3 * o4 < 0


def census(poly, cap=20000, nsamp=600, beam_cap=4000):
    rows = []
    for C in poly.classes():
        N, k = poly.N, poly.k
        vert = [j for j in range(k) if (poly.idx[j] + N) % (4 * N) in C or (poly.idx[j] + 3 * N) % (4 * N) in C]
        h = sum(1 for j in range(k) if poly.idx[j] in C or (poly.idx[j] + 2 * N) % (4 * N) in C)
        sl = poly.slots(C)
        S_pred = (k - 2) * N if N % 2 else (k - 2) * N // 2 - h
        res = {}
        for (i, d) in sl:
            res[(i, d)] = poly.prong(i, d, cap)
        perp = [(r[1], r[2]) for r in res.values() if r[0] == 'perp']
        caps = [sd for sd, r in res.items() if r[0] == 'cap']
        # reversal pairing
        paired, bad_pair = 0, 0
        for (i, d), r in res.items():
            if r[0] != 'vertex':
                continue
            rev = (r[1], (r[2] + 2 * N) % (4 * N))
            rr = res.get(rev)
            if rr is None or rr[0] != 'vertex' or rr[1] != i or (rr[2] + 2 * N) % (4 * N) != d or rev == (i, d):
                bad_pair += 1
            else:
                paired += 1
        p = paired // 2
        sprime = len(caps)
        # beam check, side by side
        beam = {}
        for j in vert:
            pts = sorted(float(s) for (jj, s) in perp if jj == j)
            samples = [(t + 0.5) / nsamp for t in range(nsamp)]
            words = [poly.beam_word(j, mp.mpf(s), beam_cap) for s in samples]
            changes = [q for q in range(nsamp - 1) if words[q] != words[q + 1]]
            capped = sum(1 for w in words if w[0] == 'cap')
            # every change brackets a prong point
            unbr = [q for q in changes if not any(samples[q] < x < samples[q + 1] for x in pts)]
            # every prong point alone in its gap separates different words
            silent = 0
            for x in pts:
                q = int(x * nsamp - 0.5)
                if 0 <= q < nsamp - 1:
                    inside = [y for y in pts if samples[q] < y < samples[q + 1]]
                    if len(inside) == 1 and words[q] == words[q + 1]:
                        silent += 1
            beam[j] = {'n_pred': len(pts) + 1, 'changes': len(changes), 'unbracketed': len(unbr),
                       'silent_points': silent, 'beam_capped': capped}
        lhs = sum(b['n_pred'] for b in beam.values())
        rows.append({'class_size': len(C), 'v': len(vert), 'h': h, 'S': len(sl), 'S_pred': S_pred,
                     'f': len(perp), 'p': p, 'sprime': sprime, 'bad_pair': bad_pair,
                     'sum_n': lhs, 'identity_rhs': len(vert) + len(sl) - 2 * p - sprime,
                     'beam': {str(j): b for j, b in beam.items()},
                     'max_prong_steps': max((r[3] for r in res.values()), default=0)})
    return rows


CASES = [
    # (label, a, N, seed)  -- triangles first (validation)
    ('tri (2,3,4)/9', [2, 3, 4], 9, 0),
    ('tri (2,5,8)/15', [2, 5, 8], 15, 0),
    ('tri (1,2,3)/6', [1, 2, 3], 6, 0),
    ('quad (2,3,4,9)/9', [2, 3, 4, 9], 9, 1),
    ('quad (2,3,4,9)/9 s2', [2, 3, 4, 9], 9, 2),
    ('quad (3,3,4,4)/7', [3, 3, 4, 4], 7, 1),
    ('quad (1,2,5,6)/7', [1, 2, 5, 6], 7, 3),
    ('quad (2,4,5,7)/9', [2, 4, 5, 7], 9, 4),
    ('quad (1,3,5,7)/8', [1, 3, 5, 7], 8, 1),
    ('quad (2,3,5,6)/8', [2, 3, 5, 6], 8, 2),
    ('quad (3,4,5,6)/9', [3, 4, 5, 6], 9, 5),
    ('quad (2,2,3,5)/6', [2, 2, 3, 5], 6, 6),
    ('pent (3,4,5,6,9)/9', [3, 4, 5, 6, 9], 9, 1),
    ('pent (4,5,5,6,7)/9', [4, 5, 5, 6, 7], 9, 2),
    ('nonconvex (1,2,9,6)/9', [1, 2, 9, 6], 9, 7),   # placeholder; a_i = N is a straight angle -> skipped
    ('nonconvex (1,1,11,5)/9', [1, 1, 11, 5], 9, 8),
    ('nonconvex (2,1,13,2)/9', [2, 1, 13, 2], 9, 9),
]


def main():
    out = []
    for label, a, N, seed in CASES:
        if N in a:
            continue
        try:
            P = Poly(a, N, seed=seed)
        except (ValueError, AssertionError) as e:
            print(f'{label}: skipped ({e})')
            continue
        rows = census(P)
        for r in rows:
            ok_S = r['S'] == r['S_pred']
            ok_id = r['sum_n'] == r['identity_rhs']
            beam_ok = all(b['unbracketed'] == 0 and b['silent_points'] == 0 for b in r['beam'].values())
            T_ok = not (r['p'] == 0 and r['sprime'] > 0)
            print(f"{label:26s} |C|={r['class_size']:3d} v={r['v']} h={r['h']} S={r['S']}/{r['S_pred']} "
                  f"f={r['f']} p={r['p']} s'={r['sprime']} badpair={r['bad_pair']} Σn={r['sum_n']} "
                  f"S-ok={ok_S} beam-ok={beam_ok} T-ok={T_ok} steps≤{r['max_prong_steps']}")
            sys.stdout.flush()
        out.append({'label': label, 'a': a, 'N': N, 'seed': seed,
                    'lengths': [mp.nstr(x, 20) for x in P.L], 'rows': rows})
    json.dump(out, open('data/s910_poly_slots.json', 'w'), indent=1)


if __name__ == '__main__':
    main()
