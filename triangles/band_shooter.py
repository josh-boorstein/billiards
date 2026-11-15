"""[Deposit] helper: the 30-digit developer of a word (used by the Theorem D certificate as an independent instrument). Paper: "Complete periodicity of perpendicular billiards in rational triangles of prime denominator". Research name: probes/s922_band_shooter.py, copied verbatim except imports (deposit names), paths (data/ and figs/ under triangles/) and this paragraph. Run through run_all.py; the usage lines below are the research repo's. Bracketed IDs and section marks (§sNNN) refer to the project's internal ledgers.

s922 (queue `s921-b`, `s921-a`): WORD-FREE search for swapped loop bands held at V0, continuous in B.

A band prong held at V0 (`s911_hband` / [NCYL-710]) is a word w whose unfolded closing direction is psi = pi/2 + s*B + c*A0
(A0 = pi/n, s, c integers) identically on an open B-interval, at fixed A = k*pi/n. A census (`probes/s912_census.py`) sees B
only at multiples of pi/N, so it misses every band narrower than that spacing (s922: the stored rows hold one slope -1
word at n = 15; this search finds two, and a slope -5 band). Here: for B on a grid of step A0/per in the canonical frame
(A < B < C, V2 = the largest angle, as in the census labelling a < b < c), and every class direction psi = pi/2 + s*B +
c*A0 entering the triangle at V2, trace the billiard path (float, vectorised) for <= L bounces; record paths ending at a
vertex; a word (final hit side stripped) seen at >= 2 grid points, not a palindrome (a palindrome retraces, Prop S), is a
swapped band CANDIDATE, then re-checked at 40 digits at its interval midpoint (every crossing proper; d psi / d B by finite
differences, integer). (Fixed before any scored run: a fwd/rev pair at one isolated B had counted as 2 points.) Frame V0 = 0, V1 = sin C, V2 = sin B e^{iA}; sides 0 = V0V1, 1 = V1V2, 2 = V2V0 (`wqm.develop`'s).

PRE-REGISTERED (s922, before any n other than 9 was searched with s = +1).
  (S+) For every odd n in 11..39, A = pi/n, s = +1, start V2: the swapped V2 -> V2 band candidates found are EXACTLY the
       `W_{q,m,j}` of [TRI-BAND-QM] with q(2m+1) = n, 0 <= j < m (each found, as itself or reversed), and nothing else --
       in particular NONE at prime n (11, 13, 17, 19, 23, 29, 31, 37).
  Could fail: any extra candidate that passes the 40-digit check (a second slope +1 family; that refutes `s921-b` as
  stated), or a predicted W not found (then the INSTRUMENT is at fault -- the theorem proves them -- and the row is
  reported, not scored).
  Scope stated in advance: bands narrower than 2 grid steps (2 A0/per) or longer than L are invisible; start vertex V2,
  held vertex V0 at A = pi/n only (A = k pi/n, k >= 2, and starts at V0, V1 are descriptive).
  Control: n = 9 (seen, s921 data): s = +1 finds only W_{3,1}; s = -1 finds the stored (5,18,22)/45 word.
  Descriptive: s = -1, -3, -5 at the same n (queue `s921-a`: the slope -1 family's m-dependence).
  PRIOR ART: `rulings.py --grep` band / swapped prong / W_{q -> [NCYL-710] (every loop seen a band, held at V0, slope
  +-1), [NCYL-712]; §s921 (a)-(c); §s915 (2) (slope -1/-5 rows, band ends); nothing word-free.
  DATA: `data_index.py --grep census` -> other strands; stored band rows `data/s911_hband.json`,
  `data/s912_census_53_99.jsonl`; output -> data/s922_band_shooter.json.
Run from the repo root:
  PYTHONPATH=.:engine:archive/scripts_2026-07:archive/scripts_2026-08:probes:data/scratch/s921 .venv/bin/python3.13 \
      probes/s922_band_shooter.py <n_lo> <n_hi> <slopes, e.g. 1 or -1,-5> [L] [per]
"""
import _setup  # noqa: F401  (deposit: lib/ on the path, cwd = triangles/)
import json
import os
import sys
import time
from collections import defaultdict

import mpmath as mp
import numpy as np

from s6_1_band_identities import word_j

SIDES = [(0, 1), (1, 2), (2, 0)]
mp.mp.dps = 40


def trace(A, Bs, psis, L, eps=1e-9):
    """Vectorised float billiard from V2 (one triangle (A, B_r) and direction psi_r per ray). Returns the end vertex
    (-1 = not within L), the step count, and the side sequence (int8, L x R; the last entry is the side holding the
    vertex hit)."""
    R = len(Bs)
    C = np.pi - A - Bs
    V = np.stack([np.zeros(R, complex), np.sin(C) + 0j, np.sin(Bs) * np.exp(1j * A)], axis=1)
    p = V[:, 2].copy()
    d = np.exp(1j * psis)
    last = np.full(R, -1)
    alive = np.ones(R, bool)
    end = np.full(R, -1)
    steps = np.zeros(R, int)
    words = np.full((L, R), -1, np.int8)
    idx = np.arange(R)
    SA = np.array([s[0] for s in SIDES]); SB = np.array([s[1] for s in SIDES])
    for k in range(L):
        ii = idx[alive]
        if len(ii) == 0:
            break
        pp, dd, VV = p[ii], d[ii], V[ii]
        best_t = np.full(len(ii), np.inf); best_s = np.full(len(ii), -1); best_u = np.zeros(len(ii))
        for s, (a, b) in enumerate(SIDES):
            P0 = VV[:, a]; E = VV[:, b] - P0
            den = dd.real * E.imag - dd.imag * E.real
            w = P0 - pp
            with np.errstate(divide='ignore', invalid='ignore'):
                t = (w.real * E.imag - w.imag * E.real) / den
                u = (w.real * dd.imag - w.imag * dd.real) / den
            ok = (t > 1e-12) & (u > -1e-7) & (u < 1 + 1e-7)
            ok &= (s == 0) if k == 0 else (last[ii] != s)
            better = ok & (t < best_t)
            best_t[better] = t[better]; best_s[better] = s; best_u[better] = u[better]
        lost = best_s < 0                      # numerical failure: no side found -> drop the ray (counted)
        bs = np.where(lost, 0, best_s)
        hit = pp + np.where(lost, 0, best_t) * dd
        words[k, ii] = best_s
        steps[ii] = k + 1
        rr = np.arange(len(ii))
        lenE = np.abs(VV[rr, SB[bs]] - VV[rr, SA[bs]])
        d0 = best_u * lenE; d1 = (1 - best_u) * lenE
        atv = (d0 < eps) | (d1 < eps) | lost
        vert = np.where(d0 < eps, SA[bs], SB[bs])
        vert = np.where(lost, -2, vert)
        end[ii[atv]] = vert[atv]
        alive[ii[atv]] = False
        nv = ~atv
        jj = ii[nv]; s_ = bs[nv]
        E = V[jj, SB[s_]] - V[jj, SA[s_]]
        e = E / np.abs(E)
        d[jj] = e * e * np.conj(d[jj])
        p[jj] = hit[nv]
        last[jj] = s_
    return end, steps, words


def develop(A, B, word):
    C = mp.pi - A - B
    V = [mp.mpc(0), mp.mpc(mp.sin(C)), mp.sin(B) * mp.expj(A)]
    S = V[2]; cps = [V]
    for ch in word:
        i, j = SIDES[int(ch)]; k = 3 - i - j
        u = (V[j] - V[i]) / abs(V[j] - V[i])
        V = list(V); V[k] = V[i] + u * u * mp.conj(V[k] - V[i]); cps.append(V)
    return S, cps


def verify(A, B, w):
    """40 digits: every crossing proper, ends at V2; returns (proper, d psi/dB)."""
    S, cps = develop(A, B, w)
    d = cps[-1][2] - S
    f = lambda z: mp.im((z - S) * mp.conj(d))
    ok = all(f(cps[k][SIDES[int(ch)][0]]) * f(cps[k][SIDES[int(ch)][1]]) < 0 for k, ch in enumerate(w))
    h = mp.mpf(10) ** -10
    ps = [mp.arg(develop(A, B + e, w)[1][-1][2] - develop(A, B + e, w)[0]) for e in (-h, h)]
    return ok, float((ps[1] - ps[0]) / (2 * h))


def family(nmax=401):
    F = {}
    for q in range(3, nmax, 2):
        for m in range(1, (nmax // q - 1) // 2 + 1):
            for j in range(m):
                w = word_j(q, m, j); F[w] = (q, m, j); F[w[::-1]] = (q, m, j)
    return F


def predicted(n):
    return sorted((q, (n // q - 1) // 2, j) for q in range(3, n + 1, 2) if n % q == 0 and n // q >= 3
                  for j in range((n // q - 1) // 2))


def search(n, s, L=2000, per=400, chunk=60000, blo=None, bhi=None):
    A0 = np.pi / n; A = A0
    Bs = np.arange(A, (np.pi - A) / 2, A0 / per)[1:]
    if blo is not None:                                  # optional window in units of A0 (descriptive runs only)
        Bs = Bs[(Bs >= blo * A0) & (Bs <= bhi * A0)]
    rays = []
    for B in Bs:
        for c in range(2 * n):
            psi = (np.pi / 2 + s * B + c * A0) % (2 * np.pi)
            x = (psi - (A - np.pi)) % (2 * np.pi)       # V2's interior arc: (A - pi, -B), width C
            if 1e-9 < x < (np.pi - A - B) - 1e-9:
                rays.append((B, psi, c))
    rays = np.array(rays)
    hits = defaultdict(list); unfinished = lost = 0
    for i0 in range(0, len(rays), chunk):
        R = rays[i0:i0 + chunk]
        end, steps, words = trace(A, R[:, 0], R[:, 1], L)
        unfinished += int((end == -1).sum()); lost += int((end == -2).sum())
        for r in np.nonzero(end == 2)[0]:
            w = ''.join(str(x) for x in words[:steps[r] - 1, r])
            if w and w != w[::-1]:
                hits[min(w, w[::-1])].append(R[r, 0] / A0)
    F = family()
    cands = []
    for w, bs in hits.items():
        bs = sorted(set(round(b, 9) for b in bs))     # a fwd/rev pair at ONE isolated B is not a band
        if len(bs) < 2:
            continue
        lo, hi = bs[0], bs[-1]
        ok, sl = verify(mp.pi / n, mp.mpf((lo + hi) / 2) * mp.pi / n, w)
        cands.append({'w': w, 'len': len(w), 'lo': lo, 'hi': hi, 'pts': len(bs), 'proper_mid': ok,
                      'slope': round(sl, 6), 'W': F.get(w)})
    cands.sort(key=lambda x: x['lo'])
    return {'n': n, 's': s, 'L': L, 'per': per, 'rays': len(rays), 'unfinished': unfinished, 'lost': lost, 'cands': cands}


if __name__ == '__main__':
    lo, hi = int(sys.argv[1]), int(sys.argv[2])
    slopes = [int(x) for x in sys.argv[3].split(',')]
    L = int(sys.argv[4]) if len(sys.argv) > 4 else 2000
    per = int(sys.argv[5]) if len(sys.argv) > 5 else 400
    out = os.environ.get('BAND_OUT', 'data/s922_band_shooter.json')   # s923: per-worker file for parallel runs (merge after)
    store = json.load(open(out)) if os.path.exists(out) else {}
    for n in range(lo, hi + 1, 2):
        for s in slopes:
            t0 = time.time()
            r = search(n, s, L, per)
            r['sec'] = round(time.time() - t0, 1)
            store[f'{n},{s},{L},{per}'] = r
            json.dump(store, open(out, 'w'), indent=1)
            real = [c for c in r['cands'] if c['proper_mid'] and abs(c['slope'] - round(c['slope'])) < 1e-4]
            found = sorted({c['W'] for c in real if c['W'] and abs(c['slope'] - 1) < 1e-4})
            extra = [c for c in real if not (c['W'] and abs(c['slope'] - s) < 1e-4)] if s == 1 else []
            line = f"n={n} s={s:+d} rays={r['rays']} unfinished={r['unfinished']} lost={r['lost']} {r['sec']}s cands={len(real)}"
            if s == 1:
                pred = predicted(n)
                line += f" | predicted {len(pred)} found {len(found)} missing {sorted(set(pred) - set(found))} EXTRA {len(extra)}"
            print(line, flush=True)
            for c in real:
                print(f"    [{c['lo']:.4f},{c['hi']:.4f}] pts {c['pts']} len {c['len']} slope {c['slope']:+.4f} W {c['W']}", flush=True)
