"""[Deposit] Proposition 6.2: the length law, checked against the unfolding and the census rows. Paper: "Complete periodicity of perpendicular billiards in rational triangles of prime denominator". Research name: probes/s942_band_cos.py, copied verbatim except imports (deposit names), paths (data/ and figs/ under triangles/) and this paragraph. Run through run_all.py; the usage lines below are the research repo's. Bracketed IDs and section marks (§sNNN) refer to the project's internal ledgers.

s942: THE BAND PRONGS' LENGTHS IN CLOSED FORM — every `W_{q,m,j}` of [TRI-BAND-QM] (b) has length
    ℓ_j = R·cos((2m − 1 − 2j)A/2),   R = 2·sin B·cos(tA)/sin(β/2)
(`A = π/n`, `n = q(2m+1)`, `β = π/q`, `B = (K + t)A`; frame `V0 = 0`, `V1 = sin C`, `V2 = sin B·e^{iA}` of §s921) — so the `m` prongs of a
band point are ONE pinned cosine run at step `A`, for every `(q, m)` and every `t` at which they are realised ([POLY-BAND-COS]).

The PROOF is by hand (general_triangles.md §s942 (2)): on §s921 (b)'s development, `E_j − S` projected on `e^{iψ_j}` is
`2 sin B cos(β/2 − x)·[sin C·e^{iβ/2}/sin(β/2) − i e^{i(β/2 − tA)}]`, `x = (1 + j)A`, and `sin C = cos(β/2 − tA)` makes the bracket
`cos(tA)/sin(β/2)`. This keeper CHECKS it:
  (I) the two root-of-unity sums of the derivation, `Σ_{i<h} ω^{−2i} = ω²/(ω − 1)`, `Σ_{i<h} ω^{2i+2} = ω/(1 − ω)` (`ω = e^{iβ}`,
      `h = (q − 1)/2`), every odd `q ≤ 101`, 40 digits;
  (N) the closed form against the 40-digit unfolding of the word, every `j < m`, at `t ∈ {0.01, 0.37, t_up(0) − 10⁻⁶}` and, for
      `j ≥ 1`, at `t` between `t_up(0)` and `t_up(j)` (prong `j` still realised: `proper` checked) — `(q, m)` over 12 pairs;
  (C) the census: every `p ≥ 2` row of `data/s942_single_run.jsonl` (102) — its lengths are exactly 2× (the `triangle_poly`
      normalisation) the formula's `p` LARGEST values at the row's `(q, m, t)`.
CONTROLS (fire): (N) with the index unreversed (`cos((2j + 1)A/2)`) must fail; (C) with `R` from the other strip vertex must fail.
PRIOR ART: §s921 (b) (the development, `t_up(j)` rising, measured 488/488), §s939 (1) ([POLY-Q5-COS]: `q = 5` only, `R_m` bounded by an
interval certificate, no closed form), §s940 (4) ([POLY-FAM-COS], measured), §s942 (1) ([POLY-COS-ALL]); grep *cos(tA)/sin*, *R_m =*,
*sin B cos(tA)* → none.
DATA: `data/s942_single_run.jsonl` (§s942 (1)).
RESULT (s942): ALL OK — (I) 50 q; (N) max error ~1e-36 on 12 (q, m), controls fire; (C) 102/102 rows, scale 2, control fires.
AUDIT (audits/s942_audit_bandcos.md, 24 V / 2 G / 0 F / 6 W): both controls near-vacuous (8d) — the main tests could still fail; a rival-shape control is suggested there.
"""
import json
import math
import sys

import mpmath as mp

import _setup  # noqa: F401  (deposit: lib/ on the path, cwd = triangles/)
import s6_1_band_identities as S

mp.mp.dps = 40


def ell(q, m, j, t, rev=True):
    n, K, g0 = S.consts(q, m)
    A, beta = mp.pi / n, mp.pi / q
    B = (K + t) * A
    k = (2 * m - 1 - 2 * j) if rev else (2 * j + 1)
    return 2 * mp.sin(B) * mp.cos(t * A) / mp.sin(beta / 2) * mp.cos(k * A / 2)


def unfold_len(q, m, j, t):
    n, K, g0 = S.consts(q, m)
    A = mp.pi / n
    w = S.word_j(q, m, j)
    c = S.unfold(A, (K + t) * A, w)
    s, E = c[0][2], c[-1][2]
    d = E - s
    f = lambda z: mp.im((z - s) * mp.conj(d))
    ok = all(f(c[k][S.SIDES[int(ch)][0]]) * f(c[k][S.SIDES[int(ch)][1]]) < 0 for k, ch in enumerate(w))
    return abs(d), ok


def t_up_j(q, m, j):
    """§s921 (b): the root of (★) with A → (1 + j)A in sin(3β/2 − ·), sin(β/2 − ·)"""
    n, K, g0 = S.consts(q, m)
    A, u, x = mp.pi / n, mp.pi / (2 * q), (1 + j) * mp.pi / n
    g = lambda t: mp.sin(3 * u - x) * mp.sin((1 - t) * A) * mp.cos(u - t * A) - mp.sin(t * A) * mp.cos(u + A - t * A) * mp.sin(u - x)
    return mp.findroot(g, (mp.mpf('0.3'), mp.mpf('0.9999999')), solver='bisect')


def check_sums():
    bad = []
    for q in range(3, 102, 2):
        w, h = mp.expj(mp.pi / q), (q - 1) // 2
        s1 = sum(w ** (-2 * i) for i in range(h))
        s2 = sum(w ** (2 * i + 2) for i in range(h))
        if abs(s1 - w ** 2 / (w - 1)) > 1e-35 or abs(s2 - w / (1 - w)) > 1e-35:
            bad.append(q)
    return bad


def check_numeric():
    worst, ctrl, bad = mp.mpf(0), mp.inf, []
    for q, m in [(3, 1), (3, 4), (3, 9), (5, 3), (5, 6), (7, 2), (9, 3), (11, 2), (13, 3), (15, 2), (21, 1), (25, 2)]:
        tu0 = S.t_up(q, m)
        for j in range(m):
            ts = [mp.mpf('0.01'), mp.mpf('0.37'), tu0 - mp.mpf('1e-6')]
            if j >= 1:
                ts.append((tu0 + t_up_j(q, m, j)) / 2)
            for t in ts:
                L, ok = unfold_len(q, m, j, t)
                if not ok:
                    bad.append((q, m, j, float(t)))
                worst = max(worst, abs(L - ell(q, m, j, t)))
                if m >= 2:
                    ctrl = min(ctrl, max(abs(unfold_len(q, m, jj, t)[0] - ell(q, m, jj, t, rev=False)) for jj in range(m)))
    return worst, ctrl, bad


def strip(abc):
    a, b, c = abc
    N = sum(abc)
    n = N // a
    for q in range(3, n + 1, 2):
        if n % q:
            continue
        m = (n // q - 1) // 2
        if m < 1 or q * (2 * m + 1) != n:
            continue
        K = (n - 2 * m - 3) // 2
        for Bv in (b, c):
            if K * a < Bv < (K + 1) * a:
                return q, m, K, Bv, (c if Bv == b else b)
    return None


def check_census():
    rows = [json.loads(l) for l in open('data/s942_single_run.jsonl')]
    bad, ctrl_pass = [], 0
    for r in rows:
        abc = tuple(r['abc'])
        st = strip(abc)
        if st is None:
            bad.append((abc, 'no strip'))
            continue
        q, m, K, Bv, other = st
        a, N = abc[0], sum(abc)
        t = mp.mpf(Bv) / a - K
        f = sorted(2 * ell(q, m, j, t) for j in range(m))
        L = sorted(r['lengths'])
        p = len(L)
        if p > m or max(abs(L[k] - f[m - p + k]) for k in range(p)) > 1e-8:
            bad.append((abc, p, m))
        # control: R built on the other non-V0 vertex as B
        Bo = mp.pi * other / N
        g = sorted(2 * 2 * mp.sin(Bo) * mp.cos(t * mp.pi / (N // a)) / mp.sin(mp.pi / (2 * q)) * mp.cos((2 * m - 1 - 2 * j) * mp.pi / (2 * (N // a))) for j in range(m))
        if max(abs(L[k] - g[m - p + k]) for k in range(p)) < 1e-8:
            ctrl_pass += 1
    return len(rows), bad, ctrl_pass


if __name__ == '__main__':
    bs = check_sums()
    print('(I) root-of-unity sums, odd q ≤ 101: bad', bs)
    worst, ctrl, bad = check_numeric()
    print('(N) max |unfold − closed form| =', mp.nstr(worst, 3), '| not proper:', bad, '| control (unreversed) min max-error =', mp.nstr(ctrl, 3))
    n, bad_c, cp = check_census()
    print(f'(C) census rows {n}: mismatches {bad_c} | control passes (should be 0): {cp}')
    ok = not bs and worst < 1e-30 and ctrl > 1e-6 and not bad and not bad_c and cp == 0
    print('ALL OK' if ok else 'FAIL')
