"""[Deposit] Appendix B, step 1: every slot walk of the two families, m = 4..30, classified into word families. Paper: "Complete periodicity of perpendicular billiards in rational triangles of prime denominator". Research name: probes/s947_q3_families.py, copied verbatim except imports (deposit names), paths (data/ and figs/ under triangles/) and this paragraph. Run through run_all.py; the usage lines below are the research repo's. Bracketed IDs and section marks (§sNNN) refer to the project's internal ledgers.

s947 `s947-a` survey: EVERY perpendicular-class slot walk of the note's Prop 8.2 families (q = 3), classified into word families.

Families F1 = (3, 6m+2, 12m+4), F2 = (3, 6m+1, 12m+5), N = 18m + 9 (polygon note Prop 8.2 / Rem 8.3; [TRI-SWAP-FAM3]).
`walk` = s939's (`s929_asym_witness.walk`, cap 40m + 100). MEASURED (m = 4 … 30, both): exactly 2m slots end at a vertex
(the m band prongs and reverses) — p = m — all others perpendicularly, none capped.
CLASSIFICATION = s939's (strands = slot residue mod 3 per vertex, t = position; segments of constant run-vector increment) ⊕ the
FIRST-LETTER pattern must agree (s939's omits it and merges different words — s947). From m = 7: 30 (F1) / 26 (F2) segments,
every segment's runs EXACTLY affine in (m, t) on m = 7 … 19, and every word on HELD-OUT m = 20 … 30 matches its fit (no word
transitions in m — unlike q = 5). Boundaries fitted as ⌊(am + b)/q⌋ on m = 7 … 19: exactly 4 per family fail held-out m = 20 … 30
(transcendental; `probes/s947_q3_graze.py` proves no lattice point lies on them).
PRIOR ART: §s939 (2) (`s939_q5_families.py`, the q = 5 survey — 30/32 segments, word transitions at m = 20), [TRI-FAM-CP] (p = m
measured to N = 303), polygon note Rem 8.3; `rulings.py --grep` *word family*, *slot walk* → none on this object. DATA:
`data_index.py --grep walk` → the q = 5 walks only.
Usage: `probes/s947_q3_families.py F1|F2` → `data/s947_q3_walks_{fam}.json`, `data/s947_q3_families_{fam}.json` (~8 min each)."""
import collections, json, sys
import _setup  # noqa: F401  (deposit: lib/ on the path, cwd = triangles/)
import mpmath as mp
import sympy as sp
import word_classify as S
from poly_slots import Poly
import slot_walk as W
mp.mp.dps = 40
ABC = {'F1': lambda m: (3, 6 * m + 2, 12 * m + 4), 'F2': lambda m: (3, 6 * m + 1, 12 * m + 5)}
MS, MS_FIT = range(4, 31), range(7, 20)


def segments(rows):
    segs = []
    for (v, res), st in S.strands(rows).items():
        cur = []
        for t, r in enumerate(st):
            if cur:
                prev = cur[-1][1]
                ok = len(r[3]) == len(prev[3]) and r[1] == prev[1] and r[4] == prev[4]
                if ok and len(cur) >= 2:
                    ok = [a - b for a, b in zip(r[3], prev[3])] == [a - b for a, b in zip(prev[3], cur[-2][1][3])]
                if not ok:
                    segs.append((v, res, cur)); cur = []
            cur.append((t, r))
        if cur:
            segs.append((v, res, cur))
    return segs


S.segments = segments


def walks(fam, m):
    abc = ABC[fam](m); N = sum(abc)
    P = Poly(list(abc), N, lengths=[mp.sin(mp.pi * x / N) for x in (abc[2], abc[0], abc[1])])
    rows = []
    for k in sorted(P.slots(P.classes()[0])):
        r = W.walk(P, *k, mp.mpf(10) ** -20, cap=40 * m + 100)
        parts = ''.join(map(str, r[5])).split('1')
        rows.append([list(k), r[0], str(r[1]), [len(x) for x in parts], ''.join(x[:1] or '-' for x in parts)])
    return rows


def main(fam):
    W_ = {m: walks(fam, m) for m in MS}
    json.dump({m: {'kinds': dict(collections.Counter(r[1] for r in W_[m])), 'rows': W_[m]} for m in MS},
              open(f'data/s947_q3_walks_{fam}.json', 'w'))
    assert all(sum(r[1] == 'vertex' for r in W_[m]) == 2 * m and all(r[1] != 'cap' for r in W_[m]) for m in MS)
    S.MS_WALK, S.MS_FIT = list(MS), MS_FIT
    counts, allaff, table = S.classify(fam, W_)
    HOLD = [m for m in MS if m >= 20]
    trans, hold = [], []
    for t in table:
        for m in HOLD:
            seg = segments(W_[m])
            if len(seg) != len(table) or len(seg[t['i']][2][0][1][3]) != t['nruns'] or any(
                    sp.Rational(t['runs'][c][0]) * m + sp.Rational(t['runs'][c][1]) * tt + sp.Rational(t['runs'][c][2]) != x
                    for tt, r in seg[t['i']][2] for c, x in enumerate(r[3]) if c in t['runs']):
                trans.append([t['i'], m]); break
        for name, j in (('t_lo', 0), ('t_hi', -1)):
            fit = t[name]
            if not (fit is not None and all((fit[0] * m + fit[1]) // fit[2] == segments(W_[m])[t['i']][2][j][0] for m in HOLD)):
                hold.append([t['i'], name, str(fit)])
    json.dump({'fam': fam, 'counts': counts, 'all_affine': allaff, 'transitions': trans, 'holdout_fail': hold, 'families': table},
              open(f'data/s947_q3_families_{fam}.json', 'w'), indent=1, default=str)
    print(fam, '| p = m, no caps (m = 4..30) | segments m≥7', sorted({counts[m] for m in MS if m >= 7}), '| all affine', allaff,
          '| transitions', trans, '| transcendental boundaries', len(hold) // 2)
    assert allaff and not trans


if __name__ == '__main__':
    main(sys.argv[1])
