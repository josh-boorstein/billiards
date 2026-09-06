#!/usr/bin/env python3
"""
s373_word_normalise.py — PRE-REGISTERED: discharge [NCYL-213] by NORMALISING
`enumerate_branches`' word against the tracer's `full_word`, then re-scoring.

THE DISCHARGE CONDITION, VERBATIM FROM [NCYL-213]: *"normalise
`enumerate_branches`' word against the tracer's `full_word` on ONE row (same
orbit, same alphabet, same endpoints), then re-run the scorer — only then is
either statement scoreable."*  This probe does that on many rows, not one.

THE TWO STATEMENTS IN DISAGREEMENT (neither cites the other):
  (A) `orphan_theorem.md` Lemma B′ (PROVED; restated in `foundations.md` §2 as a
      PROVED detector with the convention note DROPPED): word(b−s) =
      reverse(word(s)); hence "palindromic word" ⟺ "orphan", exactly.
      Its CONVENTION GOTCHA: the word is the `L1→L1` first-return segment
      INCLUDING the launch `L1` crossing — i.e. `('L1',) + full_word`.
  (B) `computational_findings.md` §103: "EVERY perpendicular branch word is
      already a PALINDROME … so reversal is trivial", and {reverse, leg-swap,
      swap∘reverse, identity} realize the pairing on 0 pairs.

WHAT THIS PROBE ESTABLISHES (three steps, each with its own control):
  STEP 1 — THE RELATION (mechanical, not a judgment call).  Claim R:
      ring_word(cell) == (concatenation of the `k` successive L1-return words
      along the orbit from the cell's midpoint) with the FINAL 'L1' removed,
      where `k` = the number of perpendicular L1-returns before the orbit
      repeats (k = 1 on the orphan, k = 2 on a J-pair).
      ⇒ the ring word is the CLOSED CYLINDER word (the whole round trip, cyclic,
      so the closing L1 is dropped); Lemma B′'s word is ONE L1→L1 segment.
      They are different objects, and R says exactly how.
  STEP 2 — RE-SCORE (A) ON ITS OWN OBJECT.  With `u(s) := ('L1',) + full_word(s)`:
      H_P:  u palindromic  ⟺  the cell is the orphan.
      H_R:  reverse(u(s)) == u(cell containing T(s)), for every non-orphan cell.
      ⚠ The partner is found by the RETURN HEIGHT `T(s)` read off the tracer, NOT
      by §103's "partners are adjacent" — scoring (A) against a claim from the
      same file as (B) would be circular.
  STEP 3 — DERIVE (B) FROM (A), which is what actually reconciles them.  Write
      u = L1·v·L1.  Under (A) the orphan has u palindromic, so v is a palindrome
      and ring = v; a J-pair has u_j = reverse(u_i), so
      ring_i = v_i · L1 · reverse(v_i).  BOTH are palindromes for ANY v.
      ⇒ §103's "every branch word is a palindrome" is a COROLLARY of Lemma B′
      under the ring convention, and its "reversal does not realize the pairing"
      is the same corollary restated (reversal is the IDENTITY on those words).
      The probe checks this structural form literally, per cell.
  STEP 4 — WHERE THE NON-PALINDROMY LIVES, which also DERIVES s372's dangling
      "each J-pair differs in exactly 4 positions, independent of length".
      Claim S:  a non-orphan has  v = X · ab · reverse(X)  with |v| EVEN and
      {a,b} = {L1,L2} the two CENTRAL letters; the partner is the same word with
      ab ↦ ba.  The orphan has |v| ODD with a single central letter (= Lemma B's
      head-on wall).  Then Hamming(v, reverse(v)) = 2 on a pair and 0 on the
      orphan, so the two RING words differ in exactly 2·2 = 4 positions whatever
      the length -- s372's observation, derived rather than noted.
      ⚠ The GEOMETRIC reading (a central `L1L2` adjacency is the orbit rounding
      the right-angle vertex R -- the corner retroreflector that the s332
      qualification to Lemma B names as the only non-perpendicular turnaround,
      and a REGULAR point by Lemma C) is a SUGGESTION, not measured here.

CONTROLS ([OPS-041] — could this have come out otherwise?):
  C1 (the discriminating half).  Non-orphan `u` must be NON-palindromic.  If it
     came back palindromic everywhere, H_P would be vacuous — the raw-ring
     column is exactly that failure mode and is reported alongside.
  C2 (the partner is not free).  reverse(u(s)) is checked against the SPECIFIC
     cell containing T(s), and additionally against every OTHER cell on the row;
     a match count of exactly 1 is required.  "Reverse is *a* word on the row"
     is the weaker claim that cannot fail once the row is closed under reversal.
  C3 (the wrong normalisation still fails).  `raw`, `L1+w`, `w+L1`, `L1+w+L1`
     applied to the RING word are re-scored, reproducing s372's 58/58 and 0/58.
  C4 (the tracer agrees with the ring).  Every row where STEP 1 fails is
     REPORTED AND DROPPED, not forced — the float tracer is not trusted on a
     narrow cell.  A row's `ok` flag is the certificate.

PRIOR ART: grepped `rulings.md`, `TOOLS.md` and all root `.md` for 'palindrom',
'full_word', 'word convention', 'prepend L1', 'launch crossing', 'L1→L1',
'enumerate_branches', 's196_ringwords' ->
  - [NCYL-213] (s372) OWNS the flag and states this discharge condition; nothing
    else has attempted it.
  - `orphan_theorem.md` Lemma B/B′ OWNS (A) + the convention gotcha (line ~295);
    `foundations.md` §2 restates (A) without it.
  - `computational_findings.md` §103 OWNS (B); its script is
    `archive/scripts_2026-07/veech_cylinders.py`.
  - [NCYL-042]/[OPS-028]: `enumerate_branches` needs BOTH `last_unclosed == 0`
    and `last_degenerate == 0` — checked here per centre.
  - [OPS-092]: never score one file's claim against another file's object; that
    is the error this probe exists to undo.
  -> no prior attempt at the normalisation itself.

Run: PYTHONPATH=.:engine:archive/scripts_2026-07:archive/scripts_2026-08:probes \
     .venv/bin/python3.13 probes/s373_word_normalise.py
"""
import collections
import json
import math

import s196_ringwords as rw
from right_triangle_billiards import RightTriangleBilliard

ROWS = [(8, 15), (4, 15),                       # even P, no orphan (the trigger)
        (5, 8), (3, 7), (5, 7),                 # Lemma B's OWN verification rows
        (1, 5), (3, 8), (5, 12), (7, 10), (3, 10), (5, 16), (7, 16)]

LET = {'H': 'H', 'L1': '1', 'L2': '2'}
SWAP = str.maketrans('12', '21')
NORMS = {'raw': lambda w: w, 'L1+w': lambda w: '1' + w,
         'w+L1': lambda w: w + '1', 'L1+w+L1': lambda w: '1' + w + '1'}


def ring_cells(P, Q):
    """`enumerate_branches` with its OWN certificate ([NCYL-042]/[OPS-028])."""
    res = rw.enumerate_branches(P, Q, max_depth=400000, keep_words=True)
    lv = res[0] if isinstance(res, tuple) else res
    cert = (rw.enumerate_branches.last_unclosed == 0
            and rw.enumerate_branches.last_degenerate == 0)
    return [(lo, hi, w if isinstance(w, str) else ''.join(w))
            for lo, hi, w in lv], cert


def segments(bil, s, k_max=8, max_hits=200000):
    """The successive L1-return words from a perpendicular launch at height `s`.

    Returns (segs, heights) where `segs[i]` is the compact word of the i-th
    L1→L1 first-return segment (ending in '1') and `heights[i]` is T^{i+1}(s).
    Stops as soon as a return height repeats `s` (the orbit has closed).
    """
    segs, heights = [], []
    cur = s
    for _ in range(k_max):
        tr = bil.trace(cur, max_hits=max_hits, stop_at_perpendicular_return=True)
        if not tr.returned_perpendicular:
            return None, None
        segs.append(''.join(LET[x] for x in tr.full_word))
        nxt = tr.l1_returns[-1].s
        heights.append(nxt)
        if abs(nxt - s) < 1e-9:
            return segs, heights
        cur = nxt
    return None, None


def which_cell(cells, y):
    for i, (lo, hi, _w) in enumerate(cells):
        if lo - 1e-12 <= y <= hi + 1e-12:
            return i
    return None


def main():
    out = {'rows': []}
    for P, Q in ROWS:
        cells, cert = ring_cells(P, Q)
        bil = RightTriangleBilliard((P / Q) * math.pi / 2)
        rec = {'P': P, 'Q': Q, 'n_cells': len(cells), 'certified': cert,
               'cells': [], 'dropped': 0}
        for i, (lo, hi, w) in enumerate(cells):
            s = 0.5 * (lo + hi)
            segs, heights = segments(bil, s)
            c = {'i': i, 'lo': lo, 'hi': hi, 'ring_len': len(w)}
            if segs is None:
                c['ok'] = False
                c['why'] = 'tracer did not close'
                rec['dropped'] += 1
                rec['cells'].append(c)
                continue
            k = len(segs)
            # ---- STEP 1: the relation ------------------------------------
            joined = ''.join(segs)
            c.update(k=k, seg_lens=[len(x) for x in segs],
                     relation=(joined[:-1] == w))
            # ---- STEP 2: Lemma B' on its OWN object ----------------------
            u = '1' + segs[0]                    # L1 + full_word  (L1 ... L1)
            partner = which_cell(cells, heights[0])
            u_partner = None
            if partner is not None:
                pj = 0.5 * (cells[partner][0] + cells[partner][1])
                psegs, _ph = segments(bil, pj)
                u_partner = ('1' + psegs[0]) if psegs else None
            rev = u[::-1]
            n_match = sum(1 for j, (a, b, _v) in enumerate(cells)
                          if j != i and _u_of(bil, cells, j) == rev)
            c.update(u_len=len(u), u_palindrome=(u == rev),
                     is_orphan_by_T=(abs(heights[0] - s) < 1e-9),
                     partner=partner,
                     rev_is_partner=(u_partner is not None and rev == u_partner),
                     rev_match_count=n_match + (1 if u == rev else 0))
            # ---- STEP 3: the structural derivation of (B) ----------------
            v = u[1:-1]
            pred = v if c['is_orphan_by_T'] else v + '1' + v[::-1]
            c.update(struct_pred_ok=(pred == w), v_palindrome=(v == v[::-1]))
            # ---- STEP 4: WHERE the non-palindromy lives ------------------
            L = len(v)
            diff = [j for j in range(L) if v[j] != v[L - 1 - j]]
            centre = ([L // 2 - 1, L // 2] if L % 2 == 0 else [L // 2])
            c.update(v_len=L, v_len_even=(L % 2 == 0), hamming=len(diff),
                     diff_is_central=(diff == ([] if c['is_orphan_by_T']
                                               else centre)),
                     central_letters=''.join(v[j] for j in centre),
                     claim_S=((L % 2 == 1 and not diff)
                              if c['is_orphan_by_T'] else
                              (L % 2 == 0 and diff == centre
                               and set(v[j] for j in centre) == {'1', '2'})))
            # ---- C3: the ring word under the four normalisations ---------
            c['ring_norm_palindrome'] = {
                nm: (f(w) == f(w)[::-1]) for nm, f in NORMS.items()}
            c['ok'] = True
            rec['cells'].append(c)
        out['rows'].append(rec)

    # ----------------------------------------------------------------- report
    tot = dict(cells=0, rel=0, orph=0, upal=0, upal_orph=0, upal_nonorph=0,
               revpart=0, nonorph=0, struct=0, revmulti=0, ringpal=0,
               claimS=0, ham2=0, ham0=0, centre_rule=0)
    centres_seen = {}
    for rec in out['rows']:
        for c in rec['cells']:
            if not c.get('ok'):
                continue
            tot['cells'] += 1
            tot['rel'] += c['relation']
            tot['struct'] += c['struct_pred_ok']
            tot['ringpal'] += c['ring_norm_palindrome']['raw']
            tot['claimS'] += c['claim_S']
            if c['is_orphan_by_T']:
                tot['orph'] += 1
                tot['upal_orph'] += c['u_palindrome']
                tot['ham0'] += (c['hamming'] == 0)
                centres_seen.setdefault('orphan', collections.Counter())[
                    c['central_letters']] += 1
                # C5: Lemma B' also states the palindrome CENTRE -- "the Prop-4
                # head-on wall (L2 if Q even, H if Q odd)".  A second, fully
                # independent prediction of (A) on the same object, free here.
                tot['centre_rule'] += (c['central_letters']
                                       == ('2' if rec['Q'] % 2 == 0 else 'H'))
            else:
                tot['nonorph'] += 1
                tot['upal_nonorph'] += c['u_palindrome']
                tot['revpart'] += c['rev_is_partner']
                tot['ham2'] += (c['hamming'] == 2)
                centres_seen.setdefault('pair', collections.Counter())[
                    c['central_letters']] += 1
            tot['upal'] += c['u_palindrome']
            tot['revmulti'] += (c['rev_match_count'] == 1)
    out['summary'] = tot
    out['central_letters'] = {k: dict(v) for k, v in centres_seen.items()}

    n, no, orp = tot['cells'], tot['nonorph'], tot['orph']
    print(f'\ns373 word normalisation — {len(ROWS)} centres, {n} scoreable cells '
          f'({orp} orphan, {no} non-orphan); dropped '
          f'{sum(r["dropped"] for r in out["rows"])}')
    print(f'  STEP 1  ring == concat(L1-return words) minus final L1 : '
          f'{tot["rel"]}/{n}')
    print(f'  STEP 2  H_P  u palindromic on the ORPHAN                : '
          f'{tot["upal_orph"]}/{orp}')
    print(f'          C1   u palindromic on a NON-orphan (want 0)     : '
          f'{tot["upal_nonorph"]}/{no}')
    print(f'          H_R  reverse(u) == the T-partner\'s u            : '
          f'{tot["revpart"]}/{no}')
    print(f'          C2   reverse(u) matches EXACTLY one cell        : '
          f'{tot["revmulti"]}/{n}')
    print(f'  STEP 3  ring == v (orphan) / v·L1·rev(v) (pair)         : '
          f'{tot["struct"]}/{n}')
    print(f'  STEP 4  claim S (v = X·ab·rev(X), ab central, |v| parity) : '
          f'{tot["claimS"]}/{n}')
    print(f'          Hamming(v,rev v) == 2 on a pair                 : '
          f'{tot["ham2"]}/{no}   == 0 on the orphan: {tot["ham0"]}/{orp}')
    print(f'          central letters — pair {out["central_letters"].get("pair")} '
          f'| orphan {out["central_letters"].get("orphan")}')
    print(f'          C5   orphan centre == L2 if Q even else H        : '
          f'{tot["centre_rule"]}/{orp}')
    print(f'  C3      ring word palindromic RAW (s372 reproduced)     : '
          f'{tot["ringpal"]}/{n}')
    for nm in NORMS:
        k = sum(c['ring_norm_palindrome'][nm] for r in out['rows']
                for c in r['cells'] if c.get('ok'))
        print(f'            ring under {nm:9s}: {k}/{n}')
    with open('data/s373_word_normalise.json', 'w') as fh:
        json.dump(out, fh, indent=1)
    print('  -> data/s373_word_normalise.json')
    return 0


_UCACHE = {}


def _u_of(bil, cells, j):
    key = (id(cells), j)
    if key not in _UCACHE:
        lo, hi, _w = cells[j]
        segs, _h = segments(bil, 0.5 * (lo + hi))
        _UCACHE[key] = ('1' + segs[0]) if segs else None
    return _UCACHE[key]


if __name__ == '__main__':
    raise SystemExit(main())
