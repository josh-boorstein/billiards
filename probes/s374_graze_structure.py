"""s374 PRE-REGISTERED -- the STRUCTURAL predictions of a proof of [NCYL-216]'s
forward direction, and the sharp risk set for its converse.

WHAT IS BEING ASKED.  s373 measured [NCYL-216] (an interior boundary is an `R`-graze iff
it separates a `J`-pair) and the s373 handoff item (f) asks whether it is PROVABLE.
Working the geometry gives a proof of the FORWARD direction and reduces the converse to
two named residuals.  This probe tests the three things that derivation PREDICTS but
which nobody has scored, and locates the one place the converse could actually fail.

  H1 (a CONSEQUENCE of the forward proof, so a genuine falsifier -- if this fails the
      proof is wrong).  NO TWO CONSECUTIVE interior boundaries are both `R`-grazes.
      Reason: an `R`-graze at `y` forces the shared return constant `b = 2y` on BOTH
      flanking cells; two `R`-grazes `y1 < y2` flanking one cell would force
      `2y1 = b = 2y2`.  ⇒ `R`-boundaries form a MATCHING: the cells tile as dominoes
      (a `J`-pair) plus singletons.

  H2 (the converse, restated as a tiling statement).  The ONLY singleton is the orphan.
      Equivalent, given H1, to `#R = #J`-pairs and to [NCYL-186]'s unproved detail
      ("the two cells of a pair are POSITIONALLY ADJACENT"), which is the actual open
      content of item (f).

  H3 ⇒⇒ **THE NAMED RISK SET, and this is the point of the probe.**  The forward proof
      uses only that the right angle `R` is a REGULAR point of the developed complex
      (cone angle `2π`) whose two side-reflections compose to `−id`.  The CONVERSE needs
      the complementary statement -- that a graze at `O` or `A` moves the first-return
      wall.  But `A` is regular exactly when `r = Q − P ∈ {1, 2}`, and `O` exactly when
      `P ∈ {1, 2}`:  a vertex with interior angle `mπ/2Q` has cone angle
      `2π·m/gcd(m, 2Q)`, and `gcd(P, 2Q) = 1` for odd `P` while `gcd(r, 2Q) = gcd(r, 2)`.
      At a REGULAR `A` the development does not branch, so an `A`-graze could keep the
      return wall and give a `T`-continuous NON-fold boundary -- which would produce a
      3-cell run and a NON-ADJACENT `J`-pair, refuting H2 there and only there.
      ⇒ So: does the `r ∈ {1,2}` (and `P = 1`) subpopulation behave differently?
      ⚠ I do NOT predict it fails -- s371 put the whole `q = 1` family (`r = 1`, all
      ring-certified) in the store and s373's `114/114` did not flinch.  The value here
      is that the risk set is NAMED and CHECKED rather than assumed away, and that a
      pass tells the proof which sub-case still needs an argument.

⚠ NON-CIRCULARITY.  Pairing is read off `chi` (RING-measured: 2 on a `J`-pair, 1 on the
orphan), never off width and never off the graze columns -- `cell_store.independent`
confirms `chi` shares no derivation path with `graze_*` or `width` ([OPS-071]).  The
equal-width row below is therefore a genuine test and not a restatement of the pairing.
⚠ COVERAGE IS QUOTED ON EVERY LINE (CLAUDE.md): `chi` sits on the ring-certified centres
only, so every score here is a small-`Q`-biased subsample whatever the totals say.

PRIOR ART: grepped `rulings.md` for 'adjacen', 'fold', 'R-graze', 'regular', 'cone angle',
'return copy', 'NCYL-216', 'NCYL-186', 'NCYL-133', 'NCYL-046' and read `orphan_theorem.md`
Lemma B/B'/C + §B'.1 ->
  [NCYL-186] states the target directly: the adjacent-equal-width-pairs picture IS proved
    from Thm D (+) Thm A / Lemma B', EXCEPT positional adjacency, which is measured
    `198/198` and unproved.  That is item (f)'s residual, so (f) is NOT a fresh question.
  [NCYL-133] supplies the forward proof's input (`R` is the only reversing vertex, it is
    REGULAR, `T` is continuous across an `R`-graze) -- asserted there inside a discharge,
    used here as the hypothesis it is.
  [NCYL-046] (s272) left EXACTLY this open: the regular-vertex-graze route is "POSSIBLE
    but NOT ESTABLISHED, and if it never operates the formula is exact".  [NCYL-216]
    establishes that it DOES operate, on ~51% of interior boundaries -- so that clause
    needs re-pointing, and H1 is the statement that says how it operates.
  [NCYL-026] (path length constant on a cell; a perpendicular `L1` hit is a VERTICAL
    developed copy at an `x` independent of `s`) is the input for the return-wall picture.
  No prior ruling scores H1 or H3; no prior ruling mentions the `r ∈ {1,2}` regular-`A`
  case at all.

Store reads only -- no tracing, no ring runs ([OPS-094]: read the JSON, never `cs.get()`).

USAGE
    PYTHONPATH=.:engine:archive/scripts_2026-07:archive/scripts_2026-08:probes \
      .venv/bin/python3.13 probes/s374_graze_structure.py
"""
import collections
import glob
import json
import os

import cell_store as cs

OUT = 'data/s374_graze_structure.json'


def rows():
    """Centres carrying BOTH the graze columns and `chi`, read straight off disk."""
    for p in sorted(glob.glob(os.path.join(cs.STORE_DIR, '*.json'))):
        try:
            rec = json.load(open(p))
        except Exception:                                      # noqa: BLE001
            continue
        if rec.get('guard') or not rec.get('complete') or not rec.get('cells'):
            continue
        cells = sorted(rec['cells'], key=lambda d: d['lo'])
        if any(c.get('chi') is None for c in cells):
            continue
        # every interior boundary must carry a graze verdict
        if any(c.get('graze_hi_V') is None for c in cells[:-1]):
            continue
        yield rec, cells


def main():
    tot = collections.Counter()
    per_r = collections.defaultdict(collections.Counter)
    census = collections.Counter()
    fails = {'H1': [], 'H2': [], 'equal_width': [], 'orphan_R': []}
    centres = []

    for rec, cells in rows():
        P, Q = rec['P'], rec['Q']
        r = Q - P
        n = len(cells)
        n_orph = sum(1 for c in cells if c['chi'] == 1)
        bnd = [cells[i]['graze_hi_V'] for i in range(n - 1)]
        nR = sum(1 for v in bnd if v == 'R')
        npair = (n - n_orph) // 2
        for v in bnd:
            census[v] += 1
            per_r[r][f'bnd_{v}'] += 1

        # --- H1: no two consecutive R boundaries -------------------------------
        h1 = all(not (bnd[i] == 'R' and bnd[i + 1] == 'R') for i in range(len(bnd) - 1))
        # --- H2: #R == #J-pairs  (given H1, == "the only singleton is the orphan")
        h2 = (nR == npair)
        # --- the domino tiling, stated directly: every non-orphan cell is in
        #     exactly one R-domino, and the orphan is in none.
        inR = [False] * n
        for i, v in enumerate(bnd):
            if v == 'R':
                inR[i] = inR[i + 1] = True
        tiling = all(inR[i] == (cells[i]['chi'] != 1) for i in range(n))
        # --- an R boundary's flanking cells have equal width (control: non-R) ---
        for i, v in enumerate(bnd):
            eq = abs((cells[i]['hi'] - cells[i]['lo'])
                     - (cells[i + 1]['hi'] - cells[i + 1]['lo'])) < 1e-9
            tot['eqw_R' if v == 'R' else 'eqw_nonR'] += int(eq)
            tot['n_R' if v == 'R' else 'n_nonR'] += 1
            if (v == 'R') != eq:
                fails['equal_width'].append((P, Q, i, v, eq))
        # --- no orphan-cell endpoint is an R graze -----------------------------
        for c in cells:
            if c['chi'] == 1:
                for side in ('graze_lo_V', 'graze_hi_V'):
                    if c.get(side) is not None:
                        tot['orphan_bnd'] += 1
                        if c[side] == 'R':
                            fails['orphan_R'].append((P, Q, side))

        tot['centres'] += 1
        tot['H1'] += h1
        tot['H2'] += h2
        tot['tiling'] += tiling
        per_r[r]['centres'] += 1
        per_r[r]['H1'] += h1
        per_r[r]['H2'] += h2
        per_r[r]['tiling'] += tiling
        if not h1:
            fails['H1'].append((P, Q, ''.join(v[0] for v in bnd)))
        if not h2:
            fails['H2'].append((P, Q, n, n_orph, nR, npair,
                                ''.join(v[0] for v in bnd)))
        centres.append(dict(P=P, Q=Q, r=r, n=n, n_orph=n_orph, nR=nR,
                            npair=npair, H1=h1, H2=h2, tiling=tiling,
                            bnd=''.join(v[0] for v in bnd)))

    print(f'independent(chi, graze_lo_V) = '
          f'{cs.independent("chi", "graze_lo_V")};  '
          f'independent(chi, width) = {cs.independent("chi", "width")}')
    print(f'\nCOVERAGE: {tot["centres"]} centres carry BOTH graze columns and `chi` '
          f'(of {len(cs.known_centres())} known).  ⚠ `chi` is ring-certified-only, '
          f'hence small-`Q` biased.')
    print(f'  interior boundaries scored: {sum(census.values())}   census '
          f'{dict(census)}')
    print(f'\nH1  no two consecutive R-boundaries : '
          f'{tot["H1"]}/{tot["centres"]} centres')
    print(f'H2  #R == #J-pairs                  : '
          f'{tot["H2"]}/{tot["centres"]} centres')
    print(f'    domino tiling (non-orphan <=> in an R-domino): '
          f'{tot["tiling"]}/{tot["centres"]} centres')
    print(f'    R-boundary flanking widths EQUAL   : '
          f'{tot["eqw_R"]}/{tot["n_R"]}')
    print(f'    CONTROL non-R flanking widths EQUAL: '
          f'{tot["eqw_nonR"]}/{tot["n_nonR"]}  (want 0)')
    print(f'    orphan-cell endpoints that are R   : '
          f'{len(fails["orphan_R"])}/{tot["orphan_bnd"]}  (want 0)')

    print('\nH3 -- THE NAMED RISK SET.  `A` is a REGULAR point iff r = Q-P in {1,2}; '
          '`O` iff P in {1,2}.')
    print(f'{"r":>4} {"centres":>8} {"H1":>6} {"H2":>6} {"tiling":>7}   '
          f'{"bnd R":>6} {"bnd A":>6} {"bnd O":>6}')
    for r in sorted(per_r):
        c = per_r[r]
        flag = '  <-- A REGULAR' if r in (1, 2) else ''
        print(f'{r:>4} {c["centres"]:>8} {c["H1"]:>6} {c["H2"]:>6} '
              f'{c["tiling"]:>7}   {c["bnd_R"]:>6} {c["bnd_A"]:>6} '
              f'{c["bnd_O"]:>6}{flag}')

    for k, v in fails.items():
        if v:
            print(f'\n⚠ {k} FAILURES ({len(v)}): {v[:12]}')

    json.dump(dict(totals=dict(tot), census=dict(census),
                   per_r={str(k): dict(v) for k, v in per_r.items()},
                   fails={k: v[:200] for k, v in fails.items()},
                   centres=centres), open(OUT, 'w'), indent=1)
    print(f'\nwrote {OUT}')


if __name__ == '__main__':
    main()
