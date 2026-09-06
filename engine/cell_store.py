#!/usr/bin/env python3
r"""cell_store.py -- THE per-branch measurement store (`future_directions.md` E1).

⇒⇒ **CONSULT THIS INSTEAD OF RE-PARTITIONING.**  One record per centre `(P,Q)`,
one row per BRANCH (cell), holding everything a session has ever wanted to know
about a cell in ONE place, so nobody assembles it from scratch again.

    from cell_store import get, cells, export_csv, recap, peel
    rec = get(5, 17)            # from disk, or measured and added
    for c in rec['cells']: ...  # width, N, T/chi, chi, area, m, gb, ...
    export_csv('data/branch_widths.csv')
    recap(caps={'budget': 300, 'max_nodes': 8_000_000})   # push the unresolved on

WHY THIS EXISTS (and why it is not another one-off).  `future_directions.md` E1
is the user's s252 proposal for exactly this store.  It sat unbuilt for ~90
sessions and its own post-mortem names the two reasons -- (i) it was filed off
the roster behind a compliance-policed trigger, and (ii) **the design crowded out
the build**.  s344 shipped the first useful subset (`ring_cache`, memoising the
partition).  This is the second: the CELL-LEVEL store, at the home E1 names.
⚠ The fan-chain layer (E1's "store one level below the hypothesis") is the part
that did the crowding out; `fan_chain` is reserved in the schema and left null in
v1 so it is unblocked without repeating the history.

** THE TWO ENGINES ARE COMPLEMENTARY, NOT REDUNDANT -- that is why rows merge. **
    glen  (`s340_glen_enum.enum_glen`)  w, N (= closure), T/chi (= fwd), nfan,
          jstar.  Float, NO depth cap, 22-38x the ring.  Never a certificate.
    ring  (`ring_cache.cells_cached`)   w (exact), T, chi, depth.  Bit-exact,
          but depth-capped, and it carries NO word length.
So `N` can only come from glen and `chi`/`T` only from the ring; a row that has
both is a merge, matched on position, and is flagged `certified`.

FIELDS, per cell
    cell_index lo hi            position on the transversal
    width                       cell width           (both engines; cross-checked)
    N                           word length, in reflections          (glen only)
    T_over_chi                  first-return length                  (glen; ring)
    chi                         multiplicity: 2 on a J-pair, 1 on the orphan (ring)
    T                           = T_over_chi * chi                   (ring only)
    area_wN = width * N         combinatorial; the quantity bounded by `8*a1`
    area_wT = width * T_over_chi  geometric; ** sum over a centre = S **
    frac_wN frac_wT             each area / the centre's TOTAL of that area --
                                so each sums to 1 over a centre, which is what
                                makes cells COMPARABLE ACROSS CENTRES (S and the
                                Kac total both grow with Q, so the raw areas are
                                not).  `frac_wT` = share of the swept sum;
                                `frac_wN` = share of the Kac mass.
    m                           rung index, or None if not a rung cell
    m_source                    how `m` was assigned
    gb hmax e_end               winding ceiling `G_B` / per-fan `H`-cap / closing
                                winding (0 or +-Q)       (glen; s357, re-sourced s383)
    nfan jstar                  O-fan count / exit vertex index      (glen only)
    is_orphan                   the cell with no equal-width twin -- so it needs
                                NO ring ([OPS-191], s431); `chi == 1` only as a
                                fallback on PARTIAL records.  `is_orphan_source`
    in_width_class              width >= t ~ r2/P (where `w*N <= 8*a1` is proved)
    fan_chain                   RESERVED, null in v1 (E1's layer)

    nopp                        Lemma-B other-branch count; 0 on the leaf class
    vertical                    the ring's vertical flag                (ring)

FIELDS, per centre -- including E1's NEGATIVE RECORD, the field it says no probe
writes ("partitioned to depth D; rungs present = {...}", not merely what was found)
    P Q r a1 a2 cf              the centre
    n_cells sumw S g_len        totals
    rungs_present               ** the negative record **
    complete                    ** did the partition CLOSE ** -- the only completion flag
    pending pending_noexit      ** resumable node state for what did NOT close **
    wpending wnoexit wsliver    the width each unresolved class accounts for
    caps attempts               what was asked of the engine, and what it bought
    kac                         the peel ledger (see THE KAC PEEL below)
    depth_read capped noexit sliver graze          what was NOT resolved
    engine certified guard store_tag

** A CAPPED CENTRE IS PARTIAL, NOT LOST (s354).  ** Its resolved cells are kept and
its unresolved subintervals are stored as glen node states, so `recap()` RESUMES
them -- cost proportional to what is left, not to the whole transversal, and the
continuation is bit-identical to never having stopped (checked against full runs at
(5,17), (3,11), (9,31), (7,31), (4,15), one of them across 87 interruptions).
⚠ Raising `budget` alone does NOT push a centre further: glen stops on four
independent caps and s353 exposed one.  See `CAPS`.

THE CERTIFICATE LIVES IN THE READER (E1's strongest stated reason to build it:
s249, s251 and s252 each re-implemented it and s252 got it wrong).  `get()` runs
`verify()` on every record it returns, cached or fresh:
    sum of widths == 1            (a partition, not a stack -- see the DO-NOT)
    sum of area_wT == S           (the swept-sum identity)
    r == Q mod P, a1 == Q // P
    width in (0, 1] for every cell

⚠⚠ **DO-NOT -- THE MISTAKE THAT PRODUCED THIS FILE.** Do NOT assemble a branch
table out of the RUNG stores (`s252_fan_chain`, `s253_wide_scan`,
`s312_offladder_props`).  Those are stacks and filtered subsets, not partitions:
the leftover band `|Y_{m*}|` and every unselected cell are absent, so the widths
sum to 0.55-0.9998 and never to 1.  s353 shipped exactly that table and the user
caught it on sight.  And do NOT merge cells ACROSS stores to repair it: different
stores hold partitions of the SAME centre at different depths, so a naive merge
double-counts (sum 4.087 at `(7,31)`).  ONE partition per centre, and `verify()`
now makes both errors impossible to ship.

⚠ `m` is the RUNG INDEX (`notation.md`), not the ladder step.  `enum_glen`'s third
argument is `a1`; on an alpha->0 ladder `Q = P*a1 + r` those coincide numerically
and the collision has bitten before.
"""
import collections
import glob
import json
import math
import os
import time

STORE_DIR = os.environ.get('CELL_STORE_DIR', 'data/cell_store')
STORE_TAG = 's354-v2'        # bump when an engine's semantics change (ring_cache (i))
_RUNG_SOURCES = ['s253_wide_scan.json', 's252_fan_chain.json',
                 's312_offladder_props.json']
_DATA = 'data'

# The four INDEPENDENT caps glen stops on.  s353 exposed only `budget`, which is
# why raising it stopped helping: at (10,23) `budget=20` stops on the clock at
# 221193 nodes / 20.3 s and `budget=120` stops on `max_nodes` at 400010 / 35.6 s,
# with `sumw` identical to the last digit.  A recap must raise all four.
CAPS = dict(budget=20.0, max_nodes=400000, max_fans=2000000, max_letters=40000000)


# --------------------------------------------------------------------------
# THE FIELD SCHEMA -- every cell column's PROVENANCE, declared in one place.
#
# ⚠⚠ WHY IT EXISTS ([OPS-071], s355, and it cost most of a session).  Scoring
# A1's rung law `#fans = min(4m+2, 4(a1-m))` against the `m` column read
# `618/620` -- which would have said the law holds everywhere, including the
# `2r > P` cone where an independent predictor shows it FAILS `80/84`.  The
# column is not a measurement: `_m_from_shape` computes `m` by INVERTING that
# law from `nfan`, so 444 of those rows COULD NOT FAIL and read `444/444`.  The
# only thing that caught it was `m_source`.  This generalises that field to the
# whole schema, so the question "is this column independent of the thing I am
# scoring?" is a function call instead of a reading of the source.
#
#     kind 'measured' -- an engine emitted it; detail names the engine.
#     kind 'derived'  -- a pure function of other STORED fields; `rederive()`
#                        owns it and detail lists its INPUTS.  Derived columns
#                        are fine BECAUSE `rederive` owns them: they cannot
#                        drift.  What is never fine is a column that DUPLICATES
#                        a measured one (`nb1 = chi*nfan` -- do not add it).
#     kind 'joined'   -- copied in from another stored table by key.
#     kind 'mixed'    -- provenance varies PER CELL; detail lists the inputs of
#                        the derived path, and a `*_source` field says which
#                        path spoke on each row.  ⚠ `derived_from` reports the
#                        derived path's inputs for these -- the CONSERVATIVE
#                        answer, which is the one a circularity check wants.
#     kind 'meta'     -- bookkeeping about the row, not data about the object.
#
# ⇒⇒ ADDING A FIELD?  DECLARE IT HERE FIRST.  `verify()` rejects a record
# carrying an undeclared cell field, so the schema cannot silently drift behind
# the data -- that is the accretion guard, and it is why this is structural
# rather than a rule someone has to remember.
CELL_FIELDS = {
    'lo':             ('measured', 'glen|ring: left endpoint'),
    'hi':             ('measured', 'glen|ring: right endpoint'),
    'width':          ('measured', 'glen|ring: hi - lo'),
    'N':              ('measured', 'glen: closure -- word length in reflections'),
    'T_over_chi':     ('measured', 'glen: fwd -- first-return length'),
    'nfan':           ('measured', 'glen: O-fan count == nb1, the FIRST-RETURN '
                                   'block count ([NCYL-182], 132/132)'),
    'jstar':          ('measured', 'glen: exit vertex index'),
    'nopp':           ('measured', 'glen: Lemma-B other-branch count'),
    'chi':            ('measured', 'ring: 2 on a J-pair, 1 on the orphan'),
    'T':              ('measured', 'ring: full-period length'),
    'depth':          ('measured', 'ring: tree depth at which the cell closed'),
    'vertical':       ('measured', 'ring: vertical-exit flag'),
    # ⚠⚠ RE-DECLARED s431 ([OPS-191], recorded s429).  It was ('derived', ('chi',))
    # and therefore RING-ONLY at 40.4% of cells -- and the ring-certified centres
    # are the small-`Q` ones, so the column was BIASED as well as partial
    # ([OPS-179]'s family).  But the orphan needs no ring: a J-pair has EQUAL
    # widths, so on a COMPLETE partition it is the unique cell whose width has no
    # twin -- a function of `width`, which is 100% covered because `verify()`
    # asserts `Σw = 1`.  'mixed' rather than plain 'derived' only because the
    # width rule needs the partition CLOSED: on the 104 partial records `chi` is
    # still the honest answer and there is no reason to throw a measurement away.
    # `is_orphan_source` says which path spoke.  Coverage 4512 -> 9109 cells.
    # ⚠ The detail tuple lists BOTH paths' inputs, not just the derived path's as
    # `m` does, and the asymmetry is deliberate: for `m` the DERIVED path is the
    # circular one, so reporting it is conservative; here it is the FALLBACK that
    # can be circular (against `chi`), so reporting only `width` would make
    # `independent('is_orphan', 'chi')` read True while 34 cells are chi-set.
    # Union = the conservative answer; recover the precision by splitting on
    # `is_orphan_source` and scoring the `width:unpaired` rows only.
    'is_orphan':      ('mixed',    ('width', 'chi')),
    'is_orphan_source': ('meta',   'which path set `is_orphan` on THIS cell: '
                                   '"width:unpaired" (the rule, COMPLETE records) '
                                   'or "ring:chi" (partial records only) -- s431'),
    'area_wN':        ('derived',  ('width', 'N')),
    'area_wT':        ('derived',  ('width', 'T_over_chi')),
    'frac_wN':        ('derived',  ('area_wN',)),
    'frac_wT':        ('derived',  ('area_wT',)),
    'cell_index':     ('derived',  ('lo',)),
    'in_width_class': ('derived',  ('width',)),
    'm':              ('mixed',    ('nfan', 'N')),
    'm_source':       ('meta',     'which path set `m` on THIS cell: '
                                   '"join:rung-stores" (independent) or '
                                   '"shape:nfan+N" (A1 INVERTED -- see [OPS-071])'),
    'm_geom':         ('derived',  ('width',)),
    'm_geom_M':       ('meta',     'the centre\'s stack top M on THIS cell\'s row '
                                   '(None outside the validated 2r<P cone)'),
    'gb':             ('measured', 'glen: G_B = max_f |E_f|, the winding '
                                   'ceiling ([WFLOOR-112]) -- s357.  ⚠ DECLARED '
                                   "'joined' from s353 to s357 and populated on "
                                   '0/5009 cells the whole time ([OPS-073]); the '
                                   'rung-table join never existed.  ⚠ SOURCE CHANGED '
                                   's383: pointwise `s357_fanword.trace_fanword` -> '
                                   "glen's own walk, same values (4992/4992 all three "
                                   'columns) -- [OPS-114] applied to `fill_fanword`'),
    'hmax':           ('measured', 'glen: max per-fan H-letter count '
                                   '= max_f |E_f - E_{f-1}| ([CHARGE-029]/[NCYL-169]).'
                                   '  ⚠ SOURCE CHANGED s383, as `gb`'),
    'e_end':          ('measured', 'glen: the winding at closure -- 0 (exact '
                                   'cancellation) or +-Q (a wrap).  ⚠ SOURCE CHANGED '
                                   's383, as `gb`'),
    'graze_lo_V':     ('measured', "glen: which vertex ('O'/'R'/'A') the graze "
                                   "that DEFINES this cell's LOWER boundary is at; "
                                   'null at the transversal end s=0.  ⚠ SOURCE CHANGED s380: '
                                   'ringwords -> glen, same values '
                                   '(3555/3555 both columns) -- [OPS-114]'),
    'graze_lo_k':     ('measured', 'glen: BOUNCES BEFORE that graze (the ray '
                                   'launched at `lo` grazes `graze_lo_V` on the step '
                                   'after its k-th reflection); null at s=0 -- s373'),
    'graze_hi_V':     ('measured', "glen: as `graze_lo_V` for the UPPER boundary; "
                                   'null at the transversal end s=1'),
    'graze_hi_k':     ('measured', 'glen: as `graze_lo_k` for the UPPER boundary; '
                                   'null at s=1'),
    'fan_chain':      ('meta',     'RESERVED (E1); null in v1'),
}
# ⚠⚠ THE GRAZE COLUMNS ARE SHARED BETWEEN ADJACENT ROWS BY CONSTRUCTION ([OPS-097]).
# A boundary belongs to TWO cells, so `graze_hi_*` of cell `i` IS `graze_lo_*` of cell
# `i+1`, always, on every centre.  That agreement is not a measurement and cannot fail
# -- do NOT score it ([OPS-041]).  They are stored per cell anyway because the cell is
# the store's key and "this cell's boundaries" is what a reader wants; the redundancy is
# ACROSS rows, not a duplicate of another column (which is the thing the measurement
# discipline actually bars).  ⚠ `graze_*_k` is NOT a proxy for `N` or `depth`: it is the
# FIRST graze along the orbit and runs far shallower (8/15 cell 0: N=53, hi-graze k=12).


def provenance(field):
    """(kind, detail) for a cell field. Raises KeyError on an undeclared one."""
    return CELL_FIELDS[field]


def derived_from(field, _seen=None):
    """The TRANSITIVE set of stored fields `field` is computed from.

    Empty for a measurement.  For `kind = 'mixed'` this reports the DERIVED
    path's inputs even though some rows come from the independent path -- the
    conservative answer, which is what a circularity check should get.
    """
    seen = set() if _seen is None else _seen
    kind, detail = CELL_FIELDS[field]
    if kind not in ('derived', 'mixed'):
        return set()
    for src in detail:
        if src in seen:
            continue
        seen.add(src)
        seen |= derived_from(src, seen)
    return seen


def independent(field, of):
    """Can `field` be scored against `of` without circularity?

    ⇒ CALL THIS BEFORE SCORING A LAW AGAINST A STORED COLUMN.  It is the
    one-line form of the check that [OPS-071] cost a session for:
        independent('m', 'nfan')  -> False   (m's shape path inverts A1 from nfan)
        independent('nfan', 'chi') -> True   (two engines, neither sees the other)
    ⚠ False does NOT mean the column is useless -- it means split on its
    `*_source` field and score only the independent rows.
    """
    return of not in derived_from(field) and field not in derived_from(of)


# --------------------------------------------------------------------------
def _path(P, Q):
    return os.path.join(STORE_DIR, f'{P}_{Q}.json')


def _cf(P, Q):
    """Continued fraction of Q/P, as the repo writes digits."""
    d, a, b = [], Q, P
    while b:
        d.append(a // b)
        a, b = b, a % b
    return d


_RUNG_LUT = None


def _rung_lut():
    """(P,Q) -> {rounded width: m}, joined from STORED rung records.

    `m` is a property of the J-pair stack; the reliable stored sources are the
    rung tables.  A width claimed by two different `m` at one centre is dropped.
    Not derived from `nfan`: (nfan-2)/4 is the IN-SCOPE branch only and warm
    rungs alias ([WFLOOR-042]).
    """
    global _RUNG_LUT
    if _RUNG_LUT is not None:
        return _RUNG_LUT
    lut, amb = {}, {}
    def add(P, Q, m, w):
        if w is None or m is None or m == '':
            return
        k, wr = (P, Q), round(float(w), 9)
        d = lut.setdefault(k, {})
        if wr in d and d[wr] != m:
            amb.setdefault(k, set()).add(wr)
        d[wr] = m
    try:
        for d_ in json.load(open(os.path.join(_DATA, 's253_wide_scan.json'))):
            add(d_['P'], d_['Q'], d_['m'], d_['U'] - d_['X'])
        for d_ in json.load(open(os.path.join(_DATA, 's252_fan_chain.json'))):
            add(d_['P'], d_['Q'], d_['m'], d_['w'])
        for d_ in json.load(open(os.path.join(_DATA, 's312_offladder_props.json'))):
            add(d_['P'], d_['a1'] * d_['P'] + d_['r'], d_.get('m'), d_['w'])
    except FileNotFoundError:
        pass
    for k, ws in amb.items():
        for w in ws:
            lut[k].pop(w, None)
    _RUNG_LUT = lut
    return lut


def _geom_rung_lut(P, Q, r, a1):
    """{rounded predicted cell width: m} for the GEOMETRIC rung stack, or None.

    ⇒⇒ THE NON-CIRCULAR RUNG INDEX, and the reason it exists (s356).  `m` above
    is `('mixed', ('nfan','N'))`: where the rung-table join has no row it falls
    back to `_m_from_shape`, which INVERTS A1 from `nfan`, so `independent('m',
    'nfan')` is False and [OPS-071] bars scoring any fan law against it.  This
    path never sees a fan count.  It reads the stack top `M` from `stack_top`
    (s257 (3), a sign test on `X_m`; `M <= a1-1` PROVED [WFLOOR-045]) and each
    rung's PREDICTED CELL WIDTH `min(V_m, X_m)` (the `want` list of s257's T5),
    and matches it against the MEASURED width -- so `independent('m_geom',
    'nfan')` is True.

    Rung `m` is a J-PAIR: the LUT maps one width to `m` and TWO cells take it.

    ⚠ SCOPE: the `2r < P` cone with `a1 >= 2`, which is where `stack_top`/`Cell`
    are validated (6916/6916, s257 (3)).  None outside it -- deliberately, rather
    than a number nobody has scored.
    ⚠ It is NOT a replacement for `m`: on the `2r < P` cone `m_geom` covers the
    RUNG cells only and is null on the residual tail, which is the honest answer
    (s356 measured the tail to hold no rung: 0/293 above the stack top).
    """
    if r == 0 or a1 < 2 or 2 * r >= P:
        return None, None
    try:
        from s254_avoidance_proof import Cell
        from s257_warmup_tiling import stack_top
    except Exception:                       # probes/ not on the path
        return None, None
    M, _ = stack_top(P, r, a1)
    lut = {}
    for m in range(M + 1):
        z = Cell(P, Q, r, a1, m)
        w = min(z.V, z.X)
        if w > 0:
            lut.setdefault(round(float(w), 9), m)
    return lut, M


def _closure_pred(a1, m):
    """The §s241/§s257 rung closure law -- the LAW, not a measurement."""
    return min(8 * (2 * m + 1) * a1 + 1 - 16 * m * m,
               16 * (a1 - m) * (m + 1) - 3)


def _m_from_shape(a1, nfan, N):
    """Rung index from the cell's own SHAPE: `nfan` gives two candidates --
    in-scope `(nfan-2)/4` and warm `a1 - nfan/4` ([WFLOOR-042]: they alias, which
    is why `nfan` alone must never be used) -- and the MEASURED word length `N`
    picks between them against `_closure_pred`.  None if 0 or 2 survive.

    ⚠ Validated against the independent rung-table join on its 168 overlapping
    rows: 166 agree.  The 2 that do not are `(11,71)`, a DEGENERATE rung where
    the measured `N = 189` departs from the law's `177` and collides with
    `closure_pred(6,3) = 189` -- so the join wins wherever it exists, and this is
    the fallback.  A degenerate rung can therefore carry a wrong shape-derived
    `m`; `m_source` records which method spoke.
    """
    if nfan is None or N is None:
        return None
    cands = set()
    if (nfan - 2) % 4 == 0:
        cands.add(int((nfan - 2) // 4))
    if nfan % 4 == 0:
        cands.add(int(a1 - nfan // 4))
    good = [m for m in cands if 0 <= m <= a1 and _closure_pred(a1, m) == int(N)]
    return good[0] if len(good) == 1 else None


# --------------------------------------------------------------------------
def _fnum(x):
    try:
        return float(x)
    except Exception:                                          # noqa: BLE001
        return None


def _inum(x):
    """As `_fnum` for a column that is an INTEGER (the winding trio, s383).

    ⚠ Kept separate on purpose: `gb`/`hmax`/`e_end` are exact small integers, they
    are compared with `==` by every consumer and by the source gate, and `e_end`'s
    whole content is that it is `0` or `+-Q` ([NCYL-188]).  Floating them would make
    the pre-s383 rows (written as ints by `fill_fanword`) and the new ones different
    JSON types for the same measurement.
    """
    try:
        return int(x)
    except Exception:                                          # noqa: BLE001
        return None


def _orphan_by_width(P, cells, rtol=1e-9):
    """(index of the unpaired-width cell | None, resolved?) -- the [OPS-191] rule.

    A J-pair has EQUAL widths (Thm D + `orphan_theorem.md` Lemma B'), so on a
    COMPLETE partition the orphan is the unique cell whose width has no twin;
    even `P` has no orphan at all ([NCYL-113]) and must show ZERO unpaired.
    ⚠ `resolved=False` is a REFUSAL, not an answer: a tolerance failure has to
    surface as "the rule declines" rather than as a wrong orphan.  Caller must
    pass `cells` already sorted, since the index is positional.
    """
    ws = sorted((float(c['width']), i) for i, c in enumerate(cells))
    paired = [False] * len(ws)
    i = 0
    while i < len(ws) - 1:
        if math.isclose(ws[i][0], ws[i + 1][0], rel_tol=rtol):
            paired[i] = paired[i + 1] = True
            i += 2
        else:
            i += 1
    lone = [k for k, p in enumerate(paired) if not p]
    if P % 2 == 0:
        return None, not lone
    if len(lone) != 1:
        return None, False
    return ws[lone[0]][1], True


def _derive(rec, base):
    """Rung index, width class, the orphan flag and the derived areas, in place."""
    P, Q, r, a1 = rec['P'], rec['Q'], rec['r'], rec['a1']
    lut = _rung_lut().get((P, Q), {})
    glut, gM = _geom_rung_lut(P, Q, r, a1)
    r2 = P - 2 * r if r else None
    t = (r2 / P) if (r2 and r2 > 0) else None
    base.sort(key=lambda c: c['lo'])
    # THE ORPHAN, FROM THE WIDTHS ALONE ([OPS-191]) -- see the `is_orphan` note in
    # CELL_FIELDS.  Gated on `complete` because a partial record's missing cells
    # can be the orphan's twin, which would make an unpaired width meaningless.
    o_idx, o_ok = ((None, False)
                   if not (rec.get('complete')
                           and all(c.get('width') is not None for c in base))
                   else _orphan_by_width(P, base))
    for i, c in enumerate(base):
        c.setdefault('chi', None)
        c.setdefault('T', None)
        c.setdefault('depth', None)
        if o_ok:
            c['is_orphan'] = (i == o_idx)
            c['is_orphan_source'] = 'width:unpaired'
        elif c.get('chi') is not None:
            c['is_orphan'] = (int(c['chi']) == 1)
            c['is_orphan_source'] = 'ring:chi'
        else:
            c['is_orphan'] = None
            c['is_orphan_source'] = None
        c['cell_index'] = i
        w, N, tc = c['width'], c['N'], c['T_over_chi']
        c['area_wN'] = (w * N) if (w is not None and N is not None) else None
        c['area_wT'] = (w * tc) if (w is not None and tc is not None) else None
        m = lut.get(round(w, 9)) if w is not None else None
        src = 'join:rung-stores' if m is not None else None
        if m is None:                       # fallback: the cell's own shape
            m = _m_from_shape(a1, c.get('nfan'), c.get('N'))
            src = 'shape:nfan+N' if m is not None else None
        c['m'] = m
        c['m_source'] = src
        c['m_geom'] = (glut.get(round(w, 9))
                       if (glut is not None and w is not None) else None)
        c['m_geom_M'] = gM
        c.setdefault('gb', None)
        c['in_width_class'] = (None if t is None else bool(w >= t))
        c.setdefault('fan_chain', None)             # RESERVED (E1); null in v1
    # ---- SHARE OF THE CENTRE'S TOTAL AREA -----------------------------------
    # `frac_wT` is the cell's share of the SWEPT SUM (denominator = S, the
    # verified `Σ w·T/χ`); `frac_wN` is its share of the KAC MASS (denominator =
    # `Σ w·N = 2μ_∂ − 1` where that is exact).  Both sum to 1 over a centre by
    # construction, which is what makes them comparable ACROSS centres -- the raw
    # areas are not, since S and the Kac total both grow with Q.
    # ⚠ On a PARTIAL record the denominator is the RESOLVED total, so the
    # fractions are shares of what has been read, not of the centre.  `complete`
    # says which you are looking at; `export_csv` drops partials by default.
    tot_wN = sum(c['area_wN'] for c in base if c['area_wN'] is not None)
    tot_wT = sum(c['area_wT'] for c in base if c['area_wT'] is not None)
    for c in base:
        c['frac_wN'] = (c['area_wN'] / tot_wN
                        if (c['area_wN'] is not None and tot_wN) else None)
        c['frac_wT'] = (c['area_wT'] / tot_wT
                        if (c['area_wT'] is not None and tot_wT) else None)
    rec['tot_wN'], rec['tot_wT'] = tot_wN, tot_wT
    rec['cells'] = base
    rec['n_cells'] = len(base)
    rec['rungs_present'] = sorted({c['m'] for c in base if c['m'] is not None})
    return rec


def measure(P, Q, max_depth=200000, ring='cached', caps=None, budget=None):
    """Measure a centre with BOTH engines and merge.  Never consults the store.

    `ring`: 'cached' reads the ring ONLY if `ring_cache` already has the row
    (disk read, free); 'compute' lets it run the enumerator, which can take
    minutes per centre; 'off' skips the ring entirely.  Default 'cached' --
    backfill must not trigger fresh ring runs, which is what made the first
    backfill attempt crawl (28 centres in 5 minutes).

    `caps`: overrides for `CAPS` (budget / max_nodes / max_fans / max_letters).
    glen's cost is linear in TOTAL LETTERS and letters are geometric in `a1`
    ([NCYL-153]), so a deep centre can run for many minutes.  ⚠ `budget=` is
    kept as a legacy alias for `caps={'budget': ...}` and is NOT enough on its
    own -- see the `CAPS` comment.

    ** A CAPPED CENTRE IS A LEGITIMATE, PARTIAL, RESUMABLE RECORD ** -- its
    resolved cells are KEPT, its unresolved subintervals are stored as `pending`
    node states, and `complete` is False.  s353 threw the partial cells away and
    stored a counter, so the only way forward was to redo the whole transversal.
    """
    from s340_glen_enum import enum_glen                       # probes/, lazy
    caps = dict(CAPS, **(caps or {}))
    if budget is not None:
        caps['budget'] = budget
    r, a1 = Q % P, Q // P
    rec = dict(P=P, Q=Q, r=r, a1=a1, a2=(P // r if r else None), cf=_cf(P, Q),
               store_tag=STORE_TAG, engine=None, certified=False, guard=None,
               depth_read=max_depth, caps=caps, cells=[], attempts=[])

    # ---- glen: the only source of the word length N -------------------------
    gcells, gst, secs = None, {}, 0.0
    gmap = {}
    try:
        t0 = time.time()
        gcells, gst = enum_glen(P, r, a1, collect_pending=True, **caps)
        # ⚠ CAPTURE IMMEDIATELY.  `last_grazes` is a module attribute; on a raise it
        # would otherwise still hold the PREVIOUS centre's map, which would be merged
        # onto this row's boundaries by float proximity -- silently, and wrongly.
        gmap = dict(enum_glen.last_grazes)
        secs = time.time() - t0
    except Exception as e:                                     # noqa: BLE001
        rec['guard'] = f'glen raised: {type(e).__name__}: {e}'
    glen_ok = bool(gst.get('ok'))

    # ---- ring: the only source of chi / T, and the certificate --------------
    rrows = None
    if ring != 'off':
        try:
            import ring_cache
            got = (ring_cache.load(P, Q, max_depth) if ring == 'cached'
                   else ring_cache.cells_cached(P, Q, max_depth=max_depth))
            if got is not None:
                _ctx, rrows, rguard = got
                if rguard:
                    rrows = None
        except Exception:                                      # noqa: BLE001
            rrows = None

    if gcells is None and rrows is None:
        # E1's NEGATIVE RECORD: a well-formed row saying what was attempted and
        # how far it got -- never a truncated dict.
        rec['guard'] = rec['guard'] or 'both engines failed'
        rec.update(n_cells=0, cells=[], sumw=None, S=None, g_len=None,
                   rungs_present=[], complete=False, pending=[],
                   pending_noexit=[], wpending=None, wnoexit=None,
                   capped=gst.get('capped'), noexit=gst.get('noexit'),
                   sliver=gst.get('sliver'), wsliver=None,
                   graze=gst.get('graze'), engine='none')
        return rec

    # The RING wins the spine when glen came back partial: a ring row is a
    # COMPLETE partition, so it turns an unresolvable centre into a resolved one
    # (glen's partial N values are then merged onto it by position).
    use_ring_spine = (gcells is None) or (not glen_ok and rrows is not None)

    if not use_ring_spine:
        base = [dict(lo=_fnum(c['lo']), hi=_fnum(c['hi']), width=_fnum(c['w']),
                     N=_fnum(c.get('closure')), T_over_chi=_fnum(c.get('fwd')),
                     nfan=_fnum(c.get('nfan')), jstar=_fnum(c.get('jstar')),
                     nopp=_fnum(c.get('nopp')),
                     # the winding trio rides along on the same glen pass (s383)
                     gb=_inum(c.get('gb')), hmax=_inum(c.get('hmax')),
                     e_end=_inum(c.get('e_end')))
                for c in gcells]
        rec['engine'] = 'glen'
        rec.update(S=_fnum(gst.get('S')), g_len=_fnum(gst.get('g_len')),
                   sumw=_fnum(gst.get('sumw')), complete=glen_ok,
                   pending=gst.get('pending') or [],
                   pending_noexit=gst.get('noexit_nodes') or [],
                   wpending=_fnum(gst.get('wpending')),
                   wnoexit=_fnum(gst.get('wnoexit')),
                   capped=gst.get('capped'), noexit=gst.get('noexit'),
                   sliver=gst.get('sliver'), wsliver=_fnum(gst.get('wsliver')),
                   graze=gst.get('graze'), maxclos=gst.get('maxclos'))
    else:
        base = [dict(lo=_fnum(d['lo']), hi=_fnum(d['hi']),
                     width=float(d['w'].to_mpf(40)), N=None,
                     T_over_chi=float(d['T'].to_mpf(40) / d['chi']),
                     nfan=None, jstar=None, nopp=None) for d in rrows]
        rec['engine'] = 'ring'
        rec.update(S=None, g_len=None, sumw=sum(c['width'] for c in base),
                   complete=True, pending=[], pending_noexit=[],
                   wpending=0.0, wnoexit=0.0, capped=None, noexit=None,
                   sliver=None, wsliver=None, graze=None, maxclos=None)
        if gcells:                    # glen's partial word lengths, by position
            for c in base:
                hit = next((g for g in gcells
                            if abs(_fnum(g['lo']) - c['lo']) < 1e-9
                            and abs(_fnum(g['hi']) - c['hi']) < 1e-9), None)
                if hit is not None:
                    c['N'] = _fnum(hit.get('closure'))
                    c['nfan'] = _fnum(hit.get('nfan'))
                    c['jstar'] = _fnum(hit.get('jstar'))
                    c['nopp'] = _fnum(hit.get('nopp'))
                    c['gb'] = _inum(hit.get('gb'))
                    c['hmax'] = _inum(hit.get('hmax'))
                    c['e_end'] = _inum(hit.get('e_end'))

    # the graze pair rides along on glen's own pass ([OPS-114]) -- no second engine,
    # no separate backfill.  A ring-spine row gets it too whenever glen also ran.
    _merge_graze(base, gmap)
    _merge_ring(rec, base, rrows)
    _derive(rec, base)
    rec['attempts'] = [_attempt(caps, gst, secs, rec, resumed=False)]
    rec['kac'] = peel(rec)
    return rec


def _merge_ring(rec, base, rrows, ring_complete=True):
    """Position-matched merge of the ring's chi / T / depth / vertical.

    ⚠⚠ `ring_complete=False` (s378) merges a PARTIAL ring row -- the cells a
    guarded run did certify.  It fills only the cells it covers (the position
    match below already skips the rest) and must NOT set `certified`.  Before
    s378 this flag was set unconditionally on any non-None `rrows`, which was
    safe only because a guarded `cells_cached` returned `rows = None`; feeding
    partial rows through the old path would have marked uncertified centres
    certified, which is the one way this change could corrupt the store.
    """
    if rrows is None:
        return
    ring = sorted(rrows, key=lambda d: d['lo'])
    rec['engine'] = 'ring' if rec['engine'] == 'ring' else 'glen+ring'
    if ring_complete:
        rec['certified'] = True
        rec.pop('ring_partial', None)
    else:
        rec['ring_partial'] = True
    for c in base:
        hit = next((d for d in ring if abs(d['lo'] - c['lo']) < 1e-9
                    and abs(d['hi'] - c['hi']) < 1e-9), None)
        if hit is None:
            continue
        c['chi'] = int(hit['chi'])
        c['T'] = float(hit['T'].to_mpf(40))
        c['depth'] = hit.get('depth')
        c['vertical'] = hit.get('vertical')
        # ⚠ `is_orphan` is NOT set here (s431).  It is owned by `_derive`, which
        # every caller of this function runs immediately afterwards, and which
        # prefers the width rule to `chi` ([OPS-191]).  Two writers would be two
        # provenances for one column -- the thing `is_orphan_source` exists to
        # make impossible.
        if c['T_over_chi'] is None:
            c['T_over_chi'] = c['T'] / c['chi']


def _attempt(caps, gst, secs, rec, resumed):
    """One row of the ATTEMPT LOG -- what was asked for, what it bought.

    The log is what makes the Kac peel predictive: two rows give the exponent in
    `next_bound`, and without it every recap re-guesses from nothing.
    """
    prev = rec.get('attempts') or []
    nodes = int(gst.get('nodes') or 0)
    return dict(caps=dict(caps), resumed=bool(resumed), secs=round(secs, 2),
                nodes=nodes,
                nodes_cum=nodes + sum(a.get('nodes', 0) for a in prev),
                capped=gst.get('capped'), noexit=gst.get('noexit'),
                n_cells_cum=rec.get('n_cells'), sumw_cum=rec.get('sumw'),
                kac_peeled_cum=_kac_peeled(rec), complete=rec.get('complete'))


# --------------------------------------------------------------------------
# THE KAC PEEL -- a progress meter with a denominator known BEFORE you enumerate.
#
# `sumw` is a bad progress meter: it is dominated by the WIDE cells, which are the
# cheap ones, so a run can sit at sumw = 0.93 with nearly all the WORK left.  Kac
# weights by exactly the thing the cost is linear in.  §s270 (4) pins the total:
#
#     <N> = sum_c w_c N_c = 2*mu_d - 1,   mu_d = s274_mu_closed.mu_closed(P, Q)
#
# EXACT at even P (the beam sweeps the whole cross-section); at odd P the swept
# fraction g_sw < 1, so the same expression is an UPPER BOUND.  Scored against the
# store at s354: 49/50 even-P records exact to 1e-9, the sole exception (8,15) at
# a1 = 1 (degenerate, outside the cone); all 136 odd-P records at ratio <= 1.000.
#
# ⚠ SO `frac` IS A LOWER BOUND ON PROGRESS WHENEVER `exact` IS FALSE -- an odd-P
# centre can be finished at frac = 0.6.  Never read `frac < 1` as "not done";
# `complete` is the only completion flag.  Kac says how much work is LEFT, not
# whether the partition closed.
def kac_total(P, Q):
    """(value, exact) -- the width-weighted total word length `sum_c w_c N_c`."""
    from s274_mu_closed import mu_closed                       # probes/, lazy
    return 2.0 * mu_closed(P, Q) - 1.0, bool(P % 2 == 0 and Q // P >= 2)


def _kac_peeled(rec):
    return sum(c['width'] * c['N'] for c in rec.get('cells', [])
               if c.get('N') is not None and c.get('width') is not None)


def peel(rec):
    """The peel ledger for one record: how much Kac mass is resolved, and the
    estimated node cost of finishing.  Pure read -- no measurement."""
    tot, exact = kac_total(rec['P'], rec['Q'])
    peeled = _kac_peeled(rec)
    covered = sum(c['width'] for c in rec.get('cells', [])
                  if c.get('N') is not None)
    frac = (peeled / tot) if tot > 0 else None
    out = dict(total=tot, exact=exact, peeled=peeled, frac=frac,
               w_covered=covered, w_pending=rec.get('wpending'),
               n_pending=len(rec.get('pending') or []),
               complete=bool(rec.get('complete')))
    out.update(next_bound(rec, frac))
    return out


ESCALATE = 4.0        # how much bigger a blind next attempt should be


def next_bound(rec, frac=None):
    """Estimate the node cost of FINISHING, and the cap for the next attempt.

    ⚠ AN ESTIMATOR, NOT A LAW, AND EVERY NUMBER IT RETURNS IS A LOWER BOUND.  The
    peel rate degrades as you go: cost per unit Kac mass is 1/w and the intervals
    still pending are the narrow ones, so a rate measured on the cells resolved so
    far is optimistic about the ones that are not.

        rate  = (Kac mass peeled) / (nodes spent), on the last pass that MOVED it
        finish_nodes = nodes_cum + (total - peeled) / rate

    ⚠⚠ A POWER-LAW FIT IN THE UNPEELED FRACTION WAS TRIED FIRST AND IS WRONG TO
    USE -- not in principle (s340 G4's geometric letters motivate it) but
    NUMERICALLY: on real logs successive passes move the peel by ~0.005 while
    nodes move by ~3%, so `log(n1/n0)/log(f0/f1)` is a ratio of two near-zero
    quantities.  Measured at (16,53): theta = 159.1, then 74.8, then a
    divide-by-zero fallback, from three passes of the same enumerator on the same
    centre.  Do not reintroduce it without conditioning on the peel actually
    having moved by a decent factor.

    `stalled` is the case worth acting on: the last pass bought nodes and NO Kac
    mass, i.e. this cap is not the binding constraint and escalating it is not the
    move -- that centre needs a different tool, not a bigger budget.
    """
    if rec.get('complete'):
        return dict(next_nodes=None, basis='complete')
    log = [a for a in (rec.get('attempts') or []) if a.get('nodes_cum')]
    tot, _ = kac_total(rec['P'], rec['Q'])
    peeled = _kac_peeled(rec)
    if frac is None:
        frac = peeled / tot if tot > 0 else 0.0
    remaining = max(tot - peeled, 0.0)
    if not log:
        return dict(next_nodes=None, basis='no attempt log')
    ncum = log[-1]['nodes_cum']
    # the last pass on which the peel actually moved
    rate, basis, stalled = None, None, False
    for i in range(len(log) - 1, 0, -1):
        dn = log[i]['nodes_cum'] - log[i - 1]['nodes_cum']
        dk = (log[i].get('kac_peeled_cum') or 0) - (log[i - 1].get('kac_peeled_cum') or 0)
        if dn > 0 and dk > 0:
            rate = dk / dn
            basis = 'marginal rate (LOWER bound)'
            stalled = (i != len(log) - 1)
            break
    if rate is None:
        if ncum > 0 and peeled > 0:
            rate, basis = peeled / ncum, 'average rate, single pass (LOWER bound)'
            stalled = len(log) > 1
        else:
            return dict(next_nodes=int(ESCALATE * max(ncum, 1)), stalled=True,
                        basis='no peel yet -- escalate blind')
    finish = int(ncum + remaining / rate)
    return dict(next_nodes=int(min(finish, ESCALATE * ncum)) if not stalled
                else int(ESCALATE * ncum),
                finish_nodes=finish, rate=rate, stalled=stalled, basis=basis)


# --------------------------------------------------------------------------
class VerifyError(AssertionError):
    pass


def verify(rec, tol=1e-8):
    """The certificate, INSIDE the reader (E1).  Raises VerifyError.

    A COMPLETE record must be a partition: the widths sum to 1.  A PARTIAL one is
    checked against the accounting identity instead --

        sum(cell widths) + wpending + wnoexit + wsliver == 1

    -- so an incomplete record is no longer un-auditable.  Before s354 the two
    drain paths (`capped`, `noexit`) discarded their intervals' WIDTH as well as
    their state, so a lost interval was undetectable and `verify` had no choice
    but to skip such records entirely.  Measured exactly 0.0 residual on partial
    runs at (5,17), (7,31), (9,31), (3,23), (10,23), (3,47).
    """
    if rec.get('guard'):
        return rec
    P, Q, cells = rec['P'], rec['Q'], rec['cells']
    if rec['r'] != Q % P or rec['a1'] != Q // P:
        raise VerifyError(f'{P}/{Q}: r/a1 inconsistent with (P,Q)')
    # THE ACCRETION GUARD -- a field may not enter the store undeclared, so the
    # provenance schema cannot drift behind the data ([OPS-071]).
    if cells:
        undeclared = set().union(*(c.keys() for c in cells)) - set(CELL_FIELDS)
        if undeclared:
            raise VerifyError(
                f'{P}/{Q}: undeclared cell field(s) {sorted(undeclared)} -- add '
                f'them to `CELL_FIELDS` with their PROVENANCE (measured / '
                f'derived / joined / mixed) before writing them')
    s = sum(c['width'] for c in cells)
    if rec.get('complete', True):
        if abs(s - 1) > tol:
            raise VerifyError(
                f'{P}/{Q}: widths sum to {s!r}, not 1 -- this is a STACK or a '
                f'merge of two partitions, not a partition (see the module DO-NOT)')
    else:
        acct = s + sum(rec.get(k) or 0.0
                       for k in ('wpending', 'wnoexit', 'wsliver'))
        if abs(acct - 1) > tol:
            raise VerifyError(
                f'{P}/{Q}: PARTIAL record does not account for the transversal: '
                f'cells {s!r} + pending {rec.get("wpending")!r} + noexit '
                f'{rec.get("wnoexit")!r} + sliver {rec.get("wsliver")!r} = {acct!r}')
    for c in cells:
        if not 0 < c['width'] <= 1 + tol:
            raise VerifyError(f'{P}/{Q}: width {c["width"]} outside (0,1]')
    if rec.get('S') is not None:
        aw = [c['area_wT'] for c in cells if c['area_wT'] is not None]
        if len(aw) == len(cells) and abs(sum(aw) - rec['S']) > 1e-6:
            raise VerifyError(f'{P}/{Q}: sum w*T/chi = {sum(aw)!r} != S = {rec["S"]!r}')
    if not rec.get('complete', True) and rec.get('kac', {}).get('exact'):
        k = rec['kac']
        if k['peeled'] > k['total'] * (1 + 1e-6):
            raise VerifyError(f'{P}/{Q}: peeled Kac mass {k["peeled"]!r} exceeds the '
                              f'exact total {k["total"]!r}')
    return rec


# --------------------------------------------------------------------------
def get(P, Q, max_depth=200000, refresh=False, ring='cached', caps=None,
        budget=None):
    """Consult the store; measure and add if absent.  Always verified."""
    p = _path(P, Q)
    if not refresh and os.path.exists(p):
        rec = json.load(open(p))
        if rec.get('store_tag') == STORE_TAG:
            return verify(rec)
    rec = measure(P, Q, max_depth=max_depth, ring=ring, caps=caps, budget=budget)
    _write(rec)
    return rec


def _write(rec):
    """Verify, then write ATOMICALLY.

    ⚠⚠ THE WRITE MUST BE ATOMIC BECAUSE READING WHILE A DAEMON WRITES IS NOW
    ROUTINE (s355).  `json.dump` straight onto the target leaves a window in
    which the file on disk is a TRUNCATED record, so a concurrent reader gets a
    `JSONDecodeError` -- or, worse, a short read that happens to parse.  s355
    worked around it by snapshotting the whole store to /tmp before every read;
    `os.replace` is atomic on POSIX and DELETES that discipline, which is the
    point: a background `recap` and a foreground read no longer need to know
    about each other.
    """
    os.makedirs(STORE_DIR, exist_ok=True)
    verify(rec)
    p = _path(rec['P'], rec['Q'])
    tmp = f'{p}.tmp{os.getpid()}'
    with open(tmp, 'w') as fh:
        json.dump(rec, fh)
    os.replace(tmp, p)                       # atomic: readers see old or new
    return rec


def cells(P, Q, **kw):
    return get(P, Q, **kw)['cells']


def known_centres():
    """(P,Q) already on disk in the store."""
    out = []
    for f in glob.glob(os.path.join(STORE_DIR, '*.json')):
        try:
            P, Q = os.path.basename(f)[:-5].split('_')
            out.append((int(P), int(Q)))
        except ValueError:
            pass
    # ⚠ s357: sort the (P,Q) TUPLES, not the filenames.  Sorting the glob sorted
    # them as STRINGS, so every consumer walked the store as 10,11,...,19,2,20,...
    # -- `(11,115)` before `(11,23)`, all of P=10..23 before P=3.  Harmless for a
    # sweep that finishes; a SAMPLING BIAS for any pass that stops early, and
    # `fill_fanword`/`recap` are exactly that ([OPS-076]).
    return sorted(out)


def stats(coverage=False):
    """Store summary; `coverage=True` adds PER-FIELD fill rates.

    ⚠⚠ QUOTE THE COVERAGE WHENEVER YOU QUOTE A SCORE OFF A COLUMN (s355).  As
    fields accrete, partial coverage is the NORM, not the exception -- and a
    partially-filled column is a biased subsample, not merely a smaller one:
    `chi` is on ~32% of cells because only the RING-CERTIFIED centres have it,
    and those are the SMALL-`Q` ones.  A law scored on `chi` is therefore scored
    on small `Q` whether or not the write-up says so.  This used to be one
    hardcoded field (`cells_with_N`), which is what a missing generic mechanism
    looks like.
    """
    ks = known_centres()
    nc = ncert = ncomp = npart = nguard = 0
    filled = {k: 0 for k in CELL_FIELDS}
    for P, Q in ks:
        try:
            rec = json.load(open(_path(P, Q)))
        except Exception:                                      # noqa: BLE001
            continue
        cells = rec.get('cells', [])
        nc += len(cells)
        ncert += bool(rec.get('certified'))
        for c in cells:
            for k in filled:
                if c.get(k) is not None:
                    filled[k] += 1
        if rec.get('guard'):
            nguard += 1
        elif rec.get('complete', True):
            ncomp += 1
        else:
            npart += 1
    out = dict(centres=len(ks), complete=ncomp, partial=npart, guarded=nguard,
               cells=nc, certified_centres=ncert, cells_with_N=filled['N'],
               tag=STORE_TAG, dir=STORE_DIR)
    if coverage:
        out['coverage'] = {
            k: dict(n=filled[k], frac=round(filled[k] / nc, 4) if nc else None,
                    kind=CELL_FIELDS[k][0])
            for k in sorted(filled, key=lambda k: -filled[k])}
    return out


_HEAD_COLS = ['P', 'Q', 'r', 'a1', 'a2']
_TAIL_COLS = ['engine', 'certified', 'complete', 'n_cells', 'S',
              'tot_wN', 'tot_wT', 'sumw', 'kac_frac', 'kac_exact', 'kac_total',
              'wpending', 'rungs_present',
              'depth_read', 'capped', 'noexit', 'sliver', 'store_tag']
# Preferred DISPLAY order only.  The authoritative list is `CELL_FIELDS`:
# anything declared there and not named here is appended automatically, so
# adding a column to the schema exports it without editing this file.  ⚠ That
# is deliberate -- the hand-maintained twin of this list is exactly how
# `vertical` came to be declared but absent (s355), and `m_geom` declared but
# unexported (s356, caught by the 45-column count not moving).
_CELL_ORDER = ('m', 'm_geom', 'm_source', 'm_geom_M', 'cell_index', 'lo', 'hi',
               'width', 'N', 'T_over_chi', 'chi', 'T', 'area_wN', 'area_wT',
               'frac_wN', 'frac_wT', 'nfan', 'jstar', 'nopp', 'depth',
               'vertical', 'is_orphan', 'in_width_class', 'gb')
_SKIP_COLS = ('fan_chain',)                       # RESERVED (E1), null in v1
_CELL_COLS = tuple(dict.fromkeys(
    [f for f in _CELL_ORDER if f in CELL_FIELDS] +
    [f for f in CELL_FIELDS if f not in _SKIP_COLS]))
CSV_COLS = _HEAD_COLS + list(_CELL_COLS) + _TAIL_COLS


def export_csv(path='data/branch_widths.csv', centres=None, include_partial=False):
    """Flat one-row-per-branch export of the whole store.

    ⚠ PARTIAL records are EXCLUDED by default.  Their cells are correct, but they
    are not a partition of the transversal, and the commonest thing anyone does
    with this file is group by centre and sum -- which is the exact error the
    module DO-NOT is about.  Pass `include_partial=True` and read the `complete`
    column if you want them.
    """
    import csv
    centres = centres or known_centres()
    n = 0
    with open(path, 'w', newline='') as fh:
        wc = csv.DictWriter(fh, fieldnames=CSV_COLS)
        wc.writeheader()
        for P, Q in sorted(centres):
            rec = verify(json.load(open(_path(P, Q))))
            if rec.get('guard'):
                continue
            if not rec.get('complete', True) and not include_partial:
                continue
            head = {k: rec.get(k) for k in
                    ('P', 'Q', 'r', 'a1', 'a2', 'engine', 'certified', 'complete',
                     'n_cells', 'S', 'tot_wN', 'tot_wT', 'sumw', 'wpending',
                     'depth_read', 'capped', 'noexit', 'sliver', 'store_tag')}
            head['rungs_present'] = ' '.join(map(str, rec.get('rungs_present') or []))
            head['kac_frac'] = (rec.get('kac') or {}).get('frac')
            head['kac_exact'] = (rec.get('kac') or {}).get('exact')
            head['kac_total'] = (rec.get('kac') or {}).get('total')
            for c in rec['cells']:
                row = dict(head)
                row.update({k: c.get(k) for k in _CELL_COLS})
                wc.writerow({k: ('' if row.get(k) is None else row.get(k))
                             for k in CSV_COLS})
                n += 1
    return n


def centres_in_existing_stores():
    """(P,Q) any prior probe has already partitioned -- the backfill worklist."""
    found = set()
    for f in sorted(glob.glob(os.path.join(_DATA, '*.json'))):
        try:
            doc = json.load(open(f))
        except Exception:                                      # noqa: BLE001
            continue
        def walk(o, depth=0):
            if depth > 4:
                return
            if isinstance(o, dict):
                c = o.get('cells')
                if isinstance(c, list) and c and isinstance(c[0], dict):
                    P, Q = o.get('P'), o.get('Q')
                    if isinstance(P, int) and isinstance(Q, int) and P >= 2:
                        found.add((P, Q))
                    return
                for v in o.values():
                    walk(v, depth + 1)
            elif isinstance(o, list):
                for v in o:
                    walk(v, depth + 1)
        walk(doc)
    return sorted(found)


def backfill(centres=None, verbose=True, ring='cached', caps=None):
    """Seed the store from every centre a prior probe already partitioned."""
    centres = centres or centres_in_existing_stores()
    ok = part = fail = skip = 0
    for P, Q in centres:
        if os.path.exists(_path(P, Q)):
            skip += 1
            continue
        try:
            rec = get(P, Q, ring=ring, caps=caps)
            if rec.get('guard'):
                fail += 1
                if verbose:
                    print(f'  guard {P}/{Q}: {rec["guard"]}')
            elif not rec.get('complete', True):
                part += 1
            else:
                ok += 1
        except Exception as e:                                 # noqa: BLE001
            fail += 1
            if verbose:
                print(f'  FAIL  {P}/{Q}: {type(e).__name__}: {e}')
    return dict(ok=ok, partial=part, failed=fail, already=skip,
                requested=len(centres))


def incomplete_centres():
    """(P,Q) with unresolved intervals -- the recap worklist."""
    out = []
    for P, Q in known_centres():
        rec = json.load(open(_path(P, Q)))
        if not rec.get('complete', True) or rec.get('guard'):
            out.append((P, Q))
    return out


def recap_one(P, Q, caps=None, verbose=True):
    """Push ONE incomplete record further by RESUMING its pending intervals.

    This is the payoff of storing node state: the cost is proportional to what is
    still unresolved, not to the whole transversal.  A record is never made worse
    -- if the resumed pass resolves nothing new, the old record stands with one
    more row in its attempt log (which is itself worth having: two rows are what
    `next_bound` needs to fit an exponent).
    """
    from s340_glen_enum import enum_glen                       # probes/, lazy
    rec = json.load(open(_path(P, Q)))
    if rec.get('complete', True) and not rec.get('guard'):
        return rec, 'already complete'
    pend = rec.get('pending') or []
    if not pend:
        # a s353-era record, or one that died before any node was stored: there
        # is nothing to resume from, so this is a fresh measurement.
        rec = measure(P, Q, caps=caps)
        _write(rec)
        if verbose:
            _report(rec, 're-measured (no stored node state)')
        return rec, 'no pending state -> full re-measure'
    caps = dict(rec.get('caps') or CAPS, **(caps or {}))
    r, a1 = Q % P, Q // P
    t0 = time.time()
    gcells, gst = enum_glen(P, r, a1, resume=pend, collect_pending=True, **caps)
    secs = time.time() - t0

    # ⚠ `gb`/`hmax`/`e_end` come back None on a cell resolved from a PRE-s383 packed
    # node -- the winding history is not in the stored node state, and glen reports
    # the loss (`wok=False`) rather than restarting the walk from 0 mid-chain.  Those
    # cells are then filled by `refresh_glen` on a from-scratch pass, like any other
    # null.  `gst['nowind']` counts them.
    new = [dict(lo=_fnum(c['lo']), hi=_fnum(c['hi']), width=_fnum(c['w']),
                N=_fnum(c.get('closure')), T_over_chi=_fnum(c.get('fwd')),
                nfan=_fnum(c.get('nfan')), jstar=_fnum(c.get('jstar')),
                nopp=_fnum(c.get('nopp')), gb=_inum(c.get('gb')),
                hmax=_inum(c.get('hmax')), e_end=_inum(c.get('e_end')))
           for c in gcells]
    # DISJOINT concatenation, not a merge: the resumed pass only ever visits
    # intervals inside the pending nodes, which are disjoint from every resolved
    # cell.  (⚠ not the (closure, fwd) merge s340's NO MERGE STEP note forbids.)
    base = [dict(c) for c in rec['cells']] + new
    rec['complete'] = bool(gst.get('ok'))
    rec['caps'] = caps
    rec['pending'] = gst.get('pending') or []
    rec['pending_noexit'] = (rec.get('pending_noexit') or []) + \
                            (gst.get('noexit_nodes') or [])
    rec['wpending'] = _fnum(gst.get('wpending'))
    rec['wnoexit'] = (rec.get('wnoexit') or 0.0) + _fnum(gst.get('wnoexit') or 0.0)
    rec['wsliver'] = (rec.get('wsliver') or 0.0) + _fnum(gst.get('wsliver') or 0.0)
    rec['capped'] = gst.get('capped')
    rec['noexit'] = (rec.get('noexit') or 0) + int(gst.get('noexit') or 0)
    rec['sumw'] = sum(c['width'] for c in base)
    rec['S'] = (sum(c['width'] * c['T_over_chi'] for c in base)
                if all(c.get('T_over_chi') is not None for c in base) else None)
    rec['g_len'] = (rec['S'] / Q) if rec['S'] is not None else None
    if rec['complete']:
        rec['guard'] = None
    _derive(rec, base)
    rec['attempts'] = (rec.get('attempts') or []) + \
                      [_attempt(caps, gst, secs, rec, resumed=True)]
    rec['kac'] = peel(rec)
    _write(rec)
    if verbose:
        _report(rec, f'resumed, {secs:.0f}s')
    return rec, 'resumed'


def _report(rec, how):
    k = rec.get('kac') or {}
    tag = 'DONE ' if rec.get('complete') else ('STALL' if k.get('stalled') else '     ')
    print(f'  {tag} ({rec["P"]:3d},{rec["Q"]:4d}) a1={rec["a1"]:2d} '
          f'cells={len(rec.get("cells") or []):3d} sumw={rec.get("sumw") or 0:.9f} '
          f'peel={k.get("frac") or 0:.4f}{"" if k.get("exact") else "*"} '
          f'pending={k.get("n_pending")} finish~{k.get("finish_nodes")}  [{how}]',
          flush=True)


def _caps_exceed(new, old):
    return any(float(new.get(k) or 0) > float(old.get(k) or 0) for k in CAPS)


def recap(caps=None, centres=None, verbose=True, auto=False, max_budget=300.0):
    """Resume every incomplete record.  `caps` overrides the stored ones.

    `auto=True` sizes EACH centre from its OWN `next_bound` estimate instead of
    one flat cap for the sweep, and scales `budget` by the same factor (clamped to
    `max_budget`) so the clock does not bind before the raised node cap does --
    the failure mode this whole change is about.  Without `auto` the estimate is
    advisory and nobody acts on it, which is how a number becomes decoration.

    STALLED centres (a pass that bought nodes and no Kac mass) are SKIPPED unless
    the caps genuinely rise above the ones they stalled under -- re-running a
    stalled centre at the same caps is guaranteed to buy nothing, and a sweep that
    does it looks like work.
    """
    centres = centres or incomplete_centres()
    done = moved = skipped = 0
    for P, Q in centres:
        rec = json.load(open(_path(P, Q)))
        old = rec.get('caps') or dict(CAPS)
        want = dict(old, **(caps or {}))
        k = rec.get('kac') or {}
        if auto and k.get('next_nodes'):
            grow = max(1.0, float(k['next_nodes']) / max(float(old.get('max_nodes') or 1), 1.0))
            want['max_nodes'] = max(int(k['next_nodes']), int(want.get('max_nodes') or 0))
            want['budget'] = min(max_budget, (old.get('budget') or CAPS['budget']) * grow)
            want['max_fans'] = max(want.get('max_fans') or 0, want['max_nodes'] * 50)
            want['max_letters'] = max(want.get('max_letters') or 0, want['max_nodes'] * 1000)
        if k.get('stalled') and not _caps_exceed(want, old):
            skipped += 1
            if verbose:
                print(f'  SKIP  ({P:3d},{Q:4d}) stalled at these caps -- raise them '
                      f'or use a different tool', flush=True)
            continue
        before = rec.get('sumw') or 0.0
        rec, _how = recap_one(P, Q, caps=want, verbose=verbose)
        if rec.get('complete'):
            done += 1
        elif (rec.get('sumw') or 0.0) > before + 1e-12:
            moved += 1
    return dict(requested=len(centres), completed=done, advanced=moved,
                skipped_stalled=skipped)


def migrate(verbose=True):
    """Upgrade s353-era records to the current tag WITHOUT re-measuring.

    A complete s353 record is already correct -- it just lacks the new fields --
    so re-running the enumerator on 186 good centres to add a flag would be pure
    waste.  Guarded s353 records carry no cells and no node state, so they are
    left for `recap_one` to re-measure from the root.
    """
    up = left = 0
    for P, Q in known_centres():
        rec = json.load(open(_path(P, Q)))
        if rec.get('store_tag') == STORE_TAG:
            continue
        if rec.get('guard') or not rec.get('cells'):
            rec.setdefault('complete', False)
            rec.setdefault('pending', [])
            rec.setdefault('pending_noexit', [])
            rec.setdefault('attempts', [])
            left += 1
        else:
            rec['complete'] = True
            rec['pending'], rec['pending_noexit'] = [], []
            rec['wpending'] = 0.0
            rec.setdefault('wnoexit', 0.0)
            rec.setdefault('caps', dict(CAPS))
            rec.setdefault('attempts', [])
            rec['kac'] = peel(rec)
            up += 1
        rec['store_tag'] = STORE_TAG
        json.dump(rec, open(_path(P, Q), 'w'))
    if verbose:
        print(f'migrate: {up} complete records upgraded in place, '
              f'{left} incomplete left for recap()')
    return up, left


def rederive(verbose=True):
    """Re-run `_derive` on every stored record IN PLACE -- no re-measurement.

    The right primitive whenever the schema gains a DERIVED column (one computed
    from fields already stored): re-measuring 237 centres to add an arithmetic
    ratio would be pure waste, and the same argument as `migrate()`.  Subsumes
    `reindex_m` (which re-runs the same rung-index logic and is kept because it
    is cited).  ⚠ Only for derived fields -- a change in ENGINE semantics needs a
    `STORE_TAG` bump and a real re-measure.
    ⚠ It is a READ-MODIFY-WRITE.  Atomic `_write` protects a concurrent READER
    but NOT against a LOST UPDATE: if a `recap` daemon advances a record between
    the read and the write, that advance is discarded.  Run it when nothing else
    is writing.  (It used a raw, non-atomic `json.dump` until s356.)
    """
    n = 0
    for P, Q in known_centres():
        rec = json.load(open(_path(P, Q)))
        if rec.get('guard') or not rec.get('cells'):
            continue
        _derive(rec, rec['cells'])
        rec['kac'] = peel(rec)
        _write(rec)
        n += 1
    if verbose:
        print(f'rederive: {n} records re-derived in place')
    return n


def refresh_ring(centres=None, verbose=True, partial=False):
    """Re-merge the RING columns into existing records -- no glen, no re-measure.

    The third member of the family, and the one the s355 coverage report is what
    surfaced: `migrate()` upgrades a record's schema, `rederive()` recomputes
    DERIVED columns, and this refreshes columns MEASURED BY THE OTHER ENGINE.
    It exists because `_merge_ring` runs only inside `measure()`, so a record
    backfilled before a ring column existed never acquires it and `recap_one`
    (which resumes glen only) never adds it either -- `vertical` was on
    `0/4967` cells for exactly this reason, and `nopp` was missing from every
    pre-s354 record.

    ⚠ `ring_cache.load` is a DISK READ (`ring='cached'` semantics): a centre
    whose ring row was never computed is skipped, not computed.  So this is
    seconds for the whole store, and the expensive tail stays a daemon job.
    """
    import ring_cache
    n_seen = n_touched = n_skipped = 0
    for P, Q in (centres or known_centres()):
        # an explicit `centres` list may name a centre with NO store record;
        # `known_centres()` cannot, which is why this path went unexercised until
        # s379 (`s378_ring_partial` queues from the RING guards, not the store, and
        # died on `(5,54)`).  Skip, do not raise.
        if not os.path.exists(_path(P, Q)):
            continue
        rec = json.load(open(_path(P, Q)))
        if rec.get('guard') or not rec.get('cells'):
            continue
        n_seen += 1
        depth = rec.get('depth_read') or 200000
        if partial:
            got = ring_cache.load(P, Q, depth, partial=True)
            if got is None:
                continue
            _ctx, rrows, _rpend, rguard = got
        else:
            got = ring_cache.load(P, Q, depth)
            if got is None:
                continue
            _ctx, rrows, rguard = got
            if rguard:
                continue
        if not rrows:
            continue
        ring_complete = rguard is None
        before = sum(1 for c in rec['cells'] for k in ('chi', 'T', 'depth',
                                                       'vertical')
                     if c.get(k) is not None)
        base = [dict(c) for c in rec['cells']]
        _merge_ring(rec, base, rrows, ring_complete=ring_complete)
        _derive(rec, base)
        after = sum(1 for c in rec['cells'] for k in ('chi', 'T', 'depth',
                                                      'vertical')
                    if c.get(k) is not None)
        if after == before:
            continue                        # nothing gained -- do NOT rewrite
        # OPTIMISTIC CONCURRENCY.  Atomic writes stop a reader seeing a torn
        # file; they do NOT stop a read-modify-write from clobbering an update
        # a background `recap` made in between.  With detached fills now routine
        # (s355), re-read and skip if the record moved under us.
        live = json.load(open(_path(P, Q)))
        if (len(live.get('cells') or []) != len(rec['cells'])
                or len(live.get('attempts') or []) != len(rec.get('attempts') or [])):
            n_skipped += 1
            if verbose:
                print(f'  ({P:3d},{Q:4d}) SKIP -- changed under us (a recap is '
                      f'running); re-run refresh_ring afterwards', flush=True)
            continue
        _write(rec)
        n_touched += 1
        if verbose:
            print(f'  ({P:3d},{Q:4d}) ring fields {before} -> {after}', flush=True)
    if verbose:
        print(f'refresh_ring: {n_seen} read, {n_touched} gained fields, '
              f'{n_skipped} skipped (changed under us)')
    return dict(seen=n_seen, touched=n_touched, skipped=n_skipped)


def _merge_graze(base, gmap, tol=1e-8):
    """Write the graze columns onto `base` from a height -> (V, k) map.

    ⇒⇒ THE GRAZE PAIR IS NOT A THIRD MEASUREMENT ([OPS-114], user-raised s380: *"there
    should be only 2 options - fast vs certified.  each should measure the same
    thing"*).  A cell boundary IS a vertex graze (`foundations.md` sec.1), so the
    engine that decided where to cut necessarily knew which vertex -- `enum_glen`
    picks the RADIUS at its split (`rho(i) = cot` at an R copy, `csc` at an A copy)
    and its two centre-splits are O copies.  It was computing the pair and dropping
    it.  This merges it onto the row on the SAME pass, so `graze_*` arrives with `N`
    instead of needing a separate third-engine backfill.

    Returns the number of boundaries matched.  A boundary with no match is left None
    -- never guessed.  `lo == 0` / `hi == 1` are transversal ends, not grazes.
    """
    if not gmap:
        return 0
    keys = sorted(gmap)

    def look(y, end):
        if end:
            return None, None
        hit = min(keys, key=lambda z: abs(z - y))
        if abs(hit - y) > tol:
            return None, None
        return gmap[hit]

    n = 0
    for c in sorted(base, key=lambda d: d['lo']):
        c['graze_lo_V'], c['graze_lo_k'] = look(c['lo'], c['lo'] < 1e-12)
        c['graze_hi_V'], c['graze_hi_k'] = look(c['hi'], c['hi'] > 1 - 1e-12)
        n += (c['graze_lo_V'] is not None) + (c['graze_hi_V'] is not None)
    return n


def fill_graze(centres=None, verbose=True, deadline=None,
               stop_file='data/s373_graze.stop', caps=None, budget=300.0):
    """Backfill the VERTEX-GRAZE columns `graze_{lo,hi}_{V,k}` on EXISTING rows.

    ⇒⇒ REWRITTEN s380: THE SOURCE IS NOW **glen**, THE FAST ENGINE, NOT A THIRD ONE.
    Until s380 this ran `s196_ringwords.enumerate_branches` -- its own docstring
    called that *"a THIRD ENGINE ON THIS RECORD"* -- to recover a pair that BOTH
    store engines already compute and discard.  Three consequences, all measured:
      * it was CAP-COUPLED: its `max_depth=400000` default ran against `ring_cache`
        rows built at `600000`, so `enumerate_branches` came back `unclosed > 0` and
        every centre was skipped as "tree uncertified" -- **0 filled / 22 skipped**,
        and silently, in ~90 min ([OPS-113]).  glen has **NO depth cap at all**, so
        that failure mode does not exist for it.
      * it was SLOW: ~3.5 min per centre of exact ring tree, against glen's ~1 s.
      * it was a per-COLUMN pass, which is why coverage was per-column (`graze` 62%)
        instead of per-row.
    The switch is gated on `probes/s380_graze_source.py`: glen vs s196 agree on the
    graze VERTEX **3555/3555** and on the DEPTH `k` **3555/3555**, 0 boundaries
    uncovered, over the 207 centres the old source had filled.  ⚠ That control CAN
    fail and did before `enum_glen._graze` kept the SHALLOWEST `k` ([OPS-097]).

    ⚠ COMPLETE records only -- unchanged.  A partial record is not a partition, so
    its cell list and a fresh run's need not correspond.  The 36 partial centres stay
    ungrazed until they complete.
    ⚠ CONTROLLED, not trusted -- also unchanged: the centre is SKIPPED unless glen
    reports `ok`, finds the same number of cells as the store, and matches every
    stored boundary to 1e-8.
    ⚠ A skip PRINTS.  s380 read 22 consecutive silent skips as one slow centre.
    ⚠ CAPS COME FROM THE RECORD, not from a default -- the same lesson as [OPS-113]'s
    `max_depth`.  A centre completed over several `recap` rounds stores the caps that
    finished it (`(3,23)`: three 45 s rounds), so a fresh run at library defaults comes
    back PARTIAL and skips.  `budget` is a wall clock, not a semantic cap, so it is the
    one knob raised here; the semantic caps are the record's own.
    ⚠ Cost class (ii) -- an engine already emits it and the store dropped it.  Schema
    addition, **no `STORE_TAG` bump**.
    """
    import time
    from s340_glen_enum import enum_glen                        # probes/, lazy
    todo = []
    for P, Q in (centres or known_centres()):
        if not os.path.exists(_path(P, Q)):
            continue                                   # [OPS-110]: skip, do not raise
        rec = json.load(open(_path(P, Q)))
        if rec.get('guard') or not rec.get('complete') or not rec.get('cells'):
            continue
        if all(c.get('graze_lo_V') is not None or c['lo'] < 1e-12
               for c in rec['cells']) and \
           all(c.get('graze_hi_V') is not None or c['hi'] > 1 - 1e-12
               for c in rec['cells']):
            continue                                   # already filled
        cost = sum(c['N'] for c in rec['cells'] if c.get('N') is not None)
        todo.append((cost, P, Q))
    todo.sort()                                        # cheapest first ([OPS-113])
    if verbose:
        print(f'fill_graze: {len(todo)} centres queued (glen source), cheapest first; '
              f'total sum-N {sum(t[0] for t in todo):,.0f}', flush=True)

    t0 = time.time()
    n_ok = n_skip = 0
    reasons = collections.Counter()
    timings = []
    for k, (_c, P, Q) in enumerate(todo, 1):
        if stop_file and os.path.exists(stop_file):
            if verbose:
                print(f'  STOP FILE -- halting after {k - 1} centres', flush=True)
            break
        if deadline is not None and time.time() - t0 > deadline:
            if verbose:
                print(f'  DEADLINE -- halting after {k - 1} centres', flush=True)
            break

        def _skip(why):
            nonlocal n_skip
            n_skip += 1
            reasons[why] += 1
            if verbose:
                print(f'  [{k:3d}/{len(todo)}] ({P:3d},{Q:4d}) SKIP -- {why}', flush=True)

        rec = json.load(open(_path(P, Q)))
        cells = rec['cells']
        r, a1 = Q % P, Q // P
        if r == 0 or a1 < 1:
            _skip('glen does not take this centre'); continue
        t_c = time.time()
        cc = dict(rec.get('caps') or CAPS)
        cc['budget'] = budget
        cc.update(caps or {})
        try:
            gcells, gst = enum_glen(P, r, a1, **cc)
        except Exception as e:                                 # noqa: BLE001
            _skip(f'raised:{type(e).__name__}'); continue
        if not gst.get('ok'):
            _skip('glen partial'); continue
        if len(gcells) != len(cells):
            _skip(f'n_cells disagree ({len(gcells)} vs {len(cells)})'); continue
        gmap = dict(enum_glen.last_grazes)

        live = json.load(open(_path(P, Q)))                    # lost-update guard
        if len(live.get('cells') or []) != len(cells):
            _skip('changed under us'); continue
        probe = [dict(lo=c['lo'], hi=c['hi']) for c in live['cells']]
        n_match = _merge_graze(probe, gmap)
        want = sum((c['lo'] > 1e-12) + (c['hi'] < 1 - 1e-12) for c in live['cells'])
        if n_match != want:
            _skip(f'boundary unmatched ({n_match}/{want})'); continue
        for c, p in zip(sorted(live['cells'], key=lambda d: d['lo']), probe):
            for f in ('graze_lo_V', 'graze_lo_k', 'graze_hi_V', 'graze_hi_k'):
                c[f] = p[f]
        _write(live)
        n_ok += 1
        dt = time.time() - t_c
        timings.append(dict(P=P, Q=Q, secs=round(dt, 2), n_cells=len(cells)))
        if verbose:
            print(f'  [{k:3d}/{len(todo)}] ({P:3d},{Q:4d}) {len(cells):4d} cells '
                  f'filled  {dt:7.1f}s', flush=True)
    if verbose:
        print(f'fill_graze: {n_ok} centres filled, {n_skip} skipped '
              f'{dict(reasons)} ({len(todo)} queued)')
    return dict(filled=n_ok, skipped=n_skip, queued=len(todo),
                reasons=dict(reasons), timings=timings)



def refresh_glen(centres=None, verbose=True, deadline=None,
                 stop_file='data/s380_glen.stop', caps=None, budget=300.0,
                 tol=1e-6, force=False):
    """Re-merge EVERY glen-sourced column from ONE glen pass -- the glen twin of
    `refresh_ring()`, and the general form of what `fill_graze` does for one column.

    ⇒⇒ WHY (user-raised s380, [OPS-114]): *"there should be only 2 options - fast vs
    certified. each should measure the same thing."*  The store's header says two
    engines; the reality was SIX writer paths, several of them one-column recovery
    passes over data a single glen run already produces.  This is the fast engine's
    ONE writer: it re-runs glen once and merges `N`, `T_over_chi`, `nfan`, `jstar`,
    `nopp`, `graze_{lo,hi}_{V,k}` AND (s383) the winding trio `gb`/`hmax`/`e_end`
    together.  A column left null by an older schema (`nopp` 42.9%, `graze_*` before
    s380, `gb`/`hmax`/`e_end` 86.7% before s383) is filled on the same pass that
    re-confirms the ones already there.

    ⚠⚠ CONTROLLED, NOT TRUSTED, AND THE CONTROL IS THE POINT: a centre is SKIPPED
    unless glen reports `ok`, finds the same `n_cells`, matches every stored boundary
    to 1e-8, AND agrees with every EXISTING non-null glen value.  A disagreement is a
    real finding -- it means the record and the engine have drifted -- so it is
    reported per centre and never written over.
    ⇒ THE AGREEMENTS ARE COUNTED TOO (`confirmed`, s383), not just the drifts.  A
    zero-drift report alone does not say whether ANYTHING was compared: a column at
    0% coverage has no stored value to disagree with, so it reads exactly like a
    column that was checked and passed.  `confirmed[f]` is the denominator that tells
    the two apart -- and when a column's SOURCE changes (`graze_*` s380, the winding
    trio s383) that denominator IS the source gate, run at full store scale.

    `force=True` (s383) queues a centre that needs NOTHING, so the pass is a pure
    re-confirmation sweep.  That is how a SOURCE CHANGE is gated at full store scale:
    the default `need` filter visits only centres with a null, which on a
    nearly-complete column is a biased handful and leaves the already-populated rows
    -- the ones written by the OLD source -- unchecked.

    ⚠ Caps come from the RECORD ([OPS-113]); only `budget`, a wall clock, is raised.
    ⚠ COMPLETE records only -- a partial record is not a partition.
    ⚠ Cost class (ii): re-measure, schema fill, **no `STORE_TAG` bump**.
    ⚠ Lost-update guard: the record is re-read immediately before the write.
    """
    import time
    from s340_glen_enum import enum_glen
    GCOLS = ('N', 'T_over_chi', 'nfan', 'jstar', 'nopp')
    # the winding trio (s383): same pass, but compared EXACTLY -- they are integers.
    WCOLS = ('gb', 'hmax', 'e_end')

    todo = []
    for P, Q in (centres or known_centres()):
        if not os.path.exists(_path(P, Q)):
            continue
        rec = json.load(open(_path(P, Q)))
        if rec.get('guard') or not rec.get('complete') or not rec.get('cells'):
            continue
        cells = rec['cells']
        need = any(c.get(f) is None for c in cells for f in GCOLS + WCOLS) or \
            any(c.get('graze_lo_V') is None and c['lo'] > 1e-12 for c in cells) or \
            any(c.get('graze_hi_V') is None and c['hi'] < 1 - 1e-12 for c in cells)
        if not (need or force):
            continue
        todo.append((sum(c['N'] for c in cells if c.get('N') is not None), P, Q))
    todo.sort()
    if verbose:
        print(f'refresh_glen: {len(todo)} centres queued, cheapest first', flush=True)

    t0 = time.time()
    n_ok = n_skip = 0
    filled = collections.Counter()
    confirmed = collections.Counter()
    reasons = collections.Counter()
    drift = []
    for k, (_c, P, Q) in enumerate(todo, 1):
        if stop_file and os.path.exists(stop_file):
            if verbose:
                print(f'  STOP FILE -- halting after {k - 1} centres', flush=True)
            break
        if deadline is not None and time.time() - t0 > deadline:
            if verbose:
                print(f'  DEADLINE -- halting after {k - 1} centres', flush=True)
            break

        def _skip(why):
            nonlocal n_skip
            n_skip += 1
            reasons[why] += 1
            if verbose:
                print(f'  [{k:3d}/{len(todo)}] ({P:3d},{Q:4d}) SKIP -- {why}', flush=True)

        rec = json.load(open(_path(P, Q)))
        cells = sorted(rec['cells'], key=lambda d: d['lo'])
        r, a1 = Q % P, Q // P
        if r == 0 or a1 < 1:
            _skip('glen does not take this centre'); continue
        cc = dict(rec.get('caps') or CAPS)
        cc['budget'] = budget
        cc.update(caps or {})
        t_c = time.time()
        try:
            gcells, gst = enum_glen(P, r, a1, **cc)
        except Exception as e:                                 # noqa: BLE001
            _skip(f'raised:{type(e).__name__}'); continue
        if not gst.get('ok'):
            _skip('glen partial'); continue
        if len(gcells) != len(cells):
            _skip(f'n_cells disagree ({len(gcells)} vs {len(cells)})'); continue
        gmap = dict(enum_glen.last_grazes)
        gcells = sorted(gcells, key=lambda d: float(d['lo']))

        # position match + AGREEMENT on everything already stored
        bad = None
        new = []
        conf = collections.Counter()      # per-centre; banked only if it WRITES
        for c, g in zip(cells, gcells):
            if abs(_fnum(g['lo']) - c['lo']) > 1e-8 or abs(_fnum(g['hi']) - c['hi']) > 1e-8:
                bad = 'boundary unmatched'; break
            got = dict(N=_fnum(g.get('closure')), T_over_chi=_fnum(g.get('fwd')),
                       nfan=_fnum(g.get('nfan')), jstar=_fnum(g.get('jstar')),
                       nopp=_fnum(g.get('nopp')), gb=_inum(g.get('gb')),
                       hmax=_inum(g.get('hmax')), e_end=_inum(g.get('e_end')))
            for f in GCOLS:
                have = c.get(f)
                if have is not None and got[f] is not None and \
                        abs(float(have) - float(got[f])) > tol * max(1.0, abs(float(have))):
                    bad = f'{f} drift ({have} vs {got[f]})'
                    drift.append(dict(P=P, Q=Q, field=f, stored=have, glen=got[f]))
                    break
                if have is not None and got[f] is not None:
                    conf[f] += 1
            # ⚠ the winding trio is EXACT: any disagreement with what
            # `fill_fanword`'s pointwise tracer wrote is a real finding, not a
            # tolerance question, so it skips the centre like any other drift.
            for f in WCOLS:
                have = c.get(f)
                if have is not None and got[f] is not None and int(have) != got[f]:
                    bad = f'{f} drift ({have} vs {got[f]})'
                    drift.append(dict(P=P, Q=Q, field=f, stored=have, glen=got[f]))
                    break
                if have is not None and got[f] is not None:
                    conf[f] += 1
            if bad:
                break
            new.append(got)
        if bad:
            _skip(bad); continue

        probe = [dict(lo=c['lo'], hi=c['hi']) for c in cells]
        _merge_graze(probe, gmap)

        live = json.load(open(_path(P, Q)))                    # lost-update guard
        if len(live.get('cells') or []) != len(cells):
            _skip('changed under us'); continue
        for c, got, pz in zip(sorted(live['cells'], key=lambda d: d['lo']), new, probe):
            for f in GCOLS + WCOLS:
                if c.get(f) is None and got[f] is not None:
                    c[f] = got[f]
                    filled[f] += 1
            for f in ('graze_lo_V', 'graze_lo_k', 'graze_hi_V', 'graze_hi_k'):
                if c.get(f) is None and pz.get(f) is not None:
                    c[f] = pz[f]
                    filled[f] += 1
        _derive(live, live['cells'])
        _write(live)
        n_ok += 1
        confirmed += conf
        if verbose:
            print(f'  [{k:3d}/{len(todo)}] ({P:3d},{Q:4d}) ok  '
                  f'{time.time() - t_c:6.1f}s', flush=True)
    if verbose:
        print(f'refresh_glen: {n_ok} centres refreshed, {n_skip} skipped '
              f'{dict(reasons)}; cells filled {dict(filled)}; '
              f'cells RE-CONFIRMED {dict(confirmed)}; {len(drift)} drifted')
    return dict(refreshed=n_ok, skipped=n_skip, queued=len(todo),
                filled=dict(filled), confirmed=dict(confirmed),
                reasons=dict(reasons), drift=drift)


def fill_rank(cells):
    """Rank a centre's cells for a RESUMABLE fill: 0 = orphan, 1 = one of each
    J-pair, 2 = the twin.  Returns a list of ranks parallel to `cells`.

    ⚠⚠ WHY A FILL IS ORDERED AT ALL ([OPS-076]).  A pass that can stop early
    (deadline, cap, kill) leaves a PARTIAL column, and the iteration order decides
    what that column is a sample OF.  Two orders are actively wrong here:
      * BY COST.  Cost is linear in `nfan`, and `nfan` is the variable `gb`/`hmax`
        exist to be scored against, so cheapest-first makes coverage a function of
        the DEPENDENT VARIABLE.  It would have deleted s357's own headline, which
        is a statement about the four LONGEST chains (`G_B = 12` at `nfan = 58724`).
      * WHOLE CENTRE AT A TIME, in store order.  That is a PREFIX, not a sample:
        "100% of the first k centres" rather than "x% of each".
    So: round-robin by RANK across all centres (see `fill_fanword`), and inside a
    centre put the row that carries information first and the redundant one last.
      rank 0, THE ORPHAN -- the only row whose `e_end` differs at all: the closure
        is a WRAP (`|E| = Q`) exactly on the orphan and an exact cancellation
        (`E = 0`) on every J-paired cell ([NCYL-188], s357).  Identified by `chi`
        where the ring has spoken, else by ODD `nfan` ([NCYL-182]'s parity law,
        `728/728`); the two agreed on every cell where both were available.
      rank 1/2, THE J-PAIR -- adjacent, equal `(width, N, nfan)` (Thm D +
        `orphan_theorem.md` Lemma B'; s356 (8), `198/198` centres).  The twin is
        pure redundancy for THESE columns: `gb`, `hmax` and `e_end` agree within
        the pair on `1097/1097` measured pairs.  It is the Lemma-B' control, so it
        is worth having -- and worth having LAST.
    """
    order = sorted(range(len(cells)), key=lambda i: cells[i]['lo'])
    rank = [1] * len(cells)
    used = set()
    for i in order:
        c = cells[i]
        orph = c.get('is_orphan')
        if orph is None and c.get('nfan') is not None:
            orph = bool(int(c['nfan']) % 2)
        if orph:
            rank[i] = 0
            used.add(i)
    prev = None
    for i in order:
        if i in used:
            prev = None
            continue
        c = cells[i]
        if prev is not None:
            a = cells[prev]
            if (a.get('nfan') == c.get('nfan') and a.get('N') == c.get('N')
                    and abs((a.get('width') or 0) - (c.get('width') or 0)) < 1e-12):
                rank[i] = 2                     # the twin goes last
                prev = None
                continue
        prev = i
    return rank


def fill_fanword(centres=None, verbose=True, max_nfan=60000, deadline=None,
                 stop_file='data/s357_fanword.stop', ranks=(0, 1, 2)):
    """Populate the WINDING columns `gb`/`hmax`/`e_end` from the POINTWISE tracer.

    ⇒⇒ ⚠⚠ **SUPERSEDED s383 AS THE ROUTINE WRITER -- reach for `refresh_glen()`.**
    This docstring already conceded the case against itself ([OPS-114] quoted it):
    *"the SAME engine that emitted `nfan`, reading three further observables off the
    same walk"*.  If it is the same walk, the store should not re-run it -- s383 made
    `enum_glen` carry the winding as three integer passengers, so `gb`/`hmax`/`e_end`
    now arrive with `N` on `measure()`'s own pass and backfill through the fast
    engine's ONE writer, `refresh_glen`.  Measured consequences, all three of which
    this function's design forced:
      * NO `max_nfan` CAP.  The trio was capped at `nfan <= 60000` because the tracer
        holds `vs`/`steps` as LISTS (~100 B/fan); the ride-along keeps three scalars,
        so the long chains -- the rows [OPS-076] warns the cap silently removes -- are
        no longer a truncated class.
      * NO `a1 >= 2` BAR.  That bar is inherited from `s356_rung_exact.scored_centres`'
        RUNG cone and is nothing to do with the winding; it left 712 cells over 32
        ring-certified `a1 = 1` centres permanently null.
      * NO SECOND PASS.  The cost of the trio is now zero: it rides on a walk that was
        happening anyway.
    ⇒ KEPT, and not deprecated to a stub, because it is the INDEPENDENT ARM OF THE
    GATE: `probes/s383_fanword_source.py` scores the ride-along against this tracer
    and that control can only run while both exist.

    Source: `s357_fanword.trace_fanword`, a POINTWISE specialisation of
    `s340_glen_enum`'s node loop -- so this is the SAME engine that emitted `nfan`,
    reading three further observables off the same walk, not a second model.  It is
    controlled against the store on `nfan`/`N`/`jstar`/`nopp` cell by cell and a
    mismatch SKIPS the cell rather than writing it (s357: 1623/1623 agreed).

    ⇒⇒ ORDER: RANK ROUNDS, AND CHEAPEST-CENTRE-FIRST INSIDE EACH ROUND ([OPS-076]).
    The two are orthogonal and both are earned:
      * The ROUNDS decide WHICH KIND of row -- 0 every centre's ORPHAN, 1 one cell of
        each J-pair, 2 the twins -- so an interrupted fill leaves a STRATIFIED sample
        rather than a prefix of the store.  `fill_rank` says why those three.
      * CHEAPEST-FIRST inside a round decides which centre, and it is nearly free
        coverage: tracing cost is linear in `nfan` and the cost distribution is
        savagely skewed -- **1% of the total buys 70% of centres and 52% of cells;
        5 centres of 237 carry 37% of it.**
    ⚠ The bias cheapest-first creates is REAL and must be quoted with any partial
    score: tau(centre cost, max `nfan`) = **0.907**, so an interrupted pass holds the
    SHORT chains.  (For `gb` itself it is mild, tau = 0.244.)  A round that RUNS TO
    COMPLETION has no order bias at all -- which is the usual case, since rank 0 is
    one cell per centre.

    ⚠ Cost class (iii), ~1 s per 1000 fans.  A cell over `max_nfan` is left null AND
    ITS `nfan` IS REPORTED (`capped_nfan`), because the truncation is a bias on the
    very column being scored: "n dropped" does not say the dropped rows were the long
    chains.  ⚠ Read-modify-write: it re-reads and skips a record that moved under it
    ([OPS-074]).  Resumable -- a cell with `gb` already set is skipped.
    """
    import mpmath as mp
    from s357_fanword import trace_fanword
    # --- PLAN.  One read per record, then work sorted by (rank, cost).
    plan = []
    capped_nfan = []
    for P, Q in (centres or known_centres()):
        # an explicit `centres` list may name a centre with NO store record;
        # `known_centres()` cannot, which is why this path went unexercised until
        # s379 (`s378_ring_partial` queues from the RING guards, not the store, and
        # died on `(5,54)`).  Skip, do not raise.
        if not os.path.exists(_path(P, Q)):
            continue
        rec = json.load(open(_path(P, Q)))
        if rec.get('guard') or not rec.get('cells'):
            continue
        if rec['a1'] < 2 or Q % P == 0:
            continue
        rk = fill_rank(rec['cells'])
        todo = collections.defaultdict(list)
        for i, (c, this) in enumerate(zip(rec['cells'], rk)):
            if c.get('gb') is not None or c.get('nfan') is None:
                continue
            if c['nfan'] > max_nfan:
                capped_nfan.append(int(c['nfan']))
                continue
            todo[this].append(i)
        for this, idx in todo.items():
            cost = sum(int(rec['cells'][i]['nfan']) for i in idx)
            plan.append((this, cost, P, Q, idx))
    plan.sort(key=lambda z: (ranks.index(z[0]) if z[0] in ranks else 99, z[1]))
    plan = [z for z in plan if z[0] in ranks]
    if verbose:
        print(f'plan: {len(plan)} (centre, rank) units, '
              f'{sum(len(z[4]) for z in plan)} cells, '
              f'{sum(z[1] for z in plan):,} fans to trace', flush=True)

    n_touched = n_cells = n_skip = n_bad = 0
    seen, rank_now = set(), None
    for rank, cost, P, Q, idx in plan:
        if stop_file and os.path.exists(stop_file):
            if verbose:
                print(f'  stop file {stop_file} -- halting', flush=True)
            break
        if deadline and time.time() > deadline:
            if verbose:
                print('  deadline -- halting', flush=True)
            break
        if rank != rank_now and verbose:
            rank_now = rank
            print(f'--- rank {rank} '
                  f'({ {0: "orphans", 1: "one per J-pair", 2: "the twins"}[rank] }),'
                  f' cheapest centre first', flush=True)
        rec = json.load(open(_path(P, Q)))
        a1, r = rec['a1'], Q % P
        seen.add((P, Q))
        got = 0
        for i in idx:
            if i >= len(rec['cells']):
                break
            c = rec['cells'][i]
            if c.get('gb') is not None or c.get('nfan') is None:
                continue
            s_mid = (mp.mpf(repr(c['lo'])) + mp.mpf(repr(c['hi']))) / 2
            t = trace_fanword(P, r, a1, s_mid, max_fans=max_nfan)
            if not t.get('closed'):
                capped_nfan.append(int(c['nfan']))
                continue
            # CONTROL -- the trace must reproduce what the store already holds.
            if (t['nfan'] != c['nfan'] or t['N'] != c['N']
                    or t['jstar'] != c['jstar']
                    or (c.get('nopp') is not None and t['nopp'] != c['nopp'])):
                n_bad += 1
                continue
            vs, steps = t['vs'], t['steps']
            e_end = vs[-1] + steps[-1]
            c['gb'] = max(max(abs(v) for v in vs), abs(e_end))
            c['hmax'] = max(abs(g) for g in steps)
            c['e_end'] = e_end
            got += 1
        if not got:
            continue
        live = json.load(open(_path(P, Q)))
        if (len(live.get('cells') or []) != len(rec['cells'])
                or len(live.get('attempts') or []) != len(rec.get('attempts') or [])):
            n_skip += 1
            if verbose:
                print(f'  ({P:3d},{Q:4d}) SKIP -- changed under us', flush=True)
            continue
        _write(rec)
        n_touched += 1
        n_cells += got
        if verbose:
            print(f'  ({P:3d},{Q:4d}) +{got}', flush=True)
    if verbose:
        cn = sorted(capped_nfan)
        print(f'fill_fanword: {len(seen)} centres touched, {n_touched} writes, '
              f'{n_cells} cells filled, {n_bad} CONTROL-FAILED, {n_skip} skipped; '
              f'{len(cn)} over cap' +
              (f' -- their `nfan` runs {cn[0]}..{cn[-1]} (median {cn[len(cn)//2]}), '
               f'so the UNFILLED rows are the LONG chains' if cn else ''), flush=True)
    return dict(seen=len(seen), touched=n_touched, cells=n_cells,
                control_failed=n_bad, skipped=n_skip,
                capped=len(capped_nfan), capped_nfan=sorted(capped_nfan))


def reindex_m(verbose=True):
    """Recompute `m`/`m_source` on every stored record in place.

    Cheap (no re-measurement) -- use it after the rung tables grow or after
    `_m_from_shape` changes, instead of a full refresh.
    """
    n = filled = 0
    for P, Q in known_centres():
        pth = _path(P, Q)
        rec = json.load(open(pth))
        if rec.get('guard'):
            continue
        lut, a1 = _rung_lut().get((P, Q), {}), rec['a1']
        for c in rec['cells']:
            m = lut.get(round(c['width'], 9))
            src = 'join:rung-stores' if m is not None else None
            if m is None:
                m = _m_from_shape(a1, c.get('nfan'), c.get('N'))
                src = 'shape:nfan+N' if m is not None else None
            c['m'], c['m_source'] = m, src
            n += 1
            filled += m is not None
        rec['rungs_present'] = sorted({c['m'] for c in rec['cells']
                                       if c['m'] is not None})
        json.dump(rec, open(pth, 'w'))
    if verbose:
        print(f'reindex_m: {filled}/{n} cells carry a rung index')
    return filled, n


if __name__ == '__main__':
    import sys
    cmd = sys.argv[1] if len(sys.argv) > 1 else ''
    if cmd == 'backfill':
        print('backfill:', backfill())
    if cmd == 'migrate':
        migrate()
    if cmd == 'reindex':
        reindex_m()
    if cmd == 'fanword':          # e.g. fanword max_nfan=1200000 hours=5
        kw = {}
        for a in sys.argv[2:]:
            k, v = a.split('=')
            if k == 'hours':
                kw['deadline'] = time.time() + float(v) * 3600
            elif k == 'ranks':
                kw['ranks'] = tuple(int(z) for z in v.split(','))
            elif k == 'stop_file':
                kw['stop_file'] = v
            else:
                kw[k] = int(float(v))
        print('fanword:', fill_fanword(**kw))
        raise SystemExit
    if cmd == 'rederive':
        rederive()
    if cmd == 'csv':
        print('rows:', export_csv())
    if cmd == 'recap':                # e.g. recap auto  |  recap budget=120 max_nodes=4e6
        kw, auto = {}, False
        for a in sys.argv[2:]:
            if a == 'auto':
                auto = True
                continue
            k, v = a.split('=')
            kw[k] = float(v) if k == 'budget' else int(float(v))
        print('recap:', recap(caps=kw or None, auto=auto))
    if cmd == 'peel':
        for P, Q in incomplete_centres():
            rec = json.load(open(_path(P, Q)))
            k = rec.get('kac') or {}
            print(f'  ({P:3d},{Q:4d}) a1={rec["a1"]:2d} cells={len(rec["cells"]):3d} '
                  f'sumw={rec.get("sumw") or 0:.6f} peel={k.get("frac") or 0:.4f}'
                  f'{"" if k.get("exact") else "*"} pending={k.get("n_pending")} '
                  f'next~{k.get("next_nodes")} [{k.get("basis")}]')
        print('  * = odd-P / a1=1: Kac denominator is an UPPER bound, so peel is a '
              'LOWER bound on progress')
    print('cell_store:', stats())
