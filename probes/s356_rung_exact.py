"""s356 (c'''), item (1a) pass 1: is the s355 rung block a LOWER bound or the EXACT
rung content of a centre?  -- and is the residual band an affine copy of a child?

THE QUESTION, AND WHY IT IS NOT s355's.  s355 (3) scored A1's rung law by asking
whether the predicted multiset is a SUB-MULTISET of the observed `nfan`.  That test
is STRUCTURALLY ONE-SIDED: a sub-multiset test can only fail by ABSENCE, so it cannot
detect UNDER-prediction.  If `stack_top` returns an `M` smaller than the number of
rungs actually present, every such centre still PASSES.  s355 (6) saw the symptom --
residual cells carrying `4m+2`-shaped values -- and correctly refused to conclude
from it, leaving DO-NOT (vi): "do NOT assume `stack_top` under-counts ... that is
untested".  This tests it.

THE INSTRUMENT -- THREE COORDINATES, NONE OF THEM THE STORE'S `m`.
    `stack_top(P,r,a1) -> (M, m*)`   the geometric stack top (s257 (3); `M <= a1-1`
                                     PROVED [WFLOOR-045])
    `turn_data(a1,m)[2]`             A1's fan count `min(4m+2, 4(a1-m))`
    `Nlaw / Nwarm`                   the s241 length law `8(2m+1)a1+1-16m^2` and the
                                     warm cap `16(a1-m)(m+1)-3` (= `_closure_pred`'s
                                     two branches), scored against MEASURED `N`
    `min(V_m, X_m)` from `Cell`      the rung's PREDICTED CELL WIDTH (s257, and the
                                     `want` list of that probe's T5), scored against
                                     MEASURED `width` -- a CONTINUUM coordinate, so
                                     unlike small integers it cannot collide by luck
For every `m` in `M+1 .. a1-1` -- the rungs A5 PERMITS and `stack_top` EXCLUDES --
ask whether a J-pair with those coordinates is present.  CONTROL: the same predicate
at `m <= M`, where the rung IS geometrically present, must fire.

⚠ DISCLOSED DENOMINATOR INFLATION.  The `N` test accepts EITHER branch (`Nlaw` or
`Nwarm`).  That is deliberate -- `_m_from_shape`'s own docstring records `(11,71)`'s
`m=3` rung as DEGENERATE, measured `N = 189` against the law's `177`, rescued only by
`closure_pred(6,3) = 189` -- but it doubles the chance of a chance hit, so it inflates
the ALTERNATIVE hypothesis and makes the negative conservative.  Mode `degen` prints
that row.

⚠ CIRCULARITY ([OPS-071], and this is what cost s355 most of a session).  Nothing here
reads `cell_store`'s `m`: `provenance('m') == ('mixed', ('nfan','N'))` and
`independent('m','nfan') is False`.  The scored columns are `nfan`, `N` and `width`,
all three `measured` and all three `independent(_, 'nfan')` -- mode `prov` asserts it.

PRIOR ART: grepped `rulings.md` + the ledgers for 'J-pair', 'twin', 'in pairs',
'adjacent cell', 'n is odd', '2m* cells'.  What came back, and it is why this probe
claims NO structural novelty for the pairing:
  * `paper.md` Thm D -- `n(P/Q)` is odd exactly when `P` is odd, PROVED, and the lone
    odd branch is the orphan.  So "pairs plus one solo" is a THEOREM, not a finding.
  * `orphan_theorem.md` Thm A / Lemma B' -- the pairing IS the time-reversal involution
    `T`: it maps a branch ISOMETRICALLY onto a different branch of the SAME width, and
    `word(b-s) = reverse(word(s))`.  Equal width, equal word length, equal block count
    all follow; the fixed point is the orphan.
  * [SELECTOR-004] / [SELECTOR-005] -- NON-RUNG J-pairs already exist as a named object
    and are known to interleave with the rungs; [SEAM-017]/[SEAM-018] -- the rungs are
    a contiguous stack of J-pairs, PROVED at the threshold `t`.
  * `data/s349_priorart_cfans.md` A1/A2/A5, D3; C3/[OPS-054] (no criterion hunt).
Mode `twin` therefore reports the pairing as a CONFIRMATION AT SCALE plus one detail
the above do not state -- that the two cells of a pair are POSITIONALLY ADJACENT.

PRE-REGISTERED:
  prov    the provenance/circularity assertions above.
  twin    every centre's cells = adjacent equal-`(width,N,nfan)` pairs + <=1 solo;
          which columns agree inside a pair; the solo is the orphan where `chi` is
          populated.  Cross-engine: the same on the RING's own partition
          (`data/s345_bounce_decomp.json`), which shares no code with glen.
  ladder  the three-coordinate test.  PREDICTED: the control fires ~100%; the
          `m > M` rate is LOW but probably not 0, because `nfan` alone collides.
  degen   the `(11,71)` degenerate rung, printed because it is what forces the
          two-branch `N` test.
  renorm  is the rescaled residual band an AFFINE copy of a stored centre's full
          partition?  Scope: affine in the store's transversal coordinate, and only
          against centres the store holds.  PREDICTED: a few hits at best.

Run:
  PYTHONPATH=.:engine:archive/scripts_2026-07:archive/scripts_2026-08:probes \
      .venv/bin/python3.13 -u probes/s356_rung_exact.py all
"""
import collections
import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'engine'))
import cell_store as CS                                            # noqa: E402
from s254_avoidance_proof import Cell                              # noqa: E402
from s257_warmup_tiling import stack_top, turn_data                # noqa: E402

OUT = 'data/s356_rung_exact.json'
WTOL = 1e-9


def Nlaw(a1, m):
    return 8 * (2 * m + 1) * a1 + 1 - 16 * m * m


def Nwarm(a1, m):
    return 16 * (a1 - m) * (m + 1) - 3


def rung_width(P, Q, r, a1, m):
    z = Cell(P, Q, r, a1, m)
    return min(z.V, z.X), z.X, z.V


def scored_centres(cone_only=True):
    """Complete store centres with full `nfan`/`N` coverage, `a1 >= 2`."""
    for P, Q in CS.known_centres():
        rec = CS.get(P, Q)
        if not rec.get('complete'):
            continue
        cs = sorted(rec['cells'], key=lambda c: c['lo'])
        if any(c.get('nfan') is None or c.get('N') is None for c in cs):
            continue
        a1, r = rec['a1'], Q % P
        if a1 < 2 or r == 0:
            continue
        if cone_only and 2 * r >= P:
            continue
        yield P, Q, a1, r, rec, cs


# ------------------------------------------------------------------- prov
def mode_prov():
    out = {}
    for f in ('nfan', 'N', 'width', 'm'):
        out[f] = dict(provenance=list(CS.provenance(f)),
                      indep_of_nfan=CS.independent(f, 'nfan'))
        print(f"  {f:6s} {CS.provenance(f)!r:70s} independent(.,'nfan')="
              f"{CS.independent(f, 'nfan')}")
    assert CS.independent('N', 'nfan') and CS.independent('width', 'nfan')
    assert not CS.independent('m', 'nfan'), "the store's `m` must stay barred here"
    print("  OK -- the three scored columns are measured and nfan-independent;"
          " the store's `m` is NOT and is not read.")
    return out


# ------------------------------------------------------------------- twin
def _twin_scan(cells, wkey, fkey, nkey):
    """(adjacency holds?, #pairs, solo index or None) for one centre's cells."""
    key = lambda c: (fkey(c), nkey(c), round(wkey(c), 12))            # noqa: E731
    cnt = collections.Counter(key(c) for c in cells)
    odd = [k for k, v in cnt.items() if v % 2]
    if len(odd) > 1:
        return False, 0, None
    cand = [i for i, c in enumerate(cells) if odd and key(c) == odd[0]] or [None]
    for si in cand:
        rest = [c for i, c in enumerate(cells) if i != si]
        if len(rest) % 2:
            continue
        if all(abs(wkey(rest[2 * t]) - wkey(rest[2 * t + 1])) < 1e-11 and
               fkey(rest[2 * t]) == fkey(rest[2 * t + 1])
               for t in range(len(rest) // 2)):
            return True, len(rest) // 2, si
    return False, 0, None


def mode_twin():
    FIELDS = ['nfan', 'N', 'jstar', 'nopp', 'T_over_chi', 'chi', 'depth']
    tot = collections.Counter()
    agree, seen = collections.Counter(), collections.Counter()
    worst = 0.0
    for P, Q, a1, r, rec, cs in scored_centres(cone_only=False):
        ok, npair, si = _twin_scan(cs, lambda c: c['width'],
                                   lambda c: int(c['nfan']), lambda c: c.get('N'))
        tot['centres'] += 1
        tot['ADJACENCY_HOLDS' if ok else 'ADJACENCY_FAILS'] += 1
        if not ok:
            continue
        tot['pairs'] += npair
        tot['solo' if si is not None else 'no_solo'] += 1
        if si is not None and cs[si].get('chi') is not None:
            tot['solo_IS_orphan' if cs[si].get('is_orphan')
                else 'solo_NOT_orphan'] += 1
        rest = [c for i, c in enumerate(cs) if i != si]
        for t in range(len(rest) // 2):
            a, b = rest[2 * t], rest[2 * t + 1]
            worst = max(worst, abs(a['width'] - b['width']))
            for f in FIELDS:
                if a.get(f) is not None and b.get(f) is not None:
                    seen[f] += 1
                    agree[f] += (a[f] == b[f])
    print(f"  glen store: {dict(tot)}")
    print(f"  worst ABSOLUTE width gap inside a pair: {worst:.3e}")
    print("  in-pair column agreement: " +
          ", ".join(f"{f} {agree[f]}/{seen[f]}" for f in FIELDS if seen[f]))
    # cross-engine: the RING's own partition, which shares no code with glen
    ring = dict(rows=0, ok=0)
    try:
        for R in json.load(open('data/s345_bounce_decomp.json')):
            cs = sorted(R['cells'], key=lambda c: c['lo'])
            ok, _, _ = _twin_scan(cs, lambda c: c['hi'] - c['lo'],
                                  lambda c: int(c['nb1']), lambda c: int(c['nb1']))
            ring['rows'] += 1
            ring['ok'] += ok
    except FileNotFoundError:
        pass
    print(f"  INDEPENDENT ENGINE (ring, data/s345_bounce_decomp.json):"
          f" adjacency {ring['ok']}/{ring['rows']} rows")
    print("  ⚠ NOT A NEW STRUCTURE -- this is Thm D + orphan_theorem Lemma B'"
          " (the time-reversal involution) confirmed at scale; see PRIOR ART.")
    return dict(tot=dict(tot), worst_width_gap=worst, ring=ring,
                agree={f: [agree[f], seen[f]] for f in FIELDS})


# ----------------------------------------------------------------- ladder
def mode_ladder():
    hit = collections.Counter()
    hits_out, per = [], []
    for P, Q, a1, r, rec, cs in scored_centres():
        M, ms = stack_top(P, r, a1)
        obs = collections.Counter((int(c['nfan']), int(c['N'])) for c in cs)
        obsf = collections.Counter(int(c['nfan']) for c in cs)
        wobs = [c['width'] for c in cs]
        row = dict(P=P, Q=Q, a1=a1, r=r, a2=P // r, M=M, n=len(cs), rungs=[])
        for m in range(0, a1):
            f = turn_data(a1, m)[2]
            wpred, X, V = rung_width(P, Q, r, a1, m)
            nb = [N for N in (Nlaw(a1, m), Nwarm(a1, m)) if obs.get((f, N), 0) >= 2]
            wmatch = sum(1 for w in wobs if abs(w - wpred) < WTOL)
            side = 'in' if m <= M else 'out'
            tag = ('BOTH' if nb else
                   'nfan_only' if obsf.get(f, 0) >= 2 else 'neither')
            hit[f'{side}:{tag}'] += 1
            hit[f'{side}:width>=2'] += (wmatch >= 2)
            if nb and side == 'out':
                hit['out:BOTH_and_width' if wmatch >= 2
                    else 'out:BOTH_but_width_FAILS'] += 1
                hits_out.append(dict(P=P, Q=Q, a1=a1, r=r, a2=P // r, M=M, m=m,
                                     nfan=f, N=nb[0], X=X, V=V, wpred=wpred,
                                     wmatch=wmatch,
                                     wmeas=sorted({round(c['width'], 10) for c in cs
                                                   if int(c['nfan']) == f})))
            row['rungs'].append(dict(m=m, side=side, nfan=f, tag=tag,
                                     wpred=wpred, wmatch=wmatch))
        per.append(row)
    tin = sum(v for k, v in hit.items() if k.startswith('in:') and 'width' not in k)
    tout = sum(v for k, v in hit.items() if k.startswith('out:')
               and 'width' not in k and 'BOTH_' not in k)
    print(f"  centres in the 2r<P cone with a full read: {len(per)}")
    print(f"  CONTROL  m <= M : (nfan,N) BOTH {hit['in:BOTH']}/{tin}"
          f" | predicted WIDTH present twice {hit['in:width>=2']}/{tin}")
    print(f"  TEST     m >  M : (nfan,N) BOTH {hit['out:BOTH']}/{tout}"
          f" | nfan only {hit['out:nfan_only']}/{tout}"
          f" | neither {hit['out:neither']}/{tout}")
    print(f"           of those {hit['out:BOTH']} joint hits, width ALSO matches"
          f" {hit['out:BOTH_and_width']}, width FAILS"
          f" {hit['out:BOTH_but_width_FAILS']}")
    print(f"  ==> rungs present above the stack top, all three coordinates:"
          f" {hit['out:BOTH_and_width']}/{tout}")
    print("  the joint (nfan,N) hits, with the width coordinate that settles them:")
    for h in hits_out:
        print(f"    ({h['P']:3d},{h['Q']:3d}) a1={h['a1']:2d} r={h['r']} M={h['M']}"
              f" m={h['m']} nfan={h['nfan']:3d} N={h['N']:4d}"
              f" X_m={h['X']:+.3e} wpred={h['wpred']:+.6f} wmeas={h['wmeas'][:3]}")
    return dict(hist=dict(hit), out_hits=hits_out, per=per, tin=tin, tout=tout)


# ------------------------------------------------------------------ degen
def mode_degen():
    P, Q, a1, r, m = 11, 71, 6, 5, 3
    rec = CS.get(P, Q)
    f = turn_data(a1, m)[2]
    got = [(int(c['nfan']), int(c['N'])) for c in rec['cells']
           if int(c['nfan']) == f]
    print(f"  ({P},{Q}) a1={a1} m={m}: A1 nfan={f}, Nlaw={Nlaw(a1, m)},"
          f" Nwarm={Nwarm(a1, m)}; measured (nfan,N) = {sorted(set(got))}")
    print("  ==> the law's N misses and the WARM branch supplies it. This row is why"
          " the ladder test accepts either branch, which inflates the ALTERNATIVE.")
    return dict(P=P, Q=Q, m=m, nfan=f, Nlaw=Nlaw(a1, m), Nwarm=Nwarm(a1, m),
                measured=sorted(set(got)))


# ----------------------------------------------------------------- renorm
def mode_renorm():
    full, resid = {}, {}
    for P, Q, a1, r, rec, cs in scored_centres(cone_only=False):
        full[(P, Q)] = sorted((c['width'], int(c['nfan'])) for c in cs)
    for P, Q, a1, r, rec, cs in scored_centres():
        M, ms = stack_top(P, r, a1)
        pool = list(cs)
        stripped = 0
        for m in range(M + 1):
            wpred, _, _ = rung_width(P, Q, r, a1, m)
            for _ in range(2):
                k = min(range(len(pool)),
                        key=lambda i: abs(pool[i]['width'] - wpred))
                if abs(pool[k]['width'] - wpred) < WTOL:
                    pool.pop(k)
                    stripped += 1
        if stripped != 2 * (M + 1):
            continue
        tot = sum(c['width'] for c in pool)
        if tot <= 0:
            continue
        resid[(P, Q)] = (sorted((c['width'] / tot, int(c['nfan'])) for c in pool),
                         tot, len(pool))
    matches = []
    for k, (rs, tot, n) in resid.items():
        for k2, fl in full.items():
            if len(fl) != n:
                continue
            if all(abs(a[0] - b[0]) < 1e-6 for a, b in zip(rs, fl)):
                matches.append([list(k), list(k2), n,
                                'W+F' if all(a[1] == b[1] for a, b in zip(rs, fl))
                                else 'W only'])
    print(f"  clean rung-strips: {len(resid)}/{len(list(scored_centres()))}"
          f" (the strip finds exactly 2(M+1) cells at the predicted widths)")
    print(f"  candidates compared against: {len(full)} stored full partitions")
    print(f"  rescaled residual == a stored full partition: {len(matches)} hits"
          f"  {matches[:8]}")
    print("  ⚠ SCOPE: AFFINE in the store's transversal coordinate, and only against"
          " centres the store holds. It does NOT bar a non-affine renormalisation.")
    return dict(n_resid=len(resid), n_cand=len(full), matches=matches,
                sizes=dict(collections.Counter(v[2] for v in resid.values())))


MODES = dict(prov=mode_prov, twin=mode_twin, ladder=mode_ladder,
             degen=mode_degen, renorm=mode_renorm)

if __name__ == '__main__':
    want = sys.argv[1:] or ['all']
    if want == ['all']:
        want = list(MODES)
    res = {}
    for w in want:
        print(f"\n=== {w} " + "=" * (66 - len(w)))
        res[w] = MODES[w]()
    json.dump(res, open(OUT, 'w'), indent=1, default=float)
    print(f"\nwrote {OUT}")
