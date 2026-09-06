"""s357_fanword.py -- PRE-REGISTERED: the FAN COUNT AS A RETURN TIME, and the step word.

THE QUESTION (the handoff's NEXT item 1, `(c‴)`(1a), since s345/s355; opened at the
user's direction).  s355 (3) + s356 (2) closed the RUNG block of a centre's `nfan`
multiset EXACTLY (A1's `min(4m+2, 4(a1-m))` over `m = 0..M`, twice, `0/293` above the
stack top).  What is left is the RESIDUAL TAIL -- the non-rung cells of the `2r < P`
cone -- where nothing predicts `nfan` and the counts explode (`(3,1)` max `nfan` runs
7, 11, 28, 132, 164, 384, 1188, 4538 at Q = 7..28 while the rung block stays [2,2]).

WHAT THIS PROBE DOES, AND WHY IT IS NOT ANOTHER LAW HUNT.  It reads the step rule of the
enumerator that PRODUCES `nfan` (`s340_glen_enum.enum_glen`) and asks what that rule makes
the count BE.  The rule is:  phi' = Om + 2Q with Om = phi + eps*jR*P, jR = j* - (j* mod 2),
and the cell CLOSES when Om = 0 (mod 2Q) ([WFLOOR-115]'s R-vertex head-on).  Since phi_0 = 0
every phi is P*(an even integer), so writing

        v_k := sum_{i<k} eps_i * (jR_i / 2)          (an integer walk)

closure reads  P*2*v_{k+1} = 0 (mod 2Q), i.e.  v = 0 (mod Q) once gcd is taken out.
⇒ `nfan` would then be the FIRST-RETURN TIME TO 0 (mod Q) of an explicit integer walk whose
steps are the signed half-R-indices -- turning "what sets #fans" into "what sets the step
word", which is the object [WFLOOR-101]'s `Sigma^(k)` and §19.4-b2's `[a1,1,a1-1,1,a1,0,...]`
already write down in CLOSED FORM for two solved families.

THE INSTRUMENT.  `trace_fanword` is `enum_glen`'s node loop specialised to a POINT s: no
interval splitting (a vertex landing exactly on s is measure zero), every comparison read
at s instead of at the endpoints.  It is NOT a new model -- it is the same arithmetic, and
the store is the control.  ⚠ The existing per-fan tracer `s313_domain_proof.full_trace` was
NOT used: it carries the out-of-domain closure defect (124/160, then 134/160 at s341, cause
still >= 2 mechanisms -- [WFLOOR-116], TOOLS.md trap (iv)), and residual cells are exactly
the off-domain ones.  `enum_glen`'s rule is the s341-amended one that reproduces the ring
51/51 odd-P + 40/40 even-P.

HYPOTHESES, FIXED BEFORE THE RUN.
  C   CONTROL.  On a stored cell, the pointwise trace from s = (lo+hi)/2 reproduces the
      store's `nfan`, `N` (= 2*letters+1), `jstar` and `nopp`.
                              PREDICT: holds.  This is the gate -- the interval enumerator
                              makes split decisions a point trace cannot, so a mismatch
                              means the specialisation is invalid and H1-H3 are void.
                              It is a 4-coordinate test on measured columns; it can fail.
  H1  `nfan` = the least k >= 1 with v_k = 0 (mod Q), and v_j != 0 (mod Q) for 0 < j < k.
                              PREDICT: holds at odd P.  At even P the modulus argument is
                              different (gcd(P,2Q) = 2, Q odd) and the reduction to
                              "v = 0 (mod Q)" uses 2 being invertible mod Q -- if the
                              algebra above is wrong this is where it shows.
  H2  PER-FAN CAP.  The step |g_i| is the fan's H-letter count, capped at a1 by
      [CHARGE-029] at `2r <= P`.  PREDICT: EXCEEDED, by exactly 1 -- that is [NCYL-169]
                              (iii), which found the cap broken "always by exactly 1" on
                              9/16 full-partition rows, all odd P with `2r <= P`.  Scored
                              here on a population ~50x larger, and separated rung vs
                              residual, which [NCYL-169] did not do.
  H3  WINDING.  `v_k` above IS [WFLOOR-112]'s PROVED winding index: iterating its own
      `Om_f = (f-1)pi + 2E_f*alpha` against phi_k = 2P*v_k + 2Q*k gives Om_{k+1} = k*pi +
      2*v_k*alpha, so E_f = v_{f-1} and G_B = max_k |v_k|.  Two questions on it:
      (a) [NCYL-169] (v) names, and explicitly leaves UN-RUN, "cap-breaking cells are
          exactly the seam cells with G_B = a1+1".  Scored here.
                              PREDICT: no prediction registered -- [NCYL-169] calls it
                              plausible and unverified, and this is the first run.
      (b) Is the closure an EXACT cancellation (v_end = 0) or a WRAP (v_end a nonzero
          multiple of Q)?     PREDICT: wraps dominate on the residual.  If instead |v| is
                              confined and every closure is v_end = 0, then `nfan` is NOT
                              a mod-Q return time in any useful sense and the explosion
                              is driven by the REAL part of the node state, not by Q.

PRIOR ART: grepped `rulings.md` (all tags), `TOOLS.md`, `cf_width_laws.md`,
`data/s349_priorart_cfans.md`, `session_state.md` for 'return time', 'first return',
'winding', 'nu-walk', 'step word', 'Sigma^(k)', 'per-fan', 'fan word', 'head-on
congruence', 'cocycle', 'residual tail', 'non-rung', 'second stack' ->
  - [WFLOOR-112] PROVES the fan phase is the accumulated winding `Om_f = (f-1)pi + 2E_f*alpha`
    and the ceiling |nu| <= a1 at both cones -- and `data/s349_priorart_cfans.md` D1 records,
    after a full read, that NO text in the repo caps the NUMBER of fans from it, in either
    direction.  H2 scores the ceiling on a new population; H3 is the un-asserted question.
  - [WFLOOR-101]/[WFLOOR-103] (A6) OWN the per-fan H-count word `Sigma^(k)` at `a2 = 1` and
    PROVE `#fans = 4(k+1)` there; `f_leg_covariance.md` §19.4-b2 owns the `a2 = 2` word.
    Both are step words for cells this probe does NOT score (it is the `2r < P` residual).
  - [WFLOOR-115] OWNS the R-vertex head-on predicate this derivation uses; [WFLOOR-102] owns
    the point-reflection hand-off.  Both PROVED; neither is re-derived here.
  - [NB-003]/C6 states the general obstruction to reading a count off a level cap ("a
    level-confined turning walk can revisit a level") -- which is precisely why H1 is posed
    as a RETURN TIME and not as a bound.
  - [OPS-054]/C3 bars hunting an arithmetic criterion over `(P,r)`; this probe reads what
    the enumerator's own step rule does, which is what [OPS-054] instructs instead.
  - [NCYL-154] already scored `closure/nfan` as the LOOSE control (0.419 vs 0.088) -- so no
    N-vs-nfan proportionality is re-run here.
  - [NCYL-185] killed the AFFINE residual-band renormalisation (0 hits); nothing found that
    traces the per-fan step word of a residual cell, at either cone.

Run:  PYTHONPATH=.:engine:probes .venv/bin/python3.13 probes/s357_fanword.py ctrl|walk|alpha|cocycle|all
"""
import collections
import json
import sys

import mpmath as mp

import cell_store as CS
from s302_model import base
from s257_warmup_tiling import stack_top
from s356_rung_exact import rung_width, scored_centres, WTOL

mp.mp.dps = 50
OUT = 'data/s357_fanword.json'
MAXFAN = 60000          # per-cell trace cap; a hit is a CAP, never a verdict


def trace_fanword(P, r, a1, s, max_fans=MAXFAN, scan_cap=None, keep_heights=False):
    """`enum_glen`'s node loop at a POINT.  Returns the per-fan step record.

    Same arithmetic as `s340_glen_enum.enum_glen`, with (lo, hi) collapsed to s:
    the interval splits (a vertex inside (lo,hi); the `cp` breakpoint) are measure
    zero at a point and drop out; every side test is read at s.

    Returns dict(closed, nfan, letters, N, jstar, nopp, steps, vs, epss, j0s, Ls).
    `steps[i]` is eps_i*(jR_i/2) and `vs[i]` the walk position BEFORE step i.

    `keep_heights=True` additionally returns `ds` (the DEVELOPED HEIGHT `d = c - s`
    BEFORE each fan) and `xs`.  ⚠ OFF by default and it must stay off for
    `cell_store.fill_fanword`: those are mpf arrays, ~100 B/fan, and the fill runs
    cells out to `nfan ~ 10^6`.  Nothing needs them to be stored anyway -- (s358)
    `d` is RECONSTRUCTIBLE from `vs`/`steps` alone, which is the point of §s358:
        c_{k+1} = c_k + 2*cot(al)*(-1)^k*sin(2*E_{k+1}*al),   E_{k+1} = vs[k]+steps[k]
    so keep them only when you want the float trace as its own control.
    """
    B = base(P, r, a1)
    Q, th, cot, csc = B['Q'], B['th'], B['cot'], B['csc']
    scan_cap = scan_cap or 64 * (a1 + 2)
    rho = lambda i: cot if i % 2 == 0 else csc
    two_q = 2 * Q

    c = mp.mpf(0)
    x = mp.mpf(0)
    phi, eps, j0, n, nfan, nopp = 0, 1, 2, 0, 0, 0
    v = 0
    steps, vs, epss, j0s, Ls = [], [], [], [], []
    ds, xs = [], []
    while nfan <= max_fans:
        if keep_heights:
            ds.append(c - s)
            xs.append(x)
        jstar = None
        i = j0
        while i - j0 <= scan_cap:
            y = c + rho(i) * mp.sin((phi + eps * i * P) * th)
            if (y <= s) if (c <= s) else (y >= s):
                jstar = i
                break
            i += 1
        if jstar is None:
            return dict(closed=False, why='noexit', nfan=nfan)
        n += jstar - j0 + 1
        ang = phi + eps * jstar * P
        jR = jstar - (jstar % 2)
        Om = phi + eps * jR * P
        vs.append(v)
        epss.append(eps)
        j0s.append(j0)
        Ls.append(jstar - j0 + 1)
        steps.append(eps * (jR // 2))
        v += eps * (jR // 2)
        nfan += 1
        if Om % two_q == 0:
            xR = x + cot * ((-1) ** ((Om // two_q) % 2))
            return dict(closed=True, nfan=nfan, letters=n, N=2 * n + 1,
                        jstar=jstar, nopp=nopp, fwd=float(1 - xR / cot),
                        steps=steps, vs=vs, epss=epss, j0s=j0s, Ls=Ls, Q=Q,
                        ds=ds, xs=xs)
        cp = 2 * (c + cot * mp.sin(Om * th)) - c
        xp = 2 * (x + cot * mp.cos(Om * th)) - x
        yj = c + rho(jstar) * mp.sin(ang * th)
        opp = (cp - s) * (yj - s) < 0
        if opp:                                        # Lemma B
            c, x, phi, eps, j0 = cp, xp, Om + two_q, -eps, 1 + (jstar % 2)
        else:                                          # [WFLOOR-113]'s other branch
            c, x, phi, eps, j0 = cp, xp, Om + two_q, eps, 2 - (jstar % 2)
            nopp += 1
        n += 1
    return dict(closed=False, why='cap', nfan=nfan)


# ----------------------------------------------------------------------------
def cone_cells(max_nfan=MAXFAN, max_centres=None):
    """Yield (P, Q, a1, r, M, cell, is_rung) over the cone centres with a clean strip."""
    seen = 0
    for P, Q, a1, r, rec, cs in scored_centres():          # 2r<P, complete, a1>=2
        M, ms = stack_top(P, r, a1)
        pool = list(cs)
        rungs = []
        for m in range(M + 1):
            wpred, _, _ = rung_width(P, Q, r, a1, m)
            for _ in range(2):
                k = min(range(len(pool)), key=lambda i: abs(pool[i]['width'] - wpred))
                if abs(pool[k]['width'] - wpred) < WTOL:
                    rungs.append(pool.pop(k))
        if len(rungs) != 2 * (M + 1):
            continue
        seen += 1
        if max_centres and seen > max_centres:
            return
        for cell in rungs:
            if cell['nfan'] and cell['nfan'] <= max_nfan:
                yield P, Q, a1, r, M, cell, True
        for cell in pool:
            if cell['nfan'] and cell['nfan'] <= max_nfan:
                yield P, Q, a1, r, M, cell, False


def _trace_cell(P, Q, a1, r, cell):
    s = (mp.mpf(repr(cell['lo'])) + mp.mpf(repr(cell['hi']))) / 2
    return trace_fanword(P, r, a1, s)


# ------------------------------------------------------------------- the sweep
def mode_sweep(max_centres=None, max_nfan=MAXFAN):
    """ONE pass: trace every traceable cone cell and store a per-cell summary row.

    Every later mode reads `data/s357_fanword.json`; the trace is the expensive part.
    """
    rows = []
    for P, Q, a1, r, M, cell, is_rung in cone_cells(max_nfan=max_nfan,
                                                    max_centres=max_centres):
        t = _trace_cell(P, Q, a1, r, cell)
        if not t['closed']:
            rows.append(dict(P=P, Q=Q, a1=a1, r=r, M=M, rung=is_rung,
                             nfan_store=int(cell['nfan']), closed=False,
                             why=t.get('why')))
            continue
        vs, steps = t['vs'], t['steps']
        vend = vs[-1] + steps[-1]
        Gb = max(abs(v) for v in vs)
        rows.append(dict(
            P=P, Q=Q, a1=a1, r=r, M=M, rung=is_rung, closed=True,
            nfan_store=int(cell['nfan']), nfan=t['nfan'],
            N_store=int(cell['N']), N=t['N'],
            jstar_store=int(cell['jstar']), jstar=t['jstar'],
            nopp_store=(int(cell['nopp']) if cell.get('nopp') is not None else None),
            nopp=t['nopp'], width=float(cell['width']),
            chi=(int(cell['chi']) if cell.get('chi') is not None else None),
            Gb=Gb, maxstep=max(abs(g) for g in steps), vend=vend,
            vend_mod=vend % Q, nzero=sum(1 for v in vs if v == 0),
            word=(steps if len(steps) <= 64 else None)))
        if len(rows) % 200 == 0:
            print(f"    ... {len(rows)} cells")
    print(f"  traced {len(rows)} cells, {sum(1 for z in rows if not z['closed'])}"
          f" not closed")
    return rows


# ------------------------------------------------------------------- C, H1, H2
def mode_ctrl(max_centres=40):
    """C + H1 + H2 in one sweep -- they read the same traces."""
    ok = collections.Counter()
    bad = []
    h1 = collections.Counter()
    amax = collections.Counter()
    pop = collections.Counter()
    for P, Q, a1, r, M, cell, is_rung in cone_cells(max_centres=max_centres):
        t = _trace_cell(P, Q, a1, r, cell)
        kind = 'rung' if is_rung else 'resid'
        pop[kind] += 1
        if not t['closed']:
            bad.append([P, Q, int(cell['nfan']), 'notclosed:' + t.get('why', '?')])
            ok['notclosed'] += 1
            continue
        hits = [(k, t[k] == int(cell[j]) if cell.get(j) is not None else None)
                for k, j in (('nfan', 'nfan'), ('N', 'N'),
                             ('jstar', 'jstar'), ('nopp', 'nopp'))]
        for k, good in hits:
            if good is None:
                ok[k + ':nodata'] += 1
            elif good:
                ok[k + ':ok'] += 1
            else:
                ok[k + ':BAD'] += 1
        if any(g is False for _, g in hits):
            bad.append([P, Q, int(cell['nfan']),
                        {k: (t[k], cell.get(j)) for k, j in
                         (('nfan', 'nfan'), ('N', 'N'), ('jstar', 'jstar'),
                          ('nopp', 'nopp'))}])
            continue
        # H1 -- first return of v to 0 mod Q
        vs, steps = t['vs'], t['steps']
        vend = vs[-1] + steps[-1]
        early = [j for j in range(1, len(vs)) if vs[j] % Q == 0]
        h1[(kind, vend % Q == 0, len(early) == 0)] += 1
        # H2 -- step alphabet against the proved ceiling a1
        amax[kind] = max(amax[kind], max(abs(g) for g in steps) - a1)
    print(f"  population: {dict(pop)}")
    print(f"  C  control vs store: {dict(ok)}")
    if bad:
        print(f"  ⚠ mismatches ({len(bad)}): {bad[:6]}")
    print(f"  H1 (kind, v_end==0 mod Q, no early return): {dict(h1)}")
    print(f"  H2 max(|step|) - a1 by kind (<=0 is the proved ceiling): {dict(amax)}")
    return dict(ctrl=dict(ok), bad=bad[:40],
                h1={str(k): v for k, v in h1.items()}, h2=dict(amax),
                pop=dict(pop))


# ------------------------------------------------------------------------- H3
def mode_cocycle(max_centres=25):
    """Is the step a function of the walk position (mod Q)?  Scored PER CENTRE."""
    per = collections.defaultdict(lambda: collections.defaultdict(dict))
    tally = collections.Counter()
    rows = []
    for P, Q, a1, r, M, cell, is_rung in cone_cells(max_centres=max_centres):
        t = _trace_cell(P, Q, a1, r, cell)
        if not t['closed']:
            continue
        for v, g, e, j0 in zip(t['vs'], t['steps'], t['epss'], t['j0s']):
            for name, key in (('v', (v % Q,)),
                              ('v,eps', (v % Q, e)),
                              ('v,eps,j0', (v % Q, e, j0))):
                d = per[(P, Q)][name]
                if key in d and d[key] != g:
                    tally[name + ':clash'] += 1
                else:
                    d.setdefault(key, g)
                    tally[name + ':ok'] += 1
    for (P, Q), d in sorted(per.items()):
        rows.append([P, Q] + [len(d[k]) for k in ('v', 'v,eps', 'v,eps,j0')])
    print(f"  H3 single-valuedness over all fans of all traced cells: {dict(tally)}")
    for name in ('v', 'v,eps', 'v,eps,j0'):
        c, k = tally[name + ':clash'], tally[name + ':ok']
        print(f"    {name:10s}  clashes {c:7d} / {c + k:7d} steps"
              f"  ({'FUNCTION' if c == 0 else 'NOT a function'})")
    print(f"  distinct keys per centre (P,Q,|v|,|v,eps|,|v,eps,j0|): {rows[:10]}")
    return dict(tally=dict(tally), rows=rows)


MODES = dict(sweep=mode_sweep, ctrl=mode_ctrl, cocycle=mode_cocycle)


def main(argv):
    want = argv[1:] or ['ctrl']
    if want == ['all']:
        want = list(MODES)
    out = {}
    for w in want:
        print(f"=== {w} ===")
        out[w] = MODES[w]()
    json.dump(out, open(OUT, 'w'), indent=1, default=str)
    print(f"wrote {OUT}")


if __name__ == '__main__':
    main(sys.argv)
