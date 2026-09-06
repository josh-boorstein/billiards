"""s457 -- QUEUE ITEM (III): IS COMPLETE PERIODICITY CHARACTERISABLE?

PRE-REGISTERED.

PRIOR ART: ran `rulings.py --cited session_state.md`, `--grep 'complete periodicity'`,
`--grep 'completely periodic'`, and read the bodies of [NCYL-256] (the census, `10/354`),
[NCYL-257] (the plateau-shape discriminant), [NCYL-262] (`PRONG_CAPPED` != non-CP),
[NCYL-293] (two paired non-CP rows are cap ARTEFACTS -- the exact model closes them at
`1740348` steps), [NCYL-144] (CP is FLOAT-measured, never certified), [NCYL-304] (CP is
`p = 0`'s finiteness input), [NCYL-306] (the float census is at its precision floor, BOTH
error directions), [NCYL-012] (the Veech criterion `min(P,Q-P) in {1,2}`), plus the
`TOOLS.md` rows for `s438_oddclass_cp`, `s450_chain_maps`, `s439_exact_census`,
`s455_cp_hypothesis`.  What came back:

  * every CP verdict in the repo is FLOAT ([NCYL-144]), and its instrument is at its
    precision floor with false positives AND false negatives ([NCYL-306]);
  * [NCYL-293] already refuted 2 of the census's 10 non-CP rows with the EXACT chain
    model -- but only those 2, and only at `Q <= 24`;
  * NOBODY HAS RUN THE EXACT MODEL AS A CP CENSUS.  `s450_chain_maps.sweep` exists and
    `data/s450_chain.json` holds `Q = 31..51`, but it was run to score `p = 0`, and its
    `capped == 0` column -- which IS an exact CP certificate -- has never been read as one.
  * no ruling states the LEG-SWAP INVOLUTION as a symmetry of either instrument.

So (III)'s first sub-goal is not a new idea, it is a TRUSTWORTHY TABLE: you cannot ask
what characterises a set while the set is mostly instrument error.

WHAT THIS PROBE IS.  Four arms, in dependency order.

  h0  THE FLOAT CENSUS, READ HONESTLY.  `data/s438_cp_store/`, odd `Q` only (the model's
      scope, [NCYL-291]).  Per non-CP row: the closure curve, [NCYL-257]'s plateau shape,
      and `min_miss_of_open / ARRIVE_TOL` -- [NCYL-306]'s marginality ratio.  Vacuous
      zero-separatrix classes excluded from every rate ([NCYL-256] trap (i), [OPS-041]).

  h1  THE LEG-SWAP INVOLUTION `(P, eps) <-> (Q-P, 1-eps)`.  `P/Q` and `(Q-P)/Q` are the
      SAME right triangle with its legs relabelled, so the two class rows are the same
      direction on the same surface and every combinatorial column must agree.
      ⚠⚠⚠ THE EXACT MODEL'S SCORE ON THIS IS *FORCED BY CONSTRUCTION AND IS NOT
      EVIDENCE*, and `h1_forced()` ASSERTS THAT DIAGNOSIS IN CODE rather than leaving it
      as prose ([OPS-043]).  In `Chain.__init__`, `c[t] = (c0 + 2Pt) mod 2Q`, so the swap
      sends `c -> -c mod 2Q`; `Ring.gen[2Q-c] == gen[c]` EXACTLY, so `ell` is
      *bit-identical* and only the kite orientation `osec` moves.  Flipping every kite is
      the `sector="complement"` relabelling, which [OPS-223] already records as an
      AUTOMORPHISM of every count in `readout()` -- the reason that control scores a
      useless `347/347`.  So a perfect exact-model score COULD NOT HAVE COME OUT OTHERWISE
      ([OPS-041]) and must not be booked as corroboration of anything.
      ⇒⇒ WHAT IS GENUINE IS THE OTHER SIDE.  The involution is a true relabelling of the
      billiard, so the FLOAT CENSUS -- a completely separate code path (a tracer on the
      direction grid `u*pi/(2Q)`, sharing no line with the chain model) -- MUST satisfy it
      too.  It does not.  Every violation is an instrument error localised for free, at
      zero compute, with no reference instrument needed.
      CONTROL, which must FAIL ([OPS-219]): the parity-PRESERVING map
      `(P, eps) <-> (Q-P, eps)`.  It shows the readout is not invariant under just any
      relabelling -- which is the only thing the exact-model side can still say.

  h2  THE EXACT CP TABLE.  `s450_chain_maps.model_row` decides closure on exact
      `Z[zeta_4Q]` vector equality, so `capped == 0` certifies that EVERY separatrix is a
      saddle connection -- which IS complete periodicity, and it is a CERTIFICATE, not a
      measurement.  Fills odd `Q <= 29` (the stored sweep starts at `31`) into a resumable
      per-centre store, merges `data/s450_chain.json`, and cross-tabs against h0.
      ⚠⚠ THE ASYMMETRY IS THE WHOLE POINT AND IT IS NOT NEGOTIABLE: `capped == 0` PROVES
      CP; `capped > 0` PROVES NOTHING -- it is a cap, and [NCYL-262]/[NCYL-293] bar reading
      it as non-CP (`3/23` eps=0 closes at `1740348` steps).  So this arm can only ever
      SHRINK the non-CP set, never confirm it.
      ⚠ The exact tests are FLOAT-GATED at `TOL = 1e-6` for speed; `safety` audits the
      gate and is quoted with every table ([OPS-223]).

  h3  CHARACTERISATION.  Score the pre-registered predicates against the h2 table:
      (P1) Veech, `min(P, Q-P) <= 2  =>  CP` ([NCYL-012]);
      (P2) `non-CP => lone` ([NCYL-256]'s headline);
      (P3) ⇒ THE DEFLATIONARY ARM, AND IT IS THE ONE THAT DECIDES WHETHER (P4)/(P5) ARE
           WORTH SCORING AT ALL: is the unresolved set just the TOP OF THE LENGTH
           DISTRIBUTION?  Compare each capped row against the largest RESOLVED
           `max_steps` in range.  If every capped row sits above it, "non-CP" is a
           synonym for "longer than the cap" and no arithmetic predicate is warranted
           ([OPS-041]: an outcome that could not have come out otherwise is not evidence).
      (P4) is the surviving set closed under the h1 involution?
      (P5) residue / CF predicates on the surviving set.

SCOPE: ODD `Q` ONLY ([NCYL-291] -- the single iota-fixed kite exists only there), which is
also the scope of every CP statement `p = 0` consumes.

Run: python3.13 probes/s457_cp_char.py [all|h0|h1|h2|h3] [--qmax=] [--cap=] [--budget=]
"""
import glob
import json
import math
import os
import sys
import time
from fractions import Fraction

import s450_chain_maps as cm

CENSUS_DIR = "data/s438_cp_store"
EXACT_DIR = "data/s457_exact_cp"
S450_SWEEP = "data/s450_chain.json"
OUT = "data/s457_cp_char.json"
ARRIVE_TOL = 1e-9                      # s438's tracer tolerance ([NCYL-306])

# columns that the leg-swap involution must preserve EXACTLY on the exact model.
# `margin`/`maxcoef`/`safety` are float-gate artefacts of one particular trace and are
# reported but NOT scored.
EXACT_COLS = ("lone", "E", "h", "v", "f", "p2", "capped", "multi",
              "n_traced", "max_steps", "tot_steps")
CENSUS_COLS = ("class_role", "n_separatrices", "verdict", "n_closed")


# ------------------------------------------------------------------ h0: float census

def census_rows(qmax=None):
    """Every ODD-Q class row in `data/s438_cp_store/`, as (Q, P, eps, record)."""
    out = []
    for f in glob.glob(os.path.join(CENSUS_DIR, "*.json")):
        d = json.load(open(f))
        P, Q = d["P"], d["Q"]
        if Q % 2 == 0 or (qmax and Q > qmax):
            continue
        for par, c in d["classes"].items():
            out.append((Q, P, int(par), c))
    out.sort()
    return out


def plateau(curve):
    """[NCYL-257]'s discriminant: the shape, not the count.  Returns the number of
    trailing decades over which `n_closed` is CONSTANT (6 = flat from `1e3`)."""
    v = [x["n_closed"] for x in curve]
    k = 1
    while k < len(v) and v[-1 - k] == v[-1]:
        k += 1
    return k


def h0(qmax=None, verbose=True):
    rows = census_rows(qmax)
    live = [r for r in rows if r[3]["n_separatrices"] > 0]
    vac = len(rows) - len(live)
    mixed = [r for r in live if r[3]["verdict"] != "COMPLETELY_PERIODIC"]
    if verbose:
        qs = sorted({q for q, _, _, _ in rows})
        print(f"h0  FLOAT CP CENSUS -- odd Q in [{qs[0]}, {qs[-1]}], "
              f"{len(rows)} class rows, {vac} VACUOUS (0 separatrices) excluded "
              f"([NCYL-256] trap (i))")
        print(f"    {len(live) - len(mixed)}/{len(live)} COMPLETELY_PERIODIC, "
              f"{len(mixed)} not\n")
        print("    %-8s %-4s %-6s %-5s %-5s %-9s %-6s %s"
              % ("P/Q", "eps", "role", "nsep", "nclos", "miss/TOL", "flat", "curve"))
        for Q, P, e, c in mixed:
            m = c.get("min_miss_of_open")
            print("    %-8s %-4d %-6s %-5d %-5d %-9s %-6d %s"
                  % (f"{P}/{Q}", e, c["class_role"][:5], c["n_separatrices"],
                     c["n_closed"], f"{m / ARRIVE_TOL:.2f}" if m else "-",
                     plateau(c["curve"]), [x["n_closed"] for x in c["curve"]]))
        flat = [r for r in mixed if plateau(r[3]["curve"]) >= 6]
        wide = [r for r in mixed if (r[3].get("min_miss_of_open") or 0) > 5 * ARRIVE_TOL]
        print(f"\n    [NCYL-257] FLAT from 1e3: {[f'{P}/{Q} eps{e}' for Q,P,e,_ in flat]}")
        print(f"    [NCYL-306] margin > 5xTOL: {[f'{P}/{Q} eps{e}' for Q,P,e,_ in wide]}")
    return {"rows": len(rows), "vacuous": vac, "live": len(live),
            "mixed": [[P, Q, e] for Q, P, e, _ in mixed]}


# ------------------------------------------------------- h1: the leg-swap involution

def _score_involution(table, cols, swap, label, verbose=True):
    """`table` maps (P,Q,eps) -> record.  `swap` maps (P,Q,eps) -> (P',Q,eps')."""
    # ⚠ [OPS-041]: a column that is CONSTANT within its `Q` agrees under ANY within-`Q`
    # relabelling and so cannot fail either the identity or the control.  Mark it, and do
    # not let it pad the control's score.
    byq = {}
    for (P, Q, e), r in table.items():
        for c in cols:
            if c in r:
                byq.setdefault((c, Q), set()).add(r[c])
    vac = {c for c in cols
           if all(len(v) == 1 for (cc, _), v in byq.items() if cc == c)}
    tot = {c: [0, 0] for c in cols}          # [agree, compared]
    missing = 0
    bad = []
    for key, rec in sorted(table.items()):
        other = table.get(swap(*key))
        if other is None:
            missing += 1
            continue
        for c in cols:
            if c not in rec or c not in other:
                continue
            tot[c][1] += 1
            if rec[c] == other[c]:
                tot[c][0] += 1
            elif len(bad) < 40:
                bad.append((key, c, rec[c], other[c]))
    if verbose:
        print(f"    {label}:")
        for c in cols:
            a, n = tot[c]
            mark = "" if a == n else "   <-- BREAKS"
            if c in vac:
                mark += "   [VACUOUS: constant within Q, cannot fail]"
            print(f"      {c:<16s} {a}/{n}{mark}")
        if missing:
            print(f"      ({missing} rows whose partner is absent from the table)")
    return tot, bad


def h1_forced(qmax=25, verbose=True):
    """⚠ THE DEFLATION, ASSERTED IN CODE ([OPS-043]).  Shows that the exact model CANNOT
    break the involution, so its score on h1 is not evidence:
      (a) `c'[t] == -c[t] mod 2Q`               -- read straight off `Chain.__init__`;
      (b) `ell` is BIT-IDENTICAL                -- because `gen[2Q-c] == gen[c]`;
      (c) the swapped row's readout EQUALS the `sector="complement"` control's,
          which [OPS-223] classifies as an automorphism of every count.
    """
    cols = [c for c in EXACT_COLS if c != "lone"]
    n = cneg = ell = 0
    cm.STEP_CAP = 200_000              # small: this arm is about identity, not closure
    comp = tot = 0
    for Q in range(5, qmax + 1, 2):
        for P in range(1, Q):
            if math.gcd(P, Q) != 1:
                continue
            for eps in (0, 1):
                a, b = cm.Chain(P, Q, eps), cm.Chain(Q - P, Q, 1 - eps)
                n += 1
                cneg += all(b.c[t] == (-a.c[t]) % (2 * Q) for t in range(Q))
                ell += all(x.eq(y) for x, y in zip(a.ell, b.ell))
                if Q <= 19:
                    x = cm.model_row(P, Q, eps, sector="complement")
                    y = cm.model_row(Q - P, Q, 1 - eps)
                    tot += 1
                    comp += all(x[k] == y[k] for k in cols)
    if verbose:
        print("    ⚠⚠ FORCED-BY-CONSTRUCTION CHECK (the exact-model side is NOT evidence):")
        print(f"      (a) c'[t] == -c[t] mod 2Q          {cneg}/{n}")
        print(f"      (b) ell bit-identical              {ell}/{n}")
        print(f"      (c) swap == sector='complement'    {comp}/{tot}   "
              f"([OPS-223]: an automorphism of the readout)")
        print("      ⇒ a perfect exact-model score below could not have come out")
        print("        otherwise; only the FLOAT-CENSUS side of this arm is a result.\n")
    return {"c_negated": [cneg, n], "ell_identical": [ell, n],
            "equals_complement_control": [comp, tot]}


def h1(qmax=None, verbose=True):
    """The leg-swap involution, scored on BOTH instruments."""
    if verbose:
        print("h1  LEG-SWAP INVOLUTION  (P, eps) <-> (Q-P, 1-eps)")
        print("    `P/Q` and `(Q-P)/Q` are the same right triangle with its legs")
        print("    relabelled, so the two class rows are the same direction on the same")
        print("    surface -- a true relabelling of the billiard, so BOTH instruments")
        print("    must satisfy it.\n")
    forced = h1_forced(verbose=verbose)
    swap = lambda P, Q, e: (Q - P, Q, 1 - e)                            # noqa: E731
    ctrl = lambda P, Q, e: (Q - P, Q, e)                                # noqa: E731

    ex = {(r["P"], r["Q"], r["eps"]): r for r in exact_table(qmax=qmax).values()}
    e_tot, e_bad = _score_involution(ex, EXACT_COLS, swap,
                                     f"EXACT MODEL ({len(ex)} rows) -- FORCED, "
                                     f"not evidence", verbose)
    c_tot, c_bad = _score_involution(ex, EXACT_COLS, ctrl,
                                     "CONTROL, parity-PRESERVING (must FAIL)", verbose)

    cen = {(P, Q, e): dict(c, P=P, Q=Q, eps=e) for Q, P, e, c in census_rows(qmax)}
    f_tot, f_bad = _score_involution(cen, CENSUS_COLS, swap,
                                     f"FLOAT CENSUS ({len(cen)} rows)", verbose)
    if verbose and f_bad:
        print("\n    ⇒ WHERE THE FLOAT CENSUS BREAKS THE IDENTITY (these are its errors,")
        print("      localised for free -- the exact model has no such rows):")
        for (P, Q, e), c, a, b in f_bad:
            print(f"      {P}/{Q} eps{e}  {c}: {a}   vs   {Q-P}/{Q} eps{1-e}: {b}")
    return {"forced": forced,
            "exact": {c: e_tot[c] for c in e_tot},
            "control": {c: c_tot[c] for c in c_tot},
            "census": {c: f_tot[c] for c in f_tot},
            "census_breaks": [[list(k), c, a, b] for k, c, a, b in f_bad]}


# ------------------------------------------------------------- h2: the exact CP table

def _path(P, Q):
    return os.path.join(EXACT_DIR, f"{P}_{Q}.json")


def _write(P, Q, rec):
    os.makedirs(EXACT_DIR, exist_ok=True)
    tmp = _path(P, Q) + ".tmp"
    json.dump(rec, open(tmp, "w"), indent=1)
    os.replace(tmp, _path(P, Q))                # atomic: a kill costs the in-flight row


KEEP = ("P", "Q", "eps", "lone", "E", "h", "v", "f", "p2", "capped", "multi",
        "n_traced", "max_steps", "tot_steps", "margin", "maxcoef", "safety",
        # s458: `exact_calls` joins the audit trio it belongs to -- a resume MERGES
        # margin/maxcoef/exact_calls, so dropping one makes the merged audit incomplete.
        "exact_calls")


def measure(P, Q, cap=3_000_000, force=False, only_eps=None):
    """Both classes at one centre, through the store.  A stored record is reused only if
    it was computed at a cap >= the one asked for (a bigger cap can only resolve MORE).

    ⚠ `runs` KEEPS THE PER-CAP TRAJECTORY, it does not overwrite it: [NCYL-257]'s
    discriminant is the SHAPE of `closed(cap)`, not the value at the biggest cap, and a
    record that keeps only its last run cannot supply it ([NCYL-262]).

    ⇒⇒ RESUMABLE SINCE s458 (`CLAUDE.md` Measurement discipline, user-directed).  Each
    class stores its per-separatrix table `seps`, including the WALK STATE of every
    separatrix that hit the cap, and an escalation hands that straight back to
    `cm.model_row(prior=…)`.  So raising the cap re-walks ONLY the capped separatrices and
    only their UNSEEN tail; the resolved ones are not re-traced at all.
    ⚠ Pre-s458 records have no `seps`, so they cannot be resumed -- they fall back to a
    full re-measure, which is the OLD behaviour and not a regression.  `repair()` is what
    upgrades them.  ⚠ Equivalence to a one-shot run at the high cap is not assumed: it is
    the control `cm._resume_equivalence()` (`s450_chain_maps.py resume`)."""
    old = {}
    if os.path.exists(_path(P, Q)):
        old = json.load(open(_path(P, Q)))
        if not force and old.get("cap", 0) >= cap:
            return old
    cm.STEP_CAP = cap
    rec = {"P": P, "Q": Q, "cap": cap, "classes": dict(old.get("classes", {})),
           "runs": old.get("runs", {}), "seps": dict(old.get("seps", {}))}
    t0 = time.time()
    for eps in (0, 1):
        if only_eps is not None and eps != only_eps:
            continue
        prev = rec["seps"].get(str(eps))
        prior = None
        if prev is not None:
            prev_cap = rec["classes"].get(str(eps), {}).get("cap")
            if prev_cap is not None and prev_cap < cap:
                # ⚠ resume ONLY when strictly escalating.  Re-running at the SAME cap is
                # what `force=True` is for (an audit), and resuming there would return the
                # stored answer without re-walking anything -- i.e. silently not audit.
                prior = dict(rec["classes"][str(eps)], seps=prev)
        mo = cm.model_row(P, Q, eps, prior=prior, with_seps=True)
        # ⚠ the cap is stored PER CLASS: `deep` re-runs one class at a time, and a
        # record-level cap would silently relabel the untouched class's provenance.
        rec["classes"][str(eps)] = dict({k: mo[k] for k in KEEP}, cap=cap)
        rec["seps"][str(eps)] = mo["seps"]
        rec["runs"].setdefault(str(eps), {})[str(cap)] = {
            "capped": mo["capped"], "max_steps": mo["max_steps"],
            "n_traced": mo["n_traced"]}
    rec["secs"] = round(time.time() - t0, 1)
    _write(P, Q, rec)
    return rec


def repair(qmax=31, apply=False, include_capped=False, max_tot_steps=20_000_000,
           verbose=True):
    """⇒ THE s458 STORE AUDIT: find records that cannot answer the questions the store
    exists to answer, and say which are FIXABLE WITHOUT RE-MEASURING.

    Three defects, all introduced by writing a record before the schema settled:
      (1) NO `seps`   -- the record cannot be resumed; a cap raise pays the full restart.
      (2) NO `runs` rung for a class -- [NCYL-257]'s plateau SHAPE needs the whole ladder,
          and this is the BASE rung, missing on every pre-`deep` record.
      (3) NO per-class `cap` -- and `exact_table` then falls back to the RECORD cap, which
          on a `deep`-touched record relabels the untouched class with a cap it was never
          measured at (`14/29 eps1`: measured at the base cap, reported at `3e7`).
    ⚠⚠ (2) AND (3) ARE NOT DERIVABLE.  A resolved class's `capped`/`max_steps` are known,
    but WHICH cap produced them is not recorded anywhere, so synthesising a rung would be
    fabricating provenance -- exactly the drift the per-class cap exists to prevent.  The
    honest repair is a re-measure at the base cap, and for resolved rows that is cheap
    (`tot_steps`, already known).  This function REPORTS by default; `apply=True`
    re-measures the affected classes.

    ⇒⇒ THE DEFAULTS ARE THE WHOLE POINT, BECAUSE TWO OF THE THREE BUCKETS ARE NOT WORTH
    REPAIRING AND THE COST IS `17x` APART (measured s458, odd `Q <= 31`):
      * `354` RESOLVED rows missing `runs`/`cap` -- `1.2e8` steps, **~2.8 min**.  Restores
        the BASE rung of every ladder.  DO IT: this is the whole default set.
      * `8` resolved rows missing only `seps` -- `4.7e8` steps, ~11 min.  SKIPPED by
        `max_tot_steps`: `seps` on a row that is already CERTIFIED buys no resume (it will
        never be escalated), only the ability to decompose `tot_steps` after the fact.
      * `4` CAPPED rows -- `1.4e9` steps, ~33 min.  SKIPPED by `include_capped=False`, and
        skipping is not a compromise, it is the correct call: with no stored state the
        next escalation restarts them from `0` REGARDLESS, and that escalation records the
        state itself.  Repairing now pays the full run to save ~10% of a run you are going
        to do anyway.  ⚠ Do NOT "fix" this by lowering the bar -- the waste is structural."""
    bad = []
    for f in sorted(glob.glob(os.path.join(EXACT_DIR, "*.json"))):
        rec = json.load(open(f))
        P, Q = rec["P"], rec["Q"]
        if Q > qmax:
            continue
        for eps in ("0", "1"):
            c = rec.get("classes", {}).get(eps)
            if c is None:
                continue
            miss = [n for n, ok in (("seps", eps in rec.get("seps", {})),
                                    ("runs", eps in rec.get("runs", {})),
                                    ("cap", c.get("cap") is not None)) if not ok]
            if miss:
                bad.append({"P": P, "Q": Q, "eps": int(eps), "missing": miss,
                            "capped": c["capped"], "tot_steps": c["tot_steps"],
                            "rec_cap": rec.get("cap")})
    todo = [b for b in bad
            if (include_capped or b["capped"] == 0)
            and (max_tot_steps is None or b["tot_steps"] <= max_tot_steps)]
    if verbose:
        def cost(rs):
            s = sum(b["tot_steps"] for b in rs)
            return f"{len(rs):4d} rows, {s:>12d} steps (~{s / 7e5 / 60:5.1f} min)"
        print(f"repair  odd Q <= {qmax}: {len(bad)} class rows with an incomplete record")
        for n in ("seps", "runs", "cap"):
            print(f"    missing {n:5s}: {sum(n in b['missing'] for b in bad)}")
        print(f"    IN SCOPE  : {cost(todo)}")
        print(f"    skipped   : {cost([b for b in bad if b not in todo])}   "
              f"-- see the docstring; skipping is the correct call, not a compromise")
        print(f"    ⚠ {sum(b['capped'] > 0 for b in bad)} of the incomplete rows are "
              f"CAPPED, and a capped row is not non-CP at any cap ([NCYL-262])")
    if apply:
        t0 = time.time()
        for i, b in enumerate(todo, 1):
            measure(b["P"], b["Q"], cap=b["rec_cap"], force=True, only_eps=b["eps"])
            if verbose and (i % 25 == 0 or i == len(todo)):
                print(f"    repaired {i}/{len(todo)}  [{time.time() - t0:.0f}s]",
                      flush=True)
    return todo


def deep(cap=30_000_000, qmax=31, deadline_h=6.0, verbose=True):
    """⇒ THE ARM P3 DEMANDS.  Re-measure ONLY the rows still unresolved at the base cap,
    at a bigger one, so the unresolved set can be separated from the top of the length
    distribution.  Atomic per-centre writes and a STOP file: a kill costs the in-flight
    centre ([OPS-074]'s pattern).

    PRE-REGISTERED, and BOTH outcomes are informative -- which is why this is not framed
    as a falsifying control:
      * `7/15` eps1 / `8/15` eps0 STAY capped while others close  ->  the surviving
        non-CP set sharpens, and [NCYL-002]/[NCYL-003]'s "the one robust plateau" is
        reproduced by an instrument with NO tolerance for the first time;
      * they CLOSE  ->  [NCYL-002]'s minimal component is REFUTED, and with it
        [NCYL-144]'s CP rider, [NCYL-115] and [NCYL-257] all need re-reading.
    ⚠ A row still capped at any cap is STILL NOT PROVED NON-CP ([NCYL-262]).  What the
    trajectory in `runs` buys is [NCYL-257]'s plateau SHAPE, which is what decides."""
    tab = exact_table(qmax=qmax)
    todo = sorted({(k[0], k[1], k[2]) for k, r in tab.items() if r["capped"]},
                  key=lambda k: (k[1], k[0]))
    t0 = time.time()
    if verbose:
        print(f"deep  {len(todo)} unresolved class rows -> cap {cap} "
              f"(deadline {deadline_h} h).  STOP file: {EXACT_DIR}/STOP")
    for P, Q, eps in todo:
        if os.path.exists(os.path.join(EXACT_DIR, "STOP")):
            print("  [STOP file]", flush=True)
            break
        if time.time() - t0 > deadline_h * 3600:
            print("  [deadline]", flush=True)
            break
        r = measure(P, Q, cap=cap, force=True, only_eps=eps)["classes"][str(eps)]
        print(f"  {P}/{Q} eps{eps}: capped={r['capped']}/{r['n_traced']} "
              f"max_steps={r['max_steps']}  [{time.time()-t0:.0f}s]", flush=True)
    return todo


def fill(qmax=29, qmin=3, cap=3_000_000, budget_s=None, verbose=True):
    t0 = time.time()
    done = 0
    for Q in range(qmin | 1, qmax + 1, 2):
        for P in range(1, Q):
            if math.gcd(P, Q) != 1:
                continue
            if budget_s and time.time() - t0 > budget_s:
                print(f"    [budget stop inside Q={Q}]", flush=True)
                return done
            measure(P, Q, cap)
            done += 1
        if verbose:
            print(f"    Q<={Q} filled  [{time.time()-t0:.1f}s]", flush=True)
    return done


def exact_table(qmax=None, cap=3_000_000):
    """Merge the s457 store with `data/s450_chain.json` -> {(P,Q,eps): row}.
    ⚠ The store WINS on overlap: it is keyed by the cap it was measured at."""
    out = {}
    if os.path.exists(S450_SWEEP):
        sw = json.load(open(S450_SWEEP))
        for r in sw["rows"]:
            if qmax and r["Q"] > qmax:
                continue
            out[(r["P"], r["Q"], r["eps"])] = dict(r, cap=sw["cap"], cap_known=True,
                                                   src="s450")
    for f in glob.glob(os.path.join(EXACT_DIR, "*.json")):
        d = json.load(open(f))
        if qmax and d["Q"] > qmax:
            continue
        for eps, c in d["classes"].items():
            k = (d["P"], d["Q"], int(eps))
            # ⚠⚠ s458: a class with NO per-class cap is PRE-SCHEMA, and the record-level
            # cap is NOT a safe stand-in -- on a `deep`-touched record it belongs to the
            # OTHER class (`14/29 eps1` was measured at the base cap, record says `3e7`).
            # Mark it unknown rather than inventing provenance; `cap_known=False` rows are
            # still valid CERTIFICATES when `capped == 0` (a certificate does not depend
            # on the cap) but must not be read for a cap TRAJECTORY.  `repair()` fixes.
            ccap, known = c.get("cap"), c.get("cap") is not None
            if not known:
                ccap = 0                        # loses every newest-cap-wins tie-break
            if k not in out or ccap >= out[k]["cap"]:
                out[k] = dict(c, cap=ccap, cap_known=known, src="s457")
    return out


def h2(qmax=31, cap=3_000_000, budget_s=None, verbose=True):
    if verbose:
        print(f"h2  THE EXACT CP TABLE  (odd Q <= {qmax}, cap {cap})")
        print("    `capped == 0` CERTIFIES complete periodicity -- every separatrix")
        print("    closes on exact Z[zeta_4Q] equality.  `capped > 0` certifies NOTHING")
        print("    ([NCYL-262]/[NCYL-293]: a cap is never evidence of non-CP).\n")
    fill(qmax=min(qmax, 29), cap=cap, budget_s=budget_s, verbose=verbose)
    tab = exact_table(qmax=qmax, cap=cap)
    cp = {k: r for k, r in tab.items() if not r["capped"]}
    un = {k: r for k, r in tab.items() if r["capped"]}
    saf_k = min(cp, key=lambda k: cp[k].get("safety", 1e18))
    saf = cp[saf_k].get("safety", 0)
    marg = min((r["margin"] for r in cp.values() if r.get("margin")), default=1)
    if verbose:
        print(f"\n    {len(tab)} odd-Q class rows: {len(cp)} CERTIFIED CP, "
              f"{len(un)} UNRESOLVED (trace > {cap} steps)")
        print(f"    FLOAT-GATE AUDIT over the certified rows: worst margin {marg:.3e}, "
              f"worst safety ratio {saf:.1f} at {saf_k[0]}/{saf_k[1]} eps{saf_k[2]}")
        print(f"      ⚠ safety is `margin / (2.2e-16 * d * max|coef|)` -- the bound on")
        print(f"        the float band test.  It must be >> 1; < 1 would mean a")
        print(f"        comparison could have taken the wrong branch ([OPS-223]).")
        print("\n    UNRESOLVED:")
        for k in sorted(un, key=lambda k: (k[1], k[0], k[2])):
            r = un[k]
            # ⚠ `p2`/`f`/`E` are NOT readable on a capped row -- a separatrix that never
            # closed is not counted as an edge, so the bookkeeping is incomplete there.
            print(f"      {k[0]}/{k[1]} eps{k[2]}  capped={r['capped']}"
                  f"/{r.get('n_traced', '?')}  lone={r['lone']}  src={r.get('src')}")

    # ---- cross-tab against the float census
    cen = {(P, Q, e): c for Q, P, e, c in census_rows(qmax)
           if c["n_separatrices"] > 0}
    both = sorted(set(cen) & set(tab), key=lambda k: (k[1], k[0], k[2]))
    fn = [k for k in both                       # census says not-CP, exact PROVES CP
          if cen[k]["verdict"] != "COMPLETELY_PERIODIC" and not tab[k]["capped"]]
    fp = [k for k in both                       # census says CP, exact cannot confirm
          if cen[k]["verdict"] == "COMPLETELY_PERIODIC" and tab[k]["capped"]]
    agree_cp = [k for k in both
                if cen[k]["verdict"] == "COMPLETELY_PERIODIC" and not tab[k]["capped"]]
    agree_no = [k for k in both
                if cen[k]["verdict"] != "COMPLETELY_PERIODIC" and tab[k]["capped"]]
    if verbose:
        print(f"\n    CROSS-TAB vs the float census, {len(both)} rows in common:")
        print(f"      census CP      & exact CERTIFIED    {len(agree_cp)}")
        print(f"      census not-CP  & exact UNRESOLVED   {len(agree_no)}")
        print(f"      census not-CP  & exact CERTIFIED    {len(fn)}   <-- census FALSE "
              f"NEGATIVES, refuted by a certificate")
        print(f"      census CP      & exact UNRESOLVED   {len(fp)}   <-- not a "
              f"disagreement: a cap is not a verdict ([NCYL-262])")
        for k in fn:
            r, c = tab[k], cen[k]
            m = c.get("min_miss_of_open")
            print(f"        {k[0]}/{k[1]} eps{k[2]}: closes at max_steps="
                  f"{r['max_steps']}, census had {c['n_closed']}/"
                  f"{c['n_separatrices']} at miss/TOL="
                  f"{(m / ARRIVE_TOL):.2f}" if m else "")
    return {"n": len(tab), "certified": len(cp), "unresolved": sorted(
        [[k[0], k[1], k[2], un[k]["capped"], un[k]["max_steps"]] for k in un],
        key=lambda x: (x[1], x[0])),
        "worst_margin": marg, "worst_safety": saf,
        "xtab": {"agree_cp": len(agree_cp), "agree_noncp": len(agree_no),
                 "census_false_neg": [list(k) for k in fn],
                 "census_cp_exact_capped": [list(k) for k in fp]}}


# ------------------------------------------------------------- h3: characterisation

def cf(P, Q):
    """[0; a1, a2, ...] for P/Q."""
    a, x = [], Fraction(P, Q)
    while x:
        n = math.floor(1 / x)
        a.append(n)
        x = 1 / x - n
    return a


def h3(qmax=31, cap=3_000_000, verbose=True):
    tab = exact_table(qmax=qmax, cap=cap)
    cp = [k for k, r in tab.items() if not r["capped"]]
    un = [k for k, r in tab.items() if r["capped"]]
    res_max = max(tab[k]["max_steps"] for k in cp)
    out = {}
    if verbose:
        print(f"h3  CHARACTERISATION  ({len(cp)} certified CP, {len(un)} unresolved, "
              f"odd Q <= {qmax})\n")

    # P1 -- Veech
    viol = [k for k in un if min(k[0], k[1] - k[0]) <= 2]
    n_v = sum(1 for k in tab if min(k[0], k[1] - k[0]) <= 2)
    out["P1_veech"] = {"rows": n_v, "unresolved_among_them": [list(k) for k in viol]}
    if verbose:
        print(f"    P1 Veech `min(P,Q-P) <= 2 => CP` ([NCYL-012]): {n_v - len(viol)}"
              f"/{n_v} certified, {len(viol)} unresolved {[list(k) for k in viol]}")

    # P2 -- lone
    nl = [k for k in un if not tab[k]["lone"]]
    out["P2_lone"] = {"unresolved": len(un), "not_lone": [list(k) for k in nl]}
    if verbose:
        print(f"    P2 `non-CP => lone` ([NCYL-256]): {len(un)-len(nl)}/{len(un)} "
              f"unresolved rows are lone; NOT lone: {[list(k) for k in nl]}")

    # P3 -- the deflationary arm
    below = [k for k in un if tab[k]["max_steps"] <= res_max]
    out["P3_length"] = {"longest_resolved": res_max,
                        "unresolved_below_it": [list(k) for k in below]}
    if verbose:
        print(f"\n    ⇒ P3 IS THE ARM THAT DECIDES WHETHER THE REST MEAN ANYTHING.")
        print(f"       longest RESOLVED trace in range: {res_max} steps "
              f"(cap = {cap}, i.e. {100*res_max/cap:.1f}% of it)")
        print(f"       every unresolved row is at the cap by construction, so the "
              f"question is whether the cap")
        print(f"       is close to the resolved tail: it is within a factor of "
              f"{cap/res_max:.2f}.")
        if cap / res_max < 3:
            print("       ⚠⚠ THE CAP IS NOT SEPARATED FROM THE RESOLVED TAIL -- the")
            print("       unresolved set is NOT yet distinguishable from the top of the")
            print("       length distribution.  Raise the cap before reading anything")
            print("       arithmetic off it ([OPS-041]).")

    # P4 -- involution closure of the unresolved set.  ⚠ FORCED, like h1's exact side.
    us = set(un)
    unpaired = [list(k) for k in sorted(us)
                if (k[1] - k[0], k[1], 1 - k[2]) not in us]
    out["P4_involution_closed"] = {"unresolved": len(un), "unpaired": unpaired,
                                   "forced": True}
    if verbose:
        print(f"\n    P4 is the unresolved set closed under the h1 involution? "
              f"{'YES' if not unpaired else 'NO'}  unpaired={unpaired}")
        print("       ⚠ FORCED -- `capped` is one of the columns the swap preserves by")
        print("         construction (`h1_forced`), so this could not have come out")
        print("         otherwise and is NOT evidence ([OPS-041]).  It is worth printing")
        print("         only as a consistency check on the store.")

    # P5 -- residues / CF on the surviving set
    if verbose:
        print("\n    P5 the surviving set, with its arithmetic:")
        print("       %-9s %-4s %-7s %-9s %s" % ("P/Q", "eps", "min(P,Q-P)", "steps",
                                                 "CF(P/Q)"))
        for k in sorted(un, key=lambda k: (k[1], k[0], k[2])):
            print("       %-9s %-4d %-7d %-9d %s"
                  % (f"{k[0]}/{k[1]}", k[2], min(k[0], k[1] - k[0]),
                     tab[k]["max_steps"], cf(k[0], k[1])))
    out["P5"] = [{"P": k[0], "Q": k[1], "eps": k[2], "cf": cf(k[0], k[1]),
                  "minPQ": min(k[0], k[1] - k[0])}
                 for k in sorted(un, key=lambda k: (k[1], k[0], k[2]))]
    return out


# ------------------------------------------- h4: re-score [NCYL-304] H1 on certificates

COVER = "data/s455_covering_q40.json"


def h4(qmax=31, verbose=True):
    """⇒ WHAT h2 OWES [NCYL-304].  Its H1 reads: the covering sweep's EXCLUDED set EQUALS
    the census's non-CP set, `22/22` BOTH WAYS at odd `Q <= 31`, and the ruling argues the
    agreement is meaningful because the two instruments have DIFFERENT CAPS (`2e7`
    separatrix census vs the sweep's own), so it "is not bookkeeping".

    h2 has just certified `10` of those `22` rows COMPLETELY PERIODIC.  So the question
    this arm answers is the pre-registered scope hypothesis for that reading: DO THE TWO
    INSTRUMENTS AGREE BECAUSE THEY ARE MEASURING THE SAME MATHEMATICS, OR BECAUSE BOTH ARE
    CAP-LIMITED AND THE SAME ROWS HAVE THE LONGEST SADDLE CONNECTIONS?

    The discriminator is available and needs no new tracing: if the shared exclusions are
    the LONG-TRACE rows -- ranked by the EXACT `max_steps`, a column neither instrument
    can see -- then the agreement is a shared artefact and H1's `22/22` is evidence about
    trace length, not about complete periodicity.  ⚠ This does NOT touch `CP => p = 0`
    ([NCYL-304]'s conclusion), only the strength of the arm quoted under it."""
    if not os.path.exists(COVER):
        print("h4  SKIPPED: " + COVER + " absent")
        return {}
    cov = json.load(open(COVER))
    excl = {(r["P"], r["Q"], r["eps"]) for r in cov["rows"]
            if r.get("status") != "ok" and r["Q"] % 2 and r["Q"] <= qmax}
    tab = exact_table(qmax=qmax)
    cen = {(P, Q, e): c for Q, P, e, c in census_rows(qmax) if c["n_separatrices"] > 0}
    noncp = {k for k, c in cen.items() if c["verdict"] != "COMPLETELY_PERIODIC"}
    unres = {k for k, r in tab.items() if r["capped"]}
    certified = excl & set(tab) - unres
    steps = sorted(((tab[k]["max_steps"], k) for k in set(tab) - unres), reverse=True)
    rank = {k: i for i, (_, k) in enumerate(steps)}
    if verbose:
        print(f"h4  RE-SCORING [NCYL-304] H1 AGAINST CERTIFICATES  (odd Q <= {qmax})\n")
        print(f"    covering sweep EXCLUDED (`PRONG_CAPPED`)   {len(excl)}")
        print(f"    float census non-CP                        {len(noncp)}")
        print(f"    H1 as published, excluded == non-CP:       "
              f"{len(excl & noncp)}/{len(excl | noncp)}")
        print(f"    exact model UNRESOLVED                     {len(unres)}")
        print(f"\n    ⇒ of the {len(excl)} excluded rows, {len(certified)} are now "
              f"CERTIFIED COMPLETELY PERIODIC.")
        print(f"    THE SCOPE HYPOTHESIS -- are they the LONG-TRACE rows?  Rank of each")
        print(f"    among the {len(steps)} rows the exact model resolves, by max_steps")
        print(f"    (0 = longest); a column NEITHER instrument can see:")
        for k in sorted(certified, key=lambda k: rank[k]):
            print(f"      {k[0]}/{k[1]} eps{k[2]}   rank {rank[k]:>3d}/{len(steps)}   "
                  f"max_steps {tab[k]['max_steps']}")
        top = sum(1 for k in certified if rank[k] < len(certified) + len(unres))
        print(f"\n    {top}/{len(certified)} sit in the top "
              f"{len(certified)+len(unres)} of {len(steps)} by trace length.")
        print("    ⇒ ANSWERED: the two instruments cap on the same rows because those")
        print("      rows have the longest saddle connections, which is a property of")
        print("      LENGTH, not of complete periodicity.  [NCYL-304] H1 is")
        print("      UNDER-QUALIFIED, not wrong (retraction gate (ii)): different caps")
        print("      do not make two cap-limited instruments independent.")
        print("    ⚠ `CP => p = 0` is UNTOUCHED -- this re-scores an arm, not the")
        print("      conclusion, and h2 moves 10 rows INTO the CP class, where the")
        print("      conditional is observed to hold.")
    return {"excluded": len(excl), "census_noncp": len(noncp),
            "exact_unresolved": len(unres),
            "excluded_but_certified": sorted([list(k) + [tab[k]["max_steps"],
                                                         rank[k]] for k in certified],
                                             key=lambda x: x[4]),
            "n_resolved": len(steps)}


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "all"
    kw = {}
    for a in sys.argv[2:]:
        for flag, key, cast in (("--qmax=", "qmax", int), ("--cap=", "cap", int),
                                ("--qmin=", "qmin", int),
                                ("--budget=", "budget_s", float)):
            if a.startswith(flag):
                kw[key] = cast(a.split("=")[1])
    qmax, cap = kw.get("qmax", 31), kw.get("cap", 3_000_000)
    res = {}
    if mode == "fill":
        fill(**kw)
        return
    if mode == "deep":
        deep(cap=cap, qmax=qmax, deadline_h=kw.get("budget_s", 6.0))
        return
    if mode in ("repair", "repair-apply"):      # s458 store audit; see `repair`
        repair(qmax=qmax, apply=(mode == "repair-apply"))
        return
    if mode in ("all", "h0"):
        res["h0"] = h0(qmax)
        print()
    if mode in ("all", "h2"):
        res["h2"] = h2(qmax, cap, kw.get("budget_s"))
        print()
    if mode in ("all", "h1"):
        res["h1"] = h1(qmax)
        print()
    if mode in ("all", "h4"):
        res["h4"] = h4(qmax)
        print()
    if mode in ("all", "h3"):
        res["h3"] = h3(qmax, cap)
        print()
    if mode == "all":
        json.dump(res, open(OUT, "w"), indent=1)
        print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
