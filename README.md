# Billiards in rational right triangles — code deposit

Exact-arithmetic code supporting four papers on right triangles — three on rational right
triangles, and one on right triangles whose angle is an irrational multiple of `π`:

- **A minimal component in a rational right triangle: an exact decomposition, and a
  conjecture of Boshernitzan** — `lib/` + `minimal_component/`
- **The orphan theorem: degenerate perpendicular orbits in rational right triangles** —
  `lib/` + `orphan/`
- **The necklace of a rational right triangle: exact counting and complete periodicity** —
  `lib/` + `neck/`
- **Non-perpendicular periodic orbits in mirror directions of right triangles** —
  `lib/` + `mirror_A/`

All four are in preparation. `MANIFEST.md` maps each script to the statement it supports.

**The repository is being rebuilt paper by paper** into a small core library (`lib/`) and
one directory per paper holding one script per statement, each named for what it checks.
All four papers are done. The legacy layout — the research code as it ran, with
session-numbered script names — remains only until a fifth paper's directory is added; it is
then removed.

## The minimal-component paper: `lib/` and `minimal_component/`

```sh
python3 -m venv .venv && .venv/bin/python -m pip install -r requirements.txt
.venv/bin/python minimal_component/run_all.py          # about 3 minutes
.venv/bin/python minimal_component/run_all.py --full   # recomputes the census, about an hour
.venv/bin/python minimal_component/figures.py          # the four figures -> minimal_component/figs/
```

No `PYTHONPATH` is needed: each script puts `lib/` on the path itself. Every script prints
the statements it checks, one line each, and exits nonzero if any fails.

`lib/` holds five modules, and everything in the paper is computed from them:

| Module | What it is |
|---|---|
| `cyclofield.py` | exact arithmetic in `ℚ(ζ_N)`; the sign of a real element is read off a certified interval enclosure, zero off the reduced polynomial |
| `triangle.py` | the triangle and the exact unfolding of an orbit along a fixed word: positions, hitting times and path lengths are affine in the launch parameter with coefficients in `ℚ(ζ_4Q)`, so a word's validity region is one interval with exact endpoints, and a tiling of a segment by such intervals is a complete certificate |
| `closure.py` | the closure certificate of §5: every separatrix of the genus-zero base walked to its end, arrival decided by exact equality in `ℤ[ζ_4Q]` |
| `develop.py` | the formal development of a word, giving a return word's translation as a function of the angle |
| `tracer.py` | a 60-digit floating-point tracer, used only to **propose** words (which the exact code then verifies or rejects) and to trace the deformed flow of §3.1 |

**What is exact and what is not.** Every statement of §§2–4 and §6 about `T₀` is decided in
`ℚ(ζ₆₀)`, with no tolerance: the interval lengths, both cuts, the seven validity intervals
and their tiling of `Σ_H̄`, the six branch intervals and their tiling of `L1`, the cylinder
table, the Kac integral and Proposition 5. The words themselves come from the tracer, but
a missed word leaves a gap in a tiling and a wrong one has an empty interval, so neither
can pass. The deformed-angle statements of §3.1 compare an exact formula with a 60-digit
trace, as the paper says. The census of §5 is exact per row and **one-sided**: a capped row
is no verdict. `minimal_component/census_record.json` holds every class row's verdict for
odd `Q ≤ 39` and even `Q ≤ 40`; `s5_census.py` recomputes any of them.

## The orphan paper: `lib/` and `orphan/`

```sh
.venv/bin/python orphan/run_all.py          # minutes; counts read from orphan/partition_record.json
.venv/bin/python orphan/run_all.py --full   # every count recomputed from scratch (hours)
.venv/bin/python orphan/figures.py          # the four figures -> orphan/figs/
```

Two modules join `lib/`:

| Module | What it is |
|---|---|
| `partition.py` | the perpendicular partition of a side at a fixed centre, **exactly**: the development is carried in `ℤ[ζ_4Q]` without division (a reflection is an index gather and a reduction modulo `Φ_4Q`), each step is ONE comparison — the beam height against the height of the copy's opposite vertex — decided exactly (a float only when its own error bound is far below the gap, a certified evaluation otherwise), and there is no width floor. The first-return map `J` is read off the development and matched endpoint for endpoint as field elements. Beams perpendicular to any side; resumable walk state. |
| `trans_surface.py` | a floating-point model of a translation surface glued from polygons, cut on separatrices and flowed; used only for the necklace paper's double regular `Q`-gon figure (`neck/figures_5_6.py`) |

**What is exact and what is not.** Every branch count, width, boundary type, orphan and
pairing is `partition.py`'s and exact. This is stronger than the paper's Appendix A, which
describes the partition tree the counts were first computed on — there side order and orbit
closure are decided numerically, under a `10⁻¹¹` width floor — and the two engines agree on
every centre where both have run (the research repository's gate: every complete centre of
its branch store, cell by cell). The surface of §§6–7 is BUILT from the gluing
(`s6_7_surface.py`), not read off the formulas it checks. The folded orbits and developments
in the figures are 60-digit traces, and each must reproduce the exact word before it is
drawn; Figure 4's leaves are float traces from three heights of the exact orphan branch, and
every claim of its caption is asserted before it is drawn. The paper cites [DR22] for the
minimal component at `8/15` (§2.1) and proves nothing about it; the necklace paper's check,
`neck/s11_cp_band.py`, reproduces it exactly with `lib/flow.py`.

`orphan/partition_record.json` holds the census of every centre the paper's ranges need,
computed by `partition.py` on a resumable step ladder. The scripts read it and re-verify a
slice of it from scratch on every run; `--full` recomputes all of it.

## The necklace paper: `lib/` and `neck/`

```sh
.venv/bin/python neck/run_all.py          # about half an hour; records re-verified on a slice
.venv/bin/python neck/run_all.py --full   # every record recomputed, every census row at its closing cap (hours)
.venv/bin/python neck/figures.py          # Figures 1-4 -> neck/figs/
.venv/bin/python neck/figures_5_6.py      # Figures 5-6 -> neck/figs/
```

Three modules join `lib/`:

| Module | What it is |
|---|---|
| `necklace.py` | the necklace readout: `closure.py`'s exact walks of every separatrix of the genus-zero base, with their incidences with the fixed set of the reversing involution recorded (the counts `h`, `v`, `f`, `D`, `E` of §§3 and 7); and the interval model of §4 |
| `flow.py` | the exact flow decomposition of Appendix A.2: a rational-triangle billiard in any direction of the `π/2d` grid decomposed into cylinders and components without periodic trajectory — a port of the decision procedure of flatsurf's `intervalxt` (Zorich-accelerated Rauzy induction, saddle-connection merges, a Boshernitzan certificate verified in integers; Kac's mass balance asserted exactly on every run). `n_flow(P, Q)` reads `n(P/Q)` off it: every perpendicular beam retraces, so a branch is one crossing of an `L1` edge by a cylinder, certified by the cylinders' exact tiling of that edge; `beam='h'` does the same for the hypotenuse |
| `partition.py`, `trans_surface.py` | as for the orphan paper (above) |

**What is exact and what is not.** Every branch count is certified twice over where both engines
have run — by `partition.py` and by `flow.py` — and the two agree on every centre. Complete
periodicity is `flow.py`'s, and both of its outcomes are certificates: a cylinder with its exact
period, or a component with a verified Boshernitzan functional. The certificates of Appendix A.3
are exact real-root isolations (sympy) on the printed conditions (`sA3_conditions.py`) and on the
mechanical prover's own (`window_prover.py`). `neck/band_record.json` holds the flow decomposition
of every coprime centre with `Q ≤ 60` (perpendicular to `L1`, and perpendicular to `H` at even `Q`
and at odd `Q ≤ 31`); `export_band_record.py` wrote it from the research repository's sweeps, and
`s11_cp_band.py` recomputes a slice of it on every run and all of it with `--full`. The figures
are drawn from exact data except Figure 5's polygon model (a float model whose closed forms are
checked to `1e-9`), and every caption claim is asserted before drawing.

## The mirror-direction paper: `lib/` and `mirror_A/`

```sh
.venv/bin/python mirror_A/run_all.py          # about a minute
.venv/bin/python mirror_A/figures.py          # the four figures -> mirror_A/figs/
```

The paper's §5 exhibits periodic cylinders at triangles with `cos 2a` algebraic and `a/π`
irrational; each is certified here **from its word alone** by
`lib/mirror.py`, on top of `lib/numfield.py`:

| Module | What it is |
|---|---|
| `numfield.py` | exact arithmetic in a real number field `ℚ(θ)`, `θ` a chosen real root of an irreducible integer polynomial, with certified signs (interval evaluation on a refining enclosure of `θ`); and its CM extension `K(w)`, `w = e^{ia}`, `t = w + 1/w = 2 cos a`, in which `Im` signs are exact |
| `mirror.py` | the triangle in `ℚ(2 cos a)(w)` and the cylinder certificate: the word unfolds to a translation exactly parallel to a mirror direction (whose class is read off, not assumed), and the open strip of parallel lines separating every crossed edge strictly is non-empty — its two ends are exact field elements, so the cylinder's width is exact; perpendicularity is decided exactly; `a/π ∉ ℚ` is decided from the minimal polynomial |

**What is exact and what is not.** Every certificate is decided in `ℚ(2 cos a)` with no
tolerance; a float is used only to propose a candidate direction index, which is then checked
exactly. The words come from the research searches' tracers and are stored in the records
(`table1_record.json`, `sqrt41_record.json`); how a word was found does not enter its
certificate. The family of §5.3 is certified exactly at its members `n ≤ 8`, which is all the
paper claims; the same script also carries the symbolic proof that it continues for every `n` (Sturm sequences over `ℚ` on `274`
rational functions of `t`). The tiles of §6 are exact (sympy, real roots). The figures are float
drawings whose captioned numbers are checked first.

## The legacy layout

`engine/`, `probes/` and `reproduce/` are the research code as it ran, kept until the next
paper's directory is added. No per-paper directory imports them. The necklace paper's
Appendix A.1 names one ingredient that lives only here, the fixed-angle backend in `ℤ[ζ_{2Q}]`
that computed §5.1's cylinder data: `engine/ring_cyc.py`, `engine/exact_state_ring.py`,
`engine/cyl_diagram.py`. The rest of this tree — the orphan paper's original partition engine,
symbolic in `cos 2α` (`engine/cos_poly.py`, `engine/exact_state.py`, `engine/word_tree*.py`,
`engine/n_exact.py`), superseded by `lib/partition.py` — is described in the research
repository's notes and is not needed to check any statement of the four papers.

⚠ That engine decides side order and orbit closure in floating point; what is exact is the value
carried at each node. And two places in it are where an obvious implementation is silently
wrong: the numbers `cos(kPπ/Q)` are `ℚ`-linearly dependent, so two equal boundaries can carry
different coefficient vectors and a coefficient-wise comparison over-counts (`n_exact.exact_key`
compares by value); and odd-`P` counts *plateau* — `n(3/29)` sits at `4` well past depth `10³`
before climbing to `19` — so only a truncation certificate (`n_exact.n_pq_certified`) is reliable.

## Running the legacy scripts

Python 3.13. The per-paper directories need no setup. In the legacy tree there is no package
and no `__init__.py` — modules import one another by
**bare name**, resolved by putting both source directories on the path. The layout is
preserved from the research repository deliberately, so that the deposited scripts are
byte-identical to the ones that produced the results rather than a repackaging of them.

```sh
python3 -m venv .venv && .venv/bin/python -m pip install -r requirements.txt
export PYTHONPATH=engine:probes
```

Two drivers reproduce what the research code computed first:

```sh
.venv/bin/python reproduce/orphan_counts.py            # ~5s;  --max-q 31 for more
.venv/bin/python reproduce/orphan_counts.py --slow     # adds the n(3/29) plateau check
.venv/bin/python reproduce/necklace_checks.py          # ~10s
```

`reproduce/orphan_counts.py` reproduces `n(2/Q) = Q − 1` on odd `Q` and, with `--slow`,
the plateau trap that motivates the certificate. `reproduce/necklace_checks.py` runs the
counting identity of §3 and the necklace prong census of §§2–3. Individual scripts run
directly (`python3 probes/<name>.py`) once `PYTHONPATH` is set; some write JSON to
`data/` and some take minutes.

## Layout

```
lib/, minimal_component/   the minimal-component paper (above)
lib/, orphan/              the orphan paper (above)
lib/, neck/                the necklace paper (above)
lib/, mirror_A/            the mirror-direction paper (above)
engine/      18 modules — legacy: the research code's exact-arithmetic core
probes/      58 modules — legacy: per-result scripts, named by the session that wrote them
reproduce/   legacy: drivers for the research code's headline numbers
data/        output directory (created on demand)
MANIFEST.md  script → paper statement
```

`engine/` and `probes/` are the transitive import closure of the legacy scripts; all 76
modules import cleanly, and nothing outside the closure is included.

## Scope, and what the code does not establish

The scripts compute; they do not prove, and the papers are explicit about which is which.
Three points matter for anyone reading the output:

- **Several of the necklace paper's statements are conditional** — on complete periodicity,
  or on the existence of an orphan branch (hypothesis (E)). Running the checks does not
  discharge either hypothesis; the paper tags each statement with its own scope.
- **Coverage is finite and is quoted where it matters.** A check reported as
  `6040/6040 class rows` means exactly that range and no more.

Script docstrings carry cross-references to the research project's internal ledgers
(bracketed IDs such as `[NCYL-283]`, and filenames such as `future_directions.md`). Those
ledgers are not part of this deposit. The docstrings are reproduced unedited rather than
rewritten, so that what is deposited is the code as it ran.

## License

MIT — see `LICENSE`.
