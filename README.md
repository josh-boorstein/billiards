# Billiards in rational right triangles — code deposit

Exact-arithmetic code supporting four papers on right triangles — three on rational right
triangles, and one on right triangles whose angle is an irrational multiple of `π`:

- **The orphan theorem: degenerate perpendicular orbits in rational right triangles** —
  `lib/` + `orphan/`
- **The necklace of a rational right triangle: exact counting and complete periodicity** —
  `lib/` + `neck/`
- **No arrivals: complete periodicity and degenerate orbits in rational right triangles** —
  `lib/` + `method/`
- **Non-perpendicular periodic orbits in mirror directions of right triangles** —
  `lib/` + `mirror_A/`

All four are in preparation. `MANIFEST.md` maps each script to the statement it supports.

The repository is a small core library (`lib/`) and one directory per paper holding one script
per statement, each named for what it checks.

```sh
python3 -m venv .venv && .venv/bin/python -m pip install -r requirements.txt
```

No `PYTHONPATH` is needed: each script puts `lib/` on the path itself. Every script prints
the statements it checks, one line each, and exits nonzero if any fails.

`lib/` begins with five modules; the paper sections below add the rest.

| Module | What it is |
|---|---|
| `cyclofield.py` | exact arithmetic in `ℚ(ζ_N)`; the sign of a real element is read off a certified interval enclosure, zero off the reduced polynomial |
| `triangle.py` | the triangle and the exact unfolding of an orbit along a fixed word: positions, hitting times and path lengths are affine in the launch parameter with coefficients in `ℚ(ζ_4Q)`, so a word's validity region is one interval with exact endpoints, and a tiling of a segment by such intervals is a complete certificate |
| `closure.py` | the closure certificate: every separatrix of the genus-zero base walked to its end, arrival decided by exact equality in `ℤ[ζ_4Q]` |
| `develop.py` | the formal development of a word, giving a return word's translation as a function of the angle |
| `tracer.py` | a 60-digit floating-point tracer, used only to **propose** words (which the exact code then verifies or rejects) |

**Notation.** The scripts write `eps` for the direction-class label that the necklace paper and
the method paper print as `𝔢` (the class `8/15 eps 0` is the paper's `8/15`, `𝔢 = 0`), and
`(E)`, `(G0b)` in older docstrings for the hypotheses the papers print as (OB) and (ZZ).

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

## The method paper: `lib/` and `method/`

```sh
.venv/bin/python method/run_all.py          # about ten minutes
.venv/bin/python method/run_all.py --full   # Proposition A.23's census recomputed to Q <= 80 (about half an hour more)
```

One script per computation the paper's *Code availability* names, built on `lib/` alone (and
`mpmath`, `sympy`, `scipy`):

- **Proposition 6.4 and §6.3** (the minimal component at `T₀ = 8/15`, and the arrival `O → A` on
  its boundary): `s6_4_*.py`, `s6_3_*.py`, on `t0.py`'s triangle and words. Every statement about
  `T₀` is decided in `ℚ(ζ₆₀)` with no tolerance — the interval lengths, both cuts, the seven
  validity intervals and their tiling, the six branch intervals and their tiling of `L1`, the
  cylinder table and the Kac mass. The words come from the tracer, but a missed word leaves a gap
  in a tiling and a wrong one has an empty interval, so neither can pass. Arms marked *(beyond
  III)* check figures the paper does not print; they come from the long-form write-up of this
  decomposition that the paper absorbed.
- **Proposition A.12 and Lemma A.15**: exact rational arithmetic, every case of the finite check.
- **Proposition A.23** (the census at `Q ≤ 80`): a functional is found by a floating-point LP
  (`scipy`) and then verified by an exact integer Bellman–Ford computation; the LP's word is never
  used. `sA_23_record.json` holds every case; the default run recomputes `Q ≤ 40` against it.
- **Theorem 6** (the counterexamples to Conjecture 2.2(a)): the words are proposed by a 50-digit
  tracer at one point of each interval; every constraint, tiling join, slope and translation is
  then an exact Laurent polynomial in `e^{iπx}` over `ℤ[i]`, each identity is checked exactly and
  each inequality by interval arithmetic (`mpmath.iv`) on the whole interval; the Galois step is
  exact (`sympy`).
- **§B.4** (the `(4, 3)` ladder): the fourteen inequalities in outward-rounded interval arithmetic
  at every `3 ≤ m ≤ 3000`, and by an exact Taylor bound with Lagrange remainder for
  `θ ≤ 1/7600`; `sB_4_record.json` holds the fourteen numerators, orders and margins. At
  `m = 1, 2` every `R`-prong is traced in interval arithmetic.

Theorems 2–5 and the rest of the paper are proofs and have no script.

## The mirror-direction paper: `lib/` and `mirror_A/`

```sh
.venv/bin/python mirror_A/run_all.py          # under two minutes
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
rational functions of `t`). Remark 6.10's period-`944` cylinder at `cos 2a = (√7 − 1)/4`, in the
class `⊥OA` inside the range `π/6 < a < π/4` of Troubetzkoy's claim (the paper's [Tro05], abstract (iii)), is certified the same way
(`tro05_record.json`). The tiles of §6 are exact (sympy, real roots). The figures are float
drawings whose captioned numbers are checked first.

## Layout

```
lib/          the core library (above)
orphan/       the orphan paper
neck/         the necklace paper
method/       the method paper
mirror_A/     the mirror-direction paper
MANIFEST.md   script → paper statement
```

## Scope, and what the code does not establish

The scripts compute; they do not prove, and the papers are explicit about which is which.
Two points matter for anyone reading the output:

- **Several of the necklace paper's statements are conditional** — on complete periodicity,
  or on the existence of an orphan branch (hypothesis (OB)). Running the checks does not
  discharge either hypothesis; the paper tags each statement with its own scope.
- **Coverage is finite and is quoted where it matters.** A check reported as
  `6040/6040 class rows` means exactly that range and no more.

Script docstrings carry cross-references to the research project's internal ledgers
(bracketed IDs such as `[NCYL-283]`, and filenames such as `future_directions.md`). Those
ledgers are not part of this deposit.

## License

MIT — see `LICENSE`.
