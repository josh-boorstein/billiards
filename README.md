# Billiards in rational right triangles — code deposit

Exact-arithmetic code supporting three papers on rational right triangles:

- **A minimal component in a rational right triangle: an exact decomposition, and a
  conjecture of Boshernitzan** — `lib/` + `minimal_component/`
- **The orphan theorem: degenerate perpendicular orbits in rational right triangles** —
  `engine/` + `probes/` (legacy layout, see below)
- **The necklace of a rational right triangle: exact counting and integrality in genus
  zero** — `engine/` + `probes/` (legacy layout)

All three are in preparation. `MANIFEST.md` maps each script to the statement it supports.

**The repository is being rebuilt paper by paper** into a small core library (`lib/`) and
one directory per paper holding one script per statement, each named for what it checks.
The minimal-component paper is done. The other two still use the legacy layout — the
research code as it ran, with session-numbered script names — and move over in turn; the
legacy tree is removed when both have.

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

## The legacy layout (orphan and necklace papers)

Every count in the orphan paper is an exact algebraic computation, never a numerical
estimate, and every one is gated by a certificate that the partition it counts is
complete. Appendix A names four ingredients, and all four are here:

| Appendix A | What it is | Where |
|---|---|---|
| A.1 | the exact partition tree, symbolic in `cos 2α` over `ℤ[½][cos 2α]` | `engine/cos_poly.py`, `engine/exact_state.py`, `engine/word_tree*.py` |
| A.2 | cyclotomic deduplication — comparison by **value** in `ℤ[ζ_{2Q}]`, never by coefficient vector | `engine/n_exact.py` (`exact_key`) |
| A.3 | the completeness certificate that gates every count | `engine/n_exact.py` (`n_pq_certified`) |
| A.4 | the fixed-angle backend and the cylinder data | `engine/ring_cyc.py`, `engine/exact_state_ring.py`, `engine/cyl_diagram.py` |

There is no floating-point tracing in the construction of the partition. Boundary
positions are converted to floating point only for display, at the very end.

**A.2 and A.3 are the two places where an obvious implementation is silently wrong.**
The numbers `cos(kPπ/Q)` are `ℚ`-linearly dependent, so two equal boundaries can carry
different coefficient vectors and a coefficient-wise comparison over-counts. And odd-`P`
counts *plateau*: `n(3/29)` sits at `4` well past depth `10³` before climbing to `19`,
so any "stable over two depths" heuristic accepts the wrong answer. Only the truncation
certificate is reliable.

## Running it

Python 3.13. There is no package and no `__init__.py` — modules import one another by
**bare name**, resolved by putting both source directories on the path. The layout is
preserved from the research repository deliberately, so that the deposited scripts are
byte-identical to the ones that produced the results rather than a repackaging of them.

```sh
python3 -m venv .venv && .venv/bin/python -m pip install -r requirements.txt
export PYTHONPATH=engine:probes
```

Two drivers reproduce the fast end of each paper:

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
engine/      18 modules — the reusable exact-arithmetic core
probes/      58 modules — the per-result scripts, named by the session that wrote them
reproduce/   drivers that regenerate the headline numbers
data/        output directory (created on demand)
MANIFEST.md  script → paper statement
```

`engine/` and `probes/` are the transitive import closure of the scripts behind the two
papers; all 76 modules import cleanly, and nothing outside the closure is included.

## Scope, and what the code does not establish

The scripts compute; they do not prove, and the papers are explicit about which is which.
Three points matter for anyone reading the output:

- **`n(2/Q) = Q − 1` (Theorem C) is VERIFIED, not proved.** It is exact on every odd
  `Q ≤ 45` and there is no proof of the general case. Running `reproduce/orphan_counts.py`
  reproduces the evidence, not a proof. What remains open is a surjectivity: whether the
  leg meets every cylinder.
- **The necklace paper's arithmetic layer rests on one unproved input, (G0b)**, and
  several of its statements are conditional on complete periodicity. Running the checks
  does not discharge either hypothesis; the paper tags each statement with its own scope.
- **Coverage is finite and is quoted where it matters.** A check reported as
  `6040/6040 class rows` means exactly that range and no more.

Script docstrings carry cross-references to the research project's internal ledgers
(bracketed IDs such as `[NCYL-283]`, and filenames such as `future_directions.md`). Those
ledgers are not part of this deposit. The docstrings are reproduced unedited rather than
rewritten, so that what is deposited is the code as it ran.

## License

MIT — see `LICENSE`.
