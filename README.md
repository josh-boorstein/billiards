# Billiards in rational right triangles — code deposit

Exact-arithmetic code supporting two papers on the perpendicular-beam partition of a
rational right triangle:

- **The orphan theorem: degenerate perpendicular orbits in rational right triangles**
- **The necklace of a rational right triangle: exact counting and integrality in genus zero**

Both are in preparation. This repository is the deposit promised in the orphan paper's
Appendix A ("Code availability"). `MANIFEST.md` maps each script to the statement it
supports.

## What the code does

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
