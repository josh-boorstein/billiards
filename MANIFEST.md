# Manifest — script to paper statement

Which script supports which statement. Each script's own docstring carries its
pre-registration, its arms and its controls: **the manuscript presents, the docstring is
the source.** Where a docstring and this table disagree, the docstring wins.

The per-paper scripts need no setup. The legacy scripts (`engine/`, `probes/`, `reproduce/`) run with `PYTHONPATH=engine:probes` set; some write JSON to `data/`.

---

## The minimal-component paper (`minimal_component/`)

Built on `lib/` alone. Each script's docstring lists the statements it checks.

| Statement | Script |
|---|---|
| §2 — the 30 directions met, the phase space and `area(B)`; the field `K`, its basis and determinant `8`; `cot α ∉ K` (§3.3) | `s2_phase_space.py` |
| §3.1 — the first return to `Σ_H̄` is a three-interval exchange: the seven words, the exact translations, `λ₁, λ₂, λ₃` and their coordinates, the boundary table | `s3_1_return_map.py` |
| §3.1 — *why boundary 5 is different*: `V(α)`, the `12°` identity, the rates `85.21` and `114.7`, the trace of the deformed flow; Figure 2's data | `s3_1_boundary5.py` |
| §3.2 — Lemma 3 (rank `3`), Lemma 4, the tower, the roof's four values | `s3_2_minimality.py` |
| §3.3 — the cylinder table, `15 g_len`, the Kac integral `μ(M)`, **Proposition 5** | `s3_3_mass_identity.py` |
| §3.4 — first-return condition, the exact tiling of `[0,1]` by seven intervals, the mass budget | `s3_4_completeness.py` |
| §4 — the second route: `n(8/15) = 6` by an exact tiling of `L1`, three equal-width pairs | `s4_branch_count.py` |
| §5 — every direction of `S¹(2d′)` is an image of a perpendicular direction (`Q ≤ 120`) | `s5_coverage.py` |
| §5, abstract, §6 — the census: the table at odd `Q ≤ 31`, the per-`Q` counts, *beyond the table*, 2.2(e) at `d ≤ 50`, `T₀`'s two class rows | `s5_census.py` (+ `census_record.json`) |
| §6 — the saddle connection `O → A`, its arrival identity, the `3`/`5`-family test, the leaves either side, its rate under deformation | `s6_connection.py` |
| Appendix — the seven return words and the six cylinder half-words | `appendix_words.py` |
| Figures 1–4 | `figures.py` |

⚠ **Scope.** The census is exact per row and one-sided: `capped == 0` proves complete
periodicity, a capped row is no verdict, and the undecided rows are listed by name. The
failure at `T₀` is not certified by the census at all; it is §3's result.

---

## The mirror-direction paper (`mirror_A/`)

*Non-perpendicular periodic orbits in mirror directions of right triangles.* Built on `lib/`
alone (`numfield.py`, `mirror.py`). Each script's docstring lists the statements it checks.

| Statement | Script |
|---|---|
| §5.1 — Examples 5.1 and 5.2: words, periods, classes, the translation `4e(2a)`, launch intervals and the width `√7/56` exactly, smallest angles; controls | `s5_1_examples.py` |
| §5.2 — Table 1: the `14` cylinders at rational `cos 2a`, `q ≤ 20`, their classes and periods; Niven; the `1 : √(4m − 1)` rows | `s5_2_rational_table.py` (+ `table1_record.json`) |
| §5.3 — the triangle `cos 2a = (1 + √41)/8`: its four cylinders, one per class, exactly; controls | `s5_3_sqrt41.py` (+ `sqrt41_record.json`) |
| §5.3 — the family `y_n → y_c`: the members `n ≤ 8` (period `30n + 6`) and the period-`30` cylinder at `y_c`, exactly — what the paper claims; beyond the paper, the symbolic checks that the family continues for every `n` (irreducibility and `a_n/π ∉ ℚ`, the alignment identity in `m`, the `274`-function clearance certificate, formal non-perpendicularity) | `s5_3_family.py` |
| §6.3 — Proposition 6.6 and Table 2: the five tiles, their end polynomials, the cover, the leg swaps, the five `F_W` and their non-degeneracy; Lemmas 6.3 and 6.4 at one angle in each tile | `s6_3_tiles.py` |
| Figures 1–4 | `figures.py` |

⚠ **Scope.** A certificate proves that the cylinder it names exists. Table 1's blank cells are
searches that found nothing, not proofs of absence, and nothing here checks them. §§2–4,
§§6.1–6.2, 6.4 and §7 are proofs in the paper and have no script.

---

## The legacy exact-arithmetic core (`engine/`)

The research code's core, kept until the next paper's directory is added; no per-paper script
imports it. The necklace paper's Appendix A.1 backend is `ring_cyc.py`, `exact_state_ring.py`
and `cyl_diagram.py`; the rest is the orphan paper's original partition engine, superseded by
`lib/partition.py`.

| Module | Role |
|---|---|
| `cos_poly.py` | `CosPoly` — an element of `ℤ[½][cos 2α]` as an integer coefficient vector |
| `cos_poly_fast.py` | numpy-backed drop-in for `cos_poly` |
| `exact_state.py` | `ExactTriangleState` — the symbolic (multi-angle) triangle state |
| `exact_state_fast.py` | drop-in for `exact_state` over `cos_poly_fast` |
| `word_tree.py` | the reflection-word partition tree |
| `word_tree_fast.py` | the same tree, faster node expansion |
| `word_tree_nofrac.py` | fully algebraic split-ordering — no float tracer in the tree |
| `n_exact.py` | **`n_pq_certified`** (A.3, the completeness certificate), **`exact_key`** (A.2, cyclotomic deduplication), `cyclotomic_poly` |
| `ring_cyc.py` | the cyclotomic ring `ℤ[ζ_{2Q}]` — `RingContext`, `RingElem` (A.4) |
| `exact_state_ring.py` | fixed-angle exact state over the ring (A.4); cost per reflection depends on `φ(2Q)`, not on depth |
| `ring_cache.py` | persistent on-disk cache for the ring backend |
| `cyclo_heights.py` | exact cylinder heights of a `π/(2Q)` direction class |
| `orphan_region.py` | `find_orphan_ring`, `orphan_width_ring`, `trace_headon`, `locate_orphan` |
| `exact_unfold.py` | exact affine unfolding of an orbit |
| `cyl_diagram.py` | `CylinderDiagram` — genus and stratum of a translation surface |
| `right_triangle_billiards.py` | plain float tracer, for cheap localization only — never for a count |
| `cell_store.py` | the per-branch measurement store |
| `zlattice.py` | exact `ℤ`-lattice arithmetic on integer row vectors |

---

## The orphan paper (`orphan/`)

Built on `lib/` alone (`partition.py` for every count). Each script's docstring lists the
statements it checks.

| Statement | Script |
|---|---|
| §4.3, Rem 5.2, Fig 3, §A.3 — the named counts `n(8/15)`, `n(5/22)`, `n(5/24)`, `n(5/9)`, `n(2/9)`, `n(3/29)`, and the `311`-letter return word of the orphan at `5/24` | `s1_named_counts.py` |
| Theorem A (Thms 4.6, 4.7), Prop 4.5, Cor 5.1 at every centre of the record; §A.3's record tallies (`598` centres, `572` certified, the `26` open, `394` / `245` / `178`) | `s4_theorem_a.py` |
| Theorem 5.3 — the three-sided refinement at the `116` centres `4 ≤ Q ≤ 19` | `s5_3_three_sided.py` |
| Thm 6.1, Cor 6.2, Props 7.1, 7.2, Rems 7.3, 7.4 — `S_α`, its hyperelliptic involution and its base, built from the gluing at every coprime `(P,Q)`, `Q ≤ 200`; §7.5's McMullen correspondence on `min(P, Q−P) = 1`, `Q = 4, …, 15` | `s6_7_surface.py` |
| §A.1, Rem A.4 — the longest passage at `11/41` (`397,153` reflections) and the size of its coefficients | `sA_passage.py` |
| Figures 1–4 | `figures.py` (+ `partition_record.json`) |
| the census record | `export_record.py` (from the research repository's census run) |

⚠ **Scope.** §§2–4 and Appendix B are proofs: Theorems 4.11, 4.18, 4.24, Prop 4.19, Cor 4.20
and Lemmas 4.21–4.23 have no script. Each script checks a statement the paper proves, over the
range it names, as Appendix A.4 says; §A.4's enumeration of §5.2's congruence at `Q ≤ 60` is
elementary and is not repeated. The statements the s740 re-cut moved from this paper to the
necklace paper are checked in `neck/` (below).

---

## The necklace paper (`neck/`)

*The necklace of a rational right triangle: exact counting and complete periodicity.* Built on
`lib/` (`closure.py` walks and `necklace.py` readout for §§3–4 and 7; `partition.py` for branch
partitions; `flow.py`, the exact flow decomposition of Appendix A.2, for every branch count and
complete-periodicity verdict on the band `Q ≤ 60`). Each script's docstring lists the statements
it checks.

| Statement | Script |
|---|---|
| Lemma 3.3, Thms 3.1, 3.4, Cor 3.5, Cor 7.8 and the swapped-edge detector after it (`D = 8 = 2 + 6` at `8/15`), Cor 7.14 at even `Q` — every census class row (odd `Q ≤ 39`, even `Q ≤ 40`) | `s3_s7_identities.py` |
| Thm 4.1, Cor 4.2 — the interval model at every odd `Q ≤ 31` class row | `s4_interval_model.py` |
| Prop 5.3(ii) at the right angle (every `R`-graze boundary a fold, its flanking cells a mirror pair, neither the orphan); no fold at `O` or `A`; Thm 10.3's `A`-graze pattern on the row `P = Q − 2` — every centre with `Q ≤ 30` | `s5_1_folds.py` |
| Rem 9.5 — `R`-graze boundaries as many as mirror pairs, at every centre of the orphan paper's record, `8/15` included | `s9_5_adjacency.py` |
| Thms 10.2, 10.3, 10.5 — the lattice-locus counts at every `Q ≤ 60` | `s10_2_lattice_counts.py` |
| Thms 10.7, 10.11 — the two walks certified by the mechanical prover; their counts and Cor 10.8 (`n(3/Q) = 2 round(Q/3) − 1`, the class completely periodic) at every centre `Q ≤ 60`; the rows `n_H(3/7)`, `n(3/4)`, `n(5/11)`, `n(5/6)` | `s10_windows.py` (+ `window_prover.py`) |
| App A.3 — the printed conditions C1–C12c and D1–D13, root-free and positive; the remarks at `Q = 7` and `Q = 11`; the prover regenerating the printed lists (all but D11) | `sA3_conditions.py` |
| §11 — all `2202` class rows on `2 ≤ Q ≤ 60` decided, `2200` completely periodic; `8/15`'s `6` cylinders and `2` certified minimal components; `⊥H` at the `373` even-`Q` centres (`n_H = Q`, Cor 7.11) and the `212` odd-`Q ≤ 31` controls; the closure certificate's `424` rows; Prop 11.1 and its second witness; Rem 10.10 (`P = 3`, `Q ≤ 60`, and `3/385`); Rem 10.12's `307` even-`P` centres | `s11_cp_band.py` (+ `band_record.json`) |
| Rem 10.12, §13, §5.1 — the named counts `n(8/15)`, `n(5/22)`, `n(5/24)`, `n(17/22)`, `n(19/24)`, `n(2/3)`, by both exact engines | `named_counts.py` |
| Figures 1–4 | `figures.py` |
| Figures 5–6 | `figures_5_6.py` (+ `band_record.json`) |
| the band record | `export_band_record.py` (from the research repository's sweeps) |

⚠ **Scope.** Every other statement of the paper is a proof and has no script. The scripts check
proved statements over the finite ranges they name, and the certificates the paper uses — the
counts of §10, the verdicts of §11, the conditions of A.3 — exactly. They do not discharge complete
periodicity or hypothesis (E) where a statement assumes them; each such statement carries its
scope on its face in the paper.

## Transitive dependencies

These are imported by the legacy scripts and are included so the legacy tree runs, but they
belong to another strand (the counting-law / fan-chain work) and support no statement in any of
the papers:

`s238_fan_anatomy`, `s242_m0fan`, `s252_chain_model`, `s253_proof_ineqs`,
`s254_avoidance_proof`, `s256_leftover`, `s257_warmup_tiling`, `s274_mu_closed`,
`s275_flow_exact`, `s302_model`, `s340_glen_enum`, `s356_rung_exact`, `s357_fanword`.
