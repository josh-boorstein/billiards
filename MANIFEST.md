# Manifest — script to paper statement

Which script supports which statement. Each script's own docstring carries its
pre-registration, its arms and its controls: **the manuscript presents, the docstring is
the source.** Where a docstring and this table disagree, the docstring wins.

The minimal-component and orphan scripts need no setup. The legacy necklace scripts run with `PYTHONPATH=engine:probes` set; some write JSON to `data/`.

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

## The exact-arithmetic core (`engine/`)

These are the four ingredients of the orphan paper's Appendix A, plus the shared
primitives the per-result scripts are built on.

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
| §1.2, §5, Rem 5.4, §9.1, §10 item 2, §A.3 — the named counts `n(5/22)`, `n(5/24)`, `n(17/22)`, `n(19/24)`, `n(8/15)`, `n(5/9)`, `n(2/9)`, `n(2/3)`, `n(3/4)`, `n(3/29)` | `s1_named_counts.py` |
| Theorem D (Thms 4.6, 4.7), Cor 5.1, Prop 4.4 | `s4_theorem_d.py` |
| Props 5.2, 5.3, 9.4 — the closed-form counts, and Prop 5.2's boundary pattern | `s5_closed_forms.py` |
| Theorem 5.5 — the three-sided refinement at the 116 centres `4 ≤ Q ≤ 19` | `s5_5_three_sided.py` |
| Thm 7.1, Props 8.1, 8.2 — `S_α`, its hyperelliptic involution and its base, built from the gluing at every coprime `(P,Q)`, `Q ≤ 200` | `s7_8_surface.py` |
| §9.1 — right-angle boundaries fold, acute ones never do; the `P = Q−2` orphan-flank pattern | `s9_1_folds.py` |
| Remark 9.6 — the `8/15` decomposition (runs the minimal-component scripts) | `s9_4_remark_9_6.py` |
| Figures 1–5 | `figures.py` (+ `partition_record.json`) |
| the census record | `export_record.py` (from the research repository's census run) |

---

## The necklace paper

Taken from the paper's own provenance table.

| § | Script |
|---|---|
| 2 — the necklace, coordinates | `probes/s447_overlap_disk.py`, `probes/s449_necklace.py`, `probes/s450_chain_maps.py` |
| 3 — the counting identity, the covering law | `probes/s444_base_graph.py`, `probes/s446_covering_identity.py` |
| 4 — the interval model, and that it closes | `probes/s451_quotient_path.py`, `probes/s455_closure_proof.py` |
| 5.1 | `probes/s453_pole_module.py`, `probes/s454_composite_q.py`, `engine/zlattice.py` |
| 5.2 | `probes/s454_self_hit.py` |
| 5.4 / 5.5 — the corner analogue | `probes/s454_self_hit.py` (`step`'s terminal `"pole"` branch), `probes/s456_neckgb.py` (`corner_start`), `probes/s492_corner_selfhit.py` (census and interiority arms), `probes/s494_lemma_interiority.py` (Lemma 5.4's arms and its control) |
| 6 — the reduction, and the dead routes | `probes/s456_neckgb.py`, `probes/s462_zz_relation.py`, `probes/s466_triv_gap.py`, `probes/s467_form_lattice.py` |
| 6.3 — the rank law | `probes/s463_rank_law.py` |
| 7 — the Veech locus | `probes/s468_ordering.py`, `probes/s469_profile_proof.py`, `probes/s470_tent_proof.py` |
| 8 — overlap | `probes/s447_overlap_disk.py`; figures `probes/s447_figure.py`; Cor. 8.2's even-`P` half `probes/s472_legswap_existence.py` |
| 9 — complete periodicity | `probes/s457_cp_char.py`, `probes/s459_corner_deep.py`, `probes/s466_cp_transfer.py` |
| **the four figures** — Fig. 1 the kite, and `B` as the necklace; Fig. 2 the `+P` gluing as the star polygon `{Q/P}`; Fig. 3 the chain `D` with its through band, fold and two degeneracies; Fig. 4 the interval model and the quotient path | `probes/s506_necklace_figures.py` |

⚠ `probes/s506_necklace_figures.py` is a figure generator and asserts nothing new: its `verify()`
re-derives every quantity its captions state — the Euler count, both cone angles, the gluing law
`j' = j ± P` from the direction arithmetic alone, the closed form for the interface measures, the
palindrome, the closure of the global coordinate, the nesting, and `|J_i| = |cos(iPπ/Q)|` — and raises
before drawing if any of them fails. Three of those checks are scored against a control that FAILS
(the `+1` and `+2` gluing orders, and a perturbed sector rule), so they are not vacuous. Two of the
four figures are DATA and two are SCHEMATIC, and each says which on its own face.

Also in the strand, cited outside the table: `probes/s464_realisable.py` (realisability),
`probes/s482_boundary_incidence.py` and `probes/s482_regular_leaves.py` (the surface-side
boundary-incidence statements).

⚠ **(G0b) is the strand's one unproved input**, and several statements above are
additionally conditional on complete periodicity. The scripts verify over stated finite
ranges; they do not discharge either hypothesis. Each conditional statement carries its
scope on its face in the paper.

---

## Transitive dependencies

These are imported by the scripts above and are included so the deposit runs, but they
belong to a third strand (the counting-law / fan-chain work) and support no statement in
either paper:

`s238_fan_anatomy`, `s242_m0fan`, `s252_chain_model`, `s253_proof_ineqs`,
`s254_avoidance_proof`, `s256_leftover`, `s257_warmup_tiling`, `s274_mu_closed`,
`s275_flow_exact`, `s302_model`, `s340_glen_enum`, `s356_rung_exact`, `s357_fanword`.
