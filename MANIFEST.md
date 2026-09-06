# Manifest — script to paper statement

Which script supports which statement. Each script's own docstring carries its
pre-registration, its arms and its controls: **the manuscript presents, the docstring is
the source.** Where a docstring and this table disagree, the docstring wins.

Run anything here with `PYTHONPATH=engine:probes` set. Some write JSON to `data/`.

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

## The orphan paper

Mapped from Appendix A and from each script's own pre-registration. The counts of §5 and
§9.3 come from the core above via `n_pq_certified`; the scripts below are the per-result
checks around them.

| Statement | Script |
|---|---|
| Theorem D — existence, uniqueness, parity, turnaround type of the orphan | `probes/orphan_theorem.py` |
| the vertex gap | `probes/s332_vertexfree.py` |
| word normalisation against the tracer's full word | `probes/s373_word_normalise.py` |
| graze structure | `probes/s374_graze_structure.py` |
| continuity of the first-return involution `T` at a cell boundary | `probes/s374_tcont.py` |
| does the limit orbit at a `J`-fold boundary reverse at the midpoint | `probes/s375_foldmid.py` |
| the residual audit behind the converse | `probes/s377_residual_audit.py` |
| §5, §9.3 — exact cylinder words and counts at large `Q` | `probes/s196_ringwords.py` |
| cylinder count and area reconciliation | `probes/s315_cylinder_count.py` |
| the odd-`P` separatrix census | `probes/s315_oddP_census.py`, `probes/s276_separatrix_census.py` |
| §7 — the unfolded surface `S_α`, identification and verification | `probes/veech_stratum.py` |
| §7 — the Apisa import, reconciled with our surface | `probes/veech_apisa.py` |
| §7, Thm 7.1 — **the even-`P` stratum**; the odd-`P` table does not apply there | `probes/s500_evenP_stratum.py` |
| §9 — cylinders of both direction classes on the Veech locus | `probes/s437_oblique_cylinders.py` |
| §9 — the separatrix / complete-periodicity census | `probes/s438_oddclass_cp.py` |
| §9 — a strip's cells of constant word, exactly | `probes/s439_exact_cells.py` |

⚠ `probes/s500_evenP_stratum.py` reports `488` rows at `Q ≤ 40`. Where the paper quotes
the Weierstrass census it quotes a far stronger independent result (`12230/12230` coprime
`(P,Q)` with `Q ≤ 200`), not this probe's coverage.

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
