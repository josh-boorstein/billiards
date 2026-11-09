# Manifest — script to paper statement

Which script supports which statement. Each script's own docstring carries its
pre-registration, its arms and its controls: **the manuscript presents, the docstring is
the source.** Where a docstring and this table disagree, the docstring wins.

The scripts need no setup: each puts `lib/` on the path itself.

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
periodicity or hypothesis (OB) where a statement assumes them; each such statement carries its
scope on its face in the paper.

---

## The method paper (`method/`)

*No arrivals: complete periodicity and degenerate orbits in rational right triangles.* Built on
`lib/` alone (`triangle.py`, `tracer.py` for §6; `mpmath`, `sympy`, `scipy` for the appendices and
Theorem 6). One script per computation of the paper's *Code availability*. Each script's
docstring lists the statements it checks.

| Statement | Script |
|---|---|
| Prop 6.4, its setting — the directions met, the phase space and its area, the field `K` | `s6_4_field.py` |
| Prop 6.4, *The return map* — the first return to `Σ_H̄` is a three-interval exchange: words, exact translations, `λ₁, λ₂, λ₃` and coordinates | `s6_4_return_map.py` |
| Prop 6.4, *The return map* — seven domains tile `[0, 1]` exactly, so there is no eighth word | `s6_4_tiling.py` |
| Prop 6.4, *Minimality*, *Unique ergodicity* — rank `3`, the circle rotation, the roof | `s6_4_minimality.py` |
| Prop 6.4, *Nothing else* — `L1` meets only the three cylinders | `s6_4_cylinders.py` |
| Prop 6.4, *Nothing else* — the cylinder traces and areas, the Kac mass, their sum `Q = 15` | `s6_4_mass.py` |
| Prop 6.4, *Code availability* — the seven return words and the six cylinder half-words | `s6_4_words.py` |
| §6.3 — the connection `O → A` on the boundary of the minimal component | `s6_3_boundary.py` |
| §6.3 — the arrival at `T₀ = 8/15`, the swapped connection it is, the reversed orbit at `7/15` | `s6_3_arrival.py` |
| Theorem 6 — on each of `I₁`, `I₂`, `I₃`: the words and their exact tiling, every constraint certified on the interval, the symmetric three-interval exchange, the Galois step (`φ(d) ≤ 20 ⟹ d ≤ 66`, least denominators `92`, `242`, `223`), the periodic orbits of `282`, `118`, `126` bounces | `s9_2_theorem6.py` |
| Prop A.12 — the certificate, `1 672` cases | `sA_12_certificate.py` |
| Lemma A.15 — the finite check at the eight 3-smooth `80 < Q < 280`, `258` good classes | `sA_15_finite_check.py` |
| Prop A.23 — the census at `Q ≤ 80`, `95 394` cases | `sA_23_census.py` (+ `sA_23_record.json`) |
| §B.4 — the fourteen inequalities of Theorem 7.11's proof at every `m ≥ 3`, and the `R`-prongs at `m = 1, 2` | `sB_4_certificate.py` (+ `sB_4_record.json`) |

⚠ **Scope.** Theorems 2–5, Proposition 9.1 and the rest of the paper are proofs and have no
script. Remarks 6.3 and 6.6 report censuses the paper calls *found by computation*; they are not
deposited. Theorem 6's words are proposed by a floating-point tracer at one point of each
interval; what is certified is that those words tile the transversal on the whole interval.

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
| §6.4 — Remark 6.10: `cos 2a = (√7 − 1)/4` lies in `π/6 < a < π/4`, `a/π ∉ ℚ`; the period-`944` primitive non-perpendicular cylinder in `⊥OA`, exactly; controls (no closure at the neighbouring rational `−2/5` in the word's frame, no one-letter change certifies) | `s6_4_troubetzkoy.py` (+ `tro05_record.json`) |
| Figures 1–4 | `figures.py` |

⚠ **Scope.** A certificate proves that the cylinder it names exists. Table 1's blank cells are
searches that found nothing, not proofs of absence, and nothing here checks them. §§2–4,
§§6.1–6.2, §6.4 apart from Remark 6.10's instance, and §7 are proofs in the paper and have no script.

---

---

## The polygon note (`polygon/`)

*Singular perpendicular orbits and complete periodicity in rational triangles.* Built on `lib/`
(`flow.py`, `beams.py`, `slots.py`; `closure.py` and `necklace.py` for the halves of isosceles
triangles at odd `N`; `sympy` for Proposition 8.2). One script per item of §9, and one for
Proposition 8.2's computer-assisted step.

| Statement | Script |
|---|---|
| Theorem 4.2 (with Lemma 4.1), Proposition 5.1, Theorem 5.2 — on every class of the `456` scalene triangles with odd `N ≤ 21` or even `N ≤ 34` (`820` classes); the control `(7, 8, 15)`: `s′ = 6`, `p = 1`, `Σ n_σ = 6` | `s4_s5_slots_and_bands.py` |
| §9, odd `N` — `1578` scalene triangles with odd `N ≤ 51`: `p ≥ 1` on `18`, all among the `626` with some `gcd(a_i, N) ≥ 3`, all with `s′ = 0`; the four with `N ≤ 21` and `Σ n_σ = N + 1`; `p = 0` at prime `N` | `s9_odd_census.py` |
| Theorem 1.4, Propositions 8.1–8.2, Remark 8.3 — the families to `N ≤ 303`: completely periodic, `p = 1`, resp. `p = m` | `s8_families.py` |
| Proposition 8.2 — the containment step: the direction identity, the harmonic form of the side functions, the exact root counts and the fourteen sign conditions | `s8_2_containment.py` |
| Theorem 6.2 (1)–(2) — the leg and base beams of the `47` isosceles triangles with even `N ≤ 30` against their halves | `s6_isosceles_beams.py` |
| Remark 8.5 — the `158` isosceles triangles with odd `N ≤ 39`: completely periodic with `p = 0` (`157` by the half's certificate, all `158` by the triangle's own decomposition) | `s8_5_isosceles_odd.py` |
| Proposition 8.6 and §9 — the `91` isosceles triangles with even `N ≤ 60`, `b < 2a`: both classes completely periodic; the controls `(7, 7, 16)`, `(4, 4, 7)` | `s8_6_isosceles_even.py` |
| Figures 1–3 | `figures.py` |

⚠ **Scope.** Everything else in the note is proved by hand and has no script, except that the case of
Proposition 8.6 in which `N` has one odd prime factor imports the method paper's Corollary 8.1, part 1,
whose computer-assisted steps (its Appendix A.4–A.6) are in `method/` (`sA_12_certificate.py`,
`sA_15_finite_check.py`, `sA_23_census.py`). The classifier's
branch positions are 40-digit evaluations (counts and kinds are exact); the figures' paths are
floating-point drawings whose content is asserted against the exact classifier first.
