---
reviewer: BERIL Automated Review (Claude, claude-sonnet-4-6)
date: 2026-05-08
project: macro_micro_nutrient_cocycling
---

# Review: Macro-Micro Nutrient Gene Co-cycling in Bacterial Pangenomes

## Summary

This is a well-executed and scientifically rigorous pangenome co-occurrence study. The research question is clearly framed, the hypothesis is falsifiable, and the analysis addresses it systematically across 27,682 GTDB species pangenomes with appropriate statistical methods (Fisher's exact test, phi coefficient, permutation null, BH FDR correction, Stouffer meta-analysis). The project goes well beyond the original 4-step plan, adding positive controls, environmental stratification, and a per-phylum forest plot. All seven notebooks have saved outputs, three publication-quality figures are present, and the REPORT.md contains a thorough treatment of results including a substantive limitations section. One important factual error appears in the README Overview section — the N×Metal statistics reported there are the broad/sensitivity-check values rather than the primary KO-only values — and should be corrected before any public sharing. Several minor display and documentation gaps round out the suggestions.

## Methodology

**Research question and hypothesis:** Both are stated with precision. RESEARCH_PLAN.md distinguishes H1 (positive association exceeding permutation null) from H0 (no greater than chance) clearly. The approach — presence/absence co-occurrence across species-level pangenomes — is appropriate for a genome-scale survey of this kind.

**Data sources:** Primary data source (`kbase_ke_pangenome`, 132.5M gene clusters, GTDB R214) is identified unambiguously. The gene family definitions table in REPORT.md specifies annotation source (KO, PFAM, gene name) for each of the 23 gene families. The distinction between KO-only N-fixation (primary) and PF00142-inclusive N-fixation (sensitivity check) is clearly explained and its impact quantified (2,746 vs. 5,872 species, OR dropping from 10.1 to 4.0).

**Statistical approach:** Multiple complementary measures are used (Jaccard, phi, OR, Fisher p, permutation Z). BH FDR correction is applied to all 72 individual-gene-pair Fisher tests (18 nutrient × 4 metal). The permutation null preserves marginal counts, which is appropriate for controlling for prevalence differences. The core/accessory analysis uses per-species Fisher tests aggregated via Stouffer meta-analysis — a defensible approach for a heterogeneous multi-species dataset.

**Phylogenetic non-independence:** Correctly flagged as Limitation 4 and as a future direction (phylogenetic independent contrasts). No attempt is made to control for it in the current analysis. Given the dataset spans all bacterial phyla, this is a genuine caveat for the group-level phi values. The phylum-level stratification and positive-control checks provide useful partial mitigation.

**Phenazine operon definition:** The ≥3-gene threshold is justified and robustness-tested (≥2: 1,125 species; ≥4: 34 species). The trade-off between sensitivity and specificity is appropriately discussed, including the known under-detection of actinomycete phenazine pathways.

**Environmental stratification:** Applying BH FDR correction across all 21 Fisher tests (7 environments × 3 pairs) is methodologically sound. The large "other" category (12,628 species, ~47% of the environmentally classified subset) is acknowledged in Limitation 8 but not quantified in the main environmental results section, which could give readers a misleading sense of environmental coverage.

## Code Quality

**SQL queries (`src/01_extract_gene_families.py`):** Queries are well-structured. PFAM ID matching uses `LIKE 'PFxxxxx.%'` to handle versioned IDs, consistent with BERDL best practice. The two-phase approach — running a main KO/gene-name query then separate PFAM queries and merging — is a reasonable workaround for the EXISTS-in-GROUP-BY limitation in Spark SQL. Subquery is skipped cleanly (`if 'EXISTS' in condition: continue`) with PFAM families handled separately.

**Gene name case sensitivity:** `ba.gene = 'pstA'` and similar exact-match predicates are used throughout. This is noted as Limitation 9. The potential for missed annotations (e.g., `PstA`, `PSTA`) is real for gene-name-based families but does not affect KO- or PFAM-based families, which cover the more novel/confident identifications. No attempt is made to estimate the magnitude of missed annotations.

**BH FDR implementation (`src/02_cooccurrence_stats.py`):** The manual BH step-up implementation is correct: sorted p-values are multiplied by `n_tests/rank`, monotonicity is enforced with a backward pass, and values are clipped to [0,1]. This matches scipy's BH behavior and is standard.

**Permutation test:** Uses a two-tailed comparison (`np.abs(null_phis) >= np.abs(observed_phi)`) with Laplace correction (`(count + 1) / (N_PERM + 1)`). Seed is fixed (`np.random.seed(42)`) for reproducibility. With 1,000 permutations, the minimum representable p is 1/1001 ≈ 0.001, clearly stated in the report.

**PFAM merge logic (minor):** In `src/01_extract_gene_families.py`, the PFAM results are merged into the main DataFrame using `pdf.loc[pfam_df.index, col] = pfam_df[col]`. If any PFAM query returns species IDs absent from the main query, `pdf.loc` would silently create new rows. This is a low-probability edge case in practice, and the subsequent `fillna(0)` would handle NaN propagation. Not critical but worth a defensive assertion.

**pitfalls.md relevance:** The pitfall on "Commit Notebooks Alongside Their Artifacts" applies here. In this project, all seven NB01–NB07 notebooks exist as committed files with outputs — the pitfall is avoided. The RESEARCH_PLAN.md clarifies that `src/*.py` scripts are the primary analysis code and notebooks are display layers, which is an unusual but valid and well-documented architecture.

**Notebooks:** All seven notebooks (NB01–NB07) have saved text outputs. NB04 generates and saves figures. NB06 shows environmental stratification but does not compute or display BH q-values (these are computed in `src/07` and reported in REPORT.md). Adding a q-value display cell to NB06 would close this gap.

## Findings Assessment

**Core results:** The three primary co-occurrence findings (P×Metal phi=0.110/OR=2.3, N×Metal phi=0.088/OR=10.1, Phz-operon×Metal 63/63) are reported consistently in REPORT.md and confirmed by NB01/NB02 outputs. The interpretation — that coupling is genomically encoded but mechanistically distinct for P-acquisition (core), N-fixation (HGT-mobile), and phenazine (clade-specific) — is well-supported by the core/accessory analysis.

**Negative associations:** The pstC/S–feoB anti-correlation (phi as low as −0.256) is an interesting finding. The report interprets it as ecological niche separation between aerobic P-scavengers and anaerobic Fe(II)-uptake specialists. This interpretation is coherent but appropriately framed as speculative.

**Statistical significance vs. effect size:** The project correctly notes that phi values are modest (0.01–0.11) despite high significance. The Phz-operon × Metal result (phi=0.014, p=9.4×10⁻³, permutation Z=2.3) is the weakest of the three primary results, and the report does a good job contextualizing it against the 91.7% background metal-gene prevalence.

**Limitations:** The nine limitations are substantive and specific. Particularly commendable: Limitation 3 (PF00142 superfamily inflation with quantitative impact), Limitation 7 (actinomycete phenazine under-detection with positive control example), and Limitation 9 (gene name case sensitivity).

**Important discrepancy — README Overview vs. REPORT.md for N×Metal:** The README Overview section states:
> "N-fixation and metal-handling genes (phi=0.107, OR=4.0, p=1.5×10⁻⁸⁷)"

These values correspond to the **N_fixation_broad** (PF00142-inclusive sensitivity check) result, not the primary KO-only result. From NB02 output and REPORT Table 1, the correct primary values are **phi=0.088, OR=10.1, p=1.3×10⁻⁷¹**. Specifically: 0.107 is the Jaccard index for N×Metal (not phi), and OR=4.0 and p=1.5×10⁻⁸⁷ are from the N_fixation_broad × Metal_handling row. The REPORT.md Table 1 and Key Findings section have the correct values throughout. The README needs to be corrected to avoid misrepresenting the primary result.

**Minor display gap in NB03:** The gene family core fraction display (NB03 cell 4) omits `pstB` and `phnD` from the `all_genes` list. These genes are included in the analytical code (`src/02`), just not displayed. Not an analytical error, but a gap in the inspection notebook.

## Suggestions

1. **[Important] Fix N×Metal statistics in README Overview.** Replace `phi=0.107, OR=4.0, p=1.5×10⁻⁸⁷` with the correct primary KO-only values: `phi=0.088, OR=10.1, p=1.3×10⁻⁷¹`. The current README cites the sensitivity-check broad values as the headline result, which could mislead readers who do not read REPORT.md.

2. **[Important] Add BH q-values to NB06 output.** NB06 currently shows raw Fisher p-values, but REPORT.md cites BH-corrected q-values for environmental stratification. Add a cell to NB06 that computes and displays q-values (or loads the FDR-corrected results if `src/07` saves them separately). This makes NB06 self-consistent with the report.

3. **[Nice-to-have] Add pstB and phnD to NB03 core fraction display.** The `all_genes` list in NB03 cell 4 omits pstB and phnD, which are included in the analysis. Adding them closes a minor documentation gap and ensures the notebook faithfully represents all 23 gene families.

4. **[Nice-to-have] Quantify the "other" environment category in the main environmental results section.** 12,628 species (~47% of environmentally classified species) fall into "other/unclassified." The REPORT cites the "other" log-OR results but does not flag the large size of this category inline. A brief parenthetical in the environmental stratification paragraph would help readers interpret the scope of the classification.

5. **[Nice-to-have] Explain forest plot p=1.0 entries.** NB07 output shows several phyla with positive log-OR but p=1.00 (e.g., Deinococcota log-OR=+4.65 [−0.155, +9.463], p=1.00). This is a known Fisher's exact test behavior for degenerate small-n contingency tables. A brief note in REPORT.md Section 7 or NB07 explaining why high log-OR can coexist with p=1.0 in small phyla would prevent confusion for readers.

6. **[Nice-to-have] Add a spot-check for gene name case sensitivity.** Limitation 9 acknowledges exact case matching for `pstA`, `pstB`, etc. A brief query in NB01 counting hits for alternative casings (e.g., `ba.gene IN ('PstA', 'PSTA')`) would either confirm the concern is negligible or quantify the potential miss rate.

7. **[Nice-to-have] Guard against PFAM merge creating extra rows in `src/01`.** Add `assert len(pdf) == initial_len` (captured before PFAM merges) after each PFAM update loop to catch any unexpected row inflation from species present in PFAM results but absent from the main query.

## Review Metadata
- **Reviewer**: BERIL Automated Review (Claude, claude-sonnet-4-6)
- **Date**: 2026-05-08
- **Scope**: README.md, RESEARCH_PLAN.md, REPORT.md, 7 notebooks (NB01–NB07), 2 source files (src/01, src/02), requirements.txt, beril.yaml, 3 figures, 12+ data files, docs/pitfalls.md
- **Note**: This review was generated by an AI system. It should be treated as advisory input, not a definitive assessment.
