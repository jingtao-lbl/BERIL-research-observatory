---
reviewer: BERIL Automated Review (Claude, claude-sonnet-4-6)
date: 2026-05-08
project: macro_micro_nutrient_cocycling
---

# Review: Macro-Micro Nutrient Gene Co-cycling in Bacterial Pangenomes

## Summary

This is a well-executed, scientifically rigorous pangenome co-occurrence study grounded in a clearly articulated mechanistic hypothesis (McRose-Newman Fe(III)-oxyhydroxide dissolution). The project successfully tests H1 across 27,682 GTDB species pangenomes using a complete 8-step pipeline with multiple statistical controls (Fisher's exact test, 1,000-permutation null, BH-FDR over 72 pairs, Stouffer meta-analysis), and the notebooks all carry saved outputs — avoiding the major pitfall of code-only submissions. The findings are genuinely interesting: significant P×Metal and N×Metal co-occurrence at genome scale, 100% overlap of phenazine operon carriers with metal-handling genes, and informative negative associations (pstC/pstS depleted vs feoB). The main technical error is a transcription mistake in the REPORT.md results table: the N×Metal Fisher p-value reported there (1.3×10⁻⁸⁷) belongs to the broad N-fixation sensitivity check, not the primary KO-only analysis (actual p=1.3×10⁻⁷¹ per NB02); the Key Findings section correctly states 10⁻⁷¹. This needs to be corrected before the report is shared. Several smaller gaps — a missing environmental-stratification figure, NaN median OR values in the core enrichment table, and the RESEARCH_PLAN's out-of-date species count — are fixable in a short revision pass.

## Methodology

**Research question:** Clearly stated and testable in both README and RESEARCH_PLAN. The explicit H1/H0 pair is a strength. The pivot from fitnessbrowser (documented in the revision history) is honest and appropriate.

**Literature grounding:** McRose & Newman (2021), Edayilam et al. (2018), Wu et al. (2022), and Dakora & Phillips (2002) are cited with sufficient specificity to connect the hypothesis to mechanism. The PhoB-regulon link and Fe(III)-oxyhydroxide surface chemistry rationale are well explained in REPORT.md.

**Approach:** The four-step analysis plan (gene family extraction → co-occurrence statistics → core/accessory enrichment → phylogenetic stratification) is logical and well-scoped. The addition of environmental stratification and a forest plot (Steps 6–8 in the final pipeline) goes beyond the plan and strengthens the submission.

**Gene family definition:** The explicit inclusion/exclusion rationale for the N-fixation PFAM (PF00142 inflates apparent diazotrophic species by 114%) and the phenazine ≥3-gene threshold (with robustness check at ≥2 and ≥4) are highlights of methodological care. The threshold sensitivity table (≥2: 1,125 species; ≥3: 63; ≥4: 34) is exactly the kind of evidence readers need.

**Reproducibility:** The README `## Reproduction` section clearly separates Spark-required steps (1 and 7) from locally reproducible ones (2–6, 8). The `requirements.txt` is present. Intermediate CSV files are committed in `data/`, so downstream notebooks can run without re-executing Spark. This is good practice.

**One gap:** RESEARCH_PLAN.md (v3 revision note) says "27,702 species pangenomes," but the actual analysis runs on 27,682. A brief note in the revision history would reconcile this and prevent reader confusion.

## Code Quality

**Spark query design (src/01):** The PFAM queries are run as separate passes (one JOIN per domain family) to avoid the correlated-subquery problem with EXISTS inside GROUP BY aggregations. This is correct and efficient. The `LIKE 'PFxxxxx.%'` versioned-ID pattern is used throughout, consistent with the pitfalls guide. The gene-name annotation approach for pst and phz genes (where KO assignments are sparse in Bakta) is the right fallback given the database's annotation coverage.

**Statistics (NB02, src/02):** The phi coefficient, Jaccard, Fisher's exact test, and permutation null are all correctly implemented. The BH-FDR correction over 72 nutrient×metal pairs is appropriate. The permutation test correctly shuffles the group membership vector while preserving marginals, and the Z-score and p-value are computed from the null distribution. One minor note: with N=1,000 permutations the minimum representable p is 1/1,001 ≈ 0.001; the code reports this as `perm_p=0.0010` — the REPORT correctly explains this convention.

**Core enrichment (NB03, src/03):** The per-species Fisher's exact test on a 2×2 (core/non-core × nutrient/metal) table followed by Stouffer meta-analysis is the appropriate design for this comparison. The Stouffer Z values (68.3 for P-genes, 3.2 for N-genes, 2.6 for Phz-genes) are clearly interpretable. **However,** the `core_enrichment_summary.csv` and NB03 output show `median_OR=NaN` for all three gene groups. This is almost certainly caused by zero cells in some 2×2 tables (producing infinite per-species ORs that make the median undefined). The NaN column is not discussed in REPORT.md and could confuse a reader inspecting the CSV. Either compute a robust central tendency (e.g., median log-OR excluding infinities) or drop the column and add a note.

**Figure code (NB04, NB07, src/05, src/08):** The multi-panel figure and forest plot code is clean and well-commented. The `TwoSlopeNorm` for the phi heatmap (centering at zero) is the right choice for a signed association metric. The forest plot correctly uses CI half-widths for the error bars, not SDs.

**Notebook organization:** All 7 notebooks follow a consistent setup → load → analyze → display pattern. All have printed outputs or displayed DataFrames — not empty code-only files. NB04 has a large base64 image output (860 KB notebook), which is expected for a multi-panel PNG.

**Auxiliary scripts:** `src/create_notebooks.py` and `src/create_enrichment_notebooks.py` are present in `src/` but are not part of the documented pipeline. They were presumably used to scaffold the notebook files. They do not cause errors but add clutter; consider moving them to a `tools/` or `dev/` subdirectory, or removing them from the committed tree.

**Known pitfall compliance:**
- PFAM versioned IDs (`LIKE 'PFxxxxx.%'`): ✓ used correctly throughout
- Notebook outputs committed: ✓ all notebooks have outputs
- Strain name collision (ENIGMA-specific): not applicable to this project

## Findings Assessment

**Conclusions supported by data:** Yes, with one error to correct. The three headline results (P×Metal phi=0.110, N×Metal OR=10.1, 100% Phz-operon×Metal) are each directly verifiable from NB01 and NB02 outputs. The core/accessory signatures (P-acquisition 70.7% core vs HMA 23.6%) are shown in NB03. The environmental stratification (strongest P×Metal in plant-associated and soil/rhizosphere) is supported by NB06.

**Critical error — N×Metal Fisher p in Results table:** REPORT.md Section 1 (group-level co-occurrence table) reports the N×Metal Fisher p as **1.3×10⁻⁸⁷**. NB02 output shows the actual value for the primary (KO-only) N-fixation definition is **1.259×10⁻⁷¹**. The value 1.3×10⁻⁸⁷ matches the broad N-fixation sensitivity check (Metal_handling × N_fixation_broad, NB02 row 11: fisher_p=1.517×10⁻⁸⁷). The Key Findings section of REPORT.md correctly states "Fisher p=1.3×10⁻⁷¹," so the table value is a transcription error — the sensitivity check p-value was accidentally placed in the primary results row. **This must be corrected.**

**Limitations acknowledged:** Yes, comprehensively. The 8-point limitations section covers presence/absence-only scoring, annotation dependence, PF00142 inflation, phylogenetic non-independence, actinomycete phenazine under-detection, and environmental classification priority. The phylogenetic non-independence caveat is particularly important and correctly flagged — the permutation test controls for marginal gene frequencies but not phylogenetic autocorrelation.

**Incomplete analysis — no figure for environmental stratification:** REPORT.md Section 6 presents a table of log-ORs by environment category but no figure. The forest plot (Figure 2) covers phyla; the environmental breakdown currently has no visual. Given that the environment-stratified result is substantive (strongest P×Metal in plant-associated and soil/rhizosphere, reversal for N×Metal in marine), a figure would improve the report.

**Positive controls:** The 8-species check (NB05) is a clear strength. The P. fluorescens nifH footnote (likely a Fer4-family ferredoxin rather than true nitrogenase) is appropriately flagged. S. coelicolor reclassification to S. anthocyanicus in GTDB R214 is documented. P. chlororaphis correctly appears as the only positive control with a full phenazine operon. All 8 species encode both P-acquisition and metal-handling genes, consistent with the genome-wide signal.

**Xenorhabdus note:** The observation that 6/63 phenazine operon carriers are Xenorhabdus (entomopathogenic nematode symbionts using phenazines for antimicrobial purposes rather than mineral dissolution) is a valuable biological nuance, correctly framed as functional convergence rather than a contradiction of the main hypothesis.

## Suggestions

1. **(Critical) Fix the N×Metal Fisher p-value in the Results table.** Change 1.3×10⁻⁸⁷ to 1.3×10⁻⁷¹ in the group-level co-occurrence table (REPORT.md Section 1). The Key Findings section is already correct. Verify against NB02 output row `N_fixation vs Metal_handling` (fisher_p=1.259×10⁻⁷¹).

2. **(Important) Resolve or annotate the `median_OR=NaN` column in `core_enrichment_summary.csv` and NB03.** Options: (a) replace with the median of finite log-ORs per group, (b) replace with fraction of species showing enrichment (already present as `frac_enriched`), or (c) drop the column entirely and add a code comment explaining that per-species ORs are frequently infinite. As-is, the NaN values will confuse any reader inspecting the CSV.

3. **(Important) Add an environmental stratification figure (NB06 / Figure 3).** A horizontal bar chart of log-OR by environment — with paired bars for P×Metal and N×Metal, error bars, and significance markers — would let readers immediately see the reversal pattern (soil/rhizosphere dominant for P×Metal, marine dominant for N×Metal). The data are in `data/env_cooccurrence.csv`. Saving to `figures/figure3_env_stratification.png` and referencing it in REPORT.md Section 6 would complete this story visually.

4. **(Minor) Update RESEARCH_PLAN.md v3 revision note** to reflect the actual species count (27,682, not 27,702). One line in the revision history suffices: "Species count updated from plan estimate (27,702) to actual query result (27,682)."

5. **(Minor) Move `src/create_notebooks.py` and `src/create_enrichment_notebooks.py`** out of `src/` or document their purpose in the README. These scaffolding scripts are not part of the analysis pipeline and their presence in `src/` implies they are, which can mislead future reproducers.

6. **(Nice-to-have) Add a prominent note in NB01** explaining that the Spark query is pre-cached and cell 3 loads from `data/species_gene_families.csv`. A reader running the notebook for the first time will not realize the output is pre-computed until they inspect the CSV modification timestamps.

7. **(Nice-to-have) Add a case-sensitivity note to the Limitations section.** The Spark queries use `ba.gene = 'pstA'` with specific capitalization. Bakta gene name capitalization is not fully standardized; a brief note that gene-name queries are case-sensitive and may miss differently-cased annotations would close a potential annotation gap not currently covered by the existing Limitation 2.

## Review Metadata

- **Reviewer**: BERIL Automated Review (Claude, claude-sonnet-4-6)
- **Date**: 2026-05-08
- **Scope**: README.md, RESEARCH_PLAN.md, REPORT.md, 7 notebooks (NB01–NB07), 18 data files, 2 figures (figure1_cooccurrence.png/pdf, forest_plot.png), requirements.txt, beril.yaml, src/ (10 scripts), docs/pitfalls.md
- **Note**: This review was generated by an AI system. It should be treated as advisory input, not a definitive assessment.
