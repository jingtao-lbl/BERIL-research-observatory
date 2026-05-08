---
reviewer: BERIL Automated Review (Claude, claude-sonnet-4-6)
date: 2026-05-08
project: macro_micro_nutrient_cocycling
---

# Review: Macro-Micro Nutrient Gene Co-cycling in Bacterial Pangenomes

## Summary

This is a well-executed, scientifically rigorous analysis of macro-micro nutrient gene co-occurrence across 27,682 GTDB bacterial pangenomes. The research question is clearly stated and testable, the statistical toolkit is appropriate and carefully applied (Fisher's exact tests with BH-FDR correction, permutation null, Stouffer meta-analysis), and the findings are significant and well-interpreted against mechanistic and ecological priors. The project is well-structured: seven notebooks with saved text/table outputs, two figures, comprehensive data files, and a thorough REPORT with explicit limitations and future directions. Three issues require attention before this project is fully reproducible by an independent researcher: a factual inconsistency in the forest plot narrative in REPORT.md Section 7, a code–artifact mismatch for the `has_N_fixation` column definition, and an incomplete Reproduction section in README.md that omits scripts 06–08.

---

## Methodology

**Research question:** Clearly stated and testable (H1/H0 in RESEARCH_PLAN.md). The pivot from Substrate C (fitnessbrowser lacks P-starvation data) to a genome-scale pangenome approach is well-motivated and documented in the revision history.

**Approach:** Sound. The four-step pipeline (extraction → co-occurrence → core/accessory → phylogenetic stratification) is appropriate for the question. The choice to define phenazine "operon carriers" at ≥3 distinct phz genes is a reasonable heuristic to distinguish true operon-bearing lineages (63 species) from the broad phzF superfamily (10,856 species), and the threshold is stated explicitly.

**Gene family definitions:** Clearly tabulated in REPORT.md. The use of both KEGG KO and PFAM PF-number anchors (with versioned wildcard matching, e.g. `LIKE 'PF09423.%'`) is correct and avoids the version-suffix pitfall. The sensitivity check for PF00142 (Fer4_NifH domain capturing the broader ferredoxin superfamily) is commendable — this distinction materially changes the N × Metal OR from 10.1 to 4.0 and is properly accounted for.

**Data sources:** Clearly identified (`kbase_ke_pangenome`, GTDB R214, `ncbi_env`). The positive-control species check (NB05) validates the extraction pipeline against eight well-characterized model organisms.

**Reproducibility (data pathway):** Steps 2–8 can be reproduced locally from the cached CSV files in `data/`. Step 1 requires on-cluster Spark access and is clearly flagged in the README. The two-tier pipeline (Spark extraction → local statistics) is a clean and sensible separation. However, the README Reproduction section only lists scripts 01–05 and omits scripts 06–08 (see Suggestions #3).

---

## Code Quality

**SQL and Spark queries (src/01, src/07):** Correct. PFAM versioned IDs are handled with `LIKE 'PFxxxxx.%'` as recommended. The `ncbi_env` join in src/07 handles `NULL` filtering explicitly. The refactored approach of running separate PFAM queries and merging them (used in src/01's second block) avoids issues with large correlated subqueries in Spark.

**Statistical implementation (src/02):** The phi coefficient and Jaccard implementations are correct. The BH-FDR step-up procedure is correctly implemented (rank-ordered, step-down monotone enforcement, clipped to [0,1]). The permutation test uses the both-sided criterion `|phi_null| ≥ |phi_obs|`, which is appropriate.

**Stouffer meta-analysis (src/03):** The per-species Fisher test setup is correct. The 2×2 table `[[nutrient_core, nutrient_noncore], [metal_core, metal_noncore]]` tests whether nutrient and metal genes have *different* core fractions within the same species, which is the right null for the claim "P-acquisition genes are more core than metal-handling genes." The sign convention (z negative when OR < 1) is consistently applied.

**Environmental classifier (src/07):** The keyword priority chain places "root" and "sediment" in `soil_terms` before `marine_terms`, so marine sediment isolates and root-zone endophytes are classified as "soil/rhizosphere" rather than "marine" or "plant-associated." Given the large "other" category (12,628 species), the misclassification risk is modest, but it is worth noting as a caveat on the environmental stratification results.

**Notebook organization:** Logical and consistent. All seven notebooks (NB01–NB07) follow a setup → load cached data → display → summarize structure, and each carries a cross-reference to its equivalent `src/xx_*.py` script. All notebooks contain saved text and table outputs — no empty code-only notebooks.

**Pitfall compliance:** Versioned PFAM IDs (`LIKE 'PFxxxxx.%'`) are used correctly throughout. No strain-name collision issues (analysis operates at GTDB species level). The `kbase_ke_pangenome` tenant prefix is used consistently in all SQL queries.

---

## Findings Assessment

**Main co-occurrence results:** All table values in REPORT.md match the NB02 notebook outputs — phi, OR, Fisher p, and permutation Z are consistent across the REPORT, NB02 text output, and `data/cooccurrence_matrix.csv`. The 100% overlap of phenazine operon carriers with metal-handling genes (63/63) is confirmed by both the NB01 summary and the NB05 positive-control lookup.

**Core/accessory analysis:** The Stouffer Z scores and per-gene core fractions in REPORT.md (pstA 73.7%, feoB 30.0%, HMA 23.6%) match the NB03 outputs precisely. The three-signature interpretation (P-acquisition core-enriched, N-fixation moderately core consistent with HGT, metal-handling bimodal) is biologically plausible and directly supported by the reported numbers.

**Forest plot narrative inconsistency (Important):** REPORT.md Section 7 states "P × Metal: Strongest effects in Cyanobacteriota (log-OR=+1.86, p=4.1×10⁻¹³), Bacillota (+1.64, p=8.0×10⁻²⁴), and Pseudomonadota (+0.94)." However, NB07's output shows Fusobacteriota with log-OR=+2.931 (95% CI: +0.872 to +4.990, p=9.64×10⁻³), which is both larger and statistically significant — larger than the Cyanobacteriota value cited as "strongest" in the REPORT. The REPORT appears to have selected phyla by significance (smallest p-value, reflecting large sample size) rather than by effect size (log-OR), but the wording "strongest effects" implies log-OR magnitude. This discrepancy needs correction.

**Environmental stratification:** The finding that P × Metal is strongest in plant-associated (log-OR=+1.51) and soil/rhizosphere (log-OR=+1.34) environments, while N × Metal peaks in marine (log-OR=+1.37) and animal-associated (log-OR=+2.27), is directly supported by the NB06 output table. The non-significance of animal-associated P × Metal (p=0.39) is correctly reported.

**Undiscussed finding in NB02:** The broad "Phenazine" group (any phz gene ≥1) shows phi=0.262 and OR=5.196 with P-acquisition in the cooccurrence matrix — notably higher than P × Metal (phi=0.110). This is visible in the NB02 output but not discussed in the REPORT. A brief note would clarify whether this reflects genuine co-occurrence of broad phzF-family enzymes with P-acquisition genes, or is driven by the ubiquity of phzF in metabolically active soil organisms regardless of operon completeness.

**Limitations:** The seven limitations in REPORT.md are specific and accurate — particularly the acknowledgment of phylogenetic non-independence and the Actinomycete phenazine under-detection. The PFAM PF00142 sensitivity analysis is an exemplary case of transparent reporting.

---

## Suggestions

1. **[Important] Correct the forest plot narrative in REPORT.md Section 7.** The statement "Strongest effects in Cyanobacteriota (log-OR=+1.86)" is contradicted by NB07, which shows Fusobacteriota at log-OR=+2.931 (p<0.01, CI does not cross zero). Either update the sentence to read "strongest effects among phyla with n≥400 species" (or equivalent qualifier), or add Fusobacteriota to the list of top positive-effect phyla. The Section 4 narrative already correctly identifies Fusobacteriota as having the highest phi (0.453) in the phylum stratification.

2. **[Important] Fix the `has_N_fixation` code–artifact mismatch.** The committed `data/species_gene_families.csv` has `has_N_fixation` = 5,872 species (broad definition, including PF00142), but the current `src/01_extract_gene_families.py` defines `has_N_fixation` using KO-only nifH+nifD (2,746 species). Re-running `src/01` would produce a different CSV than the committed artifact. The note in NB01 documents this ("has_N_fixation in this CSV uses the broad definition; src/02 recomputes KO-only for primary analysis"), and the primary analysis in src/02 is unaffected because it rebuilds N-fixation from individual `has_nifH`/`has_nifD` columns. However, a researcher re-running the pipeline would produce a silently different intermediate file. Recommended fix: update `src/01` to match the CSV (define `has_N_fixation` as the broad group and add a `has_N_fixation_ko` column for the KO-restricted definition), and add a comment in `src/01` explaining the distinction.

3. **[Important] Complete the Reproduction section in README.md.** The current instructions list only scripts 01–05, but the full analysis pipeline includes scripts 06–08. Sections 5 (positive controls), 6 (environmental stratification), and 7 (forest plot) of the REPORT, as well as Figure 2, cannot be reproduced without these steps. Append to the Reproduction section:
   ```bash
   # Step 5: Positive-control species check (runs locally from cached data)
   python src/06_positive_controls.py

   # Step 6: Environmental stratification (requires Spark for ncbi_env query)
   python src/07_environmental_stratification.py

   # Step 7: Forest plot (runs locally from cached data)
   python src/08_forest_plot.py
   ```

4. **[Nice-to-have] Add `berdl_notebook_utils` to requirements.txt.** The extraction scripts (src/01, src/07) import `from berdl_notebook_utils.setup_spark_session import get_spark_session`, which is not listed in `requirements.txt`. A researcher following the reproduction guide will hit an import error on the Spark-dependent steps without it. Consider adding it with an inline note, e.g.:
   ```
   berdl_notebook_utils  # Required for Spark steps (01, 07); not needed for local steps (02–06, 08)
   ```

5. **[Nice-to-have] Briefly discuss the broad Phenazine group phi in the REPORT.** The NB02 cooccurrence matrix shows P × broad-Phenazine (any phz ≥1, n=10,308) with phi=0.262 and OR=5.196 — substantially higher than P × Metal (phi=0.110). This is visible to any reader who examines the cooccurrence matrix output but is not addressed in the REPORT. A one-sentence note in Methods or Limitations ("the broad phzF superfamily shows phi=0.26 with P-acquisition; we do not interpret this as mechanistic evidence because phzF is prevalent across phylogenetically diverse organisms regardless of operon completeness") would preempt reviewer questions.

6. **[Nice-to-have] Note environmental classifier ordering as a caveat.** In `src/07`, terms like "root", "sediment", and "mud" appear in `soil_terms` and are matched before the marine keyword check, so marine sediment isolates are routed to "soil/rhizosphere." A single sentence in the Limitations section (e.g., "Environment classification relies on keyword matching with a fixed priority order; marine sediment and root-zone isolates may be misclassified") would complete the methodological transparency already shown for the other limitations.

---

## Review Metadata
- **Reviewer**: BERIL Automated Review (Claude, claude-sonnet-4-6)
- **Date**: 2026-05-08
- **Scope**: README.md, RESEARCH_PLAN.md, REPORT.md, 7 notebooks (NB01–NB07), 8 source scripts (src/01–08), 17 data files, 3 figures (figure1_cooccurrence.png/.pdf, forest_plot.png), requirements.txt, beril.yaml
- **Note**: This review was generated by an AI system. It should be treated as advisory input, not a definitive assessment.
