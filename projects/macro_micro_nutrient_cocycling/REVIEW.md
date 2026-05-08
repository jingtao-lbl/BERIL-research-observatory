---
reviewer: BERIL Automated Review (Claude, claude-sonnet-4-6)
date: 2026-05-08
project: macro_micro_nutrient_cocycling
---

# Review: Macro-Micro Nutrient Gene Co-cycling in Bacterial Pangenomes

## Summary

This project delivers a well-executed, genome-scale co-occurrence analysis across 27,682 GTDB species pangenomes, testing whether P-acquisition, N-fixation, and phenazine biosynthesis gene families co-occur with trace-metal-handling genes more than expected by chance. The research question is sharply defined and mechanistically grounded in the McRose-Newman Fe(III)-oxyhydroxide dissolution model. The statistical toolkit is appropriate and layered (Fisher's exact, phi coefficient, Jaccard, 1,000-permutation null, BH-FDR over 72 gene pairs, Stouffer meta-analysis, per-phylum forest plot, environmental stratification), and a thorough 9-point limitations section and positive-control species check demonstrate unusual analytical rigor. Reproducibility is well-handled: all seven notebooks have saved outputs, the Spark/local split is clearly documented in the README, a `requirements.txt` is present, and the primary analysis lives in numbered `src/*.py` scripts with the notebooks serving as display-only companions. The main gaps are a single interpretive inconsistency in the core-enrichment section (Stouffer Z sign vs. median OR for phenazine genes), an undiscussed negative phzF × feoB association that is mechanistically relevant, and the FDR-adjusted q-values for the environmental stratification not being written back to the CSV artifact. Phylogenetic non-independence remains uncontrolled, as the authors acknowledge, and is the central methodological limitation for future work.

## Methodology

**Research question and hypothesis:** Clearly stated in both README and RESEARCH_PLAN.md. H1 (positive co-occurrence exceeding a permutation null) and H0 are operationalized in testable, computable terms. The four-step plan (gene family extraction → co-occurrence → core/accessory enrichment → phylogenetic stratification) is logical and well-scoped for the available data.

**Approach:** Using species-level presence/absence across 27,682 GTDB pangenomes is an appropriate grain for this question. The dual-definition strategy for N-fixation (KO-only primary, PFAM PF00142 as sensitivity check) is thoughtfully justified — the PFAM domain captures non-diazotrophic ferredoxins and inflates the species count by 114%, so excluding it from primary analysis is correct. The ≥3 phz-gene threshold for "operon carrier" is documented with a threshold sensitivity table (≥2: 1,125 species; ≥3: 63; ≥4: 34), which supports the choice.

**Data sources:** `kbase_ke_pangenome` (132.5M gene clusters, Bakta annotations, PFAM domains, GTDB R214 taxonomy) is correctly identified as the primary source. The environmental stratification adds `ncbi_env` metadata, which is also identified and queried appropriately. The RESEARCH_PLAN.md documents the pivot from `kescience_fitnessbrowser` (no P-starvation experiments available) to pangenome presence/absence, which is a sound design decision with appropriate provenance.

**Reproducibility:** The README Reproduction section clearly distinguishes Spark steps (01, 07) from locally runnable steps (02–06, 08). The note that re-running a notebook does *not* re-query BERDL is explicit and correct. A `requirements.txt` covering numpy, pandas, scipy, matplotlib, pyspark, and berdl_notebook_utils is present. The RESEARCH_PLAN.md revision history (v1–v4) provides a useful audit trail of design decisions. One gap: the FDR-adjusted q-values computed in NB06 are not written back to `data/env_cooccurrence.csv`, meaning the CSV artifact alone does not fully reproduce the FDR results reported in REPORT.md (they must be re-derived from the inline NB06 cell).

## Code Quality

**SQL / Spark queries (src/01, src/07):** PFAM domain queries correctly use versioned `LIKE 'PF#####.%'` matching, which is the documented best practice for this schema. Separating EXISTS-based PFAM conditions into secondary queries (rather than embedding them in the main GROUP BY) is a reasonable efficiency choice given Spark's query planner behavior on large correlated subqueries. The `ncbi_env` join in src/07 filters out sentinel null values ('missing', 'not collected', etc.) appropriately.

**Statistical methods (src/02, src/03):** The BH-FDR implementation in src/02 is correct: ascending sort by p-value, q_raw = p × n / rank, then a single backward pass to enforce monotonicity. The permutation null preserves marginal species counts. The Stouffer meta-analysis aggregates per-species Fisher z-scores, which is an appropriate approach for synthesizing small within-species effects across heterogeneous species.

**Potential interpretive inconsistency — Phz_genes Stouffer Z vs. median OR (NB03, REPORT.md):** The core-enrichment summary shows that phenazine genes have a positive Stouffer Z = 2.636 (p = 0.008), while simultaneously showing median per-species OR = 0.333 and frac_enriched = 0.268 (only 26.8% of co-encoding species have OR > 1). A positive Stouffer Z conventionally means nutrient genes are more core than metal genes in the aggregate, yet median OR < 1 indicates the opposite in the typical individual species. REPORT.md states phenazine genes are "not preferentially core relative to metal genes," which aligns with median OR and frac_enriched but conflicts with the positive Stouffer Z. The most likely explanation is that a small number of species with large positive per-species z-scores (e.g., *Streptomyces* where phzF is core but HMA is accessory) dominate the Stouffer aggregate — a known vulnerability of unweighted Stouffer meta-analysis to outliers. This tension is neither acknowledged nor explained in the text, and the phrasing "not preferentially core" requires qualification.

**Undiscussed negative association — phzF × feoB (NB02, REPORT.md):** Of 72 nutrient × metal gene pairs, phzF × feoB is the fourth most negative (phi = −0.163, enrichment = 0.72×, FDR q < 10⁻¹⁶⁰). REPORT.md's "Strongest negative associations" table lists only pstC/S/A × feoB/HMA; phzF × feoB is absent from the discussion. This is notable because phzF is a phenazine gene and the mechanism under study involves phenazine-mediated Fe(III)-oxyhydroxide dissolution — an anti-correlation between phzF and direct Fe(II) uptake (FeoB) would be mechanistically interpretable and strengthens the ecological model. The omission leaves a potentially informative signal unaddressed.

**Environmental classifier (src/07):** The keyword priority order places 'sediment' and 'mud' in the `soil_terms` list, ahead of marine keywords. As documented in Limitation #8, this will misclassify marine sediment isolates as soil/rhizosphere. Given that the soil/rhizosphere group (n=3,406) and marine group (n=3,449) are similar in size, misclassification at this boundary could modestly attenuate the contrast between these environments. The magnitude is likely small given the large "other" category (n=12,624) absorbs most unclassifiable entries.

**Code organization:** The `src/*.py` scripts are logically numbered, documented with module-level docstrings, and produce named CSV outputs. The notebooks correctly load only pre-cached CSVs and do not duplicate the primary analysis logic. The correspondence between NB01–NB07 and src/01–src/08 is explicit in each notebook's markdown header.

**Pitfall compliance (docs/pitfalls.md):** PFAM versioned matching (`LIKE 'PF#####.%'`) follows the documented best practice. GTDB R214 reclassification of *S. coelicolor* to *S. anthocyanicus* is documented with the genome accession (GCF_008931305.1) in both the positive controls section and Limitation #7. Gene name case sensitivity (Limitation #9) correctly identifies that `ba.gene = 'pstA'` could miss `PstA` or `PSTA` variants, which is the relevant gotcha for name-based gene queries in Bakta. No pitfalls from docs/pitfalls.md apply directly to this project's data sources (it does not use ENIGMA strain matching, Web of Microbes action codes, cMD cohorts, or Fitness Browser KO mapping).

## Findings Assessment

**Core claims are well-supported:**
- P × Metal (phi=0.110, OR=2.3, Fisher p=1.3×10⁻⁶⁵, permutation Z=17.7): Robustly supported by both parametric and permutation tests.
- N × Metal (phi=0.088, OR=10.1, Fisher p=1.3×10⁻⁷¹, permutation Z=14.7): Strong signal; the high OR is correctly explained as near-universal metal gene prevalence in KO-defined diazotrophs (2,719/2,746 = 99.0%).
- Phz-operon × Metal (63/63, OR=∞, Fisher p=9.4×10⁻³, permutation Z=2.3): The authors correctly contextualize this: the signal is modest against a 91.7% background metal gene prevalence. The back-of-envelope calculation `0.917^63 ≈ 0.004` for the probability of 63/63 by chance is accurate and helpfully provided.
- Individual gene pairs: 64/72 FDR-significant at q<0.05. The strongest positive associations (phzG×corA 1.79×; phoD×HMA 1.72×; nifH×HMA 1.42×) are biologically interpretable and accurately reported from the pairwise_detail.csv.

**Numbers are internally consistent:** Species counts (27,682 total; 27,017 with environment metadata; 27,009 with complete annotations for environmental stratification) are consistent across README, NB01 outputs, NB06 outputs, and REPORT.md. The per-family phenazine carrier counts (24 Streptomycetaceae, 17 Pseudomonadaceae, 10 Streptosporangiaceae, 6 Enterobacteriaceae) are consistent between NB03 cell output and REPORT.md text.

**Environmental stratification (Result 6):** P×Metal significant in plant-associated and soil/rhizosphere, N×Metal strongest in marine and animal-associated — the pattern is internally consistent and biologically plausible. BH-FDR is correctly applied across all 21 tests. The q-values quoted in the REPORT table match those computed in NB06 (verified against the notebook output).

**Positive controls:** All eight model organisms encode both P-acquisition and metal-handling genes. *P. chlororaphis* (phzA/B/F/G) is the only positive control with a ≥3-gene phenazine operon, consistent with its known phenotype. The *P. fluorescens* nifH footnote (likely a divergent Fer4-family ferredoxin) is appropriately caveated and noted not to affect primary analysis.

**Limitations section:** The 9-point limitations section is one of the strongest aspects of this report. Points 1 (presence/absence only), 3 (PF00142 Fer4 superfamily inflation), 4 (phylogenetic non-independence), 6 (operon threshold heuristic), and 7 (actinomycete phenazine under-detection) directly address the most substantive weaknesses. No major limitations appear to be omitted.

**Incomplete analysis:** No cells are left with "to be filled" placeholders. The Forest Plot (NB07) and Environmental Stratification (NB06) represent genuine findings beyond the initial four-step plan, added as documented enrichments in v4.

## Suggestions

1. **(Critical) Resolve the Stouffer Z / median OR tension for phenazine core enrichment.** For the Phz_genes group (Stouffer Z=+2.636, median OR=0.333, frac_enriched=0.268), explain the apparent contradiction between a positive aggregate Z and median per-species OR < 1. Verify in src/03 whether the result is driven by outlier species (e.g., *Streptomyces* species where phzF is core but HMA/feoB are predominantly accessory) inflating the Stouffer sum. If outlier-dominated, consider winsorizing per-species z-scores or reporting both the Stouffer Z and the median OR with an explanation of why they point in different directions. As written, the phrase "not preferentially core relative to metal genes" is potentially inconsistent with the positive Stouffer Z and will confuse statistically sophisticated readers.

2. **(Critical) Discuss the phzF × feoB negative association (phi=−0.163, q<10⁻¹⁶⁰).** This is the fourth strongest negative pair in the entire 72-pair analysis and directly involves a phenazine gene — the centerpiece of the study's mechanistic model. One to two sentences would suffice: organisms using phzF-mediated phenazine biosynthesis to reductively dissolve Fe(III)-oxyhydroxides under aerobic/suboxic conditions may not co-encode direct Fe(II) uptake via FeoB, because phenazines mobilize P and metals from mineral surfaces where Fe is in the oxidized form. This negative association would reinforce rather than undermine the McRose-Newman model and should not be left out of Results 2.

3. **(Moderate) Write BH-corrected q-values to `data/env_cooccurrence.csv`.** The FDR adjustment is currently computed inline in NB06 and not saved. Add a `q_value` column to the CSV output in src/07 (or as a post-processing step in src/08 or a new script) so that the CSV artifact is self-contained and the FDR significance calls can be reproduced without re-running the notebook.

4. **(Moderate) Add a visual indicator for the nifH(Pfam) sensitivity row in Figure 1, Panel A.** The figure includes a row for `nifH(Pfam)` but its status as a sensitivity check excluded from the 72-pair FDR correction is explained only in the REPORT caption. A visual cue — dashed row border, lighter background, or "(sensitivity only)" in the y-axis label — would prevent misreading by readers who scan the heatmap before reading the caption.

5. **(Moderate) Clarify what "enrichment" means in `env_cooccurrence.csv` and NB06 display.** The column records `n_both / expected`, but values like 1.01× and 1.02× across all P×Metal environments look uniformly trivial even when the log-ORs and Fisher tests indicate real effects. Since REPORT.md correctly uses log-OR as the primary metric, consider making log-OR the primary display column in NB06 cell outputs and leading with it in the CSV, to avoid the enrichment values creating a misleading impression of negligible effect sizes.

6. **(Minor) Add a forest plot caption note about CI width reflecting sample size.** Several phyla with n<100 (Fusobacteriota n=59, Deinococcota n=52) have very wide CIs (e.g., Deinococcota log-OR=+4.65 [−0.16, +9.46], p=1.0) that can visually dominate the plot. A one-sentence caption note that "wide CIs reflect small phylum sample sizes and should not be interpreted as large effects" would help readers correctly calibrate their attention.

7. **(Minor) Elevate phylogenetic independent contrasts in Future Directions.** Of the three future directions listed, phylogenetic logistic regression is the most directly achievable with existing BERDL data (GTDB phylogenetic distances are derivable from the `gtdb_species_clade` table) and directly addresses the analysis's central methodological limitation (Limitation #4). The current ordering (Substrate C via new RB-TnSeq campaigns, Substrate B via NMDC multi-omics, PIC) buries the most tractable next step at the end of the list.

## Review Metadata

- **Reviewer**: BERIL Automated Review (Claude, claude-sonnet-4-6)
- **Date**: 2026-05-08
- **Scope**: README.md, RESEARCH_PLAN.md, REPORT.md, 7 notebooks (NB01–NB07), 8 source scripts (src/01–src/08), 10 data files, 3 figures, requirements.txt, beril.yaml, docs/pitfalls.md
- **Note**: This review was generated by an AI system. It should be treated as advisory input, not a definitive assessment.
