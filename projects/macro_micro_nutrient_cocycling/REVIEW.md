---
reviewer: BERIL Automated Review (Claude, claude-sonnet-4-6)
date: 2026-05-08
project: macro_micro_nutrient_cocycling
---

# Review: Macro-Micro Nutrient Gene Co-cycling in Bacterial Pangenomes

## Summary

This is a rigorous, well-executed genome-scale analysis of macro-nutrient and trace-metal gene co-occurrence across 27,682 GTDB species pangenomes. The project is thoroughly documented: the research question is specific and testable, the statistical approach is multi-layered (Fisher's exact test, phi coefficient, Jaccard index, permutation null, BH FDR correction, Stouffer meta-analysis, phylogenetic stratification, and environmental stratification), and the limitations section is unusually candid and detailed. The src/notebook split — where `src/*.py` scripts perform the primary analysis and notebooks display pre-computed CSVs — is unconventional but clearly documented in the README and each notebook header. All seven notebooks have saved outputs, three publication-quality figures exist in `figures/`, and the reproduction path is clearly delineated (Spark-required steps 1 and 7 vs. locally-runnable steps 2–6 and 8). The main areas for improvement are: (1) a subtle inconsistency between the notebook enrichment ratios and the log-ORs reported in the REPORT, (2) several gene families included in the analysis (pstB, phnD) are omitted from Figure 1 Panel B without explanation, (3) the 1,000-permutation limit (p_min ≈ 0.001) means the most significant signals are reported identically and conflated in the permutation column, and (4) phylogenetic non-independence remains the most consequential unaddressed confound.

---

## Methodology

**Research question and hypothesis:** Both are clearly stated and precisely operationalised. H1 (positive co-occurrence exceeding a permutation null) and H0 (no excess co-occurrence) are testable and are tested. The pivot from `kescience_fitnessbrowser` (no P-starvation experiments) to `kbase_ke_pangenome` is documented in the revision history.

**Data source:** `kbase_ke_pangenome` with 132.5M gene clusters and GTDB R214 taxonomy. Tables used (`gene_cluster`, `bakta_annotations`, `bakta_pfam_domains`, `gtdb_species_clade`) are identified. The 20-species gap between the 27,702 GTDB total and 27,682 retrieved is noted (species lacking any Bakta annotations, confirmed by the JOIN structure in `src/01`).

**Gene family definitions:** The 23 families use a sensible mix of annotation types (KEGG KO for nifH/nifD/copA/phzF/phzS, PFAM for phoD/nifH_pfam/feoB/HMA, gene name for pst/phn/corA/phzA–M). The sensitivity check comparing KO-only nifH/nifD (2,746 species) against PF00142-inclusive (5,872 species) is a highlight — it is reported transparently, its biological meaning is explained, and it correctly motivates excluding PF00142 from the primary analysis.

**Phenazine operon threshold:** The ≥3 distinct phz gene heuristic is justified empirically (≥2 gives 1,125 species; ≥4 gives 34), robustly documented in Limitation #6, and acknowledged as anchored on the canonical *Pseudomonas* operon at the expense of Actinomycete detection (Limitation #7).

**Statistical approach:** Permutation null (1,000 shuffles, seed 42) is appropriate. BH FDR correction is correctly implemented in both `src/02` (ascending-rank step-down) and `src/07` (descending step-down). Stouffer meta-analysis of per-species Fisher tests is appropriate for the core/accessory question, though the nuances of its interpretation for the Phz group (outlier-driven aggregate) are well explained in the REPORT.

**Reproducibility:** The README clearly separates Spark-required from local steps, explains the src/notebook architecture, and lists expected outputs. A reader can reproduce the local steps from the cached CSVs without cluster access. The `requirements.txt` lists all dependencies (numpy, pandas, scipy, matplotlib, pyspark, berdl_notebook_utils).

---

## Code Quality

**SQL correctness (src/01):** The main query skips PFAM gene families (`if 'EXISTS' in condition: continue`) and handles them in separate PFAM queries merged back via `pd.DataFrame.loc`. This two-pass approach is correct. PFAM IDs use `LIKE 'PFxxxxx.%'` to match versioned IDs, which is the correct pattern per the BERDL schema. Species absent from the PFAM query naturally receive 0 via `fillna(0)`.

**Potential issue — pstB/phnD omitted from Figure 1 Panel B:** `src/05_figure.py` and `NB04` both build Panel B (core fraction bars) using `gene_order = ['phoA', 'phoD_pfam', 'pstA', 'pstC', 'pstS', 'phnC', 'phnE', ...]`, silently dropping `pstB` and `phnD` from the visualisation. Both genes are included in the analysis (`groups['P_acquisition']` in `src/02`), and `NB03 cell_4` prints their core fractions. The omission is not explained in a figure caption or code comment. The REPORT and RESEARCH_PLAN describe 9 P-acquisition families but the figure shows 7. Readers comparing the figure to the text will notice the discrepancy.

**Gene name case sensitivity (Limitation #9):** Appropriately disclosed. The impact is limited to the gene-name-based families (pst, phn, phz, corA); KO and PFAM families are unaffected.

**Enrichment vs. log-OR discrepancy in NB06:** The notebook's printed summary table shows very low enrichment ratios (e.g., human-associated P×Metal: `enrich=1.02x`) alongside large log-ORs (`log-OR=+1.374`). These are not inconsistent — the enrichment ratio (n_both/expected_marginal) reflects ceiling effects when both gene families approach 100% prevalence, while log-OR captures the association structure in the 2×2 table — but a reader seeing only the notebook output will find the two metrics confusing without an explanatory note. The REPORT correctly uses log-OR for environmental comparisons, but the notebook would benefit from a markdown cell explaining why `enrich≈1.0x` does not contradict `log-OR≈1.4`.

**Permutation resolution:** With N=1,000 permutations, the minimum representable p is 1/1,001 ≈ 0.001. This is disclosed in the REPORT. However, the permutation p-values for P×Metal (0.0010), N×Metal (0.0010), and P×N (0.0010) are reported identically — all three "max out" the test. The REPORT correctly routes readers to the permutation Z-scores instead, but the permutation p column in the main results table adds little information for signals of this magnitude. A note flagging this explicitly in the REPORT table footnote or in NB02 would help prevent misreading the identical values as equivalent evidence strength.

**Stouffer meta-analysis for Phz group (src/03):** The REPORT provides an unusually detailed and correct explanation of why the Stouffer Z=2.6 for Phz genes is outlier-driven and should not be interpreted as indicating that phenazine genes are generally more core than metal genes (median OR=0.333, frac_enriched=0.268). This caveat is hard to find: it is buried deep in the Interpretation section. A brief caution note in `NB03` alongside the summary table output would improve discoverability.

**Pitfall check:** The most relevant documented pitfall (gene name case sensitivity) is directly acknowledged as Limitation #9. The project does not use ENIGMA strain matching, fitnessbrowser, or Web of Microbes, so those pitfalls do not apply.

---

## Findings Assessment

**Are conclusions supported?** Yes. The central finding — significant positive co-occurrence of P-acquisition and N-fixation genes with metal-handling genes (P×Metal: phi=0.110, OR=2.3, perm Z=17.7; N×Metal: phi=0.088, OR=10.1, perm Z=14.7) — is supported by multiple complementary tests. The 64/72 FDR-significant individual gene pairs strengthen the group-level result. The positive control species check provides external biological validation.

**Phenazine operon finding (100% overlap, 63/63):** The REPORT correctly contextualises this against the 91.7% background prevalence of metal-handling genes and provides the back-of-envelope calculation (0.917^63 ≈ 0.004). The Fisher p=9.4×10⁻³ and permutation Z=2.3 are statistically real but modest, which is honestly reported. The REPORT's caution ("statistically significant but modest in magnitude") is appropriate.

**Negative associations:** The pstC/S–feoB depletion (phi as low as −0.256) is the most concrete result in the analysis and is given appropriate ecological interpretation (aerobic/oligotrophic Pst users vs. microaerobic FeoB users). The phzF–feoB anti-correlation (phi=−0.163) is a genuinely novel observation that reinforces the McRose-Newman model without over-claiming.

**Core/accessory findings:** The three-genomic-strategy framework (P-acquisition: core-stable; Phz: clade-specific; N-fixation: HGT-flexible) is interpretively sound and maps onto known biology. The phoD exception (17.2% core, consistent with mobile-element carriage) adds specificity.

**Phylogenetic non-independence:** Limitation #4 correctly identifies this as the key unaddressed confound. The permutation test controls for marginal frequencies but not phylogenetic autocorrelation. The REPORT notes this and lists phylogenetic logistic regression as the first Future Direction. The current analysis is scientifically meaningful at genome-discovery scale, but downstream work building on this result will need to address it.

**Environmental stratification:** The finding that P×Metal is strongest in plant-associated and soil environments while N×Metal is strongest in marine and animal-associated environments is the most ecologically interpretive result and is directly testable against field data. The REPORT is appropriately cautious about the large "other" category (12,628 species, 46%) and the keyword-priority misclassification risk (Limitation #8).

**Limitations acknowledged:** All major limitations are enumerated (9 items). This is a strength.

---

## Suggestions

1. **[Figure 1, Panel B — medium priority] Document the exclusion of pstB and phnD.** Either include them in the panel or add a code comment and figure caption note explaining they were omitted for space. Readers familiar with the 9-family P-acquisition definition will notice the 2-gene gap between the text and the figure.

2. **[NB06 — medium priority] Add a markdown cell explaining the enrichment-vs-log-OR discrepancy.** A single sentence — e.g., "Enrichment ratios near 1.0 are expected when both gene groups approach universal prevalence; log-OR better captures the association structure in the full 2×2 table" — would prevent the apparent contradiction from confusing readers inspecting the notebook output.

3. **[REPORT.md / NB02 — low priority] Flag the permutation resolution saturation.** The three pairs with perm_p=0.0010 saturate the test and are indistinguishable in that column. A table footnote such as "Permutation p bounded at 1/(N+1)≈0.001 for N=1,000 permutations; the Z-scores (17.7, 14.7, 10.7) differentiate these signals" would prevent readers from treating the identical p-values as equivalent strength of evidence.

4. **[NB03 — low priority] Surface the Stouffer Phz interpretation caveat closer to the table.** The explanation of why Stouffer Z=2.6 for Phz genes is outlier-driven (median OR=0.333, frac_enriched=0.268) is currently only in the Interpretation section of REPORT.md. A one-sentence markdown cell below the `core_enrichment_summary` table in NB03 would catch readers who inspect notebooks but not the full REPORT.

5. **[src/01 — low priority] Document the two-pass PFAM strategy inline.** The `if 'EXISTS' in condition: continue` skip is not explained with a comment. A brief note explaining that PFAM subquery families are handled in separate queries and merged back (for Spark query plan reasons) would help future developers understand the architecture.

6. **[NB07 — low priority] Print n_species alongside log-OR in the top-effects summary.** Several phyla have p=1.00e+00 with large log-OR point estimates (e.g., Deinococcota, +4.654). These reflect small sample sizes or extreme cell values, but the current printout does not show n. Adding n to the per-phylum summary would clarify whether these are underpowered estimates rather than genuine large effects.

7. **[Future work — informational] Quantify the phylogenetic non-independence concern.** Limitation #4 correctly identifies the problem but gives readers no sense of scale. The Future Directions section could suggest a specific bounded approach: e.g., applying Moran's I to the species-level binary vectors using the GTDB phylogenetic distance matrix, to estimate what fraction of the observed phi=0.110 is attributable to shared ancestry before investing in full phylogenetic logistic regression.

---

## Review Metadata

- **Reviewer**: BERIL Automated Review (Claude, claude-sonnet-4-6)
- **Date**: 2026-05-08
- **Scope**: README.md, RESEARCH_PLAN.md, REPORT.md, 7 notebooks (NB01–NB07), 8 source scripts (src/01–08), 13 data files, 3 figures, requirements.txt, beril.yaml, docs/pitfalls.md
- **Note**: This review was generated by an AI system. It should be treated as advisory input, not a definitive assessment.
