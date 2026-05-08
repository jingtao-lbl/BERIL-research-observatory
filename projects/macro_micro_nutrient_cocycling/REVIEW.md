---
reviewer: BERIL Automated Review (Claude, claude-sonnet-4-6)
date: 2026-05-08
project: macro_micro_nutrient_cocycling
---

# Review: Macro-Micro Nutrient Gene Co-cycling in Bacterial Pangenomes

## Summary

This is a strong, well-documented project that addresses a clear biogeochemical hypothesis — whether bacterial pangenomes coupling macro-nutrient acquisition to metal handling reflect the biogeochemical co-mobilization of phosphate and trace metals from Fe-oxyhydroxide surfaces — using an appropriate and large dataset (27,682 GTDB species pangenomes). The eight-step analysis is logically organized, the statistical methods are sound (Fisher's exact test, BH-FDR, permutation null, Stouffer meta-analysis), and the findings are presented with commendable nuance: the effect sizes are modest but interpretable, the limitations section is unusually thorough (nine limitations), and the positive-control validation table is a highlight. The main areas for improvement are (1) the notebooks function as post-hoc data display wrappers rather than executable analysis documents, widening the gap between "notebooks have outputs" and true end-to-end reproducibility; (2) a gene name case-sensitivity bug acknowledged in Limitation 9 is not corrected in the source code; and (3) the 100% overlap finding for Phz-operon × Metal, while real, deserves more careful framing because the Fisher p=9.4×10⁻³ is only marginally significant and follows naturally from the near-universal (91.7%) baseline prevalence of metal genes.

## Methodology

**Research question and hypothesis:** The question is clearly stated and testable: positive co-occurrence (H1) vs. chance co-occurrence (H0) for three gene-group pairs across species pangenomes. The framing is appropriately grounded in mechanistic literature (McRose & Newman 2021; Edayilam et al. 2018), and the prediction that soil/rhizosphere lineages show stronger co-occurrence is specific enough to serve as a directional test.

**Approach:** The four-step plan (define gene families → co-occurrence matrix → core/accessory enrichment → phylogenetic stratification) is well-matched to the question. The addition of three supplemental steps (positive controls, environmental stratification, forest plot) strengthens the case considerably. The ≥3 distinct phz gene threshold for defining "operon carriers" is clearly motivated and the sensitivity analysis (≥2: 1,125; ≥3: 63; ≥4: 34 species) provides reassurance that the threshold does not drive the result.

**Data sources:** `kbase_ke_pangenome` is the appropriate primary source. The 27,682-vs-27,702 species count discrepancy (plan estimate vs. actual) is correctly noted and resolved. The PFAM-annotated lookup for phoD, feoB, and HMA is correct, and the decision to exclude PFAM PF00142 from primary N-fixation counts (on the grounds that it captures the broader Fer4 ferredoxin superfamily) is well-justified with quantitative evidence (5,872 vs. 2,746 species).

**Reproducibility:**  The README `## Reproduction` section clearly identifies which steps require Spark (01, 07) and which run locally on cached CSVs (02–06, 08). This is exemplary Spark/local separation. All intermediate CSV files are committed to `data/` (18 files, 19 MB total), so downstream steps can be re-run without cluster access.

One gap: the notebooks (NB01–NB07) are display wrappers, not primary analysis vehicles. NB01's code cells load from `data/species_gene_families.csv` rather than issuing Spark queries; the actual Spark query lives in `src/01_extract_gene_families.py`. This architecture is documented in the markdown cells and the README, but means that running `jupyter nbconvert --execute NB01_gene_family_extraction.ipynb` reproduces only the summary statistics display, not the data extraction itself. A reader who expects notebooks to be the analysis entry point will be surprised. Additionally, NB01 has one code cell without output (Cell 2, the `import` cell), and NB02–NB07 each have one code cell without output, which is minor.

## Code Quality

**SQL correctness:** The PFAM queries use `LIKE 'PFxxxxx.%'` versioned-ID matching, correctly avoiding the pitfall of bare-ID matching against versioned domain entries. The main gene_cluster × bakta_annotations JOIN and GROUP BY are correct.

**Case sensitivity (Limitation 9):** Gene-name-based families (pstA, pstB, pstC, pstS, phnC, phnD, phnE, corA, and all phz genes) are queried with exact case matching, e.g. `ba.gene = 'pstA'`. Bakta gene name capitalization is not fully standardized across genomes and releases. The project correctly identifies this as a limitation (Limitation 9 in REPORT.md) but does not fix it. A case-insensitive match (`LOWER(ba.gene) = 'psta'`) or `ba.gene IN ('pstA', 'PstA', 'PSTA')` for the known variants would eliminate this source of under-counting. This is the one code issue that is acknowledged but left unresolved.

**PFAM-family separation in src/01:** PFAM-based families (phoD_pfam, nifH_pfam, feoB_pfam, HMA_pfam) are detected via EXISTS subqueries, so they are skipped by the `if 'EXISTS' in condition: continue` branch in the main query builder and handled via separate PFAM queries (lines 76+). This split is functional but fragile — if a new gene family uses a different EXISTS pattern it could be silently skipped. A cleaner design would segregate the two query types at the family definition level (e.g., `{'type': 'pfam', 'id': 'PF09423'}` vs. `{'type': 'ko', 'id': 'K01077'}`), but this is a style preference and not a bug.

**Statistical methods:** Fisher's exact test is appropriate for 2×2 co-occurrence contingency tables. BH-FDR correction across 72 pairwise tests is correct. The permutation null (shuffling one group vector while holding the other fixed, 1,000 iterations, seed=42) correctly tests independence while preserving marginal prevalences. The Stouffer meta-analysis (`Z_meta = Σz / √n`) is unweighted, which treats each species as an equally informative observation; a sample-size-weighted Stouffer would better account for variation in per-species gene counts but is unlikely to change the direction of results given the large n.

**Notebook organization:** NB01 → NB07 follow the logical flow setup → query summary → analysis → visualization. Within each notebook, markdown cells provide clear context.

**Pitfall awareness:** The project correctly avoids the ENIGMA strain-name collision pitfall (not applicable here), uses versioned PFAM IDs as required, and appropriately uses `kbase_ke_pangenome` (not the deprecated `kbase_ke` prefix). No pitfalls from `docs/pitfalls.md` are violated.

## Findings Assessment

**P × Metal and N × Metal results:** The conclusions are well-supported by the reported statistics. Phi=0.110 and OR=2.3 for P × Metal, with permutation Z=17.7, constitutes a robust positive association given N=27,682. The OR=10.1 for N × Metal, driven by near-universal metal gene prevalence among diazotrophs (99.0%), is appropriately highlighted as the per-species strongest signal. The negative pstC/S–feoB associations (OR=0.64–0.72, phi as low as −0.256) are well-explained by the redox geochemistry of aerobic vs. anaerobic iron forms.

**Phz-operon × Metal:** The "100% overlap" framing deserves additional care. Given that 91.7% of all 27,682 species encode at least one metal-handling gene, the probability that all 63 phenazine operon carriers do so by chance is approximately 0.917^63 ≈ 0.004 — modest but not dramatic. The Fisher p=9.4×10⁻³ is consistent with this calculation and reflects that the effect, while real, is only marginally significant. The permutation Z=2.3 barely exceeds the conventional Z>2 threshold. The REPORT presents this as a "complete overlap" finding, which is accurate, but should note that the absolute test is whether 63/63 is surprising given the 91.7% background, not whether 63/0 is surprising in isolation. The OR=∞ label (arising from the empty n10 cell) can mislead readers into thinking the effect is larger than it is.

**Core/accessory enrichment:** The three-signature classification (P-acquisition: core-enriched; N-fixation: moderately core via HGT; metal-handling: bimodal) is well-supported by the Stouffer Z values and is the most conceptually novel finding in the analysis.

**Phylogenetic stratification:** The phylum-level distribution (positive phi in 28/34 phyla for P × Metal; negative in Myxococcota and Spirochaetota) is presented clearly. The ceiling-effect explanation for plant-associated families (100% prevalence of both groups → enrichment ratio 1.00×) is correct.

**Positive controls:** The positive-control table is a strong element. The *Pseudomonas fluorescens* nifH footnote (likely a divergent Fer4 ferredoxin, not true nitrogenase) is an honest acknowledgment of annotation ambiguity. The *S. coelicolor* → *S. anthocyanicus* GTDB reclassification is correctly identified and linked to a genome accession.

**Limitations:** The nine-limitation section is unusually comprehensive. The Actinomycete phenazine under-detection caveat (Limitation 7) is particularly well-argued and the *S. coelicolor* positive control example is apt.

**Incomplete analysis / placeholders:** No cells are left as "to be filled." All eight source scripts have corresponding data outputs. Figures 1, 2, and 3 are present and referenced correctly.

## Suggestions

1. **Fix gene name case sensitivity in src/01** *(important)*  
   Convert gene name comparisons to lowercase for all `ba.gene`-based lookups: replace `ba.gene = 'pstA'` with `LOWER(ba.gene) = 'psta'` (and similarly for pstB, pstC, pstS, phnC, phnD, phnE, corA, phzA, phzB, phzD, phzG, phzM). This eliminates the acknowledged Limitation 9 and avoids under-counting genuinely present genes. After fixing, re-running src/01 (requires Spark) and comparing result counts to the current species_gene_families.csv would quantify the impact.

2. **Clarify the Phz-operon × Metal "100% overlap" result** *(important)*  
   Add a sentence to the Results section noting that the 100% overlap is tested against the 91.7% background metal-handling prevalence and that the Fisher p=9.4×10⁻³ and permutation Z=2.3 reflect this context. Consider moving the Phz × Metal row to a separate paragraph with explicit language like: "While all 63 operon carriers encode metal-handling genes, the statistical significance of this 100% overlap (p=9.4×10⁻³, permutation Z=2.3) is modest given the near-universal background prevalence of metal genes (91.7%)." This does not change the finding — only its framing.

3. **Clarify notebook role in the Reproduction section** *(important)*  
   Add a note to README.md (or to the opening markdown cell of each notebook) explaining that NB01–NB07 are inspection/display notebooks that load pre-computed CSVs from `data/`, and that the primary analysis code is in `src/*.py`. For example: "NB01 loads the CSV produced by `src/01_extract_gene_families.py` (Spark step). Re-running NB01 does not re-query BERDL." This sets correct expectations for readers and reviewers.

4. **Add Bakta gene-name case sensitivity to docs/pitfalls.md** *(nice-to-have)*  
   Limitation 9 documents a data gotcha that affects any project querying `bakta_annotations.gene` by exact string match. This is a reusable finding that belongs in the shared pitfalls log so future projects benefit from it. A concise entry with the fix (`LOWER(ba.gene) = lower_name`) would be valuable.

5. **Cross-reference the phenazine broad vs. operon group consistently** *(nice-to-have)*  
   The REPORT introduces the "broad Phenazine group" (10,308 species with any phz gene ≥1) and the "Phenazine operon" group (63 species with ≥3 phz genes) in separate sections, but the phi=0.262 result for the broad group appears only in a note without a direct comparison to the operon group in the main results table. A brief comparison column in Table 1 (or a note in the table caption) would help readers track which Phenazine definition is being reported.

6. **Consider weighted Stouffer for core/accessory meta-analysis** *(nice-to-have)*  
   The current unweighted Stouffer meta-analysis (`Z_meta = Σz / √n`) weights each species equally. Weighting by per-species sample size (number of gene clusters tested) using `Z_meta = Σ(w_i z_i) / √(Σw_i²)` with `w_i = √n_i` would give more statistical leverage to species with large, well-sampled pangenomes. Given the large n and strong results (Stouffer Z=68.3 for P-acquisition), this is unlikely to change conclusions but would strengthen the meta-analytic claim.

## Review Metadata

- **Reviewer**: BERIL Automated Review (Claude, claude-sonnet-4-6)
- **Date**: 2026-05-08
- **Scope**: README.md, RESEARCH_PLAN.md, REPORT.md, beril.yaml, 7 notebooks (NB01–NB07), 18 data files, 4 figures, src/01–02_*.py (inspected), requirements.txt, docs/pitfalls.md
- **Note**: This review was generated by an AI system. It should be treated as advisory input, not a definitive assessment.
