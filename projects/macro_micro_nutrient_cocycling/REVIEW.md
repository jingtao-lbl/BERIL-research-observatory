---
reviewer: BERIL Automated Review (Claude, claude-sonnet-4-6)
date: 2026-05-13
project: macro_micro_nutrient_cocycling
---

# Review: Macro-Micro Nutrient Gene Co-cycling in Bacterial Pangenomes

## Summary

This is a mature, rigorous analysis that successfully tests whether bacterial pangenomes co-encode macro-nutrient acquisition and trace-metal handling genes at rates exceeding chance. Working across 27,682 GTDB species pangenomes from `kbase_ke_pangenome`, the project builds a genuinely impressive multi-layered evidentiary case: group-level co-occurrence with Fisher's exact test and permutation null, FDR-corrected pairwise analysis of 72 gene pairs, Stouffer meta-analysis of core/accessory structure, phylum- and environment-stratified effect sizes, phylogenetic logistic regression on the full 26,517-tip GTDB R214 tree, an operon-distance test that rules out physical gene linkage, and an independent validation via Wang et al. 2021. All 14 significant pairs survive phylogenetic correction (0% collapse rate). The REPORT.md is detailed, internally consistent, and thoroughly acknowledges limitations. The project's primary weaknesses are reproducibility gaps: the README's `## Reproduction` section documents only Steps 1–8 (src/01–08), leaving the entire extended analysis (src/09–15: R-based phylogenetic correction, operon-distance test, Wang validation, and their figures) unreachable by the reproduction guide, and there is no R environment specification to match the Python `requirements.txt`. These gaps do not undermine the scientific conclusions, but they would prevent an independent researcher from reproducing the full pipeline.

---

## Methodology

**Research question and hypothesis**: Both are precisely stated and testable. H1 ("presence of P-acquisition, N-fixation, and phenazine biosynthesis gene families is positively associated with metal-handling gene families, beyond a permutation null") is directly tested by the analytical pipeline. The literature context (McRose & Newman 2021, Edayilam et al. 2018) cleanly motivates the hypothesis.

**Analytical approach**: The layered strategy is well-suited to the question. Using both permutation and parametric tests as complementary checks is appropriate. The phylogenetic logistic regression choice (`phylolm::phyloglm(method="logistic_MPLE")`) is sound for binary traits on a large tree — MPLE avoids the O(n²) covariance matrix that would otherwise make the full 26,517-tip tree computationally intractable. The decision to estimate Pagel's lambda on a 5,000-tip stratified subsample (due to memory constraints) rather than the full tree is a pragmatic limitation, transparently documented.

**Gene family definitions**: The KEGG KO / gene-name / PFAM-prefix layered annotation approach is reasonable and appropriate. The distinction between primary (KO-only) and broad (PF00142) N-fixation definitions — and the explicit exclusion of the broader definition from primary analysis — is scientifically sound and properly justified. The ≥3-gene phenazine operon threshold is explained and sensitivity-checked.

**Data sources**: `kbase_ke_pangenome` (132.5M gene clusters, 27,682 species, GTDB R214) is clearly identified. Table-level details (gene_cluster, bakta_annotations, bakta_pfam_domains, gtdb_species_clade) are named.

**Reproducibility of methods**: Results in REPORT.md are cross-checked against NB01–NB07 outputs and data/CSV files — the numbers are consistent throughout.

---

## Code Quality

**SQL queries** (`src/01_extract_gene_families.py`): The CASE WHEN / GROUP BY pattern for generating the species-level presence/absence matrix is correct. The two-pass approach (gene-name / KO pass, then separate PFAM subqueries) avoids a three-way join that would be prohibitively expensive at 132.5M gene clusters. The `LIKE 'PFxxxxx.%'` pattern for PFAM versioning is consistent with the pitfalls.md guidance on versioned PFAM IDs. The `fillna(0)` call after merging PFAM results is appropriate.

**Statistical methods** (`src/02_cooccurrence_stats.py`): Fisher's exact test, phi coefficient, Jaccard, and permutation null are all implemented correctly. The BH FDR implementation (manual, step-down) is correct. The Stouffer Z meta-analysis (`src/03`) is appropriately applied as a directional aggregate test.

**R scripts** (`src/09_phylo_signal.R`, `src/10_phylo_logistic.R`): The phylogenetic signal estimation with stratified-by-phylum subsampling (src/09) is well-designed — it avoids both overrepresenting large phyla and discarding small ones. The phyloglm call uses `btol=50, log.alpha.bound=8` to improve convergence stability, and convergence is checked in the output. Both R scripts accept command-line arguments cleanly, enabling batch execution.

**Operon-distance test** (`src/11_operon_distance.py`): The `parse_gene_id` function correctly extracts the ordinal from the `{contig}_{ordinal}` BERDL gene_id format. Cross-contig pairs are excluded from distance computation but included in the same-contig fraction denominator. The permutation correctly shuffles ordinals within each species (preserving contig membership) rather than globally, which is the appropriate null.

**Notebook organization** (NB01–NB07): All seven notebooks follow a consistent pattern: markdown header describing what the notebook does and what script to run for re-analysis, data load from cached CSV, display / summary. Cells are logically ordered. NB01 and NB02 both have full text outputs. One minor code quality issue: in NB04 cell_3, `from matplotlib.patches import Patch` is imported once at the top of the cell and then imported again inside the Panel D block — a redundant import that causes no harm but could be cleaned up.

**Pitfall awareness**: The project correctly uses `LIKE 'PFxxxxx.%'` for PFAM ID matching (avoiding the fixed-version pitfall). The gene-name case sensitivity issue (e.g., `pstA` vs `PstA`) is identified as Limitation 9 in the REPORT — but no sensitivity check quantifying the potential miss rate was performed. Given that pst and phz families are gene-name-based (no KO backup for pstA/B/C/S), a simple case-insensitive comparison or a spot-check of alternative capitalizations in the database would strengthen confidence in the coverage estimates.

---

## Findings Assessment

**Conclusions supported by data**: Yes. The group-level co-occurrence table in REPORT Section 1 is reproduced verbatim from NB02 output. All phi, OR, and Fisher p values match the CSV outputs. The phylogenetic correction results (Table 5) are consistent with data/phylo/phylo_logistic.csv. The 100% phenazine operon × metal overlap is confirmed in NB01. The operon-distance result (median 910 vs null 249, Z=304) is confirmed by data/operon_distance/permutation_summary.csv and observed_dists.npy.

**Effect size interpretation**: The REPORT's discussion of "modest phi at the pan-bacterial scale" is appropriately nuanced. The three-level effect size framing (group phi → per-pair OR → environment-stratified log-OR) is well-reasoned and avoids over-claiming.

**Limitations acknowledged**: Thorough. Nine limitations are listed covering expression-level inference, annotation sensitivity, PFAM PF00142 inflation, phylogenetic subsampling, ecological vs expression inference, phenazine operon definition, actinomycete under-detection, environmental classification keyword priority, and gene name case sensitivity. The actinomycete phenazine under-detection limitation (Limitation 7) is unusually candid and well-supported by the S. coelicolor positive control.

**Negative associations**: The pstC × feoB depletion (phi=−0.256, strengthened by phylogenetic correction) is correctly interpreted as a genuine ecological niche separation rather than an artifact, and the redox-geochemical mechanism is plausibly argued.

**Incomplete or placeholder analysis**: None. All sections in REPORT.md are fully filled in with actual results.

**Wang 2021 validation**: The non-replication of the Burkholderiaceae-specific signal is clearly explained (scope difference: soil isolate enrichment vs. pan-GTDB) and is not over-claimed as a contradiction of Wang et al. 2021.

**One minor inconsistency**: REPORT.md states the P × Metal animal-associated result is "q=0.54" while the NB06 BH-corrected output shows `q=0.543`. These match. However, the REPORT's Table 4 gives `q` for soil/rhizosphere P×Metal as `8.1×10⁻⁹`, while NB06 shows `q=0.000` (rounded). The NB06 FDR computation rounds `q<0.001` to 0.000; the REPORT correctly converts this to `8.1×10⁻⁹` using the raw Fisher p — this is consistent but the derivation path is not stated.

---

## Suggestions

1. **[Critical] Extend the `## Reproduction` section in README.md to cover src/09–15.** Steps 9–15 (phylogenetic signal, phylogenetic logistic regression, operon-distance test, Wang validation, and the associated figures) are entirely absent from the reproduction guide. A reader following the README would reproduce only the first half of the analysis. The extension should specify: (a) that src/09 and src/10 require R with the `r_phylo` conda environment, (b) the command-line arguments for the R scripts (tree path, trait path, pairs CSV, output path), (c) that src/11 requires on-cluster Spark, and (d) that src/12–15 run locally from cached CSV outputs. Example sketch:
   ```bash
   # Step 9: Phylogenetic signal (R, r_phylo env, ~20 min)
   conda run -n r_phylo Rscript src/09_phylo_signal.R \
     data/phylo/bac120_r214_pruned.tree data/species_gene_families.csv \
     data/species_taxonomy.csv data/phylo/phylo_signal.csv 5000

   # Step 10: Phylogenetic logistic regression (R, r_phylo env, ~2 hr)
   conda run -n r_phylo Rscript src/10_phylo_logistic.R \
     data/phylo/bac120_r214_pruned.tree data/species_gene_families.csv \
     data/phylo/pair_definitions.csv data/phylo/phylo_logistic.csv

   # Step 11: Operon-distance test (requires Spark)
   python src/11_operon_distance.py

   # Steps 12–15: Wang validation and figures (local)
   python src/12_wang2021_validation.py
   python src/13_figure4_phylo_correction.py
   python src/14_figure5_operon_distance.py
   python src/15_figure6_wang2021.py
   ```

2. **[Important] Add an R environment specification.** The Python `requirements.txt` lists Python package versions, but the R analysis (src/09, src/10) depends on `phylolm` (v2.6.5) and `phytools` (v2.3-0) — version-sensitive packages that can produce different numerical results across versions. Add either `data/phylo/r_session_info.txt` (output of `sessionInfo()`) or an `r_requirements.txt` file listing R package versions. The REPORT already mentions these version numbers; surfacing them as a reproducible artifact closes the gap.

3. **[Important] Add display notebooks for the extended analyses (NB08+).** The existing NB01–NB07 cover src/01–08 and each loads pre-computed CSVs for display. Applying the same pattern to the extended analyses would provide inspection points for the phylogenetic correction (data/phylo/phylo_logistic.csv, data/phylo/phylo_signal.csv), operon-distance permutation (data/operon_distance/), and Wang validation (data/wang2021/). Even minimal notebooks (load CSV → display table → show figure) would allow reviewers to verify results without re-running multi-hour R jobs. This is also the pattern that the pitfalls.md entry on "Commit Notebooks Alongside Their Artifacts" explicitly recommends.

4. **[Important] Rename `figures/forest_plot.png` to `figures/figure2_forest_plot.png`** (and update REPORT.md and NB07 references accordingly) to match the numbered figure naming scheme used by all other figures. Currently figures are named figure1_*, figure3_*, figure4_*, figure5_*, figure6_* — the forest plot is the odd one out. This is a minor naming issue but can confuse readers navigating the figures/ directory.

5. **[Minor] Quantify the gene-name case sensitivity miss rate.** Limitation 9 correctly identifies that `ba.gene = 'pstA'` uses exact case matching. A one-off check of how many gene clusters in `bakta_annotations` match `LOWER(ba.gene) = 'psta'` but not `ba.gene = 'pstA'` would bound the potential miss rate. If the rate is <1%, that is worth stating explicitly to reassure readers. If it is non-trivial, a case-insensitive query (`LOWER(ba.gene) = 'psta'`) should be used in src/01.

6. **[Minor] Update `beril.yaml` status from `analysis` to `complete`.** The current status field reads `status: analysis`, which typically indicates in-progress work. Given the project is described as "Completed" in README.md and has a finalized REPORT.md, the beril.yaml field should reflect the project's actual state.

7. **[Minor] Remove the redundant `from matplotlib.patches import Patch` import in NB04 cell_3.** The import appears once at the top of the cell (line ~4) and again inside the Panel D code block. The second import is unreachable duplication.

8. **[Nice-to-have] Extend the operon-distance test to N-fixation × Metal pairs.** The current test covers only P × M gene pairs (src/11 CASE WHEN block includes pstA/B/C/S, phnC/D/E, and metal families, but omits nifH/nifD). Since N × Metal shows the strongest pan-bacterial OR (10.1), verifying that nifH/nifD are also not physically co-located with metal-handling genes would further strengthen the "ecological, not genomic linkage" conclusion. This would require a small extension to the bakta_targets temp view in src/11.

---

## Review Metadata

- **Reviewer**: BERIL Automated Review (Claude, claude-sonnet-4-6)
- **Date**: 2026-05-13
- **Scope**: README.md, RESEARCH_PLAN.md, REPORT.md, references.md, beril.yaml, requirements.txt, 7 notebooks (NB01–NB07 with cell outputs), 15 source files (src/01–15), docs/pitfalls.md, 6 figures, 18 data files/directories
- **Note**: This review was generated by an AI system. It should be treated as advisory input, not a definitive assessment.
