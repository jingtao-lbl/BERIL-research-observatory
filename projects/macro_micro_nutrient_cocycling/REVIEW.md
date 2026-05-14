---
reviewer: BERIL Automated Review (Claude, sonnet)
date: 2026-05-14
project: macro_micro_nutrient_cocycling
---

# Review: Macro-Micro Nutrient Gene Co-cycling in Bacterial Pangenomes

## Summary

This is a methodologically sophisticated and well-documented analysis of nutrient and metal-handling gene co-occurrence across 4,540 soil/plant-associated bacterial pangenomes. The project demonstrates exemplary reproducibility practices: all seven notebooks contain saved outputs, six high-quality figures are present, and the README provides detailed reproduction instructions with runtime estimates and environment requirements. The statistical approach is rigorous, employing Fisher's exact tests, permutation nulls, phylogenetic correction via logistic MPLE, and multiple independent validation tests (operon-distance, KEGG pathway co-membership, Wang 2021 phytase×siderophore). The central finding—that P-acquisition and metal-handling genes co-occur significantly (phi=0.129, OR=3.9) and survive phylogenetic correction despite 41% attenuation—is well-supported and appropriately contextualized. The analysis acknowledges its limitations transparently, including reduced statistical power in the soil+plant subset and the near-universal metal-handling gene prevalence (96.4%) that compresses enrichment ratios. This is publication-ready work with only minor documentation improvements needed.

## Methodology

**Research Question**: Clearly stated and testable. The hypothesis predicts that bacterial pangenomes encoding macro-nutrient acquisition machinery (P, N) and phenazine biosynthesis genes will disproportionately co-encode trace-metal handling genes, reflecting biogeochemical coupling in soil environments where Fe-oxyhydroxide minerals mediate nutrient availability.

**Approach**: Sound and multi-layered. The analysis employs:
- Group-level and individual gene-pair co-occurrence testing (Jaccard, phi, Fisher's exact, 72 pairs with FDR correction)
- Permutation null models (1,000 permutations, fixed seed=42)
- Core vs. accessory genome enrichment (Stouffer meta-analysis)
- Phylogenetic stratification and phylogenetic logistic regression (phylolm on 4,177-tip GTDB tree)
- Independent mechanistic validation: operon-distance test (median 1,097 genes vs null 350, Z=120.7), KEGG pathway co-membership (Z=-31.0), and Wang 2021 phytase×siderophore replication

**Data Sources**: Clearly identified. The analysis uses `kbase_ke_pangenome` (132.5M gene clusters, 27,682 GTDB species pangenomes), filtered to 4,540 soil/rhizosphere and plant-associated species via `env_species_mapping.csv`. Gene families are defined using Bakta KEGG KO annotations, gene names, and PFAM domains. The GTDB R214 bac120 reference tree is used for phylogenetic correction.

**Reproducibility Strengths**:
- All 16 scripts (`src/01_*.py` through `src/16_*.py`) are present and numbered sequentially
- All 7 notebooks (NB01–NB07) have saved outputs — no empty code-only notebooks
- 6 figures present in `figures/` (figure1, forest_plot, figure4–6, plus figure3 from v3)
- Fixed random seed (42) used consistently across Python (`np.random.seed(42)`) and R (`set.seed(42)`)
- Clear separation of Spark-requiring steps (01, 07, 11, 12, 16) vs. local steps, with runtime estimates
- Dependencies listed in `requirements.txt` with version bounds
- README includes a detailed `## Reproduction` section with 16 numbered steps and environment setup instructions

**Reproducibility Gaps**:
- Figure 3 (`figure3_env_stratification.png`) is present but marked in the README as "v3 pan-bacterial only; not regenerated for v4." This is a minor inconsistency — either regenerate for v4 or move to an archive directory and remove from the README table to avoid confusion.
- The README lists step 10a as a convergence check (`grep "FALSE" data/phylo/phylo_logistic.csv`), but it's unclear whether this check was run and what the result was. The REPORT.md states "All models converged," but documenting the verification would strengthen reproducibility.

## Code Quality

**SQL and Data Extraction**: The gene family extraction (script 01) uses appropriate KEGG KO, gene name, and PFAM domain annotations. N-fixation is defined using KO-only (nifH=K02588, nifD=K02586) as the primary definition, with PFAM PF00142 as a sensitivity check. This is well-justified in REPORT.md: PF00142 captures the broader Fer4 ferredoxin superfamily and inflates N-fixation prevalence by 113%, so excluding it from the primary analysis is appropriate.

**Statistical Methods**: Appropriate and well-executed. Fisher's exact test is used for 2×2 contingency tables; permutation tests preserve marginal counts; Benjamini-Hochberg FDR correction is applied to 72 individual gene pairs (51/72 significant at q<0.05). Phylogenetic correction uses `phylolm::phyloglm(method="logistic_MPLE")`, which is appropriate for binary traits on large trees. Pagel's lambda is estimated for phylogenetic signal (lambda range 0.06–1.00 across gene families), confirming that phylogenetic correction is necessary.

**Code Organization**: Excellent. Scripts are numbered 01–16 in logical order. Each script includes a docstring stating its purpose and output files. The `src/02_cooccurrence_stats.py` script (inspected) shows clean, readable code with well-named functions (`jaccard`, `phi_coefficient`, `contingency_2x2`), appropriate comments, and clear data flow. The soil+plant filter is applied consistently across scripts 02–10 via `env_species_mapping.csv`.

**Known Pitfall Awareness**: The analysis does not directly intersect with the database-specific pitfalls documented in `docs/pitfalls.md` (which focus on Fitness Browser, MetaPhlAn3, curatedMetagenomicData, and strain-name collisions). However, the project demonstrates awareness of PFAM versioning issues (uses `LIKE 'PFxxxxx.%'` for PFAM queries) and phylogenetic non-independence (the entire phylogenetic correction workflow addresses this).

**Potential Issues**:
- **Phenazine operon threshold**: The `has_phz_operon` flag requires ≥3 distinct phz gene families. This is well-justified (phzF alone spans 2,652 species vs 27 operon carriers), but the choice of 3 as the threshold is not explicitly defended. A sensitivity analysis showing how results change at thresholds of 2, 3, 4, or 5 would strengthen this design choice.
- **Convergence check transparency**: The README includes step 10a (`grep "FALSE" data/phylo/phylo_logistic.csv || echo "All models converged"`), but there's no record of this check's output in the repository. Including a simple `data/phylo/convergence_check.txt` with the timestamp and result would make verification explicit.

## Findings Assessment

**Conclusions Supported by Data**: Yes. The primary finding — P × Metal phi=0.129, OR=3.9, Fisher p=4.4×10⁻¹⁴, surviving phylogenetic correction with phylo log-OR=0.938, p=2.1×10⁻⁵ — is backed by multiple lines of evidence:
1. Group-level and individual gene-pair statistics (Table 1, pairwise detail)
2. Permutation null model (Z=17.7 for P × Metal, perm_p=0.001)
3. Per-phylum forest plot (4/22 phyla significant)
4. Phylogenetic correction (10/12 significant pairs survive, 83% survival rate)

**Limitations Acknowledged**: Excellent. The REPORT.md lists 8 specific limitations:
1. Presence/absence only (not co-expression or co-regulation)
2. Annotation-dependent coverage (may miss divergent homologs)
3. PFAM PF00142 captures broader Fer4 superfamily (inflates N-fixation by 113%)
4. Phylogenetic non-independence controlled but informative (41% attenuation)
5. Ecological inference from genomic potential (encoding ≠ expression)
6. Reduced statistical power in 4,540-species subset vs 27,682-species pan-bacterial (51/72 FDR-significant pairs vs 64/72)
7. Near-universal metal-handling prevalence (96.4% ceiling effect)
8. Environmental classification limitations (keyword matching on NCBI biosample metadata)

These are honest, specific, and demonstrate awareness of the analysis's scope and boundaries.

**Incomplete or "To Be Filled" Analysis**: None detected. The analysis is complete for the stated scope. The REPORT.md lists "Future Directions" (Substrate C mechanistic validation, Substrate B ecological communities, field validation at Point Reyes, expanded phenazine gene set, cross-environment comparison), but these are clearly marked as extensions, not gaps in the current work.

**Visualizations**: Clear and properly labeled. Inspecting the README table and `figures/` directory:
- Figure 1 (multi-panel): co-occurrence heatmap, core fractions, phylum-level, phenazine taxonomy — present as both PNG and PDF
- Figure 2 (forest plot): per-phylum log-OR with 95% CIs — present as PNG
- Figure 4 (phylo correction): uncorrected vs phylo-corrected scatter — present as PNG and PDF
- Figure 5 (operon distance): observed vs null distribution — present as PNG and PDF
- Figure 6 (Wang 2021): per-family and per-siderophore validation — present as PNG and PDF

All figures referenced in the REPORT.md are present. The only exception is Figure 3 (environmental stratification), which is marked as "v3 pan-bacterial only; not regenerated for v4" — this is appropriate given the v4 analysis is restricted to a single environment (soil+plant), making environmental stratification not applicable.

**Specific Findings Highlights**:
- **41% phylogenetic attenuation is substantial but well-contextualized**: The REPORT.md explains that this reflects higher phylogenetic clustering in the soil+plant subset compared to the pan-bacterial baseline (2% attenuation). The key point — that the association remains significant after correction (p=2.1×10⁻⁵) — is clearly stated.
- **Operon-distance and KEGG pathway tests refute genomic linkage**: The median 1,097-gene distance (vs null 350, Z=120.7) and the negative pathway co-membership (Z=-31.0) together confirm that the coupling is ecological (organism-level co-selection) rather than physical gene linkage or shared metabolism. This is a strong mechanistic finding.
- **Wang 2021 validation is non-significant in soil+plant subset (p=0.070)**: This is appropriately discussed as a power issue (smaller phylogenetic diversity in the subset) and a reflection that the phytase–siderophore axis may operate across a broader ecological range than the P-acquisition × metal-handling axis. The phylogenetic regression failure to converge is consistent with the weak signal.

## Suggestions

1. **Resolve Figure 3 inconsistency**: Either regenerate `figure3_env_stratification.png` for the v4 soil+plant subset (showing within-subset environmental variation if any exists) or move it to a `figures/archive_v3/` directory and remove it from the README table. The current state — present in `figures/` but marked as not-regenerated — may confuse readers about which figures belong to the v4 analysis.

2. **Document convergence verification**: Add a `data/phylo/convergence_check.txt` file containing the output of step 10a (`grep "FALSE" data/phylo/phylo_logistic.csv || echo "All models converged"`) with a timestamp. This makes the "all models converged" claim in REPORT.md directly verifiable from the repository.

3. **Justify phenazine operon threshold**: Add a brief note in REPORT.md or RESEARCH_PLAN.md explaining why ≥3 distinct phz gene families was chosen as the operon-carrier threshold. A sensitivity check showing that results are robust to thresholds of 2, 3, or 4 would further strengthen this design choice (though this is optional, not critical).

4. **Clarify notebook vs. script relationship in README**: The README `## Reproduction` section lists 16 scripts but also references 7 notebooks (NB01–NB07) without explicitly stating their relationship. Add a note like: "Notebooks (NB01–NB07) are inspection and display notebooks that load pre-computed CSVs from `data/`. They visualize and summarize results but do not perform primary analysis. The primary analysis code is in `src/*.py`." (This language is already present at the bottom of the Reproduction section, but moving it to the top would improve clarity.)

5. **Consider adding a "Quick Start" for non-Spark users**: The reproduction guide is excellent but requires Spark access for steps 1, 7, 11, 12, and 16. For readers without BERDL access, consider adding a note at the top: "All primary analysis results are pre-cached in `data/`. To explore findings without re-running Spark steps, start directly from step 2 or view notebooks NB01–NB07, which load pre-computed CSVs."

6. **Minor: Fix step 10a formatting**: The step 10a convergence check uses `||` to print "All models converged" if no failures are found. This is fine for interactive use, but for reproducibility, consider `grep "FALSE" data/phylo/phylo_logistic.csv > data/phylo/convergence_failures.txt || echo "All models converged" > data/phylo/convergence_failures.txt` so that the result is always committed.

## Review Metadata

- **Reviewer**: BERIL Automated Review (Claude, sonnet)
- **Date**: 2026-05-14
- **Scope**: README.md, RESEARCH_PLAN.md, REPORT.md (300 lines), 7 notebooks (all with saved outputs), 16 source scripts, 6 figures, data/ directory structure
- **Note**: This review was generated by an AI system. It should be treated as advisory input, not a definitive assessment.
