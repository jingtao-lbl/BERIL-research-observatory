---
reviewer: BERIL Automated Review (Claude, claude-sonnet-4-6)
date: 2026-05-07
project: macro_micro_nutrient_cocycling
---

# Review: Macro-Micro Nutrient Gene Co-cycling in Bacterial Pangenomes

## Summary

This is a well-executed, largely complete pangenomics project that tests a mechanistically motivated hypothesis: do bacterial genomes encoding macro-nutrient acquisition (P, N) and phenazine biosynthesis genes disproportionately co-encode trace-metal handling genes? The research question is specific and grounded in published biochemistry (McRose & Newman 2021), the four-step analysis is methodologically sound, and the findings are clearly supported by the data. The project is reproducible from cached CSV outputs (Steps 2–5), documentation is thorough, and limitations are honestly acknowledged — including the dominant caveat of phylogenetic non-independence. The main areas for improvement are: (1) a potentially over-broad PFAM definition for `nifH_pfam` (PF00142) that may be inflating the N-fixation × metal signal; (2) absence of multiple testing correction discussion for the 76 pairwise Fisher tests; and (3) an unreported negative Stouffer Z for the N-fixation core-enrichment analysis that has biological interest. Overall the project is close to submission-ready with minor revisions.

## Methodology

**Research question and hypothesis:** Both are clearly stated and testable. The null and alternative hypotheses (H0/H1) in RESEARCH_PLAN.md are specific and matched to the statistical tests used. The link to the McRose-Newman mechanism and cofactor biochemistry provides a compelling prior expectation for each tested pair.

**Approach:** The four-step pipeline — gene family extraction → co-occurrence statistics → core/accessory enrichment → phylogenetic stratification — is logical and well-scoped. The permutation null (1,000 shuffles preserving marginal counts) is an appropriate complement to Fisher's exact test, adding robustness to the inference. The Stouffer meta-analysis for core-enrichment is well-suited to the many-species design.

**Data sources:** The primary source (`kbase_ke_pangenome`, 132.5M gene clusters, 27,682 species) is clearly identified, along with the relevant tables (`gene_cluster`, `bakta_annotations`, `bakta_pfam_domains`, `gtdb_species_clade`). Gene family definitions are tabulated in the REPORT with annotation sources, making the query strategy auditable.

**Phenazine operon threshold:** The decision to require ≥3 distinct phz gene families to define a true "operon carrier" is well-motivated — without this filter, phzF alone (a broad superfamily) would score 10,856 species, reducing a specific biological signal to noise. This threshold is clearly documented and the 63-species result is a qualitatively distinct finding.

**Gap — PF00142 scope for nifH_pfam:** The RESEARCH_PLAN and REPORT identify nifH(Pfam) as PF00142. In PFAM, PF00142 is the "Ferredoxin" domain — a small [2Fe-2S]-binding domain found broadly in ferredoxin-like proteins, not exclusively in nitrogenase iron protein (NifH). The data bear this out: nifH by KO (K02588) covers 3,396 species, while nifH_pfam (PF00142) covers 6,632 species — nearly double. This strongly suggests PF00142 is capturing ferredoxin-like proteins beyond diazotrophs. Since the strongest individual co-occurrence signals in the analysis are nifH_pfam × feoB (phi=0.269) and nifH_pfam × HMA (phi=0.229), and these are interpreted as evidence for the Fe-Mo cofactor requirement of nitrogenase, it is important to assess how much of this signal arises from general ferredoxin-fold proteins (which also bind Fe-S clusters and would co-occur with metal-handling genes for unrelated reasons). The REPORT should either verify that PF00142 is specific to NifH in this dataset or add a caveat that the nifH_pfam results may capture ferredoxin-like proteins broadly.

**Reproducibility:** The README has a clear `## Reproduction` section that explicitly marks which step requires Spark access (Step 1 only) and which steps run locally on cached CSVs (Steps 2–5). This is excellent practice and follows the Spark/local separation guidance. Anyone with the cached data can run the full statistical analysis and figure generation without cluster access.

## Code Quality

**SQL queries (src/01_extract_gene_families.py):** Queries are clear and correctly structured. Gene name lookups use exact matching (`ba.gene = 'pstA'`), KEGG KOs use equality, and PFAM lookups use `LIKE 'PFxxxxx.%'` with version-suffix wildcarding — the recommended practice per the BERDL schema guides. The core/auxiliary/singleton breakdown per gene family is computed within the same query, which is efficient.

**Known pitfall compliance:** The PFAM LIKE pattern with trailing `.%` is correct. The project correctly avoids the "Commit Notebooks Alongside Artifacts" pitfall — notebooks NB01–NB04 are present with saved outputs (text, tables, figure). The pitfalls in docs/pitfalls.md concerning short-name collisions, WebofMicrobes binary encoding, and MetaPhlAn cross-cohort issues are not applicable to this project's data sources.

**Fragile index management in src/01 (minor bug risk):** Lines 100–103 in `src/01_extract_gene_families.py` call `pdf.set_index('gtdb_species_clade_id')` inside a column-level `for col` loop over PFAM results. The re-index only triggers on the first column (when `'gtdb_species_clade_id' in pdf.columns`), then silently uses the already-indexed DataFrame for subsequent columns. This pattern works for the current 4 PFAM families but is fragile — a future contributor who adds a PFAM family could introduce a hard-to-diagnose bug. The fix is simple: move the `set_index` call once outside the inner column loop.

```python
# Current (fragile — re-index embedded in column loop):
for col in pfam_df.columns:
    pdf = pdf.set_index('gtdb_species_clade_id') if 'gtdb_species_clade_id' in pdf.columns else pdf
    ...

# Recommended:
pdf = pdf.set_index('gtdb_species_clade_id')
for col in pfam_df.columns:
    ...
pdf = pdf.reset_index()
```

**Multiple testing correction absent:** The pairwise analysis runs 76 Fisher's exact tests (19 nutrient genes × 4 metal genes) with no Bonferroni or FDR correction. The Bonferroni threshold at α=0.05 is p < 0.00066. Nearly all reported "strong" results pass this threshold comfortably (e.g., nifH_pfam × feoB: p≈0, pstC × feoB: p≈0). However, some entries in the top-10 enrichment table in NB02 do not — notably phzM × HMA (n_both=6, p=0.13) appears in that table (sorted by enrichment ratio, not p-value) but is not statistically significant. The REPORT correctly omits this pair from its narrative, but the methods section should acknowledge the multiple testing context and state a corrected threshold.

**Statistical implementation:** phi coefficient, Jaccard, Fisher's exact, Stouffer meta-analysis, and permutation test are all correctly implemented. The permutation uses `np.random.seed(42)` for reproducibility. The Stouffer Z-score derivation from per-species Fisher p-values (converting p → Z, signing by OR direction) is standard and appropriate.

**Code organization:** Each src/ script produces clearly named output files and has a docstring listing outputs. Notebooks load pre-computed CSVs and display results cleanly. NB03 combines two src scripts (03 and 04) into one notebook — this is noted in NB03's header but creates a slight mismatch between the 5 src scripts and 4 notebooks that could confuse a new reader following the README reproduction steps.

## Findings Assessment

**Group-level co-occurrence:** The table in REPORT §Results.1 matches `cooccurrence_matrix.csv` and NB02 outputs. P × Metal (phi=0.110, OR=2.3, p=1.3×10⁻⁶⁵), N × Metal (phi=0.107, OR=4.04, p=1.5×10⁻⁸⁷), and Phz-operon × Metal (phi=0.014, OR=∞, p=9.4×10⁻³) are computed correctly and the conclusions are supported. The negative N × Phz-operon association (OR=0.19, n_both=3) is appropriately flagged as based on sparse counts.

**Individual gene pairs:** The reported top positive associations (nifH_pfam × HMA, nifH_pfam × feoB, phzG × corA, phoD × HMA, phzB × corA) and negative associations (pstC × feoB, pstS × feoB, pstC × HMA) are verified in `pairwise_detail.csv` and NB02. The biological interpretation of the pstC/S × feoB anti-correlation (redox geochemistry, aerobic vs. anaerobic niches) is plausible and well-argued.

**Core vs. accessory enrichment — undiscussed negative Stouffer Z:** The `core_enrichment_summary.csv` shows:

| Group | Stouffer Z | Direction | Discussed in REPORT? |
|-------|-----------|-----------|----------------------|
| P_genes | +68.3 (p≈0) | P-genes more core than metal genes | Yes |
| N_genes | **−4.853** (p=1×10⁻⁶) | N-genes **less** core than metal genes | **No** |
| Phz_genes | +2.636 (p=0.008) | Phz-genes more core | Briefly |

The negative Z for N_genes is an interpretable and interesting biological result: in co-occurring N+Metal species, nitrogenase genes are significantly more accessory than metal-handling genes. This is consistent with N-fixation being a conditionally acquired capacity (often on mobile genetic elements — nifH_pfam shows only 32.4% core, the lowest of any displayed gene). The REPORT's "two modes of coupling" framework (stable core vs. flexible accessory) currently places corA and pst/phn in the "stable" category and feoB/HMA/phoD in the "flexible" category, but nifH should also be in the flexible category. This result deserves a paragraph.

**Phenazine operon taxonomy — summary framing:** The REPORT states carriers are "dominated by soil Actinomycetota (35/63, 56%) and rhizosphere Pseudomonadota (27/63, 43%)." The 6 Enterobacteriaceae species are Xenorhabdus — insect-associated entomopathogenic nematode symbionts, not rhizosphere colonizers. These are within Pseudomonadota and would be included in the "27/63" figure. The detailed taxonomic section correctly identifies Xenorhabdus as "insect-associated," but the summary framing overstates the rhizosphere connection.

**Plant-associated lineage ceiling effects:** The REPORT correctly identifies that Pseudomonadaceae, Rhizobiaceae, Burkholderiaceae, Streptomycetaceae, and Xanthomonadaceae show near-universal P-acquisition and metal-handling prevalence (100%/93–100%), with enrichment ratios at ~1.00× due to ceiling effects. This is an honest acknowledgment that the phylogenetic stratification test is underpowered for families where virtually every species carries both gene sets.

**Conclusions supported:** The central conclusion — significant genomic coupling of macro-nutrient and metal-handling gene families across 27,682 bacterial pangenomes, with 100% overlap in phenazine operon carriers — is supported by the data. The mechanistic interpretations (Fe-Mo cofactor for nitrogenase, Fe-oxyhydroxide dissolution for phenazines) are appropriate inferences, not overstatements. Limitations are thoroughly and honestly enumerated.

## Suggestions

1. **Verify PF00142 specificity for nifH in this dataset (high priority).** Cross-check: what fraction of the 6,632 nifH_pfam species also have nifH by KO (K02588, 3,396 species)? If the overlap is under ~60%, PF00142 is capturing a substantial non-diazotrophic signal. If so, add a caveat to REPORT §Methods and §Results that the nifH_pfam × metal associations (phi=0.23–0.27) reflect ferredoxin-like proteins broadly, not just diazotrophs, and report the KO-based nifH × metal associations as the primary N-fixation signal.

2. **Discuss the N-fixation Stouffer Z (medium priority).** Add a paragraph in REPORT §Core vs. accessory structure noting that N_genes has Stouffer Z=−4.853 (p=1×10⁻⁶), indicating nitrogenase genes are significantly *more accessory* than co-occurring metal-handling genes in species that encode both. Interpret this as consistent with HGT-mediated acquisition of diazotrophy and update the "two modes of coupling" framework to place nifH in the "Flexible coupling" category alongside feoB/HMA/phoD.

3. **Add a multiple testing note to the Methods section (medium priority).** In REPORT §Statistical tests, state the number of individual-pair tests (76) and the Bonferroni threshold (p < 0.00066). Confirm that all narrative-featured associations pass this threshold. Note that the top-10 enrichment table in NB02 is sorted by enrichment ratio rather than p-value, so some low-n entries (e.g., phzM × HMA, p=0.13, n=6) are not statistically significant and should not be cited as evidence.

4. **Correct the summary framing of phenazine carrier ecology (low priority).** Revise "rhizosphere Pseudomonadota (27/63, 43%)" to specify that of the 26 Pseudomonadota carriers, 17 are rhizosphere Pseudomonadaceae, 6 are insect-associated Enterobacteriaceae (Xenorhabdus), 2 are soil Xanthomonadaceae, and 1 is soil Burkholderiaceae. The McRose-Newman model applies to the rhizosphere/soil fraction, not to the Xenorhabdus lineage.

5. **Fix the index management pattern in src/01 (low priority, robustness).** Move `pdf.set_index('gtdb_species_clade_id')` outside the inner `for col` loop (see Code Quality section above). This prevents a latent bug if the PFAM family list is extended in future work.

6. **Update beril.yaml status to "completed" (low priority).** `beril.yaml` shows `status: analysis` while README says `Status: Completed`. These should be consistent.

7. **Consider adding class-level stratification to the REPORT (nice-to-have).** `data/class_cooccurrence.csv` is generated by src/04 but not discussed in the REPORT. Within Pseudomonadota, phi values by class (e.g., Gammaproteobacteria vs. Alphaproteobacteria) would clarify which subclade drives the phylum-level signal and whether the Xenorhabdus-hosting Gammaproteobacteria differ from the Pseudomonas-hosting lineages.

## Review Metadata

- **Reviewer**: BERIL Automated Review (Claude, claude-sonnet-4-6)
- **Date**: 2026-05-07
- **Scope**: README.md, RESEARCH_PLAN.md, REPORT.md, 4 notebooks (NB01–NB04, all with saved outputs), 5 src/ scripts, 10 data files, 2 figures (figure1_cooccurrence.png/.pdf), requirements.txt, beril.yaml, docs/pitfalls.md
- **Note**: This review was generated by an AI system. It should be treated as advisory input, not a definitive assessment.
