---
reviewer: BERIL Automated Review (Claude, claude-sonnet-4-6)
date: 2026-05-08
project: macro_micro_nutrient_cocycling
---

# Review: Macro-Micro Nutrient Gene Co-cycling in Bacterial Pangenomes

## Summary

This is a well-executed, hypothesis-driven pangenomics study that convincingly demonstrates
statistically significant co-occurrence of P-acquisition, N-fixation, and metal-handling gene
families across 27,682 GTDB species pangenomes from `kbase_ke_pangenome`. The research question
is clearly stated, the four-step analytical plan is logical, and the REPORT is unusually complete:
all major statistics are cited numerically, limitations are thoroughly enumerated (7 points), and
future directions are specific. All 7 notebooks carry saved outputs — no empty code-only files —
and both figures are present. The main areas for improvement are: one reproducibility gap (NB02
omits a permutation result that appears in the REPORT), a potentially misleading count in NB01's
summary print, fragile PFAM-merging code in src/01, and minor issues with the environmental
stratification table and the floor-bounded permutation p-values.

---

## Methodology

**Strengths.** The hypothesis (H1: macro-nutrient and metal-handling gene co-occurrence exceeds a
permutation null) is falsifiable and operationalized concretely. The four analytical layers —
group co-occurrence → core/accessory enrichment → phylogenetic stratification → environmental
stratification — each add a distinct explanatory dimension. The positive-control species check
(NB05) is a strong design choice that grounds the genome-scale signal in known biology.
Sensitivity analyses are rigorous: N-fixation is deliberately restricted to KEGG KO-defined
nifH/nifD (2,746 species) rather than the broader PF00142 Fer4 domain (5,872 species), and the
broader definition is correctly reported as a sensitivity check. The ≥3-phz-gene operon threshold
is explained and its limitations acknowledged.

**Gap — phylogenetic non-independence.** The permutation null shuffles group membership labels
uniformly, preserving marginal species counts but not phylogenetic structure. Closely related
species share gene content by descent; inflated co-occurrence statistics in species-rich clades
(e.g., Streptomycetaceae: 24 species, all phenazine operon carriers) are not controlled for. This
is honestly acknowledged as Limitation 4. It is the most significant methodological caveat and
warrants a more prominent statement in the Results narrative, not just the Limitations section.

**Minor.** The ≥3-phz-gene threshold for "operon carrier" is a heuristic. The REPORT notes this
(Limitation 6), but does not explain why 3 was chosen over 2 or 4. A brief sentence on the
sensitivity of the 63-species count to this threshold would strengthen the justification.

---

## Code Quality

**SQL and Spark (src/01_extract_gene_families.py).** Queries are structured correctly: KO and
gene-name conditions are aggregated in a single GROUP BY pass; PFAM conditions use versioned IDs
(`LIKE 'PFxxxxx.%'`) as recommended in BERDL conventions. The two-stage approach — main query for
KO/gene-name families, then separate per-PFAM queries — is reasonable given that correlated
subqueries inside a GROUP BY can be expensive at 132M-row scale.

**Issue — fragile PFAM merging loop (src/01, lines 98–104).** The conditional
`pdf = pdf.set_index(...) if 'gtdb_species_clade_id' in pdf.columns else pdf` is executed
inside the inner column loop, meaning: the first column of the first PFAM family converts `pdf`
to index-based access; subsequent columns skip the set_index; `reset_index()` at the end of each
outer iteration restores the column for the next PFAM family. This logic is correct as written,
but its correctness depends on execution order and column-name state that is easy to break with a
minor refactor. Preferred fix: set the index once before the outer PFAM loop and call
`reset_index()` once after all PFAM families are merged.

**Statistical methods (src/02_cooccurrence_stats.py).** Phi coefficient, Jaccard, Fisher's exact,
permutation null (N=1000), and BH-FDR are all appropriate choices and correctly implemented. The
manual BH step-down implementation (lines 211–219) produces results that match the REPORT values.
Using `statsmodels.stats.multitest.multipletests` would reduce the risk of off-by-one errors in
future modifications, but the current implementation is verified correct.

**Permutation floor effect (src/02, lines 156–157, and NB02 cell_6).** With N=1000 permutations,
the minimum representable p-value under the formula `(count+1)/(N+1)` is 1/1001 ≈ 0.001. All
three strong signals (P×Metal, N×Metal, P×N) report `perm_p=0.0010`, indicating the observed phi
exceeded all 1,000 null draws. The REPORT states these p-values as exact figures (e.g.,
"permutation Z=17.7, p<0.001") without noting the floor. Adding "p<0.001 (permutation floor,
N=1,000)" would be more precise.

**Stouffer meta-analysis (src/03_core_accessory_enrichment.py, lines 92–97).** The per-species
Fisher test is on a 2×2 table of `(nutrient_core, nutrient_noncore) × (metal_core,
metal_noncore)` within each co-encoding species. This correctly tests whether the core-genome
fraction differs between the two functional gene sets within a species, consistent with what the
REPORT claims (Signatures 1–3). The sign convention (z = -z when OR < 1) is correct.

**NifH_pfam definition in NB01 summary (NB01, cell_5).** The summary print reads:
```
N + Metal: 5718
```
This uses `df['has_N_fixation']` from the CSV, which is defined as the **broad** union of
nifH+nifD+nifH_pfam (5,872 species), not the KO-only 2,746. The primary analysis in src/02 uses
KO-only (2,719 co-occurrences with Metal). A reader scanning NB01 for headline counts could cite
5,718 as the N+Metal figure when the REPORT's primary number is 2,719. The note in NB01 cell_5
("used as primary definition in src/02") is present but easy to miss. Renaming the `has_N_fixation`
column in the CSV to `has_N_fixation_broad` and adding a separate `has_N_fixation_ko` column would
eliminate the ambiguity.

**Pitfall compliance.** The known pitfalls relevant to this project were checked:
- PFAM versioned IDs (`LIKE 'PFxxxxx.%'`): ✓ correctly applied for all four PFAM families.
- kbase_ke_pangenome schema usage: no evidence of strain-name collisions (GTDB species clade IDs
  are used throughout, never short strain names).
- Notebook artifacts: notebooks are committed with outputs. The scripts in `src/` pre-date the
  notebooks and serve as the primary reproducible record; the notebooks load from cached CSVs.
  This pattern is explicitly supported by the README and is not a pitfall violation.

---

## Reproducibility

**Notebooks and outputs.** All 7 notebooks (NB01–NB07) have saved cell outputs (text tables,
printed statistics). No empty code-only notebooks. This is a clear strength.

**Figures.** Both figures are present: `figures/figure1_cooccurrence.png/.pdf` (multi-panel,
635 KB) and `figures/forest_plot.png` (386 KB). The figures are referenced from REPORT.md with
panel descriptions.

**Dependencies.** `requirements.txt` lists numpy, pandas, scipy, matplotlib, pyspark, and
berdl_notebook_utils. Versions are specified (numpy≥1.24, pandas≥2.0, etc.). The comment
distinguishing Spark-dependent from local dependencies is helpful.

**Reproduction guide.** README.md has a `## Reproduction` section listing all 8 scripts, clearly
marking steps 1 and 7 as Spark-dependent and the rest as local. This is well done.

**Reproducibility gap — NB02 vs REPORT (Phz×P permutation).** `src/02_cooccurrence_stats.py`
runs permutation tests for **5** group pairs (lines 141–145), including
`('Phenazine_operon', 'P_acquisition')`. The REPORT's Table 1 cites the result: Phz-operon × P,
permutation Z=4.0. However, **NB02 cell_6** only runs permutation tests for **4** pairs — it
omits `('Phenazine_operon', 'P_acquisition')`. The Phz×P Z=4.0 result cited in the REPORT is
therefore not visible in any notebook output; it is only verifiable by running src/02 directly.
This is the most concrete reproducibility gap in the project.

**NB03 covers two analytical steps.** NB03 is titled "Core vs. Accessory Enrichment and
Phylogenetic Stratification," corresponding to both `src/03_core_accessory_enrichment.py`
(Step 3) and `src/04_phylogenetic_stratification.py` (Step 4). The README's 8-step reproduction
guide thus has no direct 1:1 mapping to the 7 notebooks (steps 3+4 → NB03). This is a minor
structural inconsistency but not a blocking issue — the README clearly lists both scripts.

**Environmental stratification "other" category.** NB06 outputs show 12,628 of 27,017 species
(46.8%) fall in "other" (unclassifiable isolations). The REPORT's environmental stratification
table (Section 6) presents only the 6 named environments without noting the "other" category.
Since the co-occurrence analysis in NB06 does include "other" in the printed output (log-OR=+1.0,
p=5.4×10⁻³⁹), adding a brief footnote would improve transparency.

---

## Findings Assessment

**Core claims are supported.** The key statistics in the REPORT are consistent with notebook
outputs and data files:
- P×Metal: phi=0.110, OR=2.30, p=1.3×10⁻⁶⁵ (NB02 cooccurrence_matrix: 1.284×10⁻⁶⁵ ✓)
- N×Metal: phi=0.088, OR=10.1, p=1.3×10⁻⁷¹ (NB02: 1.259×10⁻⁷¹ ✓)
- Phz-operon × Metal: p=9.4×10⁻³, OR=∞ (NB02: 9.369×10⁻³ ✓)
- 63/63 phenazine operon carriers encode metal genes (NB01: 63, NB05: confirmed)

**Positive control validation.** All 8 positive-control species encode P-acquisition and
metal-handling genes. The confirmed diazotrophs (*B. diazoefficiens*, *S. meliloti*) carry nifH
and nifD. *P. chlororaphis* is the sole positive control with a complete phenazine operon (≥3 phz
genes), consistent with its known phenotype. The *S. coelicolor* → *S. anthocyanicus* GTDB R214
reclassification is correctly handled and documented.

**Negative associations are informative and correctly interpreted.** The pstC/S × feoB depletion
(phi=−0.256, enrichment=0.67×) is correctly framed as niche separation between aerobic
P-scavengers and anaerobic/microaerobic Fe²⁺ importers. The N × Phz-operon depletion (0/63
diazotroph phenazine operon carriers) is correctly attributed to taxonomic concentration of
phenazine operons in non-diazotrophic Actinomycetota and Pseudomonas lineages.

**Interpretation is appropriately hedged.** The REPORT consistently distinguishes genomic
potential from expression ("encoding a gene does not guarantee its expression"). The *Xenorhabdus*
exception (insect pathogenesis, not mineral dissolution) is explicitly noted as a case where the
same pathway serves a different ecological function.

**One internal inconsistency in the environmental stratification table.** The REPORT Table (Section
6) shows "Soil/rhizosphere: n=3,406" but NB06 cell_2 prints `soil/rhizosphere: 3409` and cell_4
prints `n= 3406` for the P×Metal test. The discrepancy (3409 vs 3406) arises because 3 species
have environment assignments but were excluded from the co-occurrence test (incomplete gene family
data). This 3-species gap is not explained in the REPORT.

---

## Suggestions

1. **(Critical — Reproducibility)** Add the `('Phenazine_operon', 'P_acquisition')` pair to NB02
   cell_6's `key_pairs` list so that the Z=4.0 result cited in the REPORT's Table 1 appears in
   notebook output. Currently this result exists only in src/02 and cannot be verified from the
   notebooks alone.

2. **(High — Code clarity)** In src/01, rename `has_N_fixation` in the output CSV to
   `has_N_fixation_broad` and the existing `has_N_fixation_ko` to `has_N_fixation`. Then update
   NB01 cell_5 to use the correct column name. This eliminates the risk that NB01's "N+Metal:
   5718" (broad) is confused with the analysis-level "N+Metal: 2719" (KO-only).

3. **(High — Code robustness)** Refactor the PFAM merging loop in src/01 (lines 96–104). Set
   the index once before the outer `for name, pfam_prefix in pfam_families.items()` loop:
   ```python
   pdf = pdf.set_index('gtdb_species_clade_id')
   for name, pfam_prefix in pfam_families.items():
       ...
       for col in pfam_df.columns:
           if col not in pdf.columns:
               pdf[col] = 0
           pdf.loc[pfam_df.index, col] = pfam_df[col]
   pdf = pdf.reset_index()
   ```
   This avoids the fragile in-loop conditional set_index pattern.

4. **(Medium — Statistics transparency)** Note the permutation floor in the REPORT and NB02. For
   example: "permutation p<0.001 (floor at N=1,000; observed phi exceeded all null draws)" for
   the three strongest pairs. This is more precise than reporting "p=0.0010" as if it were an
   exact value.

5. **(Medium — Methods transparency)** Add a brief footnote to the REPORT's environmental
   stratification table noting that 12,628 species (46.8% of the 27,017 with environment data)
   fell in the "other" (unclassifiable) category. Including the "other" row's P×Metal log-OR in
   the table would give a complete picture.

6. **(Medium — Methods justification)** Add one sentence in the Methods or Results justifying why
   the phenazine operon threshold is ≥3 genes. For example: "We confirmed that the 63-species
   count is robust to the threshold (≥2 genes: 312 species; ≥4 genes: 48 species), and chose ≥3
   as the minimum that filters broad-family phzF-only hits while capturing known Pseudomonas operon
   architectures." (Note: verify these alternate counts from the data before publishing.)

7. **(Low — REPORT clarity)** Section 6 (environmental stratification) notes "n species" with a
   3-species discrepancy for soil/rhizosphere (3,409 in the species count vs. 3,406 in the test
   table). Add a brief note (e.g., "3 species excluded from co-occurrence test due to incomplete
   gene family annotations") to prevent reader confusion.

8. **(Low — Notebook structure)** Consider splitting NB03 into NB03_core_accessory and
   NB04_phylogenetic (renumbering subsequent notebooks). Currently, step 4 (phylogenetic
   stratification from src/04) is embedded in NB03 without a section header that aligns to the
   step numbering in the README reproduction guide.

9. **(Low — Future-proofing)** Replace the manual BH correction in src/02 (lines 211–219) with
   `from statsmodels.stats.multitest import multipletests` (already available since statsmodels is
   a scipy dependency). The manual implementation is correct for this dataset but is a maintenance
   risk.

---

## Review Metadata

- **Reviewer**: BERIL Automated Review (Claude, claude-sonnet-4-6)
- **Date**: 2026-05-08
- **Scope**: README.md, RESEARCH_PLAN.md, REPORT.md, 7 notebooks (NB01–NB07), 5 src/ scripts
  reviewed (src/01–src/03, src/07), 17 data files, 2 figures (figure1_cooccurrence.png/.pdf,
  forest_plot.png), requirements.txt, beril.yaml, docs/pitfalls.md
- **Note**: This review was generated by an AI system. It should be treated as advisory input, not a definitive assessment.
