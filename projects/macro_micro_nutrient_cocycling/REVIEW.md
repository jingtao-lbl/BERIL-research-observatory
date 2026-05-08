---
reviewer: BERIL Automated Review (Claude, claude-sonnet-4-6)
date: 2026-05-08
project: macro_micro_nutrient_cocycling
---

# Review: Macro-Micro Nutrient Gene Co-cycling in Bacterial Pangenomes

## Summary

This is a well-executed, large-scale pangenome analysis testing the genomic coupling of
macro-nutrient acquisition and trace-metal handling genes across 27,682 GTDB species
pangenomes. The research question is clearly grounded in biochemical theory (McRose &
Newman 2021), the statistical approach is rigorous (Fisher's exact, BH FDR, permutation
tests, Stouffer meta-analysis), and the limitations section is unusually honest and
detailed — nine distinct caveats, including the PF00142 overcount issue, phylogenetic
non-independence, and annotation case-sensitivity. Positive-control validation against
eight model organisms, environmental stratification across six habitat categories, and a
per-phylum forest plot together give the findings strong multi-angle support. The primary
area for improvement is the non-standard separation of analysis code into `src/*.py`
scripts while relegating notebooks to display-only roles: the analytical logic is
committed and reproducible, but the notebook-centric verification workflow expected by
BERDL reviewers requires reading Python scripts rather than notebook cells.

## Methodology

**Research question:** Clearly stated and testable: do bacterial pangenomes with P- and
N-acquisition genes disproportionately co-encode trace-metal handling genes, and does
this coupling intensify in phenazine-producing lineages? The formal H1/H0 pair in
RESEARCH_PLAN.md is precise.

**Approach:** The four-step pipeline (gene extraction → co-occurrence statistics →
core/accessory enrichment → phylogenetic stratification) is well-matched to the
question. The additional enrichments (positive controls, environmental stratification,
forest plot) strengthen the findings substantially. The distinction between the Phz
"broad" group (any phz gene ≥ 1; 10,308 species) and the "operon carrier" group (≥ 3
distinct phz families; 63 species) is well-motivated and honestly interpreted: the ≥ 3
threshold is explicitly flagged as a heuristic, and the sensitivity to threshold choice
is reported (≥ 2 → 1,125 species; ≥ 4 → 34 species).

**Data sources:** Clearly identified — `kbase_ke_pangenome`, specific tables
(`gene_cluster`, `bakta_annotations`, `bakta_pfam_domains`, `gtdb_species_clade`,
`ncbi_env`), and database version (GTDB R214, Bakta annotations).

**Reproducibility:** The README `## Reproduction` section cleanly separates Spark-dependent
steps (1 and 7) from locally-runnable steps (2–6 and 8). The CSV outputs in `data/` are
committed, enabling local re-analysis without Spark access. The notebook/script
relationship is explicitly explained.

**One methodology concern — environmental stratification multiple testing:** The
environmental stratification analysis (src/07, NB06) tests 7 environments × 3 pairs = 21
Fisher tests but does not apply BH FDR correction. The REPORT correctly presents raw
p-values, but the table could mislead readers who apply a conventional α=0.05 threshold
to uncorrected tests. Adding a footnote or a corrected q-value column to Table S6 (or
the equivalent data file) would remove ambiguity.

## Code Quality

**SQL correctness:** The gene family extraction query in `src/01_extract_gene_families.py`
is correctly structured. PFAM-based families (phoD_pfam, nifH_pfam, feoB_pfam, HMA_pfam)
are extracted in separate queries to avoid correlated-subquery complexity, then merged
into the main data frame via indexed join. The two-pass approach is sound, and the
subsequent `fillna(0)` correctly handles species absent from a given PFAM query.

**Statistical implementations:**
- The BH FDR implementation in `src/02_cooccurrence_stats.py` (lines 210–218) is
  correct: sort p-values ascending, compute rank-adjusted q, enforce monotonicity
  downward, clip to [0, 1].
- The phi coefficient formula and permutation test (1,000 shuffles, seed=42, two-tailed
  |phi| comparison, Laplace-corrected p = (count ≥ |obs| + 1) / (N + 1)) are standard
  and correctly implemented.
- The Stouffer meta-analysis in `src/03_core_accessory_enrichment.py` aggregates per-
  species Fisher Z-scores. Outputs in `data/core_enrichment_summary.csv` are internally
  consistent with NB03 cell outputs.

**Gene name case sensitivity (Limitation 9 — acknowledged):** `src/01` uses exact
string matching for gene-name-based families (`ba.gene = 'pstA'`, etc.). Bakta
annotation casing is not standardized. Affected families (pstA/B/C/S, corA, phnC/D/E,
phzA/B/D/G/M) may have undercounts. KEGG KO and PFAM families are unaffected. The
practical impact is unknown without a casing audit, but a simple `LOWER(ba.gene) =
'psta'` pattern in the WHERE clause would close this gap at negligible query cost. This
is worth fixing in a follow-up run.

**Notebook organization:** Each notebook follows a clear load → display → summarize
structure. All seven notebooks have saved cell outputs (text tables, printed statistics),
enabling offline review without re-executing. NB04 saves figure1 in both PNG and PDF
formats. NB06 and NB07 also save figures and print summary tables. No empty-output
notebooks observed.

**Pitfall compliance:** The project does not trigger the main BERDL pitfalls:
- No ENIGMA strain-name matching (short-name collision pitfall not applicable).
- No Web of Microbes binary encoding.
- No cross-branch concurrent sessions evident.
- The `kbase_ke_pangenome` tenant is used correctly (tenant prefix is `kbase`).

## Findings Assessment

**Are conclusions supported?** Yes. The three headline findings (P×Metal phi=0.110
OR=2.3 p<10⁻⁶⁴; N×Metal OR=10.1 p<10⁻⁷⁰; Phz-operon 63/63 with metal genes) are
reproducible from the committed CSV data and are internally consistent across the
notebook outputs, the REPORT tables, and the figures. The permutation tests (Z > 2.3 for
all three main pairs) confirm the Fisher results are not artifacts of marginal frequency
imbalance.

**Limitations acknowledged?** Extensively. The nine-point limitations section covers:
presence/absence-only scoring (no expression data), annotation-dependent coverage, PF00142
inflation of N-fixation, phylogenetic non-independence, ecological inference from genomic
potential, phenazine threshold heuristics, actinomycete phenazine under-detection,
environmental classification keyword priority, and gene name case sensitivity. This is
one of the more thorough limitations sections in recent BERDL projects.

**One interpretive nuance to sharpen:** The REPORT (Results §3) presents the Phz_genes
core enrichment result (Stouffer Z=2.636, `core_enrichment_summary.csv` row for
`Phz_genes`) alongside P-genes (Z=68.3) and N-genes (Z=3.249). The `frac_enriched=0.268`
and `median_OR=0.333` for Phz_genes indicate that, at the per-species level, phenazine
genes are *not* preferentially in the core relative to metal genes in most co-encoding
species. This is worth a sentence in the Interpretation section — it distinguishes the
Phz case from P-acquisition (core-enriched) and is consistent with the accessory/HGT
framing of phenazine gene distribution.

**One count discrepancy to reconcile:** NB06 reports 27,017 species with environment
assignments; summing the per-environment rows yields 27,013 (3,406 + 1,134 + 3,449 +
2,059 + 3,497 + 840 + 12,628). The REPORT's text mentions "27,009" in one sentence. All
three numbers are close (within 8 of each other) and the REPORT explains that 8 species
with environment assignments lacked complete gene family annotations. The discrepancy
across these three numbers (27,017 / 27,013 / 27,009) should be unified in the REPORT
text — likely they reflect different filtering steps applied in different places.

**Figures:** Three figures are present and cover the key results:
- `figures/figure1_cooccurrence.png` — four-panel: phi heatmap, core fraction bars,
  phylum stratification, phenazine taxonomy. Clear axis labels, color coding by gene
  group, phi values annotated in each cell.
- `figures/forest_plot.png` — per-phylum log-OR for three pairs across 34 phyla.
- `figures/figure3_env_stratification.png` — environmental stratification bar chart.

Minor figure concern: Figure 2 (forest_plot.png) does not annotate sample sizes on the
x-axis bars. Phyla with very wide CIs (e.g., Deinococcota, P×Metal: CI [-0.16, +9.46])
are visually prominent but non-significant; the phylum name alone with `(n=N)` in the
y-axis label (as done for the NB03 display) would help readers calibrate quickly.

## Suggestions

1. **Apply BH FDR to environmental stratification tests (21 tests across 7 environments
   × 3 pairs).** Add a `q_value` column to `data/env_cooccurrence.csv` and update
   Table 6 in the REPORT to flag which environmental effects survive correction. Several
   of the nominally significant N×Metal effects (e.g., plant-associated p=0.064) may
   not survive; this won't change the main conclusion but will make the interpretation
   more defensible.

2. **Fix gene name case sensitivity in `src/01`.** Replace exact-match `ba.gene = 'pstA'`
   with `LOWER(ba.gene) = 'psta'` (and similarly for all other gene-name-based families).
   Re-run the full pipeline and report any changes to species counts. This is a one-line
   fix per gene family and would close Limitation 9.

3. **Add a sentence to the Interpretation section on Phz core fractions.** The `frac_enriched=0.268`
   and `median_OR=0.333` for Phz_genes (from `core_enrichment_summary.csv`) indicate
   that phenazine genes are not preferentially core relative to metal genes in most
   co-encoding species. This completes the three-signature story (P = core coupling,
   N = horizontally mobile coupling, Phz = ecologically specialized but without strong
   core enrichment) and removes a potential reader confusion when comparing the Stouffer
   Z=2.636 result to the stronger P-acquisition Z=68.3.

4. **Reconcile the three species counts in the environmental analysis** (27,017 / 27,013 /
   27,009). Choose one number as the authoritative total-with-environment count, add a
   comment in `src/07` explaining the filtering steps that reduce from that number to the
   per-environment co-occurrence test counts, and update the REPORT text to be consistent.

5. **Add sample-size annotations to Figure 2 (forest plot).** Include n= in each phylum's
   y-axis label (already done in NB07's printed output; propagate to the saved PNG). This
   is a low-effort change to `src/08_forest_plot.py` that materially improves the figure's
   standalone readability.

6. **Consider moving core analysis logic into notebooks, or annotating NB cells with
   the corresponding src/ output.** The current workflow — analysis in `src/*.py`, display
   in `NB*.ipynb` — is fully documented and reproducible, but differs from the standard
   BERDL reproducibility expectation that a reviewer can follow the analytical logic
   through notebook cells. At minimum, adding `# see src/02_cooccurrence_stats.py
   lines 50–65` cross-references in the NB cells that load pre-computed CSVs would help
   reviewers navigate between the display and the analytical source.

7. **(Future) Apply phylogenetic logistic regression** (as noted in Future Directions) to
   test co-occurrence after controlling for phylogenetic signal. The permutation test
   controls for marginal frequencies but not for shared ancestry — this is Limitation 4
   and the most important gap for the N-fixation and P-acquisition findings, where deep
   clades (Cyanobacteriota, Bacillota, Rhizobiales) may dominate the signal.

## Review Metadata

- **Reviewer**: BERIL Automated Review (Claude, claude-sonnet-4-6)
- **Date**: 2026-05-08
- **Scope**: README.md, RESEARCH_PLAN.md, REPORT.md, beril.yaml, requirements.txt,
  7 notebooks (NB01–NB07, all with saved outputs), 2 source scripts examined in full
  (src/01, src/02), 2 source scripts examined partially (src/07, src/08),
  8 data files noted, 3 figures present
- **Note**: This review was generated by an AI system. It should be treated as advisory
  input, not a definitive assessment.
