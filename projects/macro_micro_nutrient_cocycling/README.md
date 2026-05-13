# Macro-Micro Nutrient Gene Co-cycling in Bacterial Pangenomes

## Research Question

Do bacterial pangenomes encoding macro-nutrient acquisition (P, N) and phenazine biosynthesis genes disproportionately co-encode trace-metal handling genes, reflecting the biogeochemical coupling of nutrient and metal cycles through enzyme cofactor requirements and Fe-oxyhydroxide mineral dissolution?

## Status

Completed — Significant co-occurrence of macro-nutrient acquisition and metal-handling genes across 27,682 bacterial pangenomes, with all 63 phenazine operon carriers encoding both P-acquisition and metal-handling genes.

## Overview

Across 27,682 GTDB species pangenomes from `kbase_ke_pangenome`, we find significant co-occurrence of P-acquisition and metal-handling genes (phi=0.110, OR=2.3, p=1.3×10⁻⁶⁵), N-fixation and metal-handling genes (phi=0.088, OR=10.1, p=1.3×10⁻⁷¹), and complete overlap of phenazine operon carriers with metal-handling genes (63/63 species, 100%). The 63 phenazine operon carriers are concentrated in soil Actinomycetota and rhizosphere Pseudomonadota, consistent with the McRose-Newman model of phenazine-mediated Fe(III)-oxyhydroxide dissolution as a joint P and trace-metal mobilization strategy.

## Quick Links

- [RESEARCH_PLAN.md](RESEARCH_PLAN.md) — Hypothesis, approach, revision history
- [REPORT.md](REPORT.md) — Full findings and interpretation
- [figures/figure1_cooccurrence.png](figures/figure1_cooccurrence.png) — Multi-panel figure

## Reproduction

```bash
# Step 1: Extract gene families (requires on-cluster Spark access)
python src/01_extract_gene_families.py

# Step 2: Co-occurrence statistics
python src/02_cooccurrence_stats.py

# Step 3: Core vs. accessory enrichment
python src/03_core_accessory_enrichment.py

# Step 4: Phylogenetic stratification
python src/04_phylogenetic_stratification.py

# Step 5: Figure 1 (multi-panel)
python src/05_figure.py

# Step 6: Positive-control species check
python src/06_positive_controls.py

# Step 7: Environmental stratification (requires Spark for ncbi_env query)
python src/07_environmental_stratification.py

# Step 8: Forest plot (Figure 2)
python src/08_forest_plot.py

# Step 9: Phylogenetic signal — Pagel's lambda (requires R env)
conda run -n r_phylo Rscript src/09_phylo_signal.R \
  data/phylo/bac120_r214_pruned.tree \
  data/species_gene_families.csv \
  data/phylo/phylo_signal.csv

# Step 10: Phylogenetic logistic regression (requires R env)
conda run -n r_phylo Rscript src/10_phylo_logistic.R \
  data/phylo/bac120_r214_pruned.tree \
  data/species_gene_families.csv \
  data/phylo/pair_definitions.csv \
  data/phylo/phylo_logistic.csv

# Step 11: Operon-distance test (requires on-cluster Spark access)
python src/11_operon_distance.py

# Step 12: Wang 2021 phytase-siderophore validation (requires on-cluster Spark access)
python src/12_wang2021_validation.py

# Step 13: Figure 4 — phylogenetic correction scatter
python src/13_figure4_phylo_correction.py

# Step 14: Figure 5 — operon distance distribution
python src/14_figure5_operon_distance.py

# Step 15: Figure 6 — Wang 2021 family-level co-occurrence
python src/15_figure6_wang2021.py

# Step 16: KEGG pathway co-membership test (requires on-cluster Spark access)
python src/16_kegg_pathway_comembership.py
```

Steps 1, 7, 11, and 12 require access to the `kbase_ke_pangenome` tenant on BERDL. Steps 9–10 require an R environment with `phylolm` and `phytools` (create via `conda create -n r_phylo -c conda-forge r-base r-ape r-phytools r-phylolm --solver=classic`). All other steps run locally on CSV outputs in `data/`.

**Notebooks (NB01–NB07)** are inspection and display notebooks that load pre-computed CSVs from `data/`. They visualize and summarize results but do not perform primary analysis. The primary analysis code is in `src/*.py`. Re-running a notebook does not re-query BERDL — it renders cached CSV outputs.

## Authors

- Jing Tao (jingtao-lbl), Lawrence Berkeley National Laboratory
