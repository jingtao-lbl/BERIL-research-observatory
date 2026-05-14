# Macro-Micro Nutrient Gene Co-cycling in Bacterial Pangenomes

## Research Question

Do bacterial pangenomes encoding macro-nutrient acquisition (P, N) and phenazine biosynthesis genes disproportionately co-encode trace-metal handling genes, reflecting the biogeochemical coupling of nutrient and metal cycles through enzyme cofactor requirements and Fe-oxyhydroxide mineral dissolution?

## Status

Completed — Significant co-occurrence of macro-nutrient acquisition and metal-handling genes across 4,540 soil/plant-associated bacterial pangenomes (phi=0.129, OR=3.9), with the association surviving phylogenetic correction on a 4,177-tip subtree.

## Overview

This is a soil/plant-focused subset analysis of the pan-bacterial macro-micro nutrient co-cycling study (v3, 27,682 species). Restricting to 4,540 soil/rhizosphere and plant-associated GTDB species pangenomes from `kbase_ke_pangenome`, we find stronger P-acquisition × metal-handling co-occurrence (phi=0.129, OR=3.9, p=4.4×10⁻¹⁴) compared to the pan-bacterial baseline (phi=0.110, OR=2.3). The signal survives phylogenetic correction (phylo log-OR=0.938, p=2.1×10⁻⁵) though is attenuated by 41%, indicating substantial phylogenetic confounding in the soil+plant subset. All 27 phenazine operon carriers in the subset encode both P-acquisition and metal-handling genes, consistent with the 96.4% metal-handling prevalence in this subset.

The coupling is ecological, not genomic linkage or pathway integration: P × M gene pairs on the same contig are farther apart than expected by chance (median 1,097 genes vs null 350, Z=120.7), and P-acquisition and metal-handling genes share fewer KEGG pathways than expected (Z=−31.0). Co-occurrence reflects organism-level ecological co-selection, not physical linkage or shared metabolism.

The coupling is specific to alkaline phosphatases and high-affinity phosphate transporters: an independent phytase × siderophore test (Wang et al. 2021) is non-significant in the soil+plant subset (OR=1.23, p=0.070), suggesting the P × M signal is not universal across all P-mobilization strategies.

## Quick Links

- [RESEARCH_PLAN.md](RESEARCH_PLAN.md) — Hypothesis, approach, revision history
- [REPORT.md](REPORT.md) — Full findings and interpretation
- [figures/figure1_cooccurrence.png](figures/figure1_cooccurrence.png) — Multi-panel figure

## Reproduction

```bash
# Step 1: Extract gene families (requires on-cluster Spark access, ~15 min)
python src/01_extract_gene_families.py

# Step 2: Co-occurrence statistics (~2 min)
python src/02_cooccurrence_stats.py

# Step 3: Core vs. accessory enrichment (~2 min)
python src/03_core_accessory_enrichment.py

# Step 4: Phylogenetic stratification (~1 min)
python src/04_phylogenetic_stratification.py

# Step 5: Figure 1 (multi-panel, ~30 sec)
python src/05_figure.py

# Step 6: Positive-control species check (~10 sec)
python src/06_positive_controls.py

# Step 7: Environmental stratification (requires Spark, ~10 min)
python src/07_environmental_stratification.py

# Step 8: Forest plot — Figure 2 (~30 sec)
python src/08_forest_plot.py

# Step 9: Phylogenetic signal — Pagel's lambda (requires R env, ~5 min)
conda run -n r_phylo Rscript src/09_phylo_signal.R \
  data/phylo/bac120_r214_pruned.tree \
  data/species_gene_families.csv \
  data/species_taxonomy.csv \
  data/phylo/phylo_signal.csv

# Step 10: Phylogenetic logistic regression (requires R env, ~3 min)
conda run -n r_phylo Rscript src/10_phylo_logistic.R \
  data/phylo/bac120_r214_pruned.tree \
  data/species_gene_families.csv \
  data/phylo/pair_definitions.csv \
  data/phylo/phylo_logistic.csv

# Step 10a: Check for convergence failures
grep "FALSE" data/phylo/phylo_logistic.csv || echo "All models converged"

# Step 11: Operon-distance test (requires on-cluster Spark access, ~25 min)
python src/11_operon_distance.py

# Step 12: Wang 2021 phytase-siderophore validation (requires Spark, ~15 min)
python src/12_wang2021_validation.py

# Step 13: Figure 4 — phylogenetic correction scatter (~10 sec)
python src/13_figure4_phylo_correction.py

# Step 14: Figure 5 — operon distance distribution (~10 sec)
python src/14_figure5_operon_distance.py

# Step 15: Figure 6 — Wang 2021 family-level co-occurrence (~10 sec)
python src/15_figure6_wang2021.py

# Step 16: KEGG pathway co-membership test (requires Spark, ~15 min)
python src/16_kegg_pathway_comembership.py
```

Steps 1, 7, 11, 12, and 16 require access to the `kbase_ke_pangenome` tenant on BERDL. Steps 9–10 require an R environment with `phylolm` and `phytools` (create via `conda create -n r_phylo -c conda-forge r-base r-ape r-phytools r-phylolm --solver=classic`). All other steps run locally on CSV outputs in `data/`. Scripts 02–05, 08–10 apply a soil+plant species filter using `data/env_species_mapping.csv`; scripts 11, 12, and 16 filter after Spark collection.

All permutation tests and phylogenetic subsampling use a fixed random seed (42) for reproducibility. Python scripts use `np.random.seed(42)` or `np.random.default_rng(42)`; R scripts use `set.seed(42)`. Permutation p-values will differ slightly if re-run without the seed.

**Notebooks (NB01–NB07)** are inspection and display notebooks that load pre-computed CSVs from `data/`. They visualize and summarize results but do not perform primary analysis. The primary analysis code is in `src/*.py`. Re-running a notebook does not re-query BERDL — it renders cached CSV outputs.

### Figures

| Figure | Filename | Description |
|--------|----------|-------------|
| Figure 1 | `figure1_cooccurrence.png` | Multi-panel co-occurrence heatmap, core fractions, phylum-level, phenazine taxonomy |
| Figure 2 | `forest_plot.png` | Per-phylum log-OR forest plot (22 phyla, n≥20) |
| Figure 3 | `figure3_env_stratification.png` | Environmental stratification (v3 pan-bacterial only; not regenerated for v4) |
| Figure 4 | `figure4_phylo_correction.png` | Uncorrected vs phylo-corrected log-OR scatter |
| Figure 5 | `figure5_operon_distance.png` | Operon-distance distribution (observed vs null) |
| Figure 6 | `figure6_wang2021.png` | Wang 2021 per-family and per-siderophore validation |

## Authors

- Jing Tao (jingtao-lbl), Lawrence Berkeley National Laboratory
