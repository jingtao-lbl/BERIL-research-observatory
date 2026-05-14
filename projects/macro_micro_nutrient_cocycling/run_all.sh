#!/bin/bash
set -euo pipefail

# Macro-Micro Nutrient Co-cycling Pipeline Runner
# Steps 1, 7, 11, 12, 16 require on-cluster Spark access.
# Steps 9-10 require conda env r_phylo.
# All other steps run locally on cached CSVs in data/.

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

run_step() {
    local step="$1"; shift
    echo "━━━ Step $step: $* ━━━"
    "$@"
    echo ""
}

# Phase 1: Data extraction (Spark)
if [ ! -f data/species_gene_families.csv ]; then
    run_step 1 python3 src/01_extract_gene_families.py
else
    echo "━━━ Step 1: SKIP (data/species_gene_families.csv exists) ━━━"
fi

# Phase 2: Local analysis
run_step 2  python3 src/02_cooccurrence_stats.py
run_step 3  python3 src/03_core_accessory_enrichment.py
run_step 4  python3 src/04_phylogenetic_stratification.py
run_step 5  python3 src/05_figure.py
run_step 6  python3 src/06_positive_controls.py

# Phase 3: Environmental stratification (Spark)
if [ ! -f data/env_species_mapping.csv ]; then
    run_step 7 python3 src/07_environmental_stratification.py
else
    echo "━━━ Step 7: SKIP (data/env_species_mapping.csv exists) ━━━"
fi

run_step 8 python3 src/08_forest_plot.py

# Phase 4: Phylogenetic analysis (R)
run_step 9 conda run -n r_phylo Rscript src/09_phylo_signal.R \
    data/phylo/bac120_r214_pruned.tree \
    data/species_gene_families.csv \
    data/species_taxonomy.csv \
    data/phylo/phylo_signal.csv

run_step 10 conda run -n r_phylo Rscript src/10_phylo_logistic.R \
    data/phylo/bac120_r214_pruned.tree \
    data/species_gene_families.csv \
    data/phylo/pair_definitions.csv \
    data/phylo/phylo_logistic.csv

echo "━━━ Step 10a: Convergence check ━━━"
grep "FALSE" data/phylo/phylo_logistic.csv && echo "WARNING: Some models did not converge" || echo "All models converged"
echo ""

# Phase 5: Mechanistic tests (Spark)
run_step 11 python3 src/11_operon_distance.py
run_step 12 python3 src/12_wang2021_validation.py

# Phase 6: Figures
run_step 13 python3 src/13_figure4_phylo_correction.py
run_step 14 python3 src/14_figure5_operon_distance.py
run_step 15 python3 src/15_figure6_wang2021.py

# Phase 7: KEGG pathway (Spark)
run_step 16 python3 src/16_kegg_pathway_comembership.py

# Wang phylogenetic correction (requires extended CSV from step 12)
if [ -f data/species_gene_families_extended.csv ] && [ -f data/phylo/pair_definitions_wang.csv ]; then
    echo "━━━ Step 10b: Wang phylogenetic logistic regression ━━━"
    conda run -n r_phylo Rscript src/10_phylo_logistic.R \
        data/phylo/bac120_r214_pruned.tree \
        data/species_gene_families_extended.csv \
        data/phylo/pair_definitions_wang.csv \
        data/phylo/phylo_logistic_wang.csv
    echo ""
fi

echo "━━━ Pipeline complete ━━━"
echo "Figures: $(ls figures/*.png 2>/dev/null | wc -l) PNG files in figures/"
echo "Data:    $(ls data/*.csv 2>/dev/null | wc -l) CSV files in data/"
