# Macro-Micro Nutrient Co-cycling: Full Analysis and Enrichments

**Date:** May 7, 2026
**Author:** Jing Tao with Claude
**Type:** Research
**Project:** `macro_micro_nutrient_cocycling`
**Branch:** `projects/macro_micro_nutrient_cocycling`

---

## Summary

Executed the complete 4-step pangenome co-occurrence analysis (Substrate A on `kbase_ke_pangenome`) testing whether bacterial pangenomes encoding macro-nutrient acquisition (P, N) and phenazine biosynthesis genes disproportionately co-encode trace-metal handling genes. Phenazine biosynthesis was added as a 4th gene family group following the McRose & Newman (2021, Science 371:1033–1037) mechanism linking phenazine-mediated Fe(III)-oxyhydroxide dissolution to phosphate release.

Submitted for automated review via `/submit`. Addressed 6 of 7 reviewer issues (skipping issue #5 — set_index pattern). Completed three additional enrichments: (1) positive-control species check against 7 model organisms, (2) environmental stratification via `ncbi_env`, (3) per-phylum forest plot with log-odds ratios and 95% CIs.

**Key result:** P-acquisition and metal-handling genes co-occur significantly across 27,682 species (phi=0.110, OR=2.3, Fisher p=1.3×10⁻⁶⁵). N-fixation × Metal shows the strongest per-species association (OR=10.1). All 63 phenazine operon carriers encode both P-acquisition and metal-handling genes (100% overlap). Three distinct genomic signatures of macro-micro coupling: P-acquisition (core-enriched, 70.7%), N-fixation (moderately core, 63.4% — consistent with HGT), Metal-handling (bimodal — corA 73.5% core vs. HMA 23.6% accessory).

## Problem

After pivoting from Substrate C (fitnessbrowser — zero P-starvation experiments) to Substrate A (pangenome co-occurrence), needed to:

1. Define 23 gene families across 4 functional groups and extract presence/absence from 132.5M gene clusters
2. Compute pairwise co-occurrence statistics with appropriate null models
3. Test core vs. accessory genome enrichment via per-species Fisher + Stouffer meta-analysis
4. Stratify by GTDB taxonomy (phylum, class, family)
5. Validate with positive-control organisms, environmental metadata, and per-phylum effect sizes
6. Produce publication-quality figures and a full REPORT.md

## Solution

### Step 1: Gene Family Extraction (`src/01_extract_gene_families.py`)

Queried `kbase_ke_pangenome` via Spark SQL to extract 23 gene families across 4 functional groups for 27,682 species.

**Gene family definitions (final):**

| Group | Genes | Annotation source | Key pitfalls |
|-------|-------|-------------------|-------------|
| P-acquisition (9) | phoA, phoD_pfam, pstA/B/C/S, phnC/D/E | KO K01077, PFAM PF09423, gene names | pst genes have NO KO in bakta — must use gene name |
| N-fixation (2) | nifH, nifD | KO K02588, K02586 | PF00142 excluded — captures Fer4 superfamily (114% inflation) |
| Metal-handling (4) | copA, corA, feoB_pfam, HMA_pfam | KO K17686, gene name, PFAM PF02421/PF00403 | PFAM IDs versioned — must use `LIKE 'PFxxxxx.%'` |
| Phenazine (7) | phzA/B/D/F/G/S/M | KO K06998/K20940, gene names | Operon ≥3 phz genes filters 10,856→63 species |

**N-fixation definition change:** Initially included PF00142 (nifH_pfam) alongside KO-defined nifH/nifD. Reviewer flagged that PF00142 (Fer4_NifH) captures the broader Fer4 ferredoxin superfamily, inflating N-fixation species from 2,746 to 5,872 (114% increase). Restricted to KO-only for primary analysis; PF00142 reported as sensitivity check.

**Outputs:**
- `data/species_gene_families.csv` — 27,682 rows × 123 columns (presence/absence flags, counts, core/aux/singleton breakdowns per gene family)
- `data/species_taxonomy.csv` — GTDB taxonomy parsed to domain/phylum/class/order/family/genus
- `data/phenazine_operon_species.csv` — 63 species with ≥3 distinct phz gene families

### Step 2: Co-occurrence Statistics (`src/02_cooccurrence_stats.py`)

Computed Jaccard index, phi coefficient, Fisher's exact test, Benjamini-Hochberg FDR, and 1000-permutation null for all group pairs and 72 individual nutrient×metal gene pairs.

**Group-level results:**

| Pair | n_both | Jaccard | Phi | OR | Fisher p | Permutation Z |
|------|-------:|--------:|----:|---:|:--------:|--------------:|
| P × Metal | 20,683 | 0.769 | 0.110 | 2.30 | 1.3×10⁻⁶⁵ | 17.7 |
| N × Metal | 2,719 | 0.107 | 0.088 | 10.08 | 1.3×10⁻⁷¹ | 14.7 |
| Phz-operon × Metal | 63 | 0.002 | 0.014 | ∞ | 9.4×10⁻³ | 2.3 |
| P × N | 2,414 | 0.107 | 0.065 | 1.90 | 1.3×10⁻²⁹ | 10.7 |
| Phz-operon × P | 63 | 0.003 | 0.024 | ∞ | 1.5×10⁻⁶ | 4.0 |
| N × Phz-operon | 0 | 0.000 | −0.016 | 0.00 | 2.4×10⁻³ | — |

**Individual gene pair highlights (72 pairs, 64 FDR-significant at q<0.05):**
- Strongest positive: phzG × corA (1.79×, phi=0.057), phoD × HMA (1.72×, phi=0.084), nifH × HMA (1.42×, phi=0.077)
- Strongest negative: pstC × feoB (0.67×, phi=−0.256), pstS × feoB (0.64×, phi=−0.184)
- 8 non-significant pairs — all with small sample sizes (<50 co-occurrences)

**Sensitivity check (N-fixation breadth):** Including PF00142 expands N-fixation from 2,746 to 5,872 species. N × Metal phi changes from 0.088 to 0.107; OR drops from 10.1 to 4.0. The inflation includes ferredoxin-containing non-diazotrophs.

### Step 3: Core vs. Accessory Enrichment (`src/03_core_accessory_enrichment.py`)

Fisher's exact test per species comparing core/non-core fractions between nutrient and metal gene sets, combined with Stouffer meta-analysis.

**Three distinct genomic signatures:**

1. **P-acquisition — core-enriched.** 70.7% core (Stouffer meta-Z=68.3, p≈0). Individual genes: pstA 77%, pstC 76%, pstS 73%, phnC 66%, phoA 62%. Exception: phoD only 17.2% core (mobile element).

2. **N-fixation — moderately core.** 63.4% core (Stouffer meta-Z=3.2, p=0.001). Lower than P-acquisition, consistent with well-documented horizontal transfer of nif operons via symbiosis islands. **Critical fix:** Before restricting to KO-only, PF00142 inclusion caused Stouffer Z to be −4.85 (apparent accessory enrichment). After KO-only restriction, Z flipped to +3.25 — the PF00142 Fer4 ferredoxins were heavily accessory, masking the true nif signal.

3. **Metal-handling — bimodal.** corA 73.5% core (essential Mg²⁺/Co²⁺), copA 56.8% (intermediate), feoB 30.0% (accessory), HMA 23.6% (accessory). Core transport + accessory efflux pattern.

### Step 4: Phylogenetic Stratification (`src/04_phylogenetic_stratification.py`)

Stratified co-occurrence by GTDB phylum (34 with n≥50) and class (n≥20), plus plant-associated family analysis.

**Phylum-level:** Positive P × Metal phi in all 12 largest phyla except Spirochaetota (−0.085) and Chloroflexota (−0.074). Highest: Fusobacteriota (0.453), Cyanobacteriota (0.355), Bacillota (0.246).

**Plant-associated families:** Pseudomonadaceae, Rhizobiaceae, Burkholderiaceae, Streptomycetaceae, Xanthomonadaceae all show near-universal P + Metal prevalence (enrichment ratios ≈1.00 — ceiling effects).

**Phenazine operon carriers (63 species):** Concentrated in 9 families — Streptomycetaceae (24), Pseudomonadaceae (17), Streptosporangiaceae (10), Enterobacteriaceae (6 — all Xenorhabdus, insect-pathogenic), Xanthomonadaceae (2), Burkholderiaceae (1). Soil Actinomycetota (56%) + Pseudomonadota (32%) dominate.

### Figure Generation (`src/05_figure.py`)

Four-panel figure (`figures/figure1_cooccurrence.png`):
- (A) Phi coefficient heatmap for 72 nutrient × metal gene pairs
- (B) Core genome fraction bar chart per gene family
- (C) Phylum-level P × Metal co-occurrence (phi by phylum, sized by species count)
- (D) Taxonomic distribution of phenazine operon carriers (treemap by family)

### Reviewer Issues Addressed

After first `/submit`, reviewer flagged 7 issues. Addressed 6, skipped #5 (set_index pattern) per user decision:

1. **PF00142 scope** — Restricted N-fixation to KO K02588 (nifH) + K02586 (nifD). PF00142 as sensitivity check only. Re-ran NB02 and NB03.
2. **Negative Stouffer Z** — Fixed by KO-only restriction (Z flipped −4.85 → +3.25). Reframed as three-signature narrative.
3. **Multiple testing** — Applied Benjamini-Hochberg FDR to 72 Fisher tests. 64/72 significant at q<0.05. Added q-values to REPORT.md tables.
4. **Xenorhabdus** — Corrected "rhizosphere" to "insect-pathogenic (entomopathogenic nematode symbionts)."
5. *(Skipped — set_index pattern fine for project scope)*
6. **beril.yaml status** — Updated from "analysis" to "review."
7. **Orphan artifact** — `class_cooccurrence.csv` wired into analysis section 4 of REPORT.md.

### Enrichment 1: Positive-Control Species Check (`src/06_positive_controls.py`)

Looked up 8 well-characterized model organisms in GTDB R214. Streptomyces coelicolor A3(2) is absent from the dataset (not found by species name, genome accession GCF_000203835, or known synonyms S. violaceoruber / S. lividans).

**Results for 7 found species:**

| Organism | P-acquisition | N-fixation | Metal-handling | Phenazine |
|----------|---------------|------------|----------------|-----------|
| *P. fluorescens* | pstA/B/S, phnC/D/E | nifH | copA, corA | phzF |
| *P. protegens* | pstA/B/S | — | copA, corA, HMA | phzF |
| *B. diazoefficiens* | pstA/B/C/S, phnC/D/E | nifH, nifD | copA, corA | phzF, phzS |
| *S. meliloti* | phoD, pstA/B/C/S, phnC/D/E | nifH, nifD | copA, corA, feoB, HMA | phzF |
| *M. loti* | pstA/B/C/S, phnC/D/E | nifD | copA, corA, HMA | phzF |
| *S. coelicolor* | phoD, pstA/B/C/S | — | copA, corA | phzF |
| *M. extorquens* | phoA, pstA/B/C/S, phnC/D/E | — | copA, corA | phzF |
| *P. chlororaphis* | phoA, pstA/B/S | — | copA, corA | phzA/B/F/G (operon) |

**Key observations:**
- All 8 encode both P-acquisition and metal-handling — consistent with genome-wide signal
- Two known diazotrophs (B. diazoefficiens, S. meliloti) carry nifH+nifD as expected
- S. meliloti has broadest repertoire: all 4 metal families + 8/9 P families
- P. chlororaphis is only positive control with complete phenazine operon (phzA/B/F/G)
- GTDB uses reclassified genus Pseudomonas_E for fluorescens/protegens/chlororaphis
- S. coelicolor A3(2) is reclassified as s__Streptomyces_anthocyanicus in GTDB R214 (confirmed via GCF_008931305.1 in genome table)
- S. coelicolor carries only phzF (no operon) despite being a known phenazine producer — other phz genes may use different Bakta annotations

**Pitfall discovered:** Initial query used `LIKE '%coelicolor%'` which matched Cephaloticoccus (contains 'coccus' substring) instead of S. coelicolor. Fixed by using precise GTDB species prefixes. The species name `coelicolor` does not exist in GTDB R214 — it was reclassified as `anthocyanicus`. Found by searching `gtdb_metadata` for the newer RefSeq accession `RS_GCF_008931305.1`, which mapped to `s__Streptomyces_anthocyanicus`.

### Enrichment 2: Environmental Stratification (`src/07_environmental_stratification.py`)

Joined `ncbi_env` metadata to 27,009 species via genome biosample accessions (`ncbi_env.accession = genome.ncbi_biosample_id`).

**ncbi_env table structure:** Key-value store with columns: accession, attribute_name, content, display_name, harmonized_name, id, package_content. Not a columnar table. Key harmonized_name values: isolation_source (221K rows from 245K total), host (155K), env_broad_scale (72K).

**Environment classification:** Built keyword-based classifier for 6 broad categories:
- soil/rhizosphere: 3,409 species (soil, rhizosphere, rhizoplane, root, compost, sediment)
- plant-associated: 1,135 species (plant, leaf, stem, flower, endophyte, nodule, crop names)
- marine: 3,449 species (ocean, sea, coral, sponge, coastal, hydrothermal)
- freshwater/engineered: 2,059 species (lake, river, wetland, wastewater, activated sludge)
- human-associated: 3,497 species (gut, feces, blood, clinical, oral, skin)
- animal-associated: 840 species (insect, cattle, rumen, fish, poultry)
- other: 12,628 species (unclassifiable or ambiguous)

**Co-occurrence by environment:**

| Environment | n | P×Metal log-OR | P×Metal p | N×Metal log-OR | N×Metal p |
|-------------|--:|---------------:|----------:|---------------:|----------:|
| Soil/rhizosphere | 3,406 | +1.34 | 1.9×10⁻⁹ | +0.36 | 0.17 |
| Plant-associated | 1,134 | +1.51 | 1.3×10⁻⁶ | +0.94 | 0.064 |
| Marine | 3,449 | +0.26 | 7.8×10⁻³ | +1.37 | 4.2×10⁻²³ |
| Freshwater/eng. | 2,059 | +1.31 | 6.6×10⁻⁹ | +1.46 | 1.5×10⁻⁵ |
| Human-associated | 3,497 | +1.37 | 3.4×10⁻²⁰ | +0.97 | 7.1×10⁻⁶ |
| Animal-associated | 840 | +0.23 | 0.39 | +2.27 | 2.1×10⁻⁷ |

**Key findings:**
- P × Metal strongest in plant-associated (log-OR=+1.51) and soil (+1.34) — consistent with Fe-oxyhydroxide mineral surface hypothesis
- N × Metal shows different pattern: strongest in marine (+1.37) and animal-associated (+2.27); non-significant in soil (+0.36, p=0.17)
- Phz × Metal non-significant in all environments due to operon rarity (0–29 carriers per env)
- P × Metal is significant across 5/6 environments (not animal-associated) — the coupling is broadly distributed, not just soil-specific

### Enrichment 3: Per-Phylum Forest Plot (`src/08_forest_plot.py`)

Computed log-odds ratios with 95% CIs (Haldane-corrected for zero cells) for P×Metal, N×Metal, Phz×Metal across 34 GTDB phyla (n≥50 species).

**P × Metal:** Positive log-OR in 28/34 phyla. Strongest: Cyanobacteriota (+1.86, p=4.1×10⁻¹³), Bacillota (+1.64, p=8.0×10⁻²⁴), Pseudomonadota (+0.94, p=7.7×10⁻²¹). Five phyla significantly negative: Myxococcota (−2.63), Spirochaetota (−1.79), Verrucomicrobiota (−1.37), Planctomycetota (−0.98), Campylobacterota (−0.77).

**N × Metal:** Positive log-OR in 24/34 phyla. Strongest: Methanobacteriota (+5.00), Nanoarchaeota (+2.72), Cyanobacteriota (+2.11). Pattern phylogenetically broader than P × Metal — diazotrophy spans diverse lineages.

**Phz × Metal:** Only Pseudomonadota (+1.87, p=0.11) and Actinomycetota (+1.14, p=0.40) show positive trends; neither significant. All other phyla null or weakly negative. Consistent with phenazine operon concentration in just these two phyla.

**Output:** `figures/forest_plot.png` — three-panel horizontal bar chart with phyla sorted by log-OR, blue bars for p<0.05, gray for non-significant. Error bars show 95% CIs.

### Notebooks

| Notebook | Content | Status |
|----------|---------|--------|
| NB01_gene_family_extraction.ipynb | Data loading + summary stats (Spark-dependent) | Executed |
| NB02_cooccurrence_analysis.ipynb | Group + gene-pair co-occurrence, permutation test | Executed |
| NB03_core_accessory_phylogenetic.ipynb | Core/accessory, phylum, phenazine taxonomy | Executed |
| NB04_figure.ipynb | Figure 1 (4-panel) | Executed |
| NB05_positive_controls.ipynb | Positive-control species table | Executed |
| NB06_environmental_stratification.ipynb | Environmental co-occurrence | Executed |
| NB07_forest_plot.ipynb | Forest plot (Figure 2) | Executed |

**NB04 recurring issue:** `create_notebooks.py` generates NB04 with an `exec(open('05_figure.py'))` cell. The `os.path.dirname(__file__)` references in the figure script break when exec'd from notebook context (notebooks have no `__file__`). Must replace with relative `'..'` paths after notebook creation, every time.

### REPORT.md Structure (final)

1. Key Findings (5 bullets)
2. Background (McRose-Newman mechanism, Fe-oxyhydroxide biogeochemistry)
3. Hypothesis
4. Data and Methods (dataset, gene family definitions, statistical tests)
5. Results
   - 5.1 Group-level co-occurrence + sensitivity check
   - 5.2 Individual gene pair highlights (FDR-corrected)
   - 5.3 Core vs. accessory structure (three signatures)
   - 5.4 Phylogenetic distribution (phylum, class, family, phenazine taxonomy)
   - 5.5 Positive-control species check (8 model organisms)
   - 5.6 Environmental stratification (6 broad environments)
   - 5.7 Per-phylum forest plot (Figure 2)
6. Interpretation
7. Limitations (6 items including PF00142, phylogenetic non-independence)
8. Future Directions
9. References
10. Data and Reproducibility (Figure 1, Figure 2)

## Data Files

All in `projects/macro_micro_nutrient_cocycling/data/` (gitignored — must use `git add -f`):

| File | Rows | Description |
|------|-----:|-------------|
| species_gene_families.csv | 27,682 | Presence/absence + counts per gene family per species |
| species_taxonomy.csv | 27,682 | Parsed GTDB taxonomy |
| phenazine_operon_species.csv | 63 | Species with ≥3 phz genes |
| cooccurrence_matrix.csv | 21 | Group-level pairwise stats |
| pairwise_detail.csv | 72 | Gene-pair stats with FDR q-values |
| contingency_tables.txt | — | 2×2 tables for all group pairs |
| core_enrichment_summary.csv | 3 | Per-group Stouffer meta-analysis |
| core_enrichment_P_genes.csv | varies | Per-species Fisher for P×Metal |
| core_enrichment_N_genes.csv | varies | Per-species Fisher for N×Metal |
| core_enrichment_Phz_genes.csv | varies | Per-species Fisher for Phz×Metal |
| phylum_cooccurrence.csv | 34+ | Phylum-level co-occurrence stats |
| class_cooccurrence.csv | varies | Class-level co-occurrence stats |
| phenazine_operon_taxonomy.csv | 63 | Phz operon carriers with full taxonomy |
| positive_controls.csv | 8 | Model organism gene family profiles |
| env_species_mapping.csv | 27,017 | Species → primary environment assignment |
| env_cooccurrence.csv | 21 | Co-occurrence stats by environment |
| forest_plot_data.csv | 102 | Per-phylum log-OR with CIs for 3 pairs |

## Figures

| File | Description |
|------|-------------|
| figures/figure1_cooccurrence.png | 4-panel: phi heatmap, core fractions, phylum stratification, phz taxonomy |
| figures/forest_plot.png | 3-panel forest plot: per-phylum log-OR for P×Metal, N×Metal, Phz×Metal |

## Pitfalls Discovered

1. **PFAM IDs versioned in bakta:** `bakta_pfam_domains.pfam_id` stores `PF00403.26` not `PF00403`. Must use `LIKE 'PFxxxxx.%'`.
2. **Pst operon has no KEGG KO:** pstA/B/C/S must be searched by `bakta_gene_name`, not KO number.
3. **PF00142 captures Fer4 superfamily:** Using it as N-fixation proxy inflates species 114% and inverts Stouffer Z. Use KO K02588 (nifH) + K02586 (nifD) only.
4. **ncbi_env is key-value, not columnar:** Filter on `harmonized_name IN ('isolation_source', 'env_broad_scale', 'host')`. Many entries contain 'missing', 'not collected', 'not applicable' — must filter.
5. **ncbi_env → species join path:** `ncbi_env.accession = genome.ncbi_biosample_id` → `genome.gtdb_species_clade_id`.
6. **S. coelicolor reclassified in GTDB R214:** Not found under its NCBI name. Reclassified as `s__Streptomyces_anthocyanicus`. Found by searching `gtdb_metadata` for the newer RefSeq accession `RS_GCF_008931305.1`. The older accession `GCF_000203835` is not in the dataset at all.
7. **CSVs gitignored:** `.gitignore` rule `projects/*/data/**/*.csv` blocks them. Must use `git add -f`.
8. **NB04 exec() path issue:** `os.path.dirname(__file__)` breaks in notebook context. Must replace with relative paths.
9. **GTDB reclassifications:** Pseudomonas fluorescens/protegens/chlororaphis are in genus Pseudomonas_E in GTDB R214. Search queries must use GTDB names, not NCBI names.
10. **Positive control LIKE gotcha:** `LIKE '%coelicolor%'` matches Cephaloticoccus (contains 'coccus' overlap). Use precise prefix patterns like `LIKE 's__Streptomyces%coelicolor%'`.

## Remaining Work (as of 2026-05-07 end of session)

1. **Commit all enrichment changes:** src/06-08, notebooks NB05-07, figures/forest_plot.png, updated REPORT.md, data files (with `git add -f`)
2. **GitHub authentication:** `gh auth login` needed on JupyterHub before push — `git push` fails with "could not read Username"
3. **Re-submit for review:** Run `/submit` for fresh REVIEW.md incorporating enrichments
4. **RESEARCH_PLAN.md:** May need v4 update to document enrichment additions

## Files Changed

| File | Action | Description |
|------|--------|-------------|
| `src/01_extract_gene_families.py` | Created | Spark SQL extraction of 23 gene families |
| `src/02_cooccurrence_stats.py` | Created + Modified | Co-occurrence; added FDR, sensitivity check, KO-only N-fixation |
| `src/03_core_accessory_enrichment.py` | Created + Modified | Core/accessory; fixed N-fixation to KO-only |
| `src/04_phylogenetic_stratification.py` | Created | Phylum/class/family stratification |
| `src/05_figure.py` | Created | 4-panel figure |
| `src/06_positive_controls.py` | Created + Modified | Positive controls; fixed to precise GTDB matching |
| `src/07_environmental_stratification.py` | Created | Environmental stratification via ncbi_env |
| `src/08_forest_plot.py` | Created | Per-phylum forest plot |
| `src/create_notebooks.py` | Created | NB01-NB04 generation |
| `src/create_enrichment_notebooks.py` | Created | NB05-NB07 generation |
| `notebooks/NB01-NB07` | Created | All notebooks with saved outputs |
| `figures/figure1_cooccurrence.png` | Created | 4-panel main figure |
| `figures/forest_plot.png` | Created | 3-panel forest plot |
| `REPORT.md` | Created + Modified | Full report with all 7 results sections |
| `RESEARCH_PLAN.md` | Created (v1→v3) | Research plan with revision history |
| `README.md` | Created | Project overview, status: Completed |
| `REVIEW.md` | Generated | Automated review from first /submit |
| `beril.yaml` | Created | Status: review |
| `requirements.txt` | Created | numpy, pandas, scipy, matplotlib, pyspark |
| `data/*.csv` | Created | 17 data files |
| `memory/20260507g_*.md` | Created | This log |

## Verification

- All 8 analysis scripts executed successfully on-cluster via `berdl_notebook_utils.setup_spark_session.get_spark_session()`
- All 7 notebooks (NB01-NB07) executed successfully via `jupyter nbconvert --execute`
- Fisher's exact test results cross-checked: pairwise_detail.csv has 72 rows, 64 FDR-significant
- Positive controls validated: all 8 targets encode P+Metal. B. diazoefficiens carries nifH+nifD (true diazotroph), S. meliloti carries all 4 metal families, P. chlororaphis has phz operon, S. coelicolor (= S. anthocyanicus) found via GCF_008931305.1
- Environmental join: 447,570 annotations across 27,017 species (97.5% of total 27,682)
- Forest plot: 102 rows = 34 phyla × 3 pairs; figure renders correctly with CIs
- REPORT.md contains all 7 results sections plus interpretation, limitations, references
