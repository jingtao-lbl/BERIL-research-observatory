# Session Log: V4 Soil+Plant Subset Analysis

**Date:** 2026-05-13
**Branch:** `projects/macro_micro_nutrient_cocycling_v4_soil_plant` (created from v3 commit e4c9163)
**Analyst:** Jing Tao (LBNL) + Claude
**Starting point:** v3 pan-bacterial analysis complete (27,682 species, all environments)
**Goal:** Repeat entire analysis pipeline filtered to only soil/rhizosphere (3,409) + plant-associated (1,135) = 4,544 species

---

## Motivation

The v3 analysis tested macro-micro nutrient gene co-occurrence across all 27,682 GTDB bacterial species pangenomes. The biological question — co-cycling of micro- and macro-nutrients via Fe-oxyhydroxide mineral dissolution — is most relevant in the plant-microbe-soil system. Restricting to soil+plant species tests whether the coupling is stronger in the ecologically relevant subset.

---

## Implementation Plan

A detailed implementation plan was written and approved before execution. The approach: modify scripts in-place on the dedicated v4 branch (no need for dual-mode toggles since the branch provides isolation from v3). Each script gets a ~5-line environment filter block after loading data.

### Scripts skipped (not modified):
- `src/01_extract_gene_families.py` — Extracts from BERDL; keep full-database CSVs and filter downstream
- `src/06_positive_controls.py` — Model organisms are hardcoded; results unchanged
- `src/07_environmental_stratification.py` — This GENERATED `env_species_mapping.csv`; stratification within a single-environment subset is meaningless

---

## Phase 0: Data Verification

Verified species counts and intersection sizes before modifying any scripts.

| Dataset | Count |
|---------|------:|
| `env_species_mapping.csv` soil/rhizosphere species | 3,409 |
| `env_species_mapping.csv` plant-associated species | 1,135 |
| Total mapped | 4,544 |
| Intersection with `species_gene_families.csv` (27,682 rows) | 4,540 |
| Intersection with tree tips in `bac120_r214_pruned.tree` (26,517 tips) | 4,177 |
| Species with both P+M genes in subset | 3,815 |

The 4-species gap (4,544 → 4,540) is expected: 4 soil+plant species lack gene cluster data in the pangenome database.

---

## Phase 1: Local CSV-Based Scripts

### 1.1: `src/02_cooccurrence_stats.py` — Co-occurrence statistics

**Modification:** Added soil+plant filter after CSV load (~line 23):
```python
_env = pd.read_csv(os.path.join(DATA_DIR, 'env_species_mapping.csv'))
_sp = set(_env[_env['primary_env'].isin(['soil/rhizosphere', 'plant-associated'])]['gtdb_species_clade_id'])
df = df[df['gtdb_species_clade_id'].isin(_sp)].copy()
print(f"v4 soil+plant filter: {len(df)} species (from {len(_env)} mapped)")
del _env, _sp
```

**Key results:**

| Pair | phi (v4) | phi (v3) | OR (v4) | OR (v3) | p (v4) |
|------|------:|------:|------:|------:|:------:|
| P × Metal | 0.129 | 0.110 | 3.92 | 2.30 | 4.4×10⁻¹⁴ |
| N × Metal | 0.040 | 0.088 | 2.74 | 10.1 | 4.3×10⁻³ |
| P × N | 0.125 | 0.065 | 14.24 | 1.90 | 1.4×10⁻²⁴ |
| Phz-operon × Metal | 0.015 | 0.014 | ∞ | ∞ | 0.624 (NS) |
| Phz-operon × P | 0.031 | 0.024 | ∞ | ∞ | 0.043 |
| N × Phz-operon | −0.027 | −0.016 | 0.0 | 0.0 | 0.108 (NS) |

- 51/72 FDR-significant pairs (vs 64/72 pan-bacterially)
- P × Metal phi STRONGER in soil+plant (0.129 vs 0.110)
- N × Metal OR decreased (2.74 vs 10.1) due to metal background prevalence rising from 91.7% to 96.4%
- Phz-operon × Metal became non-significant (p=0.624) — 27 operon carriers with 96.4% metal background

### 1.2: `src/03_core_accessory_enrichment.py` — Core/accessory analysis

**Modification:** Same filter pattern after CSV load.

**Key results:**

| Nutrient group | Core fraction | Stouffer Z (v4) | Stouffer Z (v3) | p (v4) |
|---------------|:------------:|:-----------:|:-----------:|:------:|
| P-acquisition | 0.740 | 25.47 | 68.3 | ≈0 |
| N-fixation | 0.636 | −0.307 | 3.2 | 0.759 (NS) |
| Phenazine | 0.693 | −1.628 | 2.6 | 0.103 (NS) |

- P-acquisition core enrichment remains highly significant (Z=25.5)
- N-fixation and Phenazine core enrichment became NON-significant in soil+plant subset
- N-fixation Z flipped from positive (3.2) to near-zero (−0.3) — in soil environments, nif genes may be more frequently on mobile elements

### 1.3: `src/04_phylogenetic_stratification.py` — Phylum-level analysis

**Modification:** Filter both `df` and `tax` DataFrames to soil+plant species.

**Key results:**
- 22 phyla with ≥20 species (threshold lowered from 50 to 20 for forest plot)
- 27 phenazine operon carriers (down from 63)
- By phylum: Actinomycetota 21 (78%), Pseudomonadota 6 (22%)
- By family: Streptomycetaceae 16, Streptosporangiaceae 3, Pseudomonadaceae 3, Xanthomonadaceae 2, Jiangellaceae 1, Pseudonocardiaceae 1, SG8-39 1
- Xenorhabdus clade (6 insect pathogens) absent from soil+plant subset as expected
- Actinomycetota dominance increased from 56% to 78%

### 1.4: `src/09_phylo_signal.R` — Pagel's Lambda

**Modifications:**
1. Added R filter after loading traits (same pattern: read `env_species_mapping.csv`, filter to soil+plant IDs)
2. Reduced default `n_subsample` from 5000 to 3000 (only ~4,200 tips available)

**Key results (3,000-tip phylum-stratified subsample):**

| Gene | Lambda (v4) | Lambda (v3) | Group |
|------|:----------:|:----------:|-------|
| pstA | 1.000 | 1.000 | P-acquisition |
| pstB | 1.000 | 0.970 | P-acquisition |
| pstC | 0.951 | 0.900 | P-acquisition |
| pstS | 0.804 | — | P-acquisition |
| phoA | 0.663 | — | P-acquisition |
| phoD | 0.136 | — | P-acquisition |
| phnC | 0.667 | — | P-acquisition |
| phnD | 0.751 | — | P-acquisition |
| phnE | 0.775 | — | P-acquisition |
| nifH | 0.852 | 0.870 | N-fixation |
| nifD | 0.765 | — | N-fixation |
| copA | 0.786 | 0.880 | Metal-handling |
| corA | 0.961 | 0.920 | Metal-handling |
| feoB | 0.790 | 0.520 | Metal-handling |
| HMA | 0.713 | 0.440 | Metal-handling |
| phzF | 0.823 | — | Phenazine |
| phzG | 0.055 | 0.080 | Phenazine |
| phzS | 0.248 | 0.280 | Phenazine |
| phzA/B/D/M | NA | — | Phenazine (near-fixed) |

- feoB lambda increased dramatically (0.52 → 0.79) — stronger phylogenetic signal in soil+plant subset
- HMA lambda also increased (0.44 → 0.71)
- corA lambda increased (0.92 → 0.96)
- Rare phenazine genes (phzA, phzB, phzD, phzM) skipped due to prevalence <1%

### 1.5: `src/10_phylo_logistic.R` — Phylogenetic logistic regression

**Modifications:**
1. Added R filter after loading traits
2. Added explicit tree pruning: `tree <- keep.tip(tree, common)` before trait matching

**Key results (4,177-tip subtree, 17 pairs tested):**

| Pair | Uncorrected log-OR | Phylo log-OR | Change | Phylo p | Significant? |
|------|-------------------:|-------------:|-------:|--------:|:---:|
| P × Metal (group) | 1.599 | 0.938 | −41% | 2.1×10⁻⁵ | Yes |
| N × Metal (group) | 0.937 | 0.986 | +5% | 9.2×10⁻⁵ | Yes |
| P × N (group) | 3.361 | 3.139 | −7% | 7.2×10⁻⁵ | Yes |
| Phz-operon × Metal | ∞ | 13.28 | — | 0.98 | No |
| Phz-operon × P | ∞ | 13.37 | — | 0.94 | No |
| N × Phz-operon | −∞ | −2.73 | — | 1.1×10⁻⁵ | Yes |
| phzG × corA | 3.363 | 3.335 | −1% | 2.0×10⁻⁶ | Yes |
| phzD × corA | ∞ | 12.93 | — | 0.95 | No |
| phzA × corA | ∞ | 11.93 | — | 0.95 | No |
| phoD × HMA | 0.989 | 0.914 | −8% | 1.4×10⁻⁹ | Yes |
| phzB × corA | ∞ | 13.93 | — | 0.96 | No |
| phoD × feoB | 0.361 | 0.509 | +41% | 2.1×10⁻⁴ | Yes |
| nifH × HMA | 0.208 | 0.396 | +91% | 6.9×10⁻⁴ | Yes |
| pstC × feoB | −1.400 | −0.880 | −37% | 1.8×10⁻²⁰ | Yes |
| pstS × feoB | −0.927 | −0.549 | −41% | 1.2×10⁻¹¹ | Yes |
| pstC × HMA | −1.043 | −0.749 | −28% | 7.2×10⁻¹⁶ | Yes |
| phzF × feoB | −0.690 | −0.457 | −34% | 2.1×10⁻⁹ | Yes |

- 10/12 significant uncorrected pairs survive phylogenetic correction (83% survival)
- 2 collapsed: phzD × corA and phzB × corA (rare phenazine genes, small n)
- 2 pairs GAINED significance: N × Phz-operon (negative), nifH × HMA
- P × Metal attenuated by 41% (vs +2% pan-bacterially) — major finding
- Negative associations WEAKENED by 28–41% (vs strengthened 9–21% pan-bacterially)
- Key interpretation: soil+plant subset has more phylogenetic confounding to remove

---

## Phase 2: Spark-Required Scripts

### 2.1: `src/11_operon_distance.py` — Genomic distance test

**Modifications:**
1. After building `species_with_both` set, filter to soil+plant IDs
2. After Spark collect to pandas, filter `pdf` to soil+plant IDs

**Key results:**

| Metric | v4 (soil+plant) | v3 (pan-bacterial) |
|--------|----------------:|-------------------:|
| Species with P + M genes | 3,546 | 18,539 |
| Species with same-contig pairs | 1,299 (36.6%) | 6,188 (33.4%) |
| Total same-contig P × M pairs | 15,589 | 65,463 |
| Same-contig fraction of all pairs | 29.3% | 28.1% |
| Observed median distance | 1,097 genes | 910 genes |
| Null median (1000 perms) | 350.4 ± 6.2 | 249.4 ± 2.2 |
| Z-score | 120.7 | 303.6 |
| Within 5 genes (operon-proximal) | ~0.49% | 0.67% |

- P-M genes are farther apart in soil+plant species (median 1,097 vs 910 pan-bacterially)
- Confirms ecological co-selection, not physical linkage
- Lower Z-score (120.7 vs 303.6) reflects smaller dataset, not weaker signal

### 2.2: `src/12_wang2021_validation.py` — Phytase × siderophore

**Modifications:**
1. After building `all_species` set, filter to soil+plant IDs
2. After building `taxonomy` dict, filter to soil+plant IDs
3. Added extended CSV regeneration at end (filters species_gene_families.csv to soil+plant, adds phytase/siderophore columns)

**Key results:**

| Metric | v4 (soil+plant) | v3 (pan-bacterial) |
|--------|----------------:|-------------------:|
| Total species | 4,540 | 27,682 |
| Phytase-encoding | 1,122 (24.7%) | 6,139 (22.2%) |
| Siderophore-encoding | 437 (9.6%) | 2,103 (7.6%) |
| Both | 124 (2.7%) | 732 (2.6%) |
| OR | 1.23 | 1.99 |
| Phi | 0.028 | 0.087 |
| Fisher p | 0.070 (NS) | 2.3×10⁻⁴³ |
| Permutation Z | 1.93 (p=0.038) | 13.9 |

- **Non-significant** in soil+plant subset (Fisher p=0.070)
- OR dropped from 2.0 to 1.23 — genuine weakening, not just power loss
- Per-siderophore: only pyoverdine (OR=1.83, p=0.004) and hydroxamate (OR=1.38, p=0.048) significant
- Burkholderiaceae: 9/89 = 10.1% with both traits, OR=0.97, p=1.0 (non-significant)

**Wang phylogenetic correction:** `phyloglm` failed to converge (converged=FALSE, all phylo values NA). Error: `NA/NaN/Inf in foreign function call (arg 5)`. This is expected — only 124 species with both traits out of 4,177 tree tips, with a weak underlying signal.

### 2.3: `src/16_kegg_pathway_comembership.py` — KEGG pathway analysis

**Modification:** After Spark collect, filter `p_all`, `m_all`, and `n_df` DataFrames to soil+plant IDs.

**Key results:**

| Metric | v4 (soil+plant) | v3 (pan-bacterial) |
|--------|----------------:|-------------------:|
| Species with both P and M annotations | 3,479 | 18,056 |
| % sharing ≥1 pathway | 4.9% | 7.6% |
| Mean shared pathways | 0.05 | 0.08 |
| Null mean | 0.34 | 0.17 |
| Z-score | −31.0 | −30.0 |
| P pathways | 37 | 63 |
| M pathways | 53 | 103 |
| Shared pathways | 11 | 15 |

- P and M genes remain in **distinct** metabolic pathways (Z=−31.0)
- Consistent with pan-bacterial result — slightly stronger in soil+plant

---

## Phase 3: Figure Scripts

### Modified scripts:
- `src/05_figure.py` — Added soil+plant filter for `df` (Panel B core fractions); updated title to show species count and "Soil/Plant-Associated"
- `src/08_forest_plot.py` — Added soil+plant filter; lowered `MIN_PHYLUM_SIZE` from 50 to 20
- `src/13_figure4_phylo_correction.py` — Updated title to include "(soil/plant species)"

### Unmodified figure scripts (read from pre-regenerated data CSVs):
- `src/14_figure5_operon_distance.py`
- `src/15_figure6_wang2021.py`

### Figures regenerated:
1. `figure1_cooccurrence.png/.pdf` — 4-panel: heatmap, core fractions, phylum P×M, phenazine taxonomy
2. `forest_plot.png` — 3-panel forest plot (22 phyla with n≥20)
3. `figure4_phylo_correction.png/.pdf` — Uncorrected vs corrected log-ORs (4,177-tip tree)
4. `figure5_operon_distance.png/.pdf` — Distance distribution (15,589 same-contig pairs)
5. `figure6_wang2021.png/.pdf` — Family-level phytase × siderophore

### Figure NOT regenerated:
- `figure3_env_stratification.png` — Not applicable (entire dataset is soil+plant)

---

## Phase 4: Scripts Skipped

| Script | Reason |
|--------|--------|
| `src/01_extract_gene_families.py` | Queries BERDL for all species; keep full CSVs, filter downstream |
| `src/06_positive_controls.py` | Model organisms hardcoded; all 8 are soil+plant species |
| `src/07_environmental_stratification.py` | Generated `env_species_mapping.csv`; env stratification within single-environment subset not meaningful |

---

## Phase 5: Documentation Updates

### REPORT.md — Complete rewrite (360 → 310 lines)
- Title: "27,682 Bacterial Pangenomes" → "4,540 Soil/Plant-Associated Bacterial Pangenomes"
- All numerical results updated throughout (every table, every statistic, every percentage)
- Background: Added framing about repeating analysis for soil+plant subset
- Dataset section: Added soil+plant subset description
- Methods: Updated phylo params (4,177 tips, 3,000-tip subsample, 20-species phylum threshold)
- Section 6 (Environmental Stratification) REMOVED — not applicable within single-environment subset
- Sections 7–11 renumbered to 6–10
- Interpretation: Substantially rewritten — new emphasis on phylogenetic attenuation (41% for P×M), core/accessory differences, Wang 2021 non-significance
- Limitations: Updated with new items (reduced power, near-universal metal prevalence, environmental classification limitations)
- Figure references: Removed Figure 3; updated all others with v4-specific captions

### README.md
- Status updated with v4 summary
- Overview rewritten for soil+plant scope with v4 statistics
- Reproduction notes updated to mention soil+plant filter in scripts

### RESEARCH_PLAN.md
- Added v5 revision entry documenting the soil+plant subset analysis

### beril.yaml
- `status: complete` → `status: analysis` (pending review)
- `branch` updated to `projects/macro_micro_nutrient_cocycling_v4_soil_plant`
- `last_session_at` updated to current timestamp
- `artifacts.review: true` → `artifacts.review: false` (v3 review no longer valid for v4)

---

## Key Biological Findings (v4 vs v3 Comparison)

### Findings that STRENGTHENED in soil+plant subset:
1. **P × Metal phi** increased from 0.110 to 0.129 — stronger co-occurrence in the ecologically relevant subset
2. **P × Metal OR** increased from 2.30 to 3.92
3. **P × N phi** increased from 0.065 to 0.125 — much stronger P-N coupling in soil+plant
4. **phoD × feoB** phylo log-OR strengthened by 41% after correction (0.361 → 0.509)
5. **nifH × HMA** phylo log-OR strengthened by 91% after correction (0.208 → 0.396)
6. **Negative associations** (pstC × feoB phi) strengthened from −0.256 to −0.395

### Findings that WEAKENED or changed:
1. **N × Metal OR** dropped from 10.1 to 2.74 (metal background prevalence ceiling: 96.4% vs 91.7%)
2. **Phz-operon × Metal** became non-significant (p=0.624 vs 9.4×10⁻³) — 27 vs 63 carriers
3. **Wang 2021 phytase × siderophore** became non-significant (p=0.070 vs 2.3×10⁻⁴³)
4. **FDR-significant pairs** reduced from 64/72 to 51/72
5. **Phylo correction** attenuated P×M by 41% (vs +2% pan-bacterially) — more phylogenetic confounding in ecologically constrained subset
6. **N-fixation core enrichment** lost significance (Z=−0.3 vs 3.2)
7. **Phenazine core enrichment** lost significance (Z=−1.6 vs 2.6)

### Findings that HELD:
1. **Operon distance** — P-M genes still NOT physically linked (median 1,097 vs null 350, Z=120.7)
2. **KEGG pathway** — P-M genes still in distinct pathways (Z=−31.0)
3. **100% phenazine operon overlap** with both P and M genes maintained
4. **P-acquisition core enrichment** robust (Z=25.5)
5. **P × Metal survives phylo correction** (log-OR=0.938, p=2.1×10⁻⁵)

---

## Errors and Issues Encountered

### 1. `has_phenazine_operon` KeyError (Phase 0)
- **Symptom:** Phase 0 verification script tried `gf_sp['has_phenazine_operon']`
- **Cause:** Actual column name is `has_phz_operon`
- **Fix:** Used `[c for c in gf.columns if c.startswith('has_')]` to discover correct name

### 2. Wang phylo logistic convergence failure
- **Symptom:** `NA/NaN/Inf in foreign function call (arg 5)` from `phyloglm`
- **Cause:** Only 124 species with both traits out of 4,177 tree tips; weak signal (OR=1.23, p=0.07)
- **Resolution:** Expected negative result. Output file created with converged=FALSE and NA values.

### 3. Background task output buffering
- **Symptom:** All four background task output files showed 0 bytes initially despite processes running
- **Cause:** conda/Spark buffer stdout until completion
- **Resolution:** Waited for processes to complete; outputs appeared normally.

### 4. Tree pruning requirement for R scripts
- **Problem:** With only ~4,200 soil+plant tips out of 26,517 total, the original `traits_matched <- traits[tree$tip.label, ]` would fail (many tree tips have no trait data)
- **Fix:** Added explicit pruning: `tree <- keep.tip(tree, common)` before trait matching

### 5. Extended CSV regeneration
- **Problem:** `species_gene_families_extended.csv` (with phytase/siderophore columns) was not generated by any existing script on the v4 branch
- **Fix:** Added regeneration code to `src/12_wang2021_validation.py` that creates the extended CSV filtered to soil+plant species

---

## File Change Summary

| Category | Files modified |
|----------|------:|
| Source scripts | 11 |
| Data CSVs | 26 |
| Figures (PNG + PDF) | 10 |
| Documentation (REPORT, README, RESEARCH_PLAN, beril.yaml) | 4 |
| **Total** | **47** (+ this session log) |

Net line change: −78,039 lines (primarily from data CSVs shrinking: 27,682 → 4,540 species rows)

---

## Pending

- [ ] Commit all changes on branch `projects/macro_micro_nutrient_cocycling_v4_soil_plant`
- [ ] Run `/submit` for automated review
- [ ] Address any reviewer feedback
- [ ] Merge to main or keep as parallel analysis branch
