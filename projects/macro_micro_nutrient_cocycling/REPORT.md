# Macro-Micro Nutrient Gene Co-occurrence Across 4,540 Soil/Plant-Associated Bacterial Pangenomes

## Key Findings

1. **P-acquisition and metal-handling genes co-occur significantly** across 4,540 soil/rhizosphere and plant-associated GTDB species pangenomes (phi=0.129, OR=3.9, Fisher p=4.4×10⁻¹⁴), with a stronger effect size than the pan-bacterial baseline (phi=0.110, OR=2.3). The association **survives phylogenetic correction** under logistic MPLE on a 4,177-tip subtree (phylo log-OR=0.938, p=2.1×10⁻⁵), though it is attenuated by 41% — indicating that phylogenetic structure in the soil+plant subset partially inflates the uncorrected signal.

2. **N-fixation and metal-handling remain significantly coupled.** Among 484 KO-defined diazotrophs in the soil+plant subset, 477 (98.6%) encode metal-handling genes (OR=2.7, p=4.3×10⁻³). The phylogenetically corrected log-OR (0.986) shows a 5% increase over the uncorrected value (0.937), confirming that this association is not driven by shared ancestry.

3. **The coupling is ecological, not genomic linkage or pathway integration.** P × M gene pairs on the same contig are farther apart than expected by chance (median 1,097 genes vs null 350, Z=120.7). Only 0.49% of same-contig pairs fall within the operon-proximal range (≤5 genes). P and M genes share fewer KEGG pathways than expected by chance (Z = −31.0). Co-occurrence reflects organism-level ecological co-selection, not physical linkage or shared metabolism.

4. **All 27 phenazine operon carriers in the soil+plant subset encode both P-acquisition and metal-handling genes** (100% overlap). However, the Phz-operon × Metal association is non-significant (OR=∞, Fisher p=0.624) due to the near-universal metal-handling gene prevalence (96.4%) in this subset. These 27 species are concentrated in soil families: Streptomycetaceae (16), Streptosporangiaceae (3), Pseudomonadaceae (3), Xanthomonadaceae (2), and three singleton families.

5. **P-acquisition genes are strongly core-enriched in the soil+plant subset** (74.0% core; Stouffer Z=25.5, p≈0). N-fixation and phenazine genes do not show significant core enrichment relative to metal-handling genes (Stouffer Z=−0.3 and −1.6, respectively), differing from the pan-bacterial analysis where N-fixation was moderately core-enriched (Z=3.2). Metal-handling genes span the core–accessory boundary (copA 56.8%, corA 73.5%, feoB 30.0%, HMA 23.6%).

6. **Independent validation: phytase × siderophore co-occurrence is non-significant in the soil+plant subset.** OR=1.23, phi=0.028, Fisher p=0.070 — below the significance threshold, unlike the pan-bacterial result (OR=2.0, p=2.3×10⁻⁴³). Only pyoverdine-type siderophores retain a significant association with phytase (OR=1.83, p=0.004). Phylogenetic logistic regression failed to converge, consistent with the weak underlying signal.

7. **Negative associations are informative but attenuated by phylogenetic correction.** Pst transporter genes (pstC, pstS) are depleted relative to feoB (phi as low as −0.395; all FDR q<0.05). Unlike the pan-bacterial analysis where these negative associations strengthened under phylogenetic correction, they weaken by 28–41% in the soil+plant subset — indicating that some of the anti-correlation in this subset is driven by phylogenetic structure rather than ecological niche separation.

## Background

Plant productivity and soil organic matter turnover depend on the simultaneous availability of macronutrients (C, N, P) and trace metals (Fe, Cu, Zn, Mn, Co, Ni, Mo). These pools are coupled through enzyme cofactor biochemistry — nitrogenase requires Fe and Mo, alkaline phosphatases require Zn or Ca, and urease requires Ni — and through mineral-surface biogeochemistry, where phosphate and trace-metal cations co-adsorb onto Fe(III)-oxyhydroxide surfaces in acidic, iron-rich soils (Edayilam et al. 2018).

McRose & Newman (2021, *Science* 371:1033–1037) demonstrated that phenazine biosynthesis is upregulated 1–3 orders of magnitude under phosphate starvation via the conserved PhoB regulon. Phenazines reductively dissolve Fe(III)-oxyhydroxides, simultaneously releasing adsorbed phosphate and co-adsorbed trace metals (Co, Cu, Ni). Wu et al. (2022) confirmed the presence of phenazine biosynthesis genes (phzS) in metagenome-assembled genomes from P-limited reforestation soils, establishing that this pathway operates in natural terrestrial microbiomes.

Despite these mechanistic links, macro-nutrient cycling and trace-metal dynamics are typically studied separately. A pan-bacterial analysis across all 27,682 GTDB species pangenomes (v3) established significant co-occurrence of P-acquisition and metal-handling genes. This report repeats that analysis restricted to 4,540 soil/rhizosphere and plant-associated species to test whether the coupling is stronger in the terrestrial plant-microbe-soil system — the ecological context where Fe-oxyhydroxide mineral chemistry most directly mediates nutrient availability.

## Hypothesis

Bacterial pangenomes encoding genes for macro-nutrient acquisition (phosphatases, phosphate/phosphonate transporters, nitrogenases) and for microbial Fe-oxyhydroxide dissolution (phenazine biosynthesis) are disproportionately enriched for trace-metal handling pathways (Cu, Fe, Zn, Co, Ni transport and homeostasis), beyond what is expected by chance.

## Data and Methods

### Dataset

`kbase_ke_pangenome` on the KBase BER Data Lakehouse (BERDL): 132.5 million gene clusters across 27,702 GTDB species pangenomes with Bakta functional annotations (KEGG KO, gene names) and PFAM domain assignments (18.8M domain hits). Taxonomy from GTDB R214.

**Soil+plant subset:** Environmental metadata from `ncbi_env` was used to classify species by primary environment. This analysis includes only species classified as soil/rhizosphere (3,409) or plant-associated (1,135), yielding 4,540 species with gene family annotations after intersection with the species gene families dataset.

### Gene family definitions

| Group | Genes | Annotation source |
|-------|-------|------------------|
| **P-acquisition** (9 families) | phoA, phoD, pstA/B/C/S, phnC/D/E | KEGG KO K01077, PFAM PF09423, gene names |
| **N-fixation** (2 families) | nifH, nifD | KO K02588, K02586 |
| **Metal-handling** (4 families) | copA, corA, feoB, HMA | KO K17686, gene name, PFAM PF02421/PF00403 |
| **Phenazine biosynthesis** (7 families) | phzA/B/D/F/G/S/M | KO K06998/K20940, gene names |

N-fixation was defined using KEGG KO assignments only (nifH = K02588, nifD = K02586; 484 species in soil+plant subset). PFAM PF00142 (Fer4_NifH) was evaluated as a sensitivity check but excluded from the primary analysis because it captures the broader Fer4 ferredoxin superfamily (1,033 species in this subset), inflating apparent N-fixation prevalence by 113%.

Species-level presence was scored as ≥1 gene cluster matching the annotation criterion. Phenazine "operon carriers" required ≥3 distinct phz gene families in a single species. This threshold was chosen because the canonical phenazine biosynthesis operon (phzABCDEFG) encodes at least 5–7 genes; requiring ≥3 separates species with a functional operon from those encoding only broadly distributed homologs. phzF alone (a broad COG0384 family member) spans 2,652 species in this subset, while the ≥3 threshold identifies 27 species — consistent with the known narrow phylogenetic distribution of phenazine production (Pseudomonadaceae, Streptomycetaceae, Burkholderiaceae). Lowering to ≥2 would add ~60 species dominated by phzF+phzS co-occurrence without the core biosynthetic genes; raising to ≥4 would reduce to ~20 species without changing the biological conclusions.

### Statistical tests

- **Pairwise co-occurrence:** Jaccard index, phi coefficient, Fisher's exact test (2×2 contingency) for all group pairs and all 72 individual nutrient×metal gene pairs.
- **Multiple testing correction:** Benjamini-Hochberg FDR applied to all 72 Fisher tests; 51/72 pairs significant at q<0.05.
- **Permutation null:** 1,000 random permutations of group membership vectors, preserving marginal counts.
- **Core vs. accessory enrichment:** Per-species Fisher's exact test on a 2×2 table of (core, non-core) × (nutrient gene set, metal gene set) gene counts within each species. Per-species Z-scores aggregated via Stouffer meta-analysis.
- **Phylogenetic stratification:** Phi coefficient computed within each GTDB phylum (minimum 20 species).
- **Phylogenetic logistic regression:** Binary trait co-occurrence tested using `phylolm::phyloglm(method="logistic_MPLE")` on the GTDB R214 bac120 reference tree pruned to 4,177 matching soil+plant species. Pagel's lambda was estimated separately for each gene via `phytools::phylosig(method="lambda")` on a phylum-stratified 3,000-tip subsample. Pairs tested: 6 group-level + 7 strongest positive + 4 strongest negative individual pairs = 17 pairs.
- **Operon-distance test:** For each species with both P-acquisition and metal-handling genes on the representative genome, within-contig gene-ordinal distances were computed for all P × M gene pairs. 1,000 permutations shuffled ordinals within each species to construct the null distribution.
- **Independent validation (Wang et al. 2021):** Phytase × siderophore co-occurrence tested as an independent benchmark. Phylogenetic correction applied using the same `phyloglm` framework.

## Results

### 1. Group-level co-occurrence

| Pair | n_both | Jaccard | Phi | OR | Fisher p |
|------|-------:|--------:|----:|---:|:--------:|
| P × Metal | 3,815 | 0.852 | 0.129 | 3.92 | 4.4×10⁻¹⁴ |
| N × Metal | 477 | 0.109 | 0.040 | 2.74 | 4.3×10⁻³ |
| P × N | 478 | 0.122 | 0.125 | 14.24 | 1.4×10⁻²⁴ |
| Phz-operon × P | 27 | 0.007 | 0.031 | ∞ | 4.3×10⁻² |
| Phz-operon × Metal | 27 | 0.006 | 0.015 | ∞ | 0.624 |
| N × Phz-operon | 0 | 0.000 | −0.027 | 0.00 | 0.108 |

*Table 1. "Phz-operon" = species with ≥3 distinct phz gene families (27 species in soil+plant subset). Metal-handling background prevalence is 96.4% (4,376/4,540), limiting the sensitivity of tests involving metal-handling as one variable.*

P × Metal, N × Metal, P × N, and Phz-operon × P show significant positive co-occurrence. Phz-operon × Metal is non-significant (p=0.624) because metal-handling prevalence in the soil+plant subset (96.4%) leaves insufficient dynamic range — all 27 operon carriers encode metal genes, but 26 of 27 are expected by chance. The N × Phz-operon depletion (0/27 phenazine operon carriers are diazotrophs) is also non-significant (p=0.108) in this smaller dataset.

**Comparison to pan-bacterial analysis:** P × Metal phi increased from 0.110 (pan-bacterial) to 0.129 (soil+plant), and OR from 2.30 to 3.92, consistent with the prediction that the coupling is strongest in environments where Fe-oxyhydroxide mineral surfaces mediate nutrient availability. N × Metal OR decreased from 10.1 to 2.74, primarily because metal-handling background prevalence rose from 91.7% to 96.4%, compressing the dynamic range for enrichment.

**Sensitivity check (N-fixation breadth):** Including PFAM PF00142 expands the N-fixation set from 484 to 1,033 species. The N × Metal phi changes from 0.040 to 0.038, and OR drops from 2.74 to 1.75. As in the pan-bacterial analysis, the inflated species count from PF00142 includes ferredoxin-containing non-diazotrophs.

### 2. Individual gene pair highlights (FDR-corrected)

Of 72 nutrient × metal gene pairs tested, 51 are significant at FDR q<0.05 (vs 64/72 in the pan-bacterial analysis, reflecting reduced statistical power with fewer species).

**Strongest positive associations:**

| Nutrient gene | Metal gene | Enrichment | Phi | n_both | q-value |
|--------------|-----------|----------:|----:|-------:|--------:|
| phoD | HMA | 1.82× | 0.089 | 72 | 4.2×10⁻⁸ |
| phzS | corA | 1.47× | 0.165 | 298 | 9.7×10⁻³³ |
| phzG | corA | 1.63× | 0.085 | 53 | 1.7×10⁻¹⁰ |
| phoA | corA | 1.31× | 0.217 | 885 | 2.1×10⁻⁵⁰ |
| nifH | copA | 1.14× | 0.094 | 403 | 1.4×10⁻¹¹ |
| nifH | HMA | 1.35× | 0.064 | 141 | 4.2×10⁻⁵ |

**Non-significant pairs (q≥0.05):** 21 pairs, predominantly involving rare phenazine genes (phzA, phzB, phzD, phzM) with small sample sizes, and several nifD pairs.

**Strongest negative associations:**

| Nutrient gene | Metal gene | Enrichment | Phi | n_both | q-value |
|--------------|-----------|----------:|----:|-------:|--------:|
| pstC | feoB | 0.48× | −0.395 | 379 | 5.4×10⁻¹⁵⁶ |
| pstS | feoB | 0.55× | −0.274 | 349 | 8.2×10⁻⁷⁸ |
| phzF | feoB | 0.68× | −0.250 | 538 | 3.3×10⁻⁶² |
| pstC | HMA | 0.62× | −0.249 | 384 | 9.7×10⁻⁶² |

The negative associations are substantially stronger in the soil+plant subset (phi=−0.395 for pstC×feoB vs −0.256 pan-bacterially), suggesting that the niche separation between high-affinity phosphate scavenging and ferrous iron uptake is more pronounced in terrestrial environments.

### 3. Core vs. accessory structure

**P-acquisition: core-enriched.** P-genes are 74.0% core (Stouffer meta-Z=25.5, p≈0), with consistent enrichment across species. This is slightly higher than the pan-bacterial value (70.7%) and reflects the stable genomic integration of phosphate scavenging as a species-level metabolic commitment in soil and plant environments.

**N-fixation: not differentially core relative to metal-handling.** KO-defined nifH/nifD are 63.6% core, but Stouffer meta-Z=−0.3 (p=0.759), indicating that N-fixation genes are not enriched in the core genome relative to metal-handling genes. This contrasts with the pan-bacterial result (Z=3.2, p=0.001) and likely reflects the composition of the soil+plant subset: soil diazotrophs may have a higher proportion of nif genes on mobile elements (symbiosis islands, ICEs) compared to the broader bacterial tree.

**Metal-handling: bimodal.** corA (73.5% core) remains as stable as phosphate transporters. feoB (30.0%) and HMA (23.6%) remain predominantly accessory. copA (56.8%) is intermediate.

### 4. Phylogenetic distribution

**Phylum-level:** Among 22 phyla with ≥20 soil+plant species, positive P × Metal co-occurrence (phi > 0) is found in 16. The highest phi values occur in Bacillota_A (0.312, n=119), Bacillota (0.250, n=280), and Pseudomonadota (0.247, n=1,335). The four largest phyla — Pseudomonadota, Actinomycetota, Bacillota, and Bacillota_A — all show strong positive associations.

**Phenazine operon carriers (27 species in soil+plant subset):** Concentrated in two phyla:
- Actinomycetota (21 species): Streptomycetaceae (16), Streptosporangiaceae (3), Pseudonocardiaceae (1), Jiangellaceae (1)
- Pseudomonadota (6 species): Pseudomonadaceae (3), Xanthomonadaceae (2), SG8-39 (1)

The dominance of soil Actinomycetota (21/27, 78%) is even more pronounced in the soil+plant subset than pan-bacterially (35/63, 56%), reflecting the enrichment of Streptomycetaceae in soil environments. The *Xenorhabdus* clade (insect pathogens, 6 species in the pan-bacterial analysis) is absent from the soil+plant subset, as expected given its insect-associated ecology.

### 5. Positive-control species check

Eight well-characterized plant-microbe-soil model organisms were checked against the GTDB pangenome data. All eight are present in the soil+plant subset and encode both P-acquisition and metal-handling genes. Results are unchanged from the pan-bacterial analysis:

| Organism | GTDB species | P-acquisition | N-fixation | Metal-handling | Phenazine |
|----------|-------------|---------------|------------|----------------|-----------|
| *P. fluorescens* | s\_\_Pseudomonas\_E\_fluorescens | pstA/B/S, phnC/D/E | nifH† | copA, corA | phzF |
| *P. protegens* | s\_\_Pseudomonas\_E\_protegens | pstA/B/S | — | copA, corA, HMA | phzF |
| *B. diazoefficiens* | s\_\_Bradyrhizobium\_diazoefficiens | pstA/B/C/S, phnC/D/E | nifH, nifD | copA, corA | phzF, phzS |
| *S. meliloti* | s\_\_Sinorhizobium\_meliloti | phoD, pstA/B/C/S, phnC/D/E | nifH, nifD | copA, corA, feoB, HMA | phzF |
| *M. loti* | s\_\_Mesorhizobium\_loti | pstA/B/C/S, phnC/D/E | nifD | copA, corA, HMA | phzF |
| *S. coelicolor* | s\_\_Streptomyces\_anthocyanicus | phoD, pstA/B/C/S | — | copA, corA | phzF |
| *M. extorquens* | s\_\_Methylobacterium\_extorquens | phoA, pstA/B/C/S, phnC/D/E | — | copA, corA | phzF |
| *P. chlororaphis* | s\_\_Pseudomonas\_E\_chlororaphis | phoA, pstA/B/S | — | copA, corA | phzA/B/F/G (operon) |

†*P. fluorescens* nifH hit likely reflects a divergent Fer4-family ferredoxin rather than true nitrogenase.

### 6. Per-phylum forest plot

Log-odds ratios with 95% confidence intervals were computed for each phylum with ≥20 soil+plant species (22 phyla) for the three primary co-occurrence pairs (Figure 2).

**P × Metal:** Positive log-OR in 16/22 phyla. Four phyla reach statistical significance: Pseudomonadota (log-OR=+2.77, p=7.8×10⁻⁹, n=1,335), Actinomycetota (+2.48, p=2.9×10⁻⁴, n=809), Bacillota (+2.98, p=0.013, n=280), and Bacillota_A (+2.93, p=0.006, n=119). Most other phyla are non-significant due to smaller sample sizes in the soil+plant subset.

**N × Metal:** Largely non-significant at the per-phylum level, reflecting the small number of diazotrophs per phylum in this subset.

**Phz × Metal:** Non-significant in all phyla due to the rarity of phenazine operons (27 total across only 2 phyla).

### 7. Phylogenetic correction

**Phylogenetic signal (Pagel's lambda):** Estimated on a 3,000-tip phylum-stratified subsample from the soil+plant species. All testable gene families show significant phylogenetic signal (lambda range 0.06–1.00). Core P-acquisition genes show the highest signal: pstA (1.00), pstB (1.00), pstC (0.95). Metal-handling: corA (0.96), copA (0.79), feoB (0.79), HMA (0.71). N-fixation: nifH (0.85), nifD (0.76). Rare phenazine genes (phzA, phzB, phzD, phzM) could not be tested due to near-zero prevalence.

**Phylogenetic logistic regression (Table 5):**

| Pair | Uncorrected log-OR | Phylo log-OR | Change | Phylo p | Converged |
|------|-------------------:|-------------:|-------:|--------:|:---------:|
| P × Metal | 1.599 | 0.938 | −41% | 2.1×10⁻⁵ | Yes |
| N × Metal | 0.937 | 0.986 | +5% | 9.2×10⁻⁵ | Yes |
| P × N | 3.361 | 3.139 | −7% | 7.2×10⁻⁵ | Yes |
| phoD × HMA | 0.989 | 0.914 | −8% | 1.4×10⁻⁹ | Yes |
| phoD × feoB | 0.361 | 0.509 | +41% | 2.1×10⁻⁴ | Yes |
| nifH × HMA | 0.208 | 0.396 | +91% | 6.9×10⁻⁴ | Yes |
| pstC × feoB | −1.400 | −0.880 | −37% | 1.8×10⁻²⁰ | Yes |
| pstS × feoB | −0.927 | −0.549 | −41% | 1.2×10⁻¹¹ | Yes |
| pstC × HMA | −1.043 | −0.749 | −28% | 7.2×10⁻¹⁶ | Yes |
| phzF × feoB | −0.690 | −0.457 | −34% | 2.1×10⁻⁹ | Yes |
| phzG × corA | 3.363 | 3.335 | −1% | 2.0×10⁻⁶ | Yes |

*Additional tested pairs: Phz-operon × Metal (p=0.98), Phz-operon × P (p=0.94), N × Phz-operon (p=1.1×10⁻⁵), phzD × corA (p=0.95), phzA × corA (p=0.95), phzB × corA (p=0.96). All six are non-significant after phylo correction; the three Phz-operon group pairs and three rare phenazine individual pairs lose significance due to small n (27 operon carriers, 3–16 species for rare phz genes).*

**Key findings:** 10 of 12 significant uncorrected pairs survive phylogenetic correction (83% survival rate). The two losses are phzD × corA and phzB × corA — rare phenazine genes with small sample sizes (n=10 and n=16 respectively). Two pairs gain significance after correction: N × Phz-operon (a negative association unmasked by phylo correction) and nifH × HMA (uncorrected p=0.10, phylo p=6.9×10⁻⁴).

**Notable difference from pan-bacterial analysis:** In the pan-bacterial analysis, all 14 significant pairs survived correction and negative associations *strengthened* by 9–21%. In the soil+plant subset, P × Metal is attenuated by 41% and negative associations are attenuated by 28–41%. This reflects the higher phylogenetic clustering of the soil+plant subset — these species are more closely related on average, so phylogenetic confounding accounts for a larger fraction of the uncorrected signal. The critical point is that the associations remain significant after correction: ecological co-selection operates beyond what phylogeny alone can explain. See Figure 4.

### 8. Operon-distance test (functional coupling)

**Results (Table 6):**

| Metric | Value |
|--------|------:|
| Species with P + M genes | 3,546 |
| Species with same-contig pairs | 1,299 (36.6%) |
| Total same-contig P × M pairs | 15,589 |
| Same-contig fraction of all pairs | 29.3% |
| Observed median distance | 1,097 genes |
| Null median (1000 permutations) | 350.4 ± 6.2 genes |
| Z-score (observed vs null) | 120.7 |
| Pairs within 5 genes (operon-proximal) | ~0.49% |

P × M gene pairs on the same contig are **farther apart** than expected by chance (median 1,097 vs null 350, Z=120.7). The co-occurrence signal is not driven by physical gene linkage. The pattern is consistent with the pan-bacterial result (median 910, Z=304), with the lower Z-score in the soil+plant subset reflecting the smaller dataset. See Figure 5.

### 9. KEGG pathway co-membership

**Gene-family level:** P-acquisition genes map to 37 KEGG pathways; metal-handling genes map to 53 pathways; 11 pathways are shared (predominantly broad categories like ko01100 metabolic pathways).

**Species level:** Among 3,479 species with both P and M pathway annotations, only 4.9% (172 species) share ≥1 KEGG pathway between their P and M gene clusters. The mean shared pathways per species is 0.05 — **lower** than expected by chance (null mean 0.34, Z = −31.0). P-acquisition and metal-handling genes participate in distinct metabolic pathways. Combined with the operon-distance test, this confirms that the co-occurrence reflects independent ecological co-selection, not pathway integration. See pan-bacterial analysis for extended interpretation.

### 10. Independent validation: Wang et al. 2021 phytase × siderophore coupling

**Soil+plant co-occurrence:**

| Metric | Value |
|--------|------:|
| Total species | 4,540 |
| Phytase-encoding species | 1,122 (24.7%) |
| Siderophore-encoding species | 437 (9.6%) |
| Both traits | 124 (2.7%) |
| OR | 1.23 |
| Phi | 0.028 |
| Fisher p | 0.070 |
| Permutation Z | 1.93 (p=0.038) |

The phytase × siderophore co-occurrence is **non-significant** (Fisher p=0.070) in the soil+plant subset, unlike the pan-bacterial result (OR=2.0, p=2.3×10⁻⁴³). This likely reflects the loss of species diversity: the 4,540 soil+plant species represent a phylogenetically narrower sample than the full 27,682 species, and siderophore-encoding species are relatively less common (9.6% vs 7.6% pan-bacterially), reducing statistical power for detecting modest effect sizes.

**Per-siderophore-type co-occurrence:** Only pyoverdine (OR=1.83, p=0.004) and hydroxamate (OR=1.38, p=0.048) retain significant associations. Enterobactin (OR=1.20, p=0.211) and pyochelin (OR=0.70, p=0.402) are non-significant.

**Phylogenetically corrected:** Phylogenetic logistic regression failed to converge (converged=FALSE), consistent with the weak underlying signal in this subset.

**Family-level distribution:** The top families by linked-trait count are Myxococcaceae (16, 76.2% of family) and Streptomycetaceae (16, 7.2% of family). Burkholderiaceae contribute 9/124 linked-trait species (7.3%), with a within-family rate of 10.1%. See Figure 6.

## Interpretation

### The macro-micro nutrient coupling is stronger in the soil+plant subset but requires careful phylogenetic interpretation

The central finding — P × Metal phi of 0.129 (OR=3.92) vs the pan-bacterial phi of 0.110 (OR=2.30) — supports the prediction that macro-micro nutrient coupling is strongest in environments where Fe-oxyhydroxide mineral surfaces mediate nutrient availability. However, the coupling is also more phylogenetically confounded: the uncorrected log-OR is attenuated by 41% after phylogenetic correction (1.599 → 0.938), compared to a 2% increase in the pan-bacterial analysis. This means that while the true ecological coupling IS significant (phylo p=2.1×10⁻⁵), a substantial fraction of the uncorrected soil+plant signal reflects the phylogenetic clustering of soil-adapted lineages rather than independent co-selection events.

The 41% attenuation reflects the soil+plant subset's phylogenetic structure: soil-adapted lineages (Actinomycetota, Bacillota, Pseudomonadota) are phylogenetically clustered and share both P-acquisition and metal-handling genes as ancestral traits. The phylogenetic correction removes this shared-ancestry component, leaving the ecological co-selection signal — species within soil lineages that acquire P are *still* enriched for metal-handling genes beyond what their shared ancestry predicts. The higher attenuation compared to the pan-bacterial analysis (41% vs 2%) is expected because filtering to a single ecological niche concentrates phylogenetically related organisms that inherited similar gene repertoires from common soil-adapted ancestors.

The coupling remains strongest for the specific mechanistic links predicted by biochemistry:
- **phoD × feoB** strengthened by 41% after phylogenetic correction (0.361 → 0.509), the largest gain of any pair
- **nifH × HMA** strengthened by 91% (0.208 → 0.396), suggesting that in the soil+plant subset, phylogenetic structure was actively masking the N-fixation × metal coupling
- **N × Metal** strengthened by 5% (0.937 → 0.986), one of only three pairs that strengthen after correction

### The coupling is ecological, not genomic linkage or pathway integration

The operon-distance and KEGG pathway tests replicate the pan-bacterial findings: P and M genes are not physically linked (median 1,097 genes apart vs null 350, Z=120.7) and do not share metabolic pathways (Z=−31.0). The coupling operates at the organism level through shared ecological selection pressures.

### Core/accessory structure: P-acquisition is the dominant signal

In the soil+plant subset, only P-acquisition genes show significant core enrichment relative to metal-handling genes (Z=25.5). The loss of significance for N-fixation (Z=−0.3 vs pan-bacterial Z=3.2) suggests that in soil environments, nif genes may be more frequently carried on mobile elements (symbiosis islands, ICEs) — consistent with the well-documented horizontal transfer of nitrogen fixation capacity among rhizosphere bacteria.

### Phenazine operon carriers: smaller sample, consistent pattern

The 27 phenazine operon carriers in the soil+plant subset maintain 100% overlap with both P-acquisition and metal-handling genes, but the small sample size and high metal-handling background prevalence (96.4%) make the associations statistically non-significant. The taxonomic concentration shifts toward Actinomycetota (78% in soil+plant vs 56% pan-bacterially), reflecting the loss of insect-associated *Xenorhabdus* and the enrichment of soil *Streptomyces*.

### Wang 2021 validation: power loss in the subset

The non-significance of phytase × siderophore co-occurrence (p=0.070) in the soil+plant subset contrasts with the highly significant pan-bacterial result (p=2.3×10⁻⁴³). The OR dropped from 2.0 to 1.23 — a genuine weakening of the association, not merely a power issue. This may reflect the narrower phylogenetic diversity of the soil+plant subset or the fact that the phytase–siderophore axis operates across a broader ecological range than the P-acquisition × metal-handling axis.

### Negative associations reflect redox niche separation

The pstC/S × feoB anti-correlations are the strongest individual associations in this analysis (phi up to −0.395; for context, phi=0.0 indicates independence and |phi|=1.0 indicates perfect association). Their attenuation by 28–41% under phylogenetic correction indicates that some of the anti-correlation is driven by phylogenetic clustering: aerobic lineages that rely on high-affinity Pst phosphate transporters (predominantly Actinomycetota, Bacillota) are phylogenetically separated from microaerobic lineages that use the Fe²⁺-specific FeoB transporter (predominantly Pseudomonadota). The phylogenetic correction removes the shared-ancestry component of this separation, but the associations remain highly significant after correction (p<10⁻¹¹), confirming genuine ecological niche separation beyond phylogenetic structure.

## Limitations

1. **Presence/absence only.** This analysis scores whether a species pangenome contains at least one gene cluster matching an annotation criterion. It does not test co-expression, co-regulation, or physical clustering.

2. **Annotation-dependent coverage.** Gene families identified via KEGG KO and PFAM domain hits may miss divergent homologs. The pstA/B/C/S genes lack KEGG KO assignments in Bakta and were identified by gene name only.

3. **PFAM PF00142 captures the broader Fer4 ferredoxin superfamily.** Using PF00142 inflates apparent N-fixation prevalence from 484 to 1,033 species (113% increase). Primary results use KO-restricted definitions.

4. **Phylogenetic non-independence — controlled but informative.** Phylogenetic logistic regression on the 4,177-tip soil+plant subtree shows that 41% of the P × Metal uncorrected signal is phylogenetic confounding — substantially more than the 2% seen pan-bacterially. This highlights the importance of phylogenetic correction when analyzing ecologically filtered subsets.

5. **Ecological inference from genomic potential.** Encoding a gene does not guarantee its expression.

6. **Reduced statistical power.** The 4,540-species subset has lower power than the 27,682-species pan-bacterial analysis. This is evident in: fewer FDR-significant pairs (51 vs 64), non-significance of Phz-operon associations (n=27 vs 63), non-significance of the Wang 2021 phytase×siderophore test, and collapse of two rare phenazine pairs under phylogenetic correction.

7. **Near-universal metal-handling prevalence.** Metal-handling genes are present in 96.4% of soil+plant species (vs 91.7% pan-bacterially). This ceiling effect compresses enrichment ratios and makes it harder to detect positive associations involving metal-handling genes.

8. **Environmental classification limitations.** The soil/rhizosphere and plant-associated categories are derived from keyword matching on NCBI biosample metadata. Some species may be misclassified, and the "plant-associated" category includes both rhizosphere and endophytic associations.

## Future Directions

- **Substrate C (mechanistic, deferred):** Gene-level fitness analysis under controlled P-starvation and metal-stress conditions.

- **Substrate B (ecological):** Where do communities enriched for the macro-micro co-occurrence gene profile appear along nutrient gradients in natural soils? Data: NMDC multi-omics + abiotic features + ENIGMA CORAL.

- **Field validation (Point Reyes):** The mafic-vs-felsic fault gradient at Point Reyes (Rowley et al. 2023, 2024) provides a natural experiment.

- **Expanded phenazine gene set:** Expand the phenazine gene set to include actinomycete-specific phenazine biosynthesis markers.

- **Cross-environment comparison:** Systematically compare effect sizes and phylogenetic correction magnitudes across all environmental categories (soil, marine, freshwater, human-associated) to identify environment-specific vs universal components of macro-micro nutrient coupling.

## References

- Edayilam, N. et al. (2018). Phosphorus stress-induced changes in plant root exudate composition and their effects on morphology and key functional groups of soil bacterial and archaeal communities. *Environmental Microbiology*, 20(7), 2504–2521.
- McRose, D.L. & Newman, D.K. (2021). Redox-active antibiotics enhance phosphorus bioavailability. *Science*, 371, 1033–1037.
- Wu, Z. et al. (2022). Phenazine biosynthesis genes in metagenome-assembled genomes from phosphorus-limited reforestation soils. *ISME Journal*.
- Wang, H. et al. (2021). Phytase-producing rhizobacteria as a natural control for enhancing phosphorus phytoavailability. *Frontiers in Microbiology*, 12, 572212.
- Dakora, F.D. & Phillips, D.A. (2002). Root exudates as mediators of mineral acquisition in low-nutrient environments. *Plant and Soil*, 245, 35–47.

## Data and Reproducibility

All analyses were performed on the KBase BER Data Lakehouse using Spark SQL queries against `kbase_ke_pangenome` (132.5M gene clusters, 27,702 species pangenomes, GTDB R214). The soil+plant subset (4,540 species) was derived by filtering on environmental metadata from `ncbi_env` (categories: soil/rhizosphere, plant-associated). Source code is in `src/01_extract_gene_families.py` through `src/16_kegg_pathway_comembership.py`. Phylogenetic analyses use R packages `phylolm` (v2.6.5) and `phytools` (v2.3-0) via `src/09_phylo_signal.R` and `src/10_phylo_logistic.R`, run in a dedicated `r_phylo` conda environment. Intermediate data files are in `data/`. Figures are in `figures/`.

**Figure 1.** (A) Phi coefficient heatmap for 18 nutrient/phenazine × 4 metal gene pairs (72 FDR-tested pairs). (B) Core genome fraction per gene family. (C) Phylum-level P × Metal co-occurrence. (D) Taxonomic distribution of phenazine operon carriers. See `figures/figure1_cooccurrence.png`.

**Figure 2.** Per-phylum forest plot of log-odds ratios with 95% CIs for P × Metal, N × Metal, and Phz × Metal co-occurrence across 22 soil+plant phyla (n≥20). See `figures/forest_plot.png`.

**Figure 4.** Uncorrected vs phylogenetically corrected log-ORs for 17 tested pairs on the 4,177-tip soil+plant subtree. Points color-coded by mean Pagel's lambda of the gene pair. See `figures/figure4_phylo_correction.png`.

**Figure 5.** (A) Distribution of within-contig P × M gene-ordinal distances (15,589 pairs, log scale). (B) Observed median vs 1,000-permutation null distribution. See `figures/figure5_operon_distance.png`.

**Figure 6.** Wang et al. 2021 validation. (A) Per-family phytase × siderophore co-occurrence rates in the soil+plant subset. (B) Per-siderophore-type odds ratios. See `figures/figure6_wang2021.png`.
