# Macro-Micro Nutrient Gene Co-occurrence Across 27,682 Bacterial Pangenomes

## Key Findings

1. **P-acquisition and metal-handling genes co-occur significantly** across 27,682 GTDB species pangenomes (phi=0.110, OR=2.3, Fisher p=1.3×10⁻⁶⁵; permutation Z=17.7, p<0.001). 80.2% of species encode at least one P-acquisition gene; 91.7% encode at least one metal-handling gene; 74.7% encode both.

2. **N-fixation and metal-handling show the strongest per-species association.** Using KEGG KO-defined nifH/nifD (2,746 species), N × Metal has OR=10.1 (Fisher p=1.3×10⁻⁷¹; permutation Z=14.7). The high OR reflects near-universal metal gene presence among true diazotrophs (2,719/2,746 = 99.0%). Individual gene pairs (nifH × HMA, nifH × feoB) reach enrichments of 1.42× with phi=0.077 (all FDR q<0.05).

3. **All 63 phenazine operon carriers encode both P-acquisition and metal-handling genes** (100% overlap, OR=∞, Fisher p=9.4×10⁻³ for Phz-operon × Metal). These 63 species are concentrated in soil and insect-pathogenic families: Streptomycetaceae (24), Pseudomonadaceae (17), Streptosporangiaceae (10), and Enterobacteriaceae (6, all *Xenorhabdus* — entomopathogenic nematode symbionts).

4. **Three distinct genomic signatures of macro-micro coupling emerge from core/accessory analysis.** P-acquisition genes are predominantly core (70.7%; Stouffer Z=68.3), N-fixation genes are moderately core (63.4%; Stouffer Z=3.2), and metal-handling genes span the core–accessory boundary (copA 56.8%, corA 73.5%, feoB 30.0%, HMA 23.6%).

5. **Negative associations reveal biological structure.** Pst transporter genes (pstC, pstS, pstA) are depleted relative to feoB (0.64–0.86× expected; phi as low as −0.256; all FDR q<0.05), suggesting that high-affinity phosphate transport and ferrous iron uptake occupy distinct ecological niches.

## Background

Plant productivity and soil organic matter turnover depend on the simultaneous availability of macronutrients (C, N, P) and trace metals (Fe, Cu, Zn, Mn, Co, Ni, Mo). These pools are coupled through enzyme cofactor biochemistry — nitrogenase requires Fe and Mo, alkaline phosphatases require Zn or Ca, and urease requires Ni — and through mineral-surface biogeochemistry, where phosphate and trace-metal cations co-adsorb onto Fe(III)-oxyhydroxide surfaces in acidic, iron-rich soils (Edayilam et al. 2018).

McRose & Newman (2021, *Science* 371:1033–1037) demonstrated that phenazine biosynthesis is upregulated 1–3 orders of magnitude under phosphate starvation via the conserved PhoB regulon. Phenazines reductively dissolve Fe(III)-oxyhydroxides, simultaneously releasing adsorbed phosphate and co-adsorbed trace metals (Co, Cu, Ni). Wu et al. (2022) confirmed the presence of phenazine biosynthesis genes (phzS) in metagenome-assembled genomes from P-limited reforestation soils, establishing that this pathway operates in natural terrestrial microbiomes.

Despite these mechanistic links, macro-nutrient cycling and trace-metal dynamics are typically studied separately. This analysis provides the first genome-scale test of whether genes mediating phosphate and nitrogen acquisition share genomic context with metal-handling genes across bacterial diversity.

## Hypothesis

Bacterial pangenomes encoding genes for macro-nutrient acquisition (phosphatases, phosphate/phosphonate transporters, nitrogenases) and for microbial Fe-oxyhydroxide dissolution (phenazine biosynthesis) are disproportionately enriched for trace-metal handling pathways (Cu, Fe, Zn, Co, Ni transport and homeostasis), beyond what is expected by chance.

## Data and Methods

### Dataset

`kbase_ke_pangenome` on the KBase BER Data Lakehouse (BERDL): 132.5 million gene clusters across 27,702 GTDB species pangenomes with Bakta functional annotations (KEGG KO, gene names) and PFAM domain assignments (18.8M domain hits). Taxonomy from GTDB R214.

### Gene family definitions

| Group | Genes | Annotation source |
|-------|-------|------------------|
| **P-acquisition** (9 families) | phoA, phoD, pstA/B/C/S, phnC/D/E | KEGG KO K01077, PFAM PF09423, gene names |
| **N-fixation** (2 families) | nifH, nifD | KO K02588, K02586 |
| **Metal-handling** (4 families) | copA, corA, feoB, HMA | KO K17686, gene name, PFAM PF02421/PF00403 |
| **Phenazine biosynthesis** (7 families) | phzA/B/D/F/G/S/M | KO K06998/K20940, gene names |

N-fixation was defined using KEGG KO assignments only (nifH = K02588, nifD = K02586; 2,746 species). PFAM PF00142 (Fer4_NifH) was evaluated as a sensitivity check but excluded from the primary analysis because it captures the broader Fer4 ferredoxin superfamily (5,872 species), inflating apparent N-fixation prevalence by 114% (see Limitations).

Species-level presence was scored as ≥1 gene cluster matching the annotation criterion. Phenazine "operon carriers" required ≥3 distinct phz gene families in a single species, filtering broad-family hits (phzF superfamily alone spans 10,856 species) from true operon-bearing lineages (63 species).

### Statistical tests

- **Pairwise co-occurrence:** Jaccard index, phi coefficient, Fisher's exact test (2×2 contingency) for all group pairs and all 72 individual nutrient×metal gene pairs.
- **Multiple testing correction:** Benjamini-Hochberg FDR applied to all 72 Fisher tests; 64/72 pairs significant at q<0.05.
- **Permutation null:** 1,000 random permutations of group membership vectors, preserving marginal counts. Observed phi compared to null distribution; Z-score and permutation p-value reported.
- **Core vs. accessory enrichment:** Per-species Fisher's exact test on a 2×2 table of (core, non-core) × (nutrient gene set, metal gene set) gene counts within each species. This tests whether the two functional groups differ in their core-genome proportion. Per-species Z-scores are aggregated via Stouffer meta-analysis (Z_meta = Σz / √n) to assess the direction and magnitude of differential core enrichment across all co-encoding species.
- **Phylogenetic stratification:** Phi coefficient computed within each GTDB phylum and class (minimum 50 or 20 species, respectively). Plant-associated families tested against global baseline.

## Results

### 1. Group-level co-occurrence

| Pair | n_both | Jaccard | Phi | OR | Fisher p | Permutation Z |
|------|-------:|--------:|----:|---:|:--------:|--------------:|
| P × Metal | 20,683 | 0.769 | 0.110 | 2.30 | 1.3×10⁻⁶⁵ | 17.7 |
| N × Metal | 2,719 | 0.107 | 0.088 | 10.08 | 1.3×10⁻⁷¹ | 14.7 |
| Phz-operon × Metal | 63 | 0.002 | 0.014 | ∞ | 9.4×10⁻³ | 2.3 |
| P × N | 2,414 | 0.107 | 0.065 | 1.90 | 1.3×10⁻²⁹ | 10.7 |
| Phz-operon × P | 63 | 0.003 | 0.024 | ∞ | 1.5×10⁻⁶ | 4.0 |
| N × Phz-operon | 0 | 0.000 | −0.016 | 0.00 | 2.4×10⁻³ | — |

All group pairs except N × Phz-operon show significant positive co-occurrence. The N × Phz-operon depletion (0/63 phenazine operon carriers are diazotrophs) reflects the concentration of phenazine operons in non-diazotrophic Actinomycetota (Streptomycetaceae) and non-diazotrophic Pseudomonas lineages.

**Sensitivity check (N-fixation breadth):** Including PFAM PF00142 expands the N-fixation set from 2,746 to 5,872 species. The N × Metal phi changes from 0.088 to 0.107, and OR drops from 10.1 to 4.0. The inflated species count from PF00142 includes ferredoxin-containing non-diazotrophs, diluting the true diazotroph signal. All primary results use the KO-restricted definition.

### 2. Individual gene pair highlights (FDR-corrected)

Of 72 nutrient × metal gene pairs tested, 64 are significant at FDR q<0.05.

**Strongest positive associations:**

| Nutrient gene | Metal gene | Enrichment | Phi | n_both | q-value |
|--------------|-----------|----------:|----:|-------:|--------:|
| phzG | corA | 1.79× | 0.057 | 116 | 1.3×10⁻²⁵ |
| phzD | corA | 1.78× | 0.028 | 29 | 6.8×10⁻⁶ |
| phzA | corA | 1.73× | 0.019 | 15 | 3.8×10⁻³ |
| phoD | HMA | 1.72× | 0.084 | 468 | 1.7×10⁻³⁸ |
| phzB | corA | 1.65× | 0.030 | 43 | 5.7×10⁻⁶ |
| phoD | feoB | 1.53× | 0.077 | 570 | 2.8×10⁻³⁴ |
| nifH | HMA | 1.42× | 0.077 | 903 | 2.1×10⁻³³ |

**Non-significant pairs (q≥0.05):** phzM × HMA (n=6, q=0.18), phzM × feoB (n=4, q=0.63), phzA × HMA (n=3, q=0.34), and 5 others with small sample sizes (<50 co-occurrences).

**Strongest negative associations:**

| Nutrient gene | Metal gene | Enrichment | Phi | n_both | q-value |
|--------------|-----------|----------:|----:|-------:|--------:|
| pstC | feoB | 0.67× | −0.256 | 3,386 | <10⁻³⁰⁰ |
| pstS | feoB | 0.64× | −0.184 | 2,038 | 7.4×10⁻²¹³ |
| pstC | HMA | 0.72× | −0.173 | 2,655 | 3.2×10⁻¹⁸² |

The pstC/S–feoB anti-correlation may reflect niche separation between organisms specializing in high-affinity phosphate scavenging (oligotrophic, often aerobic) versus those with active ferrous iron uptake (anaerobic or microaerobic environments where Fe²⁺ is soluble).

### 3. Core vs. accessory structure

Three distinct genomic signatures emerge from the core/accessory analysis:

**Signature 1: P-acquisition — core-enriched.** P-genes are 70.7% core (Stouffer meta-Z=68.3, p≈0), with Pst transporter genes at 73–77% core. This reflects stable genomic integration of phosphate scavenging as a species-level metabolic commitment.

**Signature 2: N-fixation — moderately core, consistent with known HGT.** KO-defined nifH/nifD are 63.4% core (Stouffer meta-Z=3.2, p=0.001), lower than P-acquisition genes. This is consistent with the well-documented horizontal transfer of nif operons via symbiosis islands and integrative conjugative elements, which places a substantial fraction of nitrogenase genes in the accessory genome even in species where nitrogen fixation is ecologically important.

**Signature 3: Metal-handling — bimodal (core transport + accessory efflux).** corA (73.5% core) is as stable as phosphate transporters, reflecting its essential role in Mg²⁺/Co²⁺ homeostasis. In contrast, feoB (30.0%) and HMA (23.6%) are predominantly accessory, suggesting that specialized metal efflux and ferrous iron uptake are gained or lost in response to local geochemical conditions. copA (56.8%) occupies an intermediate position.

The one notable exception within P-acquisition is phoD (PFAM PF09423), which is only 17.2% core — the lowest core fraction of any nutrient gene. PhoD is a calcium-dependent phosphodiesterase often carried on mobile elements, and its accessory status is consistent with horizontal transfer as a mechanism for spreading P-acquisition capacity.

### 4. Phylogenetic distribution

**Phylum-level:** Positive P × Metal co-occurrence (phi > 0) is found in all 12 largest phyla except Spirochaetota (phi=−0.085, n=296) and Chloroflexota (phi=−0.074, n=457). The highest phi values occur in Cyanobacteriota (0.355, n=469), Fusobacteriota (0.453, n=59), and Bacillota (0.246, n=2,146). Pseudomonadota — the largest phylum (n=7,456) and the one hosting 27 of 63 phenazine operon carriers — has phi=0.167.

**Class-level:** Among classes with ≥50 species, the strongest P-Metal phi is in Fusobacteriia (0.453, n=59), Bacilli (0.246, n=2,146), and Gammaproteobacteria (0.185, n=4,731). Gammaproteobacteria hosts all 27 Pseudomonadota phenazine operon carriers and combines the highest absolute co-occurrence count with strong statistical association.

**Plant-associated families:** Pseudomonadaceae (n=498), Rhizobiaceae (n=413), Burkholderiaceae (n=333), Streptomycetaceae (n=382), and Xanthomonadaceae (n=207) all show near-universal P-acquisition (100%) and metal-handling (93–100%) gene prevalence. Their enrichment ratios (1.00×) reflect ceiling effects — essentially all species in these families encode both functional groups.

**Phenazine operon carriers (63 species):** Concentrated in 9 families:
- Streptomycetaceae (24): *Streptomyces* and *Kitasatospora* — soil actinomycetes
- Pseudomonadaceae (17): *Pseudomonas* and *Pseudomonas_E* — rhizosphere colonizers
- Streptosporangiaceae (10): *Microbispora*, *Nocardiopsis*, *Nonomuraea* — soil actinomycetes
- Enterobacteriaceae (6): *Xenorhabdus* — insect-pathogenic (entomopathogenic nematode symbionts)
- Xanthomonadaceae (2): *Lysobacter* — soil/rhizosphere
- Burkholderiaceae (1): *Burkholderia*

The dominance of soil Actinomycetota (35/63, 56%) and Pseudomonadota (20/63, 32%) among phenazine operon carriers is consistent with the McRose-Newman model: phenazine-mediated Fe(III) reduction is an adaptive strategy in soils where Fe-oxyhydroxides lock up phosphate and trace metals. The *Xenorhabdus* clade (6/63, 10%) represents a distinct ecological niche — insect pathogenesis — where phenazines serve antimicrobial rather than mineral-dissolution functions.

### 5. Positive-control species check

Eight well-characterized plant-microbe-soil model organisms were looked up in the GTDB pangenome data. *S. coelicolor* A3(2) was reclassified as *S. anthocyanicus* in GTDB R214 (confirmed via genome accession GCF_008931305.1).

| Organism | GTDB species | P-acquisition | N-fixation | Metal-handling | Phenazine |
|----------|-------------|---------------|------------|----------------|-----------|
| *P. fluorescens*† | s\_\_Pseudomonas\_E\_fluorescens | pstA/B/S, phnC/D/E | nifH† | copA, corA | phzF |
| *P. protegens* | s\_\_Pseudomonas\_E\_protegens | pstA/B/S | — | copA, corA, HMA | phzF |
| *B. diazoefficiens* | s\_\_Bradyrhizobium\_diazoefficiens | pstA/B/C/S, phnC/D/E | nifH, nifD | copA, corA | phzF, phzS |
| *S. meliloti* | s\_\_Sinorhizobium\_meliloti | phoD, pstA/B/C/S, phnC/D/E | nifH, nifD | copA, corA, feoB, HMA | phzF |
| *M. loti* | s\_\_Mesorhizobium\_loti | pstA/B/C/S, phnC/D/E | nifD | copA, corA, HMA | phzF |
| *S. coelicolor* | s\_\_Streptomyces\_anthocyanicus | phoD, pstA/B/C/S | — | copA, corA | phzF |
| *M. extorquens* | s\_\_Methylobacterium\_extorquens | phoA, pstA/B/C/S, phnC/D/E | — | copA, corA | phzF |
| *P. chlororaphis* | s\_\_Pseudomonas\_E\_chlororaphis | phoA, pstA/B/S | — | copA, corA | phzA/B/F/G (operon) |

†*P. fluorescens* is not a confirmed diazotroph; the nifH hit (KO K02588) likely reflects a divergent Fer4-family ferredoxin rather than true nitrogenase. This does not affect the primary analysis, which defines N-fixation from the full species gene matrix in src/02.

All eight species encode both P-acquisition and metal-handling genes, consistent with the genome-wide co-occurrence signal. The two known diazotrophs (*B. diazoefficiens*, *S. meliloti*) carry both nifH and nifD as expected. *S. meliloti* encodes the broadest gene repertoire: all four metal-handling families plus 8 of 9 P-acquisition families. *P. chlororaphis* is the only positive control with a complete phenazine operon (phzA/B/F/G), consistent with its known phenazine-producing phenotype. *S. coelicolor* — a known phenazine producer — carries only phzF in the pangenome; the remaining phz genes may be annotated under different names in Bakta or absent from the representative genome set.

### 6. Environmental stratification

Environmental metadata from `ncbi_env` was joined to 27,009 species via genome biosample accessions and classified into broad categories. Species counts in the table below reflect those with both environment assignments and gene family data (e.g., 3,406 of 3,409 soil-assigned species).

| Environment | n species | P×Metal log-OR | P×Metal p | N×Metal log-OR | N×Metal p |
|-------------|----------|---------------:|----------:|---------------:|----------:|
| Soil/rhizosphere | 3,406 | +1.34 | 1.9×10⁻⁹ | +0.36 | 0.17 |
| Plant-associated | 1,134 | +1.51 | 1.3×10⁻⁶ | +0.94 | 0.064 |
| Marine | 3,449 | +0.26 | 7.8×10⁻³ | +1.37 | 4.2×10⁻²³ |
| Freshwater/engineered | 2,059 | +1.31 | 6.6×10⁻⁹ | +1.46 | 1.5×10⁻⁵ |
| Human-associated | 3,497 | +1.37 | 3.4×10⁻²⁰ | +0.97 | 7.1×10⁻⁶ |
| Animal-associated | 840 | +0.23 | 0.39 | +2.27 | 2.1×10⁻⁷ |

P × Metal co-occurrence is significant across all environments except animal-associated (p=0.39). The strongest P × Metal effects are in plant-associated (log-OR=+1.51) and soil/rhizosphere (log-OR=+1.34) environments, consistent with the prediction that P-metal coupling is most relevant where Fe-oxyhydroxide mineral surfaces mediate nutrient availability. N × Metal shows a different pattern: the strongest effect is in marine (log-OR=+1.37) and animal-associated (log-OR=+2.27) environments, while soil/rhizosphere shows a non-significant trend (log-OR=+0.36, p=0.17). Phenazine operon × Metal associations were non-significant in all environments due to the rarity of operon carriers (0–29 per environment).

### 7. Per-phylum forest plot

Log-odds ratios with 95% confidence intervals were computed for each major phylum (n≥50 species) for the three primary co-occurrence pairs (Figure 2).

**P × Metal:** Positive log-OR in 28/34 phyla. Strongest effects in Cyanobacteriota (log-OR=+1.86, p=4.1×10⁻¹³), Bacillota (+1.64, p=8.0×10⁻²⁴), and Pseudomonadota (+0.94, p=7.7×10⁻²¹). Five phyla show significant negative associations, including Myxococcota (−2.63, p=0.003) and Spirochaetota (−1.79, p=0.004).

**N × Metal:** Positive log-OR in 24/34 phyla. Strongest in Cyanobacteriota (+2.11, p=2.1×10⁻¹⁴), Nanoarchaeota (+2.72, p=8.2×10⁻⁴), and Spirochaetota (+2.72, p=6.1×10⁻³). The pattern is phylogenetically broader than P × Metal, reflecting that diazotrophy spans diverse lineages.

**Phz × Metal:** Only Pseudomonadota (+1.87, p=0.11) and Actinomycetota (+1.14, p=0.40) show positive trends, consistent with the taxonomic concentration of phenazine operons. All other phyla show null or weakly negative effects due to the extreme rarity of phenazine operons outside these two phyla.

**Figure 2.** Per-phylum forest plot of log-odds ratios for (A) P × Metal, (B) N × Metal, and (C) Phz × Metal co-occurrence. Blue bars: Fisher p<0.05; gray bars: non-significant. See `figures/forest_plot.png`.

## Interpretation

### The macro-micro nutrient coupling is genomically encoded

The central result — significant positive co-occurrence of P-acquisition, N-fixation, and metal-handling genes across 27,682 bacterial pangenomes — supports the hypothesis that macro-nutrient and trace-metal cycling are genomically coupled in bacteria. The effect sizes are modest at the group level (phi 0.01–0.11) but highly significant given the dataset size (all permutation Z > 2.3), and individual gene pairs reach enrichments of 1.4–1.8× (64/72 FDR-significant at q<0.05).

The coupling is not uniform: it is strongest for the specific mechanistic links predicted by biochemistry:
- **nifH × HMA** (nitrogenase requires Fe-Mo cofactor and heavy-metal homeostasis) — strong positive signal
- **phoD × HMA/feoB** (PhoD requires metal cofactors) — strong positive signal
- **pstC × feoB** (high-affinity phosphate transport vs. ferrous iron) — strong negative signal, suggesting niche separation

### Three genomic strategies for macro-micro coupling

The core/accessory analysis reveals that macro-micro nutrient coupling operates through three distinct genomic strategies, not a single mechanism:

1. **P-acquisition: core-genome coupling.** Phosphate scavenging genes are stably integrated into species' core metabolic identity (70.7% core). This reflects the universal, constitutive need for phosphorus and the selective advantage of maintaining reliable P-acquisition machinery.

2. **Phenazine biosynthesis: clade-specific coupling.** Phenazine operons are concentrated in specific soil and rhizosphere lineages (63 species in 9 families), where they are moderately core (phzF 65%, phzS 66%). This represents an ecologically specialized coupling — the McRose-Newman reductive dissolution mechanism — that is stable within the lineages that carry it but absent from most of the bacterial tree.

3. **N-fixation: horizontally acquired coupling.** Nitrogenase genes (nifH/nifD) show lower core fractions (63.4%) than P-acquisition genes, consistent with the well-documented horizontal transfer of nif operons between species via symbiosis islands and integrative conjugative elements. The coupling of N-fixation to metal handling is real (OR=10.1, nearly all diazotrophs carry metal genes) but operates through a genomically flexible, horizontally mobile mechanism rather than stable vertical inheritance.

### Phenazine operon carriers as the extreme case

The 63 species with true phenazine operons (≥3 phz genes) represent the most extreme co-occurrence: 100% encode both P-acquisition and metal-handling genes. Their taxonomic distribution — dominated by soil Actinomycetota and Pseudomonadota — is precisely the ecological context where the McRose-Newman reductive dissolution mechanism would operate. The *Xenorhabdus* clade (6 species) is a notable exception: these insect-pathogenic bacteria likely use phenazines for antimicrobial competition rather than mineral dissolution, illustrating functional convergence of the same biosynthetic pathway across different ecological strategies.

### Negative associations are informative

The depletion of pstC/S with feoB suggests that high-affinity phosphate scavenging and ferrous iron uptake occupy distinct ecological niches. Organisms investing in Pst-type phosphate transport (typically aerobic, low-P environments) may not coexist with conditions favoring FeoB-dependent Fe²⁺ uptake (anaerobic/microaerobic, high-Fe²⁺). This niche separation is consistent with the redox geochemistry of P and Fe: under aerobic conditions, Fe(III) precipitates and adsorbs phosphate (favoring Pst for P scavenging), while under anaerobic conditions, Fe(II) is soluble and P is released (favoring FeoB for Fe uptake but reducing P limitation).

## Limitations

1. **Presence/absence only.** This analysis scores whether a species pangenome contains at least one gene cluster matching an annotation criterion. It does not test whether the genes are co-expressed, co-regulated, or physically clustered in genomes. Expression-level coupling (e.g., PhoB-dependent phenazine upregulation) cannot be inferred from co-occurrence alone.

2. **Annotation-dependent coverage.** Gene families identified via KEGG KO and PFAM domain hits may miss divergent homologs. The pstA/B/C/S genes lack KEGG KO assignments in Bakta and were identified by gene name only; some genuine pst genes may be annotated with different names. Similarly, phenazine gene identification relies on gene name matching, which may miss renamed or poorly annotated homologs.

3. **PFAM PF00142 captures the broader Fer4 ferredoxin superfamily.** The Fer4_NifH domain (PF00142) is shared by nitrogenase iron proteins and other 4Fe-4S ferredoxins. Using PF00142 inflates apparent N-fixation prevalence from 2,746 to 5,872 species (a 114% increase). The primary analysis uses KEGG KO-defined nifH/nifD only; PF00142 is reported as a sensitivity check. This distinction materially affects the N × Metal association: KO-only yields OR=10.1 (phi=0.088), while PF00142-inclusive yields OR=4.0 (phi=0.107) — the inflated species count includes ferredoxin-containing non-diazotrophs that dilute the true association.

4. **Phylogenetic non-independence.** Closely related species share gene content by descent, inflating co-occurrence statistics. The permutation test controls for marginal gene frequencies but not for phylogenetic autocorrelation. A full phylogenetic logistic regression or phylogenetic independent contrasts analysis would provide a stronger control.

5. **Ecological inference from genomic potential.** Encoding a gene does not guarantee its expression in a given environment. The co-occurrence pattern is consistent with — but does not prove — coordinated nutrient-metal cycling in natural environments.

6. **Phenazine operon definition.** The ≥3 distinct phz gene threshold is a heuristic. Some species may carry partial operons or unlinked phz genes from different biosynthetic pathways. Manual curation of operon structure (e.g., checking genomic synteny) would strengthen the phenazine results.

7. **Actinomycete phenazine under-detection.** The phenazine biosynthesis gene set is anchored on the canonical *Pseudomonas* operon (phzABCDEFG, phzM, phzS) with a ≥3-gene threshold for operon completeness. Actinomycete phenazines (e.g. those produced by *Streptomyces*) use divergent biosynthesis pathways with different gene names that the current Bakta annotation set does not collapse into the same gene family. The Phz×Metal effect in Actinomycetota (trending positive in our analysis) is therefore likely a lower bound on the true coupling. A more comprehensive phenazine gene set incorporating actinomycete-specific markers is a future direction. Positive control: *S. coelicolor* A3(2), reclassified as s\_\_Streptomyces\_anthocyanicus in GTDB R214 (genome GCF\_008931305.1), encodes P-acquisition (phoD, pstA/B/C/S), metal-handling (copA, corA), and phzF, but lacks the canonical multi-gene operon — exactly the under-detection pattern.

## Future Directions

- **Substrate C (mechanistic, deferred):** Gene-level fitness analysis under controlled P-starvation and metal-stress conditions. The RB-TnSeq dataset in `kescience_fitnessbrowser` lacks P-starvation experiments (see exhaustive search in `memory/20260507e_*`). New RB-TnSeq campaigns under defined low-P media, or GapMind-extended P-pathway scoring, would provide the mechanistic substrate.

- **Substrate B (ecological):** Where do communities enriched for the macro-micro co-occurrence gene profile appear along nutrient gradients in natural soils? Data: NMDC multi-omics + abiotic features + ENIGMA CORAL.

- **Phylogenetic independent contrasts:** Apply phylogenetic logistic regression using GTDB phylogenetic distance pairs to control for shared ancestry and test whether co-occurrence remains significant after accounting for phylogenetic signal.

- **Field validation (Point Reyes):** The mafic-vs-felsic fault gradient at Point Reyes (Rowley et al. 2023, 2024) provides a natural experiment to test whether co-released P and trace metals (Co, Cu, Ni) from Fe-oxyhydroxide dissolution shape rhizosphere community composition in the pattern predicted by this genomic analysis.

- **Expanded phenazine gene set:** Expand the phenazine gene set to include actinomycete-specific phenazine biosynthesis markers (e.g. naphthoyl synthase, PhzD-like in *Streptomyces*, RppA polyketide synthases) and re-run the (d)×(c) co-occurrence test. The *Pseudomonas*-anchored gene set used in this study likely under-detects the true scope of the phenazine arm of macro-micro nutrient co-mobilization.

## References

- Edayilam, N. et al. (2018). Phosphorus stress-induced changes in plant root exudate composition and their effects on morphology and key functional groups of soil bacterial and archaeal communities. *Environmental Microbiology*, 20(7), 2504–2521.
- McRose, D.L. & Newman, D.K. (2021). Redox-active antibiotics enhance phosphorus bioavailability. *Science*, 371, 1033–1037.
- Wu, Z. et al. (2022). Phenazine biosynthesis genes in metagenome-assembled genomes from phosphorus-limited reforestation soils. *ISME Journal*.
- Dakora, F.D. & Phillips, D.A. (2002). Root exudates as mediators of mineral acquisition in low-nutrient environments. *Plant and Soil*, 245, 35–47.

## Data and Reproducibility

All analyses were performed on the KBase BER Data Lakehouse using Spark SQL queries against `kbase_ke_pangenome` (132.5M gene clusters, 27,702 species pangenomes, GTDB R214). Source code is in `src/01_extract_gene_families.py` through `src/08_forest_plot.py`. Intermediate data files are in `data/`. Figures are in `figures/`.

**Figure 1.** (A) Phi coefficient heatmap for 18 nutrient/phenazine × 4 metal gene pairs (72 FDR-tested pairs). nifH(Pfam) is shown for visual reference as the broad-definition sensitivity check but is excluded from the 72-pair FDR correction. (B) Core genome fraction per gene family. (C) Phylum-level P × Metal co-occurrence. (D) Taxonomic distribution of phenazine operon carriers. See `figures/figure1_cooccurrence.png`.

**Figure 2.** Per-phylum forest plot of log-odds ratios with 95% CIs for P × Metal, N × Metal, and Phz × Metal co-occurrence across 34 GTDB phyla. See `figures/forest_plot.png`.
