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
- **Core vs. accessory enrichment:** Fisher's exact test per species comparing core/non-core fractions between nutrient and metal gene sets. Stouffer meta-analysis across species (Z-scores combined as Z_meta = Σz / √n).
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

## Future Directions

- **Substrate C (mechanistic, deferred):** Gene-level fitness analysis under controlled P-starvation and metal-stress conditions. The RB-TnSeq dataset in `kescience_fitnessbrowser` lacks P-starvation experiments (see exhaustive search in `memory/20260507e_*`). New RB-TnSeq campaigns under defined low-P media, or GapMind-extended P-pathway scoring, would provide the mechanistic substrate.

- **Substrate B (ecological):** Where do communities enriched for the macro-micro co-occurrence gene profile appear along nutrient gradients in natural soils? Data: NMDC multi-omics + abiotic features + ENIGMA CORAL.

- **Phylogenetic independent contrasts:** Apply phylogenetic logistic regression using GTDB phylogenetic distance pairs to control for shared ancestry and test whether co-occurrence remains significant after accounting for phylogenetic signal.

- **Field validation (Point Reyes):** The mafic-vs-felsic fault gradient at Point Reyes (Rowley et al. 2023, 2024) provides a natural experiment to test whether co-released P and trace metals (Co, Cu, Ni) from Fe-oxyhydroxide dissolution shape rhizosphere community composition in the pattern predicted by this genomic analysis.

## References

- Edayilam, N. et al. (2018). Phosphorus stress-induced changes in plant root exudate composition and their effects on morphology and key functional groups of soil bacterial and archaeal communities. *Environmental Microbiology*, 20(7), 2504–2521.
- McRose, D.L. & Newman, D.K. (2021). Redox-active antibiotics enhance phosphorus bioavailability. *Science*, 371, 1033–1037.
- Wu, Z. et al. (2022). Phenazine biosynthesis genes in metagenome-assembled genomes from phosphorus-limited reforestation soils. *ISME Journal*.
- Dakora, F.D. & Phillips, D.A. (2002). Root exudates as mediators of mineral acquisition in low-nutrient environments. *Plant and Soil*, 245, 35–47.

## Data and Reproducibility

All analyses were performed on the KBase BER Data Lakehouse using Spark SQL queries against `kbase_ke_pangenome` (132.5M gene clusters, 27,702 species pangenomes, GTDB R214). Source code is in `src/01_extract_gene_families.py` through `src/05_figure.py`. Intermediate data files are in `data/`. The multi-panel figure is at `figures/figure1_cooccurrence.png`.

**Figure 1.** (A) Phi coefficient heatmap for 72 nutrient × metal gene pairs (FDR-significant pairs shown). (B) Core genome fraction per gene family. (C) Phylum-level P × Metal co-occurrence. (D) Taxonomic distribution of phenazine operon carriers.
