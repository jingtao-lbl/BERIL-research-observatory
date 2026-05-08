"""
Enrichment 1: Positive-control species check.

Look up well-characterized rhizosphere/plant-associated model organisms
in the pangenome data and report which gene family groups each encodes.
Selects one representative clade per target organism (the base species
name without GTDB suffixes, when available).

Outputs:
  data/positive_controls.csv
"""

import sys
import os

sys.path.insert(0, os.path.dirname(__file__))
from berdl_notebook_utils.setup_spark_session import get_spark_session

spark = get_spark_session()

DATA_DIR = os.path.join(os.path.dirname(__file__), '..', 'data')

targets = {
    'Pseudomonas fluorescens': 'Pseudomonas_E_fluorescens',
    'Pseudomonas protegens': 'Pseudomonas_E_protegens',
    'Bradyrhizobium diazoefficiens': 'Bradyrhizobium_diazoefficiens',
    'Sinorhizobium meliloti': 'Sinorhizobium_meliloti',
    'Mesorhizobium loti': 'Mesorhizobium_loti',
    'Streptomyces coelicolor': 'Streptomyces_anthocyanicus',
    'Methylobacterium extorquens': 'Methylobacterium_extorquens',
    'Pseudomonas chlororaphis': 'Pseudomonas_E_chlororaphis',
}

gtdb_patterns = list(targets.values())
like_clauses = " OR ".join([f"gsc.GTDB_species LIKE 's__{p}%'" for p in gtdb_patterns])

query = f"""
SELECT gsc.gtdb_species_clade_id, gsc.GTDB_species, gsc.GTDB_taxonomy,
       gsc.no_clustered_genomes_filtered AS n_genomes
FROM kbase_ke_pangenome.gtdb_species_clade gsc
WHERE {like_clauses}
ORDER BY gsc.GTDB_species
"""

print("Searching for target species in GTDB...")
results = spark.sql(query).toPandas()
print(f"Found {len(results)} matching species clades")

import pandas as pd
gene_df = pd.read_csv(os.path.join(DATA_DIR, 'species_gene_families.csv'))

matched = gene_df.merge(results[['gtdb_species_clade_id', 'GTDB_species', 'n_genomes']],
                         on='gtdb_species_clade_id', how='inner')

print(f"Matched {len(matched)} species in gene family data")

representative = {}
for common_name, gtdb_prefix in targets.items():
    exact_name = f's__{gtdb_prefix}'
    subset = matched[matched['GTDB_species'] == exact_name]
    if len(subset) == 1:
        representative[common_name] = subset.iloc[0]
    elif len(subset) > 1:
        representative[common_name] = subset.sort_values('n_genomes', ascending=False).iloc[0]
    else:
        all_matches = matched[matched['GTDB_species'].str.startswith(exact_name)]
        if len(all_matches) > 0:
            representative[common_name] = all_matches.sort_values('n_genomes', ascending=False).iloc[0]
        else:
            print(f"  WARNING: No match for {common_name} ({gtdb_prefix})")

print(f"\nSelected {len(representative)} representative type strains:")

p_genes = ['phoA', 'phoD_pfam', 'pstA', 'pstB', 'pstC', 'pstS', 'phnC', 'phnD', 'phnE']
n_genes = ['nifH', 'nifD']
m_genes = ['copA', 'corA', 'feoB_pfam', 'HMA_pfam']
phz_genes = ['phzF', 'phzA', 'phzB', 'phzD', 'phzG', 'phzS', 'phzM']

out_rows = []
for common_name, row in representative.items():
    species = row['GTDB_species']
    clade_id = row['gtdb_species_clade_id']

    p_list = [g for g in p_genes if row.get(f'has_{g}', 0) == 1]
    n_list = [g for g in n_genes if row.get(f'has_{g}', 0) == 1]
    m_list = [g for g in m_genes if row.get(f'has_{g}', 0) == 1]
    phz_list = [g for g in phz_genes if row.get(f'has_{g}', 0) == 1]

    out_rows.append({
        'common_name': common_name,
        'GTDB_species': species,
        'clade_id': clade_id,
        'n_genomes': int(row.get('n_genomes', 0)),
        'P_acquisition': ', '.join(p_list) if p_list else 'none',
        'has_P': 1 if p_list else 0,
        'N_fixation': ', '.join(n_list) if n_list else 'none',
        'has_N': 1 if n_list else 0,
        'Metal_handling': ', '.join(m_list) if m_list else 'none',
        'has_M': 1 if m_list else 0,
        'Phenazine': ', '.join(phz_list) if phz_list else 'none',
        'has_Phz': 1 if phz_list else 0,
        'n_phz_genes': int(row.get('n_phz_genes', 0)),
        'has_phz_operon': int(row.get('has_phz_operon', 0)),
    })

    print(f"\n{common_name} -> {species} ({clade_id}, n={int(row.get('n_genomes', 0))}):")
    print(f"  P-acquisition: {', '.join(p_list) if p_list else 'NONE'}")
    print(f"  N-fixation:    {', '.join(n_list) if n_list else 'NONE'}")
    print(f"  Metal-handling: {', '.join(m_list) if m_list else 'NONE'}")
    print(f"  Phenazine:     {', '.join(phz_list) if phz_list else 'NONE'} (operon={int(row.get('has_phz_operon', 0))})")

out_df = pd.DataFrame(out_rows)
out_df.to_csv(os.path.join(DATA_DIR, 'positive_controls.csv'), index=False)
print(f"\nSaved positive_controls.csv ({len(out_df)} representative species)")
print("\nNote: S. coelicolor A3(2) is reclassified as s__Streptomyces_anthocyanicus in GTDB R214 (GCF_008931305.1).")
