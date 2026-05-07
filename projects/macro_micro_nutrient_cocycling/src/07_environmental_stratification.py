"""
Enrichment 2: Environmental stratification of co-occurrence.

Joins ncbi_env metadata (isolation_source, env_broad_scale, host) to
pangenome species via genome accessions. Groups into broad environment
categories and computes P×Metal, N×Metal, Phz×Metal co-occurrence per
environment.

Outputs:
  data/env_species_mapping.csv — species-to-environment assignments
  data/env_cooccurrence.csv — co-occurrence statistics by environment
"""

import sys
import os

sys.path.insert(0, os.path.dirname(__file__))
from berdl_notebook_utils.setup_spark_session import get_spark_session

spark = get_spark_session()

DATA_DIR = os.path.join(os.path.dirname(__file__), '..', 'data')

print("Step 1: Join ncbi_env to genomes to species clades...")

env_query = """
SELECT
    g.gtdb_species_clade_id,
    env.harmonized_name,
    LOWER(env.content) AS content
FROM kbase_ke_pangenome.ncbi_env env
JOIN kbase_ke_pangenome.genome g
    ON env.accession = g.ncbi_biosample_id
WHERE env.harmonized_name IN ('isolation_source', 'env_broad_scale', 'host')
    AND env.content IS NOT NULL
    AND LOWER(env.content) NOT IN ('missing', 'not collected', 'not applicable',
        'unknown', 'na', 'n/a', 'none', 'other', '', 'not available')
"""

print("Running join query...")
env_df = spark.sql(env_query).toPandas()
print(f"Got {len(env_df)} environment annotations across {env_df['gtdb_species_clade_id'].nunique()} species")

print("\nTop harmonized_name counts:")
print(env_df['harmonized_name'].value_counts().to_string())

soil_terms = ['soil', 'rhizosphere', 'rhizoplane', 'root', 'mycorrhiz', 'compost',
              'peat', 'sediment', 'mud', 'dirt', 'loam', 'humus', 'manure']
plant_terms = ['plant', 'leaf', 'stem', 'flower', 'fruit', 'seed', 'phyllosphere',
               'endophyte', 'epiphyte', 'xylem', 'phloem', 'nodule', 'legume',
               'rice', 'wheat', 'maize', 'corn', 'soybean', 'arabidopsis',
               'crop', 'grass', 'tree', 'wood', 'bark', 'moss']
marine_terms = ['marine', 'ocean', 'sea', 'seawater', 'coral', 'sponge',
                'coastal', 'estuar', 'tidal', 'mangrove', 'salt marsh',
                'deep sea', 'hydrothermal']
freshwater_terms = ['freshwater', 'lake', 'river', 'stream', 'pond', 'wetland',
                    'bog', 'spring', 'groundwater', 'aquifer', 'drinking water',
                    'wastewater', 'activated sludge', 'biofilm']
human_terms = ['human', 'gut', 'feces', 'stool', 'intestin', 'colon', 'oral',
               'skin', 'blood', 'urine', 'sputum', 'nasal', 'vaginal',
               'clinical', 'patient', 'hospital']
animal_terms = ['animal', 'insect', 'cattle', 'bovine', 'swine', 'pig', 'poultry',
                'chicken', 'fish', 'shrimp', 'bee', 'ant', 'termite',
                'mouse', 'rat', 'dog', 'cat', 'bird', 'rumen', 'fecal']

def classify_environment(content):
    c = str(content).lower()
    for term in soil_terms:
        if term in c:
            return 'soil/rhizosphere'
    for term in plant_terms:
        if term in c:
            return 'plant-associated'
    for term in marine_terms:
        if term in c:
            return 'marine'
    for term in freshwater_terms:
        if term in c:
            return 'freshwater/engineered'
    for term in human_terms:
        if term in c:
            return 'human-associated'
    for term in animal_terms:
        if term in c:
            return 'animal-associated'
    return 'other'

env_df['environment'] = env_df['content'].apply(classify_environment)

print("\nEnvironment classification:")
print(env_df['environment'].value_counts().to_string())

species_env = (env_df.groupby('gtdb_species_clade_id')['environment']
               .agg(lambda x: x.value_counts().index[0])
               .reset_index()
               .rename(columns={'environment': 'primary_env'}))

print(f"\nSpecies with environment assignments: {len(species_env)}")
print("Species per environment:")
print(species_env['primary_env'].value_counts().to_string())

species_env.to_csv(os.path.join(DATA_DIR, 'env_species_mapping.csv'), index=False)

import pandas as pd
import numpy as np
from scipy import stats

gene_df = pd.read_csv(os.path.join(DATA_DIR, 'species_gene_families.csv'))
merged = gene_df.merge(species_env, on='gtdb_species_clade_id', how='inner')
print(f"\nSpecies with both gene families and environment: {len(merged)}")

pairs = {
    'P x Metal': ('has_P_acquisition', 'has_metal_handling'),
    'N x Metal': ('has_N_fixation', 'has_metal_handling'),
    'Phz x Metal': ('has_phz_operon', 'has_metal_handling'),
}

out_rows = []
for env, grp in merged.groupby('primary_env'):
    n = len(grp)
    if n < 10:
        continue

    for pair_label, (col_a, col_b) in pairs.items():
        if col_a not in grp.columns or col_b not in grp.columns:
            continue

        a = grp[col_a].values.astype(int)
        b = grp[col_b].values.astype(int)

        n11 = int(np.sum((a == 1) & (b == 1)))
        n10 = int(np.sum((a == 1) & (b == 0)))
        n01 = int(np.sum((a == 0) & (b == 1)))
        n00 = int(np.sum((a == 0) & (b == 0)))

        n_a = n11 + n10
        n_b = n11 + n01
        expected = (n_a / n) * (n_b / n) * n if n > 0 else 0
        enrichment = n11 / expected if expected > 0 else 0

        if n11 == 0 or n10 == 0 or n01 == 0 or n00 == 0:
            n11_h, n10_h, n01_h, n00_h = n11+0.5, n10+0.5, n01+0.5, n00+0.5
        else:
            n11_h, n10_h, n01_h, n00_h = n11, n10, n01, n00

        log_or = np.log(n11_h * n00_h / (n10_h * n01_h))
        se = np.sqrt(1/n11_h + 1/n10_h + 1/n01_h + 1/n00_h)

        ct = np.array([[n11, n10], [n01, n00]])
        try:
            or_val, fisher_p = stats.fisher_exact(ct)
        except ValueError:
            or_val, fisher_p = 0, 1

        out_rows.append({
            'environment': env,
            'pair': pair_label,
            'n_species': n,
            'n_A': n_a,
            'n_B': n_b,
            'n_both': n11,
            'expected': round(expected, 1),
            'enrichment': round(enrichment, 3),
            'log_OR': round(log_or, 4),
            'se': round(se, 4),
            'fisher_p': fisher_p,
        })

env_cooc = pd.DataFrame(out_rows)
env_cooc.to_csv(os.path.join(DATA_DIR, 'env_cooccurrence.csv'), index=False)
print(f"\nSaved env_cooccurrence.csv ({len(env_cooc)} rows)")

print("\n=== Co-occurrence by Environment ===\n")
for pair_label in pairs:
    print(f"\n{pair_label}:")
    sub = env_cooc[env_cooc['pair'] == pair_label].sort_values('enrichment', ascending=False)
    for _, r in sub.iterrows():
        sig = '*' if r['fisher_p'] < 0.05 else ''
        print(f"  {r['environment']:>25s}: n={r['n_species']:>5d}  "
              f"both={r['n_both']:>4d}  enrich={r['enrichment']:.2f}x  "
              f"log-OR={r['log_OR']:+.3f}  p={r['fisher_p']:.2e} {sig}")
