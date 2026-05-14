"""
Enrichment 3: Per-phylum forest plot of co-occurrence effect sizes.

Computes log-odds ratios with 95% CIs for P×Metal, N×Metal, and Phz×Metal
co-occurrence within each major phylum (n>=50 species). Produces a forest
plot figure with three panels.

Outputs:
  data/forest_plot_data.csv — per-phylum log-OR and CIs
  figures/forest_plot.png — three-panel forest plot
"""

import os
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

DATA_DIR = os.path.join(os.path.dirname(__file__), '..', 'data')
FIG_DIR = os.path.join(os.path.dirname(__file__), '..', 'figures')
os.makedirs(FIG_DIR, exist_ok=True)

print("Loading data...")
df = pd.read_csv(os.path.join(DATA_DIR, 'species_gene_families.csv'))
tax = pd.read_csv(os.path.join(DATA_DIR, 'species_taxonomy.csv'))
_env = pd.read_csv(os.path.join(DATA_DIR, 'env_species_mapping.csv'))
_sp = set(_env[_env['primary_env'].isin(['soil/rhizosphere', 'plant-associated'])]['gtdb_species_clade_id'])
df = df[df['gtdb_species_clade_id'].isin(_sp)].copy()
del _env, _sp
merged = df.merge(tax[['gtdb_species_clade_id', 'phylum']], on='gtdb_species_clade_id', how='left')
print(f"v4 soil+plant filter: {len(merged)} species")

pairs = {
    'P x Metal': ('has_P_acquisition', 'has_metal_handling'),
    'N x Metal': ('has_N_fixation', 'has_metal_handling'),
    'Phz x Metal': ('has_phz_operon', 'has_metal_handling'),
}

MIN_PHYLUM_SIZE = 20

rows = []
for pair_label, (col_a, col_b) in pairs.items():
    if col_a not in merged.columns or col_b not in merged.columns:
        print(f"  Skipping {pair_label}: columns not found")
        continue

    for phylum, grp in merged.groupby('phylum'):
        n = len(grp)
        if n < MIN_PHYLUM_SIZE:
            continue

        a = grp[col_a].values.astype(int)
        b = grp[col_b].values.astype(int)

        n11 = int(np.sum((a == 1) & (b == 1)))
        n10 = int(np.sum((a == 1) & (b == 0)))
        n01 = int(np.sum((a == 0) & (b == 1)))
        n00 = int(np.sum((a == 0) & (b == 0)))

        if n11 == 0 or n10 == 0 or n01 == 0 or n00 == 0:
            n11_h = n11 + 0.5
            n10_h = n10 + 0.5
            n01_h = n01 + 0.5
            n00_h = n00 + 0.5
        else:
            n11_h, n10_h, n01_h, n00_h = n11, n10, n01, n00

        log_or = np.log(n11_h * n00_h / (n10_h * n01_h))
        se = np.sqrt(1/n11_h + 1/n10_h + 1/n01_h + 1/n00_h)
        ci_lo = log_or - 1.96 * se
        ci_hi = log_or + 1.96 * se

        _, fisher_p = stats.fisher_exact([[n11, n10], [n01, n00]])

        rows.append({
            'pair': pair_label,
            'phylum': phylum,
            'n_species': n,
            'n11': n11, 'n10': n10, 'n01': n01, 'n00': n00,
            'log_OR': round(log_or, 4),
            'se': round(se, 4),
            'ci_lo': round(ci_lo, 4),
            'ci_hi': round(ci_hi, 4),
            'fisher_p': fisher_p,
        })

forest_df = pd.DataFrame(rows)
forest_df.to_csv(os.path.join(DATA_DIR, 'forest_plot_data.csv'), index=False)
print(f"Saved forest_plot_data.csv ({len(forest_df)} rows)")

for pair_label in pairs:
    sub = forest_df[forest_df['pair'] == pair_label].sort_values('log_OR', ascending=False)
    print(f"\n{pair_label}:")
    for _, r in sub.head(10).iterrows():
        sig = '*' if r['fisher_p'] < 0.05 else ''
        print(f"  {r['phylum']:>30s}: log-OR={r['log_OR']:+.3f} [{r['ci_lo']:+.3f}, {r['ci_hi']:+.3f}]  "
              f"n={r['n_species']:>5d}  p={r['fisher_p']:.2e} {sig}")

fig, axes = plt.subplots(1, 3, figsize=(18, 10), sharey=False)

for ax, pair_label in zip(axes, pairs):
    sub = forest_df[forest_df['pair'] == pair_label].copy()
    sub = sub.sort_values('log_OR', ascending=True)

    y_pos = np.arange(len(sub))
    colors = ['#2166ac' if r['fisher_p'] < 0.05 else '#999999' for _, r in sub.iterrows()]

    ax.barh(y_pos, sub['log_OR'].values, xerr=[sub['log_OR'].values - sub['ci_lo'].values,
            sub['ci_hi'].values - sub['log_OR'].values],
            height=0.6, color=colors, edgecolor='none', alpha=0.8,
            error_kw=dict(ecolor='#333333', capsize=2, linewidth=0.8))

    ax.axvline(x=0, color='black', linewidth=0.8, linestyle='--')
    ax.set_yticks(y_pos)

    ylabels = []
    for _, r in sub.iterrows():
        p = r['phylum'].replace('p__', '')
        ylabels.append(f"{p} (n={r['n_species']})")
    ax.set_yticklabels(ylabels, fontsize=7)

    ax.set_xlabel('log Odds Ratio', fontsize=10)
    ax.set_title(pair_label, fontsize=12, fontweight='bold')

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

plt.tight_layout()
fig.savefig(os.path.join(FIG_DIR, 'forest_plot.png'), dpi=200, bbox_inches='tight')
plt.close(fig)
print(f"\nSaved figures/forest_plot.png")
