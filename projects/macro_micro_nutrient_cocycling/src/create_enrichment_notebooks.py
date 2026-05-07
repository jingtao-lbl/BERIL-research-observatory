"""Create enrichment notebooks (NB05-NB07) and execute them."""
import json
import os
import subprocess

PROJ = os.path.join(os.path.dirname(__file__), '..')
NB_DIR = os.path.join(PROJ, 'notebooks')
os.makedirs(NB_DIR, exist_ok=True)

def make_notebook(cells):
    return {
        "nbformat": 4, "nbformat_minor": 5,
        "metadata": {
            "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
            "language_info": {"name": "python", "version": "3.10.0"}
        },
        "cells": cells
    }

def md_cell(source, cid):
    return {"cell_type": "markdown", "metadata": {}, "source": source.split('\n'), "id": cid}

def code_cell(source, cid):
    return {"cell_type": "code", "metadata": {}, "source": source.split('\n'),
            "outputs": [], "execution_count": None, "id": cid}


# --- NB05: Positive Controls ---
nb05 = make_notebook([
    md_cell("# NB05: Positive-Control Species Check\n\n"
            "Validates the co-occurrence signal against 7 well-characterized model organisms.\n"
            "Streptomyces coelicolor A3(2) is absent from GTDB R214.\n\n"
            "Run equivalent script: `python src/06_positive_controls.py`", "c0"),
    code_cell("import pandas as pd\nimport os\nDATA_DIR = os.path.join('..', 'data')", "c1"),
    code_cell("pc = pd.read_csv(os.path.join(DATA_DIR, 'positive_controls.csv'))\n"
              "print(f'{len(pc)} representative species')\n"
              "pc[['common_name','GTDB_species','P_acquisition','N_fixation','Metal_handling','Phenazine','has_phz_operon']]", "c2"),
    md_cell("## Summary\n\n"
            "- All 7 targets encode both P-acquisition and metal-handling genes\n"
            "- B. diazoefficiens and S. meliloti carry nifH+nifD (known diazotrophs)\n"
            "- S. meliloti has the broadest repertoire: all 4 metal families + 8/9 P families\n"
            "- P. chlororaphis is the only positive control with a phenazine operon", "c3"),
])

# --- NB06: Environmental Stratification ---
nb06 = make_notebook([
    md_cell("# NB06: Environmental Stratification\n\n"
            "Tests whether P×Metal and N×Metal co-occurrence varies by reported isolation environment.\n"
            "Uses `kbase_ke_pangenome.ncbi_env` metadata.\n\n"
            "Run equivalent script: `python src/07_environmental_stratification.py`", "c0"),
    code_cell("import pandas as pd\nimport os\nDATA_DIR = os.path.join('..', 'data')", "c1"),
    code_cell("env_map = pd.read_csv(os.path.join(DATA_DIR, 'env_species_mapping.csv'))\n"
              "print(f'{len(env_map)} species with environment assignments')\n"
              "print()\nprint('Species per environment:')\n"
              "print(env_map['primary_env'].value_counts().to_string())", "c2"),
    md_cell("## Co-occurrence by Environment", "c3"),
    code_cell("env_cooc = pd.read_csv(os.path.join(DATA_DIR, 'env_cooccurrence.csv'))\n"
              "for pair in ['P x Metal', 'N x Metal', 'Phz x Metal']:\n"
              "    print(f'\\n{pair}:')\n"
              "    sub = env_cooc[env_cooc['pair']==pair].sort_values('enrichment', ascending=False)\n"
              "    for _, r in sub.iterrows():\n"
              "        sig = '*' if r['fisher_p'] < 0.05 else ''\n"
              "        print(f\"  {r['environment']:>25s}: n={r['n_species']:>5d}  \"\n"
              "              f\"both={r['n_both']:>4d}  enrich={r['enrichment']:.2f}x  \"\n"
              "              f\"log-OR={r['log_OR']:+.3f}  p={r['fisher_p']:.2e} {sig}\")", "c4"),
    md_cell("## Key finding\n\n"
            "P×Metal is strongest in plant-associated and soil/rhizosphere environments,\n"
            "consistent with the Fe-oxyhydroxide mineral surface hypothesis.", "c5"),
])

# --- NB07: Forest Plot ---
nb07 = make_notebook([
    md_cell("# NB07: Per-Phylum Forest Plot\n\n"
            "Log-odds ratios with 95% CIs for P×Metal, N×Metal, Phz×Metal across GTDB phyla.\n\n"
            "Run equivalent script: `python src/08_forest_plot.py`", "c0"),
    code_cell("import pandas as pd\nimport os\nDATA_DIR = os.path.join('..', 'data')\nFIG_DIR = os.path.join('..', 'figures')", "c1"),
    code_cell("forest = pd.read_csv(os.path.join(DATA_DIR, 'forest_plot_data.csv'))\n"
              "print(f'{len(forest)} phylum-pair combinations')\n"
              "print(f'Phyla represented: {forest[\"phylum\"].nunique()}')\n"
              "print(f'Pairs: {forest[\"pair\"].unique().tolist()}')", "c2"),
    md_cell("## Top effects per pair", "c3"),
    code_cell("for pair in ['P x Metal', 'N x Metal', 'Phz x Metal']:\n"
              "    print(f'\\n{pair} — top 5 phyla by log-OR:')\n"
              "    sub = forest[forest['pair']==pair].sort_values('log_OR', ascending=False).head(5)\n"
              "    for _, r in sub.iterrows():\n"
              "        sig = '*' if r['fisher_p'] < 0.05 else ''\n"
              "        print(f\"  {r['phylum']:>30s}: log-OR={r['log_OR']:+.3f} \"\n"
              "              f\"[{r['ci_lo']:+.3f}, {r['ci_hi']:+.3f}] p={r['fisher_p']:.2e} {sig}\")", "c4"),
    md_cell("## Forest Plot (Figure 2)", "c5"),
    code_cell("from IPython.display import Image\n"
              "Image(os.path.join(FIG_DIR, 'forest_plot.png'), width=900)", "c6"),
])


notebooks = {
    'NB05_positive_controls.ipynb': nb05,
    'NB06_environmental_stratification.ipynb': nb06,
    'NB07_forest_plot.ipynb': nb07,
}

for name, nb in notebooks.items():
    for cell in nb['cells']:
        if isinstance(cell.get('source'), list):
            cell['source'] = [line + '\n' for line in cell['source'][:-1]] + [cell['source'][-1]]

    path = os.path.join(NB_DIR, name)
    with open(path, 'w') as f:
        json.dump(nb, f, indent=1)
    print(f"Created {path}")

print("\nExecuting notebooks...")
for name in notebooks:
    path = os.path.join(NB_DIR, name)
    print(f"Executing {name}...")
    result = subprocess.run(
        ['jupyter', 'nbconvert', '--to', 'notebook', '--execute', '--inplace',
         '--ExecutePreprocessor.timeout=300', path],
        capture_output=True, text=True, cwd=NB_DIR
    )
    if result.returncode == 0:
        print(f"  OK: {name}")
    else:
        print(f"  FAILED: {name}")
        print(f"  stderr: {result.stderr[:500]}")
