"""
Line 3: KEGG pathway co-membership test.

Tests whether P-acquisition and metal-handling gene clusters share KEGG
pathway membership more often than expected by chance, using pathway
annotations from eggnog_mapper_annotations.
"""
import os
import csv
import numpy as np
import pandas as pd
from collections import defaultdict
from itertools import product

from berdl_notebook_utils.setup_spark_session import get_spark_session


def parse_pathways(pathway_str):
    if not pathway_str or pathway_str == "-" or pathway_str == "nan":
        return set()
    return {p.strip() for p in str(pathway_str).split(",") if p.strip() and p.strip() != "-"}


def main():
    data_dir = "projects/macro_micro_nutrient_cocycling/data"
    output_dir = os.path.join(data_dir, "kegg_pathways")
    os.makedirs(output_dir, exist_ok=True)

    spark = get_spark_session()

    # ── Step 1: Get pathway annotations for target gene clusters ──
    print("Step 1: Query KEGG pathway annotations for target gene families")

    # P-acquisition gene clusters
    p_query = """
    SELECT gc.gene_cluster_id, gc.gtdb_species_clade_id, ba.gene, ba.kegg_orthology_id,
           e.KEGG_Pathway, e.KEGG_ko
    FROM kbase_ke_pangenome.bakta_annotations ba
    JOIN kbase_ke_pangenome.gene_cluster gc ON ba.gene_cluster_id = gc.gene_cluster_id
    LEFT JOIN kbase_ke_pangenome.eggnog_mapper_annotations e ON gc.gene_cluster_id = e.query_name
    WHERE ba.gene IN ('phoA', 'pstA', 'pstB', 'pstC', 'pstS', 'phnC', 'phnD', 'phnE')
       OR ba.kegg_orthology_id = 'K01077'
    """

    p_pfam_query = """
    SELECT gc.gene_cluster_id, gc.gtdb_species_clade_id, 'phoD' as gene, NULL as kegg_orthology_id,
           e.KEGG_Pathway, e.KEGG_ko
    FROM kbase_ke_pangenome.bakta_pfam_domains bp
    JOIN kbase_ke_pangenome.gene_cluster gc ON bp.gene_cluster_id = gc.gene_cluster_id
    LEFT JOIN kbase_ke_pangenome.eggnog_mapper_annotations e ON gc.gene_cluster_id = e.query_name
    WHERE bp.pfam_id LIKE 'PF09423%'
    """

    # Metal-handling gene clusters
    m_query = """
    SELECT gc.gene_cluster_id, gc.gtdb_species_clade_id, ba.gene, ba.kegg_orthology_id,
           e.KEGG_Pathway, e.KEGG_ko
    FROM kbase_ke_pangenome.bakta_annotations ba
    JOIN kbase_ke_pangenome.gene_cluster gc ON ba.gene_cluster_id = gc.gene_cluster_id
    LEFT JOIN kbase_ke_pangenome.eggnog_mapper_annotations e ON gc.gene_cluster_id = e.query_name
    WHERE ba.gene IN ('copA', 'corA')
       OR ba.kegg_orthology_id = 'K17686'
    """

    m_pfam_query = """
    SELECT gc.gene_cluster_id, gc.gtdb_species_clade_id,
           CASE WHEN bp.pfam_id LIKE 'PF02421%' THEN 'feoB'
                WHEN bp.pfam_id LIKE 'PF00403%' THEN 'HMA'
           END as gene,
           NULL as kegg_orthology_id,
           e.KEGG_Pathway, e.KEGG_ko
    FROM kbase_ke_pangenome.bakta_pfam_domains bp
    JOIN kbase_ke_pangenome.gene_cluster gc ON bp.gene_cluster_id = gc.gene_cluster_id
    LEFT JOIN kbase_ke_pangenome.eggnog_mapper_annotations e ON gc.gene_cluster_id = e.query_name
    WHERE bp.pfam_id LIKE 'PF02421%' OR bp.pfam_id LIKE 'PF00403%'
    """

    # N-fixation gene clusters
    n_query = """
    SELECT gc.gene_cluster_id, gc.gtdb_species_clade_id, ba.gene, ba.kegg_orthology_id,
           e.KEGG_Pathway, e.KEGG_ko
    FROM kbase_ke_pangenome.bakta_annotations ba
    JOIN kbase_ke_pangenome.gene_cluster gc ON ba.gene_cluster_id = gc.gene_cluster_id
    LEFT JOIN kbase_ke_pangenome.eggnog_mapper_annotations e ON gc.gene_cluster_id = e.query_name
    WHERE ba.kegg_orthology_id IN ('K02588', 'K02586')
    """

    print("  Querying P-acquisition genes...")
    p_df = spark.sql(p_query).toPandas()
    p_pfam_df = spark.sql(p_pfam_query).toPandas()
    p_all = pd.concat([p_df, p_pfam_df], ignore_index=True)
    print(f"  P-acquisition: {len(p_all)} gene cluster rows (all species)")

    print("  Querying metal-handling genes...")
    m_df = spark.sql(m_query).toPandas()
    m_pfam_df = spark.sql(m_pfam_query).toPandas()
    m_all = pd.concat([m_df, m_pfam_df], ignore_index=True)
    print(f"  Metal-handling: {len(m_all)} gene cluster rows (all species)")

    print("  Querying N-fixation genes...")
    n_df = spark.sql(n_query).toPandas()
    print(f"  N-fixation: {len(n_df)} gene cluster rows (all species)")

    _env = pd.read_csv(os.path.join(data_dir, "env_species_mapping.csv"))
    _sp = set(_env[_env['primary_env'].isin(['soil/rhizosphere', 'plant-associated'])]['gtdb_species_clade_id'])
    p_all = p_all[p_all['gtdb_species_clade_id'].isin(_sp)].copy()
    m_all = m_all[m_all['gtdb_species_clade_id'].isin(_sp)].copy()
    n_df = n_df[n_df['gtdb_species_clade_id'].isin(_sp)].copy()
    del _env, _sp
    print(f"  After soil+plant filter: P={len(p_all)}, M={len(m_all)}, N={len(n_df)}")

    # ── Step 2: Parse pathways per gene family ──
    print("\nStep 2: Parse pathway memberships per gene family")

    def get_family_pathways(df, group_name):
        family_pathways = defaultdict(set)
        for _, row in df.iterrows():
            gene = row["gene"] if pd.notna(row["gene"]) else "unknown"
            pw = parse_pathways(row.get("KEGG_Pathway", ""))
            if pw:
                family_pathways[gene] |= pw
        return family_pathways

    p_pathways = get_family_pathways(p_all, "P")
    m_pathways = get_family_pathways(m_all, "M")
    n_pathways = get_family_pathways(n_df, "N")

    print("\n  P-acquisition pathway counts:")
    for gene, pw in sorted(p_pathways.items()):
        print(f"    {gene}: {len(pw)} pathways")

    print("\n  Metal-handling pathway counts:")
    for gene, pw in sorted(m_pathways.items()):
        print(f"    {gene}: {len(pw)} pathways")

    print("\n  N-fixation pathway counts:")
    for gene, pw in sorted(n_pathways.items()):
        print(f"    {gene}: {len(pw)} pathways")

    # ── Step 3: Gene-family-level pathway overlap ──
    print("\nStep 3: Gene-family-level pathway overlap (P × M)")

    all_p_pathways = set()
    for pw in p_pathways.values():
        all_p_pathways |= pw

    all_m_pathways = set()
    for pw in m_pathways.values():
        all_m_pathways |= pw

    all_n_pathways = set()
    for pw in n_pathways.values():
        all_n_pathways |= pw

    pm_shared = all_p_pathways & all_m_pathways
    pn_shared = all_p_pathways & all_n_pathways
    nm_shared = all_n_pathways & all_m_pathways

    print(f"  P pathways: {len(all_p_pathways)}")
    print(f"  M pathways: {len(all_m_pathways)}")
    print(f"  N pathways: {len(all_n_pathways)}")
    print(f"  P ∩ M shared: {len(pm_shared)}")
    print(f"  P ∩ N shared: {len(pn_shared)}")
    print(f"  N ∩ M shared: {len(nm_shared)}")

    if pm_shared:
        print(f"\n  Shared P × M pathways:")
        for pw in sorted(pm_shared):
            print(f"    {pw}")

    # ── Step 4: Pairwise gene-level overlap ──
    print("\nStep 4: Pairwise gene-level pathway overlap")

    pair_results = []
    for p_gene, p_pw in sorted(p_pathways.items()):
        for m_gene, m_pw in sorted(m_pathways.items()):
            shared = p_pw & m_pw
            union = p_pw | m_pw
            jaccard = len(shared) / len(union) if union else 0
            pair_results.append({
                "gene_x": p_gene,
                "gene_y": m_gene,
                "group_x": "P",
                "group_y": "M",
                "n_pathways_x": len(p_pw),
                "n_pathways_y": len(m_pw),
                "n_shared": len(shared),
                "n_union": len(union),
                "jaccard": round(jaccard, 4),
                "shared_pathways": ";".join(sorted(shared)) if shared else "",
            })
            if shared:
                print(f"  {p_gene} × {m_gene}: {len(shared)} shared (J={jaccard:.3f}): {', '.join(sorted(shared))}")

    for n_gene, n_pw in sorted(n_pathways.items()):
        for m_gene, m_pw in sorted(m_pathways.items()):
            shared = n_pw & m_pw
            union = n_pw | m_pw
            jaccard = len(shared) / len(union) if union else 0
            pair_results.append({
                "gene_x": n_gene,
                "gene_y": m_gene,
                "group_x": "N",
                "group_y": "M",
                "n_pathways_x": len(n_pw),
                "n_pathways_y": len(m_pw),
                "n_shared": len(shared),
                "n_union": len(union),
                "jaccard": round(jaccard, 4),
                "shared_pathways": ";".join(sorted(shared)) if shared else "",
            })

    # ── Step 5: Species-level pathway co-membership ──
    print("\nStep 5: Species-level pathway co-membership")

    species_p = defaultdict(lambda: defaultdict(set))
    for _, row in p_all.iterrows():
        sid = row["gtdb_species_clade_id"]
        pw = parse_pathways(row.get("KEGG_Pathway", ""))
        if pw:
            species_p[sid]["all"] |= pw

    species_m = defaultdict(lambda: defaultdict(set))
    for _, row in m_all.iterrows():
        sid = row["gtdb_species_clade_id"]
        pw = parse_pathways(row.get("KEGG_Pathway", ""))
        if pw:
            species_m[sid]["all"] |= pw

    both_species = set(species_p.keys()) & set(species_m.keys())
    print(f"  Species with P pathway annotations: {len(species_p)}")
    print(f"  Species with M pathway annotations: {len(species_m)}")
    print(f"  Species with both: {len(both_species)}")

    n_shared_any = 0
    shared_counts = []
    jaccard_vals = []

    for sid in both_species:
        p_pw = species_p[sid]["all"]
        m_pw = species_m[sid]["all"]
        shared = p_pw & m_pw
        union = p_pw | m_pw
        if shared:
            n_shared_any += 1
        shared_counts.append(len(shared))
        jaccard_vals.append(len(shared) / len(union) if union else 0)

    shared_counts = np.array(shared_counts)
    jaccard_vals = np.array(jaccard_vals)

    print(f"  Species with ≥1 shared pathway: {n_shared_any}/{len(both_species)} ({n_shared_any/len(both_species)*100:.1f}%)")
    print(f"  Mean shared pathways per species: {np.mean(shared_counts):.2f}")
    print(f"  Median shared pathways: {np.median(shared_counts):.0f}")
    print(f"  Mean Jaccard: {np.mean(jaccard_vals):.4f}")

    # ── Step 6: Permutation test ──
    print("\nStep 6: Permutation test (1000 shuffles)")

    rng = np.random.default_rng(42)
    all_pathways_pool = set()
    for pw in species_p.values():
        all_pathways_pool |= pw["all"]
    for pw in species_m.values():
        all_pathways_pool |= pw["all"]
    pathways_list = sorted(all_pathways_pool)

    obs_mean_shared = np.mean(shared_counts)
    null_means = np.empty(1000)

    species_list = list(both_species)
    p_sizes = [len(species_p[s]["all"]) for s in species_list]
    m_sizes = [len(species_m[s]["all"]) for s in species_list]

    for perm_i in range(1000):
        if perm_i % 200 == 0:
            print(f"  Permutation {perm_i}/1000...")
        perm_shared = []
        for i in range(len(species_list)):
            perm_p = set(rng.choice(pathways_list, min(p_sizes[i], len(pathways_list)), replace=False))
            perm_m = set(rng.choice(pathways_list, min(m_sizes[i], len(pathways_list)), replace=False))
            perm_shared.append(len(perm_p & perm_m))
        null_means[perm_i] = np.mean(perm_shared)

    z_score = (obs_mean_shared - np.mean(null_means)) / np.std(null_means) if np.std(null_means) > 0 else 0
    p_value = np.mean(null_means >= obs_mean_shared)

    print(f"\n  Observed mean shared pathways: {obs_mean_shared:.2f}")
    print(f"  Null mean (mean ± std): {np.mean(null_means):.2f} ± {np.std(null_means):.2f}")
    print(f"  Z-score: {z_score:.2f}")
    print(f"  Permutation p: {p_value:.4f}")

    # ── Step 7: Save results ──
    print("\nStep 7: Save results")

    pd.DataFrame(pair_results).to_csv(
        os.path.join(output_dir, "pairwise_pathway_overlap.csv"), index=False
    )

    summary = {
        "n_species_both": len(both_species),
        "n_shared_any": n_shared_any,
        "pct_shared": round(n_shared_any / len(both_species) * 100, 1) if both_species else 0,
        "mean_shared_pathways": round(obs_mean_shared, 2),
        "median_shared_pathways": float(np.median(shared_counts)),
        "mean_jaccard": round(float(np.mean(jaccard_vals)), 4),
        "null_mean": round(float(np.mean(null_means)), 2),
        "null_std": round(float(np.std(null_means)), 2),
        "z_score": round(z_score, 2),
        "perm_p": round(p_value, 4),
        "n_p_pathways": len(all_p_pathways),
        "n_m_pathways": len(all_m_pathways),
        "n_pm_shared": len(pm_shared),
    }

    with open(os.path.join(output_dir, "pathway_comembership_summary.csv"), "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["metric", "value"])
        for k, v in summary.items():
            writer.writerow([k, v])

    print(f"\nDone. Results in {output_dir}/")
    spark.stop()


if __name__ == "__main__":
    main()
