"""
Line 2: Operon-distance test using gene-interval proxy.

Uses Spark-side joins to efficiently extract gene positions for
P-acquisition and metal-handling gene clusters, then computes
distances and permutation tests in Python.
"""
import os
import csv
import numpy as np
import pandas as pd
from collections import defaultdict

from berdl_notebook_utils.setup_spark_session import get_spark_session


def parse_gene_id(gene_id):
    last_underscore = gene_id.rfind("_")
    if last_underscore == -1:
        return None, None
    contig = gene_id[:last_underscore]
    try:
        ordinal = int(gene_id[last_underscore + 1:])
    except ValueError:
        return contig, None
    return contig, ordinal


def main():
    data_dir = "projects/macro_micro_nutrient_cocycling/data"
    output_dir = os.path.join(data_dir, "operon_distance")
    os.makedirs(output_dir, exist_ok=True)

    spark = get_spark_session()

    print("Step 1: Create temp views for target gene clusters")

    spark.sql("""
    CREATE OR REPLACE TEMP VIEW bakta_targets AS
    SELECT DISTINCT gene_cluster_id,
        CASE
            WHEN kegg_orthology_id = 'K01077' THEN 'phoA'
            WHEN kegg_orthology_id = 'K17686' THEN 'copA'
            ELSE gene
        END as family,
        CASE
            WHEN gene IN ('pstA','pstB','pstC','pstS','phnC','phnD','phnE')
                 OR kegg_orthology_id = 'K01077' THEN 'P'
            ELSE 'M'
        END as func_group
    FROM kbase_ke_pangenome.bakta_annotations
    WHERE gene IN ('pstA','pstB','pstC','pstS','phnC','phnD','phnE','corA')
       OR kegg_orthology_id IN ('K01077', 'K17686')
    """)

    spark.sql("""
    CREATE OR REPLACE TEMP VIEW pfam_targets AS
    SELECT DISTINCT gene_cluster_id,
        CASE
            WHEN pfam_id LIKE 'PF09423%' THEN 'phoD'
            WHEN pfam_id LIKE 'PF02421%' THEN 'feoB'
            WHEN pfam_id LIKE 'PF00403%' THEN 'HMA'
        END as family,
        CASE
            WHEN pfam_id LIKE 'PF09423%' THEN 'P'
            ELSE 'M'
        END as func_group
    FROM kbase_ke_pangenome.bakta_pfam_domains
    WHERE pfam_id LIKE 'PF09423%'
       OR pfam_id LIKE 'PF02421%'
       OR pfam_id LIKE 'PF00403%'
    """)

    spark.sql("""
    CREATE OR REPLACE TEMP VIEW all_targets AS
    SELECT * FROM bakta_targets
    UNION ALL
    SELECT * FROM pfam_targets
    """)

    cnt = spark.sql("SELECT COUNT(DISTINCT gene_cluster_id) as cnt FROM all_targets").collect()
    print(f"  Target gene clusters: {cnt[0].cnt:,}")

    print("\nStep 2: Load species list (both P and M genes, soil+plant only)")
    import pandas as _pd
    _env = _pd.read_csv(os.path.join(data_dir, "env_species_mapping.csv"))
    soil_plant_ids = set(_env[_env['primary_env'].isin(['soil/rhizosphere', 'plant-associated'])]['gtdb_species_clade_id'])
    del _env

    species_with_both = set()
    species_accessions = {}
    with open(os.path.join(data_dir, "species_gene_families.csv")) as f:
        reader = csv.DictReader(f)
        for row in reader:
            sid = row["gtdb_species_clade_id"]
            if sid not in soil_plant_ids:
                continue
            has_p = int(row["has_P_acquisition"])
            has_m = int(row["has_metal_handling"])
            if has_p and has_m:
                acc = sid.split("--")[-1]
                species_with_both.add(sid)
                species_accessions[acc] = sid

    print(f"  {len(species_with_both)} soil+plant species with both P and metal genes")

    print("\nStep 3: Spark join — target clusters → genes → genomes")
    print("  (This joins ~200K target clusters against 1B gene rows via gene_genecluster_junction)")

    result_df = spark.sql("""
    SELECT gc.gtdb_species_clade_id as species_id,
           g.gene_id,
           g.genome_id,
           t.family,
           t.func_group
    FROM all_targets t
    JOIN kbase_ke_pangenome.gene_genecluster_junction j
        ON t.gene_cluster_id = j.gene_cluster_id
    JOIN kbase_ke_pangenome.gene g
        ON j.gene_id = g.gene_id
    JOIN kbase_ke_pangenome.gene_cluster gc
        ON t.gene_cluster_id = gc.gene_cluster_id
    """)

    print("  Collecting to pandas...")
    pdf = result_df.toPandas()
    print(f"  Got {len(pdf):,} rows (all species)")
    pdf = pdf[pdf["species_id"].isin(soil_plant_ids)].copy()
    print(f"  After soil+plant filter: {len(pdf):,} rows")

    print("\nStep 4: Filter to representative genomes and parse positions")
    pdf["rep_acc"] = pdf["species_id"].str.split("--").str[-1]
    pdf_rep = pdf[pdf["genome_id"] == pdf["rep_acc"]].copy()
    print(f"  Representative-genome rows: {len(pdf_rep):,}")

    pdf_rep["contig"] = pdf_rep["gene_id"].apply(lambda x: x[:x.rfind("_")])
    pdf_rep["ordinal"] = pdf_rep["gene_id"].apply(
        lambda x: int(x[x.rfind("_") + 1:]) if x.rfind("_") >= 0 else -1
    )
    pdf_rep = pdf_rep[pdf_rep["ordinal"] >= 0]

    print(f"  After parsing: {len(pdf_rep):,} rows with valid positions")
    print(f"  Species represented: {pdf_rep['species_id'].nunique()}")

    print("\nStep 5: Compute pairwise distances per species")

    species_results = []
    all_within_dists = []

    for sid, grp in pdf_rep.groupby("species_id"):
        p_genes = grp[grp["func_group"] == "P"][["contig", "ordinal", "family"]].values.tolist()
        m_genes = grp[grp["func_group"] == "M"][["contig", "ordinal", "family"]].values.tolist()

        if not p_genes or not m_genes:
            continue

        total_pairs = 0
        same_contig = 0
        within_dists = []

        for pc, po, pf in p_genes:
            for mc, mo, mf in m_genes:
                total_pairs += 1
                if pc == mc:
                    same_contig += 1
                    d = abs(po - mo)
                    within_dists.append(d)

        pct = same_contig / total_pairs * 100 if total_pairs > 0 else 0

        species_results.append({
            "species_id": sid,
            "n_p_genes": len(p_genes),
            "n_m_genes": len(m_genes),
            "total_pairs": total_pairs,
            "n_same_contig": same_contig,
            "pct_same_contig": round(pct, 1),
            "median_dist": np.median(within_dists) if within_dists else None,
            "min_dist": min(within_dists) if within_dists else None,
            "mean_dist": round(np.mean(within_dists), 1) if within_dists else None,
        })
        all_within_dists.extend(within_dists)

    print(f"  Species with distance data: {len(species_results)}")
    print(f"  Total same-contig pairs: {len(all_within_dists):,}")

    if all_within_dists:
        obs_median = np.median(all_within_dists)
        print(f"  Observed median distance: {obs_median:.1f} genes")
        print(f"  Observed mean distance: {np.mean(all_within_dists):.1f} genes")
        print(f"  Pairs within 5 genes: {sum(1 for d in all_within_dists if d <= 5):,}")
        print(f"  Pairs within 10 genes: {sum(1 for d in all_within_dists if d <= 10):,}")
        print(f"  Pairs within 20 genes: {sum(1 for d in all_within_dists if d <= 20):,}")
    else:
        obs_median = float("nan")

    print("\nStep 6: Permutation test (1000 shuffles)")
    n_perms = 1000
    rng = np.random.default_rng(42)

    species_gene_data = {}
    for sid, grp in pdf_rep.groupby("species_id"):
        p_idx = grp["func_group"] == "P"
        m_idx = grp["func_group"] == "M"
        if p_idx.sum() == 0 or m_idx.sum() == 0:
            continue
        contigs = grp["contig"].values
        ordinals = grp["ordinal"].values
        n_p = p_idx.sum()
        n_m = m_idx.sum()
        species_gene_data[sid] = (contigs, ordinals, n_p, n_m)

    null_medians = np.empty(n_perms)

    for perm_i in range(n_perms):
        if perm_i % 100 == 0:
            print(f"  Permutation {perm_i}/{n_perms}...")

        perm_dists = []
        for sid, (contigs, ordinals, n_p, n_m) in species_gene_data.items():
            n_total = len(ordinals)
            shuf_ord = ordinals.copy()
            rng.shuffle(shuf_ord)

            for i in range(n_p):
                for j in range(n_p, n_p + n_m):
                    if j < n_total and contigs[i] == contigs[j]:
                        perm_dists.append(abs(int(shuf_ord[i]) - int(shuf_ord[j])))

        if perm_dists:
            null_medians[perm_i] = np.median(perm_dists)
        else:
            null_medians[perm_i] = float("nan")

    valid_nulls = null_medians[~np.isnan(null_medians)]
    if len(valid_nulls) > 0 and not np.isnan(obs_median):
        p_value = np.mean(valid_nulls <= obs_median)
        z_score = (obs_median - np.mean(valid_nulls)) / np.std(valid_nulls) if np.std(valid_nulls) > 0 else 0
        print(f"\n  Observed median: {obs_median:.1f}")
        print(f"  Null median (mean ± std): {np.mean(valid_nulls):.1f} ± {np.std(valid_nulls):.1f}")
        print(f"  Z-score: {z_score:.2f}")
        print(f"  P-value (observed <= null): {p_value:.4f}")
    else:
        p_value = float("nan")
        z_score = float("nan")

    print("\nStep 7: Save results")

    results_df = pd.DataFrame(species_results)
    results_df.to_csv(os.path.join(output_dir, "species_distances.csv"), index=False)

    with open(os.path.join(output_dir, "permutation_summary.csv"), "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["metric", "value"])
        writer.writerow(["observed_median", obs_median])
        writer.writerow(["null_mean", np.mean(valid_nulls) if len(valid_nulls) > 0 else "NA"])
        writer.writerow(["null_std", np.std(valid_nulls) if len(valid_nulls) > 0 else "NA"])
        writer.writerow(["z_score", z_score])
        writer.writerow(["p_value", p_value])
        writer.writerow(["n_permutations", n_perms])
        writer.writerow(["n_species", len(species_results)])
        writer.writerow(["n_same_contig_pairs", len(all_within_dists)])

    np.save(os.path.join(output_dir, "null_medians.npy"), null_medians)
    np.save(os.path.join(output_dir, "observed_dists.npy"), np.array(all_within_dists))

    print(f"\nDone. Results in {output_dir}/")
    spark.stop()


if __name__ == "__main__":
    main()
