"""
Line 4: Wang et al. 2021 phytase-siderophore validation.

Defines phytase and siderophore biosynthesis gene families in BERDL,
runs pan-bacterial co-occurrence, stratifies within Burkholderiaceae
(Caballeronia + Paraburkholderia), and tests family-level over-representation.
"""
import os
import csv
import numpy as np
import pandas as pd
from scipy import stats
from collections import defaultdict

from berdl_notebook_utils.setup_spark_session import get_spark_session


def fisher_test(a, b, c, d):
    """2x2 Fisher's exact test. Returns OR, p-value."""
    table = np.array([[a, b], [c, d]])
    or_val, p_val = stats.fisher_exact(table)
    return or_val, p_val


def phi_coefficient(a, b, c, d):
    n = a + b + c + d
    denom = np.sqrt((a + b) * (c + d) * (a + c) * (b + d))
    if denom == 0:
        return 0.0
    return (a * d - b * c) / denom


def jaccard(a, b, c):
    denom = a + b + c
    if denom == 0:
        return 0.0
    return a / denom


def main():
    data_dir = "projects/macro_micro_nutrient_cocycling/data"
    output_dir = os.path.join(data_dir, "wang2021")
    os.makedirs(output_dir, exist_ok=True)

    spark = get_spark_session()

    # ── Step 1: Define phytase gene clusters ──
    print("Step 1: Define phytase gene families")

    phytase_queries = {
        "appA": "SELECT DISTINCT gc.gtdb_species_clade_id FROM kbase_ke_pangenome.bakta_annotations ba JOIN kbase_ke_pangenome.gene_cluster gc ON ba.gene_cluster_id = gc.gene_cluster_id WHERE ba.gene = 'appA'",
        "K01093": "SELECT DISTINCT gc.gtdb_species_clade_id FROM kbase_ke_pangenome.bakta_annotations ba JOIN kbase_ke_pangenome.gene_cluster gc ON ba.gene_cluster_id = gc.gene_cluster_id WHERE ba.kegg_orthology_id = 'K01093'",
        "K14379": "SELECT DISTINCT gc.gtdb_species_clade_id FROM kbase_ke_pangenome.bakta_annotations ba JOIN kbase_ke_pangenome.gene_cluster gc ON ba.gene_cluster_id = gc.gene_cluster_id WHERE ba.kegg_orthology_id = 'K14379'",
        "K09474": "SELECT DISTINCT gc.gtdb_species_clade_id FROM kbase_ke_pangenome.bakta_annotations ba JOIN kbase_ke_pangenome.gene_cluster gc ON ba.gene_cluster_id = gc.gene_cluster_id WHERE ba.kegg_orthology_id = 'K09474'",
        "PF02333": "SELECT DISTINCT gc.gtdb_species_clade_id FROM kbase_ke_pangenome.bakta_pfam_domains bp JOIN kbase_ke_pangenome.gene_cluster gc ON bp.gene_cluster_id = gc.gene_cluster_id WHERE bp.pfam_id LIKE 'PF02333%'",
        "PF13714": "SELECT DISTINCT gc.gtdb_species_clade_id FROM kbase_ke_pangenome.bakta_pfam_domains bp JOIN kbase_ke_pangenome.gene_cluster gc ON bp.gene_cluster_id = gc.gene_cluster_id WHERE bp.pfam_id LIKE 'PF13714%'",
        "phy": "SELECT DISTINCT gc.gtdb_species_clade_id FROM kbase_ke_pangenome.bakta_annotations ba JOIN kbase_ke_pangenome.gene_cluster gc ON ba.gene_cluster_id = gc.gene_cluster_id WHERE ba.gene IN ('phy', 'phyA', 'phyA2', 'phyC')",
    }

    phytase_species = defaultdict(set)
    all_phytase_species = set()

    for name, query in phytase_queries.items():
        rows = spark.sql(query).collect()
        species = {r.gtdb_species_clade_id for r in rows}
        phytase_species[name] = species
        all_phytase_species |= species
        print(f"  {name}: {len(species)} species")

    print(f"  Combined phytase: {len(all_phytase_species)} species")

    # ── Step 2: Define siderophore biosynthesis gene clusters ──
    print("\nStep 2: Define siderophore biosynthesis gene families")

    sid_gene_groups = {
        "enterobactin": ["entA", "entB", "entC", "entD", "entE", "entF"],
        "pyochelin": ["pchA", "pchB", "pchC", "pchD", "pchE", "pchF", "pchG"],
        "pyoverdine": ["pvdA", "pvdD", "pvdE", "pvdF", "pvdL", "pvdS"],
        "hydroxamate": ["iucA", "iucB", "iucC", "iucD"],
    }

    siderophore_species = defaultdict(set)
    all_siderophore_species = set()

    for group_name, genes in sid_gene_groups.items():
        gene_list = ",".join(f"'{g}'" for g in genes)
        query = f"""
        SELECT DISTINCT gc.gtdb_species_clade_id
        FROM kbase_ke_pangenome.bakta_annotations ba
        JOIN kbase_ke_pangenome.gene_cluster gc
            ON ba.gene_cluster_id = gc.gene_cluster_id
        WHERE ba.gene IN ({gene_list})
        """
        rows = spark.sql(query).collect()
        species = {r.gtdb_species_clade_id for r in rows}
        siderophore_species[group_name] = species
        all_siderophore_species |= species
        print(f"  {group_name} ({', '.join(genes)}): {len(species)} species")

    print(f"  Combined siderophore: {len(all_siderophore_species)} species")

    # ── Step 3: Load full species list and taxonomy ──
    print("\nStep 3: Load species list and taxonomy")

    all_species = set()
    with open(os.path.join(data_dir, "species_gene_families.csv")) as f:
        reader = csv.DictReader(f)
        for row in reader:
            all_species.add(row["gtdb_species_clade_id"])
    print(f"  Total species: {len(all_species)}")

    taxonomy = {}
    with open(os.path.join(data_dir, "species_taxonomy.csv")) as f:
        reader = csv.DictReader(f)
        for row in reader:
            sid = row["gtdb_species_clade_id"]
            taxonomy[sid] = {
                "domain": row.get("domain", ""),
                "phylum": row.get("phylum", ""),
                "class": row.get("class", ""),
                "order": row.get("order", ""),
                "family": row.get("family", ""),
                "genus": row.get("genus", ""),
            }

    # ── Step 4: Pan-bacterial co-occurrence test ──
    print("\nStep 4: Pan-bacterial phytase × siderophore co-occurrence")

    n_total = len(all_species)
    has_phyt = all_phytase_species & all_species
    has_sid = all_siderophore_species & all_species
    has_both = has_phyt & has_sid
    has_phyt_only = has_phyt - has_sid
    has_sid_only = has_sid - has_phyt
    has_neither = all_species - has_phyt - has_sid

    a = len(has_both)
    b = len(has_phyt_only)
    c = len(has_sid_only)
    d = len(has_neither)

    or_val, p_val = fisher_test(a, b, c, d)
    phi_val = phi_coefficient(a, b, c, d)
    jacc = jaccard(a, b, c)

    print(f"  Both: {a}")
    print(f"  Phytase only: {b}")
    print(f"  Siderophore only: {c}")
    print(f"  Neither: {d}")
    print(f"  OR: {or_val:.3f}")
    print(f"  Fisher p: {p_val:.2e}")
    print(f"  Phi: {phi_val:.4f}")
    print(f"  Jaccard: {jacc:.4f}")

    # Permutation test
    print("\n  Permutation test (1000 shuffles)...")
    rng = np.random.default_rng(42)
    n_phyt = len(has_phyt)
    n_sid = len(has_sid)
    species_list = list(all_species)
    null_phis = np.empty(1000)

    for i in range(1000):
        perm_phyt = set(rng.choice(species_list, n_phyt, replace=False))
        perm_sid = set(rng.choice(species_list, n_sid, replace=False))
        pa = len(perm_phyt & perm_sid)
        pb = len(perm_phyt - perm_sid)
        pc = len(perm_sid - perm_phyt)
        p_d = n_total - pa - pb - pc
        null_phis[i] = phi_coefficient(pa, pb, pc, p_d)

    perm_z = (phi_val - np.mean(null_phis)) / np.std(null_phis) if np.std(null_phis) > 0 else 0
    perm_p = np.mean(null_phis >= phi_val)
    print(f"  Permutation Z: {perm_z:.2f}")
    print(f"  Permutation p: {perm_p:.4f}")

    # ── Step 5: Per-siderophore-type co-occurrence ──
    print("\nStep 5: Per-siderophore-type co-occurrence with phytase")

    sid_type_results = []
    for stype, sp_set in siderophore_species.items():
        sp_in = sp_set & all_species
        both = has_phyt & sp_in
        phyt_only = has_phyt - sp_in
        sid_only = sp_in - has_phyt
        neither = all_species - has_phyt - sp_in

        sa, sb, sc, sd = len(both), len(phyt_only), len(sid_only), len(neither)
        s_or, s_p = fisher_test(sa, sb, sc, sd)
        s_phi = phi_coefficient(sa, sb, sc, sd)

        sid_type_results.append({
            "siderophore_type": stype,
            "n_species": len(sp_in),
            "n_both": sa,
            "OR": s_or,
            "phi": s_phi,
            "fisher_p": s_p,
        })
        print(f"  {stype}: n={len(sp_in)}, both={sa}, OR={s_or:.2f}, phi={s_phi:.4f}, p={s_p:.2e}")

    # ── Step 6: Burkholderiaceae stratification ──
    print("\nStep 6: Burkholderiaceae stratification")

    burkh_species = {s for s in all_species if taxonomy.get(s, {}).get("family") == "f__Burkholderiaceae"}
    print(f"  Burkholderiaceae species: {len(burkh_species)}")

    # Caballeronia + Paraburkholderia
    cab_para = {s for s in burkh_species
                if taxonomy.get(s, {}).get("genus", "") in ("g__Caballeronia", "g__Paraburkholderia")}
    print(f"  Caballeronia + Paraburkholderia: {len(cab_para)}")

    for label, subset in [("Burkholderiaceae", burkh_species), ("Caballeronia+Paraburkholderia", cab_para)]:
        if len(subset) < 10:
            print(f"  {label}: too few species ({len(subset)}), skipping")
            continue

        b_phyt = has_phyt & subset
        b_sid = has_sid & subset
        b_both = b_phyt & b_sid
        b_phyt_only = b_phyt - b_sid
        b_sid_only = b_sid - b_phyt
        b_neither = subset - b_phyt - b_sid

        sa2, sb2, sc2, sd2 = len(b_both), len(b_phyt_only), len(b_sid_only), len(b_neither)
        if sa2 + sb2 == 0 or sa2 + sc2 == 0 or sb2 + sd2 == 0 or sc2 + sd2 == 0:
            b_or = float("inf") if sa2 > 0 else 0
            b_p = 1.0
        else:
            b_or, b_p = fisher_test(sa2, sb2, sc2, sd2)
        b_phi = phi_coefficient(sa2, sb2, sc2, sd2)

        print(f"\n  {label} (n={len(subset)}):")
        print(f"    Both: {sa2}, Phytase only: {sb2}, Siderophore only: {sc2}, Neither: {sd2}")
        print(f"    OR: {b_or:.2f}, phi: {b_phi:.4f}, p: {b_p:.2e}")
        print(f"    Phytase prevalence: {len(b_phyt)}/{len(subset)} = {len(b_phyt)/len(subset):.1%}")
        print(f"    Siderophore prevalence: {len(b_sid)}/{len(subset)} = {len(b_sid)/len(subset):.1%}")

    # ── Step 7: Family-level over-representation ──
    print("\nStep 7: Family-level over-representation among linked-trait species")

    linked_species = has_both
    family_counts_linked = defaultdict(int)
    family_counts_all = defaultdict(int)

    for s in all_species:
        fam = taxonomy.get(s, {}).get("family", "unknown")
        family_counts_all[fam] += 1
        if s in linked_species:
            family_counts_linked[fam] += 1

    n_linked = len(linked_species)
    top_families = sorted(family_counts_linked.items(), key=lambda x: -x[1])[:20]

    print(f"\n  Linked-trait species: {n_linked}")
    print(f"  Top 20 families among linked-trait species:")

    family_results = []
    for fam, cnt in top_families:
        total = family_counts_all[fam]
        pct_linked = cnt / n_linked * 100
        pct_of_family = cnt / total * 100 if total > 0 else 0
        print(f"    {fam}: {cnt}/{n_linked} ({pct_linked:.1f}% of linked), {cnt}/{total} ({pct_of_family:.1f}% of family)")
        family_results.append({
            "family": fam,
            "n_linked": cnt,
            "n_total": total,
            "pct_of_linked": round(pct_linked, 2),
            "pct_of_family": round(pct_of_family, 2),
        })

    burkh_linked = len(linked_species & burkh_species)
    burkh_total = len(burkh_species)
    print(f"\n  Burkholderiaceae in linked-trait: {burkh_linked}/{n_linked} ({burkh_linked/n_linked*100:.1f}%)")
    print(f"  Wang 2021 reference: 251/277 = 90.6% Burkholderiaceae")

    # ── Step 8: Save results ──
    print("\nStep 8: Save results")

    summary = {
        "n_total": n_total,
        "n_phytase": len(has_phyt),
        "n_siderophore": len(has_sid),
        "n_both": a,
        "n_phytase_only": b,
        "n_siderophore_only": c,
        "n_neither": d,
        "OR": or_val,
        "phi": phi_val,
        "jaccard": jacc,
        "fisher_p": p_val,
        "perm_z": perm_z,
        "perm_p": perm_p,
    }

    with open(os.path.join(output_dir, "pan_bacterial_summary.csv"), "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["metric", "value"])
        for k, v in summary.items():
            writer.writerow([k, v])

    pd.DataFrame(sid_type_results).to_csv(
        os.path.join(output_dir, "per_siderophore_type.csv"), index=False
    )

    pd.DataFrame(family_results).to_csv(
        os.path.join(output_dir, "family_overrepresentation.csv"), index=False
    )

    # Save species-level data
    species_data = []
    for s in all_species:
        tax = taxonomy.get(s, {})
        species_data.append({
            "species_id": s,
            "has_phytase": int(s in has_phyt),
            "has_siderophore": int(s in has_sid),
            "has_both": int(s in has_both),
            "family": tax.get("family", ""),
            "genus": tax.get("genus", ""),
            "phylum": tax.get("phylum", ""),
        })

    pd.DataFrame(species_data).to_csv(
        os.path.join(output_dir, "species_phytase_siderophore.csv"), index=False
    )

    print(f"\nDone. Results in {output_dir}/")
    spark.stop()


if __name__ == "__main__":
    main()
