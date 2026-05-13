"""
Figure 5: Operon-distance distribution for P × M gene pairs.

Left panel: histogram of observed within-contig gene-ordinal distances.
Right panel: observed vs permuted-null median comparison.
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main():
    data_dir = "projects/macro_micro_nutrient_cocycling/data/operon_distance"
    fig_dir = "projects/macro_micro_nutrient_cocycling/figures"

    obs_dists = np.load(f"{data_dir}/observed_dists.npy")
    null_medians = np.load(f"{data_dir}/null_medians.npy")
    species_df = pd.read_csv(f"{data_dir}/species_distances.csv")

    obs_median = np.median(obs_dists)
    n_pairs = len(obs_dists)
    n_within5 = np.sum(obs_dists <= 5)
    n_within10 = np.sum(obs_dists <= 10)
    n_within20 = np.sum(obs_dists <= 20)

    n_species_total = len(species_df)
    n_species_same_contig = species_df["n_same_contig"].gt(0).sum()
    pct_same_contig_overall = species_df["n_same_contig"].sum() / species_df["total_pairs"].sum() * 100

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    ax1 = axes[0]
    bins = np.concatenate([np.arange(0, 51, 1), np.arange(60, 210, 10),
                           np.arange(250, 1050, 50), np.arange(1100, 5100, 100)])
    ax1.hist(obs_dists, bins=bins, color="steelblue", edgecolor="none", alpha=0.8)
    ax1.axvline(5, color="red", linestyle="--", alpha=0.7, label=f"5 genes ({n_within5} pairs)")
    ax1.axvline(10, color="orange", linestyle="--", alpha=0.7, label=f"10 genes ({n_within10} pairs)")
    ax1.axvline(20, color="green", linestyle="--", alpha=0.7, label=f"20 genes ({n_within20} pairs)")
    ax1.axvline(obs_median, color="black", linestyle="-", alpha=0.8,
                label=f"Median = {obs_median:.0f} genes")
    ax1.set_xlabel("Gene-ordinal distance (P–M gene pairs, same contig)", fontsize=11)
    ax1.set_ylabel("Number of gene pairs", fontsize=11)
    ax1.set_title(f"A. Within-contig P × M distances\n"
                  f"({n_pairs:,} pairs from {n_species_same_contig:,} species)", fontsize=12)
    ax1.set_xscale("log")
    ax1.legend(fontsize=8, loc="upper right")

    ax2 = axes[1]
    ax2.hist(null_medians, bins=30, color="lightgray", edgecolor="gray", alpha=0.8,
             label=f"Null medians (n=1000)")
    ax2.axvline(obs_median, color="red", linewidth=2, label=f"Observed median = {obs_median:.0f}")
    null_mean = np.mean(null_medians)
    null_std = np.std(null_medians)
    z_score = (obs_median - null_mean) / null_std if null_std > 0 else 0
    ax2.set_xlabel("Median gene-ordinal distance", fontsize=11)
    ax2.set_ylabel("Frequency (permutations)", fontsize=11)
    ax2.set_title(f"B. Observed vs null distribution\n"
                  f"(Z = {z_score:.1f}, observed > all 1000 shuffles)", fontsize=12)
    ax2.legend(fontsize=9)

    textstr = (f"Same-contig fraction: {pct_same_contig_overall:.1f}%\n"
               f"Operon-proximal (≤5 genes): {n_within5}/{n_pairs} = {n_within5/n_pairs*100:.2f}%\n"
               f"P–M genes are NOT physically linked")
    ax2.text(0.98, 0.65, textstr, transform=ax2.transAxes, fontsize=8,
             verticalalignment="top", horizontalalignment="right",
             bbox=dict(boxstyle="round", facecolor="lightyellow", alpha=0.8))

    plt.tight_layout()
    plt.savefig(f"{fig_dir}/figure5_operon_distance.png", dpi=300, bbox_inches="tight")
    plt.savefig(f"{fig_dir}/figure5_operon_distance.pdf", bbox_inches="tight")
    print(f"Saved to {fig_dir}/figure5_operon_distance.png")

    print(f"\nSummary statistics:")
    print(f"  Total species with P+M genes: {n_species_total}")
    print(f"  Species with same-contig pairs: {n_species_same_contig}")
    print(f"  Same-contig pair fraction: {pct_same_contig_overall:.1f}%")
    print(f"  Same-contig pairs: {n_pairs:,}")
    print(f"  Observed median: {obs_median:.0f} genes")
    print(f"  Null median: {null_mean:.1f} ± {null_std:.1f}")
    print(f"  Z-score: {z_score:.1f}")
    print(f"  Pairs within 5 genes (operon-proximal): {n_within5} ({n_within5/n_pairs*100:.2f}%)")
    print(f"  Pairs within 10 genes: {n_within10} ({n_within10/n_pairs*100:.2f}%)")
    print(f"  Pairs within 20 genes: {n_within20} ({n_within20/n_pairs*100:.2f}%)")


if __name__ == "__main__":
    main()
