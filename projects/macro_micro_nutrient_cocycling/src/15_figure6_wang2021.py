"""
Figure 6: Wang 2021 validation — phytase × siderophore co-occurrence.

Panel A: Per-family co-occurrence rates, Burkholderiaceae highlighted.
Panel B: Per-siderophore-type co-occurrence with phytase.
"""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main():
    data_dir = "projects/macro_micro_nutrient_cocycling/data/wang2021"
    fig_dir = "projects/macro_micro_nutrient_cocycling/figures"

    fam_df = pd.read_csv(f"{data_dir}/family_overrepresentation.csv")
    sid_df = pd.read_csv(f"{data_dir}/per_siderophore_type.csv")
    summary = {}
    with open(f"{data_dir}/pan_bacterial_summary.csv") as f:
        import csv
        reader = csv.reader(f)
        next(reader)
        for row in reader:
            summary[row[0]] = float(row[1])

    fig, axes = plt.subplots(1, 2, figsize=(15, 7), gridspec_kw={"width_ratios": [2, 1]})

    # Panel A: Family-level co-occurrence
    ax1 = axes[0]
    fam_df = fam_df.sort_values("pct_of_family", ascending=True)

    colors = ["#d62728" if "Burkholderiaceae" in f else "steelblue"
              for f in fam_df["family"]]
    labels = [f.replace("f__", "") for f in fam_df["family"]]

    bars = ax1.barh(range(len(fam_df)), fam_df["pct_of_family"], color=colors,
                    edgecolor="none", alpha=0.8)

    ax1.set_yticks(range(len(fam_df)))
    ax1.set_yticklabels(labels, fontsize=8)
    ax1.set_xlabel("% of family with both phytase + siderophore", fontsize=11)
    ax1.set_title("A. Phytase × siderophore co-occurrence by family\n"
                  "(red = Burkholderiaceae sensu lato)", fontsize=12)

    for i, (_, row) in enumerate(fam_df.iterrows()):
        ax1.text(row["pct_of_family"] + 0.5, i,
                 f'{row["n_linked"]}/{row["n_total"]}',
                 va="center", fontsize=7, color="gray")

    # Panel B: Per-siderophore-type
    ax2 = axes[1]
    sid_df = sid_df.sort_values("OR", ascending=True)
    colors_sid = plt.cm.Set2(np.linspace(0, 0.8, len(sid_df)))

    bars2 = ax2.barh(range(len(sid_df)), sid_df["OR"], color=colors_sid,
                     edgecolor="none", alpha=0.9)

    ax2.set_yticks(range(len(sid_df)))
    ax2.set_yticklabels(sid_df["siderophore_type"], fontsize=10)
    ax2.set_xlabel("Odds ratio (phytase co-occurrence)", fontsize=11)
    ax2.axvline(1.0, color="black", linestyle="--", alpha=0.3)
    ax2.set_title("B. Co-occurrence by siderophore type", fontsize=12)

    for i, (_, row) in enumerate(sid_df.iterrows()):
        sig = "***" if row["fisher_p"] < 0.001 else "**" if row["fisher_p"] < 0.01 else "*" if row["fisher_p"] < 0.05 else "ns"
        ax2.text(row["OR"] + 0.02, i,
                 f'n={row["n_species"]}, {sig}',
                 va="center", fontsize=9)

    plt.tight_layout()
    plt.savefig(f"{fig_dir}/figure6_wang2021.png", dpi=300, bbox_inches="tight")
    plt.savefig(f"{fig_dir}/figure6_wang2021.pdf", bbox_inches="tight")
    print(f"Saved to {fig_dir}/figure6_wang2021.png")

    print(f"\nPan-bacterial summary:")
    print(f"  OR = {summary['OR']:.3f}, phi = {summary['phi']:.4f}")
    print(f"  Fisher p = {summary['fisher_p']:.2e}")
    print(f"  Both traits: {int(summary['n_both'])}/{int(summary['n_total'])} species")
    print(f"\nPhylogenetically corrected:")

    wang_phylo = pd.read_csv("projects/macro_micro_nutrient_cocycling/data/phylo/phylo_logistic_wang.csv")
    print(f"  Uncorrected log-OR: {wang_phylo['uncorrected_logOR'].iloc[0]:.3f}")
    print(f"  Corrected log-OR: {wang_phylo['phylo_logOR'].iloc[0]:.3f}")
    print(f"  Corrected p: {wang_phylo['phylo_p'].iloc[0]:.2e}")
    reduction = (1 - abs(wang_phylo['phylo_logOR'].iloc[0]) / abs(wang_phylo['uncorrected_logOR'].iloc[0])) * 100
    print(f"  Reduction: {reduction:.0f}%")


if __name__ == "__main__":
    main()
