"""
Figure 4: Uncorrected vs phylogenetically corrected effect sizes.

17 tested pairs color-coded by mean phylogenetic signal strength
(average Pagel's lambda of the two genes in each pair).
"""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


def main():
    data_dir = "projects/macro_micro_nutrient_cocycling/data/phylo"
    fig_dir = "projects/macro_micro_nutrient_cocycling/figures"

    logistic = pd.read_csv(f"{data_dir}/phylo_logistic.csv")
    signal = pd.read_csv(f"{data_dir}/phylo_signal.csv")

    lambda_map = dict(zip(signal["gene"], signal["lambda"]))

    group_lambda = {
        "P-acquisition": np.nanmean([lambda_map.get(g, np.nan)
                                      for g in ["phoA", "pstA", "pstB", "pstC", "pstS",
                                                 "phnC", "phnD", "phnE", "phoD"]]),
        "N-fixation": np.nanmean([lambda_map.get(g, np.nan) for g in ["nifH", "nifD"]]),
        "Metal-handling": np.nanmean([lambda_map.get(g, np.nan)
                                       for g in ["copA", "corA", "feoB", "HMA"]]),
        "Phz-operon": np.nanmean([lambda_map.get(g, np.nan)
                                   for g in ["phzF", "phzA", "phzB", "phzD", "phzG", "phzS", "phzM"]]),
    }

    def get_pair_lambda(row):
        gx = row["gene_x"]
        gy = row["gene_y"]

        lx = lambda_map.get(gx, group_lambda.get(gx, np.nan))
        ly = lambda_map.get(gy, group_lambda.get(gy, np.nan))

        vals = [v for v in [lx, ly] if not np.isnan(v)]
        return np.mean(vals) if vals else 0.5

    logistic["mean_lambda"] = logistic.apply(get_pair_lambda, axis=1)

    sig_pairs = logistic[logistic["phylo_p"] < 0.05].copy()
    nonsig_pairs = logistic[logistic["phylo_p"] >= 0.05].copy()

    fig, ax = plt.subplots(figsize=(10, 8))

    cmap = plt.cm.RdYlBu_r
    norm = matplotlib.colors.Normalize(vmin=0.3, vmax=1.0)

    for _, row in sig_pairs.iterrows():
        x = row["uncorrected_logOR"]
        y = row["phylo_logOR"]
        if np.isinf(x) or np.isinf(y):
            continue
        color = cmap(norm(row["mean_lambda"]))
        ax.scatter(x, y, c=[color], s=100, edgecolors="black", linewidths=0.5, zorder=5)
        label = row["pair_label"]
        if len(label) > 20:
            label = label.replace(" (group)", "\n(group)")
        offset = (8, 8)
        if "pstC x feoB" in row["pair_label"]:
            offset = (8, -12)
        elif "pstS x feoB" in row["pair_label"]:
            offset = (-60, -12)
        elif "pstC x HMA" in row["pair_label"]:
            offset = (8, -12)
        elif "phzF x feoB" in row["pair_label"]:
            offset = (-60, 8)
        elif "P x N" in row["pair_label"]:
            offset = (8, -12)
        ax.annotate(label, (x, y), textcoords="offset points",
                    xytext=offset, fontsize=7, alpha=0.8)

    for _, row in nonsig_pairs.iterrows():
        x = row["uncorrected_logOR"]
        y = row["phylo_logOR"]
        if np.isinf(x) or np.isinf(y):
            continue
        ax.scatter(x, y, c="lightgray", s=60, marker="x",
                   edgecolors="gray", linewidths=1, zorder=4)

    lim_min = min(ax.get_xlim()[0], ax.get_ylim()[0])
    lim_max = max(ax.get_xlim()[1], ax.get_ylim()[1])
    margin = (lim_max - lim_min) * 0.1
    lim_min -= margin
    lim_max += margin

    ax.plot([lim_min, lim_max], [lim_min, lim_max], "k--", alpha=0.3, linewidth=1, zorder=1)
    ax.set_xlim(lim_min, lim_max)
    ax.set_ylim(lim_min, lim_max)

    ax.set_xlabel("Uncorrected log-OR (Fisher's exact)", fontsize=12)
    ax.set_ylabel("Phylogenetically corrected log-OR (phyloglm)", fontsize=12)
    ax.set_title("Figure 4. Effect sizes before and after phylogenetic correction\n"
                 "(color = mean Pagel's λ of gene pair)", fontsize=13)

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, shrink=0.7, pad=0.02)
    cbar.set_label("Mean Pagel's λ (phylogenetic signal)", fontsize=10)

    ax.axhline(0, color="gray", linestyle=":", alpha=0.3)
    ax.axvline(0, color="gray", linestyle=":", alpha=0.3)

    above = mpatches.Patch(color="none", label="Above diagonal = strengthened by correction")
    below = mpatches.Patch(color="none", label="Below diagonal = weakened by correction")
    ax.legend(handles=[above, below], loc="upper left", fontsize=8,
              framealpha=0.7, handlelength=0)

    ax.set_aspect("equal")
    plt.tight_layout()

    plt.savefig(f"{fig_dir}/figure4_phylo_correction.png", dpi=300, bbox_inches="tight")
    plt.savefig(f"{fig_dir}/figure4_phylo_correction.pdf", bbox_inches="tight")
    print(f"Saved to {fig_dir}/figure4_phylo_correction.png")

    print("\nPairs above diagonal (strengthened):")
    for _, row in sig_pairs.iterrows():
        x = row["uncorrected_logOR"]
        y = row["phylo_logOR"]
        if np.isinf(x) or np.isinf(y):
            continue
        if abs(y) > abs(x):
            pct = (abs(y) - abs(x)) / abs(x) * 100 if x != 0 else float("inf")
            print(f"  {row['pair_label']}: {x:.3f} -> {y:.3f} ({pct:+.0f}%)")

    print("\nPairs below diagonal (weakened):")
    for _, row in sig_pairs.iterrows():
        x = row["uncorrected_logOR"]
        y = row["phylo_logOR"]
        if np.isinf(x) or np.isinf(y):
            continue
        if abs(y) < abs(x):
            pct = (abs(y) - abs(x)) / abs(x) * 100 if x != 0 else float("inf")
            print(f"  {row['pair_label']}: {x:.3f} -> {y:.3f} ({pct:+.0f}%)")


if __name__ == "__main__":
    main()
