"""
Differential lncRNA Expression in Platinum-Resistant HGSOC
=========================================================
Biological question: which lncRNAs are expressed differently at diagnosis
between patients who go on to become platinum-resistant and those who remain
platinum-sensitive?

Method: Mann-Whitney U per lncRNA, Benjamini-Hochberg FDR correction.
Output: volcano plot with the strongest candidates labelled, plus a ranked
        table of every lncRNA tested.

Run prepare_data.py first.

Author: Himanshu Gupta, UC San Diego
"""

import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

DERIVED_DIR = os.path.join("..", "data", "derived")
FIGURE_DIR = "figures"
RESULT_DIR = "results"

# Genes this far from the origin are worth naming on the plot.
N_LABELLED = 10


def load_data():
    expr = pd.read_csv(os.path.join(DERIVED_DIR, "expr_lncrna_final.csv"), index_col=0)
    clinical = pd.read_csv(
        os.path.join(DERIVED_DIR, "clinical_hgsoc_filtered.csv"), index_col=0
    )

    print(f"lncRNAs: {expr.shape[0]}, patients: {expr.shape[1]}")
    print("Resistance tiers:")
    for tier, count in clinical["resistance_tier"].value_counts().items():
        print(f"  {tier:<20} {count}")

    return expr, clinical


def run_differential_expression(expr, clinical):
    """Compare resistant against sensitive patients, one test per lncRNA.

    Partially sensitive patients (PFI 6-12 months) are deliberately excluded:
    the contrast is only meaningful at the two extremes of platinum response.
    """
    resistant = clinical.index[clinical["resistance_tier"] == "Resistant"]
    sensitive = clinical.index[clinical["resistance_tier"] == "Sensitive"]

    print(f"\nComparing {len(resistant)} resistant vs {len(sensitive)} sensitive")

    resistant_expr = expr[resistant].to_numpy()
    sensitive_expr = expr[sensitive].to_numpy()

    # axis=1 runs one test per row, far faster than looping over 12k genes.
    _, pvalues = mannwhitneyu(
        resistant_expr, sensitive_expr, alternative="two-sided", axis=1
    )

    # Expression is already log2(TPM+1), so a difference of means is a log2
    # fold change.
    fold_changes = resistant_expr.mean(axis=1) - sensitive_expr.mean(axis=1)

    _, adj_pvalues, _, _ = multipletests(pvalues, method="fdr_bh")

    results = pd.DataFrame(
        {
            "gene_id": expr.index,
            "log2_fold_change": fold_changes,
            "pvalue": pvalues,
            "adj_pvalue": adj_pvalues,
        }
    ).sort_values("pvalue")

    return results


def summarize(results):
    n_raw = int((results["pvalue"] < 0.05).sum())
    n_fdr = int((results["adj_pvalue"] < 0.05).sum())
    n_up = int(((results["pvalue"] < 0.05) & (results["log2_fold_change"] > 0)).sum())
    n_down = int(((results["pvalue"] < 0.05) & (results["log2_fold_change"] < 0)).sum())

    print(f"\nTested:                    {len(results)}")
    print(f"Significant at raw p<0.05: {n_raw}")
    print(f"  Up in resistant:         {n_up}")
    print(f"  Down in resistant:       {n_down}")
    print(f"Significant after FDR:     {n_fdr}")
    print(f"Smallest adjusted p-value: {results['adj_pvalue'].min():.4f}")

    print("\nTop 10 candidates:")
    top = results.head(10)
    for _, row in top.iterrows():
        direction = "up" if row["log2_fold_change"] > 0 else "down"
        print(
            f"  {row['gene_id']:<20} log2FC={row['log2_fold_change']:+.3f} "
            f"({direction} in resistant)  p={row['pvalue']:.2e}"
        )

    return n_up, n_down


def plot_volcano(results, n_up, n_down, resistant_n, sensitive_n):
    """Volcano plot, with the most significant lncRNAs labelled by name."""
    significant = results["pvalue"] < 0.05
    upregulated = significant & (results["log2_fold_change"] > 0)
    downregulated = significant & (results["log2_fold_change"] < 0)

    colors = np.full(len(results), "lightgrey", dtype=object)
    colors[upregulated.to_numpy()] = "firebrick"
    colors[downregulated.to_numpy()] = "steelblue"

    fig, ax = plt.subplots(figsize=(10, 8))
    ax.scatter(
        results["log2_fold_change"],
        -np.log10(results["pvalue"]),
        c=colors,
        alpha=0.55,
        s=12,
        linewidths=0,
    )

    ax.axhline(-np.log10(0.05), color="black", linestyle="--", alpha=0.7, linewidth=1)
    ax.axvline(0, color="black", linestyle="-", alpha=0.2, linewidth=1)

    # Name the strongest hits so the plot carries information on its own.
    for _, row in results.head(N_LABELLED).iterrows():
        ax.annotate(
            row["gene_id"].split(".")[0],
            xy=(row["log2_fold_change"], -np.log10(row["pvalue"])),
            xytext=(4, 3),
            textcoords="offset points",
            fontsize=7,
            color="black",
        )

    ax.set_xlabel(
        "log2 fold change, resistant vs sensitive\n(mean difference in log2(TPM+1))",
        fontsize=11,
    )
    ax.set_ylabel("-log10(p-value)", fontsize=11)
    ax.set_title(
        "Differential lncRNA Expression in Platinum-Resistant HGSOC\n"
        f"TCGA-OV | {resistant_n} resistant vs {sensitive_n} sensitive",
        fontsize=12,
    )

    ax.legend(
        handles=[
            Patch(facecolor="firebrick", label=f"Up in resistant (n={n_up})"),
            Patch(facecolor="steelblue", label=f"Down in resistant (n={n_down})"),
            Patch(facecolor="lightgrey", label="Not significant"),
        ],
        fontsize=9,
        loc="upper left",
    )
    ax.grid(alpha=0.15)

    fig.tight_layout()
    path = os.path.join(FIGURE_DIR, "volcano_plot.png")
    fig.savefig(path, dpi=300)
    plt.close(fig)
    print(f"\nSaved {path}")


def main():
    os.makedirs(FIGURE_DIR, exist_ok=True)
    os.makedirs(RESULT_DIR, exist_ok=True)

    expr, clinical = load_data()
    results = run_differential_expression(expr, clinical)
    n_up, n_down = summarize(results)

    result_path = os.path.join(RESULT_DIR, "de_results.csv")
    results.to_csv(result_path, index=False)
    print(f"Saved {result_path}")

    plot_volcano(
        results,
        n_up,
        n_down,
        resistant_n=int((clinical["resistance_tier"] == "Resistant").sum()),
        sensitive_n=int((clinical["resistance_tier"] == "Sensitive").sum()),
    )


if __name__ == "__main__":
    main()
