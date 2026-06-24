"""
Differential lncRNA Expression in Platinum-Resistant HGSOC
=========================================================
Biological question: which lncRNAs are expressed differently at diagnosis
between patients who go on to become platinum-resistant and those who remain
platinum-sensitive?

Method: Mann-Whitney U per lncRNA, resistant against sensitive.

Run prepare_data.py first.

Author: Himanshu Gupta, UC San Diego
"""

import os

import pandas as pd
from scipy.stats import mannwhitneyu

DERIVED_DIR = os.path.join("..", "data", "derived")


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

    results = pd.DataFrame(
        {
            "gene_id": expr.index,
            "log2_fold_change": fold_changes,
            "pvalue": pvalues,
        }
    ).sort_values("pvalue")

    return results


def main():
    expr, clinical = load_data()
    results = run_differential_expression(expr, clinical)

    print(f"\nTested:                    {len(results)}")
    print(f"Significant at raw p<0.05: {int((results['pvalue'] < 0.05).sum())}")

    print("\nTop 10 candidates:")
    for _, row in results.head(10).iterrows():
        direction = "up" if row["log2_fold_change"] > 0 else "down"
        print(
            f"  {row['gene_id']:<20} log2FC={row['log2_fold_change']:+.3f} "
            f"({direction} in resistant)  p={row['pvalue']:.2e}"
        )


if __name__ == "__main__":
    main()
