"""
Kaplan-Meier Survival Analysis with a 5-lncRNA Risk Score
=========================================================
Biological question: does a composite risk score built from five lncRNAs
stratify HGSOC patients by overall survival?

Method: multivariate Cox regression on the five signature lncRNAs adjusted
        for age, then a median split of the resulting risk score.
Output: Kaplan-Meier curves for the high- and low-risk groups.

Run prepare_data.py first.

Author: Himanshu Gupta, UC San Diego
"""

import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter, KaplanMeierFitter

DERIVED_DIR = os.path.join("..", "data", "derived")
FIGURE_DIR = "figures"

# The five signature lncRNAs.
#
# Provenance: these were selected in the parent study by LASSO-Cox regression
# over the lncRNAs that were both differentially expressed (project 01) and
# associated with survival in univariate Cox models, then ranked by selection
# frequency across 1000 bootstrap resamples. They are hardcoded here so this
# analysis runs in seconds; rederiving them takes 10-20 minutes.
SIGNATURE_LNCRNAS = [
    "ENSG00000251665",
    "ENSG00000258572",  # SYNE3-AS1
    "ENSG00000241912",
    "ENSG00000260505",
    "ENSG00000287733",
]


def resolve_signature(expr):
    """Map unversioned signature IDs onto the versioned IDs in the matrix."""
    base_ids = expr.index.str.split(".").str[0]

    resolved = []
    for gene in SIGNATURE_LNCRNAS:
        matches = expr.index[base_ids == gene].tolist()
        if not matches:
            raise KeyError(f"signature lncRNA {gene} is absent from the matrix")
        resolved.append(matches[0])

    return resolved


def build_cox_data(expr, clinical, signature):
    """Assemble the modelling frame: standardized expression, age, outcome.

    Expression is z-scored so the Cox coefficients are comparable across
    lncRNAs with different dynamic ranges.
    """
    values = expr.loc[signature, clinical.index].T
    standardized = (values - values.mean()) / values.std()

    cox_data = standardized.copy()
    cox_data["age"] = clinical["age_years"]
    cox_data["os_time"] = clinical["os_time"]
    cox_data["os_event"] = clinical["os_event"]

    # Age is missing for a handful of patients; Cox needs complete cases.
    before = len(cox_data)
    cox_data = cox_data.dropna()
    print(f"Patients: {before} total, {len(cox_data)} with complete covariates")

    return cox_data


def fit_model(cox_data):
    model = CoxPHFitter()
    model.fit(cox_data, duration_col="os_time", event_col="os_event")

    print("\nMultivariate Cox model:")
    summary = model.summary[
        ["exp(coef)", "exp(coef) lower 95%", "exp(coef) upper 95%", "p"]
    ]
    for name, row in summary.iterrows():
        print(
            f"  {name:<20} HR={row['exp(coef)']:.3f} "
            f"[{row['exp(coef) lower 95%']:.3f}-{row['exp(coef) upper 95%']:.3f}]  "
            f"p={row['p']:.4f}"
        )

    return model


def stratify(model, cox_data):
    """Split the cohort at the median risk score."""
    risk_scores = model.predict_partial_hazard(cox_data)
    median = risk_scores.median()
    groups = np.where(risk_scores > median, "High", "Low")

    stratified = pd.DataFrame(
        {
            "os_time": cox_data["os_time"],
            "os_event": cox_data["os_event"],
            "risk_score": risk_scores,
            "risk_group": groups,
        }
    )

    print(f"\nHigh risk: {(groups == 'High').sum()} patients")
    print(f"Low risk:  {(groups == 'Low').sum()} patients")

    return stratified


def plot_kaplan_meier(stratified):
    high = stratified[stratified["risk_group"] == "High"]
    low = stratified[stratified["risk_group"] == "Low"]

    fig, ax = plt.subplots(figsize=(10, 6.5))

    fitted = {}
    for data, label, color in [
        (high, f"High risk (n={len(high)})", "firebrick"),
        (low, f"Low risk (n={len(low)})", "steelblue"),
    ]:
        kmf = KaplanMeierFitter()
        kmf.fit(data["os_time"], event_observed=data["os_event"], label=label)
        kmf.plot_survival_function(ax=ax, color=color, ci_show=True, linewidth=2)
        fitted[label.split()[0]] = kmf

    ax.set_title(
        "Overall Survival by 5-lncRNA Risk Score\nTCGA-OV HGSOC cohort",
        fontsize=12,
    )
    ax.set_xlabel("Time (days)", fontsize=11)
    ax.set_ylabel("Overall survival probability", fontsize=11)
    ax.set_ylim(0, 1.02)
    ax.legend(fontsize=10, loc="upper right")
    ax.grid(alpha=0.15)

    fig.tight_layout()
    path = os.path.join(FIGURE_DIR, "kaplan_meier.png")
    fig.savefig(path, dpi=300)
    plt.close(fig)
    print(f"\nSaved {path}")

    return fitted


def main():
    os.makedirs(FIGURE_DIR, exist_ok=True)

    expr = pd.read_csv(os.path.join(DERIVED_DIR, "expr_lncrna_final.csv"), index_col=0)
    clinical = pd.read_csv(
        os.path.join(DERIVED_DIR, "clinical_hgsoc_filtered.csv"), index_col=0
    )

    signature = resolve_signature(expr)
    print(f"Signature lncRNAs: {', '.join(signature)}")

    cox_data = build_cox_data(expr, clinical, signature)
    print(f"Deaths observed: {int(cox_data['os_event'].sum())}")

    model = fit_model(cox_data)
    stratified = stratify(model, cox_data)

    fitted = plot_kaplan_meier(stratified)
    for label, kmf in fitted.items():
        print(f"Median OS {label:<5} {kmf.median_survival_time_:.0f} days")


if __name__ == "__main__":
    main()
