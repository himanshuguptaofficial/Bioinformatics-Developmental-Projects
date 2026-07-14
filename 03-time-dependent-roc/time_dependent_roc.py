"""
Time-Dependent ROC: lncRNA Signature vs Clinical Variables
==========================================================
Biological question: does the 5-lncRNA risk score predict survival better
than the clinical variables an oncologist already has, namely FIGO stage and
age at diagnosis?

Method: binary survival outcome at 1, 3, and 5 years, ROC and AUC for each
        predictor. Patients censored before a timepoint are excluded at that
        timepoint because their status there is unknown.
Output: three-panel ROC figure.

Run prepare_data.py first.

Author: Himanshu Gupta, UC San Diego
"""

import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from sklearn.metrics import auc, roc_curve

DERIVED_DIR = os.path.join("..", "data", "derived")
FIGURE_DIR = "figures"

# Same signature as project 02; see that script for provenance.
SIGNATURE_LNCRNAS = [
    "ENSG00000251665",
    "ENSG00000258572",  # SYNE3-AS1
    "ENSG00000241912",
    "ENSG00000260505",
    "ENSG00000287733",
]

TIMEPOINTS = [(365, "1-Year"), (1095, "3-Year"), (1825, "5-Year")]

# FIGO substages collapse to their parent stage; the ordinal scale is what
# matters for discrimination.
FIGO_TO_ORDINAL = {
    "Stage I": 1, "Stage IA": 1, "Stage IB": 1, "Stage IC": 1,
    "Stage II": 2, "Stage IIA": 2, "Stage IIB": 2, "Stage IIC": 2,
    "Stage III": 3, "Stage IIIA": 3, "Stage IIIB": 3, "Stage IIIC": 3,
    "Stage IV": 4,
}


def compute_risk_scores(expr, clinical):
    """Refit the project 02 Cox model so this analysis stands alone."""
    base_ids = expr.index.str.split(".").str[0]
    signature = [expr.index[base_ids == gene][0] for gene in SIGNATURE_LNCRNAS]

    values = expr.loc[signature, clinical.index].T
    standardized = (values - values.mean()) / values.std()

    cox_data = standardized.copy()
    cox_data["age"] = clinical["age_years"]
    cox_data["os_time"] = clinical["os_time"]
    cox_data["os_event"] = clinical["os_event"]
    cox_data = cox_data.dropna()

    model = CoxPHFitter()
    model.fit(cox_data, duration_col="os_time", event_col="os_event")

    return model.predict_partial_hazard(cox_data)


def build_comparison_frame(clinical, risk_scores):
    """Align the three predictors on a single complete-case cohort.

    Comparing AUCs computed on different subsets would not be meaningful, so
    every predictor is evaluated on the same patients. FIGO stage is the
    limiting variable.
    """
    frame = pd.DataFrame(
        {
            "os_time": clinical["os_time"],
            "os_event": clinical["os_event"],
            "risk_score": risk_scores,
            "age": clinical["age_years"],
            "figo": clinical["figo_stage"].map(FIGO_TO_ORDINAL),
        }
    )

    before = len(frame)
    frame = frame.dropna()
    print(f"Complete-case cohort: {len(frame)} of {before} patients")
    print(f"  Deaths: {int(frame['os_event'].sum())}")
    print("  FIGO distribution:")
    for stage, count in sorted(frame["figo"].value_counts().items()):
        print(f"    stage {int(stage)}: {count}")

    return frame


def binary_outcome_at(frame, horizon):
    """Label each patient dead/alive at a horizon, dropping the unknowable.

    A patient censored before the horizon has an unknown status there, so
    including them would bias the estimate either way.
    """
    died_before = (frame["os_event"] == 1) & (frame["os_time"] <= horizon)
    followed_beyond = frame["os_time"] > horizon

    evaluable = died_before | followed_beyond
    return died_before[evaluable].astype(int), evaluable


def plot_roc_panels(frame):
    predictors = [
        ("lncRNA signature", "risk_score", "steelblue"),
        ("Age at diagnosis", "age", "forestgreen"),
        ("FIGO stage", "figo", "firebrick"),
    ]

    fig, axes = plt.subplots(1, 3, figsize=(17, 5.8))
    auc_table = {}

    for ax, (horizon, label) in zip(axes, TIMEPOINTS):
        outcome, evaluable = binary_outcome_at(frame, horizon)
        n_events = int(outcome.sum())

        print(f"\n{label} (day {horizon}): {evaluable.sum()} evaluable, {n_events} deaths")

        for name, column, color in predictors:
            scores = frame.loc[evaluable, column]
            fpr, tpr, _ = roc_curve(outcome, scores)
            roc_auc = auc(fpr, tpr)
            auc_table.setdefault(name, {})[label] = roc_auc

            ax.plot(fpr, tpr, color=color, lw=2, label=f"{name} (AUC={roc_auc:.3f})")
            print(f"  {name:<20} AUC={roc_auc:.3f}")

        ax.plot([0, 1], [0, 1], "k--", lw=1, alpha=0.6, label="Random (0.500)")
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1.02)
        ax.set_xlabel("False positive rate", fontsize=10)
        ax.set_ylabel("True positive rate", fontsize=10)
        ax.set_title(f"{label} survival (n={evaluable.sum()}, {n_events} deaths)", fontsize=11)
        ax.legend(loc="lower right", fontsize=9)
        ax.grid(alpha=0.2)

    fig.suptitle(
        "Time-Dependent ROC: 5-lncRNA Signature vs Clinical Variables\n"
        f"TCGA-OV HGSOC, n={len(frame)} with complete clinical data",
        fontsize=13,
    )
    fig.tight_layout()

    path = os.path.join(FIGURE_DIR, "time_dependent_roc.png")
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"\nSaved {path}")

    return auc_table


def main():
    os.makedirs(FIGURE_DIR, exist_ok=True)

    expr = pd.read_csv(os.path.join(DERIVED_DIR, "expr_lncrna_final.csv"), index_col=0)
    clinical = pd.read_csv(
        os.path.join(DERIVED_DIR, "clinical_hgsoc_filtered.csv"), index_col=0
    )

    risk_scores = compute_risk_scores(expr, clinical)
    frame = build_comparison_frame(clinical, risk_scores)
    auc_table = plot_roc_panels(frame)

    print("\nAUC summary:")
    header = f"{'Predictor':<20}" + "".join(f"{label:>10}" for _, label in TIMEPOINTS)
    print(header)
    print("-" * len(header))
    for name, by_time in auc_table.items():
        row = f"{name:<20}" + "".join(
            f"{by_time[label]:>10.3f}" for _, label in TIMEPOINTS
        )
        print(row)


if __name__ == "__main__":
    main()
