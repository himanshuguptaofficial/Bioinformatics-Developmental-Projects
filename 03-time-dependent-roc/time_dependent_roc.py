"""
Time-Dependent ROC for the 5-lncRNA Risk Score
==============================================
Biological question: how well does the 5-lncRNA risk score predict survival,
and does that change with the prediction horizon?

Method: binary survival outcome at 1, 3, and 5 years, ROC and AUC at each.
        Patients censored before a timepoint are excluded at that timepoint
        because their status there is unknown.

Run prepare_data.py first.

Author: Himanshu Gupta, UC San Diego
"""

import os

import pandas as pd
from lifelines import CoxPHFitter
from sklearn.metrics import auc, roc_curve

DERIVED_DIR = os.path.join("..", "data", "derived")

# Same signature as project 02; see that script for provenance.
SIGNATURE_LNCRNAS = [
    "ENSG00000251665",
    "ENSG00000258572",  # SYNE3-AS1
    "ENSG00000241912",
    "ENSG00000260505",
    "ENSG00000287733",
]

TIMEPOINTS = [(365, "1-Year"), (1095, "3-Year"), (1825, "5-Year")]


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


def binary_outcome_at(frame, horizon):
    """Label each patient dead/alive at a horizon, dropping the unknowable.

    A patient censored before the horizon has an unknown status there, so
    including them would bias the estimate either way.
    """
    died_before = (frame["os_event"] == 1) & (frame["os_time"] <= horizon)
    followed_beyond = frame["os_time"] > horizon

    evaluable = died_before | followed_beyond
    return died_before[evaluable].astype(int), evaluable


def main():
    expr = pd.read_csv(os.path.join(DERIVED_DIR, "expr_lncrna_final.csv"), index_col=0)
    clinical = pd.read_csv(
        os.path.join(DERIVED_DIR, "clinical_hgsoc_filtered.csv"), index_col=0
    )

    risk_scores = compute_risk_scores(expr, clinical)

    frame = pd.DataFrame(
        {
            "os_time": clinical["os_time"],
            "os_event": clinical["os_event"],
            "risk_score": risk_scores,
        }
    ).dropna()
    print(f"Patients with a risk score: {len(frame)}")

    for horizon, label in TIMEPOINTS:
        outcome, evaluable = binary_outcome_at(frame, horizon)
        fpr, tpr, _ = roc_curve(outcome, frame.loc[evaluable, "risk_score"])

        print(
            f"{label} (day {horizon}): {evaluable.sum()} evaluable, "
            f"{int(outcome.sum())} deaths, AUC={auc(fpr, tpr):.3f}"
        )


if __name__ == "__main__":
    main()
