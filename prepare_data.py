"""
Cohort filtering for the TCGA-OV HGSOC analyses
===============================================
Reads the raw TCGA-OV tables and cuts them down to confirmed high-grade serous
primary tumours with usable follow-up.

Inputs:
    data/TCGA_OV_TPM.csv        TCGA-OV expression matrix, TPM
    data/TCGA_OV_clinical.csv   TCGA-OV clinical table

Author: Himanshu Gupta, UC San Diego
"""

import os

import pandas as pd

DATA_DIR = "data"

EXPR_PATH = os.path.join(DATA_DIR, "TCGA_OV_TPM.csv")
CLINICAL_PATH = os.path.join(DATA_DIR, "TCGA_OV_clinical.csv")


def load_raw():
    """Read the expression matrix and the clinical table."""
    print("Loading raw data")
    expr = pd.read_csv(EXPR_PATH, index_col=0)
    print(f"  Expression: {expr.shape[0]} genes x {expr.shape[1]} samples")

    clinical = pd.read_csv(CLINICAL_PATH, index_col=0)
    print(f"  Clinical:   {len(clinical)} samples x {clinical.shape[1]} variables")

    return expr, clinical


def filter_to_hgsoc(clinical):
    """Restrict the cohort to confirmed HGSOC primary tumors with follow-up.

    Each step is printed so the attrition is auditable.
    """
    print("\nFiltering to confirmed HGSOC patients")
    print(f"  Starting samples:            {len(clinical)}")

    clinical = clinical[clinical["sample_type"] == "Primary Tumor"]
    print(f"  After primary tumor filter:  {len(clinical)}")

    clinical = clinical[
        clinical["primary_diagnosis"] == "Serous cystadenocarcinoma, NOS"
    ]
    print(f"  After HGSOC restriction:     {len(clinical)}")

    # Very short follow-up carries no survival information.
    clinical = clinical[
        (clinical["days_to_last_follow_up"] >= 30) | (clinical["days_to_death"] >= 30)
    ]
    print(f"  After minimum follow-up:     {len(clinical)}")

    clinical = clinical[
        clinical["follow_ups_disease_response"].isin(["WT-With Tumor", "TF-Tumor Free"])
    ]
    print(f"  After disease response:      {len(clinical)}")

    return clinical


def main():
    expr, clinical = load_raw()
    clinical = filter_to_hgsoc(clinical)

    matched = [s for s in clinical.index if s in expr.columns]
    print(f"\nPatients with expression data: {len(matched)}")


if __name__ == "__main__":
    main()
