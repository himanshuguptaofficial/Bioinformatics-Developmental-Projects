"""
Shared preprocessing for the TCGA-OV HGSOC analyses
===================================================
Builds the derived matrices that the project scripts read. Run this once after
downloading the raw data.

Inputs:
    data/TCGA_OV_TPM.csv        TCGA-OV expression matrix, TPM
    data/TCGA_OV_clinical.csv   TCGA-OV clinical table
    data/Table_S1_2.xlsx        Platinum-free interval, TCGA Nature 2011

Outputs:
    data/derived/clinical_hgsoc_filtered.csv   216 patients, survival + tier
    data/derived/expr_hgsoc_filtered.csv       all genes, log2(TPM+1)

Author: Himanshu Gupta, UC San Diego
"""

import os

import numpy as np
import pandas as pd

DATA_DIR = "data"
DERIVED_DIR = os.path.join(DATA_DIR, "derived")

EXPR_PATH = os.path.join(DATA_DIR, "TCGA_OV_TPM.csv")
CLINICAL_PATH = os.path.join(DATA_DIR, "TCGA_OV_clinical.csv")
PLATINUM_PATH = os.path.join(DATA_DIR, "Table_S1_2.xlsx")

# A gene must be expressed in at least this fraction of patients to be kept.
MIN_EXPRESSED_FRACTION = 0.2


def load_raw():
    """Read the three input files."""
    print("Loading raw data")
    expr = pd.read_csv(EXPR_PATH, index_col=0)
    print(f"  Expression: {expr.shape[0]} genes x {expr.shape[1]} samples")

    clinical = pd.read_csv(CLINICAL_PATH, index_col=0)
    print(f"  Clinical:   {len(clinical)} samples x {clinical.shape[1]} variables")

    platinum = pd.read_excel(PLATINUM_PATH)
    print(f"  Platinum:   {len(platinum)} patients")

    return expr, clinical, platinum


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


def assign_tier(pfi_months):
    """Three-tier platinum response, GCIG consensus definition."""
    if pfi_months < 6:
        return "Resistant"
    if pfi_months <= 12:
        return "Partially_Sensitive"
    return "Sensitive"


def merge_platinum_status(clinical, platinum):
    """Attach platinum-free interval and derive the resistance tier."""
    print("\nMerging platinum resistance status")

    # Expression barcodes look like TCGA-29-1691-01A-01R-1566-13; the platinum
    # table is keyed on the 12-character patient barcode.
    clinical = clinical.copy()
    clinical["patient_id"] = clinical.index.str[:12]

    platinum_clean = platinum[
        ["BCRPATIENTBARCODE", "PlatinumFreeInterval (mos)*", "PlatinumStatus"]
    ].copy()
    platinum_clean.columns = ["patient_id", "PFI_months", "PlatinumStatus"]

    # merge() drops the index, so restore the full barcode afterwards.
    clinical = clinical.merge(platinum_clean, on="patient_id", how="left")
    clinical = clinical.set_index("barcode")

    clinical = clinical[clinical["PlatinumStatus"].isin(["Sensitive", "Resistant"])]
    print(f"  Definitive platinum status:  {len(clinical)}")

    clinical["resistance_tier"] = clinical["PFI_months"].apply(assign_tier)
    return clinical


def add_survival_variables(clinical):
    """Build overall-survival time, event indicator, and age in years."""
    clinical = clinical.copy()

    # Dead patients contribute time to death; the rest are censored at last
    # follow-up.
    clinical["os_time"] = clinical["days_to_death"].fillna(
        clinical["days_to_last_follow_up"]
    )
    clinical["os_event"] = (clinical["vital_status"] == "Dead").astype(int)
    clinical["age_years"] = clinical["age_at_diagnosis"] / 365.25

    return clinical


def normalize_and_filter(expr, sample_ids):
    """Match samples, apply log2(TPM+1), and drop low-expression genes."""
    print("\nNormalizing expression")

    matched = [s for s in sample_ids if s in expr.columns]
    expr = expr[matched]
    print(f"  Matched samples:             {len(matched)}")

    expr_log = np.log2(expr + 1)

    min_patients = int(MIN_EXPRESSED_FRACTION * expr_log.shape[1])
    keep = (expr_log > 0).sum(axis=1) >= min_patients
    expr_log = expr_log[keep]
    print(f"  Genes after low-expr filter: {expr_log.shape[0]}")

    return expr_log, matched


def main():
    os.makedirs(DERIVED_DIR, exist_ok=True)

    expr, clinical, platinum = load_raw()

    clinical = filter_to_hgsoc(clinical)
    clinical = merge_platinum_status(clinical, platinum)
    clinical = add_survival_variables(clinical)

    expr_log, matched = normalize_and_filter(expr, clinical.index.tolist())
    clinical = clinical.loc[matched]

    print("\nWriting derived files")
    clinical.to_csv(os.path.join(DERIVED_DIR, "clinical_hgsoc_filtered.csv"))
    expr_log.to_csv(os.path.join(DERIVED_DIR, "expr_hgsoc_filtered.csv"))

    print(f"\nFinal cohort: {len(clinical)} patients")
    tiers = clinical["resistance_tier"].value_counts()
    for tier in ["Resistant", "Partially_Sensitive", "Sensitive"]:
        print(f"  {tier:<20} {tiers.get(tier, 0)}")
    print(f"Deaths observed: {clinical['os_event'].sum()}")


if __name__ == "__main__":
    main()
