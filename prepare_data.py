"""
Cohort filtering for the TCGA-OV HGSOC analyses
===============================================
Reads the raw TCGA-OV tables, cuts them down to confirmed high-grade serous
primary tumours with usable follow-up, and attaches platinum response.

Platinum-free interval is not in the GDC clinical table, so it comes from the
supplementary table of the TCGA ovarian paper (Nature 2011).

Inputs:
    data/TCGA_OV_TPM.csv        TCGA-OV expression matrix, TPM
    data/TCGA_OV_clinical.csv   TCGA-OV clinical table
    data/Table_S1_2.xlsx        Platinum-free interval, TCGA Nature 2011

Author: Himanshu Gupta, UC San Diego
"""

import os

import pandas as pd

DATA_DIR = "data"

EXPR_PATH = os.path.join(DATA_DIR, "TCGA_OV_TPM.csv")
CLINICAL_PATH = os.path.join(DATA_DIR, "TCGA_OV_clinical.csv")
PLATINUM_PATH = os.path.join(DATA_DIR, "Table_S1_2.xlsx")


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


def main():
    expr, clinical, platinum = load_raw()

    clinical = filter_to_hgsoc(clinical)
    clinical = merge_platinum_status(clinical, platinum)

    matched = [s for s in clinical.index if s in expr.columns]
    print(f"\nPatients with expression data: {len(matched)}")

    tiers = clinical["resistance_tier"].value_counts()
    for tier in ["Resistant", "Partially_Sensitive", "Sensitive"]:
        print(f"  {tier:<20} {tiers.get(tier, 0)}")


if __name__ == "__main__":
    main()
