"""
Shared preprocessing for all four analyses
==========================================
Builds the derived matrices that every project in this repository reads.
Run this once after downloading the raw data; the four project scripts then
run in seconds instead of repeating this cascade.

Inputs (see README for how to obtain each):
    data/TCGA_OV_TPM.csv        TCGA-OV expression matrix, TPM
    data/TCGA_OV_clinical.csv   TCGA-OV clinical table
    data/Table_S1_2.xlsx        Platinum-free interval, TCGA Nature 2011
    data/gencode.v44.gtf.gz     GENCODE v44 annotation

Outputs:
    data/derived/clinical_hgsoc_filtered.csv   216 patients, survival + tier
    data/derived/expr_hgsoc_filtered.csv       all genes, log2(TPM+1)
    data/derived/expr_lncrna_final.csv         lncRNAs only, log2(TPM+1)

Author: Himanshu Gupta, UC San Diego
"""

import gzip
import os

import numpy as np
import pandas as pd

DATA_DIR = "data"
DERIVED_DIR = os.path.join(DATA_DIR, "derived")

EXPR_PATH = os.path.join(DATA_DIR, "TCGA_OV_TPM.csv")
CLINICAL_PATH = os.path.join(DATA_DIR, "TCGA_OV_clinical.csv")
PLATINUM_PATH = os.path.join(DATA_DIR, "Table_S1_2.xlsx")
GTF_PATH = os.path.join(DATA_DIR, "gencode.v44.gtf.gz")

# A gene must be expressed in at least this fraction of patients to be kept.
MIN_EXPRESSED_FRACTION = 0.2


def load_raw():
    """Read the four input files."""
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


def load_lncrna_ids(gtf_path):
    """Collect Ensembl gene IDs annotated as lncRNA in GENCODE."""
    print("\nExtracting lncRNAs from GENCODE v44")

    lncrna_ids = set()
    with gzip.open(gtf_path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue

            gene_type, gene_id = "", ""
            for attr in fields[8].split(";"):
                attr = attr.strip()
                if attr.startswith("gene_type"):
                    gene_type = attr.split('"')[1]
                elif attr.startswith("gene_id"):
                    gene_id = attr.split('"')[1]

            if gene_type == "lncRNA" and gene_id:
                lncrna_ids.add(gene_id)

    print(f"  lncRNA genes in GENCODE:     {len(lncrna_ids)}")
    return lncrna_ids


def extract_lncrnas(expr_log, lncrna_ids):
    """Subset the matrix to lncRNAs, matching on unversioned Ensembl IDs."""
    # GENCODE and TCGA both carry version suffixes but not always the same one.
    lncrna_base = {gene_id.split(".")[0] for gene_id in lncrna_ids}
    expr_base = expr_log.index.str.split(".").str[0]
    expr_lncrna = expr_log[expr_base.isin(lncrna_base)]

    # Re-apply the low-expression filter within the lncRNA subset, which is
    # sparser than the transcriptome as a whole.
    min_patients = int(MIN_EXPRESSED_FRACTION * expr_lncrna.shape[1])
    keep = (expr_lncrna > 0).sum(axis=1) >= min_patients
    expr_lncrna = expr_lncrna[keep]

    print(f"  lncRNAs retained:            {expr_lncrna.shape[0]}")
    return expr_lncrna


def main():
    os.makedirs(DERIVED_DIR, exist_ok=True)

    expr, clinical, platinum = load_raw()

    clinical = filter_to_hgsoc(clinical)
    clinical = merge_platinum_status(clinical, platinum)
    clinical = add_survival_variables(clinical)

    expr_log, matched = normalize_and_filter(expr, clinical.index.tolist())
    clinical = clinical.loc[matched]

    lncrna_ids = load_lncrna_ids(GTF_PATH)
    expr_lncrna = extract_lncrnas(expr_log, lncrna_ids)

    print("\nWriting derived files")
    clinical.to_csv(os.path.join(DERIVED_DIR, "clinical_hgsoc_filtered.csv"))
    expr_log.to_csv(os.path.join(DERIVED_DIR, "expr_hgsoc_filtered.csv"))
    expr_lncrna.to_csv(os.path.join(DERIVED_DIR, "expr_lncrna_final.csv"))

    print(f"\nFinal cohort: {len(clinical)} patients")
    tiers = clinical["resistance_tier"].value_counts()
    for tier in ["Resistant", "Partially_Sensitive", "Sensitive"]:
        print(f"  {tier:<20} {tiers.get(tier, 0)}")
    print(f"Deaths observed: {clinical['os_event'].sum()}")
    print(f"lncRNAs x patients: {expr_lncrna.shape[0]} x {expr_lncrna.shape[1]}")


if __name__ == "__main__":
    main()
