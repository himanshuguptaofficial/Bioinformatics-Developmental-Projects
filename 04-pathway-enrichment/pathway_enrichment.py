"""
Genes Co-expressed with the lncRNA Signature
============================================
Biological question: lncRNAs act largely by regulating other genes, so which
protein-coding genes track the platinum-resistance signature lncRNAs across
patients?

Method: Spearman correlation between each signature lncRNA and every
        protein-coding gene, keeping |rho| >= 0.3 and p < 0.05.

Run prepare_data.py first.

Author: Himanshu Gupta, UC San Diego
"""

import gzip
import os

import pandas as pd
from scipy import stats

DERIVED_DIR = os.path.join("..", "data", "derived")
GTF_PATH = os.path.join("..", "data", "gencode.v44.gtf.gz")

# Same signature as project 02; see that script for provenance.
SIGNATURE_LNCRNAS = [
    "ENSG00000251665",
    "ENSG00000258572",  # SYNE3-AS1
    "ENSG00000241912",
    "ENSG00000260505",
    "ENSG00000287733",
]

CORRELATION_THRESHOLD = 0.3
P_THRESHOLD = 0.05


def parse_gencode(gtf_path):
    """Return protein-coding Ensembl IDs and an ID-to-symbol lookup."""
    print("Parsing GENCODE annotation")

    protein_coding = set()
    symbols = {}

    with gzip.open(gtf_path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue

            gene_id, gene_name, gene_type = "", "", ""
            for attr in fields[8].split(";"):
                attr = attr.strip()
                if attr.startswith("gene_id"):
                    gene_id = attr.split('"')[1].split(".")[0]
                elif attr.startswith("gene_name"):
                    gene_name = attr.split('"')[1]
                elif attr.startswith("gene_type"):
                    gene_type = attr.split('"')[1]

            if not gene_id:
                continue
            if gene_type == "protein_coding":
                protein_coding.add(gene_id)
            if gene_name:
                symbols[gene_id] = gene_name

    print(f"  Protein-coding genes: {len(protein_coding)}")
    return protein_coding, symbols


def find_coexpressed(expr_lncrna, expr_all, protein_coding):
    """Collect protein-coding genes correlated with any signature lncRNA."""
    base_ids = expr_lncrna.index.str.split(".").str[0]
    signature = [expr_lncrna.index[base_ids == g][0] for g in SIGNATURE_LNCRNAS]

    # Restrict the search space to genuine protein-coding genes.
    all_base = expr_all.index.str.split(".").str[0]
    coding_mask = all_base.isin(protein_coding)
    coding_expr = expr_all[coding_mask]
    coding_base = all_base[coding_mask]
    print(f"  Protein-coding genes in matrix: {len(coding_expr)}")

    union = set()
    for lncrna in signature:
        lncrna_values = expr_lncrna.loc[lncrna, coding_expr.columns]

        hits = []
        for gene_id, gene_values in zip(coding_base, coding_expr.to_numpy()):
            rho, pvalue = stats.spearmanr(lncrna_values, gene_values)
            if abs(rho) >= CORRELATION_THRESHOLD and pvalue < P_THRESHOLD:
                hits.append(gene_id)

        print(f"  {lncrna}: {len(hits)} co-expressed genes")
        union.update(hits)

    print(f"\nUnion of co-expressed genes: {len(union)}")
    return sorted(union)


def main():
    protein_coding, symbols = parse_gencode(GTF_PATH)

    print("\nLoading expression matrices")
    expr_lncrna = pd.read_csv(
        os.path.join(DERIVED_DIR, "expr_lncrna_final.csv"), index_col=0
    )
    expr_all = pd.read_csv(
        os.path.join(DERIVED_DIR, "expr_hgsoc_filtered.csv"), index_col=0
    )
    print(f"  All genes: {expr_all.shape[0]}, lncRNAs: {expr_lncrna.shape[0]}")

    coexpressed = find_coexpressed(expr_lncrna, expr_all, protein_coding)

    gene_symbols = sorted({symbols[g] for g in coexpressed if g in symbols})
    print(f"Mapped to {len(gene_symbols)} gene symbols")


if __name__ == "__main__":
    main()
