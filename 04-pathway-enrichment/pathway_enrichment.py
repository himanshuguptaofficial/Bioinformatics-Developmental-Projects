"""
Pathway Enrichment of Genes Co-expressed with the lncRNA Signature
=================================================================
Biological question: lncRNAs act largely by regulating other genes, so what
biological processes are the platinum-resistance signature lncRNAs wired
into?

Method: Spearman correlation between each signature lncRNA and every
        protein-coding gene, |rho| >= 0.3 and p < 0.05, then over-representation
        analysis of the union against KEGG and GO Biological Process.
Output: ranked enrichment tables for both libraries.

Requires network access for the Enrichr API. Run prepare_data.py first.

Author: Himanshu Gupta, UC San Diego
"""

import gzip
import os

import gseapy as gp
import numpy as np
import pandas as pd
from scipy import stats

DERIVED_DIR = os.path.join("..", "data", "derived")
GTF_PATH = os.path.join("..", "data", "gencode.v44.gtf.gz")
RESULT_DIR = "results"

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


def spearman_against_all(lncrna_expr, gene_expr):
    """Spearman rho and p-value for every (lncRNA, gene) pair.

    Spearman is Pearson on ranks, so ranking once and taking a single matrix
    product is equivalent to looping over ~19k genes five times and far
    faster.
    """
    n_samples = lncrna_expr.shape[1]

    def zscored_ranks(matrix):
        ranks = stats.rankdata(matrix, axis=1)
        centered = ranks - ranks.mean(axis=1, keepdims=True)
        return centered / np.linalg.norm(centered, axis=1, keepdims=True)

    lnc_z = zscored_ranks(lncrna_expr.to_numpy())
    gene_z = zscored_ranks(gene_expr.to_numpy())

    rho = lnc_z @ gene_z.T
    rho = np.clip(rho, -0.9999999, 0.9999999)

    # Two-sided p from the usual t approximation.
    t_stat = rho * np.sqrt((n_samples - 2) / (1 - rho**2))
    pvalues = 2 * stats.t.sf(np.abs(t_stat), df=n_samples - 2)

    return rho, pvalues


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

    lncrna_expr = expr_lncrna.loc[signature, coding_expr.columns]
    rho, pvalues = spearman_against_all(lncrna_expr, coding_expr)

    significant = (np.abs(rho) >= CORRELATION_THRESHOLD) & (pvalues < P_THRESHOLD)

    for gene, row in zip(signature, significant):
        print(f"  {gene}: {int(row.sum())} co-expressed genes")

    union = coding_base[significant.any(axis=0)]
    print(f"\nUnion of co-expressed genes: {len(union)}")

    return list(union)


def run_enrichment(gene_symbols, library, label):
    """Over-representation analysis against one Enrichr library."""
    print(f"\nRunning {label} enrichment on {len(gene_symbols)} symbols")

    enrichment = gp.enrichr(
        gene_list=gene_symbols,
        gene_sets=[library],
        organism="human",
        outdir=None,
    )
    results = enrichment.results.sort_values("Adjusted P-value")

    n_significant = int((results["Adjusted P-value"] < 0.05).sum())
    print(f"  Significant terms (adj p < 0.05): {n_significant}")
    for _, row in results.head(5).iterrows():
        print(f"    {row['Term'][:60]:<62} adj p={row['Adjusted P-value']:.2e}")

    return results


def main():
    os.makedirs(RESULT_DIR, exist_ok=True)

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

    kegg = run_enrichment(gene_symbols, "KEGG_2021_Human", "KEGG")
    go = run_enrichment(gene_symbols, "GO_Biological_Process_2021", "GO BP")

    kegg.to_csv(os.path.join(RESULT_DIR, "kegg_enrichment.csv"), index=False)
    go.to_csv(os.path.join(RESULT_DIR, "go_enrichment.csv"), index=False)
    print(f"\nWrote enrichment tables to {RESULT_DIR}/")


if __name__ == "__main__":
    main()
