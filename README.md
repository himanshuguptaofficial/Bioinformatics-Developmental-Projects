# lncRNAs and Platinum Resistance in High-Grade Serous Ovarian Cancer

Four self-contained analyses of long non-coding RNA expression in the TCGA-OV
cohort. Each one answers a single biological question and produces a single
figure. They share a dataset and build on each other, but every project runs
independently.

Author: Himanshu Gupta, UC San Diego

## Motivation
High-grade serous ovarian cancer (HGSOC) causes 70–80% of ovarian cancer
deaths. Most patients respond to first-line platinum chemotherapy, but 80–90%
relapse with resistant disease. Resistance is graded by the platinum-free
interval (PFI): under 6 months is resistant, over 12 months is sensitive, per
the GCIG consensus. These analyses ask whether lncRNAs, which are still largely
uncharacterised, carry information about that outcome.

## The four projects

| # | Project | Question | Result |
|---|---|---|---|
| 01 | [Differential expression](01-differential-expression) | Which lncRNAs differ between resistant and sensitive patients? | 554 at raw p < 0.05, none survive FDR correction |
| 02 | [Survival risk score](02-survival-risk-score) | Does a composite of five lncRNAs stratify survival? | log-rank p = 1.2 × 10⁻⁹, roughly a 2-year median survival gap |
| 03 | [Time-dependent ROC](03-time-dependent-roc) | Does it beat FIGO stage and age? | 5-year AUC 0.748 against 0.480 (age) and 0.471 (stage) |
| 04 | [Pathway enrichment](04-pathway-enrichment) | What biology is the signature in? | Overwhelmingly immune, cytokine and T-cell signaling |

Read in order they tell one story. No single lncRNA predicts platinum
resistance, but a composite of five predicts survival better than the clinical
variables currently in use, and the biology it reads out looks like the tumour
immune microenvironment rather than a cell-intrinsic resistance mechanism.

## Repository layout
```
prepare_data.py               shared preprocessing, run once
DownloadDatafromR             raw TCGA-OV download via TCGAbiolinks
data/                         all data, untracked (see below)
01-differential-expression/   script + figures/ + README
02-survival-risk-score/
03-time-dependent-roc/
04-pathway-enrichment/
```

Each project folder holds one script, the figure it produces, and a README
explaining the finding.

## Setup

### 1. Download the raw data
```bash
Rscript DownloadDatafromR
```
Pulls the TCGA-OV expression and clinical tables from the GDC along with the
GENCODE v44 annotation. Both are open access, so no login or API token is
needed. The first run caches roughly 1.7 GB of raw STAR-count files into
`GDCdata/` and takes a while. Anything already present is skipped.

One file has to be fetched by hand, because Nature serves it behind a cookie
wall. Download Supplementary Table S1.2 from
[nature.com/articles/nature10166](https://www.nature.com/articles/nature10166)
and save it as `data/Table_S1_2.xlsx`. It carries the platinum-free interval
that defines the resistant and sensitive groups.

### 2. Build the derived matrices
```bash
python prepare_data.py
```
Applies the filtering cascade all four projects depend on:

| Step | Patients |
|---|---|
| TCGA-OV samples | 434 |
| Primary tumours only | 426 |
| Confirmed HGSOC (serous cystadenocarcinoma, NOS) | 404 |
| At least 30 days follow-up | 394 |
| Non-missing disease response | 346 |
| Definitive platinum status | 216 |

The final cohort is 216 patients: 67 resistant, 55 partially sensitive and 94
sensitive, with 145 deaths observed. After log2(TPM+1) normalization and
low-expression filtering that leaves 12,290 lncRNAs. Three files are written
into `data/derived/`.

### 3. Run any project
```bash
cd 01-differential-expression
python differential_expression.py
```

## Data policy
No data is tracked in git. `data/` and `GDCdata/` are gitignored, since the raw
download is ~1.7 GB and the derived matrices add ~190 MB, all of it
reproducible from the two commands above. Figures are tracked so the results
are visible without running anything. The intermediate CSVs each project writes
to its own `results/` folder are regenerated on every run and are not tracked.

## Dependencies
R, for the data download only: `TCGAbiolinks` and `SummarizedExperiment`,
installed automatically by `DownloadDatafromR` through BiocManager.

Python 3:
```bash
pip install pandas numpy scipy statsmodels matplotlib \
            lifelines scikit-learn gseapy openpyxl
```

## Scope and limitations
These are exploratory analyses on a single public cohort. The five-lncRNA
signature was selected on the same TCGA-OV patients it is evaluated on, so the
C-index in project 02 and the AUCs in project 03 are optimistic, and external
validation would be needed before any claim of clinical utility. The immune
signal in project 04 comes from bulk RNA-seq and cannot separate
tumour-intrinsic expression from immune infiltration. Each project's README
states its own limitations in full.
