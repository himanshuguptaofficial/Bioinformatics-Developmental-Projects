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
Applies the filtering cascade all four projects depend on and writes three
files into `data/derived/`. The final cohort is 216 patients: 67 resistant, 55
partially sensitive and 94 sensitive, with 145 deaths observed.

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
