# lncRNAs and Platinum Resistance in High-Grade Serous Ovarian Cancer

Five self-contained analyses of long non-coding RNA expression in the TCGA-OV
cohort. Each one answers a single question and produces a single figure. They
share a dataset and build on each other, but every project runs independently.

Together they tell a story with a twist: projects 01–04 build a five-lncRNA
survival signature the way the biomarker literature usually does, complete
with a log-rank p of 1.2 × 10⁻⁹ — and project 05 shows that essentially all
of that performance is manufactured by the selection procedure itself. The
most useful output of this repository is not a signature. It is a quantified,
reproducible demonstration of why signatures built this way so rarely
replicate.

## Motivation
High-grade serous ovarian cancer (HGSOC) causes 70–80% of ovarian cancer
deaths. Most patients respond to first-line platinum chemotherapy, but 80–90%
relapse with resistant disease. Resistance is graded by the platinum-free
interval (PFI): under 6 months is resistant, over 12 months is sensitive, per
the GCIG consensus. These analyses ask whether lncRNAs, which are still largely
uncharacterised, carry information about that outcome — and then ask the
harder question of whether the apparent answer survives honest validation.

## The five projects

| # | Project | Question | Result |
|---|---|---|---|
| 01 | [Differential expression](01-differential-expression) | Which lncRNAs differ between resistant and sensitive patients? | 554 at raw p < 0.05, none survive FDR correction |
| 02 | [Survival risk score](02-survival-risk-score) | Does a composite of five lncRNAs stratify survival? | Apparent log-rank p = 1.2 × 10⁻⁹, C-index 0.686 — see 05 |
| 03 | [Time-dependent ROC](03-time-dependent-roc) | Does it beat FIGO stage and age? | Apparent 5-year AUC 0.748 vs 0.480 (age), 0.471 (stage) — see 05 |
| 04 | [Pathway enrichment](04-pathway-enrichment) | What biology is the signature's co-expression network in? | Overwhelmingly immune, cytokine and T-cell signaling |
| 05 | [Honest validation](05-honest-validation) | How much of 02–03 survives when selection is kept out of evaluation? | Out-of-fold C-index 0.517, permutation p = 0.33: chance |

Read in order, the arc is: no single lncRNA predicts resistance (01); a
composite of five appears to predict survival spectacularly (02) and to beat
the clinical variables oncologists already use (03); its co-expression
neighbourhood looks like tumour immunology (04); and when the entire selection
pipeline is re-run inside cross-validation folds so no patient is ever scored
by a model that saw their outcome, the prognostic performance evaporates
(05). The signature also barely re-selects itself: across 50 training folds
the pipeline picked 76 different lncRNAs, and the most stable one returned in
only 33 of 50 folds. Project 01's null was the truthful result, and 05
explains exactly how 02 and 03 turned that null into p = 10⁻⁹.

## Repository layout
```
prepare_data.py               shared preprocessing, run once
DownloadDatafromR             raw TCGA-OV download via TCGAbiolinks
data/                         all data, untracked (see below)
01-differential-expression/   script + figures/ + README
02-survival-risk-score/
03-time-dependent-roc/
04-pathway-enrichment/
05-honest-validation/
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
            lifelines scikit-learn scikit-survival gseapy openpyxl
```

## Scope and limitations
These are exploratory analyses on a single public cohort, and the single
biggest limitation is quantified rather than hand-waved: the five-lncRNA
signature was selected on the same TCGA-OV patients it is evaluated on, and
project 05 shows that this alone accounts for essentially all of the apparent
performance in projects 02 and 03 (out-of-fold C-index 0.517, permutation
p = 0.33). No suitable external cohort with both lncRNA coverage and outcomes
currently exists, so internal validation is as far as this data can go. The
immune signal in project 04 comes from bulk RNA-seq and cannot separate
tumour-intrinsic expression from immune infiltration. Each project's README
states its own limitations in full.
