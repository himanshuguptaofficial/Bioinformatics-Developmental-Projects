"""
Optimism-Corrected Validation of the Signature Selection Pipeline
=================================================================
Biological question: how much of the 5-lncRNA signature's reported
performance survives when patients are never scored by a model that has
seen their outcome?

Method: nested cross-validation that re-runs the entire selection pipeline
        (DE screen -> univariate Cox -> LASSO-Cox -> bootstrap ranking)
        inside every training fold, plus a permutation null for the whole
        procedure and a selection-stability audit.
Output: three-panel figure comparing apparent vs out-of-fold C-index,
        the permutation null, and out-of-fold Kaplan-Meier curves.

Run prepare_data.py first. Note the permutation test re-runs the full
pipeline 200 times; expect a runtime of a few hours, with progress printed
every 25 permutations.

Author: Himanshu Gupta, UC San Diego
"""

import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test
from scipy import stats
from sklearn.model_selection import StratifiedKFold
from sksurv.linear_model import CoxnetSurvivalAnalysis
from sksurv.metrics import concordance_index_censored

DERIVED_DIR = os.path.join("..", "data", "derived")
FIGURE_DIR = "figures"
RESULTS_DIR = "results"

# Selection pipeline settings, matched to the run that produced the
# signature used in projects 02 and 03.
DE_P = 0.05             # Mann-Whitney cutoff, resistant vs sensitive
COX_P = 0.05            # univariate Cox cutoff
N_SIGNATURE = 5         # size of the final signature
N_BOOTSTRAP = 40        # resamples for selection frequency inside each fold
N_SPLITS = 5
N_REPEATS = 10
N_PERM = 200
SEED = 42

# The signature as published in project 02. Selected on all 216 patients,
# which is exactly the problem this project quantifies.
PUBLISHED_SIGNATURE = [
    "ENSG00000251665", "ENSG00000258572", "ENSG00000241912",
    "ENSG00000260505", "ENSG00000287733",
]

os.makedirs(FIGURE_DIR, exist_ok=True)
os.makedirs(RESULTS_DIR, exist_ok=True)

print("Optimism-corrected validation of the selection pipeline")
print("=" * 72)

# ---------------------------------------------------------------- inputs ---

val_clin = pd.read_csv(
    os.path.join(DERIVED_DIR, "clinical_hgsoc_filtered.csv"), index_col=0
)
val_expr = pd.read_csv(
    os.path.join(DERIVED_DIR, "expr_lncrna_final.csv"), index_col=0
)
val_expr = val_expr[val_clin.index]

val_time = val_clin["os_time"].to_numpy(float)
val_event = val_clin["os_event"].to_numpy(int)
val_age = val_clin["age_years"].to_numpy(float)
val_tier = val_clin["resistance_tier"].to_numpy()
val_genes = val_expr.index.to_numpy()
val_X = val_expr.to_numpy(float)                       # genes x patients

# Age is missing for a few patients; impute at the cohort median so every
# patient stays in the fold structure. Age is a covariate, not the signal.
val_age = np.where(np.isnan(val_age), np.nanmedian(val_age), val_age)

print(f"Patients: {len(val_clin)}, lncRNAs: {val_X.shape[0]}, "
      f"deaths: {val_event.sum()}")

val_base = pd.Index([g.split(".")[0] for g in val_genes])
published_idx = np.array([int(np.where(val_base == g)[0][0])
                          for g in PUBLISHED_SIGNATURE])


def cox_score_pvalues(X, time, event):
    """Univariate Cox score test for every row of X, vectorized.

    The score test statistic for a single covariate needs only risk-set sums,
    so ranking once and using suffix cumulative sums replaces one CoxPHFitter
    call per gene. Fitting 550 genes individually takes minutes per fold and
    makes the permutation test impossible; this takes milliseconds.
    """
    order = np.argsort(time, kind="mergesort")
    t = time[order]
    e = event[order].astype(bool)
    Xo = X[:, order]

    n = len(t)
    # Risk set for an event at time t_i is every patient with time >= t_i.
    # searchsorted gives the first index of each tied time block.
    start = np.searchsorted(t, t, side="left")

    suf1 = np.zeros((X.shape[0], n + 1))
    suf2 = np.zeros((X.shape[0], n + 1))
    suf1[:, :n] = np.cumsum(Xo[:, ::-1], axis=1)[:, ::-1]
    suf2[:, :n] = np.cumsum((Xo ** 2)[:, ::-1], axis=1)[:, ::-1]
    size = (n - np.arange(n + 1)).astype(float)
    size[n] = np.nan

    k = start[e]
    s0 = size[k]
    mean = suf1[:, k] / s0
    var = suf2[:, k] / s0 - mean ** 2

    U = (Xo[:, e] - mean).sum(axis=1)
    V = var.sum(axis=1)
    V = np.where(V <= 0, np.nan, V)

    chi2 = U ** 2 / V
    return stats.chi2.sf(chi2, df=1)


def residualize_on_age(X, age):
    """Remove the linear age trend from each gene, the vectorized stand-in for
    including age as a covariate in the univariate Cox screen."""
    a = np.column_stack([np.ones_like(age), age])
    beta, *_ = np.linalg.lstsq(a, X.T, rcond=None)
    return X - (a @ beta).T


def select_signature(train, time, event, rng):
    """Run the full selection pipeline on training patients only.

    DE screen -> univariate Cox -> intersection -> LASSO-Cox -> bootstrap
    selection frequency -> top N. Returns column indices into val_X.
    """
    Xtr = val_X[:, train]
    ttr, etr = time[train], event[train]

    # 1. Differential expression, resistant vs sensitive, training patients only
    tier_tr = val_tier[train]
    res = tier_tr == "Resistant"
    sen = tier_tr == "Sensitive"
    if res.sum() < 5 or sen.sum() < 5:
        return None
    _, de_p = stats.mannwhitneyu(Xtr[:, res], Xtr[:, sen],
                                 alternative="two-sided", axis=1)
    de_hits = np.where(de_p < DE_P)[0]
    if len(de_hits) < N_SIGNATURE:
        return None

    # 2. Univariate Cox on the DE hits, age-adjusted
    resid = residualize_on_age(Xtr[de_hits], val_age[train])
    cox_p = cox_score_pvalues(resid, ttr, etr)
    keep = de_hits[np.nan_to_num(cox_p, nan=1.0) < COX_P]
    if len(keep) < N_SIGNATURE:
        return None

    # 3. LASSO-Cox with bootstrap selection frequency
    Z = Xtr[keep].T
    Z = (Z - Z.mean(axis=0)) / (Z.std(axis=0) + 1e-9)
    y = np.array(list(zip(etr.astype(bool), ttr)),
                 dtype=[("event", bool), ("time", float)])

    counts = np.zeros(len(keep))
    n_tr = len(train)
    for _ in range(N_BOOTSTRAP):
        b = rng.integers(0, n_tr, n_tr)
        if y["event"][b].sum() < 5:
            continue
        try:
            model = CoxnetSurvivalAnalysis(l1_ratio=1.0, alpha_min_ratio=0.05,
                                           n_alphas=30, max_iter=2000)
            model.fit(Z[b], y[b])
            # Take the densest solution that still keeps the model sparse.
            coefs = model.coef_
            nz = (coefs != 0).sum(axis=0)
            col = np.where(nz > 0, nz, 99)
            pick = int(np.argmin(np.abs(col - 10)))
            counts += (coefs[:, pick] != 0).astype(float)
        except Exception:
            continue

    if counts.sum() == 0:
        return None
    return keep[np.argsort(-counts)[:N_SIGNATURE]]


def fold_risk_scores(sig_idx, train, test, time, event):
    """Fit a multivariate Cox model on training patients, score the held-out."""
    Z = val_X[sig_idx].T
    mu, sd = Z[train].mean(axis=0), Z[train].std(axis=0) + 1e-9
    Z = (Z - mu) / sd

    cols = [f"g{i}" for i in range(len(sig_idx))]
    df_tr = pd.DataFrame(Z[train], columns=cols)
    df_tr["age"] = val_age[train]
    df_tr["os_time"] = time[train]
    df_tr["os_event"] = event[train]

    try:
        cph = CoxPHFitter(penalizer=0.1)
        cph.fit(df_tr, duration_col="os_time", event_col="os_event")
    except Exception:
        return None

    df_te = pd.DataFrame(Z[test], columns=cols)
    df_te["age"] = val_age[test]
    return cph.predict_partial_hazard(df_te).to_numpy()


def cross_validated_scores(time, event, seed, frozen=None, repeats=1,
                           record=None):
    """Out-of-fold risk scores. frozen=None re-runs selection inside each fold.

    record, if given, is a list that collects the signature chosen in each
    fold so selection stability can be inspected afterwards.
    """
    oof = np.full((repeats, len(time)), np.nan)
    for r in range(repeats):
        rng = np.random.default_rng(seed + r)
        skf = StratifiedKFold(n_splits=N_SPLITS, shuffle=True,
                              random_state=seed + r)
        for train, test in skf.split(np.zeros(len(time)), event):
            sig = frozen if frozen is not None else select_signature(
                train, time, event, rng)
            if sig is None:
                continue
            if record is not None:
                record.append(sig)
            pred = fold_risk_scores(sig, train, test, time, event)
            if pred is not None:
                oof[r, test] = pred
    return oof


def c_index(scores, time, event):
    ok = ~np.isnan(scores)
    if ok.sum() < 20 or event[ok].sum() < 5:
        return np.nan
    return concordance_index_censored(
        event[ok].astype(bool), time[ok], scores[ok])[0]


# ------------------------------------------------- apparent (optimistic) ---

Zp = val_X[published_idx].T
Zp = (Zp - Zp.mean(axis=0)) / Zp.std(axis=0)
df_app = pd.DataFrame(Zp, columns=[f"g{i}" for i in range(N_SIGNATURE)])
df_app["age"] = val_age
df_app["os_time"] = val_time
df_app["os_event"] = val_event
cph_app = CoxPHFitter().fit(df_app, duration_col="os_time", event_col="os_event")
apparent = c_index(cph_app.predict_partial_hazard(df_app).to_numpy(),
                   val_time, val_event)
print(f"\nApparent C-index (what project 02 reports): {apparent:.4f}")

# ------------------------------------------- honest out-of-fold estimates ---

print(f"\nNested CV, full pipeline re-selected in every fold "
      f"({N_REPEATS} repeats x {N_SPLITS} folds)...")
fold_signatures = []
oof_full = cross_validated_scores(val_time, val_event, SEED,
                                  repeats=N_REPEATS, record=fold_signatures)
c_full = np.array([c_index(oof_full[r], val_time, val_event)
                   for r in range(N_REPEATS)])
print(f"  Folds that produced a signature: {len(fold_signatures)}"
      f"/{N_REPEATS * N_SPLITS}")
print(f"  Patients scored out-of-fold: "
      f"{int((~np.isnan(oof_full)).sum(axis=1).mean())}/{len(val_time)}")
print(f"  Out-of-fold C-index: {np.nanmean(c_full):.4f} "
      f"(SD {np.nanstd(c_full):.4f})")
print(f"  Optimism: {apparent - np.nanmean(c_full):.4f}")

# Selection stability. If the pipeline were recovering real biology, the same
# lncRNAs should keep coming back when it is re-run on 80% of the cohort.
if fold_signatures:
    picked = pd.Series(np.concatenate(fold_signatures)).value_counts()
    n_folds = len(fold_signatures)
    overlap = [len(set(s) & set(published_idx)) for s in fold_signatures]
    print(f"  Distinct lncRNAs ever selected: {len(picked)} "
          f"across {n_folds} folds ({N_SIGNATURE} per fold)")
    print(f"  Most-selected lncRNA appears in {picked.iloc[0]}/{n_folds} folds")
    print(f"  Mean overlap with the published 5: {np.mean(overlap):.2f}/5")
    stability = pd.DataFrame({
        "gene_id": val_genes[picked.index.to_numpy()],
        "folds_selected": picked.to_numpy(),
        "selection_frequency": picked.to_numpy() / n_folds,
        "in_published_signature": [
            g.split(".")[0] in PUBLISHED_SIGNATURE
            for g in val_genes[picked.index.to_numpy()]],
    })
    stability.to_csv(
        os.path.join(RESULTS_DIR, "validation_selection_stability.csv"),
        index=False)

print(f"\nNested CV, signature frozen and only coefficients refit per fold...")
print("  NOTE: this arm is still contaminated. The five lncRNAs were chosen")
print("  using all 216 patients, so holding them fixed and refitting only the")
print("  coefficients leaves the selection bias untouched. It is reported to")
print("  show how much of the apparent signal survives refitting alone, not")
print("  as an unbiased estimate.")
oof_frozen = cross_validated_scores(val_time, val_event, SEED,
                                    frozen=published_idx, repeats=N_REPEATS)
c_frozen = np.array([c_index(oof_frozen[r], val_time, val_event)
                     for r in range(N_REPEATS)])
print(f"  Out-of-fold C-index: {np.nanmean(c_frozen):.4f} "
      f"(SD {np.nanstd(c_frozen):.4f})")

# --------------------------------------------------------- permutation ----
# Shuffle the (time, event) pairs and re-run the whole pipeline. This is the
# only way to know whether an out-of-fold C-index above 0.5 is meaningful for
# a procedure that screens 12k features.

print(f"\nPermutation null, {N_PERM} permutations of the full pipeline...")
perm_full = np.full(N_PERM, np.nan)
perm_frozen = np.full(N_PERM, np.nan)
rng_perm = np.random.default_rng(SEED)

for p in range(N_PERM):
    shuffle = rng_perm.permutation(len(val_time))
    pt, pe = val_time[shuffle], val_event[shuffle]
    perm_full[p] = c_index(
        cross_validated_scores(pt, pe, SEED + 1000 + p)[0], pt, pe)
    perm_frozen[p] = c_index(
        cross_validated_scores(pt, pe, SEED + 1000 + p,
                               frozen=published_idx)[0], pt, pe)
    if (p + 1) % 25 == 0:
        print(f"  {p + 1}/{N_PERM} done")

obs_full = np.nanmean(c_full)
obs_frozen = np.nanmean(c_frozen)
p_full = (1 + np.nansum(perm_full >= obs_full)) / (1 + np.sum(~np.isnan(perm_full)))
p_frozen = (1 + np.nansum(perm_frozen >= obs_frozen)) / (1 + np.sum(~np.isnan(perm_frozen)))

print(f"\n  Full pipeline:     C={obs_full:.4f}, "
      f"null mean {np.nanmean(perm_full):.4f}, permutation p={p_full:.4f}")
print(f"  Frozen signature:  C={obs_frozen:.4f}, "
      f"null mean {np.nanmean(perm_frozen):.4f}, permutation p={p_frozen:.4f}")

# ------------------------------------------- out-of-fold survival split ----
# The KM curve project 02 would show if patients were only ever scored by
# models that had never seen them.

oof_mean = np.nanmean(oof_full, axis=0)
valid = ~np.isnan(oof_mean)
groups = np.where(oof_mean[valid] > np.median(oof_mean[valid]), "High", "Low")
hi, lo = groups == "High", groups == "Low"
lr = logrank_test(val_time[valid][hi], val_time[valid][lo],
                  event_observed_A=val_event[valid][hi],
                  event_observed_B=val_event[valid][lo])
print(f"\nOut-of-fold KM split: log-rank p = {lr.p_value:.3e} "
      f"(apparent was 1.2e-09)")

# ------------------------------------------------------------- outputs ----

summary = pd.DataFrame([
    {"analysis": "Apparent (resubstitution)", "c_index": apparent,
     "sd": np.nan, "permutation_p": np.nan},
    {"analysis": "Out-of-fold, full pipeline", "c_index": obs_full,
     "sd": np.nanstd(c_full), "permutation_p": p_full},
    {"analysis": "Out-of-fold, frozen signature (still selection-biased)",
     "c_index": obs_frozen, "sd": np.nanstd(c_frozen),
     "permutation_p": p_frozen},
])
summary.to_csv(os.path.join(RESULTS_DIR, "validation_summary.csv"), index=False)
print("\n" + summary.to_string(index=False))

fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))

axes[0].bar(["Apparent", "Out-of-fold\n(full pipeline)", "Out-of-fold\n(frozen sig.)"],
            [apparent, obs_full, obs_frozen],
            yerr=[0, np.nanstd(c_full), np.nanstd(c_frozen)],
            color=["lightcoral", "steelblue", "darkseagreen"],
            edgecolor="black", capsize=5)
axes[0].axhline(0.5, color="black", linestyle="--", alpha=0.7,
                label="Chance (0.5)")
axes[0].set_ylim(0.4, 0.8)
axes[0].set_ylabel("Harrell C-index", fontsize=11)
axes[0].set_title("Optimism in the reported C-index", fontsize=12)
axes[0].legend(fontsize=9)
axes[0].grid(axis="y", alpha=0.2)

axes[1].hist(perm_full[~np.isnan(perm_full)], bins=30, color="lightgrey",
             edgecolor="white", label="Permuted outcomes")
axes[1].axvline(obs_full, color="firebrick", lw=2,
                label=f"Observed ({obs_full:.3f})")
axes[1].axvline(0.5, color="black", linestyle="--", alpha=0.6)
axes[1].set_xlabel("Out-of-fold C-index", fontsize=10)
axes[1].set_ylabel("Permutations", fontsize=10)
axes[1].set_title(f"Permutation null, full pipeline\np = {p_full:.4f}",
                  fontsize=12)
axes[1].legend(fontsize=9)
axes[1].grid(alpha=0.2)

for mask, label, color in [(hi, "High risk", "firebrick"),
                           (lo, "Low risk", "steelblue")]:
    kmf = KaplanMeierFitter()
    kmf.fit(val_time[valid][mask], event_observed=val_event[valid][mask],
            label=f"{label} (n={mask.sum()})")
    kmf.plot_survival_function(ax=axes[2], color=color, ci_show=True, linewidth=2)
axes[2].set_title(f"Out-of-fold risk groups\nlog-rank p = {lr.p_value:.2e}",
                  fontsize=12)
axes[2].set_xlabel("Time (days)", fontsize=10)
axes[2].set_ylabel("Overall survival probability", fontsize=10)
axes[2].set_ylim(0, 1.02)
axes[2].grid(alpha=0.2)

fig.suptitle("Internal validation of the 5-lncRNA signature, TCGA-OV",
             fontsize=13)
fig.tight_layout()
plt.savefig(os.path.join(FIGURE_DIR, "validation_optimism.png"), dpi=300,
            bbox_inches="tight")
plt.close()

print("\nSaved figures/validation_optimism.png, "
      "results/validation_summary.csv and "
      "results/validation_selection_stability.csv")
