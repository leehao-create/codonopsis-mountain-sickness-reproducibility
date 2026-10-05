#!/usr/bin/env python3
"""Build the frozen ML feature-prioritization figure without refitting models."""

from pathlib import Path
import re

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec


MODULE = Path(__file__).resolve().parents[1]
FINAL = MODULE / "results"
PHASE5 = MODULE / "data" / "dependencies"
OUT = MODULE / "reproduced_figures"
DATA = OUT / "plot_data"
OUT.mkdir(parents=True, exist_ok=True)
DATA.mkdir(parents=True, exist_ok=True)

COLORS = {
    "LASSO": "#0072B2",
    "Elastic Net": "#009E73",
    "RF": "#D55E00",
    "SVM-RFE": "#CC79A7",
    "Tier 1": "#A71930",
    "Tier 2": "#2A7F9E",
    "Secondary": "#AEB6BF",
}
THRESHOLD = 0.60

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7.4,
    "axes.titlesize": 8.3,
    "axes.labelsize": 7.5,
    "xtick.labelsize": 6.7,
    "ytick.labelsize": 6.7,
    "legend.fontsize": 6.7,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.55,
    "ytick.major.width": 0.55,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "savefig.facecolor": "white",
})


def clean_axis(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(False)
    ax.tick_params(length=2.5, pad=2)


def panel_label(ax, label):
    ax.text(-0.17, 1.08, label, transform=ax.transAxes, fontsize=11,
            fontweight="bold", ha="left", va="top")


def top_frequency(stability, column, n=10):
    d = stability[["gene", column]].dropna().sort_values(column, ascending=False).head(n)
    return d.sort_values(column)


def frequency_bar(ax, data, column, color, title):
    y = np.arange(len(data))
    ax.barh(y, data[column], color=color, alpha=0.88, height=0.68)
    ax.axvline(THRESHOLD, color="#777777", lw=0.7, ls=(0, (3, 2)))
    ax.set_yticks(y, data["gene"])
    ax.set_xlim(0, 1.04)
    ax.set_xlabel("Selection frequency")
    ax.set_title(title, loc="left", fontweight="bold", pad=4)
    clean_axis(ax)


stability = pd.read_csv(FINAL / "GSE103940/GSE103940_feature_stability_ranked.csv")
tuning = pd.read_csv(FINAL / "GSE103940/GSE103940_outer_fold_tuning.csv")
rf_long = pd.read_csv(PHASE5 / "ML_RF_importance_stability_long.csv")
rf_long = rf_long.loc[rf_long["dataset"].eq("GSE103940")].copy()
pred = pd.read_csv(FINAL / "GSE103940/GSE103940_final_LOSO_OOF_predictions.csv")
perf = pd.read_csv(FINAL / "GSE103940/GSE103940_final_OOF_performance.csv")

# Parse frozen fold-level tuning summaries. These are visual summaries, not refits.
lasso_tune = tuning.loc[tuning["method"].eq("LASSO")].copy()
lasso_tune["lambda_1se"] = lasso_tune["value"].str.split("/").str[1].astype(float)
lasso_tune["selected_n"] = lasso_tune["secondary"].str.extract(r"selected_n=(\d+)").astype(int)

en_tune = tuning.loc[tuning["method"].eq("Elastic Net")].copy()
en_tune["alpha"] = en_tune["value"].str.split("/").str[0].astype(float)
en_tune["lambda_1se"] = en_tune["value"].str.split("/").str[1].astype(float)
en_tune["selected_n"] = en_tune["secondary"].str.extract(r"selected_n=(\d+)").astype(int)

svm_tune = tuning.loc[tuning["method"].eq("SVM-RFE")].copy()
svm_tune["cost"] = svm_tune["value"].astype(float)
svm_tune["subset_size"] = svm_tune["secondary"].str.extract(r"subset_size=(\d+)").astype(int)
svm_tune["inner_balanced_accuracy"] = svm_tune["secondary"].str.extract(
    r"inner_balanced_accuracy=([0-9.]+)").astype(float)

rf_summary = (rf_long.groupby("gene", as_index=False)
              .agg(median_permutation_importance=("observed_importance", "median"),
                   RF_frequency=("selected", "mean")))
rf_summary = rf_summary.sort_values(["RF_frequency", "median_permutation_importance"],
                                    ascending=False).head(10)
rf_summary = rf_summary.sort_values("median_permutation_importance")

# Stable membership for the UpSet-like panel.
method_cols = {
    "LASSO": "LASSO_frequency",
    "Elastic Net": "ElasticNet_frequency",
    "RF": "RF_frequency",
    "SVM-RFE": "SVM_frequency",
}
eligible = stability.dropna(subset=list(method_cols.values())).copy()
for method, col in method_cols.items():
    eligible[method] = eligible[col] >= THRESHOLD
eligible["pattern"] = eligible[list(method_cols)].apply(
    lambda row: "|".join([m for m in method_cols if row[m]]), axis=1)
upset = (eligible.loc[eligible["pattern"].ne("")]
         .groupby("pattern", as_index=False)
         .agg(count=("gene", "size"), genes=("gene", lambda x: ";".join(x))))
upset = upset.sort_values(["count", "pattern"], ascending=[False, True]).reset_index(drop=True)

tiers = stability.loc[stability["number_of_methods_supported"].ge(3),
                      ["gene", "number_of_methods_supported"]].copy()
tiers["tier"] = np.where(tiers["number_of_methods_supported"].eq(4), "Tier 1", "Tier 2")
mean_lookup = stability.set_index("gene")[list(method_cols.values())].mean(axis=1)
tiers["mean_frequency"] = tiers["gene"].map(mean_lookup)
tiers = tiers.sort_values(["number_of_methods_supported", "mean_frequency"], ascending=[True, True])

# Save exact displayed data.
lasso_tune.to_csv(DATA / "panel_A_LASSO_nested_tuning.csv", index=False)
top_frequency(stability, "LASSO_frequency").to_csv(DATA / "panel_B_LASSO_stability.csv", index=False)
rf_summary.to_csv(DATA / "panel_C_RF_importance_stability.csv", index=False)
svm_tune.to_csv(DATA / "panel_D_SVM_RFE_nested_tuning.csv", index=False)
top_frequency(stability, "SVM_frequency").to_csv(DATA / "panel_E_SVM_stability.csv", index=False)
top_frequency(stability, "ElasticNet_frequency", 12).to_csv(DATA / "panel_F_ElasticNet_stability.csv", index=False)
upset.to_csv(DATA / "panel_G_stable_set_intersections.csv", index=False)
tiers.to_csv(DATA / "panel_H_final_tiers.csv", index=False)

# Main figure: 180 mm wide, feature prioritization only.
fig = plt.figure(figsize=(7.087, 10.1), constrained_layout=False)
gs = GridSpec(4, 2, figure=fig, left=0.105, right=0.965, top=0.965, bottom=0.075,
              hspace=0.72, wspace=0.64, height_ratios=[1.0, 1.04, 1.08, 1.15])

# A: fold-specific nested tuning instead of a leakage-prone full-data CV curve.
ax = fig.add_subplot(gs[0, 0])
x = np.arange(1, len(lasso_tune) + 1)
ax.plot(x, lasso_tune["lambda_1se"], color=COLORS["LASSO"], marker="o", ms=3, lw=1)
ax.set_xlabel("Outer LOSO fold")
ax.set_ylabel(r"Selected $\lambda$ (1se)")
ax.set_xticks([1, 3, 5, 7, 9, 11])
ax2 = ax.twinx()
ax2.plot(x, lasso_tune["selected_n"], color="#6F7A83", marker="s", ms=2.8, lw=0.8,
         ls=(0, (3, 2)))
ax2.set_ylabel("Number of selected genes", color="#59636B")
ax2.spines["top"].set_visible(False)
ax2.spines["right"].set_linewidth(0.6)
ax.set_title("LASSO nested tuning and feature selection", loc="left", fontweight="bold", pad=6)
clean_axis(ax)
panel_label(ax, "A")

# B.
ax = fig.add_subplot(gs[0, 1])
frequency_bar(ax, top_frequency(stability, "LASSO_frequency"), "LASSO_frequency",
              COLORS["LASSO"], "LASSO selection stability")
panel_label(ax, "B")

# C.
ax = fig.add_subplot(gs[1, 0])
y = np.arange(len(rf_summary))
bars = ax.barh(y, rf_summary["median_permutation_importance"], height=0.68,
               color=plt.cm.Oranges(0.25 + 0.65 * rf_summary["RF_frequency"].to_numpy()))
ax.set_yticks(y, rf_summary["gene"])
ax.set_xlabel("Median permutation importance")
ax.set_title("Random Forest permutation importance", loc="left", fontweight="bold", pad=6)
for b, f in zip(bars, rf_summary["RF_frequency"]):
    ax.text(b.get_width(), b.get_y() + b.get_height()/2, f"  {f:.2f}", va="center", fontsize=6)
clean_axis(ax)
ax.xaxis.set_major_locator(mpl.ticker.MaxNLocator(5))
panel_label(ax, "C")

# D.
ax = fig.add_subplot(gs[1, 1])
jitter = np.linspace(-0.012, 0.012, len(svm_tune))
ax.scatter(svm_tune["subset_size"], svm_tune["inner_balanced_accuracy"] + jitter,
           s=27, color=COLORS["SVM-RFE"], edgecolor="white", linewidth=0.45, zorder=3)
for size, group in svm_tune.groupby("subset_size"):
    ax.vlines(size, group["inner_balanced_accuracy"].min(), group["inner_balanced_accuracy"].max(),
              color="#B8B8B8", lw=0.7, zorder=1)
ax.set_xticks(sorted(svm_tune["subset_size"].unique()))
ax.set_ylim(0.91, 1.015)
ax.set_xlabel("Selected feature subset size")
ax.set_ylabel("Inner balanced accuracy")
ax.set_title("SVM-RFE nested optimization", loc="left", fontweight="bold", pad=6)
clean_axis(ax)
panel_label(ax, "D")

# E.
ax = fig.add_subplot(gs[2, 0])
frequency_bar(ax, top_frequency(stability, "SVM_frequency"), "SVM_frequency",
              COLORS["SVM-RFE"], "SVM-RFE selection stability")
panel_label(ax, "E")

# F.
ax = fig.add_subplot(gs[2, 1])
en_top = top_frequency(stability, "ElasticNet_frequency", 12)
frequency_bar(ax, en_top, "ElasticNet_frequency", COLORS["Elastic Net"],
              "Elastic Net selection stability")
panel_label(ax, "F")

# G: compact UpSet-like stable-set intersection summary.
sub = GridSpecFromSubplotSpec(2, 1, subplot_spec=gs[3, 0], height_ratios=[1.25, 1], hspace=0.16)
axb = fig.add_subplot(sub[0])
x = np.arange(len(upset))
axb.bar(x, upset["count"], color="#475B6B", width=0.68)
for xi, n in zip(x, upset["count"]):
    axb.text(xi, n + 0.08, str(n), ha="center", va="bottom", fontsize=6.2)
axb.set_ylabel("Genes")
axb.set_xticks([])
axb.set_title("Stable-set intersections (frequency ≥ 0.60)", loc="left", fontweight="bold", pad=6)
clean_axis(axb)
panel_label(axb, "G")
axm = fig.add_subplot(sub[1], sharex=axb)
methods_order = list(method_cols)
for xi, pattern in enumerate(upset["pattern"]):
    active = pattern.split("|")
    ys = [methods_order.index(m) for m in active]
    if len(ys) > 1:
        axm.plot([xi, xi], [min(ys), max(ys)], color="#3F4850", lw=1)
    for yi, method in enumerate(methods_order):
        active_flag = method in active
        axm.scatter(xi, yi, s=23 if active_flag else 13,
                    color=COLORS.get(method, "#666666") if active_flag else "#D9DEE2",
                    edgecolor="none", zorder=3)
axm.set_yticks(np.arange(len(methods_order)), methods_order)
axm.invert_yaxis()
axm.set_ylim(len(methods_order) - 0.35, -0.75)
axm.set_xticks(x)
pattern_labels = {
    "Elastic Net": "EN only",
    "Elastic Net|RF|SVM-RFE": "EN + RF\n+ SVM",
    "LASSO|Elastic Net|RF|SVM-RFE": "All four",
    "Elastic Net|SVM-RFE": "EN + SVM",
    "Elastic Net|RF": "EN + RF",
}
axm.set_xticklabels([pattern_labels.get(p, p.replace("|", " + ")) for p in upset["pattern"]],
                    rotation=38, ha="right", rotation_mode="anchor", fontsize=5.8)
axb.tick_params(axis="x", which="both", bottom=False, labelbottom=False)
for side in axm.spines.values():
    side.set_visible(False)
axm.tick_params(axis="both", length=0)

# H.
ax = fig.add_subplot(gs[3, 1])
y = np.arange(len(tiers))
ax.barh(y, tiers["number_of_methods_supported"],
        color=[COLORS[t] for t in tiers["tier"]], height=0.66)
ax.set_yticks(y, tiers["gene"])
ax.set_xlim(0, 4.85)
ax.set_xticks([0, 1, 2, 3, 4])
ax.set_xlabel("Algorithms with stable support")
ax.set_title("Consensus machine-learning-\nprioritized genes", loc="left",
             fontweight="bold", pad=6)
for yi, (count, tier) in enumerate(zip(tiers["number_of_methods_supported"], tiers["tier"])):
    ax.text(count + 0.18, yi, tier, va="center", fontsize=6.5)
clean_axis(ax)
panel_label(ax, "H")

for ext, kwargs in {
    "pdf": {},
    "png": {"dpi": 600},
    "tiff": {"dpi": 600, "pil_kwargs": {"compression": "tiff_lzw"}},
}.items():
    fig.savefig(OUT / f"Fig5_feature_prioritization.{ext}", **kwargs)
plt.close(fig)

# Supplementary figure: empirical OOF ROC only, never training predictions.
model_map = {
    "LASSO": "LASSO_probability",
    "Elastic Net": "ElasticNet_probability",
    "RF": "RF_probability",
    "SVM-RFE": "SVM_probability",
}
auc_map = {"LASSO": 0.942, "Elastic Net": 0.992, "RF": 0.983, "SVM-RFE": 1.000}

def roc_points(y, score):
    thresholds = np.r_[np.inf, np.sort(np.unique(score))[::-1], -np.inf]
    rows = []
    for threshold in thresholds:
        called = score >= threshold
        rows.append((np.mean(called[y == 0]), np.mean(called[y == 1])))
    return np.asarray(rows)

fig, ax = plt.subplots(figsize=(4.6, 4.25))
y = pred["observed"].to_numpy(dtype=int)
roc_rows = []
for model, col in model_map.items():
    xy = roc_points(y, pred[col].to_numpy(float))
    ax.step(xy[:, 0], xy[:, 1], where="post", lw=1.25, color=COLORS[model],
            label=f"{model} (AUC {auc_map[model]:.3f})")
    roc_rows.extend([{"model": model, "FPR": a, "TPR": b} for a, b in xy])
ax.plot([0, 1], [0, 1], color="#A7ADB2", lw=0.7, ls=(0, (3, 2)))
ax.set_xlim(0, 1)
ax.set_ylim(0, 1.02)
ax.set_aspect("equal", adjustable="box")
ax.set_xlabel("False-positive rate")
ax.set_ylabel("True-positive rate")
ax.set_title("Outer-LOSO out-of-fold performance", loc="left", fontweight="bold")
ax.legend(frameon=False, loc="lower right")
clean_axis(ax)
fig.tight_layout()
for ext, kwargs in {
    "pdf": {},
    "png": {"dpi": 600},
    "tiff": {"dpi": 600, "pil_kwargs": {"compression": "tiff_lzw"}},
}.items():
    fig.savefig(OUT / f"Supplementary_ML_OOF_ROC.{ext}", **kwargs)
plt.close(fig)
pd.DataFrame(roc_rows).to_csv(DATA / "supplementary_OOF_ROC_coordinates.csv", index=False)

print(f"Wrote figures and displayed data to {OUT}")
