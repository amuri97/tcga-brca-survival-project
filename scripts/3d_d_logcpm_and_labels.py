#!/usr/bin/env python
"""
3D-D: Normalize TCGA counts to CPM->log2(CPM+1) and align/save labels.

Inputs:
  - data_proc/tcga_counts_raw.parquet
  - data_proc/tcga_survival_join.tsv

Outputs:
  - data_proc/tcga_expr_logcpm.parquet  (genes × samples, float32)
  - data_proc/tcga_labels.tsv           (SAMPLE_ID, os_event[int8], os_time_months[float32])
"""

import pandas as pd, numpy as np
from pathlib import Path

# reload your saved counts from 3D-C
Xc = pd.read_parquet("data_proc/tcga_counts_raw.parquet")   # (59427, 1094), int32
print("counts shape:", Xc.shape, "| dtype:", Xc.dtypes.iloc[0])


# library sizes per sample (float to avoid overflow)
libs = Xc.sum(axis=0).astype(np.float64)

# guard against divide-by-zero (should be none, but be safe)
zero_libs = libs.eq(0).sum()
print("zero-size libraries:", int(zero_libs))
libs = libs.replace(0, np.nan)

# CPM = counts / libsize * 1e6
CPM = (Xc.divide(libs, axis=1) * 1e6)

# any NaNs created by zero libs? (should be 0)
print("CPM NaNs:", int(CPM.isna().sum().sum()))
# fill potential NaNs with 0 just in case
CPM = CPM.fillna(0)


X_logcpm = np.log2(CPM + 1.0).astype("float32")
print("logCPM shape:", X_logcpm.shape, "| dtype:", X_logcpm.dtypes.iloc[0])

# quick QA: no inf/NaN
nonfinite = ~np.isfinite(X_logcpm.to_numpy()).all()
print("any non-finite?", bool(nonfinite))

# optional peek
print("logCPM mean/std:", float(X_logcpm.values.mean()), float(X_logcpm.values.std()))


out_expr = Path("data_proc/tcga_expr_logcpm.parquet")
X_logcpm.to_parquet(out_expr, index=True)
print("saved:", out_expr, "| size (MB) ~", round(out_expr.stat().st_size/1e6, 1))







# load the joined table and align
joined = pd.read_csv("data_proc/tcga_survival_join.tsv", sep="\t")
joined["submitter_id"] = joined["submitter_id"].astype(str).str.upper().str.strip()

labels_tcga = (joined
    .set_index("submitter_id")
    .loc[X_logcpm.columns, ["os_event", "os_time_months"]]
    .rename_axis("SAMPLE_ID")
    .reset_index())


# tidy dtypes (optional but nice)
labels_tcga["os_event"] = labels_tcga["os_event"].astype("int8")
labels_tcga["os_time_months"] = labels_tcga["os_time_months"].astype("float32")

assert labels_tcga["SAMPLE_ID"].is_unique
assert set(labels_tcga["SAMPLE_ID"]) == set(X_logcpm.columns)


# save
out_lab = Path("data_proc/tcga_labels.tsv")
labels_tcga.to_csv(out_lab, sep="\t", index=False)
print("saved labels:", out_lab, "| shape:", labels_tcga.shape)

# safety assert: perfect alignment
assert list(labels_tcga["SAMPLE_ID"]) == list(X_logcpm.columns)
print("alignment OK ", labels_tcga.shape)

print("columns after fix:", labels_tcga.columns.tolist())
