#!/usr/bin/env python
"""
3C — Prepare METABRIC expression for modeling (portfolio-grade).

"""

# 3C-1) Pick columns to read (header-only, zero risk)


import pandas as pd, numpy as np
from pathlib import Path
mb_dir = Path("data_raw/metabric/cbioportal")
expr_candidates = sorted([p for p in mb_dir.glob("data_mrna*microarray*.txt")
                          if "zscore" not in p.name.lower()])
assert expr_candidates, f"No raw microarray expression file found in {mb_dir}"
expr_path = expr_candidates[0]
print("Using expression file:", expr_path.name)

# read header only
all_cols = list(pd.read_csv(expr_path, sep="\t", comment="#", nrows=0).columns)

# gene/meta columns often look like this
meta_like = {"HUGO_SYMBOL","GENE_SYMBOL","GENE","ENTREZ_GENE_ID","ENTREZ_GENE_IDS"}
gene_meta_cols = [c for c in all_cols if c.upper().replace(".","_").replace(" ","_") in meta_like]

# our matched sample list from 3B-3
keep_samples = pd.read_csv("data_proc/metabric_matched_samples.txt")["SAMPLE_ID"].tolist()
cols_set = set(all_cols)
keep_samples = [c for c in keep_samples if c in cols_set]  # intersection, preserves order

usecols = gene_meta_cols + keep_samples
print(f"Meta cols: {gene_meta_cols}")
print(f"Samples to load: {len(keep_samples)} / header samples {len(all_cols) - len(gene_meta_cols)}")




# 3C-2) Load just those columns (memory-friendly)

# build dtype map: strings for gene cols, float32 for expression
dtype_map = {c: "string" for c in gene_meta_cols}
for c in keep_samples: dtype_map[c] = np.float32

expr_df = pd.read_csv(
    expr_path, sep="\t", comment="#",
    usecols=usecols, dtype=dtype_map,
    na_values=["NA","na","NaN","#N/A","null","NULL"]
)
print("raw expr_df shape:", expr_df.shape)
print(expr_df.head(2))





# 3C-3) Clean gene symbols & collapse duplicates

# decide which column is the symbol
cands = [c for c in gene_meta_cols if c.lower().startswith(("hugo","gene"))]
assert cands, f"Could not find a gene symbol column among {gene_meta_cols}"
sym_col = cands[0]

# basic cleaning
expr_df[sym_col] = expr_df[sym_col].astype("string").str.strip().str.upper()
expr_df = expr_df.dropna(subset=[sym_col])
expr_df = expr_df[expr_df[sym_col] != ""]  # drop empty symbol rows

# separate numeric matrix
X = expr_df.drop(columns=[c for c in gene_meta_cols if c != sym_col]).copy()
X = X.rename(columns={sym_col: "SYMBOL"})
# groupby SYMBOL and take median across samples
Xg = (X.groupby("SYMBOL", as_index=True)
        .median(numeric_only=True))

print("post-collapse shape (genes x samples):", Xg.shape)
# drop genes that are entirely NA after collapse
Xg = Xg.dropna(how="all")
print("after dropping all-NA genes:", Xg.shape)

# tiny QA
na_frac = Xg.isna().mean(axis=1).describe()
print("per-gene NA fraction summary:\n", na_frac)

# confirm we still have it
print(Xg.shape)
print(Xg.dtypes.iloc[:5])
Xg = Xg.astype("float32")
print(Xg.dtypes.iloc[:5])  # should now show float32

mem_mb = Xg.memory_usage(deep=True).sum() / 1e6
print("in-memory size ~", round(mem_mb, 1), "MB")



# 3C-4) Save compact matrix + aligned labels (no leakage)
# write expression as Parquet (fast I/O, smaller on disk)
out_expr = Path("data_proc/metabric_expr_raw.parquet")
Xg.to_parquet(out_expr, index=True)
print("saved:", out_expr, "size (MB) ~", round(out_expr.stat().st_size / 1e6, 1))

# build labels aligned to columns order (sample axis)
samp_surv = pd.read_csv("data_proc/metabric_sample_survival.tsv", sep="\t")
labels = (samp_surv.set_index("SAMPLE_ID")
                    .loc[Xg.columns, ["PATIENT_ID","os_event","os_time_months"]]
                    .reset_index()
                    .rename(columns={"index":"SAMPLE_ID"}))

out_lab = Path("data_proc/metabric_labels.tsv")
labels.to_csv(out_lab, sep="\t", index=False)
print("saved labels:", out_lab, "shape:", labels.shape)

# final sanity: sets must match perfectly
assert list(labels["SAMPLE_ID"]) == list(Xg.columns), "Order mismatch between labels and expression" 
print("alignment OK ")





# 3C-5) (optional) quick model-readiness checks
# % missing overall (should be low)
overall_na = float(Xg.isna().mean().mean())
print("% missing overall:", round(overall_na*100, 3))

# drop genes with too many missings (e.g., >10%) — we'll hold off for now,
# but here’s the one-liner if you want to prune:
# Xg = Xg.loc[Xg.isna().mean(axis=1) <= 0.10]
