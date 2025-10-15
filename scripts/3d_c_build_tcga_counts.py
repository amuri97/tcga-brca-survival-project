#!/usr/bin/env python
"""
3D-C: Build TCGA raw counts matrix from STAR/HTSeq outputs (Primary Tumor samples).

Inputs:
  - data_proc/tcga_survival_join.tsv  (submitter_id ↔ file_name mapping)
  - data_raw/gdc_star_counts_primary/** (downloaded STAR/HTSeq files)

Outputs:
  - data_proc/tcga_counts_raw.parquet  (genes × samples, int32)

Behavior (faithful to REPL):
  - Map file_name -> actual path
  - read_star_counts() per sample (raw counts only, drop summary rows)
  - Assemble genes × samples, collapse duplicate symbols by median
  - Reorder columns strictly to join order
"""


import pandas as pd, numpy as np
from pathlib import Path


def read_star_counts(path):
    """Return Series: index = gene symbol (uppercased; fallback to gene_id), values = raw counts."""
    import pandas as pd

    kw = dict(sep="\t", comment="#", dtype=str, engine="c")
    if str(path).endswith(".gz"):
        kw["compression"] = "infer"
    df = pd.read_csv(path, **kw)

    # normalize column names
    norm = {c: c.lower().replace(".", "_").strip() for c in df.columns}
    inv  = {v: k for k, v in norm.items()}

    # columns present?
    gene_id_col   = inv.get("gene_id")
    gene_name_col = inv.get("gene_name") or inv.get("gene") or inv.get("hugo_symbol")

    # choose raw-counts column (NOT tpm/fpkm)
    for key in ["unstranded", "raw_count", "read_count", "htseq_counts", "stranded_first", "stranded_second"]:
        if key in inv:
            count_col = inv[key]
            break
    else:
        raise ValueError(f"No counts column found in {list(df.columns)}")

    # build a robust index: use gene_name if available; otherwise fall back to gene_id
    if gene_name_col is not None:
        sym = df[gene_name_col].astype("string").str.strip()
        # consider empty/NA/"NA"/"NAN"/"-" as missing
        miss = sym.isna() | (sym == "") | sym.str.upper().isin({"NA", "NAN", "NULL", "-"})
        if gene_id_col is not None:
            sym.loc[miss] = df[gene_id_col]
        else:
            sym = sym[~miss]
        index = sym.str.upper()
    elif gene_id_col is not None:
        index = df[gene_id_col].astype("string").str.strip().str.upper()
    else:
        raise ValueError("Neither gene_name nor gene_id found")

    # make the series
    s = pd.Series(df[count_col].values, index=index, name=count_col)
    print("read_star_counts:", path, "n_rows=", len(s), flush=True)


    # drop STAR/HTSeq summary rows (single pass; robust to NA)
    idx = s.index.astype("string")
    keep = ~idx.str.startswith(("__", "N_"), na=False)
    s = s[keep]

    # coerce to numeric; NaN -> 0
    s = pd.to_numeric(s, errors="coerce").fillna(0)

    # remove any residual empty/NAN index rows (defensive)
    s = s[(s.index != "") & (s.index != "NAN")]

    return s




# where your STAR counts live
base = Path("data_raw/gdc_star_counts_primary")
assert base.exists(), f"Counts folder not found: {base}"

# map basename -> full path for fast lookup
exts = (".tsv", ".tsv.gz", ".txt", ".txt.gz")
name_to_path = {
    p.name: p
    for p in base.rglob("*")
    if p.is_file() and any(p.name.lower().endswith(e) for e in exts)
}
print("Indexed files:", len(name_to_path))

# joined table from 2D-3 (submitter_id, file_name, os labels)
joined = pd.read_csv("data_proc/tcga_survival_join.tsv", sep="\t")
joined["submitter_id"] = joined["submitter_id"].astype(str).str.upper().str.strip()
joined = joined.drop_duplicates(subset=["submitter_id"]).copy()
sample_order = joined["submitter_id"].tolist()
print("Joined rows (unique samples):", len(sample_order))



series_by_sample = {}
missing_files = []

for i, row in joined.iterrows():
    sid = row["submitter_id"]
    f   = row["file_name"]
    p   = name_to_path.get(f)
    if p is None:
        missing_files.append(f)
        continue
    s = read_star_counts(p)           # <-- the function you just perfected
    series_by_sample[sid] = s
    if (i + 1) % 100 == 0:
        print(f" loaded {len(series_by_sample)}/{len(sample_order)} ...")

print("\nloaded samples:", len(series_by_sample), " | missing files:", len(missing_files))
if missing_files:
    print("example missing:", missing_files[:5])



# assemble (union of all gene indices across samples)
Xc = pd.DataFrame(series_by_sample)

# how many duplicated gene symbols?
n_dupe = int(Xc.index.duplicated().sum())
print("genes before collapse:", Xc.shape[0], "| duplicate rows:", n_dupe)

# collapse duplicates by median across duplicate-symbol rows
Xc = Xc.groupby(level=0, as_index=True).median(numeric_only=True)

# fill any gaps with 0 counts
Xc = Xc.fillna(0)

# order columns to match joined sample order (and drop any missing)
cols_loaded = [c for c in sample_order if c in Xc.columns]
Xc = Xc.reindex(columns=cols_loaded)

print("raw counts shape after collapse (genes × samples):", Xc.shape)
print("any negative counts?", bool((Xc < 0).any().any()))
print("library sizes, first 5:", list(Xc.sum(axis=0).astype(int).iloc[:5]))




out_counts = Path("data_proc/tcga_counts_raw.parquet")
Xc.astype("int32").to_parquet(out_counts, index=True)
print("saved raw counts ->", out_counts)



