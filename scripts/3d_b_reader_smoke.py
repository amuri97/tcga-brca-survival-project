#!/usr/bin/env python
"""
3D-B: Smoke-test read_star_counts on one STAR/HTSeq file and print quick QC stats.
"""
import pandas as pd
from pathlib import Path
import numpy as np


def read_star_counts(path):
    """Return Series: index = gene symbol (uppercased; fallback to gene_id), values = raw counts."""
    import pandas as pd

    kw = dict(sep="\t", comment="#", dtype=str, engine="python")
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

    # drop STAR/HTSeq summary rows
    idx = s.index.astype(str)
    s = s[~(idx.str.startswith("__") | idx.str.startswith("N_"))]

    # coerce to numeric; NaN -> 0
    s = pd.to_numeric(s, errors="coerce").fillna(0)

    # remove any residual empty/NAN index rows (defensive)
    s = s[(s.index != "") & (s.index != "NAN")]

    return s




base_smoke = Path("data_raw/gdc_star_counts_primary")
assert base_smoke.exists(), f"Counts folder not found: {base}"

# 1) map every *.tsv or *.tsv.gz file by its basename
name_to_path = {}
for p in base_smoke.rglob("*"):
    if p.is_file() and (p.suffix.lower() in [".tsv",".gz"] or p.name.endswith(".tsv.gz")):
        name_to_path[p.name] = p


joinedsmoke = pd.read_csv("data_proc/tcga_survival_join.tsv", sep="\t")

exsmoke = joinedsmoke.iloc[0]
ex_path_smoke = name_to_path.get(exsmoke["file_name"])
print("example path:", ex_path_smoke)


s_test = read_star_counts(ex_path_smoke)
print("series shape:", s_test.shape, "| nonzero:", int((s_test>0).sum()))
print("NAN rows:", int((s_test.index == "NAN").sum()))
print("head:\n", s_test.head())


