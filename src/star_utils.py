#!/usr/bin/env python
"""
Utilities for STAR/HTSeq raw counts.

read_star_counts(path) -> pd.Series
  - Chooses raw counts column (not TPM/FPKM)
  - Prefers gene symbol; falls back to gene_id for missing/blank symbols
  - Drops summary rows (__*, N_*)
  - Returns numeric counts indexed by uppercased gene symbol/id
"""

import pandas as pd, numpy as np
from pathlib import Path

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


