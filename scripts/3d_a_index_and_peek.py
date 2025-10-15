#!/usr/bin/env python
"""
3D-A: Index STAR/HTSeq count files and do a minimal 'peek' sanity check.
- Confirms files are present under data_raw/gdc_star_counts_primary/**
- Prints a few discovered files, shows first lines of one example
- Prints header columns to help confirm the raw counts column name
"""

import pandas as pd, numpy as np
from pathlib import Path

base = Path("data_raw/gdc_star_counts_primary")
assert base.exists(), f"Counts folder not found: {base}"

# 1) map every *.tsv or *.tsv.gz file by its basename
name_to_path = {}
for p in base.rglob("*"):
    if p.is_file() and (p.suffix.lower() in [".tsv",".gz"] or p.name.endswith(".tsv.gz")):
        name_to_path[p.name] = p

print("Indexed files:", len(name_to_path))

# 2) load the join table (clinical+files)
joined = pd.read_csv("data_proc/tcga_survival_join.tsv", sep="\t")
print("joined rows:", joined.shape)

# 3) pick one file from the join list to sniff
ex = joined.iloc[0]
ex_path = name_to_path.get(ex["file_name"])
print("example path:", ex_path)

# 4) peek first few lines/columns
kw = dict(sep="\t", comment="#", nrows=5, engine="python")
if str(ex_path).endswith(".gz"):
    kw["compression"] = "infer"
df_head = pd.read_csv(ex_path, **kw)
print("example columns:", df_head.columns.tolist())
print(df_head.head(3))

