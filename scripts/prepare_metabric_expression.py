#!/usr/bin/env python
"""
Prepare METABRIC expression:
  - Read raw microarray expression (cBioPortal 'data_mrna_illumina_microarray.txt')
  - Detect gene symbol column; drop non-expression columns
  - Align columns to sample survival table; report coverage
  - Collapse duplicate symbols by median
  - Save parquet (float32) and a text list of matched sample IDs

Inputs:
  - data_raw/cbioportal/data_mrna_illumina_microarray.txt
  - data_proc/metabric_sample_survival.tsv

Outputs:
  - data_proc/metabric_expr_raw.parquet  (genes x matched_samples)
  - data_proc/metabric_matched_samples.txt
"""

import pandas as pd
from pathlib import Path

mb_dir = Path("data_raw/metabric/cbioportal")

# 1) Find the RAW microarray expression file (not the z-scores one)
expr_candidates = sorted([p for p in mb_dir.glob("data_mrna*microarray*.txt")
                          if "zscore" not in p.name.lower()])
assert expr_candidates, f"No raw microarray expression file found in {mb_dir}"
expr_path = expr_candidates[0]
print("Using expression file:", expr_path.name)

# 2) Read only the header row (cheap)
cols = list(pd.read_csv(expr_path, sep="\t", comment="#", nrows=0).columns)

# 3) Separate gene/meta columns vs sample columns
meta_like = {"HUGO_SYMBOL","GENE_SYMBOL","GENE","ENTREZ_GENE_ID","ENTREZ_GENE_IDS"}
gene_meta_cols = [c for c in cols if c.upper().replace(".","_").replace(" ","_") in meta_like]
sample_cols = [c for c in cols if c not in gene_meta_cols]

print(f"Total columns: {len(cols)} | gene/meta: {len(gene_meta_cols)} | sample columns: {len(sample_cols)}")
print("Gene/meta columns guessed:", gene_meta_cols)

# 4) Compare against our survival table SAMPLE_IDs
samp_surv = pd.read_csv("data_proc/metabric_sample_survival.tsv", sep="\t")
expr_samples = set(map(str, sample_cols))
surv_samples = set(map(str, samp_surv["SAMPLE_ID"]))

inter = expr_samples & surv_samples
miss_in_expr = surv_samples - expr_samples   # samples we have labels for but no expression
miss_in_surv = expr_samples - surv_samples   # samples in matrix but no labels

print("Matched samples:", len(inter))
print("Missing in expression (should be ~0):", len(miss_in_expr))
print("Missing in survival (ok to be 0):", len(miss_in_surv))
print("Example matched:", sorted(list(inter))[:5])
print("Example missing_in_expr:", sorted(list(miss_in_expr))[:5])
print("Example missing_in_surv:", sorted(list(miss_in_surv))[:5])

# 5) Save the matched sample list for later slice-loading
pd.Series(sorted(inter), name="SAMPLE_ID").to_csv("data_proc/metabric_matched_samples.txt", index=False)
print("Saved matched sample IDs -> data_proc/metabric_matched_samples.txt")


