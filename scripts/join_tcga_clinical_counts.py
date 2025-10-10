#!/usr/bin/env python
"""
Join TCGA clinical OS labels to STAR count files for Primary Tumor samples.

Inputs:
  - data_proc/tcga_clinical_core.tsv
  - data_proc/tcga_file_metadata.tsv

Output:
  - data_proc/tcga_survival_join.tsv
"""

import numpy as np
import pandas as pd

clin = pd.read_csv("data_proc/tcga_clinical_core.tsv", sep="\t")
meta = pd.read_csv("data_proc/tcga_file_metadata.tsv", sep="\t")

def clean(x): 
    return str(x).strip().upper() if pd.notna(x) else x

# normalize ids/labels
for c in ["submitter_id","case_id"]:
    if c in clin:  clin[c]  = clin[c].map(clean)
    if c in meta:  meta[c]  = meta[c].map(clean)
if "sample_type" in meta:
    meta["sample_type"] = meta["sample_type"].map(clean)
if "sample_type_id" in meta:
    meta["sample_type_id"] = meta["sample_type_id"].astype(str).str.strip()

# --- robust PRIMARY filter ---
if "sample_type_id" in meta and meta["sample_type_id"].notna().any() and (meta["sample_type_id"]=="01").any():
    meta01 = meta[meta["sample_type_id"]=="01"].copy()
else:
    # accept "PRIMARY SOLID TUMOR" and "PRIMARY TUMOR"
    st = meta["sample_type"].str.replace(r"\s+", " ", regex=True)
    meta01 = meta[st.isin(["PRIMARY SOLID TUMOR","PRIMARY TUMOR"])].copy()

print("meta01 rows after filtering to primary:", len(meta01))

# derive patient_id from submitter barcodes (TCGA-XX-YYYY-... -> TCGA-XX-YYYY)
def to_patient(x):
    parts = str(x).split("-")
    return "-".join(parts[:3]) if len(parts)>=3 else np.nan

clin["patient_id"]  = clin["submitter_id"].map(to_patient)
meta01["patient_id"] = meta01["submitter_id"].map(to_patient)

# try UUID join first if there is overlap, else patient_id join
has_case_overlap = "case_id" in clin and "case_id" in meta01 and \
                   len(set(clin["case_id"].dropna()) & set(meta01["case_id"].dropna())) > 0

if has_case_overlap:
    joined = (clin.drop_duplicates("case_id")
                   .merge(meta01.dropna(subset=["case_id"])
                                 .drop_duplicates("case_id")[["case_id","file_id","file_name"]],
                          on="case_id", how="inner"))
    print("joined via case_id:", joined.shape)
else:
    joined = (clin.drop_duplicates("patient_id")
                   .merge(meta01.dropna(subset=["patient_id"])
                                 .drop_duplicates("patient_id")[["patient_id","file_id","file_name"]],
                          on="patient_id", how="inner"))
    print("joined via patient_id:", joined.shape)

joined.to_csv("data_proc/tcga_survival_join.tsv", sep="\t", index=False)
print("saved -> data_proc/tcga_survival_join.tsv")

