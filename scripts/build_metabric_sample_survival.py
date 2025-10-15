#!/usr/bin/env python
"""
Build METABRIC sample-level survival labels by mapping SAMPLE_ID -> PATIENT_ID and attaching patient OS.

Inputs:
  - data_raw/cbioportal/data_clinical_sample.txt
  - data_proc/metabric_clinical_core.tsv

Outputs:
  - data_proc/metabric_sample_survival.tsv (SAMPLE_ID, PATIENT_ID, os_event[int8], os_time_months[float32])
"""

import pandas as pd
from pathlib import Path

mb_dir = Path("data_raw/metabric/cbioportal")

# 1) locate sample-level clinical file
cands = sorted([p for p in mb_dir.glob("data_clinical_sample*") if p.suffix.lower() in {".txt", ".tsv"}])
assert len(cands) >= 1, f"No METABRIC sample clinical file found in {mb_dir}"
samp_path = cands[0]
print("Using sample clinical file:", samp_path.name)

# 2) read (skip lines starting with '#')
samp_raw = pd.read_csv(samp_path, sep="\t", comment="#", dtype=str)

def find_col(df, names):
    cols = {c.lower().replace(".", "_").replace(" ", "_"): c for c in df.columns}
    for n in names:
        key = n.lower().replace(".", "_").replace(" ", "_")
        if key in cols: 
            return cols[key]
    return None

col_sample  = find_col(samp_raw, ["SAMPLE_ID","SAMPLE_IDENTIFIER","SAMPLE"])
col_patient = find_col(samp_raw, ["PATIENT_ID","PATIENT_IDENTIFIER","PATIENT"])

assert col_sample and col_patient, f"Needed columns not found; got sample={col_sample}, patient={col_patient}"

sample_map = samp_raw[[col_sample, col_patient]].rename(columns={
    col_sample: "SAMPLE_ID",
    col_patient: "PATIENT_ID"
}).dropna().drop_duplicates()

# 3) bring in patient-level OS labels from 3B-1
clin = pd.read_csv("data_proc/metabric_clinical_core.tsv", sep="\t")

# 4) join sample→patient→OS
samp_surv = (sample_map.merge(clin, on="PATIENT_ID", how="inner")
                        .drop_duplicates(subset=["SAMPLE_ID"]))

samp_surv.to_csv("data_proc/metabric_sample_survival.tsv", sep="\t", index=False)
print("wrote: data_proc/metabric_sample_survival.tsv  shape:", samp_surv.shape)
print("unique patients:", samp_surv["PATIENT_ID"].nunique(), 
      "| unique samples:", samp_surv["SAMPLE_ID"].nunique())
print("os_event counts:\n", samp_surv["os_event"].value_counts(dropna=False))
