#!/usr/bin/env python
"""
Build METABRIC patient-level OS labels from cBioPortal clinical patient file.

Inputs (from brca_metabric.tar.gz):
  - data_raw/metabric/cbioportal/data_clinical_patient.txt  (typical cBioPortal format)

Outputs:
  - data_proc/metabric_clinical_core.tsv  (PATIENT_ID, os_event[int8], os_time_months[float32])

Notes:
  - OS_STATUS values like 'DECEASED'/'LIVING' or 'DECEASED:1'/'LIVING:0' are handled.
  - OS_MONTHS is used as survival time; rows with missing time are dropped.
"""

import pandas as pd, numpy as np
from pathlib import Path

mb_dir = Path("data_raw/metabric/cbioportal")
out_path = Path("data_proc/metabric_clinical_core.tsv")
out_path.parent.mkdir(parents=True, exist_ok=True)

# --- locate the patient clinical file robustly ---
cands = sorted([p for p in mb_dir.glob("data_clinical_patient*") if p.suffix.lower() in {".txt", ".tsv"}])
assert len(cands) >= 1, f"No METABRIC patient clinical file found in {mb_dir}"
clin_path = cands[0]
print("Using patient clinical file:", clin_path.name)

# --- read; skip cBioPortal header comments that start with '#'
clin_raw = pd.read_csv(clin_path, sep="\t", comment="#", dtype=str)

# helper to find columns ignoring case/dots/underscores/spaces
def find_col(df, names):
    cols = {c.lower().replace(".", "_").replace(" ", "_"): c for c in df.columns}
    for n in names:
        key = n.lower().replace(".", "_").replace(" ", "_")
        if key in cols: 
            return cols[key]
    return None

col_patient = find_col(clin_raw, ["PATIENT_ID","PATIENT_IDENTIFIER","ID","SUBJECT_ID"])
col_os_stat = find_col(clin_raw, ["OS_STATUS","OS_STATUS__DECEASED_LIVING","OVERALL_SURVIVAL_STATUS"])
col_os_mo   = find_col(clin_raw, ["OS_MONTHS","OVERALL_SURVIVAL_MONTHS"])

assert col_patient and col_os_stat and col_os_mo, (
    f"Missing required columns. Found patient={col_patient}, os_status={col_os_stat}, os_months={col_os_mo}"
)

clin = pd.DataFrame({
    "PATIENT_ID": clin_raw[col_patient].astype(str).str.strip()
})
# status → event
status = clin_raw[col_os_stat].astype(str).str.upper().str.strip()
# robust mapping: accept 'DECEASED', 'DEAD', or strings containing those
clin["os_event"] = status.str.contains("DECEASED|DEAD", regex=True).astype("int8")

# time (months)
time = pd.to_numeric(clin_raw[col_os_mo], errors="coerce")
clin["os_time_months"] = time

# drop missing/invalid time, clip zeros to tiny positive to avoid issues
clin = clin.dropna(subset=["os_time_months"]).copy()
clin["os_time_months"] = clin["os_time_months"].clip(lower=1e-6)

clin.to_csv(out_path, sep="\t", index=False)
print("wrote:", out_path, "shape:", clin.shape)
print(clin["os_event"].value_counts(dropna=False))
print(clin[["os_event","os_time_months"]].describe(include="all"))
