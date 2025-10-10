#!/usr/bin/env python
"""
Create core OS table from GDC clinical bundle.

Inputs:
  - data_raw/gdc/clinical/clinical.tsv
  - data_raw/gdc/clinical/follow_up.tsv (optional)

Output:
  - data_proc/tcga_clinical_core.t

What it does:
  - Finds columns (submitter_id/case_id/vital_status/days_to_death/last_follow_up)
  - Coerces numerics; merges best follow-up from follow_up.tsv if present
  - Creates os_event, os_time_months; drops rows missing time
"""
import pandas as pd, numpy as np
from pathlib import Path

clin_dir = Path("data_raw/gdc/clinical")
clin_path = clin_dir / "clinical.tsv"
fup_path  = clin_dir / "follow_up.tsv"
out_path  = Path("data_proc/tcga_clinical_core.tsv")
out_path.parent.mkdir(parents=True, exist_ok=True)

def find_col(df, targets):
    targets = [t.lower() for t in targets]
    for c in df.columns:
        lc = c.lower()
        if lc in targets: return c
        if any(lc.endswith("." + t) for t in targets): return c
    return None

# --- load main clinical ---
clin = pd.read_csv(clin_path, sep="\t", low_memory=False)

col_submitter = find_col(clin, ["submitter_id","case_submitter_id"])
col_case      = find_col(clin, ["case_id"])
col_vital     = find_col(clin, ["vital_status"])
col_dod       = find_col(clin, ["days_to_death"])
col_dlfu      = find_col(clin, ["days_to_last_follow_up","days_to_last_followup"])

assert col_submitter and col_vital, "key clinical columns missing"

for c in [col_dod, col_dlfu]:
    if c and c in clin.columns:
        clin[c] = pd.to_numeric(clin[c], errors="coerce")

core = pd.DataFrame({
    "submitter_id": clin[col_submitter].astype(str),
    "vital_status": clin[col_vital].astype(str).str.upper().str.strip()
})
if col_case: core["case_id"] = clin[col_case]
if col_dod:  core["days_to_death"] = clin[col_dod]
if col_dlfu: core["days_to_last_follow_up"] = clin[col_dlfu]

# --- optionally enhance with follow_up.tsv (use max per patient) ---
if fup_path.exists():
    fup = pd.read_csv(fup_path, sep="\t", low_memory=False)
    fup_submitter = find_col(fup, ["submitter_id","case_submitter_id"])
    fup_dlfu      = find_col(fup, ["days_to_last_follow_up","days_to_last_followup"])
    if fup_submitter and fup_dlfu:
        fup[fup_dlfu] = pd.to_numeric(fup[fup_dlfu], errors="coerce")
        fup_max = (fup[[fup_submitter, fup_dlfu]]
                   .groupby(fup_submitter, as_index=False)
                   .max()
                   .rename(columns={fup_submitter:"submitter_id", fup_dlfu:"fup_days_to_last_follow_up"}))
        core = core.merge(fup_max, on="submitter_id", how="left")
        core["days_to_last_follow_up"] = np.where(
            core.get("days_to_last_follow_up").notna(),
            core["days_to_last_follow_up"],
            core["fup_days_to_last_follow_up"]
        )
        core.drop(columns=["fup_days_to_last_follow_up"], inplace=True)

# --- build event flag robustly ---
vs = core["vital_status"]
core["os_event"] = (
    vs.isin(["DEAD","DECEASED"]) |
    core.get("days_to_death").fillna(0).gt(0)
).astype("int8")

# --- collapse to one row per patient ---
def agg_patient(g):
    # if any record indicates death, take max days_to_death
    event = int(g["os_event"].max())
    d_death = g.get("days_to_death")
    d_lfu   = g.get("days_to_last_follow_up")
    d_death_max = float(d_death.max()) if d_death is not None else np.nan
    d_lfu_max   = float(d_lfu.max())   if d_lfu is not None else np.nan
    time = d_death_max if event==1 and not np.isnan(d_death_max) else d_lfu_max
    return pd.Series({
        "case_id": g.get("case_id").iloc[0] if "case_id" in g else np.nan,
        "os_event": event,
        "os_time_days": time
    })

core_agg = core.groupby("submitter_id", as_index=False).apply(agg_patient)
core_agg["os_time_days"] = core_agg["os_time_days"].clip(lower=0)
core_agg = core_agg.dropna(subset=["os_time_days"]).copy()
core_agg["os_time_months"] = core_agg["os_time_days"] / 30.44

core_agg.to_csv(out_path, sep="\t", index=False)
print("wrote:", out_path, "shape:", core_agg.shape)
print(core_agg["os_event"].value_counts(dropna=False))
print(core_agg[["os_event","os_time_months"]].describe(include="all"))

core_agg = pd.read_csv("data_proc/tcga_clinical_core.tsv", sep="\t")
core_agg["os_time_days"] = core_agg["os_time_days"].clip(lower=1.0)
core_agg["os_time_months"] = core_agg["os_time_days"] / 30.44
core_agg.to_csv("data_proc/tcga_clinical_core.tsv", sep="\t", index=False)
