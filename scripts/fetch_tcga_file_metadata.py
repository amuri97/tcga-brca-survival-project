#!/usr/bin/env python
"""
Fetch GDC file metadata for STAR counts to support Primary Tumor filtering and joins.

Inputs:
  - data_raw/gdc/manifest_tcga_brca_star_counts_primary.txt

Output:
  - data_proc/tcga_file_metadata.tsv
"""
import pandas as pd, requests, json, math, pathlib

manifest_path = pathlib.Path("data_raw/gdc/manifest_tcga_brca_star_counts_primary.txt")
out_path = pathlib.Path("data_proc/tcga_file_metadata.tsv")
out_path.parent.mkdir(parents=True, exist_ok=True)

# read file UUIDs from manifest (first column is 'id')
ids = pd.read_csv(manifest_path, sep="\t")["id"].tolist()

fields = [
    "file_id","file_name","data_category","data_type","experimental_strategy",
    "analysis.workflow_type","cases.samples.sample_type","cases.samples.sample_type_id",
    "cases.case_id","cases.submitter_id"
]
url = "https://api.gdc.cancer.gov/files"
rows = []
BATCH = 200

for i in range(0, len(ids), BATCH):
    chunk = ids[i:i+BATCH]
    payload = {
        "filters": {"op":"in","content":{"field":"files.file_id","value":chunk}},
        "format": "JSON",
        "fields": ",".join(fields),
        "size": BATCH
    }
    r = requests.post(url, headers={"Content-Type":"application/json"}, data=json.dumps(payload), timeout=120)
    r.raise_for_status()
    for hit in r.json()["data"]["hits"]:
        rows.append({
            "file_id": hit.get("file_id"),
            "file_name": hit.get("file_name"),
            "data_category": hit.get("data_category"),
            "data_type": hit.get("data_type"),
            "experimental_strategy": hit.get("experimental_strategy"),
            "workflow_type": (hit.get("analysis") or {}).get("workflow_type"),
            "sample_type": ((hit.get("cases") or [{}])[0].get("samples") or [{}])[0].get("sample_type"),
            "sample_type_id": ((hit.get("cases") or [{}])[0].get("samples") or [{}])[0].get("sample_type_id"),
            "case_id": ((hit.get("cases") or [{}])[0]).get("case_id"),
            "submitter_id": ((hit.get("cases") or [{}])[0]).get("submitter_id"),
        })

meta = pd.DataFrame(rows)
meta.to_csv(out_path, sep="\t", index=False)
print("wrote", out_path, "shape", meta.shape)
print(meta.groupby(["sample_type","workflow_type","data_type"]).size().sort_values(ascending=False).head(10))
print(meta['sample_type_id'].value_counts())  # expect '01' > '11' > '06' tiny
