#!/usr/bin/env python
"""
Sanity diagnostics for Phase 2 artifacts.
"""
import pandas as pd
m = pd.read_csv("data_proc/tcga_file_metadata.tsv", sep="\t")
j = pd.read_csv("data_proc/tcga_survival_join.tsv", sep="\t")

print("metadata sample_type counts:")
print(m["sample_type"].value_counts(dropna=False).head(10))
print("\njoin shape:", j.shape)
print("\nlabels summary:")
print(j[["os_event","os_time_months"]].describe(include="all"))
