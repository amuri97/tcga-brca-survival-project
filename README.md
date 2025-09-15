# TCGA-BRCA Survival Risk Stratification (with External Validation on METABRIC)

**Goal:** Build and validate survival risk models for breast cancer using TCGA-BRCA RNA‑seq (counts ⇒ normalized) and externally validate on METABRIC (cBioPortal/MSK). This repository is structured for **reproducibility**, **interpretability**, and a clear hiring-manager story.

## Why this project?
- Demonstrates clinical awareness (censoring, survival analysis)
- Uses a large public cohort (TCGA-BRCA) and **external validation** (METABRIC)
- Balances statistical baselines (Cox) with ML (Random Survival Forests)
- Adds biological insight (top genes, pathways) in a transparent way

## Data Sources
- **TCGA-BRCA** RNA-seq counts + clinical (GDC Portal)
- **METABRIC** expression + clinical (cBioPortal/MSK)

> Raw data belongs in `data_raw/` and **should not be committed** to Git. Processed artifacts go in `data_proc/`.

## Repository Structure
```
repo/
  data_raw/            # untouched downloads (GDC, cBioPortal) — not tracked by git
  data_proc/           # processed matrices/mappings (safe to track if small)
  notebooks/
    01_cohort_and_download.ipynb
    02_preprocessing_and_eda.ipynb
    03_modeling_internal_cv.ipynb
    04_external_validation_metabric.ipynb
  src/
    preprocess/        # Python modules for cleaning/normalization/transforms
    modeling/          # model pipelines, metrics, plotting helpers
  reports/
    figs/              # generated figures
    tables/            # generated tables
  environment.yml
  LICENSE
  .gitignore
  README.md
```

## Quickstart (Conda recommended)
1. Install **Miniconda/Anaconda**.
2. Create and activate the environment:
   ```bash
   conda env create -f environment.yml
   conda activate brca-survival
   python -m ipykernel install --user --name brca-survival --display-name "Python (brca-survival)"
   ```
3. Launch JupyterLab and verify imports:
   ```bash
   jupyter lab
   # In a notebook cell:
   # from sksurv.metrics import concordance_index_censored
   ```

## Milestones
- **M1 (Data on disk):** TCGA counts + clinical; METABRIC expression + clinical
- **M2 (Preprocessing):** gene mapping, filtering, normalization, cross-platform alignment
- **M3 (Modeling):** Cox (Elastic Net), RSF; internal CV
- **M4 (External Validation):** Evaluate on METABRIC; calibration + time-dependent AUCs
- **M5 (Explain & Ship):** Interpretability, pathway hints; README/report polish

## Repro Tips
- Keep raw files in `data_raw/` only
- Save deterministic seeds for splits
- Log parameters and metrics in notebooks (or MLflow if added later)

## Author
- Amith Murikinati — Affiliation — 2025-09-12

## License
[MIT](LICENSE)
