# Data sources & filters

## GDC (TCGA-BRCA)
- Portal: https://portal.gdc.cancer.gov/
- Filters applied (Cohort Builder):
  - Project: TCGA-BRCA
  - Samples → Sample Type: Primary Tumor (01)
  - Data Category: Transcriptome Profiling
  - Data Type: Gene Expression Quantification
  - Workflow Type: STAR - Counts
  - File Format: TSV
- Note from file headers: gene-model: GENCODE v36

## cBioPortal (METABRIC)
- Portal: https://www.cbioportal.org/
- Study: Breast Cancer (METABRIC, Nature 2012 & Nat Commun 2016)
- Bundle: brca_metabric.tar.gz (stored at data_raw/cbioportal/, not tracked in Git)

