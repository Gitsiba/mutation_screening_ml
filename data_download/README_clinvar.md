# ClinVar benchmark dataset

## Purpose

This pipeline creates a reproducible, binary **germline variant-classification research benchmark** from NCBI ClinVar. It is not a cancer-detection model and must not be used for clinical decisions.

## Source and access

- Source file: [NCBI ClinVar `variant_summary.txt.gz`](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz)
- Official download guidance: https://www.ncbi.nlm.nih.gov/clinvar/docs/access/
- ClinVar data-use policy: https://www.ncbi.nlm.nih.gov/clinvar/docs/maintenance_use/
- Review-status definitions: https://www.ncbi.nlm.nih.gov/clinvar/docs/review_status/
- Attribution: ClinVar requests attribution when data are copied or redistributed. Cite Landrum et al., *Nucleic Acids Research* (2025), doi:10.1093/nar/gkae1090, and retain ClinVar accessions/IDs. Review the current NCBI policy before redistribution.

ClinVar records are submitted by multiple organizations and are not independently verified by NCBI in every case. A review-status filter improves traceability but does not make labels infallible.

## Labels and filters

The script keeps only aggregate labels `Pathogenic`, `Likely pathogenic`, `Benign`, and `Likely benign`, mapped to `pathogenic` and `benign`. It excludes uncertain/conflicting and other classification types, requires selected higher-transparency review statuses, restricts to GRCh38 single-nucleotide variants with simple A/C/G/T alleles, and collapses records to one row per `AlleleID`. Allele IDs with conflicting binary labels are removed.

The script downloads the live official file by default, records retrieval time and SHA-256, and writes a JSON manifest alongside the CSV. Because ClinVar changes weekly, the resulting row counts can change. Preserve the downloaded file and manifest to reproduce a particular run. For archival reproducibility, supply an archived source URL using `--url` if one is available.

## Run dataset preparation

```bash
python -m pip install pandas
python data_download/prepare_clinvar_dataset.py
```

Outputs (generated locally; do not commit the full dataset unless you have reviewed size and data-use implications):

- `data/raw/clinvar/variant_summary.txt.gz`
- `data/processed/clinvar_binary_snv.csv`
- `data/processed/clinvar_binary_snv.manifest.json`

## Train the baseline

```bash
python -m pip install -r requirements-ml.txt
python ml/train_clinvar_baseline.py
```

This script uses only six predeclared allele-derived features and performs a chromosome-grouped holdout split. It deliberately refuses to fall back to a random row split if it cannot find both classes in train and test. It saves `clinvar_baseline.joblib`, `metrics.json`, and `heldout_predictions.csv` under `ml/clinvar_results/`.

## Columns and leakage precautions

Candidate baseline features: `reference_allele`, `alternate_allele`, `variant_type`, `ref_gc`, `alt_gc`, `transition`. Target: `label`.

The `clinvar_*`, `clinical_significance_source`, `review_status_source`, `chromosome_source`, and `position_start_source` columns are provenance/audit fields, **not model features**. In particular, never include clinical-significance text or review status in the feature matrix: they are direct proxies for the target. The dedicated training script explicitly selects only the six features and uses chromosome as a grouping variable, not a predictor.

These simple features are deliberately limited and may yield weak performance; this pipeline establishes data provenance and a reproducible starting point, not a clinically useful predictor. A stronger study should add independently sourced, versioned biological features and evaluate on held-out genes/loci, with duplicate and homology leakage controlled. Report class balance and metrics beyond accuracy. Do not claim clinical validity.
