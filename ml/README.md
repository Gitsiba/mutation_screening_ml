# Input schema for the optional supervised ML baseline

The trainer expects a CSV with one observation per row, one trusted target column,
and candidate predictor columns. There is intentionally no synthetic dataset here:
fabricated labels would produce misleading metrics.

Example *column names only* (not example patient data):

```text
sample_id,consequence,depth,variant_allele_fraction,gene,validated_label
```

- `validated_label` is an example target name; it must be replaced with a genuinely
  curated label appropriate to the scientific question.
- Exclude IDs, patient identifiers, sample identifiers, target-derived columns, and
  any feature that leaks the label. Some common ID-like column names are excluded
  automatically; use `--drop-columns` for additional columns.
- If multiple rows come from one patient, sample, family, or cohort, do not use a
  random row split without checking for leakage. Use a group-aware validation design
  before claiming generalization.
- A VCF alone does not provide trustworthy pathogenicity or somatic labels. Variant
  FILTER status is not a substitute for a ground-truth label.
