#!/usr/bin/env python3
"""Build a documented binary ClinVar variant-classification table.

This is a research/education benchmark, not a clinical classifier. Downloads the
official NCBI ClinVar variant_summary release and retains high-review-status,
single-nucleotide germline records with unambiguous aggregate classifications.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

SOURCE_URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"
ALLOWED_LABELS = {
    "Pathogenic": "pathogenic",
    "Likely pathogenic": "pathogenic",
    "Benign": "benign",
    "Likely benign": "benign",
}
# Require evidence transparency and agreement, or expert-panel/practice-guideline review.
ACCEPTED_REVIEW_STATUSES = {
    "practice guideline",
    "reviewed by expert panel",
    "criteria provided, multiple submitters, no conflicts",
}
REQUIRED_COLUMNS = {
    "Type", "Name", "ClinicalSignificance", "ReviewStatus", "Assembly",
    "Chromosome", "Start", "Stop", "ReferenceAllele", "AlternateAllele",
    "VariationID", "AlleleID",
}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def download(url: str, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    request = urllib.request.Request(url, headers={"User-Agent": "mutation-screening-ml/1.0 (research dataset preparation)"})
    with urllib.request.urlopen(request, timeout=90) as response, destination.open("wb") as out:
        while True:
            chunk = response.read(1024 * 1024)
            if not chunk:
                break
            out.write(chunk)
    if destination.stat().st_size < 100_000:
        destination.unlink(missing_ok=True)
        raise RuntimeError("Downloaded file is unexpectedly small; refusing to process it.")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--download-dir", default="data/raw/clinvar")
    parser.add_argument("--output", default="data/processed/clinvar_binary_snv.csv")
    parser.add_argument("--url", default=SOURCE_URL, help="Official ClinVar source URL (override for archived release)")
    args = parser.parse_args()

    raw_path = Path(args.download_dir) / "variant_summary.txt.gz"
    output_path = Path(args.output)
    print(f"Downloading official ClinVar summary: {args.url}", flush=True)
    download(args.url, raw_path)
    print(f"Downloaded {raw_path.stat().st_size:,} bytes; SHA-256 {sha256_file(raw_path)}", flush=True)

    frame = pd.read_csv(raw_path, sep="\t", compression="gzip", low_memory=False, dtype=str)
    missing = sorted(REQUIRED_COLUMNS - set(frame.columns))
    if missing:
        raise ValueError(f"ClinVar file schema changed; missing required columns: {missing}")

    initial_rows = len(frame)
    # Exact accepted labels only: excludes conflicting, uncertain, risk-factor,
    # drug-response, and records not directly classified as a single variant.
    frame = frame[frame["ClinicalSignificance"].isin(ALLOWED_LABELS)].copy()
    after_label_filter = len(frame)
    frame = frame[frame["ReviewStatus"].isin(ACCEPTED_REVIEW_STATUSES)].copy()
    after_review_filter = len(frame)
    frame = frame[frame["Assembly"].eq("GRCh38")].copy()
    after_assembly_filter = len(frame)
    frame = frame[frame["Type"].eq("single nucleotide variant")].copy()
    after_type_filter = len(frame)

    # Keep simple, unambiguous SNVs only. Do not use coordinates, gene IDs,
    # clinical significance text, review status, or accession IDs as model features.
    frame["ReferenceAllele"] = frame["ReferenceAllele"].str.upper()
    frame["AlternateAllele"] = frame["AlternateAllele"].str.upper()
    frame = frame[
        frame["ReferenceAllele"].str.fullmatch(r"[ACGT]", na=False)
        & frame["AlternateAllele"].str.fullmatch(r"[ACGT]", na=False)
        & frame["ReferenceAllele"].ne(frame["AlternateAllele"])
    ].copy()

    # A variant may be represented by multiple aggregate rows. Collapse duplicate
    # allele IDs to one row and discard any allele ID with conflicting binary labels.
    frame["label"] = frame["ClinicalSignificance"].map(ALLOWED_LABELS)
    frame["AlleleID"] = frame["AlleleID"].fillna("").astype(str)
    frame = frame[frame["AlleleID"].ne("")].copy()
    label_counts = frame.groupby("AlleleID")["label"].nunique()
    conflicting_ids = set(label_counts[label_counts > 1].index)
    frame = frame[~frame["AlleleID"].isin(conflicting_ids)].copy()
    frame = frame.sort_values(["AlleleID", "ReviewStatus"], kind="stable")
    frame = frame.drop_duplicates("AlleleID", keep="first")

    # Intentionally modest, transparent features for a baseline—not sufficient
    # for clinical interpretation. No coordinates or label-proxy fields.
    result = pd.DataFrame({
        "reference_allele": frame["ReferenceAllele"],
        "alternate_allele": frame["AlternateAllele"],
        "variant_type": "SNV",
        "ref_gc": frame["ReferenceAllele"].isin(["G", "C"]).astype(int),
        "alt_gc": frame["AlternateAllele"].isin(["G", "C"]).astype(int),
        "transition": (
            frame["ReferenceAllele"].isin(["A", "G"])
            & frame["AlternateAllele"].isin(["A", "G"])
        ) | (
            frame["ReferenceAllele"].isin(["C", "T"])
            & frame["AlternateAllele"].isin(["C", "T"])
        ),
        "label": frame["label"].values,
        "clinvar_allele_id": frame["AlleleID"].values,
        "clinvar_variation_id": frame["VariationID"].values,
        "clinical_significance_source": frame["ClinicalSignificance"].values,
        "review_status_source": frame["ReviewStatus"].values,
        "chromosome_source": frame["Chromosome"].values,
        "position_start_source": frame["Start"].values,
    })
    # The source/provenance columns are kept in the prepared file for auditability.
    # The training script must explicitly select only the feature columns above.
    output_path.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(output_path, index=False)

    label_counts_out = result["label"].value_counts().to_dict()
    manifest = {
        "dataset": "NCBI ClinVar variant_summary",
        "source_url": args.url,
        "retrieved_utc": datetime.now(timezone.utc).isoformat(),
        "source_file": str(raw_path),
        "source_file_bytes": raw_path.stat().st_size,
        "source_file_sha256": sha256_file(raw_path),
        "output_file": str(output_path),
        "output_rows": len(result),
        "label_counts": label_counts_out,
        "filters": {
            "assembly": "GRCh38",
            "variant_type": "single nucleotide variant",
            "clinical_significance": list(ALLOWED_LABELS),
            "binary_mapping": ALLOWED_LABELS,
            "review_status": sorted(ACCEPTED_REVIEW_STATUSES),
            "alleles": "single A/C/G/T reference and alternate bases; reference != alternate",
            "duplicates": "one row per AlleleID; IDs with conflicting binary labels excluded",
        },
        "rows_at_stages": {
            "raw": initial_rows,
            "accepted_binary_labels": after_label_filter,
            "accepted_review_status": after_review_filter,
            "GRCh38": after_assembly_filter,
            "SNV_type": after_type_filter,
            "final_after_allele_and_duplicate_filters": len(result),
        },
        "model_features": ["reference_allele", "alternate_allele", "variant_type", "ref_gc", "alt_gc", "transition"],
        "label_column": "label",
        "limitations": [
            "ClinVar classifications are curated assertions, not perfect ground truth.",
            "This is a simple SNV-only baseline and omits sequence context and functional annotations.",
            "Random row-level validation can overestimate generalization; use gene- or locus-held-out evaluation for a stronger benchmark.",
            "Not for diagnosis, treatment decisions, or clinical use.",
        ],
        "attribution": "NCBI ClinVar; cite Landrum et al., Nucleic Acids Research (2025), doi:10.1093/nar/gkae1090, and retain ClinVar record identifiers.",
        "license_and_policy": "ClinVar is freely accessible; NCBI requests attribution when data are redistributed. Review current ClinVar data-use policy before redistribution.",
    }
    manifest_path = output_path.with_suffix(".manifest.json")
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    print(f"Prepared {len(result):,} rows -> {output_path}")
    print(f"Labels: {label_counts_out}")
    print(f"Audit manifest -> {manifest_path}")
    if len(result) == 0 or len(label_counts_out) < 2:
        print("WARNING: output is empty or has only one class; inspect source schema and filters before training.", file=sys.stderr)
        return 2
    print("Reminder: the included source/provenance columns are NOT model features.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
