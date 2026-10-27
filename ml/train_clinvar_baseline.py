#!/usr/bin/env python3
"""Train an explicitly leakage-aware baseline on the prepared ClinVar CSV.

Research/education only. Uses only the small, predeclared allele-derived feature
set and holds out chromosomes as groups. This is a demonstration baseline, not
a clinically useful or validated pathogenicity predictor.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import joblib
import pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import (
    accuracy_score,
    balanced_accuracy_score,
    classification_report,
    confusion_matrix,
    f1_score,
)
from sklearn.model_selection import GroupShuffleSplit
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder


FEATURES = [
    "reference_allele",
    "alternate_allele",
    "variant_type",
    "ref_gc",
    "alt_gc",
    "transition",
]
TARGET = "label"
GROUP = "chromosome_source"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", default="data/processed/clinvar_binary_snv.csv")
    parser.add_argument("--output-dir", default="ml/clinvar_results")
    parser.add_argument("--test-size", type=float, default=0.2)
    parser.add_argument("--random-state", type=int, default=42)
    args = parser.parse_args()

    source = Path(args.input)
    out_dir = Path(args.output_dir)
    if not source.is_file():
        parser.error(f"Input not found: {source}. Run data_download/prepare_clinvar_dataset.py first.")
    if not 0.1 <= args.test_size <= 0.5:
        parser.error("--test-size must be between 0.1 and 0.5.")

    df = pd.read_csv(source, dtype={"chromosome_source": str, "clinvar_allele_id": str})
    required = set(FEATURES + [TARGET, GROUP, "clinvar_allele_id"])
    missing = sorted(required - set(df.columns))
    if missing:
        parser.error(f"Input is missing required columns: {missing}")
    df = df.dropna(subset=FEATURES + [TARGET, GROUP, "clinvar_allele_id"]).copy()
    if df.empty:
        parser.error("No complete rows remain after filtering.")
    if df["clinvar_allele_id"].duplicated().any():
        parser.error("Duplicate ClinVar AlleleID values detected; rerun dataset preparation.")
    if df[TARGET].nunique() != 2:
        parser.error("Expected both benign and pathogenic labels.")
    if df[GROUP].nunique() < 3:
        parser.error("Need at least three distinct chromosome groups for a meaningful held-out-chromosome split.")

    # Find a deterministic group holdout where both classes occur on both sides.
    # This avoids a random row split and keeps chromosomes wholly separated.
    chosen = None
    for offset in range(200):
        splitter = GroupShuffleSplit(
            n_splits=1,
            test_size=args.test_size,
            random_state=args.random_state + offset,
        )
        train_idx, test_idx = next(splitter.split(df[FEATURES], df[TARGET], groups=df[GROUP]))
        if (
            df.iloc[train_idx][TARGET].nunique() == 2
            and df.iloc[test_idx][TARGET].nunique() == 2
        ):
            chosen = (train_idx, test_idx, args.random_state + offset)
            break
    if chosen is None:
        parser.error(
            "Could not find a chromosome-group split containing both classes in train and test. "
            "Inspect class distribution by chromosome; do not silently fall back to a row-level split."
        )
    train_idx, test_idx, split_seed = chosen
    train = df.iloc[train_idx].copy()
    test = df.iloc[test_idx].copy()

    categorical = ["reference_allele", "alternate_allele", "variant_type"]
    numeric = ["ref_gc", "alt_gc", "transition"]
    preprocessor = ColumnTransformer(
        transformers=[
            ("categorical", OneHotEncoder(handle_unknown="ignore"), categorical),
            ("numeric", "passthrough", numeric),
        ],
        remainder="drop",
    )
    model = Pipeline([
        ("features", preprocessor),
        ("classifier", RandomForestClassifier(
            n_estimators=400,
            class_weight="balanced",
            random_state=args.random_state,
            min_samples_leaf=2,
            n_jobs=-1,
        )),
    ])
    model.fit(train[FEATURES], train[TARGET])
    predictions = model.predict(test[FEATURES])
    labels = ["benign", "pathogenic"]

    out_dir.mkdir(parents=True, exist_ok=True)
    metrics = {
        "dataset": str(source),
        "target": TARGET,
        "features_used": FEATURES,
        "provenance_columns_excluded": [
            "clinvar_allele_id", "clinvar_variation_id", "clinical_significance_source",
            "review_status_source", "chromosome_source", "position_start_source",
        ],
        "split_method": "GroupShuffleSplit by chromosome_source",
        "random_state_used_for_split": split_seed,
        "train_rows": len(train),
        "test_rows": len(test),
        "train_chromosomes": sorted(train[GROUP].unique().tolist()),
        "test_chromosomes": sorted(test[GROUP].unique().tolist()),
        "train_label_counts": train[TARGET].value_counts().to_dict(),
        "test_label_counts": test[TARGET].value_counts().to_dict(),
        "accuracy": accuracy_score(test[TARGET], predictions),
        "balanced_accuracy": balanced_accuracy_score(test[TARGET], predictions),
        "macro_f1": f1_score(test[TARGET], predictions, average="macro", zero_division=0),
        "confusion_matrix_labels_benign_then_pathogenic": confusion_matrix(
            test[TARGET], predictions, labels=labels
        ).tolist(),
        "classification_report": classification_report(
            test[TARGET], predictions, labels=labels, output_dict=True, zero_division=0
        ),
        "limitations": [
            "Tiny allele-derived feature set; expected to be weak and not clinically useful.",
            "Chromosome-held-out validation is a more stringent grouping choice but is not a substitute for independent external validation.",
            "ClinVar labels are curated assertions with possible ascertainment and submission biases.",
            "Research/education only; not for diagnosis or treatment decisions.",
        ],
    }
    joblib.dump(model, out_dir / "clinvar_baseline.joblib")
    (out_dir / "metrics.json").write_text(json.dumps(metrics, indent=2, default=str) + "\n", encoding="utf-8")
    pd.DataFrame({
        "clinvar_allele_id": test["clinvar_allele_id"].values,
        "chromosome_source": test[GROUP].values,
        "true_label": test[TARGET].values,
        "predicted_label": predictions,
    }).to_csv(out_dir / "heldout_predictions.csv", index=False)
    print(json.dumps({
        "output_dir": str(out_dir),
        "train_rows": len(train),
        "test_rows": len(test),
        "train_chromosomes": metrics["train_chromosomes"],
        "test_chromosomes": metrics["test_chromosomes"],
        "accuracy": metrics["accuracy"],
        "balanced_accuracy": metrics["balanced_accuracy"],
        "macro_f1": metrics["macro_f1"],
        "warning": "Toy baseline only; do not interpret as clinical performance.",
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
