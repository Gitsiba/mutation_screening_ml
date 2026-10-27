"""Train a baseline research classifier from an explicitly labelled tabular dataset.

This script does not create labels or train on a VCF alone. Input must contain
one row per independently defined observation and a trustworthy target column.
Do not use this for clinical decisions.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
import sklearn
from sklearn.compose import ColumnTransformer
from sklearn.ensemble import RandomForestClassifier
from sklearn.impute import SimpleImputer
from sklearn.metrics import (
    accuracy_score,
    balanced_accuracy_score,
    classification_report,
    confusion_matrix,
    f1_score,
    precision_score,
    recall_score,
    roc_auc_score,
)
from sklearn.model_selection import train_test_split
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder


DEFAULT_ID_COLUMNS = {
    "id", "sample_id", "patient_id", "subject_id", "record_id",
    "variant_id", "chrom", "chromosome", "position", "pos",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Train and evaluate a baseline Random Forest on labelled tabular data."
    )
    parser.add_argument("--input", required=True, help="CSV file with features and a label column.")
    parser.add_argument("--target", required=True, help="Name of the trusted label column.")
    parser.add_argument("--output-dir", default="ml/results", help="Directory for model and metrics.")
    parser.add_argument("--drop-columns", default="", help="Comma-separated columns to exclude, e.g. sample_id,patient_id.")
    parser.add_argument("--test-size", type=float, default=0.2, help="Held-out test fraction (default: 0.2).")
    parser.add_argument("--seed", type=int, default=42, help="Random seed (default: 42).")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    source = Path(args.input)
    out_dir = Path(args.output_dir)
    if not source.is_file():
        print(f"ERROR: input CSV does not exist: {source}", file=sys.stderr)
        return 2
    if not 0.1 <= args.test_size <= 0.4:
        print("ERROR: --test-size must be between 0.1 and 0.4.", file=sys.stderr)
        return 2

    data = pd.read_csv(source, low_memory=False)
    if args.target not in data.columns:
        print(f"ERROR: target column {args.target!r} not found. Columns: {list(data.columns)}", file=sys.stderr)
        return 2
    data = data.dropna(subset=[args.target]).copy()
    if data.empty:
        print("ERROR: no rows remain after removing missing target labels.", file=sys.stderr)
        return 2
    if data[args.target].nunique() < 2:
        print("ERROR: target must contain at least two classes.", file=sys.stderr)
        return 2

    explicit_drop = {c.strip() for c in args.drop_columns.split(",") if c.strip()}
    absent = explicit_drop - set(data.columns)
    if absent:
        print(f"ERROR: requested drop columns not found: {sorted(absent)}", file=sys.stderr)
        return 2
    auto_drop = {c for c in data.columns if c.strip().lower() in DEFAULT_ID_COLUMNS}
    drop = (explicit_drop | auto_drop | {args.target}) & set(data.columns)
    X = data.drop(columns=sorted(drop)).copy()
    y = data[args.target].astype(str)
    X = X.dropna(axis=1, how="all")
    if X.shape[1] == 0:
        print("ERROR: no usable feature columns remain after exclusions.", file=sys.stderr)
        return 2

    # Remove constant columns; they provide no useful signal and can hide bad input tables.
    constant = [c for c in X.columns if X[c].nunique(dropna=True) <= 1]
    X = X.drop(columns=constant)
    if X.shape[1] == 0:
        print("ERROR: all candidate features are constant or empty.", file=sys.stderr)
        return 2

    counts = y.value_counts()
    if counts.min() < 2:
        print("ERROR: each target class needs at least 2 labelled rows for a stratified holdout split.", file=sys.stderr)
        return 2
    n_test = int(np.ceil(len(y) * args.test_size))
    if n_test < y.nunique() or (len(y) - n_test) < y.nunique():
        print("ERROR: dataset too small for a stratified train/test split with every class represented.", file=sys.stderr)
        return 2

    numeric = X.select_dtypes(include=["number", "bool"]).columns.tolist()
    categorical = [c for c in X.columns if c not in numeric]
    transformers = []
    if numeric:
        transformers.append(("numeric", Pipeline([
            ("imputer", SimpleImputer(strategy="median")),
        ]), numeric))
    if categorical:
        transformers.append(("categorical", Pipeline([
            ("imputer", SimpleImputer(strategy="most_frequent")),
            ("onehot", OneHotEncoder(handle_unknown="ignore")),
        ]), categorical))
    prep = ColumnTransformer(transformers=transformers, remainder="drop")
    model = Pipeline([
        ("preprocess", prep),
        ("classifier", RandomForestClassifier(
            n_estimators=400,
            class_weight="balanced",
            random_state=args.seed,
            n_jobs=-1,
            min_samples_leaf=2,
        )),
    ])

    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=args.test_size, random_state=args.seed, stratify=y
    )
    model.fit(X_train, y_train)
    predicted = model.predict(X_test)
    labels = sorted(y.unique().tolist())
    metrics = {
        "task": "supervised tabular classification",
        "model": "RandomForestClassifier",
        "input_file_name": source.name,
        "target_column": args.target,
        "label_counts": {str(k): int(v) for k, v in counts.items()},
        "rows_used": int(len(data)),
        "train_rows": int(len(X_train)),
        "test_rows": int(len(X_test)),
        "feature_columns_used": list(X.columns),
        "columns_excluded": sorted(drop | set(constant)),
        "numeric_feature_count": len(numeric),
        "categorical_feature_count": len(categorical),
        "random_seed": args.seed,
        "test_size": args.test_size,
        "sklearn_version": sklearn.__version__,
        "metrics": {
            "accuracy": float(accuracy_score(y_test, predicted)),
            "balanced_accuracy": float(balanced_accuracy_score(y_test, predicted)),
            "macro_f1": float(f1_score(y_test, predicted, average="macro", zero_division=0)),
            "macro_precision": float(precision_score(y_test, predicted, average="macro", zero_division=0)),
            "macro_recall": float(recall_score(y_test, predicted, average="macro", zero_division=0)),
        },
        "classification_report": classification_report(
            y_test, predicted, labels=labels, output_dict=True, zero_division=0
        ),
        "confusion_matrix": {
            "labels": labels,
            "matrix": confusion_matrix(y_test, predicted, labels=labels).tolist(),
        },
        "caveats": [
            "Exploratory baseline only; not clinically validated.",
            "Performance depends on label quality, cohort design, and absence of leakage.",
            "Random row splitting can leak information when rows from the same patient/sample are related; use group-aware splits for grouped data.",
            "Do not interpret these metrics as clinical sensitivity or specificity without a suitable validated study design.",
        ],
    }
    if len(labels) == 2 and hasattr(model, "predict_proba"):
        try:
            probs = model.predict_proba(X_test)[:, list(model.classes_).index(labels[1])]
            metrics["metrics"]["roc_auc_for_second_sorted_label"] = float(roc_auc_score((y_test == labels[1]).astype(int), probs))
            metrics["auc_positive_label"] = labels[1]
        except (ValueError, IndexError):
            pass

    out_dir.mkdir(parents=True, exist_ok=True)
    joblib.dump(model, out_dir / "baseline_model.joblib")
    with (out_dir / "metrics.json").open("w", encoding="utf-8") as handle:
        json.dump(metrics, handle, indent=2, allow_nan=False)
    predictions = pd.DataFrame({
        "row_index_from_input": X_test.index,
        "observed_label": y_test.values,
        "predicted_label": predicted,
    })
    predictions.to_csv(out_dir / "heldout_predictions.csv", index=False)
    print(f"Model saved: {out_dir / 'baseline_model.joblib'}")
    print(f"Metrics saved: {out_dir / 'metrics.json'}")
    print(f"Held-out predictions saved: {out_dir / 'heldout_predictions.csv'}")
    print(json.dumps(metrics["metrics"], indent=2))
    print("WARNING: research baseline only. Audit leakage and cohort design before reporting performance.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
