#!/usr/bin/env python3
# Random forest trainer modeled after scripts/random_forest_cmdline.R (caret + retrain)
# Usage:
#   python3 scripts/random_forest_cmdline.py \
#     --truth-file /path/to/all.truth_annotated \
#     --output-prefix /path/to/output/random_forest_caret \
#     --ntree 512 \
#     --maxnodes 0 \
#     --test-fraction 0.2
# Optional: --seed 1234 --threads 4

from __future__ import annotations

import argparse
import math
import sys
from typing import Iterable, Optional, Tuple
import xml.etree.ElementTree as ET

import joblib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import confusion_matrix, precision_recall_curve, average_precision_score
from sklearn.model_selection import GridSearchCV

TRUE_LABEL = "TRUE"
FALSE_LABEL = "FALSE"

INF_COLS = [
    "lu_gene_rate",
    "lsu_gene_rate",
    "lu_gene_rate2",
    "lsu_gene_rate2",
    "lu_gene_rate3",
    "lsu_gene_rate3",
    "lu_rate",
    "lsu_rate",
    "su_rate",
    "lsu_per_read",
    "lu_per_read",
]


def _as_bool_series(values: pd.Series) -> pd.Series:
    if values.dtype == bool:
        return values
    if np.issubdtype(values.dtype, np.number):
        return values.fillna(0).astype(float) > 0
    lowered = values.fillna("").astype(str).str.lower()
    truthy = {"true", "t", "1", "yes", "y"}
    return lowered.isin(truthy)


def sensitivity(cm: np.ndarray) -> float:
    tp = cm[1, 1]
    fn = cm[0, 1]
    return tp / (tp + fn) if (tp + fn) > 0 else 0.0


def precision(cm: np.ndarray) -> float:
    tp = cm[1, 1]
    fp = cm[1, 0]
    return tp / (tp + fp) if (tp + fp) > 0 else 0.0


def f1(cm: np.ndarray) -> float:
    tp = cm[1, 1]
    fp = cm[1, 0]
    fn = cm[0, 1]
    denom = (2 * tp + fp + fn)
    return (2 * tp) / denom if denom > 0 else 0.0


def evaluate_split(
    model: RandomForestClassifier,
    test_data: pd.DataFrame,
    train_data: pd.DataFrame,
    label_col: str = "truth",
) -> Tuple[dict, dict]:
    pred_test = model.predict(test_data.drop(columns=[label_col]))
    pred_train = model.predict(train_data.drop(columns=[label_col]))
    ct_test = confusion_matrix(
        test_data[label_col], pred_test, labels=[FALSE_LABEL, TRUE_LABEL]
    )
    ct_train = confusion_matrix(
        train_data[label_col], pred_train, labels=[FALSE_LABEL, TRUE_LABEL]
    )
    return (
        {
            "ct": ct_test,
            "sensitivity": sensitivity(ct_test),
            "precision": precision(ct_test),
            "f1": f1(ct_test),
        },
        {
            "ct": ct_train,
            "sensitivity": sensitivity(ct_train),
            "precision": precision(ct_train),
            "f1": f1(ct_train),
        },
    )


def parse_args(argv: Optional[Iterable[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--truth-file", required=True)
    parser.add_argument("--output-prefix", required=True)
    parser.add_argument("--ntree", type=int, default=256)
    parser.add_argument("--maxnodes", type=int, default=128)
    parser.add_argument("--test-fraction", type=float, default=0.2)
    parser.add_argument("--seed", type=int)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--reference-pmml")
    parser.add_argument(
        "--probability",
        action="store_true",
        help="Output prediction probabilities (0–1) instead of TRUE/FALSE labels, "
             "and save a precision-recall curve plot.",
    )
    return parser.parse_args(argv)


def load_truth_data(path: str) -> pd.DataFrame:
    try:
        df = pd.read_csv(path, sep="\t", dtype=str)
    except Exception as exc:
        raise RuntimeError(f"Failed to read truth file: {exc}") from exc
    for col in df.columns:
        # Allow numeric columns to become numeric when possible.
        try:
            df[col] = pd.to_numeric(df[col])
        except (ValueError, TypeError):
            pass
    df["dataset"] = "dataset"
    df["truth_raw"] = _as_bool_series(df["truth"])
    if "prediction" in df.columns:
        df["prediction"] = _as_bool_series(df["prediction"])
    df["truth"] = df["truth_raw"].map({True: TRUE_LABEL, False: FALSE_LABEL})
    for col in INF_COLS:
        if col in df.columns:
            values = pd.to_numeric(df[col], errors="coerce")
            values = values.replace([np.inf, -np.inf], 0).fillna(values)
            df[col] = values
    return df


def load_reference_feature_cols(path: str) -> list[str]:
    tree = ET.parse(path)
    root = tree.getroot()
    fields = []
    for elem in root.iter():
        if elem.tag.endswith("MiningField"):
            usage = elem.attrib.get("usageType", "active")
            name = elem.attrib.get("name")
            if not name:
                continue
            if usage.lower() in {"target", "predicted"}:
                continue
            fields.append(name)
    # Preserve order but deduplicate.
    seen = set()
    ordered = []
    for name in fields:
        if name not in seen:
            ordered.append(name)
            seen.add(name)
    return ordered


def load_reference_data_fields(path: str) -> list[str]:
    tree = ET.parse(path)
    root = tree.getroot()
    fields = []
    for elem in root.iter():
        if elem.tag.endswith("DataField"):
            name = elem.attrib.get("name")
            if not name or name == "truth":
                continue
            fields.append(name)
    return fields


def save_pmml(
    model: RandomForestClassifier,
    train_data: pd.DataFrame,
    label_col: str,
    pmml_path: str,
) -> bool:
    import traceback
    try:
        from sklearn2pmml import PMMLPipeline, sklearn2pmml
    except Exception as exc:
        print(f"sklearn2pmml import failed: {exc}")
        traceback.print_exc()
        return False

    pipeline = PMMLPipeline([("classifier", model)])
    pipeline.fit(train_data.drop(columns=[label_col]), train_data[label_col])
    try:
        sklearn2pmml(pipeline, pmml_path, with_repr=True)
    except Exception as exc:
        print(f"PMML export failed: {exc}")
        traceback.print_exc()
        return False
    if not normalize_pmml_double_casts(pmml_path):
        print("PMML cleanup skipped (no double() casts detected).")
    return True


def normalize_pmml_double_casts(pmml_path: str) -> bool:
    tree = ET.parse(pmml_path)
    root = tree.getroot()
    strip_pmml_namespaces(root)
    mapping = {}

    for parent in root.iter():
        for child in list(parent):
            if not child.tag.endswith("DerivedField"):
                continue
            name = child.attrib.get("name", "")
            if name.startswith("double(") and name.endswith(")"):
                base = name[len("double(") : -1]
                mapping[name] = base
                parent.remove(child)

    if not mapping:
        return False

    for elem in root.iter():
        for attr in ("name", "field"):
            if attr in elem.attrib and elem.attrib[attr] in mapping:
                elem.attrib[attr] = mapping[elem.attrib[attr]]

    tree.write(pmml_path, encoding="utf-8", xml_declaration=True)
    return True


def strip_pmml_namespaces(root: ET.Element) -> None:
    # cPMML expects unqualified element names like "PMML".
    for elem in root.iter():
        if "}" in elem.tag:
            elem.tag = elem.tag.split("}", 1)[1]
    # Remove namespace declarations left on the root, if any.
    for attr in list(root.attrib):
        if attr.startswith("xmlns"):
            del root.attrib[attr]


def find_best_threshold(
    model: RandomForestClassifier,
    test_data: pd.DataFrame,
    label_col: str,
) -> Tuple[float, float, float, float]:
    """Return (threshold, precision, recall, f1) that maximises F1 on test_data."""
    true_idx = list(model.classes_).index(TRUE_LABEL)
    scores = model.predict_proba(test_data.drop(columns=[label_col]))[:, true_idx]
    y_true = (test_data[label_col] == TRUE_LABEL).astype(int)
    prec, rec, thresholds = precision_recall_curve(y_true, scores)
    # precision_recall_curve appends a final point at recall=0 with no matching threshold.
    prec_t, rec_t = prec[:-1], rec[:-1]
    denom = prec_t + rec_t
    f1_scores = np.where(denom > 0, 2 * prec_t * rec_t / denom, 0.0)
    best_idx = int(np.argmax(f1_scores))
    return (
        float(thresholds[best_idx]),
        float(prec_t[best_idx]),
        float(rec_t[best_idx]),
        float(f1_scores[best_idx]),
    )


def plot_pr_curve(
    model: RandomForestClassifier,
    test_data: pd.DataFrame,
    label_col: str,
    output_path: str,
    best_threshold: Optional[float] = None,
) -> None:
    true_idx = list(model.classes_).index(TRUE_LABEL)
    scores = model.predict_proba(test_data.drop(columns=[label_col]))[:, true_idx]
    y_true = (test_data[label_col] == TRUE_LABEL).astype(int)
    prec, rec, thresholds = precision_recall_curve(y_true, scores)
    ap = average_precision_score(y_true, scores)

    # Identify the plot point corresponding to best_threshold (if provided).
    best_prec = best_rec = None
    if best_threshold is not None:
        prec_t, rec_t = prec[:-1], rec[:-1]
        diffs = np.abs(thresholds - best_threshold)
        bi = int(np.argmin(diffs))
        best_prec, best_rec = float(prec_t[bi]), float(rec_t[bi])

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Left: precision-recall curve
    axes[0].plot(rec, prec, lw=2, label=f"AP = {ap:.3f}")
    if best_prec is not None:
        axes[0].scatter(
            [best_rec], [best_prec], s=80, zorder=5,
            label=f"best t={best_threshold:.3f}",
        )
    axes[0].set_xlabel("Recall")
    axes[0].set_ylabel("Precision")
    axes[0].set_title("Precision-Recall Curve")
    axes[0].set_xlim([0, 1])
    axes[0].set_ylim([0, 1])
    axes[0].legend()

    # Right: precision and recall vs threshold
    axes[1].plot(thresholds, prec[:-1], lw=2, label="Precision")
    axes[1].plot(thresholds, rec[:-1], lw=2, label="Recall")
    if best_threshold is not None:
        axes[1].axvline(best_threshold, color="grey", linestyle="--", lw=1.5,
                        label=f"best t={best_threshold:.3f}")
    axes[1].set_xlabel("Decision threshold")
    axes[1].set_ylabel("Score")
    axes[1].set_title("Precision / Recall vs Threshold")
    axes[1].set_xlim([0, 1])
    axes[1].set_ylim([0, 1])
    axes[1].legend()

    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)


def main() -> int:
    opts = parse_args()
    rng = np.random.RandomState(opts.seed) if opts.seed is not None else np.random

    print("Load data")
    data = load_truth_data(opts.truth_file)

    if opts.reference_pmml:
        feature_cols = load_reference_data_fields(opts.reference_pmml)
        missing = sorted(set(feature_cols) - set(data.columns))
        if missing:
            raise RuntimeError(
                "Reference PMML fields missing from truth data: "
                + ", ".join(missing)
            )
    else:
        feature_cols = [
            col
            for col in data.columns
            if col
            not in {
                "total_hits",
                "truth",
                "truth_raw",
                "taxon",
                "prediction",
                "probability",
                "dataset",
				"taxon_name"
            }
        ]

    print("train caret on columns")
    print(feature_cols)

    max_features_end = min(25, len(feature_cols))
    if max_features_end < 1:
        raise RuntimeError("No feature columns available for training.")
    max_features_start = 1 if max_features_end < 6 else 6
    if opts.reference_pmml:
        chosen_cols = feature_cols
        chosen_mtry = len(chosen_cols)
        print(f"Reference PMML features: {chosen_mtry}")
    else:
        grid = {"max_features": list(range(max_features_start, max_features_end + 1))}
        base_model = RandomForestClassifier(
            n_estimators=opts.ntree,
            n_jobs=opts.threads,
            random_state=opts.seed,
            max_leaf_nodes=opts.maxnodes if opts.maxnodes > 0 else None,
            class_weight="balanced",
        )
        grid_search = GridSearchCV(
            estimator=base_model,
            param_grid=grid,
            cv=5,
            n_jobs=opts.threads,
        )
        grid_search.fit(data[feature_cols], data["truth"])
        chosen_mtry = int(grid_search.best_params_["max_features"])
        print(f"Mtry: {chosen_mtry}")

        importances = grid_search.best_estimator_.feature_importances_
        vi_df = pd.DataFrame(
            {"feature": feature_cols, "importance": importances}
        ).sort_values("importance", ascending=False)
        chosen_cols = vi_df["feature"].head(chosen_mtry).tolist()

    train_df = data[["truth"] + chosen_cols]
    n = len(train_df)
    test_size = max(1, int(math.floor(opts.test_fraction * n)))
    idx = rng.choice(n, size=test_size, replace=False)
    test_data = train_df.iloc[idx]
    train_data = train_df.drop(train_df.index[idx])

    print("Train forest")
    rf = RandomForestClassifier(
        n_estimators=opts.ntree,
        n_jobs=opts.threads,
        random_state=opts.seed,
        max_leaf_nodes=opts.maxnodes if opts.maxnodes > 0 else None,
        class_weight="balanced",
    )
    rf.fit(train_data.drop(columns=["truth"]), train_data["truth"])

    prefix = opts.output_prefix
    pmml_path = f"{prefix}.xml"
    varimp_path = f"{prefix}.varimp.tsv"
    rds_path = f"{prefix}.rds"
    varimp_png = f"{prefix}.varimp.png"
    pr_curve_png = f"{prefix}.pr_curve.png"

    pmml_saved = save_pmml(rf, train_data, "truth", pmml_path)

    vi_df = pd.DataFrame(
        {"feature": chosen_cols, "importance": rf.feature_importances_}
    ).sort_values("importance", ascending=False)
    vi_df.to_csv(varimp_path, sep="\t", index=False)

    plt.figure(figsize=(12, 8))
    plt.bar(vi_df["feature"], vi_df["importance"])
    plt.title("Variable Importance")
    plt.xticks(rotation=90, fontsize=8)
    plt.tight_layout()
    plt.savefig(varimp_png, dpi=150)
    plt.close()

    joblib.dump(rf, rds_path)

    eval_test, eval_train = evaluate_split(rf, test_data, train_data)
    print(
        f"Split: train={len(train_data)} test={len(test_data)} "
        f"(test fraction={opts.test_fraction:.3f})"
    )
    print(
        f"Train   sens={eval_train['sensitivity']:.4f} "
        f"prec={eval_train['precision']:.4f} "
        f"f1={eval_train['f1']:.4f}"
    )
    print(
        f"Test    sens={eval_test['sensitivity']:.4f} "
        f"prec={eval_test['precision']:.4f} "
        f"f1={eval_test['f1']:.4f}"
    )
    if opts.reference_pmml:
        print(f"Chosen features (reference PMML): {chosen_mtry}")
    else:
        print(f"Chosen mtry (grid search): {chosen_mtry}")
    if not pmml_saved:
        raise RuntimeError("PMML export failed; XML output is required.")
    print(f"Saved PMML: {pmml_path}")
    print(f"Saved varimp: {varimp_path}")
    print(f"Saved varimp plot: {varimp_png}")
    print(f"Saved model (joblib): {rds_path}")

    if opts.probability:
        best_t, best_p, best_r, best_f = find_best_threshold(rf, test_data, "truth")
        print(
            f"Best threshold (max F1 on test set): {best_t:.4f}  "
            f"prec={best_p:.4f}  recall={best_r:.4f}  f1={best_f:.4f}"
        )

        best_thresh_path = f"{prefix}.best_threshold.tsv"
        pd.DataFrame([{
            "threshold": best_t,
            "precision": best_p,
            "recall": best_r,
            "f1": best_f,
        }]).to_csv(best_thresh_path, sep="\t", index=False, float_format="%.6f")
        print(f"Saved best threshold: {best_thresh_path}")

        plot_pr_curve(rf, test_data, "truth", pr_curve_png, best_threshold=best_t)
        print(f"Saved PR curve: {pr_curve_png}")
        print("Probability mode: use model.predict_proba(X)[:, 1] for float scores (0–1).")

        threshold_path = f"{prefix}.thresholds.tsv"
        true_idx = list(rf.classes_).index(TRUE_LABEL)
        scores = rf.predict_proba(test_data.drop(columns=["truth"]))[:, true_idx]
        y_true = (test_data["truth"] == TRUE_LABEL).astype(int).to_numpy()
        thresholds = np.arange(0.0, 1.001, 0.001)
        rows = []
        for t in thresholds:
            pred = (scores >= t).astype(int)
            tp = int(((pred == 1) & (y_true == 1)).sum())
            fp = int(((pred == 1) & (y_true == 0)).sum())
            fn = int(((pred == 0) & (y_true == 1)).sum())
            prec_t = tp / (tp + fp) if (tp + fp) > 0 else 0.0
            sens_t = tp / (tp + fn) if (tp + fn) > 0 else 0.0
            denom = 2 * tp + fp + fn
            f1_t = (2 * tp) / denom if denom > 0 else 0.0
            rows.append({"threshold": round(t, 3), "precision": prec_t, "sensitivity": sens_t, "f1": f1_t})
        thresh_df = pd.DataFrame(rows)
        thresh_df.to_csv(threshold_path, sep="\t", index=False, float_format="%.6f")
        print(f"Saved threshold table: {threshold_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
