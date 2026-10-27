#!/usr/bin/env python

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr


MACRO = ["SPP1-state Macrophage", "C1QC-state Macrophage", "Other Macrophage"]
STATES = [
    "B cell", "C1QC-state Macrophage", "Dendritic cell", "Endothelial",
    "Epithelial-Normal", "Malignant", "Mast cell", "Monocyte", "NK cell",
    "Other Macrophage", "SMC/Pericyte", "SPP1-state Macrophage", "T cell",
]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def rho(x: pd.Series, y: pd.Series) -> float:
    if x.nunique() < 2 or y.nunique() < 2:
        return float("nan")
    return float(spearmanr(x.to_numpy(dtype=float), y.to_numpy(dtype=float)).statistic)


def close(a: float, b: float, tolerance: float = 1e-10) -> bool:
    if math.isnan(a) and math.isnan(b):
        return True
    return math.isclose(a, b, rel_tol=tolerance, abs_tol=tolerance)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--phase-root", required=True, type=Path)
    args = parser.parse_args()
    root = args.phase_root.resolve()
    results = root / "results_v1"
    output = root / "validation_v1"
    part = root / "validation_v1.part"
    if not results.is_dir():
        raise RuntimeError("Missing results_v1")
    if output.exists() or part.exists():
        raise RuntimeError("Refusing existing validation output")
    part.mkdir(parents=True)

    decision = json.loads((results / "PHASE86_DECISION_v1.json").read_text(encoding="utf-8"))
    counts = pd.read_csv(results / "phase86_mixture_marker_counts_v1.tsv.gz", sep="\t")
    truth = pd.read_csv(results / "phase86_mixture_truth_v1.tsv.gz", sep="\t")
    predictions = pd.read_csv(results / "phase86_mixture_predictions_v1.tsv.gz", sep="\t")
    metrics = pd.read_csv(results / "phase86_state_metrics_v1.tsv", sep="\t")
    accuracy = pd.read_csv(results / "phase86_dominant_accuracy_v1.tsv", sep="\t")
    total = pd.read_csv(results / "phase86_total_macrophage_metrics_v1.tsv", sep="\t")
    gates = pd.read_csv(results / "phase86_gate_decisions_v1.tsv", sep="\t")
    manifest = pd.read_csv(results / "FILE_MANIFEST_SHA256_v1.tsv", sep="\t")

    checks: list[dict] = []

    def record(name: str, value, passed: bool) -> None:
        checks.append({"check": name, "value": value, "pass": bool(passed)})

    record("decision_computation_pass", decision.get("status"), decision.get("status") == "PASS_COMPUTATION_COMPLETE")
    record("datasets_exact", ";".join(sorted(counts["dataset"].unique())), set(counts["dataset"]) == {"GSE207493", "GSE264586", "GSE294109"})
    record("tiers_exact", ";".join(sorted(counts["tier"].unique())), set(counts["tier"]) == {"macrophage_only", "tissue_context"})
    record("methods_exact", ";".join(sorted(predictions["method"].unique())), set(predictions["method"]) == {"frozen13", "frozen13_anchor_excluded", "lodo13", "lodo3_anchor_excluded"})
    record("counts_truth_rows_equal", len(counts), len(counts) == len(truth))
    record("prediction_rows_four_per_mixture", len(predictions), len(predictions) == 4 * len(counts))
    record("mixture_ids_unique_counts", counts["mixture_id"].nunique(), counts["mixture_id"].nunique() == len(counts))
    record("mixture_ids_unique_truth", truth["mixture_id"].nunique(), truth["mixture_id"].nunique() == len(truth))
    record("mixture_id_sets_match", len(counts), set(counts["mixture_id"]) == set(truth["mixture_id"]) == set(predictions["mixture_id"]))

    gene_columns = [column for column in counts.columns if column not in {"mixture_id", "dataset", "patient_id", "tier", "mixture_index"}]
    count_matrix = counts[gene_columns].to_numpy(dtype=float)
    record("marker_gene_count_478", len(gene_columns), len(gene_columns) == 478)
    record("mixture_counts_finite_nonnegative_integer", count_matrix.size, bool(np.isfinite(count_matrix).all() and (count_matrix >= 0).all() and np.allclose(count_matrix, np.round(count_matrix))))
    count_sums = count_matrix.sum(axis=1)
    record("mixture_read_depth_exact", f"{count_sums.min()}-{count_sums.max()}", bool(np.all(count_sums == 500000)))

    full_truth_columns = [f"true_{state}" for state in STATES]
    within_truth_columns = [f"true_within_{state}" for state in MACRO]
    record("truth_full_fractions_sum_one", float(np.max(np.abs(truth[full_truth_columns].sum(axis=1) - 1))), bool(np.allclose(truth[full_truth_columns].sum(axis=1), 1, atol=1e-12)))
    record("truth_within_fractions_sum_one", float(np.max(np.abs(truth[within_truth_columns].sum(axis=1) - 1))), bool(np.allclose(truth[within_truth_columns].sum(axis=1), 1, atol=1e-12)))

    estimate_columns = [f"estimated_within_{state}" for state in MACRO]
    estimate_matrix = predictions[estimate_columns + ["estimated_total_macrophage"]].to_numpy(dtype=float)
    within_sum = predictions[estimate_columns].sum(axis=1).to_numpy(dtype=float)
    record("estimates_finite_nonnegative", estimate_matrix.size, bool(np.isfinite(estimate_matrix).all() and (estimate_matrix >= -1e-12).all()))
    record("estimated_total_bounded", f"{predictions['estimated_total_macrophage'].min()}-{predictions['estimated_total_macrophage'].max()}", bool(((predictions["estimated_total_macrophage"] >= -1e-12) & (predictions["estimated_total_macrophage"] <= 1 + 1e-12)).all()))
    record("estimated_within_sum_zero_or_one", float(np.max(np.minimum(np.abs(within_sum), np.abs(within_sum - 1)))), bool(np.all((np.abs(within_sum) <= 1e-12) | (np.abs(within_sum - 1) <= 1e-10))))

    metric_ok = True
    for row in metrics.itertuples(index=False):
        frame = predictions.loc[
            (predictions["dataset"] == row.dataset)
            & (predictions["tier"] == row.tier)
            & (predictions["method"] == row.method)
        ]
        truth_col = frame[f"true_within_{row.state}"]
        estimate_col = frame[f"estimated_within_{row.state}"]
        difference = estimate_col.to_numpy(dtype=float) - truth_col.to_numpy(dtype=float)
        values = [
            rho(truth_col, estimate_col),
            float(np.mean(np.abs(difference))),
            float(np.sqrt(np.mean(np.square(difference)))),
            float(np.mean(estimate_col.to_numpy(dtype=float) <= 1e-12)),
        ]
        stored = [row.spearman_rho, row.mae, row.rmse, row.estimated_zero_fraction]
        metric_ok = metric_ok and all(close(float(a), float(b)) for a, b in zip(values, stored))
    record("state_metrics_recomputed", len(metrics), metric_ok)

    accuracy_ok = True
    for row in accuracy.itertuples(index=False):
        frame = predictions.loc[
            (predictions["dataset"] == row.dataset)
            & (predictions["tier"] == row.tier)
            & (predictions["method"] == row.method)
        ]
        true_matrix = frame[[f"true_within_{state}" for state in MACRO]].to_numpy(dtype=float)
        est_matrix = frame[[f"estimated_within_{state}" for state in MACRO]].to_numpy(dtype=float)
        value = float(np.mean(np.argmax(true_matrix, axis=1) == np.argmax(est_matrix, axis=1)))
        accuracy_ok = accuracy_ok and close(value, float(row.dominant_state_accuracy))
    record("dominant_accuracy_recomputed", len(accuracy), accuracy_ok)

    total_ok = True
    for row in total.itertuples(index=False):
        frame = predictions.loc[
            (predictions["dataset"] == row.dataset)
            & (predictions["tier"] == row.tier)
            & (predictions["method"] == row.method)
        ]
        values = [
            rho(frame["true_total_macrophage"], frame["estimated_total_macrophage"]),
            float(np.mean(np.abs(frame["estimated_total_macrophage"] - frame["true_total_macrophage"]))),
        ]
        total_ok = total_ok and close(values[0], float(row.spearman_rho)) and close(values[1], float(row.mae))
    record("total_macrophage_metrics_recomputed", len(total), total_ok)

    gate_ok = True
    thresholds = {
        "frozen13": dict(rho=0.70, mae=0.15, zero=0.20, accuracy=0.80, total_rho=0.70, total_mae=0.10),
        "lodo13": dict(rho=0.60, mae=0.18, zero=1.00, accuracy=0.75, total_rho=0.60, total_mae=0.12),
        "frozen13_anchor_excluded": dict(rho=0.50, mae=0.20, zero=0.50, accuracy=0.70),
        "lodo3_anchor_excluded": dict(rho=0.50, mae=0.20, zero=0.50, accuracy=0.70),
    }
    for row in gates.itertuples(index=False):
        threshold = thresholds[row.method]
        frame = metrics.loc[(metrics["method"] == row.method) & (metrics["dataset"] == row.dataset) & (metrics["tier"] == row.tier)]
        acc = float(accuracy.loc[(accuracy["method"] == row.method) & (accuracy["dataset"] == row.dataset) & (accuracy["tier"] == row.tier), "dominant_state_accuracy"].iloc[0])
        state_pass = bool((frame["spearman_rho"] >= threshold["rho"]).all() and (frame["mae"] <= threshold["mae"]).all() and (frame["estimated_zero_fraction"] <= threshold["zero"]).all())
        accuracy_pass = bool(acc >= threshold["accuracy"])
        total_pass = True
        if row.method in {"frozen13", "lodo13"} and row.tier == "tissue_context":
            current = total.loc[(total["method"] == row.method) & (total["dataset"] == row.dataset) & (total["tier"] == row.tier)].iloc[0]
            total_pass = bool(current["spearman_rho"] >= threshold["total_rho"] and current["mae"] <= threshold["total_mae"])
        gate_ok = gate_ok and bool(row.state_gate_pass) == state_pass and bool(row.accuracy_gate_pass) == accuracy_pass and bool(row.total_gate_pass) == total_pass and bool(row.group_gate_pass) == bool(state_pass and accuracy_pass and total_pass)
    record("gate_decisions_recomputed", len(gates), gate_ok)

    operational = bool(gates.loc[gates["method"] == "frozen13", "group_gate_pass"].all())
    generalization = bool(gates.loc[gates["method"] == "lodo13", "group_gate_pass"].all())
    anchor = bool(gates.loc[gates["method"].isin(["frozen13_anchor_excluded", "lodo3_anchor_excluded"]), "group_gate_pass"].all())
    record("decision_operational_recomputed", operational, operational == bool(decision["operational_gate_pass"]))
    record("decision_generalization_recomputed", generalization, generalization == bool(decision["generalization_gate_pass"]))
    record("decision_anchor_recomputed", anchor, anchor == bool(decision["anchor_gate_pass"]))
    expected_decision = "IDENTIFIABLE_FOR_FINE_STATE_BULK_USE" if operational and generalization and anchor else "NOT_IDENTIFIABLE_FOR_FINE_STATE_BULK_USE"
    record("scientific_decision_recomputed", expected_decision, decision["decision"] == expected_decision)

    record("preregistration_hash", decision["preregistration_sha256"], sha256(root / "PHASE86_PREREGISTRATION_v1.md") == decision["preregistration_sha256"])
    record("profile_hash", decision["profile_sha256"], sha256(root / "profiles_v1" / "patient_state_profiles_v1.tsv.gz") == decision["profile_sha256"])
    record("frozen_reference_hash", decision["frozen_reference_sha256"], sha256(root / "inputs_v1" / "phase83_reference_v1" / "phase83_reference_specificity_v1.tsv.gz") == decision["frozen_reference_sha256"])
    record("marker_coverage_hash", decision["marker_coverage_sha256"], sha256(root / "inputs_v1" / "phase83_reference_v1" / "phase83_marker_coverage_v1.tsv") == decision["marker_coverage_sha256"])

    manifest_ok = True
    for row in manifest.itertuples(index=False):
        path = results / row.file
        manifest_ok = manifest_ok and path.is_file() and path.stat().st_size == int(row.bytes) and sha256(path) == row.sha256
    record("result_manifest_exact", len(manifest), manifest_ok)

    forbidden = {"response", "pfs", "os", "treatment", "benefit"}
    all_columns = {column.lower() for frame in [counts, truth, predictions, metrics, accuracy, total, gates] for column in frame.columns}
    record("no_clinical_outcome_columns", ";".join(sorted(all_columns & forbidden)), not bool(all_columns & forbidden))

    check_frame = pd.DataFrame(checks)
    passed = bool(check_frame["pass"].all())
    check_frame.to_csv(part / "PHASE86_VALIDATION_CHECKS_v1.tsv", sep="\t", index=False)
    summary = {
        "status": "PASS" if passed else "FAIL",
        "checks": int(len(check_frame)),
        "passed": int(check_frame["pass"].sum()),
        "failed": int((~check_frame["pass"]).sum()),
        "decision": decision["decision"],
    }
    (part / "PHASE86_VALIDATION_STATUS_v1.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    if not passed:
        raise RuntimeError(f"Phase86 validation failed: {summary}")
    (part / "PHASE86_INDEPENDENT_VALIDATION_PASS_v1").write_text("PASS\n", encoding="utf-8")
    os.replace(part, output)
    print(json.dumps(summary))


if __name__ == "__main__":
    main()
