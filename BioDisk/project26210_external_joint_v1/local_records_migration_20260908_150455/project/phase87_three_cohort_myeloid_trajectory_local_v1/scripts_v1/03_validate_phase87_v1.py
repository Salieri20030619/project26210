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
from scipy.stats import spearmanr, wilcoxon


DATASETS = {"GSE207493", "GSE264586", "GSE294109"}
METHODS = {"slingshot", "geodesic"}
PRIMARY = {
    "P01_SPP1_minus_Monocyte_pseudotime": "positive",
    "P02_C1QC_minus_Monocyte_pseudotime": "positive",
    "P03_Other_minus_Monocyte_pseudotime": "positive",
    "P04_SPP1_anchor_excluded_marker_trend": "positive",
    "P05_C1QC_anchor_excluded_marker_trend": "positive",
    "P06_Other_anchor_excluded_marker_trend": "positive",
    "P07_Monocyte_program_trend": "negative",
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def bh(values: pd.Series) -> pd.Series:
    result = pd.Series(np.nan, index=values.index, dtype=float)
    valid = values.notna()
    p = values.loc[valid].to_numpy(dtype=float)
    if not len(p):
        return result
    order = np.argsort(p)
    ranked = p[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    restored = np.empty_like(adjusted)
    restored[order] = np.clip(adjusted, 0, 1)
    result.loc[valid] = restored
    return result


def signed_rank(values: np.ndarray) -> tuple[float, float]:
    values = values[np.isfinite(values)]
    values = values[np.abs(values) > 1e-15]
    if not len(values):
        return 0.0, 1.0
    test = wilcoxon(values, zero_method="wilcox", correction=False, alternative="two-sided", method="auto")
    return float(test.statistic), float(test.pvalue)


def close(a: float, b: float, tolerance: float = 1e-9) -> bool:
    if math.isnan(a) and math.isnan(b):
        return True
    return math.isclose(a, b, rel_tol=tolerance, abs_tol=tolerance)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--phase-root", required=True, type=Path)
    args = parser.parse_args()
    root = args.phase_root.resolve()
    inputs = root / "trajectory_inputs_v1"
    results = root / "results_v1"
    output = root / "validation_v1"
    part = root / "validation_v1.part"
    if not inputs.is_dir() or not results.is_dir() or output.exists() or part.exists():
        raise RuntimeError("Invalid Phase87 validation state")
    part.mkdir(parents=True)

    decision = json.loads((results / "PHASE87_DECISION_v1.json").read_text(encoding="utf-8"))
    cells = pd.read_csv(results / "phase87_cell_pseudotime_v1.tsv.gz", sep="\t")
    units = pd.read_csv(results / "phase87_patient_sample_effects_v1.tsv", sep="\t")
    patients = pd.read_csv(results / "phase87_patient_effects_v1.tsv", sep="\t")
    summary = pd.read_csv(results / "phase87_dataset_method_endpoint_summary_v1.tsv", sep="\t")
    unit_concordance = pd.read_csv(results / "phase87_unit_method_concordance_v1.tsv", sep="\t")
    patient_concordance = pd.read_csv(results / "phase87_patient_method_concordance_v1.tsv", sep="\t")
    concordance = pd.read_csv(results / "phase87_dataset_method_concordance_summary_v1.tsv", sep="\t")
    graph_qc = pd.read_csv(results / "phase87_geodesic_graph_qc_v1.tsv", sep="\t")
    manifest = pd.read_csv(results / "FILE_MANIFEST_SHA256_v1.tsv", sep="\t")
    input_status = pd.read_csv(inputs / "phase87_dataset_status_v1.tsv", sep="\t")
    root_qc = pd.read_csv(inputs / "phase87_root_qc_v1.tsv", sep="\t")
    program_coverage = pd.read_csv(inputs / "phase87_program_coverage_v1.tsv", sep="\t")

    checks: list[dict] = []

    def record(name: str, value, passed: bool) -> None:
        checks.append({"check": name, "value": value, "pass": bool(passed)})

    record("decision_computation_pass", decision.get("status"), decision.get("status") == "PASS_COMPUTATION_COMPLETE")
    record("datasets_exact_cells", ";".join(sorted(cells["dataset"].unique())), set(cells["dataset"]) == DATASETS)
    record("datasets_exact_status", ";".join(sorted(input_status["dataset"].unique())), set(input_status["dataset"]) == DATASETS and len(input_status) == 3)
    record("input_status_all_pass", len(input_status), input_status["status"].eq("PASS_SLINGSHOT_EXTRACTION").all())
    record("cell_id_unique", cells["cell_id"].nunique(), cells["cell_id"].nunique() == len(cells))
    record("trajectory_cell_count_conserved", len(cells), len(cells) == int(input_status["trajectory_cells"].sum()) == int(decision["trajectory_cells"]))
    record("states_exact", ";".join(sorted(cells["phase8_state"].unique())), set(cells["phase8_state"]) == {"Monocyte", "SPP1-state Macrophage", "C1QC-state Macrophage", "Other Macrophage"})
    record("one_root_cell_per_dataset", int(cells["root_cell"].sum()), cells.groupby("dataset")["root_cell"].sum().eq(1).all())
    record("root_cells_match_qc", len(root_qc), set(cells.loc[cells["root_cell"], "cell_id"]) == set(root_qc["root_cell"]))
    record("program_coverage_all_at_least_70pct", program_coverage["coverage_fraction"].min(), (program_coverage["coverage_fraction"] >= 0.70).all() and len(program_coverage) == 15)

    for method in METHODS:
        values = cells[f"{method}_pseudotime"].to_numpy(dtype=float)
        finite = np.isfinite(values)
        record(f"{method}_finite_at_least_80pct", float(finite.mean()), finite.mean() >= 0.80)
        record(f"{method}_bounded_zero_one", f"{np.nanmin(values)}-{np.nanmax(values)}", bool((values[finite] >= -1e-12).all() and (values[finite] <= 1 + 1e-12).all()))
    record("geodesic_graph_qc_three_datasets", len(graph_qc), len(graph_qc) == 3 and set(graph_qc["dataset"]) == DATASETS and (graph_qc["finite_cells"] >= 0.80 * graph_qc["cells"]).all())

    expected_endpoints = set(PRIMARY) | {"S01_Immunosuppression_no_IL10_trend"}
    record("unit_endpoint_set_exact", ";".join(sorted(units["endpoint_id"].unique())), set(units["endpoint_id"]) == expected_endpoints)
    record("unit_method_set_exact", ";".join(sorted(units["method"].unique())), set(units["method"]) == METHODS)
    record("patient_keys_unique", len(patients), not patients.duplicated(["dataset", "patient_id", "method", "endpoint_id"]).any())
    record("summary_grid_complete", len(summary), len(summary) == 3 * 2 * 8 and not summary.duplicated(["dataset", "method", "endpoint_id"]).any())

    patient_recomputed = (
        units.dropna(subset=["effect"])
        .groupby(["dataset", "patient_id", "method", "endpoint_id", "tier", "expected_direction"], sort=True)
        .agg(patient_median_effect=("effect", "median"), n_samples=("sample_id", "nunique"), total_cells=("n_cells", "sum"))
        .reset_index()
    )
    patient_columns = ["dataset", "patient_id", "method", "endpoint_id", "tier", "expected_direction", "patient_median_effect", "n_samples", "total_cells"]
    stored_patient = patients[patient_columns].sort_values(patient_columns[:6]).reset_index(drop=True)
    recomputed_patient = patient_recomputed[patient_columns].sort_values(patient_columns[:6]).reset_index(drop=True)
    patient_keys_equal = stored_patient[patient_columns[:6]].astype(str).equals(recomputed_patient[patient_columns[:6]].astype(str))
    patient_numeric_equal = (
        np.allclose(stored_patient["patient_median_effect"], recomputed_patient["patient_median_effect"], rtol=1e-9, atol=1e-9, equal_nan=True)
        and np.array_equal(stored_patient["n_samples"].to_numpy(dtype=int), recomputed_patient["n_samples"].to_numpy(dtype=int))
        and np.array_equal(stored_patient["total_cells"].to_numpy(dtype=int), recomputed_patient["total_cells"].to_numpy(dtype=int))
    )
    record("patient_effects_recomputed", len(stored_patient), patient_keys_equal and patient_numeric_equal)

    summary_ok = True
    for row in summary.itertuples(index=False):
        frame = patients.loc[
            (patients["dataset"] == row.dataset)
            & (patients["method"] == row.method)
            & (patients["endpoint_id"] == row.endpoint_id)
        ]
        values = frame["patient_median_effect"].to_numpy(dtype=float)
        statistic, p_value = signed_rank(values)
        summary_ok = summary_ok and row.n_patients == len(values)
        summary_ok = summary_ok and close(float(row.median_patient_effect), float(np.median(values)))
        summary_ok = summary_ok and close(float(row.wilcoxon_statistic), statistic)
        summary_ok = summary_ok and close(float(row.p_value), p_value)
    record("summary_effects_and_p_recomputed", len(summary), summary_ok)

    q_ok = True
    support_ok = True
    for (dataset, method), frame in summary.loc[summary["tier"] == "PRIMARY"].groupby(["dataset", "method"], sort=True):
        recalculated = bh(frame["p_value"])
        q_ok = q_ok and np.allclose(recalculated.to_numpy(dtype=float), frame["BH_q_within_dataset_method_primary7"].to_numpy(dtype=float), equal_nan=True, rtol=1e-9, atol=1e-9)
        for index, row in frame.iterrows():
            direction_ok = (row["expected_direction"] == "positive" and row["median_patient_effect"] > 0) or (row["expected_direction"] == "negative" and row["median_patient_effect"] < 0)
            expected = bool(row["n_patients"] >= 10 and row["BH_q_within_dataset_method_primary7"] < 0.05 and direction_ok)
            support_ok = support_ok and bool(row["formal_support"]) == expected
    record("primary_bh_recomputed", 6, q_ok)
    record("formal_support_recomputed", 42, support_ok)

    concordance_ok = True
    recomputed_units = []
    for (dataset, patient, sample_id), frame in cells.groupby(["dataset", "patient_id", "sample_id"], sort=True):
        x = frame["slingshot_pseudotime"].to_numpy(dtype=float)
        y = frame["geodesic_pseudotime"].to_numpy(dtype=float)
        ok = np.isfinite(x) & np.isfinite(y)
        value = float(spearmanr(x[ok], y[ok]).statistic) if ok.sum() >= 3 and np.unique(x[ok]).size >= 3 and np.unique(y[ok]).size >= 3 else float("nan")
        recomputed_units.append((dataset, patient, sample_id, value))
    recomputed_unit = pd.DataFrame(recomputed_units, columns=["dataset", "patient_id", "sample_id", "pseudotime_rho"])
    merged_unit = unit_concordance.merge(recomputed_unit, on=["dataset", "patient_id", "sample_id"], suffixes=("_stored", "_recomputed"))
    concordance_ok = len(merged_unit) == len(unit_concordance) and all(close(float(a), float(b)) for a, b in zip(merged_unit["pseudotime_rho_stored"], merged_unit["pseudotime_rho_recomputed"]))
    record("unit_method_concordance_recomputed", len(unit_concordance), concordance_ok)

    dataset_concordance_ok = True
    for row in concordance.itertuples(index=False):
        frame = patient_concordance.loc[patient_concordance["dataset"] == row.dataset]
        values = frame["patient_median_rho"].to_numpy(dtype=float)
        median = float(np.median(values))
        positive_fraction = float(np.mean(values > 0))
        gate = bool(median > 0 and positive_fraction >= 0.70)
        dataset_concordance_ok = dataset_concordance_ok and row.n_patients == len(values) and close(float(row.median_patient_rho), median) and close(float(row.positive_fraction), positive_fraction) and bool(row.concordance_gate_pass) == gate
    record("dataset_method_concordance_recomputed", len(concordance), dataset_concordance_ok)

    def support(dataset: str, method: str, endpoint: str) -> bool:
        row = summary.loc[(summary["dataset"] == dataset) & (summary["method"] == method) & (summary["endpoint_id"] == endpoint)]
        return len(row) == 1 and bool(row["formal_support"].iloc[0])

    p01_all = all(support(dataset, method, "P01_SPP1_minus_Monocyte_pseudotime") for dataset in DATASETS for method in METHODS)
    p04_external = all(support(dataset, method, "P04_SPP1_anchor_excluded_marker_trend") for dataset in {"GSE264586", "GSE294109"} for method in METHODS)
    concordance_gate = bool(len(concordance) == 3 and concordance["concordance_gate_pass"].all())
    if p01_all and p04_external and concordance_gate:
        expected_decision = "THREE_COHORT_MONOCYTE_TO_SPP1_ORDERING_REPLICATED"
    elif p01_all and concordance_gate:
        expected_decision = "STATE_ORDERING_WITHOUT_ANCHOR_EXCLUDED_PROGRAM_REPLICATION"
    else:
        expected_decision = "MONOCYTE_TO_SPP1_ORDERING_NOT_THREE_COHORT_REPLICATED"
    record("decision_p01_recomputed", p01_all, p01_all == bool(decision["p01_all_three_dual_method"]))
    record("decision_p04_recomputed", p04_external, p04_external == bool(decision["p04_two_external_dual_method"]))
    record("decision_concordance_recomputed", concordance_gate, concordance_gate == bool(decision["method_concordance_gate"]))
    record("scientific_decision_recomputed", expected_decision, decision["decision"] == expected_decision)

    record("preregistration_hash", decision["preregistration_sha256"], sha256(root / "PHASE87_PREREGISTRATION_v1.md") == decision["preregistration_sha256"])
    record("marker_table_hash", decision["marker_table_sha256"], sha256(root / "inputs_v1" / "phase83_marker_members_v1.tsv") == decision["marker_table_sha256"])
    manifest_ok = True
    for row in manifest.itertuples(index=False):
        path = results / row.file
        manifest_ok = manifest_ok and path.is_file() and path.stat().st_size == int(row.bytes) and sha256(path) == row.sha256
    record("result_manifest_exact", len(manifest), manifest_ok)

    forbidden = {"response", "pfs", "os", "treatment", "benefit"}
    columns = {column.lower() for frame in [cells, units, patients, summary, unit_concordance, patient_concordance, concordance] for column in frame.columns}
    record("no_clinical_outcome_columns", ";".join(sorted(columns & forbidden)), not bool(columns & forbidden))

    check_frame = pd.DataFrame(checks)
    passed = bool(check_frame["pass"].all())
    check_frame.to_csv(part / "PHASE87_VALIDATION_CHECKS_v1.tsv", sep="\t", index=False)
    status = {
        "status": "PASS" if passed else "FAIL",
        "checks": int(len(check_frame)),
        "passed": int(check_frame["pass"].sum()),
        "failed": int((~check_frame["pass"]).sum()),
        "decision": decision["decision"],
    }
    (part / "PHASE87_VALIDATION_STATUS_v1.json").write_text(json.dumps(status, indent=2) + "\n", encoding="utf-8")
    if not passed:
        raise RuntimeError(f"Phase87 validation failed: {check_frame.loc[~check_frame['pass']].to_dict('records')}")
    (part / "PHASE87_INDEPENDENT_VALIDATION_PASS_v1").write_text("PASS\n", encoding="utf-8")
    os.replace(part, output)
    print(json.dumps(status))


if __name__ == "__main__":
    main()
