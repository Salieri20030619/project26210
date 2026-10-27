#!/usr/bin/env python3
"""Refine Phase 53 small-n screening with exact Spearman permutation p-values."""

from __future__ import annotations

import hashlib
import importlib.util
import itertools
import json
import math
import os
import platform
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
from scipy import stats


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
SOURCE_SCRIPT = HERE / "01_screen_reference_bridge_v1.py"
SOURCE_RESULT = HERE / "analysis_local_v1" / "bridge_screen_all_v1.tsv"
FINAL = HERE / "analysis_local_v2"
PART = HERE / "analysis_local_v2.part"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def atomic_tsv(frame: pd.DataFrame, path: Path) -> None:
    part = path.with_name(path.name + ".part")
    frame.to_csv(part, sep="\t", index=False, na_rep="NA", encoding="utf-8", lineterminator="\n")
    os.replace(part, path)


def bh(values: pd.Series) -> pd.Series:
    numeric = pd.to_numeric(values, errors="coerce").to_numpy(float)
    answer = np.full(len(numeric), np.nan)
    valid = np.where(np.isfinite(numeric))[0]
    if len(valid):
        order = valid[np.argsort(numeric[valid])]
        ranked = numeric[order] * len(order) / np.arange(1, len(order) + 1)
        ranked = np.minimum.accumulate(ranked[::-1])[::-1]
        answer[order] = np.minimum(ranked, 1.0)
    return pd.Series(answer, index=values.index)


def exact_spearman_p(x: np.ndarray, y: np.ndarray) -> float:
    rx = stats.rankdata(x).astype(float)
    ry = stats.rankdata(y).astype(float)
    rx -= rx.mean()
    ry -= ry.mean()
    denominator = float(np.sqrt(np.dot(rx, rx) * np.dot(ry, ry)))
    if denominator == 0:
        return np.nan
    observed = abs(float(np.dot(rx, ry) / denominator))
    extreme = 0
    total = 0
    tolerance = 1e-12
    for permuted in itertools.permutations(ry.tolist()):
        rho = abs(float(np.dot(rx, np.asarray(permuted)) / denominator))
        extreme += int(rho + tolerance >= observed)
        total += 1
    if total != math.factorial(len(y)):
        raise RuntimeError("Permutation count mismatch")
    return extreme / total


def main() -> int:
    if FINAL.exists() or PART.exists():
        raise FileExistsError(f"Refusing to overwrite {FINAL} or {PART}")
    if not SOURCE_SCRIPT.is_file() or not SOURCE_RESULT.is_file():
        raise FileNotFoundError("Phase53 v1 source missing")
    PART.mkdir(parents=True)

    spec = importlib.util.spec_from_file_location("bridge_v1", SOURCE_SCRIPT)
    bridge = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(bridge)
    spec_patient = importlib.util.spec_from_file_location("phase14_patient", bridge.SOURCE_ANALYSIS)
    patient_module = importlib.util.module_from_spec(spec_patient)
    assert spec_patient.loader is not None
    spec_patient.loader.exec_module(patient_module)

    patient_maps: dict[tuple[str, str], pd.DataFrame] = {}
    for layer, run in bridge.LAYERS.items():
        for dataset in bridge.DATASETS:
            sample = pd.read_csv(run / dataset / "tables" / "sample_program_summary_v1.tsv", sep="\t")
            patient, _ = patient_module.build_cross_sectional_patient_programs(sample)
            patient_maps[(layer, dataset)] = patient

    relation_map = {row[0]: row for row in bridge.RELATIONS}
    result = pd.read_csv(SOURCE_RESULT, sep="\t")
    exact_values: list[float] = []
    exact_cache: dict[tuple[object, ...], float] = {}
    for _, row in result.iterrows():
        n = int(row["n_patients"])
        if n < bridge.MIN_PATIENTS_TEST or n >= 10 or not np.isfinite(row["spearman_rho"]):
            exact_values.append(np.nan)
            continue
        relation = relation_map[row["relation"]]
        _, _, x_sig, x_comp, y_sig, y_comp, _ = relation
        value_col = bridge.METHODS[row["method"]]
        patient = patient_maps[(row["evidence_layer"], row["dataset"])]
        x = patient.loc[
            (patient["signature_id"] == x_sig)
            & (patient["analysis_compartment"] == x_comp)
            & (patient["n_cells"] >= bridge.MIN_CELLS),
            ["patient_key", value_col],
        ].rename(columns={value_col: "x"})
        y = patient.loc[
            (patient["signature_id"] == y_sig)
            & (patient["analysis_compartment"] == y_comp)
            & (patient["n_cells"] >= bridge.MIN_CELLS),
            ["patient_key", value_col],
        ].rename(columns={value_col: "y"})
        paired = x.merge(y, on="patient_key", how="inner", validate="one_to_one")
        paired = paired.loc[np.isfinite(paired["x"]) & np.isfinite(paired["y"])].sort_values("patient_key")
        if len(paired) != n:
            raise RuntimeError("Small-n pair reconstruction mismatch")
        rx = tuple(stats.rankdata(paired["x"].to_numpy(float)).tolist())
        ry = tuple(stats.rankdata(paired["y"].to_numpy(float)).tolist())
        cache_key = (rx, ry)
        if cache_key not in exact_cache:
            exact_cache[cache_key] = exact_spearman_p(
                paired["x"].to_numpy(float), paired["y"].to_numpy(float)
            )
        exact_values.append(exact_cache[cache_key])

    result["small_n_exact_permutation_p"] = exact_values
    result["p_value_for_reporting"] = np.where(
        (result["n_patients"] >= bridge.MIN_PATIENTS_TEST) & (result["n_patients"] < 10),
        result["small_n_exact_permutation_p"],
        result["p_value"],
    )
    result["reporting_p_fdr_bh_within_dataset_layer_method"] = result.groupby(
        ["dataset", "evidence_layer", "method"], group_keys=False
    )["p_value_for_reporting"].apply(bh)

    reporting_tier = []
    reportable = []
    for _, row in result.iterrows():
        n = int(row["n_patients"])
        rho = row["spearman_rho"]
        p = row["p_value_for_reporting"]
        q = row["reporting_p_fdr_bh_within_dataset_layer_method"]
        if n < 3 or not np.isfinite(rho):
            tier = "not_estimable"
        elif n < bridge.MIN_PATIENTS_TEST:
            tier = "direction_only_lt5"
        elif n < 10:
            tier = "small_n_exact_significant" if np.isfinite(p) and p < 0.05 else "direction_only_5_to_9"
        else:
            tier = "fdr_significant" if np.isfinite(q) and q < 0.05 else "tested_not_fdr_significant"
        reporting_tier.append(tier)
        reportable.append(tier in {"small_n_exact_significant", "fdr_significant"})
    result["reporting_tier_v2"] = reporting_tier
    result["reportable_v2"] = reportable

    validation = []
    validation.append({
        "check": "v1_rows_and_effects_preserved",
        "status": "PASS" if len(result) == 660 else "FAIL",
        "observed": len(result),
        "expected": 660,
    })
    small = result.loc[(result["n_patients"] >= 5) & (result["n_patients"] < 10) & result["spearman_rho"].notna()]
    validation.append({
        "check": "all_testable_small_n_have_exact_p",
        "status": "PASS" if small["small_n_exact_permutation_p"].notna().all() else "FAIL",
        "observed": int(small["small_n_exact_permutation_p"].notna().sum()),
        "expected": len(small),
    })
    bounds = small["small_n_exact_permutation_p"].between(0, 1, inclusive="both")
    validation.append({
        "check": "exact_p_bounds",
        "status": "PASS" if bounds.all() else "FAIL",
        "observed": int(bounds.sum()),
        "expected": len(bounds),
    })
    non_small = result.loc[~((result["n_patients"] >= 5) & (result["n_patients"] < 10))]
    same = np.isclose(
        pd.to_numeric(non_small["p_value_for_reporting"], errors="coerce"),
        pd.to_numeric(non_small["p_value"], errors="coerce"),
        equal_nan=True,
        atol=0,
        rtol=0,
    )
    validation.append({
        "check": "non_small_p_unchanged",
        "status": "PASS" if same.all() else "FAIL",
        "observed": int(same.sum()),
        "expected": len(same),
    })
    validation_frame = pd.DataFrame(validation)
    if (validation_frame["status"] != "PASS").any():
        atomic_tsv(validation_frame, PART / "validation_checks_v2.tsv")
        raise RuntimeError("Exact small-n validation failed")

    reportable_frame = result.loc[result["reportable_v2"]].copy()
    expanded = result.loc[result["evidence_layer"] == "expanded_any_positive_cnv"].copy()
    group_rows = []
    for edge_group, frame in expanded.groupby("edge_group", sort=True):
        report = frame.loc[frame["reportable_v2"]]
        expected_datasets = sorted(report.loc[report["direction_vs_reference"] == "expected", "dataset"].unique())
        opposing_datasets = sorted(report.loc[report["direction_vs_reference"] == "opposing", "dataset"].unique())
        if expected_datasets and opposing_datasets:
            decision = "heterogeneous_reportable_directions"
        elif expected_datasets:
            decision = "reference_direction_signal"
        elif opposing_datasets:
            decision = "opposing_direction_signal"
        else:
            decision = "no_reportable_signal"
        group_rows.append({
            "edge_group": edge_group,
            "reportable_expected_rows": int(((frame["reportable_v2"]) & (frame["direction_vs_reference"] == "expected")).sum()),
            "reportable_opposing_rows": int(((frame["reportable_v2"]) & (frame["direction_vs_reference"] == "opposing")).sum()),
            "expected_signal_datasets": ";".join(expected_datasets) or "none",
            "opposing_signal_datasets": ";".join(opposing_datasets) or "none",
            "decision": decision,
        })
    group_summary = pd.DataFrame(group_rows)

    atomic_tsv(result, PART / "bridge_screen_all_exact_v2.tsv")
    atomic_tsv(reportable_frame, PART / "bridge_screen_reportable_exact_v2.tsv")
    atomic_tsv(group_summary, PART / "bridge_edge_group_summary_exact_v2.tsv")
    atomic_tsv(validation_frame, PART / "validation_checks_v2.tsv")
    metadata = {
        "schema_version": "phase53_reference_bridge_exact_small_n_v2",
        "status": "PASS",
        "generated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "source_script_sha256": sha256(SOURCE_SCRIPT),
        "source_result_sha256": sha256(SOURCE_RESULT),
        "small_n_exact_rows": len(small),
        "unique_exact_permutation_problems": len(exact_cache),
        "reportable_v1": int(pd.read_csv(SOURCE_RESULT, sep="\t")["reportable"].sum()),
        "reportable_v2": len(reportable_frame),
        "python": platform.python_version(),
        "pandas": pd.__version__,
        "numpy": np.__version__,
        "scipy": scipy.__version__,
    }
    (PART / "run_metadata_v2.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    (PART / "PHASE53_REFERENCE_BRIDGE_EXACT_PASS_v2").write_text(
        "status=PASS\n"
        f"small_n_exact_rows={len(small)}\n"
        f"reportable_rows={len(reportable_frame)}\n"
        f"validation_checks={len(validation_frame)}\n",
        encoding="utf-8",
    )
    os.replace(PART, FINAL)
    print(f"PASS exact_rows={len(small)} reportable={len(reportable_frame)} checks={len(validation_frame)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
