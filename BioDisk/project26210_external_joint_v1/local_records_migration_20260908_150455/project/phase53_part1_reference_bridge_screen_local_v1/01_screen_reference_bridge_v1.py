#!/usr/bin/env python3
"""Patient-level screen of the reference-paper CAF/TAM/T-cell bridge."""

from __future__ import annotations

import hashlib
import importlib.util
import json
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
PHASE14 = ROOT / "phase14_patient_program_validation_local_v1"
PHASE52 = ROOT / "phase52_part1_incremental_three_tier_local_v1" / "analysis_local_v1"
SOURCE_ANALYSIS = PHASE14 / "03_run_patient_level_analysis_v1.py"
SOURCE_THREE_TIER = PHASE52 / "three_tier_incremental_associations_v1.tsv"
FINAL = HERE / "analysis_local_v1"
PART = HERE / "analysis_local_v1.part"

DATASETS = ["GSE178481", "GSE202374", "GSE207493", "GSE264586", "GSE294109"]
LAYERS = {
    "strict_conservative_cnv": PHASE14 / "runs" / "20260812_final_cnv_gated_v6",
    "expanded_any_positive_cnv": PHASE14 / "runs" / "20260812_expanded_malignancy_v2",
    "provisional_tumor_like": PHASE14 / "runs" / "20260812_082000_provisional_v3",
}
METHODS = {
    "primary_gene_z_mean": "mean_primary",
    "sensitivity_mean_log_normalized": "mean_sensitivity",
}
MIN_CELLS = 20
MIN_PATIENTS_TEST = 5


RELATIONS = [
    ("B01_F1_to_POSTN_CAF", "tumor_to_caf", "f1_top100_weighted", "malignant_provisional", "postn_ecm_proxy", "fibroblast", 1),
    ("B02_TYMP_to_POSTN_CAF", "tumor_to_caf", "tymp_expression", "malignant_provisional", "postn_ecm_proxy", "fibroblast", 1),
    ("B03_POSTN_CAF_to_APOE_TAM", "caf_to_tam", "postn_ecm_proxy", "fibroblast", "apoe_tam_expression_proxy", "myeloid", 1),
    ("B04_POSTN_CAF_to_C1QC", "caf_to_tam", "postn_ecm_proxy", "fibroblast", "c1qc_state", "myeloid", 1),
    ("B05_POSTN_CAF_to_SPP1", "caf_to_tam", "postn_ecm_proxy", "fibroblast", "spp1_state_no_spp1", "myeloid", 1),
    ("B06_POSTN_CAF_to_myeloid_suppression", "caf_to_tam", "postn_ecm_proxy", "fibroblast", "myeloid_immunosuppression", "myeloid", 1),
    ("B07_APOE_TAM_to_Tdysfunction", "tam_to_tcell", "apoe_tam_expression_proxy", "myeloid", "tcell_dysfunction", "t_nk", 1),
    ("B08_APOE_TAM_to_cytotoxicity", "tam_to_tcell", "apoe_tam_expression_proxy", "myeloid", "cytotoxicity", "t_nk", -1),
    ("B09_C1QC_to_Tdysfunction", "tam_to_tcell", "c1qc_state", "myeloid", "tcell_dysfunction", "t_nk", 1),
    ("B10_C1QC_to_cytotoxicity", "tam_to_tcell", "c1qc_state", "myeloid", "cytotoxicity", "t_nk", -1),
    ("B11_SPP1_to_Tdysfunction", "tam_to_tcell", "spp1_state_no_spp1", "myeloid", "tcell_dysfunction", "t_nk", 1),
    ("B12_SPP1_to_cytotoxicity", "tam_to_tcell", "spp1_state_no_spp1", "myeloid", "cytotoxicity", "t_nk", -1),
    ("B13_APOE_TAM_to_vascular", "tam_to_vascular", "apoe_tam_expression_proxy", "myeloid", "vascular_remodeling", "myeloid", 1),
    ("B14_SPP1_to_vascular", "tam_to_vascular", "spp1_state_no_spp1", "myeloid", "vascular_remodeling", "myeloid", 1),
    ("B15_F1_to_APOE_TAM", "tumor_to_tam", "f1_top100_weighted", "malignant_provisional", "apoe_tam_expression_proxy", "myeloid", 1),
    ("B16_F1_to_C1QC", "tumor_to_tam", "f1_top100_weighted", "malignant_provisional", "c1qc_state", "myeloid", 1),
    ("B17_F1_to_SPP1", "tumor_to_tam", "f1_top100_weighted", "malignant_provisional", "spp1_state_no_spp1", "myeloid", 1),
    ("B18_TYMP_to_APOE_TAM", "tumor_to_tam", "tymp_expression", "malignant_provisional", "apoe_tam_expression_proxy", "myeloid", 1),
    ("B19_TYMP_to_C1QC", "tumor_to_tam", "tymp_expression", "malignant_provisional", "c1qc_state", "myeloid", 1),
    ("B20_TYMP_to_SPP1", "tumor_to_tam", "tymp_expression", "malignant_provisional", "spp1_state_no_spp1", "myeloid", 1),
    ("B21_APOE_TAM_to_myeloid_suppression", "tam_state_internal", "apoe_tam_expression_proxy", "myeloid", "myeloid_immunosuppression", "myeloid", 1),
    ("B22_SPP1_to_myeloid_suppression", "tam_state_internal", "spp1_state_no_spp1", "myeloid", "myeloid_immunosuppression", "myeloid", 1),
]
MALIGNANT_DEPENDENT = {row[0] for row in RELATIONS if row[3] == "malignant_provisional" or row[5] == "malignant_provisional"}
NONMALIGNANT = {row[0] for row in RELATIONS} - MALIGNANT_DEPENDENT


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


def classify(n: int, rho: float, p: float, q: float, expected: int) -> tuple[str, str]:
    if n < 3 or not np.isfinite(rho):
        return "not_estimable", "not_estimable"
    direction = "expected" if np.sign(rho) == expected else "opposing" if rho != 0 else "zero"
    if n < MIN_PATIENTS_TEST:
        return "direction_only_lt5", direction
    significant = bool(p < 0.05) if n < 10 else bool(np.isfinite(q) and q < 0.05)
    if significant:
        return ("small_n_nominal_significant" if n < 10 else "fdr_significant"), direction
    return ("direction_only_5_to_9" if n < 10 else "tested_not_fdr_significant"), direction


def analyze_relation(
    patient: pd.DataFrame,
    dataset: str,
    layer: str,
    method: str,
    value_col: str,
    relation: tuple[str, str, str, str, str, str, int],
) -> dict[str, object]:
    relation_id, edge_group, x_sig, x_comp, y_sig, y_comp, expected = relation
    x = patient.loc[
        (patient["signature_id"] == x_sig)
        & (patient["analysis_compartment"] == x_comp)
        & (patient["n_cells"] >= MIN_CELLS),
        ["patient_key", "patient_id", value_col, "n_cells"],
    ].rename(columns={value_col: "x", "n_cells": "x_cells"})
    y = patient.loc[
        (patient["signature_id"] == y_sig)
        & (patient["analysis_compartment"] == y_comp)
        & (patient["n_cells"] >= MIN_CELLS),
        ["patient_key", value_col, "n_cells"],
    ].rename(columns={value_col: "y", "n_cells": "y_cells"})
    paired = x.merge(y, on="patient_key", how="inner", validate="one_to_one")
    paired = paired.loc[np.isfinite(paired["x"]) & np.isfinite(paired["y"])].copy()
    n = len(paired)
    rho = float(stats.spearmanr(paired["x"], paired["y"]).statistic) if n >= 3 and paired["x"].nunique() > 1 and paired["y"].nunique() > 1 else np.nan
    p = float(stats.spearmanr(paired["x"], paired["y"]).pvalue) if n >= MIN_PATIENTS_TEST and np.isfinite(rho) else np.nan
    loo_values: list[float] = []
    if n >= 4 and np.isfinite(rho):
        for index in range(n):
            sub = paired.drop(paired.index[index])
            if sub["x"].nunique() > 1 and sub["y"].nunique() > 1:
                loo_values.append(float(stats.spearmanr(sub["x"], sub["y"]).statistic))
    loo_stability = float(np.mean(np.sign(loo_values) == np.sign(rho))) if loo_values else np.nan
    return {
        "dataset": dataset,
        "evidence_layer": layer,
        "method": method,
        "relation": relation_id,
        "edge_group": edge_group,
        "x_signature": x_sig,
        "x_compartment": x_comp,
        "y_signature": y_sig,
        "y_compartment": y_comp,
        "expected_direction": expected,
        "n_patients": n,
        "patient_ids": ";".join(sorted(paired["patient_id"].astype(str))),
        "spearman_rho": rho,
        "p_value": p,
        "loo_min_rho": min(loo_values) if loo_values else np.nan,
        "loo_max_rho": max(loo_values) if loo_values else np.nan,
        "loo_sign_stability_fraction": loo_stability,
    }


def main() -> int:
    if FINAL.exists() or PART.exists():
        raise FileExistsError(f"Refusing to overwrite {FINAL} or {PART}")
    required = [SOURCE_ANALYSIS, SOURCE_THREE_TIER]
    sample_sources: list[Path] = []
    for run in LAYERS.values():
        for dataset in DATASETS:
            sample_sources.append(run / dataset / "tables" / "sample_program_summary_v1.tsv")
    missing = [str(path) for path in required + sample_sources if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing frozen source(s): " + "; ".join(missing))
    PART.mkdir(parents=True)

    spec = importlib.util.spec_from_file_location("phase14_patient", SOURCE_ANALYSIS)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)

    rows: list[dict[str, object]] = []
    coverage_rows: list[dict[str, object]] = []
    for layer, run in LAYERS.items():
        for dataset in DATASETS:
            source = run / dataset / "tables" / "sample_program_summary_v1.tsv"
            sample = pd.read_csv(source, sep="\t")
            patient, _ = module.build_cross_sectional_patient_programs(sample)
            signatures = set(patient["signature_id"].astype(str))
            compartments = set(patient["analysis_compartment"].astype(str))
            coverage_rows.append(
                {
                    "dataset": dataset,
                    "evidence_layer": layer,
                    "patients": patient["patient_id"].astype(str).nunique(),
                    "postn_fibroblast_available": bool(
                        "postn_ecm_proxy" in signatures and "fibroblast" in compartments
                    ),
                    "apoe_tam_available": "apoe_tam_expression_proxy" in signatures,
                    "spp1_available": "spp1_state_no_spp1" in signatures,
                    "tcell_endpoints_available": {"tcell_dysfunction", "cytotoxicity"}.issubset(signatures),
                    "signature_count": len(signatures),
                    "compartments": ";".join(sorted(compartments)),
                }
            )
            for method, value_col in METHODS.items():
                for relation in RELATIONS:
                    rows.append(analyze_relation(patient, dataset, layer, method, value_col, relation))

    result = pd.DataFrame(rows)
    result["p_fdr_bh_within_dataset_layer_method"] = result.groupby(
        ["dataset", "evidence_layer", "method"], group_keys=False
    )["p_value"].apply(bh)
    classifications = result.apply(
        lambda row: classify(
            int(row["n_patients"]),
            float(row["spearman_rho"]) if pd.notna(row["spearman_rho"]) else np.nan,
            float(row["p_value"]) if pd.notna(row["p_value"]) else np.nan,
            float(row["p_fdr_bh_within_dataset_layer_method"]) if pd.notna(row["p_fdr_bh_within_dataset_layer_method"]) else np.nan,
            int(row["expected_direction"]),
        ),
        axis=1,
    )
    result[["reporting_tier", "direction_vs_reference"]] = pd.DataFrame(classifications.tolist(), index=result.index)
    result["reportable"] = result["reporting_tier"].isin(["small_n_nominal_significant", "fdr_significant"])
    result = result.sort_values(["dataset", "evidence_layer", "method", "relation"]).reset_index(drop=True)

    validation: list[dict[str, object]] = []
    keys = ["dataset", "evidence_layer", "method", "relation"]
    expected_rows = len(DATASETS) * len(LAYERS) * len(METHODS) * len(RELATIONS)
    validation.append({
        "check": "complete_unique_grid",
        "status": "PASS" if len(result) == expected_rows and not result.duplicated(keys).any() else "FAIL",
        "observed": f"rows={len(result)};duplicates={int(result.duplicated(keys).sum())}",
        "expected": f"rows={expected_rows};duplicates=0",
    })
    invalid_p = result.loc[(result["n_patients"] < MIN_PATIENTS_TEST) & result["p_value"].notna()]
    validation.append({
        "check": "no_tests_below_five_patients",
        "status": "PASS" if invalid_p.empty else "FAIL",
        "observed": len(invalid_p),
        "expected": 0,
    })

    # Malignant-dependent patient eligibility must be strict <= expanded <= provisional.
    malignant = result.loc[result["relation"].isin(MALIGNANT_DEPENDENT)]
    pivot = malignant.pivot(index=["dataset", "method", "relation"], columns="evidence_layer", values="n_patients")
    nested = (pivot["strict_conservative_cnv"] <= pivot["expanded_any_positive_cnv"]) & (
        pivot["expanded_any_positive_cnv"] <= pivot["provisional_tumor_like"]
    )
    validation.append({
        "check": "malignant_layer_patient_nesting",
        "status": "PASS" if bool(nested.all()) else "FAIL",
        "observed": f"passing={int(nested.sum())};rows={len(nested)}",
        "expected": f"passing={len(nested)}",
    })

    # Non-malignant bridge edges must reproduce across all malignancy layers.
    nonmal = result.loc[result["relation"].isin(NONMALIGNANT)]
    for field in ("n_patients", "spearman_rho", "p_value"):
        wide = nonmal.pivot(index=["dataset", "method", "relation"], columns="evidence_layer", values=field)
        values = wide.to_numpy(float)
        same = np.all(
            np.isclose(values[:, [0]], values, equal_nan=True, atol=1e-12, rtol=0),
            axis=1,
        )
        validation.append({
            "check": f"nonmalignant_edges_layer_invariant_{field}",
            "status": "PASS" if bool(same.all()) else "FAIL",
            "observed": f"passing={int(same.sum())};rows={len(same)}",
            "expected": f"passing={len(same)}",
        })

    # Reproduce overlapping Phase 52 full-cohort F1/TYMP-SPP1 values.
    frozen = pd.read_csv(SOURCE_THREE_TIER, sep="\t")
    mapping = {"B17_F1_to_SPP1": "H1_F1_to_SPP1", "B20_TYMP_to_SPP1": "H1_TYMP_to_SPP1"}
    compare_rows = []
    for dataset in ("GSE178481", "GSE207493"):
        for layer in LAYERS:
            for method in METHODS:
                for current_relation, frozen_relation in mapping.items():
                    current = result.loc[
                        (result["dataset"] == dataset)
                        & (result["evidence_layer"] == layer)
                        & (result["method"] == method)
                        & (result["relation"] == current_relation)
                    ].iloc[0]
                    prior = frozen.loc[
                        (frozen["dataset"] == dataset)
                        & (frozen["subset"] == "expanded_full")
                        & (frozen["evidence_layer"] == layer)
                        & (frozen["method"] == method)
                        & (frozen["relation"] == frozen_relation)
                    ].iloc[0]
                    for field in ("n_patients", "spearman_rho", "p_value"):
                        a = pd.to_numeric(pd.Series([current[field]]), errors="coerce").iloc[0]
                        b = pd.to_numeric(pd.Series([prior[field]]), errors="coerce").iloc[0]
                        compare_rows.append(bool((pd.isna(a) and pd.isna(b)) or np.isclose(a, b, atol=1e-12, rtol=0)))
    validation.append({
        "check": "phase52_core_relation_reproduction",
        "status": "PASS" if all(compare_rows) else "FAIL",
        "observed": f"passing={sum(compare_rows)};comparisons={len(compare_rows)}",
        "expected": f"passing={len(compare_rows)}",
    })

    coverage = pd.DataFrame(coverage_rows)
    g207_postn = coverage.loc[coverage["dataset"] == "GSE207493", "postn_fibroblast_available"]
    validation.append({
        "check": "gse207493_postn_unavailable_recorded",
        "status": "PASS" if not bool(g207_postn.any()) else "FAIL",
        "observed": f"available_layers={int(g207_postn.sum())}",
        "expected": "available_layers=0",
    })
    validation_frame = pd.DataFrame(validation)
    if (validation_frame["status"] != "PASS").any():
        atomic_tsv(validation_frame, PART / "validation_checks_v1.tsv")
        raise RuntimeError("Bridge screen validation failed")

    reportable = result.loc[result["reportable"]].copy()
    expanded = result.loc[result["evidence_layer"] == "expanded_any_positive_cnv"].copy()
    group_rows = []
    for edge_group, frame in expanded.groupby("edge_group", sort=True):
        report = frame.loc[frame["reportable"]]
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
            "expanded_rows": len(frame),
            "estimable_rows": int(frame["spearman_rho"].notna().sum()),
            "reportable_expected_rows": int(((frame["reportable"]) & (frame["direction_vs_reference"] == "expected")).sum()),
            "reportable_opposing_rows": int(((frame["reportable"]) & (frame["direction_vs_reference"] == "opposing")).sum()),
            "expected_signal_datasets": ";".join(expected_datasets) or "none",
            "opposing_signal_datasets": ";".join(opposing_datasets) or "none",
            "decision": decision,
        })
    group_summary = pd.DataFrame(group_rows)
    dataset_summary = result.groupby(["dataset", "evidence_layer"], as_index=False).agg(
        relations=("relation", "size"),
        estimable=("spearman_rho", lambda x: int(x.notna().sum())),
        tested=("p_value", lambda x: int(x.notna().sum())),
        reportable=("reportable", "sum"),
        reportable_expected=("direction_vs_reference", lambda x: int(((result.loc[x.index, "reportable"]) & (x == "expected")).sum())),
        reportable_opposing=("direction_vs_reference", lambda x: int(((result.loc[x.index, "reportable"]) & (x == "opposing")).sum())),
    )

    atomic_tsv(result, PART / "bridge_screen_all_v1.tsv")
    atomic_tsv(reportable, PART / "bridge_screen_reportable_v1.tsv")
    atomic_tsv(coverage, PART / "bridge_endpoint_coverage_v1.tsv")
    atomic_tsv(group_summary, PART / "bridge_edge_group_summary_v1.tsv")
    atomic_tsv(dataset_summary, PART / "bridge_dataset_summary_v1.tsv")
    atomic_tsv(validation_frame, PART / "validation_checks_v1.tsv")

    metadata = {
        "schema_version": "phase53_part1_reference_bridge_screen_v1",
        "status": "PASS",
        "generated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "datasets": DATASETS,
        "layers": list(LAYERS),
        "methods": list(METHODS),
        "relations": len(RELATIONS),
        "result_rows": len(result),
        "reportable_rows": len(reportable),
        "minimum_cells_per_patient_compartment": MIN_CELLS,
        "minimum_patients_for_test": MIN_PATIENTS_TEST,
        "source_analysis_sha256": sha256(SOURCE_ANALYSIS),
        "source_three_tier_sha256": sha256(SOURCE_THREE_TIER),
        "sample_source_sha256": {str(path.relative_to(ROOT)): sha256(path) for path in sample_sources},
        "python": platform.python_version(),
        "pandas": pd.__version__,
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "normalization_performed": False,
        "integration_performed": False,
        "clustering_performed": False,
        "annotation_performed": False,
    }
    (PART / "run_metadata_v1.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    (PART / "PHASE53_REFERENCE_BRIDGE_SCREEN_PASS_v1").write_text(
        "status=PASS\n"
        f"datasets={len(DATASETS)}\n"
        f"layers={len(LAYERS)}\n"
        f"relations={len(RELATIONS)}\n"
        f"result_rows={len(result)}\n"
        f"reportable_rows={len(reportable)}\n"
        f"validation_checks={len(validation_frame)}\n"
        "normalization_performed=0\n"
        "integration_performed=0\n"
        "clustering_performed=0\n",
        encoding="utf-8",
    )
    os.replace(PART, FINAL)
    print(f"PASS rows={len(result)} reportable={len(reportable)} checks={len(validation_frame)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
