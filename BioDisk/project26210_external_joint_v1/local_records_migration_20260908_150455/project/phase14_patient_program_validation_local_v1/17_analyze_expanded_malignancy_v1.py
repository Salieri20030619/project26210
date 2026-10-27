#!/usr/bin/env python3
"""Summarize stability and inference after adding supportive malignant cells."""

from __future__ import annotations

import json
import os
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


HERE = Path(__file__).resolve().parent
STRICT_RUN = HERE / "runs" / "20260812_final_cnv_gated_v6"
EXPANDED_RUN = HERE / "runs" / "20260812_expanded_malignancy_v2"
STRICT_ANALYSIS = HERE / "patient_analysis_final_cnv_v6"
EXPANDED_ANALYSIS = HERE / "patient_analysis_expanded_malignancy_v1"
OUT = HERE / "posthoc_expanded_malignancy_diagnostic_v1"
DATASETS = ["GSE178481", "GSE202374", "GSE207493", "GSE264586", "GSE294109"]
MIN_CELLS = 20


def atomic_tsv(frame: pd.DataFrame, path: Path) -> None:
    part = path.with_name(path.name + ".part")
    frame.to_csv(part, sep="\t", index=False, encoding="utf-8", na_rep="NA", lineterminator="\n")
    os.replace(part, path)


def atomic_json(payload: dict, path: Path) -> None:
    part = path.with_name(path.name + ".part")
    part.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    os.replace(part, path)


def selected_patient_programs(run: Path, dataset: str) -> pd.DataFrame:
    sample = pd.read_csv(run / dataset / "tables" / "sample_program_summary_v1.tsv", sep="\t")
    sample = sample[sample["analysis_compartment"] == "malignant_provisional"].copy()
    if dataset == "GSE178481":
        sample = sample[~sample["sample_id"].astype(str).str.contains("Normal", case=False, regex=False)]
    if dataset == "GSE264586":
        sample = sample[sample["timepoint"] == "Baseline"]
    rows = []
    fields = ["dataset", "patient_key", "patient_id", "signature_id", "analysis_compartment"]
    for keys, group in sample.groupby(fields, dropna=False, sort=True):
        weights = group["n_cells"].to_numpy(float)
        rows.append(
            dict(
                zip(fields, keys),
                n_cells=int(weights.sum()),
                mean_primary=float(np.average(group["mean_primary"], weights=weights)),
                mean_sensitivity=float(np.average(group["mean_sensitivity"], weights=weights)),
            )
        )
    return pd.DataFrame(rows)


def main() -> None:
    if OUT.exists():
        raise FileExistsError(OUT)
    OUT.mkdir()

    input_status = pd.read_csv(EXPANDED_RUN / "EXPANDED_MALIGNANCY_INPUT_STATUS_v1.tsv", sep="\t")
    input_status["supportive_fraction_of_expanded"] = input_status["supportive_cells"] / input_status["expanded_malignant_cells"]
    input_status["expanded_fraction_of_tumor_like"] = input_status["expanded_malignant_cells"] / input_status["tumor_like_cells"]
    atomic_tsv(input_status, OUT / "expanded_cell_and_patient_coverage_v1.tsv")

    evidence_frames = []
    for dataset in DATASETS:
        frame = pd.read_csv(EXPANDED_RUN / dataset / "tables" / "expanded_malignancy_evidence_counts_v1.tsv", sep="\t")
        evidence_frames.append(frame.groupby(["dataset", "malignancy_evidence_tier"], as_index=False)["n_cells"].sum())
    evidence = pd.concat(evidence_frames, ignore_index=True)
    atomic_tsv(evidence, OUT / "malignancy_evidence_tier_counts_v1.tsv")

    strict = pd.read_csv(STRICT_ANALYSIS / "within_dataset_patient_associations_v6.tsv", sep="\t")
    expanded = pd.read_csv(EXPANDED_ANALYSIS / "within_dataset_patient_associations_expanded_v1.tsv", sep="\t")
    malignant = expanded["x_compartment"].str.contains("malignant") | expanded["y_compartment"].str.contains("malignant")
    keys = ["dataset", "method", "relation"]
    comparison = strict[malignant].merge(expanded[malignant], on=keys, suffixes=("_strict", "_expanded"), validate="one_to_one")
    comparison["testability_change"] = np.select(
        [comparison["spearman_rho_strict"].isna() & comparison["spearman_rho_expanded"].notna(),
         comparison["spearman_rho_strict"].notna() & comparison["spearman_rho_expanded"].notna()],
        ["gained_with_expansion", "testable_both"], default="untestable_both"
    )
    comparison["rho_delta_expanded_minus_strict"] = comparison["spearman_rho_expanded"] - comparison["spearman_rho_strict"]
    comparison["same_rho_sign_when_both_tested"] = np.where(
        comparison["spearman_rho_strict"].notna() & comparison["spearman_rho_expanded"].notna(),
        np.sign(comparison["spearman_rho_strict"]) == np.sign(comparison["spearman_rho_expanded"]), np.nan
    )
    atomic_tsv(comparison, OUT / "strict_vs_expanded_malignant_associations_v1.tsv")

    tested = expanded[expanded["spearman_rho"].notna()].copy()
    tested["bootstrap_ci_excludes_zero"] = (tested["bootstrap_ci_low"] > 0) | (tested["bootstrap_ci_high"] < 0)
    nominal = tested[tested["p_value"] < 0.05].sort_values("p_value")
    atomic_tsv(nominal, OUT / "expanded_nominal_associations_v1.tsv")

    concordance_rows = []
    for dataset in DATASETS:
        high = selected_patient_programs(STRICT_RUN, dataset)
        exp = selected_patient_programs(EXPANDED_RUN, dataset)
        shared = high.merge(exp, on=["dataset", "patient_key", "patient_id", "signature_id", "analysis_compartment"], suffixes=("_strict", "_expanded"))
        shared = shared[(shared["n_cells_strict"] >= MIN_CELLS) & (shared["n_cells_expanded"] >= MIN_CELLS)]
        for signature, group in shared.groupby("signature_id", sort=True):
            for method, column in [("primary_gene_z_mean", "mean_primary"), ("sensitivity_mean_log_normalized", "mean_sensitivity")]:
                x = group[f"{column}_strict"].to_numpy(float)
                y = group[f"{column}_expanded"].to_numpy(float)
                rho = float(stats.spearmanr(x, y).statistic) if len(x) >= 3 and np.unique(x).size > 1 and np.unique(y).size > 1 else np.nan
                concordance_rows.append({
                    "dataset": dataset,
                    "signature_id": signature,
                    "method": method,
                    "shared_patients_ge20": len(group),
                    "spearman_strict_vs_expanded": rho,
                    "median_expanded_minus_strict": float(np.median(y - x)) if len(x) else np.nan,
                    "max_abs_expanded_minus_strict": float(np.max(np.abs(y - x))) if len(x) else np.nan,
                })
    concordance = pd.DataFrame(concordance_rows)
    atomic_tsv(concordance, OUT / "patient_program_strict_vs_expanded_concordance_v1.tsv")

    malignant_comparison = comparison[comparison["testability_change"] == "testable_both"]
    status = {
        "schema_version": "phase14_expanded_malignancy_diagnostic_v1",
        "status": "PASS_POSTHOC_EXPANDED_MALIGNANCY",
        "completed_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "definition": "tumor-like plus at least one positive CNV method; rebuild calls take priority for rebuild-assessed cells",
        "high_confidence_cells": int(input_status["high_confidence_cells"].sum()),
        "supportive_cells": int(input_status["supportive_cells"].sum()),
        "expanded_cells": int(input_status["expanded_malignant_cells"].sum()),
        "tumor_like_cells": int(input_status["tumor_like_cells"].sum()),
        "association_rows": int(len(expanded)),
        "association_rows_tested": int(len(tested)),
        "malignant_rows_tested": int(expanded.loc[malignant, "spearman_rho"].notna().sum()),
        "nominal_p_lt_0_05": int((tested["p_value"] < 0.05).sum()),
        "fdr_q_lt_0_05": int((tested["p_fdr_bh_within_dataset_method"] < 0.05).sum()),
        "bootstrap_ci_excludes_zero": int(tested["bootstrap_ci_excludes_zero"].sum()),
        "malignant_rows_gained_testability": int((comparison["testability_change"] == "gained_with_expansion").sum()),
        "same_sign_both_testable": int(malignant_comparison["same_rho_sign_when_both_tested"].sum()),
        "both_testable": int(len(malignant_comparison)),
        "formal_strict_v6_replaced": False,
        "warning": "The expanded tier was evaluated post hoc and is recommended as a biological working layer with strict high-confidence sensitivity analysis, not as a retrospective replacement of formal v6.",
    }
    atomic_json(status, OUT / "EXPANDED_MALIGNANCY_DIAGNOSTIC_STATUS_v1.json")
    print(json.dumps(status, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
