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


DATASETS = ["GSE207493", "GSE264586", "GSE294109"]
PRIMARY = [
    "P01_Malignant_to_SPP1_specificity",
    "P02_SPP1_to_Malignant_specificity",
    "P03_SPP1_to_Tcell_specificity",
    "P04_SPP1_to_Endothelial_specificity",
    "P05_SPP1_total_outgoing_specificity",
    "P06_SPP1_total_incoming_specificity",
]


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


def correlation(frame: pd.DataFrame) -> tuple[int, float, float]:
    clean = frame[["patient_malignant_tymp", "patient_endpoint_value"]].dropna()
    if len(clean) < 3 or clean.iloc[:, 0].nunique() < 3 or clean.iloc[:, 1].nunique() < 3:
        return len(clean), float("nan"), float("nan")
    result = spearmanr(clean.iloc[:, 0].to_numpy(dtype=float), clean.iloc[:, 1].to_numpy(dtype=float))
    return len(clean), float(result.statistic), float(result.pvalue)


def write_summary(path: Path, decision: dict, summary: pd.DataFrame, cross: pd.DataFrame) -> None:
    lines = [
        "# Phase 88 three-cohort fine-state CellChat",
        "",
        f"- Computational status: **{decision['status']}**.",
        f"- Scientific decision: **{decision['decision']}**.",
        f"- Eligible patient-sample units: {decision['patient_sample_units']}.",
        f"- Eligible cohort-specific patients: {decision['patients']}.",
        f"- CellChat interaction rows: {decision['interaction_rows']}.",
        "",
        "## Primary patient-level associations",
        "",
        "| dataset | endpoint | n | rho | p | q | formal support |",
        "|---|---|---:|---:|---:|---:|---|",
    ]
    for row in summary.loc[summary["tier"] == "PRIMARY"].sort_values(["endpoint_id", "dataset"]).itertuples(index=False):
        lines.append(
            f"| {row.dataset} | {row.endpoint_id} | {row.n_patients} | {row.spearman_rho:.3f} | "
            f"{row.p_value:.4g} | {row.BH_q_within_dataset_primary6:.4g} | {row.formal_support} |"
        )
    lines.extend([
        "",
        "## Cross-cohort decisions",
        "",
        "| endpoint | GSE207493 | GSE264586 | GSE294109 | two-external | three-cohort |",
        "|---|---|---|---|---|---|",
    ])
    for row in cross.sort_values("endpoint_id").itertuples(index=False):
        lines.append(
            f"| {row.endpoint_id} | {row.GSE207493_support} | {row.GSE264586_support} | "
            f"{row.GSE294109_support} | {row.two_external_support} | {row.three_cohort_support} |"
        )
    lines.extend([
        "",
        "## Boundary",
        "",
        "CellChat probabilities are expression-compatible communication estimates. They do not prove physical contact, secretion, receptor activation, 5′-MTA mediation, treatment response or causality. NicheNet was not rerun because the package and full frozen resource were unavailable locally.",
        "",
    ])
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--phase-root", required=True, type=Path)
    args = parser.parse_args()
    root = args.phase_root.resolve()
    runs = root / "dataset_runs_v1"
    output = root / "results_v1"
    part = root / "results_v1.part"
    if not runs.is_dir() or output.exists() or part.exists():
        raise RuntimeError("Invalid Phase88 finalization state")
    part.mkdir(parents=True)

    endpoint_frames = []
    qc_frames = []
    interaction_frames = []
    for dataset in DATASETS:
        dataset_root = runs / dataset
        status = json.loads((dataset_root / f"{dataset}_status_v1.json").read_text(encoding="utf-8"))
        if status.get("status") != "PASS_PHASE88_CELLCHAT_DATASET":
            raise RuntimeError(f"Dataset status failed: {dataset}")
        endpoint_frames.append(pd.read_csv(dataset_root / f"{dataset}_unit_endpoints_v1.tsv", sep="\t"))
        qc_frames.append(pd.read_csv(dataset_root / f"{dataset}_unit_qc_v1.tsv", sep="\t"))
        interaction_frames.append(pd.read_csv(dataset_root / f"{dataset}_cellchat_interactions_v1.tsv.gz", sep="\t"))

    units = pd.concat(endpoint_frames, ignore_index=True)
    qc = pd.concat(qc_frames, ignore_index=True)
    interactions = pd.concat(interaction_frames, ignore_index=True)
    if units.duplicated(["dataset", "unit_id", "endpoint_id"]).any():
        raise RuntimeError("Duplicate unit endpoint rows")
    if qc.duplicated(["dataset", "unit_id"]).any() or not qc["status"].eq("PASS_CELLCHAT_UNIT").all():
        raise RuntimeError("Invalid unit QC")

    patient = (
        units.groupby(["dataset", "patient_id", "endpoint_id", "tier"], sort=True)
        .agg(
            patient_malignant_tymp=("malignant_tymp_median", "median"),
            patient_endpoint_value=("endpoint_value", "median"),
            n_samples=("sample_id", "nunique"),
        )
        .reset_index()
    )
    summary_rows = []
    for key, frame in patient.groupby(["dataset", "endpoint_id", "tier"], sort=True):
        n, rho, p_value = correlation(frame)
        summary_rows.append({
            "dataset": key[0],
            "endpoint_id": key[1],
            "tier": key[2],
            "n_patients": n,
            "spearman_rho": rho,
            "p_value": p_value,
        })
    summary = pd.DataFrame(summary_rows)
    summary["BH_q_within_dataset_primary6"] = np.nan
    for dataset, index in summary.loc[summary["tier"] == "PRIMARY"].groupby("dataset").groups.items():
        summary.loc[index, "BH_q_within_dataset_primary6"] = bh(summary.loc[index, "p_value"])
    summary["formal_support"] = (
        summary["tier"].eq("PRIMARY")
        & (summary["n_patients"] >= 10)
        & (summary["spearman_rho"] > 0)
        & (summary["BH_q_within_dataset_primary6"] < 0.05)
    )
    summary["small_n_nominal_support"] = (
        summary["tier"].eq("PRIMARY")
        & summary["n_patients"].between(5, 9)
        & (summary["spearman_rho"] > 0)
        & (summary["p_value"] < 0.05)
    )

    cross_rows = []
    for endpoint in PRIMARY:
        row = {"endpoint_id": endpoint}
        support_values = {}
        for dataset in DATASETS:
            current = summary.loc[(summary["dataset"] == dataset) & (summary["endpoint_id"] == endpoint)]
            value = len(current) == 1 and bool(current["formal_support"].iloc[0])
            support_values[dataset] = value
            row[f"{dataset}_support"] = value
            row[f"{dataset}_rho"] = float(current["spearman_rho"].iloc[0]) if len(current) == 1 else float("nan")
        row["two_external_support"] = bool(support_values["GSE264586"] and support_values["GSE294109"])
        row["three_cohort_support"] = bool(row["two_external_support"] and support_values["GSE207493"])
        cross_rows.append(row)
    cross = pd.DataFrame(cross_rows)

    p01 = bool(cross.loc[cross["endpoint_id"] == PRIMARY[0], "two_external_support"].iloc[0])
    p02 = bool(cross.loc[cross["endpoint_id"] == PRIMARY[1], "two_external_support"].iloc[0])
    any_external = bool(cross["two_external_support"].any())
    if p01 and p02:
        scientific_decision = "TWO_EXTERNAL_BIDIRECTIONAL_MALIGNANT_SPP1_COMMUNICATION_SUPPORTED"
    elif any_external:
        scientific_decision = "PARTIAL_TWO_EXTERNAL_FINE_STATE_COMMUNICATION_SUPPORT"
    else:
        scientific_decision = "NO_TWO_EXTERNAL_FINE_STATE_COMMUNICATION_SUPPORT"

    decision = {
        "status": "PASS_COMPUTATION_COMPLETE",
        "decision": scientific_decision,
        "p01_two_external_support": p01,
        "p02_two_external_support": p02,
        "any_primary_two_external_support": any_external,
        "patient_sample_units": int(qc[["dataset", "unit_id"]].drop_duplicates().shape[0]),
        "patients": int(qc[["dataset", "patient_id"]].drop_duplicates().shape[0]),
        "sampled_cells": int(qc["sampled_cells"].sum()),
        "interaction_rows": int(len(interactions)),
        "preregistration_sha256": sha256(root / "PHASE88_PREREGISTRATION_v1.md"),
        "cellchat_db_sha256": sha256(root / "inputs_v1" / "cellchat_resource_v1" / "CellChatDB.human_v1.rds"),
        "nichenet_rerun": False,
    }

    units.to_csv(part / "phase88_unit_endpoints_v1.tsv", sep="\t", index=False)
    qc.to_csv(part / "phase88_unit_qc_all_v1.tsv", sep="\t", index=False)
    interactions.to_csv(part / "phase88_cellchat_interactions_all_v1.tsv.gz", sep="\t", index=False, compression="gzip")
    patient.to_csv(part / "phase88_patient_endpoints_v1.tsv", sep="\t", index=False)
    summary.to_csv(part / "phase88_dataset_endpoint_summary_v1.tsv", sep="\t", index=False)
    cross.to_csv(part / "phase88_cross_cohort_decisions_v1.tsv", sep="\t", index=False)
    with (part / "PHASE88_DECISION_v1.json").open("w", encoding="utf-8") as handle:
        json.dump(decision, handle, indent=2, ensure_ascii=False)
        handle.write("\n")
    write_summary(part / "PHASE88_CELLCHAT_SUMMARY_v1.md", decision, summary, cross)
    (part / "PHASE88_COMPUTATION_PASS_v1").write_text("PASS_COMPUTATION_COMPLETE\n", encoding="utf-8")
    manifest_rows = []
    for path in sorted(part.iterdir()):
        if path.name == "FILE_MANIFEST_SHA256_v1.tsv":
            continue
        manifest_rows.append({"file": path.name, "bytes": path.stat().st_size, "sha256": sha256(path)})
    pd.DataFrame(manifest_rows).to_csv(part / "FILE_MANIFEST_SHA256_v1.tsv", sep="\t", index=False)
    os.replace(part, output)
    print(json.dumps(decision, ensure_ascii=False))


if __name__ == "__main__":
    main()
