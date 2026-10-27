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
EXPECTED_UNITS = {"GSE207493": 16, "GSE264586": 6, "GSE294109": 14}
EXPECTED_CELLS = {"GSE207493": 22941, "GSE264586": 7347, "GSE294109": 16676}
EXPECTED_INTERACTIONS = {"GSE207493": 71725, "GSE264586": 29284, "GSE294109": 60856}
PRIMARY = [
    "P01_Malignant_to_SPP1_specificity",
    "P02_SPP1_to_Malignant_specificity",
    "P03_SPP1_to_Tcell_specificity",
    "P04_SPP1_to_Endothelial_specificity",
    "P05_SPP1_total_outgoing_specificity",
    "P06_SPP1_total_incoming_specificity",
]
EXPECTED_PREREG = "79c69ee3fe41ab8ed3219b6bd1c91b0cec9ab50b4104d0a557e0badf1706d6ae"
EXPECTED_DB = "d89941845882ca77c1e5e258b1b45385e8e1b8264afc2fa06ffcbe44ce23c4b9"


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def check(condition: bool, label: str, checks: list[dict]) -> None:
    checks.append({"check": label, "status": "PASS" if condition else "FAIL"})
    if not condition:
        raise AssertionError(label)


def bh(p_values: pd.Series) -> pd.Series:
    out = pd.Series(np.nan, index=p_values.index, dtype=float)
    keep = p_values.notna()
    p = p_values.loc[keep].to_numpy(float)
    if not len(p):
        return out
    order = np.argsort(p)
    ordered = p[order]
    adj = ordered * len(ordered) / np.arange(1, len(ordered) + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    restored = np.empty_like(adj)
    restored[order] = np.clip(adj, 0, 1)
    out.loc[keep] = restored
    return out


def correlation(frame: pd.DataFrame) -> tuple[int, float, float]:
    clean = frame[["patient_malignant_tymp", "patient_endpoint_value"]].dropna()
    if len(clean) < 3 or clean.iloc[:, 0].nunique() < 3 or clean.iloc[:, 1].nunique() < 3:
        return len(clean), math.nan, math.nan
    result = spearmanr(clean.iloc[:, 0].to_numpy(float), clean.iloc[:, 1].to_numpy(float))
    return len(clean), float(result.statistic), float(result.pvalue)


def equal_numeric(left: pd.Series, right: pd.Series) -> bool:
    return bool(np.allclose(left.to_numpy(float), right.to_numpy(float), rtol=1e-12, atol=1e-12, equal_nan=True))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--phase-root", required=True, type=Path)
    args = parser.parse_args()
    root = args.phase_root.resolve()
    result_dir = root / "results_v1"
    part = root / "validation_v1.part"
    final = root / "validation_v1"
    if part.exists() or final.exists():
        raise RuntimeError("Validation output already exists")
    part.mkdir(parents=True)
    checks: list[dict] = []

    check(digest(root / "PHASE88_PREREGISTRATION_v1.md") == EXPECTED_PREREG, "frozen preregistration hash", checks)
    check(digest(root / "inputs_v1/cellchat_resource_v1/CellChatDB.human_v1.rds") == EXPECTED_DB, "frozen CellChat DB hash", checks)
    check((root / "PHASE88_PREREGISTRATION_AMENDMENT_A_v1.md").is_file(), "amendment A preserved", checks)
    check((root / "prepared_inputs_v1_failed_attempt1/RUN_FAILURE_v1.txt").is_file(), "input failed attempt preserved", checks)
    check((root / "dataset_runs_v1/GSE264586_failed_attempt1/RUN_FAILURE_v1.txt").is_file(), "CellChat failed attempt preserved", checks)

    input_status = pd.read_csv(root / "prepared_inputs_v1/phase88_input_status_v1.tsv", sep="\t")
    check(set(input_status.dataset) == set(DATASETS), "three prepared datasets", checks)
    check(input_status.status.eq("PASS_PHASE88_INPUT_PREPARATION").all(), "input statuses PASS", checks)
    check(dict(zip(input_status.dataset, input_status.eligible_units)) == EXPECTED_UNITS, "prepared eligible unit counts", checks)
    check(dict(zip(input_status.dataset, input_status.sampled_cells)) == EXPECTED_CELLS, "prepared sampled cell counts", checks)

    sampling = pd.read_csv(root / "prepared_inputs_v1/phase88_sampling_qc_v1.tsv", sep="\t")
    check((sampling.sampled_cells <= sampling.cap).all(), "deterministic sampling cap respected", checks)
    check((sampling.sampled_cells <= sampling.available_cells).all(), "sampled cells available", checks)
    expected_groups = {"Malignant", "SPP1-state Macrophage", "C1QC-state Macrophage", "Other Macrophage", "Monocyte", "T cell", "NK cell", "Endothelial"}
    check(set(sampling.communication_group).issubset(expected_groups), "only preregistered communication groups", checks)

    meta_rows = {}
    for dataset in DATASETS:
        meta = pd.read_csv(root / f"prepared_inputs_v1/{dataset}/cellchat_cell_metadata_v1.tsv", sep="\t")
        check(meta.cell_id.notna().all() and meta.cell_id.is_unique, f"{dataset} unique nonmissing sampled cell IDs", checks)
        check(len(meta) == EXPECTED_CELLS[dataset], f"{dataset} metadata row count", checks)
        check(set(meta.communication_group).issubset(expected_groups), f"{dataset} metadata groups", checks)
        meta_rows[dataset] = len(meta)

    endpoint_parts, qc_parts, interaction_parts = [], [], []
    for dataset in DATASETS:
        run = root / "dataset_runs_v1" / dataset
        status = json.loads((run / f"{dataset}_status_v1.json").read_text(encoding="utf-8"))
        check(status["status"] == "PASS_PHASE88_CELLCHAT_DATASET", f"{dataset} formal run PASS", checks)
        check(status["units"] == EXPECTED_UNITS[dataset], f"{dataset} run unit count", checks)
        check(status["sampled_cells"] == EXPECTED_CELLS[dataset], f"{dataset} run sampled cell count", checks)
        check(status["interaction_rows"] == EXPECTED_INTERACTIONS[dataset], f"{dataset} interaction status count", checks)
        check((run / f"{dataset}_PASS_v1").is_file(), f"{dataset} PASS marker", checks)
        ep = pd.read_csv(run / f"{dataset}_unit_endpoints_v1.tsv", sep="\t")
        qc = pd.read_csv(run / f"{dataset}_unit_qc_v1.tsv", sep="\t")
        inter = pd.read_csv(run / f"{dataset}_cellchat_interactions_v1.tsv.gz", sep="\t")
        check(len(ep) == EXPECTED_UNITS[dataset] * 12, f"{dataset} 12-endpoint grid", checks)
        check(not ep.duplicated(["dataset", "unit_id", "endpoint_id"]).any(), f"{dataset} unique unit endpoints", checks)
        check(set(PRIMARY).issubset(set(ep.endpoint_id)), f"{dataset} six primary endpoints", checks)
        check(qc.status.eq("PASS_CELLCHAT_UNIT").all() and qc.unit_id.is_unique, f"{dataset} unit QC", checks)
        check(qc.sampled_cells.sum() == EXPECTED_CELLS[dataset], f"{dataset} unit sampled-cell sum", checks)
        check(len(inter) == EXPECTED_INTERACTIONS[dataset], f"{dataset} interaction row count", checks)
        check(np.isfinite(inter.prob.to_numpy(float)).all() and (inter.prob >= 0).all(), f"{dataset} finite nonnegative probabilities", checks)
        check(np.isfinite(inter.pval.to_numpy(float)).all() and inter.pval.between(0, 1).all(), f"{dataset} valid p-values", checks)
        primary_values = ep.loc[ep.endpoint_id.isin(PRIMARY), "endpoint_value"]
        check(np.isfinite(primary_values.dropna().to_numpy(float)).all(), f"{dataset} finite primary values where eligible", checks)
        if dataset != "GSE294109":
            check(primary_values.notna().all(), f"{dataset} complete primary grid", checks)
        endpoint_parts.append(ep)
        qc_parts.append(qc)
        interaction_parts.append(inter)

    units = pd.concat(endpoint_parts, ignore_index=True)
    qc_all = pd.concat(qc_parts, ignore_index=True)
    interactions = pd.concat(interaction_parts, ignore_index=True)
    check(qc_all[["dataset", "unit_id"]].drop_duplicates().shape[0] == 36, "36 eligible units", checks)
    check(len(units) == 432, "432 unit endpoint rows", checks)
    check(len(interactions) == 161865, "161865 combined interaction rows", checks)

    saved_units = pd.read_csv(result_dir / "phase88_unit_endpoints_v1.tsv", sep="\t")
    saved_qc = pd.read_csv(result_dir / "phase88_unit_qc_all_v1.tsv", sep="\t")
    check(saved_units.equals(units), "combined endpoint table exact", checks)
    check(saved_qc.equals(qc_all), "combined QC table exact", checks)

    patient = (units.groupby(["dataset", "patient_id", "endpoint_id", "tier"], sort=True)
               .agg(patient_malignant_tymp=("malignant_tymp_median", "median"),
                    patient_endpoint_value=("endpoint_value", "median"),
                    n_samples=("sample_id", "nunique")).reset_index())
    saved_patient = pd.read_csv(result_dir / "phase88_patient_endpoints_v1.tsv", sep="\t")
    keys = ["dataset", "patient_id", "endpoint_id", "tier"]
    patient = patient.sort_values(keys).reset_index(drop=True)
    saved_patient = saved_patient.sort_values(keys).reset_index(drop=True)
    check(patient[keys + ["n_samples"]].equals(saved_patient[keys + ["n_samples"]]), "patient aggregation keys exact", checks)
    check(equal_numeric(patient.patient_malignant_tymp, saved_patient.patient_malignant_tymp), "patient TYMP medians independently reproduced", checks)
    check(equal_numeric(patient.patient_endpoint_value, saved_patient.patient_endpoint_value), "patient endpoint medians independently reproduced", checks)
    check(patient[["dataset", "patient_id"]].drop_duplicates().shape[0] == 28, "28 cohort-specific patients", checks)

    rows = []
    for (dataset, endpoint, tier), frame in patient.groupby(["dataset", "endpoint_id", "tier"], sort=True):
        n, rho, p = correlation(frame)
        rows.append({"dataset": dataset, "endpoint_id": endpoint, "tier": tier, "n_patients": n, "spearman_rho": rho, "p_value": p})
    recalculated = pd.DataFrame(rows)
    recalculated["BH_q_within_dataset_primary6"] = np.nan
    for dataset, idx in recalculated.loc[recalculated.tier.eq("PRIMARY")].groupby("dataset").groups.items():
        recalculated.loc[idx, "BH_q_within_dataset_primary6"] = bh(recalculated.loc[idx, "p_value"])
    recalculated["formal_support"] = (recalculated.tier.eq("PRIMARY") & (recalculated.n_patients >= 10) &
                                       (recalculated.spearman_rho > 0) &
                                       (recalculated.BH_q_within_dataset_primary6 < 0.05))
    saved_summary = pd.read_csv(result_dir / "phase88_dataset_endpoint_summary_v1.tsv", sep="\t")
    sort_keys = ["dataset", "endpoint_id", "tier"]
    recalculated = recalculated.sort_values(sort_keys).reset_index(drop=True)
    saved_summary = saved_summary.sort_values(sort_keys).reset_index(drop=True)
    check(recalculated[sort_keys + ["n_patients"]].equals(saved_summary[sort_keys + ["n_patients"]]), "association row keys and n exact", checks)
    for column in ["spearman_rho", "p_value", "BH_q_within_dataset_primary6"]:
        check(equal_numeric(recalculated[column], saved_summary[column]), f"independent {column} reproduction", checks)
    check(recalculated.formal_support.equals(saved_summary.formal_support), "formal support independently reproduced", checks)
    check(not recalculated.loc[recalculated.tier.eq("PRIMARY"), "formal_support"].any(), "no primary formal support", checks)

    cross = pd.read_csv(result_dir / "phase88_cross_cohort_decisions_v1.tsv", sep="\t")
    check(set(cross.endpoint_id) == set(PRIMARY) and len(cross) == 6, "six cross-cohort decisions", checks)
    check(not cross.two_external_support.any() and not cross.three_cohort_support.any(), "no two-external or three-cohort support", checks)
    decision = json.loads((result_dir / "PHASE88_DECISION_v1.json").read_text(encoding="utf-8"))
    check(decision["decision"] == "NO_TWO_EXTERNAL_FINE_STATE_COMMUNICATION_SUPPORT", "scientific decision exact", checks)
    check(decision["patient_sample_units"] == 36 and decision["patients"] == 28, "decision counts exact", checks)
    check(decision["nichenet_rerun"] is False, "NicheNet not rerun recorded", checks)
    check(decision["preregistration_sha256"] == EXPECTED_PREREG and decision["cellchat_db_sha256"] == EXPECTED_DB, "decision frozen hashes exact", checks)

    manifest = pd.read_csv(result_dir / "FILE_MANIFEST_SHA256_v1.tsv", sep="\t")
    actual_files = sorted(p.name for p in result_dir.iterdir() if p.name != "FILE_MANIFEST_SHA256_v1.tsv")
    check(sorted(manifest.file.tolist()) == actual_files, "result manifest coverage exact", checks)
    for row in manifest.itertuples(index=False):
        path = result_dir / row.file
        check(path.stat().st_size == row.bytes and digest(path) == row.sha256, f"manifest exact: {row.file}", checks)

    forbidden = ("response", "benefit", "pfs", "os", "treatment")
    tables = [units, qc_all, interactions, patient, saved_summary, cross, sampling, input_status]
    check(not any(any(token in str(col).lower() for token in forbidden) for table in tables for col in table.columns), "no outcome/treatment columns", checks)

    pd.DataFrame(checks).to_csv(part / "PHASE88_VALIDATION_CHECKS_v1.tsv", sep="\t", index=False)
    validation = {
        "status": "PASS_PHASE88_INDEPENDENT_VALIDATION",
        "checks": len(checks),
        "decision": decision["decision"],
        "units": 36,
        "patients": 28,
        "interaction_rows": 161865,
        "preregistration_sha256": EXPECTED_PREREG,
        "cellchat_db_sha256": EXPECTED_DB,
    }
    (part / "PHASE88_VALIDATION_STATUS_v1.json").write_text(json.dumps(validation, indent=2) + "\n", encoding="utf-8")
    (part / "PHASE88_VALIDATION_PASS_v1").write_text("PASS_PHASE88_INDEPENDENT_VALIDATION\n", encoding="utf-8")
    os.replace(part, final)
    print(json.dumps(validation))


if __name__ == "__main__":
    main()
