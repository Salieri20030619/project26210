#!/usr/bin/env python3
"""Cell- and patient-level decomposition of the final conservative CNV gate."""

from __future__ import annotations

import json
import os
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
SYNC = ROOT / "phase11_external_cnv_consensus_local_v1" / "final_server_sync_v5"
OUT = HERE / "posthoc_cnv_gate_attrition_v2"
SOURCE = SYNC / "external_cnv_cell_consensus_final_v2.tsv.gz"
PART1 = SYNC / "phase9_all_cells_cnv_gate_v1.tsv.gz"


def atomic_tsv(frame: pd.DataFrame, path: Path) -> None:
    part = path.with_name(path.name + ".part")
    if path.exists() or part.exists():
        raise FileExistsError(path)
    frame.to_csv(part, sep="\t", index=False, encoding="utf-8", na_rep="NA", lineterminator="\n")
    os.replace(part, path)


def atomic_json(payload: dict, path: Path) -> None:
    part = path.with_name(path.name + ".part")
    if path.exists() or part.exists():
        raise FileExistsError(path)
    part.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    os.replace(part, path)


def true(series: pd.Series) -> pd.Series:
    return series.fillna(False).astype(bool)


def main() -> None:
    if OUT.exists():
        raise FileExistsError(OUT)
    OUT.mkdir(parents=True)

    external_columns = [
        "exact_cell_id",
        "infercnv_eligible", "copykat_eligible", "infercnv_call", "copykat_call",
        "external_cnv_status_incremental_v1", "external_malignant_conservative_incremental_v1",
        "external_cnv_evidence_incremental_v1",
    ]
    part1_columns = [
        "cell_id", "dataset", "patient_id", "patient_key", "phase7_fine_annotation",
        "phase9_part1_cnv_assessed", "phase9_cnv_gate_status", "phase9_malignancy_supported",
    ]
    external = pd.read_csv(SOURCE, sep="\t", usecols=external_columns, low_memory=False)
    part1 = pd.read_csv(PART1, sep="\t", usecols=part1_columns, low_memory=False).rename(columns={"cell_id": "exact_cell_id"})
    tumor = part1[part1["phase7_fine_annotation"] == "Epithelial_tumor_like"].merge(
        external, on="exact_cell_id", how="left", validate="one_to_one", indicator="external_scope"
    )
    if len(tumor) != 40302 or tumor["exact_cell_id"].duplicated().any():
        raise RuntimeError(f"Tumor-like conservation failure: rows={len(tumor)}")

    part1_assessed = true(tumor["phase9_part1_cnv_assessed"])
    part1_supported = true(tumor["phase9_malignancy_supported"])
    external_supported = true(tumor["external_malignant_conservative_incremental_v1"])
    tumor["final_supported_recomputed"] = part1_supported.where(part1_assessed, external_supported)
    tumor["gate_source"] = np.select(
        [
            part1_assessed & part1_supported,
            part1_assessed & ~part1_supported,
            ~part1_assessed & external_supported,
            ~part1_assessed & ~external_supported & (tumor["external_scope"] == "both"),
            ~part1_assessed & (tumor["external_scope"] != "both"),
        ],
        [
            "part1_supported",
            "part1_not_conservative_malignant",
            "external_supported",
            "external_not_conservative_malignant",
            "external_not_in_scope",
        ],
        default="unclassified",
    )
    tumor["gate_detail"] = np.where(
        part1_assessed,
        tumor["phase9_cnv_gate_status"].fillna("part1_status_missing"),
        tumor["external_cnv_status_incremental_v1"].fillna("external_not_in_scope"),
    )
    tumor["method_pair"] = tumor["infercnv_call"].fillna("NA") + "__" + tumor["copykat_call"].fillna("NA")

    patient = (
        tumor.groupby(["dataset", "patient_key", "patient_id", "gate_source", "gate_detail", "method_pair", "external_cnv_evidence_incremental_v1"], dropna=False)
        .size()
        .reset_index(name="n_cells")
        .sort_values(["dataset", "patient_key", "n_cells"], ascending=[True, True, False])
    )
    atomic_tsv(patient, OUT / "patient_gate_evidence_breakdown_v1.tsv")

    dataset = (
        tumor.groupby(["dataset", "gate_source", "gate_detail"], dropna=False)
        .size()
        .reset_index(name="n_cells")
        .sort_values(["dataset", "n_cells"], ascending=[True, False])
    )
    atomic_tsv(dataset, OUT / "dataset_gate_status_breakdown_v1.tsv")

    method_pairs = (
        tumor.groupby(["dataset", "method_pair", "external_cnv_status_incremental_v1"], dropna=False)
        .size()
        .reset_index(name="n_cells")
        .sort_values(["dataset", "n_cells"], ascending=[True, False])
    )
    atomic_tsv(method_pairs, OUT / "dataset_method_pair_breakdown_v1.tsv")

    patient_totals = (
        tumor.groupby(["dataset", "patient_key", "patient_id"], as_index=False)
        .agg(
            provisional_tumor_like_cells=("exact_cell_id", "size"),
            final_supported_cells=("final_supported_recomputed", "sum"),
            part1_assessed_cells=("phase9_part1_cnv_assessed", lambda x: int(true(x).sum())),
            external_supported_cells=("external_malignant_conservative_incremental_v1", lambda x: int(true(x).sum())),
            infercnv_high_cells=("infercnv_call", lambda x: int((x == "high_cnv").sum())),
            copykat_aneuploid_cells=("copykat_call", lambda x: int((x == "aneuploid").sum())),
            malignant_consensus_cells=("external_cnv_status_incremental_v1", lambda x: int((x == "malignant_consensus").sum())),
            method_discordant_cells=("external_cnv_status_incremental_v1", lambda x: int((x == "method_discordant").sum())),
            non_malignant_consensus_cells=("external_cnv_status_incremental_v1", lambda x: int((x == "non_malignant_consensus").sum())),
        )
    )
    patient_totals["final_retention_fraction"] = patient_totals["final_supported_cells"] / patient_totals["provisional_tumor_like_cells"]
    patient_totals["method_discordant_fraction"] = patient_totals["method_discordant_cells"] / patient_totals["provisional_tumor_like_cells"]
    patient_totals["non_malignant_consensus_fraction"] = patient_totals["non_malignant_consensus_cells"] / patient_totals["provisional_tumor_like_cells"]
    atomic_tsv(patient_totals.sort_values(["dataset", "patient_key"]), OUT / "patient_gate_summary_v1.tsv")

    status = {
        "schema_version": "phase14_cnv_gate_attrition_diagnostic_v2",
        "status": "PASS_POSTHOC_DIAGNOSTIC",
        "completed_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "external_rows": int(len(external)),
        "tumor_like_rows": int(len(tumor)),
        "final_supported_rows": int(tumor["final_supported_recomputed"].sum()),
        "tumor_like_in_external_scope": int((tumor["external_scope"] == "both").sum()),
        "tumor_like_not_in_external_scope": int((tumor["external_scope"] != "both").sum()),
        "gate_source_counts": {str(k): int(v) for k, v in tumor["gate_source"].value_counts().items()},
        "external_status_counts": {str(k): int(v) for k, v in tumor["external_cnv_status_incremental_v1"].fillna("NA").value_counts().items()},
        "interpretation": "Conservative support requires both inferCNV-high and CopyKAT-aneuploid when both are eligible; discordance is not equivalent to proven benign biology.",
    }
    atomic_json(status, OUT / "CNV_GATE_ATTRITION_STATUS_v1.json")


if __name__ == "__main__":
    main()
