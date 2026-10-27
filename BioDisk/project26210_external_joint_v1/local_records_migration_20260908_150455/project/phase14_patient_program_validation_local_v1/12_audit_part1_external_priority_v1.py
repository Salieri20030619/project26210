#!/usr/bin/env python3
"""Post-hoc audit of Part1-rebuild priority versus fresh external CNV calls."""

from __future__ import annotations

import json
import os
from datetime import datetime
from pathlib import Path

import pandas as pd


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
SOURCE = ROOT / "phase11_external_cnv_consensus_local_v1" / "final_server_sync_v5" / "external_cnv_cell_consensus_final_v2.tsv.gz"
OUT = HERE / "posthoc_part1_external_priority_audit_v1"


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


def main() -> None:
    if OUT.exists():
        raise FileExistsError(OUT)
    OUT.mkdir(parents=True)
    columns = [
        "dataset", "patient_key", "patient_id", "phase7_fine_annotation",
        "phase9_part1_cnv_assessed", "phase9_malignancy_supported", "phase9_cnv_gate_status",
        "external_malignant_conservative_incremental_v1", "external_cnv_status_incremental_v1",
    ]
    frame = pd.read_csv(SOURCE, sep="\t", usecols=columns, low_memory=False)
    tumor = frame[frame["phase7_fine_annotation"] == "Epithelial_tumor_like"].copy()
    assessed = tumor[tumor["phase9_part1_cnv_assessed"].fillna(False).astype(bool)].copy()
    assessed["part1_supported"] = assessed["phase9_malignancy_supported"].fillna(False).astype(bool)
    assessed["external_supported"] = assessed["external_malignant_conservative_incremental_v1"].fillna(False).astype(bool)
    assessed["priority_comparison"] = "agree_not_supported"
    assessed.loc[assessed["part1_supported"] & assessed["external_supported"], "priority_comparison"] = "agree_supported"
    assessed.loc[assessed["part1_supported"] & ~assessed["external_supported"], "priority_comparison"] = "part1_supported_external_not"
    assessed.loc[~assessed["part1_supported"] & assessed["external_supported"], "priority_comparison"] = "part1_not_external_supported"

    dataset = (
        assessed.groupby(["dataset", "priority_comparison", "phase9_cnv_gate_status", "external_cnv_status_incremental_v1"], dropna=False)
        .size().reset_index(name="n_cells")
        .sort_values(["dataset", "n_cells"], ascending=[True, False])
    )
    atomic_tsv(dataset, OUT / "dataset_part1_external_priority_v1.tsv")
    patient = (
        assessed.groupby(["dataset", "patient_key", "patient_id", "priority_comparison"], dropna=False)
        .size().reset_index(name="n_cells")
        .sort_values(["dataset", "patient_key", "n_cells"], ascending=[True, True, False])
    )
    atomic_tsv(patient, OUT / "patient_part1_external_priority_v1.tsv")
    counts = assessed["priority_comparison"].value_counts().to_dict()
    status = {
        "schema_version": "phase14_part1_external_priority_audit_v1",
        "status": "PASS_POSTHOC_AUDIT",
        "completed_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "tumor_like_cells_in_external_table": int(len(tumor)),
        "part1_assessed_tumor_like_cells": int(len(assessed)),
        "comparison_counts": {str(k): int(v) for k, v in counts.items()},
        "part1_priority_excludes_external_supported_cells": int(counts.get("part1_not_external_supported", 0)),
        "part1_priority_includes_external_not_supported_cells": int(counts.get("part1_supported_external_not", 0)),
        "interpretation": "Priority disagreements are an audit of two runs, not evidence that either method is ground truth; final v6 intentionally follows exact-cell Part1 rebuild calls.",
    }
    atomic_json(status, OUT / "PART1_EXTERNAL_PRIORITY_STATUS_v1.json")


if __name__ == "__main__":
    main()

