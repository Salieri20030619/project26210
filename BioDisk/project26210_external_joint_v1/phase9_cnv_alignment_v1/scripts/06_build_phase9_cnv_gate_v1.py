#!/usr/bin/env python3
"""Build a conservative all-cell CNV evidence gate without transferring labels to unmatched cells."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def atomic_tsv(frame: pd.DataFrame, final: Path, compression: str | None = None) -> None:
    if final.exists():
        raise FileExistsError(f"No-overwrite gate: {final}")
    part = final.with_name(final.name + ".part")
    if part.exists():
        raise FileExistsError(f"Incomplete file exists: {part}")
    frame.to_csv(part, sep="\t", index=False, na_rep="NA", compression=compression)
    os.replace(part, final)


def atomic_text(text: str, final: Path) -> None:
    if final.exists():
        raise FileExistsError(f"No-overwrite gate: {final}")
    part = final.with_name(final.name + ".part")
    part.write_text(text, encoding="utf-8")
    os.replace(part, final)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--crosswalk", required=True)
    parser.add_argument("--part1-map", required=True)
    parser.add_argument("--out-dir", required=True)
    args = parser.parse_args()

    crosswalk_path = Path(args.crosswalk).resolve()
    map_path = Path(args.part1_map).resolve()
    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    outputs = {
        "all_cells": out_dir / "phase9_all_cells_cnv_gate_v1.tsv.gz",
        "gate_summary": out_dir / "phase9_cnv_gate_summary_v1.tsv",
        "patient_summary": out_dir / "phase9_cnv_gate_by_patient_v1.tsv",
        "f1_summary": out_dir / "phase9_f1_cnv_gate_summary_v1.tsv",
        "status": out_dir / "PHASE9_CNV_GATE_STATUS",
        "metadata": out_dir / "run_metadata_v1.json",
    }
    if any(path.exists() for path in outputs.values()):
        raise FileExistsError("No-overwrite gate: one or more outputs exist")

    ext_columns = [
        "cell_id", "dataset", "object_name", "sample_id", "patient_id", "patient_key",
        "barcode_original", "src__gsm", "phase7_broad_annotation", "phase7_fine_annotation",
        "phase7_annotation_confidence", "phase8_spatial_reference_state",
        "phase8_spatial_display_label", "phase8_spatial_mapping_status",
        "phase8_spatial_annotation_confidence", "phase8_f1_projection_score",
        "phase8_f1_percentile_patient", "phase8_f1_eligible_patient",
    ]
    ext = pd.read_csv(crosswalk_path, sep="\t", dtype=str, compression="infer", usecols=ext_columns)
    direct = pd.read_csv(map_path, sep="\t", dtype=str, compression="infer")
    if len(ext) != 585119 or ext["cell_id"].duplicated().any():
        raise RuntimeError("External crosswalk count/uniqueness gate failed")
    if len(direct) != 137338 or direct["cell_id"].duplicated().any():
        raise RuntimeError("Part 1 exact map count/uniqueness gate failed")

    cnv_columns = [
        "cell_id", "part1__part1_obs_name", "part1__study", "part1__curated_patient_id",
        "part1__geo_accession", "part1__cell_type_broad_v1", "part1__cnv_assessed_v1",
        "part1__cnv_analysis_scope_v1", "part1__infercnv_call_v1",
        "part1__copykat_call_v1", "part1__malignancy_status_v1",
        "part1__malignant_conservative_v1", "part1__malignancy_evidence_v1",
        "part1__eligible_for_subclone_selection_v1",
    ]
    missing = set(cnv_columns).difference(direct.columns)
    if missing:
        raise RuntimeError(f"Exact map missing fields: {sorted(missing)}")
    merged = ext.merge(direct[cnv_columns], on="cell_id", how="left", validate="one_to_one")

    matched = merged["part1__part1_obs_name"].notna()
    assessed = merged["part1__cnv_assessed_v1"].str.lower().eq("true")
    malignant = merged["part1__malignant_conservative_v1"].str.lower().eq("true")
    epithelial = merged["phase7_broad_annotation"].eq("Epithelial")
    status = merged["part1__malignancy_status_v1"]

    gate = pd.Series("cnv_unavailable_no_identical_part1_cell", index=merged.index, dtype="string")
    gate.loc[matched & ~assessed] = "cnv_not_assessed_part1"
    gate.loc[matched & assessed & ~epithelial] = "annotation_cnv_scope_conflict_review"
    gate.loc[matched & assessed & epithelial & malignant] = "malignant_supported_part1"
    gate.loc[matched & assessed & epithelial & ~malignant & status.eq("non_malignant_consensus")] = "non_malignant_consensus_part1"
    gate.loc[matched & assessed & epithelial & ~malignant & status.eq("method_discordant")] = "cnv_method_discordant_part1"
    gate.loc[matched & assessed & epithelial & ~malignant & status.eq("indeterminate")] = "cnv_indeterminate_part1"
    gate.loc[matched & assessed & epithelial & ~malignant & ~status.isin(["non_malignant_consensus", "method_discordant", "indeterminate"])] = "cnv_assessed_other_part1"

    merged["phase9_part1_identical_cell"] = matched
    merged["phase9_part1_cnv_assessed"] = matched & assessed
    merged["phase9_current_epithelial"] = epithelial
    merged["phase9_cnv_gate_status"] = gate
    merged["phase9_malignancy_supported"] = gate.eq("malignant_supported_part1")
    merged["phase9_non_malignant_supported"] = gate.eq("non_malignant_consensus_part1")
    merged["phase9_requires_annotation_review"] = gate.eq("annotation_cnv_scope_conflict_review")
    merged["phase9_cnv_gate_version"] = "phase9_cnv_gate_v1"

    gate_summary = (
        merged.groupby(
            ["dataset", "phase7_broad_annotation", "phase8_spatial_reference_state", "phase9_cnv_gate_status"],
            observed=True,
            dropna=False,
        )
        .size()
        .reset_index(name="n_cells")
    )
    patient_summary = (
        merged.groupby(["dataset", "patient_id", "phase9_cnv_gate_status"], observed=True, dropna=False)
        .size()
        .reset_index(name="n_cells")
    )
    f1 = merged.loc[merged["phase8_spatial_reference_state"].astype(str).str.startswith("F1-")]
    f1_summary = (
        f1.groupby(["dataset", "phase8_spatial_reference_state", "phase9_cnv_gate_status"], observed=True, dropna=False)
        .size()
        .reset_index(name="n_cells")
    )

    keep = ext_columns + [column for column in cnv_columns if column != "cell_id"] + [
        "phase9_part1_identical_cell", "phase9_part1_cnv_assessed", "phase9_current_epithelial",
        "phase9_cnv_gate_status", "phase9_malignancy_supported",
        "phase9_non_malignant_supported", "phase9_requires_annotation_review",
        "phase9_cnv_gate_version",
    ]
    atomic_tsv(merged[keep], outputs["all_cells"], compression="gzip")
    atomic_tsv(gate_summary, outputs["gate_summary"])
    atomic_tsv(patient_summary, outputs["patient_summary"])
    atomic_tsv(f1_summary, outputs["f1_summary"])

    counts = {
        "all_cells": int(len(merged)),
        "identical_part1_cells": int(matched.sum()),
        "part1_cnv_assessed": int((matched & assessed).sum()),
        "current_epithelial_part1_cnv_assessed": int((matched & assessed & epithelial).sum()),
        "malignancy_supported": int(merged["phase9_malignancy_supported"].sum()),
        "non_malignant_supported": int(merged["phase9_non_malignant_supported"].sum()),
        "annotation_cnv_conflict_review": int(merged["phase9_requires_annotation_review"].sum()),
        "f1_cells": int(len(f1)),
        "f1_identical_part1": int(f1["phase9_part1_identical_cell"].sum()),
        "f1_part1_cnv_assessed": int(f1["phase9_part1_cnv_assessed"].sum()),
        "f1_malignancy_supported": int(f1["phase9_malignancy_supported"].sum()),
        "f1_non_malignant_supported": int(f1["phase9_non_malignant_supported"].sum()),
    }
    metadata = {
        "status": "PASS",
        "completed_at_utc": datetime.now(timezone.utc).isoformat(),
        "policy": "Part 1 CNV evidence is transferred only for identical cells. Malignancy support additionally requires the current Phase 7 broad class to be Epithelial. Non-epithelial assessed matches are review conflicts, not malignant calls.",
        "inputs": {
            "crosswalk": {"path": str(crosswalk_path), "bytes": crosswalk_path.stat().st_size, "sha256": sha256(crosswalk_path)},
            "part1_map": {"path": str(map_path), "bytes": map_path.stat().st_size, "sha256": sha256(map_path)},
        },
        "counts": counts,
        "outputs": {
            key: {"path": str(outputs[key]), "bytes": outputs[key].stat().st_size, "sha256": sha256(outputs[key])}
            for key in ("all_cells", "gate_summary", "patient_summary", "f1_summary")
        },
    }
    atomic_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n", outputs["metadata"])
    status_lines = ["PASS"] + [f"{key}={value}" for key, value in counts.items()]
    status_lines += [f"{key}_sha256={value['sha256']}" for key, value in metadata["outputs"].items()]
    atomic_text("\n".join(status_lines) + "\n", outputs["status"])
    print(json.dumps(counts, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
