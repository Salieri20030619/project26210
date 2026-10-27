#!/usr/bin/env python3
"""Validate terminal Phase 10 CNV units and build a conservative incremental snapshot."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


NEUTRAL_STATE = 3


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def scaled_mad(values: np.ndarray) -> float:
    median = float(np.median(values))
    return float(np.median(np.abs(values - median)) * 1.4826)


def atomic_frame(frame: pd.DataFrame, path: Path) -> None:
    part = path.with_name(path.name + ".part")
    compression = "gzip" if path.suffix == ".gz" else None
    frame.to_csv(part, sep="\t", index=False, compression=compression, lineterminator="\n")
    os.replace(part, path)


def atomic_json(path: Path, payload: dict[str, object]) -> None:
    part = path.with_name(path.name + ".part")
    part.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(part, path)


def read_mapping(path: Path) -> pd.DataFrame:
    mapping = pd.read_csv(path, sep="\t", compression="gzip", dtype=str)
    required = {
        "tool_id", "exact_cell_id", "role", "sample_id", "curated_patient_id",
        "study", "source_tissue", "cell_type_broad_v1",
    }
    if not required.issubset(mapping.columns):
        raise RuntimeError(f"Missing mapping columns: {path}")
    if not mapping["tool_id"].is_unique or not mapping["exact_cell_id"].is_unique:
        raise RuntimeError(f"Non-unique mapping IDs: {path}")
    if not set(mapping["role"]).issubset({"query", "reference"}):
        raise RuntimeError(f"Unexpected mapping role: {path}")
    return mapping


def infer_terminal_calls(
    unit: str, input_dir: Path, output_dir: Path
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, object]]:
    source_mapping = read_mapping(input_dir / "cell_id_mapping_v1.tsv.gz")
    output_mapping = pd.read_csv(output_dir / "cell_id_mapping_v2.tsv", sep="\t", dtype=str)
    if not source_mapping["tool_id"].equals(output_mapping["tool_id"]):
        raise RuntimeError(f"inferCNV tool ID order changed: {unit}")
    if not source_mapping["exact_cell_id"].equals(output_mapping["exact_cell_id"]):
        raise RuntimeError(f"inferCNV exact ID order changed: {unit}")
    run_dir = output_dir / "infercnv_run"
    groups_path = run_dir / "17_HMM_predHMMi6.leiden.hmm_mode-subclusters.cell_groupings"
    genes_path = run_dir / "17_HMM_predHMMi6.leiden.hmm_mode-subclusters.pred_cnv_genes.dat"
    post_rds = output_dir / "infercnv_object_post_run_v2.rds"
    for path in (groups_path, genes_path, post_rds, output_dir / "status.json"):
        if not path.is_file() or path.stat().st_size <= 0:
            raise RuntimeError(f"Missing terminal inferCNV output: {path}")
    status = json.loads((output_dir / "status.json").read_text(encoding="utf-8"))
    if status.get("status") != "completed":
        raise RuntimeError(f"Non-completed inferCNV status: {unit}")

    groups = pd.read_csv(groups_path, sep="\t", dtype=str).rename(columns={"cell": "tool_id"})
    groups = groups.merge(
        source_mapping[
            ["tool_id", "exact_cell_id", "role", "sample_id", "curated_patient_id", "study"]
        ],
        on="tool_id",
        how="left",
        validate="one_to_one",
    )
    if groups["exact_cell_id"].isna().any() or set(groups["tool_id"]) != set(source_mapping["tool_id"]):
        raise RuntimeError(f"Incomplete inferCNV group mapping: {unit}")
    role_counts = groups.groupby("cell_group_name", observed=True)["role"].nunique()
    if not role_counts.eq(1).all():
        raise RuntimeError(f"inferCNV group mixes roles: {unit}")

    genes = pd.read_csv(genes_path, sep="\t", dtype={"cell_group_name": str, "state": int})
    genes["abs_state_shift"] = (genes["state"] - NEUTRAL_STATE).abs()
    burden = genes.groupby("cell_group_name", observed=True).agg(
        n_cnv_genes=("gene", "size"),
        n_gain_genes=("state", lambda x: int((x > NEUTRAL_STATE).sum())),
        n_loss_genes=("state", lambda x: int((x < NEUTRAL_STATE).sum())),
        cnv_gene_burden=("abs_state_shift", "sum"),
        n_chromosomes=("chr", "nunique"),
    ).reset_index()
    roles = groups.groupby("cell_group_name", observed=True).agg(
        role=("role", "first"), n_cells=("tool_id", "size")
    ).reset_index()
    burden = roles.merge(burden, on="cell_group_name", how="left", validate="one_to_one").fillna(0)
    ref = burden.loc[burden["role"].eq("reference"), "cnv_gene_burden"].to_numpy(float)
    if ref.size == 0:
        raise RuntimeError(f"No inferCNV reference groups: {unit}")
    ref_median = float(np.median(ref))
    ref_mad = scaled_mad(ref)
    ref_q95 = float(np.quantile(ref, 0.95))
    threshold = max(ref_median + 3 * ref_mad, ref_q95)
    burden["infercnv_call"] = np.where(
        burden["role"].eq("reference"),
        "reference",
        np.where(burden["cnv_gene_burden"].gt(threshold), "high_cnv", "low_cnv"),
    )
    burden.insert(0, "infercnv_unit", unit)
    burden["high_cnv_threshold"] = threshold
    calls = groups.merge(
        burden[["cell_group_name", "cnv_gene_burden", "infercnv_call", "high_cnv_threshold"]],
        on="cell_group_name",
        how="left",
        validate="many_to_one",
    )
    calls = calls.loc[calls["role"].eq("query")].copy()
    calls.insert(0, "infercnv_unit", unit)
    if calls["exact_cell_id"].duplicated().any() or calls["infercnv_call"].isna().any():
        raise RuntimeError(f"Invalid inferCNV query calls: {unit}")
    metadata = {
        "unit_type": "infercnv",
        "unit_id": unit,
        "run_state": "completed",
        "n_query": int(len(calls)),
        "n_reference": int(source_mapping["role"].eq("reference").sum()),
        "n_genes_used": int(status.get("n_genes_used", -1)),
        "gene_coverage_qc": str(status.get("gene_coverage_qc", "missing")),
        "high_cnv_threshold": threshold,
        "rds_size_bytes": post_rds.stat().st_size,
        "validation_status": "PASS",
    }
    return calls, burden, metadata


def copy_terminal_calls(
    unit: str, input_dir: Path, output_dir: Path, run_state: str
) -> tuple[pd.DataFrame, dict[str, object]]:
    source_mapping = read_mapping(input_dir / "cell_id_mapping_v1.tsv.gz")
    query = source_mapping.loc[source_mapping["role"].eq("query")].copy()
    if run_state == "completed":
        prediction_path = output_dir / "copykat_predictions_mapped_v2.tsv"
        result_rds = output_dir / "copykat_result_v2.rds"
        for path in (prediction_path, result_rds, output_dir / "status.json"):
            if not path.is_file() or path.stat().st_size <= 0:
                raise RuntimeError(f"Missing terminal CopyKAT output: {path}")
        status = json.loads((output_dir / "status.json").read_text(encoding="utf-8"))
        if status.get("status") != "completed":
            raise RuntimeError(f"Non-completed CopyKAT status: {unit}")
        calls = pd.read_csv(prediction_path, sep="\t", dtype=str)
        calls = calls.loc[calls["role"].eq("query")].copy()
        raw = calls["copykat.pred"].fillna("")
        if not set(raw).issubset({"aneuploid", "diploid", "not.defined"}):
            raise RuntimeError(f"Unexpected CopyKAT calls: {unit}")
        calls["copykat_call"] = raw.map(
            {"aneuploid": "aneuploid", "diploid": "diploid", "not.defined": "not_defined"}
        )
        if calls["exact_cell_id"].duplicated().any() or not set(calls["exact_cell_id"]).issubset(set(query["exact_cell_id"])):
            raise RuntimeError(f"Invalid CopyKAT query mapping: {unit}")
        rds_size = result_rds.stat().st_size
        n_genes = int(status.get("n_genes_used", -1))
        gene_qc = str(status.get("gene_coverage_qc", "missing"))
    elif run_state == "not_estimable":
        calls = query.copy()
        calls["copykat.pred"] = pd.NA
        calls["copykat_call"] = "not_estimable"
        rds_size, n_genes, gene_qc = 0, -1, "not_estimable"
    else:
        raise RuntimeError(f"Unsupported CopyKAT terminal state: {run_state}")
    keep = [
        "exact_cell_id", "sample_id", "curated_patient_id", "study",
        "copykat.pred", "copykat_call",
    ]
    calls = calls[keep]
    calls.insert(0, "copykat_unit", unit)
    metadata = {
        "unit_type": "copykat",
        "unit_id": unit,
        "run_state": run_state,
        "n_query": int(len(calls)),
        "n_reference": int(source_mapping["role"].eq("reference").sum()),
        "n_genes_used": n_genes,
        "gene_coverage_qc": gene_qc,
        "high_cnv_threshold": np.nan,
        "rds_size_bytes": rds_size,
        "validation_status": "PASS",
    }
    return calls, metadata


def truth_table(infer_call: str, copy_call: str, infer_pending: bool, copy_pending: bool) -> tuple[str, bool, str, bool]:
    if infer_pending or copy_pending:
        return "pending_other_method", False, "eligible_method_not_terminal", False
    if infer_call == "high_cnv" and copy_call == "aneuploid":
        return "malignant_consensus", True, "infercnv_high_and_copykat_aneuploid", True
    if infer_call == "high_cnv" and copy_call in {"not_defined", "not_estimable", "not_eligible"}:
        return "malignant_infercnv_only", True, "infercnv_high_copykat_unavailable", True
    if infer_call == "not_eligible" and copy_call == "aneuploid":
        return "malignant_copykat_only", True, "copykat_aneuploid_infercnv_unavailable", True
    if infer_call == "low_cnv" and copy_call == "diploid":
        return "non_malignant_consensus", False, "infercnv_low_and_copykat_diploid", True
    if (infer_call == "high_cnv" and copy_call == "diploid") or (infer_call == "low_cnv" and copy_call == "aneuploid"):
        return "method_discordant", False, "infercnv_copykat_direction_conflict", True
    return "indeterminate", False, "insufficient_or_single_low_cnv_evidence", True


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--phase10-root", required=True)
    parser.add_argument("--phase9-gate", required=True)
    parser.add_argument("--out-dir", required=True)
    args = parser.parse_args()
    phase10 = Path(args.phase10_root).resolve()
    inputs = phase10 / "full_inputs_v2"
    outputs = phase10 / "results/full_v2"
    phase9_gate = Path(args.phase9_gate).resolve()
    out = Path(args.out_dir).resolve()
    if out.exists():
        raise FileExistsError(f"No-overwrite gate: {out}")
    out.mkdir(parents=True)

    validations: list[dict[str, object]] = []
    infer_calls: list[pd.DataFrame] = []
    infer_burdens: list[pd.DataFrame] = []
    copy_calls: list[pd.DataFrame] = []
    infer_eligible: list[pd.DataFrame] = []
    copy_eligible: list[pd.DataFrame] = []

    for unit_type, expected in (("infercnv", 49), ("copykat", 62)):
        unit_dirs = sorted(path for path in (inputs / unit_type).iterdir() if path.is_dir())
        if len(unit_dirs) != expected:
            raise RuntimeError(f"{unit_type} input count {len(unit_dirs)} != {expected}")
        for input_dir in unit_dirs:
            unit = input_dir.name
            mapping = read_mapping(input_dir / "cell_id_mapping_v1.tsv.gz")
            query = mapping.loc[mapping["role"].eq("query"), [
                "exact_cell_id", "sample_id", "curated_patient_id", "study"
            ]].copy()
            query[f"{unit_type}_unit"] = unit
            (infer_eligible if unit_type == "infercnv" else copy_eligible).append(query)
            output_dir = outputs / unit_type / unit
            if (output_dir / "RUN_COMPLETE").is_file():
                run_state = "completed"
            elif unit_type == "copykat" and (output_dir / "RUN_NOT_ESTIMABLE").is_file():
                run_state = "not_estimable"
            elif (output_dir / "ERROR.txt").is_file():
                run_state = "failed"
            elif output_dir.is_dir():
                run_state = "running"
            else:
                run_state = "pending"
            base_record = {
                "unit_type": unit_type,
                "unit_id": unit,
                "run_state": run_state,
                "n_query": int(len(query)),
                "n_reference": int(mapping["role"].eq("reference").sum()),
                "n_genes_used": np.nan,
                "gene_coverage_qc": "not_terminal",
                "high_cnv_threshold": np.nan,
                "rds_size_bytes": 0,
                "validation_status": "PENDING" if run_state in {"pending", "running"} else "FAIL",
            }
            if unit_type == "infercnv" and run_state == "completed":
                calls, burdens, record = infer_terminal_calls(unit, input_dir, output_dir)
                infer_calls.append(calls)
                infer_burdens.append(burdens)
                validations.append(record)
            elif unit_type == "copykat" and run_state in {"completed", "not_estimable"}:
                calls, record = copy_terminal_calls(unit, input_dir, output_dir, run_state)
                copy_calls.append(calls)
                validations.append(record)
            else:
                validations.append(base_record)

    infer_eligible_df = pd.concat(infer_eligible, ignore_index=True)
    copy_eligible_df = pd.concat(copy_eligible, ignore_index=True)
    if infer_eligible_df["exact_cell_id"].duplicated().any() or copy_eligible_df["exact_cell_id"].duplicated().any():
        raise RuntimeError("Eligible query cells duplicate within a method")
    infer_calls_df = pd.concat(infer_calls, ignore_index=True) if infer_calls else pd.DataFrame()
    copy_calls_df = pd.concat(copy_calls, ignore_index=True) if copy_calls else pd.DataFrame()
    infer_burden_df = pd.concat(infer_burdens, ignore_index=True) if infer_burdens else pd.DataFrame()
    if not infer_calls_df.empty and infer_calls_df["exact_cell_id"].duplicated().any():
        raise RuntimeError("Terminal inferCNV query calls duplicate across units")
    if not copy_calls_df.empty and copy_calls_df["exact_cell_id"].duplicated().any():
        raise RuntimeError("Terminal CopyKAT query calls duplicate across units")

    registry = infer_eligible_df.merge(
        copy_eligible_df[["exact_cell_id", "copykat_unit"]],
        on="exact_cell_id",
        how="outer",
        validate="one_to_one",
    )
    registry["infercnv_eligible"] = registry["infercnv_unit"].notna()
    registry["copykat_eligible"] = registry["copykat_unit"].notna()
    validation_df = pd.DataFrame(validations)
    state_map = validation_df.set_index(["unit_type", "unit_id"])["run_state"]
    registry["infercnv_run_state"] = [
        state_map.get(("infercnv", unit), "not_eligible") if pd.notna(unit) else "not_eligible"
        for unit in registry["infercnv_unit"]
    ]
    registry["copykat_run_state"] = [
        state_map.get(("copykat", unit), "not_eligible") if pd.notna(unit) else "not_eligible"
        for unit in registry["copykat_unit"]
    ]
    if not infer_calls_df.empty:
        infer_map = infer_calls_df.set_index("exact_cell_id")
        registry["infercnv_call"] = registry["exact_cell_id"].map(infer_map["infercnv_call"])
        registry["infercnv_group"] = registry["exact_cell_id"].map(infer_map["cell_group_name"])
        registry["infercnv_cnv_gene_burden"] = registry["exact_cell_id"].map(infer_map["cnv_gene_burden"])
        registry["infercnv_high_cnv_threshold"] = registry["exact_cell_id"].map(infer_map["high_cnv_threshold"])
    else:
        registry[["infercnv_call", "infercnv_group", "infercnv_cnv_gene_burden", "infercnv_high_cnv_threshold"]] = pd.NA
    if not copy_calls_df.empty:
        copy_map = copy_calls_df.set_index("exact_cell_id")
        registry["copykat_call"] = registry["exact_cell_id"].map(copy_map["copykat_call"])
        registry["copykat_raw_call"] = registry["exact_cell_id"].map(copy_map["copykat.pred"])
    else:
        registry[["copykat_call", "copykat_raw_call"]] = pd.NA
    registry.loc[registry["infercnv_run_state"].eq("not_eligible"), "infercnv_call"] = "not_eligible"
    registry.loc[registry["copykat_run_state"].eq("not_eligible"), "copykat_call"] = "not_eligible"
    completed_infer_missing = registry["infercnv_run_state"].eq("completed") & registry["infercnv_call"].isna()
    if completed_infer_missing.any():
        raise RuntimeError(f"Completed inferCNV cells without calls: {int(completed_infer_missing.sum())}")
    # Match the rebuilt Part 1 summary semantics: cells absent from a completed
    # CopyKAT prediction table are method-unavailable at cell level.
    completed_copy_missing = registry["copykat_run_state"].eq("completed") & registry["copykat_call"].isna()
    registry.loc[completed_copy_missing, "copykat_call"] = "not_eligible"
    registry["copykat_cell_filtered_or_uncalled"] = completed_copy_missing
    infer_pending = registry["infercnv_run_state"].isin({"pending", "running", "failed"})
    copy_pending = registry["copykat_run_state"].isin({"pending", "running", "failed"})
    integrated = [
        truth_table(str(i), str(c), bool(ip), bool(cp))
        for i, c, ip, cp in zip(
            registry["infercnv_call"].fillna("pending"),
            registry["copykat_call"].fillna("pending"),
            infer_pending,
            copy_pending,
        )
    ]
    registry["external_cnv_status_incremental_v1"] = [row[0] for row in integrated]
    registry["external_malignant_conservative_incremental_v1"] = [row[1] for row in integrated]
    registry["external_cnv_evidence_incremental_v1"] = [row[2] for row in integrated]
    registry["external_cnv_cell_finalizable_incremental_v1"] = [row[3] for row in integrated]

    gate_columns = [
        "cell_id", "dataset", "object_name", "sample_id", "patient_id", "patient_key",
        "phase7_broad_annotation", "phase7_fine_annotation", "phase7_annotation_confidence",
        "phase8_spatial_reference_state", "phase8_spatial_mapping_status",
        "phase8_spatial_annotation_confidence", "phase8_f1_projection_score",
        "phase8_f1_percentile_patient", "phase8_f1_eligible_patient",
        "part1__infercnv_call_v1", "part1__copykat_call_v1", "part1__malignancy_status_v1",
        "part1__malignant_conservative_v1", "part1__malignancy_evidence_v1",
        "phase9_part1_identical_cell", "phase9_part1_cnv_assessed", "phase9_cnv_gate_status",
        "phase9_malignancy_supported", "phase9_non_malignant_supported",
    ]
    gate = pd.read_csv(phase9_gate, sep="\t", compression="gzip", dtype=str, usecols=gate_columns)
    if not gate["cell_id"].is_unique:
        raise RuntimeError("Phase 9 gate cell IDs are not unique")
    gate = gate.loc[gate["cell_id"].isin(set(registry["exact_cell_id"]))].copy()
    registry = registry.merge(gate, left_on="exact_cell_id", right_on="cell_id", how="left", validate="one_to_one")
    if registry["cell_id"].isna().any() or not registry["phase7_broad_annotation"].eq("Epithelial").all():
        raise RuntimeError("External eligible registry is not fully current-epithelial/Phase9-mapped")

    status_counts = registry.groupby(
        ["external_cnv_status_incremental_v1", "external_cnv_cell_finalizable_incremental_v1"],
        observed=True,
    ).size().rename("n_cells").reset_index()
    unit_counts = validation_df.groupby(
        ["unit_type", "run_state", "validation_status"], observed=True
    ).size().rename("n_units").reset_index()
    control = registry.loc[
        registry["phase9_part1_cnv_assessed"].str.lower().eq("true")
        & registry["external_cnv_cell_finalizable_incremental_v1"]
    ].copy()
    control_summary = control.groupby(
        ["phase9_cnv_gate_status", "external_cnv_status_incremental_v1"], observed=True
    ).size().rename("n_cells").reset_index()

    atomic_frame(validation_df, out / "unit_validation_incremental_v1.tsv")
    atomic_frame(unit_counts, out / "unit_status_counts_incremental_v1.tsv")
    atomic_frame(infer_calls_df, out / "infercnv_query_calls_incremental_v1.tsv.gz")
    atomic_frame(infer_burden_df, out / "infercnv_group_burden_incremental_v1.tsv.gz")
    atomic_frame(copy_calls_df, out / "copykat_query_calls_incremental_v1.tsv.gz")
    atomic_frame(registry, out / "external_cnv_cell_consensus_incremental_v1.tsv.gz")
    atomic_frame(status_counts, out / "external_cnv_status_counts_incremental_v1.tsv")
    atomic_frame(control_summary, out / "part1_control_cross_table_incremental_v1.tsv")

    failures = validation_df.loc[validation_df["validation_status"].eq("FAIL")]
    metadata = {
        "schema_version": "external_cnv_incremental_snapshot_v1",
        "status": "PASS" if failures.empty else "PASS_WITH_FAILED_UNITS",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "phase10_root": str(phase10),
        "phase9_gate": str(phase9_gate),
        "phase9_gate_sha256": sha256(phase9_gate),
        "script": str(Path(__file__).resolve()),
        "script_sha256": sha256(Path(__file__).resolve()),
        "n_infercnv_units_completed": int(((validation_df["unit_type"] == "infercnv") & (validation_df["run_state"] == "completed")).sum()),
        "n_copykat_units_completed": int(((validation_df["unit_type"] == "copykat") & (validation_df["run_state"] == "completed")).sum()),
        "n_copykat_units_not_estimable": int(((validation_df["unit_type"] == "copykat") & (validation_df["run_state"] == "not_estimable")).sum()),
        "n_failed_units": int(len(failures)),
        "n_external_eligible_query_cells_union": int(len(registry)),
        "n_incrementally_finalizable_cells": int(registry["external_cnv_cell_finalizable_incremental_v1"].sum()),
        "n_part1_assessed_finalizable_controls": int(len(control)),
        "infercnv_rule": "per-unit query group burden > max(reference median + 3*scaled MAD, reference q95)",
        "malignant_conservative_rule": "Part1 rebuild truth table; an eligible but non-terminal second method forces pending_other_method",
        "interpretation": "incremental technical snapshot only; final cell labels wait for all eligible units and final independent validation",
    }
    atomic_json(out / "RUN_METADATA.json", metadata)
    manifest_rows = []
    for path in sorted(p for p in out.iterdir() if p.is_file() and p.name not in {"FILE_MANIFEST.tsv", "RUN_COMPLETE"}):
        manifest_rows.append({"file": path.name, "size_bytes": path.stat().st_size, "sha256": sha256(path)})
    atomic_frame(pd.DataFrame(manifest_rows), out / "FILE_MANIFEST.tsv")
    run_part = out / "RUN_COMPLETE.part"
    run_part.write_text(datetime.now(timezone.utc).isoformat() + "\n", encoding="utf-8")
    os.replace(run_part, out / "RUN_COMPLETE")
    print(json.dumps(metadata, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
