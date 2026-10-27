#!/usr/bin/env python3
"""Build no-expression external CNV unit manifests from validated Phase 8/9 metadata."""

from __future__ import annotations

import argparse
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


LYMPHOID = {"T_cell", "NK_cell", "B_cell"}


def atomic_tsv(frame: pd.DataFrame, final: Path) -> None:
    if final.exists():
        raise FileExistsError(f"No-overwrite gate: {final}")
    part = final.with_name(final.name + ".part")
    frame.to_csv(part, sep="\t", index=False, na_rep="NA")
    os.replace(part, final)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--crosswalk", required=True)
    parser.add_argument("--part1-map", required=True)
    parser.add_argument("--out-dir", required=True)
    args = parser.parse_args()
    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    patient_path = out_dir / "external_infercnv_patient_preflight_v1.tsv"
    sample_path = out_dir / "external_copykat_sample_preflight_v1.tsv"
    policy_path = out_dir / "external_cnv_reference_policy_v1.tsv"
    status_path = out_dir / "EXTERNAL_CNV_PREFLIGHT_STATUS"
    for path in (patient_path, sample_path, policy_path, status_path):
        if path.exists():
            raise FileExistsError(f"No-overwrite gate: {path}")

    cols = [
        "cell_id", "dataset", "patient_id", "sample_id", "object_name",
        "tissue_status", "tissue_confirmed", "phase7_broad_annotation",
        "phase7_fine_annotation", "phase8_spatial_reference_state",
    ]
    cells = pd.read_csv(args.crosswalk, sep="\t", compression="infer", dtype=str, usecols=cols)
    direct = pd.read_csv(
        args.part1_map,
        sep="\t",
        compression="infer",
        dtype=str,
        usecols=["cell_id", "part1__cnv_assessed_v1"],
    )
    if len(cells) != 585119 or cells["cell_id"].duplicated().any():
        raise RuntimeError("Crosswalk conservation gate failed")
    cells = cells.merge(direct, on="cell_id", how="left", validate="one_to_one")
    cells["part1_cnv_assessed"] = cells["part1__cnv_assessed_v1"].str.lower().eq("true")

    # Dataset-specific tumor/normal interpretation frozen from confirmed sample metadata.
    g178_normal = cells["dataset"].eq("GSE178481") & (
        cells["sample_id"].str.contains("Normal", case=False, na=False)
        | cells["object_name"].str.contains("Normal", case=False, na=False)
    )
    cells["cnv_tissue_class"] = "Tumor"
    cells.loc[g178_normal, "cnv_tissue_class"] = "Normal"
    cells["current_epithelial"] = cells["phase7_broad_annotation"].eq("Epithelial")
    cells["current_lymphoid"] = cells["phase7_broad_annotation"].isin(LYMPHOID)
    cells["query_candidate"] = cells["current_epithelial"] & cells["cnv_tissue_class"].eq("Tumor")

    policy_rows = []
    patient_rows = []
    for (dataset, patient), group in cells.groupby(["dataset", "patient_id"], sort=True, observed=True):
        query = group["query_candidate"]
        matched_normal = group["current_epithelial"] & group["cnv_tissue_class"].eq("Normal")
        same_study_normal = (
            cells["dataset"].eq(dataset)
            & cells["current_epithelial"]
            & cells["cnv_tissue_class"].eq("Normal")
        )
        tumor_lymphoid = group["current_lymphoid"] & group["cnv_tissue_class"].eq("Tumor")
        if int(matched_normal.sum()) >= 50:
            mode, n_ref = "matched_patient_normal_epithelial", int(matched_normal.sum())
        elif int(same_study_normal.sum()) >= 50:
            mode, n_ref = "same_study_normal_epithelial_pool", int(same_study_normal.sum())
        elif int(tumor_lymphoid.sum()) >= 50:
            mode, n_ref = "matched_patient_tumor_lymphoid", int(tumor_lymphoid.sum())
        else:
            mode, n_ref = "no_valid_reference", 0
        n_query = int(query.sum())
        n_assessed = int((query & group["part1_cnv_assessed"]).sum())
        n_uncovered = n_query - n_assessed
        eligible = n_query >= 50 and n_uncovered > 0 and n_ref >= 50
        patient_rows.append({
            "unit_id": f"{dataset}__{patient}",
            "dataset": dataset,
            "patient_id": patient,
            "n_tumor_epithelial": n_query,
            "n_part1_cnv_assessed_query": n_assessed,
            "n_uncovered_query": n_uncovered,
            "n_matched_normal_epithelial": int(matched_normal.sum()),
            "n_same_study_normal_epithelial": int(same_study_normal.sum()),
            "n_matched_tumor_lymphoid": int(tumor_lymphoid.sum()),
            "reference_mode": mode,
            "n_reference_candidates": n_ref,
            "eligible_infercnv": eligible,
            "status": "queued" if eligible else "ineligible",
        })
        policy_rows.append({"dataset": dataset, "patient_id": patient, "reference_mode": mode, "n_reference_candidates": n_ref})

    patient_manifest = pd.DataFrame(patient_rows)
    sample_rows = []
    for (dataset, patient, sample), group in cells.groupby(["dataset", "patient_id", "sample_id"], sort=True, observed=True):
        query = group["query_candidate"]
        if int(query.sum()) == 0:
            continue
        patient_policy = patient_manifest.loc[
            patient_manifest["unit_id"].eq(f"{dataset}__{patient}")
        ].iloc[0]
        n_query = int(query.sum())
        n_assessed = int((query & group["part1_cnv_assessed"]).sum())
        n_uncovered = n_query - n_assessed
        n_ref = int(patient_policy["n_reference_candidates"])
        eligible = n_query >= 50 and n_uncovered > 0 and n_ref >= 50
        sample_rows.append({
            "unit_id": f"{dataset}__{sample}",
            "dataset": dataset,
            "patient_id": patient,
            "sample_id": sample,
            "n_tumor_epithelial": n_query,
            "n_part1_cnv_assessed_query": n_assessed,
            "n_uncovered_query": n_uncovered,
            "reference_mode": patient_policy["reference_mode"],
            "n_reference_candidates": n_ref,
            "eligible_copykat": eligible,
            "status": "queued" if eligible else "ineligible",
        })
    sample_manifest = pd.DataFrame(sample_rows)
    policy = pd.DataFrame(policy_rows)

    atomic_tsv(patient_manifest, patient_path)
    atomic_tsv(sample_manifest, sample_path)
    atomic_tsv(policy, policy_path)
    counts = {
        "cells": len(cells),
        "datasets": int(cells["dataset"].nunique()),
        "patients": int(cells[["dataset", "patient_id"]].drop_duplicates().shape[0]),
        "infercnv_units_total": len(patient_manifest),
        "infercnv_units_eligible": int(patient_manifest["eligible_infercnv"].sum()),
        "copykat_units_total": len(sample_manifest),
        "copykat_units_eligible": int(sample_manifest["eligible_copykat"].sum()),
        "eligible_infercnv_uncovered_queries": int(patient_manifest.loc[patient_manifest["eligible_infercnv"], "n_uncovered_query"].sum()),
    }
    status = {
        "status": "PASS",
        "completed_at_utc": datetime.now(timezone.utc).isoformat(),
        "counts": counts,
        "policy": "Five datasets remain separate. Query=current epithelial cells in tumor-class samples. Reference hierarchy mirrors Part 1 rebuild: matched normal epithelial, same-study normal epithelial, then matched tumor lymphoid. Minimum query/reference=50. Existing Part 1 calls are retained and act as internal controls.",
    }
    atomic_text = json.dumps(status, indent=2, sort_keys=True) + "\n"
    part = status_path.with_name(status_path.name + ".part")
    part.write_text(atomic_text, encoding="utf-8")
    os.replace(part, status_path)
    print(json.dumps(counts, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
