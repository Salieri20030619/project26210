#!/usr/bin/env python3
import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path

import anndata as ad
import pandas as pd

PANEL = [
    "TYMP", "MTAP", "ADORA2B", "STAT3", "SPP1", "ARG1", "VEGFA", "IL10", "TIMP1",
    "GPNMB", "LGALS3", "FABP5", "CD9", "LPL", "C1QA", "C1QB", "C1QC",
    "FCN1", "VCAN", "S100A8", "S100A9", "CTSS", "LILRB1", "TGFB1", "VSIG4",
    "MARCO", "CD274", "IDO1", "MAT2A", "MAT2B", "AMD1", "SRM", "SMS",
]


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def audit_h5ad(cohort, h5ad_path, metadata_path):
    obj = ad.read_h5ad(h5ad_path, backed="r")
    meta = pd.read_csv(metadata_path, sep="\t", compression="gzip", dtype=str)
    obs_ids = pd.Index(obj.obs_names.astype(str))
    cell_ids = pd.Index(meta["cell_id"].astype(str))
    overlap = obs_ids.intersection(cell_ids)
    genes = set(map(str, obj.var_names))
    myeloid = meta.loc[meta["lineage_v1"].eq("Myeloid_DC")].copy()
    units = myeloid[["patient_key", "sample_id"]].drop_duplicates()
    return {
        "cohort": cohort,
        "h5ad": str(h5ad_path),
        "metadata": str(metadata_path),
        "n_cells_h5ad": int(obj.n_obs),
        "n_genes_h5ad": int(obj.n_vars),
        "n_metadata_rows": int(len(meta)),
        "n_exact_id_overlap": int(len(overlap)),
        "exact_id_match_fraction_h5ad": float(len(overlap) / max(1, obj.n_obs)),
        "n_myeloid_dc_cells": int(len(myeloid)),
        "n_myeloid_dc_patients": int(myeloid["patient_key"].nunique()),
        "n_myeloid_dc_patient_samples": int(len(units)),
        "panel_present": [g for g in PANEL if g in genes],
        "panel_missing": [g for g in PANEL if g not in genes],
        "obs_columns": list(map(str, obj.obs.columns)),
        "metadata_columns": list(map(str, meta.columns)),
        "h5ad_sha256": sha256(h5ad_path),
        "metadata_sha256": sha256(metadata_path),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    root60 = Path("/mnt/BioDisk/wangmingyuan/project26210_immunotherapy_scrna_phase60_v1/01_qc/cohorts")
    root67 = Path("/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1_phase67_direct_program_replication_v1/inputs_v1/metadata_v2")
    rows = []
    for cohort in ("GSE264586", "GSE294109"):
        rows.append(audit_h5ad(
            cohort,
            root60 / f"{cohort}.qc_filtered_raw_v1.h5ad",
            root67 / f"{cohort}_phase67_metadata_v1.tsv.gz",
        ))
    gse207_meta = Path("/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1_phase64_large_immune_validation_v1/results_v1/GSE207493/gse207493_cell_metadata_scores_v1.tsv.gz")
    g207 = pd.read_csv(gse207_meta, sep="\t", compression="gzip", dtype=str)
    rds_root = Path("/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase4a_10x_seurat_v1/objects")
    rds_files = sorted(rds_root.glob("GSE207493__*.rds"))
    g207_summary = {
        "cohort": "GSE207493",
        "metadata": str(gse207_meta),
        "metadata_sha256": sha256(gse207_meta),
        "n_metadata_rows": int(len(g207)),
        "metadata_columns": list(map(str, g207.columns)),
        "n_rds_files": len(rds_files),
        "rds_files": [{"path": str(p), "size": p.stat().st_size} for p in rds_files],
    }
    payload = {"external_h5ad": rows, "gse207493": g207_summary}
    (out / "phase71_input_audit_v1.json").write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    with (out / "phase71_input_audit_summary_v1.tsv").open("w", encoding="utf-8", newline="") as fh:
        fields = ["cohort", "n_cells_h5ad", "n_genes_h5ad", "n_metadata_rows", "n_exact_id_overlap", "exact_id_match_fraction_h5ad", "n_myeloid_dc_cells", "n_myeloid_dc_patients", "n_myeloid_dc_patient_samples", "panel_present", "panel_missing", "h5ad_sha256", "metadata_sha256"]
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for row in rows:
            cooked = dict(row)
            cooked["panel_present"] = ";".join(row["panel_present"])
            cooked["panel_missing"] = ";".join(row["panel_missing"])
            writer.writerow({k: cooked.get(k, "") for k in fields})
    print(json.dumps({"status": "PASS", "cohorts": len(rows), "gse207493_rds": len(rds_files)}))


if __name__ == "__main__":
    main()
