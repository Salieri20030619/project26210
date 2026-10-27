#!/usr/bin/env python3
"""Audit the rebuilt Part 1 input before inferCNV/CopyKAT analysis."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import subprocess
import sys
from datetime import datetime
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd


PROJECT = Path("/mnt/BioDisk/ludingyu/project26210")
DEFAULT_INPUT = PROJECT / (
    "results/part1_rebuild_v1/02_integration_annotation/cell_annotation_v1/"
    "part1_scvi_integrated_cell_annotated_v1.h5ad"
)
DEFAULT_GENE_ORDER = PROJECT / (
    "results/infercnv_rerun_seurat_v5/"
    "gene_order_EnsDb.Hsapiens.v86_seurat_genes.txt"
)
DEFAULT_OUTPUT = Path(
    "/mnt/BioDisk/wangmingyuan/project26210_part1_cnv_stage/"
    "malignancy_cnv_input_audit_v1"
)
DEFAULT_RSCRIPT = Path(
    "/mnt/BioDisk/ludingyu/miniconda3/envs/r453/bin/Rscript"
)

MIN_TUMOR_EPITHELIAL = 50
MIN_REFERENCE = 50
CLEAR_EPITHELIAL = "Epithelial cell"
MIXED_EPITHELIAL = "Epithelial/Myeloid mixed"
LYMPHOID_REFERENCE_LABELS = ("T cell", "NK cell", "B cell")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-h5ad", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--gene-order", type=Path, default=DEFAULT_GENE_ORDER)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--rscript", type=Path, default=DEFAULT_RSCRIPT)
    return parser.parse_args()


def sha256_file(path: Path, block_size: int = 16 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            block = handle.read(block_size)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def write_tsv(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(path, sep="\t", index=False, lineterminator="\n")


def r_capabilities(rscript: Path) -> dict[str, str | bool]:
    expression = r'''
pkgs <- c("Seurat", "Matrix", "data.table", "infercnv", "copykat")
cat("r_version\t", R.version.string, "\n", sep="")
cat("r_home\t", R.home(), "\n", sep="")
cat("rscript_resolved\t", Sys.which("Rscript"), "\n", sep="")
for (p in pkgs) {
  ok <- requireNamespace(p, quietly=TRUE)
  cat("package_", p, "_available\t", ok, "\n", sep="")
  cat("package_", p, "_version\t", if (ok) as.character(packageVersion(p)) else "NA", "\n", sep="")
}
'''
    result = subprocess.run(
        [str(rscript), "--vanilla", "-e", expression],
        check=True,
        text=True,
        capture_output=True,
    )
    values: dict[str, str | bool] = {}
    for line in result.stdout.splitlines():
        if "\t" not in line:
            continue
        key, value = line.split("\t", 1)
        value = value.strip()
        if value in {"TRUE", "FALSE"}:
            values[key] = value == "TRUE"
        else:
            values[key] = value
    values["rscript_requested"] = str(rscript)
    values["r_stderr"] = result.stderr.strip()
    return values


def reference_mode(
    matched_normal: int,
    study_normal: int,
    matched_lymphoid: int,
    all_normal: int,
) -> str:
    if matched_normal >= MIN_REFERENCE:
        return "matched_patient_normal_epithelial"
    if study_normal >= MIN_REFERENCE:
        return "same_study_normal_epithelial_pool"
    if matched_lymphoid >= MIN_REFERENCE:
        return "matched_patient_tumor_lymphoid"
    if all_normal >= MIN_REFERENCE:
        return "balanced_cross_study_normal_epithelial_pool"
    return "unavailable"


def main() -> int:
    args = parse_args()
    started = datetime.now().astimezone()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    for path in (args.input_h5ad, args.gene_order, args.rscript):
        if not path.exists():
            raise FileNotFoundError(path)

    print(f"[{datetime.now().astimezone().isoformat()}] read obs from {args.input_h5ad}")
    adata = ad.read_h5ad(args.input_h5ad, backed="r")
    obs = adata.obs.copy()
    var_names = pd.Index(adata.var_names.astype(str))
    cell_ids = pd.Index(adata.obs_names.astype(str))

    required = {
        "sample_id",
        "study",
        "curated_patient_id",
        "source_tissue",
        "cell_type_broad_v1",
        "annotation_confidence_v1",
        "annotation_ambiguous_v1",
    }
    missing = sorted(required.difference(obs.columns))
    if missing:
        raise RuntimeError(f"Missing required obs columns: {missing}")

    for col in required:
        obs[col] = obs[col].astype(str)
    obs.insert(0, "cell_id", cell_ids)

    gene_order = pd.read_csv(
        args.gene_order,
        sep="\t",
        header=None,
        names=["gene", "chromosome", "start", "end"],
        dtype={"gene": str, "chromosome": str, "start": np.int64, "end": np.int64},
    )
    gene_order["in_current_h5ad"] = gene_order["gene"].isin(var_names)
    gene_order["is_autosome"] = gene_order["chromosome"].isin(
        [f"chr{i}" for i in range(1, 23)] + [str(i) for i in range(1, 23)]
    )

    tumor_epi = obs[
        (obs["source_tissue"] == "Tumor")
        & (obs["cell_type_broad_v1"] == CLEAR_EPITHELIAL)
    ].copy()
    normal_epi = obs[
        (obs["source_tissue"] == "Normal")
        & (obs["cell_type_broad_v1"] == CLEAR_EPITHELIAL)
    ].copy()
    tumor_lymphoid = obs[
        (obs["source_tissue"] == "Tumor")
        & obs["cell_type_broad_v1"].isin(LYMPHOID_REFERENCE_LABELS)
    ].copy()
    mixed = obs[obs["cell_type_broad_v1"] == MIXED_EPITHELIAL].copy()

    matched_normal = normal_epi.groupby("curated_patient_id").size().to_dict()
    study_normal = normal_epi.groupby("study").size().to_dict()
    matched_lymphoid = tumor_lymphoid.groupby("curated_patient_id").size().to_dict()
    all_normal = len(normal_epi)

    patient_rows: list[dict[str, object]] = []
    for patient_id, group in tumor_epi.groupby("curated_patient_id", sort=True):
        studies = sorted(group["study"].unique())
        if len(studies) != 1:
            raise RuntimeError(f"Patient spans multiple studies: {patient_id}: {studies}")
        study = studies[0]
        n_matched = int(matched_normal.get(patient_id, 0))
        n_study = int(study_normal.get(study, 0))
        n_lymphoid = int(matched_lymphoid.get(patient_id, 0))
        mode = reference_mode(n_matched, n_study, n_lymphoid, all_normal)
        n_query = len(group)
        sample_ids = sorted(group["sample_id"].unique())
        patient_rows.append(
            {
                "curated_patient_id": patient_id,
                "study": study,
                "n_tumor_epithelial": n_query,
                "n_tumor_sample_ids": len(sample_ids),
                "tumor_sample_ids": ";".join(sample_ids),
                "n_matched_normal_epithelial": n_matched,
                "n_same_study_normal_epithelial": n_study,
                "n_matched_tumor_lymphoid": n_lymphoid,
                "reference_mode": mode,
                "eligible_infercnv": n_query >= MIN_TUMOR_EPITHELIAL and mode != "unavailable",
                "status": "queued" if n_query >= MIN_TUMOR_EPITHELIAL and mode != "unavailable" else "ineligible",
            }
        )
    patient_manifest = pd.DataFrame(patient_rows).sort_values(
        ["study", "curated_patient_id"], ignore_index=True
    )

    sample_rows: list[dict[str, object]] = []
    for sample_id, group in tumor_epi.groupby("sample_id", sort=True):
        patients = sorted(group["curated_patient_id"].unique())
        studies = sorted(group["study"].unique())
        if len(patients) != 1 or len(studies) != 1:
            raise RuntimeError(
                f"Tumor sample is not unique to one patient/study: {sample_id}: "
                f"patients={patients}, studies={studies}"
            )
        patient_id = patients[0]
        study = studies[0]
        n_matched = int(matched_normal.get(patient_id, 0))
        n_study = int(study_normal.get(study, 0))
        n_lymphoid = int(matched_lymphoid.get(patient_id, 0))
        mode = reference_mode(n_matched, n_study, n_lymphoid, all_normal)
        n_query = len(group)
        sample_rows.append(
            {
                "sample_id": sample_id,
                "curated_patient_id": patient_id,
                "study": study,
                "n_tumor_epithelial": n_query,
                "n_matched_normal_epithelial": n_matched,
                "n_same_study_normal_epithelial": n_study,
                "n_matched_tumor_lymphoid": n_lymphoid,
                "reference_mode": mode,
                "eligible_copykat": n_query >= MIN_TUMOR_EPITHELIAL and mode != "unavailable",
                "status": "queued" if n_query >= MIN_TUMOR_EPITHELIAL and mode != "unavailable" else "ineligible",
            }
        )
    sample_manifest = pd.DataFrame(sample_rows).sort_values(
        ["study", "curated_patient_id", "sample_id"], ignore_index=True
    )

    celltype_summary = (
        obs.groupby(["source_tissue", "cell_type_broad_v1"], observed=True)
        .size()
        .rename("n_cells")
        .reset_index()
        .sort_values(["source_tissue", "n_cells"], ascending=[True, False])
    )
    epithelial_by_study = (
        obs[obs["cell_type_broad_v1"].isin([CLEAR_EPITHELIAL, MIXED_EPITHELIAL])]
        .groupby(["study", "source_tissue", "cell_type_broad_v1"], observed=True)
        .size()
        .rename("n_cells")
        .reset_index()
        .sort_values(["study", "source_tissue", "cell_type_broad_v1"])
    )
    mixed_by_patient = (
        mixed.groupby(["study", "curated_patient_id", "source_tissue"], observed=True)
        .size()
        .rename("n_mixed_cells")
        .reset_index()
        .sort_values(["study", "curated_patient_id", "source_tissue"])
    )

    print(f"[{datetime.now().astimezone().isoformat()}] inspect HDF5 counts layer")
    with h5py.File(args.input_h5ad, "r") as handle:
        counts_group = handle["layers/counts"]
        counts_shape = tuple(int(x) for x in counts_group.attrs["shape"])
        counts_dtype = str(counts_group["data"].dtype)
        counts_nnz = int(counts_group["data"].shape[0])

    print(f"[{datetime.now().astimezone().isoformat()}] inspect R CNV environment")
    r_info = r_capabilities(args.rscript)

    checks: list[dict[str, object]] = []

    def add_check(name: str, passed: bool, observed: object, expected: object) -> None:
        checks.append(
            {
                "check": name,
                "pass": bool(passed),
                "observed": str(observed),
                "expected": str(expected),
            }
        )

    add_check("input_n_obs", adata.n_obs == 258120, adata.n_obs, 258120)
    add_check("input_n_vars", adata.n_vars == 37795, adata.n_vars, 37795)
    add_check("cell_ids_unique", not cell_ids.has_duplicates, cell_ids.has_duplicates, False)
    add_check("gene_ids_unique", not var_names.has_duplicates, var_names.has_duplicates, False)
    add_check("counts_layer_present", "counts" in adata.layers, list(adata.layers.keys()), "contains counts")
    add_check("counts_shape_matches", counts_shape == (adata.n_obs, adata.n_vars), counts_shape, (adata.n_obs, adata.n_vars))
    add_check("counts_dtype_integer", counts_dtype.startswith("int"), counts_dtype, "integer dtype")
    add_check("counts_nnz_positive", counts_nnz > 0, counts_nnz, ">0")
    add_check("source_tissue_levels", set(obs["source_tissue"]) == {"Tumor", "Normal"}, sorted(set(obs["source_tissue"])), ["Normal", "Tumor"])
    add_check("clear_epithelial_count", len(tumor_epi) + len(normal_epi) == 62652, len(tumor_epi) + len(normal_epi), 62652)
    add_check("mixed_epithelial_count", len(mixed) == 918, len(mixed), 918)
    add_check("normal_epithelial_reference_available", len(normal_epi) >= MIN_REFERENCE, len(normal_epi), f">={MIN_REFERENCE}")
    add_check("infercnv_eligible_patients", int(patient_manifest["eligible_infercnv"].sum()) >= 30, int(patient_manifest["eligible_infercnv"].sum()), ">=30")
    add_check("copykat_eligible_samples", int(sample_manifest["eligible_copykat"].sum()) >= 35, int(sample_manifest["eligible_copykat"].sum()), ">=35")
    add_check("gene_order_unique", not gene_order["gene"].duplicated().any(), int(gene_order["gene"].duplicated().sum()), 0)
    add_check("gene_order_overlap", int(gene_order["in_current_h5ad"].sum()) >= 20000, int(gene_order["in_current_h5ad"].sum()), ">=20000")
    add_check("autosomal_gene_order_overlap", int((gene_order["in_current_h5ad"] & gene_order["is_autosome"]).sum()) >= 19000, int((gene_order["in_current_h5ad"] & gene_order["is_autosome"]).sum()), ">=19000")
    for package_name in ("infercnv", "copykat"):
        key = f"package_{package_name}_available"
        add_check(key, r_info.get(key) is True, r_info.get(key), True)

    checks_frame = pd.DataFrame(checks)
    all_pass = bool(checks_frame["pass"].all())

    input_inventory = pd.DataFrame(
        [
            {
                "role": "annotated_h5ad",
                "path": str(args.input_h5ad),
                "size_bytes": args.input_h5ad.stat().st_size,
                "sha256": sha256_file(args.input_h5ad),
            },
            {
                "role": "gene_order",
                "path": str(args.gene_order),
                "size_bytes": args.gene_order.stat().st_size,
                "sha256": sha256_file(args.gene_order),
            },
        ]
    )
    software = pd.DataFrame(
        [{"component": key, "value": value} for key, value in r_info.items()]
        + [
            {"component": "python_executable", "value": sys.executable},
            {"component": "python_version", "value": sys.version.replace("\n", " ")},
            {"component": "platform", "value": platform.platform()},
            {"component": "anndata_version", "value": ad.__version__},
            {"component": "h5py_version", "value": h5py.__version__},
            {"component": "pandas_version", "value": pd.__version__},
            {"component": "numpy_version", "value": np.__version__},
        ]
    )

    write_tsv(input_inventory, args.output_dir / "cnv_input_inventory_v1.tsv")
    write_tsv(software, args.output_dir / "cnv_software_capabilities_v1.tsv")
    write_tsv(celltype_summary, args.output_dir / "celltype_by_tissue_counts_v1.tsv")
    write_tsv(epithelial_by_study, args.output_dir / "epithelial_candidates_by_study_v1.tsv")
    write_tsv(mixed_by_patient, args.output_dir / "mixed_epithelial_by_patient_v1.tsv")
    write_tsv(patient_manifest, args.output_dir / "infercnv_patient_manifest_preflight_v1.tsv")
    write_tsv(sample_manifest, args.output_dir / "copykat_sample_manifest_preflight_v1.tsv")
    write_tsv(
        gene_order,
        args.output_dir / "gene_order_compatibility_v1.tsv.gz",
    )
    write_tsv(checks_frame, args.output_dir / "malignancy_cnv_preflight_qc_v1.tsv")

    metadata = {
        "schema_version": "part1_malignancy_cnv_input_audit_v1",
        "started_at": started.isoformat(),
        "completed_at": datetime.now().astimezone().isoformat(),
        "status": "completed" if all_pass else "failed",
        "all_checks_pass": all_pass,
        "n_checks": len(checks_frame),
        "n_checks_passed": int(checks_frame["pass"].sum()),
        "input_h5ad": str(args.input_h5ad),
        "gene_order": str(args.gene_order),
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "counts_shape": counts_shape,
        "counts_dtype": counts_dtype,
        "counts_nnz": counts_nnz,
        "n_clear_epithelial": int(len(tumor_epi) + len(normal_epi)),
        "n_tumor_clear_epithelial": int(len(tumor_epi)),
        "n_normal_clear_epithelial": int(len(normal_epi)),
        "n_mixed_epithelial_excluded_primary": int(len(mixed)),
        "n_infercnv_patients_total": int(len(patient_manifest)),
        "n_infercnv_patients_eligible": int(patient_manifest["eligible_infercnv"].sum()),
        "n_copykat_samples_total": int(len(sample_manifest)),
        "n_copykat_samples_eligible": int(sample_manifest["eligible_copykat"].sum()),
        "min_tumor_epithelial": MIN_TUMOR_EPITHELIAL,
        "min_reference_cells": MIN_REFERENCE,
        "primary_candidate_definition": "source_tissue == Tumor and cell_type_broad_v1 == Epithelial cell",
        "mixed_candidate_policy": "excluded from primary analysis; retain for later sensitivity analysis",
        "reference_priority": [
            "matched_patient_normal_epithelial",
            "same_study_normal_epithelial_pool",
            "matched_patient_tumor_lymphoid_T_NK_B",
            "balanced_cross_study_normal_epithelial_pool",
        ],
        "cell_id_policy": "preserve exact H5AD IDs; later tool runs must use reversible safe surrogate IDs",
    }
    with (args.output_dir / "run_metadata_v1.json").open("w", encoding="utf-8") as handle:
        json.dump(metadata, handle, ensure_ascii=False, indent=2)
        handle.write("\n")

    marker = args.output_dir / ("RUN_COMPLETE" if all_pass else "RUN_FAILED")
    marker.write_text(datetime.now().astimezone().isoformat() + "\n", encoding="utf-8")
    print(
        "PART1_MALIGNANCY_CNV_INPUT_AUDIT_V1_"
        f"{'PASS' if all_pass else 'FAIL'} checks={len(checks_frame)} "
        f"patients={int(patient_manifest['eligible_infercnv'].sum())}/{len(patient_manifest)} "
        f"samples={int(sample_manifest['eligible_copykat'].sum())}/{len(sample_manifest)}"
    )
    return 0 if all_pass else 1


if __name__ == "__main__":
    raise SystemExit(main())
