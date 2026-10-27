#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from datetime import datetime
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import sparse


ADDED = {"N3_1", "N4_1"}


def decode_obs(handle: h5py.File, field: str) -> np.ndarray:
    values = handle[f"/obs/{field}"][()]
    category_path = f"/obs/__categories/{field}"
    if category_path in handle:
        categories = handle[category_path][()]
        categories = np.asarray([item.decode() if isinstance(item, bytes) else str(item) for item in categories], dtype=object)
        codes = np.asarray(values, dtype=int)
        answer = np.full(len(codes), "NA", dtype=object)
        valid = (codes >= 0) & (codes < len(categories))
        answer[valid] = categories[codes[valid]]
        return answer.astype(str)
    return np.asarray([item.decode() if isinstance(item, bytes) else str(item) for item in values], dtype=str)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_tsv(frame: pd.DataFrame, path: Path) -> None:
    part = path.with_name(path.name + ".part")
    frame.to_csv(part, sep="\t", index=False, na_rep="NA", encoding="utf-8", lineterminator="\n")
    os.replace(part, path)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--rebuild-root", type=Path, required=True)
    parser.add_argument("--signature", type=Path, required=True)
    args = parser.parse_args()
    h5ad = args.rebuild_root / "02_complete_cohort_build" / "gse242299_h5ad_recovery_v3" / "GSE242299_all_cells_50236_33538.h5ad"
    final = args.rebuild_root / "03_incremental_sensitivity" / "gse242299_added_normals_v1"
    part = final.with_name(final.name + ".part")
    if final.exists() or part.exists():
        raise FileExistsError("Refusing existing GSE242299 incremental output")
    part.mkdir(parents=True)

    signatures = pd.read_csv(args.signature, sep="\t")
    f1 = signatures.loc[signatures["signature_id"] == "f1_top100_weighted", ["gene", "raw_weight"]].copy()
    if len(f1) != 100 or f1["gene"].duplicated().any():
        raise RuntimeError("Frozen F1 signature is not exactly 100 unique genes")

    with h5py.File(h5ad, "r") as handle:
        shape = tuple(int(x) for x in handle["/X"].attrs["shape"])
        matrix = sparse.csr_matrix(
            (handle["/X/data"][()], handle["/X/indices"][()], handle["/X/indptr"][()]),
            shape=shape,
        )
        genes_raw = handle["/var/_index"][()]
        genes = np.asarray([item.decode() if isinstance(item, bytes) else str(item) for item in genes_raw])
        libraries = decode_obs(handle, "library_name")
        patient = decode_obs(handle, "patient")
        tissue = decode_obs(handle, "tissue")
        broad = decode_obs(handle, "broad_cell_type")
        cell_type = decode_obs(handle, "cell_type")

    if matrix.shape != (50236, 33538) or len(libraries) != matrix.shape[0]:
        raise RuntimeError("H5AD dimension or obs alignment mismatch")
    group = np.where(np.isin(libraries, sorted(ADDED)), "added2", "old29")
    if int(np.sum(group == "added2")) != 1663 or set(libraries[group == "added2"]) != ADDED:
        raise RuntimeError("Added-library cell conservation failed")

    gene_lookup = {gene: index for index, gene in enumerate(genes)}
    if "TYMP" not in gene_lookup:
        raise RuntimeError("TYMP absent from official H5AD feature index")
    available = f1.loc[f1["gene"].isin(gene_lookup)].copy()
    if len(available) < 80:
        raise RuntimeError(f"Insufficient F1 signature coverage: {len(available)}/100")
    columns = np.asarray([gene_lookup[gene] for gene in available["gene"]], dtype=int)
    weights = available["raw_weight"].to_numpy(float)
    weights = weights / weights.sum()

    totals = np.asarray(matrix.sum(axis=1)).ravel()
    if np.any(totals <= 0):
        raise RuntimeError("One or more H5AD cells have zero total counts")
    selected = matrix[:, columns].toarray().astype(float)
    lognorm = np.log1p(selected / totals[:, None] * 10000.0)
    means = lognorm.mean(axis=0)
    standard_deviation = lognorm.std(axis=0, ddof=0)
    standard_deviation[standard_deviation == 0] = 1.0
    f1_primary = ((lognorm - means) / standard_deviation) @ weights
    f1_sensitivity = lognorm @ weights
    tymp_counts = np.asarray(matrix[:, gene_lookup["TYMP"]].todense()).ravel()
    tymp_log = np.log1p(tymp_counts / totals * 10000.0)

    cell = pd.DataFrame({
        "comparison_group": group,
        "library_id": libraries,
        "patient": patient,
        "tissue": tissue,
        "broad_cell_type": broad,
        "cell_type": cell_type,
        "total_counts": totals,
        "tymp_detected": tymp_counts > 0,
        "tymp_log1p_cp10k": tymp_log,
        "f1_primary_gene_z_weighted": f1_primary,
        "f1_sensitivity_log1p_cp10k_weighted": f1_sensitivity,
    })

    composition = (
        cell.groupby(["comparison_group", "library_id", "tissue", "broad_cell_type", "cell_type"], dropna=False)
        .size().rename("cells").reset_index()
    )
    composition["fraction_within_library"] = composition["cells"] / composition.groupby(["comparison_group", "library_id"])["cells"].transform("sum")
    group_composition = cell.groupby(["comparison_group", "broad_cell_type", "cell_type"], dropna=False).size().rename("cells").reset_index()
    group_composition["fraction_within_group"] = group_composition["cells"] / group_composition.groupby("comparison_group")["cells"].transform("sum")
    score_summary = (
        cell.groupby(["comparison_group", "broad_cell_type"], dropna=False)
        .agg(
            cells=("library_id", "size"),
            libraries=("library_id", "nunique"),
            tymp_detected_cells=("tymp_detected", "sum"),
            tymp_detection_fraction=("tymp_detected", "mean"),
            tymp_mean_log1p_cp10k=("tymp_log1p_cp10k", "mean"),
            f1_primary_mean=("f1_primary_gene_z_weighted", "mean"),
            f1_sensitivity_mean=("f1_sensitivity_log1p_cp10k_weighted", "mean"),
        ).reset_index()
    )
    library_summary = (
        cell.groupby(["comparison_group", "library_id", "patient", "tissue"], dropna=False)
        .agg(
            cells=("library_id", "size"),
            broad_cell_types=("broad_cell_type", "nunique"),
            tymp_detection_fraction=("tymp_detected", "mean"),
            tymp_mean_log1p_cp10k=("tymp_log1p_cp10k", "mean"),
            f1_primary_mean=("f1_primary_gene_z_weighted", "mean"),
            f1_sensitivity_mean=("f1_sensitivity_log1p_cp10k_weighted", "mean"),
        ).reset_index()
    )

    atomic_tsv(composition, part / "added_library_celltype_composition_v1.tsv")
    atomic_tsv(group_composition, part / "old29_vs_added2_composition_v1.tsv")
    atomic_tsv(score_summary, part / "old29_vs_added2_tymp_f1_summary_v1.tsv")
    atomic_tsv(library_summary, part / "library_tymp_f1_summary_v1.tsv")
    metadata = {
        "status": "PASS",
        "generated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "h5ad": str(h5ad),
        "h5ad_sha256": sha256(h5ad),
        "signature": str(args.signature),
        "signature_sha256": sha256(args.signature),
        "old_cells": int(np.sum(group == "old29")),
        "added_cells": int(np.sum(group == "added2")),
        "added_libraries": sorted(ADDED),
        "f1_genes_available": int(len(available)),
        "normalization": "log1p(count/total_count*10000) used only for targeted TYMP/F1 score comparison",
        "integration_performed": False,
        "clustering_performed": False,
        "annotation_source": "official H5AD obs",
    }
    (part / "run_metadata_v1.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    (part / "GSE242299_INCREMENTAL_PASS_v1").write_text(
        "status=PASS\nold_libraries=29\nadded_libraries=N3_1|N4_1\nold_cells=48573\nadded_cells=1663\n"
        f"f1_genes_available={len(available)}\nintegration_performed=0\nclustering_performed=0\n",
        encoding="utf-8",
    )
    os.replace(part, final)
    print(f"PASS added_cells=1663 f1_genes={len(available)}")


if __name__ == "__main__":
    main()
