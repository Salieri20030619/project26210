#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
from datetime import datetime
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.io import mmwrite


H5AD_SHA256 = "e4779c46af13a8643e8dec584eeff4322258181f5edab3e55dcd8c20bc28a9e7"
ASSIGNMENT_SHA256 = "f8d0d01dcb76b3a19055accae42be3a2e14c6ba3a9da85db30c9592f2e635e6d"
BRIDGE_SHA256 = "cde0288a09ceabdf3c5e4d9b7aefbba5d370c2c5082b03537edb24d847250123"
SEED = 20260815
CAP = 300
GROUPS = [
    "TYMP-high Malignant", "TYMP-low Malignant", "SPP1-state Macrophage",
    "C1QC-state Macrophage", "Other Macrophage", "Endothelial",
]
TARGETED_GENES = [
    "TYMP", "HIF1A", "MST1", "MST1R", "ICAM5", "ITGAL", "ITGAM", "ITGB2",
    "CD209", "PRTG", "IGDCC4", "CA9", "VEGFA", "ADM", "ANGPTL4", "LDHA",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def atomic_tsv(frame: pd.DataFrame, path: Path, gzip_output: bool = False) -> None:
    part = Path(str(path) + ".part")
    frame.to_csv(part, sep="\t", index=False, compression="gzip" if gzip_output else None,
                 encoding="utf-8", lineterminator="\n")
    os.replace(part, path)


def atomic_json(payload: dict, path: Path) -> None:
    part = Path(str(path) + ".part")
    part.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    os.replace(part, path)


def normalize_log1p_cp10k(counts: sp.spmatrix, library: np.ndarray) -> sp.csr_matrix:
    if np.any(~np.isfinite(library)) or np.any(library <= 0):
        raise RuntimeError("Invalid raw library sizes")
    result = sp.diags(10000.0 / library.astype(float)) @ counts.astype(np.float64)
    result = result.tocsr()
    result.data = np.log1p(result.data)
    return result


def validate_counts(counts: sp.spmatrix, label: str) -> None:
    data = counts.data
    if np.any(~np.isfinite(data)) or np.any(data < 0):
        raise RuntimeError(f"{label}: negative or non-finite counts")
    if np.any(np.abs(data - np.rint(data)) > 1e-6):
        raise RuntimeError(f"{label}: non-integer counts")


def deterministic_sample(ids: pd.Series, cap: int) -> list[str]:
    values = ids.astype(str).tolist()
    if len(values) <= cap:
        return sorted(values)
    ranked = sorted(values, key=lambda value: hashlib.sha256(f"{SEED}|{value}".encode()).hexdigest())
    return ranked[:cap]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", type=Path, required=True)
    parser.add_argument("--assignments", type=Path, required=True)
    parser.add_argument("--bridge", type=Path, required=True)
    parser.add_argument("--cellchat-genes", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    args.output.mkdir(parents=True)

    hashes = {"h5ad": sha256(args.h5ad), "assignments": sha256(args.assignments), "bridge": sha256(args.bridge)}
    expected = {"h5ad": H5AD_SHA256, "assignments": ASSIGNMENT_SHA256, "bridge": BRIDGE_SHA256}
    if hashes != expected:
        raise RuntimeError(f"Input hash mismatch: {hashes}")

    assignments = pd.read_csv(args.assignments, sep="\t", dtype={"exact_cell_id": str, "patient": str})
    bridge = pd.read_csv(args.bridge, sep="\t", dtype={"exact_cell_id": str, "patient": str})
    if assignments["exact_cell_id"].duplicated().any() or bridge["exact_cell_id"].duplicated().any():
        raise RuntimeError("Duplicate cell identifiers")
    if len(assignments) != 229032 or len(bridge) != 11921:
        raise RuntimeError("Frozen assignment or bridge scope changed")
    bridge["communication_group"] = np.where(bridge["TYMP_high"].astype(bool), "TYMP-high Malignant", "TYMP-low Malignant")
    bridge_group = bridge.set_index("exact_cell_id")["communication_group"]
    malignant_patients = sorted(bridge["patient"].unique())

    meta = assignments.loc[assignments["patient"].isin(malignant_patients)].copy()
    meta["communication_group"] = meta["state_v2"]
    malignant_mask = meta["exact_cell_id"].isin(bridge_group.index)
    meta.loc[malignant_mask, "communication_group"] = meta.loc[malignant_mask, "exact_cell_id"].map(bridge_group)
    meta = meta.loc[meta["communication_group"].isin(GROUPS)].copy()
    if int(meta["communication_group"].str.contains("Malignant").sum()) != 11921:
        raise RuntimeError("Malignant group conservation failure")

    selected_ids: list[str] = []
    sample_rows = []
    for (study, patient, group), frame in meta.groupby(["study_v1", "patient", "communication_group"], sort=True):
        chosen = deterministic_sample(frame["exact_cell_id"], CAP)
        selected_ids.extend(chosen)
        sample_rows.append({"study": study, "patient": patient, "communication_group": group,
                            "n_available": len(frame), "n_sampled": len(chosen), "cap": CAP})
    selected_set = set(selected_ids)
    sampled_meta = meta.loc[meta["exact_cell_id"].isin(selected_set), ["exact_cell_id", "patient", "study_v1", "communication_group"]].copy()
    sampled_meta = sampled_meta.set_index("exact_cell_id").loc[selected_ids].reset_index()

    obj = ad.read_h5ad(args.h5ad, backed="r")
    if obj.shape != (258120, 37795) or "counts" not in obj.layers:
        raise RuntimeError("Authoritative H5AD structure changed")
    obs_index = pd.Index(obj.obs_names.astype(str))
    var_index = pd.Index(obj.var_names.astype(str))

    cellchat_tokens = pd.read_csv(args.cellchat_genes, sep="\t")["gene"].astype(str).tolist()
    expression_genes = sorted(set(cellchat_tokens + TARGETED_GENES).intersection(var_index))
    if not set(TARGETED_GENES).issubset(expression_genes):
        raise RuntimeError(f"Missing targeted genes: {sorted(set(TARGETED_GENES)-set(expression_genes))}")
    sampled_pos = obs_index.get_indexer(sampled_meta["exact_cell_id"])
    gene_pos = var_index.get_indexer(expression_genes)
    if np.any(sampled_pos < 0) or np.any(gene_pos < 0):
        raise RuntimeError("Sampled CellChat identifiers failed H5AD join")
    sampled_counts = obj.layers["counts"][sampled_pos, :][:, gene_pos]
    sampled_counts = sp.csr_matrix(sampled_counts)
    validate_counts(sampled_counts, "sampled_cellchat")
    sampled_libraries = pd.to_numeric(obj.obs.iloc[sampled_pos]["total_counts_raw"]).to_numpy(float)
    sampled_log = normalize_log1p_cp10k(sampled_counts, sampled_libraries)
    matrix_part = args.output / "cellchat_expression_genes_by_cells_v1.mtx.gz.part"
    with gzip.open(matrix_part, "wb") as handle:
        mmwrite(handle, sampled_log.T.tocoo(), field="real", precision=10)
    os.replace(matrix_part, args.output / "cellchat_expression_genes_by_cells_v1.mtx.gz")
    atomic_tsv(pd.DataFrame({"gene": expression_genes}), args.output / "cellchat_genes_v1.tsv")
    atomic_tsv(sampled_meta, args.output / "cellchat_cell_metadata_v1.tsv")
    atomic_tsv(pd.DataFrame(sample_rows), args.output / "cellchat_sampling_qc_v1.tsv")

    bridge = bridge.set_index("exact_cell_id", drop=False)
    bridge_ids = bridge["exact_cell_id"].tolist()
    bridge_pos = obs_index.get_indexer(bridge_ids)
    if np.any(bridge_pos < 0):
        raise RuntimeError("Bridge identifiers failed H5AD join")
    full_counts = sp.csr_matrix(obj.layers["counts"][bridge_pos, :])
    validate_counts(full_counts, "malignant_full")
    full_libraries = pd.to_numeric(obj.obs.iloc[bridge_pos]["total_counts_raw"]).to_numpy(float)
    full_log = normalize_log1p_cp10k(full_counts, full_libraries)

    target_pos = var_index.get_indexer(TARGETED_GENES)
    target_values = full_log[:, target_pos].toarray()
    malignant_targeted = bridge.reset_index(drop=True)[["exact_cell_id", "patient", "study_v1", "F1", "TYMP_high", "TYMP_log1p_CP10K"]].copy()
    for j, gene in enumerate(TARGETED_GENES):
        malignant_targeted[f"{gene}_log1p_CP10K"] = target_values[:, j]
    if "X_umap" in obj.obsm:
        umap = np.asarray(obj.obsm["X_umap"][bridge_pos, :])
        malignant_targeted["UMAP1"] = umap[:, 0]
        malignant_targeted["UMAP2"] = umap[:, 1]
    atomic_tsv(malignant_targeted, args.output / "malignant_targeted_expression_umap_v1.tsv.gz", gzip_output=True)

    contrast_frames = []
    contrast_qc = []
    high_flag = bridge.reset_index(drop=True)["TYMP_high"].astype(bool).to_numpy()
    patient_vec = bridge.reset_index(drop=True)["patient"].astype(str).to_numpy()
    study_vec = bridge.reset_index(drop=True)["study_v1"].astype(str).to_numpy()
    for patient in malignant_patients:
        patient_mask = patient_vec == patient
        high = np.flatnonzero(patient_mask & high_flag)
        low = np.flatnonzero(patient_mask & ~high_flag)
        eligible_contrast = len(high) >= 20 and len(low) >= 20
        contrast_qc.append({"study": study_vec[np.flatnonzero(patient_mask)[0]], "patient": patient,
                            "n_tymp_high": len(high), "n_tymp_low": len(low),
                            "eligible_high_low_contrast": eligible_contrast})
        if not eligible_contrast:
            continue
        high_mean = np.asarray(full_log[high].mean(axis=0)).ravel()
        low_mean = np.asarray(full_log[low].mean(axis=0)).ravel()
        high_detect = np.asarray((full_counts[high] > 0).mean(axis=0)).ravel()
        low_detect = np.asarray((full_counts[low] > 0).mean(axis=0)).ravel()
        contrast_frames.append(pd.DataFrame({
            "study": study_vec[np.flatnonzero(patient_mask)[0]], "patient": patient,
            "gene": var_index, "n_tymp_high": len(high), "n_tymp_low": len(low),
            "mean_log_high": high_mean, "mean_log_low": low_mean,
            "high_minus_low": high_mean - low_mean,
            "detected_fraction_high": high_detect, "detected_fraction_low": low_detect,
        }))
    contrasts = pd.concat(contrast_frames, ignore_index=True)
    atomic_tsv(contrasts, args.output / "malignant_patient_gene_contrasts_v1.tsv.gz", gzip_output=True)
    atomic_tsv(pd.DataFrame(contrast_qc), args.output / "malignant_high_low_eligibility_v1.tsv")

    group_target_rows = []
    all_target_pos = obs_index.get_indexer(meta["exact_cell_id"])
    target_counts = sp.csr_matrix(obj.layers["counts"][all_target_pos, :][:, target_pos])
    validate_counts(target_counts, "group_targeted")
    target_libraries = pd.to_numeric(obj.obs.iloc[all_target_pos]["total_counts_raw"]).to_numpy(float)
    target_log = normalize_log1p_cp10k(target_counts, target_libraries)
    meta = meta.reset_index(drop=True)
    for (study, patient, group), idx in meta.groupby(["study_v1", "patient", "communication_group"], sort=True).groups.items():
        idx_array = np.asarray(list(idx), dtype=int)
        means = np.asarray(target_log[idx_array].mean(axis=0)).ravel()
        detect = np.asarray((target_counts[idx_array] > 0).mean(axis=0)).ravel()
        for j, gene in enumerate(TARGETED_GENES):
            group_target_rows.append({"study": study, "patient": patient, "communication_group": group,
                                      "gene": gene, "n_cells": len(idx_array), "mean_log1p_CP10K": means[j],
                                      "detected_fraction": detect[j]})
    atomic_tsv(pd.DataFrame(group_target_rows), args.output / "patient_group_targeted_expression_v1.tsv.gz", gzip_output=True)

    outputs = sorted(path.name for path in args.output.iterdir() if path.is_file())
    status = {
        "status": "PASS_C01_C04_INPUT_PREPARATION",
        "completed_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "input_hashes": hashes,
        "n_assignment_cells": len(assignments), "n_malignant_cells": len(bridge),
        "n_malignant_patients": len(malignant_patients), "n_cellchat_sampled_cells": len(sampled_meta),
        "n_cellchat_expression_genes": len(expression_genes), "sample_cap_per_patient_group": CAP,
        "n_high_low_eligible_patients": int(sum(row["eligible_high_low_contrast"] for row in contrast_qc)),
        "n_patient_gene_contrast_rows": len(contrasts), "outputs": outputs,
    }
    atomic_json(status, args.output / "status.json")


if __name__ == "__main__":
    main()
