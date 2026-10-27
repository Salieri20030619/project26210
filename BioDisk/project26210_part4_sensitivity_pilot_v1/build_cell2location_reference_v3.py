#!/usr/bin/env python3
"""Build study-balanced mean-count signatures for formal cell2location mapping."""

from __future__ import annotations

import hashlib
import json
import subprocess
from datetime import datetime
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix, issparse


ROOT = Path("/mnt/BioDisk/wangmingyuan/project26210_part4_sensitivity_pilot_v1")
INPUT = Path(
    "/mnt/BioDisk/ludingyu/project26210/results/part1_rebuild_v1/03_malignancy_cnv/"
    "cnv_annotated_h5ad_v1/part1_scvi_integrated_cell_annotated_cnv_v1.h5ad"
)
ASSIGNMENTS = ROOT / "reference_v2/reference_v2_cell_assignments.tsv.gz"
OLD_SIGNATURE = Path(
    "/mnt/BioDisk/ludingyu/project26210/results/part4_spatial/ref/ref_signatures_matrix.csv"
)
MEMORY_GATE = Path("/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh")
OUT = ROOT / "cell2location_reference_v3"


def gate(stage: str, rows: list[dict[str, object]]) -> None:
    result = subprocess.run([str(MEMORY_GATE)], capture_output=True, text=True, check=False)
    line = (result.stdout or result.stderr).strip().replace("\n", " | ")
    rows.append({"stage": stage, "exit_code": result.returncode, "output": line})
    if result.returncode != 0:
        raise RuntimeError(f"Memory gate failed at {stage}: {line}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=False)
    memory_rows: list[dict[str, object]] = []
    gate("before_load", memory_rows)

    assignments = pd.read_csv(ASSIGNMENTS, sep="\t", index_col="exact_cell_id")
    assignments.index = assignments.index.astype(str)
    if assignments.index.has_duplicates or assignments["state_v2"].nunique() != 14:
        raise RuntimeError("Invalid reference-v2 assignments")

    old = pd.read_csv(OLD_SIGNATURE, index_col=0)
    genes = old.index.astype(str).tolist()
    obj = ad.read_h5ad(INPUT, backed="r")
    missing_genes = pd.Index(genes).difference(obj.var_names)
    if len(missing_genes):
        raise RuntimeError(f"Missing old-reference genes: {len(missing_genes)}")
    cell_positions = obj.obs_names.get_indexer(assignments.index)
    if np.any(cell_positions < 0):
        raise RuntimeError(f"Missing assigned cells: {int(np.sum(cell_positions < 0))}")

    gate("before_counts", memory_rows)
    gene_positions = [int(obj.var_names.get_loc(gene)) for gene in genes]
    counts = obj.layers["counts"][cell_positions, :][:, gene_positions]
    if hasattr(counts, "to_memory"):
        counts = counts.to_memory()
    if not issparse(counts):
        counts = csr_matrix(np.asarray(counts))
    counts = counts.tocsr()
    if not np.isfinite(counts.data).all() or (counts.data < 0).any():
        raise RuntimeError("Invalid raw counts")
    if not np.allclose(counts.data, np.rint(counts.data)):
        raise RuntimeError("Counts are not integer-like")

    assignments["patient_state"] = (
        assignments["patient"].astype(str) + "||" + assignments["state_v2"].astype(str)
    )
    patient_states = sorted(assignments["patient_state"].unique())
    codes = pd.Categorical(assignments["patient_state"], categories=patient_states).codes
    indicator = csr_matrix(
        (np.ones(len(assignments)), (codes, np.arange(len(assignments)))),
        shape=(len(patient_states), len(assignments)),
    )
    pseudobulk = indicator @ counts
    n_cells = np.asarray(indicator.sum(axis=1)).reshape(-1)
    patient_mean = pseudobulk.multiply((1.0 / np.maximum(n_cells, 1))[:, None]).toarray()

    meta = pd.DataFrame({"patient_state": patient_states, "n_cells": n_cells})
    meta[["patient", "state_v2"]] = meta["patient_state"].str.split(
        "||", n=1, expand=True, regex=False
    )
    patient_study = assignments.drop_duplicates("patient").set_index("patient")["study_v1"]
    meta["study"] = meta["patient"].map(patient_study)
    if meta["study"].isna().any():
        raise RuntimeError("Study mapping incomplete")

    signatures: dict[str, np.ndarray] = {}
    qc_rows: list[dict[str, object]] = []
    for state, state_rows in meta.groupby("state_v2", observed=True):
        study_profiles = []
        for _, study_rows in state_rows.groupby("study", observed=True):
            study_profiles.append(patient_mean[study_rows.index].mean(axis=0))
        profile = np.vstack(study_profiles).mean(axis=0)
        signatures[state] = profile
        qc_rows.append({
            "state_v2": state,
            "n_cells": int(state_rows["n_cells"].sum()),
            "n_patients": int(state_rows["patient"].nunique()),
            "n_studies": int(state_rows["study"].nunique()),
            "signature_sum_mean_counts": float(profile.sum()),
            "signature_nonzero_genes": int(np.sum(profile > 0)),
        })

    signature = pd.DataFrame(signatures, index=genes)
    if signature.shape != (12758, 14) or not np.isfinite(signature.to_numpy()).all():
        raise RuntimeError(f"Invalid signature matrix: {signature.shape}")
    if (signature.to_numpy() < 0).any() or signature.columns.has_duplicates:
        raise RuntimeError("Invalid signature values or identifiers")

    signature_path = OUT / "reference_v3_signatures_study_balanced_mean_counts.tsv.gz"
    signature.to_csv(signature_path, sep="\t")
    qc = pd.DataFrame(qc_rows).sort_values("state_v2")
    qc.to_csv(OUT / "reference_v3_state_qc.tsv", sep="\t", index=False)
    meta.to_csv(OUT / "reference_v3_patient_state_qc.tsv", sep="\t", index=False)
    pd.DataFrame(memory_rows).to_csv(OUT / "memory_gates_v3.tsv", sep="\t", index=False)

    metadata = {
        "status": "completed",
        "created_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "input_h5ad": str(INPUT),
        "input_h5ad_verified_sha256": "e4779c46af13a8643e8dec584eeff4322258181f5edab3e55dcd8c20bc28a9e7",
        "assignments": str(ASSIGNMENTS),
        "assignments_sha256": sha256(ASSIGNMENTS),
        "n_cells": int(len(assignments)),
        "n_genes": int(signature.shape[0]),
        "n_states": int(signature.shape[1]),
        "aggregation": "patient mean raw counts -> study mean -> equal-weight study mean",
        "purpose": "formal GLMF1/GLMF2 cell2location reference sensitivity mapping",
    }
    (OUT / "run_metadata_v3.json").write_text(
        json.dumps(metadata, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    if getattr(obj, "file", None) is not None:
        obj.file.close()

    manifest = []
    for path in sorted(OUT.iterdir()):
        if path.is_file() and path.name != "file_manifest_v3.tsv":
            manifest.append({
                "file": path.name,
                "size_bytes": path.stat().st_size,
                "sha256": sha256(path),
            })
    pd.DataFrame(manifest).to_csv(OUT / "file_manifest_v3.tsv", sep="\t", index=False)
    print(json.dumps(metadata, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
