#!/usr/bin/env python3
"""Build a 13-state reference with one total malignant state.

TYMP and F1 are retained as continuous malignant-cell annotations. This avoids
asking Cell2location to distinguish the non-identifiable TYMP-high/low pair.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import subprocess
from datetime import datetime
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy
from scipy.sparse import csr_matrix, issparse


OLD_MALIGNANT = {"F1-high Malignant", "F1-low Malignant"}
EXPECTED_STATES = sorted([
    "B cell", "C1QC-state Macrophage", "Dendritic cell", "Endothelial",
    "Epithelial-Normal", "Malignant", "Mast cell", "Monocyte", "NK cell",
    "Other Macrophage", "SMC/Pericyte", "SPP1-state Macrophage", "T cell",
])


def args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--h5ad", type=Path, required=True)
    p.add_argument("--old-assignments", type=Path, required=True)
    p.add_argument("--old-signature", type=Path, required=True)
    p.add_argument("--bridge", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--memory-gate", type=Path)
    return p.parse_args()


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def atomic_tsv(frame: pd.DataFrame, path: Path, index: bool = False) -> None:
    part = path.with_name(path.name + ".part")
    frame.to_csv(part, sep="\t", index=index, encoding="utf-8", na_rep="NA",
                 compression="gzip" if path.name.endswith(".gz") else None,
                 lineterminator="\n")
    os.replace(part, path)


def gate(path: Path | None, stage: str, rows: list[dict[str, object]]) -> None:
    if path is None:
        rows.append({"stage": stage, "exit_code": 0, "output": "not_requested"})
        return
    result = subprocess.run([str(path)], capture_output=True, text=True, check=False)
    message = (result.stdout or result.stderr).strip().replace("\n", " | ")
    rows.append({"stage": stage, "exit_code": result.returncode, "output": message})
    if result.returncode:
        raise RuntimeError(f"Memory gate failed at {stage}: {message}")


def main() -> None:
    a = args()
    for path in [a.h5ad, a.old_assignments, a.old_signature, a.bridge]:
        if not path.is_file():
            raise FileNotFoundError(path)
    if a.output.exists():
        raise FileExistsError(f"Refusing to overwrite: {a.output}")
    work = a.output.with_name(a.output.name + ".part")
    if work.exists():
        raise FileExistsError(f"Stale work directory: {work}")
    work.mkdir(parents=True)
    gates: list[dict[str, object]] = []
    gate(a.memory_gate, "before_load", gates)

    old_signature = pd.read_csv(a.old_signature, sep="\t", index_col=0)
    old_signature.index = old_signature.index.astype(str)
    if old_signature.shape != (12758, 14):
        raise RuntimeError(f"Unexpected old signature: {old_signature.shape}")
    genes = old_signature.index.tolist()

    assignments = pd.read_csv(a.old_assignments, sep="\t", index_col="exact_cell_id")
    assignments.index = assignments.index.astype(str)
    if assignments.index.has_duplicates or assignments["state_v2"].nunique() != 14:
        raise RuntimeError("Invalid frozen reference-v2 assignments")
    assignments["state_tymp_continuous_v1"] = assignments["state_v2"].replace(
        {state: "Malignant" for state in OLD_MALIGNANT}
    )
    if sorted(assignments["state_tymp_continuous_v1"].unique()) != EXPECTED_STATES:
        raise RuntimeError("Unexpected collapsed state names")

    bridge = pd.read_csv(a.bridge, sep="\t", index_col="exact_cell_id")
    bridge.index = bridge.index.astype(str)
    malignant = assignments["state_tymp_continuous_v1"].eq("Malignant")
    if int(malignant.sum()) != 11921 or set(assignments.index[malignant]) != set(bridge.index):
        raise RuntimeError("Collapsed malignant cells do not match the frozen Phase 17 bridge")
    for column in ["TYMP_raw_count", "TYMP_log1p_CP10K", "TYMP_q75_patient"]:
        assignments[column] = bridge[column].reindex(assignments.index).to_numpy()

    obj = ad.read_h5ad(a.h5ad, backed="r")
    cell_positions = obj.obs_names.get_indexer(assignments.index)
    if np.any(cell_positions < 0):
        raise RuntimeError(f"Assigned cells missing from H5AD: {int(np.sum(cell_positions < 0))}")
    missing_genes = pd.Index(genes).difference(obj.var_names)
    if len(missing_genes):
        raise RuntimeError(f"Genes missing from H5AD: {len(missing_genes)}")
    gate(a.memory_gate, "before_counts", gates)
    gene_positions = [int(obj.var_names.get_loc(g)) for g in genes]
    counts = obj.layers["counts"][cell_positions, :][:, gene_positions]
    if hasattr(counts, "to_memory"):
        counts = counts.to_memory()
    if not issparse(counts):
        counts = csr_matrix(np.asarray(counts))
    counts = counts.tocsr()
    if not np.isfinite(counts.data).all() or (counts.data < 0).any():
        raise RuntimeError("Counts contain negative or non-finite values")
    if not np.allclose(counts.data, np.rint(counts.data)):
        raise RuntimeError("Counts are not integer-like; refusing silent repair")

    assignments["patient_state"] = (
        assignments["patient"].astype(str) + "||" + assignments["state_tymp_continuous_v1"].astype(str)
    )
    patient_states = sorted(assignments["patient_state"].unique())
    codes = pd.Categorical(assignments["patient_state"], categories=patient_states).codes
    indicator = csr_matrix((np.ones(len(assignments)), (codes, np.arange(len(assignments)))),
                           shape=(len(patient_states), len(assignments)))
    n_cells = np.asarray(indicator.sum(axis=1)).reshape(-1)
    patient_mean = (indicator @ counts).multiply((1 / np.maximum(n_cells, 1))[:, None]).toarray()
    meta = pd.DataFrame({"patient_state": patient_states, "n_cells": n_cells.astype(int)})
    meta[["patient", "state_tymp_continuous_v1"]] = meta["patient_state"].str.split(
        "||", n=1, expand=True, regex=False
    )
    patient_study = assignments.drop_duplicates("patient").set_index("patient")["study_v1"]
    meta["study"] = meta["patient"].map(patient_study)
    if meta["study"].isna().any():
        raise RuntimeError("Incomplete patient-study mapping")

    profiles: dict[str, np.ndarray] = {}
    qc_rows: list[dict[str, object]] = []
    for state, rows in meta.groupby("state_tymp_continuous_v1", observed=True):
        study_means = [patient_mean[s.index].mean(axis=0) for _, s in rows.groupby("study", observed=True)]
        profile = np.vstack(study_means).mean(axis=0)
        profiles[state] = profile
        qc_rows.append({
            "state_tymp_continuous_v1": state,
            "n_cells": int(rows["n_cells"].sum()),
            "n_patients": int(rows["patient"].nunique()),
            "n_studies": int(rows["study"].nunique()),
            "signature_sum_mean_counts": float(profile.sum()),
            "signature_nonzero_genes": int(np.sum(profile > 0)),
        })
    signature = pd.DataFrame(profiles, index=genes).reindex(columns=EXPECTED_STATES)
    values = signature.to_numpy()
    if signature.shape != (12758, 13) or not np.isfinite(values).all() or (values < 0).any():
        raise RuntimeError(f"Invalid continuous-TYMP reference: {signature.shape}")

    nonmalignant = sorted(set(EXPECTED_STATES).difference({"Malignant"}))
    max_delta = float(np.abs(signature[nonmalignant].to_numpy() - old_signature[nonmalignant].to_numpy()).max())
    if not np.allclose(signature[nonmalignant], old_signature[nonmalignant], rtol=1e-10, atol=1e-12):
        raise RuntimeError(f"Non-malignant signatures drifted: {max_delta}")
    malignant_qc = next(row for row in qc_rows if row["state_tymp_continuous_v1"] == "Malignant")
    if malignant_qc["n_cells"] != 11921 or malignant_qc["n_patients"] != 20 or malignant_qc["n_studies"] != 5:
        raise RuntimeError(f"Unexpected malignant coverage: {malignant_qc}")

    validation = pd.DataFrame([
        {"check": "signature_shape", "pass": True, "observed": "12758x13", "expected": "12758x13"},
        {"check": "malignant_cells", "pass": True, "observed": malignant_qc["n_cells"], "expected": 11921},
        {"check": "malignant_patients", "pass": True, "observed": malignant_qc["n_patients"], "expected": 20},
        {"check": "malignant_studies", "pass": True, "observed": malignant_qc["n_studies"], "expected": 5},
        {"check": "nonmalignant_signatures_unchanged", "pass": max_delta <= 1e-12, "observed": max_delta, "expected": "<=1e-12"},
        {"check": "raw_counts_integer_nonnegative_finite", "pass": True, "observed": int(counts.nnz), "expected": "PASS"},
    ])

    keep = ["patient", "study_v1", "state_tymp_continuous_v1", "F1", "TYMP_raw_count",
            "TYMP_log1p_CP10K", "TYMP_q75_patient", "SPP1_program_score_v2", "C1QC_program_score_v2"]
    assignments.index.name = "exact_cell_id"
    atomic_tsv(signature, work / "reference_tymp_continuous_v1_signatures_study_balanced_mean_counts.tsv.gz", True)
    atomic_tsv(assignments[keep], work / "reference_tymp_continuous_v1_cell_assignments.tsv.gz", True)
    atomic_tsv(pd.DataFrame(qc_rows).sort_values("state_tymp_continuous_v1"), work / "reference_tymp_continuous_v1_state_qc.tsv")
    atomic_tsv(meta, work / "reference_tymp_continuous_v1_patient_state_qc.tsv")
    atomic_tsv(validation, work / "reference_tymp_continuous_v1_internal_validation.tsv")
    atomic_tsv(pd.DataFrame(gates), work / "memory_gates_tymp_continuous_v1.tsv")
    metadata = {
        "schema_version": "phase17_reference_tymp_continuous_v1", "status": "PASS_INTERNAL_VALIDATION",
        "created_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "inputs": {"h5ad": str(a.h5ad), "old_assignments": str(a.old_assignments),
                   "old_assignments_sha256": sha256(a.old_assignments), "old_signature": str(a.old_signature),
                   "old_signature_sha256": sha256(a.old_signature), "bridge": str(a.bridge),
                   "bridge_sha256": sha256(a.bridge)},
        "n_reference_cells": int(len(assignments)), "n_genes": 12758, "n_states": 13,
        "malignant_definition": "11921 frozen malignant bridge cells collapsed to one identifiable Cell2location state",
        "tymp_role": "continuous malignant exposure; not a Cell2location subtype",
        "f1_role": "continuous auxiliary covariate; not a TYMP label",
        "aggregation": "patient mean raw counts -> study mean -> equal-weight study mean",
        "failed_two_state_attempt": {"reason": "TYMP-high/low signatures not identifiable", "spearman_rho": 0.9894368522597528},
        "versions": {"python": platform.python_version(), "anndata": ad.__version__, "numpy": np.__version__,
                     "pandas": pd.__version__, "scipy": scipy.__version__},
    }
    part = work / "run_metadata_tymp_continuous_v1.json.part"
    part.write_text(json.dumps(metadata, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    os.replace(part, work / "run_metadata_tymp_continuous_v1.json")
    if getattr(obj, "file", None) is not None:
        obj.file.close()
    manifest_rows = []
    for path in sorted(work.iterdir()):
        if path.is_file() and path.name != "file_manifest_tymp_continuous_v1.tsv":
            manifest_rows.append({"file": path.name, "size_bytes": path.stat().st_size, "sha256": sha256(path)})
    atomic_tsv(pd.DataFrame(manifest_rows), work / "file_manifest_tymp_continuous_v1.tsv")
    os.replace(work, a.output)
    print(json.dumps(metadata, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
