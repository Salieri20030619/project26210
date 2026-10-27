#!/usr/bin/env python3
"""Create full, balanced and rank7-focused loom inputs for pySCENIC."""

from __future__ import annotations

import json
import subprocess
from datetime import datetime
from pathlib import Path

import loompy
import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy.sparse import csc_matrix


PROJECT = Path("/mnt/BioDisk/ludingyu/project26210")
ROOT = PROJECT / "results/part1_rebuild_v1/07_pyscenic_v1"
INPUT = ROOT / "input"
TABLES = ROOT / "tables"
REFERENCE = PROJECT / "references/pyscenic/hg38"
GATE = PROJECT / "check_memory_gate.sh"
SEED = 26210
MAX_CELLS_PER_SAMPLE = 200
MIN_FACTOR_CELLS = 5


def gate(stage: str) -> str:
    result = subprocess.run([str(GATE), stage], text=True, capture_output=True)
    output = " ".join((result.stdout + " " + result.stderr).split())
    print(output, flush=True)
    if result.returncode != 0 or "PASS" not in output:
        raise RuntimeError(f"Memory gate blocked {stage}: {output}")
    return output


def balanced_indices(meta: pd.DataFrame) -> np.ndarray:
    rng = np.random.default_rng(SEED)
    selected: list[int] = []
    factor_col = "dominant_factor_rank7_v1"
    for _, sample_df in meta.groupby("sample_id", sort=True):
        indices = sample_df.index.to_numpy()
        if len(indices) <= MAX_CELLS_PER_SAMPLE:
            selected.extend(indices.tolist())
            continue
        chosen: list[int] = []
        for _, factor_df in sample_df.groupby(factor_col, sort=True):
            proportional = int(np.floor(MAX_CELLS_PER_SAMPLE * len(factor_df) / len(sample_df)))
            quota = min(len(factor_df), max(MIN_FACTOR_CELLS, proportional))
            chosen.extend(rng.choice(factor_df.index.to_numpy(), quota, replace=False).tolist())
        chosen = list(dict.fromkeys(chosen))
        if len(chosen) > MAX_CELLS_PER_SAMPLE:
            chosen = rng.choice(np.asarray(chosen), MAX_CELLS_PER_SAMPLE, replace=False).tolist()
        elif len(chosen) < MAX_CELLS_PER_SAMPLE:
            remaining = np.setdiff1d(indices, np.asarray(chosen))
            chosen.extend(rng.choice(remaining, MAX_CELLS_PER_SAMPLE - len(chosen), replace=False).tolist())
        selected.extend(chosen)
    return np.asarray(sorted(set(selected)), dtype=int)


def create_loom(path: Path, matrix: csc_matrix, genes: np.ndarray, meta: pd.DataFrame) -> None:
    attrs = {
        "CellID": meta["exact_cell_id"].astype(str).to_numpy(dtype=object),
        "sample_id": meta["sample_id"].astype(str).to_numpy(dtype=object),
        "patient_id": meta["curated_patient_id"].astype(str).to_numpy(dtype=object),
        "study": meta["study"].astype(str).to_numpy(dtype=object),
        "dominant_factor": meta["dominant_factor_rank7_v1"].astype(str).to_numpy(dtype=object),
        "F1_score": meta["F1"].to_numpy(float),
        "TYMP_log1p_CP10K": meta["TYMP_log1p_CP10K_v1"].to_numpy(float),
    }
    loompy.create(str(path), matrix.astype(np.float32), row_attrs={"Gene": genes}, col_attrs=attrs)


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    expected = [
        INPUT / "raw_counts_genes_by_cells_v1.mtx.gz",
        INPUT / "genes_v1.tsv",
        INPUT / "cells_v1.tsv",
        INPUT / "cell_metadata_v1.tsv.gz",
        INPUT / "cnmf_rank7_model_features_v1.tsv",
    ]
    if any(not p.exists() for p in expected):
        raise FileNotFoundError([str(p) for p in expected if not p.exists()])
    outputs = [
        INPUT / "all_cells_filtered_counts_v1.loom",
        INPUT / "balanced_grn_counts_v1.loom",
        INPUT / "focused_genie3_counts_v1.loom",
    ]
    if any(p.exists() for p in outputs):
        raise RuntimeError("Refusing to overwrite existing loom input")

    gate_log = [{"stage": "pyscenic_loom_read_matrix_v1", "gate_output": gate("pyscenic_loom_read_matrix_v1") }]
    raw = csc_matrix(mmread(INPUT / "raw_counts_genes_by_cells_v1.mtx.gz"))
    genes = pd.read_csv(INPUT / "genes_v1.tsv", sep="\t")["gene"].astype(str).to_numpy()
    cells = pd.read_csv(INPUT / "cells_v1.tsv", sep="\t")["exact_cell_id"].astype(str).to_numpy()
    meta = pd.read_csv(INPUT / "cell_metadata_v1.tsv.gz", sep="\t", dtype={"exact_cell_id": str})
    if raw.shape != (37795, 11950) or raw.nnz != 35713877:
        raise RuntimeError(f"Unexpected matrix: {raw.shape}, {raw.nnz}")
    if not np.array_equal(meta["exact_cell_id"].to_numpy(str), cells):
        raise RuntimeError("Cell metadata order mismatch")

    detected = np.asarray((raw > 0).sum(axis=1)).ravel()
    totals = np.asarray(raw.sum(axis=1)).ravel()
    min_cells = max(20, int(np.ceil(0.01 * raw.shape[1])))
    min_counts = 3 * min_cells
    upper = np.char.upper(genes.astype(str))
    mito = np.char.startswith(upper, "MT-") | np.char.startswith(upper, "MTRNR")
    keep = (detected >= min_cells) & (totals >= min_counts) & (~mito)
    filtered = raw[keep, :].tocsc()
    filtered_genes = genes[keep]
    inference_idx = balanced_indices(meta)
    if len(inference_idx) < 3000 or len(inference_idx) > 5800:
        raise RuntimeError(f"Unexpected balanced GRN cell count: {len(inference_idx)}")

    gate_log.append({"stage": "pyscenic_loom_write_v1", "gate_output": gate("pyscenic_loom_write_v1")})
    create_loom(INPUT / "all_cells_filtered_counts_v1.loom", filtered, filtered_genes, meta)
    create_loom(
        INPUT / "balanced_grn_counts_v1.loom",
        filtered[:, inference_idx], filtered_genes, meta.iloc[inference_idx].reset_index(drop=True),
    )

    features = set(pd.read_csv(INPUT / "cnmf_rank7_model_features_v1.tsv", sep="\t", header=None)[0].astype(str))
    tfs = {line.strip() for line in (REFERENCE / "allTFs_hg38.txt").read_text().splitlines() if line.strip()}
    focused_keep = np.asarray([g in features or g in tfs or g in {"TYMP", "HIF1A"} for g in filtered_genes])
    create_loom(
        INPUT / "focused_genie3_counts_v1.loom",
        filtered[focused_keep, :][:, inference_idx].tocsc(),
        filtered_genes[focused_keep], meta.iloc[inference_idx].reset_index(drop=True),
    )

    qc = pd.DataFrame({
        "metric": ["n_input_genes", "n_filtered_genes", "n_all_cells", "n_balanced_grn_cells", "n_focused_genie3_genes", "max_cells_per_sample", "min_factor_cells", "seed"],
        "value": [len(genes), len(filtered_genes), len(cells), len(inference_idx), int(focused_keep.sum()), MAX_CELLS_PER_SAMPLE, MIN_FACTOR_CELLS, SEED],
    })
    qc.to_csv(TABLES / "pyscenic_input_qc_v1.tsv", sep="\t", index=False)
    infer_meta = meta.iloc[inference_idx].copy()
    infer_meta.to_csv(TABLES / "grn_inference_cells_v1.tsv.gz", sep="\t", index=False, compression="gzip")
    (
        meta.groupby(["study", "sample_id", "dominant_factor_rank7_v1"], observed=False).size().rename("n_all_cells").reset_index()
        .merge(
            infer_meta.groupby(["study", "sample_id", "dominant_factor_rank7_v1"], observed=False).size().rename("n_grn_cells").reset_index(),
            on=["study", "sample_id", "dominant_factor_rank7_v1"], how="left",
        ).fillna({"n_grn_cells": 0}).to_csv(TABLES / "grn_sampling_summary_v1.tsv", sep="\t", index=False)
    )
    pd.DataFrame(gate_log).to_csv(TABLES / "memory_gate_loom_v1.tsv", sep="\t", index=False)
    (TABLES / "loom_run_metadata_v1.json").write_text(json.dumps({
        "completed_at": datetime.now().astimezone().isoformat(), "seed": SEED,
        "max_cells_per_sample": MAX_CELLS_PER_SAMPLE, "n_grn_cells": int(len(inference_idx)),
    }, indent=2) + "\n")
    (INPUT / "LOOMS_COMPLETE").write_text(datetime.now().astimezone().isoformat() + "\n")
    print(qc.to_string(index=False), flush=True)
    print("PART1_PYSCENIC_LOOMS_V1_PASS", flush=True)


if __name__ == "__main__":
    main()
