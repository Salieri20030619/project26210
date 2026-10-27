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


SEED = 20260814
MAX_PER_PATIENT_STATE = 100
STATES = {"Monocyte", "Other Macrophage", "C1QC-state Macrophage", "SPP1-state Macrophage"}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def atomic_tsv(frame: pd.DataFrame, path: Path, compression: str | None = None) -> None:
    part = Path(str(path) + ".part")
    frame.to_csv(part, sep="\t", index=False, compression=compression, lineterminator="\n")
    os.replace(part, path)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", type=Path, required=True)
    parser.add_argument("--dpt-cell-scores", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    args.output.mkdir(parents=True)

    dpt = pd.read_csv(args.dpt_cell_scores, sep="\t", dtype={"cell_id": str})
    if len(dpt) != 50761 or dpt["cell_id"].duplicated().any() or set(dpt["state_v2"]) != STATES:
        raise RuntimeError("Unexpected validated myeloid scope")
    adata = ad.read_h5ad(args.h5ad, backed="r")
    if "counts" not in adata.layers or not adata.obs_names.is_unique:
        raise RuntimeError("Authoritative H5AD lacks raw counts or unique cell IDs")
    index = adata.obs_names.astype(str)
    selector = index.isin(set(dpt["cell_id"]))
    rows = np.flatnonzero(selector)
    cells = index[rows].to_numpy(str)
    if len(rows) != 50761:
        raise RuntimeError(f"Expected 50,761 H5AD cells, observed {len(rows)}")
    genes = adata.var_names.astype(str).to_numpy()
    counts = sp.csr_matrix(adata.layers["counts"][rows, :], dtype=np.float64)
    adata.file.close()
    if counts.shape != (50761, len(genes)) or np.any(~np.isfinite(counts.data)) or np.any(counts.data < 0) or np.any(np.abs(counts.data - np.rint(counts.data)) > 1e-8):
        raise RuntimeError("Raw-count matrix failed shape or integer-safety checks")
    meta = pd.DataFrame({"cell_id": cells}).merge(dpt, on="cell_id", how="left", validate="one_to_one")
    if meta["state_v2"].isna().any():
        raise RuntimeError("Cell metadata alignment failed")

    rng = np.random.default_rng(SEED)
    selected = []
    sampling_rows = []
    for (patient, state), frame in meta.groupby(["patient", "state_v2"], sort=True, observed=True):
        idx = frame.index.to_numpy(int)
        chosen = idx if len(idx) <= MAX_PER_PATIENT_STATE else np.sort(rng.choice(idx, MAX_PER_PATIENT_STATE, replace=False))
        selected.extend(chosen.tolist())
        sampling_rows.append({"patient": patient, "state_v2": state, "n_all_cells": len(idx), "n_discovery_cells": len(chosen)})
    selected = np.asarray(sorted(set(selected)), dtype=int)
    if len(selected) < 5000 or len(selected) >= len(meta):
        raise RuntimeError(f"Unexpected balanced discovery size: {len(selected)}")
    meta["selected_for_grn_discovery_v1"] = False
    meta.loc[selected, "selected_for_grn_discovery_v1"] = True

    with gzip.open(args.output / "raw_counts_genes_by_cells_v1.mtx.gz.part", "wb") as handle:
        mmwrite(handle, counts.T.tocsc(), field="integer")
    os.replace(args.output / "raw_counts_genes_by_cells_v1.mtx.gz.part", args.output / "raw_counts_genes_by_cells_v1.mtx.gz")
    atomic_tsv(pd.DataFrame({"gene": genes}), args.output / "genes_v1.tsv")
    atomic_tsv(pd.DataFrame({"cell_id": cells}), args.output / "cells_v1.tsv")
    atomic_tsv(meta, args.output / "cell_metadata_v1.tsv.gz", compression="gzip")
    atomic_tsv(pd.DataFrame(sampling_rows), args.output / "balanced_sampling_qc_v1.tsv")
    detected = np.asarray((counts > 0).sum(axis=0)).ravel()
    totals = np.asarray(counts.sum(axis=0)).ravel()
    atomic_tsv(pd.DataFrame({"gene": genes, "detected_cells": detected, "total_counts": totals}), args.output / "raw_gene_qc_v1.tsv.gz", compression="gzip")

    status = {
        "status": "PASS_Y03_PYSCENIC_INPUT_PREPARATION",
        "completed_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "n_cells": len(meta),
        "n_genes": len(genes),
        "nnz": int(counts.nnz),
        "n_studies": int(meta["study_v1"].nunique()),
        "n_patients": int(meta["patient"].nunique()),
        "n_discovery_cells": len(selected),
        "sampling": "maximum 100 cells per patient-by-state stratum; all smaller strata retained",
        "seed": SEED,
        "input_h5ad_sha256": sha256(args.h5ad),
        "dpt_cell_scores_sha256": sha256(args.dpt_cell_scores),
    }
    part = args.output / "status.json.part"
    part.write_text(json.dumps(status, indent=2) + "\n", encoding="utf-8")
    os.replace(part, args.output / "status.json")
    print(json.dumps(status, indent=2), flush=True)


if __name__ == "__main__":
    main()
