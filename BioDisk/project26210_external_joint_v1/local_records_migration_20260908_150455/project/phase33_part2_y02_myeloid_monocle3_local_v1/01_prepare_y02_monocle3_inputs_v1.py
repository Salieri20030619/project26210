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
from scipy import sparse
from scipy.io import mmwrite


STATES = {"Monocyte", "Other Macrophage", "C1QC-state Macrophage", "SPP1-state Macrophage"}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def atomic_tsv(frame: pd.DataFrame, path: Path, compression: str | None = None, header: bool = True) -> None:
    part = Path(str(path) + ".part")
    frame.to_csv(part, sep="\t", index=False, header=header, compression=compression, lineterminator="\n")
    os.replace(part, path)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", type=Path, required=True)
    parser.add_argument("--dpt-cell-scores", type=Path, required=True)
    parser.add_argument("--hvg", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    args.output.mkdir(parents=True)

    dpt = pd.read_csv(args.dpt_cell_scores, sep="\t", dtype={"cell_id": str})
    if len(dpt) != 50761 or dpt["cell_id"].duplicated().any() or set(dpt["state_v2"]) != STATES:
        raise RuntimeError("Unexpected validated Phase 21 DPT scope")
    hvg = pd.read_csv(args.hvg, sep="\t")["gene"].astype(str).tolist()
    if len(hvg) != 2500 or len(set(hvg)) != 2500:
        raise RuntimeError("Expected 2,500 unique frozen HVGs")

    obj = ad.read_h5ad(args.h5ad)
    if "counts" not in obj.layers or "X_scVI" not in obj.obsm or not obj.obs_names.is_unique:
        raise RuntimeError("Authoritative H5AD lacks required counts or X_scVI")
    index = obj.obs_names.astype(str)
    if not set(dpt["cell_id"]).issubset(set(index)):
        raise RuntimeError("DPT cells missing from H5AD")
    gene_index = obj.var_names.get_indexer(hvg)
    if (gene_index < 0).any():
        raise RuntimeError("Frozen HVGs missing from H5AD")
    cell_index = pd.Series(np.arange(obj.n_obs), index=index)

    study_rows = []
    for study in sorted(dpt["study_v1"].unique()):
        frame = dpt.loc[dpt["study_v1"].eq(study)].sort_values("cell_id", kind="stable").reset_index(drop=True)
        rows = cell_index.loc[frame["cell_id"]].to_numpy(int)
        counts = obj.layers["counts"][rows, :][:, gene_index]
        counts = sparse.csr_matrix(counts)
        if np.any(~np.isfinite(counts.data)) or np.any(counts.data < 0) or np.any(np.abs(counts.data - np.rint(counts.data)) > 1e-8):
            raise RuntimeError(f"{study}: counts are not finite nonnegative integers")
        latent = np.asarray(obj.obsm["X_scVI"])[rows]
        if np.any(~np.isfinite(latent)):
            raise RuntimeError(f"{study}: non-finite X_scVI")
        out = args.output / study
        out.mkdir()
        mtx_part = out / "raw_counts_genes_by_cells_v1.mtx.gz.part"
        with gzip.open(mtx_part, "wb") as handle:
            mmwrite(handle, counts.T.tocoo(), field="integer")
        os.replace(mtx_part, out / "raw_counts_genes_by_cells_v1.mtx.gz")
        atomic_tsv(pd.DataFrame({"gene": hvg}), out / "genes_v1.tsv", header=False)
        atomic_tsv(frame, out / "cell_metadata_v1.tsv.gz", compression="gzip")
        latent_frame = pd.DataFrame(latent, columns=[f"scVI_{i+1}" for i in range(latent.shape[1])])
        latent_frame.insert(0, "cell_id", frame["cell_id"])
        atomic_tsv(latent_frame, out / "scvi_latent_v1.tsv.gz", compression="gzip")
        study_rows.append({
            "study": study,
            "n_cells": len(frame),
            "n_patients": int(frame["patient"].nunique()),
            "n_monocyte": int(frame["state_v2"].eq("Monocyte").sum()),
            "n_spp1_state": int(frame["state_v2"].eq("SPP1-state Macrophage").sum()),
            "n_genes": len(hvg),
            "n_nonzero_counts": int(counts.nnz),
        })
    atomic_tsv(pd.DataFrame(study_rows), args.output / "study_input_qc_v1.tsv")
    files = []
    for path in sorted(p for p in args.output.rglob("*") if p.is_file()):
        files.append({"relative_path": str(path.relative_to(args.output)), "bytes": path.stat().st_size, "sha256": sha256(path)})
    atomic_tsv(pd.DataFrame(files), args.output / "input_manifest_v1.tsv")
    status = {
        "status": "PASS_Y02_MONOCLE3_INPUT_PREPARATION",
        "completed_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "n_cells": len(dpt),
        "n_studies": int(dpt["study_v1"].nunique()),
        "n_patients": int(dpt["patient"].nunique()),
        "n_hvg": len(hvg),
        "study_separate": True,
        "input_h5ad_sha256": sha256(args.h5ad),
        "dpt_cell_scores_sha256": sha256(args.dpt_cell_scores),
        "hvg_sha256": sha256(args.hvg),
    }
    part = args.output / "status.json.part"
    part.write_text(json.dumps(status, indent=2) + "\n", encoding="utf-8")
    os.replace(part, args.output / "status.json")
    print(json.dumps(status, indent=2), flush=True)


if __name__ == "__main__":
    main()
