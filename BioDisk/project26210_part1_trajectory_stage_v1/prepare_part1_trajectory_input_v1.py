#!/usr/bin/env python3
from __future__ import annotations

import gzip
import hashlib
import json
import subprocess
from datetime import datetime
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmwrite


SEED = 26210
PROJECT = Path("/mnt/BioDisk/ludingyu/project26210")
ROOT = PROJECT / "results/part1_rebuild_v1"
OUT = ROOT / "09_metabolism_trajectory/trajectory_v1"
INPUT = OUT / "input"
QC = OUT / "qc"
H5AD = ROOT / "03_malignancy_cnv/cnv_annotated_h5ad_v1/part1_scvi_integrated_cell_annotated_cnv_v1.h5ad"
META = ROOT / "05_cnmf/input_v1/cnmf_cell_metadata_v1.tsv.gz"
HVG = ROOT / "05_cnmf/input_v1/cnmf_hvg_genes_v1.tsv"
USAGE = ROOT / "05_cnmf/full_ranks4_8_v1/consensus_usage_rank7_v1.tsv.gz"
CLONES = ROOT / "04_subclones_tymp/cnv_subclones_v1/formal_cell_cnv_clone_assignments_v1.tsv.gz"
AUC = ROOT / "09_metabolism_trajectory/metabolism_v1/source_data/four_pathway_AUCell_cell_scores_v1.tsv.gz"
SCFEA = ROOT / "09_metabolism_trajectory/metabolism_v1/source_data/scFEA_four_axis_cell_scores_v1.tsv.gz"
GATE = PROJECT / "check_memory_gate.sh"

STATE_GENES = {
    "Base": ["PAX8", "EPCAM", "KRT8", "KRT18", "CDH1", "MIOX", "ALDOB", "FBP1", "GPX3", "SLC34A1"],
    "Hypoxia": ["CA9", "NDUFA4L2", "VEGFA", "EGLN3", "BNIP3", "ENO1", "LDHA", "SLC2A1", "ADM", "HILPDA"],
    "Cycle": ["MKI67", "TOP2A", "UBE2C", "CENPF", "CCNB1", "CCNB2", "CDK1", "TK1", "RRM2", "TYMS"],
    "TYMP_state": ["TYMP", "SPP1", "HSPB1", "LGALS1", "FABP7", "PDZK1IP1", "ANXA2", "NNMT", "S100A10", "CA9"],
}


def gate(stage: str) -> None:
    result = subprocess.run([str(GATE), stage], text=True, capture_output=True)
    QC.mkdir(parents=True, exist_ok=True)
    with (QC / "memory_gate_log_v1.tsv").open("a", encoding="utf-8") as handle:
        handle.write(stage + "\t" + str(result.returncode) + "\t" + (result.stdout + result.stderr).strip().replace("\t", "|") + "\n")
    if result.returncode != 0:
        raise RuntimeError(f"Memory gate blocked {stage}: {result.stdout}{result.stderr}")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def zscore(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    sd = np.nanstd(values, ddof=1)
    if not np.isfinite(sd) or sd == 0:
        return np.zeros(values.shape, dtype=float)
    return (values - np.nanmean(values)) / sd


def module_score(counts: sparse.csr_matrix, library_size: np.ndarray, indices: list[int]) -> np.ndarray:
    x = counts[:, indices].astype(np.float64).tocsc(copy=True)
    x = sparse.diags(1e4 / np.maximum(library_size, 1.0)) @ x
    x.data = np.log1p(x.data)
    dense = x.toarray()
    mean = dense.mean(axis=0)
    sd = dense.std(axis=0, ddof=1)
    sd[~np.isfinite(sd) | (sd == 0)] = 1.0
    return ((dense - mean) / sd).mean(axis=1)


def main() -> None:
    if (OUT / "RUN_COMPLETE").exists() or (INPUT / "INPUT_COMPLETE").exists():
        raise SystemExit("Refusing to overwrite completed trajectory input/result")
    for path in (H5AD, META, HVG, USAGE, CLONES, AUC, SCFEA):
        if not path.exists() or path.stat().st_size == 0:
            raise FileNotFoundError(path)
    INPUT.mkdir(parents=True, exist_ok=True)
    (OUT / "figures").mkdir(parents=True, exist_ok=True)
    (OUT / "tables").mkdir(parents=True, exist_ok=True)
    (OUT / "source_data").mkdir(parents=True, exist_ok=True)
    (OUT / "objects").mkdir(parents=True, exist_ok=True)
    (OUT / "logs").mkdir(parents=True, exist_ok=True)
    QC.mkdir(parents=True, exist_ok=True)

    gate("part1_trajectory_input_read_v1")
    meta = pd.read_csv(META, sep="\t", dtype={"exact_cell_id": str}).set_index("exact_cell_id", drop=False)
    if len(meta) != 11950 or not meta.index.is_unique or not meta["eligible_for_subclone_selection_v1"].astype(bool).all():
        raise ValueError("Strict 11,950-cell trajectory registry invariant failed")
    cells = meta.index.tolist()
    usage = pd.read_csv(USAGE, sep="\t", dtype={"exact_cell_id": str}).set_index("exact_cell_id")
    clones = pd.read_csv(CLONES, sep="\t", dtype={"exact_cell_id": str}).set_index("exact_cell_id")
    auc = pd.read_csv(AUC, sep="\t", dtype={"exact_cell_id": str}).set_index("exact_cell_id")
    scfea = pd.read_csv(SCFEA, sep="\t", dtype={"exact_cell_id": str}).set_index("exact_cell_id")
    for name, frame in (("usage", usage), ("clones", clones), ("auc", auc), ("scfea", scfea)):
        if set(frame.index) != set(cells):
            raise ValueError(f"{name} cell set mismatch")

    adata = ad.read_h5ad(H5AD)
    if not set(cells).issubset(set(adata.obs_names)):
        raise ValueError("Trajectory cells absent from CNV H5AD")
    adata = adata[cells].copy()
    if "counts" not in adata.layers or "X_scVI" not in adata.obsm or "X_umap" not in adata.obsm:
        raise ValueError("Required counts/X_scVI/X_umap missing")
    counts_all = sparse.csr_matrix(adata.layers["counts"])
    if counts_all.nnz == 0 or np.any(counts_all.data < 0) or not np.allclose(counts_all.data, np.rint(counts_all.data)):
        raise ValueError("Selected raw count matrix is not non-negative integer-valued")

    var_index = pd.Index(adata.var_names.astype(str))
    library_size = np.asarray(counts_all.sum(axis=1)).ravel()
    coverage_rows = []
    for state, genes in STATE_GENES.items():
        present = [gene for gene in genes if gene in var_index]
        if len(present) < 8:
            raise ValueError(f"Insufficient genes for {state}: {present}")
        indices = [int(var_index.get_loc(gene)) for gene in present]
        meta[f"state_{state}_score_v1"] = module_score(counts_all, library_size, indices)
        coverage_rows.extend({"state": state, "gene": gene, "present": gene in present} for gene in genes)

    anchor_parts = pd.DataFrame(index=meta.index)
    anchor_parts["CNV"] = pd.to_numeric(meta["infercnv_burden_over_threshold_v1"], errors="coerce").fillna(0.0)
    anchor_parts["Hypoxia"] = meta["state_Hypoxia_score_v1"]
    anchor_parts["Cycle"] = meta["state_Cycle_score_v1"]
    anchor_parts["Base_negative"] = -meta["state_Base_score_v1"]
    patient = meta["curated_patient_id"].astype(str)
    standardized = anchor_parts.groupby(patient, sort=False).transform(lambda x: zscore(x.to_numpy()))
    meta["progression_anchor_v1"] = standardized.mean(axis=1)
    root = pd.Series(False, index=meta.index)
    for _, idx in meta.groupby("curated_patient_id", sort=False).groups.items():
        n = max(5, int(np.ceil(len(idx) * 0.075)))
        chosen = meta.loc[idx, "progression_anchor_v1"].sort_values().index[:n]
        root.loc[chosen] = True
    meta["root_candidate_v1"] = root
    state_cols = [f"state_{x}_score_v1" for x in STATE_GENES]
    meta["malignant_state_v1"] = pd.DataFrame(meta[state_cols].to_numpy(), index=meta.index, columns=list(STATE_GENES)).idxmax(axis=1)

    keep_meta = [
        "exact_cell_id", "study", "sample_id", "curated_patient_id", "infercnv_burden_over_threshold_v1",
        "TYMP_raw_count_v1", "TYMP_log1p_CP10K_v1", "progression_anchor_v1", "root_candidate_v1",
        "malignant_state_v1", *state_cols,
    ]
    merged = meta[keep_meta].copy()
    merged = merged.join(usage.loc[cells, [f"F{x}" for x in range(1, 8)]])
    merged = merged.join(clones.loc[cells, ["formal_inferred_cnv_clone_v1", "clone_size_class_v1", "threshold_sensitive_patient_v1"]])
    merged = merged.join(auc.loc[cells, [c for c in auc.columns if c.startswith("AUC_")]])
    merged = merged.join(scfea.loc[cells, ["Pyrimidine_flux", "CysMet_flux", "ArgPro_flux", "Choline_proxy"]])
    if merged.isna().all(axis=0).any() or merged.index.tolist() != cells:
        raise ValueError("Merged trajectory metadata invariant failed")

    hvg = pd.read_csv(HVG, sep="\t")["gene"].astype(str).tolist()
    hvg = [gene for gene in hvg if gene in var_index]
    if len(hvg) != 2500:
        raise ValueError(f"Expected 2,500 trajectory HVGs, found {len(hvg)}")
    hvg_idx = var_index.get_indexer(hvg)
    hvg_counts = counts_all[:, hvg_idx].tocsr()
    cyto_detected = np.asarray((counts_all > 0).sum(axis=0)).ravel()
    cyto_idx = np.flatnonzero(cyto_detected >= 10)
    cyto_genes = var_index[cyto_idx].astype(str).tolist()
    cyto_counts = counts_all[:, cyto_idx].tocsr()
    if len(cyto_genes) < 8000:
        raise ValueError(f"Insufficient CytoTRACE2 gene coverage: {len(cyto_genes)}")

    gate("part1_trajectory_input_write_v1")
    mtx_path = INPUT / "raw_counts_genes_by_cells_v1.mtx"
    mmwrite(mtx_path, hvg_counts.T, field="integer")
    with mtx_path.open("rb") as source, gzip.open(str(mtx_path) + ".gz", "wb", compresslevel=6) as target:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            target.write(block)
    mtx_path.unlink()
    cyto_mtx_path = INPUT / "cytotrace_counts_genes_by_cells_v1.mtx"
    mmwrite(cyto_mtx_path, cyto_counts.T, field="integer")
    with cyto_mtx_path.open("rb") as source, gzip.open(str(cyto_mtx_path) + ".gz", "wb", compresslevel=6) as target:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            target.write(block)
    cyto_mtx_path.unlink()
    pd.Series(cells).to_csv(INPUT / "cells_v1.tsv", sep="\t", index=False, header=False)
    pd.Series(hvg).to_csv(INPUT / "genes_v1.tsv", sep="\t", index=False, header=False)
    pd.Series(cyto_genes).to_csv(INPUT / "cytotrace_genes_v1.tsv", sep="\t", index=False, header=False)
    merged.to_csv(INPUT / "cell_metadata_v1.tsv.gz", sep="\t", index=False, compression="gzip")
    latent_frame = pd.DataFrame(adata.obsm["X_scVI"], index=cells, columns=[f"scVI_{i+1}" for i in range(adata.obsm["X_scVI"].shape[1])])
    latent_frame.index.name = "exact_cell_id"
    latent_frame.to_csv(INPUT / "scvi_latent_v1.tsv.gz", sep="\t", compression="gzip")
    umap_frame = pd.DataFrame(adata.obsm["X_umap"], index=cells, columns=["UMAP1", "UMAP2"])
    umap_frame.index.name = "exact_cell_id"
    umap_frame.to_csv(INPUT / "scvi_umap_v1.tsv.gz", sep="\t", compression="gzip")
    pd.DataFrame(coverage_rows).to_csv(QC / "state_gene_coverage_v1.tsv", sep="\t", index=False)
    input_files = [INPUT / "raw_counts_genes_by_cells_v1.mtx.gz", INPUT / "cytotrace_counts_genes_by_cells_v1.mtx.gz", INPUT / "cells_v1.tsv", INPUT / "genes_v1.tsv", INPUT / "cytotrace_genes_v1.tsv", INPUT / "cell_metadata_v1.tsv.gz", INPUT / "scvi_latent_v1.tsv.gz", INPUT / "scvi_umap_v1.tsv.gz"]
    manifest = pd.DataFrame([{"file": str(x.relative_to(OUT)), "size_bytes": x.stat().st_size, "sha256": sha256(x)} for x in input_files])
    manifest.to_csv(INPUT / "file_manifest_v1.tsv", sep="\t", index=False)
    metadata = {
        "schema_version": "part1_trajectory_input_v1", "completed_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "n_cells": len(cells), "n_genes": len(hvg), "n_cytotrace_genes": len(cyto_genes), "n_patients": int(merged["curated_patient_id"].nunique()),
        "n_studies": int(merged["study"].nunique()), "n_root_candidates": int(root.sum()),
        "root_anchor": "within-patient z(CNV burden + hypoxia + cycle - base)/4; excludes TYMP,F1,metabolism",
        "embedding": "frozen Part1 scVI latent and UMAP", "random_seed": SEED,
    }
    (INPUT / "run_metadata_v1.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    (INPUT / "INPUT_COMPLETE").write_text(datetime.now().astimezone().isoformat(timespec="seconds") + "\n", encoding="utf-8")
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
