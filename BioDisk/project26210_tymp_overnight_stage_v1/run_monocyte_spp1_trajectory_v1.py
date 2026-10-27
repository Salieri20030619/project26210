#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import subprocess
from datetime import datetime
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from scipy.stats import binomtest, spearmanr


ROOT = Path("/mnt/BioDisk/ludingyu/project26210")
INPUT = ROOT / "results/part4_spatial/ref/ref_signatures.h5ad"
MEMORY_GATE = ROOT / "check_memory_gate.sh"
OUT = Path("/mnt/BioDisk/wangmingyuan/project26210_tymp_overnight_stage_v1/trajectory_v1")
CELL_TYPES = ["Monocyte", "Other Macrophage", "APOE+C1QC+ TAM", "SPP1+ TAM"]
SIGNATURES = {
    "Monocyte_program": ["FCN1", "VCAN", "S100A8", "S100A9", "CTSS", "LILRB1"],
    "SPP1_state_no_SPP1": ["GPNMB", "LGALS3", "FABP5", "CD9", "LPL"],
    "C1QC_state": ["C1QA", "C1QB", "C1QC"],
    "Vascular_remodeling": ["VEGFA", "MMP9", "MMP14", "ANGPTL4", "IL1B", "LGALS3"],
    "Immunosuppression": ["IL10", "TGFB1", "VSIG4", "MARCO", "CD274", "IDO1"],
    "MTA_context": ["MAT2A", "MAT2B", "AMD1", "SRM", "SMS"],
}


def now() -> str:
    return datetime.now().astimezone().isoformat(timespec="seconds")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def memory_gate(stage: str) -> str:
    result = subprocess.run([str(MEMORY_GATE)], capture_output=True, text=True, check=False)
    line = (result.stdout or result.stderr).strip()
    print(f"[memory_gate] {stage} {line}", flush=True)
    if result.returncode != 0:
        raise RuntimeError(line)
    return line


def zscore(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    sd = np.nanstd(values)
    return np.zeros(len(values)) if sd <= 0 or not np.isfinite(sd) else (values - np.nanmean(values)) / sd


def bh(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    order = np.argsort(values)
    ranked = values[order] * len(values) / np.arange(1, len(values) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty_like(ranked)
    out[order] = np.minimum(ranked, 1.0)
    return out


def main() -> None:
    if OUT.exists():
        raise SystemExit(f"Refusing to overwrite: {OUT}")
    tables = OUT / "tables"
    figures = OUT / "figures"
    tables.mkdir(parents=True)
    figures.mkdir(parents=True)
    started = now()
    gates = [{"time": started, "stage": "initial", "line": memory_gate("initial")}]

    backed = ad.read_h5ad(INPUT, backed="r")
    if "counts" not in backed.layers or "X_scVI" not in backed.obsm:
        raise ValueError("Reference lacks counts or X_scVI")
    mask = backed.obs["c2l_celltype"].astype(str).isin(CELL_TYPES).to_numpy()
    if mask.sum() < 1000:
        raise ValueError(f"Too few selected cells: {mask.sum()}")
    gates.append({"time": now(), "stage": "load_myeloid_subset", "line": memory_gate("load_myeloid_subset")})
    adata = backed[mask].to_memory()
    backed.file.close()
    adata.obs["trajectory_celltype"] = pd.Categorical(adata.obs["c2l_celltype"].astype(str), categories=CELL_TYPES)

    gates.append({"time": now(), "stage": "neighbors_diffusion", "line": memory_gate("neighbors_diffusion")})
    sc.pp.neighbors(adata, n_neighbors=30, n_pcs=30, use_rep="X_scVI", metric="euclidean", random_state=20260803)
    sc.tl.diffmap(adata, n_comps=20, random_state=20260803)
    monocyte = np.flatnonzero(adata.obs["trajectory_celltype"].astype(str).eq("Monocyte").to_numpy())
    latent = np.asarray(adata.obsm["X_scVI"], dtype=float)
    center = np.median(latent[monocyte], axis=0)
    root = int(monocyte[np.argmin(((latent[monocyte] - center) ** 2).sum(axis=1))])
    adata.uns["iroot"] = root
    sc.tl.dpt(adata, n_dcs=15)

    requested = sorted({gene for genes in SIGNATURES.values() for gene in genes})
    present = [gene for gene in requested if gene in adata.var_names]
    matrix = adata[:, present].layers["counts"]
    if sparse.issparse(matrix):
        matrix = matrix.toarray()
    matrix = np.asarray(matrix, dtype=float)
    totals = np.asarray(adata.layers["counts"].sum(axis=1)).reshape(-1).astype(float)
    expression = np.log1p(matrix / np.maximum(totals[:, None], 1) * 1e4)
    expression = pd.DataFrame(expression, index=adata.obs_names, columns=present)
    gene_qc = []
    for name, genes in SIGNATURES.items():
        used = [gene for gene in genes if gene in present]
        if len(used) < 3:
            raise ValueError(f"Too few genes for {name}: {used}")
        adata.obs[name] = np.column_stack([zscore(expression[gene].to_numpy()) for gene in used]).mean(axis=1)
        gene_qc.append({"signature": name, "requested": ";".join(genes), "used": ";".join(used), "n_used": len(used)})
    adata.obs["SPP1_vs_C1QC"] = zscore(adata.obs["SPP1_state_no_SPP1"]) - zscore(adata.obs["C1QC_state"])

    fields = ["sample_id", "study", "trajectory_celltype", "dpt_pseudotime", *SIGNATURES, "SPP1_vs_C1QC"]
    cell_table = adata.obs[fields].copy()
    cell_table.insert(0, "cell_id", adata.obs_names)
    cell_table.to_csv(tables / "myeloid_trajectory_cell_scores_v1.tsv.gz", sep="\t", index=False, compression="gzip")
    pd.DataFrame(gene_qc).to_csv(tables / "signature_gene_qc_v1.tsv", sep="\t", index=False)

    summary = cell_table.groupby("trajectory_celltype", observed=True).agg(
        n_cells=("cell_id", "size"), n_samples=("sample_id", "nunique"),
        median_pseudotime=("dpt_pseudotime", "median"), median_spp1_state=("SPP1_state_no_SPP1", "median"),
        median_c1qc_state=("C1QC_state", "median"), median_vascular_remodeling=("Vascular_remodeling", "median"),
        median_immunosuppression=("Immunosuppression", "median"),
    ).reset_index()
    summary.to_csv(tables / "celltype_trajectory_summary_v1.tsv", sep="\t", index=False)

    patient_rows = []
    score_names = [*SIGNATURES, "SPP1_vs_C1QC"]
    for sample, frame in cell_table.groupby("sample_id", sort=True):
        if len(frame) < 50 or frame["trajectory_celltype"].nunique() < 2:
            continue
        for score in score_names:
            result = spearmanr(frame["dpt_pseudotime"], frame[score], nan_policy="omit")
            patient_rows.append({"sample_id": sample, "study": frame["study"].mode().iat[0], "score": score,
                                 "n_cells": len(frame), "rho": float(result.statistic), "p_value": float(result.pvalue)})
    patient = pd.DataFrame(patient_rows)
    patient["p_fdr_bh_within_sample"] = patient.groupby("sample_id")["p_value"].transform(lambda x: bh(x.to_numpy()))
    patient.to_csv(tables / "sample_pseudotime_associations_v1.tsv", sep="\t", index=False)

    cohort_rows = []
    for score, frame in patient.groupby("score", sort=True):
        positive = int((frame["rho"] > 0).sum())
        cohort_rows.append({"score": score, "n_samples": len(frame), "positive_samples": positive,
                            "median_rho": float(frame["rho"].median()),
                            "sign_test_p": float(binomtest(positive, len(frame), 0.5).pvalue)})
    cohort = pd.DataFrame(cohort_rows)
    cohort["sign_test_fdr_bh"] = bh(cohort["sign_test_p"].to_numpy())
    cohort.to_csv(tables / "cohort_pseudotime_direction_summary_v1.tsv", sep="\t", index=False)

    connectivities = adata.obsp["connectivities"].tocsr()
    labels = adata.obs["trajectory_celltype"].astype(str).to_numpy()
    transition_rows = []
    for source in CELL_TYPES:
        source_mask = labels == source
        denominator = float(connectivities[source_mask].sum())
        for target in CELL_TYPES:
            weight = float(connectivities[source_mask][:, labels == target].sum())
            transition_rows.append({"source": source, "target": target, "edge_weight": weight,
                                    "source_normalized_fraction": weight / denominator if denominator else np.nan})
    pd.DataFrame(transition_rows).to_csv(tables / "trajectory_neighbor_transition_matrix_v1.tsv", sep="\t", index=False)
    pd.DataFrame(gates).to_csv(tables / "memory_gates_v1.tsv", sep="\t", index=False)

    sc.pl.umap(adata, color=["trajectory_celltype", "dpt_pseudotime", "SPP1_state_no_SPP1", "SPP1_vs_C1QC", "Vascular_remodeling"],
               ncols=3, show=False, save=False)
    import matplotlib.pyplot as plt
    plt.gcf().set_size_inches(15, 9)
    plt.tight_layout()
    plt.savefig(figures / "monocyte_spp1_trajectory_overview_v1.png", dpi=220, bbox_inches="tight")
    plt.close("all")

    metadata = {
        "status": "completed", "started_at": started, "completed_at": now(), "input": str(INPUT),
        "input_sha256": sha256(INPUT), "n_cells": adata.n_obs, "n_samples_tested": int(patient["sample_id"].nunique()),
        "root_cell": str(adata.obs_names[root]), "root_celltype": str(adata.obs["trajectory_celltype"].iat[root]),
        "method": "scVI-neighbor graph, diffusion map and DPT; sample-level Spearman direction synthesis",
        "inference_boundary": "Pseudotime is cross-sectional state ordering, not RNA velocity or lineage tracing; no spliced/unspliced layer is available.",
    }
    (OUT / "run_metadata_v1.json").write_text(json.dumps(metadata, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    (OUT / "RUN_COMPLETE").write_text(now() + "\n", encoding="utf-8")
    print(json.dumps(metadata, ensure_ascii=False, indent=2), flush=True)


if __name__ == "__main__":
    main()
