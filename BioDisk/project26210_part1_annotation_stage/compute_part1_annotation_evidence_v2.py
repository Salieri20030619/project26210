#!/usr/bin/env python3
"""Compute reproducible, marker-based evidence for Part 1 broad cell annotation."""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from datetime import datetime
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import sparse


PROJECT = Path("/mnt/BioDisk/ludingyu/project26210")
INPUT = PROJECT / "results/part1_rebuild_v1/02_integration_annotation/scvi_integration_v1/part1_scvi_integrated_v1.h5ad"
MEMORY_GATE = PROJECT / "check_memory_gate.sh"
DEFAULT_OUT = PROJECT / "results/part1_rebuild_v1/02_integration_annotation/annotation_evidence_v2"

# Positive markers only. Generic HLA, ribosomal, mitochondrial and cell-cycle
# genes are intentionally excluded because they do not establish lineage.
MARKERS: dict[str, list[str]] = {
    "T cell": ["CD3D", "CD3E", "CD3G", "TRAC", "TRBC1", "TRBC2", "LCK", "IL7R", "CD247"],
    "NK cell": ["NKG7", "GNLY", "KLRD1", "PRF1", "CTSW", "FCGR3A", "XCL1", "XCL2"],
    "B cell": ["MS4A1", "CD79A", "CD79B", "CD37", "CD22", "CD19", "CD74", "CD83"],
    "Plasma cell": ["MZB1", "JCHAIN", "SDC1", "DERL3", "TNFRSF17", "IGHG1", "IGHG3", "IGKC"],
    "Erythrocyte": ["HBB", "HBA1", "HBA2", "HBG1", "HBG2", "GYPA", "GYPB", "ALAS2", "AHSP", "SLC4A1"],
    "Macrophage": ["C1QA", "C1QB", "C1QC", "APOC1", "CD68", "CD163", "MSR1", "MRC1", "FCGR1A"],
    "Monocyte": ["S100A8", "S100A9", "FCN1", "VCAN", "CTSS", "LYZ", "LILRB1", "LGALS3"],
    "Dendritic cell": ["CLEC9A", "XCR1", "WDFY4", "BATF3", "CADM1", "FLT3", "DNASE1L3", "CD1C", "CD1E", "CLEC10A", "FCER1A"],
    "Mast cell": ["TPSAB1", "TPSB2", "CPA3", "KIT", "MS4A2", "HDC", "HPGDS"],
    "Epithelial cell": ["EPCAM", "KRT8", "KRT18", "KRT19", "KRT7", "MUC1", "KRT14", "KRT5", "PAX8"],
    "Endothelial": ["PECAM1", "VWF", "EMCN", "KDR", "ESAM", "ENG", "CDH5", "RAMP2", "PLVAP"],
    "Fibroblast": ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "COL6A1", "COL6A2", "PDGFRA", "C7", "CFD"],
    "Pericyte/SMC": ["RGS5", "PDGFRB", "CSPG4", "MCAM", "ACTA2", "MYH11", "ABCC9", "KCNJ8", "DES", "NOTCH3", "TAGLN", "CNN1"],
}

DISPLAY_MARKERS = {
    key: value[:5] for key, value in MARKERS.items()
}


def memory_gate(stage: str, records: list[dict[str, object]]) -> None:
    proc = subprocess.run([str(MEMORY_GATE)], text=True, capture_output=True)
    output = (proc.stdout or proc.stderr).strip()
    print(f"memory_gate[{stage}]: {output}", flush=True)
    records.append({"stage": stage, "returncode": proc.returncode, "output": output})
    if proc.returncode != 0:
        raise RuntimeError(f"Memory gate failed before {stage}: {output}")


def zscore_columns(values: np.ndarray) -> np.ndarray:
    means = values.mean(axis=0, keepdims=True)
    std = values.std(axis=0, keepdims=True)
    std[std < 1e-8] = 1.0
    return (values - means) / std


def save_both(fig: plt.Figure, stem: Path) -> None:
    fig.savefig(stem.with_suffix(".png"), dpi=220, bbox_inches="tight")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--chunk-size", type=int, default=4096)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    outdir: Path = args.outdir
    complete = outdir / "RUN_COMPLETE"
    if complete.exists():
        raise SystemExit(f"Refusing to overwrite completed evidence run: {complete}")
    if outdir.exists() and any(outdir.iterdir()) and not args.force:
        raise SystemExit(f"Non-empty incomplete output exists; inspect it or use --force: {outdir}")
    outdir.mkdir(parents=True, exist_ok=True)

    started = datetime.now().astimezone().isoformat()
    memory_records: list[dict[str, object]] = []
    memory_gate("before_h5ad_open", memory_records)
    data = ad.read_h5ad(INPUT, backed="r")
    if "leiden04" not in data.obs:
        raise RuntimeError("Required cluster field leiden04 is absent")
    if "X_umap" not in data.obsm:
        raise RuntimeError("Required embedding X_umap is absent")
    if not (sparse.issparse(data.X) or type(data.X).__name__ in {"_CSRDataset", "CSRDataset"}):
        raise RuntimeError(f"Expected sparse log-normalized X, got {type(data.X).__name__}")

    var_names = pd.Index(data.var_names.astype(str))
    cluster_text = data.obs["leiden04"].astype("string")
    cluster_levels = sorted(cluster_text.dropna().unique(), key=lambda x: int(x))
    cluster_to_code = {cluster: i for i, cluster in enumerate(cluster_levels)}
    codes = cluster_text.map(cluster_to_code).to_numpy(dtype=np.int32)
    if (codes < 0).any():
        raise RuntimeError("Missing Leiden04 cluster assignments")
    n_clusters = len(cluster_levels)
    n_genes = data.n_vars
    counts = np.bincount(codes, minlength=n_clusters).astype(np.int64)

    marker_rows = []
    for lineage, genes in MARKERS.items():
        for order, gene in enumerate(genes, start=1):
            marker_rows.append(
                {
                    "cell_type": lineage,
                    "gene": gene,
                    "marker_order": order,
                    "present": gene in var_names,
                    "gene_index": int(var_names.get_loc(gene)) if gene in var_names else -1,
                }
            )
    marker_panel = pd.DataFrame(marker_rows)
    marker_panel.to_csv(outdir / "annotation_marker_panel_v2.tsv", sep="\t", index=False)
    if marker_panel.groupby("cell_type")["present"].sum().min() < 3:
        raise RuntimeError("At least one cell type has fewer than three marker genes in the matrix")

    memory_gate("before_expression_scan", memory_records)
    sums = np.zeros((n_clusters, n_genes), dtype=np.float64)
    detected = np.zeros((n_clusters, n_genes), dtype=np.int64)
    x_min = np.inf
    x_max = -np.inf
    nnz = 0
    for start in range(0, data.n_obs, args.chunk_size):
        end = min(start + args.chunk_size, data.n_obs)
        x = data.X[start:end]
        if not sparse.isspmatrix_csr(x):
            x = sparse.csr_matrix(x)
        if x.data.size:
            x_min = min(x_min, float(x.data.min()))
            x_max = max(x_max, float(x.data.max()))
        nnz += int(x.nnz)
        rows = np.arange(end - start, dtype=np.int32)
        membership = sparse.csr_matrix(
            (np.ones(end - start, dtype=np.float64), (rows, codes[start:end])),
            shape=(end - start, n_clusters),
        )
        sums += (membership.T @ x).toarray()
        binary = x.copy()
        binary.data = np.ones_like(binary.data, dtype=np.int8)
        detected += (membership.T @ binary).toarray().astype(np.int64)
        if end == data.n_obs or end % 49152 < args.chunk_size:
            print(f"expression_scan: {end}/{data.n_obs} cells", flush=True)

    means = sums / counts[:, None]
    fractions = detected / counts[:, None]
    global_sums = sums.sum(axis=0)
    global_counts = counts.sum()
    other_means = (global_sums[None, :] - sums) / (global_counts - counts)[:, None]
    delta_means = means - other_means
    global_detected = detected.sum(axis=0)
    other_fractions = (global_detected[None, :] - detected) / (global_counts - counts)[:, None]
    delta_fractions = fractions - other_fractions

    marker_indices = marker_panel.loc[marker_panel["present"], "gene_index"].to_numpy(dtype=int)
    marker_gene_names = marker_panel.loc[marker_panel["present"], "gene"].tolist()
    marker_mean = means[:, marker_indices]
    marker_fraction = fractions[:, marker_indices]
    marker_zmean = zscore_columns(marker_mean)
    marker_zfrac = zscore_columns(marker_fraction)
    marker_combined = marker_zmean + 0.5 * marker_zfrac

    marker_summary_rows = []
    marker_lookup = marker_panel.loc[marker_panel["present"]].reset_index(drop=True)
    for ci, cluster in enumerate(cluster_levels):
        for gi, row in marker_lookup.iterrows():
            marker_summary_rows.append(
                {
                    "leiden04": cluster,
                    "n_cells": int(counts[ci]),
                    "cell_type": row["cell_type"],
                    "gene": row["gene"],
                    "mean_log1p_expression": float(marker_mean[ci, gi]),
                    "fraction_expressing": float(marker_fraction[ci, gi]),
                    "z_mean_across_clusters": float(marker_zmean[ci, gi]),
                    "z_fraction_across_clusters": float(marker_zfrac[ci, gi]),
                    "combined_gene_evidence": float(marker_combined[ci, gi]),
                }
            )
    marker_summary = pd.DataFrame(marker_summary_rows)
    marker_summary.to_csv(outdir / "cluster_marker_gene_evidence_v2.tsv.gz", sep="\t", index=False, compression="gzip")

    type_scores = np.zeros((n_clusters, len(MARKERS)), dtype=float)
    type_mean = np.zeros_like(type_scores)
    type_fraction = np.zeros_like(type_scores)
    score_rows = []
    types = list(MARKERS)
    for ti, lineage in enumerate(types):
        mask = marker_lookup["cell_type"].eq(lineage).to_numpy()
        type_scores[:, ti] = marker_combined[:, mask].mean(axis=1)
        type_mean[:, ti] = marker_mean[:, mask].mean(axis=1)
        type_fraction[:, ti] = marker_fraction[:, mask].mean(axis=1)
    for ci, cluster in enumerate(cluster_levels):
        order = np.argsort(type_scores[ci])[::-1]
        for rank, ti in enumerate(order, start=1):
            score_rows.append(
                {
                    "leiden04": cluster,
                    "n_cells": int(counts[ci]),
                    "rank": rank,
                    "cell_type": types[ti],
                    "marker_score": float(type_scores[ci, ti]),
                    "mean_marker_log1p_expression": float(type_mean[ci, ti]),
                    "mean_marker_fraction_expressing": float(type_fraction[ci, ti]),
                    "score_margin_from_next": float(type_scores[ci, order[0]] - type_scores[ci, order[1]]) if rank == 1 else np.nan,
                }
            )
    pd.DataFrame(score_rows).to_csv(outdir / "cluster_celltype_scores_v2.tsv", sep="\t", index=False)

    excluded = np.array(
        [
            gene.startswith(("MT-", "RPS", "RPL"))
            or gene in {
                "MALAT1", "XIST", "NEAT1", "FOS", "JUN", "JUNB", "FOSB",
                "HSP90AA1", "HSPA1A", "HSPA1B", "COX1", "COX2", "COX3",
                "ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6",
                "ATP6", "ATP8", "CYTB",
            }
            for gene in var_names
        ]
    )
    top_rows = []
    for ci, cluster in enumerate(cluster_levels):
        score = delta_means[ci] + 0.25 * np.maximum(delta_fractions[ci], 0)
        score[excluded] = -np.inf
        candidates = np.argpartition(score, -100)[-100:]
        candidates = candidates[np.argsort(score[candidates])[::-1]][:50]
        for rank, gi in enumerate(candidates, start=1):
            top_rows.append(
                {
                    "leiden04": cluster,
                    "n_cells": int(counts[ci]),
                    "rank": rank,
                    "gene": var_names[gi],
                    "specificity_score": float(score[gi]),
                    "mean_log1p_expression": float(means[ci, gi]),
                    "other_mean_log1p_expression": float(other_means[ci, gi]),
                    "delta_mean": float(delta_means[ci, gi]),
                    "fraction_expressing": float(fractions[ci, gi]),
                    "other_fraction_expressing": float(other_fractions[ci, gi]),
                    "delta_fraction": float(delta_fractions[ci, gi]),
                }
            )
    pd.DataFrame(top_rows).to_csv(outdir / "cluster_top_specific_genes_v2.tsv", sep="\t", index=False)

    cross = pd.crosstab(data.obs["leiden03"].astype("string"), cluster_text, dropna=False)
    cross.index.name = "leiden03"
    cross.to_csv(outdir / "leiden03_by_leiden04_counts_v2.tsv", sep="\t")
    composition_rows = []
    for field in ["curated_patient_id", "sample_id", "source_tissue", "study"]:
        if field not in data.obs:
            continue
        tab = pd.crosstab(cluster_text, data.obs[field].astype("string"), dropna=False)
        for cluster in tab.index:
            total = int(tab.loc[cluster].sum())
            for value, n in tab.loc[cluster][tab.loc[cluster] > 0].items():
                composition_rows.append(
                    {
                        "leiden04": str(cluster),
                        "field": field,
                        "value": str(value),
                        "n_cells": int(n),
                        "fraction_within_cluster": float(n / total),
                    }
                )
    pd.DataFrame(composition_rows).to_csv(outdir / "cluster_sample_composition_v2.tsv.gz", sep="\t", index=False, compression="gzip")

    sns.set_theme(style="white", font_scale=0.8)
    score_frame = pd.DataFrame(type_scores, index=cluster_levels, columns=types)
    fig, ax = plt.subplots(figsize=(12, 10))
    sns.heatmap(score_frame, cmap="vlag", center=0, annot=True, fmt=".1f", linewidths=0.25, ax=ax)
    ax.set_xlabel("Candidate broad cell type")
    ax.set_ylabel("Leiden04 cluster")
    ax.set_title("Marker evidence scores by Leiden04 cluster")
    save_both(fig, outdir / "cluster_celltype_score_heatmap_v2")

    display_genes = [gene for lineage in types for gene in DISPLAY_MARKERS[lineage] if gene in var_names]
    display_idx = [var_names.get_loc(gene) for gene in display_genes]
    display_z = zscore_columns(means[:, display_idx])
    fig, ax = plt.subplots(figsize=(24, 10))
    for ci in range(n_clusters):
        sizes = 10 + 90 * fractions[ci, display_idx]
        ax.scatter(np.arange(len(display_genes)), np.full(len(display_genes), ci), s=sizes, c=display_z[ci], cmap="vlag", vmin=-2.5, vmax=2.5, edgecolors="none")
    ax.set_xticks(np.arange(len(display_genes)), display_genes, rotation=90)
    ax.set_yticks(np.arange(n_clusters), cluster_levels)
    ax.set_xlabel("Canonical marker (dot size = expressing fraction; color = cluster-standardized mean)")
    ax.set_ylabel("Leiden04 cluster")
    ax.set_title("Canonical marker dot plot for annotation review")
    ax.invert_yaxis()
    save_both(fig, outdir / "cluster_marker_dotplot_v2")

    umap = np.asarray(data.obsm["X_umap"])
    fig, ax = plt.subplots(figsize=(10, 8))
    palette = sns.color_palette("husl", n_clusters)
    for ci, cluster in enumerate(cluster_levels):
        mask = codes == ci
        ax.scatter(umap[mask, 0], umap[mask, 1], s=0.35, color=palette[ci], alpha=0.65, rasterized=True, label=cluster)
        center = np.median(umap[mask], axis=0)
        ax.text(center[0], center[1], cluster, fontsize=8, weight="bold", ha="center", va="center")
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.set_title("Part 1 rebuilt scVI integration: Leiden04 clusters")
    ax.legend(title="Leiden04", bbox_to_anchor=(1.01, 1), loc="upper left", markerscale=7, ncol=2)
    save_both(fig, outdir / "umap_leiden04_v2")

    metadata = {
        "run_id": "part1_annotation_evidence_v2",
        "started_at": started,
        "completed_at": datetime.now().astimezone().isoformat(),
        "python_executable": sys.executable,
        "command": " ".join(sys.argv),
        "input_h5ad": str(INPUT),
        "input_shape": [int(data.n_obs), int(data.n_vars)],
        "cluster_field": "leiden04",
        "n_clusters": n_clusters,
        "chunk_size": args.chunk_size,
        "x_semantics": "log1p_normalized_expression_from_current_integrated_h5ad_X",
        "x_nonzero_min": float(x_min),
        "x_nonzero_max": float(x_max),
        "x_nnz": int(nnz),
        "annotation_strategy": "canonical_positive_marker_evidence_then_manual_cluster_mapping; legacy_labels_not_used_in_scores",
        "memory_gate": memory_records,
    }
    (outdir / "run_metadata_v2.json").write_text(json.dumps(metadata, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    (outdir / "RUN_COMPLETE").write_text(metadata["completed_at"] + "\n", encoding="utf-8")
    data.file.close()
    print(json.dumps(metadata, indent=2, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
