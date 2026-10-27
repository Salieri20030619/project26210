#!/usr/bin/env python3
"""Apply reviewed broad cell annotations to the rebuilt Part 1 scVI object."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
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


PROJECT = Path("/mnt/BioDisk/ludingyu/project26210")
INPUT = PROJECT / "results/part1_rebuild_v1/02_integration_annotation/scvi_integration_v1/part1_scvi_integrated_v1.h5ad"
EVIDENCE = PROJECT / "results/part1_rebuild_v1/02_integration_annotation/annotation_evidence_v2"
AUDIT = PROJECT / "results/part1_rebuild_v1/02_integration_annotation/annotation_input_audit_v1"
MEMORY_GATE = PROJECT / "check_memory_gate.sh"
DEFAULT_OUT = PROJECT / "results/part1_rebuild_v1/02_integration_annotation/cell_annotation_v1"


MAPPING = [
    ("0", "T cell", "High", False, "CD3D/CD3E/CD3G/TRAC/IL7R; T-cell score rank 1"),
    ("1", "Macrophage", "High", False, "C1QA/C1QB/C1QC/TYROBP/CD68; macrophage score rank 1"),
    ("2", "Endothelial", "High", False, "VWF/PLVAP/RAMP2/FLT1/PECAM1; endothelial score rank 1"),
    ("3", "NK cell", "High", False, "NKG7/GNLY/KLRD1/PRF1/CTSW; NK-cell score rank 1"),
    ("4", "Epithelial cell", "High", False, "KRT18/KRT8/CD24 plus epithelial score rank 1"),
    ("5", "Pericyte/SMC", "High", False, "RGS5/MYL9/TAGLN/ACTA2/PDGFRB; pericyte/SMC score rank 1"),
    ("6", "T cell", "High", False, "CD3D/CD3E/TRAC/TRBC2 with GZMK; T-cell score rank 1"),
    ("7", "Epithelial cell", "High", False, "Renal epithelial ALDOB/MIOX/FXYD2/FABP1/HPD"),
    ("8", "Macrophage", "High", False, "C1Q/APOC1 macrophage program; macrophage score rank 1"),
    ("9", "Monocyte", "High", False, "S100A8/S100A9/FCN1/VCAN/LYZ; monocyte score rank 1"),
    ("10", "Epithelial cell", "Medium", False, "SPP1/CD24 epithelial-hypoxia program; epithelial score rank 1"),
    ("11", "Epithelial cell", "Medium", False, "KRT18/KRT8/ANXA4/CD24; epithelial score rank 1"),
    ("12", "Endothelial", "High", False, "PECAM1/VWF/EMCN/KDR/ESAM; endothelial score rank 1"),
    ("13", "Epithelial cell", "Medium", False, "FXYD2/CD24 renal epithelial program with elevated mitochondrial signal"),
    ("14", "T cell", "High", False, "CD3D/CD3E/CD3G/TRAC; T-cell score rank 1"),
    ("15", "Epithelial cell", "High", False, "EPCAM/FXYD2/DEFB1/WFDC2; epithelial score rank 1"),
    ("16", "B cell", "High", False, "MS4A1/CD79A/CD37/CD74/BANK1; B-cell score rank 1"),
    ("17", "Mast cell", "High", False, "TPSAB1/TPSB2/CPA3/HPGDS; mast-cell score rank 1"),
    ("18", "T cell", "High", False, "PTPRC/CD2/GZMK/T-cell program; T-cell score rank 1"),
    ("19", "Epithelial cell", "High", False, "SPP1/CD24/FABP7 epithelial program; epithelial score rank 1"),
    ("20", "Plasma cell", "High", False, "JCHAIN/MZB1/DERL3/IGHG/IGKC; plasma-cell score rank 1"),
    ("21", "Epithelial cell", "High", False, "KRT18/CD24/CXCL14 renal epithelial program; epithelial score rank 1"),
    ("22", "Dendritic cell", "High", False, "CLEC9A/XCR1/WDFY4/BATF3/FCER1A; dendritic score rank 1"),
    ("23", "Epithelial/Myeloid mixed", "Low", True, "ALDOB/MIOX/FXYD2 co-occur with TYROBP/CST3; conflicting lineage evidence"),
    ("24", "Erythrocyte", "High", False, "HBB/HBA1/HBA2/GYPA/ALAS2; erythrocyte score rank 1"),
    ("25", "Fibroblast", "Medium", False, "DCN/COL1A1/COL1A2/CFD/CCDC80; fibroblast score rank 1; n=53"),
    ("26", "Fibroblast", "Medium", False, "COL1A1/COL1A2/COL3A1/SPARC; fibroblast score rank 1; n=39"),
]

TYPE_ORDER = [
    "Epithelial cell", "T cell", "Macrophage", "Endothelial", "NK cell",
    "Pericyte/SMC", "Monocyte", "B cell", "Plasma cell", "Mast cell",
    "Dendritic cell", "Fibroblast", "Erythrocyte", "Epithelial/Myeloid mixed",
]

COLORS = {
    "Epithelial cell": "#D55E00",
    "T cell": "#0072B2",
    "Macrophage": "#8C564B",
    "Endothelial": "#009E73",
    "NK cell": "#56B4E9",
    "Pericyte/SMC": "#CC79A7",
    "Monocyte": "#E69F00",
    "B cell": "#332288",
    "Plasma cell": "#AA4499",
    "Mast cell": "#882255",
    "Dendritic cell": "#44AA99",
    "Fibroblast": "#999933",
    "Erythrocyte": "#CC3311",
    "Epithelial/Myeloid mixed": "#777777",
}


def memory_gate(stage: str, records: list[dict[str, object]]) -> None:
    proc = subprocess.run([str(MEMORY_GATE)], text=True, capture_output=True)
    output = (proc.stdout or proc.stderr).strip()
    print(f"memory_gate[{stage}]: {output}", flush=True)
    records.append({"stage": stage, "returncode": proc.returncode, "output": output})
    if proc.returncode != 0:
        raise RuntimeError(f"Memory gate failed before {stage}: {output}")


def sha256_file(path: Path, chunk_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def save_both(fig: plt.Figure, stem: Path) -> None:
    fig.savefig(stem.with_suffix(".png"), dpi=220, bbox_inches="tight")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    outdir: Path = args.outdir
    complete = outdir / "RUN_COMPLETE"
    if complete.exists():
        raise SystemExit(f"Refusing to overwrite completed annotation: {complete}")
    if outdir.exists() and any(outdir.iterdir()):
        raise SystemExit(f"Refusing to use non-empty incomplete output directory: {outdir}")
    outdir.mkdir(parents=True, exist_ok=True)

    required = [
        INPUT,
        EVIDENCE / "RUN_COMPLETE",
        EVIDENCE / "cluster_celltype_scores_v2.tsv",
        EVIDENCE / "annotation_marker_panel_v2.tsv",
        EVIDENCE / "cluster_top_specific_genes_v2.tsv",
        AUDIT / "RUN_COMPLETE",
        AUDIT / "cluster_legacy_concordance_v1.tsv",
        AUDIT / "legacy_labels_by_current_cluster_v1.tsv.gz",
    ]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError(f"Required inputs missing: {missing}")

    started = datetime.now().astimezone().isoformat()
    memory_records: list[dict[str, object]] = []
    memory_gate("before_full_h5ad_load", memory_records)
    data = ad.read_h5ad(INPUT)
    if data.shape != (258120, 37795):
        raise RuntimeError(f"Unexpected input shape: {data.shape}")
    if not data.obs_names.is_unique or not data.var_names.is_unique:
        raise RuntimeError("Input cell or gene identifiers are not unique")

    mapping = pd.DataFrame(MAPPING, columns=["leiden04", "cell_type_broad_v1", "annotation_confidence_v1", "annotation_ambiguous_v1", "annotation_evidence_v1"])
    observed = set(data.obs["leiden04"].astype(str).unique())
    expected = set(mapping["leiden04"])
    if observed != expected:
        raise RuntimeError(f"Mapping mismatch; missing={sorted(observed - expected)}, extra={sorted(expected - observed)}")
    cluster_counts = data.obs["leiden04"].astype(str).value_counts().rename("n_cells")
    mapping = mapping.merge(cluster_counts, left_on="leiden04", right_index=True, how="left", validate="one_to_one")

    type_map = mapping.set_index("leiden04")["cell_type_broad_v1"].to_dict()
    confidence_map = mapping.set_index("leiden04")["annotation_confidence_v1"].to_dict()
    ambiguous_map = mapping.set_index("leiden04")["annotation_ambiguous_v1"].to_dict()
    evidence_map = mapping.set_index("leiden04")["annotation_evidence_v1"].to_dict()
    clusters = data.obs["leiden04"].astype(str)
    labels = clusters.map(type_map)
    confidences = clusters.map(confidence_map)
    if labels.isna().any() or confidences.isna().any():
        raise RuntimeError("Annotation mapping produced missing values")
    data.obs["cell_type_broad_v1"] = pd.Categorical(labels, categories=TYPE_ORDER, ordered=True)
    data.obs["annotation_confidence_v1"] = pd.Categorical(confidences, categories=["High", "Medium", "Low"], ordered=True)
    data.obs["annotation_ambiguous_v1"] = clusters.map(ambiguous_map).astype(bool).to_numpy()
    data.obs["annotation_evidence_v1"] = pd.Categorical(clusters.map(evidence_map))
    data.obs["annotation_cluster_v1"] = pd.Categorical(clusters, categories=sorted(expected, key=int), ordered=True)

    mapping.to_csv(outdir / "cluster_annotation_mapping_v1.tsv", sep="\t", index=False)
    cell_table = data.obs[["sample_id", "curated_patient_id", "source_tissue", "leiden04", "cell_type_broad_v1", "annotation_confidence_v1", "annotation_ambiguous_v1"]].copy()
    cell_table.index.name = "cell_id"
    cell_table.to_csv(outdir / "cell_annotation_v1.tsv.gz", sep="\t", compression="gzip")

    count_table = data.obs["cell_type_broad_v1"].value_counts(sort=False).rename("n_cells").reset_index()
    count_table["fraction"] = count_table["n_cells"] / data.n_obs
    count_table.to_csv(outdir / "celltype_counts_v1.tsv", sep="\t", index=False)
    for field in ["source_tissue", "study", "curated_patient_id", "sample_id"]:
        table = pd.crosstab(data.obs[field], data.obs["cell_type_broad_v1"], dropna=False)
        table.to_csv(outdir / f"celltype_by_{field}_counts_v1.tsv", sep="\t")

    legacy_concordance = pd.read_csv(AUDIT / "cluster_legacy_concordance_v1.tsv", sep="\t", dtype={"cluster": str})
    legacy_concordance = legacy_concordance.query("cluster_field == 'leiden04' and legacy_field == 'cell_type'")
    review = mapping.merge(legacy_concordance, left_on="leiden04", right_on="cluster", how="left", validate="one_to_one")
    review["new_matches_legacy_top"] = review["cell_type_broad_v1"].eq(review["top_legacy_value"])
    review.to_csv(outdir / "cluster_annotation_vs_legacy_v1.tsv", sep="\t", index=False)
    legacy_long = pd.read_csv(AUDIT / "legacy_labels_by_current_cluster_v1.tsv.gz", sep="\t", dtype={"cluster": str})
    legacy_long = legacy_long.query("cluster_field == 'leiden04' and legacy_field == 'cell_type'")
    crosswalk = legacy_long.merge(mapping[["leiden04", "cell_type_broad_v1"]], left_on="cluster", right_on="leiden04", how="left")
    crosswalk = crosswalk.groupby(["cell_type_broad_v1", "legacy_value"], observed=True, as_index=False)["n_cells"].sum()
    crosswalk.to_csv(outdir / "new_vs_legacy_celltype_counts_v1.tsv", sep="\t", index=False)

    sns.set_theme(style="white", font_scale=0.9)
    umap = np.asarray(data.obsm["X_umap"])
    fig, ax = plt.subplots(figsize=(11, 8.5))
    for label in TYPE_ORDER:
        mask = labels.eq(label).to_numpy()
        if not mask.any():
            continue
        ax.scatter(umap[mask, 0], umap[mask, 1], s=0.35, color=COLORS[label], alpha=0.65, rasterized=True, label=label)
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.set_title("Part 1 rebuilt broad cell annotation v1")
    ax.legend(title="Broad cell type", bbox_to_anchor=(1.01, 1), loc="upper left", markerscale=8)
    save_both(fig, outdir / "umap_cell_type_broad_v1")

    fig, ax = plt.subplots(figsize=(9, 8))
    confidence_colors = {"High": "#009E73", "Medium": "#E69F00", "Low": "#777777"}
    for level in ["High", "Medium", "Low"]:
        mask = confidences.eq(level).to_numpy()
        ax.scatter(umap[mask, 0], umap[mask, 1], s=0.35, color=confidence_colors[level], alpha=0.65, rasterized=True, label=level)
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.set_title("Broad annotation confidence v1")
    ax.legend(title="Confidence", bbox_to_anchor=(1.01, 1), loc="upper left", markerscale=8)
    save_both(fig, outdir / "umap_annotation_confidence_v1")

    plot_counts = count_table[count_table["n_cells"] > 0].sort_values("n_cells")
    fig, ax = plt.subplots(figsize=(9, 6.5))
    ax.barh(plot_counts["cell_type_broad_v1"].astype(str), plot_counts["n_cells"], color=[COLORS[x] for x in plot_counts["cell_type_broad_v1"].astype(str)])
    ax.set_xlabel("Number of retained cells")
    ax.set_ylabel("Broad cell type")
    ax.set_title("Part 1 broad cell-type composition v1")
    for i, value in enumerate(plot_counts["n_cells"]):
        ax.text(value, i, f" {value:,}", va="center", fontsize=8)
    save_both(fig, outdir / "celltype_counts_v1")

    evidence_hashes = {
        path.name: sha256_file(path)
        for path in [EVIDENCE / "annotation_marker_panel_v2.tsv", EVIDENCE / "cluster_celltype_scores_v2.tsv", EVIDENCE / "cluster_top_specific_genes_v2.tsv"]
    }
    data.uns["cell_type_broad_v1_colors"] = [COLORS[x] for x in TYPE_ORDER]
    data.uns["part1_cell_annotation_v1"] = {
        "created_at": datetime.now().astimezone().isoformat(),
        "source_h5ad": str(INPUT),
        "cluster_field": "leiden04",
        "annotation_field": "cell_type_broad_v1",
        "evidence_directory": str(EVIDENCE),
        "evidence_sha256": evidence_hashes,
        "mapping_json": json.dumps(mapping.to_dict(orient="records"), ensure_ascii=False),
        "strategy": "reviewed canonical marker evidence and cluster-specific genes; legacy annotations comparison-only",
        "legacy_labels_copied": False,
        "ambiguous_policy": "retain Epithelial/Myeloid mixed for cluster 23; do not force lineage",
    }

    memory_gate("before_annotated_h5ad_write", memory_records)
    output_h5ad = outdir / "part1_scvi_integrated_cell_annotated_v1.h5ad"
    temporary_h5ad = outdir / "part1_scvi_integrated_cell_annotated_v1.tmp.h5ad"
    data.write_h5ad(temporary_h5ad, compression="gzip", compression_opts=4)
    os.replace(temporary_h5ad, output_h5ad)
    output_sha = sha256_file(output_h5ad)
    (outdir / "part1_scvi_integrated_cell_annotated_v1.h5ad.sha256").write_text(f"{output_sha}  {output_h5ad.name}\n", encoding="utf-8")
    metadata = {
        "run_id": "part1_cell_annotation_v1",
        "started_at": started,
        "completed_at": datetime.now().astimezone().isoformat(),
        "python_executable": sys.executable,
        "command": " ".join(sys.argv),
        "input_h5ad": str(INPUT),
        "output_h5ad": str(output_h5ad),
        "output_h5ad_sha256": output_sha,
        "shape": [int(data.n_obs), int(data.n_vars)],
        "annotation_field": "cell_type_broad_v1",
        "n_cell_types": int((count_table["n_cells"] > 0).sum()),
        "n_ambiguous_cells": int(data.obs["annotation_ambiguous_v1"].sum()),
        "legacy_labels_copied": False,
        "evidence_hashes": evidence_hashes,
        "memory_gate": memory_records,
    }
    (outdir / "run_metadata_v1.json").write_text(json.dumps(metadata, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    (outdir / "RUN_COMPLETE").write_text(metadata["completed_at"] + "\n", encoding="utf-8")
    print(json.dumps(metadata, indent=2, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
