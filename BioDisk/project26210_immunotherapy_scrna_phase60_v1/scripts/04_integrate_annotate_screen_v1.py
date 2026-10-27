#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gc
import json
import math
import os
import traceback
from datetime import datetime
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from scipy import stats


SEED = 2621061
HIF1A_FROZEN = ["CA9", "VEGFA", "ADM", "ANGPTL4", "LDHA"]
MARKERS = {
    "Epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19", "KRT7", "CA9", "NDUFA4L2", "ALDOB", "FXYD2"],
    "T cell": ["CD3D", "CD3E", "TRBC1", "TRBC2", "IL7R", "CD247"],
    "NK cell": ["NKG7", "GNLY", "KLRD1", "PRF1", "GZMB"],
    "Myeloid": ["LST1", "TYROBP", "FCER1G", "CTSS", "AIF1", "LYZ"],
    "B cell": ["MS4A1", "CD79A", "CD37", "CD74", "HLA-DRA"],
    "Plasma cell": ["JCHAIN", "MZB1", "SDC1", "IGHG1", "DERL3"],
    "Endothelial": ["PECAM1", "VWF", "KDR", "EMCN", "ENG", "RAMP2"],
    "Fibroblast": ["COL1A1", "COL1A2", "DCN", "COL3A1", "LUM", "CFD"],
    "Pericyte/SMC": ["RGS5", "MCAM", "CSPG4", "ACTA2", "PDGFRB", "TAGLN"],
    "Mast cell": ["TPSAB1", "TPSB2", "KIT", "CPA3", "MS4A2"],
}


def now() -> str:
    return datetime.now().astimezone().isoformat(timespec="seconds")


def log(message: str) -> None:
    print(f"[{now()}] {message}", flush=True)


def memory_gate(stage: str, minimum_gib: float = 100.0) -> None:
    with Path("/proc/meminfo").open() as handle:
        fields = {line.split(":", 1)[0]: line.split()[1] for line in handle}
    available = float(fields["MemAvailable"]) / 1024.0 / 1024.0
    log(f"MEMORY_GATE stage={stage} available_gib={available:.2f} minimum_gib={minimum_gib:.2f}")
    if available < minimum_gib:
        raise RuntimeError(f"Memory gate failed at {stage}: {available:.2f} GiB")


def vector(obj: ad.AnnData, gene: str) -> np.ndarray:
    if gene not in obj.var_names:
        return np.full(obj.n_obs, np.nan)
    values = obj[:, gene].X
    if sp.issparse(values):
        values = values.toarray()
    return np.asarray(values, dtype=float).ravel()


def mean_score(obj: ad.AnnData, genes: list[str]) -> tuple[np.ndarray, list[str]]:
    present = [gene for gene in genes if gene in obj.var_names]
    if not present:
        return np.full(obj.n_obs, np.nan), []
    values = obj[:, present].X
    score = np.asarray(values.mean(axis=1)).ravel() if sp.issparse(values) else np.asarray(values).mean(axis=1)
    return score.astype(float), present


def f1_projection(obj: ad.AnnData, signature: pd.DataFrame) -> tuple[np.ndarray, dict]:
    signature = signature.loc[signature["factor"].eq("F1")].sort_values("gene_rank").head(100).copy()
    signature = signature.loc[signature["gene"].isin(obj.var_names)]
    if len(signature) < 60:
        raise RuntimeError(f"Insufficient F1 signature coverage: {len(signature)}/100")
    weights = signature["loading"].to_numpy(float)
    weights = weights / weights.sum()
    values = obj[:, signature["gene"].tolist()].X
    score = values @ weights
    score = np.asarray(score).ravel().astype(float)
    return score, {
        "score_name": "F1_projection_score_v1", "authority": "rank7_program_top100_genes_v1.tsv",
        "n_signature_genes": 100, "n_present_genes": int(len(signature)),
        "coverage": float(len(signature) / 100),
        "not_original_cnmf_usage": True,
        "genes_present": signature["gene"].tolist(),
    }


def annotate_clusters(obj: ad.AnnData, cluster_key: str) -> pd.DataFrame:
    score_columns = []
    coverage = []
    for label, genes in MARKERS.items():
        score, present = mean_score(obj, genes)
        column = "marker_score__" + label.replace("/", "_").replace(" ", "_")
        obj.obs[column] = score
        score_columns.append(column)
        coverage.append({"cell_type": label, "n_genes": len(genes), "n_present": len(present), "genes_present": ";".join(present)})
    cluster = obj.obs.groupby(cluster_key, observed=True)[score_columns].mean()
    z = cluster.copy()
    for column in score_columns:
        values = cluster[column].to_numpy(float)
        sd = np.std(values)
        z[column] = (values - np.mean(values)) / sd if sd > 0 else 0.0
    label_lookup = {}
    confidence_lookup = {}
    for cluster_id, row in z.iterrows():
        ordered = row.sort_values(ascending=False)
        best = ordered.index[0]
        label = best.removeprefix("marker_score__").replace("_", " ")
        label = label.replace("Pericyte SMC", "Pericyte/SMC")
        label_lookup[str(cluster_id)] = label
        confidence_lookup[str(cluster_id)] = float(ordered.iloc[0] - ordered.iloc[1])
    obj.obs["cell_type_broad_preliminary_v1"] = obj.obs[cluster_key].astype(str).map(label_lookup).astype(str)
    obj.obs["annotation_margin_preliminary_v1"] = obj.obs[cluster_key].astype(str).map(confidence_lookup).astype(float)
    table = cluster.copy()
    table.columns = [x.removeprefix("marker_score__") for x in table.columns]
    table.insert(0, "cell_type_preliminary", [label_lookup[str(x)] for x in table.index])
    table.insert(1, "annotation_margin", [confidence_lookup[str(x)] for x in table.index])
    table.index.name = "leiden_cluster"
    return table.reset_index(), pd.DataFrame(coverage)


def patient_screen(obj: ad.AnnData, cohort: str, output: Path) -> None:
    epi = obj.obs["cell_type_broad_preliminary_v1"].eq("Epithelial").to_numpy()
    rows = []
    group_columns = ["patient_id", "timepoint", "response_binary", "response_detail", "treatment"]
    obs = obj.obs.reset_index(names="cell_id")
    for keys, group in obs.groupby(group_columns, observed=True, dropna=False):
        pos = obj.obs_names.get_indexer(group["cell_id"])
        epi_pos = pos[epi[pos]]
        if len(epi_pos) >= 20:
            x = obj.obs.iloc[epi_pos]["F1_projection_score_v1"].to_numpy(float)
            y = obj.obs.iloc[epi_pos]["TYMP_log1p_CP10K"].to_numpy(float)
            rho = float(stats.spearmanr(x, y).statistic) if np.ptp(x) > 0 and np.ptp(y) > 0 else math.nan
        else:
            rho = math.nan
        patient_id, timepoint, response, detail, treatment = keys
        rows.append({
            "cohort": cohort, "patient_id": patient_id, "timepoint": timepoint,
            "response_binary": response, "response_detail": detail, "treatment": treatment,
            "n_cells_total": len(pos), "n_epithelial_preliminary": len(epi_pos),
            "epithelial_fraction": float(len(epi_pos) / len(pos)),
            "median_F1_projection_epithelial": float(np.median(obj.obs.iloc[epi_pos]["F1_projection_score_v1"])) if len(epi_pos) else math.nan,
            "median_TYMP_epithelial": float(np.median(obj.obs.iloc[epi_pos]["TYMP_log1p_CP10K"])) if len(epi_pos) else math.nan,
            "median_hif1a_downstream_epithelial": float(np.median(obj.obs.iloc[epi_pos]["hif1a_downstream_5gene_score_v1"])) if len(epi_pos) else math.nan,
            "spearman_F1_projection_TYMP_epithelial": rho,
            "eligible_preliminary_epithelial_screen": len(epi_pos) >= 20,
        })
    summary = pd.DataFrame(rows).sort_values(["patient_id", "timepoint"])
    summary.to_csv(output / f"{cohort}.preliminary_patient_endpoint_summary_v1.tsv", sep="\t", index=False)

    composition = obs.groupby(["patient_id", "timepoint", "response_binary", "cell_type_broad_preliminary_v1"], observed=True).size().rename("n_cells").reset_index()
    totals = composition.groupby(["patient_id", "timepoint"], observed=True)["n_cells"].transform("sum")
    composition["fraction"] = composition["n_cells"] / totals
    composition.to_csv(output / f"{cohort}.preliminary_patient_cell_composition_v1.tsv", sep="\t", index=False)

    primary = summary.loc[summary["eligible_preliminary_epithelial_screen"]].copy()
    if cohort == "GSE264586":
        primary = primary.loc[primary["timepoint"].eq("Baseline")]
    else:
        primary = primary.loc[primary["response_binary"].isin(["Responder", "Nonresponder"])]
    test_rows = []
    endpoints = ["spearman_F1_projection_TYMP_epithelial", "median_F1_projection_epithelial",
                 "median_TYMP_epithelial", "median_hif1a_downstream_epithelial", "epithelial_fraction"]
    for endpoint in endpoints:
        responder = primary.loc[primary["response_binary"].eq("Responder"), endpoint].dropna().to_numpy(float)
        nonresponder = primary.loc[primary["response_binary"].eq("Nonresponder"), endpoint].dropna().to_numpy(float)
        pvalue = float(stats.mannwhitneyu(responder, nonresponder, alternative="two-sided", method="exact").pvalue) if len(responder) and len(nonresponder) else math.nan
        test_rows.append({
            "cohort": cohort, "endpoint": endpoint, "n_responder": len(responder),
            "n_nonresponder": len(nonresponder), "median_responder": float(np.median(responder)) if len(responder) else math.nan,
            "median_nonresponder": float(np.median(nonresponder)) if len(nonresponder) else math.nan,
            "median_difference_R_minus_NR": float(np.median(responder) - np.median(nonresponder)) if len(responder) and len(nonresponder) else math.nan,
            "mannwhitney_exact_p": pvalue, "evidence_status": "PRELIMINARY_EPITHELIAL_SCREEN_NOT_FORMAL",
        })
    pd.DataFrame(test_rows).to_csv(output / f"{cohort}.preliminary_response_screen_v1.tsv", sep="\t", index=False)


def render_umap(obj: ad.AnnData, cohort: str, output: Path) -> None:
    rng = np.random.default_rng(SEED)
    take = np.arange(obj.n_obs)
    if len(take) > 80000:
        take = np.sort(rng.choice(take, 80000, replace=False))
    coords = np.asarray(obj.obsm["X_umap"])[take]
    labels = obj.obs.iloc[take]["cell_type_broad_preliminary_v1"].astype(str)
    categories = sorted(labels.unique())
    colors = plt.get_cmap("tab20")(np.linspace(0, 1, max(1, len(categories))))
    fig, ax = plt.subplots(figsize=(8.5, 7.0))
    for label, color in zip(categories, colors):
        mask = labels.to_numpy() == label
        ax.scatter(coords[mask, 0], coords[mask, 1], s=1.0, alpha=0.55, linewidths=0,
                   rasterized=True, label=label, color=color)
    ax.set_title(f"{cohort} preliminary broad annotation")
    ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2")
    ax.legend(frameon=False, markerscale=7, bbox_to_anchor=(1.02, 1), loc="upper left")
    fig.tight_layout()
    fig.savefig(output / f"{cohort}.preliminary_broad_annotation_umap_v1.png", dpi=240, facecolor="white")
    plt.close(fig)


def process_cohort(cohort: str, qc_root: Path, signature_path: Path, output: Path) -> dict:
    destination = output / f"{cohort}.integrated_annotated_preliminary_v1.h5ad"
    checkpoint = output / f"{cohort}.integration_status_v1.json"
    if destination.is_file() and checkpoint.is_file():
        status = json.loads(checkpoint.read_text(encoding="utf-8"))
        if status.get("status") == "PASS":
            log(f"INTEGRATION_RESUME_SKIP {cohort}")
            return status
    memory_gate(f"integration_load_{cohort}")
    source = qc_root / "cohorts" / f"{cohort}.qc_filtered_raw_v1.h5ad"
    obj = ad.read_h5ad(source)
    obj.X = obj.layers["counts"].copy()
    sc.pp.normalize_total(obj, target_sum=10000)
    sc.pp.log1p(obj)
    signature = pd.read_csv(signature_path, sep="\t")
    f1_score, f1_meta = f1_projection(obj, signature)
    obj.obs["F1_projection_score_v1"] = f1_score
    obj.obs["TYMP_log1p_CP10K"] = vector(obj, "TYMP")
    hif1a, hif1a_present = mean_score(obj, HIF1A_FROZEN)
    obj.obs["hif1a_downstream_5gene_score_v1"] = hif1a
    for gene in ["ODC1", "SRM", "PCYT2"]:
        obj.obs[f"{gene}_log1p_CP10K"] = vector(obj, gene)
    obj.uns["F1_projection_definition_v1"] = f1_meta
    obj.uns["hif1a_downstream_definition_v1"] = {"genes": HIF1A_FROZEN, "present": hif1a_present}

    log(f"HVG_PCA_START {cohort} cells={obj.n_obs} genes={obj.n_vars}")
    sc.pp.highly_variable_genes(obj, n_top_genes=min(4000, obj.n_vars - 1), flavor="seurat", batch_key="library_id")
    sc.pp.pca(obj, n_comps=50, use_highly_variable=True, zero_center=False, random_state=SEED)
    sc.external.pp.harmony_integrate(obj, key="library_id", basis="X_pca", adjusted_basis="X_pca_harmony",
                                     max_iter_harmony=20, random_state=SEED)
    sc.pp.neighbors(obj, n_neighbors=15, n_pcs=40, use_rep="X_pca_harmony", random_state=SEED)
    sc.tl.umap(obj, min_dist=0.35, spread=1.0, random_state=SEED)
    sc.tl.leiden(obj, resolution=0.6, key_added="leiden06_preliminary_v1", random_state=SEED,
                 flavor="igraph", n_iterations=2, directed=False)
    cluster_table, marker_coverage = annotate_clusters(obj, "leiden06_preliminary_v1")
    cluster_table.to_csv(output / f"{cohort}.preliminary_cluster_marker_evidence_v1.tsv", sep="\t", index=False)
    marker_coverage.to_csv(output / f"{cohort}.marker_panel_coverage_v1.tsv", sep="\t", index=False)
    patient_screen(obj, cohort, output)
    render_umap(obj, cohort, output)
    obj.uns["phase60_preliminary_annotation_boundary_v1"] = {
        "status": "PRELIMINARY_BROAD_ANNOTATION",
        "malignant_gate_complete": False,
        "response_screen_formal": False,
        "corrected_embedding_not_used_for_expression_endpoints": True,
    }
    part = Path(str(destination) + ".part")
    obj.write_h5ad(part, compression="gzip")
    os.replace(part, destination)
    status = {
        "status": "PASS", "completed_at": now(), "cohort": cohort,
        "n_cells": int(obj.n_obs), "n_genes": int(obj.n_vars),
        "n_libraries": int(obj.obs["library_id"].nunique()),
        "n_patients": int(obj.obs["patient_id"].nunique()),
        "n_clusters": int(obj.obs["leiden06_preliminary_v1"].nunique()),
        "n_preliminary_epithelial": int(obj.obs["cell_type_broad_preliminary_v1"].eq("Epithelial").sum()),
        "f1_signature_coverage": f1_meta["coverage"],
        "h5ad": str(destination),
        "evidence_boundary": "PRELIMINARY_NOT_FORMAL_UNTIL_MALIGNANT_GATING",
    }
    checkpoint.write_text(json.dumps(status, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    log("INTEGRATION_PASS " + json.dumps(status, ensure_ascii=False))
    del obj
    gc.collect()
    return status


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--qc-root", type=Path, required=True)
    parser.add_argument("--f1-signature", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    statuses = [process_cohort(cohort, args.qc_root, args.f1_signature, args.output)
                for cohort in ["GSE264586", "GSE294109"]]
    final = {"status": "PASS", "completed_at": now(), "cohorts": statuses}
    (args.output / "INTEGRATION_ANNOTATION_SCREEN_COMPLETE_v1.json").write_text(
        json.dumps(final, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(final, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    try:
        main()
    except Exception:
        traceback.print_exc()
        raise
