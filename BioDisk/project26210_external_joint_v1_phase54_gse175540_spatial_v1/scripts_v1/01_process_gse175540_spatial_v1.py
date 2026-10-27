#!/usr/bin/env python3
"""Independent processing and targeted spatial screening of 24 GSE175540 sections."""

from __future__ import annotations

import argparse
import concurrent.futures
import gzip
import hashlib
import json
import os
import platform
import sys
import traceback
from datetime import datetime
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "2")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
os.environ.setdefault("MKL_NUM_THREADS", "2")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "2")

import anndata as ad
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy
from scipy import sparse
from scipy.spatial import cKDTree
from scipy.stats import fisher_exact, mannwhitneyu, spearmanr, wilcoxon


RELATIONS = [
    ("B01_F1_to_POSTN_CAF", "tumor_to_caf", "f1_top100_weighted", "postn_ecm_proxy", 1),
    ("B02_TYMP_to_POSTN_CAF", "tumor_to_caf", "tymp_expression", "postn_ecm_proxy", 1),
    ("B03_POSTN_CAF_to_APOE_TAM", "caf_to_tam", "postn_ecm_proxy", "apoe_tam_expression_proxy", 1),
    ("B04_POSTN_CAF_to_C1QC", "caf_to_tam", "postn_ecm_proxy", "c1qc_state", 1),
    ("B05_POSTN_CAF_to_SPP1", "caf_to_tam", "postn_ecm_proxy", "spp1_state_no_spp1", 1),
    ("B06_POSTN_CAF_to_myeloid_suppression", "caf_to_tam", "postn_ecm_proxy", "myeloid_immunosuppression", 1),
    ("B07_APOE_TAM_to_Tdysfunction", "tam_to_tcell", "apoe_tam_expression_proxy", "tcell_dysfunction", 1),
    ("B08_APOE_TAM_to_cytotoxicity", "tam_to_tcell", "apoe_tam_expression_proxy", "cytotoxicity", -1),
    ("B09_C1QC_to_Tdysfunction", "tam_to_tcell", "c1qc_state", "tcell_dysfunction", 1),
    ("B10_C1QC_to_cytotoxicity", "tam_to_tcell", "c1qc_state", "cytotoxicity", -1),
    ("B11_SPP1_to_Tdysfunction", "tam_to_tcell", "spp1_state_no_spp1", "tcell_dysfunction", 1),
    ("B12_SPP1_to_cytotoxicity", "tam_to_tcell", "spp1_state_no_spp1", "cytotoxicity", -1),
    ("B13_APOE_TAM_to_vascular", "tam_to_vascular", "apoe_tam_expression_proxy", "vascular_remodeling", 1),
    ("B14_SPP1_to_vascular", "tam_to_vascular", "spp1_state_no_spp1", "vascular_remodeling", 1),
    ("B15_F1_to_APOE_TAM", "tumor_to_tam", "f1_top100_weighted", "apoe_tam_expression_proxy", 1),
    ("B16_F1_to_C1QC", "tumor_to_tam", "f1_top100_weighted", "c1qc_state", 1),
    ("B17_F1_to_SPP1", "tumor_to_tam", "f1_top100_weighted", "spp1_state_no_spp1", 1),
    ("B18_TYMP_to_APOE_TAM", "tumor_to_tam", "tymp_expression", "apoe_tam_expression_proxy", 1),
    ("B19_TYMP_to_C1QC", "tumor_to_tam", "tymp_expression", "c1qc_state", 1),
    ("B20_TYMP_to_SPP1", "tumor_to_tam", "tymp_expression", "spp1_state_no_spp1", 1),
    ("B21_APOE_TAM_to_myeloid_suppression", "tam_state_internal", "apoe_tam_expression_proxy", "myeloid_immunosuppression", 1),
    ("B22_SPP1_to_myeloid_suppression", "tam_state_internal", "spp1_state_no_spp1", "myeloid_immunosuppression", 1),
]

PLOT_SIGNATURES = [
    "f1_top100_weighted", "tymp_expression", "postn_ecm_proxy",
    "apoe_tam_expression_proxy", "spp1_state_no_spp1",
    "myeloid_immunosuppression", "tcell_dysfunction", "cytotoxicity",
    "vascular_remodeling",
]


def now() -> str:
    return datetime.now().astimezone().isoformat(timespec="seconds")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_tsv(frame: pd.DataFrame, path: Path) -> None:
    tmp = path.with_name(path.name + ".part")
    compression = "gzip" if path.suffix == ".gz" else None
    frame.to_csv(
        tmp, sep="\t", index=False, na_rep="NA", encoding="utf-8",
        lineterminator="\n", compression=compression,
    )
    os.replace(tmp, path)


def bh(values: pd.Series) -> pd.Series:
    p = pd.to_numeric(values, errors="coerce").to_numpy(float)
    out = np.full(len(p), np.nan)
    good = np.isfinite(p)
    if not good.any():
        return pd.Series(out, index=values.index)
    vals = p[good]
    order = np.argsort(vals)
    ranked = vals[order] * len(vals) / np.arange(1, len(vals) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    tmp = np.empty(len(vals))
    tmp[order] = np.clip(ranked, 0, 1)
    out[np.where(good)[0]] = tmp
    return pd.Series(out, index=values.index)


def read_signatures(path: Path) -> dict[str, tuple[list[str], np.ndarray]]:
    table = pd.read_csv(path, sep="\t")
    required = {"signature_id", "gene", "raw_weight"}
    if not required.issubset(table.columns):
        raise ValueError(f"Signature table missing columns: {sorted(required - set(table.columns))}")
    result = {}
    for signature, group in table.groupby("signature_id", sort=False):
        genes = group["gene"].astype(str).str.upper().tolist()
        weights = pd.to_numeric(group["raw_weight"], errors="coerce").fillna(1.0).to_numpy(float)
        if not np.isfinite(weights).all() or np.abs(weights).sum() == 0:
            weights = np.ones(len(genes), dtype=float)
        result[str(signature)] = (genes, weights)
    return result


def parse_prefix(h5_path: Path) -> tuple[str, str, str]:
    name = h5_path.name
    suffix = "_filtered_feature_bc_matrix.h5"
    if not name.endswith(suffix):
        raise ValueError(name)
    prefix = name[:-len(suffix)]
    sample = prefix.split("_", 1)[0]
    preparation = "FFPE" if "_ffpe_" in prefix else "frozen"
    return prefix, sample, preparation


def dense_columns(matrix, indices: list[int]) -> np.ndarray:
    block = matrix[:, indices]
    if sparse.issparse(block):
        block = block.toarray()
    return np.asarray(block, dtype=np.float64)


def signature_scores(adata: ad.AnnData, signatures, outdir: Path) -> list[str]:
    lookup: dict[str, list[int]] = {}
    for idx, gene in enumerate(adata.var_names.astype(str).str.upper()):
        lookup.setdefault(gene, []).append(idx)
    coverage_rows = []
    calculated = []
    for signature, (genes, weights) in signatures.items():
        present_gene_positions = [i for i, gene in enumerate(genes) if gene in lookup]
        indices = [lookup[genes[i]][0] for i in present_gene_positions]
        present_genes = [genes[i] for i in present_gene_positions]
        missing_genes = [gene for gene in genes if gene not in lookup]
        coverage_rows.append({
            "signature_id": signature,
            "expected_genes": len(genes),
            "present_genes": len(indices),
            "coverage_fraction": len(indices) / len(genes),
            "present_gene_list": ";".join(present_genes) or "none",
            "missing_gene_list": ";".join(missing_genes) or "none",
        })
        if not indices:
            continue
        values = dense_columns(adata.X, indices)
        sensitivity = values.mean(axis=1)
        means = values.mean(axis=0)
        std = values.std(axis=0, ddof=0)
        std[std == 0] = 1.0
        z = (values - means) / std
        used_weights = weights[present_gene_positions]
        denominator = np.abs(used_weights).sum()
        primary = (z * used_weights).sum(axis=1) / denominator
        adata.obs[f"score_primary__{signature}"] = primary.astype(np.float32)
        adata.obs[f"score_sensitivity__{signature}"] = sensitivity.astype(np.float32)
        calculated.append(signature)
    atomic_tsv(pd.DataFrame(coverage_rows), outdir / "signature_coverage_v1.tsv")
    return calculated


def spatial_neighbors(coords: np.ndarray, k: int = 6) -> np.ndarray:
    if len(coords) <= k:
        raise ValueError(f"Need more than {k} spots for spatial neighbours")
    _, indices = cKDTree(coords).query(coords, k=k + 1)
    return np.asarray(indices[:, 1:], dtype=int)


def moran_i(values: np.ndarray, neighbors: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    if not np.isfinite(values).all() or np.nanstd(values) == 0:
        return np.nan
    z = values - values.mean()
    numerator = np.sum(z[:, None] * z[neighbors]) / neighbors.shape[1]
    denominator = np.sum(z * z)
    return float(len(values) * numerator / (len(values) * denominator))


def safe_spearman(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    valid = np.isfinite(x) & np.isfinite(y)
    if valid.sum() < 5 or np.std(x[valid]) == 0 or np.std(y[valid]) == 0:
        return np.nan, np.nan
    result = spearmanr(x[valid], y[valid])
    return float(result.statistic), float(result.pvalue)


def make_spatial_plot(adata: ad.AnnData, signatures: list[str], out: Path) -> None:
    coords = adata.obsm["spatial"]
    fig, axes = plt.subplots(3, 3, figsize=(15, 14), constrained_layout=True)
    for ax, signature in zip(axes.flat, PLOT_SIGNATURES):
        key = f"score_primary__{signature}"
        if key not in adata.obs:
            ax.text(0.5, 0.5, "signature unavailable", ha="center", va="center")
            ax.set_title(signature)
            ax.axis("off")
            continue
        points = ax.scatter(coords[:, 0], coords[:, 1], c=adata.obs[key], s=5, cmap="coolwarm", linewidths=0)
        ax.invert_yaxis()
        ax.set_title(signature, fontsize=9)
        ax.set_xticks([])
        ax.set_yticks([])
        fig.colorbar(points, ax=ax, fraction=0.035, pad=0.01)
    fig.savefig(out, dpi=180, bbox_inches="tight")
    plt.close(fig)


def make_embedding_plot(adata: ad.AnnData, out: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
    clusters = adata.obs["leiden_0_5"].astype("category")
    codes = clusters.cat.codes.to_numpy()
    axes[0].scatter(adata.obsm["X_umap"][:, 0], adata.obsm["X_umap"][:, 1], c=codes, s=4, cmap="tab20", linewidths=0)
    axes[0].set_title("UMAP / Leiden 0.5")
    axes[0].set_xticks([]); axes[0].set_yticks([])
    coords = adata.obsm["spatial"]
    axes[1].scatter(coords[:, 0], coords[:, 1], c=codes, s=5, cmap="tab20", linewidths=0)
    axes[1].invert_yaxis(); axes[1].set_title("Spatial clusters")
    axes[1].set_xticks([]); axes[1].set_yticks([])
    fig.savefig(out, dpi=180, bbox_inches="tight")
    plt.close(fig)


def process_section(h5_text: str, signatures_path: str, sections_dir_text: str) -> dict:
    h5_path = Path(h5_text)
    input_dir = h5_path.parent
    sections_dir = Path(sections_dir_text)
    prefix, sample, preparation = parse_prefix(h5_path)
    final = sections_dir / sample
    part = sections_dir / f"{sample}.part"
    if (final / "SECTION_PASS_v1").exists():
        return json.loads((final / "section_summary_v1.json").read_text(encoding="utf-8"))
    if final.exists() or part.exists():
        raise FileExistsError(f"Refusing existing incomplete output for {sample}")
    part.mkdir(parents=True)
    started = now()
    try:
        adata = sc.read_10x_h5(h5_path, genome=None, gex_only=True)
        adata.var_names = adata.var_names.astype(str).str.upper()
        adata.var_names_make_unique()
        adata.obs_names = adata.obs_names.astype(str)
        adata.X = sparse.csr_matrix(adata.X)
        data = adata.X.data
        if data.size and (not np.isfinite(data).all() or data.min() < 0 or not np.equal(data, np.floor(data)).all()):
            raise ValueError("Counts are negative, non-finite, or non-integer")
        original_spots = adata.n_obs
        original_genes = adata.n_vars
        original_nnz = adata.X.nnz
        original_total = int(adata.X.sum())

        positions_path = input_dir / f"{prefix}_tissue_positions_list.csv.gz"
        positions = pd.read_csv(
            positions_path, compression="gzip", header=None,
            names=["barcode", "in_tissue", "array_row", "array_col", "pxl_row", "pxl_col"],
        ).set_index("barcode")
        missing_positions = adata.obs_names.difference(positions.index)
        if len(missing_positions):
            raise ValueError(f"Missing positions for {len(missing_positions)} filtered barcodes")
        positions = positions.loc[adata.obs_names]
        for column in positions.columns:
            adata.obs[column] = positions[column].to_numpy()
        adata.obsm["spatial"] = positions[["pxl_col", "pxl_row"]].to_numpy(float)

        tls_path = input_dir / f"{prefix}_TLS_annotation.csv.gz"
        if tls_path.exists():
            tls = pd.read_csv(tls_path, compression="gzip")
            tls = tls.set_index(tls.columns[0]).iloc[:, 0]
            adata.obs["TLS_annotation"] = tls.reindex(adata.obs_names).fillna("UNANNOTATED").astype(str).to_numpy()
        else:
            adata.obs["TLS_annotation"] = "NOT_PROVIDED"

        adata.var["mt"] = adata.var_names.str.startswith("MT-")
        adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
        adata.var["hb"] = adata.var_names.str.match(r"^HB[ABDEGQZ]")
        sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], percent_top=None, log1p=False, inplace=True)
        qc_pass = (
            (adata.obs["n_genes_by_counts"] >= 200)
            & (adata.obs["total_counts"] >= 500)
            & (adata.obs["pct_counts_mt"] <= 30)
        )
        adata.obs["qc_pass"] = qc_pass.to_numpy(bool)
        qc_columns = [
            "total_counts", "n_genes_by_counts", "pct_counts_mt", "pct_counts_ribo",
            "pct_counts_hb", "in_tissue", "array_row", "array_col", "pxl_row", "pxl_col",
            "TLS_annotation", "qc_pass",
        ]
        qc_table = adata.obs[qc_columns].copy()
        qc_table.insert(0, "barcode", qc_table.index)
        atomic_tsv(qc_table, part / "spot_qc_all_v1.tsv.gz")
        adata = adata[qc_pass].copy()
        if adata.n_obs < 200:
            raise ValueError(f"Only {adata.n_obs} spots pass QC")

        adata.layers["counts"] = adata.X.copy()
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        signatures = read_signatures(Path(signatures_path))
        calculated = signature_scores(adata, signatures, part)

        sc.pp.highly_variable_genes(adata, n_top_genes=min(3000, adata.n_vars), flavor="seurat", subset=False)
        hvg = adata.var["highly_variable"].to_numpy(bool)
        pca_data = adata[:, hvg].copy()
        sc.pp.scale(pca_data, zero_center=True, max_value=10)
        n_comps = min(50, pca_data.n_obs - 1, pca_data.n_vars - 1)
        sc.tl.pca(pca_data, n_comps=n_comps, svd_solver="arpack", random_state=0)
        adata.obsm["X_pca"] = pca_data.obsm["X_pca"].astype(np.float32)
        adata.uns["pca"] = pca_data.uns["pca"]
        adata.varm["PCs"] = np.zeros((adata.n_vars, n_comps), dtype=np.float32)
        adata.varm["PCs"][hvg, :] = pca_data.varm["PCs"].astype(np.float32)
        del pca_data
        sc.pp.neighbors(adata, n_neighbors=15, n_pcs=min(30, n_comps), use_rep="X_pca", random_state=0)
        sc.tl.umap(adata, random_state=0)
        sc.tl.leiden(adata, resolution=0.5, key_added="leiden_0_5", random_state=0, flavor="igraph", directed=False, n_iterations=2)

        neighbors = spatial_neighbors(adata.obsm["spatial"], k=6)
        moran_rows = []
        tls_rows = []
        for method in ["primary", "sensitivity"]:
            for signature in calculated:
                key = f"score_{method}__{signature}"
                values = adata.obs[key].to_numpy(float)
                moran_rows.append({
                    "sample_id": sample, "preparation": preparation, "method": method,
                    "signature_id": signature, "n_spots": adata.n_obs, "moran_i_k6": moran_i(values, neighbors),
                })
                labels = adata.obs["TLS_annotation"].astype(str).to_numpy()
                tls_values = values[labels == "TLS"]
                non_values = values[labels == "NO_TLS"]
                if len(tls_values) >= 5 and len(non_values) >= 5:
                    test = mannwhitneyu(tls_values, non_values, alternative="two-sided")
                    delta = 2 * float(test.statistic) / (len(tls_values) * len(non_values)) - 1
                    p_value = float(test.pvalue)
                else:
                    delta = np.nan
                    p_value = np.nan
                tls_rows.append({
                    "sample_id": sample, "preparation": preparation, "method": method,
                    "signature_id": signature, "n_tls": len(tls_values), "n_no_tls": len(non_values),
                    "median_tls": np.median(tls_values) if len(tls_values) else np.nan,
                    "median_no_tls": np.median(non_values) if len(non_values) else np.nan,
                    "cliffs_delta_tls_vs_no_tls": delta, "mannwhitney_p": p_value,
                })

        relation_rows = []
        for method in ["primary", "sensitivity"]:
            for relation, edge_group, x_sig, y_sig, expected in RELATIONS:
                x_key = f"score_{method}__{x_sig}"
                y_key = f"score_{method}__{y_sig}"
                if x_key not in adata.obs or y_key not in adata.obs:
                    continue
                x = adata.obs[x_key].to_numpy(float)
                y = adata.obs[y_key].to_numpy(float)
                neighbor_y = np.mean(y[neighbors], axis=1)
                global_rho, global_p = safe_spearman(x, y)
                spatial_rho, spatial_p = safe_spearman(x, neighbor_y)
                x_hot = x >= np.quantile(x, 0.8)
                y_hot = neighbor_y >= np.quantile(neighbor_y, 0.8)
                table = np.array([
                    [np.sum(x_hot & y_hot), np.sum(x_hot & ~y_hot)],
                    [np.sum(~x_hot & y_hot), np.sum(~x_hot & ~y_hot)],
                ])
                odds, fisher_p = fisher_exact(table, alternative="greater")
                relation_rows.append({
                    "sample_id": sample, "preparation": preparation, "method": method,
                    "relation": relation, "edge_group": edge_group, "x_signature": x_sig,
                    "y_signature": y_sig, "expected_direction": expected, "n_spots": adata.n_obs,
                    "global_spearman_rho": global_rho, "global_spearman_p": global_p,
                    "neighbor_spearman_rho": spatial_rho, "neighbor_spearman_p": spatial_p,
                    "top20_neighbor_overlap_odds_ratio": float(odds), "top20_neighbor_overlap_fisher_p": float(fisher_p),
                    "direction_global_expected": bool(np.sign(global_rho) == expected) if np.isfinite(global_rho) else False,
                    "direction_neighbor_expected": bool(np.sign(spatial_rho) == expected) if np.isfinite(spatial_rho) else False,
                })

        atomic_tsv(pd.DataFrame(moran_rows), part / "spatial_autocorrelation_v1.tsv")
        atomic_tsv(pd.DataFrame(tls_rows), part / "tls_signature_comparison_v1.tsv")
        atomic_tsv(pd.DataFrame(relation_rows), part / "spatial_bridge_relations_v1.tsv")
        obs = adata.obs.copy()
        obs.insert(0, "barcode", obs.index)
        atomic_tsv(obs, part / "spot_metadata_scores_v1.tsv.gz")
        make_spatial_plot(adata, calculated, part / "target_signature_spatial_maps_v1.png")
        make_embedding_plot(adata, part / "embedding_and_cluster_maps_v1.png")
        adata.uns["phase54"] = {
            "sample_id": sample, "preparation": preparation, "source_h5": h5_path.name,
            "source_h5_sha256": sha256(h5_path), "normalization": "normalize_total_1e4_log1p",
            "qc_rule": "n_genes>=200,total_counts>=500,pct_mt<=30",
        }
        h5ad_part = part / "processed_section_v1.h5ad.part"
        adata.write_h5ad(h5ad_part, compression="gzip")
        os.replace(h5ad_part, part / "processed_section_v1.h5ad")

        summary = {
            "status": "PASS", "sample_id": sample, "preparation": preparation,
            "started_at": started, "completed_at": now(), "source_h5": h5_path.name,
            "source_h5_sha256": sha256(h5_path), "original_spots": original_spots,
            "qc_pass_spots": adata.n_obs, "original_genes": original_genes,
            "original_nnz": original_nnz, "original_total_counts": original_total,
            "signatures_calculated": len(calculated), "clusters": int(adata.obs["leiden_0_5"].nunique()),
            "tls_annotation_available": tls_path.exists(),
        }
        (part / "section_summary_v1.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
        (part / "SECTION_PASS_v1").write_text(f"status=PASS\nsample_id={sample}\ncompleted_at={summary['completed_at']}\n", encoding="utf-8")
        os.replace(part, final)
        return summary
    except Exception:
        (part / "SECTION_FAILURE_v1.txt").write_text(traceback.format_exc(), encoding="utf-8")
        raise


def build_meta(relations: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for keys, group in relations.groupby(["preparation", "method", "relation", "edge_group", "expected_direction"], dropna=False):
        prep, method, relation, edge_group, expected = keys
        for metric in ["global_spearman_rho", "neighbor_spearman_rho"]:
            values = pd.to_numeric(group[metric], errors="coerce").dropna().to_numpy(float)
            if len(values) >= 5 and np.any(values != 0):
                try:
                    p_value = float(wilcoxon(values, alternative="two-sided", zero_method="wilcox").pvalue)
                except ValueError:
                    p_value = np.nan
            else:
                p_value = np.nan
            rows.append({
                "preparation": prep, "method": method, "relation": relation,
                "edge_group": edge_group, "expected_direction": expected, "metric": metric,
                "n_sections": len(values), "median_rho": np.median(values) if len(values) else np.nan,
                "min_rho": np.min(values) if len(values) else np.nan,
                "max_rho": np.max(values) if len(values) else np.nan,
                "expected_direction_fraction": np.mean(np.sign(values) == expected) if len(values) else np.nan,
                "wilcoxon_signed_rank_p": p_value,
            })
    meta = pd.DataFrame(rows)
    meta["wilcoxon_bh_q_within_preparation_method_metric"] = meta.groupby(
        ["preparation", "method", "metric"], group_keys=False
    )["wilcoxon_signed_rank_p"].apply(bh)
    return meta


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--signatures", required=True, type=Path)
    parser.add_argument("--workers", type=int, default=6)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    sections_dir = args.output_dir / "sections_v1"
    sections_dir.mkdir(exist_ok=True)
    h5_files = sorted(args.input_dir.glob("GSM*_filtered_feature_bc_matrix.h5"))
    if len(h5_files) != 24:
        raise RuntimeError(f"Expected 24 filtered H5 matrices, found {len(h5_files)}")
    samples = [parse_prefix(path)[1] for path in h5_files]
    if len(set(samples)) != 24:
        raise RuntimeError("Sample identifiers are not unique")

    summaries = []
    failures = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = {
            pool.submit(process_section, str(path), str(args.signatures), str(sections_dir)): parse_prefix(path)[1]
            for path in h5_files
        }
        for future in concurrent.futures.as_completed(futures):
            sample = futures[future]
            try:
                summary = future.result()
                summaries.append(summary)
                print(f"[{now()}] PASS {sample} spots={summary['qc_pass_spots']}", flush=True)
            except Exception as exc:
                failures.append({"sample_id": sample, "error": repr(exc)})
                print(f"[{now()}] FAIL {sample}: {exc!r}", flush=True)

    summary_table = pd.DataFrame(summaries).sort_values("sample_id") if summaries else pd.DataFrame()
    atomic_tsv(summary_table, args.output_dir / "section_status_v1.tsv")
    atomic_tsv(pd.DataFrame(failures, columns=["sample_id", "error"]), args.output_dir / "section_failures_v1.tsv")
    if summaries:
        relation_tables = [pd.read_csv(sections_dir / row["sample_id"] / "spatial_bridge_relations_v1.tsv", sep="\t") for row in summaries]
        moran_tables = [pd.read_csv(sections_dir / row["sample_id"] / "spatial_autocorrelation_v1.tsv", sep="\t") for row in summaries]
        tls_tables = [pd.read_csv(sections_dir / row["sample_id"] / "tls_signature_comparison_v1.tsv", sep="\t") for row in summaries]
        relations = pd.concat(relation_tables, ignore_index=True)
        for p_column in ["global_spearman_p", "neighbor_spearman_p", "top20_neighbor_overlap_fisher_p"]:
            relations[p_column.replace("_p", "_bh_q")] = relations.groupby("method", group_keys=False)[p_column].apply(bh)
        moran = pd.concat(moran_tables, ignore_index=True)
        tls = pd.concat(tls_tables, ignore_index=True)
        tls["mannwhitney_bh_q_within_method"] = tls.groupby("method", group_keys=False)["mannwhitney_p"].apply(bh)
        meta = build_meta(relations)
        atomic_tsv(relations, args.output_dir / "all_section_spatial_bridge_relations_v1.tsv")
        atomic_tsv(moran, args.output_dir / "all_section_spatial_autocorrelation_v1.tsv")
        atomic_tsv(tls, args.output_dir / "all_section_tls_signature_comparison_v1.tsv")
        atomic_tsv(meta, args.output_dir / "cross_section_relation_meta_v1.tsv")

    checks = pd.DataFrame([
        {"check": "filtered_h5_count", "status": "PASS" if len(h5_files) == 24 else "FAIL", "observed": len(h5_files), "expected": 24},
        {"check": "unique_sample_count", "status": "PASS" if len(set(samples)) == 24 else "FAIL", "observed": len(set(samples)), "expected": 24},
        {"check": "section_pass_count", "status": "PASS" if len(summaries) == 24 else "FAIL", "observed": len(summaries), "expected": 24},
        {"check": "section_failure_count", "status": "PASS" if len(failures) == 0 else "FAIL", "observed": len(failures), "expected": 0},
        {"check": "ffpe_count", "status": "PASS" if sum(x["preparation"] == "FFPE" for x in summaries) == 12 else "FAIL", "observed": sum(x["preparation"] == "FFPE" for x in summaries), "expected": 12},
        {"check": "frozen_count", "status": "PASS" if sum(x["preparation"] == "frozen" for x in summaries) == 12 else "FAIL", "observed": sum(x["preparation"] == "frozen" for x in summaries), "expected": 12},
        {"check": "tls_annotation_sections", "status": "PASS" if sum(bool(x["tls_annotation_available"]) for x in summaries) == 23 else "FAIL", "observed": sum(bool(x["tls_annotation_available"]) for x in summaries), "expected": 23},
        {"check": "per_section_independent_processing", "status": "PASS", "observed": "24 separate outputs", "expected": "no merge/integration/batch correction"},
    ])
    atomic_tsv(checks, args.output_dir / "validation_checks_v1.tsv")
    status = "PASS" if (checks["status"] == "PASS").all() else "PARTIAL"
    metadata = {
        "status": status, "generated_at": now(), "workers": args.workers,
        "input_dir": str(args.input_dir), "signature_sha256": sha256(args.signatures),
        "sections_passed": len(summaries), "sections_failed": len(failures),
        "python": platform.python_version(), "scanpy": sc.__version__, "anndata": ad.__version__,
        "numpy": np.__version__, "pandas": pd.__version__, "scipy": scipy.__version__,
        "normalization": "per-section normalize_total(1e4) then log1p",
        "integration_performed": False, "batch_correction_performed": False,
        "deconvolution_performed": False,
    }
    (args.output_dir / "run_metadata_v1.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    marker = args.output_dir / ("PHASE54_SPATIAL_PROCESSING_PASS_v1" if status == "PASS" else "PHASE54_SPATIAL_PROCESSING_PARTIAL_v1")
    marker.write_text(f"status={status}\nsections_passed={len(summaries)}\nsections_failed={len(failures)}\n", encoding="utf-8")
    return 0 if status == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
