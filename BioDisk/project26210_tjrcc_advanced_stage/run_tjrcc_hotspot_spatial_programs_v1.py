#!/usr/bin/env python3
"""Hotspot spatial autocorrelation, modules, and SPARK-X integration for one section."""

from __future__ import annotations

import hashlib
import importlib.metadata as metadata
import json
import math
import subprocess
import sys
import traceback
from pathlib import Path

import anndata as ad
import hotspot
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import sparse, stats


PROJECT = Path("/mnt/BioDisk/ludingyu/project26210")
RESULT_ROOT = PROJECT / "results/part4_spatial"
MEMORY_GATE = PROJECT / "check_memory_gate.sh"
PROGRAMS = {
    "Tumor-epithelial": {"CA9", "EPCAM", "KRT8", "KRT18"},
    "POSTN-ECM": {"POSTN", "COL1A1", "COL1A2", "COL3A1", "COL4A1", "FN1", "SPARC"},
    "APOE-TAM": {"APOE", "C1QA", "C1QB", "C1QC", "TREM2", "LPL"},
    "CD8-effector": {"CD8A", "NKG7", "CCL5", "PRF1", "GZMB"},
    "TGFb-response": {"TGFB1", "TGFBR1", "TGFBR2", "SMAD2", "SMAD3", "SERPINE1"},
    "TYMP-context": {"TYMP", "SPP1", "MIF", "CD74", "CXCR4", "APP"},
}


def sha256(path: Path, chunk_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def run_gate(stage: str) -> str:
    proc = subprocess.run([str(MEMORY_GATE)], check=True, text=True, capture_output=True)
    line = proc.stdout.strip()
    print(f"MEMORY_GATE stage={stage} {line}", flush=True)
    return line


def write_json(path: Path, payload: dict) -> None:
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def write_status(output_dir: Path, sample: str, status: str, note: str) -> None:
    write_json(
        output_dir / "status.json",
        {"status": status, "sample": sample, "note": note},
    )


def as_bool(series: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(series):
        return series
    return series.astype(str).str.lower().map({"true": True, "false": False}).astype(bool)


def plot_gene_concordance(table: pd.DataFrame, sample: str, output_dir: Path) -> None:
    palette = {
        "both_strict": "#7B2CBF",
        "sparkx_only": "#D95F02",
        "hotspot_only": "#1B9E77",
        "not_strict": "#BDBDBD",
    }
    fig, ax = plt.subplots(figsize=(8.2, 6.6))
    y = np.minimum(-np.log10(np.clip(table["sparkx_combined_p"].to_numpy(float), 1e-300, 1.0)), 300)
    x = table["hotspot_Z"].to_numpy(float)
    for label in ("not_strict", "sparkx_only", "hotspot_only", "both_strict"):
        mask = table["support_class"].eq(label).to_numpy()
        if mask.any():
            ax.scatter(x[mask], y[mask], s=10 if label == "not_strict" else 18, alpha=0.35 if label == "not_strict" else 0.75, c=palette[label], label=f"{label} (n={mask.sum()})", linewidths=0)
    labels = table.loc[table["support_class"].eq("both_strict")].nsmallest(15, "consensus_rank")
    for _, row in labels.iterrows():
        ax.text(row["hotspot_Z"], min(-math.log10(max(row["sparkx_combined_p"], 1e-300)), 300), row["gene"], fontsize=7)
    ax.set_xlabel("Hotspot autocorrelation Z")
    ax.set_ylabel("-log10 SPARK-X combined P")
    ax.set_title(f"{sample} gene-level spatial variability concordance")
    ax.legend(frameon=False, fontsize=8, loc="best")
    ax.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(output_dir / f"{sample}_sparkx_hotspot_concordance.png", dpi=300)
    fig.savefig(output_dir / f"{sample}_sparkx_hotspot_concordance.pdf")
    plt.close(fig)


def plot_module_maps(
    adata: ad.AnnData,
    scores: pd.DataFrame,
    module_sizes: pd.Series,
    sample: str,
    output_dir: Path,
) -> None:
    modules = module_sizes.sort_values(ascending=False).index.tolist()[:8]
    if not modules:
        fig, ax = plt.subplots(figsize=(7, 4))
        ax.text(0.5, 0.5, "No Hotspot module passed the prespecified threshold", ha="center", va="center")
        ax.axis("off")
    else:
        ncols = 2
        nrows = math.ceil(len(modules) / ncols)
        fig, axes = plt.subplots(nrows, ncols, figsize=(11, 4.8 * nrows), squeeze=False)
        spatial_key = sample if sample in adata.uns.get("spatial", {}) else next(iter(adata.uns["spatial"]))
        spatial_meta = adata.uns["spatial"][spatial_key]
        image = spatial_meta["images"]["hires"]
        scale = float(spatial_meta["scalefactors"]["tissue_hires_scalef"])
        coords = np.asarray(adata.obsm["spatial"], dtype=float) * scale
        for ax, module in zip(axes.flat, modules):
            values = scores[module].to_numpy(float)
            limit = float(np.nanpercentile(np.abs(values), 98))
            if not np.isfinite(limit) or limit == 0:
                limit = 1.0
            ax.imshow(image)
            sc = ax.scatter(coords[:, 0], coords[:, 1], c=values, s=7, cmap="coolwarm", vmin=-limit, vmax=limit, alpha=0.88, linewidths=0)
            ax.set_title(f"Module {module} (n={int(module_sizes.loc[module])})")
            ax.axis("off")
            fig.colorbar(sc, ax=ax, fraction=0.035, pad=0.01, label="Hotspot module score")
        for ax in axes.flat[len(modules):]:
            ax.axis("off")
        fig.suptitle(f"{sample} Hotspot spatial modules", fontsize=15)
        fig.tight_layout(rect=(0, 0, 1, 0.98))
    fig.savefig(output_dir / f"{sample}_hotspot_module_spatial_maps.png", dpi=300, bbox_inches="tight")
    fig.savefig(output_dir / f"{sample}_hotspot_module_spatial_maps.pdf", bbox_inches="tight")
    plt.close(fig)


def plot_module_summary(module_summary: pd.DataFrame, sample: str, output_dir: Path) -> None:
    if module_summary.empty:
        fig, ax = plt.subplots(figsize=(7, 4))
        ax.text(0.5, 0.5, "No Hotspot module passed the prespecified threshold", ha="center", va="center")
        ax.axis("off")
    else:
        plot_data = module_summary.set_index("module")[["n_genes", "n_both_strict", "distance_spearman_rho_descriptive", "tumor_minus_nontumor_median"]].copy()
        plot_data["n_genes"] = stats.zscore(plot_data["n_genes"], nan_policy="omit") if len(plot_data) > 1 else 0.0
        plot_data["n_both_strict"] = stats.zscore(plot_data["n_both_strict"], nan_policy="omit") if len(plot_data) > 1 else 0.0
        plot_data = plot_data.rename(columns={"n_genes": "Module size (z)", "n_both_strict": "Consensus genes (z)", "distance_spearman_rho_descriptive": "Distance rho", "tumor_minus_nontumor_median": "Tumor-minus-nontumor"}).fillna(0.0)
        fig, ax = plt.subplots(figsize=(8, max(3.5, 0.55 * len(plot_data) + 2)))
        sns.heatmap(plot_data, cmap="vlag", center=0, annot=True, fmt=".2f", linewidths=0.4, ax=ax, cbar_kws={"label": "Standardized or native effect"})
        ax.set_title(f"{sample} Hotspot module summary")
        ax.set_ylabel("Module")
        fig.tight_layout()
    fig.savefig(output_dir / f"{sample}_hotspot_module_summary.png", dpi=300, bbox_inches="tight")
    fig.savefig(output_dir / f"{sample}_hotspot_module_summary.pdf", bbox_inches="tight")
    plt.close(fig)


def main(sample: str, input_h5ad: Path, output_dir: Path, seed: int, boundary_path: Path) -> int:
    if not sample.replace("_", "").isalnum() or seed < 1:
        raise ValueError("Invalid sample ID or seed")
    output_dir.mkdir(parents=True, exist_ok=True)
    write_status(output_dir, sample, "running", "Hotspot and SPARK-X integration")
    memory_records = [run_gate(f"{sample.lower()}_hotspot_start")]

    spark_path = output_dir / f"{sample}_sparkx_gene_results.tsv"
    spark_meta_path = output_dir / "sparkx_run_metadata.json"
    if not spark_path.is_file() or not spark_meta_path.is_file():
        raise FileNotFoundError("Required completed SPARK-X outputs are missing")
    spark_meta = json.loads(spark_meta_path.read_text())
    if spark_meta.get("status") != "completed" or spark_meta.get("input_sha256") != sha256(input_h5ad):
        raise RuntimeError("SPARK-X metadata/input fingerprint validation failed")

    memory_records.append(run_gate(f"{sample.lower()}_hotspot_read"))
    adata = ad.read_h5ad(input_h5ad)
    if adata.n_obs < 100 or adata.n_vars < 2000:
        raise RuntimeError(f"Unexpected H5AD shape: {adata.shape}")
    if not sparse.issparse(adata.X):
        raise RuntimeError("Expected a sparse count matrix")
    matrix = adata.X.tocsr()
    if matrix.nnz < adata.n_obs or np.any(matrix.data < 0) or not np.allclose(matrix.data, np.rint(matrix.data)):
        raise RuntimeError("Count matrix integrity check failed")
    if not adata.obs_names.is_unique or not adata.var_names.is_unique or "spatial" not in adata.obsm:
        raise RuntimeError("Identifiers or spatial coordinates failed validation")

    csc = matrix.tocsc()
    detected = np.diff(csc.indptr).astype(int)
    totals = np.asarray(csc.sum(axis=0)).ravel().astype(float)
    sum_squares = np.asarray(csc.power(2).sum(axis=0)).ravel().astype(float)
    variance_numerator = sum_squares - totals**2 / adata.n_obs
    minimum_detected = max(20, math.ceil(0.01 * adata.n_obs))
    keep = (detected >= minimum_detected) & (totals >= 50) & (variance_numerator > 0)
    python_genes = adata.var_names[keep].tolist()

    spark = pd.read_csv(spark_path, sep="\t")
    if spark["gene"].duplicated().any() or set(spark["gene"]) != set(python_genes):
        raise RuntimeError("Python and R gene filters disagree")
    spark = spark.set_index("gene").loc[python_genes].reset_index()
    filtered = adata[:, python_genes].copy()
    filtered.X = filtered.X.tocsr()
    filtered.obs["total_counts_full"] = np.asarray(matrix.sum(axis=1)).ravel()

    memory_records.append(run_gate(f"{sample.lower()}_hotspot_autocorrelation"))
    hs = hotspot.Hotspot(
        filtered,
        model="danb",
        latent_obsm_key="spatial",
        umi_counts_obs_key="total_counts_full",
    )
    hs.create_knn_graph(weighted_graph=False, n_neighbors=6, approx_neighbors=False)
    autocorr = hs.compute_autocorrelations(jobs=4).rename(columns={"C": "hotspot_C", "Z": "hotspot_Z", "Pval": "hotspot_p", "FDR": "hotspot_fdr"})
    if autocorr.shape != (len(python_genes), 4) or not np.isfinite(autocorr.to_numpy(float)).all():
        raise RuntimeError(f"Hotspot autocorrelation output invalid: {autocorr.shape}")

    gene_stats = spark.merge(autocorr.reset_index(names="gene"), on="gene", how="inner", validate="one_to_one")
    gene_filter = pd.DataFrame({"gene": adata.var_names, "detected_spots": detected, "detection_fraction": detected / adata.n_obs, "total_counts": totals})
    gene_stats = gene_stats.merge(gene_filter, on="gene", how="left", validate="one_to_one")
    spark_strict = gene_stats["sparkx_by_fdr"] < 0.05
    hotspot_strict = (gene_stats["hotspot_fdr"] < 0.05) & (gene_stats["hotspot_C"] > 0)
    gene_stats["support_class"] = np.select(
        [spark_strict & hotspot_strict, spark_strict, hotspot_strict],
        ["both_strict", "sparkx_only", "hotspot_only"],
        default="not_strict",
    )
    gene_stats["sparkx_rank"] = gene_stats["sparkx_combined_p"].rank(method="min")
    gene_stats["hotspot_rank"] = gene_stats["hotspot_Z"].rank(method="min", ascending=False)
    gene_stats["consensus_rank"] = (gene_stats["sparkx_rank"] + gene_stats["hotspot_rank"]) / 2
    gene_stats = gene_stats.sort_values(["consensus_rank", "gene"]).reset_index(drop=True)

    hotspot_positive = gene_stats.loc[
        (gene_stats["hotspot_fdr"] < 0.05) & (gene_stats["hotspot_C"] > 0)
    ].sort_values(["hotspot_fdr", "hotspot_Z"], ascending=[True, False])
    if len(hotspot_positive) >= 50:
        module_genes = hotspot_positive.head(500)["gene"].tolist()
        module_selection = "top up to 500 positive Hotspot genes with FDR<0.05"
    else:
        module_genes = gene_stats.sort_values("hotspot_Z", ascending=False).head(min(200, len(gene_stats)))["gene"].tolist()
        module_selection = "prespecified fallback: top up to 200 Hotspot Z genes because fewer than 50 passed FDR<0.05"
    if len(module_genes) < 20:
        raise RuntimeError("Too few genes for Hotspot local-correlation modules")

    memory_records.append(run_gate(f"{sample.lower()}_hotspot_local_correlations"))
    local_z = hs.compute_local_correlations(module_genes, jobs=4)
    if local_z.shape != (len(module_genes), len(module_genes)) or not np.isfinite(local_z.to_numpy(float)).all():
        raise RuntimeError("Hotspot local-correlation matrix failed validation")
    modules = hs.create_modules(min_gene_threshold=15, core_only=True, fdr_threshold=0.05)
    if not modules.index.is_unique or set(modules.index) != set(module_genes):
        raise RuntimeError("Hotspot module assignment failed validation")
    module_scores = hs.calculate_module_scores()
    if not module_scores.index.equals(filtered.obs_names) or not np.isfinite(module_scores.to_numpy(float)).all():
        raise RuntimeError("Hotspot module scores failed validation")

    membership = gene_stats.set_index("gene").loc[module_genes].copy()
    membership["module"] = modules.loc[membership.index].astype(int)
    membership = membership.reset_index().sort_values(["module", "consensus_rank", "gene"])
    membership.to_csv(output_dir / f"{sample}_hotspot_module_membership.tsv", sep="\t", index=False)
    gene_stats.to_csv(output_dir / f"{sample}_sparkx_hotspot_gene_statistics.tsv", sep="\t", index=False)
    local_z.to_csv(output_dir / f"{sample}_hotspot_local_correlation_z.tsv", sep="\t")

    boundary = pd.read_csv(boundary_path, sep="\t").set_index("barcode")
    if set(boundary.index) != set(filtered.obs_names):
        raise RuntimeError("BayesSpace boundary and H5AD barcodes disagree")
    boundary = boundary.loc[filtered.obs_names]
    score_table = module_scores.copy()
    score_table.columns = [f"module_{int(column)}" for column in score_table.columns]
    score_table.insert(0, "barcode", score_table.index)
    score_table["bayesspace_cluster"] = boundary["bayesspace_cluster"].to_numpy(int)
    tumor_indicator = as_bool(boundary["computational_tumor"]).to_numpy()
    score_table["computational_tumor"] = tumor_indicator
    score_table["signed_distance_um"] = boundary["signed_distance_um"].to_numpy(float)
    score_table.to_csv(output_dir / f"{sample}_hotspot_module_scores.tsv", sep="\t", index=False)

    module_rows = []
    assigned_modules = sorted(int(value) for value in modules.unique() if int(value) != -1)
    for module in assigned_modules:
        module_members = membership.loc[membership["module"].eq(module)].sort_values("consensus_rank")
        score_column = module
        score_values = module_scores[score_column].to_numpy(float)
        rho = stats.spearmanr(score_values, boundary["signed_distance_um"].to_numpy(float)).statistic
        tumor = tumor_indicator
        overlaps = []
        member_set = set(module_members["gene"])
        for name, genes in PROGRAMS.items():
            shared = sorted(member_set & genes)
            if shared:
                overlaps.append(f"{name}:{','.join(shared)}")
        module_rows.append({
            "module": module,
            "n_genes": len(module_members),
            "n_both_strict": int(module_members["support_class"].eq("both_strict").sum()),
            "n_sparkx_strict": int((module_members["sparkx_by_fdr"] < 0.05).sum()),
            "n_hotspot_strict": int(((module_members["hotspot_fdr"] < 0.05) & (module_members["hotspot_C"] > 0)).sum()),
            "top_genes": ";".join(module_members.head(20)["gene"]),
            "candidate_program_overlaps": " | ".join(overlaps),
            "distance_spearman_rho_descriptive": float(rho),
            "tumor_median_score": float(np.median(score_values[tumor])),
            "nontumor_median_score": float(np.median(score_values[~tumor])),
            "tumor_minus_nontumor_median": float(np.median(score_values[tumor]) - np.median(score_values[~tumor])),
        })
    module_summary = pd.DataFrame(module_rows)
    module_summary.to_csv(output_dir / f"{sample}_hotspot_module_summary.tsv", sep="\t", index=False)

    memory_records.append(run_gate(f"{sample.lower()}_hotspot_plots"))
    plot_gene_concordance(gene_stats, sample, output_dir)
    if assigned_modules:
        hs.plot_local_correlations(yticklabels=False)
        correlation_figure = plt.gcf()
        correlation_figure.set_size_inches(12, 11)
        correlation_figure.suptitle(f"{sample} Hotspot local-correlation modules", y=1.01, fontsize=15)
        correlation_figure.savefig(output_dir / f"{sample}_hotspot_local_correlation_modules.png", dpi=300, bbox_inches="tight")
        correlation_figure.savefig(output_dir / f"{sample}_hotspot_local_correlation_modules.pdf", bbox_inches="tight")
        plt.close(correlation_figure)
    else:
        fig, ax = plt.subplots(figsize=(7, 4))
        ax.text(0.5, 0.5, "No Hotspot module passed the prespecified threshold", ha="center", va="center")
        ax.axis("off")
        fig.savefig(output_dir / f"{sample}_hotspot_local_correlation_modules.png", dpi=300, bbox_inches="tight")
        fig.savefig(output_dir / f"{sample}_hotspot_local_correlation_modules.pdf", bbox_inches="tight")
        plt.close(fig)
    module_sizes = modules[modules != -1].value_counts().astype(int)
    plot_module_maps(adata, module_scores, module_sizes, sample, output_dir)
    plot_module_summary(module_summary, sample, output_dir)

    run_metadata = {
        "status": "completed",
        "sample": sample,
        "input_h5ad": str(input_h5ad),
        "input_sha256": sha256(input_h5ad),
        "boundary_table": str(boundary_path),
        "boundary_sha256": sha256(boundary_path),
        "n_spots": adata.n_obs,
        "n_genes_input": adata.n_vars,
        "n_genes_tested": len(gene_stats),
        "n_sparkx_by_fdr_0_05": int((gene_stats["sparkx_by_fdr"] < 0.05).sum()),
        "n_hotspot_positive_fdr_0_05": int(((gene_stats["hotspot_fdr"] < 0.05) & (gene_stats["hotspot_C"] > 0)).sum()),
        "n_both_strict": int(gene_stats["support_class"].eq("both_strict").sum()),
        "n_module_input_genes": len(module_genes),
        "module_gene_selection": module_selection,
        "n_modules": len(assigned_modules),
        "n_assigned_module_genes": int((modules != -1).sum()),
        "parameters": {
            "gene_detection_minimum": minimum_detected,
            "gene_total_count_minimum": 50,
            "hotspot_model": "danb",
            "knn_neighbors": 6,
            "weighted_graph": False,
            "local_correlation_max_genes": 500,
            "module_min_gene_threshold": 15,
            "module_core_only": True,
            "module_fdr_threshold": 0.05,
            "jobs": 4,
            "seed": seed,
        },
        "packages": {
            "Python": sys.version.split()[0],
            "anndata": metadata.version("anndata"),
            "hotspotsc": metadata.version("hotspotsc"),
            "numpy": metadata.version("numpy"),
            "scipy": metadata.version("scipy"),
        },
        "memory_gate_records": memory_records,
        "inference_boundary": "SPARK-X and Hotspot significance is section-level spatial structure, not patient-level replication. Module-distance correlations are descriptive because spots are not independent and the boundary is computational.",
        "signed_distance_convention": "computational tumor-like side <=0; non-tumor-like side >=0",
    }
    write_json(output_dir / "hotspot_run_metadata.json", run_metadata)
    write_status(output_dir, sample, "completed", f"SPARK-X/Hotspot passed; modules={len(assigned_modules)}")
    print(
        f"SPARKX_HOTSPOT_{sample}_PASS genes={len(gene_stats)} "
        f"both={run_metadata['n_both_strict']} modules={len(assigned_modules)}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    if len(sys.argv) != 6:
        raise SystemExit("Usage: run_tjrcc_hotspot_spatial_programs_v1.py SAMPLE INPUT_H5AD OUTPUT_DIR SEED BOUNDARY_TSV")
    sample_arg = sys.argv[1]
    input_arg = Path(sys.argv[2]).resolve()
    output_arg = Path(sys.argv[3]).resolve()
    seed_arg = int(sys.argv[4])
    boundary_arg = Path(sys.argv[5]).resolve()
    try:
        raise SystemExit(main(sample_arg, input_arg, output_arg, seed_arg, boundary_arg))
    except Exception as exc:
        output_arg.mkdir(parents=True, exist_ok=True)
        write_status(output_arg, sample_arg, "failed", f"{type(exc).__name__}: {exc}")
        (output_arg / "hotspot_failure_traceback.txt").write_text(traceback.format_exc(), encoding="utf-8")
        print(f"SPARKX_HOTSPOT_{sample_arg}_FAILED: {type(exc).__name__}: {exc}", file=sys.stderr)
        raise
