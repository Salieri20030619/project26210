#!/usr/bin/env python3

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy import sparse


PROJECT_ROOT = Path("/mnt/BioDisk/ludingyu/project26210")
MEMORY_GATE = PROJECT_ROOT / "check_memory_gate.sh"
SPOT_DIAMETER_UM = 65.0


PROGRAMS = {
    "Tumor-epithelial": ["CA9", "EPCAM", "KRT8", "KRT18", "KRT19"],
    "POSTN-ECM": ["POSTN", "COL1A1", "COL1A2", "COL3A1", "COL4A1", "FN1", "SPARC"],
    "APOE-TAM": ["APOE", "C1QA", "C1QB", "C1QC", "TREM2", "LPL"],
    "CD8-effector": ["CD8A", "CD8B", "NKG7", "CCL5", "PRF1", "GZMB"],
    "TGFb-response": ["TGFB1", "TGFBR1", "TGFBR2", "SMAD2", "SMAD3", "SERPINE1", "CTGF"],
}


def memory_gate(stage: str) -> None:
    import subprocess

    print(f"MEMORY_GATE stage={stage}", flush=True)
    subprocess.run([str(MEMORY_GATE)], check=True)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def zscore(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    sd = np.nanstd(values, ddof=0)
    if not np.isfinite(sd) or sd == 0:
        return np.zeros_like(values)
    return (values - np.nanmean(values)) / sd


def program_score(adata: ad.AnnData, genes: list[str]) -> tuple[np.ndarray, list[str]]:
    present = [gene for gene in genes if gene in adata.var_names]
    if not present:
        return np.full(adata.n_obs, np.nan), []
    matrix = adata[:, present].X
    if sparse.issparse(matrix):
        matrix = matrix.toarray()
    matrix = np.asarray(matrix, dtype=float)
    totals = np.asarray(adata.X.sum(axis=1)).ravel().astype(float)
    totals[totals <= 0] = np.nan
    log_cpm = np.log1p(matrix / totals[:, None] * 1e6)
    standardized = np.column_stack([zscore(log_cpm[:, index]) for index in range(log_cpm.shape[1])])
    return np.nanmean(standardized, axis=1), present


def build_hex_neighbors(frame: pd.DataFrame) -> list[list[int]]:
    lookup = {
        (int(row.array_col), int(row.array_row)): index
        for index, row in frame.reset_index(drop=True).iterrows()
    }
    offsets = [(-2, 0), (2, 0), (-1, -1), (1, -1), (-1, 1), (1, 1)]
    result: list[list[int]] = []
    for row in frame.itertuples(index=False):
        neighbors = []
        for delta_col, delta_row in offsets:
            candidate = (int(row.array_col) + delta_col, int(row.array_row) + delta_row)
            if candidate in lookup:
                neighbors.append(lookup[candidate])
        result.append(neighbors)
    return result


def main(sample: str, analyzed_h5ad: Path, output_dir: Path) -> None:
    if not sample.replace("_", "").isalnum():
        raise ValueError("Invalid sample ID")
    assignments_path = output_dir / f"{sample}_bayesspace_assignments.tsv"
    output_dir.mkdir(parents=True, exist_ok=True)
    memory_gate("boundary_postprocess_read")
    assignments = pd.read_csv(assignments_path, sep="\t")
    if assignments.empty or assignments["barcode"].duplicated().any():
        raise RuntimeError("BayesSpace assignments failed row-count or uniqueness validation")

    adata = ad.read_h5ad(analyzed_h5ad)
    if adata.n_obs != assignments.shape[0] or not adata.obs_names.is_unique:
        raise RuntimeError("Unexpected analyzed H5AD dimensions or duplicate barcodes")
    missing = set(assignments["barcode"]) - set(adata.obs_names)
    if missing:
        raise RuntimeError(f"{len(missing)} BayesSpace barcodes are absent from analyzed H5AD")
    adata = adata[assignments["barcode"].to_numpy()].copy()

    abundance = adata.obsm["q05_cell_abundance_w_sf"]
    if not isinstance(abundance, pd.DataFrame):
        abundance = pd.DataFrame(abundance, index=adata.obs_names, columns=adata.uns["mod"]["factor_names"])
    rename = {
        column: str(column)
        .replace("q05cell_abundance_w_sf_means_per_cluster_mu_fg_", "")
        .replace("q05cell_abundance_w_sf_", "")
        for column in abundance.columns
    }
    abundance = abundance.rename(columns=rename)
    required_abundance = ["TYMP-High Malignant", "TYMP-Low Malignant", "APOE+C1QC+ TAM", "T cell"]
    if any(column not in abundance.columns for column in required_abundance):
        raise RuntimeError("Required cell2location columns are missing")

    assignments["malignant_q05"] = (
        abundance["TYMP-High Malignant"].to_numpy()
        + abundance["TYMP-Low Malignant"].to_numpy()
    )
    total_abundance = abundance.sum(axis=1).to_numpy()
    assignments["malignant_fraction"] = np.divide(
        assignments["malignant_q05"].to_numpy(),
        total_abundance,
        out=np.zeros_like(total_abundance, dtype=float),
        where=total_abundance > 0,
    )
    assignments["APOE_TAM_q05"] = abundance["APOE+C1QC+ TAM"].to_numpy()
    assignments["T_cell_q05"] = abundance["T cell"].to_numpy()

    program_genes: dict[str, list[str]] = {}
    for name, genes in PROGRAMS.items():
        score, present = program_score(adata, genes)
        assignments[name] = score
        program_genes[name] = present

    domain_summary = assignments.groupby("bayesspace_cluster", as_index=False).agg(
        n_spots=("barcode", "size"),
        malignant_q05_median=("malignant_q05", "median"),
        malignant_fraction_median=("malignant_fraction", "median"),
        tumor_epithelial_median=("Tumor-epithelial", "median"),
        postn_ecm_median=("POSTN-ECM", "median"),
        apoe_tam_program_median=("APOE-TAM", "median"),
        apoe_tam_q05_median=("APOE_TAM_q05", "median"),
        cd8_effector_median=("CD8-effector", "median"),
        t_cell_q05_median=("T_cell_q05", "median"),
    )
    domain_summary["tumor_composite"] = (
        zscore(domain_summary["malignant_fraction_median"].to_numpy())
        + zscore(domain_summary["tumor_epithelial_median"].to_numpy())
    ) / 2
    median_malignant = float(domain_summary["malignant_fraction_median"].median())
    domain_summary["tumor_like"] = (
        (domain_summary["tumor_composite"] > 0)
        & (domain_summary["malignant_fraction_median"] >= median_malignant)
    )
    if not domain_summary["tumor_like"].any():
        domain_summary.loc[domain_summary["tumor_composite"].idxmax(), "tumor_like"] = True
    if domain_summary["tumor_like"].all():
        domain_summary.loc[domain_summary["tumor_composite"].idxmin(), "tumor_like"] = False

    tumor_domains = set(domain_summary.loc[domain_summary["tumor_like"], "bayesspace_cluster"].astype(int))
    assignments["computational_tumor"] = assignments["bayesspace_cluster"].astype(int).isin(tumor_domains)
    neighbors = build_hex_neighbors(assignments)
    labels = assignments["computational_tumor"].to_numpy(dtype=bool)
    boundary = np.array(
        [any(labels[index] != labels[other] for other in adjacent) for index, adjacent in enumerate(neighbors)],
        dtype=bool,
    )
    if boundary.sum() == 0:
        raise RuntimeError("No computational tumor/non-tumor boundary was found")
    assignments["computational_boundary_spot"] = boundary

    spatial = assignments[["spatial_x_fullres", "spatial_y_fullres"]].to_numpy(dtype=float)
    tree = cKDTree(spatial[boundary])
    distance_pixels, _ = tree.query(spatial, k=1)
    spot_diameter_pixels = float(adata.uns["spatial"][sample]["scalefactors"]["spot_diameter_fullres"])
    microns_per_pixel = SPOT_DIAMETER_UM / spot_diameter_pixels
    signed_distance_um = distance_pixels * microns_per_pixel
    signed_distance_um[labels] *= -1
    assignments["signed_distance_um"] = signed_distance_um
    ordered_zones = ["Tumor core", "Inner interface", "Outer interface", "Peripheral 400-800um", "Distal >800um"]
    zone = np.select(
        [
            labels & (distance_pixels * microns_per_pixel > 400),
            labels & (distance_pixels * microns_per_pixel <= 400),
            (~labels) & (distance_pixels * microns_per_pixel <= 400),
            (~labels) & (distance_pixels * microns_per_pixel <= 800),
        ],
        ordered_zones[:4],
        default=ordered_zones[4],
    )
    assignments["distance_zone"] = pd.Categorical(zone, categories=ordered_zones, ordered=True)

    profile_columns = [
        "malignant_fraction", "Tumor-epithelial", "POSTN-ECM", "APOE-TAM", "CD8-effector",
        "TGFb-response", "APOE_TAM_q05", "T_cell_q05",
    ]
    profile = assignments.groupby("distance_zone", observed=False)[profile_columns].agg(["median", "mean", "count"])
    profile.columns = ["__".join(column) for column in profile.columns]
    profile = profile.reset_index()

    assignments.to_csv(output_dir / f"{sample}_computational_boundary_spots.tsv", sep="\t", index=False)
    domain_summary.to_csv(output_dir / f"{sample}_domain_summary.tsv", sep="\t", index=False)
    profile.to_csv(output_dir / f"{sample}_boundary_distance_profiles.tsv", sep="\t", index=False)

    memory_gate("boundary_postprocess_plot")
    hires = adata.uns["spatial"][sample]["images"]["hires"]
    hires_scale = float(adata.uns["spatial"][sample]["scalefactors"]["tissue_hires_scalef"])
    plot_x = assignments["spatial_x_fullres"].to_numpy() * hires_scale
    plot_y = assignments["spatial_y_fullres"].to_numpy() * hires_scale

    fig, axes = plt.subplots(1, 3, figsize=(18, 6), constrained_layout=True)
    for axis in axes:
        axis.imshow(hires)
        axis.set_axis_off()
    scatter = axes[0].scatter(plot_x, plot_y, c=assignments["bayesspace_cluster"], s=7, cmap="tab10", alpha=0.82)
    axes[0].set_title("BayesSpace domains")
    fig.colorbar(scatter, ax=axes[0], fraction=0.035, label="Domain")
    axes[1].scatter(plot_x[~labels], plot_y[~labels], c="#4C78A8", s=7, alpha=0.65, label="Non-tumor-like")
    axes[1].scatter(plot_x[labels], plot_y[labels], c="#E45756", s=7, alpha=0.65, label="Tumor-like")
    axes[1].scatter(plot_x[boundary], plot_y[boundary], facecolors="none", edgecolors="#111111", s=16, linewidths=0.35)
    axes[1].set_title("Computational tumor domains and boundary")
    axes[1].legend(frameon=False, loc="upper right", markerscale=2)
    vmax = np.nanpercentile(np.abs(signed_distance_um), 95)
    scatter = axes[2].scatter(plot_x, plot_y, c=signed_distance_um, s=7, cmap="coolwarm", vmin=-vmax, vmax=vmax, alpha=0.85)
    axes[2].set_title("Signed distance to computational boundary")
    fig.colorbar(scatter, ax=axes[2], fraction=0.035, label="Signed distance (um)")
    fig.savefig(output_dir / f"{sample}_bayesspace_boundary_overview.png", dpi=300, bbox_inches="tight")
    fig.savefig(output_dir / f"{sample}_bayesspace_boundary_overview.pdf", bbox_inches="tight")
    plt.close(fig)

    median_profile = assignments.groupby("distance_zone", observed=False)[profile_columns].median().reindex(ordered_zones)
    fig, axes = plt.subplots(2, 2, figsize=(12, 9), constrained_layout=True)
    plot_sets = [
        ["POSTN-ECM", "APOE-TAM", "CD8-effector", "TGFb-response"],
        ["malignant_fraction"],
        ["APOE_TAM_q05", "T_cell_q05"],
        ["Tumor-epithelial", "POSTN-ECM"],
    ]
    for axis, columns in zip(axes.flat, plot_sets):
        for column in columns:
            if column in median_profile:
                axis.plot(range(len(ordered_zones)), median_profile[column], marker="o", label=column)
        axis.set_xticks(range(len(ordered_zones)))
        axis.set_xticklabels(["Core", "Inner", "Outer", "400-800", ">800"], rotation=25, ha="right")
        axis.grid(alpha=0.2)
        axis.legend(frameon=False, fontsize=8)
    axes[0, 0].set_title("Program-score gradients")
    axes[0, 1].set_title("Malignant fraction")
    axes[1, 0].set_title("Existing cell2location proxies")
    axes[1, 1].set_title("Tumor and ECM programs")
    fig.suptitle(f"{sample} computational-boundary profiles", fontsize=14)
    fig.savefig(output_dir / f"{sample}_boundary_profile_summary.png", dpi=300, bbox_inches="tight")
    fig.savefig(output_dir / f"{sample}_boundary_profile_summary.pdf", bbox_inches="tight")
    plt.close(fig)

    metadata = {
        "status": "completed",
        "sample": sample,
        "completed_at": pd.Timestamp.now(tz="Asia/Shanghai").isoformat(),
        "assignments_sha256": sha256(assignments_path),
        "analyzed_h5ad": str(analyzed_h5ad),
        "analyzed_h5ad_sha256": sha256(analyzed_h5ad),
        "n_spots": int(assignments.shape[0]),
        "n_domains": int(assignments["bayesspace_cluster"].nunique()),
        "tumor_like_domains": sorted(int(value) for value in tumor_domains),
        "n_boundary_spots": int(boundary.sum()),
        "spot_diameter_fullres_pixels": spot_diameter_pixels,
        "microns_per_pixel": microns_per_pixel,
        "program_genes_present": program_genes,
        "boundary_definition": "BayesSpace domain composite of malignant fraction and tumor-epithelial program; not pathologist annotation",
        "inference_status": "descriptive computational-boundary pilot; spot-level zones are not biological replicates",
    }
    with (output_dir / "boundary_run_metadata.json").open("w", encoding="utf-8") as handle:
        json.dump(metadata, handle, ensure_ascii=False, indent=2)
    with (output_dir / "status.json").open("w", encoding="utf-8") as handle:
        json.dump({"status": "completed", "sample": sample, "completed_at": metadata["completed_at"]}, handle, indent=2)

    print(
        f"TJRCC_BAYESSPACE_BOUNDARY_PASS sample={sample} domains={metadata['n_domains']} "
        f"tumor_domains={metadata['tumor_like_domains']} boundary_spots={metadata['n_boundary_spots']}",
        flush=True,
    )


if __name__ == "__main__":
    if len(sys.argv) != 4:
        raise SystemExit("Usage: postprocess_tjrcc_bayesspace_boundary_v1.py SAMPLE INPUT_H5AD OUTPUT_DIR")
    main(sys.argv[1], Path(sys.argv[2]), Path(sys.argv[3]))
