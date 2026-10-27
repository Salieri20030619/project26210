#!/usr/bin/env python3

from __future__ import annotations

import hashlib
from importlib import metadata as importlib_metadata
import json
import math
import os
import subprocess
import sys
from pathlib import Path

os.environ.setdefault("POT_BACKEND_DISABLE_JAX", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import anndata as ad
import commot as ct
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse


PROJECT_ROOT = Path("/mnt/BioDisk/ludingyu/project26210")
RESULT_ROOT = PROJECT_ROOT / "results/part4_spatial"
MEMORY_GATE = PROJECT_ROOT / "check_memory_gate.sh"
SPOT_DIAMETER_UM = 65.0
COT_NITERMAX = 10000
SECRETED_THRESHOLDS = (200, 500)
CONTACT_THRESHOLD = 110
EXPECTED_CELLTYPES = {
    "APOE+C1QC+ TAM", "B cell", "Dendritic cell", "Endothelial", "Epithelial-Normal",
    "Mast cell", "Monocyte", "NK cell", "Other Macrophage", "SMC/Pericyte",
    "SPP1+ TAM", "T cell", "TYMP-High Malignant", "TYMP-Low Malignant",
}

CANDIDATES = [
    dict(axis_id="TGFB1_TGFBR1_TGFBR2", ligand="TGFB1", receptor="TGFBR1_TGFBR2", pathway="TGFb",
         signaling_type="Secreted Signaling", source_database="CellChat", source_proxy="TYMP-High Malignant", target_proxy="SMC/Pericyte"),
    dict(axis_id="POSTN_ITGAV_ITGB3", ligand="POSTN", receptor="ITGAV_ITGB3", pathway="PERIOSTIN",
         signaling_type="Secreted Signaling", source_database="CellChat", source_proxy="SMC/Pericyte", target_proxy="APOE+C1QC+ TAM"),
    dict(axis_id="POSTN_ITGAV_ITGB5", ligand="POSTN", receptor="ITGAV_ITGB5", pathway="PERIOSTIN",
         signaling_type="Secreted Signaling", source_database="CellChat", source_proxy="SMC/Pericyte", target_proxy="APOE+C1QC+ TAM"),
    dict(axis_id="CCL2_CCR2", ligand="CCL2", receptor="CCR2", pathway="CCL2",
         signaling_type="Secreted Signaling", source_database="CellChat", source_proxy="SMC/Pericyte", target_proxy="Monocyte"),
    dict(axis_id="CSF1_CSF1R", ligand="CSF1", receptor="CSF1R", pathway="CSF1",
         signaling_type="Secreted Signaling", source_database="CellChat", source_proxy="SMC/Pericyte", target_proxy="APOE+C1QC+ TAM"),
    dict(axis_id="SPP1_CD44", ligand="SPP1", receptor="CD44", pathway="SPP1",
         signaling_type="Secreted Signaling", source_database="CellChat", source_proxy="SPP1+ TAM", target_proxy="T cell"),
    dict(axis_id="SPP1_ITGAV_ITGB3", ligand="SPP1", receptor="ITGAV_ITGB3", pathway="SPP1",
         signaling_type="Secreted Signaling", source_database="CellChat", source_proxy="SPP1+ TAM", target_proxy="Endothelial"),
    dict(axis_id="MIF_CD74_CXCR4", ligand="MIF", receptor="CD74_CXCR4", pathway="MIF",
         signaling_type="Secreted Signaling", source_database="CellChat", source_proxy="TYMP-High Malignant", target_proxy="Dendritic cell"),
    dict(axis_id="MIF_CD74_CD44", ligand="MIF", receptor="CD74_CD44", pathway="MIF",
         signaling_type="Secreted Signaling", source_database="CellChat", source_proxy="TYMP-High Malignant", target_proxy="Dendritic cell"),
    dict(axis_id="VEGFA_KDR", ligand="VEGFA", receptor="KDR", pathway="VEGFA",
         signaling_type="Secreted Signaling", source_database="CellChat", source_proxy="APOE+C1QC+ TAM", target_proxy="Endothelial"),
    dict(axis_id="APOE_TREM2_TYROBP", ligand="APOE", receptor="TREM2_TYROBP", pathway="APOE_TREM2",
         signaling_type="Secreted Signaling", source_database="CellPhoneDB_v4.0", source_proxy="APOE+C1QC+ TAM", target_proxy="APOE+C1QC+ TAM"),
    dict(axis_id="APP_CD74", ligand="APP", receptor="CD74", pathway="APP",
         signaling_type="Cell-Cell Contact", source_database="CellChat", source_proxy="TYMP-High Malignant", target_proxy="Dendritic cell"),
]


def memory_gate(stage: str) -> None:
    print(f"MEMORY_GATE stage={stage}", flush=True)
    subprocess.run([str(MEMORY_GATE)], check=True)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def validate_database_rows() -> None:
    cellchat = ct.pp.ligand_receptor_database(database="CellChat", species="human", signaling_type=None)
    cpdb = ct.pp.ligand_receptor_database(database="CellPhoneDB_v4.0", species="human", signaling_type=None)
    for candidate in CANDIDATES:
        database = cellchat if candidate["source_database"] == "CellChat" else cpdb
        mask = (database.iloc[:, 0].astype(str) == candidate["ligand"]) & (
            database.iloc[:, 1].astype(str) == candidate["receptor"]
        )
        if int(mask.sum()) != 1:
            raise RuntimeError(f"Database row not uniquely resolved: {candidate['axis_id']} count={int(mask.sum())}")
        observed_type = str(database.loc[mask].iloc[0, 3])
        if observed_type != candidate["signaling_type"]:
            raise RuntimeError(f"Signaling type mismatch for {candidate['axis_id']}: {observed_type}")


def components(candidate: dict[str, str]) -> list[str]:
    return candidate["ligand"].split("_") + candidate["receptor"].split("_")


def matrix_summary(matrix: sparse.spmatrix, coords: np.ndarray) -> dict[str, float | int]:
    coo = matrix.tocoo()
    total = float(coo.data.sum())
    diagonal_mask = coo.row == coo.col
    diagonal_score = float(coo.data[diagonal_mask].sum())
    off_mask = ~diagonal_mask
    off_score = float(coo.data[off_mask].sum())
    distances = np.sqrt(np.sum((coords[coo.row] - coords[coo.col]) ** 2, axis=1))
    if total > 0:
        mean_distance = float(np.average(distances, weights=coo.data))
        order = np.argsort(distances)
        ordered_dist = distances[order]
        ordered_weight = coo.data[order]
        cumulative = np.cumsum(ordered_weight) / total
        median_distance = float(ordered_dist[np.searchsorted(cumulative, 0.5)])
        p90_distance = float(ordered_dist[np.searchsorted(cumulative, 0.9)])
        share_110 = float(coo.data[distances <= 110].sum() / total)
        share_200 = float(coo.data[distances <= 200].sum() / total)
    else:
        mean_distance = median_distance = p90_distance = share_110 = share_200 = math.nan
    return {
        "total_score": total,
        "n_edges": int(coo.nnz),
        "diagonal_score": diagonal_score,
        "diagonal_share": diagonal_score / total if total > 0 else math.nan,
        "offdiagonal_score": off_score,
        "offdiagonal_edges": int(off_mask.sum()),
        "weighted_mean_distance_um": mean_distance,
        "weighted_median_distance_um": median_distance,
        "weighted_p90_distance_um": p90_distance,
        "score_share_within_110um": share_110,
        "score_share_within_200um": share_200,
        "n_sender_spots": int(np.count_nonzero(np.asarray(matrix.sum(axis=1)).ravel())),
        "n_receiver_spots": int(np.count_nonzero(np.asarray(matrix.sum(axis=0)).ravel())),
    }


def candidate_dataframe(candidates: list[dict[str, str]]) -> pd.DataFrame:
    return pd.DataFrame([[item["ligand"], item["receptor"], item["axis_id"]] for item in candidates])


def run_model(adata: ad.AnnData, candidates: list[dict[str, str]], database_name: str, threshold: int, stage: str) -> None:
    memory_gate(stage)
    ct.tl.spatial_communication(
        adata,
        database_name=database_name,
        df_ligrec=candidate_dataframe(candidates),
        pathway_sum=True,
        heteromeric=True,
        heteromeric_rule="min",
        dis_thr=float(threshold),
        cost_type="euc",
        cot_eps_p=0.1,
        cot_rho=10.0,
        cot_nitermax=COT_NITERMAX,
        cot_weights=(0.25, 0.25, 0.25, 0.25),
        smooth=False,
    )


def plot_spatial_maps(
    sample: str,
    output_dir: Path,
    model_adata: ad.AnnData,
    analyzed: ad.AnnData,
    fullres_coords: np.ndarray,
) -> None:
    selected_axes = ["TGFB1_TGFBR1_TGFBR2", "POSTN_ITGAV_ITGB5", "MIF_CD74_CXCR4", "APOE_TREM2_TYROBP"]
    candidate_lookup = {item["axis_id"]: item for item in CANDIDATES}
    hires = analyzed.uns["spatial"][sample]["images"]["hires"]
    hires_scale = float(analyzed.uns["spatial"][sample]["scalefactors"]["tissue_hires_scalef"])
    px = fullres_coords[:, 0] * hires_scale
    py = fullres_coords[:, 1] * hires_scale
    fig, axes = plt.subplots(4, 2, figsize=(12, 20), constrained_layout=True)
    for row, axis_id in enumerate(selected_axes):
        candidate = candidate_lookup[axis_id]
        key = f"commot-candidate200-{candidate['ligand']}-{candidate['receptor']}"
        matrix = model_adata.obsp[key]
        sender = np.asarray(matrix.sum(axis=1)).ravel()
        receiver = np.asarray(matrix.sum(axis=0)).ravel()
        for column, (role, values) in enumerate((("Sender", sender), ("Receiver", receiver))):
            axis = axes[row, column]
            axis.imshow(hires)
            vmax = float(np.quantile(np.log1p(values), 0.99))
            scatter = axis.scatter(px, py, c=np.log1p(values), s=6, cmap="magma", vmin=0, vmax=max(vmax, 1e-9), alpha=0.86)
            axis.set_title(f"{axis_id}: {role} score (200 um)", fontsize=10)
            axis.set_axis_off()
            fig.colorbar(scatter, ax=axis, fraction=0.035, pad=0.01, label="log1p COMMOT marginal")
    fig.suptitle(f"{sample}: selected targeted COMMOT sender/receiver maps", fontsize=14)
    fig.savefig(output_dir / f"{sample}_commot_selected_spatial_maps_v2.png", dpi=300, bbox_inches="tight")
    fig.savefig(output_dir / f"{sample}_commot_selected_spatial_maps_v2.pdf", bbox_inches="tight")
    plt.close(fig)


def plot_summary(sample: str, output_dir: Path, summary: pd.DataFrame, hypothesis: pd.DataFrame) -> None:
    scale_order = ["contact_110um", "secreted_200um", "secreted_500um"]
    axis_order = [item["axis_id"] for item in CANDIDATES]
    total = summary.pivot(index="axis_id", columns="scale", values="total_score").reindex(axis_order)
    fig, axis = plt.subplots(figsize=(10, 7), constrained_layout=True)
    values = np.log10(total.fillna(0).to_numpy() + 1e-8)
    image = axis.imshow(values, cmap="viridis", aspect="auto")
    axis.set_xticks(range(len(total.columns)), total.columns)
    axis.set_yticks(range(len(total.index)), total.index)
    axis.set_title(f"{sample}: targeted COMMOT total transport score")
    for row in range(values.shape[0]):
        for column in range(values.shape[1]):
            if np.isfinite(total.iloc[row, column]):
                axis.text(column, row, f"{total.iloc[row, column]:.1f}", ha="center", va="center", fontsize=7,
                          color="white" if values[row, column] < np.nanmedian(values) else "black")
    fig.colorbar(image, ax=axis, fraction=0.035, label="log10(total score)")
    fig.savefig(output_dir / f"{sample}_commot_total_score_v2.png", dpi=300, bbox_inches="tight")
    fig.savefig(output_dir / f"{sample}_commot_total_score_v2.pdf", bbox_inches="tight")
    plt.close(fig)

    enrich = hypothesis.pivot(index="axis_id", columns="scale", values="proxy_log2_enrichment").reindex(axis_order)
    enrich = enrich.reindex(columns=[item for item in scale_order if item in enrich.columns])
    fig, axis = plt.subplots(figsize=(9, 7), constrained_layout=True)
    cap = max(float(np.nanpercentile(np.abs(enrich.to_numpy()), 95)), 0.5)
    image = axis.imshow(enrich.to_numpy(), cmap="coolwarm", vmin=-cap, vmax=cap, aspect="auto")
    axis.set_xticks(range(len(enrich.columns)), enrich.columns)
    axis.set_yticks(range(len(enrich.index)), enrich.index)
    axis.set_title(f"{sample}: prespecified composition-weighted source-target proxies")
    for row in range(enrich.shape[0]):
        for column in range(enrich.shape[1]):
            value = enrich.iloc[row, column]
            if np.isfinite(value):
                axis.text(column, row, f"{value:.2f}", ha="center", va="center", fontsize=7)
    fig.colorbar(image, ax=axis, fraction=0.035, label="log2 enrichment over transport marginals")
    fig.text(0.5, -0.01, "Composition-weighted mixed-spot proxy; not cell-resolved communication.", ha="center", fontsize=9)
    fig.savefig(output_dir / f"{sample}_commot_hypothesis_proxy_v2.png", dpi=300, bbox_inches="tight")
    fig.savefig(output_dir / f"{sample}_commot_hypothesis_proxy_v2.pdf", bbox_inches="tight")
    plt.close(fig)


def main(sample: str, input_h5ad: Path, boundary_path: Path, output_dir: Path) -> None:
    if not sample.replace("_", "").isalnum():
        raise ValueError("Invalid sample ID")
    sample_lower = sample.lower()
    output_dir.mkdir(parents=True, exist_ok=True)
    status_path = output_dir / "status.json"
    if status_path.exists():
        prior = json.loads(status_path.read_text())
        if prior.get("status") == "complete":
            raise RuntimeError(f"Refusing to overwrite completed output: {output_dir}")
    status_path.write_text(json.dumps({"status": "running", "sample": sample}, indent=2) + "\n")

    memory_gate(f"commot_{sample_lower}_read")
    validate_database_rows()
    mapping_path = input_h5ad
    analyzed_path = input_h5ad
    mapping = ad.read_h5ad(mapping_path)
    analyzed = ad.read_h5ad(analyzed_path)
    boundary = pd.read_csv(boundary_path, sep="\t").set_index("barcode").loc[mapping.obs_names]
    if set(mapping.obs_names) != set(analyzed.obs_names) or not mapping.obs_names.is_unique:
        raise RuntimeError("Mapping/analyzed barcode mismatch or duplicates")
    analyzed = analyzed[mapping.obs_names].copy()

    all_genes = sorted({gene for candidate in CANDIDATES for gene in components(candidate)})
    missing = sorted(set(all_genes) - set(mapping.var_names))
    present_genes = [gene for gene in all_genes if gene in mapping.var_names]
    if len(present_genes) < 2:
        raise RuntimeError(f"Too few candidate component genes are present: {present_genes}")
    raw_candidate = mapping[:, present_genes].X.tocsr().astype(float)
    library_size = np.asarray(mapping.X.sum(axis=1)).ravel().astype(float)
    if np.any(library_size <= 0) or np.any(raw_candidate.data < 0):
        raise RuntimeError("Invalid raw count matrix")
    detection = np.asarray((raw_candidate > 0).sum(axis=0)).ravel()
    detection_by_gene = {gene: 0 for gene in all_genes}
    detection_by_gene.update(dict(zip(present_genes, detection)))
    min_spots = max(20, int(math.ceil(mapping.n_obs * 0.01)))
    manifest_rows = []
    available_candidates = []
    for candidate in CANDIDATES:
        component_counts = {gene: int(detection_by_gene[gene]) for gene in components(candidate)}
        available = all(value >= min_spots for value in component_counts.values())
        row = candidate.copy()
        row.update(
            {
                "sample": sample,
                "available": available,
                "minimum_required_spots": min_spots,
                "component_detection_counts": ";".join(f"{gene}:{count}" for gene, count in component_counts.items()),
                "distance_design": "110um contact" if candidate["signaling_type"] == "Cell-Cell Contact" else "200um primary;500um sensitivity",
                "celltype_interpretation": "prespecified composition-weighted mixed-spot proxy; not cell-resolved attribution",
            }
        )
        manifest_rows.append(row)
        if available:
            available_candidates.append(candidate)
    unavailable = [item["axis_id"] for item in manifest_rows if not item["available"]]
    if len(available_candidates) < 8:
        raise RuntimeError(f"Fewer than 8 prespecified axes passed the detection gate: {unavailable}")
    pd.DataFrame(manifest_rows).to_csv(output_dir / f"{sample}_commot_candidate_manifest_v2.tsv", sep="\t", index=False)

    normalized = sparse.diags(1e4 / library_size).dot(raw_candidate).tocsr()
    normalized.data = np.log1p(normalized.data)
    model_adata = ad.AnnData(
        X=normalized,
        obs=pd.DataFrame({"library_size": library_size}, index=mapping.obs_names.copy()),
        var=pd.DataFrame(index=pd.Index(present_genes, name="gene")),
    )
    fullres_coords = np.asarray(mapping.obsm["spatial"], dtype=float)
    spot_diameter_pixels = float(analyzed.uns["spatial"][sample]["scalefactors"]["spot_diameter_fullres"])
    microns_per_pixel = SPOT_DIAMETER_UM / spot_diameter_pixels
    coords_um = fullres_coords * microns_per_pixel
    model_adata.obsm["spatial"] = coords_um

    secreted = [item for item in available_candidates if item["signaling_type"] == "Secreted Signaling"]
    contact = [item for item in available_candidates if item["signaling_type"] == "Cell-Cell Contact"]
    for threshold in SECRETED_THRESHOLDS:
        run_model(model_adata, secreted, f"candidate{threshold}", threshold,
                  f"commot_{sample_lower}_secreted_{threshold}um")
    run_model(model_adata, contact, "contact110", CONTACT_THRESHOLD, f"commot_{sample_lower}_contact_110um")

    abundance = analyzed.obsm["q05_cell_abundance_w_sf"]
    if not isinstance(abundance, pd.DataFrame):
        abundance = pd.DataFrame(abundance, index=analyzed.obs_names, columns=analyzed.uns["mod"]["factor_names"])
    abundance = abundance.rename(
        columns=lambda value: str(value)
        .replace("q05cell_abundance_w_sf_means_per_cluster_mu_fg_", "")
        .replace("q05cell_abundance_w_sf_", "")
    )
    abundance = abundance.loc[mapping.obs_names]
    abundance_values = abundance.to_numpy(dtype=float)
    celltypes = abundance.columns.tolist()
    if len(celltypes) != len(set(celltypes)) or set(celltypes) != EXPECTED_CELLTYPES:
        raise RuntimeError(f"Unexpected or duplicate abundance cell types: {celltypes}")
    if not np.all(np.isfinite(abundance_values)) or np.any(abundance_values < 0):
        raise RuntimeError("Cell-abundance matrix contains non-finite or negative values")
    abundance_row_sums = abundance_values.sum(axis=1, keepdims=True)
    if np.any(abundance_row_sums <= 0):
        raise RuntimeError("Cell-abundance matrix contains spots with non-positive total abundance")
    abundance_fraction = abundance_values / abundance_row_sums

    summary_rows = []
    celltype_rows = []
    boundary_rows = []
    spot_table = pd.DataFrame(
        {
            "barcode": mapping.obs_names,
            "x_um": coords_um[:, 0],
            "y_um": coords_um[:, 1],
            "computational_tumor": boundary["computational_tumor"].astype(bool).to_numpy(),
            "signed_distance_um": boundary["signed_distance_um"].to_numpy(dtype=float),
        }
    )
    candidate_lookup = {item["axis_id"]: item for item in CANDIDATES}
    scale_specs = [(f"candidate{threshold}", f"secreted_{threshold}um", secreted) for threshold in SECRETED_THRESHOLDS]
    scale_specs.append(("contact110", "contact_110um", contact))
    hypothesis_rows = []
    tumor_indicator = boundary["computational_tumor"].astype(bool).to_numpy()
    region_weights = np.column_stack([tumor_indicator.astype(float), (~tumor_indicator).astype(float)])
    region_names = ["computational_tumor_like", "computational_non_tumor_like"]
    for database_name, scale_label, candidates_at_scale in scale_specs:
        for candidate in candidates_at_scale:
            key = f"commot-{database_name}-{candidate['ligand']}-{candidate['receptor']}"
            matrix = model_adata.obsp[key].tocsr()
            row = {"sample": sample, "scale": scale_label, **candidate, **matrix_summary(matrix, coords_um)}
            summary_rows.append(row)
            sender = np.asarray(matrix.sum(axis=1)).ravel()
            receiver = np.asarray(matrix.sum(axis=0)).ravel()
            spot_table[f"sender__{scale_label}__{candidate['axis_id']}"] = sender
            spot_table[f"receiver__{scale_label}__{candidate['axis_id']}"] = receiver

            mass = abundance_fraction.T @ matrix.dot(abundance_fraction)
            mass_total = float(mass.sum())
            mass_fraction = mass / mass_total if mass_total > 0 else np.full_like(mass, np.nan)
            source_marginal = np.nansum(mass_fraction, axis=1)
            target_marginal = np.nansum(mass_fraction, axis=0)
            baseline = np.outer(source_marginal, target_marginal)
            enrichment = np.divide(mass_fraction, baseline, out=np.full_like(mass_fraction, np.nan), where=baseline > 0)
            flat_rank = pd.Series(mass_fraction.ravel()).rank(method="min", ascending=False).to_numpy().reshape(mass_fraction.shape)
            for source_idx, source_celltype in enumerate(celltypes):
                for target_idx, target_celltype in enumerate(celltypes):
                    celltype_rows.append(
                        {
                            "sample": sample,
                            "scale": scale_label,
                            "axis_id": candidate["axis_id"],
                            "ligand": candidate["ligand"],
                            "receptor": candidate["receptor"],
                            "source_celltype": source_celltype,
                            "target_celltype": target_celltype,
                            "weighted_mass": float(mass[source_idx, target_idx]),
                            "weighted_mass_fraction": float(mass_fraction[source_idx, target_idx]),
                            "marginal_independence_expected_fraction": float(baseline[source_idx, target_idx]),
                            "proxy_enrichment": float(enrichment[source_idx, target_idx]),
                            "proxy_log2_enrichment": float(np.log2(enrichment[source_idx, target_idx])) if enrichment[source_idx, target_idx] > 0 else math.nan,
                            "mass_rank_among_196": int(flat_rank[source_idx, target_idx]),
                            "mass_rank_fraction": float(flat_rank[source_idx, target_idx] / mass_fraction.size),
                            "interpretation_limit": "composition-weighted mixed-spot edge-mass proxy; not cell-resolved communication",
                        }
                    )
            source_idx = celltypes.index(candidate["source_proxy"])
            target_idx = celltypes.index(candidate["target_proxy"])
            hypothesis_rows.append(
                {
                    "sample": sample,
                    "scale": scale_label,
                    "axis_id": candidate["axis_id"],
                    "ligand": candidate["ligand"],
                    "receptor": candidate["receptor"],
                    "source_proxy": candidate["source_proxy"],
                    "target_proxy": candidate["target_proxy"],
                    "proxy_weighted_mass_fraction": float(mass_fraction[source_idx, target_idx]),
                    "proxy_expected_fraction": float(baseline[source_idx, target_idx]),
                    "proxy_enrichment": float(enrichment[source_idx, target_idx]),
                    "proxy_log2_enrichment": float(np.log2(enrichment[source_idx, target_idx])) if enrichment[source_idx, target_idx] > 0 else math.nan,
                    "proxy_mass_rank_among_196": int(flat_rank[source_idx, target_idx]),
                    "proxy_mass_rank_fraction": float(flat_rank[source_idx, target_idx] / mass_fraction.size),
                    "interpretation_limit": "prespecified composition-weighted mixed-spot proxy; not cell-resolved attribution or causality",
                }
            )

            region_mass = region_weights.T @ matrix.dot(region_weights)
            region_total = float(region_mass.sum())
            for source_idx, source_region in enumerate(region_names):
                for target_idx, target_region in enumerate(region_names):
                    boundary_rows.append(
                        {
                            "sample": sample,
                            "scale": scale_label,
                            "axis_id": candidate["axis_id"],
                            "source_region": source_region,
                            "target_region": target_region,
                            "edge_mass": float(region_mass[source_idx, target_idx]),
                            "edge_mass_fraction": float(region_mass[source_idx, target_idx] / region_total) if region_total > 0 else math.nan,
                            "interpretation_limit": "computational-domain summary; not pathology boundary or independent-spot inference",
                        }
                    )

    summary = pd.DataFrame(summary_rows)
    celltype_mass = pd.DataFrame(celltype_rows)
    hypothesis = pd.DataFrame(hypothesis_rows)
    boundary_mass = pd.DataFrame(boundary_rows)
    summary.to_csv(output_dir / f"{sample}_commot_transport_summary_v2.tsv", sep="\t", index=False)
    celltype_mass.to_csv(output_dir / f"{sample}_commot_celltype_edge_mass_proxy_v2.tsv", sep="\t", index=False)
    hypothesis.to_csv(output_dir / f"{sample}_commot_hypothesis_proxy_support_v2.tsv", sep="\t", index=False)
    boundary_mass.to_csv(output_dir / f"{sample}_commot_computational_domain_edge_mass_v2.tsv", sep="\t", index=False)
    spot_table.to_csv(output_dir / f"{sample}_commot_spot_sender_receiver_v2.tsv", sep="\t", index=False)
    model_adata.write_h5ad(output_dir / f"{sample}_commot_candidate_v2.h5ad", compression="gzip")

    plot_summary(sample, output_dir, summary, hypothesis)
    plot_spatial_maps(sample, output_dir, model_adata, analyzed, fullres_coords)
    run_metadata = {
        "status": "complete",
        "sample": sample,
        "n_spots": int(mapping.n_obs),
        "n_raw_genes": int(mapping.n_vars),
        "n_candidate_component_genes_prespecified": len(all_genes),
        "n_candidate_component_genes_present": len(present_genes),
        "missing_candidate_component_genes": missing,
        "n_candidate_axes_prespecified": len(CANDIDATES),
        "n_candidate_axes_modeled": len(available_candidates),
        "unavailable_axes": unavailable,
        "n_secreted_axes": len(secreted),
        "n_contact_axes": len(contact),
        "minimum_detection_spots": min_spots,
        "normalization": "library-size normalize to 1e4 then log1p; candidate genes extracted after full-library totals",
        "coordinate_units": "micrometers",
        "microns_per_pixel": microns_per_pixel,
        "secreted_distance_thresholds_um": list(SECRETED_THRESHOLDS),
        "contact_distance_threshold_um": CONTACT_THRESHOLD,
        "commot_version": importlib_metadata.version("commot"),
        "cot_nitermax": COT_NITERMAX,
        "cot_eps_p": 0.1,
        "cot_rho": 10.0,
        "cot_weights": [0.25, 0.25, 0.25, 0.25],
        "mapping_h5ad": str(mapping_path),
        "mapping_h5ad_sha256": sha256(mapping_path),
        "analyzed_h5ad": str(analyzed_path),
        "analyzed_h5ad_sha256": sha256(analyzed_path),
        "boundary_table": str(boundary_path),
        "boundary_table_sha256": sha256(boundary_path),
        "interpretation_limit": (
            "targeted collective optimal-transport scores and composition-weighted mixed-spot proxies; "
            "not direct ligand-receptor binding, cell-resolved attribution, causal direction, or patient-level replication"
        ),
    }
    (output_dir / "run_metadata_v2.json").write_text(json.dumps(run_metadata, indent=2) + "\n")
    status_path.write_text(json.dumps({"status": "complete", "sample": sample}, indent=2) + "\n")
    print(
        f"COMMOT_SAMPLE_ALL_PASS sample={sample} spots={mapping.n_obs} "
        f"axes_modeled={len(available_candidates)}/{len(CANDIDATES)} scales=3"
    )


if __name__ == "__main__":
    if len(sys.argv) != 5:
        raise SystemExit("Usage: run_tjrcc_commot_candidate_v1.py SAMPLE INPUT_H5AD BOUNDARY_TSV OUTPUT_DIR")
    sample_arg = sys.argv[1]
    input_arg = Path(sys.argv[2]).resolve()
    boundary_arg = Path(sys.argv[3]).resolve()
    output_arg = Path(sys.argv[4]).resolve()
    try:
        main(sample_arg, input_arg, boundary_arg, output_arg)
    except Exception as exc:
        output = output_arg if len(sys.argv) == 5 else RESULT_ROOT / "failed_tjrcc_commot_unknown"
        output.mkdir(parents=True, exist_ok=True)
        (output / "status.json").write_text(json.dumps({"status": "failed", "sample": sample_arg, "error": str(exc)}, indent=2) + "\n")
        raise
