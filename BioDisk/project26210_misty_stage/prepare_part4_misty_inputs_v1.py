#!/usr/bin/env python3

from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.spatial import cKDTree


PROJECT_ROOT = Path("/mnt/BioDisk/ludingyu/project26210")
RESULT_ROOT = PROJECT_ROOT / "results/part4_spatial"
OUTPUT_DIR = RESULT_ROOT / "advanced_patch_misty_inputs_v1"
MEMORY_GATE = PROJECT_ROOT / "check_memory_gate.sh"
SAMPLES = ("GLMF1", "GLMF2")
SPOT_DIAMETER_UM = 65.0
JUXTA_THRESHOLD_UM = 110.0

FIXED_PROGRAMS = (
    "Tumor-epithelial",
    "POSTN-ECM",
    "APOE-TAM",
    "CD8-effector",
    "TGFb-response",
)

CELL_IDS = {
    "APOE+C1QC+ TAM": "cell_APOE_C1QC_TAM",
    "B cell": "cell_B",
    "Dendritic cell": "cell_DC",
    "Endothelial": "cell_Endothelial",
    "Epithelial-Normal": "cell_Epithelial_Normal",
    "Mast cell": "cell_Mast",
    "Monocyte": "cell_Monocyte",
    "NK cell": "cell_NK",
    "Other Macrophage": "cell_Other_Macrophage",
    "SMC/Pericyte": "cell_SMC_Pericyte",
    "SPP1+ TAM": "cell_SPP1_TAM",
    "T cell": "cell_T",
    "TYMP-High Malignant": "cell_TYMP_High_Malignant",
    "TYMP-Low Malignant": "cell_TYMP_Low_Malignant",
}

PROGRAM_IDS = {
    "Tumor-epithelial": "program_Tumor_epithelial",
    "POSTN-ECM": "program_POSTN_ECM",
    "APOE-TAM": "program_APOE_TAM",
    "CD8-effector": "program_CD8_effector",
    "TGFb-response": "program_TGFb_response",
}


def memory_gate(stage: str) -> None:
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
    if not np.isfinite(values).all():
        raise RuntimeError("A feature contains non-finite values before standardization")
    sd = float(np.std(values, ddof=0))
    if not np.isfinite(sd) or sd <= 0:
        raise RuntimeError("A requested MISTy feature has zero or invalid variance")
    return (values - float(np.mean(values))) / sd


def gene_program_score(adata: ad.AnnData, genes: list[str]) -> tuple[np.ndarray, list[str]]:
    present = [gene for gene in genes if gene in adata.var_names]
    if not present:
        raise RuntimeError("No genes from the shared stromal program are present")
    matrix = adata[:, present].X
    if sparse.issparse(matrix):
        matrix = matrix.toarray()
    matrix = np.asarray(matrix, dtype=float)
    totals = np.asarray(adata.X.sum(axis=1)).ravel().astype(float)
    if np.any(totals <= 0):
        raise RuntimeError("Non-positive library size encountered")
    log_cpm = np.log1p(matrix / totals[:, None] * 1e6)
    standardized = np.column_stack([zscore(log_cpm[:, index]) for index in range(log_cpm.shape[1])])
    return zscore(np.mean(standardized, axis=1)), present


def load_shared_stromal_genes() -> list[str]:
    overlap_path = (
        RESULT_ROOT
        / "advanced_patch_spatial_programs_cross_sample_v2"
        / "cross_sample_module_overlap_v2.tsv"
    )
    overlap = pd.read_csv(overlap_path, sep="\t")
    row = overlap.loc[(overlap["GLMF1_module"] == 1) & (overlap["GLMF2_module"] == 1)]
    if row.shape[0] != 1 or int(row.iloc[0]["n_shared"]) != 37:
        raise RuntimeError("The expected replicated GLMF1-M1/GLMF2-M1 stromal overlap was not found")
    genes = str(row.iloc[0]["shared_genes"]).split(";")
    if len(genes) != 37 or len(set(genes)) != 37:
        raise RuntimeError("Unexpected shared stromal gene list")
    return genes


def build_sample(sample: str, stromal_genes: list[str]) -> dict[str, object]:
    sample_lower = sample.lower()
    h5ad_path = RESULT_ROOT / "analysis" / f"{sample}_analyzed.h5ad"
    boundary_path = (
        RESULT_ROOT
        / f"advanced_patch_bayesspace_{sample_lower}_v1"
        / f"{sample}_computational_boundary_spots.tsv"
    )
    adata = ad.read_h5ad(h5ad_path)
    boundary = pd.read_csv(boundary_path, sep="\t")
    if boundary["barcode"].duplicated().any() or not adata.obs_names.is_unique:
        raise RuntimeError(f"{sample}: duplicate barcodes")
    if set(boundary["barcode"]) != set(adata.obs_names):
        raise RuntimeError(f"{sample}: H5AD and boundary barcodes do not match")
    adata = adata[boundary["barcode"].to_numpy()].copy()

    abundance = adata.obsm["q05_cell_abundance_w_sf"]
    if not isinstance(abundance, pd.DataFrame):
        abundance = pd.DataFrame(
            abundance,
            index=adata.obs_names,
            columns=adata.uns["mod"]["factor_names"],
        )
    abundance = abundance.rename(columns=lambda value: str(value).split("_fg_", 1)[-1])
    if set(CELL_IDS) != set(abundance.columns):
        raise RuntimeError(f"{sample}: the common 14-cell abundance panel changed")
    if (abundance.to_numpy(dtype=float) < 0).any():
        raise RuntimeError(f"{sample}: negative q05 cell abundance")

    measurements = pd.DataFrame(index=boundary["barcode"].astype(str))
    for source_label, feature_id in CELL_IDS.items():
        measurements[feature_id] = zscore(np.log1p(abundance[source_label].to_numpy(dtype=float)))
    for source_label in FIXED_PROGRAMS:
        measurements[PROGRAM_IDS[source_label]] = zscore(boundary[source_label].to_numpy(dtype=float))
    stromal_score, present_stromal_genes = gene_program_score(adata, stromal_genes)
    if present_stromal_genes != stromal_genes:
        missing = sorted(set(stromal_genes) - set(present_stromal_genes))
        raise RuntimeError(f"{sample}: missing shared stromal genes: {missing}")
    measurements["program_shared_stromal"] = stromal_score

    spot_diameter_pixels = float(adata.uns["spatial"][sample]["scalefactors"]["spot_diameter_fullres"])
    microns_per_pixel = SPOT_DIAMETER_UM / spot_diameter_pixels
    x_um = boundary["spatial_x_fullres"].to_numpy(dtype=float) * microns_per_pixel
    y_um = boundary["spatial_y_fullres"].to_numpy(dtype=float) * microns_per_pixel
    positions = pd.DataFrame(
        {
            "barcode": boundary["barcode"].astype(str),
            "x_um": x_um - np.min(x_um),
            "y_um": y_um - np.min(y_um),
        }
    )

    coords = positions[["x_um", "y_um"]].to_numpy()
    tree = cKDTree(coords)
    nearest = tree.query(coords, k=2)[0][:, 1]
    neighbor_counts = np.array(
        [len(indices) - 1 for indices in tree.query_ball_point(coords, r=JUXTA_THRESHOLD_UM)],
        dtype=int,
    )
    if float(np.median(neighbor_counts)) != 6.0 or int(np.max(neighbor_counts)) > 6:
        raise RuntimeError(f"{sample}: 110-um geometry no longer represents the Visium six-neighbor lattice")

    matrix = measurements.to_numpy(dtype=float)
    if matrix.shape != (adata.n_obs, 20) or not np.isfinite(matrix).all():
        raise RuntimeError(f"{sample}: invalid MISTy measurement matrix")
    if np.max(np.abs(matrix.mean(axis=0))) > 1e-8 or np.max(np.abs(matrix.std(axis=0) - 1)) > 1e-8:
        raise RuntimeError(f"{sample}: standardized feature checks failed")

    measurements.insert(0, "barcode", measurements.index)
    measurement_path = OUTPUT_DIR / f"{sample}_misty_measurements_v1.tsv"
    position_path = OUTPUT_DIR / f"{sample}_misty_positions_v1.tsv"
    measurements.to_csv(measurement_path, sep="\t", index=False, float_format="%.10g")
    positions.to_csv(position_path, sep="\t", index=False, float_format="%.10g")
    return {
        "sample": sample,
        "n_spots": int(adata.n_obs),
        "n_features": int(matrix.shape[1]),
        "microns_per_pixel": microns_per_pixel,
        "nearest_neighbor_um_median": float(np.median(nearest)),
        "nearest_neighbor_um_min": float(np.min(nearest)),
        "nearest_neighbor_um_max": float(np.max(nearest)),
        "neighbors_within_110_um_min": int(np.min(neighbor_counts)),
        "neighbors_within_110_um_median": float(np.median(neighbor_counts)),
        "neighbors_within_110_um_max": int(np.max(neighbor_counts)),
        "h5ad_path": str(h5ad_path),
        "h5ad_sha256": sha256(h5ad_path),
        "boundary_path": str(boundary_path),
        "boundary_sha256": sha256(boundary_path),
        "measurements_path": str(measurement_path),
        "measurements_sha256": sha256(measurement_path),
        "positions_path": str(position_path),
        "positions_sha256": sha256(position_path),
    }


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    memory_gate("misty_input_preparation")
    stromal_genes = load_shared_stromal_genes()
    feature_rows: list[dict[str, str]] = []
    for source_label, feature_id in CELL_IDS.items():
        feature_rows.append(
            {
                "feature_id": feature_id,
                "feature_class": "cell2location_q05_abundance",
                "source_label": source_label,
                "transform": "zscore(log1p(q05 abundance)) within section",
                "genes": "",
                "selection_scope": "common prespecified 14-cell panel",
            }
        )
    for source_label in FIXED_PROGRAMS:
        feature_rows.append(
            {
                "feature_id": PROGRAM_IDS[source_label],
                "feature_class": "prespecified_gene_program",
                "source_label": source_label,
                "transform": "zscore(existing per-spot program score) within section",
                "genes": "recorded in BayesSpace boundary script",
                "selection_scope": "fixed before SPARK-X/Hotspot discovery",
            }
        )
    feature_rows.append(
        {
            "feature_id": "program_shared_stromal",
            "feature_class": "replicated_hotspot_program",
            "source_label": "GLMF1-M1 / GLMF2-M1 shared stromal-perivascular program",
            "transform": "per-gene zscore(log1p(CPM)); mean; final zscore within section",
            "genes": ";".join(stromal_genes),
            "selection_scope": "37-gene strongest reciprocal cross-section Hotspot overlap",
        }
    )
    feature_manifest_path = OUTPUT_DIR / "misty_feature_manifest_v1.tsv"
    pd.DataFrame(feature_rows).to_csv(feature_manifest_path, sep="\t", index=False)
    samples = [build_sample(sample, stromal_genes) for sample in SAMPLES]
    metadata = {
        "status": "complete",
        "method": "MISTy input preparation",
        "n_features": 20,
        "feature_composition": {
            "cell2location_q05_abundance": 14,
            "prespecified_gene_programs": 5,
            "replicated_hotspot_programs": 1,
        },
        "spatial_units": "micrometers",
        "juxtaview_threshold_um": JUXTA_THRESHOLD_UM,
        "paraview_family": "gaussian",
        "paraview_l_um": 200.0,
        "paraview_zoi_um": JUXTA_THRESHOLD_UM,
        "feature_manifest_path": str(feature_manifest_path),
        "feature_manifest_sha256": sha256(feature_manifest_path),
        "samples": samples,
    }
    metadata_path = OUTPUT_DIR / "misty_input_metadata_v1.json"
    metadata_path.write_text(json.dumps(metadata, indent=2, ensure_ascii=False) + "\n")
    print(f"MISTY_INPUT_PREPARATION_ALL_PASS samples={len(samples)} features=20")


if __name__ == "__main__":
    main()
