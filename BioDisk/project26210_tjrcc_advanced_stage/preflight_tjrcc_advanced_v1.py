#!/usr/bin/env python3

from __future__ import annotations

import hashlib
import json
import math
import subprocess
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


PROJECT_ROOT = Path("/mnt/BioDisk/ludingyu/project26210")
TJ_ROOT = PROJECT_ROOT / "external/ccRCC_spatial/TJ_RCC_figshare_24599295"
OUTPUT_ROOT = Path("/mnt/BioDisk/wangmingyuan/project26210_tjrcc_advanced_stage/preflight_v3")
MEMORY_GATE = PROJECT_ROOT / "check_memory_gate.sh"

EXPECTED_CELLTYPES = {
    "APOE+C1QC+ TAM",
    "B cell",
    "Dendritic cell",
    "Endothelial",
    "Mast cell",
    "Epithelial-Normal",
    "NK cell",
    "SMC/Pericyte",
    "SPP1+ TAM",
    "T cell",
    "TYMP-High Malignant",
    "TYMP-Low Malignant",
    "Monocyte",
    "Other Macrophage",
}

PROGRAMS = {
    "Tumor-epithelial": ["CA9", "EPCAM", "KRT8", "KRT18", "KRT19"],
    "POSTN-ECM": ["POSTN", "COL1A1", "COL1A2", "COL3A1", "COL4A1", "FN1", "SPARC"],
    "APOE-TAM": ["APOE", "C1QA", "C1QB", "C1QC", "TREM2", "LPL"],
    "CD8-effector": ["CD8A", "CD8B", "NKG7", "CCL5", "PRF1", "GZMB"],
    "TGFb-response": ["TGFB1", "TGFBR1", "TGFBR2", "SMAD2", "SMAD3", "SERPINE1", "CTGF"],
}

COMMOT_AXES = {
    "TGFB1_TGFBR1_TGFBR2": ["TGFB1", "TGFBR1", "TGFBR2"],
    "POSTN_ITGAV_ITGB3": ["POSTN", "ITGAV", "ITGB3"],
    "POSTN_ITGAV_ITGB5": ["POSTN", "ITGAV", "ITGB5"],
    "CCL2_CCR2": ["CCL2", "CCR2"],
    "CSF1_CSF1R": ["CSF1", "CSF1R"],
    "SPP1_CD44": ["SPP1", "CD44"],
    "SPP1_ITGAV_ITGB3": ["SPP1", "ITGAV", "ITGB3"],
    "MIF_CD74_CXCR4": ["MIF", "CD74", "CXCR4"],
    "MIF_CD74_CD44": ["MIF", "CD74", "CD44"],
    "VEGFA_KDR": ["VEGFA", "KDR"],
    "APOE_TREM2_TYROBP": ["APOE", "TREM2", "TYROBP"],
    "APP_CD74": ["APP", "CD74"],
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def memory_gate(sample: str) -> str:
    result = subprocess.run(
        [str(MEMORY_GATE)],
        check=True,
        capture_output=True,
        text=True,
    )
    line = result.stdout.strip()
    print(f"MEMORY_GATE sample={sample} {line}", flush=True)
    return line


def sample_files() -> list[Path]:
    files = sorted((TJ_ROOT / "cell2location_pilot_v4").glob("*/*_cell2location_pilot.h5ad"))
    files += sorted((TJ_ROOT / "cell2location_remaining_v1").glob("*/*_cell2location.h5ad"))
    return sorted(files, key=lambda path: path.name.split("_cell2location", 1)[0])


def canonical_abundance(adata: ad.AnnData) -> pd.DataFrame:
    abundance = adata.obsm["q05_cell_abundance_w_sf"]
    if not isinstance(abundance, pd.DataFrame):
        factor_names = adata.uns["mod"]["factor_names"]
        abundance = pd.DataFrame(abundance, index=adata.obs_names, columns=factor_names)
    abundance = abundance.rename(
        columns=lambda value: str(value)
        .replace("q05cell_abundance_w_sf_means_per_cluster_mu_fg_", "")
        .replace("q05cell_abundance_w_sf_", "")
    )
    return abundance


def main() -> None:
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=False)
    files = sample_files()
    if len(files) != 12:
        raise RuntimeError(f"Expected 12 TJ-RCC H5AD files, found {len(files)}")

    rows: list[dict[str, object]] = []
    gene_rows: list[dict[str, object]] = []
    axis_rows: list[dict[str, object]] = []
    gates: list[str] = []
    for path in files:
        sample = path.name.split("_cell2location", 1)[0]
        gates.append(memory_gate(sample))
        adata = ad.read_h5ad(path)
        if not sparse.isspmatrix_csr(adata.X):
            raise RuntimeError(f"{sample}: X is not CSR")
        if not adata.obs_names.is_unique or not adata.var_names.is_unique:
            raise RuntimeError(f"{sample}: duplicate spot or gene identifiers")
        data = adata.X.data
        counts_pass = bool(
            np.all(np.isfinite(data))
            and np.all(data >= 0)
            and np.all(np.abs(data - np.rint(data)) <= 1e-8)
        )
        if not counts_pass:
            raise RuntimeError(f"{sample}: X failed finite nonnegative integer count validation")
        library = np.asarray(adata.X.sum(axis=1)).ravel().astype(float)
        if np.any(library <= 0):
            raise RuntimeError(f"{sample}: non-positive spot library size")
        spatial_xy = np.asarray(adata.obsm["spatial"], dtype=float)
        if spatial_xy.shape != (adata.n_obs, 2) or not np.all(np.isfinite(spatial_xy)):
            raise RuntimeError(f"{sample}: invalid spatial coordinates")
        abundance = canonical_abundance(adata)
        if set(abundance.columns) != EXPECTED_CELLTYPES:
            raise RuntimeError(f"{sample}: unexpected abundance columns {abundance.columns.tolist()}")
        abundance_values = abundance.to_numpy(dtype=float)
        if not np.all(np.isfinite(abundance_values)) or np.any(abundance_values < 0):
            raise RuntimeError(f"{sample}: invalid q05 abundance values")
        if sample not in adata.uns["spatial"]:
            raise RuntimeError(f"{sample}: sample key absent from uns['spatial']")
        spatial_uns = adata.uns["spatial"][sample]
        scalefactors = spatial_uns["scalefactors"]
        spot_diameter = float(scalefactors["spot_diameter_fullres"])
        hires_shape = tuple(np.asarray(spatial_uns["images"]["hires"]).shape)
        if not math.isfinite(spot_diameter) or spot_diameter <= 0 or len(hires_shape) != 3:
            raise RuntimeError(f"{sample}: invalid image or spot-diameter metadata")

        detection = np.asarray((adata.X > 0).sum(axis=0)).ravel().astype(int)
        detection_by_gene = dict(zip(adata.var_names.astype(str), detection, strict=True))
        min_spots = max(20, int(math.ceil(adata.n_obs * 0.01)))
        for program, genes in PROGRAMS.items():
            present = [gene for gene in genes if gene in detection_by_gene]
            detected = [gene for gene in present if detection_by_gene[gene] >= min_spots]
            gene_rows.append(
                {
                    "sample": sample,
                    "feature": program,
                    "feature_type": "program",
                    "genes_requested": ";".join(genes),
                    "genes_present": ";".join(present),
                    "genes_passing_detection_gate": ";".join(detected),
                    "n_requested": len(genes),
                    "n_present": len(present),
                    "n_passing_detection_gate": len(detected),
                    "minimum_detection_spots": min_spots,
                }
            )
        for axis, genes in COMMOT_AXES.items():
            present = [gene for gene in genes if gene in detection_by_gene]
            component_counts = {gene: int(detection_by_gene.get(gene, 0)) for gene in genes}
            axis_pass = all(component_counts[gene] >= min_spots for gene in genes)
            axis_rows.append(
                {
                    "sample": sample,
                    "axis_id": axis,
                    "genes": ";".join(genes),
                    "component_detection_counts": ";".join(
                        f"{gene}:{component_counts[gene]}" for gene in genes
                    ),
                    "minimum_detection_spots": min_spots,
                    "pass": axis_pass,
                }
            )

        rows.append(
            {
                "sample": sample,
                "h5ad": str(path),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
                "n_spots": adata.n_obs,
                "n_genes": adata.n_vars,
                "nnz": adata.X.nnz,
                "median_umi": float(np.median(library)),
                "median_detected_genes": float(np.median(np.asarray((adata.X > 0).sum(axis=1)).ravel())),
                "counts_integer_nonnegative": counts_pass,
                "spatial_coordinates_pass": True,
                "q05_abundance_pass": True,
                "hires_image_shape": "x".join(map(str, hires_shape)),
                "spot_diameter_fullres": spot_diameter,
                "minimum_detection_spots": min_spots,
            }
        )
        del adata, data, detection, abundance, abundance_values

    samples = pd.DataFrame(rows).sort_values("sample")
    programs = pd.DataFrame(gene_rows).sort_values(["sample", "feature"])
    axes = pd.DataFrame(axis_rows).sort_values(["sample", "axis_id"])
    samples.to_csv(OUTPUT_ROOT / "tjrcc_advanced_sample_preflight_v1.tsv", sep="\t", index=False)
    programs.to_csv(OUTPUT_ROOT / "tjrcc_advanced_program_gene_preflight_v1.tsv", sep="\t", index=False)
    axes.to_csv(OUTPUT_ROOT / "tjrcc_advanced_commot_axis_preflight_v1.tsv", sep="\t", index=False)

    metadata = {
        "status": "PASS",
        "n_samples": int(samples.shape[0]),
        "total_spots": int(samples["n_spots"].sum()),
        "total_nnz": int(samples["nnz"].sum()),
        "spot_range": [int(samples["n_spots"].min()), int(samples["n_spots"].max())],
        "all_counts_pass": bool(samples["counts_integer_nonnegative"].all()),
        "all_spatial_pass": bool(samples["spatial_coordinates_pass"].all()),
        "all_q05_pass": bool(samples["q05_abundance_pass"].all()),
        "commot_axes_pass_counts": {
            axis: int(group["pass"].sum()) for axis, group in axes.groupby("axis_id")
        },
        "memory_gates": gates,
        "interpretation": (
            "Compatibility preflight only; it does not establish biological validity or completion "
            "of BayesSpace, SPARK-X/Hotspot, MISTy, or COMMOT."
        ),
    }
    (OUTPUT_ROOT / "preflight_metadata_v1.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(
        "TJRCC_ADVANCED_PREFLIGHT_PASS "
        f"samples={samples.shape[0]} spots={samples['n_spots'].sum()} "
        f"nnz={samples['nnz'].sum()} commot_pass={int(axes['pass'].sum())}/{axes.shape[0]}"
    )


if __name__ == "__main__":
    main()
