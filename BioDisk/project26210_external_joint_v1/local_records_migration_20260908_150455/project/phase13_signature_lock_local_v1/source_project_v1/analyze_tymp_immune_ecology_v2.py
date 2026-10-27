#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import os
import platform
import subprocess
import sys
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
from matplotlib.colors import TwoSlopeNorm
from scipy.sparse import issparse
from scipy.spatial import cKDTree
from scipy.stats import binomtest, rankdata, spearmanr


ROOT = Path("/mnt/BioDisk/ludingyu/project26210")
OUT = ROOT / "results/part4_spatial/tymp_immune_ecology_v2"
FIG_ORIGINAL = OUT / "figures/original"
FIG_COHORT = OUT / "figures/cohort"
TABLES = OUT / "tables"
QC = OUT / "qc"
MEMORY_GATE = ROOT / "check_memory_gate.sh"
ORIGINAL_INPUTS = {
    "GLMF1": ROOT / "results/part4_spatial/analysis/GLMF1_analyzed.h5ad",
    "GLMF2": ROOT / "results/part4_spatial/analysis/GLMF2_analyzed.h5ad",
}
COHORT_MANIFEST = (
    ROOT
    / "external/ccRCC_spatial/TJ_RCC_figshare_24599295"
    / "analysis_cell2location_full_v1/input_manifest.tsv"
)

SEED = 20260731
N_NEIGHBORS = 6
N_PERM = 9999
BLOCK_GRIDS = (4, 5)
PREFIX = "q05cell_abundance_w_sf_means_per_cluster_mu_fg_"
COHORT_PREFIX = "q05cell_abundance_w_sf_"

CELL_TARGETS = [
    "SPP1+ TAM",
    "APOE+C1QC+ TAM",
    "Other Macrophage",
    "Monocyte",
    "Dendritic cell",
    "Mast cell",
    "T cell",
    "NK cell",
    "B cell",
    "Endothelial",
]
MAP_CELL_STATES = [
    "TYMP-High Malignant",
    "SPP1+ TAM",
    "APOE+C1QC+ TAM",
    "Other Macrophage",
    "T cell",
    "NK cell",
    "Dendritic cell",
    "Endothelial",
]
CELL_SHORT = {
    "TYMP-High Malignant": "THM",
    "SPP1+ TAM": "SPP1-TAM",
    "APOE+C1QC+ TAM": "AC-TAM",
    "Other Macrophage": "Other-Mac",
    "Monocyte": "Mono",
    "Dendritic cell": "DC",
    "Mast cell": "Mast",
    "T cell": "T",
    "NK cell": "NK",
    "B cell": "B",
    "Endothelial": "Endo",
}
SIGNATURES = {
    "T-cell dysfunction": [
        "PDCD1",
        "CTLA4",
        "LAG3",
        "HAVCR2",
        "TIGIT",
        "TOX",
        "CXCL13",
        "ENTPD1",
    ],
    "Cytotoxicity": ["NKG7", "GNLY", "GZMB", "PRF1", "CTSW"],
    "Myeloid immunosuppression": [
        "SPP1",
        "CD163",
        "MRC1",
        "IL10",
        "TGFB1",
        "TREM2",
        "APOE",
        "C1QC",
        "VSIG4",
        "MARCO",
    ],
    "IFN response": ["IFNG", "STAT1", "IRF1", "CXCL9", "CXCL10", "ISG15", "IDO1"],
}
MTA_CONTEXT_GENES = [
    "TYMP",
    "MTAP",
    "MAT2A",
    "MAT2B",
    "AMD1",
    "SRM",
    "SMS",
    "ENOPH1",
    "APIP",
    "ADI1",
]
MTA_PRODUCTION_GENES = ["MAT2A", "MAT2B", "AMD1", "SRM", "SMS"]
EXPRESSION_MAPS = ["TYMP", "SPP1", "MTAP"]

SOURCE_URLS = {
    "TIDE": "https://doi.org/10.1038/s41591-018-0136-1",
    "MTA_macrophage": "https://pubmed.ncbi.nlm.nih.gov/35264604/",
    "SPP1_TAM_resistance": "https://doi.org/10.1038/s41586-024-08290-3",
}

FIGURE_ROWS: list[dict] = []
MEMORY_ROWS: list[dict] = []


def slug(text: str) -> str:
    replacements = {
        "+": "plus",
        "/": "_",
        " ": "_",
        "-": "_",
        "–": "_",
        "(": "",
        ")": "",
    }
    out = text
    for old, new in replacements.items():
        out = out.replace(old, new)
    return "".join(ch for ch in out if ch.isalnum() or ch == "_").strip("_").lower()


def sha256(path: Path, chunk: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            data = handle.read(chunk)
            if not data:
                break
            digest.update(data)
    return digest.hexdigest()


def memory_gate(stage: str, sample: str = "") -> None:
    result = subprocess.run(
        [str(MEMORY_GATE)],
        check=False,
        capture_output=True,
        text=True,
    )
    line = (result.stdout or result.stderr).strip()
    fields = {"stage": stage, "sample": sample, "exit_code": result.returncode, "line": line}
    MEMORY_ROWS.append(fields)
    print(f"[memory_gate] stage={stage} sample={sample} {line}", flush=True)
    if result.returncode != 0:
        raise RuntimeError(f"Memory gate blocked stage={stage}, sample={sample}: {line}")


def bh(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    out = np.full(values.shape, np.nan, dtype=float)
    finite = np.isfinite(values)
    if not finite.any():
        return out
    p = values[finite]
    order = np.argsort(p)
    ranked = p[order] * len(p) / np.arange(1, len(p) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adjusted = np.empty_like(ranked)
    adjusted[order] = np.minimum(ranked, 1.0)
    out[finite] = adjusted
    return out


def clean_abundance(adata: ad.AnnData) -> pd.DataFrame:
    abundance = adata.obsm["q05_cell_abundance_w_sf"]
    if not isinstance(abundance, pd.DataFrame):
        columns = adata.uns.get("mod", {}).get("factor_names")
        abundance = pd.DataFrame(abundance, index=adata.obs_names, columns=columns)
    abundance = abundance.copy()
    abundance.columns = [
        str(value).replace(PREFIX, "").replace(COHORT_PREFIX, "")
        for value in abundance.columns
    ]
    missing = sorted(set(["TYMP-High Malignant", *CELL_TARGETS]) - set(abundance.columns))
    if missing:
        raise ValueError(f"Missing abundance columns: {missing}")
    values = abundance.to_numpy(dtype=float)
    if not np.isfinite(values).all() or np.any(values < 0):
        raise ValueError("Abundance contains non-finite or negative values")
    return abundance


def neighbor_index(coords: np.ndarray) -> np.ndarray:
    _, indices = cKDTree(coords).query(coords, k=N_NEIGHBORS + 1)
    indices = np.asarray(indices[:, 1:], dtype=int)
    if indices.shape != (len(coords), N_NEIGHBORS):
        raise ValueError(f"Unexpected neighbor-index shape: {indices.shape}")
    if np.any(indices == np.arange(len(coords))[:, None]):
        raise ValueError("Self-neighbors remain after cKDTree query")
    return indices


def lag(values: np.ndarray, neighbors: np.ndarray) -> np.ndarray:
    return np.asarray(values, dtype=float)[neighbors].mean(axis=1)


def quantile_blocks(coords: np.ndarray, grid: int) -> np.ndarray:
    bins = []
    for axis in range(2):
        ranks = pd.Series(coords[:, axis]).rank(method="first")
        bins.append(pd.qcut(ranks, q=grid, labels=False).to_numpy())
    return bins[0] * grid + bins[1]


def block_summary(x: np.ndarray, y: np.ndarray, labels: np.ndarray) -> pd.DataFrame:
    frame = pd.DataFrame({"block": labels, "x": x, "y": y})
    return frame.groupby("block", as_index=False).agg(
        x=("x", "median"),
        y=("y", "median"),
        n_spots=("x", "size"),
    )


def fast_permutation_test(
    x: np.ndarray,
    y: np.ndarray,
    rng: np.random.Generator,
) -> tuple[float, float]:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    rho = float(spearmanr(x, y).statistic)
    if not np.isfinite(rho):
        return rho, np.nan
    xr = rankdata(x)
    yr = rankdata(y)
    xr = xr - xr.mean()
    yr = yr - yr.mean()
    denom = np.sqrt(np.sum(xr**2) * np.sum(yr**2))
    if denom == 0:
        return rho, np.nan
    permutations = np.vstack([rng.permutation(len(y)) for _ in range(N_PERM)])
    null = (yr[permutations] @ xr) / denom
    p_value = float((1 + np.count_nonzero(np.abs(null) >= abs(rho))) / (N_PERM + 1))
    return rho, p_value


def vectorized_spearman(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    result = spearmanr(np.asarray(x, float), np.asarray(y, float))
    return float(result.statistic), float(result.pvalue)


def normalized_required_expression(
    adata: ad.AnnData,
    required_genes: list[str],
) -> tuple[pd.DataFrame, np.ndarray]:
    genes = [gene for gene in required_genes if gene in adata.var_names]
    missing = sorted(set(required_genes) - set(genes))
    if missing:
        print(f"[gene_qc] missing={','.join(missing)}", flush=True)
    matrix = adata[:, genes].X
    if issparse(matrix):
        matrix = matrix.toarray()
    matrix = np.asarray(matrix, dtype=float)
    totals = np.asarray(adata.X.sum(axis=1)).reshape(-1).astype(float)
    if np.any(totals <= 0):
        raise ValueError("Zero-total spots found")
    log_cpm = np.log1p(matrix / totals[:, None] * 1e4)
    return pd.DataFrame(log_cpm, index=adata.obs_names, columns=genes), totals


def zscore(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    sd = values.std(ddof=0)
    if sd == 0 or not np.isfinite(sd):
        return np.zeros_like(values)
    return (values - values.mean()) / sd


def signature_scores(expression: pd.DataFrame) -> tuple[pd.DataFrame, list[dict]]:
    scores = {}
    qc_rows = []
    for name, requested in SIGNATURES.items():
        present = [gene for gene in requested if gene in expression.columns]
        if len(present) < 3:
            raise ValueError(f"Too few genes for signature {name}: {present}")
        matrix = np.column_stack([zscore(expression[gene].to_numpy()) for gene in present])
        scores[name] = matrix.mean(axis=1)
        qc_rows.append(
            {
                "signature": name,
                "requested_genes": ",".join(requested),
                "present_genes": ",".join(present),
                "missing_genes": ",".join(sorted(set(requested) - set(present))),
                "n_present": len(present),
            }
        )
    present_mta = [gene for gene in MTA_PRODUCTION_GENES if gene in expression.columns]
    matrix = np.column_stack([zscore(expression[gene].to_numpy()) for gene in present_mta])
    scores["Polyamine-MTA context"] = matrix.mean(axis=1)
    qc_rows.append(
        {
            "signature": "Polyamine-MTA context",
            "requested_genes": ",".join(MTA_PRODUCTION_GENES),
            "present_genes": ",".join(present_mta),
            "missing_genes": ",".join(sorted(set(MTA_PRODUCTION_GENES) - set(present_mta))),
            "n_present": len(present_mta),
        }
    )
    return pd.DataFrame(scores, index=expression.index), qc_rows


def save_figure(
    fig: plt.Figure,
    base_path: Path,
    group: str,
    sample: str,
    variable: str,
    plot_type: str,
    interpretation: str,
) -> None:
    base_path.parent.mkdir(parents=True, exist_ok=True)
    png = base_path.with_suffix(".png")
    pdf = base_path.with_suffix(".pdf")
    fig.savefig(png, dpi=300, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    FIGURE_ROWS.append(
        {
            "group": group,
            "sample": sample,
            "variable": variable,
            "plot_type": plot_type,
            "png": str(png),
            "pdf": str(pdf),
            "png_size_bytes": png.stat().st_size,
            "pdf_size_bytes": pdf.stat().st_size,
            "interpretation": interpretation,
        }
    )


def spatial_map(
    adata: ad.AnnData,
    sample: str,
    values: np.ndarray,
    title: str,
    colorbar_label: str,
    base_path: Path,
    group: str,
    variable: str,
    cmap: str = "magma",
) -> None:
    library = sample if sample in adata.uns["spatial"] else next(iter(adata.uns["spatial"]))
    spatial = adata.uns["spatial"][library]
    image = spatial["images"]["hires"]
    scale = float(spatial["scalefactors"]["tissue_hires_scalef"])
    coords = np.asarray(adata.obsm["spatial"], dtype=float) * scale
    values = np.asarray(values, dtype=float)
    finite = values[np.isfinite(values)]
    vmin = float(np.quantile(finite, 0.01))
    vmax = float(np.quantile(finite, 0.99))
    if vmax <= vmin:
        vmax = vmin + 1e-9
    fig, ax = plt.subplots(figsize=(5.4, 4.8))
    ax.imshow(image)
    points = ax.scatter(
        coords[:, 0],
        coords[:, 1],
        c=values,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        s=9,
        linewidths=0,
        alpha=0.92,
    )
    ax.set_title(title, fontsize=12)
    ax.set_axis_off()
    cbar = fig.colorbar(points, ax=ax, fraction=0.045, pad=0.025)
    cbar.set_label(colorbar_label, fontsize=9)
    cbar.ax.tick_params(labelsize=8)
    save_figure(
        fig,
        base_path,
        group,
        sample,
        variable,
        "spatial_map",
        "Spatial distribution; values are clipped only for color scaling at the 1st and 99th percentiles.",
    )


def association_heatmap(
    frame: pd.DataFrame,
    rows: list[str],
    columns: list[str],
    title: str,
    base_path: Path,
    group: str,
    sample: str,
    variable: str,
) -> None:
    matrix = frame.pivot(index="target", columns="grid", values="rho").reindex(
        index=rows, columns=columns
    )
    fdr = frame.pivot(index="target", columns="grid", values="fdr_bh").reindex(
        index=rows, columns=columns
    )
    fig_height = max(4.6, 0.45 * len(rows) + 1.5)
    fig, ax = plt.subplots(figsize=(5.3, fig_height))
    image = ax.imshow(
        matrix.to_numpy(dtype=float),
        cmap="coolwarm",
        norm=TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1),
        aspect="auto",
    )
    ax.set_xticks(range(len(columns)), columns)
    ax.set_yticks(range(len(rows)), [CELL_SHORT.get(x, x) for x in rows])
    ax.set_title(title, fontsize=11)
    for i in range(len(rows)):
        for j in range(len(columns)):
            value = matrix.iloc[i, j]
            q = fdr.iloc[i, j]
            star = "*" if np.isfinite(q) and q < 0.05 else ""
            ax.text(j, i, f"{value:.2f}{star}", ha="center", va="center", fontsize=8)
    cbar = fig.colorbar(image, ax=ax, fraction=0.05, pad=0.04)
    cbar.set_label("Spearman rho", fontsize=9)
    fig.text(0.5, 0.01, "*, BH FDR < 0.05", ha="center", fontsize=8)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    save_figure(
        fig,
        base_path,
        group,
        sample,
        variable,
        "association_heatmap",
        "Spatial-block effect sizes; stars indicate BH-adjusted block-permutation support.",
    )


def block_scatter(
    summary: pd.DataFrame,
    rho: float,
    p_value: float,
    fdr: float,
    title: str,
    x_label: str,
    y_label: str,
    base_path: Path,
    group: str,
    sample: str,
    variable: str,
) -> None:
    fig, ax = plt.subplots(figsize=(5.1, 4.4))
    ax.scatter(summary["x"], summary["y"], s=34, color="#3B6EA8", alpha=0.88)
    ax.set_title(title, fontsize=11)
    ax.set_xlabel(x_label, fontsize=9)
    ax.set_ylabel(y_label, fontsize=9)
    ax.text(
        0.03,
        0.97,
        f"rho={rho:.3f}\nPperm={p_value:.4g}\nFDR={fdr:.4g}",
        transform=ax.transAxes,
        va="top",
        fontsize=8,
        bbox={"facecolor": "white", "edgecolor": "0.75", "alpha": 0.85},
    )
    fig.tight_layout()
    save_figure(
        fig,
        base_path,
        group,
        sample,
        variable,
        "block_scatter",
        "Each point is a spatial-block median; inference uses 9,999 block-level permutations.",
    )


def mta_context_heatmap(
    sample: str,
    expression: pd.DataFrame,
    thm: np.ndarray,
    base_path: Path,
) -> None:
    quartile = pd.qcut(
        pd.Series(thm).rank(method="first"),
        q=4,
        labels=["Q1 low", "Q2", "Q3", "Q4 high"],
    )
    rows = []
    for gene in MTA_CONTEXT_GENES:
        if gene not in expression.columns:
            continue
        for q in quartile.cat.categories:
            rows.append(
                {
                    "gene": gene,
                    "quartile": str(q),
                    "mean_log_cpm": float(expression.loc[quartile.eq(q).to_numpy(), gene].mean()),
                }
            )
    frame = pd.DataFrame(rows)
    matrix = frame.pivot(index="gene", columns="quartile", values="mean_log_cpm").reindex(
        index=[g for g in MTA_CONTEXT_GENES if g in expression.columns],
        columns=["Q1 low", "Q2", "Q3", "Q4 high"],
    )
    scaled = matrix.apply(lambda row: pd.Series(zscore(row.to_numpy()), index=row.index), axis=1)
    fig, ax = plt.subplots(figsize=(5.5, 5.7))
    image = ax.imshow(
        scaled.to_numpy(dtype=float),
        cmap="coolwarm",
        norm=TwoSlopeNorm(vmin=-2, vcenter=0, vmax=2),
        aspect="auto",
    )
    ax.set_xticks(range(4), matrix.columns, rotation=25, ha="right")
    ax.set_yticks(range(len(matrix.index)), matrix.index)
    ax.set_title(f"{sample}: MTA-pathway context by THM abundance", fontsize=11)
    for i in range(len(matrix.index)):
        for j in range(4):
            ax.text(
                j,
                i,
                f"{matrix.iloc[i, j]:.2f}",
                ha="center",
                va="center",
                fontsize=7,
            )
    cbar = fig.colorbar(image, ax=ax, fraction=0.045, pad=0.04)
    cbar.set_label("Row z-score", fontsize=9)
    fig.text(0.5, 0.01, "Cell values show mean log1p(CPM); color is row-standardized.", ha="center", fontsize=8)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    save_figure(
        fig,
        base_path,
        "original",
        sample,
        "MTA-pathway context",
        "quartile_heatmap",
        "Transcriptomic context only; this is not a measurement or validated score of 5'-MTA accumulation.",
    )


def original_cell_tests(
    sample: str,
    abundance: pd.DataFrame,
    coords: np.ndarray,
    neighbors: np.ndarray,
) -> tuple[list[dict], dict[tuple[str, str, int], pd.DataFrame]]:
    rows = []
    scatters = {}
    relative = abundance.div(abundance.sum(axis=1).replace(0, np.nan), axis=0).fillna(0)
    sources = {
        "absolute_q05": abundance["TYMP-High Malignant"].to_numpy(dtype=float),
        "relative_abundance": relative["TYMP-High Malignant"].to_numpy(dtype=float),
    }
    for metric_i, (metric, source) in enumerate(sources.items()):
        target_frame = abundance if metric == "absolute_q05" else relative
        for target_i, target in enumerate(CELL_TARGETS):
            target_values = target_frame[target].to_numpy(dtype=float)
            target_lag = lag(target_values, neighbors)
            spot_rho, spot_p = vectorized_spearman(source, target_lag)
            rows.append(
                {
                    "sample": sample,
                    "metric": metric,
                    "source": "TYMP-High Malignant",
                    "target": target,
                    "scale": "spot_descriptive",
                    "grid": "",
                    "n_units": len(source),
                    "rho": spot_rho,
                    "p_value": spot_p,
                    "permutations": 0,
                }
            )
            for grid in BLOCK_GRIDS:
                summary = block_summary(source, target_lag, quantile_blocks(coords, grid))
                rho, p_value = fast_permutation_test(
                    summary["x"].to_numpy(),
                    summary["y"].to_numpy(),
                    np.random.default_rng(
                        SEED + (1 if sample == "GLMF2" else 0) * 100000
                        + metric_i * 10000
                        + target_i * 100
                        + grid
                    ),
                )
                rows.append(
                    {
                        "sample": sample,
                        "metric": metric,
                        "source": "TYMP-High Malignant",
                        "target": target,
                        "scale": "spatial_block_permutation",
                        "grid": f"{grid}x{grid}",
                        "n_units": len(summary),
                        "rho": rho,
                        "p_value": p_value,
                        "permutations": N_PERM,
                    }
                )
                if metric == "absolute_q05" and grid == 4:
                    scatters[(sample, target, grid)] = summary
    return rows, scatters


def original_score_tests(
    sample: str,
    abundance: pd.DataFrame,
    scores: pd.DataFrame,
    expression: pd.DataFrame,
    coords: np.ndarray,
    neighbors: np.ndarray,
) -> tuple[list[dict], dict[tuple[str, str, str, int], pd.DataFrame]]:
    rows = []
    scatters = {}
    sources = {
        "TYMP-High Malignant": abundance["TYMP-High Malignant"].to_numpy(dtype=float),
        "SPP1+ TAM": abundance["SPP1+ TAM"].to_numpy(dtype=float),
        "TYMP expression": expression["TYMP"].to_numpy(dtype=float),
    }
    targets = {
        **{name: scores[name].to_numpy(dtype=float) for name in SIGNATURES},
        "Polyamine-MTA context": scores["Polyamine-MTA context"].to_numpy(dtype=float),
        "MTAP expression": expression["MTAP"].to_numpy(dtype=float),
    }
    for source_i, (source_name, source) in enumerate(sources.items()):
        for target_i, (target_name, target) in enumerate(targets.items()):
            target_lag = lag(target, neighbors)
            spot_rho, spot_p = vectorized_spearman(source, target_lag)
            rows.append(
                {
                    "sample": sample,
                    "source": source_name,
                    "target": target_name,
                    "scale": "spot_descriptive",
                    "grid": "",
                    "n_units": len(source),
                    "rho": spot_rho,
                    "p_value": spot_p,
                    "permutations": 0,
                }
            )
            for grid in BLOCK_GRIDS:
                summary = block_summary(source, target_lag, quantile_blocks(coords, grid))
                rho, p_value = fast_permutation_test(
                    summary["x"].to_numpy(),
                    summary["y"].to_numpy(),
                    np.random.default_rng(
                        SEED + 500000
                        + (1 if sample == "GLMF2" else 0) * 100000
                        + source_i * 10000
                        + target_i * 100
                        + grid
                    ),
                )
                rows.append(
                    {
                        "sample": sample,
                        "source": source_name,
                        "target": target_name,
                        "scale": "spatial_block_permutation",
                        "grid": f"{grid}x{grid}",
                        "n_units": len(summary),
                        "rho": rho,
                        "p_value": p_value,
                        "permutations": N_PERM,
                    }
                )
                if grid == 4:
                    scatters[(sample, source_name, target_name, grid)] = summary
    return rows, scatters


def cohort_tests() -> pd.DataFrame:
    manifest = pd.read_csv(COHORT_MANIFEST, sep="\t")
    expected = set(
        [
            "sample",
            "path",
            "size_bytes",
            "sha256",
        ]
    )
    if not expected.issubset(manifest.columns):
        raise ValueError(f"Unexpected cohort manifest columns: {manifest.columns.tolist()}")
    rows = []
    for sample_i, record in manifest.sort_values("sample").reset_index(drop=True).iterrows():
        sample = str(record["sample"])
        memory_gate("cohort_sample", sample)
        path = Path(record["path"])
        if path.stat().st_size != int(record["size_bytes"]):
            raise ValueError(f"Input size changed for {sample}")
        adata = ad.read_h5ad(path)
        abundance = clean_abundance(adata)
        coords = np.asarray(adata.obsm["spatial"], dtype=float)
        neighbors = neighbor_index(coords)
        relative = abundance.div(abundance.sum(axis=1).replace(0, np.nan), axis=0).fillna(0)
        for metric_i, (metric, frame) in enumerate(
            [("absolute_q05", abundance), ("relative_abundance", relative)]
        ):
            source = frame["TYMP-High Malignant"].to_numpy(dtype=float)
            for target_i, target in enumerate(CELL_TARGETS):
                target_lag = lag(frame[target].to_numpy(dtype=float), neighbors)
                spot_rho, spot_p = vectorized_spearman(source, target_lag)
                rows.append(
                    {
                        "sample": sample,
                        "metric": metric,
                        "source": "TYMP-High Malignant",
                        "target": target,
                        "scale": "spot_descriptive",
                        "grid": "",
                        "n_units": len(source),
                        "rho": spot_rho,
                        "p_value": spot_p,
                        "permutations": 0,
                    }
                )
                for grid in BLOCK_GRIDS:
                    summary = block_summary(source, target_lag, quantile_blocks(coords, grid))
                    rho, p_value = fast_permutation_test(
                        summary["x"].to_numpy(),
                        summary["y"].to_numpy(),
                        np.random.default_rng(
                            SEED
                            + 1000000
                            + sample_i * 100000
                            + metric_i * 10000
                            + target_i * 100
                            + grid
                        ),
                    )
                    rows.append(
                        {
                            "sample": sample,
                            "metric": metric,
                            "source": "TYMP-High Malignant",
                            "target": target,
                            "scale": "spatial_block_permutation",
                            "grid": f"{grid}x{grid}",
                            "n_units": len(summary),
                            "rho": rho,
                            "p_value": p_value,
                            "permutations": N_PERM,
                        }
                    )
        del adata
    result = pd.DataFrame(rows)
    infer = result["scale"].eq("spatial_block_permutation")
    result["within_patient_fdr_bh"] = np.nan
    result.loc[infer, "within_patient_fdr_bh"] = bh(
        result.loc[infer, "p_value"].to_numpy()
    )
    return result


def cohort_summary(patient_tests: pd.DataFrame) -> pd.DataFrame:
    infer = patient_tests[patient_tests["scale"].eq("spatial_block_permutation")].copy()
    rows = []
    for (metric, target, grid), frame in infer.groupby(["metric", "target", "grid"], sort=False):
        rho = frame["rho"].to_numpy(dtype=float)
        positive = int(np.sum(rho > 0))
        negative = int(np.sum(rho < 0))
        zero = int(np.sum(rho == 0))
        nonzero = positive + negative
        p_value = (
            float(binomtest(positive, nonzero, p=0.5, alternative="two-sided").pvalue)
            if nonzero
            else 1.0
        )
        rows.append(
            {
                "metric": metric,
                "target": target,
                "grid": grid,
                "patients_tested": len(frame),
                "positive_patients": positive,
                "negative_patients": negative,
                "zero_patients": zero,
                "median_rho": float(np.median(rho)),
                "q25_rho": float(np.quantile(rho, 0.25)),
                "q75_rho": float(np.quantile(rho, 0.75)),
                "sign_test_p": p_value,
                "patients_within_test_fdr_lt_0_05": int(
                    np.sum(frame["within_patient_fdr_bh"] < 0.05)
                ),
            }
        )
    summary = pd.DataFrame(rows)
    summary["sign_test_fdr_bh"] = bh(summary["sign_test_p"].to_numpy())
    return summary


def cohort_target_plot(
    target: str,
    metric: str,
    patient_tests: pd.DataFrame,
    summary: pd.DataFrame,
) -> None:
    frame = patient_tests[
        patient_tests["scale"].eq("spatial_block_permutation")
        & patient_tests["target"].eq(target)
        & patient_tests["metric"].eq(metric)
    ].copy()
    pivot = frame.pivot(index="sample", columns="grid", values="rho").sort_index()
    fig, ax = plt.subplots(figsize=(6.0, 4.6))
    x = np.array([0, 1], dtype=float)
    for sample, row in pivot.iterrows():
        ax.plot(x, [row["4x4"], row["5x5"]], color="0.78", linewidth=0.8, zorder=1)
        ax.scatter(x, [row["4x4"], row["5x5"]], s=25, alpha=0.82, zorder=2)
    medians = pivot[["4x4", "5x5"]].median(axis=0)
    ax.scatter(x, medians, s=90, marker="D", color="#B22222", zorder=3, label="Median")
    ax.axhline(0, color="0.35", linewidth=0.8)
    ax.set_xticks(x, ["4x4", "5x5"])
    ax.set_xlim(-0.35, 1.35)
    ax.set_ylim(-1.05, 1.05)
    ax.set_ylabel("Patient spatial-block rho", fontsize=9)
    ax.set_title(
        f"THM → lag({CELL_SHORT[target]}), {metric.replace('_', ' ')}",
        fontsize=11,
    )
    annotations = []
    for grid in ["4x4", "5x5"]:
        row = summary[
            summary["metric"].eq(metric)
            & summary["target"].eq(target)
            & summary["grid"].eq(grid)
        ].iloc[0]
        annotations.append(
            f"{grid}: +{int(row.positive_patients)}/12, "
            f"median={row.median_rho:.2f}, FDR={row.sign_test_fdr_bh:.3g}"
        )
    fig.text(0.5, 0.01, "\n".join(annotations), ha="center", fontsize=8)
    fig.tight_layout(rect=(0, 0.09, 1, 1))
    save_figure(
        fig,
        FIG_COHORT / f"cohort_thm_to_{slug(target)}_{metric}",
        "cohort",
        "TJ-RCC_12",
        target,
        f"patient_paired_{metric}",
        "Each line is one patient; cohort inference is an exact two-sided sign test with BH correction.",
    )


def cohort_summary_heatmap(summary: pd.DataFrame, metric: str) -> None:
    frame = summary[summary["metric"].eq(metric)]
    matrix = frame.pivot(index="target", columns="grid", values="median_rho").reindex(
        index=CELL_TARGETS, columns=["4x4", "5x5"]
    )
    fdr = frame.pivot(index="target", columns="grid", values="sign_test_fdr_bh").reindex(
        index=CELL_TARGETS, columns=["4x4", "5x5"]
    )
    positive = frame.pivot(index="target", columns="grid", values="positive_patients").reindex(
        index=CELL_TARGETS, columns=["4x4", "5x5"]
    )
    fig, ax = plt.subplots(figsize=(5.8, 6.3))
    image = ax.imshow(
        matrix.to_numpy(dtype=float),
        cmap="coolwarm",
        norm=TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1),
        aspect="auto",
    )
    ax.set_xticks(range(2), ["4x4", "5x5"])
    ax.set_yticks(range(len(CELL_TARGETS)), [CELL_SHORT[x] for x in CELL_TARGETS])
    ax.set_title(f"TJ-RCC cohort: THM immune neighborhood ({metric.replace('_', ' ')})", fontsize=11)
    for i in range(len(CELL_TARGETS)):
        for j in range(2):
            star = "*" if fdr.iloc[i, j] < 0.05 else ""
            ax.text(
                j,
                i,
                f"{matrix.iloc[i, j]:.2f}{star}\n+{int(positive.iloc[i, j])}/12",
                ha="center",
                va="center",
                fontsize=7,
            )
    cbar = fig.colorbar(image, ax=ax, fraction=0.045, pad=0.04)
    cbar.set_label("Median patient rho", fontsize=9)
    fig.text(0.5, 0.01, "+n/12, positive patients; *, cohort sign-test FDR < 0.05.", ha="center", fontsize=8)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    save_figure(
        fig,
        FIG_COHORT / f"cohort_immune_neighborhood_summary_{metric}",
        "cohort",
        "TJ-RCC_12",
        "immune neighborhood",
        f"cohort_summary_{metric}",
        "Patient-level direction summary; patients are the biological replicates.",
    )


def write_contact_sheets() -> None:
    from PIL import Image, ImageDraw, ImageFont

    png_paths = [Path(row["png"]) for row in FIGURE_ROWS]
    chunk_size = 16
    thumb_w, thumb_h = 360, 300
    label_h = 42
    for sheet_index in range(0, len(png_paths), chunk_size):
        paths = png_paths[sheet_index : sheet_index + chunk_size]
        canvas = Image.new("RGB", (thumb_w * 4, (thumb_h + label_h) * 4), "white")
        draw = ImageDraw.Draw(canvas)
        for index, path in enumerate(paths):
            image = Image.open(path).convert("RGB")
            image.thumbnail((thumb_w - 10, thumb_h - 10))
            x = (index % 4) * thumb_w + (thumb_w - image.width) // 2
            y = (index // 4) * (thumb_h + label_h) + (thumb_h - image.height) // 2
            canvas.paste(image, (x, y))
            label = path.stem[:52]
            draw.text(
                ((index % 4) * thumb_w + 6, (index // 4) * (thumb_h + label_h) + thumb_h + 4),
                label,
                fill="black",
            )
        destination = QC / f"contact_sheet_{sheet_index // chunk_size + 1:02d}.png"
        canvas.save(destination)


def write_file_manifest() -> None:
    rows = []
    for path in sorted(OUT.rglob("*")):
        if (
            not path.is_file()
            or path.name == "file_manifest.tsv"
            or path.name in {"run.log", "verification_run.log"}
            or path.name.startswith("verification_")
        ):
            continue
        rows.append(
            {
                "relative_path": str(path.relative_to(OUT)),
                "size_bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
        )
    pd.DataFrame(rows).to_csv(OUT / "file_manifest.tsv", sep="\t", index=False)


def main() -> None:
    for directory in [FIG_ORIGINAL, FIG_COHORT, TABLES, QC]:
        directory.mkdir(parents=True, exist_ok=True)
    memory_gate("analysis_start")

    original_cell_rows = []
    original_score_rows = []
    original_spot_rows = []
    signature_qc_rows = []
    original_scatter = {}
    score_scatter = {}
    input_rows = []

    required_genes = sorted(
        set(
            EXPRESSION_MAPS
            + MTA_CONTEXT_GENES
            + [gene for genes in SIGNATURES.values() for gene in genes]
        )
    )

    for sample, path in ORIGINAL_INPUTS.items():
        memory_gate("original_sample", sample)
        input_rows.append(
            {
                "sample": sample,
                "path": str(path),
                "size_bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
        )
        adata = ad.read_h5ad(path)
        abundance = clean_abundance(adata)
        coords = np.asarray(adata.obsm["spatial"], dtype=float)
        neighbors = neighbor_index(coords)
        expression, total_counts = normalized_required_expression(adata, required_genes)
        scores, qc_rows = signature_scores(expression)
        for row in qc_rows:
            signature_qc_rows.append({"sample": sample, **row})

        for cell_state in MAP_CELL_STATES:
            spatial_map(
                adata,
                sample,
                abundance[cell_state].to_numpy(dtype=float),
                f"{sample}: {cell_state}",
                "q05 abundance",
                FIG_ORIGINAL / f"{sample}_{slug(cell_state)}_spatial",
                "original",
                cell_state,
                "viridis",
            )
        for gene in EXPRESSION_MAPS:
            spatial_map(
                adata,
                sample,
                expression[gene].to_numpy(dtype=float),
                f"{sample}: {gene} expression",
                "log1p(CPM)",
                FIG_ORIGINAL / f"{sample}_{slug(gene)}_expression_spatial",
                "original",
                gene,
                "magma",
            )
        for signature in SIGNATURES:
            spatial_map(
                adata,
                sample,
                scores[signature].to_numpy(dtype=float),
                f"{sample}: {signature} proxy",
                "mean gene z-score",
                FIG_ORIGINAL / f"{sample}_{slug(signature)}_proxy_spatial",
                "original",
                signature,
                "coolwarm",
            )
        mta_context_heatmap(
            sample,
            expression,
            abundance["TYMP-High Malignant"].to_numpy(dtype=float),
            FIG_ORIGINAL / f"{sample}_mta_pathway_context_by_thm_quartile",
        )

        for spot_index, barcode in enumerate(adata.obs_names):
            row = {
                "sample": sample,
                "barcode": barcode,
                "total_counts": total_counts[spot_index],
                "TYMP_High_Malignant_q05": abundance.iloc[spot_index][
                    "TYMP-High Malignant"
                ],
                "SPP1_TAM_q05": abundance.iloc[spot_index]["SPP1+ TAM"],
                "TYMP_log_cpm": expression.iloc[spot_index]["TYMP"],
                "SPP1_log_cpm": expression.iloc[spot_index]["SPP1"],
                "MTAP_log_cpm": expression.iloc[spot_index]["MTAP"],
            }
            for signature in scores.columns:
                row[f"score_{slug(signature)}"] = scores.iloc[spot_index][signature]
            original_spot_rows.append(row)

        cell_rows, scatters = original_cell_tests(sample, abundance, coords, neighbors)
        original_cell_rows.extend(cell_rows)
        original_scatter.update(scatters)
        score_rows, scatters = original_score_tests(
            sample, abundance, scores, expression, coords, neighbors
        )
        original_score_rows.extend(score_rows)
        score_scatter.update(scatters)
        del adata

    original_cell = pd.DataFrame(original_cell_rows)
    infer = original_cell["scale"].eq("spatial_block_permutation")
    original_cell["fdr_bh"] = np.nan
    original_cell.loc[infer, "fdr_bh"] = bh(original_cell.loc[infer, "p_value"].to_numpy())
    original_cell.to_csv(TABLES / "original_tymp_immune_spatial_tests.tsv", sep="\t", index=False)

    original_score = pd.DataFrame(original_score_rows)
    infer_score = original_score["scale"].eq("spatial_block_permutation")
    original_score["fdr_bh"] = np.nan
    original_score.loc[infer_score, "fdr_bh"] = bh(
        original_score.loc[infer_score, "p_value"].to_numpy()
    )
    original_score.to_csv(TABLES / "original_tymp_signature_spatial_tests.tsv", sep="\t", index=False)
    pd.DataFrame(original_spot_rows).to_csv(
        TABLES / "original_spot_scores.tsv.gz", sep="\t", index=False, compression="gzip"
    )
    pd.DataFrame(signature_qc_rows).to_csv(TABLES / "signature_gene_qc.tsv", sep="\t", index=False)
    pd.DataFrame(input_rows).to_csv(TABLES / "input_fingerprints.tsv", sep="\t", index=False)

    for sample in ORIGINAL_INPUTS:
        for metric in ["absolute_q05", "relative_abundance"]:
            frame = original_cell[
                original_cell["sample"].eq(sample)
                & original_cell["metric"].eq(metric)
                & original_cell["scale"].eq("spatial_block_permutation")
            ]
            association_heatmap(
                frame,
                CELL_TARGETS,
                ["4x4", "5x5"],
                f"{sample}: THM immune neighborhood ({metric.replace('_', ' ')})",
                FIG_ORIGINAL / f"{sample}_thm_immune_neighborhood_{metric}",
                "original",
                sample,
                metric,
            )

        score_frame = original_score[
            original_score["sample"].eq(sample)
            & original_score["source"].eq("TYMP-High Malignant")
            & original_score["scale"].eq("spatial_block_permutation")
        ].copy()
        score_rows_order = [*SIGNATURES.keys(), "Polyamine-MTA context", "MTAP expression"]
        score_frame = score_frame.rename(columns={"target": "target"})
        association_heatmap(
            score_frame,
            score_rows_order,
            ["4x4", "5x5"],
            f"{sample}: THM-associated transcriptomic contexts",
            FIG_ORIGINAL / f"{sample}_thm_transcriptomic_contexts",
            "original",
            sample,
            "transcriptomic contexts",
        )

    key_cell_targets = ["SPP1+ TAM", "T cell", "NK cell", "Dendritic cell", "Endothelial"]
    for sample in ORIGINAL_INPUTS:
        for target in key_cell_targets:
            summary = original_scatter[(sample, target, 4)]
            test = original_cell[
                original_cell["sample"].eq(sample)
                & original_cell["metric"].eq("absolute_q05")
                & original_cell["target"].eq(target)
                & original_cell["grid"].eq("4x4")
            ].iloc[0]
            block_scatter(
                summary,
                test["rho"],
                test["p_value"],
                test["fdr_bh"],
                f"{sample}: THM → lag({CELL_SHORT[target]})",
                "Block median THM q05",
                f"Block median lag({CELL_SHORT[target]}) q05",
                FIG_ORIGINAL / f"{sample}_thm_to_{slug(target)}_block_scatter",
                "original",
                sample,
                target,
            )

        for source, target in [
            ("SPP1+ TAM", "T-cell dysfunction"),
            ("SPP1+ TAM", "Cytotoxicity"),
            ("SPP1+ TAM", "Myeloid immunosuppression"),
            ("TYMP expression", "MTAP expression"),
            ("TYMP-High Malignant", "Polyamine-MTA context"),
        ]:
            summary = score_scatter[(sample, source, target, 4)]
            test = original_score[
                original_score["sample"].eq(sample)
                & original_score["source"].eq(source)
                & original_score["target"].eq(target)
                & original_score["grid"].eq("4x4")
            ].iloc[0]
            block_scatter(
                summary,
                test["rho"],
                test["p_value"],
                test["fdr_bh"],
                f"{sample}: {source} → lag({target})",
                f"Block median {source}",
                f"Block median lag({target})",
                FIG_ORIGINAL / f"{sample}_{slug(source)}_to_{slug(target)}_block_scatter",
                "original",
                sample,
                f"{source} to {target}",
            )

    memory_gate("cohort_analysis")
    patient_tests = cohort_tests()
    patient_tests.to_csv(TABLES / "tjrcc_patient_immune_spatial_tests.tsv", sep="\t", index=False)
    summary = cohort_summary(patient_tests)
    summary.to_csv(TABLES / "tjrcc_cohort_immune_direction_summary.tsv", sep="\t", index=False)

    for metric in ["absolute_q05", "relative_abundance"]:
        for target in CELL_TARGETS:
            cohort_target_plot(target, metric, patient_tests, summary)
        cohort_summary_heatmap(summary, metric)

    pd.DataFrame(FIGURE_ROWS).to_csv(OUT / "figure_manifest.tsv", sep="\t", index=False)

    metadata = {
        "status": "completed",
        "version": "v2",
        "created": pd.Timestamp.now(tz="Asia/Shanghai").isoformat(),
        "root": str(ROOT),
        "output": str(OUT),
        "original_inputs": {k: str(v) for k, v in ORIGINAL_INPUTS.items()},
        "cohort_manifest": str(COHORT_MANIFEST),
        "parameters": {
            "neighbors": N_NEIGHBORS,
            "block_grids": list(BLOCK_GRIDS),
            "permutations": N_PERM,
            "seed": SEED,
            "normalization": "log1p(count / spot_total_count * 1e4)",
            "signature_score": "mean within-sample z-score of available genes",
            "absolute_metric": "cell2location q05 abundance",
            "relative_metric": "q05 abundance divided by row sum across all cell states",
            "cohort_inference": "two-sided exact sign test across 12 patients; BH over all target-metric-grid tests",
        },
        "signatures": SIGNATURES,
        "mta_context_genes": MTA_CONTEXT_GENES,
        "sources": SOURCE_URLS,
        "limitations": [
            "The data do not directly measure 5'-MTA; MTA-pathway plots are transcriptomic context only.",
            "TYMP and MTAP are distinct enzymes; the analysis does not assert that TYMP directly causes MTA accumulation.",
            "Immune-state scores are spatial transcriptomic proxies and are not clinical immunotherapy-response labels.",
            "Visium spots mix cell types; signature scores are not cell-type-specific expression measurements.",
            "Spot-level associations are descriptive. Original-sample inference uses spatial-block permutations.",
            "GLMF1 and GLMF2 are analyzed separately and are not a treatment pair.",
            "TJ-RCC cohort inference treats 12 tumor cases as biological replicates and does not include R_cor/R_med.",
            "Spatial association does not prove direct communication, mediation, or causality.",
        ],
        "software": {
            "python": sys.version,
            "platform": platform.platform(),
            "anndata": ad.__version__,
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scipy": scipy.__version__,
            "matplotlib": matplotlib.__version__,
        },
        "n_figures": len(FIGURE_ROWS),
    }
    (OUT / "run_metadata.json").write_text(
        json.dumps(metadata, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    memory_gate("contact_sheet")
    write_contact_sheets()
    pd.DataFrame(MEMORY_ROWS).to_csv(OUT / "memory_gates.tsv", sep="\t", index=False)
    write_file_manifest()
    print(f"Completed: figures={len(FIGURE_ROWS)} output={OUT}", flush=True)


if __name__ == "__main__":
    main()
