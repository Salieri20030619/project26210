#!/usr/bin/env python3
"""Validate completed GLMF1/2 mappings and build descriptive spatial summaries.

This analysis deliberately treats the two sections as two biological samples. It does
not calculate spot-level between-sample p-values. Cell2location factors are reference
components, not hard per-spot cell identities.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import math
import os
import platform
import sys
from datetime import datetime
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
from scipy.spatial import cKDTree
from scipy.stats import spearmanr


SAMPLES = ("GLMF1", "GLMF2")
ABUNDANCE_KEY = "q05_cell_abundance_w_sf"
ABUNDANCE_PREFIX = "q05cell_abundance_w_sf_"
EXPECTED_EPOCH = 10_000
EXPECTED_GENES = 8_830
EXPECTED_SPOTS = {"GLMF1": 3_199, "GLMF2": 3_880}
K_NEIGHBORS = 6


def now_iso() -> str:
    return datetime.now().astimezone().isoformat(timespec="seconds")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_text(path: Path, text: str) -> None:
    part = Path(str(path) + ".part")
    part.write_text(text, encoding="utf-8")
    os.replace(part, path)


def atomic_json(path: Path, value: object) -> None:
    atomic_text(path, json.dumps(value, indent=2, ensure_ascii=False) + "\n")


def atomic_tsv(frame: pd.DataFrame, path: Path, *, index: bool = False) -> None:
    part = Path(str(path) + ".part")
    compression = "gzip" if path.suffix == ".gz" else None
    frame.to_csv(part, sep="\t", index=index, compression=compression)
    os.replace(part, path)


def clean_state(name: str) -> str:
    return str(name).removeprefix(ABUNDANCE_PREFIX)


def read_manifest(sample_dir: Path) -> pd.DataFrame:
    path = sample_dir / "file_manifest_v4.tsv"
    frame = pd.read_csv(path, sep="\t")
    if list(frame.columns) != ["file", "size_bytes", "sha256"] or len(frame) != 3:
        raise RuntimeError(f"Unexpected manifest structure: {path}")
    for row in frame.itertuples(index=False):
        candidate = sample_dir / row.file
        if not candidate.is_file():
            raise FileNotFoundError(candidate)
        if candidate.stat().st_size != int(row.size_bytes):
            raise RuntimeError(f"Size mismatch: {candidate}")
        if sha256(candidate) != row.sha256:
            raise RuntimeError(f"SHA256 mismatch: {candidate}")
    return frame


def load_sample(input_root: Path, sample: str) -> tuple[pd.DataFrame, list[str], dict]:
    sample_dir = input_root / sample
    parts = list(sample_dir.rglob("*.part"))
    if parts:
        raise RuntimeError(f"Incomplete files found for {sample}: {parts}")
    status_path = sample_dir / "status.json"
    status = json.loads(status_path.read_text(encoding="utf-8"))
    if status.get("status") != "completed" or int(status.get("epoch", -1)) != EXPECTED_EPOCH:
        raise RuntimeError(f"Completion gate failed for {sample}: {status}")
    manifest = read_manifest(sample_dir)
    h5ad_name = f"{sample}_reference_v4_cell2location.h5ad"
    h5ad_path = sample_dir / h5ad_name
    abundance_path = sample_dir / "cell_abundance_q05_v4.tsv.gz"

    adata = ad.read_h5ad(h5ad_path, backed="r")
    try:
        if adata.n_obs != EXPECTED_SPOTS[sample] or adata.n_vars != EXPECTED_GENES:
            raise RuntimeError(f"Unexpected H5AD dimensions for {sample}: {adata.shape}")
        if ABUNDANCE_KEY not in adata.obsm or "spatial" not in adata.obsm:
            raise KeyError(f"Missing abundance or spatial coordinates for {sample}")
        abundance_h5 = adata.obsm[ABUNDANCE_KEY]
        abundance_array = np.asarray(abundance_h5, dtype=float)
        h5_columns = [clean_state(x) for x in getattr(abundance_h5, "columns", [])]
        if len(h5_columns) != abundance_array.shape[1]:
            h5_columns = [str(x) for x in adata.uns["mod"]["factor_names"]]
        abundance_h5 = pd.DataFrame(
            abundance_array,
            index=adata.obs_names.astype(str),
            columns=h5_columns,
        )
        if abundance_h5.shape[1] != 14:
            raise RuntimeError(f"Unexpected H5AD abundance shape for {sample}: {abundance_h5.shape}")
        coordinates = np.asarray(adata.obsm["spatial"], dtype=float)
        spot_ids = adata.obs_names.astype(str).tolist()
    finally:
        adata.file.close()

    abundance_tsv = pd.read_csv(abundance_path, sep="\t", index_col=0)
    abundance_tsv.index = abundance_tsv.index.astype(str)
    abundance_tsv.columns = [clean_state(x) for x in abundance_tsv.columns]
    states = abundance_tsv.columns.astype(str).tolist()
    if abundance_tsv.shape != (EXPECTED_SPOTS[sample], 14):
        raise RuntimeError(f"Unexpected abundance table shape for {sample}: {abundance_tsv.shape}")
    if abundance_tsv.index.tolist() != spot_ids:
        raise RuntimeError(f"Spot order mismatch between H5AD and abundance TSV for {sample}")
    values = abundance_tsv.to_numpy(dtype=float)
    if not np.isfinite(values).all() or (values < 0).any():
        raise RuntimeError(f"Invalid abundance values for {sample}")
    if abundance_h5.columns.tolist() != states:
        abundance_h5.columns = states
    if not np.allclose(abundance_h5.to_numpy(), values, rtol=2e-6, atol=2e-7):
        raise RuntimeError(f"H5AD/TSV abundance mismatch for {sample}")
    if coordinates.shape[0] != len(spot_ids) or coordinates.shape[1] < 2:
        raise RuntimeError(f"Unexpected spatial coordinate shape for {sample}: {coordinates.shape}")
    coordinates = coordinates[:, :2]
    if not np.isfinite(coordinates).all():
        raise RuntimeError(f"Non-finite spatial coordinates for {sample}")

    result = abundance_tsv.copy()
    result.insert(0, "y", coordinates[:, 1])
    result.insert(0, "x", coordinates[:, 0])
    result.insert(0, "spot_id", spot_ids)
    result.insert(0, "sample", sample)
    provenance = {
        "sample": sample,
        "status_path": str(status_path),
        "completed_at": status["updated_at"],
        "epoch": int(status["epoch"]),
        "n_spots": len(result),
        "n_genes": EXPECTED_GENES,
        "input_h5ad": str(h5ad_path),
        "input_h5ad_sha256": status["output_sha256"],
        "manifest_sha256_verified": True,
        "manifest_rows": int(len(manifest)),
    }
    return result, states, provenance


def safe_spearman(x: np.ndarray, y: np.ndarray) -> float:
    if np.ptp(x) == 0 or np.ptp(y) == 0:
        return math.nan
    return float(spearmanr(x, y).statistic)


def correlation_long(values: np.ndarray, states: list[str], sample: str, kind: str) -> pd.DataFrame:
    matrix = np.empty((len(states), len(states)), dtype=float)
    for i in range(len(states)):
        for j in range(len(states)):
            matrix[i, j] = safe_spearman(values[:, i], values[:, j])
    return pd.DataFrame(matrix, index=states, columns=states).rename_axis("state_1").reset_index().melt(
        id_vars="state_1", var_name="state_2", value_name="spearman_rho"
    ).assign(sample=sample, correlation_type=kind)[
        ["sample", "correlation_type", "state_1", "state_2", "spearman_rho"]
    ]


def neighbor_results(
    coordinates: np.ndarray, fractions: np.ndarray, states: list[str], sample: str
) -> tuple[pd.DataFrame, pd.DataFrame]:
    k = min(K_NEIGHBORS, len(coordinates) - 1)
    neighbor_index = cKDTree(coordinates).query(coordinates, k=k + 1)[1][:, 1:]
    lag = fractions[neighbor_index].mean(axis=1)
    correlation = np.empty((len(states), len(states)), dtype=float)
    for i in range(len(states)):
        for j in range(len(states)):
            correlation[i, j] = safe_spearman(fractions[:, i], lag[:, j])
    corr_long = pd.DataFrame(correlation, index=states, columns=states).rename_axis("state_1").reset_index().melt(
        id_vars="state_1", var_name="state_2", value_name="spearman_rho"
    ).assign(sample=sample, correlation_type="neighbor_lag_fraction", k_neighbors=k)[
        ["sample", "correlation_type", "k_neighbors", "state_1", "state_2", "spearman_rho"]
    ]

    moran_rows = []
    for index, state in enumerate(states):
        vector = fractions[:, index]
        centered = vector - vector.mean()
        denominator = float(np.square(centered).sum())
        moran = math.nan if denominator == 0 else float((centered[:, None] * centered[neighbor_index]).mean(axis=1).sum() / denominator)
        moran_rows.append(
            {"sample": sample, "state": state, "k_neighbors": k, "moran_i_descriptive": moran}
        )
    return corr_long, pd.DataFrame(moran_rows)


def figure_atomic(fig: plt.Figure, path: Path, *, dpi: int = 180) -> None:
    part = Path(str(path) + ".part")
    fig.savefig(part, dpi=dpi, bbox_inches="tight", format=path.suffix.lstrip("."))
    plt.close(fig)
    os.replace(part, path)


def plot_composition(summary: pd.DataFrame, states: list[str], path: Path) -> None:
    pivot = summary.pivot(index="sample", columns="state", values="mean_fraction").reindex(columns=states)
    colors = plt.get_cmap("tab20")(np.linspace(0, 1, len(states)))
    fig, ax = plt.subplots(figsize=(11, 5.5))
    bottom = np.zeros(len(pivot))
    for state, color in zip(states, colors):
        ax.bar(pivot.index, pivot[state], bottom=bottom, label=state, color=color, width=0.65)
        bottom += pivot[state].to_numpy()
    ax.set_ylabel("Mean within-spot q05 fraction")
    ax.set_ylim(0, 1)
    ax.set_title("GLMF1/GLMF2 descriptive cell-state composition")
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False, fontsize=8)
    figure_atomic(fig, path)


def plot_spatial(sample_frame: pd.DataFrame, states: list[str], figure_dir: Path) -> None:
    sample = str(sample_frame["sample"].iloc[0])
    colors = dict(zip(states, plt.get_cmap("tab20")(np.linspace(0, 1, len(states)))))
    fig, ax = plt.subplots(figsize=(7.2, 6.5))
    for state in states:
        subset = sample_frame[sample_frame["dominant_state"] == state]
        ax.scatter(subset["x"], subset["y"], s=4, color=colors[state], label=state, linewidths=0)
    ax.invert_yaxis()
    ax.set_aspect("equal")
    ax.set_title(f"{sample}: dominant reference component per spot")
    ax.set_xlabel("spatial x")
    ax.set_ylabel("spatial y")
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False, fontsize=7, markerscale=2)
    figure_atomic(fig, figure_dir / f"{sample}_dominant_state.png")

    fig, axes = plt.subplots(4, 4, figsize=(15, 14))
    for ax, state in zip(axes.flat, states):
        values = sample_frame[state].to_numpy(dtype=float)
        upper = float(np.quantile(values, 0.99))
        if upper <= 0:
            upper = float(values.max()) if values.max() > 0 else 1.0
        scatter = ax.scatter(sample_frame["x"], sample_frame["y"], c=values, s=3, cmap="magma", vmin=0, vmax=upper, linewidths=0)
        ax.invert_yaxis()
        ax.set_aspect("equal")
        ax.set_title(state, fontsize=9)
        ax.set_xticks([])
        ax.set_yticks([])
        fig.colorbar(scatter, ax=ax, fraction=0.046, pad=0.03)
    for ax in axes.flat[len(states):]:
        ax.axis("off")
    fig.suptitle(f"{sample}: q05 abundance maps (color capped at state-specific 99th percentile)", y=0.995)
    figure_atomic(fig, figure_dir / f"{sample}_state_abundance_panels.png")


def plot_heatmaps(correlations: pd.DataFrame, states: list[str], kind: str, path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(18, 7.8), constrained_layout=True)
    for ax, sample in zip(axes, SAMPLES):
        subset = correlations[(correlations["sample"] == sample) & (correlations["correlation_type"] == kind)]
        matrix = subset.pivot(index="state_1", columns="state_2", values="spearman_rho").reindex(index=states, columns=states)
        image = ax.imshow(matrix, cmap="coolwarm", vmin=-1, vmax=1)
        ax.set_title(sample)
        ax.set_xticks(range(len(states)), states, rotation=90, fontsize=7)
        ax.set_yticks(range(len(states)), states, fontsize=7)
    fig.colorbar(image, ax=axes, shrink=0.75, label="Spearman rho")
    fig.suptitle(kind.replace("_", " "))
    figure_atomic(fig, path)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    args = parser.parse_args()
    input_root = args.input_root.resolve(strict=True)
    output_root = args.output_root.resolve()
    if output_root.exists():
        raise FileExistsError(f"Refusing to overwrite output directory: {output_root}")
    output_root.mkdir(parents=True)
    table_dir = output_root / "tables"
    figure_dir = output_root / "figures"
    table_dir.mkdir()
    figure_dir.mkdir()
    atomic_json(output_root / "RUNNING.json", {"status": "RUNNING", "started_at": now_iso()})

    frames: dict[str, pd.DataFrame] = {}
    provenance: list[dict] = []
    states_expected: list[str] | None = None
    correlations = []
    moran = []
    summary_rows = []
    long_rows = []

    for sample in SAMPLES:
        frame, states, sample_provenance = load_sample(input_root, sample)
        if states_expected is None:
            states_expected = states
        elif states != states_expected:
            raise RuntimeError("GLMF1/GLMF2 factor names or order differ")
        values = frame[states].to_numpy(dtype=float)
        totals = values.sum(axis=1)
        if (totals <= 0).any():
            raise RuntimeError(f"Non-positive total q05 abundance in {sample}")
        fractions = values / totals[:, None]
        dominant_index = values.argmax(axis=1)
        frame["total_q05_abundance"] = totals
        frame["dominant_state"] = np.asarray(states)[dominant_index]
        frame["dominant_q05_abundance"] = values[np.arange(len(frame)), dominant_index]
        frame["dominant_fraction"] = fractions[np.arange(len(frame)), dominant_index]
        frames[sample] = frame
        provenance.append(sample_provenance)

        for state_index, state in enumerate(states):
            state_values = values[:, state_index]
            state_fractions = fractions[:, state_index]
            summary_rows.append(
                {
                    "sample": sample,
                    "state": state,
                    "n_spots": len(frame),
                    "total_q05_abundance": state_values.sum(),
                    "mean_q05_abundance": state_values.mean(),
                    "median_q05_abundance": np.median(state_values),
                    "mean_fraction": state_fractions.mean(),
                    "median_fraction": np.median(state_fractions),
                    "dominant_spots": int((dominant_index == state_index).sum()),
                    "dominant_spot_fraction": float((dominant_index == state_index).mean()),
                }
            )
            long_rows.append(
                pd.DataFrame(
                    {
                        "sample": sample,
                        "spot_id": frame["spot_id"],
                        "state": state,
                        "q05_abundance": state_values,
                        "within_spot_fraction": state_fractions,
                    }
                )
            )
        correlations.append(correlation_long(fractions, states, sample, "same_spot_fraction"))
        neighbor_corr, sample_moran = neighbor_results(
            frame[["x", "y"]].to_numpy(dtype=float), fractions, states, sample
        )
        correlations.append(neighbor_corr)
        moran.append(sample_moran)

    states = states_expected or []
    summary = pd.DataFrame(summary_rows)
    correlation_table = pd.concat(correlations, ignore_index=True)
    moran_table = pd.concat(moran, ignore_index=True)
    spot_summary = pd.concat(
        [frames[s][["sample", "spot_id", "x", "y", "total_q05_abundance", "dominant_state", "dominant_q05_abundance", "dominant_fraction"]] for s in SAMPLES],
        ignore_index=True,
    )
    abundance_long = pd.concat(long_rows, ignore_index=True)
    comparison = summary.pivot(index="state", columns="sample", values="mean_fraction").reindex(states).reset_index()
    comparison["GLMF2_minus_GLMF1_mean_fraction"] = comparison["GLMF2"] - comparison["GLMF1"]
    comparison["GLMF2_over_GLMF1_mean_fraction"] = comparison["GLMF2"] / comparison["GLMF1"].replace(0, np.nan)
    sample_rank_rho = safe_spearman(comparison["GLMF1"].to_numpy(), comparison["GLMF2"].to_numpy())

    off_diagonal = correlation_table[correlation_table["state_1"] != correlation_table["state_2"]].copy()
    same_spot_unique = off_diagonal[
        (off_diagonal["correlation_type"] == "same_spot_fraction")
        & (off_diagonal["state_1"] < off_diagonal["state_2"])
    ]
    neighbor_off = off_diagonal[off_diagonal["correlation_type"] == "neighbor_lag_fraction"]
    top_rows = []
    for sample in SAMPLES:
        for kind, source in (("same_spot_fraction", same_spot_unique), ("neighbor_lag_fraction", neighbor_off)):
            subset = source[source["sample"] == sample].dropna(subset=["spearman_rho"])
            for direction, selected in (
                ("positive", subset.nlargest(15, "spearman_rho")),
                ("negative", subset.nsmallest(15, "spearman_rho")),
            ):
                selected = selected.copy()
                selected["direction"] = direction
                selected["rank_within_direction"] = range(1, len(selected) + 1)
                top_rows.append(selected)
    top_colocalization = pd.concat(top_rows, ignore_index=True)

    atomic_tsv(spot_summary, table_dir / "spot_summary.tsv.gz")
    atomic_tsv(abundance_long, table_dir / "spot_state_abundance_long.tsv.gz")
    atomic_tsv(summary, table_dir / "sample_state_summary.tsv")
    atomic_tsv(comparison, table_dir / "glmf1_glmf2_descriptive_comparison.tsv")
    atomic_tsv(correlation_table, table_dir / "state_correlations_long.tsv.gz")
    atomic_tsv(moran_table, table_dir / "spatial_autocorrelation_moran.tsv")
    atomic_tsv(top_colocalization, table_dir / "top_colocalization_pairs.tsv")
    atomic_tsv(pd.DataFrame(provenance), table_dir / "input_provenance.tsv")

    plot_composition(summary, states, figure_dir / "sample_composition.png")
    for sample in SAMPLES:
        plot_spatial(frames[sample], states, figure_dir)
    plot_heatmaps(correlation_table, states, "same_spot_fraction", figure_dir / "same_spot_correlations.png")
    plot_heatmaps(correlation_table, states, "neighbor_lag_fraction", figure_dir / "neighbor_lag_correlations.png")

    software = {
        "python": platform.python_version(),
        "anndata": importlib.metadata.version("anndata"),
        "pandas": pd.__version__,
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "matplotlib": matplotlib.__version__,
        "script_sha256": sha256(Path(__file__).resolve()),
    }
    atomic_json(output_root / "software_versions.json", software)

    generated = sorted(
        path for path in output_root.rglob("*")
        if path.is_file() and path.name not in {"RUNNING.json", "OUTPUT_MANIFEST.tsv", "PHASE12_SPATIAL_ANALYSIS_STATUS.json"}
        and not path.name.endswith(".part")
    )
    manifest = pd.DataFrame(
        [{"file": str(path.relative_to(output_root)), "size_bytes": path.stat().st_size, "sha256": sha256(path)} for path in generated]
    )
    atomic_tsv(manifest, output_root / "OUTPUT_MANIFEST.tsv")
    completed = {
        "status": "PASS",
        "completed_at": now_iso(),
        "scope": "descriptive GLMF1/GLMF2 spatial cell-state analysis",
        "n_samples": 2,
        "n_spots": int(len(spot_summary)),
        "spots_by_sample": EXPECTED_SPOTS,
        "n_states": len(states),
        "states": states,
        "mean_fraction_rank_spearman_glmf1_vs_glmf2": sample_rank_rho,
        "inference_policy": "No between-sample hypothesis test; n=2 sections. Correlations and Moran I are descriptive.",
        "interpretation_policy": "Factors are reference components; dominant state is descriptive and is not a hard cell or malignancy call.",
        "input_provenance": provenance,
    }
    atomic_json(output_root / "PHASE12_SPATIAL_ANALYSIS_STATUS.json", completed)
    (output_root / "RUNNING.json").unlink()
    print(json.dumps(completed, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
