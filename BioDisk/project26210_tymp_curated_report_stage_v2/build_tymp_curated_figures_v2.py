#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy.spatial import cKDTree


ROOT = Path("/mnt/BioDisk/ludingyu/project26210")
OUT = Path("/mnt/BioDisk/wangmingyuan/project26210_tymp_curated_report_stage_v2")
OUT.mkdir(parents=True, exist_ok=True)

IMMUNE = ROOT / "results/part4_spatial/tymp_immune_ecology_v2"
STATE = ROOT / "results/part4_spatial/tymp_spp1_state_v1"
SPOT_SCORES = IMMUNE / "tables/original_spot_scores.tsv.gz"
IMMUNE_TESTS = IMMUNE / "tables/original_tymp_immune_spatial_tests.tsv"
SIGNATURE_TESTS = IMMUNE / "tables/original_tymp_signature_spatial_tests.tsv"
PATIENT_SCORES = STATE / "tables/patient_spot_state_scores_v1.tsv.gz"
PATIENT_STATE = STATE / "tables/patient_state_associations_v1.tsv"
COHORT_STATE = STATE / "tables/cohort_state_direction_summary_v1.tsv"
PATIENT_DOWN = STATE / "tables/patient_interface_downstream_tests_v1.tsv"
COHORT_DOWN = STATE / "tables/cohort_interface_downstream_summary_v1.tsv"
H5ADS = {
    "GLMF1": ROOT / "results/part4_spatial/analysis/GLMF1_analyzed.h5ad",
    "GLMF2": ROOT / "results/part4_spatial/analysis/GLMF2_analyzed.h5ad",
}

FIG1_PNG = OUT / "Figure1_TYMP_MTA_SPP1_spatial_niche_8panel_v2.png"
FIG1_PDF = OUT / "Figure1_TYMP_MTA_SPP1_spatial_niche_8panel_v2.pdf"
FIG2_PNG = OUT / "Figure2_TYMP_SPP1_patient_spatial_atlas_6panel_v2.png"
FIG2_PDF = OUT / "Figure2_TYMP_SPP1_patient_spatial_atlas_6panel_v2.pdf"
PANEL_MANIFEST = OUT / "panel_source_manifest_v2.tsv"
SELECTED_EFFECTS = OUT / "selected_patient_effect_matrix_v2.tsv"
REPRESENTATIVE = OUT / "representative_patient_selection_v2.tsv"
METADATA = OUT / "figure_build_metadata_v2.json"

for required in [SPOT_SCORES, IMMUNE_TESTS, SIGNATURE_TESTS, PATIENT_SCORES, PATIENT_STATE, COHORT_STATE, PATIENT_DOWN, COHORT_DOWN, *H5ADS.values()]:
    if not required.is_file():
        raise FileNotFoundError(required)
for output in [FIG1_PNG, FIG1_PDF, FIG2_PNG, FIG2_PDF, PANEL_MANIFEST, SELECTED_EFFECTS, REPRESENTATIVE, METADATA]:
    if output.exists():
        raise FileExistsError(output)

plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 9,
        "axes.titlesize": 10.5,
        "axes.labelsize": 9,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "pdf.fonttype": 42,
    }
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def panel_letter(ax: plt.Axes, letter: str) -> None:
    ax.text(
        -0.035,
        1.035,
        letter,
        transform=ax.transAxes,
        fontsize=18,
        fontweight="bold",
        ha="left",
        va="top",
        color="#111111",
        bbox={"facecolor": "white", "edgecolor": "none", "pad": 1.2, "alpha": 0.92},
        zorder=100,
    )


def robust_limits(values: np.ndarray, lower: float = 0.02, upper: float = 0.98) -> tuple[float, float]:
    finite = values[np.isfinite(values)]
    low, high = np.quantile(finite, [lower, upper])
    if low == high:
        high = low + 1e-6
    return float(low), float(high)


def get_spatial(sample: str, scores: pd.DataFrame) -> tuple[pd.DataFrame, np.ndarray, np.ndarray, float]:
    adata = ad.read_h5ad(H5ADS[sample], backed="r")
    try:
        coords = np.asarray(adata.obsm["spatial"], dtype=float)
        spatial = adata.uns["spatial"][sample]
        image = np.asarray(spatial["images"]["hires"]).copy()
        scale = float(spatial["scalefactors"]["tissue_hires_scalef"])
        spot_diameter = float(spatial["scalefactors"]["spot_diameter_fullres"]) * scale
        frame = scores[scores["sample"].eq(sample)].set_index("barcode").reindex(adata.obs_names).copy()
        if frame.isna().all(axis=1).any() or len(frame) != adata.n_obs:
            raise RuntimeError(f"Barcode mismatch for {sample}")
        frame["x"] = coords[:, 0] * scale
        frame["y"] = coords[:, 1] * scale
        return frame, image, coords * scale, spot_diameter
    finally:
        adata.file.close()


def style_spatial_axis(ax: plt.Axes, image: np.ndarray, title: str) -> None:
    ax.imshow(image, alpha=0.58)
    ax.set_title(title, fontweight="bold", pad=5)
    ax.set_xlim(0, image.shape[1])
    ax.set_ylim(image.shape[0], 0)
    ax.set_aspect("equal")
    ax.axis("off")


def add_colorbar(ax: plt.Axes, artist, label: str) -> None:
    cax = inset_axes(ax, width="2.8%", height="42%", loc="lower right", borderpad=1.0)
    colorbar = plt.colorbar(artist, cax=cax)
    colorbar.ax.tick_params(labelsize=6.5, length=2)
    colorbar.set_label(label, fontsize=7)


def spatial_points(ax: plt.Axes, frame: pd.DataFrame, image: np.ndarray, column: str, title: str, cmap: str, colorbar_label: str) -> None:
    values = frame[column].to_numpy(dtype=float)
    low, high = robust_limits(values)
    style_spatial_axis(ax, image, title)
    artist = ax.scatter(
        frame["x"],
        frame["y"],
        c=np.clip(values, low, high),
        s=13,
        cmap=cmap,
        vmin=low,
        vmax=high,
        linewidths=0,
        alpha=0.94,
        rasterized=True,
    )
    add_colorbar(ax, artist, colorbar_label)


def spatial_mta_contour(ax: plt.Axes, frame: pd.DataFrame, image: np.ndarray, title: str) -> None:
    values = frame["score_polyamine_mta_context"].to_numpy(dtype=float)
    low, high = robust_limits(values)
    style_spatial_axis(ax, image, title)
    artist = ax.scatter(
        frame["x"], frame["y"], c=np.clip(values, low, high), s=11, cmap="magma", vmin=low, vmax=high,
        linewidths=0, alpha=0.78, rasterized=True,
    )
    triangulation = mtri.Triangulation(frame["x"].to_numpy(), frame["y"].to_numpy())
    triangles = triangulation.triangles
    x = frame["x"].to_numpy(); y = frame["y"].to_numpy()
    edges = np.stack(
        [
            np.hypot(x[triangles[:, 0]] - x[triangles[:, 1]], y[triangles[:, 0]] - y[triangles[:, 1]]),
            np.hypot(x[triangles[:, 1]] - x[triangles[:, 2]], y[triangles[:, 1]] - y[triangles[:, 2]]),
            np.hypot(x[triangles[:, 2]] - x[triangles[:, 0]], y[triangles[:, 2]] - y[triangles[:, 0]]),
        ],
        axis=1,
    )
    nearest = cKDTree(np.c_[x, y]).query(np.c_[x, y], k=2)[0][:, 1]
    triangulation.set_mask(edges.max(axis=1) > np.median(nearest) * 2.6)
    levels = np.quantile(values[np.isfinite(values)], [0.75, 0.90])
    ax.tricontour(triangulation, values, levels=levels, colors=["#00E5FF", "#FFF176"], linewidths=[0.8, 1.2])
    add_colorbar(ax, artist, "MTA-context z-score")
    ax.legend(
        [Line2D([0], [0], color="#00E5FF", lw=1.3), Line2D([0], [0], color="#FFF176", lw=1.5)],
        ["Q75 contour", "Q90 contour"], loc="upper right", fontsize=6.5, frameon=True, framealpha=0.82,
    )


def one_row(frame: pd.DataFrame, **filters) -> pd.Series:
    mask = pd.Series(True, index=frame.index)
    for key, value in filters.items():
        mask &= frame[key].eq(value)
    subset = frame[mask]
    if len(subset) != 1:
        raise RuntimeError(f"Expected one row for {filters}; got {len(subset)}")
    return subset.iloc[0]


def joint_niche(ax: plt.Axes, frame: pd.DataFrame, image: np.ndarray, sample: str, immune: pd.DataFrame, signature: pd.DataFrame) -> dict[str, int]:
    x = frame["x"].to_numpy(); y = frame["y"].to_numpy()
    thm = frame["TYMP_High_Malignant_q05"].to_numpy(dtype=float)
    mta = frame["score_polyamine_mta_context"].to_numpy(dtype=float)
    spp1 = frame["SPP1_TAM_q05"].to_numpy(dtype=float)
    thm_mta_core = (thm >= np.quantile(thm, 0.75)) & (mta >= np.quantile(mta, 0.75))
    spp1_high = spp1 >= np.quantile(spp1, 0.75)
    core_coords = np.c_[x[thm_mta_core], y[thm_mta_core]]
    median_nn = np.median(cKDTree(np.c_[x, y]).query(np.c_[x, y], k=2)[0][:, 1])
    adjacent = np.zeros(len(frame), dtype=bool)
    if len(core_coords):
        adjacent[spp1_high] = cKDTree(core_coords).query(np.c_[x[spp1_high], y[spp1_high]], k=1)[0] <= median_nn * 1.75
    style_spatial_axis(ax, image, f"{sample}: joint high-state niche overlay")
    ax.scatter(x, y, s=7, c="#BFC5CC", alpha=0.30, linewidths=0, rasterized=True)
    ax.scatter(x[thm_mta_core], y[thm_mta_core], s=15, c="#D81B60", alpha=0.86, linewidths=0, rasterized=True)
    ax.scatter(x[spp1_high], y[spp1_high], s=17, facecolors="none", edgecolors="#00A6A6", linewidths=0.55, alpha=0.9, rasterized=True)
    ax.scatter(x[adjacent], y[adjacent], s=24, facecolors="none", edgecolors="#FFD54F", linewidths=1.0, alpha=0.95, rasterized=True)
    tymp_mta = one_row(signature, sample=sample, source="TYMP-High Malignant", target="Polyamine-MTA context", scale="spatial_block_permutation", grid="5x5")
    spp1_mta = one_row(signature, sample=sample, source="SPP1+ TAM", target="Polyamine-MTA context", scale="spatial_block_permutation", grid="5x5")
    tymp_spp1 = one_row(immune, sample=sample, metric="absolute_q05", source="TYMP-High Malignant", target="SPP1+ TAM", scale="spatial_block_permutation", grid="5x5")
    ax.text(
        0.02,
        0.02,
        f"5×5 block FDR\nTYMP→MTA {tymp_mta.fdr_bh:.4f}\nSPP1-TAM↔MTA {spp1_mta.fdr_bh:.4f}\nTYMP→SPP1-TAM {tymp_spp1.fdr_bh:.4f}",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=6.7,
        bbox={"facecolor": "white", "edgecolor": "#666666", "alpha": 0.88, "boxstyle": "round,pad=0.25"},
    )
    ax.legend(
        handles=[
            Patch(facecolor="#D81B60", edgecolor="none", label="TYMP/MTA Q75 core"),
            Line2D([0], [0], marker="o", markerfacecolor="none", markeredgecolor="#00A6A6", linestyle="none", label="SPP1-TAM Q75"),
            Line2D([0], [0], marker="o", markerfacecolor="none", markeredgecolor="#FFD54F", linestyle="none", label="SPP1-TAM adjacent to core"),
        ],
        loc="upper right",
        fontsize=6.2,
        frameon=True,
        framealpha=0.86,
    )
    return {"core_spots": int(thm_mta_core.sum()), "spp1_high_spots": int(spp1_high.sum()), "adjacent_spp1_spots": int(adjacent.sum())}


spot_scores = pd.read_csv(SPOT_SCORES, sep="\t")
immune_tests = pd.read_csv(IMMUNE_TESTS, sep="\t")
signature_tests = pd.read_csv(SIGNATURE_TESTS, sep="\t")

spatial_inputs: dict[str, tuple[pd.DataFrame, np.ndarray, np.ndarray, float]] = {}
for sample in ["GLMF1", "GLMF2"]:
    spatial_inputs[sample] = get_spatial(sample, spot_scores)

fig, axes = plt.subplots(2, 4, figsize=(22, 11.5), constrained_layout=True)
fig.suptitle("Spatial architecture of the TYMP–MTA-context–SPP1+ macrophage niche", fontsize=20, fontweight="bold")
panel_rows: list[dict[str, str | int]] = []
joint_counts: dict[str, dict[str, int]] = {}
letters = iter("ABCDEFGH")
for row_index, sample in enumerate(["GLMF1", "GLMF2"]):
    frame, image, _, _ = spatial_inputs[sample]
    letter = next(letters); spatial_points(axes[row_index, 0], frame, image, "TYMP_High_Malignant_q05", f"{sample}: TYMP-high malignant abundance", "viridis", "q05 abundance"); panel_letter(axes[row_index, 0], letter)
    panel_rows.append({"figure": 1, "panel": letter, "type": "spatial_cell2location_map", "source": f"{H5ADS[sample]};{SPOT_SCORES}", "description": f"{sample} TYMP-high malignant q05 abundance"})
    letter = next(letters); spatial_mta_contour(axes[row_index, 1], frame, image, f"{sample}: Polyamine–MTA context"); panel_letter(axes[row_index, 1], letter)
    panel_rows.append({"figure": 1, "panel": letter, "type": "spatial_pathway_contour", "source": f"{H5ADS[sample]};{SPOT_SCORES}", "description": f"{sample} Polyamine-MTA context with Q75/Q90 contours"})
    letter = next(letters); spatial_points(axes[row_index, 2], frame, image, "SPP1_TAM_q05", f"{sample}: SPP1+ TAM abundance", "magma", "q05 abundance"); panel_letter(axes[row_index, 2], letter)
    panel_rows.append({"figure": 1, "panel": letter, "type": "spatial_cell2location_map", "source": f"{H5ADS[sample]};{SPOT_SCORES}", "description": f"{sample} SPP1+ TAM q05 abundance"})
    letter = next(letters); joint_counts[sample] = joint_niche(axes[row_index, 3], frame, image, sample, immune_tests, signature_tests); panel_letter(axes[row_index, 3], letter)
    panel_rows.append({"figure": 1, "panel": letter, "type": "spatial_joint_niche_overlay", "source": f"{H5ADS[sample]};{SPOT_SCORES};{IMMUNE_TESTS};{SIGNATURE_TESTS}", "description": f"{sample} descriptive Q75 overlay with formal 5x5 block FDR annotations"})

fig.savefig(FIG1_PNG, dpi=300, facecolor="white")
fig.savefig(FIG1_PDF, facecolor="white")
plt.close(fig)


patient_scores = pd.read_csv(PATIENT_SCORES, sep="\t")
patient_state = pd.read_csv(PATIENT_STATE, sep="\t")
cohort_state = pd.read_csv(COHORT_STATE, sep="\t")
patient_down = pd.read_csv(PATIENT_DOWN, sep="\t")
cohort_down = pd.read_csv(COHORT_DOWN, sep="\t")


def selected_state(route: str, outcome: str, grid: str) -> pd.DataFrame:
    result = patient_state[
        patient_state["route"].eq(route)
        & patient_state["adjustment"].eq("composition_qc")
        & patient_state["subset"].eq("macrophage_enriched")
        & patient_state["outcome"].eq(outcome)
        & patient_state["grid"].eq(grid)
    ][["sample", "rho"]].copy()
    if len(result) != 12:
        raise RuntimeError(f"Expected 12 rows for {route} {outcome} {grid}; got {len(result)}")
    return result


def selected_downstream(grid: str) -> pd.DataFrame:
    result = patient_down[
        patient_down["route"].eq("reference_dependent")
        & patient_down["adjustment"].eq("composition_qc")
        & patient_down["subset"].eq("macrophage_enriched")
        & patient_down["interface"].eq("TYMP_SPP1_coactivation")
        & patient_down["outcome"].eq("Invasion-EMT context")
        & patient_down["grid"].eq(grid)
    ][["sample", "rho"]].copy()
    if len(result) != 12:
        raise RuntimeError(f"Expected 12 rows for downstream {grid}; got {len(result)}")
    return result


effect_parts = []
for grid in ["4x4", "5x5"]:
    part = selected_state("reference_light", "SPP1_state_no_SPP1", grid).rename(columns={"rho": f"State {grid}"})
    effect_parts.append(part.set_index("sample"))
for grid in ["4x4", "5x5"]:
    part = selected_state("reference_dependent", "SPP1_expression", grid).rename(columns={"rho": f"Direct SPP1 {grid}"})
    effect_parts.append(part.set_index("sample"))
for grid in ["4x4", "5x5"]:
    part = selected_downstream(grid).rename(columns={"rho": f"Invasion-EMT {grid}"})
    effect_parts.append(part.set_index("sample"))
effects = pd.concat(effect_parts, axis=1)
if effects.shape != (12, 6) or effects.isna().any().any():
    raise RuntimeError(f"Unexpected effect matrix {effects.shape}")
effects["mean_selected_rho"] = effects.mean(axis=1)
effects = effects.sort_values("mean_selected_rho", ascending=False)
effects.to_csv(SELECTED_EFFECTS, sep="\t")

selection = effects[["State 5x5", "Direct SPP1 5x5", "Invasion-EMT 5x5"]].copy()
selection["minimum_selected_5x5_rho"] = selection.min(axis=1)
selection["rank"] = selection["minimum_selected_5x5_rho"].rank(method="first", ascending=False).astype(int)
selection = selection.sort_values("rank")
selection.to_csv(REPRESENTATIVE, sep="\t")
representative = str(selection.index[0])
if representative != "X49_T":
    raise RuntimeError(f"Representative patient changed unexpectedly: {representative}")


def patient_spatial(ax: plt.Axes, frame: pd.DataFrame, column: str, title: str, cmap: str, colorbar_label: str) -> None:
    values = frame[column].to_numpy(dtype=float)
    low, high = robust_limits(values)
    ax.set_facecolor("#F7F7F7")
    artist = ax.scatter(frame["spatial_x"], frame["spatial_y"], c=np.clip(values, low, high), cmap=cmap, vmin=low, vmax=high, s=10, linewidths=0, rasterized=True)
    ax.set_title(title, fontweight="bold")
    ax.set_aspect("equal")
    ax.invert_yaxis()
    ax.axis("off")
    add_colorbar(ax, artist, colorbar_label)


fig = plt.figure(figsize=(22, 12), constrained_layout=True)
grid = fig.add_gridspec(2, 4, height_ratios=[1.0, 1.25])
fig.suptitle("Patient-resolved spatial validation of the TYMP–SPP1 macrophage state", fontsize=20, fontweight="bold")
representative_scores = patient_scores[patient_scores["sample"].eq(representative)].copy()
if representative_scores.empty:
    raise RuntimeError("Representative patient scores missing")

patient_panels = [
    ("A", "neighbor_TYMP_reflight", "Neighboring TYMP state", "viridis", "state score"),
    ("B", "MTA_context", "Polyamine–MTA context", "magma", "context z-score"),
    ("C", "SPP1_state_no_SPP1", "Macrophage state without SPP1", "rocket", "state z-score"),
    ("D", "interface_coactivation", "TYMP–SPP1 state coactivation", "plasma", "coactivation"),
]
for column_index, (letter, column, label, cmap, colorbar_label) in enumerate(patient_panels):
    ax = fig.add_subplot(grid[0, column_index])
    patient_spatial(ax, representative_scores, column, f"{representative}: {label}", cmap, colorbar_label)
    panel_letter(ax, letter)
    panel_rows.append({"figure": 2, "panel": letter, "type": "representative_patient_spatial_map", "source": str(PATIENT_SCORES), "description": f"{representative} {column}; representative selected by maximum minimum selected 5x5 rho"})

heat_ax = fig.add_subplot(grid[1, 0:2])
heat_values = effects.drop(columns="mean_selected_rho")
sns.heatmap(heat_values, cmap="vlag", center=0, vmin=-0.75, vmax=0.75, annot=True, fmt=".2f", linewidths=0.45, linecolor="white", cbar_kws={"label": "Patient spatial-block rho"}, ax=heat_ax)
heat_ax.set_title("Twelve-patient effect matrix", fontweight="bold")
heat_ax.set_xlabel("Selected endpoint and block grid")
heat_ax.set_ylabel("Patient")
heat_ax.tick_params(axis="x", rotation=35)
heat_ax.tick_params(axis="y", rotation=0)
panel_letter(heat_ax, "E")
panel_rows.append({"figure": 2, "panel": "E", "type": "patient_effect_heatmap", "source": f"{PATIENT_STATE};{PATIENT_DOWN}", "description": "12 patients x six selected endpoint-grid effects"})

atlas_parent = fig.add_subplot(grid[1, 2:4])
atlas_parent.axis("off")
atlas_parent.set_title("Twelve-patient TYMP–SPP1 coactivation atlas", fontweight="bold", y=1.075, pad=0)
panel_letter(atlas_parent, "F")
atlas_grid = grid[1, 2:4].subgridspec(3, 4, wspace=0.04, hspace=0.12)
patient_order = list(effects.index)
for index, sample in enumerate(patient_order):
    ax = fig.add_subplot(atlas_grid[index // 4, index % 4])
    frame = patient_scores[patient_scores["sample"].eq(sample)]
    values = frame["interface_coactivation"].to_numpy(dtype=float)
    ax.scatter(frame["spatial_x"], frame["spatial_y"], c=values, cmap="plasma", vmin=0, vmax=1, s=2.4, linewidths=0, rasterized=True)
    ax.set_aspect("equal")
    ax.invert_yaxis()
    ax.axis("off")
    is_representative = sample == representative
    ax.text(
        0.02,
        0.90,
        f"{sample}{'  representative' if is_representative else ''}",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=7.1,
        color="white",
        fontweight="bold",
        bbox={
            "boxstyle": "round,pad=0.20",
            "facecolor": "#B7791F" if is_representative else "#263238",
            "edgecolor": "#FFD54F" if is_representative else "none",
            "linewidth": 1.0,
            "alpha": 0.92,
        },
    )
atlas_parent.text(0.98, 1.035, "Per-patient color scale: 0–1 percentile-product", transform=atlas_parent.transAxes, ha="right", va="bottom", fontsize=7.5, color="#555555")
panel_rows.append({"figure": 2, "panel": "F", "type": "twelve_patient_spatial_atlas", "source": str(PATIENT_SCORES), "description": "All 12 patients; TYMP-SPP1 coactivation percentile product"})

fig.savefig(FIG2_PNG, dpi=300, facecolor="white")
fig.savefig(FIG2_PDF, facecolor="white")
plt.close(fig)

source_hashes = {
    str(path): sha256(path)
    for path in [SPOT_SCORES, IMMUNE_TESTS, SIGNATURE_TESTS, PATIENT_SCORES, PATIENT_STATE, COHORT_STATE, PATIENT_DOWN, COHORT_DOWN, *H5ADS.values()]
}
panel_frame = pd.DataFrame(panel_rows).sort_values(["figure", "panel"])
panel_frame["source_sha256"] = panel_frame["source"].map(
    lambda value: ";".join(source_hashes[item] for item in str(value).split(";"))
)
panel_frame.to_csv(PANEL_MANIFEST, sep="\t", index=False)

metadata = {
    "status": "generated",
    "design_reference": {
        "document": str(ROOT / "分析待做（新）.docx"),
        "document_sha256": sha256(ROOT / "分析待做（新）.docx"),
        "part4_visual_types": ["spatial pathway contour", "niche/colocalization", "patient-resolved spatial atlas"],
    },
    "figure_1": {"path": str(FIG1_PNG), "pdf": str(FIG1_PDF), "panels": 8, "sha256": sha256(FIG1_PNG), "joint_overlay_counts": joint_counts},
    "figure_2": {"path": str(FIG2_PNG), "pdf": str(FIG2_PDF), "panels": 6, "sha256": sha256(FIG2_PNG), "representative_patient": representative, "selection_rule": "maximum of the minimum 5x5 rho across reference-light no-SPP1 macrophage state, reference-dependent direct SPP1, and invasion-EMT"},
    "source_hashes": source_hashes,
    "panel_manifest": str(PANEL_MANIFEST),
    "selected_effect_matrix": str(SELECTED_EFFECTS),
    "representative_selection": str(REPRESENTATIVE),
    "inference_boundary": "Spatial overlays and patient maps are descriptive visualizations; formal support remains the prespecified spatial-block and patient-direction statistics.",
}
METADATA.write_text(json.dumps(metadata, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
print(json.dumps(metadata, ensure_ascii=False, indent=2))
