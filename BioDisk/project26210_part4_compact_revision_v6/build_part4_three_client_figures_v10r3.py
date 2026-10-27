#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import importlib.util
import json
from datetime import datetime
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import TwoSlopeNorm
from matplotlib.lines import Line2D
from scipy.spatial import cKDTree


ROOT = Path("/mnt/BioDisk/ludingyu/project26210")
PART4 = ROOT / "results/part4_spatial"
OUT = Path("/mnt/BioDisk/wangmingyuan/project26210_part4_compact_revision_v6")
FIGDIR = OUT / "figures"
BASE_SCRIPT = Path(
    "/mnt/BioDisk/wangmingyuan/project26210_part4_three_figures_stage_v1/"
    "build_part4_three_client_figures_v1.py"
)

spec = importlib.util.spec_from_file_location("part4_fig_v4", BASE_SCRIPT)
base = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(base)

IMMUNE = PART4 / "tymp_immune_ecology_v2"
STATE = PART4 / "tymp_spp1_state_v1"
VASCULAR = PART4 / "tymp_vascular_niche_v1"
TRAJECTORY = PART4 / "tymp_myeloid_state_ordering_v1"
REFERENCE = PART4 / "ref/ref_signatures.h5ad"
GLMF2_H5AD = PART4 / "analysis/GLMF2_analyzed.h5ad"
NICHE_IMAGE = PART4 / "tymp_curated_validation_report_v2/Figure1_TYMP_MTA_SPP1_spatial_niche_8panel_v2.png"
METABOLIC_IMAGE = PART4 / "customer_style/figures/Fig4_customer_celltype_and_metabolic_GSVA.png"

SOURCES = {
    "immune_spot": IMMUNE / "tables/original_spot_scores.tsv.gz",
    "immune_tests": IMMUNE / "tables/original_tymp_immune_spatial_tests.tsv",
    "signature_tests": IMMUNE / "tables/original_tymp_signature_spatial_tests.tsv",
    "state_patient": STATE / "tables/patient_state_associations_v1.tsv",
    "state_cohort": STATE / "tables/cohort_state_direction_summary_v1.tsv",
    "mediation_patient": STATE / "tables/patient_exploratory_mediation_v1.tsv",
    "downstream_cohort": STATE / "tables/cohort_interface_downstream_summary_v1.tsv",
    "vascular_spot": VASCULAR / "tables/spot_mechanism_scores_v1.tsv.gz",
    "vascular_blocks": VASCULAR / "tables/sample_spatial_block_tests_v1.tsv",
    "vascular_cohort": VASCULAR / "tables/tjrcc_cohort_direction_summary_v1.tsv",
    "vascular_profiles": VASCULAR / "tables/boundary_distance_profiles_v1.tsv",
    "vascular_interface": VASCULAR / "tables/sample_interface_contrasts_v1.tsv",
    "trajectory_cells": TRAJECTORY / "tables/myeloid_trajectory_cell_scores_v1.tsv.gz",
    "trajectory_assoc": TRAJECTORY / "tables/sample_pseudotime_associations_v1.tsv",
    "trajectory_cohort": TRAJECTORY / "tables/cohort_pseudotime_direction_summary_v1.tsv",
    "trajectory_transition": TRAJECTORY / "tables/trajectory_neighbor_transition_matrix_v1.tsv",
    "trajectory_celltype": TRAJECTORY / "tables/celltype_trajectory_summary_v1.tsv",
    "reference": REFERENCE,
    "glmf2_h5ad": GLMF2_H5AD,
}

FIGURES = {
    1: FIGDIR / "Figure1_TYMP_MTA_SPP1_real_data_v10r3",
    2: FIGDIR / "Figure2_SPP1_vascular_niche_real_data_v10r3",
    3: FIGDIR / "Figure3_myeloid_continuum_real_data_v10r3",
}
MANIFEST = OUT / "panel_source_manifest_v10r3.tsv"
LEGENDS = OUT / "figure_legends_v10r3.md"
METADATA = OUT / "figure_generation_metadata_v10r3.json"

C = base.COLORS
CMAP = "RdBu_r"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def style() -> None:
    base.configure_style()
    plt.rcParams.update({
        "font.size": 9.5,
        "axes.titlesize": 11.0,
        "axes.labelsize": 9.0,
        "xtick.labelsize": 7.8,
        "ytick.labelsize": 7.8,
        "legend.fontsize": 7.5,
    })


def panel(ax: plt.Axes, letter: str) -> None:
    ax.text(
        -0.14, 1.055, letter, transform=ax.transAxes, fontsize=17,
        fontweight="bold", ha="left", va="top", color="#111827",
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.94, "pad": 0.7},
        zorder=200,
    )


def save(fig: plt.Figure, stem: Path) -> dict[str, object]:
    png = stem.with_suffix(".png")
    pdf = stem.with_suffix(".pdf")
    fig.savefig(png, dpi=360, bbox_inches="tight", facecolor="white")
    fig.savefig(pdf, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return {
        "png": str(png), "pdf": str(pdf),
        "png_size": png.stat().st_size, "pdf_size": pdf.stat().st_size,
        "png_sha256": sha256(png), "pdf_sha256": sha256(pdf),
    }


def heatmap(
    ax: plt.Axes,
    matrix: pd.DataFrame,
    title: str,
    *,
    vlim: float | None = None,
    annotate: bool = False,
    x_rotation: int = 0,
    cbar_label: str = "rho",
    number_format: str = ".2f",
) -> None:
    arr = matrix.to_numpy(dtype=float)
    if vlim is None:
        finite = np.abs(arr[np.isfinite(arr)])
        vlim = max(0.15, float(np.quantile(finite, 0.97))) if finite.size else 1.0
    image = ax.imshow(
        arr, cmap=CMAP, norm=TwoSlopeNorm(vmin=-vlim, vcenter=0, vmax=vlim),
        aspect="auto", interpolation="nearest",
    )
    ax.set_xticks(range(matrix.shape[1]), matrix.columns, rotation=x_rotation,
                  ha="right" if x_rotation else "center")
    ax.set_yticks(range(matrix.shape[0]), matrix.index)
    ax.set_title(title, fontweight="bold", pad=7)
    ax.tick_params(length=0, pad=3)
    ax.set_xticks(np.arange(-0.5, matrix.shape[1], 1), minor=True)
    ax.set_yticks(np.arange(-0.5, matrix.shape[0], 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=0.75)
    ax.tick_params(which="minor", bottom=False, left=False)
    if annotate:
        for i in range(arr.shape[0]):
            for j in range(arr.shape[1]):
                value = arr[i, j]
                if np.isfinite(value):
                    color = "white" if abs(value) > 0.56 * vlim else "#1F2937"
                    ax.text(j, i, format(value, number_format), ha="center", va="center",
                            fontsize=6.4, color=color)
    cbar = ax.figure.colorbar(image, ax=ax, fraction=0.032, pad=0.018)
    cbar.set_label(cbar_label, fontsize=7.7)
    cbar.ax.tick_params(labelsize=6.8, length=2)


def bivariate_quadrant_key(
    ax: plt.Axes,
    a_label: str,
    b_label: str,
    a_color: str,
    b_color: str,
    both_color: str,
) -> None:
    key = ax.inset_axes([0.00, -0.145, 1.00, 0.11])
    aa = np.array([0.05, .95, 0.05, .95])
    bb = np.array([0.05, 0.05, .95, .95])
    display = base.bivariate_colors(
        aa, bb, a_color=a_color, b_color=b_color, both=both_color
    ).reshape(1, 4, 3)
    key.imshow(display, aspect="auto")
    labels = ["both low", a_label, b_label, "both high"]
    for j in range(4):
        rgb = display[0, j]
        luminance = .2126 * rgb[0] + .7152 * rgb[1] + .0722 * rgb[2]
        key.text(j, 0, labels[j], ha="center", va="center", fontsize=7.4,
                 color="white" if luminance < .52 else "#111827", fontweight="bold")
    key.set_xticks([])
    key.set_yticks([])
    for spine in key.spines.values():
        spine.set_color("white")
        spine.set_linewidth(1.2)


def load_data() -> dict[str, pd.DataFrame]:
    data: dict[str, pd.DataFrame] = {}
    for name, path in SOURCES.items():
        if name in {"reference", "glmf2_h5ad"}:
            continue
        data[name] = pd.read_csv(path, sep="\t")
    return data


def add_glmf2_spatial(ax: plt.Axes, spot: pd.DataFrame) -> None:
    subset = spot[spot["sample"].eq("GLMF2")].set_index("barcode")
    backed = ad.read_h5ad(GLMF2_H5AD, backed="r")
    try:
        frame = subset.reindex(backed.obs_names).copy()
        coords = np.asarray(backed.obsm["spatial"], dtype=float)
        spatial = backed.uns["spatial"]["GLMF2"]
        image = np.asarray(spatial["images"]["hires"])
        scale = float(spatial["scalefactors"]["tissue_hires_scalef"])
    finally:
        backed.file.close()
    frame["x"] = coords[:, 0] * scale
    frame["y"] = coords[:, 1] * scale
    core = np.minimum(
        base.robust_rank(frame["TYMP_High_Malignant_q05"]),
        base.robust_rank(frame["score_polyamine_mta_context"]),
    )
    spp1 = base.robust_rank(frame["SPP1_TAM_q05"])
    colors = base.bivariate_colors(
        core, spp1, a_color=C["tumor"], b_color=C["spp1"], both="#6B3F87"
    )
    ax.imshow(image, alpha=0.60)
    ax.scatter(frame["x"], frame["y"], c=colors, s=11, lw=0, alpha=0.96, rasterized=True)
    ax.set(xlim=(0, image.shape[1]), ylim=(image.shape[0], 0))
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title("GLMF2: TYMP–MTA core and SPP1+ TAM", fontweight="bold", pad=5)
    bivariate_quadrant_key(
        ax, "core high", "SPP1 high", C["tumor"], C["spp1"], "#6B3F87"
    )


def add_glmf2_joint_niche(ax: plt.Axes, spot: pd.DataFrame) -> None:
    subset = spot[spot["sample"].eq("GLMF2")].set_index("barcode")
    backed = ad.read_h5ad(GLMF2_H5AD, backed="r")
    try:
        frame = subset.reindex(backed.obs_names).copy()
        coords = np.asarray(backed.obsm["spatial"], dtype=float)
        spatial = backed.uns["spatial"]["GLMF2"]
        image = np.asarray(spatial["images"]["hires"])
        scale = float(spatial["scalefactors"]["tissue_hires_scalef"])
    finally:
        backed.file.close()
    x, y = coords[:, 0] * scale, coords[:, 1] * scale
    tymp = frame["TYMP_High_Malignant_q05"].to_numpy(dtype=float)
    mta = frame["score_polyamine_mta_context"].to_numpy(dtype=float)
    spp1 = frame["SPP1_TAM_q05"].to_numpy(dtype=float)
    core = (tymp >= np.quantile(tymp, .75)) & (mta >= np.quantile(mta, .75))
    spp1_high = spp1 >= np.quantile(spp1, .75)
    median_nn = np.median(cKDTree(np.c_[x, y]).query(np.c_[x, y], k=2)[0][:, 1])
    adjacent = np.zeros(len(frame), dtype=bool)
    if core.any():
        distance = cKDTree(np.c_[x[core], y[core]]).query(np.c_[x[spp1_high], y[spp1_high]], k=1)[0]
        adjacent[spp1_high] = distance <= median_nn * 1.75
    ax.imshow(image, alpha=.58)
    ax.scatter(x, y, s=7, c="#BFC5CC", alpha=.28, lw=0, rasterized=True)
    ax.scatter(x[core], y[core], s=15, c="#D81B60", alpha=.87, lw=0, rasterized=True)
    ax.scatter(x[spp1_high], y[spp1_high], s=17, facecolors="none", edgecolors="#00A6A6",
               linewidths=.55, alpha=.92, rasterized=True)
    ax.scatter(x[adjacent], y[adjacent], s=25, facecolors="none", edgecolors="#F4B400",
               linewidths=1.0, alpha=.96, rasterized=True)
    ax.set(xlim=(0, image.shape[1]), ylim=(image.shape[0], 0))
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title("GLMF2: joint high-state niche overlay", fontweight="bold", pad=7)
    ax.legend(
        handles=[
            Line2D([0], [0], marker="o", markerfacecolor="#D81B60", markeredgecolor="none",
                   linestyle="none", label="TYMP/MTA Q75 core", markersize=6),
            Line2D([0], [0], marker="o", markerfacecolor="none", markeredgecolor="#00A6A6",
                   linestyle="none", label="SPP1-TAM Q75", markersize=6),
            Line2D([0], [0], marker="o", markerfacecolor="none", markeredgecolor="#F4B400",
                   linestyle="none", label="SPP1-TAM adjacent", markersize=6),
        ],
        loc="upper center", bbox_to_anchor=(.5, -.035), ncol=3, frameon=False,
        columnspacing=.8, handletextpad=.3, borderaxespad=0,
    )


def add_metabolic_pathway_maps(ax: plt.Axes) -> None:
    image = plt.imread(METABOLIC_IMAGE)
    width = image.shape[1]
    crop = image[:, int(0.445 * width):]
    ax.imshow(crop)
    ax.axis("off")
    ax.set_title("Four metabolic pathway activities across GLMF1 and GLMF2", fontweight="bold", pad=7)
    ax.text(-0.018, 0.69, "GLMF1", transform=ax.transAxes, rotation=90,
            ha="right", va="center", fontsize=8.2, fontweight="bold", color=C["ink"])
    ax.text(-0.018, 0.22, "GLMF2", transform=ax.transAxes, rotation=90,
            ha="right", va="center", fontsize=8.2, fontweight="bold", color=C["ink"])


def build_figure1(rows: list[dict[str, str]], d: dict[str, pd.DataFrame]) -> tuple[dict, dict]:
    fig = plt.figure(figsize=(20.5, 12.5))
    gs = fig.add_gridspec(2, 8, height_ratios=[1.0, 1.0])
    fig.subplots_adjust(left=0.055, right=0.975, bottom=0.075, top=0.90,
                        wspace=1.00, hspace=0.48)
    fig.suptitle(
        "TYMP–MTA–SPP1 macrophage hypothesis: spatial support and patient-level limits",
        fontsize=18.5, fontweight="bold", color=C["ink"], y=0.985,
    )

    ax = fig.add_subplot(gs[0, 0:2])
    add_glmf2_spatial(ax, d["immune_spot"])
    panel(ax, "A")
    rows.append({"figure":"1","panel":"A","type":"bivariate_spatial_map",
                 "source":str(SOURCES["immune_spot"]),
                 "description":"GLMF2 joint TYMP-MTA and SPP1 TAM spatial encoding"})

    ax = fig.add_subplot(gs[0, 2:4])
    add_glmf2_joint_niche(ax, d["immune_spot"])
    panel(ax, "B")
    rows.append({"figure":"1","panel":"B","type":"joint_niche_overlay",
                 "source":f"{SOURCES['immune_spot']};{GLMF2_H5AD}",
                 "description":"GLMF2 joint TYMP/MTA-core, SPP1-high and adjacent-spot overlay"})

    ax = fig.add_subplot(gs[0, 4:6])
    state = d["state_patient"]
    selected = state[(state["route"].eq("reference_light")) &
                     (state["adjustment"].eq("composition_qc")) &
                     (state["subset"].eq("macrophage_enriched")) &
                     (state["grid"].eq("5x5")) &
                     (state["outcome"].isin(["SPP1_state_no_SPP1", "SPP1_vs_C1QC_contrast"]))]
    patient_matrix = selected.pivot(index="outcome", columns="sample", values="rho")
    patient_matrix = patient_matrix.loc[["SPP1_state_no_SPP1", "SPP1_vs_C1QC_contrast"]]
    patient_matrix = patient_matrix[patient_matrix.iloc[0].sort_values().index]
    for sample in patient_matrix.columns:
        values = patient_matrix[sample].to_numpy(dtype=float)
        ax.plot([0, 1], values, color="#BCC5D1", lw=1.0, alpha=.78, zorder=1)
        ax.scatter(0, values[0], s=34, color=C["spp1_cell"], edgecolor="white", lw=.55, zorder=3)
        ax.scatter(1, values[1], s=34, color=C["c1qc"], edgecolor="white", lw=.55, zorder=3)
    for xpos, values, color in [
        (0, patient_matrix.iloc[0].to_numpy(dtype=float), C["spp1_cell"]),
        (1, patient_matrix.iloc[1].to_numpy(dtype=float), C["c1qc"]),
    ]:
        q1, med, q3 = np.quantile(values, [.25, .5, .75])
        ax.plot([xpos, xpos], [q1, q3], color="#1F2937", lw=5.0,
                solid_capstyle="round", alpha=.55, zorder=4)
        ax.scatter(xpos, med, s=70, marker="D", color="#1F2937",
                   edgecolor="white", lw=.7, zorder=5)
    ax.axhline(0, color="#555", ls="--", lw=.8)
    ax.set_xlim(-.35, 1.35)
    ax.set_xticks([0, 1], ["SPP1 state\n(SPP1 excluded)", "SPP1 minus\nC1QC state"])
    ax.set_ylabel("Patient-level rho")
    ax.set_title("Paired macrophage-state effects in 12 patients", fontweight="bold", pad=7)
    ax.text(.98, .03, "diamond = median; bar = IQR", transform=ax.transAxes,
            ha="right", fontsize=7.0, color=C["muted"])
    base.clean_axis(ax)
    panel(ax, "C")
    rows.append({"figure":"1","panel":"C","type":"paired_patient_effect_plot",
                 "source":str(SOURCES["state_patient"]),
                 "description":"Patient effects for SPP1-excluded and SPP1-versus-C1QC states"})

    ax = fig.add_subplot(gs[0, 6:8])
    state_cohort = d["state_cohort"]
    sc = state_cohort[(state_cohort["adjustment"].eq("composition_qc")) &
                      (state_cohort["subset"].eq("macrophage_enriched")) &
                      (state_cohort["outcome"].isin(["SPP1_state_no_SPP1", "SPP1_vs_C1QC_contrast"]))].copy()
    outcome_x = {"SPP1_state_no_SPP1":0, "SPP1_vs_C1QC_contrast":1}
    route_color = {"reference_dependent":"#7E57C2", "reference_light":"#2E75B6"}
    grid_offset = {"4x4":-.10, "5x5":.10}
    route_offset = {"reference_dependent":-.035, "reference_light":.035}
    for row in sc.itertuples(index=False):
        x = outcome_x[row.outcome] + grid_offset[row.grid] + route_offset[row.route]
        ax.errorbar(x, row.median_rho,
                    yerr=[[row.median_rho-row.q25_rho], [row.q75_rho-row.median_rho]],
                    fmt="o" if row.grid == "4x4" else "s", ms=5.4,
                    color=route_color[row.route], ecolor=route_color[row.route],
                    elinewidth=1.2, capsize=2.0, zorder=3)
    ax.axhline(0, color="#555", ls="--", lw=.8)
    ax.set_xticks([0, 1], ["SPP1 state\n(excluded)", "SPP1 minus\nC1QC"])
    ax.set_ylabel("Cohort median rho (IQR)")
    ax.set_title("State directions across route and scale", fontweight="bold", pad=7)
    ax.legend(handles=[
        Line2D([0],[0], marker="o", color="#555", label="4×4", linestyle="none", markersize=6),
        Line2D([0],[0], marker="s", color="#555", label="5×5", linestyle="none", markersize=6),
        Line2D([0],[0], marker="o", color=route_color["reference_dependent"], label="ref-dependent", linestyle="none", markersize=6),
        Line2D([0],[0], marker="o", color=route_color["reference_light"], label="ref-light", linestyle="none", markersize=6),
    ], ncol=2, frameon=False, loc="lower left", handletextpad=.3, columnspacing=.7)
    base.clean_axis(ax, grid=True)
    panel(ax, "D")
    rows.append({"figure":"1","panel":"D","type":"cohort_state_forest",
                 "source":str(SOURCES["state_cohort"]),
                 "description":"Cohort state medians and IQR across route and grid"})

    ax = fig.add_subplot(gs[1, 0:4])
    add_metabolic_pathway_maps(ax)
    panel(ax, "E")
    rows.append({"figure":"1","panel":"E","type":"metabolic_pathway_spatial_maps",
                 "source":str(METABOLIC_IMAGE),
                 "description":"Pyrimidine, Cys/Met, Arg/Pro and glycerophospholipid GSVA across GLMF1 and GLMF2"})

    ax = fig.add_subplot(gs[1, 4:6])
    down = d["downstream_cohort"].pivot(index="outcome", columns="grid", values="median_rho")
    down = down.reindex(columns=["4x4", "5x5"])
    short = {
        "Myeloid immunosuppression":"Myeloid\nimmunosuppression",
        "T-cell dysfunction":"T-cell\ndysfunction",
        "Cytotoxicity":"Cytotoxicity",
        "Hypoxia/angiogenesis":"Hypoxia /\nangiogenesis",
        "Invasion/EMT":"Invasion / EMT",
    }
    down.index = [short.get(x, x) for x in down.index]
    down = down.sort_values("5x5")
    for y, (_, row) in enumerate(down.iterrows()):
        ax.plot([row["4x4"], row["5x5"]], [y, y], color="#BFC7D2", lw=1.8, zorder=1)
        ax.scatter(row["4x4"], y, s=52, marker="o", color="#6BAED6",
                   edgecolor="white", lw=.65, zorder=3)
        ax.scatter(row["5x5"], y, s=58, marker="s", color="#D95F59",
                   edgecolor="white", lw=.65, zorder=3)
    ax.axvline(0, color="#555", ls="--", lw=.8)
    ax.set_yticks(range(len(down)), down.index)
    ax.set_xlabel("Cohort median rho")
    ax.set_title("Interface-associated downstream phenotypes", fontweight="bold", pad=7)
    ax.legend(handles=[
        Line2D([0],[0], marker="o", color="none", markerfacecolor="#6BAED6",
               markeredgecolor="white", label="4×4", markersize=7),
        Line2D([0],[0], marker="s", color="none", markerfacecolor="#D95F59",
               markeredgecolor="white", label="5×5", markersize=7),
    ], frameon=False, loc="lower right", ncol=2, handletextpad=.35, columnspacing=.8)
    base.clean_axis(ax, grid=True)
    panel(ax, "F")
    rows.append({"figure":"1","panel":"F","type":"two_scale_dumbbell",
                 "source":str(SOURCES["downstream_cohort"]),
                 "description":"Cohort median downstream associations at two grids"})

    ax = fig.add_subplot(gs[1, 6:8])
    med_all = d["mediation_patient"]
    med_all = med_all[(med_all["subset"].eq("macrophage_enriched")) &
                      (med_all["adjustment"].eq("composition_qc"))].copy()
    group_specs = [
        ("reference_dependent", "4x4", "Ref-dep\n4×4", "#9C7BC6"),
        ("reference_dependent", "5x5", "Ref-dep\n5×5", "#6E4AA5"),
        ("reference_light", "4x4", "Ref-light\n4×4", "#72A5CF"),
        ("reference_light", "5x5", "Ref-light\n5×5", "#2E75B6"),
    ]
    rng = np.random.default_rng(26211)
    for xpos, (route, grid, label, color) in enumerate(group_specs):
        values = med_all[(med_all["route"].eq(route)) & (med_all["grid"].eq(grid))]["indirect_a_times_b"].to_numpy()
        q1, medv, q3 = np.quantile(values, [.25, .5, .75])
        ax.plot([xpos, xpos], [q1, q3], color=color, lw=5, alpha=.5, solid_capstyle="round")
        ax.scatter(xpos + rng.uniform(-.11, .11, len(values)), values, s=27,
                   color=color, edgecolor="white", lw=.45, alpha=.86, zorder=3)
        ax.scatter(xpos, medv, s=62, marker="D", color="#1F2937", edgecolor="white", lw=.65, zorder=4)
    ax.axhline(0, color="#555", ls="--", lw=.8)
    ax.set_xticks(range(4), [x[2] for x in group_specs])
    ax.set_ylabel("Indirect effect a×b")
    ax.set_title("Patient mediation estimates by route and scale", fontweight="bold", pad=7)
    ax.text(.98, .03, "diamond = median; bar = IQR", transform=ax.transAxes,
            ha="right", fontsize=7.0, color=C["muted"])
    base.clean_axis(ax)
    panel(ax, "G")
    rows.append({"figure":"1","panel":"G","type":"mediation_distribution",
                 "source":str(SOURCES["mediation_patient"]),
                 "description":"Patient indirect-effect distributions across route and grid"})

    summary = {
        "patient_heatmap_patients": int(patient_matrix.shape[1]),
        "joint_niche_samples": 1,
        "metabolic_pathway_maps": 8,
        "downstream_outcomes": int(down.shape[0]),
    }
    return save(fig, FIGURES[1]), summary


RELATION_LABELS = {
    "SPP1_TAM_to_vascular_compartment":"SPP1 TAM → vascular",
    "SPP1_vascular_niche_to_T_cell_balance":"Niche → T-cell balance",
    "SPP1_vascular_niche_to_cytotoxicity":"Niche → cytotoxicity",
    "SPP1_vascular_niche_to_hypoxia_angiogenesis":"Niche → hypoxia / angio.",
    "SPP1_vascular_niche_to_invasion_EMT":"Niche → invasion / EMT",
    "TYMP_MTA_core_to_SPP1_monocyte_balance":"TYMP–MTA → SPP1 / Mono",
    "TYMP_MTA_core_to_SPP1_vascular_niche":"TYMP–MTA → niche",
}
ROUTE_TAG = {"route_independent":"", "reference_light":" [light]", "reference_dependent":" [dep]"}


def relation_display(relation: str, route: str) -> str:
    return RELATION_LABELS[relation] + ROUTE_TAG[route]


def build_figure2(rows: list[dict[str, str]], d: dict[str, pd.DataFrame]) -> tuple[dict, dict]:
    fig = plt.figure(figsize=(20.5, 12.5))
    gs = fig.add_gridspec(2, 12, height_ratios=[1.0, 1.02])
    fig.subplots_adjust(left=0.065, right=0.975, bottom=0.075, top=0.90,
                        wspace=1.05, hspace=0.48)
    fig.suptitle(
        "Patient-recurrent SPP1–vascular niche and its inverse alignment with T-cell balance",
        fontsize=18.5, fontweight="bold", color=C["ink"], y=0.985,
    )

    ax = fig.add_subplot(gs[0, 0:3])
    rep = d["vascular_spot"][d["vascular_spot"]["sample"].eq("Z43_T")].copy()
    niche = base.robust_rank(rep["spp1_vascular_niche"])
    tbal = base.robust_rank(rep["T_cell_balance"])
    colors = base.bivariate_colors(niche, tbal, a_color=C["spp1"], b_color=C["tcell"], both="#66518C")
    ax.scatter(rep["spatial_x"], rep["spatial_y"], c=colors, s=9.5, lw=0, rasterized=True)
    boundary = rep["computational_boundary_spot"].astype(bool).to_numpy()
    if boundary.any():
        ax.scatter(rep.loc[boundary, "spatial_x"], rep.loc[boundary, "spatial_y"], s=14,
                   facecolors="none", edgecolors="#222", lw=0.35, alpha=0.65, rasterized=True)
    ax.set_aspect("equal")
    ax.invert_yaxis()
    ax.axis("off")
    ax.set_title("Z43_T: niche and T-cell balance spatial map", fontweight="bold", pad=5)
    bivariate_quadrant_key(
        ax, "niche high", "T-cell high", C["spp1"], C["tcell"], "#66518C"
    )
    panel(ax, "A")
    rows.append({"figure":"2","panel":"A","type":"bivariate_spatial_map",
                 "source":str(SOURCES["vascular_spot"]),
                 "description":"Representative patient niche and T-cell balance map"})

    ax = fig.add_subplot(gs[0, 4:9])
    blocks = d["vascular_blocks"]
    b5 = blocks[(blocks["cohort"].eq("TJ-RCC")) & (blocks["grid"].eq("5x5"))].copy()
    b5["label"] = [relation_display(r, route) for r, route in zip(b5["relation"], b5["route"])]
    relation_order = (b5.groupby("label")["direction_aligned_rho"].median().sort_values().index)
    patient_order = (b5.groupby("sample")["direction_aligned_rho"].median().sort_values().index)
    matrix = b5.pivot(index="label", columns="sample", values="direction_aligned_rho")
    matrix = matrix.reindex(index=relation_order, columns=patient_order)
    heatmap(ax, matrix, "Patient-level direction-aligned effects for all prespecified relations (5×5)",
            vlim=0.75, annotate=False, x_rotation=42, cbar_label="")
    panel(ax, "B")
    rows.append({"figure":"2","panel":"B","type":"patient_relation_heatmap",
                 "source":str(SOURCES["vascular_blocks"]),
                 "description":"All prespecified 5x5 relations across 12 patients"})

    ax = fig.add_subplot(gs[0, 9:12])
    supported = [
        ("TYMP_MTA_core_to_SPP1_vascular_niche", "reference_dependent", "TYMP–MTA → niche", C["mta"]),
        ("SPP1_vascular_niche_to_T_cell_balance", "route_independent", "Niche → T-cell", C["spp1"]),
    ]
    all_xy = []
    for relation, route, label, color in supported:
        subset = blocks[(blocks["cohort"].eq("TJ-RCC")) &
                        (blocks["relation"].eq(relation)) &
                        (blocks["route"].eq(route))]
        wide = subset.pivot(index="sample", columns="grid", values="direction_aligned_rho")
        ax.scatter(wide["4x4"], wide["5x5"], s=42, color=color,
                   edgecolor="white", lw=.6, alpha=.88, label=label)
        all_xy.append(wide[["4x4", "5x5"]].to_numpy())
    xy = np.vstack(all_xy)
    lo = min(-.05, float(np.nanmin(xy)) - .05)
    hi = max(.75, float(np.nanmax(xy)) + .05)
    ax.plot([lo, hi], [lo, hi], color="#7B8794", ls="--", lw=1.0)
    ax.axhline(0, color="#B0B7C3", lw=.7)
    ax.axvline(0, color="#B0B7C3", lw=.7)
    ax.set(xlim=(lo, hi), ylim=(lo, hi))
    ax.set_xlabel("4×4 aligned rho")
    ax.set_ylabel("5×5 aligned rho")
    ax.set_title("Cross-scale patient concordance", fontweight="bold", pad=7)
    ax.legend(frameon=False, loc="lower right")
    base.clean_axis(ax, grid=True)
    panel(ax, "C")
    rows.append({"figure":"2","panel":"C","type":"cross_scale_scatter",
                 "source":str(SOURCES["vascular_blocks"]),
                 "description":"Patient concordance between 4x4 and 5x5 effects"})

    ax = fig.add_subplot(gs[1, 0:3])
    cohort = d["vascular_cohort"]
    c5 = cohort[cohort["grid"].eq("5x5")].copy()
    c5["label"] = [relation_display(r, route) for r, route in zip(c5["relation"], c5["route"])]
    c5 = c5.sort_values("median_direction_aligned_rho")
    route_color = {"route_independent":C["vascular"], "reference_light":"#D55E72", "reference_dependent":C["mta"]}
    for y, row in enumerate(c5.itertuples(index=False)):
        ax.plot([0, row.median_direction_aligned_rho], [y, y], color="#C6CDD7", lw=1.1)
        sig = row.sign_test_fdr_bh_global < 0.05
        ax.scatter(row.median_direction_aligned_rho, y,
                   s=38 + row.direction_supported_patients * 4.5,
                   color=route_color[row.route] if sig else "white",
                   edgecolor=route_color[row.route], lw=1.35, zorder=3)
    ax.axvline(0, color="#555", ls="--", lw=0.9)
    ax.set_yticks(range(len(c5)), c5["label"])
    ax.set_xlim(-0.45, 0.69)
    ax.set_xlabel("Median aligned rho")
    ax.set_title("Prespecified cohort audit", fontweight="bold", pad=7)
    ax.legend(handles=[
        Line2D([0],[0], marker="o", color="#666", markerfacecolor="#666",
               label="global FDR < 0.05", linestyle="none", markersize=6.5),
        Line2D([0],[0], marker="o", color="#666", markerfacecolor="white",
               label="global FDR ≥ 0.05", linestyle="none", markersize=6.5),
    ], frameon=False, loc="lower right", handletextpad=.35)
    base.clean_axis(ax, grid=True)
    panel(ax, "D")
    rows.append({"figure":"2","panel":"D","type":"cohort_audit_lollipop",
                 "source":str(SOURCES["vascular_cohort"]),
                 "description":"Cohort medians, support counts and global FDR"})

    ax = fig.add_subplot(gs[1, 3:6])
    profiles = d["vascular_profiles"]
    p = profiles[profiles["cohort"].eq("TJ-RCC")].copy()
    specs = [
        ("TYMP_MTA_core_reference_dependent", "TYMP–MTA core", C["tumor"]),
        ("spp1_vascular_niche", "SPP1–vascular niche", C["spp1"]),
        ("T_cell_balance", "T-cell balance", C["tcell"]),
    ]
    zones = p[["zone_order", "distance_zone"]].drop_duplicates().sort_values("zone_order")
    x = zones["zone_order"].to_numpy()
    labels = (zones["distance_zone"].str.replace("Peripheral ", "", regex=False)
              .str.replace("Distal ", "", regex=False))
    for score, label, color in specs:
        sub = p[p["score"].eq(score)]
        agg = sub.groupby("zone_order")["median_z"].agg(
            median="median", q25=lambda v: v.quantile(.25), q75=lambda v: v.quantile(.75)
        ).reindex(x)
        ax.plot(x, agg["median"], marker="o", lw=2.0, ms=4, color=color, label=label)
        ax.fill_between(x, agg["q25"], agg["q75"], color=color, alpha=.14, lw=0)
    ax.axhline(0, color="#555", ls="--", lw=.8)
    ax.set_xticks(x, labels, rotation=22, ha="right")
    ax.set_ylabel("Median within-sample z-score")
    ax.set_title("Profiles across computational interface", fontweight="bold", pad=7)
    ax.legend(frameon=False, loc="best")
    base.clean_axis(ax, grid=True)
    panel(ax, "E")
    rows.append({"figure":"2","panel":"E","type":"distance_profile",
                 "source":str(SOURCES["vascular_profiles"]),
                 "description":"Three cohort score profiles across interface zones"})

    ax = fig.add_subplot(gs[1, 6:9])
    supported_specs = [
        ("TYMP_MTA_core_to_SPP1_vascular_niche", "reference_dependent", "TYMP–MTA → niche"),
        ("SPP1_vascular_niche_to_T_cell_balance", "route_independent", "Niche → T-cell balance"),
    ]
    rng = np.random.default_rng(26210)
    positions = []
    labels = []
    colors_by_group = []
    for i, (relation, route, label) in enumerate(supported_specs):
        for j, grid in enumerate(["4x4", "5x5"]):
            vals = blocks[(blocks["cohort"].eq("TJ-RCC")) &
                          (blocks["relation"].eq(relation)) &
                          (blocks["route"].eq(route)) &
                          (blocks["grid"].eq(grid))]["direction_aligned_rho"].to_numpy()
            pos = i * 3 + j
            positions.append(pos)
            labels.append(f"{label}\n{grid}")
            color = C["mta"] if i == 0 else C["spp1"]
            colors_by_group.append(color)
            q1, med, q3 = np.quantile(vals, [.25, .5, .75])
            ax.plot([pos, pos], [q1, q3], color=color, lw=5.0, solid_capstyle="round", alpha=.45)
            ax.scatter(pos + rng.uniform(-.13, .13, len(vals)), vals, s=28, color=color,
                       edgecolor="white", lw=.45, alpha=.84, zorder=3)
            ax.scatter(pos, med, s=65, marker="D", color="#1F2937", edgecolor="white", lw=.7, zorder=4)
    ax.axhline(0, color="#555", ls="--", lw=.9)
    ax.set_xticks(positions, ["4×4", "5×5", "4×4", "5×5"])
    ax.text(.5, -.17, "TYMP–MTA → niche", transform=ax.get_xaxis_transform(),
            ha="center", va="top", fontsize=7.5, color=C["ink"])
    ax.text(3.5, -.17, "Niche → T-cell balance", transform=ax.get_xaxis_transform(),
            ha="center", va="top", fontsize=7.5, color=C["ink"])
    ax.set_ylabel("Direction-aligned patient rho")
    ax.set_title("Supported relations across two block scales", fontweight="bold", pad=7)
    ax.text(.98, .03, "diamond = median; bar = IQR", transform=ax.transAxes,
            ha="right", fontsize=7.0, color=C["muted"])
    base.clean_axis(ax)
    panel(ax, "F")
    rows.append({"figure":"2","panel":"F","type":"patient_distribution",
                 "source":str(SOURCES["vascular_blocks"]),
                 "description":"Patient distributions for two supported relations at two grids"})

    ax = fig.add_subplot(gs[1, 9:12])
    interface = d["vascular_interface"]
    interface = interface[(interface["cohort"].eq("TJ-RCC")) &
                          (interface["contrast"].eq("interface_minus_Distal_>800um")) &
                          (interface["score"].isin([
                              "TYMP_MTA_core_reference_dependent", "spp1_vascular_niche", "T_cell_balance"
                          ]))]
    interface_specs = [
        ("TYMP_MTA_core_reference_dependent", "TYMP–MTA core", C["tumor"]),
        ("spp1_vascular_niche", "SPP1–vascular niche", C["spp1"]),
        ("T_cell_balance", "T-cell balance", C["tcell"]),
    ]
    rng2 = np.random.default_rng(26212)
    for xpos, (score, label, color) in enumerate(interface_specs):
        values = interface[interface["score"].eq(score)]["difference_z"].to_numpy()
        q1, medv, q3 = np.quantile(values, [.25, .5, .75])
        ax.plot([xpos, xpos], [q1, q3], color=color, lw=5, alpha=.45, solid_capstyle="round")
        ax.scatter(xpos + rng2.uniform(-.12, .12, len(values)), values, s=27,
                   color=color, edgecolor="white", lw=.45, alpha=.86, zorder=3)
        ax.scatter(xpos, medv, s=62, marker="D", color="#1F2937", edgecolor="white", lw=.65, zorder=4)
    ax.axhline(0, color="#555", ls="--", lw=.8)
    ax.set_xticks(range(3), [x[1] for x in interface_specs], rotation=18, ha="right")
    ax.set_ylabel("Interface minus distal (>800 µm), z")
    ax.set_title("Patient interface contrasts", fontweight="bold", pad=7)
    base.clean_axis(ax)
    panel(ax, "G")
    rows.append({"figure":"2","panel":"G","type":"interface_contrast_distribution",
                 "source":str(SOURCES["vascular_interface"]),
                 "description":"Patient score contrasts between interface and distal zones"})

    summary = {"patients": int(matrix.shape[1]), "relations": int(matrix.shape[0]),
               "supported_relations": 2}
    return save(fig, FIGURES[2]), summary


def build_figure3(rows: list[dict[str, str]], d: dict[str, pd.DataFrame]) -> tuple[dict, dict]:
    cells = d["trajectory_cells"].copy()
    assoc = d["trajectory_assoc"].copy()
    cohort = d["trajectory_cohort"].copy()
    transition = d["trajectory_transition"].copy()
    celltype_summary = d["trajectory_celltype"].copy()
    backed = ad.read_h5ad(REFERENCE, backed="r")
    try:
        indexer = pd.Index(backed.obs_names).get_indexer(cells["cell_id"])
        if (indexer < 0).any():
            raise RuntimeError("Trajectory cell IDs missing from reference")
        coords = np.asarray(backed.obsm["X_umap"])[indexer]
    finally:
        backed.file.close()

    fig = plt.figure(figsize=(20.5, 12.5))
    gs = fig.add_gridspec(2, 8, height_ratios=[1.0, 1.0])
    fig.subplots_adjust(left=0.06, right=0.975, bottom=0.075, top=0.90,
                        wspace=1.08, hspace=0.48)
    fig.suptitle(
        "Myeloid cross-sectional continuum: C1QC-state gain without a proven SPP1 lineage",
        fontsize=18.5, fontweight="bold", color=C["ink"], y=0.985,
    )
    type_colors = {
        "Monocyte":C["monocyte"], "Other Macrophage":C["other_macro"],
        "APOE+C1QC+ TAM":C["apoe"], "SPP1+ TAM":C["spp1_cell"],
    }
    rng = np.random.default_rng(20260805)
    take = np.sort(rng.choice(len(cells), size=min(32000, len(cells)), replace=False))

    ax = fig.add_subplot(gs[0, 0:2])
    taken_types = cells["trajectory_celltype"].to_numpy()[take]
    for label, color in type_colors.items():
        mask = taken_types == label
        ax.scatter(coords[take][mask, 0], coords[take][mask, 1], s=3.0,
                   color=color, alpha=.58, lw=0, rasterized=True, label=label)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.set_title("Myeloid state map (32,000 cells shown)", fontweight="bold", pad=7)
    ax.legend(frameon=False, ncol=2, loc="lower left", borderaxespad=0.1,
              handletextpad=.35, columnspacing=.8, markerscale=2.5)
    base.clean_axis(ax)
    panel(ax, "A")
    rows.append({"figure":"3","panel":"A","type":"single_cell_umap",
                 "source":f"{SOURCES['reference']};{SOURCES['trajectory_cells']}",
                 "description":"Four myeloid states on the reference UMAP"})

    ax = fig.add_subplot(gs[0, 2:4])
    bins = np.linspace(cells["dpt_pseudotime"].min(), cells["dpt_pseudotime"].max(), 55)
    centers = (bins[:-1] + bins[1:]) / 2
    for offset, (label, color) in enumerate(type_colors.items()):
        vals = cells.loc[cells["trajectory_celltype"].eq(label), "dpt_pseudotime"].to_numpy()
        hist, _ = np.histogram(vals, bins=bins, density=True)
        hist = hist / (hist.max() or 1) * .78
        ax.fill_between(centers, offset, offset + hist, color=color, alpha=.68, lw=0)
        ax.plot(centers, offset + hist, color=color, lw=1.25)
    ax.set_yticks(np.arange(len(type_colors)) + .18, list(type_colors))
    ax.set_xlabel("Cross-sectional DPT")
    ax.set_title("State occupancy along DPT", fontweight="bold", pad=7)
    base.clean_axis(ax, grid=True)
    panel(ax, "B")
    rows.append({"figure":"3","panel":"B","type":"pseudotime_ridgeline",
                 "source":str(SOURCES["trajectory_cells"]),
                 "description":"DPT distributions for four annotated states"})

    score_specs = [
        ("Monocyte_program", "Monocyte program", C["monocyte"]),
        ("C1QC_state", "C1QC state", C["c1qc"]),
        ("SPP1_state_no_SPP1", "SPP1 state (excluded)", C["spp1_cell"]),
        ("SPP1_vs_C1QC", "SPP1 minus C1QC", C["mta"]),
    ]
    n_bins = 16
    cells["dpt_bin"] = pd.qcut(cells["dpt_pseudotime"], q=n_bins, labels=False, duplicates="drop")
    centers = cells.groupby("dpt_bin", observed=True)["dpt_pseudotime"].median().reindex(range(n_bins)).to_numpy()
    curve_medians: dict[str, pd.Series] = {}
    ax = fig.add_subplot(gs[0, 4:6])
    for score, label, color in score_specs:
        mean = cells.groupby("sample_id")[score].transform("mean")
        sd = cells.groupby("sample_id")[score].transform("std").replace(0, np.nan)
        cells["_z"] = (cells[score] - mean) / sd
        grouped = cells.groupby(["sample_id", "dpt_bin"], observed=True)["_z"].agg(["median", "size"]).reset_index()
        grouped = grouped[grouped["size"] >= 10]
        agg = grouped.groupby("dpt_bin")["median"].agg(
            median="median", q25=lambda x: x.quantile(.25), q75=lambda x: x.quantile(.75)
        ).reindex(range(n_bins))
        curve_medians[label] = agg["median"]
        ax.plot(centers, agg["median"], color=color, lw=2.0, label=label)
        ax.fill_between(centers, agg["q25"], agg["q75"], color=color, alpha=.12, lw=0)
    ax.axhline(0, color="#555", ls="--", lw=.8)
    ax.set_xlabel("Cross-sectional DPT")
    ax.set_ylabel("Median within-sample z-score")
    ax.set_title("Sample-balanced state-score profiles", fontweight="bold", pad=7)
    ax.legend(frameon=False, loc="best", ncol=1)
    base.clean_axis(ax, grid=True)
    panel(ax, "C")
    rows.append({"figure":"3","panel":"C","type":"sample_balanced_curves",
                 "source":str(SOURCES["trajectory_cells"]),
                 "description":"Within-sample standardized score profiles across DPT"})

    ax = fig.add_subplot(gs[0, 6:8])
    ct_order = ["Monocyte", "Other Macrophage", "APOE+C1QC+ TAM", "SPP1+ TAM"]
    ct = celltype_summary.set_index("trajectory_celltype").reindex(ct_order)
    for y, label in enumerate(ct_order):
        value = float(ct.loc[label, "median_pseudotime"])
        color = type_colors[label]
        size = 80 + 260 * np.sqrt(float(ct.loc[label, "n_cells"]) / float(ct["n_cells"].max()))
        ax.plot([0, value], [y, y], color="#C4CBD5", lw=1.4)
        ax.scatter(value, y, s=size, color=color, edgecolor="white", lw=.8, zorder=3)
    ax.set_yticks(range(4), ["Monocyte", "Other macro.", "APOE+C1QC+", "SPP1+"])
    ax.invert_yaxis()
    ax.set_xlim(0, .48)
    ax.set_xlabel("Median cross-sectional DPT")
    ax.set_title("State position and cell abundance", fontweight="bold", pad=7)
    base.clean_axis(ax, grid=True)
    panel(ax, "D")
    rows.append({"figure":"3","panel":"D","type":"state_position_bubble_lollipop",
                 "source":str(SOURCES["trajectory_celltype"]),
                 "description":"Median DPT and cell count across four annotated states"})

    ax = fig.add_subplot(gs[1, 0:2])
    curve_matrix = pd.DataFrame(curve_medians).T
    curve_matrix.index = ["Monocyte program", "C1QC state", "SPP1 state\n(SPP1 excluded)", "SPP1 minus C1QC"]
    stage_matrix = pd.DataFrame({
        "Early": curve_matrix.iloc[:, 0:4].median(axis=1),
        "Early-mid": curve_matrix.iloc[:, 4:8].median(axis=1),
        "Late-mid": curve_matrix.iloc[:, 8:12].median(axis=1),
        "Late": curve_matrix.iloc[:, 12:16].median(axis=1),
    })
    stage_values = stage_matrix.to_numpy(dtype=float)
    norm = TwoSlopeNorm(vmin=-1.05, vcenter=0, vmax=1.05)
    for i in range(stage_matrix.shape[0]):
        for j in range(stage_matrix.shape[1]):
            value = stage_values[i, j]
            ax.scatter(j, i, s=85 + 220 * abs(value), c=[value], cmap=CMAP,
                       norm=norm, edgecolor="white", lw=.8, zorder=3)
            ax.text(j, i, f"{value:.2f}", ha="center", va="center",
                    fontsize=6.7, color="white" if abs(value) > .58 else "#1F2937", zorder=4)
    ax.set_xticks(range(4), stage_matrix.columns)
    ax.set_yticks(range(4), stage_matrix.index)
    ax.invert_yaxis()
    ax.set_xlim(-.6, 3.6)
    ax.set_ylim(3.6, -.6)
    ax.set_xlabel("DPT stage (four equal-bin groups)")
    ax.set_title("State-score shifts across four DPT stages", fontweight="bold", pad=7)
    ax.grid(color=C["grid"], lw=.65, alpha=.75)
    ax.set_axisbelow(True)
    ax.tick_params(length=0)
    scalar = plt.cm.ScalarMappable(norm=norm, cmap=CMAP)
    cbar = fig.colorbar(scalar, ax=ax, fraction=.04, pad=.03)
    cbar.set_label("median z-score", fontsize=7.7)
    cbar.ax.tick_params(labelsize=6.8, length=2)
    panel(ax, "E")
    rows.append({"figure":"3","panel":"E","type":"dpt_stage_dot_matrix",
                 "source":str(SOURCES["trajectory_cells"]),
                 "description":"Four sample-balanced scores across 16 DPT bins"})

    ax = fig.add_subplot(gs[1, 2:4])
    order = list(type_colors)
    trans = transition.pivot(index="source", columns="target", values="source_normalized_fraction").reindex(index=order, columns=order)
    trans.index = ["Monocyte", "Other macro.", "APOE+C1QC+", "SPP1+"]
    trans.columns = ["Monocyte", "Other macro.", "APOE+C1QC+", "SPP1+"]
    image = ax.imshow(trans.to_numpy(), cmap="Blues", vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(4), trans.columns, rotation=32, ha="right")
    ax.set_yticks(range(4), trans.index)
    ax.set_xlabel("Neighbor state")
    ax.set_title("Shared-neighbor transition fractions", fontweight="bold", pad=7)
    for i in range(4):
        for j in range(4):
            value = trans.iloc[i, j]
            ax.text(j, i, f"{value:.3f}", ha="center", va="center",
                    fontsize=7.1, color="white" if value > .50 else "#1F2937")
    ax.set_xticks(np.arange(-.5, 4, 1), minor=True)
    ax.set_yticks(np.arange(-.5, 4, 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=.8)
    ax.tick_params(which="both", length=0)
    panel(ax, "F")
    rows.append({"figure":"3","panel":"F","type":"neighbor_transition_heatmap",
                 "source":str(SOURCES["trajectory_transition"]),
                 "description":"Source-normalized neighbor transition matrix"})

    ax = fig.add_subplot(gs[1, 4:6])
    display = {
        "Monocyte_program":"Monocyte program", "SPP1_state_no_SPP1":"SPP1 state (excluded)",
        "C1QC_state":"C1QC state", "SPP1_vs_C1QC":"SPP1 minus C1QC",
        "Vascular_remodeling":"Vascular remodeling", "Immunosuppression":"Immunosuppression",
        "MTA_context":"MTA context",
    }
    cohort["label"] = cohort["score"].map(display)
    cohort = cohort.sort_values("median_rho")
    for y, row in enumerate(cohort.itertuples(index=False)):
        sig = row.sign_test_fdr_bh < .05
        color = C["positive"] if sig else C["negative"]
        ax.plot([0, row.median_rho], [y, y], color="#C2C8D0", lw=1.3)
        ax.scatter(row.median_rho, y, s=35 + row.positive_samples * 1.1,
                   facecolor=color if sig else "white", edgecolor=color, lw=1.4, zorder=3)
        ax.text(.94, y, f"{int(row.positive_samples)}/63; q={row.sign_test_fdr_bh:.2g}",
                ha="right", va="center", fontsize=6.9)
    ax.axvline(0, color="#555", ls="--", lw=.9)
    ax.set_yticks(range(len(cohort)), cohort["label"])
    ax.set_xlim(-.54, .96)
    ax.set_xlabel("Median within-sample rho with DPT")
    ax.set_title("Direction consistency across 63 samples", fontweight="bold", pad=7)
    base.clean_axis(ax, grid=True)
    panel(ax, "G")
    rows.append({"figure":"3","panel":"G","type":"cohort_direction_lollipop",
                 "source":str(SOURCES["trajectory_cohort"]),
                 "description":"Seven score directions across 63 eligible samples"})

    ax = fig.add_subplot(gs[1, 6:8])
    assoc_specs = [
        ("Monocyte_program", "Monocyte", C["monocyte"]),
        ("C1QC_state", "C1QC", C["c1qc"]),
        ("SPP1_state_no_SPP1", "SPP1 state", C["spp1_cell"]),
        ("SPP1_vs_C1QC", "SPP1−C1QC", C["mta"]),
    ]
    rng3 = np.random.default_rng(26213)
    for xpos, (score, label, color) in enumerate(assoc_specs):
        values = assoc[assoc["score"].eq(score)]["rho"].to_numpy()
        q1, medv, q3 = np.quantile(values, [.25, .5, .75])
        ax.plot([xpos, xpos], [q1, q3], color=color, lw=5, alpha=.45, solid_capstyle="round")
        ax.scatter(xpos + rng3.uniform(-.13, .13, len(values)), values, s=19,
                   color=color, edgecolor="white", lw=.35, alpha=.62, zorder=3)
        ax.scatter(xpos, medv, s=58, marker="D", color="#1F2937", edgecolor="white", lw=.65, zorder=4)
    ax.axhline(0, color="#555", ls="--", lw=.8)
    ax.set_xticks(range(4), [x[1] for x in assoc_specs], rotation=18, ha="right")
    ax.set_ylabel("Within-sample rho with DPT")
    ax.set_title("Sample-level association distributions", fontweight="bold", pad=7)
    base.clean_axis(ax)
    panel(ax, "H")
    rows.append({"figure":"3","panel":"H","type":"sample_association_distribution",
                 "source":str(SOURCES["trajectory_assoc"]),
                 "description":"Sample-level rho distributions for four key state scores"})

    summary = {"cells": len(cells), "samples": int(cohort["n_samples"].max()),
               "displayed_cells": len(take)}
    return save(fig, FIGURES[3]), summary


def write_legends() -> None:
    LEGENDS.write_text(
        "# Figure legends (v10r3)\n\n"
        "## 图4.1  TYMP–MTA–SPP1巨噬细胞假说：空间支持与患者层面的边界\n\n"
        "**A**，GLMF2双变量空间图，颜色联合编码TYMP-high/MTA-context核心与SPP1+ TAM的切片内秩。"
        "**B**，GLMF2的联合高状态生态位叠加：粉色为TYMP/MTA Q75核心，青色为SPP1-TAM Q75，黄色为核心邻近的SPP1-TAM。"
        "**C**，12例患者中reference-light、巨噬细胞富集、composition/QC校正、5×5尺度的SPP1排除型状态及SPP1相对C1QC效应。"
        "**D**，两类状态结局在reference-dependent/reference-light路线及4×4/5×5尺度上的队列中位效应与四分位区间。"
        "**E**，GLMF1和GLMF2中嘿啶代谢、半胱氨酸/蛋氨酸代谢、精氨酸/脯氨酸代谢及甘油磷脂代谢的逐spot GSVA活性。"
        "**F**，计算界面与五类下游表型在4×4和5×5尺度的队列中位相关。"
        "**G**，两条分析路线及两个尺度下患者间接效应a×b的分布。通路颜色均按单样本单面板缩放，不表示通路间绝对强弱；空间共定位也不证明5′-MTA中介或直接通讯。\n\n"
        "## 图4.2  患者重复的SPP1–vascular niche及其与T-cell balance的反向空间关系\n\n"
        "**A**，代表性患者Z43_T的SPP1–vascular niche与T-cell balance双变量空间图。"
        "**B**，12例TJ-RCC患者中全部预设关系–路线组合的5×5方向校正效应；正值表示符合预设方向。"
        "**C**，5×5主尺度的队列审计，点大小表示方向支持患者数，实心表示全局BH FDR<0.05。"
        "**D**，TYMP–MTA core、SPP1–vascular niche和T-cell balance沿计算型肿瘤界面的患者中位轨迹及四分位区间。"
        "**E**，两条获得队列FDR支持关系在4×4和5×5尺度的患者分布；菱形为中位数，粗线为四分位区间。以上均为关联结构，不证明直接通讯、免疫排斥或治疗耐药。\n\n"
        "## 图4.3  髓系横断面连续谱以C1QC状态增强为主，而非已证明的SPP1谱系\n\n"
        "**A**，从53,794个髓系细胞中抽样展示32,000个细胞的UMAP。"
        "**B**，四类髓系状态的横断面DPT分布。"
        "**C**，16个等细胞数DPT分箱上的样本内标准化状态评分；曲线为跨样本中位数，阴影为四分位区间。"
        "**D**，同一组状态评分在DPT分箱上的热图。"
        "**E**，四类状态的source-normalized共享近邻转移矩阵；该矩阵编码局部图连通性而非谱系方向。"
        "**F**，63个合格样本中七类评分与DPT的方向一致性。DPT是横断面排序，不能替代RNA velocity、纵向取样或谱系示踪。\n",
        encoding="utf-8",
    )


def main() -> None:
    style()
    for path in SOURCES.values():
        if not path.is_file():
            raise FileNotFoundError(path)
    FIGDIR.mkdir(parents=True, exist_ok=True)
    outputs = [*(p.with_suffix(ext) for p in FIGURES.values() for ext in (".png", ".pdf")),
               MANIFEST, LEGENDS, METADATA]
    existing = [str(p) for p in outputs if p.exists()]
    if existing:
        raise FileExistsError("Refusing to overwrite existing outputs:\n" + "\n".join(existing))
    data = load_data()
    rows: list[dict[str, str]] = []
    f1, s1 = build_figure1(rows, data)
    f2, s2 = build_figure2(rows, data)
    f3, s3 = build_figure3(rows, data)
    write_legends()
    pd.DataFrame(rows).to_csv(MANIFEST, sep="\t", index=False)
    metadata = {
        "status":"generated_pending_manual_qc",
        "generated_at":datetime.now().astimezone().isoformat(timespec="seconds"),
        "design":"all panels are direct data visualizations; no workflow, network, evidence-ladder, or interpretation-box panels",
        "outputs":{"figure1":f1, "figure2":f2, "figure3":f3,
                   "legends":str(LEGENDS), "panel_manifest":str(MANIFEST)},
        "summaries":{"figure1":s1, "figure2":s2, "figure3":s3},
        "source_hashes":{name:sha256(path) for name, path in SOURCES.items()},
    }
    METADATA.write_text(json.dumps(metadata, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(metadata, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
