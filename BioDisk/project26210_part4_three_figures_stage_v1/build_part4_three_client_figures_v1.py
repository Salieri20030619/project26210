#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import math
from datetime import datetime
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import to_rgb
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Circle
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


ROOT = Path("/mnt/BioDisk/ludingyu/project26210")
PART4 = ROOT / "results/part4_spatial"
OUT = Path("/mnt/BioDisk/wangmingyuan/project26210_part4_three_figures_stage_v1")
FIGDIR = OUT / "figures"

IMMUNE = PART4 / "tymp_immune_ecology_v2"
STATE = PART4 / "tymp_spp1_state_v1"
VASCULAR = PART4 / "tymp_vascular_niche_v1"
TRAJECTORY = PART4 / "tymp_myeloid_state_ordering_v1"
DIRECTED = ROOT / "external/ccRCC_spatial/TJ_RCC_figshare_24599295/analysis_directed_communication_v1"
REFERENCE = PART4 / "ref/ref_signatures.h5ad"
GLMF2_H5AD = PART4 / "analysis/GLMF2_analyzed.h5ad"

SOURCES = {
    "immune_spot": IMMUNE / "tables/original_spot_scores.tsv.gz",
    "immune_tests": IMMUNE / "tables/original_tymp_immune_spatial_tests.tsv",
    "signature_tests": IMMUNE / "tables/original_tymp_signature_spatial_tests.tsv",
    "state_patient": STATE / "tables/patient_state_associations_v1.tsv",
    "state_cohort": STATE / "tables/cohort_state_direction_summary_v1.tsv",
    "mediation_patient": STATE / "tables/patient_exploratory_mediation_v1.tsv",
    "mediation_cohort": STATE / "tables/cohort_exploratory_mediation_summary_v1.tsv",
    "vascular_spot": VASCULAR / "tables/spot_mechanism_scores_v1.tsv.gz",
    "vascular_blocks": VASCULAR / "tables/sample_spatial_block_tests_v1.tsv",
    "vascular_cohort": VASCULAR / "tables/tjrcc_cohort_direction_summary_v1.tsv",
    "vascular_profiles": VASCULAR / "tables/boundary_distance_profiles_v1.tsv",
    "trajectory_cells": TRAJECTORY / "tables/myeloid_trajectory_cell_scores_v1.tsv.gz",
    "trajectory_assoc": TRAJECTORY / "tables/sample_pseudotime_associations_v1.tsv",
    "trajectory_cohort": TRAJECTORY / "tables/cohort_pseudotime_direction_summary_v1.tsv",
    "trajectory_transition": TRAJECTORY / "tables/trajectory_neighbor_transition_matrix_v1.tsv",
    "trajectory_celltype": TRAJECTORY / "tables/celltype_trajectory_summary_v1.tsv",
    "directed_axis": DIRECTED / "tables/dual_threshold_axis_screen.tsv",
    "reference": REFERENCE,
    "glmf2_h5ad": GLMF2_H5AD,
}

FIGURES = {
    1: FIGDIR / "Figure1_TYMP_MTA_SPP1_evidence_narrowing_v4",
    2: FIGDIR / "Figure2_SPP1_vascular_niche_Tcell_balance_v4",
    3: FIGDIR / "Figure3_myeloid_state_continuum_v4",
}
MANIFEST = OUT / "panel_source_manifest_v4.tsv"
LEGENDS = OUT / "figure_legends_v4.md"
METADATA = OUT / "figure_generation_metadata_v4.json"


COLORS = {
    "tumor": "#C33C54",
    "mta": "#7E57C2",
    "spp1": "#009E73",
    "vascular": "#0072B2",
    "tcell": "#E69F00",
    "c1qc": "#2A9D8F",
    "monocyte": "#4C78A8",
    "other_macro": "#F58518",
    "apoe": "#54A24B",
    "spp1_cell": "#E45756",
    "ink": "#253247",
    "muted": "#6B7280",
    "grid": "#D9DEE7",
    "positive": "#2A9D8F",
    "conditional": "#E9A23B",
    "negative": "#A7AFBD",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def configure_style() -> None:
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 9.2,
        "axes.titlesize": 11.5,
        "axes.labelsize": 9.5,
        "xtick.labelsize": 8.3,
        "ytick.labelsize": 8.3,
        "axes.edgecolor": "#5B6575",
        "axes.linewidth": 0.8,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "savefig.facecolor": "white",
    })


def panel_letter(ax: plt.Axes, letter: str, x: float = -0.08, y: float = 1.04) -> None:
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=18, fontweight="bold",
            color="#111827", va="top", ha="left", zorder=100,
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.92, "pad": 1.0})


def clean_axis(ax: plt.Axes, grid: bool = False) -> None:
    ax.spines[["top", "right"]].set_visible(False)
    if grid:
        ax.grid(axis="x", color=COLORS["grid"], lw=0.65, alpha=0.8)
        ax.set_axisbelow(True)


def add_box(ax: plt.Axes, xy: tuple[float, float], width: float, height: float,
            title: str, subtitle: str, face: str, edge: str, title_color: str | None = None) -> None:
    x, y = xy
    patch = FancyBboxPatch((x, y), width, height, boxstyle="round,pad=0.012,rounding_size=0.025",
                           facecolor=face, edgecolor=edge, lw=1.6)
    ax.add_patch(patch)
    ax.text(x + width / 2, y + height * 0.63, title, ha="center", va="center",
            fontsize=11.4, fontweight="bold", color=title_color or edge)
    ax.text(x + width / 2, y + height * 0.28, subtitle, ha="center", va="center",
            fontsize=8.0, color=COLORS["ink"], linespacing=1.25)


def add_arrow(ax: plt.Axes, start: tuple[float, float], end: tuple[float, float],
              color: str = "#62748A", width: float = 1.6, style: str = "-|>") -> None:
    ax.add_patch(FancyArrowPatch(start, end, arrowstyle=style, mutation_scale=13,
                                 color=color, lw=width, connectionstyle="arc3,rad=0.0"))


def save_figure(fig: plt.Figure, stem: Path) -> dict[str, str | int]:
    png = stem.with_suffix(".png")
    pdf = stem.with_suffix(".pdf")
    fig.savefig(png, dpi=450, bbox_inches="tight", facecolor="white")
    fig.savefig(pdf, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return {
        "png": str(png), "pdf": str(pdf),
        "png_size": png.stat().st_size, "pdf_size": pdf.stat().st_size,
        "png_sha256": sha256(png), "pdf_sha256": sha256(pdf),
    }


def robust_rank(values: pd.Series | np.ndarray) -> np.ndarray:
    series = pd.Series(np.asarray(values, dtype=float))
    return series.rank(method="average", pct=True).to_numpy()


def bivariate_colors(a: np.ndarray, b: np.ndarray,
                     low: str = "#EDF1F6", a_color: str = "#C33C54",
                     b_color: str = "#2E75B6", both: str = "#593B87") -> np.ndarray:
    a = np.clip(np.asarray(a, dtype=float), 0, 1)
    b = np.clip(np.asarray(b, dtype=float), 0, 1)
    c00 = np.array(to_rgb(low)); c10 = np.array(to_rgb(a_color))
    c01 = np.array(to_rgb(b_color)); c11 = np.array(to_rgb(both))
    colors = ((1-a)[:, None]*(1-b)[:, None]*c00 + a[:, None]*(1-b)[:, None]*c10
              + (1-a)[:, None]*b[:, None]*c01 + a[:, None]*b[:, None]*c11)
    return colors


def bivariate_legend(ax: plt.Axes, x_label: str, y_label: str,
                     a_color: str, b_color: str, both: str) -> None:
    inset = inset_axes(ax, width="20%", height="20%", loc="lower right", borderpad=1.0)
    grid = np.linspace(0, 1, 60)
    aa, bb = np.meshgrid(grid, grid)
    rgb = bivariate_colors(aa.ravel(), bb.ravel(), a_color=a_color, b_color=b_color, both=both).reshape(60, 60, 3)
    inset.imshow(rgb, origin="lower", aspect="auto")
    inset.set_xticks([0, 59], ["low", "high"], fontsize=5.8)
    inset.set_yticks([0, 59], ["low", "high"], fontsize=5.8)
    inset.set_xlabel(x_label, fontsize=6.2, labelpad=0)
    inset.set_ylabel(y_label, fontsize=6.2, labelpad=0)
    inset.tick_params(length=1.5, pad=1)


def build_figure1(panel_rows: list[dict[str, str]], data: dict[str, pd.DataFrame]) -> tuple[dict, dict]:
    spot = data["immune_spot"]
    immune = data["immune_tests"]
    signatures = data["signature_tests"]
    state_patient = data["state_patient"]
    state_cohort = data["state_cohort"]
    med_cohort = data["mediation_cohort"]
    directed = data["directed_axis"]

    fig = plt.figure(figsize=(19.2, 10.8), constrained_layout=True)
    gs = fig.add_gridspec(2, 6, height_ratios=[1.0, 1.08], width_ratios=[1.05, 1.05, 1.45, 1.0, 1.0, 1.0])
    fig.suptitle("Spatial evidence progressively narrows the TYMP–MTA–SPP1 macrophage hypothesis",
                 fontsize=19, fontweight="bold", color=COLORS["ink"], y=1.02)

    # A: study logic, not another data plot.
    ax = fig.add_subplot(gs[0, 0:2]); ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    add_box(ax, (0.02, 0.60), 0.27, 0.25, "GLMF1/2", "spatial discovery\n2 sections", "#FCECEF", COLORS["tumor"])
    add_box(ax, (0.365, 0.60), 0.27, 0.25, "TJ-RCC", "patient validation\n12 tumors", "#EAF3FA", COLORS["vascular"])
    add_box(ax, (0.71, 0.60), 0.27, 0.25, "Strict audit", "mediation +\n19 directed axes", "#F1EEF8", COLORS["mta"])
    add_arrow(ax, (0.29, 0.725), (0.365, 0.725)); add_arrow(ax, (0.635, 0.725), (0.71, 0.725))
    add_box(ax, (0.10, 0.13), 0.80, 0.28, "Core question",
            "Does a TYMP-high tumor region align with an MTA-context\nand a distinct SPP1-like macrophage state in human ccRCC?",
            "#F7F9FC", "#7C8798", COLORS["ink"])
    panel_letter(ax, "A", x=-0.02, y=1.01)
    panel_rows.append({"figure":"1","panel":"A","type":"study_design_schematic","source":"formal analysis design","description":"Discovery-to-validation-to-strict-audit framework"})

    # B: one high-information spatial overlay rather than three repeated maps.
    ax = fig.add_subplot(gs[0, 2:4])
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
    frame["x"] = coords[:, 0] * scale; frame["y"] = coords[:, 1] * scale
    tumor_mta = np.minimum(robust_rank(frame["TYMP_High_Malignant_q05"]), robust_rank(frame["score_polyamine_mta_context"]))
    spp1 = robust_rank(frame["SPP1_TAM_q05"])
    colors = bivariate_colors(tumor_mta, spp1, a_color=COLORS["tumor"], b_color=COLORS["spp1"], both="#6B3F87")
    ax.imshow(image, alpha=0.60)
    ax.scatter(frame["x"], frame["y"], c=colors, s=12, lw=0, alpha=0.95, rasterized=True)
    high_core = tumor_mta >= 0.75; high_spp1 = spp1 >= 0.75
    ax.scatter(frame.loc[high_core, "x"], frame.loc[high_core, "y"], s=20, facecolors="none",
               edgecolors="#FFCC33", lw=0.45, alpha=0.75, rasterized=True)
    ax.set(xlim=(0, image.shape[1]), ylim=(image.shape[0], 0)); ax.set_aspect("equal"); ax.axis("off")
    ax.set_title("GLMF2: bivariate tumor/MTA–SPP1 spatial architecture", fontweight="bold")
    bivariate_legend(ax, "TYMP–MTA core", "SPP1+ TAM", COLORS["tumor"], COLORS["spp1"], "#6B3F87")
    ax.text(0.015, 0.02, "Gold rings: TYMP–MTA core ≥ Q75", transform=ax.transAxes, fontsize=7.2,
            bbox={"facecolor":"white","edgecolor":"none","alpha":0.82,"pad":2})
    panel_letter(ax, "B", x=-0.04, y=1.02)
    panel_rows.append({"figure":"1","panel":"B","type":"bivariate_spatial_map","source":f"{SOURCES['glmf2_h5ad']};{SOURCES['immune_spot']}","description":"GLMF2 joint TYMP-MTA and SPP1 TAM spatial encoding"})

    # C: formal GLMF effects across samples and scales.
    ax = fig.add_subplot(gs[0, 4:6])
    effect_parts = []
    relation_specs = [
        ("TYMP-high → MTA context", signatures, {"source":"TYMP-High Malignant","target":"Polyamine-MTA context"}),
        ("SPP1 TAM ↔ MTA context", signatures, {"source":"SPP1+ TAM","target":"Polyamine-MTA context"}),
        ("TYMP-high → SPP1 TAM | absolute", immune, {"metric":"absolute_q05","source":"TYMP-High Malignant","target":"SPP1+ TAM"}),
        ("TYMP-high → SPP1 TAM | relative", immune, {"metric":"relative_abundance","source":"TYMP-High Malignant","target":"SPP1+ TAM"}),
    ]
    for label, source, filters in relation_specs:
        mask = source["scale"].eq("spatial_block_permutation")
        for key, value in filters.items(): mask &= source[key].eq(value)
        part = source[mask & source["grid"].isin(["4x4","5x5"])].copy()
        part["relation_label"] = label; effect_parts.append(part)
    effects = pd.concat(effect_parts, ignore_index=True)
    y_order = [x[0] for x in relation_specs]
    sample_color = {"GLMF1":"#355C7D", "GLMF2":"#D95F59"}; grid_marker={"4x4":"o","5x5":"s"}
    offsets = {("GLMF1","4x4"):-0.18,("GLMF1","5x5"):-0.06,("GLMF2","4x4"):0.06,("GLMF2","5x5"):0.18}
    for _, row in effects.iterrows():
        y = y_order.index(row["relation_label"]) + offsets[(row["sample"], row["grid"])]
        sig = float(row["fdr_bh"]) < 0.05
        ax.scatter(row["rho"], y, s=62, marker=grid_marker[row["grid"]],
                   facecolor=sample_color[row["sample"]] if sig else "white",
                   edgecolor=sample_color[row["sample"]], lw=1.4, zorder=4)
    ax.axvline(0, color="#4B5563", ls="--", lw=1)
    ax.set_yticks(range(len(y_order)), y_order); ax.invert_yaxis(); ax.set_xlim(-0.65, 0.9)
    ax.set_xlabel("Spatial-block Spearman rho"); ax.set_title("Discovery effects across sample and scale", fontweight="bold")
    clean_axis(ax, grid=True)
    handles = [Line2D([0],[0],marker="o",color="none",markerfacecolor=c,markeredgecolor=c,label=s,markersize=7) for s,c in sample_color.items()]
    handles += [Line2D([0],[0],marker=m,color="#555",markerfacecolor="white",label=g,linestyle="none",markersize=7) for g,m in grid_marker.items()]
    handles += [Line2D([0],[0],marker="o",color="#555",markerfacecolor="#555",label="FDR < 0.05",linestyle="none",markersize=7),
                Line2D([0],[0],marker="o",color="#555",markerfacecolor="white",label="FDR ≥ 0.05",linestyle="none",markersize=7)]
    ax.legend(handles=handles, ncol=3, fontsize=7.0, frameon=False, loc="lower right")
    panel_letter(ax, "C", x=-0.06, y=1.02)
    panel_rows.append({"figure":"1","panel":"C","type":"multi_scale_effect_dotplot","source":f"{SOURCES['immune_tests']};{SOURCES['signature_tests']}","description":"Four GLMF relations across two samples and two block grids"})

    # D: patient-wise relation between no-SPP1 state and SPP1-vs-C1QC contrast.
    ax = fig.add_subplot(gs[1, 0:3])
    def patient_effect(outcome: str) -> pd.DataFrame:
        return state_patient[(state_patient["route"].eq("reference_light")) &
                             (state_patient["adjustment"].eq("composition_qc")) &
                             (state_patient["subset"].eq("macrophage_enriched")) &
                             (state_patient["grid"].eq("5x5")) &
                             (state_patient["outcome"].eq(outcome))][["sample","rho"]].rename(columns={"rho":outcome})
    xy = patient_effect("SPP1_state_no_SPP1").merge(patient_effect("SPP1_vs_C1QC_contrast"), on="sample")
    ax.axhspan(ax.get_ylim()[0], 0, color="#EEF3F8", alpha=0.6); ax.axvspan(0, 1, color="#F0F8F5", alpha=0.65)
    ax.scatter(xy["SPP1_state_no_SPP1"], xy["SPP1_vs_C1QC_contrast"], s=62,
               c=np.where(xy["SPP1_state_no_SPP1"]>0, COLORS["positive"], COLORS["negative"]),
               edgecolor="white", lw=0.8, zorder=4)
    for i,row in xy.iterrows():
        dx = 4 if i % 2 == 0 else -4; ha = "left" if dx > 0 else "right"
        ax.annotate(row["sample"], (row["SPP1_state_no_SPP1"],row["SPP1_vs_C1QC_contrast"]),
                    xytext=(dx,4), textcoords="offset points", fontsize=7.0, ha=ha)
    ax.axvline(0,color="#4B5563",ls="--",lw=1); ax.axhline(0,color="#4B5563",ls="--",lw=1)
    ax.set_xlabel("No-SPP1 macrophage-state effect (rho)")
    ax.set_ylabel("SPP1-state minus C1QC-state effect (rho)")
    ax.set_title("Twelve patients: positive state trend does not imply SPP1 dominance", fontweight="bold")
    clean_axis(ax, grid=True)
    ax.text(0.98,0.04,"Reference-light · 5×5 · composition/QC adjusted",transform=ax.transAxes,ha="right",fontsize=7.3,color=COLORS["muted"])
    panel_letter(ax, "D", x=-0.05, y=1.02)
    panel_rows.append({"figure":"1","panel":"D","type":"patient_bivariate_effect_plot","source":str(SOURCES['state_patient']),"description":"No-SPP1 state versus SPP1-minus-C1QC effects across 12 patients"})

    # E: evidence ladder based on formal counts, not decorative repetition.
    ax = fig.add_subplot(gs[1, 3:6]); ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    glmf_mta = effects[effects["relation_label"].eq("TYMP-high → MTA context")]
    glmf_abs = effects[effects["relation_label"].eq("TYMP-high → SPP1 TAM | absolute")]
    ref_light = state_cohort[(state_cohort["route"].eq("reference_light")) &
                             (state_cohort["adjustment"].eq("composition_qc")) &
                             (state_cohort["subset"].eq("macrophage_enriched")) &
                             (state_cohort["outcome"].eq("SPP1_state_no_SPP1"))]
    ref_light_q = float(ref_light["sign_test_fdr_bh"].min())
    ref_light_pos = int(ref_light["positive_patients"].max())
    med_q = float(med_cohort["indirect_a_times_b_sign_fdr_bh"].min())
    directed_pass = int(directed["strict_axis_pass"].sum())
    stages = [
        ("Spatial\ncontext", f"TYMP→MTA\n{int((glmf_mta.fdr_bh<.05).sum())}/4 tests", COLORS["positive"], "supported"),
        ("Local\nco-occupancy", f"TYMP→SPP1 TAM\n{int((glmf_abs.fdr_bh<.05).sum())}/4 tests", COLORS["positive"], "supported"),
        ("Patient state\ntrend", f"{ref_light_pos}/12 positive\nglobal q={ref_light_q:.3f}", COLORS["conditional"], "conditional"),
        ("MTA\nmediation", f"indirect effect\nbest q={med_q:.3f}", COLORS["negative"], "not supported"),
        ("Directed\ncommunication", f"strict axes\n{directed_pass}/19", COLORS["negative"], "not supported"),
    ]
    widths=[0.17,0.17,0.17,0.15,0.15]; xs=[0.01,0.215,0.42,0.625,0.815]
    for idx,((title,subtitle,color,status),x,w) in enumerate(zip(stages,xs,widths)):
        add_box(ax,(x,0.32),w,0.40,title,subtitle,"white",color,color)
        ax.text(x+w/2,0.22,status.upper(),ha="center",va="center",fontsize=7.3,fontweight="bold",color=color)
        if idx < len(stages)-1: add_arrow(ax,(x+w+0.01,0.52),(xs[idx+1]-0.01,0.52),color="#7A8796",width=1.5)
    ax.text(0.5,0.91,"Evidence narrows from spatial compatibility to mechanistic stringency",ha="center",fontsize=11.5,fontweight="bold",color=COLORS["ink"])
    ax.text(0.5,0.06,"Spatial data anchor the human-tissue context; they do not independently prove 5′-MTA mediation or direct signaling.",ha="center",fontsize=8.2,color=COLORS["muted"])
    panel_letter(ax, "E", x=-0.015, y=1.02)
    panel_rows.append({"figure":"1","panel":"E","type":"evidence_narrowing_ladder","source":f"{SOURCES['immune_tests']};{SOURCES['signature_tests']};{SOURCES['state_cohort']};{SOURCES['mediation_cohort']};{SOURCES['directed_axis']}","description":"Formal evidence levels from spatial context to strict communication"})

    summary = {
        "representative_sample":"GLMF2",
        "glmf_tymp_mta_fdr_pass":int((glmf_mta.fdr_bh<.05).sum()),
        "glmf_tymp_spp1_absolute_fdr_pass":int((glmf_abs.fdr_bh<.05).sum()),
        "reference_light_positive_patients":ref_light_pos,
        "reference_light_global_fdr":ref_light_q,
        "best_mediation_fdr":med_q,
        "strict_directed_axes":directed_pass,
    }
    return save_figure(fig, FIGURES[1]), summary


def build_figure2(panel_rows: list[dict[str, str]], data: dict[str, pd.DataFrame]) -> tuple[dict, dict]:
    spot = data["vascular_spot"]
    blocks = data["vascular_blocks"]
    cohort = data["vascular_cohort"]
    profiles = data["vascular_profiles"]
    representative = "Z43_T"
    fig = plt.figure(figsize=(19.2, 11.0))
    gs = fig.add_gridspec(2, 6, height_ratios=[1.0, 1.05], width_ratios=[1.15,1.15,1.2,1.0,1.0,1.0])
    fig.subplots_adjust(left=0.055,right=0.985,bottom=0.065,top=0.885,wspace=0.55,hspace=0.18)
    fig.suptitle("A patient-recurrent SPP1–vascular niche is inversely aligned with T-cell balance",
                 fontsize=19, fontweight="bold", color=COLORS["ink"], y=0.995)

    # A: construction of the tested variables.
    ax=fig.add_subplot(gs[0,0:2]); ax.set_xlim(0,1); ax.set_ylim(0,1); ax.axis("off")
    add_box(ax,(0.03,0.62),0.25,0.22,"SPP1 TAM","spot percentile","#EAF8F3",COLORS["spp1"])
    add_box(ax,(0.37,0.70),0.25,0.18,"Endothelial","spot percentile","#EAF3FA",COLORS["vascular"])
    add_box(ax,(0.37,0.47),0.25,0.18,"Pericyte","spot percentile","#EAF3FA",COLORS["vascular"])
    add_box(ax,(0.70,0.57),0.27,0.25,"SPP1–vascular\nniche","composite spatial score","#F1EEF8",COLORS["mta"])
    add_arrow(ax,(0.28,0.73),(0.70,0.69)); add_arrow(ax,(0.62,0.79),(0.70,0.71)); add_arrow(ax,(0.62,0.56),(0.70,0.65))
    add_box(ax,(0.15,0.12),0.28,0.22,"TYMP–MTA core","reference-dependent\nand reference-light","#FCECEF",COLORS["tumor"])
    add_box(ax,(0.59,0.12),0.28,0.22,"T-cell balance","T-cell abundance\nrelative to immune context","#FFF4DF",COLORS["tcell"])
    add_arrow(ax,(0.43,0.23),(0.59,0.23),color="#9B7FB5",style="<->")
    panel_letter(ax,"A",x=-0.02,y=1.01)
    panel_rows.append({"figure":"2","panel":"A","type":"score_construction_schematic","source":str(SOURCES['vascular_spot']),"description":"Construction and relationships of vascular-niche variables"})

    # B: representative bivariate patient spatial map.
    ax=fig.add_subplot(gs[0,2:4])
    rep=spot[spot["sample"].eq(representative)].copy()
    niche=robust_rank(rep["spp1_vascular_niche"]); tbal=robust_rank(rep["T_cell_balance"])
    colors=bivariate_colors(niche,tbal,a_color=COLORS["spp1"],b_color=COLORS["tcell"],both="#66518C")
    ax.scatter(rep["spatial_x"],rep["spatial_y"],c=colors,s=8.5,lw=0,rasterized=True)
    high=rep["computational_boundary_spot"].astype(bool).to_numpy()
    if high.any():
        ax.scatter(rep.loc[high,"spatial_x"],rep.loc[high,"spatial_y"],s=14,facecolors="none",edgecolors="#222",lw=0.35,alpha=0.65,rasterized=True)
    ax.set_aspect("equal"); ax.invert_yaxis(); ax.axis("off")
    ax.set_title(f"{representative}: niche and T-cell balance\noccupy opposing regions",fontweight="bold",fontsize=10.4,pad=2)
    bivariate_legend(ax,"SPP1–vascular niche","T-cell balance",COLORS["spp1"],COLORS["tcell"],"#66518C")
    ax.text(0.02,0.02,"Representative selected by maximin 5×5 support\nfor both prespecified relations",transform=ax.transAxes,fontsize=7.0,
            bbox={"facecolor":"white","edgecolor":"none","alpha":0.85,"pad":2})
    panel_letter(ax,"B",x=-0.04,y=1.02)
    panel_rows.append({"figure":"2","panel":"B","type":"bivariate_patient_spatial_map","source":str(SOURCES['vascular_spot']),"description":"Z43_T SPP1-vascular niche and T-cell balance bivariate spatial map"})

    # C: patient dumbbells across 4x4 and 5x5 for the two supported relations.
    parent=fig.add_subplot(gs[0,4:6]); parent.axis("off")
    sub=gs[0,4:6].subgridspec(1,2,wspace=0.15)
    relations=[
        ("SPP1_vascular_niche_to_T_cell_balance","route_independent","Niche → T-cell balance\n(direction-aligned)"),
        ("TYMP_MTA_core_to_SPP1_vascular_niche","reference_dependent","TYMP–MTA → niche\n(direction-aligned)"),
    ]
    selected_rows=[]
    for j,(relation,route,title) in enumerate(relations):
        axc=fig.add_subplot(sub[0,j])
        d=blocks[(blocks["cohort"].eq("TJ-RCC"))&(blocks["relation"].eq(relation))&(blocks["route"].eq(route))].copy()
        wide=d.pivot(index="sample",columns="grid",values="direction_aligned_rho")
        wide["mean"]=wide[["4x4","5x5"]].mean(axis=1); wide=wide.sort_values("mean")
        for y,(sample,row) in enumerate(wide.iterrows()):
            axc.plot([row["4x4"],row["5x5"]],[y,y],color="#B8C0CC",lw=1.4,zorder=1)
            axc.scatter(row["4x4"],y,s=32,color="#88A9C3",marker="o",zorder=3)
            axc.scatter(row["5x5"],y,s=38,color=COLORS["mta"],marker="s",zorder=3)
        axc.axvline(0,color="#555",ls="--",lw=0.9); axc.set_yticks(range(len(wide)),wide.index if j==0 else [])
        axc.set_xlabel("Aligned rho"); axc.set_title(title,fontsize=8.6,pad=2); clean_axis(axc,grid=True)
        selected_rows.append(d)
    parent.legend(handles=[Line2D([0],[0],marker="o",color="none",markerfacecolor="#88A9C3",label="4×4",markersize=7),
                           Line2D([0],[0],marker="s",color="none",markerfacecolor=COLORS["mta"],label="5×5",markersize=7)],
                  loc="lower center",bbox_to_anchor=(0.5,-0.09),ncol=2,frameon=False,fontsize=7.5)
    panel_letter(parent,"C",x=-0.06,y=1.02)
    panel_rows.append({"figure":"2","panel":"C","type":"multi_scale_patient_dumbbell","source":str(SOURCES['vascular_blocks']),"description":"Two supported relations across 12 patients and two block grids"})

    # D: all nine prespecified 5x5 relation-route combinations in one compact audit.
    ax=fig.add_subplot(gs[1,0:2])
    c5=cohort[cohort["grid"].eq("5x5")].copy()
    label_map={
        "SPP1_TAM_to_vascular_compartment":"SPP1 TAM → vascular compartment",
        "SPP1_vascular_niche_to_T_cell_balance":"Niche → T-cell balance",
        "SPP1_vascular_niche_to_cytotoxicity":"Niche → cytotoxicity",
        "SPP1_vascular_niche_to_hypoxia_angiogenesis":"Niche → hypoxia/angiogenesis",
        "SPP1_vascular_niche_to_invasion_EMT":"Niche → invasion/EMT",
        "TYMP_MTA_core_to_SPP1_monocyte_balance":"TYMP–MTA → SPP1/Mono balance",
        "TYMP_MTA_core_to_SPP1_vascular_niche":"TYMP–MTA → niche",
    }
    route_tag={"route_independent":"","reference_light":" [ref-light]","reference_dependent":" [ref-dep]"}
    c5["label"]=[label_map[r]+route_tag[rt] for r,rt in zip(c5["relation"],c5["route"])]
    c5=c5.sort_values("median_direction_aligned_rho")
    route_color={"route_independent":COLORS["vascular"],"reference_light":"#D55E72","reference_dependent":COLORS["mta"]}
    for y,row in enumerate(c5.itertuples(index=False)):
        ax.plot([0,row.median_direction_aligned_rho],[y,y],color="#C6CDD7",lw=1.2)
        sig=row.sign_test_fdr_bh_global<0.05
        ax.scatter(row.median_direction_aligned_rho,y,s=35+row.direction_supported_patients*5,
                   color=route_color[row.route] if sig else "white",edgecolor=route_color[row.route],lw=1.5,zorder=3)
        ax.text(.64,y,f"{row.direction_supported_patients}/12; q={row.sign_test_fdr_bh_global:.3g}",
                ha="right",va="center",fontsize=6.9)
    ax.axvline(0,color="#555",ls="--",lw=0.9); ax.set_yticks(range(len(c5)),c5["label"])
    ax.set_xlim(-.44,.66)
    ax.set_xlabel("Median direction-aligned patient rho (5×5)"); ax.set_title("Prespecified cohort audit",fontweight="bold")
    clean_axis(ax,grid=True); panel_letter(ax,"D",x=-0.06,y=1.02)
    panel_rows.append({"figure":"2","panel":"D","type":"cohort_relation_audit_dotplot","source":str(SOURCES['vascular_cohort']),"description":"All nine 5x5 relation-route combinations with support counts and FDR"})

    # E: cohort median distance profiles across the computational interface.
    ax=fig.add_subplot(gs[1,2:4])
    p=profiles[profiles["cohort"].eq("TJ-RCC")].copy()
    score_specs=[("TYMP_MTA_core_reference_dependent","TYMP–MTA core",COLORS["tumor"]),
                 ("spp1_vascular_niche","SPP1–vascular niche",COLORS["spp1"]),
                 ("T_cell_balance","T-cell balance",COLORS["tcell"])]
    zones=p[["zone_order","distance_zone"]].drop_duplicates().sort_values("zone_order")
    x=zones["zone_order"].to_numpy(); labels=zones["distance_zone"].str.replace("Peripheral ","",regex=False).str.replace("Distal ","",regex=False)
    for score,label,color in score_specs:
        subp=p[p["score"].eq(score)]
        agg=subp.groupby("zone_order")["median_z"].agg(median="median",q25=lambda v:v.quantile(.25),q75=lambda v:v.quantile(.75)).reindex(x)
        ax.plot(x,agg["median"],marker="o",lw=2.2,color=color,label=label)
        ax.fill_between(x,agg["q25"],agg["q75"],color=color,alpha=.15,lw=0)
    ax.axhline(0,color="#555",ls="--",lw=.8); ax.set_xticks(x,labels,rotation=24,ha="right")
    ax.set_ylabel("Across-patient median of within-sample z-score")
    ax.set_title("Opposing profiles across a computational tumor interface",fontweight="bold")
    ax.legend(frameon=False,fontsize=7.6); clean_axis(ax,grid=True)
    ax.text(.01,.02,"Descriptive computational boundary; not a pathology-defined tumor–NAT edge.",transform=ax.transAxes,fontsize=6.9,color=COLORS["muted"])
    panel_letter(ax,"E",x=-0.06,y=1.02)
    panel_rows.append({"figure":"2","panel":"E","type":"cohort_distance_profile","source":str(SOURCES['vascular_profiles']),"description":"Three scores across computational distance zones, summarized across 12 patients"})

    # F: final association network with supported and unsupported branches.
    ax=fig.add_subplot(gs[1,4:6]); ax.set_xlim(0,1); ax.set_ylim(0,1); ax.axis("off")
    nodes={"core":(0.16,0.57,COLORS["tumor"],"TYMP–MTA\ncore"),"niche":(0.50,0.57,COLORS["spp1"],"SPP1–vascular\nniche"),
           "tcell":(0.84,0.72,COLORS["tcell"],"T-cell\nbalance"),"other":(0.84,0.35,COLORS["negative"],"Other tested\nendpoints")}
    for _,(x0,y0,color,label) in nodes.items():
        ax.add_patch(Circle((x0,y0),0.11,facecolor="white",edgecolor=color,lw=2.4)); ax.text(x0,y0,label,ha="center",va="center",fontsize=9,fontweight="bold",color=color)
    add_arrow(ax,(.27,.57),(.39,.57),color=COLORS["positive"],width=3)
    add_arrow(ax,(.61,.61),(.73,.68),color=COLORS["positive"],width=3)
    add_arrow(ax,(.61,.52),(.73,.39),color="#9AA3AF",width=1.4)
    ax.text(.33,.64,"+ 12/12\nq=0.00220",ha="center",fontsize=7.5,color=COLORS["positive"],fontweight="bold")
    ax.text(.68,.72,"− 12/12\nq=0.00220",ha="center",fontsize=7.5,color=COLORS["positive"],fontweight="bold")
    ax.text(.69,.41,"heterogeneous",ha="center",fontsize=7.2,color=COLORS["muted"])
    ax.text(.5,.93,"Association network supported by patient-level direction tests",ha="center",fontsize=11,fontweight="bold",color=COLORS["ink"])
    ax.text(.5,.10,"Edges summarize spatial associations, not direct signaling or treatment resistance.",ha="center",fontsize=7.4,color=COLORS["muted"])
    panel_letter(ax,"F",x=-0.04,y=1.02)
    panel_rows.append({"figure":"2","panel":"F","type":"evidence_network","source":str(SOURCES['vascular_cohort']),"description":"Supported and heterogeneous vascular-niche associations"})

    selected=pd.concat(selected_rows)
    supported = cohort.loc[cohort["sign_test_fdr_bh_global"] < .05, ["relation", "route"]].drop_duplicates()
    summary={"representative_patient":representative,"selection_rule":"maximum minimum 5x5 direction-aligned rho across the two supported relations",
             "supported_relations":len(supported),"patient_rows":len(selected)}
    return save_figure(fig,FIGURES[2]),summary


def build_figure3(panel_rows: list[dict[str, str]], data: dict[str, pd.DataFrame]) -> tuple[dict, dict]:
    cells=data["trajectory_cells"].copy(); assoc=data["trajectory_assoc"]; cohort=data["trajectory_cohort"]
    transition=data["trajectory_transition"]; celltype=data["trajectory_celltype"]
    backed=ad.read_h5ad(REFERENCE,backed="r")
    try:
        indexer=pd.Index(backed.obs_names).get_indexer(cells["cell_id"])
        if (indexer<0).any(): raise RuntimeError("Trajectory cell IDs missing from reference")
        coords=np.asarray(backed.obsm["X_umap"])[indexer]
    finally:
        backed.file.close()
    fig=plt.figure(figsize=(19.2,11.0),constrained_layout=True)
    gs=fig.add_gridspec(2,6,height_ratios=[1.0,1.03],width_ratios=[1.05,1.05,1.15,1.1,1.0,1.0])
    fig.suptitle("Myeloid cells form a cross-sectional continuum dominated by C1QC-state gain, not a proven SPP1 lineage",
                 fontsize=18.5,fontweight="bold",color=COLORS["ink"],y=1.02)
    rng=np.random.default_rng(20260805); take=np.sort(rng.choice(len(cells),size=min(32000,len(cells)),replace=False))
    type_colors={"Monocyte":COLORS["monocyte"],"Other Macrophage":COLORS["other_macro"],
                 "APOE+C1QC+ TAM":COLORS["apoe"],"SPP1+ TAM":COLORS["spp1_cell"]}

    # A: one UMAP anchor only.
    ax=fig.add_subplot(gs[0,0:2])
    for label,color in type_colors.items():
        mask=cells["trajectory_celltype"].to_numpy()[take]==label
        ax.scatter(coords[take][mask,0],coords[take][mask,1],s=3.0,color=color,alpha=.58,lw=0,rasterized=True,label=label)
    ax.set_xticks([]);ax.set_yticks([]);ax.set_xlabel("UMAP1");ax.set_ylabel("UMAP2");ax.set_title("Myeloid state map (32,000 cells displayed)",fontweight="bold")
    ax.legend(frameon=False,fontsize=7.4,ncol=2,loc="lower left"); clean_axis(ax); panel_letter(ax,"A",x=-.05,y=1.02)
    panel_rows.append({"figure":"3","panel":"A","type":"single_cell_umap","source":f"{SOURCES['reference']};{SOURCES['trajectory_cells']}","description":"Four myeloid states on reference UMAP"})

    # B: DPT ridgeline distributions by cell type.
    ax=fig.add_subplot(gs[0,2:4])
    bins=np.linspace(cells["dpt_pseudotime"].min(),cells["dpt_pseudotime"].max(),55)
    centers=(bins[:-1]+bins[1:])/2
    for offset,(label,color) in enumerate(type_colors.items()):
        vals=cells.loc[cells["trajectory_celltype"].eq(label),"dpt_pseudotime"].to_numpy()
        hist,_=np.histogram(vals,bins=bins,density=True); hist=hist/(hist.max() or 1)*.78
        ax.fill_between(centers,offset,offset+hist,color=color,alpha=.68,lw=0)
        ax.plot(centers,offset+hist,color=color,lw=1.3)
        ax.text(bins[0],offset+.08,label,fontsize=8.0,color=COLORS["ink"],va="bottom")
    ax.set_yticks([]);ax.set_xlabel("Cross-sectional DPT");ax.set_title("State occupancy shifts along the DPT coordinate",fontweight="bold")
    clean_axis(ax,grid=True);panel_letter(ax,"B",x=-.05,y=1.02)
    panel_rows.append({"figure":"3","panel":"B","type":"pseudotime_ridgeline","source":str(SOURCES['trajectory_cells']),"description":"DPT distributions for four myeloid states"})

    # C: sample-balanced smoothed score profiles.
    ax=fig.add_subplot(gs[0,4:6])
    n_bins=16
    cells["dpt_bin"]=pd.qcut(cells["dpt_pseudotime"],q=n_bins,labels=False,duplicates="drop")
    score_specs=[("Monocyte_program","Monocyte program",COLORS["monocyte"]),
                 ("C1QC_state","C1QC state",COLORS["c1qc"]),
                 ("SPP1_state_no_SPP1","SPP1 state (SPP1 excluded)",COLORS["spp1_cell"]),
                 ("SPP1_vs_C1QC","SPP1 minus C1QC",COLORS["mta"])]
    centers=cells.groupby("dpt_bin",observed=True)["dpt_pseudotime"].median().reindex(range(n_bins)).to_numpy()
    for score,label,color in score_specs:
        group_mean = cells.groupby("sample_id")[score].transform("mean")
        group_sd = cells.groupby("sample_id")[score].transform("std").replace(0, np.nan)
        cells["_score_z"] = (cells[score] - group_mean) / group_sd
        grouped=cells.groupby(["sample_id","dpt_bin"],observed=True)["_score_z"].agg(["median","size"]).reset_index()
        grouped=grouped[grouped["size"]>=10]
        agg=grouped.groupby("dpt_bin")["median"].agg(median="median",q25=lambda x:x.quantile(.25),q75=lambda x:x.quantile(.75)).reindex(range(n_bins))
        ax.plot(centers,agg["median"],color=color,lw=2,label=label)
        ax.fill_between(centers,agg["q25"],agg["q75"],color=color,alpha=.12,lw=0)
    ax.axhline(0,color="#555",ls="--",lw=.8);ax.set_xlabel("Cross-sectional DPT");ax.set_ylabel("Median within-sample z-score")
    ax.set_title("C1QC-state rises more consistently than SPP1-state dominance",fontweight="bold")
    ax.legend(frameon=False,fontsize=7.1);clean_axis(ax,grid=True);panel_letter(ax,"C",x=-.05,y=1.02)
    panel_rows.append({"figure":"3","panel":"C","type":"sample_balanced_state_curves","source":str(SOURCES['trajectory_cells']),"description":"Within-sample binned median scores summarized across samples"})

    # D: kNN graph connectivity network (not lineage arrows).
    ax=fig.add_subplot(gs[1,0:2]);ax.set_xlim(0,1);ax.set_ylim(0,1);ax.axis("off")
    positions={"Monocyte":(.16,.55),"Other Macrophage":(.48,.55),"APOE+C1QC+ TAM":(.78,.75),"SPP1+ TAM":(.78,.32)}
    counts=celltype.set_index("trajectory_celltype")["n_cells"].to_dict()
    cross=transition[transition["source"]!=transition["target"]]
    maxw=float(cross["edge_weight"].max())
    done=set()
    for row in cross.itertuples(index=False):
        pair=tuple(sorted([row.source,row.target]))
        if pair in done: continue
        done.add(pair); reverse=cross[(cross["source"].eq(row.target))&(cross["target"].eq(row.source))]
        weight=float(np.mean([row.edge_weight,*reverse["edge_weight"].tolist()]))
        p1=positions[row.source];p2=positions[row.target];lw=.7+6*np.sqrt(weight/maxw)
        ax.plot([p1[0],p2[0]],[p1[1],p2[1]],color="#8492A6",lw=lw,alpha=.45,zorder=1)
    for label,(x,y) in positions.items():
        color=type_colors[label];radius=.07+.06*np.sqrt(counts[label]/max(counts.values()))
        ax.add_patch(Circle((x,y),radius,facecolor=color,edgecolor="white",lw=2,zorder=3,alpha=.9))
        ax.text(x,y,label.replace(" ","\n",1),ha="center",va="center",fontsize=8,color="white",fontweight="bold",zorder=4)
        ax.text(x,y-radius-.035,f"n={counts[label]:,}",ha="center",fontsize=6.8,color=COLORS["muted"])
    ax.set_title("Shared-neighbor connectivity among annotated states",fontweight="bold")
    ax.text(.5,.05,"Line width = cross-state graph weight; no lineage direction is encoded.",ha="center",fontsize=7.2,color=COLORS["muted"])
    panel_letter(ax,"D",x=-.02,y=1.02)
    panel_rows.append({"figure":"3","panel":"D","type":"knn_connectivity_network","source":f"{SOURCES['trajectory_transition']};{SOURCES['trajectory_celltype']}","description":"Cross-state neighbor graph connectivity"})

    # E: 63-sample direction summary.
    ax=fig.add_subplot(gs[1,2:4])
    display={"Monocyte_program":"Monocyte program","SPP1_state_no_SPP1":"SPP1 state (SPP1 excluded)","C1QC_state":"C1QC state",
             "SPP1_vs_C1QC":"SPP1 minus C1QC","Vascular_remodeling":"Vascular remodeling","Immunosuppression":"Immunosuppression","MTA_context":"MTA context"}
    c=cohort.copy();c["label"]=c["score"].map(display);c=c.sort_values("median_rho")
    maxlog=max(-np.log10(c["sign_test_fdr_bh"].clip(lower=1e-12)))
    for y,row in enumerate(c.itertuples(index=False)):
        sig=row.sign_test_fdr_bh<.05;color=COLORS["positive"] if sig else COLORS["negative"]
        ax.plot([0,row.median_rho],[y,y],color="#C2C8D0",lw=1.4)
        ax.scatter(row.median_rho,y,s=35+row.positive_samples*1.2,facecolor=color if sig else "white",edgecolor=color,lw=1.5,zorder=3)
        ax.text(1.01,y,f"{int(row.positive_samples)}/63 +; q={row.sign_test_fdr_bh:.2g}",transform=ax.get_yaxis_transform(),va="center",fontsize=7.0)
    ax.axvline(0,color="#555",ls="--",lw=.9);ax.set_yticks(range(len(c)),c["label"]);ax.set_xlabel("Median within-sample rho with DPT")
    ax.set_title("Direction consistency across 63 eligible samples",fontweight="bold");clean_axis(ax,grid=True);panel_letter(ax,"E",x=-.06,y=1.02)
    panel_rows.append({"figure":"3","panel":"E","type":"cohort_direction_lollipop","source":str(SOURCES['trajectory_cohort']),"description":"Seven state scores summarized across 63 eligible samples"})

    # F: interpretation boundary.
    ax=fig.add_subplot(gs[1,4:6]);ax.set_xlim(0,1);ax.set_ylim(0,1);ax.axis("off")
    add_box(ax,(.05,.67),.90,.22,"Observed continuum","Monocyte program decreases; C1QC state rises;\nSPP1-excluded state rises weakly.","#EAF8F3",COLORS["positive"])
    add_arrow(ax,(.5,.66),(.5,.55),color="#6B7280")
    add_box(ax,(.05,.34),.90,.22,"Robust interpretation","The late DPT region is relatively C1QC-skewed,\nnot a universal SPP1-dominant endpoint.","#FFF4DF",COLORS["conditional"])
    add_arrow(ax,(.5,.33),(.5,.22),color="#6B7280")
    add_box(ax,(.05,.02),.90,.18,"Not established","No RNA velocity, lineage tracing, longitudinal sampling,\nor monocyte→SPP1+ TAM conversion proof.","#F3F4F6",COLORS["negative"])
    panel_letter(ax,"F",x=-.03,y=1.02)
    panel_rows.append({"figure":"3","panel":"F","type":"interpretation_boundary","source":f"{SOURCES['trajectory_cohort']};{SOURCES['trajectory_transition']}","description":"Observed, supported, and prohibited interpretations"})

    summary={"cells":len(cells),"samples":int(cohort["n_samples"].max()),"displayed_cells":len(take),
             "c1qc_positive_samples":int(cohort.loc[cohort.score.eq('C1QC_state'),'positive_samples'].iloc[0]),
             "spp1_positive_samples":int(cohort.loc[cohort.score.eq('SPP1_state_no_SPP1'),'positive_samples'].iloc[0])}
    return save_figure(fig,FIGURES[3]),summary


def write_legends() -> None:
    text = """# Figure legends

## Figure 1. Spatial evidence progressively narrows the TYMP–MTA–SPP1 macrophage hypothesis

**A**, Analysis design linking the two-section GLMF discovery layer, the 12-patient TJ-RCC validation layer, and strict mediation/directed-communication audits. **B**, GLMF2 bivariate spatial map. Color jointly encodes the within-section ranks of the TYMP-high/MTA-context core and SPP1+ TAM abundance; gold rings mark spots in the upper quartile of the joint TYMP–MTA score. **C**, Prespecified GLMF spatial-block effects for four relations across GLMF1/GLMF2 and 4×4/5×5 grids. Filled symbols indicate within-analysis FDR < 0.05. **D**, Patient-level reference-light effects comparing the SPP1-excluded macrophage-state score with the SPP1-state minus C1QC-state contrast in macrophage-enriched spots (5×5 grid, composition/QC adjusted). **E**, Evidence ladder showing that spatial compatibility and local co-occupancy do not remain confirmatory at the stricter mediation or directed-communication levels. Spatial results are associations and do not independently prove 5′-MTA mediation or direct signaling.

## Figure 2. A patient-recurrent SPP1–vascular niche is inversely aligned with T-cell balance

**A**, Construction of the SPP1–vascular niche and related tested variables. **B**, Representative Z43_T bivariate spatial map jointly encoding the SPP1–vascular niche and T-cell balance; the representative was selected by the maximum of the minimum 5×5 direction-aligned effect across the two prespecified supported relations. **C**, Patient-level direction-aligned effects at 4×4 and 5×5 scales for SPP1–vascular niche versus T-cell balance and reference-dependent TYMP–MTA core versus SPP1–vascular niche. **D**, Cohort audit of all nine prespecified 5×5 relation–route combinations; point size denotes direction-supported patients, fill denotes global FDR < 0.05, and annotations show support count and FDR. **E**, Across-patient median and interquartile range of three scores across a computational tumor-interface distance. This boundary is expression-derived and is not a pathology-defined tumor–NAT edge. **F**, Patient-level association network. Edge signs, support counts and FDR summarize spatial associations and do not demonstrate direct signaling, immune exclusion or treatment resistance.

## Figure 3. Myeloid cells form a cross-sectional continuum dominated by C1QC-state gain, not a proven SPP1 lineage

**A**, UMAP of 32,000 displayed cells sampled from 53,794 myeloid cells and colored by the four analysis states. **B**, Cross-sectional DPT distributions by annotated state. **C**, Sample-balanced state-score profiles across 16 equal-cell-number DPT bins. Scores were standardized within sample; medians were then calculated within sample and bin (minimum 10 cells) and summarized across samples. Ribbons show the interquartile range. **D**, Shared-neighbor graph connectivity among annotated states; node size reflects cell number and edge width reflects cross-state graph weight. No lineage direction is encoded. **E**, Median within-sample Spearman association with DPT across 63 eligible samples; point size denotes positive-sample count and filled points denote BH FDR < 0.05. **F**, Interpretation boundary distinguishing the observed continuum from unsupported lineage claims. DPT is cross-sectional and does not replace RNA velocity, longitudinal sampling or lineage tracing.
"""
    LEGENDS.write_text(text, encoding="utf-8")


def main() -> None:
    configure_style()
    for path in SOURCES.values():
        if not path.is_file(): raise FileNotFoundError(path)
    OUT.mkdir(parents=True, exist_ok=True)
    FIGDIR.mkdir(parents=True, exist_ok=True)
    expected_outputs = [
        *(path.with_suffix(ext) for path in FIGURES.values() for ext in (".png", ".pdf")),
        MANIFEST, LEGENDS, METADATA,
    ]
    existing_outputs = [str(path) for path in expected_outputs if path.exists()]
    if existing_outputs:
        raise FileExistsError("Refusing to overwrite existing outputs:\n" + "\n".join(existing_outputs))
    data = {
        "immune_spot":pd.read_csv(SOURCES["immune_spot"],sep="\t"),
        "immune_tests":pd.read_csv(SOURCES["immune_tests"],sep="\t"),
        "signature_tests":pd.read_csv(SOURCES["signature_tests"],sep="\t"),
        "state_patient":pd.read_csv(SOURCES["state_patient"],sep="\t"),
        "state_cohort":pd.read_csv(SOURCES["state_cohort"],sep="\t"),
        "mediation_cohort":pd.read_csv(SOURCES["mediation_cohort"],sep="\t"),
        "vascular_spot":pd.read_csv(SOURCES["vascular_spot"],sep="\t"),
        "vascular_blocks":pd.read_csv(SOURCES["vascular_blocks"],sep="\t"),
        "vascular_cohort":pd.read_csv(SOURCES["vascular_cohort"],sep="\t"),
        "vascular_profiles":pd.read_csv(SOURCES["vascular_profiles"],sep="\t"),
        "trajectory_cells":pd.read_csv(SOURCES["trajectory_cells"],sep="\t"),
        "trajectory_assoc":pd.read_csv(SOURCES["trajectory_assoc"],sep="\t"),
        "trajectory_cohort":pd.read_csv(SOURCES["trajectory_cohort"],sep="\t"),
        "trajectory_transition":pd.read_csv(SOURCES["trajectory_transition"],sep="\t"),
        "trajectory_celltype":pd.read_csv(SOURCES["trajectory_celltype"],sep="\t"),
        "directed_axis":pd.read_csv(SOURCES["directed_axis"],sep="\t"),
    }
    panel_rows: list[dict[str,str]]=[]
    files1,summary1=build_figure1(panel_rows,data)
    files2,summary2=build_figure2(panel_rows,data)
    files3,summary3=build_figure3(panel_rows,data)
    write_legends()
    pd.DataFrame(panel_rows).to_csv(MANIFEST,sep="\t",index=False)
    metadata={
        "status":"generated_pending_manual_qc",
        "generated_at":datetime.now().astimezone().isoformat(timespec="seconds"),
        "design_principle":"one biological theme per figure; minimize repeated plot types within a figure",
        "outputs":{"figure1":files1,"figure2":files2,"figure3":files3,"legends":str(LEGENDS),"panel_manifest":str(MANIFEST)},
        "summaries":{"figure1":summary1,"figure2":summary2,"figure3":summary3},
        "source_hashes":{name:sha256(path) for name,path in SOURCES.items()},
        "inference_boundary":"Spatial association and cross-sectional DPT do not establish direct signaling, metabolic mediation, lineage conversion, treatment resistance, or causality.",
    }
    METADATA.write_text(json.dumps(metadata,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(metadata,ensure_ascii=False,indent=2))


if __name__ == "__main__":
    main()
