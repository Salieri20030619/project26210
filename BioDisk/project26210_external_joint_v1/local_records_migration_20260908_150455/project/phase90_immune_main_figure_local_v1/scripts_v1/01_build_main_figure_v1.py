from __future__ import annotations

import gzip
import hashlib
import json
import math
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import matplotlib.patches as patches
import networkx as nx
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Wedge


SEED = 26210
RNG = np.random.default_rng(SEED)

COLORS = {
    "tumor": "#E99A2E",
    "tymp": "#C8346D",
    "monocyte": "#1E9C89",
    "spp1": "#C8445B",
    "c1qc": "#6A4C93",
    "other_mac": "#64869B",
    "tcell": "#3A78B5",
    "nk": "#55A7C6",
    "endothelial": "#55A868",
    "bcell": "#4C72B0",
    "stromal": "#CFA349",
    "mast": "#B07AA1",
    "unresolved": "#B8B8B8",
    "ink": "#20252B",
    "muted": "#66717C",
    "grid": "#D9DEE3",
    "paper": "#FFFFFF",
    "soft": "#F4F6F7",
}

STATE_COLORS = {
    "Monocyte": COLORS["monocyte"],
    "SPP1-state Macrophage": COLORS["spp1"],
    "C1QC-state Macrophage": COLORS["c1qc"],
    "Other Macrophage": COLORS["other_mac"],
}

BROAD_COLORS = {
    "Tumor / epithelial": COLORS["tumor"],
    "Myeloid": COLORS["monocyte"],
    "T cell": COLORS["tcell"],
    "NK cell": COLORS["nk"],
    "Endothelial": COLORS["endothelial"],
    "B / plasma": COLORS["bcell"],
    "Stromal": COLORS["stromal"],
    "Mast": COLORS["mast"],
    "Other": COLORS["unresolved"],
}

DIVERGING = LinearSegmentedColormap.from_list(
    "teal_white_magenta", ["#0D8C82", "#E8F5F3", "#FFFFFF", "#F7E4EC", "#B72C63"], N=256
)
TYMP_CMAP = LinearSegmentedColormap.from_list(
    "tymp", ["#E8E8E8", "#F7D6E2", "#D75A85", "#8E1649"], N=256
)
F1_CMAP = LinearSegmentedColormap.from_list(
    "f1", ["#E7E7E7", "#F6D49A", "#E99A2E", "#A84D0D"], N=256
)
SPATIAL_CMAP = LinearSegmentedColormap.from_list(
    "spatial", ["#F1F1F1", "#D9D7F2", "#7C6BB4", "#C02F67", "#F0A43A"], N=256
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def clean_axis(ax: mpl.axes.Axes) -> None:
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)


def panel_heading(ax: mpl.axes.Axes, letter: str, title: str, subtitle: str | None = None) -> None:
    ax.text(
        -0.015,
        1.045,
        letter,
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=13.2,
        fontweight="bold",
        color=COLORS["ink"],
    )
    ax.text(
        0.0,
        1.045,
        title,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=10.2,
        fontweight="bold",
        color=COLORS["ink"],
    )
    if subtitle:
        ax.text(
            1.0,
            1.045,
            subtitle,
            transform=ax.transAxes,
            ha="right",
            va="bottom",
            fontsize=6.4,
            color=COLORS["muted"],
        )


def rounded_box(ax: mpl.axes.Axes, xy, width, height, fc, ec=None, radius=0.035, lw=1.0):
    box = FancyBboxPatch(
        xy,
        width,
        height,
        boxstyle=f"round,pad=0.012,rounding_size={radius}",
        transform=ax.transAxes,
        facecolor=fc,
        edgecolor=ec if ec else fc,
        linewidth=lw,
    )
    ax.add_patch(box)
    return box


def map_broad(row: pd.Series) -> str:
    state = str(row.get("phase8_state", ""))
    broad = str(row.get("broad", ""))
    if state in {"Malignant", "Epithelial-Normal"} or broad == "Epithelial":
        return "Tumor / epithelial"
    if "Macrophage" in state or state in {"Monocyte", "Dendritic cell"} or broad == "Myeloid":
        return "Myeloid"
    if state == "T cell" or broad == "T_cell":
        return "T cell"
    if state == "NK cell" or broad == "NK_cell":
        return "NK cell"
    if state == "Endothelial" or broad == "Endothelial":
        return "Endothelial"
    if state in {"B cell", "Plasma cell"} or broad in {"B_cell", "Plasma_cell"}:
        return "B / plasma"
    if state == "SMC/Pericyte" or broad in {"Fibroblast", "Pericyte_SMC"}:
        return "Stromal"
    if state == "Mast cell" or broad == "Mast_cell":
        return "Mast"
    return "Other"


def robust_limits(values: pd.Series, low=0.02, high=0.995) -> tuple[float, float]:
    finite = values[np.isfinite(values)]
    if len(finite) == 0:
        return 0.0, 1.0
    lo, hi = np.quantile(finite, [low, high])
    if not np.isfinite(hi) or hi <= lo:
        hi = lo + 1.0
    return float(lo), float(hi)


def draw_panel_a(ax: mpl.axes.Axes) -> None:
    clean_axis(ax)
    panel_heading(
        ax,
        "A",
        "Immune-data expansion and multimodal analysis architecture",
        "Cohort-separated analyses • frozen Phase 8 labels",
    )
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    modules = [
        (0.015, 0.11, 0.20, 0.76, "Discovery atlas", "GSE207493\n19 patients • 160,710 cells", "#FFF3DD"),
        (0.265, 0.11, 0.22, 0.76, "External scRNA", "GSE264586 + GSE294109\n27 patients • 228,767 cells", "#E5F4F1"),
        (0.535, 0.11, 0.20, 0.76, "Spatial niche", "GSE175540\n12 FFPE + 12 frozen sections", "#F3EAF5"),
        (0.785, 0.11, 0.20, 0.76, "Clinical context", "CheckMate + JAVELIN\n311 + 726 RNA profiles", "#E8EEF7"),
    ]
    for x, y, w, h, title, detail, color in modules:
        rounded_box(ax, (x, y), w, h, color, "#D8DDE2", radius=0.03, lw=0.8)
        ax.text(x + w / 2, y + h - 0.11, title, transform=ax.transAxes, ha="center", va="center", fontsize=8.3, fontweight="bold")
        ax.text(x + w / 2, y + 0.10, detail, transform=ax.transAxes, ha="center", va="center", fontsize=6.3, color=COLORS["muted"], linespacing=1.25)

    # Discovery atlas glyph
    for i, color in enumerate([COLORS["tumor"], COLORS["tcell"], COLORS["monocyte"], COLORS["endothelial"], COLORS["c1qc"]]):
        theta = np.linspace(0, 2 * np.pi, 12, endpoint=False)
        rr = 0.018 + 0.004 * (i % 2)
        cx = 0.055 + 0.030 * (i % 3)
        cy = 0.48 + 0.09 * (i // 3)
        ax.scatter(cx + rr * np.cos(theta), cy + rr * np.sin(theta), s=10, color=color, alpha=0.85, transform=ax.transAxes, zorder=5)
    ax.add_patch(patches.Ellipse((0.118, 0.52), 0.13, 0.26, transform=ax.transAxes, fill=False, edgecolor="#5F6B75", lw=0.9, ls=(0, (2, 2))))

    # Three external-cohort embeddings
    for j, cx in enumerate([0.305, 0.375, 0.445]):
        n = 40
        ang = np.linspace(0, 2 * np.pi, n, endpoint=False)
        rad = 0.025 + 0.010 * np.sin(ang * (2 + j))
        xx = cx + rad * np.cos(ang)
        yy = 0.52 + 1.35 * rad * np.sin(ang)
        c = np.array([COLORS["monocyte"], COLORS["spp1"], COLORS["c1qc"], COLORS["other_mac"]])
        ax.scatter(xx, yy, s=8, c=c[np.arange(n) % 4], alpha=0.82, transform=ax.transAxes)
        ax.text(cx, 0.36, ["Discovery", "Validation 1", "Validation 2"][j], transform=ax.transAxes, ha="center", fontsize=5.4, color=COLORS["muted"])

    # Spatial tissue and honeycomb
    for j, (cx, col) in enumerate([(0.585, COLORS["tymp"]), (0.682, COLORS["c1qc"])]):
        t = np.linspace(0, 2 * np.pi, 100)
        rr = 0.060 * (1 + 0.18 * np.sin(3 * t + j))
        ax.fill(cx + rr * np.cos(t), 0.52 + 0.80 * rr * np.sin(t), color=col, alpha=0.22, transform=ax.transAxes)
        for ix in range(-2, 3):
            for iy in range(-2, 3):
                xx = cx + ix * 0.019 + (iy % 2) * 0.009
                yy = 0.52 + iy * 0.026
                if (xx - cx) ** 2 / 0.0045 + (yy - 0.52) ** 2 / 0.006 < 1:
                    ax.add_patch(patches.RegularPolygon((xx, yy), 6, radius=0.008, transform=ax.transAxes, facecolor=col, alpha=0.72, edgecolor="white", lw=0.25))
    ax.text(0.635, 0.36, "fine-state abundance\n+ tumour TYMP", transform=ax.transAxes, ha="center", fontsize=5.4, color=COLORS["muted"])

    # Clinical-profile glyph
    for row in range(2):
        for col in range(5):
            x = 0.817 + col * 0.028
            y = 0.56 - row * 0.11
            color = COLORS["tcell"] if row == 0 else COLORS["tymp"]
            ax.add_patch(patches.Circle((x, y + 0.025), 0.009, transform=ax.transAxes, facecolor=color, edgecolor="none", alpha=0.85))
            ax.add_patch(patches.FancyBboxPatch((x - 0.008, y - 0.025), 0.016, 0.043, boxstyle="round,pad=0.002", transform=ax.transAxes, facecolor=color, edgecolor="none", alpha=0.65))
    ax.text(0.885, 0.36, "response • PFS • OS\ncomposition-aware models", transform=ax.transAxes, ha="center", fontsize=5.4, color=COLORS["muted"])

    for x0, x1 in [(0.215, 0.265), (0.485, 0.535), (0.735, 0.785)]:
        arrow = FancyArrowPatch(
            (x0 + 0.008, 0.52),
            (x1 - 0.008, 0.52),
            transform=ax.transAxes,
            arrowstyle="-|>",
            mutation_scale=11,
            linewidth=1.8,
            color="#7D8790",
            connectionstyle="arc3,rad=0.0",
        )
        ax.add_patch(arrow)


def draw_umap(ax: mpl.axes.Axes, frame: pd.DataFrame, dataset: str) -> None:
    clean_axis(ax)
    order = list(BROAD_COLORS)
    for category in order:
        sub = frame[frame["plot_category"].eq(category)]
        if len(sub):
            ax.scatter(
                sub["UMAP_1"], sub["UMAP_2"], s=0.75, c=BROAD_COLORS[category], alpha=0.62,
                linewidths=0, rasterized=True
            )
    ax.set_aspect("equal", adjustable="datalim")
    status = frame.attrs.get("status", {})
    ax.set_title(
        f"{dataset}\n{status.get('patients', '?')} patients • {status.get('total_cells', len(frame)):,} cells",
        fontsize=7.2, fontweight="bold", pad=2, color=COLORS["ink"]
    )
    for category in ["Tumor / epithelial", "Myeloid", "T cell", "Endothelial"]:
        sub = frame[frame["plot_category"].eq(category)]
        if len(sub) > 150:
            x, y = sub[["UMAP_1", "UMAP_2"]].median().values
            label = {"Tumor / epithelial": "Tumour", "Myeloid": "Myeloid", "T cell": "T", "Endothelial": "Endo"}[category]
            ax.text(x, y, label, fontsize=5.2, ha="center", va="center", color=COLORS["ink"], fontweight="bold", path_effects=[pe.withStroke(linewidth=1.6, foreground="white")])


def draw_panel_b(fig: mpl.figure.Figure, spec, atlas_frames: dict[str, pd.DataFrame]) -> None:
    outer = fig.add_subplot(spec)
    clean_axis(outer)
    panel_heading(outer, "B", "Three-cohort immune atlas", "Existing UMAP coordinates • no pooled integration")
    gs = spec.subgridspec(1, 3, wspace=0.04)
    for i, dataset in enumerate(["GSE207493", "GSE264586", "GSE294109"]):
        ax = fig.add_subplot(gs[0, i])
        draw_umap(ax, atlas_frames[dataset], dataset)
    handles = [Line2D([0], [0], marker="o", color="none", markerfacecolor=color, markeredgecolor="none", markersize=4.3, label=label) for label, color in BROAD_COLORS.items()]
    outer.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.08), ncol=5, frameon=False, fontsize=5.4, handletextpad=0.3, columnspacing=0.8)


def feature_umap(ax, frame, column, cmap, title):
    clean_axis(ax)
    x, y = frame["UMAP_1"].to_numpy(), frame["UMAP_2"].to_numpy()
    values = frame[column].to_numpy(float)
    ax.scatter(x, y, s=0.55, c="#DADDE0", alpha=0.30, linewidths=0, rasterized=True)
    lo, hi = robust_limits(pd.Series(values), 0.05, 0.99)
    normed = np.clip((values - lo) / (hi - lo), 0, 1)
    mask = np.isfinite(normed) & (normed > 0.03)
    idx = np.argsort(normed[mask])
    ax.scatter(x[mask][idx], y[mask][idx], s=0.8, c=normed[mask][idx], cmap=cmap, vmin=0, vmax=1, alpha=0.80, linewidths=0, rasterized=True)
    ax.set_aspect("equal", adjustable="datalim")
    ax.set_title(title, fontsize=7.0, fontweight="bold", color=COLORS["ink"], pad=1)
    return normed


def draw_ring_glyph(ax, center, radius, fraction, color, title, lines):
    ax.add_patch(Wedge(center, radius, 0, 360, width=radius * 0.18, facecolor="#E9EDF0", edgecolor="none"))
    ax.add_patch(Wedge(center, radius, 90, 90 + 360 * fraction, width=radius * 0.18, facecolor=color, edgecolor="none"))
    ax.add_patch(patches.Circle(center, radius * 0.56, facecolor=color, alpha=0.10, edgecolor=color, lw=0.8))
    ax.text(center[0], center[1] + 0.005, title, ha="center", va="center", fontsize=6.0, fontweight="bold", color=COLORS["ink"])
    ax.text(center[0], center[1] - radius * 1.20, lines, ha="center", va="top", fontsize=5.2, color=COLORS["muted"], linespacing=1.25)


def draw_panel_c(fig, spec, atlas_frames, feature_summary, metabolic_summary, tf_summary):
    outer = fig.add_subplot(spec)
    clean_axis(outer)
    panel_heading(outer, "C", "Cellular source of TYMP and lineage-specific programs", "GSE207493 discovery cohort")
    gs = spec.subgridspec(2, 3, height_ratios=[0.69, 0.31], width_ratios=[1.0, 1.0, 1.28], hspace=0.16, wspace=0.10)
    frame = atlas_frames["GSE207493"]
    ax_tymp = fig.add_subplot(gs[0, 0])
    feature_umap(ax_tymp, frame, "TYMP_log1p", TYMP_CMAP, "TYMP expression")
    ax_f1 = fig.add_subplot(gs[0, 1])
    feature_umap(ax_f1, frame, "F1_score", F1_CMAP, "F1 projection score")

    ax_bubble = fig.add_subplot(gs[0, 2])
    show_types = ["Tumor", "Myeloid", "Dendritic", "Treg", "Exhausted_T", "CD8_T", "CD4_T", "Endothelial", "CAF", "NK", "B_cell", "Plasma"]
    fs = feature_summary.set_index("cell_type").reindex(show_types).dropna(how="all")
    metrics = ["f1_mean", "TYMP_mean", "TYMP_positive_fraction"]
    labels = ["F1 mean", "TYMP mean", "TYMP+ cells"]
    raw = fs[metrics].copy()
    scaled = (raw - raw.min()) / (raw.max() - raw.min()).replace(0, 1)
    n = fs["n_cells"].fillna(0).to_numpy()
    sizes = 12 + 100 * np.sqrt(n / max(n.max(), 1))
    for yi, cell_type in enumerate(fs.index):
        for xi, metric in enumerate(metrics):
            ax_bubble.scatter(xi, yi, s=sizes[yi], c=[TYMP_CMAP(float(scaled.loc[cell_type, metric]))], edgecolors="white", linewidths=0.45)
    ax_bubble.set_xlim(-0.55, len(metrics) - 0.45)
    ax_bubble.set_ylim(len(fs) - 0.5, -0.5)
    ax_bubble.set_xticks(range(len(metrics)), labels, rotation=35, ha="right", fontsize=5.2)
    ax_bubble.set_yticks(range(len(fs)), [x.replace("_", " ") for x in fs.index], fontsize=5.2)
    ax_bubble.grid(color=COLORS["grid"], lw=0.4, alpha=0.65)
    ax_bubble.tick_params(length=0)
    for spine in ax_bubble.spines.values():
        spine.set_visible(False)
    ax_bubble.set_title("Lineage feature landscape", fontsize=7.0, fontweight="bold", pad=1)

    ax_ring = fig.add_subplot(gs[1, :])
    clean_axis(ax_ring)
    ax_ring.set_xlim(0, 1)
    ax_ring.set_ylim(0, 1)
    contrast_to_data = {
        "F1-high tumour": ("F1_high_vs_low_Tumor", "tumor_f1", COLORS["tumor"]),
        "TYMP-high Myeloid/DC": ("TYMP_high_vs_low_MyeloidDC", "myeloid_dc_tymp", COLORS["tymp"]),
        "TYMP-high Treg": ("TYMP_high_vs_low_Treg", "treg_tymp", COLORS["c1qc"]),
    }
    for x, (title, (contrast, endpoint, color)) in zip([0.17, 0.50, 0.83], contrast_to_data.items()):
        m = metabolic_summary.loc[metabolic_summary["contrast"].eq(contrast)].iloc[0]
        t = tf_summary.loc[tf_summary["endpoint"].eq(endpoint)].iloc[0]
        draw_ring_glyph(
            ax_ring,
            (x, 0.59),
            0.16,
            float(m["adjusted_q_lt_0.05"]) / float(m["n_pathways"]),
            color,
            title,
            f"{int(m['adjusted_q_lt_0.05'])}/{int(m['n_pathways'])} pathways • {int(t['unique_cross_method_tfs'])} TFs",
        )
    ax_ring.text(0.5, 0.02, "Outer ring = cycle-adjusted metabolic FDR coverage; centre label = lineage and anchor", ha="center", va="bottom", fontsize=5.1, color=COLORS["muted"])


def shorten_label(label: str, limit=34) -> str:
    replacements = {
        " positive regulon top-weighted targets": " regulon",
        "Aerobic respiration and respiratory electron transport": "Aerobic respiration + ETC",
        "Regulation of pyruvate dehydrogenase (PDH) complex": "PDH-complex regulation",
        "Detoxification of Reactive Oxygen Species": "ROS detoxification",
        "Cellular response to hypoxia": "Hypoxia response",
        "Pentose phosphate pathway": "Pentose-phosphate pathway",
    }
    text = replacements.get(label, label)
    if len(text) > limit:
        text = text[: limit - 1].rstrip() + "…"
    return text


def select_phase67_rows(retained: pd.DataFrame) -> pd.DataFrame:
    myeloid = retained[(retained["endpoint"].eq("myeloid_dc_tymp")) & (retained["replication_grade"].eq("R2_TWO_EXTERNAL_FDR"))].copy()
    treg = retained[(retained["endpoint"].eq("treg_tymp")) & (retained["datasets_fdr_replication"].ge(1))].copy()
    selected = pd.concat([myeloid, treg], ignore_index=True)
    endpoint_order = {"myeloid_dc_tymp": 0, "treg_tymp": 1}
    family_order = {"metabolic": 0, "regulon": 1}
    selected["endpoint_order"] = selected["endpoint"].map(endpoint_order)
    selected["family_order"] = selected["family"].map(family_order)
    selected = selected.sort_values(["endpoint_order", "family_order", "label"]).reset_index(drop=True)
    return selected


def draw_program_constellation(ax, selected: pd.DataFrame) -> None:
    clean_axis(ax)
    ax.set_xlim(-1.28, 1.28)
    ax.set_ylim(-1.28, 1.28)
    core = selected[selected["endpoint"].eq("myeloid_dc_tymp")].copy()
    core["mean_rho"] = core[["GSE264586_median_adjusted_rho", "GSE294109_median_adjusted_rho"]].mean(axis=1)
    core = pd.concat([
        core[core["family"].eq("metabolic")].nlargest(7, "mean_rho"),
        core[core["family"].eq("regulon")].nlargest(7, "mean_rho"),
    ]).drop_duplicates("signature_id")
    n = len(core)
    angles = np.linspace(np.pi / 2, np.pi / 2 + 2 * np.pi, n, endpoint=False)
    ax.add_patch(patches.Circle((0, 0), 0.27, facecolor=COLORS["tymp"], edgecolor="white", lw=1.2, zorder=5))
    ax.add_patch(patches.Circle((0, 0), 0.37, facecolor=COLORS["tymp"], edgecolor="none", alpha=0.10, zorder=1))
    ax.text(0, 0.03, "TYMP", ha="center", va="center", color="white", fontsize=8.2, fontweight="bold", zorder=6)
    ax.text(0, -0.10, "Myeloid/DC", ha="center", va="center", color="white", fontsize=5.0, zorder=6)
    for angle, (_, row) in zip(angles, core.iterrows()):
        radius = 0.88 if row["family"] == "metabolic" else 1.04
        x, y = radius * np.cos(angle), radius * np.sin(angle)
        value = float(row["mean_rho"])
        edge_color = COLORS["tymp"] if value >= 0 else COLORS["monocyte"]
        ax.plot([0, x], [0, y], color=edge_color, alpha=0.28 + min(abs(value), 0.4), lw=0.7 + 4.5 * abs(value), zorder=2)
        node_color = "#E99A2E" if row["family"] == "metabolic" else "#6A4C93"
        ax.scatter([x], [y], s=62 + 240 * abs(value), c=node_color, edgecolors="white", linewidths=0.7, zorder=4)
        label = shorten_label(str(row["label"]), 20)
        ha = "left" if x >= 0 else "right"
        ax.text(x + (0.055 if x >= 0 else -0.055), y, label, ha=ha, va="center", fontsize=4.5, color=COLORS["ink"])
    ax.set_title("Replicated program constellation", fontsize=7.3, fontweight="bold", pad=2)
    ax.text(0, -1.24, "Edge width = mean |rho| across two external cohorts", ha="center", fontsize=4.9, color=COLORS["muted"])


def draw_panel_d(fig, spec, retained):
    outer = fig.add_subplot(spec)
    clean_axis(outer)
    panel_heading(outer, "D", "Cross-cohort replication of TYMP-associated immune programs", "Patient-level residual correlations • frozen signatures")
    gs = spec.subgridspec(1, 2, width_ratios=[1.62, 0.92], wspace=0.16)
    selected = select_phase67_rows(retained)
    ax = fig.add_subplot(gs[0, 0])
    matrix = selected[["GSE264586_median_adjusted_rho", "GSE294109_median_adjusted_rho"]].to_numpy(float)
    im = ax.imshow(matrix, aspect="auto", cmap=DIVERGING, vmin=-0.45, vmax=0.45, interpolation="nearest")
    labels = [shorten_label(str(x), 39) for x in selected["label"]]
    ax.set_yticks(np.arange(len(selected)), labels, fontsize=5.15)
    ax.set_xticks([0, 1], ["GSE264586\n12 patients", "GSE294109\n13 patients"], fontsize=6.0, fontweight="bold")
    ax.tick_params(axis="both", length=0)
    for i, row in selected.iterrows():
        for j, dataset in enumerate(["GSE264586", "GSE294109"]):
            rho = row[f"{dataset}_median_adjusted_rho"]
            q = row[f"{dataset}_adjusted_q_value"]
            text_color = "white" if np.isfinite(rho) and abs(rho) > 0.23 else COLORS["ink"]
            txt = "NA" if not np.isfinite(rho) else f"{rho:+.2f}" + ("*" if np.isfinite(q) and q < 0.05 else "")
            ax.text(j, i, txt, ha="center", va="center", fontsize=4.7, color=text_color, fontweight="bold" if "*" in txt else "normal")
    for spine in ax.spines.values():
        spine.set_visible(False)
    # Group strips and separators
    groups = selected.assign(group=selected["endpoint"].str.replace("_tymp", "", regex=False) + " | " + selected["family"])
    boundaries = []
    for group, sub in groups.groupby("group", sort=False):
        start, end = sub.index.min(), sub.index.max()
        color = COLORS["tymp"] if group.startswith("myeloid") else COLORS["c1qc"]
        if "metabolic" in group:
            color = COLORS["tumor"] if group.startswith("myeloid") else "#B47A9B"
        ax.add_patch(patches.Rectangle((-0.68, start - 0.5), 0.10, end - start + 1, transform=ax.transData, facecolor=color, edgecolor="none", clip_on=False))
        ax.text(-0.78, (start + end) / 2, group.replace("myeloid_dc", "Myeloid/DC").replace("treg", "Treg").replace(" | ", "\n"), ha="right", va="center", fontsize=4.8, color=COLORS["muted"], clip_on=False)
        boundaries.append(end + 0.5)
    for y in boundaries[:-1]:
        ax.axhline(y, color="white", lw=1.6)
    cax = ax.inset_axes([0.18, -0.09, 0.64, 0.025])
    cb = plt.colorbar(im, cax=cax, orientation="horizontal")
    cb.outline.set_visible(False)
    cb.ax.tick_params(labelsize=4.8, length=2)
    cb.set_label("Median adjusted Spearman rho", fontsize=5.2)
    ax.text(1.02, 1.005, "* q < 0.05", transform=ax.transAxes, ha="right", va="bottom", fontsize=5.0, color=COLORS["muted"])
    ax_net = fig.add_subplot(gs[0, 1])
    draw_program_constellation(ax_net, selected)
    ax_net.text(-1.18, 1.19, f"{sum((selected.endpoint == 'myeloid_dc_tymp'))} two-cohort Myeloid/DC signatures", fontsize=5.3, color=COLORS["tymp"], fontweight="bold")
    ax_net.text(-1.18, 1.08, f"{sum((selected.endpoint == 'treg_tymp'))} single-external Treg signatures", fontsize=5.3, color=COLORS["c1qc"], fontweight="bold")


def trajectory_path(frame: pd.DataFrame, bins=18) -> pd.DataFrame:
    values = frame.dropna(subset=["slingshot_pseudotime", "UMAP_1", "UMAP_2"]).copy()
    if len(values) < 100:
        return pd.DataFrame()
    values["pt_bin"] = pd.cut(values["slingshot_pseudotime"], bins=np.linspace(0, 1, bins + 1), include_lowest=True, labels=False)
    path = values.groupby("pt_bin", observed=True)[["UMAP_1", "UMAP_2", "slingshot_pseudotime"]].median().dropna().reset_index()
    return path


def program_dynamics(frame: pd.DataFrame, bins=20) -> pd.DataFrame:
    cols = ["MONOCYTE_PROGRAM", "SPP1_MARKER20_ANCHOR_EXCLUDED", "C1QC_MARKER20_ANCHOR_EXCLUDED", "IMMUNOSUPPRESSION_NO_IL10"]
    x = frame.dropna(subset=["slingshot_pseudotime"]).copy()
    x["pt_bin"] = pd.cut(x["slingshot_pseudotime"], bins=np.linspace(0, 1, bins + 1), include_lowest=True, labels=False)
    med = x.groupby("pt_bin", observed=True)[cols + ["slingshot_pseudotime"]].median().reset_index()
    for col in cols:
        lo, hi = med[col].min(), med[col].max()
        med[col] = (med[col] - lo) / (hi - lo) if hi > lo else 0.5
    return med


def draw_panel_e(fig, spec, atlas_frames, pseudotime, concordance):
    outer = fig.add_subplot(spec)
    clean_axis(outer)
    panel_heading(outer, "E", "Reproducible myeloid-state ordering", "Slingshot + PCA-kNN geodesic pseudotime")
    gs = spec.subgridspec(2, 3, height_ratios=[0.60, 0.40], hspace=0.12, wspace=0.08)
    datasets = ["GSE207493", "GSE264586", "GSE294109"]
    line_colors = {
        "MONOCYTE_PROGRAM": COLORS["monocyte"],
        "SPP1_MARKER20_ANCHOR_EXCLUDED": COLORS["spp1"],
        "C1QC_MARKER20_ANCHOR_EXCLUDED": COLORS["c1qc"],
        "IMMUNOSUPPRESSION_NO_IL10": COLORS["tymp"],
    }
    line_labels = {
        "MONOCYTE_PROGRAM": "Monocyte",
        "SPP1_MARKER20_ANCHOR_EXCLUDED": "SPP1-X",
        "C1QC_MARKER20_ANCHOR_EXCLUDED": "C1QC-X",
        "IMMUNOSUPPRESSION_NO_IL10": "Immunosuppression",
    }
    for j, dataset in enumerate(datasets):
        ax = fig.add_subplot(gs[0, j])
        clean_axis(ax)
        atlas = atlas_frames[dataset]
        myeloid = atlas[atlas["phase8_state"].isin(STATE_COLORS)].copy()
        joined = myeloid.merge(
            pseudotime[pseudotime["dataset"].eq(dataset)][["cell_id", "slingshot_pseudotime", "geodesic_pseudotime"]],
            on="cell_id", how="inner"
        )
        for state, color in STATE_COLORS.items():
            sub = joined[joined["phase8_state"].eq(state)]
            ax.scatter(sub["UMAP_1"], sub["UMAP_2"], s=0.75, color=color, alpha=0.62, linewidths=0, rasterized=True)
        path = trajectory_path(joined)
        if len(path) > 2:
            ax.plot(path["UMAP_1"], path["UMAP_2"], color="#1D2228", lw=1.3, alpha=0.85, zorder=6)
            ax.add_patch(FancyArrowPatch(
                (path.iloc[-3]["UMAP_1"], path.iloc[-3]["UMAP_2"]),
                (path.iloc[-1]["UMAP_1"], path.iloc[-1]["UMAP_2"]),
                arrowstyle="-|>", mutation_scale=8, color="#1D2228", lw=1.2, zorder=7
            ))
        rho = concordance.loc[concordance["dataset"].eq(dataset), "median_patient_rho"].iloc[0]
        ax.set_title(f"{dataset}  method rho={rho:.2f}", fontsize=6.5, fontweight="bold", pad=1)
        ax.set_aspect("equal", adjustable="datalim")

        axd = fig.add_subplot(gs[1, j])
        dyn = program_dynamics(pseudotime[pseudotime["dataset"].eq(dataset)])
        for col, color in line_colors.items():
            axd.plot(dyn["slingshot_pseudotime"], dyn[col], color=color, lw=1.3, label=line_labels[col])
        axd.set_xlim(0, 1)
        axd.set_ylim(-0.03, 1.03)
        axd.set_xticks([0, 0.5, 1.0])
        axd.set_yticks([0, 1], ["low", "high"])
        axd.tick_params(labelsize=4.7, length=2, color=COLORS["grid"])
        axd.grid(color=COLORS["grid"], lw=0.4, alpha=0.6)
        for spine in ["top", "right"]:
            axd.spines[spine].set_visible(False)
        axd.spines["left"].set_color(COLORS["grid"])
        axd.spines["bottom"].set_color(COLORS["grid"])
        axd.set_xlabel("Slingshot pseudotime", fontsize=5.0)
        if j == 0:
            axd.set_ylabel("scaled median", fontsize=5.0)
            axd.legend(loc="upper center", bbox_to_anchor=(0.5, -0.33), ncol=2, frameon=False, fontsize=4.7, handlelength=1.8)
    handles = [Line2D([0], [0], marker="o", lw=0, markerfacecolor=c, markeredgecolor="none", label=s, markersize=4) for s, c in STATE_COLORS.items()]
    outer.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.10), ncol=2, frameon=False, fontsize=4.9, handletextpad=0.3, columnspacing=0.8)


def read_spatial_sample(root: Path, sample: str) -> pd.DataFrame:
    meta = pd.read_csv(root / sample / "spot_metadata_v1.tsv.gz", sep="\t")
    abundance = pd.read_csv(root / sample / "fine_state_abundance_v1.tsv.gz", sep="\t")
    return meta.merge(abundance, on="spot_id", how="inner", validate="one_to_one")


def draw_spatial_map(ax, frame, feature, title, cmap=SPATIAL_CMAP):
    clean_axis(ax)
    values = frame[feature].astype(float)
    lo, hi = robust_limits(values, 0.03, 0.98)
    idx = np.argsort(values.fillna(lo).to_numpy())
    ax.scatter(frame["pxl_col"].to_numpy()[idx], -frame["pxl_row"].to_numpy()[idx], c=values.to_numpy()[idx], cmap=cmap, vmin=lo, vmax=hi, s=2.5, linewidths=0, rasterized=True)
    ax.set_aspect("equal", adjustable="datalim")
    ax.set_title(title, fontsize=5.6, pad=1, fontweight="bold")


def draw_cellchat_network(ax, interactions: pd.DataFrame) -> None:
    clean_axis(ax)
    selected = ["Malignant", "SPP1-state Macrophage", "C1QC-state Macrophage", "Other Macrophage", "T cell", "Endothelial"]
    x = interactions[
        interactions["dataset"].eq("GSE207493")
        & interactions["source"].isin(selected)
        & interactions["target"].isin(selected)
        & interactions["pval"].le(0.05)
    ].copy()
    unit_edges = x.groupby(["unit_id", "source", "target"], observed=True)["prob"].sum().reset_index()
    edge = unit_edges.groupby(["source", "target"], observed=True)["prob"].median().reset_index()
    edge = edge.sort_values("prob", ascending=False).head(20)
    angle = np.linspace(np.pi / 2, np.pi / 2 + 2 * np.pi, len(selected), endpoint=False)
    positions = {node: np.array([np.cos(a), np.sin(a)]) for node, a in zip(selected, angle)}
    node_colors = {
        "Malignant": COLORS["tumor"],
        "SPP1-state Macrophage": COLORS["spp1"],
        "C1QC-state Macrophage": COLORS["c1qc"],
        "Other Macrophage": COLORS["other_mac"],
        "T cell": COLORS["tcell"],
        "Endothelial": COLORS["endothelial"],
    }
    max_prob = max(edge["prob"].max(), 1e-9)
    for _, row in edge.iterrows():
        src, tgt = row["source"], row["target"]
        p0, p1 = positions[src], positions[tgt]
        rad = 0.18 if src != tgt else 0.35
        arrow = FancyArrowPatch(
            p0 * 0.87, p1 * 0.87, connectionstyle=f"arc3,rad={rad}", arrowstyle="-|>",
            mutation_scale=5.5, lw=0.35 + 2.2 * float(row["prob"]) / max_prob,
            color=node_colors[src], alpha=0.25 + 0.45 * float(row["prob"]) / max_prob, zorder=1
        )
        ax.add_patch(arrow)
    for node, pos in positions.items():
        ax.scatter(*pos, s=150, color=node_colors[node], edgecolors="white", linewidths=1.0, zorder=5)
        label = node.replace("-state Macrophage", "\nmacrophage").replace("Other Macrophage", "Other\nmacrophage")
        ha = "left" if pos[0] > 0.2 else "right" if pos[0] < -0.2 else "center"
        va = "bottom" if pos[1] > 0.2 else "top" if pos[1] < -0.2 else "center"
        ax.text(pos[0] + 0.06 * np.sign(pos[0]), pos[1] + 0.06 * np.sign(pos[1]), label, ha=ha, va=va, fontsize=4.6, color=COLORS["ink"])
    ax.set_xlim(-1.35, 1.35)
    ax.set_ylim(-1.35, 1.35)
    ax.set_title("Absolute CellChat network\nGSE207493 patient median", fontsize=6.0, fontweight="bold", pad=1)


def draw_endpoint_bubbles(ax, endpoint_summary: pd.DataFrame) -> None:
    primary = endpoint_summary[endpoint_summary["tier"].str.contains("PRIMARY", na=False)].copy()
    datasets = ["GSE207493", "GSE264586", "GSE294109"]
    endpoint_order = list(dict.fromkeys(primary["endpoint_id"].tolist()))
    endpoint_order = sorted(endpoint_order)[:6]
    short = {
        "P01_Malignant_to_SPP1_specificity": "Malignant → SPP1",
        "P02_SPP1_to_Malignant_specificity": "SPP1 → Malignant",
        "P03_SPP1_to_Tcell_specificity": "SPP1 → T cell",
        "P04_SPP1_to_Endothelial_specificity": "SPP1 → Endothelial",
        "P05_SPP1_total_outgoing_specificity": "SPP1 outgoing",
        "P06_SPP1_total_incoming_specificity": "SPP1 incoming",
    }
    for yi, endpoint in enumerate(endpoint_order):
        for xi, dataset in enumerate(datasets):
            row = primary[(primary["endpoint_id"].eq(endpoint)) & (primary["dataset"].eq(dataset))]
            if row.empty:
                continue
            row = row.iloc[0]
            rho = float(row["spearman_rho"]) if np.isfinite(row["spearman_rho"]) else 0
            n = float(row["n_patients"])
            ax.scatter(xi, yi, s=16 + 5.5 * n, c=[DIVERGING(np.clip((rho + 0.7) / 1.4, 0, 1))], edgecolors=COLORS["ink"] if bool(row["formal_support"]) else "white", linewidths=0.8)
            if n < 10:
                ax.text(xi, yi, f"n={int(n)}", ha="center", va="center", fontsize=3.7, color=COLORS["ink"])
    ax.set_xlim(-0.5, len(datasets) - 0.5)
    ax.set_ylim(len(endpoint_order) - 0.5, -0.5)
    ax.set_xticks(range(len(datasets)), ["GSE207493", "GSE264586", "GSE294109"], rotation=30, ha="right", fontsize=4.7)
    ax.set_yticks(range(len(endpoint_order)), [short.get(x, x.replace("_specificity", "").replace("P0", "")) for x in endpoint_order], fontsize=4.7)
    ax.grid(color=COLORS["grid"], lw=0.4, alpha=0.7)
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_title("SPP1-versus-C1QC/Other specificity", fontsize=6.0, fontweight="bold", pad=1)
    ax.text(1.0, -0.20, "bubble size = patients • outline = formal support", transform=ax.transAxes, ha="right", fontsize=4.4, color=COLORS["muted"])


def draw_panel_f(fig, spec, spatial_root, interactions, endpoint_summary):
    outer = fig.add_subplot(spec)
    clean_axis(outer)
    panel_heading(outer, "F", "Spatial immune ecology and fine-state communication", "Representative sections selected near preparation-level effect")
    gs = spec.subgridspec(2, 1, height_ratios=[0.58, 0.42], hspace=0.16)
    maps = gs[0].subgridspec(2, 4, wspace=0.03, hspace=0.08)
    samples = [("GSM5924037", "FFPE"), ("GSM5924047", "Frozen")]
    features = [
        ("score_primary__tymp_expression", "Tumour TYMP", TYMP_CMAP),
        ("q05__SPP1-state Macrophage", "SPP1 macrophage", SPATIAL_CMAP),
        ("q05__C1QC-state Macrophage", "C1QC macrophage", SPATIAL_CMAP),
        ("q05__Endothelial", "Endothelial", SPATIAL_CMAP),
    ]
    for i, (sample, prep) in enumerate(samples):
        frame = read_spatial_sample(spatial_root, sample)
        for j, (feature, title, cmap) in enumerate(features):
            ax = fig.add_subplot(maps[i, j])
            draw_spatial_map(ax, frame, feature, title if i == 0 else "", cmap)
            if j == 0:
                ax.text(-0.05, 0.5, f"{prep}\n{sample}", transform=ax.transAxes, ha="right", va="center", rotation=90, fontsize=5.2, fontweight="bold", color=COLORS["muted"])
    lower = gs[1].subgridspec(1, 2, width_ratios=[0.95, 1.15], wspace=0.22)
    ax_net = fig.add_subplot(lower[0, 0])
    draw_cellchat_network(ax_net, interactions)
    ax_bub = fig.add_subplot(lower[0, 1])
    draw_endpoint_bubbles(ax_bub, endpoint_summary)


def main() -> None:
    if len(sys.argv) != 4:
        raise SystemExit("Usage: 01_build_main_figure_v1.py <project_root> <work_root> <output_dir>")
    project_root = Path(sys.argv[1])
    work_root = Path(sys.argv[2])
    output_dir = Path(sys.argv[3])
    output_dir.mkdir(parents=True, exist_ok=False)

    mpl.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 6.2,
        "axes.titlesize": 7.0,
        "axes.labelsize": 5.5,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "axes.linewidth": 0.6,
        "savefig.facecolor": "white",
    })

    status = pd.read_csv(work_root / "atlas_extract_v1" / "atlas_extraction_status_v1.tsv", sep="\t")
    atlas_frames: dict[str, pd.DataFrame] = {}
    for _, row in status.iterrows():
        dataset = row["dataset"]
        frame = pd.read_csv(work_root / "atlas_extract_v1" / f"{dataset}_atlas_umap_sample_v1.tsv.gz", sep="\t")
        frame["plot_category"] = frame.apply(map_broad, axis=1)
        frame.attrs["status"] = row.to_dict()
        atlas_frames[dataset] = frame

    phase65 = project_root / "phase65_part1_immune_mechanistic_convergence_local_v1"
    phase67 = project_root / "phase67_part1_direct_program_replication_local_v1"
    phase87 = project_root / "phase87_three_cohort_myeloid_trajectory_local_v1"
    phase88 = project_root / "phase88_three_cohort_fine_state_cellchat_local_v1"
    feature_summary_path = phase65 / "inputs_v1" / "gse207493_celltype_feature_summary_v1.tsv"
    metabolic_summary_path = phase65 / "results_v2" / "METABOLIC_ENDPOINT_SUMMARY_v1.tsv"
    tf_summary_path = phase65 / "results_v2" / "TF_CROSS_METHOD_ENDPOINT_SUMMARY_v1.tsv"
    retained_path = phase67 / "results_v5" / "phase67_cross_dataset_retained_v1.tsv"
    pseudotime_path = phase87 / "results_v1" / "phase87_cell_pseudotime_v1.tsv.gz"
    concordance_path = phase87 / "results_v1" / "phase87_dataset_method_concordance_summary_v1.tsv"
    interactions_path = phase88 / "results_v1" / "phase88_cellchat_interactions_all_v1.tsv.gz"
    endpoint_summary_path = phase88 / "results_v1" / "phase88_dataset_endpoint_summary_v1.tsv"
    prep_meta_path = project_root / "phase69_85_server_status_sync_local_v1" / "phase78" / "analysis_v1" / "phase78_preparation_meta_v1.tsv"
    section_effects_path = project_root / "phase69_85_server_status_sync_local_v1" / "phase78" / "analysis_v1" / "phase78_section_effects_v1.tsv"

    feature_summary = pd.read_csv(feature_summary_path, sep="\t")
    metabolic_summary = pd.read_csv(metabolic_summary_path, sep="\t")
    tf_summary = pd.read_csv(tf_summary_path, sep="\t")
    retained = pd.read_csv(retained_path, sep="\t")
    pseudotime = pd.read_csv(pseudotime_path, sep="\t")
    concordance = pd.read_csv(concordance_path, sep="\t")
    interactions = pd.read_csv(interactions_path, sep="\t")
    endpoint_summary = pd.read_csv(endpoint_summary_path, sep="\t")

    fig = plt.figure(figsize=(14.0, 18.4), facecolor="white")
    outer = fig.add_gridspec(
        4, 1, height_ratios=[1.78, 4.15, 5.55, 6.00],
        left=0.045, right=0.975, top=0.978, bottom=0.030, hspace=0.21
    )
    ax_a = fig.add_subplot(outer[0])
    draw_panel_a(ax_a)
    row_bc = outer[1].subgridspec(1, 2, width_ratios=[1.08, 0.92], wspace=0.14)
    draw_panel_b(fig, row_bc[0], atlas_frames)
    draw_panel_c(fig, row_bc[1], atlas_frames, feature_summary, metabolic_summary, tf_summary)
    draw_panel_d(fig, outer[2], retained)
    row_ef = outer[3].subgridspec(1, 2, width_ratios=[1.0, 1.0], wspace=0.15)
    draw_panel_e(fig, row_ef[0], atlas_frames, pseudotime, concordance)
    draw_panel_f(fig, row_ef[1], work_root / "spatial_inputs_v1", interactions, endpoint_summary)

    fig.text(
        0.5, 0.010,
        "External immune evidence after cohort expansion | Patient-level inference; cohort-separated single-cell analyses; Phase 8 labels frozen",
        ha="center", va="bottom", fontsize=6.0, color=COLORS["muted"]
    )

    outputs = {
        "png": output_dir / "PROJECT26210_IMMUNE_MAIN_FIGURE_v1.png",
        "pdf": output_dir / "PROJECT26210_IMMUNE_MAIN_FIGURE_v1.pdf",
        "svg": output_dir / "PROJECT26210_IMMUNE_MAIN_FIGURE_v1.svg",
    }
    fig.savefig(outputs["png"], dpi=300, bbox_inches="tight", pad_inches=0.08)
    fig.savefig(outputs["pdf"], dpi=300, bbox_inches="tight", pad_inches=0.08)
    fig.savefig(outputs["svg"], dpi=300, bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)

    source_paths = [
        feature_summary_path, metabolic_summary_path, tf_summary_path, retained_path,
        pseudotime_path, concordance_path, interactions_path, endpoint_summary_path,
        prep_meta_path, section_effects_path,
        work_root / "atlas_extract_v1" / "atlas_extraction_status_v1.tsv",
    ]
    for sample in ["GSM5924037", "GSM5924047"]:
        source_paths.extend([
            work_root / "spatial_inputs_v1" / sample / "spot_metadata_v1.tsv.gz",
            work_root / "spatial_inputs_v1" / sample / "fine_state_abundance_v1.tsv.gz",
        ])
    manifest_rows = []
    for path in source_paths:
        manifest_rows.append({"path": str(path), "bytes": path.stat().st_size, "sha256": sha256(path)})
    pd.DataFrame(manifest_rows).to_csv(output_dir / "PANEL_SOURCE_MANIFEST_SHA256_v1.tsv", sep="\t", index=False)

    summary = {
        "status": "GENERATED_PENDING_VISUAL_QA",
        "seed": SEED,
        "figure_inches": [14.0, 18.4],
        "atlas_total_cells": int(status["total_cells"].sum()),
        "atlas_sampled_cells": int(status["sampled_cells"].sum()),
        "phase67_heatmap_rows": int(len(select_phase67_rows(retained))),
        "phase87_cells": int(len(pseudotime)),
        "phase88_interaction_rows": int(len(interactions)),
        "spatial_samples": ["GSM5924037", "GSM5924047"],
        "outputs": {key: {"path": str(path), "bytes": path.stat().st_size, "sha256": sha256(path)} for key, path in outputs.items()},
    }
    (output_dir / "BUILD_SUMMARY_v1.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
