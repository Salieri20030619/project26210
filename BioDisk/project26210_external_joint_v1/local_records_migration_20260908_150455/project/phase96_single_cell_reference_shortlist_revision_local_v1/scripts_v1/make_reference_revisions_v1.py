from __future__ import annotations

import importlib.util
import json
import math
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Ellipse, FancyArrowPatch, Polygon
from scipy.ndimage import gaussian_filter
from shapely.geometry import MultiPoint


BASE_SCRIPT = Path(r"C:\Users\31937\project26210_single_cell_reference_candidates_v1\make_single_cell_candidates_v1.py")
spec = importlib.util.spec_from_file_location("candidate_base", BASE_SCRIPT)
base = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(base)


def density_envelope(ax, frame: pd.DataFrame, color: str, extent: tuple[float, float, float, float]) -> None:
    """Draw disconnected local density envelopes that cover essentially all plotted cells."""
    pts = frame[["UMAP_1", "UMAP_2"]].dropna().to_numpy(float)
    if len(pts) < 35:
        return
    x0, x1, y0, y1 = extent
    hist, xedges, yedges = np.histogram2d(
        pts[:, 0], pts[:, 1], bins=(210, 210), range=((x0, x1), (y0, y1))
    )
    smooth = gaussian_filter(hist.T, sigma=1.25)
    xi = np.clip(np.searchsorted(xedges, pts[:, 0], side="right") - 1, 0, smooth.shape[1] - 1)
    yi = np.clip(np.searchsorted(yedges, pts[:, 1], side="right") - 1, 0, smooth.shape[0] - 1)
    point_density = smooth[yi, xi]
    point_density = point_density[point_density > 0]
    if not len(point_density):
        return
    level = max(float(np.quantile(point_density, 0.001)), 0.015)
    xc = (xedges[:-1] + xedges[1:]) / 2
    yc = (yedges[:-1] + yedges[1:]) / 2
    ax.contour(
        xc, yc, smooth, levels=[level], colors=[color], linewidths=0.72,
        linestyles=[(0, (2.4, 2.0))], alpha=0.82, zorder=8,
    )


def plot_f4a_fixed(atlas_frames, status, output_dir: Path) -> dict[str, str]:
    fig, axes = plt.subplots(1, 3, figsize=(14.5, 5.9))
    for ax, dataset in zip(axes, base.DATASETS):
        frame = atlas_frames[dataset]
        qx = frame["UMAP_1"].quantile([0.001, 0.999]).to_numpy()
        qy = frame["UMAP_2"].quantile([0.001, 0.999]).to_numpy()
        dx, dy = qx[1] - qx[0], qy[1] - qy[0]
        extent = (
            qx[0] - 0.03 * dx, qx[1] + 0.03 * dx,
            qy[0] - 0.03 * dy, qy[1] + 0.03 * dy,
        )
        for state in base.STATE_ORDER:
            sub = frame[frame["phase8_state"].eq(state)]
            if len(sub):
                ax.scatter(
                    sub["UMAP_1"], sub["UMAP_2"], s=1.2, c=base.STATE_COLORS[state],
                    alpha=0.58, linewidths=0, rasterized=True,
                )
                density_envelope(ax, sub, base.STATE_COLORS[state], extent)
        row = status.loc[status["dataset"].eq(dataset)].iloc[0]
        ax.set_title(
            f"{dataset}\n{int(row['patients'])} patients · {int(row['total_cells']):,} cells",
            fontweight="bold", color=base.INK, pad=7,
        )
        ax.set_xlim(extent[0], extent[1])
        ax.set_ylim(extent[2], extent[3])
        ax.set_aspect("equal", adjustable="box")
        base.clean_axis(ax)
    fig.suptitle(
        "Cohort-separated single-cell atlas with local state envelopes",
        fontsize=16, fontweight="bold", color=base.INK, y=0.985,
    )
    fig.text(
        0.5, 0.915,
        "Dashed contours are density-derived local envelopes; disconnected components are outlined separately and covariance ellipses are not used",
        ha="center", color=base.MUTED, fontsize=9.1,
    )
    handles = [
        Line2D([0], [0], marker="o", lw=0, markersize=6, markerfacecolor=base.STATE_COLORS[s],
               markeredgecolor="none", label=s.replace("-state Macrophage", " Mφ"))
        for s in base.STATE_ORDER
    ]
    fig.legend(
        handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.012), ncol=7,
        frameon=False, fontsize=8, handletextpad=0.35, columnspacing=0.9,
    )
    fig.subplots_adjust(left=0.025, right=0.99, top=0.80, bottom=0.16, wspace=0.06)
    return base.save_all(fig, output_dir, "01_f4a_local_density_envelopes_v2")


GROUPS = {
    "Tumour": ["Malignant"],
    "Myeloid": ["Monocyte", "SPP1-state Macrophage", "C1QC-state Macrophage", "Other Macrophage"],
    "Lymphoid": ["T cell", "NK cell"],
    "Vascular": ["Endothelial"],
}
GROUP_COLORS = {
    "Tumour": "#F08A82",
    "Myeloid": "#E5B36D",
    "Lymphoid": "#9DCD79",
    "Vascular": "#7AAED6",
}


def draw_group_hull(ax, points: list[tuple[float, float]], color: str) -> None:
    geom = MultiPoint(points).convex_hull.buffer(0.20, resolution=32)
    if geom.geom_type == "Polygon":
        coords = np.asarray(geom.exterior.coords)
        ax.add_patch(Polygon(coords, closed=True, facecolor=color, edgecolor=color,
                             alpha=0.15, linewidth=0.85, zorder=0))


def draw_modular_network(ax, matrix: pd.DataFrame) -> None:
    base.clean_axis(ax)
    edges = []
    for target in base.COMM_NODES:
        for source in base.COMM_NODES:
            value = float(matrix.loc[target, source])
            if source != target and value > 0:
                edges.append((source, target, value))
    edges = sorted(edges, key=lambda item: item[2], reverse=True)[:30]
    graph = nx.DiGraph()
    graph.add_nodes_from(base.COMM_NODES)
    for source, target, value in edges:
        graph.add_edge(source, target, weight=np.sqrt(value))
    layout = nx.spring_layout(
        graph.to_undirected(), seed=base.SEED, weight="weight", k=0.72, iterations=350
    )
    pos = {node: np.asarray(layout[node], dtype=float) for node in base.COMM_NODES}

    for group, nodes in GROUPS.items():
        draw_group_hull(ax, [tuple(pos[node]) for node in nodes], GROUP_COLORS[group])

    vmax = max([value for _, _, value in edges], default=1.0)
    for rank, (source, target, value) in enumerate(reversed(edges)):
        p0, p1 = pos[source], pos[target]
        sign = -1 if (base.COMM_NODES.index(source) + base.COMM_NODES.index(target)) % 2 else 1
        rad = sign * (0.05 + 0.015 * (rank % 4))
        ax.add_patch(FancyArrowPatch(
            p0, p1, arrowstyle="-|>", mutation_scale=5.0,
            connectionstyle=f"arc3,rad={rad}", shrinkA=12, shrinkB=12,
            color="#66717B", linewidth=0.35 + 1.55 * np.sqrt(value / vmax),
            alpha=0.25 + 0.32 * value / vmax, zorder=1,
        ))

    node_group = {node: group for group, nodes in GROUPS.items() for node in nodes}
    totals = matrix.sum(axis=0) + matrix.sum(axis=1)
    for node in base.COMM_NODES:
        xy = pos[node]
        size = 95 + 170 * float(totals[node] / max(totals.max(), 1e-12))
        color = GROUP_COLORS[node_group[node]]
        ax.scatter(*xy, s=size, facecolor=color, edgecolor=base.INK, linewidth=0.9, zorder=5)
        direction = xy / max(np.linalg.norm(xy), 1e-9)
        label_xy = xy + 0.12 * direction
        ax.text(
            label_xy[0], label_xy[1], base.SHORT_NODE[node],
            ha="center", va="center", fontsize=6.8, color=base.INK, zorder=7,
        )

    all_xy = np.vstack(list(pos.values()))
    xmin, ymin = all_xy.min(axis=0)
    xmax, ymax = all_xy.max(axis=0)
    ax.set_xlim(xmin - 0.48, xmax + 0.48)
    ax.set_ylim(ymin - 0.43, ymax + 0.43)
    ax.set_aspect("equal")


def plot_s5a_modular(interactions: pd.DataFrame, output_dir: Path) -> dict[str, str]:
    matrices = {}
    units = {}
    positive = []
    for dataset in base.DATASETS:
        matrices[dataset], units[dataset] = base.aggregate_communication(interactions, dataset)
        vals = matrices[dataset].to_numpy()
        positive.extend(vals[vals > 0])
    vmin = max(float(np.quantile(positive, 0.04)), 1e-7)
    vmax = float(np.quantile(positive, 0.97))
    cmap = base.LinearSegmentedColormap.from_list(
        "comm2", ["#F7FAFB", "#C7D9E7", "#688DB4", "#253C78"]
    )

    fig = plt.figure(figsize=(12.8, 12.0))
    gs = fig.add_gridspec(3, 2, width_ratios=[1.04, 0.96], hspace=0.43, wspace=0.13)
    image = None
    for i, dataset in enumerate(base.DATASETS):
        matrix = matrices[dataset]
        axh = fig.add_subplot(gs[i, 0])
        image = axh.imshow(
            matrix.to_numpy(), cmap=cmap,
            norm=base.LogNorm(vmin=vmin, vmax=max(vmax, vmin * 1.1)), aspect="equal",
        )
        axh.set_xticks(
            range(len(base.COMM_NODES)), [base.SHORT_NODE[x] for x in base.COMM_NODES],
            rotation=42, ha="right", fontsize=7.2,
        )
        axh.set_yticks(
            range(len(base.COMM_NODES)), [base.SHORT_NODE[x] for x in base.COMM_NODES],
            fontsize=7.2,
        )
        axh.set_xlabel("Sender", color=base.MUTED)
        axh.set_ylabel("Receiver", color=base.MUTED)
        axh.set_title(
            f"{dataset} · {units[dataset]} patient-sample units",
            loc="left", fontweight="bold", color=base.INK,
        )
        for spine in axh.spines.values():
            spine.set_visible(False)

        axn = fig.add_subplot(gs[i, 1])
        draw_modular_network(axn, matrix)
        axn.set_title(
            "Module-aware interaction graph", fontsize=9, color=base.MUTED, pad=2
        )

    if image is not None:
        cax = fig.add_axes([0.43, 0.035, 0.22, 0.014])
        cb = fig.colorbar(image, cax=cax, orientation="horizontal")
        cb.set_label("Patient-median CellChat probability (log scale)", fontsize=8)
        cb.ax.tick_params(labelsize=7)
    legend_handles = [
        Line2D([0], [0], marker="o", lw=0, markersize=7, markerfacecolor=GROUP_COLORS[g],
               markeredgecolor=base.INK, label=g) for g in GROUPS
    ]
    fig.legend(
        handles=legend_handles, loc="lower right", bbox_to_anchor=(0.96, 0.025),
        frameon=False, fontsize=8, ncol=2,
    )
    fig.suptitle(
        "Cohort-separated CellChat heatmaps and modular interaction graphs",
        fontsize=16, fontweight="bold", color=base.INK, y=0.995,
    )
    fig.text(
        0.5, 0.965,
        "Network positions are spring-layout modules with coloured group envelopes and weighted directed edges; no circular layout is used",
        ha="center", fontsize=9.1, color=base.MUTED,
    )
    fig.subplots_adjust(left=0.075, right=0.965, top=0.93, bottom=0.08)
    return base.save_all(fig, output_dir, "02_s5a_modular_cellchat_network_v2")


def categorical_umap(ax, frame: pd.DataFrame, column: str, palette: dict, title: str,
                     legend: bool = True, legend_cols: int = 2) -> None:
    for category, color in palette.items():
        sub = frame[frame[column].eq(category)]
        if len(sub):
            ax.scatter(
                sub["UMAP_1"], sub["UMAP_2"], s=1.25, color=color, alpha=0.64,
                linewidths=0, rasterized=True, label=str(category),
            )
    ax.set_title(title, fontweight="bold", color=base.INK, pad=5)
    ax.set_aspect("equal", adjustable="datalim")
    base.clean_axis(ax)
    if legend:
        ax.legend(
            frameon=False, fontsize=6.5, ncol=legend_cols, loc="lower center",
            bbox_to_anchor=(0.5, -0.13), handletextpad=0.25, columnspacing=0.7,
            markerscale=4,
        )


def plot_s3a_multiannotation(atlas_frames, clinical_path: Path, output_dir: Path) -> dict[str, str]:
    frame = atlas_frames["GSE294109"].copy()
    clinical = pd.read_csv(
        clinical_path, sep="\t",
        usecols=["sample_id", "patient_id", "timepoint", "treatment", "response", "tissue_confirmed"],
    ).drop_duplicates("sample_id")
    frame = frame.merge(clinical, on=["sample_id", "patient_id"], how="left", validate="many_to_one")
    frame["benefit_group"] = frame["response"].map(
        {"CR": "Benefit (CR)", "PR": "Benefit (PR)", "No_effect": "No effect", "Unknown": "Unknown"}
    ).fillna("Unknown")
    frame["tissue_short"] = frame["tissue_confirmed"].replace({
        "Primary Tumor (Multiple regions)": "Primary tumour · multi-region",
        "Primary Tumor": "Primary tumour",
        "Recurrent metastasis (vena cava)": "Recurrent metastasis",
    })

    response_palette = {
        "Benefit (CR)": "#2B8CBE", "Benefit (PR)": "#63B7AF",
        "No effect": "#C4376E", "Unknown": "#B9BDC2",
    }
    treatment_palette = {
        "Untreated": "#9E9E9E", "Ipilimumab/Nivolumab": "#6A51A3",
        "Axitinib/Pembrolizumab": "#E6550D", "Lenvatinib/Pembrolizumab": "#31A354",
        "Cabozantinib": "#D4A72C",
    }
    tissue_palette = {
        "Primary tumour · multi-region": "#4C78A8",
        "Primary tumour": "#72B7B2",
        "Recurrent metastasis": "#E45756",
    }
    patients = sorted(frame["patient_id"].dropna().unique())
    tab20 = plt.get_cmap("tab20")
    patient_palette = {patient: tab20(i % 20) for i, patient in enumerate(patients)}

    fig, axes = plt.subplots(2, 3, figsize=(13.8, 9.1))
    categorical_umap(
        axes[0, 0], frame, "phase8_state",
        {s: base.STATE_COLORS[s] for s in base.STATE_ORDER},
        "Frozen cell state", legend=False,
    )
    categorical_umap(
        axes[0, 1], frame, "benefit_group", response_palette,
        "Immunotherapy outcome", legend=True, legend_cols=2,
    )
    categorical_umap(
        axes[0, 2], frame, "treatment", treatment_palette,
        "Treatment", legend=True, legend_cols=2,
    )
    categorical_umap(
        axes[1, 0], frame, "patient_id", patient_palette,
        f"Patient of origin · {len(patients)} patients", legend=False,
    )
    categorical_umap(
        axes[1, 1], frame, "tissue_short", tissue_palette,
        "Tissue context", legend=True, legend_cols=1,
    )

    legend_ax = axes[1, 2]
    legend_ax.axis("off")
    legend_ax.text(
        0.02, 0.98, "Cell-state key", ha="left", va="top",
        fontsize=10, fontweight="bold", color=base.INK, transform=legend_ax.transAxes,
    )
    state_handles = [
        Line2D([0], [0], marker="o", lw=0, markersize=5,
               markerfacecolor=base.STATE_COLORS[state], markeredgecolor="none",
               label=state.replace("-state Macrophage", " Mφ"))
        for state in base.STATE_ORDER
    ]
    first = legend_ax.legend(
        handles=state_handles, frameon=False, fontsize=7.1, ncol=2,
        loc="upper left", bbox_to_anchor=(0.00, 0.92), handletextpad=0.25, columnspacing=0.7,
    )
    legend_ax.add_artist(first)
    legend_ax.text(
        0.02, 0.39, "Patient colour key", ha="left", va="top",
        fontsize=10, fontweight="bold", color=base.INK, transform=legend_ax.transAxes,
    )
    patient_handles = [
        Line2D([0], [0], marker="o", lw=0, markersize=5,
               markerfacecolor=patient_palette[p], markeredgecolor="none", label=p)
        for p in patients
    ]
    legend_ax.legend(
        handles=patient_handles, frameon=False, fontsize=7.0, ncol=3,
        loc="lower left", bbox_to_anchor=(0.00, 0.01), handletextpad=0.25, columnspacing=0.65,
    )

    fig.suptitle(
        "GSE294109 single-cell atlas across clinical and biological annotations",
        fontsize=16, fontweight="bold", color=base.INK, y=0.99,
    )
    fig.text(
        0.5, 0.952,
        f"Same cohort-specific UMAP coordinates in every view · {len(frame):,} visualized cells · no cross-cohort integration",
        ha="center", fontsize=9.2, color=base.MUTED,
    )
    fig.subplots_adjust(left=0.035, right=0.985, top=0.90, bottom=0.07, hspace=0.30, wspace=0.14)
    return base.save_all(fig, output_dir, "03_s3a_gse294109_multiannotation_umap")


def plot_s3g_infercnv(infercnv_png: Path, calls_tsv: Path, output_dir: Path) -> dict[str, str]:
    from PIL import Image

    source = Image.open(infercnv_png).convert("RGB")
    width, height = source.size
    top_crop = np.asarray(source.crop((
        int(0.17 * width), int(0.07 * height), int(0.965 * width), int(0.325 * height)
    )))
    bottom_crop = np.asarray(source.crop((
        int(0.17 * width), int(0.335 * height), int(0.965 * width), int(0.755 * height)
    )))
    calls = pd.read_csv(calls_tsv, sep="\t")
    unit = calls[calls["infercnv_unit"].eq("GSE156632__RCC2")].copy()
    n_query = int(len(unit))
    n_high = int(unit["infercnv_call"].eq("high_cnv").sum())
    threshold = float(unit["high_cnv_threshold"].dropna().iloc[0])

    fig = plt.figure(figsize=(14.2, 8.4))
    gs = fig.add_gridspec(
        2, 2, width_ratios=[0.055, 0.945], height_ratios=[0.34, 0.66],
        hspace=0.055, wspace=0.018,
    )
    side_top = fig.add_subplot(gs[0, 0])
    side_top.set_facecolor("#77BFB6")
    side_top.text(
        0.5, 0.5, "Normal reference", rotation=90, ha="center", va="center",
        fontsize=9, fontweight="bold", color=base.INK,
    )
    base.clean_axis(side_top)

    side_bottom = fig.add_subplot(gs[1, 0])
    side_bottom.set_facecolor("#E6A23C")
    side_bottom.text(
        0.5, 0.5, "RCC2 malignant query", rotation=90, ha="center", va="center",
        fontsize=9, fontweight="bold", color=base.INK,
    )
    base.clean_axis(side_bottom)

    ax_top = fig.add_subplot(gs[0, 1])
    ax_top.imshow(top_crop, aspect="auto")
    base.clean_axis(ax_top)
    ax_top.set_title(
        "Reference profile", loc="left", fontsize=10, fontweight="bold", color=base.INK, pad=3
    )

    ax_bottom = fig.add_subplot(gs[1, 1])
    ax_bottom.imshow(bottom_crop, aspect="auto")
    base.clean_axis(ax_bottom)
    ax_bottom.set_title(
        f"Query profile · {n_query} cells · {n_high}/{n_query} high-CNV calls",
        loc="left", fontsize=10, fontweight="bold", color=base.INK, pad=3,
    )
    ax_bottom.set_xlabel("Chromosome-ordered genomic regions", fontsize=9, color=base.INK)

    fig.suptitle(
        "Chromosome-ordered inferCNV profile — GSE156632 RCC2",
        fontsize=16, fontweight="bold", color=base.INK, y=0.985,
    )
    fig.text(
        0.5, 0.942,
        f"Existing inferCNV output reframed in the S3G visual form · high-CNV threshold = {threshold:,.1f} genes",
        ha="center", fontsize=9.2, color=base.MUTED,
    )
    fig.text(
        0.985, 0.018,
        "Red and blue indicate relative gains and losses after smoothing; this legacy malignant-CNV panel is outside the current immune-only report scope.",
        ha="right", fontsize=8.2, color=base.MUTED,
    )
    fig.subplots_adjust(left=0.035, right=0.985, top=0.90, bottom=0.08)
    return base.save_all(fig, output_dir, "04_s3g_infercnv_gse156632_rcc2")


def validate_outputs(outputs: dict[str, dict[str, str]], output_dir: Path) -> dict:
    checks = []
    for figure, paths in outputs.items():
        for fmt, text in paths.items():
            path = Path(text)
            checks.append({
                "figure": figure, "format": fmt, "exists": path.exists(),
                "bytes": path.stat().st_size if path.exists() else 0,
            })
    passed = all(row["exists"] and row["bytes"] > 5000 for row in checks)
    result = {"status": "PASS" if passed else "FAIL", "checks": checks}
    (output_dir / "REVISION_FIGURE_VALIDATION_v1.json").write_text(
        json.dumps(result, indent=2), encoding="utf-8"
    )
    if not passed:
        raise SystemExit("Revision figure validation failed")
    return result


def main() -> None:
    if len(sys.argv) != 7:
        raise SystemExit(
            "Usage: make_reference_revisions_v1.py <project_root> <atlas_dir> "
            "<clinical_sample_tsv> <infercnv_png> <infercnv_calls_tsv> <output_dir>"
        )
    project_root = Path(sys.argv[1])
    atlas_dir = Path(sys.argv[2])
    clinical_path = Path(sys.argv[3])
    infercnv_png = Path(sys.argv[4])
    calls_tsv = Path(sys.argv[5])
    output_dir = Path(sys.argv[6])
    base.setup_style()
    atlas_frames, status, pseudotime, concordance, interactions = base.load_inputs(
        project_root, atlas_dir
    )
    outputs = {
        "F4A_fixed": plot_f4a_fixed(atlas_frames, status, output_dir),
        "S5A_modular": plot_s5a_modular(interactions, output_dir),
        "S3A": plot_s3a_multiannotation(atlas_frames, clinical_path, output_dir),
        "S3G": plot_s3g_infercnv(infercnv_png, calls_tsv, output_dir),
    }
    validation = validate_outputs(outputs, output_dir)
    summary = {
        "status": validation["status"],
        "outputs": outputs,
        "review_decision": {
            "retain_from_previous_set": ["F4D myeloid topology", "S7A LR priority map"],
            "replace": ["F4A covariance ellipses", "S5A circular network"],
            "drop": ["F5D simplex"],
            "add": ["S3A multi-annotation UMAP", "S3G inferCNV"],
        },
        "boundaries": [
            "F4A uses local density envelopes with disconnected contours.",
            "S5A uses module-aware spring layouts, not circular networks.",
            "S3A is GSE294109 only and uses cohort-specific UMAP coordinates.",
            "S3G reuses a real legacy GSE156632 RCC2 inferCNV output and remains outside the immune-only report scope.",
        ],
    }
    (output_dir / "REVISION_FIGURE_SUMMARY_v1.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    print(json.dumps({
        "status": validation["status"], "figures": len(outputs), "output_dir": str(output_dir)
    }, indent=2))


if __name__ == "__main__":
    main()
