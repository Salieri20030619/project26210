from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, LogNorm
from matplotlib.lines import Line2D
from matplotlib.patches import Ellipse, FancyArrowPatch, Polygon


SEED = 26210
RNG = np.random.default_rng(SEED)
INK = "#17202A"
MUTED = "#65717D"
GRID = "#D9DEE3"
TEAL = "#178F8A"
MAGENTA = "#C22F68"
GOLD = "#E4A11B"

DATASETS = ["GSE207493", "GSE264586", "GSE294109"]
DATASET_COLORS = {"GSE207493": "#2B8CBE", "GSE264586": "#138F87", "GSE294109": "#C4376E"}
DATASET_MARKERS = {"GSE207493": "o", "GSE264586": "s", "GSE294109": "D"}

STATE_ORDER = [
    "F1-high Malignant", "F1-low Malignant", "Epithelial-Normal",
    "Monocyte", "SPP1-state Macrophage", "C1QC-state Macrophage", "Other Macrophage",
    "Dendritic cell", "T cell", "NK cell", "B cell", "Mast cell",
    "Endothelial", "SMC/Pericyte",
]
STATE_COLORS = {
    "F1-high Malignant": "#E15D44", "F1-low Malignant": "#F3A36B",
    "Epithelial-Normal": "#D6B45A", "Monocyte": "#63B7AF",
    "SPP1-state Macrophage": "#D64A78", "C1QC-state Macrophage": "#7356A5",
    "Other Macrophage": "#3C7C91", "Dendritic cell": "#84B6D9",
    "T cell": "#5B4B9A", "NK cell": "#8C6BB1", "B cell": "#B99AC8",
    "Mast cell": "#C98B68", "Endothelial": "#4CA46C", "SMC/Pericyte": "#8EBE7A",
}
MYELOID_STATES = ["Monocyte", "SPP1-state Macrophage", "C1QC-state Macrophage", "Other Macrophage"]
MAC_STATES = ["SPP1-state Macrophage", "C1QC-state Macrophage", "Other Macrophage"]
COMM_NODES = [
    "Malignant", "Monocyte", "SPP1-state Macrophage", "C1QC-state Macrophage",
    "Other Macrophage", "T cell", "NK cell", "Endothelial",
]
COMM_COLORS = {
    "Malignant": "#E15D44", "Monocyte": "#63B7AF", "SPP1-state Macrophage": "#D64A78",
    "C1QC-state Macrophage": "#7356A5", "Other Macrophage": "#3C7C91",
    "T cell": "#5B4B9A", "NK cell": "#8C6BB1", "Endothelial": "#4CA46C",
}
SHORT_NODE = {
    "Malignant": "Tumour", "Monocyte": "Monocyte", "SPP1-state Macrophage": "SPP1 Mφ",
    "C1QC-state Macrophage": "C1QC Mφ", "Other Macrophage": "Other Mφ",
    "T cell": "T cells", "NK cell": "NK", "Endothelial": "Endo",
}


def setup_style() -> None:
    mpl.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 9, "axes.titlesize": 11,
        "axes.labelsize": 9, "axes.linewidth": 0.7, "pdf.fonttype": 42,
        "ps.fonttype": 42, "svg.fonttype": "none", "savefig.facecolor": "white",
    })


def clean_axis(ax) -> None:
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)


def save_all(fig, output_dir: Path, stem: str) -> dict[str, str]:
    output_dir.mkdir(parents=True, exist_ok=True)
    paths = {}
    for ext in ["png", "pdf", "svg"]:
        path = output_dir / f"{stem}.{ext}"
        fig.savefig(path, dpi=320, bbox_inches="tight", pad_inches=0.08)
        paths[ext] = str(path)
    plt.close(fig)
    return paths


def add_covariance_halo(ax, frame: pd.DataFrame, color: str, lw: float = 0.8) -> None:
    if len(frame) < 60:
        return
    points = frame[["UMAP_1", "UMAP_2"]].dropna().to_numpy(float)
    if len(points) > 5000:
        points = points[RNG.choice(len(points), 5000, replace=False)]
    center = np.median(points, axis=0)
    centered = points - center
    dist = np.sqrt(np.sum(centered ** 2, axis=1))
    points = points[dist <= np.quantile(dist, 0.90)]
    if len(points) < 30:
        return
    cov = np.cov(points.T)
    vals, vecs = np.linalg.eigh(cov)
    order = vals.argsort()[::-1]
    vals, vecs = vals[order], vecs[:, order]
    angle = np.degrees(np.arctan2(*vecs[:, 0][::-1]))
    width, height = 4.0 * np.sqrt(np.maximum(vals, 1e-6))
    ax.add_patch(Ellipse(center, width, height, angle=angle, fill=False, edgecolor=color,
                         linewidth=lw, linestyle=(0, (2.2, 2.2)), alpha=0.72, zorder=8))


def load_inputs(project_root: Path, atlas_dir: Path):
    atlas_frames = {
        dataset: pd.read_csv(atlas_dir / f"{dataset}_atlas_umap_sample_v1.tsv.gz", sep="\t")
        for dataset in DATASETS
    }
    status = pd.read_csv(atlas_dir / "atlas_extraction_status_v1.tsv", sep="\t")
    phase87 = project_root / "phase87_three_cohort_myeloid_trajectory_local_v1" / "results_v1"
    pseudotime = pd.read_csv(phase87 / "phase87_cell_pseudotime_v1.tsv.gz", sep="\t")
    concordance = pd.read_csv(phase87 / "phase87_dataset_method_concordance_summary_v1.tsv", sep="\t")
    interactions = pd.read_csv(
        project_root / "phase88_three_cohort_fine_state_cellchat_local_v1" / "results_v1" /
        "phase88_cellchat_interactions_all_v1.tsv.gz", sep="\t"
    )
    return atlas_frames, status, pseudotime, concordance, interactions


def plot_atlas_constellation(atlas_frames, status, output_dir: Path) -> dict[str, str]:
    fig, axes = plt.subplots(1, 3, figsize=(14.5, 5.8), constrained_layout=False)
    for ax, dataset in zip(axes, DATASETS):
        frame = atlas_frames[dataset]
        for state in STATE_ORDER:
            sub = frame[frame["phase8_state"].eq(state)]
            if len(sub):
                ax.scatter(sub["UMAP_1"], sub["UMAP_2"], s=1.25, c=STATE_COLORS[state],
                           alpha=0.58, linewidths=0, rasterized=True)
                add_covariance_halo(ax, sub, STATE_COLORS[state], lw=0.55)
        row = status.loc[status["dataset"].eq(dataset)].iloc[0]
        ax.set_title(f"{dataset}\n{int(row['patients'])} patients · {int(row['total_cells']):,} cells",
                     fontweight="bold", color=INK, pad=8)
        ax.set_aspect("equal", adjustable="datalim")
        clean_axis(ax)
    fig.suptitle("Cohort-separated single-cell atlas constellation", x=0.5, y=0.985,
                 fontsize=16, fontweight="bold", color=INK)
    fig.text(0.5, 0.915,
             "Existing UMAP coordinates retained within each dataset · dashed halos summarize state topology · no pooled integration",
             ha="center", va="center", fontsize=9.2, color=MUTED)
    handles = [
        Line2D([0], [0], marker="o", lw=0, markersize=6, markerfacecolor=STATE_COLORS[s],
               markeredgecolor="none", label=s.replace("-state Macrophage", " Mφ"))
        for s in STATE_ORDER
    ]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.015), ncol=7,
               frameon=False, fontsize=8, handletextpad=0.35, columnspacing=0.9)
    fig.subplots_adjust(left=0.025, right=0.99, top=0.80, bottom=0.16, wspace=0.05)
    return save_all(fig, output_dir, "01_f4a_three_cohort_atlas_constellation")


def trajectory_path(frame: pd.DataFrame, bins: int = 22) -> pd.DataFrame:
    x = frame.dropna(subset=["slingshot_pseudotime", "UMAP_1", "UMAP_2"]).copy()
    if len(x) < 100:
        return pd.DataFrame()
    x["bin"] = pd.cut(x["slingshot_pseudotime"], bins=np.linspace(0, 1, bins + 1),
                      labels=False, include_lowest=True)
    return x.groupby("bin", observed=True)[["UMAP_1", "UMAP_2", "slingshot_pseudotime"]].median().dropna().reset_index()


def plot_myeloid_topology(atlas_frames, pseudotime, concordance, output_dir: Path) -> dict[str, str]:
    fig, axes = plt.subplots(1, 3, figsize=(13.8, 5.3))
    for ax, dataset in zip(axes, DATASETS):
        atlas = atlas_frames[dataset]
        x = atlas[atlas["phase8_state"].isin(MYELOID_STATES)].merge(
            pseudotime[pseudotime["dataset"].eq(dataset)][["cell_id", "slingshot_pseudotime"]],
            on="cell_id", how="inner",
        )
        for state in MYELOID_STATES:
            sub = x[x["phase8_state"].eq(state)]
            ax.scatter(sub["UMAP_1"], sub["UMAP_2"], s=2.1, c=STATE_COLORS[state],
                       alpha=0.64, linewidths=0, rasterized=True)
            add_covariance_halo(ax, sub, STATE_COLORS[state], lw=1.0)
        path = trajectory_path(x)
        if len(path) > 3:
            ax.plot(path["UMAP_1"], path["UMAP_2"], color=INK, lw=1.8, zorder=10)
            ax.add_patch(FancyArrowPatch(
                (path.iloc[-3]["UMAP_1"], path.iloc[-3]["UMAP_2"]),
                (path.iloc[-1]["UMAP_1"], path.iloc[-1]["UMAP_2"]),
                arrowstyle="-|>", mutation_scale=12, color=INK, lw=1.5, zorder=11
            ))
        rho = float(concordance.loc[concordance["dataset"].eq(dataset), "median_patient_rho"].iloc[0])
        n_pat = int(pseudotime.loc[pseudotime["dataset"].eq(dataset), "patient_id"].nunique())
        n_cell = int((pseudotime["dataset"].eq(dataset)).sum())
        ax.set_title(f"{dataset}\n{n_pat} patients · {n_cell:,} cells · method ρ={rho:.2f}",
                     fontweight="bold", color=INK, pad=8)
        ax.set_aspect("equal", adjustable="datalim")
        qx = x["UMAP_1"].quantile([0.005, 0.995]).to_numpy()
        qy = x["UMAP_2"].quantile([0.005, 0.995]).to_numpy()
        ax.set_xlim(qx[0] - 0.04 * (qx[1] - qx[0]), qx[1] + 0.04 * (qx[1] - qx[0]))
        ax.set_ylim(qy[0] - 0.04 * (qy[1] - qy[0]), qy[1] + 0.04 * (qy[1] - qy[0]))
        clean_axis(ax)
    fig.suptitle("Myeloid-state topology and reproducible ordering", fontsize=16,
                 fontweight="bold", color=INK, y=0.985)
    fig.text(
        0.5, 0.91,
        "Dashed halos summarize state neighborhoods; black arrows show increasing Slingshot pseudotime (cross-sectional ordering, not lineage proof)",
        ha="center", color=MUTED, fontsize=9.1
    )
    handles = [
        Line2D([0], [0], marker="o", lw=0, markersize=7, markerfacecolor=STATE_COLORS[s],
               markeredgecolor="none", label=s.replace("-state Macrophage", " Mφ"))
        for s in MYELOID_STATES
    ]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 0.01), ncol=4,
               frameon=False, fontsize=9)
    fig.subplots_adjust(left=0.03, right=0.985, top=0.79, bottom=0.16, wspace=0.08)
    return save_all(fig, output_dir, "02_f4d_myeloid_state_topology")


def ternary_xy(a, b, c):
    total = np.asarray(a) + np.asarray(b) + np.asarray(c)
    total = np.where(total == 0, np.nan, total)
    a, b, c = np.asarray(a) / total, np.asarray(b) / total, np.asarray(c) / total
    x = b + 0.5 * c
    y = (math.sqrt(3) / 2) * c
    return x, y


def ternary_segment(p1, p2):
    x1, y1 = ternary_xy(*p1)
    x2, y2 = ternary_xy(*p2)
    return [float(x1), float(x2)], [float(y1), float(y2)]


def draw_ternary(ax, data: pd.DataFrame, title: str, cohorts: list[str], show_legend: bool = False, show_vertex_labels: bool = True) -> None:
    height = math.sqrt(3) / 2
    ax.plot([0, 1, 0.5, 0], [0, 0, height, 0], color=INK, lw=1.1)
    for frac in [0.2, 0.4, 0.6, 0.8]:
        ax.plot(*ternary_segment((frac, 1-frac, 0), (frac, 0, 1-frac)), color=GRID, lw=0.55)
        ax.plot(*ternary_segment((1-frac, frac, 0), (0, frac, 1-frac)), color=GRID, lw=0.55)
        ax.plot(*ternary_segment((1-frac, 0, frac), (0, 1-frac, frac)), color=GRID, lw=0.55)
    for dataset in cohorts:
        sub = data[data["dataset"].eq(dataset)].copy()
        if sub.empty:
            continue
        x, y = ternary_xy(sub["SPP1"], sub["C1QC"], sub["Other"])
        sizes = 26 + 90 * np.sqrt(sub["n_cells"] / data["n_cells"].max())
        ax.scatter(x, y, s=sizes, c=DATASET_COLORS[dataset], marker=DATASET_MARKERS[dataset],
                   alpha=0.62, edgecolors="white", linewidths=0.55, label=dataset, zorder=4)
        med = sub[["SPP1", "C1QC", "Other"]].median()
        mx, my = ternary_xy(med["SPP1"], med["C1QC"], med["Other"])
        ax.scatter(mx, my, s=150, marker="*", c=DATASET_COLORS[dataset], edgecolors=INK,
                   linewidths=0.7, zorder=7)
        if len(sub) >= 3:
            pts = np.column_stack([x, y])
            center = pts.mean(axis=0)
            angles = np.arctan2(pts[:, 1]-center[1], pts[:, 0]-center[0])
            poly = pts[np.argsort(angles)]
            ax.add_patch(Polygon(poly, closed=True, facecolor=DATASET_COLORS[dataset],
                                 edgecolor=DATASET_COLORS[dataset], alpha=0.07, lw=0.8, zorder=1))
    if show_vertex_labels:
        ax.text(-0.035, -0.035, "SPP1 Mφ", ha="right", va="top", fontweight="bold", color=INK)
        ax.text(1.035, -0.035, "C1QC Mφ", ha="left", va="top", fontweight="bold", color=INK)
        ax.text(0.5, height + 0.04, "Other Mφ", ha="center", va="bottom", fontweight="bold", color=INK)
    ax.set_xlim(-0.08, 1.08)
    ax.set_ylim(-0.08, height + 0.09)
    ax.set_aspect("equal")
    clean_axis(ax)
    ax.set_title(title, fontweight="bold", color=INK, pad=4)
    if show_legend:
        ax.legend(frameon=False, fontsize=8, loc="lower center", bbox_to_anchor=(0.5, -0.20), ncol=3)


def plot_macrophage_simplex(pseudotime, output_dir: Path) -> dict[str, str]:
    mac = pseudotime[pseudotime["phase8_state"].isin(MAC_STATES)].copy()
    counts = mac.groupby(["dataset", "patient_id", "phase8_state"], observed=True).size().unstack(fill_value=0)
    for state in MAC_STATES:
        if state not in counts:
            counts[state] = 0
    counts = counts.reset_index()
    counts["n_cells"] = counts[MAC_STATES].sum(axis=1)
    counts = counts[counts["n_cells"].ge(30)].copy()
    counts["SPP1"] = counts["SPP1-state Macrophage"] / counts["n_cells"]
    counts["C1QC"] = counts["C1QC-state Macrophage"] / counts["n_cells"]
    counts["Other"] = counts["Other Macrophage"] / counts["n_cells"]

    fig = plt.figure(figsize=(10.8, 9.2))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.12, 0.88], hspace=0.25, wspace=0.08)
    top = fig.add_subplot(gs[0, :])
    draw_ternary(top, counts, "All eligible patients", DATASETS, show_legend=True)
    for ax, dataset in zip([fig.add_subplot(gs[1, i]) for i in range(3)], DATASETS):
        n = int(counts["dataset"].eq(dataset).sum())
        draw_ternary(ax, counts, f"{dataset} · n={n}", [dataset], show_vertex_labels=False)
    fig.suptitle("Patient-level macrophage-state composition simplex", fontsize=16,
                 fontweight="bold", color=INK, y=0.985)
    fig.text(
        0.5, 0.94,
        "Each point is one patient; coordinates are proportions within SPP1/C1QC/Other macrophage cells · star = cohort median · size = eligible cells",
        ha="center", color=MUTED, fontsize=9.2
    )
    return save_all(fig, output_dir, "03_f5d_patient_macrophage_state_simplex")


def aggregate_communication(interactions: pd.DataFrame, dataset: str) -> tuple[pd.DataFrame, int]:
    x = interactions[
        interactions["dataset"].eq(dataset) & interactions["source"].isin(COMM_NODES) &
        interactions["target"].isin(COMM_NODES) & interactions["pval"].le(0.05)
    ].copy()
    n_units = int(interactions.loc[interactions["dataset"].eq(dataset), "unit_id"].nunique())
    unit = x.groupby(["unit_id", "source", "target"], observed=True)["prob"].sum().reset_index()
    edge = unit.groupby(["source", "target"], observed=True)["prob"].median().reset_index()
    matrix = edge.pivot(index="target", columns="source", values="prob").reindex(
        index=COMM_NODES, columns=COMM_NODES
    ).fillna(0)
    return matrix, n_units


def draw_comm_network(ax, matrix: pd.DataFrame) -> None:
    clean_axis(ax)
    angles = np.linspace(np.pi/2, np.pi/2 + 2*np.pi, len(COMM_NODES), endpoint=False)
    pos = {node: np.array([math.cos(a), math.sin(a)]) for node, a in zip(COMM_NODES, angles)}
    edges = []
    for target in COMM_NODES:
        for source in COMM_NODES:
            if source != target and matrix.loc[target, source] > 0:
                edges.append((source, target, float(matrix.loc[target, source])))
    edges = sorted(edges, key=lambda z: z[2], reverse=True)[:18]
    vmax = max([v for _, _, v in edges], default=1.0)
    for source, target, value in reversed(edges):
        p0, p1 = pos[source], pos[target]
        rad = 0.18 if COMM_NODES.index(source) < COMM_NODES.index(target) else -0.18
        ax.add_patch(FancyArrowPatch(
            p0*0.82, p1*0.82, arrowstyle="-|>", mutation_scale=6.5,
            connectionstyle=f"arc3,rad={rad}", color=COMM_COLORS[source],
            lw=0.4 + 2.7*np.sqrt(value/vmax), alpha=0.24 + 0.46*value/vmax, zorder=1
        ))
    totals = matrix.sum(axis=0) + matrix.sum(axis=1)
    for node in COMM_NODES:
        xy = pos[node]
        size = 100 + 240 * float(totals[node] / max(totals.max(), 1e-12))
        ax.scatter(*xy, s=size, color=COMM_COLORS[node], edgecolors="white", linewidths=1.0, zorder=5)
        ax.text(xy[0]*1.28, xy[1]*1.28, SHORT_NODE[node], ha="center", va="center", fontsize=6.8, color=INK)
    ax.set_xlim(-1.55, 1.55)
    ax.set_ylim(-1.55, 1.55)
    ax.set_aspect("equal")


def plot_cellchat_heatmap_network(interactions, output_dir: Path) -> dict[str, str]:
    matrices = {}
    units = {}
    positive = []
    for dataset in DATASETS:
        matrices[dataset], units[dataset] = aggregate_communication(interactions, dataset)
        values = matrices[dataset].to_numpy()
        positive.extend(values[values > 0])
    vmin = max(np.quantile(positive, 0.04), 1e-7) if positive else 1e-7
    vmax = np.quantile(positive, 0.97) if positive else 1.0
    cmap = LinearSegmentedColormap.from_list("comm", ["#F7FAFB", "#C7D9E7", "#688DB4", "#253C78"])

    fig = plt.figure(figsize=(12.5, 12.0))
    gs = fig.add_gridspec(3, 2, width_ratios=[1.03, 0.97], hspace=0.43, wspace=0.16)
    image = None
    for i, dataset in enumerate(DATASETS):
        axh = fig.add_subplot(gs[i, 0])
        matrix = matrices[dataset]
        image = axh.imshow(matrix.to_numpy(), cmap=cmap,
                           norm=LogNorm(vmin=vmin, vmax=max(vmax, vmin*1.1)), aspect="equal")
        axh.set_xticks(range(len(COMM_NODES)), [SHORT_NODE[x] for x in COMM_NODES],
                       rotation=42, ha="right", fontsize=7.2)
        axh.set_yticks(range(len(COMM_NODES)), [SHORT_NODE[x] for x in COMM_NODES], fontsize=7.2)
        axh.set_xlabel("Sender", color=MUTED)
        axh.set_ylabel("Receiver", color=MUTED)
        axh.set_title(f"{dataset} · {units[dataset]} patient-sample units",
                      loc="left", fontweight="bold", color=INK)
        for spine in axh.spines.values():
            spine.set_visible(False)
        axn = fig.add_subplot(gs[i, 1])
        draw_comm_network(axn, matrix)
        axn.set_title("Top patient-median non-self edges", fontsize=9, color=MUTED, pad=2)
    if image is not None:
        cax = fig.add_axes([0.43, 0.035, 0.22, 0.014])
        cb = fig.colorbar(image, cax=cax, orientation="horizontal")
        cb.set_label("Patient-median CellChat probability (log scale)", fontsize=8)
        cb.ax.tick_params(labelsize=7)
    fig.suptitle("Cohort-separated CellChat interaction landscapes", fontsize=16,
                 fontweight="bold", color=INK, y=0.995)
    fig.text(
        0.5, 0.965,
        "Left: significant ligand-receptor probabilities aggregated to sender-receiver pairs · right: top directed edges · absolute strength is not fine-state specificity",
        ha="center", fontsize=9.1, color=MUTED
    )
    fig.subplots_adjust(left=0.075, right=0.965, top=0.93, bottom=0.08)
    return save_all(fig, output_dir, "04_s5a_cellchat_heatmap_network")


def plot_ligand_receiver_priority(interactions, output_dir: Path) -> dict[str, str]:
    targets = MAC_STATES
    x = interactions[
        interactions["source"].eq("Malignant") & interactions["target"].isin(targets)
    ].copy()
    totals = x.groupby("dataset", observed=True)["unit_id"].nunique().to_dict()
    sig = x[x["pval"].le(0.05)].copy()
    unit_lig = sig.groupby(["dataset", "unit_id", "ligand"], observed=True)["prob"].sum().reset_index()
    ligand_stats = unit_lig.groupby(["dataset", "ligand"], observed=True).agg(
        median_prob=("prob", "median"), present_units=("unit_id", "nunique")
    ).reset_index()
    ligand_stats["prevalence"] = ligand_stats.apply(
        lambda r: r["present_units"] / totals.get(r["dataset"], 1), axis=1
    )
    score = ligand_stats.groupby("ligand", observed=True).agg(
        datasets=("dataset", "nunique"),
        median_prob=("median_prob", "mean"),
        prevalence=("prevalence", "mean"),
    )
    score = score[score["datasets"].ge(2)].copy()
    score["priority"] = np.log10(score["median_prob"].clip(lower=1e-12)) + score["prevalence"]
    ligands = score.sort_values("priority", ascending=False).head(22).index.tolist()
    if len(ligands) < 12:
        ligands = ligand_stats.groupby("ligand")["median_prob"].mean().sort_values(
            ascending=False
        ).head(22).index.tolist()

    feature_unit = sig[sig["ligand"].isin(ligands)].groupby(
        ["dataset", "unit_id", "ligand", "target", "receptor"], observed=True
    )["prob"].sum().reset_index()
    target_short = {
        "SPP1-state Macrophage": "SPP1",
        "C1QC-state Macrophage": "C1QC",
        "Other Macrophage": "Other",
    }
    feature_unit["feature"] = (
        feature_unit["target"].map(target_short) + ":" + feature_unit["receptor"].astype(str)
    )
    feature_score = feature_unit.groupby("feature", observed=True)["prob"].median().sort_values(
        ascending=False
    )
    features = feature_score.head(30).index.tolist()
    matrix = feature_unit[feature_unit["feature"].isin(features)].groupby(
        ["ligand", "feature"], observed=True
    )["prob"].median().unstack(fill_value=0).reindex(index=ligands, columns=features).fillna(0)
    ligands = matrix.index.tolist()

    stats = ligand_stats[ligand_stats["ligand"].isin(ligands)].copy()
    prob = stats.pivot(index="ligand", columns="dataset", values="median_prob").reindex(
        index=ligands, columns=DATASETS
    ).fillna(0)
    prev = stats.pivot(index="ligand", columns="dataset", values="prevalence").reindex(
        index=ligands, columns=DATASETS
    ).fillna(0)
    zprob = np.log10(prob + 1e-8)
    for dataset in DATASETS:
        col = zprob[dataset]
        sd = col.std(ddof=0)
        zprob[dataset] = (col-col.mean())/sd if sd > 0 else 0
    priority = score.reindex(ligands)["priority"].fillna(
        score["priority"].min() if len(score) else 0
    )
    pnorm = (priority-priority.min()) / max(priority.max()-priority.min(), 1e-12)

    fig = plt.figure(figsize=(15.0, 8.8))
    gs = fig.add_gridspec(1, 3, width_ratios=[0.42, 0.95, 4.65], wspace=0.035)
    ax_rank = fig.add_subplot(gs[0, 0])
    y = np.arange(len(ligands))
    ax_rank.barh(y, pnorm, color=GOLD, alpha=0.78, height=0.72)
    ax_rank.set_yticks(y, ligands, fontsize=8.0)
    ax_rank.set_ylim(len(ligands)-0.5, -0.5)
    ax_rank.set_xlim(0, 1.02)
    ax_rank.set_xlabel("priority", fontsize=7.5, color=MUTED)
    ax_rank.tick_params(axis="x", labelsize=7)
    for spine in ["top", "right", "left"]:
        ax_rank.spines[spine].set_visible(False)
    ax_rank.spines["bottom"].set_color(GRID)

    ax_bub = fig.add_subplot(gs[0, 1], sharey=ax_rank)
    for j, dataset in enumerate(DATASETS):
        sizes = 14 + 120 * prev[dataset].to_numpy()
        ax_bub.scatter(
            np.full(len(ligands), j), y, s=sizes, c=zprob[dataset].to_numpy(),
            cmap="RdBu_r", vmin=-2.0, vmax=2.0, edgecolors=INK, linewidths=0.35
        )
    ax_bub.set_xticks(
        range(3), [x.replace("GSE", "") for x in DATASETS], rotation=35, ha="right", fontsize=7.5
    )
    ax_bub.tick_params(axis="y", left=False, labelleft=False)
    ax_bub.set_xlim(-0.55, 2.55)
    ax_bub.set_ylim(len(ligands)-0.5, -0.5)
    ax_bub.set_xlabel("cohort", fontsize=7.5, color=MUTED)
    ax_bub.grid(axis="x", color=GRID, lw=0.6)
    for spine in ax_bub.spines.values():
        spine.set_visible(False)

    ax_heat = fig.add_subplot(gs[0, 2], sharey=ax_rank)
    vals = matrix.to_numpy(float)
    positive = vals[vals > 0]
    lo = np.quantile(positive, 0.08) if len(positive) else 1e-7
    hi = np.quantile(positive, 0.98) if len(positive) else 1.0
    cmap = LinearSegmentedColormap.from_list(
        "lig", ["#F5F2F5", "#D5B7D5", "#9A4D9B", "#5B1E63"]
    )
    shown = np.where(vals > 0, np.log10(vals + 1e-10), np.nan)
    ax_heat.imshow(
        shown, aspect="auto", cmap=cmap,
        vmin=np.log10(max(lo, 1e-10)), vmax=np.log10(max(hi, lo*1.1))
    )
    ax_heat.set_xticks(range(len(features)), features, rotation=55, ha="right", fontsize=6.8)
    ax_heat.tick_params(axis="y", left=False, labelleft=False)
    ax_heat.set_ylim(len(ligands)-0.5, -0.5)
    ax_heat.set_xlabel("Receiver state : receptor", fontsize=8.5, color=INK)
    for spine in ax_heat.spines.values():
        spine.set_visible(False)
    ax_heat.set_xticks(np.arange(-0.5, len(features), 1), minor=True)
    ax_heat.set_yticks(np.arange(-0.5, len(ligands), 1), minor=True)
    ax_heat.grid(which="minor", color="white", linewidth=0.75)
    ax_heat.tick_params(which="minor", bottom=False, left=False)

    fig.suptitle(
        "Malignant-to-macrophage CellChat ligand-receptor priority map",
        fontsize=16, fontweight="bold", color=INK, y=0.985
    )
    fig.text(
        0.5, 0.945,
        "Left bar: cross-cohort priority · bubbles: patient prevalence (area) and within-cohort median probability (colour) · matrix: patient-median LR probability",
        ha="center", fontsize=9.2, color=MUTED
    )
    fig.text(
        0.99, 0.018,
        "CellChat expression compatibility only; this is not a NicheNet ligand-target regulatory prediction.",
        ha="right", fontsize=8.2, color=MUTED
    )
    fig.subplots_adjust(left=0.09, right=0.99, top=0.90, bottom=0.20)
    return save_all(fig, output_dir, "05_s7a_malignant_to_macrophage_lr_priority_map")


def validate_outputs(outputs: dict[str, dict[str, str]], output_dir: Path) -> dict:
    checks = []
    for figure, paths in outputs.items():
        for ext, path_text in paths.items():
            path = Path(path_text)
            checks.append({
                "figure": figure,
                "format": ext,
                "exists": path.exists(),
                "bytes": path.stat().st_size if path.exists() else 0,
            })
    passed = all(row["exists"] and row["bytes"] > 5000 for row in checks)
    result = {"status": "PASS" if passed else "FAIL", "checks": checks, "seed": SEED}
    (output_dir / "CANDIDATE_FIGURE_VALIDATION_v1.json").write_text(
        json.dumps(result, indent=2), encoding="utf-8"
    )
    if not passed:
        raise SystemExit("Output validation failed")
    return result


def main() -> None:
    if len(sys.argv) != 4:
        raise SystemExit(
            "Usage: make_single_cell_candidates_v1.py <project_root> <atlas_extract_dir> <output_dir>"
        )
    project_root = Path(sys.argv[1])
    atlas_dir = Path(sys.argv[2])
    output_dir = Path(sys.argv[3])
    setup_style()
    atlas_frames, status, pseudotime, concordance, interactions = load_inputs(
        project_root, atlas_dir
    )
    outputs = {
        "F4A_atlas": plot_atlas_constellation(atlas_frames, status, output_dir),
        "F4D_myeloid": plot_myeloid_topology(
            atlas_frames, pseudotime, concordance, output_dir
        ),
        "F5D_simplex": plot_macrophage_simplex(pseudotime, output_dir),
        "S5A_cellchat": plot_cellchat_heatmap_network(interactions, output_dir),
        "S7A_lr_map": plot_ligand_receiver_priority(interactions, output_dir),
    }
    validation = validate_outputs(outputs, output_dir)
    summary = {
        "status": validation["status"],
        "reference_forms": ["F4A", "F4D", "F5D", "S5A", "S7A"],
        "outputs": outputs,
        "notes": [
            "All cohorts remain separate in UMAP and CellChat views.",
            "Simplex points are patient-level single-cell state proportions.",
            "The S7A-inspired map uses CellChat ligand-receptor compatibility, not ligand-target prediction.",
        ],
    }
    (output_dir / "CANDIDATE_FIGURE_SUMMARY_v1.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    print(json.dumps({
        "status": validation["status"],
        "figures": len(outputs),
        "output_dir": str(output_dir),
    }, indent=2))


if __name__ == "__main__":
    main()
