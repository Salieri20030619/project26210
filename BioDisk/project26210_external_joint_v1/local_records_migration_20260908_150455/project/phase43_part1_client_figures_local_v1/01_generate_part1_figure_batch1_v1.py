from __future__ import annotations

import hashlib
import json
from pathlib import Path

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from PIL import Image


ROOT = Path(__file__).resolve().parent
SRC = ROOT / "source_remote_v1"
OUT = ROOT / "figures_v2"
QC = ROOT / "qc_v2"
OUT.mkdir(parents=True, exist_ok=True)
QC.mkdir(parents=True, exist_ok=True)

plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 9,
        "axes.titlesize": 10,
        "axes.labelsize": 9,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "axes.linewidth": 0.8,
    }
)

CELL_COLORS = {
    "Epithelial cell": "#E46C0A",
    "T cell": "#177DB3",
    "Macrophage": "#9B6A5A",
    "Endothelial": "#13A77B",
    "NK cell": "#5CB4E6",
    "Pericyte/SMC": "#CC79A7",
    "Monocyte": "#E9A900",
    "B cell": "#473A96",
    "Plasma cell": "#B75DB0",
    "Mast cell": "#9B3767",
    "Dendritic cell": "#6CBFB1",
    "Fibroblast": "#9C9E36",
    "Erythrocyte": "#D9583B",
    "Epithelial/Myeloid mixed": "#8A8A8A",
}

FACTOR_COLORS = {
    "F1": "#D55E00",
    "F2": "#0072B2",
    "F3": "#009E73",
    "F4": "#CC79A7",
    "F5": "#E69F00",
    "F6": "#6A51A3",
    "F7": "#56B4E9",
    "NA": "#BDBDBD",
}


def panel_letter(ax, letter: str) -> None:
    ax.text(
        -0.04,
        1.03,
        letter,
        transform=ax.transAxes,
        fontsize=13,
        fontweight="bold",
        ha="right",
        va="bottom",
    )


def save_figure(fig: plt.Figure, stem: str, width: float, height: float) -> None:
    fig.set_size_inches(width, height)
    fig.savefig(OUT / f"{stem}.png", dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(OUT / f"{stem}.pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def style_axes(ax) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", color="#D9D9D9", linewidth=0.5, alpha=0.7)
    ax.set_axisbelow(True)


counts = pd.read_csv(SRC / "celltype_counts_v1.tsv", sep="\t")
counts["color"] = counts["cell_type_broad_v1"].map(CELL_COLORS)

# Figure 1B: horizontal cell-count bar chart.
horizontal = counts.sort_values("n_cells", ascending=True).copy()
fig, ax = plt.subplots()
ax.barh(
    horizontal["cell_type_broad_v1"],
    horizontal["n_cells"],
    color=horizontal["color"],
    edgecolor="none",
)
for y, value in enumerate(horizontal["n_cells"]):
    ax.text(value + counts["n_cells"].max() * 0.012, y, f"{value:,}", va="center", fontsize=7.5)
ax.set_xlabel("Number of cells")
ax.set_ylabel("")
ax.set_xlim(0, counts["n_cells"].max() * 1.18)
style_axes(ax)
ax.grid(axis="x", color="#D9D9D9", linewidth=0.5, alpha=0.7)
ax.grid(axis="y", visible=False)
save_figure(fig, "Fig01B_celltype_counts_horizontal_v1", 7.0, 4.8)

# Figure 1C: vertical cell-count bar chart.
vertical = counts.sort_values("n_cells", ascending=False).copy()
fig, ax = plt.subplots()
bars = ax.bar(
    vertical["cell_type_broad_v1"],
    vertical["n_cells"],
    color=vertical["color"],
    edgecolor="none",
)
ax.set_ylabel("Number of cells")
ax.set_xlabel("")
ax.tick_params(axis="x", rotation=58)
for bar, value in zip(bars, vertical["n_cells"]):
    ax.text(
        bar.get_x() + bar.get_width() / 2,
        value + counts["n_cells"].max() * 0.012,
        f"{value / 1000:.1f}k" if value >= 1000 else f"{value}",
        ha="center",
        va="bottom",
        fontsize=6.8,
        rotation=90,
    )
ax.set_ylim(0, counts["n_cells"].max() * 1.19)
style_axes(ax)
save_figure(fig, "Fig01C_celltype_counts_vertical_v1", 8.3, 4.8)

# Figure 1A-C: cell atlas composite in the requested UMAP + two bar orientations.
umap_img = Image.open(SRC / "umap_cell_type_broad_v1.png").convert("RGB")
fig = plt.figure(constrained_layout=True)
gs = fig.add_gridspec(2, 2, width_ratios=[1.45, 1.0], height_ratios=[1.0, 1.0])
ax_a = fig.add_subplot(gs[:, 0])
ax_a.imshow(umap_img)
ax_a.axis("off")
panel_letter(ax_a, "A")

ax_b = fig.add_subplot(gs[0, 1])
ax_b.barh(
    horizontal["cell_type_broad_v1"],
    horizontal["n_cells"],
    color=horizontal["color"],
    edgecolor="none",
)
ax_b.set_xlabel("Cells")
ax_b.set_ylabel("")
ax_b.tick_params(axis="y", labelsize=6.6)
ax_b.set_xlim(0, counts["n_cells"].max() * 1.08)
style_axes(ax_b)
ax_b.grid(axis="x", color="#D9D9D9", linewidth=0.45, alpha=0.7)
ax_b.grid(axis="y", visible=False)
panel_letter(ax_b, "B")

ax_c = fig.add_subplot(gs[1, 1])
ax_c.bar(
    vertical["cell_type_broad_v1"],
    vertical["n_cells"],
    color=vertical["color"],
    edgecolor="none",
)
ax_c.set_ylabel("Cells")
ax_c.set_xlabel("")
ax_c.tick_params(axis="x", rotation=62, labelsize=5.8)
style_axes(ax_c)
panel_letter(ax_c, "C")
save_figure(fig, "Fig01_cell_landscape_composite_v1", 12.0, 7.7)

# Figure 1D: marker validation dot plot supplied by the validated annotation run.
marker_img = Image.open(SRC / "cluster_marker_dotplot_v2.png").convert("RGB")
fig, ax = plt.subplots()
ax.imshow(marker_img)
ax.axis("off")
panel_letter(ax, "D")
save_figure(fig, "Fig01D_canonical_marker_dotplot_v1", 13.8, 6.5)

# Figure 2B: CNV burden relative to the patient-specific inferCNV threshold.
calls = pd.read_csv(SRC / "tumor_epithelial_malignancy_calls_v1.tsv.gz", sep="\t")
calls = calls.loc[
    calls["malignancy_status_v1"].isin(["malignant_consensus", "non_malignant_consensus"])
].copy()
calls["group"] = calls["malignancy_status_v1"].map(
    {
        "non_malignant_consensus": "Non-malignant\nconsensus",
        "malignant_consensus": "Malignant\nconsensus",
    }
)
calls["relative_burden"] = np.log2(
    (calls["infercnv_cnv_gene_burden"] + 1)
    / (calls["infercnv_high_cnv_threshold"] + 1)
)
order = ["Non-malignant\nconsensus", "Malignant\nconsensus"]
group_colors = ["#62A88C", "#D97A4A"]

def draw_cnv_violin(ax: plt.Axes) -> None:
    arrays = [calls.loc[calls["group"] == group, "relative_burden"].dropna().to_numpy() for group in order]
    violin = ax.violinplot(arrays, positions=[0, 1], widths=0.72, showextrema=False)
    for body, color in zip(violin["bodies"], group_colors):
        body.set_facecolor(color)
        body.set_edgecolor(color)
        body.set_alpha(0.42)
    rng = np.random.default_rng(20260816)
    for index, (values, color) in enumerate(zip(arrays, group_colors)):
        plotted = values if len(values) <= 700 else rng.choice(values, 700, replace=False)
        jitter = rng.normal(index, 0.055, len(plotted))
        ax.scatter(jitter, plotted, s=5, color=color, alpha=0.32, linewidths=0)
        median = np.median(values)
        ax.plot([index - 0.22, index + 0.22], [median, median], color="#303030", linewidth=1.6)
        ax.text(index, ax.get_ylim()[1] if ax.get_ylim()[1] else 0, "")
    ax.axhline(0, color="#555555", linestyle="--", linewidth=0.8)
    ax.set_xticks([0, 1], [f"{order[i]}\n(n={len(arrays[i]):,})" for i in range(2)])
    ax.set_ylabel("log2(CNV burden / patient threshold)")
    ax.set_xlabel("")
    style_axes(ax)

fig, ax = plt.subplots()
draw_cnv_violin(ax)
save_figure(fig, "Fig02B_CNV_burden_violin_v1", 4.7, 4.7)

# Figure 2: inferCNV chromosome heatmap + consensus separation + method cross-table.
infercnv_img = Image.open(SRC / "infercnv_GSE156632_RCC2.png").convert("RGB")
cross_img = Image.open(SRC / "method_cross_table_v1.png").convert("RGB")
fig = plt.figure(constrained_layout=True)
gs = fig.add_gridspec(2, 2, height_ratios=[1.55, 1.0], width_ratios=[1.0, 1.22])
ax_a = fig.add_subplot(gs[0, :])
ax_a.imshow(infercnv_img)
ax_a.axis("off")
panel_letter(ax_a, "A")
ax_b = fig.add_subplot(gs[1, 0])
draw_cnv_violin(ax_b)
panel_letter(ax_b, "B")
ax_c = fig.add_subplot(gs[1, 1])
ax_c.imshow(cross_img)
ax_c.axis("off")
panel_letter(ax_c, "C")
save_figure(fig, "Fig02_malignancy_CNV_composite_v1", 11.5, 10.0)


def short_clone(node_id: str) -> str:
    if node_id == "NEUTRAL_ROOT":
        return "Neutral"
    return node_id.rsplit("_", 1)[-1]


edges = pd.read_csv(
    ROOT.parent / "phase32_part1_m03_cnv_clone_trees_local_v1" / "results_local_v1" / "m03_tree_edges_v1.tsv",
    sep="\t",
)
nodes = pd.read_csv(
    ROOT.parent / "phase32_part1_m03_cnv_clone_trees_local_v1" / "results_local_v1" / "m03_tree_nodes_cnmf_annotation_v1.tsv",
    sep="\t",
)


def tree_positions(graph: nx.DiGraph, root: str = "NEUTRAL_ROOT") -> dict[str, tuple[float, float]]:
    depth = nx.single_source_shortest_path_length(graph, root)
    positions: dict[str, tuple[float, float]] = {}
    for layer in sorted(set(depth.values())):
        layer_nodes = sorted(node for node, value in depth.items() if value == layer)
        if len(layer_nodes) == 1:
            y_values = [0.0]
        else:
            y_values = np.linspace((len(layer_nodes) - 1) / 2, -(len(layer_nodes) - 1) / 2, len(layer_nodes))
        for node, y_value in zip(layer_nodes, y_values):
            positions[node] = (float(layer), float(y_value))
    return positions


patients = [
    "GSE156632__RCC5",
    "GSE207493__RCC104",
    "GSE156632__RCC4",
    "GSE207493__RCC86",
]
fig, axes = plt.subplots(2, 2, constrained_layout=True)
for ax, patient, letter in zip(axes.flat, patients, ["A", "B", "C", "D"]):
    patient_edges = edges.loc[edges["infercnv_unit"] == patient].copy()
    patient_nodes = nodes.loc[nodes["infercnv_unit"] == patient].copy()
    graph = nx.DiGraph()
    graph.add_edges_from(patient_edges[["parent_node", "child_node"]].itertuples(index=False, name=None))
    positions = tree_positions(graph)
    node_lookup = patient_nodes.set_index("node_id")

    for row in patient_edges.itertuples(index=False):
        child_factor = str(node_lookup.loc[row.child_node, "dominant_cnmf_factor"])
        if child_factor == "nan":
            child_factor = "NA"
        x1, y1 = positions[row.parent_node]
        x2, y2 = positions[row.child_node]
        ax.plot([x1, x2], [y1, y2], color=FACTOR_COLORS.get(child_factor, "#777777"), linewidth=2.2)

    for node_id, (x, y) in positions.items():
        if node_id == "NEUTRAL_ROOT":
            color = FACTOR_COLORS["NA"]
            size = 170
            label = "Neutral"
        else:
            row = node_lookup.loc[node_id]
            factor = str(row["dominant_cnmf_factor"])
            color = FACTOR_COLORS.get(factor, "#777777")
            cells = int(row["clone_cells"])
            size = 150 + 36 * np.sqrt(cells)
            label = f"{short_clone(node_id)}\n{factor} · n={cells:,}"
        ax.scatter([x], [y], s=size, color=color, edgecolor="white", linewidth=1.1, zorder=3)
        ax.annotate(
            label,
            xy=(x, y),
            xytext=(0, -15),
            textcoords="offset points",
            ha="center",
            va="top",
            fontsize=6.5,
        )

    ax.set_title(patient.replace("__", " · "), pad=8)
    ax.set_xlabel("CNV-tree depth")
    ax.set_yticks([])
    ax.set_xticks(sorted({int(value[0]) for value in positions.values()}))
    y_values = [value[1] for value in positions.values()]
    ax.set_ylim(min(y_values) - 0.65, max(y_values) + 0.65)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.grid(axis="x", color="#E0E0E0", linewidth=0.55)
    ax.set_axisbelow(True)
    panel_letter(ax, letter)

legend_factors = ["F1", "F2", "F3", "F6", "NA"]
legend = [
    Line2D(
        [0],
        [0],
        marker="o",
        color="none",
        markerfacecolor=FACTOR_COLORS[factor],
        markeredgecolor="none",
        markersize=7,
        label="Neutral reference" if factor == "NA" else factor,
    )
    for factor in legend_factors
]
fig.legend(handles=legend, loc="lower center", ncol=5, frameon=False, bbox_to_anchor=(0.5, -0.02))
save_figure(fig, "Fig03_CNV_clone_trees_by_cNMF_program_v1", 10.8, 8.2)

# Reproducibility and output inventory.
records = []
for path in sorted(OUT.iterdir()):
    if path.is_file():
        records.append(
            {
                "file": path.name,
                "bytes": path.stat().st_size,
                "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            }
        )
pd.DataFrame(records).to_csv(QC / "generated_figure_manifest_v1.tsv", sep="\t", index=False)

status = {
    "status": "GENERATED_PENDING_VISUAL_QC",
    "n_files": len(records),
    "source_scope": "Part 1 rebuilt cohorts: GSE156632, GSE159115, GSE178481, GSE207493, GSE242299",
    "requested_format_reference": "分析待做（新）_reference.docx",
}
(QC / "generation_status_v1.json").write_text(
    json.dumps(status, ensure_ascii=False, indent=2), encoding="utf-8"
)
