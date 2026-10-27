from __future__ import annotations

import hashlib
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
PHASE = Path(__file__).resolve().parent
OUT = PHASE / "figures_local_v1"
QC = PHASE / "qc_local_v1"
OUT.mkdir(parents=True, exist_ok=True)
QC.mkdir(parents=True, exist_ok=True)

PH32 = ROOT / "phase32_part1_m03_cnv_clone_trees_local_v1" / "results_local_v1"
PH35 = ROOT / "phase35_part1_c01_c04_cell_communication_local_v1" / "server_results_sync_v1"

mpl.rcParams.update(
    {
        "font.family": "Arial",
        "font.size": 8.5,
        "axes.titlesize": 10.5,
        "axes.labelsize": 8.5,
        "xtick.labelsize": 7.5,
        "ytick.labelsize": 7.5,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "axes.linewidth": 0.8,
    }
)

BLUE = "#3E6FA8"
ORANGE = "#D9793D"
RED = "#B84A4A"
GRAY = "#8B9098"
BRANCH_COLORS = [
    "#5E81AC",
    "#D08770",
    "#A3BE8C",
    "#B48EAD",
    "#EBCB8B",
    "#81A1C1",
    "#BF616A",
    "#8FBCBB",
]


def save(fig: plt.Figure, stem: str) -> None:
    fig.savefig(OUT / f"{stem}.png", dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(OUT / f"{stem}.pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def panel(ax: plt.Axes, letter: str) -> None:
    ax.text(
        -0.06,
        1.04,
        letter,
        transform=ax.transAxes,
        fontsize=14,
        fontweight="bold",
        va="top",
    )


def clean(ax: plt.Axes, grid: bool = False) -> None:
    ax.spines[["top", "right"]].set_visible(False)
    if grid:
        ax.grid(axis="y", color="#DADDE1", linewidth=0.6, alpha=0.8)
        ax.set_axisbelow(True)


def short_clone(node_id: str) -> str:
    return node_id.rsplit("_", 1)[-1]


def build_tree_layout(
    patient_edges: pd.DataFrame,
) -> tuple[dict[str, tuple[float, float]], dict[str, list[str]]]:
    children: dict[str, list[str]] = {}
    for row in patient_edges.itertuples(index=False):
        children.setdefault(row.parent_node, []).append(row.child_node)
        children.setdefault(row.child_node, [])
    for node in children:
        children[node] = sorted(children[node])

    leaves: list[str] = []

    def collect_leaves(node: str) -> None:
        kids = children.get(node, [])
        if not kids:
            leaves.append(node)
            return
        for kid in kids:
            collect_leaves(kid)

    collect_leaves("NEUTRAL_ROOT")
    leaf_x = {node: float(i) for i, node in enumerate(leaves)}
    x_cache: dict[str, float] = {}

    def x_pos(node: str) -> float:
        if node in x_cache:
            return x_cache[node]
        kids = children.get(node, [])
        value = leaf_x[node] if not kids else float(np.mean([x_pos(kid) for kid in kids]))
        x_cache[node] = value
        return value

    y_cache = {"NEUTRAL_ROOT": 0.0}

    def add_y(node: str) -> None:
        for kid in children.get(node, []):
            y_cache[kid] = y_cache[node] + 1.0
            add_y(kid)

    add_y("NEUTRAL_ROOT")
    positions = {node: (x_pos(node), y_cache[node]) for node in children}
    return positions, children


def edge_event_labels(transitions: pd.DataFrame) -> dict[tuple[str, str, str], str]:
    changed = transitions.loc[transitions.state_delta_from_parent != 0].copy()
    changed["direction"] = np.where(changed.state_delta_from_parent > 0, "gain", "loss")
    grouped = (
        changed.groupby(
            ["infercnv_unit", "parent_node", "child_node", "chr_number", "direction"],
            as_index=False,
        )
        .size()
        .rename(columns={"size": "n_bins"})
    )
    labels: dict[tuple[str, str, str], str] = {}
    for key, sub in grouped.groupby(["infercnv_unit", "parent_node", "child_node"]):
        top = sub.sort_values(["n_bins", "chr_number"], ascending=[False, True]).head(1)
        labels[key] = "\n".join(
            f"chr{int(row.chr_number)} {row.direction} ({int(row.n_bins)} bins)"
            for row in top.itertuples(index=False)
        )
    return labels


def fig03_clone_trees() -> None:
    edges = pd.read_csv(PH32 / "m03_tree_edges_v1.tsv", sep="\t")
    nodes = pd.read_csv(PH32 / "m03_tree_nodes_cnmf_annotation_v1.tsv", sep="\t")
    transitions = pd.read_csv(PH32 / "m03_tree_edge_bin_transitions_v1.tsv.gz", sep="\t")
    event_lookup = edge_event_labels(transitions)

    patients = [
        ("GSE207493__RCC104", "F1 subtype"),
        ("GSE207493__RCC86", "F2 subtype"),
        ("GSE156632__RCC4", "F3 subtype"),
        ("GSE156632__RCC5", "F6 subtype"),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(10.8, 7.2))

    for ax, (patient, subtype), letter in zip(axes.flat, patients, "ABCD"):
        pe = edges.loc[edges.infercnv_unit == patient].copy()
        pn = nodes.loc[nodes.infercnv_unit == patient].copy().set_index("node_id")
        positions, _ = build_tree_layout(pe)
        node_order = [node for node in positions if node != "NEUTRAL_ROOT"]
        color_lookup = {
            node: BRANCH_COLORS[i % len(BRANCH_COLORS)] for i, node in enumerate(node_order)
        }

        max_y = max(y for _, y in positions.values())
        n_leaves = sum(1 for node in positions if node not in set(pe.parent_node))
        x_span = max(n_leaves - 1, 1)

        for idx, row in enumerate(pe.itertuples(index=False)):
            x1, y1 = positions[row.parent_node]
            x2, y2 = positions[row.child_node]
            color = color_lookup[row.child_node]
            ax.plot(
                [x1, x2],
                [y1, y2],
                color="#6E737A",
                linewidth=6.2,
                solid_capstyle="round",
                zorder=1,
            )
            ax.plot(
                [x1, x2],
                [y1, y2],
                color=color,
                linewidth=4.2,
                solid_capstyle="round",
                zorder=2,
            )
            event = event_lookup.get((patient, row.parent_node, row.child_node), "")
            if event:
                mx, my = x1 + 0.54 * (x2 - x1), y1 + 0.54 * (y2 - y1)
                if x2 > x1 + 0.05:
                    dx = 5
                elif x2 < x1 - 0.05:
                    dx = -5
                else:
                    dx = 7 if idx % 2 == 0 else -7
                ha = "left" if dx > 0 else "right"
                ax.annotate(
                    event,
                    (mx, my),
                    xytext=(dx, 0),
                    textcoords="offset points",
                    ha=ha,
                    va="center",
                    fontsize=6.2,
                    color="#30343A",
                    linespacing=1.12,
                    zorder=5,
                )

        root_x, root_y = positions["NEUTRAL_ROOT"]
        ax.scatter(
            [root_x], [root_y], s=30, color="#6E737A", edgecolor="white", linewidth=0.6, zorder=4
        )
        ax.annotate(
            "Reference",
            (root_x, root_y),
            xytext=(0, 8),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=6.5,
            color="#52575E",
        )

        for node_index, node in enumerate(node_order):
            x, y = positions[node]
            cells = int(pn.loc[node, "clone_cells"])
            color = color_lookup[node]
            size = 42 + 32 * np.sqrt(cells / max(float(pn.clone_cells.max()), 1.0))
            ax.scatter(
                [x],
                [y],
                s=size,
                color=color,
                edgecolor="white",
                linewidth=0.8,
                zorder=4,
            )
            if node not in set(pe.parent_node):
                label_offset = (0, -9)
                label_ha, label_va = "center", "top"
            else:
                label_offset = (7, 0) if node_index % 2 == 0 else (-7, 0)
                label_ha = "left" if label_offset[0] > 0 else "right"
                label_va = "center"
            ax.annotate(
                f"{short_clone(node)}  n={cells:,}",
                (x, y),
                xytext=label_offset,
                textcoords="offset points",
                ha=label_ha,
                va=label_va,
                fontsize=6.5,
                color="#25282D",
            )

        dataset, sample = patient.split("__", 1)
        ax.set_title(f"{subtype}\n{dataset} · {sample}", pad=8)
        ax.set_xlim(-0.72, x_span + 0.72)
        ax.set_ylim(max_y + 0.72, -0.50)
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)
        panel(ax, letter)

    fig.text(
        0.5,
        0.012,
        "Each branch is annotated with its largest chromosome-level CNV change (5-Mb bins).",
        ha="center",
        va="bottom",
        fontsize=7.5,
        color="#555A61",
    )
    fig.subplots_adjust(left=0.045, right=0.985, top=0.94, bottom=0.075, hspace=0.34, wspace=0.22)
    save(fig, "Fig03_CNV_clone_trees_by_cNMF_program_v4")


def feature_umap_light(
    ax: plt.Axes,
    x: np.ndarray,
    y: np.ndarray,
    value: np.ndarray,
    title: str,
) -> mpl.collections.PathCollection:
    finite = np.isfinite(x) & np.isfinite(y) & np.isfinite(value)
    positive = finite & (value > 0)
    ax.scatter(
        x[finite],
        y[finite],
        s=1.3,
        color="#D5D9DE",
        alpha=0.72,
        linewidths=0,
        rasterized=True,
        zorder=1,
    )
    values = value[positive]
    if len(values):
        hi = float(np.nanpercentile(values, 95))
        if not np.isfinite(hi) or hi <= 0:
            hi = float(np.nanmax(values))
        hi = max(hi, 1e-9)
        order = np.argsort(values)
        sc = ax.scatter(
            x[positive][order],
            y[positive][order],
            c=values[order],
            s=14,
            cmap="magma",
            vmin=0,
            vmax=hi,
            edgecolor="white",
            linewidth=0.25,
            rasterized=True,
            zorder=3,
        )
    else:
        sc = ax.scatter([], [], c=[], cmap="magma", vmin=0, vmax=1)

    n_total = int(finite.sum())
    n_positive = int(positive.sum())
    pct = 100.0 * n_positive / max(n_total, 1)
    ax.set_title(f"{title}\n{n_positive:,}/{n_total:,} cells ({pct:.2f}%)", pad=4)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_aspect("equal", adjustable="datalim")
    for spine in ax.spines.values():
        spine.set_visible(False)
    return sc


def fig14_nichenet_receptors() -> None:
    nndir = PH35 / "c03_nichenet_v1"
    c04dir = PH35 / "c02_c04_analysis_v1"
    pa = pd.read_csv(nndir / "c03_nichenet_patient_ligand_activity_v1.tsv", sep="\t")
    wt = pd.read_csv(nndir / "c03_nichenet_direct_target_weights_v1.tsv", sep="\t")
    um = pd.read_csv(c04dir / "c04_receptor_tymp_f1_umap_source_v1.tsv.gz", sep="\t")
    sm = pd.read_csv(c04dir / "c04_receptor_tymp_f1_summary_v1.tsv", sep="\t")

    fig = plt.figure(figsize=(12.2, 9.0))
    gs = fig.add_gridspec(3, 4, height_ratios=[1, 1.08, 0.95], hspace=0.52, wspace=0.32)

    ax = fig.add_subplot(gs[0, :2])
    rng = np.random.default_rng(26210)
    for i, (ligand, color) in enumerate([("MST1", BLUE), ("ICAM5", ORANGE)]):
        values = pa.loc[pa.ligand == ligand, "complete_target_pearson"].dropna().values
        parts = ax.violinplot(values, positions=[i], showmedians=True, widths=0.65)
        for body in parts["bodies"]:
            body.set_facecolor(color)
            body.set_edgecolor(color)
            body.set_alpha(0.28)
        for key in ["cmedians", "cbars", "cmins", "cmaxes"]:
            if key in parts:
                parts[key].set_color(color)
        ax.scatter(
            i + rng.normal(0, 0.055, len(values)),
            values,
            s=22,
            c=color,
            edgecolor="white",
            linewidth=0.5,
            zorder=3,
        )
    ax.axhline(0, color="#555", linewidth=0.8)
    ax.set_xticks([0, 1], ["MST1", "ICAM5"])
    ax.set_ylabel("NicheNet ligand activity (Pearson)")
    ax.set_title("Patient-level ligand activity")
    clean(ax, grid=True)
    panel(ax, "A")

    ax = fig.add_subplot(gs[0, 2:])
    ligands = ["MST1", "ICAM5"]
    targets = list(dict.fromkeys(wt.target.tolist()))
    norm = mpl.colors.Normalize(wt.nichenet_weight.min(), wt.nichenet_weight.max())
    for row in wt.itertuples(index=False):
        ax.scatter(
            targets.index(row.target),
            ligands.index(row.ligand),
            s=20 + 260 * row.target_weight_percentile,
            c=[row.nichenet_weight],
            cmap="viridis",
            norm=norm,
            edgecolor="#444",
            linewidth=0.3,
        )
    ax.set_xticks(range(len(targets)), targets, rotation=45, ha="right")
    ax.set_yticks(range(2), ligands)
    ax.set_title("Prioritized ligand–target links")
    fig.colorbar(
        mpl.cm.ScalarMappable(norm=norm, cmap="viridis"),
        ax=ax,
        fraction=0.035,
        pad=0.03,
        label="NicheNet weight",
    )
    clean(ax)
    panel(ax, "B")

    features = ["MST1R", "ITGAL_ITGB2", "PRTG", "IGDCC4"]
    for i, feature in enumerate(features):
        ax = fig.add_subplot(gs[1, i])
        sc = feature_umap_light(
            ax,
            um.UMAP1.to_numpy(),
            um.UMAP2.to_numpy(),
            um[feature].to_numpy(),
            feature,
        )
        colorbar = fig.colorbar(sc, ax=ax, fraction=0.045, pad=0.02)
        colorbar.ax.tick_params(labelsize=6)
        colorbar.set_label("log1p CP10K", fontsize=7)
        panel(ax, chr(67 + i))

    ax = fig.add_subplot(gs[2, :])
    outcomes = ["TYMP_continuous", "F1_continuous"]
    receptor_features = [
        "MST1R",
        "ITGAL",
        "ITGB2",
        "ITGAL_ITGB2",
        "CD209",
        "PRTG",
        "IGDCC4",
        "ITGAM_ITGB2_sensitivity",
    ]
    matrix = (
        sm.pivot(index="receptor_feature", columns="outcome", values="median_patient_rho")
        .reindex(receptor_features)[outcomes]
    )
    image = ax.imshow(matrix.values, cmap="RdBu_r", vmin=-0.12, vmax=0.12, aspect="auto")
    ax.set_xticks(range(2), ["TYMP", "F1"])
    ax.set_yticks(
        range(len(receptor_features)),
        [value.replace("_sensitivity", " (sens.)") for value in receptor_features],
    )
    for i, receptor in enumerate(receptor_features):
        for j, outcome in enumerate(outcomes):
            row = sm[(sm.receptor_feature == receptor) & (sm.outcome == outcome)].iloc[0]
            label = f"{row.median_patient_rho:+.2f}\nn={int(row.n_eligible_patients)}"
            if row.n_eligible_patients >= 10 and row.BH_FDR_within_outcome8features < 0.05:
                label += "*"
            color = "white" if abs(row.median_patient_rho) > 0.07 else "#202225"
            ax.text(j, i, label, ha="center", va="center", fontsize=7, color=color)
    ax.set_title("Within-patient receptor association (n < 10 shown directionally)")
    fig.colorbar(image, ax=ax, fraction=0.02, pad=0.02, label="Median Spearman rho")
    panel(ax, "G")

    fig.suptitle(
        "NicheNet ligand activity and receptor context for the TYMP/F1 axis",
        y=0.995,
        fontsize=13,
    )
    fig.subplots_adjust(left=0.07, right=0.975, top=0.92, bottom=0.065)
    save(fig, "Fig14_NicheNet_ligand_receptor_evidence_v2")


def write_qc() -> None:
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
        "revised_panels": ["Fig03", "Fig14C", "Fig14D", "Fig14E", "Fig14F"],
        "source_data_unchanged": True,
    }
    (QC / "status.json").write_text(json.dumps(status, indent=2), encoding="utf-8")


if __name__ == "__main__":
    fig03_clone_trees()
    fig14_nichenet_receptors()
    write_qc()
