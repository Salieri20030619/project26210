from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import matplotlib.patches as patches
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch


SCRIPT_DIR = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("phase90_v8_for_v11", SCRIPT_DIR / "12_build_main_figure_v8.py")
v8 = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(v8)
base = v8.base
original_panel_c = base.draw_panel_c
original_bridge = v8.draw_replication_bridge_v8


def draw_panel_c_v11(fig, spec, atlas_frames, feature_summary, metabolic_summary, tf_summary):
    before = set(fig.axes)
    original_panel_c(fig, spec, atlas_frames, feature_summary, metabolic_summary, tf_summary)
    for ax in [item for item in fig.axes if item not in before]:
        if ax.get_title() == "Lineage feature landscape":
            ax.set_xticks([0, 1, 2], ["F1", "TYMP", "TYMP+ (%)"], fontsize=5.4)
            ax.xaxis.tick_top()
            ax.tick_params(axis="x", labelrotation=0, pad=2, top=True, bottom=False)
            ax.text(
                0.98, -0.06, "colour = scaled level  •  area = cell count",
                transform=ax.transAxes, ha="right", va="top", fontsize=4.4, color=base.COLORS["muted"],
            )


def draw_replication_bridge_v11(ax, selected) -> None:
    original_bridge(ax, selected)
    for text in list(ax.texts):
        if text.get_text().startswith("Fill = rho"):
            text.remove()
    n = len(selected)
    ax.set_ylim(n + 2.0, -1.75)
    footer = FancyBboxPatch(
        (-0.82, n + 0.28), 2.35, 1.05,
        boxstyle="round,pad=0.035,rounding_size=0.08",
        facecolor="#EEF2F4", edgecolor="#CBD3D9", linewidth=0.7, zorder=2,
    )
    ax.add_patch(footer)
    ax.text(
        0.355, n + 0.80,
        "FILL = rho (teal: negative  |  magenta: positive)     AREA = −log10(q)     COLOURED LINE = same direction",
        ha="center", va="center", fontsize=5.15, fontweight="bold", color=base.COLORS["ink"], zorder=3,
    )


def compact_program(row) -> str:
    label = str(row["label"])
    if row["family"] == "regulon":
        return label.split()[0]
    replacements = {
        "Aerobic respiration and respiratory electron transport": "Aerobic respiration + ETC",
        "Cellular response to hypoxia": "Hypoxia response",
        "Detoxification of Reactive Oxygen Species": "ROS detoxification",
        "Regulation of pyruvate dehydrogenase (PDH) complex": "PDH regulation",
        "Pentose phosphate pathway": "Pentose-phosphate",
        "Respiratory electron transport": "Respiratory ETC",
    }
    return replacements.get(label, label)


def draw_dual_rail_core(ax, selected) -> None:
    base.clean_axis(ax)
    ax.set_xlim(-1.65, 1.65)
    ax.set_ylim(-1.08, 1.08)
    core = selected[selected["endpoint"].eq("myeloid_dc_tymp")].copy()
    core["mean_rho"] = core[["GSE264586_median_adjusted_rho", "GSE294109_median_adjusted_rho"]].mean(axis=1)
    metabolic = core[core["family"].eq("metabolic")].sort_values("mean_rho", ascending=False).reset_index(drop=True)
    regulon = core[core["family"].eq("regulon")].sort_values("mean_rho", ascending=False).reset_index(drop=True)

    center = FancyBboxPatch(
        (-0.24, -0.68), 0.48, 1.36,
        boxstyle="round,pad=0.02,rounding_size=0.18",
        facecolor=base.COLORS["tymp"], edgecolor="white", linewidth=1.2, zorder=5,
    )
    ax.add_patch(center)
    ax.add_patch(FancyBboxPatch(
        (-0.32, -0.76), 0.64, 1.52,
        boxstyle="round,pad=0.02,rounding_size=0.22",
        facecolor=base.COLORS["tymp"], edgecolor="none", alpha=0.10, zorder=1,
    ))
    ax.text(0, 0.16, "TYMP", ha="center", va="center", color="white", fontsize=8.5, fontweight="bold", zorder=6)
    ax.text(0, -0.02, "Myeloid/DC", ha="center", va="center", color="white", fontsize=5.0, zorder=6)
    ax.text(0, -0.22, "R2 core", ha="center", va="center", color="white", fontsize=5.0, fontweight="bold", zorder=6)

    def rail(frame, side: str, color: str) -> None:
        ys = np.linspace(0.78, -0.78, len(frame))
        node_x = -0.83 if side == "left" else 0.83
        edge_x = -0.24 if side == "left" else 0.24
        for y, (_, row) in zip(ys, frame.iterrows()):
            value = abs(float(row["mean_rho"]))
            rad = -0.10 if side == "left" else 0.10
            link = FancyArrowPatch(
                (node_x, y), (edge_x, y * 0.72), connectionstyle=f"arc3,rad={rad}",
                arrowstyle="-", linewidth=0.75 + 4.0 * value,
                color=color, alpha=0.30 + min(value, 0.4), zorder=2,
            )
            ax.add_patch(link)
            ax.scatter(node_x, y, s=38 + 210 * value, c=color, edgecolors="white", linewidths=0.7, zorder=4)
            label = compact_program(row)
            if side == "left":
                ax.text(node_x - 0.07, y, label, ha="right", va="center", fontsize=4.55, color=base.COLORS["ink"])
            else:
                ax.text(node_x + 0.07, y, label, ha="left", va="center", fontsize=4.55, color=base.COLORS["ink"])

    rail(metabolic, "left", base.COLORS["tumor"])
    rail(regulon, "right", base.COLORS["c1qc"])
    ax.text(-0.83, 0.99, "Metabolic programs", ha="center", va="center", fontsize=6.0, fontweight="bold", color=base.COLORS["tumor"])
    ax.text(0.83, 0.99, "TF regulons", ha="center", va="center", fontsize=6.0, fontweight="bold", color=base.COLORS["c1qc"])
    ax.text(0, -0.99, f"{len(metabolic)} metabolic + {len(regulon)} regulon signatures replicated in both external cohorts", ha="center", va="center", fontsize=4.8, color=base.COLORS["muted"])
    ax.set_title("Two-cohort Myeloid/DC core", fontsize=7.3, fontweight="bold", pad=2)


def draw_panel_d_v11(fig, spec, retained) -> None:
    outer = fig.add_subplot(spec)
    base.clean_axis(outer)
    base.panel_heading(outer, "D", "Cross-cohort replication of TYMP-associated immune programs", None)
    gs = spec.subgridspec(1, 2, width_ratios=[1.76, 0.84], wspace=0.10)
    selected = base.select_phase67_rows(retained)
    ax = fig.add_subplot(gs[0, 0])
    draw_replication_bridge_v11(ax, selected)
    ax_core = fig.add_subplot(gs[0, 1])
    draw_dual_rail_core(ax_core, selected)


def draw_cellchat_network_v11(ax, interactions) -> None:
    base.clean_axis(ax)
    selected = ["Malignant", "SPP1-state Macrophage", "C1QC-state Macrophage", "Other Macrophage", "T cell", "Endothelial"]
    x = interactions[
        interactions["dataset"].eq("GSE207493")
        & interactions["source"].isin(selected)
        & interactions["target"].isin(selected)
        & interactions["pval"].le(0.05)
    ].copy()
    unit_edges = x.groupby(["unit_id", "source", "target"], observed=True)["prob"].sum().reset_index()
    edge = unit_edges.groupby(["source", "target"], observed=True)["prob"].median().reset_index().sort_values("prob", ascending=False).head(20)
    angles = np.linspace(np.pi / 2, np.pi / 2 + 2 * np.pi, len(selected), endpoint=False)
    positions = {node: np.array([np.cos(a), np.sin(a)]) for node, a in zip(selected, angles)}
    node_colors = {
        "Malignant": base.COLORS["tumor"], "SPP1-state Macrophage": base.COLORS["spp1"],
        "C1QC-state Macrophage": base.COLORS["c1qc"], "Other Macrophage": base.COLORS["other_mac"],
        "T cell": base.COLORS["tcell"], "Endothelial": base.COLORS["endothelial"],
    }
    labels = {
        "Malignant": "Tumour", "SPP1-state Macrophage": "SPP1 Mph",
        "C1QC-state Macrophage": "C1QC Mph", "Other Macrophage": "Other Mph",
        "T cell": "T cell", "Endothelial": "Endo",
    }
    max_prob = max(edge["prob"].max(), 1e-9)
    for _, row in edge.iterrows():
        src, tgt = row["source"], row["target"]
        p0, p1 = positions[src], positions[tgt]
        arrow = FancyArrowPatch(
            p0 * 0.87, p1 * 0.87, connectionstyle=f"arc3,rad={0.18 if src != tgt else 0.35}",
            arrowstyle="-|>", mutation_scale=5.5,
            lw=0.35 + 2.2 * float(row["prob"]) / max_prob,
            color=node_colors[src], alpha=0.25 + 0.45 * float(row["prob"]) / max_prob, zorder=1,
        )
        ax.add_patch(arrow)
    for node, pos in positions.items():
        ax.scatter(*pos, s=145, color=node_colors[node], edgecolors="white", linewidths=1.0, zorder=5)
        offset = pos * 0.18
        ax.text(pos[0] + offset[0], pos[1] + offset[1], labels[node], ha="center", va="center", fontsize=4.8, color=base.COLORS["ink"])
    ax.set_xlim(-1.48, 1.48)
    ax.set_ylim(-1.48, 1.48)
    ax.set_title("Absolute CellChat network\nGSE207493 patient median", fontsize=6.0, fontweight="bold", pad=1)


def draw_endpoint_bubbles_v11(ax, endpoint_summary) -> None:
    primary = endpoint_summary[endpoint_summary["tier"].str.contains("PRIMARY", na=False)].copy()
    datasets = ["GSE207493", "GSE264586", "GSE294109"]
    endpoint_order = sorted(list(dict.fromkeys(primary["endpoint_id"].tolist())))[:6]
    short = {
        "P01_Malignant_to_SPP1_specificity": "Tumour → SPP1",
        "P02_SPP1_to_Malignant_specificity": "SPP1 → tumour",
        "P03_SPP1_to_Tcell_specificity": "SPP1 → T",
        "P04_SPP1_to_Endothelial_specificity": "SPP1 → Endo",
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
            ax.scatter(xi, yi, s=16 + 5.5 * n, c=[base.DIVERGING(np.clip((rho + 0.7) / 1.4, 0, 1))], edgecolors=base.COLORS["ink"] if bool(row["formal_support"]) else "white", linewidths=0.8)
            if n < 10:
                ax.text(xi, yi, f"n={int(n)}", ha="center", va="center", fontsize=3.7, color=base.COLORS["ink"])
    ax.set_xlim(-0.5, len(datasets) - 0.5)
    ax.set_ylim(len(endpoint_order) - 0.5, -0.5)
    ax.set_xticks(range(3), ["GSE207493", "GSE264586", "GSE294109"], fontsize=4.6)
    ax.xaxis.tick_top()
    ax.tick_params(axis="x", top=True, bottom=False, labelrotation=0, pad=2)
    ax.set_yticks(range(len(endpoint_order)), [short.get(x, x) for x in endpoint_order], fontsize=4.7)
    ax.grid(color=base.COLORS["grid"], lw=0.4, alpha=0.7)
    ax.tick_params(axis="y", length=0, pad=2)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_title("SPP1-versus-C1QC/Other specificity", fontsize=6.0, fontweight="bold", pad=12)
    ax.text(0.98, -0.12, "area = patients  •  outline = formal support", transform=ax.transAxes, ha="right", fontsize=4.4, color=base.COLORS["muted"])


def rehash_outputs_v11(output_dir: Path) -> None:
    summary_path = output_dir / "BUILD_SUMMARY_v1.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    renamed = {}
    for ext in ["png", "pdf", "svg"]:
        old = output_dir / f"PROJECT26210_IMMUNE_MAIN_FIGURE_v1.{ext}"
        new = output_dir / f"PROJECT26210_IMMUNE_MAIN_FIGURE_v11.{ext}"
        old.replace(new)
        renamed[ext] = {"path": str(new), "bytes": new.stat().st_size, "sha256": base.sha256(new)}
    summary["status"] = "GENERATED_PENDING_FINAL_QA_V11"
    summary["outputs"] = renamed
    (output_dir / "BUILD_SUMMARY_v11.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    summary_path.unlink()


base.draw_panel_c = draw_panel_c_v11
base.draw_panel_d = draw_panel_d_v11
base.draw_cellchat_network = draw_cellchat_network_v11
base.draw_endpoint_bubbles = draw_endpoint_bubbles_v11


if __name__ == "__main__":
    base.main()
    rehash_outputs_v11(Path(sys.argv[3]))
