from __future__ import annotations

import importlib.util
import json
import math
import sys
from pathlib import Path

import numpy as np
from matplotlib.path import Path as MplPath
from matplotlib.patches import PathPatch, Wedge


SCRIPT_DIR = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("phase90_v11_for_v12", SCRIPT_DIR / "17_build_main_figure_v11.py")
v11 = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(v11)
base = v11.base


def polar_point(angle: float, radius: float) -> tuple[float, float]:
    return radius * math.cos(angle), radius * math.sin(angle)


def ribbon_patch(src_interval, tgt_interval, color, alpha=0.30) -> PathPatch:
    s0, s1 = src_interval
    t0, t1 = tgt_interval
    outer = 0.92
    inner = 0.16
    src_arc = np.linspace(s0, s1, 7)
    tgt_arc = np.linspace(t0, t1, 7)
    vertices: list[tuple[float, float]] = [polar_point(src_arc[0], outer)]
    codes = [MplPath.MOVETO]
    for angle in src_arc[1:]:
        vertices.append(polar_point(float(angle), outer))
        codes.append(MplPath.LINETO)
    vertices.extend([
        polar_point(s1, inner),
        polar_point(t0, inner),
        polar_point(t0, outer),
    ])
    codes.extend([MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4])
    for angle in tgt_arc[1:]:
        vertices.append(polar_point(float(angle), outer))
        codes.append(MplPath.LINETO)
    vertices.extend([
        polar_point(t1, inner),
        polar_point(s0, inner),
        polar_point(s0, outer),
        (0.0, 0.0),
    ])
    codes.extend([MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4, MplPath.CLOSEPOLY])
    return PathPatch(MplPath(vertices, codes), facecolor=color, edgecolor="none", alpha=alpha, zorder=1)


def draw_cellchat_chord_v12(ax, interactions) -> None:
    base.clean_axis(ax)
    nodes = [
        "Malignant", "SPP1-state Macrophage", "C1QC-state Macrophage",
        "Other Macrophage", "T cell", "Endothelial",
    ]
    node_colors = {
        "Malignant": base.COLORS["tumor"],
        "SPP1-state Macrophage": base.COLORS["spp1"],
        "C1QC-state Macrophage": base.COLORS["c1qc"],
        "Other Macrophage": base.COLORS["other_mac"],
        "T cell": base.COLORS["tcell"],
        "Endothelial": base.COLORS["endothelial"],
    }
    labels = {
        "Malignant": "Tumour", "SPP1-state Macrophage": "SPP1 Mφ",
        "C1QC-state Macrophage": "C1QC Mφ", "Other Macrophage": "Other Mφ",
        "T cell": "T cells", "Endothelial": "Endo",
    }
    x = interactions[
        interactions["dataset"].eq("GSE207493")
        & interactions["source"].isin(nodes)
        & interactions["target"].isin(nodes)
        & interactions["pval"].le(0.05)
        & interactions["source"].ne(interactions["target"])
    ].copy()
    unit_edges = x.groupby(["unit_id", "source", "target"], observed=True)["prob"].sum().reset_index()
    edges = unit_edges.groupby(["source", "target"], observed=True)["prob"].median().reset_index()
    edges = edges.sort_values("prob", ascending=False).head(14).copy()
    edges["visual_weight"] = np.sqrt(edges["prob"].clip(lower=0))

    gap = math.radians(7.0)
    sector_span = (2 * math.pi - len(nodes) * gap) / len(nodes)
    starts = {}
    ends = {}
    angle = math.radians(90)
    for node in nodes:
        starts[node] = angle
        ends[node] = angle + sector_span
        angle += sector_span + gap

    endpoint_total = {node: 0.0 for node in nodes}
    for _, row in edges.iterrows():
        weight = float(row["visual_weight"])
        endpoint_total[row["source"]] += weight
        endpoint_total[row["target"]] += weight
    endpoint_cursor = {node: starts[node] + sector_span * 0.04 for node in nodes}
    usable_span = sector_span * 0.92
    allocations: list[tuple[tuple[float, float], tuple[float, float], str, float]] = []
    for _, row in edges.sort_values("visual_weight", ascending=False).iterrows():
        source, target = row["source"], row["target"]
        weight = float(row["visual_weight"])
        s_width = usable_span * weight / max(endpoint_total[source], 1e-12)
        t_width = usable_span * weight / max(endpoint_total[target], 1e-12)
        s_interval = (endpoint_cursor[source], endpoint_cursor[source] + s_width)
        t_interval = (endpoint_cursor[target], endpoint_cursor[target] + t_width)
        endpoint_cursor[source] += s_width
        endpoint_cursor[target] += t_width
        allocations.append((s_interval, t_interval, source, float(row["prob"])))

    for s_interval, t_interval, source, prob in sorted(allocations, key=lambda item: item[3], reverse=True):
        ax.add_patch(ribbon_patch(s_interval, t_interval, node_colors[source], alpha=0.20 + 0.20 * prob / max(edges["prob"].max(), 1e-12)))

    for node in nodes:
        theta1, theta2 = math.degrees(starts[node]), math.degrees(ends[node])
        ax.add_patch(Wedge((0, 0), 1.08, theta1, theta2, width=0.16, facecolor=node_colors[node], edgecolor="white", linewidth=1.0, zorder=5))
        mid = (starts[node] + ends[node]) / 2
        xlab, ylab = polar_point(mid, 1.29)
        ha = "left" if xlab > 0.15 else "right" if xlab < -0.15 else "center"
        va = "bottom" if ylab > 0.15 else "top" if ylab < -0.15 else "center"
        ax.text(xlab, ylab, labels[node], ha=ha, va=va, fontsize=5.2, fontweight="bold", color=base.COLORS["ink"], zorder=7)

    ax.add_patch(Wedge((0, 0), 1.08, 0, 360, width=0.01, facecolor="none", edgecolor="#D1D7DC", linewidth=0.7, zorder=6))
    ax.set_xlim(-1.48, 1.48)
    ax.set_ylim(-1.48, 1.48)
    ax.set_aspect("equal")
    ax.set_title("Absolute CellChat chord\nGSE207493 patient median", fontsize=6.0, fontweight="bold", pad=1)
    ax.text(0, -1.42, "ribbon width ∝ sqrt(median probability)  •  colour = sender", ha="center", va="bottom", fontsize=4.3, color=base.COLORS["muted"])


def rehash_outputs_v12(output_dir: Path) -> None:
    summary_path = output_dir / "BUILD_SUMMARY_v1.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    renamed = {}
    for ext in ["png", "pdf", "svg"]:
        old = output_dir / f"PROJECT26210_IMMUNE_MAIN_FIGURE_v1.{ext}"
        new = output_dir / f"PROJECT26210_IMMUNE_MAIN_FIGURE_v12.{ext}"
        old.replace(new)
        renamed[ext] = {"path": str(new), "bytes": new.stat().st_size, "sha256": base.sha256(new)}
    summary["status"] = "GENERATED_PENDING_FINAL_QA_V12"
    summary["outputs"] = renamed
    (output_dir / "BUILD_SUMMARY_v12.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    summary_path.unlink()


base.draw_cellchat_network = draw_cellchat_chord_v12


if __name__ == "__main__":
    base.main()
    rehash_outputs_v12(Path(sys.argv[3]))
