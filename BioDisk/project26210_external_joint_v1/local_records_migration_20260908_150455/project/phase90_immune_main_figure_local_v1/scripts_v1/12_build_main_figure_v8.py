from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import matplotlib.patches as patches


SCRIPT_DIR = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("phase90_v2_for_v8", SCRIPT_DIR / "02_build_main_figure_v2.py")
v2 = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(v2)
base = v2.base
original_figure_text = base.mpl.figure.Figure.text


def panel_heading_v8(ax, letter: str, title: str, subtitle: str | None = None) -> None:
    ax.text(-0.015, 1.045, letter, transform=ax.transAxes, ha="right", va="bottom", fontsize=13.2, fontweight="bold", color=base.COLORS["ink"])
    ax.text(0.0, 1.045, title, transform=ax.transAxes, ha="left", va="bottom", fontsize=10.2, fontweight="bold", color=base.COLORS["ink"])
    if subtitle and letter in {"A", "B"}:
        ax.text(1.0, 1.010, subtitle, transform=ax.transAxes, ha="right", va="bottom", fontsize=5.8, color=base.COLORS["muted"])


def compact_program_label(row) -> str:
    label = str(row["label"])
    if row["family"] == "regulon":
        return label.split()[0] + " regulon"
    replacements = {
        "Aerobic respiration and respiratory electron transport": "Aerobic respiration + ETC",
        "Cellular response to hypoxia": "Hypoxia response",
        "Detoxification of Reactive Oxygen Species": "ROS detoxification",
        "Regulation of pyruvate dehydrogenase (PDH) complex": "PDH-complex regulation",
        "Pentose phosphate pathway": "Pentose-phosphate pathway",
        "Respiratory electron transport": "Respiratory electron transport",
    }
    return replacements.get(label, label)


def draw_replication_bridge_v8(ax, selected) -> None:
    datasets = ["GSE264586", "GSE294109"]
    x_program = -0.82
    x_cohorts = [0.43, 1.22]
    ax.set_xlim(-1.10, 1.62)
    ax.set_ylim(len(selected) - 0.5, -1.75)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)

    ax.text(x_program, -1.20, "Program", ha="left", va="center", fontsize=6.0, fontweight="bold", color=base.COLORS["ink"])
    ax.text(x_cohorts[0], -1.20, "GSE264586\n12 patients", ha="center", va="center", fontsize=5.8, fontweight="bold", color=base.COLORS["ink"])
    ax.text(x_cohorts[1], -1.20, "GSE294109\n13 patients", ha="center", va="center", fontsize=5.8, fontweight="bold", color=base.COLORS["ink"])

    group_names = {
        "myeloid_dc_tymp | metabolic": ("Myeloid/DC metabolic", base.COLORS["tumor"]),
        "myeloid_dc_tymp | regulon": ("Myeloid/DC regulon", base.COLORS["tymp"]),
        "treg_tymp | metabolic": ("Treg metabolic", "#B47A9B"),
        "treg_tymp | regulon": ("Treg regulon", base.COLORS["c1qc"]),
    }
    groups = selected.assign(group=selected["endpoint"] + " | " + selected["family"])
    for group, sub in groups.groupby("group", sort=False):
        start, end = int(sub.index.min()), int(sub.index.max())
        title, color = group_names[group]
        ax.axhspan(start - 0.5, end + 0.5, color=color, alpha=0.035, zorder=0)
        ax.add_patch(patches.Rectangle((-1.075, start - 0.5), 0.13, end - start + 1, facecolor=color, edgecolor="none", zorder=1))
        ax.text(-1.01, (start + end) / 2, title, ha="center", va="center", rotation=90, fontsize=4.4, color="white", fontweight="bold", zorder=2)
        if end < len(selected) - 1:
            ax.axhline(end + 0.5, color="#DDE2E6", lw=0.8, zorder=1)

    for i, row in selected.iterrows():
        label = compact_program_label(row)
        ax.text(x_program, i, label, ha="left", va="center", fontsize=5.05, color=base.COLORS["ink"], zorder=4)
        rhos = [float(row[f"{d}_median_adjusted_rho"]) if np.isfinite(row[f"{d}_median_adjusted_rho"]) else np.nan for d in datasets]
        qs = [float(row[f"{d}_adjusted_q_value"]) if np.isfinite(row[f"{d}_adjusted_q_value"]) else np.nan for d in datasets]
        same_direction = np.all(np.isfinite(rhos)) and np.sign(rhos[0]) == np.sign(rhos[1])
        connector = base.COLORS["tymp"] if same_direction and np.nanmean(rhos) >= 0 else base.COLORS["monocyte"] if same_direction else "#CBD1D6"
        ax.plot(x_cohorts, [i, i], color=connector, lw=0.7 if not same_direction else 1.05, alpha=0.42, zorder=2)
        for x, rho, q in zip(x_cohorts, rhos, qs):
            if not np.isfinite(rho):
                ax.scatter(x, i, s=22, facecolors="white", edgecolors="#B6BDC3", linewidths=0.7, zorder=3)
                ax.text(x, i, "NA", ha="center", va="center", fontsize=3.8, color=base.COLORS["muted"], zorder=4)
                continue
            q_strength = min(max(-np.log10(max(q, 1e-12)), 0), 4) / 4 if np.isfinite(q) else 0
            size = 25 + 64 * q_strength
            color = base.DIVERGING(np.clip((rho + 0.45) / 0.90, 0, 1))
            ax.scatter(x, i, s=size, c=[color], edgecolors="white", linewidths=0.7, zorder=3)
            ax.text(x, i, f"{rho:+.2f}", ha="center", va="center", fontsize=3.8, color="white" if abs(rho) > 0.20 else base.COLORS["ink"], fontweight="bold", zorder=4)

    ax.text(
        1.58, len(selected) - 0.75,
        "Fill = rho (teal − / magenta +)   •   circle area = −log10(q)   •   coloured connector = same direction",
        ha="right", va="bottom", fontsize=4.5, color=base.COLORS["muted"],
    )


def draw_panel_d_v8(fig, spec, retained) -> None:
    outer = fig.add_subplot(spec)
    base.clean_axis(outer)
    base.panel_heading(outer, "D", "Cross-cohort replication of TYMP-associated immune programs", None)
    gs = spec.subgridspec(1, 2, width_ratios=[1.78, 0.82], wspace=0.12)
    selected = base.select_phase67_rows(retained)
    ax = fig.add_subplot(gs[0, 0])
    draw_replication_bridge_v8(ax, selected)
    ax_net = fig.add_subplot(gs[0, 1])
    base.draw_program_constellation(ax_net, selected)
    ax_net.text(-1.18, 1.19, f"{sum(selected.endpoint == 'myeloid_dc_tymp')} two-cohort Myeloid/DC signatures", fontsize=5.3, color=base.COLORS["tymp"], fontweight="bold")
    ax_net.text(-1.18, 1.08, f"{sum(selected.endpoint == 'treg_tymp')} single-external Treg signatures", fontsize=5.3, color=base.COLORS["c1qc"], fontweight="bold")


def figure_text_v8(self, x, y, s, *args, **kwargs):
    if isinstance(s, str) and s.startswith("External immune evidence after cohort expansion"):
        s = ""
    return original_figure_text(self, x, y, s, *args, **kwargs)


def rehash_outputs_v8(output_dir: Path) -> None:
    summary_path = output_dir / "BUILD_SUMMARY_v1.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    renamed = {}
    for ext in ["png", "pdf", "svg"]:
        old = output_dir / f"PROJECT26210_IMMUNE_MAIN_FIGURE_v1.{ext}"
        new = output_dir / f"PROJECT26210_IMMUNE_MAIN_FIGURE_v8.{ext}"
        old.replace(new)
        renamed[ext] = {"path": str(new), "bytes": new.stat().st_size, "sha256": base.sha256(new)}
    summary["status"] = "GENERATED_PENDING_FINAL_QA_V8"
    summary["outputs"] = renamed
    (output_dir / "BUILD_SUMMARY_v8.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    summary_path.unlink()


base.panel_heading = panel_heading_v8
base.draw_panel_d = draw_panel_d_v8
base.mpl.figure.Figure.text = figure_text_v8


if __name__ == "__main__":
    base.main()
    rehash_outputs_v8(Path(sys.argv[3]))
