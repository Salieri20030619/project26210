from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("phase90_v1", SCRIPT_DIR / "01_build_main_figure_v1.py")
base = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(base)


def draw_umap_v2(ax, frame: pd.DataFrame, dataset: str) -> None:
    base.clean_axis(ax)
    for category in base.BROAD_COLORS:
        sub = frame[frame["plot_category"].eq(category)]
        if len(sub):
            ax.scatter(
                sub["UMAP_1"], sub["UMAP_2"], s=0.75,
                c=base.BROAD_COLORS[category], alpha=0.62, linewidths=0, rasterized=True,
            )
    ax.set_aspect("equal", adjustable="datalim")
    status = frame.attrs.get("status", {})
    ax.text(
        0.5, 0.985,
        f"{dataset}  •  {status.get('patients', '?')} patients  •  {status.get('total_cells', len(frame)):,} cells",
        transform=ax.transAxes, ha="center", va="top", fontsize=6.5,
        fontweight="bold", color=base.COLORS["ink"],
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.80, pad=1.2), zorder=12,
    )
    for category in ["Tumor / epithelial", "Myeloid", "T cell", "Endothelial"]:
        sub = frame[frame["plot_category"].eq(category)]
        if len(sub) > 150:
            x, y = sub[["UMAP_1", "UMAP_2"]].median().values
            label = {"Tumor / epithelial": "Tumour", "Myeloid": "Myeloid", "T cell": "T", "Endothelial": "Endo"}[category]
            ax.text(
                x, y, label, fontsize=5.2, ha="center", va="center",
                color=base.COLORS["ink"], fontweight="bold",
                path_effects=[base.pe.withStroke(linewidth=1.6, foreground="white")],
            )


def feature_umap_v2(ax, frame, column, cmap, title):
    base.clean_axis(ax)
    x, y = frame["UMAP_1"].to_numpy(), frame["UMAP_2"].to_numpy()
    values = frame[column].to_numpy(float)
    ax.scatter(x, y, s=0.55, c="#DADDE0", alpha=0.30, linewidths=0, rasterized=True)
    lo, hi = base.robust_limits(pd.Series(values), 0.05, 0.99)
    normed = np.clip((values - lo) / (hi - lo), 0, 1)
    mask = np.isfinite(normed) & (normed > 0.03)
    idx = np.argsort(normed[mask])
    ax.scatter(
        x[mask][idx], y[mask][idx], s=0.8, c=normed[mask][idx], cmap=cmap,
        vmin=0, vmax=1, alpha=0.80, linewidths=0, rasterized=True,
    )
    ax.set_aspect("equal", adjustable="datalim")
    ax.text(
        0.5, 0.985, title, transform=ax.transAxes, ha="center", va="top",
        fontsize=6.4, fontweight="bold", color=base.COLORS["ink"],
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.80, pad=1.0), zorder=12,
    )
    return normed


def program_dynamics_v2(frame: pd.DataFrame, bins=24) -> pd.DataFrame:
    med = base._program_dynamics_v1(frame, bins=bins)
    cols = [
        "MONOCYTE_PROGRAM", "SPP1_MARKER20_ANCHOR_EXCLUDED",
        "C1QC_MARKER20_ANCHOR_EXCLUDED", "IMMUNOSUPPRESSION_NO_IL10",
    ]
    for col in cols:
        med[col] = med[col].rolling(window=5, center=True, min_periods=2).mean()
    return med


def draw_panel_e_v2(fig, spec, atlas_frames, pseudotime, concordance):
    before = set(fig.axes)
    base._draw_panel_e_v1(fig, spec, atlas_frames, pseudotime, concordance)
    for ax in [item for item in fig.axes if item not in before]:
        title = ax.get_title()
        if title.startswith("GSE") and "method rho=" in title:
            ax.set_title("")
            ax.text(
                0.5, 0.985, title, transform=ax.transAxes, ha="center", va="top",
                fontsize=6.0, fontweight="bold", color=base.COLORS["ink"],
                bbox=dict(facecolor="white", edgecolor="none", alpha=0.80, pad=1.0), zorder=12,
            )


def draw_replication_bridge(ax, selected: pd.DataFrame) -> None:
    datasets = ["GSE264586", "GSE294109"]
    y = np.arange(len(selected))
    ax.set_xlim(-1.13, 1.70)
    ax.set_ylim(len(selected) - 0.5, -0.5)
    ax.set_xticks([0.20, 1.20], ["GSE264586\n12 patients", "GSE294109\n13 patients"], fontsize=6.0, fontweight="bold")
    labels = [base.shorten_label(str(x), 30) for x in selected["label"]]
    ax.set_yticks(y, labels, fontsize=5.0)
    ax.tick_params(axis="both", length=0, pad=2)
    for spine in ax.spines.values():
        spine.set_visible(False)

    for i, row in selected.iterrows():
        rhos = [float(row[f"{d}_median_adjusted_rho"]) if np.isfinite(row[f"{d}_median_adjusted_rho"]) else np.nan for d in datasets]
        qs = [float(row[f"{d}_adjusted_q_value"]) if np.isfinite(row[f"{d}_adjusted_q_value"]) else np.nan for d in datasets]
        same_direction = np.all(np.isfinite(rhos)) and np.sign(rhos[0]) == np.sign(rhos[1])
        line_color = base.COLORS["tymp"] if same_direction and np.nanmean(rhos) >= 0 else base.COLORS["monocyte"] if same_direction else "#CBD1D6"
        ax.plot([0.20, 1.20], [i, i], color=line_color, lw=0.65 if not same_direction else 1.05, alpha=0.45, zorder=1)
        for x, rho, q in zip([0.20, 1.20], rhos, qs):
            if not np.isfinite(rho):
                ax.scatter(x, i, s=20, facecolors="white", edgecolors="#B6BDC3", linewidths=0.7, zorder=3)
                ax.text(x, i, "NA", ha="center", va="center", fontsize=3.7, color=base.COLORS["muted"], zorder=4)
                continue
            q_strength = min(max(-np.log10(max(q, 1e-12)), 0), 4) / 4 if np.isfinite(q) else 0
            size = 24 + 62 * q_strength
            color = base.DIVERGING(np.clip((rho + 0.45) / 0.90, 0, 1))
            ax.scatter(x, i, s=size, c=[color], edgecolors="white", linewidths=0.7, zorder=3)
            ax.text(x, i, f"{rho:+.2f}", ha="center", va="center", fontsize=3.8, color="white" if abs(rho) > 0.20 else base.COLORS["ink"], fontweight="bold", zorder=4)

    groups = selected.assign(group=selected["endpoint"].str.replace("_tymp", "", regex=False) + " | " + selected["family"])
    for group, sub in groups.groupby("group", sort=False):
        start, end = sub.index.min(), sub.index.max()
        color = base.COLORS["tymp"] if group.startswith("myeloid") else base.COLORS["c1qc"]
        if "metabolic" in group:
            color = base.COLORS["tumor"] if group.startswith("myeloid") else "#B47A9B"
        ax.add_patch(patches.Rectangle((-0.93, start - 0.5), 0.07, end - start + 1, transform=ax.transData, facecolor=color, edgecolor="none", clip_on=False))
        group_label = group.replace("myeloid_dc", "Myeloid/DC").replace("treg", "Treg").replace(" | ", "\n")
        ax.text(-1.00, (start + end) / 2, group_label, ha="right", va="center", fontsize=4.5, color=base.COLORS["muted"], clip_on=False)
        if end < len(selected) - 1:
            ax.axhline(end + 0.5, color="#E6EAED", lw=0.8)

    norm = base.Normalize(vmin=-0.45, vmax=0.45)
    sm = plt.cm.ScalarMappable(norm=norm, cmap=base.DIVERGING)
    cax = ax.inset_axes([0.43, -0.09, 0.47, 0.025])
    cb = plt.colorbar(sm, cax=cax, orientation="horizontal")
    cb.outline.set_visible(False)
    cb.ax.tick_params(labelsize=4.7, length=2)
    cb.set_label("Median adjusted Spearman rho", fontsize=5.0)
    ax.text(0.90, -0.09, "bubble area = -log10(q)", transform=ax.transAxes, ha="right", va="top", fontsize=4.7, color=base.COLORS["muted"])


def draw_panel_d_v2(fig, spec, retained):
    outer = fig.add_subplot(spec)
    base.clean_axis(outer)
    base.panel_heading(
        outer, "D", "Cross-cohort replication of TYMP-associated immune programs",
        "Patient-level residual correlations • frozen signatures",
    )
    gs = spec.subgridspec(1, 2, width_ratios=[1.52, 0.98], wspace=0.16)
    selected = base.select_phase67_rows(retained)
    ax = fig.add_subplot(gs[0, 0])
    draw_replication_bridge(ax, selected)
    ax_net = fig.add_subplot(gs[0, 1])
    base.draw_program_constellation(ax_net, selected)
    ax_net.text(-1.18, 1.19, f"{sum(selected.endpoint == 'myeloid_dc_tymp')} two-cohort Myeloid/DC signatures", fontsize=5.3, color=base.COLORS["tymp"], fontweight="bold")
    ax_net.text(-1.18, 1.08, f"{sum(selected.endpoint == 'treg_tymp')} single-external Treg signatures", fontsize=5.3, color=base.COLORS["c1qc"], fontweight="bold")


def rehash_outputs(output_dir: Path) -> None:
    summary_path = output_dir / "BUILD_SUMMARY_v1.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    renamed = {}
    for ext in ["png", "pdf", "svg"]:
        old = output_dir / f"PROJECT26210_IMMUNE_MAIN_FIGURE_v1.{ext}"
        new = output_dir / f"PROJECT26210_IMMUNE_MAIN_FIGURE_v2.{ext}"
        old.replace(new)
        renamed[ext] = {"path": str(new), "bytes": new.stat().st_size, "sha256": base.sha256(new)}
    summary["status"] = "GENERATED_PENDING_VISUAL_QA_V2"
    summary["outputs"] = renamed
    (output_dir / "BUILD_SUMMARY_v2.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    summary_path.unlink()


base._program_dynamics_v1 = base.program_dynamics
base._draw_panel_e_v1 = base.draw_panel_e
base.draw_umap = draw_umap_v2
base.feature_umap = feature_umap_v2
base.program_dynamics = program_dynamics_v2
base.draw_panel_e = draw_panel_e_v2
base.draw_panel_d = draw_panel_d_v2


if __name__ == "__main__":
    base.main()
    rehash_outputs(Path(sys.argv[3]))
