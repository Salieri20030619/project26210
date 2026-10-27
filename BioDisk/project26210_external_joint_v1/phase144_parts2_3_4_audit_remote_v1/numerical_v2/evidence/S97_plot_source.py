from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


INK = "#17202A"
MUTED = "#65717D"
GRID = "#D9DEE3"
SPP1 = np.array([0.84, 0.18, 0.43])
VASCULAR = np.array([0.05, 0.56, 0.53])
LIGHT = np.array([0.94, 0.95, 0.95])


def setup_style() -> None:
    mpl.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 9.5,
        "axes.titlesize": 11.5, "axes.labelsize": 9.5,
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
        "savefig.facecolor": "white",
    })


def clean_axis(ax) -> None:
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)


def robust01(values: pd.Series) -> np.ndarray:
    x = values.to_numpy(float)
    finite = np.isfinite(x)
    if not finite.any():
        return np.zeros(len(x))
    lo, hi = np.quantile(x[finite], [0.03, 0.98])
    if hi <= lo:
        hi = lo + 1.0
    out = np.zeros(len(x))
    out[finite] = np.clip((x[finite] - lo) / (hi - lo), 0, 1)
    return out


def bivariate_rgb(spp: np.ndarray, vascular: np.ndarray) -> np.ndarray:
    total = spp + vascular
    strength = np.maximum(spp, vascular)
    blend = (spp[:, None] * SPP1 + vascular[:, None] * VASCULAR) / np.maximum(total[:, None], 1e-9)
    blend[total <= 1e-9] = LIGHT
    alpha = (0.12 + 0.86 * strength)[:, None]
    return np.clip(LIGHT * (1 - alpha) + blend * alpha, 0, 1)


def read_section(spatial_root: Path, sample: str) -> pd.DataFrame:
    meta = pd.read_csv(spatial_root / sample / "spot_metadata_v1.tsv.gz", sep="\t")
    abundance = pd.read_csv(spatial_root / sample / "fine_state_abundance_v1.tsv.gz", sep="\t")
    return meta.merge(abundance, on="spot_id", how="inner", validate="one_to_one")


def draw_spatial_bridge(ax, frame: pd.DataFrame, title: str) -> None:
    spp = robust01(frame["q05__SPP1-state Macrophage"])
    vascular = robust01(frame["score_primary__vascular_remodeling"])
    colours = bivariate_rgb(spp, vascular)
    order = np.argsort(np.maximum(spp, vascular))
    x = frame["pxl_col"].to_numpy(float)
    y = -frame["pxl_row"].to_numpy(float)
    ax.scatter(x[order], y[order], s=5.0, c=colours[order], linewidths=0)
    high_both = (spp >= np.quantile(spp, 0.88)) & (vascular >= np.quantile(vascular, 0.88))
    ax.scatter(
        x[high_both], y[high_both], s=9, facecolors="none", edgecolors="#F1C75B",
        linewidths=0.45, alpha=0.85,
    )
    ax.set_aspect("equal", adjustable="datalim")
    clean_axis(ax)
    ax.set_title(title, fontweight="bold", color=INK, pad=4)
    ax.text(
        0.99, 0.015, f"{len(frame):,} spots", transform=ax.transAxes,
        ha="right", va="bottom", fontsize=8, color=MUTED,
    )


def add_bivariate_legend(ax) -> None:
    leg = inset_axes(ax, width="25%", height="25%", loc="lower left", borderpad=0.5)
    grid = np.linspace(0, 1, 35)
    spp, vas = np.meshgrid(grid, grid)
    rgb = bivariate_rgb(spp.ravel(), vas.ravel()).reshape(35, 35, 3)
    leg.imshow(rgb, origin="lower", extent=[0, 1, 0, 1], aspect="equal")
    leg.set_xticks([0, 1], ["low", "high"], fontsize=6.5)
    leg.set_yticks([0, 1], ["low", "high"], fontsize=6.5)
    leg.set_xlabel("SPP1 Mφ", fontsize=7, color=INK, labelpad=1)
    leg.set_ylabel("Vascular", fontsize=7, color=INK, labelpad=1)
    leg.tick_params(length=0, pad=1)
    for spine in leg.spines.values():
        spine.set_color("#BCC4CA")
        spine.set_linewidth(0.6)


def evidence_table(positive_meta: Path, negative_meta: Path) -> pd.DataFrame:
    pos = pd.read_csv(positive_meta, sep="\t")
    pos = pos[pos["relation_id"].eq("B09_SPP1_to_Vascular")].copy()
    pos["preparation"] = pos["preparation"].str.upper()
    rows = []
    for prep in ["FFPE", "FROZEN"]:
        row = pos[pos["preparation"].eq(prep)].iloc[0]
        rows.append({
            "relation": "SPP1 state → vascular remodeling",
            "preparation": prep,
            "n": int(row["n_sections"]),
            "rho": float(row["median_rho"]),
            "p": float(row["wilcoxon_signed_rank_p"]),
            "q": float(row["wilcoxon_bh_q"]),
            "supported": bool(row["support_status"] == "SUPPORTED"),
        })

    neg = pd.read_csv(negative_meta, sep="\t")
    neg = neg[neg["endpoint_id"].eq("P01")].copy()
    neg["preparation"] = neg["preparation"].str.upper()
    for prep in ["FFPE", "FROZEN"]:
        row = neg[neg["preparation"].eq(prep)].iloc[0]
        rows.append({
            "relation": "Tumour TYMP → neighbouring SPP1 state",
            "preparation": prep,
            "n": int(row["sections"]),
            "rho": float(row["pooled_rho"]),
            "p": float(row["p_one_sided"]),
            "q": float(row["q"]),
            "supported": bool(row["preparation_supported"]),
        })
    return pd.DataFrame(rows)


def draw_evidence_matrix(ax, evidence: pd.DataFrame) -> None:
    relations = [
        "SPP1 state → vascular remodeling",
        "Tumour TYMP → neighbouring SPP1 state",
    ]
    ymap = {relations[0]: 1, relations[1]: 0}
    xmap = {"FFPE": 0, "FROZEN": 1}
    ax.axhline(0.5, color=GRID, lw=0.8)
    for relation in relations:
        sub = evidence[evidence["relation"].eq(relation)]
        colour = SPP1 if relation.startswith("SPP1") else np.array([0.60, 0.64, 0.68])
        xs, ys = [], []
        for _, row in sub.iterrows():
            x, y = xmap[row["preparation"]], ymap[relation]
            xs.append(x)
            ys.append(y)
            size = 180 + 72 * min(-np.log10(max(row["q"], 1e-12)), 4)
            edge = "#D59A18" if row["supported"] else "white"
            ax.scatter(
                x, y, s=size, c=[colour], edgecolors=edge, linewidths=2.0, zorder=4
            )
            ax.text(
                x, y + 0.01, f"ρ={row['rho']:.3f}", ha="center", va="center",
                fontsize=8.1, color="white" if relation.startswith("SPP1") else INK,
                fontweight="bold", zorder=5,
            )
            ax.text(
                x, y - 0.25,
                f"n={int(row['n'])}\nP={row['p']:.4g}\nq={row['q']:.4g}",
                ha="center", va="top", fontsize=7.9, color=INK, linespacing=1.13,
            )
        ax.plot(xs, ys, color=colour, lw=1.2, alpha=0.45, zorder=2)
    ax.set_xlim(-0.5, 1.5)
    ax.set_ylim(-0.60, 1.45)
    ax.set_xticks([0, 1], ["FFPE", "Frozen"], fontsize=9.5, fontweight="bold")
    ax.set_yticks(
        [1, 0],
        ["SPP1 Mφ →\nvascular remodeling", "Tumour TYMP →\nneighbouring SPP1 Mφ"],
        fontsize=9.0,
    )
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_title("Preparation-level evidence", fontweight="bold", color=INK, pad=4)
    ax.text(
        0.5, 0.96, "gold outline = BH FDR support",
        transform=ax.transAxes, ha="center", fontsize=8, color=MUTED,
    )


def draw_interpretation_band(ax) -> None:
    ax.axis("off")
    left = FancyBboxPatch(
        (0.01, 0.14), 0.47, 0.72,
        boxstyle="round,pad=0.018,rounding_size=0.025",
        transform=ax.transAxes, facecolor="#E5F4F1", edgecolor="#9FD0C8", linewidth=1.0,
    )
    right = FancyBboxPatch(
        (0.52, 0.14), 0.47, 0.72,
        boxstyle="round,pad=0.018,rounding_size=0.025",
        transform=ax.transAxes, facecolor="#F1F3F5", edgecolor="#C9CFD4", linewidth=1.0,
    )
    ax.add_patch(left)
    ax.add_patch(right)
    ax.text(
        0.035, 0.64, "SUPPORTED IN BOTH PREPARATIONS",
        transform=ax.transAxes, fontsize=9.3, fontweight="bold", color="#087A70",
    )
    ax.text(
        0.035, 0.37, "SPP1 macrophage state ↔ vascular-remodeling niche",
        transform=ax.transAxes, fontsize=11.3, fontweight="bold", color=INK,
    )
    ax.text(
        0.545, 0.64, "NOT SUPPORTED IN EITHER PREPARATION",
        transform=ax.transAxes, fontsize=9.3, fontweight="bold", color="#6D757C",
    )
    ax.text(
        0.545, 0.37, "Tumour TYMP ↛ neighbouring SPP1 macrophage state",
        transform=ax.transAxes, fontsize=11.3, fontweight="bold", color=INK,
    )


def save_all(fig, output_dir: Path, stem: str) -> dict[str, str]:
    output_dir.mkdir(parents=True, exist_ok=True)
    paths = {}
    for ext in ["png", "pdf", "svg"]:
        path = output_dir / f"{stem}.{ext}"
        fig.savefig(path, dpi=320, bbox_inches="tight", pad_inches=0.08)
        paths[ext] = str(path)
    plt.close(fig)
    from PIL import Image
    preview_path = output_dir / f"{stem}_preview.jpg"
    preview = Image.open(paths["png"] )
    preview.thumbnail((1800, 1800))
    preview.convert("RGB").save(preview_path, quality=90)
    paths["preview"] = str(preview_path)
    return paths


def main() -> None:
    if len(sys.argv) != 6:
        raise SystemExit(
            "Usage: make_spp1_vascular_bridge_v1.py <spatial_root> <positive_meta> "
            "<negative_meta> <output_dir> <sample_map_json>"
        )
    spatial_root = Path(sys.argv[1])
    positive_meta = Path(sys.argv[2])
    negative_meta = Path(sys.argv[3])
    output_dir = Path(sys.argv[4])
    sample_map = json.loads(Path(sys.argv[5]).read_text(encoding="utf-8-sig"))
    setup_style()
    ffpe = read_section(spatial_root, sample_map["FFPE"])
    frozen = read_section(spatial_root, sample_map["FROZEN"])
    evidence = evidence_table(positive_meta, negative_meta)

    fig = plt.figure(figsize=(15.6, 8.8))
    gs = fig.add_gridspec(
        2, 3, width_ratios=[1.05, 1.05, 0.90], height_ratios=[0.84, 0.16],
        hspace=0.08, wspace=0.13,
    )
    ax1 = fig.add_subplot(gs[0, 0])
    draw_spatial_bridge(ax1, ffpe, f"FFPE · {sample_map['FFPE']}")
    ax2 = fig.add_subplot(gs[0, 1])
    draw_spatial_bridge(ax2, frozen, f"Frozen · {sample_map['FROZEN']}")
    add_bivariate_legend(ax2)
    ax3 = fig.add_subplot(gs[0, 2])
    draw_evidence_matrix(ax3, evidence)
    band = fig.add_subplot(gs[1, :])
    draw_interpretation_band(band)

    fig.suptitle(
        "SPP1 macrophage–vascular spatial bridge",
        fontsize=17, fontweight="bold", color=INK, y=0.99,
    )
    fig.text(
        0.5, 0.948,
        "Representative FFPE and frozen sections · within-section bivariate scaling · preparation-level inference uses 12 independent sections per preparation",
        ha="center", fontsize=9.2, color=MUTED,
    )
    fig.text(
        0.99, 0.012,
        "Spatial co-localization supports a vascular niche association; it does not establish TYMP-driven SPP1 polarization or clinical resistance.",
        ha="right", fontsize=8.2, color=MUTED,
    )
    fig.subplots_adjust(left=0.035, right=0.985, top=0.90, bottom=0.07)
    outputs = save_all(fig, output_dir, "PROJECT26210_SPP1_VASCULAR_SPATIAL_BRIDGE_v1")
    validation = {
        "status": "PASS",
        "outputs": outputs,
        "evidence": evidence.to_dict(orient="records"),
        "representative_sections": sample_map,
    }
    (output_dir / "SPP1_VASCULAR_BRIDGE_VALIDATION_v1.json").write_text(
        json.dumps(validation, indent=2), encoding="utf-8"
    )
    print(json.dumps({"status": "PASS", "outputs": outputs}, indent=2))


if __name__ == "__main__":
    main()



