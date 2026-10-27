from __future__ import annotations

import hashlib
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


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
    }
)

cell_types = [
    "T cell",
    "NK cell",
    "B cell",
    "Plasma cell",
    "Macrophage",
    "Monocyte",
    "Dendritic cell",
    "Mast cell",
    "Epithelial cell",
    "Endothelial",
    "Pericyte/SMC",
    "Fibroblast",
    "Erythrocyte",
]
markers = [
    "CD3D",
    "NKG7",
    "MS4A1",
    "MZB1",
    "C1QA",
    "S100A8",
    "CLEC9A",
    "TPSAB1",
    "EPCAM",
    "PECAM1",
    "RGS5",
    "COL1A1",
    "HBB",
]

data = pd.read_csv(SRC / "broad_celltype_marker_summary_v1.tsv", sep="\t")
data = data.loc[data["cell_type"].isin(cell_types) & data["marker"].isin(markers)].copy()
data["cell_type"] = pd.Categorical(data["cell_type"], categories=cell_types, ordered=True)
data["marker"] = pd.Categorical(data["marker"], categories=markers, ordered=True)

mean_matrix = data.pivot(index="cell_type", columns="marker", values="mean_log1p").reindex(
    index=cell_types, columns=markers
)
fraction_matrix = data.pivot(index="cell_type", columns="marker", values="fraction_expressing").reindex(
    index=cell_types, columns=markers
)
scaled = (mean_matrix - mean_matrix.mean(axis=0)) / mean_matrix.std(axis=0, ddof=0).replace(0, np.nan)
scaled = scaled.clip(-1.5, 2.5).fillna(0)

fig, ax = plt.subplots(figsize=(9.3, 6.4))
xx, yy = np.meshgrid(np.arange(len(markers)), np.arange(len(cell_types)))
scatter = ax.scatter(
    xx.ravel(),
    yy.ravel(),
    s=18 + 175 * fraction_matrix.to_numpy().ravel(),
    c=scaled.to_numpy().ravel(),
    cmap="RdBu_r",
    vmin=-1.5,
    vmax=2.5,
    edgecolors="none",
)
ax.set_xticks(np.arange(len(markers)), markers, rotation=52, ha="right")
ax.set_yticks(np.arange(len(cell_types)), cell_types)
ax.invert_yaxis()
ax.set_xlabel("Representative marker (one per broad cell type)")
ax.set_ylabel("")
ax.set_xlim(-0.7, len(markers) - 0.3)
ax.set_ylim(len(cell_types) - 0.3, -0.7)
ax.set_axisbelow(True)
ax.grid(color="#ECECEC", linewidth=0.55)
ax.spines[["top", "right"]].set_visible(False)

colorbar = fig.colorbar(scatter, ax=ax, fraction=0.035, pad=0.03)
colorbar.set_label("Scaled mean expression")
for fraction, xpos in zip([0.25, 0.50, 0.75, 1.00], [0, 1, 2, 3]):
    ax.scatter([], [], s=18 + 175 * fraction, color="#777777", alpha=0.75, label=f"{int(fraction * 100)}%")
ax.legend(
    title="Expressing cells",
    loc="upper center",
    bbox_to_anchor=(0.5, -0.18),
    ncol=4,
    frameon=False,
    columnspacing=1.8,
)
fig.tight_layout()

outputs = []
for extension, options in [("png", {"dpi": 300}), ("pdf", {})]:
    path = OUT / f"Fig01D_one_marker_per_celltype_dotplot_v2.{extension}"
    fig.savefig(path, bbox_inches="tight", facecolor="white", **options)
    outputs.append(path)
plt.close(fig)

manifest = QC / "marker_dotplot_manifest_v1.tsv"
rows = [
    {
        "file": path.name,
        "bytes": path.stat().st_size,
        "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "n_cell_types": len(cell_types),
        "n_markers": len(markers),
    }
    for path in outputs
]
pd.DataFrame(rows).to_csv(manifest, sep="\t", index=False)
