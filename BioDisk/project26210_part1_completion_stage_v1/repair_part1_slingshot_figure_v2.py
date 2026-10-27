#!/usr/bin/env python3
from __future__ import annotations

import subprocess
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


PROJECT = Path("/mnt/BioDisk/ludingyu/project26210")
ROOT = PROJECT / "results/part1_rebuild_v1/09_metabolism_trajectory/trajectory_v1"


def main() -> None:
    gate = subprocess.run(
        [str(PROJECT / "check_memory_gate.sh"), "part1_slingshot_figure_v2"],
        check=True,
        text=True,
        capture_output=True,
    )
    (ROOT / "qc/memory_gate_slingshot_figure_v2.txt").write_text(
        gate.stdout + gate.stderr, encoding="utf-8"
    )

    umap = pd.read_csv(ROOT / "input/scvi_umap_v1.tsv.gz", sep="\t", index_col=0)
    scores = pd.read_csv(
        ROOT / "source_data/trajectory_r_cell_scores_v1.tsv.gz",
        sep="\t",
        usecols=["exact_cell_id", "Slingshot_pseudotime_v1"],
        dtype={"exact_cell_id": str},
    ).set_index("exact_cell_id")
    scores = scores.loc[umap.index.astype(str)]
    values = pd.to_numeric(scores["Slingshot_pseudotime_v1"], errors="coerce").to_numpy()
    xy = umap.iloc[:, :2].to_numpy(float)
    if xy.shape != (11950, 2) or np.isfinite(values).sum() != 11950:
        raise ValueError("Slingshot v2 figure input invariant failed")

    order = np.argsort(values)
    fig, ax = plt.subplots(figsize=(5.2, 4.2))
    points = ax.scatter(
        xy[order, 0],
        xy[order, 1],
        c=values[order],
        s=1.2,
        linewidths=0,
        cmap="plasma",
        rasterized=True,
    )
    ax.set_title("Slingshot sensitivity pseudotime", fontsize=11)
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.set_aspect("equal")
    ax.spines[["top", "right"]].set_visible(False)
    colorbar = fig.colorbar(points, ax=ax, fraction=0.046, pad=0.04)
    colorbar.set_label("Pseudotime")
    fig.tight_layout()

    png = ROOT / "figures/Slingshot_pseudotime_v2.png"
    pdf = ROOT / "figures/Slingshot_pseudotime_v2.pdf"
    fig.savefig(png, dpi=400, facecolor="white")
    fig.savefig(pdf)
    plt.close(fig)
    if png.stat().st_size < 100_000 or pdf.stat().st_size < 10_000:
        raise ValueError("Slingshot v2 figure output is unexpectedly small")
    print("PART1_SLINGSHOT_FIGURE_V2_PASS")


if __name__ == "__main__":
    main()
