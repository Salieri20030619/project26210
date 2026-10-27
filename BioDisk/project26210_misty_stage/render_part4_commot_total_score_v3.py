#!/usr/bin/env python3

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path("/mnt/BioDisk/ludingyu/project26210/results/part4_spatial")
OUTPUT = ROOT / "advanced_patch_commot_figures_v3"
SAMPLES = ("GLMF1", "GLMF2")
AXIS_ORDER = [
    "TGFB1_TGFBR1_TGFBR2", "POSTN_ITGAV_ITGB3", "POSTN_ITGAV_ITGB5", "CCL2_CCR2",
    "CSF1_CSF1R", "SPP1_CD44", "SPP1_ITGAV_ITGB3", "MIF_CD74_CXCR4", "MIF_CD74_CD44",
    "VEGFA_KDR", "APOE_TREM2_TYROBP", "APP_CD74",
]
SCALE_ORDER = ["contact_110um", "secreted_200um", "secreted_500um"]


def render(sample: str) -> None:
    source = ROOT / f"advanced_patch_commot_{sample.lower()}_v2" / f"{sample}_commot_transport_summary_v2.tsv"
    summary = pd.read_csv(source, sep="\t")
    total = summary.pivot(index="axis_id", columns="scale", values="total_score").reindex(AXIS_ORDER).reindex(columns=SCALE_ORDER)
    raw = total.to_numpy(dtype=float)
    values = np.ma.masked_invalid(np.log10(raw))
    cmap = plt.colormaps["viridis"].copy()
    cmap.set_bad("white")
    finite = values.compressed()
    fig, axis = plt.subplots(figsize=(10, 7), constrained_layout=True)
    image = axis.imshow(values, cmap=cmap, aspect="auto", vmin=float(finite.min()), vmax=float(finite.max()))
    axis.set_xticks(range(len(total.columns)), total.columns)
    axis.set_yticks(range(len(total.index)), total.index)
    axis.set_title(f"{sample}: targeted COMMOT total transport score")
    median = float(np.median(finite))
    for row in range(raw.shape[0]):
        for column in range(raw.shape[1]):
            if np.isfinite(raw[row, column]):
                log_value = float(np.log10(raw[row, column]))
                axis.text(column, row, f"{raw[row, column]:.1f}", ha="center", va="center", fontsize=7,
                          color="white" if log_value < median else "black")
            else:
                axis.text(column, row, "not modeled", ha="center", va="center", fontsize=6, color="#777777")
    fig.colorbar(image, ax=axis, fraction=0.035, label="log10(total score)")
    fig.text(0.5, -0.01, "Targeted within-section scores; absolute values are not standardized across axes or sections. White cells were not modeled.",
             ha="center", fontsize=8)
    fig.savefig(OUTPUT / f"{sample}_commot_total_score_v3.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUTPUT / f"{sample}_commot_total_score_v3.pdf", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=False)
    for sample in SAMPLES:
        render(sample)
    audit = {
        "status": "complete",
        "scope": "visualization-only correction; numerical COMMOT v2 outputs unchanged",
        "supersedes": [f"{sample}_commot_total_score_v2.png/pdf" for sample in SAMPLES],
        "reason": "v2 encoded unmodeled scale-axis cells as near-zero dark-purple values",
        "correction": "v3 masks unmodeled cells white and labels them 'not modeled'",
    }
    (OUTPUT / "render_metadata_v3.json").write_text(json.dumps(audit, indent=2) + "\n")
    (OUTPUT / "status.json").write_text(json.dumps({"status": "complete"}, indent=2) + "\n")
    print("COMMOT_TOTAL_SCORE_FIGURE_V3_ALL_PASS")


if __name__ == "__main__":
    main()
