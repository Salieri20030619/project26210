#!/usr/bin/env python3
from __future__ import annotations

import json
import subprocess
from datetime import datetime
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


ROOT = Path("/mnt/BioDisk/ludingyu/project26210")
MEMORY_GATE = ROOT / "check_memory_gate.sh"
STAGE_V = Path("/mnt/BioDisk/wangmingyuan/project26210_tymp_vascular_niche_stage_v1")
STAGE_T = Path("/mnt/BioDisk/wangmingyuan/project26210_tymp_overnight_stage_v1")
VT = STAGE_V / "output_v1/tables"
TT = STAGE_T / "trajectory_v1/tables"
FIG = STAGE_T / "joint_figures_v1"
REFERENCE = ROOT / "results/part4_spatial/ref/ref_signatures.h5ad"

RELABEL = {
    "SPP1_TAM_to_vascular_compartment": "SPP1 TAM → vascular compartment",
    "SPP1_vascular_niche_to_T_cell_balance": "SPP1–vascular niche → T-cell balance",
    "SPP1_vascular_niche_to_cytotoxicity": "SPP1–vascular niche → cytotoxicity",
    "SPP1_vascular_niche_to_hypoxia_angiogenesis": "SPP1–vascular niche → hypoxia/angiogenesis",
    "SPP1_vascular_niche_to_invasion_EMT": "SPP1–vascular niche → invasion/EMT",
    "TYMP_MTA_core_to_SPP1_monocyte_balance": "TYMP–MTA core → SPP1/monocyte balance",
    "TYMP_MTA_core_to_SPP1_vascular_niche": "TYMP–MTA core → SPP1–vascular niche",
}


def save(fig: plt.Figure, stem: str) -> None:
    fig.savefig(FIG / f"{stem}.png", dpi=240, bbox_inches="tight", facecolor="white")
    fig.savefig(FIG / f"{stem}.pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main() -> None:
    if FIG.exists():
        raise SystemExit(f"Refusing to overwrite: {FIG}")
    FIG.mkdir(parents=True)
    gate = subprocess.run([str(MEMORY_GATE)], text=True, capture_output=True, check=False)
    gate_line = (gate.stdout or gate.stderr).strip()
    print(f"[memory_gate] {gate_line}", flush=True)
    if gate.returncode != 0:
        raise RuntimeError(gate_line)
    sns.set_theme(style="whitegrid", context="talk")

    block = pd.read_csv(VT / "sample_spatial_block_tests_v1.tsv", sep="\t")
    tj = block.loc[(block["cohort"] == "TJ-RCC") & (block["grid"] == "5x5")].copy()
    tj["label"] = tj["relation"].map(RELABEL)
    order = [RELABEL[key] for key in RELABEL]
    palette = {"route_independent": "#176D9C", "reference_light": "#C44E52", "reference_dependent": "#8172B2"}
    fig, ax = plt.subplots(figsize=(13.2, 9.0))
    sns.stripplot(data=tj, x="rho", y="label", order=order, hue="route", palette=palette,
                  size=7, alpha=0.76, jitter=0.20, ax=ax)
    med = tj.groupby("label")["rho"].median().reindex(order)
    for y, value in enumerate(med):
        ax.plot([value, value], [y - 0.28, y + 0.28], color="black", lw=3.2, zorder=8)
    ax.axvline(0, color="#333333", ls="--", lw=1.3)
    ax.set(xlabel="Patient-level block Spearman rho (5×5)", ylabel="",
           title="TJ-RCC vascular-niche relations across 12 patients")
    ax.legend(title="Route", loc="lower right", frameon=True, fontsize=10, title_fontsize=11)
    fig.text(0.01, 0.01, "Each dot is one patient; black bar is the median. Association does not establish causality.", fontsize=10)
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    save(fig, "vascular_niche_patient_effects_v1")

    cohort = pd.read_csv(VT / "tjrcc_cohort_direction_summary_v1.tsv", sep="\t")
    key = cohort.loc[cohort["grid"].eq("5x5")].copy()
    key["label"] = key["relation"].map(RELABEL)
    key["support"] = key["direction_supported_patients"].astype(str) + "/" + key["n_patients"].astype(str)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 7.5), gridspec_kw={"width_ratios": [1.2, 1]})
    colors = [palette.get(route, "#777777") for route in key["route"]]
    y = np.arange(len(key))
    ax1.barh(y, key["median_rho"], color=colors, alpha=0.88)
    ax1.axvline(0, color="#333333", ls="--", lw=1.2)
    ax1.set_yticks(y, key["label"])
    ax1.invert_yaxis()
    ax1.set_xlabel("Median patient rho (5×5)")
    ax1.set_title("Cohort direction")
    cells = key[["direction_supported_patients", "n_patients"]].to_numpy(float)
    sns.heatmap(cells, annot=np.array([[f"{int(a)}/{int(b)}"] for a, b in cells]), fmt="", cmap="Blues",
                vmin=0, vmax=12, cbar_kws={"label": "Direction-supported patients"},
                yticklabels=False, xticklabels=["Supported", "Total"], ax=ax2)
    ax2.set_title("Expected-direction counts")
    fig.suptitle("Primary vascular-niche mechanism audit", y=1.01)
    fig.text(0.01, 0.01, "Only the SPP1–vascular niche vs T-cell balance relation reaches global FDR < 0.05 (12/12; FDR=0.00220).", fontsize=10)
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    save(fig, "vascular_niche_cohort_summary_v1")

    # Corrected trajectory overview using stored UMAP coordinates and staged cell scores.
    cells = pd.read_csv(TT / "myeloid_trajectory_cell_scores_v1.tsv.gz", sep="\t", compression="gzip")
    backed = ad.read_h5ad(REFERENCE, backed="r")
    indexer = pd.Index(backed.obs_names).get_indexer(cells["cell_id"])
    if (indexer < 0).any():
        raise ValueError("Trajectory cells missing from reference")
    coords = np.asarray(backed.obsm["X_umap"])[indexer]
    backed.file.close()
    rng = np.random.default_rng(20260804)
    take = np.sort(rng.choice(len(cells), size=min(30000, len(cells)), replace=False))
    labels = cells["trajectory_celltype"].to_numpy()
    colors = {"Monocyte": "#4C78A8", "Other Macrophage": "#F58518",
              "APOE+C1QC+ TAM": "#54A24B", "SPP1+ TAM": "#E45756"}
    panels = [
        ("celltype", None, "Myeloid state"),
        ("dpt_pseudotime", "viridis", "Cross-sectional DPT"),
        ("SPP1_state_no_SPP1", "magma", "SPP1-state score (SPP1 excluded)"),
        ("C1QC_state", "magma", "C1QC-state score"),
        ("SPP1_vs_C1QC", "coolwarm", "SPP1-state minus C1QC-state"),
        ("Vascular_remodeling", "magma", "Vascular-remodeling score"),
    ]
    fig, axes = plt.subplots(2, 3, figsize=(16.5, 10.4))
    for ax, (field, cmap, title) in zip(axes.flat, panels):
        if field == "celltype":
            for label in colors:
                mask = labels[take] == label
                ax.scatter(coords[take][mask, 0], coords[take][mask, 1], s=2.2, alpha=0.62,
                           color=colors[label], rasterized=True, label=label)
            ax.legend(loc="upper left", bbox_to_anchor=(0, -0.02), ncol=2, frameon=False,
                      fontsize=9, markerscale=4)
        else:
            values = cells[field].to_numpy(float)[take]
            lo, hi = np.quantile(values, [0.01, 0.99])
            image = ax.scatter(coords[take, 0], coords[take, 1], c=values, s=2.1, alpha=0.72,
                               cmap=cmap, vmin=lo, vmax=hi, rasterized=True)
            fig.colorbar(image, ax=ax, fraction=0.046, pad=0.025)
        ax.set_title(title, fontsize=14)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlabel("UMAP1", fontsize=10)
        ax.set_ylabel("UMAP2", fontsize=10)
    fig.suptitle("Myeloid cross-sectional state ordering (53,794 cells; 30,000 displayed)", fontsize=18, y=0.995)
    fig.text(0.01, 0.01, "DPT uses the scVI neighbor graph. It is not RNA velocity or lineage tracing and does not prove monocyte-to-SPP1+ TAM conversion.", fontsize=10)
    fig.tight_layout(rect=(0, 0.04, 1, 0.98))
    save(fig, "myeloid_cross_sectional_state_overview_v2")

    assoc = pd.read_csv(TT / "sample_pseudotime_associations_v1.tsv", sep="\t")
    summary = pd.read_csv(TT / "cohort_pseudotime_direction_summary_v1.tsv", sep="\t")
    label_map = {
        "Monocyte_program": "Monocyte program", "SPP1_state_no_SPP1": "SPP1 state (SPP1 excluded)",
        "C1QC_state": "C1QC state", "Vascular_remodeling": "Vascular remodeling",
        "Immunosuppression": "Immunosuppression", "MTA_context": "MTA context",
        "SPP1_vs_C1QC": "SPP1-state minus C1QC-state",
    }
    assoc["label"] = assoc["score"].map(label_map)
    order_t = [label_map[key] for key in label_map]
    fig, ax = plt.subplots(figsize=(12.8, 8.2))
    sns.stripplot(data=assoc, x="rho", y="label", order=order_t, color="#4C78A8",
                  size=5.5, alpha=0.58, jitter=0.22, ax=ax)
    medians = assoc.groupby("label")["rho"].median().reindex(order_t)
    for y, value in enumerate(medians):
        ax.plot([value, value], [y - 0.30, y + 0.30], color="black", lw=3.2)
    ax.axvline(0, color="#333333", ls="--", lw=1.2)
    ax.set(xlabel="Within-sample Spearman rho with DPT", ylabel="",
           title="State-score directions across 63 eligible samples")
    for y, label in enumerate(order_t):
        row = summary.loc[summary["score"].map(label_map).eq(label)].iloc[0]
        ax.text(1.01, y, f"{int(row.positive_samples)}/63 +; FDR={row.sign_test_fdr_bh:.3g}",
                transform=ax.get_yaxis_transform(), va="center", fontsize=9)
    fig.text(0.01, 0.01, "Each dot is one sample; black bar is the median. DPT direction depends on the selected monocyte root.", fontsize=10)
    fig.tight_layout(rect=(0, 0.04, 0.86, 1))
    save(fig, "myeloid_state_patient_directions_v1")

    metadata = {
        "status": "completed",
        "completed_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "memory_gate": gate_line,
        "files": sorted(path.name for path in FIG.iterdir()),
        "inference_boundary": "All panels show spatial association or cross-sectional state ordering; none independently proves lineage conversion, 5-prime-MTA mediation, treatment resistance, or causality.",
    }
    (FIG / "figure_generation_metadata_v1.json").write_text(
        json.dumps(metadata, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(metadata, ensure_ascii=False, indent=2), flush=True)


if __name__ == "__main__":
    main()
