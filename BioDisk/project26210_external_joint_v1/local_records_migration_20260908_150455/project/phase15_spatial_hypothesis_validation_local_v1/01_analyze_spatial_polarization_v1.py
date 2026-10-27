#!/usr/bin/env python3
"""Descriptive spatial validation of single-cell-derived ecology hypotheses.

The two GLMF sections remain the biological units. Spot-level calculations describe
within-section spatial structure and never create spot-as-replicate p-values.
"""

from __future__ import annotations

import json
import math
import os
from datetime import datetime
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import spearmanr


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
SOURCE = ROOT / "phase12_glmf_spatial_analysis_local_v1" / "server_results_v1" / "tables"
OUT = HERE / "spatial_polarization_v1"
SAMPLES = ["GLMF1", "GLMF2"]
K_VALUES = [6, 12]


def atomic_tsv(frame: pd.DataFrame, path: Path) -> None:
    part = path.with_name(path.name + ".part")
    compression = "gzip" if path.suffix == ".gz" else None
    frame.to_csv(part, sep="\t", index=False, encoding="utf-8", na_rep="NA", compression=compression, lineterminator="\n")
    os.replace(part, path)


def atomic_json(payload: dict, path: Path) -> None:
    part = path.with_name(path.name + ".part")
    part.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    os.replace(part, path)


def rho(x: pd.Series | np.ndarray, y: pd.Series | np.ndarray) -> float:
    xa, ya = np.asarray(x, float), np.asarray(y, float)
    ok = np.isfinite(xa) & np.isfinite(ya)
    if ok.sum() < 3 or np.ptp(xa[ok]) == 0 or np.ptp(ya[ok]) == 0:
        return math.nan
    return float(spearmanr(xa[ok], ya[ok]).statistic)


def quantile_block(values: pd.Series, n: int = 4) -> pd.Series:
    return pd.qcut(values.rank(method="first"), n, labels=False, duplicates="drop").astype(int)


def correlation_row(frame: pd.DataFrame, sample: str, relation: str, x: str, y: str, analysis: str, mask: pd.Series, k: int | None = None) -> dict:
    sub = frame.loc[mask, [x, y]]
    return {
        "sample": sample, "relation": relation, "analysis": analysis,
        "x_metric": x, "y_metric": y, "k_neighbors": k,
        "n_spots": len(sub), "spearman_rho_descriptive": rho(sub[x], sub[y]),
    }


def main() -> None:
    if OUT.exists():
        raise FileExistsError(OUT)
    OUT.mkdir(parents=True)
    table_dir = OUT / "tables"; figure_dir = OUT / "figures"
    table_dir.mkdir(); figure_dir.mkdir()

    long = pd.read_csv(SOURCE / "spot_state_abundance_long.tsv.gz", sep="\t")
    spot = pd.read_csv(SOURCE / "spot_summary.tsv.gz", sep="\t")[["sample", "spot_id", "x", "y"]]
    wide = long.pivot(index=["sample", "spot_id"], columns="state", values="q05_abundance").reset_index()
    wide = spot.merge(wide, on=["sample", "spot_id"], how="inner", validate="one_to_one")
    required = ["F1-high Malignant", "F1-low Malignant", "C1QC-state Macrophage", "SPP1-state Macrophage", "Other Macrophage", "T cell"]
    if not all(c in wide for c in required) or len(wide) != 7079:
        raise RuntimeError("Spatial state/spot conservation failure")

    wide["malignant_total"] = wide["F1-high Malignant"] + wide["F1-low Malignant"]
    wide["macrophage_total"] = wide["C1QC-state Macrophage"] + wide["SPP1-state Macrophage"] + wide["Other Macrophage"]
    wide["f1_high_share_within_malignant"] = wide["F1-high Malignant"] / wide["malignant_total"]
    wide["c1qc_share_within_macrophage"] = wide["C1QC-state Macrophage"] / wide["macrophage_total"]
    wide["spp1_share_within_macrophage"] = wide["SPP1-state Macrophage"] / wide["macrophage_total"]
    total = wide[required].sum(axis=1)
    wide["f1_high_fraction_selected_states"] = wide["F1-high Malignant"] / total
    wide["c1qc_fraction_selected_states"] = wide["C1QC-state Macrophage"] / total
    wide["spp1_fraction_selected_states"] = wide["SPP1-state Macrophage"] / total

    relations = [
        ("F1high_to_C1QC", "f1_high_share_within_malignant", "c1qc_share_within_macrophage"),
        ("F1high_to_SPP1", "f1_high_share_within_malignant", "spp1_share_within_macrophage"),
        ("SPP1_to_T_abundance_not_dysfunction", "spp1_share_within_macrophage", "T cell"),
    ]
    correlation_rows = []
    block_rows = []
    contrast_rows = []
    sample_rows = []

    for sample in SAMPLES:
        frame = wide[wide["sample"] == sample].copy().reset_index(drop=True)
        frame["x_block"] = quantile_block(frame["x"])
        frame["y_block"] = quantile_block(frame["y"])
        frame["spatial_block_4x4"] = frame["x_block"].astype(str) + "_" + frame["y_block"].astype(str)
        tumor_q25, tumor_q50 = frame["malignant_total"].quantile([0.25, 0.5])
        macro_q25, macro_q50 = frame["macrophage_total"].quantile([0.25, 0.5])
        masks = {
            "all_spots": pd.Series(True, index=frame.index),
            "both_lineages_above_q25": (frame["malignant_total"] >= tumor_q25) & (frame["macrophage_total"] >= macro_q25),
            "both_lineages_above_median": (frame["malignant_total"] >= tumor_q50) & (frame["macrophage_total"] >= macro_q50),
        }
        sample_rows.append({
            "sample": sample, "n_spots": len(frame),
            "mean_f1_high_share": frame["f1_high_share_within_malignant"].mean(),
            "mean_c1qc_share": frame["c1qc_share_within_macrophage"].mean(),
            "mean_spp1_share": frame["spp1_share_within_macrophage"].mean(),
            "malignant_total_q25": tumor_q25, "macrophage_total_q25": macro_q25,
        })
        coords = frame[["x", "y"]].to_numpy(float)
        for relation, x, y in relations:
            for mask_name, mask in masks.items():
                correlation_rows.append(correlation_row(frame, sample, relation, x, y, f"same_spot_{mask_name}", mask))
                # Within-block and leave-one-block-out stability diagnose large-domain dependence.
                for block, indices in frame.loc[mask].groupby("spatial_block_4x4").groups.items():
                    sub = frame.loc[indices]
                    if len(sub) >= 30:
                        block_rows.append({"sample": sample, "relation": relation, "mask": mask_name,
                                           "block_analysis": "within_block", "block": block,
                                           "n_spots": len(sub), "spearman_rho_descriptive": rho(sub[x], sub[y])})
                for block in sorted(frame.loc[mask, "spatial_block_4x4"].unique()):
                    keep = mask & (frame["spatial_block_4x4"] != block)
                    block_rows.append({"sample": sample, "relation": relation, "mask": mask_name,
                                       "block_analysis": "leave_one_block_out", "block": block,
                                       "n_spots": int(keep.sum()), "spearman_rho_descriptive": rho(frame.loc[keep, x], frame.loc[keep, y])})

                selected = frame.loc[mask]
                lower, upper = selected[x].quantile([0.25, 0.75])
                lo = selected.loc[selected[x] <= lower, y]
                hi = selected.loc[selected[x] >= upper, y]
                contrast_rows.append({
                    "sample": sample, "relation": relation, "mask": mask_name,
                    "x_low_quartile_cut": lower, "x_high_quartile_cut": upper,
                    "n_low": len(lo), "n_high": len(hi), "median_y_low": lo.median(),
                    "median_y_high": hi.median(), "median_y_high_minus_low": hi.median() - lo.median(),
                })

            for k in K_VALUES:
                neighbors = cKDTree(coords).query(coords, k=k + 1)[1][:, 1:]
                lag_y = frame[y].to_numpy(float)[neighbors].mean(axis=1)
                temp = frame.copy(); temp["neighbor_lag_y"] = lag_y
                for mask_name, mask in masks.items():
                    correlation_rows.append(correlation_row(temp, sample, relation, x, "neighbor_lag_y", f"neighbor_lag_{mask_name}", mask, k))

        # Maps of internal polarization, not global composition fractions.
        fig, axes = plt.subplots(1, 3, figsize=(16, 5.2))
        for ax, metric, title in zip(
            axes,
            ["f1_high_share_within_malignant", "c1qc_share_within_macrophage", "spp1_share_within_macrophage"],
            ["F1-high share within malignant", "C1QC share within macrophages", "SPP1 share within macrophages"],
        ):
            sc = ax.scatter(frame["x"], frame["y"], c=frame[metric], s=4, cmap="viridis", vmin=0, vmax=1, linewidths=0)
            ax.invert_yaxis(); ax.set_aspect("equal"); ax.set_title(title); ax.set_xticks([]); ax.set_yticks([])
            fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.03)
        fig.suptitle(f"{sample}: within-lineage spatial polarization")
        part = figure_dir / f"{sample}_within_lineage_polarization.png.part"
        fig.savefig(part, dpi=180, bbox_inches="tight", format="png"); plt.close(fig)
        os.replace(part, figure_dir / f"{sample}_within_lineage_polarization.png")

    correlation = pd.DataFrame(correlation_rows)
    blocks = pd.DataFrame(block_rows)
    block_summary = blocks.groupby(["sample", "relation", "mask", "block_analysis"], as_index=False).agg(
        blocks=("block", "nunique"), median_rho=("spearman_rho_descriptive", "median"),
        min_rho=("spearman_rho_descriptive", "min"), max_rho=("spearman_rho_descriptive", "max"),
        positive_fraction=("spearman_rho_descriptive", lambda x: float((x > 0).mean())),
        negative_fraction=("spearman_rho_descriptive", lambda x: float((x < 0).mean())),
    )
    atomic_tsv(wide, table_dir / "spot_within_lineage_polarization_v1.tsv.gz")
    atomic_tsv(correlation, table_dir / "spatial_polarization_correlations_v1.tsv")
    atomic_tsv(blocks, table_dir / "spatial_block_stability_long_v1.tsv")
    atomic_tsv(block_summary, table_dir / "spatial_block_stability_summary_v1.tsv")
    atomic_tsv(pd.DataFrame(contrast_rows), table_dir / "polarization_quartile_contrasts_v1.tsv")
    atomic_tsv(pd.DataFrame(sample_rows), table_dir / "sample_polarization_summary_v1.tsv")

    f1c = correlation[correlation["relation"] == "F1high_to_C1QC"]
    status = {
        "schema_version": "phase15_spatial_polarization_v1",
        "status": "PASS_DESCRIPTIVE_SPATIAL_VALIDATION",
        "completed_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "samples": SAMPLES, "spots": int(len(wide)), "biological_units": 2,
        "spot_as_replicate_inference": False,
        "within_lineage_metrics": ["F1-high/(F1-high+F1-low)", "C1QC/(C1QC+SPP1+Other macrophage)"],
        "f1_c1qc_descriptive_rows": int(len(f1c)),
        "warning": "Spot correlations and block stability are descriptive. With two sections they cannot confirm a patient-level association or supply an inferential p-value.",
    }
    atomic_json(status, OUT / "SPATIAL_POLARIZATION_STATUS_v1.json")
    print(json.dumps(status, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
