#!/usr/bin/env python3
"""Score frozen programs in GLMF spatial counts and run descriptive spatial validation.

This script runs on the server. It treats GLMF1/2 as two sections and never converts
spots into independent biological replicates or emits spot-level inferential p-values.
All paths are explicit command-line inputs.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import platform
from datetime import datetime
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy
from scipy.spatial import cKDTree
from scipy.stats import rankdata, spearmanr


SAMPLES = ["GLMF1", "GLMF2"]
ABUNDANCE_KEY = "q05_cell_abundance_w_sf"
ABUNDANCE_PREFIX = "q05cell_abundance_w_sf_"
K_VALUES = [6, 12]
SIGNATURES = [
    "f1_top100_weighted", "c1qc_state", "spp1_state_no_spp1", "spp1_state_full",
    "hypoxia_angiogenesis", "tcell_dysfunction", "cytotoxicity",
]


def now_iso() -> str:
    return datetime.now().astimezone().isoformat(timespec="seconds")


def atomic_tsv(frame: pd.DataFrame, path: Path) -> None:
    part = Path(str(path) + ".part")
    compression = "gzip" if path.suffix == ".gz" else None
    frame.to_csv(part, sep="\t", index=False, encoding="utf-8", na_rep="NA", compression=compression, lineterminator="\n")
    os.replace(part, path)


def atomic_json(payload: dict, path: Path) -> None:
    part = Path(str(path) + ".part")
    part.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    os.replace(part, path)


def clean_state(name: str) -> str:
    return str(name).removeprefix(ABUNDANCE_PREFIX)


def safe_rho(x: np.ndarray, y: np.ndarray) -> float:
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 3 or np.ptp(x[ok]) == 0 or np.ptp(y[ok]) == 0:
        return math.nan
    return float(spearmanr(x[ok], y[ok]).statistic)


def residualize(y: np.ndarray, covariates: list[np.ndarray]) -> np.ndarray:
    design = np.column_stack([np.ones(len(y)), *[rankdata(np.asarray(x, float)) for x in covariates]])
    target = rankdata(np.asarray(y, float))
    coefficient, *_ = np.linalg.lstsq(design, target, rcond=None)
    return target - design @ coefficient


def qblock(values: np.ndarray) -> np.ndarray:
    return pd.qcut(pd.Series(values).rank(method="first"), 4, labels=False).to_numpy(int)


def load_sample(path: Path, signature: pd.DataFrame, sample: str) -> tuple[pd.DataFrame, list[dict]]:
    adata = ad.read_h5ad(path, backed="r")
    try:
        genes = pd.Index(adata.var_names.astype(str))
        gene_lookup = {g: i for i, g in enumerate(genes)}
        wanted = signature[signature["signature_id"].isin(SIGNATURES)]["gene"].drop_duplicates().tolist()
        present = [g for g in wanted if g in gene_lookup]
        indices = [gene_lookup[g] for g in present]
        selected = adata.X[:, indices]
        if hasattr(selected, "to_memory"):
            selected = selected.to_memory()
        if hasattr(selected, "toarray"):
            selected = selected.toarray()
        counts = np.asarray(selected, dtype=float)
        totals = adata.obs["total_counts"].to_numpy(float) if "total_counts" in adata.obs else np.asarray(adata.X.sum(axis=1)).ravel()
        lognorm = np.log1p(counts * (10000.0 / totals[:, None]))
        z = (lognorm - lognorm.mean(axis=0)) / lognorm.std(axis=0, ddof=0)
        z[:, lognorm.std(axis=0, ddof=0) == 0] = 0
        gene_position = {g: i for i, g in enumerate(present)}

        abundance = adata.obsm[ABUNDANCE_KEY]
        abundance_values = np.asarray(abundance, float)
        states = [clean_state(x) for x in getattr(abundance, "columns", [])]
        if len(states) != abundance_values.shape[1]:
            states = [str(x) for x in adata.uns["mod"]["factor_names"]]
        abund = pd.DataFrame(abundance_values, columns=states)
        coordinates = np.asarray(adata.obsm["spatial"], float)[:, :2]
        spot_ids = adata.obs_names.astype(str).tolist()
    finally:
        adata.file.close()

    frame = pd.DataFrame({"sample": sample, "spot_id": spot_ids, "x": coordinates[:, 0], "y": coordinates[:, 1], "total_counts": totals})
    overlap_rows = []
    for sig in SIGNATURES:
        definition = signature[signature["signature_id"] == sig].sort_values("gene_order")
        present_genes = [g for g in definition["gene"] if g in gene_position]
        overlap_rows.append({"sample": sample, "signature_id": sig, "defined_genes": len(definition), "present_genes": len(present_genes), "missing_genes": ";".join(g for g in definition["gene"] if g not in gene_position)})
        pos = [gene_position[g] for g in present_genes]
        if sig == "f1_top100_weighted":
            weights = definition.set_index("gene").loc[present_genes, "raw_weight"].to_numpy(float)
            frame[f"primary__{sig}"] = (lognorm[:, pos] @ weights) / weights.sum()
        else:
            frame[f"primary__{sig}"] = z[:, pos].mean(axis=1)
        frame[f"sensitivity__{sig}"] = lognorm[:, pos].mean(axis=1)

    for state in ["F1-high Malignant", "F1-low Malignant", "C1QC-state Macrophage", "SPP1-state Macrophage", "Other Macrophage", "T cell", "Endothelial"]:
        frame[state] = abund[state].to_numpy(float)
    frame["malignant_total"] = frame["F1-high Malignant"] + frame["F1-low Malignant"]
    frame["macrophage_total"] = frame["C1QC-state Macrophage"] + frame["SPP1-state Macrophage"] + frame["Other Macrophage"]
    frame["f1_high_share"] = frame["F1-high Malignant"] / frame["malignant_total"]
    frame["c1qc_share"] = frame["C1QC-state Macrophage"] / frame["macrophage_total"]
    frame["spp1_share"] = frame["SPP1-state Macrophage"] / frame["macrophage_total"]
    return frame, overlap_rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--glmf1", type=Path, required=True)
    parser.add_argument("--glmf2", type=Path, required=True)
    parser.add_argument("--signature-genes", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    args.output.mkdir(parents=True)
    signature = pd.read_csv(args.signature_genes, sep="\t")
    paths = {"GLMF1": args.glmf1, "GLMF2": args.glmf2}
    frames = []
    overlap_rows = []
    correlation_rows = []
    block_rows = []

    hypotheses = [
        ("F1polar_to_C1QC_expression", "f1_high_share", "c1qc_state", ["macrophage_total", "total_counts"], "both_lineages"),
        ("F1expression_to_C1QCpolar", "f1_top100_weighted", "c1qc_share", ["malignant_total", "total_counts"], "both_lineages"),
        ("F1polar_to_hypoxia_expression", "f1_high_share", "hypoxia_angiogenesis", ["malignant_total", "Endothelial", "total_counts"], "malignant"),
        ("SPP1polar_to_Tdysfunction_expression", "spp1_share", "tcell_dysfunction", ["T cell", "total_counts"], "macrophage_t"),
        ("SPP1polar_to_cytotoxicity_expression", "spp1_share", "cytotoxicity", ["T cell", "total_counts"], "macrophage_t"),
        ("C1QC_mapping_positive_control", "c1qc_share", "c1qc_state", ["macrophage_total", "total_counts"], "macrophage"),
        ("SPP1_mapping_positive_control", "spp1_share", "spp1_state_full", ["macrophage_total", "total_counts"], "macrophage"),
    ]

    for sample in SAMPLES:
        frame, overlaps = load_sample(paths[sample], signature, sample)
        frames.append(frame)
        overlap_rows.extend(overlaps)
        frame["x_block"] = qblock(frame["x"].to_numpy())
        frame["y_block"] = qblock(frame["y"].to_numpy())
        frame["block"] = frame["x_block"].astype(str) + "_" + frame["y_block"].astype(str)
        coords = frame[["x", "y"]].to_numpy(float)
        thresholds = {
            "malignant_q25": frame["malignant_total"].quantile(0.25), "malignant_q50": frame["malignant_total"].quantile(0.50),
            "macrophage_q25": frame["macrophage_total"].quantile(0.25), "macrophage_q50": frame["macrophage_total"].quantile(0.50),
            "t_q25": frame["T cell"].quantile(0.25), "t_q50": frame["T cell"].quantile(0.50),
        }
        for relation, xbase, ysig, covariate_names, scope in hypotheses:
            x = frame[xbase].to_numpy(float) if xbase in frame else frame[f"primary__{xbase}"].to_numpy(float)
            for method in ["primary", "sensitivity"]:
                y = frame[f"{method}__{ysig}"].to_numpy(float) if f"{method}__{ysig}" in frame else frame[ysig].to_numpy(float)
                covariates = [np.log1p(frame[name].to_numpy(float)) for name in covariate_names]
                y_resid = residualize(y, covariates)
                frame_temp = frame.copy()
                frame_temp["x_metric"] = x
                frame_temp["y_score"] = y
                frame_temp["y_residual"] = y_resid
                if scope == "both_lineages":
                    masks = {
                        "q25": (frame.malignant_total >= thresholds["malignant_q25"]) & (frame.macrophage_total >= thresholds["macrophage_q25"]),
                        "q50": (frame.malignant_total >= thresholds["malignant_q50"]) & (frame.macrophage_total >= thresholds["macrophage_q50"]),
                    }
                elif scope == "malignant":
                    masks = {"q25": frame.malignant_total >= thresholds["malignant_q25"], "q50": frame.malignant_total >= thresholds["malignant_q50"]}
                elif scope == "macrophage_t":
                    masks = {
                        "q25": (frame.macrophage_total >= thresholds["macrophage_q25"]) & (frame["T cell"] >= thresholds["t_q25"]),
                        "q50": (frame.macrophage_total >= thresholds["macrophage_q50"]) & (frame["T cell"] >= thresholds["t_q50"]),
                    }
                else:
                    masks = {"q25": frame.macrophage_total >= thresholds["macrophage_q25"], "q50": frame.macrophage_total >= thresholds["macrophage_q50"]}

                for mask_name, mask in masks.items():
                    for outcome, vector in [("raw_score", y), ("lineage_depth_residual", y_resid)]:
                        correlation_rows.append({"sample": sample, "relation": relation, "score_method": method,
                                                 "outcome": outcome, "analysis": f"same_spot_{mask_name}", "k_neighbors": np.nan,
                                                 "n_spots": int(mask.sum()), "spearman_rho_descriptive": safe_rho(x[mask], vector[mask])})
                        for block in sorted(frame.loc[mask, "block"].unique()):
                            keep = mask & (frame["block"] != block)
                            block_rows.append({"sample": sample, "relation": relation, "score_method": method,
                                               "outcome": outcome, "mask": mask_name, "left_out_block": block,
                                               "n_spots": int(keep.sum()), "spearman_rho_descriptive": safe_rho(x[keep], vector[keep])})
                        for k in K_VALUES:
                            neighbors = cKDTree(coords).query(coords, k=k + 1)[1][:, 1:]
                            lag = vector[neighbors].mean(axis=1)
                            correlation_rows.append({"sample": sample, "relation": relation, "score_method": method,
                                                     "outcome": outcome, "analysis": f"neighbor_lag_{mask_name}", "k_neighbors": k,
                                                     "n_spots": int(mask.sum()), "spearman_rho_descriptive": safe_rho(x[mask], lag[mask])})

    scores = pd.concat(frames, ignore_index=True)
    correlations = pd.DataFrame(correlation_rows)
    blocks = pd.DataFrame(block_rows)
    block_summary = blocks.groupby(["sample", "relation", "score_method", "outcome", "mask"], as_index=False).agg(
        blocks=("left_out_block", "nunique"), median_rho=("spearman_rho_descriptive", "median"),
        min_rho=("spearman_rho_descriptive", "min"), max_rho=("spearman_rho_descriptive", "max"),
        positive_fraction=("spearman_rho_descriptive", lambda x: float((x > 0).mean())),
        negative_fraction=("spearman_rho_descriptive", lambda x: float((x < 0).mean())),
    )
    atomic_tsv(scores, args.output / "spatial_program_scores_v1.tsv.gz")
    atomic_tsv(pd.DataFrame(overlap_rows), args.output / "signature_gene_overlap_v1.tsv")
    atomic_tsv(correlations, args.output / "spatial_program_correlations_v1.tsv")
    atomic_tsv(blocks, args.output / "spatial_program_leave_one_block_out_v1.tsv")
    atomic_tsv(block_summary, args.output / "spatial_program_block_stability_summary_v1.tsv")
    status = {
        "schema_version": "phase15_spatial_program_expression_v1", "status": "PASS_DESCRIPTIVE_SPATIAL_EXPRESSION_VALIDATION",
        "completed_at": now_iso(), "samples": SAMPLES, "spots": len(scores), "biological_units": 2,
        "spot_as_replicate_inference": False, "program_score_rows": len(correlations),
        "python": platform.python_version(), "anndata": ad.__version__, "numpy": np.__version__, "pandas": pd.__version__, "scipy": scipy.__version__,
        "warning": "All correlations are within-section descriptive spatial structure. Two sections cannot confirm patient-level associations.",
    }
    atomic_json(status, args.output / "SPATIAL_EXPRESSION_VALIDATION_STATUS_v1.json")
    print(json.dumps(status, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
