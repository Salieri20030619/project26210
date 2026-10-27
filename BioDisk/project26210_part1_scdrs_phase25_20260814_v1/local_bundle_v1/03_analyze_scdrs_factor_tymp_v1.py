#!/usr/bin/env python3
"""Patient-unit analysis of rebuilt scDRS scores and cNMF Factor-TYMP usage."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_text(path: Path, text: str) -> None:
    part = path.with_name(path.name + ".part")
    part.write_text(text, encoding="utf-8", newline="\n")
    os.replace(part, path)


def bh_adjust(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    order = np.argsort(values)
    ranked = values[order]
    adjusted = ranked * len(values) / np.arange(1, len(values) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    result = np.empty_like(adjusted)
    result[order] = np.clip(adjusted, 0, 1)
    return result


def patient_summary(values: np.ndarray, rng: np.random.Generator) -> dict:
    n = len(values)
    result = {
        "n_eligible_patients": int(n),
        "n_positive": int((values > 0).sum()),
        "n_negative": int((values < 0).sum()),
        "n_zero": int((values == 0).sum()),
        "median_rho": float(np.median(values)),
        "reporting_tier": "full_statistics" if n >= 10 else "direction_only",
    }
    if n >= 10:
        boot_idx = rng.integers(0, n, size=(10000, n))
        boot = np.median(values[boot_idx], axis=1)
        signs = rng.choice(np.array([-1.0, 1.0]), size=(9999, n))
        perm = np.median(signs * values[None, :], axis=1)
        observed = abs(float(np.median(values)))
        result.update(
            {
                "bootstrap_ci95_low": float(np.quantile(boot, 0.025)),
                "bootstrap_ci95_high": float(np.quantile(boot, 0.975)),
                "signflip_p_two_sided": float((1 + (np.abs(perm) >= observed).sum()) / 10000),
                "bootstrap_replicates": 10000,
                "signflip_permutations": 9999,
            }
        )
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", required=True, type=Path)
    parser.add_argument("--score", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--expected-cells", type=int, default=11950)
    parser.add_argument("--min-cells-per-patient", type=int, default=20)
    parser.add_argument("--seed", type=int, default=26210)
    args = parser.parse_args()

    if args.output.exists():
        raise FileExistsError(f"refusing existing output: {args.output}")
    part = args.output.with_name(args.output.name + ".part")
    if part.exists():
        raise FileExistsError(f"refusing stale partial output: {part}")
    part.mkdir(parents=True)

    source = ad.read_h5ad(args.h5ad, backed="r")
    required_obs = {"study", "curated_patient_id", "sample_id", "F1"}
    if not required_obs.issubset(source.obs.columns):
        raise ValueError(f"H5AD obs missing {sorted(required_obs - set(source.obs.columns))}")
    if source.n_obs != args.expected_cells or not source.obs_names.is_unique:
        raise ValueError("H5AD cell-count or uniqueness gate failed")
    obs = source.obs[list(required_obs)].copy()
    obs.index = obs.index.astype(str)

    score = pd.read_csv(args.score, sep="\t", index_col=0)
    score.index = score.index.astype(str)
    if score.index.duplicated().any() or set(score.index) != set(obs.index):
        raise ValueError("scDRS/H5AD cell identity mismatch")
    score = score.reindex(obs.index)
    needed_score = {"norm_score", "pval"}
    if not needed_score.issubset(score.columns):
        raise ValueError(f"score columns missing {sorted(needed_score - set(score.columns))}")
    norm_score = pd.to_numeric(score["norm_score"], errors="coerce").to_numpy(float)
    pval = pd.to_numeric(score["pval"], errors="coerce").to_numpy(float)
    if not np.isfinite(norm_score).all():
        raise ValueError("non-finite normalized scDRS scores")
    if not np.isfinite(pval).all() or ((pval < 0) | (pval > 1)).any():
        raise ValueError("invalid scDRS P values")

    cells = obs.copy()
    cells["scdrs_norm_score"] = norm_score
    cells["scdrs_pval"] = pval
    cells["scdrs_qval_bh_descriptive"] = bh_adjust(pval)
    if "mc_pval" in score.columns:
        cells["scdrs_mc_pval"] = pd.to_numeric(score["mc_pval"], errors="coerce").to_numpy(float)
    if "zscore" in score.columns:
        cells["scdrs_zscore"] = pd.to_numeric(score["zscore"], errors="coerce").to_numpy(float)
    if "X_umap" in source.obsm:
        umap = np.asarray(source.obsm["X_umap"])
        cells["UMAP1"] = umap[:, 0]
        cells["UMAP2"] = umap[:, 1]
    source.file.close()

    patient_rows = []
    for (study, patient), frame in cells.groupby(["study", "curated_patient_id"], sort=True):
        f1 = pd.to_numeric(frame["F1"], errors="coerce").to_numpy(float)
        sc = frame["scdrs_norm_score"].to_numpy(float)
        variable = bool(np.ptp(f1) > 0 and np.ptp(sc) > 0)
        eligible = bool(len(frame) >= args.min_cells_per_patient and variable)
        rho = float(spearmanr(f1, sc).statistic) if eligible else np.nan
        patient_rows.append(
            {
                "study": study,
                "curated_patient_id": patient,
                "n_cells": int(len(frame)),
                "f1_variable": bool(np.ptp(f1) > 0),
                "scdrs_variable": bool(np.ptp(sc) > 0),
                "eligible_primary": eligible,
                "spearman_rho_f1_vs_scdrs": rho,
                "median_f1": float(np.median(f1)),
                "median_scdrs_norm_score": float(np.median(sc)),
            }
        )
    patients = pd.DataFrame(patient_rows)
    eligible = patients.loc[patients["eligible_primary"]].copy()
    if eligible.empty:
        raise ValueError("no eligible patients")

    rng = np.random.default_rng(args.seed)
    overall = patient_summary(eligible["spearman_rho_f1_vs_scdrs"].to_numpy(float), rng)
    overall.update(
        {
            "analysis": "within-patient Spearman F1 versus scDRS normalized score",
            "min_cells_per_patient": args.min_cells_per_patient,
            "seed": args.seed,
            "n_cells": int(len(cells)),
            "n_studies": int(cells["study"].nunique()),
        }
    )

    study_rows = []
    for study, frame in eligible.groupby("study", sort=True):
        summary = patient_summary(frame["spearman_rho_f1_vs_scdrs"].to_numpy(float), rng)
        summary["study"] = study
        study_rows.append(summary)
    studies = pd.DataFrame(study_rows)

    pooled_r = float(pearsonr(pd.to_numeric(cells["F1"]).to_numpy(float), norm_score).statistic)
    descriptive = {
        "pooled_cell_pearson_r_descriptive_only": pooled_r,
        "pooled_cell_p_value_calculated": False,
        "n_cells_q_below_0_10_descriptive": int((cells["scdrs_qval_bh_descriptive"] < 0.10).sum()),
        "n_cells_q_below_0_20_descriptive": int((cells["scdrs_qval_bh_descriptive"] < 0.20).sum()),
    }

    cells.reset_index(names="exact_cell_id").to_csv(
        part / "scdrs_factor_tymp_cell_source_v1.tsv.gz",
        sep="\t",
        index=False,
        compression={"method": "gzip", "compresslevel": 6, "mtime": 0},
    )
    patients.to_csv(part / "patient_within_association_v1.tsv", sep="\t", index=False, na_rep="NA")
    studies.to_csv(part / "study_association_summary_v1.tsv", sep="\t", index=False, na_rep="NA")
    atomic_text(part / "overall_patient_association_v1.json", json.dumps(overall, indent=2) + "\n")
    atomic_text(part / "cell_level_descriptive_v1.json", json.dumps(descriptive, indent=2) + "\n")

    manifest_rows = []
    for path in sorted(p for p in part.iterdir() if p.is_file()):
        manifest_rows.append({"file": path.name, "bytes": path.stat().st_size, "sha256": sha256(path)})
    pd.DataFrame(manifest_rows).to_csv(part / "ANALYSIS_FILE_MANIFEST_v1.tsv", sep="\t", index=False)
    status = {
        "status": "PASS_ANALYSIS_COMPLETE",
        "n_cells": int(len(cells)),
        "n_patients_total": int(len(patients)),
        "n_patients_eligible": int(len(eligible)),
        "reporting_tier": overall["reporting_tier"],
        "patient_unit_primary": True,
        "pooled_cell_inference": False,
    }
    atomic_text(part / "status.json", json.dumps(status, indent=2) + "\n")
    os.replace(part, args.output)


if __name__ == "__main__":
    main()
