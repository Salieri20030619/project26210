#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse, stats
from scipy.spatial import cKDTree

SEED = 26210
N_PERMUTATIONS = 1000
MIN_ELIGIBLE = 100
MACRO_STATES = [
    "SPP1-state Macrophage",
    "C1QC-state Macrophage",
    "Other Macrophage",
]


def now() -> str:
    return datetime.now().astimezone().isoformat(timespec="seconds")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def bh(values: pd.Series) -> pd.Series:
    p = values.to_numpy(float)
    result = np.full(len(p), np.nan)
    keep = np.isfinite(p)
    observed = p[keep]
    order = np.argsort(observed)
    ranked = observed[order] * len(observed) / np.arange(1, len(observed) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adjusted = np.empty_like(observed)
    adjusted[order] = np.minimum(ranked, 1.0)
    result[np.flatnonzero(keep)] = adjusted
    return pd.Series(result, index=values.index)


def rank_residual(values: np.ndarray, covariates: list[np.ndarray]) -> np.ndarray:
    y = stats.rankdata(np.asarray(values, dtype=float), method="average")
    columns = [np.ones(len(y))]
    for covariate in covariates:
        columns.append(stats.rankdata(np.asarray(covariate, dtype=float), method="average"))
    x = np.column_stack(columns)
    beta, _, _, _ = np.linalg.lstsq(x, y, rcond=None)
    return y - x @ beta


def correlation_and_permutation(
    exposure: np.ndarray,
    outcome: np.ndarray,
    rng: np.random.Generator,
) -> tuple[float, float]:
    x = stats.rankdata(exposure, method="average")
    y = stats.rankdata(outcome, method="average")
    x = (x - x.mean()) / x.std(ddof=1)
    y = (y - y.mean()) / y.std(ddof=1)
    observed = float(np.dot(x, y) / (len(x) - 1))
    exceed = 0
    for _ in range(N_PERMUTATIONS):
        permuted = rng.permutation(x)
        rho = float(np.dot(permuted, y) / (len(x) - 1))
        exceed += rho >= observed
    return observed, (exceed + 1.0) / (N_PERMUTATIONS + 1.0)


def six_nearest(coords: np.ndarray) -> np.ndarray:
    k = min(7, len(coords))
    _, indices = cKDTree(coords).query(coords, k=k)
    if indices.ndim == 1:
        indices = indices[:, None]
    return indices[:, 1:]


def neighbour_mean(values: np.ndarray, neighbours: np.ndarray) -> np.ndarray:
    return np.asarray(values, dtype=float)[neighbours].mean(axis=1)


def extract_gene_counts(adata: ad.AnnData, gene: str) -> np.ndarray:
    genes = pd.Index(adata.var_names.astype(str).str.upper())
    if gene not in genes:
        raise RuntimeError(f"Missing gene {gene}")
    position = int(genes.get_loc(gene))
    values = adata.X[:, position]
    if sparse.issparse(values):
        values = values.toarray()
    return np.asarray(values).reshape(-1).astype(float)


def library_size(adata: ad.AnnData) -> np.ndarray:
    values = np.asarray(adata.X.sum(axis=1)).reshape(-1).astype(float)
    if not np.isfinite(values).all() or (values <= 0).any():
        raise RuntimeError("Invalid section library sizes")
    return values


def abundance_columns(frame: pd.DataFrame, prefix: str) -> dict[str, np.ndarray]:
    result = {}
    for state in [
        "Malignant", "Epithelial-Normal", "Monocyte",
        "SPP1-state Macrophage", "C1QC-state Macrophage", "Other Macrophage",
    ]:
        column = f"{prefix}__{state}"
        if column not in frame:
            raise RuntimeError(f"Missing abundance column {column}")
        result[state] = frame[column].to_numpy(float)
    return result


def analyse_endpoint(
    exposure: np.ndarray,
    exposure_covariates: list[np.ndarray],
    raw_outcome: np.ndarray,
    outcome_covariates: list[np.ndarray],
    eligible: np.ndarray,
    rng: np.random.Generator,
) -> tuple[float, float, int]:
    x = rank_residual(exposure, exposure_covariates)
    y = rank_residual(raw_outcome, outcome_covariates)
    x = x[eligible]
    y = y[eligible]
    ok = np.isfinite(x) & np.isfinite(y)
    x = x[ok]
    y = y[ok]
    if len(x) < MIN_ELIGIBLE or np.std(x) == 0 or np.std(y) == 0:
        return np.nan, np.nan, int(len(x))
    rho, permutation_p = correlation_and_permutation(x, y, rng)
    return rho, permutation_p, int(len(x))


def analyse_section(
    sample: str,
    preparation: str,
    section_dir: Path,
    sample_index: int,
) -> list[dict]:
    posterior_path = section_dir / "fine_state_posterior_v1.h5ad"
    abundance_path = section_dir / "fine_state_abundance_v1.tsv.gz"
    adata = ad.read_h5ad(posterior_path)
    abundance = pd.read_csv(abundance_path, sep="\t")
    if len(abundance) != adata.n_obs:
        raise RuntimeError(f"Abundance/spot mismatch for {sample}")
    if not np.array_equal(abundance["spot_id"].astype(str), adata.obs_names.astype(str)):
        raise RuntimeError(f"Spot order mismatch for {sample}")
    coords = np.asarray(adata.obsm["spatial"], dtype=float)
    neighbours = six_nearest(coords)
    lib = library_size(adata)
    log_lib = np.log1p(lib)
    tymp = np.log1p(extract_gene_counts(adata, "TYMP") / lib * 10000.0)
    rng = np.random.default_rng(SEED + sample_index)

    rows = []
    for prefix, endpoint_ids in [
        ("q05", ["P01", "S02", "S03", "S04"]),
        ("mean", ["S01"]),
    ]:
        values = abundance_columns(abundance, prefix)
        malignant = values["Malignant"]
        epithelial = values["Epithelial-Normal"]
        monocyte = values["Monocyte"]
        total_macro = sum(values[state] for state in MACRO_STATES)
        positive = malignant > 0
        threshold = np.median(malignant[positive]) if positive.any() else np.inf
        eligible = positive & (malignant >= threshold)
        exposure_covariates = [
            log_lib, malignant, epithelial, total_macro, monocyte
        ]
        neighbour_log_lib = neighbour_mean(log_lib, neighbours)
        neighbour_total_macro = neighbour_mean(total_macro, neighbours)
        neighbour_monocyte = neighbour_mean(monocyte, neighbours)
        common_outcome_covariates = [
            neighbour_total_macro, neighbour_monocyte, neighbour_log_lib
        ]
        for endpoint_id in endpoint_ids:
            if endpoint_id in {"P01", "S01"}:
                outcome = neighbour_mean(values["SPP1-state Macrophage"], neighbours)
                covariates = common_outcome_covariates
                label = f"{prefix}_neighbour_spp1_state"
            elif endpoint_id == "S02":
                outcome = neighbour_mean(values["C1QC-state Macrophage"], neighbours)
                covariates = common_outcome_covariates
                label = "q05_neighbour_c1qc_state"
            elif endpoint_id == "S03":
                outcome = values["SPP1-state Macrophage"]
                covariates = [total_macro, monocyte, log_lib]
                label = "q05_same_spot_spp1_state"
            elif endpoint_id == "S04":
                fraction = values["SPP1-state Macrophage"] / np.maximum(total_macro, 1e-8)
                outcome = neighbour_mean(fraction, neighbours)
                covariates = common_outcome_covariates
                label = "q05_neighbour_spp1_fraction"
            else:
                raise RuntimeError(endpoint_id)
            rho, permutation_p, n = analyse_endpoint(
                tymp, exposure_covariates, outcome, covariates, eligible, rng
            )
            rows.append({
                "sample_id": sample,
                "preparation": preparation,
                "endpoint_id": endpoint_id,
                "endpoint": label,
                "abundance_summary": prefix,
                "eligible_spots": n,
                "spearman_rho": rho,
                "permutation_p_one_sided": permutation_p,
                "malignant_positive_median": threshold,
            })
    return rows


def random_effects_meta(group: pd.DataFrame) -> dict[str, float]:
    data = group.loc[
        np.isfinite(group["spearman_rho"])
        & (group["eligible_spots"] >= MIN_ELIGIBLE)
    ].copy()
    k = len(data)
    if k < 2:
        return {
            "sections": k,
            "pooled_rho": np.nan,
            "ci_low": np.nan,
            "ci_high": np.nan,
            "p_one_sided": np.nan,
            "tau2": np.nan,
            "positive_fraction": np.nan,
        }
    rho = np.clip(data["spearman_rho"].to_numpy(float), -0.999999, 0.999999)
    fisher = np.arctanh(rho)
    variance = 1.0 / np.maximum(data["eligible_spots"].to_numpy(float) - 3.0, 1.0)
    fixed_weights = 1.0 / variance
    fixed_mean = np.sum(fixed_weights * fisher) / np.sum(fixed_weights)
    q = np.sum(fixed_weights * (fisher - fixed_mean) ** 2)
    c_value = np.sum(fixed_weights) - (
        np.sum(fixed_weights ** 2) / np.sum(fixed_weights)
    )
    tau2 = max(0.0, (q - (k - 1.0)) / c_value) if c_value > 0 else 0.0
    weights = 1.0 / (variance + tau2)
    pooled = np.sum(weights * fisher) / np.sum(weights)
    se = np.sqrt(1.0 / np.sum(weights))
    p = float(stats.norm.sf(pooled / se))
    return {
        "sections": int(k),
        "pooled_rho": float(np.tanh(pooled)),
        "ci_low": float(np.tanh(pooled - 1.96 * se)),
        "ci_high": float(np.tanh(pooled + 1.96 * se)),
        "p_one_sided": p,
        "tau2": float(tau2),
        "positive_fraction": float(np.mean(rho > 0)),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mapping-root", required=True, type=Path)
    parser.add_argument("--section-status", required=True, type=Path)
    parser.add_argument("--preregistration", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    if not (args.mapping_root / "PHASE78_MAPPING_PASS_v1").is_file():
        raise RuntimeError("Phase78 mapping is not complete")
    args.output.mkdir(parents=True)

    status = pd.read_csv(args.section_status, sep="\t")
    status = status.loc[status["status"] == "PASS"].copy()
    if len(status) != 24:
        raise RuntimeError(f"Expected 24 passed sections, found {len(status)}")
    section_rows = []
    for index, row in enumerate(status.itertuples(index=False)):
        sample = str(row.sample_id)
        section_dir = args.mapping_root / "sections_v1" / sample
        section_rows.extend(
            analyse_section(sample, str(row.preparation), section_dir, index)
        )
    sections = pd.DataFrame(section_rows)
    if len(sections) != 24 * 5:
        raise RuntimeError(f"Unexpected section endpoint rows: {len(sections)}")

    meta_rows = []
    for (preparation, endpoint_id, endpoint), group in sections.groupby(
        ["preparation", "endpoint_id", "endpoint"], sort=True
    ):
        result = random_effects_meta(group)
        meta_rows.append({
            "preparation": preparation,
            "endpoint_id": endpoint_id,
            "endpoint": endpoint,
            **result,
        })
    meta = pd.DataFrame(meta_rows)
    meta["q"] = np.nan
    for preparation in sorted(meta["preparation"].unique()):
        mask = meta["preparation"] == preparation
        meta.loc[mask, "q"] = bh(meta.loc[mask, "p_one_sided"])
    meta["preparation_supported"] = (
        (meta["pooled_rho"] > 0)
        & (meta["q"] < 0.05)
        & (meta["positive_fraction"] >= 0.60)
    )

    primary = meta.loc[meta["endpoint_id"] == "P01"].copy()
    support_by_prep = {
        row.preparation: bool(row.preparation_supported)
        for row in primary.itertuples(index=False)
    }
    dual_preparation = bool(
        support_by_prep.get("FFPE", False)
        and support_by_prep.get("frozen", False)
    )
    decision = {
        "status": "PASS",
        "completed_at": now(),
        "sections": 24,
        "section_endpoint_rows": int(len(sections)),
        "primary_support_by_preparation": support_by_prep,
        "dual_preparation_support": dual_preparation,
        "interpretation": (
            "SPATIAL_CROSS_COMPARTMENT_SUPPORTED"
            if dual_preparation
            else "SPATIAL_CROSS_COMPARTMENT_NOT_DUAL_PREPARATION_SUPPORTED"
        ),
        "preregistration": str(args.preregistration),
        "preregistration_sha256": sha256(args.preregistration),
        "mapping_metadata_sha256": sha256(args.mapping_root / "run_metadata_v1.json"),
    }

    sections.to_csv(
        args.output / "phase78_section_effects_v1.tsv", sep="\t", index=False
    )
    meta.to_csv(
        args.output / "phase78_preparation_meta_v1.tsv", sep="\t", index=False
    )
    (args.output / "PHASE78_DECISION_v1.json").write_text(
        json.dumps(decision, indent=2) + "\n", encoding="utf-8"
    )
    checks = pd.DataFrame([
        {"check": "section_count_24", "value": status["sample_id"].nunique(), "pass": status["sample_id"].nunique() == 24},
        {"check": "ffpe_count_12", "value": int((status["preparation"] == "FFPE").sum()), "pass": int((status["preparation"] == "FFPE").sum()) == 12},
        {"check": "frozen_count_12", "value": int((status["preparation"] == "frozen").sum()), "pass": int((status["preparation"] == "frozen").sum()) == 12},
        {"check": "five_endpoints_per_section", "value": len(sections), "pass": len(sections) == 120},
        {"check": "primary_meta_two_preparations", "value": len(primary), "pass": len(primary) == 2},
        {"check": "primary_effects_finite", "value": int(np.isfinite(primary["pooled_rho"]).sum()), "pass": np.isfinite(primary["pooled_rho"]).all()},
        {"check": "q_values_valid", "value": int(((meta["q"] >= 0) & (meta["q"] <= 1)).sum()), "pass": ((meta["q"] >= 0) & (meta["q"] <= 1)).all()},
    ])
    checks.to_csv(
        args.output / "phase78_validation_checks_v1.tsv", sep="\t", index=False
    )
    if not checks["pass"].all():
        raise RuntimeError("Phase78 validation failed")
    (args.output / "PHASE78_ANALYSIS_PASS_v1").write_text(
        "status=PASS\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
