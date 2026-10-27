#!/usr/bin/env python3
"""Patient-level spatial state analysis for the TYMP--MTA-context--SPP1 axis.

This analysis tests macrophage state rather than simple SPP1+ TAM abundance.
It deliberately treats MTA as a transcriptomic context proxy and does not claim
direct metabolite measurement, mediation, or causality.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from datetime import datetime
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.sparse import issparse
from scipy.spatial import cKDTree
from scipy.stats import binomtest, rankdata, spearmanr


ROOT = Path("/mnt/BioDisk/ludingyu/project26210")
MEMORY_GATE = ROOT / "check_memory_gate.sh"
DEFAULT_MANIFEST = (
    ROOT
    / "external/ccRCC_spatial/TJ_RCC_figshare_24599295"
    / "analysis_cell2location_full_v1/input_manifest.tsv"
)
DEFAULT_OUTPUT = ROOT / "results/part4_spatial/tymp_spp1_state_v1"

SEED = 20260803
N_NEIGHBORS = 6
N_PERM = 9999
BLOCK_GRIDS = (4, 5)

SPP1_STATE_FULL = ["SPP1", "GPNMB", "LGALS3", "FABP5", "CD9", "LPL"]
SPP1_STATE_NO_SPP1 = ["GPNMB", "LGALS3", "FABP5", "CD9", "LPL"]
C1QC_STATE = ["C1QA", "C1QB", "C1QC"]
PAN_MACROPHAGE = ["CD68", "CSF1R", "LST1", "AIF1", "CTSZ", "SPI1", "TYROBP", "FCER1G"]
TUMOR_CONTEXT = ["CA9", "EPCAM", "KRT8", "KRT18", "KRT19", "NDUFA4L2"]
MTA_CONTEXT = ["MAT2A", "MAT2B", "AMD1", "SRM", "SMS"]

DOWNSTREAM = {
    "Myeloid immunosuppression": ["IL10", "TGFB1", "VSIG4", "MARCO", "CD274", "IDO1"],
    "T-cell dysfunction": ["PDCD1", "CTLA4", "LAG3", "HAVCR2", "TIGIT", "TOX", "CXCL13", "ENTPD1"],
    "Cytotoxicity": ["NKG7", "GNLY", "GZMB", "PRF1", "CTSW"],
    "Hypoxia-angiogenesis": ["CA9", "VEGFA", "HIF1A", "ADM", "ANGPTL4", "LDHA"],
    "Invasion-EMT context": ["VIM", "FN1", "MMP2", "MMP9", "ITGA5", "AXL"],
}

SIGNATURES = {
    "SPP1 state full": SPP1_STATE_FULL,
    "SPP1 state no SPP1": SPP1_STATE_NO_SPP1,
    "C1QC state": C1QC_STATE,
    "Pan-macrophage": PAN_MACROPHAGE,
    "Tumor context": TUMOR_CONTEXT,
    "MTA context": MTA_CONTEXT,
    **DOWNSTREAM,
}

OUTCOMES = ["SPP1_state_no_SPP1", "SPP1_vs_C1QC_contrast", "SPP1_expression"]
ROUTES = ["reference_dependent", "reference_light"]
ADJUSTMENTS = ["composition_qc", "composition_qc_spatial"]
SUBSETS = ["all_spots", "macrophage_enriched"]

ABUNDANCE_PREFIXES = [
    "q05cell_abundance_w_sf_means_per_cluster_mu_fg_",
    "q05cell_abundance_w_sf_",
]


def now() -> str:
    return datetime.now().astimezone().isoformat(timespec="seconds")


def sha256(path: Path, chunk: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while data := handle.read(chunk):
            digest.update(data)
    return digest.hexdigest()


def bh(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    out = np.full(values.shape, np.nan)
    finite = np.isfinite(values)
    if not finite.any():
        return out
    p = values[finite]
    order = np.argsort(p)
    ranked = p[order] * len(p) / np.arange(1, len(p) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adjusted = np.empty_like(ranked)
    adjusted[order] = np.minimum(ranked, 1.0)
    out[finite] = adjusted
    return out


def zscore(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    sd = np.nanstd(values)
    if not np.isfinite(sd) or sd <= 0:
        return np.zeros_like(values)
    return (values - np.nanmean(values)) / sd


def residualize(y: np.ndarray, covariates: list[np.ndarray]) -> np.ndarray:
    y = np.asarray(y, dtype=float)
    design = [np.ones(len(y))]
    for covariate in covariates:
        design.append(zscore(np.asarray(covariate, dtype=float)))
    x = np.column_stack(design)
    beta = np.linalg.lstsq(x, y, rcond=None)[0]
    return zscore(y - x @ beta)


def memory_gate(rows: list[dict], stage: str, sample: str = "") -> None:
    result = subprocess.run([str(MEMORY_GATE)], capture_output=True, text=True, check=False)
    line = (result.stdout or result.stderr).strip()
    rows.append({"stage": stage, "sample": sample, "exit_code": result.returncode, "line": line})
    print(f"[memory_gate] stage={stage} sample={sample} {line}", flush=True)
    if result.returncode != 0:
        raise RuntimeError(f"Memory gate blocked stage={stage}, sample={sample}: {line}")


def clean_abundance(adata: ad.AnnData) -> pd.DataFrame:
    abundance = adata.obsm["q05_cell_abundance_w_sf"]
    if not isinstance(abundance, pd.DataFrame):
        names = adata.uns.get("mod", {}).get("factor_names")
        abundance = pd.DataFrame(abundance, index=adata.obs_names, columns=names)
    abundance = abundance.copy()
    columns = []
    for value in abundance.columns:
        value = str(value)
        for prefix in ABUNDANCE_PREFIXES:
            value = value.replace(prefix, "")
        columns.append(value)
    abundance.columns = columns
    required = {
        "TYMP-High Malignant",
        "TYMP-Low Malignant",
        "SPP1+ TAM",
        "APOE+C1QC+ TAM",
        "Other Macrophage",
    }
    missing = sorted(required - set(abundance.columns))
    if missing:
        raise ValueError(f"Missing abundance columns: {missing}")
    if not np.isfinite(abundance.to_numpy(dtype=float)).all():
        raise ValueError("Non-finite abundance values")
    return abundance


def required_expression(adata: ad.AnnData) -> tuple[pd.DataFrame, np.ndarray, list[dict]]:
    requested = sorted({"TYMP", "SPP1", *[g for genes in SIGNATURES.values() for g in genes]})
    present = [gene for gene in requested if gene in adata.var_names]
    matrix = adata[:, present].X
    if issparse(matrix):
        matrix = matrix.toarray()
    matrix = np.asarray(matrix, dtype=float)
    totals = np.asarray(adata.X.sum(axis=1)).reshape(-1).astype(float)
    if np.any(totals <= 0):
        raise ValueError("Zero-count spot")
    log_cpm = np.log1p(matrix / totals[:, None] * 1e4)
    expression = pd.DataFrame(log_cpm, index=adata.obs_names, columns=present)
    qc = []
    for name, genes in SIGNATURES.items():
        used = [g for g in genes if g in expression.columns]
        if len(used) < 3:
            raise ValueError(f"Too few genes for {name}: {used}")
        detection = [(matrix[:, present.index(g)] > 0).mean() for g in used]
        qc.append(
            {
                "signature": name,
                "requested_genes": ",".join(genes),
                "used_genes": ",".join(used),
                "missing_genes": ",".join(sorted(set(genes) - set(used))),
                "n_used": len(used),
                "median_detection_fraction": float(np.median(detection)),
                "min_detection_fraction": float(np.min(detection)),
            }
        )
    return expression, totals, qc


def signature_score(expression: pd.DataFrame, genes: list[str]) -> np.ndarray:
    used = [gene for gene in genes if gene in expression.columns]
    return np.column_stack([zscore(expression[gene].to_numpy()) for gene in used]).mean(axis=1)


def neighbor_index(coords: np.ndarray) -> np.ndarray:
    _, indices = cKDTree(coords).query(coords, k=N_NEIGHBORS + 1)
    indices = np.asarray(indices[:, 1:], dtype=int)
    if indices.shape != (len(coords), N_NEIGHBORS):
        raise ValueError("Unexpected neighbor index shape")
    return indices


def lag(values: np.ndarray, neighbors: np.ndarray) -> np.ndarray:
    return np.asarray(values, dtype=float)[neighbors].mean(axis=1)


def block_labels(coords: np.ndarray, grid: int) -> np.ndarray:
    bins = []
    for axis in range(2):
        ranks = pd.Series(coords[:, axis]).rank(method="first")
        bins.append(pd.qcut(ranks, q=grid, labels=False).to_numpy())
    return bins[0] * grid + bins[1]


def block_xy(x: np.ndarray, y: np.ndarray, labels: np.ndarray, mask: np.ndarray) -> pd.DataFrame:
    frame = pd.DataFrame({"block": labels[mask], "x": x[mask], "y": y[mask]})
    return frame.groupby("block", as_index=False).agg(x=("x", "median"), y=("y", "median"), n=("x", "size"))


def permutation_spearman(x: np.ndarray, y: np.ndarray, seed: int) -> tuple[float, float]:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    rho = float(spearmanr(x, y).statistic)
    if not np.isfinite(rho):
        return rho, np.nan
    xr = rankdata(x) - np.mean(rankdata(x))
    yr = rankdata(y) - np.mean(rankdata(y))
    denom = np.sqrt(np.sum(xr**2) * np.sum(yr**2))
    rng = np.random.default_rng(seed)
    permutations = np.vstack([rng.permutation(len(y)) for _ in range(N_PERM)])
    null = (yr[permutations] @ xr) / denom
    p = (1 + np.count_nonzero(np.abs(null) >= abs(rho))) / (N_PERM + 1)
    return rho, float(p)


def standardized_coefficients(x: np.ndarray, m: np.ndarray, y: np.ndarray) -> tuple[float, float, float, float]:
    x, m, y = zscore(x), zscore(m), zscore(y)
    a = float(np.linalg.lstsq(np.column_stack([np.ones(len(x)), x]), m, rcond=None)[0][1])
    bc = np.linalg.lstsq(np.column_stack([np.ones(len(x)), x, m]), y, rcond=None)[0]
    cprime, b = float(bc[1]), float(bc[2])
    return a, b, cprime, a * b


def percentile(values: np.ndarray) -> np.ndarray:
    return pd.Series(values).rank(method="average", pct=True).to_numpy(dtype=float)


def make_spatial_figure(frame: pd.DataFrame, sample: str, output: Path) -> None:
    fig, axes = plt.subplots(1, 4, figsize=(14, 3.6), constrained_layout=True)
    panels = [
        ("neighbor_THM_refdep", "Neighbor THM", "viridis"),
        ("SPP1_state_no_SPP1", "SPP1 state (no SPP1)", "magma"),
        ("SPP1_state_no_SPP1_resid_reference_dependent_composition_qc", "Adjusted SPP1 state", "coolwarm"),
        ("interface_coactivation", "TYMP-SPP1 coactivation", "plasma"),
    ]
    for ax, (column, title, cmap) in zip(axes, panels):
        scatter = ax.scatter(frame["spatial_x"], frame["spatial_y"], c=frame[column], s=4, cmap=cmap, linewidths=0)
        ax.invert_yaxis()
        ax.set_aspect("equal")
        ax.set_title(title, fontsize=9)
        ax.axis("off")
        fig.colorbar(scatter, ax=ax, fraction=0.045, pad=0.02)
    fig.suptitle(f"{sample}: TYMP-SPP1 spatial state", fontsize=11)
    fig.savefig(output.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(output.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def cohort_summary(tests: pd.DataFrame, group_columns: list[str]) -> pd.DataFrame:
    rows = []
    for keys, frame in tests.groupby(group_columns, sort=False, dropna=False):
        if not isinstance(keys, tuple):
            keys = (keys,)
        rho = frame["rho"].to_numpy(dtype=float)
        rho = rho[np.isfinite(rho)]
        positive = int(np.sum(rho > 0))
        negative = int(np.sum(rho < 0))
        nonzero = positive + negative
        p = float(binomtest(positive, nonzero, 0.5).pvalue) if nonzero else 1.0
        row = dict(zip(group_columns, keys))
        row.update(
            patients_tested=len(rho),
            positive_patients=positive,
            negative_patients=negative,
            median_rho=float(np.median(rho)),
            q25_rho=float(np.quantile(rho, 0.25)),
            q75_rho=float(np.quantile(rho, 0.75)),
            sign_test_p=p,
            patients_within_fdr_lt_0_05=int(np.sum(frame["within_patient_fdr_bh"] < 0.05)),
        )
        rows.append(row)
    result = pd.DataFrame(rows)
    result["sign_test_fdr_bh"] = bh(result["sign_test_p"].to_numpy())
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    if args.output.exists():
        raise SystemExit(f"Refusing to overwrite existing output: {args.output}")
    args.output.mkdir(parents=True)
    tables = args.output / "tables"
    figures = args.output / "figures"
    sample_figures = figures / "samples"
    for path in [tables, figures, sample_figures]:
        path.mkdir(parents=True, exist_ok=True)

    memory_rows: list[dict] = []
    memory_gate(memory_rows, "initial")
    manifest = pd.read_csv(args.manifest, sep="\t").sort_values("sample").reset_index(drop=True)
    if len(manifest) != 12 or not {"sample", "path", "size_bytes", "sha256"}.issubset(manifest.columns):
        raise ValueError("Expected a 12-sample input manifest")

    fingerprints = []
    score_frames = []
    gene_qc_rows = []
    association_rows = []
    mediation_rows = []
    downstream_rows = []

    for sample_i, record in manifest.iterrows():
        sample = str(record["sample"])
        memory_gate(memory_rows, "sample_analysis", sample)
        input_path = Path(record["path"])
        if input_path.stat().st_size != int(record["size_bytes"]):
            raise ValueError(f"{sample}: input size changed")
        actual_hash = sha256(input_path)
        if actual_hash != str(record["sha256"]):
            raise ValueError(f"{sample}: input hash changed")
        fingerprints.append({"sample": sample, "path": str(input_path), "size_bytes": input_path.stat().st_size, "sha256": actual_hash})

        adata = ad.read_h5ad(input_path)
        if not adata.obs_names.is_unique or "spatial" not in adata.obsm:
            raise ValueError(f"{sample}: invalid H5AD structure")
        abundance = clean_abundance(adata)
        expression, totals, qc = required_expression(adata)
        for row in qc:
            row["sample"] = sample
            gene_qc_rows.append(row)
        coords = np.asarray(adata.obsm["spatial"], dtype=float)
        neighbors = neighbor_index(coords)
        labels = {grid: block_labels(coords, grid) for grid in BLOCK_GRIDS}

        scores = {name: signature_score(expression, genes) for name, genes in SIGNATURES.items()}
        spp1_expression = zscore(expression["SPP1"].to_numpy(dtype=float))
        spp1_contrast = zscore(scores["SPP1 state no SPP1"] - scores["C1QC state"])
        mta_context = zscore(scores["MTA context"])
        log_umi = np.log1p(totals)
        n_genes = np.asarray((adata.X > 0).sum(axis=1)).reshape(-1)
        log_genes = np.log1p(n_genes)

        macro_columns = ["SPP1+ TAM", "APOE+C1QC+ TAM", "Other Macrophage"]
        immune_columns = [column for column in abundance.columns if column not in {"TYMP-High Malignant", "TYMP-Low Malignant", "Epithelial-Normal", "Endothelial", "SMC/Pericyte"}]
        total_macro = abundance[macro_columns].sum(axis=1).to_numpy(dtype=float)
        total_immune = abundance[immune_columns].sum(axis=1).to_numpy(dtype=float)
        total_malignant = abundance[["TYMP-High Malignant", "TYMP-Low Malignant"]].sum(axis=1).to_numpy(dtype=float)
        thm = abundance["TYMP-High Malignant"].to_numpy(dtype=float)
        refdep_cov = [np.log1p(total_macro), np.log1p(total_immune), np.log1p(total_malignant), log_umi, log_genes]
        refl_cov = [scores["Pan-macrophage"], scores["Tumor context"], log_umi, log_genes]
        sx, sy = zscore(coords[:, 0]), zscore(coords[:, 1])
        spatial_cov = [sx, sy, sx**2, sy**2, sx * sy]

        tymp_expression = zscore(expression["TYMP"].to_numpy(dtype=float))
        tymp_expression_resid = residualize(tymp_expression, refl_cov)
        exposures = {
            "reference_dependent": lag(zscore(np.log1p(thm)), neighbors),
            "reference_light": lag(tymp_expression_resid, neighbors),
        }
        raw_outcomes = {
            "SPP1_state_no_SPP1": zscore(scores["SPP1 state no SPP1"]),
            "SPP1_vs_C1QC_contrast": spp1_contrast,
            "SPP1_expression": spp1_expression,
        }
        route_cov = {"reference_dependent": refdep_cov, "reference_light": refl_cov}
        route_macro = {
            "reference_dependent": total_macro >= np.median(total_macro),
            "reference_light": scores["Pan-macrophage"] >= np.median(scores["Pan-macrophage"]),
        }
        residuals = {}
        mta_residuals = {}
        for route in ROUTES:
            for adjustment in ADJUSTMENTS:
                cov = route_cov[route] + (spatial_cov if adjustment.endswith("spatial") else [])
                mta_residuals[(route, adjustment)] = residualize(mta_context, cov)
                for outcome, values in raw_outcomes.items():
                    residuals[(route, adjustment, outcome)] = residualize(values, cov)

        for route_i, route in enumerate(ROUTES):
            for adjustment_i, adjustment in enumerate(ADJUSTMENTS):
                for subset_i, subset in enumerate(SUBSETS):
                    mask = np.ones(adata.n_obs, dtype=bool) if subset == "all_spots" else route_macro[route]
                    for outcome_i, outcome in enumerate(OUTCOMES):
                        y = residuals[(route, adjustment, outcome)]
                        for grid in BLOCK_GRIDS:
                            summary = block_xy(exposures[route], y, labels[grid], mask)
                            rho, p = permutation_spearman(
                                summary["x"].to_numpy(), summary["y"].to_numpy(),
                                SEED + sample_i * 100000 + route_i * 10000 + adjustment_i * 4000 + subset_i * 1000 + outcome_i * 100 + grid,
                            )
                            association_rows.append(
                                {"sample": sample, "route": route, "adjustment": adjustment, "subset": subset,
                                 "exposure": "neighbor_TYMP_state", "outcome": outcome, "grid": f"{grid}x{grid}",
                                 "n_spots": int(mask.sum()), "n_blocks": len(summary), "rho": rho, "p_value": p,
                                 "permutations": N_PERM}
                            )

        for route_i, route in enumerate(ROUTES):
            mask = route_macro[route]
            for grid in BLOCK_GRIDS:
                frame = pd.DataFrame({"block": labels[grid][mask], "x": exposures[route][mask],
                                      "m": mta_residuals[(route, "composition_qc")][mask],
                                      "y": residuals[(route, "composition_qc", "SPP1_state_no_SPP1")][mask]})
                block = frame.groupby("block", as_index=False).median(numeric_only=True)
                a, b, cprime, indirect = standardized_coefficients(block["x"], block["m"], block["y"])
                mediation_rows.append({"sample": sample, "route": route, "subset": "macrophage_enriched",
                                       "adjustment": "composition_qc", "grid": f"{grid}x{grid}", "n_blocks": len(block),
                                       "path_a_tymp_to_mta": a, "path_b_mta_to_spp1": b,
                                       "direct_cprime": cprime, "indirect_a_times_b": indirect})

        primary_state = residuals[("reference_dependent", "composition_qc", "SPP1_state_no_SPP1")]
        coactivation = percentile(exposures["reference_dependent"]) * percentile(primary_state)
        primary_mask = route_macro["reference_dependent"]
        for target_i, (target, values) in enumerate(DOWNSTREAM.items()):
            downstream_resid = residualize(scores[target], refdep_cov)
            for grid in BLOCK_GRIDS:
                summary = block_xy(coactivation, downstream_resid, labels[grid], primary_mask)
                rho, p = permutation_spearman(summary["x"], summary["y"], SEED + 900000 + sample_i * 10000 + target_i * 100 + grid)
                downstream_rows.append({"sample": sample, "route": "reference_dependent", "subset": "macrophage_enriched",
                                        "adjustment": "composition_qc", "interface": "TYMP_SPP1_coactivation",
                                        "outcome": target, "grid": f"{grid}x{grid}", "n_spots": int(primary_mask.sum()),
                                        "n_blocks": len(summary), "rho": rho, "p_value": p, "permutations": N_PERM})

        frame = pd.DataFrame(index=adata.obs_names)
        frame["sample"] = sample
        frame["barcode"] = adata.obs_names
        frame["spatial_x"], frame["spatial_y"] = coords[:, 0], coords[:, 1]
        frame["block_4x4"], frame["block_5x5"] = labels[4], labels[5]
        frame["total_macro_q05"] = total_macro
        frame["total_malignant_q05"] = total_malignant
        frame["pan_macrophage_score"] = scores["Pan-macrophage"]
        frame["TYMP_expression"] = tymp_expression
        frame["neighbor_THM_refdep"] = exposures["reference_dependent"]
        frame["neighbor_TYMP_reflight"] = exposures["reference_light"]
        frame["MTA_context"] = mta_context
        frame["SPP1_state_no_SPP1"] = raw_outcomes["SPP1_state_no_SPP1"]
        frame["SPP1_vs_C1QC_contrast"] = raw_outcomes["SPP1_vs_C1QC_contrast"]
        frame["SPP1_expression"] = raw_outcomes["SPP1_expression"]
        frame["macro_enriched_reference_dependent"] = route_macro["reference_dependent"]
        frame["macro_enriched_reference_light"] = route_macro["reference_light"]
        for route in ROUTES:
            for adjustment in ADJUSTMENTS:
                frame[f"MTA_context_resid_{route}_{adjustment}"] = mta_residuals[(route, adjustment)]
                for outcome in OUTCOMES:
                    frame[f"{outcome}_resid_{route}_{adjustment}"] = residuals[(route, adjustment, outcome)]
        frame["interface_coactivation"] = coactivation
        score_frames.append(frame.reset_index(drop=True))
        make_spatial_figure(frame, sample, sample_figures / f"{sample}_tymp_spp1_state_v1")
        del adata

    memory_gate(memory_rows, "cohort_synthesis")
    scores_all = pd.concat(score_frames, ignore_index=True)
    associations = pd.DataFrame(association_rows)
    associations["within_patient_fdr_bh"] = np.nan
    for sample, indices in associations.groupby("sample").groups.items():
        associations.loc[indices, "within_patient_fdr_bh"] = bh(associations.loc[indices, "p_value"].to_numpy())
    association_summary = cohort_summary(associations, ["route", "adjustment", "subset", "outcome", "grid"])

    mediation = pd.DataFrame(mediation_rows)
    mediation_summary_rows = []
    for (route, grid), frame in mediation.groupby(["route", "grid"]):
        row = {"route": route, "grid": grid, "patients_tested": len(frame)}
        for column in ["path_a_tymp_to_mta", "path_b_mta_to_spp1", "direct_cprime", "indirect_a_times_b"]:
            values = frame[column].to_numpy(dtype=float)
            positive = int(np.sum(values > 0)); negative = int(np.sum(values < 0)); nonzero = positive + negative
            row[f"{column}_positive_patients"] = positive
            row[f"{column}_median"] = float(np.median(values))
            row[f"{column}_sign_p"] = float(binomtest(positive, nonzero, 0.5).pvalue) if nonzero else 1.0
        mediation_summary_rows.append(row)
    mediation_summary = pd.DataFrame(mediation_summary_rows)
    p_columns = [c for c in mediation_summary.columns if c.endswith("_sign_p")]
    all_p = np.concatenate([mediation_summary[c].to_numpy() for c in p_columns])
    all_fdr = bh(all_p); cursor = 0
    for column in p_columns:
        mediation_summary[column.replace("_sign_p", "_sign_fdr_bh")] = all_fdr[cursor:cursor + len(mediation_summary)]
        cursor += len(mediation_summary)

    downstream_tests = pd.DataFrame(downstream_rows)
    downstream_tests["within_patient_fdr_bh"] = np.nan
    for sample, indices in downstream_tests.groupby("sample").groups.items():
        downstream_tests.loc[indices, "within_patient_fdr_bh"] = bh(downstream_tests.loc[indices, "p_value"].to_numpy())
    downstream_summary = cohort_summary(downstream_tests, ["outcome", "grid"])

    pd.DataFrame(fingerprints).to_csv(tables / "input_fingerprints_v1.tsv", sep="\t", index=False)
    pd.DataFrame(gene_qc_rows).to_csv(tables / "signature_gene_qc_v1.tsv", sep="\t", index=False)
    scores_all.to_csv(tables / "patient_spot_state_scores_v1.tsv.gz", sep="\t", index=False, compression="gzip")
    associations.to_csv(tables / "patient_state_associations_v1.tsv", sep="\t", index=False)
    association_summary.to_csv(tables / "cohort_state_direction_summary_v1.tsv", sep="\t", index=False)
    mediation.to_csv(tables / "patient_exploratory_mediation_v1.tsv", sep="\t", index=False)
    mediation_summary.to_csv(tables / "cohort_exploratory_mediation_summary_v1.tsv", sep="\t", index=False)
    downstream_tests.to_csv(tables / "patient_interface_downstream_tests_v1.tsv", sep="\t", index=False)
    downstream_summary.to_csv(tables / "cohort_interface_downstream_summary_v1.tsv", sep="\t", index=False)
    pd.DataFrame(memory_rows).to_csv(args.output / "memory_gates_v1.tsv", sep="\t", index=False)

    memory_gate(memory_rows, "cohort_figures")
    primary = associations[(associations.route == "reference_dependent") & (associations.adjustment == "composition_qc") &
                           (associations.subset == "macrophage_enriched") & (associations.outcome == "SPP1_state_no_SPP1")]
    fig, ax = plt.subplots(figsize=(8, 4.8), constrained_layout=True)
    sns.pointplot(data=primary, x="sample", y="rho", hue="grid", dodge=0.35, linestyles="none", errorbar=None, ax=ax)
    ax.axhline(0, color="black", lw=0.8); ax.set_ylabel("Spatial-block Spearman rho"); ax.set_xlabel("Patient")
    ax.set_title("Adjusted SPP1 state vs neighboring TYMP-high malignant state")
    ax.tick_params(axis="x", rotation=45)
    fig.savefig(figures / "cohort_primary_state_patient_effects_v1.png", dpi=300, bbox_inches="tight")
    fig.savefig(figures / "cohort_primary_state_patient_effects_v1.pdf", bbox_inches="tight"); plt.close(fig)

    pivot = mediation.pivot(index="sample", columns=["route", "grid"], values="indirect_a_times_b")
    fig, ax = plt.subplots(figsize=(7, 5), constrained_layout=True)
    sns.heatmap(pivot, cmap="coolwarm", center=0, annot=True, fmt=".2f", ax=ax)
    ax.set_title("Exploratory TYMP -> MTA-context -> SPP1-state indirect effect")
    fig.savefig(figures / "cohort_exploratory_mediation_heatmap_v1.png", dpi=300, bbox_inches="tight")
    fig.savefig(figures / "cohort_exploratory_mediation_heatmap_v1.pdf", bbox_inches="tight"); plt.close(fig)

    pivot = downstream_summary.pivot(index="outcome", columns="grid", values="median_rho")
    annot = downstream_summary.pivot(index="outcome", columns="grid", values="positive_patients").astype(str) + "/12"
    fig, ax = plt.subplots(figsize=(5.8, 4.2), constrained_layout=True)
    sns.heatmap(pivot, cmap="coolwarm", center=0, annot=annot, fmt="", ax=ax, cbar_kws={"label": "Median patient rho"})
    ax.set_title("TYMP-SPP1 interface and downstream states")
    fig.savefig(figures / "cohort_interface_downstream_summary_v1.png", dpi=300, bbox_inches="tight")
    fig.savefig(figures / "cohort_interface_downstream_summary_v1.pdf", bbox_inches="tight"); plt.close(fig)

    route_compare = associations[(associations.adjustment == "composition_qc") & (associations.subset == "macrophage_enriched") &
                                 (associations.outcome == "SPP1_state_no_SPP1")].pivot(index=["sample", "grid"], columns="route", values="rho").reset_index()
    fig, ax = plt.subplots(figsize=(5.2, 4.6), constrained_layout=True)
    sns.scatterplot(data=route_compare, x="reference_dependent", y="reference_light", hue="grid", s=55, ax=ax)
    lo = float(np.nanmin(route_compare[["reference_dependent", "reference_light"]])); hi = float(np.nanmax(route_compare[["reference_dependent", "reference_light"]]))
    ax.plot([lo, hi], [lo, hi], ls="--", color="grey", lw=0.8); ax.axhline(0, color="grey", lw=0.5); ax.axvline(0, color="grey", lw=0.5)
    ax.set_title("Reference-dependent vs reference-light effects")
    fig.savefig(figures / "cohort_route_concordance_v1.png", dpi=300, bbox_inches="tight")
    fig.savefig(figures / "cohort_route_concordance_v1.pdf", bbox_inches="tight"); plt.close(fig)

    metadata = {
        "status": "completed", "created_at": now(), "script": str(Path(__file__).resolve()),
        "manifest": str(args.manifest), "output": str(args.output), "samples": manifest["sample"].tolist(),
        "n_patients": 12, "n_spots": int(len(scores_all)), "seed": SEED, "neighbors": N_NEIGHBORS,
        "block_grids": list(BLOCK_GRIDS), "permutations": N_PERM, "signatures": SIGNATURES,
        "primary_analysis": {"route": "reference_dependent", "adjustment": "composition_qc",
                             "subset": "macrophage_enriched", "outcome": "SPP1_state_no_SPP1"},
        "limitations": [
            "MTA context is a transcriptomic proxy, not direct 5-prime-MTA measurement.",
            "Visium spots mix cell types; adjusted state scores are not single-cell measurements.",
            "Exploratory mediation does not establish biochemical mediation or causality.",
            "Patient direction tests treat patients, not spots or spatial blocks, as biological replicates.",
            "Reference-dependent results retain the legacy cell2location reference; reference-light results are a sensitivity analysis.",
        ],
    }
    (args.output / "run_metadata_v1.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    memory_gate(memory_rows, "file_manifest")
    pd.DataFrame(memory_rows).to_csv(args.output / "memory_gates_v1.tsv", sep="\t", index=False)
    files = []
    for path in sorted(args.output.rglob("*")):
        if path.is_file() and path.name != "file_manifest_v1.tsv":
            files.append({"path": str(path.relative_to(args.output)), "size_bytes": path.stat().st_size, "sha256": sha256(path)})
    pd.DataFrame(files).to_csv(args.output / "file_manifest_v1.tsv", sep="\t", index=False)
    print(f"TYMP_SPP1_STATE_ANALYSIS_ALL_PASS patients=12 spots={len(scores_all)} output={args.output}", flush=True)


if __name__ == "__main__":
    main()
