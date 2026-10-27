#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import subprocess
from datetime import datetime
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import binomtest, rankdata, spearmanr
from scipy.spatial import cKDTree


ROOT = Path("/mnt/BioDisk/ludingyu/project26210")
OUT = Path("/mnt/BioDisk/wangmingyuan/project26210_tymp_vascular_niche_stage_v1/output_v1")
MEMORY_GATE = ROOT / "check_memory_gate.sh"
TJ_ROOT = ROOT / "external/ccRCC_spatial/TJ_RCC_figshare_24599295"
MANIFEST = TJ_ROOT / "analysis_cell2location_full_v1/input_manifest.tsv"
STATE_SCORES = ROOT / "results/part4_spatial/tymp_spp1_state_v1/tables/patient_spot_state_scores_v1.tsv.gz"
GLMF_SPOT_SCORES = ROOT / "results/part4_spatial/tymp_immune_ecology_v2/tables/original_spot_scores.tsv.gz"
GLMF_H5AD = {
    "GLMF1": ROOT / "results/part4_spatial/analysis/GLMF1_analyzed.h5ad",
    "GLMF2": ROOT / "results/part4_spatial/analysis/GLMF2_analyzed.h5ad",
}
GLMF_BOUNDARY = {
    "GLMF1": ROOT / "results/part4_spatial/advanced_patch_bayesspace_glmf1_v1/GLMF1_computational_boundary_spots.tsv",
    "GLMF2": ROOT / "results/part4_spatial/advanced_patch_bayesspace_glmf2_v1/GLMF2_computational_boundary_spots.tsv",
}

SEED = 20260803
N_PERM = 9999
BLOCK_GRIDS = (4, 5)
N_NEIGHBORS = 6
EPS = 1e-6

SIGNATURES = {
    "Cytotoxicity": ["NKG7", "GNLY", "GZMB", "PRF1", "CTSW"],
    "Hypoxia-angiogenesis": ["CA9", "VEGFA", "HIF1A", "ADM", "ANGPTL4", "LDHA"],
    "Invasion-EMT": ["VIM", "FN1", "MMP2", "MMP9", "ITGA5", "AXL"],
}

RELATIONS = [
    ("SPP1_TAM_to_vascular_compartment", "route_independent", "spp1_percentile", "vascular_compartment", 1),
    ("SPP1_vascular_niche_to_T_cell_balance", "route_independent", "spp1_vascular_niche", "T_cell_balance", -1),
    ("SPP1_vascular_niche_to_cytotoxicity", "route_independent", "spp1_vascular_niche", "cytotoxicity_residual", -1),
    ("SPP1_vascular_niche_to_hypoxia_angiogenesis", "route_independent", "spp1_vascular_niche", "hypoxia_angiogenesis_residual", 1),
    ("SPP1_vascular_niche_to_invasion_EMT", "route_independent", "spp1_vascular_niche", "invasion_EMT_residual", 1),
    ("TYMP_MTA_core_to_SPP1_vascular_niche", "reference_light", "TYMP_MTA_core_reference_light", "spp1_vascular_niche", 1),
    ("TYMP_MTA_core_to_SPP1_monocyte_balance", "reference_light", "TYMP_MTA_core_reference_light", "SPP1_monocyte_balance", 1),
    ("TYMP_MTA_core_to_SPP1_vascular_niche", "reference_dependent", "TYMP_MTA_core_reference_dependent", "spp1_vascular_niche", 1),
    ("TYMP_MTA_core_to_SPP1_monocyte_balance", "reference_dependent", "TYMP_MTA_core_reference_dependent", "SPP1_monocyte_balance", 1),
]

ZONE_ORDER = ["Tumor core", "Inner interface", "Outer interface", "Peripheral 400-800um", "Distal >800um"]


def now() -> str:
    return datetime.now().astimezone().isoformat(timespec="seconds")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def memory_gate(rows: list[dict], stage: str, sample: str = "") -> None:
    result = subprocess.run([str(MEMORY_GATE)], capture_output=True, text=True, check=False)
    line = (result.stdout or result.stderr).strip()
    rows.append({"time": now(), "stage": stage, "sample": sample, "exit_code": result.returncode, "line": line})
    print(f"[memory_gate] {stage} {sample} {line}", flush=True)
    if result.returncode != 0:
        raise RuntimeError(f"Memory gate failed: stage={stage}; sample={sample}; {line}")


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
        return np.zeros(len(values), dtype=float)
    return (values - np.nanmean(values)) / sd


def percentile(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    return (rankdata(values, method="average") - 0.5) / len(values)


def residualize(y: np.ndarray, covariates: list[np.ndarray]) -> np.ndarray:
    y = zscore(y)
    design = [np.ones(len(y))]
    for covariate in covariates:
        design.append(zscore(covariate))
    x = np.column_stack(design)
    beta = np.linalg.lstsq(x, y, rcond=None)[0]
    return zscore(y - x @ beta)


def clean_abundance(adata: ad.AnnData) -> pd.DataFrame:
    abundance = adata.obsm["q05_cell_abundance_w_sf"]
    if not isinstance(abundance, pd.DataFrame):
        names = adata.uns.get("mod", {}).get("factor_names")
        abundance = pd.DataFrame(abundance, index=adata.obs_names, columns=names)
    abundance = abundance.copy()
    prefixes = ["q05cell_abundance_w_sf_means_per_cluster_mu_fg_", "q05cell_abundance_w_sf_"]
    cleaned = []
    for column in abundance.columns:
        value = str(column)
        for prefix in prefixes:
            value = value.replace(prefix, "")
        cleaned.append(value)
    abundance.columns = cleaned
    required = {
        "Endothelial", "SMC/Pericyte", "SPP1+ TAM", "Monocyte", "T cell",
        "TYMP-High Malignant", "TYMP-Low Malignant", "APOE+C1QC+ TAM", "Other Macrophage",
    }
    missing = sorted(required - set(abundance.columns))
    if missing:
        raise ValueError(f"Missing abundance columns: {missing}")
    if not abundance.index.equals(adata.obs_names):
        abundance = abundance.reindex(adata.obs_names)
    if abundance.isna().any().any() or not np.isfinite(abundance.to_numpy(dtype=float)).all():
        raise ValueError("Abundance matrix contains missing/non-finite values")
    return abundance


def expression_signatures(adata: ad.AnnData) -> tuple[dict[str, np.ndarray], np.ndarray, np.ndarray, list[dict]]:
    requested = sorted({gene for genes in SIGNATURES.values() for gene in genes})
    present = [gene for gene in requested if gene in adata.var_names]
    matrix = adata[:, present].X
    if sparse.issparse(matrix):
        matrix = matrix.toarray()
    matrix = np.asarray(matrix, dtype=float)
    totals = np.asarray(adata.X.sum(axis=1)).reshape(-1).astype(float)
    detected = np.asarray((adata.X > 0).sum(axis=1)).reshape(-1).astype(float)
    if np.any(totals <= 0):
        raise ValueError("Zero-count spots found")
    log_cpm = np.log1p(matrix / totals[:, None] * 1e4)
    expression = pd.DataFrame(log_cpm, index=adata.obs_names, columns=present)
    scores: dict[str, np.ndarray] = {}
    qc: list[dict] = []
    for signature, genes in SIGNATURES.items():
        used = [gene for gene in genes if gene in expression.columns]
        if len(used) < 3:
            raise ValueError(f"Too few genes for {signature}: {used}")
        scores[signature] = np.column_stack([zscore(expression[gene].to_numpy()) for gene in used]).mean(axis=1)
        detection = [(matrix[:, present.index(gene)] > 0).mean() for gene in used]
        qc.append({
            "signature": signature,
            "requested_genes": ";".join(genes),
            "used_genes": ";".join(used),
            "missing_genes": ";".join(sorted(set(genes) - set(used))),
            "n_used": len(used),
            "median_detection_fraction": float(np.median(detection)),
            "min_detection_fraction": float(np.min(detection)),
        })
    return scores, totals, detected, qc


def neighbor_lag(values: np.ndarray, coords: np.ndarray) -> np.ndarray:
    _, indices = cKDTree(coords).query(coords, k=N_NEIGHBORS + 1)
    return np.asarray(values, dtype=float)[np.asarray(indices[:, 1:], dtype=int)].mean(axis=1)


def block_labels(coords: np.ndarray, grid: int) -> np.ndarray:
    bins = []
    for axis in range(2):
        ranks = pd.Series(coords[:, axis]).rank(method="first")
        bins.append(pd.qcut(ranks, q=grid, labels=False).to_numpy())
    return bins[0] * grid + bins[1]


def block_medians(x: np.ndarray, y: np.ndarray, labels: np.ndarray, mask: np.ndarray | None = None) -> pd.DataFrame:
    if mask is None:
        mask = np.ones(len(x), dtype=bool)
    frame = pd.DataFrame({"block": labels[mask], "x": np.asarray(x)[mask], "y": np.asarray(y)[mask]})
    return frame.groupby("block", as_index=False).agg(x=("x", "median"), y=("y", "median"), n=("x", "size"))


def permutation_spearman(x: np.ndarray, y: np.ndarray, seed: int) -> tuple[float, float]:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if len(x) < 8 or np.nanstd(x) <= 0 or np.nanstd(y) <= 0:
        return np.nan, np.nan
    xr = zscore(rankdata(x))
    yr = zscore(rankdata(y))
    observed = float(np.mean(xr * yr))
    rng = np.random.default_rng(seed)
    exceed = 0
    remaining = N_PERM
    while remaining:
        batch = min(1000, remaining)
        permutations = np.vstack([rng.permutation(len(yr)) for _ in range(batch)])
        null = np.mean(xr[None, :] * yr[permutations], axis=1)
        exceed += int(np.sum(np.abs(null) >= abs(observed) - 1e-15))
        remaining -= batch
    return observed, (exceed + 1) / (N_PERM + 1)


def make_common_frame(
    sample: str,
    adata: ad.AnnData,
    source: pd.DataFrame,
    source_mode: str,
    boundary_path: Path,
) -> tuple[pd.DataFrame, list[dict]]:
    abundance = clean_abundance(adata)
    scores, totals, detected, gene_qc = expression_signatures(adata)
    coords = np.asarray(adata.obsm["spatial"], dtype=float)
    if coords.shape != (adata.n_obs, 2) or not np.isfinite(coords).all():
        raise ValueError(f"{sample}: invalid spatial coordinates")

    source = source.copy().set_index("barcode").reindex(adata.obs_names)
    if source.isna().all(axis=1).any():
        raise ValueError(f"{sample}: missing source score rows")
    boundary = pd.read_csv(boundary_path, sep="\t").set_index("barcode").reindex(adata.obs_names)
    if boundary[["signed_distance_um", "distance_zone"]].isna().any().any():
        raise ValueError(f"{sample}: missing boundary rows")

    total_immune_columns = [
        column for column in abundance.columns
        if column not in {"TYMP-High Malignant", "TYMP-Low Malignant", "Epithelial-Normal", "Endothelial", "SMC/Pericyte"}
    ]
    total_immune = abundance[total_immune_columns].sum(axis=1).to_numpy(dtype=float)
    total_malignant = abundance[["TYMP-High Malignant", "TYMP-Low Malignant"]].sum(axis=1).to_numpy(dtype=float)
    t_cell = abundance["T cell"].to_numpy(dtype=float)
    other_immune = np.maximum(total_immune - t_cell, 0)
    spp1 = abundance["SPP1+ TAM"].to_numpy(dtype=float)
    monocyte = abundance["Monocyte"].to_numpy(dtype=float)
    endothelial = abundance["Endothelial"].to_numpy(dtype=float)
    pericyte = abundance["SMC/Pericyte"].to_numpy(dtype=float)
    thm = abundance["TYMP-High Malignant"].to_numpy(dtype=float)

    p_spp1 = percentile(spp1)
    p_endo = percentile(endothelial)
    p_pericyte = percentile(pericyte)
    vascular = np.sqrt(np.maximum(p_endo * p_pericyte, EPS))
    niche = np.cbrt(np.maximum(p_spp1 * p_endo * p_pericyte, EPS))
    spp1_monocyte = residualize(
        zscore(np.log1p(spp1)) - zscore(np.log1p(monocyte)),
        [np.log1p(spp1 + monocyte + abundance["APOE+C1QC+ TAM"].to_numpy() + abundance["Other Macrophage"].to_numpy()), np.log1p(total_immune)],
    )
    t_balance = zscore(np.log1p(t_cell)) - zscore(np.log1p(other_immune))
    sx, sy = zscore(coords[:, 0]), zscore(coords[:, 1])
    spatial_cov = [sx, sy, sx**2, sy**2, sx * sy]
    quality_cov = [np.log1p(totals), np.log1p(detected), np.log1p(total_immune)]
    cytotoxicity = residualize(scores["Cytotoxicity"], [np.log1p(t_cell), *quality_cov, *spatial_cov])
    hypoxia = residualize(scores["Hypoxia-angiogenesis"], [np.log1p(total_malignant), *quality_cov, *spatial_cov])
    invasion = residualize(scores["Invasion-EMT"], [np.log1p(total_malignant), *quality_cov, *spatial_cov])

    if source_mode == "tjrcc":
        tymp_light = source["neighbor_TYMP_reflight"].to_numpy(dtype=float)
        tymp_dependent = source["neighbor_THM_refdep"].to_numpy(dtype=float)
        mta = source["MTA_context"].to_numpy(dtype=float)
    elif source_mode == "glmf":
        tymp_dependent = neighbor_lag(zscore(np.log1p(thm)), coords)
        tymp_light = tymp_dependent.copy()
        mta = source["score_polyamine_mta_context"].to_numpy(dtype=float)
    else:
        raise ValueError(source_mode)

    core_light = np.sqrt(np.maximum(percentile(tymp_light) * percentile(mta), EPS))
    core_dependent = np.sqrt(np.maximum(percentile(tymp_dependent) * percentile(mta), EPS))
    frame = pd.DataFrame({
        "sample": sample,
        "cohort": "TJ-RCC" if source_mode == "tjrcc" else "GLMF",
        "barcode": adata.obs_names,
        "spatial_x": coords[:, 0],
        "spatial_y": coords[:, 1],
        "block_4x4": block_labels(coords, 4),
        "block_5x5": block_labels(coords, 5),
        "TYMP_MTA_core_reference_light": core_light,
        "TYMP_MTA_core_reference_dependent": core_dependent,
        "MTA_context": zscore(mta),
        "spp1_percentile": p_spp1,
        "endothelial_percentile": p_endo,
        "pericyte_percentile": p_pericyte,
        "vascular_compartment": vascular,
        "spp1_vascular_niche": niche,
        "SPP1_monocyte_balance": spp1_monocyte,
        "T_cell_balance": t_balance,
        "cytotoxicity_residual": cytotoxicity,
        "hypoxia_angiogenesis_residual": hypoxia,
        "invasion_EMT_residual": invasion,
        "SPP1_TAM_q05": spp1,
        "Monocyte_q05": monocyte,
        "Endothelial_q05": endothelial,
        "SMC_Pericyte_q05": pericyte,
        "T_cell_q05": t_cell,
        "total_immune_q05": total_immune,
        "total_malignant_q05": total_malignant,
        "signed_distance_um": boundary["signed_distance_um"].to_numpy(dtype=float),
        "distance_zone": boundary["distance_zone"].astype(str).to_numpy(),
        "computational_boundary_spot": boundary["computational_boundary_spot"].astype(bool).to_numpy(),
    })
    for row in gene_qc:
        row.update({"sample": sample, "cohort": frame["cohort"].iat[0]})
    return frame, gene_qc


def association_tests(frame: pd.DataFrame, sample_index: int) -> list[dict]:
    rows = []
    cohort = frame["cohort"].iat[0]
    allowed_routes = {"route_independent", "reference_dependent"} if cohort == "GLMF" else {"route_independent", "reference_light", "reference_dependent"}
    for relation_index, (relation, route, x_col, y_col, expected) in enumerate(RELATIONS):
        if route not in allowed_routes:
            continue
        for grid in BLOCK_GRIDS:
            block = block_medians(frame[x_col].to_numpy(), frame[y_col].to_numpy(), frame[f"block_{grid}x{grid}"].to_numpy())
            rho, p_value = permutation_spearman(block["x"].to_numpy(), block["y"].to_numpy(), SEED + sample_index * 10000 + relation_index * 100 + grid)
            rows.append({
                "sample": frame["sample"].iat[0], "cohort": cohort, "relation": relation, "route": route,
                "source": x_col, "target": y_col, "expected_direction": "+" if expected > 0 else "-",
                "grid": f"{grid}x{grid}", "n_spots": len(frame), "n_blocks": len(block),
                "rho": rho, "direction_aligned_rho": rho * expected if np.isfinite(rho) else np.nan,
                "p_value": p_value, "permutations": N_PERM,
            })
    return rows


def boundary_tables(frame: pd.DataFrame) -> tuple[list[dict], list[dict]]:
    score_columns = [
        "TYMP_MTA_core_reference_light", "TYMP_MTA_core_reference_dependent", "spp1_vascular_niche",
        "SPP1_monocyte_balance", "T_cell_balance", "cytotoxicity_residual",
        "hypoxia_angiogenesis_residual", "invasion_EMT_residual",
    ]
    profile_rows = []
    contrast_rows = []
    for score in score_columns:
        standardized = zscore(frame[score].to_numpy())
        for zone in ZONE_ORDER:
            mask = frame["distance_zone"].eq(zone).to_numpy()
            if mask.any():
                profile_rows.append({
                    "sample": frame["sample"].iat[0], "cohort": frame["cohort"].iat[0], "score": score,
                    "distance_zone": zone, "zone_order": ZONE_ORDER.index(zone), "n_spots": int(mask.sum()),
                    "median_z": float(np.median(standardized[mask])), "mean_z": float(np.mean(standardized[mask])),
                })
        interface = frame["distance_zone"].isin(["Inner interface", "Outer interface"]).to_numpy()
        for reference_zone in ["Tumor core", "Distal >800um"]:
            reference = frame["distance_zone"].eq(reference_zone).to_numpy()
            if interface.any() and reference.any():
                contrast_rows.append({
                    "sample": frame["sample"].iat[0], "cohort": frame["cohort"].iat[0], "score": score,
                    "contrast": f"interface_minus_{reference_zone.replace(' ', '_')}",
                    "n_interface": int(interface.sum()), "n_reference": int(reference.sum()),
                    "interface_median_z": float(np.median(standardized[interface])),
                    "reference_median_z": float(np.median(standardized[reference])),
                    "difference_z": float(np.median(standardized[interface]) - np.median(standardized[reference])),
                })
    return profile_rows, contrast_rows


def cohort_direction_summary(patient_tests: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for keys, group in patient_tests.groupby(["relation", "route", "grid", "expected_direction"], sort=True):
        relation, route, grid, expected = keys
        aligned = group["direction_aligned_rho"].dropna().to_numpy(dtype=float)
        positive = int((aligned > 0).sum())
        p_value = binomtest(positive, len(aligned), 0.5, alternative="two-sided").pvalue if len(aligned) else np.nan
        rows.append({
            "relation": relation, "route": route, "grid": grid, "expected_direction": expected,
            "n_patients": len(aligned), "direction_supported_patients": positive,
            "median_rho": float(group["rho"].median()), "median_direction_aligned_rho": float(np.median(aligned)),
            "sign_test_p": p_value,
        })
    result = pd.DataFrame(rows)
    result["sign_test_fdr_bh_global"] = bh(result["sign_test_p"].to_numpy())
    return result


def boundary_cohort_summary(contrasts: pd.DataFrame) -> pd.DataFrame:
    expected = {
        "spp1_vascular_niche": 1,
        "SPP1_monocyte_balance": 1,
        "T_cell_balance": -1,
        "cytotoxicity_residual": -1,
        "hypoxia_angiogenesis_residual": 1,
        "invasion_EMT_residual": 1,
        "TYMP_MTA_core_reference_light": 1,
        "TYMP_MTA_core_reference_dependent": 1,
    }
    rows = []
    for (score, contrast), group in contrasts.groupby(["score", "contrast"], sort=True):
        direction = expected[score]
        values = group["difference_z"].to_numpy(dtype=float)
        aligned = values * direction
        positive = int((aligned > 0).sum())
        rows.append({
            "score": score, "contrast": contrast, "expected_direction": "+" if direction > 0 else "-",
            "n_patients": len(values), "direction_supported_patients": positive,
            "median_difference_z": float(np.median(values)), "median_direction_aligned_difference_z": float(np.median(aligned)),
            "sign_test_p": binomtest(positive, len(values), 0.5, alternative="two-sided").pvalue,
        })
    result = pd.DataFrame(rows)
    result["sign_test_fdr_bh_global"] = bh(result["sign_test_p"].to_numpy())
    return result


def main() -> None:
    if OUT.exists():
        raise SystemExit(f"Refusing to overwrite existing output: {OUT}")
    tables = OUT / "tables"
    logs = OUT / "logs"
    tables.mkdir(parents=True)
    logs.mkdir(parents=True)
    started = now()
    memory_rows: list[dict] = []
    fingerprint_rows: list[dict] = []
    gene_qc_rows: list[dict] = []
    spot_frames: list[pd.DataFrame] = []
    association_rows: list[dict] = []
    profile_rows: list[dict] = []
    contrast_rows: list[dict] = []

    memory_gate(memory_rows, "initial")
    for required in [MANIFEST, STATE_SCORES, GLMF_SPOT_SCORES, *GLMF_H5AD.values(), *GLMF_BOUNDARY.values()]:
        if not required.is_file():
            raise FileNotFoundError(required)
    manifest = pd.read_csv(MANIFEST, sep="\t").sort_values("sample").reset_index(drop=True)
    if len(manifest) != 12:
        raise ValueError(f"Expected 12 TJ-RCC cases, found {len(manifest)}")
    state_scores = pd.read_csv(STATE_SCORES, sep="\t")
    glmf_scores = pd.read_csv(GLMF_SPOT_SCORES, sep="\t")

    sample_counter = 0
    for sample, h5ad_path in GLMF_H5AD.items():
        memory_gate(memory_rows, "sample_analysis", sample)
        actual_hash = sha256(h5ad_path)
        fingerprint_rows.append({"sample": sample, "cohort": "GLMF", "path": str(h5ad_path), "size_bytes": h5ad_path.stat().st_size, "sha256": actual_hash})
        adata = ad.read_h5ad(h5ad_path)
        source = glmf_scores[glmf_scores["sample"].eq(sample)].copy()
        frame, qc = make_common_frame(sample, adata, source, "glmf", GLMF_BOUNDARY[sample])
        gene_qc_rows.extend(qc)
        spot_frames.append(frame)
        association_rows.extend(association_tests(frame, sample_counter))
        profiles, contrasts = boundary_tables(frame)
        profile_rows.extend(profiles); contrast_rows.extend(contrasts)
        sample_counter += 1
        del adata

    for record in manifest.itertuples(index=False):
        sample = str(record.sample)
        memory_gate(memory_rows, "sample_analysis", sample)
        input_path = Path(record.path)
        if input_path.stat().st_size != int(record.size_bytes) or sha256(input_path) != str(record.sha256):
            raise RuntimeError(f"Input fingerprint mismatch: {sample}")
        fingerprint_rows.append({"sample": sample, "cohort": "TJ-RCC", "path": str(input_path), "size_bytes": input_path.stat().st_size, "sha256": str(record.sha256)})
        boundary_path = TJ_ROOT / "advanced_full_v2/bayesspace" / sample / f"{sample}_computational_boundary_spots.tsv"
        if not boundary_path.is_file():
            raise FileNotFoundError(boundary_path)
        adata = ad.read_h5ad(input_path)
        source = state_scores[state_scores["sample"].eq(sample)].copy()
        frame, qc = make_common_frame(sample, adata, source, "tjrcc", boundary_path)
        gene_qc_rows.extend(qc)
        spot_frames.append(frame)
        association_rows.extend(association_tests(frame, sample_counter))
        profiles, contrasts = boundary_tables(frame)
        profile_rows.extend(profiles); contrast_rows.extend(contrasts)
        sample_counter += 1
        del adata

    memory_gate(memory_rows, "cohort_synthesis")
    spots = pd.concat(spot_frames, ignore_index=True)
    associations = pd.DataFrame(association_rows)
    associations["p_fdr_bh_within_sample"] = associations.groupby("sample")["p_value"].transform(lambda values: bh(values.to_numpy()))
    profiles = pd.DataFrame(profile_rows)
    contrasts = pd.DataFrame(contrast_rows)
    patient_tests = associations[associations["cohort"].eq("TJ-RCC")].copy()
    cohort_summary = cohort_direction_summary(patient_tests)
    boundary_patient = contrasts[contrasts["cohort"].eq("TJ-RCC")].copy()
    boundary_summary = boundary_cohort_summary(boundary_patient)

    spots.to_csv(tables / "spot_mechanism_scores_v1.tsv.gz", sep="\t", index=False, compression="gzip")
    associations.to_csv(tables / "sample_spatial_block_tests_v1.tsv", sep="\t", index=False)
    cohort_summary.to_csv(tables / "tjrcc_cohort_direction_summary_v1.tsv", sep="\t", index=False)
    profiles.to_csv(tables / "boundary_distance_profiles_v1.tsv", sep="\t", index=False)
    contrasts.to_csv(tables / "sample_interface_contrasts_v1.tsv", sep="\t", index=False)
    boundary_summary.to_csv(tables / "tjrcc_interface_direction_summary_v1.tsv", sep="\t", index=False)
    pd.DataFrame(gene_qc_rows).to_csv(tables / "signature_gene_qc_v1.tsv", sep="\t", index=False)
    pd.DataFrame(fingerprint_rows).to_csv(tables / "input_fingerprints_v1.tsv", sep="\t", index=False)
    pd.DataFrame(memory_rows).to_csv(tables / "memory_gates_v1.tsv", sep="\t", index=False)
    pd.DataFrame([{ "available_patients": 12, "patients_with_treatment_response": 0, "patients_with_survival": 0,
                    "status": "not_estimable", "reason": "No treatment-response, PFS, OS, or longitudinal therapy labels in the frozen TJ-RCC manifest." }]).to_csv(
        tables / "clinical_resistance_label_audit_v1.tsv", sep="\t", index=False
    )

    primary = cohort_summary[
        ((cohort_summary["route"].eq("route_independent")) | (cohort_summary["route"].eq("reference_light")))
        & cohort_summary["grid"].eq("5x5")
    ].copy()
    primary.to_csv(tables / "primary_5x5_mechanism_summary_v1.tsv", sep="\t", index=False)
    metadata = {
        "status": "completed",
        "started_at": started,
        "completed_at": now(),
        "samples": {"GLMF": 2, "TJ-RCC": 12},
        "spots": {key: int(value) for key, value in spots.groupby("cohort").size().items()},
        "permutations_per_spatial_test": N_PERM,
        "block_grids": ["4x4", "5x5"],
        "n_relations": len(RELATIONS),
        "signatures": SIGNATURES,
        "definitions": {
            "vascular_compartment": "geometric mean of within-sample Endothelial and SMC/Pericyte q05 abundance percentiles",
            "spp1_vascular_niche": "geometric mean of within-sample SPP1+ TAM, Endothelial, and SMC/Pericyte q05 abundance percentiles",
            "immune_exclusion": "negative association with T-cell balance and cytotoxicity residual adjusted for T-cell abundance, immune abundance, sequencing depth, detected genes, and spatial polynomial",
            "invasion_front": "pre-existing BayesSpace computational boundary; interface is Inner plus Outer interface",
            "monocyte_transition_boundary": "SPP1-monocyte balance is cross-sectional state evidence and does not prove lineage conversion",
        },
        "clinical_resistance": "not estimable: frozen data contain no response/survival/longitudinal therapy labels",
        "inference_boundary": "Spatial and cross-sectional evidence cannot independently prove 5-prime-MTA release, monocyte lineage conversion, or therapy-resistance causality.",
    }
    (OUT / "run_metadata_v1.json").write_text(json.dumps(metadata, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    (OUT / "RUN_COMPLETE").write_text(now() + "\n", encoding="utf-8")
    print(json.dumps(metadata, ensure_ascii=False, indent=2), flush=True)


if __name__ == "__main__":
    main()
