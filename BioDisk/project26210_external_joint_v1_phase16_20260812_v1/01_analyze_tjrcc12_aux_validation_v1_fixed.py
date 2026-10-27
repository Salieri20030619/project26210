#!/usr/bin/env python3
"""Pre-registered section-level validation in twelve TJRCC spatial sections."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import platform
from datetime import datetime
from pathlib import Path

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
from scipy.spatial import cKDTree
from scipy.stats import binomtest, rankdata, spearmanr, wilcoxon


SAMPLES = ["H46_T", "J38_T", "R114_T", "R29_T", "R51_T", "S15_T", "X49_T", "X98_T", "Y12_T", "Y27_T", "Y7_T", "Z43_T"]
EXPECTED_SPOTS = 34609
EXPECTED_GENES = 8830
EXPECTED_FACTORS = 14
SEED = 20260812
BOOTSTRAPS = 9999
ABUNDANCE_KEY = "q05_cell_abundance_w_sf"
ABUNDANCE_PREFIX = "q05cell_abundance_w_sf_"
SIGNATURES = [
    "f1_top100_weighted", "c1qc_state", "spp1_state_no_spp1", "spp1_state_full",
    "hypoxia_angiogenesis", "tcell_dysfunction", "cytotoxicity",
]
PRIMARY = [
    "F1polar_to_C1QC_expression",
    "F1polar_to_hypoxia_expression",
    "SPP1polar_to_Tdysfunction_expression",
]


def now_iso() -> str:
    return datetime.now().astimezone().isoformat(timespec="seconds")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


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
    x, y = np.asarray(x, float), np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 3 or np.ptp(x[ok]) == 0 or np.ptp(y[ok]) == 0:
        return math.nan
    return float(spearmanr(x[ok], y[ok]).statistic)


def residualize(y: np.ndarray, covariates: list[np.ndarray]) -> np.ndarray:
    design = np.column_stack([np.ones(len(y)), *[rankdata(np.asarray(x, float)) for x in covariates]])
    target = rankdata(np.asarray(y, float))
    coefficient, *_ = np.linalg.lstsq(design, target, rcond=None)
    return target - design @ coefficient


def bh(values: pd.Series) -> np.ndarray:
    p = values.to_numpy(float)
    order = np.argsort(p)
    ranked = p[order]
    adjusted = np.minimum.accumulate((ranked * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    out = np.empty(len(p), float)
    out[order] = np.minimum(adjusted, 1.0)
    return out


def locate_h5ad(input_root: Path, sample: str) -> Path:
    hits = sorted((input_root / sample).glob("*_reference_v4_cell2location.h5ad"))
    if len(hits) != 1:
        raise RuntimeError(f"Expected one final H5AD for {sample}, found {len(hits)}")
    return hits[0]


def validate_and_load(path: Path, status_path: Path, signature: pd.DataFrame, sample: str) -> tuple[pd.DataFrame, list[dict], dict]:
    status = json.loads(status_path.read_text(encoding="utf-8"))
    if status.get("status") != "completed" or int(status.get("epoch", -1)) != 10000:
        raise RuntimeError(f"Non-terminal sample status: {sample}")
    if int(status.get("output_size_bytes", -1)) != path.stat().st_size or status.get("output_sha256") != sha256(path):
        raise RuntimeError(f"Final H5AD hash/size mismatch: {sample}")

    adata = ad.read_h5ad(path, backed="r")
    try:
        if adata.n_vars != EXPECTED_GENES:
            raise RuntimeError(f"Unexpected gene count: {sample} {adata.n_vars}")
        genes = pd.Index(adata.var_names.astype(str))
        gene_lookup = {gene: i for i, gene in enumerate(genes)}
        wanted = signature[signature["signature_id"].isin(SIGNATURES)]["gene"].drop_duplicates().tolist()
        present = [gene for gene in wanted if gene in gene_lookup]
        indices = [gene_lookup[gene] for gene in present]
        selected = adata.X[:, indices]
        if hasattr(selected, "to_memory"):
            selected = selected.to_memory()
        if hasattr(selected, "toarray"):
            selected = selected.toarray()
        counts = np.asarray(selected, dtype=float)
        # Backed sparse datasets expose slicing but not matrix reduction methods.
        # Materialize one sample's sparse count matrix (not a dense array) before
        # computing library sizes; samples are processed serially.
        count_matrix = adata.X[:]
        totals_raw = count_matrix.sum(axis=1)
        if hasattr(totals_raw, "compute"):
            totals_raw = totals_raw.compute()
        totals = np.asarray(totals_raw, dtype=float).ravel()
        if np.any(~np.isfinite(totals)) or np.any(totals <= 0):
            raise RuntimeError(f"Invalid spatial library sizes: {sample}")
        lognorm = np.log1p(counts * (10000.0 / totals[:, None]))
        standard_deviation = lognorm.std(axis=0, ddof=0)
        z = np.zeros_like(lognorm)
        nonconstant = standard_deviation > 0
        z[:, nonconstant] = (lognorm[:, nonconstant] - lognorm[:, nonconstant].mean(axis=0)) / standard_deviation[nonconstant]
        gene_position = {gene: i for i, gene in enumerate(present)}

        abundance = adata.obsm[ABUNDANCE_KEY]
        abundance_values = np.asarray(abundance, float)
        states = [clean_state(value) for value in getattr(abundance, "columns", [])]
        if len(states) != abundance_values.shape[1]:
            states = [str(value) for value in adata.uns["mod"]["factor_names"]]
        if abundance_values.shape != (adata.n_obs, EXPECTED_FACTORS):
            raise RuntimeError(f"Unexpected abundance shape: {sample} {abundance_values.shape}")
        if np.any(~np.isfinite(abundance_values)) or np.any(abundance_values < 0):
            raise RuntimeError(f"Invalid abundance values: {sample}")
        abund = pd.DataFrame(abundance_values, columns=states)
        coordinates = np.asarray(adata.obsm["spatial"], float)[:, :2]
        spot_ids = adata.obs_names.astype(str).tolist()
    finally:
        adata.file.close()

    frame = pd.DataFrame({"sample": sample, "spot_id": spot_ids, "x": coordinates[:, 0], "y": coordinates[:, 1], "total_counts": totals})
    overlap_rows = []
    for signature_id in SIGNATURES:
        definition = signature[signature["signature_id"] == signature_id].sort_values("gene_order")
        present_genes = [gene for gene in definition["gene"] if gene in gene_position]
        missing = [gene for gene in definition["gene"] if gene not in gene_position]
        overlap_rows.append({"sample": sample, "signature_id": signature_id, "defined_genes": len(definition), "present_genes": len(present_genes), "missing_genes": ";".join(missing)})
        positions = [gene_position[gene] for gene in present_genes]
        if signature_id == "f1_top100_weighted":
            weights = definition.set_index("gene").loc[present_genes, "raw_weight"].to_numpy(float)
            frame[f"primary__{signature_id}"] = (lognorm[:, positions] @ weights) / weights.sum()
        else:
            frame[f"primary__{signature_id}"] = z[:, positions].mean(axis=1)
        frame[f"sensitivity__{signature_id}"] = lognorm[:, positions].mean(axis=1)

    required_states = ["F1-high Malignant", "F1-low Malignant", "C1QC-state Macrophage", "SPP1-state Macrophage", "Other Macrophage", "T cell", "Endothelial"]
    if not set(required_states).issubset(abund.columns):
        raise RuntimeError(f"Missing reference states: {sample}")
    for state in required_states:
        frame[state] = abund[state].to_numpy(float)
    frame["malignant_total"] = frame["F1-high Malignant"] + frame["F1-low Malignant"]
    frame["macrophage_total"] = frame["C1QC-state Macrophage"] + frame["SPP1-state Macrophage"] + frame["Other Macrophage"]
    frame["f1_high_share"] = frame["F1-high Malignant"] / frame["malignant_total"]
    frame["c1qc_share"] = frame["C1QC-state Macrophage"] / frame["macrophage_total"]
    frame["spp1_share"] = frame["SPP1-state Macrophage"] / frame["macrophage_total"]
    provenance = {"sample": sample, "path": str(path), "size_bytes": path.stat().st_size, "sha256": status["output_sha256"], "n_spots": len(frame), "n_genes": EXPECTED_GENES, "n_factors": EXPECTED_FACTORS}
    return frame, overlap_rows, provenance


def classify(row: pd.Series) -> str:
    direction = 1 if row["median_rho"] > 0 else -1 if row["median_rho"] < 0 else 0
    direction_count = row["positive_sections"] if direction > 0 else row["negative_sections"]
    reproduced = row["wilcoxon_fdr_bh"] < 0.05 and direction_count >= 8 and abs(row["median_rho"]) >= 0.10 and row["loo_direction_fraction"] == 1.0
    if reproduced:
        return "AUXILIARY_REPRODUCED"
    if row["positive_sections"] >= 8 or row["negative_sections"] >= 8:
        return "DIRECTIONAL_NOT_FDR"
    if row["positive_sections"] >= 4 and row["negative_sections"] >= 4:
        return "HETEROGENEOUS"
    return "NOT_SUPPORTED"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-root", type=Path, required=True)
    parser.add_argument("--signature-genes", type=Path, required=True)
    parser.add_argument("--registry", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    args.output.mkdir(parents=True)
    tables = args.output / "tables"; figures = args.output / "figures"
    tables.mkdir(); figures.mkdir()
    signature = pd.read_csv(args.signature_genes, sep="\t")

    frames, overlaps, provenance = [], [], []
    correlations = []
    for sample in SAMPLES:
        path = locate_h5ad(args.input_root, sample)
        frame, overlap, source = validate_and_load(path, args.input_root / sample / "status.json", signature, sample)
        frames.append(frame); overlaps.extend(overlap); provenance.append(source)
        thresholds = {
            "malignant_q25": frame["malignant_total"].quantile(0.25), "malignant_q50": frame["malignant_total"].quantile(0.50),
            "macrophage_q25": frame["macrophage_total"].quantile(0.25), "macrophage_q50": frame["macrophage_total"].quantile(0.50),
            "t_q25": frame["T cell"].quantile(0.25), "t_q50": frame["T cell"].quantile(0.50),
        }
        coordinates = frame[["x", "y"]].to_numpy(float)
        neighbor_indices = {k: cKDTree(coordinates).query(coordinates, k=k + 1)[1][:, 1:] for k in [6, 12]}
        hypotheses = [
            ("F1polar_to_C1QC_expression", "f1_high_share", "c1qc_state", ["macrophage_total", "total_counts"], "both"),
            ("F1polar_to_hypoxia_expression", "f1_high_share", "hypoxia_angiogenesis", ["malignant_total", "Endothelial", "total_counts"], "malignant"),
            ("SPP1polar_to_Tdysfunction_expression", "spp1_share", "tcell_dysfunction", ["T cell", "total_counts"], "macro_t"),
            ("SPP1polar_to_cytotoxicity_expression", "spp1_share", "cytotoxicity", ["T cell", "total_counts"], "macro_t"),
            ("C1QC_mapping_positive_control", "c1qc_share", "c1qc_state", ["macrophage_total", "total_counts"], "macrophage"),
            ("SPP1_mapping_positive_control", "spp1_share", "spp1_state_full", ["macrophage_total", "total_counts"], "macrophage"),
        ]
        for relation, x_name, y_signature, covariate_names, scope in hypotheses:
            x = frame[x_name].to_numpy(float)
            for score_method in ["primary", "sensitivity"]:
                y = frame[f"{score_method}__{y_signature}"].to_numpy(float)
                residual = residualize(y, [np.log1p(frame[name].to_numpy(float)) for name in covariate_names])
                for gate in ["q25", "q50"]:
                    if scope == "both":
                        mask = (frame.malignant_total >= thresholds[f"malignant_{gate}"]) & (frame.macrophage_total >= thresholds[f"macrophage_{gate}"])
                    elif scope == "malignant":
                        mask = frame.malignant_total >= thresholds[f"malignant_{gate}"]
                    elif scope == "macro_t":
                        mask = (frame.macrophage_total >= thresholds[f"macrophage_{gate}"]) & (frame["T cell"] >= thresholds[f"t_{gate}"])
                    else:
                        mask = frame.macrophage_total >= thresholds[f"macrophage_{gate}"]
                    for outcome, vector in [("raw_score", y), ("lineage_depth_residual", residual)]:
                        correlations.append({"sample": sample, "relation": relation, "score_method": score_method, "outcome": outcome, "analysis": f"same_spot_{gate}", "k_neighbors": np.nan, "n_spots": int(mask.sum()), "spearman_rho": safe_rho(x[mask], vector[mask])})
                        for k, neighbors in neighbor_indices.items():
                            lag = vector[neighbors].mean(axis=1)
                            correlations.append({"sample": sample, "relation": relation, "score_method": score_method, "outcome": outcome, "analysis": f"neighbor_lag_{gate}", "k_neighbors": k, "n_spots": int(mask.sum()), "spearman_rho": safe_rho(x[mask], lag[mask])})

        direct = [
            ("F1polar_to_C1QCpolar", "f1_high_share", "c1qc_share", "both"),
            ("F1polar_to_SPP1polar", "f1_high_share", "spp1_share", "both"),
            ("SPP1polar_to_T_abundance", "spp1_share", "T cell", "macro_t"),
        ]
        for relation, x_name, y_name, scope in direct:
            for gate in ["q25", "q50"]:
                if scope == "both":
                    mask = (frame.malignant_total >= thresholds[f"malignant_{gate}"]) & (frame.macrophage_total >= thresholds[f"macrophage_{gate}"])
                else:
                    mask = (frame.macrophage_total >= thresholds[f"macrophage_{gate}"]) & (frame["T cell"] >= thresholds[f"t_{gate}"])
                correlations.append({"sample": sample, "relation": relation, "score_method": "component", "outcome": "component", "analysis": f"same_spot_{gate}", "k_neighbors": np.nan, "n_spots": int(mask.sum()), "spearman_rho": safe_rho(frame.loc[mask, x_name], frame.loc[mask, y_name])})

    scores = pd.concat(frames, ignore_index=True)
    correlation_table = pd.DataFrame(correlations)
    primary = correlation_table[
        correlation_table["relation"].isin(PRIMARY)
        & (correlation_table["score_method"] == "primary")
        & (correlation_table["outcome"] == "lineage_depth_residual")
        & (correlation_table["analysis"] == "same_spot_q25")
    ].copy()
    if len(primary) != 36 or primary.groupby("relation")["sample"].nunique().ne(12).any():
        raise RuntimeError("Primary endpoint conservation failure")

    rng = np.random.default_rng(SEED)
    aggregate_rows, loo_rows = [], []
    for relation in PRIMARY:
        values = primary.loc[primary.relation == relation, "spearman_rho"].to_numpy(float)
        if np.any(~np.isfinite(values)):
            raise RuntimeError(f"Non-finite primary rho: {relation}")
        boot = np.median(values[rng.integers(0, len(values), size=(BOOTSTRAPS, len(values)))], axis=1)
        for sample in SAMPLES:
            kept = primary[(primary.relation == relation) & (primary["sample"] != sample)]["spearman_rho"].to_numpy(float)
            loo_rows.append({"relation": relation, "left_out_sample": sample, "median_rho": float(np.median(kept))})
        nonzero = values[values != 0]
        sign_p = float(binomtest(int((nonzero > 0).sum()), len(nonzero), 0.5, alternative="two-sided").pvalue) if len(nonzero) else 1.0
        try:
            wilcoxon_p = float(wilcoxon(values, alternative="two-sided", zero_method="wilcox", method="auto").pvalue)
        except ValueError:
            wilcoxon_p = 1.0
        aggregate_rows.append({
            "relation": relation, "sections": len(values), "median_rho": float(np.median(values)),
            "q25_rho": float(np.quantile(values, 0.25)), "q75_rho": float(np.quantile(values, 0.75)),
            "min_rho": float(values.min()), "max_rho": float(values.max()),
            "positive_sections": int((values > 0).sum()), "negative_sections": int((values < 0).sum()), "zero_sections": int((values == 0).sum()),
            "sign_test_p_two_sided": sign_p, "wilcoxon_p_two_sided": wilcoxon_p,
            "bootstrap_median_ci_low": float(np.quantile(boot, 0.025)), "bootstrap_median_ci_high": float(np.quantile(boot, 0.975)),
        })
    aggregate = pd.DataFrame(aggregate_rows)
    aggregate["wilcoxon_fdr_bh"] = bh(aggregate["wilcoxon_p_two_sided"])
    loo = pd.DataFrame(loo_rows)
    loo_fraction = []
    for row in aggregate.itertuples(index=False):
        values = loo.loc[loo.relation == row.relation, "median_rho"].to_numpy(float)
        direction = np.sign(row.median_rho)
        loo_fraction.append(float((np.sign(values) == direction).mean()))
    aggregate["loo_direction_fraction"] = loo_fraction
    aggregate["decision"] = aggregate.apply(classify, axis=1)

    atomic_tsv(pd.DataFrame(provenance), tables / "input_provenance_v1.tsv")
    atomic_tsv(pd.DataFrame(overlaps), tables / "signature_gene_overlap_v1.tsv")
    atomic_tsv(scores, tables / "spot_program_scores_v1.tsv.gz")
    atomic_tsv(correlation_table, tables / "section_correlations_all_v1.tsv")
    atomic_tsv(primary, tables / "primary_section_statistics_v1.tsv")
    atomic_tsv(aggregate, tables / "primary_across_section_inference_v1.tsv")
    atomic_tsv(loo, tables / "primary_leave_one_section_out_v1.tsv")

    fig, axes = plt.subplots(1, 3, figsize=(15, 5), sharey=True)
    labels = {PRIMARY[0]: "F1-C1QC", PRIMARY[1]: "F1-hypoxia", PRIMARY[2]: "SPP1-T dysfunction"}
    for axis, relation in zip(axes, PRIMARY):
        sub = primary[primary.relation == relation].set_index("sample").loc[SAMPLES].reset_index()
        colors = np.where(sub.spearman_rho >= 0, "#b2182b", "#2166ac")
        axis.bar(np.arange(len(sub)), sub.spearman_rho, color=colors)
        axis.axhline(0, color="black", linewidth=0.8)
        axis.set_xticks(np.arange(len(sub))); axis.set_xticklabels(sub["sample"], rotation=70, fontsize=8)
        axis.set_title(labels[relation]); axis.set_ylabel("section Spearman rho")
    fig.suptitle("TJRCC12 pre-registered primary spatial endpoints")
    part = figures / "primary_section_rho_v1.png.part"
    fig.savefig(part, dpi=180, bbox_inches="tight", format="png"); plt.close(fig)
    os.replace(part, figures / "primary_section_rho_v1.png")

    manifest_rows = []
    for path in sorted(p for p in args.output.rglob("*") if p.is_file() and not p.name.endswith(".part")):
        if path.name in {"OUTPUT_MANIFEST_v1.tsv", "PHASE16_TJRCC12_AUX_STATUS_v1.json"}:
            continue
        manifest_rows.append({"file": str(path.relative_to(args.output)), "size_bytes": path.stat().st_size, "sha256": sha256(path)})
    atomic_tsv(pd.DataFrame(manifest_rows), args.output / "OUTPUT_MANIFEST_v1.tsv")
    status = {
        "schema_version": "phase16_tjrcc12_aux_validation_v1", "status": "PASS_ANALYSIS_COMPLETE",
        "completed_at": now_iso(), "samples": SAMPLES, "biological_units": 12, "spots": int(len(scores)),
        "expected_spots": EXPECTED_SPOTS, "spot_as_replicate_inference": False,
        "primary_relations": PRIMARY, "primary_rows": int(len(primary)), "bootstrap_resamples": BOOTSTRAPS,
        "registry": str(args.registry), "registry_sha256": sha256(args.registry),
        "signature_genes": str(args.signature_genes), "signature_sha256": sha256(args.signature_genes),
        "python": platform.python_version(), "anndata": ad.__version__, "numpy": np.__version__, "pandas": pd.__version__, "scipy": scipy.__version__,
    }
    if len(scores) != EXPECTED_SPOTS:
        raise RuntimeError(f"Spot conservation failure: {len(scores)} != {EXPECTED_SPOTS}")
    atomic_json(status, args.output / "PHASE16_TJRCC12_AUX_STATUS_v1.json")
    print(json.dumps(status, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
