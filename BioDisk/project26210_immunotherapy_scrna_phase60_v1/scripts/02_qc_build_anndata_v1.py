#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gc
import gzip
import json
import os
import traceback
from datetime import datetime
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.io
import scipy.sparse as sp
from scipy.stats import median_abs_deviation


SEED = 2621060


def now() -> str:
    return datetime.now().astimezone().isoformat(timespec="seconds")


def log(message: str) -> None:
    print(f"[{now()}] {message}", flush=True)


def atomic_json(payload: dict, path: Path) -> None:
    part = Path(str(path) + ".part")
    part.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    os.replace(part, path)


def memory_gate(stage: str, minimum_gib: float = 80.0) -> None:
    with Path("/proc/meminfo").open() as handle:
        fields = {line.split(":", 1)[0]: line.split()[1] for line in handle}
    available = float(fields["MemAvailable"]) / 1024.0 / 1024.0
    log(f"MEMORY_GATE stage={stage} available_gib={available:.2f} minimum_gib={minimum_gib:.2f}")
    if available < minimum_gib:
        raise RuntimeError(f"Memory gate failed at {stage}: {available:.2f} GiB")


def robust_bounds(values: np.ndarray, n_mads: float, log1p: bool) -> tuple[float, float, float, float]:
    transformed = np.log1p(values) if log1p else values.astype(float)
    median = float(np.median(transformed))
    mad = float(median_abs_deviation(transformed, scale=1.0, nan_policy="omit"))
    if not np.isfinite(mad) or mad == 0:
        mad = float(np.finfo(float).eps)
    lower, upper = median - n_mads * mad, median + n_mads * mad
    if log1p:
        return float(np.expm1(lower)), float(np.expm1(upper)), median, mad
    return lower, upper, median, mad


def read_features(path: Path) -> pd.DataFrame:
    table = pd.read_csv(path, sep="\t", header=None, dtype=str, compression="gzip")
    if table.shape[1] < 2:
        raise RuntimeError(f"Invalid features file: {path}")
    table.columns = ["gene_id", "gene_symbol", "feature_type"][: table.shape[1]]
    if "feature_type" in table and table["feature_type"].eq("Gene Expression").any():
        table["keep_feature"] = table["feature_type"].eq("Gene Expression")
    else:
        table["keep_feature"] = True
    table["gene_symbol"] = table["gene_symbol"].fillna(table["gene_id"]).astype(str)
    return table


def aggregate_duplicate_symbols(counts: sp.csr_matrix, symbols: pd.Series) -> tuple[sp.csr_matrix, pd.Index]:
    codes, unique = pd.factorize(symbols.astype(str), sort=False)
    if len(unique) == len(symbols):
        return counts, pd.Index(unique, name="gene_symbol")
    mapper = sp.csr_matrix(
        (np.ones(len(codes), dtype=np.int8), (np.arange(len(codes)), codes)),
        shape=(len(codes), len(unique)),
    )
    return (counts @ mapper).tocsr(), pd.Index(unique, name="gene_symbol")


def read_library(row: pd.Series) -> ad.AnnData:
    with gzip.open(row["matrix"], "rb") as handle:
        matrix = scipy.io.mmread(handle)
    if not sp.issparse(matrix):
        matrix = sp.coo_matrix(matrix)
    counts = matrix.T.tocsr()
    counts.sort_indices()
    features = read_features(Path(row["features"]))
    if len(features) != counts.shape[1]:
        raise RuntimeError(f"Feature dimension mismatch for {row['library_id']}")
    keep = features["keep_feature"].to_numpy(bool)
    counts = counts[:, keep]
    features = features.loc[keep].reset_index(drop=True)
    counts, symbols = aggregate_duplicate_symbols(counts, features["gene_symbol"])
    barcodes = pd.read_csv(row["barcodes"], sep="\t", header=None, dtype=str, compression="gzip")[0]
    if len(barcodes) != counts.shape[0]:
        raise RuntimeError(f"Barcode dimension mismatch for {row['library_id']}")
    obs_names = pd.Index([f"{row['cohort']}__{row['gsm']}__{x}" for x in barcodes], name="cell_id")
    if obs_names.duplicated().any():
        raise RuntimeError(f"Duplicate cell IDs in {row['library_id']}")
    obs = pd.DataFrame(index=obs_names)
    for column in [
        "cohort", "gsm", "library_id", "patient_id", "timepoint", "region", "organ",
        "treatment", "response_binary", "response_detail", "primary_response_eligible",
    ]:
        obs[column] = row[column]
    return ad.AnnData(X=counts, obs=obs, var=pd.DataFrame(index=symbols))


def compute_qc(obj: ad.AnnData) -> tuple[pd.DataFrame, dict]:
    counts = obj.X.tocsr()
    total = np.asarray(counts.sum(axis=1)).ravel().astype(float)
    genes = np.diff(counts.indptr).astype(int)
    upper = obj.var_names.str.upper()
    mt_mask = np.asarray(upper.str.startswith("MT-"), dtype=bool)
    rb_mask = np.asarray(upper.str.match(r"^RP[SL][0-9]"), dtype=bool)
    hb_mask = np.asarray(upper.str.match(r"^HB[ABDEGQZ][0-9]"), dtype=bool)
    denominator = np.maximum(total, 1.0)
    mt = np.asarray(counts[:, mt_mask].sum(axis=1)).ravel() if mt_mask.any() else np.zeros(len(obj))
    rb = np.asarray(counts[:, rb_mask].sum(axis=1)).ravel() if rb_mask.any() else np.zeros(len(obj))
    hb = np.asarray(counts[:, hb_mask].sum(axis=1)).ravel() if hb_mask.any() else np.zeros(len(obj))
    qc = obj.obs.copy()
    qc["total_counts_raw"], qc["n_genes_raw"] = total, genes
    qc["pct_counts_mt_raw"] = 100.0 * mt / denominator
    qc["pct_counts_rb_raw"] = 100.0 * rb / denominator
    qc["pct_counts_hb_raw"] = 100.0 * hb / denominator
    count_low, count_high, count_med, count_mad = robust_bounds(total, 5.0, True)
    gene_low, gene_high, gene_med, gene_mad = robust_bounds(genes, 5.0, True)
    _, mt_high, mt_med, mt_mad = robust_bounds(qc["pct_counts_mt_raw"].to_numpy(float), 3.0, False)
    thresholds = {
        "count_lower_5mad": max(200.0, count_low), "count_upper_5mad_warning": count_high,
        "gene_lower_5mad": max(200.0, gene_low), "gene_upper_5mad_warning": gene_high,
        "mito_upper_3mad_and_floor10": max(10.0, mt_high),
        "log1p_count_median": count_med, "log1p_count_mad": count_mad,
        "log1p_gene_median": gene_med, "log1p_gene_mad": gene_mad,
        "mito_median": mt_med, "mito_mad": mt_mad,
    }
    qc["low_quality"] = (
        (total < thresholds["count_lower_5mad"]) | (genes < thresholds["gene_lower_5mad"])
        | (qc["pct_counts_mt_raw"].to_numpy(float) > thresholds["mito_upper_3mad_and_floor10"])
        | (total < 200) | (genes < 200) | (qc["pct_counts_mt_raw"].to_numpy(float) > 25)
    )
    qc["high_complexity_warning"] = (total > count_high) | (genes > gene_high)
    return qc, thresholds


def call_doublets(obj: ad.AnnData, seed: int) -> tuple[np.ndarray, np.ndarray, dict]:
    n_cells = obj.n_obs
    if n_cells < 100:
        return np.full(n_cells, np.nan), np.zeros(n_cells, bool), {
            "status": "not_estimable_small_library", "n_cells": n_cells, "error": ""
        }
    attempts = list(dict.fromkeys([min(20, max(5, n_cells // 25)), min(10, max(3, n_cells // 40))]))
    errors = []
    for n_pcs in attempts:
        work = ad.AnnData(X=obj.X.copy(), obs=pd.DataFrame(index=obj.obs_names), var=pd.DataFrame(index=obj.var_names))
        try:
            sc.pp.scrublet(work, expected_doublet_rate=0.05, sim_doublet_ratio=2.0,
                           n_prin_comps=n_pcs, use_approx_neighbors=True,
                           random_state=seed, verbose=False)
            score = work.obs["doublet_score"].to_numpy(float)
            predicted = work.obs["predicted_doublet"].to_numpy(bool)
            threshold = float(work.uns.get("scrublet", {}).get("threshold", np.nan))
            return score, predicted, {
                "status": "passed", "n_cells": n_cells, "n_prin_comps": n_pcs,
                "threshold": threshold, "n_predicted_doublets": int(predicted.sum()), "error": "",
            }
        except Exception as exc:
            errors.append(f"n_pcs={n_pcs}: {type(exc).__name__}: {exc}")
    return np.full(n_cells, np.nan), np.zeros(n_cells, bool), {
        "status": "not_estimable_after_retry", "n_cells": n_cells, "error": " | ".join(errors)
    }


def process_library(row: pd.Series, output: Path, number: int) -> dict:
    library = row["library_id"]
    checkpoint = output / "checkpoints" / f"{library}.json"
    h5ad_path = output / "libraries" / f"{library}.filtered_v1.h5ad"
    qc_path = output / "cell_qc" / f"{library}.cell_qc_v1.tsv.gz"
    if checkpoint.is_file() and h5ad_path.is_file() and qc_path.is_file():
        payload = json.loads(checkpoint.read_text(encoding="utf-8"))
        if payload.get("status") == "PASS":
            log(f"RESUME_SKIP {number:02d} {library}")
            return payload
    memory_gate(f"read_{library}")
    started = datetime.now().astimezone()
    obj = read_library(row)
    if obj.X.nnz == 0 or np.any(obj.X.data < 0) or not np.all(np.isfinite(obj.X.data)):
        raise RuntimeError(f"Invalid counts for {library}")
    qc, thresholds = compute_qc(obj)
    scores, predicted, doublet = call_doublets(obj, SEED + number)
    qc["scrublet_score"], qc["scrublet_predicted_doublet"] = scores, predicted
    qc["doublet_status"] = doublet["status"]
    qc["analysis_keep"] = (~qc["low_quality"]) & (~qc["scrublet_predicted_doublet"])
    qc["qc_status"] = np.select([qc["low_quality"], qc["scrublet_predicted_doublet"]],
                                 ["Low quality", "Doublet"], default="Retained")
    kept = obj[qc["analysis_keep"].to_numpy(bool)].copy()
    kept.obs = qc.loc[kept.obs_names].copy()
    kept.layers["counts"] = kept.X.copy()
    kept.uns["qc_thresholds_v1"], kept.uns["doublet_qc_v1"] = thresholds, doublet
    part_h5ad = Path(str(h5ad_path) + ".part")
    kept.write_h5ad(part_h5ad, compression="gzip")
    os.replace(part_h5ad, h5ad_path)
    part_qc = Path(str(qc_path) + ".part")
    qc.to_csv(part_qc, sep="\t", compression="gzip", index=True)
    os.replace(part_qc, qc_path)
    payload = {
        "status": "PASS", "completed_at": now(), "cohort": row["cohort"],
        "library_id": library, "gsm": row["gsm"], "patient_id": row["patient_id"],
        "n_input": int(obj.n_obs), "n_genes": int(obj.n_vars), "nnz_input": int(obj.X.nnz),
        "n_retained": int(kept.n_obs), "n_low_quality": int(qc["low_quality"].sum()),
        "n_doublet": int(qc["scrublet_predicted_doublet"].sum()),
        "retained_fraction": float(kept.n_obs / obj.n_obs),
        "doublet_status": doublet["status"], "doublet_error": doublet.get("error", ""),
        "thresholds": thresholds,
        "runtime_seconds": (datetime.now().astimezone() - started).total_seconds(),
        "h5ad": str(h5ad_path), "cell_qc": str(qc_path),
    }
    atomic_json(payload, checkpoint)
    log(f"LIBRARY_PASS {number:02d} {library} input={obj.n_obs} retained={kept.n_obs} "
        f"doublets={payload['n_doublet']} doublet_status={doublet['status']}")
    del obj, kept, qc
    gc.collect()
    return payload


def combine_cohort(cohort: str, manifest: pd.DataFrame, output: Path) -> dict:
    destination = output / "cohorts" / f"{cohort}.qc_filtered_raw_v1.h5ad"
    checkpoint = output / "cohorts" / f"{cohort}.qc_filtered_raw_v1.json"
    if destination.is_file() and checkpoint.is_file():
        payload = json.loads(checkpoint.read_text(encoding="utf-8"))
        if payload.get("status") == "PASS":
            log(f"COHORT_RESUME_SKIP {cohort}")
            return payload
    memory_gate(f"combine_{cohort}")
    objects, expected_cells = [], 0
    for library in manifest.loc[manifest["cohort"].eq(cohort), "library_id"]:
        item = ad.read_h5ad(output / "libraries" / f"{library}.filtered_v1.h5ad")
        expected_cells += item.n_obs
        objects.append(item)
    combined = ad.concat(objects, join="outer", merge="same", uns_merge="unique", index_unique=None)
    if combined.n_obs != expected_cells or combined.obs_names.duplicated().any():
        raise RuntimeError(f"Cohort cell conservation failure: {cohort}")
    if "counts" not in combined.layers:
        combined.layers["counts"] = combined.X.copy()
    combined.uns["phase60_scope_v1"] = {
        "cohort": cohort, "n_libraries": int((manifest["cohort"] == cohort).sum()),
        "patient_statistical_unit": True,
        "raw_cohort_matrices_not_merged_for_response_testing": True,
    }
    part = Path(str(destination) + ".part")
    combined.write_h5ad(part, compression="gzip")
    os.replace(part, destination)
    payload = {
        "status": "PASS", "completed_at": now(), "cohort": cohort,
        "n_cells": int(combined.n_obs), "n_genes": int(combined.n_vars), "nnz": int(combined.X.nnz),
        "n_libraries": int(combined.obs["library_id"].nunique()),
        "n_patients": int(combined.obs["patient_id"].nunique()), "h5ad": str(destination),
    }
    atomic_json(payload, checkpoint)
    log(f"COHORT_PASS {cohort} cells={combined.n_obs} genes={combined.n_vars}")
    del combined, objects
    gc.collect()
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    for subdir in ["checkpoints", "libraries", "cell_qc", "cohorts", "tables"]:
        (args.output / subdir).mkdir(parents=True, exist_ok=True)
    manifest = pd.read_csv(args.manifest, sep="\t", dtype=str)
    manifest["primary_response_eligible"] = manifest["primary_response_eligible"].str.lower().eq("true")
    results = [process_library(row, args.output, number)
               for number, (_, row) in enumerate(manifest.iterrows(), start=1)]
    summary = pd.DataFrame(results).sort_values(["cohort", "patient_id", "library_id"])
    summary.to_csv(args.output / "tables" / "library_qc_summary_v1.tsv", sep="\t", index=False)
    cohort_results = [combine_cohort(cohort, manifest, args.output) for cohort in ["GSE264586", "GSE294109"]]
    status = {
        "status": "PASS", "completed_at": now(), "n_libraries": len(results),
        "n_input": int(summary["n_input"].sum()), "n_retained": int(summary["n_retained"].sum()),
        "doublet_not_estimable_libraries": int((summary["doublet_status"] != "passed").sum()),
        "cohorts": cohort_results,
    }
    atomic_json(status, args.output / "QC_BUILD_COMPLETE_v1.json")
    log("QC_BUILD_PASS " + json.dumps(status, ensure_ascii=False))


if __name__ == "__main__":
    try:
        main()
    except Exception:
        traceback.print_exc()
        raise
