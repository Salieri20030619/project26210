#!/usr/bin/env python3
"""Bridge TYMP and rank-7 F1 in the Part 1 rebuilt frozen malignant cells."""

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
import numpy as np
import pandas as pd
import scipy
from scipy.sparse import issparse
from scipy.stats import spearmanr


EXPECTED_PRIMARY = 11950
MIN_PATIENT_CELLS = 20
MIN_QUADRANT_CELLS = 500
MIN_PATIENTS_10 = 6
MIN_STUDIES = 3
MAX_SIGNATURE_RHO = 0.98
MIN_MARKERS = 20
MIN_ABS_LOG2FC = 0.5


def now() -> str:
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


def safe_rho(x: np.ndarray, y: np.ndarray) -> float:
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 3 or np.ptp(x[ok]) == 0 or np.ptp(y[ok]) == 0:
        return math.nan
    return float(spearmanr(x[ok], y[ok]).statistic)


def clean_gene(gene: str) -> bool:
    upper = gene.upper()
    return not (upper.startswith("MT-") or upper.startswith("RPL") or upper.startswith("RPS"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", type=Path, required=True)
    parser.add_argument("--usage", type=Path, required=True)
    parser.add_argument("--reference-genes", type=Path, required=True)
    parser.add_argument("--registry", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    args.output.mkdir(parents=True)
    tables = args.output / "tables"; tables.mkdir()

    obj = ad.read_h5ad(args.h5ad, backed="r")
    required = ["eligible_for_subclone_selection_v1", "curated_patient_id", "study", "source_tissue", "total_counts_raw", "cell_type_broad_v1"]
    missing = sorted(set(required).difference(obj.obs.columns))
    if missing:
        raise RuntimeError(f"Missing H5AD metadata: {missing}")
    obs = obj.obs[required].copy()
    obs.index = obs.index.astype(str)
    eligible = obs["eligible_for_subclone_selection_v1"].astype(bool)
    if int(eligible.sum()) != EXPECTED_PRIMARY:
        raise RuntimeError(f"Primary malignant conservation failure: {eligible.sum()}")
    primary = obs.loc[eligible].copy()
    primary["exact_cell_id"] = primary.index
    primary["patient"] = primary["curated_patient_id"].astype(str)
    primary["study_v1"] = primary["study"].astype(str)

    usage = pd.read_csv(args.usage, sep="\t", usecols=["exact_cell_id", "F1"])
    if usage["exact_cell_id"].duplicated().any():
        raise RuntimeError("Duplicate F1 usage cell IDs")
    usage = usage.set_index("exact_cell_id")
    primary["F1"] = usage["F1"].reindex(primary.index).to_numpy(float)
    if primary["F1"].isna().any() or np.any(primary["F1"] < 0):
        raise RuntimeError("Incomplete or invalid F1 usage")

    if "TYMP" not in obj.var_names:
        raise RuntimeError("TYMP absent from Part 1 H5AD")
    gene_position = int(obj.var_names.get_loc("TYMP"))
    tymp_raw = obj.layers["counts"][np.flatnonzero(eligible.to_numpy()), gene_position]
    if hasattr(tymp_raw, "to_memory"):
        tymp_raw = tymp_raw.to_memory()
    if issparse(tymp_raw):
        tymp_raw = tymp_raw.toarray()
    primary["TYMP_raw_count"] = np.asarray(tymp_raw, float).ravel()
    library = pd.to_numeric(primary["total_counts_raw"], errors="coerce").to_numpy(float)
    if np.any(~np.isfinite(library)) or np.any(library <= 0):
        raise RuntimeError("Invalid total_counts_raw")
    primary["TYMP_log1p_CP10K"] = np.log1p(primary["TYMP_raw_count"].to_numpy(float) * 10000.0 / library)

    patient_counts = primary.groupby("patient", observed=True).size()
    eligible_patients = patient_counts[patient_counts >= MIN_PATIENT_CELLS].index
    primary["bridge_eligible_patient"] = primary["patient"].isin(eligible_patients)
    bridge = primary.loc[primary["bridge_eligible_patient"]].copy()
    thresholds = bridge.groupby("patient", observed=True).agg(
        n_cells=("exact_cell_id", "size"),
        f1_q75=("F1", lambda x: float(x.quantile(0.75))),
        tymp_q75=("TYMP_log1p_CP10K", lambda x: float(x.quantile(0.75))),
        tymp_zero_fraction=("TYMP_raw_count", lambda x: float((x == 0).mean())),
        tymp_mean=("TYMP_log1p_CP10K", "mean"),
        f1_mean=("F1", "mean"),
    ).reset_index()
    threshold_lookup = thresholds.set_index("patient")
    bridge["F1_q75_patient"] = bridge["patient"].map(threshold_lookup["f1_q75"])
    bridge["TYMP_q75_patient"] = bridge["patient"].map(threshold_lookup["tymp_q75"])
    bridge["F1_high"] = bridge["F1"] >= bridge["F1_q75_patient"]
    bridge["TYMP_high"] = (bridge["TYMP_log1p_CP10K"] >= bridge["TYMP_q75_patient"]) & (bridge["TYMP_raw_count"] > 0)
    bridge["quadrant"] = np.select(
        [bridge.F1_high & bridge.TYMP_high, bridge.F1_high & ~bridge.TYMP_high, ~bridge.F1_high & bridge.TYMP_high],
        ["F1high_TYMPhigh", "F1high_TYMPlow", "F1low_TYMPhigh"],
        default="F1low_TYMPlow",
    )

    correlation_rows = []
    for patient, frame in bridge.groupby("patient", observed=True):
        correlation_rows.append({"patient": patient, "n_cells": len(frame), "spearman_rho_F1_TYMP": safe_rho(frame.F1.to_numpy(float), frame.TYMP_log1p_CP10K.to_numpy(float)), "tymp_zero_fraction": float((frame.TYMP_raw_count == 0).mean())})
    correlations = pd.DataFrame(correlation_rows)

    patient_quadrants = bridge.groupby(["study_v1", "patient", "quadrant"], observed=True).agg(
        n_cells=("exact_cell_id", "size"), median_F1=("F1", "median"), median_TYMP=("TYMP_log1p_CP10K", "median")
    ).reset_index()
    all_quadrants = ["F1high_TYMPhigh", "F1high_TYMPlow", "F1low_TYMPhigh", "F1low_TYMPlow"]
    quadrant_summary_rows = []
    for quadrant in all_quadrants:
        subset = bridge[bridge.quadrant == quadrant]
        pq = patient_quadrants[patient_quadrants.quadrant == quadrant]
        quadrant_summary_rows.append({
            "quadrant": quadrant, "n_cells": len(subset), "n_patients": subset.patient.nunique(),
            "n_patients_ge10": int((pq.n_cells >= 10).sum()), "n_studies": subset.study_v1.nunique(),
            "median_TYMP": float(subset.TYMP_log1p_CP10K.median()), "median_F1": float(subset.F1.median()),
        })
    quadrant_summary = pd.DataFrame(quadrant_summary_rows)

    reference_table = pd.read_csv(args.reference_genes, sep="\t", index_col=0)
    genes = [gene for gene in reference_table.index.astype(str) if gene in obj.var_names]
    positions = [int(obj.var_names.get_loc(gene)) for gene in genes]
    counts = obj.layers["counts"][np.flatnonzero(eligible.to_numpy()), positions]
    if hasattr(counts, "to_memory"):
        counts = counts.to_memory()
    if not issparse(counts):
        counts = scipy.sparse.csr_matrix(np.asarray(counts))
    counts = counts.tocsr()
    bridge_positions = np.flatnonzero(primary["bridge_eligible_patient"].to_numpy())
    counts = counts[bridge_positions, :]
    if counts.shape[0] != len(bridge):
        raise RuntimeError("Bridge count alignment failure")

    signatures = {}
    for quadrant in all_quadrants:
        profiles = []
        for study, study_frame in bridge.groupby("study_v1", observed=True):
            study_profiles = []
            for patient, patient_frame in study_frame.groupby("patient", observed=True):
                rows = np.flatnonzero((bridge.patient.to_numpy() == patient) & (bridge.quadrant.to_numpy() == quadrant))
                if len(rows):
                    study_profiles.append(np.asarray(counts[rows].mean(axis=0)).ravel())
            if study_profiles:
                profiles.append(np.mean(study_profiles, axis=0))
        if not profiles:
            raise RuntimeError(f"No profiles for {quadrant}")
        signatures[quadrant] = np.mean(profiles, axis=0)
    signature = pd.DataFrame(signatures, index=genes)
    signature.index.name = "gene"

    correlation_matrix = signature.corr(method="spearman")
    signature_correlation_rows = []
    for i, left in enumerate(all_quadrants):
        for right in all_quadrants[i + 1:]:
            signature_correlation_rows.append({"left": left, "right": right, "spearman_rho": float(correlation_matrix.loc[left, right])})
    signature_correlations = pd.DataFrame(signature_correlation_rows)

    marker_rows = []
    clean = np.array([clean_gene(gene) for gene in genes])
    for quadrant in all_quadrants:
        other = signature.drop(columns=quadrant).mean(axis=1).to_numpy(float)
        target = signature[quadrant].to_numpy(float)
        log2fc = np.log2((target + 0.01) / (other + 0.01))
        order = np.argsort(-np.abs(log2fc))
        for rank, pos in enumerate(order[:200], 1):
            marker_rows.append({"quadrant": quadrant, "gene": genes[pos], "rank_abs_log2fc": rank, "log2fc_vs_other": float(log2fc[pos]), "clean_gene": bool(clean[pos])})
    markers = pd.DataFrame(marker_rows)
    marker_counts = markers[(markers.clean_gene) & (markers.log2fc_vs_other.abs() >= MIN_ABS_LOG2FC)].groupby("quadrant").size().reindex(all_quadrants, fill_value=0)
    quadrant_summary["markers_abs_log2fc_ge0_5"] = quadrant_summary.quadrant.map(marker_counts).astype(int)

    gates = pd.DataFrame([
        {"gate": "quadrant_cells", "result": "PASS" if (quadrant_summary.n_cells >= MIN_QUADRANT_CELLS).all() else "FAIL", "detail": quadrant_summary.set_index("quadrant").n_cells.to_dict()},
        {"gate": "quadrant_patient_coverage", "result": "PASS" if (quadrant_summary.n_patients_ge10 >= MIN_PATIENTS_10).all() else "FAIL", "detail": quadrant_summary.set_index("quadrant").n_patients_ge10.to_dict()},
        {"gate": "quadrant_study_coverage", "result": "PASS" if (quadrant_summary.n_studies >= MIN_STUDIES).all() else "FAIL", "detail": quadrant_summary.set_index("quadrant").n_studies.to_dict()},
        {"gate": "signature_pairwise_separability", "result": "PASS" if (signature_correlations.spearman_rho < MAX_SIGNATURE_RHO).all() else "FAIL", "detail": {"max_pairwise_rho": float(signature_correlations.spearman_rho.max())}},
        {"gate": "quadrant_marker_separability", "result": "PASS" if (quadrant_summary.markers_abs_log2fc_ge0_5 >= MIN_MARKERS).all() else "FAIL", "detail": quadrant_summary.set_index("quadrant").markers_abs_log2fc_ge0_5.to_dict()},
        {"gate": "positive_tymp_q75", "result": "PASS" if float((thresholds.tymp_q75 > 0).mean()) >= 0.80 else "FAIL", "detail": {"positive_fraction": float((thresholds.tymp_q75 > 0).mean()), "patients": len(thresholds)}},
    ])
    four_state_pass = bool((gates.result == "PASS").all())
    recommendation = "FOUR_QUADRANT_REFERENCE_ALLOWED" if four_state_pass else "TYMP_FOCUSED_LAYER_WITH_CONTINUOUS_F1"

    atomic_tsv(bridge.reset_index(drop=True), tables / "tymp_f1_cell_bridge_v1.tsv.gz")
    atomic_tsv(thresholds, tables / "patient_thresholds_v1.tsv")
    atomic_tsv(correlations, tables / "patient_f1_tymp_correlations_v1.tsv")
    atomic_tsv(patient_quadrants, tables / "patient_quadrant_counts_v1.tsv")
    atomic_tsv(quadrant_summary, tables / "quadrant_summary_v1.tsv")
    atomic_tsv(signature.reset_index(), tables / "quadrant_study_balanced_signatures_v1.tsv.gz")
    atomic_tsv(signature_correlations, tables / "quadrant_signature_correlations_v1.tsv")
    atomic_tsv(markers, tables / "quadrant_marker_candidates_v1.tsv")
    gates["detail"] = gates.detail.map(lambda value: json.dumps(value, ensure_ascii=False, sort_keys=True))
    atomic_tsv(gates, tables / "four_quadrant_feasibility_gates_v1.tsv")

    status = {
        "schema_version": "phase17_tymp_f1_bridge_v1", "status": "PASS_ANALYSIS_COMPLETE", "completed_at": now(),
        "primary_malignant_cells": int(len(primary)), "bridge_cells": int(len(bridge)), "eligible_patients": int(len(thresholds)),
        "median_patient_rho_F1_TYMP": float(correlations.spearman_rho_F1_TYMP.median()),
        "positive_patient_fraction_F1_TYMP": float((correlations.spearman_rho_F1_TYMP > 0).mean()),
        "four_quadrant_gate_pass": four_state_pass, "recommendation": recommendation,
        "registry": str(args.registry), "registry_sha256": sha256(args.registry),
        "inputs": {"h5ad": str(args.h5ad), "usage": str(args.usage), "reference_genes": str(args.reference_genes)},
        "python": platform.python_version(), "anndata": ad.__version__, "numpy": np.__version__, "pandas": pd.__version__, "scipy": scipy.__version__,
    }
    atomic_json(status, args.output / "PHASE17_TYMP_F1_BRIDGE_STATUS_v1.json")
    obj.file.close()
    print(json.dumps(status, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
