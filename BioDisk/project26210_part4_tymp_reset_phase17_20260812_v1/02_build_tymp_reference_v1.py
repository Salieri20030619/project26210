#!/usr/bin/env python3
"""Build the Phase 17 TYMP-focused 14-state spatial reference.

The 12 non-malignant states are rebuilt with the frozen reference-v3 rules.
Only the malignant layer changes: patients with a positive within-patient TYMP
q75 are split into TYMP-high and TYMP-low. Rank-7 F1 is retained as a
continuous annotation and is never used as a TYMP label.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import shutil
import subprocess
from datetime import datetime
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy
from scipy.sparse import csr_matrix, issparse
from scipy.stats import spearmanr


SPP1_PROGRAM = ["SPP1", "GPNMB", "LGALS3", "FABP5", "CD9", "LPL"]
C1QC_PROGRAM = ["C1QA", "C1QB", "C1QC", "APOE"]
BROAD_MAP = {
    "B cell": "B cell",
    "Dendritic cell": "Dendritic cell",
    "Endothelial": "Endothelial",
    "Mast cell": "Mast cell",
    "Monocyte": "Monocyte",
    "NK cell": "NK cell",
    "Pericyte/SMC": "SMC/Pericyte",
    "T cell": "T cell",
}
EXPECTED_STATES = sorted(
    list(BROAD_MAP.values())
    + [
        "Epithelial-Normal",
        "TYMP-high Malignant",
        "TYMP-low Malignant",
        "SPP1-state Macrophage",
        "C1QC-state Macrophage",
        "Other Macrophage",
    ]
)
OLD_MALIGNANT_STATES = {"F1-high Malignant", "F1-low Malignant"}
MIN_MALIGNANT_PER_PATIENT = 20
MIN_POSITIVE_TYMP_PATIENT_FRACTION = 0.80
TOP_MARKERS = 200


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", type=Path, required=True)
    parser.add_argument("--usage", type=Path, required=True)
    parser.add_argument("--bridge", type=Path, required=True)
    parser.add_argument("--thresholds", type=Path, required=True)
    parser.add_argument("--old-signature", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--memory-gate", type=Path)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_tsv(frame: pd.DataFrame, path: Path, *, index: bool = False) -> None:
    part = path.with_name(path.name + ".part")
    compression = "gzip" if path.name.endswith(".gz") else None
    frame.to_csv(
        part,
        sep="\t",
        index=index,
        encoding="utf-8",
        na_rep="NA",
        compression=compression,
        lineterminator="\n",
    )
    os.replace(part, path)


def atomic_json(value: dict[str, object], path: Path) -> None:
    part = path.with_name(path.name + ".part")
    part.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    os.replace(part, path)


def memory_gate(path: Path | None, stage: str, rows: list[dict[str, object]]) -> None:
    if path is None:
        rows.append({"stage": stage, "exit_code": 0, "output": "not_requested"})
        return
    result = subprocess.run([str(path)], capture_output=True, text=True, check=False)
    line = (result.stdout or result.stderr).strip().replace("\n", " | ")
    rows.append({"stage": stage, "exit_code": result.returncode, "output": line})
    if result.returncode != 0:
        raise RuntimeError(f"Memory gate failed at {stage}: {line}")


def program_score(
    counts: csr_matrix,
    positions: np.ndarray,
    gene_index: dict[str, int],
    genes: list[str],
    libraries: np.ndarray,
) -> np.ndarray:
    missing = sorted(set(genes).difference(gene_index))
    if missing:
        raise RuntimeError(f"Missing program genes: {missing}")
    matrix = counts[positions][:, [gene_index[g] for g in genes]].toarray().astype(float)
    return np.log1p(matrix / np.maximum(libraries[:, None], 1.0) * 1e4).mean(axis=1)


def main() -> None:
    args = parse_args()
    for path in [args.h5ad, args.usage, args.bridge, args.thresholds, args.old_signature]:
        if not path.is_file():
            raise FileNotFoundError(path)
    if args.output.exists():
        raise FileExistsError(f"Refusing to overwrite existing output: {args.output}")
    work = args.output.with_name(args.output.name + ".part")
    if work.exists():
        raise FileExistsError(f"Stale work directory exists; inspect it first: {work}")
    work.mkdir(parents=True)

    memory_rows: list[dict[str, object]] = []
    memory_gate(args.memory_gate, "before_load", memory_rows)

    old_signature = pd.read_csv(args.old_signature, sep="\t", index_col=0)
    old_signature.index = old_signature.index.astype(str)
    if old_signature.shape != (12758, 14):
        raise RuntimeError(f"Unexpected frozen reference shape: {old_signature.shape}")
    genes = old_signature.index.tolist()

    thresholds = pd.read_csv(args.thresholds, sep="\t")
    required_threshold_cols = {"patient", "n_cells", "f1_q75", "tymp_q75"}
    if not required_threshold_cols.issubset(thresholds.columns):
        raise RuntimeError("Patient threshold table is incomplete")
    thresholds["patient"] = thresholds["patient"].astype(str)
    if thresholds["patient"].duplicated().any() or len(thresholds) != 20:
        raise RuntimeError("Expected 20 unique bridge-eligible patients")
    positive_patients = set(thresholds.loc[thresholds["tymp_q75"] > 0, "patient"])
    if len(positive_patients) / len(thresholds) < MIN_POSITIVE_TYMP_PATIENT_FRACTION:
        raise RuntimeError("Too few patients have a positive TYMP q75")

    bridge = pd.read_csv(args.bridge, sep="\t", index_col="exact_cell_id")
    bridge.index = bridge.index.astype(str)
    required_bridge_cols = {
        "patient", "study_v1", "F1", "TYMP_raw_count", "TYMP_log1p_CP10K",
        "bridge_eligible_patient", "F1_q75_patient", "TYMP_q75_patient", "TYMP_high",
    }
    if not required_bridge_cols.issubset(bridge.columns):
        raise RuntimeError(f"Bridge columns missing: {sorted(required_bridge_cols.difference(bridge.columns))}")
    if bridge.index.has_duplicates or len(bridge) != 11921:
        raise RuntimeError(f"Unexpected bridge cells: {len(bridge)}")
    bridge["patient"] = bridge["patient"].astype(str)
    malignant_bridge = bridge.loc[bridge["patient"].isin(positive_patients)].copy()
    if malignant_bridge.empty or malignant_bridge["TYMP_high"].isna().any():
        raise RuntimeError("Invalid positive-q75 malignant bridge")
    high_check = malignant_bridge["TYMP_log1p_CP10K"] >= malignant_bridge["TYMP_q75_patient"]
    if not np.array_equal(high_check.to_numpy(), malignant_bridge["TYMP_high"].astype(bool).to_numpy()):
        raise RuntimeError("Frozen TYMP-high labels do not match positive-q75 thresholds")

    obj = ad.read_h5ad(args.h5ad, backed="r")
    missing_genes = pd.Index(genes).difference(obj.var_names)
    if len(missing_genes):
        raise RuntimeError(f"Reference genes missing from H5AD: {len(missing_genes)}")
    required_obs = [
        "cell_type_broad_v1", "eligible_for_subclone_selection_v1", "curated_patient_id",
        "study", "source_tissue", "total_counts_raw",
    ]
    missing_obs = sorted(set(required_obs).difference(obj.obs.columns))
    if missing_obs:
        raise RuntimeError(f"H5AD obs columns missing: {missing_obs}")
    obs = obj.obs[required_obs].copy()
    obs.index = obs.index.astype(str)
    obs["patient"] = obs["curated_patient_id"].astype(str)
    obs["study_v1"] = obs["study"].astype(str)
    obs["broad"] = obs["cell_type_broad_v1"].astype(str)
    obs["state_tymp_v1"] = pd.Series(pd.NA, index=obs.index, dtype="string")
    obs["F1"] = np.nan
    obs["TYMP_raw_count"] = np.nan
    obs["TYMP_log1p_CP10K"] = np.nan
    obs["TYMP_q75_patient"] = np.nan

    frozen_malignant = obs["eligible_for_subclone_selection_v1"].astype(bool)
    if int(frozen_malignant.sum()) != 11950:
        raise RuntimeError(f"Frozen malignant mother population changed: {int(frozen_malignant.sum())}")
    if not malignant_bridge.index.isin(obs.index).all():
        raise RuntimeError("Some bridge cells are absent from the authority H5AD")

    usage = pd.read_csv(args.usage, sep="\t", usecols=["exact_cell_id", "F1"]).set_index("exact_cell_id")
    usage.index = usage.index.astype(str)
    mapped_f1 = usage["F1"].reindex(obs.index)
    if int(mapped_f1.loc[frozen_malignant].notna().sum()) != 11950:
        raise RuntimeError("Incomplete F1 mapping for the frozen malignant population")
    obs["F1"] = mapped_f1.to_numpy()
    for column in ["TYMP_raw_count", "TYMP_log1p_CP10K", "TYMP_q75_patient"]:
        obs.loc[malignant_bridge.index, column] = malignant_bridge[column].to_numpy()
    obs.loc[malignant_bridge.index, "state_tymp_v1"] = np.where(
        malignant_bridge["TYMP_high"].astype(bool).to_numpy(),
        "TYMP-high Malignant",
        "TYMP-low Malignant",
    )

    for broad, state in BROAD_MAP.items():
        obs.loc[obs["broad"].eq(broad), "state_tymp_v1"] = state
    epithelial_normal = obs["broad"].eq("Epithelial cell") & obs["source_tissue"].astype(str).eq("Normal")
    obs.loc[epithelial_normal, "state_tymp_v1"] = "Epithelial-Normal"

    memory_gate(args.memory_gate, "before_counts", memory_rows)
    gene_positions = [int(obj.var_names.get_loc(gene)) for gene in genes]
    counts = obj.layers["counts"][:, gene_positions]
    if hasattr(counts, "to_memory"):
        counts = counts.to_memory()
    if not issparse(counts):
        counts = csr_matrix(np.asarray(counts))
    counts = counts.tocsr()
    if not np.isfinite(counts.data).all() or (counts.data < 0).any():
        raise RuntimeError("Raw counts contain negative or non-finite values")
    if not np.allclose(counts.data, np.rint(counts.data)):
        raise RuntimeError("Raw counts are not integer-like; refusing silent repair")
    gene_index = {gene: i for i, gene in enumerate(genes)}

    macro_pos = np.flatnonzero(obs["broad"].eq("Macrophage").to_numpy())
    libraries = pd.to_numeric(obs.iloc[macro_pos]["total_counts_raw"], errors="coerce").to_numpy(dtype=float)
    if not np.isfinite(libraries).all() or (libraries <= 0).any():
        raise RuntimeError("Invalid macrophage library sizes")
    spp1_score = program_score(counts, macro_pos, gene_index, SPP1_PROGRAM, libraries)
    c1qc_score = program_score(counts, macro_pos, gene_index, C1QC_PROGRAM, libraries)
    macro_frame = pd.DataFrame(
        {"patient": obs.iloc[macro_pos]["patient"].to_numpy(), "spp1": spp1_score, "c1qc": c1qc_score},
        index=obs.index[macro_pos],
    )
    macro_frame["spp1_pct"] = macro_frame.groupby("patient", observed=True)["spp1"].rank(method="average", pct=True)
    macro_frame["c1qc_pct"] = macro_frame.groupby("patient", observed=True)["c1qc"].rank(method="average", pct=True)
    macro_state = np.full(len(macro_frame), "Other Macrophage", dtype=object)
    spp1_high = (macro_frame["spp1_pct"] >= 0.75) & (macro_frame["spp1_pct"] >= macro_frame["c1qc_pct"])
    c1qc_high = (macro_frame["c1qc_pct"] >= 0.75) & ~spp1_high
    macro_state[spp1_high.to_numpy()] = "SPP1-state Macrophage"
    macro_state[c1qc_high.to_numpy()] = "C1QC-state Macrophage"
    obs.loc[macro_frame.index, "state_tymp_v1"] = macro_state
    obs.loc[macro_frame.index, "SPP1_program_score_v1"] = macro_frame["spp1"].to_numpy()
    obs.loc[macro_frame.index, "C1QC_program_score_v1"] = macro_frame["c1qc"].to_numpy()

    included = obs["state_tymp_v1"].notna().to_numpy()
    included_pos = np.flatnonzero(included)
    assignments = obs.iloc[included_pos].copy()
    observed_states = sorted(assignments["state_tymp_v1"].astype(str).unique())
    if observed_states != EXPECTED_STATES:
        raise RuntimeError(f"Expected 14 states, observed: {observed_states}")
    assignments["patient_state"] = assignments["patient"] + "||" + assignments["state_tymp_v1"].astype(str)
    patient_states = sorted(assignments["patient_state"].unique())
    codes = pd.Categorical(assignments["patient_state"], categories=patient_states).codes
    indicator = csr_matrix(
        (np.ones(len(assignments)), (codes, np.arange(len(assignments)))),
        shape=(len(patient_states), len(assignments)),
    )
    pseudobulk = indicator @ counts[included_pos]
    n_cells = np.asarray(indicator.sum(axis=1)).reshape(-1)
    patient_mean = pseudobulk.multiply((1.0 / np.maximum(n_cells, 1))[:, None]).toarray()
    meta = pd.DataFrame({"patient_state": patient_states, "n_cells": n_cells.astype(int)})
    meta[["patient", "state_tymp_v1"]] = meta["patient_state"].str.split("||", n=1, expand=True, regex=False)
    patient_study = assignments.drop_duplicates("patient").set_index("patient")["study_v1"]
    meta["study"] = meta["patient"].map(patient_study)
    if meta["study"].isna().any():
        raise RuntimeError("Patient-study mapping is incomplete")

    signatures: dict[str, np.ndarray] = {}
    qc_rows: list[dict[str, object]] = []
    for state, state_rows in meta.groupby("state_tymp_v1", observed=True):
        study_profiles = [
            patient_mean[study_rows.index].mean(axis=0)
            for _, study_rows in state_rows.groupby("study", observed=True)
        ]
        profile = np.vstack(study_profiles).mean(axis=0)
        signatures[state] = profile
        qc_rows.append(
            {
                "state_tymp_v1": state,
                "n_cells": int(state_rows["n_cells"].sum()),
                "n_patients": int(state_rows["patient"].nunique()),
                "n_studies": int(state_rows["study"].nunique()),
                "signature_sum_mean_counts": float(profile.sum()),
                "signature_nonzero_genes": int(np.sum(profile > 0)),
            }
        )
    signature = pd.DataFrame(signatures, index=genes)
    signature = signature.reindex(columns=EXPECTED_STATES)
    if signature.shape != (12758, 14) or not np.isfinite(signature.to_numpy()).all():
        raise RuntimeError(f"Invalid TYMP reference signature shape/values: {signature.shape}")
    if (signature.to_numpy() < 0).any():
        raise RuntimeError("Negative TYMP reference signatures")

    nonmalignant = sorted(set(EXPECTED_STATES).difference({"TYMP-high Malignant", "TYMP-low Malignant"}))
    old_nonmalignant = sorted(set(old_signature.columns).difference(OLD_MALIGNANT_STATES))
    if nonmalignant != old_nonmalignant:
        raise RuntimeError("Non-malignant state names differ from the frozen reference")
    delta = np.abs(signature[nonmalignant].to_numpy() - old_signature[nonmalignant].to_numpy())
    max_nonmalignant_delta = float(delta.max())
    if not np.allclose(signature[nonmalignant], old_signature[nonmalignant], rtol=1e-10, atol=1e-12):
        raise RuntimeError(f"Non-malignant signatures drifted: max abs delta={max_nonmalignant_delta}")

    high = signature["TYMP-high Malignant"]
    low = signature["TYMP-low Malignant"]
    tymp_high_value = float(high.loc["TYMP"])
    tymp_low_value = float(low.loc["TYMP"])
    malignant_rho = float(spearmanr(high.to_numpy(), low.to_numpy()).statistic)
    excluded = signature.index.str.upper().str.startswith(("MT-", "RPL", "RPS", "HBA", "HBB"))
    log2fc = np.log2((high + 0.01) / (low + 0.01))
    separable_markers = int(((np.abs(log2fc) >= 0.5) & ~excluded).sum())
    if not (tymp_high_value > tymp_low_value):
        raise RuntimeError("TYMP mean-count signature is not higher in TYMP-high malignant")
    if malignant_rho >= 0.98:
        raise RuntimeError(f"TYMP malignant signatures are not separable: rho={malignant_rho}")
    if separable_markers < 20:
        raise RuntimeError(f"Too few TYMP high/low markers: {separable_markers}")

    log_signature = np.log1p(signature)
    marker_rows: list[dict[str, object]] = []
    excluded_series = pd.Series(excluded, index=signature.index)
    for state in signature.columns:
        specificity = log_signature[state] - log_signature.drop(columns=state).mean(axis=1)
        top = specificity.loc[~excluded_series].nlargest(TOP_MARKERS)
        for rank, (gene, score) in enumerate(top.items(), 1):
            marker_rows.append(
                {"state_tymp_v1": state, "gene": gene, "rank": rank, "specificity_log1p_mean_counts": score}
            )

    malignant_assignments = assignments[assignments["state_tymp_v1"].isin(
        ["TYMP-high Malignant", "TYMP-low Malignant"]
    )]
    malignant_counts = malignant_assignments["state_tymp_v1"].value_counts()
    malignant_patient_counts = malignant_assignments.groupby("state_tymp_v1", observed=True)["patient"].nunique()
    if (malignant_counts < 500).any() or (malignant_patient_counts < 6).any():
        raise RuntimeError("TYMP malignant states fail minimum cells/patient coverage")

    assignments.index.name = "exact_cell_id"
    assignment_columns = [
        "patient", "study_v1", "state_tymp_v1", "F1", "TYMP_raw_count", "TYMP_log1p_CP10K",
        "TYMP_q75_patient", "SPP1_program_score_v1", "C1QC_program_score_v1",
    ]
    validation = pd.DataFrame(
        [
            {"check": "signature_shape_12758x14", "value": f"{signature.shape[0]}x{signature.shape[1]}", "pass": True},
            {"check": "positive_tymp_q75_patient_fraction", "value": len(positive_patients) / len(thresholds), "pass": True},
            {"check": "tymp_high_cells", "value": int(malignant_counts["TYMP-high Malignant"]), "pass": True},
            {"check": "tymp_low_cells", "value": int(malignant_counts["TYMP-low Malignant"]), "pass": True},
            {"check": "tymp_high_patients", "value": int(malignant_patient_counts["TYMP-high Malignant"]), "pass": True},
            {"check": "tymp_low_patients", "value": int(malignant_patient_counts["TYMP-low Malignant"]), "pass": True},
            {"check": "tymp_signature_high_gt_low", "value": f"{tymp_high_value:.9g}>{tymp_low_value:.9g}", "pass": True},
            {"check": "tymp_high_low_signature_rho_lt_0.98", "value": malignant_rho, "pass": malignant_rho < 0.98},
            {"check": "tymp_high_low_markers_abs_log2fc_ge_0.5", "value": separable_markers, "pass": separable_markers >= 20},
            {"check": "nonmalignant_signatures_unchanged", "value": max_nonmalignant_delta, "pass": max_nonmalignant_delta <= 1e-12},
            {"check": "raw_counts_integer_nonnegative_finite", "value": int(counts.nnz), "pass": True},
        ]
    )

    atomic_tsv(signature, work / "reference_tymp_v1_signatures_study_balanced_mean_counts.tsv.gz", index=True)
    atomic_tsv(pd.DataFrame(qc_rows).sort_values("state_tymp_v1"), work / "reference_tymp_v1_state_qc.tsv")
    atomic_tsv(meta, work / "reference_tymp_v1_patient_state_qc.tsv")
    atomic_tsv(assignments[assignment_columns], work / "reference_tymp_v1_cell_assignments.tsv.gz", index=True)
    atomic_tsv(pd.DataFrame(marker_rows), work / "reference_tymp_v1_marker_ranking.tsv")
    atomic_tsv(validation, work / "reference_tymp_v1_internal_validation.tsv")
    atomic_tsv(pd.DataFrame(memory_rows), work / "memory_gates_tymp_v1.tsv")
    excluded_patients = thresholds.loc[~thresholds["patient"].isin(positive_patients)].copy()
    excluded_patients["exclusion_reason"] = "TYMP q75 == 0; cannot form stable within-patient high state"
    atomic_tsv(excluded_patients, work / "reference_tymp_v1_excluded_malignant_patients.tsv")

    metadata: dict[str, object] = {
        "schema_version": "phase17_reference_tymp_v1",
        "status": "PASS_INTERNAL_VALIDATION",
        "created_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "inputs": {
            "h5ad": str(args.h5ad),
            "h5ad_verified_sha256": "e4779c46af13a8643e8dec584eeff4322258181f5edab3e55dcd8c20bc28a9e7",
            "usage": str(args.usage),
            "usage_sha256": sha256(args.usage),
            "bridge": str(args.bridge),
            "bridge_sha256": sha256(args.bridge),
            "thresholds": str(args.thresholds),
            "thresholds_sha256": sha256(args.thresholds),
            "old_signature": str(args.old_signature),
            "old_signature_sha256": sha256(args.old_signature),
        },
        "n_input_cells": int(obj.n_obs),
        "n_reference_cells": int(len(assignments)),
        "n_genes": int(signature.shape[0]),
        "n_states": int(signature.shape[1]),
        "positive_tymp_q75_patients": len(positive_patients),
        "excluded_zero_q75_patients": sorted(set(thresholds["patient"]).difference(positive_patients)),
        "malignant_definition": "Part1 rebuild eligible_for_subclone_selection_v1; patient >=20 cells; positive TYMP q75; patient-specific TYMP upper quartile",
        "f1_role": "continuous QC/covariate only; not a TYMP state label",
        "macrophage_definition": "unchanged from reference-v3: within-patient program percentiles with mutually exclusive SPP1/C1QC upper-quartile states",
        "aggregation": "patient mean raw counts -> study mean -> equal-weight study mean",
        "validation": {
            "tymp_high_mean_count": tymp_high_value,
            "tymp_low_mean_count": tymp_low_value,
            "tymp_high_low_spearman_rho": malignant_rho,
            "n_abs_log2fc_ge_0.5_markers": separable_markers,
            "max_nonmalignant_signature_abs_delta_vs_v3": max_nonmalignant_delta,
        },
        "versions": {
            "python": platform.python_version(), "anndata": ad.__version__, "numpy": np.__version__,
            "pandas": pd.__version__, "scipy": scipy.__version__,
        },
    }
    atomic_json(metadata, work / "run_metadata_tymp_v1.json")
    if getattr(obj, "file", None) is not None:
        obj.file.close()

    manifest_rows = []
    for path in sorted(work.iterdir()):
        if path.is_file() and path.name != "file_manifest_tymp_v1.tsv":
            manifest_rows.append({"file": path.name, "size_bytes": path.stat().st_size, "sha256": sha256(path)})
    atomic_tsv(pd.DataFrame(manifest_rows), work / "file_manifest_tymp_v1.tsv")
    os.replace(work, args.output)
    print(json.dumps(metadata, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
