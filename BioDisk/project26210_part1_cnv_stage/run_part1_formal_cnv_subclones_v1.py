#!/usr/bin/env python3
"""Select the audited CNV consolidation threshold and create formal clone calls."""

from __future__ import annotations

import hashlib
import json
import subprocess
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score


PROJECT = Path("/mnt/BioDisk/ludingyu/project26210")
ROOT = PROJECT / "results/part1_rebuild_v1"
PREFLIGHT = ROOT / "04_subclones_tymp/cnv_subclone_preflight_v1"
EXTENSION = ROOT / "04_subclones_tymp/cnv_subclone_threshold_extension_v1"
VERIFICATION = ROOT / "04_subclones_tymp/cnv_subclone_threshold_extension_verification_v1"
META = ROOT / "05_cnmf/input_v1/cnmf_cell_metadata_v1.tsv.gz"
OUT = ROOT / "04_subclones_tymp/cnv_subclones_v1"
MEMORY_GATE = PROJECT / "check_memory_gate.sh"
PRIMARY_MINIMUM = 20
MIN_ADJACENT_ARI = 0.90
MIN_PRIMARY_COVERAGE = 0.95


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    if OUT.exists() and any(OUT.iterdir()):
        raise RuntimeError(f"Refusing to overwrite non-empty formal clone directory: {OUT}")
    OUT.mkdir(parents=True, exist_ok=True)
    gate = subprocess.run([str(MEMORY_GATE), "formal_cnv_subclones_v1"], text=True, capture_output=True)
    gate_text = " ".join((gate.stdout + " " + gate.stderr).split())
    (OUT / "memory_gate_v1.txt").write_text(gate_text + "\n", encoding="utf-8")
    if gate.returncode != 0 or "PASS" not in gate_text:
        raise RuntimeError(f"Memory gate blocked formal clone call: {gate_text}")
    if not (VERIFICATION / "RUN_COMPLETE").is_file():
        raise RuntimeError("Independently verified threshold extension is required")

    cohort = pd.read_csv(EXTENSION / "extended_threshold_cohort_summary_v1.tsv", sep="\t")
    stability = pd.read_csv(EXTENSION / "extended_threshold_stability_v1.tsv", sep="\t")
    grid = pd.read_csv(EXTENSION / "group_to_clone_extended_grid_v1.tsv.gz", sep="\t")
    inventory = pd.read_csv(PREFLIGHT / "raw_hmm_group_profile_inventory_v1.tsv", sep="\t")
    meta = pd.read_csv(META, sep="\t", dtype={"exact_cell_id": str, "infercnv_unit": str, "infercnv_group": str})

    maximum_evaluable = int(cohort["n_evaluable_patients_two_primary_clones_min20"].max())
    plateau = cohort.loc[cohort["n_evaluable_patients_two_primary_clones_min20"].eq(maximum_evaluable), "distance_threshold"].sort_values().to_numpy(float)
    if len(plateau) % 2 != 1 or len(plateau) < 3 or not np.allclose(np.diff(plateau), 0.025):
        raise RuntimeError(f"Selection plateau is not an odd contiguous grid: {plateau}")
    selected_threshold = float(plateau[len(plateau) // 2])
    lower_threshold = float(plateau[len(plateau) // 2 - 1])
    upper_threshold = float(plateau[len(plateau) // 2 + 1])
    if not np.isclose(selected_threshold, 0.25) or not np.allclose([lower_threshold, upper_threshold], [0.225, 0.275]):
        raise RuntimeError(f"Unexpected threshold selection: {lower_threshold}, {selected_threshold}, {upper_threshold}")

    selected_cohort = cohort.loc[np.isclose(cohort["distance_threshold"], selected_threshold)].iloc[0]
    adjacent = stability.loc[
        (np.isclose(stability["threshold_a"], lower_threshold) & np.isclose(stability["threshold_b"], selected_threshold))
        | (np.isclose(stability["threshold_a"], selected_threshold) & np.isclose(stability["threshold_b"], upper_threshold))
    ].sort_values("threshold_a")
    if len(adjacent) != 2 or adjacent["cell_weighted_adjusted_rand_index"].min() < MIN_ADJACENT_ARI:
        raise RuntimeError("Selected threshold lacks required adjacent stability")
    if float(selected_cohort["fraction_cells_in_primary_clones_min20"]) < MIN_PRIMARY_COVERAGE:
        raise RuntimeError("Selected threshold lacks required primary-clone cell coverage")

    sensitivity_rows = []
    for unit, frame in grid.groupby("infercnv_unit", sort=True, observed=True):
        for threshold_a, threshold_b in ((lower_threshold, selected_threshold), (selected_threshold, upper_threshold)):
            left = frame.loc[np.isclose(frame["distance_threshold"], threshold_a)].sort_values("infercnv_group", kind="stable")
            right = frame.loc[np.isclose(frame["distance_threshold"], threshold_b)].sort_values("infercnv_group", kind="stable")
            if not left["infercnv_group"].reset_index(drop=True).equals(right["infercnv_group"].reset_index(drop=True)):
                raise RuntimeError(f"{unit}: group order mismatch in threshold sensitivity")
            weights = left["raw_group_cells"].to_numpy(int)
            ari = adjusted_rand_score(
                np.repeat(left["clone_number"].to_numpy(int), weights),
                np.repeat(right["clone_number"].to_numpy(int), weights),
            )
            sensitivity_rows.append({
                "infercnv_unit": unit,
                "threshold_a": threshold_a,
                "threshold_b": threshold_b,
                "n_raw_hmm_groups": int(len(left)),
                "n_selected_cells": int(weights.sum()),
                "cell_weighted_patient_adjusted_rand_index": float(ari),
                "patient_threshold_sensitive_ari_lt_0_8": bool(ari < 0.8),
            })
    sensitivity = pd.DataFrame(sensitivity_rows)
    sensitivity_summary = (
        sensitivity.groupby(["threshold_a", "threshold_b"], observed=True)
        .agg(
            n_patients=("infercnv_unit", "size"),
            mean_patient_ari=("cell_weighted_patient_adjusted_rand_index", "mean"),
            median_patient_ari=("cell_weighted_patient_adjusted_rand_index", "median"),
            minimum_patient_ari=("cell_weighted_patient_adjusted_rand_index", "min"),
            n_sensitive_patients_ari_lt_0_8=("patient_threshold_sensitive_ari_lt_0_8", "sum"),
        ).reset_index()
    )

    selected_groups = grid.loc[np.isclose(grid["distance_threshold"], selected_threshold)].copy()
    selected_groups = selected_groups.merge(
        inventory,
        on=["infercnv_unit", "infercnv_group"],
        how="left",
        validate="one_to_one",
        suffixes=("", "_inventory"),
    )
    if len(selected_groups) != 264 or selected_groups["n_selected_cells"].isna().any():
        raise RuntimeError("Formal group assignment merge failed")
    if not selected_groups["raw_group_cells"].equals(selected_groups["n_selected_cells"]):
        raise RuntimeError("Formal group assignment cell counts disagree")

    selected_groups["clone_size_class_v1"] = np.where(
        selected_groups["clone_cells"].ge(PRIMARY_MINIMUM),
        "primary_clone_ge20_v1",
        "minor_clone_lt20_v1",
    )
    selected_groups["formal_inferred_cnv_clone_v1"] = selected_groups["clone_id"]
    clone_summary = (
        selected_groups.groupby(["infercnv_unit", "formal_inferred_cnv_clone_v1", "clone_number", "clone_size_class_v1"], observed=True)
        .agg(
            clone_cells=("raw_group_cells", "sum"),
            n_raw_hmm_groups=("infercnv_group", "size"),
            weighted_mean_abs_bin_state_shift=("mean_abs_bin_state_shift", lambda x: float(np.average(x, weights=selected_groups.loc[x.index, "raw_group_cells"]))),
            minimum_altered_5mb_bins=("n_altered_5mb_bins", "min"),
            maximum_altered_5mb_bins=("n_altered_5mb_bins", "max"),
        ).reset_index()
    )
    if int(clone_summary["clone_cells"].sum()) != 11950:
        raise RuntimeError("Clone summary failed cell conservation")

    patient_summary = (
        clone_summary.groupby("infercnv_unit", observed=True)
        .agg(
            n_selected_cells=("clone_cells", "sum"),
            n_inferred_cnv_clones=("formal_inferred_cnv_clone_v1", "size"),
            n_primary_clones_ge20=("clone_size_class_v1", lambda x: int(x.eq("primary_clone_ge20_v1").sum())),
            n_cells_primary_clones_ge20=("clone_cells", lambda x: int(x[clone_summary.loc[x.index, "clone_size_class_v1"].eq("primary_clone_ge20_v1")].sum())),
        ).reset_index()
    )
    patient_summary["n_cells_minor_clones_lt20"] = patient_summary["n_selected_cells"] - patient_summary["n_cells_primary_clones_ge20"]
    patient_summary["fraction_cells_primary_clones_ge20"] = patient_summary["n_cells_primary_clones_ge20"] / patient_summary["n_selected_cells"]
    patient_summary["evaluable_within_patient_two_primary_clones_v1"] = patient_summary["n_primary_clones_ge20"].ge(2)

    sensitive_units = set(sensitivity.loc[sensitivity["patient_threshold_sensitive_ari_lt_0_8"], "infercnv_unit"])
    patient_summary["threshold_sensitive_patient_v1"] = patient_summary["infercnv_unit"].isin(sensitive_units)
    patient_summary["threshold_sensitivity_required_v1"] = True

    mapping = selected_groups[["infercnv_unit", "infercnv_group", "formal_inferred_cnv_clone_v1", "clone_number", "clone_cells", "clone_size_class_v1"]]
    cell_columns = [
        "exact_cell_id", "study", "curated_patient_id", "source_patient_id", "infercnv_unit", "infercnv_group",
        "copykat_call_v1", "malignancy_status_v1", "TYMP_raw_count_v1", "TYMP_detected_v1", "TYMP_log1p_CP10K_v1",
    ]
    cell_assignments = meta[cell_columns].merge(mapping, on=["infercnv_unit", "infercnv_group"], how="left", validate="many_to_one")
    cell_assignments = cell_assignments.merge(
        patient_summary[["infercnv_unit", "evaluable_within_patient_two_primary_clones_v1", "threshold_sensitive_patient_v1"]],
        on="infercnv_unit",
        how="left",
        validate="many_to_one",
    )
    cell_assignments["eligible_for_primary_within_patient_clone_analysis_v1"] = (
        cell_assignments["clone_size_class_v1"].eq("primary_clone_ge20_v1")
        & cell_assignments["evaluable_within_patient_two_primary_clones_v1"]
    )
    if len(cell_assignments) != 11950 or cell_assignments["formal_inferred_cnv_clone_v1"].isna().any():
        raise RuntimeError("Formal cell assignments failed conservation")

    selection_evidence = pd.DataFrame([{
        "selected_distance_threshold": selected_threshold,
        "selection_rule": "maximize patients with >=2 clones of >=20 cells; select center of odd contiguous maximum plateau; require both adjacent cell-weighted ARI >=0.90 and primary-clone coverage >=0.95",
        "maximum_evaluable_patients": maximum_evaluable,
        "maximum_plateau_thresholds": ",".join(f"{value:.3f}" for value in plateau),
        "lower_sensitivity_threshold": lower_threshold,
        "upper_sensitivity_threshold": upper_threshold,
        "lower_adjacent_cell_weighted_ari": float(adjacent.iloc[0]["cell_weighted_adjusted_rand_index"]),
        "upper_adjacent_cell_weighted_ari": float(adjacent.iloc[1]["cell_weighted_adjusted_rand_index"]),
        "fraction_cells_in_primary_clones_ge20": float(selected_cohort["fraction_cells_in_primary_clones_min20"]),
        "n_total_inferred_cnv_clones": int(selected_cohort["n_consolidated_clones"]),
        "n_primary_clones_ge20": int(selected_cohort["n_primary_clones_min20"]),
        "n_evaluable_patients_two_primary_clones_ge20": int(selected_cohort["n_evaluable_patients_two_primary_clones_min20"]),
        "selection_status": "selected_with_patient_level_threshold_sensitivity_warning",
    }])

    selected_groups.to_csv(OUT / "formal_group_to_cnv_clone_v1.tsv.gz", sep="\t", index=False, compression="gzip")
    clone_summary.to_csv(OUT / "formal_cnv_clone_summary_v1.tsv", sep="\t", index=False)
    patient_summary.to_csv(OUT / "formal_patient_clone_summary_v1.tsv", sep="\t", index=False)
    cell_assignments.to_csv(OUT / "formal_cell_cnv_clone_assignments_v1.tsv.gz", sep="\t", index=False, compression="gzip")
    sensitivity.to_csv(OUT / "patient_threshold_sensitivity_v1.tsv", sep="\t", index=False)
    sensitivity_summary.to_csv(OUT / "patient_threshold_sensitivity_summary_v1.tsv", sep="\t", index=False)
    selection_evidence.to_csv(OUT / "selected_threshold_evidence_v1.tsv", sep="\t", index=False)

    metadata = {
        "schema_version": "part1_formal_inferred_cnv_subclones_v1",
        "completed_at": datetime.now().astimezone().isoformat(),
        "selected_distance_threshold": selected_threshold,
        "primary_clone_minimum_cells": PRIMARY_MINIMUM,
        "n_cells": 11950,
        "n_patients": 26,
        "n_total_inferred_cnv_clones": int(len(clone_summary)),
        "n_primary_clones_ge20": int(clone_summary["clone_size_class_v1"].eq("primary_clone_ge20_v1").sum()),
        "n_evaluable_patients": int(patient_summary["evaluable_within_patient_two_primary_clones_v1"].sum()),
        "n_threshold_sensitive_patients": int(patient_summary["threshold_sensitive_patient_v1"].sum()),
        "mandatory_sensitivity_thresholds": [lower_threshold, upper_threshold],
        "interpretation": "inferred CNV clone; not mutation-validated genetic clone",
        "cnmf_distinction": "cNMF factors are transcriptional programs and must not be relabeled as CNV clones",
        "dependency_changes": False,
    }
    (OUT / "run_metadata_v1.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    files = sorted(path for path in OUT.iterdir() if path.name not in {"file_manifest_v1.tsv", "RUN_COMPLETE"})
    pd.DataFrame({
        "file": [path.name for path in files],
        "size_bytes": [path.stat().st_size for path in files],
        "sha256": [sha256(path) for path in files],
    }).to_csv(OUT / "file_manifest_v1.tsv", sep="\t", index=False)
    (OUT / "RUN_COMPLETE").write_text(datetime.now().astimezone().isoformat() + "\n", encoding="utf-8")
    print(selection_evidence.to_string(index=False))
    print(sensitivity_summary.to_string(index=False))
    print(patient_summary[["n_selected_cells", "n_inferred_cnv_clones", "n_primary_clones_ge20", "evaluable_within_patient_two_primary_clones_v1", "threshold_sensitive_patient_v1"]].sum(numeric_only=True).to_string())
    print("PART1_FORMAL_CNV_SUBCLONES_V1_PASS")


if __name__ == "__main__":
    main()
