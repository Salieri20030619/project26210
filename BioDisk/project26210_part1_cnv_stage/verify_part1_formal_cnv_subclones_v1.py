#!/usr/bin/env python3
"""Independent terminal verification for formal inferred CNV clone calls v1."""

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
EXTENSION = ROOT / "04_subclones_tymp/cnv_subclone_threshold_extension_v1"
SOURCE = ROOT / "04_subclones_tymp/cnv_subclones_v1"
META = ROOT / "05_cnmf/input_v1/cnmf_cell_metadata_v1.tsv.gz"
OUT = ROOT / "04_subclones_tymp/cnv_subclones_verification_v1"
MEMORY_GATE = PROJECT / "check_memory_gate.sh"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    if OUT.exists() and any(OUT.iterdir()):
        raise RuntimeError(f"Refusing to overwrite non-empty formal verification directory: {OUT}")
    OUT.mkdir(parents=True, exist_ok=True)
    checks: list[dict[str, object]] = []

    def check(name: str, passed: bool, detail: str) -> None:
        checks.append({"check": name, "status": "PASS" if passed else "FAIL", "detail": detail})
        if not passed:
            pd.DataFrame(checks).to_csv(OUT / "validation_results_v1.tsv", sep="\t", index=False)
            raise RuntimeError(f"Verification failed: {name}: {detail}")

    gate = subprocess.run([str(MEMORY_GATE), "formal_cnv_subclones_verification_v1"], text=True, capture_output=True)
    gate_text = " ".join((gate.stdout + " " + gate.stderr).split())
    (OUT / "memory_gate_v1.txt").write_text(gate_text + "\n", encoding="utf-8")
    check("memory_gate", gate.returncode == 0 and "PASS" in gate_text, gate_text)
    check("source_complete", (SOURCE / "RUN_COMPLETE").is_file(), str(SOURCE))

    manifest = pd.read_csv(SOURCE / "file_manifest_v1.tsv", sep="\t")
    for row in manifest.itertuples(index=False):
        path = SOURCE / row.file
        check(f"file_exists::{row.file}", path.is_file(), str(path))
        check(f"file_size::{row.file}", path.stat().st_size == row.size_bytes, f"observed={path.stat().st_size}; expected={row.size_bytes}")
        check(f"file_sha256::{row.file}", sha256(path) == row.sha256, row.sha256)

    metadata = json.loads((SOURCE / "run_metadata_v1.json").read_text(encoding="utf-8"))
    check("selected_threshold_metadata", metadata["selected_distance_threshold"] == 0.25, str(metadata["selected_distance_threshold"]))
    check("interpretation_guard", metadata["interpretation"] == "inferred CNV clone; not mutation-validated genetic clone", metadata["interpretation"])
    check("cnmf_distinction_guard", "transcriptional programs" in metadata["cnmf_distinction"], metadata["cnmf_distinction"])
    check("mandatory_sensitivity", np.allclose(metadata["mandatory_sensitivity_thresholds"], [0.225, 0.275]), str(metadata["mandatory_sensitivity_thresholds"]))

    extension = pd.read_csv(EXTENSION / "group_to_clone_extended_grid_v1.tsv.gz", sep="\t")
    expected_groups = extension.loc[np.isclose(extension["distance_threshold"], 0.25)].sort_values(["infercnv_unit", "infercnv_group"], kind="stable")
    groups = pd.read_csv(SOURCE / "formal_group_to_cnv_clone_v1.tsv.gz", sep="\t").sort_values(["infercnv_unit", "infercnv_group"], kind="stable")
    clones = pd.read_csv(SOURCE / "formal_cnv_clone_summary_v1.tsv", sep="\t")
    patients = pd.read_csv(SOURCE / "formal_patient_clone_summary_v1.tsv", sep="\t")
    cells = pd.read_csv(SOURCE / "formal_cell_cnv_clone_assignments_v1.tsv.gz", sep="\t", dtype={"exact_cell_id": str})
    sensitivity = pd.read_csv(SOURCE / "patient_threshold_sensitivity_v1.tsv", sep="\t")
    sensitivity_summary = pd.read_csv(SOURCE / "patient_threshold_sensitivity_summary_v1.tsv", sep="\t")
    evidence = pd.read_csv(SOURCE / "selected_threshold_evidence_v1.tsv", sep="\t")
    meta = pd.read_csv(META, sep="\t", dtype={"exact_cell_id": str})

    check("formal_group_rows", len(groups) == 264, str(len(groups)))
    check("formal_group_unique", not groups[["infercnv_unit", "infercnv_group"]].duplicated().any(), "unit/group")
    for column in ["clone_number", "clone_id", "clone_cells", "raw_group_cells"]:
        observed = groups[column].reset_index(drop=True)
        expected = expected_groups[column].reset_index(drop=True)
        check(f"formal_group_exact_extension::{column}", observed.equals(expected), f"mismatches={int((observed != expected).sum())}")
    expected_class = np.where(groups["clone_cells"].ge(20), "primary_clone_ge20_v1", "minor_clone_lt20_v1")
    check("clone_size_class", np.array_equal(groups["clone_size_class_v1"], expected_class), str(groups["clone_size_class_v1"].value_counts().to_dict()))

    check("clone_count", len(clones) == 159, str(len(clones)))
    check("primary_clone_count", clones["clone_size_class_v1"].eq("primary_clone_ge20_v1").sum() == 48, str(clones["clone_size_class_v1"].value_counts().to_dict()))
    check("clone_cell_conservation", int(clones["clone_cells"].sum()) == 11950, str(clones["clone_cells"].sum()))
    group_clone_counts = groups.groupby("formal_inferred_cnv_clone_v1")["raw_group_cells"].sum().sort_index()
    declared_clone_counts = clones.set_index("formal_inferred_cnv_clone_v1")["clone_cells"].sort_index()
    check("clone_counts_recomputed", group_clone_counts.equals(declared_clone_counts), f"n={len(group_clone_counts)}")

    check("patient_count", len(patients) == 26 and not patients["infercnv_unit"].duplicated().any(), str(len(patients)))
    check("patient_cell_conservation", int(patients["n_selected_cells"].sum()) == 11950, str(patients["n_selected_cells"].sum()))
    check("evaluable_patient_count", int(patients["evaluable_within_patient_two_primary_clones_v1"].sum()) == 14, str(patients["evaluable_within_patient_two_primary_clones_v1"].sum()))
    check("sensitive_patient_count", int(patients["threshold_sensitive_patient_v1"].sum()) == 5, str(patients["threshold_sensitive_patient_v1"].sum()))
    check("primary_cell_count", int(patients["n_cells_primary_clones_ge20"].sum()) == 11604, str(patients["n_cells_primary_clones_ge20"].sum()))
    check("minor_cell_count", int(patients["n_cells_minor_clones_lt20"].sum()) == 346, str(patients["n_cells_minor_clones_lt20"].sum()))

    check("cell_rows", len(cells) == 11950 and not cells["exact_cell_id"].duplicated().any(), str(len(cells)))
    check("cell_identity", set(cells["exact_cell_id"]) == set(meta["exact_cell_id"]), f"formal={len(cells)}; meta={len(meta)}")
    expected_cell_map = cells.groupby("formal_inferred_cnv_clone_v1").size().sort_index()
    check("cell_clone_counts", expected_cell_map.equals(declared_clone_counts), f"n={len(expected_cell_map)}")
    expected_eligible = cells["clone_size_class_v1"].eq("primary_clone_ge20_v1") & cells["evaluable_within_patient_two_primary_clones_v1"]
    check("cell_primary_eligibility", expected_eligible.equals(cells["eligible_for_primary_within_patient_clone_analysis_v1"]), str(expected_eligible.sum()))

    recomputed_sensitivity = []
    for unit, frame in extension.groupby("infercnv_unit", sort=True, observed=True):
        for threshold_a, threshold_b in ((0.225, 0.25), (0.25, 0.275)):
            left = frame.loc[np.isclose(frame["distance_threshold"], threshold_a)].sort_values("infercnv_group", kind="stable")
            right = frame.loc[np.isclose(frame["distance_threshold"], threshold_b)].sort_values("infercnv_group", kind="stable")
            weights = left["raw_group_cells"].to_numpy(int)
            ari = adjusted_rand_score(np.repeat(left["clone_number"].to_numpy(int), weights), np.repeat(right["clone_number"].to_numpy(int), weights))
            recomputed_sensitivity.append((unit, threshold_a, threshold_b, ari))
    recomputed = pd.DataFrame(recomputed_sensitivity, columns=["infercnv_unit", "threshold_a", "threshold_b", "ari"])
    declared = sensitivity[["infercnv_unit", "threshold_a", "threshold_b", "cell_weighted_patient_adjusted_rand_index"]].rename(columns={"cell_weighted_patient_adjusted_rand_index": "ari"})
    compared = recomputed.merge(declared, on=["infercnv_unit", "threshold_a", "threshold_b"], suffixes=("_recomputed", "_declared"), validate="one_to_one")
    check("patient_sensitivity_rows", len(compared) == 52, str(len(compared)))
    check("patient_sensitivity_recomputed", np.allclose(compared["ari_recomputed"], compared["ari_declared"], atol=1e-12, rtol=0), str(np.max(np.abs(compared["ari_recomputed"] - compared["ari_declared"]))))
    check("sensitivity_summary_medians", np.allclose(sensitivity_summary["median_patient_ari"], [1.0, 1.0]), str(sensitivity_summary["median_patient_ari"].tolist()))

    check("selection_evidence_rows", len(evidence) == 1, str(len(evidence)))
    check("selection_plateau", evidence.iloc[0]["maximum_plateau_thresholds"] == "0.225,0.250,0.275", str(evidence.iloc[0]["maximum_plateau_thresholds"]))
    check("selection_adjacent_ari", evidence.iloc[0]["lower_adjacent_cell_weighted_ari"] >= 0.90 and evidence.iloc[0]["upper_adjacent_cell_weighted_ari"] >= 0.90, str(evidence.iloc[0].to_dict()))
    check("selection_coverage", evidence.iloc[0]["fraction_cells_in_primary_clones_ge20"] >= 0.95, str(evidence.iloc[0]["fraction_cells_in_primary_clones_ge20"]))

    results = pd.DataFrame(checks)
    results.to_csv(OUT / "validation_results_v1.tsv", sep="\t", index=False)
    summary = {
        "schema_version": "part1_formal_cnv_subclones_verification_v1",
        "completed_at": datetime.now().astimezone().isoformat(),
        "n_checks": int(len(results)),
        "n_pass": int(results["status"].eq("PASS").sum()),
        "n_fail": int(results["status"].eq("FAIL").sum()),
        "n_cells": 11950,
        "n_inferred_cnv_clones": 159,
        "n_primary_clones_ge20": 48,
        "n_evaluable_patients": 14,
        "dependency_changes": False,
    }
    (OUT / "verification_summary_v1.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    files = sorted(path for path in OUT.iterdir() if path.name not in {"file_manifest_v1.tsv", "RUN_COMPLETE"})
    pd.DataFrame({
        "file": [path.name for path in files],
        "size_bytes": [path.stat().st_size for path in files],
        "sha256": [sha256(path) for path in files],
    }).to_csv(OUT / "file_manifest_v1.tsv", sep="\t", index=False)
    (OUT / "RUN_COMPLETE").write_text(datetime.now().astimezone().isoformat() + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))
    print("PART1_FORMAL_CNV_SUBCLONES_VERIFICATION_V1_PASS")


if __name__ == "__main__":
    main()
