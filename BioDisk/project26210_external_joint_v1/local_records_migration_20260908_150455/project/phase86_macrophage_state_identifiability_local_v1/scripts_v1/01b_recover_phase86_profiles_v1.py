#!/usr/bin/env python

from __future__ import annotations

import argparse
import hashlib
import os
import shutil
from pathlib import Path

import numpy as np
import pandas as pd


EXPECTED_DATASETS = {"GSE207493", "GSE264586", "GSE294109"}
EXPECTED_STATES = {
    "B cell", "C1QC-state Macrophage", "Dendritic cell", "Endothelial",
    "Epithelial-Normal", "Malignant", "Mast cell", "Monocyte", "NK cell",
    "Other Macrophage", "SMC/Pericyte", "SPP1-state Macrophage", "T cell",
}
CORE_FILES = [
    "patient_state_profiles_v1.tsv.gz",
    "patient_state_profile_counts_v1.tsv",
    "profile_extraction_status_v1.tsv",
    "profile_gene_coverage_v1.tsv",
]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--phase-root", required=True, type=Path)
    args = parser.parse_args()
    root = args.phase_root.resolve()
    failed = root / "profiles_v1_failed_attempt1"
    output = root / "profiles_v1"
    part = root / "profiles_v1.part"
    if not failed.is_dir() or output.exists() or part.exists():
        raise RuntimeError("Invalid recovery directory state")
    for name in CORE_FILES + ["RUN_FAILURE_v1.txt"]:
        if not (failed / name).is_file():
            raise RuntimeError(f"Missing failed-attempt artifact: {name}")

    profiles = pd.read_csv(failed / CORE_FILES[0], sep="\t")
    counts = pd.read_csv(failed / CORE_FILES[1], sep="\t")
    status = pd.read_csv(failed / CORE_FILES[2], sep="\t")
    coverage = pd.read_csv(failed / CORE_FILES[3], sep="\t")
    metadata_columns = ["dataset", "patient_id", "state", "n_cells"]
    genes = [column for column in profiles.columns if column not in metadata_columns]
    matrix = profiles[genes].to_numpy(dtype=float)

    checks = []

    def record(name: str, value, passed: bool) -> None:
        checks.append({"check": name, "value": value, "pass": bool(passed)})

    record("datasets_exact", ";".join(sorted(profiles["dataset"].unique())), set(profiles["dataset"]) == EXPECTED_DATASETS)
    record("states_expected", profiles["state"].nunique(), set(profiles["state"]).issubset(EXPECTED_STATES) and set(profiles["state"]) == EXPECTED_STATES)
    record("profile_key_unique", len(profiles), not profiles.duplicated(["dataset", "patient_id", "state"]).any())
    record("common_genes_478", len(genes), len(genes) == 478)
    record("counts_finite_nonnegative_integer", matrix.size, bool(np.isfinite(matrix).all() and (matrix >= 0).all() and np.allclose(matrix, np.round(matrix))))
    record("cell_counts_positive_integer", profiles["n_cells"].sum(), bool((profiles["n_cells"] > 0).all() and np.allclose(profiles["n_cells"], np.round(profiles["n_cells"]))))
    record("coverage_three_rows", len(coverage), len(coverage) == 3 and set(coverage["dataset"]) == EXPECTED_DATASETS)
    record("coverage_at_least_95pct", coverage["present_genes"].min(), bool((coverage["present_genes"] >= 0.95 * coverage["panel_genes"]).all()))
    record("status_three_pass_rows", len(status), len(status) == 3 and set(status["dataset"]) == EXPECTED_DATASETS and status["status"].eq("PASS_PROFILE_EXTRACTION").all())
    record("source_bytes_exact", status["source_bytes"].sum(), set(status["source_bytes"].astype(int)) == {1923312451, 941486138, 1462703849})

    recomputed = profiles.groupby(["dataset", "state"], sort=True).agg(patients=("patient_id", "nunique"), cells=("n_cells", "sum")).reset_index()
    stored = counts.sort_values(["dataset", "state"]).reset_index(drop=True)
    recomputed = recomputed.sort_values(["dataset", "state"]).reset_index(drop=True)
    record("state_count_table_recomputed", len(stored), stored.equals(recomputed))
    record("failed_marker_present", "RUN_FAILURE_v1.txt", (failed / "RUN_FAILURE_v1.txt").is_file())

    check_frame = pd.DataFrame(checks)
    if not check_frame["pass"].all():
        raise RuntimeError(f"Recovery validation failed: {check_frame.loc[~check_frame['pass']].to_dict('records')}")

    part.mkdir(parents=True)
    for name in CORE_FILES:
        shutil.copy2(failed / name, part / name)
    check_frame.to_csv(part / "PROFILE_RECOVERY_VALIDATION_v1.tsv", sep="\t", index=False)
    metadata = pd.DataFrame([
        {"record_key": "status", "value": "PASS_PROFILE_EXTRACTION_RECOVERED"},
        {"record_key": "recovery_source", "value": str(failed)},
        {"record_key": "failure_boundary", "value": "AFTER_CORE_PROFILE_WRITE_BEFORE_FORMAL_DIRECTORY"},
        {"record_key": "common_platform_genes", "value": len(genes)},
        {"record_key": "datasets", "value": profiles["dataset"].nunique()},
        {"record_key": "profile_rows", "value": len(profiles)},
        {"record_key": "total_cells", "value": int(profiles["n_cells"].sum())},
        {"record_key": "recovery_checks", "value": len(check_frame)},
    ])
    metadata.to_csv(part / "profile_run_metadata_v1.tsv", sep="\t", index=False)
    (part / "RECOVERY_RECORD_v1.md").write_text(
        "# Phase 86 profile extraction recovery\n\n"
        "The initial R run completed all patient-state count aggregation and safety checks, then failed while constructing the final run-metadata table because `data.table(key=...)` treated `key` as a set-key argument. The failed directory is preserved unchanged. This recovery independently reloaded and validated every written core table before copying it into the formal non-overwriting `profiles_v1` directory. No Seurat object was reopened and no count was recomputed or modified during recovery.\n",
        encoding="utf-8",
    )
    (part / "PHASE86_PROFILE_EXTRACTION_PASS_v1").write_text("PASS_PROFILE_EXTRACTION_RECOVERED\n", encoding="utf-8")
    manifest = []
    for path in sorted(part.iterdir()):
        if path.name == "PROFILE_FILE_MANIFEST_SHA256_v1.tsv":
            continue
        manifest.append({"file": path.name, "bytes": path.stat().st_size, "sha256": sha256(path)})
    pd.DataFrame(manifest).to_csv(part / "PROFILE_FILE_MANIFEST_SHA256_v1.tsv", sep="\t", index=False)
    os.replace(part, output)
    print(f"PASS_PROFILE_EXTRACTION_RECOVERED profiles={len(profiles)} genes={len(genes)} checks={len(check_frame)}")


if __name__ == "__main__":
    main()
