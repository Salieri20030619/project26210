#!/usr/bin/env python3
"""Independent verification of patient-level CNV clone/TYMP/cNMF association v1."""

from __future__ import annotations

import hashlib
import json
import subprocess
from datetime import datetime
from pathlib import Path

import matplotlib.image as mpimg
import numpy as np
import pandas as pd
from scipy.stats import binomtest, wilcoxon


PROJECT = Path("/mnt/BioDisk/ludingyu/project26210")
ROOT = PROJECT / "results/part1_rebuild_v1"
SOURCE = ROOT / "04_subclones_tymp/cnv_clone_tymp_cnmf_association_v1"
CLONES = ROOT / "04_subclones_tymp/cnv_subclones_v1"
EXTENSION = ROOT / "04_subclones_tymp/cnv_subclone_threshold_extension_v1"
CNMF = ROOT / "05_cnmf/full_ranks4_8_v1/consensus_usage_rank7_v1.tsv.gz"
OUT = ROOT / "04_subclones_tymp/cnv_clone_tymp_cnmf_association_verification_v1"
MEMORY_GATE = PROJECT / "check_memory_gate.sh"
THRESHOLDS = (0.225, 0.250, 0.275)
FACTORS = tuple(f"F{i}" for i in range(1, 8))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def signed(values: np.ndarray) -> tuple[float, float, float]:
    nonzero = values[~np.isclose(values, 0)]
    if len(nonzero) == 0:
        return 1.0, 1.0, 0.0
    wp = float(wilcoxon(nonzero, alternative="two-sided", zero_method="wilcox", method="auto").pvalue)
    sp = float(binomtest(int(np.sum(nonzero > 0)), len(nonzero), 0.5, alternative="two-sided").pvalue)
    return wp, sp, float(np.median(values))


def main() -> None:
    if OUT.exists() and any(OUT.iterdir()):
        raise RuntimeError(f"Refusing to overwrite non-empty association verification: {OUT}")
    OUT.mkdir(parents=True, exist_ok=True)
    checks: list[dict[str, object]] = []

    def check(name: str, passed: bool, detail: str) -> None:
        checks.append({"check": name, "status": "PASS" if passed else "FAIL", "detail": detail})
        if not passed:
            pd.DataFrame(checks).to_csv(OUT / "validation_results_v1.tsv", sep="\t", index=False)
            raise RuntimeError(f"Verification failed: {name}: {detail}")

    gate = subprocess.run([str(MEMORY_GATE), "cnv_clone_tymp_cnmf_association_verification_v1"], text=True, capture_output=True)
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
    check("patient_statistical_unit", metadata["statistical_unit"] == "patient", metadata["statistical_unit"])
    check("preselection_warning", any("preselected" in x for x in metadata["interpretation_limits"]), str(metadata["interpretation_limits"]))
    check("causality_guard", any("causal" in x for x in metadata["interpretation_limits"]), str(metadata["interpretation_limits"]))

    cells = pd.read_csv(CLONES / "formal_cell_cnv_clone_assignments_v1.tsv.gz", sep="\t", dtype={"exact_cell_id": str})
    usage = pd.read_csv(CNMF, sep="\t", dtype={"exact_cell_id": str})
    grid = pd.read_csv(EXTENSION / "group_to_clone_extended_grid_v1.tsv.gz", sep="\t")
    declared_clones = pd.read_csv(SOURCE / "primary_clone_tymp_cnmf_summary_v1.tsv.gz", sep="\t")
    declared_contrasts = pd.read_csv(SOURCE / "patient_factor_high_clone_tymp_contrasts_v1.tsv", sep="\t")
    declared_summary = pd.read_csv(SOURCE / "cohort_factor_contrast_summary_v1.tsv", sep="\t")
    declared_heterogeneity = pd.read_csv(SOURCE / "patient_clone_program_heterogeneity_v1.tsv", sep="\t")
    base = cells[["exact_cell_id", "study", "infercnv_unit", "infercnv_group", "TYMP_detected_v1", "TYMP_log1p_CP10K_v1"]].merge(usage, on="exact_cell_id", validate="one_to_one")
    check("input_cell_identity", len(base) == 11950 and not base["exact_cell_id"].duplicated().any(), str(len(base)))
    check("usage_simplex", np.allclose(base[list(FACTORS)].sum(axis=1), 1.0, atol=1e-6, rtol=0), str(base[list(FACTORS)].sum(axis=1).agg(["min", "max"]).to_dict()))

    recomputed_clone_rows = []
    recomputed_contrasts = []
    expected_patient_thresholds = 0
    for threshold in THRESHOLDS:
        mapping = grid.loc[np.isclose(grid["distance_threshold"], threshold), ["infercnv_unit", "infercnv_group", "clone_id", "clone_cells"]]
        frame = base.merge(mapping, on=["infercnv_unit", "infercnv_group"], validate="many_to_one")
        primary_count = mapping.loc[mapping["clone_cells"].ge(20)].groupby("infercnv_unit")["clone_id"].nunique()
        eligible_units = set(primary_count.loc[primary_count >= 2].index)
        eligible = frame.loc[frame["clone_cells"].ge(20) & frame["infercnv_unit"].isin(eligible_units)]
        expected_patient_thresholds += len(eligible_units)
        check(f"evaluable_patients::{threshold}", len(eligible_units) == 14, str(len(eligible_units)))
        if np.isclose(threshold, 0.25):
            check("formal_primary_cells", len(eligible) <= 11604 and len(eligible) > 0, str(len(eligible)))
        for (study, unit, clone_id), clone in eligible.groupby(["study", "infercnv_unit", "clone_id"], sort=True, observed=True):
            row = {"distance_threshold": threshold, "study": study, "infercnv_unit": unit, "clone_id": clone_id, "n_cells": len(clone), "TYMP_median_log1p_CP10K": float(clone["TYMP_log1p_CP10K_v1"].median())}
            row.update({f"{factor}_median_usage": float(clone[factor].median()) for factor in FACTORS})
            recomputed_clone_rows.append(row)
    recomputed_clones = pd.DataFrame(recomputed_clone_rows)
    key = ["distance_threshold", "study", "infercnv_unit", "clone_id"]
    compared_clones = recomputed_clones.merge(declared_clones, on=key, suffixes=("_recomputed", "_declared"), validate="one_to_one")
    check("clone_summary_identity", len(compared_clones) == len(recomputed_clones) == len(declared_clones), str(len(compared_clones)))
    check("clone_cell_counts", (compared_clones["n_cells_recomputed"] == compared_clones["n_cells_declared"]).all(), str(compared_clones[["n_cells_recomputed", "n_cells_declared"]].sum().to_dict()))
    for column in ["TYMP_median_log1p_CP10K", *[f"{factor}_median_usage" for factor in FACTORS]]:
        check(f"clone_metric::{column}", np.allclose(compared_clones[f"{column}_recomputed"], compared_clones[f"{column}_declared"], atol=1e-12, rtol=0), column)

    for (threshold, study, unit), patient_clones in recomputed_clones.groupby(["distance_threshold", "study", "infercnv_unit"], sort=True, observed=True):
        for factor in FACTORS:
            ordered = patient_clones.sort_values([f"{factor}_median_usage", "clone_id"], ascending=[False, True], kind="stable")
            delta = float(ordered.iloc[0]["TYMP_median_log1p_CP10K"] - ordered.iloc[1:]["TYMP_median_log1p_CP10K"].mean())
            recomputed_contrasts.append({"distance_threshold": threshold, "study": study, "infercnv_unit": unit, "factor": factor, "factor_high_clone_id": ordered.iloc[0]["clone_id"], "patient_TYMP_delta_factor_high_vs_other_clones": delta})
    recomputed_contrasts = pd.DataFrame(recomputed_contrasts)
    contrast_key = ["distance_threshold", "study", "infercnv_unit", "factor"]
    compared = recomputed_contrasts.merge(declared_contrasts, on=contrast_key, suffixes=("_recomputed", "_declared"), validate="one_to_one")
    check("contrast_rows", len(compared) == expected_patient_thresholds * len(FACTORS), str(len(compared)))
    check("factor_high_clone", (compared["factor_high_clone_id_recomputed"] == compared["factor_high_clone_id_declared"]).all(), "exact")
    check("patient_deltas", np.allclose(compared["patient_TYMP_delta_factor_high_vs_other_clones_recomputed"], compared["patient_TYMP_delta_factor_high_vs_other_clones_declared"], atol=1e-12, rtol=0), str(np.max(np.abs(compared["patient_TYMP_delta_factor_high_vs_other_clones_recomputed"] - compared["patient_TYMP_delta_factor_high_vs_other_clones_declared"]))))

    for row in declared_summary.itertuples(index=False):
        values = recomputed_contrasts.loc[np.isclose(recomputed_contrasts["distance_threshold"], row.distance_threshold) & recomputed_contrasts["factor"].eq(row.factor), "patient_TYMP_delta_factor_high_vs_other_clones"].to_numpy(float)
        wp, sp, median = signed(values)
        check(f"summary_n::{row.distance_threshold}::{row.factor}", len(values) == row.n_patients, str(len(values)))
        check(f"summary_median::{row.distance_threshold}::{row.factor}", abs(median - row.median_patient_tymp_delta) <= 1e-12, f"recomputed={median}")
        check(f"summary_wilcoxon::{row.distance_threshold}::{row.factor}", abs(wp - row.wilcoxon_two_sided_p) <= 1e-12, f"recomputed={wp}")
        check(f"summary_sign::{row.distance_threshold}::{row.factor}", abs(sp - row.sign_test_two_sided_p) <= 1e-12, f"recomputed={sp}")

    check("heterogeneity_patient_thresholds", len(declared_heterogeneity) == expected_patient_thresholds, str(len(declared_heterogeneity)))
    check("formal_F1_rows", len(declared_contrasts.loc[np.isclose(declared_contrasts["distance_threshold"], 0.25) & declared_contrasts["factor"].eq("F1")]) == 14, "n=14")
    figures = [
        "factor_contrast_threshold_sensitivity_v1.png",
        "selected_F1_patient_tymp_deltas_v1.png",
        "selected_F1_threshold_patient_sensitivity_v1.png",
        "formal_primary_clone_tymp_landscape_v1.png",
    ]
    for name in figures:
        image = mpimg.imread(SOURCE / name)
        check(f"figure_shape::{name}", image.ndim in (2, 3) and image.shape[0] >= 800 and image.shape[1] >= 1000, str(image.shape))
        check(f"figure_finite::{name}", np.isfinite(image).all(), name)

    results = pd.DataFrame(checks)
    results.to_csv(OUT / "validation_results_v1.tsv", sep="\t", index=False)
    summary = {
        "schema_version": "part1_cnv_clone_tymp_cnmf_association_verification_v1",
        "completed_at": datetime.now().astimezone().isoformat(),
        "n_checks": int(len(results)),
        "n_pass": int(results["status"].eq("PASS").sum()),
        "n_fail": int(results["status"].eq("FAIL").sum()),
        "n_patient_thresholds": int(expected_patient_thresholds),
        "n_patient_factor_contrasts": int(len(recomputed_contrasts)),
        "manual_figure_review_pending": True,
        "dependency_changes": False,
    }
    (OUT / "verification_summary_v1.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    files = sorted(path for path in OUT.iterdir() if path.name not in {"file_manifest_v1.tsv", "RUN_COMPLETE"})
    pd.DataFrame({"file": [p.name for p in files], "size_bytes": [p.stat().st_size for p in files], "sha256": [sha256(p) for p in files]}).to_csv(OUT / "file_manifest_v1.tsv", sep="\t", index=False)
    (OUT / "RUN_COMPLETE").write_text(datetime.now().astimezone().isoformat() + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))
    print("PART1_CNV_CLONE_TYMP_CNMF_ASSOCIATION_VERIFICATION_V1_PASS")


if __name__ == "__main__":
    main()
