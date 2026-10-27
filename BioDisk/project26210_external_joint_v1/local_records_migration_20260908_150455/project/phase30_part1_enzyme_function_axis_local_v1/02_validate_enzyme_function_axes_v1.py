#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def add(checks: list[dict], name: str, passed: bool, detail: str) -> None:
    checks.append({"check": name, "status": "PASS" if passed else "FAIL", "detail": detail})


def main() -> None:
    root = Path(__file__).resolve().parent
    analysis = root / "results_local_v1"
    output = root / "validation_local_v1"
    part = root / "validation_local_v1.part"
    if output.exists() or part.exists():
        raise FileExistsError("refusing to overwrite validation output")
    if not analysis.is_dir():
        raise FileNotFoundError(analysis)
    part.mkdir(parents=True)
    checks: list[dict] = []

    required = [
        "input_manifest_v1.tsv", "enzyme_axis_reaction_registry_v1.tsv",
        "enzyme_axis_component_qc_v1.tsv", "enzyme_axis_feature_registry_v1.tsv",
        "enzyme_axis_cell_scores_v1.tsv.gz", "enzyme_axis_patient_effects_v1.tsv",
        "enzyme_axis_endpoint_summary_v1.tsv", "enzyme_axis_integrated_evidence_v1.tsv",
        "enzyme_axis_joint_terminal_summary_v1.tsv", "status.json",
    ]
    missing = [x for x in required if not (analysis / x).is_file()]
    add(checks, "required_outputs", not missing, f"missing={missing}")

    manifest = pd.read_csv(analysis / "input_manifest_v1.tsv", sep="\t")
    hashes_ok = len(manifest) == 6
    for row in manifest.itertuples(index=False):
        p = root / "input_snapshot_local_v1" / row.file
        hashes_ok &= p.is_file() and sha256(p) == row.sha256
    add(checks, "input_hashes", bool(hashes_ok), f"n_inputs={len(manifest)}")

    reactions = pd.read_csv(analysis / "enzyme_axis_reaction_registry_v1.tsv", sep="\t")
    expected = {"M_45", "M_164", "M_157", "M_138", "M_140"}
    reaction_ok = set(reactions["module"]) == expected
    reaction_ok &= "TYMP" in reactions.loc[reactions["module"] == "M_45", "module_genes"].iloc[0].split(";")
    reaction_ok &= "RRM2B" in reactions.loc[reactions["module"] == "M_157", "module_genes"].iloc[0].split(";")
    add(checks, "reaction_mapping", bool(reaction_ok), f"modules={sorted(reactions['module'].tolist())}")

    scores = pd.read_csv(analysis / "enzyme_axis_cell_scores_v1.tsv.gz", sep="\t")
    cell_ok = len(scores) == 11950 and scores["exact_cell_id"].nunique() == 11950
    add(checks, "cell_identity", cell_ok, f"rows={len(scores)}; unique={scores['exact_cell_id'].nunique()}")
    patient_count = scores["curated_patient_id"].nunique()
    add(checks, "patient_registry", patient_count == 26, f"n_patients_total={patient_count}")

    required_numeric = [c for c in scores if c.startswith("z_") or c.endswith("axis_primary")
                        or c.endswith("axis_broad_sensitivity")]
    finite_ok = np.isfinite(scores[required_numeric].to_numpy(dtype=float)).all()
    add(checks, "finite_axis_values", bool(finite_ok), f"n_numeric_axis_columns={len(required_numeric)}")

    formula_diffs = {
        "TYMP_axis_primary": np.max(np.abs(scores["TYMP_axis_primary"] - scores[["z_TYMP_RNA", "z_TYMP_M45_flux", "z_TYMP_M164_flux"]].mean(axis=1))),
        "TYMP_axis_broad_sensitivity": np.max(np.abs(scores["TYMP_axis_broad_sensitivity"] - scores[["z_TYMP_RNA", "z_TYMP_Pyrimidine_flux"]].mean(axis=1))),
        "RRM2B_axis_primary": np.max(np.abs(scores["RRM2B_axis_primary"] - scores[["z_RRM2B_RNA", "z_RRM_M157_flux"]].mean(axis=1))),
        "RRM2B_axis_broad_sensitivity": np.max(np.abs(scores["RRM2B_axis_broad_sensitivity"] - scores[["z_RRM2B_RNA", "z_RRM_M138_flux", "z_RRM_M140_flux", "z_RRM_M157_flux"]].mean(axis=1))),
    }
    add(checks, "axis_formulas", max(formula_diffs.values()) < 1e-12,
        json.dumps(formula_diffs, sort_keys=True))

    patient_means = scores.groupby("curated_patient_id")[required_numeric].mean().abs().to_numpy()
    centered_ok = float(np.nanmax(patient_means)) < 1e-10
    add(checks, "within_patient_centering", centered_ok,
        f"max_abs_patient_mean={float(np.nanmax(patient_means)):.3g}")

    features = pd.read_csv(analysis / "enzyme_axis_feature_registry_v1.tsv", sep="\t")
    scope_ok = len(features) == 12 and not features["is_VIPER_activity"].astype(bool).any()
    scope_ok &= not features["is_direct_biochemical_activity"].astype(bool).any()
    add(checks, "interpretation_boundary", bool(scope_ok),
        f"n_features={len(features)}; viper_rows={features['is_VIPER_activity'].astype(bool).sum()}")

    detail = pd.read_csv(analysis / "enzyme_axis_patient_effects_v1.tsv", sep="\t")
    expected_endpoints = {"F1_continuous", "F1_dominant_group", "F1_patient_quartile",
                          "trajectory_terminal_shift"}
    detail_ok = set(detail["endpoint"]) == expected_endpoints
    detail_ok &= set(detail.loc[detail["endpoint"] == "trajectory_terminal_shift", "method"]) == {"Monocle3", "Slingshot", "scTour"}
    add(checks, "endpoint_coverage", bool(detail_ok),
        f"endpoints={sorted(detail['endpoint'].unique())}")

    threshold_ok = True
    threshold_ok &= (detail.loc[detail["endpoint"] == "F1_continuous", "n_cells"] >= 20).all()
    grouped = detail[detail["endpoint"] == "F1_dominant_group"]
    threshold_ok &= (grouped["n_low"] >= 10).all() and (grouped["n_high"] >= 10).all()
    other = detail[detail["endpoint"].isin(["F1_patient_quartile", "trajectory_terminal_shift"])]
    threshold_ok &= (other["n_low"] >= 5).all() and (other["n_high"] >= 5).all()
    add(checks, "cell_count_thresholds", bool(threshold_ok), f"n_patient_effect_rows={len(detail)}")

    summary = pd.read_csv(analysis / "enzyme_axis_endpoint_summary_v1.tsv", sep="\t")
    grouped_detail = detail.groupby(["feature", "feature_class", "endpoint", "method"], dropna=False)["effect"]
    recomputed = grouped_detail.agg(n="count", median="median").reset_index()
    merged = summary.merge(recomputed, on=["feature", "feature_class", "endpoint", "method"], validate="one_to_one")
    summary_ok = len(merged) == len(summary)
    summary_ok &= (merged["n_eligible_patients"] == merged["n"]).all()
    summary_ok &= np.allclose(merged["median_patient_effect"], merged["median"], atol=1e-12)
    add(checks, "summary_recalculation", bool(summary_ok), f"n_summary_rows={len(summary)}")

    tier_ok = ((summary["n_eligible_patients"] >= 10) ==
               (summary["reporting_tier"] == "FDR_AT_LEAST_10")).all()
    tier_ok &= ((summary["n_eligible_patients"].between(3, 9)) ==
                (summary["reporting_tier"] == "NOMINAL_P_3_TO_9")).all()
    add(checks, "small_n_reporting_rule", bool(tier_ok),
        str(summary["reporting_tier"].value_counts().to_dict()))

    evidence = pd.read_csv(analysis / "enzyme_axis_integrated_evidence_v1.tsv", sep="\t")
    bool_cols = ["F1_continuous", "F1_dominant_group", "F1_patient_quartile",
                 "trajectory_convergence_2_of_3"]
    counts = evidence[bool_cols].astype(bool).sum(axis=1)
    expected_tier = np.where(counts == 4, "HIGH_SUPPORT",
                             np.where(counts >= 3, "SECONDARY_SUPPORT", "LIMITED_OR_UNSUPPORTED"))
    evidence_ok = set(evidence["axis"]) == {"TYMP_axis_primary", "RRM2B_axis_primary"}
    evidence_ok &= (counts == evidence["n_supported_endpoints_of_4"]).all()
    evidence_ok &= (expected_tier == evidence["evidence_tier"]).all()
    add(checks, "integrated_evidence_rule", bool(evidence_ok), evidence.to_json(orient="records"))

    joint = pd.read_csv(analysis / "enzyme_axis_joint_terminal_summary_v1.tsv", sep="\t")
    joint_ok = len(joint) == 3 and set(joint["method"]) == {"Monocle3", "Slingshot", "scTour"}
    joint_ok &= joint["fraction_both_positive"].between(0, 1).all()
    add(checks, "joint_terminal_summary", bool(joint_ok), f"n_methods={len(joint)}")

    analysis_status = json.loads((analysis / "status.json").read_text(encoding="utf-8"))
    status_ok = analysis_status["status"] == "PASS_ANALYSIS_COMPLETE_AWAITING_INDEPENDENT_VALIDATION"
    status_ok &= analysis_status["figures_created"] is False
    add(checks, "analysis_status", bool(status_ok), analysis_status["status"])

    figure_files = [p.name for p in analysis.rglob("*") if p.suffix.lower() in {".png", ".pdf", ".svg", ".jpg", ".jpeg"}]
    add(checks, "no_figures", not figure_files, f"figure_files={figure_files}")

    checks_frame = pd.DataFrame(checks)
    checks_frame.to_csv(part / "verification_checks_v1.tsv", sep="\t", index=False)
    passed = int((checks_frame["status"] == "PASS").sum())
    status = {
        "status": "PASS_INDEPENDENT_VALIDATION" if passed == len(checks) else "FAIL_VALIDATION",
        "checks_passed": passed,
        "checks_total": len(checks),
        "analysis_directory_sha256": {
            p.name: sha256(p) for p in sorted(analysis.iterdir()) if p.is_file()
        },
    }
    (part / "status.json").write_text(json.dumps(status, indent=2) + "\n", encoding="utf-8")
    os.replace(part, output)
    if passed != len(checks):
        raise RuntimeError(f"validation failed: {passed}/{len(checks)}")


if __name__ == "__main__":
    main()

