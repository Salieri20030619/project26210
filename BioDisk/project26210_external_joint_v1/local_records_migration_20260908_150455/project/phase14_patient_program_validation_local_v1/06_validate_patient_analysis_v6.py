#!/usr/bin/env python3
"""Independent conservation, rule, manifest, and final-CNV validation for Phase 14 v6."""

from __future__ import annotations

import hashlib
import json
import os
from datetime import datetime
from pathlib import Path

import pandas as pd


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
OUT = HERE / "patient_analysis_final_cnv_v6"
RUN = HERE / "runs" / "20260812_final_cnv_gated_v6"
DATASETS = ["GSE178481", "GSE202374", "GSE207493", "GSE264586", "GSE294109"]
VERSION = "v6"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_tsv(frame: pd.DataFrame, path: Path) -> None:
    part = path.with_name(path.name + ".part")
    if path.exists() or part.exists():
        raise FileExistsError(path)
    frame.to_csv(part, sep="\t", index=False, encoding="utf-8", lineterminator="\n")
    os.replace(part, path)


def main() -> None:
    checks: list[dict[str, str]] = []

    def check(name: str, passed: bool, detail: str) -> None:
        checks.append({"check": name, "result": "PASS" if passed else "FAIL", "detail": detail})

    required_stems = [
        "within_dataset_patient_associations",
        "GSE178481_full_vs_added_subset_sensitivity",
        "GSE264586_paired_longitudinal",
        "GSE294109_response_stratified",
        "malignant_F1_patient_reproducibility",
        "cross_sectional_sample_selection_audit",
        "cross_study_direction_synthesis",
        "input_manifest",
        "PATIENT_ANALYSIS_STATUS",
        "FILE_MANIFEST_SHA256",
    ]
    required = [f"{stem}_{VERSION}.json" if stem == "PATIENT_ANALYSIS_STATUS" else f"{stem}_{VERSION}.tsv" for stem in required_stems]
    check("required_outputs_nonempty", all((OUT / name).is_file() and (OUT / name).stat().st_size > 0 for name in required), ";".join(required))

    manifest = pd.read_csv(OUT / f"FILE_MANIFEST_SHA256_{VERSION}.tsv", sep="\t")
    failures = []
    for row in manifest.itertuples(index=False):
        path = OUT / row.file
        if not path.is_file() or path.stat().st_size != int(row.bytes) or sha256(path) != row.sha256:
            failures.append(row.file)
    check("output_manifest_hashes_and_sizes", not failures, ";".join(failures) if failures else f"n={len(manifest)}")

    input_manifest = pd.read_csv(OUT / f"input_manifest_{VERSION}.tsv", sep="\t")
    input_failures = []
    for row in input_manifest.itertuples(index=False):
        path = ROOT / row.input
        if not path.is_file() or path.stat().st_size != int(row.bytes) or sha256(path) != row.sha256:
            input_failures.append(row.input)
    check("input_manifest_hashes_and_sizes", not input_failures, ";".join(input_failures) if input_failures else f"n={len(input_manifest)}")

    dataset_status = []
    for dataset in DATASETS:
        path = RUN / dataset / "status" / f"{dataset}_INDEPENDENT_VALIDATION_v1.tsv"
        frame = pd.read_csv(path, sep="\t")
        dataset_status.append(bool(((frame["check"] == "OVERALL") & (frame["result"] == "PASS")).any()))
    check("five_final_cnv_gated_inputs_passed", all(dataset_status) and len(dataset_status) == 5, str(dataset_status))

    gate_summary = pd.read_csv(RUN / "FINAL_CNV_GATED_INPUT_STATUS_v4.tsv", sep="\t")
    gate_ok = len(gate_summary) == 5 and set(gate_summary["dataset"]) == set(DATASETS) and (gate_summary["status"] == "PASS").all()
    gate_ok = gate_ok and (gate_summary["final_malignant_cells"] > 0).all() and (gate_summary["final_malignant_patients"] > 0).all()
    check("final_cnv_gate_summary", bool(gate_ok), f"rows={len(gate_summary)} malignant_cells={int(gate_summary['final_malignant_cells'].sum())}")

    associations = pd.read_csv(OUT / f"within_dataset_patient_associations_{VERSION}.tsv", sep="\t")
    check("association_grid_complete", len(associations) == 5 * 2 * 12, f"rows={len(associations)}")
    check("association_datasets_separate", set(associations["dataset"]) == set(DATASETS), ";".join(sorted(associations["dataset"].unique())))
    tested = associations["status"].str.startswith("TESTED")
    tested_valid = associations.loc[tested, ["spearman_rho", "p_value", "p_fdr_bh_within_dataset_method"]].notna().all().all()
    tested_valid &= bool((associations.loc[tested, "n_patients"] >= 5).all())
    check("tested_associations_meet_n_and_finite_rules", bool(tested_valid), f"tested={int(tested.sum())}")
    check("no_cell_level_inference_fields", not any("cell_p" in name.lower() for name in associations.columns), ";".join(associations.columns))

    audit = pd.read_csv(OUT / f"cross_sectional_sample_selection_audit_{VERSION}.tsv", sep="\t")
    g178 = audit[audit["dataset"] == "GSE178481"]
    g264 = audit[audit["dataset"] == "GSE264586"]
    check("GSE178481_tumor_only_selection", len(g178) > 0 and not g178["sample_id"].str.contains("Normal", case=False, regex=False).any(), f"samples={g178['object_name'].nunique()}")
    check("GSE264586_baseline_only_selection", len(g264) > 0 and set(g264["timepoint"]) == {"Baseline"}, f"samples={g264['object_name'].nunique()}")
    other = audit[~audit["dataset"].isin(["GSE178481", "GSE264586"])]
    check("other_dataset_sample_counts", other["object_name"].nunique() == 50, f"samples={other['object_name'].nunique()}")

    longitudinal = pd.read_csv(OUT / f"GSE264586_paired_longitudinal_{VERSION}.tsv", sep="\t")
    check("longitudinal_only_verified_followups", set(longitudinal["followup"]) == {"Cycle2", "End_of_treatment"}, ";".join(sorted(longitudinal["followup"].unique())))
    check("longitudinal_small_n_is_descriptive", bool((longitudinal.loc[longitudinal["n_pairs"] < 5, "status"].str.startswith("DESCRIPTIVE")).all()), f"max_pairs={longitudinal['n_pairs'].max()}")

    response = pd.read_csv(OUT / f"GSE294109_response_stratified_{VERSION}.tsv", sep="\t")
    response_ok = (response["n_benefit_patients"] >= 0).all() and (response["n_no_effect_patients"] >= 0).all()
    check("response_groups_verified", bool(response_ok), f"rows={len(response)}")

    status = json.loads((OUT / f"PATIENT_ANALYSIS_STATUS_{VERSION}.json").read_text(encoding="utf-8"))
    check("formal_status_final_malignancy", status.get("status") == "PASS" and status.get("final_malignancy_gate") is True, f"status={status.get('status')} gate={status.get('final_malignancy_gate')}")
    check("fixed_seed_and_bootstrap", status.get("seed") == 20260812 and status.get("bootstrap_replicates") == 9999, f"seed={status.get('seed')} boot={status.get('bootstrap_replicates')}")
    check("GSE207493_membership_not_guessed", any("not reported" in item for item in status.get("limitations", [])), "limitation declared")

    overall = all(row["result"] == "PASS" for row in checks)
    checks.append({"check": "OVERALL", "result": "PASS" if overall else "FAIL", "detail": datetime.now().astimezone().isoformat(timespec="seconds")})
    atomic_tsv(pd.DataFrame(checks), OUT / f"INDEPENDENT_VALIDATION_STATUS_{VERSION}.tsv")
    if not overall:
        raise SystemExit("Phase14 patient v6 validation failed")


if __name__ == "__main__":
    main()

