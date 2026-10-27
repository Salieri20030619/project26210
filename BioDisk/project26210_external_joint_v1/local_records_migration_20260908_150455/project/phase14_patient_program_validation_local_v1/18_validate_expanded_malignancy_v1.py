#!/usr/bin/env python3
"""Independent structural and inference validation for expanded malignancy analysis."""

from __future__ import annotations

import hashlib
import json
import os
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
RUN = HERE / "runs" / "20260812_expanded_malignancy_v2"
ANALYSIS = HERE / "patient_analysis_expanded_malignancy_v1"
DIAG = HERE / "posthoc_expanded_malignancy_diagnostic_v1"
OUT = HERE / "posthoc_expanded_malignancy_validation_v1"
DATASETS = ["GSE178481", "GSE202374", "GSE207493", "GSE264586", "GSE294109"]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_tsv(frame: pd.DataFrame, path: Path) -> None:
    part = path.with_name(path.name + ".part")
    frame.to_csv(part, sep="\t", index=False, encoding="utf-8", lineterminator="\n")
    os.replace(part, path)


def atomic_json(payload: dict, path: Path) -> None:
    part = path.with_name(path.name + ".part")
    part.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    os.replace(part, path)


def main() -> None:
    if OUT.exists():
        raise FileExistsError(OUT)
    OUT.mkdir()
    rows: list[dict[str, str]] = []

    def check(name: str, passed: bool, detail: str) -> None:
        rows.append({"check": name, "result": "PASS" if passed else "FAIL", "detail": detail})

    required = [
        RUN / "EXPANDED_MALIGNANCY_INPUT_STATUS_v1.tsv",
        ANALYSIS / "PATIENT_ANALYSIS_STATUS_expanded_v1.json",
        ANALYSIS / "within_dataset_patient_associations_expanded_v1.tsv",
        ANALYSIS / "FILE_MANIFEST_SHA256_expanded_v1.tsv",
        DIAG / "EXPANDED_MALIGNANCY_DIAGNOSTIC_STATUS_v1.json",
        DIAG / "strict_vs_expanded_malignant_associations_v1.tsv",
        DIAG / "patient_program_strict_vs_expanded_concordance_v1.tsv",
    ]
    for dataset in DATASETS:
        required.extend([
            RUN / dataset / "tables" / "patient_program_summary_v1.tsv",
            RUN / dataset / "tables" / "sample_program_summary_v1.tsv",
            RUN / dataset / "tables" / "expanded_malignancy_evidence_counts_v1.tsv",
            RUN / dataset / "status" / f"{dataset}_INDEPENDENT_VALIDATION_v1.tsv",
        ])
    check("required_outputs_nonempty", all(p.is_file() and p.stat().st_size > 0 for p in required), f"n={len(required)}")

    run_status = pd.read_csv(RUN / "EXPANDED_MALIGNANCY_INPUT_STATUS_v1.tsv", sep="\t")
    counts_ok = (
        len(run_status) == 5
        and set(run_status["dataset"]) == set(DATASETS)
        and (run_status["status"] == "PASS_POSTHOC_EXPANDED_MALIGNANCY").all()
        and int(run_status["tumor_like_cells"].sum()) == 40302
        and int(run_status["high_confidence_cells"].sum()) == 16645
        and int(run_status["supportive_cells"].sum()) == 17775
        and int(run_status["expanded_malignant_cells"].sum()) == 34420
        and (run_status["high_confidence_cells"] + run_status["supportive_cells"] == run_status["expanded_malignant_cells"]).all()
    )
    check("expanded_cell_counts_and_nesting", bool(counts_ok), "tumor_like=40302 high=16645 supportive=17775 expanded=34420")

    validations_ok = True
    for dataset in DATASETS:
        frame = pd.read_csv(RUN / dataset / "status" / f"{dataset}_INDEPENDENT_VALIDATION_v1.tsv", sep="\t")
        validations_ok &= bool(((frame["check"] == "OVERALL") & (frame["result"] == "PASS")).any())
        validations_ok &= bool(((frame["check"] == "part1_rebuild_priority") & (frame["result"] == "PASS")).any())
        validations_ok &= bool(((frame["check"] == "tier_nesting") & (frame["result"] == "PASS")).any())
    check("five_dataset_rebuild_priority_and_nesting", validations_ok, "5/5 dataset validators PASS")

    manifest = pd.read_csv(ANALYSIS / "FILE_MANIFEST_SHA256_expanded_v1.tsv", sep="\t")
    failures = []
    for row in manifest.itertuples(index=False):
        path = ANALYSIS / row.file
        if not path.is_file() or path.stat().st_size != int(row.bytes) or sha256(path) != row.sha256:
            failures.append(row.file)
    check("analysis_output_manifest_hashes", not failures, ";".join(failures) if failures else f"n={len(manifest)}")

    status = json.loads((ANALYSIS / "PATIENT_ANALYSIS_STATUS_expanded_v1.json").read_text(encoding="utf-8"))
    status_ok = (
        status.get("status") == "PASS_POSTHOC_EXPANDED_MALIGNANCY"
        and status.get("malignancy_gate_label") == "high_confidence_plus_supportive_rebuild_priority"
        and status.get("association_rows") == 120
        and status.get("association_rows_tested") == 84
        and status.get("minimum_cells") == 20
        and status.get("minimum_patients") == 5
        and status.get("bootstrap_replicates") == 9999
    )
    check("analysis_status_and_frozen_parameters", status_ok, "rows=120 tested=84 min_cells=20 min_patients=5 boot=9999")

    assoc = pd.read_csv(ANALYSIS / "within_dataset_patient_associations_expanded_v1.tsv", sep="\t")
    tested = assoc["spearman_rho"].notna()
    nominal = int((assoc.loc[tested, "p_value"] < 0.05).sum())
    fdr = int((assoc.loc[tested, "p_fdr_bh_within_dataset_method"] < 0.05).sum())
    ci_excludes = int(((assoc.loc[tested, "bootstrap_ci_low"] > 0) | (assoc.loc[tested, "bootstrap_ci_high"] < 0)).sum())
    grid_ok = len(assoc) == 120 and set(assoc["dataset"]) == set(DATASETS) and int(tested.sum()) == 84 and nominal == 4 and fdr == 0 and ci_excludes == 4
    check("association_grid_and_inference_counts", grid_ok, f"rows={len(assoc)} tested={tested.sum()} nominal={nominal} fdr={fdr} ci_excludes={ci_excludes}")

    diag = json.loads((DIAG / "EXPANDED_MALIGNANCY_DIAGNOSTIC_STATUS_v1.json").read_text(encoding="utf-8"))
    diag_ok = (
        diag.get("status") == "PASS_POSTHOC_EXPANDED_MALIGNANCY"
        and diag.get("high_confidence_cells") == 16645
        and diag.get("supportive_cells") == 17775
        and diag.get("expanded_cells") == 34420
        and diag.get("association_rows_tested") == 84
        and diag.get("malignant_rows_tested") == 54
        and diag.get("malignant_rows_gained_testability") == 18
        and diag.get("formal_strict_v6_replaced") is False
    )
    check("diagnostic_status_reconciled", diag_ok, json.dumps(diag, ensure_ascii=False, sort_keys=True))

    strict_assoc = pd.read_csv(HERE / "patient_analysis_final_cnv_v6" / "within_dataset_patient_associations_v6.tsv", sep="\t")
    strict_fdr = int((strict_assoc["p_fdr_bh_within_dataset_method"] < 0.05).sum())
    check("formal_strict_v6_unchanged", strict_fdr == 0, f"strict_v6_fdr={strict_fdr}; expanded_is_posthoc=true")

    all_pass = all(row["result"] == "PASS" for row in rows)
    completed = datetime.now().astimezone().isoformat(timespec="seconds")
    rows.append({"check": "OVERALL", "result": "PASS" if all_pass else "FAIL", "detail": completed})
    atomic_tsv(pd.DataFrame(rows), OUT / "EXPANDED_MALIGNANCY_VALIDATION_STATUS_v1.tsv")
    payload = {
        "schema_version": "phase14_expanded_malignancy_validation_v1",
        "status": "PASS" if all_pass else "FAIL",
        "completed_at": completed,
        "checks_passed": sum(row["result"] == "PASS" for row in rows),
        "checks_total": len(rows),
        "formal_strict_v6_replaced": False,
    }
    atomic_json(payload, OUT / "EXPANDED_MALIGNANCY_VALIDATION_SUMMARY_v1.json")
    print(json.dumps(payload, ensure_ascii=False, indent=2))
    if not all_pass:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
