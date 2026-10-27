#!/usr/bin/env python3
"""Independent consistency checks for the Phase 14 zero-FDR diagnostic."""

from __future__ import annotations

import hashlib
import json
import os
import sys
from datetime import datetime
from pathlib import Path
from zoneinfo import ZoneInfo

import pandas as pd


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent
OUT = HERE / "posthoc_zero_fdr_validation_v4"

FORMAL = HERE / "patient_analysis_final_cnv_v6"
DIAG = HERE / "posthoc_zero_fdr_diagnostic_v1"
ATTR = HERE / "posthoc_cnv_gate_attrition_v2"
TIERS = HERE / "posthoc_malignancy_gate_tiers_v3"
FDR = HERE / "posthoc_fdr_family_and_key_patients_v1"
SUMMARY = HERE / "ZERO_FDR_DIAGNOSTIC_SUMMARY_v1.md"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_atomic(path: Path, text: str) -> None:
    part = path.with_suffix(path.suffix + ".part")
    part.write_text(text, encoding="utf-8", newline="")
    os.replace(part, path)


def add_check(rows: list[dict[str, str]], name: str, passed: bool, detail: str) -> None:
    rows.append({"check": name, "result": "PASS" if passed else "FAIL", "detail": detail})


def main() -> int:
    if OUT.exists():
        raise RuntimeError(f"Refusing to overwrite existing validation directory: {OUT}")
    OUT.mkdir(parents=False)

    required = [
        FORMAL / "PATIENT_ANALYSIS_STATUS_v6.json",
        FORMAL / "INDEPENDENT_VALIDATION_STATUS_v6.tsv",
        FORMAL / "within_dataset_patient_associations_v6.tsv",
        DIAG / "DIAGNOSTIC_STATUS_v1.json",
        DIAG / "all_tested_associations_ranked_v1.tsv",
        DIAG / "malignant_retention_summary_v1.tsv",
        DIAG / "exploratory_cross_study_correlation_meta_v1.tsv",
        ATTR / "CNV_GATE_ATTRITION_STATUS_v1.json",
        ATTR / "patient_gate_summary_v1.tsv",
        TIERS / "GATE_TIER_SENSITIVITY_STATUS_v1.json",
        TIERS / "gate_tier_summary_v1.tsv",
        TIERS / "malignant_relations_by_gate_tier_v1.tsv",
        FDR / "FDR_FAMILY_AUDIT_STATUS_v1.json",
        FDR / "planned_family_bh_audit_v1.tsv",
        SUMMARY,
    ]
    checks: list[dict[str, str]] = []
    add_check(
        checks,
        "required_outputs_nonempty",
        all(path.is_file() and path.stat().st_size > 0 for path in required),
        f"n={len(required)}",
    )

    formal_status = json.loads((FORMAL / "PATIENT_ANALYSIS_STATUS_v6.json").read_text(encoding="utf-8"))
    formal_validation = pd.read_csv(FORMAL / "INDEPENDENT_VALIDATION_STATUS_v6.tsv", sep="\t")
    associations = pd.read_csv(FORMAL / "within_dataset_patient_associations_v6.tsv", sep="\t")
    diagnostic = json.loads((DIAG / "DIAGNOSTIC_STATUS_v1.json").read_text(encoding="utf-8"))
    attrition = json.loads((ATTR / "CNV_GATE_ATTRITION_STATUS_v1.json").read_text(encoding="utf-8"))
    tier_status = json.loads((TIERS / "GATE_TIER_SENSITIVITY_STATUS_v1.json").read_text(encoding="utf-8"))
    tier_summary = pd.read_csv(TIERS / "gate_tier_summary_v1.tsv", sep="\t")
    fdr_status = json.loads((FDR / "FDR_FAMILY_AUDIT_STATUS_v1.json").read_text(encoding="utf-8"))
    retention = pd.read_csv(DIAG / "malignant_retention_summary_v1.tsv", sep="\t")

    add_check(
        checks,
        "formal_v6_passed_independent_validation",
        formal_status.get("status") == "PASS"
        and bool(formal_status.get("final_malignancy_gate"))
        and ((formal_validation["check"] == "OVERALL") & (formal_validation["result"] == "PASS")).any(),
        f"status={formal_status.get('status')} gate={formal_status.get('final_malignancy_gate')}",
    )
    tested = associations[associations["status"] == "TESTED"].copy()
    formal_fdr = int((pd.to_numeric(tested["p_fdr_bh_within_dataset_method"], errors="coerce") < 0.05).sum())
    nominal = int((pd.to_numeric(tested["p_value"], errors="coerce") < 0.05).sum())
    ci_excludes = int(
        (
            (pd.to_numeric(tested["bootstrap_ci_low"], errors="coerce") > 0)
            | (pd.to_numeric(tested["bootstrap_ci_high"], errors="coerce") < 0)
        ).sum()
    )
    add_check(
        checks,
        "formal_association_counts_reconciled",
        len(associations) == 120 and len(tested) == 66 and nominal == 5 and ci_excludes == 2 and formal_fdr == 0,
        f"rows={len(associations)} tested={len(tested)} nominal={nominal} ci_excludes={ci_excludes} fdr={formal_fdr}",
    )
    add_check(
        checks,
        "diagnostic_status_reconciled",
        diagnostic.get("status") == "PASS_POSTHOC_DIAGNOSTIC"
        and diagnostic.get("tested_associations") == len(tested)
        and diagnostic.get("nominal_p_lt_0_05") == nominal
        and diagnostic.get("bootstrap_ci_excludes_zero") == ci_excludes
        and diagnostic.get("lost_testability_after_final_cnv_gate") == 36,
        json.dumps(diagnostic, ensure_ascii=False, sort_keys=True),
    )
    add_check(
        checks,
        "cnv_gate_cell_counts_reconciled",
        attrition.get("status") == "PASS_POSTHOC_DIAGNOSTIC"
        and attrition.get("tumor_like_rows") == 40302
        and attrition.get("final_supported_rows") == 16645
        and attrition.get("tumor_like_in_external_scope") == 40118
        and attrition.get("tumor_like_not_in_external_scope") == 184
        and int(retention["provisional_cells"].sum()) == 40302
        and int(retention["final_cells"].sum()) == 16645,
        "phenotype=40302 final=16645 external_scope=40118 outside=184",
    )

    tier_map = tier_summary.set_index("gate_tier").to_dict(orient="index")
    expected_tiers = {
        "phenotype_candidate": (72, 6, 2),
        "any_method_positive": (54, 4, 0),
        "strict_conservative": (36, 4, 0),
    }
    tier_ok = tier_status.get("status") == "PASS_POSTHOC_SENSITIVITY"
    for name, values in expected_tiers.items():
        row = tier_map.get(name, {})
        tier_ok = tier_ok and tuple(int(row.get(field, -1)) for field in ("rows_tested", "nominal_p_lt_0_05", "fdr_q_lt_0_05")) == values
    add_check(checks, "gate_tier_sensitivity_reconciled", tier_ok, json.dumps(expected_tiers, sort_keys=True))

    add_check(
        checks,
        "planned_fdr_family_audit_reconciled",
        fdr_status.get("status") == "PASS_POSTHOC_AUDIT"
        and fdr_status.get("reported_fdr_significant") == 0
        and fdr_status.get("planned_family_fdr_significant") == 0
        and fdr_status.get("planned_family_sensitivity_uses_12_hypotheses") is True,
        "reported=0 planned_12_with_untested_as_1=0",
    )

    summary_text = " ".join(SUMMARY.read_text(encoding="utf-8").split())
    precedence_terms = [
        "part1_rebuild_v1",
        "Original Part 1 clustering and annotation are not used",
        "relabel cells",
        "do not reverse rebuild priority",
        "no Phase 14 v6 result passes its prespecified FDR gate",
    ]
    add_check(
        checks,
        "summary_declares_rebuild_priority_and_unchanged_formal_result",
        all(term in summary_text for term in precedence_terms),
        ";".join(precedence_terms),
    )

    all_pass = all(row["result"] == "PASS" for row in checks)
    completed = datetime.now(ZoneInfo("Asia/Taipei")).isoformat(timespec="seconds")
    checks.append({"check": "OVERALL", "result": "PASS" if all_pass else "FAIL", "detail": completed})
    checks_df = pd.DataFrame(checks)
    write_atomic(OUT / "ZERO_FDR_DIAGNOSTIC_VALIDATION_STATUS_v1.tsv", checks_df.to_csv(sep="\t", index=False, lineterminator="\n"))

    manifest_rows = []
    for path in required:
        manifest_rows.append(
            {
                "relative_path": path.relative_to(PROJECT).as_posix(),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
        )
    manifest_df = pd.DataFrame(manifest_rows)
    write_atomic(OUT / "ZERO_FDR_DIAGNOSTIC_FILE_MANIFEST_v1.tsv", manifest_df.to_csv(sep="\t", index=False, lineterminator="\n"))

    status = {
        "schema_version": "phase14_zero_fdr_diagnostic_validation_v4",
        "status": "PASS" if all_pass else "FAIL",
        "completed_at": completed,
        "checks_passed": sum(row["result"] == "PASS" for row in checks),
        "checks_total": len(checks),
        "formal_result_unchanged": formal_fdr == 0,
        "cnv_precedence": "part1_rebuild_v1 exact-cell calls retain priority",
        "original_part1_clustering_annotation_used": False,
        "python": sys.version,
        "pandas": pd.__version__,
    }
    write_atomic(OUT / "ZERO_FDR_DIAGNOSTIC_VALIDATION_SUMMARY_v1.json", json.dumps(status, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps(status, indent=2, ensure_ascii=False))
    return 0 if all_pass else 1


if __name__ == "__main__":
    raise SystemExit(main())
