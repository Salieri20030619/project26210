#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path

import pandas as pd


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def add(rows: list[dict], name: str, passed: bool, detail: str) -> None:
    rows.append({"check": name, "status": "PASS" if passed else "FAIL", "detail": detail})


def main() -> None:
    root = Path(__file__).resolve().parent
    output = root / "terminal_validation_local_v1"
    part = root / "terminal_validation_local_v1.part"
    if output.exists() or part.exists():
        raise FileExistsError("refusing to overwrite Phase 30 terminal validation")
    part.mkdir(parents=True)
    rows: list[dict] = []

    independent = json.loads((root / "validation_local_v1" / "status.json").read_text(encoding="utf-8"))
    add(rows, "independent_validation", independent.get("status") == "PASS_INDEPENDENT_VALIDATION"
        and independent.get("checks_passed") == 17 and independent.get("checks_total") == 17,
        f"{independent.get('checks_passed')}/{independent.get('checks_total')}")

    analysis_hashes_ok = True
    for name, expected in independent["analysis_directory_sha256"].items():
        path = root / "results_local_v1" / name
        analysis_hashes_ok &= path.is_file() and sha256(path) == expected
    add(rows, "validated_analysis_hashes", bool(analysis_hashes_ok),
        f"n_files={len(independent['analysis_directory_sha256'])}")

    integration_status = json.loads((root / "all_axis_integration_local_v2" / "status.json").read_text(encoding="utf-8"))
    integration = pd.read_csv(root / "all_axis_integration_local_v2" /
                              "enzyme_axis_integrated_evidence_all_axes_v2.tsv", sep="\t")
    add(rows, "all_axis_integration_status",
        integration_status.get("status") == "PASS_ALL_AXIS_INTEGRATION" and len(integration) == 4,
        f"status={integration_status.get('status')}; n_axes={len(integration)}")

    expected = {
        "TYMP_axis_primary": (0, "LIMITED_OR_UNSUPPORTED", 0, False),
        "RRM2B_axis_primary": (0, "LIMITED_OR_UNSUPPORTED", 0, False),
        "TYMP_axis_broad_sensitivity": (3, "SECONDARY_SUPPORT", 0, False),
        "RRM2B_axis_broad_sensitivity": (0, "LIMITED_OR_UNSUPPORTED", 3, True),
    }
    integration_ok = set(integration["axis"]) == set(expected)
    for row in integration.itertuples(index=False):
        e = expected[row.axis]
        integration_ok &= (row.n_positive_endpoints_of_4, row.positive_evidence_tier,
                           row.n_negative_endpoints_of_4, bool(row.opposing_evidence)) == e
    add(rows, "all_axis_integrated_results", bool(integration_ok),
        integration[["axis", "n_positive_endpoints_of_4", "positive_evidence_tier",
                     "n_negative_endpoints_of_4", "opposing_evidence"]].to_json(orient="records"))

    joint = pd.read_csv(root / "results_local_v1" / "enzyme_axis_joint_terminal_summary_v1.tsv", sep="\t")
    joint_ok = joint.set_index("method")["n_both_positive"].to_dict() == {
        "Monocle3": 1, "Slingshot": 4, "scTour": 5
    }
    add(rows, "joint_terminal_counts", bool(joint_ok),
        joint.set_index("method")["n_both_positive"].to_dict().__str__())

    record = (root / "RUN_RECORD_v1.md").read_text(encoding="utf-8")
    record_ok = all(x in record for x in ["SECONDARY_SUPPORT", "not VIPER activity",
                                           "not direct biochemical enzyme assays", "17/17"])
    add(rows, "run_record_scope", record_ok, f"characters={len(record)}")

    register = pd.read_csv(root / "PHASE30_COMPLETION_REGISTER_v1.tsv", sep="\t")
    register_ok = len(register) == 1 and register.iloc[0]["phase30_status"] == "COMPLETE_VALIDATED"
    add(rows, "completion_register", bool(register_ok), f"rows={len(register)}")

    failed_preserved = (root / "failed_attempt_v1" / "input_manifest_v1.tsv").is_file()
    add(rows, "failed_attempt_preserved", failed_preserved, "failed_attempt_v1/input_manifest_v1.tsv")

    figures = [str(p.relative_to(root)) for p in root.rglob("*")
               if p.suffix.lower() in {".png", ".pdf", ".svg", ".jpg", ".jpeg"}]
    add(rows, "no_figures", not figures, f"figures={figures}")

    frame = pd.DataFrame(rows)
    frame.to_csv(part / "terminal_checks_v1.tsv", sep="\t", index=False)
    passed = int((frame["status"] == "PASS").sum())
    files = [
        root / "PHASE30_ENZYME_AXIS_REGISTRY_v1.md",
        root / "RUN_RECORD_v1.md",
        root / "PHASE30_COMPLETION_REGISTER_v1.tsv",
        root / "results_local_v1" / "enzyme_axis_endpoint_summary_v1.tsv",
        root / "results_local_v1" / "enzyme_axis_joint_terminal_summary_v1.tsv",
        root / "validation_local_v1" / "status.json",
        root / "all_axis_integration_local_v2" / "enzyme_axis_integrated_evidence_all_axes_v2.tsv",
        root / "all_axis_integration_local_v2" / "status.json",
    ]
    pd.DataFrame([{"file": str(p.relative_to(root)), "size_bytes": p.stat().st_size,
                   "sha256": sha256(p)} for p in files]).to_csv(
        part / "PHASE30_VALIDATION_MANIFEST_v1.tsv", sep="\t", index=False)
    status = {"status": "PASS_PHASE30_TERMINAL_VALIDATION" if passed == len(rows) else "FAIL",
              "checks_passed": passed, "checks_total": len(rows)}
    (part / "status.json").write_text(json.dumps(status, indent=2) + "\n", encoding="utf-8")
    os.replace(part, output)
    if passed != len(rows):
        raise RuntimeError(f"terminal validation failed: {passed}/{len(rows)}")


if __name__ == "__main__":
    main()

