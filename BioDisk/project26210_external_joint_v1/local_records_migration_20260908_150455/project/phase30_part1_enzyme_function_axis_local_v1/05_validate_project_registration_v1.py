#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import os
import re
from pathlib import Path

import pandas as pd


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    root = Path(__file__).resolve().parent
    project = root.parent
    output = root / "registration_validation_local_v2"
    part = root / "registration_validation_local_v2.part"
    if output.exists() or part.exists():
        raise FileExistsError("refusing to overwrite registration validation")
    part.mkdir(parents=True)
    checks = []

    status_text = (project / "PROJECT_STATUS.md").read_text(encoding="utf-8")
    checks.append({"check": "project_status", "status": "PASS" if all(x in status_text for x in
                   ["Phase 30 Part 1 TYMP/RRM2B enzyme-function-axis addendum",
                    "PASS_PHASE30_TERMINAL_VALIDATION", "SECONDARY_SUPPORT", "(3/4)"]) else "FAIL",
                   "detail": "Phase 30 status block and key conclusion"})

    manifest = pd.read_csv(project / "RESULTS_MANIFEST.tsv", sep="\t", dtype=str)
    ids = [
        "phase30_enzyme_axis_reaction_registry_local_v1",
        "phase30_enzyme_axis_endpoint_summary_local_v1",
        "phase30_enzyme_axis_all_axis_integration_local_v2",
        "phase30_enzyme_axis_independent_validation_local_v1",
        "phase30_enzyme_axis_terminal_validation_local_v1",
    ]
    subset = manifest[manifest["result_id"].isin(ids)]
    checks.append({"check": "manifest_rows", "status": "PASS" if len(subset) == 5 and subset["result_id"].nunique() == 5 else "FAIL",
                   "detail": f"rows={len(subset)}; unique={subset['result_id'].nunique()}"})
    paths_ok = True
    hashes_ok = True
    for row in subset.itertuples(index=False):
        path = Path(row.absolute_path)
        paths_ok &= path.is_file()
        match = re.search(r"SHA256 ([0-9a-fA-F]{64})", row.validation_evidence)
        hashes_ok &= bool(match) and path.is_file() and sha256(path).lower() == match.group(1).lower()
    checks.append({"check": "manifest_paths", "status": "PASS" if paths_ok else "FAIL", "detail": "all five paths exist"})
    checks.append({"check": "manifest_hashes", "status": "PASS" if hashes_ok else "FAIL", "detail": "all five SHA256 values match"})

    register = pd.read_csv(root / "PHASE30_COMPLETION_REGISTER_v2.tsv", sep="\t")
    register_ok = len(register) == 1 and register.iloc[0]["phase30_status"] == "COMPLETE_VALIDATED"
    register_ok &= register.iloc[0]["validation"] == "17/17 independent and 9/9 terminal checks PASS"
    checks.append({"check": "completion_register_v2", "status": "PASS" if register_ok else "FAIL",
                   "detail": register.iloc[0]["validation"]})

    frame = pd.DataFrame(checks)
    frame.to_csv(part / "registration_checks_v1.tsv", sep="\t", index=False)
    passed = int((frame["status"] == "PASS").sum())
    status = {"status": "PASS_PROJECT_REGISTRATION_VALIDATION" if passed == len(frame) else "FAIL",
              "checks_passed": passed, "checks_total": len(frame)}
    (part / "status.json").write_text(json.dumps(status, indent=2) + "\n", encoding="utf-8")
    os.replace(part, output)
    if passed != len(frame):
        raise RuntimeError(f"registration validation failed: {passed}/{len(frame)}")


if __name__ == "__main__":
    main()
