#!/usr/bin/env python3
from __future__ import annotations

import csv
import hashlib
import json
from datetime import datetime
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parent
CORE = ROOT / "server_results_sync_v1"
RUN_METADATA = ROOT / "server_run_metadata_sync_v1"
ANALYSIS = ROOT / "analysis_local_v1"
VALIDATION = ROOT / "independent_validation_local_v2"
OUT = ROOT / "terminal_validation_local_v1"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    if OUT.exists():
        raise FileExistsError(OUT)
    OUT.mkdir()
    analysis_status = json.loads((ANALYSIS / "status.json").read_text(encoding="utf-8"))
    validation_status = json.loads((VALIDATION / "status.json").read_text(encoding="utf-8"))
    input_status = json.loads((RUN_METADATA / "input_status_v1.json").read_text(encoding="utf-8"))
    loom_status = json.loads((RUN_METADATA / "loom_status_v1.json").read_text(encoding="utf-8"))
    summary = pd.read_csv(ANALYSIS / "y03_candidate_endpoint_summary_v1.tsv", sep="\t")
    integrated = pd.read_csv(ANALYSIS / "y03_candidate_integrated_evidence_v1.tsv", sep="\t").set_index("tf")
    core_log = (RUN_METADATA / "logs" / "pyscenic_core_v1.log").read_text(encoding="utf-8", errors="replace")
    checks = [
        ("input", input_status["status"] == "PASS_Y03_PYSCENIC_INPUT_PREPARATION", input_status["status"]),
        ("loom", loom_status["status"] == "PASS_Y03_PYSCENIC_LOOM_BUILD", loom_status["status"]),
        ("core", "Y03_PYSCENIC_CORE_PASS" in core_log and (CORE / "RUN_COMPLETE").is_file(), "terminal PASS log plus completion marker"),
        ("analysis", analysis_status["status"] == "PASS_Y03_CANDIDATE_ANALYSIS_COMPLETE_PRESENTATION_PENDING", analysis_status["status"]),
        ("independent_validation", validation_status["status"] == "PASS_Y03_INDEPENDENT_VALIDATION" and validation_status["n_pass"] == validation_status["n_checks"] == 29, f"{validation_status['n_pass']}/{validation_status['n_checks']}"),
        ("scope", input_status["n_cells"] == 50761 and input_status["n_patients"] == 40 and input_status["n_studies"] == 5, "50761 cells; 40 patients; 5 studies"),
        ("candidate_endpoint_grid", len(summary) == 12 and set(summary["tf"]) == {"JUN", "JUNB", "CEBPD"}, f"{len(summary)} rows; 3 candidates"),
        ("integrated_result", integrated.loc["CEBPD", "n_supported_endpoints_of_4"] == 4 and integrated.loc["JUNB", "n_supported_endpoints_of_4"] == 4 and integrated.loc["JUN", "n_supported_endpoints_of_4"] == 1, "CEBPD 4/4; JUNB 4/4; JUN 1/4"),
        ("two_sided", analysis_status["two_sided_testing"] is True and analysis_status["direction_prespecified"] is False, "two-sided; no prespecified direction"),
        ("no_figures", analysis_status["figures_created"] is False and not any(ROOT.rglob("*.png")) and not any(ROOT.rglob("*.pdf")) and not any(ROOT.rglob("*.svg")), "no figures"),
    ]
    failures = [name for name, ok, _ in checks if not ok]
    with (OUT / "terminal_checks_v1.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["check", "pass", "detail"])
        writer.writerows((name, str(ok), detail) for name, ok, detail in checks)
    files = []
    for directory in (CORE, RUN_METADATA, ANALYSIS, VALIDATION):
        for path in sorted(directory.rglob("*")):
            if path.is_file():
                files.append({"relative_path": str(path.relative_to(ROOT)), "bytes": path.stat().st_size, "sha256": sha256(path)})
    pd.DataFrame(files).to_csv(OUT / "validated_file_manifest_v1.tsv", sep="\t", index=False, lineterminator="\n")
    payload = {
        "status": "PASS_PHASE34_Y03_TERMINAL_VALIDATION" if not failures else "FAIL_PHASE34_Y03_TERMINAL_VALIDATION",
        "timestamp_local": datetime.now().astimezone().isoformat(timespec="seconds"),
        "checks_passed": len(checks) - len(failures),
        "checks_total": len(checks),
        "component_checks": "29/29",
        "failed_checks": failures,
        "presentation_status": "HEATMAP_PENDING_USER_REVIEW",
    }
    (OUT / "status.json").write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    if failures:
        raise RuntimeError(f"Terminal validation failed: {failures}")


if __name__ == "__main__":
    main()
