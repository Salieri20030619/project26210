#!/usr/bin/env python3
from __future__ import annotations

import csv
import hashlib
import json
from datetime import datetime
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parent
ANALYSIS = ROOT / "analysis_local_v1"
VALIDATION = ROOT / "validation_local_v1"
INTEGRATION = ROOT / "integration_local_v1"
SERVER = ROOT / "server_results_sync_v1"
OUT = ROOT / "terminal_validation_local_v1"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    if OUT.exists(): raise FileExistsError(OUT)
    OUT.mkdir()
    v = json.loads((VALIDATION / "status.json").read_text())
    a = json.loads((ANALYSIS / "status.json").read_text())
    r = json.loads((SERVER / "status.json").read_text())
    i = json.loads((INTEGRATION / "status.json").read_text())
    conclusions = pd.read_csv(INTEGRATION / "y02_cross_method_endpoint_conclusions_v1.tsv", sep="\t").set_index("endpoint")
    checks = [
        ("remote_compute", r["status"] == "PASS_Y02_MONOCLE3_COMPUTE", r["status"]),
        ("analysis", a["status"] == "PASS_Y02_MONOCLE3_ANALYSIS_COMPLETE_PRESENTATION_PENDING", a["status"]),
        ("independent_validation", v["status"] == "PASS_INDEPENDENT_VALIDATION" and v["checks_passed"] == v["checks_total"] == 13, f"{v['checks_passed']}/{v['checks_total']}"),
        ("integration", i["status"] == "PASS_Y02_CROSS_METHOD_INTEGRATION", i["status"]),
        ("scope", r["n_cells"] == 50761 and r["n_studies"] == 5 and r["n_patients"] == 40, "50761 cells; 5 studies; 40 patients"),
        ("five_rds", len(list(SERVER.glob("*_monocle3_cds_v1.rds"))) == 5, str(len(list(SERVER.glob("*_monocle3_cds_v1.rds"))))),
        ("state_endpoint_positive", conclusions.loc["SPP1_minus_Monocyte_pseudotime", "cross_method_conclusion"] == "CROSS_METHOD_POSITIVE_DIRECTION", conclusions.loc["SPP1_minus_Monocyte_pseudotime", "cross_method_conclusion"]),
        ("spp1_program_positive", conclusions.loc["spp1_state_no_spp1", "cross_method_conclusion"] == "CROSS_METHOD_POSITIVE_DIRECTION", conclusions.loc["spp1_state_no_spp1", "cross_method_conclusion"]),
        ("spp1_vs_c1qc_negative", conclusions.loc["spp1_vs_c1qc", "cross_method_conclusion"] == "CROSS_METHOD_NEGATIVE_DIRECTION", conclusions.loc["spp1_vs_c1qc", "cross_method_conclusion"]),
        ("two_sided", a["two_sided_testing"] is True and a["direction_prespecified"] is False, "two-sided; no direction prespecified"),
        ("interpretation_guard", "not lineage tracing" in a["interpretation"], a["interpretation"]),
        ("no_figures", a["figures_created"] is False and r["figures_created"] is False and not any(ROOT.rglob("*.png")) and not any(ROOT.rglob("*.pdf")), "no figures"),
    ]
    failed = [name for name, ok, _ in checks if not ok]
    with (OUT / "terminal_checks_v1.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n"); writer.writerow(["check", "pass", "detail"]); writer.writerows((name, str(ok), detail) for name, ok, detail in checks)
    files = []
    for directory in (SERVER, ANALYSIS, VALIDATION, INTEGRATION):
        for path in sorted(directory.iterdir()):
            if path.is_file(): files.append({"relative_path":str(path.relative_to(ROOT)), "bytes":path.stat().st_size, "sha256":sha256(path)})
    pd.DataFrame(files).to_csv(OUT / "validated_file_manifest_v1.tsv", sep="\t", index=False, lineterminator="\n")
    payload = {"status":"PASS_PHASE33_Y02_TERMINAL_VALIDATION" if not failed else "FAIL_PHASE33_Y02_TERMINAL_VALIDATION", "timestamp_local":datetime.now().astimezone().isoformat(timespec="seconds"), "checks_passed":len(checks)-len(failed), "checks_total":len(checks), "component_checks":"13/13", "failed_checks":failed, "presentation_status":"FIGURES_PENDING_USER_REVIEW"}
    (OUT / "status.json").write_text(json.dumps(payload, indent=2)+"\n", encoding="utf-8")
    if failed: raise RuntimeError(f"Terminal validation failed: {failed}")


if __name__ == "__main__": main()
