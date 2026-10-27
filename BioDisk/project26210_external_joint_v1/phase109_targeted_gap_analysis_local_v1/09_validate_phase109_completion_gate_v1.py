from __future__ import annotations

import csv
import hashlib
import json
import re
from datetime import datetime
from pathlib import Path


ROOT = Path(r"D:\project26210_external_joint_v1")
PHASE = ROOT / "phase109_targeted_gap_analysis_local_v1"
NOW = datetime.now().astimezone().isoformat(timespec="seconds")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


checks: list[dict[str, object]] = []


def check(name: str, passed: bool, detail: str) -> None:
    checks.append({"check": name, "passed": bool(passed), "detail": detail})


markers = {
    "PROJECT_STATUS.md": "## Phase 109 targeted supplementation",
    "DECISIONS.md": "## Phase 109 targeted-gap decisions",
    "CODEX_HANDOFF.md": "## Phase 109 current breakpoint",
}
for filename, marker in markers.items():
    text = (ROOT / filename).read_text(encoding="utf-8-sig")
    count = text.count(marker)
    check(f"unique marker {filename}", count == 1, f"count={count}")

analysis = json.loads((PHASE / "analysis_v3" / "PHASE109_ANALYSIS_STATUS_v1.json").read_text(encoding="utf-8"))
validation = json.loads((PHASE / "independent_validation_v1" / "PHASE109_INDEPENDENT_VALIDATION_STATUS_v1.json").read_text(encoding="utf-8"))
state = json.loads((PHASE / "PHASE109_STATE_UPDATE_STATUS_v1.json").read_text(encoding="utf-8"))
check("formal analysis PASS", analysis.get("status") == "PASS_ANALYSIS_COMPLETE", str(analysis.get("status")))
check("independent validation 54/54", validation.get("checks_passed") == 54 and validation.get("checks_total") == 54, f"{validation.get('checks_passed')}/{validation.get('checks_total')}")
check("state update PASS", state.get("status") == "PASS_PHASE109_STATE_UPDATE", str(state.get("status")))
check("primary abundance support count", analysis.get("primary_abundance_supported") == 0, str(analysis.get("primary_abundance_supported")))
check("fine-state CellChat support count", analysis.get("fine_state_cellchat_supported") == 2, str(analysis.get("fine_state_cellchat_supported")))
check("direct MTA remains untested", analysis.get("direct_mta_tested") is False, str(analysis.get("direct_mta_tested")))
check("causal MTA-SPP1 remains untested", analysis.get("causal_mta_to_spp1_tested") is False, str(analysis.get("causal_mta_to_spp1_tested")))

manifest_path = ROOT / "RESULTS_MANIFEST.tsv"
with manifest_path.open("r", encoding="utf-8-sig", newline="") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    fieldnames = reader.fieldnames
    rows = list(reader)
check("manifest IDs globally unique before final rows", len(rows) == len({r["result_id"] for r in rows}), f"rows={len(rows)} unique={len({r['result_id'] for r in rows})}")

validator_path = PHASE / "09_validate_phase109_completion_gate_v1.py"
output_path = PHASE / "PHASE109_COMPLETION_GATE_VALIDATION_v1.json"
pre_status = "PASS" if all(item["passed"] for item in checks) else "FAIL"
pre_payload = {
    "status": "PASS_PHASE109_COMPLETION_GATE" if pre_status == "PASS" else "FAIL_PHASE109_COMPLETION_GATE",
    "generated_at": NOW,
    "checks_passed": sum(bool(item["passed"]) for item in checks),
    "checks_total": len(checks),
    "checks": checks,
}
output_path.write_text(json.dumps(pre_payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
if pre_status != "PASS":
    raise RuntimeError(json.dumps(pre_payload, ensure_ascii=False, indent=2))

extra = [
    ("phase109_state_update_status_v1", "state_update_json", PHASE / "PHASE109_STATE_UPDATE_STATUS_v1.json"),
    ("phase109_completion_validator_v1", "formal_validation_script", validator_path),
    ("phase109_completion_gate_validation_v1", "state_update_validation_json", output_path),
]
existing = {r["result_id"] for r in rows}
if any(result_id in existing for result_id, _, _ in extra):
    raise RuntimeError("Completion rows already exist; refusing non-append-only rerun")
with manifest_path.open("a", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
    for result_id, artifact_type, path in extra:
        writer.writerow(
            {
                "result_id": result_id,
                "stage": "phase109_project_state_update",
                "artifact_type": artifact_type,
                "absolute_path": str(path),
                "generated_at": NOW,
                "validation_status": "PASS",
                "validation_evidence": f"SHA256 {sha256(path).upper()}; Phase109 mandatory four-file completion gate",
                "notes": "Phase109 append-only state synchronization evidence",
            }
        )

with manifest_path.open("r", encoding="utf-8-sig", newline="") as handle:
    final_rows = list(csv.DictReader(handle, delimiter="\t"))
phase_rows = [r for r in final_rows if r["stage"] in {"phase109_targeted_gap_analysis", "phase109_project_state_update"}]
final_checks: list[dict[str, object]] = []


def fcheck(name: str, passed: bool, detail: str) -> None:
    final_checks.append({"check": name, "passed": bool(passed), "detail": detail})


fcheck("manifest IDs globally unique", len(final_rows) == len({r["result_id"] for r in final_rows}), f"rows={len(final_rows)}")
fcheck("Phase109 manifest row count", len(phase_rows) == 24, f"rows={len(phase_rows)}")
for row in phase_rows:
    path = Path(row["absolute_path"])
    fcheck(f"path exists {row['result_id']}", path.is_file(), str(path))
    match = re.search(r"SHA256 ([0-9A-F]{64})", row["validation_evidence"])
    actual = sha256(path).upper() if path.is_file() else "MISSING"
    fcheck(f"hash matches {row['result_id']}", bool(match) and match.group(1) == actual, actual)

if not all(item["passed"] for item in final_checks):
    raise RuntimeError(json.dumps({"status": "FAIL_PHASE109_FINAL_MANIFEST_CHECK", "checks": final_checks}, ensure_ascii=False, indent=2))

receipt = {
    "status": "PASS_PHASE109_COMPLETION_GATE",
    "pre_checks": f"{pre_payload['checks_passed']}/{pre_payload['checks_total']}",
    "final_checks": f"{sum(bool(item['passed']) for item in final_checks)}/{len(final_checks)}",
    "phase109_manifest_rows": len(phase_rows),
    "four_file_sha256": {name: sha256(ROOT / name) for name in ["PROJECT_STATUS.md", "DECISIONS.md", "RESULTS_MANIFEST.tsv", "CODEX_HANDOFF.md"]},
}
print(json.dumps(receipt, ensure_ascii=False, indent=2))
