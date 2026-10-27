from __future__ import annotations

import csv
import hashlib
import json
import shutil
from datetime import datetime, timezone, timedelta
from pathlib import Path


ROOT = Path(r"D:\project26210_external_joint_v1")
PHASE = ROOT / "phase117_23figure_word_package_correction_local_v1"
MANIFEST = ROOT / "RESULTS_MANIFEST.tsv"
SCRIPT_DST = PHASE / "scripts_v1" / "05_finalize_phase117_state_gate_v2.py"
VALIDATION = PHASE / "validation_v2" / "PHASE117_FINAL_STATE_VALIDATION_v2.json"
BACKUP = PHASE / "state_backup_before_phase117_v1"
DELIV = PHASE / "deliverable_v2"
NOW = datetime.now(timezone(timedelta(hours=8))).isoformat(timespec="seconds")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest().upper()


shutil.copy2(Path(__file__), SCRIPT_DST)
rows = list(csv.DictReader(MANIFEST.open("r", encoding="utf-8-sig", newline=""), delimiter="\t"))
if not any(r["result_id"] == "phase117_state_finalize_script_v2" for r in rows):
    with MANIFEST.open("a", encoding="utf-8", newline="") as fh:
        csv.writer(fh, delimiter="\t", lineterminator="\n").writerow([
            "phase117_state_finalize_script_v2",
            "phase117_23figure_word_package_correction",
            "formal_state_validation_script",
            str(SCRIPT_DST),
            NOW,
            "PASS",
            f"Phase117 corrected final state gate; SHA256 {sha256(SCRIPT_DST)}",
            "Final manifest schema, unique-ID, path, marker, backup and package validation",
        ])

rows = list(csv.DictReader(MANIFEST.open("r", encoding="utf-8-sig", newline=""), delimiter="\t"))
ids = [r["result_id"] for r in rows]
phase_rows = [r for r in rows if r["stage"] == "phase117_23figure_word_package_correction"]
package_validation = json.loads((PHASE / "validation_v2" / "PACKAGE_VALIDATION_v2.json").read_text(encoding="utf-8"))
delivery_validation = json.loads((PHASE / "INSTALL_AND_DELIVERY_VALIDATION_v2.json").read_text(encoding="utf-8"))
checks = {
    "manifest_schema": list(rows[0].keys()) == ["result_id", "stage", "artifact_type", "absolute_path", "generated_at", "validation_status", "validation_evidence", "notes"],
    "manifest_unique_ids": len(ids) == len(set(ids)),
    "phase117_row_count_18": len(phase_rows) == 18,
    "phase117_paths_exist": all(Path(r["absolute_path"]).exists() for r in phase_rows),
    "phase117_status_pass": all(r["validation_status"] == "PASS" for r in phase_rows),
    "state_backup_four_files": all((BACKUP / n).exists() for n in ["PROJECT_STATUS.md", "DECISIONS.md", "RESULTS_MANIFEST.tsv", "CODEX_HANDOFF.md"]),
    "project_status_marker": "Phase 117" in (ROOT / "PROJECT_STATUS.md").read_text(encoding="utf-8"),
    "decisions_marker": "Phase 117" in (ROOT / "DECISIONS.md").read_text(encoding="utf-8"),
    "handoff_marker": "Phase 117" in (ROOT / "CODEX_HANDOFF.md").read_text(encoding="utf-8"),
    "package_validation_18_of_18": package_validation["status"] == "PASS" and package_validation["n_checks"] == 18 and package_validation["n_failed"] == 0,
    "word_has_23_inline_images": package_validation["word_inline_images"] == 23,
    "word_rendered_pages_28": package_validation["rendered_pages"] == 28,
    "no_tsv_json_in_user_package": not any(p.suffix.lower() in {".tsv", ".json"} for p in DELIV.rglob("*") if p.is_file()),
    "downloads_hash_match_44": delivery_validation["downloads_exact_hash_matches"] == 44,
    "zip_entry_count_44": delivery_validation["zip_file_entries"] == 44,
}
result = {
    "status": "PASS" if all(checks.values()) else "FAIL",
    "generated_at": NOW,
    "n_checks": len(checks),
    "n_passed": sum(checks.values()),
    "checks": checks,
}
VALIDATION.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
if result["status"] != "PASS":
    raise RuntimeError(json.dumps(result, ensure_ascii=False, indent=2))
print(json.dumps(result, ensure_ascii=False, indent=2))
