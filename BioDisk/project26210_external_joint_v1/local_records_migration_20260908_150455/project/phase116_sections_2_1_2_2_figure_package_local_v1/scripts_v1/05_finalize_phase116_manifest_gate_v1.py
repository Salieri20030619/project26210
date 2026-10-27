from __future__ import annotations

import csv
import hashlib
import json
import re
import shutil
from pathlib import Path


ROOT = Path(r"D:\project26210_external_joint_v1")
PHASE = ROOT / "phase116_sections_2_1_2_2_figure_package_local_v1"
SCRIPTS = PHASE / "scripts_v1"
VALID = PHASE / "validation_v1"
BACKUP = PHASE / "state_backup_before_phase116_manifest_finalize_v1"
MANIFEST = ROOT / "RESULTS_MANIFEST.tsv"
STATE = ["PROJECT_STATUS.md", "DECISIONS.md", "RESULTS_MANIFEST.tsv", "CODEX_HANDOFF.md"]
HEADER = ["result_id", "stage", "artifact_type", "absolute_path", "generated_at", "validation_status", "validation_evidence", "notes"]
STAMP = "2026-09-02T17:15:00+08:00"


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest().upper()


def row_for(result_id, artifact_type, path, note):
    return {
        "result_id": result_id,
        "stage": "phase116_project_state_update",
        "artifact_type": artifact_type,
        "absolute_path": str(path),
        "generated_at": STAMP,
        "validation_status": "PASS",
        "validation_evidence": f"Phase116 final state consistency gate; SHA256 {sha(path)}",
        "notes": note,
    }


def main():
    if BACKUP.exists():
        raise FileExistsError(BACKUP)
    BACKUP.mkdir(parents=True)
    for name in STATE:
        shutil.copy2(ROOT / name, BACKUP / name)

    with MANIFEST.open("r", encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        rows = list(reader)
        if reader.fieldnames != HEADER:
            raise RuntimeError(reader.fieldnames)

    for row in rows:
        if not row["result_id"].strip('"').startswith("phase116_"):
            continue
        path = Path(row["absolute_path"])
        if path.is_file() and "SHA256 " not in row["validation_evidence"]:
            row["validation_evidence"] = row["validation_evidence"].rstrip() + f"; SHA256 {sha(path)}"

    add = [
        row_for("phase116_state_completion_script_v1", "formal_state_script", SCRIPTS / "04_complete_phase116_state_v1.py", "Creates backup and updates all four project-state files"),
        row_for("phase116_manifest_finalize_script_v1", "formal_state_validation_script", SCRIPTS / "05_finalize_phase116_manifest_gate_v1.py", "Adds per-artifact hashes and final state gate"),
    ]
    ids = [r["result_id"].strip('"') for r in rows]
    for item in add:
        if item["result_id"] in ids:
            raise RuntimeError(f"duplicate id {item['result_id']}")
        rows.append(item)
        ids.append(item["result_id"])

    with MANIFEST.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=HEADER, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)

    checks = {}
    checks["manifest_schema"] = True
    checks["manifest_unique_ids"] = len(ids) == len(set(ids))
    phase_rows = [r for r in rows if r["result_id"].strip('"').startswith("phase116_")]
    checks["phase116_paths_exist"] = all(Path(r["absolute_path"]).is_file() for r in phase_rows)
    hash_ok = True
    for row in phase_rows:
        match = re.search(r"SHA256 ([0-9A-Fa-f]{64})", row["validation_evidence"])
        if not match or sha(Path(row["absolute_path"])) != match.group(1).upper():
            hash_ok = False
    checks["phase116_registered_hashes_match"] = hash_ok
    checks["project_status_complete"] = "PASS_PACKAGE_COMPLETE" in (ROOT / "PROJECT_STATUS.md").read_text(encoding="utf-8")
    checks["decision_boundary_recorded"] = "Reject the wording that TYMP has been proved" in (ROOT / "DECISIONS.md").read_text(encoding="utf-8")
    checks["handoff_resume_point_recorded"] = "## Phase 116 2.1–2.2 figure package current breakpoint" in (ROOT / "CODEX_HANDOFF.md").read_text(encoding="utf-8")
    checks["pre_finalize_backup_complete"] = all((BACKUP / x).is_file() for x in STATE)
    checks["package_validation_pass"] = json.loads((VALID / "PACKAGE_VALIDATION_v1.json").read_text(encoding="utf-8"))["status"] == "PASS"
    gate = {
        "status": "PASS" if all(checks.values()) else "FAIL",
        "n_checks": len(checks),
        "n_failed": sum(not x for x in checks.values()),
        "checks": checks,
        "phase116_manifest_rows_before_gate_self_registration": len(phase_rows),
        "state_hashes_before_gate_self_registration": {name: sha(ROOT / name) for name in STATE},
        "remote_mirror": "/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1",
        "remote_sync_required_before_completion_report": True,
    }
    gate_path = VALID / "PHASE116_FINAL_STATE_VALIDATION_v1.json"
    gate_path.write_text(json.dumps(gate, ensure_ascii=False, indent=2), encoding="utf-8")
    if gate["status"] != "PASS":
        raise RuntimeError(json.dumps(gate, ensure_ascii=False))

    self_row = row_for("phase116_final_state_validation_v1", "state_update_validation_json", gate_path, "Final manifest schema, unique-ID, path, hash, marker and backup validation")
    if self_row["result_id"] in ids:
        raise RuntimeError("duplicate final validation id")
    with MANIFEST.open("a", encoding="utf-8", newline="") as f:
        csv.DictWriter(f, fieldnames=HEADER, delimiter="\t", lineterminator="\n").writerow(self_row)

    with MANIFEST.open("r", encoding="utf-8-sig", newline="") as f:
        final_rows = list(csv.DictReader(f, delimiter="\t"))
    final_ids = [r["result_id"].strip('"') for r in final_rows]
    if len(final_ids) != len(set(final_ids)):
        raise RuntimeError("post-registration duplicate IDs")
    print(json.dumps({"status": "PASS", "checks": len(checks), "phase116_rows": sum(x.startswith('phase116_') for x in final_ids), "gate": str(gate_path)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
