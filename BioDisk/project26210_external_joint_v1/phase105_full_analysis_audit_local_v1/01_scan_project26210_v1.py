from __future__ import annotations

import csv
import json
import re
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path


ROOT = Path(r"D:\project26210_external_joint_v1")
OUT = Path(r"C:\Users\31937\project26210_phase105_full_analysis_audit_v1")
OUT.mkdir(parents=True, exist_ok=True)


def write_tsv(path: Path, rows: list[dict[str, object]], columns: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def phase_number(text: str) -> int | None:
    match = re.search(r"(?i)phase[_ -]?(\d+)", text)
    return int(match.group(1)) if match else None


with (ROOT / "RESULTS_MANIFEST.tsv").open("r", encoding="utf-8-sig", newline="") as handle:
    manifest = list(csv.DictReader(handle, delimiter="\t"))

manifest_by_phase: dict[int, list[dict[str, str]]] = defaultdict(list)
for row in manifest:
    number = phase_number(row["stage"])
    if number is not None:
        manifest_by_phase[number].append(row)

status_text = (ROOT / "PROJECT_STATUS.md").read_text(encoding="utf-8-sig")
decision_text = (ROOT / "DECISIONS.md").read_text(encoding="utf-8-sig")
status_lines = status_text.splitlines()

heading_rows: list[dict[str, object]] = []
status_heading_counts: Counter[int] = Counter()
for index, line in enumerate(status_lines, start=1):
    if line.startswith("## "):
        number = phase_number(line)
        if number is not None:
            status_heading_counts[number] += 1
        heading_rows.append({"line_number": index, "phase": number or "", "heading": line[3:]})
write_tsv(OUT / "project_status_headings_v1.tsv", heading_rows, ["line_number", "phase", "heading"])

decision_heading_counts: Counter[int] = Counter()
for line in decision_text.splitlines():
    if line.startswith("## "):
        number = phase_number(line)
        if number is not None:
            decision_heading_counts[number] += 1

directories_by_phase: dict[int, list[Path]] = defaultdict(list)
directory_rows: list[dict[str, object]] = []
for directory in sorted((item for item in ROOT.iterdir() if item.is_dir()), key=lambda p: p.name.lower()):
    number = phase_number(directory.name)
    if number is None:
        continue
    directories_by_phase[number].append(directory)
    files = [item for item in directory.rglob("*") if item.is_file()]
    latest = max((item.stat().st_mtime for item in files), default=directory.stat().st_mtime)
    summary_like = [
        item
        for item in files
        if re.search(r"(?i)(summary|status|decision|addendum|validation|register|audit)", item.name)
    ]
    registered = [
        row
        for row in manifest
        if row["absolute_path"].lower().startswith(str(directory).lower() + "\\")
        or row["absolute_path"].lower() == str(directory).lower()
    ]
    directory_rows.append(
        {
            "phase": number,
            "directory": directory.name,
            "absolute_path": str(directory),
            "file_count": len(files),
            "total_size_bytes": sum(item.stat().st_size for item in files),
            "latest_file_time": datetime.fromtimestamp(latest).astimezone().isoformat(timespec="seconds"),
            "summary_like_file_count": len(summary_like),
            "summary_like_examples": ";".join(str(item.relative_to(directory)) for item in summary_like[:12]),
            "manifest_rows_under_directory": len(registered),
            "directory_registration": "REGISTERED" if registered else "NO_MANIFEST_PATH",
        }
    )
write_tsv(
    OUT / "top_level_phase_directory_inventory_v1.tsv",
    directory_rows,
    [
        "phase",
        "directory",
        "absolute_path",
        "file_count",
        "total_size_bytes",
        "latest_file_time",
        "summary_like_file_count",
        "summary_like_examples",
        "manifest_rows_under_directory",
        "directory_registration",
    ],
)

path_audit_rows: list[dict[str, object]] = []
for row in manifest:
    value = row["absolute_path"]
    if re.match(r"^[A-Za-z]:\\", value):
        location = "LOCAL_WINDOWS"
        exists = Path(value).is_file()
        check = "EXISTS" if exists else "MISSING_LOCAL_PATH"
    elif value.startswith("/mnt/"):
        location = "REMOTE_LINUX"
        exists = ""
        check = "REMOTE_CHECK_REQUIRED"
    else:
        location = "OTHER"
        exists = Path(value).exists()
        check = "EXISTS" if exists else "UNRESOLVED_PATH"
    path_audit_rows.append(
        {
            "result_id": row["result_id"],
            "stage": row["stage"],
            "phase": phase_number(row["stage"]) or "",
            "artifact_type": row["artifact_type"],
            "location": location,
            "path_check": check,
            "validation_status": row["validation_status"],
            "absolute_path": value,
        }
    )
write_tsv(
    OUT / "manifest_path_audit_v1.tsv",
    path_audit_rows,
    [
        "result_id",
        "stage",
        "phase",
        "artifact_type",
        "location",
        "path_check",
        "validation_status",
        "absolute_path",
    ],
)

presence_rows: list[dict[str, object]] = []
for number in range(1, 105):
    dirs = directories_by_phase.get(number, [])
    manifest_rows = manifest_by_phase.get(number, [])
    file_count = sum(row["file_count"] for row in directory_rows if row["phase"] == number)
    flags: list[str] = []
    if dirs and not manifest_rows:
        flags.append("DIRECTORY_WITHOUT_MANIFEST_STAGE")
    if dirs and status_heading_counts[number] == 0:
        flags.append("DIRECTORY_WITHOUT_STATUS_HEADING")
    if manifest_rows and not dirs:
        flags.append("MANIFEST_WITHOUT_TOP_LEVEL_LOCAL_DIRECTORY")
    if any(row["validation_status"].upper().startswith(("FAIL", "BLOCK")) for row in manifest_rows):
        flags.append("HAS_FAILURE_OR_BLOCKED_RECORD")
    presence_rows.append(
        {
            "phase": number,
            "top_level_directory_count": len(dirs),
            "top_level_directories": ";".join(item.name for item in dirs),
            "file_count": file_count,
            "manifest_row_count": len(manifest_rows),
            "manifest_statuses": ";".join(sorted({row["validation_status"] for row in manifest_rows})),
            "project_status_heading_count": status_heading_counts[number],
            "decision_heading_count": decision_heading_counts[number],
            "audit_flags": ";".join(flags) if flags else "NONE",
        }
    )
write_tsv(
    OUT / "phase_presence_matrix_v1.tsv",
    presence_rows,
    [
        "phase",
        "top_level_directory_count",
        "top_level_directories",
        "file_count",
        "manifest_row_count",
        "manifest_statuses",
        "project_status_heading_count",
        "decision_heading_count",
        "audit_flags",
    ],
)

stage_rows: list[dict[str, object]] = []
for stage in sorted({row["stage"] for row in manifest}):
    group = [row for row in manifest if row["stage"] == stage]
    stage_rows.append(
        {
            "stage": stage,
            "phase": phase_number(stage) or "",
            "rows": len(group),
            "artifact_types": ";".join(sorted({row["artifact_type"] for row in group})),
            "validation_statuses": ";".join(sorted({row["validation_status"] for row in group})),
            "local_paths": sum(re.match(r"^[A-Za-z]:\\", row["absolute_path"]) is not None for row in group),
            "remote_paths": sum(row["absolute_path"].startswith("/mnt/") for row in group),
            "missing_local_paths": sum(
                re.match(r"^[A-Za-z]:\\", row["absolute_path"]) is not None
                and not Path(row["absolute_path"]).is_file()
                for row in group
            ),
        }
    )
write_tsv(
    OUT / "manifest_stage_summary_v1.tsv",
    stage_rows,
    [
        "stage",
        "phase",
        "rows",
        "artifact_types",
        "validation_statuses",
        "local_paths",
        "remote_paths",
        "missing_local_paths",
    ],
)

summary = {
    "status": "PASS_SCAN_COMPLETE",
    "generated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
    "authoritative_root": str(ROOT),
    "manifest_rows": len(manifest),
    "manifest_unique_ids": len({row["result_id"] for row in manifest}),
    "manifest_stages": len({row["stage"] for row in manifest}),
    "top_level_phase_directories": len(directory_rows),
    "phases_with_top_level_directories": sorted(directories_by_phase),
    "phases_with_manifest_rows": sorted(manifest_by_phase),
    "directory_without_manifest_count": sum(row["directory_registration"] == "NO_MANIFEST_PATH" for row in directory_rows),
    "missing_local_manifest_paths": sum(row["path_check"] == "MISSING_LOCAL_PATH" for row in path_audit_rows),
    "remote_manifest_paths_requiring_server_check": sum(row["path_check"] == "REMOTE_CHECK_REQUIRED" for row in path_audit_rows),
    "presence_flags": Counter(
        flag
        for row in presence_rows
        for flag in str(row["audit_flags"]).split(";")
        if flag != "NONE"
    ),
}
(OUT / "scan_summary_v1.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
print(json.dumps(summary, ensure_ascii=False, indent=2))
