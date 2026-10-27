from __future__ import annotations

import csv
import hashlib
import json
from datetime import datetime
from pathlib import Path


ROOT = Path(__file__).resolve().parent


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


checks: list[dict[str, object]] = []


def check(name: str, passed: bool, detail: str) -> None:
    checks.append({"check": name, "passed": bool(passed), "detail": detail})


required = [
    "PHASE110_TARGETED_SUPPLEMENT_PLAN_v1.md",
    "CURRENT_ASSET_UTILIZATION_MATRIX_v1.tsv",
    "CHAIN_GAP_ACTION_MATRIX_v1.tsv",
    "PLANNED_MODULE_REGISTRY_v1.tsv",
    "GSE175540_24_SECTION_REUSE_PLAN_v1.md",
    "SOURCE_REGISTRY_v1.tsv",
]
for name in required:
    check(f"required file {name}", (ROOT / name).is_file(), str(ROOT / name))

tables = {
    "CURRENT_ASSET_UTILIZATION_MATRIX_v1.tsv": ("asset_id", 14),
    "CHAIN_GAP_ACTION_MATRIX_v1.tsv": ("edge_id", 9),
    "PLANNED_MODULE_REGISTRY_v1.tsv": ("module_id", 8),
    "SOURCE_REGISTRY_v1.tsv": ("source_id", 12),
}
for name, (id_col, expected_rows) in tables.items():
    path = ROOT / name
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    check(f"row count {name}", len(rows) == expected_rows, f"rows={len(rows)} expected={expected_rows}")
    ids = [row[id_col] for row in rows]
    check(f"unique IDs {name}", len(ids) == len(set(ids)), f"rows={len(ids)} unique={len(set(ids))}")
    check(f"complete rows {name}", all(all(value != "" for value in row.values()) for row in rows), "no blank fields")

plan = (ROOT / "PHASE110_TARGETED_SUPPLEMENT_PLAN_v1.md").read_text(encoding="utf-8")
reuse = (ROOT / "GSE175540_24_SECTION_REUSE_PLAN_v1.md").read_text(encoding="utf-8")
for token in ["72,952", "12 FFPE", "12 frozen", "M110A", "M110B", "M110C", "M110D", "section is the inferential unit", "No new biological analysis"]:
    check(f"plan token {token}", token.lower() in plan.lower(), token)
for token in ["24 processed", "23 TLS", "Do not rerun cell2location", "genuinely missing"]:
    check(f"reuse token {token}", token.lower() in reuse.lower(), token)
check("causal boundary retained", "cannot close metabolite or causal edges" in plan or "cannot close metabolite or causal edges" in reuse, "metabolite/causal boundary")
check("Phase109 reclassified not erased", "Phase109 is retained" in plan, "Phase109 retained as completed gap test")

status = "PASS_PHASE110_TARGETED_SUPPLEMENT_PLAN_VALIDATION" if all(item["passed"] for item in checks) else "FAIL_PHASE110_TARGETED_SUPPLEMENT_PLAN_VALIDATION"
payload = {
    "status": status,
    "generated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
    "checks_passed": sum(bool(item["passed"]) for item in checks),
    "checks_total": len(checks),
    "analysis_executed": False,
    "planned_first_block": ["M110A", "M110B", "M110C", "M110D"],
    "checks": checks,
}
(ROOT / "PHASE110_PLANNING_VALIDATION_v1.json").write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

manifest_rows = []
for path in sorted(ROOT.iterdir()):
    if path.is_file() and path.name != "FILE_MANIFEST_SHA256_v1.tsv":
        manifest_rows.append({"file": path.name, "bytes": path.stat().st_size, "sha256": sha256(path)})
with (ROOT / "FILE_MANIFEST_SHA256_v1.tsv").open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=["file", "bytes", "sha256"], delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(manifest_rows)

print(json.dumps({"status": status, "checks_passed": payload["checks_passed"], "checks_total": len(checks), "files_hashed": len(manifest_rows)}, ensure_ascii=False, indent=2))
if status.startswith("FAIL"):
    raise SystemExit(1)
