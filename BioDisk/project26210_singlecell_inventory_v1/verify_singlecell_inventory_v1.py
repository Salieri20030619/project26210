#!/usr/bin/env python3
"""Verify the staged single-cell inventory and emit reproducibility records."""

from __future__ import annotations

import csv
import hashlib
import json
from datetime import datetime
from pathlib import Path


ROOT = Path(__file__).resolve().parent


def read_tsv(name: str) -> list[dict[str, str]]:
    with (ROOT / name).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


checks: list[tuple[str, bool, str]] = []


def check(name: str, passed: bool, detail: str) -> None:
    checks.append((name, passed, detail))


inventory = json.loads((ROOT / "singlecell_h5ad_metadata_inventory_v1.json").read_text())
joint = json.loads((ROOT / "joint_reference_coverage_metadata_v1.json").read_text())
current = json.loads((ROOT / "current_joint_cell_pool_metadata_v1.json").read_text())

check("inventory_status", inventory.get("status") == "complete", str(inventory.get("status")))
check("joint_status", joint.get("status") == "complete", str(joint.get("status")))
check("current_status", current.get("status") == "complete", str(current.get("status")))
check("current_cells", current.get("n_cells") == 258120, str(current.get("n_cells")))
check("current_genes", current.get("n_genes") == 37795, str(current.get("n_genes")))
check("reference_cells", joint.get("n_cells") == 230618, str(joint.get("n_cells")))
check("reference_patients", joint.get("n_patients") == 40, str(joint.get("n_patients")))
check("reference_samples", joint.get("n_samples") == 75, str(joint.get("n_samples")))
check("reference_studies", joint.get("n_studies") == 5, str(joint.get("n_studies")))

pool = read_tsv("current_joint_cell_pool_coverage_v1.tsv")
pool_by_name = {row["joint_group_v1"]: row for row in pool}
for cell_pool, expected in {
    "Macrophage": 42546,
    "T cell": 63071,
    "Endothelial": 30025,
    "Conservative malignant": 11950,
}.items():
    observed = int(pool_by_name[cell_pool]["n_cells"])
    check(f"pool_{cell_pool}", observed == expected, f"observed={observed}; expected={expected}")

c2l = read_tsv("c2l_celltype_patient_sample_coverage_v1.tsv")
c2l_by_name = {row["c2l_celltype"]: row for row in c2l}
for cell_type, expected in {
    "SPP1+ TAM": 11290,
    "APOE+C1QC+ TAM": 6030,
    "TYMP-High Malignant": 1720,
    "TYMP-Low Malignant": 4730,
}.items():
    observed = int(c2l_by_name[cell_type]["n_cells"])
    check(f"c2l_{cell_type}", observed == expected, f"observed={observed}; expected={expected}")

candidates = read_tsv("candidate_dataset_priority_v1.tsv")
check("candidate_rows", len(candidates) == 8, f"observed={len(candidates)}; expected=8")
check(
    "candidate_p0_metadata",
    any(row["accession"] == "GSE294109 + GSE294006" and row["priority"] == "P0_metadata_first" for row in candidates),
    "GSE294109 + GSE294006 must be P0_metadata_first",
)
check(
    "candidate_p0_clinical",
    any(row["accession"] == "GSE264586" and row["priority"] == "P0_if_clinical_story" for row in candidates),
    "GSE264586 must be P0_if_clinical_story",
)

report = (ROOT / "SINGLECELL_INVENTORY_AND_DATASET_DECISION_v1.md").read_text(encoding="utf-8")
for phrase in (
    "现有本地数据足以启动",
    "旧版细分标签不能直接作为新文章的最终参考",
    "空间中的 SPP1 信号不能自动归因于 TAM",
    "没有下载任何外部数据",
):
    check(f"report_phrase_{phrase[:12]}", phrase in report, phrase)

check_path = ROOT / "inventory_validation_checks_v1.tsv"
with check_path.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(["check", "status", "detail"])
    for name, passed, detail in checks:
        writer.writerow([name, "PASS" if passed else "FAIL", detail])

failed = [name for name, passed, _ in checks if not passed]

manifest_path = ROOT / "inventory_file_manifest_v1.tsv"
excluded = {manifest_path.name, "inventory_verification_summary_v1.json"}
files = sorted(path for path in ROOT.iterdir() if path.is_file() and path.name not in excluded)
with manifest_path.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(["file", "size_bytes", "sha256"])
    for path in files:
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        writer.writerow([path.name, path.stat().st_size, digest])

summary = {
    "status": "complete" if not failed else "failed",
    "verified_at": datetime.now().astimezone().isoformat(timespec="seconds"),
    "n_checks": len(checks),
    "n_passed": len(checks) - len(failed),
    "n_failed": len(failed),
    "failed_checks": failed,
    "checks": str(check_path),
    "manifest": str(manifest_path),
}
(ROOT / "inventory_verification_summary_v1.json").write_text(
    json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
)
print(json.dumps(summary, ensure_ascii=False, indent=2))
raise SystemExit(1 if failed else 0)
