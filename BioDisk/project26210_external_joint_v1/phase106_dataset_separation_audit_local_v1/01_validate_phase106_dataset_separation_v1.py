from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path


ROOT = Path(r"C:\Users\31937\project26210_phase106_dataset_separation_audit_v1")


def tsv(name: str) -> list[dict[str, str]]:
    with (ROOT / name).open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


checks: list[dict[str, object]] = []


def check(name: str, passed: bool, detail: object = "") -> None:
    checks.append({"check": name, "passed": bool(passed), "detail": str(detail)})


required = [
    "PHASE106_DATASET_SEPARATION_RULE_v1.md",
    "DATASET_UNIVERSE_REGISTRY_v1.tsv",
    "PHASE_TO_DATASET_ROUTE_v1.tsv",
    "DATASET_SPECIFIC_CLAIM_MATRIX_v1.tsv",
    "NONCOMPARABILITY_AND_CONFLICT_REGISTER_v1.tsv",
    "CLIENT_RESPONSE_BLUEPRINT_v2.md",
    "PHASE106_DATASET_SEPARATION_AUDIT_SUMMARY_v1.md",
]
for name in required:
    path = ROOT / name
    check(f"required_{name}", path.is_file() and path.stat().st_size > 0, path)

registry = tsv("DATASET_UNIVERSE_REGISTRY_v1.tsv")
registry_ids = [row["universe_id"] for row in registry]
check("registry_rows_16", len(registry) == 16, len(registry))
check("registry_ids_unique", len(registry_ids) == len(set(registry_ids)))
check("registry_original_part1_exact", any(row["universe_id"] == "D01" and row["datasets_or_samples"] == "GSE156632;GSE159115;GSE178481;GSE207493;GSE242299" for row in registry))
check("registry_immune_two_exact", any(row["universe_id"] == "D02" and row["datasets_or_samples"] == "GSE264586;GSE294109" for row in registry))
check("registry_immune_bridge_exact", any(row["universe_id"] == "D03" and row["datasets_or_samples"] == "GSE207493;GSE264586;GSE294109" for row in registry))
check("registry_spatial_three_separate", all(key in registry_ids for key in ["D04", "D05", "D06"]))
check("registry_metabolite_separate", all(key in registry_ids for key in ["D09", "D10", "D11", "D12"]))
check("registry_client_separate", any(row["universe_id"] == "D14" and "Do not request, import, analyse, redraw" in row["forbidden_cross_use"] for row in registry))

routes = tsv("PHASE_TO_DATASET_ROUTE_v1.tsv")
route_ids = [row["route_id"] for row in routes]
check("route_rows_44", len(routes) == 44, len(routes))
check("route_ids_unique", len(route_ids) == len(set(route_ids)))
check("route_phase105_superseded", any(row["phase_scope"] == "Phase105" and row["cross_universe_status"] == "SUPERSEDED_CROSS_DATASET_BLUEPRINT" for row in routes))
check("route_phase106_no_pooling", any(row["phase_scope"] == "Phase106" and row["cross_universe_status"] == "NO_BIOLOGICAL_POOLING" for row in routes))
check("route_nci_separate", any(row["phase_scope"] == "Phase102-103" and row["primary_universe_ids"] == "D10" for row in routes))
check(
    "route_spatial_separate",
    all(any(key in row["primary_universe_ids"].split(";") for row in routes) for key in ["D04", "D05", "D06"])
    and all(
        row["cross_universe_status"] in {"SINGLE_UNIVERSE", "SINGLE_SPATIAL_UNIVERSE", "TWO_SPATIAL_COHORTS_KEEP_SEPARATE", "CROSS_SOURCE_ANALYSIS_EXPLICIT", "PREPARATION_STRATIFIED_SINGLE_UNIVERSE", "MULTI_OUTPUT_KEEP_SEPARATE"}
        for row in routes
        if any(key in row["primary_universe_ids"].split(";") for key in ["D04", "D05", "D06"])
    ),
)

claims = tsv("DATASET_SPECIFIC_CLAIM_MATRIX_v1.tsv")
claim_ids = [row["claim_id"] for row in claims]
check("claim_rows_27", len(claims) == 27, len(claims))
check("claim_ids_unique", len(claim_ids) == len(set(claim_ids)))
for row in claims:
    for universe in row["universe_id"].split(";"):
        check(f"claim_{row['claim_id']}_universe_{universe}", universe in registry_ids)
check("claim_part1_not_mta", any(row["claim_id"] == "C001" and "Measured MTA claims" in row["must_not_be_combined_with"] for row in claims))
check("claim_spatial_edges_separate", any(row["claim_id"] == "C015" for row in claims) and any(row["claim_id"] == "C016" for row in claims))
check("claim_nci_exploratory", any(row["claim_id"] == "C022" and row["evidence_decision"] == "EXPLORATORY_DIRECTION_ONLY" for row in claims))
check("claim_client_no_analysis", any(row["claim_id"] == "C027" and row["evidence_decision"] == "OUTSIDE_PROJECT_ANALYSIS" for row in claims))

conflicts = tsv("NONCOMPARABILITY_AND_CONFLICT_REGISTER_v1.tsv")
conflict_ids = [row["issue_id"] for row in conflicts]
check("conflict_rows_12", len(conflicts) == 12, len(conflicts))
check("conflict_ids_unique", len(conflict_ids) == len(set(conflict_ids)))
check("conflict_baseline_perturbation", any(row["issue_id"] == "N08" and row["classification"] == "BASELINE_VERSUS_PERTURBATION" for row in conflicts))
check("conflict_same_accession", any(row["issue_id"] == "N05" and row["classification"] == "OVERLAPPING_DATASET_DIFFERENT_ANALYSIS_UNIVERSE" for row in conflicts))
check("conflict_spatial_edges", any(row["issue_id"] == "N06" and row["classification"] == "DIFFERENT_EDGES_WITHIN_SAME_UNIVERSE" for row in conflicts))

rule = (ROOT / "PHASE106_DATASET_SEPARATION_RULE_v1.md").read_text(encoding="utf-8")
blueprint = (ROOT / "CLIENT_RESPONSE_BLUEPRINT_v2.md").read_text(encoding="utf-8")
summary = (ROOT / "PHASE106_DATASET_SEPARATION_AUDIT_SUMMARY_v1.md").read_text(encoding="utf-8")
check("rule_no_rescue", "must not be used to fill a missing edge, rescue a negative result" in rule)
check("rule_phase105_superseded", "supersede the cross-dataset narrative in the Phase 105 blueprint" in rule)
check("blueprint_immediate_reply_separates", "不再把原第一部分单细胞、后续免疫治疗单细胞、空间转录组和NCI/细胞系结果串成同一条证据链" in blueprint)
check("blueprint_module_a_only", "Use Module A only" in blueprint)
check("blueprint_phase105_superseded", "v2 supersedes" in blueprint)
check("summary_16_universes", "16 explicit data/presentation universes" in summary)
check("summary_27_claims", "27 dataset-specific claims" in summary)
check("summary_no_unified_direction", "must not force a unified direction" in summary)
check("summary_d01_target", "target is D01 original Part 1 single-cell only" in summary)

failed = [row for row in checks if not row["passed"]]
payload = {
    "status": "PASS_PHASE106_DATASET_SEPARATION_VALIDATION" if not failed else "FAIL_PHASE106_DATASET_SEPARATION_VALIDATION",
    "checks_passed": len(checks) - len(failed),
    "checks_total": len(checks),
    "failed_checks": failed,
}
(ROOT / "PHASE106_DATASET_SEPARATION_VALIDATION_v1.json").write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

manifest_rows = []
for path in sorted(item for item in ROOT.iterdir() if item.is_file() and item.name != "FILE_MANIFEST_SHA256_v1.tsv"):
    manifest_rows.append({"file": path.name, "size_bytes": path.stat().st_size, "sha256": sha256(path)})
with (ROOT / "FILE_MANIFEST_SHA256_v1.tsv").open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=["file", "size_bytes", "sha256"], delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(manifest_rows)

print(json.dumps(payload, ensure_ascii=False, indent=2))
raise SystemExit(0 if not failed else 1)
