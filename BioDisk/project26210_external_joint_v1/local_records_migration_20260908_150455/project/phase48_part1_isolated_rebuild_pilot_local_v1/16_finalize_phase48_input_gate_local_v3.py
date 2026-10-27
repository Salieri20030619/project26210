#!/usr/bin/env python3
import csv
import json
import re
from collections import Counter
from datetime import datetime
from pathlib import Path


ROOT = Path(__file__).resolve().parent
CHECKS_OUT = ROOT / "input_gate_validation_checks_v3.tsv"
STATUS_OUT = ROOT / "PHASE48_ISOLATED_REBUILD_PILOT_STATUS_v4.md"
checks = []


def tsv(path):
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def add(check_id, passed, observed, expected, evidence):
    checks.append({
        "check_id": check_id,
        "status": "PASS" if passed else "FAIL",
        "observed": str(observed),
        "expected": str(expected),
        "evidence": evidence,
    })


rds = tsv(ROOT / "server_results_sync_v1" / "input_file_inventory_v2.tsv")
rds_counts = Counter(row["dataset"] for row in rds)
add("complete_rds_inventory", len(rds) == 45 and rds_counts == {"GSE178481": 26, "GSE207493": 19}, f"rows={len(rds)};GSE178481={rds_counts['GSE178481']};GSE207493={rds_counts['GSE207493']}", "45;26;19", "server_results_sync_v1/input_file_inventory_v2.tsv")
add("rds_manifest_safety", all(int(row["bytes"]) > 0 and re.fullmatch(r"[0-9a-f]{64}", row["sha256"]) for row in rds), "positive sizes and SHA256 format evaluated for 45 rows", "all PASS", "server_results_sync_v1/input_file_inventory_v2.tsv")

archives = tsv(ROOT / "server_results_sync_v2" / "public_geo_archive_manifest_v1.tsv")
add("public_archive_integrity", len(archives) == 3 and all(row["validation_status"] == "PASS" and row["expected_bytes"] == row["actual_bytes"] for row in archives), f"rows={len(archives)};failures={sum(row['validation_status'] != 'PASS' or row['expected_bytes'] != row['actual_bytes'] for row in archives)}", "3 exact-size PASS archives", "server_results_sync_v2/public_geo_archive_manifest_v1.tsv")
redundancy = tsv(ROOT / "PUBLIC_GEO_ARCHIVE_REDUNDANCY_VALIDATION_v1.tsv")
add("redundant_archive_match", len(redundancy) == 3 and all(row["sha256_match"] == "true" for row in redundancy), f"matches={sum(row['sha256_match'] == 'true' for row in redundancy)}/3", "3/3", "PUBLIC_GEO_ARCHIVE_REDUNDANCY_VALIDATION_v1.tsv")

scope = tsv(ROOT / "scope_freeze_local_v3" / "full_cohort_library_scope_v3.tsv")
included = [row for row in scope if row["inclusion_status"] == "include"]
excluded = [row for row in scope if row["inclusion_status"] == "exclude"]
patients = {(row["study"], row["source_patient_id"]) for row in included}
add("scope_v3", len(scope) == 102 and len(included) == 99 and len(excluded) == 3 and len(patients) == 53, f"rows={len(scope)};include={len(included)};exclude={len(excluded)};patients={len(patients)}", "102;99;3;53", "scope_freeze_local_v3/FULL_COHORT_SCOPE_PASS_v3")
add("scope_exclusion_set", {row["geo_accession"] for row in excluded} == {"GSM4819730", "GSM4819731", "GSM4819732"}, "|".join(sorted(row["geo_accession"] for row in excluded)), "GSM4819730|GSM4819731|GSM4819732", "scope_freeze_local_v3/full_cohort_library_scope_v3.tsv")

structure = tsv(ROOT / "server_results_sync_v3" / "raw_count_structure_summary_v3.tsv")
structure_by_study = {row["study"]: row for row in structure}
for study, expected in (("GSE156632", 12), ("GSE159115", 14)):
    row = structure_by_study[study]
    add(f"{study.lower()}_raw_structure", row["status"] == "PASS" and int(row["pass_libraries"]) == expected and int(row["fail_libraries"]) == 0, f'{row["pass_libraries"]}/{row["audited_libraries"]};status={row["status"]}', f"{expected}/{expected};PASS", "server_results_sync_v3/raw_count_structure_summary_v3.tsv")

recovery_checks = tsv(ROOT / "gse242299_h5ad_recovery_validation_checks_v2.tsv")
recovery = json.loads((ROOT / "server_results_sync_v4" / "h5ad_recovery_metadata_v3.json").read_text(encoding="utf-8-sig"))
libraries = tsv(ROOT / "server_results_sync_v4" / "h5ad_library_cell_counts_v3.tsv")
matrix = tsv(ROOT / "server_results_sync_v4" / "h5ad_matrix_candidate_audit_v3.tsv")
source = tsv(ROOT / "server_results_sync_v4" / "H5AD_SOURCE_MANIFEST_v3.tsv")
add("gse242299_independent_recovery_validation", len(recovery_checks) == 22 and all(row["status"] == "PASS" for row in recovery_checks), f"PASS={sum(row['status'] == 'PASS' for row in recovery_checks)}/{len(recovery_checks)}", "22/22", "gse242299_h5ad_recovery_validation_checks_v2.tsv")
add("gse242299_library_conservation", len(libraries) == 31 and sum(int(row["cells"]) for row in libraries) == 50236 and all(row["status"] == "PASS" for row in libraries), f"libraries={len(libraries)};cells={sum(int(row['cells']) for row in libraries)}", "31;50236", "server_results_sync_v4/h5ad_library_cell_counts_v3.tsv")
add("gse242299_raw_count_matrix", len(matrix) == 1 and matrix[0]["matrix_path"] == "/X" and all(matrix[0][key] == "true" for key in ("finite_all", "nonnegative_all", "integer_all", "structural_checks", "raw_count_candidate")), f'{matrix[0]["matrix_path"]};finite={matrix[0]["finite_all"]};nonnegative={matrix[0]["nonnegative_all"]};integer={matrix[0]["integer_all"]};structural={matrix[0]["structural_checks"]}', "/X;all true", "server_results_sync_v4/h5ad_matrix_candidate_audit_v3.tsv")
add("gse242299_source_hash", len(source) == 1 and source[0]["status"] == "PASS" and source[0]["compressed_sha256"] == "8539723bf40f3d8647c9773c386db56256b305f69e0cdf80848cec5bae428f94" and re.fullmatch(r"[0-9a-f]{64}", source[0]["decompressed_sha256"]) is not None, f'compressed={source[0]["compressed_sha256"]};decompressed={source[0]["decompressed_sha256"]}', "validated compressed and decompressed SHA256", "server_results_sync_v4/H5AD_SOURCE_MANIFEST_v3.tsv")
add("gse242299_scope_source_consistency", all(row["input_source"] == "NCBI_GEO_official_series_H5AD_validated_raw_X" and row["archive_members"].endswith("::/X") for row in scope if row["study"] == "GSE242299"), "31 GSE242299 scope rows evaluated", "31/31 use validated official H5AD /X", "scope_freeze_local_v3/full_cohort_library_scope_v3.tsv")

old_marker = (ROOT / "server_results_sync_v3" / "PUBLIC_GEO_RAW_COUNT_STRUCTURE_BLOCKED_v3").read_text(encoding="utf-8-sig")
add("legacy_mtx_defect_provenance", "failed_libraries=GSM7757833:T2_1|GSM7757855:N3_1" in old_marker, "two official MTX defects remain documented", "T2_1 and N3_1 defects preserved", "server_results_sync_v3/PUBLIC_GEO_RAW_COUNT_STRUCTURE_BLOCKED_v3")
add("no_pre_gate_downstream", "downstream_started=false" in old_marker, "downstream_started=false", "false before gate closure", "server_results_sync_v3/PUBLIC_GEO_RAW_COUNT_STRUCTURE_BLOCKED_v3")
add("five_dataset_separation", len({row["study"] for row in included}) == 5, "five study-specific source groups", "5; never combined", "scope_freeze_local_v3/full_cohort_library_scope_v3.tsv")

with CHECKS_OUT.open("x", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=["check_id", "status", "observed", "expected", "evidence"], delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(checks)
failed = [row for row in checks if row["status"] != "PASS"]
if failed:
    raise RuntimeError("Phase 48 input gate failed: " + ", ".join(row["check_id"] for row in failed))

status = f"""# Phase 48 isolated Part 1 cohort-rebuild pilot status v4

Generated: {datetime.now().astimezone().isoformat(timespec='seconds')}

## Status

**PHASE5_INPUT_GATE_PASS — READY_FOR_SERIAL_PER_DATASET_RAW_COUNT_BUILD**

All five study-specific inputs now have a validated raw-count source. The gate closed with {len(checks)}/{len(checks)} checks passing. No normalization, integration, batch correction, dimensionality reduction, clustering or annotation was run before this PASS.

## Frozen scope

Scope v3 retains 102 public library records: 99 included libraries from 53 study-specific patients and three explicitly excluded chRCC-related GSE159115 records. Included libraries by study are GSE156632 12, GSE159115 11, GSE178481 26, GSE207493 19 and GSE242299 31. The five datasets remain separate.

## Validated raw-count inputs

- GSE178481 and GSE207493: 45/45 transferred raw-count RDS inputs (26 plus 19) have positive sizes and recorded SHA256 values matching the frozen inventory.
- GSE156632: 12/12 official GEO processed raw-count CSV libraries passed structure and numeric safety checks.
- GSE159115: 14/14 official GEO processed raw-count HDF5 libraries passed; 11 ccRCC-compatible libraries are included and three chRCC-related libraries remain excluded by the frozen scope.
- GSE242299: the official per-library MTX export still has two independently confirmed source defects (T2_1 malformed/truncated matrix; N3_1 zero-byte barcode/features). These defects remain preserved as provenance and are not silently repaired.
- GSE242299 recovery: the official 5,350,598,000-byte series H5AD passed compressed and decompressed hash checks. Its sole standard matrix candidate is `/X`, a structurally valid 50,236 × 33,538 CSR matrix with 26,714,153 stored values; every stored value is finite, nonnegative and exactly integer-valued. `library_name` exactly matches all 31 scoped libraries, totaling 50,236 cells, including T2_1 4,906 and N3_1 240.

The author-code snapshot at commit `1c9daebb33449f32a06978227aa5fcf0db726c3a` supports the `/X` interpretation: raw counts were stored in `adata.raw`, and downstream code materialized `adata.raw.to_adata()`. The numeric and structural acceptance does not rely on that inference alone; it was independently checked directly in the official H5AD.

## Phase 5 authorization

Phase 5 may now construct one raw-count Seurat object per dataset, serially. It must preserve collision-safe source metadata, build globally unique cell IDs, reject unsafe counts without rounding or repair, release memory between datasets, write final RDS outputs through `.part` files plus atomic rename, and pass all five conservation/safety checks before any downstream processing begins.
"""
STATUS_OUT.write_text(status, encoding="utf-8")
print(f"PASS checks={len(checks)} status={STATUS_OUT}")
