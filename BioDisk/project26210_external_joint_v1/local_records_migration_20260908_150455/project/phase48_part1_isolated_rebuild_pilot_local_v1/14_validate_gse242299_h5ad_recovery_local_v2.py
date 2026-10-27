#!/usr/bin/env python3
import csv
import json
import re
from datetime import datetime
from pathlib import Path


PHASE_ROOT = Path(__file__).resolve().parent
SYNC_ROOT = PHASE_ROOT / "server_results_sync_v4"
CHECKS_PATH = PHASE_ROOT / "gse242299_h5ad_recovery_validation_checks_v2.tsv"
STATUS_PATH = PHASE_ROOT / "GSE242299_H5AD_RECOVERY_LOCAL_VALIDATION_PASS_v2.md"
checks = []


def add_check(name, passed, observed, expected):
    checks.append({
        "check": name,
        "status": "PASS" if passed else "FAIL",
        "observed": str(observed),
        "expected": str(expected),
    })


def read_tsv(name):
    with (SYNC_ROOT / name).open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


required = [
    "GSE242299_H5AD_RECOVERY_PASS_v3",
    "H5AD_SOURCE_MANIFEST_v3.tsv",
    "h5ad_dataset_inventory_v3.tsv",
    "h5ad_library_cell_counts_v3.tsv",
    "h5ad_matrix_candidate_audit_v3.tsv",
    "h5ad_obs_candidate_columns_v3.tsv",
    "h5ad_recovery_metadata_v3.json",
    "11_finalize_gse242299_h5ad_recovery_v3.log",
]
for name in required:
    add_check(f"required_file:{name}", (SYNC_ROOT / name).is_file(), "evaluated", "present")

marker = (SYNC_ROOT / required[0]).read_text(encoding="utf-8-sig")
metadata = json.loads((SYNC_ROOT / "h5ad_recovery_metadata_v3.json").read_text(encoding="utf-8-sig"))
matrix = read_tsv("h5ad_matrix_candidate_audit_v3.tsv")
libraries = read_tsv("h5ad_library_cell_counts_v3.tsv")
source = read_tsv("H5AD_SOURCE_MANIFEST_v3.tsv")
log = (SYNC_ROOT / "11_finalize_gse242299_h5ad_recovery_v3.log").read_text(encoding="utf-8-sig")

add_check("marker_status", bool(re.search(r"(?m)^status=PASS\s*$", marker)), marker.splitlines()[0], "status=PASS")
add_check("metadata_status", metadata["status"] == "PASS", metadata["status"], "PASS")
add_check("dimensions", (metadata["n_obs"], metadata["n_vars"]) == (50236, 33538), f'{metadata["n_obs"]}x{metadata["n_vars"]}', "50236x33538")
add_check("library_column", metadata["selected_library_column"] == "library_name" and metadata["exact_library_columns"] == ["library_name"], metadata["exact_library_columns"], "library_name as sole exact 31-library field")
add_check("selected_count_matrix", metadata["selected_count_matrix_path"] == "/X" and metadata["count_matrix_rule"] == "accept_unique_integer_X_when_raw_and_layers_absent", f'{metadata["selected_count_matrix_path"]};{metadata["count_matrix_rule"]}', "/X; unique integer X with raw/layers absent")
add_check("author_code_commit", metadata["author_code_commit"] == "1c9daebb33449f32a06978227aa5fcf0db726c3a", metadata["author_code_commit"], "1c9daebb33449f32a06978227aa5fcf0db726c3a")

add_check("matrix_row_count", len(matrix) == 1, len(matrix), 1)
if len(matrix) == 1:
    m = matrix[0]
    matrix_pass = (
        m["matrix_path"] == "/X" and m["storage"] == "csr_matrix"
        and int(m["shape_rows"]) == 50236 and int(m["shape_columns"]) == 33538
        and int(m["stored_values"]) == 26714153
        and all(m[key] == "true" for key in ("finite_all", "nonnegative_all", "integer_all", "structural_checks", "raw_count_candidate"))
        and float(m["minimum"]) == 1 and float(m["maximum"]) == 1251
    )
    add_check("matrix_safety_and_structure", matrix_pass, f'{m["matrix_path"]};{m["shape_rows"]}x{m["shape_columns"]};nnz={m["stored_values"]};range={m["minimum"]}-{m["maximum"]}', "/X CSR; exact dimensions/nnz; all safety checks true; range 1-1251")

library_ids = [row["library_id"] for row in libraries]
cell_sum = sum(int(row["cells"]) for row in libraries)
cell_lookup = {row["library_id"]: int(row["cells"]) for row in libraries}
add_check("library_count_and_uniqueness", len(libraries) == 31 and len(set(library_ids)) == 31, f"rows={len(libraries)};unique={len(set(library_ids))}", "rows=31;unique=31")
add_check("library_cell_conservation", cell_sum == 50236 and all(row["status"] == "PASS" for row in libraries), f"sum={cell_sum};failures={sum(row['status'] != 'PASS' for row in libraries)}", "sum=50236;failures=0")
add_check("affected_library_recovery", cell_lookup.get("T2_1") == 4906 and cell_lookup.get("N3_1") == 240, f'T2_1={cell_lookup.get("T2_1")};N3_1={cell_lookup.get("N3_1")}', "T2_1=4906;N3_1=240")

add_check("source_manifest_row_count", len(source) == 1, len(source), 1)
if len(source) == 1:
    s = source[0]
    source_pass = (
        s["status"] == "PASS" and int(s["decompressed_size_bytes"]) == 7345987312
        and bool(re.fullmatch(r"[0-9a-f]{64}", s["decompressed_sha256"]))
        and int(s["compressed_size_bytes"]) == 5350598000
        and s["compressed_sha256"] == "8539723bf40f3d8647c9773c386db56256b305f69e0cdf80848cec5bae428f94"
    )
    add_check("source_manifest_integrity", source_pass, f'decompressed={s["decompressed_size_bytes"]}:{s["decompressed_sha256"]};compressed={s["compressed_size_bytes"]}:{s["compressed_sha256"]};status={s["status"]}', "validated sizes and hashes; PASS")
add_check("atomic_finalize_log", "finished_at=" in log and not any(token in log for token in ("Traceback", "RuntimeError", "BLOCKED")), "finish/error scan evaluated", "finished_at present; no error token")

with CHECKS_PATH.open("x", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=["check", "status", "observed", "expected"], delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(checks)
failed = [row for row in checks if row["status"] != "PASS"]
if failed:
    raise RuntimeError("GSE242299 H5AD independent validation failed: " + ", ".join(row["check"] for row in failed))

source_row = source[0]
status = f"""# GSE242299 official H5AD recovery — independent local validation

- Status: **PASS**
- Generated: {datetime.now().astimezone().isoformat(timespec='seconds')}
- Checks: {len(checks)}/{len(checks)} PASS
- Matrix: `/X`, CSR, 50,236 cells × 33,538 genes, 26,714,153 stored nonzero values
- Count safety: finite, nonnegative, exactly integer-valued, structurally valid
- Library field: `library_name`, exact 31-library match; cell total 50,236
- Recovered defective-MTX libraries: T2_1 = 4,906 cells; N3_1 = 240 cells
- Official compressed source SHA256: `{source_row['compressed_sha256']}`
- Decompressed H5AD SHA256: `{source_row['decompressed_sha256']}`
- Author-code snapshot commit: `{metadata['author_code_commit']}`
- Large H5AD remains on `sgs10`; only small audit evidence was synchronized locally.
"""
STATUS_PATH.write_text(status, encoding="utf-8")
print(f"PASS checks={len(checks)} status={STATUS_PATH}")
