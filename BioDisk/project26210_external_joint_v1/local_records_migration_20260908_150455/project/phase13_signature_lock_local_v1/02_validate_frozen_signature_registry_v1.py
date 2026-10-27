#!/usr/bin/env python3
"""Independent structural/hash validation for the Phase 13 registry."""

from __future__ import annotations

import csv
import hashlib
import json
import os
from datetime import datetime
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
REG = HERE / "registry_v1"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def atomic_tsv(path: Path, fieldnames: list[str], rows: list[dict]) -> None:
    part = path.with_name(path.name + ".part")
    if path.exists() or part.exists():
        raise FileExistsError(f"Refusing existing output: {path}")
    with part.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(part, path)


def main() -> None:
    expected = [
        REG / "frozen_signature_registry_v1.tsv",
        REG / "frozen_signature_genes_v1.tsv",
        REG / "frozen_analysis_rules_v1.json",
        REG / "source_manifest_v1.tsv",
        REG / "run_metadata_v1.json",
    ]
    checks: list[dict] = []

    def check(name: str, passed: bool, detail: str) -> None:
        checks.append({"check": name, "result": "PASS" if passed else "FAIL", "detail": detail})

    check("all_registry_files_readable", all(path.is_file() and path.stat().st_size > 0 for path in expected), ";".join(path.name for path in expected))
    registry = read_tsv(expected[0])
    genes = read_tsv(expected[1])
    sources = read_tsv(expected[3])
    rules = json.loads(expected[2].read_text(encoding="utf-8"))

    signature_ids = [row["signature_id"] for row in registry]
    check("signature_ids_unique", len(signature_ids) == len(set(signature_ids)), f"n={len(signature_ids)}")
    gene_pairs = [(row["signature_id"], row["gene"]) for row in genes]
    check("genes_unique_within_signature", len(gene_pairs) == len(set(gene_pairs)), f"n={len(gene_pairs)}")
    counts = {signature_id: 0 for signature_id in signature_ids}
    for row in genes:
        counts[row["signature_id"]] = counts.get(row["signature_id"], 0) + 1
    count_match = all(counts.get(row["signature_id"]) == int(row["n_genes"]) for row in registry)
    check("registry_gene_counts_conserved", count_match, json.dumps(counts, sort_keys=True))

    f1 = [row for row in genes if row["signature_id"] == "f1_top100_weighted"]
    f1_ok = len(f1) == 100 and all(float(row["raw_weight"]) >= 0 for row in f1) and len({row["gene"] for row in f1}) == 100
    check("f1_top100_nonnegative_unique", f1_ok, f"n={len(f1)}")
    check("outcome_blind_lock_declared", rules.get("outcome_blind") is True, str(rules.get("outcome_blind")))
    check("patient_primary_unit_declared", "patient" in rules.get("primary_unit", ""), rules.get("primary_unit", ""))
    check("datasets_separate_declared", rules.get("datasets_are_separate") is True, str(rules.get("datasets_are_separate")))

    source_failures = []
    for row in sources:
        path = ROOT / row["path_from_project_root"]
        if not path.is_file() or path.stat().st_size != int(row["bytes"]) or sha256(path).lower() != row["sha256"].lower():
            source_failures.append(row["path_from_project_root"])
    check("source_hashes_and_sizes_match", not source_failures, ";".join(source_failures) if source_failures else f"n={len(sources)}")

    overall = all(row["result"] == "PASS" for row in checks)
    checks.append({"check": "OVERALL", "result": "PASS" if overall else "FAIL", "detail": datetime.now().astimezone().isoformat(timespec="seconds")})
    atomic_tsv(REG / "VALIDATION_STATUS_v1.tsv", ["check", "result", "detail"], checks)
    if not overall:
        raise SystemExit("Phase13 registry validation failed")


if __name__ == "__main__":
    main()
