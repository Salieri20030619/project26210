#!/usr/bin/env python3
"""Create the authoritative Phase 53 finalization manifest with final paths."""

from __future__ import annotations

import hashlib
import os
from pathlib import Path

import pandas as pd


HERE = Path(__file__).resolve().parent
FINAL = HERE / "finalization_local_v2"
PART = HERE / "finalization_local_v2.part"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_tsv(frame: pd.DataFrame, path: Path) -> None:
    staging = path.with_name(path.name + ".part")
    frame.to_csv(staging, sep="\t", index=False, encoding="utf-8", lineterminator="\n")
    os.replace(staging, path)


def main() -> int:
    if FINAL.exists() or PART.exists():
        raise FileExistsError("Refusing to overwrite existing Phase 53 finalization v2")
    PART.mkdir(parents=True)
    stages = [
        ("initial_screen", HERE / "analysis_local_v1" / "validation_checks_v1.tsv", 8),
        ("exact_small_n", HERE / "analysis_local_v2" / "validation_checks_v2.tsv", 4),
        ("deduplication_overlap_weighting", HERE / "final_evidence_local_v1" / "validation_checks_v1.tsv", 10),
        ("method_consolidation", HERE / "consolidated_findings_local_v1" / "validation_checks_v1.tsv", 8),
    ]
    rows = []
    for stage, path, expected in stages:
        table = pd.read_csv(path, sep="\t")
        passed = int((table["status"] == "PASS").sum())
        rows.append({
            "stage": stage,
            "status": "PASS" if len(table) == expected and passed == expected else "FAIL",
            "passed_checks": passed,
            "total_checks": len(table),
            "expected_checks": expected,
            "validation_file": str(path.relative_to(HERE)).replace("\\", "/"),
            "sha256": sha256(path),
        })
    summary = pd.DataFrame(rows)
    if (summary["status"] != "PASS").any():
        atomic_tsv(summary, PART / "validation_summary_v2.tsv")
        raise RuntimeError("Phase 53 finalization v2 validation failed")
    atomic_tsv(summary, PART / "validation_summary_v2.tsv")
    (PART / "PHASE53_REFERENCE_BRIDGE_SCREEN_PASS_v2").write_text(
        "status=PASS\nvalidation=30/30\nreportable_findings=10\nindependent_reportable_findings=5\nauthoritative_manifest=v2\n",
        encoding="utf-8",
    )

    selected = []
    excluded_dirs = {"__pycache__", "signature_overlap_audit_local_v1.part", "finalization_local_v1", "finalization_local_v2.part"}
    for path in HERE.rglob("*"):
        if not path.is_file() or any(part in excluded_dirs for part in path.parts) or path.suffix == ".pyc":
            continue
        selected.append(path)
    selected.extend([PART / "validation_summary_v2.tsv", PART / "PHASE53_REFERENCE_BRIDGE_SCREEN_PASS_v2"])
    manifest_rows = []
    for path in sorted(set(selected), key=lambda item: str(item)):
        if PART in path.parents:
            rel = Path(FINAL.name) / path.relative_to(PART)
        else:
            rel = path.relative_to(HERE)
        manifest_rows.append({
            "relative_path": str(rel).replace("\\", "/"),
            "bytes": path.stat().st_size,
            "sha256": sha256(path),
        })
    atomic_tsv(pd.DataFrame(manifest_rows), PART / "FILE_MANIFEST_SHA256_v2.tsv")
    os.replace(PART, FINAL)
    print(f"PASS stages={len(summary)} checks={summary['passed_checks'].sum()}/{summary['total_checks'].sum()} files={len(manifest_rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
