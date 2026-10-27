#!/usr/bin/env python3
"""Validate Phase 53 stage markers and write the final SHA256 manifest."""

from __future__ import annotations

import hashlib
import os
from pathlib import Path

import pandas as pd


HERE = Path(__file__).resolve().parent
FINAL = HERE / "finalization_local_v1"
PART = HERE / "finalization_local_v1.part"


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
        raise FileExistsError("Refusing to overwrite existing Phase 53 finalization")
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
        status = "PASS" if len(table) == expected and passed == expected else "FAIL"
        rows.append({
            "stage": stage,
            "status": status,
            "passed_checks": passed,
            "total_checks": len(table),
            "expected_checks": expected,
            "validation_file": str(path.relative_to(HERE)),
            "sha256": sha256(path),
        })
    summary = pd.DataFrame(rows)
    required = [
        HERE / "PHASE53_REFERENCE_BRIDGE_SCREEN_STATUS_v1.md",
        HERE / "INTERPRETATION_v1.md",
        HERE / "REPRODUCIBILITY_v1.md",
        HERE / "signature_overlap_audit_local_v2" / "SIGNATURE_OVERLAP_AUDIT_PASS_v1",
        HERE / "final_evidence_local_v1" / "PHASE53_FINAL_EVIDENCE_PASS_v1",
        HERE / "consolidated_findings_local_v1" / "PHASE53_CONSOLIDATED_FINDINGS_PASS_v1",
    ]
    if (summary["status"] != "PASS").any() or not all(path.exists() for path in required):
        atomic_tsv(summary, PART / "validation_summary_v1.tsv")
        raise RuntimeError("Phase 53 finalization validation failed")
    atomic_tsv(summary, PART / "validation_summary_v1.tsv")
    (PART / "PHASE53_REFERENCE_BRIDGE_SCREEN_PASS_v1").write_text(
        "status=PASS\nvalidation=30/30\nreportable_findings=10\nindependent_reportable_findings=5\n",
        encoding="utf-8",
    )

    selected = []
    for path in HERE.rglob("*"):
        if not path.is_file() or "__pycache__" in path.parts or path.suffix == ".pyc":
            continue
        if "signature_overlap_audit_local_v1.part" in path.parts:
            continue
        if "finalization_local_v1.part" in path.parts:
            continue
        selected.append(path)
    selected.extend([PART / "validation_summary_v1.tsv", PART / "PHASE53_REFERENCE_BRIDGE_SCREEN_PASS_v1"])
    manifest = pd.DataFrame([
        {
            "relative_path": str(path.relative_to(HERE if HERE in path.parents else PART)).replace("\\", "/")
            if HERE in path.parents else path.name,
            "bytes": path.stat().st_size,
            "sha256": sha256(path),
        }
        for path in sorted(set(selected), key=lambda item: str(item))
    ])
    atomic_tsv(manifest, PART / "FILE_MANIFEST_SHA256_v1.tsv")
    os.replace(PART, FINAL)
    print(f"PASS stages={len(summary)} checks={summary['passed_checks'].sum()}/{summary['total_checks'].sum()} files={len(manifest)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
