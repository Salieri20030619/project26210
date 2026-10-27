#!/usr/bin/env python3
"""Validate the compact server sync and build the local Phase 54 manifest."""

from __future__ import annotations

import hashlib
import os
from pathlib import Path

import pandas as pd


HERE = Path(__file__).resolve().parent
SYNC = HERE / "server_summary_sync_v1"
FINAL = HERE / "local_finalization_v1"
PART = HERE / "local_finalization_v1.part"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_tsv(frame: pd.DataFrame, path: Path) -> None:
    temporary = path.with_name(path.name + ".part")
    frame.to_csv(temporary, sep="\t", index=False, encoding="utf-8", lineterminator="\n")
    os.replace(temporary, path)


def main() -> int:
    if FINAL.exists() or PART.exists():
        raise FileExistsError("Refusing to overwrite local finalization")
    PART.mkdir(parents=True)
    remote_manifest = pd.read_csv(
        SYNC / "postrun_validation_v1" / "compact_deliverables_manifest_v1.tsv", sep="\t"
    )
    rows = []
    for record in remote_manifest.itertuples(index=False):
        local = SYNC / Path(record.relative_path).name
        observed = sha256(local) if local.exists() else "MISSING"
        rows.append({
            "relative_path_on_server": record.relative_path,
            "local_file": str(local.relative_to(HERE)).replace("\\", "/"),
            "status": "PASS" if observed == record.sha256 else "FAIL",
            "expected_sha256": record.sha256,
            "observed_sha256": observed,
        })
    transfer = pd.DataFrame(rows)
    postrun = pd.read_csv(
        SYNC / "postrun_validation_v1" / "postrun_validation_checks_v1.tsv", sep="\t"
    )
    checks = pd.DataFrame([
        {"check": "compact_transfer_hashes", "status": "PASS" if (transfer["status"] == "PASS").all() and len(transfer) == 10 else "FAIL", "observed": int((transfer["status"] == "PASS").sum()), "expected": 10},
        {"check": "server_postrun_validation", "status": "PASS" if (postrun["status"] == "PASS").all() and len(postrun) == 12 else "FAIL", "observed": int((postrun["status"] == "PASS").sum()), "expected": 12},
        {"check": "section_status_rows", "status": "PASS" if len(pd.read_csv(SYNC / "section_status_v1.tsv", sep="\t")) == 24 else "FAIL", "observed": len(pd.read_csv(SYNC / "section_status_v1.tsv", sep="\t")), "expected": 24},
        {"check": "failure_table_empty", "status": "PASS" if len(pd.read_csv(SYNC / "section_failures_v1.tsv", sep="\t")) == 0 else "FAIL", "observed": len(pd.read_csv(SYNC / "section_failures_v1.tsv", sep="\t")), "expected": 0},
    ])
    atomic_tsv(transfer, PART / "compact_transfer_hash_validation_v1.tsv")
    atomic_tsv(checks, PART / "local_sync_validation_checks_v1.tsv")
    if (checks["status"] != "PASS").any():
        raise RuntimeError("Local compact sync validation failed")
    (PART / "PHASE54_LOCAL_SYNC_PASS_v1").write_text(
        "status=PASS\ntransfer_hashes=10/10\nserver_postrun_validation=12/12\nsections=24/24\n",
        encoding="utf-8",
    )

    selected = []
    for path in HERE.rglob("*"):
        if not path.is_file() or "__pycache__" in path.parts or "local_finalization_v1.part" in path.parts:
            continue
        selected.append(path)
    selected.extend([
        PART / "compact_transfer_hash_validation_v1.tsv",
        PART / "local_sync_validation_checks_v1.tsv",
        PART / "PHASE54_LOCAL_SYNC_PASS_v1",
    ])
    manifest_rows = []
    for path in sorted(set(selected), key=lambda item: str(item)):
        if PART in path.parents:
            relative = Path(FINAL.name) / path.relative_to(PART)
        else:
            relative = path.relative_to(HERE)
        manifest_rows.append({
            "relative_path": str(relative).replace("\\", "/"),
            "bytes": path.stat().st_size,
            "sha256": sha256(path),
        })
    atomic_tsv(pd.DataFrame(manifest_rows), PART / "FILE_MANIFEST_SHA256_v1.tsv")
    os.replace(PART, FINAL)
    print(f"PASS transfer=10/10 postrun=12/12 files={len(manifest_rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

