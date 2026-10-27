#!/usr/bin/env python3
"""Audit Phase 17 GLMF/TJ-RCC mapping queue after all jobs finish."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from datetime import datetime
from pathlib import Path

import pandas as pd


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--root", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    if a.output.exists():
        raise FileExistsError(a.output)
    a.output.mkdir(parents=True)
    groups = {
        "GLMF": (a.root / "glmf12_cell2location_tymp_continuous_v1", 2),
        "TJRCC": (a.root / "tjrcc12_cell2location_tymp_continuous_v1", 12),
    }
    rows = []
    for cohort, (root, expected) in groups.items():
        statuses = sorted(root.glob("*/status.json"))
        for path in statuses:
            value = json.loads(path.read_text(encoding="utf-8"))
            h5ad = Path(value["output_h5ad"]) if value.get("output_h5ad") else None
            rows.append({"cohort": cohort, "sample": value.get("sample", path.parent.name),
                         "status": value.get("status"), "epoch": value.get("epoch"),
                         "output_exists": bool(h5ad and h5ad.is_file()),
                         "output_sha256_matches": bool(h5ad and h5ad.is_file() and sha256(h5ad) == value.get("output_sha256")),
                         "status_file": str(path)})
        if len(statuses) != expected:
            rows.append({"cohort": cohort, "sample": "__cohort_count__", "status": "failed",
                         "epoch": 0, "output_exists": False, "output_sha256_matches": False,
                         "status_file": f"expected {expected}, found {len(statuses)}"})
    table = pd.DataFrame(rows)
    passed = len(table) == 14 and table["status"].eq("completed").all() and table["output_exists"].all() and table["output_sha256_matches"].all()
    table.to_csv(a.output / "mapping_queue_audit_v1.tsv", sep="\t", index=False)
    status = {"schema_version": "phase17_mapping_queue_audit_v1", "status": "PASS" if passed else "FAIL",
              "completed_at": datetime.now().astimezone().isoformat(timespec="seconds"),
              "n_rows": len(table), "n_completed": int(table["status"].eq("completed").sum()),
              "n_failed": int((~table["status"].eq("completed")).sum())}
    (a.output / "PHASE17_MAPPING_QUEUE_STATUS_v1.json").write_text(json.dumps(status, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(status, ensure_ascii=False, indent=2))
    if not passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
