#!/usr/bin/env python3
"""Create review packets for current-annotation versus rebuild-CNV disagreements."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def atomic_tsv(frame: pd.DataFrame, final: Path, compression: str | None = None) -> None:
    if final.exists():
        raise FileExistsError(f"No-overwrite gate: {final}")
    part = final.with_name(final.name + ".part")
    frame.to_csv(part, sep="\t", index=False, na_rep="NA", compression=compression)
    os.replace(part, final)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gate", required=True)
    parser.add_argument("--out-dir", required=True)
    args = parser.parse_args()
    gate_path = Path(args.gate).resolve()
    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    cell_path = out_dir / "annotation_cnv_conflict_cells_v1.tsv.gz"
    summary_path = out_dir / "annotation_cnv_conflict_summary_v1.tsv"
    status_path = out_dir / "ANNOTATION_CNV_CONFLICT_REVIEW_STATUS"
    for path in (cell_path, summary_path, status_path):
        if path.exists():
            raise FileExistsError(f"No-overwrite gate: {path}")

    cells = pd.read_csv(gate_path, sep="\t", compression="infer", dtype=str)
    if len(cells) != 585119 or cells["cell_id"].duplicated().any():
        raise RuntimeError("All-cell gate conservation failed")
    conflict = cells.loc[cells["phase9_requires_annotation_review"].str.lower().eq("true")].copy()
    if len(conflict) != 696:
        raise RuntimeError(f"Expected 696 conflicts, found {len(conflict)}")
    conflict["review_priority"] = "medium_scope_disagreement"
    malignant = conflict["part1__malignant_conservative_v1"].str.lower().eq("true")
    conflict.loc[malignant, "review_priority"] = "high_part1_conservative_malignant"
    summary = (
        conflict.groupby(
            [
                "dataset", "phase7_broad_annotation", "phase7_fine_annotation",
                "phase8_spatial_reference_state", "part1__cell_type_broad_v1",
                "part1__malignancy_status_v1", "part1__malignant_conservative_v1",
                "review_priority",
            ],
            observed=True,
            dropna=False,
        )
        .size()
        .reset_index(name="n_cells")
        .sort_values(["review_priority", "n_cells"], ascending=[True, False])
    )
    atomic_tsv(conflict, cell_path, compression="gzip")
    atomic_tsv(summary, summary_path)
    counts = {
        "all_cells": len(cells),
        "conflict_cells": len(conflict),
        "part1_conservative_malignant_conflicts": int(malignant.sum()),
        "datasets": int(conflict["dataset"].nunique()),
        "patients": int(conflict[["dataset", "patient_id"]].drop_duplicates().shape[0]),
    }
    status = {
        "status": "PASS",
        "completed_at_utc": datetime.now(timezone.utc).isoformat(),
        "counts": counts,
        "input": {"path": str(gate_path), "bytes": gate_path.stat().st_size, "sha256": sha256(gate_path)},
        "outputs": {
            "cells": {"path": str(cell_path), "bytes": cell_path.stat().st_size, "sha256": sha256(cell_path)},
            "summary": {"path": str(summary_path), "bytes": summary_path.stat().st_size, "sha256": sha256(summary_path)},
        },
        "interpretation": "Conflicts are review flags. They are not automatically converted to malignant or deleted.",
    }
    part = status_path.with_name(status_path.name + ".part")
    part.write_text(json.dumps(status, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(part, status_path)
    print(json.dumps(counts, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
