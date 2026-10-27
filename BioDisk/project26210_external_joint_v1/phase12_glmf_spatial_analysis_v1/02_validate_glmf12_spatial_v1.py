#!/usr/bin/env python3
"""Independent structural and conservation checks for Phase 12 outputs."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd


EXPECTED_SPOTS = {"GLMF1": 3_199, "GLMF2": 3_880}
EXPECTED_STATES = 14


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_json(path: Path, value: object) -> None:
    part = Path(str(path) + ".part")
    part.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    os.replace(part, path)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-root", type=Path, required=True)
    args = parser.parse_args()
    root = args.output_root.resolve(strict=True)
    status = json.loads((root / "PHASE12_SPATIAL_ANALYSIS_STATUS.json").read_text(encoding="utf-8"))
    checks: list[dict] = []

    def check(name: str, passed: bool, observed: object) -> None:
        checks.append({"check": name, "pass": bool(passed), "observed": observed})

    check("formal_status_pass", status.get("status") == "PASS", status.get("status"))
    check("no_part_files", not list(root.rglob("*.part")), len(list(root.rglob("*.part"))))
    manifest = pd.read_csv(root / "OUTPUT_MANIFEST.tsv", sep="\t")
    manifest_ok = True
    for row in manifest.itertuples(index=False):
        path = root / row.file
        manifest_ok &= path.is_file() and path.stat().st_size == int(row.size_bytes) and sha256(path) == row.sha256
    check("manifest_files_size_sha256", manifest_ok, len(manifest))

    spot = pd.read_csv(root / "tables/spot_summary.tsv.gz", sep="\t")
    counts = spot.groupby("sample").size().to_dict()
    check("spot_counts", counts == EXPECTED_SPOTS, counts)
    check("spot_ids_unique_within_sample", not spot.duplicated(["sample", "spot_id"]).any(), int(spot.duplicated(["sample", "spot_id"]).sum()))
    numeric_spot = spot[["x", "y", "total_q05_abundance", "dominant_q05_abundance", "dominant_fraction"]].to_numpy(float)
    check("spot_values_finite", np.isfinite(numeric_spot).all(), bool(np.isfinite(numeric_spot).all()))
    check("spot_values_nonnegative", (spot[["total_q05_abundance", "dominant_q05_abundance", "dominant_fraction"]].to_numpy(float) >= 0).all(), float(numeric_spot[:, 2:].min()))

    long = pd.read_csv(root / "tables/spot_state_abundance_long.tsv.gz", sep="\t")
    states = sorted(long["state"].astype(str).unique())
    check("state_count", len(states) == EXPECTED_STATES, len(states))
    check("long_row_count", len(long) == sum(EXPECTED_SPOTS.values()) * EXPECTED_STATES, len(long))
    check("long_keys_unique", not long.duplicated(["sample", "spot_id", "state"]).any(), int(long.duplicated(["sample", "spot_id", "state"]).sum()))
    values = long[["q05_abundance", "within_spot_fraction"]].to_numpy(float)
    check("long_values_finite_nonnegative", np.isfinite(values).all() and (values >= 0).all(), [float(values.min()), float(values.max())])
    fraction_sums = long.groupby(["sample", "spot_id"], sort=False)["within_spot_fraction"].sum().to_numpy()
    check("within_spot_fractions_sum_one", np.allclose(fraction_sums, 1.0, atol=2e-6), float(np.max(np.abs(fraction_sums - 1))))

    summary = pd.read_csv(root / "tables/sample_state_summary.tsv", sep="\t")
    moran = pd.read_csv(root / "tables/spatial_autocorrelation_moran.tsv", sep="\t")
    corr = pd.read_csv(root / "tables/state_correlations_long.tsv.gz", sep="\t")
    check("summary_rows", len(summary) == 2 * EXPECTED_STATES, len(summary))
    check("moran_rows", len(moran) == 2 * EXPECTED_STATES, len(moran))
    check("correlation_rows", len(corr) == 2 * 2 * EXPECTED_STATES * EXPECTED_STATES, len(corr))
    check("summary_state_set", sorted(summary["state"].astype(str).unique()) == states, len(summary["state"].unique()))

    expected_figures = [
        "sample_composition.png",
        "GLMF1_dominant_state.png",
        "GLMF1_state_abundance_panels.png",
        "GLMF2_dominant_state.png",
        "GLMF2_state_abundance_panels.png",
        "same_spot_correlations.png",
        "neighbor_lag_correlations.png",
    ]
    figure_sizes = {name: (root / "figures" / name).stat().st_size if (root / "figures" / name).is_file() else 0 for name in expected_figures}
    check("figures_present_nonempty", all(size > 0 for size in figure_sizes.values()), figure_sizes)

    failures = [item for item in checks if not item["pass"]]
    result = {
        "status": "PASS" if not failures else "FAIL",
        "validated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "checks_passed": len(checks) - len(failures),
        "checks_total": len(checks),
        "checks": checks,
    }
    atomic_json(root / "INDEPENDENT_VALIDATION_STATUS.json", result)
    print(json.dumps(result, ensure_ascii=False, indent=2))
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
