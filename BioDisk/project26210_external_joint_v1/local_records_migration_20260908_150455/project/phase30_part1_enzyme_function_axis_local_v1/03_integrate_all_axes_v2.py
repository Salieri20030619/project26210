#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path

import pandas as pd


EXPECTED_SUMMARY_SHA256 = "c55ef2bbd85b3c7203f9417ca100710377e28a38dc1f5112a82e2fb09affb45f"
AXES = {
    "TYMP_axis_primary": "primary",
    "RRM2B_axis_primary": "primary",
    "TYMP_axis_broad_sensitivity": "broad_sensitivity",
    "RRM2B_axis_broad_sensitivity": "broad_sensitivity",
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    root = Path(__file__).resolve().parent
    source = root / "results_local_v1" / "enzyme_axis_endpoint_summary_v1.tsv"
    output = root / "all_axis_integration_local_v2"
    part = root / "all_axis_integration_local_v2.part"
    if output.exists() or part.exists():
        raise FileExistsError("refusing to overwrite all-axis integration output")
    if sha256(source) != EXPECTED_SUMMARY_SHA256:
        raise ValueError("validated endpoint-summary hash mismatch")
    part.mkdir(parents=True)
    summary = pd.read_csv(source, sep="\t")
    rows = []
    for axis, role in AXES.items():
        d = summary[summary["feature"] == axis]
        if len(d) != 6:
            raise ValueError(f"incomplete endpoint rows for {axis}")
        endpoint_positive = {
            endpoint: bool(d.loc[d["endpoint"] == endpoint, "positive_support"].iloc[0])
            for endpoint in ["F1_continuous", "F1_dominant_group", "F1_patient_quartile"]
        }
        endpoint_negative = {
            endpoint: bool(d.loc[d["endpoint"] == endpoint, "negative_support"].iloc[0])
            for endpoint in ["F1_continuous", "F1_dominant_group", "F1_patient_quartile"]
        }
        trajectory_positive = int(d.loc[d["endpoint"] == "trajectory_terminal_shift",
                                        "positive_support"].sum())
        trajectory_negative = int(d.loc[d["endpoint"] == "trajectory_terminal_shift",
                                        "negative_support"].sum())
        positive_flags = [*endpoint_positive.values(), trajectory_positive >= 2]
        negative_flags = [*endpoint_negative.values(), trajectory_negative >= 2]
        n_positive = sum(positive_flags)
        n_negative = sum(negative_flags)
        tier = "HIGH_SUPPORT" if n_positive == 4 else "SECONDARY_SUPPORT" if n_positive >= 3 else "LIMITED_OR_UNSUPPORTED"
        rows.append({
            "axis": axis,
            "axis_role": role,
            "F1_continuous_positive": endpoint_positive["F1_continuous"],
            "F1_dominant_group_positive": endpoint_positive["F1_dominant_group"],
            "F1_patient_quartile_positive": endpoint_positive["F1_patient_quartile"],
            "trajectory_positive_convergence_2_of_3": trajectory_positive >= 2,
            "n_positive_trajectory_methods": trajectory_positive,
            "n_positive_endpoints_of_4": n_positive,
            "positive_evidence_tier": tier,
            "F1_continuous_negative": endpoint_negative["F1_continuous"],
            "F1_dominant_group_negative": endpoint_negative["F1_dominant_group"],
            "F1_patient_quartile_negative": endpoint_negative["F1_patient_quartile"],
            "trajectory_negative_convergence_2_of_3": trajectory_negative >= 2,
            "n_negative_trajectory_methods": trajectory_negative,
            "n_negative_endpoints_of_4": n_negative,
            "opposing_evidence": n_negative >= 3,
        })
    integrated = pd.DataFrame(rows)
    integrated.to_csv(part / "enzyme_axis_integrated_evidence_all_axes_v2.tsv", sep="\t", index=False)
    status = {
        "status": "PASS_ALL_AXIS_INTEGRATION",
        "source_summary_sha256": EXPECTED_SUMMARY_SHA256,
        "n_axes": len(integrated),
        "tier_rule": "4_of_4_high; 3_of_4_secondary; 0_to_2_limited_or_unsupported",
    }
    (part / "status.json").write_text(json.dumps(status, indent=2) + "\n", encoding="utf-8")
    os.replace(part, output)


if __name__ == "__main__":
    main()

