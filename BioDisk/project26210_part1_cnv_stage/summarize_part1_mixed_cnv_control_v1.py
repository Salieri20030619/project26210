#!/usr/bin/env python3
"""Resolve inferCNV rerun confounding with matched no-add 4-thread controls."""

from __future__ import annotations

import hashlib
import json
from datetime import datetime
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt

import summarize_part1_mixed_cnv_sensitivity_v1 as base


PROJECT = Path("/mnt/BioDisk/ludingyu/project26210")
ROOT = PROJECT / "results/part1_rebuild_v1/03_malignancy_cnv"
SENS = ROOT / "mixed_cnv_sensitivity_v1"
CONTROL = SENS / "control_no_add_4threads"
SUMMARY_V1 = SENS / "summary_v1"
OUT = SENS / "summary_v2_with_control"
PRIMARY = ROOT / "malignancy_cnv_summary_v1/infercnv_query_cell_calls_v1.tsv.gz"
UNITS = ["GSE156632__RCC3", "GSE156632__RCC5"]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def compare_calls(
    unit: str,
    left_name: str,
    left: pd.DataFrame,
    right_name: str,
    right: pd.DataFrame,
) -> tuple[dict[str, object], pd.DataFrame]:
    l = left.loc[left["infercnv_unit"].eq(unit), ["exact_cell_id", "infercnv_call"]].copy()
    r = right.loc[right["infercnv_unit"].eq(unit), ["exact_cell_id", "infercnv_call"]].copy()
    merged = l.merge(r, on="exact_cell_id", suffixes=("_left", "_right"), validate="one_to_one")
    if len(merged) != len(l) or len(r) != len(l):
        raise RuntimeError(f"{unit}: ID registry differs for {left_name} vs {right_name}")
    same = merged["infercnv_call_left"].eq(merged["infercnv_call_right"])
    summary = {
        "unit_id": unit,
        "left_run": left_name,
        "right_run": right_name,
        "n_clear_query": len(merged),
        "n_exact_call_agreement": int(same.sum()),
        "n_call_changed": int((~same).sum()),
        "exact_call_agreement_fraction": float(same.mean()),
    }
    transition = merged.groupby(
        ["infercnv_call_left", "infercnv_call_right"], observed=True
    ).size().rename("n_cells").reset_index().rename(
        columns={"infercnv_call_left": "left_call", "infercnv_call_right": "right_call"}
    )
    transition.insert(0, "right_run", right_name)
    transition.insert(0, "left_run", left_name)
    transition.insert(0, "unit_id", unit)
    return summary, transition


def main() -> None:
    if not (CONTROL / "RUN_COMPLETE").is_file():
        raise RuntimeError("No-add inferCNV controls are incomplete")
    if OUT.exists() and any(OUT.iterdir()):
        raise RuntimeError(f"Refusing to overwrite non-empty control summary: {OUT}")
    OUT.mkdir(parents=True, exist_ok=True)
    gates: list[dict[str, str]] = []
    previous_root = base.SENS
    base.SENS = CONTROL
    try:
        control_parts = []
        burden_parts = []
        for unit in UNITS:
            cells, burden = base.read_infer(unit, gates)
            control_parts.append(cells)
            burden_parts.append(burden)
    finally:
        base.SENS = previous_root
    control = pd.concat(control_parts, ignore_index=True)
    control_burden = pd.concat(burden_parts, ignore_index=True)
    primary = pd.read_csv(PRIMARY, sep="\t", dtype=str)
    plus_one = pd.read_csv(
        SUMMARY_V1 / "infercnv_query_cell_calls_mixed_sensitivity_v1.tsv.gz", sep="\t", dtype=str
    )
    plus_one_clear = plus_one.loc[plus_one["cell_type_broad_v1"].eq("Epithelial cell")].copy()

    comparisons = []
    transitions = []
    for unit in UNITS:
        for left_name, left, right_name, right in [
            ("primary_8threads", primary, "control_no_add_4threads", control),
            ("control_no_add_4threads", control, "plus_one_mixed_4threads", plus_one_clear),
            ("primary_8threads", primary, "plus_one_mixed_4threads", plus_one_clear),
        ]:
            summary, transition = compare_calls(unit, left_name, left, right_name, right)
            comparisons.append(summary)
            transitions.append(transition)
    comparison = pd.DataFrame(comparisons)
    transition = pd.concat(transitions, ignore_index=True)

    threshold_rows = []
    for label, burden in [
        ("control_no_add_4threads", control_burden),
        (
            "plus_one_mixed_4threads",
            pd.read_csv(SUMMARY_V1 / "infercnv_group_burden_mixed_sensitivity_v1.tsv", sep="\t"),
        ),
    ]:
        for unit in UNITS:
            subset = burden.loc[burden["infercnv_unit"].eq(unit)]
            first = subset.iloc[0]
            threshold_rows.append(
                {
                    "run": label,
                    "unit_id": unit,
                    "n_query_groups": int((subset["role"] == "query").sum()),
                    "n_reference_groups": int((subset["role"] == "reference").sum()),
                    "reference_median_cnv_gene_burden": float(first["reference_median_cnv_gene_burden"]),
                    "reference_mad_cnv_gene_burden": float(first["reference_mad_cnv_gene_burden"]),
                    "reference_q95_cnv_gene_burden": float(first["reference_q95_cnv_gene_burden"]),
                    "high_cnv_threshold": float(first["high_cnv_threshold"]),
                }
            )
    thresholds = pd.DataFrame(threshold_rows)

    mixed = pd.read_csv(SUMMARY_V1 / "tumor_mixed_cnv_calls_v1.tsv", sep="\t")
    copy_stability = pd.read_csv(SUMMARY_V1 / "primary_clear_call_stability_v1.tsv", sep="\t")
    copy_stability = copy_stability.loc[copy_stability["method"].eq("copykat")].set_index("unit_id")
    patient_to_copy = {"GSE156632__RCC3": "RCC3t", "GSE156632__RCC5": "RCC5t"}
    control_effect = comparison.loc[
        comparison["left_run"].eq("control_no_add_4threads")
        & comparison["right_run"].eq("plus_one_mixed_4threads")
    ].set_index("unit_id")
    mixed["infercnv_control_to_plus_one_clear_agreement_fraction_v2"] = mixed["curated_patient_id"].map(
        control_effect["exact_call_agreement_fraction"]
    )
    mixed["copykat_primary_to_plus_one_clear_agreement_fraction_v2"] = [
        float(copy_stability.loc[patient_to_copy[patient], "exact_call_agreement_fraction"])
        for patient in mixed["curated_patient_id"]
    ]
    mixed["infercnv_sensitivity_reliability_v2"] = np.where(
        mixed["infercnv_control_to_plus_one_clear_agreement_fraction_v2"].ge(0.95),
        "stable_clear_background",
        "unstable_clear_background",
    )
    mixed["copykat_sensitivity_reliability_v2"] = np.where(
        mixed["copykat_primary_to_plus_one_clear_agreement_fraction_v2"].ge(0.95),
        "stable_clear_background",
        "unstable_clear_background",
    )
    mixed["combined_sensitivity_interpretation_v2"] = np.where(
        mixed["infercnv_sensitivity_reliability_v2"].eq("stable_clear_background"),
        "two_method_support_in_plus_one_run_with_stable_clear_background",
        "copykat_support_stable_but_infercnv_plus_one_call_has_unstable_clear_background",
    )
    mixed["eligible_for_primary_or_subclone_selection_v2"] = False
    mixed["exclusion_reason_v2"] = (
        "mixed lineage; only one tumor cell per patient/sample; sensitivity result is not part of primary malignancy calls"
    )

    comparison.to_csv(OUT / "infercnv_paired_control_call_stability_v2.tsv", sep="\t", index=False, lineterminator="\n")
    transition.to_csv(OUT / "infercnv_paired_control_call_transitions_v2.tsv", sep="\t", index=False, lineterminator="\n")
    thresholds.to_csv(OUT / "infercnv_paired_control_thresholds_v2.tsv", sep="\t", index=False, lineterminator="\n")
    control.to_csv(OUT / "infercnv_control_query_cell_calls_v2.tsv.gz", sep="\t", index=False, compression="gzip", lineterminator="\n")
    control_burden.to_csv(OUT / "infercnv_control_group_burden_v2.tsv", sep="\t", index=False, lineterminator="\n")
    mixed.to_csv(OUT / "tumor_mixed_cnv_calls_with_control_v2.tsv", sep="\t", index=False, lineterminator="\n")
    pd.DataFrame(gates).to_csv(OUT / "memory_gate_log_v2.tsv", sep="\t", index=False, lineterminator="\n")

    fig, ax = plt.subplots(figsize=(10.5, 4.8))
    order = comparison.copy()
    order["label"] = order["unit_id"].str.replace("GSE156632__", "", regex=False) + "\n" + order["left_run"].str.replace("_", " ") + "\nvs " + order["right_run"].str.replace("_", " ")
    colors = np.where(order["left_run"].eq("control_no_add_4threads"), "#d95f02", "#1b9e77")
    ax.bar(np.arange(len(order)), order["exact_call_agreement_fraction"], color=colors)
    ax.set_xticks(np.arange(len(order)), order["label"], rotation=35, ha="right", fontsize=7)
    ax.set_ylim(0, 1.03)
    ax.set_ylabel("Exact clear-cell inferCNV call agreement")
    ax.set_title("Matched no-add controls separate rerun/thread effects from +1 mixed-cell effects")
    for i, value in enumerate(order["exact_call_agreement_fraction"]):
        ax.text(i, min(value + 0.02, 1.0), f"{value:.1%}", ha="center", fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / "infercnv_paired_control_stability_v2.png", dpi=300)
    fig.savefig(OUT / "infercnv_paired_control_stability_v2.pdf")
    plt.close(fig)

    metadata = {
        "schema_version": "part1_mixed_cnv_sensitivity_control_summary_v2",
        "completed_at": datetime.now().astimezone().isoformat(),
        "design": "primary 8-thread vs no-add 4-thread vs plus-one-mixed 4-thread",
        "n_control_units": 2,
        "mixed_included_in_primary": False,
        "mixed_eligible_for_subclone_selection": False,
    }
    (OUT / "run_metadata_v2.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    manifest = []
    for path in sorted(OUT.iterdir()):
        if path.is_file() and path.name not in {"file_manifest_v2.tsv", "RUN_COMPLETE"}:
            manifest.append({"file": path.name, "size_bytes": path.stat().st_size, "sha256": sha256(path)})
    pd.DataFrame(manifest).to_csv(OUT / "file_manifest_v2.tsv", sep="\t", index=False, lineterminator="\n")
    (OUT / "RUN_COMPLETE").write_text(datetime.now().astimezone().isoformat() + "\n", encoding="utf-8")
    print(comparison.to_string(index=False))
    print(mixed[["exact_cell_id", "infercnv_call", "copykat_call", "infercnv_sensitivity_reliability_v2", "combined_sensitivity_interpretation_v2", "eligible_for_primary_or_subclone_selection_v2"]].to_string(index=False))
    print("PART1_MIXED_CNV_CONTROL_SUMMARY_V2_PASS")


if __name__ == "__main__":
    main()
