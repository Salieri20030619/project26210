#!/usr/bin/env python3
"""Summarize the isolated mixed-lineage CNV sensitivity reruns."""

from __future__ import annotations

import hashlib
import json
import subprocess
from datetime import datetime
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt


PROJECT = Path("/mnt/BioDisk/ludingyu/project26210")
ROOT = PROJECT / "results/part1_rebuild_v1/03_malignancy_cnv"
SENS = ROOT / "mixed_cnv_sensitivity_v1"
OUT = SENS / "summary_v1"
PRIMARY = ROOT / "malignancy_cnv_summary_v1"
MEMORY_GATE = PROJECT / "check_memory_gate.sh"
NEUTRAL_STATE = 3
UNITS = {
    "infercnv": ["GSE156632__RCC3", "GSE156632__RCC5"],
    "copykat": ["RCC3t", "RCC5t"],
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def memory_gate(stage: str, records: list[dict[str, str]]) -> None:
    result = subprocess.run([str(MEMORY_GATE)], text=True, capture_output=True)
    output = " ".join((result.stdout + " " + result.stderr).split())
    records.append({"stage": stage, "exit_status": str(result.returncode), "gate_output": output})
    print(output, flush=True)
    if result.returncode != 0 or "PASS" not in output:
        raise RuntimeError(f"Memory gate blocked {stage}: {output}")


def scaled_mad(values: np.ndarray) -> float:
    median = float(np.median(values))
    return float(np.median(np.abs(values - median)) * 1.4826)


def read_infer(unit: str, gates: list[dict[str, str]]) -> tuple[pd.DataFrame, pd.DataFrame]:
    memory_gate(f"summarize_mixed_infercnv_{unit}", gates)
    unit_dir = SENS / "infercnv" / unit
    mapping = pd.read_csv(unit_dir / "cell_id_mapping_v2.tsv", sep="\t", dtype=str)
    groups = pd.read_csv(
        unit_dir / "infercnv_run/17_HMM_predHMMi6.leiden.hmm_mode-subclusters.cell_groupings",
        sep="\t", dtype=str,
    ).rename(columns={"cell": "tool_id"})
    genes = pd.read_csv(
        unit_dir / "infercnv_run/17_HMM_predHMMi6.leiden.hmm_mode-subclusters.pred_cnv_genes.dat",
        sep="\t", dtype={"cell_group_name": str, "state": int, "gene": str, "chr": str},
    )
    groups = groups.merge(
        mapping[
            [
                "tool_id", "exact_cell_id", "role", "sample_id", "curated_patient_id",
                "study", "source_tissue", "cell_type_broad_v1",
            ]
        ],
        on="tool_id", how="left", validate="one_to_one",
    )
    if groups["exact_cell_id"].isna().any() or groups["tool_id"].duplicated().any():
        raise RuntimeError(f"{unit}: incomplete inferCNV mapping")
    if not (groups.groupby("cell_group_name", observed=True)["role"].nunique() == 1).all():
        raise RuntimeError(f"{unit}: HMM group mixes query and reference")
    genes["abs_state_shift"] = (genes["state"] - NEUTRAL_STATE).abs()
    burden = genes.groupby("cell_group_name", observed=True).agg(
        n_cnv_genes=("gene", "size"),
        n_gain_genes=("state", lambda x: int((x > NEUTRAL_STATE).sum())),
        n_loss_genes=("state", lambda x: int((x < NEUTRAL_STATE).sum())),
        cnv_gene_burden=("abs_state_shift", "sum"),
        n_chromosomes=("chr", "nunique"),
    ).reset_index()
    group_roles = groups.groupby("cell_group_name", observed=True).agg(
        role=("role", "first"), n_cells=("tool_id", "size")
    ).reset_index()
    burden = group_roles.merge(burden, on="cell_group_name", how="left", validate="one_to_one")
    numeric = ["n_cnv_genes", "n_gain_genes", "n_loss_genes", "cnv_gene_burden", "n_chromosomes"]
    burden[numeric] = burden[numeric].fillna(0)
    ref = burden.loc[burden["role"].eq("reference"), "cnv_gene_burden"].to_numpy(float)
    if ref.size == 0:
        raise RuntimeError(f"{unit}: no reference HMM groups")
    ref_median = float(np.median(ref))
    ref_mad = scaled_mad(ref)
    ref_q95 = float(np.quantile(ref, 0.95))
    threshold = max(ref_median + 3 * ref_mad, ref_q95)
    burden["infercnv_call"] = np.where(
        burden["role"].eq("reference"), "reference",
        np.where(burden["cnv_gene_burden"].gt(threshold), "high_cnv", "low_cnv"),
    )
    burden.insert(0, "infercnv_unit", unit)
    burden["reference_median_cnv_gene_burden"] = ref_median
    burden["reference_mad_cnv_gene_burden"] = ref_mad
    burden["reference_q95_cnv_gene_burden"] = ref_q95
    burden["high_cnv_threshold"] = threshold
    cells = groups.merge(
        burden[["cell_group_name", "cnv_gene_burden", "infercnv_call", "high_cnv_threshold"]],
        on="cell_group_name", how="left", validate="many_to_one",
    )
    cells.insert(0, "infercnv_unit", unit)
    query = cells.loc[cells["role"].eq("query")].copy()
    if query["infercnv_call"].isna().any() or query["exact_cell_id"].duplicated().any():
        raise RuntimeError(f"{unit}: invalid inferCNV query calls")
    return query, burden


def read_copykat(unit: str, gates: list[dict[str, str]]) -> pd.DataFrame:
    memory_gate(f"summarize_mixed_copykat_{unit}", gates)
    calls = pd.read_csv(
        SENS / "copykat" / unit / "copykat_predictions_mapped_v2.tsv", sep="\t", dtype=str
    )
    calls = calls.loc[calls["role"].eq("query")].copy()
    allowed = {"aneuploid", "diploid", "not.defined"}
    if not set(calls["copykat.pred"].dropna().unique()).issubset(allowed):
        raise RuntimeError(f"{unit}: unexpected CopyKAT call")
    calls["copykat_call"] = calls["copykat.pred"].map(
        {"aneuploid": "aneuploid", "diploid": "diploid", "not.defined": "not_defined"}
    )
    calls.insert(0, "copykat_unit", unit)
    if calls["exact_cell_id"].isna().any() or calls["exact_cell_id"].duplicated().any():
        raise RuntimeError(f"{unit}: invalid CopyKAT mapping")
    return calls


def integrated_status(infer_call: str, copy_call: str) -> tuple[str, bool, str]:
    if infer_call == "high_cnv" and copy_call == "aneuploid":
        return "malignant_consensus", True, "infercnv_high_and_copykat_aneuploid"
    if infer_call == "high_cnv" and copy_call in {"not_defined", "not_estimable", "not_eligible"}:
        return "malignant_infercnv_only", True, "infercnv_high_copykat_unavailable"
    if infer_call == "not_eligible" and copy_call == "aneuploid":
        return "malignant_copykat_only", True, "copykat_aneuploid_infercnv_unavailable"
    if infer_call == "low_cnv" and copy_call == "diploid":
        return "non_malignant_consensus", False, "infercnv_low_and_copykat_diploid"
    if (infer_call == "high_cnv" and copy_call == "diploid") or (
        infer_call == "low_cnv" and copy_call == "aneuploid"
    ):
        return "method_discordant", False, "infercnv_copykat_direction_conflict"
    return "indeterminate", False, "insufficient_or_single_low_cnv_evidence"


def stability_table(
    method: str,
    units: list[str],
    primary: pd.DataFrame,
    sensitivity: pd.DataFrame,
    call_column: str,
    unit_column: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    summaries: list[dict[str, object]] = []
    transitions: list[pd.DataFrame] = []
    for unit in units:
        p = primary.loc[primary[unit_column].eq(unit), ["exact_cell_id", call_column]].copy()
        s = sensitivity.loc[
            sensitivity[unit_column].eq(unit) & sensitivity["cell_type_broad_v1"].eq("Epithelial cell"),
            ["exact_cell_id", call_column],
        ].copy()
        merged = p.merge(s, on="exact_cell_id", suffixes=("_primary", "_sensitivity"), validate="one_to_one")
        if len(merged) != len(p) or len(s) != len(p):
            raise RuntimeError(f"{method} {unit}: clear query ID set changed")
        same = merged[f"{call_column}_primary"].eq(merged[f"{call_column}_sensitivity"])
        summaries.append(
            {
                "method": method,
                "unit_id": unit,
                "n_primary_clear_query": len(p),
                "n_sensitivity_clear_query": len(s),
                "n_exact_call_agreement": int(same.sum()),
                "n_call_changed": int((~same).sum()),
                "exact_call_agreement_fraction": float(same.mean()),
            }
        )
        transition = merged.groupby(
            [f"{call_column}_primary", f"{call_column}_sensitivity"], observed=True
        ).size().rename("n_cells").reset_index()
        transition.insert(0, "unit_id", unit)
        transition.insert(0, "method", method)
        transitions.append(transition.rename(columns={
            f"{call_column}_primary": "primary_call",
            f"{call_column}_sensitivity": "sensitivity_call",
        }))
    return pd.DataFrame(summaries), pd.concat(transitions, ignore_index=True)


def plot_summary(mixed: pd.DataFrame, stability: pd.DataFrame, output: Path) -> None:
    call_score = {
        "low_cnv": 0, "diploid": 0, "not_defined": 1,
        "high_cnv": 2, "aneuploid": 2,
    }
    matrix = np.array(
        [
            [call_score.get(x, 1) for x in mixed["infercnv_call"]],
            [call_score.get(x, 1) for x in mixed["copykat_call"]],
        ], dtype=float,
    )
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 3.8), gridspec_kw={"width_ratios": [1.25, 1]})
    axes[0].imshow(matrix, vmin=0, vmax=2, cmap=matplotlib.colors.ListedColormap(["#2166ac", "#bdbdbd", "#b2182b"]), aspect="auto")
    for i, column in enumerate(["infercnv_call", "copykat_call"]):
        for j, value in enumerate(mixed[column]):
            axes[0].text(j, i, str(value).replace("_", "\n"), ha="center", va="center", fontsize=8, color="white" if matrix[i, j] != 1 else "black")
    axes[0].set_yticks([0, 1], ["inferCNV", "CopyKAT"])
    axes[0].set_xticks(range(len(mixed)), mixed["sample_id"], rotation=20, ha="right")
    axes[0].set_title("Tumor mixed-cell sensitivity calls")
    labels = stability["method"] + "\n" + stability["unit_id"].str.replace("GSE156632__", "", regex=False)
    values = stability["exact_call_agreement_fraction"].to_numpy(float)
    axes[1].bar(np.arange(len(values)), values, color=["#4c78a8", "#4c78a8", "#f58518", "#f58518"])
    axes[1].set_xticks(np.arange(len(values)), labels, rotation=30, ha="right", fontsize=8)
    axes[1].set_ylim(0, 1.02)
    axes[1].set_ylabel("Exact clear-cell call agreement")
    axes[1].set_title("Primary vs +1 mixed-cell rerun")
    for i, value in enumerate(values):
        axes[1].text(i, min(value + 0.02, 1.0), f"{value:.1%}", ha="center", va="bottom", fontsize=8)
    fig.tight_layout()
    fig.savefig(output.with_suffix(".png"), dpi=300)
    fig.savefig(output.with_suffix(".pdf"))
    plt.close(fig)


def main() -> None:
    if not (SENS / "RUN_COMPLETE").is_file():
        raise RuntimeError("Sensitivity reruns are not complete")
    OUT.mkdir(parents=True, exist_ok=True)
    gates: list[dict[str, str]] = []
    infer_parts: list[pd.DataFrame] = []
    burden_parts: list[pd.DataFrame] = []
    for unit in UNITS["infercnv"]:
        cells, burden = read_infer(unit, gates)
        infer_parts.append(cells)
        burden_parts.append(burden)
    infer = pd.concat(infer_parts, ignore_index=True)
    burden = pd.concat(burden_parts, ignore_index=True)
    copykat = pd.concat([read_copykat(unit, gates) for unit in UNITS["copykat"]], ignore_index=True)

    primary_infer = pd.read_csv(PRIMARY / "infercnv_query_cell_calls_v1.tsv.gz", sep="\t", dtype=str)
    primary_copy = pd.read_csv(PRIMARY / "copykat_query_cell_calls_v1.tsv.gz", sep="\t", dtype=str)
    infer_stability, infer_transition = stability_table(
        "infercnv", UNITS["infercnv"], primary_infer, infer, "infercnv_call", "infercnv_unit"
    )
    copy_stability, copy_transition = stability_table(
        "copykat", UNITS["copykat"], primary_copy, copykat, "copykat_call", "copykat_unit"
    )
    stability = pd.concat([infer_stability, copy_stability], ignore_index=True)
    transitions = pd.concat([infer_transition, copy_transition], ignore_index=True)

    mixed_infer = infer.loc[infer["cell_type_broad_v1"].eq("Epithelial/Myeloid mixed")].copy()
    mixed_copy = copykat.loc[copykat["cell_type_broad_v1"].eq("Epithelial/Myeloid mixed")].copy()
    if len(mixed_infer) != 2 or len(mixed_copy) != 2:
        raise RuntimeError(f"Unexpected mixed calls: infer={len(mixed_infer)} copykat={len(mixed_copy)}")
    keep_infer = [
        "exact_cell_id", "sample_id", "curated_patient_id", "study", "source_tissue",
        "cell_type_broad_v1", "infercnv_unit", "cell_group_name", "cnv_gene_burden",
        "high_cnv_threshold", "infercnv_call",
    ]
    keep_copy = ["exact_cell_id", "copykat_unit", "copykat.pred", "copykat_call"]
    mixed = mixed_infer[keep_infer].merge(
        mixed_copy[keep_copy], on="exact_cell_id", how="outer", validate="one_to_one"
    )
    if mixed[["infercnv_call", "copykat_call"]].isna().any().any():
        raise RuntimeError("Mixed cell lacks a method call")
    integrated = [integrated_status(i, c) for i, c in zip(mixed["infercnv_call"], mixed["copykat_call"])]
    mixed["malignancy_status_mixed_sensitivity_v1"] = [x[0] for x in integrated]
    mixed["malignant_conservative_mixed_sensitivity_v1"] = [x[1] for x in integrated]
    mixed["malignancy_evidence_mixed_sensitivity_v1"] = [x[2] for x in integrated]
    mixed["included_in_primary_malignancy_summary_v1"] = False

    registry = pd.read_csv(SENS / "inputs/mixed_cell_registry_v1.tsv", sep="\t", dtype=str)
    scope = registry.groupby(["source_tissue", "sensitivity_scope_v1"], observed=True).size().rename("n_cells").reset_index()
    if scope["n_cells"].sum() != 918:
        raise RuntimeError("Mixed registry total changed")
    normal = registry.loc[registry["source_tissue"].eq("Normal")].copy()
    normal["cnv_sensitivity_status_v1"] = "not_applicable_normal_tissue"

    infer.to_csv(OUT / "infercnv_query_cell_calls_mixed_sensitivity_v1.tsv.gz", sep="\t", index=False, compression="gzip", lineterminator="\n")
    burden.to_csv(OUT / "infercnv_group_burden_mixed_sensitivity_v1.tsv", sep="\t", index=False, lineterminator="\n")
    copykat.to_csv(OUT / "copykat_query_cell_calls_mixed_sensitivity_v1.tsv.gz", sep="\t", index=False, compression="gzip", lineterminator="\n")
    mixed.to_csv(OUT / "tumor_mixed_cnv_calls_v1.tsv", sep="\t", index=False, lineterminator="\n")
    normal.to_csv(OUT / "normal_mixed_scope_v1.tsv.gz", sep="\t", index=False, compression="gzip", lineterminator="\n")
    scope.to_csv(OUT / "mixed_scope_summary_v1.tsv", sep="\t", index=False, lineterminator="\n")
    stability.to_csv(OUT / "primary_clear_call_stability_v1.tsv", sep="\t", index=False, lineterminator="\n")
    transitions.to_csv(OUT / "primary_clear_call_transitions_v1.tsv", sep="\t", index=False, lineterminator="\n")
    pd.DataFrame(gates).to_csv(OUT / "memory_gate_log_v1.tsv", sep="\t", index=False, lineterminator="\n")
    plot_summary(mixed, stability, OUT / "mixed_cnv_sensitivity_summary_v1")

    metadata = {
        "schema_version": "part1_mixed_cnv_sensitivity_summary_v1",
        "completed_at": datetime.now().astimezone().isoformat(),
        "n_mixed_total": 918,
        "n_tumor_mixed_rerun": 2,
        "n_normal_mixed_not_tumor_queries": 916,
        "n_infercnv_units": 2,
        "n_copykat_units": 2,
        "primary_main_summary_unchanged": True,
        "mixed_calls_included_in_primary": False,
        "stability_min_exact_call_agreement_fraction": float(stability["exact_call_agreement_fraction"].min()),
    }
    (OUT / "run_metadata_v1.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    manifest_rows = []
    for path in sorted(OUT.iterdir()):
        if path.is_file() and path.name not in {"file_manifest_v1.tsv", "RUN_COMPLETE"}:
            manifest_rows.append({"file": path.name, "size_bytes": path.stat().st_size, "sha256": sha256(path)})
    pd.DataFrame(manifest_rows).to_csv(OUT / "file_manifest_v1.tsv", sep="\t", index=False, lineterminator="\n")
    (OUT / "RUN_COMPLETE").write_text(datetime.now().astimezone().isoformat() + "\n", encoding="utf-8")
    print(mixed[["exact_cell_id", "infercnv_call", "copykat_call", "malignancy_status_mixed_sensitivity_v1"]].to_string(index=False))
    print(stability.to_string(index=False))
    print("PART1_MIXED_CNV_SENSITIVITY_SUMMARY_V1_PASS")


if __name__ == "__main__":
    main()
