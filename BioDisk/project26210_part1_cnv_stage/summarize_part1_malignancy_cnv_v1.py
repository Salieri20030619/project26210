#!/usr/bin/env python3
from __future__ import annotations

import csv
import gzip
import hashlib
import json
import subprocess
from datetime import datetime
from pathlib import Path

import anndata as ad
import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt


PROJECT = Path("/mnt/BioDisk/ludingyu/project26210")
ROOT = PROJECT / "results/part1_rebuild_v1/03_malignancy_cnv"
INPUT_H5AD = PROJECT / "results/part1_rebuild_v1/02_integration_annotation/cell_annotation_v1/part1_scvi_integrated_cell_annotated_v1.h5ad"
INFER_ROOT = ROOT / "infercnv_v2/patients"
COPY_ROOT = ROOT / "copykat_v2/samples"
COPY_INPUT_ROOT = ROOT / "cnv_unit_inputs_v1/copykat"
TERMINAL_TABLE = ROOT / "cnv_v2_terminal_recovery_v1/resolved_unit_status_v1.tsv"
OUT = ROOT / "malignancy_cnv_summary_v1"
MEMORY_GATE = PROJECT / "check_memory_gate.sh"
NEUTRAL_STATE = 3


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
        raise RuntimeError(f"Memory gate blocked: {stage}: {output}")


def scaled_mad(values: np.ndarray) -> float:
    median = float(np.median(values))
    return float(np.median(np.abs(values - median)) * 1.4826)


def read_infer_unit(unit: str, gates: list[dict[str, str]]) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, object]]:
    memory_gate(f"summarize_infercnv_{unit}", gates)
    unit_dir = INFER_ROOT / unit
    run_dir = unit_dir / "infercnv_run"
    mapping = pd.read_csv(unit_dir / "cell_id_mapping_v2.tsv", sep="\t", dtype=str)
    groups = pd.read_csv(
        run_dir / "17_HMM_predHMMi6.leiden.hmm_mode-subclusters.cell_groupings",
        sep="\t",
        dtype=str,
    )
    genes = pd.read_csv(
        run_dir / "17_HMM_predHMMi6.leiden.hmm_mode-subclusters.pred_cnv_genes.dat",
        sep="\t",
        dtype={"cell_group_name": str, "state": int, "gene": str, "chr": str},
    )
    groups = groups.rename(columns={"cell": "tool_id"})
    groups = groups.merge(
        mapping[["tool_id", "exact_cell_id", "role", "sample_id", "curated_patient_id", "study"]],
        on="tool_id",
        how="left",
        validate="one_to_one",
    )
    if groups["exact_cell_id"].isna().any() or groups["tool_id"].duplicated().any():
        raise RuntimeError(f"{unit}: incomplete inferCNV cell mapping")
    role_counts = groups.groupby("cell_group_name", observed=True)["role"].nunique()
    if not (role_counts == 1).all():
        raise RuntimeError(f"{unit}: inferCNV group mixes query and reference roles")

    genes["abs_state_shift"] = (genes["state"] - NEUTRAL_STATE).abs()
    burden = genes.groupby("cell_group_name", observed=True).agg(
        n_cnv_genes=("gene", "size"),
        n_gain_genes=("state", lambda x: int((x > NEUTRAL_STATE).sum())),
        n_loss_genes=("state", lambda x: int((x < NEUTRAL_STATE).sum())),
        cnv_gene_burden=("abs_state_shift", "sum"),
        mean_abs_state_shift=("abs_state_shift", "mean"),
        max_abs_state_shift=("abs_state_shift", "max"),
        n_chromosomes=("chr", "nunique"),
    ).reset_index()
    group_roles = groups.groupby("cell_group_name", observed=True).agg(
        role=("role", "first"),
        n_cells=("tool_id", "size"),
    ).reset_index()
    burden = group_roles.merge(burden, on="cell_group_name", how="left", validate="one_to_one")
    numeric = [
        "n_cnv_genes", "n_gain_genes", "n_loss_genes", "cnv_gene_burden",
        "mean_abs_state_shift", "max_abs_state_shift", "n_chromosomes",
    ]
    burden[numeric] = burden[numeric].fillna(0)
    ref_scores = burden.loc[burden["role"].eq("reference"), "cnv_gene_burden"].to_numpy(float)
    if ref_scores.size == 0:
        raise RuntimeError(f"{unit}: no reference HMM groups")
    ref_median = float(np.median(ref_scores))
    ref_mad = scaled_mad(ref_scores)
    ref_q95 = float(np.quantile(ref_scores, 0.95))
    threshold = max(ref_median + 3 * ref_mad, ref_q95)
    burden["infercnv_call"] = np.where(
        burden["role"].eq("reference"),
        "reference",
        np.where(burden["cnv_gene_burden"].gt(threshold), "high_cnv", "low_cnv"),
    )
    burden.insert(0, "infercnv_unit", unit)
    burden["reference_median_cnv_gene_burden"] = ref_median
    burden["reference_mad_cnv_gene_burden"] = ref_mad
    burden["reference_q95_cnv_gene_burden"] = ref_q95
    burden["high_cnv_threshold"] = threshold

    cell_calls = groups.merge(
        burden[["cell_group_name", "cnv_gene_burden", "infercnv_call", "high_cnv_threshold"]],
        on="cell_group_name",
        how="left",
        validate="many_to_one",
    )
    cell_calls.insert(0, "infercnv_unit", unit)
    query = cell_calls.loc[cell_calls["role"].eq("query")].copy()
    if query["infercnv_call"].isna().any() or query["exact_cell_id"].duplicated().any():
        raise RuntimeError(f"{unit}: invalid inferCNV query calls")
    metadata = {
        "infercnv_unit": unit,
        "n_query_cells": int(len(query)),
        "n_reference_cells": int((cell_calls["role"] == "reference").sum()),
        "n_query_groups": int((burden["role"] == "query").sum()),
        "n_reference_groups": int((burden["role"] == "reference").sum()),
        "reference_median": ref_median,
        "reference_mad": ref_mad,
        "reference_q95": ref_q95,
        "high_cnv_threshold": threshold,
    }
    return query, burden, metadata


def read_copy_unit(unit: str, resolved_status: str, gates: list[dict[str, str]]) -> pd.DataFrame:
    memory_gate(f"summarize_copykat_{unit}", gates)
    if resolved_status == "completed":
        calls = pd.read_csv(COPY_ROOT / unit / "copykat_predictions_mapped_v2.tsv", sep="\t", dtype=str)
        calls = calls.loc[calls["role"].eq("query")].copy()
        raw = calls["copykat.pred"].fillna("")
        allowed = {"aneuploid", "diploid", "not.defined"}
        if not set(raw.unique()).issubset(allowed):
            raise RuntimeError(f"{unit}: unexpected CopyKAT calls: {sorted(set(raw.unique()))}")
        calls["copykat_call"] = raw.map({"aneuploid": "aneuploid", "diploid": "diploid", "not.defined": "not_defined"})
        keep = ["exact_cell_id", "sample_id", "curated_patient_id", "study", "copykat.pred", "copykat_call"]
        calls = calls[keep]
    elif resolved_status == "not_estimable":
        with gzip.open(COPY_INPUT_ROOT / unit / "cell_id_mapping_v1.tsv.gz", "rt") as handle:
            calls = pd.read_csv(handle, sep="\t", dtype=str)
        calls = calls.loc[calls["role"].eq("query")].copy()
        calls["copykat.pred"] = pd.NA
        calls["copykat_call"] = "not_estimable"
        calls = calls[["exact_cell_id", "sample_id", "curated_patient_id", "study", "copykat.pred", "copykat_call"]]
    else:
        raise RuntimeError(f"{unit}: unsupported resolved status {resolved_status}")
    calls.insert(0, "copykat_unit", unit)
    if calls["exact_cell_id"].isna().any() or calls["exact_cell_id"].duplicated().any():
        raise RuntimeError(f"{unit}: invalid CopyKAT query mapping")
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
    if (infer_call == "high_cnv" and copy_call == "diploid") or (infer_call == "low_cnv" and copy_call == "aneuploid"):
        return "method_discordant", False, "infercnv_copykat_direction_conflict"
    return "indeterminate", False, "insufficient_or_single_low_cnv_evidence"


def plot_cross_table(table: pd.DataFrame, out_prefix: Path) -> None:
    infer_order = ["high_cnv", "low_cnv", "not_eligible"]
    copy_order = ["aneuploid", "diploid", "not_defined", "not_estimable", "not_eligible"]
    matrix = table.pivot(index="infercnv_call", columns="copykat_call", values="n_cells").reindex(index=infer_order, columns=copy_order, fill_value=0).fillna(0)
    values = matrix.to_numpy(float)
    fig, ax = plt.subplots(figsize=(8.2, 4.8))
    image = ax.imshow(np.log10(values + 1), cmap="Blues", aspect="auto")
    for i in range(values.shape[0]):
        for j in range(values.shape[1]):
            ax.text(j, i, f"{int(values[i, j]):,}", ha="center", va="center", fontsize=9)
    ax.set_xticks(range(len(copy_order)), [x.replace("_", " ") for x in copy_order], rotation=30, ha="right")
    ax.set_yticks(range(len(infer_order)), [x.replace("_", " ") for x in infer_order])
    ax.set_xlabel("CopyKAT call")
    ax.set_ylabel("inferCNV call")
    ax.set_title("Tumor epithelial CNV method cross-table")
    fig.colorbar(image, ax=ax, label="log10(n cells + 1)")
    fig.tight_layout()
    fig.savefig(out_prefix.with_suffix(".png"), dpi=300)
    fig.savefig(out_prefix.with_suffix(".pdf"))
    plt.close(fig)


def plot_patient_status(patient: pd.DataFrame, out_prefix: Path) -> None:
    order = [
        "malignant_consensus", "malignant_infercnv_only", "malignant_copykat_only",
        "non_malignant_consensus", "method_discordant", "indeterminate",
    ]
    colors = {
        "malignant_consensus": "#b2182b",
        "malignant_infercnv_only": "#ef8a62",
        "malignant_copykat_only": "#fddbc7",
        "non_malignant_consensus": "#2166ac",
        "method_discordant": "#762a83",
        "indeterminate": "#bdbdbd",
    }
    wide = patient.pivot(index="curated_patient_id", columns="malignancy_status_v1", values="fraction").reindex(columns=order, fill_value=0).fillna(0)
    wide = wide.assign(_malignant=wide[[x for x in order if x.startswith("malignant_")]].sum(axis=1)).sort_values("_malignant", ascending=False).drop(columns="_malignant")
    fig, ax = plt.subplots(figsize=(14, 7.5))
    bottom = np.zeros(len(wide))
    x = np.arange(len(wide))
    for status in order:
        values = wide[status].to_numpy(float)
        ax.bar(x, values, bottom=bottom, color=colors[status], width=0.85, label=status.replace("_", " "))
        bottom += values
    ax.set_xticks(x, wide.index, rotation=70, ha="right", fontsize=7)
    ax.set_ylim(0, 1)
    ax.set_ylabel("Fraction of tumor epithelial cells")
    ax.set_xlabel("Curated patient")
    ax.set_title("Conservative malignancy evidence by patient")
    ax.legend(frameon=False, ncol=3, loc="upper center", bbox_to_anchor=(0.5, 1.16), fontsize=8)
    fig.tight_layout()
    fig.savefig(out_prefix.with_suffix(".png"), dpi=300)
    fig.savefig(out_prefix.with_suffix(".pdf"))
    plt.close(fig)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    gates: list[dict[str, str]] = []
    memory_gate("read_annotated_h5ad_obs", gates)
    adata = ad.read_h5ad(INPUT_H5AD, backed="r")
    required = ["study", "sample_id", "curated_patient_id", "source_tissue", "cell_type_broad_v1"]
    missing = [column for column in required if column not in adata.obs.columns]
    if missing:
        raise RuntimeError(f"Annotated H5AD missing columns: {missing}")
    epithelial = adata.obs.loc[adata.obs["cell_type_broad_v1"].astype(str).eq("Epithelial cell"), required].copy()
    epithelial.insert(0, "exact_cell_id", epithelial.index.astype(str))
    adata.file.close()
    if len(epithelial) != 62652 or epithelial["exact_cell_id"].duplicated().any():
        raise RuntimeError(f"Unexpected epithelial registry: n={len(epithelial)}")
    tumor = epithelial.loc[epithelial["source_tissue"].astype(str).eq("Tumor")].copy()
    normal = epithelial.loc[epithelial["source_tissue"].astype(str).eq("Normal")].copy()
    if len(tumor) != 38670 or len(normal) != 23982:
        raise RuntimeError(f"Unexpected epithelial tissue counts: tumor={len(tumor)} normal={len(normal)}")

    terminal = pd.read_csv(TERMINAL_TABLE, sep="\t", dtype=str)
    infer_units = terminal.loc[(terminal["method"] == "infercnv") & (terminal["resolved_status"] == "completed"), "unit_id"].tolist()
    copy_terminal = terminal.loc[terminal["method"] == "copykat", ["unit_id", "resolved_status"]]
    if len(infer_units) != 34 or len(copy_terminal) != 41:
        raise RuntimeError("Terminal table unit counts do not match 34 inferCNV / 41 CopyKAT")

    infer_cells = []
    infer_groups = []
    infer_params = []
    for unit in sorted(infer_units):
        cells, groups, params = read_infer_unit(unit, gates)
        infer_cells.append(cells)
        infer_groups.append(groups)
        infer_params.append(params)
    infer_cells_df = pd.concat(infer_cells, ignore_index=True)
    infer_groups_df = pd.concat(infer_groups, ignore_index=True)
    if infer_cells_df["exact_cell_id"].duplicated().any():
        raise RuntimeError("inferCNV query cells are duplicated across patient units")

    copy_cells = []
    for row in copy_terminal.sort_values("unit_id").itertuples(index=False):
        copy_cells.append(read_copy_unit(row.unit_id, row.resolved_status, gates))
    copy_cells_df = pd.concat(copy_cells, ignore_index=True)
    if copy_cells_df["exact_cell_id"].duplicated().any():
        raise RuntimeError("CopyKAT query cells are duplicated across sample units")

    tumor_ids = set(tumor["exact_cell_id"])
    if not set(infer_cells_df["exact_cell_id"]).issubset(tumor_ids):
        raise RuntimeError("inferCNV query calls include non-tumor/non-epithelial cells")
    if not set(copy_cells_df["exact_cell_id"]).issubset(tumor_ids):
        raise RuntimeError("CopyKAT query calls include non-tumor/non-epithelial cells")

    infer_map = infer_cells_df.set_index("exact_cell_id")
    copy_map = copy_cells_df.set_index("exact_cell_id")
    tumor["infercnv_call"] = tumor["exact_cell_id"].map(infer_map["infercnv_call"]).fillna("not_eligible")
    tumor["infercnv_unit"] = tumor["exact_cell_id"].map(infer_map["infercnv_unit"])
    tumor["infercnv_group"] = tumor["exact_cell_id"].map(infer_map["cell_group_name"])
    tumor["infercnv_cnv_gene_burden"] = pd.to_numeric(tumor["exact_cell_id"].map(infer_map["cnv_gene_burden"]), errors="coerce")
    tumor["infercnv_high_cnv_threshold"] = pd.to_numeric(tumor["exact_cell_id"].map(infer_map["high_cnv_threshold"]), errors="coerce")
    tumor["copykat_call"] = tumor["exact_cell_id"].map(copy_map["copykat_call"]).fillna("not_eligible")
    tumor["copykat_raw_call"] = tumor["exact_cell_id"].map(copy_map["copykat.pred"])
    tumor["copykat_unit"] = tumor["exact_cell_id"].map(copy_map["copykat_unit"])

    integrated = [integrated_status(i, c) for i, c in zip(tumor["infercnv_call"], tumor["copykat_call"])]
    tumor["malignancy_status_v1"] = [x[0] for x in integrated]
    tumor["malignant_conservative_v1"] = [x[1] for x in integrated]
    tumor["malignancy_evidence_v1"] = [x[2] for x in integrated]
    if tumor["malignancy_status_v1"].isna().any():
        raise RuntimeError("Integrated malignancy status contains NA")

    cross = tumor.groupby(["infercnv_call", "copykat_call"], observed=True).size().rename("n_cells").reset_index()
    status_counts = tumor.groupby(["malignancy_status_v1", "malignant_conservative_v1"], observed=True).size().rename("n_cells").reset_index()
    patient_status = tumor.groupby(["curated_patient_id", "malignancy_status_v1"], observed=True).size().rename("n_cells").reset_index()
    patient_status["fraction"] = patient_status["n_cells"] / patient_status.groupby("curated_patient_id", observed=True)["n_cells"].transform("sum")
    patient_summary = tumor.groupby("curated_patient_id", observed=True).agg(
        n_tumor_epithelial=("exact_cell_id", "size"),
        n_infercnv_high=("infercnv_call", lambda x: int((x == "high_cnv").sum())),
        n_infercnv_low=("infercnv_call", lambda x: int((x == "low_cnv").sum())),
        n_infercnv_not_eligible=("infercnv_call", lambda x: int((x == "not_eligible").sum())),
        n_copykat_aneuploid=("copykat_call", lambda x: int((x == "aneuploid").sum())),
        n_copykat_diploid=("copykat_call", lambda x: int((x == "diploid").sum())),
        n_copykat_not_defined=("copykat_call", lambda x: int((x == "not_defined").sum())),
        n_copykat_not_estimable=("copykat_call", lambda x: int((x == "not_estimable").sum())),
        n_copykat_not_eligible=("copykat_call", lambda x: int((x == "not_eligible").sum())),
        n_malignant_conservative=("malignant_conservative_v1", "sum"),
        n_method_discordant=("malignancy_status_v1", lambda x: int((x == "method_discordant").sum())),
    ).reset_index()
    patient_summary["fraction_malignant_conservative"] = patient_summary["n_malignant_conservative"] / patient_summary["n_tumor_epithelial"]
    sample_summary = tumor.groupby(["curated_patient_id", "sample_id"], observed=True).agg(
        n_tumor_epithelial=("exact_cell_id", "size"),
        n_malignant_conservative=("malignant_conservative_v1", "sum"),
        n_copykat_aneuploid=("copykat_call", lambda x: int((x == "aneuploid").sum())),
        n_copykat_diploid=("copykat_call", lambda x: int((x == "diploid").sum())),
        n_copykat_unavailable=("copykat_call", lambda x: int(x.isin(["not_defined", "not_estimable", "not_eligible"]).sum())),
    ).reset_index()
    sample_summary["fraction_malignant_conservative"] = sample_summary["n_malignant_conservative"] / sample_summary["n_tumor_epithelial"]

    definitive = tumor.loc[tumor["infercnv_call"].isin(["high_cnv", "low_cnv"]) & tumor["copykat_call"].isin(["aneuploid", "diploid"])]
    agreement = int(((definitive["infercnv_call"] == "high_cnv") == (definitive["copykat_call"] == "aneuploid")).sum())
    coverage = pd.DataFrame([
        {"metric": "all_epithelial_cells", "value": len(epithelial)},
        {"metric": "tumor_epithelial_cells", "value": len(tumor)},
        {"metric": "normal_epithelial_cells", "value": len(normal)},
        {"metric": "infercnv_eligible_query_cells", "value": int((tumor["infercnv_call"] != "not_eligible").sum())},
        {"metric": "copykat_completed_or_not_estimable_query_cells", "value": int((tumor["copykat_call"] != "not_eligible").sum())},
        {"metric": "both_methods_definitive_cells", "value": len(definitive)},
        {"metric": "both_methods_definitive_agreement_cells", "value": agreement},
        {"metric": "both_methods_definitive_agreement_fraction", "value": agreement / len(definitive) if len(definitive) else np.nan},
        {"metric": "malignant_conservative_cells", "value": int(tumor["malignant_conservative_v1"].sum())},
        {"metric": "patients_with_malignant_conservative_cells", "value": int((patient_summary["n_malignant_conservative"] > 0).sum())},
    ])

    memory_gate("write_malignancy_summary_outputs", gates)
    infer_cells_df.to_csv(OUT / "infercnv_query_cell_calls_v1.tsv.gz", sep="\t", index=False, compression="gzip")
    infer_groups_df.to_csv(OUT / "infercnv_group_burden_v1.tsv", sep="\t", index=False)
    pd.DataFrame(infer_params).to_csv(OUT / "infercnv_threshold_parameters_v1.tsv", sep="\t", index=False)
    copy_cells_df.to_csv(OUT / "copykat_query_cell_calls_v1.tsv.gz", sep="\t", index=False, compression="gzip")
    tumor.to_csv(OUT / "tumor_epithelial_malignancy_calls_v1.tsv.gz", sep="\t", index=False, compression="gzip")
    cross.to_csv(OUT / "method_cross_table_v1.tsv", sep="\t", index=False)
    status_counts.to_csv(OUT / "malignancy_status_counts_v1.tsv", sep="\t", index=False)
    patient_status.to_csv(OUT / "patient_status_composition_v1.tsv", sep="\t", index=False)
    patient_summary.to_csv(OUT / "patient_malignancy_summary_v1.tsv", sep="\t", index=False)
    sample_summary.to_csv(OUT / "sample_malignancy_summary_v1.tsv", sep="\t", index=False)
    coverage.to_csv(OUT / "method_coverage_summary_v1.tsv", sep="\t", index=False)
    pd.DataFrame(gates).to_csv(OUT / "memory_gate_log_v1.tsv", sep="\t", index=False)
    plot_cross_table(cross, OUT / "method_cross_table_v1")
    plot_patient_status(patient_status, OUT / "patient_malignancy_composition_v1")

    output_files = sorted(path for path in OUT.iterdir() if path.is_file() and path.name not in {"file_manifest_v1.tsv", "RUN_COMPLETE", "run_metadata_v1.json"})
    with (OUT / "file_manifest_v1.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["path", "size_bytes", "sha256"], delimiter="\t")
        writer.writeheader()
        for path in output_files:
            writer.writerow({"path": str(path), "size_bytes": path.stat().st_size, "sha256": sha256(path)})
    status_dict = dict(zip(status_counts["malignancy_status_v1"], status_counts["n_cells"].astype(int)))
    metadata = {
        "schema_version": "part1_malignancy_cnv_summary_v1",
        "created_at": datetime.now().astimezone().isoformat(),
        "input_h5ad": str(INPUT_H5AD),
        "input_h5ad_sha256": sha256(INPUT_H5AD),
        "n_epithelial": len(epithelial),
        "n_tumor_epithelial": len(tumor),
        "n_normal_epithelial": len(normal),
        "infercnv_units": len(infer_units),
        "copykat_units_completed": int((copy_terminal["resolved_status"] == "completed").sum()),
        "copykat_units_not_estimable": int((copy_terminal["resolved_status"] == "not_estimable").sum()),
        "infercnv_rule": "per-unit query group burden > max(reference median + 3*scaled MAD, reference q95)",
        "malignant_conservative_rule": "consensus high/aneuploid OR inferCNV high with CopyKAT unavailable OR CopyKAT aneuploid with inferCNV unavailable",
        "malignancy_status_counts": status_dict,
        "both_methods_definitive_n": len(definitive),
        "both_methods_definitive_agreement_n": agreement,
        "both_methods_definitive_agreement_fraction": agreement / len(definitive) if len(definitive) else None,
        "interpretation": "cell counts and agreement are descriptive; cells are not independent biological replicates",
    }
    (OUT / "run_metadata_v1.json").write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    (OUT / "RUN_COMPLETE").write_text(datetime.now().astimezone().isoformat() + "\n")
    print(
        "PART1_MALIGNANCY_CNV_SUMMARY_V1_ALL_PASS "
        f"tumor_epithelial={len(tumor)} malignant_conservative={int(tumor['malignant_conservative_v1'].sum())} "
        f"definitive_overlap={len(definitive)} agreement={agreement}",
        flush=True,
    )


if __name__ == "__main__":
    main()
