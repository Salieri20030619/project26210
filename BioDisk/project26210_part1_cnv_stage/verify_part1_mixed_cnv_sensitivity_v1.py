#!/usr/bin/env python3
"""Independent terminal verification for the isolated mixed-cell CNV sensitivity analysis."""

from __future__ import annotations

import gzip
import hashlib
import json
import subprocess
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
from PIL import Image
from scipy.io import mmread


PROJECT = Path("/mnt/BioDisk/ludingyu/project26210")
ROOT = PROJECT / "results/part1_rebuild_v1/03_malignancy_cnv"
SENS = ROOT / "mixed_cnv_sensitivity_v1"
SUMMARY = SENS / "summary_v1"
CONTROL_SUMMARY = SENS / "summary_v2_with_control"
VERIFY = ROOT / "mixed_cnv_sensitivity_verification_v1"
PRIMARY_INPUT = ROOT / "cnv_unit_inputs_v1"
PRIMARY_SUMMARY = ROOT / "malignancy_cnv_summary_v1"
MEMORY_GATE = PROJECT / "check_memory_gate.sh"
UNITS = [
    ("infercnv", "GSE156632__RCC3", 1508, 200, "RCC3t_GAGGTGATCAGGCGAA.1"),
    ("infercnv", "GSE156632__RCC5", 2381, 200, "RCC5t_ACTTGTTGTTCAGCGC.1"),
    ("copykat", "RCC3t", 1508, 1000, "RCC3t_GAGGTGATCAGGCGAA.1"),
    ("copykat", "RCC5t", 2381, 1000, "RCC5t_ACTTGTTGTTCAGCGC.1"),
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_matrix(path: Path) -> sp.csc_matrix:
    with gzip.open(path, "rb") as handle:
        return sp.csc_matrix(mmread(handle))


def integrated(infer: str, copy: str) -> tuple[str, bool]:
    if infer == "high_cnv" and copy == "aneuploid":
        return "malignant_consensus", True
    if infer == "high_cnv" and copy in {"not_defined", "not_estimable", "not_eligible"}:
        return "malignant_infercnv_only", True
    if infer == "not_eligible" and copy == "aneuploid":
        return "malignant_copykat_only", True
    if infer == "low_cnv" and copy == "diploid":
        return "non_malignant_consensus", False
    if (infer == "high_cnv" and copy == "diploid") or (infer == "low_cnv" and copy == "aneuploid"):
        return "method_discordant", False
    return "indeterminate", False


def main() -> None:
    if VERIFY.exists() and any(VERIFY.iterdir()):
        raise RuntimeError(f"Refusing to overwrite non-empty verification directory: {VERIFY}")
    VERIFY.mkdir(parents=True, exist_ok=True)
    checks: list[dict[str, str]] = []
    gates: list[dict[str, str]] = []

    def check(name: str, passed: bool, evidence: str) -> None:
        checks.append({"check": name, "status": "PASS" if passed else "FAIL", "evidence": evidence})
        if not passed:
            raise RuntimeError(f"Verification failed: {name}: {evidence}")

    def gate(stage: str) -> None:
        result = subprocess.run([str(MEMORY_GATE)], text=True, capture_output=True)
        output = " ".join((result.stdout + " " + result.stderr).split())
        gates.append({"stage": stage, "exit_status": str(result.returncode), "gate_output": output})
        print(output, flush=True)
        check(f"memory_gate_{stage}", result.returncode == 0 and "PASS" in output, output)

    gate("start")
    check("sensitivity_run_complete", (SENS / "RUN_COMPLETE").is_file(), str(SENS / "RUN_COMPLETE"))
    check("summary_complete", (SUMMARY / "RUN_COMPLETE").is_file(), str(SUMMARY / "RUN_COMPLETE"))
    check("control_summary_complete", (CONTROL_SUMMARY / "RUN_COMPLETE").is_file(), str(CONTROL_SUMMARY / "RUN_COMPLETE"))
    registry = pd.read_csv(SENS / "inputs/mixed_cell_registry_v1.tsv", sep="\t", dtype=str)
    check("mixed_total", len(registry) == 918, str(len(registry)))
    check("mixed_tumor", int(registry["source_tissue"].eq("Tumor").sum()) == 2, str(int(registry["source_tissue"].eq("Tumor").sum())))
    check("mixed_normal", int(registry["source_tissue"].eq("Normal").sum()) == 916, str(int(registry["source_tissue"].eq("Normal").sum())))

    for method, unit, n_query, n_ref, mixed_id in UNITS:
        gate(f"verify_input_{method}_{unit}")
        primary_dir = PRIMARY_INPUT / method / unit
        sensitivity_dir = SENS / "inputs" / method / unit
        p_map = pd.read_csv(primary_dir / "cell_id_mapping_v1.tsv.gz", sep="\t", dtype=str)
        s_map = pd.read_csv(sensitivity_dir / "cell_id_mapping_v1.tsv.gz", sep="\t", dtype=str)
        check(f"{method}_{unit}_mapping_plus_one", len(s_map) == len(p_map) + 1, f"{len(p_map)}->{len(s_map)}")
        check(f"{method}_{unit}_mapping_prefix_exact", s_map.iloc[:-1].equals(p_map), "all primary rows exact")
        check(f"{method}_{unit}_mixed_last", s_map.iloc[-1]["exact_cell_id"] == mixed_id and s_map.iloc[-1]["role"] == "query", str(s_map.iloc[-1].to_dict()))
        check(f"{method}_{unit}_query_count", int(s_map["role"].eq("query").sum()) == n_query + 1, str(int(s_map["role"].eq("query").sum())))
        check(f"{method}_{unit}_reference_count", int(s_map["role"].eq("reference").sum()) == n_ref, str(int(s_map["role"].eq("reference").sum())))
        p_gene = pd.read_csv(primary_dir / "gene_order_v1.tsv.gz", sep="\t", dtype=str)
        s_gene = pd.read_csv(sensitivity_dir / "gene_order_v1.tsv.gz", sep="\t", dtype=str)
        check(f"{method}_{unit}_gene_order_exact", s_gene.equals(p_gene), f"n={len(s_gene)}")
        p_matrix = read_matrix(primary_dir / "counts_genes_by_cells_v1.mtx.gz")
        s_matrix = read_matrix(sensitivity_dir / "counts_genes_by_cells_v1.mtx.gz")
        check(f"{method}_{unit}_matrix_shape_plus_one", s_matrix.shape == (p_matrix.shape[0], p_matrix.shape[1] + 1), f"{p_matrix.shape}->{s_matrix.shape}")
        diff = s_matrix[:, :-1] - p_matrix
        diff.eliminate_zeros()
        check(f"{method}_{unit}_primary_matrix_exact", diff.nnz == 0, f"difference_nnz={diff.nnz}")
        check(f"{method}_{unit}_mixed_counts_nonzero", s_matrix[:, -1].nnz > 0, f"nnz={s_matrix[:, -1].nnz}")
        output_dir = SENS / method / unit
        check(f"{method}_{unit}_run_complete", (output_dir / "RUN_COMPLETE").is_file(), str(output_dir))
        status = json.loads((output_dir / "status.json").read_text(encoding="utf-8"))
        check(f"{method}_{unit}_status", status.get("status") == "completed", str(status.get("status")))
        run_summary = pd.read_csv(output_dir / "run_summary_v2.tsv", sep="\t", dtype=str).set_index("parameter")["value"]
        query_key = "n_query_cells" if method == "infercnv" else "n_query_input"
        thread_key = "threads" if method == "infercnv" else "n_cores"
        check(f"{method}_{unit}_run_query_count", int(run_summary[query_key]) == n_query + 1, str(run_summary[query_key]))
        check(f"{method}_{unit}_low_resource", int(run_summary[thread_key]) == 4, str(run_summary[thread_key]))
        gate_table = pd.read_csv(output_dir / "memory_gate_log.tsv", sep="\t", dtype=str)
        check(f"{method}_{unit}_all_internal_gates", gate_table["gate_output"].str.contains("PASS", regex=False).all(), f"n={len(gate_table)}")

    mixed = pd.read_csv(SUMMARY / "tumor_mixed_cnv_calls_v1.tsv", sep="\t")
    check("mixed_call_rows", len(mixed) == 2 and not mixed["exact_cell_id"].duplicated().any(), str(len(mixed)))
    check("mixed_expected_ids", set(mixed["exact_cell_id"]) == {row[4] for row in UNITS}, ",".join(mixed["exact_cell_id"]))
    for row in mixed.itertuples(index=False):
        status, malignant = integrated(row.infercnv_call, row.copykat_call)
        check(f"mixed_{row.sample_id}_integrated_status", row.malignancy_status_mixed_sensitivity_v1 == status, f"{row.infercnv_call}/{row.copykat_call}->{status}")
        check(f"mixed_{row.sample_id}_integrated_bool", bool(row.malignant_conservative_mixed_sensitivity_v1) == malignant, str(malignant))
        check(f"mixed_{row.sample_id}_not_primary", not bool(row.included_in_primary_malignancy_summary_v1), "false")

    normal = pd.read_csv(SUMMARY / "normal_mixed_scope_v1.tsv.gz", sep="\t", dtype=str)
    check("normal_scope_rows", len(normal) == 916, str(len(normal)))
    check("normal_scope_status", normal["cnv_sensitivity_status_v1"].eq("not_applicable_normal_tissue").all(), "916/916")
    stability = pd.read_csv(SUMMARY / "primary_clear_call_stability_v1.tsv", sep="\t")
    check("stability_rows", len(stability) == 4, str(len(stability)))
    check("stability_denominators", set(stability["n_primary_clear_query"]) == {1508, 2381}, str(sorted(stability["n_primary_clear_query"].tolist())))
    check("stability_reconciles", (stability["n_exact_call_agreement"] + stability["n_call_changed"] == stability["n_primary_clear_query"]).all(), "4/4")
    check("stability_fraction_range", stability["exact_call_agreement_fraction"].between(0, 1).all(), str(stability["exact_call_agreement_fraction"].tolist()))

    control_stability = pd.read_csv(CONTROL_SUMMARY / "infercnv_paired_control_call_stability_v2.tsv", sep="\t")
    check("control_stability_rows", len(control_stability) == 6, str(len(control_stability)))
    primary_vs_control = control_stability.loc[
        control_stability["left_run"].eq("primary_8threads")
        & control_stability["right_run"].eq("control_no_add_4threads")
    ]
    check(
        "primary_vs_no_add_control_exact",
        len(primary_vs_control) == 2
        and primary_vs_control["exact_call_agreement_fraction"].eq(1.0).all(),
        primary_vs_control[["unit_id", "exact_call_agreement_fraction"]].to_dict("records").__str__(),
    )
    control_vs_plus_one = control_stability.loc[
        control_stability["left_run"].eq("control_no_add_4threads")
        & control_stability["right_run"].eq("plus_one_mixed_4threads")
    ].set_index("unit_id")
    infer_stability = stability.loc[stability["method"].eq("infercnv")].set_index("unit_id")
    check(
        "control_vs_plus_one_matches_primary_comparison",
        set(control_vs_plus_one.index) == set(infer_stability.index)
        and np.allclose(
            control_vs_plus_one.loc[infer_stability.index, "exact_call_agreement_fraction"],
            infer_stability["exact_call_agreement_fraction"],
            rtol=0,
            atol=0,
        ),
        control_vs_plus_one["exact_call_agreement_fraction"].to_dict().__str__(),
    )
    check(
        "infercnv_plus_one_background_unstable",
        control_vs_plus_one["exact_call_agreement_fraction"].lt(0.95).all(),
        control_vs_plus_one["exact_call_agreement_fraction"].to_dict().__str__(),
    )
    copy_stability = stability.loc[stability["method"].eq("copykat")]
    check(
        "copykat_plus_one_background_stable",
        len(copy_stability) == 2 and copy_stability["exact_call_agreement_fraction"].ge(0.95).all(),
        copy_stability[["unit_id", "exact_call_agreement_fraction"]].to_dict("records").__str__(),
    )
    mixed_control = pd.read_csv(CONTROL_SUMMARY / "tumor_mixed_cnv_calls_with_control_v2.tsv", sep="\t")
    check("mixed_control_rows", len(mixed_control) == 2, str(len(mixed_control)))
    check(
        "mixed_control_unstable_infercnv",
        mixed_control["infercnv_sensitivity_reliability_v2"].eq("unstable_clear_background").all(),
        mixed_control["infercnv_sensitivity_reliability_v2"].tolist().__str__(),
    )
    check(
        "mixed_control_stable_copykat",
        mixed_control["copykat_sensitivity_reliability_v2"].eq("stable_clear_background").all(),
        mixed_control["copykat_sensitivity_reliability_v2"].tolist().__str__(),
    )
    check(
        "mixed_ineligible_for_primary_or_subclone",
        not mixed_control["eligible_for_primary_or_subclone_selection_v2"].astype(bool).any(),
        "2/2 false",
    )

    primary_infer = pd.read_csv(PRIMARY_SUMMARY / "infercnv_query_cell_calls_v1.tsv.gz", sep="\t", dtype=str)
    sensitivity_infer = pd.read_csv(SUMMARY / "infercnv_query_cell_calls_mixed_sensitivity_v1.tsv.gz", sep="\t", dtype=str)
    primary_copy = pd.read_csv(PRIMARY_SUMMARY / "copykat_query_cell_calls_v1.tsv.gz", sep="\t", dtype=str)
    sensitivity_copy = pd.read_csv(SUMMARY / "copykat_query_cell_calls_mixed_sensitivity_v1.tsv.gz", sep="\t", dtype=str)
    for method, unit_column, call_column, primary, sensitivity in [
        ("infercnv", "infercnv_unit", "infercnv_call", primary_infer, sensitivity_infer),
        ("copykat", "copykat_unit", "copykat_call", primary_copy, sensitivity_copy),
    ]:
        for unit in [row[1] for row in UNITS if row[0] == method]:
            p = primary.loc[primary[unit_column].eq(unit), ["exact_cell_id", call_column]]
            s = sensitivity.loc[sensitivity[unit_column].eq(unit) & sensitivity["cell_type_broad_v1"].eq("Epithelial cell"), ["exact_cell_id", call_column]]
            merged = p.merge(s, on="exact_cell_id", suffixes=("_p", "_s"), validate="one_to_one")
            exact = int(merged[f"{call_column}_p"].eq(merged[f"{call_column}_s"]).sum())
            recorded = stability.loc[(stability["method"] == method) & (stability["unit_id"] == unit)].iloc[0]
            check(f"{method}_{unit}_stability_independent", len(merged) == int(recorded["n_primary_clear_query"]) and exact == int(recorded["n_exact_call_agreement"]), f"n={len(merged)} exact={exact}")

    figure = SUMMARY / "mixed_cnv_sensitivity_summary_v1.png"
    with Image.open(figure) as image:
        width, height = image.size
    check("figure_resolution", width >= 2000 and height >= 800, f"{width}x{height}")
    manifest = pd.read_csv(SUMMARY / "file_manifest_v1.tsv", sep="\t", dtype=str)
    for row in manifest.itertuples(index=False):
        path = SUMMARY / row.file
        check(f"summary_hash_{row.file}", path.is_file() and sha256(path) == row.sha256, row.sha256)
    control_manifest = pd.read_csv(CONTROL_SUMMARY / "file_manifest_v2.tsv", sep="\t", dtype=str)
    for row in control_manifest.itertuples(index=False):
        path = CONTROL_SUMMARY / row.file
        check(f"control_summary_hash_{row.file}", path.is_file() and sha256(path) == row.sha256, row.sha256)
    primary_manifest = pd.read_csv(PRIMARY_SUMMARY / "file_manifest_v1.tsv", sep="\t", dtype=str)
    for row in primary_manifest.itertuples(index=False):
        path = Path(row.path)
        check(f"primary_unchanged_{path.name}", path.is_file() and sha256(path) == row.sha256, row.sha256)
    for script in [
        ROOT / "prepare_part1_mixed_cnv_sensitivity_inputs_v1.py",
        ROOT / "run_part1_mixed_cnv_sensitivity_v1.sh",
        ROOT / "summarize_part1_mixed_cnv_sensitivity_v1.py",
        ROOT / "run_part1_mixed_cnv_sensitivity_summary_v1.sh",
        ROOT / "summarize_part1_mixed_cnv_control_v1.py",
        ROOT / "run_part1_mixed_cnv_control_summary_v2.sh",
    ]:
        text = script.read_text(encoding="utf-8")
        check(f"part4_not_referenced_{script.name}", "part4" not in text.lower(), "no Part 4 path/reference")

    pd.DataFrame(checks).to_csv(VERIFY / "verification_checks_v1.tsv", sep="\t", index=False, lineterminator="\n")
    pd.DataFrame(gates).to_csv(VERIFY / "memory_gate_log_v1.tsv", sep="\t", index=False, lineterminator="\n")
    critical = []
    for path in sorted([p for p in SENS.rglob("*") if p.is_file()]):
        if path.name.endswith((".rds", ".jpeg")) or "infercnv_run" in path.parts or "copykat_run" in path.parts:
            continue
        critical.append({"path": str(path.resolve()), "size_bytes": path.stat().st_size, "sha256": sha256(path)})
    pd.DataFrame(critical).to_csv(VERIFY / "critical_file_manifest_v1.tsv", sep="\t", index=False, lineterminator="\n")
    metadata = {
        "schema_version": "part1_mixed_cnv_sensitivity_verification_v1",
        "completed_at": datetime.now().astimezone().isoformat(),
        "n_checks": len(checks),
        "n_pass": sum(row["status"] == "PASS" for row in checks),
        "n_fail": sum(row["status"] == "FAIL" for row in checks),
        "n_critical_files_hashed": len(critical),
    }
    (VERIFY / "run_metadata_v1.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    (VERIFY / "RUN_COMPLETE").write_text(datetime.now().astimezone().isoformat() + "\n", encoding="utf-8")
    print(json.dumps(metadata, indent=2))
    print("PART1_MIXED_CNV_SENSITIVITY_VERIFICATION_V1_PASS")


if __name__ == "__main__":
    main()
