#!/usr/bin/env python3
from __future__ import annotations

import csv
import hashlib
import json
import subprocess
from datetime import datetime
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from PIL import Image


PROJECT = Path("/mnt/BioDisk/ludingyu/project26210")
ROOT = PROJECT / "results/part1_rebuild_v1/03_malignancy_cnv"
INPUT_H5AD = PROJECT / "results/part1_rebuild_v1/02_integration_annotation/cell_annotation_v1/part1_scvi_integrated_cell_annotated_v1.h5ad"
SUMMARY = ROOT / "malignancy_cnv_summary_v1"
OUT = ROOT / "malignancy_cnv_summary_verification_v1"
MEMORY_GATE = PROJECT / "check_memory_gate.sh"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str)


def expected_status(infer_call: str, copy_call: str) -> tuple[str, bool]:
    if infer_call == "high_cnv" and copy_call == "aneuploid":
        return "malignant_consensus", True
    if infer_call == "high_cnv" and copy_call in {"not_defined", "not_estimable", "not_eligible"}:
        return "malignant_infercnv_only", True
    if infer_call == "not_eligible" and copy_call == "aneuploid":
        return "malignant_copykat_only", True
    if infer_call == "low_cnv" and copy_call == "diploid":
        return "non_malignant_consensus", False
    if (infer_call == "high_cnv" and copy_call == "diploid") or (infer_call == "low_cnv" and copy_call == "aneuploid"):
        return "method_discordant", False
    return "indeterminate", False


def normalized_rows(frame: pd.DataFrame, columns: list[str]) -> list[tuple[str, ...]]:
    work = frame[columns].fillna("").astype(str)
    return sorted(tuple(row) for row in work.itertuples(index=False, name=None))


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    checks: list[dict[str, str]] = []

    def check(name: str, passed: bool, detail: str) -> None:
        checks.append({"check": name, "status": "PASS" if passed else "FAIL", "detail": detail})

    gate = subprocess.run([str(MEMORY_GATE)], text=True, capture_output=True)
    gate_output = " ".join((gate.stdout + " " + gate.stderr).split())
    check("verification:memory_gate", gate.returncode == 0 and "PASS" in gate_output, gate_output)
    check("summary:run_complete", (SUMMARY / "RUN_COMPLETE").is_file(), str(SUMMARY / "RUN_COMPLETE"))

    calls_path = SUMMARY / "tumor_epithelial_malignancy_calls_v1.tsv.gz"
    calls = pd.read_csv(calls_path, sep="\t", low_memory=False)
    check("calls:row_count", len(calls) == 38670, f"n={len(calls)}")
    check("calls:unique_exact_ids", calls["exact_cell_id"].nunique() == 38670, f"unique={calls['exact_cell_id'].nunique()}")
    check("calls:tumor_only", set(calls["source_tissue"].astype(str)) == {"Tumor"}, str(sorted(calls["source_tissue"].astype(str).unique())))
    check("calls:epithelial_only", set(calls["cell_type_broad_v1"].astype(str)) == {"Epithelial cell"}, str(sorted(calls["cell_type_broad_v1"].astype(str).unique())))
    check("calls:no_na_status", not calls[["infercnv_call", "copykat_call", "malignancy_status_v1", "malignancy_evidence_v1"]].isna().any().any(), "status columns")

    adata = ad.read_h5ad(INPUT_H5AD, backed="r")
    source = adata.obs
    expected_ids = set(source.index[source["cell_type_broad_v1"].astype(str).eq("Epithelial cell") & source["source_tissue"].astype(str).eq("Tumor")].astype(str))
    normal_count = int((source["cell_type_broad_v1"].astype(str).eq("Epithelial cell") & source["source_tissue"].astype(str).eq("Normal")).sum())
    adata.file.close()
    check("calls:exact_source_id_set", set(calls["exact_cell_id"].astype(str)) == expected_ids, f"source={len(expected_ids)};calls={len(calls)}")
    check("source:normal_epithelial_count", normal_count == 23982, f"n={normal_count}")

    predicted = [expected_status(i, c) for i, c in zip(calls["infercnv_call"], calls["copykat_call"])]
    predicted_status = [x[0] for x in predicted]
    predicted_bool = [x[1] for x in predicted]
    recorded_bool = calls["malignant_conservative_v1"].astype(str).str.lower().map({"true": True, "false": False})
    check("calls:integrated_status_truth_table", predicted_status == calls["malignancy_status_v1"].tolist(), "independent truth-table replay")
    check("calls:conservative_boolean_truth_table", predicted_bool == recorded_bool.tolist(), "independent truth-table replay")

    infer_cells = pd.read_csv(SUMMARY / "infercnv_query_cell_calls_v1.tsv.gz", sep="\t", low_memory=False)
    infer_groups = pd.read_csv(SUMMARY / "infercnv_group_burden_v1.tsv", sep="\t", low_memory=False)
    infer_params = pd.read_csv(SUMMARY / "infercnv_threshold_parameters_v1.tsv", sep="\t")
    check("infercnv:unit_count", infer_cells["infercnv_unit"].nunique() == 34 and len(infer_params) == 34, f"cells_units={infer_cells['infercnv_unit'].nunique()};params={len(infer_params)}")
    check("infercnv:query_cell_count", len(infer_cells) == 38547 and infer_cells["exact_cell_id"].nunique() == 38547, f"n={len(infer_cells)};unique={infer_cells['exact_cell_id'].nunique()}")
    param_index = infer_params.set_index("infercnv_unit")
    for unit, unit_groups in infer_groups.groupby("infercnv_unit", observed=True):
        refs = unit_groups.loc[unit_groups["role"].eq("reference"), "cnv_gene_burden"].astype(float).to_numpy()
        median = float(np.median(refs))
        mad = float(np.median(np.abs(refs - median)) * 1.4826)
        q95 = float(np.quantile(refs, 0.95))
        threshold = max(median + 3 * mad, q95)
        recorded = float(param_index.loc[unit, "high_cnv_threshold"])
        check(f"infercnv:{unit}:threshold", np.isclose(recorded, threshold, rtol=0, atol=1e-8), f"recorded={recorded};recomputed={threshold}")
        query_groups = unit_groups.loc[unit_groups["role"].eq("query")]
        expected_group_calls = np.where(query_groups["cnv_gene_burden"].astype(float).gt(threshold), "high_cnv", "low_cnv")
        check(f"infercnv:{unit}:group_calls", np.array_equal(expected_group_calls, query_groups["infercnv_call"].to_numpy()), f"groups={len(query_groups)}")

    infer_compare = calls.loc[calls["infercnv_call"].ne("not_eligible"), ["exact_cell_id", "infercnv_call", "infercnv_unit", "infercnv_group"]]
    infer_source = infer_cells.rename(columns={"cell_group_name": "infercnv_group"})[["exact_cell_id", "infercnv_call", "infercnv_unit", "infercnv_group"]]
    check("infercnv:integrated_mapping_exact", normalized_rows(infer_compare, list(infer_compare.columns)) == normalized_rows(infer_source, list(infer_source.columns)), f"integrated={len(infer_compare)};source={len(infer_source)}")

    copy_cells = pd.read_csv(SUMMARY / "copykat_query_cell_calls_v1.tsv.gz", sep="\t", low_memory=False)
    check("copykat:query_cell_count", len(copy_cells) == 38191 and copy_cells["exact_cell_id"].nunique() == 38191, f"n={len(copy_cells)};unique={copy_cells['exact_cell_id'].nunique()}")
    check("copykat:unit_count", copy_cells["copykat_unit"].nunique() == 41, f"units={copy_cells['copykat_unit'].nunique()}")
    copy_compare = calls.loc[calls["copykat_call"].ne("not_eligible"), ["exact_cell_id", "copykat_call", "copykat_unit"]]
    copy_source = copy_cells[["exact_cell_id", "copykat_call", "copykat_unit"]]
    check("copykat:integrated_mapping_exact", normalized_rows(copy_compare, list(copy_compare.columns)) == normalized_rows(copy_source, list(copy_source.columns)), f"integrated={len(copy_compare)};source={len(copy_source)}")
    check("copykat:not_estimable_units", set(copy_cells.loc[copy_cells["copykat_call"].eq("not_estimable"), "copykat_unit"]) == {"T2_1", "T2_2", "T3_1"}, str(sorted(copy_cells.loc[copy_cells['copykat_call'].eq('not_estimable'), 'copykat_unit'].unique())))

    cross_expected = calls.groupby(["infercnv_call", "copykat_call"], observed=True).size().rename("n_cells").reset_index()
    cross_recorded = pd.read_csv(SUMMARY / "method_cross_table_v1.tsv", sep="\t")
    check("summary:cross_table_exact", normalized_rows(cross_expected.astype(str), ["infercnv_call", "copykat_call", "n_cells"]) == normalized_rows(cross_recorded.astype(str), ["infercnv_call", "copykat_call", "n_cells"]), f"rows={len(cross_expected)}")

    status_expected = calls.groupby(["malignancy_status_v1", "malignant_conservative_v1"], observed=True).size().rename("n_cells").reset_index()
    status_recorded = pd.read_csv(SUMMARY / "malignancy_status_counts_v1.tsv", sep="\t")
    check("summary:status_counts_exact", normalized_rows(status_expected.astype(str), list(status_expected.columns)) == normalized_rows(status_recorded.astype(str), list(status_recorded.columns)), f"rows={len(status_expected)}")
    status_map = dict(zip(status_recorded["malignancy_status_v1"], status_recorded["n_cells"].astype(int)))
    expected_status_map = {
        "malignant_consensus": 10844,
        "malignant_infercnv_only": 1106,
        "method_discordant": 22216,
        "non_malignant_consensus": 498,
        "indeterminate": 4006,
    }
    check("summary:expected_status_values", status_map == expected_status_map, f"actual={status_map}")

    patient = pd.read_csv(SUMMARY / "patient_malignancy_summary_v1.tsv", sep="\t")
    sample = pd.read_csv(SUMMARY / "sample_malignancy_summary_v1.tsv", sep="\t")
    check("summary:patient_count", len(patient) == 38 and patient["n_tumor_epithelial"].sum() == 38670, f"patients={len(patient)};cells={patient['n_tumor_epithelial'].sum()}")
    check("summary:sample_count", len(sample) == 52 and sample["n_tumor_epithelial"].sum() == 38670, f"samples={len(sample)};cells={sample['n_tumor_epithelial'].sum()}")
    check("summary:conservative_cell_count", int(recorded_bool.sum()) == 11950 and int(patient["n_malignant_conservative"].sum()) == 11950, f"cells={int(recorded_bool.sum())};patient_sum={int(patient['n_malignant_conservative'].sum())}")
    check("summary:patients_with_conservative_calls", int((patient["n_malignant_conservative"] > 0).sum()) == 26, f"n={int((patient['n_malignant_conservative'] > 0).sum())}")

    definitive = calls["infercnv_call"].isin(["high_cnv", "low_cnv"]) & calls["copykat_call"].isin(["aneuploid", "diploid"])
    agreement = ((calls.loc[definitive, "infercnv_call"] == "high_cnv") == (calls.loc[definitive, "copykat_call"] == "aneuploid"))
    check("summary:definitive_overlap", int(definitive.sum()) == 33558, f"n={int(definitive.sum())}")
    check("summary:definitive_agreement", int(agreement.sum()) == 11342 and np.isclose(float(agreement.mean()), 0.33798200131116274), f"n={int(agreement.sum())};fraction={float(agreement.mean())}")

    gate_log = pd.read_csv(SUMMARY / "memory_gate_log_v1.tsv", sep="\t")
    check("summary:all_memory_gates_pass", len(gate_log) == 77 and gate_log["exit_status"].astype(int).eq(0).all() and gate_log["gate_output"].str.contains("PASS", regex=False).all(), f"rows={len(gate_log)}")

    for name in ["method_cross_table_v1", "patient_malignancy_composition_v1"]:
        png = SUMMARY / f"{name}.png"
        pdf = SUMMARY / f"{name}.pdf"
        readable = False
        dimensions = ""
        try:
            with Image.open(png) as image:
                image.verify()
            with Image.open(png) as image:
                dimensions = f"{image.width}x{image.height}"
                readable = image.width >= 1200 and image.height >= 800
        except Exception as exc:
            dimensions = str(exc)
        check(f"figure:{name}:readable_png", readable, dimensions)
        check(f"figure:{name}:pdf_exists", pdf.is_file() and pdf.stat().st_size > 10000, f"size={pdf.stat().st_size if pdf.exists() else -1}")

    manifest_path = SUMMARY / "file_manifest_v1.tsv"
    manifest = read_tsv(manifest_path)
    manifest_ok = True
    for row in manifest.itertuples(index=False):
        path = Path(row.path)
        if not path.is_file() or str(path.stat().st_size) != row.size_bytes or sha256(path) != row.sha256:
            manifest_ok = False
            break
    check("summary:file_manifest_integrity", manifest_ok, f"files={len(manifest)}")

    metadata = json.loads((SUMMARY / "run_metadata_v1.json").read_text())
    check("summary:input_h5ad_hash", metadata.get("input_h5ad_sha256") == sha256(INPUT_H5AD), metadata.get("input_h5ad_sha256", ""))
    check("summary:no_biological_independence_claim", "not independent biological replicates" in metadata.get("interpretation", ""), metadata.get("interpretation", ""))

    checks_file = OUT / "verification_checks_v1.tsv"
    with checks_file.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["check", "status", "detail"], delimiter="\t")
        writer.writeheader()
        writer.writerows(checks)

    critical = sorted({
        calls_path,
        SUMMARY / "infercnv_query_cell_calls_v1.tsv.gz",
        SUMMARY / "infercnv_group_burden_v1.tsv",
        SUMMARY / "infercnv_threshold_parameters_v1.tsv",
        SUMMARY / "copykat_query_cell_calls_v1.tsv.gz",
        SUMMARY / "method_cross_table_v1.tsv",
        SUMMARY / "malignancy_status_counts_v1.tsv",
        SUMMARY / "patient_malignancy_summary_v1.tsv",
        SUMMARY / "sample_malignancy_summary_v1.tsv",
        SUMMARY / "method_coverage_summary_v1.tsv",
        SUMMARY / "memory_gate_log_v1.tsv",
        SUMMARY / "method_cross_table_v1.png",
        SUMMARY / "patient_malignancy_composition_v1.png",
        SUMMARY / "run_metadata_v1.json",
        SUMMARY / "RUN_COMPLETE",
        ROOT / "summarize_part1_malignancy_cnv_v1.py",
    })
    with (OUT / "critical_file_manifest_v1.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["path", "size_bytes", "sha256"], delimiter="\t")
        writer.writeheader()
        for path in critical:
            writer.writerow({"path": str(path), "size_bytes": path.stat().st_size, "sha256": sha256(path)})

    failed = [row for row in checks if row["status"] != "PASS"]
    run_metadata = {
        "schema_version": "part1_malignancy_cnv_summary_verification_v1",
        "created_at": datetime.now().astimezone().isoformat(),
        "checks_total": len(checks),
        "checks_passed": len(checks) - len(failed),
        "checks_failed": len(failed),
        "critical_files_hashed": len(critical),
        "manual_visual_review_separate": True,
    }
    (OUT / "run_metadata_v1.json").write_text(json.dumps(run_metadata, indent=2, sort_keys=True) + "\n")
    if failed:
        (OUT / "RUN_FAILED").write_text(datetime.now().astimezone().isoformat() + "\n")
        raise SystemExit(f"PART1_MALIGNANCY_CNV_SUMMARY_VERIFICATION_FAILED checks={len(checks)} failed={len(failed)}")
    (OUT / "RUN_COMPLETE").write_text(datetime.now().astimezone().isoformat() + "\n")
    print(f"PART1_MALIGNANCY_CNV_SUMMARY_VERIFICATION_ALL_PASS checks={len(checks)} files={len(critical)}")


if __name__ == "__main__":
    main()
