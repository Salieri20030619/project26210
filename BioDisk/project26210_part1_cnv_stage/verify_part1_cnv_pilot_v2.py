#!/usr/bin/env python3
"""Independent structural and identity verification for CNV pilot v2."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from datetime import datetime
from pathlib import Path

import anndata as ad
import pandas as pd


PROJECT = Path("/mnt/BioDisk/ludingyu/project26210")
INPUT_H5AD = PROJECT / (
    "results/part1_rebuild_v1/02_integration_annotation/cell_annotation_v1/"
    "part1_scvi_integrated_cell_annotated_v1.h5ad"
)
DEFAULT_PILOT = Path(
    "/mnt/BioDisk/wangmingyuan/project26210_part1_cnv_stage/cnv_pilot_v2"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pilot-dir", type=Path, default=DEFAULT_PILOT)
    return parser.parse_args()


def sha256_file(path: Path, block_size: int = 16 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_size):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    args = parse_args()
    pilot = args.pilot_dir
    verify_dir = pilot / "verification_v1"
    verify_dir.mkdir(parents=True, exist_ok=True)

    infer_input = pilot / "inputs/infercnv/GSE156632__RCC1"
    copy_input = pilot / "inputs/copykat/RCC1t"
    infer_out = pilot / "infercnv/GSE156632__RCC1"
    copy_out = pilot / "copykat/RCC1t"
    hmm_group_file = infer_out / (
        "infercnv_run/17_HMM_predHMMi6.leiden.hmm_mode-subclusters.cell_groupings"
    )
    hmm_regions_file = infer_out / (
        "infercnv_run/17_HMM_predHMMi6.leiden.hmm_mode-subclusters.pred_cnv_regions.dat"
    )

    adata = ad.read_h5ad(INPUT_H5AD, backed="r")
    current_ids = set(adata.obs_names.astype(str))
    current_obs = adata.obs[["sample_id", "curated_patient_id", "source_tissue", "cell_type_broad_v1"]].copy()
    current_obs.index = adata.obs_names.astype(str)

    infer_map_in = pd.read_csv(infer_input / "cell_id_mapping_v1.tsv.gz", sep="\t")
    copy_map_in = pd.read_csv(copy_input / "cell_id_mapping_v1.tsv.gz", sep="\t")
    infer_map_out = pd.read_csv(infer_out / "cell_id_mapping_v1.tsv", sep="\t")
    copy_map_out = pd.read_csv(copy_out / "cell_id_mapping_v1.tsv", sep="\t")
    hmm_groups = pd.read_csv(hmm_group_file, sep="\t")
    hmm_regions = pd.read_csv(hmm_regions_file, sep="\t")
    predictions = pd.read_csv(copy_out / "copykat_predictions_mapped_v1.tsv", sep="\t")

    infer_status = json.loads((infer_out / "status.json").read_text())
    copy_status = json.loads((copy_out / "status.json").read_text())
    input_status = json.loads((pilot / "inputs/run_metadata_v1.json").read_text())
    log_text = (pilot / "cnv_pilot_v2.log").read_text(errors="replace")

    checks: list[dict[str, object]] = []

    def check(name: str, passed: bool, observed: object, expected: object) -> None:
        checks.append({"check": name, "pass": bool(passed), "observed": str(observed), "expected": str(expected)})

    expected_files = [
        pilot / "inputs/RUN_COMPLETE",
        infer_out / "RUN_COMPLETE",
        copy_out / "RUN_COMPLETE",
        infer_out / "infercnv_object_post_run_v1.rds",
        copy_out / "copykat_result_v1.rds",
        hmm_group_file,
        hmm_regions_file,
        copy_out / "copykat_predictions_mapped_v1.tsv",
    ]
    check("all_required_files_present", all(x.is_file() for x in expected_files), sum(x.is_file() for x in expected_files), len(expected_files))
    check("input_preparation_status", input_status.get("status") == "completed", input_status.get("status"), "completed")
    check("infercnv_status", infer_status.get("status") == "completed", infer_status.get("status"), "completed")
    check("copykat_status", copy_status.get("status") == "completed", copy_status.get("status"), "completed")
    check("no_error_files", not any(pilot.rglob("ERROR.txt")), len(list(pilot.rglob("ERROR.txt"))), 0)
    check("pilot_terminal_token", "PART1_CNV_PILOT_V2_ALL_PASS" in log_text, "PART1_CNV_PILOT_V2_ALL_PASS" in log_text, True)
    check("infercnv_terminal_token", "PART1_INFERCNV_UNIT_V1_PASS" in log_text, "PART1_INFERCNV_UNIT_V1_PASS" in log_text, True)
    check("copykat_terminal_token", "PART1_COPYKAT_UNIT_V1_PASS" in log_text, "PART1_COPYKAT_UNIT_V1_PASS" in log_text, True)
    check("no_traceback", "Traceback" not in log_text, "Traceback" in log_text, False)
    check("no_r_hard_error", "执行已停止" not in log_text and "Execution halted" not in log_text, "hard error present" if ("执行已停止" in log_text or "Execution halted" in log_text) else "absent", "absent")

    safe_pattern = re.compile(r"^[QR][0-9]{7}$")
    for label, mapping, expected_query, expected_reference in (
        ("infercnv", infer_map_in, 69, 200),
        ("copykat", copy_map_in, 69, 277),
    ):
        check(f"{label}_tool_ids_unique", not mapping["tool_id"].duplicated().any(), mapping["tool_id"].duplicated().sum(), 0)
        check(f"{label}_exact_ids_unique", not mapping["exact_cell_id"].duplicated().any(), mapping["exact_cell_id"].duplicated().sum(), 0)
        check(f"{label}_safe_id_format", mapping["tool_id"].map(lambda x: bool(safe_pattern.fullmatch(str(x)))).all(), "all" if mapping["tool_id"].map(lambda x: bool(safe_pattern.fullmatch(str(x)))).all() else "not_all", "all")
        check(f"{label}_all_ids_in_h5ad", set(mapping["exact_cell_id"]).issubset(current_ids), len(set(mapping["exact_cell_id"]) - current_ids), 0)
        check(f"{label}_query_count", int((mapping["role"] == "query").sum()) == expected_query, int((mapping["role"] == "query").sum()), expected_query)
        check(f"{label}_reference_count", int((mapping["role"] == "reference").sum()) == expected_reference, int((mapping["role"] == "reference").sum()), expected_reference)
        joined = mapping.set_index("exact_cell_id").join(current_obs, rsuffix="_h5ad")
        query = joined[joined["role"] == "query"]
        reference = joined[joined["role"] == "reference"]
        check(f"{label}_query_is_tumor_epithelial", ((query["source_tissue_h5ad"] == "Tumor") & (query["cell_type_broad_v1_h5ad"] == "Epithelial cell")).all(), int(((query["source_tissue_h5ad"] == "Tumor") & (query["cell_type_broad_v1_h5ad"] == "Epithelial cell")).sum()), len(query))
        check(f"{label}_reference_is_normal_epithelial", ((reference["source_tissue_h5ad"] == "Normal") & (reference["cell_type_broad_v1_h5ad"] == "Epithelial cell")).all(), int(((reference["source_tissue_h5ad"] == "Normal") & (reference["cell_type_broad_v1_h5ad"] == "Epithelial cell")).sum()), len(reference))

    check("infercnv_output_mapping_exact", set(infer_map_in["exact_cell_id"]) == set(infer_map_out["exact_cell_id"]), len(set(infer_map_in["exact_cell_id"]) ^ set(infer_map_out["exact_cell_id"])), 0)
    check("copykat_output_mapping_exact", set(copy_map_in["exact_cell_id"]) == set(copy_map_out["exact_cell_id"]), len(set(copy_map_in["exact_cell_id"]) ^ set(copy_map_out["exact_cell_id"])), 0)
    check("hmm_groupings_all_cells", set(hmm_groups["cell"]) == set(infer_map_in["tool_id"]), len(set(hmm_groups["cell"]) ^ set(infer_map_in["tool_id"])), 0)
    check("hmm_query_cells_complete", set(hmm_groups.loc[hmm_groups["cell"].str.startswith("Q"), "cell"]) == set(infer_map_in.loc[infer_map_in["role"] == "query", "tool_id"]), len(set(hmm_groups.loc[hmm_groups["cell"].str.startswith("Q"), "cell"]) ^ set(infer_map_in.loc[infer_map_in["role"] == "query", "tool_id"])), 0)
    check("hmm_regions_present", len(hmm_regions) > 0, len(hmm_regions), ">0")
    check("hmm_states_valid", set(hmm_regions["state"]).issubset(set(range(1, 7))), sorted(hmm_regions["state"].unique()), "subset 1..6")

    check("copykat_prediction_count", len(predictions) == 346, len(predictions), 346)
    check("copykat_predictions_tool_unique", not predictions["tool_id"].duplicated().any(), predictions["tool_id"].duplicated().sum(), 0)
    check("copykat_predictions_exact_mapped", predictions["exact_cell_id"].notna().all(), predictions["exact_cell_id"].isna().sum(), 0)
    check("copykat_all_input_cells_predicted", set(predictions["tool_id"]) == set(copy_map_in["tool_id"]), len(set(predictions["tool_id"]) ^ set(copy_map_in["tool_id"])), 0)
    reference_prediction = predictions[predictions["role"] == "reference"]
    check("copykat_reference_all_diploid", (reference_prediction["copykat.pred"] == "diploid").all(), int((reference_prediction["copykat.pred"] == "diploid").sum()), len(reference_prediction))
    check("copykat_query_prediction_nonmissing", predictions.loc[predictions["role"] == "query", "copykat.pred"].notna().all(), predictions.loc[predictions["role"] == "query", "copykat.pred"].isna().sum(), 0)

    gate_files = [
        pilot / "inputs/memory_gate_log_v1.tsv",
        infer_out / "memory_gate_log.tsv",
        copy_out / "memory_gate_log.tsv",
    ]
    gate_frames = [pd.read_csv(path, sep="\t") for path in gate_files]
    gate_text = "\n".join(path.read_text() for path in gate_files)
    check("memory_gate_files_present", all(path.is_file() for path in gate_files), sum(path.is_file() for path in gate_files), 3)
    check("memory_gates_all_pass", all("PASS" in line for line in gate_text.splitlines()[1:] if line.strip()), "all PASS" if all("PASS" in line for line in gate_text.splitlines()[1:] if line.strip()) else "not all PASS", "all PASS")
    check("memory_gate_record_count", sum(len(frame) for frame in gate_frames) == 10, sum(len(frame) for frame in gate_frames), 10)

    copy_summary = (
        predictions.groupby(["role", "copykat.pred"], dropna=False)
        .size()
        .rename("n_cells")
        .reset_index()
    )
    infer_summary = (
        hmm_groups.assign(role=hmm_groups["cell"].str[0].map({"Q": "query", "R": "reference"}))
        .groupby(["role", "cell_group_name"], dropna=False)
        .size()
        .rename("n_cells")
        .reset_index()
    )
    copy_summary.to_csv(verify_dir / "copykat_prediction_summary_v1.tsv", sep="\t", index=False)
    infer_summary.to_csv(verify_dir / "infercnv_group_summary_v1.tsv", sep="\t", index=False)

    checks_frame = pd.DataFrame(checks)
    checks_frame.to_csv(verify_dir / "verification_checks_v1.tsv", sep="\t", index=False)
    all_pass = bool(checks_frame["pass"].all())

    file_rows = []
    for path in sorted(pilot.rglob("*")):
        if not path.is_file() or verify_dir in path.parents:
            continue
        file_rows.append(
            {
                "relative_path": str(path.relative_to(pilot)),
                "size_bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
        )
    pd.DataFrame(file_rows).to_csv(verify_dir / "file_manifest_v1.tsv", sep="\t", index=False)

    metadata = {
        "schema_version": "part1_cnv_pilot_v2_verification_v1",
        "status": "completed" if all_pass else "failed",
        "all_checks_pass": all_pass,
        "n_checks": len(checks_frame),
        "n_checks_passed": int(checks_frame["pass"].sum()),
        "n_files_hashed": len(file_rows),
        "verified_at": datetime.now().astimezone().isoformat(),
        "pilot_dir": str(pilot),
        "copykat_query_prediction_counts": predictions.loc[predictions["role"] == "query", "copykat.pred"].value_counts().to_dict(),
        "infercnv_query_subclusters": int(hmm_groups.loc[hmm_groups["cell"].str.startswith("Q"), "cell_group_name"].nunique()),
        "warnings": [
            "pilot contains one patient/sample and is a technical validation, not a cohort biological conclusion",
            "timedatectl/system-bus and Matrix coercion deprecation warnings did not interrupt either tool",
        ],
    }
    (verify_dir / "run_metadata_v1.json").write_text(json.dumps(metadata, indent=2, ensure_ascii=False) + "\n")
    marker = verify_dir / ("RUN_COMPLETE" if all_pass else "RUN_FAILED")
    marker.write_text(datetime.now().astimezone().isoformat() + "\n")
    print(f"PART1_CNV_PILOT_V2_VERIFICATION_V1_{'PASS' if all_pass else 'FAIL'} checks={len(checks_frame)} files={len(file_rows)}")
    return 0 if all_pass else 1


if __name__ == "__main__":
    raise SystemExit(main())
