#!/usr/bin/env python3
from __future__ import annotations

import csv
import hashlib
import json
import os
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parent
INPUT = ROOT / "inputs_v1"
RESULTS = ROOT / "results_v2"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


checks: list[dict[str, str]] = []


def check(name: str, condition: bool, detail: str) -> None:
    checks.append({"check": name, "status": "PASS" if condition else "FAIL", "detail": detail})


expected_inputs = {
    "IMMUNE_EVIDENCE_REGISTER_v1.tsv",
    "IMMUNE_EVIDENCE_GRADING_RULES_v1.tsv",
    "PHASE63_VALIDATION_v1.tsv",
    "large_cohort_convergence_grading_v2.tsv",
    "gse207493_cell_cycle_data_gate_v1.tsv",
    "gse207493_cell_cycle_adjusted_effect_summary_v1.tsv",
    "gse207493_metabolic_result_counts_v1.tsv",
    "gse207493_metabolic_contrast_summary_retained_v1.tsv",
    "gse207493_metabolic_correlation_summary_retained_v1.tsv",
    "gse207493_celltype_feature_summary_v1.tsv",
    "gse207493_author_composition_metrics_v1.tsv",
    "de_novo_tf_integrated_retained_v1.tsv",
    "de_novo_tf_cross_method_FDR_v1.tsv",
    "de_novo_integrated_status_v1.json",
    "de_novo_integrated_manifest_v1.tsv",
    "gse207493_de_novo_tf_regulon_record_v1.md",
    "de_novo_scenic_status_v1.json",
    "de_novo_collectri_viper_status_v1.json",
}
check("expected_inputs", expected_inputs.issubset({path.name for path in INPUT.iterdir()}), f"expected={len(expected_inputs)}")

phase63_validation = read_tsv(INPUT / "PHASE63_VALIDATION_v1.tsv")
check("phase63_no_fail", not any(row["status"] == "FAIL" for row in phase63_validation), str(Counter(row["status"] for row in phase63_validation)))

composition = {row["metric"]: float(row["value"]) for row in read_tsv(INPUT / "gse207493_author_composition_metrics_v1.tsv")}
check("gse207493_cells", composition["total_cells"] == 102723, str(composition["total_cells"]))
check("author_composition_concordance", composition["pearson_r_across_9_broad_types"] > 0.99, str(composition["pearson_r_across_9_broad_types"]))

cycle_gate = read_tsv(INPUT / "gse207493_cell_cycle_data_gate_v1.tsv")
check("cell_cycle_gates", len(cycle_gate) == 5 and all(row["status"] == "PASS" for row in cycle_gate), f"rows={len(cycle_gate)}")
cycle = read_tsv(INPUT / "gse207493_cell_cycle_adjusted_effect_summary_v1.tsv")
expected_cycle = {"F1_Treg_vs_otherT", "F1_Tumor_vs_nonTumor", "TYMP_MyeloidDC_vs_other", "TYMP_Treg_vs_otherT"}
check("cell_cycle_contrasts", {row["contrast"] for row in cycle} == expected_cycle, str(sorted(row["contrast"] for row in cycle)))
check("cell_cycle_all_fdr", all(float(row["adjusted_q_value"]) < 0.05 for row in cycle), "4/4 q<0.05")

met = {row["contrast"]: row for row in read_tsv(INPUT / "gse207493_metabolic_result_counts_v1.tsv")}
check("metabolic_contrasts", set(met) == {"F1_high_vs_low_Tumor", "TYMP_high_vs_low_MyeloidDC", "TYMP_high_vs_low_Treg"}, str(sorted(met)))
check("metabolic_fdr_counts", [int(met[key]["adjusted_q_lt_0.05"]) for key in ["F1_high_vs_low_Tumor", "TYMP_high_vs_low_MyeloidDC", "TYMP_high_vs_low_Treg"]] == [22, 17, 10], "22,17,10")

with (INPUT / "de_novo_integrated_status_v1.json").open(encoding="utf-8") as handle:
    integrated = json.load(handle)
with (INPUT / "de_novo_scenic_status_v1.json").open(encoding="utf-8") as handle:
    scenic = json.load(handle)
with (INPUT / "de_novo_collectri_viper_status_v1.json").open(encoding="utf-8") as handle:
    viper = json.load(handle)
check("de_novo_status", integrated["status"] == "PASS_DE_NOVO_TF_REGULON_PRIMARY_V1", integrated["status"])
check("de_novo_no_old_primary", integrated["old_part1_primary_used"] is False and scenic["old_part1_regulons_used"] is False and viper["old_part1_network_used"] is False, "all_false")
check("de_novo_dimensions", integrated["n_cells"] == 39912 and integrated["n_samples"] == 19 and integrated["n_integrated_rows"] == 4992 and integrated["n_retained_rows"] == 2826, json.dumps(integrated, sort_keys=True))

tf_cross = read_tsv(INPUT / "de_novo_tf_cross_method_FDR_v1.tsv")
counts = Counter(row["endpoint"] for row in tf_cross)
check("cross_method_total", len(tf_cross) == 176, str(len(tf_cross)))
check("cross_method_endpoint_counts", counts == Counter({"tumor_f1": 95, "myeloid_dc_tymp": 67, "treg_tymp": 14}), str(dict(counts)))
check("cross_method_grade", all(row["evidence_grade"] == "A_cross_method_FDR" for row in tf_cross), "all A_cross_method_FDR")
check("cross_method_fdr", all(float(row["q_SCENIC"]) < 0.05 and float(row["q_VIPER"]) < 0.05 for row in tf_cross), "both methods q<0.05")
check("cross_method_direction", all(row["same_direction"].lower() == "true" for row in tf_cross), "all same_direction")

remote_manifest = {row["file"]: row for row in read_tsv(INPUT / "de_novo_integrated_manifest_v1.tsv")}
manifest_mapping = {
    "de_novo_tf_cross_method_FDR_v1.tsv": "de_novo_tf_cross_method_FDR_v1.tsv",
    "de_novo_tf_integrated_retained_v1.tsv": "de_novo_tf_integrated_retained_v1.tsv",
    "gse207493_de_novo_tf_regulon_record_v1.md": "gse207493_de_novo_tf_regulon_record_v1.md",
    "status_v1.json": "de_novo_integrated_status_v1.json",
}
for remote_name, local_name in manifest_mapping.items():
    path = INPUT / local_name
    expected = remote_manifest[remote_name]
    check(f"remote_hash_{remote_name}", path.stat().st_size == int(expected["size_bytes"]) and sha256(path) == expected["sha256"], local_name)

evidence = read_tsv(RESULTS / "IMMUNE_MECHANISTIC_EVIDENCE_REGISTER_v2.tsv")
check("phase65_evidence_rows", len(evidence) == 15, str(len(evidence)))
check("phase65_evidence_ids", len({row["evidence_id"] for row in evidence}) == 15, "unique")
check("phase65_grade_counts", Counter(row["global_grade"] for row in evidence) == Counter({"A": 3, "B": 4, "C1": 6, "D": 2}), str(dict(Counter(row["global_grade"] for row in evidence))))
direct = [row for row in evidence if row["domain"] == "direct_clinical"]
check("direct_predictors_not_promoted", len(direct) == 2 and all(row["global_grade"] == "D" for row in direct), "F1/TYMP both D")

tf_summary = read_tsv(RESULTS / "TF_CROSS_METHOD_ENDPOINT_SUMMARY_v1.tsv")
check("tf_summary_rows", len(tf_summary) == 3, str(len(tf_summary)))
check("tf_summary_closure", sum(int(row["cross_method_fdr_rows"]) for row in tf_summary) == 176, str(sum(int(row["cross_method_fdr_rows"]) for row in tf_summary)))
frozen_lists = {}
for line in (INPUT / "gse207493_de_novo_tf_regulon_record_v1.md").read_text(encoding="utf-8").splitlines():
    for endpoint in ("tumor_f1", "myeloid_dc_tymp", "treg_tymp"):
        prefix = f"- {endpoint}: "
        if line.startswith(prefix):
            frozen_lists[endpoint] = ",".join(value.strip() for value in line[len(prefix):].split(","))
check("tf_frozen_display_lists", all(row["top_15_unique_tfs"] == frozen_lists[row["endpoint"]] for row in tf_summary), "all endpoints match frozen Phase64 record")
met_summary = read_tsv(RESULTS / "METABOLIC_ENDPOINT_SUMMARY_v1.tsv")
check("metabolic_summary_rows", len(met_summary) == 3, str(len(met_summary)))

summary_text = (RESULTS / "PART1_IMMUNE_SERIES_CONVERGENCE_v2.md").read_text(encoding="utf-8")
check("summary_boundary", "仍未验证 F1 或 TYMP 是独立的免疫治疗疗效标志物" in summary_text and "不能写" in summary_text, "clinical boundary present")
check("supersession_record", (RESULTS / "SUPERSESSION_RECORD_v2.md").exists(), "results_v2 supersession documented")

validation_path = RESULTS / "PHASE65_VALIDATION_v1.tsv"
with validation_path.open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=["check", "status", "detail"], delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(checks)

failures = [row for row in checks if row["status"] != "PASS"]
manifest_rows = []
for path in sorted(RESULTS.iterdir()):
    if path.is_file() and path.name not in {"OUTPUT_MANIFEST_SHA256_v1.tsv", "PHASE65_PASS_v1"}:
        manifest_rows.append({"file": path.name, "size_bytes": path.stat().st_size, "sha256": sha256(path)})
with (RESULTS / "OUTPUT_MANIFEST_SHA256_v1.tsv").open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=["file", "size_bytes", "sha256"], delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(manifest_rows)

if failures:
    raise SystemExit("Phase65 validation failed: " + ", ".join(row["check"] for row in failures))
(RESULTS / "PHASE65_PASS_v1").write_text(f"status=PASS\nchecks={len(checks)}\n", encoding="utf-8")
print(json.dumps({"status": "PASS", "checks": len(checks), "files_hashed": len(manifest_rows)}, indent=2))
