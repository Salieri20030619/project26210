from __future__ import annotations

import csv
import hashlib
import json
from datetime import datetime
from pathlib import Path


ROOT = Path(r"D:\project26210_external_joint_v1")
PHASE = ROOT / "phase109_targeted_gap_analysis_local_v1"
NOW = datetime.now().astimezone().isoformat(timespec="seconds")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def append_once(path: Path, marker: str, text: str) -> None:
    current = path.read_text(encoding="utf-8-sig")
    if marker in current:
        raise RuntimeError(f"Refusing duplicate state section {marker} in {path}")
    with path.open("a", encoding="utf-8", newline="") as handle:
        if current and not current.endswith("\n"):
            handle.write("\n")
        handle.write("\n" + text.strip() + "\n")


required = [
    PHASE / "PHASE109_PREREGISTRATION_v1.md",
    PHASE / "analysis_v3" / "PHASE109_ANALYSIS_STATUS_v1.json",
    PHASE / "analysis_v3" / "abundance_primary_results_v1.tsv",
    PHASE / "analysis_v3" / "cellchat_fine_state_specificity_summary_v1.tsv",
    PHASE / "independent_validation_v1" / "PHASE109_INDEPENDENT_VALIDATION_STATUS_v1.json",
    PHASE / "final_synthesis_v1" / "PHASE109_CHAIN_EDGE_UPDATE_v1.tsv",
    PHASE / "final_synthesis_v1" / "PHASE109_TARGETED_SUPPLEMENT_DECISION_v1.md",
]
missing = [str(p) for p in required if not p.is_file()]
if missing:
    raise FileNotFoundError("Missing Phase109 required files: " + "; ".join(missing))

analysis_status = json.loads((PHASE / "analysis_v3" / "PHASE109_ANALYSIS_STATUS_v1.json").read_text(encoding="utf-8"))
validation_status = json.loads((PHASE / "independent_validation_v1" / "PHASE109_INDEPENDENT_VALIDATION_STATUS_v1.json").read_text(encoding="utf-8"))
if analysis_status["status"] != "PASS_ANALYSIS_COMPLETE":
    raise RuntimeError("Formal Phase109 analysis is not PASS")
if validation_status["status"] != "PASS_PHASE109_INDEPENDENT_VALIDATION" or validation_status["checks_passed"] != 54:
    raise RuntimeError("Independent Phase109 validation is not 54/54 PASS")

status_section = """
## Phase 109 targeted supplementation of unresolved chain edges (2026-08-31)

- Status: **PASS_ANALYSIS_COMPLETE; PASS_PHASE109_INDEPENDENT_VALIDATION 54/54**. The analysis was preregistered before result inspection and retained 100,000 within-study permutations with seed 262109.
- Original Part1 patient-level bridge: 20 patients met the frozen eligibility rule. Zero of four preregistered tumour-axis-to-SPP1 abundance/polarization endpoints were supported. TYMP-axis estimates were negative; TYMP RNA estimates were weak/inconsistent; all BH q values were at least 0.67747. Broad-axis and F1 sensitivities did not pass the frozen gate.
- Internal fine-state communication: the existing Phase35 patient interaction table allowed exact reaggregation into SPP1/C1QC/Other macrophage states without rerunning CellChat. Four endpoints were tested across 14 patients. Malignant-to-SPP1 and SPP1-to-malignant specificity were not supported.
- SPP1-to-endothelial specificity was supported in 14/14 patients (median 0.21065; BH q=0.00305), and endothelial-to-SPP1 specificity was supported in 13/14 patients (median 0.14642; BH q=0.00305). These are computational communication-compatibility results, not ligand binding, induction or causality.
- Frozen Phase70 synthesis retains GSE175540 SPP1-state-to-vascular support in FFPE and frozen preparations. It triangulates the downstream interface but is not pooled with the CellChat patients.
- Direct TYMP-to-measured-MTA and MTA-to-SPP1 induction were not tested in Phase109. They remain unsupported/pending client perturbational validation; no client data or figure was imported or anticipated.
- Two pre-formal implementation failures are preserved under `analysis_v1` and `analysis_v2`. They occurred after numerical computation during support-label construction and are not accepted result sets. `analysis_v3` is the sole formal output.
- Authoritative package: `phase109_targeted_gap_analysis_local_v1`. Current breakpoint: update the unified report around the supported downstream SPP1/endothelial-vascular interface while explicitly showing the upstream TYMP/MTA/SPP1 breaks; do not continue same-family endpoint expansion to seek a preferred upstream result.
"""

decision_section = """
## Phase 109 targeted-gap decisions (2026-08-31)

- Accept the frozen Phase109 primary result: current original Part1 data do not support a general tumour TYMP/F1-to-SPP1-state abundance, polarization or malignant-SPP1 communication edge.
- Accept bidirectional SPP1-state/endothelial fine-state communication as supported computational compatibility: SPP1-to-endothelial 14/14 positive and endothelial-to-SPP1 13/14 positive, both BH q=0.00305.
- Triangulate, but do not pool, the fine-state CellChat result with the previously validated GSE175540 SPP1-state-to-vascular spatial result.
- Do not describe CellChat as direct ligand binding, protein validation, induction or causality.
- Keep direct TYMP-to-measured-MTA as current-data-do-not-support and MTA-to-SPP1 as pending-client-validation. Client perturbational data and figures remain outside project analysis.
- Do not add post-result endpoints or relax thresholds to rescue the upstream chain. The next report revision must preserve unsupported edges in their biological positions and use `PHASE109_CHAIN_EDGE_UPDATE_v1.tsv` as the current edge ledger.
"""

handoff_section = """
## Phase 109 current breakpoint (2026-08-31)

- Completion: formal Phase109 analysis PASS; independent recomputation 54/54 PASS.
- Authoritative directory: `D:\\project26210_external_joint_v1\\phase109_targeted_gap_analysis_local_v1`.
- Formal result directory: `analysis_v3` only. `analysis_v1` and `analysis_v2` are preserved implementation-failure audits and must not be cited as results.
- Primary upstream result: 0/4 patient-level TYMP/F1-to-SPP1 abundance/polarization endpoints supported in 20 eligible original Part1 patients; 0/2 malignant-SPP1 fine-state communication endpoints supported in 14 CellChat patients.
- New supported downstream result: SPP1-to-endothelial specificity 14/14 positive, median 0.21065, q=0.00305; endothelial-to-SPP1 13/14 positive, median 0.14642, q=0.00305.
- Frozen spatial result: GSE175540 SPP1-state-to-vascular remodeling remains supported in both FFPE and frozen preparations; TYMP-to-SPP1 remains unsupported spatially.
- Use `final_synthesis_v1/PHASE109_CHAIN_EDGE_UPDATE_v1.tsv` and `final_synthesis_v1/PHASE109_TARGETED_SUPPLEMENT_DECISION_v1.md` for the next unified report revision.
- Direct measured MTA and causal MTA-to-SPP1 remain outside the available data and pending client validation. Do not request, import, predict or redraw client-owned validation outputs unless explicitly asked.
- Immediate next action: reconstruct the unified client Results around component evidence, placing the supported SPP1/endothelial-vascular interface downstream of, but not as a rescue for, the unsupported TYMP/MTA/SPP1 links.
"""

append_once(ROOT / "PROJECT_STATUS.md", "## Phase 109 targeted supplementation", status_section)
append_once(ROOT / "DECISIONS.md", "## Phase 109 targeted-gap decisions", decision_section)
append_once(ROOT / "CODEX_HANDOFF.md", "## Phase 109 current breakpoint", handoff_section)

# Extend the package hash inventory with the state updater itself.
package_manifest = PHASE / "FILE_MANIFEST_SHA256_v1.tsv"
with package_manifest.open("r", encoding="utf-8-sig", newline="") as handle:
    package_rows = list(csv.DictReader(handle, delimiter="\t"))
script_path = PHASE / "08_update_phase109_project_state_v1.py"
script_rel = script_path.relative_to(PHASE).as_posix()
if not any(row["relative_path"] == script_rel for row in package_rows):
    package_rows.append({"relative_path": script_rel, "size_bytes": str(script_path.stat().st_size), "sha256": sha256(script_path)})
with package_manifest.open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=["relative_path", "size_bytes", "sha256"], delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(sorted(package_rows, key=lambda row: row["relative_path"]))

manifest = ROOT / "RESULTS_MANIFEST.tsv"
with manifest.open("r", encoding="utf-8-sig", newline="") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    fieldnames = reader.fieldnames
    rows = list(reader)
if fieldnames != ["result_id", "stage", "artifact_type", "absolute_path", "generated_at", "validation_status", "validation_evidence", "notes"]:
    raise RuntimeError(f"Unexpected RESULTS_MANIFEST header: {fieldnames}")

artifacts = [
    ("phase109_preregistration_v1", "preregistration_md", PHASE / "PHASE109_PREREGISTRATION_v1.md"),
    ("phase109_input_audit_status_v1", "input_audit_json", PHASE / "input_audit_remote_sync_v1" / "PHASE109_INPUT_AUDIT_STATUS_v1.json"),
    ("phase109_feasibility_status_v1", "feasibility_json", PHASE / "feasibility_v1" / "PHASE109_FEASIBILITY_STATUS_v1.json"),
    ("phase109_analysis_script_v1", "formal_script", PHASE / "02_run_phase109_targeted_analysis_v1.py"),
    ("phase109_analysis_status_v1", "analysis_status_json", PHASE / "analysis_v3" / "PHASE109_ANALYSIS_STATUS_v1.json"),
    ("phase109_abundance_primary_v1", "evidence_table_tsv", PHASE / "analysis_v3" / "abundance_primary_results_v1.tsv"),
    ("phase109_abundance_sensitivity_v1", "evidence_table_tsv", PHASE / "analysis_v3" / "abundance_exposure_sensitivity_v1.tsv"),
    ("phase109_cellchat_patient_specificity_v1", "evidence_table_tsv", PHASE / "analysis_v3" / "cellchat_patient_fine_state_specificity_v1.tsv"),
    ("phase109_cellchat_specificity_summary_v1", "evidence_table_tsv", PHASE / "analysis_v3" / "cellchat_fine_state_specificity_summary_v1.tsv"),
    ("phase109_spatial_frozen_synthesis_v1", "evidence_table_tsv", PHASE / "analysis_v3" / "spatial_frozen_edge_synthesis_v1.tsv"),
    ("phase109_analysis_summary_v1", "analysis_document_md", PHASE / "analysis_v3" / "PHASE109_TARGETED_ANALYSIS_SUMMARY_v1.md"),
    ("phase109_independent_validator_v1", "formal_validation_script", PHASE / "06_validate_phase109_targeted_v1.py"),
    ("phase109_independent_validation_v1", "validation_json", PHASE / "independent_validation_v1" / "PHASE109_INDEPENDENT_VALIDATION_STATUS_v1.json"),
    ("phase109_failure_audit_analysis_v1", "failure_audit_md", PHASE / "analysis_v1" / "FAILURE_AUDIT_v1.md"),
    ("phase109_failure_audit_analysis_v2", "failure_audit_md", PHASE / "analysis_v2" / "FAILURE_AUDIT_v2.md"),
    ("phase109_chain_edge_update_v1", "evidence_table_tsv", PHASE / "final_synthesis_v1" / "PHASE109_CHAIN_EDGE_UPDATE_v1.tsv"),
    ("phase109_targeted_supplement_decision_v1", "analysis_document_md", PHASE / "final_synthesis_v1" / "PHASE109_TARGETED_SUPPLEMENT_DECISION_v1.md"),
    ("phase109_targeted_gap_summary_png_v1", "figure_png", PHASE / "final_synthesis_v1" / "PHASE109_TARGETED_GAP_SUMMARY_v1.png"),
    ("phase109_targeted_gap_summary_pdf_v1", "figure_pdf", PHASE / "final_synthesis_v1" / "PHASE109_TARGETED_GAP_SUMMARY_v1.pdf"),
    ("phase109_file_manifest_v1", "file_manifest_tsv", package_manifest),
    ("phase109_state_update_script_v1", "formal_script", script_path),
]
existing_ids = {row["result_id"] for row in rows}
new_rows = []
for result_id, artifact_type, path in artifacts:
    if result_id in existing_ids:
        raise RuntimeError(f"Duplicate result_id refused: {result_id}")
    if not path.is_file():
        raise FileNotFoundError(path)
    new_rows.append(
        {
            "result_id": result_id,
            "stage": "phase109_targeted_gap_analysis",
            "artifact_type": artifact_type,
            "absolute_path": str(path),
            "generated_at": NOW,
            "validation_status": "PASS",
            "validation_evidence": f"SHA256 {sha256(path).upper()}; Phase109 formal analysis and independent validation 54/54 PASS",
            "notes": "Targeted unresolved-edge supplementation under the unified chain and dataset-separation rules",
        }
    )
with manifest.open("a", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
    writer.writerows(new_rows)

completion = {
    "status": "PASS_PHASE109_STATE_UPDATE",
    "generated_at": NOW,
    "manifest_rows_added": len(new_rows),
    "analysis_status": analysis_status["status"],
    "independent_validation": f"{validation_status['checks_passed']}/{validation_status['checks_total']}",
    "four_file_markers": {
        "PROJECT_STATUS.md": "## Phase 109 targeted supplementation",
        "DECISIONS.md": "## Phase 109 targeted-gap decisions",
        "RESULTS_MANIFEST.tsv": "phase109_chain_edge_update_v1",
        "CODEX_HANDOFF.md": "## Phase 109 current breakpoint",
    },
}
(PHASE / "PHASE109_STATE_UPDATE_STATUS_v1.json").write_text(json.dumps(completion, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
print(json.dumps(completion, ensure_ascii=False, indent=2))
