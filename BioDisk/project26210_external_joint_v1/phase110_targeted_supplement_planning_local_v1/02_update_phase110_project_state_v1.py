from __future__ import annotations

import csv
import hashlib
import json
from datetime import datetime
from pathlib import Path


ROOT = Path(r"D:\project26210_external_joint_v1")
PHASE = ROOT / "phase110_targeted_supplement_planning_local_v1"
NOW = datetime.now().astimezone().isoformat(timespec="seconds")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def append_once(path: Path, marker: str, section: str) -> None:
    current = path.read_text(encoding="utf-8-sig")
    if marker in current:
        raise RuntimeError(f"Duplicate marker refused: {marker}")
    with path.open("a", encoding="utf-8", newline="") as handle:
        if current and not current.endswith("\n"):
            handle.write("\n")
        handle.write("\n" + section.strip() + "\n")


validation = json.loads((PHASE / "PHASE110_PLANNING_VALIDATION_v1.json").read_text(encoding="utf-8"))
if validation.get("status") != "PASS_PHASE110_TARGETED_SUPPLEMENT_PLAN_VALIDATION" or validation.get("checks_passed") != 33:
    raise RuntimeError("Phase110 planning validation is not 33/33 PASS")

status_section = """
## Phase 110 targeted supplementation planning and 24-section spatial reuse audit (2026-09-01)

- Status: **PASS_PHASE110_TARGETED_SUPPLEMENT_PLAN_VALIDATION 33/33; PLANNING ONLY; NO NEW BIOLOGICAL ANALYSIS EXECUTED**.
- Phase109 is retained as a valid completed local gap test but is no longer treated as the whole targeted-supplementation task. The corrected task is a project-wide asset-to-gap-to-algorithm/data plan.
- The GSE175540 24-section asset is intact on `sgs10-codex`: 24 Phase54 processed H5AD files, 24 Phase78 fine-state abundance tables, 24 posterior H5AD files, 12 FFPE and 12 frozen sections, 72,952 QC-passing spots, tissue images for all sections and official TLS annotations for 23 sections.
- Prior GSE175540 use was substantial but narrow: preprocessing, fixed signatures, Moran's I, simple TLS comparisons, coarse/fine cell2location, targeted pairwise spatial relations, an MTA-response proxy and one tumour-TYMP-to-neighbour-SPP1 test. Cross-section fine-state niche discovery, spatially gated ligand-receptor validation, TLS boundary gradients and niche-specific transcriptomics remain unperformed.
- The first planned execution block is M110A-D: immutable input/metadata audit; FFPE-discovered/frozen-transferred fine-state niche atlas; a maximum ten-pair SPP1-endothelial spatial ligand-receptor test selected independently from Phase35; and official-TLS-distance immune gradients.
- M110E niche-specific paired section pseudobulk is next. Conditional TYMP heterogeneity, histology embeddings and new dataset acquisition are secondary or conditional.
- No rerun of cell2location is planned unless input validation invalidates an existing posterior. Sections, not spots, remain the inferential units; FFPE and frozen remain separate preparation families.
- Direct TYMP-to-measured-MTA and MTA-to-SPP1 induction cannot be closed by added spatial algorithms. Public perturbation data must pass a strict direct-exposure/biological-replicate/SPP1-readout gate; otherwise these edges remain client-owned validation tasks.
- Authoritative plan: `phase110_targeted_supplement_planning_local_v1/PHASE110_TARGETED_SUPPLEMENT_PLAN_v1.md`.
"""

decision_section = """
## Phase 110 targeted-supplementation planning decisions (2026-09-01)

- Correct the scope: targeted supplementation starts with a complete asset and result inventory, not immediate execution of another isolated association analysis.
- Retain Phase109 as completed evidence, but use Phase110 as the current planning breakpoint.
- Prioritize the existing GSE175540 24-section data before downloading another broad dataset.
- Approve M110A-D as the recommended first executable block after explicit execution approval and frozen preregistration; this planning phase itself does not authorize analysis.
- Learn spatial niches without TYMP, TLS or clinical outcomes, then transfer FFPE-derived centroids unchanged to frozen sections. Reverse transfer is sensitivity only.
- Cap the primary spatial ligand-receptor family at five SPP1-to-endothelial and five endothelial-to-SPP1 pairs, selected from the independent Phase35 patient table before spatial outcome inspection.
- Use official TLS geometry for distance and boundary tests. Do not replace it with a transcript-derived TLS label in the primary analysis.
- Do not rerun cell2location, reopen TYMP-MTA baseline correlations, pool FFPE with frozen, treat spots as biological replicates or add datasets solely to seek a preferred direction.
- Add GSE210038/041, PRJNA705464, SCP1288 or GSE222703 only for a fixed phenotype and a declared validation role after M110B-E.
"""

handoff_section = """
## Phase 110 current planning breakpoint (2026-09-01)

- Completion: planning validation 33/33 PASS; no new biological analysis executed.
- Authoritative directory: `D:\\project26210_external_joint_v1\\phase110_targeted_supplement_planning_local_v1`.
- Read first: `PHASE110_TARGETED_SUPPLEMENT_PLAN_v1.md`, `CHAIN_GAP_ACTION_MATRIX_v1.tsv`, `CURRENT_ASSET_UTILIZATION_MATRIX_v1.tsv` and `GSE175540_24_SECTION_REUSE_PLAN_v1.md`.
- Correct interpretation of Phase109: a valid completed local gap test, not the complete targeted-supplementation programme.
- Highest-value available data: GSE175540 24 sections with 24/24 processed H5AD, 24/24 fine-state q05/mean posteriors, all tissue images and 23 official TLS maps. Existing Cell2location outputs should be reused.
- Recommended first execution block after approval: M110A input/metadata audit, M110B fine-state niche transfer, M110C capped LR spatial validation and M110D TLS-distance immune gradients. M110E paired section pseudobulk follows only after niche labels are frozen.
- Causal upstream gaps remain outside computational closure: direct TYMP-to-MTA and MTA-to-SPP1 require eligible perturbation data or client experiments.
- Do not start M110F-G or acquire new cohorts before M110A-D results define a fixed target.
"""

append_once(ROOT / "PROJECT_STATUS.md", "## Phase 110 targeted supplementation planning", status_section)
append_once(ROOT / "DECISIONS.md", "## Phase 110 targeted-supplementation planning decisions", decision_section)
append_once(ROOT / "CODEX_HANDOFF.md", "## Phase 110 current planning breakpoint", handoff_section)

manifest = ROOT / "RESULTS_MANIFEST.tsv"
with manifest.open("r", encoding="utf-8-sig", newline="") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    fields = reader.fieldnames
    rows = list(reader)
expected = ["result_id", "stage", "artifact_type", "absolute_path", "generated_at", "validation_status", "validation_evidence", "notes"]
if fields != expected:
    raise RuntimeError(f"Unexpected manifest header: {fields}")

artifacts = [
    ("phase110_targeted_supplement_plan_v1", "planning_document_md", PHASE / "PHASE110_TARGETED_SUPPLEMENT_PLAN_v1.md"),
    ("phase110_asset_utilization_matrix_v1", "planning_table_tsv", PHASE / "CURRENT_ASSET_UTILIZATION_MATRIX_v1.tsv"),
    ("phase110_chain_gap_action_matrix_v1", "planning_table_tsv", PHASE / "CHAIN_GAP_ACTION_MATRIX_v1.tsv"),
    ("phase110_module_registry_v1", "planning_table_tsv", PHASE / "PLANNED_MODULE_REGISTRY_v1.tsv"),
    ("phase110_gse175540_reuse_plan_v1", "planning_document_md", PHASE / "GSE175540_24_SECTION_REUSE_PLAN_v1.md"),
    ("phase110_source_registry_v1", "source_registry_tsv", PHASE / "SOURCE_REGISTRY_v1.tsv"),
    ("phase110_plan_validator_v1", "formal_validation_script", PHASE / "01_validate_phase110_plan_v1.py"),
    ("phase110_planning_validation_v1", "validation_json", PHASE / "PHASE110_PLANNING_VALIDATION_v1.json"),
    ("phase110_file_manifest_v1", "file_manifest_tsv", PHASE / "FILE_MANIFEST_SHA256_v1.tsv"),
    ("phase110_state_update_script_v1", "formal_script", PHASE / "02_update_phase110_project_state_v1.py"),
]
existing = {row["result_id"] for row in rows}
new_rows = []
for result_id, artifact_type, path in artifacts:
    if result_id in existing:
        raise RuntimeError(f"Duplicate manifest ID refused: {result_id}")
    if not path.is_file():
        raise FileNotFoundError(path)
    new_rows.append({
        "result_id": result_id,
        "stage": "phase110_targeted_supplement_planning",
        "artifact_type": artifact_type,
        "absolute_path": str(path),
        "generated_at": NOW,
        "validation_status": "PASS",
        "validation_evidence": f"SHA256 {sha256(path).upper()}; Phase110 plan validation 33/33 PASS; no analysis executed",
        "notes": "Project-wide gap inventory and targeted algorithm/dataset supplementation plan",
    })
with manifest.open("a", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writerows(new_rows)

state = {
    "status": "PASS_PHASE110_STATE_UPDATE",
    "generated_at": NOW,
    "manifest_rows_added": len(new_rows),
    "planning_validation": "33/33",
    "analysis_executed": False,
}
(PHASE / "PHASE110_STATE_UPDATE_STATUS_v1.json").write_text(json.dumps(state, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
print(json.dumps(state, ensure_ascii=False, indent=2))
