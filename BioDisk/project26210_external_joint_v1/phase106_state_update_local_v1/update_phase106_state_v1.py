from __future__ import annotations

import csv
import hashlib
import re
from datetime import datetime
from pathlib import Path


STAGE = Path(r"C:\Users\31937\project26210_phase106_state_update_v1")
ROOT = Path(r"D:\project26210_external_joint_v1")
AUDIT = ROOT / "phase106_dataset_separation_audit_local_v1"


STATUS = r"""

## Phase 106 dataset-separation and claim-routing correction (2026-08-31)

- Status: **PASS_PHASE106_DATASET_SEPARATION_VALIDATION 69/69**. This phase corrects the remaining cross-dataset conflation in the Phase 105 client blueprint. No biological analysis was rerun or modified.
- Phase 105 scan counts, artifact inventory and missing-registration findings remain valid. Its `CLIENT_RESPONSE_BLUEPRINT_v1.md` is superseded because it still linked original Part 1 single-cell, immunotherapy single-cell, spatial and public metabolite/cell-line results too closely.
- The project is now routed into 16 independent data/presentation universes, 44 phase routes, 27 dataset-specific claims and 12 non-comparability/conflict rules. A result cannot fill a missing edge, rescue a negative result or upgrade causality in another universe.
- Original Part 1 single-cell is D01 only: GSE156632, GSE159115, GSE178481, GSE207493 and GSE242299. Its F1/TYMP, macrophage subdivision, DPT/Monocle3, pySCENIC and broad CellChat conclusions remain inside D01 and its module-specific patient subsets.
- Immunotherapy single-cell is separate: D02 uses GSE264586/GSE294109 for clinical-context analyses; D03 uses the explicitly frozen GSE207493/GSE264586/GSE294109 three-cohort bridge. Overlap of GSE207493 does not make D01 and D03 results interchangeable because the cell subsets, cohort partners, endpoints and multiplicity families differ.
- Spatial results remain separated into D04 GLMF1/2, D05 TJ-RCC12 and D06 GSE175540. Sections/spots are not single-cell patients. GSE175540 supports the SPP1-state-to-vascular edge but not the upstream tumour-TYMP-to-neighbouring-SPP1 edge; the downstream positive result cannot complete the chain.
- Public measured-metabolite evidence is separated into D09 UnitedMet paired tumour RNA/MTA, D10 NCI60, D11 CCL180 and D12 CCLE. UnitedMet and NCI60 conclusions are not direct replication because they differ in system, endpoint and sample size. CCL180/CCLE adjacent metabolites are not MTA substitutes.
- Client unpublished validation is D14 and remains outside project analysis. Public observational baseline data cannot replace, predict or pre-empt the client's perturbation experiment or figure.
- For the current screenshot request, the client-facing target is D01 original Part 1 single-cell only. D02-D14 may be described later as separate modules or parallel evidence columns, never as pooled support for a single mechanism.
- Authoritative Phase 106 files are under `phase106_dataset_separation_audit_local_v1`. `CLIENT_RESPONSE_BLUEPRINT_v2.md` supersedes the Phase 105 v1 blueprint.
- Current breakpoint: reply to the client using the dataset-separated wording, then decide whether the requested deliverable is D01-only or a multi-module report with strict dataset boundaries. The unresolved D01 fine-state communication question must not be filled with D03 or spatial results.
"""


DECISIONS = r"""

## 2026-08-31 - Phase 106 dataset and conclusion separation decision

- Supersede the Phase 105 cross-dataset client blueprint with `phase106_dataset_separation_audit_local_v1/CLIENT_RESPONSE_BLUEPRINT_v2.md`. Preserve Phase 105 scan and registration-gap evidence.
- Do not combine original Part 1 single-cell, immunotherapy single-cell, GLMF, TJ-RCC, GSE175540, UnitedMet, NCI60, CCL180, CCLE, external bulk clinical or client unpublished validation into one statistical or causal evidence chain.
- Require every reported conclusion to name its dataset universe, biological unit, endpoint and evidence decision. A shared gene or verbal mechanism is not sufficient for direct comparison.
- Treat GSE207493 in D01 and D03 as overlapping accession use under different analysis universes. Do not transfer results without a frozen sample/cell intersection and new statistical contract.
- Keep GLMF1/2, TJ-RCC12 and GSE175540 as distinct spatial cohorts. Do not pool their spots or sections as one patient series.
- Keep UnitedMet tumour tissue and NCI60 cell-line associations separate. Retain NCI60 renal n=8 as exploratory direction only and do not call it validation or contradiction of UnitedMet.
- Record CCL180 and CCLE direct MTA as not evaluable. Adjacent metabolites and near-mass ions cannot rescue or replace direct MTA evidence.
- Keep the client's validation experiment and figure as a separate unavailable perturbational universe. Do not request, import, analyse, redraw, predict or replace it without explicit authorization.
- For the current client response, use D01 original Part 1 only. Other completed analyses may be mentioned as separate future modules, not as support rows for D01.
"""


HANDOFF = r"""

## Phase 106 dataset-separated current state (2026-08-31)

- Completion: Phase106 dataset-separation validation 69/69 PASS.
- Authoritative directory: `D:\project26210_external_joint_v1\phase106_dataset_separation_audit_local_v1`.
- Phase105 blueprint status: superseded for cross-dataset mixing; Phase105 scan counts and missing-registration audit remain valid.
- Current routing:
  - D01 original Part1 scRNA: GSE156632/GSE159115/GSE178481/GSE207493/GSE242299.
  - D02 immunotherapy clinical scRNA: GSE264586/GSE294109.
  - D03 three-cohort immune bridge: GSE207493/GSE264586/GSE294109 under a separate contract.
  - D04 GLMF spatial; D05 TJ-RCC spatial; D06 GSE175540 spatial.
  - D09 UnitedMet paired tumour RNA/MTA.
  - D10 NCI60; D11 CCL180; D12 CCLE.
  - D14 client unpublished validation; outside project analysis.
- Do not derive one unified TYMP-MTA-SPP1 conclusion across these universes. Use `DATASET_SPECIFIC_CLAIM_MATRIX_v1.tsv` and `NONCOMPARABILITY_AND_CONFLICT_REGISTER_v1.tsv`.
- Immediate client reply target: D01 original Part1 only. State that immune-treatment scRNA, spatial and NCI/cell-line analyses are separate and will not be mixed into this section.
- D01-specific unresolved item: internal five-study SPP1/C1QC/Other fine-state communication. Do not substitute D03 Phase88 or spatial results.
- Client validation data/figure remain client-owned; no project intervention.
- No Phase106 background job remains. Historical stopped Phase17 wrapper is unchanged and is not a live scientific task.
"""


def read(path: Path) -> tuple[str, bool]:
    raw = path.read_bytes()
    return raw.decode("utf-8-sig"), raw.startswith(b"\xef\xbb\xbf")


def write(path: Path, text: str, bom: bool) -> None:
    path.write_text(text, encoding="utf-8-sig" if bom else "utf-8", newline="")


def append_once(path: Path, marker: str, addendum: str) -> None:
    text, bom = read(path)
    count = text.count(marker)
    if count > 1:
        raise RuntimeError(f"duplicate marker: {marker}")
    if count == 0:
        newline = "\r\n" if "\r\n" in text else "\n"
        write(path, text.rstrip("\r\n") + addendum.replace("\n", newline) + newline, bom)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest().upper()


append_once(STAGE / "PROJECT_STATUS.md", "## Phase 106 dataset-separation and claim-routing correction (2026-08-31)", STATUS)
append_once(STAGE / "DECISIONS.md", "## 2026-08-31 - Phase 106 dataset and conclusion separation decision", DECISIONS)
append_once(STAGE / "CODEX_HANDOFF.md", "## Phase 106 dataset-separated current state (2026-08-31)", HANDOFF)

manifest_path = STAGE / "RESULTS_MANIFEST.tsv"
text, bom = read(manifest_path)
reader = csv.DictReader(text.splitlines(), delimiter="\t")
header = reader.fieldnames
rows = list(reader)
if header is None:
    raise RuntimeError("manifest header missing")
ids = {row["result_id"] for row in rows}
timestamp = datetime.now().astimezone().isoformat(timespec="seconds")


def kind(path: Path) -> str:
    return {".py": "formal_script", ".tsv": "routing_table_tsv", ".json": "validation_json", ".md": "audit_document_md"}.get(path.suffix, "audit_artifact")


for path in sorted(item for item in AUDIT.iterdir() if item.is_file()):
    stem = re.sub(r"[^a-z0-9]+", "_", path.stem.lower()).strip("_")
    result_id = f"phase106_dataset_separation_{stem}_v1"
    if result_id in ids:
        continue
    row = {
        "result_id": result_id,
        "stage": "phase106_dataset_separation_audit",
        "artifact_type": kind(path),
        "absolute_path": str(path),
        "generated_at": timestamp,
        "validation_status": "PASS",
        "validation_evidence": f"SHA256 {sha256(path)}; Phase106 dataset-separation validation 69/69 PASS",
        "notes": "Dataset-universe routing, dataset-specific claims and non-comparability controls",
    }
    values = [row[column] for column in header]
    text = text.rstrip("\r\n") + "\r\n" + "\t".join(values) + "\r\n"
    ids.add(result_id)
write(manifest_path, text, bom)
print("PHASE106_STATE_STAGING_UPDATE_PASS")
