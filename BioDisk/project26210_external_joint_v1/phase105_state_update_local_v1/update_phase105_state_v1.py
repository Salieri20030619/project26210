from __future__ import annotations

import csv
import hashlib
import re
from datetime import datetime
from pathlib import Path


STAGE = Path(r"C:\Users\31937\project26210_phase105_state_update_v1")
ROOT = Path(r"D:\project26210_external_joint_v1")
AUDIT = ROOT / "phase105_full_analysis_audit_local_v1"


STATUS = r"""

## Phase 105 full-project analysis audit and client-scope re-lock (2026-08-31)

- Status: **PASS_FULL_PROJECT_SCAN; PASS_PHASE105_AUDIT_VALIDATION 56/56**. The audit scanned the authoritative D-drive project and `sgs10-codex` recovery workspace without rerunning or modifying any biological analysis.
- Coverage: 574 manifest rows with 574 unique IDs and 103 stage labels; 91 local top-level phase directories; 81 remote `26210` directories; zero missing registered local paths; 7/7 required remote evidence files present.
- State-record gaps were confirmed: 16 top-level phase directories have no registered artifact path, 15 numeric phases have a local directory without a manifest stage, and 12 phases have a directory without a dedicated status heading. Phase 21, Phase 57 and Phase 58 are completed analyses directly relevant to the current client request but were not previously registered in the long-term manifest.
- The client scope is now frozen to the internal five-study single-cell cohort: GSE156632, GSE159115, GSE178481, GSE207493 and GSE242299. External validation cohorts, spatial/bulk/public-metabolomics/cell-line modules and the client's unpublished validation experiments do not drive the immediate client-facing main line.
- The required macrophage-first SPP1 chain is already substantially complete. Phase 8 subdivided 71,045 macrophage-pool cells into Other Macrophage 40,429, SPP1-state 15,322 and C1QC-state 15,294. Phase 21 and Phase 33 used 50,761 myeloid cells from 40 patients in five study-separated graphs and showed consistent cross-sectional Monocyte-to-SPP1-state ordering. Phase 34 completed internal five-study myeloid pySCENIC.
- A key internal negative is preserved: the SPP1-versus-C1QC contrast did not increase along Phase 21 DPT in any of the five studies. SPP1-state ordering therefore does not imply uniform replacement of the C1QC programme.
- The remaining analysis gap is explicit: Phase 35 internal CellChat used a broad macrophage grouping, whereas Phase 88 fine-state CellChat used later external cohorts. No formal internal-five-study SPP1/C1QC/Other fine-state-specific communication analysis was found. Broad Phase 35 communication and external Phase 88 cannot be relabelled as internal SPP1-state-specific communication.
- Existing Phase 22/23, 40–46 and 89–100 reports are not acceptable bases for incremental revision under the new request. The next client-facing content must be redesigned from the internal evidence ledger in the order atlas -> myeloid extraction -> macrophage subdivision -> SPP1-state definition -> patient distribution -> DPT/Monocle3 ordering -> myeloid pySCENIC -> communication/function.
- The client's own validation experiments and future figures remain outside project analysis. No client data were requested, imported, reanalysed, redrawn or pre-emptively interpreted in Phase 105.
- Authoritative Phase 105 outputs are under `phase105_full_analysis_audit_local_v1`. Current breakpoint: decide whether to run the identified internal fine-state communication gap before building the redesigned client content.
"""


DECISIONS = r"""

## 2026-08-31 - Phase 105 internal-single-cell redesign decision

- Use only GSE156632, GSE159115, GSE178481, GSE207493 and GSE242299 as the main evidence scope for the immediate client-facing single-cell redesign. Keep studies separate and patients as the inferential units.
- Introduce SPP1 only after explicit myeloid extraction and macrophage subdivision into SPP1-state, C1QC-state and Other Macrophage. Do not define an SPP1 subgroup directly from the unsplit immune atlas.
- Use Phase 8, Phase 21, Phase 33 and Phase 34 as the validated internal macrophage/SPP1 main line. Use Phase 57–58 to disclose internal cohort heterogeneity and GSE207493 weighting.
- Preserve the internal negative result that SPP1-versus-C1QC does not uniformly increase along the state ordering. Do not describe the trajectory as C1QC replacement, lineage conversion, RNA velocity or causal polarization.
- Treat Phase 35 as broad macrophage communication only. Do not claim selective SPP1-state communication until a formal internal-five-study fine-state-specific analysis is completed.
- Keep Phase 60–104 external, spatial, bulk, perturbation, public-metabolomics and cell-line results in the internal audit layer. They may constrain claims but must not structure the immediate client response unless the client later asks for them.
- Do not use the Phase 91–100 reports or Phase 95–97 candidate figures as the new report base. Rebuild the content architecture from the Phase 105 evidence ledger.
- Keep the client's unpublished validation experiment and figure outside the analysis workspace. Do not request, import, reanalyse, redraw or infer those results without new explicit authorization.
"""


HANDOFF = r"""

## Phase 105 full-project audit and current client task (2026-08-31)

- Completion: `PASS_FULL_PROJECT_SCAN`; independent audit validation 56/56 PASS.
- Authoritative audit directory: `D:\project26210_external_joint_v1\phase105_full_analysis_audit_local_v1`.
- Frozen client requirement: internal single-cell cohorts only; macrophage subdivision before SPP1; external cohorts do not drive the main line; client validation data/figures remain client-owned and outside analysis.
- Internal cohort: GSE156632, GSE159115, GSE178481, GSE207493 and GSE242299.
- Reusable main-line evidence:
  - Phase 8 macrophage subdivision: 71,045 cells; Other 40,429; SPP1-state 15,322; C1QC-state 15,294.
  - Phase 21 DPT and Phase 33 Monocle3: 50,761 myeloid cells, 40 patients, five separate study graphs; SPP1-state later than Monocyte in all five studies.
  - Phase 34 internal myeloid pySCENIC: CEBPD/JUNB state-associated regulon evidence with explicit non-causal boundary.
  - Phase 57–58: internal F1–TYMP direction and GSE207493 stability/weighting audit.
- Required negative/boundary: Phase 21 SPP1-versus-C1QC does not rise in any study; do not claim uniform C1QC replacement or lineage conversion.
- Unresolved gap: no internal-five-study SPP1/C1QC/Other fine-state-specific communication analysis. Phase 35 is broad macrophage; Phase 88 is external. Decide whether to run the internal gap before drafting the final redesigned content.
- Do not resume from Phase 91–100 reports. Use `CLIENT_RESPONSE_BLUEPRINT_v1.md`, `CLIENT_INTERNAL_MAINLINE_EVIDENCE_v1.tsv` and `ANALYSIS_MODULE_LEDGER_v1.tsv`.
- No current Phase 105 background job remains. The only server process matching project26210 is the historical stopped Phase 17 wrapper with an existing completion marker.
"""


def read_text(path: Path) -> tuple[str, bool]:
    raw = path.read_bytes()
    return raw.decode("utf-8-sig"), raw.startswith(b"\xef\xbb\xbf")


def write_text(path: Path, text: str, bom: bool) -> None:
    path.write_text(text, encoding="utf-8-sig" if bom else "utf-8", newline="")


def append_once(path: Path, marker: str, addendum: str) -> None:
    text, bom = read_text(path)
    count = text.count(marker)
    if count > 1:
        raise RuntimeError(f"duplicate marker: {marker}")
    if count == 0:
        newline = "\r\n" if "\r\n" in text else "\n"
        write_text(path, text.rstrip("\r\n") + addendum.replace("\n", newline) + newline, bom)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest().upper()


append_once(STAGE / "PROJECT_STATUS.md", "## Phase 105 full-project analysis audit and client-scope re-lock (2026-08-31)", STATUS)
append_once(STAGE / "DECISIONS.md", "## 2026-08-31 - Phase 105 internal-single-cell redesign decision", DECISIONS)
append_once(STAGE / "CODEX_HANDOFF.md", "## Phase 105 full-project audit and current client task (2026-08-31)", HANDOFF)

manifest_path = STAGE / "RESULTS_MANIFEST.tsv"
text, bom = read_text(manifest_path)
reader = csv.DictReader(text.splitlines(), delimiter="\t")
header = reader.fieldnames
rows = list(reader)
if header is None:
    raise RuntimeError("manifest header missing")
ids = {row["result_id"] for row in rows}
timestamp = datetime.now().astimezone().isoformat(timespec="seconds")


def artifact_type(path: Path) -> str:
    if path.suffix == ".py":
        return "formal_script"
    if path.suffix == ".json":
        return "audit_status_json"
    if path.suffix == ".tsv":
        return "audit_table_tsv"
    if path.suffix == ".txt":
        return "process_snapshot_txt"
    return "audit_document_md"


for path in sorted(item for item in AUDIT.iterdir() if item.is_file()):
    stem = re.sub(r"[^a-z0-9]+", "_", path.stem.lower()).strip("_")
    result_id = f"phase105_full_audit_{stem}_v1"
    if result_id in ids:
        continue
    row = {
        "result_id": result_id,
        "stage": "phase105_full_analysis_audit",
        "artifact_type": artifact_type(path),
        "absolute_path": str(path),
        "generated_at": timestamp,
        "validation_status": "PASS",
        "validation_evidence": f"SHA256 {sha256(path)}; Phase105 audit validation 56/56 PASS",
        "notes": "Full-project scan, scope reconciliation and internal-single-cell client redesign evidence",
    }
    values = [row[column] for column in header]
    if any("\t" in value or "\r" in value or "\n" in value for value in values):
        raise RuntimeError(f"invalid TSV value: {result_id}")
    text = text.rstrip("\r\n") + "\r\n" + "\t".join(values) + "\r\n"
    ids.add(result_id)
write_text(manifest_path, text, bom)

print("PHASE105_STATE_STAGING_UPDATE_PASS")
