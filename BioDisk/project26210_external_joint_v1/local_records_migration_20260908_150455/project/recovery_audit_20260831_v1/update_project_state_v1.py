from __future__ import annotations

import csv
import hashlib
from datetime import datetime
from pathlib import Path


STAGE = Path(r"C:\Users\31937\project26210_recovery_audit_20260831_v1")
AUTHORITATIVE = Path(r"D:\project26210_external_joint_v1")


STATUS_ADDENDUM = r"""

## Recovery audit and authoritative breakpoint (2026-08-31)

- A cross-workspace recovery audit confirmed that the project has progressed through **Phase 104**, not Phase 99. Phase 100–104 artifacts, status sections, decisions and manifest rows were checked against the authoritative local root `D:\project26210_external_joint_v1`.
- Phase 100 produced the final six-page TYMP–MTA immune public-validation brief v3; its DOCX/PDF hashes and page-level review remain PASS. Phase 101 recovered and identity-audited the public NCI60, CCL180 and CCLE cell-line resources used by the subsequent modules.
- Phase 102 preprocessing is complete with 14/14 internal checks and 15/15 independent checks. Phase 103 analysis is complete with 10/10 internal checks, 20/20 independent checks and final PNG/PDF visual review. Phase 104 is complete with 12/12 internal checks, 18/18 independent checks and final PNG/PDF visual review.
- The latest scientific breakpoint is therefore the Phase 104 decision: direct MTA is not evaluable in CCL180 or CCLE; zero of 12 adjacent-metabolite relations passed both prespecified primary and tissue-adjusted FDR gates. The NCI60 renal TYMP–MTA direction remains an n=8 exploratory candidate and is not a 786-0 perturbation result.
- `CODEX_HANDOFF.md` was newly created at the authoritative root because no current handoff file existed there. It records the Phase 104 breakpoint, validated evidence, operational paths and permitted next actions.
- The four key state files are mirrored to `/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1` on `sgs10-codex` for recovery. The local D-drive versions remain authoritative because manifest artifact paths are local Windows paths; the server copies are synchronized state mirrors, not an independent result installation.
- This recovery module performed no new biological analysis and did not alter any Phase 100–104 source artifact or scientific result.
"""


DECISION_ADDENDUM = r"""

## 2026-08-31 - project recovery and authoritative-state decision

- Treat `D:\project26210_external_joint_v1` as the authoritative current project root and Phase 104 as the latest completed scientific phase.
- Treat `/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1` on `sgs10-codex` as a synchronized recovery mirror for the four key state files. Do not reinterpret Windows manifest paths as server-local analysis paths.
- Do not resume from the older `/mnt/BioDisk/ludingyu/project26210` records or the Phase 98/99 workspace without first reconciling them to the Phase 104 state.
- Preserve Phase 100–104 source artifacts and their completed validation records unchanged. The recovery audit adds project-state documentation only.
- Keep the current 786-0 boundary explicit: the public NCI60 observation is a single baseline cell line with all three MTA repeats at the observed floor, not a TYMP-knockdown result. Direct perturbation evidence with matched controls is still required for a causal claim.
"""


SUMMARY = """# project26210 recovery audit summary

- Audit date: 2026-08-31.
- Authoritative root: `D:\\project26210_external_joint_v1`.
- Latest completed scientific module: Phase 104.
- Phase 102 independent validation: 15/15 PASS.
- Phase 103 independent validation: 20/20 PASS; final figure review PASS.
- Phase 104 independent validation: 18/18 PASS; final figure review PASS.
- Before this recovery update, the installed `PROJECT_STATUS.md`, `DECISIONS.md` and `RESULTS_MANIFEST.tsv` were byte-identical to the Phase 104 state-update candidates.
- The missing `CODEX_HANDOFF.md` was reconstructed from the current Phase 100–104 source summaries, validation records, decisions and manifest.
- No source analysis artifact or scientific result was modified.
"""


def read_text_preserve(path: Path) -> tuple[str, bool]:
    raw = path.read_bytes()
    bom = raw.startswith(b"\xef\xbb\xbf")
    return raw.decode("utf-8-sig"), bom


def write_text_preserve(path: Path, text: str, bom: bool) -> None:
    encoding = "utf-8-sig" if bom else "utf-8"
    path.write_text(text, encoding=encoding, newline="")


def append_once(path: Path, marker: str, addendum: str) -> None:
    text, bom = read_text_preserve(path)
    if text.count(marker) > 1:
        raise RuntimeError(f"duplicate marker before update: {marker}")
    if marker not in text:
        ending = "\r\n" if "\r\n" in text else "\n"
        normalized = addendum.replace("\n", ending)
        write_text_preserve(path, text.rstrip("\r\n") + normalized + ending, bom)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest().upper()


append_once(
    STAGE / "PROJECT_STATUS.md",
    "## Recovery audit and authoritative breakpoint (2026-08-31)",
    STATUS_ADDENDUM,
)
append_once(
    STAGE / "DECISIONS.md",
    "## 2026-08-31 - project recovery and authoritative-state decision",
    DECISION_ADDENDUM,
)
(STAGE / "RECOVERY_AUDIT_SUMMARY_v1.md").write_text(SUMMARY, encoding="utf-8", newline="\n")

manifest_path = STAGE / "RESULTS_MANIFEST.tsv"
manifest_text, manifest_bom = read_text_preserve(manifest_path)
reader = csv.DictReader(manifest_text.splitlines(), delimiter="\t")
rows = list(reader)
ids = {row["result_id"] for row in rows}
timestamp = datetime.now().astimezone().isoformat(timespec="seconds")

new_rows = [
    {
        "result_id": "project_state_handoff_phase104_recovery_v1",
        "stage": "project_recovery_20260831",
        "artifact_type": "codex_handoff_md",
        "absolute_path": str(AUTHORITATIVE / "CODEX_HANDOFF.md"),
        "generated_at": timestamp,
        "validation_status": "PASS",
        "validation_evidence": f"SHA256 {sha256(STAGE / 'CODEX_HANDOFF.md')}; Phase 100-104 recovery handoff",
        "notes": "Current authoritative breakpoint, evidence boundaries and operational recovery instructions",
    },
    {
        "result_id": "project_state_recovery_audit_summary_v1",
        "stage": "project_recovery_20260831",
        "artifact_type": "recovery_audit_summary_md",
        "absolute_path": str(AUTHORITATIVE / "recovery_audit_20260831_v1" / "RECOVERY_AUDIT_SUMMARY_v1.md"),
        "generated_at": timestamp,
        "validation_status": "PASS",
        "validation_evidence": f"SHA256 {sha256(STAGE / 'RECOVERY_AUDIT_SUMMARY_v1.md')}; source validation statuses reconciled",
        "notes": "Read-only progress audit; no biological analysis artifact changed",
    },
]

header = reader.fieldnames
if header is None:
    raise RuntimeError("manifest header missing")
for row in new_rows:
    if row["result_id"] not in ids:
        fields = [row[column] for column in header]
        if any("\t" in value or "\n" in value or "\r" in value for value in fields):
            raise RuntimeError(f"invalid manifest field: {row['result_id']}")
        manifest_text = manifest_text.rstrip("\r\n") + "\r\n" + "\t".join(fields) + "\r\n"
        ids.add(row["result_id"])
write_text_preserve(manifest_path, manifest_text, manifest_bom)

print("PROJECT_STATE_STAGING_UPDATE_PASS")
