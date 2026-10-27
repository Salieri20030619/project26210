from __future__ import annotations

import csv
import hashlib
import json
import shutil
from datetime import datetime, timezone, timedelta
from pathlib import Path


ROOT = Path(r"D:\project26210_external_joint_v1")
PHASE = ROOT / "phase117_23figure_word_package_correction_local_v1"
DELIV = PHASE / "deliverable_v2"
DOWNLOAD = Path(r"C:\Users\31937\Downloads\项目26210_2.1-2.2图件证据包_23图完整版_20260902_v2")
ZIP_PATH = Path(str(DOWNLOAD) + ".zip")
BACKUP = PHASE / "state_backup_before_phase117_v1"
STATE_FILES = ["PROJECT_STATUS.md", "DECISIONS.md", "RESULTS_MANIFEST.tsv", "CODEX_HANDOFF.md"]
NOW = datetime.now(timezone(timedelta(hours=8))).isoformat(timespec="seconds")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest().upper()


def append_text(path: Path, text: str) -> None:
    with path.open("a", encoding="utf-8", newline="") as fh:
        fh.write("\n" + text.strip() + "\n")


if BACKUP.exists():
    raise FileExistsError(f"Refusing to overwrite backup: {BACKUP}")
for name in STATE_FILES:
    if not (ROOT / name).exists():
        raise FileNotFoundError(ROOT / name)
for p in [DELIV, DOWNLOAD, ZIP_PATH]:
    if not p.exists():
        raise FileNotFoundError(p)

BACKUP.mkdir(parents=True)
for name in STATE_FILES:
    shutil.copy2(ROOT / name, BACKUP / name)

summary = PHASE / "PHASE117_PACKAGE_SUMMARY_v2.md"
summary.write_text(
    """# Phase 117 2.1–2.2图件包修正版摘要\n\n"
    "- 本阶段纠正Phase116交付的可读性问题，不改变任何生物学分析结果。\n"
    "- Word正文实际嵌入23张PNG图，共28页；不是仅列出图件文件。\n"
    "- 18张PDF仅作为独立矢量备份，置于单独目录。\n"
    "- 面向主管的交付目录不含TSV或JSON；关键统计量和证据边界直接写入Word。\n"
    "- 正式校验：包结构18/18 PASS，Word 23张内嵌图，28/28页渲染，a11y high/medium/low=0/0/0，Downloads 44/44哈希一致，ZIP 44项。\n"
    "- 结论边界保持不变：当前数据支持TYMP RNA与cNMF rank-7 F1转录状态相关，不支持宣称TYMP促进肿瘤进入更恶性状态。\n"
    """,
    encoding="utf-8",
)

append_text(
    ROOT / "PROJECT_STATUS.md",
    f"""## Phase 117 2.1–2.2 supervisor package correction (2026-09-02)\n\n"
    "- Status: **PASS_CORRECTED_PACKAGE; WORD 23/23 EMBEDDED FIGURES; WORD RENDER 28/28 PAGES; FORMAL VALIDATION 18/18 PASS; A11Y 0/0/0; DOWNLOADS 44/44 HASH MATCH; ZIP 44 ENTRIES**.\n"
    "- Scope: corrected packaging and readability only; no biological analysis was rerun and no result changed. Phase116 remains an audit record but its user-facing package is superseded by this Phase117 v2 delivery.\n"
    "- Inputs and biological sample sizes are inherited: 11,950 malignant epithelial cells; patient-level F1/TYMP inference up to n=20; clone analysis n=14. Cells are not treated as independent biological replicates.\n"
    "- Deliverable design: the Word document now actually embeds all 23 PNG figures across 28 pages. The 18 PDFs are separated as vector backups. No TSV or JSON is included in the supervisor-facing package; all key statistical explanations are readable directly in Word.\n"
    "- Principal evidence retained: F1–TYMP is positive in 17/20 patients (median Spearman rho=0.115, P=0.0153, 18-endpoint BH q=0.0197); GSE207493 is 10/10 positive (median rho=0.147) but contributes 99.66% of F1-dominant cells. F1 versus inferCNV burden is null (median rho=-0.033, q=0.923). Primary clone contrast at CNV threshold 0.25 has median delta=0.200, 9/14 positive, P=0.021, q=0.0735. Trajectory methods conflict and F1-high clone pathway ORA is FDR-negative.\n"
    "- Evidence boundary: current data support an association between TYMP RNA and the cNMF rank-7 F1 transcriptional state; they do not establish that TYMP drives transition to a more malignant state.\n"
    f"- Authoritative package: `{PHASE}`. Downloads delivery: `{DOWNLOAD}` and adjacent full-name ZIP.\n"
    "- Current breakpoint: supervisor discussion and evidence collection. Active project26210 background jobs: none.""",
)

append_text(
    ROOT / "DECISIONS.md",
    """## Phase 117 corrected 23-figure package decisions (2026-09-02)\n\n"
    "- Accepted: a supervisor evidence package must be readable without opening TSV/JSON files; the Word document itself is the primary explanation artifact.\n"
    "- Accepted: all 23 PNG figures corresponding to the first malignant-epithelial sections 2.1 and 2.2 are embedded in the Word document; 18 PDFs are retained only as separately labelled vector backups.\n"
    "- Superseded for user delivery: the Phase116 v1 package, because its Word document embedded only four figures and exposed raw TSV/JSON files that were not directly usable by the user. Phase116 remains preserved for audit and provenance.\n"
    "- No biological interpretation changed. Retain the positive F1–TYMP association with cohort imbalance/heterogeneity and retain the F1–CNV null, clone-threshold sensitivity, cross-method trajectory conflict and FDR-negative pathway ORA.\n"
    "- Rejected claim remains: TYMP has been proven to promote entry into a more malignant tumour state. The allowed claim remains TYMP RNA association with the cNMF rank-7 F1 transcriptional state.\n"
    "- Current action remains supervisor review; no new computation is authorized by this packaging correction.""",
)

append_text(
    ROOT / "CODEX_HANDOFF.md",
    f"""## Phase 117 corrected 23-figure package current breakpoint (2026-09-02)\n\n"
    "- Completed: non-overwriting v2 supervisor package correcting Phase116 readability. Word contains exactly 23 embedded PNG figures and renders to 28 pages; 18 PDFs are separated as vector backups; the user package contains no TSV or JSON.\n"
    "- Validation: package 18/18 PASS; Word inline images 23; 28/28 rendered pages with no blank/clipped page; accessibility high=0, medium=0, low=0; authoritative-to-Downloads 44/44 exact hash matches; ZIP 44 entries.\n"
    f"- Authoritative package: `{PHASE}`. Delivery folder: `{DOWNLOAD}`. Full-name ZIP: `{ZIP_PATH}`.\n"
    "- Exact resume point: discuss the evidence with the supervisor. If more evidence is requested, first clarify whether the target is association synthesis, perturbational/longitudinal causal evidence or a different malignant-state definition.\n"
    "- Active jobs: none. Do not restart downloads. Claims that must not be made: TYMP-driven malignant-state transition, a consistent trajectory direction across methods, CNV validation of F1, or FDR-significant F1-high clone pathway enrichment.""",
)

validation_path = PHASE / "validation_v2" / "PHASE117_FINAL_STATE_VALIDATION_v2.json"
validation_path.parent.mkdir(exist_ok=True)
validation_path.write_text("{}\n", encoding="utf-8")

artifacts = [
    ("phase117_summary_v2", "package_summary", summary, "Phase117 correction summary"),
    ("phase117_deliverable_dir_v2", "delivered_package", DELIV, "23 embedded PNG figures; 18 separate PDF backups; no TSV/JSON"),
    ("phase117_docx_v2", "delivered_docx", DELIV / "项目26210_2.1-2.2图件与统计说明_23图完整版_v2.docx", "23 embedded images; 28 pages"),
    ("phase117_pdf_v2", "delivered_pdf", DELIV / "项目26210_2.1-2.2图件与统计说明_23图完整版_v2.pdf", "Word fixed-format render; 28 pages"),
    ("phase117_file_manifest_v2", "file_manifest", PHASE / "FILE_MANIFEST_SHA256_v2.tsv", "All 44 user-package files with SHA256"),
    ("phase117_package_validation_v2", "validation_json", PHASE / "validation_v2" / "PACKAGE_VALIDATION_v2.json", "18/18 PASS; 23 inline images; 28 rendered pages"),
    ("phase117_a11y_v2", "accessibility_validation_json", PHASE / "validation_v2" / "a11y_audit_v2.json", "High=0, medium=0, low=0"),
    ("phase117_render_status_v2", "render_validation_json", PHASE / "word_render_v2" / "render_status.json", "Microsoft Word render PASS; 28 pages"),
    ("phase117_install_validation_v2", "delivery_validation_json", PHASE / "INSTALL_AND_DELIVERY_VALIDATION_v2.json", "Downloads exact 44/44; ZIP 44 entries"),
    ("phase117_build_script_v2", "formal_document_script", PHASE / "scripts_v1" / "01_build_23figure_word_v2.py", "Builds Word with all 23 embedded PNG figures"),
    ("phase117_validate_script_v2", "formal_validation_script", PHASE / "scripts_v1" / "02_validate_23figure_package_v2.py", "Validates package, OOXML and render"),
    ("phase117_install_script_v2", "formal_delivery_script", PHASE / "scripts_v1" / "03_install_phase117_package_v2.py", "Non-overwriting package installation and delivery"),
    ("phase117_zip_fix_script_v2", "recovery_script", PHASE / "scripts_v1" / "03b_fix_phase117_zip_name_v2.py", "Corrects dotted-directory ZIP naming without changing contents"),
    ("phase117_state_script_v2", "formal_state_script", PHASE / "scripts_v1" / "04_complete_phase117_state_v2.py", "Updates mandatory four state files and validates gate"),
    ("phase117_download_dir_v2", "delivered_package", DOWNLOAD, "Exact 44-file Downloads copy"),
    ("phase117_download_zip_v2", "delivered_zip", ZIP_PATH, "44-entry send-ready ZIP"),
    ("phase117_final_state_validation_v2", "state_update_validation_json", validation_path, "Mandatory four-file completion gate"),
]

manifest = ROOT / "RESULTS_MANIFEST.tsv"
with manifest.open("a", encoding="utf-8", newline="") as fh:
    writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
    for rid, atype, path, notes in artifacts:
        if not path.exists():
            raise FileNotFoundError(path)
        evidence = "Phase117 formal correction validation"
        if path.is_file() and path != validation_path:
            evidence += f"; SHA256 {sha256(path)}"
        writer.writerow([rid, "phase117_23figure_word_package_correction", atype, str(path), NOW, "PASS", evidence, notes])

rows = list(csv.DictReader(manifest.open("r", encoding="utf-8-sig", newline=""), delimiter="\t"))
ids = [r["result_id"] for r in rows]
phase_rows = [r for r in rows if r["stage"] == "phase117_23figure_word_package_correction"]
checks = {
    "manifest_schema": list(rows[0].keys()) == ["result_id", "stage", "artifact_type", "absolute_path", "generated_at", "validation_status", "validation_evidence", "notes"],
    "manifest_unique_ids": len(ids) == len(set(ids)),
    "phase117_row_count_17": len(phase_rows) == 17,
    "phase117_paths_exist": all(Path(r["absolute_path"]).exists() for r in phase_rows),
    "phase117_status_pass": all(r["validation_status"] == "PASS" for r in phase_rows),
    "state_backup_four_files": all((BACKUP / n).exists() for n in STATE_FILES),
    "project_status_marker": "Phase 117" in (ROOT / "PROJECT_STATUS.md").read_text(encoding="utf-8"),
    "decisions_marker": "Phase 117" in (ROOT / "DECISIONS.md").read_text(encoding="utf-8"),
    "handoff_marker": "Phase 117" in (ROOT / "CODEX_HANDOFF.md").read_text(encoding="utf-8"),
    "word_has_23_inline_images": json.loads((PHASE / "validation_v2" / "PACKAGE_VALIDATION_v2.json").read_text(encoding="utf-8"))["details"]["word_inline_images"] == 23,
    "no_tsv_json_in_user_package": not any(p.suffix.lower() in {".tsv", ".json"} for p in DELIV.rglob("*") if p.is_file()),
    "downloads_hash_match_44": json.loads((PHASE / "INSTALL_AND_DELIVERY_VALIDATION_v2.json").read_text(encoding="utf-8"))["downloads_exact_hash_matches"] == 44,
}
result = {"status": "PASS" if all(checks.values()) else "FAIL", "generated_at": NOW, "n_checks": len(checks), "n_passed": sum(checks.values()), "checks": checks}
validation_path.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
if result["status"] != "PASS":
    raise RuntimeError(json.dumps(result, ensure_ascii=False, indent=2))

shutil.copy2(Path(__file__), PHASE / "scripts_v1" / Path(__file__).name)
print(json.dumps(result, ensure_ascii=False, indent=2))
