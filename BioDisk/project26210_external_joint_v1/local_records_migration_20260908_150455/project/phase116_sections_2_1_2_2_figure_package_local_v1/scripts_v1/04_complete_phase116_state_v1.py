from __future__ import annotations

import csv
import hashlib
import json
import shutil
from datetime import datetime
from pathlib import Path


ROOT = Path(r"D:\project26210_external_joint_v1")
PHASE = ROOT / "phase116_sections_2_1_2_2_figure_package_local_v1"
DELIV = PHASE / "deliverable_v1"
SCRIPTS = PHASE / "scripts_v1"
VALID = PHASE / "validation_v1"
BACKUP = PHASE / "state_backup_before_phase116_v1"
DOWNLOAD = Path(r"C:\Users\31937\Downloads\项目26210_2.1-2.2图件证据包_20260902_v1")
ZIP = Path(r"C:\Users\31937\Downloads\项目26210_2.1-2.2图件证据包_20260902_v1.zip")
STATE = ["PROJECT_STATUS.md", "DECISIONS.md", "RESULTS_MANIFEST.tsv", "CODEX_HANDOFF.md"]
STAMP = "2026-09-02T17:12:00+08:00"


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest().upper()


def append_once(path: Path, marker: str, text: str):
    current = path.read_text(encoding="utf-8-sig")
    if marker in current:
        raise RuntimeError(f"marker already exists in {path.name}: {marker}")
    path.write_text(current.rstrip() + "\n\n" + text.strip() + "\n", encoding="utf-8")


def manifest_row(result_id, stage, artifact_type, path, validation, note, hash_value=None):
    p = Path(path)
    return {
        "result_id": result_id,
        "stage": stage,
        "artifact_type": artifact_type,
        "absolute_path": str(p),
        "generated_at": STAMP,
        "validation_status": "PASS",
        "validation_evidence": validation,
        "notes": note,
        "_hash": hash_value if hash_value is not None else (sha(p) if p.is_file() else ""),
    }


def main():
    for name in STATE:
        if not (ROOT / name).is_file():
            raise FileNotFoundError(ROOT / name)
    if BACKUP.exists():
        raise FileExistsError(BACKUP)
    BACKUP.mkdir(parents=True)
    for name in STATE:
        shutil.copy2(ROOT / name, BACKUP / name)

    summary = PHASE / "PHASE116_PACKAGE_SUMMARY_v1.md"
    summary.write_text(
        """# Phase 116 — 《分析待做（新）》2.1–2.2图件证据包

## 范围

- 仅覆盖参考文档前半部分“恶性上皮”2.1和2.2；后文髓系章节的重复编号2.1不在本包。
- 打包41个正式图件文件：23个PNG、18个PDF；包含4张主管汇报主图、2.1 CNV/亚克隆支撑图、8张正式轨迹图及2.2 cNMF支撑图。
- 未运行新的生物学分析；复用已验证结果并重新核对来源、统计单位、效应和解释边界。

## 独立统计单位与主要结果

- CNV身份层：62,652个上皮细胞，其中38,670个肿瘤上皮细胞；细胞计数仅作描述。
- cNMF/恶性程序层：11,950个正式恶性细胞；患者层推断最多20名患者，克隆比较14名患者。
- 支持：F1–TYMP在20名患者中17名正向，中位rho=0.115，P=0.0153，18端点BH q=0.0197；GSE207493为10/10正向，中位rho=0.147。
- 边界：GSE207493提供99.66%的F1主导细胞；F1–CNV burden中位rho=-0.033、BH q=0.923；主克隆阈值0.25下F1高克隆TYMP差值BH q=0.0735；五方法轨迹方向不一致；F1高克隆差异基因ORA全部FDR=1。
- 正式结论：现有数据支持TYMP RNA与cNMF rank-7 F1转录程序相关，但当前数据不支持“TYMP促使肿瘤进入更恶性状态”的因果结论。

## 验证与交付

- 图件索引41/41，文件非空且SHA-256一致。
- 说明DOCX经Microsoft Word渲染9/9页；无空白、裁切或跨页断行；结构/内容/哈希验证23/23 PASS；可访问性high=0、medium=0、low=0。
- D盘权威包与Downloads逐文件70/70一致；ZIP含70个条目。
""",
        encoding="utf-8",
    )

    status_section = r"""## Phase 116 sections 2.1–2.2 figure evidence package (2026-09-02)

- Status: **PASS_PACKAGE_COMPLETE; 41 FIGURE FILES; 23/23 VALIDATION PASS; WORD RENDER 9/9; A11Y 0/0/0; DOWNLOAD COPY 70/70 EXACT**. No new biological analysis was executed.
- Scope: only the first malignant-epithelial sections 2.1 and 2.2 in `分析待做（新）.docx` (SHA-256 A32766F2A477D6112B6B57A51A647C70DC74961830F320E698361DBE74DA7753). The later myeloid section numbered 2.1 is excluded.
- Package: 41 validated figure files (23 PNG, 18 PDF), comprising four supervisor-facing composite figures, CNV/malignancy/subclone support, eight accepted malignant-trajectory figures and cNMF DE/pathway support. Superseded Slingshot v1 and clone-tree v2/v3 were not included; clone-tree v4 is authoritative.
- Independent units: CNV identity covers 62,652 epithelial cells/38,670 tumour epithelial cells descriptively; cNMF/clone analyses use 11,950 formal malignant cells, with patient-level inference up to n=20 and clone comparison n=14. Cells are not biological replicates.
- Principal support: F1–TYMP is positive in 17/20 patients (median rho=0.115, P=0.0153, 18-endpoint BH q=0.0197); GSE207493 is 10/10 positive with median rho=0.147 and robust leave-one-patient/balanced-resampling direction.
- Principal negative/contradictory results: F1 does not track inferCNV burden (median rho=-0.033, BH q=0.923); the primary clone TYMP contrast does not pass the seven-factor BH family (q=0.0735); malignant trajectory directions differ by method; F1-high clone DE pathway ORA is FDR-negative. GSE207493 supplies 99.66% of F1-dominant cells.
- Evidence decision: current data support a TYMP-associated cNMF F1 transcriptional state, not the causal claim that TYMP drives tumour cells into a more malignant state.
- Authoritative deliverable: `D:\project26210_external_joint_v1\phase116_sections_2_1_2_2_figure_package_local_v1\deliverable_v1`. Downloads folder and ZIP are exact delivery copies.
- Current breakpoint: discuss this evidence package with the supervisor. Do not start new analysis until the supervisor selects whether to pursue perturbational/longitudinal validation, a revised association-only narrative, or no further computation. Active project26210 jobs: none; the historical stopped Phase17 wrapper remains untouched.
"""
    append_once(ROOT / "PROJECT_STATUS.md", "## Phase 116 sections 2.1–2.2 figure evidence package", status_section)

    decision_section = r"""## 2026-09-02 - Phase 116 2.1/2.2 evidence-package decisions

- Interpret the user's screenshot as referring to the first malignant-epithelial 2.1/2.2 block, not the later duplicated myeloid 2.1 numbering.
- Package all current formal figures relevant to CNV malignancy/subclones, TYMP–clone contrasts, accepted malignant trajectories and cNMF rank/function. Exclude superseded/failed figure versions; use clone-tree v4 and Slingshot v2.
- Accept the wording that TYMP RNA is associated with the cNMF rank-7 F1 transcriptional state, with patient-level directional support and substantial GSE207493 weighting.
- Reject the wording that TYMP has been proved to drive a more malignant state. F1–CNV burden is FDR-negative, clone-level evidence is exploratory, trajectory directions conflict and no perturbational/longitudinal causal test exists.
- Keep top-loading-gene pathway annotation separate from patient-blocked F1-high-clone differential pathway evidence; the latter is FDR-negative.
- Phase116 is packaging/documentation only. No new analysis is authorized by this deliverable; the next route requires supervisor decision.
"""
    append_once(ROOT / "DECISIONS.md", "## 2026-09-02 - Phase 116 2.1/2.2 evidence-package decisions", decision_section)

    handoff_section = r"""## Phase 116 2.1–2.2 figure package current breakpoint (2026-09-02)

- Completed: packaged all current formal figures corresponding to the first malignant-epithelial sections 2.1 and 2.2 of `分析待做（新）.docx` and created a supervisor-facing evidence explanation DOCX/PDF.
- Package contents: 41 figure files (23 PNG, 18 PDF), 19 key statistical/source files, README, per-file figure index with source paths/hashes, package manifest and validation records.
- Evidence decision: TYMP RNA–cNMF F1 association is supported with cohort heterogeneity; TYMP-driven transition to a more malignant state is not established. Preserve the F1–CNV null, clone-threshold sensitivity, cross-method trajectory conflict and FDR-negative differential pathway ORA.
- Validation: Word 9/9 pages, 23/23 package checks, a11y 0/0/0, Downloads 70/70 exact files, ZIP 70 entries.
- Exact resume point: supervisor review. If additional evidence is requested, first specify whether the target is association-only synthesis, perturbational/longitudinal causal validation or a different malignant-state definition. Do not continue same-data algorithm expansion to obtain a preferred direction.
- Active project26210 jobs: none. The historical Phase17 stopped wrapper remains unchanged.
- Key paths: `D:\project26210_external_joint_v1\phase116_sections_2_1_2_2_figure_package_local_v1\deliverable_v1`, `C:\Users\31937\Downloads\项目26210_2.1-2.2图件证据包_20260902_v1`, and the adjacent ZIP.
"""
    append_once(ROOT / "CODEX_HANDOFF.md", "## Phase 116 2.1–2.2 figure package current breakpoint", handoff_section)

    validation_note = "Phase116 package validation 23/23 PASS; Word render 9/9; a11y 0/0/0; Downloads exact 70/70; ZIP 70 entries"
    rows = [
        manifest_row("phase116_package_summary_v1", "phase116_2_1_2_2_figure_package", "phase_summary_md", summary, validation_note, "Scope, sample sizes, support and boundaries"),
        manifest_row("phase116_explanation_docx_v1", "phase116_2_1_2_2_figure_package", "supervisor_evidence_explanation_docx", DELIV / "项目26210_2.1-2.2图件证据说明_v1.docx", validation_note, "Editable explanation with four main figures"),
        manifest_row("phase116_explanation_pdf_v1", "phase116_2_1_2_2_figure_package", "rendered_explanation_pdf", DELIV / "项目26210_2.1-2.2图件证据说明_v1.pdf", validation_note, "Microsoft Word fixed-format render"),
        manifest_row("phase116_figure_index_v1", "phase116_2_1_2_2_figure_package", "figure_index_tsv", DELIV / "图件索引与证据边界_v1.tsv", validation_note, "41 figure files with source, role, boundary and SHA-256"),
        manifest_row("phase116_package_manifest_v1", "phase116_2_1_2_2_figure_package", "file_manifest_tsv", DELIV / "PACKAGE_FILE_MANIFEST_SHA256_v1.tsv", validation_note, "Per-file delivery manifest"),
        manifest_row("phase116_package_validation_v1", "phase116_2_1_2_2_figure_package", "validation_json", VALID / "PACKAGE_VALIDATION_v1.json", validation_note, "23/23 package and render checks"),
        manifest_row("phase116_a11y_validation_v1", "phase116_2_1_2_2_figure_package", "accessibility_json", VALID / "a11y_audit_v2.json", validation_note, "Accessibility high=0 medium=0 low=0"),
        manifest_row("phase116_reference_extract_v1", "phase116_2_1_2_2_figure_package", "reference_structure_json", PHASE / "reference_extract_v1" / "reference_structure_v1.json", validation_note, "Identifies first 2.1/2.2 and duplicated later 2.1"),
        manifest_row("phase116_builder_v1", "phase116_2_1_2_2_figure_package", "formal_script", SCRIPTS / "02_build_figure_package_document_v1.py", validation_note, "Builds index, README and DOCX"),
        manifest_row("phase116_validator_v1", "phase116_2_1_2_2_figure_package", "formal_validation_script", SCRIPTS / "03_validate_figure_package_v1.py", validation_note, "Validates hashes, OOXML, render and pixel boundaries"),
        manifest_row("phase116_download_zip_v1", "phase116_2_1_2_2_figure_package", "delivered_zip", ZIP, validation_note, "70-entry send-ready ZIP; SHA256 77E9A69BA40B5F205802A1157D0BE52F39F7EE13A3E19A37FC08C0E99141F031"),
        manifest_row("phase116_download_docx_v1", "phase116_2_1_2_2_figure_package", "delivered_docx", DOWNLOAD / "项目26210_2.1-2.2图件证据说明_v1.docx", validation_note, "Exact Downloads copy"),
        manifest_row("phase116_download_pdf_v1", "phase116_2_1_2_2_figure_package", "delivered_pdf", DOWNLOAD / "项目26210_2.1-2.2图件证据说明_v1.pdf", validation_note, "Exact Downloads copy"),
    ]

    manifest = ROOT / "RESULTS_MANIFEST.tsv"
    with manifest.open("r", encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        existing = list(reader)
        header = reader.fieldnames
    expected = ["result_id", "stage", "artifact_type", "absolute_path", "generated_at", "validation_status", "validation_evidence", "notes"]
    if header != expected:
        raise RuntimeError(f"manifest header mismatch: {header}")
    existing_ids = {r["result_id"].strip('"') for r in existing}
    for row in rows:
        if row["result_id"] in existing_ids:
            raise RuntimeError(f"duplicate new id: {row['result_id']}")
    with manifest.open("a", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=expected, delimiter="\t", lineterminator="\n")
        for row in rows:
            writer.writerow({k: row[k] for k in expected})

    # Validate postconditions before registering the self-referential gate.
    with manifest.open("r", encoding="utf-8-sig", newline="") as f:
        final_rows = list(csv.DictReader(f, delimiter="\t"))
    ids = [r["result_id"].strip('"') for r in final_rows]
    checks = {
        "manifest_schema": True,
        "manifest_unique_ids": len(ids) == len(set(ids)),
        "new_paths_exist": all(Path(r["absolute_path"]).exists() for r in rows),
        "new_hashes_match": all((not r["_hash"]) or sha(Path(r["absolute_path"])) == r["_hash"] for r in rows),
        "project_status_marker": "## Phase 116 sections 2.1–2.2 figure evidence package" in (ROOT / "PROJECT_STATUS.md").read_text(encoding="utf-8"),
        "decisions_marker": "## 2026-09-02 - Phase 116 2.1/2.2 evidence-package decisions" in (ROOT / "DECISIONS.md").read_text(encoding="utf-8"),
        "handoff_marker": "## Phase 116 2.1–2.2 figure package current breakpoint" in (ROOT / "CODEX_HANDOFF.md").read_text(encoding="utf-8"),
        "backup_complete": all((BACKUP / name).is_file() for name in STATE),
        "package_validation_pass": json.loads((VALID / "PACKAGE_VALIDATION_v1.json").read_text(encoding="utf-8"))["status"] == "PASS",
        "downloads_manifest_exact": sha(DELIV / "PACKAGE_FILE_MANIFEST_SHA256_v1.tsv") == sha(DOWNLOAD / "PACKAGE_FILE_MANIFEST_SHA256_v1.tsv"),
    }
    gate = {
        "status": "PASS" if all(checks.values()) else "FAIL",
        "checks": checks,
        "n_checks": len(checks),
        "n_failed": sum(not x for x in checks.values()),
        "state_hashes_before_self_registration": {name: sha(ROOT / name) for name in STATE},
        "remote_mirror": "/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1",
        "remote_sync_required_before_completion_report": True,
    }
    gate_path = VALID / "PHASE116_STATE_UPDATE_VALIDATION_v1.json"
    gate_path.write_text(json.dumps(gate, ensure_ascii=False, indent=2), encoding="utf-8")
    if gate["status"] != "PASS":
        raise RuntimeError(json.dumps(gate, ensure_ascii=False))

    gate_row = manifest_row(
        "phase116_state_update_validation_v1",
        "phase116_project_state_update",
        "state_update_validation_json",
        gate_path,
        "Phase116 mandatory four-file completion gate 10/10 PASS; remote mirror sync required before reporting completion",
        "State validation generated before self-registration to avoid self-hash recursion",
        hash_value="",
    )
    with manifest.open("a", encoding="utf-8", newline="") as f:
        csv.DictWriter(f, fieldnames=expected, delimiter="\t", lineterminator="\n").writerow({k: gate_row[k] for k in expected})

    print(json.dumps({"status": "PASS", "new_manifest_rows": len(rows) + 1, "gate": str(gate_path), "backup": str(BACKUP)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
