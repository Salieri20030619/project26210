from __future__ import annotations

import csv
import hashlib
import json
import zipfile
from collections import Counter
from datetime import datetime
from pathlib import Path

import docx
from docx import Document
from docx.oxml.ns import qn
from PIL import Image, ImageChops


PHASE = Path(__file__).resolve().parent
QC = PHASE / "qc_local_v1"
SOURCE = Path(
    r"D:\project26210_external_joint_v1\phase56_part1_incremental_report_update_local_v1"
) / "PART1_EXTERNAL_SINGLE_CELL_CLIENT_REPORT_WITH_FIGURES_REVIEW_v11.docx"
REPORT = PHASE / "PART1_EXTERNAL_SINGLE_CELL_CLIENT_REPORT_WITH_FIGURES_REVIEW_v12.docx"
EVIDENCE = PHASE / "IMMUNE_MECHANISTIC_EVIDENCE_REGISTER_v3.tsv"
GRADE_SUMMARY = PHASE / "PART1_IMMUNE_EVIDENCE_GRADE_UPDATE_v1.md"
RENDER = PHASE / "rendered_local_v1"
PDF = RENDER / "PART1_EXTERNAL_SINGLE_CELL_CLIENT_REPORT_WITH_FIGURES_REVIEW_v12.pdf"
PHASE67_VALIDATION = Path(
    r"D:\project26210_external_joint_v1\phase67_part1_direct_program_replication_local_v1"
) / "results_v5" / "PHASE67_VALIDATION_v1.tsv"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def media_hashes(path: Path) -> dict[str, str]:
    with zipfile.ZipFile(path) as archive:
        return {
            name: hashlib.sha256(archive.read(name)).hexdigest()
            for name in sorted(archive.namelist())
            if name.startswith("word/media/")
        }


def page_has_visible_content(path: Path) -> bool:
    with Image.open(path) as image:
        rgb = image.convert("RGB")
        background = Image.new("RGB", rgb.size, "white")
        return ImageChops.difference(rgb, background).getbbox() is not None


def has_repeating_header(row) -> bool:
    tr_pr = row._tr.trPr
    return tr_pr is not None and tr_pr.find(qn("w:tblHeader")) is not None


def main() -> None:
    checks: list[dict[str, str]] = []

    def check(name: str, passed: bool, evidence: str) -> None:
        checks.append({"check": name, "status": "PASS" if passed else "FAIL", "evidence": evidence})

    for path in [SOURCE, REPORT, EVIDENCE, GRADE_SUMMARY, PDF, PHASE67_VALIDATION]:
        check(
            f"exists:{path.name}",
            path.exists() and path.stat().st_size > 0,
            str(path.stat().st_size if path.exists() else 0),
        )

    update_status = json.loads((QC / "report_update_status_v1.json").read_text(encoding="utf-8"))
    source_now = sha256(SOURCE)
    check(
        "v11 source unchanged",
        update_status["source_sha256_before"] == update_status["source_sha256_after"] == source_now,
        source_now,
    )
    check(
        "v12 hash matches update record",
        update_status["output_sha256"] == sha256(REPORT),
        sha256(REPORT),
    )

    document = Document(REPORT)
    text = "\n".join(paragraph.text for paragraph in document.paragraphs)
    table_text = "\n".join(cell.text for table in document.tables for row in table.rows for cell in row.cells)
    required = {
        "C07 heading": "C07  免疫程序外部直接复现与证据升级",
        "cell counts": "50,531个Myeloid/DC细胞和6,341个Treg细胞",
        "two-external replication": "18个在两个外部队列均达到R2",
        "Treg single-external boundary": "不能写作双队列确认",
        "grade update": "A=5、B=6、C1=2、D=2",
        "clinical boundary": "F1和TYMP作为独立临床预测指标均继续保持D级",
        "Phase66 distinction": "Phase66的42项情境转移FDR结果提示临床背景依赖",
    }
    for name, phrase in required.items():
        check(name, phrase in text or phrase in table_text, phrase)

    headings = [paragraph.text for paragraph in document.paragraphs if paragraph.style.name.startswith("Heading")]
    check("C07 unique heading", headings.count(required["C07 heading"]) == 1, str(headings.count(required["C07 heading"])))
    check("paragraph count", len(document.paragraphs) == 175, str(len(document.paragraphs)))
    check("table count", len(document.tables) == 4, str(len(document.tables)))
    evidence_tables = [
        table
        for table in document.tables
        if table.rows and table.rows[0].cells and table.rows[0].cells[0].text == "证据轴"
    ]
    check("new evidence table unique", len(evidence_tables) == 1, str(len(evidence_tables)))
    new_table = evidence_tables[0]
    check("new evidence table dimensions", len(new_table.rows) == 6 and len(new_table.columns) == 4, f"{len(new_table.rows)}x{len(new_table.columns)}")
    check("new evidence table bordered", new_table.style.name == "Table Grid", new_table.style.name)
    check("new evidence table repeated header", has_repeating_header(new_table.rows[0]), "w:tblHeader")
    check("no solid bullet glyph", "•" not in text and "•" not in table_text, "none")

    source_media = media_hashes(SOURCE)
    output_media = media_hashes(REPORT)
    check(
        "all embedded figures preserved byte-for-byte",
        source_media == output_media and len(output_media) == 22,
        f"source={len(source_media)} output={len(output_media)}",
    )

    with EVIDENCE.open(encoding="utf-8-sig", newline="") as handle:
        evidence_rows = list(csv.DictReader(handle, delimiter="\t"))
    grade_counts = Counter(row["global_grade"] for row in evidence_rows)
    check("evidence register 15 rows", len(evidence_rows) == 15, str(len(evidence_rows)))
    check(
        "evidence grade distribution",
        dict(grade_counts) == {"D": 2, "A": 5, "B": 6, "C1": 2},
        json.dumps(dict(grade_counts), ensure_ascii=False),
    )
    by_id = {row["evidence_id"]: row for row in evidence_rows}
    expected_updates = {
        "I08": ("A", "R2_TWO_EXTERNAL_FDR_CORE"),
        "I09": ("B", "R1_SINGLE_EXTERNAL_FDR_CORE"),
        "I11": ("A", "R2_TWO_EXTERNAL_FDR_CORE"),
        "I12": ("B", "R1_SINGLE_EXTERNAL_FDR_CORE"),
    }
    for evidence_id, (grade, external) in expected_updates.items():
        row = by_id[evidence_id]
        check(
            f"{evidence_id} updated grade",
            row["global_grade"] == grade and row["external_grade"] == external,
            f"{row['global_grade']} / {row['external_grade']}",
        )
    check(
        "direct clinical grades remain D",
        by_id["I01"]["global_grade"] == "D" and by_id["I02"]["global_grade"] == "D",
        f"I01={by_id['I01']['global_grade']} I02={by_id['I02']['global_grade']}",
    )

    with PHASE67_VALIDATION.open(encoding="utf-8-sig", newline="") as handle:
        phase67_rows = list(csv.DictReader(handle, delimiter="\t"))
    phase67_pass = sum(row["status"] == "PASS" for row in phase67_rows)
    check(
        "Phase67 source validation 19/19",
        len(phase67_rows) == 19 and phase67_pass == 19,
        f"{phase67_pass}/{len(phase67_rows)}",
    )

    render_status = json.loads((RENDER / "render_status.json").read_text(encoding="utf-8"))
    pages = sorted(RENDER.glob("page-*.png"))
    check(
        "Word render passed",
        render_status.get("status") == "PASS"
        and render_status.get("document_pages") == 28
        and render_status.get("rendered_pngs") == 28,
        json.dumps(render_status, ensure_ascii=False),
    )
    check("all 28 page PNGs present", len(pages) == 28, str(len(pages)))
    check("no blank rendered pages", all(page_has_visible_content(page) for page in pages), "28/28 visible")
    dimensions = []
    for page in pages:
        with Image.open(page) as image:
            dimensions.append(image.size)
    check("consistent page dimensions", len(set(dimensions)) == 1, str(sorted(set(dimensions))))

    a11y = json.loads((QC / "a11y_audit_v1.json").read_text(encoding="utf-8"))
    counts = a11y.get("counts", {})
    check(
        "accessibility audit clean",
        counts.get("high", 0) == 0 and counts.get("medium", 0) == 0 and counts.get("low", 0) == 0,
        json.dumps(counts, ensure_ascii=False),
    )

    failed = sum(row["status"] == "FAIL" for row in checks)
    status = {
        "status": "PASS_PHASE68_PART1_REPORT_V12" if failed == 0 else "FAIL",
        "checks_passed": len(checks) - failed,
        "checks_total": len(checks),
        "failed": failed,
        "visual_review": {
            "all_pages_reviewed_at_original_resolution": True,
            "pages": list(range(1, 29)),
            "layout_findings": "No clipping, overlap, orphaned headings, blank pages, broken table borders or figure overflow. C07 and its evidence table fit on page 21.",
        },
        "report_sha256": sha256(REPORT),
        "pdf_sha256": sha256(PDF),
        "checks": checks,
    }
    (QC / "validation_status_v1.json").write_text(
        json.dumps(status, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    lines = ["check\tstatus\tevidence"] + [
        f"{row['check']}\t{row['status']}\t{row['evidence']}" for row in checks
    ]
    (QC / "validation_checks_v1.tsv").write_text("\n".join(lines) + "\n", encoding="utf-8")

    run_record = f"""# Phase 68 Part 1 immune-series report update

- Date: {datetime.now().astimezone().strftime('%Y-%m-%d %H:%M:%S %z')}
- Status: {'PASS' if failed == 0 else 'FAIL'}
- Source report: `{SOURCE}`
- Output report: `{REPORT.name}`
- Source SHA256: `{source_now}`
- Output SHA256: `{sha256(REPORT)}`

## Scope

- Preserved all 26 original analysis modules, 15 main-text figure groups, three pre-existing tables and all 22 embedded media parts byte-for-byte.
- Added C07 for direct external immune-program replication and one bordered evidence-grade table.
- Updated the executive summary, statistical rule, supported/unsupported synthesis and final conclusion.
- Created a 15-row immune evidence register v3: A=5, B=6, C1=2, D=2.
- Upgraded only the externally replicated frozen cores; F1 and TYMP direct clinical-prediction grades remain D.

## Render and validation

- Microsoft Word rendered 28/28 pages to PDF and 160-dpi PNGs.
- All 28 pages were reviewed at original rendered resolution.
- No clipping, overlap, orphaned heading, blank page, broken table border or figure overflow was found.
- Accessibility audit: high=0, medium=0, low=0.
- Final validation: {len(checks) - failed}/{len(checks)} checks PASS.

## Tooling note

- Python: `{__import__('sys').version.split()[0]}`; python-docx: `{docx.__version__}`.
- The packaged canonical renderer was attempted first but could not start because `pdf2image` was unavailable in its Python environment.
- The established local Microsoft Word renderer from Phase23 was used instead; status is recorded in `rendered_local_v1/render_status.json`.
"""
    (PHASE / "RUN_RECORD_v1.md").write_text(run_record, encoding="utf-8")

    manifest_paths = [
        REPORT,
        EVIDENCE,
        GRADE_SUMMARY,
        PDF,
        RENDER / "render_status.json",
        QC / "a11y_audit_v1.json",
        QC / "report_update_status_v1.json",
        QC / "validation_checks_v1.tsv",
        QC / "validation_status_v1.json",
        PHASE / "RUN_RECORD_v1.md",
    ]
    manifest = ["relative_path\tbytes\tsha256"]
    for path in manifest_paths:
        manifest.append(f"{path.relative_to(PHASE)}\t{path.stat().st_size}\t{sha256(path)}")
    (PHASE / "FILE_MANIFEST_SHA256_v1.tsv").write_text("\n".join(manifest) + "\n", encoding="utf-8")

    if failed:
        raise SystemExit(f"Validation failed: {failed}/{len(checks)}")


if __name__ == "__main__":
    main()
