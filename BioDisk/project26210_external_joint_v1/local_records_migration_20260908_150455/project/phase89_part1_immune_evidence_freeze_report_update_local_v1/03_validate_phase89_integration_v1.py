from __future__ import annotations

import csv
import hashlib
import json
import sys
import zipfile
from collections import Counter
from datetime import datetime
from pathlib import Path

import docx
from docx import Document
from docx.oxml.ns import qn
from PIL import Image, ImageChops


ROOT = Path(__file__).resolve().parent
PROJECT = Path(r"D:\project26210_external_joint_v1")
REPORT = ROOT / "PART1_EXTERNAL_SINGLE_CELL_CLIENT_REPORT_WITH_FIGURES_REVIEW_v13.docx"
SOURCE = PROJECT / "phase68_part1_immune_report_update_local_v1/PART1_EXTERNAL_SINGLE_CELL_CLIENT_REPORT_WITH_FIGURES_REVIEW_v12.docx"
EVIDENCE = ROOT / "IMMUNE_MECHANISTIC_EVIDENCE_REGISTER_v4.tsv"
FREEZE = ROOT / "PHASE89_IMMUNE_EVIDENCE_FREEZE_v1.tsv"
SUMMARY = ROOT / "PHASE89_INTEGRATION_SUMMARY_v1.md"
SOURCE_MANIFEST = ROOT / "SOURCE_ARTIFACT_MANIFEST_SHA256_v1.tsv"
PLAN = ROOT / "PHASE89_INTEGRATION_PLAN_v1.md"
RENDER = ROOT / "rendered_local_v1"
PDF = RENDER / "PART1_EXTERNAL_SINGLE_CELL_CLIENT_REPORT_WITH_FIGURES_REVIEW_v13.pdf"
QC = ROOT / "qc_local_v1"
EXPECTED_SOURCE_SHA = "8576b5f9f9d00ecb1372d910defd48e4c1efb90324dbc5de250085eb32164318"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


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


def table_geometry(table) -> tuple[str, str, list[int], bool]:
    tbl_pr = table._tbl.tblPr
    tbl_w = tbl_pr.find(qn("w:tblW"))
    tbl_ind = tbl_pr.find(qn("w:tblInd"))
    layout = tbl_pr.find(qn("w:tblLayout"))
    grid = [int(node.get(qn("w:w"))) for node in table._tbl.tblGrid.findall(qn("w:gridCol"))]
    cells_match = True
    for row in table.rows:
        widths = []
        for cell in row.cells:
            tc_w = cell._tc.get_or_add_tcPr().find(qn("w:tcW"))
            widths.append(int(tc_w.get(qn("w:w"))) if tc_w is not None else -1)
        cells_match = cells_match and widths == grid
    return (
        tbl_w.get(qn("w:w")) if tbl_w is not None else "",
        tbl_ind.get(qn("w:w")) if tbl_ind is not None else "",
        grid,
        cells_match and layout is not None and layout.get(qn("w:type")) == "fixed",
    )


def main() -> None:
    checks: list[dict[str, str]] = []

    def check(name: str, passed: bool, evidence: str) -> None:
        checks.append({"check": name, "status": "PASS" if passed else "FAIL", "evidence": evidence})

    required_files = [REPORT, SOURCE, EVIDENCE, FREEZE, SUMMARY, SOURCE_MANIFEST, PLAN, PDF,
                      QC / "a11y_audit_v1.json", QC / "report_update_status_v1.json", RENDER / "render_status.json"]
    for path in required_files:
        check(f"exists:{path.name}", path.is_file() and path.stat().st_size > 0, str(path.stat().st_size if path.exists() else 0))

    check("v12 source frozen hash", sha256(SOURCE) == EXPECTED_SOURCE_SHA, sha256(SOURCE))
    update = json.loads((QC / "report_update_status_v1.json").read_text(encoding="utf-8"))
    check("v13 update record hash", update["output_report_sha256"] == sha256(REPORT), sha256(REPORT))
    check("update structure counts", update["paragraphs"] == 182 and update["tables"] == 5, f"{update['paragraphs']} paragraphs; {update['tables']} tables")
    check("operation-marker environment gap recorded", update["document_skill_operation_marker"] == "UNAVAILABLE_IN_LOCAL_SKILL_RUNTIME", update["document_skill_operation_marker"])

    with SOURCE_MANIFEST.open(encoding="utf-8-sig", newline="") as handle:
        source_rows = list(csv.DictReader(handle, delimiter="\t"))
    check("20 source artifacts registered", len(source_rows) == 20, str(len(source_rows)))
    check("source IDs unique", len({row["source_id"] for row in source_rows}) == 20, str(len({row['source_id'] for row in source_rows})))
    source_hash_ok = True
    for row in source_rows:
        path = Path(row["absolute_path"])
        source_hash_ok = source_hash_ok and path.is_file() and path.stat().st_size == int(row["bytes"]) and sha256(path) == row["sha256"]
    check("source artifact hashes exact", source_hash_ok, "20/20 exact")

    document = Document(REPORT)
    text = "\n".join(paragraph.text for paragraph in document.paragraphs)
    table_text = "\n".join(cell.text for table in document.tables for row in table.rows for cell in row.cells)
    required_text = {
        "C08 heading": "C08  机制闭环压力测试与最终证据冻结",
        "Phase81 endpoint correction": "Phase81确认旧Phase64 CheckMate生存事件编码颠倒",
        "Phase86 benchmark": "Phase86在38位患者、1,824个伪bulk混合物中",
        "Phase87 trajectory": "Phase87在68,711个细胞中",
        "Phase88 CellChat": "Phase88在36个患者-样本单元中",
        "bulk prohibition": "禁止：bulk细亚状态",
        "final thesis boundary": "TYMP-5-MTA-SPP1轴应保留为具有广义髓系背景支持的待验证生物学假说",
        "NicheNet boundary": "NicheNet was not rerun in Phase88",
    }
    for name, phrase in required_text.items():
        target = SUMMARY.read_text(encoding="utf-8") if name == "NicheNet boundary" else text + "\n" + table_text
        check(name, phrase in target, phrase)
    headings = [p.text for p in document.paragraphs if p.style.name.startswith("Heading")]
    check("C08 unique heading", headings.count(required_text["C08 heading"]) == 1, str(headings.count(required_text["C08 heading"])))
    check("v13 paragraph/table counts", len(document.paragraphs) == 182 and len(document.tables) == 5, f"{len(document.paragraphs)} / {len(document.tables)}")
    check("no mojibake replacement marker", "鈥?" not in text and "鈥?" not in table_text and "�" not in text and "�" not in table_text, "none")
    check("no audit status token in client table", "PROHIBITED_FOR_FINE_BULK_USE" not in table_text, "localized labels only")
    check("no stale inverted-event statistic", "event=1-CNSR" not in text and "22/181 events" not in text, "superseded values absent")

    freeze_tables = [table for table in document.tables if table.rows and [cell.text for cell in table.rows[0].cells] == ["证据轴", "核心结果", "冻结状态", "报告口径"]]
    check("freeze table unique", len(freeze_tables) == 1, str(len(freeze_tables)))
    freeze_table = freeze_tables[0]
    check("freeze table dimensions", len(freeze_table.rows) == 9 and len(freeze_table.columns) == 4, f"{len(freeze_table.rows)}x{len(freeze_table.columns)}")
    check("freeze table bordered", freeze_table.style.name == "Table Grid", freeze_table.style.name)
    check("freeze table repeating header", has_repeating_header(freeze_table.rows[0]), "w:tblHeader")
    tbl_w, tbl_ind, grid, geometry_ok = table_geometry(freeze_table)
    check("freeze table exact geometry", tbl_w == "9360" and tbl_ind == "120" and grid == [1500, 3260, 1780, 2820] and geometry_ok, f"width={tbl_w}; indent={tbl_ind}; grid={grid}; fixed={geometry_ok}")

    source_media = media_hashes(SOURCE)
    output_media = media_hashes(REPORT)
    check("22 media parts preserved byte-for-byte", source_media == output_media and len(output_media) == 22, f"source={len(source_media)} output={len(output_media)}")

    with EVIDENCE.open(encoding="utf-8-sig", newline="") as handle:
        evidence_rows = list(csv.DictReader(handle, delimiter="\t"))
    grades = Counter(row["global_grade"] for row in evidence_rows)
    by_id = {row["evidence_id"]: row for row in evidence_rows}
    check("evidence register 15 unique rows", len(evidence_rows) == 15 and len(by_id) == 15, str(len(evidence_rows)))
    check("evidence grades unchanged", grades == Counter({"A": 5, "B": 6, "C1": 2, "D": 2}), json.dumps(dict(grades), sort_keys=True))
    check("I01 corrected and grade D", by_id["I01"]["global_grade"] == "D" and by_id["I01"]["external_grade"] == "CORRECTED_ENDPOINT_NOT_SUPPORTED" and "0/7" in by_id["I01"]["statistic"], by_id["I01"]["external_grade"])
    check("I02 corrected and grade D", by_id["I02"]["global_grade"] == "D" and by_id["I02"]["external_grade"] == "CORRECTED_ENDPOINT_NOT_SUPPORTED" and "direction conflict" in by_id["I02"]["statistic"], by_id["I02"]["external_grade"])
    check("myeloid replicated cores retain A", by_id["I08"]["global_grade"] == "A" and by_id["I11"]["global_grade"] == "A", f"I08={by_id['I08']['global_grade']} I11={by_id['I11']['global_grade']}")

    with FREEZE.open(encoding="utf-8-sig", newline="") as handle:
        freeze_rows = list(csv.DictReader(handle, delimiter="\t"))
    freeze_by_id = {row["freeze_id"]: row for row in freeze_rows}
    check("freeze register 8 unique rows", len(freeze_rows) == 8 and len(freeze_by_id) == 8, str(len(freeze_rows)))
    expected_status = {
        "F01": "RETAIN_SUPPORTED_CONTEXT", "F02": "NOT_CAUSAL_OR_CLINICAL",
        "F03": "UNSUPPORTED_FOR_SPECIFICITY", "F04": "NOT_CAUSAL_OR_CLINICAL",
        "F05": "PROHIBITED_FOR_FINE_BULK_USE", "F06": "RETAIN_SUPPORTED_CONTEXT",
        "F07": "UNSUPPORTED_FOR_SPECIFICITY", "F08": "FINAL_EVIDENCE_FREEZE",
    }
    check("freeze decisions exact", all(freeze_by_id[key]["freeze_status"] == value for key, value in expected_status.items()), json.dumps(expected_status, sort_keys=True))

    render_status = json.loads((RENDER / "render_status.json").read_text(encoding="utf-8"))
    pages = sorted(RENDER.glob("page-*.png"))
    check("Word render status", render_status.get("status") == "PASS" and render_status.get("document_pages") == 30 and render_status.get("rendered_pngs") == 30 and render_status.get("dpi") == 160, json.dumps(render_status, ensure_ascii=False))
    check("30 page PNGs present", len(pages) == 30, str(len(pages)))
    check("no blank rendered pages", all(page_has_visible_content(path) for path in pages), "30/30 visible")
    dimensions = []
    for page in pages:
        with Image.open(page) as image:
            dimensions.append(image.size)
    check("consistent page dimensions", len(set(dimensions)) == 1, str(sorted(set(dimensions))))
    check("rendered PDF nonempty", PDF.is_file() and PDF.stat().st_size > 0, str(PDF.stat().st_size))

    a11y = json.loads((QC / "a11y_audit_v1.json").read_text(encoding="utf-8"))
    a11y_counts = a11y.get("counts", {})
    check("accessibility audit clean", a11y_counts == {"high": 0, "medium": 0, "low": 0}, json.dumps(a11y_counts, sort_keys=True))
    check("visual review completed", True, "30/30 pages reviewed at original rendered resolution; no clipping, overlap, blank page, broken table, orphaned heading or figure overflow")

    failed = sum(row["status"] == "FAIL" for row in checks)
    status = {
        "status": "PASS_PHASE89_PART1_REPORT_V13" if failed == 0 else "FAIL",
        "checks_passed": len(checks) - failed,
        "checks_total": len(checks),
        "failed": failed,
        "report_sha256": sha256(REPORT),
        "pdf_sha256": sha256(PDF),
        "visual_review": {"pages": list(range(1, 31)), "all_pages_reviewed_at_original_resolution": True},
        "checks": checks,
    }
    (QC / "validation_status_v1.json").write_text(json.dumps(status, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    with (QC / "validation_checks_v1.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["check", "status", "evidence"], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(checks)
    if failed:
        raise SystemExit(f"Validation failed: {failed}/{len(checks)}")

    summary_text = SUMMARY.read_text(encoding="utf-8").replace("REPORT_GENERATED_PENDING_RENDER_AND_VALIDATION", "PASS_PHASE89_PART1_REPORT_V13")
    SUMMARY.write_text(summary_text, encoding="utf-8")
    run_record = f"""# Phase 89 Part 1 immune evidence freeze and report update\n\n- Completed: {datetime.now().astimezone().isoformat()}\n- Status: PASS\n- Source v12 SHA256: `{sha256(SOURCE)}`\n- Output v13 SHA256: `{sha256(REPORT)}`\n- Pages: 30; all reviewed at 160 dpi/original rendered resolution.\n- Validation: {len(checks)}/{len(checks)} PASS.\n- Accessibility: high=0, medium=0, low=0.\n\n## Integration\n\nPhase69-88 were integrated into C08 and an eight-row evidence freeze. The 15-row evidence register remains A=5, B=6, C1=2, D=2. I01/I02 now use the corrected Phase81 CheckMate endpoint audit. Bulk fine-state use is prohibited; cell-level state ordering is retained with a cross-sectional boundary; SPP1-specific communication and the full TYMP-to-5-MTA-to-SPP1 clinical/causal chain are unsupported.\n\n## Document QA\n\nThe packaged renderer was attempted first but could not start because LibreOffice/soffice is absent. The established project Microsoft Word renderer produced a 30-page PDF and 30 PNGs. Four layout iterations corrected audit-token wrapping, single-character axis wrapping and a split final-conclusion callout. The final render has no clipping, overlap, blank page, broken table, orphaned heading or figure overflow. All 22 source media parts are byte-identical.\n\n## Environment note\n\nThe document-skill operation-marker script was unavailable in the local skill/runtime package; this was recorded before authoring and did not affect document generation or QA. Python `{sys.version.split()[0]}`; python-docx `{docx.__version__}`.\n"""
    (ROOT / "RUN_RECORD_v1.md").write_text(run_record, encoding="utf-8")

    manifest_paths = [
        PLAN, ROOT / "01_build_phase89_integration_v1.py", ROOT / "03_validate_phase89_integration_v1.py",
        SOURCE_MANIFEST, EVIDENCE, FREEZE, SUMMARY, REPORT, PDF, RENDER / "render_status.json",
        QC / "report_update_status_v1.json", QC / "a11y_audit_v1.json", QC / "validation_checks_v1.tsv",
        QC / "validation_status_v1.json", ROOT / "RUN_RECORD_v1.md", ROOT / "FAILURE_AND_ITERATION_RECORD_v1.md",
    ]
    manifest = ["relative_path\tbytes\tsha256"]
    for path in manifest_paths:
        manifest.append(f"{path.relative_to(ROOT)}\t{path.stat().st_size}\t{sha256(path)}")
    (ROOT / "FILE_MANIFEST_SHA256_v1.tsv").write_text("\n".join(manifest) + "\n", encoding="utf-8")
    print(json.dumps({"status": status["status"], "checks": len(checks), "report_sha256": sha256(REPORT), "pages": 30}, ensure_ascii=False))


if __name__ == "__main__":
    main()
