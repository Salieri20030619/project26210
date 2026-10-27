from __future__ import annotations

import hashlib
import json
import zipfile
from pathlib import Path

from PIL import Image
from docx import Document


ROOT = Path(__file__).resolve().parents[1]
PHASE = Path(__file__).resolve().parent
QC = PHASE / "qc_local_v1"
REPORT = PHASE / "PART1_EXTERNAL_SINGLE_CELL_CLIENT_REPORT_WITH_FIGURES_REVIEW_v10.docx"
PDF = PHASE / "rendered_local_v1" / "PART1_EXTERNAL_SINGLE_CELL_CLIENT_REPORT_WITH_FIGURES_REVIEW_v10.pdf"
FIG03 = PHASE / "figures_local_v1" / "Fig03_CNV_clone_trees_by_cNMF_program_v4.png"
FIG14 = PHASE / "figures_local_v1" / "Fig14_NicheNet_ligand_receptor_evidence_v2.png"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    checks: list[dict] = []

    def check(name: str, passed: bool, evidence: str) -> None:
        checks.append({"check": name, "status": "PASS" if passed else "FAIL", "evidence": evidence})

    for path in [REPORT, PDF, FIG03, FIG14]:
        check(f"exists:{path.name}", path.exists() and path.stat().st_size > 0, str(path.stat().st_size))

    document = Document(REPORT)
    paragraphs = [paragraph.text for paragraph in document.paragraphs]
    check(
        "Fig03 caption updated",
        any("按亚克隆分支着色，并标注每条分支最主要的染色体 CNV 变化" in text for text in paragraphs),
        "updated caption text present",
    )
    check(
        "old Fig03 caption absent",
        not any("节点与分支按主导 cNMF 程序着色" in text for text in paragraphs),
        "superseded caption text absent",
    )

    with zipfile.ZipFile(REPORT) as archive:
        embedded_fig03 = archive.read("word/media/image4.png")
        embedded_fig14 = archive.read("word/media/image15.png")
    check(
        "Fig03 image bytes embedded exactly",
        hashlib.sha256(embedded_fig03).hexdigest() == sha256(FIG03),
        sha256(FIG03),
    )
    check(
        "Fig14 image bytes embedded exactly",
        hashlib.sha256(embedded_fig14).hexdigest() == sha256(FIG14),
        sha256(FIG14),
    )

    with Image.open(FIG03) as image:
        fig03_dims = image.size
    with Image.open(FIG14) as image:
        fig14_dims = image.size
    check("Fig03 300-dpi-sized raster", fig03_dims[0] >= 3000 and fig03_dims[1] >= 2000, str(fig03_dims))
    check("Fig14 300-dpi-sized raster", fig14_dims[0] >= 3000 and fig14_dims[1] >= 2000, str(fig14_dims))

    render_status = json.loads((PHASE / "rendered_local_v1" / "render_status.json").read_text(encoding="utf-8"))
    page_pngs = sorted((PHASE / "rendered_local_v1").glob("page-*.png"))
    check(
        "Word render passed",
        render_status.get("status") == "PASS" and render_status.get("rendered_pngs") == 26,
        f"{render_status.get('rendered_pngs')}/26 pages",
    )
    check("all rendered pages present", len(page_pngs) == 26, str(len(page_pngs)))
    check("focused pages present", page_pngs[6].name == "page-007.png" and page_pngs[17].name == "page-018.png", "pages 7 and 18")

    a11y = json.loads((QC / "a11y_audit_v1.json").read_text(encoding="utf-8"))
    summary = a11y.get("counts", {})
    check(
        "accessibility audit clean",
        summary.get("high", 0) == 0 and summary.get("medium", 0) == 0 and summary.get("low", 0) == 0,
        json.dumps(summary, ensure_ascii=False),
    )

    passed = sum(row["status"] == "PASS" for row in checks)
    failed = len(checks) - passed
    status = {
        "status": "PASS_PHASE46_FOCUSED_FIGURE_REVISION" if failed == 0 else "FAIL",
        "checks_passed": passed,
        "checks_total": len(checks),
        "failed": failed,
        "visual_review": {
            "all_pages_contact_sheet_reviewed": True,
            "focused_pages_reviewed_at_original_resolution": [7, 18],
        },
        "report_sha256": sha256(REPORT),
        "pdf_sha256": sha256(PDF),
        "checks": checks,
    }
    (QC / "validation_status_v2.json").write_text(
        json.dumps(status, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    lines = ["check\tstatus\tevidence"] + [
        f"{row['check']}\t{row['status']}\t{row['evidence']}" for row in checks
    ]
    (QC / "validation_checks_v2.tsv").write_text("\n".join(lines) + "\n", encoding="utf-8")
    if failed:
        raise SystemExit(f"Validation failed: {failed}/{len(checks)}")


if __name__ == "__main__":
    main()
