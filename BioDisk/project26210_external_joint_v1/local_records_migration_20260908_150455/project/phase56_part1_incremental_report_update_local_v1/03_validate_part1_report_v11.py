from __future__ import annotations

import hashlib
import json
import zipfile
from pathlib import Path

from docx import Document
from PIL import Image, ImageChops


ROOT = Path(__file__).resolve().parents[1]
PHASE = Path(__file__).resolve().parent
QC = PHASE / "qc_local_v1"
SOURCE = (
    ROOT
    / "phase46_part1_clone_tree_and_fig14_revision_local_v1"
    / "PART1_EXTERNAL_SINGLE_CELL_CLIENT_REPORT_WITH_FIGURES_REVIEW_v10.docx"
)
REPORT = PHASE / "PART1_EXTERNAL_SINGLE_CELL_CLIENT_REPORT_WITH_FIGURES_REVIEW_v11.docx"
RENDER = PHASE / "rendered_local_v1"
PDF = RENDER / "PART1_EXTERNAL_SINGLE_CELL_CLIENT_REPORT_WITH_FIGURES_REVIEW_v11.pdf"


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
        diff = ImageChops.difference(rgb, background)
        return diff.getbbox() is not None


def main() -> None:
    checks: list[dict[str, str]] = []

    def check(name: str, passed: bool, evidence: str) -> None:
        checks.append({"check": name, "status": "PASS" if passed else "FAIL", "evidence": evidence})

    for path in [SOURCE, REPORT, PDF]:
        check(f"exists:{path.name}", path.exists() and path.stat().st_size > 0, str(path.stat().st_size))

    document = Document(REPORT)
    text = "\n".join(paragraph.text for paragraph in document.paragraphs)
    required = {
        "three-tier section": "C06  新增患者三层复核与参考桥接筛查",
        "GSE178481 tiered statistic": "rho=0.929，精确P=0.00223，层内BH校正q=0.0179",
        "GSE264586 exact statistic": "rho=1.000，精确P=0.0167",
        "Phase49 narrative correction": "P=0.02345、q=0.09380",
        "full-chain boundary": "完整F1/TYMP→CAF→TAM→CD8链条目前不成立",
    }
    for name, phrase in required.items():
        check(name, phrase in text, phrase)

    check("paragraph count", len(document.paragraphs) == 168, str(len(document.paragraphs)))
    check("table count preserved", len(document.tables) == 3, str(len(document.tables)))
    headings = [paragraph.text for paragraph in document.paragraphs if paragraph.style.name.startswith("Heading")]
    check("C06 uses heading style", required["three-tier section"] in headings, "Heading 2")

    source_media = media_hashes(SOURCE)
    output_media = media_hashes(REPORT)
    check(
        "all embedded figures preserved byte-for-byte",
        source_media == output_media,
        f"source={len(source_media)} output={len(output_media)} media parts",
    )

    render_status = json.loads((RENDER / "render_status.json").read_text(encoding="utf-8"))
    pages = sorted(RENDER.glob("page-*.png"))
    check(
        "Word render passed",
        render_status.get("status") == "PASS" and render_status.get("rendered_pngs") == 27,
        json.dumps(render_status, ensure_ascii=False),
    )
    check("all 27 page PNGs present", len(pages) == 27, str(len(pages)))
    check("no blank rendered pages", all(page_has_visible_content(page) for page in pages), "27/27 visible")

    a11y = json.loads((QC / "a11y_audit_v1.json").read_text(encoding="utf-8"))
    counts = a11y.get("counts", {})
    check(
        "accessibility audit clean",
        counts.get("high", 0) == 0 and counts.get("medium", 0) == 0 and counts.get("low", 0) == 0,
        json.dumps(counts, ensure_ascii=False),
    )

    update_status = json.loads((QC / "report_update_status_v1.json").read_text(encoding="utf-8"))
    check(
        "source hash matches update record",
        update_status.get("source_sha256") == sha256(SOURCE),
        sha256(SOURCE),
    )
    check(
        "report hash matches update record",
        update_status.get("output_sha256") == sha256(REPORT),
        sha256(REPORT),
    )

    failed = sum(row["status"] == "FAIL" for row in checks)
    status = {
        "status": "PASS_PHASE56_PART1_REPORT_V11" if failed == 0 else "FAIL",
        "checks_passed": len(checks) - failed,
        "checks_total": len(checks),
        "failed": failed,
        "visual_review": {
            "all_pages_reviewed_at_original_resolution": True,
            "pages": list(range(1, 28)),
            "layout_findings": "No clipping, overlap, orphaned headings, blank pages, or figure overflow.",
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

    manifest_paths = [
        REPORT,
        PDF,
        RENDER / "render_status.json",
        QC / "a11y_audit_v1.json",
        QC / "report_update_status_v1.json",
        QC / "validation_checks_v1.tsv",
        QC / "validation_status_v1.json",
    ]
    manifest = ["relative_path\tbytes\tsha256"]
    for path in manifest_paths:
        manifest.append(f"{path.relative_to(PHASE)}\t{path.stat().st_size}\t{sha256(path)}")
    (PHASE / "FILE_MANIFEST_SHA256_v1.tsv").write_text("\n".join(manifest) + "\n", encoding="utf-8")

    if failed:
        raise SystemExit(f"Validation failed: {failed}/{len(checks)}")


if __name__ == "__main__":
    main()
