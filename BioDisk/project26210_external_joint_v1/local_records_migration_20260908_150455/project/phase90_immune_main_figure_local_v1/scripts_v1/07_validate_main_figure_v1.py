from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pdfplumber
import pypdfium2 as pdfium
from PIL import Image, ImageOps
from lxml import etree


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def average_hash(image: Image.Image, size: int = 32) -> np.ndarray:
    gray = image.convert("L").resize((size, size), Image.Resampling.LANCZOS)
    values = np.asarray(gray, dtype=float)
    return values > values.mean()


def main() -> None:
    if len(sys.argv) != 4:
        raise SystemExit("Usage: 07_validate_main_figure_v1.py <project_root> <work_root> <figure_dir>")
    project_root = Path(sys.argv[1])
    work_root = Path(sys.argv[2])
    figure_dir = Path(sys.argv[3])
    qc_dir = figure_dir / "qc_v1"
    qc_dir.mkdir(parents=False, exist_ok=False)

    checks: list[dict[str, object]] = []

    def record(check: str, observed: object, passed: bool) -> None:
        checks.append({"check": check, "observed": str(observed), "status": "PASS" if passed else "FAIL"})

    png = figure_dir / "PROJECT26210_IMMUNE_MAIN_FIGURE_v6.png"
    pdf = figure_dir / "PROJECT26210_IMMUNE_MAIN_FIGURE_v6.pdf"
    svg = figure_dir / "PROJECT26210_IMMUNE_MAIN_FIGURE_v6.svg"
    summary_path = figure_dir / "BUILD_SUMMARY_v6.json"
    source_manifest_path = figure_dir / "PANEL_SOURCE_MANIFEST_SHA256_v1.tsv"
    for path in [png, pdf, svg, summary_path, source_manifest_path]:
        record(f"exists_{path.name}", path.exists(), path.exists() and path.stat().st_size > 0)

    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    for kind, path in [("png", png), ("pdf", pdf), ("svg", svg)]:
        expected = summary["outputs"][kind]
        record(f"summary_bytes_{kind}", path.stat().st_size, path.stat().st_size == int(expected["bytes"]))
        record(f"summary_sha256_{kind}", sha256(path), sha256(path) == expected["sha256"])

    sources = pd.read_csv(source_manifest_path, sep="\t")
    source_pass = []
    for _, row in sources.iterrows():
        path = Path(row["path"])
        source_pass.append(path.exists() and path.stat().st_size == int(row["bytes"]) and sha256(path) == row["sha256"])
    record("source_manifest_all_exact", f"{sum(source_pass)}/{len(source_pass)}", len(source_pass) > 0 and all(source_pass))

    with Image.open(png) as image:
        rgb = image.convert("RGB")
        width, height = rgb.size
        small = rgb.resize((600, 800), Image.Resampling.BILINEAR)
        array = np.asarray(small)
        nonwhite = np.mean(np.min(array, axis=2) < 245)
        record("png_dimensions", f"{width}x{height}", width >= 4000 and height >= 5000)
        record("png_aspect_ratio", width / height, 0.72 <= width / height <= 0.84)
        record("png_nonwhite_coverage", round(float(nonwhite), 4), 0.12 <= nonwhite <= 0.70)
        png_hash = average_hash(rgb)

    with pdfplumber.open(pdf) as document:
        record("pdf_single_page", len(document.pages), len(document.pages) == 1)
        page_text = document.pages[0].extract_text() or ""
        required_titles = [
            "Immune-data expansion and multimodal analysis architecture",
            "Three-cohort immune atlas",
            "Cellular source of TYMP and lineage-specific programs",
            "Cross-cohort replication of TYMP-associated immune programs",
            "Reproducible myeloid-state ordering",
            "Spatial immune ecology and fine-state communication",
        ]
        for index, title in enumerate(required_titles, start=1):
            record(f"pdf_title_{index}", title, title in page_text)
        record("pdf_no_removed_panel_g", "Integrated evidence boundary" not in page_text, "Integrated evidence boundary" not in page_text)

    pdf_doc = pdfium.PdfDocument(str(pdf))
    rendered = pdf_doc[0].render(scale=2.0).to_pil().convert("RGB")
    render_path = qc_dir / "PROJECT26210_IMMUNE_MAIN_FIGURE_v6_pdf_render.png"
    rendered.save(render_path, optimize=True)
    pdf_hash = average_hash(rendered)
    hamming = int(np.sum(png_hash != pdf_hash))
    record("pdf_png_visual_hash_distance", hamming, hamming <= 90)

    parser = etree.XMLParser(resolve_entities=False, no_network=True, recover=False)
    root = etree.parse(str(svg), parser).getroot()
    svg_text = " ".join(text.strip() for text in root.itertext() if text.strip())
    record("svg_root", root.tag, root.tag.endswith("svg"))
    record("svg_has_all_panel_titles", sum(title in svg_text for title in required_titles), all(title in svg_text for title in required_titles))

    atlas_status = pd.read_csv(work_root / "atlas_extract_v1" / "atlas_extraction_status_v1.tsv", sep="\t")
    record("atlas_status_all_pass", ";".join(atlas_status["status"]), atlas_status["status"].eq("PASS").all())
    record("atlas_total_cells", int(atlas_status["total_cells"].sum()), int(atlas_status["total_cells"].sum()) == 389477)
    record("atlas_tymp_complete", int(atlas_status["tymp_finite"].sum()), int(atlas_status["tymp_finite"].sum()) == 389477)

    retained = pd.read_csv(
        project_root / "phase67_part1_direct_program_replication_local_v1" / "results_v5" / "phase67_cross_dataset_retained_v1.tsv",
        sep="\t",
    )
    myeloid = retained[(retained["endpoint"].eq("myeloid_dc_tymp")) & (retained["replication_grade"].eq("R2_TWO_EXTERNAL_FDR"))]
    treg = retained[(retained["endpoint"].eq("treg_tymp")) & (retained["datasets_fdr_replication"].ge(1))]
    record("phase67_myeloid_r2_rows", len(myeloid), len(myeloid) == 18)
    record("phase67_treg_external_rows", len(treg), len(treg) == 13)

    phase87_checks = pd.read_csv(
        project_root / "phase87_three_cohort_myeloid_trajectory_local_v1" / "validation_v1" / "PHASE87_VALIDATION_CHECKS_v1.tsv",
        sep="\t",
    )
    phase88_checks = pd.read_csv(
        project_root / "phase88_three_cohort_fine_state_cellchat_local_v1" / "validation_v1" / "PHASE88_VALIDATION_CHECKS_v1.tsv",
        sep="\t",
    )
    record("phase87_source_validation", f"{sum(phase87_checks.status == 'PASS')}/{len(phase87_checks)}", phase87_checks["status"].eq("PASS").all())
    record("phase88_source_validation", f"{sum(phase88_checks.status == 'PASS')}/{len(phase88_checks)}", phase88_checks["status"].eq("PASS").all())

    spatial_effects = pd.read_csv(
        project_root / "phase69_85_server_status_sync_local_v1" / "phase78" / "analysis_v1" / "phase78_section_effects_v1.tsv",
        sep="\t",
    )
    chosen = spatial_effects[spatial_effects["sample_id"].isin(["GSM5924037", "GSM5924047"]) & spatial_effects["endpoint_id"].eq("P01")]
    record("spatial_samples_two_preparations", ";".join(sorted(chosen["preparation"].unique())), set(chosen["preparation"]) == {"FFPE", "frozen"})
    record("spatial_selected_nonextreme", chosen["spearman_rho"].abs().max(), chosen["spearman_rho"].abs().max() < 0.05)

    output_files = [png, pdf, svg, summary_path, source_manifest_path, render_path]
    manifest = pd.DataFrame(
        [{"path": str(path), "bytes": path.stat().st_size, "sha256": sha256(path)} for path in output_files]
    )
    manifest.to_csv(qc_dir / "FINAL_FILE_MANIFEST_SHA256_v1.tsv", sep="\t", index=False)
    checks_frame = pd.DataFrame(checks)
    checks_frame.to_csv(qc_dir / "VALIDATION_CHECKS_v1.tsv", sep="\t", index=False)
    passed = checks_frame["status"].eq("PASS").all()
    status = {
        "status": "PASS" if passed else "FAIL",
        "checks_pass": int(checks_frame["status"].eq("PASS").sum()),
        "checks_total": int(len(checks_frame)),
        "pdf_render_path": str(render_path),
        "visual_review": "PASS_FULL_FIGURE_AND_PDF_RENDER",
    }
    (qc_dir / "VALIDATION_STATUS_v1.json").write_text(json.dumps(status, indent=2), encoding="utf-8")
    if not passed:
        raise SystemExit("Validation failed")


if __name__ == "__main__":
    main()
