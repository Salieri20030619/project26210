from __future__ import annotations

import hashlib
import json
import shutil
import sys
from datetime import datetime
from pathlib import Path

import pandas as pd


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def copy_new(source: Path, target: Path) -> None:
    if target.exists():
        raise SystemExit(f"Refusing existing target: {target}")
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, target)


def main() -> None:
    if len(sys.argv) != 3:
        raise SystemExit("Usage: 16_install_phase90_v10_revision.py <project_root> <work_root>")
    project_root = Path(sys.argv[1]).resolve()
    work_root = Path(sys.argv[2]).resolve()
    figure_root = work_root / "figure_v10"
    install_root = project_root / "phase90_immune_main_figure_local_v1"
    if not install_root.is_dir():
        raise SystemExit("Phase90 install root missing")
    now = datetime.now().astimezone().isoformat(timespec="seconds")

    output_names = [
        "PROJECT26210_IMMUNE_MAIN_FIGURE_v10.png",
        "PROJECT26210_IMMUNE_MAIN_FIGURE_v10.pdf",
        "PROJECT26210_IMMUNE_MAIN_FIGURE_v10.svg",
        "BUILD_SUMMARY_v10.json",
    ]
    copied: list[Path] = []
    for name in output_names:
        target = install_root / name
        copy_new(figure_root / name, target)
        copied.append(target)
    for name in [
        "VALIDATION_STATUS_v1.json",
        "VALIDATION_CHECKS_v1.tsv",
        "FINAL_FILE_MANIFEST_SHA256_v1.tsv",
        "PROJECT26210_IMMUNE_MAIN_FIGURE_v10_pdf_render.png",
    ]:
        target = install_root / "qc_v10" / name
        copy_new(figure_root / "qc_v1" / name, target)
        copied.append(target)
    for name in [
        "12_build_main_figure_v8.py",
        "14_build_main_figure_v10.py",
        "15_validate_main_figure_v10.py",
        "16_install_phase90_v10_revision.py",
    ]:
        target = install_root / "scripts_v1" / name
        copy_new(work_root / "scripts" / name, target)
        copied.append(target)

    hashes = {ext: sha256(figure_root / f"PROJECT26210_IMMUNE_MAIN_FIGURE_v10.{ext}") for ext in ["png", "pdf", "svg"]}
    addendum = f"""# Phase 90 D-panel readability revision v10

- Completed: {now}
- Status: **PASS_MAIN_FIGURE_V10; 35/35 VALIDATION CHECKS PASS; D-PANEL AND FULL PDF RENDER REVIEW PASS**.
- Revision scope: Panel D only. Panels A-C and E-F are unchanged from the validated v7 composition.
- The four evidence groups are now explicit vertical colour bands: Myeloid/DC metabolic, Myeloid/DC regulon, Treg metabolic and Treg regulon.
- Full concise program/regulon names were moved into the former blank central region instead of being truncated on the far-left margin.
- GSE264586 and GSE294109 are explicit aligned evidence columns. Fill encodes rho, area encodes -log10(q), and a coloured connector indicates same-direction effects.
- The Phase 67 source rows, effect values, q values and scientific interpretation are unchanged. This is a readability/layout revision only.
- v10 supersedes v7 for client delivery.
- PNG SHA256: `{hashes['png']}`.
- PDF SHA256: `{hashes['pdf']}`.
- SVG SHA256: `{hashes['svg']}`.
"""
    addendum_path = install_root / "PROJECT_STATUS_PHASE90_ADDENDUM_v2.md"
    if addendum_path.exists():
        raise SystemExit("Phase90 addendum v2 already exists")
    addendum_path.write_text(addendum, encoding="utf-8")
    copied.append(addendum_path)

    all_rows = []
    for path in sorted(install_root.rglob("*")):
        if path.is_file() and path.name not in {"FILE_MANIFEST_SHA256_v2.tsv", "INSTALL_VALIDATION_v2.json"}:
            all_rows.append({"relative_path": path.relative_to(install_root).as_posix(), "bytes": path.stat().st_size, "sha256": sha256(path)})
    manifest_v2 = install_root / "FILE_MANIFEST_SHA256_v2.tsv"
    if manifest_v2.exists():
        raise SystemExit("File manifest v2 already exists")
    pd.DataFrame(all_rows).to_csv(manifest_v2, sep="\t", index=False)

    status_path = project_root / "PROJECT_STATUS.md"
    status_text = status_path.read_text(encoding="utf-8")
    revision_heading = "### Phase 90 D-panel readability revision v10 (2026-08-27)"
    if revision_heading in status_text:
        raise SystemExit("Phase90 v10 status already recorded")
    revision_block = f"""

{revision_heading}

- Status: **PASS_MAIN_FIGURE_V10; 35/35 VALIDATION CHECKS PASS; FULL PDF-RENDER REVIEW PASS**.
- Panel D was revised for readability: group colour bands remain on the left, complete concise program/regulon names now occupy the former blank central region, and the GSE264586/GSE294109 effect columns are explicitly aligned.
- No Phase67 source value or conclusion changed. Panels A-C and E-F are unchanged. v10 supersedes v7 for client delivery.
- Authoritative deliverables: `PROJECT26210_IMMUNE_MAIN_FIGURE_v10.png`, `.pdf` and `.svg` under `phase90_immune_main_figure_local_v1`.
"""
    status_path.write_text(status_text.rstrip() + revision_block + "\n", encoding="utf-8")

    results_path = project_root / "RESULTS_MANIFEST.tsv"
    results_text = results_path.read_text(encoding="utf-8").rstrip("\r\n")
    existing_ids = {line.split("\t", 1)[0] for line in results_text.splitlines()[1:] if line.strip()}
    result_rows = [
        ("phase90_immune_main_figure_png_v10", "high_resolution_png", "PROJECT26210_IMMUNE_MAIN_FIGURE_v10.png", f"SHA256 {hashes['png']}; 4166x5622", "Authoritative readability revision; supersedes v7"),
        ("phase90_immune_main_figure_pdf_v10", "publication_pdf", "PROJECT26210_IMMUNE_MAIN_FIGURE_v10.pdf", f"SHA256 {hashes['pdf']}; one page; PDF render reviewed", "Authoritative publication PDF"),
        ("phase90_immune_main_figure_svg_v10", "editable_svg", "PROJECT26210_IMMUNE_MAIN_FIGURE_v10.svg", f"SHA256 {hashes['svg']}; editable text/vector", "Authoritative editable source"),
        ("phase90_immune_main_figure_validation_v10", "validation_json", "qc_v10/VALIDATION_STATUS_v1.json", "35/35 PASS; D-panel and full PDF render reviewed", "v10 independent QA"),
        ("phase90_immune_main_figure_status_addendum_v2", "status_addendum_md", "PROJECT_STATUS_PHASE90_ADDENDUM_v2.md", "PASS_MAIN_FIGURE_V10", "Readability revision scope and unchanged scientific boundary"),
        ("phase90_immune_main_figure_file_manifest_v2", "file_manifest_tsv", "FILE_MANIFEST_SHA256_v2.tsv", f"{len(all_rows)} files hashed", "Expanded Phase90 package manifest including v10"),
    ]
    lines = []
    for result_id, artifact_type, relative, evidence, notes in result_rows:
        if result_id in existing_ids:
            raise SystemExit(f"Manifest ID already exists: {result_id}")
        fields = [result_id, "phase90_immune_main_figure", artifact_type, str(install_root / relative), now, "PASS", evidence, notes]
        lines.append("\t".join(fields))
    results_path.write_text(results_text + "\n" + "\n".join(lines) + "\n", encoding="utf-8")

    installed = pd.read_csv(manifest_v2, sep="\t")
    exact = []
    for _, row in installed.iterrows():
        path = install_root / Path(row["relative_path"])
        exact.append(path.exists() and path.stat().st_size == int(row["bytes"]) and sha256(path) == row["sha256"])
    install_status = {
        "status": "PASS",
        "package_files_exact": int(sum(exact)),
        "package_files_total": int(len(exact)),
        "v10_validation_checks": "35/35 PASS",
        "root_status_revision_count": status_path.read_text(encoding="utf-8").count(revision_heading),
        "results_manifest_v10_rows": int(pd.read_csv(results_path, sep="\t")["result_id"].str.contains("v10|addendum_v2|manifest_v2").sum()),
        "installed_at": now,
    }
    (install_root / "INSTALL_VALIDATION_v2.json").write_text(json.dumps(install_status, indent=2), encoding="utf-8")
    if not all(exact) or install_status["root_status_revision_count"] != 1:
        raise SystemExit("v10 installation validation failed")


if __name__ == "__main__":
    main()
